//***********************************************************************
// csatlakozó face osztályok header
// Creation date:  2018. 08. 31.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef FACE_CSATLAKOZO_HEADER
#define	FACE_CSATLAKOZO_HEADER
//***********************************************************************


//***********************************************************************
#include "cella_es_face.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
class dc_csatlakozo_face : public dc_keresztface {
//***********************************************************************
public:
    //***********************************************************************
    dc_csatlakozo_face(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_keresztface(p_tok, p_cella, p_adat_face, is_el) {}
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_keresztface::init();
    }
    //***********************************************************************
    void update_UT_IP_NR() override {
    //***********************************************************************
        //if (*p_tok->p_UTi > rvt(10))
        //    *p_tok->p_UTi = rvt(10);
        //if (*p_tok->p_UTi < rvt(-10))
        //    *p_tok->p_UTi = rvt(-10);
        p_tok->emlekek.get_akt().UT = p_tok->emlekek.get_kiindulo().UT + *p_tok->p_UTi * akt_sim.alfa;
        //p_tok->IP.get_akt() = p_tok->IP.get_kiindulo() + *p_tok->p_IPi * akt_sim.alfa; // nem tudom, hogy ennek van-e bármi értelme
    }
    //***********************************************************************
    void update_UT_IP_SA() override {
    //***********************************************************************
        p_tok->emlekek.get_akt().UT = *p_tok->p_UTi;
        p_tok->emlekek.get_akt().IP = *p_tok->p_IPi; // nem tudom, hogy ennek van-e bármi értelme
    }
    //***********************************************************************
    void update_hibaaram_NR_kereszt(rvt & sum_el_hiba, rvt & sum_th_hiba) override {
    //***********************************************************************
        rvt IP = p_tok->emlekek.get_akt().IP;
        *p_tok->p_Ji = IP;     // !! betesszük a teljes áramot az inhomogén áram vektorba, így két cella összevonásakor pont a hibaáram marad itt.
        csatlakozo_aramok_dc[csatlakozo_index] = IP;
        if (is_el) {
            sum_el_hiba -= IP; // ki kell vonni az áramot, mert a külsõ kapocsnak megfelelõ az iránya-
        }
        else {
            sum_th_hiba -= IP; // ki kell vonni az áramot, mert a külsõ kapocsnak megfelelõ az iránya 
        }
    }
    //***********************************************************************
    virtual void debug_write(::std::ofstream & fs, bool is_endl) const override { if (is_endl) fs << "csatlak  "; dc_keresztface::debug_write(fs, false); if (is_endl) fs << ::std::endl; }
    //***********************************************************************
};


//***********************************************************************
class dc_csatlakozo_face_ellenallas_elth_nem_disszipalo_const : public dc_csatlakozo_face {
// keresztface, paramétereit a p_adat_face-bõl veszi
//***********************************************************************
    rvt G; // csak inicializáláskor kap értéket: nem jó, ha emlékezõben van.
public:
    //***********************************************************************
    dc_csatlakozo_face_ellenallas_elth_nem_disszipalo_const(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_csatlakozo_face(p_tok, p_cella, p_adat_face, is_el), G{ rvt { g_min } } {}
    //***********************************************************************
    face_azonosito get_azon() const { return fa_cd_r_kdn; }
    //***********************************************************************
    bool is_konst() const { return true; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_csatlakozo_face::init();
        rvt fajlagos_vez = is_el 
            ? p_anyag->elvez.get_value(p_cella->get_Tc_megtartando(), ertek_t()).ertek
            : p_anyag->thvez.get_value(p_cella->get_Tc_megtartando(), ertek_t()).ertek;
        G = fajlagos_vez * A_per_L;
    }
    //***********************************************************************
    void update_peremfeltetel() override {}
    //***********************************************************************
    void update_uj_lepeshez() override {
    // kapacitások frissítése, konstansok admittanciamátrixokba írása, 
    // konstans középponti admittanciák növelése
    //***********************************************************************
        if (is_el) {
            *p_tok->p_Yii = G;
            *p_tok->p_Yie = -G;
            *p_tok->p_Yei = -G;
            p_cella->konst_Yee_dc += G;
        }
        else {
            *p_tok->p_Yii = G;
            *p_tok->p_Yit = -G;
            *p_tok->p_Yti = -G;
            p_cella->konst_Ytt_dc += G;
        }
    }
    //***********************************************************************
    void update_nemlinearis_parameterek() override {} // konstans
    //***********************************************************************
    void update_admittanciamatrix_SA() override {} // konstans
    //***********************************************************************
    void update_Jakobi_NR_nincs_csatolt() override {} // konstans
    //***********************************************************************
    void update_Jakobi_NR_van_csatolt() override {} // konstans
    //***********************************************************************
    void agaram_szamitasa() override {
    //***********************************************************************
        rvt UTc = is_el ? p_cella->el_center_face_dc.emlekek.get_akt().UT : p_cella->th_center_face_dc.emlekek.get_akt().UT;
        rvt V = p_tok->emlekek.get_akt().UT - UTc;
        p_tok->emlekek.get_akt().IP = G * V;
    }
};


//***********************************************************************
class dc_csatlakozo_face_ellenallas_elth_nem_disszipalo_Tfuggo : public dc_csatlakozo_face {
// keresztface, paramétereit a p_adat_face-bõl veszi
//***********************************************************************
    struct emlekek_strukt {
        rvt G, dG_per_dT;
    };
    emlek<emlekek_strukt> emlekek;
public:
    //***********************************************************************
    dc_csatlakozo_face_ellenallas_elth_nem_disszipalo_Tfuggo(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_csatlakozo_face(p_tok, p_cella, p_adat_face, is_el) { emlekek.get_akt().G = rvt(g_min); }
    //***********************************************************************
    face_azonosito get_azon() const { return fa_cd_r_tdn; }
    //***********************************************************************
    bool is_konst() const { return false; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_csatlakozo_face::init();
        ertek_t fajlagos_ertek = is_el
            ? p_anyag->elvez.get_value(p_cella->get_Tc_megtartando(), p_cella->emlekek.get_megtartando().H)
            : p_anyag->thvez.get_value(p_cella->get_Tc_megtartando(), p_cella->emlekek.get_megtartando().H);
        emlekek_strukt & eml = emlekek.get_akt();
        eml.G = fajlagos_ertek.ertek * A_per_L;
        eml.dG_per_dT = fajlagos_ertek.derivalt * A_per_L;
    }
    //***********************************************************************
    void update_peremfeltetel() override {}
    //***********************************************************************
    void update_uj_lepeshez() override {} // nem konstans
    //***********************************************************************
    void update_nemlinearis_parameterek() override {
    //***********************************************************************
        ertek_t fajlagos_ertek = is_el
            ? p_anyag->elvez.get_value(p_cella->get_Tc_akt(), p_cella->emlekek.get_akt().H)
            : p_anyag->thvez.get_value(p_cella->get_Tc_akt(), p_cella->emlekek.get_akt().H);
        emlekek_strukt & eml = emlekek.get_akt();
        eml.G = fajlagos_ertek.ertek * A_per_L;
        eml.dG_per_dT = fajlagos_ertek.derivalt * A_per_L;
    }
    //***********************************************************************
    void update_admittanciamatrix_SA() override {
    //***********************************************************************
        rvt Gelozo = emlekek.get_elozo().G;
        if (is_el) {
            *p_tok->p_Yii = Gelozo;
            *p_tok->p_Yie = -Gelozo;
            *p_tok->p_Yei = -Gelozo;
            *p_cella->p_Yee_dc += Gelozo;
        }
        else {
            *p_tok->p_Yii = Gelozo;
            *p_tok->p_Yit = -Gelozo;
            *p_tok->p_Yti = -Gelozo;
            *p_cella->p_Ytt_dc += Gelozo;
        }
    }
    //***********************************************************************
    void set_cella_is_szimm() {
    //***********************************************************************
        if (is_el)
            return; // csak termikusnál nemszimmetrikus
        // TODO: kikapcsolom a Newton-Raphson algoritmust termikus szimulációnál => iteráció. Eddigi tapasztalat szerint nem hoz a konyhára.
        // p_cella->set_to_NR_nincs_csatolt_nemszimm();
        p_cella->set_to_NR_van_csatolt_nemszimm();
    }
    //***********************************************************************
    void update_Jakobi_NR_nincs_csatolt() override {
    //***********************************************************************
        const emlekek_strukt & eml = emlekek.get_elozo();
        rvt Gelozo = eml.G;
        if (is_el) {
            *p_tok->p_Yii = Gelozo;
            *p_tok->p_Yie = -Gelozo;
            *p_tok->p_Yei = -Gelozo;
            *p_cella->p_Yee_dc += Gelozo;
        }
        else {
            // TODO: kikapcsolom a Newton-Raphson algoritmust termikus szimulációnál => iteráció. Eddigi tapasztalat szerint nem hoz a konyhára.
            // rvt dT = p_tok->UT.get_elozo() - p_cella->th_center_face_dc.UT.get_elozo(); // a teljes esõ hõmérsékletet kell venni
            rvt dT = 0;
            *p_tok->p_Yii = Gelozo;
            *p_tok->p_Yit = eml.dG_per_dT*dT - Gelozo;
            *p_tok->p_Yti = -Gelozo;
            *p_cella->p_Ytt_dc += Gelozo - eml.dG_per_dT*dT;
        }
    }
    //***********************************************************************
    void update_Jakobi_NR_van_csatolt() override {
    //***********************************************************************
        const emlekek_strukt & eml = emlekek.get_elozo();
        rvt Gelozo = eml.G;
        if (is_el) { // nincs disszipáció
            *p_tok->p_Yii = Gelozo;
            *p_tok->p_Yie = -Gelozo;
            *p_tok->p_Yei = -Gelozo;
            *p_cella->p_Yee_dc += Gelozo;
        }
        else {
            rvt dT = p_tok->emlekek.get_elozo().UT - p_cella->th_center_face_dc.emlekek.get_elozo().UT; // a teljes esõ hõmérsékletet kell venni, nem a hibahõmérsékletet!!!
            *p_tok->p_Yii = Gelozo;
            *p_tok->p_Yit = eml.dG_per_dT*dT - Gelozo;
            *p_tok->p_Yti = -Gelozo;
            *p_cella->p_Ytt_dc += Gelozo - eml.dG_per_dT*dT;
        }
    }
    //***********************************************************************
    void agaram_szamitasa() override {
    //***********************************************************************
        rvt UTc = is_el ? p_cella->el_center_face_dc.emlekek.get_akt().UT : p_cella->th_center_face_dc.emlekek.get_akt().UT;
        rvt V = p_tok->emlekek.get_akt().UT - UTc;
        p_tok->emlekek.get_akt().IP = emlekek.get_akt().G * V;
    }
};


//***********************************************************************
class dc_csatlakozo_face_ellenallas_el_disszipalo_const : public dc_csatlakozo_face {
// keresztface, paramétereit a p_adat_face-bõl veszi
//***********************************************************************
    rvt G; // csak inicializáláskor kap értéket: nem jó, ha emlékezõben van.
public:
    //***********************************************************************
    dc_csatlakozo_face_ellenallas_el_disszipalo_const(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_csatlakozo_face(p_tok, p_cella, p_adat_face, is_el), G{ rvt{ g_min } } {}
    //***********************************************************************
    face_azonosito get_azon() const { return fa_cd_r_ked; }
    //***********************************************************************
    bool is_konst() const { return true; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_csatlakozo_face::init();
        rvt fajlagos_vez = p_anyag->elvez.get_value(p_cella->get_Tc_megtartando(), ertek_t()).ertek;
        G = fajlagos_vez * A_per_L;
    }
    //***********************************************************************
    void update_peremfeltetel() override {}
    //***********************************************************************
    void update_uj_lepeshez() override {
    // kapacitások frissítése, konstansok admittanciamátrixokba írása, 
    // konstans középponti admittanciák növelése
    //***********************************************************************
        *p_tok->p_Yii = G;
        *p_tok->p_Yie = -G;
        *p_tok->p_Yei = -G;
        p_cella->konst_Yee_dc += G;
    }
    //***********************************************************************
    void update_nemlinearis_parameterek() override {} // konstans
    //***********************************************************************
    void update_admittanciamatrix_SA() override {} // konstans
    //***********************************************************************
    void set_cella_is_szimm() {
    //***********************************************************************
        p_cella->set_to_NR_van_csatolt_nemszimm();
    }
    //***********************************************************************
    void update_Jakobi_NR_nincs_csatolt() override {} // konstans
    //***********************************************************************
    void update_Jakobi_NR_van_csatolt() override {  // csak a disszipáció nem konstans
    //***********************************************************************
        rvt V = p_tok->emlekek.get_elozo().UT - p_cella->el_center_face_dc.emlekek.get_elozo().UT; // a teljes esõ feszültséget kell venni, nem a hibafeszültséget!!!
        rvt dP_per_dV = 2*G*V;
        *p_tok->p_Yti = -dP_per_dV; // nem áramfüggõ ellenállásnál -dP/dV = -d(GV^2)/dV = -2GV
        *p_cella->p_Yte_dc += dP_per_dV;
    }
    //***********************************************************************
    void agaram_szamitasa() override {
    //***********************************************************************
        rvt Uc = p_cella->el_center_face_dc.emlekek.get_akt().UT;
        rvt V = p_tok->emlekek.get_akt().UT - Uc;
        rvt I = G * V;
        p_tok->emlekek.get_akt().IP = I;
        *p_cella->p_sum_disszipacio_dc += V * I;
    }
};


//***********************************************************************
class dc_csatlakozo_face_ellenallas_el_disszipalo_Tfuggo : public dc_csatlakozo_face {
// keresztface, paramétereit a p_adat_face-bõl veszi
//***********************************************************************
    struct emlekek_strukt {
        rvt G, dG_per_dT;
    };
    emlek<emlekek_strukt> emlekek;
public:
    //***********************************************************************
    dc_csatlakozo_face_ellenallas_el_disszipalo_Tfuggo(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_csatlakozo_face(p_tok, p_cella, p_adat_face, is_el) { emlekek.get_akt().G = rvt(g_min); }
    //***********************************************************************
    face_azonosito get_azon() const { return fa_cd_r_ted; }
    //***********************************************************************
    bool is_konst() const { return false; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_csatlakozo_face::init();
        ertek_t fajlagos_ertek = p_anyag->elvez.get_value(p_cella->get_Tc_megtartando(), p_cella->emlekek.get_megtartando().H);
        emlekek_strukt & eml = emlekek.get_akt();
        eml.G = fajlagos_ertek.ertek * A_per_L;
        eml.dG_per_dT = fajlagos_ertek.derivalt * A_per_L;
    }
    //***********************************************************************
    void update_peremfeltetel() override {}
    //***********************************************************************
    void update_uj_lepeshez() override {} // nem konstans
    //***********************************************************************
    void update_nemlinearis_parameterek() override {
    //***********************************************************************
        ertek_t fajlagos_ertek = p_anyag->elvez.get_value(p_cella->get_Tc_akt(), p_cella->emlekek.get_akt().H);
        emlekek_strukt & eml = emlekek.get_akt();
        eml.G = fajlagos_ertek.ertek * A_per_L;
        eml.dG_per_dT = fajlagos_ertek.derivalt * A_per_L;
    }
    //***********************************************************************
    void update_admittanciamatrix_SA() override {
    //***********************************************************************
        rvt Gelozo = emlekek.get_elozo().G;
        *p_tok->p_Yii = Gelozo;
        *p_tok->p_Yie = -Gelozo;
        *p_tok->p_Yei = -Gelozo;
        *p_cella->p_Yee_dc += Gelozo;
    }
    //***********************************************************************
    void set_cella_is_szimm() {
    //***********************************************************************
        p_cella->set_to_NR_van_csatolt_nemszimm();
    }
    //***********************************************************************
    void update_Jakobi_NR_nincs_csatolt() override {
    //***********************************************************************
        rvt Gelozo = emlekek.get_elozo().G;
        *p_tok->p_Yii = Gelozo;
        *p_tok->p_Yie = -Gelozo;
        *p_tok->p_Yei = -Gelozo;
        *p_cella->p_Yee_dc += Gelozo;
    }
    //***********************************************************************
    void update_Jakobi_NR_van_csatolt() override {
    //***********************************************************************
        const emlekek_strukt & eml = emlekek.get_elozo();
        rvt Gelozo = eml.G;
        rvt V = p_tok->emlekek.get_elozo().UT - p_cella->el_center_face_dc.emlekek.get_elozo().UT; // a teljes esõ feszültséget kell venni, nem a hibafeszültséget!!!
        rvt dP_per_dV = 2 * Gelozo*V;
        *p_tok->p_Yii = Gelozo;
        *p_tok->p_Yie = -Gelozo;
        *p_tok->p_Yit = eml.dG_per_dT*V;
        *p_tok->p_Yei = -Gelozo;
        *p_tok->p_Yti = -dP_per_dV; // nem áramfüggõ ellenállásnál -dP/dV = -d(GV^2)/dV = -2GV
        *p_cella->p_Yee_dc += Gelozo;
        *p_cella->p_Yet_dc -= eml.dG_per_dT*V;
        *p_cella->p_Yte_dc += dP_per_dV;
        // A tt tagot a cella::fw_klaszter_1_update_Jakobi_dc adja hozzá
    }
    //***********************************************************************
    void agaram_szamitasa() override {
    //***********************************************************************
        const emlekek_strukt & eml = emlekek.get_akt();
        rvt Uc = p_cella->el_center_face_dc.emlekek.get_akt().UT;
        rvt V = p_tok->emlekek.get_akt().UT - Uc;
        rvt I = eml.G * V;
        p_tok->emlekek.get_akt().IP = I;
        *p_cella->p_sum_disszipacio_dc += V * I;
        *p_cella->p_sum_ddissz_per_dT_dc += eml.dG_per_dT * V * V;
    }
};

//***********************************************************************
class dc_csatlakozo_face_junction_el_disszipalo_Tfuggo : public dc_csatlakozo_face {
// keresztface, paramétereit a p_adat_face-bõl veszi
//***********************************************************************
    struct emlekek_strukt {
        rvt G_dioda, G_ellenallas, dG_ellenallas_per_dT, dI_dioda_per_dT, g_diff, g_diff_dioda, g_diff_ellenallas;
        rvt G, P_dissz_P_eso;
        rvt bU1, bUc, bUx, bTc;
        rvt U_kozep;
        LED_model_result_pack pack;
    };
    emlek<emlekek_strukt> emlekek;
    rvt Aface_per_Ajunction_full, Ajunction_per_Aface_full;
    rvt Aface_per_Ajunction_1, Ajunction_per_Aface_1;
public:
    //***********************************************************************
    dc_csatlakozo_face_junction_el_disszipalo_Tfuggo(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_csatlakozo_face(p_tok, p_cella, p_adat_face, is_el), Aface_per_Ajunction_full{ rvt() }, Ajunction_per_Aface_full{ rvt() }, 
        Aface_per_Ajunction_1{ rvt() }, Ajunction_per_Aface_1{ rvt() } { emlekek.get_akt().G = rvt(g_min); }
    //***********************************************************************
    face_azonosito get_azon() const { return fa_cd_j_ted; }
    //***********************************************************************
    bool is_konst() const { return false; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_csatlakozo_face::init();

        rvt Tc = p_cella->get_Tc_megtartando();
        emlekek_strukt & eml = emlekek.get_akt();

        // junction paraméterek

        if (p_junction == nullptr)
            throw hiba("dc_csatlakozo_face_junction_el_disszipalo_Tfuggo::init", "missing junction");
        //***********************************************************************
        // !!!!                                                                    !!!!
        // !!!! A COB LED-nél a LED-ek számával szorozzuk                          !!!!
        // !!!!                                                                    !!!!
        Aface_per_Ajunction_full = A / p_junction->area;
        Aface_per_Ajunction_1 = Aface_per_Ajunction_full * 12/*24*/;
        //***********************************************************************
        Ajunction_per_Aface_full = 1 / Aface_per_Ajunction_full;
        Ajunction_per_Aface_1 = 1 / Aface_per_Ajunction_1;
        //rvt fajlagos_diff_vezetes_j = rvt();
        //ertek_t fajlagos_ertek_j = p_junction->el_egyenlet.get_G_junction(Tc, rvt(), fajlagos_diff_vezetes_j);
        p_junction->el_egyenlet.junction_calc_el_from_U(2.1, Tc, 0.05, eml.pack, true);
        eml.G_dioda = eml.pack.G * Aface_per_Ajunction_1;
        eml.g_diff_dioda = eml.pack.dI_per_dVj * Aface_per_Ajunction_1;
        eml.dI_dioda_per_dT = eml.pack.dI_per_dT * Aface_per_Ajunction_1;
        eml.P_dissz_P_eso = 1; // A melegítõ (esõ-sugárzott) teljesítmény és a junctionon esõ teljesítmény aránya. Kezdetben 1.
        //P_dissz_P_eso.get_elozo() = 1;

        // ellenállás paraméterek

        ertek_t fajlagos_ertek_r = p_anyag->elvez.get_value(Tc, p_cella->emlekek.get_megtartando().H);
        eml.G_ellenallas = fajlagos_ertek_r.ertek * A_per_L;
        eml.g_diff_ellenallas = eml.G_ellenallas;
        eml.dG_ellenallas_per_dT = fajlagos_ertek_r.derivalt * A_per_L;

        // teljes face paraméterei

        eml.G = replusz(eml.G_ellenallas, eml.G_dioda);
        eml.g_diff = replusz(eml.g_diff_ellenallas, eml.g_diff_dioda);
        //printf("%g, %g, %g, %g, %g, %g\n", G, G_dioda, G_ellenallas, g_diff, g_diff_dioda, g_diff_ellenallas);
    }
    //***********************************************************************
    void update_peremfeltetel() override {}
    //***********************************************************************
    void update_uj_lepeshez() override {} // nem konstans
    //***********************************************************************
    void update_nemlinearis_parameterek() override {
    //***********************************************************************

        rvt Tc = p_cella->get_Tc_akt();
        emlekek_strukt & eml = emlekek.get_akt();

        eml.bTc = Tc;
        eml.bU1 = p_tok->emlekek.get_akt().UT;
        eml.bUx = eml.U_kozep;
        eml.bUc = p_cella->el_center_face_dc.emlekek.get_akt().UT;

        // junction paraméterek

        //rvt Vj = p_tok->UT.get_elozo() - U_kozep.get_elozo();
        rvt Vj = p_tok->emlekek.get_akt().UT - eml.U_kozep;
//        rvt fajlagos_diff_vezetes_j = rvt();
//        ertek_t fajlagos_ertek_j = p_junction->el_egyenlet.get_G_junction(Tc, Vj, fajlagos_diff_vezetes_j);

        rvt Uc = p_cella->el_center_face_dc.emlekek.get_akt().UT;
        rvt U1 = p_tok->emlekek.get_akt().UT;
        rvt Iface = (U1 - Uc) * eml.G; // a friss U-val számoljuk újra I-t, hogy konvergáljon. Az elõzõ I-vel nem konvergál.
        rvt I = eml.pack.get_votma() ? (Iface * Ajunction_per_Aface_1) : eml.pack.IF;
        //rvt I = pack.get_akt().get_votma() ? abs(p_tok->IP.get_akt() * Ajunction_per_Aface_1) : pack.get_akt().IF;
//printf("votma1=%u\n", (uns)pack.get_akt().get_votma());
        bool is_votma = eml.pack.get_votma();
        //printf("1 %g\n", p_tok->IP.get_akt());
        //printf("%g %g\n", I, Ajunction_per_Aface_1);
        p_junction->el_egyenlet.junction_calc_el(Vj, Tc, I/*pack.get_akt().IF*/, eml.pack);
        if (is_votma && !eml.pack.get_votma())
            eml.pack.set_votma(true);
//printf("votma2=%u\n", (uns)pack.get_akt().get_votma());
        //        p_junction->el_egyenlet.junction_calc_el_from_U(Vj, Tc, pack.get_akt().IF, pack.get_akt());
        eml.G_dioda = eml.pack.G * Aface_per_Ajunction_1;
        eml.g_diff_dioda = eml.pack.dI_per_dVj * Aface_per_Ajunction_1;
        eml.dI_dioda_per_dT = eml.pack.dI_per_dT * Aface_per_Ajunction_1;
        //P_dissz_P_eso.get_akt() = P_dissz_P_eso.get_elozo();

        // ellenállás paraméterek

        ertek_t fajlagos_ertek_r = p_anyag->elvez.get_value(Tc, p_cella->emlekek.get_akt().H);
        eml.G_ellenallas = fajlagos_ertek_r.ertek * A_per_L;
        eml.g_diff_ellenallas = eml.G_ellenallas;
        eml.dG_ellenallas_per_dT = fajlagos_ertek_r.derivalt * A_per_L;

        // teljes face paraméterei

        eml.G = replusz(eml.G_ellenallas, eml.G_dioda);
        eml.g_diff = replusz(eml.g_diff_ellenallas, eml.g_diff_dioda);
        //printf("%g, %g\n", G.get_akt(), g_diff.get_akt());
    }
    //***********************************************************************
    void update_admittanciamatrix_SA() override {
    //***********************************************************************
        rvt Gelozo = emlekek.get_elozo().G;
        *p_tok->p_Yii = Gelozo;
        *p_tok->p_Yie = -Gelozo;
        *p_tok->p_Yei = -Gelozo;
        *p_cella->p_Yee_dc += Gelozo;
    }
    //***********************************************************************
    void set_cella_is_szimm() {
    //***********************************************************************
        p_cella->set_to_NR_van_csatolt_nemszimm();
    }
    //***********************************************************************
    void update_Jakobi_NR_nincs_csatolt() override {
    //***********************************************************************
        rvt Gdelozo = emlekek.get_elozo().g_diff;
        *p_tok->p_Yii = Gdelozo;
        *p_tok->p_Yie = -Gdelozo;
        *p_tok->p_Yei = -Gdelozo;
        *p_cella->p_Yee_dc += Gdelozo;
    }
    //***********************************************************************
    void update_Jakobi_NR_van_csatolt() override {
    //***********************************************************************
        const emlekek_strukt & eml = emlekek.get_elozo();
        rvt U1 = p_tok->emlekek.get_elozo().UT;
        rvt Ux = eml.U_kozep;
        rvt Uc = p_cella->el_center_face_dc.emlekek.get_elozo().UT;
        rvt Gdelozo = eml.g_diff;
        if(U1 != eml.bU1)printf("%g, %g, %g, %g\n", U1 - eml.bU1, Uc, eml.bUc, p_cella->el_center_face_dc.emlekek.get_akt().UT);
        rvt Vj = U1 - Ux;
        rvt Vr = Ux - Uc;
//        if (p_cella->cella_index == 22493)
//            printf("\nUc=%g, Ux=%g, U1=%g, Vj=%g, Vr=%g\n", Uc, Ux, U1, Vj, Vr);
        rvt egy_per_nevezo = 1 / (eml.g_diff_dioda + eml.g_diff_ellenallas);
//printf("%g\n", Gdelozo);
        *p_tok->p_Yii = Gdelozo;
        *p_tok->p_Yie = -Gdelozo;
        *p_tok->p_Yei = -Gdelozo;
        *p_cella->p_Yee_dc += Gdelozo;
        // valamiért nem eszi meg a tényleges et tagot
        // ha a rendes et tagot használom, szétbassza, ha a -1-szeresét, akkor csak lassítja kissé a konvergenciát
        rvt et = Vr * eml.dG_ellenallas_per_dT; // (g_diff_ellenallas.get_elozo() * dI_dioda_per_dT.get_elozo() + g_diff_dioda.get_elozo() * Vr * dG_ellenallas_per_dT.get_elozo()) * egy_per_nevezo;
        //printf("gr=%g, gj=%g, g=%g, et=%g\n", g_diff_ellenallas, g_diff_dioda, g_diff, et);
        *p_tok->p_Yit = et;
        *p_cella->p_Yet_dc -= et;
        *p_tok->p_Yti = 0;
        // rvt dP_per_dV = 2*G*V;
        // *p_tok->p_Yti = -dP_per_dV; // nem áramfüggõ ellenállásnál -dP/dV = -d(GV^2)/dV = -2GV
        // *p_cella->p_Yte_dc += dP_per_dV;
        // A tt tagot a cella::fw_klaszter_1_update_Jakobi_dc adja hozzá
        //rvt dPd_per_dVj = g_diff_dioda.get_elozo()*Vj*P_dissz_P_eso.get_elozo();
        rvt dPd_per_dVj = eml.g_diff_dioda*Vj;
        rvt dPd_per_dVr = 2 * eml.G_ellenallas * Vr * Vr;
        rvt te = (eml.G_ellenallas * dPd_per_dVj + eml.g_diff_dioda * dPd_per_dVr) * egy_per_nevezo;
        *p_tok->p_Yti = -te;
        *p_cella->p_Yte_dc += te;
    }
    //***********************************************************************
    void agaram_szamitasa() override {
    //***********************************************************************
        emlekek_strukt & eml = emlekek.get_akt();
        rvt Uc = p_cella->el_center_face_dc.emlekek.get_akt().UT;
        rvt U1 = p_tok->emlekek.get_akt().UT;
        rvt Ux = eml.U_kozep;
        rvt V = U1 - Uc; 
//        V = V > 3 ? 3 : V; // !!!! a konvergencia szétszállt, ezért a korlátozás. Ha lehet ennél nagyogg junction fesz, akkor módosítani kell.
        rvt I = eml.G * V;
        p_tok->emlekek.get_akt().IP = I;
        //printf("2 %g\n", p_tok->IP.get_akt());
        //printf("%g, %g, %g\n", I*Ajunction_per_Aface_1, pack.get_akt().IF, abs(I)*Ajunction_per_Aface_1 - pack.get_akt().IF);
        //if (I != 0)
        //    pack.get_akt().IF = abs(I)*Ajunction_per_Aface;
        rvt Tc = p_cella->get_Tc_akt();
        Tc = Tc < 0 ? 0 : Tc;

        // sugárzott teljesítmény

        p_cella->emlekek.get_akt().A_LED = p_cella->emlekek.get_akt().A_LED + A;

        //rvt junction_dF_per_dVj;
        //ertek_t junction_F = p_junction->rad.get_rad_flux(Tc, I*Ajunction_per_Aface, g_diff_dioda.get_akt()*Ajunction_per_Aface, junction_dF_per_dVj);
        p_junction->rad.junction_calc_rad(Tc, eml.pack);
        rvt rad = eml.pack.Fi_e * Aface_per_Ajunction_1;
        //printf("rad=%g, V*I=%g, V=%g, I=%g, Uc=%g, U1=%g, Ux=%g, Tc=%g\n", rad, V*I, V, I, Uc, U1, Ux, Tc);
        if (rad > V*I) {
            rad = V*I;
        }
        p_tok->emlekek.get_akt().rad = rad;
        p_cella->emlekek.get_akt().sum_rad = p_cella->emlekek.get_akt().sum_rad + rad;

        // luminancia

        p_junction->lum.junction_calc_lum(Tc, eml.pack);
        rvt lum = eml.pack.Fi_nu * Aface_per_Ajunction_1;
        p_tok->emlekek.get_akt().lum = lum;
        p_cella->emlekek.get_akt().sum_lum = p_cella->emlekek.get_akt().sum_lum + lum;

        // disszipáció

        *p_cella->p_sum_disszipacio_dc += V * I - rad;
        eml.P_dissz_P_eso = 1 - rad / (V * I + 1e-20);
        //*p_cella->p_sum_ddissz_per_dT_dc += dG_per_dT * V * V;
        rvt Vj = U1 - Ux;
        rvt Vr = Ux - Uc;
        rvt dP_per_dT = eml.dI_dioda_per_dT*Vj + eml.dG_ellenallas_per_dT * Vr * Vr;
        // Ha -= van, ami kijön, akkor szétbassza. Ha +=, akkor csak lassítja a konvergenciát.
        // *p_cella->p_sum_ddissz_per_dT_dc += dP_per_dT;
    }
    //***********************************************************************
    void update_UT_IP_NR() override {
    //***********************************************************************
        rvt dUc0 = *p_cella->el_center_face_dc.p_UTi;
        rvt dUc = dUc0 * akt_sim.alfa;
        rvt dU10 = *p_tok->p_UTi;
        rvt dU1 = dU10 * akt_sim.alfa;
        rvt Uc = p_cella->el_center_face_dc.emlekek.get_kiindulo().UT + dUc;
        rvt U1 = p_tok->emlekek.get_kiindulo().UT + dU1;
        p_tok->emlekek.get_akt().UT = U1;
        const emlekek_strukt & eml = emlekek.get_kiindulo();
        emlekek.get_akt().U_kozep = (U1*eml.G_dioda + Uc*eml.G_ellenallas) / (eml.G_dioda + eml.G_ellenallas);
        //p_tok->IP.get_akt() = p_tok->IP.get_kiindulo() + *p_tok->p_IPi * akt_sim.alfa; // ez hibásan mûködik, nem szabad !!!
//        if (p_cella->cella_index == 22493)
//            printf("\ndUc0=%g, dUc=%g, Uc=%g, dU10=%g, dU1=%g, U1=%g, U_kozep=%g\n", dUc0, dUc, Uc, dU10, dU1, U1, U_kozep.get_akt());
    }
    //***********************************************************************
    void update_UT_IP_SA() override {
    //***********************************************************************
        rvt Uc = *p_cella->el_center_face_dc.p_UTi;
        rvt U1 = *p_tok->p_UTi;
        p_tok->emlekek.get_akt().UT = U1;
        const emlekek_strukt & eml = emlekek.get_kiindulo();
        emlekek.get_akt().U_kozep = (U1*eml.G_dioda + Uc*eml.G_ellenallas) / (eml.G_dioda + eml.G_ellenallas);
        p_tok->emlekek.get_akt().IP = *p_tok->p_IPi; // nem tudom, hogy ennek van-e bármi értelme
    }
/*
    //***********************************************************************
    void print_G()const override { 
        rvt Tc = p_cella->get_Tc_akt();
        rvt Vj = p_tok->UT.get_elozo() - U_kozep.get_elozo();
        rvt fajlagos_diff_vezetes_j = rvt();
        printf("\n");
        printf("Tc_elozo=%g, Tc_akt=%g, Vj_elozo=%g, Vj_akt=%g\n", p_cella->th_center_face_dc.UT.get_elozo(), Tc, Vj, p_tok->UT.get_akt() - U_kozep.get_akt());
        ertek_t fajlagos_ertek_j = p_junction->el_egyenlet.get_G_junction(Tc, Vj, fajlagos_diff_vezetes_j, true);
        printf("Gj=%g, gj=%g\n", G_dioda.get_akt(), g_diff_dioda.get_akt());
        //printf("G=%g, Gj=%g, Gr=%g, g=%g, gj=%g, gr=%g\n", G, G_dioda, G_ellenallas, g_diff, g_diff_dioda, g_diff_ellenallas);
        printf("\n");
    }
*/
};

/* Úgy látom, ezek maradhatnak.
    //***********************************************************************
    void update_hibaaram_NR_kereszt(rvt & sum_el_hiba, rvt & sum_th_hiba) override {
    //***********************************************************************
        rvt IP = p_tok->IP.get_akt();
        *p_tok->p_Ji = IP;     // !! betesszük a teljes áramot az inhomogén áram vektorba, így két cella összevonásakor pont a hibaáram marad itt.
        csatlakozo_aramok_dc[csatlakozo_index] = IP;
        if (is_el) {
            sum_el_hiba -= IP; // ki kell vonni az áramot, mert a külsõ kapocsnak megfelelõ az iránya-
        }
        else {
            sum_th_hiba -= IP; // ki kell vonni az áramot, mert a külsõ kapocsnak megfelelõ az iránya 
        }
    }
    //***********************************************************************
    void debug_write(::std::ofstream & fs, bool is_endl) const override { 
    //***********************************************************************
        if (is_endl) fs << "csatlak  "; 
        os_face_core_dc::debug_write(fs, false); 
        fs << " G=" << ::std::setw(12) << G;
        if (is_endl) fs << ::std::endl;
    }
*/

}

#endif
