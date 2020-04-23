//***********************************************************************
// perem face osztályok header
// Creation date:  2018. 08. 31.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef FACE_PEREM_HEADER
#define	FACE_PEREM_HEADER
//***********************************************************************


//***********************************************************************
#include "cella_es_face.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
class dc_peremface : public dc_keresztface {
//***********************************************************************
protected:
    const adat_peremfeltetel * p_perem;
public:
    //***********************************************************************
    dc_peremface(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_keresztface(p_tok, p_cella, p_adat_face, is_el), p_perem{ nullptr } {}
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_keresztface::init();
        p_perem = &akt_sim.p_akt_sim->perem[akt_sim.akt_perem_index[p_adat_face->perem_index]];
        // az open perem face ezt nem hívja meg, úgyhogy ha valami miatt kéne, ne felejtsd el
    }
    //***********************************************************************

    //***********************************************************************
    void update_UT_IP_NR() override {} // nincs mit frissíteni, a peremfeszt/T-t az agaram_szamitasa fogja számolni
    //***********************************************************************
    void update_UT_IP_SA() override {} // nincs mit frissíteni, a peremfeszt/T-t az agaram_szamitasa fogja számolni
    //***********************************************************************
    void update_hibaaram_NR_kereszt(rvt & sum_el_hiba, rvt & sum_th_hiba) override { // Specializálni kell, ha módosítja az inhomogén áramot!!
    //***********************************************************************
        rvt IP = p_tok->emlekek.get_akt().IP;
        if (is_el) {
            sum_el_hiba -= IP; // ki kell vonni az áramot, mert a külsõ kapocsnak megfelelõ az iránya-
        }
        else {
            sum_th_hiba -= IP; // ki kell vonni az áramot, mert a külsõ kapocsnak megfelelõ az iránya 
        }
    }
    //***********************************************************************
    virtual void debug_write(::std::ofstream & fs, bool is_endl) const override { if (is_endl) fs << "perem    "; dc_keresztface::debug_write(fs, false); if (is_endl) fs << ::std::endl; }
    //***********************************************************************
};


//***********************************************************************
class dc_perem_HTC_face_ellenallas_elth_nem_disszipalo_konst : public dc_peremface {
// keresztface, paramétereit a p_adat_face-bõl veszi
//***********************************************************************
    rvt Ghtc, Ycc, kulso_arany, G;
public:
    //***********************************************************************
    dc_perem_HTC_face_ellenallas_elth_nem_disszipalo_konst(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_peremface(p_tok, p_cella, p_adat_face, is_el), Ghtc{ rvt() }, Ycc{ rvt() }, kulso_arany{ rvt() }, G{ rvt{ g_min } } {}
    //***********************************************************************
    face_azonosito get_azon() const { return fa_rd_r_kdn; }
    //***********************************************************************
    bool is_konst() const { return true; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_peremface::init();
        rvt fajlagos_vez = is_el 
            ? p_anyag->elvez.get_value(p_cella->get_Tc_megtartando(), ertek_t()).ertek
            : p_anyag->thvez.get_value(p_cella->get_Tc_megtartando(), ertek_t()).ertek;
        G = fajlagos_vez * A_per_L;
        update_peremfeltetel();
    }
    //***********************************************************************
    void update_peremfeltetel() override {
    //***********************************************************************
        Ghtc = p_perem->htc * A;
        Ycc = replusz(G, Ghtc);
        kulso_arany = G / (G + Ghtc);
    }
    //***********************************************************************
    void update_uj_lepeshez() override {
    // kapacitások frissítése, konstansok admittanciamátrixokba írása, 
    // konstans középponti admittanciák növelése
    // a perem csak a középpontit változtatja
    //***********************************************************************
        if (is_el) {
            p_cella->konst_Yee_dc += Ycc;
        }
        else {
            p_cella->konst_Ytt_dc += Ycc;
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
        rvt UT = UTc * kulso_arany;
        rvt V = UT - UTc;
        p_tok->emlekek.get_akt().UT = UT;
        p_tok->emlekek.get_akt().IP = G * V;
    }
};


//***********************************************************************
class dc_perem_HTC_face_ellenallas_elth_nem_disszipalo_Tfuggo : public dc_peremface {
// keresztface, paramétereit a p_adat_face-bõl veszi
//***********************************************************************
    struct emlekek_strukt {
        rvt Ycc, dG_per_dT, G;
    };
    emlek<emlekek_strukt> emlekek;
    rvt Ghtc;
public:
    //***********************************************************************
    dc_perem_HTC_face_ellenallas_elth_nem_disszipalo_Tfuggo(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_peremface(p_tok, p_cella, p_adat_face, is_el), Ghtc{ rvt() } { emlekek.get_akt().G = rvt(g_min); }
    //***********************************************************************
    face_azonosito get_azon() const { return fa_rd_r_tdn; }
    //***********************************************************************
    bool is_konst() const { return false; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_peremface::init();
        ertek_t fajlagos_ertek = is_el
            ? p_anyag->elvez.get_value(p_cella->get_Tc_megtartando(), p_cella->emlekek.get_megtartando().H)
            : p_anyag->thvez.get_value(p_cella->get_Tc_megtartando(), p_cella->emlekek.get_megtartando().H);
        emlekek_strukt & eml = emlekek.get_akt();
        eml.G = fajlagos_ertek.ertek * A_per_L;
        eml.dG_per_dT = fajlagos_ertek.derivalt * A_per_L;
        update_peremfeltetel();
    }
    //***********************************************************************
    void update_peremfeltetel() override {
    //***********************************************************************
        emlekek_strukt & eml = emlekek.get_akt();
        Ghtc = p_perem->htc * A;
        eml.Ycc = replusz(eml.G, Ghtc);
    }
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
        eml.Ycc = replusz(eml.G, Ghtc);
    }
    //***********************************************************************
    void update_admittanciamatrix_SA() override {
    //***********************************************************************
        const emlekek_strukt & eml = emlekek.get_elozo();
        if (is_el) {
            *p_cella->p_Yee_dc += eml.Ycc;
        }
        else {
            *p_cella->p_Ytt_dc += eml.Ycc;
        }
    }
    //***********************************************************************
    void update_Jakobi_NR_nincs_csatolt() override {
    //***********************************************************************
        const emlekek_strukt & eml = emlekek.get_elozo();
        if (is_el) {
            *p_cella->p_Yee_dc += eml.Ycc;
        }
        else {
            rvt Tc = p_cella->th_center_face_dc.emlekek.get_elozo().UT; // a teljes esõ hõmérsékletet kell venni
            rvt arany = Ghtc / (Ghtc + eml.G);
            *p_cella->p_Ytt_dc += eml.Ycc + eml.dG_per_dT*Tc*arany*arany;
        }
    }
    //***********************************************************************
    void update_Jakobi_NR_van_csatolt() override {
    //***********************************************************************
        const emlekek_strukt & eml = emlekek.get_elozo();
        if (is_el) { // nincs disszipáció
            *p_cella->p_Yee_dc += eml.Ycc;
        }
        else {
            rvt Tc = p_cella->th_center_face_dc.emlekek.get_elozo().UT; // a teljes esõ hõmérsékletet kell venni
            rvt arany = Ghtc / (Ghtc + eml.G);
            *p_cella->p_Ytt_dc += eml.Ycc + eml.dG_per_dT*Tc*arany*arany;
        }
    }
    //***********************************************************************
    void agaram_szamitasa() override {
    //***********************************************************************
        rvt UTc = is_el ? p_cella->el_center_face_dc.emlekek.get_akt().UT : p_cella->th_center_face_dc.emlekek.get_akt().UT;
        const emlekek_strukt & eml = emlekek.get_akt();
        rvt UT = UTc * eml.G / (eml.G + Ghtc);
        rvt V = UT - UTc;
        p_tok->emlekek.get_akt().UT = UT;
        p_tok->emlekek.get_akt().IP = eml.G * V;
    }
};


//***********************************************************************
class dc_perem_UT_face_ellenallas_elth_nem_disszipalo_konst : public dc_peremface {
// keresztface, paramétereit a p_adat_face-bõl veszi
//***********************************************************************
    rvt UT, G;
public:
    //***********************************************************************
    dc_perem_UT_face_ellenallas_elth_nem_disszipalo_konst(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_peremface(p_tok, p_cella, p_adat_face, is_el), UT{ rvt() }, G{ rvt{ g_min } } {}
    //***********************************************************************
    face_azonosito get_azon() const { return fa_ud_r_kdn; }
    //***********************************************************************
    bool is_konst() const { return true; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_peremface::init();
        rvt fajlagos_vez = is_el
            ? p_anyag->elvez.get_value(p_cella->get_Tc_megtartando(), ertek_t()).ertek
            : p_anyag->thvez.get_value(p_cella->get_Tc_megtartando(), ertek_t()).ertek;
        G = fajlagos_vez * A_per_L;
        update_peremfeltetel();
    }
    //***********************************************************************
    void update_peremfeltetel() override {
    //***********************************************************************
        UT = p_perem->UT_tolt;
        p_tok->emlekek.get_akt().UT = UT;
    }
    //***********************************************************************
    void update_uj_lepeshez() override {
    // kapacitások frissítése, konstansok admittanciamátrixokba írása, 
    // konstans középponti admittanciák növelése
    // a perem csak a középpontit változtatja
    //***********************************************************************
        if (is_el) {
            p_cella->konst_Yee_dc += G;
        }
        else {
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
        rvt V = UT - UTc;
        p_tok->emlekek.get_akt().UT = UT;
        p_tok->emlekek.get_akt().IP = G * V;
    }
};


//***********************************************************************
class dc_perem_UT_face_ellenallas_elth_nem_disszipalo_Tfuggo : public dc_peremface {
// keresztface, paramétereit a p_adat_face-bõl veszi
//***********************************************************************
    struct emlekek_strukt {
        rvt dG_per_dT, G;
    };
    emlek<emlekek_strukt> emlekek;
    rvt UT;
public:
    //***********************************************************************
    dc_perem_UT_face_ellenallas_elth_nem_disszipalo_Tfuggo(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_peremface(p_tok, p_cella, p_adat_face, is_el), UT{ rvt() } { emlekek.get_akt().G = rvt(g_min); }
    //***********************************************************************
    face_azonosito get_azon() const { return fa_ud_r_tdn; }
    //***********************************************************************
    bool is_konst() const { return false; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_peremface::init();
        ertek_t fajlagos_ertek = is_el
            ? p_anyag->elvez.get_value(p_cella->get_Tc_megtartando(), p_cella->emlekek.get_megtartando().H)
            : p_anyag->thvez.get_value(p_cella->get_Tc_megtartando(), p_cella->emlekek.get_megtartando().H);
        emlekek_strukt & eml = emlekek.get_akt();
        eml.G = fajlagos_ertek.ertek * A_per_L;
        eml.dG_per_dT = fajlagos_ertek.derivalt * A_per_L;
        update_peremfeltetel();
    }
    //***********************************************************************
    void update_peremfeltetel() override {
    //***********************************************************************
        UT = p_perem->UT_tolt;
        p_tok->emlekek.get_akt().UT = UT;
    }
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
        const emlekek_strukt & eml = emlekek.get_elozo();
        if (is_el) {
            *p_cella->p_Yee_dc += eml.G;
        }
        else {
            *p_cella->p_Ytt_dc += eml.G;
        }
    }
    //***********************************************************************
    void update_Jakobi_NR_nincs_csatolt() override {
    //***********************************************************************
        const emlekek_strukt & eml = emlekek.get_elozo();
        if (is_el) {
            *p_cella->p_Yee_dc += eml.G;
        }
        else {
            rvt dT = p_tok->emlekek.get_elozo().UT - p_cella->th_center_face_dc.emlekek.get_elozo().UT; // a teljes esõ hõmérsékletet kell venni
            *p_cella->p_Ytt_dc += eml.G - eml.dG_per_dT*dT;
        }
    }
    //***********************************************************************
    void update_Jakobi_NR_van_csatolt() override {
    //***********************************************************************
        const emlekek_strukt & eml = emlekek.get_elozo();
        if (is_el) { // nincs disszipáció
            *p_cella->p_Yee_dc += eml.G;
        }
        else {
            rvt dT = p_tok->emlekek.get_elozo().UT - p_cella->th_center_face_dc.emlekek.get_elozo().UT; // a teljes esõ hõmérsékletet kell venni, nem a hibahõmérsékletet!!!
            *p_cella->p_Ytt_dc += eml.G - eml.dG_per_dT*dT;
        }
    }
    //***********************************************************************
    void agaram_szamitasa() override {
    //***********************************************************************
        rvt UTc = is_el ? p_cella->el_center_face_dc.emlekek.get_akt().UT : p_cella->th_center_face_dc.emlekek.get_akt().UT;
        rvt V = UT - UTc;
        p_tok->emlekek.get_akt().UT = UT;
        p_tok->emlekek.get_akt().IP = emlekek.get_akt().G * V;
    }
};


//***********************************************************************
class dc_perem_U_face_ellenallas_el_disszipalo_konst : public dc_peremface {
// keresztface, paramétereit a p_adat_face-bõl veszi
//***********************************************************************
    rvt U, G;
public:
    //***********************************************************************
    dc_perem_U_face_ellenallas_el_disszipalo_konst(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_peremface(p_tok, p_cella, p_adat_face, is_el), U{ rvt() }, G{ rvt{ g_min } } {}
    //***********************************************************************
    face_azonosito get_azon() const { return fa_ud_r_ked; }
    //***********************************************************************
    bool is_konst() const { return true; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_peremface::init();
        rvt fajlagos_vez = p_anyag->elvez.get_value(p_cella->get_Tc_megtartando(), ertek_t()).ertek;
        G = fajlagos_vez * A_per_L;
        update_peremfeltetel();
    }
    //***********************************************************************
    void update_peremfeltetel() override {
    //***********************************************************************
        U = p_perem->UT_tolt;
        p_tok->emlekek.get_akt().UT = U;
    }
    //***********************************************************************
    void update_uj_lepeshez() override {
    // kapacitások frissítése, konstansok admittanciamátrixokba írása, 
    // konstans középponti admittanciák növelése
    // a perem csak a középpontit változtatja
    //***********************************************************************
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
    void update_Jakobi_NR_van_csatolt() override { // csak a disszipáció nem konstans
    //***********************************************************************
        rvt V = p_tok->emlekek.get_elozo().UT - p_cella->el_center_face_dc.emlekek.get_elozo().UT; // a teljes esõ feszültséget kell venni, nem a hibafeszültséget!!!
        *p_cella->p_Yte_dc += 2*G*V;
    }
    //***********************************************************************
    void agaram_szamitasa() override {
    //***********************************************************************
        rvt Uc = p_cella->el_center_face_dc.emlekek.get_akt().UT;
        rvt V = U - Uc;
        rvt I = G * V;
        p_tok->emlekek.get_akt().UT = U;
        p_tok->emlekek.get_akt().IP = I;
        *p_cella->p_sum_disszipacio_dc += V * I;
    }
};


//***********************************************************************
class dc_perem_U_face_ellenallas_el_disszipalo_Tfuggo : public dc_peremface {
// keresztface, paramétereit a p_adat_face-bõl veszi
//***********************************************************************
    struct emlekek_strukt {
        rvt dG_per_dT, G;
    };
    emlek<emlekek_strukt> emlekek;
    rvt U;
public:
    //***********************************************************************
    dc_perem_U_face_ellenallas_el_disszipalo_Tfuggo(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_peremface(p_tok, p_cella, p_adat_face, is_el), U{ rvt() } { emlekek.get_akt().G = rvt(g_min); }
    //***********************************************************************
    face_azonosito get_azon() const { return fa_ud_r_ted; }
    //***********************************************************************
    bool is_konst() const { return false; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_peremface::init();
        ertek_t fajlagos_ertek = p_anyag->elvez.get_value(p_cella->get_Tc_megtartando(), ertek_t());
        emlekek_strukt & eml = emlekek.get_akt();
        eml.G = fajlagos_ertek.ertek * A_per_L;
        eml.dG_per_dT = fajlagos_ertek.derivalt * A_per_L;
        update_peremfeltetel();
    }
    //***********************************************************************
    void update_peremfeltetel() override {
    //***********************************************************************
        U = p_perem->UT_tolt;
        p_tok->emlekek.get_akt().UT = U;
    }
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
        *p_cella->p_Yee_dc += emlekek.get_elozo().G;
    }
    //***********************************************************************
    void set_cella_is_szimm() {
    //***********************************************************************
        p_cella->set_to_NR_van_csatolt_nemszimm();
    }
    //***********************************************************************
    void update_Jakobi_NR_nincs_csatolt() override {
    //***********************************************************************
        *p_cella->p_Yee_dc += emlekek.get_elozo().G;
    }
    //***********************************************************************
    void update_Jakobi_NR_van_csatolt() override {
    //***********************************************************************
        rvt dU = p_tok->emlekek.get_elozo().UT - p_cella->el_center_face_dc.emlekek.get_elozo().UT; // a teljes esõ feszültséget kell venni, nem a hibafeszültséget!!!
        const emlekek_strukt & eml = emlekek.get_elozo();
        *p_cella->p_Yee_dc += eml.G;
        *p_cella->p_Yet_dc -= eml.dG_per_dT*dU;
        *p_cella->p_Yte_dc += 2* eml.G*dU;
        // A tt tagot a cella::fw_klaszter_1_update_Jakobi_dc adja hozzá
    }
    //***********************************************************************
    void agaram_szamitasa() override {
    //***********************************************************************
        const emlekek_strukt & eml = emlekek.get_akt();
        rvt Uc = p_cella->el_center_face_dc.emlekek.get_akt().UT;
        rvt V = U - Uc;
        rvt I = eml.G * V;
        p_tok->emlekek.get_akt().UT = U;
        p_tok->emlekek.get_akt().IP = I;
        *p_cella->p_sum_disszipacio_dc += V * I;
        *p_cella->p_sum_ddissz_per_dT_dc += eml.dG_per_dT * V * V;
    }
};


//***********************************************************************
class dc_perem_open_face_ellenallas_elth : public dc_peremface {
// keresztface, paramétereit a p_adat_face-bõl veszi
// open esetén mindegy, hogy lineáris vagy nemlineáris lenne, úgysem folyik áram, és így nem is disszipál
//***********************************************************************
public:
    //***********************************************************************
    dc_perem_open_face_ellenallas_elth(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_peremface(p_tok, p_cella, p_adat_face, is_el) {}
    //***********************************************************************
    face_azonosito get_azon() const { return fa_od_r_kdn; }
    //***********************************************************************
    bool is_konst() const { return true; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        // nem kell a peremface init, úgyhogy azt nem hívjuk
    }
    //***********************************************************************
    void update_peremfeltetel() override {}
    //***********************************************************************
    void update_uj_lepeshez() override {} // nem változtatja az adm.mx-ot
    //***********************************************************************
    void update_nemlinearis_parameterek() override {} // 0
    //***********************************************************************
    void update_admittanciamatrix_SA() override {} // 0
    //***********************************************************************
    void update_Jakobi_NR_nincs_csatolt() override {} // 0
    //***********************************************************************
    void update_Jakobi_NR_van_csatolt() override {} // 0
    //***********************************************************************
    void agaram_szamitasa() override {
    //***********************************************************************
        rvt UTc = is_el ? p_cella->el_center_face_dc.emlekek.get_akt().UT : p_cella->th_center_face_dc.emlekek.get_akt().UT;
        p_tok->emlekek.get_akt().UT = UTc;
        p_tok->emlekek.get_akt().IP = rvt();
    }
};


}

#endif
