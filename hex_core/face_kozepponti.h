//***********************************************************************
// középponti face osztályok header
// Creation date:  2018. 08. 31.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef FACE_KOZEPPONTI_HEADER
#define	FACE_KOZEPPONTI_HEADER
//***********************************************************************


//***********************************************************************
#include "cella_es_face.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
class dc_kozepponti_face : public os_face_core_dc {
//***********************************************************************
protected:
    rvt I_gen, V_per_Vstruct;
    const adat_struktura * p_struktura;
public:
    //***********************************************************************
    dc_kozepponti_face(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :os_face_core_dc(p_tok, p_cella, p_adat_face, is_el), I_gen{ rvt() }, V_per_Vstruct{ rvt() }, p_struktura{ nullptr } {}
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        p_anyag = &akt_sim.p_akt_sim->anyagok[akt_sim.akt_anyag_index[p_cella->p_adatcella->anyag_index]]; // cella anyag !!!
        p_struktura = &akt_sim.akt_strukturak[p_cella->p_adatcella->struktura_index];
        V_per_Vstruct = akt_sim.p_akt_sim->meretek[p_cella->p_adatcella->volume_index] / p_struktura->volume;
    }
    //***********************************************************************
    void update_peremfeltetel() override {}
    //***********************************************************************
    void update_UT_IP_NR() override {
    //***********************************************************************
        //if (*p_tok->p_UTi > rvt(10))
        //    *p_tok->p_UTi = rvt(10);
        //if (*p_tok->p_UTi < rvt(-10))
        //    *p_tok->p_UTi = rvt(-10);
        p_tok->emlekek.get_akt().UT = p_tok->emlekek.get_kiindulo().UT + *p_tok->p_UTi * akt_sim.alfa;
    }
    //***********************************************************************
    void update_UT_IP_SA() override {
    //***********************************************************************
        p_tok->emlekek.get_akt().UT = *p_tok->p_UTi;
    }
    //***********************************************************************
    void update_hibaaram_NR_kereszt(rvt & sum_el_hiba, rvt & sum_th_hiba) override {}
    //***********************************************************************
    virtual void debug_write(::std::ofstream & fs, bool is_endl) const override { if (is_endl) fs << "kozep    "; os_face_core_dc::debug_write(fs, false); fs << " I_gen=" << ::std::setw(12) << I_gen; if (is_endl) fs << ::std::endl; }
    //***********************************************************************
};


//***********************************************************************
class dc_kozepponti_face_szakadas : public dc_kozepponti_face {
// középponti face
//***********************************************************************
public:
    //***********************************************************************
    dc_kozepponti_face_szakadas(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_kozepponti_face(p_tok, p_cella, p_adat_face, is_el) {}
    //***********************************************************************
    face_azonosito get_azon() const { return fa_kd_x_kdn; }
    //***********************************************************************
    bool is_konst() const { return true; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_kozepponti_face::init();
        // G=0, J=0 mindig, úgyhogy semmi tennivaló
    }
    //***********************************************************************
    void update_gerj() override {}
    //***********************************************************************
    void update_uj_lepeshez() override {} // nincs tennivaló
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
        p_tok->emlekek.get_akt().IP = rvt();
    }
    //***********************************************************************
    void update_inhom_SA_kozepponti(rvt dissz) override {
    //***********************************************************************
        if (is_el)
            *p_cella->p_Je_dc = rvt();
        else
            *p_cella->p_Jt_dc = -dissz;
    }
    //***********************************************************************
    void update_hibaaram_NR_kozepponti(rvt & sum_el_hiba, rvt & sum_th_hiba) override {
    //***********************************************************************
        if (is_el) {
            *p_cella->p_Je_dc = sum_el_hiba;
        }
        else {
            *p_cella->p_Jt_dc = sum_th_hiba;
        }
    }
};


//***********************************************************************
class dc_kozepponti_face_IP_generator_konst : public dc_kozepponti_face {
// középponti face, paraméterei az p_cella-ból veszi
//***********************************************************************
public:
    //***********************************************************************
    dc_kozepponti_face_IP_generator_konst(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_kozepponti_face(p_tok, p_cella, p_adat_face, is_el) {}
    //***********************************************************************
    face_azonosito get_azon() const { return fa_kd_i_kdn; }
    //***********************************************************************
    bool is_konst() const { return true; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_kozepponti_face::init();
        update_gerj();
    }
    //***********************************************************************
    void update_gerj() override {
    //***********************************************************************
        if (is_el) {
            igaz_e_hiba(p_struktura->p_el_gerj == nullptr || p_struktura->p_el_gerj->tipus != gt_I, "dc_kozepponti_face_IP_generator_konst::init", "excit type");
            I_gen = p_struktura->p_el_gerj->ertek_tolt * V_per_Vstruct;
        }
        else {
            igaz_e_hiba(p_struktura->p_th_gerj == nullptr || p_struktura->p_th_gerj->tipus != gt_P, "dc_kozepponti_face_IP_generator_konst::init", "excit type");
            I_gen = p_struktura->p_th_gerj->ertek_tolt * V_per_Vstruct;
        }
    }
    //***********************************************************************
    void update_uj_lepeshez() override {} // nincs tennivaló
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
        p_tok->emlekek.get_akt().IP = -I_gen;
    }
    //***********************************************************************
    void update_inhom_SA_kozepponti(rvt dissz) override {
    //***********************************************************************
        if (is_el)
            *p_cella->p_Je_dc = -I_gen;
        else
            *p_cella->p_Jt_dc = -I_gen - dissz;
    }
    //***********************************************************************
    void update_hibaaram_NR_kozepponti(rvt & sum_el_hiba, rvt & sum_th_hiba) override {
    //***********************************************************************
        if (is_el) {
            sum_el_hiba -= I_gen;
            *p_cella->p_Je_dc = sum_el_hiba;
        }
        else {
            sum_th_hiba -= I_gen;
            *p_cella->p_Jt_dc = sum_th_hiba;
        }
    }
};


//***********************************************************************
class dc_kozepponti_face_Cth_IP_generator_konst : public dc_kozepponti_face {
// középponti face, paraméterei az p_cella-ból veszi
// hõkapacitás ill. hõkapacitás + disszipáció gerjesztés
//***********************************************************************
    rvt T_start; // Az elõzõ valódi lépésben kialakult hõmérséklet
    rvt Gcth, inhom;
public:
    //***********************************************************************
    dc_kozepponti_face_Cth_IP_generator_konst(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_kozepponti_face(p_tok, p_cella, p_adat_face, is_el), T_start{ rvt() }, Gcth{ rvt() }, inhom{ rvt() } {}
    //***********************************************************************
    face_azonosito get_azon() const { return fa_kd_ci_ktn; }
    //***********************************************************************
    bool is_konst() const { return true; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        igaz_e_hiba(is_el, "dc_kozepponti_face_Cth_IP_generator_konst::init", "electrical central face with Cth");
        dc_kozepponti_face::init();
        update_gerj();
    }
    //***********************************************************************
    void update_dt() override {
    //***********************************************************************
        // Gcth = C / dt
        rvt fajlagos_Cth = p_anyag->Cth.get_value(p_cella->get_Tc_megtartando(), ertek_t()).ertek; // p_anyag a középponti face-eknél a cella anyagára van állítva
        rvt Cth = rvt(fajlagos_Cth * akt_sim.p_akt_sim->meretek[p_cella->p_adatcella->volume_index] * zaj_szorzo);
        Gcth = Cth / akt_sim.value; // ez a face csak tranziens módban jöhet létre, amikor akt_sim.value nem 0
    }
    //***********************************************************************
    void update_gerj() override {
    //***********************************************************************
        igaz_e_hiba(p_struktura->p_th_gerj != nullptr && p_struktura->p_th_gerj->tipus != gt_P, "dc_kozepponti_face_Cth_IP_generator_konst::init", "excit type");
        igaz_e_hiba(akt_sim.tipus != alt_trans, "dc_kozepponti_face_Cth_IP_generator_konst::init", "this face type can be used in transient mode");
        I_gen = p_struktura->p_th_gerj==nullptr ? rvt() : p_struktura->p_th_gerj->ertek_tolt * V_per_Vstruct;
        update_dt();
    }
    //***********************************************************************
    void update_uj_lepeshez() override { 
    //***********************************************************************
        T_start = p_tok->emlekek.get_megtartando().UT;
        p_cella->konst_Ytt_dc += Gcth;
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
        rvt T = p_tok->emlekek.get_akt().UT;
        inhom = -I_gen - T_start * Gcth;
        p_tok->emlekek.get_akt().IP = inhom + T * Gcth;
    }
    //***********************************************************************
    void update_inhom_SA_kozepponti(rvt dissz) override {
    //***********************************************************************
        *p_cella->p_Jt_dc = inhom - dissz;
    }
    //***********************************************************************
    void update_hibaaram_NR_kozepponti(rvt & sum_el_hiba, rvt & sum_th_hiba) override {
    //***********************************************************************
        sum_th_hiba += p_tok->emlekek.get_akt().IP;
        *p_cella->p_Jt_dc = sum_th_hiba;
    }
    //***********************************************************************
    virtual void debug_write(::std::ofstream & fs, bool is_endl) const override { 
    //***********************************************************************
        if (is_endl) fs << "Cth_cons "; 
        dc_kozepponti_face::debug_write(fs, false);
        fs << " T_start=" << ::std::setw(12) << T_start;
        fs << " Gcth=" << ::std::setw(12) << Gcth;
        if (is_endl) fs << ::std::endl;
    }
    //***********************************************************************
};


//***********************************************************************
class dc_kozepponti_face_Cth_IP_generator_Tfuggo : public dc_kozepponti_face {
// középponti face, paraméterei az p_cella-ból veszi
// hõkapacitás ill. hõkapacitás + disszipáció gerjesztés
//***********************************************************************
    struct emlekek_strukt {
        rvt Gcth, inhom;
        ertek_t CthAB, CBm; // Ha fázisváltó
        rvt dG_per_dT;
    };
    emlek<emlekek_strukt> emlekek;
    rvt T_start; // Az elõzõ valódi lépésben kialakult hõmérséklet
    rvt fajlagos_Cth_szorzoja;
public:
    //***********************************************************************
    dc_kozepponti_face_Cth_IP_generator_Tfuggo(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_kozepponti_face(p_tok, p_cella, p_adat_face, is_el), fajlagos_Cth_szorzoja{ rvt() } {}
    //***********************************************************************
    face_azonosito get_azon() const { return fa_kd_ci_ttn; }
    //***********************************************************************
    bool is_konst() const { return false; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        igaz_e_hiba(is_el, "dc_kozepponti_face_Cth_IP_generator_Tfuggo::init", "electrical central face with Cth");
        igaz_e_hiba(akt_sim.tipus != alt_trans, "dc_kozepponti_face_Cth_IP_generator_Tfuggo::init", "this face type can be used in transient mode");
        dc_kozepponti_face::init();
        if (p_anyag->fazisvaltas.is_fazisvalto()) {
            ertek_t Cth = p_anyag->Cth.get_value(rvt(), p_anyag->fazisvaltas.H(rvt(), rvt()));
            emlekek.get_akt().CthAB = emlekek.get_elozo_to_overwrite().CthAB = emlekek.get_prevprev_to_overwrite().CthAB
                = emlekek.get_megtartando_to_overwrite().CthAB = emlekek.get_mentendo_to_overwrite().CthAB
                = emlekek.get_kiindulo_to_overwrite().CthAB = Cth;
            emlekek.get_akt().CBm = emlekek.get_elozo_to_overwrite().CBm = emlekek.get_prevprev_to_overwrite().CBm
                = emlekek.get_megtartando_to_overwrite().CBm = emlekek.get_mentendo_to_overwrite().CBm
                = emlekek.get_kiindulo_to_overwrite().CBm = Cth;
        }
        update_gerj();
    }
    //***********************************************************************
    void update_dt() override {
    //***********************************************************************
        fajlagos_Cth_szorzoja = rvt(akt_sim.p_akt_sim->meretek[p_cella->p_adatcella->volume_index] / akt_sim.value * zaj_szorzo); // ez a face csak tranziens módban jöhet létre, amikor akt_sim.value nem 0
    }
    //***********************************************************************
    void update_gerj() override {
    //***********************************************************************
        igaz_e_hiba(p_struktura->p_th_gerj != nullptr && p_struktura->p_th_gerj->tipus != gt_P, "dc_kozepponti_face_Cth_IP_generator_Tfuggo::init", "excit type");
        I_gen = p_struktura->p_th_gerj == nullptr ? rvt() : p_struktura->p_th_gerj->ertek_tolt * V_per_Vstruct;
        update_dt();
    }
    //***********************************************************************
    void update_uj_lepeshez() override { T_start = p_tok->emlekek.get_megtartando().UT; }
    //***********************************************************************
    void update_nemlinearis_parameterek() override {
    //***********************************************************************
        // p_anyag a középponti face-eknél a cella anyagára van állítva
        rvt Tc_akt = p_cella->get_Tc_akt();
        ertek_t H_akt = p_cella->emlekek.get_akt().H;
        ertek_t fajlagos_Cth = p_anyag->Cth.get_value(Tc_akt, H_akt);
        emlekek_strukt & eml = emlekek.get_akt();
        if (p_anyag->fazisvaltas.is_fazisvalto()) {
            eml.CthAB = fajlagos_Cth;
            ertek_t akt_CBm = p_anyag->fazisvaltas.CBm(p_cella->get_Tc_megtartando(), Tc_akt
                , p_cella->emlekek.get_megtartando().H.ertek, H_akt.ertek
                , emlekek.get_megtartando().CthAB.ertek, fajlagos_Cth.ertek, emlekek.get_megtartando().CBm);
            eml.CBm = akt_CBm;
            eml.Gcth = akt_CBm.ertek * fajlagos_Cth_szorzoja;
            eml.dG_per_dT = akt_CBm.derivalt * fajlagos_Cth_szorzoja;
        }
        else {
            eml.Gcth = fajlagos_Cth.ertek * fajlagos_Cth_szorzoja;
            eml.dG_per_dT = fajlagos_Cth.derivalt * fajlagos_Cth_szorzoja;
        }
        eml.inhom = emlekek.get_elozo().inhom;
    }
    //***********************************************************************
    void update_admittanciamatrix_SA() override {
    //***********************************************************************
        *p_cella->p_Ytt_dc += emlekek.get_elozo().Gcth;
    }
    //***********************************************************************
    void update_Jakobi_NR_nincs_csatolt() override {
    //***********************************************************************
        *p_cella->p_Ytt_dc += emlekek.get_elozo().Gcth;
        // TODO: elõször rakjuk össze, utána megnézhetjük, hogy a hõkapacitás NR-osítása mûködik-e (a felsõ sort kell cserélni az alsó kettõre)
//        rvt dT = p_tok->UT.get_elozo() - p_cella->th_center_face_dc.UT.get_elozo(); // a teljes esõ hõmérsékletet kell venni
//        *p_cella->p_Ytt_dc += Gcth - dG_per_dT*dT;
    }
    //***********************************************************************
    void update_Jakobi_NR_van_csatolt() override {
    //***********************************************************************
        *p_cella->p_Ytt_dc += emlekek.get_elozo().Gcth;
        // Jelenleg szimmetrikus ettõl még az admittanciamátrix, de ha ez változna, akkor a nemsimm beállítandó
        // p_cella->set_to_nemsimm();
        // TODO: elõször rakjuk össze, utána megnézhetjük, hogy a hõkapacitás NR-osítása mûködik-e (a felsõ sort kell cserélni az alsó kettõre)
//        rvt dT = p_tok->UT.get_elozo() - p_cella->th_center_face_dc.UT.get_elozo(); // a teljes esõ hõmérsékletet kell venni, nem a hibahõmérsékletet!!!
//        *p_cella->p_Ytt_dc += Gcth - dG_per_dT*dT;
    }
    //***********************************************************************
    void agaram_szamitasa() override {
    //***********************************************************************
        rvt T = p_tok->emlekek.get_akt().UT;
        emlekek_strukt & eml = emlekek.get_akt();
        eml.inhom = -I_gen - T_start * eml.Gcth;
        p_tok->emlekek.get_akt().IP = eml.inhom + T * eml.Gcth;
    }
    //***********************************************************************
    void update_inhom_SA_kozepponti(rvt dissz) override {
    //***********************************************************************
        *p_cella->p_Jt_dc = emlekek.get_elozo().inhom - dissz;
    }
    //***********************************************************************
    void update_hibaaram_NR_kozepponti(rvt & sum_el_hiba, rvt & sum_th_hiba) override {
    //***********************************************************************
        sum_th_hiba += p_tok->emlekek.get_akt().IP;
        *p_cella->p_Jt_dc = sum_th_hiba;
    }
    //***********************************************************************
    virtual void debug_write(::std::ofstream & fs, bool is_endl) const override { 
    //***********************************************************************
        if (is_endl) fs << "Cth_fugg "; 
        dc_kozepponti_face::debug_write(fs, false);
        const emlekek_strukt & eml = emlekek.get_akt();
        fs << " T_start=" << ::std::setw(12) << T_start;
        fs << " Gcth=" << ::std::setw(12) << eml.Gcth;
        fs << " dG_per_dT=" << ::std::setw(12) << eml.dG_per_dT;
        fs << " CBm=" << ::std::setw(12) << eml.CBm.ertek;
        fs << " CBmd=" << ::std::setw(12) << eml.CBm.derivalt;
        if (is_endl) fs << ::std::endl;
    }
    //***********************************************************************
};


//***********************************************************************
class dc_kozepponti_face_UT_generator_konst : public dc_kozepponti_face {
// középponti face, paraméterei az p_cella-ból veszi
// fesz v. hõm. gerjesztés
//***********************************************************************
    rvt UT;
public:
    //***********************************************************************
    dc_kozepponti_face_UT_generator_konst(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :dc_kozepponti_face(p_tok, p_cella, p_adat_face, is_el), UT{ rvt() } {}
    //***********************************************************************
    face_azonosito get_azon() const { return fa_kd_u_kdn; }
    //***********************************************************************
    bool is_konst() const { return true; }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        dc_kozepponti_face::init();
        update_gerj();
        // Ycc = g_max;
    }
    //***********************************************************************
    void update_gerj() override {
    //***********************************************************************
        if (is_el) {
            igaz_e_hiba(p_struktura->p_el_gerj == nullptr || p_struktura->p_el_gerj->tipus != gt_U, "dc_kozepponti_face_UT_generator_konst::init", "excit type");
            UT = p_struktura->p_el_gerj->ertek_tolt;
            I_gen = UT * rvt(g_max);
        }
        else {
            igaz_e_hiba(p_struktura->p_th_gerj == nullptr || p_struktura->p_th_gerj->tipus != gt_T, "dc_kozepponti_face_UT_generator_konst::init", "excit type");
            UT = p_struktura->p_th_gerj->ertek_tolt;
            I_gen = UT * rvt(g_max);
        }
    }
    //***********************************************************************
    void update_uj_lepeshez() override {
    // kapacitások frissítése, konstansok admittanciamátrixokba írása, 
    // konstans középponti admittanciák növelése
    // a perem csak a középpontit változtatja
    //***********************************************************************
        if (is_el) {
            p_cella->konst_Yee_dc += rvt(g_max);
        }
        else {
            p_cella->konst_Ytt_dc += rvt(g_max);
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
        p_tok->emlekek.get_akt().IP = -I_gen;
        // p_tok->UT = UT;
    }
    //***********************************************************************
    void update_inhom_SA_kozepponti(rvt dissz) override {
    //***********************************************************************
        if (is_el)
            *p_cella->p_Je_dc = -I_gen;
        else
            *p_cella->p_Jt_dc = -I_gen - dissz;
    }
    //***********************************************************************
    void update_hibaaram_NR_kozepponti(rvt & sum_el_hiba, rvt & sum_th_hiba) override {
    //***********************************************************************
        rvt dI = (UT - p_tok->emlekek.get_akt().UT) * rvt(g_max);
        if (is_el) {
            sum_el_hiba -= dI; // I_gen;
            *p_cella->p_Je_dc = sum_el_hiba;
            //sum_el_hiba = rvt();
        }
        else {
            sum_th_hiba -= dI; // I_gen;
            *p_cella->p_Jt_dc = sum_th_hiba;
            //sum_th_hiba = rvt();
        }
    }
};


    

}

#endif
