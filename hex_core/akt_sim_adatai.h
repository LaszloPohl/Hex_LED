//***********************************************************************
// aktu�lis szimul�ci� adatai oszt�ly header
// Creation date:  2018. 08. 10.
// Creator:        Pohl L�szl�
//***********************************************************************


//***********************************************************************
#ifndef AKT_SIM_ADATAI_HEADER
#define	AKT_SIM_ADATAI_HEADER
//***********************************************************************


//***********************************************************************
#include "bemenet.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
class vezerlo;
//***********************************************************************


//***********************************************************************
class akt_sim_adatai {
//***********************************************************************
    //***********************************************************************
    void clear() {
    //***********************************************************************
        tipus = alt_dc;
        value = I0 = max_hiba = Tamb = rvt();
        alfa = rvt(1);
        max_iter = 0;
        is_uj_fa_kell = is_uj_cellaszerkezet_kell = is_uj_cellaszerkezet_elso_kor = is_uj_facek_letrehozasa_kell = is_tamb_update_kell = false;
        is_gerj_update_kell = is_del_all_prev = is_del_all_fa = is_peremfelt_update_kell = is_force_face_update = false;
        is_subiter_dt_changed = is_iter_dt_csokkento = false;
        is_ignore_error = false;
        fa_adat = fat_double;
        solver_type = st_sunred;
        p_akt_sim = nullptr;
        p_akt_eredm = nullptr;
        akt_perem_index.clear();
        akt_anyag_index.clear();
        akt_junct_index.clear();
        akt_perem.clear();
        akt_strukturak.clear();
    }
    //***********************************************************************
    void clear_nem_oroklodok() {
    //***********************************************************************
        tipus = alt_dc;
        value = rvt();
        alfa = rvt(1);
        is_uj_fa_kell = is_uj_cellaszerkezet_kell = is_uj_cellaszerkezet_elso_kor = is_uj_facek_letrehozasa_kell = is_tamb_update_kell = false;
        is_gerj_update_kell = is_del_all_prev = is_del_all_fa = is_peremfelt_update_kell = is_force_face_update = false;
        is_subiter_dt_changed = is_iter_dt_csokkento = false;
        is_ignore_error = false;
    }
    //***********************************************************************
    void clear_for_uj_iteracio() {
    //***********************************************************************
        is_uj_fa_kell = is_uj_cellaszerkezet_kell = is_uj_cellaszerkezet_elso_kor = is_uj_facek_letrehozasa_kell = is_tamb_update_kell = false;
        is_gerj_update_kell = is_del_all_prev = is_del_all_fa = is_peremfelt_update_kell = is_force_face_update = false;
        is_subiter_dt_changed = is_iter_dt_csokkento = false;
    }
    //***********************************************************************
    void reset_akt_perem_index() {
    //***********************************************************************
        akt_perem_index.set_size(p_akt_sim->perem.size());
        for (uns i = 0; i < akt_perem_index.size(); i++)
            akt_perem_index[i] = i;
    }
    //***********************************************************************
    void reset_akt_perem() {
    //***********************************************************************
        akt_perem = p_akt_sim->perem;
    }
    //***********************************************************************
    void reset_akt_strukturak() {
    //***********************************************************************
        akt_strukturak = p_akt_sim->strukturak;
    }
    //***********************************************************************
    void reset_akt_anyag_index() {
    //***********************************************************************
        akt_anyag_index.set_size(p_akt_sim->anyagok.size());
        for (uns i = 0; i < akt_anyag_index.size(); i++)
            akt_anyag_index[i] = i;
    }
    //***********************************************************************
    void reset_akt_junct_index() {
    //***********************************************************************
        akt_junct_index.set_size(p_akt_sim->junctionok.size());
        for (uns i = 0; i < akt_junct_index.size(); i++)
            akt_junct_index[i] = i;
    }
    //***********************************************************************
    void del_all_gerj() {
    //***********************************************************************
        is_gerj_update_kell = true;
        for (uns i = 0; i < akt_strukturak.size(); i++)
            akt_strukturak[i].p_el_gerj = akt_strukturak[i].p_th_gerj = nullptr;
    }
public:
    //***********************************************************************
    analizis_lepes_tipus tipus;             // alt_dc, alt_trans, alt_ac
    rvt value;                              // tranziens l�p�sk�z / AC frekvencia
    rvt I0, max_hiba, Tamb;
    rvt alfa;                               // A vez�rl� �ll�tja be a fesz�lts�gek �jrasz�mol�sa el�tt, default 1.0
    uns max_iter;
    solver_tipus solver_type;
    bool is_ignore_error;                   //
    bool is_uj_fa_kell;                     // vez�rl�, cella �s redukcios_fa figyeli
    bool is_uj_cellaszerkezet_kell;         // vez�rl� �s cella figyeli
    bool is_uj_cellaszerkezet_elso_kor;     // vez�rl�� figyeli, uj_cellaszerkezet eset�n k�tszer kell futtatni a pre_cellafeldolgoz�st, els� k�rben csak foglalunk. (Sug�rk�vet�s elsz�llt, mert nem l�tez� cell�kat c�mezt�nk.)
    bool is_uj_facek_letrehozasa_kell;      // cella figyeli
    bool is_tamb_update_kell;               // vezerlo �s cella figyeli, cell�n�l update-elni kell mindent �s is_del_all_prev is legyen
    bool is_gerj_update_kell;               // cella figyeli, a gerjeszt�st tartalmaz� face-t kell �jragy�rtani
    bool is_peremfelt_update_kell;          // cella figyeli, a peremfelt�tel �rt�ke v�ltozik, friss�teni kell a peremre kapcsol�d� face-eket
    bool is_del_all_prev;                   // cella figyeli
    bool is_del_all_fa;                     // vezerlo figyeli
    bool is_force_face_update;              // cella figyeli, akkor true, ha v�ltozik a dt tranziens szimul�ci�n�l
    bool is_subiter_dt_changed;             // cella figyeli, ha iter�ci� k�zben v�ltozik a dt tranziens szimul�ci�n�l
    bool is_iter_dt_csokkento;              // cella figyeli, ha id�ben visszal�p�nk tranziens szimul�ci�n�l
    fa_adat_tipus fa_adat;                  // fa_adat_tipus, fat_double, fat_float, fat_double_ended_float
    adat_szimulacio * p_akt_sim;
    const vektor<adat_eredm> * p_akt_eredm; // const vektor<adat_eredm> *, mik a ki�rand� eredm�nyek.
    //const vezerlo * p_vezerlo;
    vektor<uns> akt_perem_index;            // alapb�l [i] = i
    vektor<uns> akt_anyag_index;            // alapb�l [i] = i
    vektor<uns> akt_junct_index;            // alapb�l [i] = i
    vektor<adat_peremfeltetel> akt_perem;   // a peremfelt�telek m�solata, ebben lehet m�dos�tani
    vektor<adat_struktura> akt_strukturak;  // a strukt�r�k m�solata, ebben lehet m�dos�tani a gerjeszt�sek c�m�t (m�solat az�rt kell, mert p_akt_sim konstans pointer)
    //***********************************************************************


    //***********************************************************************
    akt_sim_adatai() /*: p_vezerlo{ nullptr }*/ { clear(); }
    //***********************************************************************
    void set_akt_anal_step(adat_szimulacio * p_uj_sim, const adat_analizis_lepes & uj_lepes); // p_uj_sim = nullptr, ha most nem v�ltozik
    void set_akt_iter() { clear_for_uj_iteracio(); } // ha ez nem az els� iter�ci�, akkor h�vhat� �s h�vand�
    //***********************************************************************
    void set_Tamb() {
    //***********************************************************************
        for (uns i = 1; i < akt_perem.size(); i++) {
            akt_perem[i].set_Tamb(Tamb);
        }
    }
    //***********************************************************************
};


}

#endif
