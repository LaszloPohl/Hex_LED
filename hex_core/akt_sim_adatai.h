//***********************************************************************
// aktuális szimuláció adatai osztály header
// Creation date:  2018. 08. 10.
// Creator:        Pohl László
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
    rvt value;                              // tranziens lépésköz / AC frekvencia
    rvt I0, max_hiba, Tamb;
    rvt alfa;                               // A vezérlõ állítja be a feszültségek újraszámolása elõtt, default 1.0
    uns max_iter;
    solver_tipus solver_type;
    bool is_ignore_error;                   //
    bool is_uj_fa_kell;                     // vezérlõ, cella és redukcios_fa figyeli
    bool is_uj_cellaszerkezet_kell;         // vezérlõ és cella figyeli
    bool is_uj_cellaszerkezet_elso_kor;     // vezérlõõ figyeli, uj_cellaszerkezet esetén kétszer kell futtatni a pre_cellafeldolgozást, elsõ körben csak foglalunk. (Sugárkövetés elszállt, mert nem létezõ cellákat címeztünk.)
    bool is_uj_facek_letrehozasa_kell;      // cella figyeli
    bool is_tamb_update_kell;               // vezerlo és cella figyeli, cellánál update-elni kell mindent és is_del_all_prev is legyen
    bool is_gerj_update_kell;               // cella figyeli, a gerjesztést tartalmazó face-t kell újragyártani
    bool is_peremfelt_update_kell;          // cella figyeli, a peremfeltétel értéke változik, frissíteni kell a peremre kapcsolódó face-eket
    bool is_del_all_prev;                   // cella figyeli
    bool is_del_all_fa;                     // vezerlo figyeli
    bool is_force_face_update;              // cella figyeli, akkor true, ha változik a dt tranziens szimulációnál
    bool is_subiter_dt_changed;             // cella figyeli, ha iteráció közben változik a dt tranziens szimulációnál
    bool is_iter_dt_csokkento;              // cella figyeli, ha idõben visszalépünk tranziens szimulációnál
    fa_adat_tipus fa_adat;                  // fa_adat_tipus, fat_double, fat_float, fat_double_ended_float
    adat_szimulacio * p_akt_sim;
    const vektor<adat_eredm> * p_akt_eredm; // const vektor<adat_eredm> *, mik a kiírandó eredmények.
    //const vezerlo * p_vezerlo;
    vektor<uns> akt_perem_index;            // alapból [i] = i
    vektor<uns> akt_anyag_index;            // alapból [i] = i
    vektor<uns> akt_junct_index;            // alapból [i] = i
    vektor<adat_peremfeltetel> akt_perem;   // a peremfeltételek másolata, ebben lehet módosítani
    vektor<adat_struktura> akt_strukturak;  // a struktúrák másolata, ebben lehet módosítani a gerjesztések címét (másolat azért kell, mert p_akt_sim konstans pointer)
    //***********************************************************************


    //***********************************************************************
    akt_sim_adatai() /*: p_vezerlo{ nullptr }*/ { clear(); }
    //***********************************************************************
    void set_akt_anal_step(adat_szimulacio * p_uj_sim, const adat_analizis_lepes & uj_lepes); // p_uj_sim = nullptr, ha most nem változik
    void set_akt_iter() { clear_for_uj_iteracio(); } // ha ez nem az elsõ iteráció, akkor hívható és hívandó
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
