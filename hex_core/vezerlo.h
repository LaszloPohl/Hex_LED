//***********************************************************************
// szimul�ci� f� vez�rl� oszt�lya header
// Creation date:  2018. 08. 04.
// Creator:        Pohl L�szl�
//***********************************************************************


//***********************************************************************
#ifndef VEZERLO_HEADER
#define	VEZERLO_HEADER
//***********************************************************************


//***********************************************************************
#include "akt_sim_adatai.h"
#include "redukcios_fa.hpp"
#include "Parhuzamos_Feldolgozo.h"
#include "cella_es_face.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
// *******************   Glob�lis v�ltoz�k    ***************************
extern bemenet bemeneti_adatok;
extern akt_sim_adatai akt_sim, prev_sim;
extern redukcios_fa<float, false>  dcfa_f;
extern redukcios_fa<double, false> dcfa_d;
extern redukcios_fa<::std::complex<float>, true>  acfa_f;
extern redukcios_fa<::std::complex<double>, true> acfa_d;
extern vektor<os_cella*> cellak; // a 0 index� dummy!
extern vektor<rvt> csatlakozo_aramok_dc; // a 0 index� dummy!
extern vektor<iter_csomopont> iter_csomopontok_dc; // a 0 index� dummy!
extern iter_solver<float> dc_itsolver_f;
extern iter_solver<double> dc_itsolver_d;
extern Para_engine feldolgozo;
//***********************************************************************
//***********************************************************************


//***********************************************************************
class vezerlo{
//***********************************************************************

    //***********************************************************************
    bool is_bemeneti_adatok_inicializalva;
    //***********************************************************************

    //***********************************************************************
    vezerlo() :is_bemeneti_adatok_inicializalva{ false } {}
    //***********************************************************************
    vezerlo(const vezerlo&) = delete;
    //***********************************************************************
    vezerlo & operator=(const vezerlo &) = delete;
    //***********************************************************************
    void set_ossz_lepesszam(); // megsz�molja, hogy h�ny anal�zis l�p�s van �sszesen, �s be�ll�tja az �sszl�p�ssz�mot
    //***********************************************************************
    void update_Tamb() {
    //***********************************************************************
        bemeneti_adatok.set_Tamb(akt_sim.Tamb);
        akt_sim.set_Tamb();
    }
    //***********************************************************************
    void run_step_settings();
    void run_fw_threads(uns akt_fa, fa_adat_tipus akt_fa_tipus);
    void run_apply_externals(uns gyoker_index);
    void run_bw_threads(uns akt_fa, fa_adat_tipus akt_fa_tipus, uns gyoker_index);
    void run_pre_threads(sumhiba_tipus & kezdo_hiba, const vektor<klaszter_tartomany> & cella_kt);
    void run_post1_threads(sumhiba_tipus & akt_hiba, const vektor<klaszter_tartomany> & cella_kt);
    void run_post2_threads(rvt alfa, sumhiba_tipus & akt_hiba, const vektor<klaszter_tartomany> & cella_kt);
    void run_lepes_threads(szal_altipus szat);
    void run_belso_mentes(uns akt_anal_index, rvt sum_time);
    void debug_write(const char * fajlnev);
    //***********************************************************************
public:
    //***********************************************************************
    ~vezerlo() { /*szal_feladat_adatai adat; adat.beallitas_szaltorleshez(); feldolgozo.betesz(adat);*/ }
    //***********************************************************************
    void init_from_file(const ::std::string & fajlnev) { bemeneti_adatok.felepit_fajlbol(fajlnev); is_bemeneti_adatok_inicializalva = true; }
    //***********************************************************************
    void run_simulation(); // mindent ez csin�l, a mem�riafoglal�st�l a szimul�ci�kon �t a ment�sekig
    //***********************************************************************

    //***********************************************************************
    // k�v�lr�l vez�relhet� szimul�ci�
    //***********************************************************************
private:
    //***********************************************************************
    ::std::string sim_step_szoveg;
    uns akt_sim_index, akt_anal_index, sum_iter;
    rvt sum_time; // a tranziens szimul�ci� aktu�lis id�pontja
    rvt normal_lepes_eleje;
    rvt normal_dt;
    rvt normal_lepes_vege;
    fa_adat_tipus akt_fa_tipus;
    sumhiba_tipus akt_hiba, prev_hiba;
    uns analizis_lepes_iteracioszama;
    uns iteracioszam_stoppig, akt_szimulacio_iteracioszam;
    uns plusz_lepes;
    rvt akt_iterstop_hiba; // ami alapj�n az iter�ci�t meg�ll�tjuk
    enum iterfajta{ iter_normal, iter_dt_csokkento, iter_dt_novelo, iter_dt_stabilizalo }; 
                // dt_csokkento: t�l nagy volt a dT, visszal�p�nk id�ben
                // dt_stabilizalo: olyan, mint a norm�l l�p�s, csak nem az anal intervallum v�g�ig jutunk el, ez�rt ut�na intervalluon bel�l l�p�nk
                // dt_novelo: stikeres stabiliz�ci� ut�n elfogadjuk az aktu�lis id�pontot �s tov�bbl�p�nk, dt-t is n�velj�k. Megpr�b�ljuk esetleg a marad�k intervallumot egyszerre?
    iterfajta allapot;
    rvt akt_lepes_eleje;
    bool is_elso_iter;
    //***********************************************************************
public:
    //***********************************************************************
    void run_10_elokeszites(); // Az els� szimul�ci�s l�p�s el�tt fut.
    bool run_20_start_next_step(); // A bemeneti adatokb�l bet�lti a k�vetkez� l�p�st. Ha az el�z� szimul�ci� k�sz volt, akkor a k�v szim. els� l�p�s�t. Retunr: false, nincs t�bb l�p�s.
    void run_2A_start_external_step(const adat_analizis_lepes & lepes); // Ha k�v�lr�l, pl. guib�l adjuk meg a k�vetkez� l�p�st. (Ha van bemeneti adat, akkor az ak�r folytathat� is ez ut�n a l�p�s ut�n.)
    bool run_30_start_iteration(bool is_first_cycle); // K�tszer kell megh�vni: el�bb true-val, majd az iter�ci�k j�nnek, azt�n false-szal. Ha az is_first_cycle false, �s nem duplav�g� float, akkor nem csin�l semmit, �s false-t ad vissza, ekkor ut�na nem h�vjuk az iter�ci�kat.
    bool run_40_next_iteration(bool is_eloiteracio); // Lefuttat egy iter�ci�s l�p�st, �s visszat�r�sben jelzi, hogy k�sz van-e, vagy kell m�g iter�ci�.
    void run_80_stop_step(uns mentes_index); // bels� menti az eredm�nyt a mentes_index azonos�t�j� map-ekbe, ha external step, pl. >1M, �s a txt-ben esetleg lehetne jel�lni.
    void run_90_befejezes(); // Ha lefutott az �sszes szimul�ci�, akkor az ut�munka.
    //***********************************************************************

    //***********************************************************************
    static vezerlo & get_vezerlo() {
    // singleton lesz a vez�rl�!
    //***********************************************************************
        static vezerlo az_egyetlen_peldany;
        return az_egyetlen_peldany;
    }
    //***********************************************************************
};



}

#endif
