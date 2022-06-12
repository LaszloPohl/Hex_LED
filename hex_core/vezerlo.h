//***********************************************************************
// szimuláció fõ vezérlõ osztálya header
// Creation date:  2018. 08. 04.
// Creator:        Pohl László
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
// *******************   Globális változók    ***************************
extern bemenet bemeneti_adatok;
extern akt_sim_adatai akt_sim, prev_sim;
extern redukcios_fa<float, false>  dcfa_f;
extern redukcios_fa<double, false> dcfa_d;
extern redukcios_fa<::std::complex<float>, true>  acfa_f;
extern redukcios_fa<::std::complex<double>, true> acfa_d;
extern vektor<os_cella*> cellak; // a 0 indexû dummy!
extern vektor<rvt> csatlakozo_aramok_dc; // a 0 indexû dummy!
extern vektor<iter_csomopont> iter_csomopontok_dc; // a 0 indexû dummy!
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
    void set_ossz_lepesszam(); // megszámolja, hogy hány analízis lépés van összesen, és beállítja az összlépésszámot
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
    void run_simulation(); // mindent ez csinál, a memóriafoglalástól a szimulációkon át a mentésekig
    //***********************************************************************

    //***********************************************************************
    // kívülrõl vezérelhetõ szimuláció
    //***********************************************************************
private:
    //***********************************************************************
    ::std::string sim_step_szoveg;
    uns akt_sim_index, akt_anal_index, sum_iter;
    rvt sum_time; // a tranziens szimuláció aktuális idõpontja
    rvt normal_lepes_eleje;
    rvt normal_dt;
    rvt normal_lepes_vege;
    fa_adat_tipus akt_fa_tipus;
    sumhiba_tipus akt_hiba, prev_hiba;
    uns analizis_lepes_iteracioszama;
    uns iteracioszam_stoppig, akt_szimulacio_iteracioszam;
    uns plusz_lepes;
    rvt akt_iterstop_hiba; // ami alapján az iterációt megállítjuk
    enum iterfajta{ iter_normal, iter_dt_csokkento, iter_dt_novelo, iter_dt_stabilizalo }; 
                // dt_csokkento: túl nagy volt a dT, visszalépünk idõben
                // dt_stabilizalo: olyan, mint a normál lépés, csak nem az anal intervallum végéig jutunk el, ezért utána intervalluon belül lépünk
                // dt_novelo: stikeres stabilizáció után elfogadjuk az aktuális idõpontot és továbblépünk, dt-t is növeljük. Megpróbáljuk esetleg a maradék intervallumot egyszerre?
    iterfajta allapot;
    rvt akt_lepes_eleje;
    bool is_elso_iter;
    //***********************************************************************
public:
    //***********************************************************************
    void run_10_elokeszites(); // Az elsõ szimulációs lépés elõtt fut.
    bool run_20_start_next_step(); // A bemeneti adatokból betölti a következõ lépést. Ha az elõzõ szimuláció kész volt, akkor a köv szim. elsõ lépését. Retunr: false, nincs több lépés.
    void run_2A_start_external_step(const adat_analizis_lepes & lepes); // Ha kívülrõl, pl. guiból adjuk meg a következõ lépést. (Ha van bemeneti adat, akkor az akár folytatható is ez után a lépés után.)
    bool run_30_start_iteration(bool is_first_cycle); // Kétszer kell meghívni: elõbb true-val, majd az iterációk jönnek, aztán false-szal. Ha az is_first_cycle false, és nem duplavégû float, akkor nem csinál semmit, és false-t ad vissza, ekkor utána nem hívjuk az iterációkat.
    bool run_40_next_iteration(bool is_eloiteracio); // Lefuttat egy iterációs lépést, és visszatérésben jelzi, hogy kész van-e, vagy kell még iteráció.
    void run_80_stop_step(uns mentes_index); // belsõ menti az eredményt a mentes_index azonosítójú map-ekbe, ha external step, pl. >1M, és a txt-ben esetleg lehetne jelölni.
    void run_90_befejezes(); // Ha lefutott az összes szimuláció, akkor az utómunka.
    //***********************************************************************

    //***********************************************************************
    static vezerlo & get_vezerlo() {
    // singleton lesz a vezérlõ!
    //***********************************************************************
        static vezerlo az_egyetlen_peldany;
        return az_egyetlen_peldany;
    }
    //***********************************************************************
};



}

#endif
