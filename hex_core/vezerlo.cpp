//***********************************************************************
// aktu�lis szimul�ci� adatai + szimul�ci� f� vez�rl� oszt�lya cpp
// Creation date:  2018. 08. 04.
// Creator:        Pohl L�szl�
//***********************************************************************


//***********************************************************************
#include "vezerlo.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
// *******************   Glob�lis v�ltoz�k    ***************************
// *************** (nesze neked, programoz�s :D) ************************
bemenet bemeneti_adatok;
akt_sim_adatai akt_sim, prev_sim;
redukcios_fa<float, false>  dcfa_f;
redukcios_fa<double, false> dcfa_d;
redukcios_fa<::std::complex<float>, true>  acfa_f;
redukcios_fa<::std::complex<double>, true> acfa_d;
vektor<os_cella*> cellak; // a 0 index� dummy!
vektor<rvt> csatlakozo_aramok_dc; // a 0 index� dummy!
Para_engine feldolgozo;
//***********************************************************************
//***********************************************************************


//***********************************************************************
void vezerlo::set_ossz_lepesszam(){
//***********************************************************************
    uns db = 0;
    for (uns i = 1; i < bemeneti_adatok.szimulaciok.size(); i++)
        db += bemeneti_adatok.szimulaciok[i].analizis_lepesek.size() - 1;
    akt_lepes::set_ossz_lepesszam(db);
    akt_lepes::set_aktualis_lepes_neve("the simulation");
}


//***********************************************************************
void vezerlo::run_fw_threads(uns akt_fa, fa_adat_tipus akt_fa_tipus){
//***********************************************************************
    szal_feladat_adatai fw_klaszter_sablon, fw_egyedi_sablon;
    fw_klaszter_sablon.beallitas_fw_klaszternek();
    fw_egyedi_sablon.beallitas_fw_egyedinek();
    fw_egyedi_sablon.is_ac = fw_klaszter_sablon.is_ac = akt_sim.tipus == alt_ac;
    fw_egyedi_sablon.is_float = fw_klaszter_sablon.is_float = akt_fa_tipus == fat_float;

    // fw klaszterek ind�t�sa

    fw_klaszter_sablon.melyik_fa = akt_fa;

    const vektor<klaszter_tartomany> & kt = (akt_fa == 1)
        ? akt_sim.p_akt_sim->fa_klaszter_tartomanyok_1
        : akt_sim.p_akt_sim->fa_klaszter_tartomanyok_2;

    for (uns szal_index = 0; szal_index < kt.size(); szal_index++) {
        fw_klaszter_sablon.set_fa_kezdoindex_es_utolso_index(
            kt.unsafe(szal_index).klaszter_kezdoindex,
            kt.unsafe(szal_index).klaszter_utolso_index,
            (akt_fa == 1) ? 1 : akt_sim.p_akt_sim->gyoker_1_index + 1,
            1
            ); // a m�sodik fa az els� gy�k�releme ut�n kezd�dik
        feldolgozo.betesz(fw_klaszter_sablon, false);
    }

    // fw egyediek ind�t�sa

    fw_egyedi_sablon.melyik_fa = akt_fa;

    const vektor<adat_egyedi> & efe = (akt_fa == 1)
        ? akt_sim.p_akt_sim->egyedi_fa_elemek_1
        : akt_sim.p_akt_sim->egyedi_fa_elemek_2;

    for (uns szal_index = 0; szal_index < efe.size(); szal_index++) {
        fw_egyedi_sablon.set_fa_kezdoindex_es_utolso_index(
            efe.unsafe(szal_index).index,
            efe.unsafe(szal_index).index,
            (akt_fa == 1) ? 1 : akt_sim.p_akt_sim->gyoker_1_index + 1,
            efe.unsafe(szal_index).threads
            ); // a m�sodik fa az els� gy�k�relemeut�n kezd�dik
        feldolgozo.betesz(fw_egyedi_sablon, false);
    }

    // v�r a fw v�g�re

    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_fw_klaszter);
    feldolgozo.torli_a_keszeket(szt_fw_klaszter);
    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_fw_egyedi);
    feldolgozo.torli_a_keszeket(szt_fw_egyedi);
}


//***********************************************************************
void vezerlo::run_bw_threads(uns akt_fa, fa_adat_tipus akt_fa_tipus, uns gyoker_index){
//***********************************************************************
    szal_feladat_adatai bw_klaszter_sablon, bw_egyedi_sablon;
    bw_klaszter_sablon.beallitas_bw_klaszternek();
    bw_egyedi_sablon.beallitas_bw_egyedinek();
    bw_egyedi_sablon.is_ac = bw_klaszter_sablon.is_ac = akt_sim.tipus == alt_ac;
    bw_egyedi_sablon.is_float = bw_klaszter_sablon.is_float = akt_fa_tipus == fat_float;

    // bw egyediek ind�t�sa

    const vektor<adat_egyedi> & efe = (akt_fa == 1)
        ? akt_sim.p_akt_sim->egyedi_fa_elemek_1
        : akt_sim.p_akt_sim->egyedi_fa_elemek_2;

    bw_egyedi_sablon.melyik_fa = akt_fa;
    for (uns szal_index = 0; szal_index < efe.size(); szal_index++) {
        bw_egyedi_sablon.set_fa_kezdoindex_es_utolso_index(
            efe.unsafe(szal_index).index,
            efe.unsafe(szal_index).index,
            (akt_fa == 1) ? 1 : akt_sim.p_akt_sim->gyoker_1_index + 1,
            efe.unsafe(szal_index).threads
            ); // a m�sodik fa az els� gy�k�relemeut�n kezd�dik
        bw_egyedi_sablon.is_felteteles = gyoker_index != efe.unsafe(szal_index).index; // a legfels�t felt�tel n�lk�l elind�tja
        feldolgozo.betesz(bw_egyedi_sablon, false);
    }

    // bw klaszterek ind�t�sa

    const vektor<klaszter_tartomany> & kt = (akt_fa == 1)
        ? akt_sim.p_akt_sim->fa_klaszter_tartomanyok_1
        : akt_sim.p_akt_sim->fa_klaszter_tartomanyok_2;

    bw_klaszter_sablon.melyik_fa = akt_fa;
    for (uns szal_index = 0; szal_index < kt.size(); szal_index++) {
        bw_klaszter_sablon.set_fa_kezdoindex_es_utolso_index(
            kt.unsafe(szal_index).klaszter_kezdoindex,
            kt.unsafe(szal_index).klaszter_utolso_index,
            (akt_fa == 1) ? 1 : akt_sim.p_akt_sim->gyoker_1_index + 1,
            1
            ); // a m�sodik fa az els� gy�k�relemeut�n kezd�dik
        feldolgozo.betesz(bw_klaszter_sablon, false);
    }

    // v�r a bw v�g�re

    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_bw_egyedi);
    feldolgozo.torli_a_keszeket(szt_bw_egyedi);
    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_bw_klaszter);
    feldolgozo.torli_a_keszeket(szt_bw_klaszter);
}


//***********************************************************************
void vezerlo::run_step_settings(){
//***********************************************************************
    if (akt_sim.is_tamb_update_kell)
        update_Tamb();

    if (akt_sim.is_del_all_fa) {
        dcfa_f.clear();
        dcfa_d.clear();
        acfa_f.clear();
        acfa_d.clear();
    }

    if (akt_sim.is_uj_cellaszerkezet_kell) {
        // befejez�d�tt a ment�s? k�l�nben meg kell v�rni a v�g�t, csak ut�na lehet resize-olni a cell�kat.
        feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_belso_mentes);
        feldolgozo.torli_a_keszeket(szt_belso_mentes);
        emlekek.clear();
        for (uns i = 0; i < cellak.size(); i++)
            delete cellak.unsafe(i);
        cellak.set_size(akt_sim.p_akt_sim->cellak.size());
        for (uns i = 0; i < cellak.size(); i++)
            cellak.unsafe(i) = nullptr;
        csatlakozo_aramok_dc.set_size(akt_sim.p_akt_sim->csatlakozo_db + 1);
        csatlakozo_aramok_dc.zero();
        // A pre_init h�v�sa �tker�lt a cellaklaszterbe => nem lehet olyan face, amelyiknek a p�rja m�sik cellaklaszterben van! (Ett�l t�bb helyen f�gg a helyes m�k�d�s.)
    }

    if (akt_sim.is_uj_fa_kell) {
        if (akt_sim.tipus == alt_ac) { // AC
            if (akt_sim.fa_adat != fat_double) acfa_f.init(); // float vagy double v�g� float
            if (akt_sim.fa_adat != fat_float)  acfa_d.init(); // double vagy double v�g� float
        }
        else { // DC vagy tranziens
            if (akt_sim.fa_adat != fat_double) dcfa_f.init(); // float vagy double v�g� float
            if (akt_sim.fa_adat != fat_float)  dcfa_d.init(); // double vagy double v�g� float
        }
    }

    // ha is_del_all_prev vagy is_tamb_update_kell van, akkor 4-et kell l�ptetni az eml�kez� t�rat (ekkor a mentend� �s a megtartand� marad csak, a t�bbi 0 lesz a clear h�v�s miatt)

    if (akt_sim.is_del_all_prev || akt_sim.is_tamb_update_kell) { /*emlekezo::leptet1(); emlekezo::leptet1(); emlekezo::leptet1(); emlekezo::leptet1();*/ }
}


//***********************************************************************
void vezerlo::run_pre_threads(sumhiba_tipus & kezdo_hiba, const vektor<klaszter_tartomany> & cella_kt){
// pre cella klaszterek ind�t�sa
//***********************************************************************
    szal_feladat_adatai cellafeldolgozo_klaszter_sablon_pre;
    cellafeldolgozo_klaszter_sablon_pre.beallitas_cellafeldolgozo_klaszternek_pre();
    cellafeldolgozo_klaszter_sablon_pre.is_ac = akt_sim.tipus == alt_ac;
    akt_sim.alfa = rvt(1);
    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_cellafeldolgozo_klaszter);
    for (uns szal_index = 0; szal_index < cella_kt.size(); szal_index++) {
        cellafeldolgozo_klaszter_sablon_pre.set_cellatomb_kezdoindex_es_utolso_index(
            cella_kt.unsafe(szal_index).klaszter_kezdoindex, cella_kt.unsafe(szal_index).klaszter_utolso_index);
        feldolgozo.betesz(cellafeldolgozo_klaszter_sablon_pre, false);
    }
    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_cellafeldolgozo_klaszter);
    kezdo_hiba = feldolgozo.get_sumhiba(szt_cellafeldolgozo_klaszter);
    feldolgozo.torli_a_keszeket(szt_cellafeldolgozo_klaszter);
}


//***********************************************************************
void vezerlo::run_lepes_threads(szal_altipus szat) {
// pre cella klaszterek ind�t�sa
//***********************************************************************
    size_t szalszam = bemeneti_adatok.cputhreads;
    size_t db = emlekek.get_size();
    szal_feladat_adatai lepo_klaszter_sablon;
    lepo_klaszter_sablon.beallitas_lepo_klaszternek(szat);
    lepo_klaszter_sablon.is_ac = akt_sim.tipus == alt_ac;
    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_lepo_klaszter);
    for (uns szal_index = 0; szal_index < bemeneti_adatok.cputhreads; szal_index++) {
        lepo_klaszter_sablon.set_cellatomb_kezdoindex_es_utolso_index(
            (uns)((db * szal_index) / szalszam), (uns)((db * (szal_index + 1)) / szalszam));
        feldolgozo.betesz(lepo_klaszter_sablon, false);
    }
    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_lepo_klaszter);
    feldolgozo.torli_a_keszeket(szt_lepo_klaszter);
}


//***********************************************************************
void vezerlo::run_apply_externals(uns gyoker_index) {
// externalok �rv�nyes�t�se
//***********************************************************************
    if (akt_sim.p_akt_sim->fa_elemek[gyoker_index].A != 0)
        throw hiba(1, "%u non-reducted node remained after forwsubs", akt_sim.p_akt_sim->fa_elemek[gyoker_index].A);
}


//***********************************************************************
void vezerlo::run_post1_threads(sumhiba_tipus & akt_hiba, const vektor<klaszter_tartomany> & cella_kt){
// post1 cella klaszterek ind�t�sa
//***********************************************************************
    szal_feladat_adatai cellafeldolgozo_klaszter_sablon_post1;
    cellafeldolgozo_klaszter_sablon_post1.beallitas_cellafeldolgozo_klaszternek_post();
    cellafeldolgozo_klaszter_sablon_post1.is_ac = akt_sim.tipus == alt_ac;
    akt_sim.alfa = rvt(1);
    for (uns szal_index = 0; szal_index < cella_kt.size(); szal_index++) {
        cellafeldolgozo_klaszter_sablon_post1.set_cellatomb_kezdoindex_es_utolso_index(
            cella_kt.unsafe(szal_index).klaszter_kezdoindex, cella_kt.unsafe(szal_index).klaszter_utolso_index);
        feldolgozo.betesz(cellafeldolgozo_klaszter_sablon_post1, false);
    }
    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_cellafeldolgozo_klaszter);
    // Val�sz�n�leg jobb k�t fa eset�n �sszeadni a hib�kat, mint csak a m�sodikat n�zni.
    sumhiba_tipus sumhiba = feldolgozo.get_sumhiba(szt_cellafeldolgozo_klaszter);
    sumhiba.to_sum(akt_hiba);
    feldolgozo.torli_a_keszeket(szt_cellafeldolgozo_klaszter);
}


//***********************************************************************
void vezerlo::run_post2_threads(rvt alfa, sumhiba_tipus & akt_hiba, const vektor<klaszter_tartomany> & cella_kt){
// post2 cella klaszterek ind�t�sa
//***********************************************************************
    szal_feladat_adatai cellafeldolgozo_klaszter_sablon_post2;
    cellafeldolgozo_klaszter_sablon_post2.beallitas_cellafeldolgozo_klaszternek_post();
    cellafeldolgozo_klaszter_sablon_post2.is_ac = akt_sim.tipus == alt_ac;
    akt_sim.alfa = alfa;
    for (uns szal_index = 0; szal_index < cella_kt.size(); szal_index++) {
        cellafeldolgozo_klaszter_sablon_post2.set_cellatomb_kezdoindex_es_utolso_index(
            cella_kt.unsafe(szal_index).klaszter_kezdoindex, cella_kt.unsafe(szal_index).klaszter_utolso_index);
        feldolgozo.betesz(cellafeldolgozo_klaszter_sablon_post2, false);
    }
    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_cellafeldolgozo_klaszter);
    sumhiba_tipus sumhiba = feldolgozo.get_sumhiba(szt_cellafeldolgozo_klaszter);
    sumhiba.to_sum(akt_hiba);
    feldolgozo.torli_a_keszeket(szt_cellafeldolgozo_klaszter);
}


//***********************************************************************
void vezerlo::run_belso_mentes(uns akt_anal_index, rvt sum_time){
//***********************************************************************
    szal_feladat_adatai belso_mentes_sablon;
    belso_mentes_sablon.beallitas_belso_mentesnek();
    belso_mentes_sablon.is_ac = akt_sim.tipus == alt_ac;
    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_belso_mentes); // csak akkor kezdi menteni az �j l�p�s eredm�ny�t, ha k�sz a r�gi
    feldolgozo.torli_a_keszeket(szt_belso_mentes);
    run_lepes_threads(szat_megtart_akt); //emlekek.megtartando_az_aktualis();
    run_lepes_threads(szat_ment_akt); //emlekek.mentendo_az_aktualis();
    // ment�s1 => ha k�sz passzolja az eredm�nyt tov�bb, ami ment�s2-h�z ker�l, aki t�nyleg menti; Tamb-ot is adja �t, mert v�ltozhat!
    // Tamb legyeb a szal_feladat_adatai param�tere, �s kiolvas�skor azonnal adja a h�m�rs�klethez!
    // �j cellastrukt�ra l�trehoz�sa el�tt be kell fejezni a ment�st
    belso_mentes_sablon.Tamb = akt_sim.Tamb;
    belso_mentes_sablon.p_akt_eredm = akt_sim.p_akt_eredm;
    belso_mentes_sablon.eredm_utvonal = akt_sim.p_akt_sim->eredm_utvonal;
    belso_mentes_sablon.analizis_tipus = akt_sim.tipus;
    belso_mentes_sablon.analizis_value = akt_sim.tipus == alt_dc ? rvt() : akt_sim.tipus == alt_ac ? akt_sim.value : sum_time;
    belso_mentes_sablon.is_use_commas = akt_sim.p_akt_sim->is_use_commas;
    belso_mentes_sablon.akt_anal_index = akt_anal_index;
    feldolgozo.betesz(belso_mentes_sablon, false);
}


//***********************************************************************
void vezerlo::debug_write(const char * fajlnev){
//***********************************************************************
    ::std::ofstream fs(fajlnev);
    if (fs.is_open()) {
        for (uns i = 1; i < cellak.size(); i++) {
            fs << "cell " << i << ":" << ::std::endl;
            cellak.unsafe(i)->debug_write(fs);
        }
        for (uns i = 0; i < dcfa_d.fa_1.size(); i++) { // 0-t�l indexel
            fs << "tree " << (i + 1) << ":" << ::std::endl;
            dcfa_d.fa_1.unsafe(i).debug_write(fs);
        }
        fs.close();
    }
}

#define JUNCTION_PRECALC

//***********************************************************************
void vezerlo::run_simulation/*_darabokbol*/() {
//***********************************************************************
    run_10_elokeszites();
//    uns maxlepes = 500;
    while (run_20_start_next_step()) {
#ifdef JUNCTION_PRECALC
        LED_model_result_pack * el_pack = akt_sim.p_akt_sim->junctionok[1].el_egyenlet.set_junction_pre_calc(/*0.25*/0.05, 0, true, nullptr);
        akt_sim.p_akt_sim->junctionok[1].rad.set_junction_pre_calc(/*0.25*/0.05, 0, false, el_pack);
#endif
        if (run_30_start_iteration(true)) {
#ifdef JUNCTION_PRECALC

            while (run_40_next_iteration()) {
//                maxlepes--;
//                if (maxlepes == 0)
//                    break;
            }
            printf("***\n");
            akt_sim.p_akt_sim->junctionok[1].el_egyenlet.set_junction_normal_calc();
            akt_sim.p_akt_sim->junctionok[1].rad.set_junction_normal_calc();
            while (run_40_next_iteration()) {
                //                maxlepes--;
                //                if (maxlepes == 0)
                //                    break;
            }
#endif
            while (run_40_next_iteration()) {
                //                maxlepes--;
                //                if (maxlepes == 0)
                //                    break;
            }
        }
        if (run_30_start_iteration(false)) {
            while (run_40_next_iteration())
                ;
        }
        run_80_stop_step(akt_anal_index);
    }
    run_90_befejezes();
}


//***********************************************************************
void vezerlo::run_10_elokeszites() {
// Az els� szimul�ci�s l�p�s el�tt fut.
//***********************************************************************
    set_hiba_hol h("vezerlo::run_simulation");

    if (!is_bemeneti_adatok_inicializalva)
        throw hiba(1, "No input data.");
    set_ossz_lepesszam();

    feldolgozo.set_szalszam(bemeneti_adatok.cputhreads);

    akt_sim_index = 1;
    akt_anal_index = 0;
    sum_iter = 0;
    sum_time = rvt();
}


//***********************************************************************
bool vezerlo::run_20_start_next_step() {
// A bemeneti adatokb�l bet�lti a k�vetkez� l�p�st. Ha az el�z� szimul�ci�
// k�sz volt, akkor a k�v szim. els� l�p�s�t. Return: false, nincs t�bb l�p�s.
//***********************************************************************

    // Az aktu�lis szimul�ci�s l�p�s kiv�laszt�sa

    if (akt_sim_index >= bemeneti_adatok.szimulaciok.size())
        return false;
    akt_anal_index++;
    while (akt_anal_index >= bemeneti_adatok.szimulaciok[akt_sim_index].analizis_lepesek.size()) {
        akt_sim_index++;
        if (akt_sim_index >= bemeneti_adatok.szimulaciok.size())
            return false;
        akt_anal_index = 1;
        sum_time = rvt();
    }

    // A szimul�ci�s l�p�s be�ll�t�sa 

    sim_step_szoveg = "sim: " + ::std::to_string(akt_sim_index) + "/" + ::std::to_string(bemeneti_adatok.szimulaciok.size() - 1)
        + ", step: " + ::std::to_string(akt_anal_index) + "/"
        + ::std::to_string(bemeneti_adatok.szimulaciok[akt_sim_index].analizis_lepesek.size() - 1);
    if (akt_anal_index < 10 || akt_anal_index > bemeneti_adatok.szimulaciok[akt_sim_index].analizis_lepesek.size() - 3)
        most("run_simulation start (" + sim_step_szoveg + ")");

    prev_sim = akt_sim;

    akt_sim.set_akt_anal_step(
        akt_anal_index == 1 ? &bemeneti_adatok.szimulaciok[akt_sim_index] : nullptr,
        bemeneti_adatok.szimulaciok[akt_sim_index].analizis_lepesek[akt_anal_index]);

    akt_lepes::inc_akt_lepesszam();

    normal_lepes_eleje = sum_time;
    if (akt_sim.tipus == alt_trans)
        sum_time += akt_sim.value; // dt
    else
        sum_time = rvt();
    normal_dt = akt_sim.value;
    normal_lepes_vege = sum_time;

    run_step_settings(); // Akt step alapj�n be�ll�t�sok, pl. fa t�rl�se
    akt_fa_tipus = (akt_sim.fa_adat != fat_double) ? fat_float : fat_double;
    //bool is_iteracio_vege = false;

    // pre cella klaszterek ind�t�sa

    sumhiba_tipus kezdo_hiba;
    run_pre_threads(kezdo_hiba, akt_sim.p_akt_sim->cella_klaszter_tartomanyok);

/*
    for (uns i = 1; i < cellak.size(); i++)
        cellak.unsafe(i).randomize_zaj();
*/
    // iter�ci� 

    if (akt_anal_index < 3)most("run_simulation pre cella cluster is done");
    prev_hiba.zero();
    akt_hiba.zero();
    akt_hiba.sum_IP_hiba = kezdo_hiba.sum_IP_hiba;
    analizis_lepes_iteracioszama = 0;
    akt_sim.is_subiter_dt_changed = false;

    return true;
}


//***********************************************************************
void vezerlo::run_2A_start_external_step(const adat_analizis_lepes & lepes) {
// Ha k�v�lr�l, pl. guib�l adjuk meg a k�vetkez� l�p�st. (Ha van bemeneti
// adat, akkor az ak�r folytathat� is ez ut�n a l�p�s ut�n.)
//***********************************************************************
    TODO("vezerlo::run_2A_start_external_step");
}


//***********************************************************************
bool vezerlo::run_30_start_iteration(bool is_first_cycle) {
// K�tszer kell megh�vni: el�bb true-val, majd az iter�ci�k j�nnek, azt�n
// false-szal. Ha az is_first_cycle false, �s nem duplav�g� float, akkor
// nem csin�l semmit, �s false-t ad vissza, ekkor ut�na nem h�vjuk az iter�ci�kat.
//***********************************************************************
    if (!is_first_cycle) {
        if (akt_sim.fa_adat == fat_double_ended_float && akt_fa_tipus == fat_float)
            akt_fa_tipus = fat_double;
        else
            return false;
    }
    if (akt_anal_index < 3)most("run_simulation start fw");
    iteracioszam_stoppig = 0;
    plusz_lepes = 0;
    akt_iterstop_hiba = rvt(); // ami alapj�n az iter�ci�t meg�ll�tjuk
    allapot = iter_normal;
    akt_lepes_eleje = normal_lepes_eleje;
    // normal_lepes_eleje, normal_dt, normal_lepes_vege
    is_elso_iter = true;
    return true;
}


//***********************************************************************
bool vezerlo::run_40_next_iteration() {
// Lefuttat egy iter�ci�s l�p�st, �s visszat�r�sben jelzi, hogy kell-e m�g iter�ci�.
//***********************************************************************
    char itertip = '-';
    if (allapot == iter_dt_csokkento)itertip = 'v';
    if (allapot == iter_dt_novelo)itertip = '^';
    if (allapot == iter_dt_stabilizalo)itertip = '#';
    analizis_lepes_iteracioszama++;
    sum_iter++;
    sumhiba_tipus dummy_hiba;
    const vektor<klaszter_tartomany> & cella_kt = akt_sim.p_akt_sim->cella_klaszter_tartomanyok;
    switch (allapot) {
        case iter_normal:
            run_lepes_threads(szat_kiind_akt_es_elore);
            iteracioszam_stoppig++;
            akt_lepes::set_aktualis_lepes_neve(::std::string("the simulation (")
                + sim_step_szoveg + ", iter: " + ::std::to_string(iteracioszam_stoppig) + ") normal");
            prev_hiba = akt_hiba;
            break;
        case iter_dt_stabilizalo:
            run_lepes_threads(szat_kiind_akt_es_elore);
            iteracioszam_stoppig++;
            akt_lepes::set_aktualis_lepes_neve(::std::string("the simulation (")
                + sim_step_szoveg + ", iter: " + ::std::to_string(iteracioszam_stoppig) + ") stabilizing");
            prev_hiba = akt_hiba;
            break;
        case iter_dt_csokkento:
            akt_sim.value /= 2 + 4 * akt_hiba.max_T_hiba;  //*= akt_sim.value<1e-10 ? 0.1 : 0.5; // fel�re cs�kkentj�k a l�p�sk�zt
            akt_sim.is_subiter_dt_changed = true; // �jrasz�moland�k a kapacit�sok
            akt_sim.is_iter_dt_csokkento = true; // ne sz�molja �jra a fesz�lts�geket, mert az akt a kiindul�
            sum_time = akt_lepes_eleje + akt_sim.value;
            run_lepes_threads(szat_vissza); //emlekek.visszalep();
            run_post1_threads(dummy_hiba, cella_kt);
            run_lepes_threads(szat_elore); //emlekek.leptet();
            plusz_lepes++;
            iteracioszam_stoppig = 1;
            akt_lepes::set_aktualis_lepes_neve(::std::string("the simulation (")
                + sim_step_szoveg + ", iter: " + ::std::to_string(iteracioszam_stoppig) + ") decreasing");
            break;
        case iter_dt_novelo:
            run_lepes_threads(szat_megtart_akt); //emlekek.megtartando_az_aktualis();
            akt_lepes_eleje = sum_time; // el�z� l�p�s elfogad�sa
            if (akt_lepes_eleje + (akt_sim.value<1e-11 ? rvt(20) : rvt(4)) * akt_sim.value > normal_lepes_vege) {
                akt_sim.value = normal_lepes_vege - akt_lepes_eleje;
                sum_time = normal_lepes_vege;
                allapot = iter_normal;
            }
            else {
                akt_sim.value *= (akt_sim.value<1e-11 ? rvt(10) : rvt(2));
                sum_time = akt_lepes_eleje + akt_sim.value;
                allapot = iter_dt_stabilizalo;
            }
            akt_sim.is_subiter_dt_changed = true;
            akt_sim.is_iter_dt_csokkento = true; //nem sz�mohat fesz�lts�geket
            run_post1_threads(dummy_hiba, cella_kt);
            run_lepes_threads(szat_kiind_akt_es_elore);
            plusz_lepes++;
            iteracioszam_stoppig = 1;
            akt_lepes::set_aktualis_lepes_neve(::std::string("the simulation (")
                + sim_step_szoveg + ", iter: " + ::std::to_string(iteracioszam_stoppig) + ") increasing");
            prev_hiba = akt_hiba;
            is_elso_iter = true;
            break;
    }

    if (analizis_lepes_iteracioszama > 1)
        akt_sim.set_akt_iter();

    akt_hiba.zero();
    for (uns akt_fa = 1; akt_fa <= akt_sim.p_akt_sim->db_fa; akt_fa++) { // az �sszes f�t megcsin�lja
        const uns gyoker_index = (akt_fa == 1) ? akt_sim.p_akt_sim->gyoker_1_index : akt_sim.p_akt_sim->gyoker_2_index;
        run_fw_threads(akt_fa, akt_fa_tipus);                                                   if (akt_anal_index < 3)most("run_simulation fw is done");
        run_apply_externals(gyoker_index);
        run_bw_threads(akt_fa, akt_fa_tipus, gyoker_index);                                     if (akt_anal_index < 3)most("run_simulation bw is done");
        run_post1_threads(akt_hiba, cella_kt); if (akt_anal_index < 3)most("run_simulation post cella cluster is done");
    }

    // Konvergencia jav�t�sa

    if (akt_sim.tipus == alt_trans) {
        const rvt csokkento_max_T_hiba = rvt(1); // mondjuk max. 1 K ugr�st enged�nk meg.
        if (akt_sim.value>1e-15 && (akt_hiba.max_T_hiba > (is_elso_iter ? (5 * csokkento_max_T_hiba) : csokkento_max_T_hiba) || abs((akt_hiba.max_T_hiba-prev_hiba.max_T_hiba)/akt_hiba.max_T_hiba) < 1e-5)) {
            allapot = iter_dt_csokkento;
        }
        else {
            if (allapot == iter_dt_csokkento)
                allapot = iter_dt_stabilizalo;
        }
    }

    if ((allapot==iter_normal || allapot==iter_dt_stabilizalo) && akt_sim.p_akt_sim->nemlin_tipus != nmt_klasszik_iteracio && iteracioszam_stoppig > 1) { // Az els� iter�ci� ut�n a prev_hib�khoz nincs �rtelme hasonl�tani (ak�r �j l�p�s, ak�r float/double v�lt�s)
        // Az IP hiba a hib�k �sszege, az UT hiba a hib�k �tlaga
        if (akt_hiba.sum_IP_hiba >= akt_sim.max_hiba || (prev_hiba.sum_IP_hiba >= akt_sim.max_hiba && akt_hiba.sum_IP_hiba >= 0.1*prev_hiba.sum_IP_hiba)) {
            if (akt_hiba.sum_IP_hiba < 100 * prev_hiba.sum_IP_hiba) {
                rvt alfa = prev_hiba.sum_IP_hiba / (prev_hiba.sum_IP_hiba + akt_hiba.sum_IP_hiba);
                if (alfa < 0.25)
                    alfa = 0.25;
                run_lepes_threads(szat_elore); //emlekek.leptet();
                sumhiba_tipus uj_hiba;
                                
                run_post2_threads(alfa, uj_hiba, cella_kt); // post2 cella klaszterek ind�t�sa
                                
                if (uj_hiba.sum_IP_hiba > akt_hiba.sum_IP_hiba) {
                    run_lepes_threads(szat_vissza); //emlekek.visszalep();
                }
                else {
                    akt_hiba = uj_hiba;
                }
                // A gyakorlat alapj�n lehet m�g finom�tani, nem teljesen k�veti a trtr-t. Ill. l�sd m�g az iteracio.pdf-et.
            }
            else { // akt hiba >= 100*el�z� hiba
                rvt alfa = rvt(1);
                const uns sorozatszam = 3;
                for (uns i = 0; i < sorozatszam; i++) {
                    alfa *= rvt(0.4);
                    run_lepes_threads(szat_elore); //emlekek.leptet();
                    sumhiba_tipus uj_hiba;

                    run_post2_threads(alfa, uj_hiba, cella_kt); // post2 cella klaszterek ind�t�sa
                                    
                    if (uj_hiba.sum_IP_hiba > rvt(2)*prev_hiba.sum_IP_hiba)
                        if (i < sorozatszam - 1) {
                            run_lepes_threads(szat_vissza); //emlekek.visszalep();
                        }
                        else {
                            akt_hiba = uj_hiba;
                        }
                    else {
                        akt_hiba = uj_hiba;
                        break;
                    }
                }
            }
        }
    }
    // akt_iterstop_hiba = akt_IP_hiba < akt_sim.max_hiba*1.0e-6 ? akt_IP_hiba : akt_IP_hiba + akt_UT_hiba;
    akt_iterstop_hiba = akt_hiba.sum_IP_hiba + akt_hiba.max_UT_hiba;
    if (allapot == iter_dt_stabilizalo && (akt_iterstop_hiba < akt_sim.max_hiba || iteracioszam_stoppig >= akt_sim.max_iter))
        allapot = iter_dt_novelo;
    // ide esetleg lehet tenni egy opci�s ment�st, ha n�zni akarjuk a konvergenci�t

    //rvt atlag_IP_hiba = sqrt(akt_IP_hiba) / (akt_center_csomopont_db + csatlakozo_aramok_dc.size() - 1);
    printf("anal=%-4u iter=%-3u sum_iter=%-5u sum_time=%-12.6g dt=%-9.3g ", akt_anal_index, iteracioszam_stoppig, sum_iter, sum_time, akt_sim.value);
    printf("iter_err=%-9.3g IP_err=%-9.3g VT_err=%-9.3g T_err=%-9.3g alfa=%.3g", akt_iterstop_hiba, akt_hiba.sum_IP_hiba, akt_hiba.max_UT_hiba, akt_hiba.max_T_hiba, akt_sim.alfa);
//                    printf("VT_err=%-9.3g %c T_err=%-9.3g alfa=%.3g", akt_UT_hiba, itertip, akt_T_hiba, akt_sim.alfa);
/*
    uns kiirando_T_index = 22436;
    if (cellak.size() > kiirando_T_index) {
        if (cellak[kiirando_T_index].is_th)printf(" %.3gK", cellak[kiirando_T_index].th_center_face_dc.emlekek.get_akt().UT); // 10273
        if (cellak[kiirando_T_index].is_el)printf(" %.3gV", cellak[kiirando_T_index].el_center_face_dc.emlekek.get_akt().UT); // 10273
        printf("\n");
    }
    else
*/        printf("\n");

    //printf("Tc=%g\tHm=%g\tHa=%g\n", cellak[kiirando_T_index].th_center_face_dc.UT.get_akt(), cellak[kiirando_T_index].H.get_megtartando().ertek, cellak[kiirando_T_index].H.get_akt().ertek);
/*                    if (sum_iter == 3)
        debug_write("jo.txt");

    if (sum_iter == 298)
        debug_write("ti_298.txt");
    if (sum_iter == 299)
        debug_write("ti_299.txt");
    if (sum_iter == 300)
        debug_write("ti_300.txt");
    if (sum_iter == 301)
        debug_write("ti_301.txt");
*/
    is_elso_iter = false;
    return (allapot != iter_normal || (akt_iterstop_hiba >= akt_sim.max_hiba && iteracioszam_stoppig < akt_sim.max_iter));
}


//***********************************************************************
void vezerlo::run_80_stop_step(uns mentes_index) {
// bels� menti az eredm�nyt a mentes_index azonos�t�j� map-ekbe, ha
// external step, pl. >1M, �s a txt-ben esetleg lehetne jel�lni.
//***********************************************************************
    run_belso_mentes(mentes_index, sum_time); // Ha k�ls�, akkor mentes_index speci�lis, pl. 1000000-t�l induljon
    if (mentes_index < 3)
        most(analizis_lepes_iteracioszama == 1
            ? "run_simulation sim step is done"
            : ::std::string("run_simulation sim step is done (iter: ") + ::std::to_string(analizis_lepes_iteracioszama) + ")");
}


//***********************************************************************
void vezerlo::run_90_befejezes() {
// Ha lefutott az �sszes szimul�ci�, akkor az ut�munka.
//***********************************************************************
    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_belso_mentes);
    feldolgozo.torli_a_keszeket(szt_belso_mentes);
    feldolgozo.var_mig_van_varo_vagy_futo_feladat(szt_fajlba_mentes);
    feldolgozo.torli_a_keszeket(szt_fajlba_mentes);
    akt_lepes::set_ossz_lepesszam(0);
    log_print::log_print(1, "the simulation: OK                                                \n");
}


#include <time.h>

//***********************************************************************
void akt_sim_adatai::set_akt_anal_step(adat_szimulacio * p_uj_sim, const adat_analizis_lepes & uj_lepes) {
//***********************************************************************
    set_hiba_hol h("akt_sim_adatai::run_simulation");

    clear_nem_oroklodok();
    if (tipus != uj_lepes.tipus || value != uj_lepes.value)
        is_force_face_update = true; // a kapacit�sok miatt a face-ek friss�t�se sz�ks�ges
    tipus = uj_lepes.tipus;
    value = uj_lepes.value;

    // �j szimul�ci� eset�n alaphelyzetbe �ll�t�s

    if (p_uj_sim != nullptr) {
        p_akt_sim = p_uj_sim;
        reset_akt_perem_index();
        reset_akt_perem(); // ez m�dos�that�
        reset_akt_anyag_index();
        reset_akt_junct_index();
        reset_akt_strukturak(); // ez m�dos�that�, tartalmazza a gerjeszt�sek c�m�t is

        if (p_akt_sim->is_zaj_random) {
            uns initnum = (uns)time(nullptr) % RAND_MAX;
            srand(initnum);
        }
        else {
            srand(1976);
        }

        // mi �j?

        if (!p_akt_sim->is_KM) { // �j anyagok
            is_uj_facek_letrehozasa_kell = true;
        }
        if (!p_akt_sim->is_KC) { // �j strukt�r�k
            is_uj_facek_letrehozasa_kell = true;
        }
        if (!p_akt_sim->is_KJ) { // �j junctionok
            is_uj_facek_letrehozasa_kell = true;
        }
        if (!p_akt_sim->is_KP) { // �j peremkelt�telek
            is_uj_facek_letrehozasa_kell = true;
            is_tamb_update_kell = true; // ?
            is_peremfelt_update_kell = true;
        }
        // !p_akt_sim->is_KA : nincs tennival� // �j anal�zis l�p�sek
        if (!p_akt_sim->is_KL) { // �j cellastrukt�ra
            is_uj_fa_kell = true;
            is_uj_cellaszerkezet_kell = true;
            is_uj_facek_letrehozasa_kell = true;
            is_tamb_update_kell = true;
            is_gerj_update_kell = true;
            is_peremfelt_update_kell = true;
            is_del_all_prev = true;
            is_del_all_fa = true; // t�r�lj�k a r�gi f�kat, nehogy bezavarjanak pl. adatt�pus m�dos�t�sn�l
        }
        if (!p_akt_sim->is_KD) { // �j m�retek
            is_uj_facek_letrehozasa_kell = true;
            is_peremfelt_update_kell = true;
        }
        if (!p_akt_sim->is_KF) { // �j fa
            is_uj_fa_kell = true;
        }
        if (p_akt_sim->fa_adat != fa_adat) { // ha m�s adatt�pussal alarunk dolgozni
            fa_adat = p_akt_sim->fa_adat;
            // nem felt�tlen�l kell �j fa, �gyhogy ezt kiveszem: is_uj_fa_kell = true;
            // a r�gi f�kat nem kell t�r�lni, mert lehet, hogy haszn�ljuk m�g
        }
        if (!p_akt_sim->is_KR) {
            p_akt_eredm = nullptr;
        }
    }

    // be�ll�t�sok alkalmaz�sa

    for (uns i = 0; i < uj_lepes.beallitasok.size(); i++) {
        switch (uj_lepes.beallitasok[i].tipus) {
        case abt_I0: I0 = uj_lepes.beallitasok[i].rvt_ertek; break;
        case abt_max_error: max_hiba = uj_lepes.beallitasok[i].rvt_ertek; break;
        case abt_max_iter:  max_iter = uj_lepes.beallitasok[i].uns_ertek; break;
        case abt_del_all_prev: is_del_all_prev = true; break;
        case abt_del_all_gerj: del_all_gerj(); break;
        case abt_del_fa: is_del_all_fa = true; is_uj_fa_kell = true; break;
        case abt_tamb: Tamb = uj_lepes.beallitasok[i].rvt_ertek; is_tamb_update_kell = true; break;
        case abt_change_bn:
            akt_perem_index[uj_lepes.beallitasok[i].uns_ertek] = uj_lepes.beallitasok[i].uns_ertek_2;
            is_uj_facek_letrehozasa_kell = true;
            break;
        case abt_change_bv:
            switch (akt_perem[uj_lepes.beallitasok[i].uns_ertek].tipus) {
            case pt_u:
            case pt_t:
                akt_perem[uj_lepes.beallitasok[i].uns_ertek].UT = uj_lepes.beallitasok[i].rvt_ertek;
                break;
            case pt_htc:
                akt_perem[uj_lepes.beallitasok[i].uns_ertek].htc = uj_lepes.beallitasok[i].rvt_ertek;
                break;
            case pt_thtc:
                akt_perem[uj_lepes.beallitasok[i].uns_ertek].UT = uj_lepes.beallitasok[i].rvt_ertek_2;
                akt_perem[uj_lepes.beallitasok[i].uns_ertek].htc = uj_lepes.beallitasok[i].rvt_ertek;
                break;
            }
            is_peremfelt_update_kell = true;
            break;
        case abt_change_jn:
            akt_junct_index[uj_lepes.beallitasok[i].uns_ertek] = uj_lepes.beallitasok[i].uns_ertek_2;
            is_uj_facek_letrehozasa_kell = true;
            break;
        case abt_change_mn:
            akt_anyag_index[uj_lepes.beallitasok[i].uns_ertek] = uj_lepes.beallitasok[i].uns_ertek_2;
            is_uj_facek_letrehozasa_kell = true;
            break;
        case abt_ct:
            switch (uj_lepes.beallitasok[i].uns_ertek) {
                case 1:
                    //if (fa_adat != fat_double)
                    //    is_uj_fa_kell = true;
                    fa_adat = fat_double;
                    break;
                case 2:
                    //if (fa_adat != fat_float)
                    //    is_uj_fa_kell = true;
                    fa_adat = fat_float;
                    break;
                case 3:
                    //if (fa_adat != fat_double_ended_float)
                    //    is_uj_fa_kell = true;
                    fa_adat = fat_double_ended_float;
                    break;
                default:
                    throw hiba(1, "abt_ct: unknown tree type (%u)", uj_lepes.beallitasok[i].uns_ertek);
            }
            break;
        }
    }

    // gerjeszt�sek be�ll�t�sa (nem cser�lhet� fel a be�ll�t�sok alkalmaz�s�val, mert ott szerepelhet a gerjeszt�sek t�rl�se)

    for (uns i = 0; i < uj_lepes.gerjesztesek.size(); i++) {
        is_gerj_update_kell = true;
        switch (uj_lepes.gerjesztesek[i].tipus) {
            case gt_U:
            case gt_I: akt_strukturak[uj_lepes.gerjesztesek[i].struktura_index].p_el_gerj = &uj_lepes.gerjesztesek[i]; break;
            case gt_P:
            case gt_T: akt_strukturak[uj_lepes.gerjesztesek[i].struktura_index].p_th_gerj = &uj_lepes.gerjesztesek[i]; break;
            case gt_el_none: akt_strukturak[uj_lepes.gerjesztesek[i].struktura_index].p_el_gerj = nullptr; break;
            case gt_th_none: akt_strukturak[uj_lepes.gerjesztesek[i].struktura_index].p_th_gerj = nullptr; break;
            default: throw hiba(1, "unknown excitation type");
        }
    }

    // eredm�nyek be�ll�t�sa: csak akkor, ha van megadva �j eredm�nyk�r�s

    if (uj_lepes.eredmenyek.size()>0)
        p_akt_eredm = &uj_lepes.eredmenyek;
}


}
