//***********************************************************************
// párhuzamos feldolgozó cpp
// Creation date:  2018. 08. 09.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#include "Parhuzamos_Feldolgozo.h"
#include "vezerlo.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
template<typename adattipus> 
void parallel_seged<adattipus>::run_job(hivando_matrixfuggveny<adattipus> & fv) {
//***********************************************************************
    szal_feladat_adatai adatok;
    adatok.beallitas_matrixmuveletnek(fv);
    feldolgozo.betesz(adatok, fv.s_is_priority);
}


//***********************************************************************
const szaltipus_jellemzoi szaltipus_jellemzoi_tomb[szal_tipus::szt_N]
= {// valid,                        egyszalu,   azonnal_torolheto,  felteteles, ertesitendo_felteteles
    { szt_nem_csinal_semmit,        false,      false,              false,      szt_N },
    { szt_cellafeldolgozo_klaszter, false,      false,              false,      szt_N },
    { szt_lepo_klaszter,            false,      true,               false,      szt_N },
    { szt_fw_klaszter,              false,      true,               false,      szt_fw_egyedi },
    { szt_fw_egyedi,                false,      true,               true,       szt_N },
    { szt_bw_egyedi,                false,      true,               true,       szt_bw_klaszter },
    { szt_bw_klaszter,              false,      true,               true,       szt_N },
    { szt_matrix,                   false,      true,               false,      szt_N },
    { szt_belso_mentes,             true,       true,               false,      szt_N },
    { szt_fajlba_mentes,            true,       true,               false,      szt_N },
};
//***********************************************************************


//***********************************************************************
void run_fw_klaszter(szal_feladat_adatai * padat) {
// A klaszter mátrix szálszáma mindig 1
//***********************************************************************
//printf("\nfw klaszter start\n");
    const uns utolso = padat->klaszter_utolso_index;
    const vektor<adat_fa_elem> & fa_elemek = akt_sim.p_akt_sim->fa_elemek;
    if (padat->is_ac) {
        if (padat->is_float) {
            if (padat->melyik_fa == 1) {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                for (; i_be <= utolso; i_be++, i_fa++)
                    acfa_f.fa_1[i_fa].forward(1, fa_elemek[i_be], 1, acfa_f.fa_1); // egyelõre az unsafe helyett maradjon az [], hátha nem jó
            }
            else {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                for (; i_be <= utolso; i_be++, i_fa++)
                    acfa_f.fa_2[i_fa].forward(1, fa_elemek[i_be], 2, acfa_f.fa_2);
            }
        }
        else {
            if (padat->melyik_fa == 1) {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                for (; i_be <= utolso; i_be++, i_fa++)
                    acfa_d.fa_1[i_fa].forward(1, fa_elemek[i_be], 1, acfa_d.fa_1);
            }
            else {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                for (; i_be <= utolso; i_be++, i_fa++)
                    acfa_d.fa_2[i_fa].forward(1, fa_elemek[i_be], 2, acfa_d.fa_2);
            }
        }
    }
    else {
        if (padat->is_float) {
            if (padat->melyik_fa == 1) {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                for (; i_be <= utolso; i_be++, i_fa++)
                    dcfa_f.fa_1[i_fa].forward(1, fa_elemek[i_be], 1, dcfa_f.fa_1); // egyelõre az unsafe helyett maradjon az [], hátha nem jó
            }
            else {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                for (; i_be <= utolso; i_be++, i_fa++)
                    dcfa_f.fa_2[i_fa].forward(1, fa_elemek[i_be], 2, dcfa_f.fa_2);
            }
        }
        else {
            if (padat->melyik_fa == 1) {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                for (; i_be <= utolso; i_be++, i_fa++)
                    dcfa_d.fa_1[i_fa].forward(1, fa_elemek[i_be], 1, dcfa_d.fa_1);
            }
            else {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                for (; i_be <= utolso; i_be++, i_fa++)
                    dcfa_d.fa_2[i_fa].forward(1, fa_elemek[i_be], 2, dcfa_d.fa_2);
            }
        }
    }
//printf("\nfw klaszter vege\n");
}


//***********************************************************************
void run_fw_egyedi(szal_feladat_adatai * padat) {
//***********************************************************************
//printf("\nfw egyedi start\n");
    const vektor<adat_fa_elem> & fa_elemek = akt_sim.p_akt_sim->fa_elemek;
    if (padat->is_ac) {
        if (padat->is_float) {
            if (padat->melyik_fa == 1) {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                acfa_f.fa_1[i_fa].forward(padat->al_szalszam, fa_elemek[i_be], 1, acfa_f.fa_1); // egyelõre az unsafe helyett maradjon az [], hátha nem jó
            }
            else {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                acfa_f.fa_2[i_fa].forward(padat->al_szalszam, fa_elemek[i_be], 2, acfa_f.fa_2);
            }
        }
        else {
            if (padat->melyik_fa == 1) {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                acfa_d.fa_1[i_fa].forward(padat->al_szalszam, fa_elemek[i_be], 1, acfa_d.fa_1);
            }
            else {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                acfa_d.fa_2[i_fa].forward(padat->al_szalszam, fa_elemek[i_be], 2, acfa_d.fa_2);
            }
        }
    }
    else {
        if (padat->is_float) {
            if (padat->melyik_fa == 1) {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                dcfa_f.fa_1[i_fa].forward(padat->al_szalszam, fa_elemek[i_be], 1, dcfa_f.fa_1); // egyelõre az unsafe helyett maradjon az [], hátha nem jó
            }
            else {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                dcfa_f.fa_2[i_fa].forward(padat->al_szalszam, fa_elemek[i_be], 2, dcfa_f.fa_2);
            }
        }
        else {
            if (padat->melyik_fa == 1) {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                dcfa_d.fa_1[i_fa].forward(padat->al_szalszam, fa_elemek[i_be], 1, dcfa_d.fa_1);
            }
            else {
                uns i_be = padat->klaszter_kezdoindex, i_fa = padat->faelem_kezdoindex;
                dcfa_d.fa_2[i_fa].forward(padat->al_szalszam, fa_elemek[i_be], 2, dcfa_d.fa_2);
            }
        }
    }
//printf("\nfw egyedi vege\n");
}


//***********************************************************************
void run_bw_egyedi(szal_feladat_adatai * padat) {
//***********************************************************************
    const vektor<adat_fa_elem> & fa_elemek = akt_sim.p_akt_sim->fa_elemek;
    if (padat->is_ac) {
        if (padat->is_float) {
            if (padat->melyik_fa == 1) acfa_f.fa_1[padat->faelem_utolso_index].backward(padat->al_szalszam); // egyelõre az unsafe helyett maradjon az [], hátha nem jó
            else                       acfa_f.fa_2[padat->faelem_utolso_index].backward(padat->al_szalszam);
        }
        else {
            if (padat->melyik_fa == 1) acfa_d.fa_1[padat->faelem_utolso_index].backward(padat->al_szalszam);
            else                       acfa_d.fa_2[padat->faelem_utolso_index].backward(padat->al_szalszam);
        }
    }
    else {
        if (padat->is_float) {
            if (padat->melyik_fa == 1) dcfa_f.fa_1[padat->faelem_utolso_index].backward(padat->al_szalszam); // egyelõre az unsafe helyett maradjon az [], hátha nem jó
            else                       dcfa_f.fa_2[padat->faelem_utolso_index].backward(padat->al_szalszam);
        }
        else {
            if (padat->melyik_fa == 1) dcfa_d.fa_1[padat->faelem_utolso_index].backward(padat->al_szalszam);
            else                       dcfa_d.fa_2[padat->faelem_utolso_index].backward(padat->al_szalszam);
        }
    }
}


//***********************************************************************
void run_bw_klaszter(szal_feladat_adatai * padat) {
// A klaszter mátrix szálszáma mindig 1
//***********************************************************************
    const uns elso = padat->klaszter_kezdoindex;
    const vektor<adat_fa_elem> & fa_elemek = akt_sim.p_akt_sim->fa_elemek;
    if (padat->is_ac) {
        if (padat->is_float) {
            if (padat->melyik_fa == 1) {
                uns i_be = padat->klaszter_utolso_index, i_fa = padat->faelem_utolso_index;
                for (; i_be >= elso; i_be--, i_fa--) // kihasználjuk, hogy a klaszter kezdõindexe nem lehet 1-nél kisebb
                    acfa_f.fa_1.unsafe(i_fa).backward(1);
            }
            else {
                uns i_be = padat->klaszter_utolso_index, i_fa = padat->faelem_utolso_index;
                for (; i_be >= elso; i_be--, i_fa--)
                    acfa_f.fa_2.unsafe(i_fa).backward(1);
            }
        }
        else {
            if (padat->melyik_fa == 1) {
                uns i_be = padat->klaszter_utolso_index, i_fa = padat->faelem_utolso_index;
                for (; i_be >= elso; i_be--, i_fa--)
                    acfa_d.fa_1.unsafe(i_fa).backward(1);
            }
            else {
                uns i_be = padat->klaszter_utolso_index, i_fa = padat->faelem_utolso_index;
                for (; i_be >= elso; i_be--, i_fa--)
                    acfa_d.fa_2.unsafe(i_fa).backward(1);
            }
        }
    }
    else {
        if (padat->is_float) {
            if (padat->melyik_fa == 1) {
                uns i_be = padat->klaszter_utolso_index, i_fa = padat->faelem_utolso_index;
                for (; i_be >= elso; i_be--, i_fa--)
                    dcfa_f.fa_1.unsafe(i_fa).backward(1);
            }
            else {
                uns i_be = padat->klaszter_utolso_index, i_fa = padat->faelem_utolso_index;
                for (; i_be >= elso; i_be--, i_fa--)
                    dcfa_f.fa_2.unsafe(i_fa).backward(1);
            }
        }
        else {
            if (padat->melyik_fa == 1) {
                uns i_be = padat->klaszter_utolso_index, i_fa = padat->faelem_utolso_index;
                for (; i_be >= elso; i_be--, i_fa--)
                    dcfa_d.fa_1.unsafe(i_fa).backward(1);
            }
            else {
                uns i_be = padat->klaszter_utolso_index, i_fa = padat->faelem_utolso_index;
                for (; i_be >= elso; i_be--, i_fa--)
                    dcfa_d.fa_2.unsafe(i_fa).backward(1);
            }
        }
    }
}


//***********************************************************************
void run_matrixmuvelet(szal_feladat_adatai * padat) {
//***********************************************************************
    switch (padat->matrix_muvelet) {
        case szal_feladat_adatai::mmt_dc_f: padat->matrix_adat_dc_f.run_job();  break;
        case szal_feladat_adatai::mmt_dc_d: padat->matrix_adat_dc_d.run_job();  break;
        case szal_feladat_adatai::mmt_ac_f: padat->matrix_adat_ac_f.run_job();  break;
        case szal_feladat_adatai::mmt_ac_d: padat->matrix_adat_ac_d.run_job();  break;
        default: throw hiba("run_matrixmuvelet", "unknown matrix_muvelet type");
    }
}


//***********************************************************************
void run_cellafeldolgozo_klaszter(szal_feladat_adatai * padat) {
//***********************************************************************
    if (padat->is_ac) {
        TODO("run_cellafeldolgozo_klaszter: AC");
    }
    else { // dc + transi

        // elõkészítés, ha ez ûj lépés

        rvt max_abs_Uc = rvt(), max_abs_Tc = rvt();
        if (padat->altipus == szat_pre) {
            if (akt_sim.is_uj_cellaszerkezet_kell) {
                for (uns i = padat->klaszter_kezdoindex; i <= padat->klaszter_utolso_index; i++) {
                    os_cella ** pp_cella = &cellak.unsafe(i);
                    if (*pp_cella == nullptr) { // létrehozzuk a cellát
                        if (akt_sim.p_akt_sim->cellak.unsafe(i).cella_tipus == ct_faces_cella) {
                            *pp_cella = new faces_cella;
                        }
                        else
                            TODO("run_cellafeldolgozo_klaszter => new");
                    }
                    (*pp_cella)->pre_init(i);
                }
                for (uns i = padat->klaszter_kezdoindex; i <= padat->klaszter_utolso_index; i++) {
                    os_cella & akt_cella = *cellak.unsafe(i);
                    akt_cella.init_dc(); // nem kell törölni a prev-et, mert a vezérlõ set_size-ot hívott, ami mindent törölt. Hívja a resize matrixokat.
                    akt_cella.halozat_foglalas_cellaban_dc();
                    akt_cella.peremfeltetel_update_dc();
                    akt_cella.gerj_update_dc();
                    akt_cella.update_for_uj_lepes_dc();
                    akt_cella.set_is_szimm();
                }
            }
            else { // ha nem kell új cellaszerkezet, de új lépés
                for (uns i = padat->klaszter_kezdoindex; i <= padat->klaszter_utolso_index; i++) {
                    os_cella & akt_cella = *cellak.unsafe(i);
                    if (akt_sim.is_uj_fa_kell)
                        akt_cella.resize_matrixok_dc();
                    if (akt_sim.is_del_all_prev || akt_sim.is_tamb_update_kell) {
                        akt_cella.del_all_prev_dc(); // dc-t is és ac-t is töröljük
                        akt_cella.del_all_prev_ac();
                    }
                    if (akt_sim.is_uj_facek_letrehozasa_kell || akt_sim.is_force_face_update) {
                        if (akt_sim.is_uj_facek_letrehozasa_kell) // új anyag, stb.
                            akt_cella.reset_H();
                        akt_cella.halozat_foglalas_cellaban_dc();
                        akt_cella.peremfeltetel_update_dc();
                        akt_cella.gerj_update_dc();
                        akt_cella.set_is_szimm();
                    }
                    else {
                        if (akt_sim.is_peremfelt_update_kell)
                            akt_cella.peremfeltetel_update_dc();
                        if (akt_sim.is_gerj_update_kell) {
                            akt_cella.gerj_update_dc();
                            akt_cella.set_is_szimm();
                        }
                    }
                    akt_cella.update_for_uj_lepes_dc();
                }
            }
        }
        else { // post típusú cellafeldolgozás: feszültségek/hõmérsékletek számítása, UT hibák számítása, dt frissítés, ha kell
            if (!akt_sim.is_iter_dt_csokkento) {
                if (akt_sim.p_akt_sim->nemlin_tipus == nmt_klasszik_iteracio) {
                    rvt UT_hiba = rvt(), U_hiba = rvt(), T_hiba = rvt();
                    for (uns i = padat->klaszter_kezdoindex; i <= padat->klaszter_utolso_index; i++) {
                        os_cella & akt_cella = *cellak.unsafe(i);
                        rvt akt_T_hiba = rvt();
                        rvt akt_hiba = akt_cella.cella_klaszter_2_SA_UT_dc(akt_T_hiba);
                        rvt akt_U_hiba = akt_hiba - akt_T_hiba;
                        if (akt_hiba > UT_hiba)
                            UT_hiba = akt_hiba;
                        if (akt_T_hiba > T_hiba)
                            T_hiba = akt_T_hiba;
                        if (akt_U_hiba > U_hiba)
                            U_hiba = akt_U_hiba;
                    }
                    padat->sumhiba.max_UT_hiba = UT_hiba;
                    padat->sumhiba.max_T_hiba = T_hiba;
                    padat->sumhiba.max_U_hiba = U_hiba;
                }
                else {
                    rvt UT_hiba = rvt(), U_hiba = rvt(), T_hiba = rvt();
                    for (uns i = padat->klaszter_kezdoindex; i <= padat->klaszter_utolso_index; i++) {
                        os_cella & akt_cella = *cellak.unsafe(i);
                        rvt Uc, Tc, akt_T_hiba = rvt();
                        rvt akt_hiba = akt_cella.cella_klaszter_2_NR_UT_dc(Uc, Tc, akt_T_hiba);
                        rvt akt_U_hiba = akt_hiba - akt_T_hiba;
                        if (akt_hiba > UT_hiba)
                            UT_hiba = akt_hiba;
                        if (akt_T_hiba > T_hiba)
                            T_hiba = akt_T_hiba;
                        if (akt_U_hiba > U_hiba)
                            U_hiba = akt_U_hiba;
                        Uc = abs(Uc);
                        Tc = abs(Tc);
                        if (Uc > max_abs_Uc)
                            max_abs_Uc = Uc;
                        if (Tc > max_abs_Tc)
                            max_abs_Tc = Tc;
                    }
                    padat->sumhiba.max_UT_hiba = UT_hiba;
                    padat->sumhiba.max_T_hiba = T_hiba;
                    padat->sumhiba.max_U_hiba = U_hiba;
                    padat->sumhiba.max_abs_Uc = max_abs_Uc;
                    padat->sumhiba.max_abs_Tc = max_abs_Tc;
                }
            }
            if (akt_sim.is_subiter_dt_changed) {
                for (uns i = padat->klaszter_kezdoindex; i <= padat->klaszter_utolso_index; i++) {
                    os_cella & akt_cella = *cellak.unsafe(i);
                    akt_cella.update_dt_dc();
                    akt_cella.update_for_uj_lepes_dc();
                }
            }
        }

        // paraméterek frissítése

        for (uns i = padat->klaszter_kezdoindex; i <= padat->klaszter_utolso_index; i++) {
            os_cella & akt_cella = *cellak.unsafe(i);
            akt_cella.cella_klaszter_3_parameterek_frissitese_dc();
            //              if (i == 22493)
            //                  cellak.unsafe(i).print_G();
        }

        // ágáramok + dissz számítása (szum dissz nullázás az elõzõ lépésben)

        for (uns i = padat->klaszter_kezdoindex; i <= padat->klaszter_utolso_index; i++) {
            os_cella & akt_cella = *cellak.unsafe(i);
            akt_cella.cella_klaszter_4_agaramok_szamitasa_dc();
        }

        // hibaáramok számítása, J feltöltése

        if (akt_sim.p_akt_sim->nemlin_tipus == nmt_klasszik_iteracio) {
            for (uns i = padat->klaszter_kezdoindex; i <= padat->klaszter_utolso_index; i++) {
                os_cella & akt_cella = *cellak.unsafe(i);
                akt_cella.cella_klaszter_5_SA_fill_J_dc();
            }
        }
        else {
            rvt IP_hiba = rvt();
            rvt sum_egyik = rvt(), sum_masik = rvt();
            for (uns i = padat->klaszter_kezdoindex; i <= padat->klaszter_utolso_index; i++) {
                os_cella & akt_cella = *cellak.unsafe(i);
                rvt akt_hiba = akt_cella.cella_klaszter_5_NR_hibaaramok_szamitasa_dc();
                IP_hiba += akt_hiba*akt_hiba;
                //                  if (akt_hiba > 1e-14)
                //                      printf("%u: %g\n", i, akt_hiba);
            }
            padat->sumhiba.sum_IP_hiba = IP_hiba;
        } 
    }
}


//***********************************************************************
void run_lepo_klaszter(szal_feladat_adatai * padat) {
//***********************************************************************
    switch (padat->altipus) {
        case szat_elore:
            emlekek.leptet(padat->klaszter_kezdoindex, padat->klaszter_utolso_index);
            break;
        case szat_vissza:
            emlekek.visszalep(padat->klaszter_kezdoindex, padat->klaszter_utolso_index);
            break;
        case szat_megtart_akt:
            emlekek.megtartando_az_aktualis(padat->klaszter_kezdoindex, padat->klaszter_utolso_index);
            break;
        case szat_ment_akt:
            emlekek.mentendo_az_aktualis(padat->klaszter_kezdoindex, padat->klaszter_utolso_index);
            break;
        case szat_kiind_akt:
            emlekek.kiindulo_az_aktualis(padat->klaszter_kezdoindex, padat->klaszter_utolso_index);
            break;
        case szat_kiind_akt_es_elore:
            emlekek.kiindulo_az_aktualis_es_leptet(padat->klaszter_kezdoindex, padat->klaszter_utolso_index);
            break;
        default:
            throw hiba("run_lepo_klaszter", "Unknown altipus (%u)", padat->altipus);
    }
}


//***********************************************************************
struct egy_adat_dc {
//***********************************************************************
    bool is_el, is_th;
    float el_adat, th_adat;
    float diss, blue, yellow, lum; // sum-nál a kimenõ van a blue-ban és a yellow-ban, és a lum-ban a bemenõ blue
    egy_adat_dc() :is_el{ false }, is_th{ false }, el_adat{ 0.0f }, th_adat{ 0.0f },
        diss{ 0 }, blue{ 0 }, yellow{ 0 }, lum{ 0 } {}
    void set_el(rvt adat) { is_el = true; el_adat = (float)adat; }
    void set_th(rvt adat) { is_th = true; th_adat = (float)adat; }
};


//***********************************************************************
struct egy_cella_dc {
//***********************************************************************
    egy_adat_dc center;
    vektor<float> face_adatok; // itt csak az érték, a kiolvasónak kell tudnia
};


//***********************************************************************
struct mentendo_adatok {
// TODO: AC
//***********************************************************************
    eredmeny_tipus tipus; // et_c_pontprobe, et_f_pontprobe, et_c_map, et_f_map
    mit_ment_tipus mit_ment; // mm_UT, mm_IP, mm_RL
    analizis_lepes_tipus analizis_tipus; // alt_dc, alt_trans, alt_ac
    float analizis_value; // idõpont vagy frekvencia
    egy_adat_dc probe_dc;
    vektor<egy_adat_dc> c_map_dc; // középpont mentésekor
    vektor<egy_cella_dc> f_map_dc; // minden face adat mentésekor
    mentendo_adatok() {}
    vektor<float> nyers_map;
    void set_nyers_map();
private:
    uns calc_map_size();
};


//***********************************************************************
uns mentendo_adatok::calc_map_size() {
//***********************************************************************
    uns hossz = 0;
    if (tipus == et_c_map) {
        if (mit_ment == mm_UT || mit_ment == mm_IP) {
            for (uns i = 1; i < c_map_dc.size(); i++) {
                const egy_adat_dc & akt_adat = c_map_dc.unsafe(i);
                if (akt_adat.is_el)
                    hossz++;
                if (akt_adat.is_th)
                    hossz++;
            }
        }
        else if (mit_ment == mm_RL) { // lustaságból a csak termikus cellákhoz is kiírja a rad/lum értéket. majd egyszer éredemes lenne optimalizálni
            for (uns i = 1; i < c_map_dc.size(); i++) { // a c_map már eleve 6*-os méretû
                const egy_adat_dc & akt_adat = c_map_dc.unsafe(i);
                if (akt_adat.is_el)
                    hossz++;
                if (akt_adat.is_th)
                    hossz++;
            }
        }
/*        else if (mit_ment == mm_RL) { // csak az elektromos cellákhoz ír rad/lum értéket. igazából még ezen is lehetne optimalizálni, mert elég lenne csak a junctont tartalmazó cellákat kiírni.
            for (uns i = 1; i < c_map_dc.size(); i++) {
                const egy_adat_dc & akt_adat = c_map_dc.unsafe(i);
                if (akt_adat.is_el)
                    hossz += 2;
            }
        }
*/
    }
    else if (tipus == et_f_map) {
        if (mit_ment == mm_UT || mit_ment == mm_IP) {
            for (uns i = 1; i < f_map_dc.size(); i++) {
                const egy_cella_dc & akt_cella = f_map_dc.unsafe(i);
                if (akt_cella.center.is_el)
                    hossz++;
                if (akt_cella.center.is_th)
                    hossz++;
                hossz += akt_cella.face_adatok.size();
            }
        }
        else if (mit_ment == mm_RL) {
            TODO("mentendo_adatok::calc_map_size() => mit_ment == mm_RL");
        }
/*        else if (mit_ment == mm_RL) {
            for (uns i = 1; i < f_map_dc.size(); i++) {
                const egy_cella_dc & akt_cella = f_map_dc.unsafe(i);
                if (akt_cella.center.is_el)
                    hossz += 2;
                hossz += akt_cella.face_adatok.size();
            }
        }
*/
    }
    return hossz;
}


//***********************************************************************
void mentendo_adatok::set_nyers_map() {
//***********************************************************************
    nyers_map.set_size(calc_map_size());
    uns hossz = 0;
    if (tipus == et_c_map) {
        if (mit_ment == mm_UT || mit_ment == mm_IP) {
            for (uns i = 1; i < c_map_dc.size(); i++) {
                const egy_adat_dc & akt_adat = c_map_dc.unsafe(i);
                if (akt_adat.is_el)
                    nyers_map.unsafe(hossz++) = akt_adat.el_adat;
                if (akt_adat.is_th)
                    nyers_map.unsafe(hossz++) = akt_adat.th_adat;
            }
        }
        else if (mit_ment == mm_RL) {
            for (uns i = 1; i < c_map_dc.size(); i++) { // a c_map-ben egymás után jön a rada, radb, luma, lumb, dissa, dissb érték, de mindegyik elektromos-termikus párként
                const egy_adat_dc & akt_adat = c_map_dc.unsafe(i);
                if (akt_adat.is_el)
                    nyers_map.unsafe(hossz++) = akt_adat.el_adat;
                if (akt_adat.is_th) // a c_map a termikus cellákra 0-t tartalmaz, de hogy ne kelljen külön leíró fájl, ezt is kiírjuk
                    nyers_map.unsafe(hossz++) = akt_adat.th_adat;
            }
        }
    }
    else if (tipus == et_f_map) {
        if (mit_ment == mm_UT || mit_ment == mm_IP) {
            for (uns i = 1; i < f_map_dc.size(); i++) {
                const egy_cella_dc & akt_cella = f_map_dc.unsafe(i);
                if (akt_cella.center.is_el)
                    nyers_map.unsafe(hossz++) = akt_cella.center.el_adat;
                if (akt_cella.center.is_th)
                    nyers_map.unsafe(hossz++) = akt_cella.center.th_adat;
                for (uns j = 1; j < akt_cella.face_adatok.size(); j++)
                    nyers_map.unsafe(hossz++) = akt_cella.face_adatok.unsafe(j);
            }
        }
        else if (mit_ment == mm_RL) {
            TODO("mentendo_adatok::set_nyers_map CRL map");
        }
    }
}


//***********************************************************************
void run_belso_mentes(szal_feladat_adatai * padat) {
//***********************************************************************
    set_hiba_hol h("run_belso_mentes");
    if (padat->p_akt_eredm == nullptr)
        return;
    
    szal_feladat_adatai fajlba_mentes;
    fajlba_mentes.beallitas_fajlba_mentesnek();
    fajlba_mentes.p_mentendo_adatok = new lista<mentendo_adatok>; // fel kell majd szabadítani
    fajlba_mentes.eredm_utvonal = padat->eredm_utvonal;
    fajlba_mentes.analizis_tipus = padat->analizis_tipus;
    fajlba_mentes.analizis_value = padat->analizis_value;
    fajlba_mentes.is_use_commas = padat->is_use_commas;
    fajlba_mentes.akt_anal_index = padat->akt_anal_index;
    rvt Tamb = padat->Tamb;

    for (uns i = 0; i < padat->p_akt_eredm->size(); i++) {
        mentendo_adatok & akt = fajlba_mentes.p_mentendo_adatok->push_back();
        akt.tipus = padat->p_akt_eredm->unsafe(i).tipus; // pontprobe/imagemap
        akt.mit_ment = padat->p_akt_eredm->unsafe(i).mit_ment;
        akt.analizis_tipus = padat->analizis_tipus;
        akt.analizis_value = (float)padat->analizis_value; // tranziensnél nem dt, hanem teljes idõ
        if (akt.analizis_tipus == alt_ac) {
            TODO("run_belso_mentes: AC");
        }
        else{
            switch (akt.tipus) {
                case et_c_pontprobe: {
                    if (cellak.unsafe(padat->p_akt_eredm->unsafe(i).cella_index)->get_tipus() == ct_faces_cella) {
                        const faces_cella & akt_cella = *static_cast<faces_cella*>(cellak[padat->p_akt_eredm->unsafe(i).cella_index]);
                        switch (akt.mit_ment) {
                            case mm_UT: {
                                if (akt_cella.is_el)
                                    akt.probe_dc.set_el(akt_cella.el_center_face_dc.emlekek.get_mentendo().UT);
                                if (akt_cella.is_th)
                                    akt.probe_dc.set_th(akt_cella.th_center_face_dc.emlekek.get_mentendo().UT + Tamb);
                                float A_top = (float)(akt_cella.A_top * 4 * 3.14159265358);
                                akt.probe_dc.diss   = (float)akt_cella.emlekek.get_mentendo().sum_diss;
                                akt.probe_dc.blue   = (float)akt_cella.emlekek.get_mentendo().sum_blue / A_top;
                                akt.probe_dc.yellow = (float)akt_cella.emlekek.get_mentendo().sum_yellow / A_top;
                                akt.probe_dc.lum    = (float)akt_cella.emlekek.get_mentendo().sum_lum / A_top;
                            }
                            break;
                            case mm_IP: {
                                if (akt_cella.is_el)
                                    akt.probe_dc.set_el(akt_cella.el_center_face_dc.emlekek.get_mentendo().IP);
                                if (akt_cella.is_th)
                                    akt.probe_dc.set_th(akt_cella.th_center_face_dc.emlekek.get_mentendo().IP);
                            }
                            break;
                            default:
                                throw hiba(1, "unsupported probe type");
                            break;
                        }
                    }
                    else
                        TODO("run_belso_mentes => 1");
                }
                break;
                case et_f_pontprobe: {
                    if (cellak.unsafe(padat->p_akt_eredm->unsafe(i).cella_index)->get_tipus() == ct_faces_cella) {
                        const faces_cella & akt_cella = *static_cast<faces_cella*>(cellak[padat->p_akt_eredm->unsafe(i).cella_index]);
                        const face_tok_dc & akt_tok = akt_cella.faces_dc[padat->p_akt_eredm->unsafe(i).face_index];
                        switch (akt.mit_ment) {
                            case mm_UT: {
                                if (akt_tok.p_core->is_El())
                                    akt.probe_dc.set_el(akt_tok.emlekek.get_mentendo().UT);
                                else
                                    akt.probe_dc.set_th(akt_tok.emlekek.get_mentendo().UT + Tamb);
                            }
                            break;
                            case mm_IP: {
                                if (akt_tok.p_core->is_El())
                                    akt.probe_dc.set_el(akt_tok.emlekek.get_mentendo().IP);
                                else
                                    akt.probe_dc.set_th(akt_tok.emlekek.get_mentendo().IP);
                            }
                            break;
                            default:
                                throw hiba(1, "unsupported probe type");
                            break;
                        }
                    }
                    else
                        TODO("run_belso_mentes => 2");
                }
                break;
                case et_c_map: {
                    switch (akt.mit_ment) {
                        case mm_UT: {
                            akt.c_map_dc.set_size(cellak.size());
                            for (uns i = 1; i < cellak.size(); i++) {
                                if (cellak.unsafe(i)->get_tipus() == ct_faces_cella) {
                                    const faces_cella & akt_cella = *static_cast<faces_cella*>(cellak.unsafe(i));
                                    egy_adat_dc & akt_adat = akt.c_map_dc.unsafe(i);
                                    if (akt_cella.is_el)
                                        akt_adat.set_el(akt_cella.el_center_face_dc.emlekek.get_mentendo().UT);
                                    if (akt_cella.is_th)
                                        akt_adat.set_th(akt_cella.th_center_face_dc.emlekek.get_mentendo().UT + Tamb);
                                }
                                else
                                    TODO("run_belso_mentes => 3");
                            }
                        }
                        break;
                        case mm_IP: {
                            akt.c_map_dc.set_size(cellak.size());
                            for (uns i = 1; i < cellak.size(); i++) {
                                if (cellak.unsafe(i)->get_tipus() == ct_faces_cella) {
                                    const faces_cella & akt_cella = *static_cast<faces_cella*>(cellak.unsafe(i));
                                    egy_adat_dc & akt_adat = akt.c_map_dc.unsafe(i);
                                    if (akt_cella.is_el)
                                        akt_adat.set_el(akt_cella.el_center_face_dc.emlekek.get_mentendo().IP);
                                    if (akt_cella.is_th)
                                        akt_adat.set_th(akt_cella.th_center_face_dc.emlekek.get_mentendo().IP);
                                }
                                else
                                    TODO("run_belso_mentes => 4");
                            }
                        }
                        break;
                        case mm_RL: {
                            akt.c_map_dc.set_size(cellak.size()*6); // 6*-es méret a 6 féle adatnak
                            //printf("//////////////// %u\n", akt.c_map_dc.size());
                            //rvt sum_R = 0, sum_A = 0, sum_rad = 0;
                            for (uns i = 1; i < cellak.size(); i++) {
                                if (cellak.unsafe(i)->get_tipus() == ct_faces_cella) {
                                    const faces_cella & akt_cella = *static_cast<faces_cella*>(cellak.unsafe(i));
                                    rvt A_LED = akt_cella.emlekek.get_mentendo().A_LED * 4 * rvt(3.14159265358); // !!! meg van szorozva 4 pi-vel
                                    rvt V_cella = akt_sim.p_akt_sim->meretek[akt_cella.p_adatcella->volume_index];
                                    rvt A_top = akt_cella.A_top * 4 * rvt(3.14159265358); // !!! meg van szorozva 4 pi-vel
                                    egy_adat_dc & akt_adat_R_telj = akt.c_map_dc.unsafe(6 * i);
                                    egy_adat_dc & akt_adat_R_suruseg = akt.c_map_dc.unsafe(6 * i + 1);
                                    egy_adat_dc & akt_adat_L_telj = akt.c_map_dc.unsafe(6 * i + 2);
                                    egy_adat_dc & akt_adat_L_suruseg = akt.c_map_dc.unsafe(6 * i + 3);
                                    egy_adat_dc & akt_adat_dissz_telj = akt.c_map_dc.unsafe(6 * i + 4);
                                    egy_adat_dc & akt_adat_dissz_suruseg = akt.c_map_dc.unsafe(6 * i + 5);
                                    rvt dissz = akt_cella.emlekek.get_mentendo().sum_diss;
                                    if (A_LED > 0) {
                                        if (akt_cella.is_el || akt_cella.p_cella_anyag->is_fenypor) {
                                            rvt rad = akt_sim.p_akt_sim->fenyutak.is_write_output ? akt_cella.emlekek.get_mentendo().sum_blue_out : akt_cella.emlekek.get_mentendo().sum_rad;
                                            rvt lum = akt_sim.p_akt_sim->fenyutak.is_write_output ? akt_cella.emlekek.get_mentendo().sum_yellow_out : akt_cella.emlekek.get_mentendo().sum_lum;
                                            akt_adat_R_telj.set_el(rad);
                                            akt_adat_R_suruseg.set_el(rad / A_LED);
                                            akt_adat_L_telj.set_el(lum);
                                            akt_adat_L_suruseg.set_el(lum / A_LED);
                                            akt_adat_dissz_telj.set_el(dissz);
                                            akt_adat_dissz_suruseg.set_el(dissz / V_cella);
                                        }
                                        if (akt_cella.is_th) {
                                            rvt blue = akt_cella.emlekek.get_mentendo().sum_blue;
                                            rvt yllw = akt_cella.emlekek.get_mentendo().sum_yellow;
                                            akt_adat_R_telj.set_th(blue);
                                            akt_adat_R_suruseg.set_th(blue / A_top);
                                            akt_adat_L_telj.set_th(yllw);
                                            akt_adat_L_suruseg.set_th(yllw / A_top);
                                            akt_adat_dissz_telj.set_th(dissz);
                                            akt_adat_dissz_suruseg.set_th(dissz / V_cella);
                                        }
                                    }
                                    else { // 0-kat mentünk, ha nem OUT fénypor üzemmód van
                                        if (akt_cella.is_el || akt_cella.p_cella_anyag->is_fenypor) {
                                            rvt blue = akt_sim.p_akt_sim->fenyutak.is_write_output ? akt_cella.emlekek.get_mentendo().sum_blue_out : 0;
                                            rvt yllw = akt_sim.p_akt_sim->fenyutak.is_write_output ? akt_cella.emlekek.get_mentendo().sum_yellow_out : 0;
                                            akt_adat_R_telj.set_el(blue);
                                            akt_adat_R_suruseg.set_el(blue / A_top);
                                            akt_adat_L_telj.set_el(yllw);
                                            akt_adat_L_suruseg.set_el(yllw / A_top);
                                            akt_adat_dissz_telj.set_el(dissz);
                                            akt_adat_dissz_suruseg.set_el(dissz / V_cella);
                                            //rvt by = blue + yllw;
                                            //sum_R += by;
                                            //sum_A += by != 0 ? A_top : 0;
                                            //sum_rad += by / A_top;
                                            //if (by != 0)
                                            //    printf("A = %g, A*4pi = %g, P = %g, rad = %g\n", akt_cella.A_top, A_top, by, by / A_top);
                                        }
                                        if (akt_cella.is_th) {
                                            rvt blue = akt_cella.emlekek.get_mentendo().sum_blue;
                                            rvt yllw = akt_cella.emlekek.get_mentendo().sum_yellow;
                                            akt_adat_R_telj.set_th(blue);
                                            akt_adat_R_suruseg.set_th(blue / A_top);
                                            akt_adat_L_telj.set_th(yllw);
                                            akt_adat_L_suruseg.set_th(yllw / A_top);
                                            akt_adat_dissz_telj.set_th(dissz);
                                            akt_adat_dissz_suruseg.set_th(dissz / V_cella);
                                        }
                                    }
                                }
                                else
                                    TODO("run_belso_mentes => 5");
                            }
                            //printf("\n\nsumR = %g, sum_A = %g, sum_rad = %g\n\n", sum_R, sum_A, sum_rad);
                        }
                        break;
                        default:
                            throw hiba(1, "unsupported map type");
                        break;
                    }
                }
                break;
                case et_f_map: {
                    akt.f_map_dc.set_size(cellak.size());
                    switch (akt.mit_ment) {
                        case mm_UT: {
                            for (uns i = 1; i < cellak.size(); i++) {
                                if (cellak.unsafe(i)->get_tipus() == ct_faces_cella) {
                                    const faces_cella & akt_cella = *static_cast<faces_cella*>(cellak.unsafe(i));
                                    egy_cella_dc & akt_ki = akt.f_map_dc.unsafe(i);
                                    if (akt_cella.is_el)
                                        akt_ki.center.set_el(akt_cella.el_center_face_dc.emlekek.get_mentendo().UT);
                                    if (akt_cella.is_th)
                                        akt_ki.center.set_th(akt_cella.th_center_face_dc.emlekek.get_mentendo().UT + Tamb);
                                    akt_ki.face_adatok.set_size(akt_cella.faces_dc.size());
                                    for (uns j = 0; j < akt_cella.faces_dc.size(); j++)
                                        akt_ki.face_adatok.unsafe(j) = (float)(akt_cella.faces_dc.unsafe(j).emlekek.get_mentendo().UT + akt_cella.faces_dc.unsafe(j).p_core->is_El() ? rvt() : +Tamb);
                                }
                                else
                                    TODO("run_belso_mentes => 6");
                            }
                        }
                        break;
                        case mm_IP: {
                            for (uns i = 1; i < cellak.size(); i++) {
                                if (cellak.unsafe(i)->get_tipus() == ct_faces_cella) {
                                    const faces_cella & akt_cella = *static_cast<faces_cella*>(cellak.unsafe(i));
                                    egy_cella_dc & akt_ki = akt.f_map_dc.unsafe(i);
                                    if (akt_cella.is_el)
                                        akt_ki.center.set_el(akt_cella.el_center_face_dc.emlekek.get_mentendo().IP);
                                    if (akt_cella.is_th)
                                        akt_ki.center.set_th(akt_cella.th_center_face_dc.emlekek.get_mentendo().IP);
                                    akt_ki.face_adatok.set_size(akt_cella.faces_dc.size());
                                    for (uns j = 0; j < akt_cella.faces_dc.size(); j++)
                                        akt_ki.face_adatok.unsafe(j) = (float)akt_cella.faces_dc.unsafe(j).emlekek.get_mentendo().IP;
                                }
                                else
                                    TODO("run_belso_mentes => 7");
                            }
                        }
                        break;
                        default:
                            throw hiba(1, "unsupported map type");
                        break;
                    }
                }
                break;
            }
        }
    }
    if (padat->analizis_tipus != alt_ac) {
        mentendo_adatok & akt = fajlba_mentes.p_mentendo_adatok->push_back();
        akt.tipus = et_sum_value;
        akt.mit_ment = mm_RL;
        akt.analizis_tipus = padat->analizis_tipus;
        akt.analizis_value = 0;
        rvt sum_be_blue = 0, sum_ki_blue = 0, sum_ki_yellow = 0, sum_ki_diss = 0;
        for (uns i = 1; i < cellak.size(); i++) {
            if (cellak.unsafe(i)->get_tipus() == ct_faces_cella) {
                const faces_cella & akt_cella = *static_cast<faces_cella*>(cellak.unsafe(i));
                for (uns j = 1; j < akt_cella.faces_dc.size(); j++) {
                    sum_be_blue += akt_cella.faces_dc[j].emlekek.get_mentendo().rad;
                }
                sum_ki_blue += akt_cella.emlekek.get_mentendo().sum_blue_out;
                sum_ki_yellow += akt_cella.emlekek.get_mentendo().sum_yellow_out;
                sum_ki_diss += akt_cella.emlekek.get_mentendo().sum_diss_out;
            }
        }
        akt.probe_dc.blue   = (float)sum_ki_blue;
        akt.probe_dc.yellow = (float)sum_ki_yellow;
        akt.probe_dc.lum    = (float)sum_be_blue;
        akt.probe_dc.diss   = (float)sum_ki_diss;
    }
    feldolgozo.betesz(fajlba_mentes, false);
}


//***********************************************************************
void szamkiiro(FILE *fp, float ertek, bool is_use_commas) {
//***********************************************************************
    char s[50];
    sprintf_s(s, 50, "\t%g", ertek);
    if (is_use_commas) {
        char *ppont = strchr(s, '.');
        if (ppont != nullptr)
            *ppont = ',';
    }
    fprintf(fp, s);
}


//***********************************************************************
struct maphead {
//***********************************************************************
    char azon[8];
    uns adatmeret; // bájtban
    uns tomorites_tipus; // 0 = nyers
    maphead() :azon{ "hexmap" }, adatmeret{ 0 }, tomorites_tipus{ 0 } {}
};


//***********************************************************************
void run_fajlba_mentes(szal_feladat_adatai * padat) {
//***********************************************************************
    if (padat->p_mentendo_adatok == nullptr)
        return;

    set_hiba_hol h("run_fajlba_mentes");

    // probe-ok

    uns db = 0;
    for (void * it = padat->p_mentendo_adatok->get_it(); it != nullptr; padat->p_mentendo_adatok->inc_it(it)) {
        const mentendo_adatok & adat = padat->p_mentendo_adatok->get_akt(it);
        if (adat.tipus == et_c_pontprobe || adat.tipus == et_f_pontprobe)
            db++;
    }

    if (db > 0) {
        static ::std::string elozo_utvonal;
        bool is_felulir = elozo_utvonal != padat->eredm_utvonal;
        if (is_felulir) {
            elozo_utvonal = padat->eredm_utvonal;
        }
        FILE *fp;
        if (fopen_s(&fp, (padat->eredm_utvonal + "probe_results.txt").c_str(), is_felulir ? "wt" : "at") != 0)
            throw hiba(1, "cannot open %s to write", (padat->eredm_utvonal + "probe_results.txt").c_str());
        // TODO: fejléc kiírása
        fprintf(fp, "%u", padat->akt_anal_index);
        if (padat->analizis_tipus == alt_dc)
            fprintf(fp, "\tDC");
        else
            szamkiiro(fp, (float)padat->analizis_value, padat->is_use_commas);
        for (void * it = padat->p_mentendo_adatok->get_it(); it != nullptr; padat->p_mentendo_adatok->inc_it(it)) {
            const mentendo_adatok & adat = padat->p_mentendo_adatok->get_akt(it);
            if (adat.tipus == et_c_pontprobe || adat.tipus == et_f_pontprobe) {
                if (adat.probe_dc.is_el)
                    szamkiiro(fp, adat.probe_dc.el_adat, padat->is_use_commas);
                else fprintf(fp, "\t-");
                if (adat.probe_dc.is_th)
                    szamkiiro(fp, adat.probe_dc.th_adat, padat->is_use_commas);
                else fprintf(fp, "\t-");
                szamkiiro(fp, adat.probe_dc.diss, padat->is_use_commas);
                szamkiiro(fp, adat.probe_dc.blue, padat->is_use_commas);
                szamkiiro(fp, adat.probe_dc.yellow, padat->is_use_commas);
                szamkiiro(fp, adat.probe_dc.lum, padat->is_use_commas);
                fprintf(fp, "\n");
            }
            else if (adat.tipus == et_sum_value) {
                if (adat.probe_dc.lum > 0) {
                    fprintf(fp, "\n");
                    fprintf(fp, "blue_in=");
                    szamkiiro(fp, adat.probe_dc.lum, padat->is_use_commas);
                    fprintf(fp, "\nblue_out=");
                    szamkiiro(fp, adat.probe_dc.blue, padat->is_use_commas);
                    fprintf(fp, " out/in=");
                    szamkiiro(fp, adat.probe_dc.blue / adat.probe_dc.lum, padat->is_use_commas);
                    fprintf(fp, "\nyellow_out=");
                    szamkiiro(fp, adat.probe_dc.yellow, padat->is_use_commas);
                    fprintf(fp, " out/in=");
                    szamkiiro(fp, adat.probe_dc.yellow / adat.probe_dc.lum, padat->is_use_commas);
                    fprintf(fp, "\ndiss_out=");
                    szamkiiro(fp, adat.probe_dc.diss, padat->is_use_commas);
                    fprintf(fp, " out/in=");
                    szamkiiro(fp, adat.probe_dc.diss / adat.probe_dc.lum, padat->is_use_commas);
                    fprintf(fp, "\nsum out/in=");
                    szamkiiro(fp, (adat.probe_dc.blue + adat.probe_dc.yellow + adat.probe_dc.diss) / adat.probe_dc.lum, padat->is_use_commas);
                    fprintf(fp, "\n");
                }
            }
        }
        fprintf(fp, "\n");
        fclose(fp);
    }

    // map-ek

    for (void * it = padat->p_mentendo_adatok->get_it(); it != nullptr; padat->p_mentendo_adatok->inc_it(it)) {
        mentendo_adatok & adat = padat->p_mentendo_adatok->get_akt(it);
        if (adat.tipus != et_c_map && adat.tipus != et_f_map)
            continue;
        ::std::string tipusnev = adat.tipus == et_c_map ? "c" : "f";
        switch (adat.mit_ment) {
            case mm_UT: tipusnev += "vt"; break;
            case mm_IP: tipusnev += "ip"; break;
            case mm_RL: tipusnev += "rl"; break;
            default: throw hiba(1, "impossible adat.mit_ment value");
        }
        ::std::string fajlnev = padat->eredm_utvonal + "hexres" + ::std::to_string(padat->akt_anal_index) + '.' + tipusnev;
        adat.set_nyers_map();
        char fej[128] = { 0 };
        maphead & headstrukt = *((maphead*)fej);
        headstrukt = maphead();
        headstrukt.adatmeret = adat.nyers_map.size()*sizeof(float);
        // TODO: rendes fejléc
        FILE *fp;
        if (fopen_s(&fp, fajlnev.c_str(), "wb") != 0)
            throw hiba(1, "cannot open %s to write", (padat->eredm_utvonal + "probe_results.txt").c_str());
        fwrite(fej, 1, 128, fp);
        fwrite(&adat.nyers_map.unsafe(0), sizeof(float), adat.nyers_map.size(), fp); // TODO: jelenleg a nyers map-et menti, de késõbb a tömörítettet
        fclose(fp);
    }

    delete padat->p_mentendo_adatok;
}

//***********************************************************************
void szalfuttato_fuggveny(szal_tipus mit_futtat) {
//***********************************************************************
    try {
        szal_feladat_adatai * padat;
        bool is_torolhetok_a_szalak;
        //printf("%u\n", mit_futtat);
        do {
            while((padat = feldolgozo.kivesz(mit_futtat))==nullptr)
                ;
            //log_print("szal be %u\n", (uns)(padat->tipus));
            is_torolhetok_a_szalak = padat->is_torolhetok_a_szalak;
            if(!is_torolhetok_a_szalak){
                switch (padat->tipus) {
                    case szt_nem_csinal_semmit:                                                 break;
                    case szt_cellafeldolgozo_klaszter:  run_cellafeldolgozo_klaszter(padat);    break;
                    case szt_lepo_klaszter:             run_lepo_klaszter(padat);               break;
                    case szt_fw_klaszter:               run_fw_klaszter(padat);                 break;
                    case szt_fw_egyedi:                 run_fw_egyedi(padat);                   break;
                    case szt_bw_egyedi:                 run_bw_egyedi(padat);                   break;
                    case szt_bw_klaszter:               run_bw_klaszter(padat);                 break;
                    case szt_matrix:                    run_matrixmuvelet(padat);               break;
                    case szt_belso_mentes:              run_belso_mentes(padat);                break;
                    case szt_fajlba_mentes:             run_fajlba_mentes(padat);               break;
                    default:
                        throw hiba("szalfuttato_fuggveny", "unsupported thread job type (%u)", (uns)padat->tipus);
                }
            }
            feldolgozo.kesz_egy_munka(padat);
            feldolgozo.torli_a_keszeket(szt_fajlba_mentes);
            //log_print("szal ki (%u)\n", mit_futtat);
            //feldolgozo.print_status();
        } while (!is_torolhetok_a_szalak);
        szal_feladat_adatai dummy;
        dummy.beallitas_szaltorleshez(mit_futtat);
        feldolgozo.betesz(dummy, false); // betesz mégegyet, hogy a következõ szálat is kinyírja
    }
    catch (hiba h) {
        printf("\n%s\n", h.what());
        abort();
    }
}


//***********************************************************************
inline bool szal_feladat_adatai::is_fw_egyedi_indithato() const {
//***********************************************************************
    const vektor<adat_fa_elem> & fa_elemek = akt_sim.p_akt_sim->fa_elemek;
    if (is_ac) {
        if (is_float) {
            if (melyik_fa == 1) return acfa_f.fa_1[faelem_kezdoindex].is_fw_indithato(fa_elemek[klaszter_kezdoindex], 1, acfa_f.fa_1);
            else                return acfa_f.fa_2[faelem_kezdoindex].is_fw_indithato(fa_elemek[klaszter_kezdoindex], 2, acfa_f.fa_2);
        }
        else {
            if (melyik_fa == 1) return acfa_d.fa_1[faelem_kezdoindex].is_fw_indithato(fa_elemek[klaszter_kezdoindex], 1, acfa_d.fa_1);
            else                return acfa_d.fa_2[faelem_kezdoindex].is_fw_indithato(fa_elemek[klaszter_kezdoindex], 2, acfa_d.fa_2);
        }
    }
    else {
        if (is_float) {
            if (melyik_fa == 1) return dcfa_f.fa_1[faelem_kezdoindex].is_fw_indithato(fa_elemek[klaszter_kezdoindex], 1, dcfa_f.fa_1);
            else                return dcfa_f.fa_2[faelem_kezdoindex].is_fw_indithato(fa_elemek[klaszter_kezdoindex], 2, dcfa_f.fa_2);
        }
        else {
            if (melyik_fa == 1) return dcfa_d.fa_1[faelem_kezdoindex].is_fw_indithato(fa_elemek[klaszter_kezdoindex], 1, dcfa_d.fa_1);
            else                return dcfa_d.fa_2[faelem_kezdoindex].is_fw_indithato(fa_elemek[klaszter_kezdoindex], 2, dcfa_d.fa_2);
        }
    }
}


//***********************************************************************
bool szal_feladat_adatai::is_bw_indithato() const {
//***********************************************************************
    if (is_ac) {
        if (is_float) {
            if (melyik_fa == 1) return acfa_f.fa_1[faelem_utolso_index].is_bw_indithato();
            else                return acfa_f.fa_2[faelem_utolso_index].is_bw_indithato();
        }
        else {
            if (melyik_fa == 1) return acfa_d.fa_1[faelem_utolso_index].is_bw_indithato();
            else                return acfa_d.fa_2[faelem_utolso_index].is_bw_indithato();
        }
    }
    else {
        if (is_float) {
            if (melyik_fa == 1) return dcfa_f.fa_1[faelem_utolso_index].is_bw_indithato();
            else                return dcfa_f.fa_2[faelem_utolso_index].is_bw_indithato();
        }
        else {
            if (melyik_fa == 1) return dcfa_d.fa_1[faelem_utolso_index].is_bw_indithato();
            else                return dcfa_d.fa_2[faelem_utolso_index].is_bw_indithato();
        }
    }
}


//***********************************************************************
bool szal_feladat_adatai::is_indithato() const{
//***********************************************************************
    if (!is_felteteles) // bw egyedi elsõ eleme nem feltételes
        return true;
    switch (tipus) {
        case szt_fw_egyedi: return is_fw_egyedi_indithato();
        case szt_bw_egyedi:
        case szt_bw_klaszter: return is_bw_indithato();
    }
    return true;
}


//***********************************************************************
void Para_engine::kesz_egy_munka(szal_feladat_adatai * p_feladat) {
//***********************************************************************
    ::std::unique_lock<::std::mutex> lock(lezaro_mutex);

    if (p_feladat->p_is_kesz != nullptr)
        *p_feladat->p_is_kesz = true;

    list_elem * kesz_listaelem = szalak[p_feladat->tipus].in_progress_list.pop_elem(p_feladat);
    if (kesz_listaelem == nullptr)
        throw hiba("Para_engine::kesz_egy_munka", "impossibility: kesz_listaelem==nullptr");

    // mivel töröljük p_feladatot, adatait tárolni kell

    ::std::condition_variable * p_kesz_egy_feladat = p_feladat->p_kesz_egy_feladat;
    szal_tipus tipus = p_feladat->tipus; 

    osszes_nincs_kesz_db--;
    in_progress_threads_db--;
    if (is_log_szalak) {
        thread_log_elem log;
        log.thread_id = kesz_listaelem->adatok.thread_id;
        log.in_progress_threads_db = in_progress_threads_db;
        log.free_threads_db = szalak[tipus].waiting_list.getElementNum() + szalak[tipus].waiting_priority_list.getElementNum();
        log.tipus = tipus;
        log.mi_tortenik = thread_log_elem::tlmt_most_van_vege;
        log.matrix_fuggveny_tipus = tipus == szt_matrix ? kesz_listaelem->adatok.matrix_fv_tipus : 0;
        log.ido = std::chrono::system_clock::now();
        thread_log_lista.push_back(log);
    }
    if (szalak[tipus].is_azonnal_torolheto)
        delete kesz_listaelem;
    else
        szalak[tipus].finished_list.push_back(kesz_listaelem);

    lock.unlock();
    szalak[tipus].van_varakozo.notify_all(); // saját típusú indulhat, ha tud
    if (szalak[tipus].ertesitendo_felteteles != szt_N) // más típusú indulhat, ha tud
        szalak[szalak[tipus].ertesitendo_felteteles].van_varakozo.notify_all();
    szalak[tipus].van_futo.notify_all(); // kész ven-e az összes ellenõrzése indulhat
    van_futo.notify_one(); // kész ven-e az összes ellenõrzése indulhat

    // A hívó folyamat értesítése arról, hogy kész a feladat (mátrixmûveletnél használjuk, de lehetne máshol is)

    if (p_kesz_egy_feladat != nullptr) {
        p_kesz_egy_feladat->notify_all();
    }
}


//***********************************************************************
szal_feladat_adatai * Para_engine::kivesz(szal_tipus tipus) {
//***********************************************************************
    ::std::unique_lock<::std::mutex> lock(lezaro_mutex);

    szalak[tipus].van_varakozo.wait(lock, [this, tipus]() {
        bool felt = (!szalak[tipus].is_felteteles || szalak[tipus].waiting_list.is_indithato() || szalak[tipus].waiting_priority_list.is_indithato());
        //printf("felt = %s (%u)\n", felt ? "true" : "false", tipus);
        return (szalak[tipus].waiting_list.getElementNum() + szalak[tipus].waiting_priority_list.getElementNum()) > 0 && felt;
    });

    // Elõször a priority listából veszi ki, ha abból nem lehet, akkor a normálból

    list_elem * kivett_listaelem = szalak[tipus].waiting_priority_list.pop_indithato();
    if (kivett_listaelem == nullptr)
        kivett_listaelem = szalak[tipus].waiting_list.pop_indithato();

    if (kivett_listaelem == nullptr)
        throw hiba("Para_engine::kivesz", "impossibility: kivett_listaelem==nullptr");

    szalak[tipus].in_progress_list.push_back(kivett_listaelem);
    in_progress_threads_db++;
    if (is_log_szalak) {
        thread_log_elem log;
        log.thread_id = kivett_listaelem->adatok.thread_id;
        log.in_progress_threads_db = in_progress_threads_db;
        log.free_threads_db = szalak[tipus].waiting_list.getElementNum() + szalak[tipus].waiting_priority_list.getElementNum();
        log.tipus = tipus;
        log.mi_tortenik = thread_log_elem::tlmt_most_indul;
        log.matrix_fuggveny_tipus = tipus == szt_matrix ? kivett_listaelem->adatok.matrix_fv_tipus : 0;
        log.ido = std::chrono::system_clock::now();
        thread_log_lista.push_back(log);
    }
    //osszes_feldolgozasra_varo_munka_db--;

    //kivett_listaelem->adatok.all = szal_feladat_adatai::a_feldolgozas_alatt;

    lock.unlock();
    return &kivett_listaelem->adatok;
}


//***********************************************************************
void Para_engine::set_szalszam(uns szal_db) {
//***********************************************************************
    //szal_db = 1;
    if (szalak[1].szalak.size() > 0)
        throw hiba("Para_engine::set_szalszam", "this function can be called only once during the program run");
    if (szal_db == 0)
        throw hiba("Para_engine::set_szalszam", "thread number must be > 0");
    for (szal_tipus i = szal_tipus(1); i < szt_N; i = szal_tipus(i + 1)) {
        //uns akt_db = szal_db;
        //if (i == szt_fw_klaszter || i == szt_matrix) {
        //    akt_db = 4 * szal_db;
        //}
        if(szaltipus_jellemzoi_tomb[i].valid != i)
            throw hiba("Para_engine::set_szalszam", "invalid szaltipus_jellemzoi_tomb[%u]", i);
        szalak[i].is_felteteles = szaltipus_jellemzoi_tomb[i].is_felteteles;
        szalak[i].is_azonnal_torolheto = szaltipus_jellemzoi_tomb[i].is_azonnal_torolheto;
        szalak[i].ertesitendo_felteteles = szaltipus_jellemzoi_tomb[i].ertesitendo_felteteles;
        szalak[i].szalak.set_size(szaltipus_jellemzoi_tomb[i].is_egyszalu ? 1 : szal_db);
        for (uns j = 0; j < szalak[i].szalak.size(); j++)
            szalak[i].szalak[j] = ::std::thread{ szalfuttato_fuggveny, i };
    }
}


//***********************************************************************
void Para_engine::betesz(const szal_feladat_adatai & be, bool is_priority) {
//***********************************************************************
    ::std::unique_lock<::std::mutex> lock(lezaro_mutex);

    if (is_priority)
        szalak[be.tipus].waiting_priority_list.push_back(be, akt_thread_id);
    else
        szalak[be.tipus].waiting_list.push_back(be, akt_thread_id);
    akt_thread_id++;
    osszes_nincs_kesz_db++;
    //osszes_feldolgozasra_varo_munka_db++;

    lock.unlock();
    szalak[be.tipus].van_varakozo.notify_all();
}


//***********************************************************************
void Para_engine::var_mig_van_varo_vagy_futo_feladat_barmilyen() {
//***********************************************************************
    ::std::unique_lock<::std::mutex> lock(lezaro_mutex);

    van_futo.wait(lock, [this]() {return osszes_nincs_kesz_db == 0; });
}


//***********************************************************************
void Para_engine::var_mig_van_varo_vagy_futo_feladat(szal_tipus tipus) {
//***********************************************************************
    ::std::unique_lock<::std::mutex> lock(lezaro_mutex);

    szalak[tipus].van_futo.wait(lock, [this, tipus]() {
        return szalak[tipus].waiting_list.getElementNum() == 0 && szalak[tipus].waiting_priority_list.getElementNum() == 0 && szalak[tipus].in_progress_list.getElementNum() == 0;
    });
}


//***********************************************************************
void Para_engine::torli_a_keszeket(szal_tipus tipus) {
//***********************************************************************
    list_elem * kesz_listaelem;
    while((kesz_listaelem = szalak[tipus].finished_list.pop_front()))
        delete kesz_listaelem;
}


//***********************************************************************
Para_engine::~Para_engine() {
//***********************************************************************
    szal_feladat_adatai adat; 
    for (szal_tipus i = szal_tipus(1); i < szt_N; i = szal_tipus(i + 1)) {
        adat.beallitas_szaltorleshez(i); 
        betesz(adat, false);
        for (uns j = 0; j < szalak[i].szalak.size(); j++)
            szalak[i].szalak[j].join();
    }
    if (is_log_szalak) {
        for (void *it = thread_log_lista.get_it(); it != nullptr; thread_log_lista.inc_it(it)) {
            thread_log_elem & akt = thread_log_lista.get_akt(it);
            if (akt.mi_tortenik == thread_log_elem::tlmt_most_indul) {
                for (void *jt = it; jt != nullptr; thread_log_lista.inc_it(jt)) {
                    thread_log_elem & cs = thread_log_lista.get_akt(jt);
                    if (cs.mi_tortenik == thread_log_elem::tlmt_most_van_vege && cs.thread_id == akt.thread_id) {
                        std::chrono::duration<double> diff = cs.ido - akt.ido;
                        akt.id_time = 0;
                        cs.id_time = diff.count();
                        break;
                    }
                }
            }
        }
        FILE * fp;
        const char * szaltipusnev[szt_N] = { "szt_nem_csinal_semmit", "szt_cellafeldolgozo_klaszter", "szt_lepo_klaszter", "szt_fw_klaszter", "szt_fw_egyedi",
            "szt_bw_egyedi", "szt_bw_klaszter", "szt_matrix",  "szt_belso_mentes", "szt_fajlba_mentes" };
        if (fopen_s(&fp, "thread_log.txt", "wt") != 0) {
            printf("\n\nError: Cannot open thread_log.txt\n");
            return;
        }
        thread_log_elem prev;
        bool is_nem_elso = false;
        fprintf(fp, "time\tid\tdb\tevent\tdt\tid_time\ttipus\tmx_fv_tipus\n");
        for (void *it = thread_log_lista.get_it(); it != nullptr; thread_log_lista.inc_it(it)) {
            const thread_log_elem & akt = thread_log_lista.get_akt(it);
            std::chrono::duration<double> diff = akt.ido - indulasi_ido;
            if (is_nem_elso) {
                fprintf(fp, "%.0f\t", diff.count() * 1000000.0);
                fprintf(fp, "0\t");
                fprintf(fp, "%u\t", prev.in_progress_threads_db);
                fprintf(fp, "%u\n", prev.free_threads_db);
            }
            fprintf(fp, "%.0f\t", diff.count() * 1000000.0);
            fprintf(fp, "%u\t", akt.thread_id);
            fprintf(fp, "%u\t", akt.in_progress_threads_db);
            fprintf(fp, "%u\t", akt.free_threads_db);
            fprintf(fp, "%s\t", akt.mi_tortenik==thread_log_elem::tlmt_most_indul ? "start" : "stop");
            if (is_nem_elso) {
                diff = akt.ido - prev.ido;
                fprintf(fp, "%.0f\t", diff.count() * 1000000.0);
            }
            else {
                fprintf(fp, "0\t");
            }
            fprintf(fp, "%.0f\t", akt.id_time * 1000000.0);
            fprintf(fp, "%s\t", szaltipusnev[akt.tipus]);
            fprintf(fp, "%u\n", akt.matrix_fuggveny_tipus);
            prev = akt;
            is_nem_elso = true;
        }
        fclose(fp);
    }
}


}

