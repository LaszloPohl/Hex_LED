//***********************************************************************
// cella osztály és face osztályok cpp
// Creation date:  2018. 08. 10.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#include "cella_es_face.h"
#include "vezerlo.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************



//***********************************************************************
void faces_cella::pre_init_egyedi() {
//***********************************************************************
    const adat_cella & adatcella = *p_adatcella;
    kapcsolodo_cella = adatcella.kapcsolodo_cella_index == 0 ? nullptr : static_cast<faces_cella*>(cellak[adatcella.kapcsolodo_cella_index]);
    faces_dc.set_size(adatcella.facek.size());
    faces_ac.set_size(adatcella.facek.size());
    if (is_el && kapcsolodo_cella != nullptr) { // szétválasztott elektromos-termikus celláknál a termikus cellába gyûjtjük
        p_sum_disszipacio_dc = &kapcsolodo_cella->sum_disszipacio_dc;
        p_sum_ddissz_per_dT_dc = &kapcsolodo_cella->sum_ddissz_per_dT_dc;
    }
    //return (is_el ? 1 : 0) + (is_th ? 1 : 0); center csomópont db
}


//***********************************************************************
void faces_cella::init_dc() {
//***********************************************************************
    const adat_cella & adatcella = *p_adatcella;

    resize_matrixok_dc();

    // face párja beállítása

    if (kapcsolodo_cella != nullptr) { // csak akkor lehet, ha elektrotermikus szimulációnál ez tisztán elektromos vagy termikus cella
        if (is_el) el_center_face_dc.p_parja = &kapcsolodo_cella->th_center_face_dc;
        if (is_th) th_center_face_dc.p_parja = &kapcsolodo_cella->el_center_face_dc;
        for (uns i = 1; i < faces_dc.size(); i++)
            faces_dc.unsafe(i).p_parja = adatcella.facek.unsafe(i).parja == 0 ? nullptr : &kapcsolodo_cella->faces_dc[adatcella.facek.unsafe(i).parja];
    }
    else {
        if (is_el && is_th) {
            el_center_face_dc.p_parja = &th_center_face_dc;
            th_center_face_dc.p_parja = &el_center_face_dc;
            for (uns i = 1; i < faces_dc.size(); i++)
                faces_dc.unsafe(i).p_parja = adatcella.facek.unsafe(i).parja == 0 ? nullptr : &faces_dc[adatcella.facek.unsafe(i).parja];
        }
        else {
            el_center_face_dc.p_parja = nullptr;
            th_center_face_dc.p_parja = nullptr;
            for (uns i = 1; i < faces_dc.size(); i++)
                faces_dc.unsafe(i).p_parja = nullptr;
        }
    }
}


//***********************************************************************
void faces_cella::resize_matrixok_dc() {
//***********************************************************************
    const adat_cella & adatcella = *p_adatcella;
    is_2_fa = akt_sim.p_akt_sim->db_fa == 2;
    uns A_el = 0, A_th = 0, A_elth_1_fa = 0; // hány külsõ el/th face van? (centroid is az)
    uns P_el = 0, P_th = 0; // hány perem el/th face van?
    for (uns i = 1; i < adatcella.facek.size(); i++) {
        const adat_face & aface = adatcella.facek.unsafe(i);
        if (aface.is_perem) {
            if (aface.is_el)P_el++;
            else            P_th++;
        }
        else if (aface.is_kulso) {
            if (aface.tipus == ft_centroid) {
                if (aface.is_el) { el_centroid_index_bemeneti = i; el_centroid_index_matrix = is_2_fa ? A_el : A_elth_1_fa; }
                else             { th_centroid_index_bemeneti = i; th_centroid_index_matrix = is_2_fa ? A_th : A_elth_1_fa; }
            }
            if (aface.is_el)A_el++;
            else            A_th++;
            A_elth_1_fa++;
        }
    }

    is_perem_face = (P_el + P_th != 0);
    uns uj_A_1 = is_2_fa ? A_el : A_elth_1_fa; // (A_el + A_th)
    uns uj_A_2 = is_2_fa ? A_th : 0;
    uns B_el = is_el ? (el_centroid_index_bemeneti > 0 ? 0 : 1) : 0; // ha van elektromos része a cellának, és nem centroid az elektromos csomópont, akkor a B-ben 1 sor/oszlop tartozik hozzá, egyébként 0.
    uns B_th = is_th ? (th_centroid_index_bemeneti > 0 ? 0 : 1) : 0;
    uns uj_B_1 = is_2_fa ? B_el : (B_el + B_th);
    uns uj_B_2 = is_2_fa ? B_th : 0;

    bool is_kell_set_size = false;
    if (A_1 != uj_A_1 || A_2 != uj_A_2 || B_1 != uj_B_1 || B_2 != uj_B_2) {
        is_kell_set_size = true;
        A_1 = uj_A_1;
        A_2 = uj_A_2;
        B_1 = uj_B_1;
        B_2 = uj_B_2;
    }

    if(is_kell_set_size){

        // mátrixok foglalása

        dc_yred_1.set_size(A_1, A_1);
        dc_yred_2.set_size(A_2, A_2);
        dc_jred_1.set_size(A_1);
        dc_jred_2.set_size(A_2);
        dc_UA_1.set_size(A_1);
        dc_UA_2.set_size(A_2);
        dc_IA_1.set_size(A_1);
        dc_IA_2.set_size(A_2);

        dc_ya_1.set_size(A_1, A_1);
        dc_ya_2.set_size(A_2, A_2);
        dc_xa_1.set_size(B_1, A_1);
        dc_xa_2.set_size(B_2, A_2);
        dc_xb_1.set_size(A_1, B_1);
        dc_xb_2.set_size(A_2, B_2);
        dc_yb_1.set_size(B_1, B_1);
        dc_yb_2.set_size(B_2, B_2);
        dc_ja_1.set_size(A_1);
        dc_ja_2.set_size(A_2);
        dc_jb_1.set_size(B_1);
        dc_jb_2.set_size(B_2);

        dc_nzbxa_1.set_size(B_1, A_1);
        dc_nzbxa_2.set_size(B_2, A_2);
        dc_nzbjb_1.set_size(B_1);
        dc_nzbjb_2.set_size(B_2);
        dc_ub_1.set_size(B_1);
        dc_ub_2.set_size(B_2);

        // pointerek beállítása

        feltoltes_tipus = ft_semmi;
        if (p_adatcella->mezotipus == mt_elektromos) {
            feltoltes_tipus = ft_el;
            if (el_centroid_index_bemeneti == 0) {
                uns csatlakozo_index = 0;
                for (uns i = 1; i < faces_dc.size(); i++) {
                    face_tok_dc & akt_tok = faces_dc.unsafe(i);
                    switch (p_adatcella->facek.unsafe(i).tipus) {
                        case ft_csatlakozo: 
                            akt_tok.p_Yii = &dc_ya_1.get_elem_unsafe(csatlakozo_index, csatlakozo_index);
                            akt_tok.p_Yie = &dc_xb_1.get_elem_unsafe(csatlakozo_index); // egy oszlopa van
                            akt_tok.p_Yei = &dc_xa_1.get_elem_unsafe(csatlakozo_index); // egy sora van
                            akt_tok.p_Ji  = &dc_ja_1.unsafe(csatlakozo_index);
                            akt_tok.p_UTi = &dc_UA_1.unsafe(csatlakozo_index);
                            akt_tok.p_IPi = &dc_IA_1.unsafe(csatlakozo_index);
                            csatlakozo_index++;
                            break;
                       // más face típus esetén nincs a face tokban beállítandó
                    }
                }
                p_Yee_dc = &dc_yb_1.get_elem_unsafe(0);
                p_Je_dc  = &dc_jb_1.unsafe(0);
                el_center_face_dc.p_UTi = &dc_ub_1.unsafe(0);
            }
            else { // van el centroid
                uns csatlakozo_index = 0;
                for (uns i = 1; i < faces_dc.size(); i++) {
                    face_tok_dc & akt_tok = faces_dc.unsafe(i);
                    switch (p_adatcella->facek.unsafe(i).tipus) {
                        case ft_csatlakozo:
                            akt_tok.p_Yii = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         csatlakozo_index);
                            akt_tok.p_Yie = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         el_centroid_index_matrix);
                            akt_tok.p_Yei = &dc_ya_1.get_elem_unsafe(el_centroid_index_matrix, csatlakozo_index);
                            akt_tok.p_Ji  = &dc_ja_1.unsafe(csatlakozo_index);
                            akt_tok.p_UTi = &dc_UA_1.unsafe(csatlakozo_index);
                            akt_tok.p_IPi = &dc_IA_1.unsafe(csatlakozo_index);
                            csatlakozo_index++;
                            break;
                        case ft_centroid:
                            egyforma_e_hiba(csatlakozo_index, el_centroid_index_matrix, "cella::resize_matrixok_dc() el csatlakozo_index!=centroid_index");
                            // A centroid face nem állít Yt ill J-t, mert csak egy alteregója a középsõ csomópontnak, ezrt nem kell beállítani a peremétereit
                            break;
                            // más face típus esetén nincs a face tokban beállítandó
                    }
                }
                p_Yee_dc = &dc_ya_1.get_elem_unsafe(el_centroid_index_matrix, el_centroid_index_matrix);
                p_Je_dc  = &dc_ja_1.unsafe(el_centroid_index_matrix);
                el_center_face_dc.p_UTi = &dc_UA_1.unsafe(el_centroid_index_matrix);
            }
        }
        else if (p_adatcella->mezotipus == mt_termikus) {
            feltoltes_tipus = ft_th;
            if (th_centroid_index_bemeneti == 0) {
                uns csatlakozo_index = 0;
                for (uns i = 1; i < faces_dc.size(); i++) {
                    face_tok_dc & akt_tok = faces_dc.unsafe(i);
                    switch (p_adatcella->facek.unsafe(i).tipus) {
                        case ft_csatlakozo:
                            akt_tok.p_Yii = &dc_ya_1.get_elem_unsafe(csatlakozo_index, csatlakozo_index);
                            akt_tok.p_Yit = &dc_xb_1.get_elem_unsafe(csatlakozo_index); // egy oszlopa van
                            akt_tok.p_Yti = &dc_xa_1.get_elem_unsafe(csatlakozo_index); // egy sora van
                            akt_tok.p_Ji  = &dc_ja_1.unsafe(csatlakozo_index);
                            akt_tok.p_UTi = &dc_UA_1.unsafe(csatlakozo_index);
                            akt_tok.p_IPi = &dc_IA_1.unsafe(csatlakozo_index);
                            csatlakozo_index++;
                            break;
                            // más face típus esetén nincs a face tokban beállítandó
                    }
                }
                p_Ytt_dc = &dc_yb_1.get_elem_unsafe(0);
                p_Jt_dc  = &dc_jb_1.unsafe(0);
                th_center_face_dc.p_UTi = &dc_ub_1.unsafe(0);
            }
            else { // van th centroid
                uns csatlakozo_index = 0;
                for (uns i = 1; i < faces_dc.size(); i++) {
                    face_tok_dc & akt_tok = faces_dc.unsafe(i);
                    switch (p_adatcella->facek.unsafe(i).tipus) {
                        case ft_csatlakozo:
                            akt_tok.p_Yii = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         csatlakozo_index);
                            akt_tok.p_Yit = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         th_centroid_index_matrix);
                            akt_tok.p_Yti = &dc_ya_1.get_elem_unsafe(th_centroid_index_matrix, csatlakozo_index);
                            akt_tok.p_Ji  = &dc_ja_1.unsafe(csatlakozo_index);
                            akt_tok.p_UTi = &dc_UA_1.unsafe(csatlakozo_index);
                            akt_tok.p_IPi = &dc_IA_1.unsafe(csatlakozo_index);
                            csatlakozo_index++;
                            break;
                        case ft_centroid:
                            egyforma_e_hiba(csatlakozo_index, th_centroid_index_matrix, "cella::resize_matrixok_dc() th csatlakozo_index!=centroid_index");
                            // A centroid face nem állít Yt ill J-t, mert csak egy alteregója a középsõ csomópontnak, ezrt nem kell beállítani a peremétereit
                            break;
                            // más face típus esetén nincs a face tokban beállítandó
                    }
                }
                p_Ytt_dc = &dc_ya_1.get_elem_unsafe(th_centroid_index_matrix, th_centroid_index_matrix);
                p_Jt_dc  = &dc_ja_1.unsafe(th_centroid_index_matrix);
                th_center_face_dc.p_UTi = &dc_UA_1.unsafe(th_centroid_index_matrix);
            }
        }
        else if (p_adatcella->mezotipus == mt_elektrotermikus) {
            if (is_2_fa) { // két fa esetén nincsenek vegyes Y-ok vagy J-k, azaz ugyanazokat a pointereket kell feltölteni, mint egyteres szimulációnál

                feltoltes_tipus = ft_elth_nincs_csatolt;

                // 1. fa = elektromos

                if (el_centroid_index_bemeneti == 0) {
                    uns csatlakozo_index = 0;
                    for (uns i = 1; i < p_adatcella->facek.size(); i++)
                        if (p_adatcella->facek.unsafe(i).is_el) {
                            face_tok_dc & akt_tok = faces_dc.unsafe(i);
                            switch (p_adatcella->facek.unsafe(i).tipus) {
                                case ft_csatlakozo:
                                    akt_tok.p_Yii = &dc_ya_1.get_elem_unsafe(csatlakozo_index, csatlakozo_index);
                                    akt_tok.p_Yie = &dc_xb_1.get_elem_unsafe(csatlakozo_index); // egy oszlopa van
                                    akt_tok.p_Yei = &dc_xa_1.get_elem_unsafe(csatlakozo_index); // egy sora van
                                    akt_tok.p_Ji  = &dc_ja_1.unsafe(csatlakozo_index);
                                    akt_tok.p_UTi = &dc_UA_1.unsafe(csatlakozo_index);
                                    akt_tok.p_IPi = &dc_IA_1.unsafe(csatlakozo_index);
                                    csatlakozo_index++;
                                    break;
                                    // más face típus esetén nincs a face tokban beállítandó
                            }
                        }
                    p_Yee_dc = &dc_yb_1.get_elem_unsafe(0);
                    p_Je_dc  = &dc_jb_1.unsafe(0);
                    el_center_face_dc.p_UTi = &dc_ub_1.unsafe(0);
                }
                else { // van el centroid
                    uns csatlakozo_index = 0;
                    for (uns i = 1; i < p_adatcella->facek.size(); i++)
                        if (p_adatcella->facek.unsafe(i).is_el) {
                            face_tok_dc & akt_tok = faces_dc.unsafe(i);
                            switch (p_adatcella->facek.unsafe(i).tipus) {
                                case ft_csatlakozo:
                                    akt_tok.p_Yii = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         csatlakozo_index);
                                    akt_tok.p_Yie = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         el_centroid_index_matrix);
                                    akt_tok.p_Yei = &dc_ya_1.get_elem_unsafe(el_centroid_index_matrix, csatlakozo_index);
                                    akt_tok.p_Ji  = &dc_ja_1.unsafe(csatlakozo_index);
                                    akt_tok.p_UTi = &dc_UA_1.unsafe(csatlakozo_index);
                                    akt_tok.p_IPi = &dc_IA_1.unsafe(csatlakozo_index);
                                    csatlakozo_index++;
                                    break;
                                case ft_centroid:
                                    egyforma_e_hiba(csatlakozo_index, el_centroid_index_matrix, "cella::resize_matrixok_dc() el csatlakozo_index!=centroid_index");
                                    // A centroid face nem állít Yt ill J-t, mert csak egy alteregója a középsõ csomópontnak, ezrt nem kell beállítani a peremétereit
                                    break;
                                    // más face típus esetén nincs a face tokban beállítandó
                            }
                        }
                    p_Yee_dc = &dc_ya_1.get_elem_unsafe(el_centroid_index_matrix, el_centroid_index_matrix);
                    p_Je_dc  = &dc_ja_1.unsafe(el_centroid_index_matrix);
                    el_center_face_dc.p_UTi = &dc_UA_1.unsafe(el_centroid_index_matrix);
                }

                // 2. fa = termikus

                if (th_centroid_index_bemeneti == 0) {
                    uns csatlakozo_index = 0;
                    for (uns i = 1; i < p_adatcella->facek.size(); i++)
                        if (!p_adatcella->facek.unsafe(i).is_el) {
                            face_tok_dc & akt_tok = faces_dc.unsafe(i);
                            switch (p_adatcella->facek.unsafe(i).tipus) {
                                case ft_csatlakozo:
                                    akt_tok.p_Yii = &dc_ya_2.get_elem_unsafe(csatlakozo_index, csatlakozo_index);
                                    akt_tok.p_Yit = &dc_xb_2.get_elem_unsafe(csatlakozo_index); // egy oszlopa van
                                    akt_tok.p_Yti = &dc_xa_2.get_elem_unsafe(csatlakozo_index); // egy sora van
                                    akt_tok.p_Ji  = &dc_ja_2.unsafe(csatlakozo_index);
                                    akt_tok.p_UTi = &dc_UA_2.unsafe(csatlakozo_index);
                                    akt_tok.p_IPi = &dc_IA_2.unsafe(csatlakozo_index);
                                    csatlakozo_index++;
                                    break;
                                    // más face típus esetén nincs a face tokban beállítandó
                            }
                        }
                    p_Ytt_dc = &dc_yb_2.get_elem_unsafe(0);
                    p_Jt_dc  = &dc_jb_2.unsafe(0);
                    th_center_face_dc.p_UTi = &dc_ub_2.unsafe(0);
                }
                else { // van th centroid
                    uns csatlakozo_index = 0;
                    for (uns i = 1; i < p_adatcella->facek.size(); i++)
                        if (!p_adatcella->facek.unsafe(i).is_el) {
                            face_tok_dc & akt_tok = faces_dc.unsafe(i);
                            switch (p_adatcella->facek.unsafe(i).tipus) {
                                case ft_csatlakozo:
                                    akt_tok.p_Yii = &dc_ya_2.get_elem_unsafe(csatlakozo_index,         csatlakozo_index);
                                    akt_tok.p_Yit = &dc_ya_2.get_elem_unsafe(csatlakozo_index,         th_centroid_index_matrix);
                                    akt_tok.p_Yti = &dc_ya_2.get_elem_unsafe(th_centroid_index_matrix, csatlakozo_index);
                                    akt_tok.p_Ji  = &dc_ja_2.unsafe(csatlakozo_index);
                                    akt_tok.p_UTi = &dc_UA_2.unsafe(csatlakozo_index);
                                    akt_tok.p_IPi = &dc_IA_2.unsafe(csatlakozo_index);
                                    csatlakozo_index++;
                                    break;
                                case ft_centroid:
                                    egyforma_e_hiba(csatlakozo_index, th_centroid_index_matrix, "cella::resize_matrixok_dc() th csatlakozo_index!=centroid_index");
                                    // A centroid face nem állít Yt ill J-t, mert csak egy alteregója a középsõ csomópontnak, ezrt nem kell beállítani a peremétereit
                                    break;
                                    // más face típus esetén nincs a face tokban beállítandó
                            }
                        }
                    p_Ytt_dc = &dc_ya_2.get_elem_unsafe(th_centroid_index_matrix, th_centroid_index_matrix);
                    p_Jt_dc  = &dc_ja_2.unsafe(th_centroid_index_matrix);
                    th_center_face_dc.p_UTi = &dc_UA_2.unsafe(th_centroid_index_matrix);
                }
            }
            else { // 1 fa van
                feltoltes_tipus = ft_elth_van_csatolt;
                if (el_centroid_index_bemeneti == 0 && th_centroid_index_bemeneti == 0) { // nincs centroid
                    uns csatlakozo_index = 0;
                    for (uns i = 1; i < p_adatcella->facek.size(); i++) {
                        face_tok_dc & akt_tok = faces_dc.unsafe(i);
                        switch (p_adatcella->facek.unsafe(i).tipus) {
                            case ft_csatlakozo: 
                                akt_tok.p_Yii = &dc_ya_1.get_elem_unsafe(csatlakozo_index, csatlakozo_index);
                                akt_tok.p_Yie = &dc_xb_1.get_elem_unsafe(csatlakozo_index, 0); // két oszlopa van
                                akt_tok.p_Yit = &dc_xb_1.get_elem_unsafe(csatlakozo_index, 1); // két oszlopa van
                                akt_tok.p_Yei = &dc_xa_1.get_elem_unsafe(0, csatlakozo_index); // két sora van
                                akt_tok.p_Yti = &dc_xa_1.get_elem_unsafe(1, csatlakozo_index); // két sora van
                                akt_tok.p_Ji  = &dc_ja_1.unsafe(csatlakozo_index);
                                akt_tok.p_UTi = &dc_UA_1.unsafe(csatlakozo_index);
                                akt_tok.p_IPi = &dc_IA_1.unsafe(csatlakozo_index);
                                csatlakozo_index++;
                                break;
                           // más face típus esetén nincs a face tokban beállítandó
                        }
                    }
                    p_Yee_dc = &dc_yb_1.get_elem_unsafe(0); // 0,0
                    p_Yet_dc = &dc_yb_1.get_elem_unsafe(1); // 0,1
                    p_Yte_dc = &dc_yb_1.get_elem_unsafe(2); // 1,0
                    p_Ytt_dc = &dc_yb_1.get_elem_unsafe(3); // 1,1
                    p_Je_dc  = &dc_jb_1.unsafe(0);
                    p_Jt_dc  = &dc_jb_1.unsafe(1);
                    el_center_face_dc.p_UTi = &dc_ub_1.unsafe(0);
                    th_center_face_dc.p_UTi = &dc_ub_1.unsafe(1);
                }
                else if (th_centroid_index_bemeneti == 0) { // elektromos centroid van
                    uns csatlakozo_index = 0;
                    for (uns i = 1; i < p_adatcella->facek.size(); i++) {
                        face_tok_dc & akt_tok = faces_dc.unsafe(i);
                        switch (p_adatcella->facek.unsafe(i).tipus) {
                            case ft_csatlakozo:
                                akt_tok.p_Yii = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         csatlakozo_index);
                                akt_tok.p_Yie = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         el_centroid_index_matrix);
                                akt_tok.p_Yit = &dc_xb_1.get_elem_unsafe(csatlakozo_index); // egy oszopa van
                                akt_tok.p_Yei = &dc_ya_1.get_elem_unsafe(el_centroid_index_matrix, csatlakozo_index);
                                akt_tok.p_Yti = &dc_xa_1.get_elem_unsafe(csatlakozo_index); // egy sora van
                                akt_tok.p_Ji  = &dc_ja_1.unsafe(csatlakozo_index);
                                akt_tok.p_UTi = &dc_UA_1.unsafe(csatlakozo_index);
                                akt_tok.p_IPi = &dc_IA_1.unsafe(csatlakozo_index);
                                csatlakozo_index++;
                                break;
                            case ft_centroid:
                                egyforma_e_hiba(csatlakozo_index, el_centroid_index_matrix, "cella::resize_matrixok_dc() el csatlakozo_index!=centroid_index");
                                // A centroid face nem állít Yt ill J-t, mert csak egy alteregója a középsõ csomópontnak, ezrt nem kell beállítani a peremétereit
                                break;
                                // más face típus esetén nincs a face tokban beállítandó
                        }
                    }
                    p_Yee_dc = &dc_ya_1.get_elem_unsafe(el_centroid_index_matrix, el_centroid_index_matrix);
                    p_Yet_dc = &dc_xb_1.get_elem_unsafe(el_centroid_index_matrix); // egy oszopa van
                    p_Yte_dc = &dc_xa_1.get_elem_unsafe(el_centroid_index_matrix); // egy sora van
                    p_Ytt_dc = &dc_yb_1.get_elem_unsafe(0);
                    p_Je_dc  = &dc_ja_1.unsafe(el_centroid_index_matrix);
                    p_Jt_dc  = &dc_jb_1.unsafe(0);
                    el_center_face_dc.p_UTi = &dc_UA_1.unsafe(el_centroid_index_matrix);
                    th_center_face_dc.p_UTi = &dc_ub_1.unsafe(0);
                }
                else if (el_centroid_index_bemeneti == 0) { // termikus centroid van
                    uns csatlakozo_index = 0;
                    for (uns i = 1; i < p_adatcella->facek.size(); i++) {
                        face_tok_dc & akt_tok = faces_dc.unsafe(i);
                        switch (p_adatcella->facek.unsafe(i).tipus) {
                            case ft_csatlakozo:
                                akt_tok.p_Yii = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         csatlakozo_index);
                                akt_tok.p_Yie = &dc_xb_1.get_elem_unsafe(csatlakozo_index); // egy oszopa van
                                akt_tok.p_Yit = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         th_centroid_index_matrix);
                                akt_tok.p_Yei = &dc_xa_1.get_elem_unsafe(csatlakozo_index); // egy sora van
                                akt_tok.p_Yti = &dc_ya_1.get_elem_unsafe(th_centroid_index_matrix, csatlakozo_index);
                                akt_tok.p_Ji  = &dc_ja_1.unsafe(csatlakozo_index);
                                akt_tok.p_UTi = &dc_UA_1.unsafe(csatlakozo_index);
                                akt_tok.p_IPi = &dc_IA_1.unsafe(csatlakozo_index);
                                csatlakozo_index++;
                                break;
                            case ft_centroid:
                                egyforma_e_hiba(csatlakozo_index, th_centroid_index_matrix, "cella::resize_matrixok_dc() th csatlakozo_index!=centroid_index");
                                // A centroid face nem állít Yt ill J-t, mert csak egy alteregója a középsõ csomópontnak, ezrt nem kell beállítani a peremétereit
                                break;
                                // más face típus esetén nincs a face tokban beállítandó
                        }
                    }
                    p_Yee_dc = &dc_yb_1.get_elem_unsafe(0);
                    p_Yet_dc = &dc_xa_1.get_elem_unsafe(th_centroid_index_matrix); // egy sora van
                    p_Yte_dc = &dc_xb_1.get_elem_unsafe(th_centroid_index_matrix); // egy oszopa van
                    p_Ytt_dc = &dc_ya_1.get_elem_unsafe(th_centroid_index_matrix, th_centroid_index_matrix);
                    p_Je_dc  = &dc_jb_1.unsafe(0);
                    p_Jt_dc  = &dc_ja_1.unsafe(th_centroid_index_matrix);
                    el_center_face_dc.p_UTi = &dc_ub_1.unsafe(0);
                    th_center_face_dc.p_UTi = &dc_UA_1.unsafe(th_centroid_index_matrix);
                }
                else { // elektromos és termikus centroid is van
                    uns csatlakozo_index = 0;
                    for (uns i = 1; i < p_adatcella->facek.size(); i++) {
                        face_tok_dc & akt_tok = faces_dc.unsafe(i);
                        switch (p_adatcella->facek.unsafe(i).tipus) {
                            case ft_csatlakozo:
                                akt_tok.p_Yii = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         csatlakozo_index);
                                akt_tok.p_Yie = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         el_centroid_index_matrix);
                                akt_tok.p_Yit = &dc_ya_1.get_elem_unsafe(csatlakozo_index,         th_centroid_index_matrix);
                                akt_tok.p_Yei = &dc_ya_1.get_elem_unsafe(el_centroid_index_matrix, csatlakozo_index);
                                akt_tok.p_Yti = &dc_ya_1.get_elem_unsafe(th_centroid_index_matrix, csatlakozo_index);
                                akt_tok.p_Ji  = &dc_ja_1.unsafe(csatlakozo_index);
                                akt_tok.p_UTi = &dc_UA_1.unsafe(csatlakozo_index);
                                akt_tok.p_IPi = &dc_IA_1.unsafe(csatlakozo_index);
                                csatlakozo_index++;
                                break;
                            case ft_centroid:
                                if (p_adatcella->facek.unsafe(i).is_el)
                                    egyforma_e_hiba(csatlakozo_index, el_centroid_index_matrix, "cella::resize_matrixok_dc() el csatlakozo_index!=centroid_index");
                                else
                                    egyforma_e_hiba(csatlakozo_index, th_centroid_index_matrix, "cella::resize_matrixok_dc() th csatlakozo_index!=centroid_index");
                                // A centroid face nem állít Yt ill J-t, mert csak egy alteregója a középsõ csomópontnak, ezrt nem kell beállítani a peremétereit
                                break;
                                // más face típus esetén nincs a face tokban beállítandó
                        }
                    }
                    p_Yee_dc = &dc_ya_1.get_elem_unsafe(el_centroid_index_matrix, el_centroid_index_matrix);
                    p_Yet_dc = &dc_ya_1.get_elem_unsafe(el_centroid_index_matrix, th_centroid_index_matrix);
                    p_Yte_dc = &dc_ya_1.get_elem_unsafe(th_centroid_index_matrix, el_centroid_index_matrix);
                    p_Ytt_dc = &dc_ya_1.get_elem_unsafe(th_centroid_index_matrix, th_centroid_index_matrix);
                    p_Je_dc  = &dc_ja_1.unsafe(el_centroid_index_matrix);
                    p_Jt_dc  = &dc_ja_1.unsafe(th_centroid_index_matrix);
                    el_center_face_dc.p_UTi = &dc_UA_1.unsafe(el_centroid_index_matrix);
                    th_center_face_dc.p_UTi = &dc_UA_1.unsafe(th_centroid_index_matrix);
                }
            }
        }
        else throw hiba("cella::resize_matrixok_dc()", "unknown cell type (th/el/elth/?)");
    }
}


//***********************************************************************
void os_cella::reset_H() {
//***********************************************************************
    emlekek.get_akt().H = ertek_t();
    emlekek.get_elozo_to_overwrite().H = ertek_t();
    emlekek.get_prevprev_to_overwrite().H = ertek_t();
    emlekek.get_megtartando_to_overwrite().H = ertek_t();
    if (p_cella_anyag!=nullptr && p_cella_anyag->fazisvaltas.is_fazisvalto()) {
        emlekek.get_akt().H = p_cella_anyag->fazisvaltas.H(rvt(), rvt()); // H = 0 az elõzõ állapot, és ambient hõmérsékletre áll be
    }
}


//***********************************************************************
void faces_cella::halozat_foglalas_cellaban_dc() {
//***********************************************************************
    const adat_cella & adatcella = *p_adatcella;
    set_hiba_hol h("faces_cella::halozat_foglalas_cellaban_dc");

    is_lin = true;
    p_cella_anyag = &akt_sim.p_akt_sim->anyagok[akt_sim.akt_anyag_index[adatcella.anyag_index]];

    if (p_cella_anyag->fazisvaltas.is_fazisvalto()) { 
        is_lin = false;
    }

    // csatlakozó és perem face-ek

    A_cella = 0; // Nem használjuk sehol, ezért nem viszem át az õsbe
    A_top = 1; // ha nincs top, akkor 1-gyel osszunk
    for (uns i = 1; i < faces_dc.size(); i++) {

        const adat_face & akt_face = adatcella.facek[i];

        // milyen típus kell ?

        face_azonosito azon = fa_none;
        switch (akt_face.tipus) {
            case ft_csatlakozo:
            case ft_normal_perem:
            case ft_spec_perem: {

                if (akt_face.tipus == ft_spec_perem)
                    TODO("cella::halozat_foglalas_cellaban_dc ft_spec_perem");

                // bool is_homersekletfuggo = false;
                // bool is_junction = false;
                
                // TODO: seebeck
                // belsõ rész típusának meghatározása
                // 1. anyag

                azon = fa_belso_R;

                const adat_anyag & akt_anyag = akt_sim.p_akt_sim->anyagok[akt_sim.akt_anyag_index[akt_face.anyag_index]];
                if (akt_face.is_el) {
                    if (!akt_anyag.elvez.is_exists())
                        throw hiba(1, "material (%u) electrical conduction is undefined but required", akt_face.anyag_index);
                    if (!akt_anyag.elvez.is_const())
                        azon = face_azonosito(azon | fa_hom_fuggo);
                    if (is_elektrotermikus_szimulacio)
                        azon = face_azonosito(azon | fa_disszipalo | fa_el_specific);

                }
                else {
                    if (!akt_anyag.thvez.is_exists())
                        throw hiba(1, "material (%u) thermal conduction is undefined but required", akt_face.anyag_index);
                    if (!akt_anyag.thvez.is_const())
                        azon = face_azonosito(azon | fa_hom_fuggo);
                }

                // 2. junction

                if(akt_face.is_el && akt_sim.akt_junct_index[akt_face.junction_index] != 0) // nem csinálunk termikus face-t a junctionnak. miért kéne? de ha mégis kell, szabad.
                    azon = face_azonosito(azon | fa_belso_J | fa_hom_fuggo); // a junction-t mindig hõmérsékletfüggõnek tekintjük (fölösleges külön tárgyalni a nem hõmérsékletfüggõ esetet)

                // perem típusa

                if (akt_face.tipus == ft_normal_perem) {
                    uns akt_perem_index = akt_sim.akt_perem_index[akt_face.perem_index];
                    if (akt_perem_index == 0) {
                        azon = fa_od_r_kdn; // open peremnél mindegy, hogy hõmérsékletfüggõ-e ill. disszipáló-e, nem is térspecifikus, vagy van-e juction
                    }
                    else {
                        switch (akt_sim.p_akt_sim->perem[akt_perem_index].tipus) {
                            case pt_u:      azon = face_azonosito(azon | fa_peremre | fa_perem_UT); break;
                            case pt_t:      azon = face_azonosito(azon | fa_peremre | fa_perem_UT); break;
                            case pt_htc:    azon = face_azonosito(azon | fa_peremre | fa_perem_HTC); break;
                            case pt_thtc:   azon = face_azonosito(azon | fa_peremre | fa_perem_HTC_T); break;
                            default: throw hiba(1, "unknown boundary type (%u)", uns(akt_sim.p_akt_sim->perem[akt_perem_index].tipus));
                        }
                    }
                }
                else if (akt_face.tipus == ft_spec_perem) {
                    TODO("Spec perem");
                }
                else { // csatlakozó
                    azon = face_azonosito(azon | fa_csatlakozo);
                }

            }
            break;
            case ft_centroid: {
                azon = fa_Cd__d;
            }
            break;
            default: throw hiba(1, "unknown face type (%u)", uns(akt_face.tipus));
        }
        // ha nem ilyen típusú, akkor újat foglal

        if (faces_dc[i].p_core == nullptr || faces_dc[i].p_core->get_azon() != azon) {
            delete faces_dc[i].p_core;
            faces_dc[i].p_core = os_face_core_dc::alloc_new(azon, &faces_dc[i], this, &akt_face, akt_face.is_el);
        }

        // beállítja

        faces_dc[i].p_core->init();
        if ((azon & fa_hom_fuggo) != 0 || (azon & fa_I_fuggo) != 0 || (azon & fa_belso_J) != 0)
            is_lin = false;
    }
}


//***********************************************************************
void faces_cella::gerj_update_dc() {
//***********************************************************************
// változhat a gerjesztés típusa is, ezért itt foglaljuk a középponti face-eket

    const adat_struktura & akt_stukt = akt_sim.akt_strukturak[p_adatcella->struktura_index];
    if (is_el) {

        // milyen típus kell ?
        
        face_azonosito azon;
        gerjesztes_tipus tip = akt_stukt.p_el_gerj == nullptr ? gt_el_none : akt_stukt.p_el_gerj->tipus;
        switch (tip) {
            case gt_el_none:    azon = fa_kd_x_kdn; break;
            case gt_I:          azon = fa_kd_i_kdn; break;
            case gt_U:          azon = fa_kd_u_kdn; break;
            default: throw hiba("cella::gerj_update_dc", "bad electrical excitation type");
        }

        // ha nem ilyen típusú, akkor újat foglal

        if (el_center_face_dc.p_core == nullptr || el_center_face_dc.p_core->get_azon() != azon) {
            delete el_center_face_dc.p_core;
            el_center_face_dc.p_core = os_face_core_dc::alloc_new(azon, &el_center_face_dc, this, nullptr, true);
        }

        // beállítja

        el_center_face_dc.p_core->init();

        // frissítés

        el_center_face_dc.p_core->update_gerj();
    }
    if (is_th) {
        
        // milyen típus kell ?

        face_azonosito azon;
        const adat_anyag & akt_anyag = akt_sim.p_akt_sim->anyagok[akt_sim.akt_anyag_index[p_adatcella->anyag_index]];
        gerjesztes_tipus tip = akt_stukt.p_th_gerj == nullptr ? gt_th_none : akt_stukt.p_th_gerj->tipus;
        switch (tip) {
            case gt_th_none:    azon = akt_sim.tipus == alt_trans ? (akt_anyag.Cth.is_const() ? fa_kd_ci_ktn : fa_kd_ci_ttn) : fa_kd_x_kdn; break;
            case gt_P:          azon = akt_sim.tipus == alt_trans ? (akt_anyag.Cth.is_const() ? fa_kd_ci_ktn : fa_kd_ci_ttn) : fa_kd_i_kdn; break;
            case gt_T:          azon = fa_kd_u_kdn; break;
            default: throw hiba("cella::gerj_update_dc", "bad thermal excitation type");
        }
        
        // ha nem ilyen típusú, akkor újat foglal

        if (th_center_face_dc.p_core == nullptr || th_center_face_dc.p_core->get_azon() != azon) {
            delete th_center_face_dc.p_core;
            th_center_face_dc.p_core = os_face_core_dc::alloc_new(azon, &th_center_face_dc, this, nullptr, false);
        }

        // beállítja

        th_center_face_dc.p_core->init();
        if ((azon & fa_hom_fuggo) != 0)
            is_lin = false;

        // frissítés

        th_center_face_dc.p_core->update_gerj();
    }

    // TODO: ha megengedett lenne a keresztface-eknél is a gerjesztés, akkor itt végig kellene járni az összeset, és meghívni az update_gerj-üket
}


//***********************************************************************
void faces_cella::set_is_szimm() {
//***********************************************************************
    is_NR_nincs_csatolt_nemszimm = is_NR_van_csatolt_nemszimm = is_nemszimm = false;
    if (akt_sim.tipus == alt_ac)
        return; // AC csak szimmetrikus lehet
    if (akt_sim.p_akt_sim->nemlin_tipus == nmt_klasszik_iteracio)
        return; // A szukcesszív approximáció csak szimmetrikus lehet

    // Newton-Raphson van

    if (konst_Yet_dc != konst_Yte_dc)
        set_to_NR_van_csatolt_nemszimm(); // TODO: update_for_uj_lepes_dc állíthatja be a két konst-ot, de csak is_uj_cellaszerkezet_kell 
                                          // esetén hívjuk set_is_szimm()-et a update_for_uj_lepes_dc után, így ha esetleg a jövõben
                                          // lenne egymástól különbözõ konst_Yet_dc vagy konst_Yte_dc, meg kell nézni, hogy okoz-e
                                          // problémát, hogy nem új cellaszerkezet esetén itt a korábbi értéket vizsgáljuk.
                                          // (Egyébként mivel ez a fv nagyon gyors, be is lehet rakni minden iterációba, de fölöslegesen minek?)
    for (uns i = 1; i < faces_dc.size(); i++)
        faces_dc.unsafe(i).p_core->set_cella_is_szimm();
    if (is_el)el_center_face_dc.p_core->set_cella_is_szimm();
    if (is_th)th_center_face_dc.p_core->set_cella_is_szimm();

    // A cella nemszimmetriájának beállítása (ide csak akkor jutunk, ha DC NR a szimuláció)

    if (feltoltes_tipus == ft_elth_van_csatolt)
        is_nemszimm = is_NR_van_csatolt_nemszimm;
    else
        is_nemszimm = is_NR_nincs_csatolt_nemszimm;
}


//***********************************************************************
void faces_cella::update_dt_dc() {
//***********************************************************************
    for (uns i = 1; i < faces_dc.size(); i++)
        faces_dc.unsafe(i).p_core->update_dt();
    if (is_el)el_center_face_dc.p_core->update_dt();
    if (is_th)th_center_face_dc.p_core->update_dt();
}


//***********************************************************************
void faces_cella::peremfeltetel_update_dc() {
// itt nem változhat a perem típusa, ezért a halozat_foglalas_cellaban_dc-ben marad a
// nem középponti face-ek foglalása
//***********************************************************************
    if(is_perem_face) // ha a cellában nincs perem,akkor nem is fut a ciklus
        for (uns i = 1; i < faces_dc.size(); i++)
            faces_dc.unsafe(i).p_core->update_peremfeltetel();
}


//***********************************************************************
void faces_cella::del_all_prev_dc() {
//***********************************************************************
    for (uns i = 0; i < faces_dc.size(); i++) {
        faces_dc[i].emlekek.clear_megtartando_is();
    }
    el_center_face_dc.emlekek.clear_megtartando_is();
    th_center_face_dc.emlekek.clear_megtartando_is();
    // kiinduló állapotba hozzuk a hiszterézist
    emlekek.get_akt().H = ertek_t();
    emlekek.get_elozo_to_overwrite().H = ertek_t();
    emlekek.get_prevprev_to_overwrite().H = ertek_t();
    emlekek.get_megtartando_to_overwrite().H = ertek_t();
}


//***********************************************************************
void faces_cella::update_for_uj_lepes_dc(){
// A hiszterézis állapotát és a kapacitások feszültségét az emlékezõ tár
// automatikusan kezeli, itt nincs dolgunk velük.
// ha face foglalás vagy gerj update volt, akkor be kell írni az 
// admittancia mátrixban minden értéket, hogy a konstansok aktuálisak 
// legyenek (ekkor elõbb 0-za a mátrixokat és vektorokat)
//***********************************************************************
    if (akt_sim.is_uj_facek_letrehozasa_kell || akt_sim.is_gerj_update_kell || akt_sim.is_force_face_update || akt_sim.is_subiter_dt_changed) {
        is_change_adm = is_change_j = true;
        dc_ya_1.zero_nembiztos();
        dc_ya_2.zero_nembiztos();
        dc_xa_1.zero_nembiztos();
        dc_xa_2.zero_nembiztos();
        dc_xb_1.zero_nembiztos();
        dc_xb_2.zero_nembiztos();
        dc_yb_1.zero_nembiztos();
        dc_yb_2.zero_nembiztos();
        dc_ja_1.zero();
        dc_ja_2.zero();
        dc_jb_1.zero();
        dc_jb_2.zero();
        konst_Yee_dc = konst_Yet_dc = konst_Yte_dc = konst_Ytt_dc = konst_Je_dc = konst_Jt_dc = rvt();
        for (uns i = 1; i < faces_dc.size(); i++)
            faces_dc.unsafe(i).p_core->update_uj_lepeshez();
        if(is_el)
            el_center_face_dc.p_core->update_uj_lepeshez();
        if(is_th)
            th_center_face_dc.p_core->update_uj_lepeshez();
    }
}


//***********************************************************************
void faces_cella::fw_dc(uns melyik_fa) {
// TODO: Ha DC-AC-DC váltás során AC-ben változik valami, akkor nem 
// biztos, hogy ez a felépítés jól fog mûködni.
// is_change_adm és is_change_j beállítása, ha változott
//***********************************************************************
    fw_klaszter_1_update_Jakobi_dc();
    // TODO: itt kellene az is_lin-t figyelembe venni, de majd egyszer, most mindig redukljunk
    if (melyik_fa == 1) {
        if (B_1 == 1) { // dc_xat_1 és dc_nzbxat_1 érvénytelen, mert nem számoljuk ki!
            dc_nzbxa_1.math_1_ninv_mul(dc_yb_1, dc_xa_1);
            dc_yred_1.math_1_add_mul(dc_ya_1, dc_xb_1, dc_nzbxa_1);
            math_1x1_mul(dc_nzbjb_1, dc_yb_1, dc_jb_1);
            math_1_add_mul_jred(dc_jred_1, dc_ja_1, dc_xb_1, dc_nzbjb_1);
        }
        else if (B_1 == 2) {
            if(is_nemszimm) dc_nzbxa_1.math_2_ninv_mul(dc_yb_1, dc_xa_1);
            else            dc_nzbxa_1.math_2_ninv_mul_symm(dc_yb_1, dc_xa_1);
            dc_yred_1.math_2_add_mul(dc_ya_1, dc_xb_1, dc_nzbxa_1);
            math_2x2_mul(dc_nzbjb_1, dc_yb_1, dc_jb_1);
            math_2_add_mul_jred(dc_jred_1, dc_ja_1, dc_xb_1, dc_nzbjb_1);
        }
        else { // centroidnál tán lehet 0
            throw hiba("cella::fw_dc", "B size 1 or 2 expected, %u arrived", B_1);
        }
    }
    else {
        if (B_2 == 1) { // dc_xat_1 és dc_nzbxat_1 érvénytelen, mert nem számoljuk ki!
            dc_nzbxa_2.math_1_ninv_mul(dc_yb_2, dc_xa_2);
            dc_yred_2.math_1_add_mul(dc_ya_2, dc_xb_2, dc_nzbxa_2);
            math_1x1_mul(dc_nzbjb_2, dc_yb_2, dc_jb_2);
            math_1_add_mul_jred(dc_jred_2, dc_ja_2, dc_xb_2, dc_nzbjb_2);
        }
        else if (B_2 == 2) {
            if (is_nemszimm) dc_nzbxa_2.math_2_ninv_mul(dc_yb_2, dc_xa_2);
            else             dc_nzbxa_2.math_2_ninv_mul_symm(dc_yb_2, dc_xa_2);
            dc_yred_2.math_2_add_mul(dc_ya_2, dc_xb_2, dc_nzbxa_2);
            math_2x2_mul(dc_nzbjb_2, dc_yb_2, dc_jb_2);
            math_2_add_mul_jred(dc_jred_2, dc_ja_2, dc_xb_2, dc_nzbjb_2);
        }
        else {
            throw hiba("cella::fw_dc", "B size 1 or 2 expected, %u arrived", B_2);
        }
    }
    is_change_adm = is_change_j = true; // TODO: a tényleges helyzetet kell elõállítani
}


//***********************************************************************
void faces_cella::bw_dc(uns melyik_fa) {
//***********************************************************************
    if (melyik_fa == 1) {
        if (B_1 == 1)        math_1_add_mul_ub(dc_ub_1, dc_nzbjb_1, dc_nzbxa_1, dc_UA_1);
        else if (B_1 == 2)   math_2_add_mul_ub(dc_ub_1, dc_nzbjb_1, dc_nzbxa_1, dc_UA_1);
        else throw hiba("cella::bw_dc", "B size 1 or 2 expected, %u arrived", B_1);
        if (akt_sim.p_akt_sim->nemlin_tipus == nmt_el_th_newton) {
            dc_UA_2.zero();
            dc_IA_2.zero();
            dc_ub_2.zero();
        }
    }
    else {
        if (B_2 == 1)        math_1_add_mul_ub(dc_ub_2, dc_nzbjb_2, dc_nzbxa_2, dc_UA_2);
        else if (B_2 == 2)   math_2_add_mul_ub(dc_ub_2, dc_nzbjb_2, dc_nzbxa_2, dc_UA_2);
        else throw hiba("cella::bw_dc", "B size 1 or 2 expected, %u arrived", B_2);
        if (akt_sim.p_akt_sim->nemlin_tipus == nmt_el_th_newton) {
            dc_UA_1.zero();
            dc_IA_1.zero();
            dc_ub_1.zero();
        }
    }
    // TODO: a face-tokokban a dUT-kbe és dIP-kbe kell tenni az értékeket. Kell ez?
    // TODO: csak akkor számoljuk újra az UT stb. értékeket, ha volt bw
    is_change_adm = false;
    is_change_j = false;
}


//***********************************************************************
void faces_cella::print_G() const {
//***********************************************************************
    for (uns i = 1; i < faces_dc.size(); i++)
        faces_dc[i].p_core->print_G();
}


//***********************************************************************
void faces_cella::debug_write(::std::ofstream & fs) const {
//***********************************************************************
    for (uns i = 1; i < faces_dc.size(); i++) {
        fs << "face " << ::std::setfill('0') << ::std::setw(2) << i <<": ";
        faces_dc.unsafe(i).p_core->debug_write(fs, true);
    }
    if(is_el) {
        fs << "face EC: ";
        el_center_face_dc.p_core->debug_write(fs, true);
    }
    if (is_th) {
        fs << "face TC: ";
        th_center_face_dc.p_core->debug_write(fs, true);
    }
    fs << "JA=\n";
    dc_ja_1.debug_write(fs);
    fs << "JB=\n";
    dc_jb_1.debug_write(fs);
    fs << "UA=\n";
    dc_UA_1.debug_write(fs);
    fs << "UB=\n";
    dc_ub_1.debug_write(fs);
    fs << "YA=\n";
    dc_ya_1.debug_write(fs);
    fs << "XA=\n";
    dc_xa_1.debug_write(fs);
    fs << "XB=\n";
    dc_xb_1.debug_write(fs);
    fs << "YB=\n";
    dc_yb_1.debug_write(fs);
    fs << "YRED=\n";
    dc_yred_1.debug_write(fs);
    fs << "JRED=\n";
    dc_jred_1.debug_write(fs);
}


//***********************************************************************
void faces_cella::randomize_zaj() {
//***********************************************************************
    for (uns i = 1; i < faces_dc.size(); i++){
        faces_dc.unsafe(i).p_core->randomize_zaj();
    }
    if (is_el)
        el_center_face_dc.p_core->randomize_zaj();
    if (is_th)
        th_center_face_dc.p_core->randomize_zaj();
}


//***********************************************************************
rvt faces_cella::cella_klaszter_2_NR_UT_dc(rvt & Uc, rvt & Tc, rvt & T_hiba) {
// A feszültséghibák MAXIMUMÁT! adja vissza
// Ha a középpontra fesz.gen kapcsolódik, 0-t adjon vissza! TODO: ezt egyelõre nem teszem bele. Lássuk a gyakorlatban, hogy kell-e.
//***********************************************************************
    rvt sum_hiba = rvt();
    for (uns i = 1; i < faces_dc.size(); i++)
        faces_dc.unsafe(i).p_core->update_UT_IP_NR();
    if (is_el) {
        el_center_face_dc.p_core->update_UT_IP_NR();
        sum_hiba += abs(*el_center_face_dc.p_UTi * akt_sim.alfa);
    }
    if (is_th) {
        th_center_face_dc.p_core->update_UT_IP_NR();
        sum_hiba += T_hiba = abs(*th_center_face_dc.p_UTi * akt_sim.alfa);
    }
    else {
        T_hiba = rvt();
    }
    Uc = is_el ? el_center_face_dc.emlekek.get_akt().UT : rvt();
    Tc = is_th ? th_center_face_dc.emlekek.get_akt().UT : rvt();
    return sum_hiba;
}


//***********************************************************************
rvt faces_cella::cella_klaszter_2_SA_UT_dc(rvt & T_hiba) {
// A feszültséghibák MAXIMUMÁT! adja vissza
// Ha a középpontra fesz.gen kapcsolódik, 0-t adjon vissza! TODO: ezt egyelõre nem teszem bele. Lássuk a gyakorlatban, hogy kell-e.
//***********************************************************************
    rvt sum_hiba = rvt();
    for (uns i = 1; i < faces_dc.size(); i++)
        faces_dc.unsafe(i).p_core->update_UT_IP_SA();
    if (is_el) {
        el_center_face_dc.p_core->update_UT_IP_SA();
        sum_hiba += abs(el_center_face_dc.emlekek.get_akt().UT - el_center_face_dc.emlekek.get_kiindulo().UT);
    }
    if (is_th) {
        th_center_face_dc.p_core->update_UT_IP_SA();
        sum_hiba += T_hiba = abs(th_center_face_dc.emlekek.get_akt().UT - th_center_face_dc.emlekek.get_kiindulo().UT);
    }
    else {
        T_hiba = rvt();
    }
    return sum_hiba;
}


//***********************************************************************
void faces_cella::cella_klaszter_3_parameterek_frissitese_dc() {
// Ha van hiszterézis, itt számoljuk újra, a középponti hõmérséklet meg-
// határozása után.
// A nem konstans face-ek a cella H-jával kell hívják a tulajdonságokat
//***********************************************************************
    if (p_cella_anyag->fazisvaltas.is_fazisvalto()) { // az akt-ba kellett kerülnie a most kiszámított középponti hõmérsékletnek.
        emlekek.get_akt().H = p_cella_anyag->fazisvaltas.H(get_Tc_akt(), emlekek.get_megtartando().H.ertek);
    }
    for (uns i = 1; i < faces_dc.size(); i++)
        faces_dc.unsafe(i).p_core->update_nemlinearis_parameterek();
    if(is_el)
        el_center_face_dc.p_core->update_nemlinearis_parameterek();
    if(is_th)
        th_center_face_dc.p_core->update_nemlinearis_parameterek();
    sum_disszipacio_dc = rvt();
    sum_ddissz_per_dT_dc = rvt();
    emlekek.get_akt().sum_lum = rvt();
    emlekek.get_akt().sum_rad = rvt();
    emlekek.get_akt().A_LED = rvt();
}


//***********************************************************************
void faces_cella::cella_klaszter_4_agaramok_szamitasa_dc() {
// disszipáció és a deriváltja is. Amelyik face disszipáló fajta, az kiszámolja a disszipációt is. 
// Csak akkor foglaljunk disszipáló face-t, ha az tényleg disszipál!
// Peremeknél a külsõ UT is!
//***********************************************************************
    for (uns i = 1; i < faces_dc.size(); i++)
        faces_dc.unsafe(i).p_core->agaram_szamitasa();
    if (is_el)
        el_center_face_dc.p_core->agaram_szamitasa();
    if (is_th) {
        th_center_face_dc.p_core->agaram_szamitasa();
        if (is_besugarzo) {
            for (uns i = 0; i < p_adatcella->besugarzo_cella.size(); i++) {
                const adat_besugarzo_cella & be = p_adatcella->besugarzo_cella[i];
                sum_disszipacio_dc += cellak[be.cella_index]->emlekek.get_elozo().sum_rad * be.arany.get_value(th_center_face_dc.emlekek.get_elozo().UT, emlekek.get_elozo().H).ertek;
                //printf("Pbe=%gW, arany=%g\n", cellak[be.cella_index].sum_rad.get_elozo(), be.arany.get_value(th_center_face_dc.UT.get_elozo(), H.get_elozo()).ertek);
            }
        }

        // a cellán kilépõ kék és sárga fény kiszámítása, ill a belõlük eredõ disszipáció figyelembe vétele

        rvt ossz_kek = 0;
        rvt ossz_sarga = 0;
        rvt ossz_uj_sarga = 0;
        rvt ossz_dissz = 0;
        rvt ossz_ki_kek = 0;
        rvt ossz_ki_sarga = 0;
        const rvt d_mul_b = akt_sim.p_akt_sim->dlppb;
        const rvt d_mul_y = akt_sim.p_akt_sim->dlppy;
        for (uns i = 0; i < p_adatcella->fenyutak.kek_indexek.size(); i++) {
            const auto & index = p_adatcella->fenyutak.kek_indexek[i];
            const auto & uut = akt_sim.p_akt_sim->fenyutak.kek_fenyutak[index.path_index];
            
            // junction cella

            rvt akt_kek = 0;
            if (cellak[uut.cells[0].cella_index]->get_tipus() == ct_faces_cella)
                akt_kek = static_cast<faces_cella*>(cellak[uut.cells[0].cella_index])->faces_dc[uut.face_index].emlekek.get_elozo().rad;
            else
                TODO("faces_cella::cella_klaszter_4_agaramok_szamitasa_dc() => utvonal ut");
            rvt akt_sarga = 0;
            akt_kek *= uut.K0; // a junction által sugárzott kék ekkora része jut erre a path-ra

            // köztes cellák

            for (uns j = 1; j < index.cella_index; j++) {
                const auto & kc = uut.cells[j];
                rvt be_kek = akt_kek;
                rvt be_sarga = akt_sarga;
                rvt d_akt = uut.d_mul * kc.ud;

                rvt T_koztes = 0;
                if(cellak[kc.cella_index]->get_tipus() == ct_faces_cella)
                    T_koztes = static_cast<faces_cella*>(cellak[kc.cella_index])->get_Tc_elozo(); // másik szálon is lehet
                else
                    TODO("faces_cella::cella_klaszter_4_agaramok_szamitasa_dc() => utvonal kc");
                rvt blue_abs_coeff_koztes = akt_sim.p_akt_sim->fenyutak.fenyut_tulajdonsagok[kc.blue_absorption_coeff_index].get_value(T_koztes, ertek_t()).ertek * d_mul_b;
                rvt conv_efficienc_koztes = akt_sim.p_akt_sim->fenyutak.fenyut_tulajdonsagok[kc.conversion_efficiency_index].get_value(T_koztes, ertek_t()).ertek;
                rvt yllw_abs_coeff_koztes = akt_sim.p_akt_sim->fenyutak.fenyut_tulajdonsagok[kc.yllw_absorption_coeff_index_sig / 4].get_value(T_koztes, ertek_t()).ertek * d_mul_y;
                rvt re_conv_effici_koztes = akt_sim.p_akt_sim->fenyutak.fenyut_tulajdonsagok[kc.re_conversion_efficiency_index].get_value(T_koztes, ertek_t()).ertek;

                rvt ki_kek = be_kek * exp(-blue_abs_coeff_koztes * d_akt);
                rvt elnyelt_kek = be_kek - ki_kek;
                rvt uj_ki_sarga = elnyelt_kek * conv_efficienc_koztes;
                rvt sarga_tovabb = be_sarga * exp(-yllw_abs_coeff_koztes * d_akt);
                rvt elnyelt_sarga = be_sarga - sarga_tovabb;
                rvt re_sarga = elnyelt_sarga * re_conv_effici_koztes;

                akt_kek = ki_kek * kc.uK / 32768;
                akt_sarga = uj_ki_sarga * kc.uUS / 32768 + sarga_tovabb * kc.uRS / 32768 + re_sarga * kc.uUS / 32768;
            }

            // aktuális cella

            const auto & ac = uut.cells[index.cella_index];
            rvt d_akt = uut.d_mul * ac.ud;
            rvt T_akt = get_Tc_akt();
            rvt akt_blue_abs_coeff = akt_sim.p_akt_sim->fenyutak.fenyut_tulajdonsagok[ac.blue_absorption_coeff_index].get_value(T_akt, ertek_t()).ertek * d_mul_b;
            rvt akt_conv_efficienc = akt_sim.p_akt_sim->fenyutak.fenyut_tulajdonsagok[ac.conversion_efficiency_index].get_value(T_akt, ertek_t()).ertek;
            rvt akt_yllw_abs_coeff = akt_sim.p_akt_sim->fenyutak.fenyut_tulajdonsagok[ac.yllw_absorption_coeff_index_sig / 4].get_value(T_akt, ertek_t()).ertek * d_mul_y;
            rvt akt_re_conv_effici = akt_sim.p_akt_sim->fenyutak.fenyut_tulajdonsagok[ac.re_conversion_efficiency_index].get_value(T_akt, ertek_t()).ertek;

            rvt akt_ki_kek = akt_kek * exp(-akt_blue_abs_coeff * d_akt);
            rvt akt_elnyelt_kek = akt_kek - akt_ki_kek;
            rvt akt_ki_uj_sarga = akt_elnyelt_kek * akt_conv_efficienc;
            rvt akt_sarga_tovabb = akt_sarga * exp(-akt_yllw_abs_coeff * d_akt);
            rvt akt_elnyelt_sarga = akt_sarga - akt_sarga_tovabb;
            rvt akt_re_sarga = akt_elnyelt_sarga * akt_re_conv_effici;

            if ((ac.yllw_absorption_coeff_index_sig & 3) == full_significant) ossz_kek += akt_ki_kek;
            else if ((ac.yllw_absorption_coeff_index_sig & 3) == half_significant) ossz_kek += 0.5*akt_ki_kek;
            if ((ac.yllw_absorption_coeff_index_sig & 3) == full_significant) ossz_sarga += akt_ki_uj_sarga * ac.uUS / 32768 + akt_sarga_tovabb * ac.uRS / 32768 + akt_re_sarga * ac.uUS / 32768;
            else if ((ac.yllw_absorption_coeff_index_sig & 3) == half_significant) ossz_sarga += 0.5*(akt_ki_uj_sarga * ac.uUS / 32768 + akt_sarga_tovabb* ac.uRS / 32768 + akt_re_sarga * ac.uUS / 32768);
            ossz_uj_sarga += akt_ki_uj_sarga + akt_re_sarga;
            ossz_dissz += akt_elnyelt_kek - akt_ki_uj_sarga + akt_elnyelt_sarga - akt_re_sarga;
            if (index.cella_index == uut.cells.size() - 1) {
                ossz_ki_kek += akt_ki_kek;
                //ossz_ki_kek += akt_ki_kek * ac.uK / 32768;
                //ossz_ki_sarga += akt_sarga_tovabb;
                //ossz_ki_sarga += akt_ki_uj_sarga + akt_sarga_tovabb + akt_re_sarga;
                ossz_ki_sarga += akt_ki_uj_sarga * ac.uUS / 32768 + akt_sarga_tovabb * ac.uRS / 32768 + akt_re_sarga * ac.uUS / 32768;
            }
        }
        for (uns i = 0; i < p_adatcella->fenyutak.sarga_indexek.size(); i++) {
            const auto & index = p_adatcella->fenyutak.sarga_indexek[i];
            const auto & uut = akt_sim.p_akt_sim->fenyutak.sarga_fenyutak[index.path_index];
            //printf("%u\n", akt_sim.p_akt_sim->fenyutak.sarga_fenyutak.size());

            // source cella

            rvt akt_sarga = 0;
            if (cellak[uut.cells[0].cella_index]->get_tipus() == ct_faces_cella)
                akt_sarga = static_cast<faces_cella*>(cellak[uut.cells[0].cella_index])->emlekek.get_elozo().sum_new_yellow;
            else
                TODO("faces_cella::cella_klaszter_4_agaramok_szamitasa_dc() => utvonal ut 2");
            akt_sarga *= uut.R0; // a source által sugárzott új sárga ekkora része jut erre a path-ra

            // köztes cellák

            for (uns j = 1; j < index.cella_index; j++) {
                const auto & kc = uut.cells[j];
                rvt be_sarga = akt_sarga;
                rvt d_akt = uut.d_mul * kc.ud;

                rvt T_koztes = 0;
                if (cellak[kc.cella_index]->get_tipus() == ct_faces_cella)
                    T_koztes = static_cast<faces_cella*>(cellak[kc.cella_index])->get_Tc_elozo(); // másik szálon is lehet
                else
                    TODO("faces_cella::cella_klaszter_4_agaramok_szamitasa_dc() => utvonal kc 2");
                
                rvt yllw_abs_coeff_koztes = akt_sim.p_akt_sim->fenyutak.fenyut_tulajdonsagok[kc.yllw_absorption_coeff_index_sig / 4].get_value(T_koztes, ertek_t()).ertek * d_mul_y;
                rvt regi_ki_sarga = be_sarga * exp(-yllw_abs_coeff_koztes *d_akt);
                akt_sarga = regi_ki_sarga * kc.uRS / 32768; // a reemitted sárga nem kerül be az aktuális sárga sugárba
            }

            // aktuális cella

            const auto & ac = uut.cells[index.cella_index];
            rvt d_akt = uut.d_mul * ac.ud;
            rvt T_akt = get_Tc_akt();
            rvt akt_yllw_abs_coeff = akt_sim.p_akt_sim->fenyutak.fenyut_tulajdonsagok[ac.yllw_absorption_coeff_index_sig / 4].get_value(T_akt, ertek_t()).ertek * d_mul_y;
            rvt akt_re_conv_effici = akt_sim.p_akt_sim->fenyutak.fenyut_tulajdonsagok[ac.re_conversion_efficiency_index].get_value(T_akt, ertek_t()).ertek;

            rvt akt_sarga_tovabb = akt_sarga * exp(-akt_yllw_abs_coeff * d_akt);
            rvt akt_elnyelt_sarga = akt_sarga - akt_sarga_tovabb;
            rvt akt_re_sarga = akt_elnyelt_sarga * akt_re_conv_effici;

            if ((ac.yllw_absorption_coeff_index_sig & 3) == full_significant) ossz_sarga += akt_sarga_tovabb;// +akt_re_sarga;
            else if ((ac.yllw_absorption_coeff_index_sig & 3) == half_significant) ossz_sarga += 0.5*(akt_sarga_tovabb/* + akt_re_sarga*/);
            ossz_uj_sarga += akt_re_sarga; // a reemitted sárga növeli a cella sugárzott teljesítményét
            ossz_dissz += akt_elnyelt_sarga - akt_re_sarga;
//            if ((ac.yllw_absorption_coeff_index_sig & 3) == full_significant && index.cella_index != uut.cells.size() - 1) {
//                printf("cella_index = %u, size = %u, sig = %u, s=%g\n", index.cella_index, uut.cells.size(), (ac.yllw_absorption_coeff_index_sig & 3), akt_sarga_tovabb);
//            }
            if (index.cella_index == uut.cells.size() - 1) {
                //ossz_ki_sarga += akt_sarga_tovabb;
                //ossz_ki_sarga += akt_sarga_tovabb * ac.uRS / 32768 + akt_re_sarga;
                ossz_ki_sarga += akt_sarga_tovabb;// +akt_re_sarga;
            }
        }
        emlekek.get_akt().sum_blue = emlekek.get_akt().sum_rad + ossz_kek; // a junction cella kék fényt sugároz
        emlekek.get_akt().sum_yellow = ossz_sarga;
        emlekek.get_akt().sum_new_yellow = ossz_uj_sarga;
        emlekek.get_akt().sum_blue_out = ossz_ki_kek;
        emlekek.get_akt().sum_yellow_out = ossz_ki_sarga;
        emlekek.get_akt().sum_diss_out = ossz_dissz;
        sum_disszipacio_dc += ossz_dissz;
    }
}


//***********************************************************************
rvt faces_cella::cella_klaszter_5_NR_hibaaramok_szamitasa_dc() {
// Az áramhibák ÖSSZEGÉT! adja vissza
// Ha a középpontra fesz.gen kapcsolódik, 0-t adjon vissza!
// is_change_j beállítása, ha változik
//***********************************************************************
    rvt sum_el_hiba = rvt(), sum_th_hiba = rvt();
    dc_ja_1.zero();
    dc_ja_2.zero();
    dc_jb_1.zero();
    dc_jb_2.zero();
    // TODO: ha van centroid, akkor annak mennyi az árama? 
    for (uns i = 1; i < faces_dc.size(); i++)
        faces_dc.unsafe(i).p_core->update_hibaaram_NR_kereszt(sum_el_hiba, sum_th_hiba);
    if (is_el) {
        el_center_face_dc.p_core->update_hibaaram_NR_kozepponti(sum_el_hiba, sum_th_hiba);
    }
    if (is_th) {
        emlekek.get_akt().sum_diss = sum_disszipacio_dc;
        sum_th_hiba -= sum_disszipacio_dc; // a disszipációt itt veszük figyelembe
        th_center_face_dc.p_core->update_hibaaram_NR_kozepponti(sum_el_hiba, sum_th_hiba);
    }
    //printf("%g\t%g\t", sum_el_hiba, sum_th_hiba);
    return abs(sum_el_hiba) + abs(sum_th_hiba);
}


//***********************************************************************
void faces_cella::cella_klaszter_5_SA_fill_J_dc() {
// is_change_j beállítása, ha változik
//***********************************************************************
    for (uns i = 1; i < faces_dc.size(); i++)
        faces_dc.unsafe(i).p_core->update_inhom_SA_kereszt();
    if (is_el)
        el_center_face_dc.p_core->update_inhom_SA_kozepponti(rvt());
    if (is_th)
        th_center_face_dc.p_core->update_inhom_SA_kozepponti(sum_disszipacio_dc);
}


//***********************************************************************
void faces_cella::fw_klaszter_1_update_Jakobi_dc() {
// A mátrixok konstans elemeit az update_for_uj_lepes tölti fel, a
// többit itt állítjuk be.
//***********************************************************************
    switch (feltoltes_tipus) {
        case ft_el:
            *p_Yee_dc = konst_Yee_dc;
            break;
        case ft_th:
            *p_Ytt_dc = konst_Ytt_dc - sum_ddissz_per_dT_dc; // különválasztott elektromos és termikus cellánál lehet termikus cella esetében is disszipáció
            break;
        case ft_elth_nincs_csatolt:
            *p_Yee_dc = konst_Yee_dc;
            *p_Ytt_dc = konst_Ytt_dc - sum_ddissz_per_dT_dc;
            break;
        case ft_elth_van_csatolt:
            *p_Yee_dc = konst_Yee_dc;
            *p_Yet_dc = konst_Yet_dc;
            *p_Yte_dc = konst_Yte_dc;
            *p_Ytt_dc = konst_Ytt_dc - sum_ddissz_per_dT_dc;
            break;
        default:
            throw hiba("fw_klaszter_1_update_Jakobi_dc", "unknown feltoltes_tipus");
            break;
    }
    if (akt_sim.p_akt_sim->nemlin_tipus == nmt_klasszik_iteracio) {
        for (uns i = 1; i < faces_dc.size(); i++)
            faces_dc.unsafe(i).p_core->update_admittanciamatrix_SA();
        if (is_el)
            el_center_face_dc.p_core->update_admittanciamatrix_SA();
        if (is_th)
            th_center_face_dc.p_core->update_admittanciamatrix_SA();
    }
    else if (feltoltes_tipus == ft_elth_van_csatolt) {
        for (uns i = 1; i < faces_dc.size(); i++)
            faces_dc.unsafe(i).p_core->update_Jakobi_NR_van_csatolt();
        if (is_el)
            el_center_face_dc.p_core->update_Jakobi_NR_van_csatolt();
        if (is_th)
            th_center_face_dc.p_core->update_Jakobi_NR_van_csatolt();
    }
    else {
        for (uns i = 1; i < faces_dc.size(); i++)
            faces_dc.unsafe(i).p_core->update_Jakobi_NR_nincs_csatolt();
        if (is_el)
            el_center_face_dc.p_core->update_Jakobi_NR_nincs_csatolt();
        if (is_th)
            th_center_face_dc.p_core->update_Jakobi_NR_nincs_csatolt();
    }
}

/*
//***********************************************************************
void cella::fw_klaszter_1_NR_update_Jakobi_Y_dc(bool force) {
//       if (force || !is_lin)...
//***********************************************************************
    TODO("fw_klaszter_1_NR_fill_Jakobi_Y_dc");
}


//***********************************************************************
void cella::fw_klaszter_1_SA_update_Y_dc(bool force) {
//       if (force || !is_lin)...
//***********************************************************************
    TODO("fw_klaszter_1_SA_fill_Y_dc");
}
*/

//***********************************************************************
void faces_cella::init_ac() {
//***********************************************************************
    TODO("AC init");
}


//***********************************************************************
void faces_cella::face_foglalas_ac() {
//***********************************************************************
    is_lin = true; // AC esetben mindig lin, ne változtasd!
    TODO("AC face foglalas");
}


//***********************************************************************
void faces_cella::gerj_update_ac() {
//***********************************************************************
    TODO("AC gerj update");
}


//***********************************************************************
void faces_cella::peremfeltetel_update_ac() {
//***********************************************************************
    TODO("AC perem update");
}


//***********************************************************************
void faces_cella::del_all_prev_ac() {
//***********************************************************************
    for (uns i = 0; i < faces_ac.size(); i++) {
        faces_ac[i].bele->emlekek.clear_megtartando_is();
    }
    el_center_face_ac.bele->emlekek.clear_megtartando_is();
    th_center_face_ac.bele->emlekek.clear_megtartando_is();
}


//***********************************************************************
void faces_cella::update_ac(bool force) {
//***********************************************************************
    TODO("AC update");
}


//***********************************************************************
void faces_cella::fw_ac(uns melyik_fa) {
// Ha AC-DC-AC váltás során DC-ben változik valami, akkor nem biztos, hogy
// ez a felépítés jól fog mûködni.
//***********************************************************************
    TODO("fw_ac");
}


//***********************************************************************
void faces_cella::bw_ac(uns melyik_fa) {
// A fa backwardja nem hívja meg, ha nem kell backward!
//***********************************************************************
    TODO("AC BW");
}


//***********************************************************************
void th6_cella::resize_matrixok_dc() {
//***********************************************************************
    
    // Ha a fa változik, de a cella nem, attól még kellhet újrafoglalás
    
    is_2_fa = akt_sim.p_akt_sim->db_fa == 2;
    uns uj_A_1 = is_2_fa ? 0 : 6; // (A_el + A_th)
    uns uj_A_2 = is_2_fa ? 6 : 0;
    uns uj_B_1 = is_2_fa ? 0 : 1;
    uns uj_B_2 = is_2_fa ? 1 : 0;

    bool is_kell_set_size = false;
    if (A_1 != uj_A_1 || A_2 != uj_A_2 || B_1 != uj_B_1 || B_2 != uj_B_2) {
        is_kell_set_size = true;
        A_1 = uj_A_1;
        A_2 = uj_A_2;
        B_1 = uj_B_1;
        B_2 = uj_B_2;
    }
    if (is_kell_set_size) {

        // mátrixok foglalása

        dc_yred_1.set_size(A_1, A_1);
        dc_yred_2.set_size(A_2, A_2);
        dc_jred_1.set_size(A_1);
        dc_jred_2.set_size(A_2);
        dc_UA_1.set_size(A_1);
        dc_UA_2.set_size(A_2);
        dc_IA_1.set_size(A_1);
        dc_IA_2.set_size(A_2);
    }
}


//***********************************************************************
void th6_cella::halozat_foglalas_cellaban_dc() {
//***********************************************************************
    const adat_cella & adatcella = *p_adatcella;
    set_hiba_hol h("th6_cella::halozat_foglalas_cellaban_dc");

    is_lin = true;
    p_cella_anyag = &akt_sim.p_akt_sim->anyagok[akt_sim.akt_anyag_index[adatcella.anyag_index]];

    if (p_cella_anyag->fazisvaltas.is_fazisvalto()) {
        is_lin = false;
    }

    A_top = 0;
    for (meret_t i = 0; i < p_adatcella->facek.size(); i++) {
        if (p_adatcella->facek.unsafe(i).oldal == pit_top) {
            A_top = akt_sim.p_akt_sim->meretek[p_adatcella->facek.unsafe(i).A_index]; // Pontosan egy top-ja lehet
        }
    }
}

}

