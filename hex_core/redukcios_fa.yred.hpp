//***********************************************************************
// redukciós fa template header
// Creation date:  2018. 08. 09.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef REDUKCIOS_FA_HEADER
#define	REDUKCIOS_FA_HEADER
//***********************************************************************


//***********************************************************************
#include "matrix.hpp"
#include "akt_sim_adatai.h"
#include "cella_es_face.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
template<bool B> struct template_seged {};
//***********************************************************************


//***********************************************************************
template<typename adattipus, bool is_ac> class redukcios_fa_elem {
//***********************************************************************

    //***********************************************************************
    const adat_fa_elem * p_faelem_adat;
    //***********************************************************************
    cella * p_cella;
    //***********************************************************************
    redukcios_fa_elem *p_bal, *p_jobb;
    //***********************************************************************
    uns melyik_fa;
    //***********************************************************************
    bool is_kell_redukcio, is_kell_fw;
    bool is_fw_kesz, is_bw_mehet, is_szimm;
    //***********************************************************************
    matrix<adattipus> YB_NZB, XA, XB;   // fw be
    vektor<adattipus> JA, JB;               // fw be
    matrix<adattipus> XAT, NZBXA, NZBXAT;   // fw temp
    vektor<adattipus> NZBJB;                // fw temp
    matrix<adattipus> YRED;                 // fw ki: nulladik szinten meg kell õrizze az elõzõ iterációban betöltött értékeket (frissítés miatt)!
    vektor<adattipus> JRED;                 // fw ki: nulladik szinten meg kell õrizze az elõzõ iterációban betöltött értékeket (frissítés miatt)!
    vektor<adattipus> UA;                   // bw be
    vektor<adattipus> IA, UB;               // bw ki, IA-t csak az elemi cellákhoz kell kiszámolni
    //***********************************************************************
    void init_fa_elem(const adat_fa_elem & faelemadat, uns melyikfa, vektor<redukcios_fa_elem<adattipus, is_ac> > & sajat_taroloja); // mátrixoknak memóriafoglalás + nullzás + egyéb beállítások
    //***********************************************************************
    void YA_YB_XA_XB_feltoltese();
    void JA_JB_feltoltese();
    void math_reduce_symm();
    void math_reduce_nonsymm();
    void math_forward();
    void load_from_cella();
    void load_from_cella_belso(template_seged<true>);
    void load_from_cella_belso(template_seged<false>);
    void store_to_cella();
    void store_to_cella_belso(template_seged<true>);
    void store_to_cella_belso(template_seged<false>);
    void vissza_copy();
public:
    //***********************************************************************
    redukcios_fa_elem() :p_faelem_adat{ nullptr }, melyik_fa{ 0 }, p_cella{ nullptr }, p_bal{ nullptr }, p_jobb{ nullptr },
        is_kell_redukcio{ true }, is_kell_fw{ true }, is_fw_kesz{ false }, is_bw_mehet{ false }, is_szimm{ false } {}
    //***********************************************************************
    void forward(const adat_fa_elem & faelemadat, uns melyikfa, vektor<redukcios_fa_elem<adattipus, is_ac> > & sajat_taroloja); // redukció, ha kell + fw
    //***********************************************************************
    void backward();        // bw + cellák eredményeinek kiszámítása
    //***********************************************************************

    //***********************************************************************
    bool is_fw_indithato(const adat_fa_elem & faelemadat, uns melyikfa, vektor<redukcios_fa_elem<adattipus, is_ac> > & sajat_taroloja) {
    //***********************************************************************
        if (p_faelem_adat == nullptr)
            init_fa_elem(faelemadat, melyikfa, sajat_taroloja);
        if (p_cella != nullptr) return true;
        return p_bal->is_fw_kesz && p_jobb->is_fw_kesz;
    }
    //***********************************************************************

    //***********************************************************************
    bool is_bw_indithato()const { return is_bw_mehet; }
    //***********************************************************************
    void debug_write(::std::ofstream & fs) const {
    //***********************************************************************
        fs << "YB_NZB=\n";
        YB_NZB.debug_write(fs);
        fs << "YRED=\n";
        YRED.debug_write(fs);
        fs << "JRED=\n";
        JRED.debug_write(fs);
    }

};


//***********************************************************************
template<typename adattipus, bool is_ac> class redukcios_fa {
//***********************************************************************

public:
    //***********************************************************************
    vektor<redukcios_fa_elem<adattipus, is_ac> > fa_1, fa_2; // a 0 indexû is hasznos!
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        fa_1.set_size(akt_sim.p_akt_sim->gyoker_1_index); // A 0 indexû nem dummy, így a sim->gyoker_1_index pont a szükséges méret !
        if (akt_sim.p_akt_sim->db_fa > 1) {
            fa_2.set_size(akt_sim.p_akt_sim->gyoker_2_index - akt_sim.p_akt_sim->gyoker_1_index);
        }
    }

    //***********************************************************************
    void clear() { fa_1.clear(); fa_2.clear(); }
    //***********************************************************************
};


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::init_fa_elem(const adat_fa_elem & faelemadat, uns melyikfa, vektor<redukcios_fa_elem<adattipus, is_ac> > & sajat_taroloja){
//***********************************************************************
    p_faelem_adat = &faelemadat;
    melyik_fa = melyikfa; // 1 vagy 2
    if (p_faelem_adat->cella_index != 0) { // cellára kapcsolódik
        p_cella = &cellak[p_faelem_adat->cella_index];
        is_szimm = !p_cella->is_nonsymmetrical();
        // itt nem kell inicializálni a cellát, mert a cella fw-je úgyis inicializálja, ha kell
        if (is_szimm) YRED.set_size_szimm(p_faelem_adat->A);
        else          YRED.set_size(p_faelem_adat->A, p_faelem_adat->A);
        JRED.set_size(p_faelem_adat->A);
        UA.set_size(p_faelem_adat->A);
        IA.set_size(p_faelem_adat->A);
    }
    else { // két faelemre kapcsolódik
        p_bal  = &sajat_taroloja[p_faelem_adat->bal_elem_indexe - 1]; // A saját tárolója 0-tól indexelõdik, a faelem_adat 1-tõl adja
        p_jobb = &sajat_taroloja[p_faelem_adat->jobb_elem_indexe - 1]; // A saját tárolója 0-tól indexelõdik, a faelem_adat 1-tõl adja
        is_szimm = p_bal->is_szimm && p_jobb->is_szimm;
        if (is_szimm) {
            YRED.set_size_szimm_and_zero(p_faelem_adat->A);
        }
        else {
            YRED.set_size_and_zero(p_faelem_adat->A, p_faelem_adat->A);
        }
        YB_NZB.set_size_and_zero(p_faelem_adat->B, p_faelem_adat->B);
        XA.set_size_and_zero(p_faelem_adat->B, p_faelem_adat->A);
        XB.set_size_and_zero(p_faelem_adat->A, p_faelem_adat->B);
        JA.set_size_and_zero(p_faelem_adat->A);
        JB.set_size_and_zero(p_faelem_adat->B);
        NZBXA.set_size(p_faelem_adat->B, p_faelem_adat->A);
        if (p_faelem_adat->B < 3) { // 1 és 2 esetén nem használjuk ezeket a segédmátrixokat
            XAT.clear();
            NZBXAT.clear();
        }
        else {
            if(is_szimm)    XAT.clear();
            else            XAT.set_size(p_faelem_adat->A, p_faelem_adat->B);
            NZBXAT.set_size(p_faelem_adat->A, p_faelem_adat->B);
        }
        NZBJB.set_size(p_faelem_adat->B);
        JRED.set_size(p_faelem_adat->A);
        UA.set_size(p_faelem_adat->A);
        UB.set_size(p_faelem_adat->B);
    }
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::YA_YB_XA_XB_feltoltese(){
//***********************************************************************
    
    // Ellenõrizzük, hogy a két bemenõ és az aktuális faelem szimmetrikussága megfelelõ-e, és ha nem, módosítjuk a fát.
    // Elvileg ez nagyon ritka eset

    if (is_szimm != p_bal->is_szimm && p_jobb->is_szimm) {
        is_szimm = !is_szimm;
        if (is_szimm) {
            YRED.set_size_szimm_and_zero(p_faelem_adat->A);
            XAT.clear();
        }
        else {
            YRED.set_size_and_zero(p_faelem_adat->A, p_faelem_adat->A);
            if (p_faelem_adat->B > 2) 
                XAT.set_size(p_faelem_adat->A, p_faelem_adat->B);
        }
    }

    const vektor<adat_masolando> & mas = p_faelem_adat->masolando;
    for (uns i = 0; i < mas.size(); i++) {
        const adat_masolando & sor = mas.unsafe(i);
        for (uns j = 0; j < mas.size(); j++) {
            const adat_masolando & osz = mas.unsafe(j);
            if (sor.is_masolando) { 
                if (sor.is_balbol) { // másolandónál csak az egyikbõl másol 
                    if (osz.is_masolando) {
                        if (osz.is_balbol)   YRED.submatrix_copy(p_bal->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.db, osz.db); // mindkettõ másolandó, YA-ba megy a blokk
                        else                 ; // mindkettõ másolandó, de nem ugyanaz a forrás, így ilyenkor nincs mûvelet
                    }
                    else {
                        if (osz.is_centroid) YRED.submatrix_copy(p_bal->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.db, osz.db); // sor másolandó, oszlop centroid: másolás YA-ba
                        else                 XB.submatrix_copy(p_bal->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.db, osz.db); // sor másolandó, oszlop redukálandó: XB-be másolja
                    }
                }
                else { // sor: másolandó jobból
                    if (osz.is_masolando) {
                        if (osz.is_balbol)   ; // mindkettõ másolandó, de nem ugyanz a forrás, így ilyenkor nincs mûvelet
                        else                 YRED.submatrix_copy(p_jobb->YRED, sor.hova, osz.hova, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // mindkettõ másolandó, YA-ba megy a blokk
                    }
                    else {
                        if (osz.is_centroid) YRED.submatrix_copy(p_jobb->YRED, sor.hova, osz.hova, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // sor másolandó, oszlop centroid: másolás YA-ba
                        else                 XB.submatrix_copy(p_jobb->YRED, sor.hova, osz.hova, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // sor másolandó, oszlop redukálandó: XB-be másolja
                    }
                }
            }
            else { // a sor nem másolandó
                if (sor.is_centroid) { 
                    if (osz.is_masolando) {
                        if (osz.is_balbol)   YRED.submatrix_copy(p_bal->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.db, osz.db); // sor centroid, oszlop másolandó: YA-ba másolja
                        else                 YRED.submatrix_copy(p_jobb->YRED, sor.hova, osz.hova, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // sor centroid, oszlop másolandó: YA-ba másolja
                    }
                    else { // oszlop nem másolandó
                        if (osz.is_centroid) YRED.submatrix_add(p_bal->YRED, p_jobb->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // mindkettõ centroid: YA-ba teszi a két blokk összegét
                        else                 XB.submatrix_add(p_bal->YRED, p_jobb->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // sor centroid, oszlop redukálandó: XB-be teszi a két blokk összegét
                    }
                }
                else { // sor: redukálandó
                    if (osz.is_masolando) {
                        if (osz.is_balbol)   XA.submatrix_copy(p_bal->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.db, osz.db); // sor redukálandó, oszlop másolandó: XA-ba másolja
                        else                 XA.submatrix_copy(p_jobb->YRED, sor.hova, osz.hova, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // sor redukálandó, oszlop másolandó: XA-ba másolja
                    }
                    else {
                        if (osz.is_centroid) XA.submatrix_add(p_bal->YRED, p_jobb->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // sor redukálandó, oszlop cetroid: XA-ba teszi a két blokk összegét
                        else                 YB_NZB.submatrix_add(p_bal->YRED, p_jobb->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // mindkettõ redukálandó: YB_NZB-be teszi a két blokk összegét
                    }
                }
            }
        }
    }
    printf("\n/*****/\n");
    p_bal->YRED.print_z();
    p_jobb->YRED.print_z();
    YRED.print_z();
    printf("\n/*****/\n");
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::JA_JB_feltoltese(){
//***********************************************************************
    const vektor<adat_masolando> & mas = p_faelem_adat->masolando;
    for (uns j = 0; j < mas.size(); j++) {
        const adat_masolando & osz = mas.unsafe(j);
        if (osz.is_masolando) {
            if (osz.is_balbol)  JA.subvektor_copy(p_bal->JRED, osz.hova, osz.bal_from, osz.db);
            else                JA.subvektor_copy(p_jobb->JRED, osz.hova, osz.jobb_from, osz.db); 
        }
        else {
            if (osz.is_centroid) JA.subvektor_add(p_bal->JRED, p_jobb->JRED, osz.hova, osz.bal_from, osz.jobb_from, osz.db);
            else                 JB.subvektor_add(p_bal->JRED, p_jobb->JRED, osz.hova, osz.bal_from, osz.jobb_from, osz.db);
        }
    }
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::math_reduce_symm(){
//***********************************************************************
    if (YB_NZB.get_col() == 1) {
        NZBXA.math_1_ninv_mul(YB_NZB, XA);
        YRED.math_1_add_mul_symm(YRED, XB, NZBXA);
    }
    else if (YB_NZB.get_col() == 2) {
        NZBXA.math_2_ninv_mul_symm(YB_NZB, XA);
        YRED.math_2_add_mul_symm(YRED, XB, NZBXA);
    }
    else {
        YB_NZB.math_symm_ninv_of_nonsymm();
        NZBXA.math_mul_t(YB_NZB, XB);
        NZBXAT.transp(NZBXA);
        YRED.math_add_mul_t_symm(YRED, XB, NZBXAT);
    }
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::math_reduce_nonsymm(){
//***********************************************************************
    if (YB_NZB.get_col() == 1) {
        NZBXA.math_1_ninv_mul(YB_NZB, XA);
        YRED.math_1_add_mul(YRED, XB, NZBXA);
    }
    else if (YB_NZB.get_col() == 2) {
        NZBXA.math_2_ninv_mul(YB_NZB, XA);
        YRED.math_2_add_mul(YRED, XB, NZBXA);
    }
    else {
        YB_NZB.math_ninv_np();
        XAT.transp(XA);
        NZBXA.math_mul_t(YB_NZB, XAT);
        NZBXAT.transp(NZBXA);
        YRED.math_add_mul_t(YRED, XB, NZBXAT);
    }
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::math_forward(){
//***********************************************************************
    if (YB_NZB.get_col() == 1) {
        math_1x1_mul(NZBJB, YB_NZB, JB);
        math_1_add_mul_jred(JRED, JA, XB, NZBJB);
    }
    else if (YB_NZB.get_col() == 2) {
        math_2x2_mul(NZBJB, YB_NZB, JB);
        math_2_add_mul_jred(JRED, JA, XB, NZBJB);
    }
    else {
        math_mul(NZBJB, YB_NZB, JB);
        math_add_mul(JRED, JA, XB, NZBJB);
    }
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::load_from_cella(){
// betölti a cellából YRED-et és JRED-et frissítéssel
//***********************************************************************
    load_from_cella_belso(template_seged<is_ac>{});
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::load_from_cella_belso(template_seged<true>) {
// AC eset
//***********************************************************************
    p_cella->fw_ac(melyik_fa);

    // Ellenõrizzük, hogy a cella és a faelem szimmetrikussága azonos-e, és ha nem, módosítjuk a fát.
    // Elvileg ez nagyon ritka eset

    if (p_cella->is_nonsymmetrical() == is_szimm) {
        is_szimm = !is_szimm;
        if (is_szimm) YRED.set_size_szimm(p_faelem_adat->A);
        else          YRED.set_size(p_faelem_adat->A, p_faelem_adat->A);
        if (!p_cella->is_adm_changed())
            throw hiba("redukcios_fa_elem::load_from_cella_belso", "symm/nonsymm change without adm mx change");
    }

    if (melyik_fa == 1) {
        egyforma_e_hiba(JRED.size(), p_cella->ac_jred_1.size(), "redukcios_fa_elem::load_from_cella JRED.size != jred.size");
        bool is_Y_updated = akt_sim.is_uj_fa_kell, is_J_updated = akt_sim.is_uj_fa_kell;
        if (p_cella->is_adm_changed()) { // csak akkor frissít, ha a cella változott
            if (is_szimm) {
                for (meret_t row = 0; row < YRED.get_row(); row++)
                    for (meret_t col = row; col < YRED.get_col(); col++)
                        is_Y_updated = YRED.frissit_unsafe(row, col, p_cella->ac_yred_1.get_elem_unsafe(row, col)) || is_Y_updated;
            }
            else {
                for (meret_t i = 0; i < YRED.size(); i++)
                    is_Y_updated = YRED.frissit_unsafe(i, p_cella->ac_yred_1.get_elem_unsafe(i)) || is_Y_updated;
            }
        }
        is_kell_redukcio = is_Y_updated;
        if (p_cella->is_j_changed()) { // csak akkor frissít, ha a cella változott
            for (meret_t i = 0; i < JRED.size(); i++)
                is_J_updated = JRED.frissit_unsafe(i, p_cella->ac_jred_1.unsafe(i)) || is_J_updated;
        }
        is_kell_fw = is_J_updated || is_Y_updated;
    }
    else {
        egyforma_e_hiba(JRED.size(), p_cella->ac_jred_2.size(), "redukcios_fa_elem::load_from_cella JRED.size != jred.size");
        bool is_Y_updated = akt_sim.is_uj_fa_kell, is_J_updated = akt_sim.is_uj_fa_kell;
        if (p_cella->is_adm_changed()) { // csak akkor frissít, ha a cella változott
            if (is_szimm) {
                for (meret_t row = 0; row < YRED.get_row(); row++)
                    for (meret_t col = row; col < YRED.get_col(); col++)
                        is_Y_updated = YRED.frissit_unsafe(row, col, p_cella->ac_yred_2.get_elem_unsafe(row, col)) || is_Y_updated;
            }
            else {
                for (meret_t i = 0; i < YRED.size(); i++)
                    is_Y_updated = YRED.frissit_unsafe(i, p_cella->ac_yred_2.get_elem_unsafe(i)) || is_Y_updated;
            }
        }
        is_kell_redukcio = is_Y_updated;
        if (p_cella->is_j_changed()) { // csak akkor frissít, ha a cella változott
            for (meret_t i = 0; i < JRED.size(); i++)
                is_J_updated = JRED.frissit_unsafe(i, p_cella->ac_jred_2.unsafe(i)) || is_J_updated;
        }
        is_kell_fw = is_J_updated || is_Y_updated;
    }
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::load_from_cella_belso(template_seged<false>) {
// DC eset
//***********************************************************************
    p_cella->fw_dc(melyik_fa);

    // Ellenõrizzük, hogy a cella és a faelem szimmetrikussága azonos-e, és ha nem, módosítjuk a fát.
    // Elvileg ez nagyon ritka eset

    if (p_cella->is_nonsymmetrical() == is_szimm) {
        is_szimm = !is_szimm;
        if (is_szimm) YRED.set_size_szimm(p_faelem_adat->A);
        else          YRED.set_size(p_faelem_adat->A, p_faelem_adat->A);
        if (!p_cella->is_adm_changed())
            throw hiba("redukcios_fa_elem::load_from_cella_belso", "symm/nonsymm change without adm mx change");
    }

    if (melyik_fa == 1) {
        egyforma_e_hiba(JRED.size(), p_cella->dc_jred_1.size(), "redukcios_fa_elem::load_from_cella JRED.size != jred.size");
        bool is_Y_updated = akt_sim.is_uj_fa_kell, is_J_updated = akt_sim.is_uj_fa_kell;
        if (p_cella->is_adm_changed()) { // csak akkor frissít, ha a cella változott
            if (is_szimm) {
                for (meret_t row = 0; row < YRED.get_row(); row++) 
                    for (meret_t col = row; col < YRED.get_col(); col++)
                        is_Y_updated = YRED.frissit_unsafe(row, col, adattipus(p_cella->dc_yred_1.get_elem_unsafe(row, col))) || is_Y_updated;
            }
            else {
                for (meret_t i = 0; i < YRED.size(); i++)
                    is_Y_updated = YRED.frissit_unsafe(i, adattipus(p_cella->dc_yred_1.get_elem_unsafe(i))) || is_Y_updated;
            }
        }
        is_kell_redukcio = is_Y_updated;
        if (p_cella->is_j_changed()) { // csak akkor frissít, ha a cella változott
            for (meret_t i = 0; i < JRED.size(); i++)
                is_J_updated = JRED.frissit_unsafe(i, adattipus(p_cella->dc_jred_1.unsafe(i))) || is_J_updated;
        }
        is_kell_fw = is_J_updated || is_Y_updated;
    }
    else {
        egyforma_e_hiba(JRED.size(), p_cella->dc_jred_2.size(), "redukcios_fa_elem::load_from_cella JRED.size != jred.size");
        bool is_Y_updated = akt_sim.is_uj_fa_kell, is_J_updated = akt_sim.is_uj_fa_kell;
        if (p_cella->is_adm_changed()) { // csak akkor frissít, ha a cella változott
            if (is_szimm) {
                for (meret_t row = 0; row < YRED.get_row(); row++)
                    for (meret_t col = row; col < YRED.get_col(); col++)
                        is_Y_updated = YRED.frissit_unsafe(row, col, adattipus(p_cella->dc_yred_2.get_elem_unsafe(row, col))) || is_Y_updated;
            }
            else {
                for (meret_t i = 0; i < YRED.size(); i++)
                    is_Y_updated = YRED.frissit_unsafe(i, adattipus(p_cella->dc_yred_2.get_elem_unsafe(i))) || is_Y_updated;
            }
        }
        is_kell_redukcio = is_Y_updated;
        if (p_cella->is_j_changed()) { // csak akkor frissít, ha a cella változott
            for (meret_t i = 0; i < JRED.size(); i++)
                is_J_updated = JRED.frissit_unsafe(i, adattipus(p_cella->dc_jred_2.unsafe(i))) || is_J_updated;
        }
        is_kell_fw = is_J_updated || is_Y_updated;
    }
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::store_to_cella() {
// Ez mindenképp másol, de a backsubs csak akkor hívja, ha kell update
//***********************************************************************
    store_to_cella_belso(template_seged<is_ac>{});
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::store_to_cella_belso(template_seged<true>) {
// Ez mindenképp másol, de a backsubs csak akkor hívja, ha kell update
//***********************************************************************
    if (melyik_fa == 1) {
        egyforma_e_hiba(UA.size(), p_cella->ac_UA_1.size(), "redukcios_fa_elem::store_to_cella UA size");
        for (meret_t i = 0; i < UA.size(); i++)
            p_cella->ac_UA_1.unsafe(i) = UA.unsafe(i);
        for (meret_t i = 0; i < IA.size(); i++)
            p_cella->ac_IA_1.unsafe(i) = IA.unsafe(i);
    }
    else {
        egyforma_e_hiba(UA.size(), p_cella->ac_UA_2.size(), "redukcios_fa_elem::store_to_cella UA size");
        for (meret_t i = 0; i < UA.size(); i++)
            p_cella->ac_UA_2.unsafe(i) = UA.unsafe(i);
        for (meret_t i = 0; i < IA.size(); i++)
            p_cella->ac_IA_2.unsafe(i) = IA.unsafe(i);
    }
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::store_to_cella_belso(template_seged<false>) {
// Ez mindenképp másol, de a backsubs csak akkor hívja, ha kell update
//***********************************************************************
    if (melyik_fa == 1) {
        egyforma_e_hiba(UA.size(), p_cella->dc_UA_1.size(), "redukcios_fa_elem::store_to_cella UA size");
        for (meret_t i = 0; i < UA.size(); i++)
            p_cella->dc_UA_1.unsafe(i) = rvt(UA.unsafe(i));
        for (meret_t i = 0; i < IA.size(); i++)
            p_cella->dc_IA_1.unsafe(i) = rvt(IA.unsafe(i));
    }
    else {
        egyforma_e_hiba(UA.size(), p_cella->dc_UA_2.size(), "redukcios_fa_elem::store_to_cella UA size");
        for (meret_t i = 0; i < UA.size(); i++)
            p_cella->dc_UA_2.unsafe(i) = rvt(UA.unsafe(i));
        for (meret_t i = 0; i < IA.size(); i++)
            p_cella->dc_IA_2.unsafe(i) = rvt(IA.unsafe(i));
    }
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::vissza_copy() {
//***********************************************************************
    const vektor<adat_masolando> & mas = p_faelem_adat->masolando;
    for (uns j = 0; j < mas.size(); j++) {
        const adat_masolando & osz = mas.unsafe(j);
        if (osz.is_masolando) {
            if (osz.is_balbol)  
                p_bal ->UA.subvektor_copy(UA, osz.bal_from,  osz.hova, osz.db);
            else                
                p_jobb->UA.subvektor_copy(UA, osz.jobb_from, osz.hova, osz.db);
        }
        else {
            if (osz.is_centroid) {
                p_bal ->UA.subvektor_copy(UA, osz.bal_from,  osz.hova, osz.db);
                p_jobb->UA.subvektor_copy(UA, osz.jobb_from, osz.hova, osz.db);
            }
            else {
                p_bal ->UA.subvektor_copy(UB, osz.bal_from,  osz.hova, osz.db);
                p_jobb->UA.subvektor_copy(UB, osz.jobb_from, osz.hova, osz.db);
            }
        }
    }
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::forward(const adat_fa_elem & faelemadat, uns melyikfa, vektor<redukcios_fa_elem<adattipus, is_ac> > & sajat_taroloja){
//***********************************************************************
    if (p_faelem_adat == nullptr)
        init_fa_elem(faelemadat, melyikfa, sajat_taroloja);

    if (p_cella != nullptr) { // cellára kapcsolódik
        load_from_cella();
    }
    else { // két faelemre kapcsolódik

        // redukció

        is_kell_redukcio = p_bal->is_kell_redukcio || p_jobb->is_kell_redukcio;
        if (is_kell_redukcio) {
            YA_YB_XA_XB_feltoltese();
            if (is_szimm)   math_reduce_symm();
            else            math_reduce_nonsymm();
        }

        // forwsubs

        is_kell_fw = p_bal->is_kell_fw || p_jobb->is_kell_fw;
        if (is_kell_fw) {
            JA_JB_feltoltese();
            math_forward();
        }
        //printf("\n\n");
        //YRED.print();
        //printf("\n\n");
        //JRED.print();
        //printf("\n\n");
        //getchar();
    }
    is_bw_mehet = false;
    is_fw_kesz = true;
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::backward(){
//***********************************************************************
    if (p_faelem_adat->cella_index != 0) { // cellára kapcsolódik
        if(is_szimm)    math_add_mul_symm(IA, JRED, YRED, UA);
        else            math_add_mul(IA, JRED, YRED, UA);
        store_to_cella();
        if (is_ac)  p_cella->bw_ac(melyik_fa);
        else        p_cella->bw_dc(melyik_fa);
    }
    else { // két faelemre kapcsolódik
        if (UB.size() == 1)      math_1_add_mul_ub(UB, NZBJB, NZBXA, UA);
        else if (UB.size() == 2) math_2_add_mul_ub(UB, NZBJB, NZBXA, UA);
        else                     math_add_mul(UB, NZBJB, NZBXA, UA);
        vissza_copy();
        p_bal->is_fw_kesz = false;
        p_bal->is_bw_mehet = true;
        p_jobb->is_fw_kesz = false;
        p_jobb->is_bw_mehet = true;
    }
    is_fw_kesz = false;
}


}

#endif