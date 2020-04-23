//***********************************************************************
// redukci�s fa template header
// Creation date:  2018. 08. 09.
// Creator:        Pohl L�szl�
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
    matrix<adattipus> YRED;                 // fw ki: nulladik szinten meg kell �rizze az el�z� iter�ci�ban bet�lt�tt �rt�keket (friss�t�s miatt)!
    vektor<adattipus> JRED;                 // fw ki: nulladik szinten meg kell �rizze az el�z� iter�ci�ban bet�lt�tt �rt�keket (friss�t�s miatt)!
    vektor<adattipus> UA;                   // bw be
    vektor<adattipus> IA, UB;               // bw ki, IA-t csak az elemi cell�khoz kell kisz�molni
    //***********************************************************************
    void init_fa_elem(const adat_fa_elem & faelemadat, uns melyikfa, vektor<redukcios_fa_elem<adattipus, is_ac> > & sajat_taroloja); // m�trixoknak mem�riafoglal�s + nullz�s + egy�b be�ll�t�sok
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
    void forward(const adat_fa_elem & faelemadat, uns melyikfa, vektor<redukcios_fa_elem<adattipus, is_ac> > & sajat_taroloja); // redukci�, ha kell + fw
    //***********************************************************************
    void backward();        // bw + cell�k eredm�nyeinek kisz�m�t�sa
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
    vektor<redukcios_fa_elem<adattipus, is_ac> > fa_1, fa_2; // a 0 index� is hasznos!
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        fa_1.set_size(akt_sim.p_akt_sim->gyoker_1_index); // A 0 index� nem dummy, �gy a sim->gyoker_1_index pont a sz�ks�ges m�ret !
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
    if (p_faelem_adat->cella_index != 0) { // cell�ra kapcsol�dik
        p_cella = &cellak[p_faelem_adat->cella_index];
        is_szimm = !p_cella->is_nonsymmetrical();
        // itt nem kell inicializ�lni a cell�t, mert a cella fw-je �gyis inicializ�lja, ha kell
        if (is_szimm) YRED.set_size_szimm(p_faelem_adat->A);
        else          YRED.set_size(p_faelem_adat->A, p_faelem_adat->A);
        JRED.set_size(p_faelem_adat->A);
        UA.set_size(p_faelem_adat->A);
        IA.set_size(p_faelem_adat->A);
    }
    else { // k�t faelemre kapcsol�dik
        p_bal  = &sajat_taroloja[p_faelem_adat->bal_elem_indexe - 1]; // A saj�t t�rol�ja 0-t�l indexel�dik, a faelem_adat 1-t�l adja
        p_jobb = &sajat_taroloja[p_faelem_adat->jobb_elem_indexe - 1]; // A saj�t t�rol�ja 0-t�l indexel�dik, a faelem_adat 1-t�l adja
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
        if (p_faelem_adat->B < 3) { // 1 �s 2 eset�n nem haszn�ljuk ezeket a seg�dm�trixokat
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
    
    // Ellen�rizz�k, hogy a k�t bemen� �s az aktu�lis faelem szimmetrikuss�ga megfelel�-e, �s ha nem, m�dos�tjuk a f�t.
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
                if (sor.is_balbol) { // m�soland�n�l csak az egyikb�l m�sol 
                    if (osz.is_masolando) {
                        if (osz.is_balbol)   YRED.submatrix_copy(p_bal->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.db, osz.db); // mindkett� m�soland�, YA-ba megy a blokk
                        else                 ; // mindkett� m�soland�, de nem ugyanaz a forr�s, �gy ilyenkor nincs m�velet
                    }
                    else {
                        if (osz.is_centroid) YRED.submatrix_copy(p_bal->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.db, osz.db); // sor m�soland�, oszlop centroid: m�sol�s YA-ba
                        else                 XB.submatrix_copy(p_bal->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.db, osz.db); // sor m�soland�, oszlop reduk�land�: XB-be m�solja
                    }
                }
                else { // sor: m�soland� jobb�l
                    if (osz.is_masolando) {
                        if (osz.is_balbol)   ; // mindkett� m�soland�, de nem ugyanz a forr�s, �gy ilyenkor nincs m�velet
                        else                 YRED.submatrix_copy(p_jobb->YRED, sor.hova, osz.hova, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // mindkett� m�soland�, YA-ba megy a blokk
                    }
                    else {
                        if (osz.is_centroid) YRED.submatrix_copy(p_jobb->YRED, sor.hova, osz.hova, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // sor m�soland�, oszlop centroid: m�sol�s YA-ba
                        else                 XB.submatrix_copy(p_jobb->YRED, sor.hova, osz.hova, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // sor m�soland�, oszlop reduk�land�: XB-be m�solja
                    }
                }
            }
            else { // a sor nem m�soland�
                if (sor.is_centroid) { 
                    if (osz.is_masolando) {
                        if (osz.is_balbol)   YRED.submatrix_copy(p_bal->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.db, osz.db); // sor centroid, oszlop m�soland�: YA-ba m�solja
                        else                 YRED.submatrix_copy(p_jobb->YRED, sor.hova, osz.hova, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // sor centroid, oszlop m�soland�: YA-ba m�solja
                    }
                    else { // oszlop nem m�soland�
                        if (osz.is_centroid) YRED.submatrix_add(p_bal->YRED, p_jobb->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // mindkett� centroid: YA-ba teszi a k�t blokk �sszeg�t
                        else                 XB.submatrix_add(p_bal->YRED, p_jobb->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // sor centroid, oszlop reduk�land�: XB-be teszi a k�t blokk �sszeg�t
                    }
                }
                else { // sor: reduk�land�
                    if (osz.is_masolando) {
                        if (osz.is_balbol)   XA.submatrix_copy(p_bal->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.db, osz.db); // sor reduk�land�, oszlop m�soland�: XA-ba m�solja
                        else                 XA.submatrix_copy(p_jobb->YRED, sor.hova, osz.hova, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // sor reduk�land�, oszlop m�soland�: XA-ba m�solja
                    }
                    else {
                        if (osz.is_centroid) XA.submatrix_add(p_bal->YRED, p_jobb->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // sor reduk�land�, oszlop cetroid: XA-ba teszi a k�t blokk �sszeg�t
                        else                 YB_NZB.submatrix_add(p_bal->YRED, p_jobb->YRED, sor.hova, osz.hova, sor.bal_from, osz.bal_from, sor.jobb_from, osz.jobb_from, sor.db, osz.db); // mindkett� reduk�land�: YB_NZB-be teszi a k�t blokk �sszeg�t
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
// bet�lti a cell�b�l YRED-et �s JRED-et friss�t�ssel
//***********************************************************************
    load_from_cella_belso(template_seged<is_ac>{});
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::load_from_cella_belso(template_seged<true>) {
// AC eset
//***********************************************************************
    p_cella->fw_ac(melyik_fa);

    // Ellen�rizz�k, hogy a cella �s a faelem szimmetrikuss�ga azonos-e, �s ha nem, m�dos�tjuk a f�t.
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
        if (p_cella->is_adm_changed()) { // csak akkor friss�t, ha a cella v�ltozott
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
        if (p_cella->is_j_changed()) { // csak akkor friss�t, ha a cella v�ltozott
            for (meret_t i = 0; i < JRED.size(); i++)
                is_J_updated = JRED.frissit_unsafe(i, p_cella->ac_jred_1.unsafe(i)) || is_J_updated;
        }
        is_kell_fw = is_J_updated || is_Y_updated;
    }
    else {
        egyforma_e_hiba(JRED.size(), p_cella->ac_jred_2.size(), "redukcios_fa_elem::load_from_cella JRED.size != jred.size");
        bool is_Y_updated = akt_sim.is_uj_fa_kell, is_J_updated = akt_sim.is_uj_fa_kell;
        if (p_cella->is_adm_changed()) { // csak akkor friss�t, ha a cella v�ltozott
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
        if (p_cella->is_j_changed()) { // csak akkor friss�t, ha a cella v�ltozott
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

    // Ellen�rizz�k, hogy a cella �s a faelem szimmetrikuss�ga azonos-e, �s ha nem, m�dos�tjuk a f�t.
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
        if (p_cella->is_adm_changed()) { // csak akkor friss�t, ha a cella v�ltozott
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
        if (p_cella->is_j_changed()) { // csak akkor friss�t, ha a cella v�ltozott
            for (meret_t i = 0; i < JRED.size(); i++)
                is_J_updated = JRED.frissit_unsafe(i, adattipus(p_cella->dc_jred_1.unsafe(i))) || is_J_updated;
        }
        is_kell_fw = is_J_updated || is_Y_updated;
    }
    else {
        egyforma_e_hiba(JRED.size(), p_cella->dc_jred_2.size(), "redukcios_fa_elem::load_from_cella JRED.size != jred.size");
        bool is_Y_updated = akt_sim.is_uj_fa_kell, is_J_updated = akt_sim.is_uj_fa_kell;
        if (p_cella->is_adm_changed()) { // csak akkor friss�t, ha a cella v�ltozott
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
        if (p_cella->is_j_changed()) { // csak akkor friss�t, ha a cella v�ltozott
            for (meret_t i = 0; i < JRED.size(); i++)
                is_J_updated = JRED.frissit_unsafe(i, adattipus(p_cella->dc_jred_2.unsafe(i))) || is_J_updated;
        }
        is_kell_fw = is_J_updated || is_Y_updated;
    }
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::store_to_cella() {
// Ez mindenk�pp m�sol, de a backsubs csak akkor h�vja, ha kell update
//***********************************************************************
    store_to_cella_belso(template_seged<is_ac>{});
}


//***********************************************************************
template<typename adattipus, bool is_ac>
inline void redukcios_fa_elem<adattipus, is_ac>::store_to_cella_belso(template_seged<true>) {
// Ez mindenk�pp m�sol, de a backsubs csak akkor h�vja, ha kell update
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
// Ez mindenk�pp m�sol, de a backsubs csak akkor h�vja, ha kell update
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

    if (p_cella != nullptr) { // cell�ra kapcsol�dik
        load_from_cella();
    }
    else { // k�t faelemre kapcsol�dik

        // redukci�

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
    if (p_faelem_adat->cella_index != 0) { // cell�ra kapcsol�dik
        if(is_szimm)    math_add_mul_symm(IA, JRED, YRED, UA);
        else            math_add_mul(IA, JRED, YRED, UA);
        store_to_cella();
        if (is_ac)  p_cella->bw_ac(melyik_fa);
        else        p_cella->bw_dc(melyik_fa);
    }
    else { // k�t faelemre kapcsol�dik
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