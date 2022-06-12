//***********************************************************************
// vektor template header
// Creation date:  2018. 08. 08.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef MATRIX_HEADER
#define	MATRIX_HEADER
//***********************************************************************


//***********************************************************************
#include "vektor.hpp"
#include <thread>
#include <mutex>
#include <chrono>
#include <atomic>
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
//#define MIN_MUL_ROW 64
//#define MIN_MUL_COL 64
//#define MIN_MUL_COM 64
//***********************************************************************


//***********************************************************************
template<typename adattipus> class parallel_seged;
template<typename adattipus> class matrix;
//***********************************************************************


//***********************************************************************
template<typename adattipus> struct hivando_matrixfuggveny {
//***********************************************************************

    //***********************************************************************
    enum fuggveny_tipus {
    //***********************************************************************
        fvt_none,
        fvt_math_add_mul_t_biztos,
        fvt_math_add_mul_t_symm,
        fvt_math_mul_t_biztos,
        fvt_math_ninv_np,
        fvt_math_nmul_t_biztos,
        fvt_math_sub_mul_t_biztos,
        fvt_math_sub_mul_t_symm_in_nonsymm,
        fvt_math_symm_ninv_of_nonsymm,
        fvt_math_bigblock_ninv_np = 100, // A bigblock azonosÍtója 100-nál kezdõdik
        fvt_math_bigblock_nonsymm_add_mul_t,
        fvt_math_bigblock_nonsymm_mul_t,
        fvt_math_bigblock_nonsymm_nmul_t,
        fvt_math_bigblock_nonsymm_sub_mul_t,
        fvt_math_bigblock_symm_add_mul_t,
        fvt_math_bigblock_symm_in_nonsymm_ninv_np,
        fvt_math_bigblock_symm_in_nosymm_sub_mul_t
    };

    //***********************************************************************
    matrix<adattipus> *p_dest;
    const matrix<adattipus> *p_src1, *p_src2, *p_src3;
    fuggveny_tipus tipus;
    ::std::condition_variable * p_kesz_egy_feladat;
    std::atomic<bool> * p_is_kesz;
    uns s_osztasszam; 
    bool s_is_priority;
    bool s_is_symmetrize_needed;
    //***********************************************************************

    //***********************************************************************
    hivando_matrixfuggveny() :p_dest{ nullptr }, p_src1{ nullptr }, p_src2{ nullptr }, p_src3{ nullptr },
        tipus{ fvt_none }, p_kesz_egy_feladat{ nullptr }, p_is_kesz{ nullptr }, s_osztasszam{ 0 }, 
        s_is_priority{ false }, s_is_symmetrize_needed{ false } {}
    //***********************************************************************

    //***********************************************************************
    void math_add_mul_t_biztos(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        const matrix<adattipus> & c, const matrix<adattipus> & a, const matrix<adattipus> & b_t, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        s_is_priority = is_priority;
        p_dest = &dest;
        p_src1 = &c;
        p_src2 = &a;
        p_src3 = &b_t;
        tipus = fvt_math_add_mul_t_biztos;
    }

    //***********************************************************************
    void math_add_mul_t_symm(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        const matrix<adattipus> & ya, const matrix<adattipus> & xb, const matrix<adattipus> & nzbxat, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        s_is_priority = is_priority;
        p_dest = &dest;
        p_src1 = &ya;
        p_src2 = &xb;
        p_src3 = &nzbxat;
        tipus = fvt_math_add_mul_t_symm;
    }

    //***********************************************************************
    void math_bigblock_ninv_np(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        uns osztasszam, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        s_is_priority = is_priority;
        p_dest = &dest;
        s_osztasszam = osztasszam;
        tipus = fvt_math_bigblock_ninv_np;
    }
    //***********************************************************************
    void math_bigblock_nonsymm_add_mul_t(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        uns osztasszam, const matrix<adattipus> & src1, const matrix<adattipus> & src2, const matrix<adattipus> & src3, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        p_dest = &dest;
        s_osztasszam = osztasszam;
        s_is_priority = is_priority;
        p_src1 = &src1;
        p_src2 = &src2;
        p_src3 = &src3;
        tipus = fvt_math_bigblock_nonsymm_add_mul_t;
    }

    //***********************************************************************
    void math_bigblock_nonsymm_mul_t(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        uns osztasszam, const matrix<adattipus> & src1, const matrix<adattipus> & src2, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        p_dest = &dest;
        s_osztasszam = osztasszam;
        s_is_priority = is_priority;
        p_src1 = &src1;
        p_src2 = &src2;
        tipus = fvt_math_bigblock_nonsymm_mul_t;
    }

    //***********************************************************************
    void math_bigblock_nonsymm_nmul_t(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        uns osztasszam, const matrix<adattipus> & src1, const matrix<adattipus> & src2, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        p_dest = &dest;
        s_osztasszam = osztasszam;
        s_is_priority = is_priority;
        p_src1 = &src1;
        p_src2 = &src2;
        tipus = fvt_math_bigblock_nonsymm_nmul_t;
    }

    //***********************************************************************
    void math_bigblock_nonsymm_sub_mul_t(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        uns osztasszam, const matrix<adattipus> & src1, const matrix<adattipus> & src2, const matrix<adattipus> & src3, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        p_dest = &dest;
        s_osztasszam = osztasszam;
        s_is_priority = is_priority;
        p_src1 = &src1;
        p_src2 = &src2;
        p_src3 = &src3;
        tipus = fvt_math_bigblock_nonsymm_sub_mul_t;
    }

    //***********************************************************************
    void math_bigblock_symm_add_mul_t(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        uns osztasszam, const matrix<adattipus> & src1, const matrix<adattipus> & src2, const matrix<adattipus> & src3, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        p_dest = &dest;
        s_osztasszam = osztasszam;
        s_is_priority = is_priority;
        p_src1 = &src1;
        p_src2 = &src2;
        p_src3 = &src3;
        tipus = fvt_math_bigblock_symm_add_mul_t;
    }

    //***********************************************************************
    void math_bigblock_symm_in_nonsymm_ninv_np(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        uns osztasszam, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        s_is_priority = is_priority;
        p_dest = &dest;
        s_osztasszam = osztasszam;
        tipus = fvt_math_bigblock_symm_in_nonsymm_ninv_np;
    }

    //***********************************************************************
    void math_bigblock_symm_in_nosymm_sub_mul_t(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        uns osztasszam, const matrix<adattipus> & src1,
        const matrix<adattipus> & src2, const matrix<adattipus> & src3, bool is_symmetrize_needed, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        s_is_priority = is_priority;
        p_dest = &dest;
        p_src1 = &src1;
        p_src2 = &src2;
        p_src3 = &src3;
        s_osztasszam = osztasszam;
        s_is_symmetrize_needed = is_symmetrize_needed;
        tipus = fvt_math_bigblock_symm_in_nosymm_sub_mul_t;
    }

    //***********************************************************************
    void math_mul_t_biztos(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        const matrix<adattipus> & a, const matrix<adattipus> & b_t, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        s_is_priority = is_priority;
        p_dest = &dest;
        p_src1 = &a;
        p_src2 = &b_t;
        tipus = fvt_math_mul_t_biztos;
    }

    //***********************************************************************
    void math_ninv_np(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        s_is_priority = is_priority;
        p_dest = &dest;
        tipus = fvt_math_ninv_np;
    }

    //***********************************************************************
    void math_nmul_t_biztos(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        const matrix<adattipus> & a, const matrix<adattipus> & b_t, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        s_is_priority = is_priority;
        p_dest = &dest;
        p_src1 = &a;
        p_src2 = &b_t;
        tipus = fvt_math_nmul_t_biztos;
    }

    //***********************************************************************
    void math_sub_mul_t_biztos(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        const matrix<adattipus> & c, const matrix<adattipus> & a, const matrix<adattipus> & b_t, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        s_is_priority = is_priority;
        p_dest = &dest;
        p_src1 = &c;
        p_src2 = &a;
        p_src3 = &b_t;
        tipus = fvt_math_sub_mul_t_biztos;
    }

    //***********************************************************************
    void math_sub_mul_t_symm_in_nonsymm(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest,
        const matrix<adattipus> & c, const matrix<adattipus> & a, const matrix<adattipus> & b, bool is_symmetrize_needed, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        s_is_priority = is_priority;
        p_dest = &dest;
        p_src1 = &c;
        p_src2 = &a;
        p_src3 = &b;
        s_is_symmetrize_needed = is_symmetrize_needed;
        tipus = fvt_math_sub_mul_t_symm_in_nonsymm;
    }

    //***********************************************************************
    void math_symm_ninv_of_nonsymm(std::atomic<bool> & is_kesz, ::std::condition_variable & kesz_egy_feladat, matrix<adattipus> & dest, bool is_priority) {
    //***********************************************************************
        is_kesz = false;
        p_is_kesz = &is_kesz;
        p_kesz_egy_feladat = &kesz_egy_feladat;
        s_is_priority = is_priority;
        p_dest = &dest;
        tipus = fvt_math_symm_ninv_of_nonsymm;
    }

    //***********************************************************************
    void run_job();
    //***********************************************************************
};


//***********************************************************************
template<typename adattipus> class matrix {
//***********************************************************************
    
    //***********************************************************************
    vektor<adattipus> t;
    vektor< vektor<adattipus> > sorok;
    parallel_seged<adattipus> *p_parallel_seged;
    meret_t row, col;
    bool is_szimm;
    //***********************************************************************

    //***********************************************************************
    void set_sorok() {
    //***********************************************************************
        sorok.set_size(row);
        if (is_szimm) {
            meret_t kezd = 0, db = row;
            for (meret_t i = 0; i < row; i++) {
                sorok.unsafe(i).rafektet(t, kezd - i, col); // így az m[i][j] valóban az i,j indexû elemet adja, de j>=i kötelezõ! 
                kezd += db;
                db--;
            }
        }
        else {
            for (meret_t i = 0; i < row; i++)
                sorok.unsafe(i).rafektet(t, i*col, col);
        }
    }

    //***********************************************************************
    matrix(const matrix&) = delete;
    matrix(const matrix&&) = delete;
    void operator=(const matrix&&) = delete;
    //***********************************************************************
public:
    
    //***********************************************************************
    const matrix & operator=(const matrix & masik) {
    // Nem ráfektetésbiztos.
    //***********************************************************************
        t = masik.t;
        row = masik.row;
        col = masik.col;
        is_szimm = masik.is_szimm;
        set_sorok();
        return *this;
    }
    //***********************************************************************
    matrix() :p_parallel_seged{ nullptr }, row{ 0 }, col{ 0 }, is_szimm{ false } {}
    //***********************************************************************
    ~matrix();
    //***********************************************************************
    meret_t get_col() const { return col; }
    //***********************************************************************
    meret_t get_row() const { return row; }
    //***********************************************************************
    bool get_is_szimm() const { return is_szimm; }
    //***********************************************************************
    meret_t size() const { return t.size(); }
    //***********************************************************************
    void clear() {
    //***********************************************************************
        t.clear();
        sorok.clear();
        row = col = 0;
        is_szimm = false;
        delete p_parallel_seged;
        p_parallel_seged = nullptr;
    }
    //***********************************************************************
    bool frissit_unsafe(meret_t i, const adattipus & uj) { return t.frissit_unsafe(i, uj); }
    //***********************************************************************
    bool frissit_unsafe(meret_t row, meret_t col, const adattipus & uj) { return sorok.unsafe(row).frissit_unsafe(col, uj); }
    //***********************************************************************
    void set_size(meret_t uj_row, meret_t uj_col) { clear(); row = uj_row; col = uj_col; is_szimm = false; t.set_size(row*col); set_sorok(); }
    //***********************************************************************
    void set_size_szimm(meret_t uj_rowcol) { clear(); row = uj_rowcol; col = uj_rowcol; is_szimm = true; t.set_size((uj_rowcol*(uj_rowcol + 1)) / 2); set_sorok(); }
    //***********************************************************************
    void set_size_and_zero(meret_t uj_r, meret_t uj_c) { set_size(uj_r, uj_c); zero_nembiztos(); }
    //***********************************************************************
    void set_size_szimm_and_zero(meret_t uj_rowcol) { set_size_szimm(uj_rowcol); zero_nembiztos(); }
    //***********************************************************************
    void rafektet(matrix & masik, uns start_row, uns start_col, uns row_db, uns col_db, bool is_symm) {
    // t üres lesz, így azok a tagfv-ek, amelyek t-t használják, nem mûködnek
    // a sorok használatával címezhetõk a cellák
    //***********************************************************************
        kisebb_e_hiba(start_row + row_db - 1, masik.row, "matrix::rafektet start_row + row_db");
        kisebb_e_hiba(start_col + col_db - 1, masik.col, "matrix::rafektet start_col + col_db");
        clear();
        sorok.set_size(row_db);
        row = row_db;
        col = col_db;
        is_szimm = is_symm;
        for (uns i = 0; i < row_db; i++)
            sorok.unsafe(i).rafektet(masik.sorok.unsafe(start_row + i), start_col, col_db); // (vektor & masik, meret_t kezdet, meret_t db)
    }
    //***********************************************************************
    // Nem ráfektetésbiztos.
    void math_add_nembiztos(const matrix & a, const matrix & b) { t.math_add(a.t, b.t); }
    //***********************************************************************
    // Nem ráfektetésbiztos.
    void math_sub_nembiztos(const matrix & a, const matrix & b) { t.math_sub(a.t, b.t); }
    //***********************************************************************
    // Nem ráfektetésbiztos.
    void math_neg_nembiztos() { t.math_neg(); }
    //***********************************************************************
    vektor<adattipus> & operator[](meret_t i) { return sorok[i]; }
    //***********************************************************************
    const vektor<adattipus> & operator[](meret_t i) const { return sorok[i]; }
    //***********************************************************************
    const adattipus & get_elem_unsafe(meret_t i) const { return t.unsafe(i); }
    //***********************************************************************
    adattipus & get_elem_unsafe(meret_t i) { return t.unsafe(i); }
    //***********************************************************************
    const adattipus & get_elem_unsafe(meret_t y, meret_t x) const { return sorok.unsafe(y).unsafe(x); }
    //***********************************************************************
    adattipus & get_elem_unsafe(meret_t row, meret_t col) { return sorok.unsafe(row).unsafe(col); }
    //***********************************************************************
    const adattipus & get_elem(meret_t row, meret_t col) const {
    //***********************************************************************
        if (is_szimm) {// így az m[i][j] valóban az i,j indexû elemet adja, de j>=i kötelezõ!
            return (col < row) ? sorok.unsafe(col).unsafe(row) : sorok.unsafe(row).unsafe(col);
        }
        else {
            return sorok.unsafe(row).unsafe(col);
        }
    }

    //***********************************************************************
    void debug_write(::std::ofstream & fs) const{
    //***********************************************************************
        for (meret_t i = 0; i < row; i++){
            sorok[i].debug_write(fs);
        }            
    }
    
    //***********************************************************************
    void print() const{
    //***********************************************************************
        for (meret_t i = 0; i < row; i++){
            sorok[i].print(is_szimm ? i : 0);
            ::std::cout << ::std::endl;
        }            
        ::std::cout << ::std::endl;
    }
    
    //***********************************************************************
    void print_z() const{
    //***********************************************************************
        for (meret_t i = 0; i < row; i++){
            sorok[i].print_z(is_szimm ? i : 0);
            ::std::cout << ::std::endl;
        }            
        ::std::cout << ::std::endl;
    }
    
    //***********************************************************************
    void print_size(char be = ' ', char end = '\n') const {
    //***********************************************************************
        printf("%c(%u,%u:%u)%c", be, row, col, t.size(), end);
        //::std::cout << be << '(' << row << ',' << col << ':' << t.size() << ')' << end;
    }
    
    //***********************************************************************
    void math_mul_t_nembiztos(const matrix & a, const matrix & b_t) {
    // Nem ráfektetésbiztos.
    //***********************************************************************
        egyforma_e_hiba(a.row, row,     "math_mul_t row row");
        egyforma_e_hiba(b_t.row, col,   "math_mul_t row col");
        egyforma_e_hiba(a.col, b_t.col, "math_mul_t col col");
        igaz_e_hiba(is_szimm || a.is_szimm || b_t.is_szimm, "matrix::math_mul_t", "symmetrical matrix not allowed");

        cuns ni = row, nj = col, nk = a.col;
        cuns di = row % 4, dj = col % 4, dk = a.col % 4;
        cuns hi = ni - di, hj = nj - dj, hk = nk - dk;
        adattipus *d0 = &t.unsafe(0), *d1 = d0 + nj, *d2 = d0 + 2 * nj, *d3 = d0 + 3 * nj;
        const adattipus *a0 = &a.t.unsafe(0), *a1 = a0 + nk, *a2 = a0 + 2 * nk, *a3 = a0 + 3 * nk;
        for (uns i = 0; i < hi; i += 4) {
            const adattipus *b0 = &b_t.t.unsafe(0), *b1 = b0 + nk, *b2 = b0 + 2 * nk, *b3 = b0 + 3 * nk;
            for (uns j = 0; j < hj; j += 4) {
                adattipus d[16] = { adattipus() };
                for (uns k = 0; k < hk; k += 4, a0 += 4, a1 += 4, a2 += 4, a3 += 4, b0 += 4, b1 += 4, b2 += 4, b3 += 4) {
                    d[0]  += a0[0] * b0[0] + a0[1] * b0[1] + a0[2] * b0[2] + a0[3] * b0[3];
                    d[1]  += a0[0] * b1[0] + a0[1] * b1[1] + a0[2] * b1[2] + a0[3] * b1[3];
                    d[2]  += a0[0] * b2[0] + a0[1] * b2[1] + a0[2] * b2[2] + a0[3] * b2[3];
                    d[3]  += a0[0] * b3[0] + a0[1] * b3[1] + a0[2] * b3[2] + a0[3] * b3[3];

                    d[4]  += a1[0] * b0[0] + a1[1] * b0[1] + a1[2] * b0[2] + a1[3] * b0[3];
                    d[5]  += a1[0] * b1[0] + a1[1] * b1[1] + a1[2] * b1[2] + a1[3] * b1[3];
                    d[6]  += a1[0] * b2[0] + a1[1] * b2[1] + a1[2] * b2[2] + a1[3] * b2[3];
                    d[7]  += a1[0] * b3[0] + a1[1] * b3[1] + a1[2] * b3[2] + a1[3] * b3[3];

                    d[8]  += a2[0] * b0[0] + a2[1] * b0[1] + a2[2] * b0[2] + a2[3] * b0[3];
                    d[9]  += a2[0] * b1[0] + a2[1] * b1[1] + a2[2] * b1[2] + a2[3] * b1[3];
                    d[10] += a2[0] * b2[0] + a2[1] * b2[1] + a2[2] * b2[2] + a2[3] * b2[3];
                    d[11] += a2[0] * b3[0] + a2[1] * b3[1] + a2[2] * b3[2] + a2[3] * b3[3];

                    d[12] += a3[0] * b0[0] + a3[1] * b0[1] + a3[2] * b0[2] + a3[3] * b0[3];
                    d[13] += a3[0] * b1[0] + a3[1] * b1[1] + a3[2] * b1[2] + a3[3] * b1[3];
                    d[14] += a3[0] * b2[0] + a3[1] * b2[1] + a3[2] * b2[2] + a3[3] * b2[3];
                    d[15] += a3[0] * b3[0] + a3[1] * b3[1] + a3[2] * b3[2] + a3[3] * b3[3];
                }
                for (uns k = 0; k < dk; k++) {
                    d[0]  += a0[k] * b0[k];	d[1]  += a0[k] * b1[k];	d[2]  += a0[k] * b2[k];	d[3]  += a0[k] * b3[k];
                    d[4]  += a1[k] * b0[k];	d[5]  += a1[k] * b1[k];	d[6]  += a1[k] * b2[k];	d[7]  += a1[k] * b3[k];
                    d[8]  += a2[k] * b0[k];	d[9]  += a2[k] * b1[k];	d[10] += a2[k] * b2[k];	d[11] += a2[k] * b3[k];
                    d[12] += a3[k] * b0[k];	d[13] += a3[k] * b1[k];	d[14] += a3[k] * b2[k];	d[15] += a3[k] * b3[k];
                }
                d0[0] = d[0];  d0[1] = d[1];  d0[2] = d[2];  d0[3] = d[3];
                d1[0] = d[4];  d1[1] = d[5];  d1[2] = d[6];  d1[3] = d[7];
                d2[0] = d[8];  d2[1] = d[9];  d2[2] = d[10]; d2[3] = d[11];
                d3[0] = d[12]; d3[1] = d[13]; d3[2] = d[14]; d3[3] = d[15];
                d0 += 4, d1 += 4, d2 += 4, d3 += 4;
                a0 -= hk, a1 -= hk, a2 -= hk, a3 -= hk, b0 += dk + 3 * nk, b1 += dk + 3 * nk, b2 += dk + 3 * nk, b3 += dk + 3 * nk;
            }
            for (uns j = 0; j < dj; j++) {
                d0[j] = d1[j] = d2[j] = d3[j] = 0;
                for (uns k = 0; k < nk; k++) {
                    d0[j] += a0[k] * b0[j*nk + k]; d1[j] += a1[k] * b0[j*nk + k]; d2[j] += a2[k] * b0[j*nk + k]; d3[j] += a3[k] * b0[j*nk + k];
                }
            }
            d0 += dj + 3 * nj; d1 += dj + 3 * nj; d2 += dj + 3 * nj; d3 += dj + 3 * nj;
            a0 += 4 * nk; a1 += 4 * nk; a2 += 4 * nk; a3 += 4 * nk;
        }
        const adattipus *b0 = &b_t.t.unsafe(0);
        for (uns i = 0; i < di; i++)
            for (uns j = 0; j < nj; j++) {
                d0[i*nj + j] = 0;
                for (uns k = 0; k < nk; k++)
                    d0[i*nj + j] += a0[i*nk + k] * b0[j*nk + k];
            }
    }

    //***********************************************************************
    void math_mul_t_biztos(const matrix & a, const matrix & b_t) {
    // kb 4%-kal lassabb, mint a pointeres verzió, de elvileg ráfektetett mátrixon is megy
    //***********************************************************************
        egyforma_e_hiba(a.row, row,     "math_mul_t row row");
        egyforma_e_hiba(b_t.row, col,   "math_mul_t row col");
        egyforma_e_hiba(a.col, b_t.col, "math_mul_t col col");
        igaz_e_hiba(is_szimm || a.is_szimm || b_t.is_szimm, "matrix::math_mul_t", "symmetrical matrix not allowed");

        cuns ni = row, nj = col, nk = a.col;
        cuns di = row % 4, dj = col % 4, dk = a.col % 4;
        cuns hi = ni - di, hj = nj - dj, hk = nk - dk;
        for (uns i = 0; i < hi; i += 4) {
            vektor<adattipus> & dest_sor_i0 = sorok.unsafe(i + 0);
            vektor<adattipus> & dest_sor_i1 = sorok.unsafe(i + 1);
            vektor<adattipus> & dest_sor_i2 = sorok.unsafe(i + 2);
            vektor<adattipus> & dest_sor_i3 = sorok.unsafe(i + 3);
            const vektor<adattipus> & a_sor_i0 = a.sorok.unsafe(i + 0);
            const vektor<adattipus> & a_sor_i1 = a.sorok.unsafe(i + 1);
            const vektor<adattipus> & a_sor_i2 = a.sorok.unsafe(i + 2);
            const vektor<adattipus> & a_sor_i3 = a.sorok.unsafe(i + 3);
            for (uns j = 0; j < hj; j += 4) {
                const vektor<adattipus> & b_sor_j0 = b_t.sorok.unsafe(j + 0);
                const vektor<adattipus> & b_sor_j1 = b_t.sorok.unsafe(j + 1);
                const vektor<adattipus> & b_sor_j2 = b_t.sorok.unsafe(j + 2);
                const vektor<adattipus> & b_sor_j3 = b_t.sorok.unsafe(j + 3);
                adattipus d[16] = { adattipus() };
                for (uns k = 0; k < hk; k += 4) {
                    d[0]  += a_sor_i0.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j0.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j0.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[1]  += a_sor_i0.unsafe(k + 0) * b_sor_j1.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j1.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j1.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[2]  += a_sor_i0.unsafe(k + 0) * b_sor_j2.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j2.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j2.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[3]  += a_sor_i0.unsafe(k + 0) * b_sor_j3.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j3.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j3.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);

                    d[4]  += a_sor_i1.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i1.unsafe(k + 1) * b_sor_j0.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j0.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[5]  += a_sor_i1.unsafe(k + 0) * b_sor_j1.unsafe(k + 0)
                           + a_sor_i1.unsafe(k + 1) * b_sor_j1.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j1.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[6]  += a_sor_i1.unsafe(k + 0) * b_sor_j2.unsafe(k + 0)
                           + a_sor_i1.unsafe(k + 1) * b_sor_j2.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j2.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[7]  += a_sor_i1.unsafe(k + 0) * b_sor_j3.unsafe(k + 0)
                           + a_sor_i1.unsafe(k + 1) * b_sor_j3.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j3.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);

                    d[8]  += a_sor_i2.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i2.unsafe(k + 1) * b_sor_j0.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j0.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[9]  += a_sor_i2.unsafe(k + 0) * b_sor_j1.unsafe(k + 0)
                           + a_sor_i2.unsafe(k + 1) * b_sor_j1.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j1.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[10] += a_sor_i2.unsafe(k + 0) * b_sor_j2.unsafe(k + 0)
                           + a_sor_i2.unsafe(k + 1) * b_sor_j2.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j2.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[11] += a_sor_i2.unsafe(k + 0) * b_sor_j3.unsafe(k + 0)
                           + a_sor_i2.unsafe(k + 1) * b_sor_j3.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j3.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);

                    d[12] += a_sor_i3.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i3.unsafe(k + 1) * b_sor_j0.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j0.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[13] += a_sor_i3.unsafe(k + 0) * b_sor_j1.unsafe(k + 0)
                           + a_sor_i3.unsafe(k + 1) * b_sor_j1.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j1.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[14] += a_sor_i3.unsafe(k + 0) * b_sor_j2.unsafe(k + 0)
                           + a_sor_i3.unsafe(k + 1) * b_sor_j2.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j2.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[15] += a_sor_i3.unsafe(k + 0) * b_sor_j3.unsafe(k + 0)
                           + a_sor_i3.unsafe(k + 1) * b_sor_j3.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j3.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);
                }
                for (uns k = hk; k < nk; k++) {
                    d[0]  += a_sor_i0.unsafe(k) * b_sor_j0.unsafe(k);
                    d[1]  += a_sor_i0.unsafe(k) * b_sor_j1.unsafe(k);
                    d[2]  += a_sor_i0.unsafe(k) * b_sor_j2.unsafe(k);
                    d[3]  += a_sor_i0.unsafe(k) * b_sor_j3.unsafe(k);

                    d[4]  += a_sor_i1.unsafe(k) * b_sor_j0.unsafe(k);
                    d[5]  += a_sor_i1.unsafe(k) * b_sor_j1.unsafe(k);
                    d[6]  += a_sor_i1.unsafe(k) * b_sor_j2.unsafe(k);
                    d[7]  += a_sor_i1.unsafe(k) * b_sor_j3.unsafe(k);

                    d[8]  += a_sor_i2.unsafe(k) * b_sor_j0.unsafe(k);
                    d[9]  += a_sor_i2.unsafe(k) * b_sor_j1.unsafe(k);
                    d[10] += a_sor_i2.unsafe(k) * b_sor_j2.unsafe(k);
                    d[11] += a_sor_i2.unsafe(k) * b_sor_j3.unsafe(k);

                    d[12] += a_sor_i3.unsafe(k) * b_sor_j0.unsafe(k);
                    d[13] += a_sor_i3.unsafe(k) * b_sor_j1.unsafe(k);
                    d[14] += a_sor_i3.unsafe(k) * b_sor_j2.unsafe(k);
                    d[15] += a_sor_i3.unsafe(k) * b_sor_j3.unsafe(k);
                }
                dest_sor_i0.unsafe(j + 0) = d[0];
                dest_sor_i0.unsafe(j + 1) = d[1];
                dest_sor_i0.unsafe(j + 2) = d[2];
                dest_sor_i0.unsafe(j + 3) = d[3];

                dest_sor_i1.unsafe(j + 0) = d[4];
                dest_sor_i1.unsafe(j + 1) = d[5];
                dest_sor_i1.unsafe(j + 2) = d[6];
                dest_sor_i1.unsafe(j + 3) = d[7];

                dest_sor_i2.unsafe(j + 0) = d[8];
                dest_sor_i2.unsafe(j + 1) = d[9];
                dest_sor_i2.unsafe(j + 2) = d[10];
                dest_sor_i2.unsafe(j + 3) = d[11];

                dest_sor_i3.unsafe(j + 0) = d[12];
                dest_sor_i3.unsafe(j + 1) = d[13];
                dest_sor_i3.unsafe(j + 2) = d[14];
                dest_sor_i3.unsafe(j + 3) = d[15];
            }
            for (uns j = hj; j < nj; j++) {
                const vektor<adattipus> & b_sor_j0 = b_t.sorok.unsafe(j);
                adattipus d[4] = { adattipus() };
                for (uns k = 0; k < nk; k++) {
                    d[0] += a_sor_i0.unsafe(k) * b_sor_j0.unsafe(k);
                    d[1] += a_sor_i1.unsafe(k) * b_sor_j0.unsafe(k);
                    d[2] += a_sor_i2.unsafe(k) * b_sor_j0.unsafe(k);
                    d[3] += a_sor_i3.unsafe(k) * b_sor_j0.unsafe(k);
                }
                dest_sor_i0.unsafe(j) = d[0];
                dest_sor_i1.unsafe(j) = d[1];
                dest_sor_i2.unsafe(j) = d[2];
                dest_sor_i3.unsafe(j) = d[3];
            }
        }
        for (uns i = hi; i < ni; i++) {
            vektor<adattipus> & dest_sor_i = sorok.unsafe(i);
            const vektor<adattipus> & a_sor_i = a.sorok.unsafe(i);
            for (uns j = 0; j < nj; j++) {
                const vektor<adattipus> & b_sor_j = b_t.sorok.unsafe(j);
                adattipus d = adattipus();
                for (uns k = 0; k < nk; k++)
                    d += a_sor_i.unsafe(k) * b_sor_j.unsafe(k);
                dest_sor_i.unsafe(j) = d;
            }
        }
    }

    //***********************************************************************
    void math_nmul_t_biztos(const matrix & a, const matrix & b_t) {
    // ráfektetett mátrixon is megy, a szorzat -1-szeresét adja
    //***********************************************************************
        egyforma_e_hiba(a.row, row,     "math_nmul_t row row");
        egyforma_e_hiba(b_t.row, col,   "math_nmul_t row col");
        egyforma_e_hiba(a.col, b_t.col, "math_nmul_t col col");
        igaz_e_hiba(is_szimm || a.is_szimm || b_t.is_szimm, "matrix::math_nmul_t", "symmetrical matrix not allowed");

        cuns ni = row, nj = col, nk = a.col;
        cuns di = row % 4, dj = col % 4, dk = a.col % 4;
        cuns hi = ni - di, hj = nj - dj, hk = nk - dk;
        for (uns i = 0; i < hi; i += 4) {
            vektor<adattipus> & dest_sor_i0 = sorok.unsafe(i + 0);
            vektor<adattipus> & dest_sor_i1 = sorok.unsafe(i + 1);
            vektor<adattipus> & dest_sor_i2 = sorok.unsafe(i + 2);
            vektor<adattipus> & dest_sor_i3 = sorok.unsafe(i + 3);
            const vektor<adattipus> & a_sor_i0 = a.sorok.unsafe(i + 0);
            const vektor<adattipus> & a_sor_i1 = a.sorok.unsafe(i + 1);
            const vektor<adattipus> & a_sor_i2 = a.sorok.unsafe(i + 2);
            const vektor<adattipus> & a_sor_i3 = a.sorok.unsafe(i + 3);
            for (uns j = 0; j < hj; j += 4) {
                const vektor<adattipus> & b_sor_j0 = b_t.sorok.unsafe(j + 0);
                const vektor<adattipus> & b_sor_j1 = b_t.sorok.unsafe(j + 1);
                const vektor<adattipus> & b_sor_j2 = b_t.sorok.unsafe(j + 2);
                const vektor<adattipus> & b_sor_j3 = b_t.sorok.unsafe(j + 3);
                adattipus d[16] = { adattipus() };
                for (uns k = 0; k < hk; k += 4) {
                    d[0]  += a_sor_i0.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j0.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j0.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[1]  += a_sor_i0.unsafe(k + 0) * b_sor_j1.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j1.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j1.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[2]  += a_sor_i0.unsafe(k + 0) * b_sor_j2.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j2.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j2.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[3]  += a_sor_i0.unsafe(k + 0) * b_sor_j3.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j3.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j3.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);

                    d[4]  += a_sor_i1.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i1.unsafe(k + 1) * b_sor_j0.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j0.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[5]  += a_sor_i1.unsafe(k + 0) * b_sor_j1.unsafe(k + 0)
                           + a_sor_i1.unsafe(k + 1) * b_sor_j1.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j1.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[6]  += a_sor_i1.unsafe(k + 0) * b_sor_j2.unsafe(k + 0)
                           + a_sor_i1.unsafe(k + 1) * b_sor_j2.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j2.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[7]  += a_sor_i1.unsafe(k + 0) * b_sor_j3.unsafe(k + 0)
                           + a_sor_i1.unsafe(k + 1) * b_sor_j3.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j3.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);

                    d[8]  += a_sor_i2.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i2.unsafe(k + 1) * b_sor_j0.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j0.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[9]  += a_sor_i2.unsafe(k + 0) * b_sor_j1.unsafe(k + 0)
                           + a_sor_i2.unsafe(k + 1) * b_sor_j1.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j1.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[10] += a_sor_i2.unsafe(k + 0) * b_sor_j2.unsafe(k + 0)
                           + a_sor_i2.unsafe(k + 1) * b_sor_j2.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j2.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[11] += a_sor_i2.unsafe(k + 0) * b_sor_j3.unsafe(k + 0)
                           + a_sor_i2.unsafe(k + 1) * b_sor_j3.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j3.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);

                    d[12] += a_sor_i3.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i3.unsafe(k + 1) * b_sor_j0.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j0.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[13] += a_sor_i3.unsafe(k + 0) * b_sor_j1.unsafe(k + 0)
                           + a_sor_i3.unsafe(k + 1) * b_sor_j1.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j1.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[14] += a_sor_i3.unsafe(k + 0) * b_sor_j2.unsafe(k + 0)
                           + a_sor_i3.unsafe(k + 1) * b_sor_j2.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j2.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[15] += a_sor_i3.unsafe(k + 0) * b_sor_j3.unsafe(k + 0)
                           + a_sor_i3.unsafe(k + 1) * b_sor_j3.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j3.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);
                }
                for (uns k = hk; k < nk; k++) {
                    d[0]  += a_sor_i0.unsafe(k) * b_sor_j0.unsafe(k);
                    d[1]  += a_sor_i0.unsafe(k) * b_sor_j1.unsafe(k);
                    d[2]  += a_sor_i0.unsafe(k) * b_sor_j2.unsafe(k);
                    d[3]  += a_sor_i0.unsafe(k) * b_sor_j3.unsafe(k);

                    d[4]  += a_sor_i1.unsafe(k) * b_sor_j0.unsafe(k);
                    d[5]  += a_sor_i1.unsafe(k) * b_sor_j1.unsafe(k);
                    d[6]  += a_sor_i1.unsafe(k) * b_sor_j2.unsafe(k);
                    d[7]  += a_sor_i1.unsafe(k) * b_sor_j3.unsafe(k);

                    d[8]  += a_sor_i2.unsafe(k) * b_sor_j0.unsafe(k);
                    d[9]  += a_sor_i2.unsafe(k) * b_sor_j1.unsafe(k);
                    d[10] += a_sor_i2.unsafe(k) * b_sor_j2.unsafe(k);
                    d[11] += a_sor_i2.unsafe(k) * b_sor_j3.unsafe(k);

                    d[12] += a_sor_i3.unsafe(k) * b_sor_j0.unsafe(k);
                    d[13] += a_sor_i3.unsafe(k) * b_sor_j1.unsafe(k);
                    d[14] += a_sor_i3.unsafe(k) * b_sor_j2.unsafe(k);
                    d[15] += a_sor_i3.unsafe(k) * b_sor_j3.unsafe(k);
                }
                dest_sor_i0.unsafe(j + 0) = -d[0];
                dest_sor_i0.unsafe(j + 1) = -d[1];
                dest_sor_i0.unsafe(j + 2) = -d[2];
                dest_sor_i0.unsafe(j + 3) = -d[3];

                dest_sor_i1.unsafe(j + 0) = -d[4];
                dest_sor_i1.unsafe(j + 1) = -d[5];
                dest_sor_i1.unsafe(j + 2) = -d[6];
                dest_sor_i1.unsafe(j + 3) = -d[7];

                dest_sor_i2.unsafe(j + 0) = -d[8];
                dest_sor_i2.unsafe(j + 1) = -d[9];
                dest_sor_i2.unsafe(j + 2) = -d[10];
                dest_sor_i2.unsafe(j + 3) = -d[11];

                dest_sor_i3.unsafe(j + 0) = -d[12];
                dest_sor_i3.unsafe(j + 1) = -d[13];
                dest_sor_i3.unsafe(j + 2) = -d[14];
                dest_sor_i3.unsafe(j + 3) = -d[15];
            }
            for (uns j = hj; j < nj; j++) {
                const vektor<adattipus> & b_sor_j0 = b_t.sorok.unsafe(j);
                adattipus d[4] = { adattipus() };
                for (uns k = 0; k < nk; k++) {
                    d[0] += a_sor_i0.unsafe(k) * b_sor_j0.unsafe(k);
                    d[1] += a_sor_i1.unsafe(k) * b_sor_j0.unsafe(k);
                    d[2] += a_sor_i2.unsafe(k) * b_sor_j0.unsafe(k);
                    d[3] += a_sor_i3.unsafe(k) * b_sor_j0.unsafe(k);
                }
                dest_sor_i0.unsafe(j) = -d[0];
                dest_sor_i1.unsafe(j) = -d[1];
                dest_sor_i2.unsafe(j) = -d[2];
                dest_sor_i3.unsafe(j) = -d[3];
            }
        }
        for (uns i = hi; i < ni; i++) {
            vektor<adattipus> & dest_sor_i = sorok.unsafe(i);
            const vektor<adattipus> & a_sor_i = a.sorok.unsafe(i);
            for (uns j = 0; j < nj; j++) {
                const vektor<adattipus> & b_sor_j = b_t.sorok.unsafe(j);
                adattipus d = adattipus();
                for (uns k = 0; k < nk; k++)
                    d += a_sor_i.unsafe(k) * b_sor_j.unsafe(k);
                dest_sor_i.unsafe(j) = -d;
            }
        }
    }

    //***********************************************************************
    void math_add_mul_t_nembiztos(const matrix & c, const matrix & a, const matrix & b_t) {
    // Nem ráfektetésbiztos.
    //***********************************************************************
        egyforma_e_hiba(c.col, col,     "math_add_mul_t col");
        egyforma_e_hiba(c.row, row,     "math_add_mul_t row");
        egyforma_e_hiba(a.row, row,     "math_add_mul_t row row");
        egyforma_e_hiba(b_t.row, col,   "math_add_mul_t row col");
        egyforma_e_hiba(a.col, b_t.col, "math_add_mul_t col col");
        igaz_e_hiba(is_szimm || a.is_szimm || b_t.is_szimm, "matrix::math_add_mul_t", "symmetrical matrix not allowed");

        cuns ni = row, nj = col, nk = a.col;
        cuns di = row % 4, dj = col % 4, dk = a.col % 4;
        cuns hi = ni - di, hj = nj - dj, hk = nk - dk;
        adattipus *d0 = &t.unsafe(0), *d1 = d0 + nj, *d2 = d0 + 2 * nj, *d3 = d0 + 3 * nj;
        const adattipus *c0 = &c.t.unsafe(0), *c1 = c0 + nj, *c2 = c0 + 2 * nj, *c3 = c0 + 3 * nj;
        const adattipus *a0 = &a.t.unsafe(0), *a1 = a0 + nk, *a2 = a0 + 2 * nk, *a3 = a0 + 3 * nk;
        for (uns i = 0; i < hi; i += 4) {
            const adattipus *b0 = &b_t.t.unsafe(0), *b1 = b0 + nk, *b2 = b0 + 2 * nk, *b3 = b0 + 3 * nk;
            for (uns j = 0; j < hj; j += 4) {
                adattipus d[16] = { adattipus() };
                for (uns k = 0; k < hk; k += 4, a0 += 4, a1 += 4, a2 += 4, a3 += 4, b0 += 4, b1 += 4, b2 += 4, b3 += 4) {
                    d[0] += a0[0] * b0[0] + a0[1] * b0[1] + a0[2] * b0[2] + a0[3] * b0[3];
                    d[1] += a0[0] * b1[0] + a0[1] * b1[1] + a0[2] * b1[2] + a0[3] * b1[3];
                    d[2] += a0[0] * b2[0] + a0[1] * b2[1] + a0[2] * b2[2] + a0[3] * b2[3];
                    d[3] += a0[0] * b3[0] + a0[1] * b3[1] + a0[2] * b3[2] + a0[3] * b3[3];

                    d[4] += a1[0] * b0[0] + a1[1] * b0[1] + a1[2] * b0[2] + a1[3] * b0[3];
                    d[5] += a1[0] * b1[0] + a1[1] * b1[1] + a1[2] * b1[2] + a1[3] * b1[3];
                    d[6] += a1[0] * b2[0] + a1[1] * b2[1] + a1[2] * b2[2] + a1[3] * b2[3];
                    d[7] += a1[0] * b3[0] + a1[1] * b3[1] + a1[2] * b3[2] + a1[3] * b3[3];

                    d[8] += a2[0] * b0[0] + a2[1] * b0[1] + a2[2] * b0[2] + a2[3] * b0[3];
                    d[9] += a2[0] * b1[0] + a2[1] * b1[1] + a2[2] * b1[2] + a2[3] * b1[3];
                    d[10] += a2[0] * b2[0] + a2[1] * b2[1] + a2[2] * b2[2] + a2[3] * b2[3];
                    d[11] += a2[0] * b3[0] + a2[1] * b3[1] + a2[2] * b3[2] + a2[3] * b3[3];

                    d[12] += a3[0] * b0[0] + a3[1] * b0[1] + a3[2] * b0[2] + a3[3] * b0[3];
                    d[13] += a3[0] * b1[0] + a3[1] * b1[1] + a3[2] * b1[2] + a3[3] * b1[3];
                    d[14] += a3[0] * b2[0] + a3[1] * b2[1] + a3[2] * b2[2] + a3[3] * b2[3];
                    d[15] += a3[0] * b3[0] + a3[1] * b3[1] + a3[2] * b3[2] + a3[3] * b3[3];
                }
                for (uns k = 0; k < dk; k++) {
                    d[0] += a0[k] * b0[k];	d[1] += a0[k] * b1[k];	d[2] += a0[k] * b2[k];	d[3] += a0[k] * b3[k];
                    d[4] += a1[k] * b0[k];	d[5] += a1[k] * b1[k];	d[6] += a1[k] * b2[k];	d[7] += a1[k] * b3[k];
                    d[8] += a2[k] * b0[k];	d[9] += a2[k] * b1[k];	d[10] += a2[k] * b2[k];	d[11] += a2[k] * b3[k];
                    d[12] += a3[k] * b0[k];	d[13] += a3[k] * b1[k];	d[14] += a3[k] * b2[k];	d[15] += a3[k] * b3[k];
                }
                d0[0] = d[0]  + c0[0]; d0[1] = d[1]  + c0[1]; d0[2] = d[2]  + c0[2]; d0[3] = d[3]  + c0[3];
                d1[0] = d[4]  + c1[0]; d1[1] = d[5]  + c1[1]; d1[2] = d[6]  + c1[2]; d1[3] = d[7]  + c1[3];
                d2[0] = d[8]  + c2[0]; d2[1] = d[9]  + c2[1]; d2[2] = d[10] + c2[2]; d2[3] = d[11] + c2[3];
                d3[0] = d[12] + c3[0]; d3[1] = d[13] + c3[1]; d3[2] = d[14] + c3[2]; d3[3] = d[15] + c3[3];
                d0 += 4, d1 += 4, d2 += 4, d3 += 4;
                c0 += 4, c1 += 4, c2 += 4, c3 += 4;
                a0 -= hk, a1 -= hk, a2 -= hk, a3 -= hk, b0 += dk + 3 * nk, b1 += dk + 3 * nk, b2 += dk + 3 * nk, b3 += dk + 3 * nk;
            }
            for (uns j = 0; j < dj; j++) {
                d0[j] = c0[j]; d1[j] = c1[j]; d2[j] = c2[j]; d3[j] = c3[j];
                for (uns k = 0; k < nk; k++) {
                    d0[j] += a0[k] * b0[j*nk + k]; d1[j] += a1[k] * b0[j*nk + k]; d2[j] += a2[k] * b0[j*nk + k]; d3[j] += a3[k] * b0[j*nk + k];
                }
            }
            d0 += dj + 3 * nj; d1 += dj + 3 * nj; d2 += dj + 3 * nj; d3 += dj + 3 * nj;
            c0 += dj + 3 * nj; c1 += dj + 3 * nj; c2 += dj + 3 * nj; c3 += dj + 3 * nj;
            a0 += 4 * nk; a1 += 4 * nk; a2 += 4 * nk; a3 += 4 * nk;
        }
        const adattipus *b0 = &b_t.t.unsafe(0);
        for (uns i = 0; i < di; i++)
            for (uns j = 0; j < nj; j++) {
                d0[i*nj + j] = c0[i*nj + j];
                for (uns k = 0; k < nk; k++)
                    d0[i*nj + j] += a0[i*nk + k] * b0[j*nk + k];
            }
    }

    //***********************************************************************
    void math_add_mul_t_biztos(const matrix & c, const matrix & a, const matrix & b_t) {
    // kb 4%-kal lassabb, mint a pointeres verzió, de elvileg ráfektetett mátrixon is megy
    //***********************************************************************
        egyforma_e_hiba(a.row, row,     "math_add_mul_t_ row row");
        egyforma_e_hiba(b_t.row, col,   "math_add_mul_t_ row col");
        egyforma_e_hiba(a.col, b_t.col, "math_add_mul_t_ col col");
        igaz_e_hiba(is_szimm || a.is_szimm || b_t.is_szimm, "matrix::math_add_mul_t_", "symmetrical matrix not allowed");

        cuns ni = row, nj = col, nk = a.col;
        cuns di = row % 4, dj = col % 4, dk = a.col % 4;
        cuns hi = ni - di, hj = nj - dj, hk = nk - dk;
        for (uns i = 0; i < hi; i += 4) {
            vektor<adattipus> & dest_sor_i0 = sorok.unsafe(i + 0);
            vektor<adattipus> & dest_sor_i1 = sorok.unsafe(i + 1);
            vektor<adattipus> & dest_sor_i2 = sorok.unsafe(i + 2);
            vektor<adattipus> & dest_sor_i3 = sorok.unsafe(i + 3);
            const vektor<adattipus> & a_sor_i0 = a.sorok.unsafe(i + 0);
            const vektor<adattipus> & a_sor_i1 = a.sorok.unsafe(i + 1);
            const vektor<adattipus> & a_sor_i2 = a.sorok.unsafe(i + 2);
            const vektor<adattipus> & a_sor_i3 = a.sorok.unsafe(i + 3);
            const vektor<adattipus> & c_sor_i0 = c.sorok.unsafe(i + 0);
            const vektor<adattipus> & c_sor_i1 = c.sorok.unsafe(i + 1);
            const vektor<adattipus> & c_sor_i2 = c.sorok.unsafe(i + 2);
            const vektor<adattipus> & c_sor_i3 = c.sorok.unsafe(i + 3);
            for (uns j = 0; j < hj; j += 4) {
                const vektor<adattipus> & b_sor_j0 = b_t.sorok.unsafe(j + 0);
                const vektor<adattipus> & b_sor_j1 = b_t.sorok.unsafe(j + 1);
                const vektor<adattipus> & b_sor_j2 = b_t.sorok.unsafe(j + 2);
                const vektor<adattipus> & b_sor_j3 = b_t.sorok.unsafe(j + 3);
                adattipus d[16] = { adattipus() };
                for (uns k = 0; k < hk; k += 4) {
                    d[0]  += a_sor_i0.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j0.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j0.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[1]  += a_sor_i0.unsafe(k + 0) * b_sor_j1.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j1.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j1.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[2]  += a_sor_i0.unsafe(k + 0) * b_sor_j2.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j2.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j2.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[3]  += a_sor_i0.unsafe(k + 0) * b_sor_j3.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j3.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j3.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);

                    d[4]  += a_sor_i1.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i1.unsafe(k + 1) * b_sor_j0.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j0.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[5]  += a_sor_i1.unsafe(k + 0) * b_sor_j1.unsafe(k + 0)
                           + a_sor_i1.unsafe(k + 1) * b_sor_j1.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j1.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[6]  += a_sor_i1.unsafe(k + 0) * b_sor_j2.unsafe(k + 0)
                           + a_sor_i1.unsafe(k + 1) * b_sor_j2.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j2.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[7]  += a_sor_i1.unsafe(k + 0) * b_sor_j3.unsafe(k + 0)
                           + a_sor_i1.unsafe(k + 1) * b_sor_j3.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j3.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);

                    d[8]  += a_sor_i2.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i2.unsafe(k + 1) * b_sor_j0.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j0.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[9]  += a_sor_i2.unsafe(k + 0) * b_sor_j1.unsafe(k + 0)
                           + a_sor_i2.unsafe(k + 1) * b_sor_j1.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j1.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[10] += a_sor_i2.unsafe(k + 0) * b_sor_j2.unsafe(k + 0)
                           + a_sor_i2.unsafe(k + 1) * b_sor_j2.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j2.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[11] += a_sor_i2.unsafe(k + 0) * b_sor_j3.unsafe(k + 0)
                           + a_sor_i2.unsafe(k + 1) * b_sor_j3.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j3.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);

                    d[12] += a_sor_i3.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i3.unsafe(k + 1) * b_sor_j0.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j0.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[13] += a_sor_i3.unsafe(k + 0) * b_sor_j1.unsafe(k + 0)
                           + a_sor_i3.unsafe(k + 1) * b_sor_j1.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j1.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[14] += a_sor_i3.unsafe(k + 0) * b_sor_j2.unsafe(k + 0)
                           + a_sor_i3.unsafe(k + 1) * b_sor_j2.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j2.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[15] += a_sor_i3.unsafe(k + 0) * b_sor_j3.unsafe(k + 0)
                           + a_sor_i3.unsafe(k + 1) * b_sor_j3.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j3.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);
                }
                for (uns k = hk; k < nk; k++) {
                    d[0]  += a_sor_i0.unsafe(k) * b_sor_j0.unsafe(k);
                    d[1]  += a_sor_i0.unsafe(k) * b_sor_j1.unsafe(k);
                    d[2]  += a_sor_i0.unsafe(k) * b_sor_j2.unsafe(k);
                    d[3]  += a_sor_i0.unsafe(k) * b_sor_j3.unsafe(k);

                    d[4]  += a_sor_i1.unsafe(k) * b_sor_j0.unsafe(k);
                    d[5]  += a_sor_i1.unsafe(k) * b_sor_j1.unsafe(k);
                    d[6]  += a_sor_i1.unsafe(k) * b_sor_j2.unsafe(k);
                    d[7]  += a_sor_i1.unsafe(k) * b_sor_j3.unsafe(k);

                    d[8]  += a_sor_i2.unsafe(k) * b_sor_j0.unsafe(k);
                    d[9]  += a_sor_i2.unsafe(k) * b_sor_j1.unsafe(k);
                    d[10] += a_sor_i2.unsafe(k) * b_sor_j2.unsafe(k);
                    d[11] += a_sor_i2.unsafe(k) * b_sor_j3.unsafe(k);

                    d[12] += a_sor_i3.unsafe(k) * b_sor_j0.unsafe(k);
                    d[13] += a_sor_i3.unsafe(k) * b_sor_j1.unsafe(k);
                    d[14] += a_sor_i3.unsafe(k) * b_sor_j2.unsafe(k);
                    d[15] += a_sor_i3.unsafe(k) * b_sor_j3.unsafe(k);
                }
                dest_sor_i0.unsafe(j + 0) = d[0] + c_sor_i0.unsafe(j + 0);
                dest_sor_i0.unsafe(j + 1) = d[1] + c_sor_i0.unsafe(j + 1);
                dest_sor_i0.unsafe(j + 2) = d[2] + c_sor_i0.unsafe(j + 2);
                dest_sor_i0.unsafe(j + 3) = d[3] + c_sor_i0.unsafe(j + 3);

                dest_sor_i1.unsafe(j + 0) = d[4] + c_sor_i1.unsafe(j + 0);
                dest_sor_i1.unsafe(j + 1) = d[5] + c_sor_i1.unsafe(j + 1);
                dest_sor_i1.unsafe(j + 2) = d[6] + c_sor_i1.unsafe(j + 2);
                dest_sor_i1.unsafe(j + 3) = d[7] + c_sor_i1.unsafe(j + 3);

                dest_sor_i2.unsafe(j + 0) = d[8] + c_sor_i2.unsafe(j + 0);
                dest_sor_i2.unsafe(j + 1) = d[9] + c_sor_i2.unsafe(j + 1);
                dest_sor_i2.unsafe(j + 2) = d[10] + c_sor_i2.unsafe(j + 2);
                dest_sor_i2.unsafe(j + 3) = d[11] + c_sor_i2.unsafe(j + 3);

                dest_sor_i3.unsafe(j + 0) = d[12] + c_sor_i3.unsafe(j + 0);
                dest_sor_i3.unsafe(j + 1) = d[13] + c_sor_i3.unsafe(j + 1);
                dest_sor_i3.unsafe(j + 2) = d[14] + c_sor_i3.unsafe(j + 2);
                dest_sor_i3.unsafe(j + 3) = d[15] + c_sor_i3.unsafe(j + 3);
            }
            for (uns j = hj; j < nj; j++) {
                const vektor<adattipus> & b_sor_j0 = b_t.sorok.unsafe(j);
                adattipus d[4] = { adattipus() };
                for (uns k = 0; k < nk; k++) {
                    d[0] += a_sor_i0.unsafe(k) * b_sor_j0.unsafe(k);
                    d[1] += a_sor_i1.unsafe(k) * b_sor_j0.unsafe(k);
                    d[2] += a_sor_i2.unsafe(k) * b_sor_j0.unsafe(k);
                    d[3] += a_sor_i3.unsafe(k) * b_sor_j0.unsafe(k);
                }
                dest_sor_i0.unsafe(j) = d[0] + c_sor_i0.unsafe(j);
                dest_sor_i1.unsafe(j) = d[1] + c_sor_i1.unsafe(j);
                dest_sor_i2.unsafe(j) = d[2] + c_sor_i2.unsafe(j);
                dest_sor_i3.unsafe(j) = d[3] + c_sor_i3.unsafe(j);
            }
        }
        for (uns i = hi; i < ni; i++) {
            vektor<adattipus> & dest_sor_i = sorok.unsafe(i);
            const vektor<adattipus> & a_sor_i = a.sorok.unsafe(i);
            const vektor<adattipus> & c_sor_i = c.sorok.unsafe(i);
            for (uns j = 0; j < nj; j++) {
                const vektor<adattipus> & b_sor_j = b_t.sorok.unsafe(j);
                adattipus d = adattipus();
                for (uns k = 0; k < nk; k++)
                    d += a_sor_i.unsafe(k) * b_sor_j.unsafe(k);
                dest_sor_i.unsafe(j) = d + c_sor_i.unsafe(j);
            }
        }
    }

    //***********************************************************************
    void math_sub_mul_t_biztos(const matrix & c, const matrix & a, const matrix & b_t) {
    // ráfektetett mátrixon is megy, *this = c - a * b_t 
    //***********************************************************************
        egyforma_e_hiba(a.row, row,     "math_sub_mul_t row row");
        egyforma_e_hiba(b_t.row, col,   "math_sub_mul_t row col");
        egyforma_e_hiba(a.col, b_t.col, "math_sub_mul_t col col");
        igaz_e_hiba(is_szimm || a.is_szimm || b_t.is_szimm, "matrix::math_sub_mul_t", "symmetrical matrix not allowed");

        cuns ni = row, nj = col, nk = a.col;
        cuns di = row % 4, dj = col % 4, dk = a.col % 4;
        cuns hi = ni - di, hj = nj - dj, hk = nk - dk;
        for (uns i = 0; i < hi; i += 4) {
            vektor<adattipus> & dest_sor_i0 = sorok.unsafe(i + 0);
            vektor<adattipus> & dest_sor_i1 = sorok.unsafe(i + 1);
            vektor<adattipus> & dest_sor_i2 = sorok.unsafe(i + 2);
            vektor<adattipus> & dest_sor_i3 = sorok.unsafe(i + 3);
            const vektor<adattipus> & a_sor_i0 = a.sorok.unsafe(i + 0);
            const vektor<adattipus> & a_sor_i1 = a.sorok.unsafe(i + 1);
            const vektor<adattipus> & a_sor_i2 = a.sorok.unsafe(i + 2);
            const vektor<adattipus> & a_sor_i3 = a.sorok.unsafe(i + 3);
            const vektor<adattipus> & c_sor_i0 = c.sorok.unsafe(i + 0);
            const vektor<adattipus> & c_sor_i1 = c.sorok.unsafe(i + 1);
            const vektor<adattipus> & c_sor_i2 = c.sorok.unsafe(i + 2);
            const vektor<adattipus> & c_sor_i3 = c.sorok.unsafe(i + 3);
            for (uns j = 0; j < hj; j += 4) {
                const vektor<adattipus> & b_sor_j0 = b_t.sorok.unsafe(j + 0);
                const vektor<adattipus> & b_sor_j1 = b_t.sorok.unsafe(j + 1);
                const vektor<adattipus> & b_sor_j2 = b_t.sorok.unsafe(j + 2);
                const vektor<adattipus> & b_sor_j3 = b_t.sorok.unsafe(j + 3);
                adattipus d[16] = { adattipus() };
                for (uns k = 0; k < hk; k += 4) {
                    d[0]  += a_sor_i0.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j0.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j0.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[1]  += a_sor_i0.unsafe(k + 0) * b_sor_j1.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j1.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j1.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[2]  += a_sor_i0.unsafe(k + 0) * b_sor_j2.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j2.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j2.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[3]  += a_sor_i0.unsafe(k + 0) * b_sor_j3.unsafe(k + 0) 
                           + a_sor_i0.unsafe(k + 1) * b_sor_j3.unsafe(k + 1) 
                           + a_sor_i0.unsafe(k + 2) * b_sor_j3.unsafe(k + 2) 
                           + a_sor_i0.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);

                    d[4]  += a_sor_i1.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i1.unsafe(k + 1) * b_sor_j0.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j0.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[5]  += a_sor_i1.unsafe(k + 0) * b_sor_j1.unsafe(k + 0)
                           + a_sor_i1.unsafe(k + 1) * b_sor_j1.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j1.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[6]  += a_sor_i1.unsafe(k + 0) * b_sor_j2.unsafe(k + 0)
                           + a_sor_i1.unsafe(k + 1) * b_sor_j2.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j2.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[7]  += a_sor_i1.unsafe(k + 0) * b_sor_j3.unsafe(k + 0)
                           + a_sor_i1.unsafe(k + 1) * b_sor_j3.unsafe(k + 1)
                           + a_sor_i1.unsafe(k + 2) * b_sor_j3.unsafe(k + 2)
                           + a_sor_i1.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);

                    d[8]  += a_sor_i2.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i2.unsafe(k + 1) * b_sor_j0.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j0.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[9]  += a_sor_i2.unsafe(k + 0) * b_sor_j1.unsafe(k + 0)
                           + a_sor_i2.unsafe(k + 1) * b_sor_j1.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j1.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[10] += a_sor_i2.unsafe(k + 0) * b_sor_j2.unsafe(k + 0)
                           + a_sor_i2.unsafe(k + 1) * b_sor_j2.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j2.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[11] += a_sor_i2.unsafe(k + 0) * b_sor_j3.unsafe(k + 0)
                           + a_sor_i2.unsafe(k + 1) * b_sor_j3.unsafe(k + 1)
                           + a_sor_i2.unsafe(k + 2) * b_sor_j3.unsafe(k + 2)
                           + a_sor_i2.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);

                    d[12] += a_sor_i3.unsafe(k + 0) * b_sor_j0.unsafe(k + 0) 
                           + a_sor_i3.unsafe(k + 1) * b_sor_j0.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j0.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j0.unsafe(k + 3);
                    d[13] += a_sor_i3.unsafe(k + 0) * b_sor_j1.unsafe(k + 0)
                           + a_sor_i3.unsafe(k + 1) * b_sor_j1.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j1.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j1.unsafe(k + 3);
                    d[14] += a_sor_i3.unsafe(k + 0) * b_sor_j2.unsafe(k + 0)
                           + a_sor_i3.unsafe(k + 1) * b_sor_j2.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j2.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j2.unsafe(k + 3);
                    d[15] += a_sor_i3.unsafe(k + 0) * b_sor_j3.unsafe(k + 0)
                           + a_sor_i3.unsafe(k + 1) * b_sor_j3.unsafe(k + 1)
                           + a_sor_i3.unsafe(k + 2) * b_sor_j3.unsafe(k + 2)
                           + a_sor_i3.unsafe(k + 3) * b_sor_j3.unsafe(k + 3);
                }
                for (uns k = hk; k < nk; k++) {
                    d[0]  += a_sor_i0.unsafe(k) * b_sor_j0.unsafe(k);
                    d[1]  += a_sor_i0.unsafe(k) * b_sor_j1.unsafe(k);
                    d[2]  += a_sor_i0.unsafe(k) * b_sor_j2.unsafe(k);
                    d[3]  += a_sor_i0.unsafe(k) * b_sor_j3.unsafe(k);

                    d[4]  += a_sor_i1.unsafe(k) * b_sor_j0.unsafe(k);
                    d[5]  += a_sor_i1.unsafe(k) * b_sor_j1.unsafe(k);
                    d[6]  += a_sor_i1.unsafe(k) * b_sor_j2.unsafe(k);
                    d[7]  += a_sor_i1.unsafe(k) * b_sor_j3.unsafe(k);

                    d[8]  += a_sor_i2.unsafe(k) * b_sor_j0.unsafe(k);
                    d[9]  += a_sor_i2.unsafe(k) * b_sor_j1.unsafe(k);
                    d[10] += a_sor_i2.unsafe(k) * b_sor_j2.unsafe(k);
                    d[11] += a_sor_i2.unsafe(k) * b_sor_j3.unsafe(k);

                    d[12] += a_sor_i3.unsafe(k) * b_sor_j0.unsafe(k);
                    d[13] += a_sor_i3.unsafe(k) * b_sor_j1.unsafe(k);
                    d[14] += a_sor_i3.unsafe(k) * b_sor_j2.unsafe(k);
                    d[15] += a_sor_i3.unsafe(k) * b_sor_j3.unsafe(k);
                }
                dest_sor_i0.unsafe(j + 0) = c_sor_i0.unsafe(j + 0) - d[0];
                dest_sor_i0.unsafe(j + 1) = c_sor_i0.unsafe(j + 1) - d[1];
                dest_sor_i0.unsafe(j + 2) = c_sor_i0.unsafe(j + 2) - d[2];
                dest_sor_i0.unsafe(j + 3) = c_sor_i0.unsafe(j + 3) - d[3];

                dest_sor_i1.unsafe(j + 0) = c_sor_i1.unsafe(j + 0) - d[4];
                dest_sor_i1.unsafe(j + 1) = c_sor_i1.unsafe(j + 1) - d[5];
                dest_sor_i1.unsafe(j + 2) = c_sor_i1.unsafe(j + 2) - d[6];
                dest_sor_i1.unsafe(j + 3) = c_sor_i1.unsafe(j + 3) - d[7];

                dest_sor_i2.unsafe(j + 0) = c_sor_i2.unsafe(j + 0) - d[8];
                dest_sor_i2.unsafe(j + 1) = c_sor_i2.unsafe(j + 1) - d[9];
                dest_sor_i2.unsafe(j + 2) = c_sor_i2.unsafe(j + 2) - d[10];
                dest_sor_i2.unsafe(j + 3) = c_sor_i2.unsafe(j + 3) - d[11];

                dest_sor_i3.unsafe(j + 0) = c_sor_i3.unsafe(j + 0) - d[12];
                dest_sor_i3.unsafe(j + 1) = c_sor_i3.unsafe(j + 1) - d[13];
                dest_sor_i3.unsafe(j + 2) = c_sor_i3.unsafe(j + 2) - d[14];
                dest_sor_i3.unsafe(j + 3) = c_sor_i3.unsafe(j + 3) - d[15];
            }
            for (uns j = hj; j < nj; j++) {
                const vektor<adattipus> & b_sor_j0 = b_t.sorok.unsafe(j);
                adattipus d[4] = { adattipus() };
                for (uns k = 0; k < nk; k++) {
                    d[0] += a_sor_i0.unsafe(k) * b_sor_j0.unsafe(k);
                    d[1] += a_sor_i1.unsafe(k) * b_sor_j0.unsafe(k);
                    d[2] += a_sor_i2.unsafe(k) * b_sor_j0.unsafe(k);
                    d[3] += a_sor_i3.unsafe(k) * b_sor_j0.unsafe(k);
                }
                dest_sor_i0.unsafe(j) = c_sor_i0.unsafe(j) - d[0];
                dest_sor_i1.unsafe(j) = c_sor_i1.unsafe(j) - d[1];
                dest_sor_i2.unsafe(j) = c_sor_i2.unsafe(j) - d[2];
                dest_sor_i3.unsafe(j) = c_sor_i3.unsafe(j) - d[3];
            }
        }
        for (uns i = hi; i < ni; i++) {
            vektor<adattipus> & dest_sor_i = sorok.unsafe(i);
            const vektor<adattipus> & a_sor_i = a.sorok.unsafe(i);
            const vektor<adattipus> & c_sor_i = c.sorok.unsafe(i);
            for (uns j = 0; j < nj; j++) {
                const vektor<adattipus> & b_sor_j = b_t.sorok.unsafe(j);
                adattipus d = adattipus();
                for (uns k = 0; k < nk; k++)
                    d += a_sor_i.unsafe(k) * b_sor_j.unsafe(k);
                dest_sor_i.unsafe(j) = c_sor_i.unsafe(j) - d;
            }
        }
    }

    //***********************************************************************
    void transp(const matrix & a) {
    //***********************************************************************
        egyforma_e_hiba(a.col, row, "matrix::transp col-row");
        egyforma_e_hiba(a.row, col, "matrix::transp row-col");
        igaz_e_hiba(is_szimm, "matrix::transp", "symmetrical matrix not allowed");
        for (meret_t i = 0; i < row; i++)
            for (meret_t j = 0; j < col; j++)
                sorok.unsafe(i).unsafe(j) = a.sorok.unsafe(j).unsafe(i);
    }

    
    //***********************************************************************
    void copy(const matrix & a) {
    // ráfektetettre is mûködik
    //***********************************************************************
        egyforma_e_hiba(a.col, col, "matrix::copy col-col");
        egyforma_e_hiba(a.row, row, "matrix::copy row-row");
        igaz_e_hiba(is_szimm, "matrix::copy", "symmetrical matrix not allowed");
        for (meret_t i = 0; i < row; i++)
            for (meret_t j = 0; j < col; j++)
                sorok.unsafe(i).unsafe(j) = a.sorok.unsafe(i).unsafe(j);
    }

    
    //***********************************************************************
    void submatrix_copy(const matrix & src, uns start_dest_row, uns start_dest_col, uns start_src_row, uns start_src_col, uns db_row, uns db_col){
    //***********************************************************************
        if (db_col == 0 || db_row == 0)
            return;
        kisebb_e_hiba(start_dest_row + db_row - 1, row,     "matrix::submatrix_copy start_dest_row + db_row");
        kisebb_e_hiba(start_src_row  + db_row - 1, src.row, "matrix::submatrix_copy start_src_row  + db_row");
        if (is_szimm) {
            if (start_dest_row >= start_dest_col + db_col) // azaz start_dest_row > start_dest_col + db_col - 1
                ; // ha a cél teljes egészében az alsó háromszögbe esik, akkor eldobjuk
            else {
                if (src.is_szimm) {  // szimm => szimm
                    if (start_dest_row + db_row - 1 <= start_dest_col) { // a cél teljesen fõátló fölötti
                        if (start_src_row + db_row - 1 <= start_src_col) { // a forrás teljesen fõátló fölötti
                            for (uns i = 0; i < db_row; i++)
                                sorok.unsafe(start_dest_row + i).subvektor_copy(src.sorok.unsafe(start_src_row + i), start_dest_col, start_src_col, db_col);
                        }
                        else if (start_src_row >= start_src_col + db_col) { // a forrás teljesen fõátló alatti (elvileg ilyen nem fordul elõ)
                            uns src_row = start_src_row, dest_row = start_dest_row;
                            for (uns i = 0; i < db_row; i++, src_row++, dest_row++) {
                                uns src_col = start_src_col, dest_col = start_dest_col;
                                for (uns j = 0; j < db_col; j++, src_col++, dest_col++) {
                                    sorok.unsafe(dest_row).unsafe(dest_col) = src.sorok.unsafe(src_col).unsafe(src_row); // dest(i,j) = src(j,i)
                                }
                            }
                        }
                        else { // a forrást metszi a fõátló (elvileg ilyen nem fordul elõ)
                            uns src_row = start_src_row, dest_row = start_dest_row;
                            for (uns i = 0; i < db_row; i++, src_row++, dest_row++) {
                                uns dest_col = start_dest_col;
                                cuns ig = start_src_col + db_col;
                                for (uns src_col = start_src_col; src_col < src_row && src_col < ig; src_col++, dest_col++) {
                                    sorok.unsafe(dest_row).unsafe(dest_col) = src.sorok.unsafe(src_col).unsafe(src_row); // dest(i,j) = src(j,i)
                                }
                                for (uns src_col = start_src_col > src_row ? start_src_col : src_row; src_col < ig; src_col++, dest_col++) {
                                    sorok.unsafe(dest_row).unsafe(dest_col) = src.sorok.unsafe(src_row).unsafe(src_col); // dest(i,j) = src(i,j)
                                }
                            }
                        }
                    }
                    else { // a célt metszi a fõátló
                        if (start_src_row + db_row - 1 <= start_src_col) { // a forrás teljesen fõátló fölötti (elvileg ilyen nem fordul elõ)
                            uns src_row = start_src_row;
                            for (uns dest_row = start_dest_row; dest_row < start_dest_row + db_row; dest_row++, src_row++) {
                                uns dest_col = dest_row > start_dest_col ? dest_row : start_dest_col;
                                uns src_col = start_src_col + dest_col - start_dest_col;
                                for (; dest_col < start_dest_col + db_col; dest_col++, src_col++) {
                                    sorok.unsafe(dest_row).unsafe(dest_col) = src.sorok.unsafe(src_row).unsafe(src_col); // dest(i,j) = src(i,j)
                                }
                            }
                        }
                        else if (start_src_row >= start_src_col + db_col) { // a forrás teljesen fõátló alatti (elvileg ilyen nem fordul elõ)
                            uns src_row = start_src_row;
                            for (uns dest_row = start_dest_row; dest_row < start_dest_row + db_row; dest_row++, src_row++) {
                                uns dest_col = dest_row > start_dest_col ? dest_row : start_dest_col;
                                uns src_col = start_src_col + dest_col - start_dest_col;
                                for (; dest_col < start_dest_col + db_col; dest_col++, src_col++) {
                                    sorok.unsafe(dest_row).unsafe(dest_col) = src.sorok.unsafe(src_col).unsafe(src_row); // dest(i,j) = src(j,i)
                                }
                            }
                        }
                        else { // a forrást metszi a fõátló
                            uns src_row = start_src_row;
                            for (uns dest_row = start_dest_row; dest_row < start_dest_row + db_row; dest_row++, src_row++) {
                                uns dest_col = dest_row > start_dest_col ? dest_row : start_dest_col;
                                uns src_col = start_src_col + dest_col - start_dest_col;
                                for (; dest_col < start_dest_col + db_col; dest_col++, src_col++) {
                                    if (src_col >= src_row)
                                        sorok.unsafe(dest_row).unsafe(dest_col) = src.sorok.unsafe(src_row).unsafe(src_col); // dest(i,j) = src(i,j)
                                    else
                                        sorok.unsafe(dest_row).unsafe(dest_col) = src.sorok.unsafe(src_col).unsafe(src_row); // dest(i,j) = src(j,i)
                                }
                            }
                        }
                    }
                }
                else { // nemszimm => szimm (lehet ilyen?)
                    if (start_dest_row + db_row - 1 <= start_dest_col) { // a cél teljesen fõátló fölötti
                        for (uns i = 0; i < db_row; i++)
                            sorok.unsafe(start_dest_row + i).subvektor_copy(src.sorok.unsafe(start_src_row + i), start_dest_col, start_src_col, db_col);
                    }
                    else { // fõátló alatti rész is van a célban
                        uns src_row = start_src_row, dest_row = start_dest_row;
                        for (uns i = 0; i < db_row; i++, src_row++, dest_row++) {
                            uns dest_col = start_dest_col >= dest_row ? start_dest_col : dest_row;
                            uns src_col = start_src_col + dest_col - start_dest_col;
                            for (; dest_col < start_dest_col + db_col; src_col++, dest_col++) {
                                sorok.unsafe(dest_row).unsafe(dest_col) = src.sorok.unsafe(src_row).unsafe(src_col); // dest(i,j) = src(i,j)
                            }
                        }
                    }
                }
            }
        }
        else {
            if (src.is_szimm) { // szimm => nemszimm
                if (start_src_row + db_row - 1 <= start_src_col) { // a forrás teljesen fõátló fölötti
                    for (uns i = 0; i < db_row; i++)
                        sorok.unsafe(start_dest_row + i).subvektor_copy(src.sorok.unsafe(start_src_row + i), start_dest_col, start_src_col, db_col);
                }
                else if (start_src_row >= start_src_col + db_col) { // a forrás teljesen fõátló alatti
                    uns src_row = start_src_row, dest_row = start_dest_row;
                    for (uns i = 0; i < db_row; i++, src_row++, dest_row++) {
                        uns src_col = start_src_col, dest_col = start_dest_col;
                        for (uns j = 0; j < db_col; j++, src_col++, dest_col++) {
                            sorok.unsafe(dest_row).unsafe(dest_col) = src.sorok.unsafe(src_col).unsafe(src_row); // dest(i,j) = src(j,i)
                        }
                    }
                }
                else { // a forrást metszi a fõátló
                    uns src_row = start_src_row, dest_row = start_dest_row;
                    for (uns i = 0; i < db_row; i++, src_row++, dest_row++) {
                        uns src_col = start_src_col, dest_col = start_dest_col;
                        for (uns j = 0; j < db_col; j++, src_col++, dest_col++) {
                            if (src_col >= src_row)
                                sorok.unsafe(dest_row).unsafe(dest_col) = src.sorok.unsafe(src_row).unsafe(src_col); // dest(i,j) = src(i,j)
                            else
                                sorok.unsafe(dest_row).unsafe(dest_col) = src.sorok.unsafe(src_col).unsafe(src_row); // dest(i,j) = src(j,i)
                        }
                    }
                }
            }
            else { // nemszimm => nemszimm
                for (uns i = 0; i < db_row; i++)
                    sorok.unsafe(start_dest_row + i).subvektor_copy(src.sorok.unsafe(start_src_row + i), start_dest_col, start_src_col, db_col);
            }
        }
    }


    //***********************************************************************
    void submatrix_puszegyenlo(const matrix & src, uns start_dest_row, uns start_dest_col, uns start_src_row, uns start_src_col, uns db_row, uns db_col){
    //***********************************************************************
        if (db_col == 0 || db_row == 0)
            return;
        kisebb_e_hiba(start_dest_row + db_row - 1, row,     "matrix::submatrix_puszegyenlo start_dest_row + db_row");
        kisebb_e_hiba(start_src_row  + db_row - 1, src.row, "matrix::submatrix_puszegyenlo start_src_row  + db_row");
        if (is_szimm) {
            if (start_dest_row >= start_dest_col + db_col) // azaz start_dest_row > start_dest_col + db_col - 1
                ; // ha a cél teljes egészében az alsó háromszögbe esik, akkor eldobjuk
            else {
                if (src.is_szimm) {  // szimm => szimm
                    if (start_dest_row + db_row - 1 <= start_dest_col) { // a cél teljesen fõátló fölötti
                        if (start_src_row + db_row - 1 <= start_src_col) { // a forrás teljesen fõátló fölötti
                            for (uns i = 0; i < db_row; i++)
                                sorok.unsafe(start_dest_row + i).subvektor_pluszegyenlo(src.sorok.unsafe(start_src_row + i), start_dest_col, start_src_col, db_col);
                        }
                        else if (start_src_row >= start_src_col + db_col) { // a forrás teljesen fõátló alatti (elvileg ilyen nem fordul elõ)
                            uns src_row = start_src_row, dest_row = start_dest_row;
                            for (uns i = 0; i < db_row; i++, src_row++, dest_row++) {
                                uns src_col = start_src_col, dest_col = start_dest_col;
                                for (uns j = 0; j < db_col; j++, src_col++, dest_col++) {
                                    sorok.unsafe(dest_row).unsafe(dest_col) += src.sorok.unsafe(src_col).unsafe(src_row); // dest(i,j) = src(j,i)
                                }
                            }
                        }
                        else { // a forrást metszi a fõátló (elvileg ilyen nem fordul elõ)
                            uns src_row = start_src_row, dest_row = start_dest_row;
                            for (uns i = 0; i < db_row; i++, src_row++, dest_row++) {
                                uns dest_col = start_dest_col;
                                cuns ig = start_src_col + db_col;
                                for (uns src_col = start_src_col; src_col < src_row && src_col < ig; src_col++, dest_col++) {
                                    sorok.unsafe(dest_row).unsafe(dest_col) += src.sorok.unsafe(src_col).unsafe(src_row); // dest(i,j) = src(j,i)
                                }
                                for (uns src_col = start_src_col > src_row ? start_src_col : src_row; src_col < ig; src_col++, dest_col++) {
                                    sorok.unsafe(dest_row).unsafe(dest_col) += src.sorok.unsafe(src_row).unsafe(src_col); // dest(i,j) = src(i,j)
                                }
                            }
                        }
                    }
                    else { // a célt metszi a fõátló
                        if (start_src_row + db_row - 1 <= start_src_col) { // a forrás teljesen fõátló fölötti (elvileg ilyen nem fordul elõ)
                            uns src_row = start_src_row;
                            for (uns dest_row = start_dest_row; dest_row < start_dest_row + db_row; dest_row++, src_row++) {
                                uns dest_col = dest_row > start_dest_col ? dest_row : start_dest_col;
                                uns src_col = start_src_col + dest_col - start_dest_col;
                                for (; dest_col < start_dest_col + db_col; dest_col++, src_col++) {
                                    sorok.unsafe(dest_row).unsafe(dest_col) += src.sorok.unsafe(src_row).unsafe(src_col); // dest(i,j) = src(i,j)
                                }
                            }
                        }
                        else if (start_src_row >= start_src_col + db_col) { // a forrás teljesen fõátló alatti (elvileg ilyen nem fordul elõ)
                            uns src_row = start_src_row;
                            for (uns dest_row = start_dest_row; dest_row < start_dest_row + db_row; dest_row++, src_row++) {
                                uns dest_col = dest_row > start_dest_col ? dest_row : start_dest_col;
                                uns src_col = start_src_col + dest_col - start_dest_col;
                                for (; dest_col < start_dest_col + db_col; dest_col++, src_col++) {
                                    sorok.unsafe(dest_row).unsafe(dest_col) += src.sorok.unsafe(src_col).unsafe(src_row); // dest(i,j) = src(j,i)
                                }
                            }
                        }
                        else { // a forrást metszi a fõátló
                            uns src_row = start_src_row;
                            for (uns dest_row = start_dest_row; dest_row < start_dest_row + db_row; dest_row++, src_row++) {
                                uns dest_col = dest_row > start_dest_col ? dest_row : start_dest_col;
                                uns src_col = start_src_col + dest_col - start_dest_col;
                                for (; dest_col < start_dest_col + db_col; dest_col++, src_col++) {
                                    if (src_col >= src_row)
                                        sorok.unsafe(dest_row).unsafe(dest_col) += src.sorok.unsafe(src_row).unsafe(src_col); // dest(i,j) = src(i,j)
                                    else
                                        sorok.unsafe(dest_row).unsafe(dest_col) += src.sorok.unsafe(src_col).unsafe(src_row); // dest(i,j) = src(j,i)
                                }
                            }
                        }
                    }
                }
                else { // nemszimm => szimm (lehet ilyen?)
                    if (start_dest_row + db_row - 1 <= start_dest_col) { // a cél teljesen fõátló fölötti
                        for (uns i = 0; i < db_row; i++)
                            sorok.unsafe(start_dest_row + i).subvektor_pluszegyenlo(src.sorok.unsafe(start_src_row + i), start_dest_col, start_src_col, db_col);
                    }
                    else { // fõátló alatti rész is van a célban
                        uns src_row = start_src_row, dest_row = start_dest_row;
                        for (uns i = 0; i < db_row; i++, src_row++, dest_row++) {
                            uns dest_col = start_dest_col >= dest_row ? start_dest_col : dest_row;
                            uns src_col = start_src_col + dest_col - start_dest_col;
                            for (; dest_col < start_dest_col + db_col; src_col++, dest_col++) {
                                sorok.unsafe(dest_row).unsafe(dest_col) += src.sorok.unsafe(src_row).unsafe(src_col); // dest(i,j) = src(i,j)
                            }
                        }
                    }
                }
            }
        }
        else {
            if (src.is_szimm) { // szimm => nemszimm
                if (start_src_row + db_row - 1 <= start_src_col) { // a forrás teljesen fõátló fölötti
                    for (uns i = 0; i < db_row; i++)
                        sorok.unsafe(start_dest_row + i).subvektor_pluszegyenlo(src.sorok.unsafe(start_src_row + i), start_dest_col, start_src_col, db_col);
                }
                else if (start_src_row >= start_src_col + db_col) { // a forrás teljesen fõátló alatti
                    uns src_row = start_src_row, dest_row = start_dest_row;
                    for (uns i = 0; i < db_row; i++, src_row++, dest_row++) {
                        uns src_col = start_src_col, dest_col = start_dest_col;
                        for (uns j = 0; j < db_col; j++, src_col++, dest_col++) {
                            sorok.unsafe(dest_row).unsafe(dest_col) += src.sorok.unsafe(src_col).unsafe(src_row); // dest(i,j) = src(j,i)
                        }
                    }
                }
                else { // a forrást metszi a fõátló
                    uns src_row = start_src_row, dest_row = start_dest_row;
                    for (uns i = 0; i < db_row; i++, src_row++, dest_row++) {
                        uns src_col = start_src_col, dest_col = start_dest_col;
                        for (uns j = 0; j < db_col; j++, src_col++, dest_col++) {
                            if (src_col >= src_row)
                                sorok.unsafe(dest_row).unsafe(dest_col) += src.sorok.unsafe(src_row).unsafe(src_col); // dest(i,j) = src(i,j)
                            else
                                sorok.unsafe(dest_row).unsafe(dest_col) += src.sorok.unsafe(src_col).unsafe(src_row); // dest(i,j) = src(j,i)
                        }
                    }
                }
            }
            else { // nemszimm => nemszimm
                for (uns i = 0; i < db_row; i++)
                    sorok.unsafe(start_dest_row + i).subvektor_pluszegyenlo(src.sorok.unsafe(start_src_row + i), start_dest_col, start_src_col, db_col);
            }
        }
    }


    //***********************************************************************
    void submatrix_add(const matrix & src_1, const matrix & src_2, uns start_dest_row, uns start_dest_col, uns start_src_1_row, uns start_src_1_col, uns start_src_2_row, uns start_src_2_col, uns db_row, uns db_col){
    //***********************************************************************
        if (db_col == 0 || db_row == 0)
            return;
        kisebb_e_hiba(start_dest_row + db_row - 1, row,         "matrix::submatrix_add start_dest_row + db_row");
        kisebb_e_hiba(start_src_1_row  + db_row - 1, src_1.row, "matrix::submatrix_add start_src_1_row  + db_row");
        kisebb_e_hiba(start_src_2_row  + db_row - 1, src_2.row, "matrix::submatrix_add start_src_2_row  + db_row");

        // a leggyakoribb eset, ha semmi sem szimmetrikus

        if (!is_szimm && !src_1.is_szimm && !src_2.is_szimm) {
            for (uns i = 0; i < db_row; i++)
                sorok.unsafe(start_dest_row + i).subvektor_add(src_1.sorok.unsafe(start_src_1_row + i), src_2.sorok.unsafe(start_src_2_row + i), start_dest_col, start_src_1_col, start_src_2_col, db_col);
            return;
        }

        // Valami szimmetrikus: túl komplex lenne kibontani. Remélem, így is elég gyors lesz.

        submatrix_copy(src_1, start_dest_row, start_dest_col, start_src_1_row, start_src_1_col, db_row, db_col);
        submatrix_puszegyenlo(src_2, start_dest_row, start_dest_col, start_src_2_row, start_src_2_col, db_row, db_col);
    }


    //***********************************************************************
    friend inline void math_mul(vektor<adattipus> & cel, const matrix & src1, const vektor<adattipus> & src2) {
    //***********************************************************************
        egyforma_e_hiba(cel.size(), src1.get_row(), "vektor math_mul size");
        igaz_e_hiba(src1.is_szimm, "matrix::math_mul vektor", "symmetrical matrix not allowed");
        for (meret_t i = 0; i < cel.size(); i++)
            cel.unsafe(i) = math_mul(src1.sorok.unsafe(i), src2);
    }


    //***********************************************************************
    friend inline void math_add_mul(vektor<adattipus> & cel, const vektor<adattipus> & addando, const matrix & src1, const vektor<adattipus> & src2) {
    //***********************************************************************
        egyforma_e_hiba(cel.size(), addando.size(), "vektor math_add_mul size");
        egyforma_e_hiba(cel.size(), src1.get_row(), "vektor math_mul size");
        igaz_e_hiba(src1.is_szimm, "matrix::math_add_mul vektor", "symmetrical matrix not allowed");
        for (meret_t i = 0; i < cel.size(); i++)
            cel.unsafe(i) = addando.unsafe(i) + math_mul(src1.sorok.unsafe(i), src2);
    }


    //***********************************************************************
    friend inline void math_add_mul_symm(vektor<adattipus> & cel, const vektor<adattipus> & addando, const matrix & src1, const vektor<adattipus> & src2) {
    //***********************************************************************
        egyforma_e_hiba(cel.size(), addando.size(), "vektor math_add_mul size");
        egyforma_e_hiba(cel.size(), src1.get_row(), "vektor math_mul size");
        igaz_e_hiba(!src1.is_szimm, "matrix::math_add_mul_symm vektor", "nonsymmetrical matrix not allowed");
        for (meret_t row = 0; row < cel.size(); row++) {
            adattipus sum = addando.unsafe(row);
            for (meret_t col = 0; col < row; col++) {
                sum += src1.sorok.unsafe(col).unsafe(row) * src2.unsafe(col); // (j,i)
            }
            for (meret_t col = row; col < cel.size(); col++) {
                sum += src1.sorok.unsafe(row).unsafe(col) * src2.unsafe(col); // (i,j)
            }
            cel.unsafe(row) = sum;
        }
    }


    //****************************************************************
    void math_ninv_np_() {
    // fõelemkiválasztás nélkül
    // Ez nem ráfektetésbiztos
    //****************************************************************
        egyforma_e_hiba(row, col, "matrix::math_inv_np");
        igaz_e_hiba(is_szimm, "matrix::math_ninv_np_", "symmetrical matrix not allowed");

        if (row == 1) {
            t.unsafe(0) = adattipus(-1) / t.unsafe(0); 
            return;
        }
        if (row == 2) {
            adattipus * const p = &t.unsafe(0);
            adattipus p0 = abs(p[0]) < 1e-20f ? 1e20f : adattipus(1.0f) / p[0];
            adattipus p2 = -p[2] * p0;
            adattipus p1 = p[1] * p0;
            adattipus p3 = p[3] + p2 * p[1];
            adattipus oszto2 = abs(p3) < 1e-20f ? 1e20f : adattipus(1.0f) / p3;
            adattipus C2 = -p1 * oszto2;
            p[0] = -(p0 + C2 * p2);
            p[1] = -C2;
            p[2] = -p2 * oszto2;
            p[3] = -oszto2;
            return;
        }

        for (meret_t i = 0; i < row; i++){
            vektor<adattipus> & sor_i = sorok.unsafe(i);
            adattipus oszto = abs(sor_i.unsafe(i)) < 1e-20f ? 1e20f : adattipus(1.0f) / sor_i.unsafe(i);
            meret_t j;
            for (j = 0; j < i; j++) {
                vektor<adattipus> & sor_j = sorok.unsafe(j);
                adattipus C = sor_j.unsafe(i) * oszto;
                meret_t k;
                for (k = 0; k < i; k++) sor_j.unsafe(k) -= C * sor_i.unsafe(k);
                sor_j.unsafe(k) = C; // nem negált invnél = -C;
                for (k++; k < row; k++) sor_j.unsafe(k) -= C * sor_i.unsafe(k);
            }
            for (j++; j < row; j++) {
                vektor<adattipus> & sor_j = sorok.unsafe(j);
                adattipus C = sor_j.unsafe(i) * oszto;
                meret_t k;
                for (k = 0; k < i; k++) sor_j.unsafe(k) -= C * sor_i.unsafe(k);
                sor_j.unsafe(k) = C; // nem negált invnél = -C;
                for (k++; k < row; k++) sor_j.unsafe(k) -= C * sor_i.unsafe(k);
            }
            for (j = 0; j < i; j++) sor_i.unsafe(j) *= oszto;
            sor_i.unsafe(j) = -oszto; // nem negált invnél = oszto;
            for (j++; j < row; j++) sor_i.unsafe(j) *= oszto;
            //if(i==2)return;
        }
    }

    //****************************************************************
    void math_ninv_np_blokk_2x2() {
    // fõelemkiválasztás nélkül
    //****************************************************************
        egyforma_e_hiba(row, col, "matrix::math_ninv_np_blokk_2x2");
        igaz_e_hiba(is_szimm, "matrix::math_ninv_np_blokk_2x2", "symmetrical matrix not allowed");

        if (row == 1) {
            sorok.unsafe(0).unsafe(0) = adattipus(-1) / sorok.unsafe(0).unsafe(0);
            return;
        }
        if (row == 2) {
            math_ninv_2x2_fv(&sorok.unsafe(0).unsafe(0), &sorok.unsafe(1).unsafe(0));
            return;
        }

        if (row == 4) {
            math_ninv_4x4_fv(&sorok.unsafe(0).unsafe(0), &sorok.unsafe(1).unsafe(0), &sorok.unsafe(2).unsafe(0), &sorok.unsafe(3).unsafe(0));
            return;
        }

        cuns drow = row % 2;
        cuns hrow = row - drow;
        cuns nrow = row;
        for (meret_t i = 0; i < hrow; i += 2){
            vektor<adattipus> & sor_i0 = sorok.unsafe(i);
            vektor<adattipus> & sor_i1 = sorok.unsafe(i + 1);
            
            // (i,i)-nél lévõ elem inverze
            const adattipus p0 = abs(sor_i0.unsafe(i)) < 1e-20f ? 1e20f : adattipus(1.0f) / sor_i0.unsafe(i);
            const adattipus p1o = sor_i0.unsafe(i + 1);
            const adattipus p1 = p1o * p0;
            const adattipus p2 = -sor_i1.unsafe(i) * p0;
            const adattipus p3 = sor_i1.unsafe(i+1) + p2 * p1o;
            const adattipus oszto2 = abs(p3) < 1e-20f ? 1e20f : adattipus(1.0f) / p3;
            const adattipus C2 = -p1 * oszto2;
            const adattipus a0 = p0 + C2 * p2; // a0...a3 az inverz
            const adattipus a1 = C2;
            const adattipus a2 = p2 * oszto2;
            const adattipus a3 = oszto2;

            meret_t j;
            for (j = 0; j < i; j += 2) {
                vektor<adattipus> & sor_j0 = sorok.unsafe(j);
                vektor<adattipus> & sor_j1 = sorok.unsafe(j + 1);
                const adattipus C0 = sor_j0.unsafe(i) * a0 + sor_j0.unsafe(i + 1) * a2;
                const adattipus C1 = sor_j0.unsafe(i) * a1 + sor_j0.unsafe(i + 1) * a3;
                const adattipus C2 = sor_j1.unsafe(i) * a0 + sor_j1.unsafe(i + 1) * a2;
                const adattipus C3 = sor_j1.unsafe(i) * a1 + sor_j1.unsafe(i + 1) * a3;
                meret_t k;
                for (k = 0; k < i; k += 2) {
                    sor_j0.unsafe(k)     -= C0 * sor_i0.unsafe(k)     + C1 * sor_i1.unsafe(k);
                    sor_j0.unsafe(k + 1) -= C0 * sor_i0.unsafe(k + 1) + C1 * sor_i1.unsafe(k + 1);
                    sor_j1.unsafe(k)     -= C2 * sor_i0.unsafe(k)     + C3 * sor_i1.unsafe(k);
                    sor_j1.unsafe(k + 1) -= C2 * sor_i0.unsafe(k + 1) + C3 * sor_i1.unsafe(k + 1);
                }
                sor_j0.unsafe(k)     = C0; // nem negált invnél = -C;
                sor_j0.unsafe(k + 1) = C1;
                sor_j1.unsafe(k)     = C2;
                sor_j1.unsafe(k + 1) = C3;
                for (k += 2; k < hrow; k += 2) {
                    sor_j0.unsafe(k)     -= C0 * sor_i0.unsafe(k)     + C1 * sor_i1.unsafe(k);
                    sor_j0.unsafe(k + 1) -= C0 * sor_i0.unsafe(k + 1) + C1 * sor_i1.unsafe(k + 1);
                    sor_j1.unsafe(k)     -= C2 * sor_i0.unsafe(k)     + C3 * sor_i1.unsafe(k);
                    sor_j1.unsafe(k + 1) -= C2 * sor_i0.unsafe(k + 1) + C3 * sor_i1.unsafe(k + 1);
                }
                for (k = hrow; k < nrow; k++) {
                    sor_j0.unsafe(k)     -= C0 * sor_i0.unsafe(k)     + C1 * sor_i1.unsafe(k);
                    sor_j1.unsafe(k)     -= C2 * sor_i0.unsafe(k)     + C3 * sor_i1.unsafe(k);
                }
            }
            for (j+=2; j < hrow; j += 2) {
                vektor<adattipus> & sor_j0 = sorok.unsafe(j);
                vektor<adattipus> & sor_j1 = sorok.unsafe(j + 1);
                const adattipus C0 = sor_j0.unsafe(i) * a0 + sor_j0.unsafe(i + 1) * a2;
                const adattipus C1 = sor_j0.unsafe(i) * a1 + sor_j0.unsafe(i + 1) * a3;
                const adattipus C2 = sor_j1.unsafe(i) * a0 + sor_j1.unsafe(i + 1) * a2;
                const adattipus C3 = sor_j1.unsafe(i) * a1 + sor_j1.unsafe(i + 1) * a3;
                meret_t k;
                for (k = 0; k < i; k += 2) {
                    sor_j0.unsafe(k)     -= C0 * sor_i0.unsafe(k)     + C1 * sor_i1.unsafe(k);
                    sor_j0.unsafe(k + 1) -= C0 * sor_i0.unsafe(k + 1) + C1 * sor_i1.unsafe(k + 1);
                    sor_j1.unsafe(k)     -= C2 * sor_i0.unsafe(k)     + C3 * sor_i1.unsafe(k);
                    sor_j1.unsafe(k + 1) -= C2 * sor_i0.unsafe(k + 1) + C3 * sor_i1.unsafe(k + 1);
                }
                sor_j0.unsafe(k)     = C0; // nem negált invnél = -C;
                sor_j0.unsafe(k + 1) = C1;
                sor_j1.unsafe(k)     = C2;
                sor_j1.unsafe(k + 1) = C3;
                for (k += 2; k < hrow; k += 2) {
                    sor_j0.unsafe(k)     -= C0 * sor_i0.unsafe(k)     + C1 * sor_i1.unsafe(k);
                    sor_j0.unsafe(k + 1) -= C0 * sor_i0.unsafe(k + 1) + C1 * sor_i1.unsafe(k + 1);
                    sor_j1.unsafe(k)     -= C2 * sor_i0.unsafe(k)     + C3 * sor_i1.unsafe(k);
                    sor_j1.unsafe(k + 1) -= C2 * sor_i0.unsafe(k + 1) + C3 * sor_i1.unsafe(k + 1);
                }
                for (k = hrow; k < nrow; k++) {
                    sor_j0.unsafe(k)     -= C0 * sor_i0.unsafe(k)     + C1 * sor_i1.unsafe(k);
                    sor_j1.unsafe(k)     -= C2 * sor_i0.unsafe(k)     + C3 * sor_i1.unsafe(k);
                }
            }
            for (j = hrow; j < nrow; j++) {
                vektor<adattipus> & sor_j0 = sorok.unsafe(j);
                const adattipus C0 = sor_j0.unsafe(i) * a0 + sor_j0.unsafe(i + 1) * a2;
                const adattipus C1 = sor_j0.unsafe(i) * a1 + sor_j0.unsafe(i + 1) * a3;
                meret_t k;
                for (k = 0; k < i; k += 2) {
                    sor_j0.unsafe(k)     -= C0 * sor_i0.unsafe(k)     + C1 * sor_i1.unsafe(k);
                    sor_j0.unsafe(k + 1) -= C0 * sor_i0.unsafe(k + 1) + C1 * sor_i1.unsafe(k + 1);
                }
                sor_j0.unsafe(k)     = C0; // nem negált invnél = -C;
                sor_j0.unsafe(k + 1) = C1;
                for (k += 2; k < hrow; k += 2) {
                    sor_j0.unsafe(k)     -= C0 * sor_i0.unsafe(k)     + C1 * sor_i1.unsafe(k);
                    sor_j0.unsafe(k + 1) -= C0 * sor_i0.unsafe(k + 1) + C1 * sor_i1.unsafe(k + 1);
                }
                for (k = hrow; k < nrow; k++) {
                    sor_j0.unsafe(k)     -= C0 * sor_i0.unsafe(k)     + C1 * sor_i1.unsafe(k);
                }
            }
            for (j = 0; j < i; j += 2) {
                const adattipus b0 = a0 * sor_i0.unsafe(j)     + a1 * sor_i1.unsafe(j);
                const adattipus b1 = a0 * sor_i0.unsafe(j + 1) + a1 * sor_i1.unsafe(j + 1);
                const adattipus b2 = a2 * sor_i0.unsafe(j)     + a3 * sor_i1.unsafe(j);
                const adattipus b3 = a2 * sor_i0.unsafe(j + 1) + a3 * sor_i1.unsafe(j + 1);
                sor_i0.unsafe(j)     = b0;
                sor_i0.unsafe(j + 1) = b1;
                sor_i1.unsafe(j)     = b2;
                sor_i1.unsafe(j + 1) = b3;
            }
            sor_i0.unsafe(j)     = -a0; // nem negált invnél = ai;
            sor_i0.unsafe(j + 1) = -a1;
            sor_i1.unsafe(j)     = -a2;
            sor_i1.unsafe(j + 1) = -a3;
            for (j += 2; j < hrow; j += 2) {
                const adattipus b0 = a0 * sor_i0.unsafe(j)     + a1 * sor_i1.unsafe(j);
                const adattipus b1 = a0 * sor_i0.unsafe(j + 1) + a1 * sor_i1.unsafe(j + 1);
                const adattipus b2 = a2 * sor_i0.unsafe(j)     + a3 * sor_i1.unsafe(j);
                const adattipus b3 = a2 * sor_i0.unsafe(j + 1) + a3 * sor_i1.unsafe(j + 1);
                sor_i0.unsafe(j)     = b0;
                sor_i0.unsafe(j + 1) = b1;
                sor_i1.unsafe(j)     = b2;
                sor_i1.unsafe(j + 1) = b3;
            }
            for (j = hrow; j < nrow; j++) {
                const adattipus b0 = a0 * sor_i0.unsafe(j)     + a1 * sor_i1.unsafe(j);
                const adattipus b2 = a2 * sor_i0.unsafe(j)     + a3 * sor_i1.unsafe(j);
                sor_i0.unsafe(j)     = b0;
                sor_i1.unsafe(j)     = b2;
            }
        }
        for (meret_t i = hrow; i < nrow; i++) {
            vektor<adattipus> & sor_i = sorok.unsafe(i);
            adattipus oszto = abs(sor_i.unsafe(i)) < 1e-20f ? 1e20f : adattipus(1.0f) / sor_i.unsafe(i);
            meret_t j;
            for (j = 0; j < i; j++) {
                vektor<adattipus> & sor_j = sorok.unsafe(j);
                adattipus C = sor_j.unsafe(i) * oszto;
                meret_t k;
                for (k = 0; k < i; k++) sor_j.unsafe(k) -= C * sor_i.unsafe(k);
                sor_j.unsafe(k) = C; // nem negált invnél = -C;
                for (k++; k < row; k++) sor_j.unsafe(k) -= C * sor_i.unsafe(k);
            }
            for (j++; j < row; j++) {
                vektor<adattipus> & sor_j = sorok.unsafe(j);
                adattipus C = sor_j.unsafe(i) * oszto;
                meret_t k;
                for (k = 0; k < i; k++) sor_j.unsafe(k) -= C * sor_i.unsafe(k);
                sor_j.unsafe(k) = C; // nem negált invnél = -C;
                for (k++; k < row; k++) sor_j.unsafe(k) -= C * sor_i.unsafe(k);
            }
            for (j = 0; j < i; j++) sor_i.unsafe(j) *= oszto;
            sor_i.unsafe(j) = -oszto; // nem negált invnél = oszto;
            for (j++; j < row; j++) sor_i.unsafe(j) *= oszto;
        }
    }

    //****************************************************************
    static inline void math_ninv_2x2_fv(adattipus * be0, adattipus * be1) {
    // 1 db 4x4-es mátrixot invertál, bemenet a 4 sor elsõ elemének címe
    //****************************************************************
        adattipus p0 = abs(be0[0]) < 1e-20f ? 1e20f : adattipus(1.0f) / be0[0];
        adattipus p2 = -be1[0] * p0;
        adattipus p1 = be0[1] * p0;
        adattipus p3 = be1[1] + p2 * be0[1];
        adattipus oszto2 = abs(p3) < 1e-20f ? 1e20f : adattipus(1.0f) / p3;
        adattipus C2 = -p1 * oszto2;
        be0[0] = -(p0 + C2 * p2);
        be0[1] = -C2;
        be1[0] = -p2 * oszto2;
        be1[1] = -oszto2;
    }

    //****************************************************************
    static inline void math_inv_2x2_fv(adattipus * be0, adattipus * be1) {
    // 1 db 4x4-es mátrixot invertál, bemenet a 4 sor elsõ elemének címe
    //****************************************************************
        adattipus p0 = abs(be0[0]) < 1e-20f ? 1e20f : adattipus(1.0f) / be0[0];
        adattipus p2 = -be1[0] * p0;
        adattipus p1 = be0[1] * p0;
        adattipus p3 = be1[1] + p2 * be0[1];
        adattipus oszto2 = abs(p3) < 1e-20f ? 1e20f : adattipus(1.0f) / p3;
        adattipus C2 = -p1 * oszto2;
        be0[0] = p0 + C2 * p2;
        be0[1] = C2;
        be1[0] = p2 * oszto2;
        be1[1] = oszto2;
    }

    //****************************************************************
    static inline void math_inv_4x4_fv(adattipus * be0, adattipus * be1, adattipus * be2, adattipus * be3) {
    // 1 db 4x4-es mátrixot invertál, bemenet a 4 sor elsõ elemének címe
    //****************************************************************
        // 0. sor
        adattipus a00 = abs(be0[0]) < 1e-20f ? 1e20f : adattipus(1.0f) / be0[0];
        adattipus a01 = a00 * be0[1];
        adattipus a02 = a00 * be0[2];
        adattipus a03 = a00 * be0[3];
        adattipus a10 =        - be1[0] * a00;
        adattipus a11 = be1[1] - be1[0] * a01;
        adattipus a12 = be1[2] - be1[0] * a02;
        adattipus a13 = be1[3] - be1[0] * a03;
        adattipus a20 =        - be2[0] * a00;
        adattipus a21 = be2[1] - be2[0] * a01;
        adattipus a22 = be2[2] - be2[0] * a02;
        adattipus a23 = be2[3] - be2[0] * a03;
        adattipus a30 =        - be3[0] * a00;
        adattipus a31 = be3[1] - be3[0] * a01;
        adattipus a32 = be3[2] - be3[0] * a02;
        adattipus a33 = be3[3] - be3[0] * a03;
        // 1. sor
        a11 = abs(a11) < 1e-20f ? 1e20f : adattipus(1.0f) / a11;
        a10 *= a11;
        a12 *= a11;
        a13 *= a11;
        a00 -= a01 * a10;
        a02 -= a01 * a12;
        a03 -= a01 * a13;
        a01 *= -a11;
        a20 -= a21 * a10;
        a22 -= a21 * a12;
        a23 -= a21 * a13;
        a21 *= -a11;
        a30 -= a31 * a10;
        a32 -= a31 * a12;
        a33 -= a31 * a13;
        a31 *= -a11;
        // 2. sor
        a22 = abs(a22) < 1e-20f ? 1e20f : adattipus(1.0f) / a22;
        a20 *= a22;
        a21 *= a22;
        a23 *= a22;
        a00 -= a02 * a20;
        a01 -= a02 * a21;
        a03 -= a02 * a23;
        a02 *= -a22;
        a10 -= a12 * a20;
        a11 -= a12 * a21;
        a13 -= a12 * a23;
        a12 *= -a22;
        a30 -= a32 * a20;
        a31 -= a32 * a21;
        a33 -= a32 * a23;
        a32 *= -a22;
        // 3. sor
        a33 = abs(a33) < 1e-20f ? 1e20f : adattipus(1.0f) / a33;
        a30 *= a33;
        a31 *= a33;
        a32 *= a33;
        be0[0] = a00 - a03 * a30;
        be0[1] = a01 - a03 * a31;
        be0[2] = a02 - a03 * a32;
        be0[3] = -a03 * a33;
        be1[0] = a10 - a13 * a30;
        be1[1] = a11 - a13 * a31;
        be1[2] = a12 - a13 * a32;
        be1[3] = -a13 * a33;
        be2[0] = a20 - a23 * a30;
        be2[1] = a21 - a23 * a31;
        be2[2] = a22 - a23 * a32;
        be2[3] = -a23 * a33;
        be3[0] = a30;
        be3[1] = a31;
        be3[2] = a32;
        be3[3] = a33;
    }

    //****************************************************************
    static inline void math_ninv_4x4_fv(adattipus * be0, adattipus * be1, adattipus * be2, adattipus * be3) {
    // 1 db 4x4-es mátrixot invertál, bemenet a 4 sor elsõ elemének címe
    //****************************************************************
        // 0. sor
        adattipus a00 = abs(be0[0]) < 1e-20f ? 1e20f : adattipus(1.0f) / be0[0];
        adattipus a01 = a00 * be0[1];
        adattipus a02 = a00 * be0[2];
        adattipus a03 = a00 * be0[3];
        adattipus a10 =        - be1[0] * a00;
        adattipus a11 = be1[1] - be1[0] * a01;
        adattipus a12 = be1[2] - be1[0] * a02;
        adattipus a13 = be1[3] - be1[0] * a03;
        adattipus a20 =        - be2[0] * a00;
        adattipus a21 = be2[1] - be2[0] * a01;
        adattipus a22 = be2[2] - be2[0] * a02;
        adattipus a23 = be2[3] - be2[0] * a03;
        adattipus a30 =        - be3[0] * a00;
        adattipus a31 = be3[1] - be3[0] * a01;
        adattipus a32 = be3[2] - be3[0] * a02;
        adattipus a33 = be3[3] - be3[0] * a03;
        // 1. sor
        a11 = abs(a11) < 1e-20f ? 1e20f : adattipus(1.0f) / a11;
        a10 *= a11;
        a12 *= a11;
        a13 *= a11;
        a00 -= a01 * a10;
        a02 -= a01 * a12;
        a03 -= a01 * a13;
        a01 *= -a11;
        a20 -= a21 * a10;
        a22 -= a21 * a12;
        a23 -= a21 * a13;
        a21 *= -a11;
        a30 -= a31 * a10;
        a32 -= a31 * a12;
        a33 -= a31 * a13;
        a31 *= -a11;
        // 2. sor
        a22 = abs(a22) < 1e-20f ? 1e20f : adattipus(1.0f) / a22;
        a20 *= a22;
        a21 *= a22;
        a23 *= a22;
        a00 -= a02 * a20;
        a01 -= a02 * a21;
        a03 -= a02 * a23;
        a02 *= -a22;
        a10 -= a12 * a20;
        a11 -= a12 * a21;
        a13 -= a12 * a23;
        a12 *= -a22;
        a30 -= a32 * a20;
        a31 -= a32 * a21;
        a33 -= a32 * a23;
        a32 *= -a22;
        // 3. sor
        a33 = abs(a33) < 1e-20f ? 1e20f : adattipus(1.0f) / a33;
        a30 *= a33;
        a31 *= a33;
        a32 *= a33;
        be0[0] = a03 * a30 - a00;
        be0[1] = a03 * a31 - a01;
        be0[2] = a03 * a32 - a02;
        be0[3] = a03 * a33;
        be1[0] = a13 * a30 - a10;
        be1[1] = a13 * a31 - a11;
        be1[2] = a13 * a32 - a12;
        be1[3] = a13 * a33;
        be2[0] = a23 * a30 - a20;
        be2[1] = a23 * a31 - a21;
        be2[2] = a23 * a32 - a22;
        be2[3] = a23 * a33;
        be3[0] = -a30;
        be3[1] = -a31;
        be3[2] = -a32;
        be3[3] = -a33;
    }

    //****************************************************************
    void math_ninv_np/*_blokk_4x4*/() {
    // fõelemkiválasztás nélkül
    //****************************************************************
        egyforma_e_hiba(row, col, "matrix::math_ninv_np");
        igaz_e_hiba(is_szimm, "matrix::math_ninv_np", "symmetrical matrix not allowed");

        if (row == 1) {
            sorok.unsafe(0).unsafe(0) = adattipus(-1) / sorok.unsafe(0).unsafe(0);
            return;
        }
        if (row == 2) {
            math_ninv_2x2_fv(&sorok.unsafe(0).unsafe(0), &sorok.unsafe(1).unsafe(0));
            return;
        }
        if (row == 4) {
            math_ninv_4x4_fv(&sorok.unsafe(0).unsafe(0), &sorok.unsafe(1).unsafe(0), &sorok.unsafe(2).unsafe(0), &sorok.unsafe(3).unsafe(0));
            return;
        }

        cuns drow = row % 4;
        cuns hrow = row - drow;
        cuns nrow = row;
        for (meret_t i = 0; i < hrow; i += 4){
            vektor<adattipus> & sor_i0 = sorok.unsafe(i);
            vektor<adattipus> & sor_i1 = sorok.unsafe(i + 1);
            vektor<adattipus> & sor_i2 = sorok.unsafe(i + 2);
            vektor<adattipus> & sor_i3 = sorok.unsafe(i + 3);

            // (i,i)-nél lévõ elem inverze
            adattipus a0[4], a1[4], a2[4], a3[4];
            a0[0] = sor_i0.unsafe(i);    a0[1] = sor_i0.unsafe(i + 1);    a0[2] = sor_i0.unsafe(i + 2);    a0[3] = sor_i0.unsafe(i + 3);
            a1[0] = sor_i1.unsafe(i);    a1[1] = sor_i1.unsafe(i + 1);    a1[2] = sor_i1.unsafe(i + 2);    a1[3] = sor_i1.unsafe(i + 3);
            a2[0] = sor_i2.unsafe(i);    a2[1] = sor_i2.unsafe(i + 1);    a2[2] = sor_i2.unsafe(i + 2);    a2[3] = sor_i2.unsafe(i + 3);
            a3[0] = sor_i3.unsafe(i);    a3[1] = sor_i3.unsafe(i + 1);    a3[2] = sor_i3.unsafe(i + 2);    a3[3] = sor_i3.unsafe(i + 3);
            math_inv_4x4_fv(a0, a1, a2, a3); // a0...a3 az inverz

            meret_t j;
            for (j = 0; j < i; j += 4) {
                vektor<adattipus> & sor_j0 = sorok.unsafe(j);
                vektor<adattipus> & sor_j1 = sorok.unsafe(j + 1);
                vektor<adattipus> & sor_j2 = sorok.unsafe(j + 2);
                vektor<adattipus> & sor_j3 = sorok.unsafe(j + 3);
                const adattipus C00 = sor_j0.unsafe(i) * a0[0] + sor_j0.unsafe(i + 1) * a1[0] + sor_j0.unsafe(i + 2) * a2[0] + sor_j0.unsafe(i + 3) * a3[0];
                const adattipus C01 = sor_j0.unsafe(i) * a0[1] + sor_j0.unsafe(i + 1) * a1[1] + sor_j0.unsafe(i + 2) * a2[1] + sor_j0.unsafe(i + 3) * a3[1];
                const adattipus C02 = sor_j0.unsafe(i) * a0[2] + sor_j0.unsafe(i + 1) * a1[2] + sor_j0.unsafe(i + 2) * a2[2] + sor_j0.unsafe(i + 3) * a3[2];
                const adattipus C03 = sor_j0.unsafe(i) * a0[3] + sor_j0.unsafe(i + 1) * a1[3] + sor_j0.unsafe(i + 2) * a2[3] + sor_j0.unsafe(i + 3) * a3[3];
                const adattipus C10 = sor_j1.unsafe(i) * a0[0] + sor_j1.unsafe(i + 1) * a1[0] + sor_j1.unsafe(i + 2) * a2[0] + sor_j1.unsafe(i + 3) * a3[0];
                const adattipus C11 = sor_j1.unsafe(i) * a0[1] + sor_j1.unsafe(i + 1) * a1[1] + sor_j1.unsafe(i + 2) * a2[1] + sor_j1.unsafe(i + 3) * a3[1];
                const adattipus C12 = sor_j1.unsafe(i) * a0[2] + sor_j1.unsafe(i + 1) * a1[2] + sor_j1.unsafe(i + 2) * a2[2] + sor_j1.unsafe(i + 3) * a3[2];
                const adattipus C13 = sor_j1.unsafe(i) * a0[3] + sor_j1.unsafe(i + 1) * a1[3] + sor_j1.unsafe(i + 2) * a2[3] + sor_j1.unsafe(i + 3) * a3[3];
                const adattipus C20 = sor_j2.unsafe(i) * a0[0] + sor_j2.unsafe(i + 1) * a1[0] + sor_j2.unsafe(i + 2) * a2[0] + sor_j2.unsafe(i + 3) * a3[0];
                const adattipus C21 = sor_j2.unsafe(i) * a0[1] + sor_j2.unsafe(i + 1) * a1[1] + sor_j2.unsafe(i + 2) * a2[1] + sor_j2.unsafe(i + 3) * a3[1];
                const adattipus C22 = sor_j2.unsafe(i) * a0[2] + sor_j2.unsafe(i + 1) * a1[2] + sor_j2.unsafe(i + 2) * a2[2] + sor_j2.unsafe(i + 3) * a3[2];
                const adattipus C23 = sor_j2.unsafe(i) * a0[3] + sor_j2.unsafe(i + 1) * a1[3] + sor_j2.unsafe(i + 2) * a2[3] + sor_j2.unsafe(i + 3) * a3[3];
                const adattipus C30 = sor_j3.unsafe(i) * a0[0] + sor_j3.unsafe(i + 1) * a1[0] + sor_j3.unsafe(i + 2) * a2[0] + sor_j3.unsafe(i + 3) * a3[0];
                const adattipus C31 = sor_j3.unsafe(i) * a0[1] + sor_j3.unsafe(i + 1) * a1[1] + sor_j3.unsafe(i + 2) * a2[1] + sor_j3.unsafe(i + 3) * a3[1];
                const adattipus C32 = sor_j3.unsafe(i) * a0[2] + sor_j3.unsafe(i + 1) * a1[2] + sor_j3.unsafe(i + 2) * a2[2] + sor_j3.unsafe(i + 3) * a3[2];
                const adattipus C33 = sor_j3.unsafe(i) * a0[3] + sor_j3.unsafe(i + 1) * a1[3] + sor_j3.unsafe(i + 2) * a2[3] + sor_j3.unsafe(i + 3) * a3[3];
                meret_t k;
                for (k = 0; k < i; k += 4) {
                    sor_j0.unsafe(k)     -= C00 * sor_i0.unsafe(k)     + C01 * sor_i1.unsafe(k)     + C02 * sor_i2.unsafe(k)     + C03 * sor_i3.unsafe(k);
                    sor_j0.unsafe(k + 1) -= C00 * sor_i0.unsafe(k + 1) + C01 * sor_i1.unsafe(k + 1) + C02 * sor_i2.unsafe(k + 1) + C03 * sor_i3.unsafe(k + 1);
                    sor_j0.unsafe(k + 2) -= C00 * sor_i0.unsafe(k + 2) + C01 * sor_i1.unsafe(k + 2) + C02 * sor_i2.unsafe(k + 2) + C03 * sor_i3.unsafe(k + 2);
                    sor_j0.unsafe(k + 3) -= C00 * sor_i0.unsafe(k + 3) + C01 * sor_i1.unsafe(k + 3) + C02 * sor_i2.unsafe(k + 3) + C03 * sor_i3.unsafe(k + 3);
                    sor_j1.unsafe(k)     -= C10 * sor_i0.unsafe(k)     + C11 * sor_i1.unsafe(k)     + C12 * sor_i2.unsafe(k)     + C13 * sor_i3.unsafe(k);
                    sor_j1.unsafe(k + 1) -= C10 * sor_i0.unsafe(k + 1) + C11 * sor_i1.unsafe(k + 1) + C12 * sor_i2.unsafe(k + 1) + C13 * sor_i3.unsafe(k + 1);
                    sor_j1.unsafe(k + 2) -= C10 * sor_i0.unsafe(k + 2) + C11 * sor_i1.unsafe(k + 2) + C12 * sor_i2.unsafe(k + 2) + C13 * sor_i3.unsafe(k + 2);
                    sor_j1.unsafe(k + 3) -= C10 * sor_i0.unsafe(k + 3) + C11 * sor_i1.unsafe(k + 3) + C12 * sor_i2.unsafe(k + 3) + C13 * sor_i3.unsafe(k + 3);
                    sor_j2.unsafe(k)     -= C20 * sor_i0.unsafe(k)     + C21 * sor_i1.unsafe(k)     + C22 * sor_i2.unsafe(k)     + C23 * sor_i3.unsafe(k);
                    sor_j2.unsafe(k + 1) -= C20 * sor_i0.unsafe(k + 1) + C21 * sor_i1.unsafe(k + 1) + C22 * sor_i2.unsafe(k + 1) + C23 * sor_i3.unsafe(k + 1);
                    sor_j2.unsafe(k + 2) -= C20 * sor_i0.unsafe(k + 2) + C21 * sor_i1.unsafe(k + 2) + C22 * sor_i2.unsafe(k + 2) + C23 * sor_i3.unsafe(k + 2);
                    sor_j2.unsafe(k + 3) -= C20 * sor_i0.unsafe(k + 3) + C21 * sor_i1.unsafe(k + 3) + C22 * sor_i2.unsafe(k + 3) + C23 * sor_i3.unsafe(k + 3);
                    sor_j3.unsafe(k)     -= C30 * sor_i0.unsafe(k)     + C31 * sor_i1.unsafe(k)     + C32 * sor_i2.unsafe(k)     + C33 * sor_i3.unsafe(k);
                    sor_j3.unsafe(k + 1) -= C30 * sor_i0.unsafe(k + 1) + C31 * sor_i1.unsafe(k + 1) + C32 * sor_i2.unsafe(k + 1) + C33 * sor_i3.unsafe(k + 1);
                    sor_j3.unsafe(k + 2) -= C30 * sor_i0.unsafe(k + 2) + C31 * sor_i1.unsafe(k + 2) + C32 * sor_i2.unsafe(k + 2) + C33 * sor_i3.unsafe(k + 2);
                    sor_j3.unsafe(k + 3) -= C30 * sor_i0.unsafe(k + 3) + C31 * sor_i1.unsafe(k + 3) + C32 * sor_i2.unsafe(k + 3) + C33 * sor_i3.unsafe(k + 3);
                }
                sor_j0.unsafe(k)     = C00;    sor_j0.unsafe(k + 1) = C01;    sor_j0.unsafe(k + 2) = C02;    sor_j0.unsafe(k + 3) = C03; // nem negált invnél = -C;
                sor_j1.unsafe(k)     = C10;    sor_j1.unsafe(k + 1) = C11;    sor_j1.unsafe(k + 2) = C12;    sor_j1.unsafe(k + 3) = C13;
                sor_j2.unsafe(k)     = C20;    sor_j2.unsafe(k + 1) = C21;    sor_j2.unsafe(k + 2) = C22;    sor_j2.unsafe(k + 3) = C23;
                sor_j3.unsafe(k)     = C30;    sor_j3.unsafe(k + 1) = C31;    sor_j3.unsafe(k + 2) = C32;    sor_j3.unsafe(k + 3) = C33;
                for (k += 4; k < hrow; k += 4) {
                    sor_j0.unsafe(k)     -= C00 * sor_i0.unsafe(k)     + C01 * sor_i1.unsafe(k)     + C02 * sor_i2.unsafe(k)     + C03 * sor_i3.unsafe(k);
                    sor_j0.unsafe(k + 1) -= C00 * sor_i0.unsafe(k + 1) + C01 * sor_i1.unsafe(k + 1) + C02 * sor_i2.unsafe(k + 1) + C03 * sor_i3.unsafe(k + 1);
                    sor_j0.unsafe(k + 2) -= C00 * sor_i0.unsafe(k + 2) + C01 * sor_i1.unsafe(k + 2) + C02 * sor_i2.unsafe(k + 2) + C03 * sor_i3.unsafe(k + 2);
                    sor_j0.unsafe(k + 3) -= C00 * sor_i0.unsafe(k + 3) + C01 * sor_i1.unsafe(k + 3) + C02 * sor_i2.unsafe(k + 3) + C03 * sor_i3.unsafe(k + 3);
                    sor_j1.unsafe(k)     -= C10 * sor_i0.unsafe(k)     + C11 * sor_i1.unsafe(k)     + C12 * sor_i2.unsafe(k)     + C13 * sor_i3.unsafe(k);
                    sor_j1.unsafe(k + 1) -= C10 * sor_i0.unsafe(k + 1) + C11 * sor_i1.unsafe(k + 1) + C12 * sor_i2.unsafe(k + 1) + C13 * sor_i3.unsafe(k + 1);
                    sor_j1.unsafe(k + 2) -= C10 * sor_i0.unsafe(k + 2) + C11 * sor_i1.unsafe(k + 2) + C12 * sor_i2.unsafe(k + 2) + C13 * sor_i3.unsafe(k + 2);
                    sor_j1.unsafe(k + 3) -= C10 * sor_i0.unsafe(k + 3) + C11 * sor_i1.unsafe(k + 3) + C12 * sor_i2.unsafe(k + 3) + C13 * sor_i3.unsafe(k + 3);
                    sor_j2.unsafe(k)     -= C20 * sor_i0.unsafe(k)     + C21 * sor_i1.unsafe(k)     + C22 * sor_i2.unsafe(k)     + C23 * sor_i3.unsafe(k);
                    sor_j2.unsafe(k + 1) -= C20 * sor_i0.unsafe(k + 1) + C21 * sor_i1.unsafe(k + 1) + C22 * sor_i2.unsafe(k + 1) + C23 * sor_i3.unsafe(k + 1);
                    sor_j2.unsafe(k + 2) -= C20 * sor_i0.unsafe(k + 2) + C21 * sor_i1.unsafe(k + 2) + C22 * sor_i2.unsafe(k + 2) + C23 * sor_i3.unsafe(k + 2);
                    sor_j2.unsafe(k + 3) -= C20 * sor_i0.unsafe(k + 3) + C21 * sor_i1.unsafe(k + 3) + C22 * sor_i2.unsafe(k + 3) + C23 * sor_i3.unsafe(k + 3);
                    sor_j3.unsafe(k)     -= C30 * sor_i0.unsafe(k)     + C31 * sor_i1.unsafe(k)     + C32 * sor_i2.unsafe(k)     + C33 * sor_i3.unsafe(k);
                    sor_j3.unsafe(k + 1) -= C30 * sor_i0.unsafe(k + 1) + C31 * sor_i1.unsafe(k + 1) + C32 * sor_i2.unsafe(k + 1) + C33 * sor_i3.unsafe(k + 1);
                    sor_j3.unsafe(k + 2) -= C30 * sor_i0.unsafe(k + 2) + C31 * sor_i1.unsafe(k + 2) + C32 * sor_i2.unsafe(k + 2) + C33 * sor_i3.unsafe(k + 2);
                    sor_j3.unsafe(k + 3) -= C30 * sor_i0.unsafe(k + 3) + C31 * sor_i1.unsafe(k + 3) + C32 * sor_i2.unsafe(k + 3) + C33 * sor_i3.unsafe(k + 3);
                }
                for (k = hrow; k < nrow; k++) {
                    sor_j0.unsafe(k)     -= C00 * sor_i0.unsafe(k)     + C01 * sor_i1.unsafe(k)     + C02 * sor_i2.unsafe(k)     + C03 * sor_i3.unsafe(k);
                    sor_j1.unsafe(k)     -= C10 * sor_i0.unsafe(k)     + C11 * sor_i1.unsafe(k)     + C12 * sor_i2.unsafe(k)     + C13 * sor_i3.unsafe(k);
                    sor_j2.unsafe(k)     -= C20 * sor_i0.unsafe(k)     + C21 * sor_i1.unsafe(k)     + C22 * sor_i2.unsafe(k)     + C23 * sor_i3.unsafe(k);
                    sor_j3.unsafe(k)     -= C30 * sor_i0.unsafe(k)     + C31 * sor_i1.unsafe(k)     + C32 * sor_i2.unsafe(k)     + C33 * sor_i3.unsafe(k);
                }
            }
            for (j+=4; j < hrow; j += 4) {
                vektor<adattipus> & sor_j0 = sorok.unsafe(j);
                vektor<adattipus> & sor_j1 = sorok.unsafe(j + 1);
                vektor<adattipus> & sor_j2 = sorok.unsafe(j + 2);
                vektor<adattipus> & sor_j3 = sorok.unsafe(j + 3);
                const adattipus C00 = sor_j0.unsafe(i) * a0[0] + sor_j0.unsafe(i + 1) * a1[0] + sor_j0.unsafe(i + 2) * a2[0] + sor_j0.unsafe(i + 3) * a3[0];
                const adattipus C01 = sor_j0.unsafe(i) * a0[1] + sor_j0.unsafe(i + 1) * a1[1] + sor_j0.unsafe(i + 2) * a2[1] + sor_j0.unsafe(i + 3) * a3[1];
                const adattipus C02 = sor_j0.unsafe(i) * a0[2] + sor_j0.unsafe(i + 1) * a1[2] + sor_j0.unsafe(i + 2) * a2[2] + sor_j0.unsafe(i + 3) * a3[2];
                const adattipus C03 = sor_j0.unsafe(i) * a0[3] + sor_j0.unsafe(i + 1) * a1[3] + sor_j0.unsafe(i + 2) * a2[3] + sor_j0.unsafe(i + 3) * a3[3];
                const adattipus C10 = sor_j1.unsafe(i) * a0[0] + sor_j1.unsafe(i + 1) * a1[0] + sor_j1.unsafe(i + 2) * a2[0] + sor_j1.unsafe(i + 3) * a3[0];
                const adattipus C11 = sor_j1.unsafe(i) * a0[1] + sor_j1.unsafe(i + 1) * a1[1] + sor_j1.unsafe(i + 2) * a2[1] + sor_j1.unsafe(i + 3) * a3[1];
                const adattipus C12 = sor_j1.unsafe(i) * a0[2] + sor_j1.unsafe(i + 1) * a1[2] + sor_j1.unsafe(i + 2) * a2[2] + sor_j1.unsafe(i + 3) * a3[2];
                const adattipus C13 = sor_j1.unsafe(i) * a0[3] + sor_j1.unsafe(i + 1) * a1[3] + sor_j1.unsafe(i + 2) * a2[3] + sor_j1.unsafe(i + 3) * a3[3];
                const adattipus C20 = sor_j2.unsafe(i) * a0[0] + sor_j2.unsafe(i + 1) * a1[0] + sor_j2.unsafe(i + 2) * a2[0] + sor_j2.unsafe(i + 3) * a3[0];
                const adattipus C21 = sor_j2.unsafe(i) * a0[1] + sor_j2.unsafe(i + 1) * a1[1] + sor_j2.unsafe(i + 2) * a2[1] + sor_j2.unsafe(i + 3) * a3[1];
                const adattipus C22 = sor_j2.unsafe(i) * a0[2] + sor_j2.unsafe(i + 1) * a1[2] + sor_j2.unsafe(i + 2) * a2[2] + sor_j2.unsafe(i + 3) * a3[2];
                const adattipus C23 = sor_j2.unsafe(i) * a0[3] + sor_j2.unsafe(i + 1) * a1[3] + sor_j2.unsafe(i + 2) * a2[3] + sor_j2.unsafe(i + 3) * a3[3];
                const adattipus C30 = sor_j3.unsafe(i) * a0[0] + sor_j3.unsafe(i + 1) * a1[0] + sor_j3.unsafe(i + 2) * a2[0] + sor_j3.unsafe(i + 3) * a3[0];
                const adattipus C31 = sor_j3.unsafe(i) * a0[1] + sor_j3.unsafe(i + 1) * a1[1] + sor_j3.unsafe(i + 2) * a2[1] + sor_j3.unsafe(i + 3) * a3[1];
                const adattipus C32 = sor_j3.unsafe(i) * a0[2] + sor_j3.unsafe(i + 1) * a1[2] + sor_j3.unsafe(i + 2) * a2[2] + sor_j3.unsafe(i + 3) * a3[2];
                const adattipus C33 = sor_j3.unsafe(i) * a0[3] + sor_j3.unsafe(i + 1) * a1[3] + sor_j3.unsafe(i + 2) * a2[3] + sor_j3.unsafe(i + 3) * a3[3];
                meret_t k;
                for (k = 0; k < i; k += 4) {
                    sor_j0.unsafe(k)     -= C00 * sor_i0.unsafe(k)     + C01 * sor_i1.unsafe(k)     + C02 * sor_i2.unsafe(k)     + C03 * sor_i3.unsafe(k);
                    sor_j0.unsafe(k + 1) -= C00 * sor_i0.unsafe(k + 1) + C01 * sor_i1.unsafe(k + 1) + C02 * sor_i2.unsafe(k + 1) + C03 * sor_i3.unsafe(k + 1);
                    sor_j0.unsafe(k + 2) -= C00 * sor_i0.unsafe(k + 2) + C01 * sor_i1.unsafe(k + 2) + C02 * sor_i2.unsafe(k + 2) + C03 * sor_i3.unsafe(k + 2);
                    sor_j0.unsafe(k + 3) -= C00 * sor_i0.unsafe(k + 3) + C01 * sor_i1.unsafe(k + 3) + C02 * sor_i2.unsafe(k + 3) + C03 * sor_i3.unsafe(k + 3);
                    sor_j1.unsafe(k)     -= C10 * sor_i0.unsafe(k)     + C11 * sor_i1.unsafe(k)     + C12 * sor_i2.unsafe(k)     + C13 * sor_i3.unsafe(k);
                    sor_j1.unsafe(k + 1) -= C10 * sor_i0.unsafe(k + 1) + C11 * sor_i1.unsafe(k + 1) + C12 * sor_i2.unsafe(k + 1) + C13 * sor_i3.unsafe(k + 1);
                    sor_j1.unsafe(k + 2) -= C10 * sor_i0.unsafe(k + 2) + C11 * sor_i1.unsafe(k + 2) + C12 * sor_i2.unsafe(k + 2) + C13 * sor_i3.unsafe(k + 2);
                    sor_j1.unsafe(k + 3) -= C10 * sor_i0.unsafe(k + 3) + C11 * sor_i1.unsafe(k + 3) + C12 * sor_i2.unsafe(k + 3) + C13 * sor_i3.unsafe(k + 3);
                    sor_j2.unsafe(k)     -= C20 * sor_i0.unsafe(k)     + C21 * sor_i1.unsafe(k)     + C22 * sor_i2.unsafe(k)     + C23 * sor_i3.unsafe(k);
                    sor_j2.unsafe(k + 1) -= C20 * sor_i0.unsafe(k + 1) + C21 * sor_i1.unsafe(k + 1) + C22 * sor_i2.unsafe(k + 1) + C23 * sor_i3.unsafe(k + 1);
                    sor_j2.unsafe(k + 2) -= C20 * sor_i0.unsafe(k + 2) + C21 * sor_i1.unsafe(k + 2) + C22 * sor_i2.unsafe(k + 2) + C23 * sor_i3.unsafe(k + 2);
                    sor_j2.unsafe(k + 3) -= C20 * sor_i0.unsafe(k + 3) + C21 * sor_i1.unsafe(k + 3) + C22 * sor_i2.unsafe(k + 3) + C23 * sor_i3.unsafe(k + 3);
                    sor_j3.unsafe(k)     -= C30 * sor_i0.unsafe(k)     + C31 * sor_i1.unsafe(k)     + C32 * sor_i2.unsafe(k)     + C33 * sor_i3.unsafe(k);
                    sor_j3.unsafe(k + 1) -= C30 * sor_i0.unsafe(k + 1) + C31 * sor_i1.unsafe(k + 1) + C32 * sor_i2.unsafe(k + 1) + C33 * sor_i3.unsafe(k + 1);
                    sor_j3.unsafe(k + 2) -= C30 * sor_i0.unsafe(k + 2) + C31 * sor_i1.unsafe(k + 2) + C32 * sor_i2.unsafe(k + 2) + C33 * sor_i3.unsafe(k + 2);
                    sor_j3.unsafe(k + 3) -= C30 * sor_i0.unsafe(k + 3) + C31 * sor_i1.unsafe(k + 3) + C32 * sor_i2.unsafe(k + 3) + C33 * sor_i3.unsafe(k + 3);
                }
                sor_j0.unsafe(k)     = C00;    sor_j0.unsafe(k + 1) = C01;    sor_j0.unsafe(k + 2) = C02;    sor_j0.unsafe(k + 3) = C03; // nem negált invnél = -C;
                sor_j1.unsafe(k)     = C10;    sor_j1.unsafe(k + 1) = C11;    sor_j1.unsafe(k + 2) = C12;    sor_j1.unsafe(k + 3) = C13;
                sor_j2.unsafe(k)     = C20;    sor_j2.unsafe(k + 1) = C21;    sor_j2.unsafe(k + 2) = C22;    sor_j2.unsafe(k + 3) = C23;
                sor_j3.unsafe(k)     = C30;    sor_j3.unsafe(k + 1) = C31;    sor_j3.unsafe(k + 2) = C32;    sor_j3.unsafe(k + 3) = C33;
                for (k += 4; k < hrow; k += 4) {
                    sor_j0.unsafe(k)     -= C00 * sor_i0.unsafe(k)     + C01 * sor_i1.unsafe(k)     + C02 * sor_i2.unsafe(k)     + C03 * sor_i3.unsafe(k);
                    sor_j0.unsafe(k + 1) -= C00 * sor_i0.unsafe(k + 1) + C01 * sor_i1.unsafe(k + 1) + C02 * sor_i2.unsafe(k + 1) + C03 * sor_i3.unsafe(k + 1);
                    sor_j0.unsafe(k + 2) -= C00 * sor_i0.unsafe(k + 2) + C01 * sor_i1.unsafe(k + 2) + C02 * sor_i2.unsafe(k + 2) + C03 * sor_i3.unsafe(k + 2);
                    sor_j0.unsafe(k + 3) -= C00 * sor_i0.unsafe(k + 3) + C01 * sor_i1.unsafe(k + 3) + C02 * sor_i2.unsafe(k + 3) + C03 * sor_i3.unsafe(k + 3);
                    sor_j1.unsafe(k)     -= C10 * sor_i0.unsafe(k)     + C11 * sor_i1.unsafe(k)     + C12 * sor_i2.unsafe(k)     + C13 * sor_i3.unsafe(k);
                    sor_j1.unsafe(k + 1) -= C10 * sor_i0.unsafe(k + 1) + C11 * sor_i1.unsafe(k + 1) + C12 * sor_i2.unsafe(k + 1) + C13 * sor_i3.unsafe(k + 1);
                    sor_j1.unsafe(k + 2) -= C10 * sor_i0.unsafe(k + 2) + C11 * sor_i1.unsafe(k + 2) + C12 * sor_i2.unsafe(k + 2) + C13 * sor_i3.unsafe(k + 2);
                    sor_j1.unsafe(k + 3) -= C10 * sor_i0.unsafe(k + 3) + C11 * sor_i1.unsafe(k + 3) + C12 * sor_i2.unsafe(k + 3) + C13 * sor_i3.unsafe(k + 3);
                    sor_j2.unsafe(k)     -= C20 * sor_i0.unsafe(k)     + C21 * sor_i1.unsafe(k)     + C22 * sor_i2.unsafe(k)     + C23 * sor_i3.unsafe(k);
                    sor_j2.unsafe(k + 1) -= C20 * sor_i0.unsafe(k + 1) + C21 * sor_i1.unsafe(k + 1) + C22 * sor_i2.unsafe(k + 1) + C23 * sor_i3.unsafe(k + 1);
                    sor_j2.unsafe(k + 2) -= C20 * sor_i0.unsafe(k + 2) + C21 * sor_i1.unsafe(k + 2) + C22 * sor_i2.unsafe(k + 2) + C23 * sor_i3.unsafe(k + 2);
                    sor_j2.unsafe(k + 3) -= C20 * sor_i0.unsafe(k + 3) + C21 * sor_i1.unsafe(k + 3) + C22 * sor_i2.unsafe(k + 3) + C23 * sor_i3.unsafe(k + 3);
                    sor_j3.unsafe(k)     -= C30 * sor_i0.unsafe(k)     + C31 * sor_i1.unsafe(k)     + C32 * sor_i2.unsafe(k)     + C33 * sor_i3.unsafe(k);
                    sor_j3.unsafe(k + 1) -= C30 * sor_i0.unsafe(k + 1) + C31 * sor_i1.unsafe(k + 1) + C32 * sor_i2.unsafe(k + 1) + C33 * sor_i3.unsafe(k + 1);
                    sor_j3.unsafe(k + 2) -= C30 * sor_i0.unsafe(k + 2) + C31 * sor_i1.unsafe(k + 2) + C32 * sor_i2.unsafe(k + 2) + C33 * sor_i3.unsafe(k + 2);
                    sor_j3.unsafe(k + 3) -= C30 * sor_i0.unsafe(k + 3) + C31 * sor_i1.unsafe(k + 3) + C32 * sor_i2.unsafe(k + 3) + C33 * sor_i3.unsafe(k + 3);
                }
                for (k = hrow; k < nrow; k++) {
                    sor_j0.unsafe(k)     -= C00 * sor_i0.unsafe(k)     + C01 * sor_i1.unsafe(k)     + C02 * sor_i2.unsafe(k)     + C03 * sor_i3.unsafe(k);
                    sor_j1.unsafe(k)     -= C10 * sor_i0.unsafe(k)     + C11 * sor_i1.unsafe(k)     + C12 * sor_i2.unsafe(k)     + C13 * sor_i3.unsafe(k);
                    sor_j2.unsafe(k)     -= C20 * sor_i0.unsafe(k)     + C21 * sor_i1.unsafe(k)     + C22 * sor_i2.unsafe(k)     + C23 * sor_i3.unsafe(k);
                    sor_j3.unsafe(k)     -= C30 * sor_i0.unsafe(k)     + C31 * sor_i1.unsafe(k)     + C32 * sor_i2.unsafe(k)     + C33 * sor_i3.unsafe(k);
                }
            }
            for (j = hrow; j < nrow; j++) {
                vektor<adattipus> & sor_j0 = sorok.unsafe(j);
                const adattipus C00 = sor_j0.unsafe(i) * a0[0] + sor_j0.unsafe(i + 1) * a1[0] + sor_j0.unsafe(i + 2) * a2[0] + sor_j0.unsafe(i + 3) * a3[0];
                const adattipus C01 = sor_j0.unsafe(i) * a0[1] + sor_j0.unsafe(i + 1) * a1[1] + sor_j0.unsafe(i + 2) * a2[1] + sor_j0.unsafe(i + 3) * a3[1];
                const adattipus C02 = sor_j0.unsafe(i) * a0[2] + sor_j0.unsafe(i + 1) * a1[2] + sor_j0.unsafe(i + 2) * a2[2] + sor_j0.unsafe(i + 3) * a3[2];
                const adattipus C03 = sor_j0.unsafe(i) * a0[3] + sor_j0.unsafe(i + 1) * a1[3] + sor_j0.unsafe(i + 2) * a2[3] + sor_j0.unsafe(i + 3) * a3[3];
                meret_t k;
                for (k = 0; k < i; k += 4) {
                    sor_j0.unsafe(k)     -= C00 * sor_i0.unsafe(k)     + C01 * sor_i1.unsafe(k)     + C02 * sor_i2.unsafe(k)     + C03 * sor_i3.unsafe(k);
                    sor_j0.unsafe(k + 1) -= C00 * sor_i0.unsafe(k + 1) + C01 * sor_i1.unsafe(k + 1) + C02 * sor_i2.unsafe(k + 1) + C03 * sor_i3.unsafe(k + 1);
                    sor_j0.unsafe(k + 2) -= C00 * sor_i0.unsafe(k + 2) + C01 * sor_i1.unsafe(k + 2) + C02 * sor_i2.unsafe(k + 2) + C03 * sor_i3.unsafe(k + 2);
                    sor_j0.unsafe(k + 3) -= C00 * sor_i0.unsafe(k + 3) + C01 * sor_i1.unsafe(k + 3) + C02 * sor_i2.unsafe(k + 3) + C03 * sor_i3.unsafe(k + 3);
                }
                sor_j0.unsafe(k) = C00;    sor_j0.unsafe(k + 1) = C01;    sor_j0.unsafe(k + 2) = C02;    sor_j0.unsafe(k + 3) = C03; // nem negált invnél = -C;
                for (k += 4; k < hrow; k += 4) {
                    sor_j0.unsafe(k)     -= C00 * sor_i0.unsafe(k)     + C01 * sor_i1.unsafe(k)     + C02 * sor_i2.unsafe(k)     + C03 * sor_i3.unsafe(k);
                    sor_j0.unsafe(k + 1) -= C00 * sor_i0.unsafe(k + 1) + C01 * sor_i1.unsafe(k + 1) + C02 * sor_i2.unsafe(k + 1) + C03 * sor_i3.unsafe(k + 1);
                    sor_j0.unsafe(k + 2) -= C00 * sor_i0.unsafe(k + 2) + C01 * sor_i1.unsafe(k + 2) + C02 * sor_i2.unsafe(k + 2) + C03 * sor_i3.unsafe(k + 2);
                    sor_j0.unsafe(k + 3) -= C00 * sor_i0.unsafe(k + 3) + C01 * sor_i1.unsafe(k + 3) + C02 * sor_i2.unsafe(k + 3) + C03 * sor_i3.unsafe(k + 3);
                }
                for (k = hrow; k < nrow; k++) {
                    sor_j0.unsafe(k)     -= C00 * sor_i0.unsafe(k)     + C01 * sor_i1.unsafe(k)     + C02 * sor_i2.unsafe(k)     + C03 * sor_i3.unsafe(k);
                }
            }
            for (j = 0; j < i; j += 4) {
                const adattipus b00 = a0[0] * sor_i0.unsafe(j)     + a0[1] * sor_i1.unsafe(j)     + a0[2] * sor_i2.unsafe(j)     + a0[3] * sor_i3.unsafe(j);
                const adattipus b01 = a0[0] * sor_i0.unsafe(j + 1) + a0[1] * sor_i1.unsafe(j + 1) + a0[2] * sor_i2.unsafe(j + 1) + a0[3] * sor_i3.unsafe(j + 1);
                const adattipus b02 = a0[0] * sor_i0.unsafe(j + 2) + a0[1] * sor_i1.unsafe(j + 2) + a0[2] * sor_i2.unsafe(j + 2) + a0[3] * sor_i3.unsafe(j + 2);
                const adattipus b03 = a0[0] * sor_i0.unsafe(j + 3) + a0[1] * sor_i1.unsafe(j + 3) + a0[2] * sor_i2.unsafe(j + 3) + a0[3] * sor_i3.unsafe(j + 3);
                const adattipus b10 = a1[0] * sor_i0.unsafe(j)     + a1[1] * sor_i1.unsafe(j)     + a1[2] * sor_i2.unsafe(j)     + a1[3] * sor_i3.unsafe(j);
                const adattipus b11 = a1[0] * sor_i0.unsafe(j + 1) + a1[1] * sor_i1.unsafe(j + 1) + a1[2] * sor_i2.unsafe(j + 1) + a1[3] * sor_i3.unsafe(j + 1);
                const adattipus b12 = a1[0] * sor_i0.unsafe(j + 2) + a1[1] * sor_i1.unsafe(j + 2) + a1[2] * sor_i2.unsafe(j + 2) + a1[3] * sor_i3.unsafe(j + 2);
                const adattipus b13 = a1[0] * sor_i0.unsafe(j + 3) + a1[1] * sor_i1.unsafe(j + 3) + a1[2] * sor_i2.unsafe(j + 3) + a1[3] * sor_i3.unsafe(j + 3);
                const adattipus b20 = a2[0] * sor_i0.unsafe(j)     + a2[1] * sor_i1.unsafe(j)     + a2[2] * sor_i2.unsafe(j)     + a2[3] * sor_i3.unsafe(j);
                const adattipus b21 = a2[0] * sor_i0.unsafe(j + 1) + a2[1] * sor_i1.unsafe(j + 1) + a2[2] * sor_i2.unsafe(j + 1) + a2[3] * sor_i3.unsafe(j + 1);
                const adattipus b22 = a2[0] * sor_i0.unsafe(j + 2) + a2[1] * sor_i1.unsafe(j + 2) + a2[2] * sor_i2.unsafe(j + 2) + a2[3] * sor_i3.unsafe(j + 2);
                const adattipus b23 = a2[0] * sor_i0.unsafe(j + 3) + a2[1] * sor_i1.unsafe(j + 3) + a2[2] * sor_i2.unsafe(j + 3) + a2[3] * sor_i3.unsafe(j + 3);
                const adattipus b30 = a3[0] * sor_i0.unsafe(j)     + a3[1] * sor_i1.unsafe(j)     + a3[2] * sor_i2.unsafe(j)     + a3[3] * sor_i3.unsafe(j);
                const adattipus b31 = a3[0] * sor_i0.unsafe(j + 1) + a3[1] * sor_i1.unsafe(j + 1) + a3[2] * sor_i2.unsafe(j + 1) + a3[3] * sor_i3.unsafe(j + 1);
                const adattipus b32 = a3[0] * sor_i0.unsafe(j + 2) + a3[1] * sor_i1.unsafe(j + 2) + a3[2] * sor_i2.unsafe(j + 2) + a3[3] * sor_i3.unsafe(j + 2);
                const adattipus b33 = a3[0] * sor_i0.unsafe(j + 3) + a3[1] * sor_i1.unsafe(j + 3) + a3[2] * sor_i2.unsafe(j + 3) + a3[3] * sor_i3.unsafe(j + 3);
                sor_i0.unsafe(j)     = b00;   sor_i0.unsafe(j + 1) = b01;    sor_i0.unsafe(j + 2) = b02;    sor_i0.unsafe(j + 3) = b03;
                sor_i1.unsafe(j)     = b10;   sor_i1.unsafe(j + 1) = b11;    sor_i1.unsafe(j + 2) = b12;    sor_i1.unsafe(j + 3) = b13;
                sor_i2.unsafe(j)     = b20;   sor_i2.unsafe(j + 1) = b21;    sor_i2.unsafe(j + 2) = b22;    sor_i2.unsafe(j + 3) = b23;
                sor_i3.unsafe(j)     = b30;   sor_i3.unsafe(j + 1) = b31;    sor_i3.unsafe(j + 2) = b32;    sor_i3.unsafe(j + 3) = b33;
            }
            sor_i0.unsafe(j) = -a0[0];   sor_i0.unsafe(j + 1) = -a0[1];    sor_i0.unsafe(j + 2) = -a0[2];    sor_i0.unsafe(j + 3) = -a0[3]; // nem negált invnél = ai[j];
            sor_i1.unsafe(j) = -a1[0];   sor_i1.unsafe(j + 1) = -a1[1];    sor_i1.unsafe(j + 2) = -a1[2];    sor_i1.unsafe(j + 3) = -a1[3];
            sor_i2.unsafe(j) = -a2[0];   sor_i2.unsafe(j + 1) = -a2[1];    sor_i2.unsafe(j + 2) = -a2[2];    sor_i2.unsafe(j + 3) = -a2[3];
            sor_i3.unsafe(j) = -a3[0];   sor_i3.unsafe(j + 1) = -a3[1];    sor_i3.unsafe(j + 2) = -a3[2];    sor_i3.unsafe(j + 3) = -a3[3];
            for (j += 4; j < hrow; j += 4) {
                const adattipus b00 = a0[0] * sor_i0.unsafe(j)     + a0[1] * sor_i1.unsafe(j)     + a0[2] * sor_i2.unsafe(j)     + a0[3] * sor_i3.unsafe(j);
                const adattipus b01 = a0[0] * sor_i0.unsafe(j + 1) + a0[1] * sor_i1.unsafe(j + 1) + a0[2] * sor_i2.unsafe(j + 1) + a0[3] * sor_i3.unsafe(j + 1);
                const adattipus b02 = a0[0] * sor_i0.unsafe(j + 2) + a0[1] * sor_i1.unsafe(j + 2) + a0[2] * sor_i2.unsafe(j + 2) + a0[3] * sor_i3.unsafe(j + 2);
                const adattipus b03 = a0[0] * sor_i0.unsafe(j + 3) + a0[1] * sor_i1.unsafe(j + 3) + a0[2] * sor_i2.unsafe(j + 3) + a0[3] * sor_i3.unsafe(j + 3);
                const adattipus b10 = a1[0] * sor_i0.unsafe(j)     + a1[1] * sor_i1.unsafe(j)     + a1[2] * sor_i2.unsafe(j)     + a1[3] * sor_i3.unsafe(j);
                const adattipus b11 = a1[0] * sor_i0.unsafe(j + 1) + a1[1] * sor_i1.unsafe(j + 1) + a1[2] * sor_i2.unsafe(j + 1) + a1[3] * sor_i3.unsafe(j + 1);
                const adattipus b12 = a1[0] * sor_i0.unsafe(j + 2) + a1[1] * sor_i1.unsafe(j + 2) + a1[2] * sor_i2.unsafe(j + 2) + a1[3] * sor_i3.unsafe(j + 2);
                const adattipus b13 = a1[0] * sor_i0.unsafe(j + 3) + a1[1] * sor_i1.unsafe(j + 3) + a1[2] * sor_i2.unsafe(j + 3) + a1[3] * sor_i3.unsafe(j + 3);
                const adattipus b20 = a2[0] * sor_i0.unsafe(j)     + a2[1] * sor_i1.unsafe(j)     + a2[2] * sor_i2.unsafe(j)     + a2[3] * sor_i3.unsafe(j);
                const adattipus b21 = a2[0] * sor_i0.unsafe(j + 1) + a2[1] * sor_i1.unsafe(j + 1) + a2[2] * sor_i2.unsafe(j + 1) + a2[3] * sor_i3.unsafe(j + 1);
                const adattipus b22 = a2[0] * sor_i0.unsafe(j + 2) + a2[1] * sor_i1.unsafe(j + 2) + a2[2] * sor_i2.unsafe(j + 2) + a2[3] * sor_i3.unsafe(j + 2);
                const adattipus b23 = a2[0] * sor_i0.unsafe(j + 3) + a2[1] * sor_i1.unsafe(j + 3) + a2[2] * sor_i2.unsafe(j + 3) + a2[3] * sor_i3.unsafe(j + 3);
                const adattipus b30 = a3[0] * sor_i0.unsafe(j)     + a3[1] * sor_i1.unsafe(j)     + a3[2] * sor_i2.unsafe(j)     + a3[3] * sor_i3.unsafe(j);
                const adattipus b31 = a3[0] * sor_i0.unsafe(j + 1) + a3[1] * sor_i1.unsafe(j + 1) + a3[2] * sor_i2.unsafe(j + 1) + a3[3] * sor_i3.unsafe(j + 1);
                const adattipus b32 = a3[0] * sor_i0.unsafe(j + 2) + a3[1] * sor_i1.unsafe(j + 2) + a3[2] * sor_i2.unsafe(j + 2) + a3[3] * sor_i3.unsafe(j + 2);
                const adattipus b33 = a3[0] * sor_i0.unsafe(j + 3) + a3[1] * sor_i1.unsafe(j + 3) + a3[2] * sor_i2.unsafe(j + 3) + a3[3] * sor_i3.unsafe(j + 3);
                sor_i0.unsafe(j)     = b00;   sor_i0.unsafe(j + 1) = b01;    sor_i0.unsafe(j + 2) = b02;    sor_i0.unsafe(j + 3) = b03;
                sor_i1.unsafe(j)     = b10;   sor_i1.unsafe(j + 1) = b11;    sor_i1.unsafe(j + 2) = b12;    sor_i1.unsafe(j + 3) = b13;
                sor_i2.unsafe(j)     = b20;   sor_i2.unsafe(j + 1) = b21;    sor_i2.unsafe(j + 2) = b22;    sor_i2.unsafe(j + 3) = b23;
                sor_i3.unsafe(j)     = b30;   sor_i3.unsafe(j + 1) = b31;    sor_i3.unsafe(j + 2) = b32;    sor_i3.unsafe(j + 3) = b33;
            }
            for (j = hrow; j < nrow; j++) {
                const adattipus b00 = a0[0] * sor_i0.unsafe(j)     + a0[1] * sor_i1.unsafe(j)     + a0[2] * sor_i2.unsafe(j)     + a0[3] * sor_i3.unsafe(j);
                const adattipus b10 = a1[0] * sor_i0.unsafe(j)     + a1[1] * sor_i1.unsafe(j)     + a1[2] * sor_i2.unsafe(j)     + a1[3] * sor_i3.unsafe(j);
                const adattipus b20 = a2[0] * sor_i0.unsafe(j)     + a2[1] * sor_i1.unsafe(j)     + a2[2] * sor_i2.unsafe(j)     + a2[3] * sor_i3.unsafe(j);
                const adattipus b30 = a3[0] * sor_i0.unsafe(j)     + a3[1] * sor_i1.unsafe(j)     + a3[2] * sor_i2.unsafe(j)     + a3[3] * sor_i3.unsafe(j);
                sor_i0.unsafe(j)     = b00;
                sor_i1.unsafe(j)     = b10;
                sor_i2.unsafe(j)     = b20;
                sor_i3.unsafe(j)     = b30;
            }
        }
        for (meret_t i = hrow; i < nrow; i++) {
            vektor<adattipus> & sor_i = sorok.unsafe(i);
            adattipus oszto = abs(sor_i.unsafe(i)) < 1e-20f ? 1e20f : adattipus(1.0f) / sor_i.unsafe(i);
            meret_t j;
            for (j = 0; j < i; j++) {
                vektor<adattipus> & sor_j = sorok.unsafe(j);
                adattipus C = sor_j.unsafe(i) * oszto;
                meret_t k;
                for (k = 0; k < i; k++) sor_j.unsafe(k) -= C * sor_i.unsafe(k);
                sor_j.unsafe(k) = C; // nem negált invnél = -C;
                for (k++; k < row; k++) sor_j.unsafe(k) -= C * sor_i.unsafe(k);
            }
            for (j++; j < row; j++) {
                vektor<adattipus> & sor_j = sorok.unsafe(j);
                adattipus C = sor_j.unsafe(i) * oszto;
                meret_t k;
                for (k = 0; k < i; k++) sor_j.unsafe(k) -= C * sor_i.unsafe(k);
                sor_j.unsafe(k) = C; // nem negált invnél = -C;
                for (k++; k < row; k++) sor_j.unsafe(k) -= C * sor_i.unsafe(k);
            }
            for (j = 0; j < i; j++) sor_i.unsafe(j) *= oszto;
            sor_i.unsafe(j) = -oszto; // nem negált invnél = oszto;
            for (j++; j < row; j++) sor_i.unsafe(j) *= oszto;
        }
    }

    //****************************************************************
    void math_inv_p(bool neg) {
    // fõelemkiválasztással
    // Ez nem ráfektetésbiztos
    //****************************************************************
        egyforma_e_hiba(row, col, "matrix::math_inv_p");
        igaz_e_hiba(is_szimm, "matrix::math_ninv_p", "symmetrical matrix not allowed");

        const meret_t S2=row+row;
        meret_t *x = new meret_t[S2], *y = x + row;
        for (meret_t i = 0; i < S2; i++) x[i] = ~0;

        for (meret_t i = 0; i < row; i++) {

            // Pivot kiválasztása

            vektor<adattipus> & sor_i = sorok.unsafe(i);
            double diff = 0.0;
            meret_t V = ~0;
            meret_t j;
            for (j = 0; j < row; j++) if (x[j] == ~0) {
                double temp = abs(sor_i.unsafe(j));
                if (temp>diff) { diff = temp; V = j; }//v-edik oszlopot választjuk
            }
            if ((V == ~0) || (diff == 0))
                throw hiba("matrix::math_inv_p", "singular matrix");
            x[V] = i;
            y[i] = V;

            // Elemcsere

            adattipus A = adattipus(1.0) / sor_i.unsafe(V);

            for (j = 0; j < i; j++) {
                vektor<adattipus> & sor_j = sorok.unsafe(j);
                adattipus C = -sor_j.unsafe(V) * A;
                meret_t k;
                for (k = 0; k < V; k++)sor_j.unsafe(k) += C*sor_i.unsafe(k);
                sor_j.unsafe(k) = C;
                for (k++; k < row; k++)sor_j.unsafe(k) += C*sor_i.unsafe(k);
            }
            for (j++; j < row; j++) {
                vektor<adattipus> & sor_j = sorok.unsafe(j);
                adattipus C = -sor_j.unsafe(V) * A;
                meret_t k;
                for (k = 0; k < V; k++)sor_j.unsafe(k) += C*sor_i.unsafe(k);
                sor_j.unsafe(k) = C;
                for (k++; k < row; k++)sor_j.unsafe(k) += C*sor_i.unsafe(k);
            }
            for (j = 0;j<V;j++)sor_i.unsafe(j) *= A;
            sor_i.unsafe(j) = A;
            for (j++; j < row; j++)sor_i.unsafe(j) *= A;
        }

        // Sorok és oszlopok sorrendbe rakása

        adattipus * N = &t[0];
        meret_t i;
        for (i = 0; i < row; i++) {
            meret_t j;
            for (j = i; y[j] != i; j++)
                ;
            if (i != j) {
                vektor<adattipus> & sor_i = sorok.unsafe(i);
                vektor<adattipus> & sor_j = sorok.unsafe(j);
                for (meret_t k = 0; k < row; k++) { adattipus temp = sor_i.unsafe(k); sor_i.unsafe(k) = sor_j.unsafe(k); sor_j.unsafe(k) = temp; }
                y[j] = y[i];
            }
        }
        for (i = 0;i < row; i++) {
            meret_t j;
            for (j = i; x[j] != i; j++)
                ;
            if (i != j) {
                for (meret_t k = 0; k < row; k++) { 
                    vektor<adattipus> & sor_k = sorok.unsafe(k); 
                    adattipus temp  = sor_k.unsafe(i);
                    sor_k.unsafe(i) = sor_k.unsafe(j);
                    sor_k.unsafe(j) = temp; 
                }
                x[j] = x[i];
            }
        }

        delete[] x;
        if (neg) math_neg_nembiztos();
    }

    //***********************************************************************
    void zero_nembiztos() {
    // Ez nem ráfektetésbiztos
    //***********************************************************************
        t.zero();
    }

    //***********************************************************************
    void egyseg() {
    //***********************************************************************
        igaz_e_hiba(is_szimm, "matrix::egyseg", "symmetrical matrix not allowed");
        t.zero();
        for (meret_t i = 0; i < row; i++)
            sorok.unsafe(i).unsafe(i) = 1.0;
    }

    //***********************************************************************
    void math_1_ninv_mul(matrix & yb, const matrix & xa) {
    // yb 1x1-es, ezt ninvertálja és zb*xa-t magába teszi
    //***********************************************************************
        const adattipus nzb = adattipus(-1) / yb.sorok.unsafe(0).unsafe(0);
        yb.sorok.unsafe(0).unsafe(0) = nzb;
        for (uns i = 0; i < col; i++)
            sorok.unsafe(0).unsafe(i) = nzb*xa.sorok.unsafe(0).unsafe(i);
    }

    //***********************************************************************
    void math_2_ninv_mul(matrix & yb, const matrix & xa) {
    // yb 2x2-es, ezt ninvertálja és zb*xa-t magába teszi
    //***********************************************************************
        igaz_e_hiba(yb.is_szimm, "matrix::math_2_ninv_mul", "symmetrical matrix not allowed");
        adattipus * be0 = &yb.sorok.unsafe(0).unsafe(0);
        adattipus * be1 = &yb.sorok.unsafe(1).unsafe(0);
        const adattipus p0 = abs(be0[0]) < 1e-20f ? 1e20f : adattipus(1.0f) / be0[0];
        const adattipus p2 = -be1[0] * p0;
        const adattipus p1 = be0[1] * p0;
        const adattipus p3 = be1[1] + p2 * be0[1];
        const adattipus oszto2 = abs(p3) < 1e-20f ? 1e20f : adattipus(1.0f) / p3;
        const adattipus C2 = -p1 * oszto2;
        const adattipus nzb00 = be0[0] = -(p0 + C2 * p2);
        const adattipus nzb01 = be0[1] = -C2;
        const adattipus nzb10 = be1[0] = -p2 * oszto2;
        const adattipus nzb11 = be1[1] = -oszto2;

        vektor<adattipus> & nzbxa0 = sorok.unsafe(0);
        vektor<adattipus> & nzbxa1 = sorok.unsafe(1);
        const vektor<adattipus> & xa0 = xa.sorok.unsafe(0);
        const vektor<adattipus> & xa1 = xa.sorok.unsafe(1);
        for (uns i = 0; i < col; i++) {
            nzbxa0.unsafe(i) = nzb00 * xa0.unsafe(i) + nzb01 * xa1.unsafe(i);
            nzbxa1.unsafe(i) = nzb10 * xa0.unsafe(i) + nzb11 * xa1.unsafe(i);
        }
    }

    //***********************************************************************
    void math_2_ninv_mul_symm(matrix & yb, const matrix & xa) {
    // yb 2x2-es szimmetrikus teljesként tárolt, ezt ninvertálja és zb*xa-t magába teszi
    //***********************************************************************
        igaz_e_hiba(yb.is_szimm, "matrix::math_2_ninv_mul", "symmetrical matrix not allowed");
        adattipus * be0 = &yb.sorok.unsafe(0).unsafe(0);
        adattipus * be1 = &yb.sorok.unsafe(1).unsafe(0);
        const adattipus p0 = abs(be0[0]) < 1e-20f ? 1e20f : adattipus(1.0f) / be0[0];
        const adattipus p2 = -be0[1] * p0;
        const adattipus p3 = be1[1] + p2 * be0[1];
        const adattipus oszto2 = abs(p3) < 1e-20f ? 1e20f : adattipus(1.0f) / p3;
        const adattipus C2 = -p2 * oszto2;
        const adattipus nzb00 = be0[0] = C2 * p2 - p0;
        const adattipus nzb01 = be0[1] = be1[0] = C2;
        const adattipus nzb11 = be1[1] = -oszto2;

        vektor<adattipus> & nzbxa0 = sorok.unsafe(0);
        vektor<adattipus> & nzbxa1 = sorok.unsafe(1);
        const vektor<adattipus> & xa0 = xa.sorok.unsafe(0);
        const vektor<adattipus> & xa1 = xa.sorok.unsafe(1);
        for (uns i = 0; i < col; i++) {
            nzbxa0.unsafe(i) = nzb00 * xa0.unsafe(i) + nzb01 * xa1.unsafe(i);
            nzbxa1.unsafe(i) = nzb01 * xa0.unsafe(i) + nzb11 * xa1.unsafe(i);
        }
    }

    //***********************************************************************
    void math_1_add_mul(const matrix & ya, const matrix & xb, const matrix & nzbxa) {
    // xb 1 oszlop, nzbxa 1 sor
    //***********************************************************************
        igaz_e_hiba(ya.is_szimm || is_szimm, "matrix::math_1_add_mul", "symmetrical matrix not allowed");
        for (uns i = 0; i < row; i++) {
            const adattipus xbe = xb.sorok.unsafe(i).unsafe(0);
            for (uns j = 0; j < col; j++)
                sorok.unsafe(i).unsafe(j) = ya.sorok.unsafe(i).unsafe(j) + xbe * nzbxa.sorok.unsafe(0).unsafe(j);
        }
    }

    //***********************************************************************
    void math_1_add_mul_symm(const matrix & ya, const matrix & xb, const matrix & nzbxa) {
    // xb 1 oszlop, nzbxa 1 sor
    //***********************************************************************
        igaz_e_hiba(!ya.is_szimm || !is_szimm, "matrix::math_1_add_mul_symm", "nonsymmetrical matrix not allowed");
        for (uns i = 0; i < row; i++) {
            const adattipus xbe = xb.sorok.unsafe(i).unsafe(0);
            for (uns j = i; j < col; j++)
                sorok.unsafe(i).unsafe(j) = ya.sorok.unsafe(i).unsafe(j) + xbe * nzbxa.sorok.unsafe(0).unsafe(j);
        }
    }

    //***********************************************************************
    void math_2_add_mul(const matrix & ya, const matrix & xb, const matrix & nzbxa) {
    // xb 2 oszlop, nzbxa 2 sor
    //***********************************************************************
        igaz_e_hiba(ya.is_szimm || is_szimm, "matrix::math_2_add_mul", "symmetrical matrix not allowed");
        const vektor<adattipus> & nzbxa0 = nzbxa.sorok.unsafe(0);
        const vektor<adattipus> & nzbxa1 = nzbxa.sorok.unsafe(1);
        for (uns i = 0; i < row; i++) {
            const adattipus xb0 = xb.sorok.unsafe(i).unsafe(0);
            const adattipus xb1 = xb.sorok.unsafe(i).unsafe(1);
            for (uns j = 0; j < col; j++)
                sorok.unsafe(i).unsafe(j) = ya.sorok.unsafe(i).unsafe(j) + xb0 * nzbxa0.unsafe(j) + xb1 * nzbxa1.unsafe(j);
        }
    }

    //***********************************************************************
    void math_2_add_mul_symm(const matrix & ya, const matrix & xb, const matrix & nzbxa) {
    // xb 2 oszlop, nzbxa 2 sor
    //***********************************************************************
        igaz_e_hiba(!ya.is_szimm || !is_szimm, "matrix::math_2_add_mul_symm", "symmetrical matrix not allowed");
        const vektor<adattipus> & nzbxa0 = nzbxa.sorok.unsafe(0);
        const vektor<adattipus> & nzbxa1 = nzbxa.sorok.unsafe(1);
        for (uns i = 0; i < row; i++) {
            const adattipus xb0 = xb.sorok.unsafe(i).unsafe(0);
            const adattipus xb1 = xb.sorok.unsafe(i).unsafe(1);
            for (uns j = i; j < col; j++)
                sorok.unsafe(i).unsafe(j) = ya.sorok.unsafe(i).unsafe(j) + xb0 * nzbxa0.unsafe(j) + xb1 * nzbxa1.unsafe(j);
        }
    }

    //***********************************************************************
    friend inline void math_1x1_mul(vektor<adattipus> & nzbjb, const matrix & nzb, const vektor<adattipus> & jb) {
    // minden 1-es méretû
    //***********************************************************************
        nzbjb.unsafe(0) = nzb.sorok.unsafe(0).unsafe(0) * jb.unsafe(0);
    }

    //***********************************************************************
    friend inline void math_2x2_mul(vektor<adattipus> & nzbjb, const matrix & nzb, const vektor<adattipus> & jb) {
    // nzb 2x2-es méretû
    //***********************************************************************
        igaz_e_hiba(nzb.is_szimm, "math_2x2_mul", "symmetrical matrix not allowed");
        nzbjb.unsafe(0) = nzb.sorok.unsafe(0).unsafe(0) * jb.unsafe(0) + nzb.sorok.unsafe(0).unsafe(1) * jb.unsafe(1);
        nzbjb.unsafe(1) = nzb.sorok.unsafe(1).unsafe(0) * jb.unsafe(0) + nzb.sorok.unsafe(1).unsafe(1) * jb.unsafe(1);
    }

    //***********************************************************************
    friend inline void math_1_add_mul_jred(vektor<adattipus> & jred, const vektor<adattipus> & ja, const matrix & xb, const vektor<adattipus> & nzbjb) {
    // nzbjb 1-es
    //***********************************************************************
        const adattipus a = nzbjb.unsafe(0);
        for (meret_t i = 0; i < jred.size(); i++)
            jred.unsafe(i) = ja.unsafe(i) + xb.sorok.unsafe(0).unsafe(i) * a;
    }

    //***********************************************************************
    friend inline void math_2_add_mul_jred(vektor<adattipus> & jred, const vektor<adattipus> & ja, const matrix & xb, const vektor<adattipus> & nzbjb) {
    // nzbjb 2-es
    //***********************************************************************
        const adattipus a0 = nzbjb.unsafe(0);
        const adattipus a1 = nzbjb.unsafe(1);
        meret_t j = 0;
        for (meret_t i = 0; i < jred.size(); i++, j+=2)
            jred.unsafe(i) = ja.unsafe(i) + xb.sorok.unsafe(0).unsafe(j) * a0 + xb.sorok.unsafe(0).unsafe(j + 1) * a1;
    }

    //***********************************************************************
    friend inline void math_1_add_mul_ub(vektor<adattipus> & ub, const vektor<adattipus> & nzbjb, const matrix & nzbxa, const vektor<adattipus> & UA) {
    // ub 1-es
    //***********************************************************************
        ub.unsafe(0) = nzbjb.unsafe(0) + math_mul(nzbxa.sorok.unsafe(0), UA);
    }

    //***********************************************************************
    friend inline void math_2_add_mul_ub(vektor<adattipus> & ub, const vektor<adattipus> & nzbjb, const matrix & nzbxa, const vektor<adattipus> & UA) {
    // ub 2-es
    //***********************************************************************
        ub.unsafe(0) = nzbjb.unsafe(0) + math_mul(nzbxa.sorok.unsafe(0), UA);
        ub.unsafe(1) = nzbjb.unsafe(1) + math_mul(nzbxa.sorok.unsafe(1), UA);
    }
    
    //***********************************************************************
    void math_symm_ninv_of_nonsymm(){
    // A szimmetrikus mátrixot nemszimmetrikus formában tároljuk, ráfektetett is lehet
    //***********************************************************************
        igaz_e_hiba(is_szimm, "math_symm_ninv_of_nonsymm_rafektetett", "only nonsymmetrically stored symmetrical matrix allowed");
#define NEWKIVALTO 32

        if (row == 1) {
            sorok.unsafe(0).unsafe(0) = adattipus(-1) / sorok.unsafe(0).unsafe(0);
            return;
        }
        if (row == 2) {
            adattipus * const ps0 = &sorok.unsafe(0).unsafe(0);
            adattipus * const ps1 = &sorok.unsafe(1).unsafe(0);
            adattipus p0 = abs(ps0[0]) < 1e-20f ? 1e20f : adattipus(1.0f) / ps0[0];
            adattipus p2 = -ps0[1] * p0;
            adattipus p3 = ps1[1] + p2 * ps0[1];
            adattipus oszto2 = abs(p3) < 1e-20f ? 1e20f : adattipus(1.0f) / p3;
            adattipus C2 = -p2 * oszto2;
            ps0[0] = C2 * p2 - p0;
            ps1[0] = ps0[1] = C2;
            ps1[1] = -oszto2;
            return;
        }

        adattipus dum[2 * NEWKIVALTO], *b0 = row > NEWKIVALTO ? new adattipus[2 * row] : dum, *bs = b0 + row;
        uns i, j, k;

        for (i = 0; i < row; i++) {
            vektor<adattipus> & sor_i = sorok.unsafe(i);
            adattipus & pivot = sor_i.unsafe(i);
            adattipus oszto = abs(pivot) < 1e-20f ? 1e20f : adattipus(1.0f) / pivot;
            for (k = 0; k < i; ++k) { adattipus & v = sorok.unsafe(k).unsafe(i); b0[k] = v;	    bs[k] = v *= oszto; }
                                                                                 b0[k] = pivot;	bs[k] = pivot = -oszto;
            for (k++; k < row; k++) { adattipus & v = sor_i.unsafe(k);           b0[k] = v;	    bs[k] = v *= oszto; }

            for (j = 0; j < i; j++) {//j=0-tól i-1-ig
                vektor<adattipus> & sor_j = sorok.unsafe(j);
                const adattipus x = b0[j]; // i. sor j. eleme
                const adattipus x2 = sor_j.unsafe(i); // j. sor i. eleme
                for (k = j; k < row; k++) sor_j.unsafe(k) -= x*bs[k];
                sor_j.unsafe(i) = x2;
            }
            for (j++; j < row; j++) {//j=i+1-tõl S-1-ig
                vektor<adattipus> & sor_j = sorok.unsafe(j);
                const adattipus x = b0[j];
                for (k = j; k < row; k++) sor_j.unsafe(k) -= x*bs[k];//k!=i mindig igaz, mert j=i+1-tõl indul a ciklus
            }
        }
        if (row > NEWKIVALTO)delete[] b0;
        symmetrize_from_upper();
    }

    //***********************************************************************
    void symmetrize_from_upper() {
    // Egy nemszimmetrikus mátrix felsõ háromszögét az alsóba másolja
    //***********************************************************************
        igaz_e_hiba(is_szimm, "symmetrize_from_upper", "only nonsymmetrically stored symmetrical matrix allowed");
        egyforma_e_hiba(row, col, "matrix::symmetrize_from_upper");
        for (uns i = 1; i < col; i++)
            for (uns j = 0; j < i; j++) {
                sorok.unsafe(i).unsafe(j) = sorok.unsafe(j).unsafe(i);
            }
    }

    //***********************************************************************
    void math_add_mul_t_symm(const matrix & ya, const matrix & xb, const matrix & nzbxat){
    //***********************************************************************
        igaz_e_hiba(!is_szimm || !ya.is_szimm, "matrix::math_add_mul_t_symm", "symmetrical matrix required");
        egyforma_e_hiba(ya.col, col, "math_add_mul_t_symm col");
        egyforma_e_hiba(xb.row, row, "math_add_mul_t_symm row row");
        egyforma_e_hiba(nzbxat.row, col, "math_add_mul_t_symm row col");
        egyforma_e_hiba(xb.col, nzbxat.col, "math_add_mul_t_symm col col");

        cuns ni = row, nj = col, nk = xb.col;
        cuns di = row % 4, dj = col % 4, dk = xb.col % 4;
        cuns hi = ni - di, hj = nj - dj, hk = nk - dk;
        for (uns i = 0; i < hi; i += 4) {
            const vektor<adattipus> & xb_sor_i0 = xb.sorok.unsafe(i + 0);
            const vektor<adattipus> & xb_sor_i1 = xb.sorok.unsafe(i + 1);
            const vektor<adattipus> & xb_sor_i2 = xb.sorok.unsafe(i + 2);
            const vektor<adattipus> & xb_sor_i3 = xb.sorok.unsafe(i + 3);
            const vektor<adattipus> & nzbxat_sor_i0 = nzbxat.sorok.unsafe(i + 0);
            const vektor<adattipus> & nzbxat_sor_i1 = nzbxat.sorok.unsafe(i + 1);
            const vektor<adattipus> & nzbxat_sor_i2 = nzbxat.sorok.unsafe(i + 2);
            const vektor<adattipus> & nzbxat_sor_i3 = nzbxat.sorok.unsafe(i + 3);
            adattipus haromszog[10] = { adattipus() };
            for (uns k = 0; k < hk; k += 4) {
                haromszog[0] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_i0.unsafe(k + 0) 
                              + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_i0.unsafe(k + 1) 
                              + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_i0.unsafe(k + 2) 
                              + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_i0.unsafe(k + 3);
                haromszog[1] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_i1.unsafe(k + 0) 
                              + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_i1.unsafe(k + 1) 
                              + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_i1.unsafe(k + 2) 
                              + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_i1.unsafe(k + 3);
                haromszog[2] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_i2.unsafe(k + 0) 
                              + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_i2.unsafe(k + 1) 
                              + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_i2.unsafe(k + 2) 
                              + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_i2.unsafe(k + 3);
                haromszog[3] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_i3.unsafe(k + 0) 
                              + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_i3.unsafe(k + 1) 
                              + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_i3.unsafe(k + 2) 
                              + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_i3.unsafe(k + 3);

                haromszog[4] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_i1.unsafe(k + 0) 
                              + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_i1.unsafe(k + 1) 
                              + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_i1.unsafe(k + 2) 
                              + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_i1.unsafe(k + 3);
                haromszog[5] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_i2.unsafe(k + 0) 
                              + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_i2.unsafe(k + 1) 
                              + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_i2.unsafe(k + 2) 
                              + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_i2.unsafe(k + 3);
                haromszog[6] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_i3.unsafe(k + 0) 
                              + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_i3.unsafe(k + 1) 
                              + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_i3.unsafe(k + 2) 
                              + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_i3.unsafe(k + 3);

                haromszog[7] += xb_sor_i2.unsafe(k + 0) * nzbxat_sor_i2.unsafe(k + 0) 
                              + xb_sor_i2.unsafe(k + 1) * nzbxat_sor_i2.unsafe(k + 1) 
                              + xb_sor_i2.unsafe(k + 2) * nzbxat_sor_i2.unsafe(k + 2) 
                              + xb_sor_i2.unsafe(k + 3) * nzbxat_sor_i2.unsafe(k + 3);
                haromszog[8] += xb_sor_i2.unsafe(k + 0) * nzbxat_sor_i3.unsafe(k + 0) 
                              + xb_sor_i2.unsafe(k + 1) * nzbxat_sor_i3.unsafe(k + 1) 
                              + xb_sor_i2.unsafe(k + 2) * nzbxat_sor_i3.unsafe(k + 2) 
                              + xb_sor_i2.unsafe(k + 3) * nzbxat_sor_i3.unsafe(k + 3);

                haromszog[9] += xb_sor_i3.unsafe(k + 0) * nzbxat_sor_i3.unsafe(k + 0) 
                              + xb_sor_i3.unsafe(k + 1) * nzbxat_sor_i3.unsafe(k + 1) 
                              + xb_sor_i3.unsafe(k + 2) * nzbxat_sor_i3.unsafe(k + 2) 
                              + xb_sor_i3.unsafe(k + 3) * nzbxat_sor_i3.unsafe(k + 3);
            }
            for (uns k = hk; k < nk; k++) {
                haromszog[0] += xb_sor_i0.unsafe(k) * nzbxat_sor_i0.unsafe(k);
                haromszog[1] += xb_sor_i0.unsafe(k) * nzbxat_sor_i1.unsafe(k);
                haromszog[2] += xb_sor_i0.unsafe(k) * nzbxat_sor_i2.unsafe(k);
                haromszog[3] += xb_sor_i0.unsafe(k) * nzbxat_sor_i3.unsafe(k);

                haromszog[4] += xb_sor_i1.unsafe(k) * nzbxat_sor_i1.unsafe(k);
                haromszog[5] += xb_sor_i1.unsafe(k) * nzbxat_sor_i2.unsafe(k);
                haromszog[6] += xb_sor_i1.unsafe(k) * nzbxat_sor_i3.unsafe(k);

                haromszog[7] += xb_sor_i2.unsafe(k) * nzbxat_sor_i2.unsafe(k);
                haromszog[8] += xb_sor_i2.unsafe(k) * nzbxat_sor_i3.unsafe(k);

                haromszog[9] += xb_sor_i3.unsafe(k) * nzbxat_sor_i3.unsafe(k);

            }

            vektor<adattipus> & sor_i0 = sorok.unsafe(i + 0);
            vektor<adattipus> & sor_i1 = sorok.unsafe(i + 1);
            vektor<adattipus> & sor_i2 = sorok.unsafe(i + 2);
            vektor<adattipus> & sor_i3 = sorok.unsafe(i + 3);
            const vektor<adattipus> & ya_sor_i0 = ya.sorok.unsafe(i + 0);
            const vektor<adattipus> & ya_sor_i1 = ya.sorok.unsafe(i + 1);
            const vektor<adattipus> & ya_sor_i2 = ya.sorok.unsafe(i + 2);
            const vektor<adattipus> & ya_sor_i3 = ya.sorok.unsafe(i + 3);

            sor_i0.unsafe(i + 0) = haromszog[0] + ya_sor_i0.unsafe(i + 0);
            sor_i0.unsafe(i + 1) = haromszog[1] + ya_sor_i0.unsafe(i + 1);
            sor_i0.unsafe(i + 2) = haromszog[2] + ya_sor_i0.unsafe(i + 2);
            sor_i0.unsafe(i + 3) = haromszog[3] + ya_sor_i0.unsafe(i + 3);
            
            sor_i1.unsafe(i + 1) = haromszog[4] + ya_sor_i1.unsafe(i + 1);
            sor_i1.unsafe(i + 2) = haromszog[5] + ya_sor_i1.unsafe(i + 2);
            sor_i1.unsafe(i + 3) = haromszog[6] + ya_sor_i1.unsafe(i + 3);

            sor_i2.unsafe(i + 2) = haromszog[7] + ya_sor_i2.unsafe(i + 2);
            sor_i2.unsafe(i + 3) = haromszog[8] + ya_sor_i2.unsafe(i + 3);

            sor_i3.unsafe(i + 3) = haromszog[9] + ya_sor_i3.unsafe(i + 3);

            for (uns j = i + 4; j < hj; j += 4) {
                const vektor<adattipus> & nzbxat_sor_j0 = nzbxat.sorok.unsafe(j + 0);
                const vektor<adattipus> & nzbxat_sor_j1 = nzbxat.sorok.unsafe(j + 1);
                const vektor<adattipus> & nzbxat_sor_j2 = nzbxat.sorok.unsafe(j + 2);
                const vektor<adattipus> & nzbxat_sor_j3 = nzbxat.sorok.unsafe(j + 3);
                adattipus negyzet[16] = { adattipus() };
                for (uns k = 0; k < hk; k += 4) {
                    negyzet[0] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_j0.unsafe(k + 0) 
                                + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_j0.unsafe(k + 1) 
                                + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_j0.unsafe(k + 2) 
                                + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_j0.unsafe(k + 3);
                    negyzet[1] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_j1.unsafe(k + 0) 
                                + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_j1.unsafe(k + 1) 
                                + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_j1.unsafe(k + 2) 
                                + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_j1.unsafe(k + 3);
                    negyzet[2] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_j2.unsafe(k + 0) 
                                + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_j2.unsafe(k + 1) 
                                + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_j2.unsafe(k + 2) 
                                + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_j2.unsafe(k + 3);
                    negyzet[3] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_j3.unsafe(k + 0) 
                                + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_j3.unsafe(k + 1) 
                                + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_j3.unsafe(k + 2) 
                                + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_j3.unsafe(k + 3);

                    negyzet[4] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_j0.unsafe(k + 0) 
                                + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_j0.unsafe(k + 1)
                                + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_j0.unsafe(k + 2)
                                + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_j0.unsafe(k + 3);
                    negyzet[5] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_j1.unsafe(k + 0)
                                + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_j1.unsafe(k + 1)
                                + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_j1.unsafe(k + 2)
                                + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_j1.unsafe(k + 3);
                    negyzet[6] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_j2.unsafe(k + 0)
                                + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_j2.unsafe(k + 1)
                                + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_j2.unsafe(k + 2)
                                + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_j2.unsafe(k + 3);
                    negyzet[7] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_j3.unsafe(k + 0)
                                + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_j3.unsafe(k + 1)
                                + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_j3.unsafe(k + 2)
                                + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_j3.unsafe(k + 3);

                    negyzet[8] += xb_sor_i2.unsafe(k + 0) * nzbxat_sor_j0.unsafe(k + 0) 
                                + xb_sor_i2.unsafe(k + 1) * nzbxat_sor_j0.unsafe(k + 1)
                                + xb_sor_i2.unsafe(k + 2) * nzbxat_sor_j0.unsafe(k + 2)
                                + xb_sor_i2.unsafe(k + 3) * nzbxat_sor_j0.unsafe(k + 3);
                    negyzet[9] += xb_sor_i2.unsafe(k + 0) * nzbxat_sor_j1.unsafe(k + 0)
                                + xb_sor_i2.unsafe(k + 1) * nzbxat_sor_j1.unsafe(k + 1)
                                + xb_sor_i2.unsafe(k + 2) * nzbxat_sor_j1.unsafe(k + 2)
                                + xb_sor_i2.unsafe(k + 3) * nzbxat_sor_j1.unsafe(k + 3);
                    negyzet[10]+= xb_sor_i2.unsafe(k + 0) * nzbxat_sor_j2.unsafe(k + 0)
                                + xb_sor_i2.unsafe(k + 1) * nzbxat_sor_j2.unsafe(k + 1)
                                + xb_sor_i2.unsafe(k + 2) * nzbxat_sor_j2.unsafe(k + 2)
                                + xb_sor_i2.unsafe(k + 3) * nzbxat_sor_j2.unsafe(k + 3);
                    negyzet[11]+= xb_sor_i2.unsafe(k + 0) * nzbxat_sor_j3.unsafe(k + 0)
                                + xb_sor_i2.unsafe(k + 1) * nzbxat_sor_j3.unsafe(k + 1)
                                + xb_sor_i2.unsafe(k + 2) * nzbxat_sor_j3.unsafe(k + 2)
                                + xb_sor_i2.unsafe(k + 3) * nzbxat_sor_j3.unsafe(k + 3);

                    negyzet[12]+= xb_sor_i3.unsafe(k + 0) * nzbxat_sor_j0.unsafe(k + 0) 
                                + xb_sor_i3.unsafe(k + 1) * nzbxat_sor_j0.unsafe(k + 1)
                                + xb_sor_i3.unsafe(k + 2) * nzbxat_sor_j0.unsafe(k + 2)
                                + xb_sor_i3.unsafe(k + 3) * nzbxat_sor_j0.unsafe(k + 3);
                    negyzet[13]+= xb_sor_i3.unsafe(k + 0) * nzbxat_sor_j1.unsafe(k + 0)
                                + xb_sor_i3.unsafe(k + 1) * nzbxat_sor_j1.unsafe(k + 1)
                                + xb_sor_i3.unsafe(k + 2) * nzbxat_sor_j1.unsafe(k + 2)
                                + xb_sor_i3.unsafe(k + 3) * nzbxat_sor_j1.unsafe(k + 3);
                    negyzet[14]+= xb_sor_i3.unsafe(k + 0) * nzbxat_sor_j2.unsafe(k + 0)
                                + xb_sor_i3.unsafe(k + 1) * nzbxat_sor_j2.unsafe(k + 1)
                                + xb_sor_i3.unsafe(k + 2) * nzbxat_sor_j2.unsafe(k + 2)
                                + xb_sor_i3.unsafe(k + 3) * nzbxat_sor_j2.unsafe(k + 3);
                    negyzet[15]+= xb_sor_i3.unsafe(k + 0) * nzbxat_sor_j3.unsafe(k + 0)
                                + xb_sor_i3.unsafe(k + 1) * nzbxat_sor_j3.unsafe(k + 1)
                                + xb_sor_i3.unsafe(k + 2) * nzbxat_sor_j3.unsafe(k + 2)
                                + xb_sor_i3.unsafe(k + 3) * nzbxat_sor_j3.unsafe(k + 3);
                }
                for (uns k = hk; k < nk; k++) {
                    negyzet[0]  += xb_sor_i0.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyzet[1]  += xb_sor_i0.unsafe(k) * nzbxat_sor_j1.unsafe(k);
                    negyzet[2]  += xb_sor_i0.unsafe(k) * nzbxat_sor_j2.unsafe(k);
                    negyzet[3]  += xb_sor_i0.unsafe(k) * nzbxat_sor_j3.unsafe(k);

                    negyzet[4]  += xb_sor_i1.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyzet[5]  += xb_sor_i1.unsafe(k) * nzbxat_sor_j1.unsafe(k);
                    negyzet[6]  += xb_sor_i1.unsafe(k) * nzbxat_sor_j2.unsafe(k);
                    negyzet[7]  += xb_sor_i1.unsafe(k) * nzbxat_sor_j3.unsafe(k);

                    negyzet[8]  += xb_sor_i2.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyzet[9]  += xb_sor_i2.unsafe(k) * nzbxat_sor_j1.unsafe(k);
                    negyzet[10] += xb_sor_i2.unsafe(k) * nzbxat_sor_j2.unsafe(k);
                    negyzet[11] += xb_sor_i2.unsafe(k) * nzbxat_sor_j3.unsafe(k);

                    negyzet[12] += xb_sor_i3.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyzet[13] += xb_sor_i3.unsafe(k) * nzbxat_sor_j1.unsafe(k);
                    negyzet[14] += xb_sor_i3.unsafe(k) * nzbxat_sor_j2.unsafe(k);
                    negyzet[15] += xb_sor_i3.unsafe(k) * nzbxat_sor_j3.unsafe(k);
                }

                sor_i0.unsafe(j + 0) = negyzet[0]  + ya_sor_i0.unsafe(j + 0);
                sor_i0.unsafe(j + 1) = negyzet[1]  + ya_sor_i0.unsafe(j + 1);
                sor_i0.unsafe(j + 2) = negyzet[2]  + ya_sor_i0.unsafe(j + 2);
                sor_i0.unsafe(j + 3) = negyzet[3]  + ya_sor_i0.unsafe(j + 3);

                sor_i1.unsafe(j + 0) = negyzet[4]  + ya_sor_i1.unsafe(j + 0);
                sor_i1.unsafe(j + 1) = negyzet[5]  + ya_sor_i1.unsafe(j + 1);
                sor_i1.unsafe(j + 2) = negyzet[6]  + ya_sor_i1.unsafe(j + 2);
                sor_i1.unsafe(j + 3) = negyzet[7]  + ya_sor_i1.unsafe(j + 3);

                sor_i2.unsafe(j + 0) = negyzet[8]  + ya_sor_i2.unsafe(j + 0);
                sor_i2.unsafe(j + 1) = negyzet[9]  + ya_sor_i2.unsafe(j + 1);
                sor_i2.unsafe(j + 2) = negyzet[10] + ya_sor_i2.unsafe(j + 2);
                sor_i2.unsafe(j + 3) = negyzet[11] + ya_sor_i2.unsafe(j + 3);

                sor_i3.unsafe(j + 0) = negyzet[12] + ya_sor_i3.unsafe(j + 0);
                sor_i3.unsafe(j + 1) = negyzet[13] + ya_sor_i3.unsafe(j + 1);
                sor_i3.unsafe(j + 2) = negyzet[14] + ya_sor_i3.unsafe(j + 2);
                sor_i3.unsafe(j + 3) = negyzet[15] + ya_sor_i3.unsafe(j + 3);
            }
            for (uns j = hj; j < nj; j++) {
                const vektor<adattipus> & nzbxat_sor_j0 = nzbxat.sorok.unsafe(j);
                adattipus negyes[4] = { adattipus() };
                for (uns k = 0; k < nk; k++) {
                    negyes[0] += xb_sor_i0.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyes[1] += xb_sor_i1.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyes[2] += xb_sor_i2.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyes[3] += xb_sor_i3.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                }
                sor_i0.unsafe(j) = negyes[0] + ya_sor_i0.unsafe(j);
                sor_i1.unsafe(j) = negyes[1] + ya_sor_i1.unsafe(j);
                sor_i2.unsafe(j) = negyes[2] + ya_sor_i2.unsafe(j);
                sor_i3.unsafe(j) = negyes[3] + ya_sor_i3.unsafe(j);
            }
	    }
        for (uns i = hi; i < ni; i++) {
            const vektor<adattipus> & xb_sor_i = xb.sorok.unsafe(i);
            for (uns j = i; j < nj; j++) {
                const vektor<adattipus> & nzbxat_sor_j = nzbxat.sorok.unsafe(j);
                adattipus sum = adattipus();
                for (uns k = 0; k < nk; k++) {
                    sum += xb_sor_i.unsafe(k) * nzbxat_sor_j.unsafe(k);
                }
                sorok.unsafe(i).unsafe(j) = sum + ya.sorok.unsafe(i).unsafe(j);
            }
        }
    }

    //***********************************************************************
    void math_sub_mul_t_symm_in_nonsymm(const matrix & c, const matrix & a, const matrix & b, bool is_symmetrize_needed){
    //***********************************************************************
        igaz_e_hiba(is_szimm || c.is_szimm, "matrix::math_sub_mul_t_symm_in_nonsymm", "nonsymmetrical matrix required");
        egyforma_e_hiba(c.col, col, "math_sub_mul_t_symm_in_nonsymm col");
        egyforma_e_hiba(a.row, row, "math_sub_mul_t_symm_in_nonsymm row row");
        egyforma_e_hiba(b.row, col, "math_sub_mul_t_symm_in_nonsymm row col");
        egyforma_e_hiba(a.col, b.col, "math_sub_mul_t_symm_in_nonsymm col col");

        cuns ni = row, nj = col, nk = a.col;
        cuns di = row % 4, dj = col % 4, dk = a.col % 4;
        cuns hi = ni - di, hj = nj - dj, hk = nk - dk;
        for (uns i = 0; i < hi; i += 4) {
            const vektor<adattipus> & xb_sor_i0 = a.sorok.unsafe(i + 0);
            const vektor<adattipus> & xb_sor_i1 = a.sorok.unsafe(i + 1);
            const vektor<adattipus> & xb_sor_i2 = a.sorok.unsafe(i + 2);
            const vektor<adattipus> & xb_sor_i3 = a.sorok.unsafe(i + 3);
            const vektor<adattipus> & nzbxat_sor_i0 = b.sorok.unsafe(i + 0);
            const vektor<adattipus> & nzbxat_sor_i1 = b.sorok.unsafe(i + 1);
            const vektor<adattipus> & nzbxat_sor_i2 = b.sorok.unsafe(i + 2);
            const vektor<adattipus> & nzbxat_sor_i3 = b.sorok.unsafe(i + 3);
            adattipus haromszog[10] = { adattipus() };
            for (uns k = 0; k < hk; k += 4) {
                haromszog[0] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_i0.unsafe(k + 0) 
                              + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_i0.unsafe(k + 1) 
                              + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_i0.unsafe(k + 2) 
                              + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_i0.unsafe(k + 3);
                haromszog[1] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_i1.unsafe(k + 0) 
                              + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_i1.unsafe(k + 1) 
                              + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_i1.unsafe(k + 2) 
                              + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_i1.unsafe(k + 3);
                haromszog[2] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_i2.unsafe(k + 0) 
                              + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_i2.unsafe(k + 1) 
                              + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_i2.unsafe(k + 2) 
                              + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_i2.unsafe(k + 3);
                haromszog[3] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_i3.unsafe(k + 0) 
                              + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_i3.unsafe(k + 1) 
                              + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_i3.unsafe(k + 2) 
                              + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_i3.unsafe(k + 3);

                haromszog[4] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_i1.unsafe(k + 0) 
                              + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_i1.unsafe(k + 1) 
                              + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_i1.unsafe(k + 2) 
                              + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_i1.unsafe(k + 3);
                haromszog[5] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_i2.unsafe(k + 0) 
                              + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_i2.unsafe(k + 1) 
                              + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_i2.unsafe(k + 2) 
                              + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_i2.unsafe(k + 3);
                haromszog[6] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_i3.unsafe(k + 0) 
                              + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_i3.unsafe(k + 1) 
                              + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_i3.unsafe(k + 2) 
                              + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_i3.unsafe(k + 3);

                haromszog[7] += xb_sor_i2.unsafe(k + 0) * nzbxat_sor_i2.unsafe(k + 0) 
                              + xb_sor_i2.unsafe(k + 1) * nzbxat_sor_i2.unsafe(k + 1) 
                              + xb_sor_i2.unsafe(k + 2) * nzbxat_sor_i2.unsafe(k + 2) 
                              + xb_sor_i2.unsafe(k + 3) * nzbxat_sor_i2.unsafe(k + 3);
                haromszog[8] += xb_sor_i2.unsafe(k + 0) * nzbxat_sor_i3.unsafe(k + 0) 
                              + xb_sor_i2.unsafe(k + 1) * nzbxat_sor_i3.unsafe(k + 1) 
                              + xb_sor_i2.unsafe(k + 2) * nzbxat_sor_i3.unsafe(k + 2) 
                              + xb_sor_i2.unsafe(k + 3) * nzbxat_sor_i3.unsafe(k + 3);

                haromszog[9] += xb_sor_i3.unsafe(k + 0) * nzbxat_sor_i3.unsafe(k + 0) 
                              + xb_sor_i3.unsafe(k + 1) * nzbxat_sor_i3.unsafe(k + 1) 
                              + xb_sor_i3.unsafe(k + 2) * nzbxat_sor_i3.unsafe(k + 2) 
                              + xb_sor_i3.unsafe(k + 3) * nzbxat_sor_i3.unsafe(k + 3);
            }
            for (uns k = hk; k < nk; k++) {
                haromszog[0] += xb_sor_i0.unsafe(k) * nzbxat_sor_i0.unsafe(k);
                haromszog[1] += xb_sor_i0.unsafe(k) * nzbxat_sor_i1.unsafe(k);
                haromszog[2] += xb_sor_i0.unsafe(k) * nzbxat_sor_i2.unsafe(k);
                haromszog[3] += xb_sor_i0.unsafe(k) * nzbxat_sor_i3.unsafe(k);

                haromszog[4] += xb_sor_i1.unsafe(k) * nzbxat_sor_i1.unsafe(k);
                haromszog[5] += xb_sor_i1.unsafe(k) * nzbxat_sor_i2.unsafe(k);
                haromszog[6] += xb_sor_i1.unsafe(k) * nzbxat_sor_i3.unsafe(k);

                haromszog[7] += xb_sor_i2.unsafe(k) * nzbxat_sor_i2.unsafe(k);
                haromszog[8] += xb_sor_i2.unsafe(k) * nzbxat_sor_i3.unsafe(k);

                haromszog[9] += xb_sor_i3.unsafe(k) * nzbxat_sor_i3.unsafe(k);

            }

            vektor<adattipus> & sor_i0 = sorok.unsafe(i + 0);
            vektor<adattipus> & sor_i1 = sorok.unsafe(i + 1);
            vektor<adattipus> & sor_i2 = sorok.unsafe(i + 2);
            vektor<adattipus> & sor_i3 = sorok.unsafe(i + 3);
            const vektor<adattipus> & ya_sor_i0 = c.sorok.unsafe(i + 0);
            const vektor<adattipus> & ya_sor_i1 = c.sorok.unsafe(i + 1);
            const vektor<adattipus> & ya_sor_i2 = c.sorok.unsafe(i + 2);
            const vektor<adattipus> & ya_sor_i3 = c.sorok.unsafe(i + 3);

            sor_i0.unsafe(i + 0) = -haromszog[0] + ya_sor_i0.unsafe(i + 0);
            sor_i0.unsafe(i + 1) = -haromszog[1] + ya_sor_i0.unsafe(i + 1);
            sor_i0.unsafe(i + 2) = -haromszog[2] + ya_sor_i0.unsafe(i + 2);
            sor_i0.unsafe(i + 3) = -haromszog[3] + ya_sor_i0.unsafe(i + 3);
            
            sor_i1.unsafe(i + 1) = -haromszog[4] + ya_sor_i1.unsafe(i + 1);
            sor_i1.unsafe(i + 2) = -haromszog[5] + ya_sor_i1.unsafe(i + 2);
            sor_i1.unsafe(i + 3) = -haromszog[6] + ya_sor_i1.unsafe(i + 3);

            sor_i2.unsafe(i + 2) = -haromszog[7] + ya_sor_i2.unsafe(i + 2);
            sor_i2.unsafe(i + 3) = -haromszog[8] + ya_sor_i2.unsafe(i + 3);

            sor_i3.unsafe(i + 3) = -haromszog[9] + ya_sor_i3.unsafe(i + 3);

            for (uns j = i + 4; j < hj; j += 4) {
                const vektor<adattipus> & nzbxat_sor_j0 = b.sorok.unsafe(j + 0);
                const vektor<adattipus> & nzbxat_sor_j1 = b.sorok.unsafe(j + 1);
                const vektor<adattipus> & nzbxat_sor_j2 = b.sorok.unsafe(j + 2);
                const vektor<adattipus> & nzbxat_sor_j3 = b.sorok.unsafe(j + 3);
                adattipus negyzet[16] = { adattipus() };
                for (uns k = 0; k < hk; k += 4) {
                    negyzet[0] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_j0.unsafe(k + 0) 
                                + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_j0.unsafe(k + 1) 
                                + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_j0.unsafe(k + 2) 
                                + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_j0.unsafe(k + 3);
                    negyzet[1] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_j1.unsafe(k + 0) 
                                + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_j1.unsafe(k + 1) 
                                + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_j1.unsafe(k + 2) 
                                + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_j1.unsafe(k + 3);
                    negyzet[2] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_j2.unsafe(k + 0) 
                                + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_j2.unsafe(k + 1) 
                                + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_j2.unsafe(k + 2) 
                                + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_j2.unsafe(k + 3);
                    negyzet[3] += xb_sor_i0.unsafe(k + 0) * nzbxat_sor_j3.unsafe(k + 0) 
                                + xb_sor_i0.unsafe(k + 1) * nzbxat_sor_j3.unsafe(k + 1) 
                                + xb_sor_i0.unsafe(k + 2) * nzbxat_sor_j3.unsafe(k + 2) 
                                + xb_sor_i0.unsafe(k + 3) * nzbxat_sor_j3.unsafe(k + 3);

                    negyzet[4] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_j0.unsafe(k + 0) 
                                + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_j0.unsafe(k + 1)
                                + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_j0.unsafe(k + 2)
                                + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_j0.unsafe(k + 3);
                    negyzet[5] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_j1.unsafe(k + 0)
                                + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_j1.unsafe(k + 1)
                                + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_j1.unsafe(k + 2)
                                + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_j1.unsafe(k + 3);
                    negyzet[6] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_j2.unsafe(k + 0)
                                + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_j2.unsafe(k + 1)
                                + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_j2.unsafe(k + 2)
                                + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_j2.unsafe(k + 3);
                    negyzet[7] += xb_sor_i1.unsafe(k + 0) * nzbxat_sor_j3.unsafe(k + 0)
                                + xb_sor_i1.unsafe(k + 1) * nzbxat_sor_j3.unsafe(k + 1)
                                + xb_sor_i1.unsafe(k + 2) * nzbxat_sor_j3.unsafe(k + 2)
                                + xb_sor_i1.unsafe(k + 3) * nzbxat_sor_j3.unsafe(k + 3);

                    negyzet[8] += xb_sor_i2.unsafe(k + 0) * nzbxat_sor_j0.unsafe(k + 0) 
                                + xb_sor_i2.unsafe(k + 1) * nzbxat_sor_j0.unsafe(k + 1)
                                + xb_sor_i2.unsafe(k + 2) * nzbxat_sor_j0.unsafe(k + 2)
                                + xb_sor_i2.unsafe(k + 3) * nzbxat_sor_j0.unsafe(k + 3);
                    negyzet[9] += xb_sor_i2.unsafe(k + 0) * nzbxat_sor_j1.unsafe(k + 0)
                                + xb_sor_i2.unsafe(k + 1) * nzbxat_sor_j1.unsafe(k + 1)
                                + xb_sor_i2.unsafe(k + 2) * nzbxat_sor_j1.unsafe(k + 2)
                                + xb_sor_i2.unsafe(k + 3) * nzbxat_sor_j1.unsafe(k + 3);
                    negyzet[10]+= xb_sor_i2.unsafe(k + 0) * nzbxat_sor_j2.unsafe(k + 0)
                                + xb_sor_i2.unsafe(k + 1) * nzbxat_sor_j2.unsafe(k + 1)
                                + xb_sor_i2.unsafe(k + 2) * nzbxat_sor_j2.unsafe(k + 2)
                                + xb_sor_i2.unsafe(k + 3) * nzbxat_sor_j2.unsafe(k + 3);
                    negyzet[11]+= xb_sor_i2.unsafe(k + 0) * nzbxat_sor_j3.unsafe(k + 0)
                                + xb_sor_i2.unsafe(k + 1) * nzbxat_sor_j3.unsafe(k + 1)
                                + xb_sor_i2.unsafe(k + 2) * nzbxat_sor_j3.unsafe(k + 2)
                                + xb_sor_i2.unsafe(k + 3) * nzbxat_sor_j3.unsafe(k + 3);

                    negyzet[12]+= xb_sor_i3.unsafe(k + 0) * nzbxat_sor_j0.unsafe(k + 0) 
                                + xb_sor_i3.unsafe(k + 1) * nzbxat_sor_j0.unsafe(k + 1)
                                + xb_sor_i3.unsafe(k + 2) * nzbxat_sor_j0.unsafe(k + 2)
                                + xb_sor_i3.unsafe(k + 3) * nzbxat_sor_j0.unsafe(k + 3);
                    negyzet[13]+= xb_sor_i3.unsafe(k + 0) * nzbxat_sor_j1.unsafe(k + 0)
                                + xb_sor_i3.unsafe(k + 1) * nzbxat_sor_j1.unsafe(k + 1)
                                + xb_sor_i3.unsafe(k + 2) * nzbxat_sor_j1.unsafe(k + 2)
                                + xb_sor_i3.unsafe(k + 3) * nzbxat_sor_j1.unsafe(k + 3);
                    negyzet[14]+= xb_sor_i3.unsafe(k + 0) * nzbxat_sor_j2.unsafe(k + 0)
                                + xb_sor_i3.unsafe(k + 1) * nzbxat_sor_j2.unsafe(k + 1)
                                + xb_sor_i3.unsafe(k + 2) * nzbxat_sor_j2.unsafe(k + 2)
                                + xb_sor_i3.unsafe(k + 3) * nzbxat_sor_j2.unsafe(k + 3);
                    negyzet[15]+= xb_sor_i3.unsafe(k + 0) * nzbxat_sor_j3.unsafe(k + 0)
                                + xb_sor_i3.unsafe(k + 1) * nzbxat_sor_j3.unsafe(k + 1)
                                + xb_sor_i3.unsafe(k + 2) * nzbxat_sor_j3.unsafe(k + 2)
                                + xb_sor_i3.unsafe(k + 3) * nzbxat_sor_j3.unsafe(k + 3);
                }
                for (uns k = hk; k < nk; k++) {
                    negyzet[0]  += xb_sor_i0.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyzet[1]  += xb_sor_i0.unsafe(k) * nzbxat_sor_j1.unsafe(k);
                    negyzet[2]  += xb_sor_i0.unsafe(k) * nzbxat_sor_j2.unsafe(k);
                    negyzet[3]  += xb_sor_i0.unsafe(k) * nzbxat_sor_j3.unsafe(k);

                    negyzet[4]  += xb_sor_i1.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyzet[5]  += xb_sor_i1.unsafe(k) * nzbxat_sor_j1.unsafe(k);
                    negyzet[6]  += xb_sor_i1.unsafe(k) * nzbxat_sor_j2.unsafe(k);
                    negyzet[7]  += xb_sor_i1.unsafe(k) * nzbxat_sor_j3.unsafe(k);

                    negyzet[8]  += xb_sor_i2.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyzet[9]  += xb_sor_i2.unsafe(k) * nzbxat_sor_j1.unsafe(k);
                    negyzet[10] += xb_sor_i2.unsafe(k) * nzbxat_sor_j2.unsafe(k);
                    negyzet[11] += xb_sor_i2.unsafe(k) * nzbxat_sor_j3.unsafe(k);

                    negyzet[12] += xb_sor_i3.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyzet[13] += xb_sor_i3.unsafe(k) * nzbxat_sor_j1.unsafe(k);
                    negyzet[14] += xb_sor_i3.unsafe(k) * nzbxat_sor_j2.unsafe(k);
                    negyzet[15] += xb_sor_i3.unsafe(k) * nzbxat_sor_j3.unsafe(k);
                }

                sor_i0.unsafe(j + 0) = -negyzet[0]  + ya_sor_i0.unsafe(j + 0);
                sor_i0.unsafe(j + 1) = -negyzet[1]  + ya_sor_i0.unsafe(j + 1);
                sor_i0.unsafe(j + 2) = -negyzet[2]  + ya_sor_i0.unsafe(j + 2);
                sor_i0.unsafe(j + 3) = -negyzet[3]  + ya_sor_i0.unsafe(j + 3);

                sor_i1.unsafe(j + 0) = -negyzet[4]  + ya_sor_i1.unsafe(j + 0);
                sor_i1.unsafe(j + 1) = -negyzet[5]  + ya_sor_i1.unsafe(j + 1);
                sor_i1.unsafe(j + 2) = -negyzet[6]  + ya_sor_i1.unsafe(j + 2);
                sor_i1.unsafe(j + 3) = -negyzet[7]  + ya_sor_i1.unsafe(j + 3);

                sor_i2.unsafe(j + 0) = -negyzet[8]  + ya_sor_i2.unsafe(j + 0);
                sor_i2.unsafe(j + 1) = -negyzet[9]  + ya_sor_i2.unsafe(j + 1);
                sor_i2.unsafe(j + 2) = -negyzet[10] + ya_sor_i2.unsafe(j + 2);
                sor_i2.unsafe(j + 3) = -negyzet[11] + ya_sor_i2.unsafe(j + 3);

                sor_i3.unsafe(j + 0) = -negyzet[12] + ya_sor_i3.unsafe(j + 0);
                sor_i3.unsafe(j + 1) = -negyzet[13] + ya_sor_i3.unsafe(j + 1);
                sor_i3.unsafe(j + 2) = -negyzet[14] + ya_sor_i3.unsafe(j + 2);
                sor_i3.unsafe(j + 3) = -negyzet[15] + ya_sor_i3.unsafe(j + 3);
            }
            for (uns j = hj; j < nj; j++) {
                const vektor<adattipus> & nzbxat_sor_j0 = b.sorok.unsafe(j);
                adattipus negyes[4] = { adattipus() };
                for (uns k = 0; k < nk; k++) {
                    negyes[0] += xb_sor_i0.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyes[1] += xb_sor_i1.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyes[2] += xb_sor_i2.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                    negyes[3] += xb_sor_i3.unsafe(k) * nzbxat_sor_j0.unsafe(k);
                }
                sor_i0.unsafe(j) = -negyes[0] + ya_sor_i0.unsafe(j);
                sor_i1.unsafe(j) = -negyes[1] + ya_sor_i1.unsafe(j);
                sor_i2.unsafe(j) = -negyes[2] + ya_sor_i2.unsafe(j);
                sor_i3.unsafe(j) = -negyes[3] + ya_sor_i3.unsafe(j);
            }
	    }
        for (uns i = hi; i < ni; i++) {
            const vektor<adattipus> & xb_sor_i = a.sorok.unsafe(i);
            for (uns j = i; j < nj; j++) {
                const vektor<adattipus> & nzbxat_sor_j = b.sorok.unsafe(j);
                adattipus sum = adattipus();
                for (uns k = 0; k < nk; k++) {
                    sum += xb_sor_i.unsafe(k) * nzbxat_sor_j.unsafe(k);
                }
                sorok.unsafe(i).unsafe(j) = -sum + c.sorok.unsafe(i).unsafe(j);
            }
        }
        if(is_symmetrize_needed)
            symmetrize_from_upper();
    }

    //***********************************************************************
    void math_bigblock_ninv_np(uns osztasszam, bool is_priority);
    void math_bigblock_symm_in_nonsymm_ninv_np(uns osztasszam, bool is_priority);
    void math_bigblock_nonsymm_mul_t(uns osztasszam, const matrix & src1, const matrix & src2, bool is_priority);
    void math_bigblock_nonsymm_nmul_t(uns osztasszam, const matrix & src1, const matrix & src2, bool is_priority);
    void math_bigblock_nonsymm_add_mul_t(uns osztasszam, const matrix & src1, const matrix & src2, const matrix & src3, bool is_priority);
    void math_bigblock_nonsymm_sub_mul_t(uns osztasszam, const matrix & src1, const matrix & src2, const matrix & src3, bool is_priority);
    void math_bigblock_symm_add_mul_t(uns osztasszam, const matrix & src1, const matrix & src2, const matrix & src3, bool is_priority);
    void math_bigblock_symm_in_nosymm_sub_mul_t(uns osztasszam, const matrix & src1, const matrix & src2, const matrix & src3, bool is_symmetrize_needed, bool is_priority);
    //***********************************************************************
};

using namespace std::chrono_literals;

//***********************************************************************
template<typename adattipus> class parallel_seged {
//***********************************************************************
    enum para_feladat_tipus{ pft_inv, pft_mul_t, pft_addsub_mul_t };
    //***********************************************************************
    matrix<adattipus> A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P;
    matrix<adattipus> Q1, Q2, Q3, Q4, R1, R2, R3, S1, S2, S3, S4, T1, T2, T3;
    matrix<adattipus> *p_dest;
    const matrix<adattipus> *p_src1, *p_src2, *p_src3;
    uns szalszam;
    para_feladat_tipus tipus;
    bool is_symm, is_vizszintes;
    ::std::mutex lezaro_mutex;
    ::std::condition_variable kesz_egy_feladat;
    std::atomic<bool> is_kesz_1{ false }, is_kesz_2[6]{ false }, is_kesz_3[3]{ false },
        is_kesz_4{ false }, is_kesz_5[3]{ false }, is_kesz_6[3]{ false };
    //***********************************************************************

public:
    //***********************************************************************
    void var_egy_esemenyre(const std::atomic<bool> & is_kesz) {
    //***********************************************************************
        if (is_kesz)
            return;
        ::std::unique_lock<::std::mutex> lock(lezaro_mutex);
        kesz_egy_feladat.wait(lock, [] { return true;});
//        if(kesz_egy_feladat.wait_for(lock, 100ms) == std::cv_status::timeout)
//            printf("\nTime out.\n");
    }

    //***********************************************************************
    void run_job(hivando_matrixfuggveny<adattipus> & fv);
    // A parhuzamos_feldolgozo-ban van definiálva
    //***********************************************************************

    //***********************************************************************
    parallel_seged(uns szalszam, matrix<adattipus> & NZB, bool is_nonsymm_inv) : szalszam{ szalszam }, tipus{ pft_inv }, 
        p_src1{ nullptr }, p_src2{ nullptr }, p_src3{ nullptr }, is_vizszintes{ false } {
    // invertáláshoz állítja be
    //***********************************************************************
        egyforma_e_hiba(NZB.get_row(), NZB.get_col(), "parallel_seged::parallel_seged (inv)");
        igaz_e_hiba(NZB.get_is_szimm(), "parallel_seged::parallel_seged (inv)", "symmetrical matrix not allowed");
        uns COL = NZB.get_col();
        uns n = (COL + 2) / 4;
        uns m = COL - 3 * n;
        A.rafektet(NZB, 0 * n, 0 * n, n, n, false);
        B.rafektet(NZB, 0 * n, 1 * n, n, n, false);
        C.rafektet(NZB, 0 * n, 2 * n, n, n, false);
        D.rafektet(NZB, 0 * n, 3 * n, n, m, false);
        E.rafektet(NZB, 1 * n, 0 * n, n, n, false);
        F.rafektet(NZB, 1 * n, 1 * n, n, n, false);
        G.rafektet(NZB, 1 * n, 2 * n, n, n, false);
        H.rafektet(NZB, 1 * n, 3 * n, n, m, false);
        I.rafektet(NZB, 2 * n, 0 * n, n, n, false);
        J.rafektet(NZB, 2 * n, 1 * n, n, n, false);
        K.rafektet(NZB, 2 * n, 2 * n, n, n, false);
        L.rafektet(NZB, 2 * n, 3 * n, n, m, false);
        M.rafektet(NZB, 3 * n, 0 * n, m, n, false);
        N.rafektet(NZB, 3 * n, 1 * n, m, n, false);
        O.rafektet(NZB, 3 * n, 2 * n, m, n, false);
        P.rafektet(NZB, 3 * n, 3 * n, m, m, false);
        Q1.set_size(n, n);
        Q2.set_size(n, n);
        Q3.set_size(m, n);
        S1.set_size(n, m);
        S2.set_size(n, m);
        S3.set_size(n, m);
        if (is_nonsymm_inv) {
            Q4.set_size(n, n);
            S4.set_size(m, m);
            R1.set_size(n, n);
            R2.set_size(n, n);
            R3.set_size(m, n);
            T1.set_size(n, m);
            T2.set_size(n, m);
            T3.set_size(n, m);
        }
        p_dest = &NZB;
        is_symm = false;
    }

    //***********************************************************************
    parallel_seged(uns szalszam, matrix<adattipus> & dest, const matrix<adattipus> & src1,
        const matrix<adattipus> & src2) : szalszam{ szalszam }, tipus{ pft_mul_t }, p_src3{ nullptr }, is_vizszintes{ false } {
    // szorzáshoz állítja be
    //***********************************************************************
        p_dest = &dest;
        p_src1 = &src1;
        p_src2 = &src2;
        is_symm = false;
        if (szalszam < 2)
            return;

        // Ha legalább 2 szál van, akkor itt csak kettéosztjuk, rekurzív lesz a szálakra bontás
        // A cél hosszabbik oldalát osztjuk fel a szálak arányában
        // A szal_1-es blokk megy külön szálon, a szal_2-es itt, ezért a szal_1-esnél eggyel kevesebb szállal számolunk, ha az rekurzív.

        cuns szal_1 = (szalszam < 6) ?       1        : ((szalszam - 1) / 2);
        cuns szal_2 = (szalszam < 6) ? (szalszam - 1) : (szalszam - szal_1 - 1);
        cuns szal_ossz = szal_1 + szal_2;

        if (dest.get_col() > dest.get_row()) {

            // vízszintesen osztjuk ketté

            is_vizszintes = true;
            cuns db_1 = (((dest.get_col() * szal_1 / szal_ossz) + 2) / 4) * 4;
            cuns db_2 = dest.get_col() - db_1;

            // src1-et szorozzuk src2 felsõ ill. alsó felével, tehát src1-re nem kell ráfektetni
            // const_cast, hogy ne kelljen még const ráfektetõket is csinálni

            A.rafektet(dest, 0, 0, dest.get_row(), db_1, false);
            B.rafektet(dest, 0, db_1, dest.get_row(), db_2, false);
            C.rafektet(const_cast<matrix<adattipus> &>(src2), 0, 0, db_1, src2.get_col(), false);
            D.rafektet(const_cast<matrix<adattipus> &>(src2), db_1, 0, db_2, src2.get_col(), false);
        }
        else {

            // függõlegesen osztjuk ketté

            is_vizszintes = false;
            cuns db_1 = (((dest.get_row() * szal_1 / szal_ossz) + 2) / 4) * 4;
            cuns db_2 = dest.get_row() - db_1;

            // src1 felsõ ill. alsó felét szorozzuk src2-vel, tehát src2-re nem kell ráfektetni
            // const_cast, hogy ne kelljen még const ráfektetõket is csinálni

            A.rafektet(dest, 0, 0, db_1, dest.get_col(), false);
            B.rafektet(dest, db_1, 0, db_2, dest.get_col(), false);
            C.rafektet(const_cast<matrix<adattipus> &>(src1), 0, 0, db_1, src1.get_col(), false);
            D.rafektet(const_cast<matrix<adattipus> &>(src1), db_1, 0, db_2, src1.get_col(), false);
        }
    }

    //***********************************************************************
    parallel_seged(uns szalszam, bool is_symm, bool is_symm_in_nosymm, matrix<adattipus> & dest, const matrix<adattipus> & src1,
        const matrix<adattipus> & src2, const matrix<adattipus> & src3) : szalszam{ szalszam }, 
        is_symm{ is_symm }, tipus{ pft_addsub_mul_t }, is_vizszintes{ false } {
    // add_mul-hoz vagy sub_mul-hoz állítja be
    // dest = src1 +/- src2 * src3
    //***********************************************************************
        p_dest = &dest;
        p_src1 = &src1;
        p_src2 = &src2;
        p_src3 = &src3;
        if (szalszam < 2)
            return;
        if (is_symm) {
            // Ha 2 szál van, akkor 3 részre bontjuk: 1. szálhoz a felsõ és az alsó háromszög, 2. szálhoz a négyzet tartozik
            // Ha >2 szál van, akkor 4 részre osztjuk: a négyzetet is félbe vágjuk. >4 szálnál rekurzívan kell majd hívni
            // 2 szálnál nincs D, H, mert C, G az egész blokkot fedi, de L, M és O mindig van
            // --------   --------     -----   -----
            // \A |C|D|   \E |G|H|     |   |   |   |
            //  \------ =  \------  +  | I |   | K |
            //   \    |     \    |  -  |___| * |___|, O = az L és az M egyben (B kiszámításához kell, ill. 2 szálnál C-hez)
            //    \ B |      \ F |     |   |   | L |
            //     \  |       \  |     | J |   |---|
            //      \ |        \ |     |   |   | M |
            //       \|         \|     -----   -----

            cuns db_1 = ((dest.get_col() + 2) / 4) * 2;
            cuns db_2 = dest.get_col() - db_1;
            cuns db_3 = ((db_2 + 2) / 4) * 2;
            cuns db_4 = db_2 - db_3;

            A.rafektet(dest, 0, 0, db_1, db_1, true);
            B.rafektet(dest, db_1, db_1, db_2, db_2, true);
            E.rafektet(const_cast<matrix<adattipus> &>(src1), 0, 0, db_1, db_1, true);
            F.rafektet(const_cast<matrix<adattipus> &>(src1), db_1, db_1, db_2, db_2, true);

            I.rafektet(const_cast<matrix<adattipus> &>(src2),           0, 0, db_1, src2.get_col(), false);
            J.rafektet(const_cast<matrix<adattipus> &>(src2),        db_1, 0, db_2, src2.get_col(), false);
            K.rafektet(const_cast<matrix<adattipus> &>(src3),           0, 0, db_1, src3.get_col(), false);
            L.rafektet(const_cast<matrix<adattipus> &>(src3),        db_1, 0, db_3, src3.get_col(), false);
            M.rafektet(const_cast<matrix<adattipus> &>(src3), db_1 + db_3, 0, db_4, src3.get_col(), false);
            O.rafektet(const_cast<matrix<adattipus> &>(src3),        db_1, 0, db_2, src3.get_col(), false);

            if (szalszam == 2) {
                C.rafektet(dest, 0, db_1, db_1, db_2, false);
                G.rafektet(const_cast<matrix<adattipus> &>(src1), 0, db_1, db_1, db_2, false);
            }
            else {
                C.rafektet(dest, 0,        db_1, db_1, db_3, false);
                D.rafektet(dest, 0, db_1 + db_3, db_1, db_4, false);
                G.rafektet(const_cast<matrix<adattipus> &>(src1), 0,        db_1, db_1, db_3, false);
                H.rafektet(const_cast<matrix<adattipus> &>(src1), 0, db_1 + db_3, db_1, db_4, false);
            }
        }
        else  if (is_symm_in_nosymm) {
            // Ha 2 szál van, akkor 3 részre bontjuk: 1. szálhoz a felsõ és az alsó háromszög, 2. szálhoz a négyzet tartozik
            // Ha >2 szál van, akkor 4 részre osztjuk: a négyzetet is félbe vágjuk. >4 szálnál rekurzívan kell majd hívni
            // 2 szálnál nincs D, H, mert C, G az egész blokkot fedi, de L, M és O mindig van
            // --------   --------     -----   -----
            // \A |C|D|   \E |G|H|     |   |   |   |
            //  \------ =  \------  +  | I |   | K |
            //   \    |     \    |  -  |___| * |___|, O = az L és az M egyben (B kiszámításához kell, ill. 2 szálnál C-hez)
            //    \ B |      \ F |     |   |   | L |
            //     \  |       \  |     | J |   |---|
            //      \ |        \ |     |   |   | M |
            //       \|         \|     -----   -----

            cuns db_1 = ((dest.get_col() + 2) / 4) * 2;
            cuns db_2 = dest.get_col() - db_1;
            cuns db_3 = ((db_2 + 2) / 4) * 2;
            cuns db_4 = db_2 - db_3;

            A.rafektet(dest, 0, 0, db_1, db_1, false);
            B.rafektet(dest, db_1, db_1, db_2, db_2, false);
            E.rafektet(const_cast<matrix<adattipus> &>(src1), 0, 0, db_1, db_1, false);
            F.rafektet(const_cast<matrix<adattipus> &>(src1), db_1, db_1, db_2, db_2, false);

            I.rafektet(const_cast<matrix<adattipus> &>(src2), 0, 0, db_1, src2.get_col(), false);
            J.rafektet(const_cast<matrix<adattipus> &>(src2), db_1, 0, db_2, src2.get_col(), false);
            K.rafektet(const_cast<matrix<adattipus> &>(src3), 0, 0, db_1, src3.get_col(), false);
            L.rafektet(const_cast<matrix<adattipus> &>(src3), db_1, 0, db_3, src3.get_col(), false);
            M.rafektet(const_cast<matrix<adattipus> &>(src3), db_1 + db_3, 0, db_4, src3.get_col(), false);
            O.rafektet(const_cast<matrix<adattipus> &>(src3), db_1, 0, db_2, src3.get_col(), false);

            if (szalszam == 2) {
                C.rafektet(dest, 0, db_1, db_1, db_2, false);
                G.rafektet(const_cast<matrix<adattipus> &>(src1), 0, db_1, db_1, db_2, false);
            }
            else {
                C.rafektet(dest, 0, db_1, db_1, db_3, false);
                D.rafektet(dest, 0, db_1 + db_3, db_1, db_4, false);
                G.rafektet(const_cast<matrix<adattipus> &>(src1), 0, db_1, db_1, db_3, false);
                H.rafektet(const_cast<matrix<adattipus> &>(src1), 0, db_1 + db_3, db_1, db_4, false);
            }
        }
        else {

            // Ha legalább 2 szál van, akkor itt csak kettéosztjuk, rekurzív lesz a szálakra bontás
            // A cél hosszabbik oldalát osztjuk fel a szálak arányában

            cuns szal_1 = (szalszam < 6) ?       1        : ((szalszam - 1) / 2);
            cuns szal_2 = (szalszam < 6) ? (szalszam - 1) : (szalszam - szal_1 - 1);
            cuns szal_ossz = szal_1 + szal_2;

            if (dest.get_col() > dest.get_row()) { // Ha felosztott mátrixra hívjuk, akkor lehet ilyen is

                // vízszintesen osztjuk ketté

                is_vizszintes = true;
                cuns db_1 = (((dest.get_col() * szal_1 / szal_ossz) + 2) / 4) * 4;
                cuns db_2 = dest.get_col() - db_1;

                // src2-t szorozzuk src3 felsõ ill. alsó felével, tehát src2-re nem kell ráfektetni
                // const_cast, hogy ne kelljen még const ráfektetõket is csinálni

                A.rafektet(dest, 0, 0, dest.get_row(), db_1, false);
                B.rafektet(dest, 0, db_1, dest.get_row(), db_2, false);
                C.rafektet(const_cast<matrix<adattipus> &>(src1), 0, 0, src1.get_row(), db_1, false);
                D.rafektet(const_cast<matrix<adattipus> &>(src1), 0, db_1, src1.get_row(), db_2, false);
                E.rafektet(const_cast<matrix<adattipus> &>(src3), 0, 0, db_1, src3.get_col(), false);
                F.rafektet(const_cast<matrix<adattipus> &>(src3), db_1, 0, db_2, src3.get_col(), false);
            }
            else {

                // függõlegesen osztjuk ketté

                is_vizszintes = false;
                cuns db_1 = (((dest.get_row() * szal_1 / szal_ossz) + 2) / 4) * 4;
                cuns db_2 = dest.get_row() - db_1;

                // src2 felsõ ill. alsó felét szorozzuk src3-mal, tehát src3-ra nem kell ráfektetni
                // const_cast, hogy ne kelljen még const ráfektetõket is csinálni

                A.rafektet(dest, 0, 0, db_1, dest.get_col(), false);
                B.rafektet(dest, db_1, 0, db_2, dest.get_col(), false);
                C.rafektet(const_cast<matrix<adattipus> &>(src1), 0, 0, db_1, src1.get_col(), false);
                D.rafektet(const_cast<matrix<adattipus> &>(src1), db_1, 0, db_2, src1.get_col(), false);
                E.rafektet(const_cast<matrix<adattipus> &>(src2), 0, 0, db_1, src2.get_col(), false);
                F.rafektet(const_cast<matrix<adattipus> &>(src2), db_1, 0, db_2, src2.get_col(), false);
            }
        }
    }

    //***********************************************************************
    //******** Ha külön szálként futtatott mûvelet "rekurzív", azaz  ********
    //******** parallel_seged-et tartalmaz (matrix::...bigblock...), ********
    //******** akkor számára 1-gyel kevesebb szálszámot kell megadni ********
    //******** Ha a bigblock hívás nem külön szálon töténik, akkor   ********
    //******** nem kell levonni a rendelkezésre álló szálak számából ********
    //***********************************************************************

    //***********************************************************************
    void nonsymm_in_nonsymm_ninv(bool is_priority) {
    // az invertálandó mátrix nemszimmetrikus
    //***********************************************************************
        if (tipus != pft_inv)
            throw hiba("parallel_seged::nonsymm_in_nonsymm_ninv", "tipus != pft_inv");

        if (szalszam < 2) {
            p_dest->math_ninv_np();
            return;
        }

        cuns min_szalszam = 4;
        cuns min_meret = 32;
        const bool felbontando = (szalszam >= min_szalszam) && (A.get_col() / 4 >= min_meret);

        if (felbontando) {
//printf("nonsymm_in_nonsymm_ninv %u %u\n", p_dest->get_row(), p_dest->get_col());

            hivando_matrixfuggveny<adattipus> fv;

            cuns al_szalszam_2 = szalszam / 2 > 0 ? szalszam / 2 : 1;
            cuns al_szalszam_3 = szalszam / 3 > 0 ? szalszam / 3 : 1;
            cuns al_szalszam_4 = szalszam / 4 > 0 ? szalszam / 4 : 1;
            cuns al_szalszam_6 = szalszam / 6 > 0 ? szalszam / 6 : 1;
            cuns al_szalszam_8 = szalszam / 8 > 0 ? szalszam / 8 : 1;

            // 1

            is_kesz_1 = true;  fv.math_bigblock_ninv_np(is_kesz_1, kesz_egy_feladat, A, szalszam - 1, true);  run_job(fv);

            Q1.transp(B);
            Q2.transp(C);
            Q3.transp(D);
            R1.copy(E);
            R2.copy(I);
            R3.copy(M);

            while (!is_kesz_1) {
                var_egy_esemenyre(is_kesz_1);
            }

            // 2

            Q4.transp(A);

            is_kesz_2[0] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[0], kesz_egy_feladat, E, al_szalszam_4, R1, Q4, true);         run_job(fv); // ! 2[0]
            is_kesz_2[1] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[1], kesz_egy_feladat, I, al_szalszam_6, R2, Q4, is_priority);  run_job(fv);
            is_kesz_2[2] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[2], kesz_egy_feladat, M, al_szalszam_6, R3, Q4, is_priority);  run_job(fv);
            is_kesz_2[3] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[3], kesz_egy_feladat, B, al_szalszam_8,  A, Q1, is_priority);  run_job(fv);
            is_kesz_2[4] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[4], kesz_egy_feladat, C, al_szalszam_8,  A, Q2, is_priority);  run_job(fv);
            is_kesz_2[5] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[5], kesz_egy_feladat, D, al_szalszam_8,  A, Q3, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_3_elindult = false;
            bool is_4_elindult = false;
            bool is_5_elindult = false;
            bool is_6_elindult = false;
            bool is_2_kesz = false;
            bool is_3_kesz = false;
            bool is_5_kesz = false;
            bool is_6_kesz = false;


            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 3

                if (!is_3_elindult && is_kesz_2[0]) {
                    is_kesz_3[0] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[0], kesz_egy_feladat, F, al_szalszam_4, F, E, Q1, true);         run_job(fv);
                    is_kesz_3[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[1], kesz_egy_feladat, G, al_szalszam_8, G, E, Q2, is_priority);  run_job(fv);
                    is_kesz_3[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[2], kesz_egy_feladat, H, al_szalszam_8, H, E, Q3, is_priority);  run_job(fv);
                    is_3_elindult = true;
                }

                // 4

                if (!is_4_elindult && is_3_elindult && is_kesz_3[0]) {
                    is_kesz_4 = true;  fv.math_bigblock_ninv_np(is_kesz_4, kesz_egy_feladat, F, al_szalszam_4, true);  run_job(fv);
                    is_4_elindult = true;
                }

                // 5

                if (!is_5_elindult && is_kesz_2[1]) {
                    is_kesz_5[0] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[0], kesz_egy_feladat, J, al_szalszam_8, J, I, Q1, is_priority);  run_job(fv);
                    is_kesz_5[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[1], kesz_egy_feladat, K, al_szalszam_8, K, I, Q2, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[2], kesz_egy_feladat, L, al_szalszam_8, L, I, Q3, is_priority);  run_job(fv);
                    is_5_elindult = true;
                }

                // 6

                if (!is_6_elindult && is_kesz_2[2]) {
                    is_kesz_6[0] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[0], kesz_egy_feladat, N, al_szalszam_8, N, M, Q1, is_priority);  run_job(fv);
                    is_kesz_6[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[1], kesz_egy_feladat, O, al_szalszam_8, O, M, Q2, is_priority);  run_job(fv);
                    is_kesz_6[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[2], kesz_egy_feladat, P, al_szalszam_8, P, M, Q3, is_priority);  run_job(fv);
                    is_6_elindult = true;
                }

                // készek összeszeddése

                if (!is_2_kesz) {
                    if (is_3_elindult && is_5_elindult && is_6_elindult) // hogy ne az atomicot kelljen vizsgálni
                        if (is_kesz_2[3] && is_kesz_2[4] && is_kesz_2[5])
                            is_2_kesz = true;
                }

                if (!is_3_kesz && is_3_elindult) {
                    if (is_4_elindult && is_kesz_3[1] && is_kesz_3[2])
                        is_3_kesz = true;
                }

                if (!is_5_kesz && is_5_elindult) {
                    if (is_kesz_5[0] && is_kesz_5[1] && is_kesz_5[2])
                        is_5_kesz = true;
                }

                if (!is_6_kesz && is_6_elindult) {
                    if (is_kesz_6[0] && is_kesz_6[1] && is_kesz_6[2])
                        is_6_kesz = true;
                }

                if (is_2_kesz && is_3_kesz && is_5_kesz && is_6_kesz)
                    is_kesz_1 = true;
            }

            // 7

            Q1.transp(E);
            Q2.transp(G);
            Q3.transp(H);
            R1.copy(B);
            R2.copy(J);
            R3.copy(N);

            // 8

            while (!is_kesz_4) {
                var_egy_esemenyre(is_kesz_4);
            }

            Q4.transp(F);

            is_kesz_2[1] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[1], kesz_egy_feladat, J, al_szalszam_4, R2, Q4, true);         run_job(fv); // ! 2[1]
            is_kesz_2[0] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[0], kesz_egy_feladat, B, al_szalszam_6, R1, Q4, is_priority);  run_job(fv);
            is_kesz_2[2] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[2], kesz_egy_feladat, N, al_szalszam_6, R3, Q4, is_priority);  run_job(fv);
            is_kesz_2[3] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[3], kesz_egy_feladat, E, al_szalszam_8,  F, Q1, is_priority);  run_job(fv);
            is_kesz_2[4] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[4], kesz_egy_feladat, G, al_szalszam_8,  F, Q2, is_priority);  run_job(fv);
            is_kesz_2[5] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[5], kesz_egy_feladat, H, al_szalszam_8,  F, Q3, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_9_elindult = false;
            bool is_10_elindult = false;
            bool is_11_elindult = false;
            bool is_12_elindult = false;
            bool is_8_kesz = false;
            bool is_9_kesz = false;
            bool is_11_kesz = false;
            bool is_12_kesz = false;


            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 9

                if (!is_9_elindult && is_kesz_2[1]) {  // ! 2[1]
                    is_kesz_3[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[1], kesz_egy_feladat, K, al_szalszam_4, K, J, Q2, true);         run_job(fv); // ! 3[1]
                    is_kesz_3[0] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[0], kesz_egy_feladat, I, al_szalszam_8, I, J, Q1, is_priority);  run_job(fv);
                    is_kesz_3[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[2], kesz_egy_feladat, L, al_szalszam_8, L, J, Q3, is_priority);  run_job(fv);
                    is_9_elindult = true;
                }

                // 10

                if (!is_10_elindult && is_9_elindult && is_kesz_3[1]) {  // ! 3[1]
                    is_kesz_4 = true;  fv.math_bigblock_ninv_np(is_kesz_4, kesz_egy_feladat, K, al_szalszam_4, true);  run_job(fv);
                    is_10_elindult = true;
                }

                // 11

                if (!is_11_elindult && is_kesz_2[0]) { // !
                    is_kesz_5[0] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[0], kesz_egy_feladat, A, al_szalszam_8, A, B, Q1, is_priority);  run_job(fv);
                    is_kesz_5[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[1], kesz_egy_feladat, C, al_szalszam_8, C, B, Q2, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[2], kesz_egy_feladat, D, al_szalszam_8, D, B, Q3, is_priority);  run_job(fv);
                    is_11_elindult = true;
                }

                // 12

                if (!is_12_elindult && is_kesz_2[2]) { // !
                    is_kesz_6[0] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[0], kesz_egy_feladat, M, al_szalszam_8, M, N, Q1, is_priority);  run_job(fv);
                    is_kesz_6[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[1], kesz_egy_feladat, O, al_szalszam_8, O, N, Q2, is_priority);  run_job(fv);
                    is_kesz_6[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[2], kesz_egy_feladat, P, al_szalszam_8, P, N, Q3, is_priority);  run_job(fv);
                    is_12_elindult = true;
                }

                // készek összeszeddése

                if (!is_8_kesz) {
                    if (is_9_elindult && is_11_elindult && is_12_elindult) // hogy ne az atomicot kelljen vizsgálni
                        if (is_kesz_2[3] && is_kesz_2[4] && is_kesz_2[5])
                            is_8_kesz = true;
                }

                if (!is_9_kesz && is_9_elindult) {
                    if (is_10_elindult && is_kesz_3[0] && is_kesz_3[2]) // !
                        is_9_kesz = true;
                }

                if (!is_11_kesz && is_11_elindult) {
                    if (is_kesz_5[0] && is_kesz_5[1] && is_kesz_5[2])
                        is_11_kesz = true;
                }

                if (!is_12_kesz && is_12_elindult) {
                    if (is_kesz_6[0] && is_kesz_6[1] && is_kesz_6[2])
                        is_12_kesz = true;
                }

                if (is_8_kesz && is_9_kesz && is_11_kesz && is_12_kesz)
                    is_kesz_1 = true;
            }

            // 13

            Q1.transp(I);
            Q2.transp(J);
            Q3.transp(L);
            R1.copy(C);
            R2.copy(G);
            R3.copy(O);

            while (!is_kesz_4) {
                var_egy_esemenyre(is_kesz_4);
            }

            // 14

            Q4.transp(K);

            is_kesz_2[2] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[2], kesz_egy_feladat, O, al_szalszam_4, R3, Q4, true);         run_job(fv); // ! 2[2]
            is_kesz_2[0] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[0], kesz_egy_feladat, C, al_szalszam_6, R1, Q4, is_priority);  run_job(fv);
            is_kesz_2[1] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[1], kesz_egy_feladat, G, al_szalszam_6, R2, Q4, is_priority);  run_job(fv);
            is_kesz_2[3] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[3], kesz_egy_feladat, I, al_szalszam_8, K,  Q1, is_priority);  run_job(fv);
            is_kesz_2[4] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[4], kesz_egy_feladat, J, al_szalszam_8, K,  Q2, is_priority);  run_job(fv);
            is_kesz_2[5] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[5], kesz_egy_feladat, L, al_szalszam_8, K,  Q3, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_15_elindult = false;
            bool is_16_elindult = false;
            bool is_17_elindult = false;
            bool is_18_elindult = false;
            bool is_14_kesz = false;
            bool is_15_kesz = false;
            bool is_17_kesz = false;
            bool is_18_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 15

                if (!is_15_elindult && is_kesz_2[2]) {  // ! 2[2]
                    is_kesz_3[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[2], kesz_egy_feladat, P, al_szalszam_4, P, O, Q3, true);         run_job(fv); // ! 3[2]
                    is_kesz_3[0] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[0], kesz_egy_feladat, M, al_szalszam_8, M, O, Q1, is_priority);  run_job(fv);
                    is_kesz_3[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[1], kesz_egy_feladat, N, al_szalszam_8, N, O, Q2, is_priority);  run_job(fv);
                    is_15_elindult = true;
                }

                // 16

                if (!is_16_elindult && is_15_elindult && is_kesz_3[2]) {  // ! 3[2]
                    is_kesz_4 = true;  fv.math_bigblock_ninv_np(is_kesz_4, kesz_egy_feladat, P, al_szalszam_4, true);  run_job(fv);
                    is_16_elindult = true;
                }

                // 17

                if (!is_17_elindult && is_kesz_2[0]) { // !
                    is_kesz_5[0] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[0], kesz_egy_feladat, A, al_szalszam_8, A, C, Q1, is_priority);  run_job(fv);
                    is_kesz_5[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[1], kesz_egy_feladat, B, al_szalszam_8, B, C, Q2, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[2], kesz_egy_feladat, D, al_szalszam_8, D, C, Q3, is_priority);  run_job(fv);
                    is_17_elindult = true;
                }

                // 18

                if (!is_18_elindult && is_kesz_2[1]) { // !
                    is_kesz_6[0] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[0], kesz_egy_feladat, E, al_szalszam_8, E, G, Q1, is_priority);  run_job(fv);
                    is_kesz_6[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[1], kesz_egy_feladat, F, al_szalszam_8, F, G, Q2, is_priority);  run_job(fv);
                    is_kesz_6[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[2], kesz_egy_feladat, H, al_szalszam_8, H, G, Q3, is_priority);  run_job(fv);
                    is_18_elindult = true;
                }

                // készek összeszeddése

                if (!is_14_kesz) {
                    if (is_15_elindult && is_17_elindult && is_18_elindult) // hogy ne az atomicot kelljen vizsgálni
                        if (is_kesz_2[3] && is_kesz_2[4] && is_kesz_2[5])
                            is_14_kesz = true;
                }

                if (!is_15_kesz && is_15_elindult) {
                    if (is_16_elindult && is_kesz_3[0] && is_kesz_3[1]) // !
                        is_15_kesz = true;
                }

                if (!is_17_kesz && is_17_elindult) {
                    if (is_kesz_5[0] && is_kesz_5[1] && is_kesz_5[2])
                        is_17_kesz = true;
                }

                if (!is_18_kesz && is_18_elindult) {
                    if (is_kesz_6[0] && is_kesz_6[1] && is_kesz_6[2])
                        is_18_kesz = true;
                }

                if (is_14_kesz && is_15_kesz && is_17_kesz && is_18_kesz)
                    is_kesz_1 = true;
            }

            // 19

            S1.transp(M);
            S2.transp(N);
            S3.transp(O);
            T1.copy(D);
            T2.copy(H);
            T3.copy(L);

            // 20

            while (!is_kesz_4) {
                var_egy_esemenyre(is_kesz_4);
            }

            S4.transp(P);

            is_kesz_2[0] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[0], kesz_egy_feladat, D, al_szalszam_4, T1, S4, is_priority);  run_job(fv);
            is_kesz_2[1] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[1], kesz_egy_feladat, H, al_szalszam_4, T2, S4, is_priority);  run_job(fv);
            is_kesz_2[2] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[2], kesz_egy_feladat, L, al_szalszam_4, T3, S4, is_priority);  run_job(fv);
            is_kesz_2[3] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[3], kesz_egy_feladat, M, al_szalszam_8, P,  S1, is_priority);  run_job(fv);
            is_kesz_2[4] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[4], kesz_egy_feladat, N, al_szalszam_8, P,  S2, is_priority);  run_job(fv);
            is_kesz_2[5] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[5], kesz_egy_feladat, O, al_szalszam_8, P,  S3, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_21_elindult = false;
            bool is_22_elindult = false;
            bool is_23_elindult = false;
            bool is_20_kesz = false;
            bool is_21_kesz = false;
            bool is_22_kesz = false;
            bool is_23_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 21

                if (!is_21_elindult && is_kesz_2[0]) {
                    is_kesz_3[0] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[0], kesz_egy_feladat, A, al_szalszam_8, A, D, S1, is_priority);  run_job(fv);
                    is_kesz_3[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[1], kesz_egy_feladat, B, al_szalszam_8, B, D, S2, is_priority);  run_job(fv);
                    is_kesz_3[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[2], kesz_egy_feladat, C, al_szalszam_8, C, D, S3, is_priority);  run_job(fv);
                    is_21_elindult = true;
                }

                // 22

                if (!is_22_elindult && is_kesz_2[1]) {
                    is_kesz_5[0] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[0], kesz_egy_feladat, E, al_szalszam_8, E, H, S1, is_priority);  run_job(fv);
                    is_kesz_5[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[1], kesz_egy_feladat, F, al_szalszam_8, F, H, S2, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[2], kesz_egy_feladat, G, al_szalszam_8, G, H, S3, is_priority);  run_job(fv);
                    is_22_elindult = true;
                }

                // 23

                if (!is_23_elindult && is_kesz_2[2]) {
                    is_kesz_6[0] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[0], kesz_egy_feladat, I, al_szalszam_8, I, L, S1, is_priority);  run_job(fv);
                    is_kesz_6[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[1], kesz_egy_feladat, J, al_szalszam_8, J, L, S2, is_priority);  run_job(fv);
                    is_kesz_6[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[2], kesz_egy_feladat, K, al_szalszam_8, K, L, S3, is_priority);  run_job(fv);
                    is_23_elindult = true;
                }

                // készek összeszeddése

                if (!is_20_kesz) {
                    if (is_21_elindult && is_22_elindult && is_23_elindult) // hogy ne az atomicot kelljen vizsgálni
                        if (is_kesz_2[3] && is_kesz_2[4] && is_kesz_2[5])
                            is_20_kesz = true;
                }

                if (!is_21_kesz && is_21_elindult) {
                    if (is_kesz_3[0] && is_kesz_3[1] && is_kesz_3[2])
                        is_21_kesz = true;
                }

                if (!is_22_kesz && is_22_elindult) {
                    if (is_kesz_5[0] && is_kesz_5[1] && is_kesz_5[2])
                        is_22_kesz = true;
                }

                if (!is_23_kesz && is_23_elindult) {
                    if (is_kesz_6[0] && is_kesz_6[1] && is_kesz_6[2])
                        is_23_kesz = true;
                }

                if (is_20_kesz && is_21_kesz && is_22_kesz && is_23_kesz)
                    is_kesz_1 = true;
            }
        }
        else {

            hivando_matrixfuggveny<adattipus> fv;

            // 1

            is_kesz_1 = true;  fv.math_ninv_np(is_kesz_1, kesz_egy_feladat, A, true);  run_job(fv);

            Q1.transp(B);
            Q2.transp(C);
            Q3.transp(D);
            R1.copy(E);
            R2.copy(I);
            R3.copy(M);

            while (!is_kesz_1) {
                var_egy_esemenyre(is_kesz_1);
            }

            // 2

            Q4.transp(A);

            is_kesz_2[0] = true;  fv.math_nmul_t_biztos(is_kesz_2[0], kesz_egy_feladat, E, R1, Q4, true);         run_job(fv); // ! 2[0]
            is_kesz_2[1] = true;  fv.math_nmul_t_biztos(is_kesz_2[1], kesz_egy_feladat, I, R2, Q4, is_priority);  run_job(fv);
            is_kesz_2[2] = true;  fv.math_nmul_t_biztos(is_kesz_2[2], kesz_egy_feladat, M, R3, Q4, is_priority);  run_job(fv);
            is_kesz_2[3] = true;  fv.math_nmul_t_biztos(is_kesz_2[3], kesz_egy_feladat, B,  A, Q1, is_priority);  run_job(fv);
            is_kesz_2[4] = true;  fv.math_nmul_t_biztos(is_kesz_2[4], kesz_egy_feladat, C,  A, Q2, is_priority);  run_job(fv);
            is_kesz_2[5] = true;  fv.math_nmul_t_biztos(is_kesz_2[5], kesz_egy_feladat, D,  A, Q3, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_3_elindult = false;
            bool is_4_elindult = false;
            bool is_5_elindult = false;
            bool is_6_elindult = false;
            bool is_2_kesz = false;
            bool is_3_kesz = false;
            bool is_5_kesz = false;
            bool is_6_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 3

                if (!is_3_elindult && is_kesz_2[0]) {
                    is_kesz_3[0] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[0], kesz_egy_feladat, F, F, E, Q1, true);         run_job(fv);
                    is_kesz_3[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[1], kesz_egy_feladat, G, G, E, Q2, is_priority);  run_job(fv);
                    is_kesz_3[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[2], kesz_egy_feladat, H, H, E, Q3, is_priority);  run_job(fv);
                    is_3_elindult = true;
                }

                // 4

                if (!is_4_elindult && is_3_elindult && is_kesz_3[0]) {
                    is_kesz_4 = true;  fv.math_ninv_np(is_kesz_4, kesz_egy_feladat, F, true);  run_job(fv);
                    is_4_elindult = true;
                }

                // 5

                if (!is_5_elindult && is_kesz_2[1]) {
                    is_kesz_5[0] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[0], kesz_egy_feladat, J, J, I, Q1, is_priority);  run_job(fv);
                    is_kesz_5[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[1], kesz_egy_feladat, K, K, I, Q2, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[2], kesz_egy_feladat, L, L, I, Q3, is_priority);  run_job(fv);
                    is_5_elindult = true;
                }

                // 6

                if (!is_6_elindult && is_kesz_2[2]) {
                    is_kesz_6[0] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[0], kesz_egy_feladat, N, N, M, Q1, is_priority);  run_job(fv);
                    is_kesz_6[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[1], kesz_egy_feladat, O, O, M, Q2, is_priority);  run_job(fv);
                    is_kesz_6[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[2], kesz_egy_feladat, P, P, M, Q3, is_priority);  run_job(fv);
                    is_6_elindult = true;
                }

                // készek összeszeddése

                if (!is_2_kesz) {
                    if (is_3_elindult && is_5_elindult && is_6_elindult) // hogy ne az atomicot kelljen vizsgálni
                        if (is_kesz_2[3] && is_kesz_2[4] && is_kesz_2[5])
                            is_2_kesz = true;
                }

                if (!is_3_kesz && is_3_elindult) {
                    if (is_4_elindult && is_kesz_3[1] && is_kesz_3[2])
                        is_3_kesz = true;
                }

                if (!is_5_kesz && is_5_elindult) {
                    if (is_kesz_5[0] && is_kesz_5[1] && is_kesz_5[2])
                        is_5_kesz = true;
                }

                if (!is_6_kesz && is_6_elindult) {
                    if (is_kesz_6[0] && is_kesz_6[1] && is_kesz_6[2])
                        is_6_kesz = true;
                }

                if (is_2_kesz && is_3_kesz && is_5_kesz && is_6_kesz)
                    is_kesz_1 = true;
            }

            // 7

            Q1.transp(E);
            Q2.transp(G);
            Q3.transp(H);
            R1.copy(B);
            R2.copy(J);
            R3.copy(N);

            // 8

            while (!is_kesz_4) {
                var_egy_esemenyre(is_kesz_4);
            }

            Q4.transp(F);

            is_kesz_2[1] = true;  fv.math_nmul_t_biztos(is_kesz_2[1], kesz_egy_feladat, J, R2, Q4, true);         run_job(fv); // ! 2[1]
            is_kesz_2[0] = true;  fv.math_nmul_t_biztos(is_kesz_2[0], kesz_egy_feladat, B, R1, Q4, is_priority);  run_job(fv);
            is_kesz_2[2] = true;  fv.math_nmul_t_biztos(is_kesz_2[2], kesz_egy_feladat, N, R3, Q4, is_priority);  run_job(fv);
            is_kesz_2[3] = true;  fv.math_nmul_t_biztos(is_kesz_2[3], kesz_egy_feladat, E,  F, Q1, is_priority);  run_job(fv);
            is_kesz_2[4] = true;  fv.math_nmul_t_biztos(is_kesz_2[4], kesz_egy_feladat, G,  F, Q2, is_priority);  run_job(fv);
            is_kesz_2[5] = true;  fv.math_nmul_t_biztos(is_kesz_2[5], kesz_egy_feladat, H,  F, Q3, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_9_elindult = false;
            bool is_10_elindult = false;
            bool is_11_elindult = false;
            bool is_12_elindult = false;
            bool is_8_kesz = false;
            bool is_9_kesz = false;
            bool is_11_kesz = false;
            bool is_12_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 9

                if (!is_9_elindult && is_kesz_2[1]) {  // ! 2[1]
                    is_kesz_3[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[1], kesz_egy_feladat, K, K, J, Q2, true);         run_job(fv); // ! 3[1]
                    is_kesz_3[0] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[0], kesz_egy_feladat, I, I, J, Q1, is_priority);  run_job(fv);
                    is_kesz_3[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[2], kesz_egy_feladat, L, L, J, Q3, is_priority);  run_job(fv);
                    is_9_elindult = true;
                }

                // 10

                if (!is_10_elindult && is_9_elindult && is_kesz_3[1]) {  // ! 3[1]
                    is_kesz_4 = true;  fv.math_ninv_np(is_kesz_4, kesz_egy_feladat, K, true);  run_job(fv);
                    is_10_elindult = true;
                }

                // 11

                if (!is_11_elindult && is_kesz_2[0]) { // !
                    is_kesz_5[0] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[0], kesz_egy_feladat, A, A, B, Q1, is_priority);  run_job(fv);
                    is_kesz_5[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[1], kesz_egy_feladat, C, C, B, Q2, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[2], kesz_egy_feladat, D, D, B, Q3, is_priority);  run_job(fv);
                    is_11_elindult = true;
                }

                // 12

                if (!is_12_elindult && is_kesz_2[2]) { // !
                    is_kesz_6[0] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[0], kesz_egy_feladat, M, M, N, Q1, is_priority);  run_job(fv);
                    is_kesz_6[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[1], kesz_egy_feladat, O, O, N, Q2, is_priority);  run_job(fv);
                    is_kesz_6[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[2], kesz_egy_feladat, P, P, N, Q3, is_priority);  run_job(fv);
                    is_12_elindult = true;
                }

                // készek összeszeddése

                if (!is_8_kesz) {
                    if (is_9_elindult && is_11_elindult && is_12_elindult) // hogy ne az atomicot kelljen vizsgálni
                        if (is_kesz_2[3] && is_kesz_2[4] && is_kesz_2[5])
                            is_8_kesz = true;
                }

                if (!is_9_kesz && is_9_elindult) {
                    if (is_10_elindult && is_kesz_3[0] && is_kesz_3[2]) // !
                        is_9_kesz = true;
                }

                if (!is_11_kesz && is_11_elindult) {
                    if (is_kesz_5[0] && is_kesz_5[1] && is_kesz_5[2])
                        is_11_kesz = true;
                }

                if (!is_12_kesz && is_12_elindult) {
                    if (is_kesz_6[0] && is_kesz_6[1] && is_kesz_6[2])
                        is_12_kesz = true;
                }

                if (is_8_kesz && is_9_kesz && is_11_kesz && is_12_kesz)
                    is_kesz_1 = true;
            }

            // 13

            Q1.transp(I);
            Q2.transp(J);
            Q3.transp(L);
            R1.copy(C);
            R2.copy(G);
            R3.copy(O);

            while (!is_kesz_4) {
                var_egy_esemenyre(is_kesz_4);
            }

            // 14

            Q4.transp(K);

            is_kesz_2[2] = true;  fv.math_nmul_t_biztos(is_kesz_2[2], kesz_egy_feladat, O, R3, Q4, true);         run_job(fv); // ! 2[2]
            is_kesz_2[0] = true;  fv.math_nmul_t_biztos(is_kesz_2[0], kesz_egy_feladat, C, R1, Q4, is_priority);  run_job(fv);
            is_kesz_2[1] = true;  fv.math_nmul_t_biztos(is_kesz_2[1], kesz_egy_feladat, G, R2, Q4, is_priority);  run_job(fv);
            is_kesz_2[3] = true;  fv.math_nmul_t_biztos(is_kesz_2[3], kesz_egy_feladat, I,  K, Q1, is_priority);  run_job(fv);
            is_kesz_2[4] = true;  fv.math_nmul_t_biztos(is_kesz_2[4], kesz_egy_feladat, J,  K, Q2, is_priority);  run_job(fv);
            is_kesz_2[5] = true;  fv.math_nmul_t_biztos(is_kesz_2[5], kesz_egy_feladat, L,  K, Q3, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_15_elindult = false;
            bool is_16_elindult = false;
            bool is_17_elindult = false;
            bool is_18_elindult = false;
            bool is_14_kesz = false;
            bool is_15_kesz = false;
            bool is_17_kesz = false;
            bool is_18_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 15

                if (!is_15_elindult && is_kesz_2[2]) {  // ! 2[2]
                    is_kesz_3[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[2], kesz_egy_feladat, P, P, O, Q3, true);         run_job(fv); // ! 3[2]
                    is_kesz_3[0] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[0], kesz_egy_feladat, M, M, O, Q1, is_priority);  run_job(fv);
                    is_kesz_3[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[1], kesz_egy_feladat, N, N, O, Q2, is_priority);  run_job(fv);
                    is_15_elindult = true;
                }

                // 16

                if (!is_16_elindult && is_15_elindult && is_kesz_3[2]) {  // ! 3[2]
                    is_kesz_4 = true;  fv.math_ninv_np(is_kesz_4, kesz_egy_feladat, P, true);  run_job(fv);
                    is_16_elindult = true;
                }

                // 17

                if (!is_17_elindult && is_kesz_2[0]) { // !
                    is_kesz_5[0] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[0], kesz_egy_feladat, A, A, C, Q1, is_priority);  run_job(fv);
                    is_kesz_5[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[1], kesz_egy_feladat, B, B, C, Q2, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[2], kesz_egy_feladat, D, D, C, Q3, is_priority);  run_job(fv);
                    is_17_elindult = true;
                }

                // 18

                if (!is_18_elindult && is_kesz_2[1]) { // !
                    is_kesz_6[0] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[0], kesz_egy_feladat, E, E, G, Q1, is_priority);  run_job(fv);
                    is_kesz_6[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[1], kesz_egy_feladat, F, F, G, Q2, is_priority);  run_job(fv);
                    is_kesz_6[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[2], kesz_egy_feladat, H, H, G, Q3, is_priority);  run_job(fv);
                    is_18_elindult = true;
                }

                // készek összeszeddése

                if (!is_14_kesz) {
                    if (is_15_elindult && is_17_elindult && is_18_elindult) // hogy ne az atomicot kelljen vizsgálni
                        if (is_kesz_2[3] && is_kesz_2[4] && is_kesz_2[5])
                            is_14_kesz = true;
                }

                if (!is_15_kesz && is_15_elindult) {
                    if (is_16_elindult && is_kesz_3[0] && is_kesz_3[1]) // !
                        is_15_kesz = true;
                }

                if (!is_17_kesz && is_17_elindult) {
                    if (is_kesz_5[0] && is_kesz_5[1] && is_kesz_5[2])
                        is_17_kesz = true;
                }

                if (!is_18_kesz && is_18_elindult) {
                    if (is_kesz_6[0] && is_kesz_6[1] && is_kesz_6[2])
                        is_18_kesz = true;
                }

                if (is_14_kesz && is_15_kesz && is_17_kesz && is_18_kesz)
                    is_kesz_1 = true;
            }

            // 19

            S1.transp(M);
            S2.transp(N);
            S3.transp(O);
            T1.copy(D);
            T2.copy(H);
            T3.copy(L);

            // 20

            while (!is_kesz_4) {
                var_egy_esemenyre(is_kesz_4);
            }

            S4.transp(P);

            is_kesz_2[0] = true;  fv.math_nmul_t_biztos(is_kesz_2[0], kesz_egy_feladat, D, T1, S4, is_priority);  run_job(fv);
            is_kesz_2[1] = true;  fv.math_nmul_t_biztos(is_kesz_2[1], kesz_egy_feladat, H, T2, S4, is_priority);  run_job(fv);
            is_kesz_2[2] = true;  fv.math_nmul_t_biztos(is_kesz_2[2], kesz_egy_feladat, L, T3, S4, is_priority);  run_job(fv);
            is_kesz_2[3] = true;  fv.math_nmul_t_biztos(is_kesz_2[3], kesz_egy_feladat, M, P,  S1, is_priority);  run_job(fv);
            is_kesz_2[4] = true;  fv.math_nmul_t_biztos(is_kesz_2[4], kesz_egy_feladat, N, P,  S2, is_priority);  run_job(fv);
            is_kesz_2[5] = true;  fv.math_nmul_t_biztos(is_kesz_2[5], kesz_egy_feladat, O, P,  S3, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_21_elindult = false;
            bool is_22_elindult = false;
            bool is_23_elindult = false;
            bool is_20_kesz = false;
            bool is_21_kesz = false;
            bool is_22_kesz = false;
            bool is_23_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 21

                if (!is_21_elindult && is_kesz_2[0]) {
                    is_kesz_3[0] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[0], kesz_egy_feladat, A, A, D, S1, is_priority);  run_job(fv);
                    is_kesz_3[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[1], kesz_egy_feladat, B, B, D, S2, is_priority);  run_job(fv);
                    is_kesz_3[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[2], kesz_egy_feladat, C, C, D, S3, is_priority);  run_job(fv);
                    is_21_elindult = true;
                }

                // 22

                if (!is_22_elindult && is_kesz_2[1]) {
                    is_kesz_5[0] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[0], kesz_egy_feladat, E, E, H, S1, is_priority);  run_job(fv);
                    is_kesz_5[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[1], kesz_egy_feladat, F, F, H, S2, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[2], kesz_egy_feladat, G, G, H, S3, is_priority);  run_job(fv);
                    is_22_elindult = true;
                }

                // 23

                if (!is_23_elindult && is_kesz_2[2]) {
                    is_kesz_6[0] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[0], kesz_egy_feladat, I, I, L, S1, is_priority);  run_job(fv);
                    is_kesz_6[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[1], kesz_egy_feladat, J, J, L, S2, is_priority);  run_job(fv);
                    is_kesz_6[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[2], kesz_egy_feladat, K, K, L, S3, is_priority);  run_job(fv);
                    is_23_elindult = true;
                }

                // készek összeszeddése

                if (!is_20_kesz) {
                    if (is_21_elindult && is_22_elindult && is_23_elindult) // hogy ne az atomicot kelljen vizsgálni
                        if (is_kesz_2[3] && is_kesz_2[4] && is_kesz_2[5])
                            is_20_kesz = true;
                }

                if (!is_21_kesz && is_21_elindult) {
                    if (is_kesz_3[0] && is_kesz_3[1] && is_kesz_3[2])
                        is_21_kesz = true;
                }

                if (!is_22_kesz && is_22_elindult) {
                    if (is_kesz_5[0] && is_kesz_5[1] && is_kesz_5[2])
                        is_22_kesz = true;
                }

                if (!is_23_kesz && is_23_elindult) {
                    if (is_kesz_6[0] && is_kesz_6[1] && is_kesz_6[2])
                        is_23_kesz = true;
                }

                if (is_20_kesz && is_21_kesz && is_22_kesz && is_23_kesz)
                    is_kesz_1 = true;
            }
        }
    }
    //***********************************************************************
    void symm_in_nonsymm_ninv(bool is_priority) {
    // az invertálandó mátrix nemszimmetrikus mátrixban tárolt szimmetrikus
    //***********************************************************************
        if (tipus != pft_inv)
            throw hiba("parallel_seged::symm_in_nonsymm_ninv", "tipus != pft_inv");

        if (szalszam < 2) {
            p_dest->math_symm_ninv_of_nonsymm();
            return;
        }

        cuns min_szalszam = 4;
        cuns min_meret = 32;
        const bool felbontando = (szalszam >= min_szalszam) && (A.get_col() / 4 >= min_meret);

        if (felbontando) {
//printf("symm_in_nonsymm_ninv %u %u\n", p_dest->get_row(), p_dest->get_col());

            cuns al_inv_szalszam = szalszam / 4;
            cuns al_mul_szalszam = szalszam / 2;

            hivando_matrixfuggveny<adattipus> fv;

            cuns al_szalszam_2 = szalszam / 2 > 0 ? szalszam / 2 : 1;
            cuns al_szalszam_4 = szalszam / 4 > 0 ? szalszam / 4 : 1;

            // 1

            is_kesz_1 = true;  fv.math_bigblock_symm_in_nonsymm_ninv_np(is_kesz_1, kesz_egy_feladat, A, szalszam - 1, true);  run_job(fv);

            Q1.transp(B);
            Q2.transp(C);
            Q3.transp(D);

            while (!is_kesz_1) {
                var_egy_esemenyre(is_kesz_1);
            }

            // 2

            is_kesz_2[0] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[0], kesz_egy_feladat, E, al_szalszam_2, Q1, A, true);         run_job(fv); // ! 2[0]
            is_kesz_2[1] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[1], kesz_egy_feladat, I, al_szalszam_4, Q2, A, is_priority);  run_job(fv);
            is_kesz_2[2] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[2], kesz_egy_feladat, M, al_szalszam_4, Q3, A, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_3_elindult = false;
            bool is_4_elindult = false;
            bool is_5_elindult = false;
            bool is_6_elindult = false;
            bool is_3_kesz = false;
            bool is_5_kesz = false;
            bool is_6_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 3

                if (!is_3_elindult && is_kesz_2[0]) {
                    is_kesz_3[0] = true;  fv.math_bigblock_symm_in_nosymm_sub_mul_t(is_kesz_3[0], kesz_egy_feladat, F, al_szalszam_2, F, E, Q1, false, true);  run_job(fv);
                    is_kesz_3[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[1], kesz_egy_feladat, G, al_szalszam_4, G, E, Q2, is_priority);  run_job(fv);
                    is_kesz_3[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[2], kesz_egy_feladat, H, al_szalszam_4, H, E, Q3, is_priority);  run_job(fv);
                    B.transp(E);
                    is_3_elindult = true;
                }

                // 4

                if (!is_4_elindult && is_3_elindult && is_kesz_3[0]) {
                    is_kesz_4 = true;  fv.math_bigblock_symm_in_nonsymm_ninv_np(is_kesz_4, kesz_egy_feladat, F, al_szalszam_4, true);  run_job(fv);
                    is_4_elindult = true;
                }

                // 5

                if (!is_5_elindult && is_kesz_2[1]) {
                    is_kesz_5[1] = true;  fv.math_bigblock_symm_in_nosymm_sub_mul_t(is_kesz_5[1], kesz_egy_feladat, K, al_szalszam_4, K, I, Q2, true, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[2], kesz_egy_feladat, L, al_szalszam_4, L, I, Q3, is_priority);  run_job(fv);
                    C.transp(I);
                    is_5_elindult = true;
                }

                // 6

                if (!is_6_elindult && is_kesz_2[2]) {
                    is_kesz_6[2] = true;  fv.math_bigblock_symm_in_nosymm_sub_mul_t(is_kesz_6[2], kesz_egy_feladat, P, al_szalszam_4, P, M, Q3, true, is_priority);  run_job(fv);
                    D.transp(M);
                    is_6_elindult = true;
                }

                // készek összeszeddése

                if (!is_3_kesz && is_3_elindult) {
                    if (is_4_elindult && is_kesz_3[1] && is_kesz_3[2])
                        is_3_kesz = true;
                }

                if (!is_5_kesz && is_5_elindult) {
                    if (is_kesz_5[1] && is_kesz_5[2])
                        is_5_kesz = true;
                }

                if (!is_6_kesz && is_6_elindult) {
                    if (is_kesz_6[2])
                        is_6_kesz = true;
                }

                if (is_3_kesz && is_5_kesz && is_6_kesz)
                    is_kesz_1 = true;
            }

            // 7

            Q1.copy(B);
            Q2.transp(G);
            Q3.transp(H);

            while (!is_kesz_4) {
                var_egy_esemenyre(is_kesz_4);
            }

            // 8

            is_kesz_2[1] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[1], kesz_egy_feladat, J, al_szalszam_2, Q2, F, true);         run_job(fv); // ! 2[1]
            is_kesz_2[0] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[0], kesz_egy_feladat, B, al_szalszam_4, Q1, F, is_priority);  run_job(fv);
            is_kesz_2[2] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[2], kesz_egy_feladat, N, al_szalszam_4, Q3, F, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_9_elindult = false;
            bool is_10_elindult = false;
            bool is_11_elindult = false;
            bool is_12_elindult = false;
            bool is_9_kesz = false;
            bool is_11_kesz = false;
            bool is_12_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 9

                if (!is_9_elindult && is_kesz_2[1]) {  // ! 2[1]
                    is_kesz_3[1] = true;  fv.math_bigblock_symm_in_nosymm_sub_mul_t(is_kesz_3[1], kesz_egy_feladat, K, al_szalszam_2, K, J, Q2, false, true);         run_job(fv); // ! 3[1]
                    is_kesz_3[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[2], kesz_egy_feladat, L, al_szalszam_4, L, J, Q3, is_priority);  run_job(fv);
                    G.transp(J);
                    is_9_elindult = true;
                }

                // 10

                if (!is_10_elindult && is_9_elindult && is_kesz_3[1]) {  // ! 3[1]
                    is_kesz_4 = true;  fv.math_bigblock_symm_in_nonsymm_ninv_np(is_kesz_4, kesz_egy_feladat, K, al_szalszam_4, true);  run_job(fv);
                    is_10_elindult = true;
                }

                // 11

                if (!is_11_elindult && is_kesz_2[0]) { // !
                    is_kesz_5[0] = true;  fv.math_bigblock_symm_in_nosymm_sub_mul_t(is_kesz_5[0], kesz_egy_feladat, A, al_szalszam_4, A, B, Q1, true, is_priority);  run_job(fv);
                    is_kesz_5[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[1], kesz_egy_feladat, C, al_szalszam_4, C, B, Q2, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[2], kesz_egy_feladat, D, al_szalszam_4, D, B, Q3, is_priority);  run_job(fv);
                    is_11_elindult = true;
                }

                // 12

                if (!is_12_elindult && is_kesz_2[2]) { // !
                    is_kesz_6[2] = true;  fv.math_bigblock_symm_in_nosymm_sub_mul_t(is_kesz_6[2], kesz_egy_feladat, P, al_szalszam_4, P, N, Q3, true, is_priority);  run_job(fv);
                    H.transp(N);
                    is_12_elindult = true;
                }

                // készek összeszeddése

                if (!is_9_kesz && is_9_elindult) {
                    if (is_10_elindult && is_kesz_3[2]) // !
                        is_9_kesz = true;
                }

                if (!is_11_kesz && is_11_elindult) {
                    if (is_kesz_5[0] && is_kesz_5[1] && is_kesz_5[2])
                        is_11_kesz = true;
                }

                if (!is_12_kesz && is_12_elindult) {
                    if (is_kesz_6[2])
                        is_12_kesz = true;
                }

                if (is_9_kesz && is_11_kesz && is_12_kesz)
                    is_kesz_1 = true;
            }

            // 13

            Q1.copy(C);
            Q2.copy(G);
            Q3.transp(L);

            while (!is_kesz_4) {
                var_egy_esemenyre(is_kesz_4);
            }

            // 14

            is_kesz_2[2] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[2], kesz_egy_feladat, O, al_szalszam_2, Q3, K, true);         run_job(fv); // ! 2[2]
            is_kesz_2[0] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[0], kesz_egy_feladat, C, al_szalszam_4, Q1, K, is_priority);  run_job(fv);
            is_kesz_2[1] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[1], kesz_egy_feladat, G, al_szalszam_4, Q2, K, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_15_elindult = false;
            bool is_16_elindult = false;
            bool is_17_elindult = false;
            bool is_18_elindult = false;
            bool is_15_kesz = false;
            bool is_17_kesz = false;
            bool is_18_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 15

                if (!is_15_elindult && is_kesz_2[2]) {  // ! 2[2]
                    is_kesz_3[2] = true;  fv.math_bigblock_symm_in_nosymm_sub_mul_t(is_kesz_3[2], kesz_egy_feladat, P, al_szalszam_2, P, O, Q3, false, true);  run_job(fv); // ! 3[2]
                    L.transp(O);
                    is_15_elindult = true;
                }

                // 16

                if (!is_16_elindult && is_15_elindult && is_kesz_3[2]) {  // ! 3[2]
                    is_kesz_4 = true;  fv.math_bigblock_symm_in_nonsymm_ninv_np(is_kesz_4, kesz_egy_feladat, P, al_szalszam_4, true);  run_job(fv);
                    is_16_elindult = true;
                }

                // 17

                if (!is_17_elindult && is_kesz_2[0]) { // !
                    is_kesz_5[0] = true;  fv.math_bigblock_symm_in_nosymm_sub_mul_t(is_kesz_5[0], kesz_egy_feladat, A, al_szalszam_4, A, C, Q1, true, is_priority);  run_job(fv);
                    is_kesz_5[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[1], kesz_egy_feladat, B, al_szalszam_4, B, C, Q2, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[2], kesz_egy_feladat, D, al_szalszam_4, D, C, Q3, is_priority);  run_job(fv);
                    is_17_elindult = true;
                }

                // 18

                if (!is_18_elindult && is_kesz_2[1]) { // !
                    is_kesz_6[1] = true;  fv.math_bigblock_symm_in_nosymm_sub_mul_t(is_kesz_6[1], kesz_egy_feladat, F, al_szalszam_4, F, G, Q2, true, is_priority);  run_job(fv);
                    is_kesz_6[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_6[2], kesz_egy_feladat, H, al_szalszam_4, H, G, Q3, is_priority);  run_job(fv);
                    is_18_elindult = true;
                }

                // készek összeszeddése

                if (!is_15_kesz && is_15_elindult) {
                    if (is_16_elindult) // !
                        is_15_kesz = true;
                }

                if (!is_17_kesz && is_17_elindult) {
                    if (is_kesz_5[0] && is_kesz_5[1] && is_kesz_5[2])
                        is_17_kesz = true;
                }

                if (!is_18_kesz && is_18_elindult) {
                    if (is_kesz_6[1] && is_kesz_6[2])
                        is_18_kesz = true;
                }

                if (is_15_kesz && is_17_kesz && is_18_kesz)
                    is_kesz_1 = true;
            }

            // 19

            S1.copy(D);
            S2.copy(H);
            S3.copy(L);

            while (!is_kesz_4) {
                var_egy_esemenyre(is_kesz_4);
            }

            // 20

            is_kesz_2[0] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[0], kesz_egy_feladat, D, al_szalszam_2, S1, P, is_priority);  run_job(fv);
            is_kesz_2[1] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[1], kesz_egy_feladat, H, al_szalszam_4, S2, P, is_priority);  run_job(fv);
            is_kesz_2[2] = true;  fv.math_bigblock_nonsymm_nmul_t(is_kesz_2[2], kesz_egy_feladat, L, al_szalszam_4, S3, P, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_21_elindult = false;
            bool is_22_elindult = false;
            bool is_23_elindult = false;
            bool is_21_kesz = false;
            bool is_22_kesz = false;
            bool is_23_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 21

                if (!is_21_elindult && is_kesz_2[0]) {
                    is_kesz_3[0] = true;  fv.math_bigblock_symm_in_nosymm_sub_mul_t(is_kesz_3[0], kesz_egy_feladat, A, al_szalszam_4, A, D, S1, false, is_priority);  run_job(fv);
                    is_kesz_3[1] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[1], kesz_egy_feladat, B, al_szalszam_4, B, D, S2, is_priority);  run_job(fv);
                    is_kesz_3[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_3[2], kesz_egy_feladat, C, al_szalszam_4, C, D, S3, is_priority);  run_job(fv);
                    is_21_elindult = true;
                }

                // 22

                if (!is_22_elindult && is_kesz_2[1]) {
                    is_kesz_5[1] = true;  fv.math_bigblock_symm_in_nosymm_sub_mul_t(is_kesz_5[1], kesz_egy_feladat, F, al_szalszam_4, F, H, S2, false, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_5[2], kesz_egy_feladat, G, al_szalszam_4, G, H, S3, is_priority);  run_job(fv);
                    is_22_elindult = true;
                }

                // 23

                if (!is_23_elindult && is_kesz_2[2]) {
                    is_kesz_6[2] = true;  fv.math_bigblock_symm_in_nosymm_sub_mul_t(is_kesz_6[2], kesz_egy_feladat, K, al_szalszam_4, K, L, S3, false, is_priority);  run_job(fv);
                    is_23_elindult = true;
                }

                // készek összeszeddése

                if (!is_21_kesz && is_21_elindult) {
                    if (is_kesz_3[0] && is_kesz_3[1] && is_kesz_3[2])
                        is_21_kesz = true;
                }

                if (!is_22_kesz && is_22_elindult) {
                    if (is_kesz_5[1] && is_kesz_5[2])
                        is_22_kesz = true;
                }

                if (!is_23_kesz && is_23_elindult) {
                    if (is_kesz_6[2])
                        is_23_kesz = true;
                }

                if (is_21_kesz && is_22_kesz && is_23_kesz)
                    is_kesz_1 = true;
            }
        }
        else {

            hivando_matrixfuggveny<adattipus> fv;

            // 1

            is_kesz_1 = true;  fv.math_symm_ninv_of_nonsymm(is_kesz_1, kesz_egy_feladat, A, true);  run_job(fv);

            Q1.transp(B);
            Q2.transp(C);
            Q3.transp(D);

            while (!is_kesz_1) {
                var_egy_esemenyre(is_kesz_1);
            }

            // 2

            is_kesz_2[0] = true;  fv.math_nmul_t_biztos(is_kesz_2[0], kesz_egy_feladat, E, Q1, A, true);         run_job(fv); // ! 2[0]
            is_kesz_2[1] = true;  fv.math_nmul_t_biztos(is_kesz_2[1], kesz_egy_feladat, I, Q2, A, is_priority);  run_job(fv);
            is_kesz_2[2] = true;  fv.math_nmul_t_biztos(is_kesz_2[2], kesz_egy_feladat, M, Q3, A, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_3_elindult = false;
            bool is_4_elindult = false;
            bool is_5_elindult = false;
            bool is_6_elindult = false;
            bool is_3_kesz = false;
            bool is_5_kesz = false;
            bool is_6_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 3

                if (!is_3_elindult && is_kesz_2[0]) {
                    is_kesz_3[0] = true;  fv.math_sub_mul_t_symm_in_nonsymm(is_kesz_3[0], kesz_egy_feladat, F, F, E, Q1, false, true);  run_job(fv);
                    is_kesz_3[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[1], kesz_egy_feladat, G, G, E, Q2, is_priority);  run_job(fv);
                    is_kesz_3[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[2], kesz_egy_feladat, H, H, E, Q3, is_priority);  run_job(fv);
                    B.transp(E);
                    is_3_elindult = true;
                }

                // 4

                if (!is_4_elindult && is_3_elindult && is_kesz_3[0]) {
                    is_kesz_4 = true;  fv.math_symm_ninv_of_nonsymm(is_kesz_4, kesz_egy_feladat, F, true);  run_job(fv);
                    is_4_elindult = true;
                }

                // 5

                if (!is_5_elindult && is_kesz_2[1]) {
                    is_kesz_5[1] = true;  fv.math_sub_mul_t_symm_in_nonsymm(is_kesz_5[1], kesz_egy_feladat, K, K, I, Q2, true, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[2], kesz_egy_feladat, L, L, I, Q3, is_priority);  run_job(fv);
                    C.transp(I);
                    is_5_elindult = true;
                }

                // 6

                if (!is_6_elindult && is_kesz_2[2]) {
                    is_kesz_6[2] = true;  fv.math_sub_mul_t_symm_in_nonsymm(is_kesz_6[2], kesz_egy_feladat, P, P, M, Q3, true, is_priority);  run_job(fv);
                    D.transp(M);
                    is_6_elindult = true;
                }

                // készek összeszeddése

                if (!is_3_kesz && is_3_elindult) {
                    if (is_4_elindult && is_kesz_3[1] && is_kesz_3[2])
                        is_3_kesz = true;
                }

                if (!is_5_kesz && is_5_elindult) {
                    if (is_kesz_5[1] && is_kesz_5[2])
                        is_5_kesz = true;
                }

                if (!is_6_kesz && is_6_elindult) {
                    if (is_kesz_6[2])
                        is_6_kesz = true;
                }

                if (is_3_kesz && is_5_kesz && is_6_kesz)
                    is_kesz_1 = true;
            }

            // 7

            Q1.copy(B);
            Q2.transp(G);
            Q3.transp(H);

            // 8

            while (!is_kesz_4) {
                var_egy_esemenyre(is_kesz_4);
            }
            is_kesz_2[1] = true;  fv.math_nmul_t_biztos(is_kesz_2[1], kesz_egy_feladat, J, Q2, F, true);         run_job(fv); // ! 2[1]
            is_kesz_2[0] = true;  fv.math_nmul_t_biztos(is_kesz_2[0], kesz_egy_feladat, B, Q1, F, is_priority);  run_job(fv);
            is_kesz_2[2] = true;  fv.math_nmul_t_biztos(is_kesz_2[2], kesz_egy_feladat, N, Q3, F, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_9_elindult = false;
            bool is_10_elindult = false;
            bool is_11_elindult = false;
            bool is_12_elindult = false;
            bool is_9_kesz = false;
            bool is_11_kesz = false;
            bool is_12_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 9

                if (!is_9_elindult && is_kesz_2[1]) {  // ! 2[1]
                    is_kesz_3[1] = true;  fv.math_sub_mul_t_symm_in_nonsymm(is_kesz_3[1], kesz_egy_feladat, K, K, J, Q2, false, true);  run_job(fv); // ! 3[1]
                    is_kesz_3[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[2], kesz_egy_feladat, L, L, J, Q3, is_priority);  run_job(fv);
                    G.transp(J);
                    is_9_elindult = true;
                }

                // 10

                if (!is_10_elindult && is_9_elindult && is_kesz_3[1]) {  // ! 3[1]
                    is_kesz_4 = true;  fv.math_symm_ninv_of_nonsymm(is_kesz_4, kesz_egy_feladat, K, true);  run_job(fv);
                    is_10_elindult = true;
                }

                // 11

                if (!is_11_elindult && is_kesz_2[0]) { // !
                    is_kesz_5[0] = true;  fv.math_sub_mul_t_symm_in_nonsymm(is_kesz_5[0], kesz_egy_feladat, A, A, B, Q1, true, is_priority);  run_job(fv);
                    is_kesz_5[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[1], kesz_egy_feladat, C, C, B, Q2, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[2], kesz_egy_feladat, D, D, B, Q3, is_priority);  run_job(fv);
                    is_11_elindult = true;
                }

                // 12

                if (!is_12_elindult && is_kesz_2[2]) { // !
                    is_kesz_6[2] = true;  fv.math_sub_mul_t_symm_in_nonsymm(is_kesz_6[2], kesz_egy_feladat, P, P, N, Q3, true, is_priority);  run_job(fv);
                    H.transp(N);
                    is_12_elindult = true;
                }

                // készek összeszeddése

                if (!is_9_kesz && is_9_elindult) {
                    if (is_10_elindult && is_kesz_3[2]) // !
                        is_9_kesz = true;
                }

                if (!is_11_kesz && is_11_elindult) {
                    if (is_kesz_5[0] && is_kesz_5[1] && is_kesz_5[2])
                        is_11_kesz = true;
                }

                if (!is_12_kesz && is_12_elindult) {
                    if (is_kesz_6[2])
                        is_12_kesz = true;
                }

                if (is_9_kesz && is_11_kesz && is_12_kesz)
                    is_kesz_1 = true;
            }

            // 13

            Q1.copy(C);
            Q2.copy(G);
            Q3.transp(L);

            while (!is_kesz_4) {
                var_egy_esemenyre(is_kesz_4);
            }

            // 14

            is_kesz_2[2] = true;  fv.math_nmul_t_biztos(is_kesz_2[2], kesz_egy_feladat, O, Q3, K, true);         run_job(fv); // ! 2[2]
            is_kesz_2[0] = true;  fv.math_nmul_t_biztos(is_kesz_2[0], kesz_egy_feladat, C, Q1, K, is_priority);  run_job(fv);
            is_kesz_2[1] = true;  fv.math_nmul_t_biztos(is_kesz_2[1], kesz_egy_feladat, G, Q2, K, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_15_elindult = false;
            bool is_16_elindult = false;
            bool is_17_elindult = false;
            bool is_18_elindult = false;
            bool is_15_kesz = false;
            bool is_17_kesz = false;
            bool is_18_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 15

                if (!is_15_elindult && is_kesz_2[2]) {  // ! 2[2]
                    is_kesz_3[2] = true;  fv.math_sub_mul_t_symm_in_nonsymm(is_kesz_3[2], kesz_egy_feladat, P, P, O, Q3, false, true);  run_job(fv); // ! 3[2]
                    L.transp(O);
                    is_15_elindult = true;
                }

                // 16

                if (!is_16_elindult && is_15_elindult && is_kesz_3[2]) {  // ! 3[2]
                    is_kesz_4 = true;  fv.math_symm_ninv_of_nonsymm(is_kesz_4, kesz_egy_feladat, P, true);  run_job(fv);
                    is_16_elindult = true;
                }

                // 17

                if (!is_17_elindult && is_kesz_2[0]) { // !
                    is_kesz_5[0] = true;  fv.math_sub_mul_t_symm_in_nonsymm(is_kesz_5[0], kesz_egy_feladat, A, A, C, Q1, true, is_priority);  run_job(fv);
                    is_kesz_5[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[1], kesz_egy_feladat, B, B, C, Q2, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[2], kesz_egy_feladat, D, D, C, Q3, is_priority);  run_job(fv);
                    is_17_elindult = true;
                }

                // 18

                if (!is_18_elindult && is_kesz_2[1]) { // !
                    is_kesz_6[1] = true;  fv.math_sub_mul_t_symm_in_nonsymm(is_kesz_6[1], kesz_egy_feladat, F, F, G, Q2, true, is_priority);  run_job(fv);
                    is_kesz_6[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_6[2], kesz_egy_feladat, H, H, G, Q3, is_priority);  run_job(fv);
                    is_18_elindult = true;
                }

                // készek összeszeddése

                if (!is_15_kesz && is_15_elindult) {
                    if (is_16_elindult) // !
                        is_15_kesz = true;
                }

                if (!is_17_kesz && is_17_elindult) {
                    if (is_kesz_5[0] && is_kesz_5[1] && is_kesz_5[2])
                        is_17_kesz = true;
                }

                if (!is_18_kesz && is_18_elindult) {
                    if (is_kesz_6[1] && is_kesz_6[2])
                        is_18_kesz = true;
                }

                if (is_15_kesz && is_17_kesz && is_18_kesz)
                    is_kesz_1 = true;
            }

            // 19

            S1.copy(D);
            S2.copy(H);
            S3.copy(L);

            while (!is_kesz_4) {
                var_egy_esemenyre(is_kesz_4);
            }

            // 20

            is_kesz_2[0] = true;  fv.math_nmul_t_biztos(is_kesz_2[0], kesz_egy_feladat, D, S1, P, is_priority);  run_job(fv);
            is_kesz_2[1] = true;  fv.math_nmul_t_biztos(is_kesz_2[1], kesz_egy_feladat, H, S2, P, is_priority);  run_job(fv);
            is_kesz_2[2] = true;  fv.math_nmul_t_biztos(is_kesz_2[2], kesz_egy_feladat, L, S3, P, is_priority);  run_job(fv);

            is_kesz_1 = false;
            bool is_21_elindult = false;
            bool is_22_elindult = false;
            bool is_23_elindult = false;
            bool is_21_kesz = false;
            bool is_22_kesz = false;
            bool is_23_kesz = false;

            while (!is_kesz_1) {

                var_egy_esemenyre(is_kesz_1);

                // 21

                if (!is_21_elindult && is_kesz_2[0]) {
                    is_kesz_3[0] = true;  fv.math_sub_mul_t_symm_in_nonsymm(is_kesz_3[0], kesz_egy_feladat, A, A, D, S1, false, is_priority);  run_job(fv);
                    is_kesz_3[1] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[1], kesz_egy_feladat, B, B, D, S2, is_priority);  run_job(fv);
                    is_kesz_3[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_3[2], kesz_egy_feladat, C, C, D, S3, is_priority);  run_job(fv);
                    is_21_elindult = true;
                }

                // 22

                if (!is_22_elindult && is_kesz_2[1]) {
                    is_kesz_5[1] = true;  fv.math_sub_mul_t_symm_in_nonsymm(is_kesz_5[1], kesz_egy_feladat, F, F, H, S2, false, is_priority);  run_job(fv);
                    is_kesz_5[2] = true;  fv.math_sub_mul_t_biztos(is_kesz_5[2], kesz_egy_feladat, G, G, H, S3, is_priority);  run_job(fv);
                    is_22_elindult = true;
                }

                // 23

                if (!is_23_elindult && is_kesz_2[2]) {
                    is_kesz_6[2] = true;  fv.math_sub_mul_t_symm_in_nonsymm(is_kesz_6[2], kesz_egy_feladat, K, K, L, S3, false, is_priority);  run_job(fv);
                    is_23_elindult = true;
                }

                // készek összeszeddése

                if (!is_21_kesz && is_21_elindult) {
                    if (is_kesz_3[0] && is_kesz_3[1] && is_kesz_3[2])
                        is_21_kesz = true;
                }

                if (!is_22_kesz && is_22_elindult) {
                    if (is_kesz_5[1] && is_kesz_5[2])
                        is_22_kesz = true;
                }

                if (!is_23_kesz && is_23_elindult) {
                    if (is_kesz_6[2])
                        is_23_kesz = true;
                }

                if (is_21_kesz && is_22_kesz && is_23_kesz)
                    is_kesz_1 = true;
            }
        }

        p_dest->symmetrize_from_upper();
    }

    //***********************************************************************
    void nonsymm_mul_t(bool is_priority) {
    //***********************************************************************
        if (tipus != pft_mul_t)
            throw hiba("parallel_seged::nonsymm_mul_t", "tipus != pft_mul_t");

        if (szalszam < 2) {
            p_dest->math_mul_t_biztos(*p_src1, *p_src2);
            return;
        }
//printf("nonsymm_mul_t %u %u\n", p_dest->get_row(), p_dest->get_col());
        cuns szal_1 = (szalszam < 6) ?       1        : ((szalszam - 1) / 2);
        cuns szal_2 = (szalszam < 6) ? (szalszam - 1) : (szalszam - szal_1 - 1);

        is_kesz_1 = true;
        hivando_matrixfuggveny<adattipus> fv;
        if (is_vizszintes) {
            if (szal_1 < 2) fv.math_mul_t_biztos(is_kesz_1, kesz_egy_feladat, A, *p_src1, C, is_priority);
            else            fv.math_bigblock_nonsymm_mul_t(is_kesz_1, kesz_egy_feladat, A, szal_1, *p_src1, C, is_priority);
            run_job(fv);
            if (szal_2 < 2) B.math_mul_t_biztos(*p_src1, D);
            else            B.math_bigblock_nonsymm_mul_t(szal_2, *p_src1, D, is_priority);
        }
        else {
            if (szal_1 < 2) fv.math_mul_t_biztos(is_kesz_1, kesz_egy_feladat, A, C, *p_src2, is_priority);
            else            fv.math_bigblock_nonsymm_mul_t(is_kesz_1, kesz_egy_feladat, A, szal_1, C, *p_src2, is_priority);
            run_job(fv);
            if (szal_2 < 2) B.math_mul_t_biztos(D, *p_src2);
            else            B.math_bigblock_nonsymm_mul_t(szal_2, D, *p_src2, is_priority);
        }
        while (!is_kesz_1) {
            var_egy_esemenyre(is_kesz_1);
        }
    }

    //***********************************************************************
    void nonsymm_nmul_t(bool is_priority) {
    //***********************************************************************
        if (tipus != pft_mul_t)
            throw hiba("parallel_seged::nonsymm_nmul_t", "tipus != pft_mul_t");

        if (szalszam < 2) {
            p_dest->math_nmul_t_biztos(*p_src1, *p_src2);
            return;
        }
//printf("nonsymm_nmul_t %u %u\n", p_dest->get_row(), p_dest->get_col());

        cuns szal_1 = (szalszam < 6) ?       1        : ((szalszam - 1) / 2);
        cuns szal_2 = (szalszam < 6) ? (szalszam - 1) : (szalszam - szal_1 - 1);

        is_kesz_1 = true;
        hivando_matrixfuggveny<adattipus> fv;
        if (is_vizszintes) {
            if (szal_1 < 2) fv.math_nmul_t_biztos(is_kesz_1, kesz_egy_feladat, A, *p_src1, C, is_priority);
            else            fv.math_bigblock_nonsymm_nmul_t(is_kesz_1, kesz_egy_feladat, A, szal_1, *p_src1, C, is_priority);
            run_job(fv);
            if (szal_2 < 2) B.math_nmul_t_biztos(*p_src1, D);
            else            B.math_bigblock_nonsymm_nmul_t(szal_2, *p_src1, D, is_priority);
        }
        else {
            if (szal_1 < 2) fv.math_nmul_t_biztos(is_kesz_1, kesz_egy_feladat, A, C, *p_src2, is_priority);
            else            fv.math_bigblock_nonsymm_nmul_t(is_kesz_1, kesz_egy_feladat, A, szal_1, C, *p_src2, is_priority);
            run_job(fv);
            if (szal_2 < 2) B.math_nmul_t_biztos(D, *p_src2);
            else            B.math_bigblock_nonsymm_nmul_t(szal_2, D, *p_src2, is_priority);
        }
        while (!is_kesz_1) {
            var_egy_esemenyre(is_kesz_1);
        }
    }

    //***********************************************************************
    void nonsymm_add_mul_t(bool is_priority) {
    //***********************************************************************
        if (tipus != pft_addsub_mul_t)
            throw hiba("parallel_seged::nonsymm_add_mul_t", "tipus != pft_addsub_mul_t");

        if (szalszam < 2) {
            p_dest->math_add_mul_t_biztos(*p_src1, *p_src2, *p_src3);
            return;
        }
//printf("nonsymm_add_mul_t %u %u\n", p_dest->get_row(), p_dest->get_col());

        cuns szal_1 = (szalszam < 6) ?       1        : ((szalszam - 1) / 2);
        cuns szal_2 = (szalszam < 6) ? (szalszam - 1) : (szalszam - szal_1 - 1);

        is_kesz_1 = true;
        hivando_matrixfuggveny<adattipus> fv;
        if (is_vizszintes) {
            if (szal_1 < 2) fv.math_add_mul_t_biztos(is_kesz_1, kesz_egy_feladat, A, C, *p_src2, E, is_priority);
            else            fv.math_bigblock_nonsymm_add_mul_t(is_kesz_1, kesz_egy_feladat, A, szal_1, C, *p_src2, E, is_priority);
            run_job(fv);
            if (szal_2 < 2) B.math_add_mul_t_biztos(D, *p_src2, F);
            else            B.math_bigblock_nonsymm_add_mul_t(szal_2, D, *p_src2, F, is_priority);
        }
        else {
            if (szal_1 < 2) fv.math_add_mul_t_biztos(is_kesz_1, kesz_egy_feladat, A, C, E, *p_src3, is_priority);
            else            fv.math_bigblock_nonsymm_add_mul_t(is_kesz_1, kesz_egy_feladat, A, szal_1, C, E, *p_src3, is_priority);
            run_job(fv);
            if (szal_2 < 2) B.math_add_mul_t_biztos(D, F, *p_src3);
            else            B.math_bigblock_nonsymm_add_mul_t(szal_2, D, F, *p_src3, is_priority);
        }
        while (!is_kesz_1) {
            var_egy_esemenyre(is_kesz_1);
        }
    }

    //***********************************************************************
    void nonsymm_sub_mul_t(bool is_priority) {
    //***********************************************************************
        if (tipus != pft_addsub_mul_t)
            throw hiba("parallel_seged::nonsymm_sub_mul_t", "tipus != pft_addsub_mul_t");

        if (szalszam < 2) {
            p_dest->math_sub_mul_t_biztos(*p_src1, *p_src2, *p_src3);
            return;
        }
//printf("nonsymm_sub_mul_t %u %u\n", p_dest->get_row(), p_dest->get_col());

        cuns szal_1 = (szalszam < 6) ?       1        : ((szalszam - 1) / 2);
        cuns szal_2 = (szalszam < 6) ? (szalszam - 1) : (szalszam - szal_1 - 1);

        is_kesz_1 = true;
        hivando_matrixfuggveny<adattipus> fv;
        if (is_vizszintes) {
            if (szal_1 < 2) fv.math_sub_mul_t_biztos(is_kesz_1, kesz_egy_feladat, A, C, *p_src2, E, is_priority);
            else            fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_1, kesz_egy_feladat, A, szal_1, C, *p_src2, E, is_priority);
            run_job(fv);
            if (szal_2 < 2) B.math_sub_mul_t_biztos(D, *p_src2, F);
            else            B.math_bigblock_nonsymm_sub_mul_t(szal_2, D, *p_src2, F, is_priority);
        }
        else {
            if (szal_1 < 2) fv.math_sub_mul_t_biztos(is_kesz_1, kesz_egy_feladat, A, C, E, *p_src3, is_priority);
            else            fv.math_bigblock_nonsymm_sub_mul_t(is_kesz_1, kesz_egy_feladat, A, szal_1, C, E, *p_src3, is_priority);
            run_job(fv);
            if (szal_2 < 2) B.math_sub_mul_t_biztos(D, F, *p_src3);
            else            B.math_bigblock_nonsymm_sub_mul_t(szal_2, D, F, *p_src3, is_priority);
        }
        while (!is_kesz_1) {
            var_egy_esemenyre(is_kesz_1);
        }
    }

    //***********************************************************************
    void symm_add_mul_t(bool is_priority) {
    //***********************************************************************
        if (tipus != pft_addsub_mul_t)
            throw hiba("parallel_seged::symm_add_mul_t", "tipus != pft_addsub_mul_t");

        if (szalszam < 2) {
            p_dest->math_add_mul_t_symm(*p_src1, *p_src2, *p_src3);
        }
        else if (szalszam == 2) {
            A.math_add_mul_t_symm(E, I, K);
            C.math_add_mul_t_biztos(G, I, O);
            B.math_add_mul_t_symm(F, J, O);
        }
        else if (szalszam < 5) {
            A.math_add_mul_t_symm(E, I, K);
            C.math_add_mul_t_biztos(G, I, L);
            D.math_add_mul_t_biztos(H, I, M);
            B.math_add_mul_t_symm(F, J, O);
        }
        else {
//printf("symm_add_mul_t %u %u\n", p_dest->get_row(), p_dest->get_col());
            cuns al_szal = (szalszam + 3) / 4;
            A.math_bigblock_symm_add_mul_t(al_szal, E, I, K, is_priority);
            C.math_bigblock_nonsymm_add_mul_t(al_szal, G, I, L, is_priority);
            D.math_bigblock_nonsymm_add_mul_t(al_szal, H, I, M, is_priority);
            B.math_bigblock_symm_add_mul_t(al_szal, F, J, O, is_priority);
        }
    }

    //***********************************************************************
    void symm_in_nosymm_sub_mul_t(bool is_symmetrize_needed, bool is_priority) {
    //***********************************************************************
        if (tipus != pft_addsub_mul_t)
            throw hiba("parallel_seged::symm_sub_mul_t", "tipus != pft_addsub_mul_t");

        if (szalszam < 2) {
            p_dest->math_sub_mul_t_symm_in_nonsymm(*p_src1, *p_src2, *p_src3, false);
        }
        else if (szalszam == 2) {
            A.math_sub_mul_t_symm_in_nonsymm(E, I, K, false);
            C.math_sub_mul_t_biztos(G, I, O);
            B.math_sub_mul_t_symm_in_nonsymm(F, J, O, false);
        }
        else if (szalszam < 5) {
            A.math_sub_mul_t_symm_in_nonsymm(E, I, K, false);
            C.math_sub_mul_t_biztos(G, I, L);
            D.math_sub_mul_t_biztos(H, I, M);
            B.math_sub_mul_t_symm_in_nonsymm(F, J, O, false);
        }
        else {
//printf("symm_in_nosymm_sub_mul_t %u %u\n", p_dest->get_row(), p_dest->get_col());
            cuns al_szal = (szalszam + 3) / 4;
            A.math_bigblock_symm_in_nosymm_sub_mul_t(al_szal, E, I, K, false, is_priority);
            C.math_bigblock_nonsymm_sub_mul_t(al_szal, G, I, L, is_priority);
            D.math_bigblock_nonsymm_sub_mul_t(al_szal, H, I, M, is_priority);
            B.math_bigblock_symm_in_nosymm_sub_mul_t(al_szal, F, J, O, false, is_priority);
        }
        if (is_symmetrize_needed)
            p_dest->symmetrize_from_upper();
    }

    //***********************************************************************
    bool check_inv(matrix<adattipus> & d) const {
    //***********************************************************************
        return tipus == pft_inv && p_dest == &d;
    }

    //***********************************************************************
    bool check_mul_t(matrix<adattipus> & d, const matrix<adattipus> & c1, const matrix<adattipus> & c2) const {
    //***********************************************************************
        return tipus == pft_mul_t && p_dest == &d && p_src1 == &c1 && p_src2 == &c2;
    }

    //***********************************************************************
    bool check_addsub_mul_t(matrix<adattipus> & d, const matrix<adattipus> & c1, const matrix<adattipus> & c2,
        const matrix<adattipus> & c3) const {
    //***********************************************************************
        return tipus == pft_addsub_mul_t && p_dest == &d && p_src1 == &c1 && p_src2 == &c2 && p_src3 == &c3;
    }

};


//***********************************************************************
template<typename adattipus> matrix<adattipus>::~matrix() {
//***********************************************************************
    delete p_parallel_seged;
};


//***********************************************************************
template<typename adattipus> void matrix<adattipus>::math_bigblock_ninv_np(uns osztasszam, bool is_priority) {
//***********************************************************************
    if (p_parallel_seged == nullptr || !p_parallel_seged->check_inv(*this)) {
        delete p_parallel_seged;
        p_parallel_seged = new parallel_seged<adattipus>(osztasszam, *this, true);
    }
    p_parallel_seged->nonsymm_in_nonsymm_ninv(is_priority);
};


//***********************************************************************
template<typename adattipus> void matrix<adattipus>::math_bigblock_symm_in_nonsymm_ninv_np(uns osztasszam, bool is_priority) {
//***********************************************************************
    if (p_parallel_seged == nullptr || !p_parallel_seged->check_inv(*this)) {
        delete p_parallel_seged;
        p_parallel_seged = new parallel_seged<adattipus>(osztasszam, *this, false);
    }
    p_parallel_seged->symm_in_nonsymm_ninv(is_priority);
};


//***********************************************************************
template<typename adattipus> void matrix<adattipus>::math_bigblock_nonsymm_mul_t(uns osztasszam, const matrix & src1, 
    const matrix & src2, bool is_priority) {
//***********************************************************************
    if (row == 0 || col == 0)
        return;
    if (p_parallel_seged == nullptr || !p_parallel_seged->check_mul_t(*this, src1, src2)) {
        delete p_parallel_seged;
        p_parallel_seged = new parallel_seged<adattipus>(osztasszam, *this, src1, src2);
    }
    p_parallel_seged->nonsymm_mul_t(is_priority);
};


//***********************************************************************
template<typename adattipus> void matrix<adattipus>::math_bigblock_nonsymm_nmul_t(uns osztasszam, const matrix & src1, 
    const matrix & src2, bool is_priority) {
//***********************************************************************
    if (row == 0 || col == 0)
        return;
    if (p_parallel_seged == nullptr || !p_parallel_seged->check_mul_t(*this, src1, src2)) {
        delete p_parallel_seged;
        p_parallel_seged = new parallel_seged<adattipus>(osztasszam, *this, src1, src2);
    }
    p_parallel_seged->nonsymm_nmul_t(is_priority);
};


//***********************************************************************
template<typename adattipus> void matrix<adattipus>
    ::math_bigblock_nonsymm_add_mul_t(uns osztasszam, const matrix & src1, const matrix & src2, const matrix & src3, bool is_priority) {
//***********************************************************************
    if (row == 0 || col == 0)
        return;
    if (p_parallel_seged == nullptr || !p_parallel_seged->check_addsub_mul_t(*this, src1, src2, src3)) {
        delete p_parallel_seged;
        p_parallel_seged = new parallel_seged<adattipus>(osztasszam, false, false, *this, src1, src2, src3);
    }
    p_parallel_seged->nonsymm_add_mul_t(is_priority);
};


//***********************************************************************
template<typename adattipus> void matrix<adattipus>
    ::math_bigblock_nonsymm_sub_mul_t(uns osztasszam, const matrix & src1, const matrix & src2, const matrix & src3, bool is_priority) {
//***********************************************************************
    if (row == 0 || col == 0)
        return;
    if (p_parallel_seged == nullptr || !p_parallel_seged->check_addsub_mul_t(*this, src1, src2, src3)) {
        delete p_parallel_seged;
        p_parallel_seged = new parallel_seged<adattipus>(osztasszam, false, false, *this, src1, src2, src3);
    }
    p_parallel_seged->nonsymm_sub_mul_t(is_priority);
};


//***********************************************************************
template<typename adattipus> void matrix<adattipus>
    ::math_bigblock_symm_add_mul_t(uns osztasszam, const matrix & src1, const matrix & src2, const matrix & src3, bool is_priority) {
//***********************************************************************
    if (row == 0 || col == 0)
        return;
    if (p_parallel_seged == nullptr || !p_parallel_seged->check_addsub_mul_t(*this, src1, src2, src3)) {
        delete p_parallel_seged;
        p_parallel_seged = new parallel_seged<adattipus>(osztasszam, true, false, *this, src1, src2, src3);
    }
    p_parallel_seged->symm_add_mul_t(is_priority);
};


//***********************************************************************
template<typename adattipus> void matrix<adattipus>
    ::math_bigblock_symm_in_nosymm_sub_mul_t(uns osztasszam, const matrix & src1, const matrix & src2, const matrix & src3, 
        bool is_symmetrize_needed, bool is_priority) {
//***********************************************************************
    if (row == 0 || col == 0)
        return;
    if (p_parallel_seged == nullptr || !p_parallel_seged->check_addsub_mul_t(*this, src1, src2, src3)) {
        delete p_parallel_seged;
        p_parallel_seged = new parallel_seged<adattipus>(osztasszam, false, true, *this, src1, src2, src3);
    }
    p_parallel_seged->symm_in_nosymm_sub_mul_t(is_symmetrize_needed, is_priority);
};


//***********************************************************************
template<typename adattipus> 
void hivando_matrixfuggveny<adattipus>::run_job() {
//***********************************************************************
    switch (tipus) {
        case fvt_math_add_mul_t_biztos:
            p_dest->math_add_mul_t_biztos(*p_src1, *p_src2, *p_src3);
            break;
        case fvt_math_add_mul_t_symm:
            p_dest->math_add_mul_t_symm(*p_src1, *p_src2, *p_src3);
            break;
        case fvt_math_bigblock_ninv_np:
            p_dest->math_bigblock_ninv_np(s_osztasszam, s_is_priority);
            break;
        case fvt_math_bigblock_nonsymm_add_mul_t:
            p_dest->math_bigblock_nonsymm_add_mul_t(s_osztasszam, *p_src1, *p_src2, *p_src3, s_is_priority);
            break;
        case fvt_math_bigblock_nonsymm_mul_t:
            p_dest->math_bigblock_nonsymm_mul_t(s_osztasszam, *p_src1, *p_src2, s_is_priority);
            break;
        case fvt_math_bigblock_nonsymm_nmul_t:
            p_dest->math_bigblock_nonsymm_nmul_t(s_osztasszam, *p_src1, *p_src2, s_is_priority);
            break;
        case fvt_math_bigblock_nonsymm_sub_mul_t:
            p_dest->math_bigblock_nonsymm_sub_mul_t(s_osztasszam, *p_src1, *p_src2, *p_src3, s_is_priority);
            break;
        case fvt_math_bigblock_symm_add_mul_t:
            p_dest->math_bigblock_symm_add_mul_t(s_osztasszam, *p_src1, *p_src2, *p_src3, s_is_priority);
            break;
        case fvt_math_bigblock_symm_in_nonsymm_ninv_np:
            p_dest->math_bigblock_symm_in_nonsymm_ninv_np(s_osztasszam, s_is_priority);
            break;
        case fvt_math_bigblock_symm_in_nosymm_sub_mul_t:
            p_dest->math_bigblock_symm_in_nosymm_sub_mul_t(s_osztasszam, *p_src1, *p_src2, *p_src3, s_is_symmetrize_needed, s_is_priority);
            break;
        case fvt_math_mul_t_biztos:
            p_dest->math_mul_t_biztos(*p_src1, *p_src2);
            break;
        case fvt_math_ninv_np:
            p_dest->math_ninv_np();
            break;
        case fvt_math_nmul_t_biztos:
            p_dest->math_nmul_t_biztos(*p_src1, *p_src2);
            break;
        case fvt_math_sub_mul_t_biztos:
            p_dest->math_sub_mul_t_biztos(*p_src1, *p_src2, *p_src3);
            break;
        case fvt_math_sub_mul_t_symm_in_nonsymm:
            p_dest->math_sub_mul_t_symm_in_nonsymm(*p_src1, *p_src2, *p_src3, s_is_symmetrize_needed);
            break;
        case fvt_math_symm_ninv_of_nonsymm:
            p_dest->math_symm_ninv_of_nonsymm();
            break;
        default:
            throw hiba("hivando_matrixfuggveny::run_job", "unknown job type");
    }
}


}

#endif