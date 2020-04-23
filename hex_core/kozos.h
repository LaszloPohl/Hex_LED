//***********************************************************************
// közös típusok, deklarációk header
// Creation date:  2018. 06. 27.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef KOZOS_HEADER
#define	KOZOS_HEADER
//***********************************************************************

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
//#include <sstream>
#include <cmath>
#include <complex>
#include "hiba.h"


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
typedef unsigned uns;
typedef const unsigned cuns;
typedef double rvt; // real value type
typedef uns meret_t;
enum nonlin_tipus { nt_konstans, nt_lin, nt_szakszok, nt_fazisvalto, nt_exp, nt_mizs, nt_diode_1, nt_erno, nt_polinom, nt_jani_1, nt_inv }; // nt_konstans, nt_lin, nt_szakszok, nt_fazisvalto, nt_exp, nt_mizs, nt_diode_1, nt_erno, nt_inv
enum mezo_tipus { mt_elektromos, mt_termikus, mt_elektrotermikus };
enum nemlin_megoldo_tipus{ nmt_klasszik_iteracio, nmt_el_th_newton };  // ha két fa van, akkor külön newton, egyébként egyesített newton
enum peremtipus { pt_open, pt_u, pt_t, pt_htc, pt_thtc };
enum analizis_lepes_tipus { alt_dc, alt_trans, alt_ac };
enum gerjesztes_tipus { gt_el_none, gt_th_none, gt_U, gt_I, gt_T, gt_P };
enum eredmeny_tipus { et_c_pontprobe, et_f_pontprobe, et_c_map, et_f_map, et_sum_value };
enum mit_ment_tipus { mm_UT, mm_IP, mm_RL };
enum analizis_beallitas_tipus { abt_I0, abt_max_error, abt_max_iter, abt_del_all_prev, abt_del_all_gerj, abt_del_fa, abt_tamb,
    abt_change_bn, abt_change_bv, abt_change_jn, abt_change_mn, abt_ct }; // abt_del_all_prev: az elõzõ lépésben számolt össze U/T/I/P 0-zása
    // bn: boundary number, bv: boundary value, jn: junction number, mn: material number
enum cellatipus { ct_none, ct_faces_cella, ct_th6 };
enum face_tipus { ft_none, ft_csatlakozo, ft_normal_perem, ft_spec_perem, ft_centroid }; // centroid = kivezetett középponti
enum peremirany_tipus{ pit_west, pit_east, pit_south, pit_north, pit_bottom, pit_top, pit_none };
enum fa_adat_tipus { fat_double, fat_float, fat_double_ended_float };
enum szal_tipus { szt_nem_csinal_semmit, szt_cellafeldolgozo_klaszter, szt_lepo_klaszter, szt_fw_klaszter, szt_fw_egyedi, szt_bw_egyedi, szt_bw_klaszter, szt_matrix, szt_belso_mentes, szt_fajlba_mentes, szt_N }; // szt_N: hány száltípus van, tömbmérethez
enum szal_altipus { szat_semmi, szat_pre, szat_post, szat_elore, szat_vissza, szat_megtart_akt, szat_ment_akt, szat_kiind_akt, szat_kiind_akt_es_elore }; // szat_semmi, szat_pre, szat_post, szat_elore, szat_vissza, szat_megtart_akt, szat_ment_akt, szat_kiind_akt, szat_kiind_akt_es_elore
enum illesztes_tipusa { it_none, it_lin, it_strong }; // it_none, it_lin, it_strong
#define absT 273.15
#define g_max 1.0e+20
#define g_min 1.0e-20
//***********************************************************************
template<typename T>inline T sqr(const T & a) { return a*a; }
//***********************************************************************
template<typename T>inline T replusz(const T & a, const T & b) { return a * b / (a + b); }
//***********************************************************************
//***********************************************************************
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#define fpu_error(x) (isinf(x) || isnan(x))
#define FPU_ERROR_CHECK
#ifdef FPU_ERROR_CHECK
#define chech_fpu_error(x) if(fpu_error(x)){ printf("\n\n%s\n", __func__); getchar(); }
#define chech_fpu_error_name(x,n) if(fpu_error(x)){ printf("\n\n%s: %s\n", __func__, n); getchar(); }
#else
#define chech_fpu_error(x) 
#define chech_fpu_error_name(x,n)
#endif
//***********************************************************************


//***********************************************************************
struct szaltipus_jellemzoi {
//***********************************************************************
    szal_tipus valid;   // ha nem szal_tipus::szt_N, ki van töltve, Para_engine::set_szalszam ellenõrzi
    bool is_egyszalu;   // ha true, 1 szál, ha false, N szálú
    bool is_azonnal_torolheto; // ha false, számítási eredményt õriz, csak ennek feldolgozása után törölhetõ, egyébként azonnal, ha kész
    bool is_felteteles; // true esetén csak akkor indítható, ha teljesül a feltétel, egyébként várnia kell
    szal_tipus ertesitendo_felteteles; // ha kész ez a feladat, akkor kell-e értesítenie egy másik típusú szálcsoportot, hogy nézze meg, indítható-e
    szaltipus_jellemzoi(szal_tipus valid = szal_tipus::szt_N, bool is_egyszalu = false, bool is_azonnal_torolheto = false,
        bool is_felteteles = false, szal_tipus ertesitendo_felteteles = szal_tipus::szt_N)
        :valid{ valid }, is_egyszalu{ is_egyszalu }, is_azonnal_torolheto{ is_azonnal_torolheto }, is_felteteles{ is_felteteles },
         ertesitendo_felteteles{ ertesitendo_felteteles } {}
};
//***********************************************************************
extern const szaltipus_jellemzoi szaltipus_jellemzoi_tomb[szal_tipus::szt_N];
//***********************************************************************


//***********************************************************************
struct klaszter_tartomany {
//***********************************************************************
    uns klaszter_kezdoindex, klaszter_utolso_index;
    klaszter_tartomany() :klaszter_kezdoindex{ 0 }, klaszter_utolso_index{ 0 } {}
};
//***********************************************************************


//***********************************************************************
struct Itrio{
//***********************************************************************
    rvt I, dI_per_dU, dI_per_dT;
};
//***********************************************************************


//***********************************************************************
void most(const ::std::string & mi_tortent);
void kiirt_ido_forma(bool masodpercben); // ha true, akkor nem mutat ms-ban
void esemenyek_kiirasa();
void utolso_esemeny_kiirasa();
//***********************************************************************


//***********************************************************************
class akt_lepes {
//***********************************************************************
    static ::std::string aktualis_lepes_neve;
    static uns aktualis_lepes_szazalek, ossz_lepesszam, aktualis_lepes_szama;
    static bool akt_lepes_kiir;
public:
    //***********************************************************************
    static void set_aktualis_lepes_neve(const ::std::string & nev) {
   //***********************************************************************
        aktualis_lepes_neve = nev;
        aktualis_lepes_szazalek = 0;
    }

    //***********************************************************************
    static void get_aktualis_lepes_szazalek(const ::std::string * & nev, uns & szazalek, uns & akt_lepes, uns & ossz_lepes) {
        //***********************************************************************
        nev = &aktualis_lepes_neve;
        szazalek = aktualis_lepes_szazalek;
        akt_lepes = aktualis_lepes_szama;
        ossz_lepes = ossz_lepesszam;
    }

    //***********************************************************************
    static void set_aktualis_lepes_szazalek(uns szazalek) { aktualis_lepes_szazalek = szazalek; }
    //***********************************************************************
    static void set_ossz_lepesszam(uns n) { ossz_lepesszam = n; }
    //***********************************************************************
    static void inc_ossz_lepesszam(uns n) { ossz_lepesszam++; }
    //***********************************************************************
    static void inc_akt_lepesszam() { aktualis_lepes_szama++; }
    //***********************************************************************
    static void set_akt_lepes_kiir(bool is_kiir) { akt_lepes_kiir = is_kiir; }
    //***********************************************************************
    static bool get_akt_lepes_kiir() { return akt_lepes_kiir; }
    //***********************************************************************
};
//***********************************************************************


}

#endif
