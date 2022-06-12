//***********************************************************************
// cella oszt�ly �s face oszt�lyok header
// Creation date:  2018. 08. 10.
// Creator:        Pohl L�szl�
//***********************************************************************


//***********************************************************************
#ifndef CELLA_ES_FACE_HEADER
#define	CELLA_ES_FACE_HEADER
//***********************************************************************


//***********************************************************************
#include "matrix.hpp"
#include "akt_sim_adatai.h"
#include "emlekezo_tar.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
extern akt_sim_adatai akt_sim, prev_sim;
//***********************************************************************


//***********************************************************************
class os_cella {
//***********************************************************************
protected:
    bool is_besugarzo; // bool, ha van olyan cella, amelyiknek a f�ny�t ez valamilyen ar�nyban h�v� alak�tja
    bool is_elektrotermikus_szimulacio; // elektrotermikus szimul�ci� eset�n lehet olyan cella, ahol nem elektrotermikus szimul�ci� t�rt�nik, ez a cell�ra vonatkozik!
    bool is_change_adm, is_change_j; // bool. bw t�rli. kell-e a cell�ban redukci� ill fw? Ha false, akkor a fa nem friss�ti mag�t!
    bool is_2_fa; // 1 vagy 2 fa van-e, azaz k�l�n kell-e szedni a k�t teret
    bool is_lin; // ha kapacit�s van, akkor is lin, ha lin
    uns A_1, A_2; // YA m�trixok m�rete
    uns B_1, B_2; // YB m�trixok m�rete
public:
    //***********************************************************************
    rvt A_top; // F�nysz�m�t�shoz sz�ks�ges hack. Ha van top oldal defini�lva, akkor az inicializ�ci� sor�n be kell �ll�tani. Default 1!
    //***********************************************************************
    matrix<rvt> dc_yred_1, dc_yred_2;
    vektor<rvt> dc_jred_1, dc_jred_2;
    vektor<rvt> dc_UA_1, dc_UA_2, dc_IA_1, dc_IA_2;
    vektor<uns> dc_iter_csp_index;
    matrix<::std::complex<rvt> > ac_yred_1, ac_yred_2;
    vektor<::std::complex<rvt> > ac_jred_1, ac_jred_2;
    vektor<::std::complex<rvt> > ac_UA_1, ac_UA_2, ac_IA_1, ac_IA_2;
    //***********************************************************************
public:
    uns cella_index; // debug c�lb�l
    const adat_cella * p_adatcella;
    const adat_anyag * p_cella_anyag; // nem minden cellat�pushoz van hozz�rendelt anyag, pl. komplex_cella, ekkor nullptr
    bool is_el, is_th; // ha van elektromos/termikus face, akkor igaz
    bool is_iter_csp_frissitendo;

    //***********************************************************************
    struct cella_emlek {
    //***********************************************************************
        rvt sum_rad, sum_lum, sum_diss, A_LED;  // a cella teljes sug�rzott teljes�tm�nye/f�ny�rama/f�nysug�rz� fel�lete
        rvt sum_blue, sum_yellow;               // a cell�b�l kil�p� �sszes k�k ill. s�rga f�ny
        rvt sum_new_yellow;                     // a kil�p� �sszes s�rg�b�l ennyi az �j (s�rga utakhoz)
        rvt sum_blue_out, sum_yellow_out, sum_diss_out; // a cell�ban v�gz�d� sugarakb�l mennyi j�n ki (ezek �sszega addja a teljes k�k ill s�rga kimen� f�nyt a LED-b�l)
        ertek_t H;                              // Az aktu�lis HB (�s HA = megtartando) hiszter�zis, valamint a meredeks�ge
    };
    //***********************************************************************

    //***********************************************************************
    emlek<cella_emlek> emlekek;
    //***********************************************************************

    //***********************************************************************
    void pre_init(uns Cella_Index) {
    //***********************************************************************
        cella_index = Cella_Index;
        p_adatcella = &akt_sim.p_akt_sim->cellak[cella_index];
        const adat_cella & adatcella = *p_adatcella;
        is_el = adatcella.mezotipus == mt_elektromos || adatcella.mezotipus == mt_elektrotermikus;
        is_th = adatcella.mezotipus == mt_termikus || adatcella.mezotipus == mt_elektrotermikus;
        is_besugarzo = adatcella.besugarzo_cella.size() > 0;
        is_elektrotermikus_szimulacio = is_el && is_th;
        pre_init_egyedi();
    }

    //***********************************************************************
    os_cella() :is_lin{ false }, cella_index{ 0 }, is_el{ false }, is_th{ false }, p_adatcella{ nullptr },
        is_elektrotermikus_szimulacio{ false }, is_besugarzo{ false }, p_cella_anyag{ nullptr }, 
        A_1{ 0 }, A_2{ 0 }, B_1{ 0 }, B_2{ 0 }, is_2_fa{ false }, is_change_adm{ true }, is_change_j{ true },
        A_top{ 1 }, is_iter_csp_frissitendo{ false } {}
    void reset_H(); // �j szimul�ci� �s l�nyeges v�ltoz�s eset�n
    bool is_adm_changed()const { return is_change_adm; }
    bool is_j_changed()const { return is_change_j; }
    //***********************************************************************

    //***********************************************************************
    virtual void debug_write(::std::ofstream & fs) const = 0;
    virtual cellatipus get_tipus()const = 0;
    virtual void pre_init_egyedi() = 0;
    virtual bool is_nonsymmetrical()const = 0;
    //***********************************************************************

    //***********************************************************************
    //                             *****  DC  *****
    //***********************************************************************
    virtual void bw_dc(uns melyik_fa) = 0;
    virtual void fw_dc(uns melyik_fa) = 0;
    virtual void init_dc() = 0; // mem�riafoglal�s
    virtual void resize_matrixok_dc() = 0; // kisz�m�tja a m�trixok m�ret�t, �s lefoglalja azokat, �j fa eset�n is
    virtual void halozat_foglalas_cellaban_dc() = 0; // is_lin-t is ez �ll�tja, eredetileg face_foglalas volt
    virtual void gerj_update_dc() = 0;
    virtual void set_is_szimm() = 0; // A h�l�zat (face-ek) foglal�sa ut�n
    virtual void update_dt_dc() = 0; // csak akkor kell k�l�n megh�vni, ha iter�ci� k�zben v�ltozik a dt
    virtual void peremfeltetel_update_dc() = 0;
    virtual void del_all_prev_dc() = 0;
    virtual void update_for_uj_lepes_dc() = 0; // pl. kapacit�sok fesz�lts�ge�nek friss�t�se
    virtual void update_iter_csp() = 0;
    // void update_for_uj_lepes_ac(); // AC-n�l nincs bels� iter�ci�
    //***********************************************************************
    // NR = Newton-Raphson
    // SA = Successive Approximation
    virtual rvt  cella_klaszter_2_NR_UT_dc(rvt & Uc, rvt & Tc, rvt & T_hiba) = 0; // Ui = kiindulo_Ui + dUi, vissza: szumma dUi hiba �s a cella Uc �s Tc �rt�ke
    virtual rvt  cella_klaszter_2_SA_UT_dc(rvt & T_hiba) = 0; // csak hib�t sz�mol, vissza: szumma Ui-el�z�_Ui
    virtual void cella_klaszter_3_parameterek_frissitese_dc() = 0; // h�m�rs�kletf�gg� vezet�sek, junction vezet�se �gfesz�lts�gb�l
    virtual void cella_klaszter_4_agaramok_szamitasa_dc() = 0; // disszip�ci� is
    virtual rvt  cella_klaszter_5_NR_hibaaramok_szamitasa_dc() = 0; // J felt�lt�se, vissza: cellak�z�pponti I+P hiba�ram
    virtual void cella_klaszter_5_SA_fill_J_dc() = 0; //
    virtual void fw_klaszter_1_update_Jakobi_dc() = 0; // NR eset�ben a Jakobi, SA eset�ben az admittancia m�trixot t�lti fel a nemline�ris �rt�kekkel
    //***********************************************************************

    //***********************************************************************
    //                             *****  AC  *****
    //***********************************************************************
    virtual void bw_ac(uns melyik_fa) = 0;
    virtual void fw_ac(uns melyik_fa) = 0;
    virtual void init_ac() = 0; // mem�riafoglal�s
    virtual void face_foglalas_ac() = 0;
    virtual void gerj_update_ac() = 0;
    virtual void peremfeltetel_update_ac() = 0;
    virtual void del_all_prev_ac() = 0;
    virtual void update_ac(bool force) = 0;
    //***********************************************************************
};


//***********************************************************************
class th6_cella : public os_cella {
//***********************************************************************
private:
public:
    void debug_write(::std::ofstream & fs) const override { TODO("th6_cella::debug_write"); }
    cellatipus get_tipus()const override { return ct_th6; }
    void pre_init_egyedi() override {}
    bool is_nonsymmetrical()const override { return false; }
    //***********************************************************************
    //                             *****  DC  *****
    //***********************************************************************
    void init_dc() override { resize_matrixok_dc(); }
    void resize_matrixok_dc() override;
    void halozat_foglalas_cellaban_dc() override;
    void bw_dc(uns melyik_fa) override { TODO("th6_cella::bw_dc"); }
    void fw_dc(uns melyik_fa) override { TODO("th6_cella::fw_dc"); }
    void gerj_update_dc() override { TODO("th6_cella::gerj_update_dc"); }
    void set_is_szimm() override { TODO("th6_cella::set_is_szimm"); }
    void update_dt_dc() override { TODO("th6_cella::update_dt_dc"); }
    void peremfeltetel_update_dc() override { TODO("th6_cella::peremfeltetel_update_dc"); }
    void del_all_prev_dc() override { TODO("th6_cella::del_all_prev_dc"); }
    void update_for_uj_lepes_dc() override { TODO("th6_cella::update_for_uj_lepes_dc"); }
    //***********************************************************************
    rvt  cella_klaszter_2_NR_UT_dc(rvt & Uc, rvt & Tc, rvt & T_hiba) override { TODO("th6_cella::cella_klaszter_2_NR_UT_dc"); }
    rvt  cella_klaszter_2_SA_UT_dc(rvt & T_hiba) override { TODO("th6_cella::cella_klaszter_2_SA_UT_dc"); }
    void cella_klaszter_3_parameterek_frissitese_dc() override { TODO("th6_cella::cella_klaszter_3_parameterek_frissitese_dc"); }
    void cella_klaszter_4_agaramok_szamitasa_dc() override { TODO("th6_cella::cella_klaszter_4_agaramok_szamitasa_dc"); }
    rvt  cella_klaszter_5_NR_hibaaramok_szamitasa_dc() override { TODO("th6_cella::cella_klaszter_5_NR_hibaaramok_szamitasa_dc"); }
    void cella_klaszter_5_SA_fill_J_dc() override { TODO("th6_cella::cella_klaszter_5_SA_fill_J_dc"); }
    void fw_klaszter_1_update_Jakobi_dc() override { TODO("th6_cella::fw_klaszter_1_update_Jakobi_dc"); }
    //***********************************************************************

    //***********************************************************************
    //                             *****  AC  *****
    //***********************************************************************
    void bw_ac(uns melyik_fa) override { TODO("th6_cella::bw_ac"); }
    void fw_ac(uns melyik_fa) override { TODO("th6_cella::fw_ac"); }
    void init_ac() override { TODO("th6_cella::init_ac"); }
    void face_foglalas_ac() override { TODO("th6_cella::face_foglalas_ac"); }
    void gerj_update_ac() override { TODO("th6_cella::gerj_update_ac"); }
    void peremfeltetel_update_ac() override { TODO("th6_cella::peremfeltetel_update_ac"); }
    void del_all_prev_ac() override { TODO("th6_cella::del_all_prev_ac"); }
    void update_ac(bool force) override { TODO("th6_cella::update_ac"); }
    //***********************************************************************

    //***********************************************************************
    th6_cella::th6_cella() {}
    //***********************************************************************
};


//***********************************************************************
class os_face_core_dc;
class os_face_core_ac;
extern vektor<rvt> csatlakozo_aramok_dc; // a 0 index� dummy!
extern vektor<iter_csomopont> iter_csomopontok_dc; // a 0 index� dummy!
//***********************************************************************


//***********************************************************************
struct face_tok_dc {
//***********************************************************************

    //***********************************************************************
    struct face_fok_dc_emlek {
    //***********************************************************************
        rvt UT, IP; // a face cser�l�dhet, ez�rt itt kell t�rolni a kor�bbi �rt�keket (ez a teljes UT vagy IP)
        rvt rad, lum;
    };

    os_face_core_dc * p_core;
    face_tok_dc * p_parja; // elektrotermikusn�l a m�sik
    //rvt dUT, dIP; // Newton-Raphson eset�n ide ker�l a bw-b�l visszakapott eredm�ny
    //egytarolo<emlekezo::tar<rvt> > dissz, rad, lum; // ezek kellenek majd
    rvt *p_Yii, *p_Yie, *p_Yit, *p_Yei, *p_Yti, *p_Ji, *p_UTi, *p_IPi;
    emlek<face_fok_dc_emlek> emlekek;
    //***********************************************************************
    face_tok_dc() : p_core{ nullptr }, p_parja{ nullptr }, p_Yii{ nullptr }, p_Yie{ nullptr },
        p_Yit{ nullptr }, p_Yei{ nullptr }, p_Yti{ nullptr }, p_Ji{ nullptr }, 
        p_UTi{ nullptr }, p_IPi{ nullptr }/*, dUT{ rvt() }, dIP{ rvt() }*/
        {}
    //***********************************************************************
    ~face_tok_dc();
    //***********************************************************************
    void debug_write(::std::ofstream & fs) const {
    //***********************************************************************
        fs << "UT_akt=" << ::std::left << ::std::setfill(' ') << ::std::setw(12) << emlekek.get_akt().UT << " UT_ki=" << ::std::setw(12) << emlekek.get_kiindulo().UT << " UT_meg=" << ::std::setw(12) << emlekek.get_megtartando().UT;
        fs << " IP_akt=" << ::std::setw(12) << emlekek.get_akt().IP << " IP_ki=" << ::std::setw(12) << emlekek.get_kiindulo().IP << " IP_meg=" << ::std::setw(12) << emlekek.get_megtartando().IP;
    }
};


//***********************************************************************
struct ac_tok_bele {
//***********************************************************************

    //***********************************************************************
    struct face_fok_ac_emlek {
    //***********************************************************************
        ::std::complex<rvt> UT, IP;
    };

    emlek<face_fok_ac_emlek> emlekek;
    //egytarolo<emlek<::std::complex<rvt> > > dissz, rad, lum;
    ::std::complex<rvt> Yi, Yc, Jc, Jx;
    ac_tok_bele() : Yi{ rvt() }, Yc{ rvt() }, Jc{ rvt() }, Jx{ rvt() } {}
};


//***********************************************************************
struct face_tok_ac {
//***********************************************************************
    os_face_core_ac * p_core;
    face_tok_ac * p_parja; // elektrotermikusn�l a m�sik
    egytarolo<ac_tok_bele> bele;
    face_tok_ac() : p_core{ nullptr }, p_parja{ nullptr } {}
    ~face_tok_ac();
};

//***********************************************************************
class dc_kozepponti_face;
class szal_feladat_adatai;
//***********************************************************************
class faces_cella : public os_cella {
//***********************************************************************
    friend class dc_kozepponti_face;
    //friend void run_cellafeldolgozo_klaszter(szal_feladat_adatai *);
    //***********************************************************************
    bool is_nemszimm; // ha nemszimmetrikusak a m�trixok, true
    bool is_NR_nincs_csatolt_nemszimm; // Van-e olyan face, amelyn�l NR nem csatolt (egyteres) szimul�ci�n�l nemszimmetrikus az Y
    bool is_NR_van_csatolt_nemszimm; // Van-e olyan face, amelyn�l NR csatolt (k�tteres) szimul�ci�n�l nemszimmetrikus az Y
    bool is_perem_face; // van-e perem_face a cell�ban
    enum { ft_semmi, ft_el, ft_th, ft_elth_nincs_csatolt, ft_elth_van_csatolt } feltoltes_tipus; // mely m�trixelemeket kell felt�lteni
public: // Hogy a face-ek hozz�f�rhessenek

    rvt sum_disszipacio_dc, sum_ddissz_per_dT_dc; // ide gy�jtj�k, mert t�bb elektromos cella is egy termikusra adhatja.
    rvt *p_sum_disszipacio_dc, *p_sum_ddissz_per_dT_dc; // ha k�l�nv�lasztott elektrotermikus van, akkor ez a m�sik cell�ra mutat
    rvt konst_Yee_dc, konst_Yet_dc, konst_Yte_dc, konst_Ytt_dc, konst_Je_dc, konst_Jt_dc; // Az adott �rt�k konstans r�sze, a face::update_uj_lepeshez() �rja +=-vel, ha face foglal�s vagy gerj update volt
    rvt *p_Yee_dc, *p_Yet_dc, *p_Yte_dc, *p_Ytt_dc, *p_Je_dc, *p_Jt_dc; // �gy kevesebb pointer kell a face tokokba
    face_tok_dc el_center_face_dc, th_center_face_dc;
    vektor<face_tok_dc> faces_dc;   // face_tok_dc vektor, a 0 index� dummy
    rvt A_cella; // A keresztface init-ek t�ltik fel
private:
    faces_cella * kapcsolodo_cella; // van, ahol felt�telezz�k, hogy ha ez nem nullptr, akkor a m�sik t�rbeli p�rj�ra mutat, teh�t m�s c�lra ne haszn�ld!
    uns el_centroid_index_bemeneti, th_centroid_index_bemeneti; // 1-t�l indexel�dik, a 0 dummy
    uns el_centroid_index_matrix, th_centroid_index_matrix; // vigy�zz, 0-t�l indexel�dik, NEM dummy
    //uns P_1, P_2; // premfacek sz�ma
    //***********************************************************************
    matrix<rvt> dc_ya_1, dc_ya_2;
    matrix<rvt> dc_xa_1, dc_xa_2;
    matrix<rvt> dc_xb_1, dc_xb_2;
    matrix<rvt> dc_yb_1, dc_yb_2;
    vektor<rvt> dc_ja_1, dc_ja_2;
    vektor<rvt> dc_jb_1, dc_jb_2;
    //***********************************************************************
    matrix<rvt> dc_nzbxa_1, dc_nzbxa_2;
    vektor<rvt> dc_nzbjb_1, dc_nzbjb_2;
    vektor<rvt> dc_ub_1, dc_ub_2;
    //***********************************************************************
public:
    //***********************************************************************
    void init_dc() override; // mem�riafoglal�s
    void resize_matrixok_dc() override; // kisz�m�tja a m�trixok m�ret�t, �s lefoglalja azokat, �j fa eset�n is
    void halozat_foglalas_cellaban_dc() override; // is_lin-t is ez �ll�tja
    void gerj_update_dc() override;
    void set_is_szimm() override; // A face-ek foglal�sa ut�n
    void update_dt_dc() override; // csak akkor kell k�l�n megh�vni, ha iter�ci� k�zben v�ltozik a dt
    void peremfeltetel_update_dc() override;
    void del_all_prev_dc() override;
    void update_for_uj_lepes_dc() override; // pl. kapacit�sok fesz�lts�ge�nek friss�t�se
    void update_iter_csp() override;
    // void update_for_uj_lepes_ac(); // AC-n�l nincs bels� iter�ci�
    //***********************************************************************
    // NR = Newton-Raphson
    // SA = Successive Approximation
    rvt  cella_klaszter_2_NR_UT_dc(rvt & Uc, rvt & Tc, rvt & T_hiba) override; // Ui = kiindulo_Ui + dUi, vissza: szumma dUi hiba �s a cella Uc �s Tc �rt�ke
    rvt  cella_klaszter_2_SA_UT_dc(rvt & T_hiba) override; // csak hib�t sz�mol, vissza: szumma Ui-el�z�_Ui
    void cella_klaszter_3_parameterek_frissitese_dc() override; // h�m�rs�kletf�gg� vezet�sek, junction vezet�se �gfesz�lts�gb�l
    void cella_klaszter_4_agaramok_szamitasa_dc() override; // disszip�ci� is
    rvt  cella_klaszter_5_NR_hibaaramok_szamitasa_dc() override; // J felt�lt�se, vissza: cellak�z�pponti I+P hiba�ram
    void cella_klaszter_5_SA_fill_J_dc() override; //
    void fw_klaszter_1_update_Jakobi_dc() override; // NR eset�ben a Jakobi, SA eset�ben az admittancia m�trixot t�lti fel a nemline�ris �rt�kekkel
    //***********************************************************************
    //                             *****  AC  *****
    //***********************************************************************
private:
    ::std::complex<rvt> sum_disszipacio_ac;
    vektor<face_tok_ac> faces_ac;   // face_tok_ac vektor, a 0 index� dummy
    face_tok_ac el_center_face_ac, th_center_face_ac;
    ref_vektor<face_tok_ac> perem_faces_ac_1, perem_faces_ac_2; // face_tok_ac vektor, a 0 index� dummy
    ref_vektor<face_tok_ac> kulso_faces_ac_1, kulso_faces_ac_2; // face_tok_ac vektor, a 0 index� dummy
public:
    void init_ac() override; // mem�riafoglal�s
    void face_foglalas_ac() override;
    void gerj_update_ac() override;
    void peremfeltetel_update_ac() override;
    void del_all_prev_ac() override;
    void update_ac(bool force) override;
    //***********************************************************************
public:
    //***********************************************************************
    faces_cella() :is_nemszimm{ false }, is_NR_nincs_csatolt_nemszimm{ false }, is_NR_van_csatolt_nemszimm{ false},
        kapcsolodo_cella{ nullptr },
        el_centroid_index_bemeneti{ 0 }, th_centroid_index_bemeneti{ 0 },
        sum_disszipacio_dc{ rvt() }, sum_disszipacio_ac{ rvt() }, is_perem_face { false },
        p_Yee_dc{ nullptr }, p_Yet_dc{ nullptr }, p_Yte_dc{ nullptr }, p_Ytt_dc{ nullptr }, p_Je_dc{ nullptr }, p_Jt_dc{ nullptr },
        el_centroid_index_matrix{ 0 }, th_centroid_index_matrix{ 0 }, konst_Yee_dc{ rvt() }, konst_Yet_dc{ rvt() }, konst_Yte_dc{ rvt() },
        konst_Ytt_dc{ rvt() }, konst_Je_dc{ rvt() }, konst_Jt_dc{ rvt() },
        feltoltes_tipus{ ft_semmi }, sum_ddissz_per_dT_dc{ rvt() }, p_sum_disszipacio_dc{ &sum_disszipacio_dc},
        p_sum_ddissz_per_dT_dc{ &sum_ddissz_per_dT_dc } {}
    //***********************************************************************
    void pre_init_egyedi() override;
    void fw_ac(uns melyik_fa) override;
    void bw_ac(uns melyik_fa) override;
    void fw_dc(uns melyik_fa) override;
    void bw_dc(uns melyik_fa) override;
    void set_to_NR_nincs_csatolt_nemszimm() { is_NR_nincs_csatolt_nemszimm = true; } // A face h�vja, ha nemszimmetrikus m�trixot hoz l�tre
    void set_to_NR_van_csatolt_nemszimm() { is_NR_van_csatolt_nemszimm = true; } // A face h�vja, ha nemszimmetrikus m�trixot hoz l�tre
    bool is_nonsymmetrical()const override { return is_nemszimm; }
    //***********************************************************************

    //***********************************************************************
    rvt get_Tc_akt() const {
    //***********************************************************************
        if (is_th)
            return th_center_face_dc.emlekek.get_akt().UT;
        if (kapcsolodo_cella != nullptr && kapcsolodo_cella->is_th)
            return kapcsolodo_cella->get_Tc_akt();
        return rvt();
    }
    //***********************************************************************
    rvt get_Tc_elozo() const {
    //***********************************************************************
        if (is_th)
            return th_center_face_dc.emlekek.get_elozo().UT;
        if (kapcsolodo_cella != nullptr && kapcsolodo_cella->is_th)
            return kapcsolodo_cella->get_Tc_elozo();
        return rvt();
    }
    //***********************************************************************
    rvt get_Tc_megtartando() const {
    //***********************************************************************
        if (is_th)
            return th_center_face_dc.emlekek.get_megtartando().UT;
        if (kapcsolodo_cella != nullptr && kapcsolodo_cella->is_th)
            return kapcsolodo_cella->get_Tc_megtartando();
        return rvt();
    }
    //***********************************************************************
    void print_G() const;
    //***********************************************************************
    void debug_write(::std::ofstream & fs) const override;
    //***********************************************************************
    rvt get_prev_IP_hiba() {
    //***********************************************************************
        rvt hibaI = rvt();
        uns db = 0;
        if (is_el) {
            hibaI = abs(*p_Je_dc);
            db = 1;
        }
        if (is_th) {
            hibaI += abs(*p_Jt_dc);
            db++;
        }
        return hibaI / db;
    }
    //***********************************************************************
    void randomize_zaj();
    //***********************************************************************
    cellatipus get_tipus()const override { return ct_faces_cella; }
    //***********************************************************************
};


//***********************************************************************
// jel�l�s fa_AB_bels�_CDE
// A: C/c/k/o/u/r/w: Centroid/csatlakoz�/k�zponti/open/UT/HTC/nem ambientre csatlakoz� HTC
// B: a/d: AC/DC
// bels�:
//      r: ellen�ll�s
//      i: �ramgener�tor
//      u: fesz�lts�ggener�tor
//      c: kapacit�s
//      j: junction
//      s: seebeck tag => ennek m�g nincs defini�lva a konstansa
//      x: szakad�s
// C: k/t/e/d: konstans/h�m�rs�kletf�gg�/�ramf�gg�/mindkett�t�l f�gg�
// D: e/t/d: elektromos/termikus/mindk�t t�rben alkalmazhat�
// E: d/n: disszip�l�/nem disszip�l�
//***********************************************************************
enum face_azonosito { 
    fa_none             = 0,
    fa_csatlakozo       = 0x00000001,
    fa_peremre          = 0x00000002,
    fa_kozepponti       = 0x00000004,
    fa_centroid         = 0x00000008,
    fa_perem_open       = 0x00000010,
    fa_perem_UT         = 0x00000020,
    fa_perem_HTC        = 0x00000040,
    fa_perem_HTC_T      = 0x00000060,
    fa_DC               = 0x00000000, 
    fa_AC               = 0x00000100, 
    fa_nem_hom_fuggo    = 0,
    fa_hom_fuggo        = 0x00000200,
    fa_nem_I_fuggo      = 0,
    fa_I_fuggo          = 0x00000400,
    fa_nem_disszipalo   = 0,
    fa_disszipalo       = 0x00000800,
    fa_nem_el_specific  = 0,
    fa_el_specific      = 0x00001000,
    fa_nem_th_specific  = 0,
    fa_th_specific      = 0x00002000,
    fa_belso_R          = 0x00010000,
    fa_belso_C          = 0x00020000,
    fa_belso_IP         = 0x00040000,
    fa_belso_UT         = 0x00080000,
    fa_belso_J          = 0x00100000, // A J(unction) az R|C|I|U-val sorba kapcsol�dik.

    fa_cd_r_kdn         = fa_csatlakozo | fa_DC | fa_belso_R, // dc_csatlakozo_face_ellenallas_elth_nem_disszipalo_const
    fa_od_r_kdn         = fa_peremre    | fa_DC | fa_perem_open | fa_belso_R, // dc_perem_open_face_ellenallas_elth
    fa_cd_r_ked         = fa_csatlakozo | fa_DC | fa_disszipalo | fa_el_specific | fa_belso_R, // dc_csatlakozo_face_ellenallas_el_disszipalo_const
    fa_kd_i_kdn         = fa_kozepponti | fa_DC | fa_belso_IP, // dc_kozepponti_face_IP_generator_konst
    fa_kd_ci_ktn        = fa_kozepponti | fa_DC | fa_th_specific | fa_belso_C | fa_belso_IP, // dc_kozepponti_face_Cth_IP_generator_konst
    fa_kd_u_kdn         = fa_kozepponti | fa_DC | fa_belso_UT, // dc_kozepponti_face_UT_generator_konst
    fa_cd_r_tdn         = fa_csatlakozo | fa_DC | fa_hom_fuggo | fa_belso_R, // dc_csatlakozo_face_ellenallas_elth_nem_disszipalo_Tfuggo
    fa_cd_r_ted         = fa_csatlakozo | fa_DC | fa_hom_fuggo | fa_disszipalo | fa_el_specific | fa_belso_R, // dc_csatlakozo_face_ellenallas_el_disszipalo_Tfuggo
    fa_cd_j_ted         = fa_csatlakozo | fa_DC | fa_hom_fuggo | fa_disszipalo | fa_el_specific | fa_belso_R | fa_belso_J, // dc_csatlakozo_face_junction_el_disszipalo_Tfuggo
    fa_kd_ci_ttn        = fa_kozepponti | fa_DC | fa_hom_fuggo | fa_th_specific | fa_belso_C | fa_belso_IP, // dc_kozepponti_face_Cth_IP_generator_Tfuggo
    fa_rd_r_kdn         = fa_peremre    | fa_DC | fa_perem_HTC | fa_belso_R, // dc_perem_HTC_face_ellenallas_elth_nem_disszipalo_konst
    fa_rd_r_tdn         = fa_peremre    | fa_DC | fa_perem_HTC | fa_hom_fuggo | fa_belso_R, // dc_perem_HTC_face_ellenallas_elth_nem_disszipalo_Tfuggo
    fa_ud_r_kdn         = fa_peremre    | fa_DC | fa_perem_UT  | fa_belso_R, // dc_perem_UT_face_ellenallas_elth_nem_disszipalo_konst
    fa_ud_r_tdn         = fa_peremre    | fa_DC | fa_perem_UT  | fa_hom_fuggo | fa_belso_R, // dc_perem_UT_face_ellenallas_elth_nem_disszipalo_Tfuggo
    fa_ud_r_ked         = fa_peremre    | fa_DC | fa_perem_UT  | fa_disszipalo | fa_el_specific | fa_belso_R, // dc_perem_U_face_ellenallas_el_disszipalo_konst
    fa_ud_r_ted         = fa_peremre    | fa_DC | fa_perem_UT  | fa_hom_fuggo | fa_disszipalo | fa_el_specific | fa_belso_R, // dc_perem_U_face_ellenallas_el_disszipalo_Tfuggo
    fa_kd_x_kdn         = fa_kozepponti | fa_DC, // dc_kozepponti_face_szakadas
    fa_Cd__d            = fa_centroid   | fa_DC
};
//***********************************************************************


//***********************************************************************
class os_face_core_ac {
//***********************************************************************
protected:
    //***********************************************************************
    bool is_el;
    face_tok_ac * p_tok;
    faces_cella * p_cella;
    const adat_face * p_adat_face;
    //***********************************************************************
public:
    //***********************************************************************
    os_face_core_ac(face_tok_ac * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :p_tok{ p_tok }, p_cella{ p_cella }, p_adat_face{ p_adat_face }, is_el{ is_el } {}
    //***********************************************************************
    virtual ~os_face_core_ac() {}
    //***********************************************************************
    virtual face_azonosito get_azon() const = 0;
    //***********************************************************************
    // mindig konstans
    //***********************************************************************
};


//***********************************************************************
class os_face_core_dc {
//***********************************************************************
    // ugyanazt a face-t haszn�ljuk dc �s tranziens esetben is
protected:
    //***********************************************************************
    bool is_el, is_kell_rand;
    uns csatlakozo_index;
    face_tok_dc * p_tok;
    faces_cella * p_cella;
    const adat_face * p_adat_face;
    const adat_anyag * p_anyag;
    const adat_junction * p_junction;
    double zaj_szorzo; // default: 1,  a float itt nem el�g pontos
    //***********************************************************************
    void randomize_zaj_belso() {
    //***********************************************************************
        if (!is_kell_rand)
            return;
        uns zaj_index = p_adat_face == nullptr ? p_cella->p_adatcella->zaj_index : p_adat_face->zaj_index;
        if (zaj_index == 0) { // a k�zponti zaj �rt�k sz�m�t
            if (akt_sim.p_akt_sim->zaj == rvt()) { // nincs zaj
                zaj_szorzo = 1.0;
            }
            else {
                zaj_szorzo = 1.0 + akt_sim.p_akt_sim->zaj*((2.0*rand() / RAND_MAX) - 1.0);
            }
        }
        else {
            zaj_szorzo = 1.0 + akt_sim.p_akt_sim->meretek[p_adat_face->zaj_index];
        }
    }
public:
    //***********************************************************************
    os_face_core_dc(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :p_tok{ p_tok }, p_cella{ p_cella }, p_adat_face{ p_adat_face }, p_anyag{ nullptr }, p_junction{ nullptr },
        is_el{ is_el }, is_kell_rand{ true }, csatlakozo_index{ p_adat_face==nullptr ? 0 : p_adat_face->csatlakozo_index },
        zaj_szorzo{ 1 } {}
    //***********************************************************************
    bool is_junction()const { return p_junction != nullptr; }
    //***********************************************************************
    uns get_csatlakozo_index()const { return csatlakozo_index; }
    //***********************************************************************
    void randomize_zaj() {
    //***********************************************************************
        if (!is_kell_rand)
            return;
        randomize_zaj_belso();
        is_kell_rand = false;
    }
    //***********************************************************************
    virtual ~os_face_core_dc() {}
    //***********************************************************************
    virtual face_azonosito get_azon() const = 0;
    //***********************************************************************
    virtual bool is_konst() const = 0;
    //***********************************************************************
    virtual void init() = 0;
    //***********************************************************************
    virtual void set_cella_is_szimm() {}; // amelyik face nemszimmetrikuss� teszi az admittanciam�trixot, az be�ll�tja ezt
    //***********************************************************************
    virtual void update_uj_lepeshez() = 0; // kapacit�sok friss�t�se, konstansok admittanciam�trixokba �r�sa, konstans k�z�pponti admittanci�k n�vel�se
    //***********************************************************************
    virtual void update_admittanciamatrix_SA() = 0; // a nem konstans elemekkel friss�ti az admittanciam�trixot
    //***********************************************************************
    virtual void update_Jakobi_NR_nincs_csatolt() = 0; // a nem konstans elemekkel friss�ti a Jakobi m�trixot (kisz�m�tja a deriv�ltakat is)
    //***********************************************************************
    virtual void update_Jakobi_NR_van_csatolt() = 0; // a nem konstans elemekkel friss�ti a Jakobi m�trixot (kisz�m�tja a deriv�ltakat is)
    //***********************************************************************
    virtual void update_nemlinearis_parameterek() = 0; // fesz�lts�gek/h�m�rs�kletek alapj�n a face be�ll�t�sa
    //***********************************************************************
    virtual void update_UT_IP_NR() = 0; // csatlakoz�, perem �s k�z�pponti face-eken elt�r�, de csak a peremn�l egyedi, a m�sik kett�n�l egys�ges
    //***********************************************************************
    virtual void update_UT_IP_SA() = 0; // csatlakoz�, perem �s k�z�pponti face-eken elt�r�, de csak a peremn�l egyedi, a m�sik kett�n�l egys�ges
    //***********************************************************************
    virtual void agaram_szamitasa() = 0; // disszip�ci� �s peremekn�l a k�ls� UT is
    //***********************************************************************
    void update_inhom_SA_kereszt() {} // TODO: ha lesz olyan keresztface, ami m�dos�tja az inhomog�n �ramot, akkor virtu�liss� kell tenni
    //***********************************************************************
    virtual void update_inhom_SA_kozepponti(rvt dissz) {}
    //***********************************************************************
    virtual void update_hibaaram_NR_kereszt(rvt & sum_el_hiba, rvt & sum_th_hiba) = 0;
    //***********************************************************************
    virtual void update_hibaaram_NR_kozepponti(rvt & sum_el_hiba, rvt & sum_th_hiba) {}
    //***********************************************************************
    virtual void update_dt() {} // ha kapacit�s van, �s v�ltozik a dt tranziens szimul�ci�n�l
    //***********************************************************************
    virtual void update_gerj() = 0; // egyel�re csak a k�z�pponti face-ekn�l
    //***********************************************************************
    virtual void update_peremfeltetel() = 0; // peremface-ekn�l
    //***********************************************************************
    static os_face_core_dc * alloc_new(face_azonosito tipus, face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el);
    //***********************************************************************
    bool is_El()const { return is_el; }
    //***********************************************************************
    virtual void print_G()const {}
    //***********************************************************************
    virtual void debug_write(::std::ofstream & fs, bool is_endl) const { 
    //***********************************************************************
        fs << (is_el ? "el " : "th ");
        p_tok->debug_write(fs);
        if (is_endl)
            fs << ::std::endl;
    }
    //***********************************************************************
};


//***********************************************************************
class dc_centroid : public os_face_core_dc {
//***********************************************************************
public:
    //***********************************************************************
    dc_centroid(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :os_face_core_dc(p_tok, p_cella, p_adat_face, is_el) {}
    //***********************************************************************

    //***********************************************************************
    face_azonosito get_azon() const override { return fa_Cd__d; }
    //***********************************************************************
    bool is_konst() const override { return true; }
    //***********************************************************************
    void init() override {}
    //***********************************************************************
    void update_uj_lepeshez() override {}
    //***********************************************************************
    void update_admittanciamatrix_SA() override {}
    //***********************************************************************
    void update_Jakobi_NR_nincs_csatolt() override {}
    //***********************************************************************
    void update_Jakobi_NR_van_csatolt() override {}
    //***********************************************************************
    void update_nemlinearis_parameterek() override {}
    //***********************************************************************
    void update_UT_IP_NR() override {}
    //***********************************************************************
    void update_UT_IP_SA() override {}
    //***********************************************************************
    void agaram_szamitasa() override {}
    //***********************************************************************
    void update_gerj() override {}
    //***********************************************************************
    void update_peremfeltetel() override {}
    //***********************************************************************
    void update_hibaaram_NR_kereszt(rvt & sum_el_hiba, rvt & sum_th_hiba) override { TODO("Centroid �ram�t meg kell hat�rozni."); }
    //***********************************************************************
    virtual void debug_write(::std::ofstream & fs, bool is_endl) const override { fs << "centroid "; os_face_core_dc::debug_write(fs, false); if (is_endl) fs << ::std::endl; }
    //***********************************************************************
};


//***********************************************************************
class dc_keresztface : public os_face_core_dc {
//***********************************************************************
protected:
    rvt A_per_L, A;
    // emlekezo::tar<rvt> G; ne legyen az �sben, mert konstansokn�l nem friss�tj�k mindig
public:
    //***********************************************************************
    dc_keresztface(face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el)
        :os_face_core_dc(p_tok, p_cella, p_adat_face, is_el), A_per_L{ rvt() }, A{ rvt() } { /*G = rvt(g_min);*/ }
    //***********************************************************************

    //***********************************************************************
    void init() {
    //***********************************************************************
        A = akt_sim.p_akt_sim->meretek[p_adat_face->A_index];
        p_cella->A_cella += A;
        if (p_adat_face->oldal == pit_top)
            p_cella->A_top = A;
        A_per_L = rvt(A / akt_sim.p_akt_sim->meretek[p_adat_face->L_index] * zaj_szorzo);
        p_anyag = &akt_sim.p_akt_sim->anyagok[akt_sim.akt_anyag_index[p_adat_face->anyag_index]];
        p_junction = p_adat_face->junction_index==0 ? nullptr : &akt_sim.p_akt_sim->junctionok[akt_sim.akt_junct_index[p_adat_face->junction_index]];
    }
    //***********************************************************************
    void update_gerj() override {}
    //***********************************************************************
    //void print_G()const override { printf("%g\n", G); }
    //***********************************************************************
    virtual void debug_write(::std::ofstream & fs, bool is_endl) const override { if (is_endl) fs << "kereszt  "; os_face_core_dc::debug_write(fs, false); fs << " G=TODO" /*<< ::std::setw(12) << G.get_akt()*/; if (is_endl) fs << ::std::endl; }
    //***********************************************************************
};


}


//***********************************************************************
#include "face_csatlakozo.h"
#include "face_perem.h"
#include "face_kozepponti.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
inline os_face_core_dc * os_face_core_dc::alloc_new(face_azonosito tipus, face_tok_dc * p_tok, faces_cella * p_cella, const adat_face * p_adat_face, bool is_el) {
//***********************************************************************
    os_face_core_dc * pret = nullptr;
    switch (tipus) {
        case fa_Cd__d:      pret = new dc_centroid(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_cd_r_kdn:   pret = new dc_csatlakozo_face_ellenallas_elth_nem_disszipalo_const(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_od_r_kdn:   pret = new dc_perem_open_face_ellenallas_elth(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_cd_r_ked:   pret = new dc_csatlakozo_face_ellenallas_el_disszipalo_const(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_kd_i_kdn:   pret = new dc_kozepponti_face_IP_generator_konst(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_kd_ci_ktn:  pret = new dc_kozepponti_face_Cth_IP_generator_konst(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_kd_u_kdn:   pret = new dc_kozepponti_face_UT_generator_konst(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_cd_r_tdn:   pret = new dc_csatlakozo_face_ellenallas_elth_nem_disszipalo_Tfuggo(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_cd_r_ted:   pret = new dc_csatlakozo_face_ellenallas_el_disszipalo_Tfuggo(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_cd_j_ted:   pret = new dc_csatlakozo_face_junction_el_disszipalo_Tfuggo(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_kd_ci_ttn:  pret = new dc_kozepponti_face_Cth_IP_generator_Tfuggo(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_rd_r_kdn:   pret = new dc_perem_HTC_face_ellenallas_elth_nem_disszipalo_konst(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_rd_r_tdn:   pret = new dc_perem_HTC_face_ellenallas_elth_nem_disszipalo_Tfuggo(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_ud_r_kdn:   pret = new dc_perem_UT_face_ellenallas_elth_nem_disszipalo_konst(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_ud_r_tdn:   pret = new dc_perem_UT_face_ellenallas_elth_nem_disszipalo_Tfuggo(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_ud_r_ked:   pret = new dc_perem_U_face_ellenallas_el_disszipalo_konst(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_ud_r_ted:   pret = new dc_perem_U_face_ellenallas_el_disszipalo_Tfuggo(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_kd_x_kdn:   pret = new dc_kozepponti_face_szakadas(p_tok, p_cella, p_adat_face, is_el); break;
        case fa_none: throw hiba("os_face_core_dc::alloc_new", "none type");
        default: throw hiba("os_face_core_dc::alloc_new", "unknown type");
    }
    if (pret->get_azon() != tipus)
        throw hiba("os_face_core_dc::alloc_new", "byd type binding");
    return pret;
}


//***********************************************************************
inline face_tok_dc::~face_tok_dc() { delete p_core; }
inline face_tok_ac::~face_tok_ac() { delete p_core; }
//***********************************************************************


}

#endif
