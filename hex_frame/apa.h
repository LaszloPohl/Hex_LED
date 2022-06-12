//***********************************************************************
// apa class header
// Creation date:  2009. 07. 11.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef VSUN3_APA_HEADER
#define	VSUN3_APA_HEADER
//***********************************************************************


//***********************************************************************
#include "tipusok.h"
#include "gfunc.h"
#include "sugarak.h"
//***********************************************************************


//***********************************************************************
class meret_tomb_tipus{
//***********************************************************************
    std::vector<dbl> t;
    struct meret_index_par { dbl meret; size_t index; };
    meret_index_par last_ten[10];
    size_t lt_index;
public:
    //***********************************************************************
    meret_tomb_tipus() {
    //***********************************************************************
        clear();
    }

    //***********************************************************************
    size_t get_index(dbl meret) {
    //***********************************************************************
        for (size_t i = 0; i < 10; i++)
            if (last_ten[i].meret == meret)
                return last_ten[i].index;
        for (size_t i = 0; i < t.size(); i++)
            if (t[i] == meret) {
                last_ten[lt_index].meret = meret;
                last_ten[lt_index].index = i;
                lt_index = (lt_index + 1) % 10;
                return i;
            }
        t.push_back(meret);
        last_ten[lt_index].meret = meret;
        last_ten[lt_index].index = t.size() - 1;
        lt_index = (lt_index + 1) % 10;
        return t.size() - 1;
    }

    //***********************************************************************
    void clear() {
    //***********************************************************************
        lt_index = 0;
        t.clear();
        t.push_back(-1.0);
        for (size_t i = 0; i < 10; i++)
            last_ten[i] = meret_index_par{ -1.0, 0 };
    }

    //***********************************************************************
    const std::vector<dbl> & get_vector()const { return t; }
    //***********************************************************************
};


//***********************************************************************
struct real_cell_res{ // potenciálok ill. áramok visszaadására szolgál
//***********************************************************************
    dbl t[BASIC_SIDES]; // EXTERNAL=CENTER módon használjuk helytakarékosság miatt
    real_cell_res(){for(uns i=0;i<BASIC_SIDES;i++)t[i]=nulla;}
};


//***********************************************************************
struct GTpar{
//***********************************************************************
    dbl G[3],T; // dbl G,T;
    bool is_F;
    bool is_isotrop()const { return G[0] == G[1] && G[1] == G[2]; }
    GTpar() :is_F{ false } {}
};


//***********************************************************************
struct material;
struct v6eredm;
//***********************************************************************


//***********************************************************************
struct lsfit_adat {
//***********************************************************************
    dbl U, T, I; // A T mindig hõmérséklet, a másik kettõ bármi lehet. U-ra (és T-re, ha trió) illeszt, az eredmény I
};


//***********************************************************************
struct lsfit {
//***********************************************************************
    bool is_trio; // duó vagy trió? azaz egy vagy két adatra kell fittelni a harmadikat
    dbl multiplier;
    lsfit_veg_tipus start, stop;
    lsfit_egyenlet_tipus egyenlet;
    uns egyenlet_unspar_1, egyenlet_unspar_2; // polinomnál az U és T polinom foka
    tomb<lsfit_adat> meresek;
    bool read(const PLString & fajlnev);
    lsfit() :is_trio{ true }, multiplier{ 1 }, start{ lsv_none }, stop{ lsv_none }, egyenlet{ lse_polinom },
        egyenlet_unspar_1{ 0 }, egyenlet_unspar_2{ 0 } {}
};


//***********************************************************************
#define GGSIZE 8
//***********************************************************************
struct vezetes{
//***********************************************************************
    monlintipus tipus,semitip;  // enum monlintipus{nlt_lin,nlt_exp,nlt_diode,nlt_quadratic,nlt_erno,nlt_szakaszok,nlt_mizs1};
    dbl g[3],gg[GGSIZE],T;           // dbl g[3],gg[GGSIZE],T;
    tomb<GTpar> szakaszok;      // tomb<GZpar> szakaszok;
    lsfit ls_fit_adatok;        // ha legkisebb négyzet illesztõ
    bool specific;              // bool, Ernõ típusnál ha a mátrix 1 m2-re vonatkozik, azaz a félvezetõ felületével szorzandó
    bool is_resistivity;        // ha true, a felhasználó fajlagos ellenállást adott meg vezetés helyett, a hõmérsékletfüggés függgvényét erre kell alkalmazni
    bool is_his;
    dbl his_value_1,his_value_2; // A hiszterézises tartományban az érték. Ha value_2=0, akkor kapacitás, azaz a visszaadott érték value_1/(T_max-T_min), egyébként vezetés, ami lineárisan skálázódik min és max között
    // nlt_erno esetén g lesz B×d, gg pedig M×d
    vezetes(dbl default_g=nulla):tipus(nlt_lin),semitip(nlt_lin),T(nulla),specific(false),is_resistivity(false)
        , his_value_2(nulla), is_his(false) {
        g[0] = g[1] = g[2] = default_g; 
        for (uns i = 0; i < GGSIZE; i++) 
            gg[i] = nulla;
    }
    dbl get_ertek(dbl T)const;
    bool is_szakaszok_isotrop()const {
        if (tipus != nlt_szakaszok)return true;
        for (uns i = 0; i < szakaszok.size(); i++)
            if (!szakaszok[i].is_isotrop())
                return false;
        return true;
    }
    bool is_isotrop()const { return g[0] == g[1] && g[1] == g[2] && is_szakaszok_isotrop(); }
    void write_const(FILE *fp)const {
        if (tipus != nlt_lin)
            throw hiba("vezetes::write_const", "only const value is acceptable");
        fprintf(fp, "%g;\n",g[0]);
    }
    void write_normal(FILE *fp, uns irany, bool is_eol = true) const{ // irány: 0=x, 1=y, 2=z
        if(is_resistivity)
            fprintf(fp, "I");
        else
            fprintf(fp, "N");
        //if (is_his)
        //    fprintf(fp, "H%g;%g;", his_value_1, his_value_2);
        switch (tipus) {
            case nlt_lin:
                fprintf(fp, "C%g;", g[irany]);
                if (is_eol)fprintf(fp, "\n");
                break;
            case nlt_linearis:
                fprintf(fp, "L%g;%g;", gg[0], gg[1]);
                if (is_eol)fprintf(fp, "\n");
                break;
            case nlt_exp:
                fprintf(fp, "X%g;%g;", g[irany], gg[0]);
                if (is_eol)fprintf(fp, "\n");
                break;
            case nlt_quadratic:
                hiba("vezetes::write_normal", "quadratic property is not supported");
                break;
            case nlt_szakaszok: {
                bool is_fazis_broken = false;
                for (uns i = 0; i < szakaszok.size(); i++)
                    if (szakaszok[i].is_F)
                        is_fazis_broken = true;
                fprintf(fp, "%c%u", is_fazis_broken ? 'F' : 'B', szakaszok.size());
                for (uns i = 0; i < szakaszok.size(); i++) {
                    if (szakaszok[i].is_F)
                        fprintf(fp, "F%g;", szakaszok[i].G[irany]);
                    else
                        fprintf(fp, "T%g;%g;", szakaszok[i].T, szakaszok[i].G[irany]);
                }
                if (is_eol)fprintf(fp, "\n");
                break;
            }
            case nlt_mizs1:
                fprintf(fp, "M");
                for (uns i = 0; i < 7; i++)
                    fprintf(fp, "%g;", gg[i]);
                if (is_eol)fprintf(fp, "\n");
                break;
            default:
                throw hiba("vezetes::write_normal", "unsupported nonlin type");
                break;
        }
    }
    void write_semi(FILE *fp) { 
        switch (semitip) {
        case nlt_exp:
            fprintf(fp, "NX%g;1.0;\n", g[0]);
            break;
        case nlt_diode:
            fprintf(fp, "ND%g;%g;%g;\n", g[0], g[1], g[2]);
            break;
        case nlt_quadratic:
            hiba("vezetes::write_semi", "quadratic property is not supported");
            break;
        case nlt_erno:
            fprintf(fp, "NE%g;%g;%g;%g;%g;%g;\n", g[0], g[1], g[2], gg[0], gg[1], gg[2]);
            break;
        case nlt_lsfit: {
                if (ls_fit_adatok.egyenlet == lse_polinom) {
                    char eleje = 'N', vege = 'N';
                    switch (ls_fit_adatok.start) {
                        case lsv_lin:    eleje = 'L'; break;
                        case lsv_strong: eleje = 'S'; break;
                    }
                    switch (ls_fit_adatok.stop) {
                        case lsv_lin:    vege = 'L'; break;
                        case lsv_strong: vege = 'S'; break;
                    }
                    if (ls_fit_adatok.is_trio) {
                        fprintf(fp, "NI%g;TPV%uT%u%c%c%uM", ls_fit_adatok.multiplier, ls_fit_adatok.egyenlet_unspar_1, ls_fit_adatok.egyenlet_unspar_2, eleje, vege, ls_fit_adatok.meresek.size());
                        for (uns i = 0; i < ls_fit_adatok.meresek.size(); i++)
                            fprintf(fp, "%g;%g;%g;", ls_fit_adatok.meresek[i].U, ls_fit_adatok.meresek[i].T, ls_fit_adatok.meresek[i].I);
                        fprintf(fp, "\n");
                    }
                    else { // duó
                        fprintf(fp, "NI%g;DPV%u%c%c%uM", ls_fit_adatok.multiplier, ls_fit_adatok.egyenlet_unspar_1, eleje, vege, ls_fit_adatok.meresek.size());
                        for (uns i = 0; i < ls_fit_adatok.meresek.size(); i++)
                            fprintf(fp, "%g;%g;", ls_fit_adatok.meresek[i].U, ls_fit_adatok.meresek[i].I);
                        fprintf(fp, "\n");
                    }
                }
                else
                    throw hiba("vezetes::write_semi", "unsupported fit equation type");
            }
            break;
        default:
            throw hiba("vezetes::write_semi", "unsupported nonlin type");
            break;
        }
    }
    void write_sem_rad_lum(FILE *fp) {
        switch (semitip) {
        case nlt_lin:
            fprintf(fp, "NC%g;\n", g[0]);
            break;
        case nlt_linearis:
            fprintf(fp, "NL%g;%g;\n", gg[0], gg[1]);
            break;
        case nlt_erno:
            fprintf(fp, "NE%g;%g;%g;%g;%g;%g;\n", g[0], g[1], g[2], gg[0], gg[1], gg[2]);
            break;
        case nlt_lsfit: {
                if (ls_fit_adatok.egyenlet == lse_polinom) {
                    char eleje = 'N', vege = 'N';
                    switch (ls_fit_adatok.start) {
                        case lsv_lin:    eleje = 'L'; break;
                        case lsv_strong: eleje = 'S'; break;
                    }
                    switch (ls_fit_adatok.stop) {
                        case lsv_lin:    vege = 'L'; break;
                        case lsv_strong: vege = 'S'; break;
                    }
                    if (ls_fit_adatok.is_trio) {
                        fprintf(fp, "NI%g;TPV%uT%u%c%c%uM", ls_fit_adatok.multiplier, ls_fit_adatok.egyenlet_unspar_1, ls_fit_adatok.egyenlet_unspar_2, eleje, vege, ls_fit_adatok.meresek.size());
                        for (uns i = 0; i < ls_fit_adatok.meresek.size(); i++)
                            fprintf(fp, "%g;%g;%g;", ls_fit_adatok.meresek[i].U, ls_fit_adatok.meresek[i].T, ls_fit_adatok.meresek[i].I);
                        fprintf(fp, "\n");
                    }
                    else { // duó
                        fprintf(fp, "NI%g;DPV%u%c%c%uM", ls_fit_adatok.multiplier, ls_fit_adatok.egyenlet_unspar_1, eleje, vege, ls_fit_adatok.meresek.size());
                        for (uns i = 0; i < ls_fit_adatok.meresek.size(); i++)
                            fprintf(fp, "%g;%g;", ls_fit_adatok.meresek[i].U, ls_fit_adatok.meresek[i].I);
                        fprintf(fp, "\n");
                    }
                }
                else
                    throw hiba("vezetes::write_sem_rad_lum", "unsupported fit equation type");
            }
            break;
        default:
            throw hiba("vezetes::write_sem_rad_lum", "unsupported nonlin type, only const and erno is supported");
            break;
        }
    }
};

//***********************************************************************
inline dbl vezetes::get_ertek(dbl T)const { // irany: 0=x, 1=y, 2=z
// Hõmérsékletfüggõ vezetést ad vissza
//***********************************************************************
    if (is_resistivity) { // ha a g mátrixban nem vezetés, hanem ellenállás van
        switch (tipus) {
            case nlt_lin: return egy / g[0];
            case nlt_linearis: return egy / (gg[0] + T * gg[1]);
            case nlt_exp: {
                cd szorzat = gg[0] * (T - 25.0);
                cd kitevo = szorzat > 7.0 ? 7.0 : szorzat < -7.0 ? -7.0 : szorzat;
                return egy / (g[0] * exp(kitevo));
            }
            case nlt_diode: throw hiba("vezetes::get_ertek", "diode equation is not applicable");
            case nlt_quadratic: return egy / (g[0] * (gg[0] * T * T + gg[1] * T + gg[2]));
            case nlt_szakaszok: {
                cu32 n = szakaszok.size();
                dbl G0, G1, T0, T1;
                if (n == 0)throw hiba("vezetes::get_ertek", "0 vertice polygon");
                if (n == 1)return egy / szakaszok[0].G[0];
                if (T <= szakaszok[0].T) {
                    return egy / szakaszok[0].G[0];
                }
                else if (T >= szakaszok[n - 1].T) {
                    return egy / szakaszok[n - 1].G[0];
                }
                else {
                    u32 i = 1;
                    while (T > szakaszok[i].T) i++;
                    G0 = szakaszok[i - 1].G[0];  G1 = szakaszok[i].G[0];
                    T0 = szakaszok[i - 1].T;     T1 = szakaszok[i].T;
                    return egy / (G0 + (G1 - G0) * (T - T0) / (T1 - T0));
                }
            }
            case nlt_mizs1:
                return egy / fn_mizs1(T, gg[0], gg[1], gg[2], gg[3], gg[4], gg[5], gg[6]);
            }
    }
    else {
        switch (tipus) {
            case nlt_lin: return g[0];
            case nlt_linearis: return gg[0] + T * gg[1];
            case nlt_exp: {
                cd szorzat = gg[0] * (T - 25.0);
                cd kitevo = szorzat>7.0 ? 7.0 : szorzat<-7.0 ? -7.0 : szorzat;
                return g[0] * exp(kitevo);
            }
            case nlt_diode: throw hiba("vezetes::get_ertek", "diode equation is not applicable");
            case nlt_quadratic: return g[0] * (gg[0] * T*T + gg[1] * T + gg[2]);
            case nlt_szakaszok: {
                cu32 n = szakaszok.size();
                dbl G0, G1, T0, T1, dC = nulla;
                if (n == 0)throw hiba("vezetes::get_ertek", "0 vertice polygon");
                if (n == 1)return szakaszok[0].G[0];
                if (T <= szakaszok[0].T) {
                    return dC + szakaszok[0].G[0];
                }
                else if (T >= szakaszok[n - 1].T) {
                    return dC + szakaszok[n - 1].G[0];
                }
                else {
                    u32 i = 1;
                    while (T > szakaszok[i].T) i++;
                    G0 = szakaszok[i - 1].G[0];  G1 = szakaszok[i].G[0];
                    T0 = szakaszok[i - 1].T;     T1 = szakaszok[i].T;
                    return dC + G0 + (G1 - G0) * (T - T0) / (T1 - T0);
                }
            }
            case nlt_mizs1:
                return fn_mizs1(T, gg[0], gg[1], gg[2], gg[3], gg[4], gg[5], gg[6]);
            }
    }
    throw hiba("vezetes::get_ertek", "unexpected control path (%u)", tipus);
}



//***********************************************************************
class vezetes_tomb_tipus{
//***********************************************************************
    std::vector<vezetes> t;
    struct vezetes_index_par { vezetes ertek; size_t index; };
    vezetes_index_par last_ten[10];
    size_t lt_index;
    //***********************************************************************
    bool is_egyforma(const vezetes & a, const vezetes & b) {
    //***********************************************************************
        if (a.tipus != b.tipus || a.is_resistivity != b.is_resistivity)
            return false;
        switch (a.tipus) {
            case nlt_lin:       return a.g[0] == b.g[0];
            case nlt_linearis:  return a.gg[0] == b.gg[0] && a.gg[1] == b.gg[1];
            case nlt_exp:       return a.gg[0] == b.gg[0] && a.g[0] == b.g[0] && a.g[1] == b.g[1] && a.g[2] == b.g[2];
            case nlt_szakaszok: {
                if (a.szakaszok.size() != b.szakaszok.size())
                    return false;
                for (uns i = 0; i < a.szakaszok.size(); i++) {
                    if (a.szakaszok[i].is_F != b.szakaszok[i].is_F)
                        return false;
                    if (a.szakaszok[i].G[0] != b.szakaszok[i].G[0] || a.szakaszok[i].G[1] != b.szakaszok[i].G[1] || a.szakaszok[i].G[2] != b.szakaszok[i].G[2])
                        return false;
                    if (!a.szakaszok[i].is_F && a.szakaszok[i].T != b.szakaszok[i].T)
                            return false;
                }
                return true;
            }
            default:
                throw hiba("vezetes_tomb_tipus::is_egyforma", "unsupported nonlin type");
                break;
        }
    }
public:
    //***********************************************************************
    vezetes_tomb_tipus() {
    //***********************************************************************
        clear();
    }

    //***********************************************************************
    uns get_index(const vezetes & uj_ertek) {
    //***********************************************************************
        for (size_t i = 0; i < 10; i++)
            if (is_egyforma(last_ten[i].ertek, uj_ertek))
                return (uns)last_ten[i].index;
        for (size_t i = 0; i < t.size(); i++)
            if (is_egyforma(t[i], uj_ertek)) {
                last_ten[lt_index].ertek = uj_ertek;
                last_ten[lt_index].index = i;
                lt_index = (lt_index + 1) % 10;
                return (uns)i;
            }
        t.push_back(uj_ertek);
        last_ten[lt_index].ertek = uj_ertek;
        last_ten[lt_index].index = t.size() - 1;
        lt_index = (lt_index + 1) % 10;
        return (uns)(t.size() - 1);
    }

    //***********************************************************************
    void clear() {
    //***********************************************************************
        lt_index = 0;
        t.clear();
        t.push_back(-1.0);
        for (size_t i = 0; i < 10; i++)
            last_ten[i] = vezetes_index_par{ -1.0, 0 };
    }

    //***********************************************************************
    const std::vector<vezetes> & get_vector()const { return t; }
    //***********************************************************************
};


//***********************************************************************
struct material{
//***********************************************************************
    PLString nev;               // PLString 
    bool is_el,is_th,is_his;    // bool, létezik-e a tér, van-e hiszterézis
    bool is_fenypor;
    bool is_phase_change_energy;// bool, volt-e definiálva
    char output_side;           // char, csak akkor menti az össz fényhez, ha ezt az oldalt metszi a sugár
    vezetes elvez;              // vezetes elvez;
    vezetes thvez;              // vezetes
    vezetes Ce,Cth,S;           // vezetesd
    vezetes D;                  // vezetes, 0..1, az esõ elektromos telj mekkora része fût
    vezetes emissivity;         // vezetés, 0..1, def 1.
    vezetes light_blue_absorption_coeff;    // vezetés, def 1
    vezetes light_conversion_efficiency;    // vezetés, 0..1, def 1
    vezetes light_yellow_absorption_coeff;  // vezetés, def 1
    vezetes light_re_conversion_efficiency; // vezetés, def 1
    vezetes light_vertical_light_conversion; // vezetés, 0..1, def 0, ideiglenes
    dbl reflectivity, yellow_correction; // dbl, 0..1, def 0
    uns direction_mode; // uns, 1, 2 3 vagy 4, def 1
    uns ray_per_cell_dir; // uns, csak páratlan lehet
    bool is_sarga_szetmegy; // bool, ha a sárga elválhat a kéktõl, és mindenfelé mehet
    dbl cut_level; // dbl, default 0.01, sok sugárnál eldobja az ez alattiakat
    dbl d_light_powder;
    dbl his_T_min, his_T_max, his_T_width_fele; // A hiszterézis sávjának alsó és felsõ középhõmérséklete, valamit a sáv szélessége (a szélesség fele-fele van a közép két oldalán)
    dbl his_T1,his_T2,his_T3,his_T4,his_delta,his_inv_delta;
    dbl phase_change_energy;
    uns anyagindex[3];

    material() :is_el(false), is_th(false), is_his(false), is_fenypor{ false }, is_phase_change_energy{ false }, D(egy),
        emissivity(egy), his_T_min(nulla), his_T_max(egy), his_T_width_fele(nulla), phase_change_energy{ 0.0 },
        light_blue_absorption_coeff{ egy }, light_conversion_efficiency{ egy }, light_yellow_absorption_coeff{ egy },
        light_vertical_light_conversion{ nulla }, direction_mode{ 1 }, is_sarga_szetmegy{ false }, cut_level{ 0.01 }, d_light_powder{ 0 },
        ray_per_cell_dir{ 1 }, output_side{ 'T' }, light_re_conversion_efficiency{ egy }, reflectivity{ 0 }, yellow_correction{ 0 } {}
    bool is_lin_el()const {
        if (elvez.tipus != nlt_lin)
            return false;
        if (S.tipus != nlt_lin)
            return false;
        return true;
    }
    bool is_lin_th()const {
        if (thvez.tipus != nlt_lin)
            return false;
        if (Cth.tipus != nlt_lin)
            return false;
        return true;
    }
    bool is_isotrop()const{
        if (is_el) {
            if (!elvez.is_isotrop())
                return false;
            if (!S.is_isotrop())
                return false;
            if (!D.is_isotrop())
                return false;
            if (!emissivity.is_isotrop())
                return false;
        }
        if (is_th) {
            if (!thvez.is_isotrop())
                return false;
            if (!Cth.is_isotrop())
                return false;
        }
        return true;
    }
    void write(FILE * fp, uns irany, uns index) { // irány: 0=x, 1=y, 2=z
        if (irany == 0)
            anyagindex[2] = anyagindex[1] = anyagindex[0] = index;
        else
            anyagindex[irany] = index;
        fprintf(fp, "BM%u\n", index);
        if (is_his) {
            if (!is_phase_change_energy)
                throw hiba("material::write", "phase_change_energy is missing");
            fprintf(fp, "F%g;H%g;%g;%g;\n", phase_change_energy, his_T_min, his_T_max, his_T_width_fele);
        }
        fprintf(fp, "PGT");
        thvez.write_normal(fp, irany);
        fprintf(fp, "PGE");
        elvez.write_normal(fp, irany);
        fprintf(fp, "PCT");
        Cth.write_normal(fp, irany);
        fprintf(fp, "PS");
        S.write_normal(fp, irany);
        fprintf(fp, "PD");
        D.write_normal(fp, irany);
        fprintf(fp, "PE");
        emissivity.write_normal(fp, irany);
        if (is_fenypor)
            fprintf(fp, "LP\n");
        fprintf(fp, "EM%u\n", index);
    }
};


//***********************************************************************
struct semiconductor{
//***********************************************************************
    static uns db;              // ststic uns, a létezõ félvezetõk darabszáma, az azonosító meghatározásához kell.
    uns azon;                   // uns, a félvezetõ azonosítója, a full_I_semi tömb indexeléséhez kell
    uns col2;                   // uns, a másik szín
    vezetes par;                // vezetes par, par.g[0]=Vt, par.g[1]=I0, vigyázat I0 1 négyzetméterre!!!, 3 par esetén másodfokú egyenlet, gg a hõm függés
    vezetes D;                  // vezetes, 0..1, az esõ elektromos telj mekkora része fût
    vezetes R;                  // vezetes, >=1, a sugárzott teljesítmény hányszorosa vonódjon le a számított disszipációból Dissz=U^2*G*D-F*R, F a radianciatomb megfelelõ eleme
    vezetes rad,lum;            // vezetes, az adott átmenethez tartozó radiancia ill luminancia értékeket tartalmazza
    dbl As;                     // dbl, mekkora felületen érintkezik a két anyag egymással
    uns index;
    semiconductor():col2(0),As(nulla),D(egy),R(egy){azon=db++;}
};


//***********************************************************************
struct color{
//***********************************************************************
    PLString nev;               // PLString
    bool is;                    // bool, ha van definiálva, true
    SzinTipus tipus;            // SzinTipus: SzinNormal,SzinBoundary,SzinUchannel
    dbl terfogat;               // dbl
    material * pmat;            // material *, peremnél NULL, a színhez tartozó anyagra mutat
    MezoTipus field;            // MezoTipus : FieldEl,FieldTherm,FieldElTherm, a tényleges mezõ típusával ÉS kapcsolatban van
    tomb<semiconductor> tsemi;  // tomb<semiconductor>, mely anyagok felé félvezetõ
    uns index;
    color() :is(false), tipus(SzinNormal), pmat(NULL), terfogat(0.0), index(0), field{ FieldElTherm } {}
};


struct model;


//***********************************************************************
struct z_a_hengerhez{
//***********************************************************************
    dbl ertek;      // dbl, a z irányú vastagság, ha nem henger, egyébként nem definiált
    bool henger_e;  // bool, ha henger koordinátarendszer, true
    z_a_hengerhez():ertek(nulla),henger_e(false){}
    dbl get(dbl sugar)const{ return henger_e ? 2*M_PI*sugar : ertek; }
};


//***********************************************************************
struct high_res_region_struct {
//***********************************************************************
    uns x1, y1, z1, x2, y2, z2;
    uns x_res, y_res, z_res;
};
//***********************************************************************


//***********************************************************************
struct model{
//***********************************************************************
    PLString fileName,name;     // PLString
    PLString mod_nev;           // a modellfájl neve kiterjesztés nélkül
    uns simdb;                  // uns, hány szimulációs fájl tartozik a modellhez
    uns dim;                    // uns, 1, 2 v 3.
    uns x_res,y_res,z_res;      // uns, 1..2048
    bool kerekcoord;            // bool, cartesian: false, sphere v. cylindrical: true
    bool is_35;                 // bool, ha a vsun 3.5 solver fusson
    bool is_half;               // bool, LED szimulációnál a W irányú falra nem küld sugarat
    dbl r_min,r_max;            // dbl, ha kerekcoord==true
    tomb<z_a_hengerhez> x_pit,y_pit,z_pit;  // tomb<z_a_hengerhez>, a henger_e csak z irányban lehet beállítva
    tomb<dbl> x_hossz,y_hossz,z_hossz;// tomb<dbl>, kétszer hosszabb, mint a pit, minden cella közepének és tetejének a távolságát is tárolja, 2*x a közép, 2*x+1 a tetõ
    tomb<bitmap> tbmp;          // tomb<bitmap>
    vezetes amb_emiss;          // vezetes, a környezet emissziós tényezõje, 0..1, default=1
    dbl general_reflectivity;   // dbl, az anyagok elnyelik vagy visszaverik a LED fényét, 0..1, default 0
    tomb<material> tmat;        // tomb<material>
    color tcolor[colmax+1];     // color tcolor[colmax]
    dbl A_semi_full;            // dbl, az összes félvezetõ felülete
    uns coupled_azon_start;     // 
    tomb<high_res_region_struct> righ_res_regions; // tomb<high_res_region_struct>
    model() :simdb(0), A_semi_full(nulla), general_reflectivity{ 0 }, is_half{ false }
    { tcolor[colmax].is = true; tcolor[colmax].tipus = SzinBoundary; amb_emiss.g[0] = amb_emiss.g[1] = amb_emiss.g[2] = egy; }
    // Ide kellenek, mert a félvezetõnél innen könnyebb számolni, meg a felületet is itt ismeri.
    //bool isCellExist(uns x,uns y,uns z)const; // a megadott koordinátákkal létezik-e cella (ha belsõ perem, akkor nem)
    void read(PLString path);   // beolvassa a fileName nevû modelfájlt
};


//***********************************************************************
struct convection{
//***********************************************************************
    bool is_defined;
    PLString nev;               // PLString
    vezetes radiation;          // vezetes, a sugárzásra számolt HTC mekkora részét vegye figyelembe, default=1
    ConvTipus vertical_tipus,lower_tipus,upper_tipus; // ConvTipus
    vezetes vertical_value,lower_value,upper_value;   // vezetes
    irany axis;                 // melyik tengely mentén forog, X_IRANY vagy Y_IRANY vagy Z_IRANY
    dbl angle;                  // dbl, fokban!
    ConvTipus edge;             // ConvHTC,ConvEdgeWei_H,ConvEdgeWei_I,ConvEdgeWei_HI, ha HTC, akkor nem számol élt.
    convection():is_defined(false),radiation(egy),vertical_tipus(ConvHTC),lower_tipus(ConvHTC)
        ,upper_tipus(ConvHTC),vertical_value(10.0),lower_value(10.0),upper_value(10.0),axis(X_IRANY),angle(90.0),edge(ConvHTC){}
};


//***********************************************************************
struct boundary{
//***********************************************************************
    PeremTipus tipus;           // PeremTipus: PeremOpen,PeremU,PeremR,PeremRU. A nonuniformot most nem implementálom, majd ha szükség lesz rá.
    dbl value, value2;          // dbl, a value2 csak conv_temp peremfeltétel esetén kap értéket, ez a csatolt hõmérséklet
    uns v6_index;               //
    bool is_special;            // ha helyfüggõ a peremfeltétel (convection és conv map esetén)
    convection conv;            // convection
    tomb2d<dbl> conv_map;       // tomb2d<dbl>
    const tomb3d<real_cell_res> *Ttomb,*Ptomb; // const tomb3d<real_cell_res> *
    dbl A,TA,P,Agyujt,TAgyujt,Pgyujt,ambiT; // dbl, A: felület, TA: hõmérséklet*terület, P teljesítmény, gyujt: a masina::get_all_matrix az ertek() segítségével ezekbe gyûjti az aktuális értékeket
    dbl rHTCA,rHTCAgyujt,cHTCA,cHTCAgyujt;
    boundary() :value(nulla), value2(nulla), Ttomb(NULL), Ptomb(NULL), ambiT(nulla), v6_index(0), is_special{ false } {}
};


//***********************************************************************
struct csomag{
//***********************************************************************
    uns szin;                   // uns
    boundary el[BASIC_SIDES];   // boundary el[BASIC_SIDES];
    boundary th[BASIC_SIDES];   // boundary th[BASIC_SIDES];
    csomag():szin(colmax){}     
};


//***********************************************************************
struct excitation{
//***********************************************************************
    bool is;                    // bool
    GerjTipus tipus;            // GerjSemmi,GerjU,GerjI
    dbl ertek;                  // dbl
    excitation():is(false){}
};


//***********************************************************************
struct probeT{
//***********************************************************************
    PLString cimke;             // PLString
    Oldal oldal;                // Oldal: WEST=1,EAST=2,SOUTH=3,NORTH=4,BOTTOM=5,TOP=6,CENTER=19
    uns x,y,z;                  // uns, map esetén x=0..2, ahol 0=x, 1=y, 2=z, az érték pedig y-ba kerül
    uns x2,y2,z2;               // uns
    uns current_type;           // uns, if current probe => 0=color, 1=volume, 2=side
};


//***********************************************************************
struct excitation_2{
//***********************************************************************
    bool is_el;                 // elektromos vagy termikus
    uns color_index;            // uns
    GerjTipus tipus;            // GerjSemmi,GerjU,GerjI
    dbl ertek;                  // dbl
    excitation_2():color_index(colmax){}
};


//***********************************************************************
struct change_time{
//***********************************************************************
    dbl time;                   // dbl
    dbl timestep;               // dbl, ha 0, akkor az adott idõpontban nem változik
    tomb<excitation_2> excit;   // tomb<excitation_2>, a megváltozó gerjesztések
    change_time():timestep(nulla){}
};
    

//***********************************************************************
struct time_and_change{
//***********************************************************************
    i32 change_index;           // i32, Ha az idõponthoz nem tartozik change: -1, egyébként a ctrl tömb indexe.
    dbl time;                   // Szimulációs idõpont.
    time_and_change():change_index(-1),time(nulla){}
};

//***********************************************************************
struct analysis{
//***********************************************************************
    PLString nev;               // PLString
    AnalizisTipus tipus;        // AnalizisTipus {AnalDC,AnalNDC,AnalAC,AnalLinTran,AnalLogTran,AnalBode,AnalIdo,AnalCtrl}
    dbl from,to;                // dbl
    uns step;                   // uns
    uns ndc_maxiter;            // uns
    dbl relhiba,ndc_I0;         // dbl
    tomb<change_time> ctrl;     // tomb<change_time>
    tomb<time_and_change> times;// tomb<time_and_change>, controlled analízisnél ebbe tároljuk el a szimulációs idõpontokat
    
    //****************************************************************
    void fill_times(){          // a times tömbot feltölti a beállított adatok alapján
    //****************************************************************
        uns i_ctrl = 0, i_times = 0;
        dbl akt_time = nulla; // csak a lépésköz szerint változik, a change_time szerint nem
        dbl akt_step = from;
        times.clear();
        time_and_change tc;
        
        for( i_times = 0; akt_time <= to + akt_step / 1000.0; i_times++ ){
            dbl dt = akt_step / 1000.0;
            if( ctrl.size() > i_ctrl && ctrl[i_ctrl].time < akt_time + dt ){
                // Ha a change_time kisebb az akt-nál, vagy egybeesik vele, majdnem ugyanazt kell csinálni.
                tc.change_index = i_ctrl;
                tc.time = ctrl[i_ctrl].time;
                times.puffer_add( i_times, tc );
                if(ctrl[i_ctrl].timestep > nulla){
                    akt_step = ctrl[i_ctrl].timestep;
                    akt_time = ctrl[i_ctrl].time;
                }
                if(ctrl[i_ctrl].time >= akt_time - dt) // egybeesik a change_time és az akt_time
                    akt_time += akt_step;
                 i_ctrl++;
           }
            else{
                tc.change_index = -1;
                tc.time = akt_time;
                times.puffer_add( i_times, tc );
                akt_time += akt_step;
            }
            if (i_times == 0){ // Betesz egy elsõ szimulációs idõpontot
                tc.change_index = -1;
                tc.time = akt_step*1e-6;
                times.puffer_add(1, tc);
                i_times++;
            }
        }
        times.resize(i_times);
    }
    analysis():step(0),ndc_maxiter(10),from(nulla),to(nulla),relhiba(0.005),ndc_I0(0.1){}
};

struct simulation;

//***********************************************************************
class uchannel{
//***********************************************************************
    friend class Gyuri_uchannel;
    const color * pszin;        // const color, pointer a modell::tcolor tömb elemére
    simulation * psim;          // az õt tartalmazó szimuláció címe
    uns szinindex;              // uns, a modell::tcolor tömb elemének indexe
    uns n;                      // uns, ennyi részre van osztva a csatorna
    CsatornaTipus tipus;        // CsatornaTipus: CsatRect, CsatCirc, CsatTrap
    bool is_auto_wall_temp;     // bool
    bool is_reverse;            // bool, a csatorna melyik végén fújunk be?
    dbl flow_rate;              // dbl, l/h
    dbl fixed_wall_temp;        // dbl, if(!is_auto_wall_temp) 
    dbl fluid_temp;             // dbl
    dbl width_bottom, width_top;// dbl
    dbl roughness;              // dbl
    dbl height;                 // dbl
    dbl density, dynamic_visc;  // dbl
    dbl spec_heat, heat_cond;   // dbl
    PLString nev;               // PLString
    Vec3d start, stop;          // Vec3d, a csatorna eleje és vége [m], már figyelembe van véve is_reverse
    dbl length;                 // dbl
public:
    void init1(simulation *psim, uns szinindex, const color *pszin){ this->psim = psim; this->szinindex = szinindex; this->pszin = pszin; }
    void init2(PLString path, PLString fajlnev, const PLString &cimke);
    void read(PLString path, PLString fajlnev);
    void set_start_stop_length();
    void set_uchannel_in_bmp();
    const PLString & getLabel()const{
        if (pszin == nullptr)
            throw hiba("uchannel::getLabel", "pszin==nullptr (using uninitialized object)");
        return pszin->nev;
    }
    uchannel() :pszin(nullptr), szinindex(colmax), n(8), tipus(CsatTrap), is_auto_wall_temp(true), is_reverse(false), 
        flow_rate(1.0), fixed_wall_temp(60.0), fluid_temp(nulla), width_bottom(250e-6), width_top(350e-6), height(67e-6),
        density(1.1614), dynamic_visc(1.84e-5), spec_heat(1005.0), heat_cond(0.0261), nev("undefined-channel-label"),
        length(nulla), roughness(1.0e-6){}
};


//***********************************************************************
class Gyuri_uchannel{
//***********************************************************************
    struct Rladder{
        double resistance_along;
        double resistance_cross;
        double heat_exchange;
        double Twall;
        double Tout;
    };

    struct Segment{
        double mass_flow_rate;
        double whole_area;
        double Nusselt_number;
        double local_heat_transfer_coeff;
    };

    struct Channel{
        //Channel geometries and properties
        double length;
        double height;
        double width_top;
        double width_bottom;
        double roughness;
        enum geometry { rect, circ, trap } cross_section_type; // 0 Rectangular, 1 circle, 2 trapezoid, etc.

        //Fluid properties, mechanical & thermal parameters
//        PLString fluid_mat;
        double density;
        double specific_heat; // J/(Kg*K)
        double heat_conductivity;
        double dynamic_viscosity;
        double fluid_inlet_temperature;

        //Flow properties
        double mass_flow_rate;
        double volumetric_flow_rate;
        double velocity_of_flow;
        enum type_of_flow { lam, turb, mixed } laminar_or_turbulent; 

        //Calculated values
        double area;
        double perimeter;
        double side_ratio;
        double hydraulic_diameter; //4*A/P
        double Prandtl_number;
        double Reynolds_number;
        double friction_factor;
        double head_loss;
        double press_drop;
        double press_drop_dV;
        double average_Nusselt_number;
        double global_heat_transfer_coefficient;
        double global_thermal_conductivity;

        Segment * segments;
        Rladder * stages;
        Channel() :segments(nullptr), stages(nullptr){}
        ~Channel(){ delete[] segments; delete[] stages; }
    };

};


//***********************************************************************
class powermap {
//***********************************************************************
    simulation *psim;
    bool van_e;

    enum Hol { top, volume };
    Hol hol;

    bool is_exact, is_stop;
    u32 x, y, z; // u32, resolution
    tomb<dbl> x_pitch, y_pitch, z_pitch;

    struct map_t {
        double t;
        tomb3d<dbl> pmap;
        void read(srfajl &fajl, uns &sor, u32 x, u32 y, u32 z);
    };
    tomb<map_t> tombok;

public:
    void init(simulation *psim) { this->psim = psim; }
    void read(PLString path, PLString fajlnev);
    bool is_exists()const { return van_e; }
    void rescale() { // a beolvasott felbontású mapet átskálázza a modellben megadott felbontáshoz
        if (!is_exact)
            throw hiba("powermap::rescale", "not exact powermap is not supported"); 
    }
    void get_map(dbl t, dbl * map); // t idõpontban érvényes disszipációtérképet adja (Wattban, nem sûrûségben)
    powermap() :psim(nullptr), van_e(false), hol(top), is_exact(false), is_stop(false),
        x(1), y(1), z(1) {}
};


//***********************************************************************
struct auto_transi {
//***********************************************************************
    bool is_V, is_T, is_no_plus_step_data;
    dbl V_max_rel, V_min_dt;
    dbl T_max_rel, T_min_dt;
    uns V_max_plus_steps, T_max_plus_steps; // ennyi plusz sikeres lépést enged két elõre megadott idõpont között; ha 0, bármennyit
    auto_transi() :is_V(false), is_T(false), is_no_plus_step_data(false), V_max_plus_steps(0), T_max_plus_steps(0) {}
};


class masina;

//***********************************************************************
struct simulation{
//***********************************************************************
    enum autotr_state{atr_normal, atr_auto, atr_auto_fixpont, atr_hisz_back};
    PLString fileName,name;     // PLString
    PLString fimdesc_name;
    PLString sim_nev;           // a fájl neve kiterjesztés nélkül
    model * pmodel;             // model * pmodel;
    MezoTipus mezo;             // MezoTipus mezo: FieldEl,FieldTherm,FieldElTherm
    uns mezo_szamitasi_mod;     // elektrotermikusnál, 1: hagyományos kétteres iteráció, 2: Newton-Raphson kétteres, 3: Newton-Raphson egyteres
    bool is_lin,is_no_semi,is_no_seebeck,is_no_peltier,is_peltier_center; // bool
    bool is_no_thomson, is_no_joule, is_no_accelerator, is_fim_txt; // bool
    bool is_vesszo,is_nofim,aramot,is_always_quad,is_always_strassen_mul;	// bool
    uns valostipus;             // uns: 0: double, 1: quad_float; elthermnél: 1: el=quad, 2: th=quad, 3: mindkettõ quad
    uns el_nonlin_subiter;      // uns, default 1
    uns ndc_miniter;            // uns, default 0
    uns cpu_threads;            // uns
    uns hdd;                    // uns
    uns optimize;               // uns, egyelõre nem definiált, hogy mire akarom használni
    uns FIM_res_xy,FIM_res_z;   // uns, FIM kényszerítése adott felbontásra (nem lehet kisebb, mint a tényleges)
    uns FIM_diff_x, FIM_diff_y, FIM_diff_z; // uns, eltolás a FIM fájlon belül.
    tomb<convection> tconv;     // tomb<convection> tconv;
    csomag normalperem;         // csomag, a színt nem használjuk
    tomb<csomag> tinner;        // tomb<csomag>, innerindexszel indexelni!!!
    tomb<boundary> peremlista_menteshez; // a beolvasott peemeket ide is bemásoljuk, a v6core számára innen mentünk
    uns innerindex[colmax];     // uns[colmax], segédtömb a tinner kereséshez
    powermap pmap;              // powermap, külsõ disszipációtérkép alkalmazása esetén
    tomb<dbl> prev_T_map, prev_U_map, akt_T_map, akt_U_map; // cellaközéppontok hõmérséklete/feszültsége az elõzõ és az aktuális tranziens/controlled lépés után, az automatikus léptetéshez használjuk. A run_transi/run_controlled resize-zolja
    auto_transi auto_tra;       // auto_transi, automatikus tranziens lépésköz paraméterei
    tomb<uchannel> tucha;       // tomb<uchannel>, a mikrocsatornák paraméterei
    tomb<csomag> tucha_perem[colmax]; // tomb<csomag> tucha_perem[colmax], minden uchannel színhez van csomagtömb, egy tömb annyi elemû, ahány darabra osztott a csatorna
    excitation texcitE[colmax]; // excitation texcitE[colmax], elektromos gerj
    tomb<excitation> mulE[colmax]; // tomb<excitation> mulE[colmax], ha több gerjesztés van
    excitation texcitT[colmax]; // excitation texcitT[colmax], hõ gerj
    tomb<excitation> mulT[colmax]; // tomb<excitation> mulT[colmax], ha több gerjesztés van
    uns db_temp;                // uns, beolvasáskor használt ideiglenes változó
    uns index_temp[4];          // uns[4], beolvasáskor használt ideiglenes változó
    dbl ambiT;                  // dbl
    tomb<dbl> mulAmbiT;         // tomb<dbl>, ha több külsõ hõmérséklet van megadva
    tomb<probeT> tproV;         // tomb<probeT>
    tomb<probeT> tproT;         // tomb<probeT>
    tomb<probeT> tproC;         // tomb<probeT>, áram probe, elektromos és hõáramra egyaránt
    tomb<probeT> tproM;         // tomb<probeT>, map
    tomb<probeT> tproS;         // tomb<probeT>, section
    tomb<analysis> tanal;       // tomb<analysis> tanal;
    simulation():pmodel(NULL),is_lin(false),is_no_semi(false),is_fim_txt(false),is_no_seebeck(false),is_no_peltier(false),
        is_peltier_center(false),is_no_thomson(false),is_no_joule(false),is_no_accelerator(false),el_nonlin_subiter(1),ndc_miniter(0),valostipus(0),is_always_strassen_mul(false),
        cpu_threads(32767),hdd(0),optimize(1),ambiT(nulla),is_vesszo(false),is_nofim(false),aramot(true),is_always_quad(false),
        FIM_res_xy(0), FIM_res_z(0), FIM_diff_x(0), FIM_diff_y(0), FIM_diff_z(0), mezo_szamitasi_mod{ 3 } {}
    void read(PLString path);   // beolvassa a fileName nevû szimulációs fájlt
    const boundary * get_perem(uns x, uns y, uns z, Oldal oldal, bool is_el)const;
    const csomag & get_inner_perem(uns x, uns y, uns z)const;
};

//***********************************************************************
struct t_modell_face_adat {
//***********************************************************************
    uns kulso_el_db, kulso_th_db;
    uns csatlakozo_index_el, csatlakozo_index_th;
    dbl A, L;
    uns anyag_index; // a megfelelõ irányú anyag
    uns el_perem_index, th_perem_index;
    uns th_perem_x, th_perem_y;
    char th_perem_c; // WESNBT
    bool is_el_perem, is_th_perem;
    uns junction_index;
    uns face_index_el, face_index_th;
    tomb2d<t_modell_face_adat> belso_facek; // ha az oldal több face-re van osztva
    t_modell_face_adat() :csatlakozo_index_el{ 0 }, csatlakozo_index_th{ 0 }, junction_index { 0 }, kulso_el_db{ 0 }, kulso_th_db{ 0 },
        face_index_el{ 0 }, face_index_th{ 0 }, th_perem_c{ 0 } {}
    void face2face(uns xx, uns yy, bool is_el, bool is_th); // saját adatai alapján létrehozza a belso face-eket
    void write_face(FILE *fp, bool is_el, bool is_th, uns & start_face_index, meret_tomb_tipus & meret_tomb, uns cella_anyag_index, uns cella_tipus);
    uns facetipus_azonosito(bool is_el, bool is_th, uns cella_anyag_index);
    void set_face_indexek(bool is_el, bool is_th, uns & start_face_index);
    void addSideToCurrentProbe(v6eredm& eredm, uns cella_index) const;
};


//***********************************************************************
struct ketpont {
//***********************************************************************
    dbl x0, x1, y0, y1, z0, z1;
    ketpont() :x0{ 0 }, x1{ 0 }, y0{ 0 }, y1{ 0 }, z0{ 0 }, z1{ 0 } {}
};

//***********************************************************************
struct t_modell_cella {
//***********************************************************************
    bool is_cella;  // Ha tényleg cella, nem perem/csatorna
    bool is_el, is_th;
    t_modell_face_adat face_adat[BASIC_SIDES];
    bool is_nonlin_el, is_nonlin_th;
    tomb3d<t_modell_cella> belso_cellak;
    uns color_index;
    uns anyag_index; // csak az irányfüggetlen tulajdonságok (=hõkapacitás) kiszámításához (amúgy az x irányú anyagot rendeljük hozzá)
    material * pmat; // fényporhoz kell
    bool is_junction; // ha van a cella valamelyik oldalán junction
    tomb<uns> fenyforras_cella_index; // ha van fényforrás párja, akkor az indexe ebben van. Egyelõre csak a vele egy oszlopban lévõk fényét kapja.
    uns junction_bottom_face, junction_top_face; // Ha top vagy bottom face junction, akkor ez az indexe: 2 lépésben beállítva
    //tomb<fenyut> kek_fenyutak, sarga_fenyutak; // ha ez fénypor cella, és figyelembe veszünk ide érkezõ fényt, akkor a fényutak
    dbl V; // térfogat
    ketpont sarkok; // csak nem high_res celláknál, fényporhoz kell
    dbl sum_P_kek, sum_P_sarga; // fényporhoz kell
    uns cella_index; // 0: ha nincs beállítva, >0 a tényleges cellánál
    t_modell_cella() :is_nonlin_el{ false }, is_nonlin_th{ false }, is_cella{ true }, cella_index{ 0 }, is_junction{ false }, 
        pmat{ nullptr }, junction_bottom_face{ 0 }, junction_top_face{ 0 }, sum_P_kek{ 0 }, sum_P_sarga{ 0 } {}

    void set_face_el(Oldal o, bool is_el_perem, uns el_perem_index, uns kulso_el_db, bool is_nonlin = false) { 
        face_adat[o].is_el_perem = is_el_perem;  
        face_adat[o].el_perem_index = el_perem_index;
        face_adat[o].kulso_el_db = kulso_el_db;
        if (is_nonlin)
            is_nonlin_el = true;
    }

    void set_face_th(Oldal o, bool is_th_perem, uns th_perem_index, uns kulso_th_db, uns x, uns y, char c, bool is_nonlin = false) {
        face_adat[o].is_th_perem = is_th_perem;
        face_adat[o].th_perem_index = th_perem_index;
        face_adat[o].th_perem_x = x;
        face_adat[o].th_perem_y = y;
        face_adat[o].th_perem_c = c;
        face_adat[o].kulso_th_db = kulso_th_db;
        if (is_nonlin)
            is_nonlin_th = true;
    }

    void mul_face(uns x_mul, uns y_mul, uns z_mul) {
        face_adat[WEST].kulso_el_db *= x_mul;
        face_adat[WEST].kulso_th_db *= x_mul;
        face_adat[EAST].kulso_el_db *= x_mul;
        face_adat[EAST].kulso_th_db *= x_mul;
        face_adat[SOUTH].kulso_el_db *= y_mul;
        face_adat[SOUTH].kulso_th_db *= y_mul;
        face_adat[NORTH].kulso_el_db *= y_mul;
        face_adat[NORTH].kulso_th_db *= y_mul;
        face_adat[BOTTOM].kulso_el_db *= z_mul;
        face_adat[BOTTOM].kulso_th_db *= z_mul;
        face_adat[TOP].kulso_el_db *= z_mul;
        face_adat[TOP].kulso_th_db *= z_mul;
    }

    void set_belso_cellak_mul(uns x_mul, uns y_mul, uns z_mul);
    void set_egy_belso_cella(const t_modell_cella & tulaj);
    uns get_szomszed_color(model & aktMod, uns x, uns y, uns z, Oldal old) {
        switch (old) {
            case WEST: x--; break;
            case EAST: x++; break;
            case SOUTH: y--; break;
            case NORTH: y++; break;
            case BOTTOM: z--; break;
            case TOP: z++; break;
        }
        return aktMod.tbmp[z].getpixel_also(x, y);
    }
    void set_cella_index(uns & startindex) {
        if (!is_cella)
            return;
        if (belso_cellak.size() > 0) {
            for (uns z = 0; z < belso_cellak.z_size(); z++)
                for (uns y = 0; y < belso_cellak.y_size(); y++)
                    for (uns x = 0;x < belso_cellak.x_size(); x++)
                        belso_cellak.getref(x, y, z).set_cella_index(startindex);
        }
        else
            cella_index = startindex++;
    }
    void write_cella(FILE *fp, FILE *descfp, meret_tomb_tipus & meret_tomb, uns cx, uns cy, uns cz); // ha tényleg cella, kiírja magát
    void set_face_indexek();
    void addSideToCurrentProbe(v6eredm& eredm, Oldal oldal) const;
private:
    uns write_faces(FILE *fp, meret_tomb_tipus & meret_tomb, uns cella_tipus); // a write_cella hívja
    uns cellatipus_azonosito(); // a write_cella hívja
};

//***********************************************************************
class doboz {
//***********************************************************************
public:
    bool is_lin;
    uns x1, x2, y1, y2, z1, z2;
    doboz(uns x1, uns x2, uns y1, uns y2, uns z1, uns z2) :x1{ x1 }, x2{ x2 }, y1{ y1 }, y2{ y2 }, z1{ z1 }, z2{ z2 }, is_lin{ false } {}
    doboz() :x1{ 0 }, x2{ 0 }, y1{ 0 }, y2{ 0 }, z1{ 0 }, z2{ 0 }, is_lin{ false } {}
    doboz(const doboz & m) :is_lin{ m.is_lin }, x1{ m.x1 }, x2{ m.x2 }, y1{ m.y1 }, y2{ m.y2 }, z1{ m.z1 }, z2{ m.z2 } {}
};

//***********************************************************************
class red_fa {
//***********************************************************************
public:
    bool is_alcella;
    doboz meret;
    doboz almeret; // ha alcella
    uns level;
    uns index; // kiíráskor az elágazás indexe, 1-tõl
    uns blokk_kezdo_index; // az ehhez a csomóponthoz tartozó részfa összes eleme egybefüggõ indextartományon van, ennek a tartományak az elsõ eleme ez, az utolsó elem az aktuális, azaz az indexszel jelölt.
    red_fa *bal, *jobb;
    const t_modell_cella * cella;
    tomb<mit_hova_masol> mit_hova;
    uns kozos_be1, kozos_be2, kozos_db;
    uns A0, B0;     // ha elemi cella, akkor a külsõ és redukálandó csomópontok száma. centroidnál figyelni a számolásra (A0++, B0--)!
    tol_db oldalak[BASIC_SIDES];
    Oldal kozos_oldal_1; // A bal cella mely oldala a közös? Csak EAST, NORTH és TOP lehet. (Ha nem így lenne, akkor a jobb-bal megcserélésével elõállítható ez az állapot.)
    red_fa(doboz meret, uns level) :bal{ nullptr }, jobb{ nullptr }, is_alcella{ false },
        meret{ meret }, level{ level }, cella{ nullptr }, A0{ 0 }, B0{ 0 } {}
    void optimized_add_mit_hova(const mit_hova_masol & miho);
     static red_fa * build_tree_uj(MezoTipus mt, const tomb3d<t_modell_cella> & r, doboz befoglalo_d, const tomb<Oldal> & iranyok, uns level = 0);
     static red_fa * build_tree_old(MezoTipus mt, const tomb3d<t_modell_cella> & r, doboz d, uns level = 0);
     static red_fa * build_subtree_uj(MezoTipus mt, const t_modell_cella & c, doboz befoglalo_d, const tomb<Oldal> & iranyok, uns level);
     static red_fa * build_subtree_old(MezoTipus mt, const t_modell_cella & c, doboz d, uns level);
     static void decrease_level(red_fa * gy);
     static void optimize_tree(red_fa * & gy);
     static void tol_ig_db_szamol(MezoTipus mt, red_fa * gy);
     static void feloszt(MezoTipus mt, const tomb3d<t_modell_cella> & r, const doboz & d, doboz & ki_1, doboz & ki_2, Oldal & kozos_oldal_1);
     static bool szukito(MezoTipus mt, const tomb3d<t_modell_cella> & r, doboz & d);
     static void iranyszamolo_uj(const doboz & be_meret, tomb<Oldal> & iranyok);
     static uns set_indexek(red_fa * gy, uns & start_index);
     static void write_tree(FILE *fp, red_fa * gy);
};


//***********************************************************************
struct v6gerj {
//***********************************************************************
    char tipus, ter; // tipus: 'U', 'I', 'T', 'P', 'N'; ter: 'E', 'T'
    dbl ertek;
    uns color_index;
};


//***********************************************************************
struct v6eredm {
//***********************************************************************
    struct CurrIndex {
        uns cella_index;
        uns face_index;
        CurrIndex(uns ci = 0, uns fi = 0) :cella_index{ ci }, face_index{ fi }{}
    };
    char eredm_fajta;           // 'P' vagy 'M' vagy 'C'
    char probe_map_fajta;       // 'E', 'T', 'F', 'I', 'L', 'R'
    uns probe_cella_index;      // hányas cella
    uns probe_face_index;       // ha 'F' a probe, akkor a face indexe
    tomb<CurrIndex> elCurrentProbe; // faces of the electrical current probe if this is a current probe
    tomb<CurrIndex> thCurrentProbe; // faces of the thermal current probe if this is a current probe
    v6eredm() {}
};


//***********************************************************************
struct v6anal {
//***********************************************************************
    char tipus; // 'D', 'A', 'S', 'T': DC, AC, step, timeconst
    dbl ertek; // 'A': Hz, 'S': s
    bool is_del_excits; // törli az összes korábbi gerjesztést, természetesen az aktuálisakat nem
    bool is_reset; // minden kiindulási hõmérséklet, fesz, áram 0, azaz mint az elsõ analízis elõtt
    uns maxiter; // Ha 0, akkor nem változik
    dbl maxhiba; // Ha 0, akkor nem változik
    dbl I0; // Ha 0, akkor nem változik
    tomb<v6gerj> gerj;
    tomb<v6eredm> eredm;
    v6anal() :tipus{ 'D' }, ertek{ 0 }, is_del_excits{ false }, is_reset{ false }, maxiter{ 0 }, maxhiba{ 0 }, I0{ 0 } {}
    void reset() { ertek = 0; is_del_excits = false; is_reset = false; maxhiba = 0; maxiter = 0; I0 = 0; gerj.clear(); eredm.clear(); } // csak az analízisítpus marad
};

//***********************************************************************
class apa{
//***********************************************************************
    PLString path,projFile;     // projFile tartalmazza az útvonalat is
    PLString proj_nev;          // a fájl neve kiterjesztés nélkül
    tomb<model> tmodels;        // tomb<model> tmodels;
    tomb<simulation> tsim;      // tomb<simulation> tsim;
    tomb3d<t_modell_cella> modell_racs;  // az aktuális rács
    red_fa *modell_fa_el, *modell_fa_th, *modell_fa_elth; // törölni kell használat elõtt!
    tomb<v6anal> aktAnalizisek; // törölni kell használat elõtt!
    meret_tomb_tipus meret_tomb; // a térfogatok, felületek és hosszak gyûjteménye
    sugarfeldolgozo sugar_feldolgozo; // a fénypor sugarak számításához, index_elemi_cells hívja

    //***********************************************************************
    void face_par_index_beallito(t_modell_face_adat & face_1, t_modell_face_adat & face_2, uns & akt_csatlakozo_index) {
    //***********************************************************************
        const char * fvnev = "apa::face_par_index_beallito";
        if (face_1.kulso_el_db == 1) {
            if (face_2.kulso_el_db != 1)throw hiba(fvnev, "face_2.kulso_el_db != 1 (%u)", face_2.kulso_el_db);
            face_1.csatlakozo_index_el = akt_csatlakozo_index++;
            face_2.csatlakozo_index_el = akt_csatlakozo_index++;
        }
        if (face_1.kulso_th_db == 1) {
            if (face_2.kulso_th_db != 1)throw hiba(fvnev, "face_2.kulso_th_db != 1 (%u)", face_2.kulso_th_db);
            face_1.csatlakozo_index_th = akt_csatlakozo_index++;
            face_2.csatlakozo_index_th = akt_csatlakozo_index++;
        }
        if (face_1.kulso_el_db>1)throw hiba(fvnev, "face_1.kulso_el_db>1 (%u)", face_1.kulso_el_db);
        if (face_2.kulso_el_db>1)throw hiba(fvnev, "face_2.kulso_el_db>1 (%u)", face_2.kulso_el_db);
        if (face_1.kulso_th_db>1)throw hiba(fvnev, "face_1.kulso_th_db>1 (%u)", face_1.kulso_th_db);
        if (face_2.kulso_th_db>1)throw hiba(fvnev, "face_2.kulso_el_db>1 (%u)", face_2.kulso_th_db);
    }
    //***********************************************************************
    void face_index_beallito(t_modell_cella & cella_1, t_modell_cella & cella_2, uns & akt_csatlakozo_index, irany ir) {
    //***********************************************************************
        Oldal oldal_1, oldal_2;
        uns dx = 0, dy = 0, dz = 0;
        uns xx_1 = 0, xx_2 = 0, yy_1 = 0, yy_2 = 0, max_xx = 0, max_yy = 0;
        if (ir == X_IRANY) {
            oldal_1 = EAST;
            oldal_2 = WEST;
            dx = 1;
            xx_1 = cella_1.belso_cellak.y_size();
            xx_2 = cella_2.belso_cellak.y_size();
            yy_1 = cella_1.belso_cellak.z_size();
            yy_2 = cella_2.belso_cellak.z_size();
        }
        else if (ir == Y_IRANY) {
            oldal_1 = NORTH;
            oldal_2 = SOUTH;
            dy = 1;
            xx_1 = cella_1.belso_cellak.x_size();
            xx_2 = cella_2.belso_cellak.x_size();
            yy_1 = cella_1.belso_cellak.z_size();
            yy_2 = cella_2.belso_cellak.z_size();
        }
        else {
            oldal_1 = TOP;
            oldal_2 = BOTTOM;
            dz = 1;
            xx_1 = cella_1.belso_cellak.x_size();
            xx_2 = cella_2.belso_cellak.x_size();
            yy_1 = cella_1.belso_cellak.y_size();
            yy_2 = cella_2.belso_cellak.y_size();
        }
        max_xx = xx_1 > xx_2 ? xx_1 : xx_2;
        max_yy = yy_1 > yy_2 ? yy_1 : yy_2;
        uns face_per_oldal_x_1 = xx_1 == 0 ? 0 : max_xx / xx_1;
        uns face_per_oldal_x_2 = xx_2 == 0 ? 0 : max_xx / xx_2;
        uns face_per_oldal_y_1 = yy_1 == 0 ? 0 : max_yy / yy_1;
        uns face_per_oldal_y_2 = yy_2 == 0 ? 0 : max_yy / yy_2;
        bool is_f1_osztott = (xx_1 != max_xx || yy_1 != max_yy);
        bool is_f2_osztott = (xx_2 != max_xx || yy_2 != max_yy);
        if (cella_1.belso_cellak.size() == 0 && cella_2.belso_cellak.size() == 0)
            face_par_index_beallito(cella_1.face_adat[oldal_1], cella_2.face_adat[oldal_2], akt_csatlakozo_index);
        else if (cella_1.belso_cellak.size() == 0) { // A cella_2 High.res, a cella_1 nem. cella_2 belsõ celláinak nem lehet belsõ face-e
            for (uns i = 0; i < yy_2; i++)
                for (uns j = 0; j < xx_2; j++) { // cella_2 WEST/SOUTH/BOTTOM
                    if (ir == X_IRANY)face_par_index_beallito(cella_1.face_adat[oldal_1].belso_facek.getref(j, i), cella_2.belso_cellak.getref(0, j, i).face_adat[oldal_2], akt_csatlakozo_index);
                    if (ir == Y_IRANY)face_par_index_beallito(cella_1.face_adat[oldal_1].belso_facek.getref(j, i), cella_2.belso_cellak.getref(j, 0, i).face_adat[oldal_2], akt_csatlakozo_index);
                    if (ir == Z_IRANY)face_par_index_beallito(cella_1.face_adat[oldal_1].belso_facek.getref(j, i), cella_2.belso_cellak.getref(j, i, 0).face_adat[oldal_2], akt_csatlakozo_index);
                }
        }
        else if (cella_2.belso_cellak.size() == 0) { // A cella_1 High.res, a cella_2 nem. cella_1 belsõ celláinak nem lehet belsõ face-e
            for (uns i = 0; i < yy_1; i++)
                for (uns j = 0; j < xx_1; j++) { // cella_1 EAST/NORTH/TOP
                    if (ir == X_IRANY)face_par_index_beallito(cella_1.belso_cellak.getref(cella_1.belso_cellak.x_size() - 1, j, i).face_adat[oldal_1], cella_2.face_adat[oldal_2].belso_facek.getref(j, i), akt_csatlakozo_index);
                    if (ir == Y_IRANY)face_par_index_beallito(cella_1.belso_cellak.getref(j, cella_1.belso_cellak.y_size() - 1, i).face_adat[oldal_1], cella_2.face_adat[oldal_2].belso_facek.getref(j, i), akt_csatlakozo_index);
                    if (ir == Z_IRANY)face_par_index_beallito(cella_1.belso_cellak.getref(j, i, cella_1.belso_cellak.z_size() - 1).face_adat[oldal_1], cella_2.face_adat[oldal_2].belso_facek.getref(j, i), akt_csatlakozo_index);
                }
        }
        else { // Mindkét cella High.res
            for (uns i = 0; i < max_yy; i++) {
                uns c1_belso_y = i / face_per_oldal_y_1;
                uns f1_belso_y = i % face_per_oldal_y_1;
                uns c2_belso_y = i / face_per_oldal_y_2;
                uns f2_belso_y = i % face_per_oldal_y_2;
                for (uns j = 0; j < max_xx; j++) {
                    uns c1_belso_x = j / face_per_oldal_x_1;
                    uns f1_belso_x = j % face_per_oldal_x_1;
                    uns c2_belso_x = j / face_per_oldal_x_2;
                    uns f2_belso_x = j % face_per_oldal_x_2;
                    if (ir == X_IRANY) {
                        t_modell_cella & c1 = cella_1.belso_cellak.getref(cella_1.belso_cellak.x_size() - 1, c1_belso_x, c1_belso_y);
                        t_modell_cella & c2 = cella_2.belso_cellak.getref(0, c2_belso_x, c2_belso_y);
                        t_modell_face_adat & face_1 = is_f1_osztott ? c1.face_adat[oldal_1].belso_facek.getref(f1_belso_x, f1_belso_y) : c1.face_adat[oldal_1];
                        t_modell_face_adat & face_2 = is_f2_osztott ? c2.face_adat[oldal_2].belso_facek.getref(f2_belso_x, f2_belso_y) : c2.face_adat[oldal_2];
                        face_par_index_beallito(face_1, face_2, akt_csatlakozo_index);
                    }
                    else if (ir == Y_IRANY) {
                        t_modell_cella & c1 = cella_1.belso_cellak.getref(c1_belso_x, cella_1.belso_cellak.y_size() - 1, c1_belso_y);
                        t_modell_cella & c2 = cella_2.belso_cellak.getref(c2_belso_x, 0, c2_belso_y);
                        t_modell_face_adat & face_1 = is_f1_osztott ? c1.face_adat[oldal_1].belso_facek.getref(f1_belso_x, f1_belso_y) : c1.face_adat[oldal_1];
                        t_modell_face_adat & face_2 = is_f2_osztott ? c2.face_adat[oldal_2].belso_facek.getref(f2_belso_x, f2_belso_y) : c2.face_adat[oldal_2];
                        face_par_index_beallito(face_1, face_2, akt_csatlakozo_index);
                    }
                    else {
                        t_modell_cella & c1 = cella_1.belso_cellak.getref(c1_belso_x, c1_belso_y, cella_1.belso_cellak.z_size() - 1);
                        t_modell_cella & c2 = cella_2.belso_cellak.getref(c2_belso_x, c2_belso_y, 0);
                        t_modell_face_adat & face_1 = is_f1_osztott ? c1.face_adat[oldal_1].belso_facek.getref(f1_belso_x, f1_belso_y) : c1.face_adat[oldal_1];
                        t_modell_face_adat & face_2 = is_f2_osztott ? c2.face_adat[oldal_2].belso_facek.getref(f2_belso_x, f2_belso_y) : c2.face_adat[oldal_2];
                        face_par_index_beallito(face_1, face_2, akt_csatlakozo_index);
                    }
                }
            }
        }
    }
public:
    apa(const char * ProjectFile);
    ~apa() { del_modell_fa(modell_fa_el); del_modell_fa(modell_fa_th); del_modell_fa(modell_fa_elth); }
    void write_v6sim();
    void write_akt_sim(FILE *fp, simulation & aktSim, analysis & aktAnal, const PLString & leiras, uns akt_n);
    void write_materials(FILE *fp, simulation & aktSim);
    void write_colors(FILE *fp, simulation & aktSim);
    uns write_junctions(FILE *fp, simulation & aktSim);
    void write_special_boundaries(FILE *fp, simulation & aktSim);
    void write_boundary_conditions(FILE *fp, simulation & aktSim);
    uns build_modell_racs(simulation & aktSim, uns & csatlakozo_face_db);
    void build_modell_fa(simulation & aktSim);
    void build_analizisek(simulation & aktSim);
    // TODO: fix hõmérsékletû/feszültségû gerjesztés mint perem 
    // TODO: mindenütt WESNBT az oldalsorrend. Ha a redukciós fában WSBENT-ra lehetne cserélni, akkor 10 helyett 8 másolással 
    // megoldható lenne a szétmásolás. (Nem egy nagy nyereség idõben...)
    void write_analizisek(FILE *fp, simulation & aktSim, uns junction_db);
    void index_and_write_elemi_cells(FILE *fp, simulation & aktSim, uns cellaszam, uns csatlakozo_db);
    void write_elemi_cells(FILE *fp, simulation & aktSim, uns cellaszam, uns csatlakozo_db);
    void write_meretek(FILE *fp);
    void write_modell_tree(FILE *fp, simulation & aktSim);
    void del_modell_fa(red_fa * & fa);
};


#endif
