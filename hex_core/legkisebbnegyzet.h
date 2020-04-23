//***********************************************************************
// legkisebb négyzet számoló header
// Creation date:  2019. 02. 14.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef LEGKISEBBNEGYZET_HEADER
#define	LEGKISEBBNEGYZET_HEADER
//***********************************************************************


//***********************************************************************
#include "kozos.h"
#include "hiba.h"
#include "matrix.hpp"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
struct trio {
//***********************************************************************
    rvt U, T, I; // nem feltétlenül U, T, I, hanem két paraméter meg egy érték
    trio(rvt U = 0, rvt T = 0, rvt I = 0) :U(U), T(T), I(I) {}
};


//***********************************************************************
struct polinom {
//***********************************************************************
    matrix<rvt> pars; // sor: U hatványai, oszlop: T hatványai
    rvt min_U, max_U; // ebben a tartományban illeszt a polinomra
    rvt min_T, max_T; // ebben a tartományban illeszt a polinomra
    illesztes_tipusa ill_eleje, ill_vege;
    rvt szorzo; // a polinom nem használja, csak tárolja! a visszaadott eredmények felhasználásánál kell figelembe venni
    polinom() : min_U{ 0 }, max_U{ 0 }, min_T{ 0 }, max_T{ 0 }, ill_eleje{ it_none }, ill_vege{ it_none }, szorzo{ 1 } {}
private:
    //***********************************************************************
    rvt get_I_poli(rvt U, rvt T) const{
    //***********************************************************************
        rvt ret = 0, nU = 1;
        for (uns i = 0; i < pars.get_row(); i++, nU *= U) {
            rvt nT = 1;
            for (uns j = 0; j < pars.get_col(); j++, nT *= T)
                ret += pars[i][j] * nU * nT;
        }
        return ret;
    }
    //***********************************************************************
    rvt get_dI_per_dU_poli(rvt U, rvt T) const{
    //***********************************************************************
        rvt ret = 0, nU = 1;
        for (uns i = 1; i < pars.get_row(); i++, nU *= U) {
            rvt nT = 1;
            for (uns j = 0; j < pars.get_col(); j++, nT *= T)
                ret += pars[i][j] * nU * nT * i;
        }
        return ret;
    }
    //***********************************************************************
    rvt get_dI_per_dT_poli(rvt U, rvt T) const{
    //***********************************************************************
        rvt ret = 0, nU = 1;
        for (uns i = 0; i < pars.get_row(); i++, nU *= U) {
            rvt nT = 1;
            for (uns j = 1; j < pars.get_col(); j++, nT *= T)
                ret += pars[i][j] * nU * nT * j;
        }
        return ret;
    }
    //***********************************************************************
    Itrio get_Itrio_poli(rvt U, rvt T) const{
    //***********************************************************************
        Itrio ret;
        ret.I = get_I_poli(U, T);
        ret.dI_per_dU = get_dI_per_dU_poli(U, T);
        ret.dI_per_dT = get_dI_per_dT_poli(U, T);
        return ret;
    }
    //***********************************************************************
    rvt get_I_lin_1(rvt U, rvt T) const {
    // Az origóba menõ egyenes a min_U pontból
    //***********************************************************************
        return get_I_poli(min_U, T) * U / min_U;
    }
    //***********************************************************************
    rvt get_dI_per_dU_lin_1(rvt U, rvt T) const {
    //***********************************************************************
        return get_I_poli(min_U, T) / min_U;
    }
    //***********************************************************************
    rvt get_dI_per_dT_lin_1(rvt U, rvt T) const {
    //***********************************************************************
        return get_dI_per_dT_poli(min_U, T) * U / min_U;
    }
    //***********************************************************************
    Itrio get_Itrio_lin_1(rvt U, rvt T) const{
    //***********************************************************************
        Itrio ret;
        rvt I0 = get_I_poli(min_U, T);
        ret.dI_per_dU = I0 / min_U;
        ret.I = ret.dI_per_dU * U;
        ret.dI_per_dT = get_dI_per_dT_lin_1(U, T);
        return ret;
    }
    //***********************************************************************
    rvt get_I_lin_2(rvt U, rvt T) const {
    // A max_U pontból az ott érvényes meredekséggel tovább menõ egyenes
    //***********************************************************************
        rvt I0 = get_I_poli(max_U, T);
        rvt dIdU0 = get_dI_per_dU_poli(max_U, T);
        rvt b = I0 - dIdU0 * max_U;
        return dIdU0 * U + b;
    }
    //***********************************************************************
    rvt get_dI_per_dU_lin_2(rvt U, rvt T) const {
    //***********************************************************************
        return get_dI_per_dU_poli(max_U, T);
    }
    //***********************************************************************
    rvt get_dI_per_dT_lin_2(rvt U, rvt T) const {
    //***********************************************************************
        return get_dI_per_dT_poli(max_U, T);
    }
    //***********************************************************************
    Itrio get_Itrio_lin_2(rvt U, rvt T) const{
    //***********************************************************************
        Itrio ret;
        ret.I = get_I_lin_2(U, T);
        ret.dI_per_dU = get_dI_per_dU_lin_2(U, T);
        ret.dI_per_dT = get_dI_per_dT_lin_2(U, T);
        return ret;
    }
    //***********************************************************************
    rvt get_I_strong_1(rvt U, rvt T) const {
    // I(U,T) = I0(U0, T)*(U/U0)^p
    // dI/dU = I(U,T) * p / U => p = U0/I0*dI/dU
    //***********************************************************************
        rvt I0 = get_I_poli(min_U, T);
        rvt dIdU0 = get_dI_per_dU_poli(min_U, T);
        rvt hatvany = I0 / min_U * dIdU0;
        if (hatvany < 6)
            hatvany = 6;
        return I0*pow(U / min_U, hatvany);
    }
    //***********************************************************************
    rvt get_dI_per_dU_strong_1(rvt U, rvt T) const {
    //***********************************************************************
        rvt I0 = get_I_poli(min_U, T);
        rvt dIdU0 = get_dI_per_dU_poli(min_U, T);
        rvt hatvany = I0 / min_U * dIdU0;
        if (hatvany < 6)
            hatvany = 6;
        rvt I = I0*pow(U / min_U, hatvany);
        return I * hatvany / U;
    }
    //***********************************************************************
    rvt get_dI_per_dT_strong_1(rvt U, rvt T) const {
    //***********************************************************************
        rvt I0 = get_I_poli(min_U, T);
        rvt dIdU0 = get_dI_per_dU_poli(min_U, T);
        rvt hatvany = I0 / min_U * dIdU0;
        if (hatvany < 6)
            hatvany = 6;
        return get_dI_per_dT_poli(min_U, T)*pow(U / min_U, hatvany);
    }
    //***********************************************************************
    Itrio get_Itrio_strong_1(rvt U, rvt T) const{
    //***********************************************************************
        Itrio ret;
        rvt I0 = get_I_poli(min_U, T);
        rvt dIdU0 = get_dI_per_dU_poli(min_U, T);
        rvt hatvany = I0 / min_U * dIdU0;
        if (hatvany < 6)
            hatvany = 6;
        rvt szorzo = pow(U / min_U, hatvany);
        ret.I = I0 * szorzo;
        ret.dI_per_dU = ret.I * hatvany / U;
        ret.dI_per_dT = get_dI_per_dT_poli(min_U, T) * szorzo;
        return ret;
    }
public:
/*
    //***********************************************************************
    rvt get_I(rvt U, rvt T, bool is_abs = true) const {
    //***********************************************************************
        set_hiba_hol h("polinom::get_I");
        rvt absU = is_abs ? abs(U) : U;
        if (absU < min_U) {
            if (T < min_T) {

            }
            else if (T > max_T) {

            }
            else {

            }
        }
        else if (absU > max_U) {
            if (T < min_T) {

            }
            else if (T > max_T) {

            }
            else {

            }
        }
        else {
            if (T < min_T) {

            }
            else if (T > max_T) {

            }
            else {
                return get_I_poli(absU, T);
            }
        }
    }
*/
    //***********************************************************************
    Itrio get_Itrio(rvt U, rvt T, bool is_abs = true) const {
    //***********************************************************************
        set_hiba_hol h("polinom::get_Itrio");
        rvt absU = is_abs ? abs(U) : U;
        Itrio ret;
        ret.dI_per_dT = 0; // 9-bõl 6 tartományban ez a helyzet
        if (absU < min_U) {
            rvt G = get_I_poli(min_U, T < min_T ? min_T : T > max_T ? max_T : T) / min_U;
            ret.I = G * absU;
            ret.dI_per_dU = G;
        }
        else if (absU > max_U) {
            rvt G = get_I_poli(max_U, T < min_T ? min_T : T > max_T ? max_T : T) / max_U;
            ret.I = G * absU;
            ret.dI_per_dU = G;
        }
        else {
            if (T < min_T) {
                ret.I = get_I_poli(absU, min_T);
                ret.dI_per_dU = ret.I / absU;//get_dI_per_dU_poli(absU, min_T);
            }
            else if (T > max_T) {
                ret.I = get_I_poli(absU, max_T);
                ret.dI_per_dU = ret.I / absU;//get_dI_per_dU_poli(absU, max_T);
            }
            else {
                ret.I = get_I_poli(absU, T);
                ret.dI_per_dU = get_dI_per_dU_poli(absU, T);
                ret.dI_per_dT = get_dI_per_dT_poli(absU, T);
            }
        }
        return ret;
    }

    //***********************************************************************
    rvt calc_U(rvt I, rvt T) {
    // Felezõ módszer, csak elõszámításhoz
    //***********************************************************************
        rvt ya = get_Itrio(min_U, T).I - I;
        rvt yb = get_Itrio(max_U, T).I - I;
        if (ya*yb > 0) {
            throw hiba("polinom::calc_U", "Azonos elojel");
        }
        rvt U1 = ya < 0 ? min_U : max_U;
        rvt U2 = ya < 0 ? max_U : min_U;
        while (abs((U1 - U2) / U2) > 1.0e-04) {
            rvt Uc = (U1 + U2)*0.5;
            rvt Ic = get_Itrio(Uc, T).I - I;
            if (Ic < 0)
                U1 = Uc;
            else
                U2 = Uc;
        }
        return U1;
    }
/*
    //***********************************************************************
    rvt get_I(rvt U, rvt T, bool is_abs = true) const {
    //***********************************************************************
        set_hiba_hol h("polinom::get_I");
        rvt absU = is_abs ? abs(U) : U;
        if (absU < min_U) {
            switch (ill_eleje) {
                case it_none:
                    return get_I_poli(absU, T);
                    break;
                case it_lin:
                    return get_I_lin_1(absU, T);
                    break;
                case it_strong:
                    return get_I_strong_1(absU, T);
                    break;
                default:
                    throw hiba(1, "program error: unknown fitting type");
            }
        }
        else if (absU > max_U) {
            switch (ill_vege) {
                case it_none:
                    return get_I_poli(absU, T);
                    break;
                case it_lin:
                    return get_I_lin_2(absU, T);
                    break;
                case it_strong:
                    throw hiba(1, "strong fitting not supported on the high end");
                    break;
                default:
                    throw hiba(1, "program error: unknown fitting type");
            }
        }
        else {
            return get_I_poli(absU, T);
        }
    }
    //***********************************************************************
    Itrio get_Itrio(rvt U, rvt T, bool is_abs = true) const {
    //***********************************************************************
        set_hiba_hol h("polinom::get_Itrio");
        rvt absU = is_abs ? abs(U) : U;
        if (absU < min_U) {
            switch (ill_eleje) {
                case it_none:
                    return get_Itrio_poli(absU, T);
                    break;
                case it_lin:
                    return get_Itrio_lin_1(absU, T);
                    break;
                case it_strong:
                    return get_Itrio_strong_1(absU, T);
                    break;
                default:
                    throw hiba(1, "program error: unknown fitting type");
            }
        }
        else if (absU > max_U) {
            switch (ill_vege) {
                case it_none:
                    return get_Itrio_poli(absU, T);
                    break;
                case it_lin:
                    return get_Itrio_lin_2(absU, T);
                    break;
                case it_strong:
                    throw hiba(1, "strong fitting not supported on the high end");
                    break;
                default:
                    throw hiba(1, "program error: unknown fitting type");
            }
        }
        else {
            return get_Itrio_poli(absU, T);
        }
    }
*/
};


//***********************************************************************
void lkn_illeszto(matrix<rvt> & eredmeny, uns U_hatvany_sor, uns T_hatvany_oszlop, const vektor<trio> & vuti);
//***********************************************************************


//***********************************************************************
struct illesztendo_adatok {
//***********************************************************************
    enum lsfit_egyenlet_tipus { lse_polinom };
    vektor<trio> vuti;              // illesztett esetben a mérési pontok (U-ra és T-re illeszti I-t), ha esetleg kellene, megõrzi
    rvt  szorzo;
    uns egyenlet_unspar_1, egyenlet_unspar_2; // polinomnál az U és T polinom foka
    bool is_trio; // duó vagy trió? azaz egy vagy két adatra kell fittelni a harmadikat
    illesztes_tipusa start, stop;
    lsfit_egyenlet_tipus egyenlet;
    //***********************************************************************
    illesztendo_adatok() :szorzo{ 1 }, egyenlet_unspar_1{ 0 }, egyenlet_unspar_2{ 0 }, is_trio{ false },
        start{ it_none }, stop{ it_none }, egyenlet{ lse_polinom } {}
    //***********************************************************************
    void beolvas_fajlbol();
    //***********************************************************************
    void set_polinom_by_fitting(polinom & result) {
    //***********************************************************************
        lkn_illeszto(result.pars, egyenlet_unspar_1, egyenlet_unspar_2, vuti);
        result.ill_eleje = start;
        result.ill_vege = stop;
        result.min_U = result.max_U = vuti[0].U;
        result.min_T = result.max_T = vuti[0].T;
        for (uns i = 0; i < vuti.size(); i++) {
            if (vuti[i].U < result.min_U)
                result.min_U = vuti[i].U;
            if (vuti[i].U > result.max_U)
                result.max_U = vuti[i].U;
            if (vuti[i].T < result.min_T)
                result.min_T = vuti[i].T;
            if (vuti[i].T > result.max_T)
                result.max_T = vuti[i].T;
        }
    }
    //***********************************************************************
};


//***********************************************************************
class lkn_2_par {
//***********************************************************************
    vektor<trio> vuti;
public:
    void clear() { vuti.clear(); }
    void push_back(const trio & uti) { vuti.push_back(uti); }
    void calc_parameters(matrix<rvt> & eredmeny, uns U_hatvany_sor, uns T_hatvany_oszlop) const {
        lkn_illeszto(eredmeny, U_hatvany_sor, T_hatvany_oszlop, vuti);
    }
};


}
#endif
