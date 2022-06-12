//***********************************************************************
// bemeneti adatok header
// Creation date:  2018. 06. 28.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef BEMENET_HEADER
#define	BEMENET_HEADER
//***********************************************************************


//***********************************************************************
#include "kozos.h"
#include "hiba.h"
#include "vektor.hpp"
#include "adattipusok.h"
#include "legkisebbnegyzet.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
class fazisvalto {
//***********************************************************************
    bool is;
    rvt TH1, TH2, SZ, energia;              // alapmennyiségek
    rvt TH1_tolt, TH2_tolt;                 // tolt alapmennyiségek
    rvt TK_tolt, TV_tolt, TN_tolt, TM, ITM; // származtatott mennyiségek
public:
    //***********************************************************************
    fazisvalto() :is{ false } {}
    //***********************************************************************
    bool is_fazisvalto()const { return is; }
    //***********************************************************************
    void set_params(rvt Th1, rvt Th2, rvt Sz, rvt Energia) { is = true; TH1 = Th1; TH2 = Th2; SZ = Sz; energia = Energia; }
    //***********************************************************************
    rvt get_F1()const { return TH1 - SZ; }
    //***********************************************************************
    rvt get_F2()const { return TH2 + SZ; }
    //***********************************************************************
    rvt get_L1()const { return TH1 + SZ; }
    //***********************************************************************
    rvt get_L2()const { return TH2 - SZ; }
    //***********************************************************************

    //***********************************************************************
    void set_Tamb(rvt Tamb) {
    //***********************************************************************
        TH1_tolt = TH1 - Tamb;
        TH2_tolt = TH2 - Tamb;
        TK_tolt = TH1_tolt - SZ;
        TV_tolt = TH2_tolt + SZ;
        TN_tolt = TH1_tolt + SZ;
        TM = TH2_tolt - TH1_tolt;
        ITM = rvt(1.0) / TM;
    }

    //***********************************************************************
    ertek_t H(rvt TB_tolt, rvt HA) const{
    //***********************************************************************
        if (TB_tolt <= TK_tolt)
            return ertek_t(rvt(), rvt());
        if (TB_tolt >= TV_tolt)
            return ertek_t(rvt(1.0), rvt());
        double TC_tolt = TH1_tolt + HA * TM;
        if (TB_tolt < TC_tolt - SZ)
            return ertek_t((TB_tolt - TK_tolt) * ITM, ITM);
        if (TB_tolt <= TC_tolt + SZ)
            return ertek_t(HA, rvt());
        return ertek_t((TB_tolt - TN_tolt) * ITM, ITM);
    }

    //***********************************************************************
    ertek_t CBm(rvt TA, rvt TB, rvt HA, rvt HB, rvt CthA, rvt CthB, const ertek_t & prev_CBm) const {
    //***********************************************************************
        if (TA == TB)
            return ertek_t((HA == rvt() || HA == rvt(1)) ? CthA : (CthA + energia * ITM), rvt());
        rvt nev = 1 / (TB - TA);
        rvt uj_CBm = energia * abs((HB - HA) * nev) + (CthA + CthB) * rvt(0.5);
        rvt deriv = (uj_CBm - prev_CBm.ertek) * nev; // TODO: érdemes lenne megnézni, hogy milyen értékek fordulnak elõ a gyakorlatban, kell-e korlátozni?
/*
        if (abs(deriv)>1.0e+007) {
            printf("uj_CBm=%g\n", uj_CBm);
            printf("prev_CBm.ertek=%g\n", prev_CBm.ertek);
            printf("TB=%g\n", TB);
            printf("TA=%g\n", TA);
            printf("nev=%g\n", nev);
            printf("deriv=%g\n", deriv);
            getchar();
        }
*/
        return ertek_t(uj_CBm, deriv);
    }

    //***********************************************************************
    void beolvas_fajlbol();
    //***********************************************************************
};


//***********************************************************************
struct LED_model_result_pack {
//***********************************************************************
    double G, Pdiss, Fi_e, Fi_nu;
    double dI_per_dVj; // differenciális vezetés
    double dI_per_dT;
    double dF_per_dVj; // F=Fi_e
    double dF_per_dT;  // F=Fi_e
    double IF, VF, I_rad, VF_rad, V_rad;
    bool votma_x;
    void set_votma(bool uj) { votma_x = uj; /*printf("set(%u)", (uns)uj);*/ }
    bool get_votma()const { return votma_x; }
    LED_model_result_pack() :G{ 1 }, Pdiss{ 0 }, Fi_e{ 1 }, Fi_nu{ 1 }, dI_per_dVj{ 0 }, dI_per_dT{ 0 },
        dF_per_dVj{ 0 }, dF_per_dT{ 0 }, IF{ 1 }, VF{ 1 }, I_rad{ 0.1 }, VF_rad{ 1 }, V_rad{ 1 } { set_votma(false); }
};


//***********************************************************************
class multi_domain_LED_model {
//***********************************************************************
    // type-1 paraméterek
    double Tref, VT;
    double I0e, me, Rse;
    double I0rad, mrad, Rrrad;
    double ael, bel, cel, del, eel, fel;
    double arad, brad, crad, drad, erad, frad;
    double akap, bkap, ckap, dkap, ekap, fkap, gkap, hkap, ikap;
    //***********************************************************************
    double calc_VT(double T) const { return 1.380649e-23 * (T + 273.15) / 1.60217662e-19; }
    //***********************************************************************
    double calc_VF_type_1(double IF, double Tj, double & VF_rad) const {
    //***********************************************************************
        if (IF <= 0)
            throw hiba("multi_domain_LED_model::calc_VF_type_1", "IF <= 0 (%g)", IF);
        IF = abs(IF);
        double V_dioda = me * VT * log(IF / I0e + 1);
        double dVF_el = (ael * IF * IF + bel * IF + cel) * (Tj * Tj - Tref * Tref) +
            (del * IF * IF + eel * IF + fel) * (Tj - Tref);
        double V_R = Rse * IF;
        VF_rad = V_dioda + dVF_el;
        return dVF_el + V_dioda + V_R;
    }
    //***********************************************************************
    double calc_VF_type_1_print(double IF, double Tj, double & VF_rad) const {
    //***********************************************************************
        if (IF <= 0)
            throw hiba("multi_domain_LED_model::calc_VF_type_1", "IF <= 0 (%g)", IF);
        IF = abs(IF);
        double V_dioda = me * VT * log(IF / I0e + 1);
        double dVF_el = (ael * IF * IF + bel * IF + cel) * (Tj * Tj - Tref * Tref) +
            (del * IF * IF + eel * IF + fel) * (Tj - Tref);
        double V_R = Rse * IF;
        VF_rad = V_dioda + dVF_el;
        static int db = 1000;
        db--;
        if (db == 0) {
            printf("VR=%g, Vd=%g, dV=%g, IF=%g\n", V_R, V_dioda, dVF_el, IF);
            printf("me=%g, VT=%g, I0e=%g, IF=%g\n", me, VT, I0e, IF);
            db = 1000;
        }
        return dVF_el + V_dioda + V_R;
    }
    //***********************************************************************
    double calc_VF_rad_type_1(double IF_rad, double Tj, double & V_rad) const {
    //***********************************************************************
        if (IF_rad <= 0)
            throw hiba("multi_domain_LED_model::calc_VF_rad_type_1", "IF_rad <= 0 (%g)", IF_rad);
        double Vrad_dioda = mrad * VT * log(IF_rad / I0rad + 1);
        double dVF_rad = (arad * IF_rad * IF_rad + brad * IF_rad + crad) * (Tj * Tj - Tref * Tref) +
            (drad * IF_rad * IF_rad + erad * IF_rad + frad) * (Tj - Tref);
        double V_R = Rrrad * IF_rad;
        V_rad = Vrad_dioda + dVF_rad;
        return dVF_rad + Vrad_dioda + V_R;
    }
    //***********************************************************************
    double iterate_IF_type_1(double IF_prev, double VF, double Tj) const {
    //***********************************************************************
        if (IF_prev <= 0)
            throw hiba("multi_domain_LED_model::iterate_IF_type_1", "IF_prev <= 0 (%g)", IF_prev);
        double V_dioda = me * VT * log(IF_prev / I0e + 1);
        double dVF_el = (ael * IF_prev * IF_prev + bel * IF_prev + cel) * (Tj * Tj - Tref * Tref) +
            (del * IF_prev * IF_prev + eel * IF_prev + fel) * (Tj - Tref);
        double IF = (VF - V_dioda - dVF_el) / Rse;
        // printf("\tIF = %.15g -> %.15g (%f%%)\n", IF_prev, IF, 100 * (IF - IF_prev) / IF_prev);
        return IF;
    }
    //***********************************************************************
    double iterate_I_rad_type_1(double IF_rad_prev, double VF_rad, double Tj) const {
    //***********************************************************************
        if (IF_rad_prev <= 0)
            throw hiba("multi_domain_LED_model::iterate_I_rad_type_1", "IF_rad_prev <= 0 (%g)", IF_rad_prev);
        double Vrad_dioda = mrad * VT * log(IF_rad_prev / I0rad + 1);
        double dVF_rad = (arad * IF_rad_prev * IF_rad_prev + brad * IF_rad_prev + crad) * (Tj * Tj - Tref * Tref) +
            (drad * IF_rad_prev * IF_rad_prev + erad * IF_rad_prev + frad) * (Tj - Tref);
        return (VF_rad - Vrad_dioda - dVF_rad) / Rrrad;
    }
    //***********************************************************************
    double Calc_IF_felezve(double I_start, double V_cel, double Tj, double & VF_rad) const{
    //***********************************************************************
        double VF_rad_calc;
        double y1 = calc_VF_type_1(I_start, Tj, VF_rad_calc) - V_cel;
        double I1 = I_start;
        if (abs(y1) > 1e-9) {
            double I2 = I_start * 0.5;
            double y2 = calc_VF_type_1(I2, Tj, VF_rad_calc) - V_cel;
            if (y1*y2 > 0) {
                double szorzo = abs(y2) > abs(y1) ? 2 : 0.5;
                I2 = I_start;
                uns it_num = 0;
                do {
                    I2 *= szorzo;
                    y2 = calc_VF_type_1(I2, Tj, VF_rad_calc) - V_cel;
                    if (y1*y2 < 0)
                        break;
                    if (abs(y2) < abs(y1))
                        I1 = I2;
                    else throw hiba("Calc_IF_felezve", "divergent");
                    it_num++;
                    if (it_num > 1000) {
                        throw hiba("Calc_IF_felezve", "iter 1: I1 = %g, I2 = %g, y1 = %g, y2 = %g\n", I1, I2, y1, y2);
                    }
                } while (true);
            }
            // A gyök I1 és I2 között van
            uns it_num = 0;
            double I3, y3;
            do {
                I3 = 0.5*(I1 + I2);
                y3 = calc_VF_type_1(I3, Tj, VF_rad_calc) - V_cel;
                if (y1*y3 > 0) {
                    I1 = I3;
                    y1 = y3;
                }
                else {
                    I2 = I3;
                    y2 = y3;
                }
                it_num++;
                if (it_num > 1000) {
                    throw hiba("Calc_IF_felezve", "iter 2: I1 = %g, I2 = %g, I3= %g, y1 = %g, y2 = %g, y3=%g\n", I1, I2, I3, y1, y2, y3);
                }
            } while (abs(y3) > 1e-9);
            I1 = I3;
        }
        return I1;
    }
    //***********************************************************************
    double calc_IF_type_1(double IF_prev, double VF, double Tj, double & VF_rad) const {
    //***********************************************************************
        double VF_rad_calc;
        double VF_calc = calc_VF_type_1(IF_prev, Tj, VF_rad_calc);
        double I_start = IF_prev;
        uns it_num = 0;
        while (abs(VF_calc - VF) > 1e-9) {
            IF_prev = iterate_IF_type_1(IF_prev, VF, Tj);
            if (IF_prev < 0 || it_num > 500) {
                IF_prev = Calc_IF_felezve(I_start, VF, Tj, VF_rad_calc);
            }
            VF_calc = calc_VF_type_1(IF_prev, Tj, VF_rad_calc);
            it_num++;
            if (it_num > 1000) {
                throw hiba("calc_IF_type_1", "VF_calc = %g, VF = %g, IF_prev = %g\n", VF_calc, VF, IF_prev);
            }
        }
        //calc_VF_type_1_print(IF_prev, Tj, VF_rad_calc); // csak a kiírás miatt
        VF_rad = VF_rad_calc;
        return IF_prev;
    }
    //***********************************************************************
    double Calc_I_rad_felezve(double I_start, double V_cel, double Tj, double & V_rad) const{
    //***********************************************************************
        double V_rad_calc;
        double y1 = calc_VF_rad_type_1(I_start, Tj, V_rad_calc) - V_cel;
        double I1 = I_start;
        if (abs(y1) > 1e-9) {
            double I2 = I_start * 0.5;
            double y2 = calc_VF_rad_type_1(I2, Tj, V_rad_calc) - V_cel;
            if (y1*y2 > 0) {
                double szorzo = abs(y2) > abs(y1) ? 2 : 0.5;
                I2 = I_start;
                uns it_num = 0;
                do {
                    I2 *= szorzo;
                    y2 = calc_VF_rad_type_1(I2, Tj, V_rad_calc) - V_cel;
                    if (y1*y2 < 0)
                        break;
                    if (abs(y2) < abs(y1))
                        I1 = I2;
                    else throw hiba("Calc_I_rad_felezve", "divergent");
                    it_num++;
                    if (it_num > 100000) {
                        throw hiba("Calc_I_rad_felezve", "iter 1: I1 = %g, I2 = %g, y1 = %g, y2 = %g\n", I1, I2, y1, y2);
                    }
                } while (true);
            }
            // A gyök I1 és I2 között van
            uns it_num = 0;
            double I3, y3;
            do {
                I3 = 0.5*(I1 + I2);
                y3 = calc_VF_rad_type_1(I3, Tj, V_rad_calc) - V_cel;
                if (y1*y3 > 0) {
                    I1 = I3;
                    y1 = y3;
                }
                else {
                    I2 = I3;
                    y2 = y3;
                }
                it_num++;
                if (it_num > 100000) {
                    throw hiba("Calc_I_rad_felezve", "iter 2: I1 = %g, I2 = %g, I3= %g, y1 = %g, y2 = %g, y3=%g\n", I1, I2, I3, y1, y2, y3);
                }
            } while (abs(y3) > 1e-9);
            I1 = I3;
        }
        return I1;
    }
    //***********************************************************************
    double calc_I_rad_type_1(double IF_rad_prev, double VF_rad, double Tj, double & V_rad) const {
    //***********************************************************************
        double V_rad_calc;
        double VF_rad_calc = calc_VF_rad_type_1(IF_rad_prev, Tj, V_rad_calc);
        double I_start = IF_rad_prev;
        uns it_num = 0;
        while (abs(VF_rad_calc - VF_rad) > 1e-9) {
            double uj_IF_rad_prev = iterate_I_rad_type_1(IF_rad_prev, VF_rad, Tj);
            if (uj_IF_rad_prev < 0 || abs(uj_IF_rad_prev / IF_rad_prev) > 10 || abs(IF_rad_prev / uj_IF_rad_prev) > 10) {
                IF_rad_prev = Calc_I_rad_felezve(I_start, VF_rad, Tj, V_rad_calc);
            }
            else IF_rad_prev = uj_IF_rad_prev;
            VF_rad_calc = calc_VF_rad_type_1(IF_rad_prev, Tj, V_rad_calc);
            it_num++;
            if (it_num > 10000000) {
                throw hiba("calc_I_rad_type_1", "VF_rad_calc = %g, VF_rad = %g, IF_rad_prev = %g\n", VF_rad_calc, VF_rad, IF_rad_prev);
            }
        }
        V_rad = V_rad_calc;
        return IF_rad_prev;
    }
    //***********************************************************************
    double calc_Fi_e(double VF, double IF, double I_rad, double Tj) const {
    //***********************************************************************
        double VF_rad, V_rad;
        calc_IF_type_1(IF, VF, Tj, VF_rad);
        I_rad = calc_I_rad_type_1(I_rad, VF_rad, Tj, V_rad);
        return I_rad*V_rad;
    }
    //***********************************************************************
    double calc_dFi_e_per_dVF(double VF, double IF, double I_rad, double Tj) const {
    //***********************************************************************
        return (calc_Fi_e(VF + 1e-6, IF, I_rad, Tj) - calc_Fi_e(VF - 1e-6, IF, I_rad, Tj)) / 2e-6;
    }
    //***********************************************************************
    double calc_dFi_e_per_dT(double VF, double IF, double I_rad, double Tj) const {
    //***********************************************************************
        return (calc_Fi_e(VF, IF, I_rad, Tj + +1e-6) - calc_Fi_e(VF, IF, I_rad, Tj - 1e-6)) / 2e-6;
    }

public:
    //***********************************************************************
    multi_domain_LED_model() {
    //***********************************************************************
        Tref = 70;
        VT = calc_VT(Tref);
        I0e = 7.32187956118988e-19;
        me = 2.29275643623816;//2.275;//
        Rse = 1.532;//1.25;//
        I0rad = 4.29260164799285e-23;
        mrad = 1.83561771204239;
        Rrrad = 0.899;
        ael = -4.20320084575734e-05;
        bel = 4.69533080515079e-05;
        cel = 8.3781624233249e-07;
        del = 0.0100986814367841;
        eel = -0.0125056006514783;
        fel = -0.00107358797611507;
        arad = 0.000038198290216593;
        brad = 0.0000225698585460935;
        crad = 2.86052740648461e-06;
        drad = -0.0122727605875887;
        erad = -0.00747414076628783;
        frad = -0.00132565719838538;
        akap = 0;
        bkap = 0;
        ckap = 0;
        dkap = 0;
        ekap = 0;
        fkap = 0;
        gkap = 0;
        hkap = 0;
        ikap = 1000;
    }
/*
    //***********************************************************************
    multi_domain_LED_model() {
    //***********************************************************************
        Tref = 70;
        VT = calc_VT(Tref);
        I0e = 9.4235e-15;
        me = 2.993840735;
        Rse = 0.456051822;
        I0rad = 1.95675e-14;
        mrad = 3.132431602;
        Rrrad = 0.051179859;
        ael = -1.6221E-05;
        bel = 2.49301E-05;
        cel = 6.03936E-06;
        del = 0.003537067;
        eel = -0.006962476;
        fel = -0.002481474;
        arad = -1.28415E-05;
        brad = 2.24149E-05;
        crad = 7.93926E-06;
        drad = 0.003486183;
        erad = -0.006874837;
        frad = -0.002486097;
        akap = -0.00189312;
        bkap = 0.159623633;
        ckap = 2.689638531;
        dkap = 0.004088094;
        ekap = -0.439392729;
        fkap = 1.008384099;
        gkap = -0.001998127;
        hkap = 0.299603441;
        ikap = 35.0516427;
    }
*/
    //***********************************************************************
    // A junction feszültsége alapján számolja ki a visszaadott értékeket
    void calc_el_from_U(double VF, double Tj, double start_I, LED_model_result_pack & ret) const {
    // A junction fesz alapján számolja ki a visszaadott értékeket egy kezdõ áramból iterálva
    //***********************************************************************
        VF = abs(VF);
        if (VF == 0) {
            VF = calc_VF_type_1(0.1, Tj, ret.VF_rad);
        }
        if (VF < 1e-9)VF = 1e-9; // 1 nV alatt nem számol
        if (VF > 4)VF = 2.9;
        ret.VF = VF;
        ret.IF = calc_IF_type_1(start_I, VF, Tj, ret.VF_rad);
        //printf("From U: U=%g, start I = %g, uj IF = %g, Tj = %g\n", VF, start_I, ret.IF, Tj);
        ret.G = ret.IF / ret.VF;
        //ret.dI_per_dVj = 1 / (Rse + me*VT / (ret.IF + I0e)
        //    + (2 * ael * ret.IF + bel) * (Tj * Tj - Tref * Tref) + (2 * del * ret.IF + eel) * (Tj - Tref));
        ret.dI_per_dVj = ret.G;
        double dV_per_dT = 2 * (ael * ret.IF * ret.IF + bel * ret.IF + cel) * Tj + (del * ret.IF * ret.IF + eel * ret.IF + fel);
        ret.dI_per_dT = ret.dI_per_dVj * dV_per_dT;
        chech_fpu_error_name(ret.VF, "ret.VF");
        chech_fpu_error_name(ret.IF, "ret.IF");
        chech_fpu_error_name(ret.G, "ret.G");
        chech_fpu_error_name(ret.dI_per_dVj, "ret.dI_per_dVj");
        chech_fpu_error_name(ret.dI_per_dT, "ret.dI_per_dT");
    }
    //***********************************************************************
    void calc_el_from_I(double IF, double Tj, LED_model_result_pack & ret) const {
    // A junction árama alapján számolja ki a visszaadott értékeket
    //***********************************************************************
        IF = abs(IF);
        if (IF == 0)IF = 1;
        if (IF < 1e-9)IF = 1e-9; // 1 nA alatt nem számol
        ret.IF = IF;
        ret.VF = calc_VF_type_1(IF, Tj, ret.VF_rad);
        //printf("From I: I=%g, uj VF=%g, Tj = %g\n", IF, ret.VF, Tj);
        ret.G = ret.IF / ret.VF;
        //ret.dI_per_dVj = 1 / (Rse + me*VT/(ret.IF+I0e) 
        //    + (2 * ael * ret.IF + bel) * (Tj * Tj - Tref * Tref) + (2 * del * ret.IF + eel) * (Tj - Tref));
        ret.dI_per_dVj = ret.G;
        double dV_per_dT = 2 * (ael * ret.IF * ret.IF + bel * ret.IF + cel) * Tj + (del * ret.IF * ret.IF + eel * ret.IF + fel);
        ret.dI_per_dT = ret.dI_per_dVj * dV_per_dT;
        chech_fpu_error_name(ret.VF, "ret.VF");
        chech_fpu_error_name(ret.IF, "ret.IF");
        chech_fpu_error_name(ret.G, "ret.G");
        chech_fpu_error_name(ret.dI_per_dVj, "ret.dI_per_dVj");
        chech_fpu_error_name(ret.dI_per_dT, "ret.dI_per_dT");
    }
    //***********************************************************************
    void calc_rad(double Tj, LED_model_result_pack & ret) const {
    //***********************************************************************
        double Ibe = ret.I_rad;
        //ret.I_rad = 0.053102915895065522;
        //ret.VF_rad = 2.7611807211808559;
        //Tj = 61.50496099839561026101;
        ret.I_rad = calc_I_rad_type_1(ret.I_rad, ret.VF_rad, Tj, ret.V_rad);
        ret.Fi_e = ret.I_rad*ret.V_rad;
        if (fpu_error(ret.Fi_e)) {
            printf("\n%.20f\n", Tj);
        }
        ret.Pdiss = ret.IF*ret.VF - ret.Fi_e;
        ret.dF_per_dVj = calc_dFi_e_per_dVF(ret.VF, ret.IF, ret.I_rad, Tj);
        ret.dF_per_dT = 0; //calc_dFi_e_per_dT(ret.VF, ret.IF, ret.I_rad, Tj);
        chech_fpu_error_name(ret.I_rad, "ret.I_rad");
        chech_fpu_error_name(ret.Fi_e, "ret.Fi_e");
        chech_fpu_error_name(ret.Pdiss, "ret.Pdiss");
        chech_fpu_error_name(ret.dF_per_dVj, "ret.dF_per_dVj");
        chech_fpu_error_name(ret.dF_per_dT, "ret.dF_per_dT");
    }
    //***********************************************************************
    void calc_lum(double Tj, LED_model_result_pack & ret) const {
    //***********************************************************************
        double K = (akap*Tj*Tj + bkap*Tj + ckap)*ret.IF*ret.IF 
                 + (dkap*Tj*Tj + ekap*Tj + fkap)*ret.IF 
                 + (gkap*Tj*Tj + hkap*Tj + ikap);
        ret.Fi_nu = ret.Fi_e*K;
    }
};


//***********************************************************************
class tulajdonsag {
//***********************************************************************
    //zárttá kell alakítani, legyen is paramétere, set_konstans(...), set_lin, ....
    bool is;
    nonlin_tipus tipus, inv_tipus;  // nt_konstans, nt_lin, nt_szakszok, nt_fazisvalto, nt_exp, nt_mizs, nt_diode_1, nt_erno, nt_polinom, nt_inv
    rvt ertek;                      // alapérték, tolt: Tamb-ra normalizált
    rvt seged;                      // exp kitevõ konstansa, lin meredeksége, tolt: Tamb-ra normalizált
    vektor<rvt> par;                // típusfüggõ segédparaméterek, tolt: Tamb-ra normalizált
    broken_line szakasz;            // nt_szakszok esetén a törtvonal
    fazis_broken_line fazis;        // nt_fazisvalto esetén a dupla törtvonal

    rvt ertek_tolt;                 // alapérték, tolt: Tamb-ra normalizált
    rvt seged_tolt;                 // exp kitevõ konstansa, lin meredeksége, tolt: Tamb-ra normalizált
    vektor<rvt> par_tolt;           // típusfüggõ segédparaméterek, tolt: Tamb-ra normalizált
    // polinom ill. illesztett polinom paraméterei
    bool is_illesztett;             // ha false, akkor eleve polinomként adott, ha true, akkor mérési pontokként adott, amibõl beolvasás után illesztett polinomot
    illesztendo_adatok ill_adat;    // illesztett esetben a mérési pontok (U-ra és T-re illeszti I-t), ha esetleg kellene, megõrzi
    polinom poli;                   // ha polinommal adott a tulajdonság
    // junction elõszámításhoz
    bool is_junction_pre_calc;      // elõszámító lépés junction esetén, konstans vezetést/fényt ad vissza
    //***********************************************************************
    LED_model_result_pack junction_pre_pack; // ha elõszámító lépés fut, ezt adja vissza a calc
    multi_domain_LED_model jani_junction_model;

    //***********************************************************************
    double mizs(double T_tolt) const{
    // a driválás miatt nem rvt, hanem double
    //***********************************************************************
        double Tabs = T_tolt + seged_tolt; // seged_tolt = Tamb + 273.15
        double Tbc = sqr(Tabs - par[1]);   // par[2] lenne a hatvány, de fixáljuk c=2-re
        double signa = Tabs > par[1] ? -par[0] : par[0];
        double tag1 = signa * Tbc / (Tbc + par[3]) + par[0] + par[5];
        return 0.01*(tag1*exp(par[4] / Tabs) + par[6]); // 0.01* az ohm cm miatt
    }
    //***********************************************************************
    ertek_t get_value_konstans(rvt T_tolt, const ertek_t & H) const { return ertek_t(ertek_tolt, rvt()); }
    //***********************************************************************
    ertek_t get_value_lin(rvt T_tolt, const ertek_t & H) const { return ertek_t(ertek_tolt + seged_tolt * T_tolt, seged_tolt); }
    //***********************************************************************
    ertek_t get_value_exp(rvt T_tolt, const ertek_t & H) const { return exp_par(T_tolt, ertek_tolt, seged_tolt); }
    //***********************************************************************
    ertek_t get_value_szakaszok(rvt T_tolt, const ertek_t & H) const { return szakasz.value(T_tolt); }
    //***********************************************************************
    ertek_t get_value_fazisvalto(rvt T_tolt, const ertek_t & H) const { return fazis.value(T_tolt, H); }
    //***********************************************************************

    //***********************************************************************
    ertek_t get_value_mizs(rvt T_tolt, const ertek_t & H) const {
    //***********************************************************************
        double f = mizs(T_tolt); // a deriválás miatt használunk double-t
        return ertek_t(rvt(f), rvt((mizs(T_tolt + 1.0e-6) - f) * 1.0e+6));
    }

    //***********************************************************************
    ertek_t get_value_erno(rvt T_tolt, const ertek_t & H) const {
    //***********************************************************************
        throw hiba("get_value_erno", "erno is not resistance");
    }

    //***********************************************************************
    ertek_t get_value_polinom(rvt T_tolt, const ertek_t & H) const {
    //***********************************************************************
        throw hiba("get_value_polinom", "polinom is not resistance");
    }

    //***********************************************************************
    ertek_t (tulajdonsag::*fvp_get_value)(rvt T_tolt, const ertek_t & H) const;
    ertek_t (tulajdonsag::*fvp_get_inv_value)(rvt T_tolt, const ertek_t & H) const;
    //***********************************************************************

    //***********************************************************************
    ertek_t get_value_inv(rvt T_tolt, const ertek_t & H) const {
    //***********************************************************************
        ertek_t t = (this->*fvp_get_inv_value)(T_tolt, H);
        rvt f = rvt(1.0) / t.ertek;
        return ertek_t(f, -t.derivalt * f * f);
    }

    //***********************************************************************
    ertek_t get_G_junction(rvt T_tolt, rvt V_junction, rvt & diff_vezetes, bool is_print=false) const{
    // A pn átmenet teljes felületének vezetését adja, a face vezetése ez szorozva
    // a face felület és a teljes junction felület arányával.
    // A diff_vezetes a teljes pn átmenet dI/dU differenciális vezetése. (Nem a G deriváltja !!)
    // A derivált a dI/dT
    //***********************************************************************
        if (is_junction_pre_calc) {
            diff_vezetes = junction_pre_pack.dI_per_dVj;
            return ertek_t(junction_pre_pack.G, junction_pre_pack.dI_per_dT);
        }
        ertek_t ret;
        if (tipus == nt_erno) {
            rvt T = T_tolt + seged_tolt; // seged_tolt-ba a T_amb került
            if (T > 130)
                T = 130;
            if (T < 20)
                T = 20;
            rvt b = par.unsafe(0) + par.unsafe(1)*T + par.unsafe(2)*T*T;
            rvt m = par.unsafe(3) + par.unsafe(4)*T + par.unsafe(5)*T*T;
            if (V_junction == 0)
                V_junction = pow(0.1 / b, 1 / m);
            V_junction = abs(V_junction);
            //if (V_junction < 0.1)
            //    V_junction = 0.1;
            b = abs(b);
            rvt U_az_m_miusz_egyediken = pow(V_junction, m - 1);
            ret.ertek = b * U_az_m_miusz_egyediken; // I = b*U^m => G = I/U = b*u^(m-1)
            diff_vezetes = ret.ertek*m; // dI/dU = m*b*U^(m-1)
            ret.derivalt = U_az_m_miusz_egyediken*((par.unsafe(1) + 2 * par.unsafe(2)*T)*V_junction + b*m*(par.unsafe(4) + 2 * par.unsafe(5)*T));
            if (is_print)
                printf("T=%g, V=%g, b=%g, m=%g, I=%g\n", T, V_junction, b, m, ret.ertek*V_junction);
        }
        else if (tipus == nt_polinom) {
            rvt T = T_tolt + seged_tolt; // seged_tolt-ba a T_amb került
            if (V_junction == 0)
                V_junction = 1e-3;//0.5*(poli.min_U + poli.max_U);
/*            V_junction = abs(V_junction);
            if (V_junction > poli.max_U * 1.3)
                V_junction = 1.3 * poli.max_U;
            if (T > 150)
                T = 150;
*/            Itrio trio;
            trio = poli.get_Itrio(V_junction, T, true);
            ret.ertek = abs(poli.szorzo * trio.I / V_junction);//0.063;0.221; //
            ret.derivalt = poli.szorzo * trio.dI_per_dT;//0;0;//
            diff_vezetes = poli.szorzo * trio.dI_per_dU;//0.063;0.221;//
            //printf("T_tolt=%g, seged_tolt=%g\n", T_tolt, seged_tolt);
            //printf("T=%g, V=%g, I=%g, G=%g, dI/dT=%g, dI/dU=%g\n", T, V_junction, poli.szorzo * trio.I, ret.ertek, ret.derivalt, diff_vezetes);
        }
        else {
            throw hiba("get_G_junction", "unsupported junction type");
        }
        return ret;
    }

    //***********************************************************************
    ertek_t get_rad_flux(rvt T_tolt, rvt I_junction, rvt g_diff_j, rvt & dF_per_dVj) const{
    // Bemenet a face árama felszorozva, mintha az egész junction árama volna, 
    // mert a karakterisztika az egész junctionre vonatkozik.
    // Azaz I_junction = Iface * Ajunction / Aface.
    // Bemenet még a face hõmérséklete, valamint a differenciális vezetése.
    // Ez a fluxus feszültségszerinti deriváltjához kell.
    // Vissza: a fluxus és hõmérséklet szerinti deriváltja returnnel, valamint
    // a fluxus feszültség szerinti deriváltja, paraméterként.
    //***********************************************************************
        if (is_junction_pre_calc) {
            dF_per_dVj = junction_pre_pack.dF_per_dVj;
            return ertek_t(junction_pre_pack.Fi_e, junction_pre_pack.dF_per_dT);
        }
        ertek_t ret;
        if (tipus == nt_erno) {
            rvt T = T_tolt + seged_tolt; // seged_tolt-ba a T_amb került
            if (T > 130)
                T = 130;
            if (T < 20)
                T = 20;
            rvt b = par.unsafe(0) + par.unsafe(1)*T + par.unsafe(2)*T*T;
            rvt m = par.unsafe(3) + par.unsafe(4)*T + par.unsafe(5)*T*T;
            rvt db_per_dT = par.unsafe(1) + 2 * par.unsafe(2) * T;
            rvt dm_per_dT = par.unsafe(4) + 2 * par.unsafe(5) * T;
            I_junction = abs(I_junction);
            ret.ertek = m * I_junction + b; // F = m(T) * I(U) + b(T)
            if (ret.ertek < 0) {
                ret.ertek = 0;
                dF_per_dVj = 0;
                ret.derivalt = 0;
            }
            else {
                dF_per_dVj = g_diff_j * m; // dF/dVj = m * dI/dU = m * g_diff
                ret.derivalt = dm_per_dT * I_junction + db_per_dT;
            }
        }
        else if (tipus == nt_polinom) {
            rvt T = T_tolt + seged_tolt; // seged_tolt-ba a T_amb került
            I_junction = abs(I_junction);
            if (T > 150)
                T = 150;
            Itrio trio;
            trio = poli.get_Itrio(I_junction, T, true);
            ret.ertek = poli.szorzo * trio.I;
            ret.derivalt = poli.szorzo * trio.dI_per_dT;
            dF_per_dVj = poli.szorzo * trio.dI_per_dU * g_diff_j; // trio.dI_per_dU valójában dF/dI, így dF/dV = dF/dI*dI/dV
        }
        else if (tipus == nt_konstans) {
            ret.ertek = ertek_tolt*I_junction; // C*I
            ret.derivalt = rvt();
            dF_per_dVj = ertek_tolt * g_diff_j; // C*dI/dU
        }
        else {
            throw hiba("get_rad_flux", "unsupported junction type");
        }
        return ret;
    }
public:
    //***********************************************************************
    rvt get_lum_flux(rvt T_tolt, rvt I_junction) const{
    // Bemenet a face árama felszorozva, mintha az egész junction árama volna, 
    // mert a karakterisztika az egész junctionre vonatkozik.
    // Azaz I_junction = Iface * Ajunction / Aface.
    // Bemenet még a face hõmérséklete, valamint a differenciális vezetése.
    // Ez a fluxus feszültségszerinti deriváltjához kell.
    // Vissza: a fluxus és hõmérséklet szerinti deriváltja returnnel, valamint
    // a fluxus feszültség szerinti deriváltja, paraméterként.
    //***********************************************************************
        if (is_junction_pre_calc) {
            return junction_pre_pack.Fi_nu;
        }
        rvt ret;
        if (tipus == nt_erno) {
            rvt T = T_tolt + seged_tolt; // seged_tolt-ba a T_amb került

            if (T > 130)
                T = 130;
            if (T < 20)
                T = 20;
            rvt b = par.unsafe(0) + par.unsafe(1)*T + par.unsafe(2)*T*T;
            rvt m = par.unsafe(3) + par.unsafe(4)*T + par.unsafe(5)*T*T;
            I_junction = abs(I_junction);

            ret = m * I_junction + b; // F = m(T) * I(U) + b(T)
            if (ret < 0) {
                ret = 0;
            }
            //printf("%g, %g, %g, %g, %g, %g\n", par.unsafe(0), par.unsafe(1), par.unsafe(2), par.unsafe(3), par.unsafe(4), par.unsafe(5));
            ///*if(I_junction / ret>50)*/printf("T=%g, I=%g, m=%g, b=%g, ret=%g, I/ret=%g\n", T, I_junction, m, b, ret, I_junction/ret);
        }
        else if (tipus == nt_polinom) {
            rvt T = T_tolt + seged_tolt; // seged_tolt-ba a T_amb került
            I_junction = abs(I_junction);
            if (T > 150)
                T = 150;
            Itrio trio;
            trio = poli.get_Itrio(I_junction, T, true);
            ret = poli.szorzo * trio.I;
        }
        else if (tipus == nt_konstans) {
            ret = ertek_tolt*I_junction; // C*I
        }
        else {
            throw hiba("get_lum_flux", "unsupported junction type");
        }
        return ret;
    }

public:
    //***********************************************************************
    tulajdonsag(rvt kezdoertek = rvt()) :is{ false }, tipus { nt_konstans }, inv_tipus{ nt_konstans },
        ertek{ kezdoertek }, ertek_tolt{ kezdoertek }, seged{ rvt() }, seged_tolt{ rvt() },
        is_junction_pre_calc{ false } {}
    //***********************************************************************
    bool is_exists()const { return is; }
    //***********************************************************************
    bool is_const()const { return tipus == nt_konstans; }
    //***********************************************************************
    void set_to_konstans(rvt Ertek) { is = true; tipus = nt_konstans; ertek = Ertek; seged = rvt(); fvp_get_value = &tulajdonsag::get_value_konstans; }
    //***********************************************************************
    void set_to_lin(rvt b, rvt m) { is = true; tipus = nt_lin; ertek = b; seged = m; fvp_get_value = &tulajdonsag::get_value_lin; }
    //***********************************************************************
    void set_to_exp(rvt A, rvt gamma) { is = true; tipus = nt_exp; ertek = A; seged = gamma; fvp_get_value = &tulajdonsag::get_value_exp; }
    //***********************************************************************
    broken_line & set_to_szakaszok() { is = true; tipus = nt_szakszok; fvp_get_value = &tulajdonsag::get_value_szakaszok; return szakasz; }
    //***********************************************************************
    fazis_broken_line & set_to_fazisvalto() { is = true; tipus = nt_fazisvalto; fvp_get_value = &tulajdonsag::get_value_fazisvalto; return fazis; }
    //***********************************************************************
    void set_to_erno(rvt a, rvt b, rvt c, rvt d, rvt e, rvt f) { 
    //***********************************************************************
        is = true; tipus = nt_erno; 
        par.set_size(6);
        par.unsafe(0) = a; par.unsafe(1) = b; par.unsafe(2) = c;
        par.unsafe(3) = d; par.unsafe(4) = e; par.unsafe(5) = f;
        fvp_get_value = &tulajdonsag::get_value_erno;
    }
    //***********************************************************************
    void set_to_polinom() { 
    //***********************************************************************
        is = true; tipus = nt_jani_1;// nt_polinom; 
        fvp_get_value = &tulajdonsag::get_value_polinom;
    }
    //***********************************************************************
    void set_to_diode_1(rvt a, rvt b, rvt c) { is = true; tipus = nt_diode_1; TODO("set_to_diode_1"); }
    //***********************************************************************

    //***********************************************************************
    void set_to_mizs(rvt a, rvt b, rvt d, rvt e, rvt f, rvt g) {
    //***********************************************************************
        is = true; tipus = nt_mizs;
        par.set_size(7);
        par.unsafe(0) = a; par.unsafe(1) = b; par.unsafe(2) = 2; par.unsafe(3) = d; 
        par.unsafe(4) = e; par.unsafe(5) = f; par.unsafe(6) = g;
        fvp_get_value = &tulajdonsag::get_value_mizs;
    }
    //***********************************************************************

    //***********************************************************************
    void set_to_inv() {  // a set_to_xxxx UTÁN hívandó
    //***********************************************************************
        inv_tipus = tipus; 
        tipus = nt_inv; 
        fvp_get_inv_value = fvp_get_value;
        fvp_get_value = &tulajdonsag::get_value_inv;
    }

    //***********************************************************************
    ertek_t get_value(rvt T_tolt, const ertek_t & H) const{
    //***********************************************************************
        return (this->*fvp_get_value)(T_tolt, H);
    }

    //***********************************************************************
    // A junction feszültsége alapján számolja ki a visszaadott értékeket
    void junction_calc_el_from_U(double VF, double Tj_tolt, double start_I, LED_model_result_pack & ret, bool from_I_if_possible) const {
    //***********************************************************************
        if (is_junction_pre_calc) {
            bool is_votma = ret.get_votma();
            ret = junction_pre_pack;
            if (is_votma && !ret.get_votma())
                ret.set_votma(true);
            return;
        }
        switch (tipus) {
            case nt_erno:
            case nt_polinom: {
                    rvt diff_vez;
                    ertek_t G = get_G_junction(Tj_tolt, VF, diff_vez);
                    ret.VF = VF;
                    ret.G = G.ertek;
                    ret.IF = ret.G * ret.VF;
                    ret.dI_per_dVj = diff_vez;
                    ret.dI_per_dT = G.derivalt;
                }
            break;
            case nt_jani_1:
                if(from_I_if_possible)
                    jani_junction_model.calc_el_from_I(start_I, Tj_tolt + seged_tolt, ret);
                else
                    jani_junction_model.calc_el_from_U(VF, Tj_tolt + seged_tolt, start_I, ret);
                break;
            default:
                throw hiba("calc_el_from_U", "unsupported junction type");
        }
    }
    //***********************************************************************
    void junction_calc_el(double VF, double Tj_tolt, double start_I, LED_model_result_pack & ret) const {
    // Ha a modell támogatja, a junction árama alapján számolja ki a visszaadott értékeket, egyébként a fszültsége alapján
    //***********************************************************************
        if (is_junction_pre_calc) {
            bool is_votma = ret.get_votma();
            ret = junction_pre_pack;
            if (is_votma && !ret.get_votma())
                ret.set_votma(true);
            return;
        }
        switch (tipus) {
            case nt_erno:
            case nt_polinom: {
                    rvt diff_vez;
                    ertek_t G = get_G_junction(Tj_tolt, VF, diff_vez);
                    ret.VF = VF;
                    ret.G = G.ertek;
                    ret.IF = ret.G * ret.VF;
                    ret.dI_per_dVj = diff_vez;
                    ret.dI_per_dT = G.derivalt;
                }
            break;
            case nt_jani_1:
                if (ret.get_votma()) {
                    jani_junction_model.calc_el_from_I(start_I, Tj_tolt + seged_tolt, ret);
                }
                else {
                    jani_junction_model.calc_el_from_U(VF, Tj_tolt + seged_tolt, start_I, ret);
                    ret.set_votma(true);
                }
                break;
            default:
                throw hiba("calc_el_from_U", "unsupported junction type");
        }
    }
    //***********************************************************************
    void junction_calc_rad(double Tj_tolt, LED_model_result_pack & ret) const {
    //***********************************************************************
        if (is_junction_pre_calc) {
            bool is_votma = ret.get_votma();
            ret = junction_pre_pack;
            if (is_votma && !ret.get_votma())
                ret.set_votma(true);
            return;
        }
        switch (tipus) {
            case nt_erno:
            case nt_polinom: {
                    rvt dF_per_dVj;
                    ertek_t F = get_rad_flux(Tj_tolt, ret.IF, ret.dI_per_dVj, dF_per_dVj);
                    ret.I_rad = 0;
                    ret.Fi_e = F.ertek;
                    ret.Pdiss = ret.VF * ret.IF - ret.Fi_e;
                    ret.dF_per_dVj = dF_per_dVj;
                    ret.dF_per_dT = F.derivalt;
            }
            break;
            case nt_jani_1:
                jani_junction_model.calc_rad(Tj_tolt + seged_tolt, ret);
                break;
            default:
                throw hiba("calc_rad_lum", "unsupported junction type");
        }
    }
    //***********************************************************************
    void junction_calc_lum(double Tj_tolt, LED_model_result_pack & ret) const {
    //***********************************************************************
        if (is_junction_pre_calc) {
            bool is_votma = ret.get_votma();
            ret = junction_pre_pack;
            if (is_votma && !ret.get_votma())
                ret.set_votma(true);
            return;
        }
        switch (tipus) {
            case nt_erno:
            case nt_polinom: {
                    ret.Fi_nu = get_lum_flux(Tj_tolt, ret.IF);
            }
            break;
            case nt_jani_1:
                jani_junction_model.calc_lum(Tj_tolt + seged_tolt, ret);
                break;
            default:
                throw hiba("calc_rad_lum", "unsupported junction type");
        }
    }

    //***********************************************************************
    void set_junction_normal_calc() { is_junction_pre_calc = false; }
    LED_model_result_pack * set_junction_pre_calc(rvt I, rvt T_tolt, bool is_el, LED_model_result_pack * el_pack_rad_hoz) {
    // visszaadja az elektromos pakkot a rad kiszámításához. A rad nullptr-t ad.
    //***********************************************************************
        if (tipus == nt_polinom) {
            is_junction_pre_calc = false;
            if (is_el) {
                rvt U = poli.calc_U(I, T_tolt + seged_tolt);
                rvt junction_pre_diff_vez = I / U;
                junction_pre_pack.IF = I;
                junction_pre_pack.VF = U;
                junction_pre_pack.G = junction_pre_diff_vez;
                junction_pre_pack.dI_per_dVj = junction_pre_diff_vez;
                junction_pre_pack.dI_per_dT = 0;
                is_junction_pre_calc = true;
                return &junction_pre_pack;
            }
            else {
                rvt dummy;
                ertek_t junction_pre_rad_flux = get_rad_flux(T_tolt, I, 0, dummy);
                rvt junction_pre_lum_flux = get_lum_flux(T_tolt, I);
                
                if (el_pack_rad_hoz == nullptr)
                    throw hiba("set_junction_pre_calc", "el_pack_rad_hoz == nullptr");

                el_pack_rad_hoz->Fi_e = junction_pre_rad_flux.ertek;
                el_pack_rad_hoz->Fi_nu = junction_pre_lum_flux;
                el_pack_rad_hoz->Pdiss = el_pack_rad_hoz->IF * el_pack_rad_hoz->VF - el_pack_rad_hoz->Fi_e;
                el_pack_rad_hoz->I_rad = el_pack_rad_hoz->VF_rad = el_pack_rad_hoz->V_rad = 0;
                el_pack_rad_hoz->dF_per_dVj = el_pack_rad_hoz->dF_per_dT = 0;
                junction_pre_pack = *el_pack_rad_hoz;
                is_junction_pre_calc = true;
                return nullptr;
            }
        }
        else if (tipus == nt_jani_1) {
            is_junction_pre_calc = false;
            if (is_el) {
                jani_junction_model.calc_el_from_I(I/* / 12*/, T_tolt + seged_tolt, junction_pre_pack);
                junction_pre_pack.dI_per_dVj = junction_pre_pack.G;
                junction_pre_pack.dI_per_dT = 0;
                is_junction_pre_calc = true;
                return &junction_pre_pack;
            }
            else {
                jani_junction_model.calc_rad(T_tolt + seged_tolt, *el_pack_rad_hoz);
                jani_junction_model.calc_lum(T_tolt + seged_tolt, *el_pack_rad_hoz);
                junction_pre_pack = *el_pack_rad_hoz;
                is_junction_pre_calc = true;
                return nullptr;
            }
        }
        else {
            throw hiba("set_junction_pre_calc", "unsupported junction type");
        }
    }

    //***********************************************************************
    void set_Tamb(rvt Tamb) {
    //***********************************************************************
        nonlin_tipus tip = tipus == nt_inv ? inv_tipus : tipus;
        switch (tip) {
            case nt_konstans:
                ertek_tolt = ertek;
                seged_tolt = seged;
                break;
            case nt_lin: // return ertek + seged * (T - 25.0)
                ertek_tolt = ertek + seged * (Tamb - rvt(25.0));
                seged_tolt = seged;
                break;
            case nt_exp: // return ertek*exp(seged * (T - 25.0))
                ertek_tolt = ertek * exp(seged * (Tamb - rvt(25.0)));
                seged_tolt = seged;
                break;
            case nt_szakszok:
                szakasz.set_Tamb(Tamb);
                break;
            case nt_fazisvalto:
                fazis.set_Tamb(Tamb);
                break;
            case nt_mizs:
                seged_tolt = Tamb + rvt(absT);
                break;
            case nt_erno:
                seged_tolt = Tamb;
                break;
            case nt_polinom:
                seged_tolt = Tamb;
                break;
            case nt_jani_1:
                seged_tolt = Tamb;
                break;
            case nt_diode_1:
                TODO("set_Tamb: nt_diode_1");
                break;
            default:
                TODO("set_Tamb: unknown type");
        }
    }

    //***********************************************************************
    void beolvas_fajlbol(const fazisvalto & fv);
    //***********************************************************************
    void beolvas_fajlbol(); // junction-hoz
    //***********************************************************************
};


//***********************************************************************
struct adat_anyag {
//***********************************************************************

    //***********************************************************************
    tulajdonsag elvez, thvez;       // tulajdonsag
    tulajdonsag Cth, epsz_rel, S;   // tulajdonsag
    tulajdonsag D;                  // tulajdonsag, 0..1, az esõ elektromos telj mekkora része fût
    tulajdonsag emissivity;         // tulajdonsag, 0..1, def 1.
    fazisvalto fazisvaltas;         // fazisvalto
    bool is_fenypor;
    //***********************************************************************

    //***********************************************************************
    void set_Tamb(rvt Tamb) {
    //***********************************************************************
        elvez.set_Tamb(Tamb);
        thvez.set_Tamb(Tamb);
        Cth.set_Tamb(Tamb);
        epsz_rel.set_Tamb(Tamb);
        S.set_Tamb(Tamb);
        D.set_Tamb(Tamb);
        emissivity.set_Tamb(Tamb);
        fazisvaltas.set_Tamb(Tamb);
    }

    //***********************************************************************
    void beolvas_fajlbol(uns hanyas);
    adat_anyag() :is_fenypor{ false } {}
    //***********************************************************************
};


//***********************************************************************
struct adat_junction {
//***********************************************************************
    rvt area;
    tulajdonsag el_egyenlet;            // tulajdonsag
    tulajdonsag D;                      // tulajdonsag
    tulajdonsag R;                      // tulajdonsag, Dissz=U^2*G*D-F*R, F a számíott sug. telj.
    tulajdonsag rad, lum;               // tulajdonsag, rad, lum egyenlete, konst vagy erno lehet
    adat_junction() :area{ rvt() }, D{ rvt(1) }, R{ rvt(1) } {}

    //***********************************************************************
    void set_Tamb(rvt Tamb) {
    //***********************************************************************
        el_egyenlet.set_Tamb(Tamb);
        D.set_Tamb(Tamb);
        R.set_Tamb(Tamb);
        rad.set_Tamb(Tamb);
        lum.set_Tamb(Tamb);
    }

    //***********************************************************************
    void beolvas_fajlbol(uns hanyas);
    //***********************************************************************
};


//***********************************************************************
struct adat_peremfeltetel {
//***********************************************************************
    peremtipus tipus;   // pt_open, pt_u, pt_t, pt_htc, pt_thtc
    rvt UT, UT_tolt;    // rvt
    rvt htc;            // rvt
    void set_Tamb(rvt Tamb) { UT_tolt = (tipus == pt_t || tipus == pt_thtc) ? UT - Tamb : UT; }
    adat_peremfeltetel() : tipus{ pt_open }, UT { rvt() }, UT_tolt{ rvt() }, htc{ rvt() } {}
};


//***********************************************************************
struct adat_gerj {
//***********************************************************************
    gerjesztes_tipus tipus; // gt_none, gt_U, gt_I, gt_T, gt_P
    uns struktura_index;    // melyik struktúrára adandó
    rvt ertek, ertek_tolt;  // inicializáláskor értéket, kilovasáskor toltat
    void set_Tamb(rvt Tamb) { ertek_tolt = (tipus == gt_T) ? ertek - Tamb : ertek; }
    adat_gerj() :struktura_index{ 0 }, ertek{ rvt() }, ertek_tolt{ rvt() } {}
};


//***********************************************************************
struct adat_eredm {
//***********************************************************************
    struct CurrIndex {
        uns cella_index;
        uns face_index;
        CurrIndex(uns ci = 0, uns fi = 0) :cella_index{ ci }, face_index{ fi }{}
    };
    eredmeny_tipus tipus; // et_c_pontprobe, et_f_pontprobe, et_currentprobe, et_c_map, et_f_map
    mit_ment_tipus mit_ment; // mm_UT, MM_IP, mm_lum, mm_rad
    uns cella_index; // pontprobe esetén a cella indexe
    uns face_index; // pontprobe és face esetén a face indexe
    vektor<CurrIndex> currentProbe; // faces of the current probe if this is a current probe
    adat_eredm() :cella_index{ 0 }, face_index{ 0 } {}
};

//***********************************************************************
struct adat_struktura {
//***********************************************************************
    rvt volume;
    const adat_gerj *p_el_gerj, *p_th_gerj;  // adat_gerj *, a vezérlõ állítja a másolatában, a bemenetben nincs funkciója
    adat_struktura() :volume{ rvt() }, p_el_gerj{ nullptr }, p_th_gerj{ nullptr } {}
};


//***********************************************************************
struct adat_anal_beall {
//***********************************************************************
    analizis_beallitas_tipus tipus; // abt_I0, abt_max_error, abt_max_iter, abt_del_all_prev, abt_del_all_gerj
    uns uns_ertek, uns_ertek_2;
    rvt rvt_ertek, rvt_ertek_2;
    adat_anal_beall() :uns_ertek{ 0 }, rvt_ertek{ rvt() }, uns_ertek_2{ 0 }, rvt_ertek_2{ rvt() } {}
    void clear() { uns_ertek = uns_ertek_2 = 0; rvt_ertek = rvt_ertek_2 = rvt(); }
    // abt_I0            uns_ertek: -,           uns_ertek_2: -,  rvt_ertek: I0,        rvt_ertek_2: - 
    // abt_max_error     uns_ertek: -,           uns_ertek_2: -,  rvt_ertek: max error, rvt_ertek_2: -
    // abt_max_iter      uns_ertek: max iter,    uns_ertek_2: -,  rvt_ertek: -,         rvt_ertek_2: -
    // abt_del_all_prev  uns_ertek: -,           uns_ertek_2: -,  rvt_ertek: -,         rvt_ertek_2: -
    // abt_del_all_gerj  uns_ertek: -,           uns_ertek_2: -,  rvt_ertek: -,         rvt_ertek_2: -
    // abt_tamb          uns_ertek: -,           uns_ertek_2: -,  rvt_ertek: Tamb,      rvt_ertek_2: -
    // abt_change_bn     uns_ertek: perem,       uns_ertek_2: uj, rvt_ertek: -,         rvt_ertek_2: -
    // abt_change_bv     uns_ertek: perem,       uns_ertek_2: -,  rvt_ertek: új érték,  rvt_ertek_2: HTC + T peremnél a T
    // abt_change_jn     uns_ertek: junction,    uns_ertek_2: uj, rvt_ertek: -,         rvt_ertek_2: -
    // abt_change_mn     uns_ertek: anyag,       uns_ertek_2: uj, rvt_ertek: -,         rvt_ertek_2: -
    // abt_ct            uns_ertek: 1/2/3,       uns_ertek_2: -,  rvt_ertek: -,         rvt_ertek_2: -
};


//***********************************************************************
struct adat_analizis_lepes {
//***********************************************************************
    analizis_lepes_tipus tipus;             // alt_dc, alt_trans, alt_ac
    rvt value;                              // tranziens lépésköz / AC frekvencia
    vektor<adat_gerj> gerjesztesek;         // adat_gerj vektor, a 0 indexû is érvényes!
    vektor<adat_anal_beall> beallitasok;    // adat_anal_beall vektor, a 0 indexû is érvényes!
    vektor<adat_eredm> eredmenyek;          // adat_eredme vektor, a 0 indexû is érvényes!
    bool is_ignore_error;                   // bool, hagyja jóvá a lépést a vezérlõ akkor is, ha hiba van

    //***********************************************************************
    void beolvas_fajlbol(uns hanyas);
    //***********************************************************************
    adat_analizis_lepes() :value{ rvt() }, is_ignore_error{ false } {}
    //***********************************************************************
    void set_Tamb(rvt Tamb) {
    //***********************************************************************
        for (uns i = 0; i < gerjesztesek.size(); i++) {
            gerjesztesek[i].set_Tamb(Tamb);
        }
    }
};


//***********************************************************************
struct adat_cella;
//***********************************************************************


//***********************************************************************
struct adat_face {
//***********************************************************************
    bool is_el;                             // bool
    bool is_perem;                          // bool, ha a típus ft_normal_perem vagy ft_spec_perem
    bool is_kulso;                          // bool, ha a típus ft_csatlakozo vagy ft_centroid
    face_tipus tipus;                       // face_tipus, ft_csatlakozo, ft_normal_perem, ft_spec_perem, ft_kozepponti, ft_centroid
    uns parja;                              // uns, is_parja_peremre-vel együtt érvényes. nincs=0. ha el, akkor a termikus párja és fordítva. ha a cella egyterû, akkor másik cellában is lehet
    uns anyag_index;                        // uns
    uns zaj_index;                          // rvt, 0, akkor az általános zaj érték szerint kerül beállításra
    uns junction_index;                     // uns, nincs=0
    uns perem_index;                        // uns, open=0
    uns A_index, L_index;                   // uns
    uns csatlakozo_index;
    peremirany_tipus peremirany;            // peremirany_tipus, ha special perem, pit_west, pit_east, pit_south, pit_north, pit_bottom, pit_top
    peremirany_tipus oldal;                 // peremirany_tipus, ha WESNBT oldal, akkor jelzi, ha nem olyan, akkor pit_none
    uns perem_x, perem_y;                   // uns, ha special perem
    adat_face() : is_el{ false }, is_perem{ false }, is_kulso{ false }, tipus{ ft_none }, parja{ 0 }, anyag_index{ 0 }, zaj_index{ 0 }, 
        junction_index{ 0 }, perem_index{ 0 }, A_index{ 0 }, L_index{ 0 }, csatlakozo_index{ 0 }, peremirany{ pit_west }, oldal{ pit_none }, perem_x{ 0 }, perem_y{ 0 } {}
    //***********************************************************************
    void beolvas_fajlbol(uns hanyas, const adat_cella & cella);
    //***********************************************************************
};


//***********************************************************************
struct adat_besugarzo_cella {
//***********************************************************************
    uns cella_index;                        // uns, a sugárzást kibocsátó cella indexe
    tulajdonsag arany;                      // tulajdonsag, a kibocsátott sugárzás mekkora része alakul itt hõvé
};

/*
//***********************************************************************
struct light_path {
//***********************************************************************
    //***********************************************************************
    struct koztes_cella {
    //***********************************************************************
        uns cella_index;
        rvt d;                                // rvt, akt cella vastagsága a fény irányában
        uns blue_absorption_coeff_index;      // uns, kék_ki=kék_be*exp(-ez*d)
        uns conversion_efficiency_index;      // uns, sárga=elnyelt_kék*ez, a többi P_diss
        uns yllw_absorption_coeff_index;      // uns, sárga_ki=sárga_be*exp(-ez*d)+új_sárga
        rvt K;                                // rvt, kimenõ kékjének ezen az útvonalon figyelembe veendõ aránya
        rvt US;                               // rvt, új sárgájának ezen az útvonalon figyelembe veendõ aránya 
        rvt RS;                               // rvt, kimenõ régi sárgájának ezen az útvonalon figyelembe veendõ aránya
    };
    uns junction_cella_index;
    uns junction_cella_face_index;
    rvt junction_cella_K;                     // rvt, forrás face-bõl milyen arányban jön ezen az úton
    rvt akt_d;                                // rvt, akt cella vastagsága a fény irányában
    uns akt_blue_absorption_coeff_index;      // uns, kék_ki=kék_be*exp(-ez*d)
    uns akt_conversion_efficiency_index;      // uns, sárga=elnyelt_kék*ez, a többi P_diss
    uns akt_yllw_absorption_coeff_index;      // uns, sárga_ki=sárga_be*exp(-ez*d)+új_sárga
    vektor<koztes_cella> koztes_cellak;       // vektor<koztes_cella>, 0 indexû is érvényes
};


//***********************************************************************
struct yellow_light_path {
//***********************************************************************
    //***********************************************************************
    struct koztes_cella {
    //***********************************************************************
        uns cella_index;
        rvt d;                                // rvt, akt cella vastagsága a fény irányában
        uns yllw_absorption_coeff_index;      // uns, sárga_ki=sárga_be*exp(-ez*d)+új_sárga
        rvt RS;                               // rvt, kimenõ régi sárgájának ezen az útvonalon figyelembe veendõ aránya
    };
    uns source_cella_index;
    rvt source_cella_K;                       // rvt, forrás cellából milyen arányban jön ezen az úton
    rvt akt_d;                                // rvt, akt cella vastagsága a fény irányában
    uns akt_yllw_absorption_coeff_index;      // uns, sárga_ki=sárga_be*exp(-ez*d)+új_sárga
    vektor<koztes_cella> koztes_cellak;       // vektor<koztes_cella>, 0 indexû is érvényes
};
*/


//***********************************************************************
struct adat_cella {
//***********************************************************************
    //***********************************************************************
    struct fenyutak_t {
    //***********************************************************************
#pragma pack(2)
        //***********************************************************************
        struct path_index_t{
        //***********************************************************************
            uns path_index;                 // uns, melyik path
            unsigned short cella_index;     // uns, path-on belül ennek a fénypor cellának az indexe
            path_index_t() :path_index{ 0 }, cella_index{ 0 } {}
        };
#pragma pack()
        //***********************************************************************
        uns kek_db, sarga_db;               // uns, a cellához kapcsolódó fényutak száma (a vektor átméretezéséhez számoljuk)
        vektor<path_index_t> kek_indexek;   // vektor<path_index_t>, a 0 indexû is érvényes
        vektor<path_index_t> sarga_indexek; // vektor<path_index_t>, a 0 indexû is érvényes
        fenyutak_t() :kek_db{ 0 }, sarga_db{ 0 } {}
    };
    //***********************************************************************
    mezo_tipus mezotipus;                       // mt_elektromos, mt_termikus, mt_elektrotermikus
    cellatipus cella_tipus;                     // ct_faces_cella, ct_th6
    uns struktura_index;                        // uns
    uns zaj_index;                              // uns, a középponti face-ek zaja
    uns anyag_index;                            // uns
    uns kapcsolodo_cella_index;                 // uns, ha eltherm szimulációban különválasztjuk az elektromos cellát a termikustól, akkor a párja
    uns volume_index;                           // uns
    vektor<adat_face> facek;                    // adat_face vektor, a 0 indexû dummy
    vektor<adat_besugarzo_cella> besugarzo_cella; // adat_besugarzo_cella vektor, ha size>0, akkor van olyan cella, ahonnan itt disszipálódik valami
    fenyutak_t fenyutak;                        // fenyutak_t
    adat_cella() :cella_tipus{ ct_none }, struktura_index { 0 }, zaj_index{ 0 }, anyag_index{ 0 }, kapcsolodo_cella_index{ 0 },
        volume_index{ 0 } {}
    //***********************************************************************
    void beolvas_fajlbol(uns hanyas);
    //***********************************************************************
    void set_Tamb(rvt Tamb) {
    //***********************************************************************
        for (uns i = 0; i < besugarzo_cella.size(); i++) {
            besugarzo_cella[i].arany.set_Tamb(Tamb);
        }
    }
};


//***********************************************************************
struct adat_masolando {
//***********************************************************************
    bool is_masolando;                      // bool, ha true, akkor egy mátrixból másolunk, ha false, akkor két mátrixból összeadunk (redukció vagy centroid)
    bool is_centroid;                       // bool, ha is_masolando=false, akkor ha is_centroid=true, akkor centroid (A-ba megy), ha false, akkor redukálandó (B-be megy)
                                            // megj.: mivel másolandóból sokkal több van, mint redukálandóból, ezért az is_centroid vizsgálata csak utóbbinál szükséges
    bool is_balbol;                         // bool, másolandó esetén a forrás a bal vagy a jobb?
    uns hova;                               // A-ban vagy B-ben hová másol (másolandó vagy centroid: A-ban, redukálandó: B-ben), 0-tól indul
    uns bal_from;                           // uns, kezdõindex a bal ágon lévõ mátrixban (0-tól indul)
    uns jobb_from;                          // uns, kezdõindex a jobb ágon lévõ mátrixban (0-tól indul)
    uns db;                                 // uns, a másolandó/összeadandó csomópontok száma (centroidnál mindig 1!)
    adat_masolando() :is_masolando{ false }, is_centroid{ false }, is_balbol{ false }, hova{ 0 }, bal_from{ 0 }, jobb_from{ 0 }, db{ 0 } {}
};


//***********************************************************************
struct adat_fa_elem{
//***********************************************************************
    uns cella_index;                        // uns, ha levél elem, akkor a kapcsolódó cella indexe, ha belsõ elem, akkor 0
    uns tartomany_kezdoindex;               // uns, a tartomany_kezdoindex-tõl a this-ig vannak ebben a tartományban a faelemek
    uns bal_elem_indexe;                    // uns, a bal ágon lévõ bemenõ faelem indexe
    uns jobb_elem_indexe;                   // uns, a jobb ágon lévõ bemenõ faelem indexe
    uns A, B;                               // uns, A és B mátrix mérete (a centroid A-ban van)
    vektor<adat_masolando> masolando;       // adat_masolando vektor, a 0 indexû is érvényes!
    adat_fa_elem() : cella_index{ 0 }, tartomany_kezdoindex{ 0 }, bal_elem_indexe{ 0 }, jobb_elem_indexe{ 0 }, A{ 0 }, B{ 0 } {}
    //***********************************************************************
    void beolvas_fajlbol(uns hanyas);
    //***********************************************************************
};


//***********************************************************************
struct adat_egyedi {
//***********************************************************************
    uns index;
    uns threads;    // a redukciónál ennyi szálra bonttható a mátrixszorzás/invertálás
    adat_egyedi() :index{ 0 }, threads{ 1 } {}
    adat_egyedi(uns index, uns threads) :index{ index }, threads{ threads } {}
};


//***********************************************************************
enum significance { not_significant, half_significant, full_significant };
//***********************************************************************


//***********************************************************************
struct light_path_blue_t {
//***********************************************************************
#pragma pack(2)
    //***********************************************************************
    struct cella {
    //***********************************************************************
        uns cella_index;
        unsigned short ud;                              // rvt, akt cella vastagsága a fény irányában
        unsigned short uK;                              // rvt, kimenõ kékjének ezen az útvonalon figyelembe veendõ aránya
        unsigned short uUS;                             // rvt, új és reemitted sárgájának ezen az útvonalon figyelembe veendõ aránya
        unsigned short uRS;                             // rvt, kimenõ régi sárgájának ezen az útvonalon figyelembe veendõ aránya
        unsigned short blue_absorption_coeff_index;     // uns, kék_ki=kék_be*exp(-ez*d)
        unsigned short conversion_efficiency_index;     // uns, sárga=elnyelt_kék*ez, a többi P_diss
        unsigned short yllw_absorption_coeff_index_sig; // uns, sárga_tovább=sárga_be*exp(-ez*d), alsó két bit: significant
                                                        //      elnyelt_sárga=sárga_be-sárga_tovább
        unsigned short re_conversion_efficiency_index;  // uns, re_sárga=elnyelt_sárga*ez, a többi P_diss
                                                        //      sárga_ki=sárga_tovább+új_sárga+re_sárga
        cella() :ud{ 0 }, uK{ 0 }, uUS{ 0 }, uRS{ 0 }, cella_index{ 0 }, blue_absorption_coeff_index{ 0 }, conversion_efficiency_index{ 0 },
            yllw_absorption_coeff_index_sig{ full_significant }, re_conversion_efficiency_index{ 0 } {}
    };
#pragma pack()
    vektor<cella> cells;                      // vektor<cella>, a 0 indexû a junction vagy fénypor kiinduló cella
    uns face_index;                           // uns, csak kék sugárnál, akkor is csak a junction cellánál
    float d_mul, K0;                          // float, max dossz a pathban per 32768
    light_path_blue_t() : face_index{ 0 }, d_mul{ 0 }, K0{ 0 } {}
    //***********************************************************************
    void beolvas_fajlbol();
    //***********************************************************************
};


//***********************************************************************
struct light_path_yellow_t {
//***********************************************************************
#pragma pack(2)
//***********************************************************************
    struct cella {
    //***********************************************************************
        uns cella_index;
        unsigned short ud;                      // rvt, akt cella vastagsága a fény irányában
        unsigned short uRS;                     // rvt, kimenõ régi sárgájának ezen az útvonalon figyelembe veendõ aránya
        unsigned short yllw_absorption_coeff_index_sig; // uns, sárga_tovább=sárga_be*exp(-ez*d), alsó két bit: significant
                                                        //      elnyelt_sárga=sárga_be-sárga_tovább
        unsigned short re_conversion_efficiency_index;  // uns, re_sárga=elnyelt_sárga*ez, a többi P_diss
                                                        //      sárga_ki=sárga_tovább+új_sárga+re_sárga
        cella() :ud{ 0 }, uRS{ 0 }, cella_index{ 0 }, yllw_absorption_coeff_index_sig{ full_significant }, re_conversion_efficiency_index{ 0 } {}
    };
#pragma pack()
    vektor<cella> cells;                      // vektor<cella>, a 0 indexû a junction vagy fénypor kiinduló cella
    float d_mul, R0;                          // float, max dossz a pathban per 32768
    light_path_yellow_t() : d_mul{ 0 }, R0{ 0 } {}
    //***********************************************************************
    void beolvas_fajlbol();
    //***********************************************************************
};


//***********************************************************************
struct light_paths_t {
//***********************************************************************
    vektor<light_path_blue_t> kek_fenyutak;         // vektor<light_path_blue_t>, a 0 indexû dummy
    vektor<light_path_yellow_t> sarga_fenyutak;     // vektor<light_path_yellow_t>, a 0 indexû dummy
    vektor<tulajdonsag> fenyut_tulajdonsagok;       // tulajdonsag vektor, a 0 indexû dummy
    bool is_write_output;                           // bool, ha nem az összegzett mennyiségek, hanem a sugárvégi mennyiségek írandók
    //***********************************************************************
    void beolvas_fajlbol();
    void cellak_beallitasa(vektor<adat_cella> & cellak);
    //***********************************************************************
    light_paths_t() :is_write_output{ false } {}
    //***********************************************************************
};


//***********************************************************************
struct adat_szimulacio {
//***********************************************************************
    bool is_KM, is_KC, is_KJ, is_KP;                // bool, az elõzõ szimulációban definiált anyagok, struktúrák, junctionök, peremfeltételek,
    bool is_KA, is_KL, is_KD, is_KF, is_KR;         // bool, analízis lépések, cellák, méretek, fák érvényesek.
    bool is_joint_cells;                            // bool, ha van a cellák között kapcsolódó
    bool is_use_commas;                             // bool, ha a szöveges eredményfájlban tizedesvesszõt használunk
    bool is_zaj_random;                             // bool, ha true, akkor srand(time), egyébként srand(1976)
    rvt zaj;                                        // rvt, deafult: 0
    rvt dlppb, dlppy;                                // rvt, deafult: 1
    fa_adat_tipus fa_adat;                          // fa_adat_tipus, fat_double, fat_float, fat_double_ended_float
    ::std::string neve;                             // A szimuláció neve
    ::std::string eredm_utvonal;                    // ide menti az eredményeket
    mezo_tipus mezo;                                // mt_elektromos, mt_termikus, mt_elektrotermikus
    nemlin_megoldo_tipus nemlin_tipus;              // nmt_klasszik_iteracio, nmt_el_th_newton
    solver_tipus solver_type;                       // st_sunred, st_iter
    vektor<adat_anyag> anyagok;                     // adat_anyag vektor, a 0 indexû dummy
    vektor<adat_struktura> strukturak;              // adat_struktura vektor, a 0 indexû dummy
    vektor<adat_junction> junctionok;               // adat_junction vektor, a 0 indexû dummy
    vektor<adat_peremfeltetel> perem;               // adat_peremfeltetel vektor, a 0 indexû dummy, de a 0 index az opent jelenti
    vektor<adat_analizis_lepes> analizis_lepesek;   // adat_analizis_lepes vektor, a 0 indexû dummy
    vektor<adat_cella> cellak;                      // adat_cella vektor, a 0 indexû dummy
    light_paths_t fenyutak;                         // light_paths_t
    vektor<rvt> meretek;                            // rvt vektor, a 0 indexû dummy (a volume_indexhez, A_indexhez, L_indexhez tartozó méretek)
    vektor<klaszter_tartomany> cella_klaszter_tartomanyok; // klaszter_tartomany vektor, a 0 indexû is hasznos! A cellak vektoron belüli tartományt adja!
    uns csatlakozo_db;                              // uns, a csatlakozó csomópontok száma
    uns db_fa;                                      // uns, 1 vagy 2
    uns gyoker_1_index, gyoker_2_index;             // uns, a két fa gyökerének indexe a fa_elemek tömbben
    vektor<adat_fa_elem> fa_elemek;                 // adat_fa_elem vektor, a 0 indexû dummy
    vektor<klaszter_tartomany> fa_klaszter_tartomanyok_1, fa_klaszter_tartomanyok_2; // klaszter_tartomany vektor, a 0 indexû is hasznos! A fa_elemek vektoron belüli tartományt adja!
    vektor<adat_egyedi> egyedi_fa_elemek_1, egyedi_fa_elemek_2; // adat_egyedi vektor, a 0 indexû is hasznos! A fa_elemek vektoron belüli indexet és a megengedett szálszámot adja!
    adat_szimulacio() :csatlakozo_db{ 0 }, db_fa{ 0 }, gyoker_1_index{ 0 }, gyoker_2_index{ 0 }, is_KA{ false }, is_KC{ false }, is_KD{ false }, is_KF{ false }, is_KR{ false },
        is_KJ{ false }, is_KL{ false }, is_KM{ false }, is_KP{ false }, is_joint_cells{ false }, is_zaj_random{ false }, fa_adat{ fat_double }, is_use_commas{ false }, zaj{ 0 }, dlppb{ 1 }, dlppy{ 1 } {}
    //***********************************************************************
    void beolvas_fajlbol(uns hanyas, adat_szimulacio * elozo, uns cputhreads);
    //***********************************************************************
    void kapcsolodo_cellak_beallitasa();
    //***********************************************************************
    void cellatartomany_feltolto(vektor<klaszter_tartomany> & klaszter_tartomanyok, uns szal_kell);
    //***********************************************************************
    void fa_rekurziv_klasztertartomany_feltolto(uns kezdoindex, vektor<klaszter_tartomany> & klaszter_tartomanyok, vektor<adat_egyedi> & egyedi_fa_elemek, uns szal_kell);
    //***********************************************************************
    void set_Tamb(rvt Tamb) {
    //***********************************************************************
        for (uns i = 1; i < anyagok.size(); i++) {
            anyagok[i].set_Tamb(Tamb);
        }
        for (uns i = 1; i < junctionok.size(); i++) {
            junctionok[i].set_Tamb(Tamb);
        }
        for (uns i = 1; i < perem.size(); i++) {
            perem[i].set_Tamb(Tamb);
        }
        for (uns i = 1; i < analizis_lepesek.size(); i++) {
            analizis_lepesek[i].set_Tamb(Tamb);
        }
        for (uns i = 1; i < cellak.size(); i++) {
            cellak[i].set_Tamb(Tamb);
        }
        for (uns i = 1; i < fenyutak.fenyut_tulajdonsagok.size(); i++) {
            fenyutak.fenyut_tulajdonsagok[i].set_Tamb(Tamb);
        }
    }
};


//***********************************************************************
struct bemenet{
//***********************************************************************
    //***********************************************************************
    vektor<adat_szimulacio> szimulaciok; // adat_szimulacio vektor, a 0 indexû dummy
    uns cputhreads;
    ::std::string utvonal;
    //***********************************************************************
    bemenet() : cputhreads{ 8 } {}
    //***********************************************************************
    void felepit_fajlbol(const ::std::string & fajlnev);  // beolvassa a szimulációkat egy fájlból és felépíti a bemeneti adatszerkezetet
    //***********************************************************************
    void set_Tamb(rvt Tamb) {
    //***********************************************************************
        for (uns i = 1; i < szimulaciok.size(); i++) {
            szimulaciok[i].set_Tamb(Tamb);
        }
    }
};


//***********************************************************************
struct iter_csomopont {
//***********************************************************************
    uns cella_1_index, cella_2_index;
    uns cella_1_sajat, cella_2_sajat; // Az admittanciamátrix mely sora (=soron belül mely érték)
    iter_csomopont() :cella_1_index{ 0 }, cella_2_index{ 0 }, cella_1_sajat{ 0 }, cella_2_sajat{ 0 } {}
};


}
#endif
