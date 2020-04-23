//***********************************************************************
// bemeneti adatok header
// Creation date:  2018. 07. 17.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef ADATTIPUSOK_HEADER
#define	ADATTIPUSOK_HEADER
//***********************************************************************


//***********************************************************************
#include "kozos.h"
#include "hiba.h"
#include "vektor.hpp"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
struct tol_db {
//***********************************************************************
    meret_t tol, db;  // kezdõindex és innen kezdve hány darab
};


//***********************************************************************
struct mit_hova_masol {
//***********************************************************************
    bool is_be1;  // true: a be1-rõl kell másolni, false: a be2-rõl kell másolni (0. szinten az alapcella végzi a beállítást, ott nincs mit_hova_masol)
    uns honnan, hova, hanyat; // hány csomópont van az adott oldalon, és az hol kezdõdik a forrás és a cél mátrixban
};


//***********************************************************************
struct ertek_t {
//***********************************************************************
    rvt ertek;
    rvt derivalt;
    ertek_t() :ertek{ rvt() }, derivalt{ rvt() } {}
    ertek_t(rvt ertek, rvt derivalt) :ertek{ ertek }, derivalt{ derivalt } {}
};


//***********************************************************************
struct adatpar {
//***********************************************************************
    rvt T, T_tolt, value; // T: tényleges hõmérséklet, T_tolt: T-Tamb, value: érték
    adatpar(rvt T = rvt(), rvt value = rvt()) : T{ T }, T_tolt{ T }, value{ value } {}
    void eltol(rvt Tamb) { T_tolt = T - Tamb; }
};


//***********************************************************************
class fazisvalto;
//***********************************************************************


//***********************************************************************
class broken_line {
//***********************************************************************
    vektor<adatpar> pontok; // a 0 indexû is érvényes
public:

    //***********************************************************************
    void add_pont(rvt T, rvt value) { // Utána set_Tamb futtatandó
    //***********************************************************************
        pontok.push_back(adatpar(T, value));
    }

    //***********************************************************************
    void set_Tamb(rvt Tamb) {
    //***********************************************************************
        for (uns i = 0; i < pontok.size(); i++)
            pontok[i].eltol(Tamb);
    }

    //***********************************************************************
    ertek_t value(rvt T_tolt) const{
    //***********************************************************************
        if (T_tolt <= pontok.first().T_tolt)
            return ertek_t(pontok.first().value, rvt());
        if (T_tolt >= pontok.last().T_tolt)
            return ertek_t(pontok.last().value, rvt());
        uns masodik = 1;
        while (pontok[masodik].T_tolt < T_tolt)
            masodik++;
        // G0 + (G1 - G0) * (T - T0) / (T1 - T0)
        const adatpar & a0 = pontok[masodik - 1];
        const adatpar & a1 = pontok[masodik];
        rvt m = (a1.value - a0.value) / (a1.T_tolt - a0.T_tolt);
        return ertek_t(a0.value + m * (T_tolt - a0.T_tolt), m);
    }

    //***********************************************************************
    adatpar & operator[](meret_t i) { return pontok[i]; }
    //***********************************************************************
    const adatpar & operator[](meret_t i) const{ return pontok[i]; }
    //***********************************************************************
    void set_size(meret_t n) { pontok.set_size(n); }
    //***********************************************************************
    void beolvas_fajlbol(const fazisvalto & fv);
    //***********************************************************************
};


//***********************************************************************
class fazis_broken_line {
//***********************************************************************
    broken_line gorbe_0, gorbe_1;
public:
    //***********************************************************************
    void add_pont_0(rvt T, rvt value) { gorbe_0.add_pont(T, value); }
    //***********************************************************************
    void add_pont_1(rvt T, rvt value) { gorbe_1.add_pont(T, value); }
    //***********************************************************************
    void set_Tamb(rvt Tamb) { gorbe_0.set_Tamb(Tamb); gorbe_1.set_Tamb(Tamb); }
    //***********************************************************************

    //***********************************************************************
    ertek_t value(rvt T_tolt, const ertek_t & H) const { // H.ertek 0...1, a fázisátalakulás mértéke; H.derivalt a fázisátalakulás meredeksége (vízszintes szakaszon 0!)
    //***********************************************************************
        if (H.ertek == rvt()) return gorbe_0.value(T_tolt);
        if (H.ertek == rvt(1.0)) return gorbe_1.value(T_tolt);
        ertek_t v0 = gorbe_0.value(T_tolt);
        ertek_t v1 = gorbe_1.value(T_tolt);
        rvt szorzo = rvt(1.0) - H.ertek;
        rvt ertek = szorzo * v0.ertek + H.ertek * v1.ertek;
        rvt derivalt = H.derivalt * (v1.ertek - v0.ertek) + szorzo * v0.derivalt + H.ertek * v1.derivalt;
/*
        if (abs(derivalt) > 1.0e+007) {
            printf("ertek=%g\n", ertek);
            printf("derivalt=%g\n", derivalt);
            printf("v0.ertek=%g\n", v0.ertek);
            printf("v0.derivalt=%g\n", v0.derivalt);
            printf("v1.ertek=%g\n", v1.ertek);
            printf("v1.derivalt=%g\n", v1.derivalt);
            printf("H.ertek=%g\n", H.ertek);
            printf("H.derivalt=%g\n\n", H.derivalt);
            getchar();
        }
*/
        return ertek_t(ertek, derivalt);
    }

    //***********************************************************************
    void beolvas_fajlbol(const fazisvalto & fv);
    //***********************************************************************
};


//***********************************************************************
inline ertek_t exp_par(rvt T_tolt, rvt A, rvt B) {
// y  = A*exp(B*T)
// y' = A*B*exp(B*T)
//***********************************************************************
    rvt BT = B * T_tolt;
    if (fabs(BT) <= 7.0) {
        rvt result = A * exp(BT);
        return ertek_t(result, B * result);
    }
    // y  = A*exp( 7)*(B*T+1-7)
    // y' = A*B*exp( 7)
    // y  = A*exp(-7)*(B*T+1+7)
    // y' = A*B*exp(-7)
    if (BT > rvt()) {
        static const rvt exp7 = rvt(exp(7.0));
        const rvt Aexp7 = A * exp7;
        return ertek_t(Aexp7 * (BT - rvt(6.0)), Aexp7*B);
    }
    else {
        static const rvt expm7 = rvt(exp(-7.0));
        const rvt Aexpm7 = A * expm7;
        return ertek_t(Aexpm7*(BT + rvt(8.0)), Aexpm7*B);
    }
}
//***********************************************************************

}
#endif

