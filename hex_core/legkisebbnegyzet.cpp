//***********************************************************************
#include "legkisebbnegyzet.h"
#include "matrix.hpp"
#include "bemenet.h"
#include "fajlolvasas_segito_rutinok.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
void lkn_illeszto(matrix<rvt> & eredmeny, uns U_hatvany_sor, uns T_hatvany_oszlop, const vektor<trio> & vuti) {
// Ha U_hatvany vagy T_hatvany pl. 2, akkor az 3 db szorzót jelent (T^0, T^1, T^2)
//***********************************************************************
    U_hatvany_sor++; // hány darab szorzó kell
    T_hatvany_oszlop++;
    eredmeny.set_size(U_hatvany_sor, T_hatvany_oszlop);
    matrix<rvt> M, Mt, MtM;
    vektor<rvt> B, MtB, A;
    const uns sor = vuti.size();
    const uns oszl = U_hatvany_sor*T_hatvany_oszlop;
    M.set_size(sor, oszl);
    Mt.set_size(oszl, sor);
    MtM.set_size(oszl, oszl);
    B.set_size(sor);
    MtB.set_size(oszl);
    A.set_size(oszl);
    for (uns i = 0; i < sor; i++) {
        trio uti_i = vuti[i];
        B[i] = uti_i.I;
        rvt U_a_j_ediken = 1;
        for (uns j = 0; j < U_hatvany_sor; j++, U_a_j_ediken *= uti_i.U) {
            rvt T_a_k_adikon = 1;
            for (uns k = 0; k < T_hatvany_oszlop; k++, T_a_k_adikon *= uti_i.T) {
                M[i][j * T_hatvany_oszlop + k] = U_a_j_ediken * T_a_k_adikon;
            }
        }
    }
    Mt.transp(M);
    MtM.math_mul_t_nembiztos(Mt, Mt);
    MtM.math_inv_p(false);
    math_mul(MtB, Mt, B);
    math_mul(A, MtM, MtB);
    for (uns i = 0; i < U_hatvany_sor; i++)
        for (uns j = 0; j < T_hatvany_oszlop; j++)
            eredmeny[i][j] = A[i*T_hatvany_oszlop + j];
    printf("\nfit matrix:\n");
    eredmeny.print();
    printf("\nU\tT\tI\tI_fit:\n");
    for (uns i = 0; i < sor; i++) {
        trio uti_i = vuti[i];
        rvt U_a_j_ediken = 1;
        rvt I_szam = 0;
        for (uns j = 0; j < U_hatvany_sor; j++, U_a_j_ediken *= uti_i.U) {
            rvt T_a_k_adikon = 1;
            for (uns k = 0; k < T_hatvany_oszlop; k++, T_a_k_adikon *= uti_i.T) {
                I_szam += A[j*T_hatvany_oszlop + k] * U_a_j_ediken * T_a_k_adikon;
            }
        }
        printf("%g %g %g %g %g%%\n", uti_i.U, uti_i.T, uti_i.I, I_szam, 100 * abs(uti_i.I - I_szam) / uti_i.I);
    }
}


//***********************************************************************
void illesztendo_adatok::beolvas_fajlbol() {
//***********************************************************************
    int ch;
    set_hiba_hol h("illesztendo_adatok::beolvas_fajlbol");

    // szorzó

    szorzo = fajl::get_rvt("fit type junction property multiplier value", true);

    // Trió/Duó

    switch (ch = fajl::get_char("T(rio)/D(uo)")) {
        case 'T': 
            is_trio = true;
            break;
        case 'D': 
            is_trio = false;
            break;
        default:  throw hiba(1, "unknown fit type parameter (%c) (T(rio)/D(uo) expected)", ch);
    }

    // Illesztendõ típus
    // Jelenleg csak polinom támogatott
    // érték és hõmérséklet foka

    fajl::check_text("P");
    egyenlet = lse_polinom;
    fajl::check_text("V");
    egyenlet_unspar_1 = fajl::get_uns("polinom value dimension");
    if (is_trio) {
        fajl::check_text("T");
        egyenlet_unspar_2 = fajl::get_uns("polinom temperature dimension");
    }

    // eleje és vége illesztése

    switch (ch = fajl::get_char("N(one)/L(inear)/S(trong)")) {
        case 'N': 
            start = it_none;
            break;
        case 'L': 
            start = it_lin;
            break;
        case 'S':
            start = it_strong;
            break;
        default:  throw hiba(1, "unknown fit start parameter (%c) (N(one)/L(inear)/S(trong) expected)", ch);
    }

    switch (ch = fajl::get_char("N(one)/L(inear)")) {
        case 'N': 
            stop = it_none;
            break;
        case 'L': 
            stop = it_lin;
            break;
        default:  throw hiba(1, "unknown fit stop parameter (%c) (N(one)/L(inear) expected)", ch);
    }

    // mért pontok

    uns n = fajl::get_uns("number of measured data");
    fajl::check_text("M");
    trio uti;
    for (uns i = 0; i < n; i++) {
        uti.U = fajl::get_rvt("measured U value", true);
        if(is_trio)uti.T = fajl::get_rvt("measured T value", true);
        uti.I = fajl::get_rvt("measured I value", true);
        vuti.push_back(uti);
    }
}


}