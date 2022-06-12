
#include "vezerlo.h"

#include <thread>
#include <chrono>

using namespace ns_v6core;

void kiir_akt_lepes(){
    for (;;){
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        if (akt_lepes::get_akt_lepes_kiir()) {
            const std::string * ps;
            uns szazalek, akt_lepesszam, ossz_lepesszam;
            akt_lepes::get_aktualis_lepes_szazalek(ps, szazalek, akt_lepesszam, ossz_lepesszam);
            if (szazalek == 0)
                log_print(1, "%-60s\r", ps->c_str());
            else
                log_print(1, "%s: %u%%%-50s\r", ps->c_str(), szazalek,"");
        }
    }
}

#include "matrix.hpp"

int main(int n, const char ** params) {
    // printf("%u, %u\n", (uns)sizeof(light_path_blue_t::cella), (uns)sizeof(light_path_yellow_t::cella));
    // return 0;
    time(&prog_start_time);
    try {
/*
        multi_domain_LED_model mod;
        LED_model_result_pack pack;
        mod.calc_form_I(0.5, 25, pack);
        printf("U = %g V\nI = %g A\n", pack.VF, pack.IF);
        mod.calc_form_U(0.7, 25, 1e1, pack);
        printf("U = %g V\nI = %g A\n", pack.VF, pack.IF);
        return 0;
*/
/*
        matrix<double> eredeti, inverz, egyseg, nulla;

        cuns N = 27;
        eredeti.set_size(N, N);
        inverz.set_size(N, N);
        egyseg.set_size(N, N); // egyseg.set_size_szimm(N);
        nulla.set_size_and_zero(N, N);

//        for (uns i = 0; i < N; i++)
//            for (uns j = i; j < N; j++)
//                eredeti[i][j] = eredeti[j][i] = 1.0 / (double)(rand() + 1);

        for (uns i = 0; i < N; i++)
            for (uns j = 0; j < N; j++)
                eredeti[i][j] = 1.0 / (double)(rand() + 1);

        // inverz = eredeti;
        inverz.transp(eredeti);
        //inverz.math_symm_ninv_of_nonsymm();
        //inverz.math_bigblock_symm_in_nonsymm_ninv_np();
        inverz.math_ninv_np_();

//        egyseg.math_add_mul_t_symm(nulla, eredeti, inverz);
//        egyseg.math_add_mul_t_biztos(nulla, eredeti, inverz);
        egyseg.math_bigblock_nonsymm_add_mul_t(2, nulla, eredeti, inverz, false); // !! NEM LEHET, mert nincsenek szálak foglalva !!
        eredeti.print_z();
        inverz.print_z();
        egyseg.print_z();

        return 0;
*/
/*        
        matrix<double> mm, mn, mo, mp, mq;
        cuns N = 4;
        mm.set_size(N, N);
        mn.set_size(N, N);
        mo.set_size(N, N);
        mp.set_size(N, N);
        mq.set_size(N, N);
*/
        /*
        mm[0][0] = 4;
        mm[0][1] = 2;
        mm[0][2] = 1;
        mm[1][0] = 3;
        mm[1][1] = 1;
        mm[1][2] = 1;
        mm[2][0] = 2;
        mm[2][1] = 1;
        mm[2][2] = 0;
        */
/*
        for (uns i = 0; i < N; i++)
            for (uns j = 0; j < N; j++)
                mm[i][j] = 1.0/(double)(rand() + 1);
        printf("Be:\n");
        mm.print_z();printf("\n");
        mn.transp(mm);
        mm.math_ninv_np_();
        printf("Jo inv:\n");
        mm.print_z();printf("\n");
        mp.transp(mn);
        mp.math_ninv_np();
        printf("Rossz inv:\n");
        mp.print_z();printf("\n");
        mo.math_mul_t(mp, mn);
        printf("Egysegmatrix:\n");
        mo.print_z();printf("\n");
        mq.math_sub(mm, mp);
        printf("Jo-rossz:\n");
        mq.print_z();printf("\n");
*/
/*        
        matrix<double> mp;
        mp.set_size(2, 2);
        mp[0][0] = 4;
        mp[0][1] = 3;
        mp[1][0] = 3;
        mp[1][1] = 1;
        mp.math_ninv_np();
        mp.math_symm_ninv_of_nonsymm();
        mp.print_z();printf("\n");
        

        return 0;
*/

        std::thread(kiir_akt_lepes).detach();
        kiirt_ido_forma(true);
        most("start");
        if (n < 2)
            throw hiba("main", "missing sim file name. Usage:\n\t%s <sim file name.v6sim>", params[0]);
        // log_print::set_fajlba(true);
        vezerlo::get_vezerlo().init_from_file(params[1]);
        vezerlo::get_vezerlo().run_simulation();

        most("stop");
        akt_lepes::set_akt_lepes_kiir(false);
        //utolso_esemeny_kiirasa();
        fflush(stdout);
        esemenyek_kiirasa();
        //getchar();
        for (uns i = 0; i < cellak.size(); i++)
            delete cellak.unsafe(i);
    }
    catch (hiba h) {
        printf("\n%s\n", h.what());
    }
    //getchar();
    return 0;
}

#include "legkisebbnegyzet.h"

int main_x() {
    lkn_2_par x;
    //x.push_back(trio(0, 30, 0));
    //x.push_back(trio(0, 50, 0));
    //x.push_back(trio(0, 85, 0));
    //x.push_back(trio(0, 110, 0));
    x.push_back(trio(2.59, 85, 0.02));
    x.push_back(trio(2.614, 70, 0.02));
    x.push_back(trio(2.618, 85, 0.03));
    x.push_back(trio(2.644, 70, 0.03));
    x.push_back(trio(2.647, 50, 0.02));
    x.push_back(trio(2.674, 85, 0.06));
    x.push_back(trio(2.677, 50, 0.03));
    x.push_back(trio(2.685, 30, 0.02));
    x.push_back(trio(2.701, 70, 0.06));
    x.push_back(trio(2.718, 30, 0.03));
    x.push_back(trio(2.726, 85, 0.1));
    x.push_back(trio(2.74, 50, 0.06));
    x.push_back(trio(2.755, 70, 0.1));
    x.push_back(trio(2.785, 30, 0.06));
    x.push_back(trio(2.798, 50, 0.1));
    x.push_back(trio(2.847, 30, 0.1));
    x.push_back(trio(2.932, 85, 0.35));
    x.push_back(trio(2.969, 70, 0.35));
    x.push_back(trio(3.026, 50, 0.35));
    x.push_back(trio(3.026, 85, 0.5));
    x.push_back(trio(3.066, 70, 0.5));
    x.push_back(trio(3.074, 110, 0.7));
    x.push_back(trio(3.096, 30, 0.35));
    x.push_back(trio(3.13, 50, 0.5));
    x.push_back(trio(3.138, 85, 0.7));
    x.push_back(trio(3.183, 70, 0.7));
    x.push_back(trio(3.209, 30, 0.5));
    x.push_back(trio(3.22, 110, 1));
    x.push_back(trio(3.257, 50, 0.7));
    x.push_back(trio(3.293, 85, 1));
    x.push_back(trio(3.346, 70, 1));
    matrix<rvt> par;
    x.calc_parameters(par, 3, 2);
    return 0;
}