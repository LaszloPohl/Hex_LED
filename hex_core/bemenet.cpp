//***********************************************************************
// bemeneti adatok cpp
// Creation date:  2018. 06. 28.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#include "bemenet.h"
#include "fajlolvasas_segito_rutinok.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//****************************************************************
::std::string get_utvonal(const ::std::string& fajlnev){
// bemenet: útvonal + fájlnév
//****************************************************************
    // ::std::string utvonal;
    size_t per = fajlnev.rfind('/');
    size_t bs = fajlnev.rfind('\\');
    if (per == ::std::string::npos && bs == ::std::string::npos)
        return "./";
    if (per == ::std::string::npos)
        return fajlnev.substr(0, bs + 1);
    if (bs == ::std::string::npos)
        return fajlnev.substr(0, per + 1);
    return fajlnev.substr(0, per > bs ? per + 1 : bs + 1);
}

#include <direct.h>
//***********************************************************************
void bemenet::felepit_fajlbol(const ::std::string & fajlnev){
//***********************************************************************
    using namespace std;
    set_hiba_hol h("bemenet::felepit_fajlbol");
    utvonal = get_utvonal(fajlnev);
    _mkdir((utvonal + "hex_res").c_str());
    log_print("Path: %s\n", utvonal.c_str());
    fajl f(fajlnev);
    log_print("Reading: %s\n", fajlnev.c_str());
    fajl::check_text("V6SIMFILE");
    fajl::check_text("CPUTHREADS");
    cputhreads = fajl::get_uns("CPUTHREADS");
    fajl::check_text("NS");
    szimulaciok.set_size(1 + fajl::get_uns("NS simulation number"));
    for (uns i = 1; i < szimulaciok.size(); i++) {
        szimulaciok[i].beolvas_fajlbol(i, i == 1 ? nullptr : &szimulaciok[i - 1], cputhreads);
        std::string eredm_ut = utvonal + "hex_res/" + szimulaciok[i].neve;
        szimulaciok[i].eredm_utvonal = eredm_ut + '/';
        _mkdir((eredm_ut).c_str());

    }
    fajl::check_text("EV");
}


//***********************************************************************
void light_path_blue_t::beolvas_fajlbol() {
//***********************************************************************
    struct prev_values_t {
        const light_path_blue_t * p_path;
        uns C0, D;
        unsigned short B, E, Y, M, uK, uU, uR, F0;
        uns K0;
    };
    static prev_values_t prevs;

    // junction cella
    
    fajl::check_text("P");
    uns D0 = fajl::get_uns("D0");
    if (D0 == 0) D0 = prevs.D;
    else prevs.D = D0;
    d_mul = (float)(D0 * 1e-9 / 32768);

    char ch = fajl::get_char("blue path first char");
    
    uns C0;
    if (ch == 'C') {
        prevs.C0 = C0 = fajl::get_uns("C0 index");
        ch = fajl::get_char("blue path char");
    }
    else C0 = prevs.C0;
    
    if (ch == 'F') {
        prevs.F0 = face_index = fajl::get_uns("F0 index");
        ch = fajl::get_char("blue path char");
    }
    else face_index = prevs.F0;
    
    uns uK0;
    if (ch == 'K') {
        prevs.K0 = uK0 = fajl::get_uns("K0");
        ch = fajl::get_char("blue path char");
    }
    else uK0 = prevs.K0;
    
    K0 = (float)(uK0 * 1e-9f);
    
    if (ch != 'N')
        throw hiba("light_path_blue_t::beolvas_fajlbol", "blue light path, N expected, %c arrived", ch);
    cells.set_size(1 + fajl::get_uns("N blue path cell dimension number"));
    
    cells[0].cella_index = C0;
    cells[0].uK = 0; // ezt nem szabad figyelembe venni a számolásnál, hanem K0-t kell

    // többi cella

    for (uns i = 1; i < cells.size(); i++) {
        ch = fajl::get_char("C_i");
        if (ch == 'C' || ch == 'S' || ch == 'H') {
            if (ch == 'S') cells[i].yllw_absorption_coeff_index_sig = full_significant;
            else if (ch == 'H') cells[i].yllw_absorption_coeff_index_sig = half_significant;
            else cells[i].yllw_absorption_coeff_index_sig = not_significant;
            cells[i].cella_index = fajl::get_uns("C_i index");
            ch = fajl::get_char("B_i");
        }
        else {
            cells[i].yllw_absorption_coeff_index_sig = 3 & prevs.p_path->cells[i].yllw_absorption_coeff_index_sig;
            cells[i].cella_index = prevs.p_path->cells[i].cella_index;
        }
        
        if (ch == 'B') {
            prevs.B = cells[i].blue_absorption_coeff_index = fajl::get_uns("B_i");
            ch = fajl::get_char("E_i");
        }
        else cells[i].blue_absorption_coeff_index = prevs.B;
        
        if (ch == 'E') {
            prevs.E = cells[i].conversion_efficiency_index = fajl::get_uns("E_i");
            ch = fajl::get_char("Y_i");
        }
        else cells[i].conversion_efficiency_index = prevs.E;

        if (ch == 'Y') {
            cells[i].yllw_absorption_coeff_index_sig += prevs.Y = 4 * fajl::get_uns("Y_i");
            ch = fajl::get_char("M_i");
        }
        else cells[i].yllw_absorption_coeff_index_sig += prevs.Y;

        if (ch == 'M') {
            prevs.M = cells[i].re_conversion_efficiency_index = fajl::get_uns("M_i");
            ch = fajl::get_char("K_i");
        }
        else cells[i].re_conversion_efficiency_index = prevs.M;

        if (ch == 'K') {
            float akt_K = (float)fajl::get_rvt("K_i", true);
            prevs.uK = cells[i].uK = (uns)(akt_K * 32768 + 0.5);
            ch = fajl::get_char("U_i");
        }
        else cells[i].uK = prevs.uK;

        if (ch == 'U') {
            prevs.uU = cells[i].uUS = (uns)(fajl::get_rvt("U_i", true) * 32768 + 0.5);
            ch = fajl::get_char("R_i");
        }
        else cells[i].uUS = prevs.uU;

        if (ch == 'R') {
            prevs.uR = cells[i].uRS = (uns)(fajl::get_rvt("R_i", true) * 32768 + 0.5);
            ch = fajl::get_char("I_i");
        }
        else cells[i].uRS = prevs.uR;

        if (ch != 'D')
            throw hiba("light_path_blue_t::beolvas_fajlbol", "blue light path, D expected, %c arrived", ch);
        uns ud = fajl::get_uns("d_i");
        if (ud == 0) ud = 32768;
        
        cells[i].ud = ud;//(float)(ud * d_max);
    }
    prevs.p_path = this;
}


//***********************************************************************
void light_path_yellow_t::beolvas_fajlbol() {
//***********************************************************************
    struct prev_values_t {
        const light_path_yellow_t * p_path;
        uns C0, D;
        unsigned short Y, M, uR;
        uns R0;
    };
    static prev_values_t prevs;

    // forrás cella

    fajl::check_text("P");
    uns D0 = fajl::get_uns("D0");
    if (D0 == 0) D0 = prevs.D;
    else prevs.D = D0;
    d_mul = (float)(D0 * 1e-9 / 32768);

    char ch = fajl::get_char("blue path first char");

    uns C0;
    if (ch == 'C') {
        prevs.C0 = C0 = fajl::get_uns("C0 index");
        ch = fajl::get_char("blue path char");
    }
    else C0 = prevs.C0;

    uns uR0;
    if (ch == 'R') {
        prevs.R0 = uR0 = fajl::get_uns("R0");
        ch = fajl::get_char("blue path char");
    }
    else uR0 = prevs.R0;

    R0 = uR0 * 1e-9f;

    if (ch != 'N')
        throw hiba("light_path_blue_t::beolvas_fajlbol", "yellow light path, N expected, %c arrived", ch);
    cells.set_size(1 + fajl::get_uns("N blue path cell dimension number"));

    cells[0].cella_index = C0;
    cells[0].uRS = 0; // ezt nem szabad figyelembe venni a számolásnál, hanem R0-t kell

    // többi cella

    for (uns i = 1; i < cells.size(); i++) {
        ch = fajl::get_char("C_i");
        if (ch == 'C' || ch == 'S' || ch == 'H') {
            if (ch == 'S') cells[i].yllw_absorption_coeff_index_sig = full_significant;
            else if (ch == 'H') cells[i].yllw_absorption_coeff_index_sig = half_significant;
            else cells[i].yllw_absorption_coeff_index_sig = not_significant;
            cells[i].cella_index = fajl::get_uns("C_i index");
            ch = fajl::get_char("B_i");
        }
        else {
            cells[i].yllw_absorption_coeff_index_sig = 3 & prevs.p_path->cells[i].yllw_absorption_coeff_index_sig;
            cells[i].cella_index = prevs.p_path->cells[i].cella_index;
        }

        if (ch == 'Y') {
            cells[i].yllw_absorption_coeff_index_sig += prevs.Y = 4 * fajl::get_uns("Y_i");
            ch = fajl::get_char("D_i");
        }
        else cells[i].yllw_absorption_coeff_index_sig += prevs.Y;

        //if (cells[i].yllw_absorption_coeff_index_sig > 12)
        //    printf("bbb");

        if (ch == 'M') {
            prevs.M = cells[i].re_conversion_efficiency_index = fajl::get_uns("M_i");
            ch = fajl::get_char("D_i");
        }
        else cells[i].re_conversion_efficiency_index = prevs.M;

        if (ch == 'R') {
            float akt_R = (float)fajl::get_rvt("R_i", true);
            prevs.uR = cells[i].uRS = (uns)(akt_R * 32768 + 0.5);
            ch = fajl::get_char("X_i");
        }
        else cells[i].uRS = prevs.uR;

        if (ch != 'D')
            throw hiba("light_path_blue_t::beolvas_fajlbol", "blue light path, D expected, %c arrived", ch);
        uns ud = fajl::get_uns("d_i");
        if (ud == 0) ud = 32768;

        cells[i].ud = ud;//(float)(ud * d_max);
    }
    prevs.p_path = this;
}


//***********************************************************************
void light_paths_t::beolvas_fajlbol() {
//***********************************************************************
    char ch1 = fajl::get_char("WO or NB");
    char ch2 = fajl::get_char("WO or NB");

    is_write_output = (ch1 == 'W'&&ch2 == 'O');
    if (is_write_output) {
        ch1 = fajl::get_char("WO or NB");
        ch2 = fajl::get_char("WO or NB");
    }

    // kék utak

    if (ch1 != 'N' || ch2 != 'B')
        throw hiba("light_paths_t::beolvas_fajlbol", "NB expected %c%c arrived", ch1, ch2);
    kek_fenyutak.set_size(1 + fajl::get_uns("NB dimensions number"));
    if (kek_fenyutak.size() > 1) {
        akt_lepes::set_aktualis_lepes_neve("blue paths read");
        uns egyszazalek = kek_fenyutak.size() / 100;
        for (uns i = 1; i < kek_fenyutak.size(); i++) {
            kek_fenyutak[i].beolvas_fajlbol();
            if (i % egyszazalek == 0)
                akt_lepes::set_aktualis_lepes_szazalek(i / egyszazalek);
        }
        akt_lepes::set_aktualis_lepes_szazalek(100);
        log_print::log_print(1, "blue paths read: OK                  \n");
    }

    // sárga utak

    fajl::check_text("NY");
    sarga_fenyutak.set_size(1 + fajl::get_uns("NY dimensions number"));
    if (sarga_fenyutak.size() > 1) {
        akt_lepes::set_aktualis_lepes_neve("yellow paths read");
        uns egyszazalek = sarga_fenyutak.size() / 100;
        for (uns i = 1; i < sarga_fenyutak.size(); i++) {
            sarga_fenyutak[i].beolvas_fajlbol();
            if (i % egyszazalek == 0)
                akt_lepes::set_aktualis_lepes_szazalek(i / egyszazalek);
        }
        akt_lepes::set_aktualis_lepes_szazalek(100);
        log_print::log_print(1, "yellow paths read: OK                  \n");
    }

    // tulajdonságok

    fajl::check_text("NR");
    fenyut_tulajdonsagok.set_size(1 + fajl::get_uns("NR dimensions number"));
    for (uns i = 1; i < fenyut_tulajdonsagok.size(); i++) {
        fajl::check_text("R");
        uns uns_temp = fajl::get_uns("lightpath property R");
        if (uns_temp != i)
            throw hiba(1, "R%u expected, R%u arrived", i, uns_temp);
        fajl::check_text("=");
        fenyut_tulajdonsagok[i].beolvas_fajlbol();
    }
}


//***********************************************************************
void light_paths_t::cellak_beallitasa(vektor<adat_cella> & cellak) {
//***********************************************************************
    uns egyszazalek;

    // darabszámítás

    akt_lepes::set_aktualis_lepes_neve("blue paths cell count");
    egyszazalek = kek_fenyutak.size() / 100;
    for (uns i = 1; i < kek_fenyutak.size(); i++) {
        const light_path_blue_t & akt_ut = kek_fenyutak[i];
        for (uns j = 1; j < akt_ut.cells.size(); j++) {
            cellak[akt_ut.cells[j].cella_index].fenyutak.kek_db++;
        }
        if (i % egyszazalek == 0)
            akt_lepes::set_aktualis_lepes_szazalek(i / egyszazalek);
    }

    akt_lepes::set_aktualis_lepes_neve("yellow paths cell count");
    egyszazalek = sarga_fenyutak.size() / 100;
    for (uns i = 1; i < sarga_fenyutak.size(); i++) {
        const light_path_yellow_t & akt_ut = sarga_fenyutak[i];
        for (uns j = 1; j < akt_ut.cells.size(); j++) {
            cellak[akt_ut.cells[j].cella_index].fenyutak.sarga_db++;
        }
        if (i % egyszazalek == 0)
            akt_lepes::set_aktualis_lepes_szazalek(i / egyszazalek);
    }

    // átméretezés, a darabokat visszaállítja 0-ra, mert az indexek beállításánál innen tudjuk, hogy hányast kell beállítani

    akt_lepes::set_aktualis_lepes_neve("cell path resize");
    egyszazalek = cellak.size() / 100;
    for (uns i = 1; i < cellak.size(); i++) {
        adat_cella::fenyutak_t & akt_fu = cellak[i].fenyutak;
        if (akt_fu.kek_db > 0)
            akt_fu.kek_indexek.set_size(akt_fu.kek_db);
        akt_fu.kek_db = 0;
        if (akt_fu.sarga_db > 0)
            akt_fu.sarga_indexek.set_size(akt_fu.sarga_db);
        akt_fu.sarga_db = 0;
        if (i % egyszazalek == 0)
            akt_lepes::set_aktualis_lepes_szazalek(i / egyszazalek);
    }

    // indexek beállítása

    akt_lepes::set_aktualis_lepes_neve("set blue paths cell indices");
    egyszazalek = kek_fenyutak.size() / 100;
    for (uns i = 1; i < kek_fenyutak.size(); i++) {
        const light_path_blue_t & akt_ut = kek_fenyutak[i];
        for (uns j = 1; j < akt_ut.cells.size(); j++) {
            adat_cella & akt_cella = cellak[akt_ut.cells[j].cella_index];
            auto & akt_index = akt_cella.fenyutak.kek_indexek[akt_cella.fenyutak.kek_db++];
            akt_index.path_index = i;
            akt_index.cella_index = j;
        }
        if (i % egyszazalek == 0)
            akt_lepes::set_aktualis_lepes_szazalek(i / egyszazalek);
    }

    akt_lepes::set_aktualis_lepes_neve("set yellow paths cell indices");
    egyszazalek = sarga_fenyutak.size() / 100;
    for (uns i = 1; i < sarga_fenyutak.size(); i++) {
        const light_path_yellow_t & akt_ut = sarga_fenyutak[i];
        for (uns j = 1; j < akt_ut.cells.size(); j++) {
            adat_cella & akt_cella = cellak[akt_ut.cells[j].cella_index];
            auto & akt_index = akt_cella.fenyutak.sarga_indexek[akt_cella.fenyutak.sarga_db++];
            akt_index.path_index = i;
            akt_index.cella_index = j;
        }
        if (i % egyszazalek == 0)
            akt_lepes::set_aktualis_lepes_szazalek(i / egyszazalek);
    }
    log_print::log_print(1, "setting cells: OK                  \n");
}


//***********************************************************************
void adat_szimulacio::beolvas_fajlbol(uns hanyas, adat_szimulacio * elozo, uns cputhreads){
//***********************************************************************
    set_hiba_hol h("adat_szimulacio::beolvas_fajlbol");
    
    // BS
    
    akt_lepes::set_aktualis_lepes_neve("starting read sim");
    fajl::check_text("BS");
    uns uns_temp = fajl::get_uns("BS");
    if (uns_temp != hanyas)
        throw hiba(1, "BS%u expected, BS%u arrived", hanyas, uns_temp);
    neve = fajl::get_quoted_text("BSx name is missing");
    log_print("Reading simulation %u: \"%s\"\n", hanyas, neve.c_str());
    
    // keep

    int ch;
    do{
        switch (ch = fajl::get_char("keep")) {
            case 'K':  while ((ch = fajl::get_char_or_eol("keep ?")) != 0) {
                    switch (ch) {
                        case 'M': is_KM = true; break;
                        case 'C': is_KC = true; break;
                        case 'J': is_KJ = true; break;
                        case 'P': is_KP = true; break;
                        case 'A': is_KA = true; break;
                        case 'L': is_KL = true; break;
                        case 'D': is_KD = true; break;
                        case 'F': is_KF = true; break;
                        case 'R': is_KR = true; break;
                        default:
                            throw hiba(1, "keep type");
                    }
                }
                ch = 'K';
                break;
            case 'T': 
                switch (ch = fajl::get_char("TFLOAT/TDOUBLE")) {
                    case 'F': fa_adat = fat_float;  break;
                    case 'D': fa_adat = fat_double; break;
                    case 'E': fa_adat = fat_double_ended_float; break;
                    default:  throw hiba(1, "unknown tree data type TF/TD/TE (T%c)", ch);
                }
                ch = 'T';
                break;
            case 'Z':
                switch (ch = fajl::get_char("noise D/R")) {
                    case 'D': is_zaj_random = false; break;
                    case 'R': is_zaj_random = true;  break;
                    default:  throw hiba(1, "unknown noise type D/R (Z%c)", ch);
                }
                zaj = fajl::get_rvt("noise", true);
                ch = 'Z';
                break;
            case 'F': break;
            default:  throw hiba(1, "unknown sim starting parameter (%c)", ch);
        }
    } while (ch != 'F');
    
    // mezõ típus
    
    enum all { all_alap, all_el, all_th, all_elth };
    all a = all_alap;
    nemlin_tipus = nmt_klasszik_iteracio;
    solver_type = st_sunred;//st_iter;//
    while ((ch = fajl::get_char_or_eol("simulation field type")) != 0) {
        switch (ch) {
            case '1': nemlin_tipus = nmt_klasszik_iteracio; break;
            case '2': nemlin_tipus = nmt_el_th_newton; break;
            case 'E': a = all_el; break;
            case 'T': a = (a == all_el) ? all_elth : all_th; break;
            case 'S': solver_type = st_sunred; break;
            case 'I': solver_type = st_iter; break;
            default:
                throw hiba(1, "simulation field type format (%c)", ch);
        }
    }
    switch (a) {
        case all_el:   mezo = mt_elektromos; break;
        case all_th:   mezo = mt_termikus; break;
        case all_elth: mezo = mt_elektrotermikus; break;
        default:
            throw hiba(1, "simulation field type format");
    }

    // anyagok

    if (is_KM) {
        if (elozo == nullptr)
            throw hiba(1, "keep materials is defined but no previous simulation");
        anyagok.rafektet(elozo->anyagok);
    }
    else {
        fajl::check_text("NM");
        anyagok.set_size(1 + fajl::get_uns("NM material number"));
        for (uns i = 1; i < anyagok.size(); i++)
            anyagok[i].beolvas_fajlbol(i);
    }

    // struktúrák

    if (is_KC) {
        if (elozo == nullptr)
            throw hiba(1, "keep structures is defined but no previous simulation");
        strukturak.rafektet(elozo->strukturak);
    }
    else {
        fajl::check_text("NC");
        strukturak.set_size(1 + fajl::get_uns("NC structure number"));
        for (uns i = 1; i < strukturak.size(); i++) {
            fajl::check_text("C", "structure");
            uns uns_temp = fajl::get_uns("structure Ci");
            if (uns_temp != i)
                throw hiba(1, "C%u expected, C%u arrived", i, uns_temp);
            fajl::check_text("V", "structure");
            strukturak[i].volume = fajl::get_rvt("structure volume", false);
        }
    }

    // junctionok

    if (is_KJ) {
        if (elozo == nullptr)
            throw hiba(1, "keep junctions is defined but no previous simulation");
        junctionok.rafektet(elozo->junctionok);
    }
    else {
        fajl::check_text("NJ");
        junctionok.set_size(1 + fajl::get_uns("NJ junction number"));
        for (uns i = 1; i < junctionok.size(); i++)
            junctionok[i].beolvas_fajlbol(i);
    }

    // peremfeltételek

    if (is_KP) {
        if (elozo == nullptr)
            throw hiba(1, "keep boundary conditions is defined but no previous simulation");
        perem.rafektet(elozo->perem);
    }
    else {
        fajl::check_text("NP");
        perem.set_size(1 + fajl::get_uns("NP boundary condition number"));
        for (uns i = 1; i < perem.size(); i++) {
            fajl::check_text("P", "boundary condition");
            uns uns_temp = fajl::get_uns("boundary condition Pi");
            if (uns_temp != i)
                throw hiba(1, "P%u expected, P%u arrived", i, uns_temp);
            int ch;
            switch (ch = fajl::get_char("boundary condition type")) {
                case 'T': perem[i].tipus = pt_t; perem[i].UT = fajl::get_rvt("boundary condition T", false); break;
                case 'U': perem[i].tipus = pt_u; perem[i].UT = fajl::get_rvt("boundary condition U", false); break;
                case 'H': perem[i].tipus = pt_htc; perem[i].htc = fajl::get_rvt("boundary condition H", false); break;
                case 'K': 
                    perem[i].tipus = pt_thtc; 
                    perem[i].htc = fajl::get_rvt("boundary condition K", false); 
                    fajl::check_text("T", "boundary condition K T");
                    perem[i].UT = fajl::get_rvt("boundary condition K T", false);
                    break;
                case 'M': TODO("map boundary condition");  break;
                case 'S': TODO("special boundary condition");  break;
                default:  throw hiba(1, "unknown junction property N(ormal)/I(nverz) (%c)", ch);
            }
        }
    }

    // cellák

    if (is_KL) {
        if (elozo == nullptr)
            throw hiba(1, "keep cells is defined but no previous simulation");
        cellak.rafektet(elozo->cellak);
        is_joint_cells = elozo->is_joint_cells;
    }
    else {
        akt_lepes::set_aktualis_lepes_neve("cells data read");
        fajl::check_text("NL");
        cellak.set_size(1 + fajl::get_uns("NL cells number"));
        fajl::check_text("C");
        csatlakozo_db = fajl::get_uns("connection number");
        if (csatlakozo_db % 2 == 1)
            throw hiba(1, "connection number is odd, it must be even");
        uns egyszazalek = cellak.size() / 100;
        bool is_joint = false;
        for (uns i = 1; i < cellak.size(); i++) {
            cellak.unsafe(i).beolvas_fajlbol(i);
            if (cellak.unsafe(i).kapcsolodo_cella_index != 0)
                is_joint = true;
            if (i % egyszazalek == 0)
                akt_lepes::set_aktualis_lepes_szazalek(i / egyszazalek);
        }
        is_joint_cells = is_joint;
        if(is_joint)
            kapcsolodo_cellak_beallitasa();
        cella_klaszter_tartomanyok.clear();
#ifdef _DEBUG
        cellatartomany_feltolto(cella_klaszter_tartomanyok, 1);
#else
printf("\n cellatartomany_feltolto cputhreads:%u\n", cputhreads);
        cellatartomany_feltolto(cella_klaszter_tartomanyok, cputhreads);
#endif
        log_print::log_print(1, "cells data read: OK                  \n");
    }

    // fenyut_tulajdonsagok

    if (is_KL) {
        if (elozo == nullptr)
            throw hiba(1, "keep cells is defined but no previous simulation");
        fenyutak.kek_fenyutak.rafektet(elozo->fenyutak.kek_fenyutak);
        fenyutak.sarga_fenyutak.rafektet(elozo->fenyutak.sarga_fenyutak);
        fenyutak.fenyut_tulajdonsagok.rafektet(elozo->fenyutak.fenyut_tulajdonsagok);
    }
    else {
        fenyutak.beolvas_fajlbol();
        fenyutak.cellak_beallitasa(cellak);
    }

    // méretek

    if (is_KD) {
        if (elozo == nullptr)
            throw hiba(1, "keep dimensions is defined but no previous simulation");
        meretek.rafektet(elozo->meretek);
    }
    else {
        fajl::check_text("ND");
        meretek.set_size(1 + fajl::get_uns("ND dimensions number"));
        for (uns i = 1; i < meretek.size(); i++) {
            fajl::check_text("D");
            uns_temp = fajl::get_uns("dimensions data D");
            if (uns_temp != i)
                throw hiba(1, "D%u expected, D%u arrived", i, uns_temp);
            fajl::check_text("=");
            meretek[i] = fajl::get_rvt("dimension size", false);
        }
    }

    // fa

    if (is_KF) {
        if (elozo == nullptr)
            throw hiba(1, "keep trees is defined but no previous simulation");
        fa_elemek.rafektet(elozo->fa_elemek);
        // biztos, hogy itt elég a fa elemek ráfektetése, és nem kell pl. a klaszter tartományokat is másolni?
    }
    else {
        akt_lepes::set_aktualis_lepes_neve("tree data read");
        fajl::check_text("NF");
        db_fa = fajl::get_uns("tree number NF");
        fajl::check_text("B");
        fa_elemek.set_size(1 + fajl::get_uns("NF B tree element number"));
        fajl::check_text("G1G");
        gyoker_1_index = fajl::get_uns("G1G");
        if (db_fa == 2) {
            fajl::check_text("G2G");
            gyoker_2_index = fajl::get_uns("G2G");
        }
        if (db_fa != 1 && db_fa != 2)
            throw hiba(1, "1 or 2 tree is supported, %u arrived", db_fa);
        uns egyszazalek = fa_elemek.size() / 100;
        for (uns i = 1; i < fa_elemek.size(); i++) {
            switch (ch = fajl::get_char("tree element")) {
                case 'L': 
                    uns_temp = fajl::get_uns("leaf element");
                    if (uns_temp != i)
                        throw hiba(1, "L%u expected, L%u arrived", i, uns_temp);
                    fajl::check_text("C");
                    fa_elemek[i].cella_index = fajl::get_uns("leaf element C");
                    fajl::check_text("A");
                    fa_elemek[i].A = fajl::get_uns("leaf element A");
                    fajl::check_text("B");
                    fa_elemek[i].B = fajl::get_uns("leaf element B");
                    break;
                case 'B': fa_elemek[i].beolvas_fajlbol(i); break;
                default:  throw hiba(1, "unknown tree element type (%c)", ch);
            }
            if (i % egyszazalek == 0)
                akt_lepes::set_aktualis_lepes_szazalek(i / egyszazalek);
        }
        fa_klaszter_tartomanyok_1.clear();
        fa_klaszter_tartomanyok_2.clear();
        egyedi_fa_elemek_1.clear();
        egyedi_fa_elemek_2.clear();
        fa_rekurziv_klasztertartomany_feltolto(gyoker_1_index, fa_klaszter_tartomanyok_1, egyedi_fa_elemek_1, cputhreads);
        if (db_fa == 2) {
            fa_rekurziv_klasztertartomany_feltolto(gyoker_2_index, fa_klaszter_tartomanyok_2, egyedi_fa_elemek_2, cputhreads);
        }
        log_print::log_print(1, "tree data read: OK                  \n");
    }

    // analízis lépések

    if (is_KA) {
        if (elozo == nullptr)
            throw hiba(1, "keep analysis steps is defined but no previous simulation");
        analizis_lepesek.rafektet(elozo->analizis_lepesek);
    }
    else {
        int ch = fajl::get_char("analysis");
        dlppb = dlppy = 1;
        if (ch == 'D') {
            fajl::check_text("LPP"); 
            ch = fajl::get_char("DLPPB");
            if (ch == 'B') dlppb = fajl::get_rvt("DLPPB", true);
            else if (ch == 'Y') dlppy = fajl::get_rvt("DLPPY", true);
            ch = fajl::get_char("analysis");
        }
        if (ch == 'D') {
            fajl::check_text("LPP");
            ch = fajl::get_char("DLPPB");
            if (ch == 'B') dlppb = fajl::get_rvt("DLPPB", true);
            else if (ch == 'Y') dlppy = fajl::get_rvt("DLPPY", true);
            ch = fajl::get_char("analysis");
        }
        switch (ch) {
            case 'U': fajl::check_text("C"); is_use_commas = true; fajl::check_text("NA"); break;
            case 'N': is_use_commas = false; fajl::check_text("A"); break;
            default:  throw hiba(1, "unknown analysis start command (%c)", ch);
        }
        vektor<adat_analizis_lepes> elo_analizis_lepesek;
        elo_analizis_lepesek.set_size(1 + fajl::get_uns("NA analysis steps number"));
        analizis_lepesek.set_size(1);
        uns plusz_lepes_db = 0;
        const rvt hibalepes_ido = 1e-10;
        for (uns i = 1; i < elo_analizis_lepesek.size(); i++) {
            elo_analizis_lepesek[i].beolvas_fajlbol(i);
            analizis_lepesek.push_back(elo_analizis_lepesek[i]);
/*
            // tranziensnél gerjesztésugráshoz plusz lépés
            if (i > 1 && elo_analizis_lepesek[i].tipus == alt_trans && elo_analizis_lepesek[i].gerjesztesek.size() > 0) {

                analizis_lepesek.last().value = hibalepes_ido;
                analizis_lepesek.last().is_ignore_error = true;
                analizis_lepesek.push_back(elo_analizis_lepesek[i]);

                analizis_lepesek.last().value = hibalepes_ido;
                analizis_lepesek.last().is_ignore_error = true;
                analizis_lepesek.push_back(elo_analizis_lepesek[i]);

                analizis_lepesek.last().value -= 2*hibalepes_ido;
                analizis_lepesek.last().gerjesztesek.clear();
            }
*/
        }
    }

    // ES

    fajl::check_text("ES");
    uns_temp = fajl::get_uns("ES");
    if (uns_temp != hanyas)
        throw hiba(1, "ES%u expected, ES%u arrived", hanyas, uns_temp);
    most("simualtion data read OK");
}


//***********************************************************************
void adat_szimulacio::kapcsolodo_cellak_beallitasa() {
//***********************************************************************
    set_hiba_hol h("adat_szimulacio::kapcsolodo_cellak_beallitasa");
    for (uns i = 1; i < cellak.size(); i++) {
        if (cellak.unsafe(i).kapcsolodo_cella_index != 0) {
            
            // elõször termikus joint cella kell
            
            if (cellak.unsafe(i).mezotipus == mt_termikus) // az elektotermikusakat kiszûrtük beolvasásnál, az nem lehet joint
                throw hiba(1, "thermal joint cell without a preceding joint electrical cell (%u)", i);
            
            // megkeressük az utána következõ elsõ nem elektromos cellát

            uns j;
            for (j = i + 1; j < cellak.size() && cellak.unsafe(j).mezotipus == mt_elektromos; j++)
                ;
            if (j >= cellak.size())
                throw hiba(1, "electrical joint cell without a following joint electrical/electrothermal cell (%u)", i);

            // beállítjuk az összes elektromos cella indexét i-tõl j-ig, i-t is léptetjük

            for (; i < j; i++) {
                if (cellak.unsafe(i).kapcsolodo_cella_index == 0)
                    throw hiba(1, "non-joint electrical cell among joint cells (%u)", i);
                cellak.unsafe(i).kapcsolodo_cella_index = j;
            }

            // beállítjuk a termikus cella indexét, ha kapcsolódó (elektrotermikus nem lehet kapcsolódó, azaz nála a kapcsolodo_cella_index==0)

            if (cellak.unsafe(i).kapcsolodo_cella_index != 0)
                cellak.unsafe(i).kapcsolodo_cella_index = i - 1;
        }
    }
}


//***********************************************************************
void adat_szimulacio::cellatartomany_feltolto(vektor<klaszter_tartomany>& klaszter_tartomanyok, uns szal_kell) {
// Célszerû lenne összetett cellákat használni, amik akár kompakt modellt
// is tartalmazhatnak, ezzel megoldható lenne a több elektromos egy ter-
// mikusban és hasonló problémák, nem kéne ez a joint cell móka.
//***********************************************************************
    uns start = 0, stop = 0;
    klaszter_tartomany tart;printf("\n\n\nCellaszam = %u\n\n\n", cellak.size() - 1);
    for (uns i = 1; i <= szal_kell; i++) {
        start = stop + 1;
        stop = start + (cellak.size() - start) / (szal_kell - i + 1) - 1; // a cellak_size az utolsónál eggyel nagyobb
        while (stop < cellak.size() && cellak.unsafe(stop).kapcsolodo_cella_index != 0)
            stop++;
        if (i < szal_kell && stop >= cellak.size() - 1)
            throw hiba("adat_szimulacio::cellatartomany_feltolto", "too long joint cell row"); // ez igazából warning is lehetne. Azt jelenti, hogy a joint cellák miatt nem tudunk annyi szálat csinálni, amennyit szeretnénk.
        tart.klaszter_kezdoindex = start;
        tart.klaszter_utolso_index = stop;
        klaszter_tartomanyok.push_back(tart);
    }
}


//***********************************************************************
void adat_szimulacio::fa_rekurziv_klasztertartomany_feltolto(uns kezdoindex, vektor<klaszter_tartomany>& klaszter_tartomanyok, vektor<adat_egyedi> & egyedi_fa_elemek, uns szal_kell) {
//***********************************************************************
    set_hiba_hol h("adat_szimulacio::rekurziv_klasztertartomany_feltolto");

    if (szal_kell == 0)
        throw hiba(1, "szal_kell == 0");
    if (szal_kell == 1) {
        klaszter_tartomany tart;
        tart.klaszter_kezdoindex = fa_elemek[kezdoindex].tartomany_kezdoindex;
        tart.klaszter_utolso_index = kezdoindex;
        klaszter_tartomanyok.push_back(tart);
        return;
    }
    if (fa_elemek[kezdoindex].cella_index != 0)
        throw hiba(1, "unballanced tree");
    uns alszal = (szal_kell + 1) / 2; // így a megengedettnél több szálra bonthatja, de kiegyensúlyozottabb => 2021.03.15.: 32 szál lesz 10 helyett. Ez jó?
    fa_rekurziv_klasztertartomany_feltolto(fa_elemek[kezdoindex].bal_elem_indexe, klaszter_tartomanyok, egyedi_fa_elemek, alszal); // szal_kell - szal_kell / 2
    fa_rekurziv_klasztertartomany_feltolto(fa_elemek[kezdoindex].jobb_elem_indexe, klaszter_tartomanyok, egyedi_fa_elemek, alszal); // szal_kell / 2
    egyedi_fa_elemek.push_back(adat_egyedi{ kezdoindex, szal_kell });
}


//***********************************************************************
void adat_anyag::beolvas_fajlbol(uns hanyas){
//***********************************************************************
    set_hiba_hol h("adat_anyag::beolvas_fajlbol (BM"+uns_to_string(hanyas)+")");

    // BM

    fajl::check_text("BM");
    uns uns_temp = fajl::get_uns("BM");
    if (uns_temp != hanyas)
        throw hiba(1, "BM%u expected, BM%u arrived", hanyas, uns_temp);
    
    // tulajdonsagok

    bool is_vege = false;
    is_fenypor = false;
    do {
        int ch;
        switch (ch = fajl::get_char("material")) {
            case 'F':
                fazisvaltas.beolvas_fajlbol();
                break;
            case 'P':
                switch (ch = fajl::get_char("material P?")) {
                    case 'G':
                        switch (ch = fajl::get_char("material PG?")) {
                            case 'E': elvez.beolvas_fajlbol(fazisvaltas); break;
                            case 'T': thvez.beolvas_fajlbol(fazisvaltas); break;
                            default:  throw hiba(1, "unknown material property (PG%c)", ch);
                        }
                        break;
                    case 'C':
                        switch (ch = fajl::get_char("material PC?")) {
                            case 'T': Cth.beolvas_fajlbol(fazisvaltas); break;
                            default:  throw hiba(1, "unknown material property (PC%c)", ch);
                        }
                        break;
                    case 'S': S.beolvas_fajlbol(fazisvaltas); break;
                    case 'D': D.beolvas_fajlbol(fazisvaltas); break;
                    case 'E': emissivity.beolvas_fajlbol(fazisvaltas); break;
                    default:  throw hiba(1, "unknown material property (P%c)", ch);
                }
                break;
            case 'E':
                fajl::check_text("M");
                uns_temp = fajl::get_uns("EM");
                if (uns_temp != hanyas)
                    throw hiba(1, "EM%u expected, EM%u arrived", hanyas, uns_temp);
                is_vege = true;
                break;
            case 'L':
                fajl::check_text("P");
                is_fenypor = true;
                break;
            default:
                throw hiba(1, "unknown material property (%c)", ch);
        }
    } while (!is_vege);
}


//***********************************************************************
void fazisvalto::beolvas_fajlbol(){
//***********************************************************************
    set_hiba_hol h("fazisvalto::beolvas_fajlbol");
    is = true;
    energia = fajl::get_rvt("material phase change", true);
    int ch = fajl::get_char("material phase change N/H");
    if (ch != 'N' && ch != 'H')
        throw hiba(1, "unknown phase change type (N/H expected, %c arrived)", ch);
    TH1 = fajl::get_rvt("material phase change TH1", true);
    TH2 = fajl::get_rvt("material phase change TH2", true);
    SZ = ch=='H' ? fajl::get_rvt("material phase change SZ", true) : rvt();
}


//***********************************************************************
void tulajdonsag::beolvas_fajlbol(const fazisvalto & fv){
//***********************************************************************
    set_hiba_hol h("tulajdonsag::beolvas_fajlbol");
    int ch;
    bool is_inverz;
    rvt tpars[7];

    // normál / inverz

    switch (ch = fajl::get_char("material property N/I")) {
        case 'N': is_inverz = false; break;
        case 'I': is_inverz = true;  break;
        default:  throw hiba(1, "unknown material property N(ormal)/I(nverz) (%c)", ch);
    }

    // fajta és értékek

    switch (ch = fajl::get_char("material property type (const/exp/...)")) {
        case 'C': 
            set_to_konstans(fajl::get_rvt("material property constant value",true)); 
            break;
        case 'L': 
            for (uns i = 0; i < 2; i++)
                tpars[i] = fajl::get_rvt("material property linear value", true); 
            set_to_lin(tpars[0], tpars[1]); // tpars[0]=b, tpars[1]=m
            break;
        case 'X':
            for (uns i = 0; i < 2; i++)
                tpars[i] = fajl::get_rvt("material property exp value", true);
            set_to_exp(tpars[0], tpars[1]);
            break;
        case 'D':
            for (uns i = 0; i < 3; i++)
                tpars[i] = fajl::get_rvt("material property diode value", true);
            set_to_diode_1(tpars[0], tpars[1], tpars[2]);
            break;
        case 'M':
            for (uns i = 0; i < 7; i++)
                tpars[i] = fajl::get_rvt("material property mizs value", true);
            if (tpars[2] != 2.0)
                throw hiba(1, "mizs function c must be 2, %g found", tpars[2]);
            set_to_mizs(tpars[0], tpars[1], tpars[3], tpars[4], tpars[5], tpars[6]);
            break;
        case 'B': {
                broken_line & bl = set_to_szakaszok();
                bl.beolvas_fajlbol(fv);
            }
            break;
        case 'F': {
                fazis_broken_line & fbl = set_to_fazisvalto();
                fbl.beolvas_fajlbol(fv);
            }
            break;
        default:  throw hiba(1, "unknown material property type (%c) instead of C(onst)/L(inear)/...", ch);
    }
    if (is_inverz)
        set_to_inv();
}

//***********************************************************************
void tulajdonsag::beolvas_fajlbol(){
//***********************************************************************
    set_hiba_hol h("tulajdonsag::beolvas_fajlbol");
    int ch;
    bool is_inverz;
    rvt tpars[7];

    // normál / inverz

    switch (ch = fajl::get_char("junction property N/I")) {
        case 'N': is_inverz = false; break;
        case 'I': is_inverz = true;  break;
        default:  throw hiba(1, "unknown junction property N(ormal)/I(nverz) (%c)", ch);
    }

    // fajta és értékek

    switch (ch = fajl::get_char("junction property type (const/exp/...)")) {
        case 'C': 
            set_to_konstans(fajl::get_rvt("junction property constant value",true)); 
            break;
        case 'L': 
            for (uns i = 0; i < 2; i++)
                tpars[i] = fajl::get_rvt("junction property linear value", true); 
            set_to_lin(tpars[0], tpars[1]);
            break;
        case 'X':
            for (uns i = 0; i < 2; i++)
                tpars[i] = fajl::get_rvt("junction property exp value", true);
            set_to_exp(tpars[0], tpars[1]);
            break;
        case 'D':
            for (uns i = 0; i < 3; i++)
                tpars[i] = fajl::get_rvt("junction property diode value", true);
            set_to_diode_1(tpars[0], tpars[1], tpars[2]);
            break;
        case 'M':
            for (uns i = 0; i < 7; i++)
                tpars[i] = fajl::get_rvt("junction property mizs value", true);
            if (tpars[2] != 2.0)
                throw hiba(1, "mizs function c must be 2, %g found", tpars[2]);
            set_to_mizs(tpars[0], tpars[1], tpars[3], tpars[4], tpars[5], tpars[6]);
            break;
        case 'E':
            for (uns i = 0; i < 6; i++)
                tpars[i] = fajl::get_rvt("junction property erno value", true);
            set_to_erno(tpars[0], tpars[1], tpars[2], tpars[3], tpars[4], tpars[5]);
            break;
        case 'I':
            ill_adat.beolvas_fajlbol();
            ill_adat.set_polinom_by_fitting(poli);
            set_to_polinom();
            break;
        default:  throw hiba(1, "unknown junction property (%c) instead of E(rno)/D(iode)/...", ch);
    }
    if (is_inverz)
        set_to_inv();
}


//***********************************************************************
void adat_junction::beolvas_fajlbol(uns hanyas){
//***********************************************************************
    set_hiba_hol h("adat_junction::beolvas_fajlbol (BJ"+uns_to_string(hanyas)+")");

    // BJ

    fajl::check_text("BJ");
    uns uns_temp = fajl::get_uns("BJ");
    if (uns_temp != hanyas)
        throw hiba(1, "BJ%u expected, BJ%u arrived", hanyas, uns_temp);
    
    // tulajdonsagok

    bool is_vege = false;
    do {
        int ch;
        switch (ch = fajl::get_char("junction")) {
            case 'P':
                switch (ch = fajl::get_char("junction P?")) {
                    case 'G':
                        switch (ch = fajl::get_char("junction PG?")) {
                            case 'J': el_egyenlet.beolvas_fajlbol(); break;
                            default:  throw hiba(1, "unknown junction property (PG%c)", ch);
                        }
                        break;
                    case 'D':
                        switch (ch = fajl::get_char("junction PD?")) {
                            case 'C': D.beolvas_fajlbol(); break;
                            default:  throw hiba(1, "unknown junction property (PD%c)", ch);
                        }
                        break;
                    case 'R':
                        switch (ch = fajl::get_char("junction PR?")) {
                            case 'C': R.beolvas_fajlbol(); break;
                            case 'A': rad.beolvas_fajlbol(); break;
                            default:  throw hiba(1, "unknown junction property (PR%c)", ch);
                        }
                        break;
                    case 'L':
                        switch (ch = fajl::get_char("junction PL?")) {
                            case 'U': lum.beolvas_fajlbol(); break;
                            default:  throw hiba(1, "unknown junction property (PL%c)", ch);
                        }
                        break;
                    case 'A':
                        switch (ch = fajl::get_char("junction PA?")) {
                            case 'S': area = fajl::get_rvt("junction PAS", false); break;
                            default:  throw hiba(1, "unknown junction property (PA%c)", ch);
                        }
                        break;
                    default:  throw hiba(1, "unknown junction property (P%c)", ch);
                }
                break;
            case 'E':
                fajl::check_text("J", "end junction");
                uns_temp = fajl::get_uns("EJ");
                if (uns_temp != hanyas)
                    throw hiba(1, "EJ%u expected, EJ%u arrived", hanyas, uns_temp);
                is_vege = true;
                break;
            default:
                throw hiba(1, "unknown junction property (%c)", ch);
        }
    } while (!is_vege);
}


//***********************************************************************
void adat_analizis_lepes::beolvas_fajlbol(uns hanyas){
//***********************************************************************
    set_hiba_hol h("adat_analizis_lepes::beolvas_fajlbol (BA" + uns_to_string(hanyas) + ")");

    // BA

    fajl::check_text("BA");
    uns uns_temp = fajl::get_uns("BA");
    if (uns_temp != hanyas)
        throw hiba(1, "BA%u expected, BA%u arrived", hanyas, uns_temp);

    // jellemzõk

    bool is_vege = false;
    do {
        int ch;
        adat_anal_beall beall;
        adat_gerj gerj;
        adat_eredm eredm;
        switch (ch = fajl::get_char("analysis")) {
            beall.clear();
            case 'T':
                fajl::check_text("AMB", "Tamb");
                beall.tipus = abt_tamb;
                beall.rvt_ertek = fajl::get_rvt("analysis TAMB", false);
                beallitasok.push_back(beall);
                break;
            case 'C':
                switch (ch = fajl::get_char("analysis change")) {
                    case 'B':
                        switch (ch = fajl::get_char("analysis change")) {
                            case 'N':
                                beall.tipus = abt_change_bn;
                                beall.uns_ertek   = fajl::get_uns("analysis CBN");
                                fajl::check_text("T", "CBNT");
                                beall.uns_ertek_2 = fajl::get_uns("analysis CBNT");
                                beallitasok.push_back(beall);
                                break;
                            case 'V':
                                beall.tipus = abt_change_bv;
                                beall.uns_ertek = fajl::get_uns("analysis CBV");
                                fajl::check_text("T", "CBVT");
                                beall.rvt_ertek = fajl::get_rvt("analysis CBVT", false);
                                switch (ch = fajl::get_char_or_eol("analysis CBVT")) {
                                    case 0:
                                        beallitasok.push_back(beall);
                                        break;
                                    case ';':
                                        beall.rvt_ertek_2 = fajl::get_rvt("analysis CBVT2", false);
                                        beallitasok.push_back(beall);
                                        break;
                                    default:  throw hiba(1, "unknown analysis change property last character (CBV%uT%g%c)", beall.uns_ertek, beall.rvt_ertek, ch);
                                }
                                break;
                            default:  throw hiba(1, "unknown analysis change property (CB%c)", ch);
                        }
                        break;
                    case 'M':
                        fajl::check_text("N", "CMN");
                        beall.tipus = abt_change_mn;
                        beall.uns_ertek   = fajl::get_uns("analysis CMN");
                        fajl::check_text("T", "CMNT");
                        beall.uns_ertek_2 = fajl::get_uns("analysis CMNT");
                        beallitasok.push_back(beall);
                        break;
                    case 'J':
                        fajl::check_text("N", "CJN");
                        beall.tipus = abt_change_jn;
                        beall.uns_ertek = fajl::get_uns("analysis CJN");
                        fajl::check_text("T", "CJNT");
                        beall.uns_ertek_2 = fajl::get_uns("analysis CJNT");
                        beallitasok.push_back(beall);
                        break;
                    case 'T':
                        beall.tipus = abt_ct;
                        switch (ch = fajl::get_char("analysis CT?")) {
                            case 'D': beall.uns_ertek = 1; break;
                            case 'F': beall.uns_ertek = 2; break;
                            case 'E': beall.uns_ertek = 3; break;
                            default:  throw hiba(1, "unknown analysis property (CT%c)", ch);
                        }
                        beallitasok.push_back(beall);
                        break;
                    default:  throw hiba(1, "unknown analysis change property (C%c)", ch);
                }
                break;
            case 'P':
                switch (ch = fajl::get_char("analysis P?")) {
                    case 'C':
                        beall.tipus = abt_I0;
                        beall.rvt_ertek = fajl::get_rvt("analysis PC", false);
                        break;
                    case 'E':
                        beall.tipus = abt_max_error;
                        beall.rvt_ertek = fajl::get_rvt("analysis PE", false);
                        break;
                    case 'I':
                        beall.tipus = abt_max_iter;
                        beall.uns_ertek = fajl::get_uns("analysis PI");
                        break;
                    default:  throw hiba(1, "unknown analysis property (P%c)", ch);
                }
                beallitasok.push_back(beall);
                break;
            case 'A':
                switch (ch = fajl::get_char("analysis P?")) {
                    case 'D':
                        tipus = alt_dc;
                        break;
                    case 'S':
                        tipus = alt_trans;
                        value = fajl::get_rvt("analysis step time", false);
                        break;
                    case 'A':
                        tipus = alt_ac;
                        value = fajl::get_rvt("analysis ac frequency", false);
                        break;
                    default:  throw hiba(1, "unknown analysis type (A%c)", ch);
                }
                break;
            case 'D':
                switch (ch = fajl::get_char("analysis D?")) {
                    case 'A': beall.tipus = abt_del_all_prev; break;
                    case 'X': beall.tipus = abt_del_all_gerj; break;
                    case 'T': beall.tipus = abt_del_fa;       break;
                    default:  throw hiba(1, "unknown analysis property (D%c)", ch);
                }
                beallitasok.push_back(beall);
                break;
            case 'X':
                gerj.struktura_index = fajl::get_uns("analysis excitation structure index");
                switch (ch = fajl::get_char("analysis X ?")) {
                    case 'E':
                        switch (ch = fajl::get_char("analysis P?")) {
                            case 'U':
                                gerj.tipus = gt_U;
                                gerj.ertek = fajl::get_rvt("analysis excitation EU value", false);
                                gerjesztesek.push_back(gerj);
                                break;
                            case 'I':
                                gerj.tipus = gt_I;
                                gerj.ertek = fajl::get_rvt("analysis excitation EI value", false);
                                gerjesztesek.push_back(gerj);
                                break;
                            case 'N':
                                gerj.tipus = gt_el_none;
                                gerj.ertek = rvt();
                                gerjesztesek.push_back(gerj);
                                break;
                            default:  throw hiba(1, "unknown analysis excitation type (E%c)", ch);
                        }
                        break;
                    case 'T':
                        switch (ch = fajl::get_char("analysis P?")) {
                            case 'T':
                                gerj.tipus = gt_T;
                                gerj.ertek = fajl::get_rvt("analysis excitation TT value", false);
                                gerjesztesek.push_back(gerj);
                                break;
                            case 'P':
                                gerj.tipus = gt_P;
                                gerj.ertek = fajl::get_rvt("analysis excitation TP value", false);
                                gerjesztesek.push_back(gerj);
                                break;
                            case 'N':
                                gerj.tipus = gt_th_none;
                                gerj.ertek = rvt();
                                gerjesztesek.push_back(gerj);
                                break;
                            default:  throw hiba(1, "unknown analysis excitation type (T%c)", ch);
                        }
                        break;
                    default:  throw hiba(1, "unknown analysis property (D%c)", ch);
                }
                break;
            case 'R':
                switch (ch = fajl::get_char("analysis R?")) {
                    case 'P':
                        eredm.cella_index = fajl::get_uns("analysis RP cella index");
                        switch (ch = fajl::get_char("analysis RP data type (U/T/F)")) {
                            case 'C':
                                eredm.tipus = et_c_pontprobe;
                                eredm.mit_ment = mm_UT;
                                eredmenyek.push_back(eredm);
                                break;
                            case 'F':
                                eredm.tipus = et_f_pontprobe;
                                eredm.mit_ment = mm_UT;
                                eredm.face_index = fajl::get_uns("analysis RP face index");
                                eredmenyek.push_back(eredm);
                                break;
                            default:  throw hiba(1, "unknown analysis RP data type (RP%c)", ch);
                        }
                        break;
                    case 'C':
                        eredm.tipus = et_currentprobe;
                        eredm.mit_ment = mm_IP;
                        eredm.currentProbe.set_size(fajl::get_uns("analysis RC face number"));
                        fajl::check_text(";", "analysis RC face number");
                        for (uns i = 0; i < eredm.currentProbe.size(); i++) {
                            eredm.currentProbe[i].cella_index = fajl::get_uns("analysis RC cell index");
                            fajl::check_text("F", "analysis RC cell and index separator");
                            eredm.currentProbe[i].face_index = fajl::get_uns("analysis RC face index");
                            fajl::check_text(";", "analysis RC face index");
                        }
                        eredmenyek.push_back(eredm);
                        break;
                    case 'M':
                        ch = fajl::get_char("analysis RP where type (C/F)");
                        if (ch == 'C')
                            eredm.tipus = et_c_map;
                        else if (ch == 'F')
                            eredm.tipus = et_f_map;
                        else throw hiba(1, "unknown analysis RP where type (RM%c)", ch);
                        switch (ch = fajl::get_char("analysis RP data type (V/I/R)")) {
                            case 'V':
                                eredm.mit_ment = mm_UT;
                                eredmenyek.push_back(eredm);
                                break;
                            case 'I':
                                eredm.mit_ment = mm_IP;
                                eredmenyek.push_back(eredm);
                                break;
                            case 'R':
                                eredm.mit_ment = mm_RL;
                                eredmenyek.push_back(eredm);
                                break;
                            default:  throw hiba(1, "unknown analysis RP data type (RM%c%c)", eredm.tipus == et_c_map ? 'C' : 'F', ch);
                        }
                        break;
                    default:  throw hiba(1, "unknown analysis excitation type (T%c)", ch);
                }
                break;
            case 'E':
                fajl::check_text("A", "end analysis");
                uns_temp = fajl::get_uns("EA");
                if (uns_temp != hanyas)
                    throw hiba(1, "EA%u expected, EA%u arrived", hanyas, uns_temp);
                is_vege = true;
                break;
            default:
                throw hiba(1, "unknown analysis property (%c)", ch);
        }
    } while (!is_vege);
}


//***********************************************************************
void adat_cella::beolvas_fajlbol(uns hanyas){
//***********************************************************************
    set_hiba_hol h("adat_cella::beolvas_fajlbol (BL" + uns_to_string(hanyas) + ")");

    // BL

    fajl::check_text("BL");
    cella_tipus = ct_faces_cella; // !!!
    uns uns_temp = fajl::get_uns("BL");
    if (uns_temp != hanyas)
        throw hiba(1, "BL%u expected, BL%u arrived", hanyas, uns_temp);

    // jellemzõk

    fajl::check_text("C");
    struktura_index = fajl::get_uns("BL Cx");
    fajl::check_text("M");
    anyag_index = fajl::get_uns("BL Mx");
    fajl::check_text("V");
    volume_index = fajl::get_uns("BL Vx");
    int ch = fajl::get_char("cell field type");
    if (ch == 'F') {
        switch (ch = fajl::get_char("cell field type")) {
        case 'E':
            //ch = fajl::get_char_or_eol("cell field type");
            switch (ch = fajl::get_char_or_eol("cell field type")) {
            case 0:
                mezotipus = mt_elektromos;
                break;
            case 'T':
                mezotipus = mt_elektrotermikus;
                break;
            default:  throw hiba(1, "unknown field type (FE%c)", ch);
            }
            //tipus = mt_elektromos;
            break;
        case 'T':
            mezotipus = mt_termikus;
            break;
        default:  throw hiba(1, "unknown field type (F%c)", ch);
        }
        if ((ch = fajl::get_char("joint cell (K) or noise (Z) or cell face count (NF)")) == 'K') {
            if (mezotipus == mt_elektromos || mezotipus == mt_termikus)
                kapcsolodo_cella_index = 1; // a kapcsolodo_cellak_beallitasa()-ban állítjuk be a ténylegeset
            else throw hiba(1, "electrothermal cell (%u) is defined as joint", hanyas);
            ch = fajl::get_char("noise (Z) or cell face count (NF)");
        }

        // center face zaj

        if (ch == 'Z') {
            zaj_index = fajl::get_uns("noise index");
            ch = fajl::get_char("cell face count (NF)");
        }

        // face-ek

        if (ch != 'N')
            throw hiba(1, "cell face count (N) expected, %c arrived", ch);
        fajl::check_text("F", "cell face count (NF)");
        facek.set_size(1 + fajl::get_uns("NF face count"));
    }
    else if (ch == 'T') {
        uns tipus = fajl::get_uns("cell type");
        if (tipus == 1) {
            cella_tipus = ct_th6;
            facek.set_size(1 + 6);
        }
        else {
            throw hiba(1, "unknown cell type (%u)", tipus);
        }
    }
    else
        throw hiba(1, "cell type (T) or field type (F) expected, %c arrived", ch);
    for (uns i = 1; i < facek.size(); i++)
        facek[i].beolvas_fajlbol(i, *this);

    // L (besugárzó cella)

    while ((ch = fajl::get_char("radiating cell (L), light path (NP) or end (EL)")) == 'L') {
        if (mezotipus == mt_elektromos)
            throw hiba(1, "electrical cell (%u) is defined as drain of radiation", hanyas);
        adat_besugarzo_cella be;
        be.cella_index = fajl::get_uns("radiating cell index");
        fazisvalto fv;
        be.arany.beolvas_fajlbol(fv);
        besugarzo_cella.push_back(be);
    }

    // EL

    fajl::check_text("L");
    uns_temp = fajl::get_uns("EL");
    if (uns_temp != hanyas)
        throw hiba(1, "EL%u expected, EL%u arrived", hanyas, uns_temp);
}


//***********************************************************************
void adat_face::beolvas_fajlbol(uns hanyas, const adat_cella & cella){
//***********************************************************************
    set_hiba_hol h("adat_face::beolvas_fajlbol");

    // típus

    int ch;
    switch (ch = fajl::get_char("face")) {
        case 'J':
            tipus = ft_csatlakozo;
            is_perem = false;
            is_kulso = true;
            break;
        case 'P':
            tipus = ft_normal_perem;
            is_kulso = false;
            is_perem = true;
            break;
        case 'C':
            tipus = ft_centroid;
            is_perem = false;
            is_kulso = true;
            break;
        case 'S':
            tipus = ft_spec_perem;
            is_kulso = false;
            is_perem = true;
            break;
        default:  throw hiba(1, "unknown face type (%c)", ch);
    }

    // index

    uns uns_temp = fajl::get_uns("face index");
    if (uns_temp != hanyas)
        throw hiba(1, "%c%u expected, %c%u arrived", ch, hanyas, ch, uns_temp);

    // el/th

    switch (ch = fajl::get_char("face E/T")) {
        case 'E':
            is_el = true;
            break;
        case 'T':
            is_el = false;
            break;
        default:  throw hiba(1, "unknown face type, E/T expected, %c arrived", ch);
    }

    // csatlakozó index

    csatlakozo_index = fajl::get_uns("connetion index");

    if (tipus == ft_centroid)
        goto sor_vege;

    // dual

    if ((ch = fajl::get_char("face noise/material/junction/area/dual")) == 'D') {
        parja = hanyas - 1;
        anyag_index = cella.facek[hanyas - 1].anyag_index;
        A_index = cella.facek[hanyas - 1].A_index;
        L_index = cella.facek[hanyas - 1].L_index;
        goto sor_vege;
    }

    // zaj

    if (ch == 'Z') {
        anyag_index = fajl::get_uns("noise index");
        ch = fajl::get_char("face material/junction/area/dual");
    }

    // anyag

    if (ch == 'M') {
        anyag_index = fajl::get_uns("face material");
        ch = fajl::get_char("cell face junction/area");
    }
    else {
        anyag_index = cella.anyag_index;
    }

    // junction

    if (ch == 'J') {
        junction_index = fajl::get_uns("face junction");
        ch = fajl::get_char("cell face area");
    }

    // area index

    switch (ch) {
        case 'A': oldal = pit_none; break;
        case 'W': oldal = pit_west; break;
        case 'E': oldal = pit_east; break;
        case 'S': oldal = pit_south; break;
        case 'N': oldal = pit_north; break;
        case 'B': oldal = pit_bottom; break;
        case 'T': oldal = pit_top; break;
        default: throw hiba(1, "cell face area (AWESNBT) expected, %c arrived", ch);
    }
    A_index = fajl::get_uns("face area index");

    // lenght index

    if ((ch = fajl::get_char("face length index")) == 'L') {
        L_index = fajl::get_uns("face length index");
    }
    else throw hiba(1, "cell face length (L) expected, %c arrived", ch);

    // duális

    switch (ch = fajl::get_char("other field face type")) {
        case 'X':
            parja = 0;
            break;
        case 'F':
            parja = fajl::get_uns("the other field face index");
            break;
        default:  throw hiba(1, "unknown other field face type (%c), F/X expected", ch);
    }

    // peremindex

    switch (ch = fajl::get_char_or_eol("face boundary condition index or EOL")) {
        case 0:
            goto sor_vege;
            break;
        case 'I':
            perem_index = fajl::get_uns("face boundary condition index");
            break;
        default:  throw hiba(1, "unknown face property after L index (%c), I was expected", ch);
    }

    // speciális perem koordináták

    switch (ch = fajl::get_char_or_eol("face special boundary condition X or EOL")) {
        case 0:
            goto sor_vege;
            break;
        case 'X':
            perem_x = fajl::get_uns("face special boundary condition X");
            fajl::check_text("Y", "face special boundary condition Y");
            perem_y = fajl::get_uns("face special boundary condition Y");
            switch (ch = fajl::get_char("face special boundary condition direction")) {
                case 'E':
                    peremirany = pit_east;
                    break;
                case 'W':
                    peremirany = pit_west;
                    break;
                case 'S':
                    peremirany = pit_south;
                    break;
                case 'N':
                    peremirany = pit_north;
                    break;
                case 'B':
                    peremirany = pit_bottom;
                    break;
                case 'T':
                    peremirany = pit_top;
                    break;
                default:  throw hiba(1, "unknown face special boundary condition direction (%c), WESNBT was expected", ch);
            }
            break;
        default:  throw hiba(1, "unknown face property after L index (%c), I was expected", ch);
    }

sor_vege:
    return;
}


//***********************************************************************
void adat_fa_elem::beolvas_fajlbol(uns hanyas){
//***********************************************************************
    set_hiba_hol h("adat_fa_elem::beolvas_fajlbol (BF" + uns_to_string(hanyas) + ")");

    // BL

    fajl::check_text("F");
    uns uns_temp = fajl::get_uns("BF");
    if (uns_temp != hanyas)
        throw hiba(1, "BF%u expected, BF%u arrived", hanyas, uns_temp);

    // jellemzõk

    fajl::check_text("F","BF F");
    tartomany_kezdoindex = fajl::get_uns("BF Fx");
    fajl::check_text("L", "BF L");
    bal_elem_indexe = fajl::get_uns("BF Lx");
    fajl::check_text("R", "BF R");
    jobb_elem_indexe = fajl::get_uns("BF Rx");

    // Másolandó blokkok

    fajl::check_text("NM");
    masolando.set_size(fajl::get_uns("NM from-to-n blocks"));
    fajl::check_text("A");
    A = fajl::get_uns("MN A");
    fajl::check_text("B");
    B = fajl::get_uns("MN B");
    adat_masolando akt;
    int ch;
    for (uns i = 0; i < masolando.size(); i++) {
        switch (ch = fajl::get_char("from-to-n block type")) {
            case 'A':
                masolando[i].is_masolando = true;
                masolando[i].is_centroid = false;
                masolando[i].hova = fajl::get_uns("to");
                switch (ch = fajl::get_char("from-to-n block A source")) {
                    case 'E':
                        masolando[i].is_balbol = true;
                        masolando[i].bal_from = fajl::get_uns("from");
                        break;
                    case 'K':
                        masolando[i].is_balbol = false;
                        masolando[i].jobb_from = fajl::get_uns("from");
                        break;
                    default:  throw hiba(1, "unknown from-to-n block A source type (%c), E/K was expected", ch);
                }
                fajl::check_text("N");
                masolando[i].db = fajl::get_uns("n");
                break;
            case 'B':
                masolando[i].is_masolando = false;
                masolando[i].is_centroid = false;
                masolando[i].is_balbol = false;
                masolando[i].hova = fajl::get_uns("to");
                fajl::check_text("E");
                masolando[i].bal_from = fajl::get_uns("from left");
                fajl::check_text("K");
                masolando[i].jobb_from = fajl::get_uns("from right");
                fajl::check_text("N");
                masolando[i].db = fajl::get_uns("n");
                break;
            case 'C':
                masolando[i].is_masolando = false;
                masolando[i].is_centroid = true;
                masolando[i].is_balbol = false;
                masolando[i].hova = fajl::get_uns("to");
                fajl::check_text("E");
                masolando[i].bal_from = fajl::get_uns("from left");
                fajl::check_text("K");
                masolando[i].jobb_from = fajl::get_uns("from right");
                masolando[i].db = 1;
                break;
            default:  throw hiba(1, "unknown face special boundary condition direction (%c), WESNBT was expected", ch);
        }
    }
        
    // EF

    fajl::check_text("EF");
    uns_temp = fajl::get_uns("EF");
    if (uns_temp != hanyas)
        throw hiba(1, "EF%u expected, EF%u arrived", hanyas, uns_temp);
}

}