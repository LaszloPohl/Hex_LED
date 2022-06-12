//***********************************************************************
// párhuzamos feldolgozó header
// Creation date:  2018. 08. 09.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef PARHUZAMOS_FELDOLGOZO_HEADER
#define	PARHUZAMOS_FELDOLGOZO_HEADER
//***********************************************************************


//***********************************************************************
#include "kozos.h"
#include <thread>
#include <mutex>
#include "vektor.hpp"
#include "matrix.hpp"
#include "bemenet.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
class vezerlo;
struct adat_eredm;
struct mentendo_adatok;
class Parhuzamos_Feldolgozo;
//***********************************************************************


//***********************************************************************
struct sumhiba_tipus {
//***********************************************************************
    //***********************************************************************
    rvt max_UT_hiba, max_U_hiba, max_T_hiba;
    rvt sum_IP_hiba;
    rvt max_abs_Uc, max_abs_Tc;
    //***********************************************************************
    sumhiba_tipus() :max_UT_hiba{ 0 }, max_U_hiba{ 0 }, max_T_hiba{ 0 },
        sum_IP_hiba{ 0 }, max_abs_Uc{ 0 }, max_abs_Tc{ 0 } {}
    //***********************************************************************
    void zero() {
    //***********************************************************************
        max_UT_hiba = max_U_hiba = max_T_hiba = sum_IP_hiba = max_abs_Uc = max_abs_Tc = rvt();
    }
    //***********************************************************************
    void to_sum(sumhiba_tipus & sum) {
    //***********************************************************************
        if (max_UT_hiba > sum.max_UT_hiba) 
            sum.max_UT_hiba = max_UT_hiba;
        if (max_U_hiba > sum.max_U_hiba)
            sum.max_U_hiba = max_U_hiba;
        if (max_T_hiba > sum.max_T_hiba)
            sum.max_T_hiba = max_T_hiba;
        sum.sum_IP_hiba += sum_IP_hiba;
        if (max_abs_Uc > sum.max_abs_Uc)
            sum.max_abs_Uc = max_abs_Uc;
        if (max_abs_Tc > sum.max_abs_Tc)
            sum.max_abs_Tc = max_abs_Tc;
    }
};


//***********************************************************************
class szal_feladat_adatai {
//***********************************************************************

    //***********************************************************************
    bool is_fw_egyedi_indithato() const;
    bool is_bw_indithato() const;
    //***********************************************************************

public:
    //***********************************************************************
    enum allapot { a_ures, a_feldolgozasra_var, a_feldolgozas_alatt, a_kesz };
    allapot all;
    bool is_felteteles, is_torolhetok_a_szalak;
    szal_tipus tipus;   // szt_nem_csinal_semmit, szt_cellafeldolgozo_klaszter, szt_fw_klaszter, szt_fw_egyedi, szt_bw_egyedi, szt_bw_klaszter, 
                        // szt_matrix, szt_belso_mentes, szt_fajlba_mentes, szt_N
    uns thread_id; // A Para_engine::betesz függvény állítja be
    // az összes tartomány típus futtatásához szükséges adatok
    uns klaszter_kezdoindex, klaszter_utolso_index; // egyedinél a klaszter_kezdoindexet használjuk, az adat_szimulacio::fa_elemek-en belüli index, nem a dcfa_d/dcfa_f/acfa_d/acfa_f indexe !
    uns faelem_kezdoindex, faelem_utolso_index; // az elem indexe a dcfa_d/dcfa_f/acfa_d/acfa_f-n belül
    uns al_szalszam; // ennyi szálra osztható a redukció
    bool is_ac, is_float; // Az aktuális fa típusa
    uns melyik_fa; // 1 vagy 2
    // cella klaszternél
    szal_altipus altipus;
    rvt atlag_IP_hiba, max_Uc, max_Tc;
    // belsõ mentéshez
    rvt Tamb;
    const vektor<adat_eredm> * p_akt_eredm; // const vektor<adat_eredm> *, mik a kiírandó eredmények.
    // belsõ mentéshez és fájlba mentéshez
    ::std::string eredm_utvonal;
    lista<mentendo_adatok> * p_mentendo_adatok;  // a belsõ mentés hozza létre, õ passzolja a fájlba mentésnek
    analizis_lepes_tipus analizis_tipus;
    rvt analizis_value; // idõ vagy frekvencia
    uns akt_anal_index;
    bool is_use_commas;
    // mátrixmûveletek
    enum matrix_muvelet_tipus{ mmt_none, mmt_dc_f, mmt_dc_d, mmt_ac_f, mmt_ac_d };
    matrix_muvelet_tipus matrix_muvelet;
    uns matrix_fv_tipus;
    hivando_matrixfuggveny< float > matrix_adat_dc_f;
    hivando_matrixfuggveny< double > matrix_adat_dc_d;
    hivando_matrixfuggveny< ::std::complex<float> > matrix_adat_ac_f;
    hivando_matrixfuggveny< ::std::complex<double> > matrix_adat_ac_d;
    ::std::condition_variable * p_kesz_egy_feladat;
    std::atomic<bool> * p_is_kesz;
    // visszaadott adatok
    sumhiba_tipus sumhiba;
    //***********************************************************************

    //***********************************************************************
    szal_feladat_adatai() :all{ a_ures }, is_felteteles{ false }, is_torolhetok_a_szalak{ false }, tipus{ szt_nem_csinal_semmit },
        klaszter_kezdoindex{ 0 }, klaszter_utolso_index{ 0 }, faelem_kezdoindex{ 0 }, faelem_utolso_index{ 0 }, al_szalszam{ 1 },
        altipus{ szat_semmi }, p_akt_eredm{ nullptr }, p_mentendo_adatok{ nullptr }, is_use_commas{ false }, akt_anal_index{ 0 },
        atlag_IP_hiba{ rvt() }, max_Uc{ rvt() }, max_Tc{ rvt() }, p_kesz_egy_feladat{ nullptr }, p_is_kesz{ nullptr }, 
        matrix_muvelet{ mmt_none }, thread_id{ 0 }, matrix_fv_tipus{ 0 } {}
    //***********************************************************************
    bool is_indithato() const; // lekérdezi a bemenõ faelemeket, hogy készen vannak-e, azaz elindítható-e ez a feladat
    //***********************************************************************
    
    //***********************************************************************
    void set_fa_kezdoindex_es_utolso_index(uns adat_kezdoi, uns adat_utolsoi, uns eltolas, uns sub_threads) {
    //***********************************************************************
        klaszter_kezdoindex = adat_kezdoi; klaszter_utolso_index = adat_utolsoi; faelem_kezdoindex = adat_kezdoi - eltolas; faelem_utolso_index = adat_utolsoi - eltolas; al_szalszam = sub_threads;
    }
    //***********************************************************************
    void set_cellatomb_kezdoindex_es_utolso_index(uns adat_kezdoi, uns adat_utolsoi) {
    //***********************************************************************
        klaszter_kezdoindex = adat_kezdoi; klaszter_utolso_index = adat_utolsoi;
    }
    //***********************************************************************
    void beallitas_szaltorleshez(szal_tipus tip) { all = a_feldolgozasra_var; is_felteteles = false;  is_torolhetok_a_szalak = true;  tipus = tip; }
    void beallitas_fw_klaszternek() { all = a_feldolgozasra_var; is_felteteles = false; is_torolhetok_a_szalak = false; tipus = szt_fw_klaszter; }
    void beallitas_fw_egyedinek()   { all = a_feldolgozasra_var; is_felteteles = true;  is_torolhetok_a_szalak = false; tipus = szt_fw_egyedi; }
    void beallitas_bw_egyedinek()   { all = a_feldolgozasra_var; is_felteteles = true;  is_torolhetok_a_szalak = false; tipus = szt_bw_egyedi;   }
    void beallitas_bw_klaszternek() { all = a_feldolgozasra_var; is_felteteles = true;  is_torolhetok_a_szalak = false; tipus = szt_bw_klaszter; }
    void beallitas_lepo_klaszternek(szal_altipus szat) {all = a_feldolgozasra_var; is_felteteles = false; is_torolhetok_a_szalak = false; tipus = szt_lepo_klaszter; altipus = szat; }
    void beallitas_cellafeldolgozo_klaszternek_pre()  { all = a_feldolgozasra_var; is_felteteles = false; is_torolhetok_a_szalak = false; tipus = szt_cellafeldolgozo_klaszter; sumhiba.zero(); altipus = szat_pre; }
    void beallitas_cellafeldolgozo_klaszternek_post() { all = a_feldolgozasra_var; is_felteteles = false; is_torolhetok_a_szalak = false; tipus = szt_cellafeldolgozo_klaszter; sumhiba.zero(); altipus = szat_post; }
    void beallitas_belso_mentesnek() { all = a_feldolgozasra_var;  is_felteteles = false; is_torolhetok_a_szalak = false; tipus = szt_belso_mentes; analizis_value = rvt(); }
    void beallitas_fajlba_mentesnek() { all = a_feldolgozasra_var; is_felteteles = false; is_torolhetok_a_szalak = false; tipus = szt_fajlba_mentes; analizis_value = rvt(); }
    //***********************************************************************

    //***********************************************************************
    void beallitas_matrixmuveletnek(hivando_matrixfuggveny< float > & fv) {
    //***********************************************************************
        all = a_feldolgozasra_var; 
        is_felteteles = false; 
        is_torolhetok_a_szalak = false; 
        tipus = szt_matrix; 
        matrix_muvelet = mmt_dc_f;
        matrix_adat_dc_f = fv;
        p_kesz_egy_feladat = fv.p_kesz_egy_feladat;
        p_is_kesz = fv.p_is_kesz;
        matrix_fv_tipus = (uns)fv.tipus;
    }

    //***********************************************************************
    void beallitas_matrixmuveletnek(hivando_matrixfuggveny< double > & fv) {
    //***********************************************************************
        all = a_feldolgozasra_var; 
        is_felteteles = false; 
        is_torolhetok_a_szalak = false; 
        tipus = szt_matrix; 
        matrix_muvelet = mmt_dc_d;
        matrix_adat_dc_d = fv;
        p_kesz_egy_feladat = fv.p_kesz_egy_feladat;
        p_is_kesz = fv.p_is_kesz;
        matrix_fv_tipus = (uns)fv.tipus;
    }

    //***********************************************************************
    void beallitas_matrixmuveletnek(hivando_matrixfuggveny< ::std::complex<float> > & fv) {
    //***********************************************************************
        all = a_feldolgozasra_var; 
        is_felteteles = false; 
        is_torolhetok_a_szalak = false; 
        tipus = szt_matrix; 
        matrix_muvelet = mmt_ac_f;
        matrix_adat_ac_f = fv;
        p_kesz_egy_feladat = fv.p_kesz_egy_feladat;
        p_is_kesz = fv.p_is_kesz;
        matrix_fv_tipus = (uns)fv.tipus;
    }

    //***********************************************************************
    void beallitas_matrixmuveletnek(hivando_matrixfuggveny< ::std::complex<double> > & fv) {
    //***********************************************************************
        all = a_feldolgozasra_var; 
        is_felteteles = false; 
        is_torolhetok_a_szalak = false; 
        tipus = szt_matrix; 
        matrix_muvelet = mmt_ac_d;
        matrix_adat_ac_d = fv;
        p_kesz_egy_feladat = fv.p_kesz_egy_feladat;
        p_is_kesz = fv.p_is_kesz;
        matrix_fv_tipus = (uns)fv.tipus;
    }
};


//***********************************************************************
void szalfuttato_fuggveny(szal_tipus mit_futtat);
//***********************************************************************


//***********************************************************************
extern vektor<rvt> csatlakozo_aramok_dc; // a 0 indexû dummy!
extern vektor<iter_csomopont> iter_csomopontok_dc; // a 0 indexû dummy!
//***********************************************************************


//***********************************************************************
class Para_engine {
//***********************************************************************
    //***********************************************************************
    struct list_elem {
    //***********************************************************************
        szal_feladat_adatai adatok;
        list_elem *next;
        list_elem(const szal_feladat_adatai & adatok) :adatok{ adatok }, next{ nullptr } {}
        list_elem(const szal_feladat_adatai & adatok, uns thread_id) :adatok{ adatok }, next{ nullptr } { this->adatok.thread_id = thread_id; }
    };
    //***********************************************************************
    class para_list {
    //***********************************************************************
        list_elem *first, *last;
        uns akt_db; // hány elemû a lista
        //***********************************************************************
        list_elem * get_indithato() {
        //***********************************************************************
            for (list_elem *p = first; p != nullptr; p = p->next)
                if (p->adatok.is_indithato()) 
                    return p;
            return nullptr;
        }
    public:
        //***********************************************************************
        para_list() :first{ nullptr }, last{ nullptr }, akt_db{ 0 } {}
        //***********************************************************************
        void push_back(list_elem * elem) { // elem->next nullptr kell legyen
        //***********************************************************************
            elem->next = nullptr;
            if (first == nullptr) {
                first = last = elem;
                akt_db = 1;
            }
            else {
                last->next = elem;
                last = elem;
                akt_db++;
            }
        }
        //***********************************************************************
        void push_back(const szal_feladat_adatai & adatok) {
        //***********************************************************************
            push_back(new list_elem{ adatok });
        }
        //***********************************************************************
        void push_back(const szal_feladat_adatai & adatok, uns thread_id) {
        //***********************************************************************
            push_back(new list_elem{ adatok, thread_id });
        }
        //***********************************************************************
        list_elem * pop_front() {
        //***********************************************************************
            if (first == nullptr)
                return nullptr; // akt_db nem változik, marad 0
            list_elem * p = first;
            first = first->next;
            if (first == nullptr)
                last = nullptr;
            akt_db--;
            return p;
        }
        //***********************************************************************
        list_elem * pop_elem(const szal_feladat_adatai * p_feladat) {
        //***********************************************************************
            list_elem *p, *prev;
            for (p = first, prev = nullptr; p != nullptr; prev = p, p = p->next)
                if (&p->adatok == p_feladat) {
                    if (prev == nullptr) { // first
                        if (first == last)
                            last = nullptr;
                        first = first->next; // = p->next
                    }
                    else {
                        if (p == last) {
                            last = prev;
                        }
                        prev->next = p->next;
                    }
                    akt_db--;
                    return p;
                }
            return p; // nullptr
        }
        //***********************************************************************
        list_elem * pop_indithato() {
        //***********************************************************************
            list_elem * p = get_indithato();
            if (p == nullptr)
                return nullptr;
            p = pop_elem(&p->adatok);
            return p;
        }
        //***********************************************************************
        uns getElementNum() const { return akt_db; }
        //***********************************************************************
        bool is_indithato() { return get_indithato() != nullptr; }
        //***********************************************************************
        sumhiba_tipus get_sumhiba() const {
        //***********************************************************************
            sumhiba_tipus sumhiba;
            for (list_elem *p = first; p != nullptr; p = p->next) {
                p->adatok.sumhiba.to_sum(sumhiba);
            }
            return sumhiba;
        }
    };
    //***********************************************************************
    struct adott_tipusu_szalak {
    //***********************************************************************
        para_list waiting_list, waiting_priority_list; // indításra várakozó feladatok, száltípus szerint, a priority elemek elõbb indulnak
        para_list in_progress_list;  // futó feladatok, száltípus szerint
        para_list finished_list;     // kész feladatok, száltípus szerint: ide csak akkor kerül, ha nem azonnal törölhetõ
        vektor<::std::thread> szalak;
        ::std::condition_variable van_varakozo, van_futo;
        bool is_felteteles, is_azonnal_torolheto;
        szal_tipus ertesitendo_felteteles;
        adott_tipusu_szalak() :is_felteteles{ false }, is_azonnal_torolheto{ false }, ertesitendo_felteteles{ szt_N } {}
        void print_status(uns i)const { 
            if(waiting_list.getElementNum() + waiting_priority_list.getElementNum() + in_progress_list.getElementNum() + finished_list.getElementNum() > 0)
                printf("%u: var: %u + %u, fut: %u, kesz: %u; ", i, waiting_list.getElementNum(), waiting_priority_list.getElementNum(), in_progress_list.getElementNum(), finished_list.getElementNum());
        }
    };
    //***********************************************************************
    struct thread_log_elem {
    //***********************************************************************
        uns thread_id;
        uns in_progress_threads_db, free_threads_db;
        szal_tipus tipus;
        enum thread_log_mit_tortenik { tlmt_most_indul, tlmt_most_van_vege };
        thread_log_mit_tortenik mi_tortenik;
        uns matrix_fuggveny_tipus; // ide betesszük az enum uns értékét
        std::chrono::time_point<std::chrono::system_clock> ido;
        double id_time;
        thread_log_elem() :id_time{ 0 } {}
    };
    //***********************************************************************
    uns osszes_nincs_kesz_db;
    uns akt_thread_id; // A betesz függvény beteszi a szál adatába, és növeli az értékét
    uns in_progress_threads_db; // Ennyi szál aktív jelenleg, azaz ennyi van valamelyik in_progress_list-ben
    bool is_log_szalak;
    lista<thread_log_elem> thread_log_lista;
    std::chrono::time_point<std::chrono::system_clock> indulasi_ido;
    ::std::mutex lezaro_mutex;
    adott_tipusu_szalak szalak[szt_N];
    ::std::condition_variable van_futo;
    //***********************************************************************
    void kesz_egy_munka(szal_feladat_adatai * p_feladat);
    szal_feladat_adatai * kivesz(szal_tipus tipus);
    friend void szalfuttato_fuggveny(szal_tipus mit_futtat);
    //***********************************************************************
public:
    //***********************************************************************
    Para_engine() :osszes_nincs_kesz_db{ 0 }, akt_thread_id{ 1000 }, in_progress_threads_db{ 0 }, is_log_szalak{ false }
        { indulasi_ido = std::chrono::system_clock::now(); }
    //***********************************************************************
    void print_status()const { printf("nincs_kesz: %u, ", osszes_nincs_kesz_db); for (uns i = 1; i < szt_N; i++) szalak[i].print_status(i); printf("\n"); }
    void set_szalszam(uns szal_db); // beállítja a szálak számát, összesen egyszer hívható a program futása során
    //***********************************************************************
    void betesz(const szal_feladat_adatai & be, bool is_priority);
    //***********************************************************************
    void var_mig_van_varo_vagy_futo_feladat_barmilyen();
    //***********************************************************************
    void var_mig_van_varo_vagy_futo_feladat(szal_tipus tipus);
    //bool is_fut_e_meg(szal_tipus tipus) const;
    //***********************************************************************
    sumhiba_tipus get_sumhiba(szal_tipus tipus) const {
    //***********************************************************************
        sumhiba_tipus sumhiba = szalak[tipus].finished_list.get_sumhiba();
        rvt sum_csat = 0;
        for (uns i = 1; i < csatlakozo_aramok_dc.size(); i += 2) {
            rvt akt_hiba = csatlakozo_aramok_dc.unsafe(i) + csatlakozo_aramok_dc.unsafe(i + 1);
            sum_csat += akt_hiba*akt_hiba;
        }
        sumhiba.sum_IP_hiba += sum_csat;
        return sumhiba;
    }
    //***********************************************************************
    void torli_a_keszeket(szal_tipus tipus);
    //***********************************************************************
    //uns get_nincs_kesz_db(szal_tipus tipus) { return nincs_kesz_db[tipus]; }
    //***********************************************************************
    ~Para_engine();
    //***********************************************************************
};


// 3. párhuzamos feldolgozó

//***********************************************************************
class spinlock {
// gyors "mutex", csak rövid várakozáshoz, mivel végtelen ciklusban vár
//***********************************************************************
    std::atomic_flag is_locked = ATOMIC_FLAG_INIT;
public:
    void lock() { while (is_locked.test_and_set(std::memory_order_acquire)); }
    void unlock() { is_locked.clear(std::memory_order_release); }
};


//***********************************************************************
class nagyfelbontasu_utemezo { // high grain sceduler
//***********************************************************************

    //***********************************************************************
    struct feladat {
    //***********************************************************************
        uns tole_fuggo_blokk_1;
        uns tole_fuggo_blokk_2;
        feladat() :tole_fuggo_blokk_1{ 0 }, tole_fuggo_blokk_2{ 0 } {}
    };

    //***********************************************************************
    struct blokk {
    //***********************************************************************
        uns db;
        uns feladat_index; // a feladat vektorban hányas indexûnél kezdõdik ez a blokk
        uns next;
        uns akt_fuggoseg;
        uns eredeti_fuggoseg;
        bool is_priority;
        blokk() :db{ 0 }, feladat_index{ 0 }, next{ 0 }, akt_fuggoseg{ 0 }, eredeti_fuggoseg{ 0 }, is_priority{ false } {}
    };

    //***********************************************************************
    spinlock kivesz_slock; // az aktuális (kivevõ) listát védi, csak a kivesz lokkolhatja, hogy ne legyen deadlock
    spinlock kesz_slock; // a kész és a priority listát védi, alapból a kész lockolja, de a kivesz is lockolhatja
    uns akt_indithato;
    uns hanyadik_jon;
    uns kesz_lista_eleje;
    uns kesz_lista_vege;
    uns priority_lista_eleje;
    uns priority_lista_vege;
    vektor<feladat> feladatok;  // a 0 indexû dummy
    vektor<blokk> blokkok;      // a 0 indexû dummy
    //***********************************************************************

    //***********************************************************************
    bool belso_kesz(uns fuggo_blokk_num) {
    // lezárt kesz_slock-nál hívható
    //***********************************************************************
        bool ret = false;
        blokk & fuggo_blokk = blokkok[fuggo_blokk_num];
        if (fuggo_blokk.akt_fuggoseg == 0) {
            kivesz_slock.unlock();
            throw hiba("nagyfelbontasu_utemezo::belso_kesz", "fuggo_blokk => akt_fuggoseg == 0");
        }
        fuggo_blokk.akt_fuggoseg--;
        if (fuggo_blokk.akt_fuggoseg == 0) { // ha megszûnt a blokk függõsége, beteszi a priority vagy a kész listába
            fuggo_blokk.next = 0;
            ret = true;
            if (fuggo_blokk.is_priority) { // a priority listába teszi
                if (priority_lista_vege == 0) {
                    priority_lista_eleje = priority_lista_vege = fuggo_blokk_num;
                }
                else {
                    blokkok[priority_lista_vege].next = fuggo_blokk_num;
                    priority_lista_vege = fuggo_blokk_num;
                }
            }
            else { // a kész listába teszi
                if (kesz_lista_vege == 0) {
                    kesz_lista_eleje = kesz_lista_vege = fuggo_blokk_num;
                }
                else {
                    blokkok[kesz_lista_vege].next = fuggo_blokk_num;
                    kesz_lista_vege = fuggo_blokk_num;
                }
            }
        }
        return ret;
    }

public:
    //***********************************************************************
    nagyfelbontasu_utemezo() :akt_indithato{ 0 }, hanyadik_jon{ 0 }, kesz_lista_eleje{ 0 }, kesz_lista_vege{ 0 },
        priority_lista_eleje{ 0 }, priority_lista_vege{ 0 } {}
    //***********************************************************************

    //***********************************************************************
    uns kivesz() {
    //***********************************************************************
        kivesz_slock.lock();

        // ha van fontos feladat, azt az aktuális (kivevõ) lista elejére tesszük

        if (priority_lista_eleje != 0) { 
            kesz_slock.lock();
            blokkok[priority_lista_vege].next = akt_indithato;
            akt_indithato = priority_lista_eleje;
            priority_lista_eleje = priority_lista_vege = 0;
            kesz_slock.unlock();
        }
        // ha nincs futtatható feladat

        if (akt_indithato == 0) { 
            kivesz_slock.unlock();
            return 0;
        }

        // ha az aktuális blokk elindítva, megyünk a következõre

        blokk * p_akt_blokk = &blokkok[akt_indithato];
        if (p_akt_blokk->db < hanyadik_jon) { 
            p_akt_blokk->akt_fuggoseg = p_akt_blokk->eredeti_fuggoseg;
            hanyadik_jon = 1;
            akt_indithato = p_akt_blokk->next;
            if (akt_indithato == 0) { // átteszük a készeket
                kesz_slock.lock();
                akt_indithato = kesz_lista_eleje;
                kesz_lista_eleje = kesz_lista_vege = 0;
                kesz_slock.unlock();
                if (akt_indithato == 0) { // nincs futtatható feladat
                    kivesz_slock.unlock();
                    return 0;
                }
            }
            p_akt_blokk = &blokkok[akt_indithato];
        }

        // kivesszük az aktuális feladatot
        
        uns ret = p_akt_blokk->feladat_index + hanyadik_jon - 1; // A hanyadik_jon 1-nél indul, ezért a -1
        hanyadik_jon++;
        kivesz_slock.unlock();
        return ret;
    }

    //***********************************************************************
    bool kesz(uns kesz_feladat) { // !! A hívó feladata: ha van váró thread, akkor notify
    // return: tett-e be feladatot a futtatható listára (ekkor kell a notify)
    //***********************************************************************
        kivesz_slock.lock();
        bool ret = false;
        const feladat & akt_feladat = feladatok[kesz_feladat];

        if (akt_feladat.tole_fuggo_blokk_1 != 0) {
            ret = belso_kesz(akt_feladat.tole_fuggo_blokk_1) || ret;
        }

        if (akt_feladat.tole_fuggo_blokk_2 != 0) {
            ret = belso_kesz(akt_feladat.tole_fuggo_blokk_2) || ret;
        }

        kivesz_slock.unlock();
        return ret;
    }
    //***********************************************************************
    //*****   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! *****//
    // Most legyen külön ütemezõ a cellaklaszternek és a fw-bw-nek minden esetre,
    // a következõ verzióban a csatlakozó cella megszüntetendõ, helyette komplex cella, ami kompakt modellt is tartalmazhat,
    // SPICE struktúraleírás lehetõségével, a netlistából automatikus admittanciamátrix generálással.
    //***********************************************************************

};


}

#endif
