//***********************************************************************
// vektor template header
// Creation date:  2018. 01. 06.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef VEKTOR_HEADER
#define	VEKTOR_HEADER
//***********************************************************************


//***********************************************************************
#include "kozos.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
#define epszilon 1e-8
//***********************************************************************


//***********************************************************************
template<typename adattipus> class vektor {
//***********************************************************************
    adattipus *tomb;
    meret_t n;
    bool torlendo;
    //***********************************************************************
    vektor(const vektor&) = delete;
    vektor(const vektor&&) = delete;
    void operator=(const vektor&&) = delete;
    //***********************************************************************
public:
    //***********************************************************************
    vektor() :tomb{ nullptr }, n{ 0 }, torlendo(false) {}
    //***********************************************************************
    // vektor(adattipus * p, meret_t size) :tomb{ p }, n{ size }, torlendo(false) {}
    //***********************************************************************
    // vektor(const adattipus * p, meret_t size) :tomb{ const_cast<adattipus*>(p) }, n{ size }, torlendo(false) {}
    //***********************************************************************
    // vektor(vektor & v, meret_t size, meret_t from = 0) :tomb{ v.tomb + from }, n{ size }, torlendo(false) {kisebb_e_hiba(from+size, v.n+1, "vektor ctor");}
    //***********************************************************************
    // vektor(const vektor & v, meret_t size, meret_t from = 0) :tomb{ const_cast<adattipus*>(v.tomb + from) }, n{ size }, torlendo(false) {kisebb_e_hiba(from + size, v.n + 1, "vektor ctor");}
    //***********************************************************************
    ~vektor() { clear_tagvaltozobeallitas_nelkul(); }
    //***********************************************************************
    meret_t size()const { return n; }
    //***********************************************************************
    void clear_tagvaltozobeallitas_nelkul() { if (torlendo)delete[]tomb; }
    //***********************************************************************
    void clear() { clear_tagvaltozobeallitas_nelkul(); tomb = nullptr; n = 0; torlendo = false; }
    //***********************************************************************
    void set_size(meret_t new_size) { clear_tagvaltozobeallitas_nelkul(); n = new_size; tomb = (n == 0) ? nullptr : new adattipus[n]; torlendo = true; }
    //***********************************************************************
    void set_size_and_zero(meret_t new_size) { set_size(new_size); zero(); }
    //***********************************************************************
    adattipus & operator[](meret_t i) { kisebb_e_hiba(i, n, "vektor::operator[]"); return tomb[i]; }
    //***********************************************************************
    const adattipus & operator[](meret_t i)const { kisebb_e_hiba(i, n, "vektor::operator[]const"); return tomb[i]; }
    //***********************************************************************
    adattipus & unsafe(meret_t i) { return tomb[i]; }
    //***********************************************************************
    const adattipus & unsafe(meret_t i)const { return tomb[i]; }
    //***********************************************************************
    const adattipus & first()const { return operator[](0); }
    //***********************************************************************
    const adattipus & last()const { return operator[](n-1); }
    //***********************************************************************
    adattipus & last() { return operator[](n - 1); }
    //***********************************************************************
        
    //***********************************************************************
    vektor & operator=(const vektor& masik){
    //***********************************************************************
        set_size(masik.n);
        for (uns i = 0; i < n; i++)
            tomb[i] = masik.tomb[i];
        return *this;
    }


    //***********************************************************************
    void debug_write(::std::ofstream & fs) const{
    //***********************************************************************
        for (meret_t i = 0; i < n; i++) {
            if (i % 128 == 0) {
                fs << (i == 0 ? "   " : "\n>> ");
            }
            fs << ::std::setw(12) << tomb[i] << ' ';
        }
        fs << ::std::endl;
    }


    //***********************************************************************
    void print(uns kezd = 0) const{
    //***********************************************************************
        for (meret_t i = 0; i < kezd; i++)
            ::std::cout << "- ";
        for (meret_t i = kezd; i < n; i++)
            ::std::cout << tomb[i] << ' ';
    }


    //***********************************************************************
    void print_z(uns kezd = 0) const{
    //***********************************************************************
        for (meret_t i = 0; i < kezd; i++)
            ::std::cout << "- ";
        for (meret_t i = kezd; i < n; i++)
            ::std::cout << (abs(tomb[i]) < 1e-10 ? 0 : tomb[i]) << ' ';
    }


    //***********************************************************************
    void rafektet(vektor & masik){
    // a this vektor a másikra mutasson.
    //***********************************************************************
        clear_tagvaltozobeallitas_nelkul();
        tomb = masik.tomb;
        n = masik.n;
        torlendo = false;
    }


    //***********************************************************************
    void rafektet(vektor & masik, meret_t kezdet, meret_t db){
    // a this vektor a másikra mutasson.
    //***********************************************************************
        clear_tagvaltozobeallitas_nelkul();
        n = db;
        torlendo = false;
        tomb = (db == 0) ? nullptr : (masik.tomb + kezdet);
    }

/*
    //***********************************************************************
    bool frissit(meret_t i, const adattipus & uj){ // return: false, ha nincs frissítve (egyformák)
    //***********************************************************************
        kisebb_e_hiba(i, n, "vektor::frissit"); 
        if (abs(tomb[i] - uj) <= epszilon*abs(uj) + 1.0e-10) // < nem jó, mert regi=uj=0 esetén false-t ad
            return false;
        tomb[i] = uj;
        return true;
    }
*/

    //***********************************************************************
    bool frissit_unsafe(meret_t i, const adattipus & uj){ // return: false, ha nincs frissítve (egyformák)
    //***********************************************************************
        double egyik = abs(tomb[i] - uj);
        double masik = epszilon*abs(uj) + 1.0e-10; // < nem jó, mert regi=uj=0 esetén false-t ad
        if (egyik == 0)
            return false;
        //if (egyik <= masik) // egyelõre kiveszem a frissítést, így mindenképpen másol
        //    return false;
        tomb[i] = uj;
        return true;
    }
    //***********************************************************************


    //***********************************************************************
    void subvektor_copy(const vektor & src, uns start_dest, uns start_src, uns db){
    //***********************************************************************
        if (db == 0)
            return;
        kisebb_e_hiba(start_dest + db - 1, n,     "vektor::subvektor_copy start_dest + db");
        kisebb_e_hiba(start_src  + db - 1, src.n, "vektor::subvektor_copy start_src + db");
        for (uns i = 0; i < db; i++)
            unsafe(start_dest + i) = src.unsafe(start_src + i);
    }


    //***********************************************************************
    void subvektor_pluszegyenlo(const vektor & src, uns start_dest, uns start_src, uns db){
    //***********************************************************************
        if (db == 0)
            return;
        kisebb_e_hiba(start_dest + db - 1, n,     "vektor::subvektor_pluszegyenlo start_dest + db");
        kisebb_e_hiba(start_src  + db - 1, src.n, "vektor::subvektor_pluszegyenlo start_src + db");
        for (uns i = 0; i < db; i++)
            unsafe(start_dest + i) += src.unsafe(start_src + i);
    }


    //***********************************************************************
    bool subvektor_frissit(const vektor & src, uns start_dest, uns start_src, uns db){
    //***********************************************************************
        if (db == 0)
            return false;
        kisebb_e_hiba(start_dest + db - 1, n,     "vektor::subvektor_frissit start_dest + db");
        kisebb_e_hiba(start_src  + db - 1, src.n, "vektor::subvektor_frissit start_src + db");
        bool is_valtozott = false;
        for (uns i = 0; i < db; i++)
            is_valtozott = frissit_unsafe(start_dest + i, src.unsafe(start_src + i)) || is_valtozott;
        return is_valtozott;
    }

/*
    //***********************************************************************
    bool frissit(const vektor & src){
    //***********************************************************************
        egyforma_e_hiba(n, src.n, "vektor::frissit n!=src.n");
        bool is_valtozott = false;
        for (uns i = 0; i < n; i++)
            is_valtozott = frissit_unsafe(i, src.unsafe(i)) || is_valtozott;
        return is_valtozott;
    }
*/

    //***********************************************************************
    void subvektor_add(const vektor & src_1, const vektor & src_2, uns start_dest, uns start_src_1, uns start_src_2, uns db){
    //***********************************************************************
        if (db == 0)
            return;
        kisebb_e_hiba(start_dest  + db - 1, n,       "vektor::subvektor_copy start_dest + db");
        kisebb_e_hiba(start_src_1 + db - 1, src_1.n, "vektor::subvektor_copy start_src_1 + db");
        kisebb_e_hiba(start_src_2 + db - 1, src_2.n, "vektor::subvektor_copy start_src_2 + db");
        for (uns i = 0; i < db; i++)
            unsafe(start_dest + i) = src_1.unsafe(start_src_1 + i) + src_2.unsafe(start_src_2 + i);
    }


    //***********************************************************************
    void push_back(const adattipus & uj){
    //***********************************************************************
        adattipus *t2 = new adattipus[n + 1];
        for (uns i = 0; i < n; i++)
            t2[i] = tomb[i];
        if (torlendo)
            delete[]tomb;
        tomb = t2;
        tomb[n] = uj;
        n++;
        torlendo = true;
    }


    //***********************************************************************
    void zero(){
    //***********************************************************************
        for (uns i = 0; i < n;i++)
            tomb[i] = adattipus();
    }


    //***********************************************************************
    void math_add(const vektor & v1, const vektor & v2) {
    //***********************************************************************
        egyforma_e_hiba(v1.n, v2.n, "vektor::add v1.n!=v2.n");
        egyforma_e_hiba(this->n, v2.n, "vektor::add this->n!=input.n ");
        // TODO: érdemes kis vektorokra kézzel kifejteni?
        for (meret_t i = 0; i < n; i++)
            tomb[i] = v1.tomb[i] + v2.tomb[i];
    }
    //***********************************************************************


    //***********************************************************************
    void math_sub(const vektor & v1, const vektor & v2) {
    //***********************************************************************
        egyforma_e_hiba(v1.n, v2.n, "vektor::sub v1.n!=v2.n");
        egyforma_e_hiba(this->n, v2.n, "vektor::sub this->n!=input.n ");
        // TODO: érdemes kis vektorokra kézzel kifejteni?
        for (meret_t i = 0; i < n; i++)
            tomb[i] = v1.tomb[i] - v2.tomb[i];
    }
    //***********************************************************************


    //***********************************************************************
    void math_neg() {
    //***********************************************************************
        for (meret_t i = 0; i < n; i++)
            tomb[i] = -tomb[i];
    }
    //***********************************************************************


    //***********************************************************************
    friend inline adattipus math_mul(const vektor & v1, const vektor & v2) {
    //***********************************************************************
        egyforma_e_hiba(v1.n, v2.n, "vektor mul v1.n!=v2.n");
        // TODO: érdemes kis vektorokra kézzel kifejteni?
        adattipus sum = adattipus();
        for (meret_t i = 0; i < v1.n; i++)
            sum += v1.tomb[i] * v2.tomb[i];
        return sum;
    }
    //***********************************************************************
};


//***********************************************************************
template<typename adattipus> class ref_vektor {
//***********************************************************************

    //***********************************************************************
    vektor<adattipus *> ptomb;
    //***********************************************************************
public:
    //***********************************************************************
    meret_t size()const { return ptomb.size(); }
    //***********************************************************************
    void set_ref(meret_t i, adattipus & p) { ptomb[i] = &p; }
    //***********************************************************************
    void set_pointer(meret_t i, adattipus * p) { ptomb[i] = p; }
    //***********************************************************************
    void clear() { ptomb.clear(); }
    //***********************************************************************
    void set_size(meret_t new_size) { ptomb.set_size(new_size); }
    //***********************************************************************
    adattipus & operator[](meret_t i) { return *ptomb[i]; }
    //***********************************************************************
    const adattipus & operator[](meret_t i)const { return *ptomb[i]; }
    //***********************************************************************
};


//***********************************************************************
template<typename adattipus> class lista {
//***********************************************************************
    struct listaelem {
        listaelem *prev, *next;
        adattipus adat;
        listaelem() :prev{ nullptr }, next{ nullptr } {}
        listaelem(const adattipus & uj) :prev{ nullptr }, next{ nullptr }, adat(uj) {}
    };
    listaelem start, stop;
public:
    lista() { start.next = &stop; stop.prev = &start; }
    ~lista() { for (listaelem * it = start.next->next; it != nullptr; it = it->next) delete it->prev; }
    adattipus &push_back(const adattipus & uj) { listaelem *p = new listaelem(uj); p->prev = stop.prev; p->next = &stop; stop.prev->next = p; stop.prev = p; return p->adat; }
    adattipus &push_back() { listaelem *p = new listaelem; p->prev = stop.prev; p->next = &stop; stop.prev->next = p; stop.prev = p; return p->adat; }
    adattipus &push_front(const adattipus & uj) { listaelem *p = new listaelem(uj); p->prev = &start; p->next = start.next; start.next->prev = p; start.next = p;  return p->adat; }
    adattipus &push_front() { listaelem *p = new listaelem; p->prev = &start; p->next = start.next; start.next->prev = p; start.next = p;  return p->adat; }
    adattipus pop_front() { listaelem *p = start.next; adattipus vissza = p->adat; p->next->prev = &start; start.next = p->next; delete p; return vissza; }
    void *get_it() { return start.next == &stop ? nullptr : start.next; }
    const void *get_it() const { return start.next == &stop ? nullptr : start.next; }
    int inc_it(void *& it) { listaelem *p = (listaelem*)it; if (p->next == &stop) it = nullptr; else it = p->next; return 0; }
    int inc_it(const void *& it)const { const listaelem *p = (const listaelem*)it; if (p->next == &stop) it = nullptr; else it = p->next; return 0; }
    adattipus & get_akt(void * it) { listaelem *p = (listaelem*)it; return p->adat; }
    const adattipus & get_akt(const void * it) const { const listaelem *p = (const listaelem*)it; return p->adat; }
    bool is_empty()const { return start.next == &stop; }
};


}

#endif

