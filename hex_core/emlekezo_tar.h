//***********************************************************************
// emlékezõ tároló header
// Creation date:  2018. 08. 11.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef EMLEKEZO_TAR_HEADER
#define	EMLEKEZO_TAR_HEADER
//***********************************************************************


//***********************************************************************
#include "kozos.h"
#include <vector>
#include <thread>
#include <mutex>
#include "Parhuzamos_Feldolgozo.h"
//***********************************************************************

//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
template <typename T>
class egytarolo{
//***********************************************************************
    T * egy;
    egytarolo(const egytarolo &) = delete;
    egytarolo &operator=(const egytarolo &) = delete;
public:
    egytarolo(): egy{nullptr}{}
    ~egytarolo() { delete egy; }
    void alloc() { if (egy == nullptr)egy = new T; }
    bool is_allocated()const { return egy != nullptr; }
    //T & operator=(const T& be) { *egy = be; return *egy; }
    T * operator->() { return egy; }
    const T * operator->() const { return egy; }
    T & operator *() { return *egy; }
    const T & operator *() const { return *egy; }
};


//***********************************************************************
class os_emlek {
//***********************************************************************
    os_emlek(const os_emlek &) = delete;
    os_emlek & operator=(const os_emlek &) = delete;
public:
    //***********************************************************************
    virtual void leptet() = 0;
    virtual void visszalep() = 0;
    virtual void megtartando_az_aktualis() = 0;
    virtual void mentendo_az_aktualis() = 0;
    virtual void kiindulo_az_aktualis() = 0;
    os_emlek() {}
    //***********************************************************************
};


//***********************************************************************
template<typename T>
class emlek : public os_emlek {
//***********************************************************************
    size_t indexx; // Az emlekek vektorban
    T akt, prev, prevprev;
    T megtartando, mentendo, kiindulo;
public:
    //***********************************************************************
    emlek();
//    ~emlek();
    //***********************************************************************
    void leptet() override { prevprev = prev; prev = akt; }
    //***********************************************************************
    void visszalep() override { akt = prev; prev = prevprev; }
    //***********************************************************************
    void megtartando_az_aktualis() override { megtartando = akt; }
    //***********************************************************************
    void mentendo_az_aktualis() override { mentendo = akt; }
    //***********************************************************************
    void kiindulo_az_aktualis() override { kiindulo = akt; }
    //***********************************************************************
    T & get_akt() { return akt; }
    const T & get_akt() const { return akt; }
    const T & get_elozo() const { return prev; }
    T & get_elozo_to_overwrite() { return prev; }
    T & get_prevprev_to_overwrite() { return prevprev; }
    const T & get_megtartando() const { return megtartando; }
    T & get_megtartando_to_overwrite() { return megtartando; }
    const T & get_mentendo() const { return mentendo; }
    T & get_mentendo_to_overwrite() { return mentendo; }
    const T & get_kiindulo() const { return kiindulo; }
    T & get_kiindulo_to_overwrite() { return kiindulo; }
    void clear() { akt = prev = prevprev = T(); }
    void clear_megtartando_is() { clear(); megtartando = T(); }
    void init_all(const T & value) { akt = prev = prevprev = megtartando = mentendo = kiindulo = value; }
    void operator=(const T & masik) { get_akt() = masik; }
    //***********************************************************************
};


//***********************************************************************
class emlekek_ore {
//***********************************************************************
    //***********************************************************************
    std::vector<os_emlek*> tomb;
    ::std::mutex or_mutex;
    spinlock sl;
    //***********************************************************************
public:
    //***********************************************************************
    size_t push_back(os_emlek * uj) { 
    // A spinlock Intelen lassabb, AMD-n gyorsabb
    //***********************************************************************
        //::std::unique_lock<::std::mutex> lock(or_mutex);
		sl.lock();
        tomb.push_back(uj);
		size_t siz = tomb.size() - 1;
		sl.unlock();
        return siz; 
    }
    //***********************************************************************
    // void torol(size_t indexx) { tomb[indexx] = nullptr; }
    void clear() { tomb.clear(); }
    //***********************************************************************
    void leptet(size_t from = 0, size_t to = 0) {
    //***********************************************************************
        const size_t n = to == 0 ? tomb.size() : to;
        for (size_t i = from; i < n; i++)
            tomb[i]->leptet();
    }
    //***********************************************************************
    void visszalep(size_t from = 0, size_t to = 0) {
    //***********************************************************************
        const size_t n = to == 0 ? tomb.size() : to;
        for (size_t i = from; i < n; i++)
            tomb[i]->visszalep();
    }
    //***********************************************************************
    void megtartando_az_aktualis(size_t from = 0, size_t to = 0) {
    //***********************************************************************
        const size_t n = to == 0 ? tomb.size() : to;
        for (size_t i = from; i < n; i++)
            tomb[i]->megtartando_az_aktualis();
    }
    //***********************************************************************
    void mentendo_az_aktualis(size_t from = 0, size_t to = 0) {
    //***********************************************************************
        const size_t n = to == 0 ? tomb.size() : to;
        for (size_t i = from; i < n; i++)
            tomb[i]->mentendo_az_aktualis();
    }
    //***********************************************************************
    void kiindulo_az_aktualis(size_t from = 0, size_t to = 0) {
    //***********************************************************************
        const size_t n = to == 0 ? tomb.size() : to;
        for (size_t i = from; i < n; i++)
            tomb[i]->kiindulo_az_aktualis();
    }
    //***********************************************************************
    void kiindulo_az_aktualis_es_leptet(size_t from = 0, size_t to = 0) {
    //***********************************************************************
        const size_t n = to == 0 ? tomb.size() : to;
        for (size_t i = from; i < n; i++) {
            os_emlek * akt = tomb[i];
            akt->kiindulo_az_aktualis();
            akt->leptet();
        }
    }
    //***********************************************************************
    size_t get_size()const { return tomb.size(); }
    //***********************************************************************
};


//***********************************************************************
extern emlekek_ore emlekek;
//***********************************************************************


//***********************************************************************
template<typename T>
inline emlek<T>::emlek() :akt{ T() }, prev{ T() }, prevprev{ T() }, megtartando{ T() }, mentendo{ T() }, kiindulo{ T() } {
//***********************************************************************
    indexx = emlekek.push_back(this);
}


/*
//***********************************************************************
template<typename T>
inline emlek<T>::~emlek() {
//***********************************************************************
    emlekek.torol(indexx);
}
*/

}

#endif

