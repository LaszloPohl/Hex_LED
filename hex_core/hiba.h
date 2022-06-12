//***********************************************************************
// hiba header
// Creation date:  2009. 07. 11.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef PL_HIBA_HEADER
#define	PL_HIBA_HEADER
//***********************************************************************
#include <cstdio>
#include <cstdarg>
#include <string>
#include <fstream>
#include <sstream>
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
class hiba {
//***********************************************************************
    char t[1024];
    static ::std::string hol; // tömeges hibaellenõrzéskor ne kelljen mindig megadni a helyet, a függvény elején beállítjuk
    friend class set_hiba_hol;
public:

    //***********************************************************************
    hiba(const char * hely, const char * formatum, ...) {
    //***********************************************************************
        va_list p;
        char s[1024];

        va_start(p, formatum);
        sprintf_s(s, 1024, "Error: %s => %s\n", hely, formatum);
        vsprintf_s(t, 1024, s, p);
        va_end(p);
    }
    //***********************************************************************

    //***********************************************************************
    hiba(int, const char * formatum, ...) {
    //***********************************************************************
        va_list p;
        char s[1024];

        va_start(p, formatum);
        sprintf_s(s, 1024, "Error: %s => %s\n", hol.c_str(), formatum);
        vsprintf_s(t, 1024, s, p);
        va_end(p);
    }
    //***********************************************************************

    //***********************************************************************
    const char * what()const { return t; }
    //***********************************************************************

};
//***********************************************************************

//***********************************************************************
#define hibaellenorzes
#define todo_hiba
//***********************************************************************


//***********************************************************************
template<typename T>
inline void egyforma_e_hiba(T szam_1, T szam_2, const char * hely) {
//***********************************************************************
#ifdef hibaellenorzes
    if (szam_1 != szam_2) {
        std::stringstream ss;
        ss << szam_1 << "!=" << szam_2;
        throw hiba(hely, ss.str().c_str());
    }
#endif
}


//***********************************************************************
template<typename T>
inline void kisebb_e_hiba(T szam_1, T szam_2, const char * hely) {
//***********************************************************************
#ifdef hibaellenorzes
    if (szam_1 >= szam_2) {
        std::stringstream ss;
        ss << szam_1 << ">=" << szam_2;
        throw hiba(hely, ss.str().c_str());
    }
#endif
}


//***********************************************************************
template<typename T>
inline void nagyobb_e_hiba(T szam_1, T szam_2, const char * hely) {
//***********************************************************************
#ifdef hibaellenorzes
    if (szam_1 <= szam_2) {
        std::stringstream ss;
        ss << szam_1 << "<=" << szam_2;
        throw hiba(hely, ss.str().c_str());
    }
#endif
}


//***********************************************************************
inline void igaz_e_hiba(bool is_hiba, const char * hely, const char * mi) {
//***********************************************************************
#ifdef hibaellenorzes
    if (is_hiba) {
        throw hiba(hely, mi);
    }
#endif
}


//***********************************************************************
inline void TODO(const char * mi) {
//***********************************************************************
#ifdef todo_hiba
    throw hiba("TODO", mi);
#endif
}


//***********************************************************************
inline void fail_to_read_error_check(::std::ifstream & fs, const char *mit) {
//***********************************************************************
    if(!fs.good())
        throw hiba(1, "cannot read %s", mit);
}


//***********************************************************************
class set_hiba_hol {
//***********************************************************************
    ::std::string old_hiba;
public:
    set_hiba_hol(::std::string hol) : old_hiba{ hiba::hol } { hiba::hol = hol; }
    ~set_hiba_hol() { hiba::hol = old_hiba; }
};
//***********************************************************************


//***********************************************************************
extern time_t prog_start_time;
//***********************************************************************


//***********************************************************************
class log_print {
//***********************************************************************
    static bool is_kepernyore;
    static bool is_fajlba;
public:
    static void set_kepernyore(bool is) { is_kepernyore = is; }
    static void set_fajlba(bool is) { is_fajlba = is; }
    log_print(const char * formatum, ...);
    log_print(int, const char * formatum, ...);
};

}

//***********************************************************************
#endif
//***********************************************************************

