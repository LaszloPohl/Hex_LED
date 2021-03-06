//***********************************************************************
// f�jlolvas�st seg�t� rutinok
// Creation date:  2018. 07. 27.
// Creator:        Pohl L�szl�
//***********************************************************************


//***********************************************************************
#ifndef FOSR_HEADER
#define	FOSR_HEADER
//***********************************************************************


//***********************************************************************
#include "kozos.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
inline ::std::string uns_to_string(uns szam){
//***********************************************************************
    ::std::stringstream ss;
    ss << szam;
    return ss.str();
}


//***********************************************************************
class fajl {
//***********************************************************************
    static ::std::ifstream * pfs;
public:

    //***********************************************************************
    fajl(const ::std::string & fajlnev) {
    //***********************************************************************
        pfs = new ::std::ifstream(fajlnev);
        if (!pfs->is_open())
            throw hiba(1, "Cannot open \"%s\" to read.", fajlnev.c_str());
    }


    //***********************************************************************
    ~fajl() { pfs->close(); delete pfs; pfs = nullptr; }
    //***********************************************************************


    //***********************************************************************
    static int get_char(const char * mit) {
    //***********************************************************************
        int ch = 0;
        while (isspace(ch = pfs->get()) && pfs->good())
            ;
        if (ch == '#') {
            while (pfs->get() != '\n' && pfs->good())
                ;
            if (pfs->good())
                ch = get_char(mit);
        }
        if (!pfs->good())
            throw hiba(1, "cannot read %s", mit);
        return ch;
    }
    //***********************************************************************


    //***********************************************************************
    static int get_char_or_eol(const char * mit) {
    //***********************************************************************
        int ch = 0;
        while ((ch = pfs->get()) != '\n' && isspace(ch) && pfs->good())
            ;
        if (ch == '#') {
            while ((ch = pfs->get()) != '\n' && pfs->good())
                ;
        }
        if (!pfs->good())
            throw hiba(1, "cannot read %s", mit);
        return ch == '\n' ? 0 : ch;
    }
    //***********************************************************************


    //***********************************************************************
    static void check_text(const char * mit, const char * mit2 = "") {
    //***********************************************************************
        for (uns i = 0; mit[i] != '\0'; i++)
            if (get_char(mit) != mit[i])
                throw hiba(1, "cannot read %s %s", mit2, mit);
    }
    //***********************************************************************


    //***********************************************************************
    static uns get_uns(const char * mit) {
    // a megjegyz�st nem kezeli
    //***********************************************************************
        uns adat;
        (*pfs) >> adat;
        if (!pfs->good())
            throw hiba(1, "cannot read %s", mit);
        return adat;
    }
    //***********************************************************************


    //***********************************************************************
    static rvt get_rvt(const char * mit, bool is_pontosvesszo) {
    // a megjegyz�st nem kezeli
    //***********************************************************************
        rvt adat;
        (*pfs) >> adat;
        if (!pfs->good())
            throw hiba(1, "cannot read %s", mit);
        if (is_pontosvesszo)
            check_text(";", mit);
        return adat;
    }
    //***********************************************************************


    //***********************************************************************
    static ::std::string get_quoted_text(const char * mit) {
    //***********************************************************************
        ::std::string adat;
        check_text("\"", mit);
        int ch;
        while ((ch = get_char(mit)) != '\"')
            adat += char(ch);
        return adat;
    }
    //***********************************************************************


    //***********************************************************************
    static void eldobja_a_sor_veget() {
    //***********************************************************************
        while (pfs->get() != '\n' && pfs->good())
            ;
    }

};

}

#endif

