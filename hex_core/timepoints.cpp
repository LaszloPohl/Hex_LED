//***********************************************************************
// futásidõ számítása cpp
// Creation date:  2018. 01. 14.
// Creator:        Pohl László
//***********************************************************************

//***********************************************************************
#include <chrono>
#include <string>
#include <vector>
#include <iostream>
#include "kozos.h"
//***********************************************************************


//***********************************************************************
using namespace ::std;
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
class Time_points {
//***********************************************************************

    //***********************************************************************
    struct idopont {
    //***********************************************************************
        chrono::time_point<chrono::system_clock> ido;
        string mi_tortent;
    };
    //***********************************************************************

    //***********************************************************************
    vector<idopont> tomb;
    //***********************************************************************
    bool csak_masodpercben;
    //***********************************************************************

    //***********************************************************************
    void idokonverter(double t, char * puffer) {
    //***********************************************************************
        char ctomb[128];
        if (t < 1.0 && !csak_masodpercben) {
            sprintf_s(ctomb, 128, "%9.3f ms,", t*1000.0);
        }
        else {
            unsigned h = (unsigned)(t / 3600);
            unsigned m = (unsigned)(t / 60) - h * 60;
            double s = t - 60 * m - 3600 * h;
            unsigned is = (unsigned)s;
            unsigned ms = (unsigned)((s - is)*1000.0 + 0.5);
            if (t < 60) {
                if (is == 0 && ms == 0)
                    sprintf_s(ctomb, 128, "        0  s,");
                else if(ms == 0)
                    sprintf_s(ctomb, 128, "   %2u      s,", is);
                else if(is == 0)
                    sprintf_s(ctomb, 128, "      %3u ms,", ms);
                else
                    sprintf_s(ctomb, 128, "   %2u.%03u  s,", is, ms);
            }
            else if (t < 3600)
                sprintf_s(ctomb, 128, "%2u:%02u.%03u  s,", m, is, ms);
            else
                sprintf_s(ctomb, 128, "%2u:%02u:%02u.%03u  s,", h, m, is, ms);
        }
        sprintf_s(puffer, 128, "%-14s", ctomb);
    }
public:

    //***********************************************************************
    Time_points() :csak_masodpercben{ false } { add("Program started"); }
    //***********************************************************************

    //***********************************************************************
    void add(const string & mi_tortent) {
    //***********************************************************************
        idopont t;
        t.ido = chrono::system_clock::now();
        t.mi_tortent = mi_tortent;
        tomb.push_back(t);
    }

    //***********************************************************************
    void print() {
    //***********************************************************************
        cout << "                                                              " << endl;
        chrono::time_point<chrono::system_clock> start = tomb[0].ido;
        char ido[128], form[128];
        uns maxhossz = 0;
        for (uns i = 1; i < tomb.size();i++)
            if (tomb[i].mi_tortent.length()>maxhossz)
                maxhossz = (uns)tomb[i].mi_tortent.length();
        sprintf_s(form, 128, "  %%-%us", maxhossz + 2);
        for (uns i = 1; i < tomb.size(); i++) {
            sprintf_s(ido, 128, form, (tomb[i].mi_tortent + ": ").c_str());
            cout << ido;
            chrono::duration<double> teljes_ido_eddig = tomb[i].ido - start;
            chrono::duration<double> dt = tomb[i].ido - tomb[i - 1].ido;
            idokonverter(dt.count(), ido);
            cout << ido;
            idokonverter(teljes_ido_eddig.count(), ido);
            cout << "full time: " << ido << endl;
        }
        cout << endl;
    }

    //***********************************************************************
    void print_last() {
    //***********************************************************************
        chrono::time_point<chrono::system_clock> start = tomb[0].ido;
        char ido[128], form[128];
        uns maxhossz = 0;
        if (tomb.size()>1) {
            maxhossz = (uns)tomb.back().mi_tortent.length();
            sprintf_s(form, 128, "  %%-%us", maxhossz + 2);
            sprintf_s(ido, 128, form, (tomb.back().mi_tortent + ": ").c_str());
            cout << ido;
            chrono::duration<double> teljes_ido_eddig = tomb.back().ido - start;
            chrono::duration<double> dt = tomb.back().ido - tomb[tomb.size() - 2].ido;
            idokonverter(dt.count(), ido);
            cout << ido;
            idokonverter(teljes_ido_eddig.count(), ido);
            cout << "full time: " << ido << endl;
        }
    }

    //***********************************************************************
    void kiirt_ido_forma(bool masodpercben) {
    //***********************************************************************
        csak_masodpercben = masodpercben;
    }
};
//***********************************************************************

//***********************************************************************
Time_points timepoints;
//***********************************************************************
::std::string akt_lepes::aktualis_lepes_neve;
//***********************************************************************
uns akt_lepes::aktualis_lepes_szazalek = 0, akt_lepes::ossz_lepesszam = 0, akt_lepes::aktualis_lepes_szama = 0;
//***********************************************************************
bool akt_lepes::akt_lepes_kiir = true;
//***********************************************************************

//***********************************************************************
void most(const ::std::string & mi_tortent) {
//***********************************************************************
    timepoints.add(mi_tortent);
}

//***********************************************************************
void esemenyek_kiirasa() {
//***********************************************************************
    timepoints.print();
}

//***********************************************************************
void utolso_esemeny_kiirasa() {
//***********************************************************************
    timepoints.print_last();
}

//***********************************************************************
void kiirt_ido_forma(bool masodpercben) {
//***********************************************************************
    timepoints.kiirt_ido_forma(masodpercben);
}

}