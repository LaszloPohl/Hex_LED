//***********************************************************************
// hiba cpp
// Creation date:  2018. 07. 27.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#include "hiba.h"
#include <ctime>
#include "fajlolvasas_segito_rutinok.h"
//***********************************************************************


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
::std::string hiba::hol;
bool log_print::is_kepernyore = true;
bool log_print::is_fajlba = false;
::std::ifstream * fajl::pfs = nullptr;
//***********************************************************************


//***********************************************************************
log_print::log_print(const char * formatum, ...){
//***********************************************************************
    if (is_kepernyore) {
        va_list p;
        va_start(p, formatum);
        vprintf(formatum, p);
        va_end(p);
    }
    if (is_fajlba) {
        va_list p;
        char ido[1024];
        char logline[1024];
        FILE *fp;
        if(fopen_s(&fp, "srlog.log", "at")!=0)
            return;
        time_t ltime;
        time(&ltime);
        ctime_s(ido, 1024, &ltime);
        ido[24] = '\0';
        va_start(p, formatum);
        sprintf_s(logline, 1024, "[%s] %s\n", ido + 4, formatum);
        if (formatum[0])vfprintf(fp, logline, p);
        else fprintf(fp, "\n");
        va_end(p);
        fclose(fp);
    }
}

time_t prog_start_time;

//***********************************************************************
log_print::log_print(int, const char * formatum, ...){
//***********************************************************************
    char ido[1024];
    char logline[1024];
    time_t ltime;
    time(&ltime);
    ctime_s(ido, 1024, &ltime);
    uns dtime = (uns)(ltime - prog_start_time);
    if (is_kepernyore) {
        ido[19] = '\0';
        const std::string * ps;
        uns szazalek, akt_lepesszam, ossz_lepesszam;
        akt_lepes::get_aktualis_lepes_szazalek(ps, szazalek, akt_lepesszam, ossz_lepesszam);
        if (ossz_lepesszam > 0)
            sprintf_s(logline, 1024, "[%s => %02u:%02u:%02u] step %u/%u of %s", ido + 11, dtime / 3600, dtime % 3600 / 60, dtime % 60, akt_lepesszam, ossz_lepesszam, formatum);
        else
            sprintf_s(logline, 1024, "[%s => %02u:%02u:%02u] %s", ido + 11, dtime / 3600, dtime % 3600 / 60, dtime % 60, formatum);
        va_list p;
        va_start(p, formatum);
        if (formatum[0])vprintf(logline, p);
        else printf("\n");
        va_end(p);
    }
    if (is_fajlba) {
        ido[24] = '\0';
        sprintf_s(logline, 1024, "[%s => %02u:%02u:%02u] %s", ido + 4, dtime / 3600, dtime % 3600 / 60, dtime % 60, formatum);
        va_list p;
        FILE *fp;
        if(fopen_s(&fp, "srlog.log", "at")!=0 || fp==nullptr)
            return;
        va_start(p, formatum);
        if (formatum[0])vfprintf(fp, logline, p);
        else fprintf(fp, "\n");
        va_end(p);
        fclose(fp);
    }
}

}