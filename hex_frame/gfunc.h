//***********************************************************************
// general function header
// Creation date:  2009. 07. 12.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef VSUN3_GFUNC_HEADER
#define	VSUN3_GFUNC_HEADER
//***********************************************************************


//***********************************************************************
#include "tipusok.h"
//***********************************************************************


//***********************************************************************
inline double replusz(const double a,const double b){return a==nulla ? nulla : b==nulla ? nulla : a==g_max ? b : b==g_max ? a : a*b/(a+b);}
inline dbl negyzet(dbl a){return a*a;}
void GetTime(PLString & ret);
PLString GetTmpName();
void logprint(const char * s,...);
const char * getlogline(); // GUI-ból hívható
void resetAnalLevel(unsigned level);
void setAnalStepDb(unsigned level,unsigned db);
void setAnalKeszStep(unsigned level,unsigned n,const char * kovlepesnev="");
unsigned getAnalProgress(unsigned level);
inline dbl fn_mizs1(dbl T, dbl a, dbl b, dbl c, dbl d, dbl e, dbl f, dbl g) {
    T += absT;
    cd Tbc = pow(T - b, c);
    cd signa = T > b ? -a : a;
    cd tag1 = signa * Tbc / (Tbc + d) + a + f;
    return 0.01*(tag1*exp(e / T) + g);
}
//***********************************************************************


//***********************************************************************
#endif
//***********************************************************************
