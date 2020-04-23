//***********************************************************************
// hiba header
// Creation date:  2009. 07. 11.
// Creator:        Pohl L�szl�
//***********************************************************************


//***********************************************************************
#ifndef PL_HIBA_HEADER
#define	PL_HIBA_HEADER
//***********************************************************************
#include <stdio.h>
#include <stdarg.h>
//***********************************************************************


//***********************************************************************
class hiba{
//***********************************************************************
    char t[1024];
public:
    hiba(const char * hely,const char * formatum,...){
	    va_list p;
        char s[1024];

	    va_start(p,formatum);
        sprintf(s,"Error: %s => %s",hely,formatum);
	    vsprintf(t,s,p);	
	    va_end(p);
    }
    const char * what()const{return t;}
};
//***********************************************************************


//***********************************************************************
#endif
//***********************************************************************