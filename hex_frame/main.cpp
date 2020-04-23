//***********************************************************************
// main
// Creation date:  2009. 07. 11.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#include "tipusok.h"
#include "apa.h"
#include "gfunc.h"
#include <iostream>
//***********************************************************************


//***********************************************************************
char hibaUzenet[1024];
const char * getHiba(){return hibaUzenet;}
//***********************************************************************


//***********************************************************************
bool run(const char * Project){
//***********************************************************************
    try{
        apa a(Project);
        a.write_v6sim();
    }
    catch(const hiba h){
        strcpy(hibaUzenet,h.what());
        return false;
    }
    return true;
}

char idobelyeg[100];

//***********************************************************************
int main(int argc,char **argv){ // normál sunred
//***********************************************************************
    time_t ltime;
    time(&ltime);
    strftime(idobelyeg, 100, "%y%m%d_%H%M%S", localtime(&ltime));
    if (argc < 2)
        printf("Error:\nStart: v6_frame_v2 some.proj\n");
    else
        if(!run(argv[1]))printf("\n%s\n",getHiba());
    return 0;
}

