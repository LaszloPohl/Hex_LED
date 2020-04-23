//***********************************************************************
// general function source
// Creation date:  2009. 07. 12.
// Creator:        Pohl L�szl�
//***********************************************************************


//***********************************************************************
#include "gfunc.h"
//***********************************************************************


//***********************************************************************
u32 VSUN_CPU_Thread=8; // csak a beolvas� �ll�thatja
u32 VSUN_Akt_Thread=0; // nodered �ll�tja be, a matrix n�velheti, szalak �s matrix cs�kkenti
bool parhuzamosanfut=false;
const char *OldalNev[]={"CENTER","WEST","EAST","SOUTH","NORTH","BOTTOM","TOP"};
//***********************************************************************


//***********************************************************************
void GetTime(PLString & ret){
//***********************************************************************
    time_t ltime;
    time(&ltime);
    ret=ctime(&ltime)+4;
    ret.trunc(20);
}


//***********************************************************************
PLString GetTmpName(){
//***********************************************************************
	static u32 n=0;
	char to[64];
	sprintf(to,"tmp%05u.tmp",n++);
	return PLString(to);
}


static clock_t start=clock();
static clock_t gyujt[10],gys;

void GyujtStart(){gys=clock();}
void GyujtStop(cu32 n){gyujt[n]+=clock()-gys;}
void GyzjtPrint(cu32 n){printf("%lu: %.3f seconds\n",n,(double)gyujt[n]/CLOCKS_PER_SEC);}

//***********************************************************************
void ResetRunTimeClock(){
//***********************************************************************
	start=clock();
}

//***********************************************************************
double GetRunTime(PLString & ret){
//***********************************************************************
	clock_t finish=clock();
	char s[128];
	double duration=(double)(finish-start)/CLOCKS_PER_SEC;
	sprintf(s,"%.3f seconds",duration);
	ret=s;
	start=finish;
	return duration;
}

//***********************************************************************
void printTime(const char * muvelet,const char * projektnev,const char *res, uns db){
//***********************************************************************
	PLString s;
	double ido=GetRunTime(s);
	double ido1=ido/db;
	char dbido[256];
	sprintf(dbido,"%.4g %s",ido1<1.0?ido1*1000.0:ido1,ido1<1.0?"ms":"s");
	if(db==1)printf("* %-30s = %s\n",muvelet,s.c_str());
	else printf("* %-30s = %u x %s = %s\n",muvelet,db,dbido,s.c_str());
    FILE *fp=fopen("vsunruntime.txt","at");
    if(fp!=NULL){
        fprintf(fp,"%s, %s, %s = %u x %s = %s\n",projektnev,res,muvelet,db,dbido,s.c_str());
    }
    fclose(fp);
}

//***********************************************************************
char logline[1024];
const char * getlogline(){return logline;}
void logprint(const char * formatum,...){
//***********************************************************************
	va_list p;
    FILE *fp=fopen("srlog.log","at");
    if(fp==NULL)return;
    PLString ido;
    GetTime(ido);
	va_start(p,formatum);
    sprintf(logline,"[%s] %s\n",ido.c_str(),formatum);
	if(formatum[0])vfprintf(fp,logline,p);
    else fprintf(fp,"\n");
	va_end(p);
    fclose(fp);
}


//***********************************************************************
// [0]: h�ny l�p�s egy dolog (nodered, forwsubs, backsubs)
// [1]: h�ny dolog egy DC/AC anal�zis
// [2]: h�ny DC/AC anal�zis az anal�zis
// [3]: h�ny anal�zis a szimul�ci�
// [4]: h�ny szimul�ci� a projekt
unsigned StepDb[5];
unsigned KeszStep[5]; // h�ny step van k�sz
PLString AktFolyamatNev;
bool ConsoleText=false;
bool ConDebug=false;
ProgressTipus proti=PT_Egyenletes;
bool is_simple_anal_progress = true;

//***********************************************************************
void kiirAnal(){
//***********************************************************************
    if(!ConsoleText)return;
    if(ConDebug){
        char s[60];
        sprintf(s,"%s: %2u/%u               ",AktFolyamatNev.c_str(),KeszStep[0],StepDb[0]);
        s[22]=0;
        printf("\r* sim: %u/%u, anal: %u/%u, subanal: %u/%u, step: %u/%u - %s",
            KeszStep[4],StepDb[4],KeszStep[3],StepDb[3],KeszStep[2],StepDb[2],KeszStep[1],StepDb[1],s);
    }
    else{
//        ConDebug=true;
//        kiirAnal();
//        ConDebug=false;
        if (is_simple_anal_progress) {
            printf("\r* [%u/%u] (%u-%u=%u)                 ", KeszStep[2], StepDb[2], StepDb[2], KeszStep[2], StepDb[2]-KeszStep[2]);
        }
        else {
            cuns anpr4 = getAnalProgress(1);
            char progress[80];
            sprintf(progress, "[%u/%u]", KeszStep[2], StepDb[2]);
            printf("\r* %-9s %3u%% ", progress, anpr4);
            cuns n = anpr4 > 100 ? 50 : anpr4 * 50 / 100;
            for (uns i = 0;i < n;i++)progress[i] = '#';
            progress[n] = 0;
            printf("%s", progress);
        }
    }
//    if(StepDb[4]>1)printf(", sim: %3u%%",getAnalProgress(3));
//    if(StepDb[3]>1)printf(", anal: %3u%%, sub anal %3u%%, step %3u%% of %.15s ",getAnalProgress(2),getAnalProgress(1),getAnalProgress(0),AktFolyamatNev.c_str());
//    else printf(", anal: %3u%%, step %3u%% of %-15s ",getAnalProgress(1),getAnalProgress(0),AktFolyamatNev.c_str());
}


//***********************************************************************
void resetAnalLevel(unsigned level){
//***********************************************************************
    if(level>4)throw hiba("resetAnalLevel","level>4 (%u>4)",level);
    for(uns i=0;i<=level;i++){
        StepDb[i]=1;
        KeszStep[i]=0;
    }
}


//***********************************************************************
void setAnalStepDb(unsigned level,unsigned db){
//***********************************************************************
    if(level>4)throw hiba("setAnalStepDb","level>4 (%u>4)",level);
    StepDb[level]=db;
    KeszStep[level]=0;
}


//***********************************************************************
void setAnalKeszStep(unsigned level,unsigned n,const char * kovlepesnev){
// h�ny step van k�sz
//***********************************************************************
    if(level>4)throw hiba("setAnalKeszStep","level>4 (%u>4)",level);
//    if(level>0)resetAnalLevel(level-1);
    if(*kovlepesnev!=0)AktFolyamatNev=kovlepesnev;
    KeszStep[level]=n;
    if(level<=2)kiirAnal();
}


//***********************************************************************
unsigned getAnalProgress(unsigned level){
// h�ny sz�zal�kn�l tart
//***********************************************************************
    if(level>4)level=4;
    cuns also = (level==0) ? 0 : getAnalProgress(level-1);
    if(level==1)switch(proti){
        case PT_NoderedElott: return (also+4)/10;
        case PT_Nodered: return 10+also*80/100;
        case PT_NoderedUtan: return 90+(also+4)/10;
    } // k�bgy�k�sen emelkedjen
    return KeszStep[level]==StepDb[level]?100:(100*KeszStep[level]*KeszStep[level]*KeszStep[level]+also*StepDb[level]*StepDb[level])/(StepDb[level]*StepDb[level]*StepDb[level]);
}
//***********************************************************************

