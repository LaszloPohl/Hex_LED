//***********************************************************************
// t�pusok header
// Creation date:  2009. 07. 11.
// Creator:        Pohl L�szl�
//***********************************************************************


//***********************************************************************
#ifndef VSUN3_TIPUSOK_HEADER
#define	VSUN3_TIPUSOK_HEADER
//***********************************************************************
#define vsundebugmode
//***********************************************************************
#pragma warning(disable : 4996)
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include "isspace.h"
#include "PLString.h"
#include "hiba.h"
#include "srfajl.h"
#include "listaestomb.h"
#include "bmp.h"
//***********************************************************************
typedef signed char			i8;
typedef unsigned char		u8;
typedef short				i16;
typedef unsigned short		u16;
typedef int					i32;
typedef unsigned			u32;
typedef const i8			ci8;
typedef const u8			cu8;
typedef const i16			ci16;
typedef const u16			cu16;
typedef const i32			ci32;
typedef const u32			cu32;
typedef double	    		dbl;
typedef const double		cd;
typedef unsigned            uns;
typedef const unsigned      cuns;
//***********************************************************************
enum ProgressTipus{PT_NoderedElott,PT_Nodered,PT_NoderedUtan,PT_Egyenletes};
extern u32 VSUN_CPU_Thread,VSUN_Akt_Thread;
extern bool ConsoleText,ConDebug,is_simple_anal_progress;
extern bool is_transi_step_as_dc;
extern ProgressTipus proti;
extern char hibaUzenet[1024];
extern char idobelyeg[100];
extern FILE *logfile;
//***********************************************************************
cd nulla=0.0;
cd egy=1.0;
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
cuns colmax=256;
cd g_max=1.0e+030;
cd g_min=egy/g_max;
cd absT=273.15;
cd szigma=5.6704e-8;
cd his_megengedett_arany=0.1;
cd his_szamitott_arany=0.5*his_megengedett_arany;
//cd xUt=0.026;
enum Oldal{	EXTERNAL=0,WEST=1,EAST=2,SOUTH=3,NORTH=4,BOTTOM=5,TOP=6,//Basic oldalak //Tilos megv�ltoztatni a sorrendet!!
			WS=7,WN=8,WB=9,WT=10,ES=11,EN=12,EB=13,ET=14,SB=15,ST=16,NB=17,NT=18,CENTER=19};// Extended oldalak //Tilos megv�ltoztatni a sorrendet!!
extern const char *OldalNev[];
cuns BASIC_SIDES=        7; //oldalak sz�ma (kocka 6 oldala+externals)
cuns EXTENDED_SIDES=	20; //oldalak sz�ma (kocka 6 oldala+externals+12 �lk�z�pi "oldal")(A CENTER-t is belesz�m�tom, 0. l�p�sben kell!)
cu32 NEWKIVALTO=		32;
enum monlintipus{nlt_lin,nlt_linearis,nlt_exp,nlt_diode,nlt_quadratic,nlt_erno,nlt_szakaszok,nlt_mizs1,nlt_lsfit};
enum lsfit_veg_tipus { lsv_none, lsv_lin, lsv_strong };
enum lsfit_egyenlet_tipus{ lse_polinom };
enum MezoTipus {FieldEl,FieldTherm,FieldElTherm}; // FieldEl,FieldTherm,FieldElTherm
enum SzinTipus{SzinNormal,SzinBoundary,SzinUchannel};
enum CsatornaTipus{CsatRect, CsatCirc, CsatTrap};
enum ConvTipus {ConvHTC,ConvUpper_1,ConvLower_1,ConvVertical_1,ConvVerticalChurchill_P_1,ConvVerticalChurchill_P_2,
                ConvVerticalChurchill_T,ConvVerticalChurchill_C,ConvVerticalLee_T,ConvVerticalLee_P,ConvVerticalMihajev,
                ConvYovanovichMin,ConvYovanovichMax,ConvWei,ConvEdgeWei_H,ConvEdgeWei_I,ConvEdgeWei_HI,ConvEdgeWei_HH};
enum PeremTipus {PeremOpen, PeremV, PeremT, PeremR,PeremRU}; // RU=nem ambient h�m�rs�klethez konvekci�
enum GerjTipus {GerjSemmi,GerjU,GerjI};
enum ProbeTipus {ProbeV,ProbeT,ProbeC,ProbeM,ProbeS};// M=map, S=section
enum AnalizisTipus {AnalDC,AnalNDC,AnalAC,AnalLinTran,AnalLogTran,AnalBode,AnalIdo,AnalCtrl};
enum irany{X_IRANY,Y_IRANY,Z_IRANY,N_IRANY};//N_IRANY=nincs ir�ny
enum CompactTipus{CompactMatrix};

struct stackedMaterial {
    std::string name;
    double gTh, cTh; // gTh: fajl. h�vezet�s, Si eset�n 156.3, cTh: fajl. h�kapac, Si eset�n 1.596e+006
};
//***********************************************************************
struct Vec3d{
    dbl x, y, z;
    Vec3d() :x(nulla), y(nulla), z(nulla){}
    Vec3d(dbl X, dbl Y, dbl Z) :x(X), y(Y), z(Z){}
    dbl hossz()const { return sqrt(x*x + y*y + z*z); }
    dbl get_vertical_angle()const { return fabs(acos(z / sqrt(x*x + y*y + z*z))); } // A z ir�ny� egys�gvektorral bez�rt sz�g
    Vec3d operator+(const Vec3d &b)const{ return Vec3d(x + b.x, y + b.y, z + b.z); }
    Vec3d operator-(const Vec3d &b)const{ return Vec3d(x - b.x, y - b.y, z - b.z); }
    dbl   operator*(const Vec3d &b)const{ return x*b.x + y*b.y + z*b.z; }
    Vec3d operator*(dbl c)const { return Vec3d{ x*c,y*c,z*c }; }
    friend Vec3d operator*(dbl c, const Vec3d & v) { return v*c; }
    friend dbl Abs(const Vec3d &v){ return sqrt(v*v); }
    friend void swap(Vec3d &a, Vec3d &b){ dbl t; t = a.x; a.x = b.x; b.x = t; t = a.y; a.y = b.y; b.y = t; t = a.z; a.z = b.z; b.z = t; }
    inline friend dbl vegyesszorzat(const Vec3d &a, const Vec3d &b, const Vec3d &c) {
        return a.x*(b.y*c.z - b.z*c.y) - a.y*(b.x*c.z - b.z*c.x) + a.z*(b.x*c.y - b.y*c.x);
    }
};
inline dbl terszog(const Vec3d &a, const Vec3d &b, const Vec3d &c) {
    dbl ah = a.hossz();
    dbl bh = b.hossz();
    dbl ch = c.hossz();
    dbl vsz = vegyesszorzat(a, b, c);
    dbl tag1 = ah*bh*ch;
    dbl tag2 = ah*(b*c);
    dbl tag3 = bh*(a*c);
    dbl tag4 = ch*(a*b);
    return 2*fabs(atan2(vsz, tag1 + tag2 + tag3 + tag4));
}
//***********************************************************************
struct tol_db { // kezd�index �s innen kezdve h�ny darab
    uns tol, db; 
    tol_db() :tol{ 0 }, db{ 0 } {} 
};
struct mit_hova_masol {
    bool is_be1;  // true: a be1-r�l kell m�solni, false: a be2-r�l kell m�solni
    uns honnan, hova, hanyat; // h�ny csom�pont van az adott oldalon, �s az hol kezd�dik a forr�s �s a c�l m�trixban
    mit_hova_masol() :is_be1{ false }, honnan{ 0 }, hova{ 0 }, hanyat{ 0 } {}
};

//***********************************************************************
void most(const std::string & mi_tortent);
void esemenyek_kiirasa();
//***********************************************************************

#endif