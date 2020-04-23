//***********************************************************************
// Vector SUNRED bitmap class header
// Creation date:	2004. 04. 21.
// Creator:			Pohl László
// Owner:			freeware
//***********************************************************************


//***********************************************************************
#ifndef VSUN_BITMAP_HEADER
#define	VSUN_BITMAP_HEADER
//***********************************************************************


//***********************************************************************
class bitmap;
//***********************************************************************


//***********************************************************************
//#ifndef _WINDOWS_
//***********************************************************************


//***********************************************************************
#define BYTEpl				unsigned char
#define WORDpl				unsigned short int
#define DWORDpl				unsigned long int
#define LONGpl				long int
//***********************************************************************


//***********************************************************************
//#endif
//***********************************************************************


//***********************************************************************
#pragma pack(push)
#pragma pack(1)
//***********************************************************************


//***********************************************************************
typedef struct tagBITMAPFILEHEADERx {
//***********************************************************************
        WORDpl    bfType;
        DWORDpl   bfSize;
        WORDpl    bfReserved1;
        WORDpl    bfReserved2;
        DWORDpl   bfOffBits;
} BITMAPFILEHEADERx;


//***********************************************************************
typedef struct tagBITMAPINFOHEADERx {
//***********************************************************************
        DWORDpl      biSize;
        LONGpl       biWidth;
        LONGpl       biHeight;
        WORDpl       biPlanes;
        WORDpl       biBitCount;
        DWORDpl      biCompression;
        DWORDpl      biSizeImage;
        LONGpl       biXPelsPerMeter;
        LONGpl       biYPelsPerMeter;
        DWORDpl      biClrUsed;
        DWORDpl      biClrImportant;
} BITMAPINFOHEADERx;


//***********************************************************************
typedef struct tagRGBQUADx {
//***********************************************************************
        BYTEpl    rgbBlue;
        BYTEpl    rgbGreen;
        BYTEpl    rgbRed;
        BYTEpl    rgbReserved;
} RGBQUADx;


//***********************************************************************
#pragma pack(pop)
//***********************************************************************


//***********************************************************************
#include "PLString.h"
//***********************************************************************


//***********************************************************************
class bitmap
//***********************************************************************
{
	BITMAPFILEHEADERx FH;
	BITMAPINFOHEADERx IH;
	RGBQUADx pal[256];//paletta, mindig annyi eleme van,
							 //ahany bites a szinmelyseg
	DWORDpl *m;
    unsigned m_size; // a kép mérete pixelben
	PLString FileName;
    enum konstansok{also_bit=8,szinszam=1<<also_bit,also_maszk=szinszam-1};
public:
	bitmap(const char *filename):m_size(0){m=0;load(filename);}
    bitmap(unsigned y=0,unsigned x=0,unsigned BitCount=24):m_size(0),m(NULL){resize(y,x,BitCount);}
	~bitmap(){if(m!=NULL)delete [] m;m=NULL;}

	void resize(unsigned y=0,unsigned x=0,unsigned BitCount=24); // ugyanaz, mint a konstruktor
	void load(const char *filename);
	void save(const char *filename);
    void load(const PLString & filename){load(filename.c_str());}
    void save(const PLString & filename){save(filename.c_str());}

	void to24bit();//Convert <24 bit -> 24bit bmp

	const RGBQUADx& getpal(int n)const{return pal[n];}
	unsigned getxsize()const{return IH.biWidth;}
	unsigned getysize()const{return IH.biHeight;}
	unsigned getbitcount()const{return IH.biBitCount;}
    DWORDpl getpixel_also(unsigned x, unsigned y)const{ return m[IH.biWidth*y + x] & also_maszk; }
    DWORDpl getpixel_also(unsigned n)const{ return m[n] & also_maszk; }
    DWORDpl getpixel_felso(unsigned x, unsigned y)const{ return m[IH.biWidth*y + x] >> also_bit; }
    DWORDpl getpixel_felso(unsigned n)const{ return m[n] >> also_bit; }
    DWORDpl getpixel_full(unsigned x, unsigned y)const{ return m[IH.biWidth*y + x]; }
    DWORDpl getpixel_full(unsigned n)const{ return m[n]; }
    void setpixel_also(unsigned x, unsigned y, DWORDpl value){ m[IH.biWidth*y + x] = (m[IH.biWidth*y + x] & ~also_maszk) | (value & also_maszk); }
    void setpixel_felso(unsigned x, unsigned y, DWORDpl value){ m[IH.biWidth*y + x] = (m[IH.biWidth*y + x] & also_maszk) | (value << also_bit); }
    //	DWORDpl & getpixelref(unsigned x,unsigned y){return m[IH.biWidth*y+x];}
//	DWORDpl & getpixelref(unsigned n){return m[n];}
//	DWORDpl * getm(){return m;}
	const PLString & GetFileName(){return FileName;}
    bool find_also(unsigned szin)const{
        for(unsigned i=0; i<m_size; i++)
            if( (m[i] & also_maszk) == szin ) return true;
        return false;
    }
			
	void destroy(){if(m!=NULL)delete [] m;m=NULL;}
	void kiir()const;
};


#endif
