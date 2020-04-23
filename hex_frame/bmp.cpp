//***********************************************************************
// Vector SUNRED bitmap class source
// Creation date:	2004. 04. 21.
// Creator:			Pohl László
// Owner:			freeware
//***********************************************************************


//***********************************************************************
#ifndef VSUN_BITMAP_SOURCE
#define	VSUN_BITMAP_SOURCE
//***********************************************************************


//***********************************************************************
#include "bmp.h"
#include "hiba.h"
//***********************************************************************


//***********************************************************************
const RGBQUADx alappaletta[]=
//***********************************************************************
{{0,0,0,0},		{0,0,244,0},	{0,244,0,0},	{0,244,244,0},	{244,0,0,0},	{244,0,244,0},	{244,244,0,0},	{255,255,255,0},
{85,85,85,0},	{113,113,198,0},{113,198,113,0},{56,142,142,0},	{198,113,113,0},{142,56,142,0},	{142,142,56,0},	{44,44,44,0},
{96,44,44,0},	{148,44,44,0},	{200,44,44,0},	{252,44,44,0},	{44,96,44,0},	{96,96,44,0},	{148,96,44,0},	{200,96,44,0},
{252,96,44,0},	{44,148,44,0},	{96,148,44,0},	{148,148,44,0},	{200,148,44,0},	{252,148,44,0},	{44,200,44,0},	{96,200,44,0},
{148,200,44,0},	{200,200,44,0},	{252,200,44,0},	{44,252,44,0},	{96,252,44,0},	{148,252,44,0},	{200,252,44,0},	{252,252,44,0},
{44,44,96,0},	{96,44,96,0},	{148,44,96,0},	{200,44,96,0},	{252,44,96,0},	{44,96,96,0},	{96,96,96,0},	{148,96,96,0},
{200,96,96,0},	{252,96,96,0},	{44,148,96,0},	{96,148,96,0},	{148,148,96,0},	{200,148,96,0},	{252,148,96,0},	{44,200,96,0},
{96,200,96,0},	{148,200,96,0},	{200,200,96,0},	{252,200,96,0},	{44,252,96,0},	{96,252,96,0},	{148,252,96,0},	{200,252,96,0},
{252,252,96,0},	{44,44,148,0},	{96,44,148,0},	{148,44,148,0},	{200,44,148,0},	{252,44,148,0},	{44,96,148,0},	{96,96,148,0},
{148,96,148,0},	{200,96,148,0},	{252,96,148,0},	{44,148,148,0},	{96,148,148,0},	{148,148,148,0},{200,148,148,0},{252,148,148,0},
{44,200,148,0},	{96,200,148,0},	{148,200,148,0},{200,200,148,0},{252,200,148,0},{44,252,148,0},	{96,252,148,0},	{148,252,148,0},
{200,252,148,0},{252,252,148,0},{44,44,200,0},	{96,44,200,0},	{148,44,200,0},	{200,44,200,0},	{252,44,200,0},	{44,96,200,0},
{96,96,200,0},	{148,96,200,0},	{200,96,200,0},	{252,96,200,0},	{44,148,200,0},	{96,148,200,0},	{148,148,200,0},{200,148,200,0},
{252,148,200,0},{44,200,200,0},	{96,200,200,0},	{148,200,200,0},{200,200,200,0},{252,200,200,0},{44,252,200,0},	{96,252,200,0},
{148,252,200,0},{200,252,200,0},{252,252,200,0},{44,44,252,0},	{96,44,252,0},	{148,44,252,0},	{200,44,252,0},	{252,44,252,0},
{44,96,252,0},	{96,96,252,0},	{148,96,252,0},	{200,96,252,0},	{252,96,252,0},	{44,148,252,0},	{96,148,252,0},	{148,148,252,0},
{200,148,252,0},{252,148,252,0},{44,200,252,0},	{96,200,252,0},	{148,200,252,0},{200,200,252,0},{252,200,252,0},{44,252,252,0},
{96,252,252,0},	{148,252,252,0},{200,252,252,0},{255,255,255,0},{0,0,0,0},		{4,4,4,0},		{8,8,8,0},		{12,12,12,0},
{16,16,16,0},	{20,20,20,0},	{24,24,24,0},	{28,28,28,0},	{32,32,32,0},	{36,36,36,0},	{40,40,40,0},	{44,44,44,0},
{48,48,48,0},	{52,52,52,0},	{56,56,56,0},	{60,60,60,0},	{64,64,64,0},	{68,68,68,0},	{72,72,72,0},	{76,76,76,0},
{80,80,80,0},	{84,84,84,0},	{88,88,88,0},	{92,92,92,0},	{96,96,96,0},	{100,100,100,0},{104,104,104,0},{108,108,108,0},
{112,112,112,0},{116,116,116,0},{120,120,120,0},{124,124,124,0},{128,128,128,0},{132,132,132,0},{136,136,136,0},{140,140,140,0},
{144,144,144,0},{148,148,148,0},{152,152,152,0},{156,156,156,0},{160,160,160,0},{164,164,164,0},{168,168,168,0},{172,172,172,0},
{176,176,176,0},{180,180,180,0},{184,184,184,0},{188,188,188,0},{192,192,192,0},{196,196,196,0},{200,200,200,0},{204,204,204,0},
{208,208,208,0},{212,212,212,0},{216,216,216,0},{220,220,220,0},{224,224,224,0},{228,228,228,0},{232,232,232,0},{236,236,236,0},
{240,240,240,0},{244,244,244,0},{248,248,248,0},{252,252,252,0},{255,255,255,0},{90,90,50,0},	{0,0,0,0},		{0,0,0,0},
{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},
{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},
{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},
{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},
{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},
{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0},		{0,0,0,0}};
//***********************************************************************


//***********************************************************************
void bitmap::resize(unsigned y,unsigned x,unsigned BitCount)
//***********************************************************************
{
    const unsigned meret = x * y;
	FH.bfOffBits=BitCount==24?54:54+(1<<BitCount)*4;
	FH.bfReserved1=0;
	FH.bfReserved2=0;
	FH.bfSize=((x*3+3)&(0xfffffffc))*y+54;
	FH.bfType=0x4d42;
	IH.biBitCount=(unsigned short)BitCount;
	IH.biClrImportant=0;
	IH.biClrUsed=0;
	IH.biCompression=0;
	IH.biHeight=y;
	IH.biPlanes=1;
	IH.biSize=40;
	IH.biSizeImage=((x*3+3)&(0xfffffffc))*y;
	IH.biWidth=x;
	IH.biXPelsPerMeter=3800;
	IH.biYPelsPerMeter=3800;
	switch(BitCount){
		case 1:
			FH.bfSize=((x/8+3)&(0xfffffffc))*y+FH.bfOffBits;
			IH.biSizeImage=((x/8+3)&(0xfffffffc))*y;
			break;
		case 4:
			FH.bfSize=((x/2+3)&(0xfffffffc))*y+FH.bfOffBits;
			IH.biSizeImage=((x/2+3)&(0xfffffffc))*y;
			break;
		case 8:
			FH.bfSize=((x+3)&(0xfffffffc))*y+FH.bfOffBits;
			IH.biSizeImage=((x+3)&(0xfffffffc))*y;
			break;
		case 24:
			FH.bfSize=((x*3+3)&(0xfffffffc))*y+FH.bfOffBits;
			IH.biSizeImage=((x*3+3)&(0xfffffffc))*y;
			break;
		default:throw hiba("bitmap::bitmap","unknown bitdepth (%u)",BitCount);
	}
    delete [] m;
    m = NULL;
	if(IH.biSizeImage!=0){
		m=new DWORDpl[meret];
		if(m==0)throw hiba("bitmap::bitmap","allocation fault");
	}
	for(unsigned i=0;i<meret;i++)m[i]=0;
	for(unsigned i=0;i<256;i++)pal[i]=alappaletta[i];
}


//***********************************************************************
void bitmap::kiir()const
//***********************************************************************
{
	printf("x:%d\ny:%d\n\n",IH.biWidth,IH.biHeight);
	for(int i=0;i<IH.biHeight;i++)
	{
		for(int j=0;j<IH.biWidth;j++)printf("%02X ",getpixel_full(j,i));
		printf("\n\t    ");
	}
    printf("\n-------------\n\t    ");
    unsigned db=1,ertek=getpixel_full(0,0);
    for(int i=0;i<IH.biHeight;i++){
        for(int j=(i?0:1);j<IH.biWidth;j++){
            if(ertek==getpixel_full(j,i)/*&&j!=0*/)db++;
            else{
                if(db!=1)printf("%u*%u ",ertek,db);
                else printf("%u ",ertek);
                ertek=getpixel_full(j,i);
                db=1;
//        		if(j==0)printf("\n\t    ");
            }
        }
    }
    if(db!=1)printf("%u*%u ",ertek,db);
    else printf("%u ",ertek);
    printf("\n-------------\n");
}


//***********************************************************************
void bitmap::load(const char *filename)
//***********************************************************************
{
	FILE *fp;
	BYTEpl btemp;
	unsigned  xlength,i,j,k,IT1,IT2;

	FileName=filename;
	if((fp=fopen(filename,"rb"))==0)throw hiba("bitmap::load","cannot open file %s",filename);
	if(fread(&FH,sizeof(FH),1,fp)==0)throw hiba("bitmap::load","data read error");
	if(fread(&IH,sizeof(IH),1,fp)==0)throw hiba("bitmap::load","data read error");
	if(FH.bfType!=0x4d42)throw hiba("bitmap::load","file is not bitmap");
	if(IH.biCompression!=0)throw hiba("bitmap::load","cannot use compressed file");
	if(m!=0)delete [] m;
	m=new DWORDpl[m_size=IH.biWidth*IH.biHeight];
	if(m==0)throw hiba("bitmap::load","allocation fault");
	if(IH.biBitCount<9)
		if(fread(pal,4,1<<IH.biBitCount,fp)==0)
			throw hiba("bitmap::load","data read error");
	fseek(fp,FH.bfOffBits,SEEK_SET);
	switch(IH.biBitCount)
	{
		case  1:xlength=((IH.biWidth+31)>>3)&(0xfffffffc);
				IT1=0;
				for(i=0;i<unsigned(IH.biHeight);i++)
				{
					IT2=0;
					for(j=0;j<xlength;j++)
					{
						if(fread(&btemp,1,1,fp)==0)
							throw hiba("bitmap::load","data read error");
						for(k=0;k<8;k++)
							if(IT2<unsigned(IH.biWidth))
							{
								m[IT1]=btemp/128;
								IT2++;
								IT1++;
								btemp<<=1;
							}
					}
				}
				break;
		case  4:xlength=((IH.biWidth+7)>>1)&(0xfffffffc);
				IT1=0;
				for(i=0;i<unsigned(IH.biHeight);i++)
				{
					IT2=0;
					for(j=0;j<xlength;j++)
					{
						if(fread(&btemp,1,1,fp)==0)
							throw hiba("bitmap::load","data read error");
						if(IT2<unsigned(IH.biWidth))
						{
							m[IT1]=btemp>>4;
							IT2++;
							IT1++;
						}
						if(IT2<unsigned(IH.biWidth))
						{
							m[IT1]=btemp%16;
							IT2++;
							IT1++;
						}
					}
				}
				break;
		case  8:xlength=(IH.biWidth+3)&(0xfffffffc);
				IT1=0;
				for(i=0;i<unsigned(IH.biHeight);i++)
				{
					IT2=0;
					for(j=0;j<xlength;j++)
					{
						if(fread(&btemp,1,1,fp)==0)
							throw hiba("bitmap::load","data read error");
						if(IT2<unsigned(IH.biWidth))
						{
							m[IT1]=btemp;
							IT2++;
							IT1++;
						}
					}
				}
				break;
		case 24:xlength=(IH.biWidth*3+3)&(0xfffffffc);
				IT1=0;
				for(i=0;i<unsigned(IH.biHeight);i++)
				{
					IT2=0;
					for(j=0;j<xlength;j++)
					{
						if(fread(&btemp,1,1,fp)==0)
							throw hiba("bitmap::load","data read error");
						if(IT2<unsigned(IH.biWidth))
						{
							m[IT1]=((DWORDpl)btemp)<<16;
							if(fread(&btemp,1,1,fp)==0)
								throw hiba("bitmap::load","data read error");
							m[IT1]+=((DWORDpl)btemp)<<8;
							if(fread(&btemp,1,1,fp)==0)
								throw hiba("bitmap::load","data read error");
							m[IT1]+=btemp;
							IT2++;
							j+=2;
							IT1++;
						}
					}
				}
				break;
		default:throw hiba("bitmap::load","illegal bit count in the given file");
	}
	fclose(fp);
	return;
}

//***********************************************************************
void bitmap::save(const char *filename)
//***********************************************************************
{
	FILE *fp;
	BYTEpl btemp;
	unsigned  xlength,i,j,k,IT1,IT2;

	FileName=filename;
	if(m==0)throw hiba("bitmap::save","cannot save zero size images");
	if((fp=fopen(filename,"wb"))==0)throw hiba("bitmap::save","cannot open file");
	if(fwrite(&FH,sizeof(FH),1,fp)==0)throw hiba("bitmap::save","data write error");
	if(fwrite(&IH,sizeof(IH),1,fp)==0)throw hiba("bitmap::save","data write error");
	if(IH.biBitCount<9)
		if(fwrite(pal,4,1<<IH.biBitCount,fp)==0)
			throw hiba("bitmap::save","data write error");
	fseek(fp,FH.bfOffBits,SEEK_SET);
	switch(IH.biBitCount)
	{
		case  1:xlength=((IH.biWidth+31)>>3)&(0xfffffffc);
				IT1=0;
				for(i=0;i<unsigned(IH.biHeight);i++)
				{
					IT2=0;
					btemp=0;
					for(j=0;j<xlength;j++)
					{
						for(k=0;k<8;k++)
							if(IT2<unsigned(IH.biWidth))
							{
								btemp>>=1;
								btemp+=(BYTEpl)(m[IT1]<<7);
								IT2++;
								IT1++;
							}
						if(fwrite(&btemp,1,1,fp)==0)
							throw hiba("bitmap::save","data write error");
					}
				}
				break;
		case  4:xlength=((IH.biWidth+7)>>1)&(0xfffffffc);
				IT1=0;
				for(i=0;i<unsigned(IH.biHeight);i++)
				{
					IT2=0;
					btemp=0;
					for(j=0;j<xlength;j++)
					{
						btemp=0;
						if(IT2<unsigned(IH.biWidth))
						{
							btemp=(BYTEpl)(m[IT1]<<4);
							IT2++;
							IT1++;
						}
						if(IT2<unsigned(IH.biWidth))
						{
							btemp+=(BYTEpl)(m[IT1]);
							IT2++;
							IT1++;
						}
						if(fwrite(&btemp,1,1,fp)==0)
							throw hiba("bitmap::save","data write error");
					}
				}
				break;
		case  8:xlength=(IH.biWidth+3)&(0xfffffffc);
				IT1=0;
				for(i=0;i<unsigned(IH.biHeight);i++)
				{
					btemp=0;
					IT2=0;
					for(j=0;j<xlength;j++)
					{
						btemp=0;
						if(IT2<unsigned(IH.biWidth))
						{
							btemp=(BYTEpl)(m[IT1]);
							IT2++;
							IT1++;
						}
						if(fwrite(&btemp,1,1,fp)==0)
							throw hiba("bitmap::save","write read error");
					}
				}
				break;
		case 24:xlength=(IH.biWidth*3+3)&(0xfffffffc);
				IT1=0;
				for(i=0;i<unsigned(IH.biHeight);i++)
				{
					btemp=0;
					IT2=0;
					for(j=0;j<xlength;j++)
					{
						if(IT2<unsigned(IH.biWidth))
						{
							btemp=(BYTEpl)(m[IT1]>>16);
							if(fwrite(&btemp,1,1,fp)==0)
								throw hiba("bitmap::save","data write error");
							btemp=(BYTEpl)((m[IT1]>>8)&255);
							if(fwrite(&btemp,1,1,fp)==0)
								throw hiba("bitmap::save","data write error");
							btemp=(BYTEpl)(m[IT1]&255);
							IT2++;
							j+=2;
							IT1++;
						}
						if(fwrite(&btemp,1,1,fp)==0)
							throw hiba("bitmap::save","data write error");
					}
				}
				break;
		default:throw hiba("bitmap::save","illegal bit count in the given file");
	}
	fclose(fp);
	return;
}


//***********************************************************************
void bitmap::to24bit()
//***********************************************************************
{
	if(IH.biBitCount==24)return;
	IH.biBitCount=24;
	FH.bfOffBits=54;
	const int fut=IH.biHeight*IH.biWidth;
	for(int i=0;i<fut;i++)
		m[i]=(pal[m[i]].rgbRed)+(pal[m[i]].rgbGreen<<8)
			+(pal[m[i]].rgbBlue<<16);
	return;
}


#endif
