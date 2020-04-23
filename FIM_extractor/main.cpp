//***********************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
//***********************************************************************


//***********************************************************************
#pragma warning(disable : 4996)
#pragma pack(push)
#pragma pack(8)
//***********************************************************************
typedef struct { int kgrid; /* gridsize=2^kgrid       */
		 int kresol;        /* field size =2^kresol   */
		 int nlay;          /* # of the layers        */
		 int analtype;      /* analysis type 0...4    */
		 int nstep;         /* sequence No.           */
		 double t;          /* time value [s]         */
		 double dt;         /* arriving time-step     */
		 double f;          /* frequency [Hz]         */
		 double angle;      /* for timeconst [rad]    */
		 int douim;         /* 0/1=normal/double im   */
		 int sizeofdata;    /* double/complex (8/16)  */
		 } resultheadtype;
//***********************************************************************
#pragma pack(pop)
//***********************************************************************

//***********************************************************************
struct fim_tipus{
//***********************************************************************
    double *re,*im;
    unsigned x,y,z;
    bool komplex;
    fim_tipus():re(NULL),im(NULL),x(0),y(0),z(0),komplex(false){}
    void free(){delete [] re; delete [] im;re=im=NULL;}
};


//***********************************************************************
char * get_fajlnev(char * be){
//***********************************************************************
// az elejérõl és a végérõl leszedi a szóközöket
    size_t n=strlen(be);
    for(size_t i=n-1;isspace(be[i]);i--)be[i]=0;
    for(size_t i=0; i<n; i++)
        if(!isspace(be[i]))return be+i;
    return be;
}

//***********************************************************************
fim_tipus get_fim(char * fimfajl_neve){
//***********************************************************************
    fim_tipus fim;
    FILE * fp=fopen(fimfajl_neve,"rb");
    if(fp==NULL){
        printf("Error: cannot open %s\n",fimfajl_neve);
        exit(1);
    }
	resultheadtype fej;
	fread(&fej,sizeof(fej),1,fp);
	fim.y = fim.x = 1 << fej.kresol;
	fim.z = fej.nlay;
    if(fej.sizeofdata==16)fim.komplex=true;
    else if(fej.sizeofdata!=8){
        printf("Error: fim data size is %d\n",fej.sizeofdata);
        exit(1);
    }
    bool dupla=fej.douim==1;
    fim.re=new double[fim.x*fim.y*fim.z];
    if(fim.komplex)fim.im=new double[fim.x*fim.y*fim.z];
    unsigned xy = fim.x * fim.y;
    unsigned retegmeret= (fim.komplex ? 2 : 1) * (dupla ? xy*2-1 : xy);
    double * reteg=new double[retegmeret];
    for(unsigned z=0; z<fim.z; z++){
        if(fread(reteg,sizeof(double),retegmeret,fp)!=retegmeret){
            printf("Error: cannot read layer %u data\n",z);
            exit(1);
        }
        if(dupla){
            if(fim.komplex){
                for(unsigned j=0; j<fim.y; j++)
                    for(unsigned i=0; i<fim.x; i++){
                        fim.re[xy*z+j*fim.x+i]=reteg[(j*fim.x+i)*4];
                        fim.im[xy*z+j*fim.x+i]=reteg[(j*fim.x+i)*4+1];
                    }
           }
            else{ // valos
                for(unsigned j=0; j<fim.y; j++)
                    for(unsigned i=0; i<fim.x; i++){
                        fim.re[xy*z+j*fim.x+i]=reteg[(j*fim.x+i)*2];
                    }
            }
        }
        else{ // szimpla
            if(fim.komplex){
                for(unsigned i=0; i<xy; i++){
                    fim.re[xy*z+i]=reteg[i*2];
                    fim.im[xy*z+i]=reteg[i*2+1];
                }
           }
            else{ // valos
                for(unsigned i=0; i<xy; i++)
                    fim.re[xy*z+i]=reteg[i];
            }
        }
    }
    fclose(fp);
    delete [] reteg;
    return fim;
}

//***********************************************************************
void kov_sor(char temp[],char s[],FILE*fp,const char fajlnev[]){
//***********************************************************************
    if(fgets(temp,1024,fp)==NULL){
        printf("Error: unexpected end in %s\n",fajlnev);
        exit(1);
    }
    if(sscanf(temp,"%s",s)!=1){
        printf("Error: unexpected end in %s\n",fajlnev);
        exit(1);
    }
}

//***********************************************************************
bool vesszo=false;
//***********************************************************************

//***********************************************************************
const char * double2string(double x,const unsigned eltol=0){
//***********************************************************************
    static char s[1024];
    sprintf(s+eltol,"%g",x);
    if(vesszo)for(int i=0;s[i];i++)if(s[i+eltol]=='.')s[i+eltol]=',';
    return s+eltol;
}


//***********************************************************************
int main(int argc, char ** argv){
//***********************************************************************
    const char * fajlnev = argc>1 ? argv[1] : "fim_extractor_ctrl.txt";
    FILE * fp = fopen(fajlnev,"rt");
    if(fp==NULL){
        printf("Error: cannot open %s\n",fajlnev);
        exit(1);
    }
    char kinev[1024],temp[1024],s[1024];
    if(fgets(temp,1024,fp)==NULL){
        printf("Error: cannot read FIM_FILE= in %s\n",fajlnev);
        exit(1);
    }
    if(sscanf(temp,"%s",s)!=1){
        printf("Error: FIM_FILE= is missing in %s\n",fajlnev);
        exit(1);
    }
    if(strcmpi(s,"FIM_FILE=")!=0){
        printf("Error: FIM_FILE= is missing in %s, %s found\n",fajlnev,s);
        exit(1);
    }
    fim_tipus fim=get_fim(get_fajlnev(temp+9));
    if(fgets(temp,1024,fp)==NULL){
        printf("Error: cannot read OUT_FILE= in %s\n",fajlnev);
        exit(1);
    }
    if(sscanf(temp,"%s",s)!=1){
        printf("Error: OUT_FILE= is missing in %s\n",fajlnev);
        exit(1);
    }
    if(strcmpi(s,"OUT_FILE=")!=0){
        printf("Error: OUT_FILE= is missing in %s, %s found\n",fajlnev,s);
        exit(1);
    }
    strcpy(kinev,get_fajlnev(temp+9));

    kov_sor(temp,s,fp,fajlnev);
    if(strcmpi(s,"USE_COMMAS")==0){
        vesszo=true;
        kov_sor(temp,s,fp,fajlnev);
    }
    remove(kinev);
    while(strcmpi(s,"END")!=0){

        FILE *ki=fopen(kinev,"at");
        if(ki==NULL){
            printf("Error: cannot open %s\n",kinev);
            exit(1);
        }
        printf("%s\n",temp);
        fprintf(ki,"\n%s",temp);
        if(strcmpi(s,"PROBE")==0){
            unsigned x,y,z;
            if(sscanf(temp+6,"%u %u %u",&x,&y,&z)!=3){
                printf("Error: cannot read PROBE x y z koordinates: %s\n",temp);
                exit(1);
            }
            if(fim.komplex) fprintf(ki,"%s\t%s\n",double2string(fim.re[fim.x*fim.y*z+fim.x*y+x]),double2string(fim.im[fim.x*fim.y*z+fim.x*y+x],200));
            else fprintf(ki,"%s\n",double2string(fim.re[fim.x*fim.y*z+fim.x*y+x]));
        }
        else if(strcmpi(s,"X-SECTION")==0){
            unsigned y,z;
            if(sscanf(temp+10,"%u %u",&y,&z)!=2){
                printf("Error: cannot read X-SECTION y z koordinates: %s\n",temp);
                exit(1);
            }
            if(fim.komplex) for(unsigned x=0; x<fim.x; x++)fprintf(ki,"%u\t%s\t%s\n",x,double2string(fim.re[fim.x*fim.y*z+fim.x*y+x]),double2string(fim.im[fim.x*fim.y*z+fim.x*y+x],200));
            else for(unsigned x=0; x<fim.x; x++)fprintf(ki,"%u\t%s\n",x,double2string(fim.re[fim.x*fim.y*z+fim.x*y+x]));
        }
        else if(strcmpi(s,"Y-SECTION")==0){
            unsigned x,z;
            if(sscanf(temp+10,"%u %u",&x,&z)!=2){
                printf("Error: cannot read Y-SECTION x z koordinates: %s\n",temp);
                exit(1);
            }
            if(fim.komplex) for(unsigned y=0; y<fim.y; y++)fprintf(ki,"%u\t%s\t%s\n",y,double2string(fim.re[fim.x*fim.y*z+fim.x*y+x]),double2string(fim.im[fim.x*fim.y*z+fim.x*y+x],200));
            else for(unsigned y=0; y<fim.y; y++)fprintf(ki,"%u\t%s\n",y,double2string(fim.re[fim.x*fim.y*z+fim.x*y+x]));
        }
        else if(strcmpi(s,"Z-SECTION")==0){
            unsigned x,y;
            if(sscanf(temp+10,"%u %u",&x,&y)!=2){
                printf("Error: cannot read Z-SECTION x y koordinates: %s\n",temp);
                exit(1);
            }
            if(fim.komplex) for(unsigned z=0; z<fim.z; z++)fprintf(ki,"%u\t%s\t%s\n",z,double2string(fim.re[fim.x*fim.y*z+fim.x*y+x]),double2string(fim.im[fim.x*fim.y*z+fim.x*y+x],200));
            else for(unsigned z=0; z<fim.z; z++)fprintf(ki,"%u\t%s\n",z,double2string(fim.re[fim.x*fim.y*z+fim.x*y+x]));
        }
        else if(strcmpi(s,"LAYER")==0){
            unsigned z;
            if(sscanf(temp+6,"%u",&z)!=1){
                printf("Error: cannot read LAYER z koordinate: %s\n",temp);
                exit(1);
            }
            if(fim.komplex) 
                for(unsigned x=0; x<fim.x; x++, fprintf(ki,"\n"))
                    for(unsigned y=0; y<fim.y; y++)
                        fprintf(ki,"%s\t%s\t",double2string(fim.re[fim.x*fim.y*z+fim.x*y+x]),double2string(fim.im[fim.x*fim.y*z+fim.x*y+x],200));
            else for(unsigned x=0; x<fim.x; x++, fprintf(ki,"\n"))
                    for(unsigned y=0; y<fim.y; y++)
                        fprintf(ki,"%s\t",double2string(fim.re[fim.x*fim.y*z+fim.x*y+x]));
        }
        fclose(ki);
        kov_sor(temp,s,fp,fajlnev);
    }

    fclose(fp);
    fim.free();
    return 0;
}