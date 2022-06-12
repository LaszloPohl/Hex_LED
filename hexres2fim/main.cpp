#include <stdio.h>
#include <stdarg.h>

#include "PLString.h"
typedef unsigned            uns;
typedef unsigned            u32;
typedef const unsigned      cuns;
typedef double              dbl;


//****************************************************************
template<class C> class tomb3d{
//****************************************************************
    C * m;
    unsigned x,y,z,n;
public:
    //****************************************************************
    tomb3d(){m=NULL;x=y=z=n=0;}
    ~tomb3d(){if(m)delete [] m;}

    unsigned x_size()const{return x;}
    unsigned y_size()const{return y;}
    unsigned z_size()const{return z;}
    unsigned   size()const{return n;}

    //****************************************************************


    //****************************************************************
    C & getref(unsigned X,unsigned Y,unsigned Z){
    //****************************************************************
#ifdef vsundebugmode
        if(X>=x||Y>=y||Z>=z){
            throw hiba("tomb3d::getref","X>=x||Y>=y||Z>=z (%u>=%u||%u>=%u||%u>=%u)",X,x,Y,y,Z,z);
        }
#endif
        return m[Z*(x*y)+Y*x+X];
    }


    //****************************************************************
    C & getref(unsigned N){
    //****************************************************************
#ifdef vsundebugmode
        if(N>=n)throw hiba("tomb3d::getref","N>n (%u>=%u)",N,n);
#endif
        return m[N];
    }


    //****************************************************************
    const C & getconstref(unsigned X,unsigned Y,unsigned Z)const{
    //****************************************************************
#ifdef vsundebugmode
        if(X>=x||Y>=y||Z>=z)throw hiba("tomb3d::getconstref","X>=x||Y>=y||Z>=z (%u>=%u||%u>=%u||%u>=%u)",X,x,Y,y,Z,z);
#endif
        return m[Z*(x*y)+Y*x+X];
    }


     //****************************************************************
    const C & getconstref(unsigned N)const{
    //****************************************************************
#ifdef vsundebugmode
        if(N>=n)throw hiba("tomb3d::getref","N>n (%u>=%u)",N,n);
#endif
        return m[N];
    }


   //****************************************************************
    C get(unsigned X,unsigned Y,unsigned Z){
    //****************************************************************
#ifdef vsundebugmode
        if(X>=x||Y>=y||Z>=z)throw hiba("tomb3d::get","X>=x||Y>=y||Z>=z (%u>=%u||%u>=%u||%u>=%u)",X,x,Y,y,Z,z);
#endif
        return m[Z*(x*y)+Y*x+X];
    }


    //****************************************************************
    tomb3d(const tomb3d & t){
    //****************************************************************
        m=new C[n=t.n];
        if(m==0)throw hiba("tomb::tomb(tomb&)","alloc failed");
        x=t.x; y=t.y; z=t.z;
        for(unsigned i=0;i<n;i++)m[i]=t.m[i];
    }


    //****************************************************************
    void resize(unsigned new_x,unsigned new_y,unsigned new_z){
    //****************************************************************
        unsigned n2=new_x*new_y*new_z;
        x=new_x;
        y=new_y;
        z=new_z;
        if(n==n2&&m!=NULL)return;
        delete [] m;
        m=new C[n=n2];
        if(m==NULL)throw hiba("tomb3d::resize","alloc failed");
    }

    //****************************************************************
    tomb3d & operator=(const tomb3d & a){
    //****************************************************************
        if(this!=&a){
            if(m)delete [] m;
            m=new C[n=(x=a.x)*(y=a.y)*(z=a.z)];
            if(m==0)throw hiba("tomb3d::operator=","alloc failed");
            for(unsigned i=0;i<n;i++)m[i]=a.m[i];
        }
        return *this;
    }
};


//****************************************************************
template<class C> class tomb
//****************************************************************
{
    C * m;
    unsigned n,findindex;
public:
    //****************************************************************
    tomb():n(0),m(NULL),findindex(0){}
    ~tomb(){delete [] m;m=NULL;}

    unsigned size()const{return n;}
    unsigned getfindindex()const{return findindex;}
    void clear(){delete [] m; m=NULL;n=findindex=0;}
    void swap(tomb<C> &masik) { 
        C * ptemp = m; m = masik.m; masik.m = ptemp;
        unsigned utemp = n; n = masik.n; masik.n = utemp;
        utemp = findindex; findindex = masik.findindex; masik.findindex = utemp;
    }
    void zero() {
        for (unsigned i = 0; i < n; i++)
            m[i] = C();
    }
    //****************************************************************


    //****************************************************************
    tomb(const tomb & t){
    //****************************************************************
        m=new C[n=t.n];
        if(m==0)throw hiba("tomb::tomb(tomb&)","alloc failed");
        for(unsigned i=0;i<n;i++)m[i]=t.m[i];
    }


    //****************************************************************
    C * add(C elem,unsigned db=1){
    // ha többet tesz be, akkor az ezek közül utolsó címét adja vissza
    //****************************************************************
        C *m2=new C[n+db];
        if(m2==0)throw hiba("tomb::add","alloc failed");
        if(m){for(unsigned i=0;i<n;i++)m2[i]=m[i];delete [] m;}
        while(db--)m2[n++]=elem;
        m=m2;
        return &m2[n-1];
    }


    //****************************************************************
    void add(const tomb<C> &t){
    //****************************************************************
        C *m2=new C[n+t.n];
        if(m2==0)throw hiba("tomb::add","alloc failed");
        if(m){for(unsigned i=0;i<n;i++)m2[i]=m[i];delete [] m;}
        m = m2;
        for(unsigned i=0; i<t.n; i++) m[i+n] = t.m[i];
        n += t.n;
    }


    //****************************************************************
    void puffer_add(unsigned index, C elem){
    // Ha a tömb után eggyel akarunk új elemet betenni, akkor 100
    // új elemnek foglal helyet
    //****************************************************************
        if(index==n)resize(n+100);
        (*this)[index] = elem;
    }


    //****************************************************************
    void resize(unsigned newsize){
    //****************************************************************
        if(n==newsize)return;
        C *m2=new C[newsize];
        if(m2==0)throw hiba("tomb::resize","alloc failed");
        unsigned masol=(newsize<n)?newsize:n;
        if(m){for(unsigned i=0;i<masol;i++)m2[i]=m[i];delete [] m;}
        n=newsize;
        m=m2;
    }


    //****************************************************************
    void decsize(){// eggyel csökkenti a tömb méretét, gyorsabb, mint a resize
    //****************************************************************
        if(n)n--;
        else throw hiba("tomb::decsize","n==0");
    }


    //****************************************************************
    unsigned incsize(){resize(n+1); return n;}
    //****************************************************************


    //****************************************************************
    C & getLast(){if(n>0)return m[n-1];else throw hiba("tomb::getLast","n==0");}
    const C & getLast()const{if(n>0)return m[n-1];else throw hiba("tomb::getLast","n==0");}
    //****************************************************************


    //****************************************************************
    C & operator[](unsigned index){
    //****************************************************************
//#ifdef vsundebugmode
        if(index>=n)
            throw hiba("tomb::operator[]","index>=n (%u>=%u)",index,n);
//#endif
        return m[index];
    }


    //****************************************************************
    const C & operator[](unsigned index)const{
    //****************************************************************
//#ifdef vsundebugmode
        if(index>=n)
            throw hiba("tomb::operator[]","index>=n (%u>=%u)",index,n);
//#endif
        return m[index];
    }


    //****************************************************************
    tomb & operator=(const tomb & a){
    //****************************************************************
        if(this!=&a){
            if(m)delete [] m;
            m=new C[n=a.n];
            if(m==0)throw hiba("tomb::operator=","alloc failed");
            for(unsigned i=0;i<n;i++)m[i]=a.m[i];
        }
        return *this;
    }
};


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
inline uns lk2hatvany(uns ertek){
//***********************************************************************
    uns n=1,step=1u;
    for(;step<<=1;n++);
    for(u32 i=0;i<n;i++)if(ertek<=(1u<<i))return i;
    throw hiba("legkisebb2hatvanyamibebelefer()","impossibility");
}


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
struct cella_adat {
//***********************************************************************
    bool is_el, is_th, is_ph;
    uns c, x, y, z, f;
    cella_adat() :is_el{ false }, is_th{ false }, is_ph{ false }, c { 0 }, x{ 0 }, y{ 0 }, z{ 0 }, f{ 0 } {}
};


//***********************************************************************
struct fimdesc {
//***********************************************************************
    uns cellnum, x, y, z;
    uns FIM_res_xy, FIM_res_z, FIM_diff_x, FIM_diff_y, FIM_diff_z;
    uns anal_db;
    tomb<cella_adat> cellak;
    void read_from_file(const PLString & fajlnev);
};


//***********************************************************************
struct maphead {
//***********************************************************************
    char azon[8];
    uns adatmeret; // bájtban
    uns tomorites_tipus; // 0 = nyers
    maphead() :azon{ "hexmap" }, adatmeret{ 0 }, tomorites_tipus{ 0 } {}
};


//***********************************************************************
struct hexres {
//***********************************************************************
    tomb<float> dc_ertekek;
    bool is_face_data;
    bool read(const PLString & fajlnev, bool is_face_dat);
};


//***********************************************************************
void get_sim_names(const PLString & hex_file, tomb<PLString> & sim_names) {
//***********************************************************************
    char tomb[1024];
    FILE *fp;
    if (fopen_s(&fp, hex_file.c_str(), "rt") != 0)
        throw hiba("get_sim_names", "cannot open %s to read", hex_file.c_str());
    int sim_db = -1;
    while (sim_db < 0 && fgets(tomb, 1024, fp) != nullptr) {
        PLString elso_szo = PLString(tomb).get_first_word().UpCase();
        if (elso_szo[0] == 'N' && elso_szo[1] == 'S')
            if (!elso_szo.toint(sim_db, 2))
                throw hiba("get_sim_names", "cannot convert NS? to number (%s)", elso_szo.c_str());
    }
    sim_names.resize(sim_db);
    int akt_sim_index = 0;
    while (fgets(tomb, 1024, fp) != nullptr && akt_sim_index < sim_db) {
        PLString sor(tomb);
        PLString elso_szo = sor.get_first_word().UpCase();
        if (elso_szo[0] == 'B' && elso_szo[1] == 'S') {
            int eleje = sor.find('\"');
            int vege = sor.findr('\"');
            if (eleje == vege)
                throw hiba("get_sim_names", "unknown BS line format (%s)", tomb);
            sim_names[akt_sim_index] = sor.substr(eleje + 1, vege - eleje - 1);
            akt_sim_index++;
        }
    }
    fclose(fp);
    //sim_names[0] = "Marci_LED_modell_190401_152718";
    if(akt_sim_index < sim_db)
        throw hiba("get_sim_names", "not enough simulation (%d<%d)", akt_sim_index, sim_db);
}


//***********************************************************************
struct sim_adat {
//***********************************************************************
    uns FIM_res_xy;
    uns FIM_res_z;
    uns FIM_diff_x;
    uns FIM_diff_y;
    uns FIM_diff_z;
    tomb3d<float> map;
};


//***********************************************************************
struct sim_adatok {
//***********************************************************************
    sim_adat el, th, rad_a, rad_b, lum_a, lum_b, diss_a, diss_b;
    sim_adat blue_a, blue_b, yllw_a, yllw_b, diss_c, diss_d;
    bool is_el, is_th;
    void build(const fimdesc & desc, const hexres & res);
    void build_rad_lum(const fimdesc & desc, const hexres & res);
};


//***********************************************************************
void write_fim(const PLString & fim_file_name, const sim_adat & sim){
//***********************************************************************
    cuns x = sim.map.x_size(), y = sim.map.y_size(), z = sim.map.z_size();

    uns wxy=(1u<<lk2hatvany(x)>1u<<lk2hatvany(y))?1u<<lk2hatvany(x):1u<<lk2hatvany(y);
    uns wz=1u<<lk2hatvany(z);

    if (sim.FIM_res_xy > wxy)
        wxy = sim.FIM_res_xy;
    if (sim.FIM_res_z > wz)
        wz = sim.FIM_res_z;

    cuns dx = sim.FIM_diff_x + x > wxy ? 0 : sim.FIM_diff_x;
    cuns dy = sim.FIM_diff_y + y > wxy ? 0 : sim.FIM_diff_y;
    cuns dz = sim.FIM_diff_z + z > wz  ? 0 : sim.FIM_diff_z;

    resultheadtype fej;
    memset(&fej,0,sizeof(fej));

    fej.kgrid=lk2hatvany(wxy*wz);
    fej.kresol=lk2hatvany(wxy);
    fej.nlay=wz;

    fej.analtype = 0; // 0: dc, 1: ac
    fej.sizeofdata = sizeof(double); // ac: 2*

    fej.nstep = 1;
    fej.t=0;
    fej.dt = 0;
    fej.f = 0;// par.is_ac? par.f : 0.0;
    fej.angle=0.0;
    fej.douim=0;

    // fájl írása

    if(false){ //
        /*
        const tomb3d<dcomplex_cell_res> * ptemp=t_dcmap[is_el?par.u_azon:par.T_azon];
        if(!sim.is_nofim){
            cuns db=2*wxy*wxy*wz+wxy*wxy;
            cuns Ldb=wxy*wxy;
            cuns Hdb=Ldb*wz;
            uns N=0;
            dcomplex * tor=new dcomplex[db];
            dcomplex * fim2=tor+Hdb;
            dcomplex * fim3=fim2+Hdb;

            for (uns i = 0;i < db;i++)tor[i] = nulla;
            for(uns i=0;i<z;i++){
                uns layelt=Ldb*(i+dz);
                for(u32 j=0;j<y;j++){
                    uns sorelt=wxy*(j+dy);
                    for(uns k=0;k<x;k++,N++){
                        tor[layelt+sorelt+k+dx] =ptemp->getconstref(N).t[EXTERNAL];
                        fim2[layelt+sorelt+k+dx]=ptemp->getconstref(N).t[BOTTOM];
                        fim3[sorelt+k+dx]       =ptemp->getconstref(N).t[TOP]; // csak az utolsó rétegnél kellene
                    }
                }
            }

            PLString FajlNev=nev+"SRE"+t+(is_el?".FIMe":".FIM");
            if((fp=fopen(FajlNev.c_str(),"wb"))==0)throw hiba("masina::writefim()","cannot open %s",FajlNev.c_str());
            fwrite(&fej,sizeof(fej),1,fp);
            fwrite(tor,sizeof(dcomplex),Hdb+Hdb+Ldb,fp);
            fclose(fp);

            delete [] tor;
        }
        */
    }
    else {
        uns N = 0;
        cuns db = 2 * wxy*wxy*wz + wxy*wxy;
        cuns Ldb = wxy*wxy;
        cuns Hdb = Ldb*wz;
        cuns db_i = 8 * Hdb;
        dbl * tor = new dbl[db];
        dbl * fim2 = tor + Hdb;
        dbl * fim3 = fim2 + Hdb;
        for (uns i = 0;i < db;i++)tor[i] = sim.map.getconstref(0);
        for (uns i = 0;i < z;i++) {
            uns layelt = Ldb*(i + dz);
            for (u32 j = 0;j < y;j++) {
                uns sorelt = wxy*(j + dy);
                for (uns k = 0;k < x;k++, N++) {
                    tor[layelt + sorelt + k + dx] = sim.map.getconstref(N); // center
                    fim2[layelt + sorelt + k + dx] = sim.map.getconstref(N); // bottom
                    fim3[sorelt + k + dx] = sim.map.getconstref(N);  // top, csak az utolsó rétegnél kellene
                }
            }
        }

        // FIM

        FILE *fp;
        if ((fp = fopen(fim_file_name.c_str(), "wb")) == 0)throw hiba("masina::writefim()", "cannot open %s", fim_file_name.c_str());
        fwrite(&fej, sizeof(fej), 1, fp);
        fwrite(tor, sizeof(double), Hdb + Hdb + Ldb, fp);
        fclose(fp);
        printf("File created: %s\n", fim_file_name.c_str());

        delete[] tor;
    }
}


//***********************************************************************
int main(int par_db, char ** pars) {
//***********************************************************************
    try {
        tomb<PLString> sim_names;
        // kiolvassa a hex fájlból a szimulációk neveit
        if (par_db < 2)
            throw hiba("main", "missing command line argument: hex file name");
        get_sim_names(pars[1], sim_names);
        PLString path = getPath(pars[1]);
        fimdesc desc;
        hexres res;
        sim_adatok adat;
        for (uns i = 0; i < sim_names.size(); i++) {
            desc.read_from_file(path + sim_names[i] + ".fimdesc");
            for (uns j = 1; j <= desc.anal_db; j++) {
                printf("**** %u ****\n", j);
                if (res.read(path + "hex_res/" + sim_names[i] + "/hexres" + j + ".cvt", false)) {
                    adat.build(desc, res);
                    char szam[100];
                    sprintf(szam, "%04u", j);
                    if (adat.is_el)
                        write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fime", adat.el);
                    if (adat.is_th)
                        write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim", adat.th);
                }
                printf("**** %u ****\n", j);
                if (res.read(path + "hex_res/" + sim_names[i] + "/hexres" + j + ".crl", false)) {
                    adat.build_rad_lum(desc, res);
                    char szam[100];
                    sprintf(szam, "%04u", j);
                    write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim_rada", adat.rad_a);
                    write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim_radb", adat.rad_b);
                    write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim_luma", adat.lum_a);
                    write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim_lumb", adat.lum_b);
                    write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim_dissa", adat.diss_a);
                    write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim_dissb", adat.diss_b);
                    write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim_bluea", adat.blue_a);
                    write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim_blueb", adat.blue_b);
                    write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim_yllwa", adat.yllw_a);
                    write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim_yllwb", adat.yllw_b);
                    write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim_dissc", adat.diss_c);
                    write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + ".fim_dissd", adat.diss_d);
                }
                printf("**** %u ****\n", j);
                if (res.read(path + "hex_res/" + sim_names[i] + "/hexres" + j + ".fvt", true)) {
                    adat.build(desc, res);
                    char szam[100];
                    sprintf(szam, "%04u", j);
                    if (adat.is_el)
                        write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + "f.fime", adat.el);
                    if (adat.is_th)
                        write_fim(path + "hex_res/" + sim_names[i] + "/sre" + szam + "f.fim", adat.th);
                }
                // a rad-lum-f még TODO, a build_rad_lum-ot megfelelõen kell elkészíteni
            }
        }
    }
    catch (hiba h) {
        printf("\n%s\n", h.what());
    }
}


//***********************************************************************
void fimdesc::read_from_file(const PLString & fajlnev) {
//***********************************************************************
    cellak.clear();
    FILE *fp;
    if (fopen_s(&fp, fajlnev.c_str(), "rt") != 0)
        throw hiba("fimdesc::read_from_file", "cannot open %s to read", fajlnev.c_str());
    char tomb[1024];

    // 1. sor

    if (fgets(tomb, 1024, fp) == nullptr)
        throw hiba("fimdesc::read_from_file", "cannot read 1st line of %s", fajlnev.c_str());
    if(sscanf(tomb, "%u%u%u%u", &cellnum, &x, &y, &z)!=4)
        throw hiba("fimdesc::read_from_file", "no 4 data in 1st line of %s (%s)", fajlnev.c_str(), tomb);

    // 2. sor

    if (fgets(tomb, 1024, fp) == nullptr)
        throw hiba("fimdesc::read_from_file", "cannot read 2nd line of %s", fajlnev.c_str());
    if (sscanf(tomb, "%u%u%u%u%u", &FIM_res_xy, &FIM_res_z, &FIM_diff_x, &FIM_diff_y, &FIM_diff_z) != 5)
        throw hiba("fimdesc::read_from_file", "no 5 data in 2nd line of %s (%s)", fajlnev.c_str(), tomb);

    // cellasorok

    cellak.resize(cellnum);
    for (uns i = 0; i < cellnum; i++) {
        char terek[1024];
        if (fgets(tomb, 1024, fp) == nullptr)
            throw hiba("fimdesc::read_from_file", "cannot read cell %u line from %s", i+1, fajlnev.c_str());
        if (sscanf(tomb, "%s%u%u%u%u%u", terek, &cellak[i].c, &cellak[i].x, &cellak[i].y, &cellak[i].z, &cellak[i].f) != 6)
            throw hiba("fimdesc::read_from_file", "no 6 data in cell %u line of %s (%s)", i+1, fajlnev.c_str(), tomb);
        cellak[i].is_el = terek[0] == 'E' || terek[1] == 'E';
        cellak[i].is_ph = terek[0] == 'P' || terek[1] == 'P';
        cellak[i].is_th = terek[0] == 'T' || terek[1] == 'T';
    }

    // analízisszám

    if (fgets(tomb, 1024, fp) == nullptr)
        throw hiba("fimdesc::read_from_file", "cannot read last line of %s", fajlnev.c_str());
    if (sscanf(tomb, "%u", &anal_db) != 1)
        throw hiba("fimdesc::read_from_file", "no analysis number in last line of %s (%s)", fajlnev.c_str(), tomb);
}


//***********************************************************************
bool hexres::read(const PLString & fajlnev, bool is_face_dat) {
//***********************************************************************
    char fej[128] = { 0 };
    maphead & headstrukt = *((maphead*)fej);
    FILE *fp;
    is_face_data = is_face_dat;
    if (fopen_s(&fp, fajlnev.c_str(), "rb") != 0)
        return false;
    fread(fej, 1, 128, fp);
    if(PLString(headstrukt.azon)!="hexmap")
        throw hiba("hexres::read", "%s is not a hexres file", fajlnev.c_str());
    if (headstrukt.tomorites_tipus != 0)
        throw hiba("hexres::read", "compression type is not 0 in %s", fajlnev.c_str());
    dc_ertekek.clear();
    dc_ertekek.resize(headstrukt.adatmeret / sizeof(float));
    if (fread(&dc_ertekek[0], sizeof(float), dc_ertekek.size(), fp) != dc_ertekek.size())
        throw hiba("hexres::read", "file is corrupt (%s)", fajlnev.c_str());
    fclose(fp);
    return true;
}


//***********************************************************************
void sim_adatok::build(const fimdesc & desc, const hexres & res) {
//***********************************************************************
    is_el = is_th = false;
    el.FIM_res_xy = th.FIM_res_xy = desc.FIM_res_xy;
    el.FIM_res_z  = th.FIM_res_z  = desc.FIM_res_z;
    el.FIM_diff_x = th.FIM_diff_x = desc.FIM_diff_x;
    el.FIM_diff_y = th.FIM_diff_y = desc.FIM_diff_y;
    el.FIM_diff_z = el.FIM_diff_z = desc.FIM_diff_z;
    el.map.resize(desc.x, desc.y, desc.z);
    th.map.resize(desc.x, desc.y, desc.z);
    float el_filler[16] = { 0 }, th_filler[16] = { 0 };
    uns db = 0;
    uns ertek_index = 0;
    bool first_el = true, first_th = true;
    for (uns i = 0; i < desc.cellnum; i++) {
        if (desc.cellak[i].z>15) {
            printf("Error: max 16 layers are supported. z=%u found\n", desc.cellak[i].z);
            return;
        }
        if (desc.cellak[i].is_el) {
            if (is_el) {
                if (first_el) {
                    first_el = false;
                    for (uns j = 0; j < 16; j++)
                        el_filler[j] = res.dc_ertekek[ertek_index];
                }
                else el_filler[desc.cellak[i].z] = res.dc_ertekek[ertek_index];
            }
            is_el = true;
            ertek_index++;
        }
        if (desc.cellak[i].is_th) {
            if (is_th) {
                if (first_th) {
                    first_th = false;
                    for (uns j = 0; j < 16; j++)
                        th_filler[j] = res.dc_ertekek[ertek_index]; //40;
                    //th_filler[15] = 60;
                }
                else th_filler[desc.cellak[i].z] = res.dc_ertekek[ertek_index];
            }
            is_th = true;
            ertek_index++;
        }
        if (res.is_face_data)
            ertek_index += desc.cellak[i].f; // átugorja a face adatokat
    }
    for (uns i = 0; i < desc.z; i++)
        for (uns j = 0; j < desc.y; j++)
            for (uns k = 0; k < desc.x; k++) {
                el.map.getref(k, j, i) = 0* /*0;*/el_filler[i];
                th.map.getref(k, j, i) = /*i>9 || i<8 ? (i==15?60:40) : */th_filler[i];
            }
    ertek_index = 0;
    for (uns i = 0; i < desc.cellnum; i++) {
        if (desc.cellak[i].is_el) {
            el.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
            is_el = true;
        }
        if (desc.cellak[i].is_th) {
            float T = res.dc_ertekek[ertek_index++];
            //T = T < 40 ? 40 : T;
            //if (desc.cellak[i].z < 6)T = 32;
            th.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = /*desc.cellak[i].z>9|| desc.cellak[i].z<8 ? 40 : */T;
            is_th = true;
        }
        if(res.is_face_data)
            ertek_index += desc.cellak[i].f; // átugorja a face adatokat
    }
}



//***********************************************************************
void sim_adatok::build_rad_lum(const fimdesc & desc, const hexres & res) {
//***********************************************************************
    rad_a.FIM_res_xy = rad_b.FIM_res_xy = lum_a.FIM_res_xy = lum_b.FIM_res_xy = diss_a.FIM_res_xy = diss_b.FIM_res_xy = desc.FIM_res_xy;
    rad_a.FIM_res_z  = rad_b.FIM_res_z  = lum_a.FIM_res_z  = lum_b.FIM_res_z  = diss_a.FIM_res_z  = diss_b.FIM_res_z  = desc.FIM_res_z;
    rad_a.FIM_diff_x = rad_b.FIM_diff_x = lum_a.FIM_diff_x = lum_b.FIM_diff_x = diss_a.FIM_diff_x = diss_b.FIM_diff_x = desc.FIM_diff_x;
    rad_a.FIM_diff_y = rad_b.FIM_diff_y = lum_a.FIM_diff_y = lum_b.FIM_diff_y = diss_a.FIM_diff_y = diss_b.FIM_diff_y = desc.FIM_diff_y;
    rad_a.FIM_diff_z = rad_b.FIM_diff_z = lum_a.FIM_diff_z = lum_b.FIM_diff_z = diss_a.FIM_diff_z = diss_b.FIM_diff_z = desc.FIM_diff_z;
    blue_a.FIM_res_xy = blue_b.FIM_res_xy = yllw_a.FIM_res_xy = yllw_b.FIM_res_xy = diss_c.FIM_res_xy = diss_d.FIM_res_xy = desc.FIM_res_xy;
    blue_a.FIM_res_z  = blue_b.FIM_res_z  = yllw_a.FIM_res_z  = yllw_b.FIM_res_z  = diss_c.FIM_res_z  = diss_d.FIM_res_z  = desc.FIM_res_z;
    blue_a.FIM_diff_x = blue_b.FIM_diff_x = yllw_a.FIM_diff_x = yllw_b.FIM_diff_x = diss_c.FIM_diff_x = diss_d.FIM_diff_x = desc.FIM_diff_x;
    blue_a.FIM_diff_y = blue_b.FIM_diff_y = yllw_a.FIM_diff_y = yllw_b.FIM_diff_y = diss_c.FIM_diff_y = diss_d.FIM_diff_y = desc.FIM_diff_y;
    blue_a.FIM_diff_z = blue_b.FIM_diff_z = yllw_a.FIM_diff_z = yllw_b.FIM_diff_z = diss_c.FIM_diff_z = diss_d.FIM_diff_z = desc.FIM_diff_z;
    rad_a.map.resize(desc.x, desc.y, desc.z);
    rad_b.map.resize(desc.x, desc.y, desc.z);
    lum_a.map.resize(desc.x, desc.y, desc.z);
    lum_b.map.resize(desc.x, desc.y, desc.z);
    diss_a.map.resize(desc.x, desc.y, desc.z);
    diss_b.map.resize(desc.x, desc.y, desc.z);
    blue_a.map.resize(desc.x, desc.y, desc.z);
    blue_b.map.resize(desc.x, desc.y, desc.z);
    yllw_a.map.resize(desc.x, desc.y, desc.z);
    yllw_b.map.resize(desc.x, desc.y, desc.z);
    diss_c.map.resize(desc.x, desc.y, desc.z);
    diss_d.map.resize(desc.x, desc.y, desc.z);
    uns db = 0;
    for (uns i = 0; i < desc.z; i++)
        for (uns j = 0; j < desc.y; j++)
            for (uns k = 0; k < desc.x; k++) {
                rad_a.map.getref(k, j, i) = 0;
                rad_b.map.getref(k, j, i) = 0;
                lum_a.map.getref(k, j, i) = 0;
                lum_b.map.getref(k, j, i) = 0;
                diss_a.map.getref(k, j, i) = 0;
                diss_b.map.getref(k, j, i) = 0;
                blue_a.map.getref(k, j, i) = 0;
                blue_b.map.getref(k, j, i) = 0;
                yllw_a.map.getref(k, j, i) = 0;
                yllw_b.map.getref(k, j, i) = 0;
                diss_c.map.getref(k, j, i) = 0;
                diss_d.map.getref(k, j, i) = 0;
            }
    uns ertek_index = 0;
    for (uns i = 0; i < desc.cellnum; i++) {
        if (desc.cellak[i].is_el || desc.cellak[i].is_ph) {
            rad_a.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
        }
        if (desc.cellak[i].is_th) {
            blue_a.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
        }
        if (desc.cellak[i].is_el || desc.cellak[i].is_ph) {
            rad_b.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
        }
        if (desc.cellak[i].is_th) {
            blue_b.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
        }
        if (desc.cellak[i].is_el || desc.cellak[i].is_ph) {
            lum_a.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
        }
        if (desc.cellak[i].is_th) {
            yllw_a.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
        }
        if (desc.cellak[i].is_el || desc.cellak[i].is_ph) {
            lum_b.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
        }
        if (desc.cellak[i].is_th) {
            yllw_b.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
        }
        if (desc.cellak[i].is_el || desc.cellak[i].is_ph) {
            diss_a.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
        }
        if (desc.cellak[i].is_th) {
            diss_c.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
        }
        if (desc.cellak[i].is_el || desc.cellak[i].is_ph) {
            diss_b.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
        }
        if (desc.cellak[i].is_th) {
            diss_d.map.getref(desc.cellak[i].x, desc.cellak[i].y, desc.cellak[i].z) = res.dc_ertekek[ertek_index++];
        }
        //printf("%u\n", ertek_index);
        if(res.is_face_data)
            ertek_index += desc.cellak[i].f; // átugorja a face adatokat => nincs jól kezelve!!
    }
}
