//***********************************************************************
// Vector SUNRED PLList and tomb class header
// Creation date:    2004. 04. 23.
// Creator:            Pohl László
//***********************************************************************


//***********************************************************************
#ifndef VSUN2_LISTAESTOMB_HEADER
#define VSUN2_LISTAESTOMB_HEADER
//***********************************************************************


//****************************************************************
#include "PLString.h"
#include "hiba.h"
//****************************************************************

//****************************************************************
class hizotomb{
// 3D tömb double értékek tárolására
// a set megnöveli a tömb méretét, ha szükséges
// a get 0-t ad vissza, ha a koordináták kilógnak a tömbbõl
//****************************************************************
    double *m;
    unsigned x,y,z,n;
public:
    //****************************************************************
    hizotomb():x(0),y(0),z(0),n(0),m(NULL){}
    //****************************************************************
    ~hizotomb(){free();}
    //****************************************************************
    void free(){delete [] m; m = NULL; x = y = z = n = 0;}
    //****************************************************************
    hizotomb(unsigned x,unsigned y,unsigned z):x(x),y(y),z(z),n(x*y*z){
    //****************************************************************
        m=new double[n];
        for(unsigned i=0;i<n;i++)m[i]=0.0;
    }
    //****************************************************************
    hizotomb(const hizotomb & be){ m=NULL; *this=be; }
    //****************************************************************
    hizotomb & operator=(const hizotomb & be){
    //****************************************************************
        if(&be!=this){
            delete [] m;
            x=be.x;
            y=be.y;
            z=be.z;
            n=be.n;
            m=new double[n];
            for(unsigned i=0;i<n;i++)m[i]=be.m[i];
        }
        return *this;
    }
    //****************************************************************
    void expand(unsigned x,unsigned y=0,unsigned z=0){
    //****************************************************************
        if( x > this->x || y > this->y || z > this->z ){
            resize( x > this->x ? x : this->x, y > this->y ? y : this->y, z > this->z ? z : this->z );
        }
    }
    //****************************************************************
    void resize(unsigned x,unsigned y=0,unsigned z=0){
    //****************************************************************
        if( x != this->x || y != this->y || z != this->z ){
            unsigned n2 = x * y * z;
            double *m2 = new double[n2];
            for(unsigned i=0; i<n2; i++)m2[i]=0.0;
            const unsigned zz = z < this->z ? z : this->z ;
            const unsigned yy = y < this->y ? y : this->y ;
            const unsigned xx = x < this->x ? x : this->x ;
            for(unsigned k=0; k<zz; k++)
                for(unsigned j=0; j<yy; j++)
                    for(unsigned i=0; i<xx; i++)
                        m2[ i + j*x + k*x*y ] = m[ i + j*this->x + k*this->x*this->y ];
            delete [] m;
            m = m2;
            this->x = x;
            this->y = y;
            this->z = z;
            n = n2;
        }
    }
    //****************************************************************
    void set(double value,unsigned x,unsigned y=0,unsigned z=0){
    //****************************************************************
        expand(x+1,y+1,z+1);
        m[ x + y*this->x + z*this->x*this->y ] = value;
    }
    //****************************************************************
    void inc(double value,unsigned x,unsigned y=0,unsigned z=0){
    //****************************************************************
        expand(x+1,y+1,z+1);
        m[ x + y*this->x + z*this->x*this->y ] += value;
    }
    //****************************************************************
    double get(unsigned x,unsigned y=0,unsigned z=0)const{
    //****************************************************************
        if( x >= this->x || y >= this->y || z >= this->z ) return 0.0;
        return m[ x + y*this->x + z*this->x*this->y ];
    }
    //****************************************************************
    double sum(){ // a tömb elemeinek összegét adja.
    //****************************************************************
        double s=0.0;
        for(unsigned i=0; i<n; i++) s += m[ i ];
        return s;
    }
    //****************************************************************
    unsigned size()const{return n;}
    //****************************************************************
    void neg(){for(unsigned i=0; i<n; i++) if(m[i]!=0.0)m[i]=1.0/m[i];}
    //****************************************************************
    double sulyozott_osszeg(const hizotomb & suly){
    //****************************************************************
        double s=0.0;
        for(unsigned k=0; k<z; k++)
            for(unsigned j=0; j<y; j++)
                for(unsigned i=0; i<x; i++){
                    double v1 = get(i,j,k);
                    double v2 = suly.get(i,j,k);
                    s += v1 * v2;
                    if( v1 != 0.0 || v2 != 0.0 )
                        printf("%g * %g (%g)\n",v1,v2,s);
                }
        return s;
    }
    //****************************************************************
    void max_suruseg(const hizotomb & terulet,double & max_p,double & max_area){
    // ebben van a fényteljesítmény, paraméterként kapja a területeket
    // visszaadja azt a teljesítmény-terület párt, amelyik hányadosa maximális
    //****************************************************************
        double max = 0.0;
        for(unsigned k=0; k<z; k++)
            for(unsigned j=0; j<y; j++)
                for(unsigned i=0; i<x; i++){
                    double v1 = get(i,j,k);
                    double v2 = terulet.get(i,j,k);
                    if( v2 > 0.0 && v1/v2 > max ){
                        max = v1 / v2;
                        max_p = v1;
                        max_area = v2;
                    }
                }
    }
    //****************************************************************
    double sum_cut(double cut){ // A cut-nál nagyobb elemek összegét számolja
    //****************************************************************
        double s = 0.0;
        for(unsigned i = 0; i < n; i++)
            if( m[ i ] > cut ) s += m[ i ];
        return s;
    }
    //****************************************************************
    double sum_cut(const hizotomb & p,double cut){
    // teljesítménysûrûség alapján vág
    //****************************************************************
        double s=0.0;
        for(unsigned k=0; k<z; k++)
            for(unsigned j=0; j<y; j++)
                for(unsigned i=0; i<x; i++){
                    double v1 = p.get(i,j,k);
                    double v2 = get(i,j,k); // terület
                    if( v2 > 0.0 && v1/v2 > cut )
                        s += v2;
                }
        return s;
    }
    //****************************************************************
};


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
    bool is_exists(unsigned X, unsigned Y, unsigned Z) {
    //****************************************************************
        return X < x && Y < y && Z < z;
    }

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

    //****************************************************************
    void save(FILE * fp){
    //****************************************************************
        unsigned teszt='3';
        if(fwrite(&teszt,sizeof(unsigned),1,fp)!=1)throw hiba("tomb3d::save","teszt write failed");
        if(fwrite(&x,sizeof(unsigned),1,fp)!=1)throw hiba("tomb3d::save","x write failed");
        if(fwrite(&y,sizeof(unsigned),1,fp)!=1)throw hiba("tomb3d::save","y write failed");
        if(fwrite(&z,sizeof(unsigned),1,fp)!=1)throw hiba("tomb3d::save","z write failed");
        if(fwrite(&n,sizeof(unsigned),1,fp)!=1)throw hiba("tomb3d::save","n write failed");
        if(fwrite(m,sizeof(C),n,fp)!=n)throw hiba("tomb3d::save","m write failed");
    }

    //****************************************************************
    void load(FILE * fp){
    //****************************************************************
        unsigned teszt;
        if(fread(&teszt,sizeof(unsigned),1,fp)!=1)throw hiba("tomb3d::load","teszt read failed");
        if(teszt!='3')throw hiba("tomb3d::load","not a tomb3d to read");
        if(fread(&x,sizeof(unsigned),1,fp)!=1)throw hiba("tomb3d::load","x read failed");
        if(fread(&y,sizeof(unsigned),1,fp)!=1)throw hiba("tomb3d::load","y read failed");
        if(fread(&z,sizeof(unsigned),1,fp)!=1)throw hiba("tomb3d::load","z read failed");
        if(fread(&n,sizeof(unsigned),1,fp)!=1)throw hiba("tomb3d::load","n read failed");
        delete [] m;
        m=new C[n];
        if(fread(m,sizeof(C),n,fp)!=n)throw hiba("tomb3d::load","m read failed");
    }

    //****************************************************************
    void store(FILE * fp){
    //****************************************************************
        save(fp);
        delete [] m;
        x=y=z=n=0;
        m=NULL;
    }
};


//****************************************************************
template<class C> class tomb2d{
//****************************************************************
    C * m;
    unsigned x,y,n;
public:
    //****************************************************************
    tomb2d(){m=0;x=y=n=0;}
    void free(){delete [] m;m=0;x=y=n=0;}
    ~tomb2d(){free();}

    unsigned x_size()const{return x;}
    unsigned y_size()const{return y;}
    unsigned size()const { return n; }

    //****************************************************************


    //****************************************************************
    C & operator[](unsigned index){
    //****************************************************************
//#ifdef vsundebugmode
        if(index>=n)throw hiba("tomb2d::operator[]","index>=n (%u>=%u)",index,n);
//#endif
        return m[index];
    }


    //****************************************************************
    C & getref(unsigned X,unsigned Y){
    //****************************************************************
#ifdef vsundebugmode
        if(X>=x||Y>=y)
            throw hiba("tomb2d::getref","X>=x||Y>=y (%u>=%u||%u>=%u)",X,x,Y,y);
#endif
        return m[Y*x+X];
    }


    //****************************************************************
    const C & getconstref(unsigned X,unsigned Y)const{
    //****************************************************************
#ifdef vsundebugmode
        if(X>=x||Y>=y)throw hiba("tomb2d::getconstref","X>=x||Y>=y (%u>=%u||%u>=%u)",X,x,Y,y);
#endif
        return m[Y*x+X];
    }


    //****************************************************************
    C get(unsigned X,unsigned Y){
    //****************************************************************
#ifdef vsundebugmode
        if(X>=x||Y>=y)throw hiba("tomb2d::get","X>=x||Y>=y (%u>=%u||%u>=%u)",X,x,Y,y);
#endif
        return m[Y*x+X];
    }


    //****************************************************************
    tomb2d(const tomb2d & t){
    //****************************************************************
        m=new C[n=t.n];
        if(m==0)throw hiba("tomb2d::tomb2d(tomb2d&)","alloc failed");
        x=t.x; y=t.y;
        for(unsigned i=0;i<n;i++)m[i]=t.m[i];
    }


    //****************************************************************
    void resize(unsigned new_x,unsigned new_y=1){
    //****************************************************************
        if(m)delete [] m;
        m=new C[n=(x=new_x)*(y=new_y)];
        if(m==0)throw hiba("tomb2d::resize","alloc failed");
    }

    //****************************************************************
    tomb2d & operator=(const tomb2d & a){
    //****************************************************************
        if(this!=&a){
            if(m)delete [] m;
            m=new C[n=(x=a.x)*(y=a.y)];
            if(m==0)throw hiba("tomb2d::operator=","alloc failed");
            for(unsigned i=0;i<n;i++)m[i]=a.m[i];
        }
        return *this;
    }

    //****************************************************************
    void save(FILE * fp){
    //****************************************************************
        unsigned teszt='2';
        if(fwrite(&teszt,sizeof(unsigned),1,fp)!=1)throw hiba("tomb2d::save","teszt write failed");
        if(fwrite(&x,sizeof(unsigned),1,fp)!=1)throw hiba("tomb2d::save","x write failed");
        if(fwrite(&y,sizeof(unsigned),1,fp)!=1)throw hiba("tomb2d::save","y write failed");
        if(fwrite(&n,sizeof(unsigned),1,fp)!=1)throw hiba("tomb2d::save","n write failed");
        if(fwrite(m,sizeof(C),n,fp)!=n)throw hiba("tomb2d::save","m write failed");
    }

    //****************************************************************
    void load(FILE * fp){
    //****************************************************************
        unsigned teszt;
        if(fread(&teszt,sizeof(unsigned),1,fp)!=1)throw hiba("tomb2d::load","teszt read failed");
        if(teszt!='2')throw hiba("tomb2d::load","not a tomb2d to read");
        if(fread(&x,sizeof(unsigned),1,fp)!=1)throw hiba("tomb2d::load","x read failed");
        if(fread(&y,sizeof(unsigned),1,fp)!=1)throw hiba("tomb2d::load","y read failed");
        if(fread(&n,sizeof(unsigned),1,fp)!=1)throw hiba("tomb2d::load","n read failed");
        delete [] m;
        m=new C[n];
        if(fread(m,sizeof(C),n,fp)!=n)throw hiba("tomb2d::load","m read failed");
    }

    //****************************************************************
    void store(FILE * fp){
    //****************************************************************
        save(fp);
        delete [] m;
        m=NULL;
    }
};


//****************************************************************
template<class C> class tomb
//****************************************************************
{
    C * m;
    unsigned n,findindex,allocated;
public:
    //****************************************************************
    tomb() :n(0), m(NULL), findindex(0), allocated{ 0 } {}
    ~tomb(){delete [] m;m=NULL;}

    unsigned size()const{return n;}
    unsigned getfindindex()const{return findindex;}
    void clear() { delete[] m; m = NULL; n = findindex = 0; allocated = 0; }
    void swap(tomb<C> &masik) { 
        C * ptemp = m; m = masik.m; masik.m = ptemp;
        unsigned utemp = n; n = masik.n; masik.n = utemp;
        utemp = findindex; findindex = masik.findindex; masik.findindex = utemp;
        utemp = allocated; allocated = masik.allocated; masik.allocated = utemp;
    }
    void zero() {
        for (unsigned i = 0; i < n; i++)
            m[i] = C();
    }
    unsigned get_index(const C * p)const { return (unsigned)(p - m); }
    //****************************************************************


    //****************************************************************
    tomb(const tomb & t){
    //****************************************************************
        m=new C[n=t.n];
        allocated = n;
        findindex = t.findindex;
        if(m==0)throw hiba("tomb::tomb(tomb&)","alloc failed");
        for(unsigned i=0;i<n;i++)m[i]=t.m[i];
    }


    //****************************************************************
    C * add(C elem,unsigned db=1){
    // ha többet tesz be, akkor az ezek közül utolsó címét adja vissza
    //****************************************************************
        C *m2=new C[n+db];
        allocated = n + db;
        if(m2==0)throw hiba("tomb::add","alloc failed");
        if(m){for(unsigned i=0;i<n;i++)m2[i]=m[i];delete [] m;}
        while(db--)m2[n++]=elem;
        m=m2;
        return &m2[n-1];
    }


    //****************************************************************
    void add(const tomb<C> &t){
    //****************************************************************
        if (t.n == 0)
            return;
        C *m2=new C[n+t.n];
        allocated = n + t.n;
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
        allocated = newsize;
        if(m2==0)throw hiba("tomb::resize","alloc failed");
        unsigned masol=(newsize<n)?newsize:n;
        if(m){for(unsigned i=0;i<masol;i++)m2[i]=m[i];delete [] m;}
        n=newsize;
        m=m2;
    }

    //****************************************************************
    void reallocate(unsigned min_size) {
    //****************************************************************
        if (allocated >= min_size) {
            n = min_size;
            return;
        }
        allocated = 3 * allocated / 2 > min_size ? 3 * allocated / 2 : min_size;
        C * m2 = new C[allocated];
        for (unsigned i = 0; i < n; i++)
            m2[i] = m[i];
        n = min_size;
        delete[] m;
        m = m2;
    }


    //****************************************************************
    C * add_with_reallocate(C elem){
    //****************************************************************
        reallocate(n + 1);
        m[n-1] = elem;
        return &m[n-1];
    }


    //****************************************************************
    C * add_with_reallocate(){
    //****************************************************************
        reallocate(n + 1);
        return &m[n-1];
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
        if (n == 0 && a.n == 0)
            return *this;
        if(this!=&a){
            if (n != a.n) {
                delete[] m;
                m = new C[n = a.n];
            }
            if(m==0)throw hiba("tomb::operator=","alloc failed");
            for(unsigned i=0;i<n;i++)m[i]=a.m[i];
        }
        return *this;
    }


    //****************************************************************
    bool findUnsorted(C elem){
    // ha létezik az == operátor a C osztályhoz, akkor használható csak
    // nem rendezett tömb
    //****************************************************************
        for(findindex=0;findindex<n;findindex++)if(m[findindex]==elem)return true;
        return false;
    }


    //****************************************************************
    bool findSorted(C elem){
    // ha létezik az == operátor a C osztályhoz, akkor használható csak
    // rendezett tömb
    //****************************************************************
        for(findindex=0;findindex<n&&m[findindex]<=elem;findindex++)if(m[findindex]==elem)return true;
        return false;
    }


    //****************************************************************
    void insertIfNotExists(C elem){
    // ha létezik már az elem, akkor nem szúrja be (a tömb rendezett)
    //****************************************************************
        if(findSorted(elem))return;

        C *m2=new C[n+1];
        allocated = n + 1;
        unsigned i=0;

        if(m2==0)throw hiba("tomb::insertIfNotExists","alloc failed");
        if(m){for(;i<findindex;i++)m2[i]=m[i];}
        m2[i]=elem;
        if(m){for(;i<n;i++)m2[i+1]=m[i];delete [] m;}
        n++;m=m2;
    }


    //****************************************************************
    void insertAlways(C elem){
    // beszúrja a megfelelõ helyre akkor is, ha már van ilyen elem
    //****************************************************************
        findSorted(elem);

        C *m2=new C[n+1];
        allocated = n + 1;
        unsigned i=0;

        if(m2==0)throw hiba("tomb::insertAlways","alloc failed");
        if(m){for(;i<findindex;i++)m2[i]=m[i];}
        m2[i]=elem;
        if(m){for(;i<n;i++)m2[i+1]=m[i];delete [] m;}
        n++;m=m2;
    }


    //****************************************************************
    void insert(C elem,unsigned pozicio){
    // beszúrja a megfelelõ helyre
    //****************************************************************
        if(pozicio>n)throw hiba("tomb::insert()","inserting out of array");

        C *m2=new C[n+1];
        allocated = n + 1;
        unsigned i=0;

        if(m2==0)throw hiba("tomb::insert","alloc failed");
        if(m){for(;i<pozicio;i++)m2[i]=m[i];}
        m2[i]=elem;
        if(m){for(;i<n;i++)m2[i+1]=m[i];delete [] m;}
        n++;m=m2;
    }


    //****************************************************************
    void removeElem(C elem){
    //****************************************************************
        unsigned i=0;
        for(;i<n;i++)if(m[i]==elem)break;
        if(m[i]!=elem)throw hiba("tomb::removeElem","Elem missing");
        removeN(i);
    }


    //****************************************************************
    void removeN(unsigned N){
    //****************************************************************
        if(N>=n)throw hiba("tomb::removeN","trying remove elemet out of array");
        if(n==1){delete [] m;m=0;n=findindex=0;return;}
        C *m2=new C[n-1]; if(m2==0)throw hiba("tomb::insert","alloc failed");
        allocated = n - 1;
        unsigned i; for(i=0;i<N;i++)m2[i]=m[i];
        for(++i;i<n;i++)m2[i-1]=m[i];
        delete m; n--; m=m2;
    }

    //****************************************************************
    void save(FILE * fp){
    //****************************************************************
        unsigned teszt='1';
        if(fwrite(&teszt,sizeof(unsigned),1,fp)!=1)throw hiba("tomb::save","teszt write failed");
        if(fwrite(&findindex,sizeof(unsigned),1,fp)!=1)throw hiba("tomb::save","findindex write failed");
        if(fwrite(&n,sizeof(unsigned),1,fp)!=1)throw hiba("tomb::save","n write failed");
        if(fwrite(m,sizeof(C),n,fp)!=n)throw hiba("tomb::save","m write failed");
    }

    //****************************************************************
    void load(FILE * fp){
    //****************************************************************
        unsigned teszt;
        if(fread(&teszt,sizeof(unsigned),1,fp)!=1)throw hiba("tomb::load","teszt read failed");
        if(teszt!='1')throw hiba("tomb::load","not a tomb to read");
        if(fread(&findindex,sizeof(unsigned),1,fp)!=1)throw hiba("tomb::load","findindex read failed");
        if(fread(&n,sizeof(unsigned),1,fp)!=1)throw hiba("tomb::load","n read failed");
        delete [] m;
        m=new C[n];
        allocated = n;
        if(fread(m,sizeof(C),n,fp)!=n)throw hiba("tomb::load","m read failed");
    }

    //****************************************************************
    void store(FILE * fp){
    //****************************************************************
        save(fp);
        delete [] m;
        m=NULL;
    }
};

//****************************************************************
template<class C> class asszoc{
//****************************************************************
    struct ex{
        C elem;
        PLString s;
        ex & operator=(ex & src){if(&src!=this){elem=src.elem;s=src.s;}return *this;}
        ex(ex & src){elem=src.elem;s=src.s;}
        ex(){}
    };
    ex * m;
    unsigned n,findindex;
public:
    //****************************************************************
    asszoc(){m=0;n=findindex=0;}
    ~asszoc(){if(m)delete [] m;}

    unsigned size()const{return n;}
    unsigned getfindindex()const{return findindex;}
    void clear(){delete [] m; m=0;n=0;}
    //****************************************************************


    //****************************************************************
    asszoc(const asszoc & t){
    //****************************************************************
        m=new ex[n=t.n];
        if(m==0)throw hiba("asszoc::asszoc(asszoc&)","alloc failed");
        for(unsigned i=0;i<n;i++)m[i]=t.m[i];
    }


    //****************************************************************
    C * add(const char * s,C elem){return add(PLString(s),elem);}
    C * add(PLString s,C elem){
    //****************************************************************
        ex *m2=new ex[n+1];
        if(m2==0)throw hiba("asszoc::add","alloc failed");
        if(m){for(unsigned i=0;i<n;i++)m2[i]=m[i];delete [] m;}
        m2[n].elem=elem;
        m2[n].s=s;
        n++;
        m=m2;
        return &m2[n-1].elem;
    }


    //****************************************************************
    void resize(unsigned newsize){
    //****************************************************************
        ex *m2=new ex[newsize];
        if(m2==0)throw hiba("asszoc::resize","alloc failed");
        unsigned masol=(newsize<n)?newsize:n;
        if(m){for(unsigned i=0;i<masol;i++)m2[i]=m[i];delete [] m;}
        n=newsize;
        m=m2;
    }


    //****************************************************************
    C & operator[](unsigned index){
    //****************************************************************
//#ifdef vsundebugmode
        if(index>=n)throw hiba("asszoc::operator[]","index>=n (%u>=%u)",index,n);
//#endif
        return m[index].elem;
    }


    //****************************************************************
    bool index(const char * s,C & ret){return index(PLString(s),ret);}
    bool index(PLString s,C & ret){
    //****************************************************************
        for(unsigned i=0;i<n;i++)if(s.UpCase()==m[i].s.UpCase()){ret=m[i].elem;return true;}
        return false;
    }


    //****************************************************************
    asszoc & operator=(const asszoc & a){
    //****************************************************************
        if(this!=&a){
            if(m)delete [] m;
            m=new ex[n=a.n];
            if(m==0)throw hiba("asszoc::operator=","alloc failed");
            for(unsigned i=0;i<n;i++)m[i]=a.m[i];
        }
        return *this;
    }


    //****************************************************************
    bool findUnsorted(C elem){
    // ha létezik az == operátor a C osztályhoz, akkor használható csak
    // nem rendezett tömb
    //****************************************************************
        for(findindex=0;findindex<n;findindex++)if(m[findindex].elem==elem)return true;
        return false;
    }


    //****************************************************************
    bool findUnsorted(C elem,PLString & key){
    // ha létezik az == operátor a C osztályhoz, akkor használható csak
    // nem rendezett tömb
    // key-ben adja vissza a kulcsot
    //****************************************************************
        for(findindex=0;findindex<n;findindex++)if(m[findindex].elem==elem){key=m[findindex].s; return true;}
        return false;
    }


    //****************************************************************
    bool findSorted(C elem){
    // ha létezik az == operátor a C osztályhoz, akkor használható csak
    // rendezett tömb
    //****************************************************************
        for(findindex=0;findindex<n&&m[findindex].elem<=elem;findindex++)if(m[findindex].elem==elem)return true;
        return false;
    }


    //****************************************************************
    void removeN(unsigned N){
    //****************************************************************
        if(N>=n)throw hiba("asszoc::removeN","trying remove elemet out of array");
        if(n==1){delete [] m;m=0;n=findindex=0;return;}
        ex *m2=new ex[n-1]; if(m2==0)throw hiba("asszoc::insert","alloc failed");
        unsigned i; for(i=0;i<N;i++)m2[i]=m[i];
        for(++i;i<n;i++)m2[i-1]=m[i];
        delete m; n--; m=m2;
    }

    //****************************************************************
    PLString & key(unsigned index){
    //****************************************************************
        if(index>=n)throw hiba("asszoc::key","index error");
        return m[index].s;
    }
};

#define PL_MAX_VEREM_MERET    256

//****************************************************************
template<class C> class PLList
//****************************************************************
{
    //************************************************************
    struct DataUnit
    //************************************************************
    {
        C            D;
        DataUnit    *Prev,*Next;
    };

    DataUnit    *First,*Last,*Act;
    DataUnit    *Verem[PL_MAX_VEREM_MERET];
    int            VeremElemSzam;
    C            Dummy;
    DataUnit    DDummy;
public:


    //************************************************************
    PLList&    operator=(PLList &L)
    //************************************************************
    {
        if(this!=&L)
        {
            Free();
            L.Act=L.First;
            while(L.Act!=0)
            {
                Dummy=L.Act->D;
                Add(Dummy);
                L.Act=L.Act->Next;
            }
            L.Act=L.Last;
        }
        return *this;
    }


    //************************************************************
    PLList(PLList& L)
    //************************************************************
    {
        First=Last=Act=0;
        L.Act=L.First;
        while(L.Act!=0)
        {
            Dummy=L.Act->D;
            Add(Dummy);
            L.Act=L.Act->Next;
        }
        L.Act=L.Last;
    }


    //************************************************************
    PLList(){First=Last=Act=0;VeremElemSzam=0;}
    ~PLList(){Free();}


    C        GetFirst(){Act=First;if(Act!=0)return Act->D;else return Dummy;}
    C        GetLast(){Act=Last;if(Act!=0)return Act->D;else return Dummy;}
    C        GetAct(){if(Act==0)return Dummy;return Act->D;}
    C&        GetFirstRef(){Act=First;if(Act!=0)return Act->D;else return Dummy;}
    C&        GetLastRef(){Act=Last;if(Act!=0)return Act->D;else return Dummy;}
    C&        GetActRef(){if(Act==0)return Dummy;return Act->D;}
    C*        GetFirstPoi(){return &GetFirstRef();}
    C*        GetLastPoi(){return &GetLastRef();}
    C*        GetActPoi(){return &GetActRef();}
    C*        GetNextPoi(){return &GetNextRef();}
    bool    Begin(){return Act==First;}
    void    SetBegin(){Act=First;}
    bool    End(){
        if(Last)return Act==Last;
        else return true;
    }
    void    SetEnd(){Act=Last;}
    void    SetNext(){if((Act!=0)&&(Act->Next!=0))Act=Act->Next;}
    void    StepBack(){if(!Begin())Act=Act->Prev;}
    void    Free(){while(First!=0){Act=First->Next;delete First;First=Act;}First=Last=Act=0;}
    bool    Empty(){return (Last==0)?true:false;}
    int        Count(){int i=0;PushAct();for(Act=First;Act!=0;i++,Act=Act->Next);PopAct();return i;}
    void    ResetToDDummy(){DDummy.Prev=0;DDummy.Next=First;Act=&DDummy;}
    void    PushAct(){if(VeremElemSzam<PL_MAX_VEREM_MERET)Verem[VeremElemSzam++]=Act;else throw hiba("PLList::PopAct()","PLList Verem overflowed.");}
    void    PopAct(){if(VeremElemSzam>0)Act=Verem[--VeremElemSzam];else throw hiba("PLList::PopAct()","PLList Verem underflowed.");}

    friend class TMLib;
    //************************************************************


    //************************************************************
    C    *Add(C &Data)
    //************************************************************
    {
        DataUnit *DU=new DataUnit;
        DU->D=Data;
        DU->Next=0;
        if(Last==0){
            DU->Prev=0;
            First=DU;
        }
        else{
            DU->Prev=Last;
            Last->Next=DU;
        }
        Last=Act=DU;
        return &(Act->D);
    }


    //************************************************************
    void    Insert(C &Data)//a kijelölt mögé szúrja be
    //************************************************************
    {
        DataUnit *DU=new DataUnit;
        DU->D=Data;
        if(End())Add(Data);//Ha nincs elem a listában, akkor is mûködik,mert Last==Act==0
        else
        {
            DU->Prev=Act;
            DU->Next=Act->Next;
            Act->Next=DU;
            DU->Next->Prev=DU;
            Act=DU;
        }
    }


    //************************************************************
    void    InsertFirst(C &Data)
    //************************************************************
    {
        DataUnit *DU=new DataUnit;
        DU->D=Data;
        if(End())Add(Data);//Ha nincs elem a listában, akkor is mûködik,mert Last==Act==0
        else
        {
            DU->Prev=0;
            DU->Next=First;
            First->Prev=DU;
            First=Act=DU;
        }
    }


    //************************************************************
    void    DelAct()
    //************************************************************
    {
        if(Act==0)return;
        if(Act->Prev!=0)Act->Prev->Next=Act->Next;
        if(Act->Next!=0)Act->Next->Prev=Act->Prev;
        if(Act==First)First=Act->Next;
        if(Act==Last)Last=Act->Prev;
        DataUnit *DU=Act;
        if(Act->Prev!=0)Act=Act->Prev;
        else Act=Act->Next;
        delete DU;
    }


    //************************************************************
    void    DelLast()
    //************************************************************
    {
        Act=Last;
        DelAct();
    }


    //************************************************************
    C        GetNext()
    //************************************************************
    {
        if(Act==0)return Dummy;
        if(Act->Next==0)return Act->D;
        Act=Act->Next;
        return Act->D;
    }


    //************************************************************
    C&        GetNextRef()
    //************************************************************
    {
        if(Act==0)return Dummy;
        if(Act->Next==0)return Act->D;
        Act=Act->Next;
        return Act->D;
    }

};


/****************************************************************/
inline bool StringKereso(PLList<PLString> &Lista,PLString &Elem)
//ha a Listában szerepel az Elem, akkor true
/****************************************************************/
{
    Lista.PushAct();
    Lista.ResetToDDummy();
    while(!Lista.End())
    {
        if(Lista.GetNextRef().UpCase()==Elem.UpCase())
        {
            Lista.PopAct();
            return true;
        }
    }
    Lista.PopAct();
    return false;
}


#endif
