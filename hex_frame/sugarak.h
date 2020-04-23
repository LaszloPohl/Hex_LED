//***********************************************************************
// sugarak header
// Creation date:  2020. 02. 06.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef HEX_FRAME_SUGARAK_HEADER
#define	HEX_FRAME_SUGARAK_HEADER
//***********************************************************************


//***********************************************************************
#include "tipusok.h"
//***********************************************************************


//***********************************************************************
// Az alap a sugar_cella struktúrák listája.
// A lista felépítését a junction cellához tartozó sugar_cella foglalásá-
// val kezdjük.
// Ennek van egy indulo_sugarak listája, amit feltöltünk a junction face
// és a fénypor cellák peremfaceinek középpontját összekötõ kék sugarak-
// kal.
// Töröljük a gyenge sugarakat.
// A maradék kék sugárhoz létrehozzuk a sugar_cella listákat, melyek kez-
// dõeleme a cellasugarak listában van, utolsó elemük a junction cella.
// Létrehozzuk a sárga sugarakat a kék sugar_cella-kban.
// Töröljük a gyenge sugarakat.
// Létrehozzuk a sárga sugar_cella listákat.
//***********************************************************************


//***********************************************************************
struct t_modell_cella;
//***********************************************************************


//***********************************************************************
struct kulso_oldal {
//***********************************************************************
    Vec3d p1, p2, p3, p4, c;
    const t_modell_cella * p_modell_cella;
    uns cella_x, cella_y, cella_z;
    bool is_oldalso;
    //***********************************************************************
    kulso_oldal() :p_modell_cella{ nullptr }, cella_x{ 0 }, cella_y{ 0 }, cella_z{ 0 }, is_oldalso{ false } {}
    //***********************************************************************
    void set(dbl x1, dbl y1, dbl z1, dbl x2, dbl y2, dbl z2, const t_modell_cella & modell_cella, uns x, uns y, uns z, bool oldalso_e) {
    //***********************************************************************
        p_modell_cella = &modell_cella;
        is_oldalso = oldalso_e;
        cella_x = x;
        cella_y = y;
        cella_z = z;
        p1.x = x1;
        p1.y = y1;
        p1.z = z1;
        p3.x = x2;
        p3.y = y2;
        p3.z = z2;
        if (z1 == z2) {
            p2.x = x1;
            p2.y = y2;
            p2.z = z1;
            p4.x = x2;
            p4.y = y1;
            p4.z = z1;
        }
        else if (y1 == y2) {
            p2.x = x1;
            p2.y = y1;
            p2.z = z2;
            p4.x = x2;
            p4.y = y1;
            p4.z = z1;
        }
        else {
            p2.x = x1;
            p2.y = y1;
            p2.z = z2;
            p4.x = x1;
            p4.y = y2;
            p4.z = z1;
        }
        c = 0.5*(p1 + p3);
    }
};


//***********************************************************************
struct egyenes {
//***********************************************************************
    Vec3d p1, p2;
    //***********************************************************************
    void set(const Vec3d & P1, const Vec3d & P2) { p1 = P1; p2 = P2; }
    //***********************************************************************
    Vec3d get_x_metszet(dbl x)const {
    //***********************************************************************
        return Vec3d{ x, (x - p1.x)*(p2.y - p1.y) / (p2.x - p1.x) + p1.y, (x - p1.x)*(p2.z - p1.z) / (p2.x - p1.x) + p1.z }; 
    }
    //***********************************************************************
    Vec3d get_y_metszet(dbl y)const {
    //***********************************************************************
        return Vec3d{ (y - p1.y)*(p2.x - p1.x) / (p2.y - p1.y) + p1.x, y, (y - p1.y)*(p2.z - p1.z) / (p2.y - p1.y) + p1.z };
    }
    //***********************************************************************
    Vec3d get_z_metszet(dbl z)const {
    //***********************************************************************
        return Vec3d{ (z - p1.z)*(p2.x - p1.x) / (p2.z - p1.z) + p1.x, (z - p1.z)*(p2.y - p1.y) / (p2.z - p1.z) + p1.y, z };
    }
    //***********************************************************************
    bool get_x_metszet(dbl x, dbl y0, dbl y1, dbl z0, dbl z1, Vec3d & metszespont) const{
    // y0<=y1 és z0<=z1 kell legyen
    //***********************************************************************
        if (p2.x == p1.x)
            return false;
        Vec3d mp = get_x_metszet(x);
        if (mp.y<y0 || mp.y>y1 || mp.z<z0 || mp.z>z1)
            return false;
        metszespont = mp;
        return true;
    }
    //***********************************************************************
    bool get_y_metszet(dbl y, dbl x0, dbl x1, dbl z0, dbl z1, Vec3d & metszespont) const{
    // x0<=x1 és z0<=z1 kell legyen
    //***********************************************************************
        if (p2.y == p1.y)
            return false;
        Vec3d mp = get_y_metszet(y);
        if (mp.x<x0 || mp.x>x1 || mp.z<z0 || mp.z>z1)
            return false;
        metszespont = mp;
        return true;
    }
    //***********************************************************************
    bool get_z_metszet(dbl z, dbl x0, dbl x1, dbl y0, dbl y1, Vec3d & metszespont) const{
    // x0<=x1 és y0<=y1 kell legyen
    //***********************************************************************
        if (p2.z == p1.z)
            return false;
        Vec3d mp = get_z_metszet(z);
        if (mp.x<x0 || mp.x>x1 || mp.y<y0 || mp.y>y1)
            return false;
        metszespont = mp;
        return true;
    }
    //***********************************************************************
    Vec3d get_p2_kozelpont()const {
    //***********************************************************************
        return (p2.z != p1.z) 
            ? get_z_metszet(p2.z + (p1.z - p2.z)*1e-9)
            : ((p2.y != p1.y) 
                ? get_y_metszet(p2.y + (p1.y - p2.y)*1e-9)
                : get_x_metszet(p2.x + (p1.x - p2.x)*1e-9)); 
    }
    //***********************************************************************
    dbl get_p1_p2_distance()const {
    //***********************************************************************
        Vec3d d = p2 - p1;
        return d.hossz();
    }
    //***********************************************************************
    dbl get_vertical_angle()const {
    //***********************************************************************
        Vec3d d = p2 - p1;
        return d.get_vertical_angle();
    }
};


//***********************************************************************
struct sugar_adat {
//***********************************************************************
    uns src_cella_index;
    unsigned short face_index;
    unsigned short x0, y0, z0; // a kiinduló junction/fénypor cella koordinátái
    unsigned short x1, y1, z1; // a végzõdõ fénypor cella koordinátái
    egyenes e;
    dbl P; // kéknél a felület és a látószög szorzata, sárgánál még a bemenõ sugár arányával is szorozva
    dbl A_kek; // kéknél a junction cella felülete
    dbl P_kek_sarga; // A sugár bemenete a sugárzó cella teljesítményének mekkora részét jelképezi
    bool is_oldalso;
    //***********************************************************************
    sugar_adat *prev, *next;
    //***********************************************************************
    sugar_adat() :prev{ nullptr }, next{ nullptr }, A_kek{ 0 }, is_oldalso{ false } {}
    //***********************************************************************
    void remove() {
    // Ne közvetlenül hívjuk, hanem a sugar_lista!
    //***********************************************************************
        prev->next = next;
        next->prev = prev;
    }

};


//***********************************************************************
struct sugar_cella {
//***********************************************************************
    uns cella_index;
    unsigned short x, y, z, face_index;
    float K_P, d;
    char dir1, dir2;
    bool is_vertical, is_oldalsohoz;
    sugar_cella() :K_P{ 0 }, d{ 0 }, cella_index{ 0 }, face_index{ 0 }, x{ 0 }, y{ 0 }, 
        z{ 0 }, dir1{ 'x' }, dir2{ 'x' }, is_vertical{ false }, is_oldalsohoz{ false } {}
};


//***********************************************************************
struct sugar_lista {
//***********************************************************************
    //***********************************************************************
    sugar_adat start, stop;
    sugar_cella kiindulo_cella;
    uns db;
    //***********************************************************************
    sugar_lista() : db{ 0 } { start.next = &stop; stop.prev = &start; }
    //***********************************************************************
    ~sugar_lista() {
    //***********************************************************************
        free();
    }
    //***********************************************************************
    sugar_lista & operator=(sugar_lista & masik) {
    //***********************************************************************
        start = masik.start;
        stop = masik.stop;
        start.next->prev = &start;
        stop.prev->next = &stop;
        kiindulo_cella = masik.kiindulo_cella;
        db = masik.db;
        masik.start.next = &masik.stop;
        masik.stop.prev = &masik.start;
        masik.db = 0;
        return *this;
    }
    //***********************************************************************
    void free() {
    //***********************************************************************
        sugar_adat * pakt = start.next;
        while (pakt->next != nullptr) {
            sugar_adat * temp = pakt->next;
            delete pakt;
            pakt = temp;
        }
        start.next = &stop; stop.prev = &start;
    }
    //***********************************************************************
    sugar_adat * add(const sugar_adat & uj) {
    // a lista végére teszi
    //***********************************************************************
        sugar_adat *p_uj = new sugar_adat{ uj };
        p_uj->next = &stop;
        p_uj->prev = stop.prev;
        stop.prev->next = p_uj;
        stop.prev = p_uj;
        db++;
        return p_uj;
    }
    //***********************************************************************
    void add(const sugar_adat & uj, uns db, bool is_fele, const Vec3d & dp1i, const Vec3d & dp1j, const Vec3d & dp2i, const Vec3d & dp2j) {
    // a lista végére teszi
    // db: irányonként ennyi sugár
    // dpxy: x: kezdõ vagy végpont, y: i vagy j irányban lépve
    //***********************************************************************
        sugar_adat akt = uj;
        cuns n = db*db;
        akt.A_kek /= n;
        akt.P /= n;
        akt.P_kek_sarga /= n;
        int idb = (db - 1) / 2;
        uns iii = 0;
        for (int i = -idb; i <= idb; i++) {
            for (int j = -idb; j <= idb; j++, iii++) 
                if (!is_fele || iii % 2 == 0) {
                    akt.e.p1 = (double)i*dp1i + (double)j*dp1j + uj.e.p1;
                    akt.e.p2 = (double)i*dp2i + (double)j*dp2j + uj.e.p2;
                    add(akt);
                }
        }
    }
    //***********************************************************************
    void remove(sugar_adat * p) {
    //***********************************************************************
        p->remove();
        db--;
    }
    //***********************************************************************
    void recount() {
    //***********************************************************************
        db = 0;
        for (auto it = start.next; it != &stop; db++, it = it->next)
            ;
    }
    //***********************************************************************
    void normalize_P_kek_sarga() {
    //***********************************************************************
        dbl sum = 0;
        for (auto it = start.next; it != &stop; it = it->next)
            sum += it->P_kek_sarga;
        dbl isum = 1 / sum;
        for (auto it = start.next; it != &stop; it = it->next) {
            it->P_kek_sarga *= isum;
        }
    }
    //***********************************************************************
    void recalculate_P_kek(dbl & sum_dP, dbl & sum_P) {
    // kék sugárnál P = A_kek * P_kek_sarga
    // paraméterek a súlyozott hossz kiszámításához
    //***********************************************************************
        dbl dP = 0, P = 0;
        for (auto it = start.next; it != &stop; it = it->next) {
            it->P = it->A_kek * it->P_kek_sarga;
            P += it->P;
            dP += it->P * it->e.get_p1_p2_distance();
        }
        sum_dP += dP;
        sum_P += P;
    }
};


//***********************************************************************
struct global_sugar_tomb {
//***********************************************************************
    tomb<sugar_adat*> t;
    //***********************************************************************
    void clear() { t.clear(); }
    //***********************************************************************
    void add(sugar_lista & sl) {
    //***********************************************************************
        uns start_index = t.size();
        uns n = start_index + sl.db;
        t.reallocate(n);
        sugar_adat * it = sl.start.next;
        for (uns i = start_index; i < n; i++, it=it->next) {
            t[i] = it;
        }
    }
    //***********************************************************************
    void add(tomb<sugar_adat*> mt, uns mstart) {
    //***********************************************************************
        uns start_index = t.size();
        uns n = start_index + mt.size() - mstart;
        t.reallocate(n);
        for (uns i = start_index; i < n; i++, mstart++) {
            t[i] = mt[mstart];
        }
    }
    //***********************************************************************
};


//***********************************************************************
struct cella_lista {
//***********************************************************************
    //***********************************************************************
    tomb<sugar_cella> cellatomb;
    //***********************************************************************
    cella_lista * next;
    //***********************************************************************
    cella_lista() :next{ nullptr } {}
    //***********************************************************************
    sugar_cella & add() { cellatomb.reallocate(cellatomb.size() + 1); return cellatomb.getLast(); }
    //***********************************************************************
};


//***********************************************************************
struct cellasugar_lista {
//***********************************************************************
    //***********************************************************************
    cella_lista * lista;
    //***********************************************************************
    cellasugar_lista() :lista{ nullptr } {}
    //***********************************************************************
    ~cellasugar_lista() { clear(); }
    //***********************************************************************
    cella_lista * add() { cella_lista * uj = new cella_lista; uj->next = lista; lista = uj; return uj; }
    //***********************************************************************
    void clear() {
    //***********************************************************************
        // lista felszabadítása
        for (cella_lista *pi = lista; pi != nullptr;) {
            cella_lista *temp = pi->next;
            delete pi;
            pi = temp;
        }
        lista = nullptr;
    }
};


//***********************************************************************
struct simulation;
struct fenyut;
class vezetes_tomb_tipus;
//***********************************************************************


//***********************************************************************
class sugarfeldolgozo {
//***********************************************************************
    tomb<kulso_oldal> kulso_oldalak;
    tomb<sugar_lista> kek_sugarlista, sarga_sugarlista;
    global_sugar_tomb global_kek, global_sarga;
    uns dir_mode, ray_per_cell_dir;
    dbl cut_level;
    dbl d_light_powder, dlppb, dlppy;
    dbl blue_tenyleges, yellow_tenyleges;
    dbl yellow_correction;
    bool is_sarga_szetmegy;
    bool is_top_junction; 
    char output_side;
    uns junction_layer;
    //***********************************************************************
    void clear() {
    //***********************************************************************
        kulso_oldalak.clear();
        global_kek.clear();
        global_sarga.clear();
    }
    //***********************************************************************
    dbl get_blue_weighted_length_and_recalculate_blue_P() {
    //***********************************************************************
        dbl sum_dP = 0, sum_P = 0;
        for (uns i=0; i<kek_sugarlista.size(); i++) {
            kek_sugarlista[i].recalculate_P_kek(sum_dP, sum_P);
        }
        return sum_dP / sum_P;
    }
    //***********************************************************************
    dbl get_yellow_weighted_length_and_recalculate_yellow_P() {
    //***********************************************************************
        dbl sum_dP = 0, sum_P = 0;
        for (uns i=0; i<sarga_sugarlista.size(); i++) {
            sarga_sugarlista[i].recalculate_P_kek(sum_dP, sum_P);
        }
        return sum_dP / sum_P;
    }
    //***********************************************************************
    void kulso_oldalak_gyujtese(const simulation & aktSim, const tomb3d<t_modell_cella> & modell_racs);
    dbl kek_sugarak_gyujtese(const simulation & aktSim, const tomb3d<t_modell_cella> & modell_racs);
    dbl sarga_sugarak_gyujtese(const simulation & aktSim, const tomb3d<t_modell_cella> & modell_racs);
    dbl build_blue_cell_rays(const simulation & aktSim, tomb3d<t_modell_cella> & modell_racs, FILE * fp, vezetes_tomb_tipus & vt, uns & db_cell);
    dbl build_yellow_cell_rays(const simulation & aktSim, tomb3d<t_modell_cella> & modell_racs, FILE * fp, vezetes_tomb_tipus & vt, uns & db_cell);
    dbl get_cell_list(bool is_blue, bool is_half, tomb<sugar_cella> & uj_lista, sugar_adat * p_sugar, tomb3d<t_modell_cella> & modell_racs, dbl P0);
    void write_one_path(const simulation & aktSim, const tomb3d<t_modell_cella> & modell_racs, FILE * fp, vezetes_tomb_tipus & vt, const tomb<sugar_cella> & cellatomb, bool is_blue, bool is_first);
    //***********************************************************************
public:
    sugarfeldolgozo() :dir_mode{ 0 }, ray_per_cell_dir{ 1 }, cut_level{ 0 }, d_light_powder{ 0 }, is_sarga_szetmegy{ false }, is_top_junction{ false },
        junction_layer{ 0 }, dlppb{ 0 }, dlppy{ 0 }, output_side{ 'T' }, blue_tenyleges{ 1 }, yellow_tenyleges{ 1 }, yellow_correction{ 0 } { }
    void build_and_write(const simulation & aktSim, tomb3d<t_modell_cella> & modell_racs, uns dir_mode, char output_side, 
        uns ray_per_cell_dir, bool is_sarga_szetmegy, dbl cut_level, dbl d_light_powder, bool is_top_junction, uns junction_layer, 
        FILE * fp, vezetes_tomb_tipus & vt, dbl blue_tenyleges, dbl yellow_tenyleges, dbl yellow_correction);
    dbl get_d_light_powder_proportional_blue()const { return dlppb; }
    dbl get_d_light_powder_proportional_yellow()const { return dlppy; }
};


//***********************************************************************
#endif
//***********************************************************************
