#include "adattipusok.h"
#include "bemenet.h"
#include "fajlolvasas_segito_rutinok.h"


//***********************************************************************
namespace ns_v6core {
//***********************************************************************


//***********************************************************************
void fazis_broken_line::beolvas_fajlbol(const fazisvalto & fv){
//***********************************************************************
    set_hiba_hol h("fazis_broken_line::beolvas_fajlbol");
    uns db = fajl::get_uns("phase broken line length");
    bool volt_mar_F = false;
    for (uns i = 0; i < db; i++) {
        int ch;
        rvt T, v;
        switch (ch = fajl::get_char("phase broken line T/F/L")) {
            case 'T':  
                T = fajl::get_rvt("phase broken line T temp", true); 
                v = fajl::get_rvt("phase broken line T value", true);
                if (volt_mar_F)
                    gorbe_1.add_pont(T, v);
                else
                    gorbe_0.add_pont(T, v);
                break;
            case 'F':
                T = volt_mar_F ? fv.get_F2() : fv.get_F1();
                v = fajl::get_rvt("phase broken line F value", true);
                if (volt_mar_F)
                    gorbe_1.add_pont(T, v);
                else
                    gorbe_0.add_pont(T, v);
                volt_mar_F = true;
                break;
            case 'L':
                T = volt_mar_F ? fv.get_L2() : fv.get_L1();
                v = fajl::get_rvt("phase broken line L value", true);
                if (volt_mar_F)
                    gorbe_1.add_pont(T, v);
                else
                    gorbe_0.add_pont(T, v);
                volt_mar_F = true;
                break;
            default:  throw hiba(1, "unknown phase broken line property T/F/L (%c)", ch);
        }
    }
}


//***********************************************************************
void broken_line::beolvas_fajlbol(const fazisvalto & fv){
//***********************************************************************
    set_hiba_hol h("broken_line::beolvas_fajlbol");
    uns db = fajl::get_uns("broken line length");
    bool volt_mar_F = false;
    for (uns i = 0; i < db; i++) {
        int ch;
        rvt T, v;
        switch (ch = fajl::get_char("broken line T/F/L")) {
            case 'T':  
                T = fajl::get_rvt("broken line T temp", true); 
                v = fajl::get_rvt("broken line T value", true);
                add_pont(T, v);
                break;
            case 'F':
                if(!fv.is_fazisvalto())
                    throw hiba(1, "F property in broken line type parameter of non phase change material", ch);
                T = volt_mar_F ? fv.get_F2() : fv.get_F1();
                v = fajl::get_rvt("broken line F value", true);
                add_pont(T, v);
                volt_mar_F = true;
                break;
            case 'L':
                if (!fv.is_fazisvalto())
                    throw hiba(1, "L property in broken line type parameter of non phase change material", ch);
                T = volt_mar_F ? fv.get_L2() : fv.get_L1();
                v = fajl::get_rvt("broken line L value", true);
                add_pont(T, v);
                volt_mar_F = true;
                break;
            default:  throw hiba(1, "unknown broken line property T/F/L (%c)", ch);
        }
    }
}

}