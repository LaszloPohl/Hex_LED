//***********************************************************************
// apa class source
// Creation date:  2009. 07. 11.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#include "apa.h"
#include "gfunc.h"
#include "sugarak.h"
//***********************************************************************

//***********************************************************************
uns semiconductor::db=0;
//***********************************************************************
//***********************************************************************
apa::apa(const char * ProjectFile) 
    :projFile(ProjectFile), modell_fa_el{ nullptr }, modell_fa_th{ nullptr }, modell_fa_elth{ nullptr } {
//***********************************************************************
    path = getPath(projFile);
    proj_nev = getFileNameWithoutExtension(projFile);
    srfajl fajl;
    resetAnalLevel(4);
    logprint("Open project: %s",ProjectFile);
    fajl.open(projFile);
    if(fajl.lines().size()<1||fajl.lines()[0][0].LowCase()!="vsun3-project")
        hiba("apa::apa","hibas projekt fajl");
    else{
        
        // modellfájlok nevei
        
        uns db=0,m=0;
        for(uns i=1;i<fajl.lines().size();i++)if(fajl.lines()[i][0].LowCase()=="model")db++;
        if(!db)throw hiba("apa::apa()","no models are defined in %s",ProjectFile);
        tmodels.clear();
        tmodels.resize(db);
        db=0;
        for(uns i=1;i<fajl.lines().size();i++)if(fajl.lines()[i][0].LowCase()=="model"){
            tmodels[db].fileName=fajl.lines()[i][1];
            tmodels[db].mod_nev = getFileNameWithoutExtension(tmodels[db].fileName);
            tmodels[db].simdb=0;
            db++;
        }
        
        // szimulációs fájlok nevei
        
        db=0;
        for(uns i=1;i<fajl.lines().size();i++)if(fajl.lines()[i][0].LowCase()=="simulation")db++;
        if(!db)throw hiba("apa::apa()","no simulations are defined in %s",ProjectFile);
        tsim.clear();
        tsim.resize(db);
        db=0;
        for(uns i=1;i<fajl.lines().size();i++)if(fajl.lines()[i][0].LowCase()=="simulation"){
            if(m==0)throw hiba("apa::apa()","simulation before model in %s",ProjectFile);
            tsim[db].fileName=fajl.lines()[i][1];
            tsim[db].sim_nev = getFileNameWithoutExtension(tsim[db].fileName);
            tmodels[m-1].simdb++;
            tsim[db].pmodel=&tmodels[m-1];
            db++;
        }
        else if(fajl.lines()[i][0].LowCase()=="model")m++;

        // model fájlok olvasása

        for(uns i=0;i<tmodels.size();i++)tmodels[i].read(path);

        // szimulációs fájlok olvasása

        for(uns i=0;i<tsim.size();i++)tsim[i].read(path);
        uns sumdb=0;
        for(uns i=0;i<tsim.size();i++){ // hány szimuláció van, a többszörös gerjesztéseket figyelembe véve
            uns sim_darab=1,vane=0;
            for(uns j=0;j<4;j++)tsim[i].index_temp[j]=colmax;
            for(uns j=0;j<colmax;j++){
                if(tsim[i].mulE[j].size()>0){
                    vane++;
                    sim_darab*=tsim[i].mulE[j].size();
                    if(tsim[i].index_temp[0]==colmax)tsim[i].index_temp[0]=j;
                    else tsim[i].index_temp[1]=j;
                }
                if(tsim[i].mulT[j].size()>0){
                    vane++;
                    sim_darab*=tsim[i].mulT[j].size();
                    if(tsim[i].index_temp[2]==colmax)tsim[i].index_temp[2]=j;
                    else tsim[i].index_temp[3]=j;
                }
                if(vane>2)throw hiba("apa::apa()","maximum 2 multiple excitation is allowed in %s",ProjectFile);
            }
            if(vane){ sumdb+=sim_darab; tsim[i].db_temp=sim_darab; }
            else { sumdb++; tsim[i].db_temp=1; }
        }
        tomb<simulation> sim_temp;
        sim_temp.resize(sumdb);
        uns j=0;
        for(uns i=0;i<tsim.size();i++){// többszörös szimulációk létrehozása
            cuns db=tsim[i].db_temp;
            cuns * ix=tsim[i].index_temp;
            for(uns k=0;k<db;k++)sim_temp[j+k]=tsim[i];
            if(db>1){ // van többszörös
                if(ix[0]!=colmax){
                    if(tsim[i].mulE[ix[0]].size()==db){ // egy elektromos gerjesztés értéke változik
                        for(uns k=0;k<db;k++){
                            sim_temp[j+k].texcitE[ix[0]]=tsim[i].mulE[ix[0]][k];
                            char s[100];
                            sprintf(s,"_%g",sim_temp[j+k].texcitE[ix[0]].ertek);
                            sim_temp[j+k].name+=s;
                        }
                    }
                    else if(ix[1]!=colmax){ // két elektomos multi gerjesztés van
                        uns db1=tsim[i].mulE[ix[0]].size();
                        uns db2=tsim[i].mulE[ix[1]].size();
                        if(db1*db2!=db)throw hiba("apa::apa()","Impossibility happend, program error. (db)");
                        for(uns k=0;k<db1;k++)
                            for(uns L=0;L<db2;L++){
                                sim_temp[j+k*db2+L].texcitE[ix[0]]=tsim[i].mulE[ix[0]][k];
                                sim_temp[j+k*db2+L].texcitE[ix[1]]=tsim[i].mulE[ix[1]][L];
                                char s[100];
                                sprintf(s,"_%g_%g",tsim[i].mulE[ix[0]][k].ertek,tsim[i].mulE[ix[1]][L].ertek);
                                sim_temp[j+k*db2+L].name+=s;
                            }
                    }
                    else if(ix[2]!=colmax){ // vegyes gerjesztés van
                        uns db1=tsim[i].mulE[ix[0]].size();
                        uns db2=tsim[i].mulT[ix[2]].size();
                        if(db1*db2!=db)throw hiba("apa::apa()","Impossibility happend, program error. (db)");
                        for(uns k=0;k<db1;k++)
                            for(uns L=0;L<db2;L++){
                                sim_temp[j+k*db2+L].texcitE[ix[0]]=tsim[i].mulE[ix[0]][k];
                                sim_temp[j+k*db2+L].texcitT[ix[2]]=tsim[i].mulT[ix[2]][L];
                                char s[100];
                                sprintf(s,"_%g_%g",tsim[i].mulE[ix[0]][k].ertek,tsim[i].mulT[ix[2]][L].ertek);
                                sim_temp[j+k*db2+L].name+=s;
                            }
                    }
                    else throw hiba("apa::apa()","Impossibility happend, program error. (ix2)");
               }
                else if(ix[2]!=colmax){ // csak termikus gerj van
                    if(tsim[i].mulT[ix[2]].size()==db){ // egy termikus gerjesztés értéke változik
                        for(uns k=0;k<db;k++){
                            sim_temp[j+k].texcitT[ix[2]]=tsim[i].mulT[ix[2]][k];
                            char s[100];
                            sprintf(s,"_%g",sim_temp[j+k].texcitT[ix[0]].ertek);
                            sim_temp[j+k].name+=s;
                        }
                    }
                    else if(ix[3]!=colmax){ // két termikus multi gerjesztés van
                        uns db1=tsim[i].mulT[ix[2]].size();
                        uns db2=tsim[i].mulT[ix[3]].size();
                        if(db1*db2!=db)throw hiba("apa::apa()","Impossibility happend, program error. (db)");
                        for(uns k=0;k<db1;k++)
                            for(uns L=0;L<db2;L++){
                                sim_temp[j+k*db2+L].texcitT[ix[2]]=tsim[i].mulT[ix[2]][k];
                                sim_temp[j+k*db2+L].texcitT[ix[3]]=tsim[i].mulT[ix[3]][L];
                                char s[100];
                                sprintf(s,"_%g_%g",tsim[i].mulT[ix[2]][k].ertek,tsim[i].mulT[ix[3]][L].ertek);
                                sim_temp[j+k*db2+L].name+=s;
                            }
                    }
                    else throw hiba("apa::apa()","Impossibility happend, program error. (ix3)");
                }
                else throw hiba("apa::apa()","Impossibility happend, program error. (nix)");
            }
            j+=db;
        }
        tsim=sim_temp;
        sumdb=tsim.size();
        for(uns i=0;i<tsim.size();i++) // több ambient hõmérséklet van-e?
            sumdb += tsim[i].mulAmbiT.size();
        if(sumdb!=tsim.size()){
            sim_temp.resize(sumdb);
            j=0;
            for(uns i=0;i<tsim.size();i++){// többszörös szimulációk létrehozása különbözõ ambient hõmérsékletre
                cuns db=tsim[i].mulAmbiT.size()+1;
                for(uns k=0;k<db;k++){
                    sim_temp[j+k]=tsim[i];
                    if(db>1){
                        char s[100];
                        if( k == 0 )
                            sprintf(s,"_%gC",sim_temp[j+k].ambiT);
                        else 
                            sprintf(s,"_%gC",sim_temp[j+k].mulAmbiT[k-1]);
                        sim_temp[j+k].name+=s;
                        if( k > 0 )
                            sim_temp[j+k].ambiT=sim_temp[j+k].mulAmbiT[k-1];
                    }
                }
                j+=db;
            }
        }
        tsim=sim_temp;
    }
    logprint("Opening done");
}


//***********************************************************************
void apa::write_v6sim(){
//***********************************************************************
    PLString v6_nev = path + proj_nev + ".hex";
    FILE *fp = fopen(v6_nev.c_str(), "wt");
    if (fp == nullptr)
        throw hiba("apa::write_v6sim", "cannot open %s", v6_nev.c_str());
    printf("creating file: %s\n", v6_nev.c_str());
    fprintf(fp, "V6SIMFILE\n");
    if (tsim.size() > 0) {
        fprintf(fp, "CPUTHREADS%u\n", tsim[0].cpu_threads);
    }
    uns analizisszam = 0;
    for (uns i = 0; i < tsim.size(); i++) {
        analizisszam += tsim[i].tanal.size();
    }
    fprintf(fp, "NS%u\n", analizisszam);
    uns akt_n = 0;
    for (uns i = 0; i < tsim.size(); i++) {
        bool is_uj_model = i == 0 || tsim[i].pmodel != tsim[i - 1].pmodel;
        for (uns j = 0; j < tsim[i].tanal.size(); j++) {
            PLString leiras = proj_nev;
            if (tmodels.size()>1)
                leiras = leiras + "_(" + tsim[i].pmodel->name + ")";
            if (tsim.size()>tmodels.size())
                leiras = leiras + "_(" + tsim[i].name + ")";
            if (tsim[i].tanal.size()>1)
                leiras = leiras + "_(" + tsim[i].tanal[j].nev + ")";
            leiras.replace('-', '_');
            leiras.replace(' ', '_');
            leiras.replace('\"', '_');
            akt_n++;
            write_akt_sim(fp, tsim[i], tsim[i].tanal[j], leiras, akt_n);
        }
    }
    fprintf(fp, "EV\n");
    fclose(fp);
}

FILE *logfile = nullptr;

//***********************************************************************
void apa::write_akt_sim(FILE *fp, simulation & aktSim, analysis & aktAnal, const PLString & leiras, uns akt_n){
// beleír egy teljes modellt + szimulációt a hex fájlba
//***********************************************************************
    aktSim.fimdesc_name = leiras + "_" + idobelyeg;
    fprintf(fp, "BS%u\"%s\"\n", akt_n, aktSim.fimdesc_name.c_str());
    switch (aktSim.mezo) {
        case FieldEl:
            fprintf(fp, "FE%u\n", aktSim.mezo_szamitasi_mod > 1 ? 2 : 1);
            break;
        case FieldTherm:
            fprintf(fp, "FT%u\n", aktSim.mezo_szamitasi_mod > 1 ? 2 : 1);
            break;
        case FieldElTherm:
            fprintf(fp, "FET%u\n", aktSim.mezo_szamitasi_mod > 1 ? 2 : 1);
            break;
    }
    logfile = fopen((path + aktSim.fimdesc_name+".ide").c_str(), "wt");
    write_materials(fp, aktSim);
    write_colors(fp, aktSim);
    uns junction_db = write_junctions(fp, aktSim);
    // write_special_boundaries(fp, aktSim); // Lehet, hogy inkább ki kellene fejteni, face-enként menteni.
    write_boundary_conditions(fp, aktSim);
    most("cell grid generation starts");
    printf("cell grid generation starts\n");
    uns csatlakozo_db = 0;
    uns cellaszam = build_modell_racs(aktSim, csatlakozo_db);
    printf("cell indexing and path finding starts\n");
    index_and_write_elemi_cells(fp, aktSim, cellaszam, csatlakozo_db);
    most("cells written, buildig tree starts");
    printf("cells written, buildig tree starts\n");
    write_meretek(fp);
    build_modell_fa(aktSim);
    write_modell_tree(fp, aktSim);
    if (sugar_feldolgozo.get_d_light_powder_proportional_blue() > 0) {
        fprintf(fp, "DLPPB%.6g;\n", sugar_feldolgozo.get_d_light_powder_proportional_blue());
    }
    if (sugar_feldolgozo.get_d_light_powder_proportional_yellow() > 0) {
        fprintf(fp, "DLPPY%.6g;\n", sugar_feldolgozo.get_d_light_powder_proportional_yellow());
    }
    write_analizisek(fp, aktSim, junction_db);
    fprintf(fp, "ES%u\n", akt_n);
    most("sim write is ready");
    printf("sim write is ready\n");
    esemenyek_kiirasa();
    fclose(logfile);
    logfile = nullptr;
}


//***********************************************************************
void apa::write_materials(FILE *fp, simulation & aktSim){
//***********************************************************************
    uns db = 0;
    for (uns i = 0; i < aktSim.pmodel->tmat.size(); i++)
        if (aktSim.pmodel->tmat[i].is_isotrop())
            db++;
        else db += 3;
    fprintf(fp, "NM%u\n", db);
    uns anyagindex = 1;
    for (uns i = 0; i < aktSim.pmodel->tmat.size(); i++)
        if (aktSim.pmodel->tmat[i].is_isotrop()) {
            aktSim.pmodel->tmat[i].write(fp, 0, anyagindex);
            anyagindex++;
        }
        else {
            aktSim.pmodel->tmat[i].write(fp, 0, anyagindex);
            anyagindex++;
            aktSim.pmodel->tmat[i].write(fp, 1, anyagindex);
            anyagindex++;
            aktSim.pmodel->tmat[i].write(fp, 2, anyagindex);
            anyagindex++;
        }
    if (anyagindex != db + 1)
        hiba("apa::write_materials", "db+1!=anyagindex");
}


//***********************************************************************
void apa::write_colors(FILE * fp, simulation & aktSim){
//***********************************************************************
    uns db = 0;
    for (uns i = 0; i < colmax; i++)
        if (aktSim.pmodel->tcolor[i].is && aktSim.pmodel->tcolor[i].tipus == SzinNormal && aktSim.pmodel->tcolor[i].terfogat>0.0) {
            db++;
            aktSim.pmodel->tcolor[i].index = db;
        }
        else aktSim.pmodel->tcolor[i].index = 0;
    fprintf(fp, "NC%u\n", db);
    for (uns i = 0; i < colmax; i++)
        if (aktSim.pmodel->tcolor[i].is && aktSim.pmodel->tcolor[i].tipus == SzinNormal && aktSim.pmodel->tcolor[i].terfogat>0.0) {
            fprintf(fp, "C%uV%g\n", aktSim.pmodel->tcolor[i].index, aktSim.pmodel->tcolor[i].terfogat);
            //fprintf(fp, "PM%u;%u;%u;\n", aktSim.pmodel->tcolor[i].pmat->anyagindex[0], aktSim.pmodel->tcolor[i].pmat->anyagindex[1], aktSim.pmodel->tcolor[i].pmat->anyagindex[2]);
        }
}

//***********************************************************************
uns apa::write_junctions(FILE * fp, simulation & aktSim){
//***********************************************************************
    uns db = 0;
    for (uns i = 0; i < colmax; i++)
        if (aktSim.pmodel->tcolor[i].is && aktSim.pmodel->tcolor[i].tipus == SzinNormal && aktSim.pmodel->tcolor[i].terfogat>0.0) {
            for (uns j = 0; j < aktSim.pmodel->tcolor[i].tsemi.size(); j++)
                aktSim.pmodel->tcolor[i].tsemi[j].index = db + j + 1;
            db+= aktSim.pmodel->tcolor[i].tsemi.size();
        }
    fprintf(fp, "NJ%u\n", db);
    for (uns i = 0; i < colmax; i++)
        if (aktSim.pmodel->tcolor[i].is && aktSim.pmodel->tcolor[i].tipus == SzinNormal && aktSim.pmodel->tcolor[i].terfogat>0.0) {
            for (uns j = 0; j < aktSim.pmodel->tcolor[i].tsemi.size(); j++) {
                fprintf(fp, "BJ%u\n", aktSim.pmodel->tcolor[i].tsemi[j].index);
                fprintf(fp, "PGJ");
                aktSim.pmodel->tcolor[i].tsemi[j].par.write_semi(fp);
                fprintf(fp, "PDC");
                aktSim.pmodel->tcolor[i].tsemi[j].D.write_normal(fp, 0);
                fprintf(fp, "PRC");
                aktSim.pmodel->tcolor[i].tsemi[j].R.write_normal(fp, 0);
                fprintf(fp, "PRA");
                aktSim.pmodel->tcolor[i].tsemi[j].rad.write_sem_rad_lum(fp);
                fprintf(fp, "PLU");
                aktSim.pmodel->tcolor[i].tsemi[j].lum.write_sem_rad_lum(fp);
                fprintf(fp, "PAS%g\n", aktSim.pmodel->tcolor[i].tsemi[j].As);
                fprintf(fp, "EJ%u\n", aktSim.pmodel->tcolor[i].tsemi[j].index);
            }
        }
    return db;
}


//***********************************************************************
void write_bou_tip(FILE * fp, ConvTipus tipus) {
//***********************************************************************
    switch (tipus) {
    case ConvHTC:
        fprintf(fp, "HTC");
        break;
    case ConvVertical_1:
        fprintf(fp, "VE1");
        break;
    case ConvUpper_1:
        fprintf(fp, "U1P");
        break;
    case ConvLower_1:
        fprintf(fp, "L1W");
        break;
    case ConvVerticalChurchill_P_1:
        fprintf(fp, "C1P");
        break;
    case ConvVerticalChurchill_P_2:
        fprintf(fp, "C2P");
        break;
    case ConvVerticalChurchill_T:
        fprintf(fp, "C3T");
        break;
    case ConvVerticalChurchill_C:
        fprintf(fp, "C4C");
        break;
    case ConvVerticalLee_T:
        fprintf(fp, "LET");
        break;
    case ConvVerticalLee_P:
        fprintf(fp, "LEP");
        break;
    case ConvVerticalMihajev:
        fprintf(fp, "MIH");
        break;
    case ConvYovanovichMin:
        fprintf(fp, "YMI");
        break;
    case ConvYovanovichMax:
        fprintf(fp, "YMA");
        break;
    case ConvWei:
        fprintf(fp, "WEI");
        break;
    default:
        throw hiba("write_bou_tip", "unsupported boundary type");
    }
}



//***********************************************************************
void apa::write_special_boundaries(FILE * fp, simulation & aktSim) {
//***********************************************************************
    fprintf(fp, "NB%u\n", aktSim.tconv.size() + aktSim.tucha.size());
    for (uns i = 0; i < aktSim.tconv.size(); i++) {
        fprintf(fp, "BB%u\n", i + 1);
        fprintf(fp, "TC\n");
        fprintf(fp, "RA");
        aktSim.tconv[i].radiation.write_const(fp);
        fprintf(fp, "CV");
        write_bou_tip(fp, aktSim.tconv[i].vertical_tipus);
        aktSim.tconv[i].vertical_value.write_const(fp);
        fprintf(fp, "CU");
        write_bou_tip(fp, aktSim.tconv[i].upper_tipus);
        aktSim.tconv[i].upper_value.write_const(fp);
        fprintf(fp, "CL");
        write_bou_tip(fp, aktSim.tconv[i].lower_tipus);
        aktSim.tconv[i].lower_value.write_const(fp);
        if (aktSim.tconv[i].angle != 0) {
            if (aktSim.tconv[i].axis == X_IRANY)
                fprintf(fp, "XX\n");
            else if (aktSim.tconv[i].axis == Y_IRANY)
                fprintf(fp, "XY\n");
            else throw hiba("apa::write_special_boundaries", "axis cannot be Z");
            fprintf(fp, "A%g\n", aktSim.tconv[i].angle);
        }
        if (aktSim.tconv[i].edge!=ConvHTC) {
            switch (aktSim.tconv[i].edge) {
                case ConvEdgeWei_H:
                    fprintf(fp, "CEW_H\n");
                    break;
                case ConvEdgeWei_I:
                    fprintf(fp, "CEW_I\n");
                    break;
                case ConvEdgeWei_HH:
                    fprintf(fp, "CEWHH\n");
                    break;
                case ConvEdgeWei_HI:
                    fprintf(fp, "CEWHI\n");
                    break;
                default:
                    throw hiba("apa::write_special_boundaries", "unsupported Wei boundary type");
            }
        }
        fprintf(fp, "EB%u\n", i + 1);
    }
    for (uns i = aktSim.tconv.size(); i < aktSim.tconv.size() + aktSim.tucha.size(); i++) {
        fprintf(fp, "BB%u\n", i + 1);
        fprintf(fp, "TG\n");
        fprintf(fp, "TODO:uchannel_mentese");
    }
}

//***********************************************************************
void apa::write_boundary_conditions(FILE * fp, simulation & aktSim) {
//***********************************************************************
    fprintf(fp, "NP%u\n", aktSim.peremlista_menteshez.size());
    for (uns i = 0; i < aktSim.peremlista_menteshez.size(); i++) {
        fprintf(fp, "P%u", aktSim.peremlista_menteshez[i].v6_index);
        switch (aktSim.peremlista_menteshez[i].tipus) {
            case PeremOpen:
                throw hiba("apa::write_boundary_conditions", "PeremOpen nem lehetne a tombben...");
                break;
            case PeremV:
                fprintf(fp, "U%g\n", aktSim.peremlista_menteshez[i].value);
                break;
            case PeremT:
                fprintf(fp, "T%g\n", aktSim.peremlista_menteshez[i].value);
                break;
            case PeremR:
                if (aktSim.peremlista_menteshez[i].conv.is_defined) {
                    fprintf(fp, "SC\n"); // speciális konvektív (+rad), TODO
                }
                else if (aktSim.peremlista_menteshez[i].conv_map.x_size()>0) {
                    fprintf(fp, "M\n"); // HTC map, TODO
                }
                else 
                    fprintf(fp, "H%g\n", aktSim.peremlista_menteshez[i].value);
                break;
            case PeremRU:
                fprintf(fp, "K%gT%g\n", aktSim.peremlista_menteshez[i].value, aktSim.peremlista_menteshez[i].value2);
                break;
        }
    }
}


//***********************************************************************
struct perembeallito {
//***********************************************************************
    Oldal oldal;
    bool is_rendes_perem;
    char c;
    uns x, y;
    const color *szomszed;
};


//***********************************************************************
uns apa::build_modell_racs(simulation & aktSim, uns & csatlakozo_face_db) {
//***********************************************************************
    const char * fvnev = "apa::build_modell_racs()";
    model & aktMod = *aktSim.pmodel;
    cuns x_res = aktMod.x_res;
    cuns y_res = aktMod.y_res;
    cuns z_res = aktMod.z_res;
    modell_racs.resize(x_res, y_res, z_res);

    perembeallito b[6];
    b[0].oldal = WEST;
    b[0].c = 'W';
    b[1].oldal = EAST;
    b[1].c = 'E';
    b[2].oldal = SOUTH;
    b[2].c = 'S';
    b[3].oldal = NORTH;
    b[3].c = 'N';
    b[4].oldal = BOTTOM;
    b[4].c = 'B';
    b[5].oldal = TOP;
    b[5].c = 'T';

    // 1. kör: cellák lényegének létrehozása, anyag, méret hozzárendelés, nagyfelbontásúra bontás, ahol kell

    for (uns z = 0; z < z_res; z++)
        for (uns y = 0; y < y_res; y++)
            for (uns x = 0; x < x_res; x++) {

                cuns cella_szin = aktMod.tbmp[z].getpixel_also(x, y);
                const color & cella_color = aktMod.tcolor[cella_szin];
                t_modell_cella & aktCella = modell_racs.getref(x, y, z);

                if (cella_color.tipus != SzinNormal) {
                    aktCella.is_cella = false;
                    continue;
                }

                aktCella.is_cella = true;
                aktCella.color_index = cella_color.index;
                aktCella.anyag_index = cella_color.pmat->anyagindex[0];
                aktCella.pmat = cella_color.pmat;
                aktCella.is_nonlin_el = !cella_color.pmat->is_lin_el();
                aktCella.is_nonlin_th = !cella_color.pmat->is_lin_th();

                // méretek

                cd x_pit = aktMod.x_pit[x].get(nulla);
                cd y_pit = aktMod.y_pit[y].get(nulla);
                cd z_pit = aktMod.z_pit[z].get(aktMod.x_hossz[2 * x]);
                cd A_x = y_pit*z_pit;
                cd A_y = x_pit*z_pit;
                cd A_z = x_pit*y_pit;

                aktCella.V = A_x*x_pit;

                aktCella.sarkok.x0 = (x == 0) ? 0 : aktMod.x_hossz[2 * x - 1];
                aktCella.sarkok.x1 = aktMod.x_hossz[2 * x + 1];
                aktCella.sarkok.y0 = (y == 0) ? 0 : aktMod.y_hossz[2 * y - 1];
                aktCella.sarkok.y1 = aktMod.y_hossz[2 * y + 1];
                aktCella.sarkok.z0 = (z == 0) ? 0 : aktMod.z_hossz[2 * z - 1];
                aktCella.sarkok.z1 = aktMod.z_hossz[2 * z + 1];
                aktCella.sarkok.x0 = round(aktCella.sarkok.x0*1e12)*1e-12;
                aktCella.sarkok.x1 = round(aktCella.sarkok.x1*1e12)*1e-12;
                aktCella.sarkok.y0 = round(aktCella.sarkok.y0*1e12)*1e-12;
                aktCella.sarkok.y1 = round(aktCella.sarkok.y1*1e12)*1e-12;
                aktCella.sarkok.z0 = round(aktCella.sarkok.z0*1e12)*1e-12;
                aktCella.sarkok.z1 = round(aktCella.sarkok.z1*1e12)*1e-12;

                aktCella.face_adat[WEST].A = aktCella.face_adat[EAST].A = A_x;
                aktCella.face_adat[SOUTH].A = aktCella.face_adat[NORTH].A = A_y;
                aktCella.face_adat[BOTTOM].A = aktCella.face_adat[TOP].A = A_z;
                aktCella.face_adat[WEST].L = aktCella.face_adat[EAST].L = 0.5*x_pit;
                aktCella.face_adat[SOUTH].L = aktCella.face_adat[NORTH].L = 0.5*y_pit;
                aktCella.face_adat[BOTTOM].L = aktCella.face_adat[TOP].L = 0.5*z_pit;
                aktCella.face_adat[WEST].anyag_index = aktCella.face_adat[EAST].anyag_index = cella_color.pmat->anyagindex[0];
                aktCella.face_adat[SOUTH].anyag_index = aktCella.face_adat[NORTH].anyag_index = cella_color.pmat->anyagindex[1];
                aktCella.face_adat[BOTTOM].anyag_index = aktCella.face_adat[TOP].anyag_index = cella_color.pmat->anyagindex[2];

                aktCella.face_adat[WEST].th_perem_c = 'W';
                aktCella.face_adat[EAST].th_perem_c = 'E';
                aktCella.face_adat[SOUTH].th_perem_c = 'S';
                aktCella.face_adat[NORTH].th_perem_c = 'N';
                aktCella.face_adat[BOTTOM].th_perem_c = 'B';
                aktCella.face_adat[TOP].th_perem_c = 'T';

                aktCella.is_el = (cella_color.field != FieldTherm) && (aktSim.mezo != FieldTherm);
                aktCella.is_th = (cella_color.field != FieldEl) && (aktSim.mezo != FieldEl);

                // perem-e vagy másik cellához kapcsolódik?

                b[0].is_rendes_perem = x == 0;
                b[1].is_rendes_perem = x == aktMod.x_res - 1;
                b[2].is_rendes_perem = y == 0;
                b[3].is_rendes_perem = y == aktMod.y_res - 1;
                b[4].is_rendes_perem = z == 0;
                b[5].is_rendes_perem = z == aktMod.z_res - 1;

                b[0].szomszed = b[0].is_rendes_perem ? &cella_color : &(aktMod.tcolor[aktMod.tbmp[z].getpixel_also(x - 1, y)]);
                b[1].szomszed = b[1].is_rendes_perem ? &cella_color : &(aktMod.tcolor[aktMod.tbmp[z].getpixel_also(x + 1, y)]);
                b[2].szomszed = b[2].is_rendes_perem ? &cella_color : &(aktMod.tcolor[aktMod.tbmp[z].getpixel_also(x, y - 1)]);
                b[3].szomszed = b[3].is_rendes_perem ? &cella_color : &(aktMod.tcolor[aktMod.tbmp[z].getpixel_also(x, y + 1)]);
                b[4].szomszed = b[4].is_rendes_perem ? &cella_color : &(aktMod.tcolor[aktMod.tbmp[z - 1].getpixel_also(x, y)]);
                b[5].szomszed = b[5].is_rendes_perem ? &cella_color : &(aktMod.tcolor[aktMod.tbmp[z + 1].getpixel_also(x, y)]);

                b[0].x = y; b[0].y = z;
                b[1].x = y; b[1].y = z;
                b[2].x = x; b[2].y = z;
                b[3].x = x; b[3].y = z;
                b[4].x = x; b[4].y = y;
                b[5].x = x; b[5].y = y;

                for (uns i = 0; i < 6; i++) {
                    if (b[i].is_rendes_perem || b[i].szomszed->tipus != SzinNormal) { // mellette biztosan perem van
                        // a mikrocsatorna perem nincs kezelve sem nemlinearitás, sem egyéb szempontból
                        if (aktCella.is_el) {
                            const boundary * pb = aktSim.get_perem(x, y, z, b[i].oldal, true);
                            aktCella.set_face_el(b[i].oldal, true, pb->v6_index, 0, pb->conv.is_defined);
                        }
                        else aktCella.set_face_el(b[i].oldal, false, 0, 0);
                        if (aktCella.is_th) {
                            const boundary * pb = aktSim.get_perem(x, y, z, b[i].oldal, false);
                            aktCella.set_face_th(b[i].oldal, true, pb->v6_index, 0, pb->is_special ? b[i].x : ~0, b[i].y, b[i].c, pb->conv.is_defined);
                        }
                        else aktCella.set_face_th(b[i].oldal, false, 0, 0, ~0, b[i].y, b[i].c);
                    }
                    else { // a b[i].oldal mellett nem perem és nem belsõ perem van, de nem szimulált tér még lehet
                        if (aktCella.is_el) {
                            //if (b[i].szomszed->field == FieldTherm) // a szomszéd cella elektromosan nem szimulált, szigetelõ!
                            if (b[i].szomszed->field == FieldTherm || aktSim.mezo == FieldTherm) // a szomszéd cella elektromosan nem szimulált, szigetelõ!
                                aktCella.set_face_el(b[i].oldal, true, 0, 0);
                            else { // igazi szomszédja van
                                aktCella.set_face_el(b[i].oldal, false, 0, 1);
                                for (uns j = 0; j < cella_color.tsemi.size(); j++)
                                    if (cella_color.tsemi[j].col2 == aktCella.get_szomszed_color(aktMod, x, y, z, b[i].oldal)) {
                                        aktCella.face_adat[b[i].oldal].junction_index = cella_color.tsemi[j].index;
                                        aktCella.is_nonlin_el = true; // nehéz lenne olyan junctiont találni, aminek kostans a vezetése...
                                        aktCella.is_junction = true;
                                        if (b[i].oldal == BOTTOM)   aktCella.junction_bottom_face = 1000;
                                        if (b[i].oldal == TOP)      aktCella.junction_top_face = 1001;
                                    }
                            }
                        }
                        else aktCella.set_face_el(b[i].oldal, false, 0, 0); // nem szimuláljuk elektromosan
                        if (aktCella.is_th) {
                            if (b[i].szomszed->field == FieldEl || aktSim.mezo == FieldEl) // a szomszéd cella termikusan nem szimulált, szigetelõ!
                                aktCella.set_face_th(b[i].oldal, true, 0, 0, ~0, b[i].y, b[i].c);
                            else aktCella.set_face_th(b[i].oldal, false, 0, 1, ~0, b[i].y, b[i].c); // igazi szomszédja van
                        }
                        else aktCella.set_face_th(b[i].oldal, false, 0, 0, ~0, b[i].y, b[i].c); // nem szimuláljuk termikusan
                    }
                }

                // high res region, 1. kör

                uns x_mul = 0, y_mul = 0, z_mul = 0;
                for (uns i = 0; i < aktMod.righ_res_regions.size(); i++) {
                    if (x >= aktMod.righ_res_regions[i].x1 && x <= aktMod.righ_res_regions[i].x2
                        && y >= aktMod.righ_res_regions[i].y1 && y <= aktMod.righ_res_regions[i].y2
                        && z >= aktMod.righ_res_regions[i].z1 && z <= aktMod.righ_res_regions[i].z2) {
                        x_mul = aktMod.righ_res_regions[i].x_res; // ha átfednek, akkor az utolsó marad, ahogy kell
                        y_mul = aktMod.righ_res_regions[i].y_res;
                        z_mul = aktMod.righ_res_regions[i].z_res;
                    }
                }

                aktCella.set_belso_cellak_mul(x_mul, y_mul, z_mul);
            }

    // 2. kör: szomszédos cellák alapján csomópontszám illesztése
    //         egyszerre állítja az oldalpárt
    // beállítja a csatlakozó face indexeket is: két csatlakozó face indexe két egymás utáni érték: egy páratlan szám és az utána következõ páros szám
    // a centroid is ugyanolyan indexet kap, mint a sima csatlakozó

    uns akt_csatlakozo_index = 1;
    for (uns z = 0; z < z_res; z++)
        for (uns y = 0; y < y_res; y++)
            for (uns x = 0; x < x_res; x++) {
                t_modell_cella & aktCella = modell_racs.getref(x, y, z);

                if (!aktCella.is_cella)
                    continue;

                cuns akt_bx = aktCella.belso_cellak.x_size();
                cuns akt_by = aktCella.belso_cellak.y_size();
                cuns akt_bz = aktCella.belso_cellak.z_size();
                cuns akt_s = aktCella.belso_cellak.size();

                // belsõ cellák belsõ face-einek indexelése

                if (akt_s > 0) {
                    for (uns bz = 0; bz < akt_bz; bz++)
                        for (uns by = 0; by < akt_by; by++)
                            for (uns bx = 0; bx < akt_bx; bx++) {
                                t_modell_cella & akt_bc = aktCella.belso_cellak.getref(bx, by, bz);
                                if (bx < akt_bx - 1) // van belsõ szomszédja x irányban
                                    face_par_index_beallito(akt_bc.face_adat[EAST], aktCella.belso_cellak.getref(bx + 1, by, bz).face_adat[WEST], akt_csatlakozo_index);
                                if (by < akt_by - 1) // van belsõ szomszédja y irányban
                                    face_par_index_beallito(akt_bc.face_adat[NORTH], aktCella.belso_cellak.getref(bx, by + 1, bz).face_adat[SOUTH], akt_csatlakozo_index);
                                if (bz < akt_bz - 1) // van belsõ szomszédja z irányban
                                    face_par_index_beallito(akt_bc.face_adat[TOP], aktCella.belso_cellak.getref(bx, by, bz + 1).face_adat[BOTTOM], akt_csatlakozo_index);
                            }
                }

                // akt.EAST + szomszéd.WEST kiegyenlítése

                if (x != x_res - 1 && modell_racs.getref(x + 1, y, z).is_cella) { // Nem a szélén van és mellette nem belsõ perem van, hanem cella.
                    t_modell_cella & szomszedCella = modell_racs.getref(x + 1, y, z);
                    cuns sz_bx = szomszedCella.belso_cellak.x_size();
                    cuns sz_by = szomszedCella.belso_cellak.y_size();
                    cuns sz_bz = szomszedCella.belso_cellak.z_size();
                    cuns sz_s = szomszedCella.belso_cellak.size();
                    if (akt_s != 0 || sz_s != 0) { // Ha legalább az egyik High.res.region
                        if (akt_s == 0) { // A szomszéd High.res, az akt nem.
                            aktCella.face_adat[EAST].face2face(sz_by, sz_bz, aktCella.is_el, aktCella.is_th);
                        }
                        else if (sz_s == 0) { // Az akt High.res, a szomszéd nem
                            szomszedCella.face_adat[WEST].face2face(akt_by, akt_bz, szomszedCella.is_el, szomszedCella.is_th);
                        }
                        else if (akt_by != sz_by || akt_bz != sz_bz) {
                            // mindkettõ high.res, de eltérõ a felbontásuk
                            if (akt_by % sz_by != 0 && sz_by % akt_by != 0)
                                throw hiba("apa::build_modell_racs", "not compatible high res regions in y direction: (%u,%u,%u EAST): %u vs %u", x, y, z, akt_by, sz_by);
                            if (akt_bz % sz_bz != 0 && sz_bz % akt_bz != 0)
                                throw hiba("apa::build_modell_racs", "not compatible high res regions in z direction: (%u,%u,%u EAST): %u vs %u", x, y, z, akt_bz, sz_bz);
                            cuns max_by = akt_by > sz_by ? akt_by : sz_by;
                            cuns max_bz = akt_bz > sz_bz ? akt_bz : sz_bz;
                            if (akt_by != max_by || akt_bz != max_bz) {
                                cuns dy = max_by / akt_by;
                                cuns dz = max_bz / akt_bz;
                                for (uns j = 0; j < akt_bz; j++)
                                    for (uns i = 0; i < akt_by; i++)
                                        aktCella.belso_cellak.getref(akt_bx - 1, i, j).face_adat[EAST].face2face(dy, dz, aktCella.is_el, aktCella.is_th);
                                aktCella.face_adat[EAST].kulso_el_db *= dy*dz;
                                aktCella.face_adat[EAST].kulso_th_db *= dy*dz;
                            }
                            if (sz_by != max_by || sz_bz != max_bz) {
                                cuns dy = max_by / sz_by;
                                cuns dz = max_bz / sz_bz;
                                for (uns j = 0; j < sz_bz; j++)
                                    for (uns i = 0; i < sz_by; i++)
                                        szomszedCella.belso_cellak.getref(0, i, j).face_adat[WEST].face2face(dy, dz, szomszedCella.is_el, szomszedCella.is_th);
                                szomszedCella.face_adat[WEST].kulso_el_db *= dy*dz;
                                szomszedCella.face_adat[WEST].kulso_th_db *= dy*dz;
                            }
                        } // else mindkettõ egyforma felbontású high.res => nothing to do
                    }
                    face_index_beallito(aktCella, szomszedCella, akt_csatlakozo_index, X_IRANY);
                }

                // akt.NORTH + szomszéd.SOUTH kiegyenlítése

                if (y != y_res - 1 && modell_racs.getref(x, y + 1, z).is_cella) { // Nem a szélén van és mellette nem belsõ perem van, hanem cella.
                    t_modell_cella & szomszedCella = modell_racs.getref(x, y + 1, z);
                    cuns sz_bx = szomszedCella.belso_cellak.x_size();
                    cuns sz_by = szomszedCella.belso_cellak.y_size();
                    cuns sz_bz = szomszedCella.belso_cellak.z_size();
                    cuns sz_s = szomszedCella.belso_cellak.size();
                    if (akt_s != 0 || sz_s != 0) { // Ha legalább az egyik High.res.region
                        if (akt_s == 0) { // A szomszéd High.res, az akt nem.
                            aktCella.face_adat[NORTH].face2face(sz_bx, sz_bz, aktCella.is_el, aktCella.is_th);
                        }
                        else if (sz_s == 0) { // Az akt High.res, a szomszéd nem
                            szomszedCella.face_adat[SOUTH].face2face(akt_bx, akt_bz, szomszedCella.is_el, szomszedCella.is_th);
                        }
                        else if (akt_bx != sz_bx || akt_bz != sz_bz) {
                            // mindkettõ high.res, de eltérõ a felbontásuk
                            if (akt_bx % sz_bx != 0 && sz_bx % akt_bx != 0)
                                throw hiba("apa::build_modell_racs", "not compatible high res regions in x direction: (%u,%u,%u NORTH): %u vs %u", x, y, z, akt_bx, sz_bx);
                            if (akt_bz % sz_bz != 0 && sz_bz % akt_bz != 0)
                                throw hiba("apa::build_modell_racs", "not compatible high res regions in z direction: (%u,%u,%u NORTH): %u vs %u", x, y, z, akt_bz, sz_bz);
                            cuns max_bx = akt_bx > sz_bx ? akt_bx : sz_bx;
                            cuns max_bz = akt_bz > sz_bz ? akt_bz : sz_bz;
                            if (akt_bx != max_bx || akt_bz != max_bz) {
                                cuns dx = max_bx / akt_bx;
                                cuns dz = max_bz / akt_bz;
                                for (uns j = 0; j < akt_bz; j++)
                                    for (uns i = 0; i < akt_bx; i++)
                                        aktCella.belso_cellak.getref(i, akt_by - 1, j).face_adat[NORTH].face2face(dx, dz, aktCella.is_el, aktCella.is_th);
                                aktCella.face_adat[NORTH].kulso_el_db *= dx*dz;
                                aktCella.face_adat[NORTH].kulso_th_db *= dx*dz;
                            }
                            if (sz_bx != max_bx || sz_bz != max_bz) {
                                cuns dx = max_bx / sz_bx;
                                cuns dz = max_bz / sz_bz;
                                for (uns j = 0; j < sz_bz; j++)
                                    for (uns i = 0; i < sz_bx; i++)
                                        szomszedCella.belso_cellak.getref(i, 0, j).face_adat[SOUTH].face2face(dx, dz, szomszedCella.is_el, szomszedCella.is_th);
                                szomszedCella.face_adat[SOUTH].kulso_el_db *= dx*dz;
                                szomszedCella.face_adat[SOUTH].kulso_th_db *= dx*dz;
                            }
                        }// else mindkettõ egyforma felbontású high.res => nothing to do
                    }
                    face_index_beallito(aktCella, szomszedCella, akt_csatlakozo_index, Y_IRANY);
                }

                // akt.TOP + szomszéd.BOTTOM kiegyenlítése

                if (z != z_res - 1 && modell_racs.getref(x, y, z + 1).is_cella) { // Nem a szélén van és mellette nem belsõ perem van, hanem cella.
                    t_modell_cella & szomszedCella = modell_racs.getref(x, y, z + 1);
                    cuns sz_bx = szomszedCella.belso_cellak.x_size();
                    cuns sz_by = szomszedCella.belso_cellak.y_size();
                    cuns sz_bz = szomszedCella.belso_cellak.z_size();
                    cuns sz_s = szomszedCella.belso_cellak.size();
                    if (akt_s != 0 || sz_s != 0) { // Ha legalább az egyik High.res.region
                        if (akt_s == 0) { // A szomszéd High.res, az akt nem.
                            aktCella.face_adat[TOP].face2face(sz_bx, sz_by, aktCella.is_el, aktCella.is_th);
                        }
                        else if (sz_s == 0) { // Az akt High.res, a szomszéd nem
                            szomszedCella.face_adat[BOTTOM].face2face(akt_bx, akt_by, szomszedCella.is_el, szomszedCella.is_th);
                        }
                        else if (akt_bx != sz_bx || akt_by != sz_by) {
                            // mindkettõ high.res, de eltérõ a felbontásuk
                            if (akt_bx % sz_bx != 0 && sz_bx % akt_bx != 0)
                                throw hiba("apa::build_modell_racs", "not compatible high res regions in x direction: (%u,%u,%u TOP): %u vs %u", x, y, z, akt_bx, sz_bx);
                            if (akt_by % sz_by != 0 && sz_by % akt_by != 0)
                                throw hiba("apa::build_modell_racs", "not compatible high res regions in z direction: (%u,%u,%u TOP): %u vs %u", x, y, z, akt_by, sz_by);
                            cuns max_bx = akt_bx > sz_bx ? akt_bx : sz_bx;
                            cuns max_by = akt_by > sz_by ? akt_by : sz_by;
                            if (akt_bx != max_bx || akt_by != max_by) {
                                cuns dx = max_bx / akt_bx;
                                cuns dy = max_by / akt_by;
                                for (uns j = 0; j < akt_by; j++)
                                    for (uns i = 0; i < akt_bx; i++)
                                        aktCella.belso_cellak.getref(i, j, akt_bz - 1).face_adat[TOP].face2face(dx, dy, aktCella.is_el, aktCella.is_th);
                                aktCella.face_adat[TOP].kulso_el_db *= dx*dy;
                                aktCella.face_adat[TOP].kulso_th_db *= dx*dy;
                            }
                            if (sz_bx != max_bx || sz_by != max_by) {
                                cuns dx = max_bx / sz_bx;
                                cuns dy = max_by / sz_by;
                                for (uns j = 0; j < sz_by; j++)
                                    for (uns i = 0; i < sz_bx; i++)
                                        szomszedCella.belso_cellak.getref(i, j, 0).face_adat[BOTTOM].face2face(dx, dy, szomszedCella.is_el, szomszedCella.is_th);
                                szomszedCella.face_adat[BOTTOM].kulso_el_db *= dx*dy;
                                szomszedCella.face_adat[BOTTOM].kulso_th_db *= dx*dy;
                            }
                        }// else mindkettõ egyforma felbontású high.res => nothing to do
                    }
                    face_index_beallito(aktCella, szomszedCella, akt_csatlakozo_index, Z_IRANY);
                }
            }
    uns startindex = 1;
    for (uns z = 0; z < z_res; z++)
        for (uns y = 0; y < y_res; y++)
            for (uns x = 0; x < x_res; x++) {
                if (x == 35 && y == 30 && z == 5)
                    printf("\nCellaindex=%u (%u)\n", startindex, z);
                modell_racs.getref(x, y, z).set_cella_index(startindex);
            }

    printf("%u cells are indexed\n", startindex - 1);
    csatlakozo_face_db = akt_csatlakozo_index - 1;
    return startindex - 1;
}


//***********************************************************************
void sugarfeldolgozo::kulso_oldalak_gyujtese(const simulation & aktSim, const tomb3d<t_modell_cella>& modell_racs) {
//***********************************************************************
    model & aktMod = *aktSim.pmodel;
    cuns x_res = aktMod.x_res;
    cuns y_res = aktMod.y_res;
    cuns z_res = aktMod.z_res;
    //**********************************************
    // Tükrözõdõ középvonalhoz (half_reflex) ezt true-ra állítjuk
    bool is_not_half = !aktMod.is_half;
    //bool is_not_half = true;
    //**********************************************
    kulso_oldal uj_oldal;
    bool hatar;
    kulso_oldalak.clear();
    for (uns z = 0; z < z_res; z++)
        for (uns y = 0; y < y_res; y++)
            for (uns x = 0; x < x_res; x++) {
                const t_modell_cella & aktCella = modell_racs.getconstref(x, y, z);
                if (aktCella.is_cella && aktCella.pmat->is_fenypor) {
                    
                    // bottom
                    
                    if (dir_mode > 3 || !is_top_junction) {
                        hatar = false;
                        if (z == 0)hatar = true;
                        else {
                            if (!modell_racs.getconstref(x, y, z - 1).is_cella)hatar = true;
                            else if (!modell_racs.getconstref(x, y, z - 1).pmat->is_fenypor)hatar = true;
                        }
                        if (hatar) {
                            uj_oldal.set(aktCella.sarkok.x0, aktCella.sarkok.y0, aktCella.sarkok.z0, aktCella.sarkok.x1, aktCella.sarkok.y1, aktCella.sarkok.z0, aktCella, x, y, z, false);
                            kulso_oldalak.add(uj_oldal);
                        }
                    }

                    if (dir_mode > 2) {
                        // west
                        hatar = false;
                        if (x == 0)hatar = true;
                        else {
                            if (!modell_racs.getconstref(x - 1, y, z).is_cella)hatar = true;
                            else if (!modell_racs.getconstref(x - 1, y, z).pmat->is_fenypor)hatar = true;
                        }
                        if (hatar) {
                            if (is_not_half || x != 0) {
                                uj_oldal.set(aktCella.sarkok.x0, aktCella.sarkok.y0, aktCella.sarkok.z0, aktCella.sarkok.x0, aktCella.sarkok.y1, aktCella.sarkok.z1, aktCella, x, y, z, true);
                                kulso_oldalak.add(uj_oldal);
                            }
                        }
                        // east
                        hatar = false;
                        if (x == x_res - 1)hatar = true;
                        else {
                            if (!modell_racs.getconstref(x + 1, y, z).is_cella)hatar = true;
                            else if (!modell_racs.getconstref(x + 1, y, z).pmat->is_fenypor)hatar = true;
                        }
                        if (hatar) {
                            uj_oldal.set(aktCella.sarkok.x1, aktCella.sarkok.y0, aktCella.sarkok.z0, aktCella.sarkok.x1, aktCella.sarkok.y1, aktCella.sarkok.z1, aktCella, x, y, z, true);
                            kulso_oldalak.add(uj_oldal);
                        }
                        // south
                        hatar = false;
                        if (y == 0)hatar = true;
                        else {
                            if (!modell_racs.getconstref(x, y - 1, z).is_cella)hatar = true;
                            else if (!modell_racs.getconstref(x, y - 1, z).pmat->is_fenypor)hatar = true;
                        }
                        if (hatar) {
                            uj_oldal.set(aktCella.sarkok.x0, aktCella.sarkok.y0, aktCella.sarkok.z0, aktCella.sarkok.x1, aktCella.sarkok.y0, aktCella.sarkok.z1, aktCella, x, y, z, true);
                            kulso_oldalak.add(uj_oldal);
                        }
                        // north
                        hatar = false;
                        if (y == y_res - 1)hatar = true;
                        else {
                            if (!modell_racs.getconstref(x, y + 1, z).is_cella)hatar = true;
                            else if (!modell_racs.getconstref(x, y + 1, z).pmat->is_fenypor)hatar = true;
                        }
                        if (hatar) {
                            uj_oldal.set(aktCella.sarkok.x0, aktCella.sarkok.y1, aktCella.sarkok.z0, aktCella.sarkok.x1, aktCella.sarkok.y1, aktCella.sarkok.z1, aktCella, x, y, z, true);
                            kulso_oldalak.add(uj_oldal);
                        }
                    }

                    // top

                    if (dir_mode > 3 || is_top_junction) {
                        hatar = false;
                        if (z == z_res - 1)hatar = true;
                        else {
                            const t_modell_cella & mc = modell_racs.getconstref(x, y, z + 1);
                            if (!mc.is_cella)hatar = true;
                            else if (!mc.pmat->is_fenypor)hatar = true;
                        }
                        if (hatar) {
                            uj_oldal.set(aktCella.sarkok.x0, aktCella.sarkok.y0, aktCella.sarkok.z1, aktCella.sarkok.x1, aktCella.sarkok.y1, aktCella.sarkok.z1, aktCella, x, y, z, false);
                            kulso_oldalak.add(uj_oldal);
                        }
                    }
                }
            }
}


//***********************************************************************
inline bool is_in_the_cell(const Vec3d & pont, const t_modell_cella & cella) {
//***********************************************************************
    dbl x0, x1, y0, y1, z0, z1;
    if (cella.sarkok.x0 <= cella.sarkok.x1) { x0 = cella.sarkok.x0; x1 = cella.sarkok.x1; }
    else { x0 = cella.sarkok.x1; x1 = cella.sarkok.x0; }
    if (cella.sarkok.y0 <= cella.sarkok.y1) { y0 = cella.sarkok.y0; y1 = cella.sarkok.y1; }
    else { y0 = cella.sarkok.y1; y1 = cella.sarkok.y0; }
    if (cella.sarkok.z0 <= cella.sarkok.z1) { z0 = cella.sarkok.z0; z1 = cella.sarkok.z1; }
    else { z0 = cella.sarkok.z1; z1 = cella.sarkok.z0; }
    return pont.x >= x0 && pont.x <= x1 && pont.y >= y0 && pont.y <= y1 && pont.z >= z0 && pont.z <= z1;
}


//***********************************************************************
int sugar_cmp(const void *p1, const void *p2) {
//***********************************************************************
    const sugar_adat ** a1 = (const sugar_adat **)p1;
    const sugar_adat ** a2 = (const sugar_adat **)p2;
    if ((*a1)->P < (*a2)->P) return -1;
    if ((*a1)->P > (*a2)->P) return  1;
    return 0;
}


//***********************************************************************
void _unsprintin(uns szam, char *s, bool is_comma) {
//***********************************************************************
    if (szam < 1000)
        sprintf(s, "%u", szam);
    else {
        _unsprintin(szam / 1000, s, is_comma);
        sprintf(s + strlen(s), "%c%03u", is_comma?',':' ', szam % 1000);
    }
}


//***********************************************************************
char * unsprint(uns szam, bool is_comma = false) {
//***********************************************************************
    static char s[10][100];
    static uns index = 0;
    _unsprintin(szam, s[index % 10], is_comma);
    return s[(index++) % 10];
}


//***********************************************************************
dbl sugarfeldolgozo::kek_sugarak_gyujtese(const simulation & aktSim, const tomb3d<t_modell_cella> & modell_racs) {
// végigmegy az összes junctionon, kiszámítja az összes kék sugarat és beteszi a listákba, kidobja a legkisebb teljesítményûeket
// vissza: súlyozott kék sugár hossz
//***********************************************************************
    model & aktMod = *aktSim.pmodel;
    cuns x_res = aktMod.x_res;
    cuns y_res = aktMod.y_res;
    cuns z_res = aktMod.z_res;
    uns sum_db_betett = 0;
    cuns x_min = 0, x_max = x_res;
    cuns y_min = 0, y_max = y_res;
    cuns z_min = 0, z_max = z_res;
    //////////////////////////////////////////////////////////////////////
    bool is_paros = (ray_per_cell_dir % 2) == 0;
    cuns ray_per_dir_per_cell = dir_mode == 1 ? 1 : (is_paros ? (ray_per_cell_dir + 1) : ray_per_cell_dir); // ray_per_dir_per_cell*ray_per_dir_per_cell db sugár indul a cellából, csak páratlan lehet
    //////////////////////////////////////////////////////////////////////
    cuns ray_per_cell = is_paros ? (ray_per_dir_per_cell*ray_per_dir_per_cell / 2 + 1) : ray_per_dir_per_cell*ray_per_dir_per_cell;
    printf("gather blue rays\n");
    printf("%u ray / cell\n", ray_per_cell);
    uns junctdb = 0, lpdb = 0;
    for (uns z = z_min; z < z_max; z++)
        for (uns y = y_min; y < y_max; y++)
            for (uns x = x_min; x < x_max; x++) {
                const t_modell_cella & aktCella = modell_racs.getconstref(x, y, z);
                if (aktCella.is_cella && aktCella.is_junction)
                    junctdb++;
                if (aktCella.is_cella && aktCella.pmat->is_fenypor)
                    lpdb++;
                if (aktCella.is_cella && aktCella.is_junction // ha van fõ sugárzó irányû junction (top/bottom)
                    && ((is_top_junction && aktCella.junction_top_face != 0) || (!is_top_junction && aktCella.junction_bottom_face != 0))) {
                    sugar_adat akt_sugar;
                    // A kiindulópontot kicsit eltoljuk, hogy ne pont élben metssze az átlós sugár a cellákat
                    Vec3d p1{ 0.5*(aktCella.sarkok.x0 + aktCella.sarkok.x1 + 1e-10), 0.5*(aktCella.sarkok.y0 + aktCella.sarkok.y1 + 3e-10), 1e-10 + is_top_junction ? aktCella.sarkok.z1 : aktCella.sarkok.z0 };
                    Vec3d dp1i = { (aktCella.sarkok.x1 - aktCella.sarkok.x0) / ray_per_dir_per_cell, 0, 0 };
                    Vec3d dp1j = { 0, (aktCella.sarkok.y1 - aktCella.sarkok.y0) / ray_per_dir_per_cell, 0 };
                    akt_sugar.x0 = x;
                    akt_sugar.y0 = y;
                    akt_sugar.z0 = z;
                    akt_sugar.e.p1 = p1;
                    akt_sugar.src_cella_index = aktCella.cella_index;
                    akt_sugar.face_index = is_top_junction ? aktCella.junction_top_face : aktCella.junction_bottom_face;

                    kek_sugarlista.reallocate(kek_sugarlista.size() + 1);
                    sugar_lista & akt_sugar_lista = kek_sugarlista.getLast();
                    sugar_cella & uj_sugar_cella = akt_sugar_lista.kiindulo_cella;
                    uj_sugar_cella.cella_index = akt_sugar.src_cella_index;
                    uj_sugar_cella.face_index = akt_sugar.face_index;

                    dbl sum_P = 0; // egy cellából kilõtt sugarak összteljesítménye
                    uns db_in = 0;
                    dbl A_junction = 1; //fabs((aktCella.sarkok.x1 - aktCella.sarkok.x0)*(aktCella.sarkok.y1 - aktCella.sarkok.y0));
                    // egyformának tekintjük a junction cellák területét, mert ez a számítás úgysem a ténylegest adja
                    akt_sugar.A_kek = A_junction;
                    for (uns i = 0; i < kulso_oldalak.size(); i++) {
                        const kulso_oldal & cel_oldal = kulso_oldalak[i];
                        if (dir_mode > 1 || (cel_oldal.cella_x == x && cel_oldal.cella_y == y)) {
                            akt_sugar.x1 = cel_oldal.cella_x;
                            akt_sugar.y1 = cel_oldal.cella_y;
                            akt_sugar.z1 = cel_oldal.cella_z;
                            akt_sugar.e.p2 = cel_oldal.c;
                            Vec3d dp2i = (cel_oldal.p2 - cel_oldal.p1) * (1.0 / ray_per_dir_per_cell);
                            akt_sugar.is_oldalso = cel_oldal.is_oldalso;
                            Vec3d dp2j = (cel_oldal.p4 - cel_oldal.p1) * (1.0 / ray_per_dir_per_cell);
                            bool is_nem_vissza = true;
                            if (is_top_junction)is_nem_vissza = akt_sugar.e.p2.z > akt_sugar.e.p1.z;
                            else is_nem_vissza = akt_sugar.e.p2.z < akt_sugar.e.p1.z;
                            if (is_nem_vissza && is_in_the_cell(akt_sugar.e.get_p2_kozelpont(), *cel_oldal.p_modell_cella)) {
                                dbl P_akt = dir_mode==1 ? 1.0 : (terszog(cel_oldal.p1 - p1, cel_oldal.p2 - p1, cel_oldal.p3 - p1)
                                    + terszog(cel_oldal.p3 - p1, cel_oldal.p4 - p1, cel_oldal.p1 - p1));
                                akt_sugar.P_kek_sarga = P_akt;
                                akt_sugar.P = A_junction*P_akt;
                                sum_P += P_akt;
                                akt_sugar_lista.add(akt_sugar, ray_per_dir_per_cell, is_paros, dp1i, dp1j, dp2i, dp2j);
                                db_in += ray_per_cell;
                            }
                        }
                    }
                    sum_db_betett += db_in;
                    global_kek.add(akt_sugar_lista);
                    printf(".");
                    // A sugár belép vagy kilép a cél cellából? => get_p2_kozelpont benne van-e a cél cellában?
                    // Az emiatt eldobott teljesítményt összegyûjtjük, és a végén felszorozzuk az egyes P_kek_sarga-kat
                    // =>nem! a sugarak által képviselt teljesítményt kell összeszámolni, és ezzel felszorozni a teljesítményt,
                    // mert nem biztos, hogy minden térszögben van fénypor
                }
            }
    printf("\njunction=%u, lp=%u, ko=%u\n", junctdb, lpdb, kulso_oldalak.size());
    fprintf(logfile, "junction=%u\nphosphor_cell=%u\nexternal_side=%u\n", junctdb, lpdb, kulso_oldalak.size());
    printf("\nblue sort started\r");
    uns i = 0;
    if (dir_mode > 1) {
        if (global_kek.t.size() > 1) qsort(&global_kek.t[0], global_kek.t.size(), sizeof(sugar_adat*), sugar_cmp);
        printf("blue remove       \r");
        dbl sum_P = 0; // egy cellából kilõtt sugarak összteljesítménye
        for (uns i = 0; i < global_kek.t.size(); i++)
            sum_P += global_kek.t[i]->P;
        dbl P_kidobhato = sum_P * cut_level;
        dbl P_kidobott = 0; // sum_P-t a P_kek_sarga alapján számoljuk
        for (i = 0; P_kidobott < P_kidobhato && i < global_kek.t.size(); i++) {
            P_kidobott += global_kek.t[i]->P;
            global_kek.t[i]->remove();
            delete global_kek.t[i];
        }
    }
    for (uns i = 0; i < kek_sugarlista.size(); i++) {
        kek_sugarlista[i].recount();
        kek_sugarlista[i].normalize_P_kek_sarga();
    }
    printf("blue remove finished\n");
    dbl d_weighted = get_blue_weighted_length_and_recalculate_blue_P();
    printf("all blue / after remove (removed) = %s / %s (%s)\nd_blue_weighted = %g\n",
        unsprint(sum_db_betett), unsprint(global_kek.t.size() - i), unsprint(i), d_weighted);
    fprintf(logfile, "all blue=%s\nafter remove=%s\n(removed=%s)\n",
        unsprint(sum_db_betett,true), unsprint(global_kek.t.size() - i,true), unsprint(i,true));
    return d_weighted;
}


//***********************************************************************
dbl sugarfeldolgozo::sarga_sugarak_gyujtese(const simulation & aktSim, const tomb3d<t_modell_cella> & modell_racs) {
// végigmegy az összes fénypor cellán, kiszámítja az összes sárga sugarat 
// és beteszi a listákba, kidobja a legkisebb teljesítményûeket
// vissza: súlyozott sárga sugár hossz
//***********************************************************************
    model & aktMod = *aktSim.pmodel;
    cuns x_res = aktMod.x_res;
    cuns y_res = aktMod.y_res;
    cuns z_res = aktMod.z_res;
    uns sum_db_betett = 0;
    cuns x_min = 0, x_max = x_res;
    cuns y_min = 0, y_max = y_res;
    cuns z_min = 0, z_max = z_res;
    uns dbdb = 0;
    // Az átlagos cellánkénti kék teljesítmény meghatározása
    dbl sum_sum_P_kek = 0;
    uns sum_sum_P_kek_db = 0;
    for (uns z = z_min; z < z_max; z++)
        for (uns y = y_min; y < y_max; y++)
            for (uns x = x_min; x < x_max; x++) {
                const t_modell_cella & aktCella = modell_racs.getconstref(x, y, z);
                if (aktCella.is_cella && aktCella.pmat->is_fenypor) {
                    sum_sum_P_kek += aktCella.sum_P_kek;
                    sum_sum_P_kek_db++;
                }
            }
    dbl avg_sum_P_kek = sum_sum_P_kek_db>0 ? sum_sum_P_kek / sum_sum_P_kek_db : 1;
    printf("gather yellow rays\n");
    for (uns z = z_min; z < z_max; z++)
        for (uns y = y_min; y < y_max; y++)
            for (uns x = x_min; x < x_max; x++) {
                const t_modell_cella & aktCella = modell_racs.getconstref(x, y, z);
                if (aktCella.is_cella && aktCella.pmat->is_fenypor/* && aktCella.sum_P_kek > 0*/) {
                    sugar_adat akt_sugar;
                    // A kiindulópontot kicsit eltoljuk, hogy ne pont élben metssze az átlós sugár a cellákat
                    Vec3d p1{ 0.5*(aktCella.sarkok.x0 + aktCella.sarkok.x1 + 1e-10), 0.5*(aktCella.sarkok.y0 + aktCella.sarkok.y1 + 3e-10), 0.5*(aktCella.sarkok.z0 + aktCella.sarkok.z1 + 2e-10) };
                    akt_sugar.x0 = x;
                    akt_sugar.y0 = y;
                    akt_sugar.z0 = z;
                    akt_sugar.e.p1 = p1;
                    akt_sugar.src_cella_index = aktCella.cella_index;
                    akt_sugar.face_index = 0;

                    sarga_sugarlista.reallocate(sarga_sugarlista.size() + 1);
                    sugar_lista & akt_sugar_lista = sarga_sugarlista.getLast();
                    sugar_cella & uj_sugar_cella = akt_sugar_lista.kiindulo_cella;
                    uj_sugar_cella.cella_index = akt_sugar.src_cella_index;
                    uj_sugar_cella.face_index = 0;

                    dbl sum_P = 0; // egy cellából kilõtt sugarak összteljesítménye
                    uns db_in = 0;
                    akt_sugar.A_kek = 1;//aktCella.sum_P_kek;
                    //printf("%g, %g\n", aktCella.sum_P_kek, avg_sum_P_kek);
                    //akt_sugar.A_kek = A_junction;
                    for (uns i = 0; i < kulso_oldalak.size(); i++) {
                        const kulso_oldal & cel_oldal = kulso_oldalak[i];
                        if (dir_mode > 1 || (cel_oldal.cella_x == x && cel_oldal.cella_y == y)) {
                            akt_sugar.x1 = cel_oldal.cella_x;
                            akt_sugar.y1 = cel_oldal.cella_y;
                            akt_sugar.z1 = cel_oldal.cella_z;
                            akt_sugar.e.p2 = cel_oldal.c;
                            akt_sugar.is_oldalso = cel_oldal.is_oldalso;
                            if (is_in_the_cell(akt_sugar.e.get_p2_kozelpont(), *cel_oldal.p_modell_cella)) {
                                dbl P_akt = dir_mode==1 ? 1.0 : (terszog(cel_oldal.p1 - p1, cel_oldal.p2 - p1, cel_oldal.p3 - p1)
                                    + terszog(cel_oldal.p3 - p1, cel_oldal.p4 - p1, cel_oldal.p1 - p1));
                                akt_sugar.P_kek_sarga = P_akt;
                                akt_sugar.P = (aktCella.sum_P_kek + avg_sum_P_kek) * P_akt;
                                //akt_sugar.P = P_akt;
                                //akt_sugar.P = A_junction*P_akt;
                                sum_P += P_akt;
                                akt_sugar_lista.add(akt_sugar);
                                db_in++;
                            }
                        }
                    }
                    sum_db_betett += db_in;
                    global_sarga.add(akt_sugar_lista);
                    if (dbdb++ % 10 == 0) printf(".");
                }
            }
    printf("\nyellow sort started\r");
    uns i = 0;
    if (dir_mode > 1) {
        if (global_sarga.t.size() > 1) qsort(&global_sarga.t[0], global_sarga.t.size(), sizeof(sugar_adat*), sugar_cmp);
        printf("yellow remove       \r");
        dbl sum_P = 0; // egy cellából kilõtt sugarak összteljesítménye
        for (uns i = 0; i < global_sarga.t.size(); i++)
            sum_P += global_sarga.t[i]->P;
        dbl P_kidobhato = sum_P * cut_level;
        dbl P_kidobott = 0; // sum_P-t a P_kek_sarga alapján számoljuk
        for (i = 0; P_kidobott < P_kidobhato && i < global_sarga.t.size(); i++) {
            P_kidobott += global_sarga.t[i]->P;
            global_sarga.t[i]->remove();
            delete global_sarga.t[i];
        }
    }
    for (uns i = 0; i < sarga_sugarlista.size(); i++) {
        sarga_sugarlista[i].recount();
        sarga_sugarlista[i].normalize_P_kek_sarga();
    }
    printf("yellow remove finished\n");
    dbl d_weighted = get_yellow_weighted_length_and_recalculate_yellow_P();
    printf("all yellow / after remove (removed) = %s / %s (%s)\nd_yellow_weighted = %g\n",
        unsprint(sum_db_betett), unsprint(global_sarga.t.size() - i), unsprint(i), d_weighted);
    fprintf(logfile, "all yellow=%s\nafter remove=%s\n(removed=%s)\n",
        unsprint(sum_db_betett,true), unsprint(global_sarga.t.size() - i,true), unsprint(i,true));
    return d_weighted;
}


//***********************************************************************
bool metszesszamolo(const egyenes & e, const ketpont & c, dbl & d, char & dir1, char & dir2, Vec3d & mp1, Vec3d & mp2, uns & itt_irany, uns kivant_irany) {
// kivant_irany: 1: x, 2: y, 3: z
//***********************************************************************
    d = 0;
    bool is_elso = true, is_ketto = false;
    itt_irany = 0;
    uns db = 0;
    if (e.get_x_metszet(c.x0, c.y0, c.y1, c.z0, c.z1, mp1)) {
        is_elso = false;
        dir1 = 'W';
        itt_irany |= 1;
        db++;
    }
    if (e.get_x_metszet(c.x1, c.y0, c.y1, c.z0, c.z1, is_elso ? mp1 : mp2)) {
        is_ketto = is_elso ? false : true;
        if (is_elso)dir1 = 'E';
        else dir2 = 'E';
        is_elso = false;
        itt_irany |= 2;
        db++;
    }
    if (e.get_y_metszet(c.y0, c.x0, c.x1, c.z0, c.z1, is_elso ? mp1 : mp2)) {
        is_ketto = is_elso ? false : true;
        if (is_elso)dir1 = 'S';
        else dir2 = 'S';
        is_elso = false;
        itt_irany |= 4;
        db++;
    }
    if (e.get_y_metszet(c.y1, c.x0, c.x1, c.z0, c.z1, is_elso ? mp1 : mp2)) {
        is_ketto = is_elso ? false : true;
        if (is_elso)dir1 = 'N';
        else dir2 = 'N';
        is_elso = false;
        itt_irany |= 8;
        db++;
    }
    if (e.get_z_metszet(c.z0, c.x0, c.x1, c.y0, c.y1, is_elso ? mp1 : mp2)) {
        is_ketto = is_elso ? false : true;
        if (is_elso)dir1 = 'B';
        else dir2 = 'B';
        is_elso = false;
        itt_irany |= 16;
        db++;
    }
    if (e.get_z_metszet(c.z1, c.x0, c.x1, c.y0, c.y1, is_elso ? mp1 : mp2)) {
        is_ketto = is_elso ? false : true;
        if (is_elso)dir1 = 'T';
        else dir2 = 'T';
        is_elso = false;
        itt_irany |= 32;
        db++;
    }
    if (db > 2) {
        ketpont kp = c;
        if (kivant_irany != 1) kp.x0 -= 5e-10;
        if (kivant_irany != 2) kp.y0 -= 5e-10;
        if (kivant_irany != 3) kp.z0 -= 5e-10;
        if (kivant_irany != 1) kp.x1 += 5e-10;
        if (kivant_irany != 2) kp.x1 += 5e-10;
        if (kivant_irany != 3) kp.x1 += 5e-10;
        //printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
        if (kivant_irany > 0)
            metszesszamolo(e, kp, d, dir1, dir2, mp1, mp2, itt_irany, 0);
        //printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
    }
    if (!is_ketto)
        return false;
    d = (mp1 - mp2).hossz();
    return true;
}


//***********************************************************************
dbl sugarfeldolgozo::get_cell_list(bool is_blue, bool is_half, tomb<sugar_cella> & uj_tomb, sugar_adat * p_sugar, tomb3d<t_modell_cella> & modell_racs, dbl P0) {
// sugár alapján cellalistát hoz létre
//***********************************************************************
    dbl sugarhossz = p_sugar->e.get_p1_p2_distance();
    dbl min_d = sugarhossz * cut_level
        / ((p_sugar->x1 - p_sugar->x0)*(p_sugar->x1 - p_sugar->x0)
            + (p_sugar->y1 - p_sugar->y0)*(p_sugar->y1 - p_sugar->y0)
            + (p_sugar->z1 - p_sugar->z0)*(p_sugar->z1 - p_sugar->z0)
            );
    dbl sum_d = 0;
    bool is_vertical = p_sugar->e.get_vertical_angle() < 0.1;
    Vec3d mp1, mp2;
    uns itt_irany;

    // ha egyformák, 0

    int dx = 0, dy = 0, dz = 0;
    if (p_sugar->x1 > p_sugar->x0) dx = 1;
    else if (p_sugar->x1 < p_sugar->x0) dx = -1;
    if (p_sugar->y1 > p_sugar->y0) dy = 1;
    else if (p_sugar->y1 < p_sugar->y0) dy = -1;
    if (p_sugar->z1 > p_sugar->z0) dz = 1;
    else if (p_sugar->z1 < p_sugar->z0) dz = -1;

    // Ha ugyanaz a kiinduló és a cél

    if (dx == 0 && dy == 0 && dz == 0) {
        t_modell_cella & akt_cella = modell_racs.getref(p_sugar->x0, p_sugar->y0, p_sugar->z0);
        dbl d = 0;
        char dir1, dir2;
        metszesszamolo(p_sugar->e, akt_cella.sarkok, d, dir1, dir2, mp1, mp2, itt_irany, 0);
        switch (dir1) {
            case 'W': 
            case 'E': 
                if (p_sugar->e.p1.x > p_sugar->e.p2.x) { // nyugatra megy
                    dx = -1;
                    if (dir1 == 'E')dir1 = 'X'; // keletre nem megy
                    if (dir2 == 'E')dir2 = 'X';
                }
                else { // keletre megy
                    dx = 1;
                    if (dir1 == 'W')dir1 = 'X';
                    if (dir2 == 'W')dir2 = 'X';
                }
                break;
            case 'S': 
            case 'N': 
                if (p_sugar->e.p1.y > p_sugar->e.p2.y) { // délre megy
                    dy = -1;
                    if (dir1 == 'N')dir1 = 'X';
                    if (dir2 == 'N')dir2 = 'X';
                }
                else { // északra megy
                    dy = 1; 
                    if (dir1 == 'S')dir1 = 'X';
                    if (dir2 == 'S')dir2 = 'X';
                }
                break;
            case 'B': 
            case 'T': 
                if (p_sugar->e.p1.z > p_sugar->e.p2.z) { // le megy
                    dz = -1;
                    if (dir1 == 'T')dir1 = 'X';
                    if (dir2 == 'T')dir2 = 'X';
                }
                else { // fel megy
                    dz = 1;
                    if (dir1 == 'B')dir1 = 'X';
                    if (dir2 == 'B')dir2 = 'X';
                }
                break;
        }
        if (akt_cella.is_cella && akt_cella.pmat->is_fenypor) {
            uj_tomb.reallocate(uj_tomb.size() + 1);
            sugar_cella & uj = uj_tomb.getLast();
            uj.cella_index = akt_cella.cella_index;
            uj.x = p_sugar->x0;
            uj.y = p_sugar->y0;
            uj.z = p_sugar->z0;
            uj.d = (float)(0.5*d);
            uj.K_P = (float)1;
            uj.dir1 = dir1;
            uj.dir2 = dir2;
            uj.is_vertical = is_vertical;
            uj.is_oldalsohoz = p_sugar->is_oldalso;
            sum_d += 0.5*d;
            if (is_blue) akt_cella.sum_P_kek += P0;
            else akt_cella.sum_P_sarga += P0;
        }
        min_d = 0;
    }
    else { // ha nem ugyanaz, de foszfor cella, akkor a foszfor közepétõl a sugarat hozzáadjuk
        t_modell_cella & akt_cella = modell_racs.getref(p_sugar->x0, p_sugar->y0, p_sugar->z0);
        if (akt_cella.is_cella && akt_cella.pmat->is_fenypor) {
            dbl d = 0;
            char dir1, dir2;
            metszesszamolo(p_sugar->e, akt_cella.sarkok, d, dir1, dir2, mp1, mp2, itt_irany, 0);
            switch (dir1) {
            case 'W':
            case 'E':
                if (p_sugar->e.p1.x > p_sugar->e.p2.x) { // nyugatra megy
                    if (dir1 == 'E')dir1 = 'X'; // keletre nem megy
                    if (dir2 == 'E')dir2 = 'X';
                }
                else { // keletre megy
                    if (dir1 == 'W')dir1 = 'X';
                    if (dir2 == 'W')dir2 = 'X';
                }
                break;
            case 'S':
            case 'N':
                if (p_sugar->e.p1.y > p_sugar->e.p2.y) { // délre megy
                    if (dir1 == 'N')dir1 = 'X';
                    if (dir2 == 'N')dir2 = 'X';
                }
                else { // északra megy
                    if (dir1 == 'S')dir1 = 'X';
                    if (dir2 == 'S')dir2 = 'X';
                }
                break;
            case 'B':
            case 'T':
                if (p_sugar->e.p1.z > p_sugar->e.p2.z) { // le megy
                    if (dir1 == 'T')dir1 = 'X';
                    if (dir2 == 'T')dir2 = 'X';
                }
                else { // fel megy
                    if (dir1 == 'B')dir1 = 'X';
                    if (dir2 == 'B')dir2 = 'X';
                }
                break;
            }
            uj_tomb.reallocate(uj_tomb.size() + 1);
            sugar_cella & uj = uj_tomb.getLast();
            uj.cella_index = akt_cella.cella_index;
            uj.x = p_sugar->x0;
            uj.y = p_sugar->y0;
            uj.z = p_sugar->z0;
            uj.d = (float)(0.5*d);
            uj.K_P = (float)1;
            uj.dir1 = dir1;
            uj.dir2 = dir2;
            uj.is_vertical = is_vertical;
            uj.is_oldalsohoz = p_sugar->is_oldalso;
            sum_d += 0.5*d;
            if (is_blue) akt_cella.sum_P_kek += P0;
            else akt_cella.sum_P_sarga += P0;
        }
    }

    // ellépegetünk a kiinduló cellából a végsõbe

    uns akt_x = p_sugar->x0, akt_y = p_sugar->y0, akt_z = p_sugar->z0;
    while (akt_x != p_sugar->x1 || akt_y != p_sugar->y1 || akt_z != p_sugar->z1) {
        dbl d = 0;
        char dir1, dir2;
        if (dz != 0 && akt_z != p_sugar->z1 && modell_racs.is_exists(akt_x, akt_y, akt_z + dz) && metszesszamolo(p_sugar->e, modell_racs.getconstref(akt_x, akt_y, akt_z + dz).sarkok, d, dir1, dir2, mp1, mp2, itt_irany, 3)) {
            akt_z += dz;
        }
        else if (dy != 0 && akt_y != p_sugar->y1 && modell_racs.is_exists(akt_x, akt_y + dy, akt_z) && metszesszamolo(p_sugar->e, modell_racs.getconstref(akt_x, akt_y + dy, akt_z).sarkok, d, dir1, dir2, mp1, mp2, itt_irany, 2)) {
            akt_y += dy;
        }
        else if (dx != 0 && akt_x != p_sugar->x1 && modell_racs.is_exists(akt_x + dx, akt_y, akt_z) && metszesszamolo(p_sugar->e, modell_racs.getconstref(akt_x + dx, akt_y, akt_z).sarkok, d, dir1, dir2, mp1, mp2, itt_irany, 1)) {
            akt_x += dx;
        }
        else {
            printf("\n\nError: No neighbour is light cell path.\n");
            exit(1);
        }
        t_modell_cella & akt_cella = modell_racs.getref(akt_x, akt_y, akt_z);
        if (d > min_d && akt_cella.is_cella && akt_cella.pmat->is_fenypor) {
            Vec3d kozepe = 0.5*(mp1 + mp2);
            dbl eddig_hossz = (kozepe - p_sugar->e.p1).hossz();
            dbl hosszarany = eddig_hossz / sugarhossz;
            uj_tomb.reallocate(uj_tomb.size() + 1);
            sugar_cella & uj = uj_tomb.getLast();
            uj.cella_index = akt_cella.cella_index;
            uj.x = akt_x;
            uj.y = akt_y;
            uj.z = akt_z;
            uj.d = (float)d;
            uj.K_P = (float)1;
            uj.dir1 = dir1;
            uj.dir2 = dir2;
            uj.is_vertical = is_vertical;
            uj.is_oldalsohoz = p_sugar->is_oldalso;
            sum_d += d;
            //if (is_blue) akt_cella.sum_P_kek += P0;
            if (is_blue) akt_cella.sum_P_kek += P0*(hosszarany + 0.2);
            else akt_cella.sum_P_sarga += P0;
        }
    }

    // visszaverõdés számítása

    egyenes aktualis_egyenes = p_sugar->e;
    bool is_vege = false;
    bool half_reflex = false;
    uns reflexiok_szama = 0;
    cuns max_reflexio = 5;
    dbl akt_P_szorzo = 1;
    dx = dy = dz = 0;
    if (aktualis_egyenes.p1.x < aktualis_egyenes.p2.x) dx =  1;
    if (aktualis_egyenes.p1.x > aktualis_egyenes.p2.x) dx = -1;
    if (aktualis_egyenes.p1.y < aktualis_egyenes.p2.y) dy =  1;
    if (aktualis_egyenes.p1.y > aktualis_egyenes.p2.y) dy = -1;
    if (aktualis_egyenes.p1.z < aktualis_egyenes.p2.z) dz =  1;
    if (aktualis_egyenes.p1.z > aktualis_egyenes.p2.z) dz = -1;

    while (!is_vege) {
        
        // a következõ metszett cella meghatározása
        
        dbl d = 0;
        char dir1, dir2;
        int akt_dx = 0, akt_dy = 0, akt_dz = 0;
        half_reflex = false;
        if (dz != 0 && modell_racs.is_exists(akt_x, akt_y, akt_z + dz) && metszesszamolo(aktualis_egyenes, modell_racs.getconstref(akt_x, akt_y, akt_z + dz).sarkok, d, dir1, dir2, mp1, mp2, itt_irany, 3)) {
            akt_dz = dz;
        }
        else if (dy != 0 && modell_racs.is_exists(akt_x, akt_y + dy, akt_z) && metszesszamolo(aktualis_egyenes, modell_racs.getconstref(akt_x, akt_y + dy, akt_z).sarkok, d, dir1, dir2, mp1, mp2, itt_irany, 2)) {
            akt_dy = dy;
        }
        else if (dx != 0 && modell_racs.is_exists(akt_x + dx, akt_y, akt_z) && metszesszamolo(aktualis_egyenes, modell_racs.getconstref(akt_x + dx, akt_y, akt_z).sarkok, d, dir1, dir2, mp1, mp2, itt_irany, 1)) {
            akt_dx = dx;
        }
        else if (is_half && dx < 0 && akt_x == 0) { // ha a West oldalon tükrözõdik
            half_reflex = true;
            if (!metszesszamolo(aktualis_egyenes, modell_racs.getconstref(akt_x, akt_y, akt_z).sarkok, d, dir1, dir2, mp1, mp2, itt_irany, 1))
                is_vege = true;
            if (dir1 != 'W' && dir2 != 'W')
                is_vege = true;
            else {
                if (dir1 == 'W')
                    dir1 = 'E';
                if (dir2 == 'W')
                    dir2 = 'E';
            }
            akt_dx = dx;
        }
        else {
            is_vege = true;
        }
        
        if (!is_vege) {
            t_modell_cella & kov_cella = half_reflex 
                ? modell_racs.getref(akt_x, akt_y, akt_z) // dummy
                : modell_racs.getref(akt_x + akt_dx, akt_y + akt_dy, akt_z + akt_dz);
            if (!kov_cella.is_cella)
                is_vege = true;
            else {
                dbl prev_K_P = 1;
                if (!half_reflex && kov_cella.pmat->is_fenypor) {
                    // ha az fénypor => továbblépés oda
                    akt_x += akt_dx, akt_y += akt_dy, akt_z += akt_dz;
                }
                else {
                    // ha az nem fénypor => hány visszaverõdés volt > max vagy reflexió=0? ha igen, vége. ha nem, egyenes tükrözése
                    dbl akt_K_P = half_reflex ? 1.0 : kov_cella.pmat->reflectivity;
                    if (reflexiok_szama >= max_reflexio || akt_K_P == 0)
                        is_vege = true;
                    else {
                        prev_K_P = akt_K_P;
                        reflexiok_szama++;
                        uns kiv_ir = 0;
                        if (akt_dx > 0) { // a W oldalon van a metszet
                            dx = -dx;
                            if (dir1 == 'W') {
                                aktualis_egyenes.p1.x = 2 * mp1.x - aktualis_egyenes.p1.x;
                                aktualis_egyenes.p2.x = 2 * mp1.x - aktualis_egyenes.p2.x;
                            }
                            else if (dir2 == 'W') {
                                aktualis_egyenes.p1.x = 2 * mp2.x - aktualis_egyenes.p1.x;
                                aktualis_egyenes.p2.x = 2 * mp2.x - aktualis_egyenes.p2.x;
                            }
                            else
                                break; // throw hiba("sugarfeldolgozo::get_cell_list", "No W intersection.");
                            kiv_ir = 1;
                        }
                        else if (akt_dx < 0) { // az E oldalon van a metszet
                            dx = -dx;
                            if (dir1 == 'E') {
                                aktualis_egyenes.p1.x = 2 * mp1.x - aktualis_egyenes.p1.x;
                                aktualis_egyenes.p2.x = 2 * mp1.x - aktualis_egyenes.p2.x;
                            }
                            else if (dir2 == 'E') {
                                aktualis_egyenes.p1.x = 2 * mp2.x - aktualis_egyenes.p1.x;
                                aktualis_egyenes.p2.x = 2 * mp2.x - aktualis_egyenes.p2.x;
                            }
                            else
                                break; // throw hiba("sugarfeldolgozo::get_cell_list", "No E intersection.");
                            kiv_ir = 1;
                        }
                        else if (akt_dy > 0) { // a S oldalon van a metszet
                            dy = -dy;
                            if (dir1 == 'S') {
                                aktualis_egyenes.p1.y = 2 * mp1.y - aktualis_egyenes.p1.y;
                                aktualis_egyenes.p2.y = 2 * mp1.y - aktualis_egyenes.p2.y;
                            }
                            else if (dir2 == 'S') {
                                aktualis_egyenes.p1.y = 2 * mp2.y - aktualis_egyenes.p1.y;
                                aktualis_egyenes.p2.y = 2 * mp2.y - aktualis_egyenes.p2.y;
                            }
                            else
                                break; // throw hiba("sugarfeldolgozo::get_cell_list", "No S intersection.");
                            kiv_ir = 2;
                        }
                        else if (akt_dy < 0) { // a N oldalon van a metszet
                            dy = -dy;
                            if (dir1 == 'N') {
                                aktualis_egyenes.p1.y = 2 * mp1.y - aktualis_egyenes.p1.y;
                                aktualis_egyenes.p2.y = 2 * mp1.y - aktualis_egyenes.p2.y;
                            }
                            else if (dir2 == 'N') {
                                aktualis_egyenes.p1.y = 2 * mp2.y - aktualis_egyenes.p1.y;
                                aktualis_egyenes.p2.y = 2 * mp2.y - aktualis_egyenes.p2.y;
                            }
                            else
                                break; // throw hiba("sugarfeldolgozo::get_cell_list", "No N intersection.");
                            kiv_ir = 2;
                        }
                        else if (akt_dz > 0) { // a B oldalon van a metszet
                            dz = -dz;
                            if (dir1 == 'B') {
                                aktualis_egyenes.p1.z = 2 * mp1.z - aktualis_egyenes.p1.z;
                                aktualis_egyenes.p2.z = 2 * mp1.z - aktualis_egyenes.p2.z;
                            }
                            else if (dir2 == 'B') {
                                aktualis_egyenes.p1.z = 2 * mp2.z - aktualis_egyenes.p1.z;
                                aktualis_egyenes.p2.z = 2 * mp2.z - aktualis_egyenes.p2.z;
                            }
                            else
                                break; // throw hiba("sugarfeldolgozo::get_cell_list", "No B intersection.");
                            kiv_ir = 3;
                        }
                        else { // a T oldalon van a metszet
                            dz = -dz;
                            if (dir1 == 'T') {
                                aktualis_egyenes.p1.z = 2 * mp1.z - aktualis_egyenes.p1.z;
                                aktualis_egyenes.p2.z = 2 * mp1.z - aktualis_egyenes.p2.z;
                            }
                            else if (dir2 == 'T') {
                                aktualis_egyenes.p1.z = 2 * mp2.z - aktualis_egyenes.p1.z;
                                aktualis_egyenes.p2.z = 2 * mp2.z - aktualis_egyenes.p2.z;
                            }
                            else
                                break; // throw hiba("sugarfeldolgozo::get_cell_list", "No T intersection.");
                            kiv_ir = 3;
                        }
                        if (!metszesszamolo(aktualis_egyenes, modell_racs.getconstref(akt_x, akt_y, akt_z).sarkok, d, dir1, dir2, mp1, mp2, itt_irany, kiv_ir)) {
                            break; // 
                            printf("\n\nError: Reflected light problem.\n");
                            exit(1);
                        }
                        is_vertical = aktualis_egyenes.get_vertical_angle() < 0.1;
                    }
                }

                // mindkettõ => cella berakása (tükrözésnél ugyanazt)

                if (!is_vege) {
                    t_modell_cella & akt_cella = modell_racs.getref(akt_x, akt_y, akt_z);
                    if (d > min_d && akt_cella.is_cella && akt_cella.pmat->is_fenypor) {
                        uj_tomb.getLast().K_P *= (float)prev_K_P; // !! Az elõzõ cella kimenete
                        akt_P_szorzo *= prev_K_P;
                        uj_tomb.reallocate(uj_tomb.size() + 1);
                        sugar_cella & uj = uj_tomb.getLast();
                        uj.cella_index = akt_cella.cella_index;
                        uj.x = akt_x;
                        uj.y = akt_y;
                        uj.z = akt_z;
                        uj.d = (float)d;
                        uj.K_P = (float)1;
                        uj.dir1 = dir1;
                        uj.dir2 = dir2;
                        uj.is_vertical = is_vertical;
                        sum_d += d;
                        uj.is_oldalsohoz = p_sugar->is_oldalso;
                        if (is_blue) akt_cella.sum_P_kek += P0 * akt_P_szorzo; // !!
                        else akt_cella.sum_P_sarga += P0 * akt_P_szorzo;
                    }
                }
            }
        }

    };

    //const sugar_cella & last = uj_tomb.getLast();
    //if (last.z == 7 && last.x < 20 && last.x > 0 && last.y>20 && last.y < 45) {
    //    printf("\n%u, %u, %u\n", last.x, last.y, last.z);
    //}

    return sum_d;
}


//***********************************************************************
struct sulyhossz {
//***********************************************************************
    dbl hossz, suly;
};


//***********************************************************************
dbl get_dlpp(const tomb<sulyhossz> & sh, dbl kivant_arany, dbl tenyleges_arany, dbl sum_P, dbl d_light_powder) {
// kivant_arany : ilyen kék arányt akarunk kapni a szimulációban
// tenyleges_arany : a mod fájlban lévõ fénypor alfa ekkora arányt ad a valóságban
//***********************************************************************
    if (tenyleges_arany > 0.98 || tenyleges_arany > 0.98)
        return 1;
    dbl min = 1, max = 1000000;
    dbl dlpp_uj = 0;
    while (max - min > 1e-4) {
        dbl alfa = 0.5*(min + max);
        dbl sum = 0;
        for (uns i = 0; i < sh.size(); i++) {
            sum += sh[i].suly*exp(-alfa*sh[i].hossz);
        }
        dbl akt_arany = sum / sum_P;
        //printf("\nki=%f, alfa=%f, v=%g", per, alfa, 1/v);
        if (akt_arany > kivant_arany)
            min = alfa;
        else
            max = alfa;
        dlpp_uj = -alfa*d_light_powder / log(tenyleges_arany);
    }
    return dlpp_uj;
}


//***********************************************************************
dbl sugarfeldolgozo::build_blue_cell_rays(const simulation & aktSim, tomb3d<t_modell_cella>& modell_racs, FILE * fp, vezetes_tomb_tipus & vt, uns & db_cell) {
// végigmegy a cellasugarak listán. Ahol a cellalista egy elemû, azaz junction,
// ott végigmegy az indulo_sugarak-on, és mindegyik alapján gyárt egy 
// cellasugarat, melyek utolsó cellája a junction cella
// az elsõ cellasugarat a junction cella helyére köti be
//***********************************************************************
    dbl sum_dP = 0, sum_P = 0;
    printf("write blue cell rays\n");
    tomb<sugar_cella> cellatomb;
    tomb<sulyhossz> sh;
    bool is_first = true;
    //**********************************************
    // Tükrözõdõ középvonalhoz (half_reflex) ezt true-ra (vagy aktSim.pmodel->is_half-ra) állítjuk
    bool is_half = false; // aktSim.pmodel->is_half;
    //bool is_half = true; // aktSim.pmodel->is_half;
    //**********************************************
    uns cell_db = 0;
    for (uns i = 0; i < kek_sugarlista.size(); i++) {
        sugar_lista & akt_sugarlista = kek_sugarlista[i];
        akt_sugarlista.kiindulo_cella.K_P = 1;
        printf(".");
        for (auto jp = akt_sugarlista.start.next; jp != &akt_sugarlista.stop; jp = jp->next) {
            sum_P += jp->P;
            cellatomb.reallocate(1);
            cellatomb[0] = akt_sugarlista.kiindulo_cella;
            cellatomb[0].K_P = (float)jp->P_kek_sarga;
            dbl sugarhossz = get_cell_list(true, is_half, cellatomb, jp, modell_racs, jp->P);
            sum_dP += sugarhossz * jp->P;
            sulyhossz s;
            s.hossz = sugarhossz;
            s.suly = jp->P;
            sh.add_with_reallocate(s);
            //printf("%g\t%g\t%g\t%g\t%g\n", sugarhossz, jp->P, sugarhossz * jp->P, sum_dP, sum_P);
            write_one_path(aktSim, modell_racs, fp, vt, cellatomb, true, is_first);
            cell_db += cellatomb.size();
            is_first = false;
        }
    }
    kek_sugarlista.clear();
    dlppy = dlppb = get_dlpp(sh, blue_tenyleges, blue_tenyleges, sum_P, d_light_powder);//d_light_powder * sum_P / sum_dP;
    //dlppy *= 0.5;
    printf("\nd_light_powder = %g, d_blue_weighted = %g, dlpp = %g\n", d_light_powder, sum_dP / sum_P, dlppb);
    fprintf(logfile, "d_light_powder=%g\nd_blue_weighted=%g\ndlpp=%g\n", d_light_powder*1e6, 1e6*sum_dP / sum_P, dlppb);
    db_cell = cell_db;
    return sum_dP / sum_P;
}


//***********************************************************************
dbl sugarfeldolgozo::build_yellow_cell_rays(const simulation & aktSim, tomb3d<t_modell_cella>& modell_racs, FILE * fp, vezetes_tomb_tipus & vt, uns & db_cell) {
// végigmegy a cellasugarak listán. Ahol a cellalista egy elemû, azaz junction,
// ott végigmegy az indulo_sugarak-on, és mindegyik alapján gyárt egy 
// cellasugarat, melyek utolsó cellája a junction cella
// az elsõ cellasugarat a junction cella helyére köti be
//***********************************************************************
    dbl sum_dP = 0, sum_P = 0;
    printf("write yellow cell rays\n");
    tomb<sugar_cella> cellatomb;
    tomb<sulyhossz> sh;
    bool is_first = true;
    //**********************************************
    // Tükrözõdõ középvonalhoz (half_reflex) ezt true-ra (vagy aktSim.pmodel->is_half-ra) állítjuk
    bool is_half = false; // aktSim.pmodel->is_half;
    //bool is_half = true; // aktSim.pmodel->is_half;
    //**********************************************
    uns dbdb = 0;
    uns cell_db = 0;
    for (uns i = 0; i < sarga_sugarlista.size(); i++) {
        sugar_lista & akt_sugarlista = sarga_sugarlista[i];
        akt_sugarlista.kiindulo_cella.K_P = 1;
        if (dbdb++ % 10 == 0) printf(".");
        uns hanyadik = 0;
        for (auto jp = akt_sugarlista.start.next; jp != &akt_sugarlista.stop; jp = jp->next) {
            sum_P += jp->P;
            cellatomb.reallocate(1);
            cellatomb[0] = akt_sugarlista.kiindulo_cella;
            cellatomb[0].K_P = (float)jp->P_kek_sarga;
            //if (cellatomb[0].K_P > 0.5)
            //    printf("%g\n", cellatomb[0].K_P);
//            if (i == 341 && jp->src_cella_index == 15080 && hanyadik==1)
//                printf("*");
            dbl sugarhossz = get_cell_list(false, is_half, cellatomb, jp, modell_racs, jp->P);
            sum_dP += sugarhossz * jp->P;
            sulyhossz s;
            s.hossz = sugarhossz;
            s.suly = jp->P;
            sh.add_with_reallocate(s);
            write_one_path(aktSim, modell_racs, fp, vt, cellatomb, false, is_first);
            cell_db += cellatomb.size();
            is_first = false;
            hanyadik++;
        }
    }
    uns sh_index = 0;
    sum_dP = 0;
    sum_P = 0;
    for (uns i = 0; i < sarga_sugarlista.size(); i++) {
        sugar_lista & akt_sugarlista = sarga_sugarlista[i];
        for (auto jp = akt_sugarlista.start.next; jp != &akt_sugarlista.stop; jp = jp->next) {
            const t_modell_cella & aktCella = modell_racs.getconstref(jp->x0, jp->y0, jp->z0);
            dbl suly = jp->P * aktCella.sum_P_sarga;
            sum_P += suly;
            sum_dP += sh[sh_index].hossz * suly;
            sh[sh_index].suly = suly;
            sh_index++;
        }
    }
    sarga_sugarlista.clear();
    dbl kivant_arany = yellow_tenyleges + yellow_correction;
    if (kivant_arany > 1)
        kivant_arany = 1;
    dlppy = get_dlpp(sh, kivant_arany, yellow_tenyleges, sum_P, d_light_powder);
    //dlppy = /*0.5 * */d_light_powder * sum_P / sum_dP; // 0.5*, mert átlagosan feleakkora utat tesz meg a sárgában, mint a kék
    printf("\nd_light_powder = %g, d_yellow_weighted = %g, dlppy = %g\n", d_light_powder, sum_dP / sum_P, dlppy);
    fprintf(logfile, "d_light_powder=%g\nd_yellow_weighted=%g\ndlppy=%g\n", 1e6*d_light_powder, 1e6*sum_dP / sum_P, dlppy);
    db_cell = cell_db;
    return sum_dP / sum_P;
}


//***********************************************************************
void sugarfeldolgozo::build_and_write(const simulation & aktSim, tomb3d<t_modell_cella> & modell_racs, uns dir_mode, char output_side, 
    uns ray_per_cell_dir, bool is_sarga_szetmegy, dbl cut_level, dbl d_light_powder, bool is_top_junction, uns junction_layer, FILE * fp, 
    vezetes_tomb_tipus & vt, dbl blue_tenyleges, dbl yellow_tenyleges, dbl yellow_correction) {
// is_top_junction: top vagy bottom a világítás iránya? Amerre több junction van, az lesz.
// junction layer: ahol a legtöbb junction van, az a z irányú index. Ez alatt/fölött nem megy ki kék sugár, sárga mehet.
// cellaútvonal meghatározása: melyik oldalon megy ki az aktuális cellából, az dönti el, hogy mi lesz a köv cella
// súlyozott kék hossz meghatározása
// a metszett cellák közül is eldobjuk az alig metszetteket
// ha az utolsó elõtti cella külsõ, akkor eldobjuk a sugarat (kívülrõl megy be a fény)
//***********************************************************************
    vt.clear();
    if (dir_mode == 0) {
        fprintf(fp, "NB0\nNY0\n");
        return;
    }
    this->dir_mode = dir_mode;
    this->ray_per_cell_dir = ray_per_cell_dir;
    this->cut_level = cut_level;
    this->d_light_powder = d_light_powder;
    this->is_top_junction = is_top_junction;
    this->junction_layer = junction_layer;
    this->is_sarga_szetmegy = is_sarga_szetmegy;
    this->output_side = output_side;
    this->blue_tenyleges = blue_tenyleges;
    this->yellow_tenyleges = yellow_tenyleges;
    this->yellow_correction = yellow_correction;
    clear();
    kulso_oldalak_gyujtese(aktSim, modell_racs);
    dbl d_blue_weighted_1 = kek_sugarak_gyujtese(aktSim,modell_racs); // elõzetes súlyozott hossz a kék sugarak hossza alapján
    uns sugar_db = 0;
    for (uns i = 0; i < kek_sugarlista.size();i++)
        sugar_db += kek_sugarlista[i].db;
    if(output_side=='O')
        fprintf(fp, "WO\n");
    fprintf(fp, "NB%u\n", sugar_db);
    uns db_cell = 0;
    dbl d_blue_weighted_2 = build_blue_cell_rays(aktSim, modell_racs, fp, vt, db_cell); // tényleges súlyozott hossz a fénypor cellákon átmenõ hosszak alapján
    printf("blue ray cells = %s\n", unsprint(db_cell));
    fprintf(logfile, "blue ray cells=%s\n", unsprint(db_cell,true));
    if (is_sarga_szetmegy) {
        // csak akkor számoljuk ki a sárga sugarakat, ha kell
        dbl d_yellow_weighted_1 = sarga_sugarak_gyujtese(aktSim, modell_racs); // elõzetes súlyozott hossz a kék sugarak hossza alapján
        sugar_db = 0;
        for (uns i = 0; i < sarga_sugarlista.size();i++)
            sugar_db += sarga_sugarlista[i].db;
        fprintf(fp, "NY%u\n", sugar_db);
        db_cell = 0;
        dbl d_yellow_weighted_2 = build_yellow_cell_rays(aktSim, modell_racs, fp, vt, db_cell); // tényleges súlyozott hossz a fénypor cellákon átmenõ hosszak alapján
        fprintf(logfile, "yellow ray cells=%s\n", unsprint(db_cell,true));
    }
    else {
        fprintf(fp, "NY0\n");
    }
}


//***********************************************************************
void apa::build_modell_fa(simulation & aktSim){
//***********************************************************************
    del_modell_fa(modell_fa_el);
    del_modell_fa(modell_fa_th);
    del_modell_fa(modell_fa_elth);
    doboz d;
    d.x1 = d.y1 = d.z1 = 0;
    d.x2 = aktSim.pmodel->x_res - 1;
    d.y2 = aktSim.pmodel->y_res - 1;
    d.z2 = aktSim.pmodel->z_res - 1;
    tomb<Oldal> iranyok;
    red_fa::iranyszamolo_uj(d, iranyok);
    switch (aktSim.mezo) {
        case FieldEl:
            modell_fa_el = red_fa::build_tree_uj(FieldEl, modell_racs, d, iranyok, 0);
            red_fa::optimize_tree(modell_fa_el);
            red_fa::tol_ig_db_szamol(FieldEl, modell_fa_el);
            break;
        case FieldTherm:
            modell_fa_th = red_fa::build_tree_uj(FieldTherm, modell_racs, d, iranyok, 0);
            red_fa::optimize_tree(modell_fa_th);
            red_fa::tol_ig_db_szamol(FieldTherm, modell_fa_th);
            break;
        case FieldElTherm: {
                switch (aktSim.mezo_szamitasi_mod) {
                    case 1:
                    case 2:
                        modell_fa_el = red_fa::build_tree_uj(FieldEl, modell_racs, d, iranyok, 0);
                        red_fa::optimize_tree(modell_fa_el);
                        red_fa::tol_ig_db_szamol(FieldEl, modell_fa_el);
                        modell_fa_th = red_fa::build_tree_uj(FieldTherm, modell_racs, d, iranyok, 0);
                        red_fa::optimize_tree(modell_fa_th);
                        red_fa::tol_ig_db_szamol(FieldTherm, modell_fa_th);
                        break;
                    case 3:
                        modell_fa_elth = red_fa::build_tree_uj(FieldElTherm, modell_racs, d, iranyok, 0);
                        red_fa::optimize_tree(modell_fa_elth);
                        red_fa::tol_ig_db_szamol(FieldElTherm, modell_fa_elth);
                        break;
                }
            }
            break;
    }
}


//***********************************************************************
void apa::build_analizisek(simulation & aktSim){
//***********************************************************************
    aktAnalizisek.clear();
    for (uns i = 0; i < aktSim.tanal.size(); i++) {
        v6anal aktAnal;
        const analysis & anal = aktSim.tanal[i];
        if (aktAnalizisek.size() == 0) {
            aktAnal.maxiter = 100;
            aktAnal.maxhiba = 1e-4;
            aktAnal.I0 = 0.1;
            
            for (uns j = 0; j < colmax; j++) {
                if (aktSim.texcitE[j].is) {
                    v6gerj gerj;
                    const excitation & x = aktSim.texcitE[j];
                    gerj.ter = 'E';
                    gerj.tipus = 'N';
                    if (x.tipus == GerjU)
                        gerj.tipus = 'U';
                    if (x.tipus == GerjI)
                        gerj.tipus = 'I';
                    gerj.ertek = x.ertek;
                    gerj.color_index = aktSim.pmodel->tcolor[j].index;
                    aktAnal.gerj.add(gerj);
                }
                if (aktSim.texcitT[j].is) {
                    v6gerj gerj;
                    const excitation & x = aktSim.texcitT[j];
                    gerj.ter = 'T';
                    gerj.tipus = 'N';
                    if (x.tipus == GerjU)
                        gerj.tipus = 'T';
                    if (x.tipus == GerjI)
                        gerj.tipus = 'P';
                    gerj.ertek = x.ertek;
                    gerj.color_index = aktSim.pmodel->tcolor[j].index;
                    aktAnal.gerj.add(gerj);
                }
            }
            // feszültség probe-ok
            for (uns j = 0; j < aktSim.tproV.size(); j++) {
                const probeT & akt_probe = aktSim.tproV[j];
                const t_modell_cella & akt_cella = modell_racs.getconstref(akt_probe.x, akt_probe.y, akt_probe.z);
                v6eredm eredm;
                eredm.eredm_fajta = 'P';
                if(!akt_cella.is_cella)
                    throw hiba("apa::build_analizisek", "probe is set on a no-cell (%u, %u, %u)", akt_probe.x, akt_probe.y, akt_probe.z);
                if (akt_cella.belso_cellak.size() > 0) {
                    eredm.probe_cella_index = akt_cella.belso_cellak.getconstref(
                        akt_cella.belso_cellak.x_size() / 2, 
                        akt_cella.belso_cellak.y_size() / 2, 
                        akt_cella.belso_cellak.z_size() / 2).cella_index;
                }
                else {
                    eredm.probe_cella_index = akt_cella.cella_index;
                }
                eredm.probe_map_fajta = 'V'; // egyelõre mindig center értéket kér, oldalsót nem tud                
                aktAnal.eredm.add(eredm);
            }
            // hõmérséklet probe-ok
            for (uns j = 0; j < aktSim.tproT.size(); j++) {
                const probeT & akt_probe = aktSim.tproT[j];
                const t_modell_cella & akt_cella = modell_racs.getconstref(akt_probe.x, akt_probe.y, akt_probe.z);
                v6eredm eredm;
                eredm.eredm_fajta = 'P';
                if (!akt_cella.is_cella)
                    throw hiba("apa::build_analizisek", "probe is set on a no-cell (%u, %u, %u)", akt_probe.x, akt_probe.y, akt_probe.z);
                if (akt_cella.belso_cellak.size() > 0) {
                    eredm.probe_cella_index = akt_cella.belso_cellak.getconstref(
                        akt_cella.belso_cellak.x_size() / 2,
                        akt_cella.belso_cellak.y_size() / 2,
                        akt_cella.belso_cellak.z_size() / 2).cella_index;
                }
                else {
                    eredm.probe_cella_index = akt_cella.cella_index;
                }
                eredm.probe_map_fajta = 'V'; // egyelõre mindig center értéket kér, oldalsót nem tud
                bool megvan = false;
                for (uns k = 0; k < aktAnal.eredm.size(); k++) {
                    if (aktAnal.eredm[k].eredm_fajta == eredm.eredm_fajta
                        && aktAnal.eredm[k].probe_cella_index == eredm.probe_cella_index
                        && aktAnal.eredm[k].probe_map_fajta == eredm.probe_map_fajta)
                        megvan = true;
                }
                if(!megvan)
                    aktAnal.eredm.add(eredm);
            }
            if (!aktSim.is_nofim) {
                v6eredm eredm;
                eredm.eredm_fajta = 'M';
                eredm.probe_map_fajta = 'V';
                aktAnal.eredm.add(eredm);
            }
        }
        switch (anal.tipus) {
            case AnalDC: 
                aktAnal.tipus = 'D';
                aktAnalizisek.add(aktAnal);
                aktAnal.reset();
                break;
            case AnalNDC:
                aktAnal.tipus = 'D';
                aktAnal.maxiter = anal.ndc_maxiter;
                aktAnal.maxhiba = anal.relhiba;
                aktAnal.I0 = anal.ndc_I0;
                aktAnalizisek.add(aktAnal);
                aktAnal.reset();
                break;
            case AnalAC:
                aktAnal.tipus = 'A';
                aktAnal.ertek = anal.from;
                aktAnalizisek.add(aktAnal);
                aktAnal.reset();
                break;
            case AnalLogTran:
            case AnalBode:
            case AnalIdo: {
                    if (anal.tipus == AnalLogTran)
                        aktAnal.tipus = 'S';
                    else if (anal.tipus == AnalBode)
                        aktAnal.tipus = 'A';
                    else aktAnal.tipus = 'T';
                    cd ddf = pow(10.0, 1.0 / anal.step);
                    dbl aktTimeFreq, prevTimeFreq = 0;
                    for (aktTimeFreq = anal.from; aktTimeFreq <= anal.to + aktTimeFreq * 1.0e-6; aktTimeFreq *= ddf) {
                        aktAnal.ertek = aktTimeFreq - prevTimeFreq;
                        aktAnalizisek.add(aktAnal);
                        aktAnal.reset();
                        prevTimeFreq = aktTimeFreq;
                    }
                }
                break;
            case AnalLinTran:
                aktAnal.tipus = 'S';
                for (dbl aktTime = anal.from; aktTime <= anal.to + aktTime * 1.0e-6; aktTime += anal.from) {
                    aktAnal.ertek = anal.from;
                    aktAnalizisek.add(aktAnal);
                    aktAnal.reset();
                }
                break;
            case AnalCtrl: {
                    aktAnal.tipus = 'S';
                    dbl aktStep = anal.from;
                    dbl limit = anal.to;
                    uns next_change_index = 0;
                    dbl aktTime = 0;
                    while (aktTime < limit + aktStep*1e-6) { // nincs több változás
                        dbl prevTime = aktTime;
                        bool noadd = false;
                        if (next_change_index >= anal.ctrl.size()) {
                            aktTime += aktStep;
                            aktAnal.ertek = (aktTime < limit ? aktTime : limit) - prevTime;
                            if (aktAnal.ertek<aktStep*1e-3)
                                noadd = true;
                        }
                        else { // van még változás
                            if (aktTime + aktStep*1.000001 > anal.ctrl[next_change_index].time) { // változás van
                                const change_time & ct = anal.ctrl[next_change_index];
                                aktTime = ct.time;
                                if (aktTime >= limit + aktStep*1e-6) { // timelimit után lenne a következõ lépés a változással együtt
                                    aktAnal.ertek = limit - prevTime;
                                }
                                else { // tényleg változás van
                                    aktAnal.ertek = aktTime - prevTime;
                                    if (ct.timestep > 0)
                                        aktStep = ct.timestep;
                                    for (uns j = 0; j < ct.excit.size(); j++) {
                                        v6gerj gerj;
                                        const excitation_2 & x = ct.excit[j];
                                        gerj.ter = x.is_el ? 'E' : 'T';
                                        switch (x.tipus) {
                                            case GerjU: gerj.tipus = x.is_el ? 'U' : 'T'; gerj.ertek = x.ertek; break;
                                            case GerjI: gerj.tipus = x.is_el ? 'I' : 'P'; gerj.ertek = x.ertek; break;
                                            case GerjSemmi: gerj.tipus = 'N'; gerj.ertek = 0; break;
                                        }
                                        gerj.color_index = aktSim.pmodel->tcolor[x.color_index].index;
                                        aktAnal.gerj.add(gerj);
                                    }
                                }
                                next_change_index++;
                            }
                            else { // most nincs változás
                                aktTime += aktStep;
                                aktAnal.ertek = (aktTime < limit ? aktTime : limit) - prevTime;
                                if (aktAnal.ertek<aktStep*1e-3)
                                    noadd = true;
                            }
                        }
                        if(!noadd)
                            aktAnalizisek.add(aktAnal);
                        aktAnal.reset();
                    }
                }
                break;
        }
    }
}


//***********************************************************************
void apa::write_analizisek(FILE * fp, simulation & aktSim, uns junction_db){
//***********************************************************************
    build_analizisek(aktSim);
    FILE *descfp;
    if (fopen_s(&descfp, (path + aktSim.fimdesc_name + ".fimdesc").c_str(), "at") != 0)
        throw hiba("apa::write_analizisek", "cannot open %s to append", (path + aktSim.fimdesc_name + ".fimdesc").c_str());
    fprintf(descfp, "%u\n", aktAnalizisek.size());
    fclose(descfp);
    if (aktSim.is_vesszo)
        fprintf(fp, "UC\n");
    fprintf(fp, "NA%u\n", aktAnalizisek.size());
    for (uns i = 0; i < aktAnalizisek.size(); i++) {
        fprintf(fp, "BA%u\n", i+1);
        const v6anal & a = aktAnalizisek[i];
        for (uns j = 0; j < a.gerj.size(); j++) {
            fprintf(fp, "X%u%c%c", a.gerj[j].color_index, a.gerj[j].ter, a.gerj[j].tipus);
            if (a.gerj[j].tipus != 'N')
                fprintf(fp, "%g", a.gerj[j].ertek);
            fprintf(fp, "\n");
        }
        if(i==0)
            fprintf(fp, "TAMB%g\n", aktSim.ambiT);
        fprintf(fp, "A%c", a.tipus);
        if (a.tipus != 'D')
            fprintf(fp, "%g", a.ertek);
        fprintf(fp, "\n");
        if (a.I0 != 0)
            fprintf(fp, "PC%g\n", a.I0);
        if (a.maxhiba != 0)
            fprintf(fp, "PE%g\n", a.maxhiba);
        if (a.maxiter != 0)
            fprintf(fp, "PI%u\n", a.maxiter);
        if(a.is_del_excits)
            fprintf(fp, "DX\n");
        if (a.is_reset)
            fprintf(fp, "DA\n");
        uns probe_db = a.eredm.size();
        if (probe_db > 0) {
            for (uns j = 0; j < probe_db; j++) {
                const v6eredm & er = a.eredm[j];
                switch (er.eredm_fajta) {
                    case 'P':
                        if (er.probe_map_fajta == 'F') {
                            fprintf(fp, "RP%uF%u\n", er.probe_cella_index, er.probe_face_index);
                        }
                        else {
                            fprintf(fp, "RP%uC\n", er.probe_cella_index);
                        }
                        break;
                    case 'M': 
                        fprintf(fp, "RMC%c\n", er.probe_map_fajta);
                        if (junction_db > 0 && er.probe_map_fajta != 'R')
                            fprintf(fp, "RMCR\n");
                        break;
                }
            }
        }
        fprintf(fp, "EA%u\n", i+1);
    }
}


//***********************************************************************
void sugarfeldolgozo::write_one_path(const simulation & aktSim, const tomb3d<t_modell_cella>& modell_racs, FILE * fp, vezetes_tomb_tipus & vt, const tomb<sugar_cella> & cellatomb, bool is_blue, bool is_first) {
//***********************************************************************
    struct prev_values_t {
        tomb<sugar_cella> cellatomb;
        uns C0, B, E, Y, M, P;
        unsigned short F0;
        float KR0;
        dbl U, R, K;
        void reset() {
            cellatomb.clear();
            C0 = B = E = Y = M = P = 0;
            F0 = 0;
            KR0 = 0;
            U = R = K = 0;
        }
    };
    static prev_values_t prevs;
    if(is_first) prevs.reset();

    if (cellatomb.size() == 0)
        return;
    
    dbl d_max = 0;
    for (uns i = 1; i < cellatomb.size(); i++)
        if (cellatomb[i].d > d_max)
            d_max = cellatomb[i].d;

    uns P_max = (uns)((d_max + 5e-10)*1e9);
    if (is_first || prevs.P != P_max) { fprintf(fp, "P%u", P_max); prevs.P = P_max; }
    else fprintf(fp, "P0");
    if (is_first || prevs.C0 != cellatomb[0].cella_index) fprintf(fp, "C%u", prevs.C0 = cellatomb[0].cella_index);
    if (is_blue) {
        if (is_first || prevs.F0 !=  cellatomb[0].face_index)  fprintf(fp, "F%u",  prevs.F0 = cellatomb[0].face_index);
        if (is_first || prevs.KR0 != cellatomb[0].K_P)         fprintf(fp, "K%u", (uns)(((prevs.KR0 = cellatomb[0].K_P)+5e-10)*1e9));
    }
    else {
        if (is_first || prevs.KR0 != cellatomb[0].K_P)         fprintf(fp, "R%u", (uns)(((prevs.KR0 = cellatomb[0].K_P) + 5e-10)*1e9));
    }
    
    fprintf(fp, "N%u", cellatomb.size() - 1);
    //*******************************************************
    bool csak_az_oldalsohoz_meno_sugarakat_irja_ki = false; // !!!
    //*******************************************************
    for (uns i = 1; i < cellatomb.size(); i++) {
        const sugar_cella & akt = cellatomb[i];
        
        char significant = (akt.dir1 == output_side || akt.dir2 == output_side) ? 'S' : 'C';
        if (output_side == 'O' && (akt.dir1 == 'T' || akt.dir2 == 'T')) // OUT típusnál a normál fájlokba a TOP-nak megfelelõ íródik
            significant = 'S';
        if (csak_az_oldalsohoz_meno_sugarakat_irja_ki && !akt.is_oldalsohoz)significant = 'C';
        if (output_side == 'A') significant = 'S';
        else if (output_side == 'H') {
            significant = (akt.dir1 == 'T' || akt.dir2 == 'T') ? 'S' : 'H';
        }
        else if (output_side == 'V') significant = akt.is_vertical ? 'S' : 'C';

        if (!is_first && prevs.cellatomb.size() > i) {
            const sugar_cella * p_prev_cella = &prevs.cellatomb[i];
            char prev_significant = (p_prev_cella->dir1 == output_side || p_prev_cella->dir2 == output_side) ? 'S' : 'C';
            if (output_side == 'O' && (p_prev_cella->dir1 == 'T' || p_prev_cella->dir2 == 'T')) // OUT típusnál a normál fájlokba a TOP-nak megfelelõ íródik
                prev_significant = 'S';
            if (csak_az_oldalsohoz_meno_sugarakat_irja_ki && !p_prev_cella->is_oldalsohoz)prev_significant = 'C';
            if (output_side == 'A') prev_significant = 'S';
            else if (output_side == 'H') {
                prev_significant = (p_prev_cella->dir1 == 'T' || p_prev_cella->dir2 == 'T') ? 'S' : 'H';
            }
            else if (output_side == 'V') prev_significant = p_prev_cella->is_vertical ? 'S' : 'C';
            if (significant != prev_significant || p_prev_cella->cella_index != akt.cella_index) {
                fprintf(fp, "%c%u", significant, akt.cella_index);
            }
        }
        else {
            fprintf(fp, "%c%u", significant, akt.cella_index);
        }
        
        const t_modell_cella & akt_cella = modell_racs.getconstref(akt.x, akt.y, akt.z);
        
        if (is_blue) {
            uns B = vt.get_index(akt_cella.pmat->light_blue_absorption_coeff);
            uns E = vt.get_index(akt_cella.pmat->light_conversion_efficiency);
            if (is_first || prevs.B != B) fprintf(fp, "B%u", prevs.B = B);
            if (is_first || prevs.E != E) fprintf(fp, "E%u", prevs.E = E);
        }
        
        uns Y = vt.get_index(akt_cella.pmat->light_yellow_absorption_coeff);
        if (is_first || prevs.Y != Y) fprintf(fp, "Y%u", prevs.Y = Y);
        
        uns M = vt.get_index(akt_cella.pmat->light_re_conversion_efficiency);
        if (is_first || prevs.M != M) fprintf(fp, "M%u", prevs.M = M);

        if (is_blue) {
            dbl K = akt.K_P;
            if (is_first || prevs.K != K) fprintf(fp, "K%g;", prevs.K = K);
            dbl U = is_sarga_szetmegy ? 0.0 : 1.0; // Az új sárga nem csökken a reflexió miatt, ezért 1 marad
            if (is_first || prevs.U != U) fprintf(fp, "U%g;", prevs.U = U);
        }
        
        dbl R = is_blue ? (is_sarga_szetmegy ? 0.0 : akt.K_P) : akt.K_P;
        if (is_first || prevs.R != R) fprintf(fp, "R%g;", prevs.R = R);
        
        uns akt_D = (uns)((akt.d / d_max) * 32768 + 0.5);
        if (akt_D == 0) akt_D = 1;
        if (akt_D == 32768) akt_D = 0;
        fprintf(fp, "D%u", akt_D);
    }
    prevs.cellatomb = cellatomb;
    fprintf(fp, "\n");
}


//***********************************************************************
void apa::write_elemi_cells(FILE * fp, simulation & aktSim, uns cellaszam, uns csatlakozo_db){
//***********************************************************************
    meret_tomb.clear();
    fprintf(fp, "NL%uC%u\n", cellaszam, csatlakozo_db);
    FILE *descfp;
    if (fopen_s(&descfp, (path + aktSim.fimdesc_name + ".fimdesc").c_str(), "wt") != 0)
        throw hiba("apa::write_elemi_cells", "cannot open %s to write", (path + aktSim.fimdesc_name + ".fimdesc").c_str());
    fprintf(descfp, "%u\t%u\t%u\t%u\n", cellaszam, modell_racs.x_size(), modell_racs.y_size(), modell_racs.z_size());
    fprintf(descfp, "%u\t%u\t%u\t%u\t%u\n", aktSim.FIM_res_xy, aktSim.FIM_res_z, aktSim.FIM_diff_x, aktSim.FIM_diff_y, aktSim.FIM_diff_z);
    for (uns z = 0; z < modell_racs.z_size(); z++) {
        for (uns y = 0; y < modell_racs.y_size(); y++) {
            for (uns x = 0; x < modell_racs.x_size(); x++) {
                modell_racs.getref(x, y, z).write_cella(fp, descfp, meret_tomb, x, y, z);
            }
        }
    }
    fclose(descfp);
}


//***********************************************************************
void apa::index_and_write_elemi_cells(FILE *fp, simulation & aktSim, uns cellaszam, uns csatlakozo_db){
//***********************************************************************
    for (uns z = 0; z < modell_racs.z_size(); z++)
        for (uns y = 0; y < modell_racs.y_size(); y++)
            for (uns x = 0; x < modell_racs.x_size(); x++) {
                modell_racs.getref(x, y, z).set_face_indexek();
            }

    model & aktMod = *aktSim.pmodel;
    cuns x_res = aktMod.x_res;
    cuns y_res = aktMod.y_res;
    cuns z_res = aktMod.z_res;
    
    // fénypor vertical_light_conversion beállítása

    uns dir_mode = 0; // egyféle mód megengedett az egész modellben. Ha több van, az elsõt veszi figyelembe.
    uns ray_per_cell_dir = 1;
    bool is_sarga_szetmegy; // vannak-e külön sárga sugarak
    char output_side = 'T';
    dbl cut_level = 0;
    dbl d_light_powder = 0;
    uns top_sugarzo = 0; // hány top és hány bottom irányban sugárzó junction van?
    uns bottom_sugarzo = 0;
    dbl legfelso_fenypor = 0, legalso_fenypor = 0; // Az alap távolság automatikus meghatározásához
    material *pmat = nullptr;
    tomb<uns> db_junction; // hány junction van az adott rétegben?
    dbl yellow_correction = 0;
    db_junction.resize(z_res);
    for (uns z = 0; z < z_res; z++)
        db_junction[z] = 0;

    for (uns z = 0; z < z_res; z++)
        for (uns y = 0; y < y_res; y++)
            for (uns x = 0; x < x_res; x++) {
                t_modell_cella & akt_cella = modell_racs.getref(x, y, z);
                // az anyagból megnézni, hogy ez fénypor - e
                if (akt_cella.is_cella && akt_cella.pmat != nullptr && akt_cella.pmat->is_fenypor && akt_cella.pmat->light_vertical_light_conversion.g[0] == 0) {
                    if (dir_mode == 0) {
                        dir_mode = akt_cella.pmat->direction_mode;
                        ray_per_cell_dir = akt_cella.pmat->ray_per_cell_dir;
                        cut_level = akt_cella.pmat->cut_level;
                        output_side = akt_cella.pmat->output_side;
                        d_light_powder = akt_cella.pmat->d_light_powder;
                        is_sarga_szetmegy = akt_cella.pmat->is_sarga_szetmegy;
                        legalso_fenypor = akt_cella.sarkok.z0;
                        legfelso_fenypor = akt_cella.sarkok.z1;
                        pmat = akt_cella.pmat;
                        yellow_correction = akt_cella.pmat->yellow_correction;
                    }
                    else
                        legfelso_fenypor = akt_cella.sarkok.z1;
                }
                if (akt_cella.is_cella && akt_cella.is_junction) {
                    if (akt_cella.junction_bottom_face > 0)
                        bottom_sugarzo++;
                    if (akt_cella.junction_top_face > 0)
                        top_sugarzo++;
                    db_junction[z]++;
                }
                if (akt_cella.is_cella && akt_cella.pmat != nullptr && akt_cella.pmat->is_fenypor && akt_cella.pmat->light_vertical_light_conversion.g[0] != 0) {
                    // ha fénypor, akkor megkeresni a hozzá tartozó sugárzó cellát, indexét eltárolni; ha high_res, elkaszálni;
                    if(akt_cella.belso_cellak.size()!=0)
                        throw hiba("apa::build_modell_racs", "light powder in high res region not supported: (%u,%u,%u)", x, y, z);
                    for (uns i = 0; i < z_res; i++) {
                        t_modell_cella & sugarzo_cella = modell_racs.getref(x, y, i);
                        if (sugarzo_cella.is_junction) {
                            if (sugarzo_cella.belso_cellak.size() != 0)
                                throw hiba("apa::build_modell_racs", "light source in high res region with light powder is not supported: (%u,%u,%u)", x, y, i);
                            akt_cella.fenyforras_cella_index.add(sugarzo_cella.cella_index);
                        }
                    }
                }
            }

    printf("writing cells starts\n");
    write_elemi_cells(fp, aktSim, cellaszam, csatlakozo_db);

    // fénypor útvonalak beállítása (lehetne együtt a vertical_light_conversion résszel, de ha azt törlöm, egyszerûbb így, szétszedve)
    // => vigyázz, a dir_mode-ot is az nézi meg.

    vezetes_tomb_tipus vt;
    if (dir_mode > 0) {
        uns junction_layer = 0;
        for (uns z = 1; z < z_res; z++)
            if (db_junction[z]>db_junction[junction_layer])
                junction_layer = z;
        if (d_light_powder == 0) {
            if (top_sugarzo >= bottom_sugarzo) { // felfelé sugároz
                dbl junc = aktMod.z_hossz[2 * junction_layer + 1];
                d_light_powder = legfelso_fenypor - (junc > legalso_fenypor? junc : legalso_fenypor);
            }
            else {
                dbl junc = aktMod.z_hossz[2 * junction_layer - 1];
                dbl teto = junc < legfelso_fenypor ? junc : legfelso_fenypor;
                d_light_powder = junction_layer == 0 ? 0 : teto - legalso_fenypor;
            }
        }
        dbl kek_ki = 0;
        dbl sarga_ki = 0;
        if (pmat != nullptr) {
            dbl alfa_blue = pmat->light_blue_absorption_coeff.get_ertek(50);
            dbl alfa_yellow = pmat->light_yellow_absorption_coeff.get_ertek(50);
            dbl eta_conv = pmat->light_conversion_efficiency.get_ertek(50);
            dbl eta_re = pmat->light_re_conversion_efficiency.get_ertek(50);
            kek_ki = exp(-alfa_blue*d_light_powder);
            dbl sarga_arany = exp(-alfa_yellow*d_light_powder);
            dbl sarga_re = (1 - sarga_arany)*eta_re;
            sarga_ki = sarga_arany + sarga_re;
        }
        sugar_feldolgozo.build_and_write(aktSim, modell_racs, dir_mode, output_side, ray_per_cell_dir, is_sarga_szetmegy,
            cut_level, d_light_powder, top_sugarzo >= bottom_sugarzo, junction_layer, fp, vt, kek_ki, sarga_ki, yellow_correction);
    }

    // fényút tulajdonságok

    fprintf(fp, "NR%zu\n", vt.get_vector().size() - 1);
    for (size_t i = 1; i < vt.get_vector().size(); i++) {
        fprintf(fp, "R%u=", (uns)i);
        vt.get_vector()[i].write_normal(fp, 0, true);
    }
}


//***********************************************************************
void apa::write_meretek(FILE * fp){
//***********************************************************************
    fprintf(fp, "ND%zu\n", meret_tomb.get_vector().size() - 1);
    for (size_t i = 1; i < meret_tomb.get_vector().size(); i++)
        fprintf(fp, "D%zu=%e\n", i, meret_tomb.get_vector()[i]);
}


//***********************************************************************
void apa::write_modell_tree(FILE * fp, simulation & aktSim){
//***********************************************************************
    red_fa *modell_fa_1 = nullptr, *modell_fa_2 = nullptr;
    switch (aktSim.mezo) {
        case FieldEl:
            modell_fa_1 = modell_fa_el;
            break;
        case FieldTherm:
            modell_fa_1 = modell_fa_th;
            break;
        case FieldElTherm: {
                switch (aktSim.mezo_szamitasi_mod) {
                    case 1:
                    case 2:
                        modell_fa_1 = modell_fa_el;
                        modell_fa_2 = modell_fa_th;
                        break;
                    case 3:
                        modell_fa_1 = modell_fa_elth;
                        break;
                }
            }
            break;
    }
    uns start_index = 1;
    uns db1 = red_fa::set_indexek(modell_fa_1, start_index);
    uns db2 = red_fa::set_indexek(modell_fa_2, start_index);
    if (db1 + db2 != start_index - 1)
        throw hiba("apa::write_modell_tree", "db1 + db2 != start_index - 1");
    if (db2 != 0)
        fprintf(fp, "NF2B%u\n", db1 + db2);
    else
        fprintf(fp, "NF1B%u\n", db1);
    fprintf(fp, "G1G%u\n", db1);
    if (db2 != 0)
        fprintf(fp, "G2G%u\n", db1 + db2);
    red_fa::write_tree(fp, modell_fa_1);
    red_fa::write_tree(fp, modell_fa_2);
}


//***********************************************************************
void apa::del_modell_fa(red_fa * & fa){
//***********************************************************************
    if (fa == nullptr)
        return;
    del_modell_fa(fa->bal);
    del_modell_fa(fa->jobb);
    delete fa;
    fa = nullptr;
}


//***********************************************************************
void t_modell_face_adat::face2face(uns xx, uns yy, bool is_el, bool is_th){
//***********************************************************************
    belso_facek.resize(xx, yy);
    cuns xy = xx*yy;
    if (kulso_el_db > 1 || kulso_th_db > 1)
        throw hiba("t_modell_face_adat::face2face", "kulso_el_db > 1 || kulso_th_db > 1 (%u, %u)", kulso_el_db, kulso_th_db);
    kulso_el_db = is_el_perem || !is_el ? 0 : xy;
    kulso_th_db = is_th_perem || !is_th ? 0 : xy;
    for (uns y = 0; y < yy; y++)
        for (uns x = 0; x < xx; x++) {
            t_modell_face_adat & aktFace = belso_facek.getref(x, y);
            aktFace.A = A / xy;
            aktFace.anyag_index = anyag_index;
            aktFace.el_perem_index = el_perem_index;
            aktFace.is_el_perem = is_el_perem;
            aktFace.is_th_perem = is_th_perem;
            aktFace.junction_index = junction_index;
            aktFace.kulso_el_db = is_el_perem || !is_el ? 0 : 1;
            aktFace.kulso_th_db = is_th_perem || !is_th ? 0 : 1;
            aktFace.L = L;
            aktFace.th_perem_index = th_perem_index;
            aktFace.th_perem_x = th_perem_x;
            aktFace.th_perem_y = th_perem_y;
            aktFace.th_perem_c = th_perem_c;
        }
}

//***********************************************************************
void t_modell_face_adat::write_face(FILE * fp, bool is_el, bool is_th, uns & sfi, meret_tomb_tipus & meret_tomb, 
    uns cella_anyag_index, uns cella_tipus){
// sfi = start face index
//***********************************************************************
    if (belso_facek.size() > 0) {
        for (uns y = 0; y < belso_facek.y_size(); y++)
            for (uns x = 0; x < belso_facek.x_size(); x++)
                belso_facek.getref(x, y).write_face(fp, is_el, is_th, sfi, meret_tomb, cella_anyag_index, cella_anyag_index);
    }
    else {
        if (is_el && is_th) { // TODO: és nem centroid egyik sem
            bool is_dual = !is_el_perem && !is_th_perem && junction_index==0;

            // Elektromos

            if (face_index_el != sfi)
                throw hiba("write_face", "face_index_el != sfi (%u!=%u)", face_index_el, sfi);

                                                    fprintf(fp, "%c%uE%u", is_el_perem ? 'P' : 'J', sfi, csatlakozo_index_el);
            if (cella_anyag_index != anyag_index)   fprintf(fp, "M%u", anyag_index);
            if (junction_index > 0)                 fprintf(fp, "J%u", junction_index);
                                                    fprintf(fp, "%c%zuL%zuF%u", (th_perem_c == 0 ? 'A' : th_perem_c), meret_tomb.get_index(A), meret_tomb.get_index(L), sfi + 1);
            if (is_el_perem)                        fprintf(fp, "I%u", el_perem_index);
                                                    fprintf(fp, "\n");
            sfi++;

            // Termikus

            if (face_index_th != sfi)
                throw hiba("write_face", "face_index_th != sfi (%u!=%u)", face_index_th, sfi);

            if (!is_th_perem && is_dual) {
                fprintf(fp, "J%uT%uD\n", sfi, csatlakozo_index_th);
            }
            else {
                                                                fprintf(fp, "%c%uT%u", is_th_perem ? ((th_perem_index > 0 && th_perem_x != ~0) ? 'S' : 'P') : 'J', sfi, csatlakozo_index_th);
                if (cella_anyag_index != anyag_index)           fprintf(fp, "M%u", anyag_index);
                                                                fprintf(fp, "%c%zuL%zuF%u", (th_perem_c == 0 ? 'A' : th_perem_c), meret_tomb.get_index(A), meret_tomb.get_index(L), sfi - 1);
                if (is_th_perem) {
                                                                fprintf(fp, "I%u", th_perem_index);
                    if (th_perem_index > 0 && th_perem_x != ~0) fprintf(fp, "X%uY%u%c", th_perem_x, th_perem_y, th_perem_c);
                }
                                                                fprintf(fp, "\n");
            }
            sfi++;
        }
        else { //  vagy csak elektromos, vagy csak termikus a cella, azaz nincs párja
            if (is_el && face_index_el != sfi)
                throw hiba("write_face", "face_index_el != sfi (%u!=%u)", face_index_el, sfi);
            if (is_th && face_index_th != sfi)
                throw hiba("write_face", "face_index_th != sfi (%u!=%u)", face_index_th, sfi);

            // TODO: ha centroid, azt a rövidet kiírni
            bool is_spec_th_perem = is_th && is_th_perem && th_perem_index > 0 && th_perem_x != ~0;
            char perembetu = is_spec_th_perem ? 'S' : ((is_el && is_el_perem) || (is_th && is_th_perem) ? 'P' : 'J');
                                                                    fprintf(fp, "%c%u%c%u", perembetu, sfi, is_el ? 'E' : 'T', is_el ? csatlakozo_index_el : csatlakozo_index_th);
            if (cella_anyag_index != anyag_index)                   fprintf(fp, "M%u", anyag_index);
            if (is_el && junction_index > 0)                        fprintf(fp, "J%u", junction_index);
                                                                    fprintf(fp, "%c%zuL%zuX", (th_perem_c == 0 ? 'A' : th_perem_c), meret_tomb.get_index(A), meret_tomb.get_index(L));
            if ((is_el && is_el_perem) || (is_th && is_th_perem))   fprintf(fp, "I%u", is_el ? el_perem_index : th_perem_index);
            if (is_spec_th_perem)                                   fprintf(fp, "X%uY%u%c", th_perem_x, th_perem_y, th_perem_c);
                                                                    fprintf(fp, "\n");
            sfi++;
        }
    }
}


//***********************************************************************
uns t_modell_face_adat::facetipus_azonosito(bool is_el, bool is_th, uns cella_anyag_index) {
// 1: th, 2: el, 3: elth
// 4: th_perem, 8: el_perem
// 16: el_junction
// 32: más anyag
//***********************************************************************
    if (belso_facek.size() > 0)
        return 0;
    uns tipus = 0;
    if (is_th)
        tipus |= 1;
    if(is_el)
        tipus |= 2;
    if (is_th) {
        if (is_th_perem)
            tipus |= 4;
    }
    if (is_el) {
        if (is_el_perem)
            tipus |= 8;
        if (junction_index > 0)
            tipus |= 16;
    }
    if (cella_anyag_index != anyag_index)
        tipus |= 32;
    return tipus;
}


//***********************************************************************
void t_modell_face_adat::set_face_indexek(bool is_el, bool is_th, uns & sfi){
// sfi = start face index
//***********************************************************************
    if (belso_facek.size() > 0) {
        for (uns y = 0; y < belso_facek.y_size(); y++)
            for (uns x = 0; x < belso_facek.x_size(); x++)
                belso_facek.getref(x, y).set_face_indexek(is_el, is_th, sfi);
    }
    else {
        if (is_el && is_th) { // TODO: és nem centroid egyik sem

            // Elektromos

            face_index_el = sfi;
            sfi++;

            // Termikus

            face_index_th = sfi;
            sfi++;
        }
        else { //  vagy csak elektromos, vagy csak termikus a cella, azaz nincs párja
            // TODO: ha centroid, azt a rövidet kiírni
            if (is_el)face_index_el = sfi;
            if (is_th)face_index_th = sfi;
            sfi++;
        }
    }
}


//***********************************************************************
void t_modell_cella::set_belso_cellak_mul(uns x_mul, uns y_mul, uns z_mul) {
//***********************************************************************
    if (x_mul == 0)
        return;

    belso_cellak.resize(x_mul, y_mul, z_mul);
    mul_face(y_mul*z_mul, x_mul*z_mul, x_mul*y_mul);

    for (uns z = 0; z < belso_cellak.z_size(); z++)
        for (uns y = 0; y < belso_cellak.y_size(); y++)
            for (uns x = 0; x < belso_cellak.x_size(); x++)
                belso_cellak.getref(x, y, z).set_egy_belso_cella(*this);

    if (face_adat[WEST].is_el_perem || face_adat[WEST].is_th_perem) {
        for (uns z = 0; z < belso_cellak.z_size(); z++)
            for (uns y = 0; y < belso_cellak.y_size(); y++) {
                if (face_adat[WEST].is_el_perem)belso_cellak.getref(0, y, z).set_face_el(WEST, true, face_adat[WEST].el_perem_index, 0);
                if (face_adat[WEST].is_th_perem)belso_cellak.getref(0, y, z).set_face_th(WEST, true, face_adat[WEST].th_perem_index, 0, face_adat[WEST].th_perem_x, face_adat[WEST].th_perem_y, face_adat[WEST].th_perem_c);
            }
    }
    if (face_adat[EAST].is_el_perem || face_adat[EAST].is_th_perem) {
        for (uns z = 0; z < belso_cellak.z_size(); z++)
            for (uns y = 0; y < belso_cellak.y_size(); y++) {
                if (face_adat[EAST].is_el_perem)belso_cellak.getref(belso_cellak.x_size() - 1, y, z).set_face_el(EAST, true, face_adat[EAST].el_perem_index, 0);
                if (face_adat[EAST].is_th_perem)belso_cellak.getref(belso_cellak.x_size() - 1, y, z).set_face_th(EAST, true, face_adat[EAST].th_perem_index, 0, face_adat[EAST].th_perem_x, face_adat[EAST].th_perem_y, face_adat[EAST].th_perem_c);
            }
    }
    if (face_adat[SOUTH].is_el_perem || face_adat[SOUTH].is_th_perem) {
        for (uns z = 0; z < belso_cellak.z_size(); z++)
            for (uns x = 0; x < belso_cellak.x_size(); x++) {
if (face_adat[SOUTH].is_el_perem)belso_cellak.getref(x, 0, z).set_face_el(SOUTH, true, face_adat[SOUTH].el_perem_index, 0);
if (face_adat[SOUTH].is_th_perem)belso_cellak.getref(x, 0, z).set_face_th(SOUTH, true, face_adat[SOUTH].th_perem_index, 0, face_adat[SOUTH].th_perem_x, face_adat[SOUTH].th_perem_y, face_adat[SOUTH].th_perem_c);
            }
    }
    if (face_adat[NORTH].is_el_perem || face_adat[NORTH].is_th_perem) {
        for (uns z = 0; z < belso_cellak.z_size(); z++)
            for (uns x = 0; x < belso_cellak.x_size(); x++) {
                if (face_adat[NORTH].is_el_perem)belso_cellak.getref(x, belso_cellak.y_size() - 1, z).set_face_el(NORTH, true, face_adat[NORTH].el_perem_index, 0);
                if (face_adat[NORTH].is_th_perem)belso_cellak.getref(x, belso_cellak.y_size() - 1, z).set_face_th(NORTH, true, face_adat[NORTH].th_perem_index, 0, face_adat[NORTH].th_perem_x, face_adat[NORTH].th_perem_y, face_adat[NORTH].th_perem_c);
            }
    }
    if (face_adat[BOTTOM].is_el_perem || face_adat[BOTTOM].is_th_perem) {
        for (uns y = 0; y < belso_cellak.y_size(); y++)
            for (uns x = 0; x < belso_cellak.x_size(); x++) {
                if (face_adat[BOTTOM].is_el_perem)belso_cellak.getref(x, y, 0).set_face_el(BOTTOM, true, face_adat[BOTTOM].el_perem_index, 0);
                if (face_adat[BOTTOM].is_th_perem)belso_cellak.getref(x, y, 0).set_face_th(BOTTOM, true, face_adat[BOTTOM].th_perem_index, 0, face_adat[BOTTOM].th_perem_x, face_adat[BOTTOM].th_perem_y, face_adat[BOTTOM].th_perem_c);
            }
    }
    if (face_adat[TOP].is_el_perem || face_adat[TOP].is_th_perem) {
        for (uns y = 0; y < belso_cellak.y_size(); y++)
            for (uns x = 0; x < belso_cellak.x_size(); x++) {
                if (face_adat[TOP].is_el_perem)belso_cellak.getref(x, y, belso_cellak.z_size() - 1).set_face_el(TOP, true, face_adat[TOP].el_perem_index, 0);
                if (face_adat[TOP].is_th_perem)belso_cellak.getref(x, y, belso_cellak.z_size() - 1).set_face_th(TOP, true, face_adat[TOP].th_perem_index, 0, face_adat[TOP].th_perem_x, face_adat[TOP].th_perem_y, face_adat[TOP].th_perem_c);
            }
    }
    for (uns z = 0; z < belso_cellak.z_size(); z++)
        for (uns y = 0; y < belso_cellak.y_size(); y++) {
            belso_cellak.getref(0, y, z).face_adat[WEST].junction_index = face_adat[WEST].junction_index;
            belso_cellak.getref(belso_cellak.x_size() - 1, y, z).face_adat[EAST].junction_index = face_adat[EAST].junction_index;
        }
    for (uns z = 0; z < belso_cellak.z_size(); z++)
        for (uns x = 0; x < belso_cellak.x_size(); x++) {
            belso_cellak.getref(x, 0, z).face_adat[SOUTH].junction_index = face_adat[SOUTH].junction_index;
            belso_cellak.getref(x, belso_cellak.y_size() - 1, z).face_adat[NORTH].junction_index = face_adat[NORTH].junction_index;
        }
    for (uns y = 0; y < belso_cellak.y_size(); y++)
        for (uns x = 0; x < belso_cellak.x_size(); x++) {
            belso_cellak.getref(x, y, 0).face_adat[BOTTOM].junction_index = face_adat[BOTTOM].junction_index;
            belso_cellak.getref(x, y, belso_cellak.z_size() - 1).face_adat[TOP].junction_index = face_adat[TOP].junction_index;
        }
}

//***********************************************************************
void t_modell_cella::set_egy_belso_cella(const t_modell_cella & tulaj) {
// minden belsõ cellát ugyanarra állít, azaz nem figyeli, hogy peremen van-e
//***********************************************************************
    uns x_mul = tulaj.belso_cellak.x_size();
    uns y_mul = tulaj.belso_cellak.y_size();
    uns z_mul = tulaj.belso_cellak.z_size();
    is_cella = tulaj.is_cella;
    is_el = tulaj.is_el;
    is_th = tulaj.is_th;
    is_nonlin_el = tulaj.is_nonlin_el;
    is_nonlin_th = tulaj.is_nonlin_th;
    color_index = tulaj.color_index;
    anyag_index = tulaj.anyag_index;
    pmat = tulaj.pmat;
    V = tulaj.V / (x_mul*y_mul*z_mul);
    for (uns i = 1; i < BASIC_SIDES; i++) {
        face_adat[i].kulso_el_db = is_el ? 1 : 0;
        face_adat[i].kulso_th_db = is_th ? 1 : 0;
        face_adat[i].is_el_perem = face_adat[i].is_th_perem = false;
        face_adat[i].el_perem_index = face_adat[i].th_perem_index = 0;
        face_adat[i].th_perem_x = tulaj.face_adat[i].th_perem_x;
        face_adat[i].th_perem_y = tulaj.face_adat[i].th_perem_y;
        face_adat[i].th_perem_c = tulaj.face_adat[i].th_perem_c;
        face_adat[i].anyag_index = tulaj.face_adat[i].anyag_index;
    }
    face_adat[WEST].L = tulaj.face_adat[WEST].L / x_mul;
    face_adat[WEST].A = tulaj.face_adat[WEST].A / (y_mul*z_mul);
    face_adat[EAST].L = tulaj.face_adat[EAST].L / x_mul;
    face_adat[EAST].A = tulaj.face_adat[EAST].A / (y_mul*z_mul);
    face_adat[SOUTH].L = tulaj.face_adat[SOUTH].L / y_mul;
    face_adat[SOUTH].A = tulaj.face_adat[SOUTH].A / (x_mul*z_mul);
    face_adat[NORTH].L = tulaj.face_adat[NORTH].L / y_mul;
    face_adat[NORTH].A = tulaj.face_adat[NORTH].A / (x_mul*z_mul);
    face_adat[BOTTOM].L = tulaj.face_adat[BOTTOM].L / z_mul;
    face_adat[BOTTOM].A = tulaj.face_adat[BOTTOM].A / (x_mul*y_mul);
    face_adat[TOP].L = tulaj.face_adat[TOP].L / z_mul;
    face_adat[TOP].A = tulaj.face_adat[TOP].A / (x_mul*y_mul);
}


//***********************************************************************
void t_modell_cella::write_cella(FILE * fp, FILE *descfp, meret_tomb_tipus & meret_tomb, uns cx, uns cy, uns cz) {
//***********************************************************************
    if (!is_cella)
        return;
    if (belso_cellak.size() > 0) {
        for (uns z = 0; z < belso_cellak.z_size(); z++)
            for (uns y = 0; y < belso_cellak.y_size(); y++)
                for (uns x = 0;x < belso_cellak.x_size(); x++)
                    belso_cellak.getref(x, y, z).write_cella(fp, descfp, meret_tomb, cx, cy, cz);
    }
    else {
        fprintf(fp, "BL%uC%uM%uV%zu#(%u,%u,%u)\n", cella_index, color_index, anyag_index, meret_tomb.get_index(V), cx, cy, cz);
        uns cella_tipus = 0; // cellatipus_azonosito(); // egyelõre maradjunk a faces celláknál
        if (cella_tipus == 1) {
            fprintf(fp, "T1\n");
        }
        else {
            fprintf(fp, "F");
            if (is_el)fprintf(fp, "E");
            if (is_th)fprintf(fp, "T");
            fprintf(fp, "\n");
        }
        uns db_face = write_faces(fp, meret_tomb, cella_tipus);
        if (fenyforras_cella_index.size() > 0) {
            for (uns i = 0; i < fenyforras_cella_index.size(); i++) {
                fprintf(fp, "L%u", fenyforras_cella_index[i]);
                pmat->light_vertical_light_conversion.write_normal(fp, 0);
            }
        }
        fprintf(fp, "EL%u\n", cella_index);
        if (is_el || (pmat != nullptr && pmat->is_fenypor))fprintf(descfp, is_el ? "E" : "P");
        if (is_th)fprintf(descfp, "T");
        fprintf(descfp, "\t%u\t%u\t%u\t%u\t%u\n", cella_index, cx, cy, cz, db_face);
    }
}


//***********************************************************************
uns t_modell_cella::write_faces(FILE * fp, meret_tomb_tipus & meret_tomb, uns cella_tipus){
// Elemi cella hívhatja
//***********************************************************************

    // face-ek számolása

    uns db_face = 0;

    // external: i=0 esetén 

    for (uns i = 1; i < BASIC_SIDES; i++) {
        if (face_adat[i].belso_facek.size() > 0) {
            for (uns y = 0; y < face_adat[i].belso_facek.y_size(); y++)
                for (uns x = 0; x < face_adat[i].belso_facek.x_size(); x++) {
                    if (is_el) db_face++;
                    if (is_th) db_face++;
                }
        }
        else {
            if (is_el) db_face++;
            if (is_th) db_face++;
        }
    }

    // External face-ek: TODO
    // az external a face_adat[0] lenne?

    if (cella_tipus != 1)
        fprintf(fp, "NF%u\n", db_face);

    //uns akt_extern_index = 1;
    //uns akt_csatlakozo_index = 1;
    //uns akt_perem_index = 1;
    uns akt_face_index = 1;

    // externalok:
    // face_adat[0].write_face(fp, is_el, is_th, akt_extern_index, akt_perem_index);

    for (uns i = 1; i < BASIC_SIDES; i++) {
        face_adat[i].write_face(fp, is_el, is_th, akt_face_index, meret_tomb, anyag_index, cella_tipus);
    }
    return db_face;
}


//***********************************************************************
uns t_modell_cella::cellatipus_azonosito() {
// 1: pontosan 6 normál csatlakozó termikus face-e van, semmi egyéb
//***********************************************************************
    if (belso_cellak.size() > 0)
        return 0;
    uns db_1 = 0, db_egyeb = 0;
    for (uns i = 1; i < BASIC_SIDES; i++) {
        uns akt_face_tipus = face_adat[i].facetipus_azonosito(is_el, is_th, anyag_index);
        if (akt_face_tipus == 0)
            return 0;
        if (akt_face_tipus == 1)
            db_1++;
        else
            db_egyeb++;
    }
    
    if (db_egyeb > 0)
        return 0;

    if (db_1 == 6)
        return 1;

    return 0;
}


//***********************************************************************
void t_modell_cella::set_face_indexek(){
//***********************************************************************
    if (!is_cella)
        return;
    if (belso_cellak.size() > 0) {
        for (uns z = 0; z < belso_cellak.z_size(); z++)
            for (uns y = 0; y < belso_cellak.y_size(); y++)
                for (uns x = 0;x < belso_cellak.x_size(); x++)
                    belso_cellak.getref(x, y, z).set_face_indexek();
    }

    uns akt_face_index = 1;

    for (uns i = 1; i < BASIC_SIDES; i++) {
        face_adat[i].set_face_indexek(is_el, is_th, akt_face_index);
    }
    if (junction_bottom_face != 0) {
        junction_bottom_face = face_adat[BOTTOM].face_index_el;
        if (junction_bottom_face == 0)
            throw hiba("set_face_indexek", "junction_bottom_face == 0");
    }
    if (junction_top_face != 0) {
        junction_top_face = face_adat[TOP].face_index_el;
        if (junction_top_face == 0)
            throw hiba("set_face_indexek", "junction_top_face == 0");
    }
}


//***********************************************************************
const csomag & simulation::get_inner_perem(uns x, uns y, uns z)const {
//***********************************************************************
    uns szin = pmodel->tbmp[z].getpixel_also(x, y);
    const color &aktszin = pmodel->tcolor[szin];
    if (!aktszin.is)
        throw hiba("simulation::get_inner_perem", "color not exists (%u)", szin);
    if (aktszin.tipus == SzinBoundary)
        return tinner[innerindex[szin]];
    else if (aktszin.tipus == SzinUchannel) {
        uns tucha_index = pmodel->tbmp[z].getpixel_felso(x, y);
        if (tucha_perem[szin].size() <= tucha_index)
            throw hiba("simulation::get_inner_perem", "program error: tucha_perem[%u].size()<=tucha_index (%u<=%u)", szin, tucha_perem[szin].size(), tucha_index);
        return tucha_perem[szin][tucha_index];
    }
    else
        throw hiba("simulation::get_inner_perem", "program error: impossible boundary type", szin);
}


//***********************************************************************
const boundary * simulation::get_perem(uns x, uns y, uns z, Oldal oldal, bool is_el)const {
//***********************************************************************
    const boundary * b = NULL;
    if (is_el) {
        switch (oldal) {
        case WEST:   b = x == 0 ? &normalperem.el[WEST] : &get_inner_perem(x - 1, y, z).el[EAST];   break;
        case EAST:   b = x == pmodel->x_res - 1 ? &normalperem.el[EAST] : &get_inner_perem(x + 1, y, z).el[WEST];   break;
        case SOUTH:  b = y == 0 ? &normalperem.el[SOUTH] : &get_inner_perem(x, y - 1, z).el[NORTH];  break;
        case NORTH:  b = y == pmodel->y_res - 1 ? &normalperem.el[NORTH] : &get_inner_perem(x, y + 1, z).el[SOUTH];  break;
        case BOTTOM: b = z == 0 ? &normalperem.el[BOTTOM] : &get_inner_perem(x, y, z - 1).el[TOP];    break;
        case TOP:    b = z == pmodel->z_res - 1 ? &normalperem.el[TOP] : &get_inner_perem(x, y, z + 1).el[BOTTOM]; break;
        default: throw hiba("simulation::get_perem_u", "Unknown side.");
        }
    }
    else {
        switch (oldal) {
        case WEST:   b = x == 0 ? &normalperem.th[WEST] : &get_inner_perem(x - 1, y, z).th[EAST];   break;
        case EAST:   b = x == pmodel->x_res - 1 ? &normalperem.th[EAST] : &get_inner_perem(x + 1, y, z).th[WEST];   break;
        case SOUTH:  b = y == 0 ? &normalperem.th[SOUTH] : &get_inner_perem(x, y - 1, z).th[NORTH];  break;
        case NORTH:  b = y == pmodel->y_res - 1 ? &normalperem.th[NORTH] : &get_inner_perem(x, y + 1, z).th[SOUTH];  break;
        case BOTTOM: b = z == 0 ? &normalperem.th[BOTTOM] : &get_inner_perem(x, y, z - 1).th[TOP];    break;
        case TOP:    b = z == pmodel->z_res - 1 ? &normalperem.th[TOP] : &get_inner_perem(x, y, z + 1).th[BOTTOM]; break;
        default: throw hiba("simulation::get_perem_u", "Unknown side.");
        }
    }
    return b;
}


//***********************************************************************
enum sorazon{az_unknown,az_nev,az_dim,az_xres,az_yres,az_zres,az_coor,
             az_rmin,az_rmax,az_xpit,az_ypit,az_zpit,az_highres,az_bmp,
             az_mater,az_mat,az_color,az_col,az_semicon,az_semi,az_coupled,
             az_end,az_field,az_calc_mod,az_lin,az_nosemi,az_noseebeck,az_nopeltier,
             az_peltier_center,az_nothomson,az_no_joule,az_noaccelerator,
             az_el_nonlin_subiter,az_comp,az_always_strassen_mul,az_cpu,
             az_hdd,az_bound,az_optimize,az_bou,az_excita,az_fimtxt,
             az_commas,az_nofim,az_ambem,az_excit,az_ambient,az_powermap,
             az_uchannel,az_no_plus_step_data, az_general_reflectivity,
             az_probe,az_pro,az_analy,az_convection,az_conv,az_extnods,
             az_footpr,az_admmx,az_inhom,az_cm_temp,az_always_quad,
             az_timestep,az_timelimit,az_change_time,az_electrical,
             az_thermal,az_ndc_miniter,az_FIM_res,az_FIM_diff,az_half,
             az_auto_transi_steps_V,az_auto_transi_steps_T,az_auto_transi_steps};
//***********************************************************************


//***********************************************************************
sorazon azonosit(const PLString & szo){
// szo csupa kisbetûs kell legyen
//***********************************************************************
    if(szo=="name")return az_nev;
    if(szo=="dimension")return az_dim;
    if(szo=="x-resolution")return az_xres;
    if(szo=="y-resolution")return az_yres;
    if(szo=="z-resolution")return az_zres;
    if(szo=="coordinates")return az_coor;
    if(szo=="r-min")return az_rmin;
    if(szo=="r-max")return az_rmax;
    if(szo=="x-pitch")return az_xpit;
    if(szo=="y-pitch")return az_ypit;
    if(szo=="z-pitch")return az_zpit;
    if (szo == "high-res-region")return az_highres;
    if (szo == "bitmap")return az_bmp;
    if(szo=="ambient-emissivity")return az_ambem;
    if (szo == "general-reflectivity")return az_general_reflectivity;
    if (szo == "half")return az_half;
    if(szo=="material")return az_mater;
    if(szo=="mat")return az_mat;
    if(szo=="color")return az_color;
    if(szo=="col")return az_col;
    if(szo=="semiconductor")return az_semicon;
    if(szo=="semi")return az_semi;
    if(szo=="coupled_model")return az_coupled;
    if(szo=="end")return az_end;
    if(szo=="field")return az_field;
    if (szo == "calc_mode")return az_calc_mod;
    if(szo=="linear")return az_lin;
    if(szo=="no_semicond")return az_nosemi;
    if(szo=="no_seebeck")return az_noseebeck;
    if(szo=="no_peltier")return az_nopeltier;
    if(szo=="peltier_center")return az_peltier_center;
    if(szo=="no_thomson")return az_nothomson;
    if(szo=="no_joule")return az_no_joule;
    if(szo=="no_accelerator")return az_noaccelerator;
    if(szo=="cm_temperature")return az_cm_temp;
    if (szo == "el_nonlin_subiter")return az_el_nonlin_subiter;
    if (szo == "ndc_miniter")return az_ndc_miniter;
    if(szo=="computation")return az_comp;
    if(szo=="always_quad")return az_always_quad;
    if(szo=="always_strassen_mul")return az_always_strassen_mul;
    if(szo=="cpu-threads")return az_cpu;
    if(szo=="hdd-cache")return az_hdd;
    if(szo=="optimize")return az_optimize;
    if(szo=="convection")return az_convection;
    if(szo=="conv")return az_conv;
    if(szo=="boundary")return az_bound;
    if(szo=="bou")return az_bou;
    if(szo=="excitation")return az_excita;
    if(szo=="excit")return az_excit;
    if(szo=="ambient")return az_ambient;
    if(szo=="ambient-temperature")return az_ambient;
    if (szo == "powermap")return az_powermap;
    if (szo == "uchannel")return az_uchannel;
    if(szo == "probe")return az_probe;
    if(szo == "map_to_txt")return az_fimtxt;
    if(szo=="no_images")return az_nofim;
    if (szo == "use_commas")return az_commas;
    if (szo == "fim_res")return az_FIM_res;
    if (szo == "fim_diff")return az_FIM_diff;
    if(szo=="pro")return az_pro;
    if (szo == "no_plus_step_data")return az_no_plus_step_data;
    if (szo == "auto_transi_steps_v")return az_auto_transi_steps_V;
    if (szo == "auto_transi_steps_t")return az_auto_transi_steps_T;
    if (szo == "auto_transi_steps")return az_auto_transi_steps;
    if(szo=="analysis")return az_analy;
    if(szo=="external_nodes")return az_extnods;
    if(szo=="footprint")return az_footpr;
    if(szo=="adm")return az_admmx;
    if(szo=="inhom")return az_inhom;
    if(szo=="timestep")return az_timestep;
    if(szo=="timelimit")return az_timelimit;
    if(szo=="change_time")return az_change_time;
    if(szo=="electrical")return az_electrical;
    if(szo=="thermal")return az_thermal;
    if(szo!="")printf("# unknown parameter is ignored: %s\n",szo.c_str());
    return az_unknown;
}


//***********************************************************************
enum matazon{maz_unknown, maz_phase_change_energy,maz_heatcond,maz_heatres,maz_heatcap,maz_conduct,
             maz_resistivity,maz_capac,maz_seebeck,maz_disscoeff,
             maz_emissivity, maz_light_blue_absorption_coeff, maz_light_conversion_efficiency, maz_light_yellow_absorption_coeff,
             maz_vertical_light_conversion, maz_direction_mode, maz_cut_level, maz_d_light_powder, maz_ray_per_cell_dir,
             maz_output_side, maz_re_conversion_efficiency, maz_reflectivity, maz_yellow_correction
};
//***********************************************************************


//***********************************************************************
matazon azon_mat(const PLString & szo){
// szo csupa kisbetûs kell legyen
//***********************************************************************
    if (szo == "phase_change_energy")return maz_phase_change_energy;
    if(szo=="heatcond")return maz_heatcond;
    if(szo=="heatres")return maz_heatres;
    if(szo=="heatcap")return maz_heatcap;
    if(szo=="conduct")return maz_conduct;
    if(szo=="resistivity")return maz_resistivity;
    if(szo=="capacity")return maz_capac;
    if(szo=="seebeck")return maz_seebeck;
    if(szo=="dissip_coeff")return maz_disscoeff;
    if (szo == "emissivity")return maz_emissivity;
    if (szo == "blue_absorption_coeff")return maz_light_blue_absorption_coeff;
    if (szo == "conversion_efficiency")return maz_light_conversion_efficiency;
    if (szo == "yellow_absorption_coeff")return maz_light_yellow_absorption_coeff;
    if (szo == "yellow_correction")return maz_yellow_correction;
    if (szo == "re_conversion_efficiency")return maz_re_conversion_efficiency;
    if (szo == "vertical_light_conversion")return maz_vertical_light_conversion;
    if (szo == "direction_mode")return maz_direction_mode;
    if (szo == "d_light_powder")return maz_d_light_powder;
    if (szo == "ray_per_cell_dir")return maz_ray_per_cell_dir;
    if (szo == "output_side")return maz_output_side;
    if (szo == "cut_level")return maz_cut_level;
    if (szo == "reflectivity")return maz_reflectivity;
    if(szo!="")printf("# unknown meterial parameter is ignored: %s\n",szo.c_str());
    return maz_unknown;
}


//***********************************************************************
enum convazon{caz_unknown,caz_rad,caz_vertical,caz_upper,caz_lower,
              caz_axis,caz_angle,caz_edge};
//***********************************************************************


//***********************************************************************
convazon azon_conv(const PLString & szo){
// szo csupa kisbetûs kell legyen
//***********************************************************************
    if(szo=="radiation")return caz_rad;
    if(szo=="free-vertical")return caz_vertical;
    if(szo=="free-upper")return caz_upper;
    if(szo=="free-lower")return caz_lower;
    if(szo=="axis")return caz_axis;
    if(szo=="angle")return caz_angle;
    if(szo=="edge")return caz_edge;
    if(szo!="")printf("# unknown convection parameter is ignored: %s\n",szo.c_str());
    return caz_unknown;
}


//***********************************************************************
bool pitchconverter(const tomb<PLString> & be,tomb<z_a_hengerhez> & ki,tomb<dbl> & hossz){
// be[1]-tõl megy, nem 0-tól!
//***********************************************************************
//    uns j=0;
    ki.clear();
    for(uns i=1;i<be.size();i++){
        int n=be[i].find('*');
        z_a_hengerhez d;
        if(!be[i].todouble(d.ertek))return false;
        if(n!=npos){
            if(!be[i].toint(n,n+1))return false;
            ki.add(d,n);
        }
        else ki.add(d);
    }
    hossz.clear();
    hossz.add(ki[0].ertek*0.5);
    hossz.add(ki[0].ertek);
    for(u32 i=1;i<ki.size();i++){
        hossz.add(hossz.getLast()+ki[i].ertek*0.5);
        hossz.add(hossz.getLast()+ki[i].ertek*0.5); // ez már az elõzõ sorban kapotthoz adja
    }
    return true;
}


//***********************************************************************
bool pitchconverter(const tomb<PLString> & be, tomb<dbl> & ki) {
// be[1]-tõl megy, nem 0-tól!
//***********************************************************************
//    uns j=0;
    ki.clear();
    for (uns i = 1;i<be.size();i++) {
        int n = be[i].find('*');
        dbl d;
        if (!be[i].todouble(d))return false;
        if (n != npos) {
            if (!be[i].toint(n, n + 1))return false;
            ki.add(d, n);
        }
        else ki.add(d);
    }
    return true;
}


//***********************************************************************
double interpolalo(tomb2d<dbl> &t, uns x_res, uns y_res, dbl x, dbl y){
//***********************************************************************
    double uj_x = x * ( t.x_size() - 1 ) / ( x_res - 1) ;
    double uj_y = y * ( t.y_size() - 1 ) / ( y_res - 1) ;
    uns also_x = (uns)uj_x;
    uns also_y = (uns)uj_y;
    uns felso_x = ( also_x >= t.x_size() - 1 ) ? also_x : also_x + 1;
    uns felso_y = ( also_y >= t.y_size() - 1 ) ? also_y : also_y + 1;
    double arany_x = uj_x - also_x;
    double arany_y = uj_y - also_y;
    return (1.0 - arany_x) * (1.0 - arany_y) * t.get( also_x,  also_y  )+
           (      arany_x) * (1.0 - arany_y) * t.get( felso_x, also_y  )+
           (1.0 - arany_x) * (      arany_y) * t.get( also_x,  felso_y )+
           (      arany_x) * (      arany_y) * t.get( felso_x, felso_y );
}


//***********************************************************************
bool konv_map_konverter(const PLString fajlnev, tomb2d<dbl> & cel, uns x_res, uns y_res){
// Az x_res és y_res a cél kívánt mérete, ekkorára skálázza fel a beolvasott mapot
//***********************************************************************
    const char * fvnev="konv_map_konverter()";
    srfajl fajl;
    fajl.open(fajlnev);
    if(fajl.lines().size()<1||fajl.lines()[0][0].LowCase()!="vsun3-convection-map")
        throw hiba(fvnev,"first line is not vsun3-convection-map in %s",fajlnev.c_str());
    if(fajl.lines()[1].size()<2)throw hiba(fvnev,"incomplete line 2 in file: %s",fajlnev.c_str());
    if(fajl.lines()[2].size()<2)throw hiba(fvnev,"incomplete line 3 in file: %s",fajlnev.c_str());
    if(fajl.lines()[1][0].LowCase()!="x-size")throw hiba(fvnev,"line 2 must be \"X-SIZE=<number>\" in %s",fajlnev.c_str());
    if(fajl.lines()[2][0].LowCase()!="y-size")throw hiba(fvnev,"line 3 must be \"Y-SIZE=<number>\" in %s",fajlnev.c_str());
    uns x,y;
    if(!fajl.lines()[1][1].tounsigned(x))throw hiba(fvnev,"X-SIZE is not a number (%s) in %s",fajl.lines()[1][1].c_str(),fajlnev.c_str());
    if(!fajl.lines()[2][1].tounsigned(y))throw hiba(fvnev,"Y-SIZE is not a number (%s) in %s",fajl.lines()[2][1].c_str(),fajlnev.c_str());

    cu32 n=x*y;
    tomb<dbl> t;
    t.clear();
    for(uns i=3; t.size()<n; i++){
        if(fajl.lines().size()<=i)
            throw hiba(fvnev,"unexpected end of file in %s",fajlnev.c_str());
        if(fajl.lines()[i][0].LowCase()=="end")
            throw hiba(fvnev,"less convection value than X-SIZE*Y-SIZE (%u<%u) in %s",t.size(),n,fajlnev.c_str());
        for(uns j=0; j<fajl.lines()[i].size() && t.size()<n; j++){
            int m=fajl.lines()[i][j].find('*');
            dbl d;
            if(!fajl.lines()[i][j].todouble(d))
                throw hiba("konv_map_konverter","cannot be converted to double (%s) in %s",fajl.lines()[i][j].c_str(),fajlnev.c_str());
            if(m!=npos){
                if(!fajl.lines()[i][j].toint(m,m+1))
                    throw hiba("konv_map_konverter","cannot be converted to int (%s) in %s",fajl.lines()[i][j].c_str(),fajlnev.c_str());
                t.add(d,m);
            }
            else t.add(d);
        }
    }

    tomb2d<dbl> t_2d;
    t_2d.free();
    t_2d.resize(x,y);
    uns k=0;
    for(uns i=0; i<y; i++)
        for(uns j=0; j<x; j++,k++)
            t_2d.getref(j,i)=t[k];

    cel.free();
    cel.resize(x_res,y_res);
    for(uns i=0; i<y_res; i++)
        for(uns j=0; j<x_res; j++){
            double d=interpolalo(t_2d,x_res,y_res,j,y_res-i-1)*2.1839; // !!!!!! Fejjel lefelé és szorozva egy konstanssal! Ha rendes bemenetet kap, ezt kivenni!
            cel.getref(j,i)=d;
        }

    return true;
}


//***********************************************************************
uns function_azonosito(const PLString &nev, monlintipus &tipus) {
// a név alapján beállítja a típust és visszaadja a beolvasandó valós 
// értékek számát, melyek a gg-be kerülnek
//***********************************************************************
    if (nev == "mizs1") {
        tipus = nlt_mizs1;
        return 7;
    }
    if (nev == "linear") {
        tipus = nlt_linearis;
        return 2;
    }
    else
        throw hiba("function_azonosito", "unknown material function name (%s)", nev.c_str());
    return 0;
}


//***********************************************************************
bool R_converter(const tomb<PLString> & sor,u32 startindex,vezetes & dest,const material & m,bool is_diode=false){
//***********************************************************************
    const PLString s=sor[startindex].LowCase();
    if (s == "function") {
        cuns db = function_azonosito(sor[startindex + 1], dest.tipus);
        if (db > GGSIZE)
            throw("R_converter", "program error: function parameter number>GGSIZE, GGSIZE increase is required");
        for (uns i = 0;i < db; i++)
            if (!sor[startindex + 2 + i].todouble(dest.gg[i]))
                return false;
        if (dest.tipus == nlt_mizs1) {
            dest.g[0] = dest.g[1] = dest.g[2] = fn_mizs1(25.0, dest.gg[0], dest.gg[1], dest.gg[2], dest.gg[3], dest.gg[4], dest.gg[5], dest.gg[6]);
        }
        if (dest.tipus == nlt_linearis) {
            dest.g[0] = dest.g[1] = dest.g[2] = dest.gg[0] + 25 * dest.gg[1];
        }
        dest.T = nulla;
    }
    else if(s[0]!='['){
        int y=s.find('|');
        int z=(y==npos)?npos:s.find('|',y+1);
        int g=s.find('#');
        int g2=(g==npos)?npos:s.find('|',g+1);
        int g3=(g2==npos)?npos:s.find('|',g2+1);
        if(g!=npos){
            if(y!=npos&&y>g)y=npos;
            if(z!=npos&&z>g)z=npos;
        }

        if(!s.todouble(dest.g[0]))return false;
        dest.g[2]=dest.g[1]=dest.g[0];
        dest.T=nulla;

        if(y!=npos)if(!s.todouble(dest.g[1],y+1))return false;
        if(z!=npos)if(!s.todouble(dest.g[2],z+1))return false;
        if(g!=npos){if(!s.todouble(dest.gg[0],g+1))return false;}
        else dest.gg[0]=nulla;
        if(g2!=npos)if(!s.todouble(dest.gg[1],g2+1))return false;
        if(g3!=npos)if(!s.todouble(dest.gg[2],g3+1))return false;
    
        if(g==npos)dest.tipus=nlt_lin;
        else if(g2==npos)dest.tipus=nlt_exp;
        else if(g3==npos)dest.tipus=nlt_diode;
        else dest.tipus = nlt_quadratic;

        if(y==npos)dest.semitip=nlt_exp;
        else if(z==npos){
            dest.semitip=nlt_diode;
            dest.g[2]=0.1;
        }
        else{
            if(is_diode){
                dest.semitip=nlt_diode;
            }
            else dest.semitip=nlt_quadratic;
        }
    }
    else{ // poligonnal megadott paraméter esetén
        dest.tipus = nlt_szakaszok;
        while(startindex<sor.size()-1){
            if(sor[startindex][0]!='[')return false;
            u32 poliindex=dest.szakaszok.size();
            dest.szakaszok.resize(dest.szakaszok.size()+1);

            // '[' átugrása és T beolvasása

            bool normal=true;
            if(m.is_his){
                bool f=true;
                if(sor[startindex].LowCase()=="[f" || sor[startindex+1].LowCase()=="f")
                    normal=false;
                else if(sor[startindex].LowCase()=="[h" || sor[startindex+1].LowCase()=="h"){
                    normal=false;
                    f=false;
                }
                if(!normal){
                    dest.is_his = true;
                    if(sor[startindex]=="[")
                        startindex++;
                    startindex++;
                    if(!sor[startindex].todouble(dest.his_value_1))return false;
                    if(sor[startindex][sor[startindex].Length()-1]!=']'&& sor[startindex+1]!="]"){
                        startindex++;
                        if(!sor[startindex].todouble(dest.his_value_2))return false;
                    }
                    startindex++;
                    if(startindex<sor.size() && sor[startindex]=="]")startindex++; // ha nem volt szóköz, akkor automatikusan átugrotta
                    if(dest.his_value_2==nulla)
                        dest.szakaszok.resize(dest.szakaszok.size()-1);
                    else{
                        dest.szakaszok[poliindex].T = m.his_T_min-m.his_T_width_fele;
                        dest.szakaszok[poliindex].G[0] = dest.szakaszok[poliindex].G[1] =
                            dest.szakaszok[poliindex].G[2] = dest.his_value_1;
                        dest.szakaszok[poliindex].is_F = true;
                        dest.szakaszok.resize(dest.szakaszok.size()+1);
                        poliindex++;
                        dest.szakaszok[poliindex].T = m.his_T_max+m.his_T_width_fele;
                        dest.szakaszok[poliindex].G[0] = dest.szakaszok[poliindex].G[1] =
                            dest.szakaszok[poliindex].G[2] = dest.his_value_2;
                        dest.szakaszok[poliindex].is_F = true;
                    }
                }
            }
            if(normal){
                if(sor[startindex]=="["){// szóköz volt a [ után
                    startindex++;
                    if(!sor[startindex].todouble(dest.szakaszok[poliindex].T))return false;
                }
                else{ // a [ közvetlenül a szám elõtt van
                    if(!sor[startindex].todouble(dest.szakaszok[poliindex].T,1))return false;
                }

                // G beolvasása és ] átugrása

                startindex++;
                if(startindex>=sor.size())return false;
                int y=sor[startindex].find('|');
                int z=(y==npos)?npos:sor[startindex].find('|',y+1);

                if(!sor[startindex].todouble(dest.szakaszok[poliindex].G[0]))return false;
                dest.szakaszok[poliindex].G[2] = dest.szakaszok[poliindex].G[1] = dest.szakaszok[poliindex].G[0];

                if(y!=npos)if(!sor[startindex].todouble(dest.szakaszok[poliindex].G[1],y+1))return false;
                if(z!=npos)if(!sor[startindex].todouble(dest.szakaszok[poliindex].G[2],z+1))return false;

                startindex++;
                if(startindex<sor.size() && sor[startindex]=="]")startindex++; // ha nem volt szóköz, akkor automatikusan átugrotta
            }
        }
    }

    return true;
}


//***********************************************************************
bool Erno_converter(const PLString & fajlnev,vezetes & dest){
// a fájlban erno modell van
//***********************************************************************
    const char * fvnev="Erno_converter()";
    dbl A=1.0;
    uns index=0;
    srfajl fajl;
    fajl.open(fajlnev);
    if(fajl.lines().size()<1||fajl.lines()[0][0].LowCase()!="vsun3-ernomodel")
        throw hiba(fvnev,"first line is not vsun3-ernomodel in %s",fajlnev.c_str());
    double d,m[3][3],b[3][3];
    if(fajl.lines().size()<5)throw hiba(fvnev,"incomplete file: %s",fajlnev.c_str());
    if(fajl.lines()[1].size()<2)throw hiba(fvnev,"incomplete line 2 in file: %s",fajlnev.c_str());
    if(fajl.lines()[2].size()<1)throw hiba(fvnev,"incomplete line 3 in file: %s",fajlnev.c_str());
    if(fajl.lines()[2][0].LowCase()=="specific"){ // fajlagos (1 m2-re vonatkozó)
        dest.specific=true;
        index++;
    }
    if(fajl.lines()[2+index][0].LowCase()=="a"){ // van arányszorzó
        if(!fajl.lines()[2+index][1].todouble(d))throw hiba(fvnev,"line 3 must be \"A=<number>\", not a number in %s",fajlnev.c_str());
        A=d;
        index++;
    }
    if(fajl.lines()[2+index].size()<10)throw hiba(fvnev,"incomplete line %u in file: %s",3+index,fajlnev.c_str());
    if(fajl.lines()[3+index].size()<10)throw hiba(fvnev,"incomplete line %u in file: %s",4+index,fajlnev.c_str());
    if(fajl.lines()[1][0].LowCase()!="d")throw hiba(fvnev,"line 2 must be \"d=<number>\" in %s",fajlnev.c_str());
    if(!fajl.lines()[1][1].todouble(d))throw hiba(fvnev,"line 2 must be \"d=<number>\", not a number in %s",fajlnev.c_str());
    if(fajl.lines()[2+index][0].LowCase()!="b")throw hiba(fvnev,"line %u must be \"B=<9 number>\" in %s",3+index,fajlnev.c_str());
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)
        if(!fajl.lines()[2+index][i*3+j+1].todouble(b[i][j]))throw hiba(fvnev,"in line %u not a number (%s) in %s",3+index,fajl.lines()[2+index][i*3+j].c_str(),fajlnev.c_str());
    if(fajl.lines()[3+index][0].LowCase()!="m")throw hiba(fvnev,"line %u must be \"M=<9 number>\" in %s",4+index,fajlnev.c_str());
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)
        if(!fajl.lines()[3+index][i*3+j+1].todouble(m[i][j]))throw hiba(fvnev,"in line %u not a number (%s) in %s",4+index,fajl.lines()[2][i*3+j].c_str(),fajlnev.c_str());

    for(int i=0;i<3;i++)dest.g[i] =(b[i][0]+b[i][1]*d+b[i][2]*d*d)*A;
    for(int i=0;i<3;i++)dest.gg[i]=m[i][0]+m[i][1]*d+m[i][2]*d*d;
    dest.tipus=dest.semitip=nlt_erno;
    return true;
}


//***********************************************************************
bool lsfit::read(const PLString & fajlnev) {
// a fájlban a legkisebb négyzetek módszerével illesztendõ adatok vannak
//***********************************************************************
    const char * fvnev = "lsfit::read()";
    uns index = 0;
    srfajl fajl;
    fajl.open(fajlnev);
    if (fajl.lines().size()<1 || fajl.lines()[0][0].LowCase() != "hex-least-square")
        throw hiba(fvnev, "first line is not hex-least-square in %s", fajlnev.c_str());

    // duo|trio

    index++;
    if (fajl.lines()[index][0].LowCase() == "duo")
        is_trio = false;
    else if (fajl.lines()[index][0].LowCase() == "trio")
        is_trio = true;
    else throw hiba(fvnev, "line %u must be duo|trio in file: %s, %s found", index, fajlnev.c_str(), fajl.lines()[index][0].c_str());
    index++;

    // multiplier (opc)

    multiplier = 1;
    if (fajl.lines()[index][0].LowCase() == "multiplier") {
        if(!fajl.lines()[index][1].todouble(multiplier))
            throw hiba(fvnev, "multiplier value is not a number in line %u in file: %s, %s found", index, fajlnev.c_str(), fajl.lines()[index][1].c_str());
        index++;
    }

    // start (opc)

    start = lsv_none;
    if (fajl.lines()[index][0].LowCase() == "start") {
        if (fajl.lines()[index][1].LowCase() == "none")
            start = lsv_none;
        else if (fajl.lines()[index][1].LowCase() == "lin")
            start = lsv_lin;
        else if (fajl.lines()[index][1].LowCase() == "strong")
            start = lsv_strong;
        else throw hiba(fvnev, "start must be none|lin|strong in file: %s, %s found", fajlnev.c_str(), fajl.lines()[index][0].c_str());
        index++;
    }

    // stop (opc)

    stop = lsv_none;
    if (fajl.lines()[index][0].LowCase() == "stop") {
        if (fajl.lines()[index][1].LowCase() == "none")
            stop = lsv_none;
        else if (fajl.lines()[index][1].LowCase() == "lin")
            stop = lsv_lin;
        else if (fajl.lines()[index][1].LowCase() == "strong")
            stop = lsv_strong;
        else throw hiba(fvnev, "stop must be none|lin|strong in file: %s, %s found", fajlnev.c_str(), fajl.lines()[index][0].c_str());
        index++;
    }

    // fit

    if (fajl.lines()[index][0].LowCase() == "fit") {
        if (fajl.lines()[index][1].LowCase() == "poli") {
            egyenlet = lse_polinom;
            if(!fajl.lines()[index][2].tounsigned(egyenlet_unspar_1) || !fajl.lines()[index][3].tounsigned(egyenlet_unspar_2))
                throw hiba(fvnev, "bad fit=poli parameters (%s,%s) in line %u in file: %s", 
                    fajl.lines()[index][2].c_str(), fajl.lines()[index][3].c_str(), index, fajlnev.c_str());
        }
        else throw hiba(fvnev, "unsupported fitting type (%s) in line %u in file: %s", fajl.lines()[index][1].c_str(), index, fajlnev.c_str());
    }
    else throw hiba(fvnev, "missing fit in line %u in file: %s, %s found", index, fajlnev.c_str(), fajl.lines()[index][0].c_str());

    // measures

    index++;
    if (fajl.lines()[index][0].LowCase() != "measures")
        throw hiba(fvnev, "missing measures in line %u in file: %s, %s found", index, fajlnev.c_str(), fajl.lines()[index][0].c_str());
    for (index++; index < fajl.lines().size() && fajl.lines()[index][0].LowCase() != "end"; index++) {
        lsfit_adat akt;
        if (is_trio) {
            if(fajl.lines()[index].size()<3)
                throw hiba(fvnev, "3 measures data needed, %u found in line %u in file: %s", fajl.lines()[index].size(), index, fajlnev.c_str());
            if(!fajl.lines()[index][0].todouble(akt.U) || !fajl.lines()[index][1].todouble(akt.T) || !fajl.lines()[index][2].todouble(akt.I))
                throw hiba(fvnev, "wrong measures data in line %u in file: %s, %s, %s, %s found", index, fajlnev.c_str(),
                    fajl.lines()[index][0].c_str(), fajl.lines()[index][1].c_str(), fajl.lines()[index][2].c_str());
            meresek.add(akt);
        }
        else { // duó
            if (fajl.lines()[index].size()<2)
                throw hiba(fvnev, "2 measures data needed, %u found in line %u in file: %s", fajl.lines()[index].size(), index, fajlnev.c_str());
            if (!fajl.lines()[index][0].todouble(akt.U) || !fajl.lines()[index][1].todouble(akt.I))
                throw hiba(fvnev, "wrong measures data in line %u in file: %s, %s, %s, %s found", index, fajlnev.c_str(),
                    fajl.lines()[index][0].c_str(), fajl.lines()[index][1].c_str());
            akt.T = 0;
            meresek.add(akt);
        }
    }
    if (fajl.lines()[index][0].LowCase() != "end")
        throw hiba(fvnev, "\"end\" is missing in line %u in file: %s, %s found", index, fajlnev.c_str(), fajl.lines()[index][0].c_str());
    return true;
}


//***********************************************************************
bool VR_converter(const PLString & s,vezetes & dest){
// vari-condhoz, ahol az elsõ a hõmérséklet, aztán a vezetés
//***********************************************************************
    int T=s.find('|');
    int y=(T==npos)?npos:s.find('|',T+1);
    int z=(y==npos)?npos:s.find('|',y+1);
    int g=s.find('#');
    int g2=(g==npos)?npos:s.find('|',g+1);
    int g3=(g2==npos)?npos:s.find('|',g2+1);
    if(g!=npos){
        if(T!=npos&&T>g)T=npos;
        if(y!=npos&&y>g)y=npos;
        if(z!=npos&&z>g)z=npos;
    }

    if(T==npos)return false;
    if(!s.todouble(dest.T))return false;
    if(!s.todouble(dest.g[0],T+1))return false;
    dest.g[2]=dest.g[1]=dest.g[0];

    if(y!=npos)if(!s.todouble(dest.g[1],y+1))return false;
    if(z!=npos)if(!s.todouble(dest.g[2],z+1))return false;
    if(g!=npos){if(!s.todouble(dest.gg[0],g+1))return false;}
    else dest.gg[0]=nulla;
    if(g2!=npos)if(!s.todouble(dest.gg[1],g2+1))return false;
    if(g3!=npos)if(!s.todouble(dest.gg[2],g3+1))return false;
    
    if(g==npos)dest.tipus=nlt_lin;
    else if(g2==npos)dest.tipus=nlt_exp;
    else if(g3==npos)dest.tipus=nlt_diode;
    else dest.tipus=nlt_quadratic;

    return true;
}


//***********************************************************************
Oldal oldalazonosíto(PLString azon){
//***********************************************************************
    Oldal oldal=EXTERNAL;
    if(azon=="west")oldal=WEST;
    else if(azon=="east")oldal=EAST;
    else if(azon=="south")oldal=SOUTH;
    else if(azon=="north")oldal=NORTH;
    else if(azon=="top")oldal=TOP;
    else if(azon=="bottom")oldal=BOTTOM;
    return oldal;
}


//***********************************************************************
void model::read(PLString path){
//***********************************************************************
    const char * fvnev="model::read()";
    logprint("Read %s",fileName.c_str());
    srfajl fajl;
    fajl.open(path+fileName);
    if(fajl.lines().size()<1||(fajl.lines()[0][0].LowCase()!="vsun3-model"
        && fajl.lines()[0][0].LowCase() != "vsun35-model"))
        throw hiba(fvnev,"first line is not vsun3-model or vsun35-model in %s (%s)",fileName.c_str(), fajl.lines()[0][0].LowCase().c_str());
    is_35 = fajl.lines()[0][0].LowCase() == "vsun35-model";
    uns i=0;
    sorazon az=az_unknown;

    // name

    try{
        do{
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_nev && az!=az_unknown)throw hiba(fvnev,"name is missing in line %u in %s",i,fileName.c_str());
        }while(az!=az_nev);
        name=fajl.lines()[i][1];

        // dimension

        do{
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_dim && az!=az_unknown)throw hiba(fvnev,"dimension is missing in line %u in %s",i,fileName.c_str());
        }while(az!=az_dim);
        if(!fajl.lines()[i][1].tounsigned(dim))throw hiba(fvnev,"dimension is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());

        // x-resolution

        do{
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_xres && az!=az_unknown)throw hiba(fvnev,"x-resolution is missing in line %u in %s",i,fileName.c_str());
        }while(az!=az_xres);
        if(!fajl.lines()[i][1].tounsigned(x_res))throw hiba(fvnev,"x-resolution is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());

        // y-resolution

        if(dim>1){
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_yres && az!=az_unknown)throw hiba(fvnev,"y-resolution is missing in line %u in %s",i,fileName.c_str());
            }while(az!=az_yres);
            if(!fajl.lines()[i][1].tounsigned(y_res))throw hiba(fvnev,"y-resolution is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());
        }
        else y_res=1;

        // z-resolution

        if(dim>2){
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_zres && az!=az_unknown)throw hiba(fvnev,"z-resolution is missing in line %u in %s",i,fileName.c_str());
            }while(az!=az_zres);
            if(!fajl.lines()[i][1].tounsigned(z_res))throw hiba(fvnev,"z-resolution is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());
        }
        else z_res=1;

        // coordinates

        do{
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az==az_yres||az==az_zres)continue;
            if(az!=az_coor && az!=az_unknown && az!=az_xpit)throw hiba(fvnev,"coordinates and x-pitch is missing in line %u in %s",i,fileName.c_str());
        }while(az!=az_coor && az!=az_xpit);
        kerekcoord=false;
        if(az==az_coor){
            PLString s=fajl.lines()[i][1].LowCase();
            if(s=="cartesian")kerekcoord=false;
            else if(s=="cylindrical" && dim==2)kerekcoord=true;
            else if(s=="sphere" && dim==1)kerekcoord=true;
            else throw hiba(fvnev,"unknown or not applicable coordinates (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());
        }

        // r-min, r-max

        if(kerekcoord){
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_rmin && az!=az_unknown)throw hiba(fvnev,"r-min is missing in line %u in %s",i,fileName.c_str());
            }while(az!=az_rmin);
            if(!fajl.lines()[i][1].todouble(r_min))throw hiba(fvnev,"r-min is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());
            
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_rmax && az!=az_unknown)throw hiba(fvnev,"r-max is missing in line %u in %s",i,fileName.c_str());
            }while(az!=az_rmax);
            if(!fajl.lines()[i][1].todouble(r_max))throw hiba(fvnev,"r-max is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());
            
            // x_pit és x_hossz feltöltése

            x_pit.clear();
            z_a_hengerhez d;
            d.ertek = (r_max - r_min) / x_res;
            x_pit.add(d,x_res);
            x_hossz.clear();
            x_hossz.add(x_pit[0].ertek*0.5 + r_min);
            x_hossz.add(x_pit[0].ertek + r_min);
            for(u32 i=1;i<x_pit.size();i++){
                x_hossz.add(x_hossz.getLast()+x_pit[i].ertek*0.5);
                x_hossz.add(x_hossz.getLast()+x_pit[i].ertek*0.5); // ez már az elõzõ sorban kapotthoz adja
            }
        }
        else{

        // x-pitch

            if(az!=az_xpit){
                do{
                    i++;
                    az=azonosit(fajl.lines()[i][0].LowCase());
                    if(az!=az_xpit && az!=az_unknown)throw hiba(fvnev,"x-pitch is missing in line %u in %s",i,fileName.c_str());
                }while(az!=az_xpit);
            }
            if(!pitchconverter(fajl.lines()[i],x_pit,x_hossz))throw hiba(fvnev,"x-pitch data corrupt in line %u in %s",i,fileName.c_str());
            if(x_pit.size()<x_res)throw hiba(fvnev,"fewer x-pitch data than resolution (%u<%u) in line %u in %s",x_pit.size(),x_res,i,fileName.c_str());
            if(x_pit.size()>x_res)logprint("warning: model::read() => more x-pitch data than resolution (%u>%u) in line %u in %s",x_pit.size(),x_res,i,fileName.c_str());
        }

        // y-pitch

        if(dim>1){
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_ypit && az!=az_unknown)throw hiba(fvnev,"y-pitch is missing in line %u in %s",i,fileName.c_str());
            }while(az!=az_ypit);
            if(!pitchconverter(fajl.lines()[i],y_pit,y_hossz))throw hiba(fvnev,"y-pitch data corrupt in line %u in %s",i,fileName.c_str());
            if(y_pit.size()<y_res)throw hiba(fvnev,"fewer y-pitch data than resolution (%u<%u) in line %u in %s",y_pit.size(),y_res,i,fileName.c_str());
            if(y_pit.size()>y_res)logprint("warning: model::read() => more y-pitch data than resolution (%u>%u) in line %u in %s",y_pit.size(),y_res,i,fileName.c_str());
        }

        // z-pitch

        if(dim>2 || !kerekcoord){
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_zpit && az!=az_unknown)throw hiba(fvnev,"z-pitch is missing in line %u in %s",i,fileName.c_str());
            }while(az!=az_zpit);
            if(!pitchconverter(fajl.lines()[i],z_pit,z_hossz))throw hiba(fvnev,"z-pitch data corrupt in line %u in %s",i,fileName.c_str());
            if(z_pit.size()<z_res)throw hiba(fvnev,"fewer z-pitch data than resolution (%u<%u) in line %u in %s",z_pit.size(),z_res,i,fileName.c_str());
            if(z_pit.size()>z_res)logprint("warning: model::read() => more z-pitch data than resolution (%u>%u) in line %u in %s",z_pit.size(),z_res,i,fileName.c_str());
        }
        else{
            z_pit.clear();
            z_a_hengerhez d;
            d.henger_e = true;
            z_pit.add(d);
            z_hossz.clear();
            z_hossz.add(0.0);
            z_hossz.add(0.0);
        }

		sprintf(hibaUzenet,"%ux%ux%u",x_res,y_res,z_res);

        // high-res-region

        do {
            i++;
            az = azonosit(fajl.lines()[i][0].LowCase());
            if (az == az_highres) {
                if (fajl.lines()[i].size() < 10)
                    throw hiba(fvnev, "too few parameters in line %u (high-res-region) in %s", i, fileName.c_str());
                righ_res_regions.incsize();
                high_res_region_struct & be = righ_res_regions.getLast();

                if (!fajl.lines()[i][1].tounsigned(be.x_res))
                    throw hiba(fvnev, "high-res-region x resolution is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
                if (!fajl.lines()[i][2].tounsigned(be.y_res))
                    throw hiba(fvnev, "high-res-region y resolution is not a number (%s) in %s", fajl.lines()[i][2].c_str(), fileName.c_str());
                if (!fajl.lines()[i][3].tounsigned(be.z_res))
                    throw hiba(fvnev, "high-res-region z resolution is not a number (%s) in %s", fajl.lines()[i][3].c_str(), fileName.c_str());

                if (!fajl.lines()[i][4].tounsigned(be.x1))
                    throw hiba(fvnev, "high-res-region x1 coord is not a number (%s) in %s", fajl.lines()[i][4].c_str(), fileName.c_str());
                if (!fajl.lines()[i][5].tounsigned(be.y1))
                    throw hiba(fvnev, "high-res-region y1 coord is not a number (%s) in %s", fajl.lines()[i][5].c_str(), fileName.c_str());
                if (!fajl.lines()[i][6].tounsigned(be.z1))
                    throw hiba(fvnev, "high-res-region z1 coord is not a number (%s) in %s", fajl.lines()[i][6].c_str(), fileName.c_str());

                if (!fajl.lines()[i][7].tounsigned(be.x2))
                    throw hiba(fvnev, "high-res-region x2 coord is not a number (%s) in %s", fajl.lines()[i][7].c_str(), fileName.c_str());
                if (!fajl.lines()[i][8].tounsigned(be.y2))
                    throw hiba(fvnev, "high-res-region y2 coord is not a number (%s) in %s", fajl.lines()[i][8].c_str(), fileName.c_str());
                if (!fajl.lines()[i][9].tounsigned(be.z2))
                    throw hiba(fvnev, "high-res-region z2 coord is not a number (%s) in %s", fajl.lines()[i][9].c_str(), fileName.c_str());

                if (be.x_res < 1 || be.x_res>100 || be.y_res < 1 || be.y_res>100 || be.z_res < 1 || be.z_res>100)
                    throw hiba(fvnev, "high-res-region: invalid resolution (%u,%u,%u) in %s", be.x_res, be.y_res, be.z_res, fileName.c_str());
                if (be.x1 > be.x2 || be.y1 > be.y2 || be.z1 > be.z2)
                    throw hiba(fvnev, "high-res-region: invalid region (%u,%u,%u) > (%u,%u,%u) in %s", be.x1, be.y1, be.z1, be.x2, be.y2, be.z2, fileName.c_str());
                if (be.x1 >= x_res || be.y1 >= y_res || be.z1 >= z_res)
                    throw hiba(fvnev, "high-res-region: invalid region (%u,%u,%u) >= (%u,%u,%u) in %s", be.x1, be.y1, be.z1, x_res, y_res, z_res, fileName.c_str());
                if (be.x2 >= x_res || be.y2 >= y_res || be.z2 >= z_res)
                    throw hiba(fvnev, "high-res-region: invalid region (%u,%u,%u) >= (%u,%u,%u) in %s", be.x1, be.y1, be.z1, x_res, y_res, z_res, fileName.c_str());
            }
        } while (az == az_highres||az == az_unknown); // minimum 1 material szükséges
        i--;

        // bitmap

        tbmp.clear();
        tbmp.resize(z_res);
        for(uns j=0;j<z_res;j++){
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_bmp && az!=az_unknown)throw hiba(fvnev,"bitmap is missing in line %u in %s",i,fileName.c_str());
            }while(az!=az_bmp);
            tbmp[j].load(path+fajl.lines()[i][1]);
            if(tbmp[j].getxsize()<x_res)
                throw hiba(fvnev,"bitmap x size < x resolution (%u<%u) in %s in %s",tbmp[j].getxsize(),x_res,fajl.lines()[i][1].c_str(),fileName.c_str());
            if(tbmp[j].getysize()<y_res)
                throw hiba(fvnev,"bitmap y size < y resolution (%u<%u) in %s in %s",tbmp[j].getysize(),y_res,fajl.lines()[i][1].c_str(),fileName.c_str());
        }

        // material vagy ambient-emissivity vagy general-reflectivity

        tmat.clear();
        material m_dummy;
        do{
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_mater && az!=az_unknown && az!=az_bmp && az!=az_ambem && az!=az_general_reflectivity && az!=az_half)
                throw hiba(fvnev,"ambient-emissivity or general-reflectivity or half or material is missing in line %u in %s",i,fileName.c_str());
        }while(az!=az_mater&&az!=az_ambem&&az!=az_general_reflectivity&&az!=az_half); // minimum 1 material szükséges

        if(az==az_ambem){
            if(!R_converter(fajl.lines()[i],1,amb_emiss,m_dummy))
                throw hiba(fvnev,"ambient-emissivity value wrong in line %u in %s",i,fileName.c_str());
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_mater && az!=az_unknown && az != az_general_reflectivity && az != az_half)
                    throw hiba(fvnev,"material or general-reflectivity or half is missing in line %u in %s",i,fileName.c_str());
            }while(az!=az_mater&&az != az_general_reflectivity && az != az_half); // minimum 1 material szükséges
        }
        
        general_reflectivity = 0;
        if (az == az_general_reflectivity) {
            if (!fajl.lines()[i][1].todouble(general_reflectivity))
                throw hiba(fvnev, "ambient-emissivity value wrong in line %u in %s", i, fileName.c_str());
            do {
                i++;
                az = azonosit(fajl.lines()[i][0].LowCase());
                if (az != az_mater && az != az_unknown && az != az_half)
                    throw hiba(fvnev, "material or half is missing in line %u in %s", i, fileName.c_str());
            } while (az != az_mater && az != az_half); // minimum 1 material szükséges
        }

        if (az == az_half) {
            is_half = true;
            do {
                i++;
                az = azonosit(fajl.lines()[i][0].LowCase());
                if (az != az_mater && az != az_unknown)
                    throw hiba(fvnev, "material is missing in line %u in %s", i, fileName.c_str());
            } while (az != az_mater); // minimum 1 material szükséges
        }

        // material és mat

        while(az==az_mater){
            tmat.resize(tmat.size()+1);
            material & aktmat=tmat[tmat.size()-1];
            aktmat.nev=fajl.lines()[i][1].LowCase();
            aktmat.reflectivity = general_reflectivity;
            if(fajl.lines()[i].size() > 2){ // Fázisváltás vagy hiszterézis
                uns melyik=0;
                if(fajl.lines()[i][2].LowCase()=="f")
                    melyik = 1;
                else if(fajl.lines()[i][2].LowCase()=="h")
                    melyik = 2;
                if(melyik==0)
                    throw hiba(fvnev,"unknown identifier (%s found) in line %u (material) in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                aktmat.is_his=true;
                if(fajl.lines()[i].size() < 5)
                    throw hiba(fvnev,"too few parameters in line %u (material) in %s",i,fileName.c_str());
                if(!fajl.lines()[i][3].todouble(aktmat.his_T_min))
                    throw hiba(fvnev,"not a number (%s found) in line %u (material) in %s",fajl.lines()[i][3].c_str(),i,fileName.c_str());
                if(!fajl.lines()[i][4].todouble(aktmat.his_T_max))
                    throw hiba(fvnev,"not a number (%s found) in line %u (material) in %s",fajl.lines()[i][4].c_str(),i,fileName.c_str());
                if(melyik==2){
                    if(fajl.lines()[i].size() < 6)
                        throw hiba(fvnev,"too few parameters in line %u (material) in %s",i,fileName.c_str());
                    if(!fajl.lines()[i][5].todouble(aktmat.his_T_width_fele))
                        throw hiba(fvnev,"not a number (%s found) in line %u (material) in %s",fajl.lines()[i][5].c_str(),i,fileName.c_str());
                    aktmat.his_T_width_fele*=0.5;
                }
                // aktmat.set_his_derivaltak();
            }
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_mat && az!=az_mater && az!=az_color && az!=az_unknown)
                    throw hiba(fvnev,"material or mat or color is missing (%s found)in line %u in %s",fajl.lines()[i][0].c_str(),i,fileName.c_str());
                if(az==az_mat){
                    matazon maz=azon_mat(fajl.lines()[i][1].LowCase());
                    //uns kelldb=3;
                    //if(fajl.lines()[i].size()!=kelldb)logprint("warning: model::read() => mat arguments may be in wrong format in line %u in %s",i,fileName.c_str());
                    switch(maz){
                        case maz_phase_change_energy:
                            aktmat.is_phase_change_energy = true;
                            if (!fajl.lines()[i][2].todouble(aktmat.phase_change_energy))
                                throw hiba(fvnev, "phase_change_energy is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
                            break;
                        case maz_heatcond:
                        case maz_heatres:
                            aktmat.is_th=true;
                            if(!R_converter(fajl.lines()[i],2,aktmat.thvez,aktmat))
                                throw hiba(fvnev,"heatcond value wrong in line %u in %s",i,fileName.c_str());
                            aktmat.thvez.is_resistivity = ( maz == maz_heatres );
                            break;
                        case maz_heatcap:
                            if(!R_converter(fajl.lines()[i],2,aktmat.Cth,aktmat))
                                throw hiba(fvnev,"heatcap value wrong in line %u in %s",i,fileName.c_str());
                            break;
                        case maz_conduct:
                        case maz_resistivity:
                            aktmat.is_el=true;
                            if(!R_converter(fajl.lines()[i],2,aktmat.elvez,aktmat))
                                throw hiba(fvnev,"conduct value wrong in line %u in %s",i,fileName.c_str());
                            aktmat.elvez.is_resistivity = ( maz == maz_resistivity );
                            break;
                        case maz_capac:
                            if(!R_converter(fajl.lines()[i],2,aktmat.Ce,aktmat))
                                throw hiba(fvnev,"capacity value wrong in line %u in %s",i,fileName.c_str());
                            break;
                        case maz_seebeck:
                            if(!R_converter(fajl.lines()[i],2,aktmat.S,aktmat))
                                throw hiba(fvnev,"seebeck value wrong in line %u in %s",i,fileName.c_str());
                            break;
                        case maz_disscoeff:
                            if(!R_converter(fajl.lines()[i],2,aktmat.D,aktmat))
                                throw hiba(fvnev,"dissip_coeff value wrong in line %u in %s",i,fileName.c_str());
                            break;
                        case maz_emissivity:
                            if(!R_converter(fajl.lines()[i],2,aktmat.emissivity,aktmat))
                                throw hiba(fvnev,"emissivity value wrong in line %u in %s",i,fileName.c_str());
                            break;
                        case maz_reflectivity:
                            if (!fajl.lines()[i][2].todouble(aktmat.reflectivity))
                                throw hiba(fvnev, "reflectivity is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
                            break;
                        case maz_yellow_correction:
                            if (!fajl.lines()[i][2].todouble(aktmat.yellow_correction))
                                throw hiba(fvnev, "yellow_correction is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
                            break;
                        case maz_light_blue_absorption_coeff:
                            aktmat.is_fenypor = true;
                            if (!R_converter(fajl.lines()[i], 2, aktmat.light_blue_absorption_coeff, aktmat))
                                throw hiba(fvnev, "absorption value wrong in line %u in %s", i, fileName.c_str());
                            break;
                        case maz_light_conversion_efficiency:
                            aktmat.is_fenypor = true;
                            if (!R_converter(fajl.lines()[i], 2, aktmat.light_conversion_efficiency, aktmat))
                                throw hiba(fvnev, "conversion_efficiency value wrong in line %u in %s", i, fileName.c_str());
                            break;
                        case maz_light_yellow_absorption_coeff:
                            aktmat.is_fenypor = true;
                            if (!R_converter(fajl.lines()[i], 2, aktmat.light_yellow_absorption_coeff, aktmat))
                                throw hiba(fvnev, "yellow_absorption_coeff value wrong in line %u in %s", i, fileName.c_str());
                            break;
                        case maz_re_conversion_efficiency:
                            aktmat.is_fenypor = true;
                            if (!R_converter(fajl.lines()[i], 2, aktmat.light_re_conversion_efficiency, aktmat))
                                throw hiba(fvnev, "re_conversion_efficiency value wrong in line %u in %s", i, fileName.c_str());
                            break;
                        case maz_vertical_light_conversion:
                            aktmat.is_fenypor = true;
                            if (!R_converter(fajl.lines()[i], 2, aktmat.light_vertical_light_conversion, aktmat))
                                throw hiba(fvnev, "light_vertical_light_conversion value wrong in line %u in %s", i, fileName.c_str());
                            break;
                        case maz_direction_mode:
                            aktmat.is_fenypor = true;
                            if (!fajl.lines()[i][2].tounsigned(aktmat.direction_mode))
                                throw hiba(fvnev, "direction_mode is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
                            if ((aktmat.direction_mode & 8) != 0) {
                                aktmat.direction_mode -= 8;
                                aktmat.is_sarga_szetmegy = true;
                            }
                            else aktmat.is_sarga_szetmegy = false;
                            break;
                        case maz_ray_per_cell_dir:
                            aktmat.is_fenypor = true;
                            if (!fajl.lines()[i][2].tounsigned(aktmat.ray_per_cell_dir))
                                throw hiba(fvnev, "ray_per_cell_dir is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
                            break;
                        case maz_output_side:
                            aktmat.is_fenypor = true;
                            aktmat.output_side = fajl.lines()[i][2][0];
                            if(aktmat.output_side != 'E' && aktmat.output_side != 'W' && aktmat.output_side != 'S' &&
                                aktmat.output_side != 'N' && aktmat.output_side != 'B' && aktmat.output_side != 'T' &&
                                aktmat.output_side != 'A' && aktmat.output_side != 'V' && aktmat.output_side != 'H' &&
                                aktmat.output_side != 'O')
                                throw hiba(fvnev, "output_side is unknown (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
                            break;
                        case maz_cut_level:
                            aktmat.is_fenypor = true;
                            if (!fajl.lines()[i][2].todouble(aktmat.cut_level))
                                throw hiba(fvnev, "dut_level is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
                            break;
                        case maz_d_light_powder:
                            aktmat.is_fenypor = true;
                            if (!fajl.lines()[i][2].todouble(aktmat.d_light_powder))
                                throw hiba(fvnev, "d_light_powder is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
                            break;
                    }
                }
            }while(az==az_unknown || az==az_mat);
        } // legalább 1 color kell

        // color

        while(az==az_color){
            uns szin;
            if(!fajl.lines()[i][1].tounsigned(szin))
                throw hiba(fvnev,"color value wrong in line %u in %s",i,fileName.c_str());
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_col && az!=az_color && az!=az_semicon && az!=az_end && az!=az_coupled && az!=az_unknown)
                    throw hiba(fvnev,"color or col or semiconductor or coupled_model or end is missing (%s found)in line %u in %s",fajl.lines()[i][0].c_str(),i,fileName.c_str());
                if(az==az_col){
                    PLString aktpar=fajl.lines()[i][2].LowCase();
                    PLString kisb=fajl.lines()[i][1].LowCase();
                    if(kisb=="type"){
                        tcolor[szin].is=true;
                        if(aktpar=="normal")tcolor[szin].tipus=SzinNormal;
                        else if (aktpar == "boundary")tcolor[szin].tipus = SzinBoundary;
                        else if (aktpar == "uchannel")tcolor[szin].tipus = SzinUchannel;
                        else logprint("warning: model::read() => unrecognized color type (%s) in line %u in %s", aktpar.c_str(), i, fileName.c_str());
                    }
                    else if(kisb=="material"){
                        uns j;
                        for(j=0;j<tmat.size() && tmat[j].nev!=aktpar;j++);
                        if(j<tmat.size())tcolor[szin].pmat=&tmat[j];
                        else throw hiba(fvnev,"missing material (%s) in color definition in line %u in %s",aktpar.c_str(),i,fileName.c_str());
                    }
                    else if (kisb == "label") {
                        tcolor[szin].nev = aktpar;
                    }
                    else if (kisb == "field") {
                        if(aktpar.LowCase()=="el")
                            tcolor[szin].field=FieldEl;
                        else if (aktpar.LowCase() == "th")
                            tcolor[szin].field = FieldTherm;
                        else tcolor[szin].field = FieldElTherm;
                    }
                }
            }while(az==az_unknown || az==az_col);
        }// semiconductor vagy coupled_model vagy end a vége

        // a bmp színeinek ellenõrzése, térfogatszámítás, tseged feltöltése

        // tseged.resize( x_res * y_res * z_res );
        uns j4 = 0;
        for(uns j1=0;j1<z_res;j1++)
            for(uns j2=0;j2<y_res;j2++)
                for(uns j3=0;j3<x_res;j3++,j4++){
                    color & cakt=tcolor[tbmp[j1].getpixel_also(j3,j2)];
                    if(!cakt.is)
                        throw hiba(fvnev,"undefined color (%u) in bitmap in %s",tbmp[j1].getpixel_also(j3,j2),fileName.c_str());
                    else{
                        cakt.terfogat+=x_pit[j3].get(nulla)*y_pit[j2].get(nulla)*z_pit[j1].get(x_hossz[2*j3]);
                        // tseged[j4].mater = cakt.pmat;
                    }
                }

        // semiconductor

        while(az==az_semicon){
           uns col1,col2;
            if(!fajl.lines()[i][1].tounsigned(col1))
                throw hiba(fvnev,"semiconductor color1 value wrong in line %u in %s",i,fileName.c_str());
            if(!fajl.lines()[i][2].tounsigned(col2))
                throw hiba(fvnev,"semiconductor color2 value wrong in line %u in %s",i,fileName.c_str());
            if(col1==col2)
                throw hiba(fvnev," semiconductor colors must be different (%u) in line %u in %s",col1,i,fileName.c_str());
            if(!tcolor[col1].is||!tcolor[col2].is)
                throw hiba(fvnev," undefined semiconductor color (%u or %u) in line %u in %s",col1,col2,i,fileName.c_str());
            tcolor[col1].tsemi.resize(tcolor[col1].tsemi.size()+1);
            semiconductor & aktsemi=tcolor[col1].tsemi[tcolor[col1].tsemi.size()-1];
            aktsemi.col2=col2;
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_semi && az!=az_semicon && az!=az_coupled && az!=az_end && az!=az_unknown)
                    throw hiba(fvnev,"semiconductor or semi or coupled_model or end is missing (%s found)in line %u in %s",fajl.lines()[i][0].c_str(),i,fileName.c_str());
                if(az==az_semi){
                    PLString aktpar=fajl.lines()[i][2].LowCase();
                    PLString kisb=fajl.lines()[i][1].LowCase();
                    if(kisb=="model"){
                        if(aktpar=="equation"){
                            if(!R_converter(fajl.lines()[i],3,aktsemi.par,m_dummy))
                                throw hiba(fvnev,"equation value wrong in line %u in %s",i,fileName.c_str());
                        }
                        else if(aktpar=="diode"){
                            if(!R_converter(fajl.lines()[i],3,aktsemi.par,m_dummy,true))
                                throw hiba(fvnev,"equation value wrong in line %u in %s",i,fileName.c_str());
                        }
                        else if(aktpar=="erno"){
                            if(!Erno_converter(path+fajl.lines()[i][3],aktsemi.par))
                                throw hiba(fvnev,"erno model wrong in line %u in %s",i,fileName.c_str());
                        }
                        else if (aktpar == "lsfit") {
                            aktsemi.par.semitip = nlt_lsfit;
                            if (!aktsemi.par.ls_fit_adatok.read(path + fajl.lines()[i][3]))
                                throw hiba(fvnev, "LSfit model wrong in line %u in %s", i, fileName.c_str());
                        }
                        else throw hiba(fvnev,"unknown semiconductor model type (%s) in line %u in %s",aktpar.c_str(),i,fileName.c_str());
                    }
                    else if(kisb=="dissip_coeff"){
                        if(!R_converter(fajl.lines()[i],2,aktsemi.D,m_dummy))
                            throw hiba(fvnev,"dissip_coeff value wrong in line %u in %s",i,fileName.c_str());
                    }
                    else if(kisb=="radiation_coeff"){
                        if(!R_converter(fajl.lines()[i],2,aktsemi.R,m_dummy))
                            throw hiba(fvnev,"radiation_coeff value wrong in line %u in %s",i,fileName.c_str());
                    }
                    else if(kisb=="radiance"){
                        if(aktpar=="equation"){
                            if(!R_converter(fajl.lines()[i],3,aktsemi.rad,m_dummy))
                                throw hiba(fvnev,"radiance equation value wrong in line %u in %s",i,fileName.c_str());
                        }
                        else if(aktpar=="erno"){
                            if(!Erno_converter(path+fajl.lines()[i][3],aktsemi.rad))
                                throw hiba(fvnev,"erno radiance model wrong in line %u in %s",i,fileName.c_str());
                        }
                        else if (aktpar == "lsfit") {
                            aktsemi.rad.semitip = nlt_lsfit;
                            if (!aktsemi.rad.ls_fit_adatok.read(path + fajl.lines()[i][3]))
                                throw hiba(fvnev, "LSfit model wrong in line %u in %s", i, fileName.c_str());
                        }
                        else throw hiba(fvnev,"unknown semiconductor radiance type (%s) in line %u in %s",aktpar.c_str(),i,fileName.c_str());
                    }
                    else if(kisb=="luminance"){
                        if(aktpar=="equation"){
                            if(!R_converter(fajl.lines()[i],3,aktsemi.lum,m_dummy))
                                throw hiba(fvnev,"luminance equation value wrong in line %u in %s",i,fileName.c_str());
                        }
                        else if(aktpar=="erno"){
                            if(!Erno_converter(path+fajl.lines()[i][3],aktsemi.lum))
                                throw hiba(fvnev,"erno luminance model wrong in line %u in %s",i,fileName.c_str());
                        }
                        else if (aktpar == "lsfit") {
                            aktsemi.lum.semitip = nlt_lsfit;
                            if (!aktsemi.lum.ls_fit_adatok.read(path + fajl.lines()[i][3]))
                                throw hiba(fvnev, "LSfit model wrong in line %u in %s", i, fileName.c_str());
                        }
                        else throw hiba(fvnev,"unknown semiconductor luminance type (%s) in line %u in %s",aktpar.c_str(),i,fileName.c_str());
                    }
                     else if(kisb=="equation"){
                        throw hiba(fvnev,"semiconductor farmat is deprecated, use model=equation=... instead of equation=... in line %u in %s",i,fileName.c_str());
                    }
                     else if(kisb!="")printf("# unknown semiconductor parameter is ignored: %s\n",kisb.c_str());
               }
            }while(az==az_unknown || az==az_semi);
        }

        // félvezetõ felületek számítása

        A_semi_full=nulla;
        for(uns ix=0;ix<colmax;ix++)
            if(tcolor[ix].is&&tcolor[ix].tipus==SzinNormal)
                for(uns iy=0;iy<tcolor[ix].tsemi.size();iy++){
                    semiconductor & sakt=tcolor[ix].tsemi[iy];
                    sakt.As=nulla;
                    if(sakt.col2==colmax){ // külsõ peremhez félvezet
                        // alja
                        for(uns j2=0;j2<y_res;j2++)
                            for(uns j3=0;j3<x_res;j3++)
                                if(tbmp[0].getpixel_also(j3,j2)==ix)sakt.As+=x_pit[j3].get(nulla)*y_pit[j2].get(nulla);
                        // teteje
                        for(uns j2=0;j2<y_res;j2++)
                            for(uns j3=0;j3<x_res;j3++)
                                if(tbmp[z_res-1].getpixel_also(j3,j2)==ix)sakt.As+=x_pit[j3].get(nulla)*y_pit[j2].get(nulla);
                        // west
                        for(uns j1=0;j1<z_res;j1++)
                            for(uns j2=0;j2<y_res;j2++)
                                if(tbmp[j1].getpixel_also(0,j2)==ix)sakt.As+=z_pit[j1].get(x_hossz[0])*y_pit[j2].get(nulla);
                        // east
                        for(uns j1=0;j1<z_res;j1++)
                            for(uns j2=0;j2<y_res;j2++)
                                if(tbmp[j1].getpixel_also(x_res-1,j2)==ix)sakt.As+=z_pit[j1].get(x_hossz[2*x_res-2])*y_pit[j2].get(nulla);
                        // south
                        for(uns j1=0;j1<z_res;j1++)
                            for(uns j3=0;j3<x_res;j3++)
                                if(tbmp[j1].getpixel_also(j3,0)==ix)sakt.As+=z_pit[j1].get(x_hossz[j3*2])*x_pit[j3].get(nulla);
                        // north
                        for(uns j1=0;j1<z_res;j1++)
                            for(uns j3=0;j3<x_res;j3++)
                                if(tbmp[j1].getpixel_also(j3,y_res-1)==ix)sakt.As+=z_pit[j1].get(x_hossz[j3*2])*x_pit[j3].get(nulla);
                    }
                    else{ // belül félvezet
                        cuns szin1=ix,szin2=sakt.col2;
                        // z irányban
                        for(uns j1=0;j1<z_res-1;j1++)
                            for(uns j2=0;j2<y_res;j2++)
                                for(uns j3=0;j3<x_res;j3++){
                                    cuns s1=tbmp[j1].getpixel_also(j3,j2);
                                    cuns s2=tbmp[j1+1].getpixel_also(j3,j2);
                                    if( (s1==szin1&&s2==szin2) || (s1==szin2&&s2==szin1) )sakt.As+=y_pit[j2].get(nulla)*x_pit[j3].get(nulla);
                                }
                        // y irányban
                        for(uns j1=0;j1<z_res;j1++)
                            for(uns j2=0;j2<y_res-1;j2++)
                                for(uns j3=0;j3<x_res;j3++){
                                    cuns s1=tbmp[j1].getpixel_also(j3,j2);
                                    cuns s2=tbmp[j1].getpixel_also(j3,j2+1);
                                    if( (s1==szin1&&s2==szin2) || (s1==szin2&&s2==szin1) )sakt.As+=z_pit[j1].get(x_hossz[j3*2])*x_pit[j3].get(nulla);
                                }
                        // x irányban
                        for(uns j1=0;j1<z_res;j1++)
                            for(uns j2=0;j2<y_res;j2++)
                                for(uns j3=0;j3<x_res-1;j3++){
                                    cuns s1=tbmp[j1].getpixel_also(j3,j2);
                                    cuns s2=tbmp[j1].getpixel_also(j3+1,j2);
                                    if( (s1==szin1&&s2==szin2) || (s1==szin2&&s2==szin1) )sakt.As+=y_pit[j2].get(nulla)*z_pit[j1].get(x_hossz[j3*2]);
                                }
                    }
                    A_semi_full += sakt.As;
        }
        
        // coupled_model

        uns azon_start = 0;
        while(az==az_coupled){
            hiba("model::read", "kampakt modell csatolas nem tamogatott");
        }
        
    }
    catch(const hiba & h){
        PLString s=h.what();
        if(s.find("tomb::operator[]")==npos)throw;
        throw hiba(fvnev,"missing parameters or unexpected end of file in line %u in %s",i,fileName.c_str());
    }
}


//***********************************************************************
bool bouvalue(const tomb<PLString> & be,boundary & cel,const tomb <convection> & tconv
    ,const PLString & path, uns x_res, uns y_res, tomb<boundary> & peremlista_menteshez, bool isel){
//***********************************************************************
    const PLString azon=be[2].LowCase();
    cel.is_special = false;
    if(azon=="open"||azon=="float"||azon=="adiabatic"){
        cel.tipus=PeremOpen;
    }
    else if(azon=="isothermal"||azon=="equipot"){
        cel.tipus=isel ? PeremV : PeremT;
        if(!be[3].todouble(cel.value))return false;
    }
    else if(azon=="convection"){
        cel.tipus=PeremR;
        if(!be[3].todouble(cel.value)){
            PLString nev=be[3].LowCase();
            for(u32 i=0;i<tconv.size();i++)
                if (tconv[i].nev == nev){
                    cel.conv = tconv[i];
                    cel.is_special = true;
                    goto visszater;
                }
            return false;
        }
        else cel.conv.is_defined = false;
    }
    else if (azon == "conv_temp"){
        cel.tipus = PeremRU;
        if (!be[3].todouble(cel.value)){
            PLString nev = be[3].LowCase();
            for (u32 i = 0; i<tconv.size(); i++)
                if (tconv[i].nev == nev){
                    cel.is_special = true;
                    cel.conv = tconv[i];
                    if (!be[4].todouble(cel.value2))
                        return false;
                    goto visszater;
                }
            return false;
        }
        else{
            cel.conv.is_defined = false;
            if (!be[4].todouble(cel.value2))
                return false;
        }
    }
    else if (azon == "convection-map"){
        cel.tipus=PeremR;
        cel.is_special = true;
        if(!konv_map_konverter(path+be[3],cel.conv_map,x_res,y_res))return false;
    }
    else return false;
visszater:
    if (cel.tipus == PeremOpen) {
        cel.value = nulla;
        cel.v6_index = 0;
    }
    else {
        cel.v6_index = peremlista_menteshez.size() + 1;
        peremlista_menteshez.add(cel);
    }
    return true;
}


//***********************************************************************
bool bouolvas(const tomb<PLString> & be,csomag & cs,bool & isel,const tomb <convection> & tconv
    ,const PLString & path, uns x_res, uns y_res, tomb<boundary> & peremlista_menteshez){
//***********************************************************************
    const PLString azon=be[1].LowCase();
    if(azon=="electrical"){
        isel=true;
        if(!bouvalue(be,cs.el[0],tconv,path,x_res,y_res, peremlista_menteshez,isel))return false;
        for(uns j=1;j<BASIC_SIDES;j++){
            cs.el[j]=cs.el[0];
        }
    }
    else if(azon=="thermal"){
        isel=false;
        if(!bouvalue(be,cs.th[0],tconv,path,x_res,y_res, peremlista_menteshez,isel))return false;
        for(uns j=1;j<BASIC_SIDES;j++){
            cs.th[j]=cs.th[0];
        }
    }
    else if(azon=="internal"){
        throw hiba("bouolvas()","azon==\"internal\"");
    }
    else{
        Oldal oldal=oldalazonosíto(azon);
        if(oldal==EXTERNAL)return false;
        if(isel){
            if(!bouvalue(be,cs.el[oldal],tconv,path,x_res,y_res, peremlista_menteshez,isel))return false;
        }
        else{
            if(!bouvalue(be,cs.th[oldal],tconv,path,x_res,y_res, peremlista_menteshez,isel))return false;
        }
    }

    return true;
}


//***********************************************************************
bool excitvalue(const tomb<PLString> & be,excitation & ertek,tomb<excitation> & mul){
// Ha több gerjesztés érték szerepel, akkor az elsõ kerül az ertek-be, a többi a mul-ba
//***********************************************************************
    const PLString azon=be[2].LowCase();
    if(azon=="voltage"||azon=="temp"||azon=="isothermal"){
        ertek.tipus=GerjU;
        if(!be[3].todouble(ertek.ertek))return false;
    }
    else if(azon=="current"||azon=="power"){
        ertek.tipus=GerjI;
        if(!be[3].todouble(ertek.ertek))return false;
    }
    else throw hiba("excitvalue","unknown excitation type (%s)",azon.c_str());
    if(be.size()>4){
        mul.resize(be.size()-3);
        for(uns i=0;i<mul.size();i++){
            mul[i]=ertek;
            if(!be[i+3].todouble(mul[i].ertek))return false;
        }
    }
    return true;
}


//***********************************************************************
bool excitvalue_controlled(const tomb<PLString> & be,excitation_2 & ertek){
// Controlled analízis fájlhoz
//***********************************************************************
    const PLString azon=be[2].LowCase();
    if(azon=="voltage"||azon=="temp"||azon=="isothermal"){
        ertek.tipus=GerjU;
        if(!be[3].todouble(ertek.ertek))return false;
    }
    else if(azon=="current"||azon=="power"){
        ertek.tipus=GerjI;
        if(!be[3].todouble(ertek.ertek))return false;
    }
    else if(azon=="none")
        ertek.tipus = GerjSemmi;
    else throw hiba("excitvalue_controlled","unknown excitation type (%s)",azon.c_str());
    return true;
}


//***********************************************************************
void read_ctrl(const PLString & fileName,analysis & aktanal, model * pmodel, simulation & sim){
//***********************************************************************
    const char * fvnev="simulation::read_ctrl()";
    logprint("Read %s",fileName.c_str());
    srfajl fajl;
    fajl.open(fileName);
    if(fajl.lines().size()<1||fajl.lines()[0][0].LowCase()!="vsun3-controlled")
        throw hiba(fvnev,"first line is not vsun3-controlled in %s",fileName.c_str());
    uns i=0;
    sorazon az=az_unknown;

    // timestep

    do{
        i++;
        az=azonosit(fajl.lines()[i][0].LowCase());
        if(az!=az_timestep && az!=az_unknown)
            throw hiba(fvnev,"timestep is missing in line %u in %s",i,fileName.c_str());
    }while(az!=az_timestep);
    if(!fajl.lines()[i][1].todouble(aktanal.from)) // az anal.step uns, ezért ebbe tesszük
        throw hiba(fvnev,"timestep is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());

    // timelimit

    do{
        i++;
        az=azonosit(fajl.lines()[i][0].LowCase());
        if(az!=az_timelimit && az!=az_unknown)
            throw hiba(fvnev,"timelimit is missing in line %u in %s",i,fileName.c_str());
    }while(az!=az_timelimit);
    if(!fajl.lines()[i][1].todouble(aktanal.to))
        throw hiba(fvnev,"timelimit is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());

    // change_time

    do{
        i++;
        az=azonosit(fajl.lines()[i][0].LowCase());
        if(az!=az_change_time && az!=az_end && az!=az_unknown)
            throw hiba(fvnev,"change_time is missing in line %u in %s",i,fileName.c_str());
    }while(az!=az_change_time && az!=az_end);

    aktanal.ctrl.clear();
    uns darab=0;
    while(az!=az_end){

        // change_time

        if(az==az_change_time){
            darab++;
            aktanal.ctrl.resize(darab);
            if(!fajl.lines()[i][1].todouble(aktanal.ctrl[darab-1].time))
                throw hiba(fvnev,"change_time is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());
            if(darab>1 && aktanal.ctrl[darab-1].time <= aktanal.ctrl[darab-2].time)
                throw hiba(fvnev,"change_time have to be increasing (%s<=%g) in line %u in %s",
                                  fajl.lines()[i][1].c_str(),aktanal.ctrl[darab-2].time,i,fileName.c_str());
        }

        // electrical vagy thermal

        else if(az==az_electrical || az==az_thermal){
            cuns index=aktanal.ctrl[darab-1].excit.size();
            aktanal.ctrl[darab-1].excit.resize(index+1);

            excitation_2 & aktexcit = aktanal.ctrl[darab - 1].excit[index];
            
            aktexcit.is_el = az==az_electrical;
            
            PLString kisb=fajl.lines()[i][1].LowCase();
            uns szin=colmax;
            for(uns j=0;j<colmax && szin==colmax;j++){
                const color & cc=pmodel->tcolor[j];
                if(cc.is && cc.tipus==SzinNormal && cc.nev==kisb)szin=j;
            }
            if(szin==colmax)
                throw hiba(fvnev,"missing color label (%s) in excitation definition in line %u in %s",fajl.lines()[i][1].c_str(),i,fileName.c_str());
            
            aktexcit.color_index=szin;
            if(!excitvalue_controlled(fajl.lines()[i], aktexcit))
                throw hiba(fvnev,"missing or bad value in excitation definition (%s) in line %u in %s",fajl.lines()[i][1].c_str(),i,fileName.c_str());
        }
        else if(az==az_timestep){
            if(!fajl.lines()[i][1].todouble(aktanal.ctrl[darab-1].timestep))
                throw hiba(fvnev,"timestep is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());
        }
        i++;
        az=azonosit(fajl.lines()[i][0].LowCase());
        if(az!=az_change_time && az!=az_end && az!=az_electrical && az!=az_thermal && az!=az_timestep && az!=az_unknown)
            throw hiba(fvnev,"%s found in line %u in %s",fajl.lines()[i][0].c_str(),i,fileName.c_str());
    }
}


//***********************************************************************
void uchannel::read(PLString path, PLString fajlnev){
//***********************************************************************
    const char * fvnev = "uchannel::init2()";
    logprint("Read %s", fajlnev.c_str());
    srfajl fajl;
    fajl.open(path + fajlnev);

    for (uns i = 0; i < fajl.lines().size(); i++){
        PLString azon = fajl.lines()[i][0].LowCase();
        if (azon == "reverse"){
            is_reverse = true;
            continue;
        }
        if (azon == "segment_number"){
            if (!fajl.lines()[i][1].tounsigned(n))
                throw hiba(fvnev, "segment_number is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "flow_rate"){
            if (!fajl.lines()[i][1].todouble(flow_rate))
                throw hiba(fvnev, "flow_rate is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "fixed_wall_temp"){
            if (!fajl.lines()[i][1].todouble(fixed_wall_temp))
                throw hiba(fvnev, "fixed_wall_temp is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            is_auto_wall_temp = false;
            continue;
        }
        if (azon == "auto_wall_temp"){
            is_auto_wall_temp = true;
            continue;
        }
        if (azon == "fluid_temp"){
            if (!fajl.lines()[i][1].todouble(fluid_temp))
                throw hiba(fvnev, "fluid_temp is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "type"){
            if (fajl.lines()[i][1].LowCase() == "rect")
                tipus = CsatRect;
            else if (fajl.lines()[i][1].LowCase() == "circ")
                tipus = CsatCirc;
            else if (fajl.lines()[i][1].LowCase() == "trap")
                tipus = CsatTrap;
            else
                throw hiba(fvnev, "invalid uchannel type (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "width_bottom"){
            if (!fajl.lines()[i][1].todouble(width_bottom))
                throw hiba(fvnev, "width_bottom is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "width_top"){
            if (!fajl.lines()[i][1].todouble(width_top))
                throw hiba(fvnev, "width_top is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "roughness") {
            if (!fajl.lines()[i][1].todouble(roughness))
                throw hiba(fvnev, "roughness is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "height"){
            if (!fajl.lines()[i][1].todouble(height))
                throw hiba(fvnev, "height is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "density"){
            if (!fajl.lines()[i][1].todouble(density))
                throw hiba(fvnev, "density is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "dynamic_viscosity"){
            if (!fajl.lines()[i][1].todouble(dynamic_visc))
                throw hiba(fvnev, "dynamic_viscosity is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "specific_heat"){
            if (!fajl.lines()[i][1].todouble(spec_heat))
                throw hiba(fvnev, "specific_heat is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "heat_conductivity"){
            if (!fajl.lines()[i][1].todouble(heat_cond))
                throw hiba(fvnev, "heat_conductivity is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        throw hiba(fvnev, "unknown identifier (%s) in %s", fajl.lines()[i][0].c_str(), fajlnev.c_str());
    }
}


//***********************************************************************
void uchannel::set_start_stop_length(){
//***********************************************************************
    if (psim->pmodel->kerekcoord)
        throw hiba("uchannel::set_start_stop_length", "cylindrical coordinata system is mot supported");
    if (szinindex == colmax)
        throw hiba("uchannel::set_start_stop_length", "program error: szinindex == colmax");

    struct pont{ uns x, y, z; pont() :x(0), y(0), z(0){} };
    pont i, min, max;
    pont xmin_ymin, xmin_ymax, xmin_zmin, xmin_zmax;
    pont xmax_ymin, xmax_ymax, xmax_zmin, xmax_zmax;
    pont ymin_xmin, ymin_xmax, ymin_zmin, ymin_zmax;
    pont ymax_xmin, ymax_xmax, ymax_zmin, ymax_zmax;
    pont zmin_xmin, zmin_xmax, zmin_ymin, zmin_ymax;
    pont zmax_xmin, zmax_xmax, zmax_ymin, zmax_ymax;

    bool is_first = true;
    cuns x_res = psim->pmodel->x_res, y_res = psim->pmodel->y_res, z_res = psim->pmodel->z_res;

    // befoglaló téglatest meghatározása pixelben
    // a nyolc csúcshoz legközelebbi pixel meghatározása

    for (i.z = 0; i.z < z_res; i.z++){
        const bitmap & bm = psim->pmodel->tbmp[i.z];
        for (i.y = 0; i.y < y_res; i.y++)
            for (i.x = 0; i.x < x_res; i.x++)
                if (bm.getpixel_also(i.x, i.y) == szinindex)
                    if (is_first){
                        min = max = i;
                        xmin_ymin = xmin_ymax = xmin_zmin = xmin_zmax = xmax_ymin = xmax_ymax = xmax_zmin = xmax_zmax = i;
                        ymin_xmin = ymin_xmax = ymin_zmin = ymin_zmax = ymax_xmin = ymax_xmax = ymax_zmin = ymax_zmax = i;
                        zmin_xmin = zmin_xmax = zmin_ymin = zmin_ymax = zmax_xmin = zmax_xmax = zmax_ymin = zmax_ymax = i;
                        is_first = false;
                    }
                    else{
                        if (i.x < min.x){
                            min.x = i.x;
                            xmin_ymin = xmin_ymax = xmin_zmin = xmin_zmax = i;
                        }
                        else if (i.x == min.x){
                            if (i.y < xmin_ymin.y) xmin_ymin = i;
                            if (i.y > xmin_ymax.y) xmin_ymax = i;
                            if (i.z < xmin_zmin.z) xmin_zmin = i;
                            if (i.z > xmin_zmax.z) xmin_zmax = i;
                        }

                        if (i.x > max.x){
                            max.x = i.x;
                            xmax_ymin = xmax_ymax = xmax_zmin = xmax_zmax = i;
                        }
                        else if (i.x == max.x){
                            if (i.y < xmax_ymin.y) xmax_ymin = i;
                            if (i.y > xmax_ymax.y) xmax_ymax = i;
                            if (i.z < xmax_zmin.z) xmax_zmin = i;
                            if (i.z > xmax_zmax.z) xmax_zmax = i;
                        }

                        if (i.y < min.y){
                            min.y = i.y;
                            ymin_xmin = ymin_xmax = ymin_zmin = ymin_zmax = i;
                        }
                        else if (i.y == min.y){
                            if (i.x < ymin_xmin.x) ymin_xmin = i;
                            if (i.x > ymin_xmax.x) ymin_xmax = i;
                            if (i.z < ymin_zmin.z) ymin_zmin = i;
                            if (i.z > ymin_zmax.z) ymin_zmax = i;
                        }

                        if (i.y > max.y){
                            max.y = i.y;
                            ymax_xmin = ymax_xmax = ymax_zmin = ymax_zmax = i;
                        }
                        else if (i.y == max.y){
                            if (i.x < ymax_xmin.x) ymax_xmin = i;
                            if (i.x > ymax_xmax.x) ymax_xmax = i;
                            if (i.z < ymax_zmin.z) ymax_zmin = i;
                            if (i.z > ymax_zmax.z) ymax_zmax = i;
                        }

                        if (i.z < min.z){
                            min.z = i.z;
                            zmin_xmin = zmin_xmax = zmin_ymin = zmin_ymax = i;
                        }
                        else if (i.z == min.z){
                            if (i.x < zmin_xmin.x) zmin_xmin = i;
                            if (i.x > zmin_xmax.x) zmin_xmax = i;
                            if (i.y < zmin_ymin.y) zmin_ymin = i;
                            if (i.y > zmin_ymax.y) zmin_ymax = i;
                        }

                        if (i.z > max.z){
                            max.z = i.z;
                            zmax_xmin = zmax_xmax = zmax_ymin = zmax_ymax = i;
                        }
                        else if (i.z == max.z){
                            if (i.x < zmax_xmin.x) zmax_xmin = i;
                            if (i.x > zmax_xmax.x) zmax_xmax = i;
                            if (i.y < zmax_ymin.y) zmax_ymin = i;
                            if (i.y > zmax_ymax.y) zmax_ymax = i;
                        }
                    }
    }
 
    // a csatorna kezdõ- és végpontjának meghatározása

    const tomb<dbl> & x_hossz = psim->pmodel->x_hossz;
    const tomb<dbl> & y_hossz = psim->pmodel->y_hossz;
    const tomb<dbl> & z_hossz = psim->pmodel->z_hossz;
    dbl W = x_hossz[2 * max.x] - x_hossz[2 * min.x];
    dbl H = y_hossz[2 * max.y] - y_hossz[2 * min.y];
    dbl L = z_hossz[2 * max.z] - z_hossz[2 * min.z];

    if (W >= H && W >= L){ // x irányú az áramlás
        if (H >= L){ // másodlagos irány az y
            
            // meg kell nézni, hogy fent vagy lent indul-e, és az egyenes kezdõpontja beefoglaló kocka megfelelõ pixelének oldalközepe lesz
            
            if (y_hossz[2 * xmin_ymin.y] - y_hossz[2 * min.y] <= y_hossz[2 * max.y] - y_hossz[2 * xmin_ymax.y]){ // az alsó pixelbõl indul az egyenes
                start.x = (min.x == 0) ? nulla : x_hossz[2 * min.x - 1];
                start.y = y_hossz[2 * min.y];
                start.z = z_hossz[2 * xmin_ymin.z];
            }
            else{ // a felsõ pixelbõl indul az egyenes
                start.x = (min.x == 0) ? nulla : x_hossz[2 * min.x - 1];
                start.y = y_hossz[2 * max.y];
                start.z = z_hossz[2 * xmin_ymax.z];
            }

            // az egyenes végpontja

            if (y_hossz[2 * xmax_ymin.y] - y_hossz[2 * min.y] <= y_hossz[2 * max.y] - y_hossz[2 * xmax_ymax.y]){ // az alsó pixelben végzõdik az egyenes
                stop.x = x_hossz[2 * max.x + 1];
                stop.y = y_hossz[2 * min.y];
                stop.z = z_hossz[2 * xmax_ymin.z];
            }
            else{ // a felsõ pixelben végzõdik az egyenes
                stop.x = x_hossz[2 * max.x + 1];
                stop.y = y_hossz[2 * max.y];
                stop.z = z_hossz[2 * xmax_ymax.z];
            }
        }
        else{ // másodlagos irány a z
            if (z_hossz[2 * xmin_zmin.z] - z_hossz[2 * min.z] <= z_hossz[2 * max.z] - z_hossz[2 * xmin_zmax.z]){ // az alsó pixelbõl indul az egyenes
                start.x = (min.x == 0) ? nulla : x_hossz[2 * min.x - 1];
                start.y = y_hossz[2 * xmin_zmin.y];
                start.z = z_hossz[2 * min.z];
            }
            else{ // a felsõ pixelbõl indul az egyenes
                start.x = (min.x == 0) ? nulla : x_hossz[2 * min.x - 1];
                start.y = y_hossz[2 * xmin_zmax.y];
                start.z = z_hossz[2 * max.z];
            }
            if (z_hossz[2 * xmax_zmin.z] - z_hossz[2 * min.z] <= z_hossz[2 * max.z] - z_hossz[2 * xmax_zmax.z]){ // az alsó pixelben végzõdik az egyenes
                stop.x = x_hossz[2 * max.x + 1];
                stop.y = y_hossz[2 * xmax_zmin.y];
                stop.z = z_hossz[2 * min.z];
            }
            else{ // a felsõ pixelben végzõdik az egyenes
                stop.x = x_hossz[2 * max.x + 1];
                stop.y = y_hossz[2 * xmax_zmax.y];
                stop.z = z_hossz[2 * max.z];
            }
        }
    }
    else if (H >= L){ // y irányú az áramlás
        if (W >= L){ // másodlagos irány az x
            if (x_hossz[2 * ymin_xmin.x] - x_hossz[2 * min.x] <= x_hossz[2 * max.x] - x_hossz[2 * ymin_xmax.x]){ // az alsó pixelbõl indul az egyenes
                start.x = x_hossz[2 * min.x];
                start.y = (min.y == 0) ? nulla : y_hossz[2 * min.y - 1];
                start.z = z_hossz[2 * ymin_xmin.z];
            }
            else{ // a felsõ pixelbõl indul az egyenes
                start.x = x_hossz[2 * max.x];
                start.y = (min.y == 0) ? nulla : y_hossz[2 * min.y - 1];
                start.z = z_hossz[2 * ymin_xmax.z];
            }
            if (x_hossz[2 * ymax_xmin.x] - x_hossz[2 * min.x] <= x_hossz[2 * max.x] - x_hossz[2 * ymax_xmax.x]){ // az alsó piyelben végzõdik az egyenes
                stop.x = x_hossz[2 * min.x];
                stop.y = y_hossz[2 * max.y + 1];
                stop.z = z_hossz[2 * ymax_xmin.z];
            }
            else{ // a felsõ pixelben végzõdik az egyenes
                stop.x = x_hossz[2 * max.x];
                stop.y = y_hossz[2 * max.y + 1];
                stop.z = z_hossz[2 * ymax_xmax.z];
            }
        }
        else{ // másodlagos irány a z
            if (z_hossz[2 * ymin_zmin.z] - z_hossz[2 * min.z] <= z_hossz[2 * max.z] - z_hossz[2 * ymin_zmax.z]){ // az alsó pixelbõl indul az egyenes
                start.x = x_hossz[2 * ymin_zmin.x];
                start.y = (min.y == 0) ? nulla : y_hossz[2 * min.y - 1];
                start.z = z_hossz[2 * min.z];
            }
            else{ // a felsõ pixelbõl indul az egyenes
                start.x = x_hossz[2 * ymin_zmax.x];
                start.y = (min.y == 0) ? nulla : y_hossz[2 * min.y - 1];
                start.z = z_hossz[2 * max.z];
            }
            if (z_hossz[2 * ymax_zmin.z] - z_hossz[2 * min.z] <= z_hossz[2 * max.z] - z_hossz[2 * ymax_zmax.z]){ // az alsó pixelben végzõdik az egyenes
                stop.x = x_hossz[2 * ymax_zmin.x];
                stop.y = y_hossz[2 * max.y + 1];
                stop.z = z_hossz[2 * min.z];
            }
            else{ // a felsõ pixelben végzõdik az egyenes
                stop.x = x_hossz[2 * ymax_zmax.x];
                stop.y = y_hossz[2 * max.y + 1];
                stop.z = z_hossz[2 * max.z];
            }
        }
    }
    else{ // z irányú az áramlás
        if (W >= H){ // másodlagos irány az x
            if (x_hossz[2 * zmin_xmin.x] - x_hossz[2 * min.x] <= x_hossz[2 * max.x] - x_hossz[2 * zmin_xmax.x]){ // az alsó pixelbõl indul az egyenes
                start.x = x_hossz[2 * min.x];
                start.y = y_hossz[2 * zmin_xmin.y];
                start.z = (min.z == 0) ? nulla : z_hossz[2 * min.z - 1];
            }
            else{ // a felsõ pixelbõl indul az egyenes
                start.x = x_hossz[2 * max.x];
                start.y = y_hossz[2 * zmin_xmax.y];
                start.z = (min.z == 0) ? nulla : z_hossz[2 * min.z - 1];
            }
            if (x_hossz[2 * zmax_xmin.x] - x_hossz[2 * min.x] <= x_hossz[2 * max.x] - x_hossz[2 * zmax_xmax.x]){ // az alsó pizelben végyõdik az egyenes
                stop.x = x_hossz[2 * min.x];
                stop.y = y_hossz[2 * zmax_xmin.y];
                stop.z = z_hossz[2 * max.z + 1];
            }
            else{ // a felsõ pixelben végyõdik az egyenes
                stop.x = x_hossz[2 * max.x];
                stop.y = y_hossz[2 * zmax_xmax.y];
                stop.z = z_hossz[2 * max.z + 1];
            }
        }
        else{ // másodlagos irány a y
            if (y_hossz[2 * zmin_ymin.y] - y_hossz[2 * min.y] <= y_hossz[2 * max.y] - y_hossz[2 * zmin_ymax.y]){ // az alsó pixelbõl indul az egyenes
                start.x = x_hossz[2 * zmin_ymin.x];
                start.y = y_hossz[2 * min.y];
                start.z = (min.z == 0) ? nulla : z_hossz[2 * min.z - 1];
            }
            else{ // a felsõ pixelbõl indul az egyenes
                start.x = x_hossz[2 * zmin_ymax.x];
                start.y = y_hossz[2 * max.y];
                start.z = (min.z == 0) ? nulla : z_hossz[2 * min.z - 1];
            }
            if (y_hossz[2 * zmax_ymin.y] - y_hossz[2 * min.y] <= y_hossz[2 * max.y] - y_hossz[2 * zmax_ymax.y]){ // az alsó pixelben végyõdik az egyenes
                stop.x = x_hossz[2 * zmax_ymin.x];
                stop.y = y_hossz[2 * min.y];
                stop.z = z_hossz[2 * max.z + 1];
            }
            else{ // a felsõ pixelben végyõdik az egyenes
                stop.x = x_hossz[2 * zmax_ymax.x];
                stop.y = y_hossz[2 * max.y];
                stop.z = z_hossz[2 * max.z + 1];
            }
        }
    }
    // Mindig x növekedése irányába folyik. Ha függõleges, akkor y növekedése irányába.
    if (fabs(start.x - stop.x) < 1e-9){
        if (start.y > stop.y)
            swap(start, stop);
    }
    else if (start.x > stop.x)
        swap(start, stop);

    if (is_reverse)
        swap(start, stop);

    // hossz kiszámítása
    
    length = Abs(stop - start);
}


//***********************************************************************
void uchannel::set_uchannel_in_bmp(){
//***********************************************************************
    // Minden szinindex színû pixel középpontháról eldönti, hogy a 
    // csatorna mely szakaszához tartozik az n-bõl, és beállítja ennek
    // indexét (azaz tucha_perem[szinindex][i]-hez rendeli)

    // ALGORITMUS:
    // 1. Ha van két vektorunk: e és f (közös pontból: az origóból indulnak),
    //    akkor skaláris szorzatuk e*f=|e|*|f|*cos(alfa), alfa az általuk
    //    bezárt szög. Ez az egyik vektor merõleges vetületét jelenti a másik
    //    vektorra. Vagyis |f|*cos(alfa)/|e| megmondja, hogy az f által kije-
    //    lölt ponthoz legközelebbi pont az origó-e szakaszon hol van, milyen
    //    arányban osztja ketté a szakaszt. (Értéke lehet negatív nagy 1-nél
    //    nagyobb is, ha az f pont nem az origo-e szakasz mellett van.
    // 2. e = stop - start, f = cellakozep - start
    // 3. e*f = |e|*|f|*cos(alfa) => |f|*cos(alfa) = e*f/|e|
    //    A vetítési arány tehát: 
    //    va = |f|*cos(alfa)/|e| = e*f/|e|^2 = (e*f)/(e*e)
    //    if(va<0.0) va = 0; if(va>=1.0) va = 1.0 - 1.0e-010; (a -1e-10 azért kell, hogy max. n-1-re kerekítsen.
    //    A keresett index = (int)(n*va);

    tomb<uns> letezo_indexek;

    cuns x_res = psim->pmodel->x_res, y_res = psim->pmodel->y_res, z_res = psim->pmodel->z_res;
    const tomb<dbl> & x_hossz = psim->pmodel->x_hossz;
    const tomb<dbl> & y_hossz = psim->pmodel->y_hossz;
    const tomb<dbl> & z_hossz = psim->pmodel->z_hossz;
    const Vec3d e = stop - start;
    cd egy_per_hossz_negyzet = egy / (e*e);
    cd majdnemegy = egy - 1.0e-010;

    for (uns z = 0; z < z_res; z++){
        bitmap & bm = psim->pmodel->tbmp[z];
        for (uns y = 0; y < y_res; y++)
            for (uns x = 0; x < x_res; x++)
                if (bm.getpixel_also(x, y) == szinindex){
                    const Vec3d cellakozep(x_hossz[2 * x], y_hossz[2 * y], z_hossz[2 * z]);
                    const Vec3d f = cellakozep - start;
                    cd va = e*f*egy_per_hossz_negyzet;
                    cd vacut = va < nulla ? nulla : va >= egy ? majdnemegy : va;
                    cuns index = (uns)(n*vacut);
                    bm.setpixel_felso(x, y, index);
                    letezo_indexek.insertIfNotExists(index);
                }
    }

    // Ha a felhasználó túl nagy n-t adott meg, töröljük a felesleges szegmenseket, és n értékét a megfelelõ értékre csökkentjük

    for (uns z = 0; z < z_res; z++){
        bitmap & bm = psim->pmodel->tbmp[z];
        for (uns y = 0; y < y_res; y++)
            for (uns x = 0; x < x_res; x++)
                if (bm.getpixel_also(x, y) == szinindex)
                    if (letezo_indexek.findSorted(bm.getpixel_felso(x, y)))
                        bm.setpixel_felso(x, y, letezo_indexek.getfindindex());
                    else
                        throw hiba("uchannel::set_uchannel_in_bmp", "Impossible index");
    }
    n = letezo_indexek.size();
}

//***********************************************************************
void uchannel::init2(PLString path, PLString fajlnev, const PLString &cimke){
//***********************************************************************
// beolvassa a fájlt, az alapján beállítja a paramétereket
// meghatározza a csatorna kezdõ és végpontját, kiszámolja a csatornahosszt, felosztja n részre
// beállítja a modell tbmp tömbjében a csatornához tartozó színû pixelek felsõ értékének a csatornarész sorszámát
// átméretezi tucha_perem tömb szinindexedik elemét, a vezetéseknek kezdõértéket ad

    // kezdeti értékek

    nev = cimke;
    n = 8;
    tipus = CsatTrap;
    is_auto_wall_temp = true;
    is_reverse = false;
    flow_rate = 1.0;
    fixed_wall_temp = 60.0;
    fluid_temp = 0.0;
    width_bottom = 250e-6;
    width_top = 350e-6;
    roughness = 1.0e-6;
    height = 67e-6;
    density = 1.1614;
    dynamic_visc = 1.84e-5;
    spec_heat = 1005.0;
    heat_cond = 0.0261;

    read(path, fajlnev);
    set_start_stop_length();
    set_uchannel_in_bmp();

    psim->tucha_perem[szinindex].resize(n);
}

//***********************************************************************
void powermap::map_t::read(srfajl & fajl, uns & sor, u32 x, u32 y, u32 z){
//***********************************************************************
    pmap.resize(x, y, z);
    uns sorindex = 0;
    for (uns k = 0; k < z; k++)
        for (uns j = 0; j < y; j++)
            for (uns i = 0; i < x; i++, sorindex++) {
                if (sorindex >= fajl.lines()[sor].size()) {
                    sorindex = 0;
                    sor++;
                }
                if (!fajl.lines()[sor][sorindex].todouble(pmap.getref(i, j, k)))
                    throw hiba("powermap::map_t::read", "power value is not a number (%s) (t=%g)", fajl.lines()[sor][sorindex].c_str(),t);
            }
    sor++;
}


//***********************************************************************
void powermap::read(PLString path, PLString fajlnev) {
//***********************************************************************
    const char * fvnev = "powermap::read()";
    logprint("Read %s", fajlnev.c_str());
    srfajl fajl;
    fajl.open(path + fajlnev);

    for (uns i = 0; i < fajl.lines().size(); i++) {
        PLString azon = fajl.lines()[i][0].LowCase();
        if (azon == "vsun3-power-map")
            continue;
        if (azon == "top") {
            hol = Hol::top;
            continue;
        }
        if (azon == "volume") {
            hol = Hol::volume;
            continue;
        }
        if (azon == "exact") {
            is_exact = true;
            continue;
        }
        if (azon == "x-resolution") {
            if (!fajl.lines()[i][1].tounsigned(x))
                throw hiba(fvnev, "segment_number is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "y-resolution") {
            if (!fajl.lines()[i][1].tounsigned(y))
                throw hiba(fvnev, "segment_number is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "z-resolution") {
            if (!fajl.lines()[i][1].tounsigned(z))
                throw hiba(fvnev, "segment_number is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
            continue;
        }
        if (azon == "x-pitch") {
            if (!pitchconverter(fajl.lines()[i], x_pitch))
                throw hiba(fvnev, "x-pitch data corrupt in line %u in %s", i, fajlnev.c_str());
            if (x_pitch.size()<x)
                throw hiba(fvnev, "fewer x-pitch data than resolution (%u<%u) in line %u in %s", 
                    x_pitch.size(), x, i, fajlnev.c_str());
            if (x_pitch.size()>x)
                logprint("warning: model::read() => more x-pitch data than resolution (%u>%u) in line %u in %s", 
                    x_pitch.size(), x, i, fajlnev.c_str());
            continue;
        }
        if (azon == "y-pitch") {
            if (!pitchconverter(fajl.lines()[i], y_pitch))
                throw hiba(fvnev, "y-pitch data corrupt in line %u in %s", i, fajlnev.c_str());
            if (y_pitch.size()<y)
                throw hiba(fvnev, "fewer y-pitch data than resolution (%u<%u) in line %u in %s", 
                    y_pitch.size(), y, i, fajlnev.c_str());
            if (y_pitch.size()>y)
                logprint("warning: model::read() => more y-pitch data than resolution (%u>%u) in line %u in %s", 
                    y_pitch.size(), y, i, fajlnev.c_str());
            continue;
        }
        if (azon == "z-pitch") {
            if (!pitchconverter(fajl.lines()[i], z_pitch))
                throw hiba(fvnev, "z-pitch data corrupt in line %u in %s", i, fajlnev.c_str());
            if (z_pitch.size()<z)
                throw hiba(fvnev, "fewer z-pitch data than resolution (%u<%u) in line %u in %s", 
                    z_pitch.size(), z, i, fajlnev.c_str());
            if (z_pitch.size()>z)
                logprint("warning: model::read() => more z-pitch data than resolution (%u>%u) in line %u in %s", 
                    z_pitch.size(), z, i, fajlnev.c_str());
            continue;
        }
        if (azon == "t") {
            uns db = 0;
            for (uns j = i; j < fajl.lines().size(); j++) {
                PLString az = fajl.lines()[j][0].LowCase();
                if (az == "t")
                    db++;
            }
            tombok.resize(db);
            if (hol == Hol::top)
                z = 1;
            for (uns j = 0; j < db; j++) {
                PLString az = fajl.lines()[i][0].LowCase();
                if(az!="t")
                    throw hiba(fvnev, "identifier should be \"t\" instead of %s in %u. iteration of %u in line %u in %s", az.c_str(), j, db, i, fajlnev.c_str());
                if (!fajl.lines()[i][1].todouble(tombok[j].t))
                    throw hiba(fvnev, "t is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fajlnev.c_str());
                i++;
                tombok[j].read(fajl, i, x, y, z);
            }
            continue;
        }
        if (azon == "stop") {
            is_stop = true;
            continue;
        }
        if (azon == "end") {
            if (!is_exact) {
                if (x_pitch.size()<x)throw hiba(fvnev, "x-pitch is missing in %s", fajlnev.c_str());
                if (y_pitch.size()<y)throw hiba(fvnev, "y-pitch is missing in %s", fajlnev.c_str());
                if (hol==Hol::volume && z_pitch.size()<z)throw hiba(fvnev, "z-pitch is missing in %s", fajlnev.c_str());
            }
            break;
        }
        throw hiba(fvnev, "unknown identifier (%s) in %s", fajl.lines()[i][0].c_str(), fajlnev.c_str());
    }
    van_e = true;
}


//***********************************************************************
void powermap::get_map(dbl t, dbl * map){
// A map psim->pmodel->x_res * psim->pmodel->y_res * psim->pmodel->z_res
// méretû kell legyen. Top map esetén az alsóbbakat nullázza.
//***********************************************************************
    if (psim->pmodel->x_res != x)
        throw hiba("powermap::get_map", "model and powermap x-resolution is not equal (%u!=%u)", psim->pmodel->x_res, x);
    if (psim->pmodel->y_res != y)
        throw hiba("powermap::get_map", "model and powermap y-resolution is not equal (%u!=%u)", psim->pmodel->y_res, y);
    if (hol == Hol::volume && psim->pmodel->z_res != z)
        throw hiba("powermap::get_map", "model and powermap z-resolution is not equal (%u!=%u)", psim->pmodel->z_res, z);

    // Az elsõ t-ig nincs disszipáció, vagy ha stop be van állítva, akkor utolsó t után nincs disszipáció

    if (tombok[0].t >= t + 1e-12 || (is_stop && tombok[tombok.size()-1].t < t - 1e-12)) {
        cuns db = x*y*psim->pmodel->z_res;
        for (uns i = 0; i < db; i++)
            map[i] = nulla;
        return;
    }    
    
    // amikortól a disszipáció definiálva van, hányadik mape-et használja

    uns hanyadik = 0;
    for (uns i = 1; i < tombok.size(); i++)
        if (tombok[i].t >= t + 1e-12) {
            hanyadik = i - 1;
            break;
        }
        else if (i == tombok.size() - 1)
            hanyadik = i;

    switch (hol) {
        case Hol::top: {
            int n = 0;
            for (uns k = 0; k < psim->pmodel->z_res - 1; k++)
                for (uns j = 0; j < y; j++)
                    for (uns i = 0; i < x; i++, n++) {
                        map[n] = nulla;
                    }
            cuns eltolas = x * y * (psim->pmodel->z_res - 1);
            for (uns k = psim->pmodel->z_res - 1; k < psim->pmodel->z_res; k++)
                for (uns j = 0; j < y; j++)
                    for (uns i = 0; i < x; i++, n++) {
                        map[n] = tombok[hanyadik].pmap.getconstref(n - eltolas);
                    }
        }
            break;
        case Hol::volume: {
                int n = 0;
                for (uns k = 0; k < z; k++)
                    for (uns j = 0; j < y; j++)
                        for (uns i = 0; i < x; i++, n++) {
                            map[n] = tombok[hanyadik].pmap.getconstref(n);
                        }
            }
            break;
    }
}

#include <thread>

//***********************************************************************
void simulation::read(PLString path){
//***********************************************************************
    const char * fvnev="simulation::read()";
    logprint("Read %s",fileName.c_str());
    srfajl fajl;
    fajl.open(path+fileName);
    if(fajl.lines().size()<1||fajl.lines()[0][0].LowCase()!="vsun3-simulation")
        throw hiba(fvnev,"first line is not vsun3-model in %s",fileName.c_str());
    uns i=0;
    sorazon az=az_unknown;

    // name

    try{
        do{
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_nev && az!=az_unknown)throw hiba(fvnev,"name is missing in line %u in %s",i,fileName.c_str());
        }while(az!=az_nev);
        name=fajl.lines()[i][1];

        // field

        do{
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_field&& az!=az_unknown)throw hiba(fvnev,"field is missing in line %u in %s",i,fileName.c_str());
        }while(az!=az_field);
        PLString kisb=fajl.lines()[i][1].LowCase();
        if(kisb=="electrical")mezo=FieldEl;
        else if(kisb=="thermal")mezo=FieldTherm;
        else if(kisb=="eltherm")mezo=FieldElTherm;
        else throw hiba(fvnev,"unknown field type (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());

        // computation

        bool cango=false;
        do{ // az!=az_comp && az!=az_bound && az!=az_convection && az!=az_cpu && az!=az_hdd && az!=az_optimize
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            switch(az){
                case az_calc_mod:  
                    if (!fajl.lines()[i][1].tounsigned(mezo_szamitasi_mod))
                        throw hiba(fvnev, "calc_mode is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
                    if (mezo_szamitasi_mod < 1 || mezo_szamitasi_mod > 3)
                        throw hiba(fvnev, "calc_mode must be 1..3 but %u was set in %s", mezo_szamitasi_mod, fileName.c_str());
                    break;
                case az_lin:       is_lin=true;         break;
                case az_nosemi:    is_no_semi=true;     break;
                case az_noseebeck: is_no_seebeck=true;  break;
                case az_nopeltier: is_no_peltier=true;  break;
                case az_peltier_center: is_peltier_center=true;  break;
                case az_nothomson: is_no_thomson=true;  break;
                case az_no_joule:  is_no_joule=true;    break;
                // A read_ctrl true-ra állíthatja az is_no_excit_as_bou-t.
                case az_noaccelerator: is_no_accelerator=true;  break;
                case az_cm_temp:   aramot=false;        break;
                case az_el_nonlin_subiter:
                    if(!fajl.lines()[i][1].tounsigned(el_nonlin_subiter))
                        throw hiba(fvnev,"el_nonlin_subiter is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());
                    break;
                case az_ndc_miniter:
                    if (!fajl.lines()[i][1].tounsigned(ndc_miniter))
                        throw hiba(fvnev, "ndc_miniter is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
                    break;
                case az_comp:
                    if(!fajl.lines()[i][1].tounsigned(valostipus))
                        throw hiba(fvnev,"computation is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());
                    break;
                case az_always_strassen_mul: is_always_strassen_mul=true;  break;
                case az_always_quad: is_always_quad=true;  break;
                case az_bound:      cango=true;   break;
                case az_convection: cango=true;   break;
                case az_cpu:
                    if (!fajl.lines()[i][1].tounsigned(cpu_threads))
                        if (fajl.lines()[i][1].LowCase() == "auto")
                            cpu_threads = 32767;
                        else
                            throw hiba(fvnev,"unrecognizable cpu-threads value (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());
                    break;
                case az_hdd:
                    if(!fajl.lines()[i][1].tounsigned(hdd))
                        throw hiba(fvnev,"hdd-cache is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());
                    break;
                case az_optimize:
                    if(!fajl.lines()[i][1].tounsigned(optimize))
                        throw hiba(fvnev,"optimize is not a number (%s) in %s",fajl.lines()[i][1].c_str(),fileName.c_str());
                    break;
                case az_unknown:                break;
                default: throw hiba(fvnev,"linear or no_semicond or no_seebeck or computation or cpu_threads or hdd-cache or optimize or boundary is missing in line %u in %s",i,fileName.c_str());
            }
        }while(!cango);

        uns tenyleges_szalszam = std::thread::hardware_concurrency();
        if (cpu_threads == 32767) {
            cpu_threads = tenyleges_szalszam;
            cpu_threads = cpu_threads == 0 ? 16 : cpu_threads;
        }

        // convection

        material m_dummy;
        bool isel;
        while(az!=az_bound){
            switch(az){
                case az_unknown: break;
                case az_convection:{
                        cu32 n=tconv.incsize();
                        tconv[n-1].nev=fajl.lines()[i][1].LowCase();
                        tconv[n-1].is_defined=true;
                        if(tconv[n-1].nev=="")throw hiba(fvnev,"missing convection name in line %u in %s",i,fileName.c_str());
                    }
                    break;
                case az_conv:{
                        cu32 n=tconv.size();
                        switch(azon_conv(fajl.lines()[i][1].LowCase())){
                            case caz_unknown: break;
                            case caz_rad:
                                if(!R_converter(fajl.lines()[i],2,tconv.getLast().radiation,m_dummy))
                                    throw hiba(fvnev,"conv>radiation missing multiplier value in line %u in %s",i,fileName.c_str());
                                break;
                            case caz_vertical:
                                if(fajl.lines()[i][2].LowCase()=="htc")             tconv[n-1].vertical_tipus=ConvHTC;
                                else if(fajl.lines()[i][2].LowCase()=="vertical-1") tconv[n-1].vertical_tipus=ConvVertical_1; 
                                else if(fajl.lines()[i][2].LowCase()=="churchill-p1") tconv[n-1].vertical_tipus=ConvVerticalChurchill_P_1; 
                                else if(fajl.lines()[i][2].LowCase()=="churchill-p2") tconv[n-1].vertical_tipus=ConvVerticalChurchill_P_2; 
                                else if(fajl.lines()[i][2].LowCase()=="churchill-t") tconv[n-1].vertical_tipus=ConvVerticalChurchill_T; 
                                else if(fajl.lines()[i][2].LowCase()=="churchill-c") tconv[n-1].vertical_tipus=ConvVerticalChurchill_C; 
                                else if(fajl.lines()[i][2].LowCase()=="lee-t") tconv[n-1].vertical_tipus=ConvVerticalLee_T;
                                else if(fajl.lines()[i][2].LowCase()=="lee-p") tconv[n-1].vertical_tipus=ConvVerticalLee_P;
                                else if(fajl.lines()[i][2].LowCase()=="mihajev") tconv[n-1].vertical_tipus=ConvVerticalMihajev;
                                else if(fajl.lines()[i][2].LowCase()=="wei") tconv[n-1].vertical_tipus=ConvWei;
                                else throw hiba(fvnev,"conv>free-vertical unknown type \"%s\" in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                                
                                if(!R_converter(fajl.lines()[i],3,tconv[n-1].vertical_value,m_dummy))
                                    throw hiba(fvnev,"conv>free-vertical multiplier value in line %u in %s",i,fileName.c_str());
                                break;
                            case caz_lower:
                                if(fajl.lines()[i][2].LowCase()=="htc")          tconv[n-1].lower_tipus=ConvHTC;
                                else if(fajl.lines()[i][2].LowCase()=="lower-1") tconv[n-1].lower_tipus=ConvLower_1;
                                else if(fajl.lines()[i][2].LowCase()=="yovanovich-min") tconv[n-1].lower_tipus=ConvYovanovichMin;
                                else if(fajl.lines()[i][2].LowCase()=="yovanovich-max") tconv[n-1].lower_tipus=ConvYovanovichMax;
                                else if(fajl.lines()[i][2].LowCase()=="wei") tconv[n-1].lower_tipus=ConvWei;
                                else throw hiba(fvnev,"conv>free-lower unknown type \"%s\" in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                                
                                if(!R_converter(fajl.lines()[i],3,tconv[n-1].lower_value,m_dummy))
                                    throw hiba(fvnev,"conv>free-lower multiplier value in line %u in %s",i,fileName.c_str());
                                break;
                            case caz_upper:
                                if(fajl.lines()[i][2].LowCase()=="htc")          tconv[n-1].upper_tipus=ConvHTC;
                                else if(fajl.lines()[i][2].LowCase()=="upper-1") tconv[n-1].upper_tipus=ConvUpper_1;
                                else if(fajl.lines()[i][2].LowCase()=="wei") tconv[n-1].upper_tipus=ConvWei;
                                else throw hiba(fvnev,"conv>free-upper unknown type \"%s\" in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                                
                                if(!R_converter(fajl.lines()[i],3,tconv[n-1].upper_value,m_dummy))
                                    throw hiba(fvnev,"conv>rfree-upper multiplier value in line %u in %s",i,fileName.c_str());
                                break;
                            case caz_axis:
                                if(fajl.lines()[i][2].LowCase()=="x")       tconv[n-1].axis=X_IRANY;
                                else if(fajl.lines()[i][2].LowCase()=="y")  tconv[n-1].axis=Y_IRANY;
                                else if(fajl.lines()[i][2].LowCase()=="z")  tconv[n-1].axis=Z_IRANY;
                                else throw hiba(fvnev,"conv>axis unknown axis \"%s\" in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                                break;
                            case caz_angle:
                                if(!fajl.lines()[i][2].todouble(tconv[n-1].angle))
                                    throw hiba(fvnev,"conv>angle not a number \"%s\" in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                                if(tconv[n-1].angle<=-360 || tconv[n-1].angle>=360)
                                    throw hiba(fvnev,"conv>angle is not in -360<\"%s\"<360 range in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                                if(tconv[n-1].angle<0)tconv[n-1].angle+=360;
                                break;
                            case caz_edge:
                                if(fajl.lines()[i][2].LowCase()=="wei-h")        tconv[n-1].edge=ConvEdgeWei_H;
                                else if(fajl.lines()[i][2].LowCase()=="wei-i")   tconv[n-1].edge=ConvEdgeWei_I;
                                else if(fajl.lines()[i][2].LowCase()=="wei-hi")  tconv[n-1].edge=ConvEdgeWei_HI;
                                else if(fajl.lines()[i][2].LowCase()=="wei-hh")  tconv[n-1].edge=ConvEdgeWei_HH;
                                else throw hiba(fvnev,"conv>edge unknown type \"%s\" in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                                break;
                        }
                    }
                    break;
            }
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_bound && az!=az_convection && az!=az_conv && az!=az_unknown)
                throw hiba(fvnev,"convection or conv or boundary is missing in line %u in %s",i,fileName.c_str());
        }

        // boundary, electrical

        if(mezo==FieldEl||mezo==FieldElTherm){
            if(fajl.lines()[i][1].LowCase()!="electrical")
                throw hiba(fvnev,"missing electrical boundary in line %u in %s",i,fileName.c_str());
            if(!bouolvas(fajl.lines()[i],normalperem,isel,tconv,path,pmodel->x_res,pmodel->y_res,peremlista_menteshez))
                throw hiba(fvnev,"electrical boundary type (%s) or value is wrong in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_bound && az!=az_bou && az!=az_excita && az!=az_ambient && az!=az_unknown)
                    throw hiba(fvnev,"boundary or bou or excitation or ambient is missing in line %u in %s",i,fileName.c_str());
            }while(az==az_unknown);
            while(az==az_bou){
                if(!bouolvas(fajl.lines()[i],normalperem,isel,tconv,path,pmodel->x_res,pmodel->y_res, peremlista_menteshez))
                    throw hiba(fvnev,"electrical boundary type (%s) or value is wrong in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                do{
                    i++;
                    az=azonosit(fajl.lines()[i][0].LowCase());
                    if(az!=az_bound && az!=az_bou && az!=az_excita && az!=az_ambient && az!=az_unknown)
                        throw hiba(fvnev,"boundary or bou or excitation or ambient is missing in line %u in %s",i,fileName.c_str());
                }while(az==az_unknown);
            }
            if(pmodel->dim < 3){
                normalperem.el[BOTTOM].tipus = normalperem.el[TOP].tipus = PeremOpen;
            }
        }

        // boundary, thermal

        if(mezo==FieldTherm||mezo==FieldElTherm){
            while(az!=az_bound){
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_bound && az!=az_unknown)
                    throw hiba(fvnev,"thermal boundary is missing in line %u in %s",i,fileName.c_str());
            }
            if(fajl.lines()[i][1].LowCase()!="thermal")
                throw hiba(fvnev,"missing thermal boundary in line %u in %s",i,fileName.c_str());
            if(!bouolvas(fajl.lines()[i],normalperem,isel,tconv,path,pmodel->x_res,pmodel->y_res, peremlista_menteshez))
                throw hiba(fvnev,"thermal boundary type (%s) or value is wrong in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_bound && az!=az_bou && az!=az_excita && az!=az_ambient && az!=az_unknown)
                    throw hiba(fvnev,"boundary or bou or excitation or ambient is missing in line %u in %s",i,fileName.c_str());
            }while(az==az_unknown);
            while(az==az_bou){
                if(!bouolvas(fajl.lines()[i],normalperem,isel,tconv,path,pmodel->x_res,pmodel->y_res, peremlista_menteshez))
                    throw hiba(fvnev,"thermal boundary type (%s) or value is wrong in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                do{
                    i++;
                    az=azonosit(fajl.lines()[i][0].LowCase());
                    if(az!=az_bound && az!=az_bou && az!=az_excita && az!=az_ambient && az!=az_unknown)
                        throw hiba(fvnev,"boundary or bou or excitation or ambient is missing in line %u in %s",i,fileName.c_str());
                }while(az==az_unknown);
            }
            if(pmodel->dim < 3){
                normalperem.th[BOTTOM].tipus = normalperem.th[TOP].tipus = PeremOpen;
            }
        }

        // boundary, internal

        tinner.clear();
        uns darab=0;
        for(uns j=0;j<colmax;j++)if(pmodel->tcolor[j].is && pmodel->tcolor[j].tipus==SzinBoundary)darab++;
        for(uns j=0;j<colmax;j++)innerindex[j]=colmax;
        tinner.resize(darab);
        for(;darab>0;darab--){
            while(az!=az_bound){
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_bound && az!=az_unknown)
                    throw hiba(fvnev,"internal boundary is missing in line %u in %s",i,fileName.c_str());
            }
            if(fajl.lines()[i][1].LowCase()!="internal")
                throw hiba(fvnev,"missing internal boundary in line %u in %s",i,fileName.c_str());
            const PLString cimke=fajl.lines()[i][2].LowCase();
            for(uns j=0;j<colmax && tinner[darab-1].szin==colmax;j++)
                if(pmodel->tcolor[j].is && pmodel->tcolor[j].tipus==SzinBoundary && cimke==pmodel->tcolor[j].nev)tinner[darab-1].szin=j;
            if(tinner[darab-1].szin==colmax)
                throw hiba(fvnev,"undefined internal boundary label (%s) in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
            
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_bound && az!=az_bou && az!=az_excita && az!=az_ambient && az!=az_unknown)
                    throw hiba(fvnev,"boundary or bou or excitation or ambient is missing in line %u in %s",i,fileName.c_str());
            }while(az==az_unknown);
            bool voltel=false,voltth=false;
            while(az==az_bou){
                if(fajl.lines()[i][1].LowCase()=="electrical")voltel=true;
                if(fajl.lines()[i][1].LowCase()=="thermal")voltth=true;
                if(!bouolvas(fajl.lines()[i],tinner[darab-1],isel,tconv,path,pmodel->x_res,pmodel->y_res, peremlista_menteshez))
                    throw hiba(fvnev,"internal boundary type (%s) or value is wrong in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                do{
                    i++;
                    az=azonosit(fajl.lines()[i][0].LowCase());
                    if(az!=az_bound && az!=az_bou && az!=az_excita && az!=az_ambient && az!=az_unknown)
                        throw hiba(fvnev,"boundary or bou or excitation or ambient is missing in line %u in %s",i,fileName.c_str());
                }while(az==az_unknown);
            }
            if((mezo==FieldEl||mezo==FieldElTherm)&&!voltel)
                throw hiba(fvnev,"missing electrical internal boundary definition in line %u in %s",i,fileName.c_str());
            if((mezo==FieldTherm||mezo==FieldElTherm)&&!voltth)
                throw hiba(fvnev,"missing thermal internal boundary definition in line %u in %s",i,fileName.c_str());
        }
        for(uns j=0;j<tinner.size();j++)innerindex[tinner[j].szin]=j;

        // excitation

        while(az!=az_excita && az!=az_ambient){
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if (az==az_bou || az==az_bound)
                throw hiba(fvnev, "unnecessary boundary conditin in line %u in %s (internal boundary is turned off?)", i, fileName.c_str());
            if (az != az_excita && az != az_ambient && az != az_unknown)
                throw hiba(fvnev,"excitation or ambient is missing in line %u in %s",i,fileName.c_str());
        }
        while(az==az_excita){
            bool isel=false;
            const PLString kisb=fajl.lines()[i][1].LowCase();
            if(kisb=="electrical")isel=true;
            else if(kisb=="thermal")isel=false;
            else throw hiba(fvnev,"excitation type is wrong (%s) in line %u in %s",fajl.lines()[i][1].c_str(),i,fileName.c_str());
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_excit && az!=az_excita && az!=az_ambient && az!=az_unknown)
                    throw hiba(fvnev,"excitation or excit or ambient is missing (%s found)in line %u in %s",fajl.lines()[i][0].c_str(),i,fileName.c_str());
                if(az==az_excit){
                    PLString kisb=fajl.lines()[i][1].LowCase();
                    uns szin=colmax;
                    for(uns j=0;j<colmax && szin==colmax;j++){
                        const color & cc=pmodel->tcolor[j];
                        if(cc.is && cc.tipus==SzinNormal && cc.nev==kisb)szin=j;
                    }
                    if(szin==colmax)
                        throw hiba(fvnev,"missing color label (%s) in excit definition in line %u in %s",fajl.lines()[i][1].c_str(),i,fileName.c_str());
                    if(isel){
                        texcitE[szin].is=true;
                        if(!excitvalue(fajl.lines()[i],texcitE[szin],mulE[szin]))                        
                            throw hiba(fvnev,"missing or bad excitation value in excit definition (%s) in line %u in %s",fajl.lines()[i][1].c_str(),i,fileName.c_str());
                    }
                    else{
                        texcitT[szin].is=true;
                        if(!excitvalue(fajl.lines()[i],texcitT[szin],mulT[szin]))                        
                            throw hiba(fvnev,"missing or bad excitation value in excit definition (%s) in line %u in %s",fajl.lines()[i][1].c_str(),i,fileName.c_str());
                    }
                }
            }while(az==az_unknown || az==az_excit);
        }

        // ambient

        while(az!=az_ambient){
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_ambient && az!=az_unknown)
                throw hiba(fvnev,"ambient is missing in line %u in %s",i,fileName.c_str());
        }
        ambiT=nulla;
        dbl dummy;
        if(mezo==FieldEl?false
            :mezo==FieldTherm?!fajl.lines()[i][1].todouble(ambiT)
            :mezo==FieldElTherm?!fajl.lines()[i][1].to2double(dummy,ambiT,true)
            :true)
            throw hiba(fvnev,"incorrect ambient value (%s) in line %u in %s",fajl.lines()[i][1].c_str(),i,fileName.c_str());
        if(fajl.lines()[i].size()>2){
            // több külsõ hõmérsékleten kell futtatni
            mulAmbiT.resize(fajl.lines()[i].size()-2);
            for(uns j=0; j<mulAmbiT.size(); j++){
                if(!fajl.lines()[i][j+2].todouble(mulAmbiT[j]))
                    throw hiba(fvnev,"incorrect ambient value (%s) in line %u in %s",fajl.lines()[i][1].c_str(),i,fileName.c_str());
            }
        }

        // powermap

        i++;
        az = azonosit(fajl.lines()[i][0].LowCase());
        while (az == az_unknown || az == az_powermap) {
            if (az == az_powermap) {
                pmap.init(this);
                pmap.read(path, fajl.lines()[i][1]);
                pmap.rescale();
            }
            i++;
            az = azonosit(fajl.lines()[i][0].LowCase());
        }

        // uchannel

        //      tucha tömb inicializálása

        uns uchadb = 0, ix;
        for (ix = 0; ix < colmax; ix++)
            if (pmodel->tcolor[ix].is && pmodel->tcolor[ix].tipus == SzinUchannel)
                uchadb++;

        tucha.resize(uchadb);

        for (ix = uchadb = 0; ix < colmax; ix++)
            if (pmodel->tcolor[ix].is && pmodel->tcolor[ix].tipus == SzinUchannel){
                tucha[uchadb].init1(this, ix, &pmodel->tcolor[ix]);
                uchadb++;
            }

        //      uchannel fájl nevének beolvasása és tucha tömb elemeinek inicializálása

        uchadb = 0;
        while (az == az_unknown || az == az_uchannel){
            if (az == az_uchannel){
                const PLString cimke = fajl.lines()[i][1].LowCase();
                uns j;
                for (j = 0; j < tucha.size() && tucha[j].getLabel() != cimke; j++);
                if (j == tucha.size())
                    throw hiba(fvnev, "unknown uchannel label (%s) in line %u in %s", cimke.c_str(), i, fileName.c_str());
                tucha[j].init2(path, fajl.lines()[i][2], cimke);
                uchadb++;
            }
            i++;
            az = azonosit(fajl.lines()[i][0].LowCase());
        }
        if (uchadb != tucha.size())
            throw hiba(fvnev, "missing uchannel definition (%u!=%u) in line %u in %s", uchadb, tucha.size(), i, fileName.c_str());

        // probe

        while(az!=az_probe && az!=az_fimtxt && az != az_commas && az != az_FIM_res && az != az_FIM_diff && az!=az_nofim 
            && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps 
            && az != az_analy && az!=az_end && az != az_no_plus_step_data){
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_probe && az!=az_analy && az!=az_fimtxt && az != az_commas && az != az_FIM_res 
                && az != az_FIM_diff && az!=az_nofim && az!=az_end && az!=az_unknown && az != az_no_plus_step_data
                && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps)
                throw hiba(fvnev,"probe or map_to_txt or FIM_res or FIM_diff or use_commas or no_images or auto_transi_steps or analysis or end is missing in line %u in %s",i,fileName.c_str());
        }
        tproV.clear();
        tproT.clear();
        while(az==az_probe){
            const PLString kisb=fajl.lines()[i][1].LowCase();
            ProbeTipus pr;
            if(kisb=="voltage")pr=ProbeV;
            else if(kisb=="temperature")pr=ProbeT;
            else if(kisb=="current")pr=ProbeC;
            else if(kisb=="map")pr=ProbeM;
            else if(kisb=="section")pr=ProbeS;
            else throw hiba(fvnev,"unknown probe type (%s) in line %u in %s",fajl.lines()[i][1].c_str(),i,fileName.c_str());
            uns aktindex=0;
            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_pro && az!=az_probe && az!=az_analy && az!=az_fimtxt && az != az_commas && az != az_FIM_res 
                    && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps
                    && az != az_FIM_diff && az!=az_nofim && az!=az_end && az!=az_unknown && az != az_no_plus_step_data)
                    throw hiba(fvnev,"pro or probe or auto_transi_steps or analysis or end is missing (%s found)in line %u in %s",fajl.lines()[i][0].c_str(),i,fileName.c_str());
                if(az==az_pro){
                    if(pr==ProbeV)tproV.resize((aktindex=tproV.size())+1);
                    else if(pr==ProbeT)tproT.resize((aktindex=tproT.size())+1);
                    else if(pr==ProbeC)tproC.resize((aktindex=tproC.size())+1);
                    else if(pr==ProbeM)tproM.resize((aktindex=tproM.size())+1);
                    else if(pr==ProbeS)tproS.resize((aktindex=tproS.size())+1);
                    if(pr==ProbeV)tproV[aktindex].cimke=fajl.lines()[i][1].LowCase();
                    else if(pr==ProbeT)tproT[aktindex].cimke=fajl.lines()[i][1].LowCase();
                    else if(pr==ProbeC)tproC[aktindex].cimke=fajl.lines()[i][1].LowCase();
                    else if(pr==ProbeM)tproM[aktindex].cimke=fajl.lines()[i][1].LowCase();
                    else if(pr==ProbeS)tproS[aktindex].cimke=fajl.lines()[i][1].LowCase();
                    PLString aktpar=fajl.lines()[i][2].LowCase();
                    if(pr==ProbeV||pr==ProbeT||pr==ProbeC||pr==ProbeM){
                        Oldal dal;
                        if(aktpar=="center")dal=CENTER;
                        else if(aktpar=="west")dal=WEST;
                        else if(aktpar=="east")dal=EAST;
                        else if(aktpar=="south")dal=SOUTH;
                        else if(aktpar=="north")dal=NORTH;
                        else if(aktpar=="top")dal=TOP;
                        else if(aktpar=="bottom")dal=BOTTOM;
                        else throw hiba(fvnev,"unknown probe place (%s) in line %u in %s (center, west, etc. needed)",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                        if(pr==ProbeV){
                            tproV[aktindex].oldal=dal;
                            if(!fajl.lines()[i][3].to3uns(tproV[aktindex].x,tproV[aktindex].y,tproV[aktindex].z))
                                throw hiba(fvnev,"bad probe valuse (%s) in line %u in %s",fajl.lines()[i][3].c_str(),i,fileName.c_str());
                        }
                        else if(pr==ProbeT){
                            tproT[aktindex].oldal=dal;
                            if(!fajl.lines()[i][3].to3uns(tproT[aktindex].x,tproT[aktindex].y,tproT[aktindex].z))
                                throw hiba(fvnev,"bad probe valuse (%s) in line %u in %s",fajl.lines()[i][3].c_str(),i,fileName.c_str());
                        }
                        else if(pr==ProbeC){
                            tproC[aktindex].oldal=dal;
                            if(!fajl.lines()[i][3].to3uns(tproC[aktindex].x,tproC[aktindex].y,tproC[aktindex].z))
                                throw hiba(fvnev,"bad probe valuse (%s) in line %u in %s",fajl.lines()[i][3].c_str(),i,fileName.c_str());
                            if(!fajl.lines()[i][4].to3uns(tproC[aktindex].x2,tproC[aktindex].y2,tproC[aktindex].z2))
                                throw hiba(fvnev,"bad probe valuse (%s) in line %u in %s",fajl.lines()[i][4].c_str(),i,fileName.c_str());
                        }
                        else if(pr==ProbeM){
                            tproM[aktindex].oldal=dal;
                            if(fajl.lines()[i][3].LowCase()=="x")tproM[aktindex].x=0;
                            else if(fajl.lines()[i][3].LowCase()=="y")tproM[aktindex].x=1;
                            else if(fajl.lines()[i][3].LowCase()=="z")tproM[aktindex].x=2;
                            else throw hiba(fvnev,"bad map probe direction (%s) in line %u in %s",fajl.lines()[i][3].c_str(),i,fileName.c_str());
                            if(!fajl.lines()[i][4].tounsigned(tproM[aktindex].y))
                                throw hiba(fvnev,"bad map probe valuse (%s) in line %u in %s",fajl.lines()[i][4].c_str(),i,fileName.c_str());
                        }
                    }
                    else if(pr==ProbeS){
                        if(aktpar=="x")     tproS[aktindex].oldal=WEST;
                        else if(aktpar=="y")tproS[aktindex].oldal=SOUTH;
                        else if(aktpar=="z")tproS[aktindex].oldal=BOTTOM;
                        else throw hiba(fvnev,"bad section probe type (%s) in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
                        if(!fajl.lines()[i][3].to2uns(tproS[aktindex].x,tproS[aktindex].y)) // mindig az x-y-ba teszi a koordinátákat, akkor is, ha nem az!
                            throw hiba(fvnev,"bad probe valuse (%s) in line %u in %s",fajl.lines()[i][3].c_str(),i,fileName.c_str());
                    }
                }
            }while(az==az_unknown || az==az_pro);
        }

        // map_to_txt

        while(az!=az_fimtxt && az!=az_commas && az != az_FIM_res && az != az_FIM_diff && az!=az_nofim && az!=az_analy 
            && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps 
            && az!=az_end && az != az_no_plus_step_data){
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_fimtxt && az != az_commas && az != az_FIM_res && az != az_FIM_diff && az!=az_nofim 
                && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps 
                && az!=az_analy && az!=az_end && az!=az_unknown && az != az_no_plus_step_data)
                throw hiba(fvnev,"map_to_txt or FIM_res or FIM_diff or use_commas or no_images or auto_transi_steps or analysis or end is missing in line %u in %s",i,fileName.c_str());
        }
        if(az==az_fimtxt)is_fim_txt=true;

        // FIM_res

        FIM_res_xy = FIM_res_z = 0;
        while (az != az_commas && az != az_FIM_res && az != az_FIM_diff && az != az_nofim && az != az_analy 
            && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps 
            && az != az_end && az != az_no_plus_step_data) {
            i++;
            az = azonosit(fajl.lines()[i][0].LowCase());
            if (az != az_commas && az != az_FIM_res && az != az_FIM_diff && az != az_nofim && az != az_analy && az != az_no_plus_step_data
                && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps && az != az_end && az != az_unknown)
                throw hiba(fvnev, "map_to_txt or FIM_res or FIM_diff or use_commas or no_images or auto_transi_steps or analysis or end is missing in line %u in %s", i, fileName.c_str());
        }
        if (az == az_FIM_res) {
            if (!fajl.lines()[i][1].tounsigned(FIM_res_xy))
                throw hiba(fvnev, "FIM_res_xy is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
            if (!fajl.lines()[i][2].tounsigned(FIM_res_z))
                throw hiba(fvnev, "FIM_res_z is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
        }

        // FIM_diff

        FIM_diff_x = FIM_diff_y = FIM_diff_z = 0;
        while (az != az_commas && az != az_FIM_diff && az != az_nofim && az != az_analy && az != az_end
            && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps
            && az != az_no_plus_step_data) {
            i++;
            az = azonosit(fajl.lines()[i][0].LowCase());
            if (az != az_commas && az != az_FIM_diff && az != az_nofim && az != az_analy && az != az_end && az != az_no_plus_step_data
                && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps && az != az_unknown)
                throw hiba(fvnev, "map_to_txt or FIM_diff or use_commas or no_images or auto_transi_steps or analysis or end is missing in line %u in %s", i, fileName.c_str());
        }
        if (az == az_FIM_diff) {
            if (!fajl.lines()[i][1].tounsigned(FIM_diff_x))
                throw hiba(fvnev, "FIM_diff_x is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
            if (!fajl.lines()[i][2].tounsigned(FIM_diff_y))
                throw hiba(fvnev, "FIM_diff_y is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
            if (!fajl.lines()[i][3].tounsigned(FIM_diff_z))
                throw hiba(fvnev, "FIM_diff_z is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
        }

        // use_commas

        while(az!=az_commas && az!=az_nofim && az!=az_analy && az!=az_end && az != az_no_plus_step_data
            && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps){
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_commas && az!=az_nofim && az!=az_analy && az!=az_end && az!=az_unknown && az != az_no_plus_step_data
                && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps)
                throw hiba(fvnev,"use_commas or no_images or auto_transi_steps or analysis or end is missing in line %u in %s",i,fileName.c_str());
        }
		if(az==az_commas)is_vesszo=true;

        // no_images

        while(az!=az_nofim && az!=az_analy && az!=az_end && az != az_auto_transi_steps_V 
            && az != az_auto_transi_steps_T && az != az_auto_transi_steps && az != az_no_plus_step_data){
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_nofim && az!=az_analy && az!=az_end && az!=az_unknown && az != az_no_plus_step_data
                && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps)
                throw hiba(fvnev,"no_images or auto_transi_steps or analysis or end is missing in line %u in %s",i,fileName.c_str());
        }
		if(az==az_nofim)is_nofim=true;

        // no_plus_step_data

        while (az != az_analy && az != az_end && az != az_auto_transi_steps_V
            && az != az_auto_transi_steps_T && az != az_auto_transi_steps && az != az_no_plus_step_data) {
            i++;
            az = azonosit(fajl.lines()[i][0].LowCase());
            if (az != az_analy && az != az_end && az != az_unknown && az != az_no_plus_step_data
                && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps)
                throw hiba(fvnev, "no_images or auto_transi_steps or analysis or end is missing in line %u in %s", i, fileName.c_str());
        }
        if (az == az_no_plus_step_data)auto_tra.is_no_plus_step_data = true;

        //  az_auto_transi_steps

        while (az != az_analy && az != az_end && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps) {
            i++;
            az = azonosit(fajl.lines()[i][0].LowCase());
            if (az != az_analy && az != az_end && az != az_unknown
                && az != az_auto_transi_steps_V && az != az_auto_transi_steps_T && az != az_auto_transi_steps)
                throw hiba(fvnev, "auto_transi_steps or analysis or end is missing in line %u in %s", i, fileName.c_str());
        }
        if (az == az_auto_transi_steps) {
            if (auto_tra.is_T || auto_tra.is_V)
                throw hiba(fvnev, "auto_transi_steps_T/auto_transi_steps_V/auto_transi_steps redefinition with auto_transi_steps in line %u in %s", i, fileName.c_str());
            if (!fajl.lines()[i][1].todouble(auto_tra.V_max_rel))
                throw hiba(fvnev, "auto_transi_steps par1 is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
            if (!fajl.lines()[i][2].todouble(auto_tra.V_min_dt))
                throw hiba(fvnev, "auto_transi_steps par2 is not a number (%s) in %s", fajl.lines()[i][2].c_str(), fileName.c_str());
            if (!fajl.lines()[i][3].tounsigned(auto_tra.V_max_plus_steps))
                throw hiba(fvnev, "auto_transi_steps par3 is not a number (%s) in %s", fajl.lines()[i][3].c_str(), fileName.c_str());
            auto_tra.T_max_rel = auto_tra.V_max_rel;
            auto_tra.T_min_dt = auto_tra.V_min_dt;
            auto_tra.T_max_plus_steps = auto_tra.V_max_plus_steps;
            auto_tra.is_T = auto_tra.is_V = true;
        }
        else  if (az == az_auto_transi_steps_V) {
            if (auto_tra.is_V)
                throw hiba(fvnev, "auto_transi_steps_V/auto_transi_steps redefinition with auto_transi_steps_V in line %u in %s", i, fileName.c_str());
            if (!fajl.lines()[i][1].todouble(auto_tra.V_max_rel))
                throw hiba(fvnev, "auto_transi_steps_V par1 is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
            if (!fajl.lines()[i][2].todouble(auto_tra.V_min_dt))
                throw hiba(fvnev, "auto_transi_steps_V par2 is not a number (%s) in %s", fajl.lines()[i][2].c_str(), fileName.c_str());
            if (!fajl.lines()[i][3].tounsigned(auto_tra.V_max_plus_steps))
                throw hiba(fvnev, "auto_transi_steps_V par3 is not a number (%s) in %s", fajl.lines()[i][3].c_str(), fileName.c_str());
            auto_tra.is_V = true;
        }
        else  if (az == az_auto_transi_steps_T) {
            if (auto_tra.is_T)
                throw hiba(fvnev, "auto_transi_steps_T/auto_transi_steps redefinition with auto_transi_steps_T in line %u in %s", i, fileName.c_str());
            if (!fajl.lines()[i][1].todouble(auto_tra.T_max_rel))
                throw hiba(fvnev, "auto_transi_steps_T par1 is not a number (%s) in %s", fajl.lines()[i][1].c_str(), fileName.c_str());
            if (!fajl.lines()[i][2].todouble(auto_tra.T_min_dt))
                throw hiba(fvnev, "auto_transi_steps_T par2 is not a number (%s) in %s", fajl.lines()[i][2].c_str(), fileName.c_str());
            if (!fajl.lines()[i][3].tounsigned(auto_tra.T_max_plus_steps))
                throw hiba(fvnev, "auto_transi_steps_T par3 is not a number (%s) in %s", fajl.lines()[i][3].c_str(), fileName.c_str());
            auto_tra.is_T = true;
        }

        
        // analysis

        while(az!=az_analy && az!=az_end){
            i++;
            az=azonosit(fajl.lines()[i][0].LowCase());
            if(az!=az_analy && az!=az_end && az!=az_unknown)
                throw hiba(fvnev,"analysis or end is missing (%s found) in line %u in %s", fajl.lines()[i][0].c_str(), i,fileName.c_str());
        }
        tanal.clear();
        while(az==az_analy){
            tanal.resize(tanal.size()+1);
            analysis & aktanal=tanal[tanal.size()-1];
            aktanal.nev=fajl.lines()[i][1];
            const PLString kisb=fajl.lines()[i][2].LowCase();
            if(kisb=="dc")aktanal.tipus=AnalDC;
            else if(kisb=="ndc")aktanal.tipus=AnalNDC;
            else if(kisb=="ac")aktanal.tipus=AnalAC;
            else if(kisb=="lintransient")aktanal.tipus=AnalLinTran;
            else if(kisb=="logtransient")aktanal.tipus=AnalLogTran;
            else if(kisb=="bode")aktanal.tipus=AnalBode;
            else if(kisb=="timeconst")aktanal.tipus=AnalIdo;
            else if(kisb=="controlled")aktanal.tipus=AnalCtrl;
            else throw hiba(fvnev,"bad analysis type (%s) in line %u in %s",fajl.lines()[i][2].c_str(),i,fileName.c_str());
            if(aktanal.tipus==AnalAC){
                if(!fajl.lines()[i][3].todouble(aktanal.from))
                    throw hiba(fvnev,"bad analysis format (%s) in line %u in %s",fajl.lines()[i][3].c_str(),i,fileName.c_str());
            }
            else if(aktanal.tipus==AnalLinTran){
                if(!fajl.lines()[i][3].to2double(aktanal.from,aktanal.to))
                    throw hiba(fvnev,"bad analysis format (%s) in line %u in %s",fajl.lines()[i][3].c_str(),i,fileName.c_str());
            }
            else if(aktanal.tipus==AnalLogTran||aktanal.tipus==AnalBode||aktanal.tipus==AnalIdo){
                if(!fajl.lines()[i][3].to2double1uns(aktanal.from,aktanal.to,aktanal.step))
                    throw hiba(fvnev,"bad analysis format (%s) in line %u in %s",fajl.lines()[i][3].c_str(),i,fileName.c_str());
            }
            else if(aktanal.tipus==AnalNDC){
                if(fajl.lines()[i].size()>3&&!fajl.lines()[i][3].to1uns2double(aktanal.ndc_maxiter,aktanal.relhiba,aktanal.ndc_I0))
                    throw hiba(fvnev,"bad analysis format (%s) in line %u in %s",fajl.lines()[i][3].c_str(),i,fileName.c_str());
                // aktanal.ndc_I0 /= (pmodel->A_semi_full==nulla) ? 1.0 : pmodel->A_semi_full; // az áram az egész félvezetõ felületre érvényes
            }
            else if(aktanal.tipus==AnalCtrl){
                read_ctrl(path+fajl.lines()[i][3],aktanal,pmodel,*this);
            }

            do{
                i++;
                az=azonosit(fajl.lines()[i][0].LowCase());
                if(az!=az_analy && az!=az_end && az!=az_unknown)
                    throw hiba(fvnev,"analysis or end is missing in line %u in %s",i,fileName.c_str());
            }while(az!=az_analy && az!=az_end);
        }
    }
    catch(const hiba & h){
        PLString s=h.what();
        if(s.find("tomb::operator[]")==npos)throw;
        throw hiba(fvnev,"missing parameters or unexpected end of file in line %u in %s",i,fileName.c_str());
    }
}


//***********************************************************************
bool red_fa::szukito(MezoTipus mt, const tomb3d<t_modell_cella> & r, doboz & d){
// Az mt a szimuláció típusát jelenti, de elth esetén a konkrét cella lehet csak el vagy th!
// return: false, ha nincs cella a térfogatban
//***********************************************************************
    bool kell_el = mt == FieldEl || mt == FieldElTherm;
    bool kell_th = mt == FieldTherm || mt == FieldElTherm;
    uns x, y, z;

    // WEST

    for (x = d.x1; x <= d.x2; x++) {
        for (y = d.y1; y <= d.y2; y++)
            for (z = d.z1; z <= d.z2; z++) {
                const t_modell_cella & c = r.getconstref(x, y, z);
                if (c.is_cella && ((c.is_el && kell_el) || (c.is_th && kell_th)))
                    goto ki_west;
            }
    }
ki_west:
    if (x > d.x2)
        return false;
    d.x1 = x;

    // EAST

    for (x = d.x2; x >= d.x1; x--) {
        for (y = d.y1; y <= d.y2; y++)
            for (z = d.z1; z <= d.z2; z++) {
                const t_modell_cella & c = r.getconstref(x, y, z);
                if (c.is_cella && ((c.is_el && kell_el) || (c.is_th && kell_th)))
                    goto ki_east;
            }
    }
ki_east:
    d.x2 = x;

    // SOUTH

    for (y = d.y1; y <= d.y2; y++) {
        for (x = d.x1; x <= d.x2; x++)
            for (z = d.z1; z <= d.z2; z++) {
                const t_modell_cella & c = r.getconstref(x, y, z);
                if (c.is_cella && ((c.is_el && kell_el) || (c.is_th && kell_th)))
                    goto ki_south;
            }
    }
ki_south:
    d.y1 = y;

    // NORTH

    for (y = d.y2; y >= d.y1; y--) {
        for (x = d.x1; x <= d.x2; x++)
            for (z = d.z1; z <= d.z2; z++) {
                const t_modell_cella & c = r.getconstref(x, y, z);
                if (c.is_cella && ((c.is_el && kell_el) || (c.is_th && kell_th)))
                    goto ki_north;
            }
    }
ki_north:
    d.y2 = (uns)y;

    // BOTTOM

    for (z = d.z1; z <= d.z2; z++) {
        for (y = d.y1; y <= d.y2; y++)
            for (x = d.x1; x <= d.x2; x++) {
                const t_modell_cella & c = r.getconstref(x, y, z);
                if (c.is_cella && ((c.is_el && kell_el) || (c.is_th && kell_th)))
                    goto ki_bottom;
            }
    }
ki_bottom:
    d.z1 = z;

    // TOP

    for (z = d.z2; z >= d.z1; z--) {
        for (y = d.y1; y <= d.y2; y++)
            for (x = d.x1; x <= d.x2; x++) {
                const t_modell_cella & c = r.getconstref(x, y, z);
                if (c.is_cella && ((c.is_el && kell_el) || (c.is_th && kell_th)))
                    goto ki_top;
            }
    }
ki_top:
    d.z2 = z;
    return true;
}


//***********************************************************************
void red_fa::iranyszamolo_uj(const doboz & be_meret, tomb<Oldal>& iranyok){
// Hogy jó sorrendben legyenek a face-ek a redukció során, azonos szinten mindig
// ugyanolyan irányban történik az összevonás.
// Azonos méret esetén a z irányú összevonás a preferált, mert ott általában azonos az anyag.
//***********************************************************************
    dbl dx = be_meret.x2 - be_meret.x1 + 1;
    dbl dy = be_meret.y2 - be_meret.y1 + 1;
    dbl dz = be_meret.z2 - be_meret.z1 + 1;
    iranyok.clear();
    while (dx >= 1 || dy >= 1 || dz >= 1) {
        if (dx > dy) {
            if (dx > dz) { // dx a legnagyobb
                iranyok.add(EAST);
                dx /= 2;
            }
            else { // dz a legnagyobb
                iranyok.add(TOP);
                dz /= 2;
            }
        }
        else if (dy > dz) { // dy a legnagyobb
            iranyok.add(NORTH);
            dy /= 2;
        }
        else { // dz a legnagyobb
            iranyok.add(TOP);
            dz /= 2;
        }
    }
}


//***********************************************************************
uns red_fa::set_indexek(red_fa * gy, uns & start_index){
//***********************************************************************
    if (gy == nullptr)
        return 0;
    gy->blokk_kezdo_index = start_index;
    uns db = 1 + set_indexek(gy->bal, start_index);
    db += set_indexek(gy->jobb, start_index);
    gy->index = start_index;
    start_index++;
    return db;
}

//***********************************************************************
void red_fa::write_tree(FILE * fp, red_fa * gy){
//***********************************************************************
    if (gy == nullptr)
        return;
    write_tree(fp, gy->bal);
    write_tree(fp, gy->jobb);
    if (gy->bal == nullptr) { // elemi cella
        if (gy->cella == nullptr)
            throw hiba("red_fa::write_tree", "gy->cella == nullptr");
        fprintf(fp, "L%uC%uA%uB%u\n", gy->index, gy->cella->cella_index, gy->A0, gy->B0);
    }
    else{ // belsõ csomópont
        fprintf(fp, "BF%uF%uL%uR%u\n", gy->index, gy->blokk_kezdo_index, gy->bal->index, gy->jobb->index);
        fprintf(fp, "NM%uA%uB%u\n", gy->mit_hova.size() + 1, gy->mit_hova.size() > 0 ? gy->mit_hova.getLast().hova + gy->mit_hova.getLast().hanyat : 0, gy->kozos_db); // A közös miatt a +1
        fprintf(fp, "B0E%uK%uN%u\n", gy->kozos_be1, gy->kozos_be2, gy->kozos_db);
        for (uns i = 0; i < gy->mit_hova.size(); i++)
            fprintf(fp, "A%u%c%uN%u\n", gy->mit_hova[i].hova, (gy->mit_hova[i].is_be1 ? 'E' : 'K'), gy->mit_hova[i].honnan, gy->mit_hova[i].hanyat);
        fprintf(fp, "EF%u\n", gy->index);
    }
}


//***********************************************************************
inline void red_fa::optimized_add_mit_hova(const mit_hova_masol & uj){
//***********************************************************************
    if (mit_hova.size() == 0) {
        mit_hova.add(uj);
        return;
    }
    mit_hova_masol & elozo = mit_hova.getLast();
    if (elozo.is_be1 == uj.is_be1 && elozo.honnan + elozo.hanyat == uj.honnan && elozo.hova + elozo.hanyat == uj.hova) {
        elozo.hanyat += uj.hanyat;
    }
    else {
        mit_hova.add(uj);
    }
}


//***********************************************************************
red_fa * red_fa::build_tree_uj(MezoTipus mt, const tomb3d<t_modell_cella> & r, doboz befoglalo_d, const tomb<Oldal> & iranyok, uns level) {
// Az mt a szimuláció típusát jelenti, de elth esetén a konkrét cella lehet csak el vagy th!
// A face-ek sorrendjét jól generálja=>adott szinten azonos irányú összevonás, mindig ugyanannál a koordinátánál kell legyen a felosztás, így
// a befoglaló doboz alapján számoljuk, nem a tányleges doboz alapján.
//***********************************************************************
    red_fa * gy = nullptr;
    doboz dm = befoglalo_d;
    if (!szukito(mt, r, dm))
        return nullptr; // Ha nincs cella a térfogatban, nullptr-t ad vissza.
    if (dm.x1 != dm.x2 || dm.y1 != dm.y2 || dm.z1 != dm.z2) { // tovább kell bontani
        gy = new red_fa(dm, level);
        doboz ki_1 = befoglalo_d, ki_2 = befoglalo_d;
        gy->kozos_oldal_1 = iranyok[level];
        switch (gy->kozos_oldal_1) {
            case EAST:
                ki_1.x2 = (befoglalo_d.x1 + befoglalo_d.x2) / 2;
                ki_2.x1 = (befoglalo_d.x1 + befoglalo_d.x2) / 2 + 1;
                break;
            case NORTH:
                ki_1.y2 = (befoglalo_d.y1 + befoglalo_d.y2) / 2;
                ki_2.y1 = (befoglalo_d.y1 + befoglalo_d.y2) / 2 + 1;
                break;
            case TOP:
                ki_1.z2 = (befoglalo_d.z1 + befoglalo_d.z2) / 2;
                ki_2.z1 = (befoglalo_d.z1 + befoglalo_d.z2) / 2 + 1;
                break;
            default:
                throw hiba("red_fa::build_tree_uj", "impossible reduction direction");
        }
        gy->bal = build_tree_uj(mt, r, ki_1, iranyok, level + 1);
        gy->jobb = build_tree_uj(mt, r, ki_2, iranyok, level + 1);
    }
    else { // egy cellánk van: felbontás alcellákra
        const t_modell_cella & c = r.getconstref(dm.x1, dm.y1, dm.z1);
        if (c.belso_cellak.size() == 0) {
            gy = new red_fa(dm, level);
            gy->cella = &c;
        }
        else {
            doboz belso_doboz(0, c.belso_cellak.x_size() - 1, 0, c.belso_cellak.y_size() - 1, 0, c.belso_cellak.z_size() - 1);
            tomb<Oldal> belso_iranyok;
            iranyszamolo_uj(belso_doboz, belso_iranyok);
            gy = build_subtree_uj(mt, c, belso_doboz, belso_iranyok, level);
            optimize_tree(gy);
        }
    }
    return gy;
}


//***********************************************************************
red_fa * red_fa::build_tree_old(MezoTipus mt, const tomb3d<t_modell_cella> & r, doboz d, uns level) {
// Az mt a szimuláció típusát jelenti, de elth esetén a konkrét cella lehet csak el vagy th!
//***********************************************************************
    red_fa * gy = nullptr;
    if (!szukito(mt, r, d))
        return nullptr; // Ha nincs cella a térfogatban, nullptr-t ad vissza.
    if (d.x1 != d.x2 || d.y1 != d.y2 || d.z1 != d.z2) { // tovább kell bontani
        gy = new red_fa(d, level);
        doboz ki_1, ki_2;
        feloszt(mt, r, d, ki_1, ki_2, gy->kozos_oldal_1);
        gy->bal = build_tree_old(mt, r, ki_1, level + 1);
        gy->jobb = build_tree_old(mt, r, ki_2, level + 1);
    }
    else { // egy cellánk van: felbontás alcellákra
        const t_modell_cella & c = r.getconstref(d.x1, d.y1, d.z1);
        if (c.belso_cellak.size() == 0) {
            gy = new red_fa(d, level);
            gy->cella = &c;
        }
        else {
            gy = build_subtree_old(mt, c, doboz(0, c.belso_cellak.x_size() - 1, 0, c.belso_cellak.y_size() - 1, 0, c.belso_cellak.z_size() - 1), level);
        }
    }
    return gy;
}


//***********************************************************************
red_fa * red_fa::build_subtree_old(MezoTipus mt, const t_modell_cella & c, doboz d, uns level){
// Az mt a szimuláció típusát jelenti, de elth esetén a konkrét cella lehet csak el vagy th!
//***********************************************************************
    red_fa * gy = new red_fa(d, level);
    gy->is_alcella = true;
    if (d.x1 != d.x2 || d.y1 != d.y2 || d.z1 != d.z2) { // tovább kell bontani
        doboz ki_1, ki_2;
        feloszt(mt, c.belso_cellak, d, ki_1, ki_2, gy->kozos_oldal_1);
        gy->bal = build_subtree_old(mt, c, ki_1, level + 1);
        gy->jobb = build_subtree_old(mt, c, ki_2, level + 1);
    }
    else { // elemi cella
        gy->cella = &c.belso_cellak.getconstref(d.x1, d.y1, d.z1);
    }
    return gy;
}


//***********************************************************************
red_fa * red_fa::build_subtree_uj(MezoTipus mt, const t_modell_cella & c, doboz befoglalo_d, const tomb<Oldal> & iranyok, uns level) {
// Az mt a szimuláció típusát jelenti, de elth esetén a konkrét cella lehet csak el vagy th!
//***********************************************************************
    red_fa * gy = new red_fa(befoglalo_d, level);
    gy->is_alcella = true;
    if (befoglalo_d.x1 != befoglalo_d.x2 || befoglalo_d.y1 != befoglalo_d.y2 || befoglalo_d.z1 != befoglalo_d.z2) { // tovább kell bontani
        doboz ki_1 = befoglalo_d, ki_2 = befoglalo_d;
        gy->kozos_oldal_1 = iranyok[level];
        switch (gy->kozos_oldal_1) {
            case EAST:
                ki_1.x2 = (befoglalo_d.x1 + befoglalo_d.x2) / 2;
                ki_2.x1 = (befoglalo_d.x1 + befoglalo_d.x2) / 2 + 1;
                break;
            case NORTH:
                ki_1.y2 = (befoglalo_d.y1 + befoglalo_d.y2) / 2;
                ki_2.y1 = (befoglalo_d.y1 + befoglalo_d.y2) / 2 + 1;
                break;
            case TOP:
                ki_1.z2 = (befoglalo_d.z1 + befoglalo_d.z2) / 2;
                ki_2.z1 = (befoglalo_d.z1 + befoglalo_d.z2) / 2 + 1;
                break;
            default:
                throw hiba("red_fa::build_subtree_uj", "impossible reduction direction");
        }
        gy->bal = build_subtree_uj(mt, c, ki_1, iranyok, level + 1);
        gy->jobb = build_subtree_uj(mt, c, ki_2, iranyok, level + 1);
    }
    else { // elemi cella
        gy->cella = &c.belso_cellak.getconstref(befoglalo_d.x1, befoglalo_d.y1, befoglalo_d.z1);
    }
    return gy;
}


//***********************************************************************
void red_fa::decrease_level(red_fa * gy){
//***********************************************************************
    if (gy == nullptr)
        return;
    decrease_level(gy->bal);
    decrease_level(gy->jobb);
    gy->level--;
}


//***********************************************************************
void red_fa::optimize_tree(red_fa *& gy){
//***********************************************************************
    if (gy == nullptr)
        return;
    if (gy->bal == nullptr && gy->jobb == nullptr)
        return;
    optimize_tree(gy->bal);
    optimize_tree(gy->jobb);
    if (gy->bal == nullptr) {
        // printf("opt: bal %u\n",gy->level);
        decrease_level(gy->jobb);
        red_fa *temp = gy;
        gy = gy->jobb;
        delete temp;
        return;
    }
    if (gy->jobb == nullptr) {
        // printf("opt: jobb %u\n", gy->level);
        decrease_level(gy->bal);
        red_fa *temp = gy;
        gy = gy->bal;
        delete temp;
        return;
    }
}


//***********************************************************************
void red_fa::tol_ig_db_szamol(MezoTipus mt, red_fa * gy){
// Az mt a szimuláció típusát jelenti, de elth esetén a konkrét cella lehet csak el vagy th!
//***********************************************************************
    if (gy == nullptr)
        return;
    if (gy->bal == nullptr || gy->jobb == nullptr) { // elemi cella
        if (gy->bal != gy->jobb)
            throw hiba("red_fa::tol_ig_db_szamol", "impossible tree (bal or jobb is nullptr not both)");
        if (gy->cella == nullptr)
            throw hiba("red_fa::tol_ig_db_szamol", "impossible tree (cella is nullptr)");
        if (gy->cella->belso_cellak.size()!=0)
            throw hiba("red_fa::tol_ig_db_szamol", "impossible tree (belso_cellak.size()!=0)");
        uns tol = 0;
        for (uns i = 1; i < BASIC_SIDES; i++) { // TODO: ha van external, azt kezelni kéne, most 0-nak tekinti
            switch (mt) {
                case FieldEl:
                    if (!gy->cella->is_el)
                        throw hiba("red_fa::tol_ig_db_szamol", "FieldEl: !gy->cella->is_el");
                    gy->oldalak[i].tol = tol;
                    tol += gy->oldalak[i].db = gy->cella->face_adat[i].kulso_el_db;
                    break;
                case FieldTherm:
                    if (!gy->cella->is_th)
                        throw hiba("red_fa::tol_ig_db_szamol", "FieldTherm: !gy->cella->is_th");
                    gy->oldalak[i].tol = tol;
                    tol += gy->oldalak[i].db = gy->cella->face_adat[i].kulso_th_db;
                    break;
                case FieldElTherm:
                    if (!gy->cella->is_el && gy->cella->face_adat[i].kulso_el_db > 0)
                        throw hiba("red_fa::tol_ig_db_szamol", "FieldElTherm: !gy->cella->is_el && gy->cella->face_adat[i].kulso_el_db > 0 (%u)", gy->cella->face_adat[i].kulso_el_db);
                    if (!gy->cella->is_th && gy->cella->face_adat[i].kulso_th_db > 0)
                        throw hiba("red_fa::tol_ig_db_szamol", "FieldElTherm: !gy->cella->is_th && gy->cella->face_adat[i].kulso_th_db > 0 (%u)", gy->cella->face_adat[i].kulso_th_db);
                    gy->oldalak[i].tol = tol;
                    tol += gy->oldalak[i].db = gy->cella->face_adat[i].kulso_el_db + gy->cella->face_adat[i].kulso_th_db;
                    break;
            }
        }
        gy->A0 = tol; // centroidnál ez ++
        gy->B0 = (gy->cella->is_el ? 1 : 0) + (gy->cella->is_th ? 1 : 0);// ez szerintem szar, mert a cella típusa számít, nem a szimulációé: mt == FieldElTherm ? 2 : 1; // centroidnál ez --
    }
    else { // összevonás
        tol_ig_db_szamol(mt, gy->bal);
        tol_ig_db_szamol(mt, gy->jobb);
        mit_hova_masol akt;
        uns hova = 0;

        // EXTERNAL

        if (gy->bal->oldalak[EXTERNAL].db != 0) {
            akt.is_be1 = true;
            akt.honnan = gy->bal->oldalak[EXTERNAL].tol;
            akt.hanyat = gy->bal->oldalak[EXTERNAL].db;
            akt.hova = hova;
            if (gy->oldalak[EXTERNAL].db == 0)
                gy->oldalak[EXTERNAL].tol = hova;
            gy->oldalak[EXTERNAL].db += akt.hanyat;
            hova += akt.hanyat;
            gy->optimized_add_mit_hova(akt);
        }
        if (gy->jobb->oldalak[EXTERNAL].db != 0) {
            akt.is_be1 = false;
            akt.honnan = gy->jobb->oldalak[EXTERNAL].tol;
            akt.hanyat = gy->jobb->oldalak[EXTERNAL].db;
            akt.hova = hova;
            if (gy->oldalak[EXTERNAL].db == 0)
                gy->oldalak[EXTERNAL].tol = hova;
            gy->oldalak[EXTERNAL].db += akt.hanyat;
            hova += akt.hanyat;
            gy->optimized_add_mit_hova(akt);
        }

        // WEST és EAST

        if (gy->bal->oldalak[WEST].db != 0) {
            akt.is_be1 = true;
            akt.honnan = gy->bal->oldalak[WEST].tol;
            akt.hanyat = gy->bal->oldalak[WEST].db;
            akt.hova = hova;
            if (gy->oldalak[WEST].db == 0)
                gy->oldalak[WEST].tol = hova;
            gy->oldalak[WEST].db += akt.hanyat;
            hova += akt.hanyat;
            gy->optimized_add_mit_hova(akt);
        }
        if (gy->kozos_oldal_1 == EAST) { // A bal cella EAST oldala a jobb cella WEST oldalához csatlakozik
            if(gy->bal->oldalak[EAST].db!=gy->jobb->oldalak[WEST].db)
                throw hiba("red_fa::tol_ig_db_szamol", "gy->bal->oldalak[EAST].db!=gy->jobb->oldalak[WEST].db (%u!=%u)", gy->bal->oldalak[EAST].db, gy->jobb->oldalak[WEST].db);
            gy->kozos_be1 = gy->bal->oldalak[EAST].tol;
            gy->kozos_be2 = gy->jobb->oldalak[WEST].tol;
            gy->kozos_db = gy->bal->oldalak[EAST].db;
        }
        else {
            if (gy->jobb->oldalak[WEST].db != 0) {
                akt.is_be1 = false;
                akt.honnan = gy->jobb->oldalak[WEST].tol;
                akt.hanyat = gy->jobb->oldalak[WEST].db;
                akt.hova = hova;
                if (gy->oldalak[WEST].db == 0)
                    gy->oldalak[WEST].tol = hova;
                gy->oldalak[WEST].db += akt.hanyat;
                hova += akt.hanyat;
                gy->optimized_add_mit_hova(akt);
            }
            if (gy->bal->oldalak[EAST].db != 0) {
                akt.is_be1 = true;
                akt.honnan = gy->bal->oldalak[EAST].tol;
                akt.hanyat = gy->bal->oldalak[EAST].db;
                akt.hova = hova;
                if (gy->oldalak[EAST].db == 0)
                    gy->oldalak[EAST].tol = hova;
                gy->oldalak[EAST].db += akt.hanyat;
                hova += akt.hanyat;
                gy->optimized_add_mit_hova(akt);
            }
        }
        if (gy->jobb->oldalak[EAST].db != 0) {
            akt.is_be1 = false;
            akt.honnan = gy->jobb->oldalak[EAST].tol;
            akt.hanyat = gy->jobb->oldalak[EAST].db;
            akt.hova = hova;
            if (gy->oldalak[EAST].db == 0)
                gy->oldalak[EAST].tol = hova;
            gy->oldalak[EAST].db += akt.hanyat;
            hova += akt.hanyat;
            gy->optimized_add_mit_hova(akt);
        }

        // SOUTH és NORTH

        if (gy->bal->oldalak[SOUTH].db != 0) {
            akt.is_be1 = true;
            akt.honnan = gy->bal->oldalak[SOUTH].tol;
            akt.hanyat = gy->bal->oldalak[SOUTH].db;
            akt.hova = hova;
            if (gy->oldalak[SOUTH].db == 0)
                gy->oldalak[SOUTH].tol = hova;
            gy->oldalak[SOUTH].db += akt.hanyat;
            hova += akt.hanyat;
            gy->optimized_add_mit_hova(akt);
        }
        if (gy->kozos_oldal_1 == NORTH) { // A bal cella NORTH oldala a jobb cella SOUTH oldalához csatlakozik
            if (gy->bal->oldalak[NORTH].db != gy->jobb->oldalak[SOUTH].db)
                throw hiba("red_fa::tol_ig_db_szamol", "gy->bal->oldalak[NORTH].db!=gy->jobb->oldalak[SOUTH].db (%u!=%u)", gy->bal->oldalak[NORTH].db, gy->jobb->oldalak[SOUTH].db);
            gy->kozos_be1 = gy->bal->oldalak[NORTH].tol;
            gy->kozos_be2 = gy->jobb->oldalak[SOUTH].tol;
            gy->kozos_db = gy->bal->oldalak[NORTH].db;
        }
        else {
            if (gy->jobb->oldalak[SOUTH].db != 0) {
                akt.is_be1 = false;
                akt.honnan = gy->jobb->oldalak[SOUTH].tol;
                akt.hanyat = gy->jobb->oldalak[SOUTH].db;
                akt.hova = hova;
                if (gy->oldalak[SOUTH].db == 0)
                    gy->oldalak[SOUTH].tol = hova;
                gy->oldalak[SOUTH].db += akt.hanyat;
                hova += akt.hanyat;
                gy->optimized_add_mit_hova(akt);
            }
            if (gy->bal->oldalak[NORTH].db != 0) {
                akt.is_be1 = true;
                akt.honnan = gy->bal->oldalak[NORTH].tol;
                akt.hanyat = gy->bal->oldalak[NORTH].db;
                akt.hova = hova;
                if (gy->oldalak[NORTH].db == 0)
                    gy->oldalak[NORTH].tol = hova;
                gy->oldalak[NORTH].db += akt.hanyat;
                hova += akt.hanyat;
                gy->optimized_add_mit_hova(akt);
            }
        }
        if (gy->jobb->oldalak[NORTH].db != 0) {
            akt.is_be1 = false;
            akt.honnan = gy->jobb->oldalak[NORTH].tol;
            akt.hanyat = gy->jobb->oldalak[NORTH].db;
            akt.hova = hova;
            if (gy->oldalak[NORTH].db == 0)
                gy->oldalak[NORTH].tol = hova;
            gy->oldalak[NORTH].db += akt.hanyat;
            hova += akt.hanyat;
            gy->optimized_add_mit_hova(akt);
        }

        // BOTTOM és TOP

        if (gy->bal->oldalak[BOTTOM].db != 0) {
            akt.is_be1 = true;
            akt.honnan = gy->bal->oldalak[BOTTOM].tol;
            akt.hanyat = gy->bal->oldalak[BOTTOM].db;
            akt.hova = hova;
            if (gy->oldalak[BOTTOM].db == 0)
                gy->oldalak[BOTTOM].tol = hova;
            gy->oldalak[BOTTOM].db += akt.hanyat;
            hova += akt.hanyat;
            gy->optimized_add_mit_hova(akt);
        }
        if (gy->kozos_oldal_1 == TOP) { // A bal cella TOP oldala a jobb cella BOTTOM oldalához csatlakozik
            if (gy->bal->oldalak[TOP].db != gy->jobb->oldalak[BOTTOM].db)
                throw hiba("red_fa::tol_ig_db_szamol", "gy->bal->oldalak[TOP].db!=gy->jobb->oldalak[BOTTOM].db (%u!=%u)", gy->bal->oldalak[TOP].db, gy->jobb->oldalak[BOTTOM].db);
            gy->kozos_be1 = gy->bal->oldalak[TOP].tol;
            gy->kozos_be2 = gy->jobb->oldalak[BOTTOM].tol;
            gy->kozos_db = gy->bal->oldalak[TOP].db;
        }
        else {
            if (gy->jobb->oldalak[BOTTOM].db != 0) {
                akt.is_be1 = false;
                akt.honnan = gy->jobb->oldalak[BOTTOM].tol;
                akt.hanyat = gy->jobb->oldalak[BOTTOM].db;
                akt.hova = hova;
                if (gy->oldalak[BOTTOM].db == 0)
                    gy->oldalak[BOTTOM].tol = hova;
                gy->oldalak[BOTTOM].db += akt.hanyat;
                hova += akt.hanyat;
                gy->optimized_add_mit_hova(akt);
            }
            if (gy->bal->oldalak[TOP].db != 0) {
                akt.is_be1 = true;
                akt.honnan = gy->bal->oldalak[TOP].tol;
                akt.hanyat = gy->bal->oldalak[TOP].db;
                akt.hova = hova;
                if (gy->oldalak[TOP].db == 0)
                    gy->oldalak[TOP].tol = hova;
                gy->oldalak[TOP].db += akt.hanyat;
                hova += akt.hanyat;
                gy->optimized_add_mit_hova(akt);
            }
        }
        if (gy->jobb->oldalak[TOP].db != 0) {
            akt.is_be1 = false;
            akt.honnan = gy->jobb->oldalak[TOP].tol;
            akt.hanyat = gy->jobb->oldalak[TOP].db;
            akt.hova = hova;
            if (gy->oldalak[TOP].db == 0)
                gy->oldalak[TOP].tol = hova;
            gy->oldalak[TOP].db += akt.hanyat;
            hova += akt.hanyat;
            gy->optimized_add_mit_hova(akt);
        }
    }
}


//***********************************************************************
void red_fa::feloszt(MezoTipus mt, const tomb3d<t_modell_cella>& r, const doboz & d, doboz & ki_1, doboz & ki_2, Oldal & kozos_oldal_1){
// Az 1. (azaz bal) cella közös oldala az EAST, NORT vagy TOP lehet. 
// (Ha nem így adódna, akkor ebben a fv-ben meg kell cserélni a két kimenetet.)
// Az mt a szimuláció típusát jelenti, de elth esetén a konkrét cella lehet csak el vagy th!
//***********************************************************************
    // az intelligens verziónál figyelj a d.is_lin-re is, ha true, akkor nem kell vizsgálni, hogy van-e nemlineáris
    ki_1 = d;
    ki_2 = d;
    if (d.x2 - d.x1 > d.y2 - d.y1) {
        if (d.x2 - d.x1 > d.z2 - d.z1) {
            ki_1.x2 = (d.x1 + d.x2) / 2;
            ki_2.x1 = (d.x1 + d.x2) / 2 + 1;
            kozos_oldal_1 = EAST;
        }
        else {
            ki_1.z2 = (d.z1 + d.z2) / 2;
            ki_2.z1 = (d.z1 + d.z2) / 2 + 1;
            kozos_oldal_1 = TOP;
        }
    }
    else if (d.y2 - d.y1 > d.z2 - d.z1) {
        ki_1.y2 = (d.y1 + d.y2) / 2;
        ki_2.y1 = (d.y1 + d.y2) / 2 + 1;
        kozos_oldal_1 = NORTH;
    }
    else {
        ki_1.z2 = (d.z1 + d.z2) / 2;
        ki_2.z1 = (d.z1 + d.z2) / 2 + 1;
        kozos_oldal_1 = TOP;
    }
}
