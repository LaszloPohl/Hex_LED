//***********************************************************************
// hiba header
// Creation date:  2009. 07. 11.
// Creator:        Pohl László
//***********************************************************************


//***********************************************************************
#ifndef PL_SRFAJL_HEADER
#define	PL_SRFAJL_HEADER
//***********************************************************************
#include "PLString.h"
#include "listaestomb.h"
//***********************************************************************


//***********************************************************************
class srfajl{
//***********************************************************************
    tomb< tomb<PLString> > sorok;
    //***********************************************************************
    void szokozcsere(PLString & s){
    //***********************************************************************
        const unsigned n=s.Length();
        for(unsigned i=0;i<n;i++)
            if(plisspace(s[i])||s[i]==','||s[i]==';'||s[i]=='='||s[i]=='<'
                ||s[i]=='>' || s[i] == ':' || s[i] == '(' || s[i] == ')')s[i]=' ';
    }
    //***********************************************************************
    int egysegszamolo(const PLString & s){
    // nem kezdõdhet szóközzel!
    //***********************************************************************
        const unsigned n=s.Length();
        unsigned i,db,in=s[0]!='\"';
        for(i=db=1;i<n;i++){
            if(in && s[i-1]==' ' && s[i]!=' ')db++;
            if(s[i]=='\"')in=1-in;
        }
        return in?int(db):-1;
    }
    //***********************************************************************
    bool explode(const PLString & s,tomb<PLString> & t){
    // nem kezdõdhet szóközzel! Csak szóköz elválasztó van benne.
    //***********************************************************************
        const int db=egysegszamolo(s);
        const unsigned n=s.Length();
        if(db<1)return false;
        t.clear();
        t.resize(db);
        unsigned i=0,j,in=s[0]!='\"';
        enum all{kezd,normal,idezet};
        all a=kezd;
        for(j=0;j<n;j++){
            switch(a){
                case kezd: if(in){t[i]=s[j]; a=normal;}else a=idezet; break;
                case normal: 
                    if(s[j-1]==' ' && s[j]!=' ')i++;
                    if(s[j]!='\"'){if(s[j]!=' ')t[i]+=s[j];}
                    else a=idezet; 
                    break;
                case idezet: if(s[j]!='\"')t[i]+=s[j]; else a=normal; break;
            }
        }
        return true;
    }
    //***********************************************************************
    void rm_unsupported_version(PLString & s, const char *version){
    //***********************************************************************
        if (s.find(version) == 0)
            s.trunc();
    }
    //***********************************************************************
    void rm_supported_version_head(PLString & s, const char *version){
    //***********************************************************************
        if (s.find(version) == 0) {
            size_t db = strlen(version);
            for (size_t i = 0; i < db; i++)s[i] = ' ';
        }
    }
public:
    //***********************************************************************
    const tomb< tomb<PLString> > &lines(){return sorok;}
    //***********************************************************************
    void open(const PLString nev){
    //***********************************************************************
        unsigned sor=0;
        char s[1024];
        PLString ps;
        FILE *fp=fopen(nev.c_str(),"rt");
        if(fp==NULL)throw hiba("srfajl::open()","Cannot open %s to read",nev.c_str());
        sorok.clear();
        sorok.resize(8);
        while(fgets(s,1024,fp)){
            PLString ideiglenes=s;
            ideiglenes.rm_comment();
            rm_unsupported_version(ideiglenes, "V3");
            rm_unsupported_version(ideiglenes, "v3");
            rm_supported_version_head(ideiglenes, "V6");
            rm_supported_version_head(ideiglenes, "v6");
            if(ideiglenes.isspace())continue;
            ps += ideiglenes;
            //ps+=s;
            //ps.rm_comment();
            szokozcsere(ps);
            ps=trim(ps);
            if(ps[ps.Length()-1]!='\\'){
                if((sor+1)%8==0)sorok.resize(sor+9);
                if(!explode(ps,sorok[sor]))throw hiba("srfajl::open()","missing \" in line %u in %s",sor,nev.c_str());
                sor++;
                ps.trunc();
            }
            else ps[ps.Length()-1]=' ';
        }
        sorok.resize(sor);
		fclose(fp);
        return;
    }
};


#endif
