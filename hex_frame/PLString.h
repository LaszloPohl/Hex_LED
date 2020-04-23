//****************************************************************
// PLString osztály
// 2002-2009
// By Pohl László
// freeware
//****************************************************************

#if !defined PLSTRINGCLASS
#define PLSTRINGCLASS

//****************************************************************
#pragma warning(disable : 4996)
#include <string.h>
#include "isspace.h"
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
//****************************************************************

const int npos=-1;
const double floatvalue=123456789e123; // Ha a PLString todouble tagja FLOAT stringet kap, ezt az értéket helyettesíti be

//****************************************************************
class PLString{
//****************************************************************
    char     *t;    // char *
    unsigned len;   // int
public:
    PLString(const char *s="");
    PLString(const PLString& s);
    ~PLString(){delete [] t;}

    char&    operator[](size_t i)const{return t[i];}
    PLString operator+(const PLString& s)const;
    PLString operator+(const char *s)const{return *this+PLString(s);}
    friend PLString operator+(const char *s,const PLString s2){return PLString(s)+s2;}
    PLString operator+(char c)const;
    PLString&operator+=(char c);
    PLString&operator+=(const PLString& s){return *this=*this+s;}
    PLString operator+(int a)const;
    PLString operator+(unsigned a)const;
    PLString operator+(long a)const;
    PLString operator+(unsigned long a)const;
    PLString operator+(double a)const;
    PLString&operator=(const char c);
    PLString&operator=(const PLString& s);
    PLString&operator=(char *s){*this=PLString(s);return *this;}
    PLString&operator=(const char *s){*this=PLString(s);return *this;}
    bool     operator==(const PLString& s)const{return strcmp(t,s.t)==0;}
    bool     operator!=(const PLString& s)const{return !(*this==s);}
    bool     operator==(const char *s)const{return strcmp(t,s)==0;}
    bool     operator!=(const char *s)const{return !(*this==s);}
    bool     operator<(const PLString& s)const{return strcmp(t,s.t)<0;}
    bool     operator>(const PLString& s)const{return strcmp(t,s.t)>0;}
        
    PLString substr(unsigned Start=0,unsigned n=npos)const;
    int      Length()const{return len;}
    PLString UpCase()const{//A string nagybetûs változatát    adja vissza, az eredeti nem változik
        PLString s=*this;
        for(unsigned i=0;i<len;i++)s.t[i]=(char)toupper(s.t[i]);
        return s;
    }
    PLString LowCase()const{//A string nagybetûs változatát adja vissza, az eredeti nem változik
        PLString s=*this;
        for(unsigned i=0;i<len;i++)s.t[i]=(char)tolower(s.t[i]);
        return s;
    }
    PLString trim()const{
        PLString s=*this;
        for(unsigned i=s.Length()-1;i>0&&plisspace(s[i]);s[i--]=0);
        for(unsigned i=0;s[i]!=0;i++)if(!plisspace(s[i]))return PLString(s.t+i);
        return PLString("");
    }
    void    trunc(unsigned n=0){if(n<len){len=n;t[n]=0;}}
    void    rm_comment(){int n; if((n=find("//"))!=npos)trunc(n);}

    int     find(const char c,unsigned start=0)const;
    int     find(const char *c)const;
    int     find(const PLString &s)const;
    int     findr(const char c)const;
    int     findr(const char *c)const;
    int     findr(const PLString &s)const;
    
    bool todouble(double & d,unsigned start=0)const{
        if(npos==start)return false;
        if(sscanf(t+start,"%lg",&d)!=1){
            if(UpCase()=="FLOAT"){
                d=floatvalue;
                return true;
            }
            return false;
        }
        else return true;
    }
    void replace(char mit,char mire){
        for(unsigned i=0;i<len;i++)if(t[i]==mit)t[i]=mire;
    }
//    bool    todouble(double & d)const{d=atof(t);if((d==0.0)&&(t[0]!='0')){if((*this).UpCase()=="FLOAT"){d=floatvalue;return true;}return false;}return true;}
    bool    toint(int & i,int start=0)const{if(npos==start)return false;i=atoi(t+start);if((i==0)&&(t[0+start]!='0'))return false;return true;}
    bool    tounsigned(unsigned & u,int start=0)const{if(npos==start)return false;u=unsigned(atoi(t+start));if((u==0)&&(t[0+start]!='0'))return false;return true;}
    bool    tolong(long & l)const{l=atol(t);if((l==0)&&(t[0]!='0'))return false;return true;}

    bool    to2double(double & d1, double & d2,bool egyislehet=false)const{int x=find('|'); if(egyislehet&&x==npos){ bool b=todouble(d1); d2=d1; return b;}else return (x==npos||!todouble(d1)||!todouble(d2,x+1))?false:true;}
    bool    to2uns(unsigned & u1, unsigned & u2)const{int x=find('|'); return (x==npos||!tounsigned(u1)||!tounsigned(u2,x+1))?false:true;}
    bool    to3uns(unsigned & u1, unsigned & u2, unsigned & u3)const{int x=find('|'),y=find('|',x+1); return (x==npos||y==npos||!tounsigned(u1)||!tounsigned(u2,x+1)||!tounsigned(u3,y+1))?false:true;}
    bool    to2double1uns(double & d1, double & d2, unsigned & u)const{int x=find('|'),y=find('|',x+1); return (x==npos||y==npos||!todouble(d1)||!todouble(d2,x+1)||!tounsigned(u,y+1))?false:true;}
    // a to1uns2double megengedi, hogy elhagyjuk bármelyik értéket, akkor nem változtatja meg a megfelelõ paramétert
    bool    to1uns2double(unsigned & u,double & d1, double & d2)const{int x=find('|'),y=find('|',x+1); return ((len==0||tounsigned(u))&&(x==npos||todouble(d1,x+1))&&(y==npos||todouble(d2,y+1))?true:false);}
    char  * c_str()const{return t;}
    bool    isspace(){for(unsigned i=0; i<len; i++)if(!::isspace(t[i]))return false; return true;}

    friend  PLString operator +(int a,PLString& s);
    friend  PLString getPath(const PLString& s);
    friend  PLString getFileName(const PLString& s);
    friend  PLString getFileNameWithoutExtension(const PLString& s);
    friend  PLString getExtension(const PLString& s);
    friend  PLString trim(const PLString & s);
};
#endif
