//****************************************************************
// PLString osztály
// 2002-2009
// By Pohl László
// freeware
//****************************************************************


//***********************************************************************
#ifndef PLSTRING_SOURCE
#define PLSTRING_SOURCE
//***********************************************************************


//***********************************************************************
#include "PLString.h"
//***********************************************************************

#define MIN_STRING_SIZE 8

//****************************************************************
PLString::PLString(const char *s){
//****************************************************************
    t=new char[(len=unsigned(strlen(s)))+1+MIN_STRING_SIZE];
    strcpy(t,s);
}


//****************************************************************
PLString::PLString(const PLString& s){
//****************************************************************
    t=new char[(len=s.len)+1+MIN_STRING_SIZE];
    strcpy(t,s.t);
}


//****************************************************************
PLString PLString::operator +(const PLString& s)const{
//****************************************************************
    char *psum=new char[len+s.len+1];
    strcpy(psum,t);
    strcpy(psum+len,s.t);
    PLString sum(psum);
    delete [] psum;
    return sum;
}


//****************************************************************
PLString PLString::operator +(char c)const{
//****************************************************************
    char *psum=new char[len+2];
    strcpy(psum,t);
    psum[len]=c;
    psum[len+1]='\0';
    PLString sum(psum);
    delete [] psum;
    return sum;
}


//****************************************************************
PLString & PLString::operator +=(char c){
//****************************************************************
    if(len%MIN_STRING_SIZE==0){
        char *psum=new char[len+1+MIN_STRING_SIZE];
        strcpy(psum,t);
        delete [] t;
        t=psum;
   }
    t[len++]=c;
    t[len]=0;
    return *this;
}


//****************************************************************
PLString PLString::operator +(int a)const{
//****************************************************************
    char c[256];
    return *this+PLString(itoa(a,c,10));
}


//****************************************************************
PLString PLString::operator +(unsigned a)const{
//****************************************************************
    char c[256];
    return *this+PLString(itoa(a,c,10));
}


//****************************************************************
PLString PLString::operator +(long a)const{
//****************************************************************
    char c[256];
    return *this+PLString(ltoa(a,c,10));
}


//****************************************************************
PLString PLString::operator +(unsigned long a)const{
//****************************************************************
    char c[256];
    return *this+PLString(ltoa(a,c,10));
}


//****************************************************************
PLString PLString::operator +(double a)const{
//****************************************************************
    char c[256];
    sprintf(c,"%.10g",double(a));
    return *this+c;
}


//****************************************************************
PLString  operator +(int a, const PLString &s){
//****************************************************************
    char c[256];
    return PLString(itoa(a,c,10))+s;
}


//****************************************************************
PLString& PLString::operator =(const PLString& s){
//****************************************************************
    if(this!=&s){
        delete [] t;
        t=new char[(len=s.len)+1+MIN_STRING_SIZE];
        strcpy(t,s.t);
    }
    return *this;
}


//****************************************************************
PLString & PLString::operator =(const char c){
//****************************************************************
    delete [] t;
    t=new char[(len=1)+1+MIN_STRING_SIZE];
    t[0]=c;
    t[1]=0;
    return *this;
}


//****************************************************************
PLString PLString::substr(unsigned Start,unsigned n)const{
//****************************************************************
    if(Start>=len)return PLString("");
    if((n==npos)||(n>len-Start))n=len-Start;
    char *c=new char[n+1];
    for(unsigned i=Start;i<Start+n;i++)c[i-Start]=t[i];
    c[n]='\0';
    PLString s(c);
    delete [] c;
    return s;
}


//****************************************************************
int PLString::find(const char c,unsigned start)const{
//****************************************************************
    if(start==npos)return npos;
    for(unsigned i=start;i<len;i++)if(t[i]==c)return i;
    return npos;
}

           
//****************************************************************
int PLString::findr(const char c)const{
//****************************************************************
    for(int i=len-1;i>=0;i--)if(t[i]==c)return i;
    return npos;
}

           
//****************************************************************
int PLString::find(const char *c)const{
//****************************************************************
    int clen=int(strlen(c));
    for(int i=0;i<int(len)-clen;i++)if(t[i]==c[0]){
        bool b=true;
        for(int j=0;j<clen;j++)if(t[i+j]!=c[j])b=false;
        if(b)return i;
    }
    return npos;
}

 
//****************************************************************
int PLString::findr(const char *c)const{
//****************************************************************
    int clen=int(strlen(c));
    for(int i=len-clen-1;i>=0;i--)if(t[i]==c[0]){
        bool b=true;
        for(int j=0;j<clen;j++)if(t[i+j]!=c[j])b=false;
        if(b)return i;
    }
    return npos;
}

 
//****************************************************************
int PLString::find(const PLString &s)const{
//****************************************************************
    for(int i=0;i<=int(len)-int(s.len);i++)if(t[i]==s[0]){
        bool b=true;
        for(unsigned j=0;j<s.len;j++)if(t[i+j]!=s[j])b=false;
        if(b)return i;
    }
    return npos;
}


//****************************************************************
int PLString::findr(const PLString &s)const{
//****************************************************************
    for(int i=len-s.len;i>=0;i--)if(t[i]==s[0]){
        bool b=true;
        for(unsigned j=0;j<s.len;j++)if(t[i+j]!=s[j])b=false;
        if(b)return i;
    }
    return npos;
}


//****************************************************************
PLString getPath(const PLString& s){
// bemenet: útvonal+fájlnév
//****************************************************************
    PLString s2;
    int per=s.findr('/');
    int bs=s.findr('\\');
    if(per==npos&&bs==npos){s2="./";return s2;}
    if(per==npos){s2=s.substr(0,bs+1);return s2;}
    if(bs==npos){s2=s.substr(0,per+1);return s2;}
    s2=s.substr(0,per>bs?per+1:bs+1);
    return s2;
}


//****************************************************************
PLString getFileName(const PLString& s){
// bemenet: útvonal+fájlnév
//****************************************************************
    PLString s2;
    int per=s.findr('/');
    int bs=s.findr('\\');
    if(per==npos&&bs==npos)return s;
    if(per==npos){s2=s.substr(bs+1);return s2;}
    if(bs==npos){s2=s.substr(per+1);return s2;}
    s2=s.substr(per>bs?per+1:bs+1);
    return s2;
}


//****************************************************************
PLString getFileNameWithoutExtension(const PLString& s){
// bemenet: útvonal+fájlnév
//****************************************************************
    PLString s2=getFileName(s);
    int pont=s2.findr('.');
    if(pont==npos)return s;
    return s2.substr(0,pont);
}


//****************************************************************
PLString getExtension(const PLString& s){
// bemenet: útvonal+fájlnév
//****************************************************************
    PLString s2;
    int pont=s.findr('.');
    if(pont==npos)return s2;
    return s.substr(pont+1);
}


//****************************************************************
PLString trim(const PLString& s){
// bemenet: levágja a whitespace-eket
//****************************************************************
    int start=0,stop=s.Length()-1;
    while(plisspace(s.t[start])&&start<stop)start++;
    while(plisspace(s.t[stop])&&start<stop)stop--;
    if(start>=stop)return PLString();
    if(start==0&&stop==s.Length()-1)return s;
    return s.substr(start,stop-start+1);
}


#endif
