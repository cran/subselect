#ifndef VSMAC
#define VSMAC

#define NONE     1
#define GCD      2
#define RV       3 
#define MCB2     4 
#define DEFAULT  5 

extern const double EPS;
extern const unsigned VC;
extern const unsigned CORR;
extern const unsigned CHS;
extern const unsigned ALL;

extern unsigned p,p0,fp,lp,hdf,m,ms,df,q,daind,pcrt,hrcrt,aiccrt;
extern long unsigned cntg,cntp,fpcnt,fpcnt0,fpcnt1;
extern double   *lbnd,*ubnd,trs,trs2;
extern unsigned *prvks,*dmyv,*actv;  
extern char indRM[],indGCD[],indRV[];
extern char cmpcrMC2[],cmpcrGCD[],cmpcrGCD1[],cmpcrRV[];

extern unsigned n;
extern double *Fl,d0,w0,u0,v0,c0,t0,ic0,*vc0;
extern unsigned *Flp;
extern unsigned tmpcnt;

#endif
