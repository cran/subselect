#ifndef SSCMA
#define SSCMA

#include <string>

namespace leapsnbnds  {


					//  SSCMA Special types

typedef  short int  vind;    	// Integer type used to index variables. Should be able to run 
                             	// from 0 to the largest possible number of variables under comparison. 
typedef  double     real;    	// Floating point type used to represent real numbers.                        

enum   pcskept {given,firstk};  //  Enumeration used to identify the set of PCs considered by the GCD criterion

			          	
					//  SSCMA Global declarations


//  Constants shared by several files

const unsigned SRC	= 1;
const unsigned INV      = 0;

const unsigned MINIMZ   = 0;
const unsigned MAXIMZ   = 1;

const int  GCD      = 1;
const int  RV       = 2;
const int  MCB2     = 3;

const int NOTFOUND = 99;


//  Variables shared by several files


extern vind  p,q,fp,lp,mindim,lastvar,flst,*actv;  
extern int pcrt;    
extern long unsigned ms,sbsetind,maxsbst,*sbsetcnt;                                                  
extern real  *lbnd,*ubnd,trs,trs2;
extern double btime,maxtime; 
extern pcskept pcsets;
extern std::string memmsg;

#ifdef COUNTING 
extern long unsigned cntg,fpcnt,fpcnt1;		//  Floating point operations counters   
#endif                

//  Functions shared by several files

void crtwrksp(void);
void savfull(void);
void savfrst(void);
bool prcksp(vind,vind,vind,vind,vind);
void msg(const std::string&);
void errmsg(const std::string&);

}

#endif
