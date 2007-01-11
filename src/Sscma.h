#ifndef SSCMA
#define SSCMA

#include <string>

namespace extendedleaps {


	/*  Extended Leaps special types  */

typedef  short int  vind;    /* Integer type used to index variables. Should be able to run
                                from 0 to the largest possible number of variables under comparison  */
typedef  double     real;    /* Floating point type used to represent real numbers.                  */

enum   pcskept {given,firstk};  /*  Enumeration used to identify the set of PCs considered by the GCD criterion  */

	/*   Extended Leaps global declarations  */

/*  Constants shared by more than two files  */

const int MINIMZ   = 0;	/* Comparision criterion is to be minimized        */
const int MAXIMZ   = 1;	/* Comparision criterion is to be maximized        */

/*  Variables shared by more than two files  */

extern vind  p,q,fp,lp,mindim,lastvar,flst,*actv;  
extern int ms;                                                     
extern real  *lbnd,*ubnd;   
extern double btime,maxtime; 
extern pcskept pcsets;
extern std::string memmsg;

#ifdef COUNTING 
extern int cntg,fpcnt,fpcnt1;	/*  Floating point operation counters  */ 
#endif                

/*  Functions shared by more than two files   */

void crtwrksp(void);
void savfull(void);
void savfrst(void);
bool prcksp(vind,vind,vind,vind,vind);
void msg(const std::string&);
void errmsg(const std::string&);

}

#endif
