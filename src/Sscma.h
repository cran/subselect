#ifndef SSCMA
#define SSCMA

#include <string>
#include "ErrMReals.h" 

using namespace ErrMReals; 

namespace extendedleaps {


	/*  Extended Leaps special types  */

typedef  short int  vind;    /* Integer type used to index variables. Should be able to run
                                from 0 to the largest possible number of variables under comparison  */

// typedef  double     real;    // Floating point type used to represent real numbers.
typedef  ErrMReals::errmonitreal<double>	real;   // Floating point type used to represent real numbers	


enum pcskept {given,firstk}; 		 /*  Enumeration used to identify the set of PCs considered by the GCD criterion  */
enum direction {forward,backward};	 /*  Enumeration used to specifiy the specification of the criterion updates */

	/*   Extended Leaps global declarations  */

/*  Constants shared by more than two files  */

const int MINIMZ   = 0;	/* Comparision criterion is to be minimized        */
const int MAXIMZ   = 1;	/* Comparision criterion is to be maximized        */


/*  Variables shared by more than two files  */
 
extern vind  p,q,fp,lp,mindim,lastvar,flst,*actv;  
extern long unsigned ms;                                                     
extern double  *lbnd,*ubnd;
extern double btime,maxtime,numtol;
extern pcskept pcsets;
extern std::string memmsg;

#ifdef COUNTING 
extern int cntg,fpcnt,fpcnt1;	/*  Floating point operation counters  */ 
#endif                

/*  Functions shared by more than two files   */

bool Leaps_Search(vind frwind0,vind bckind0,vind fvind,vind lvind,vind nvfrwd,vind nvbckwrd);	 
bool Rev_Leaps_Search(vind frwind0,vind bckind0,vind fvind,vind lvind,vind nvfrwd,vind nvbckwrd);
bool Forward_BreadthF_Search(vind frwind0,vind nvfrwd,vind fvind,vind lvind);
void crtwrksp(void);
void savfull(void);
void savfrst(void);
void msg(const std::string&);
void errmsg(const std::string&);

}

#endif
