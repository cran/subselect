#include "SRI_sscma.h"

using namespace extendedleaps;

extern "C" 
int leaps(double *S,double *S2,double *Si,double *Segval,double *Segvct,
	double *E,double *Ei,double *Hegvct,double *HegvctTinv,double *HegvctEinv,
	double *wilksval,double *bartpival,double *lawhotval,double *ccr12val,int* r,
	int *kmin,int *kmax,int *nsol,
	int *exclude,int *include,int *nexclude,int *ninclude,
	char **criterion,int *fixed,int *pcindices,int *nbindices,
	 int *dim,double *timelimit,int *found,
	int *subsets,double *values,double *bvalues,int *bsets)
{
	return callsscma(S,S2,Si,Segval,Segvct,E,Ei,Hegvct,HegvctTinv,HegvctEinv,
			*wilksval,*bartpival,*lawhotval,*ccr12val,*r,
		  	*kmin,*kmax,*nsol,
			exclude,include,*nexclude,*ninclude,
			*criterion,*fixed,pcindices,*nbindices,*dim,*timelimit,found,
			subsets,values,bvalues,bsets);
}
