#include "SR_sscma.h"
#include "Sscma.h"

extern "C" 
int leaps(double *S,double *S2,double *Si,double *Segval,double *Segvct,
	int *kmin,int *kmax,int *nsol,
	int *exclude,int *include,int *nexclude,int *ninclude,
	char **criterion,int *fixed,int *pcindices,int *nbindices,
	 int *dim,double *timelimit,int *found,
	int *subsets,double *values,double *bvalues,int *bsets)
{
	return callsscma(S,S2,Si,Segval,Segvct,
		  	*kmin,*kmax,*nsol,
			exclude,include,*nexclude,*ninclude,
			*criterion,fixed,pcindices,*nbindices,*dim,*timelimit,found,
			subsets,values,bvalues,bsets);
}
