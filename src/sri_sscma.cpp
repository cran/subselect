// SRI_sscma.cpp : Defines the entry point for the DLL application.
//

#include "SRI_sscma.h"
#include "SR_sscma.h"

extern "C"  
int leaps(double *mat,int *kmin,int *kmax,int *nsol,
		  int *exclude,int *include,int *nexclude,int *ninclude,
		  char **criterion,int *pcindices,int *nbindices,
		  int *dim,double *timelimit,int *found,
		  int *subsets,double *values,double *bvalues,int *bsets)
{
	return callsscma(mat,*kmin,*kmax,*nsol,exclude,include,*nexclude,*ninclude,
		         *criterion,pcindices,*nbindices,*dim,*timelimit,
			 found,subsets,values,bvalues,bsets);
}

