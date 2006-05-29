extern "C" 
int leaps(double *S,double *S2,double *Si,double *Segval,double *Segvct,
	double *E,double *Ei,double *Hegvct,double *HegvctTinv,double *HegvctEinv,
	double *wilksval,double *bartpival,double *lawhotval,double *ccr12val,int* r,
	int *kmin,int *kmax,int *nsol,
	int *exclude,int *include,int *nexclude,int *ninclude,
	char **criterion,int *fixed,int *pcindices,int *nbindices,
	 int *dim,double *timelimit,int *found,
	int *subsets,double *values,double *bvalues,int *bsets);	

namespace extendedleaps {
		
int callsscma(double *S,double *S2,double *Si,double *Segval,double *Segvct,
	double *E,double *Ei,double *Hegvct,double *HegvctTinv,double *HegvctEinv,
	double wilksval,double bartpival,double lawhotval,double ccr12val,int r,
	int kmin,int kmax,int nsol,int *out,int *in,int nout,int nin,
	char *cmpcr,int fixed,int *pcind,int nind,int nvar,double timelimit,
	int *found,int *subs,double *subsv,double *bestsv,int *bests);
}

