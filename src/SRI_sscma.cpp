#include "ErrMReals.h" 
#include "SRI_sscma.h"

using namespace extendedleaps;

extern "C" 
SEXP eleaps(SEXP S,SEXP S2,SEXP Si,SEXP Segval,SEXP Segvct,
	SEXP E,SEXP Ei,SEXP Hegvct,SEXP HegvctTinv,SEXP HegvctEinv,
	SEXP wilksval,SEXP bartpival,SEXP lawhotval,SEXP ccr12val,
	SEXP  r,SEXP kmin,SEXP kmax,SEXP nsol,
	SEXP exclude,SEXP include,SEXP nexclude,SEXP ninclude,
	SEXP criterion,SEXP fixed,SEXP pcindices,SEXP nbindices,
	SEXP dim,SEXP timelimit,SEXP maxaperr,SEXP onlyforward)
{
  SEXP subsets,values,bestsets,bestvalues,dimsub,dimval,dimbsets,ans;
  
  Rf_protect(r = AS_INTEGER(r));
  Rf_protect(kmin = AS_INTEGER(kmin));
  Rf_protect(kmax = AS_INTEGER(kmax));
  Rf_protect(nsol = AS_INTEGER(nsol));
  Rf_protect(nexclude = AS_INTEGER(nexclude));
  Rf_protect(ninclude = AS_INTEGER(ninclude));
  Rf_protect(fixed = AS_INTEGER(fixed));
  Rf_protect(nbindices = AS_INTEGER(nbindices));
  Rf_protect(dim = AS_INTEGER(dim));
  Rf_protect(onlyforward = AS_INTEGER(onlyforward));
  
	bool optimal,nomemory;
	int r1 = INTEGER(r)[0]; 
	int kmin1 = INTEGER(kmin)[0]; 
	int kmax1 = INTEGER(kmax)[0]; 
	int nsol1 = INTEGER(nsol)[0];
  int nexclude1 = INTEGER(nexclude)[0]; 
  int ninclude1 = INTEGER(ninclude)[0]; 
  int fixed1 = INTEGER(fixed)[0]; 
  int nbindices1 = INTEGER(nbindices)[0]; 
  int dim1 = INTEGER(dim)[0]; 
  int klength = kmax1 - kmin1 + 1;
	int checkcolinearity = INTEGER(onlyforward)[0];  
	
  Rf_protect(wilksval = AS_NUMERIC(wilksval));
	Rf_protect(bartpival = AS_NUMERIC(bartpival));
	Rf_protect(lawhotval = AS_NUMERIC(lawhotval));
	Rf_protect(ccr12val = AS_NUMERIC(ccr12val));
	Rf_protect(timelimit = AS_NUMERIC(timelimit));
	Rf_protect(maxaperr = AS_NUMERIC(maxaperr));
	
	double wilksval1 =  REAL(wilksval)[0];
	double bartpival1 = REAL(bartpival)[0];
	double lawhotval1 = REAL(lawhotval)[0];
	double ccr12val1 = REAL(ccr12val)[0];
	double timelimit1 = REAL(timelimit)[0];
	double maxaperr1 = REAL(maxaperr)[0];

	Rf_protect(criterion = AS_CHARACTER(criterion));
	const char* criterion1 = CHAR(STRING_ELT(criterion,0));
		
	if (!checkcolinearity) ErrMReals::errmonitreal<double>::dropec = true;   
	else ErrMReals::errmonitreal<double>::dropec = false;   

  Rf_protect(subsets = Rf_allocVector(INTSXP,nsol1*kmax1*klength));
  Rf_protect(values = Rf_allocVector(REALSXP,nsol1*klength));
  Rf_protect(bestsets = Rf_allocVector(INTSXP,kmax1*klength));
  Rf_protect(bestvalues = Rf_allocVector(REALSXP,klength));

	int retcode = extendedleaps::callsscma(
	  REAL(S),REAL(S2),REAL(Si),REAL(Segval),REAL(Segvct),
	  REAL(E),REAL(Ei),REAL(Hegvct),REAL(HegvctTinv),REAL(HegvctEinv),
	  wilksval1,bartpival1,lawhotval1,ccr12val1,
	  r1,kmin1,kmax1,nsol1,
	  INTEGER(exclude),INTEGER(include),nexclude1,ninclude1,
    criterion1,fixed1,INTEGER(pcindices),nbindices1,
	  dim1,timelimit1,maxaperr1,checkcolinearity,
	  INTEGER(subsets),REAL(values),REAL(bestvalues),INTEGER(bestsets),
	  false);

	if (retcode == 4) nomemory = true;
	else nomemory = false;
	if (retcode==0 || retcode==2)  optimal = true;
	else optimal = false;
	if (retcode==2 || retcode==3)  {
		Rprintf("\nWarning: Because of numerical problems caused by strong multicolinearity\n");
		Rprintf("some subsets were excluded from the analysis.\n");
		Rprintf("You can try to increase the number of subsets to be compared by reducing the value\n");
		Rprintf("of the function argument maxaperr but the numerical accuracy of results may be compromised\n\n");
	}

  Rf_protect(dimsub = Rf_allocVector(INTSXP,3));
	Rf_protect(dimval = Rf_allocVector(INTSXP,2));
	Rf_protect(dimbsets = Rf_allocVector(INTSXP,2));
	Rf_protect(ans = NEW_LIST(6));
	  
  INTEGER(dimsub)[0] = nsol1;
  INTEGER(dimsub)[1] = kmax1;
  INTEGER(dimsub)[2] = klength;
  SET_DIM(subsets,dimsub); 

  INTEGER(dimval)[0] = nsol1;
  INTEGER(dimval)[1] = klength;
  SET_DIM(values,dimval); 
  	
  INTEGER(dimbsets)[0] = klength;
  INTEGER(dimbsets)[1] = kmax1; 
  SET_DIM(bestsets,dimbsets); 

 	SET_VECTOR_ELT(ans, 0, subsets);
  SET_VECTOR_ELT(ans, 1, values);
  SET_VECTOR_ELT(ans, 2, bestvalues);
  SET_VECTOR_ELT(ans, 3, bestsets);
  SET_VECTOR_ELT(ans, 4, Rf_ScalarInteger(optimal));
  SET_VECTOR_ELT(ans, 5, Rf_ScalarInteger(nomemory));

	UNPROTECT(25);
	
  return(ans);
}


