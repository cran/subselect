#include <cassert>
#include <cmath>
#include "Sscma.h"
#include "Vsmabo.h"
#include "Qform.h"

using namespace leapsnbnds;

#ifdef COUNTING 
extern long unsigned fpcnt1;
#endif

qfauxmem::qfauxmem(vind nparcels)
  :   r(nparcels)
{
	tmpv.reserve(r);
}

qfauxmem::~qfauxmem()  {  }

qfdata::qfdata(vind totalnv,vind partialnv,vind nparcels,qfauxmem* mem,real criterion)
  :  p(totalnv), k(partialnv), r(nparcels), auxmem(mem), crt(criterion)
{
	ve.assign(r,vector<real>(k));
	e = new symtwodarray(k);
}

qfdata::~qfdata()
{
	delete e;
}

real qfdata::updatecrt(vind *,vind *flist,vind var,vind fvarind) const
{
	assert(flist[var-1] >= fvarind && flist[var-1]  < p);

	real *tv = getauxmem()->gettmpv();
	vind varind = flist[var-1]-fvarind;
	real newcrt = crt,e1 = (*e)(varind,varind);

	for (vind i=0;i<r;i++) {
		tv[i] = ve[i][varind]/e1;
		newcrt += tv[i]*ve[i][varind];
	}
	#ifdef COUNTING 
	fpcnt1 += 2*r;
	#endif
	return newcrt;
}

void qfdata::pivot(vind vp,vind v1,vind vl,vind fvarind,real newcrt,vind *list,vind *,subsetdata * newdtgpnt,bool last)
{
	vind pivotind;
	if (list) pivotind = list[vp-fvarind-1];
	else pivotind = vp-fvarind-1;
	real pivotval = (*e)(pivotind,pivotind);
	real *tv = getauxmem()->gettmpv();

	qfdata *newdata = static_cast<qfdata *>(newdtgpnt);    //Attention: newdtgtpnt MUST point to qfdata object !!!

	// For safety, in debug mode use the alternative code with dynamic_cast and assert
	
//	qfdata *newdata = dynamic_cast<qfdata *>(newdtgpnt);    
//	assert(newdata);

	newdata->crt = newcrt;
	if (!last)  {
		symatpivot(list,pivotval,*e,*(newdata->e),vp-fvarind,v1-fvarind,vl-fvarind);
		for (vind j=0;j<r;j++) 
			vectorpivot(list,ve[j],newdata->ve[j],*e,tv[j],vp-fvarind,v1-fvarind,vl-fvarind); 
	}
}
