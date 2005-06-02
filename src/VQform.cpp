#include <cassert>
#include <cmath>
#include <vector>
#include "Sscma.h"
#include "Vsmabo.h"
#include "Qform.h"
#include "VQform.h"

using namespace leapsnbnds;

#ifdef COUNTING  
extern long unsigned fpcnt1;
#endif

vqfdata::vqfdata(vind totalnv,vind partialnv,vind nparcels,vqfauxmem* mem,vqfgdata* data,vind sbdimension,real criterion)
  :  qfdata(totalnv,partialnv,nparcels,0,criterion), auxmem(mem), gdt(data), dim(sbdimension)
{
	real vc0 = 0.;
	if (gdt->getdirection() == backward) vc0 = 1.;
	vc.assign(r,vc0);
}

vqfdata::~vqfdata()  {  }

real vqfdata::updatecrt(vind *,vind *flist,vind var,vind fvarind) const
{
	assert(flist[var-1] >= fvarind && flist[var-1]  < p);

	vind maxk=0,varind=flist[var-1]-fvarind;
	real inc,newcrt,e1=(*e)(varind,varind);
	real *tv=getauxmem()->gettmpv(),*newvc=getauxmem()->gettmpvc();
	dir d = gdt->getdirection();
	switch (d)  {
		case forward:
			maxk = dim+1;
			if (maxk > r) maxk = r;
			newcrt = crt + vc[dim];
			break;
		case backward:
			maxk = dim-1;
			if (maxk > r) maxk = r;
			if (r > dim-1) newcrt = crt - vc[dim-1];
			else newcrt = crt;
			break;
	}
	for (vind j=0;j<maxk;j++) {
		newvc[j] = vc[j];
		tv[j] = ve[j][varind]/e1;
		newvc[j] += (inc = tv[j]*ve[j][varind]);
		newcrt += inc;
	} 
	#ifdef COUNTING  
	fpcnt1 += 2*maxk;
	#endif
	return newcrt;
}

void vqfdata::pivot(vind vp,vind v1,vind vl,vind fvarind,real newcrt,vind *list,vind *,subsetdata * newdtgpnt,bool last)
{
	vind pivotind,maxk=0;
	if (list) pivotind = list[vp-fvarind-1];
	else pivotind = vp-fvarind-1;
	real pivotval = (*e)(pivotind,pivotind);
	real *tv=getauxmem()->gettmpv(),*newvc=getauxmem()->gettmpvc();
	dir d = getgdata()->getdirection();

	vqfdata *newdata = static_cast<vqfdata *>(newdtgpnt);    //Attention: newdtgtpnt MUST point to vqfdata object !!!

	// For safety, in debug mode use the alternative code with dynamic_cast and assert
	
//	vqfdata *newdata = dynamic_cast<vqfdata *>(newdtgpnt);    
//	assert(newdata);

	newdata->crt = newcrt;
	switch (d)  {
		case forward:
			maxk = (newdata->dim = dim+1) + vl-vp;
			if (maxk > r) maxk = r;
			break;
		case backward:
			maxk = newdata->dim = dim-1;
			if (maxk > r) maxk = r;
			break;
	}
	{ for (vind j=0;j<maxk;j++) {
		if (j < newdata->dim) newdata->vc[j] = newvc[j];
		else {
			tv[j] = ve[j][pivotind]/pivotval;
			newdata->vc[j] = vc[j] + tv[j]*ve[j][pivotind];
			#ifdef COUNTING  
			fpcnt += 2*maxk;
			#endif
		}
	} }
	if (!last)  {
		symatpivot(list,pivotval,*e,(*newdata->e),vp-fvarind,v1-fvarind,vl-fvarind);
		{ for (vind j=0;j<maxk;j++) 
			vectorpivot(list,ve[j],newdata->ve[j],*e,tv[j],vp-fvarind,v1-fvarind,vl-fvarind); }
	}
}
