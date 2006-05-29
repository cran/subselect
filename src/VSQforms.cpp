#include <cassert>
#include <cmath>
#include <vector>
#include "Sscma.h"
#include "Vsmabo.h"
#include "Qforms.h"
#include "VSQforms.h"

namespace extendedleaps {

#ifdef COUNTING  
extern long unsigned fpcnt1;
#endif

vsqfdata::vsqfdata(vind tnv,vind nvtopiv,vind nparcels,real vc0,real sum)
  :  sqfdata(tnv,nvtopiv,nparcels,sum)
{
	vc.reserve(nparcels);
	vc.assign(nparcels,vc0);
}

vsqfdata::vsqfdata(vind tnv,vind nvtopiv,vind nparcels,const vector<real>& ovc,real sum)
  :  sqfdata(tnv,nvtopiv,nparcels,sum)
{
	vc.reserve(nparcels);
	vc.assign(ovc.begin(),ovc.end());
}

vsqfdata::~vsqfdata()  {  }

void  vsqfdata::setvc(double* x)  
{ 
	for (vind j=0;j<r;j++) vc[j] = x[j];
}

inline real vsqfdata::updatesum(direction d,mindices& mmind,vind var,vind dim,partialvsqfdata* pdt) const
{ 
	if (mmind.direct()) return updatesum(d,(*(mmind.idpm()))[var-1],dim,pdt); 
	else return updatesum(d,(*(mmind.iipm()))[var-1],dim,pdt); 
}

inline void vsqfdata::pivot(direction d,mindices& mmind,vind vp,vind t,vind dim,partialvsqfdata* pdt,vsqfdata* fdt,bool last)
{ 
	if (mmind.direct()) pivot(d,*(mmind.idpm()),vp,t,dim,pdt,fdt,last); 
	else pivot(d,*(mmind.iipm()),vp,t,dim,pdt,fdt,last); 
}

real vsqfdata::updatesum(direction d,vind varind,vind dim,partialvsqfdata* newdata) const
{
	vind maxk=0;
	real inc,newsum,e1=(*e)(varind,varind);
	real *tv = newdata->gettmpv(),*newvc=newdata->gettmpvc();
	switch (d)  {
		case forward:
			maxk = dim+1;
			if (maxk > r) maxk = r;
			newsum = qfsum() + vc[dim];
			break;
		case backward:
			maxk = dim-1;
			if (maxk > r) maxk = r;
			if (r > dim-1) newsum = qfsum() - vc[dim-1];
			else newsum = qfsum();
			break;
	}
	for (vind j=0;j<maxk;j++) {
		newvc[j] = vc[j];
		tv[j] = ve[j][varind]/e1;
		newvc[j] += (inc = tv[j]*ve[j][varind]);
		newsum += inc;
	} 
	newdata->setpivotval(e1);
	newdata->setsum(newsum);

	#ifdef COUNTING  
	fpcnt1 += 2*maxk;
	#endif
	return newsum;
}

template<accesstp tp> 
void vsqfdata::pivot(direction d,lagindex<tp>& prtmmit,vind vp,vind t,vind dim,partialvsqfdata* newpdata,vsqfdata* newfdata,bool last)
{
	vind pivotind,newdim,maxk=0;
	pivotind = prtmmit[vp-1];
	real pivotval = newpdata->getpivotval();
	real *tv = newpdata->gettmpv();

	switch (d)  {
		case forward:
			maxk = (newdim = dim+1) + t;
			if (maxk > r) maxk = r;
			break;
		case backward:
			maxk = newdim = dim-1;
			if (maxk > r) maxk = r;
			break;
	}
	{ for (vind j=newdim;j<maxk;j++) {
		tv[j] = ve[j][pivotind]/pivotval;
		newfdata->vc[j] = vc[j] + tv[j]*ve[j][pivotind];
	} }
	#ifdef COUNTING  
	if (maxk > newdim) fpcnt += 2*(maxk - newfdata->dim);
	#endif

	symatpivot(prtmmit,pivotval,*e,*(newfdata->e),vp,t);
	for (vind j=0;j<maxk;j++) 
		vectorpivot(prtmmit,ve[j],newfdata->ve[j],*e,tv[j],vp,t); 
}

}
