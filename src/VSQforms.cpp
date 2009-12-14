#include <cassert>
#include <cmath>
#include <vector>
#include "Sscma.h"
#include "Vsmabo.h"
#include "Qforms.h"
#include "VSQforms.h"

namespace extendedleaps {

#ifdef COUNTING  
extern int fpcnt1;
#endif

partialvsqfdata::partialvsqfdata(vind nparcels,real vc0)
  :  partialsqfdata(nparcels)
{
	tmpvc.reserve(nparcels);
	tmpvc.assign(nparcels,vc0);
}

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

void  vsqfdata::setvc(real* x,vind nparcels)  
{ 
	for (vind j=0;j<nparcels;j++) vc[j] = x[j];
}

real vsqfdata::updatesum(direction dir,mindices& mmind,vind var,vind dim,partialvsqfdata* pdt) const
{ 
	if (mmind.direct()) return updatesum(dir,(*(mmind.idpm()))[var-1],dim,pdt); 
	else return updatesum(dir,(*(mmind.iipm()))[var-1],dim,pdt); 
}

void vsqfdata::pivot(direction dir,mindices& mmind,vind vp,vind t,vind dim,partialvsqfdata* pdt,vsqfdata* fdt,bool last)
{ 
	if (mmind.direct()) pivot(dir,*(mmind.idpm()),vp,t,dim,pdt,fdt,last); 
	else pivot(dir,*(mmind.iipm()),vp,t,dim,pdt,fdt,last); 
}

real vsqfdata::updatesum(direction dir,vind varind,vind dim,partialvsqfdata* newdata) const
{
	vind maxk=0;
	real inc,newsum=0.,e1=(*e)(varind,varind);
	real *tv = newdata->gettmpv(),*newvc=newdata->gettmpvc();
	switch (dir)  {
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

void vsqfdata::pivot(direction dir,lagindex<d>& prtmmit,vind vp,vind t,vind dim,partialvsqfdata* newpdata,vsqfdata* newfdata,bool last)
{
	vind pivotind,newdim=0,maxk=0;
	pivotind = prtmmit[vp-1];
	real pivotval = newpdata->getpivotval();
	real *tv = newpdata->gettmpv();

	switch (dir)  {
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

void vsqfdata::pivot(direction dir,lagindex<i>& prtmmit,vind vp,vind t,vind dim,partialvsqfdata* newpdata,vsqfdata* newfdata,bool last)
{
	vind pivotind,newdim=0,maxk=0;
	pivotind = prtmmit[vp-1];
	real pivotval = newpdata->getpivotval();
	real *tv = newpdata->gettmpv();

	switch (dir)  {
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
