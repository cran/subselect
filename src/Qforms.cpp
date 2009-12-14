#include <cassert>
#include <cmath>
#include "Sscma.h"
#include "Vsmabo.h"
#include "Qforms.h"

namespace extendedleaps {

#ifdef COUNTING 
extern int fpcnt1;
#endif

partialqfdata::partialqfdata(vind nparcels)
  :   r(nparcels)
{
	tmpv.reserve(r);
}

qfdata::qfdata(vind tnv,vind nvtopiv,vind nparcels)
  :  p(tnv), k(nvtopiv), r(nparcels)
{
	ve.assign(r,vector<real>(k));
	e = new symtwodarray(k);
}

qfdata::~qfdata()
{
	delete e;
}

 void qfdata::pivot(direction dir,mindices& mmind,vind vp,vind t,partialqfdata* pdt,qfdata* fdt,bool last)
{ 
	if (mmind.direct()) pivot(*(mmind.idpm()),vp,t,pdt,fdt,last); 
	else pivot(*(mmind.iipm()),vp,t,pdt,fdt,last); 
}

void qfdata::pivot(lagindex<d>& prtmmit,vind vp,vind t,partialqfdata* newpdata,qfdata* newfdata,bool last)
{
	symatpivot(prtmmit,newpdata->getpivotval(),*e,*(newfdata->e),vp,t);
		for (vind j=0;j<r;j++) 
			vectorpivot(prtmmit,ve[j],newfdata->ve[j],*e,(newpdata->gettmpv())[j],vp,t); 
}

void qfdata::pivot(lagindex<i>& prtmmit,vind vp,vind t,partialqfdata* newpdata,qfdata* newfdata,bool last)
{
	symatpivot(prtmmit,newpdata->getpivotval(),*e,*(newfdata->e),vp,t);
		for (vind j=0;j<r;j++) 
			vectorpivot(prtmmit,ve[j],newfdata->ve[j],*e,(newpdata->gettmpv())[j],vp,t); 
}

real sqfdata::updatesum(mindices& mmind,vind var,partialsqfdata* pdt) const
{ 
	if (mmind.direct()) return updatesum((*(mmind.idpm()))[var-1],pdt); 
	else return updatesum((*(mmind.iipm()))[var-1],pdt); 
}
		
real sqfdata::updatesum(vind varind,partialsqfdata* newdata) const
{
	real *tv = newdata->gettmpv();
	real newsum = sum_,e1 = (*e)(varind,varind);

	for (vind i=0;i<r;i++) {
		tv[i] = ve[i][varind]/e1;
		newsum += tv[i]*ve[i][varind];
	}
	#ifdef COUNTING 
	fpcnt1 += 2*r;
	#endif

	newdata->setpivotval(e1);
	newdata->setsum(newsum);
	return newsum;
}

real singleqfdata::updatecrt(direction dir,mindices& mmind,vind var,partialdata* pdt) const
{
	partialsingleqfdata *newdata = static_cast<partialsingleqfdata *>(pdt);    
	
	/* Attention: pdt MUST point to partialsingleqfdata object !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert    */
	
/*	partialsingleqfdata *newdata = dynamic_cast<partialsingleqfdata *>(pdt);
	assert(newdata);                                                             */
	
	return qf->updatesum(mmind,var,newdata->pqf);  
}

void singleqfdata::pivot(direction dir,mindices& mmind,vind vp,vind t,partialdata* pdt,subsetdata* fdt,bool last)
{	
	partialsingleqfdata* newpdata = static_cast<partialsingleqfdata *>(pdt);    
	singleqfdata* newfdata = static_cast<singleqfdata *>(fdt);    
	
	/*  Attention: pdt and fdt MUST point to partialsingleqfdata and singleqfdata objects !!!
	    For safety, in debug mode use the alternative code with dynamic_cast and assert         */
	
/*	partialsingleqfdata* newpdata = dynamic_cast<partialsingleqfdata *>(pdt);
	singleqfdata* newfdata = dynamic_cast<singleqfdata *>(fdt);
	assert(newpdata && newfdata);                                */

	qf->pivot(dir,mmind,vp,t,newpdata->pqf,newfdata->qf,last);  
}

}
