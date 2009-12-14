#include <cassert>
#include <cmath>
#include "Sscma.h"
#include "Vsmabo.h"
#include "RMcrt.h"

namespace extendedleaps {

#ifdef COUNTING  
extern int fpcnt1;
#endif

partialrmdata::partialrmdata(vind nvariables)
  :   p(nvariables)
{
	tmpv.reserve(p);
}

rmdata::rmdata(vind lastvariab,vind nvtopiv,vind tnv,rmgdata *data,const deque<bool>& active,real criterion)
  :  lastv(lastvariab), p(tnv), k(nvtopiv), crt(criterion), varin(active), e(0), gdt(data)
{
	if (k > 0) try {
		e = new symtwodarray(k);
		ovct.assign(p,0);
		{ for (vind i=0;i<p;i++) {
			if (i+k >= lastv) ovct[i] = new matvectarray(k,e,i-(lastv-k));
			else ovct[i] = new matvectarray(k,0,0);
		} }
	}
	catch (std::bad_alloc)   {
		delete e;
		{ for (unsigned i=0;i<ovct.size();i++) delete ovct[i]; }
		throw;
	}
}

rmdata::~rmdata()
{
	{ for (unsigned i=0;i<ovct.size();i++) delete ovct[i]; }
	delete e;
}

real rmdata::updatecrt(direction dir,mindices& mmind,vind var,partialdata* pdt) const
{ 
	if (mmind.direct()) return updatecrt(dir,(*(mmind.idfm())),var,(*(mmind.idpm()))[var-1],pdt); 
	else return updatecrt(dir,*(mmind.iifm()),var,(*(mmind.iipm()))[var-1],pdt); 
}
		
 void rmdata::pivot(direction dir,mindices& mmind,vind vp,vind t,
						   partialdata* pdt,subsetdata* fdt,bool last)
{ 
	if (mmind.direct()) pivot(dir,*(mmind.idpm()),*(mmind.idfm()),vp,t,pdt,fdt,last); 
	else pivot(dir,*(mmind.iipm()),*(mmind.iifm()),vp,t,pdt,fdt,last); 
}

real rmdata::updatecrt(direction dir,itindex<d>& fmmind,vind var,vind varind,partialdata* newdtpnt) const
{
	partialrmdata *newdata = static_cast<partialrmdata *>(newdtpnt);    
	
	/*  Attention: newdtpnt MUST point to partialrmdata object !!!
	    For safety, in debug mode use the alternative code with dynamic_cast and assert    */
	
/*	partialrmdata *newdata = dynamic_cast<partialrmdata *>(newdtpnt);
	assert(newdata);                                                       */
	
	real *tv = newdata->gettmpv();
	real newcrt=crt,e1 = (*e)(varind,varind);

	if (dir == forward) newcrt -= e1;
	else newcrt -= 1./e1;
	fmmind.reset();
	for (vind i=0;i<p;fmmind++,i++) {
		if (!varin[i] && (i!=var-1) ) {
			tv[i] = (*ovct[fmmind()])[varind]/e1;
			newcrt -= tv[i] * (*ovct[fmmind()])[varind];
			#ifdef COUNTING  
			fpcnt1 += 2;
			#endif
	}
	}
	newdata->setpivotval(e1);
	newdata->setcrt(newcrt);
	return newcrt;
}

real rmdata::updatecrt(direction dir,itindex<i>& fmmind,vind var,vind varind,partialdata* newdtpnt) const
{
	partialrmdata *newdata = static_cast<partialrmdata *>(newdtpnt);    
	
	/*  Attention: newdtpnt MUST point to partialrmdata object !!!
	    For safety, in debug mode use the alternative code with dynamic_cast and assert    */
	
/*	partialrmdata *newdata = dynamic_cast<partialrmdata *>(newdtpnt);
	assert(newdata);                                                       */
	
	real *tv = newdata->gettmpv();
	real newcrt=crt,e1 = (*e)(varind,varind);

	if (dir == forward) newcrt -= e1;
	else newcrt -= 1./e1;
	fmmind.reset();
	for (vind i=0;i<p;fmmind++,i++) {
		if (!varin[i] && (i!=var-1) ) {
			tv[i] = (*ovct[fmmind()])[varind]/e1;
			newcrt -= tv[i] * (*ovct[fmmind()])[varind];
			#ifdef COUNTING  
			fpcnt1 += 2;
			#endif
	}
	}
	newdata->setpivotval(e1);
	newdata->setcrt(newcrt);
	return newcrt;
}

void rmdata::pivot(direction dir,lagindex<d>& prtmmit,itindex<d>& fmmind,vind vp,vind t,partialdata* newpdtpnt,subsetdata* newfdtpnt,bool last)
{
	partialrmdata* newpdata = static_cast<partialrmdata *>(newpdtpnt);    
	rmdata* newfdata = static_cast<rmdata *>(newfdtpnt);    
	
	/*  Attention: newpdtpnt and newfdttpnt MUST point to partialrmdata and rmdata objects !!!
	    For safety, in debug mode use the alternative code with dynamic_cast and assert           */
	
/*	partialrmdata* newpdata = dynamic_cast<partialrmdata *>(newpdtpnt);
	rmdata* newfdata = dynamic_cast<rmdata *>(newfdtpnt);
	assert(newpdata && newfdata);                                              */

	real pivotval = newpdata->getpivotval();
	real *tv = newpdata->gettmpv();

	{ for (vind i=0;i<p;i++)  
		if (i+1 != vp) newfdata->varin[i] = varin[i]; }
	if (dir == backward) newfdata->varin[vp-1] = false;
	else newfdata->varin[vp-1] = true;
	symatpivot(prtmmit,pivotval,*e,*(newfdata->e),vp,t);
	fmmind.reset();
	{ for (vind i=0;i<vp;fmmind++,i++)  
		if (i+1 != vp && !newfdata->varin[i])  {
			vectorpivot(prtmmit,*ovct[fmmind()],*newfdata->ovct[i],*e,tv[i],vp,t); 
			newfdata->ovct[i]->switchtoowndata();
	} }
	if (dir == backward) {
		prtmmit.reset(vp);
		for (vind j=vp;j<vp+t;prtmmit++,j++) 
			newfdata->ovct[vp-1]->setvalue(j-vp,-(*ovct[fmmind[vp-1]])[prtmmit()]/pivotval); 
		#ifdef COUNTING  
		fpcnt += t;
		#endif
		newfdata->ovct[vp-1]->switchtoowndata();
	}
	fmmind.reset(vp+t);
	{ for (vind i=vp+t;i<p;fmmind++,i++)  
		if (!newfdata->varin[i])  {
			vectorpivot(prtmmit,*ovct[fmmind()],*newfdata->ovct[i],*e,tv[i],vp,t); 
			newfdata->ovct[i]->switchtoowndata();
	} }
}

void rmdata::pivot(direction dir,lagindex<i>& prtmmit,itindex<i>& fmmind,vind vp,vind t,partialdata* newpdtpnt,subsetdata* newfdtpnt,bool last)
{
	partialrmdata* newpdata = static_cast<partialrmdata *>(newpdtpnt);    
	rmdata* newfdata = static_cast<rmdata *>(newfdtpnt);    
	
	/*  Attention: newpdtpnt and newfdttpnt MUST point to partialrmdata and rmdata objects !!!
	    For safety, in debug mode use the alternative code with dynamic_cast and assert           */
	
/*	partialrmdata* newpdata = dynamic_cast<partialrmdata *>(newpdtpnt);
	rmdata* newfdata = dynamic_cast<rmdata *>(newfdtpnt);
	assert(newpdata && newfdata);                                              */

	real pivotval = newpdata->getpivotval();
	real *tv = newpdata->gettmpv();

	{ for (vind i=0;i<p;i++)  
		if (i+1 != vp) newfdata->varin[i] = varin[i]; }
	if (dir == backward) newfdata->varin[vp-1] = false;
	else newfdata->varin[vp-1] = true;
	symatpivot(prtmmit,pivotval,*e,*(newfdata->e),vp,t);
	fmmind.reset();
	{ for (vind i=0;i<vp;fmmind++,i++)  
		if (i+1 != vp && !newfdata->varin[i])  {
			vectorpivot(prtmmit,*ovct[fmmind()],*newfdata->ovct[i],*e,tv[i],vp,t); 
			newfdata->ovct[i]->switchtoowndata();
	} }
	if (dir == backward) {
		prtmmit.reset(vp);
		for (vind j=vp;j<vp+t;prtmmit++,j++) 
			newfdata->ovct[vp-1]->setvalue(j-vp,-(*ovct[fmmind[vp-1]])[prtmmit()]/pivotval); 
		#ifdef COUNTING  
		fpcnt += t;
		#endif
		newfdata->ovct[vp-1]->switchtoowndata();
	}
	fmmind.reset(vp+t);
	{ for (vind i=vp+t;i<p;fmmind++,i++)  
		if (!newfdata->varin[i])  {
			vectorpivot(prtmmit,*ovct[fmmind()],*newfdata->ovct[i],*e,tv[i],vp,t); 
			newfdata->ovct[i]->switchtoowndata();
	} }
}

}
