#include <cmath>
#include "Sscma.h"
#include "Vsmabo.h"
#include "SpecialArrays.h"
#include "Qforms.h"
#include "MStcrt.h"

namespace extendedleaps {

wilksdata::wilksdata(vind nv,vind tnv,vind nvtopiv,vind hr,real wst)
  :  nvar(nv), p(tnv), k(nvtopiv), hrank(hr), wilksst(wst)
{
	emat = new symtwodarray(k);
	tmat = new symtwodarray(k);
}

wilksdata::~wilksdata(void)
{ 
	delete emat; 
	delete tmat; 
}

const real wilksdata::indice(void)	const	
{
	if (hrank < nvar) return 1. - pow(wilksst,1./hrank); 
	else return 1. - pow(wilksst,1./nvar); 
} 

void  wilksdata::getpdata(partialdata* pd)  
{ 
	partialwilksdata *pdaswilks = static_cast<partialwilksdata *>(pd);    
	
	/* Attention: pd MUST point to partialwilksdata object !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert    */
	
/*	partialwilksdata *pdaswilks = dynamic_cast<partialwilksdata *>(pd);
	assert(pdaswilks);                                                      */

	wilksst = pdaswilks->getcrt();
	nvar = pdaswilks->nvar;
}

real wilksdata::updatecrt(direction d,mindices& mmind,vind var,partialdata* pdt) const
{ 
	if (mmind.direct()) return updatecrt(d,(*(mmind.idpm()))[var-1],pdt); 
	else return updatecrt(d,(*(mmind.iipm()))[var-1],pdt); 
}
   
		
void wilksdata::pivot(direction,mindices& mmind,vind vp,vind t,partialdata* pdt,subsetdata* fdt,bool last)
{ 
	if (mmind.direct()) pivot(*(mmind.idpm()),vp,t,pdt,fdt,last); 
	else pivot(*(mmind.iipm()),vp,t,pdt,fdt,last); 
}

real wilksdata::updatecrt(direction d,vind varind,partialdata* newdtpnt) const  
{  
	partialwilksdata *newdata = static_cast<partialwilksdata *>(newdtpnt);    
	
	/* Attention: newdtpnt MUST point to partialwilksdata object !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert    */
	
/*	partialwilksdata *newdata = dynamic_cast<partialwilksdata *>(pdt);
	assert(newdata);                                                       */

	if (d==forward) newdata->nvar=nvar+1 ; 
	else newdata->nvar=nvar-1; 
	real e1 = (*emat)(varind,varind);
	real t1 = (*tmat)(varind,varind);
	real newwilksst = wilksst * (e1/t1);

	#ifdef COUNTING 
	fpcnt1 += 2;
	#endif

	newdata->setepivot(e1);
	newdata->settpivot(t1);
	newdata->setcrt(newwilksst);

	return newwilksst;
} 

void wilksdata::pivot(lagindex<d>& prtmmit,vind vp,vind t,partialdata* newpdtpnt,subsetdata* newfdtpnt,bool last)
{	
	partialwilksdata* newpdata = static_cast<partialwilksdata *>(newpdtpnt);    
	wilksdata* newfdata = static_cast<wilksdata *>(newfdtpnt);    
	
	/*  Attention: newpdtpnt and newfdtpnt MUST point to partialwilksdata and wilksdata objects !!!
	    For safety, in debug mode use the alternative code with dynamic_cast and assert             */
	
/*	partialwilksdata* newpdata = dynamic_cast<partialwilksdata *>(newpdtpnt);
	wilksdata* newfdata = dynamic_cast<wilksdata *>(newfdtpnt);
	assert(newpdata && newfdata);                                              */

	symatpivot(prtmmit,newpdata->getepivot(),*emat,*(newfdata->emat),vp,t);
	symatpivot(prtmmit,newpdata->gettpivot(),*tmat,*(newfdata->tmat),vp,t);
} 

void wilksdata::pivot(lagindex<i>& prtmmit,vind vp,vind t,partialdata* newpdtpnt,subsetdata* newfdtpnt,bool last)
{	
	partialwilksdata* newpdata = static_cast<partialwilksdata *>(newpdtpnt);    
	wilksdata* newfdata = static_cast<wilksdata *>(newfdtpnt);    
	
	/*  Attention: newpdtpnt and newfdtpnt MUST point to partialwilksdata and wilksdata objects !!!
	    For safety, in debug mode use the alternative code with dynamic_cast and assert             */
	
/*	partialwilksdata* newpdata = dynamic_cast<partialwilksdata *>(newpdtpnt);
	wilksdata* newfdata = dynamic_cast<wilksdata *>(newfdtpnt);
	assert(newpdata && newfdata);                                              */

	symatpivot(prtmmit,newpdata->getepivot(),*emat,*(newfdata->emat),vp,t);
	symatpivot(prtmmit,newpdata->gettpivot(),*tmat,*(newfdata->tmat),vp,t);
} 

partialtracedata::partialtracedata(vind nvars,vind hrank)		
{
	nvar = nvars;
	pqf = new partialsqfdata(hrank); 
}

  const real partialtracedata::getcrt(void) const	
{ 
	return pqf->getsum(); 
}

tracedata::tracedata(vind nv,vind tnv,vind nvtopiv,vind hr,real crt)
	: hrank(hr)
{
	nvar = nv;
	sqf = new sqfdata(tnv,nvtopiv,hr,crt);		
}

		
const real tracedata::criterion(void) const	
{ 
	return sqf->qfsum();  
}

void tracedata::setcriterion(real c)			
{
	sqf->setqfsum(c); 
}

void  tracedata::getpdata(partialdata* pd)  
{ 
	partialtracedata *pdastracest = static_cast<partialtracedata *>(pd);    
	
	/* Attention: pd MUST point to partialtracedata object !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert     */
	
/*	partialtracedata *pdasfgcd = dynamic_cast<partialtracedata *>(pd);
	assert(pdasfgcd);                                                    */

	setcriterion(pdastracest->getcrt());
	nvar = pdastracest->nvar;
}

real tracedata::updatecrt(direction d,mindices& mmind,vind var,partialdata* pdt) const
{  
	partialtracedata *newdata = static_cast<partialtracedata *>(pdt);    
	
	/* Attention: newdtpnt MUST point to partialtracedata object !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert     */
	
/*	partialtracedata *newdata = dynamic_cast<partialtracedata *>(pdt);
	assert(newdata);                                                    */

	if (d==forward) newdata->nvar=nvar+1 ; 
	else newdata->nvar=nvar-1; 
	return sqf->updatesum(mmind,var,newdata->pqf);  
} 

void tracedata::pivot(direction d,mindices& mmind,vind vp,vind t,partialdata* pdt,subsetdata* fdt,bool last)
{	
	partialtracedata* newpdata = static_cast<partialtracedata *>(pdt);    
	tracedata* newfdata = static_cast<tracedata *>(fdt);    
	
	/* Attention: pdt and fdt MUST point to partialtracedata and gcddata objects !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert    */
	
/*	partialtracedata* newpdata = dynamic_cast<partialtracedata *>(pdt);
	tracedata* newfdata = dynamic_cast<tracedata *>(fdt);
	assert(newpdata && newfdata);                                  */

	sqf->pivot(d,mmind,vp,t,newpdata->pqf,newfdata->sqf,last);  
} 

const real bartpistdata::indice(void)	const	
{
	if (hrank < nvar) return criterion()/hrank; 
	else return criterion()/nvar; 
} 

const real lawlhotstdata::indice(void)	const	
{
	if (hrank < nvar) return criterion()/(criterion()+hrank); 
	else return criterion()/(criterion()+nvar); 
} 

}
