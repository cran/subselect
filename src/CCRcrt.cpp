#include <cmath>
#include "Sscma.h"
#include "Vsmabo.h"
#include "CCRcrt.h"

namespace extendedleaps {

partialccrdata::partialccrdata(vind nvars,vind hrank)
{
	nvar = nvars;
	bptmpv.reserve(hrank);
}

partialccrdata::partialccrdata(vind nvars,vind hrank,real r2,real w,real bp)
  :	  ccr12(r2), wilksst(w), bartpist(bp)
{
	nvar = nvars;
	bptmpv.reserve(hrank);
}

ccrdata::ccrdata(vind nv,vind tnv,vind nvtopiv,vind hr,real w,real bp,real r2)
  :  p(tnv), k(nvtopiv), hrank(hr), nvar(nv), ccr12(r2), wilksst(w), bartpist(bp) 
{
	htinv.assign(hrank,vector<real>(k));
	emat = new symtwodarray(k);
	tmat = new symtwodarray(k);
}

ccrdata::~ccrdata(void)
{ 
	delete emat; 
	delete tmat; 
}

void  ccrdata::getpdata(partialdata* pd)  
{ 
	partialccrdata *pdasccr = static_cast<partialccrdata *>(pd);    
	
	/* Attention: pd MUST point to partialccrdata object !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert  */
	
/*	partialccrdata *pdasccr = dynamic_cast<partialccr *>(pd);	*/
/*	assert(pdasccr);                                                */

	ccr12 = pdasccr->getcrt();
	wilksst = pdasccr->getwilks();
	bartpist = pdasccr->getbartpi();
	nvar = pdasccr->nvar;

}

real ccrdata::updatecrt(direction dir,mindices& mmind,vind var,partialdata* pdt) const
{ 
	if (mmind.direct()) return updatecrt(dir,(*(mmind.idpm()))[var-1],pdt); 
	else return updatecrt(dir,(*(mmind.iipm()))[var-1],pdt); 
}

 void ccrdata::pivot(direction dir,mindices& mmind,vind vp,vind t,partialdata* pdt,subsetdata* fdt,bool last)
{ 
	if (mmind.direct()) pivot(*(mmind.idpm()),vp,t,pdt,fdt,last); 
	else pivot(*(mmind.iipm()),vp,t,pdt,fdt,last); 
}

void ccrdata::updatest(real& newwilksst,real& newbartpist,vind varind,partialccrdata* newdata) const  
{  
	real e1 = (*emat)(varind,varind);
	real t1 = (*tmat)(varind,varind);
	real *tv = newdata->getbptmpv();

	newwilksst = wilksst * (e1/t1);
	newbartpist = bartpist;
	for (vind i=0;i<hrank;i++) {
		tv[i] = htinv[i][varind]/t1;
		newbartpist += tv[i]*htinv[i][varind];
	}

	#ifdef COUNTING 
	fpcnt1 += 4;
	#endif

	newdata->setepivot(e1);
	newdata->settpivot(t1);
	newdata->setwilks(newwilksst);
	newdata->setbartpi(newbartpist);

	return;
} 

void ccrdata::pivot(lagindex<d>& prtmmit,vind vp,vind t,partialdata* newpdtpnt,subsetdata* newfdtpnt,bool last)
{	
	partialccrdata* newpdata = static_cast<partialccrdata *>(newpdtpnt);    
	ccrdata* newfdata = static_cast<ccrdata *>(newfdtpnt);    
	
	/* Attention: newpdtpnt and newfdtpnt MUST point to partialccrdata and ccrdata objects !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert     */
	
/*	partialccrdata* newpdata = dynamic_cast<partialccrdata *>(newpdtpnt);
	ccrdata* newfdata = dynamic_cast<ccrdata *>(newfdtpnt);                                
	assert(newpdata && newfdata);                                                          */

	symatpivot(prtmmit,newpdata->getepivot(),*emat,*(newfdata->emat),vp,t);
	symatpivot(prtmmit,newpdata->gettpivot(),*tmat,*(newfdata->tmat),vp,t);
	for (vind j=0;j<hrank;j++) 
		vectorpivot(prtmmit,htinv[j],newfdata->htinv[j],*tmat,(newpdata->getbptmpv())[j],vp,t); 
} 

void ccrdata::pivot(lagindex<i>& prtmmit,vind vp,vind t,partialdata* newpdtpnt,subsetdata* newfdtpnt,bool last)
{	
	partialccrdata* newpdata = static_cast<partialccrdata *>(newpdtpnt);    
	ccrdata* newfdata = static_cast<ccrdata *>(newfdtpnt);    
	
	/* Attention: newpdtpnt and newfdtpnt MUST point to partialccrdata and ccrdata objects !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert     */
	
/*	partialccrdata* newpdata = dynamic_cast<partialccrdata *>(newpdtpnt);
	ccrdata* newfdata = dynamic_cast<ccrdata *>(newfdtpnt);                                
	assert(newpdata && newfdata);                                                          */

	symatpivot(prtmmit,newpdata->getepivot(),*emat,*(newfdata->emat),vp,t);
	symatpivot(prtmmit,newpdata->gettpivot(),*tmat,*(newfdata->tmat),vp,t);
	for (vind j=0;j<hrank;j++) 
		vectorpivot(prtmmit,htinv[j],newfdata->htinv[j],*tmat,(newpdata->getbptmpv())[j],vp,t); 
} 

real rnk2ccrdata::updatecrt(direction dir,vind varind,partialdata* newdtpnt) const  
{  
	
	partialccrdata *newdata = static_cast<partialccrdata *>(newdtpnt);    
	
	/* Attention: newdtpnt MUST point to partialccrdata object !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert   */
	
/*	partialccrdata *newdata = dynamic_cast<partialccrdata *>(pdt);
	assert(newdata);                                                   */

	real newwilksst,newbartpist,newccr12;

	updatest(newwilksst,newbartpist,varind,newdata);
	if (dir==forward) newdata->nvar=nvar+1 ; 
	else newdata->nvar=nvar-1;
	
	if (newdata->nvar == 1) newccr12 = newbartpist; 
 	else  {
		newccr12 = 0.5 * ( newbartpist +  sqrt(newbartpist*newbartpist -4.*(newbartpist+newwilksst-1.)) );

		#ifdef COUNTING 
		fpcnt1 += 3;
		#endif
	}
	

	newdata->setcrt(newccr12);
	return newccr12;
} 

}
