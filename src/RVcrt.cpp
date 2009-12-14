#include <cassert>
#include <cmath>
#include "Sscma.h"
#include "Vsmabo.h"
#include "RVcrt.h"

namespace extendedleaps {

#ifdef COUNTING  
extern int fpcnt1;
#endif

rvgdata::rvgdata(vind nvariables)
  :   p(nvariables)
{
	s2 = new symtwodarray(p);
}

rvgdata::~rvgdata()
{
	delete s2;
}

partialrvdata::partialrvdata(vind nvariables)
  :   p(nvariables)
{
	tmpv.reserve(p);
	cndv.reserve(p);
	vin.resize(p);
	m1t.assign(p,vector<real>(p));
}

rvdata::rvdata(vind lastvariab,vind nvtopiv,vind tnv,rvgdata *data,const deque<bool>& active,vind *origvarlist,real criterion)
  :  lastv(lastvariab), p(tnv), k(nvtopiv), crt(criterion), varin(active), orgvar(origvarlist), e(0), gdt(data)
{
	try {
		if (k > 0)  {
			e = new symtwodarray(k);
			ivct.assign(p,0);
			for (vind i=0;i<p;i++) {
				if (i+k >= lastv) ivct[i] = new matvectarray(k,e,i-(lastv-k));
				else ivct[i] = new matvectarray(k,0,0);
			}
		}
		s2m1.assign(p,vector<real>(p));
	}
	catch (std::bad_alloc)   {
		delete e;
		{ for (unsigned i=0;i<ivct.size();i++) delete ivct[i]; }
		throw;
	}
}

rvdata::~rvdata()
{
	for (unsigned i=0;i<ivct.size();i++) delete ivct[i]; 
	delete e;
}

void  rvdata::getpdata(partialdata* pd)  
{ 
	partialrvdata *pdasrv = static_cast<partialrvdata *>(pd);    
	
	/* Attention: pd MUST point to partialrvdata object !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert    */

/*	partialrvdata *pdasrv = dynamic_cast<partialrvdata *>(pd);
	assert(pdasrv);                                               */

	setcriterion(pdasrv->getcrt());
	{ for (vind j=0;j<p;j++) varin[j] = pdasrv->vin[j]; }
	for (vind i=0;i<p;i++) if (varin[i])
		for (vind j=0;j<p;j++)  if (varin[j]) { 
			real tmp = pdasrv->m1t[i][j]; 
			s2m1[i][j] = tmp; 
		}
}

 real rvdata::updatecrt(direction dir,mindices& mmind,vind var,partialdata* pdt) const
{ 
	if (mmind.direct()) return updatecrt(dir,*(mmind.idpm()),*(mmind.idfm()),var,pdt); 
	else return updatecrt(dir,*(mmind.iipm()),*(mmind.iifm()),var,pdt); 
}
		
void rvdata::pivot(direction dir,mindices& mmind,vind vp,vind t,
						   partialdata* pdt,subsetdata* fdt,bool last)
{ 
	if (mmind.direct()) pivot(dir,*(mmind.idpm()),*(mmind.idfm()),vp,t,pdt,fdt,last); 
	else pivot(dir,*(mmind.iipm()),*(mmind.iifm()),vp,t,pdt,fdt,last); 
}

real rvdata::updatecrt(direction dir,lagindex<d>& prtmmit,itindex<d>& fmmind,vind var,partialdata* newdtpnt) const
{
	partialrvdata *newdata = static_cast<partialrvdata *>(newdtpnt);    
	
	/* Attention: newdtpnt MUST point to partialrvdata object !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert    */
	
/*	partialrvdata *newdata = dynamic_cast<partialrvdata *>(newdtpnt);
	assert(newdata);                                                         */
	
	vind varind = prtmmit[var-1];
	real newcrt,e1 = (*e)(varind,varind);
	real *cv = newdata->getcndv();
	deque<bool>& vin = newdata->vin;

	vin = varin;
	if (dir == forward) vin[var-1] = true;
	else vin[var-1] = false;
	fmmind.reset();
	for (vind i=0;i<p;fmmind++,i++)	
		if (vin[i] && (i!=var-1) )  cv[i] = (*ivct[fmmind()])[varind]/e1;
	if (dir == forward) cv[var-1] = 1./e1;
	cmpts2sm1(prtmmit,fmmind,newdata,newdata->getm1t(),orgvar,var,&vin[0],&vin[0]);
	newcrt = frobenius(newdata->getm1t(),&vin[0]);
	#ifdef COUNTING  
	fpcnt1 += p;
	#endif

	newdata->setpivotval(e1);
	newdata->setcrt(newcrt);
	return newcrt;
}


real rvdata::updatecrt(direction dir,lagindex<i>& prtmmit,itindex<i>& fmmind,vind var,partialdata* newdtpnt) const
{
	partialrvdata *newdata = static_cast<partialrvdata *>(newdtpnt);    
	
	/* Attention: newdtpnt MUST point to partialrvdata object !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert    */
	
/*	partialrvdata *newdata = dynamic_cast<partialrvdata *>(newdtpnt);
	assert(newdata);                                                         */
	
	vind varind = prtmmit[var-1];
	real newcrt,e1 = (*e)(varind,varind);
	real *cv = newdata->getcndv();
	deque<bool>& vin = newdata->vin;

	vin = varin;
	if (dir == forward) vin[var-1] = true;
	else vin[var-1] = false;
	fmmind.reset();
	for (vind i=0;i<p;fmmind++,i++)	
		if (vin[i] && (i!=var-1) )  cv[i] = (*ivct[fmmind()])[varind]/e1;
	if (dir == forward) cv[var-1] = 1./e1;
	cmpts2sm1(prtmmit,fmmind,newdata,newdata->getm1t(),orgvar,var,&vin[0],&vin[0]);
	newcrt = frobenius(newdata->getm1t(),&vin[0]);
	#ifdef COUNTING  
	fpcnt1 += p;
	#endif

	newdata->setpivotval(e1);
	newdata->setcrt(newcrt);
	return newcrt;
}

void rvdata::pivot(direction dir,lagindex<d>& prtmmit,itindex<d>& fmmind,vind vp,vind t,partialdata* newpdtpnt,subsetdata* newfdtpnt,bool last)
{
	vind pivotind,fpivotind = fmmind[vp-1];              
	pivotind = prtmmit[vp-1];

	partialrvdata* pdata = static_cast<partialrvdata *>(newpdtpnt);    
	rvdata* newdata = static_cast<rvdata *>(newfdtpnt);    
	
	/* Attention: pdtpnt and newdttpnt MUST point to partialrvdata and rvdata objects !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert       */
	
/*	partialrvdata* pdata = dynamic_cast<partialrvdata *>(newpdtpnt);
	rvdata* newdata = dynamic_cast<rvdata *>(newfdtpnt);
	assert(pdata && newdata);                                  */

	real pivotval = pdata->getpivotval();
	real *cv = pdata->getcndv();
	deque<bool>& colin = pdata->vin;

	symatpivot(prtmmit,pivotval,*e,*(newdata->e),vp,t);
	fmmind.reset();
	for (vind i=0;i<vp;fmmind++,i++)  
	if (newdata->varin[i])  {
		vectorpivot(prtmmit,*ivct[fmmind()],*newdata->ivct[i],*e,cv[i],vp,t); 
		newdata->ivct[i]->switchtoowndata();
	} 
	if (dir == forward)  {
		prtmmit.reset(vp);
		for (vind j=vp;j<vp+t;prtmmit++,j++)   
			newdata->ivct[vp-1]->setvalue(j-vp,-(*ivct[fpivotind])[prtmmit()]/pivotval);  
		#ifdef COUNTING  
		fpcnt += t;
		#endif
		newdata->ivct[vp-1]->switchtoowndata();
	}
	fmmind.reset(vp+t);
	{ for (vind i=vp+t;i<p;fmmind++,i++)  
		if (newdata->varin[i])  {
			vectorpivot(prtmmit,*ivct[fmmind()],*newdata->ivct[i],*e,cv[i],vp,t); 
			newdata->ivct[i]->switchtoowndata();
		} 
	}

	{ for (vind j=0;j<p;j++)
		if (j+1 > vp && j+1 <= vp+t && !colin[j]) colin[j] = true;
		else colin[j] = false;
	}
	cmpts2sm1(prtmmit,fmmind,pdata,newdata->s2m1,orgvar,vp,&(newdata->varin[0]),&colin[0]);
}


void rvdata::pivot(direction dir,lagindex<i>& prtmmit,itindex<i>& fmmind,vind vp,vind t,partialdata* newpdtpnt,subsetdata* newfdtpnt,bool last)
{
	vind pivotind,fpivotind = fmmind[vp-1];              
	pivotind = prtmmit[vp-1];

	partialrvdata* pdata = static_cast<partialrvdata *>(newpdtpnt);    
	rvdata* newdata = static_cast<rvdata *>(newfdtpnt);    
	
	/* Attention: pdtpnt and newdttpnt MUST point to partialrvdata and rvdata objects !!!
	   For safety, in debug mode use the alternative code with dynamic_cast and assert       */
	
/*	partialrvdata* pdata = dynamic_cast<partialrvdata *>(newpdtpnt);
	rvdata* newdata = dynamic_cast<rvdata *>(newfdtpnt);
	assert(pdata && newdata);                                  */

	real pivotval = pdata->getpivotval();
	real *cv = pdata->getcndv();
	deque<bool>& colin = pdata->vin;

	symatpivot(prtmmit,pivotval,*e,*(newdata->e),vp,t);
	fmmind.reset();
	for (vind i=0;i<vp;fmmind++,i++)  
	if (newdata->varin[i])  {
		vectorpivot(prtmmit,*ivct[fmmind()],*newdata->ivct[i],*e,cv[i],vp,t); 
		newdata->ivct[i]->switchtoowndata();
	} 
	if (dir == forward)  {
		prtmmit.reset(vp);
		for (vind j=vp;j<vp+t;prtmmit++,j++)   
			newdata->ivct[vp-1]->setvalue(j-vp,-(*ivct[fpivotind])[prtmmit()]/pivotval);  
		#ifdef COUNTING  
		fpcnt += t;
		#endif
		newdata->ivct[vp-1]->switchtoowndata();
	}
	fmmind.reset(vp+t);
	{ for (vind i=vp+t;i<p;fmmind++,i++)  
		if (newdata->varin[i])  {
			vectorpivot(prtmmit,*ivct[fmmind()],*newdata->ivct[i],*e,cv[i],vp,t); 
			newdata->ivct[i]->switchtoowndata();
		} 
	}

	{ for (vind j=0;j<p;j++)
		if (j+1 > vp && j+1 <= vp+t && !colin[j]) colin[j] = true;
		else colin[j] = false;
	}
	cmpts2sm1(prtmmit,fmmind,pdata,newdata->s2m1,orgvar,vp,&(newdata->varin[0]),&colin[0]);
}

void rvdata::cmpts2sm1(lagindex<d>&,itindex<d>&,partialrvdata* pdata,twodarray& outmat,vind* orgvlst,vind vp,bool* rowlst,bool* collst) const
{
	vind fvarind=lastv-k,pivotind=vp-fvarind-1;
	real *tv=pdata->gettmpv(),*fl=pdata->getcndv();

	for (vind j=0;j<p;j++) if (collst[j] ) {
		tv[j] = 0.;
		for (vind a=0;a<p;a++) { 
			if ( !rowlst[a] || a+1 == vp ) continue;
			tv[j] += -(*ivct[a])[pivotind]*gdt->gets2(orgvlst[a],orgvlst[j]);  
			#ifdef COUNTING  
			fpcnt++;
			#endif
		}
	} 

	if (rowlst[vp-1]) {
		for (vind i=0;i<p;i++)  {  
			if ( !rowlst[i] || i+1 == vp ) continue;
			for (vind j=0;j<p;j++)  if (collst[j] ) {
				outmat[i][j] = 
					s2m1[i][j] +  fl[i] * ( gdt->gets2(orgvlst[vp-1],orgvlst[j]) - tv[j] );
				#ifdef COUNTING  
				fpcnt++;
				#endif
			}
		} 
		for (vind j=0;j<p;j++)  if (collst[j] ) {
			outmat[vp-1][j] = 
				fl[vp-1] * ( gdt->gets2(orgvlst[vp-1],orgvlst[j]) - tv[j] );
			#ifdef COUNTING  
			fpcnt++;
			#endif
		} 
	}

	else {
		for (vind i=0;i<p;i++)  {  
			if (!rowlst[i] ) continue;
			for (vind j=0;j<p;j++)  if (collst[j] ) {
				outmat[i][j] =  s2m1[i][j] + 
					  (*ivct[i])[pivotind] * gdt->gets2(orgvlst[j],orgvlst[vp-1]) - fl[i]*tv[j];
				#ifdef COUNTING  
				fpcnt += 2;
				#endif
			}
		}
	}
}

void rvdata::cmpts2sm1(lagindex<i>& prtmmit,itindex<i>& fmmind,partialrvdata* pdata,twodarray& outmat,vind* orgvlst,vind vp,bool* rowlst,bool* collst) const
{
	real *tv=pdata->gettmpv(),*fl=pdata->getcndv();
	vind inrowi;
	vind pivotind=prtmmit[vp-1];                
	itindex<i>& rowind = fmmind;
	itindex<i> colind(fmmind);

	for (vind j=0;j<p;j++) if (collst[j] ) {
		tv[j] = 0.;
		rowind.reset();
		for (vind a=0;a<p;rowind++,a++) { 
			if ( !rowlst[a] || a+1 == vp ) continue;
			tv[j] += -(*ivct[rowind()])[pivotind]*gdt->gets2(orgvlst[a],orgvlst[j]);  
			#ifdef COUNTING  
			fpcnt++;
			#endif
		}
	} 

	if (rowlst[vp-1]) {
		rowind.reset();
		for (vind i=0;i<p;rowind++,i++)  {  
			if ( !rowlst[i] || i+1 == vp ) continue;
			inrowi = rowind();
			colind.reset();
			for (vind j=0;j<p;colind++,j++)  if (collst[j] ) {
				outmat[i][j] =
					s2m1[inrowi][colind()] +  fl[i] * ( gdt->gets2(orgvlst[vp-1],orgvlst[j]) - tv[j] );
				#ifdef COUNTING  
				fpcnt++;
				#endif
			}
		} 
		for (vind j=0;j<p;j++)  if (collst[j] ) {
			outmat[vp-1][j] = fl[vp-1] * ( gdt->gets2(orgvlst[vp-1],orgvlst[j]) - tv[j] );
			#ifdef COUNTING  
			fpcnt++;
			#endif
		} 
	}

	else {
		rowind.reset();
		for (vind i=0;i<p;rowind++,i++)  if (rowlst[i] ) {
			inrowi = rowind();
			colind.reset();
			for (vind j=0;j<p;colind++,j++)  if (collst[j] ) {
				outmat[i][j] = s2m1[inrowi][colind()] +
					  (*ivct[inrowi])[pivotind] * gdt->gets2(orgvlst[j],orgvlst[vp-1]) - fl[i]*tv[j];
				#ifdef COUNTING  
				fpcnt += 2;
				#endif
			}
		}
	}
}

real rvdata::frobenius(twodarray& m,bool *inlst) const
{
	real tmp = 0.;

	for (vind i=0;i<p;i++)  if (inlst[i]) {
		tmp += pow(m[i][i],2);   
		#ifdef COUNTING  
		fpcnt++;
		#endif
		for (vind j=0;j<i;j++)  if (inlst[j])  { 
			tmp += 2*m[i][j]*m[j][i];  
			#ifdef COUNTING  
			fpcnt++;
			#endif
		}
	}
	return tmp;
}

}
