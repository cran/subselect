#include <cassert>
#include <cmath>
#include "Sscma.h"
#include "Vsmabo.h"
#include "RMcrt.h"

using namespace leapsnbnds;

#ifdef COUNTING  
extern long unsigned fpcnt1;
#endif

rmgdata::rmgdata(vind nvariables)
  :   p(nvariables)
{
	tmpv.reserve(p);
}

rmdata::rmdata(vind lastvariab,vind partialnv,vind numbvar,rmgdata *data,const deque<bool>& active,real criterion)
  :  lastv(lastvariab), k(partialnv), p(numbvar), gdt(data), crt(criterion), e(0), varin(active)
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
		{ for (vind i=0;i<ovct.size();i++) delete ovct[i]; }
		throw;
	}
}

rmdata::~rmdata()
{
	for (vind i=0;i<ovct.size();i++) delete ovct[i]; 
	delete e;
}

real rmdata::updatecrt(vind *,vind *flist,vind var,vind fvarind) const
{

	assert(flist[var-1] >= fvarind && flist[var-1]  < p);

	vind irowi,varind = flist[var-1]-fvarind;
	real newcrt=crt,e1 = (*e)(varind,varind);

	if (!varin[var-1]) newcrt -= e1;
	else newcrt -= 1./e1;
	{ for (vind i=0;i<var-1;i++) if (!varin[i]) {
		irowi = flist[i] ;
		newcrt -= pow((*ovct[irowi])[varind],2)/e1;
		#ifdef COUNTING  
		fpcnt1 += 2;
		#endif
	}}
	{ for (vind i=var;i<p;i++) if (!varin[i]) {
		irowi = flist[i] ;
		newcrt -= pow((*ovct[irowi])[varind],2)/e1;
		#ifdef COUNTING  
		fpcnt1 += 2;
		#endif
	}}

	return newcrt;
}

void rmdata::pivot(vind vp,vind v1,vind vl,vind fvarind,real newcrt,vind* plist,vind* flist,subsetdata * newdtgpnt,bool last)
{
	vind irowi,icoli,pivotind;
	if (plist) pivotind = plist[vp-fvarind-1];
	else pivotind = vp-fvarind-1;
	real tmp,pivotval = (*e)(pivotind,pivotind);
	real *tv = getgdata()->gettmpv();

	rmdata *newdata = static_cast<rmdata *>(newdtgpnt);    //Attention: newdtgtpnt MUST point to rmdata object !!!

	// For safety, in debug mode use the alternative code with dynamic_cast and assert
	
//	rmdata *newdata = dynamic_cast<rmdata *>(newdtgpnt);    
//	assert(newdata);

	newdata->crt = newcrt;

	{ for (vind j=0;j<p;j++)  
		if (j+1 != vp) newdata->varin[j] = varin[j]; }
	if (varin[vp-1]) newdata->varin[vp-1] = false;
	else newdata->varin[vp-1] = true;
	{  for (vind j=0;j<v1-1;j++) 
		if (j+1 != vp && !newdata->varin[j])  {
			irowi = flist[j];
			tmp = (*ovct[irowi])[pivotind];
			tv[j] = tmp/pivotval;
			#ifdef COUNTING  
			fpcnt ++;
			#endif
		}
	}
	{ for (vind j=vl;j<p;j++) 
		if (!newdata->varin[j])  {
			irowi = flist[j];
			tmp = (*ovct[irowi])[pivotind];
			tv[j] = tmp/pivotval;
			#ifdef COUNTING  
			fpcnt ++;
			#endif
		}
	}

	if (!last)  {
		symatpivot(plist,pivotval,*e,*(newdata->e),vp-fvarind,v1-fvarind,vl-fvarind);
		{ for (vind j=0;j<v1-1;j++)  
			if (j+1 != vp && !newdata->varin[j])  {
				irowi = flist[j];
				vectorpivot(plist,*ovct[irowi],*newdata->ovct[j],*e,tv[j],vp-fvarind,v1-fvarind,vl-fvarind); 
				newdata->ovct[j]->switchtoowndata();
		} }
		if (!newdata->varin[vp-1]) {
			{ for (vind j=v1-1;j<vl;j++) {
				irowi = flist[vp-1];
				if (plist) icoli = plist[j-fvarind];
				else icoli = j-fvarind;
				newdata->ovct[vp-1]->setvalue(j+1-v1,-(*ovct[irowi])[icoli]/pivotval); 
			} }
			#ifdef COUNTING  
			fpcnt += (vl-v1+1);
			#endif
			newdata->ovct[vp-1]->switchtoowndata();
		}
		{ for (vind j=vl;j<p;j++)  
			if (!newdata->varin[j])  {
				irowi = flist[j];
				vectorpivot(plist,*ovct[irowi],*newdata->ovct[j],*e,tv[j],vp-fvarind,v1-fvarind,vl-fvarind); 
				newdata->ovct[j]->switchtoowndata();
		} }
	}
}
