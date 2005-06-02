#include <cassert>
#include <cmath>
#include "Sscma.h"
#include "Vsmabo.h"
#include "RVcrt.h"

using namespace leapsnbnds;

#ifdef COUNTING  
extern long unsigned fpcnt1;
#endif

rvgdata::rvgdata(vind nvariables)
  :   p(nvariables)
{
	tmpv.reserve(p);
	cndv.reserve(p);
	s2 = new symtwodarray(p);
	m1t.assign(p,vector<real>(p));
}

rvgdata::~rvgdata()
{
	delete s2;
}

rvdata::rvdata(vind lastvariab,vind partialnv,vind numbvar,rvgdata *data,const deque<bool>& active,vind *origvarlist,real criterion)
  :  lastv(lastvariab), k(partialnv), p(numbvar), gdt(data), varin(active), orgvar(origvarlist), crt(criterion), e(0)
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
		{ for (vind i=0;i<ivct.size();i++) delete ivct[i]; }
		throw;
	}
}

rvdata::~rvdata()
{
	for (vind i=0;i<ivct.size();i++) delete ivct[i]; 
	delete e;
}

real rvdata::updatecrt(vind *plist,vind *flist,vind var,vind fvarind) const
{
	assert(flist[var-1] >= fvarind && flist[var-1]  < p);

	vind irowi,varind; 
	if (plist) varind = plist[var-fvarind-1];
	else varind = var-fvarind-1;
	real newcrt,e1 = (*e)(varind,varind);
	real *cv=gdt->getcndv();
	static deque<bool> vin(p);

	vin = varin;
	if (varin[var-1]) vin[var-1] = false;
	else vin[var-1] = true;
	{  for (vind i=0;i<var-1;i++)	if (vin[i])  {
		irowi = flist[i];
		cv[i] = (*ivct[irowi])[varind]/e1;  
	} } 
	{  for (vind i=var;i<p;i++)  if (vin[i])  {
		irowi = flist[i];
		cv[i] = (*ivct[irowi])[varind]/e1;  
	} }
	if (vin[var-1]) cv[var-1] = 1./e1;
	cmpts2sm1(plist,flist,gdt->getm1t(),&vin[0],orgvar,var,var+1,var);
	newcrt = frobenius(gdt->getm1t(),&vin[0]);
	#ifdef COUNTING  
	fpcnt1 += p;
	#endif
	return newcrt;
}

void rvdata::pivot(vind vp,vind v1,vind vl,vind fvarind,real newcrt,vind* plist,vind* flist,subsetdata* newdtgpnt,bool last)
{
	vind inrowi,incoli,pivotind,fpivotind = flist[vp-1];
	if (plist) pivotind = plist[vp-fvarind-1];
	else pivotind = vp-fvarind-1;
	real pivotval = (*e)(pivotind,pivotind);
	real *cv = getgdata()->getcndv();

	rvdata *newdata = static_cast<rvdata *>(newdtgpnt);    // Attention: newdtgtpnt MUST point to rvdata object !!!

	// For safety, in debug mode use the alternative code with dynamic_cast and assert
	
//	rvdata *newdata = dynamic_cast<rvdata *>(newdtgpnt);    
//	assert(newdata);

	newdata->crt = newcrt;

	{ for (vind j=0;j<p;j++)  
		if (j+1 != vp) newdata->varin[j] = varin[j]; }
	if (varin[vp-1]) newdata->varin[vp-1] = false;
	else newdata->varin[vp-1] = true;
	{  for (vind j=0;j<v1-1;j++) 
		if (j+1 != vp && newdata->varin[j])  {
			inrowi = flist[j];
			cv[j] = (*ivct[inrowi])[pivotind]/pivotval;
			#ifdef COUNTING  
			fpcnt ++;
			#endif
		}
	}
	{ for (vind j=vl;j<p;j++) 
		if (newdata->varin[j])  {
			inrowi = flist[j];
			cv[j] = (*ivct[inrowi])[pivotind]/pivotval;
			#ifdef COUNTING  
			fpcnt ++;
			#endif
		}
	}
	
	if (!last)  {
		symatpivot(plist,pivotval,*e,*(newdata->e),vp-fvarind,v1-fvarind,vl-fvarind);
		for (vind j=0;j<v1-1;j++)  
			if (j+1 != vp && newdata->varin[j])  {
				inrowi = flist[j];
				vectorpivot(plist,*ivct[inrowi],*newdata->ivct[j],*e,cv[j],vp-fvarind,v1-fvarind,vl-fvarind); 
				newdata->ivct[j]->switchtoowndata();
			} 
		if (newdata->varin[vp-1]) {
			for (vind j=v1-1;j<vl;j++)  { 
				if (plist) incoli = plist[j-fvarind];
				else incoli = j-fvarind;
				newdata->ivct[vp-1]->setvalue(j+1-v1,-(*ivct[fpivotind])[incoli]/pivotval);  
			} 
			#ifdef COUNTING  
			fpcnt += (vl-v1+1);
			#endif
			newdata->ivct[vp-1]->switchtoowndata();
		}

		{ for (vind j=vl;j<p;j++)  
			if (newdata->varin[j])  {
				inrowi = flist[j];
				vectorpivot(plist,*ivct[inrowi],*newdata->ivct[j],*e,cv[j],vp-fvarind,v1-fvarind,vl-fvarind); 
				newdata->ivct[j]->switchtoowndata();
			} 
		}
	}

	if (newdata->varin[vp-1])  cv[vp-1] = 1./pivotval;
	{ for (vind j=v1-1;j<vl;j++)  
		if (newdata->varin[j]) {
			if (!plist) inrowi = j-fvarind;
			else inrowi = plist[j-fvarind];
			cv[j] = (*e)(inrowi,pivotind)/pivotval;
			#ifdef COUNTING  
			fpcnt += (vl-v1+1);
			#endif
		}
	}
	cmpts2sm1(plist,flist,newdata->s2m1,&(newdata->varin[0]),orgvar,vp,v1,vl);

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

void rvdata::cmpts2sm1(vind *plist,vind *flist,twodarray& outmat,bool *inlst,vind *orgvlst,vind vp,vind v1,vind vl) const
{
	real val,*tv=gdt->gettmpv(),*fl=gdt->getcndv();
	vind pivotind,fpivotind,inrowi,incoli,fvarind=lastv-k;
	if (plist) pivotind = plist[vp-fvarind-1];
	else pivotind = vp-fvarind-1;
	fpivotind = flist[vp-1];

	{ for (vind j=0;j<p;j++) {
		if ( (j+1 < v1 || j+1 > vl) && !inlst[j]) continue;
		incoli = flist[j];
		tv[j] = 0.;
		for (vind a=0;a<p;a++) 
			if (inlst[a] && a+1 != vp)  {
				inrowi = flist[a];
				tv[j] += -(*ivct[inrowi])[pivotind]*gdt->gets2(orgvlst[a],orgvlst[j]);  
				#ifdef COUNTING  
				fpcnt++;
				#endif
			}
	} }

	if (inlst[vp-1]) {
		{ for (vind i=0;i<p;i++)  
			if (inlst[i] && i+1 != vp) {
				inrowi = flist[i];
				for (vind j=0;j<vl;j++)  {
					if (j+1 < v1 && !inlst[j]) continue;
					incoli = flist[j];
					val = s2m1[inrowi][incoli] +  fl[i] * ( gdt->gets2(orgvlst[vp-1],orgvlst[j]) - tv[j] );
					outmat[i][j] = val;
					#ifdef COUNTING  
					fpcnt++;
					#endif
				}
		} }
		{ for (vind j=0;j<vl;j++)  {
			if (j+1 < v1 && !inlst[j]) continue;
			incoli = flist[j];
			val = fl[vp-1] * ( gdt->gets2(orgvlst[vp-1],orgvlst[j]) - tv[j] );
			outmat[vp-1][j] = val;
			#ifdef COUNTING  
			fpcnt++;
			#endif
		} }
  }

  else {
	for (vind i=0;i<p;i++)  
		if (inlst[i]) {
			inrowi = flist[i];
			for (vind j=0;j<p;j++)  {
				if ( (j+1 < v1 || j+1 > vl) && !inlst[j] ) continue;
				incoli = flist[j];
				val = s2m1[inrowi][incoli] + 
					  (*ivct[inrowi])[pivotind] * gdt->gets2(orgvlst[j],orgvlst[vp-1]) - fl[i]*tv[j];
				outmat[i][j] = val;
				#ifdef COUNTING  
				fpcnt += 2;
				#endif
			}
		}
  }
}
