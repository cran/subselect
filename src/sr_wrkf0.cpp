#include <math.h>
#include <assert.h>
#include "fullsqmi.h"
#include "fullmati.h"
#include "lagmvct.h"
#include "SR_vsmao.h"
#include "SR_vsda.h"
#include "vsmac.h"
#include "vsmav.h"

unsigned lastvar;
extern unsigned gcrt,flsti,flsts,maxcmp,*cmpl;  
extern fullvct **cmpv;

#ifndef NULL

#define NULL ((void *) 0)

// class  {
//	public:
//		template<class T>
//		operator T*() const  { return 0;  }
//	private:
//		void operator&() const;
// }    NULL;

#define NULLPOINT

#endif  


void 	crtwrksp()
{
	unsigned i;

	fpcnt0 = 0;
	switch (pcrt) {
		case MCB2:
			trs = S->v(1,1);  for(i=2;i<=p;i++) trs += S->v(i,i);
			if ((SW = new wrkspace(SRC,PTR,p,0,S,static_cast<symatrix *>(NULL),
                             static_cast<fullvct **>(NULL))) == static_cast<wrkspace *>(NULL) ) 
	          	  prmtend("SW");
			if ((IW = new wrkspace(INV,PTR,p,0,S,static_cast<symatrix *>(NULL),
                             static_cast<fullvct **>(NULL))) == static_cast<wrkspace *>(NULL) ) 
	          	  prmtend("IW");
			fpcnt0 = SQR(p)*(p+1)/2;
			break;
		case GCD:
			S->diag(maxcmp);
			for (i=0;i<q;i++) 
				scamult(sqrt(S->egval->v(cmpl[i])),S->egvct[cmpl[i]-1],cmpv[i]);
			if ( (SW = new wrkspace(SRC,V,p,q,S,static_cast<symatrix *>(NULL),cmpv) )
				== static_cast<wrkspace *>(NULL) ) prmtend("SW");			
			if ( (IW = new wrkspace(INV,V,p,q,S,static_cast<symatrix *>(NULL),cmpv) )
				== static_cast<wrkspace *>(NULL) ) prmtend("IW");			
			fpcnt0 = SQR(p)*(p+1)/2 + q*p*(p+1);
			break;
		case RV:
			multmat(S,S,S2);
			trs2 = S2->v(1,1);  for(i=2;i<=p;i++) trs2 += S2->v(i,i);
			if ( (SW = new wrkspace(SRC,NS2SM1,p,0,S,static_cast<symatrix *>(NULL),
                             static_cast<fullvct **>(NULL))) == static_cast<wrkspace *>(NULL) ) 
	          	  prmtend("SW");
			if ( (IW = new wrkspace(INV,NS2SM1,p,0,S,static_cast<symatrix *>(NULL),
                             static_cast<fullvct **>(NULL))) == static_cast<wrkspace *>(NULL) ) 
	          	  prmtend("IW");
			fpcnt0 = CUBE(p)+(SQR(p)+p)/2;
			break;
	}
}

void fsort()
{
	unsigned i,j,var;
	unsigned *iind,*sind,*fiind,*fsind;
	kspace *il,*sl,*ilt,*slt;

	il = IW->wrklst[flsti];
	sl = SW->wrklst[flsts];
	il->sort(il->srtv,fp+lp+1,p);
	lastvar = il->srtv[p-1];

	for (i=0;i<=flsts;i++)  {
		slt = SW->wrklst[i];
		for (j=fp+lp;j<slt->p;j++)  slt->srtv[j] = il->srtv[j];
	}
	for (i=0;i<=flsti-1;i++)  {
		ilt = IW->wrklst[i];
		for (j=fp+lp;j<ilt->p;j++)  ilt->srtv[j] = il->srtv[j];
	}

	iind = new unsigned[p-fp-lp];
	sind = new unsigned[p-fp-lp];
	fiind = new unsigned[p-fp-lp];
	fsind = new unsigned[p-fp-lp];
	for (i=0;i<p-fp-lp;i++)  {
		var = il->srtv[fp+lp+i];
		iind[i] = (fiind[i] = il->varp[var-1]) - lp;
		sind[i] = (fsind[i] = sl->varp[var-1]) - fp;
	}
	il->asgvar(fp+1,fp+lp+1,p-fp-lp,iind,fiind);
	sl->asgvar(lp+1,fp+lp+1,p-fp-lp,sind,fsind);
	delete[] iind;
	delete[] sind;
	delete[] fiind;
	delete[] fsind;

	for (i=0;i<=flsts;i++)  {
		slt = SW->wrklst[i];
		for (j=fp+lp;j<slt->p;j++)  slt->varp[slt->srtv[j]-1] = j+1;
	}
	for (i=0;i<=flsti;i++)  {
		ilt = IW->wrklst[i];
		for (j=fp+lp;j<ilt->p;j++)  ilt->varp[ilt->srtv[j]-1] = j+1;
	}
}
