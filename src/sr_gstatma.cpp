#include <assert.h>
#include <math.h>
#include <iostream.h>
#include <stdlib.h>
#include "fullsqmi.h"
#include "fullmati.h"
#include "SR_vsda.h"
#include "lagmvct.h"
#include "SR_vsmao.h"
#include "vsmac.h"

extern unsigned gcrt,g,mttp,nx,ny,flsts,flst,lastvar,mindim,*inactv;
extern psbsq      *bsts;  


void savfrst()
{
	unsigned i;
	double mcbc,rv,gcd;
	kspace *fset;
	sbset  *st=NULL;

	assert(flst >= 0 && flst < p+1);
	fset = SW->wrklst[flst];
	for (i=lp;i<p-1;i++) {
		actv[i-lp] = fset->srtv[i];
		assert(actv[i-lp] > 0 && actv[i-lp] <= p);
	}
	assert(lastvar > 0 && lastvar <= p);  
	actv[p-lp-1] = lastvar;
	switch (pcrt) {
		case MCB2:
			mcbc = fset->c;
			st = csbset(fp,0,actv,cmpcrMC2,indRM,mcbc);
			break;
		case GCD:
			gcd = fset->v;
			if (q < 10) st = csbset(fp,q,actv,cmpcrGCD,indGCD,gcd);
			else st = csbset(fp,q,actv,cmpcrGCD1,indGCD,gcd);
			break;
		case RV:
			rv = fset->c;
			st = csbset(fp,p,actv,cmpcrRV,indRV,rv);
			break;
	}
	bsts[fp-1]->inserte(csbsetqe(st));
	cntg++;
	return;
}

void savfull(void)
{
	unsigned i;
	double mcbc,gcd,trsm1sg,trsm1s2sqr;
	sbset *st=NULL;
	kspace *wlst; 

	assert(flst >= 0 && flst <= p);
	wlst = IW->wrklst[flst];
	for (i=0;i<fp;i++) {
		actv[i] = wlst->srtv[i];
		assert(actv[i] > 0 && actv[i] <= p);
	}
	for (i=fp+lp;i<p-1;i++) {
		actv[i-lp] = wlst->srtv[i];
		assert(actv[i-lp] > 0 && actv[i-lp] <= p);
	}
	assert(lastvar > 0 && lastvar <= p);  
	actv[p-lp-1] = lastvar;
//  if (lp == 0)  {
//		cout << "VARS:   ";
//		for (i=0;i<p-lp;i++) {
//			if (!(i%VARPL) && i) cout << "\n        ";
//			cout << "  " << varn[actv[i]-1];
//		}
//	}
	switch (pcrt) {
		case MCB2:
			mcbc = wlst->c;
//			if (mttp == VC && lp == 0) cout << "\nTr S = " << trs << "\n";
			st = csbset(p-lp,0,actv,cmpcrMC2,indRM,mcbc);
			break;
		case GCD:
//			cout << "\n";
			trsm1sg = wlst->c;
			gcd = trsm1sg/sqrt(p*q);
//			cout << "\nGCD[G(" << q << "),X] = " << gcd << "\n";
			if (q < 10) st = csbset(p-lp,q,actv,cmpcrGCD,indGCD,trsm1sg);
			else st = csbset(p-lp,q,actv,cmpcrGCD1,indGCD,trsm1sg);
			break;
		case RV:
//			if (lp == 0)
//				if (mttp == CORR) cout << "\nTr R^2 = " << trs2 << "\n";
//				else cout << "\nTr S^2 = " << trs2 << "\n";
			trsm1s2sqr = wlst->c;
			st = csbset(p-lp,p,actv,cmpcrRV,indRV,trsm1s2sqr);
			break;
	}
	bsts[p-lp-mindim]->inserte(csbsetqe(st));
	cntg++;
	return;
}


