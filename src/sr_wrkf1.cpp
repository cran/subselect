#include <assert.h>
#include <time.h>
#include "fullsqmi.h"
#include "fullmati.h"
#include "SR_vsda.h"
#include "lagmvct.h"
#include "SR_vsmao.h"
#include "vsmac.h"

void prcksp1(wrkspace *,unsigned,unsigned,unsigned,unsigned);
void getactv(wrkspace *,unsigned,unsigned,unsigned);
void showcnt(unsigned,unsigned,unsigned);

extern unsigned gcrt,pcrttp,lastvar,mindim,maxdim;  
extern double   crub,*bndl;
extern psbsq    *bsts;  

int cmp(const void *a,const void *b)
{
	unsigned ai,bi;

	ai = *(unsigned *)a - 1;
	bi = *(unsigned *)b - 1;
	if ( pcrt != MCB2 ) {
		if (Fl[ai] > Fl[bi]) return 1;
		else if (Fl[ai] < Fl[bi]) return -1;
			else return 0;
	}
	else {
		if (Fl[ai] < Fl[bi]) return 1;
		else if (Fl[ai] > Fl[bi]) return -1;
			else return 0;
	}
}

void actvcnv(unsigned p,unsigned p1,unsigned *v0,unsigned *v1)
{
	unsigned i,j,k;

	for (i=0,j=0,k=1;i<p1;i++,k++)  while (k<v0[i])  v1[j++] = k++;
	while (k<=p) v1[j++] = k++;
	return;
}

void getactv(wrkspace *w,unsigned t,unsigned k1,unsigned nv)
{
	unsigned i;
	kspace *wlst = w->wrklst[k1]; 

	if (t == INV)  {
		if (wlst-> p == p)  {  
			actvcnv(p,p-nv,wlst->var,actv);
			for (i=0;i<nv;i++) actv[i] = wlst->srtv[actv[i]-1];
		}
		else { 
			actvcnv(p-1,p-nv,wlst->var,actv);
			for (i=0;i<nv-1;i++) actv[i] = wlst->srtv[actv[i]-1];
			actv[nv-1] = lastvar;
		}
	}
	else
		for (i=0;i<nv;i++) actv[i] = wlst->srtv[wlst->var[i]-1];
	return;
}

void prcksp1(wrkspace *w,unsigned tree,unsigned k0,unsigned k1,unsigned nv,unsigned vi)
{

	double crt,bnd=NO;
	sbset *st=static_cast<sbset *>(NULL);
	sbsetqe *sq;

	assert(k0 > k1);
	if (k1 == 0)  w->pivot(vi,vi+1,vi,k0,0);
	else w->pivot(vi,vi+1,p-1,k0,k1);
	if (nv < mindim || nv > maxdim) return;
	++cntg;
	crt = w->wrklst[k1]->c;  
	if (pcrttp == MAX && crt < lbnd[nv-mindim]) return;  
	if (pcrttp == MIN && crt > ubnd[nv-mindim]) return;
	
	getactv(w,tree,k1,nv);
	switch (pcrt) {
		case MCB2:
			st = csbset(nv,0,actv,cmpcrMC2,indRM,crt);
			break;
		case GCD:
			if (q < 10) st = csbset(nv,q,actv,cmpcrGCD,indGCD,crt);
			else st = csbset(nv,q,actv,cmpcrGCD1,indGCD,crt);
			break;
		case RV:
			st = csbset(nv,p,actv,cmpcrRV,indRV,crt);
			break;
	}
	sq = csbsetqe(st);
	bnd = bsts[nv-mindim]->inserte(sq);
	if (bnd == NO) dsbsetqe(sq);
	if (bsts[nv-mindim]->full())
		if (pcrttp == MAX) lbnd[nv-mindim] = bnd;
		else ubnd[nv-mindim] = bnd;
	return;
}

int prcksp(unsigned k,unsigned ks,unsigned nvs,unsigned nvi,unsigned fv)
{
	unsigned i,ks0;
	static unsigned k1;
	unsigned minnvs,minnvi,maxnvs,maxnvi,minnv,maxnv;  
	double   bnd = NO;                              //  ,currtime;
	extern double btime,maxtime;

	if (p-fv > 15) 
		if (clock()-btime > maxtime) return FALSE;

	maxnvs = nvs+k; 
	maxnvi = nvi-1;
	for (minnvs=nvs+1,minnvi=nvi-k,i=0;i<k;minnvs++,minnvi++,i++)  {
		if (i == 0) ks0 = ks;
		else ks0 = k1;
		k1 = p-1-fv-i;
		if (k1 > 0) prvks[k1-1] = ks0;
		if (maxnvs >= mindim && minnvs <= maxdim) 
			if (minnvs == maxdim) prcksp1(SW,SRC,ks0,0,minnvs,fv+i);
			else prcksp1(SW,SRC,ks0,k1,minnvs,fv+i);
		if (maxnvi >= mindim && minnvi <= maxdim) 
			if (maxnvi == mindim) prcksp1(IW,INV,k,0,maxnvi,fv+i);
			else prcksp1(IW,INV,k,k1,maxnvi,fv+i);
	}
	for (i=1;i<k;i++)   {
		minnv = nvs+k-i;
		maxnv = nvi-2;
		if (minnv > maxnv) minnv = maxnv;
		bnd = IW->wrklst[i]->c;
		if (maxnv>=mindim && minnv<=maxdim && !leap(pcrttp,bnd,minnv,maxnv)          )
			if ( prcksp(i,prvks[i-1],nvs+k-i-1,nvi-1,p-i) == FALSE )
				return FALSE;
	}
	return TRUE;
}


char leap(unsigned dir,double crt,unsigned minv,unsigned maxv)
{
	unsigned i;
	char l;

	for (l=TRUE,i=minv;l&&i<=maxv;i++) {
		if (dir == MAX && crt > lbnd[i-mindim]) l = FALSE; 
		else if (dir == MIN && crt < ubnd[i-mindim]) l = FALSE;
	}
	return l;
}

