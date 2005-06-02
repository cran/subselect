#include <cassert>
#include <ctime>
#include "Sscma.h"
#include "Subsets.h"
#include "Vsmabo.h"

using namespace leapsnbnds;

const real   NOBND = 1E99;

void prcksp1(wrkspace *,vind,vind,vind,vind);
void getactv(wrkspace *,vind,vind,vind);
void actvcnv(vind,vind,vind *,vind *);
bool leap(vind,real,const real *,vind,vind);

#ifdef COUNTING
void showcnt(int,long unsigned,vind);
#endif

extern vind maxdim,*prvks;    
extern int pcrttp;
extern enum {given,firstk} pcsets;
extern real   crub,*bndl;
extern psbstlist *bsts;

int cmp(const void *a,const void *b)
{
	vind ai = *(vind *)a - 1;
	vind bi = *(vind *)b - 1;
	if (pcrttp == MAXIMZ)  {
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

void actvcnv(vind p,vind p1,vind *v0,vind *v1)
{

	vind i,j,k;

	for (i=0,j=0,k=1;i<p1;i++,k++)  while (k<v0[i])  v1[j++] = k++;
	while (k<=p) v1[j++] = k++;
}

void getactv(wrkspace *w,vind t,vind k1,vind nv)
{
	subset *wlst = &w->subsetat(k1+1); 

	if (t == INV)  {
		if (wlst->getp() == p)  {  
			actvcnv(p,p-nv,wlst->getvar(),actv);
			{ for (vind i=0;i<nv;i++) actv[i] = wlst->getithvar(actv[i]-1)+1; }
		}
		else { 
			actvcnv(p-1,p-nv,wlst->getvar(),actv);
			{ for (vind i=0;i<nv-1;i++) actv[i] = wlst->getithvar(actv[i]-1)+1; }
			actv[nv-1] = lastvar;
		}
	}
	else
		{ for (vind i=0;i<nv;i++) actv[i] = wlst->getithvar(wlst->getvar(i+1)-1)+1; }
}

void prcksp1(wrkspace *w,vind tree,vind k0,vind k1,vind nv,vind vi)
{
	real crt;              
	sbset *prevsbset,*st;
	sbstlist::iterator ptprevsbset;
	
	assert(k0 > k1);
	if (k1 == 0)  w->pivot(vi,vi+1,vi,k0,0);
	else w->pivot(vi,vi+1,p-1,k0,k1);
	if (nv < mindim || nv > maxdim) return;
	#ifdef COUNTING
	++cntg;
	#endif
	crt = w->subsetat(k1+1).getdata().criterion();  

	if (pcrttp == MAXIMZ && crt < lbnd[nv-mindim]) return;	// Check if new subset is better than any of the
	if (pcrttp == MINIMZ && crt > ubnd[nv-mindim]) return;	// sets in current best set list for this dimension
	
	getactv(w,tree,k1,nv);
	st = csbset(nv,actv,crt);	
	psbstlist curlist = bsts[nv-mindim];

	if (sbsetcnt[nv-mindim] == ms)  {				// Remove and discard worst subset saved	   
		prevsbset = *(ptprevsbset=curlist->begin());
		curlist->erase(ptprevsbset);
		dsbset(prevsbset);
	}
	else sbsetcnt[nv-mindim]++;
	curlist->insert(st);					// Insert new subset in best sets list
	if (sbsetcnt[nv-mindim] == ms)				   
		if (pcrttp == MAXIMZ) lbnd[nv-mindim] = (*curlist->begin())->crt();
		else ubnd[nv-mindim] = (*curlist->begin())->crt();
	return;
}

bool leapsnbnds::prcksp(vind k,vind ks,vind nvs,vind nvi,vind fv)
{
	vind ks0;
	static vind k1;
	vind maxnvs,maxnvi,minnv,maxnv;  
	real   bnd = NOBND;                                     
	const  real*	cbnd;                         

	if (p-fv > 10) {
		if (clock()-btime > maxtime) return false;
	} 

	maxnvs = nvs+k; 
	maxnvi = nvi-1;
	{ for (vind minnvs=nvs+1,minnvi=nvi-k,i=0;i<k;minnvs++,minnvi++,i++)  {
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
	} }
	{ for (vind i=1;i<k;i++)   {
		minnv = nvs+k-i;
		maxnv = nvi-2;
		if (minnv <= maxdim && maxnv >= mindim) {
			if (minnv < mindim) minnv = mindim;
			if (maxnv > maxdim) maxnv = maxdim;
			if (minnv > maxnv) minnv = maxnv;
			bnd = IW->subsetat(i+1).getdata().criterion();  
			cbnd = IW->subsetat(i+1).getdata().getbnds();
			if (!leap(pcrttp,bnd,cbnd,minnv,maxnv) )
				if (!prcksp(i,prvks[i-1],nvs+k-i-1,nvi-1,p-i)) return false;
		}
	} }
	return true;
}


bool leap(vind dir,real crt,const real *crtcrr,vind minv,vind maxv)
{
	vind i;
	bool l;

	for (l=true,i=maxv;l&&i>=minv;i--) {
		if (crtcrr && i < maxv) crt -= crtcrr[i];	// Remove ith+1 parcel when using GCD with i (variable) PCs 
		if (dir == MAXIMZ && crt > lbnd[i-mindim]) l = false; 
		else if (dir == MINIMZ && crt < ubnd[i-mindim]) l = false;
	}
	return l;
}