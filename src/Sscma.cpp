#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include "Sscma.h"
#include "Subsets.h"
#include "Vsmabo.h"

namespace leapsnbnds  {

vind		p,q,fp,lp,mindim,flst,lastvar,*actv;  
real		trs,trs2,*lbnd,*ubnd;
int		pcrt; 
long unsigned   ms,*sbsetcnt;                      
pcskept		pcsets;
std::string 	memmsg("Unable to find enough memory to run leaps with this data\n");


#ifdef COUNTING 
long unsigned cntg,fpcnt,fpcnt1;	//  Floating point operations counters   
#endif                

}

using namespace leapsnbnds;

subsetdata *idata,*fulldata;
globaldata *gidata,*gfulldata;
auxmemory *dataam;
vind ndim,maxdim,*prvks,*dmyv,*cmpl,*Flp,*ivlst,*ovlst;                
extern vind flsts,flsti;
int	pcrttp;                                         
wrkspace   *SW,*IW;   
real  c0,v0,*vc0,*Fl;                               				                   
long unsigned m0,maxsbqe,maxcmp;

#ifdef COUNTING
long unsigned  cntp,fpcnt0;   
#endif

psbst*		sbsarr;
psbstlist*	bsts;

void resetvar(void);
int getpcrt(char *,int *);
void initvlist(int *,int *,int *,int,int,int);
void cleanlists(void);
void fillres(vind,vind,long unsigned,int *,int *,real *,real *);
void saveset(psbstlist,int *,real *l,long unsigned,vind);
int trivialcmp(const void *,const void *);
void matasvcttranspose(int,int,int *);
void asgmemory(void);
void cleanup(void);

void fsort(void);

void resetvar()
{
	actv = prvks = dmyv = cmpl = Flp = ivlst = ovlst= 0;
	idata = fulldata = 0;
	gidata = gfulldata = 0;
	dataam = 0;
	lbnd = ubnd = vc0 = Fl = 0;
	c0 = v0 = 0.;                   
	sbsarr = 0;
	bsts = 0;
	sbsetind = 0;
	#ifdef COUNTING
	cntp = 0;   
	#endif
}

int getpcrt(char *st,int *fixed)
{
	if ( !strncmp(st,"rm",2) || !strncmp(st,"RM",2) || !strncmp(st,"Rm",2) )        return MCB2;    
	else if ( !strncmp(st,"rv",2) || !strncmp(st,"RV",2) || !strncmp(st,"Rv",2) )   return RV;
	else if ( !strncmp(st,"gcd",3) || !strncmp(st,"GCD",3) || !strncmp(st,"Gcd",3)) 
		{
			assert(fixed);
			if (*fixed) pcsets = given;
			else pcsets = firstk;
			return GCD;  
		}
	else return NOTFOUND;
}

void initvlist(int *ilist,int *olist,int *pcind,int ni,int no,int nind)
{
	try {
		if ( ni > 0) ivlst = new vind[ni];
		if ( no > 0) ovlst = new vind[no];
		if (pcrt != GCD) q = 0;
		else {
			if (pcsets == firstk) q = maxcmp = maxdim;
			else if ( (q=nind) == 0) {
				cleanlists();
				errmsg("Criterion GCD requires a non-empty list of S eigenvectors\n");
			}
			cmpl =  new vind[q];
		}
	}

	catch (std::bad_alloc)   {
		cleanlists();
		errmsg(memmsg);
	}

	fp = ni;
	{ for (int i=0;i<ni;i++) ivlst[i] = ilist[i]; }
	lp = no;
	{ for (int i=0;i<no;i++) ovlst[i] = olist[i]; }
	if (q > 0)  {
		{ for (int i=1;i<=(int)q;i++)  {
			if (pcsets == firstk) cmpl[i-1] = i;
			else {
				if (i == 1) maxcmp = cmpl[0] = pcind[0];
				else if ( (cmpl[i-1] = pcind[i-1]) > maxcmp )   maxcmp = cmpl[i-1];
			}
		} }
	}
}

void cleanlists() {
	delete[] ivlst;
	delete[] ovlst;
	delete[] cmpl;
}


bool sscma(bool fullwrksp,subsetdata *nullsetdt,subsetdata *fullsetdt)
{
	SW = new wrkspace(fullwrksp,SRC,p,nullsetdt);
	IW = new wrkspace(fullwrksp,INV,p,fullsetdt);
	flst = flsts;
	#ifdef COUNTING  
	fpcnt1 = 0;
	#endif
	if (p > fp+lp+1) fsort();
	else lastvar = IW->subsetat(flsti+1).getithvar(p-1)+1;
	assert(lastvar > 0 && lastvar <= p);
	if (fp > 0 && fp == mindim) savfrst();
	if (maxdim == p-lp) savfull();
	#ifdef COUNTING  
	fpcnt = 0;
	#endif
	if (p <= fp+lp+1 || prcksp(flst,flst,fp,p-lp,fp+lp+1) )  return  true;
	else return false;
}

void fsort()
{
	vind var,*iind,*sind;                  
	vind sskipvar=0,iskipvar=0;
	subset *il,*sl,*ilt,*slt;

	il = &IW->subsetat(flsti+1);
	sl = &SW->subsetat(flsts+1);

	il->sort(fp+lp+1,p);
	lastvar = il->getithvar(p-1)+1;

	{ if (SW != NULL) for (vind i=1;i<=flsts+1;i++)  {
		slt = &SW->subsetat(i);
		for (vind j=fp+lp;j<p;j++)  slt->setithvar(j,il->getithvar(j));
	} }
	{ for (vind i=1;i<=flsti;i++)  {
		ilt = &IW->subsetat(i);
		for (vind j=fp+lp;j<p;j++)  ilt->setithvar(j,il->getithvar(j));
	} }

	iind = new vind[p-fp-lp];
	sind = new vind[p-fp-lp];

	for (vind i=0;i<p-fp-lp;i++)  {
		var = il->getithvar(fp+lp+i);
		if (lp == 0 && fp > 0) {
			iskipvar = fp;
			iind[i] = var+1;
		}
		else iind[i] = il->getvarp(var)-fp-lp+1;
		if (fp == 0 && lp > 0) {
			sskipvar = lp;
			sind[i] = var+1;
		}
		else sind[i] = sl->getvarp(var)-fp-lp+1;
	}

	il->asgvar(iskipvar,p-fp-lp,iind);
	sl->asgvar(sskipvar,p-fp-lp,sind);

	delete[] iind;
	delete[] sind;

	if (SW != NULL) 
		{	for (unsigned i=1;i<=flsts+1;i++)  {
				slt = &SW->subsetat(i);
				for (unsigned j=fp+lp;j<slt->getp();j++)  slt->setvarp(slt->getithvar(j),j);
			}
	}
	{ for (vind i=1;i<=flsti+1;i++)  {
		ilt = &IW->subsetat(i);
		for (vind j=fp+lp;j<ilt->getp();j++)  ilt->setvarp(ilt->getithvar(j),j);
	} }
}

void fillres(vind fk,vind nk,long unsigned ns,int *bst,int *st,real *bvl,real *vl)
{
	vind i,lk=fk+nk-1;
	int *stmatp;
	
	for (i=0;i<nk;i++)  {
		stmatp = &(st[i*ns*lk]);
		saveset(bsts[i],stmatp,&(vl[i*ns]),ns,lk);
		matasvcttranspose(ns,lk,stmatp);
		saveset(bsts[i],&(bst[i*lk]),&(bvl[i]),1,lk);
	}
	matasvcttranspose(nk,lk,bst);
	return;
}

void saveset(psbstlist pset,int *bvar,real *bcrtval,long unsigned nel,vind dim)
{
	int i=0,j,*var;

	for (sbstlist::reverse_iterator qep=pset->rbegin();i<nel&&qep!=pset->rend();++qep)  {
		var = &(bvar[(i++)*dim]);
		for (j=0;j<(*qep)->nvar();j++) var[j] = (*qep)->actvar()[j];  
		qsort(static_cast<void *>(var),(*qep)->nvar(),sizeof(*var),trivialcmp);
		for (j=(*qep)->nvar();j<dim;j++) var[j] = 0; 
		switch (pcrt)  {
			case MCB2:
				*bcrtval++ = sqrt(1. - (*qep)->crt()/trs);
				break;
			case GCD:
				if (pcsets == given) *bcrtval++ = (*qep)->crt()/sqrt(q*(*qep)->nvar());
				else if (pcsets == firstk) *bcrtval++ = (*qep)->crt()/(*qep)->nvar();
				break;
			case RV:
				*bcrtval++ = sqrt((*qep)->crt()/trs2);
				break;
		}
	}
	for (i=pset->size();i<nel;i++)  {
		for (vind j=0;j<dim;j++) bvar[i*dim+j] = 0;
		*bcrtval++ = 0.;
	} 
}

int trivialcmp(const void *a,const void *b)
{
	int ai,bi;

	ai = *(int*)a;
	bi = *(int*)b;
	if (ai > bi)  return  1;
	else if (ai < bi)  return -1;
		 else return 0;
}

void matasvcttranspose(int m,int n,int *data)
{
	vind mn=m*n;
	int  *tmp=0;
	
	try { tmp = new int[mn]; }
	catch (std::bad_alloc)   {
		cleanup();
		errmsg(memmsg);
	}
	{ for (vind i=0;i<m;i++)
		for (vind j=0;j<n;j++)  tmp[i+j*m] = data[i*n+j]; }
	{ for (vind i=0;i<mn;i++) data[i] = tmp[i]; }
	delete[] tmp;
}

void asgmemory()
{
	try  {

		actv = new vind[p];
		dmyv = new vind[p];
		Fl = new real[p];
		Flp = new vind[p];

		if (ms > 0) {
			bsts = new psbstlist[ndim];
			sbsetcnt =new long unsigned[ndim];
			{for (vind i=0;i<ndim;i++) sbsetcnt[i] = 0; }
		}
		if (ndim == p-fp-lp+1) maxsbst = maxsbqe = ms*(ndim-1)+2;
		else maxsbst = maxsbqe = ms*ndim+3;
		sbsarr = new psbst[maxsbst];
		{ for (long unsigned i=0;i<maxsbst;i++) sbsarr[i]=new sbset(i,p);  }

		if ( pcrt == MCB2 ) {
			pcrttp = MINIMZ;
			if (ms > 0) 
				{ for (vind i=0;i<ndim;i++) bsts[i] = new sbstlist(sbstsort::descending);  }
				if (ndim == p-fp-lp+1) { 
					ubnd = new real[ndim-1];
					{ for (vind i=0;i<ndim-1;i++) 
						ubnd[i] = std::numeric_limits<real>::infinity();  } 
				}
				else {
					ubnd = new real[ndim];
					{ for (vind i=0;i<ndim;i++) 
						ubnd[i] = std::numeric_limits<real>::infinity();;  }
				}
			lbnd = 0;
		}
		else  {
			pcrttp = MAXIMZ;
			if (ms > 0) {
				{  for (vind i=0;i<ndim;i++) bsts[i] = new sbstlist(sbstsort::ascending);  }
				if (ndim == p-fp-lp+1) { 
					lbnd = new real[ndim-1];
					{ for (vind i=0;i<ndim-1;i++) lbnd[i] = 0.;  }
				}
				else {
					lbnd = new real[ndim];
					{ for (vind i=0;i<ndim;i++) lbnd[i] = 0.;  }
				}
			}
			ubnd = 0;
		}
		if ( pcrt == GCD && pcsets == firstk)  vc0 = new real[q];  
		else vc0 = 0;
		prvks = new vind[p-1];
	}

	catch (std::bad_alloc)   {
		cleanup();
		errmsg(memmsg);
	}
}

void cleanup(void)
{
	long unsigned i;
	
	delete   SW;
	delete   IW;
	delete[] actv;
	delete[] dmyv;
	delete[] Fl;
	delete[] Flp;
	if (bsts != 0) {
		for (i=0;i<ndim;i++) delete bsts[i];
		delete[] bsts;
	}
	delete[] sbsetcnt; 
	delete[] lbnd;
	delete[] ubnd;
	delete[] prvks;
	if ( sbsarr != 0) {
		for (i=0;i< maxsbst;i++) delete sbsarr[i];
		delete[] sbsarr;
	}
	delete[] vc0;
	delete idata;
	delete fulldata;
	delete gidata;
	delete gfulldata;
	delete dataam;
	cleanlists();
}
