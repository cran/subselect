#include <iomanip.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "fullsqmi.h"
#include "fullmati.h"
#include "lagmvct.h"
#include "SR_vsmao.h"
#include "SR_vsda.h"
#include "vsmac.h"
#include "vsmav.h"

const unsigned MAT = 999999;

typedef lagvct*  plagvct;
typedef cndvct*  pcndvct;
extern  unsigned g,back,*frcv;
unsigned flsts,flsti;

double csymip(unsigned,symatrixb *,symatrixb *,unsigned *,unsigned);
double ctrsq(unsigned,fwrksqm *,unsigned *,unsigned);
void initzero(fwrksqm *);
void us2sm1(fwrksqm *,fwrksqm *,cndvct **,unsigned *,unsigned *,unsigned,unsigned,unsigned,double *,double *);

double cndvct::v(unsigned  i)
{
	if (mt != NULL)  return mt->v(mr,i);
	else return fullvct::v(i-lag());
}

kspace::kspace(unsigned crt,unsigned tnv,unsigned pnv,unsigned rh,unsigned nvo)
  :  ind(crt), p(tnv), k(pnv), r(rh), p0(nvo), nvar(0)
{
	unsigned i;

	assgnmem();
	for (i=0;i<p;i++)  srtv[i] = varp[i] = i+1;
}

kspace::kspace(unsigned *ivar,unsigned crt,unsigned tnv,unsigned pnv,unsigned rh,unsigned nvo)
  :  ind(crt), p(tnv), k(pnv), r(rh), p0(nvo), nvar(0)
{
	unsigned i;

	assgnmem();
	for (i=0;i<tnv;i++)  {
		srtv[i]  = ivar[i];
		varp[ivar[i]-1] = i+1;
	}
}

void kspace::assgnmem()
{
	unsigned i;

	if (k<p) var = new unsigned[p-k];
	else var = NULL;
	srtv = new unsigned[p];
	varp = new unsigned[p0];
	if (ind == V)  vst = NULL;
	else vst = new unsigned[p0];
	if (ind != NS2SM1) s2m1 = NULL;
	else s2m1 = new fwrksqm(p);
	if (k > 0)  {
		if (ind != V) ve = NULL;
		else {
			ve = new plagvct[r];
			for (i=0;i<r;i++) ve[i] = new lagvct(k,p);
		}
		e = new kmatrix(k,p);
		if (ind != PTR)  ovct = NULL;
		else {
			ovct = new pcndvct[p0];
			for (i=0;i<p0;i++) ovct[i] = new cndvct(k,p,e,i+1);
		}
		if (ind != NS2SM1)  ivct = NULL;
		else {
			ivct = new pcndvct[p0];
			for (i=0;i<p0;i++)  ivct[i] = new cndvct(k,p,e,i+1);
		}
	}
	else  {
		ve = NULL;
		e = NULL;
		ovct = ivct = NULL;
	}
	switch (ind) {
		case V:
			c = v = UNKNW;
			if (e != NULL) e->det = NOTDF;
			break;
		case PTR:
			c = UNKNW;
			if (e != NULL) e->det = UNKNW;
			v = NOTDF;
			break;
		case NS2SM1:
			c = s2m1->det = UNKNW;
			if (e != NULL) e->det = UNKNW;
			v = NOTDF;
			break;
	}
}

kspace::~kspace()
{
	unsigned i;

	delete[] var;
	delete[] srtv;
	delete[] varp;
	delete e;
	delete s2m1;
	if (ve != NULL)   {
		for (i=0;i<r;i++) delete ve[i];
		delete[] ve;
	}
	delete[] vst;
	if (ovct != NULL)   {
		for (i=0;i<p0;i++) delete ovct[i];
		delete[] ovct;
	}
	if (ivct != NULL)   {
		for (i=0;i<p0;i++) delete ivct[i];
		delete[] ivct;
	}
}

void kspace::mvtofront(unsigned *l,unsigned nele)
{
	unsigned i,j;

	assert(l && nele > 0 && nele <= p);
	for (i=0;i<nele;i++) {
		e->swpvar(i+1,e->getvarp(l[i]));
		if (ind == NS2SM1) {
			s2m1->swprows(i+1,s2m1->getrowp(l[i]));
			s2m1->swpcols(i+1,s2m1->getcolp(l[i]));
			for (j=0;j<p0;j++) ivct[j]->swprows(i+1,ivct[j]->getrowp(l[i]));
		}
		if (ind == PTR) 
			for (j=0;j<p0;j++) ovct[j]->swprows(i+1,ovct[j]->getrowp(l[i]));
		if (ve != NULL) 
			for (j=0;j<r;j++) ve[j]->swprows(i+1,ve[j]->getrowp(l[i]));
	}
}

void kspace::asgvar(unsigned flagv,unsigned ffullv,unsigned nv,unsigned *lagv,unsigned *fullv)
{
	unsigned i;

	assert(lagv && flagv > 0 && nv > 0 && flagv+nv-1 <= k );
	assert( ind != NS2SM1 || (fullv && ffullv > 0 && ffullv+nv-1 <= p) );
	e->asgvar(flagv,nv,lagv);
	switch (ind)  {
		case V:
			for (i=0;i<r;i++)  ve[i]->asgrows(flagv,nv,lagv);
			break;
		case PTR:
			for (i=0;i<p0;i++)  ovct[i]->asgrows(flagv,nv,lagv);
			break;
		case NS2SM1:
			for (i=0;i<p0;i++)  ivct[i]->asgrows(flagv,nv,lagv);
			s2m1->asgrows(ffullv,nv,fullv);
			s2m1->asgcols(ffullv,nv,fullv);
			break;
	}
}

void kspace::sort(unsigned *l,unsigned fv,unsigned lv)
{
	unsigned i;

	assert(l && fv > 0 && lv > fv && lv <= p);
	for (i=0;i<=lv-fv;i++)  {
		Fl[i] = Ftolv(fv+i);
		Flp[srtv[fv+i-1]-1] = i+1;
		dmyv[i] = i+1;
	}
	qsort((void *)dmyv,lv-fv+1,sizeof(*dmyv),cmp);
	for (i=fv;i<=lv;i++) dmyv[i-fv] = l[dmyv[i-fv]+fv-2];
	for (i=fv;i<=lv;i++) l[i-1] = dmyv[i-fv];
}

double   *cv,*tv;
fwrksqm  *m1t;

double kspace::Ftolv(unsigned vi)
{
	unsigned i;
	double   dv,dt,e1,tt;

	assert(vi >= 1 && vi <= p);
	switch (ind)  {
		case V:
			dv = 0.;
			e1 = e->v(vi,vi);
			for (i=0;i<r;i++) dv += pow(ve[i]->v(vi),2)/e1;
			fpcnt1 += 2*r;
			return v0+dv;
		case PTR:
			e1 = e->v(vi,vi);
			if (vst[vi-1] == VOUT) dt = -e1;
			else dt = -1./e1;
			for (i=1;i<vi;i++) if (vst[i-1] == VOUT) {
				dt -= pow(ovct[i-1]->v(vi),2)/e1;
				fpcnt1 += 2;
			}
			for (i=vi+1;i<=p;i++) if (vst[i-1] == VOUT) {
				dt -= pow(e->v(i,vi),2)/e1;
				fpcnt1 += 2;
			}
			return c0 + dt;
		case NS2SM1:
			vst[vi-1] = VOUT;
			e1 = e->v(vi,vi);
			for (i=1;i<vi;i++)	if (vst[i-1] == VIN) 
				cv[i-1] = ivct[i-1]->v(vi)/e1; 
			cv[vi-1] = 1./e1;
			for (i=vi+1;i<=p;i++)	cv[i-1] = e->v(i,vi)/e1;
			us2sm1(s2m1,m1t,ivct,vst,srtv,vi,vi+1,vi,cv,tv);
			tt = ctrsq(p,m1t,vst,VIN);
			vst[vi-1] = VIN;
			fpcnt1 += p;
			return tt;
	}
	return 0;
}

wrkspace::wrkspace(unsigned tp,unsigned crt,unsigned nv,unsigned rh,symatrix *m1,symatrix *m2,fullvct **v)
   :  ind(crt), p (nv), r(rh), E(m1), T(m2),  hv(v)
{
	unsigned j,*tlst=NULL;                
	int j1;
	kspace *ispc;
	extern unsigned *ivlst,*ovlst;

	if (ind == V) ve = new double[r];
	else ve = NULL;
	if (ind == PTR || ind == NS2SM1) cndv = new double[p];
	else cndv = NULL;
	if (ind == NS2SM1) tmpv = new double[p];
	else tmpv = NULL;

	switch (tp)  {
		case SRC:
			flsts = p-lp-1;
			wrklst = new pkspc[flsts+1];
			if (fp+lp > 0)  {
				tlst   = new unsigned[p];
				frontlsts(ovlst,ivlst,lp,fp,tlst);
				ispc = wrklst[flsts] = new kspace(tlst,ind,p,p,r,p);
			}
			else ispc = wrklst[flsts] = new kspace(ind,p,p,r,p);
			cpmatdt(E,ispc->e);	
			if (ind == NS2SM1) initzero(ispc->s2m1);
			if (ispc->ve != NULL) for (j=0;j<r;j++) cpvectdt(hv[j],ispc->ve[j]);
			if (lp+fp > 0) ispc->mvtofront(tlst,fp+lp);
			if (ispc->vst != NULL)  for (j=0;j<p;j++) ispc->vst[j] = VOUT;
			if (ind == V) ispc->v = 0.;
			if (ind == PTR) {
				ispc->c = ispc->e->v(1,1);
				for (j=2;j<=p;j++)  	ispc->c += ispc->e->v(j,j);
			}
			for (j=1;j<=fp;j++)  {
				if (tlst == NULL) wrklst[flsts-j] = new kspace(ind,p,p-j,r,p);
				else wrklst[flsts-j] = new kspace(tlst,ind,p,p-j,r,p);
				if (flsts > j) pivot(lp+j,lp+j+1,p,flsts+1-j,flsts-j);
				else pivot(lp+j,lp+j+1,p-1,flsts+1-j,0);
				delete wrklst[flsts+1-j];
			}
			nwl = (flsts -= fp) + 1;
			break;
		case INV:
			flsti = p-fp-1;
			wrklst = new pkspc[flsti+1];
			if (fp+lp > 0)  {
				tlst   = new unsigned[p];
				frontlsts(ivlst,ovlst,fp,lp,tlst);
				ispc = wrklst[flsti] = new kspace(tlst,ind,p,p,r,p);
			}
			else ispc = wrklst[flsti] = new kspace(ind,p,p,r,p);
			cpmatdt(E,ispc->e);
			ispc->e->siminvrp();
			E->det = -(ispc->e->det = -1./ispc->e->det);
			if (ind == NS2SM1) cpmatdt(E,ispc->s2m1);
			if (ispc->ve != NULL) for (j=0;j<r;j++) rightmult(ispc->e,hv[j],ispc->ve[j]);
			if (ispc->vst != NULL)  for (j=0;j<p;j++) ispc->vst[j] = VIN;
			if (ind == V) {
				ispc->v = 0.;
				for (j=0;j<r;j++) ispc->v -= vectprod(hv[j],ispc->ve[j]);
			}
			if (lp+fp > 0) ispc->mvtofront(tlst,fp+lp);
			if (ind == V)      v0 = ispc->c = ispc->v;
			if (ind == PTR)    c0 = ispc->c = 0.;
			if (ind == NS2SM1) c0 = ispc->c = trs2;
			for (j=1;j<=lp;j++)  {
				if (tlst == NULL) wrklst[flsti-j] = new kspace(ind,p,p-j,r,p);
				else wrklst[flsti-j] = new kspace(tlst,ind,p,p-j,r,p);
				pivot(fp+j,fp+j+1,p,flsti+1-j,flsti-j);
				delete wrklst[flsti+1-j];
			}
			nwl = (flsti -= lp) + 1;
			break;
	}
	for (j1=nwl-2;j1>=0;j1--) {
		if (ind != NS2SM1)
			if (tlst == NULL) wrklst[j1] = new kspace(ind,p-1,j1,r,p);
			else wrklst[j1] = new kspace(tlst,ind,p-1,j1,r,p);
		else  
			if (tlst == NULL) wrklst[j1] = new kspace(ind,p,j1+1,r,p);
			else wrklst[j1] = new kspace(tlst,ind,p,j1+1,r,p);
	}
	delete[] tlst;
}

void wrkspace::frontlsts(unsigned *l1,unsigned *l2,unsigned nl1,unsigned nl2,unsigned *ol)
{
	unsigned i,tmp,*elep;
	
	elep = new unsigned[p];
	for (i=0;i<p;i++)  ol[i] = elep[i] = i+1;
	for (i=0;i<nl1;i++) {
		tmp = ol[i];
		ol[i] = l1[i];
		ol[elep[l1[i]-1]-1] = tmp;
		elep[tmp-1] = elep[l1[i]-1];
		elep[l1[i]-1] = i+1;
	}
	for (i=0;i<nl2;i++) {
		tmp = ol[nl1+i];
		ol[nl1+i] = l2[i];
		ol[elep[l2[i]-1]-1] = tmp;
		elep[tmp-1] = elep[l2[i]-1];
		elep[l2[i]-1] = nl1+i+1;
	}
	delete[] elep;
}

wrkspace::~wrkspace(void)
{
	unsigned j;

	delete[] ve;
	delete[] cndv;
	delete[] tmpv;
	for (j=0;j<nwl;j++)  delete wrklst[j];
	delete[] wrklst;
}

void wrkspace::pivot(unsigned vp,unsigned v1,unsigned vl,unsigned li,unsigned lo)
{
	unsigned j;
	double e;
	kspace *isi,*iso;

	isi = wrklst[li];
	iso = wrklst[lo];
	for (j=0;j<isi->nvar;j++) iso->var[j] = isi->var[j];
	iso->var[isi->nvar] = vp;
	iso->nvar = isi->nvar+1;
	if (isi->vst != NULL)  {
		for (j=1;j<=p;j++)  if (j != vp) iso->vst[j-1] = isi->vst[j-1];
		if (isi->vst[vp-1] == VIN) iso->vst[vp-1] = VOUT;
		else iso->vst[vp-1] = VIN;
	}
	e = isi->e->v(vp,vp);
	if (ind == V)  {
		iso->v = isi->v;
		for (j=0;j<r;j++) {
			ve[j] = isi->ve[j]->v(vp)/e;
			iso->v += ve[j]*isi->ve[j]->v(vp);
		}
		fpcnt += 2*r;
	}
	if (ind == PTR) {
		if (iso->vst[vp-1] == VOUT) { iso->c = isi->c - 1/e; fpcnt++; }
		else iso->c = isi->c - e;
		for (j=1;j<v1;j++) if (j != vp && iso->vst[j-1] == VOUT)  {
			cndv[j-1] = isi->ovct[j-1]->v(vp)/e;
			iso->c -=  cndv[j-1]*isi->ovct[j-1]->v(vp);
			fpcnt += 2;
		}
		for (j=vl+1;j<=p;j++) if (iso->vst[j-1] == VOUT)  {
			cndv[j-1] = isi->ovct[j-1]->v(vp)/e;
			iso->c -=  cndv[j-1]*isi->ovct[j-1]->v(vp);
			fpcnt += 2;
		}
	}
	if (ind == NS2SM1) {
		for (j=1;j<v1;j++) if (j != vp && iso->vst[j-1] == VIN)
			{ cndv[j-1] = isi->ivct[j-1]->v(vp)/e;  fpcnt ++;  }
		for (j=vl+1;j<=p;j++) if (iso->vst[j-1] == VIN)
		  { cndv[j-1] = isi->ivct[j-1]->v(vp)/e;  fpcnt ++;  }
	}

	if (lo != 0)  {
		pivot(isi->e,iso->e,vp,v1,vl);
		if (ind == V)
			for (j=0;j<r;j++) pivot(isi->ve[j],iso->ve[j],isi->e,ve[j],vp,v1,vl);
		if (ind == PTR)  {
			for (j=1;j<v1;j++)  if (j != vp && iso->vst[j-1] == VOUT)  {
				pivot(isi->ovct[j-1],iso->ovct[j-1],isi->e,cndv[j-1],vp,v1,vl);
				iso->ovct[j-1]->mt = NULL;
			}
			if (iso->vst[vp-1] == VOUT) {
				for (j=v1;j<=vl;j++) iso->ovct[vp-1]->put(-isi->ovct[vp-1]->v(j)/e,j);
				iso->ovct[vp-1]->mt = NULL;
			}
			for (j=vl+1;j<=p;j++)  if (iso->vst[j-1] == VOUT) {
				pivot(isi->ovct[j-1],iso->ovct[j-1],isi->e,cndv[j-1],vp,v1,vl);
				iso->ovct[j-1]->mt = NULL;
			}
			for (j=v1;j<=vl;j++)  if (iso->vst[j-1] == VOUT)
				iso->c += ( iso->e->v(j,j) - isi->e->v(j,j) );
		}
		if (ind == NS2SM1)  {
			for (j=1;j<v1;j++)  if (j != vp && iso->vst[j-1] == VIN)  {
				pivot(isi->ivct[j-1],iso->ivct[j-1],isi->e,cndv[j-1],vp,v1,vl);
				iso->ivct[j-1]->mt = NULL;
			}
			if (iso->vst[vp-1] == VIN) {
				for (j=v1;j<=vl;j++) iso->ivct[vp-1]->put(-isi->ivct[vp-1]->v(j)/e,j);
				fpcnt += (vl-v1+1);
				iso->ivct[vp-1]->mt = NULL;
			}
			for (j=vl+1;j<=p;j++)  if (iso->vst[j-1] == VIN) {
				pivot(isi->ivct[j-1],iso->ivct[j-1],isi->e,cndv[j-1],vp,v1,vl);
				iso->ivct[j-1]->mt = NULL;
			}
		}
	}
	if (ind == V) iso->c = iso->v;
	if (ind == NS2SM1) {
		if (iso->vst[vp-1] == VIN)  cndv[vp-1] = 1./e;
		for (j=v1;j<=vl;j++)  if (iso->vst[j-1] == VIN) {
			cndv[j-1] = isi->e->v(j,vp)/e;
			fpcnt += (vl-v1+1);
		}
		us2sm1(isi->s2m1,iso->s2m1,isi->ivct,iso->vst,iso->srtv,vp,v1,vl,cndv,tmpv);
		iso->c = ctrsq(p,iso->s2m1,iso->vst,VIN);
	}

	return;
}

void wrkspace::pivot(symatrixb *im,symatrixb *om,unsigned vp,unsigned v1,unsigned vl)
{
	unsigned i,j;
	double t,t1;

	t = 1./im->v(vp,vp);
	if (vl > v1) for (i=0;i<vl-v1+1;i++)  om->vi[i] = om->iv[i] = i+1;
	for (i=v1;i<=vl;i++)  {
		t1 = im->v(i,vp)*t;
		for (j=i;j<=vl;j++) om->put(im->v(i,j)-t1*im->v(vp,j),i,j);
	}
	fpcnt += (vl-v1+1)*(vl-v1+4)/2;
	return;
}

void wrkspace::cpivot(symatrixb *im,symatrixb *om,unsigned vp,unsigned v1,unsigned vl,unsigned cnd,unsigned *lst)
{
	unsigned j,k;
	double t,t1=0.;

	t = im->v(vp,vp);
	for (j=1;j<=p;j++) if ( (j >= v1 && j <= vl) || lst[j-1] == cnd)  {
		if (j == vp) { om->put(-1/t,j,j); fpcnt++; }
		else {
			t1 = im->v(j,vp)/t;
			om->put(im->v(j,j)-t1*im->v(vp,j),j,j);
			fpcnt += 2;
		}
		for (k=j+1;k<=p;k++)  if ( (k >= v1 && k <= vl) || lst[k-1] == cnd) {
			if (j == vp) { om->put(-im->v(j,k)/t,j,k); fpcnt++; }
			if (k == vp) om->put(-t1,j,k);
			if (j != vp && k != vp)
				{  om->put(im->v(j,k)-t1*im->v(vp,k),j,k); fpcnt++;  }
		}
	}
}


void wrkspace::pivot(fullvct *iv,fullvct *ov,kmatrix *im,double t1,unsigned vp,unsigned v1,unsigned vl)
{
	unsigned j;

	if (vl >= v1) for (j=0;j<vl-v1+1;j++)  ov->ri[j] = ov->ir[j] = j+1;
	for (j=v1;j<=vl;j++) ov->put(iv->v(j)-t1*im->v(vp,j),j);
	fpcnt += vl-v1+1;
	return;
}

double csymip(unsigned p,symatrixb *m1,symatrixb *m2,unsigned *lst,unsigned cnd)
{
	double tmp;
	unsigned i,j;

	tmp = 0.;
	for (i=1;i<=p;i++)  if (lst[i-1] == cnd)  {
		tmp += m1->v(i,i)*m2->v(i,i);  fpcnt++;
		for (j=i+1;j<=p;j++)  if (lst[j-1] == cnd)
			{  tmp += 2*m1->v(i,j)*m2->v(i,j);	 fpcnt++;  }
	}
	return tmp;
}

double ctrsq(unsigned p,fwrksqm *m,unsigned *lst,unsigned cnd)
{
	double tmp;
	unsigned i,j;

	tmp = 0.;
	for (i=1;i<=p;i++)  if (lst[i-1] == cnd) {
		tmp += pow(m->v(i,i),2);   fpcnt++;
		for (j=1;j<i;j++)  if (lst[j-1] == cnd)
			{ tmp += 2*m->v(i,j)*m->v(j,i);  fpcnt++; }
	}
	return tmp;
}

void initzero(fwrksqm *m)
{
	unsigned i,j;

	for (i=1;i<=m->n;i++)
		for (j=1;j<=m->p;j++) m->put(0.,i,j);
}

void us2sm1(fwrksqm *m1,fwrksqm *m2,cndvct **wcol,unsigned *lst,unsigned *slst,
				unsigned u,unsigned u1,unsigned ul,double *fl,double *tmpv)
{
	unsigned a,i,j;
	double val;

	for (j=1;j<=p;j++) {
		if ( (j < u1 || j > ul) && lst[j-1] != VIN) continue;
		tmpv[j-1] = 0.;
		for (a=1;a<=p;a++) if (lst[a-1] == VIN && a != u)
		 { tmpv[j-1] += -wcol[a-1]->v(u)*S2->v(slst[a-1],slst[j-1]);  fpcnt++; }
	}

	if (lst[u-1] == VIN) {
		for (i=1;i<=p;i++)  if (lst[i-1] == VIN && i != u) {
			for (j=1;j<=ul;j++)  {
				if (j < u1 && lst[j-1] != VIN) continue;
				val = m1->v(i,j) +  fl[i-1] * (S2->v(slst[j-1],slst[u-1])-tmpv[j-1]);
				m2->put(val,i,j);
				fpcnt++;
			}
		}
		for (j=1;j<=ul;j++)  {
			if (j < u1 && lst[j-1] != VIN) continue;
			val = fl[u-1] * (S2->v(slst[u-1],slst[j-1])-tmpv[j-1]);
			m2->put(val,u,j);
			fpcnt++;
		}
  }

	if (lst[u-1] == VOUT) {
		for (i=1;i<=p;i++) if (lst[i-1] == VIN){
			for (j=1;j<=p;j++)  {
				if ( (j < u1 || j > ul) && lst[j-1] != VIN) continue;
				val = m1->v(i,j) + wcol[i-1]->v(u)*S2->v(slst[u-1],slst[j-1]) - fl[i-1]*tmpv[j-1];
				m2->put(val,i,j);
				fpcnt += 2;
			}
		}
	}

	return;
}

