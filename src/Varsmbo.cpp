#include <cstdlib>
#include <cassert>
#include "Sscma.h"
#include "Vsmabo.h"

using namespace leapsnbnds;

vind flsts,flsti;
int  cmp(const void *,const void *);

subset::subset(vind tnv,vind pnv,subsetdata *id,bool pdt,vind nvo)
  :  lstvarpm(tnv-1), frstvarpm(tnv-pnv), data(id), privatedata(pdt), p(nvo), k(0), var(0), pmemorypos(0)
{
	assgnmem();
	for (vind i=0;i<p;i++)  
		fmemorypos[i] = orgvarind[i] =  orgvarpos[i] = i; 
	if (id) id->setorgvarl(orgvarind);                    // Attach list of original variable indicators to subsetdata
}

subset::subset(vind * const ivar,vind tnv,vind pnv,subsetdata *id,bool pdt,vind nvo)
  :  lstvarpm(tnv-1), frstvarpm(tnv-pnv), data(id), privatedata(pdt), p(nvo), k(0), var(0), pmemorypos(0)
{
	assgnmem();
	for (vind i=0;i<p;i++)  {
		orgvarind[i] = ivar[i]-1;
		orgvarpos[orgvarind[i]] = i;
		fmemorypos[i] = i;
	}
	if (id) id->setorgvarl(orgvarind);          // Attach list of original variable indicators to subsetdata
}

void subset::assgnmem()
{
	if (frstvarpm) var = new vind[frstvarpm];
	orgvarind = new vind[p];
	orgvarpos = new vind[p];
	fmemorypos = new vind[p];
}

subset::~subset()
{
	delete[] var;
	delete[] orgvarind;
	delete[] orgvarpos;
	delete[] fmemorypos;
	delete[] pmemorypos;
	if (privatedata) delete data;
}

void subset::copyvar(subset & newsp)
{
	for (vind i=0;i<k;i++)  newsp.var[i] = var[i];  
}

void subset::sort(vind fv,vind lv)
{
	assert(fv > 0 && lv > fv && lv <= p);

	for (vind i=0;i<=lv-fv;i++)  {
		Fl[i] = data->updatecrt(pmemorypos,fmemorypos,fv+i,frstvarpm);
		Flp[orgvarind[fv+i-1]] = i+1;
		dmyv[i] = i+1;
	}
	qsort((void *)dmyv,lv-fv+1,sizeof(*dmyv),cmp);
	{ for (vind i=fv;i<=lv;i++) dmyv[i-fv] = orgvarind[dmyv[i-fv]+fv-2]; }
	{ for (vind i=fv;i<=lv;i++) orgvarind[i-1] = dmyv[i-fv]; }
}

void subset::reorder(vind *list)
{
	vind lag = frstvarpm + (p-lstvarpm) - 1;

	if (!pmemorypos) pmemorypos = new vind[lstvarpm-frstvarpm+1];
	for (vind i=0;i<p;i++)  {
		fmemorypos[i] = list[i] - 1;
		if (i >= lag) pmemorypos[i-lag] = list[i] - lag -1;
	} 
}

void subset::asgvar(vind fvar,vind nv,vind *list)
{
	vind lag = frstvarpm + (p-lstvarpm) - 1;

	if (!pmemorypos) pmemorypos = new vind[lstvarpm-frstvarpm+1];
	for (vind i=0;i<nv;i++)  {
		pmemorypos[fvar+i] = list[i] -1;
		fmemorypos[lag+fvar+i] = lag + list[i] -1;
	}
}

void subset::pivot(vind vp,vind v1,vind vl,subset *newsp,bool last)
{  
	real newcrt = data->updatecrt(pmemorypos,fmemorypos,vp,frstvarpm);
	data->pivot(vp,v1,vl,frstvarpm,newcrt,pmemorypos,fmemorypos,newsp->data,last); 
	for (vind i=0;i<lstvarpm;i++) newsp->fmemorypos[i] = i;
}


wrkspace::wrkspace(bool pvar,vind tp,vind nv,subsetdata *data0)
	:  full(pvar), p(nv)
{
	vind j,lstind=0,nvattop=0,nvatbot=0;       
	vind *vattop=NULL,*vatbot=NULL,*tlst=NULL;
	int j1;
	subset *ispc;
	subsetdata *newdata;
	extern vind *ivlst,*ovlst;

	switch (tp)  {
		case SRC:
			lstind = flsts = p-lp-1;
			nvattop = fp;
			nvatbot = lp;
			vattop = ivlst;
			vatbot = ovlst;
			break;
		case INV:
			lstind = flsti = p-fp-1;
			nvattop = lp;
			nvatbot = fp;
			vattop = ovlst;
			vatbot = ivlst;
			break;
	}
	wrklst = new pkspc[lstind+1];
	if (fp+lp > 0)  {
		tlst   = new vind[p];
		frontlsts(vatbot,vattop,nvatbot,nvattop,tlst);
		ispc = wrklst[lstind] = new subset(tlst,p,p,data0,false,p);
	}
	else ispc = wrklst[lstind] = new subset(p,p,data0,false,p);
	if (lp+fp > 0) ispc->reorder(tlst);
	for (j=1;j<=nvattop;j++)  {
		newdata = data0->crcopy(p,p-nvatbot-j);
		if (!tlst) wrklst[lstind-j] = new subset(p,p-nvatbot-j,newdata,true,p);
		else wrklst[lstind-j] = new subset(tlst,p,p-nvatbot-j,newdata,true,p);
		if (lstind > j) pivot(nvatbot+j,nvatbot+j+1,p,lstind+1-j,lstind-j);
		else pivot(nvatbot+j,nvatbot+j+1,nvatbot+j,lstind+1-j,0);
		delete wrklst[lstind+1-j];
	}
	if (tp == SRC) nwl = (flsts -= fp) + 1;
	if (tp == INV) nwl = (flsti -= lp) + 1;
	for (j1=nwl-2;j1>=0;j1--) {
		newdata = data0->crcopy(p-1,j1);
		if (!tlst) wrklst[j1] = new subset(p-1,j1,newdata,true,p);
		else wrklst[j1] = new subset(tlst,p-1,j1,newdata,true,p);
	}
	delete[] tlst;
}

void wrkspace::frontlsts(vind *l1,vind *l2,vind nl1,vind nl2,vind *ol)
{
	vind tmp,*elep;
	
	elep = new vind[p];
	{ for (vind i=0;i<p;i++)  ol[i] = elep[i] = i+1; }
	{ for (vind i=0;i<nl1;i++) {
		tmp = ol[i];
		ol[i] = l1[i];
		ol[elep[l1[i]-1]-1] = tmp;
		elep[tmp-1] = elep[l1[i]-1];
		elep[l1[i]-1] = i+1;
	} }
	{ for (vind i=0;i<nl2;i++) {
		tmp = ol[nl1+i];
		ol[nl1+i] = l2[i];
		ol[elep[l2[i]-1]-1] = tmp;
		elep[tmp-1] = elep[l2[i]-1];
		elep[l2[i]-1] = nl1+i+1;
	} }
	delete[] elep;
}

wrkspace::~wrkspace(void)
{
	for (vind j=0;j<nwl;j++)  delete wrklst[j];
	delete[] wrklst;
}

void wrkspace::pivot(vind vp,vind v1,vind vl,vind li,vind lo)
{
	vind isonvar;
	subset* isi = wrklst[li];
	subset* iso = wrklst[lo];

	isi->copyvar(*iso);
	iso->setnvar(isonvar=isi->getnvar()+1);
	iso->setvar(isonvar,vp);
	isi->pivot(vp,v1,vl,iso,(lo==0));
}
