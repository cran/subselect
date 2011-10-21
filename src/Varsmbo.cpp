#include <cstdlib>
#include <cassert>
#include <vector>
#include <limits>
#include "Sscma.h"
#include "Vsmabo.h"

using namespace std;

namespace extendedleaps {

const double INF = std::numeric_limits<double>::infinity(); 
vind *dmyv,*Flp;  
double *Fl;	 

extern std::vector<partialdata *> pdata;
extern vind flsts,flsti;
extern bool numericalprob;

int  cmp(const void *,const void *);
int  revcmp(const void *,const void *);


SRCwrkspace::SRCwrkspace(bool pivotall,vind nv,subsetdata *data0,vind* ivlst,vind* ovlst)
{
	pvall_ = pivotall;
	if (pivotall) flsts = nv-lp;
	else flsts = nv-lp-1;
	initwrkspace(pivotall,nv,data0,flsts,fp,lp,ivlst,ovlst);
	flsts -= fp;
}

bool SRCwrkspace::pivot(vind vp,vind t,vind li,vind lo,double acpbound)
{
	subset *isi=0,*iso=0;

	pivotinit(isi,iso,vp,li,lo);
	return isi->pivot(forward,vp,t,iso,(lo==0),acpbound);
}

INVwrkspace::INVwrkspace(bool pivotall,vind nv,subsetdata *data0,vind* ivlst,vind* ovlst)
{
	pvall_ = pivotall;
	if (pivotall) flsti = nv-fp;
	else flsti = nv-fp-1;
	initwrkspace(pivotall,nv,data0,flsti,lp,fp,ovlst,ivlst);
	flsti -= lp;
}

bool INVwrkspace::pivot(vind vp,vind t,vind li,vind lo,double acpbound)
{
	subset *isi=0,*iso=0;

	pivotinit(isi,iso,vp,li,lo);
	return isi->pivot(backward,vp,t,iso,(lo==0),acpbound);
}

void wrkspace::initwrkspace(bool pivotall,vind nv,subsetdata *data0,vind lstind,vind nvattop,vind nvatbot,vind* vattop,vind* vatbot)
{
	vind lastv;
	vind* tlst=0;
	subset *ispc;
	subsetdata *newdata;
	double trivbnd(-INF);

	maxim = data0->max();
	if ( !max() ) trivbnd *= -1;
	wrklst = new pkspc[lstind+1];
	p = nv;
	if (pivotall)  {
		lastv = nv;
		nwl = p-fp-lp+1;
	}
	else  {
		lastv = nv-1;
		nwl = p-fp-lp;
	}
	if (fp+lp > 0)  {
		tlst   = new vind[p];
		frontlsts(vatbot,vattop,nvatbot,nvattop,tlst);
		ispc = wrklst[lstind] = new subset(tlst,p,p,data0,false,p);
	}
	else ispc = wrklst[lstind] = new subset(p,p,data0,false,p);
	if (lp+fp > 0) ispc->reorder(tlst);
	for (vind j=1;j<=nvattop;j++)  {
		newdata = data0->crcopy(p,p-nvatbot-j);
		if (!tlst) wrklst[lstind-j] = new subset(p,p-nvatbot-j,newdata,true,p);
		else wrklst[lstind-j] = new subset(tlst,p,p-nvatbot-j,newdata,true,p);
		if (lstind > j) pivot(nvatbot+j,p-nvatbot-j,lstind+1-j,lstind-j,trivbnd);
		else pivot(nvatbot+j,0,lstind+1-j,0,trivbnd);
		delete wrklst[lstind+1-j];
	}
	for (int j1=nwl-2;j1>=0;j1--) {
		newdata = data0->crcopy(lastv,j1); 
		if (!tlst) wrklst[j1] = new subset(lastv,j1,newdata,true,p);
		else wrklst[j1] = new subset(tlst,lastv,j1,newdata,true,p);
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

void wrkspace::pivotinit(subset*& isi,subset*& iso,vind vp,vind li,vind lo)
{
	vind isonvar;
	isi = wrklst[li];
	iso = wrklst[lo];

	isi->copyvar(*iso);
	iso->setnvar(isonvar=isi->getnvar()+1);
	iso->setvar(isonvar,vp);
}

subset::subset(vind nvp,vind pnv,subsetdata *id,bool pdt,vind tnv)
   :  p(tnv), t(pnv), k(0), var(0), frstvarpm(nvp-pnv),  pmemorypos(0), memii(0), data(id), privatedata(pdt), nxtres(0) 
{
	assgnmem();
	for (vind i=0;i<p;i++)  
		fmemorypos[i] = orgvarind[i] =  orgvarpos[i] = i;
	if (id) id->setorgvarl(orgvarind);                    // Attach list of original variable indicators to subsetdata
}

subset::subset(vind * const ivar,vind nvp,vind pnv,subsetdata *id,bool pdt,vind tnv)
  :  p(tnv), t(pnv), k(0), var(0), frstvarpm(nvp-pnv),  pmemorypos(0), memii(0), data(id), privatedata(pdt), nxtres(0) 
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
	memii = new mindices(p,p-frstvarpm,frstvarpm,fmemorypos);
}

subset::~subset()
{
	delete[] var;
	delete[] orgvarind;
	delete[] orgvarpos;
	delete[] fmemorypos;
	delete[] pmemorypos;
	if (privatedata) delete data;
	delete memii;
	delete[] nxtres;
}

void subset::copyvar(subset & newsp)
{
	for (vind i=0;i<k;i++)  newsp.var[i] = var[i];  
}

void subset::sort(direction dir,bool reverse,vind fv,vind lv)
{
	bool reliable(true);
	double trivbnd(-INF);

	assert(fv > 0 && lv > fv && lv <= p);
	if (!data->max())  trivbnd *= -1;
	if (!nxtres) nxtres = new knownres[t];
	for (vind i=0;i<=lv-fv;i++)  {
		Fl[i] = data->updatecrt(dir,*memii,fv+i,pdata[i+1],reliable,numtol,trivbnd);
		if (!reliable)	Fl[i] = trivbnd; 
 		Flp[orgvarind[fv+i-1]] = i+1;
		dmyv[i] = i+1;
		nxtres[i].criterion = Fl[i];
		nxtres[i].pres = pdata[i+1];
	}
	if (reverse) qsort((void *)dmyv,lv-fv+1,sizeof(*dmyv),revcmp);
	else {
		qsort((void *)dmyv,lv-fv+1,sizeof(*dmyv),cmp);
		int tmp=dmyv[0];
		for (int a=0;a<lv-fv;a++) dmyv[a] = dmyv[a+1];
		dmyv[lv-fv] = tmp;
	}

	{ for (vind i=fv;i<=lv;i++) dmyv[i-fv] = orgvarind[dmyv[i-fv]+fv-2]; }
	{ for (vind i=fv;i<=lv;i++) orgvarind[i-1] = dmyv[i-fv]; }
}


void subset::reorder(vind *l)
{
	vind lag = p - t;

	if (!pmemorypos) 
		memii->asgnpmmiid(new lagindex<i>(t,frstvarpm,pmemorypos = new vind[t]));
	for (vind i=0;i<p;i++)  {
		fmemorypos[i] = l[i] - 1;
		if (i >= lag) pmemorypos[i-lag] = l[i]-lag-1;
	} 
}

void subset::asgvar(vind fvar,vind nv,vind *list)
{
	vind lag = p - t;

	if (!pmemorypos) 
		memii->asgnpmmiid(new lagindex<i>(t,frstvarpm,pmemorypos = new vind[t]));
	for (vind i=0;i<nv;i++)  {
		pmemorypos[fvar+i] = list[i] -1;
		fmemorypos[lag+fvar+i] = lag + list[i] -1;
	}
}

bool subset::pivot(direction d,vind vp,vind t,subset *newsp,bool last,double acpbound)
{
	bool reliable=true;
	partialdata* pdt;
	vind vporgi;

	if ( nopivot() )  { 
		newsp->forbidpivot(); 	
		numericalprob = true;
		return false;
	}
	if (nxtres)  { 
		if (memii->direct()) vporgi = (*(memii->idpm()))[vp-1];  
		else vporgi = (*(memii->iipm()))[vp-1];  
		pdt = nxtres[vporgi].pres;
	}
	else data->updatecrt(d,*memii,vp,pdt=pdata[0],reliable,numtol,acpbound);

	if (!reliable) {
		newsp->forbidpivot(); 	
		numericalprob = true;
		return false;
	}
	newsp->getdatap()->getpdata(pdt);
	if (!last) data->pivot(d,*memii,vp,t,pdt,newsp->data,last,reliable,INF); 
	
	return true;
}

mindices::mindices(vind szfm,vind szpm,vind pmemlag)
{
	idfm_ = new itindex<d>(szfm);
	idpm_ = new lagindex<d>(szpm,pmemlag);
	iifm_ = 0;
	iipm_ = 0;
}

mindices::mindices(vind szfm,vind szpm,vind pmemlag,vind* fmmlst)
{
	idfm_ = new itindex<d>(szfm);
	idpm_ = new lagindex<d>(szpm,pmemlag);
	iifm_ = new itindex<i>(szfm,fmmlst);
	iipm_ = 0;
}

mindices::mindices(vind szfm,vind szpm,vind pmemlag,vind* fmmlst,vind* pmmlst)
{
	idfm_ = new itindex<d>(szfm);
	idpm_ = new lagindex<d>(szpm,pmemlag);
	iifm_ = new itindex<i>(szfm,fmmlst);
	iipm_ = new lagindex<i>(szpm,pmemlag,pmmlst);
}

mindices::~mindices(void)
{
	delete idfm_;
	delete idpm_;
	delete iifm_;
	delete iipm_;
}

}
