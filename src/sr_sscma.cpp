#include <assert.h>
#include <iostream.h>
#include <string.h>
#include <time.h>
#include "fullsqmi.h"
#include "fullmati.h"
#include "SR_vsda.h"
#include "lagmvct.h"
#include "SR_vsmao.h"
#include "vsmac.h"
#include "vsmav.h"

const double      EPS        = 1E-12;
const unsigned    UNKNOWNINT = 99999;
const unsigned    VC         =     0;
const unsigned    CORR       =     1;
const unsigned    CHS        =     2;
const unsigned    MAXVNL     =    20;
const unsigned    NOTFOUND   =    99;

fullvct    **cmpv;
symatrix   *S,*S2;    
wrkspace   *SW,*IW;   
unsigned p,fp,lp,q,ms,m0,pcrt,pcrttp,maxsbst,maxsbqe,flst,mttp,maxcmp,mindim,maxdim,ndim; 
unsigned *prvks,*dmyv,*actv,*cmpl,*Flp,*ivlst,*ovlst;   
double   *lbnd,*ubnd,*bndl,*Fl,c0,v0,trs,trs2,maxtime;  
long unsigned cntg,cntp,fpcnt,fpcnt0,fpcnt1;
extern unsigned flsts,flsti,lastvar;
extern int sbsetind,sbqeind;
double     btime;
psbsq      *bsts;              
psbst      *sbsarr;
psbstqe    *sbqearr;

unsigned getpcrt(char *);
void initvlist(int *,int *,int *,int,int,int);
void compS(double *,unsigned);
int sscma(void);
void fillres(unsigned,unsigned,unsigned,int *,int *,double *,double *);

char cmpcrMC2[] =  "McCabe's 2nd Criterion:\nTr[S(out.in)]";
char cmpcrGCD[] =  "Tr [(S_in)^-1 SG( )_in]";
char cmpcrGCD1[]=  "Tr [(S_in)^-1 SG(  )_in]";
char cmpcrRV[] =   "Tr {[(S_in)^-1 (S^2)_in]^2}";
char indGCD[]  =   "GCD";
char indRM[]   =   "RM";
char indRV[]   =   "RV";

void matasvcttranspose(unsigned,unsigned,int *);
void initpvar(void);

int callsscma(double *mat,int kmin,int kmax,int nsol,
			  int *out,int *in,int nout,int nin,
			  char *cmpcr,int *pcind,int nind,int nvar,double timelimit,
			  int *found,int *subs,double *subsv,double *bestsv,int *bests)
{
	btime = clock();
	initpvar();
	p = nvar;
	ms = nsol;
	if (kmin > nin) mindim = kmin;
	else mindim = nin;
	if (kmax < nvar-nout) maxdim = kmax;
	else maxdim = nvar - nout;
	ndim = maxdim-mindim+1;
	if ( (pcrt = getpcrt(cmpcr)) == NOTFOUND )  
		errmsg("The Comparison criterion supplied is not supported\n");
	initvlist(in,out,pcind,nin,nout,nind);
	maxtime = timelimit*CLOCKS_PER_SEC;
	asgmemory();
	compS(mat,nvar);
	if ( (*found = sscma()) == TRUE)
		fillres(mindim,ndim,nsol,bests,subs,bestsv,subsv);
	cleanup();
	return 0;
}

void initpvar()
{
	cmpv = NULL;
    	S = S2 = NULL;
	SW = IW = NULL;
	prvks = dmyv = actv = ivlst = ovlst = cmpl = Flp = NULL;  
	lbnd = ubnd = bndl = Fl = NULL;  
    	bsts = NULL;    
    	sbqearr = NULL;  
    	sbsarr = NULL;
	mttp = VC;
	cntg = cntp = 0;
	sbsetind = sbqeind = 0;
}

unsigned getpcrt(char *st)
{
	if ( !strncmp(st,"rm",2) || !strncmp(st,"RM",2) || !strncmp(st,"Rm",2) )        return MCB2;    
	else if ( !strncmp(st,"rv",2) || !strncmp(st,"RV",2) || !strncmp(st,"Rv",2) )   return RV;
	else if ( !strncmp(st,"gcd",3) || !strncmp(st,"GCD",3) || !strncmp(st,"Gcd",3) ) return GCD;  
	else return NOTFOUND;
}

void initvlist(int *ilist,int *olist,int *pcind,int ni,int no,int nind)
{
	int i,j,vn0;  

	if ( ni > 0  && (ivlst = new unsigned[ni]) ==  NULL)  prmtend("ivlst");
	if ( no > 0  && (ovlst = new unsigned[no]) ==  NULL)  prmtend("ovlst");

	fp = ni;
	for (i=0;i<ni;i++) ivlst[i] = ilist[i];
	lp = no;
	for (i=0;i<no;i++) ovlst[i] = olist[i];

	if (pcrt != GCD) q = 0;
	else {
		if ( (q=nind) == 0) 
			errmsg("Criterion GCD requires a non-empty list of S eigenvectors\n");
		if ( (cmpl =  new unsigned[q]) == NULL) prmtend("cmpl");
		if ( (cmpv = new (fullvct *[q])) == NULL) prmtend("cmpv");
		for (i=1;i<=(int)q;i++)  {
			if ( (cmpv[i-1] = new fullvct(p)) == NULL) prmtend("fullvct");
			if (i == 1) maxcmp = cmpl[0] = pcind[0];
			else if ( (cmpl[i-1] = pcind[i-1]) > maxcmp )
				maxcmp = cmpl[i-1];
		}
	}
}

void compS(double *data,unsigned ncol)
{
	unsigned i,j;

	for (i=1;i<=p;i++)
		for (j=1;j<=i;j++)  S->put(data[(i-1)*ncol+j-1],i,j);
}


int sscma()
{
	crtwrksp();
	flst = flsts;
	fpcnt1 = 0;
	if (p > fp+lp+1) fsort();
	else lastvar = IW->wrklst[flsti]->srtv[p-1];
	assert(lastvar > 0 && lastvar <= p);
	if (fp > 0 && fp == mindim) savfrst();
	if (maxdim == p-lp) savfull();
	fpcnt = 0;
	if (p <= fp+lp+1 || prcksp(flst,flst,fp,p-lp,fp+lp+1) == TRUE)  return TRUE;
	else return FALSE;
}

void fillres(unsigned fk,unsigned nk,unsigned ns,int *bst,int *st,double *bvl,double *vl)
{
	unsigned i,lk=fk+nk-1;
	int *stmatp;
	
	for (i=0;i<nk;i++)  {
		stmatp = &(st[i*ns*lk]);
		bsts[i]->save(stmatp,&(vl[i*ns]),ns,lk,-1);
		matasvcttranspose(ns,lk,stmatp);
		bsts[i]->save(&(bst[i*lk]),&(bvl[i]),1,lk,-1);
	}
	matasvcttranspose(nk,lk,bst);
	return;
}

void matasvcttranspose(unsigned m,unsigned n,int *data)
{
	unsigned i,j,mn=m*n;
	int  *tmp;
	
	tmp = new int[mn];
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)  tmp[i+j*m] = data[i*n+j];
	for (i=0;i<mn;i++) data[i] = tmp[i];
	delete[] tmp;
}

void asgmemory()
{
	unsigned i;
	extern double *cv,*tv;
	extern fwrksqm  *m1t;

	if ( (S = new symatrix(p)) == NULL) prmtend("S");
	if ( pcrt == RV && (S2 = new symatrix(p)) == NULL ) prmtend("S2");
	if ( (actv = new unsigned[p]) == NULL) prmtend("actv");
	if ( (dmyv = new unsigned[p]) == NULL) prmtend("dmyv");
	if ( (Fl = new double[p]) == NULL) prmtend("Fl");
	if ( (Flp = new unsigned[p]) == NULL) prmtend("Flp");
	if (ms > 0 && (bsts = new psbsq[ndim]) == NULL) prmtend("bsts");
	if (ndim == p-fp-lp) maxsbst = maxsbqe = ms*(ndim-1)+3;
	else maxsbst = maxsbqe = ms*ndim+2;
	if ( (sbsarr = new psbst[maxsbst]) == NULL) prmtend("psbst");
	if ( (sbqearr = new psbstqe[maxsbqe]) == NULL) prmtend("sbqearr");
	for (i=0;i<maxsbst;i++) if ((sbsarr[i]=new sbset(i,p)) == NULL) prmtend("sbsarr-i");
	for (i=0;i<maxsbqe;i++) if ((sbqearr[i]=new sbsetqe(i)) == NULL) prmtend("sbqearr-i");
	if ( pcrt == MCB2 ) {
		pcrttp = MIN;
		if (ms > 0) 
			if (ndim == p-fp-lp) { 
				for (i=0;i<ndim-1;i++) if ( (bsts[i] = new subsetq(-1,ms)) == NULL) prmtend("bsts-i");
				if ( (bsts[ndim-1] = new subsetq(-1,1)) == NULL) prmtend("bsts-full");
				if ( (ubnd = new double[ndim-1]) == NULL) prmtend("ubnd");
				for (i=0;i<ndim-1;i++) ubnd[i] = BIG;
			}
			else {
				for (i=0;i<ndim;i++) if ( (bsts[i] = new subsetq(-1,ms)) == NULL) prmtend("bsts-i");
				if ( (ubnd = new double[ndim]) == NULL) prmtend("ubnd");
				for (i=0;i<ndim;i++) ubnd[i] = BIG;
			}
		lbnd = NULL;
	}
	else  {
		pcrttp = MAX;
		if (ms > 0) {
			if (ndim == p-fp-lp) { 
				for (i=0;i<ndim-1;i++) if ( (bsts[i] = new subsetq(1,ms)) == NULL) prmtend("bsts-i");
				if ( (bsts[ndim-1] = new subsetq(1,1)) == NULL) prmtend("bsts-full");
				if ( (lbnd = new double[ndim-1]) == NULL) prmtend("lbnd");
				for (i=0;i<ndim-1;i++) lbnd[i] = 0.;
			}
			else {
				for (i=0;i<ndim;i++) if ( (bsts[i] = new subsetq(1,ms)) == NULL) prmtend("bsts-i");
				if ( (lbnd = new double[ndim]) == NULL) prmtend("lbnd");
				for (i=0;i<ndim;i++) lbnd[i] = 0.;
			}
		}
		ubnd = NULL;
	}
	if (pcrt == RV) {
		if ( (m1t = new fwrksqm(p)) == NULL) prmtend("m1t");
		if ( (cv = new double[p]) == NULL) prmtend("cv");
		if ( (tv = new double[p]) == NULL) prmtend("tv");
	}
	else {
		m1t = NULL;
		cv = tv = NULL;
		bndl = NULL;    
	}
	if ( (prvks=new unsigned[p-1]) == NULL) prmtend("unsigned");
	return;
}

void cleanup(void)
{
	unsigned i;
	
	extern double *cv,*tv;
	extern fwrksqm  *m1t;

	delete   SW;
	delete   IW;
	delete S;
	delete S2;
	delete[] actv;
	delete[] dmyv;
	delete[] Fl;
	delete[] Flp;
	if (bsts != NULL) {
		for (i=0;i<ndim;i++) delete bsts[i];
		delete[] bsts;
	}
	delete[] lbnd;
	delete[] ubnd;
	delete[] bndl;
	delete[] prvks;
	if ( sbsarr != NULL) {
		for (i=0;i< maxsbst;i++) delete sbsarr[i];
		delete[] sbsarr;
	}
	if ( sbqearr != NULL) {
		for (i=0;i< maxsbqe;i++) delete sbqearr[i];
		delete[] sbqearr;
	}
	delete m1t;
	delete[] cv;
	delete[] tv;
	delete[] ivlst;
	delete[] ovlst;
	delete[] cmpl;
	if (cmpv !=  NULL)  {
		for (i=0;i<q;i++) delete cmpv[i];
		delete[] cmpv;
	}
	return;
}





