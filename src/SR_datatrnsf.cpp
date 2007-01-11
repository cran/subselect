#include <ctime>
#include <cmath>
#include <deque>
#include <string>
#include "SRI_sscma.h"
#include "Sscma.h"
#include "Subsets.h"
#include "Vsmabo.h"
#include "Qforms.h"
#include "VSQforms.h"
#include "GCDcrt.h"
#include "RMcrt.h"
#include "RVcrt.h"
#include "MStcrt.h"
#include "CCRcrt.h"
#include "Rnk3CCRcrt.h"

namespace extendedleaps {

extern double	btime,maxtime;

extern const int  GCD      = 1;
extern const int  RV       = 2;
extern const int  MCB2     = 3;
extern const int  TAU      = 4;	
extern const int  XI       = 5;	
extern const int  ZETA     = 6;	
extern const int  CCR1     = 7;	
extern const int  NOTFOUND = 99;

extern int pcrt;    

bool fwrkspace = false;

using std::deque;

extern subsetdata *idata,*fulldata;
extern globaldata *gidata,*gfulldata;
extern vector<partialdata *> pdata;
extern int *prvks,*cmpl,*ivlst,*ovlst;    
extern vind ndim,maxdim;      
extern real  c0,v0,*vc0;                    

void trnsfdwst(double *S,double *Sinv,double *E,double *Einv,double wstval,int hrank);
void trnsfdtrst(double *M,double *Minv,double *Hegvct,double *HegvctMinv,double trval,
				int hrank);
void trnsfdccr(double *S,double *Sinv,double *E,double *Einv,
			   double *Hegvct,double *HegvctTinv,double *HegvctEinv,
			   double ccr12,double wstval,double bartpival,double lawhotval,int hrank);
void trnsfdgcd(double *S,double *Sinv,double *Segval,double *Segvct,int npcs);
void trnsfdrm(double *S,double *Sinv);
void trnsfdrv(double *S,double *Sinv,double *Ssqr);

void resetvar(void);
int getpcrt(char* st,bool fixed);
void initvlist(int *,int *,int *,int,int,int);
void fillres(vind fk,vind nk,int ns,int* bst,int* st,real* bvl,real* vl);
void asgmemory(void);
void cleanup(void);

int callsscma(double* S,double* S2,double* Si,double* Segval,double* Segvct,
	  double* E,double* Ei,double* Hegvct,double* HegvctTinv,double* HegvctEinv,
	  double wilksval,double bartpival,double lawhotval,double ccr12val,int r,
	  int kmin,int kmax,int nsol,int* out,int* in,int nout,int nin,
	  char* cmpcr,int fixed,int* pcind,int nind,int nvar,double timelimit,
	  int* found,int* subs,double* subsv,double* bestsv,int* bests)
{
	resetvar();
	btime = clock();
	maxtime = timelimit*CLOCKS_PER_SEC;
	p = nvar;
	ms = nsol;
	if (kmin > nin) mindim = kmin;
	else mindim = nin;
	if (kmax < nvar-nout) maxdim = kmax;
	else maxdim = nvar - nout;
	ndim = maxdim-mindim+1;
	if ( (pcrt = getpcrt(cmpcr,static_cast<bool>(fixed))) == NOTFOUND )  {  
		std::string st1("The Comparison criterion suplied, "),st2(cmpcr),st3(", is not supported\n");
		errmsg(st1+st2+st3);
	}
	initvlist(in,out,pcind,nin,nout,nind);
	asgmemory();

	switch (pcrt)  {
		case TAU:    
			trnsfdwst(S,Si,E,Ei,wilksval,r);
			break;
		case XI: 
			trnsfdtrst(S,Si,Hegvct,HegvctTinv,bartpival,r);
			break;
		case ZETA: 
			trnsfdtrst(E,Ei,Hegvct,HegvctEinv,lawhotval,r);
			break;
		case CCR1:    
			trnsfdccr(S,Si,E,Ei,Hegvct,HegvctTinv,HegvctEinv,ccr12val,wilksval,bartpival,lawhotval,r);
			break;
		case GCD: 
			trnsfdgcd(S,Si,Segval,Segvct,q);
			break;
		case MCB2: 
			trnsfdrm(S,Si);
			break;
		case RV: 
			trnsfdrv(S,Si,S2);
			break;
	}
	if ( (*found = sscma(fwrkspace,idata,fulldata)) == true)
		fillres(mindim,ndim,nsol,bests,subs,bestsv,subsv);
	cleanup();
	return 0;
}

void trnsfdwst(double *S,double *Sinv,double *E,double *Einv,double wstval,int hrank)
{
	wilksdata  *idataaswilks=0,*fulldataaswilks=0;
			
	try  {
		pdata.reserve(p+1);
		{ for (int j=0;j<=p;j++) pdata[j] = 0; }
		for (int j=0;j<=p;j++) pdata[j] = new partialwilksdata(p,0.);
		idataaswilks = static_cast<wilksdata *>(idata = new wilksdata(0,p,p,hrank,1.));
		fulldataaswilks = static_cast<wilksdata *>(fulldata = new wilksdata(p,p,p,hrank,wstval));  
	}
	catch (std::bad_alloc)   {
		cleanup();
		errmsg(memmsg);
	}
	{ for (int i=0;i<p;i++)   
		for (int j=0;j<=i;j++) {
			idataaswilks->setematcoef(i,j,E[j*p+i]); 
			idataaswilks->settmatcoef(i,j,S[j*p+i]); 
			fulldataaswilks->setematcoef(i,j,-Einv[j*p+i]); 
			fulldataaswilks->settmatcoef(i,j,-Sinv[j*p+i]); 
		}
	}
}

void trnsfdtrst(double *M,double *Minv,double *Hegvct,double *HegvctMinv,double trval,int hrank)
{
	tracedata  *idataastrst=0,*fulldataastrst=0;
			
	try  {
		pdata.reserve(p+1);
		{ for (int j=0;j<=p;j++) pdata[j] = 0; }
		for (int j=0;j<=p;j++) pdata[j] = new partialtracedata(p,hrank);
		if (pcrt == XI)  {
			idataastrst = static_cast<tracedata *>(idata = new bartpistdata(0,p,p,hrank,0.));
			fulldataastrst= static_cast<tracedata *>(fulldata = new bartpistdata(p,p,p,hrank,v0=trval));  
		}
		else if (pcrt == ZETA)  {
			idataastrst = static_cast<tracedata *>(idata = new lawlhotstdata(0,p,p,hrank,0.));
			fulldataastrst= static_cast<tracedata *>(fulldata = new lawlhotstdata(p,p,p,hrank,v0=trval));  
		}
	}
	catch (std::bad_alloc)   {
		cleanup();
		errmsg(memmsg);
	}
	{ for (int i=0;i<p;i++)   
		for (int j=0;j<=i;j++) {
			idataastrst->getqfdata()->setcoefmatel(i,j,M[j*p+i]); 
			fulldataastrst->getqfdata()->setcoefmatel(i,j,-Minv[j*p+i]); 
		}
	}
	{ for (int i=0;i<hrank;i++)  {
		for (int j=0;j<p;j++) {
			idataastrst->getqfdata()->setvectel(i,j,Hegvct[i*p+j]); 
			fulldataastrst->getqfdata()->setvectel(i,j,-HegvctMinv[i*p+j]); 
		}
	} }
	return;
}


void trnsfdccr(double *S,double *Sinv,double *E,double *Einv,
			   double *Hegvct,double *HegvctTinv,double *HegvctEinv,
			   double ccr12,double wstval,double bartpival,double lawhotval,int hrank)
{
	singleqfdata  *idataassgqf=0,*fulldataassgqf=0;
	ccrdata  *idataasccr=0,*fulldataasccr=0;
	rnk3ccrdata *idataasrnk3ccr=0,*fulldataasrnk3ccr=0;
		
	try  {
		pdata.reserve(p+1);
		{ for (int j=0;j<=p;j++) pdata[j] = 0; }
		if (hrank == 1)  {
			for (int j=0;j<=p;j++) pdata[j] = new partialsingleqfdata();
			idataassgqf = 
				static_cast<singleqfdata *>(idata = new singleqfdata(p,p,0.));
			fulldataassgqf = 
				static_cast<singleqfdata *>(fulldata = new singleqfdata(p,p,ccr12));
		}
		else if (hrank == 2)  {
			for (int j=0;j<=p;j++) pdata[j] = new partialccrdata(0,hrank);
			idataasccr = 
				static_cast<ccrdata *>(idata = new rnk2ccrdata(0,p,p,1.,0.,0.));
			fulldataasccr = 
				static_cast<ccrdata *>(fulldata = new rnk2ccrdata(p,p,p,wstval,bartpival,ccr12));
		}
		else if (hrank == 3)  {
			for (int j=0;j<=p;j++) pdata[j] = new partialrnk3ccrdata(0,hrank);
			idataasrnk3ccr = static_cast<rnk3ccrdata *>( idataasccr = 
					static_cast<ccrdata *>(idata = new rnk3ccrdata(0,p,p,1.,0.,0.,0.)) );
			fulldataasrnk3ccr = static_cast<rnk3ccrdata *>( fulldataasccr = 
				static_cast<ccrdata *>(fulldata = new rnk3ccrdata(p,p,p,wstval,bartpival,lawhotval,ccr12)) );
		}
	}
	catch (std::bad_alloc)   {
		cleanup();
		errmsg(memmsg);
	}
	if (hrank == 1) {
		for (int i=0;i<p;i++)  { 
			for (int j=0;j<=i;j++) {
				idataassgqf->getqfdata()->setcoefmatel(i,j,S[j*p+i]); 
				fulldataassgqf->getqfdata()->setcoefmatel(i,j,-Sinv[j*p+i]); 
			}
		}
		for (int j=0;j<p;j++) {
			idataassgqf->getqfdata()->setvectel(0,j,Hegvct[j]); 
			fulldataassgqf->getqfdata()->setvectel(0,j,-HegvctTinv[j]); 
		}
	}
	else  {
		{ for (int i=0;i<p;i++)   
			for (int j=0;j<=i;j++) {
				idataasccr->settmatcoef(i,j,S[j*p+i]); 
				fulldataasccr->settmatcoef(i,j,-Sinv[j*p+i]); 
				idataasccr->setematcoef(i,j,E[j*p+i]); 
				fulldataasccr->setematcoef(i,j,-Einv[j*p+i]);
			}
		}
		{ for (int i=0;i<hrank;i++)  {
			for (int j=0;j<p;j++) {
				idataasccr->sethtinvel(i,j,Hegvct[i*p+j]); 
				fulldataasccr->sethtinvel(i,j,-HegvctTinv[i*p+j]); 
			}
		} }
		if (hrank == 3) { 
			for (int i=0;i<3;i++)  {
				for (int j=0;j<p;j++) {
					idataasrnk3ccr->setheinvel(i,j,Hegvct[i*p+j]); 
					fulldataasrnk3ccr->setheinvel(i,j,-HegvctEinv[i*p+j]); 
				}
			}
		}
	}
	return;
}

void trnsfdgcd(double *S,double *Sinv,double *Segval,double *Segvct,int npcs)
{
	real srtegval;
	gcddata  *idataasgcd=0,*fulldataasgcd=0;
			
	try  {
		pdata.reserve(p+1);
		for (int j=0;j<=p;j++) pdata[j] = 0;
		switch (pcsets)  {
			case (given): {
					for (int j=0;j<=p;j++) pdata[j] = new partialfgcddata(p,npcs);
					idataasgcd = static_cast<gcddata *>(idata = new fgcddata(0,p,p,npcs,0.));
					fulldataasgcd = static_cast<gcddata *>(fulldata = new fgcddata(p,p,p,npcs,v0=npcs));  
				}
				break;
			case (firstk): {
					{ for (int j=0;j<=p;j++) pdata[j] = new partialvgcddata(p,p); }
					idataasgcd = static_cast<gcddata *>(idata = new vgcddata(0,p,p,0.,0.));
					fulldataasgcd = static_cast<gcddata *>(fulldata = new vgcddata(p,p,p,1.,v0=p));  
					for (int j=0;j<npcs;j++) vc0[j] = 0.;  
				}
				break;
		}
	}
	catch (std::bad_alloc)   {
		cleanup();
		errmsg(memmsg);
	}
	{ for (int i=0;i<p;i++)   
		for (int j=0;j<=i;j++) {
			idataasgcd->getqfdata()->setcoefmatel(i,j,S[j*p+i]); 
			fulldataasgcd->getqfdata()->setcoefmatel(i,j,-Sinv[j*p+i]); 
		}
	}
	{ for (int i=0;i<npcs;i++)  {
		srtegval = sqrt(Segval[i]);
		for (int j=0;j<p;j++) {
			idataasgcd->getqfdata()->setvectel(i,j,srtegval*Segvct[i*p+j]); 
			fulldataasgcd->getqfdata()->setvectel(i,j,-Segvct[i*p+j]/srtegval); 
		}
	} } 
}

void trnsfdrm(double *S,double *Sinv)
{
	deque<bool>	avars(p,false); 
	rmdata *idataasrmdt=0,*fulldataasrmdt=0;

	real trs = S[0];
	{ for (int i=1;i<p;i++) trs += S[i*p+i]; }
	try  {
		pdata.reserve(p+1);
		{ for (int j=0;j<=p;j++) pdata[j] = 0; }
		for (int j=0;j<=p;j++) pdata[j] = new partialrmdata(p);
		rmgdata *gdataasrmdt = static_cast<rmgdata *>(gidata = new rmgdata(p));
		idataasrmdt = static_cast<rmdata *>(idata = new rmdata(p,p,p,gdataasrmdt,avars,trs));
		gdataasrmdt->settrs(trs);
		avars.assign(p,true);
		fulldataasrmdt = static_cast<rmdata *>(fulldata = new rmdata(p,p,p,gdataasrmdt,avars,c0=0.));
	}
	catch (std::bad_alloc)   {
		cleanup();
		errmsg(memmsg);
	}
	for (int i=0;i<p;i++)
		for (int j=0;j<=i;j++) {
			idataasrmdt->setcoefmatel(i,j,S[j*p+i]); 
			fulldataasrmdt->setcoefmatel(i,j,-Sinv[j*p+i]); 
		}
}

void trnsfdrv(double *S,double *Sinv,double *Ssqr)
{
	deque<bool>	avars(p,false); 
	rvdata *idataasrvdt=0,*fulldataasrvdt=0;

	real trs2 = Ssqr[0];
	for (int i=1;i<p;i++) trs2 += Ssqr[i*p+i];
	try  {
		pdata.reserve(p+1);
		{ for (int j=0;j<=p;j++) pdata[j] = 0; }
		for (int j=0;j<=p;j++) pdata[j] = new partialrvdata(p);
		rvgdata *gdataasrvdt = static_cast<rvgdata *>(gidata = new rvgdata(p));
		{ for (int i=0;i<p;i++)
			for (int j=0;j<=i;j++) gdataasrvdt->sets2(i,j,Ssqr[j*p+i]);  }
		idataasrvdt =  static_cast<rvdata *>(idata = new rvdata(p,p,p,gdataasrvdt,avars,0,0.));
		avars.assign(p,true);
		fulldataasrvdt = static_cast<rvdata *>(fulldata = new rvdata(p,p,p,gdataasrvdt,avars,0,c0=trs2));
		gdataasrvdt->settrs2(trs2);		
	}
	catch (std::bad_alloc)   {
		cleanup();
		errmsg(memmsg);
	}
	{ for (int i=0;i<p;i++)
		for (int j=0;j<=i;j++) {
			idataasrvdt->setcoefmatel(i,j,S[j*p+i]); 
			idataasrvdt->sets2m1(i,j,0.); 
			idataasrvdt->sets2m1(j,i,0.); 
			fulldataasrvdt->setcoefmatel(i,j,-Sinv[j*p+i]); 
			fulldataasrvdt->sets2m1(i,j,S[j*p+i]); 
			fulldataasrvdt->sets2m1(j,i,S[j*p+i]); 
		}
	} 
	fwrkspace = true; 
}

}
