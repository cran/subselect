#include <ctime>
#include <cmath>
#include <deque>
#include <cstring> 
#include "SR_sscma.h"
#include "Sscma.h"
#include "Subsets.h"
#include "Vsmabo.h"
#include "Qform.h"
#include "VQform.h"
#include "RMcrt.h"
#include "RVcrt.h"


namespace leapsnbnds  {

double	btime,maxtime; 

}
bool fwrkspace = false;

using namespace leapsnbnds;
using std::deque;

extern subsetdata *idata,*fulldata;
extern globaldata *gidata,*gfulldata;
extern auxmemory  *dataam;
extern vind ndim,maxdim,*prvks,*cmpl,*ivlst,*ovlst;      
extern real  c0,v0,*vc0;                    


void transferdata(real *,real *,real *,real *,real *);
void resetvar(void);
int getpcrt(char *,int *);
void initvlist(int *,int *,int *,int,int,int);
void fillres(vind,vind,long unsigned,int *,int *,real *,real *);
void asgmemory(void);
void cleanup(void);

void transferdata(double *S,double *Ssqr,double *Sinv,double *Segval,double *Segvct)
{

	switch (pcrt)  {
		case GCD:    {
			real    sqrtlmdi;
			qfdata  *idataasqf=0,*fulldataasqf=0;
			
			try  {
				switch (pcsets)  {
					case (given): {
							qfauxmem* dataamasqf = static_cast<qfauxmem *>(dataam = new qfauxmem(q));
							idataasqf = static_cast<qfdata *>(idata = new qfdata(p,p,q,dataamasqf,0.));
							fulldataasqf = static_cast<qfdata *>(fulldata = new qfdata(p,p,q,dataamasqf,v0=q));  
						}
						break;
					case (firstk):  {
							vqfauxmem* dataamasvqf = static_cast<vqfauxmem *>(dataam = new vqfauxmem(q));
							vqfgdata* gidataasqf = static_cast<vqfgdata *>(gidata = new vqfgdata(forward));
							vqfgdata* gfulldataasqf = static_cast<vqfgdata *>(gfulldata = new vqfgdata(backward));
							idataasqf = 
								static_cast<qfdata *>(idata = new vqfdata(p,p,q,dataamasvqf,gidataasqf,0,0.));
							fulldataasqf = 
								static_cast<qfdata *>(fulldata = new vqfdata(p,p,q,dataamasvqf,gfulldataasqf,p,v0=q));  
							for (vind j=0;j<q;j++) vc0[j] = 0.;  
						}
						break;
					}
			}
			catch (std::bad_alloc)   {
				cleanup();
				errmsg(memmsg);
			}

			{ for (vind i=0;i<p;i++)   
				for (vind j=0;j<=i;j++) {
					idataasqf->setcoefmatel(i,j,S[j*p+i]); 
					fulldataasqf->setcoefmatel(i,j,-Sinv[j*p+i]); 
				}
			}
			{ for (vind i=0;i<q;i++)  {
				sqrtlmdi = sqrt(Segval[i]);
				for (vind j=0;j<p;j++) {
					idataasqf->setvectel(i,j,sqrtlmdi*Segvct[i*p+j]); 
					fulldataasqf->setvectel(i,j,-Segvct[i*p+j]/sqrtlmdi);
				}
			} } }
			break;
		case MCB2:  {
				deque<bool>	avars(p,false); 
				rmdata *idataasrmdt=0,*fulldataasrmdt=0;
				trs = S[0];
				{ for (vind i=1;i<p;i++) trs += S[i*p+i]; }
				try  {
					rmgdata *gdataasrmdt = static_cast<rmgdata *>(gidata = new rmgdata(p));
					idataasrmdt = static_cast<rmdata *>(idata = new rmdata(p,p,p,gdataasrmdt,avars,trs));
					avars.assign(p,true);
					fulldataasrmdt = static_cast<rmdata *>(fulldata = new rmdata(p,p,p,gdataasrmdt,avars,c0=0.));
				}
				catch (std::bad_alloc)   {
					cleanup();
					errmsg(memmsg);
				}

				for (vind i=0;i<p;i++)
					for (vind j=0;j<=i;j++) {
						idataasrmdt->setcoefmatel(i,j,S[j*p+i]); 
						fulldataasrmdt->setcoefmatel(i,j,-Sinv[j*p+i]); 
					}
			}
			break;
		case RV:  {
				deque<bool>	avars(p,false); 
				rvdata *idataasrvdt=0,*fulldataasrvdt=0;
				trs2 = Ssqr[0];
				for (vind i=1;i<p;i++) trs2 += Ssqr[i*p+i];
			
				try  {
					rvgdata *gdataasrvdt = static_cast<rvgdata *>(gidata = new rvgdata(p));
					{ for (vind i=0;i<p;i++)
						for (vind j=0;j<=i;j++) gdataasrvdt->sets2(i,j,Ssqr[j*p+i]);
					}
					idataasrvdt =  static_cast<rvdata *>(idata = new rvdata(p,p,p,gdataasrvdt,avars,0,0.));
					avars.assign(p,true);
					fulldataasrvdt = 
						static_cast<rvdata *>(fulldata = new rvdata(p,p,p,gdataasrvdt,avars,0,c0=trs2));
				}
				catch (std::bad_alloc)   {
					cleanup();
					errmsg(memmsg);
				}

				{ for (vind i=0;i<p;i++)
					for (vind j=0;j<=i;j++) {
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
			break;
	}
}


int callsscma(double *S,double *S2,double *Si,double *Segval,double *Segvct,
			  int kmin,int kmax,int nsol,int *out,int *in,int nout,int nin,
			  char *cmpcr,int *fixed,int *pcind,int nind,int nvar,double timelimit,
			  int *found,int *subs,double *subsv,double *bestsv,int *bests)
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
	if ( (pcrt = getpcrt(cmpcr,fixed)) == NOTFOUND )  
		errmsg("The Comparison criterion suplied is not supported\n");
	initvlist(in,out,pcind,nin,nout,nind);
	asgmemory();
	transferdata(S,S2,Si,Segval,Segvct);
	if ( (*found = sscma(fwrkspace,idata,fulldata)) == true)
		fillres(mindim,ndim,nsol,bests,subs,bestsv,subsv);
	cleanup();
	return 0;
}
