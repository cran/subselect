// #include <strstrea.h>
#include <stdlib.h>
#include <math.h>
#include "fullsqmi.h"
#include "fullmati.h"
#include "SR_vsda.h"
#include "vsmac.h"

const double UNABLE    = 1E99;

int trivialcmp(const void *,const void *);


int sbsetind,sbqeind;
extern unsigned maxsbst,maxsbqe;
extern sbset    **sbsarr;
extern sbsetqe  **sbqearr;

sbsetqe *csbsetqe(sbset *s)
{

	sbsetqe  *sqe;
	if ( ++sbqeind >= (int)maxsbqe)  errmsg("Memory management problems");
	sqe = sbqearr[sbqeind-1];
	sqe->assgnset(s);
	sqe->assgnpset(NULL);
	sqe->assgnnset(NULL);
	return sqe;
}

void dsbsetqe(sbsetqe *sqe)
{
	if ( --sbqeind < 0)  errmsg("Memory management problems");
	dsbset(sqe->set());
	(sbqearr[sqe->pos()]=sbqearr[sbqeind])->assgnpos(sqe->pos());
	(sbqearr[sbqeind]=sqe)->assgnpos(sbqeind);
	return;
}

sbset *csbset(unsigned n,unsigned hrnk,unsigned *v,char *cn,char *in,double c)
{
	unsigned i;
	sbset    *s;

	if ( ++sbsetind > (int)maxsbst)  errmsg("Memory management problems");
	s = sbsarr[sbsetind-1];
	s->newset = TRUE;
	s->nvar_ = n;   s->r_ = hrnk;
	for (i=0;i<n;i++) s->actvar_[i] = v[i];
	s->crtn_ = cn;
	s->indn_ = in;
	s->crt_ = c;
	return s;
}

void dsbset(sbset *s)
{
	if ( --sbsetind < 0)  errmsg("Memory management problems");
	(sbsarr[s->pos]=sbsarr[sbsetind])->pos = s->pos;
	(sbsarr[sbsetind]=s)->pos = sbsetind;
	return;
}

sbset::sbset(unsigned p,unsigned n)
:  pos(p), nvar_(n)
{
	actvar_ = new unsigned[nvar_];
}

sbset::~sbset(void)
{
	delete[] actvar_;
}

/*
void sbset::show(const char *st)
{
	unsigned i;
	double ind;

	cout << st << "\nVARS:  ";
	for (i=0;i<nvar_;i++) {
		if (!(i%VARPL) && i) cout << "\n       ";
		cout << "  " << varn[actvar_[i]-1];
	}
	switch (pcrt)  {
		case MCB2:
			ind = sqrt(1. - crt_/trs);
			cout << "\n" << indn_ << " = " << ind << "\n";
			cout << "(" << crtn_ << " = " << crt_ << ")\n";
			break;
		case GCD:
			ind = crt_/sqrt(q*nvar_);
			cout << "\n" << indn_ << " = " << ind << "\n";
			cout << "(" << crtn_ << " = " << crt_ << ")\n";
			break;
		case RV:
			ind = sqrt(crt_/trs2);
			cout << "\n" << indn_ << " = " << ind << "\n";
			cout << "(" << crtn_ << " = " << crt_ << ")\n";
			break;
	}
}
*/

void sbset::save(int *var,double *crtval,unsigned dim)
{
	unsigned i;

	for (i=0;i<nvar_;i++) var[i] = actvar_[i];
	qsort(static_cast<void *>(var),nvar_,sizeof(*var),trivialcmp);
	for (i=nvar_;i<dim;i++) var[i] = 0;
	switch (pcrt)  {
		case MCB2:
			*crtval = sqrt(1. - crt_/trs);
			break;
		case GCD:
			*crtval = crt_/sqrt(q*nvar_);
			break;
		case RV:
			*crtval = sqrt(crt_/trs2);
			break;
	}
}

int trivialcmp(const void *a,const void *b)
{
	unsigned ai,bi;

	ai = *(unsigned *)a;
	bi = *(unsigned *)b;
	if (ai > bi)  return  1;
	else if (ai < bi)  return -1;
		 else return 0;
}


sbsetqe::sbsetqe(unsigned p)    
 :  pos_(p)  
{ }

sbsetqe::sbsetqe(unsigned p,sbset *sbs)  	
 :  pos_(p), set_(sbs), nxtset_(NULL), prvset_(NULL), elegl_(FALSE)  
{ }

char cmpset(sbset *e1,sbset *e2)
{
	if ( e1->crt() > e2->crt() )  return  1;
	if ( e1->crt() < e2->crt() )  return -1;
	else return 0;
}

subsetq::subsetq(int d,unsigned m)
 :  dir_(d), maxele_(m), nele_(0), all_(TRUE)
{
	pmin_ = pmax_ = NULL;
	if (m > 0) full_ = FALSE;
	else full_ = TRUE;
}

void subsetq::delalle(void)
{
	sbsetqe *fele,*elep;

	if (dir_ == 1)  fele = pmin_;
	else fele = pmax_;
	if (fele != NULL)  {
		if ( (elep=fele->nxtset()) != NULL )
			for (;elep!=NULL;elep=elep->nxtset())  
				;
		if (dir_ == 1)  delel(pmax_);
		if (dir_ == -1) delel(pmin_);
	}
}

double subsetq::inserte(sbsetqe *newel)
{
	sbsetqe **felep,*elep;

	if (full_)  {
		if (all_) all_ = FALSE;
		if  ( (dir_ ==  1 && cmpset(newel->set(),pmin_->set()) != 1 )  ||
			  (dir_ == -1 && cmpset(newel->set(),pmax_->set()) != -1)  ) return UNABLE;
	}
	if (dir_ == 1)  felep = &pmin_;
	else felep = &pmax_;
	if ( *felep == NULL ||
		 (dir_ ==  1 && cmpset(newel->set(),pmin_->set()) != 1 )  ||
		 (dir_ == -1 && cmpset(newel->set(),pmax_->set()) != -1)
		)
	 {
		newel->assgnnset(*felep);
		if (*felep != NULL) (*felep)->assgnpset(newel);
		*felep = newel;
	}
	else  {
		elep = *felep;
		while ( elep->nxtset() != NULL &&
				  (   (dir_ == 1  && cmpset(newel->set(),elep->nxtset()->set())  == 1 ) ||
					  (dir_ == -1 && cmpset(newel->set(),elep->nxtset()->set())  == -1)   )
				 ) elep = elep->nxtset();
		newel->assgnnset(elep->nxtset());
		newel->assgnpset(elep);
		if (elep->nxtset() != NULL) elep->nxtset()->assgnpset(newel);
		elep->assgnnset(newel);
		if (full_) {
			*felep = (*felep)->nxtset();
			dsbsetqe((*felep)->prvset());
			(*felep)->assgnpset(NULL);
		}
	}
	if (!full_ && ++nele_ == maxele_) full_ = TRUE;
	if (newel->nxtset() == NULL)  {
		if (dir_ == 1)  pmax_ = newel;
		if (dir_ == -1) pmin_ = newel;
	}

	return (*felep)->set()->crt();
}

void subsetq::delel(sbsetqe *e)
{
	if (e == pmin_ && dir_ == 1) pmin_ = e->nxtset();
	else if (e == pmax_ && dir_ == -1) pmax_ = e->nxtset();
		  else e->prvset()->assgnnset(e->nxtset());
	if (e == pmax_ && dir_ == 1) pmax_ = e->prvset();
	else if (e == pmin_ && dir_ == -1) pmin_ = e->prvset();
		  else e->nxtset()->assgnpset(e->prvset());
	dsbsetqe(e);
	nele_--;
	full_ = FALSE;
}

/*
void subsetq::show(const char *st,unsigned nel,int d1)
{
	unsigned i;
	sbsetqe *fele,*lele,*qep;

	cout << st << "\n";
	if (d1 == 1)  {
		if (dir_ == 1) fele = pmin_;
		else fele = pmax_;
		for (i=0,qep=fele;i<nel&&qep!=NULL;qep=qep->nxtset(),i++)
			qep->set()->show("");
	}
	if (d1 == -1) {
		if (dir_ == 1)  lele = pmax_;
		else lele = pmin_;
		for (i=0,qep=lele;i<nel&&qep!=NULL;qep=qep->prvset(),i++) {
			qep->set()->show("");
		}
	}
	cout << "\n";
}
*/

void subsetq::save(int *bvar,double *bcrtval,unsigned nel,unsigned dim,int d1)
{
	unsigned i;
	sbsetqe *fele,*lele,*qep;

	if (d1 == 1)  {
		if (dir_ == 1)  fele = pmin_;
		else fele = pmax_;
		for (i=0,qep=fele;i<nel&&qep!=NULL;qep=qep->nxtset(),i++)
			qep->set()->save(&(bvar[i*dim]),bcrtval++,dim);
	}
	if (d1 == -1) {
		if (dir_ == 1)  lele = pmax_;
		else lele = pmin_;
		for (i=0,qep=lele;i<nel&&qep!=NULL;qep=qep->prvset(),i++) {
			qep->set()->save(&(bvar[i*dim]),bcrtval++,dim);
		}
	}
}
