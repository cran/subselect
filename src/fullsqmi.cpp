#include <math.h>
#include <assert.h>
#include "fullsqmi.h"

const double    BIG = 9999E99;

fwrksqm::fwrksqm(unsigned dim)
	:  wrksqmat(dim)
{
	unsigned i,j;

	data = new double[n*p];
	ri = new unsigned[n];
	ir = new unsigned[n];
	for (i=0;i<n;i++) ri[i] = ir[i] = i+1;
	ci = new unsigned[p];
	ic = new unsigned[p];
	for (j=0;j<p;j++) ci[j] = ic[j] = j+1;
}

fwrksqm::~fwrksqm(void)
{
	delete[] data;
	delete[] ri;
	delete[] ci;
	delete[] ir;
	delete[] ic;
}

void fwrksqm::rewind(void)
{
	unsigned i,j;

	for (i=0;i<n;i++) ri[i] = ir[i] = i+1;
	for (j=0;j<p;j++) ci[j] = ic[j] = j+1;
}

double fwrksqm::v(unsigned i,unsigned j)
{
	assert(i <= n && j <= p);
	return data[p*(getrow(i)-1)+getcol(j)-1];
}

void fwrksqm::put(double x,unsigned i,unsigned j)
{
	assert(i <= n && j <= p);
	data[p*(getrow(i)-1)+getcol(j)-1] = x;
}

void fwrksqm::swprows(unsigned row1,unsigned row2)
{
	unsigned t1,t2;

	assert(row1 <= n && row2 <= n);
	ir[(t1=ri[row1-1])-1] = row2;
	ir[(t2=ri[row2-1])-1] = row1;
	ri[row1-1] = t2;
	ri[row2-1] = t1;
}

void fwrksqm::swpcols(unsigned col1,unsigned col2)
{
	unsigned t1,t2;

	assert(col1 <= p && col2 <= p);
	ic[(t1=ci[col1-1])-1] = col2;
	ic[(t2=ci[col2-1])-1] = col1;
	ci[col1-1] = t2;
	ci[col2-1] = t1;
}

void fwrksqm::asgrows(unsigned frow,unsigned nrows,unsigned *rows)
{
	unsigned i,*t;

	assert(frow + nrows - 1 <= n);
	t = new unsigned[nrows];
	for (i=0;i<nrows;i++)  ir[(t[i]=ri[rows[i]-1])-1] = frow+i;
	for (i=0;i<nrows;i++)  ri[frow+i-1] = t[i];
	delete[] t;
}

void fwrksqm::asgcols(unsigned fcol,unsigned ncols,unsigned *cols)
{
	unsigned j,*t;

	assert(fcol + ncols - 1 <= p);
	t = new unsigned[ncols];
	for (j=0;j<ncols;j++)  ic[(t[j]=ci[cols[j]-1])-1] = fcol+j;
	for (j=0;j<ncols;j++)  ci[fcol+j-1] = t[j];
	delete[] t;
}

unsigned fwrksqm::getrow(unsigned rowi)
{
	assert(rowi <= n);
	return ri[rowi-1];
}

unsigned fwrksqm::getcol(unsigned coli)
{
	assert(coli <= p);
	return ci[coli-1];
}

unsigned fwrksqm::getrowp(unsigned orow)
{
	assert(orow <= n);
	return ir[orow-1];
}

unsigned fwrksqm::getcolp(unsigned ocol)
{
	assert(ocol <= p);
	return ic[ocol-1];
}

void fwrksqm::multcol(double sca,unsigned col,unsigned fele,unsigned lele)
{
	unsigned i;

	assert(col <= p && fele <= n && lele <= n);
	if (lele == 0) lele = n;
	for (i=fele;i<=lele;i++)  put(sca*v(i,col),i,col);
}

void fwrksqm::coloper(unsigned col1,double sca,unsigned col2,unsigned frow,unsigned lrow)
{
	unsigned i;

	assert(col1 <= p && col2 <= p && frow <= n && lrow <= n);
	if (lrow == 0) lrow = n;
	for (i=frow;i<=lrow;i++)  put(v(i,col1)+sca*v(i,col2),i,col1);
}

void fwrksqm::multrow(double sca,unsigned row,unsigned fele,unsigned lele)
{
	unsigned j;

	assert(row <= n && fele <= p && lele <= p);
	if (lele == 0) lele = p;
	for (j=fele;j<=lele;j++)  put(sca*v(row,j),row,j);
}

void fwrksqm::rowoper(unsigned row1,double sca,unsigned row2,unsigned fcol,unsigned lcol)
{
	unsigned j;

	assert(row1 <= n && row2 <= n && fcol <= p && lcol <= p);
	if (lcol == 0) lcol = p;
	for (j=fcol;j<=lcol;j++)  put(v(row1,j)+sca*v(row2,j),row1,j);
}

void fwrksqm::maxele(unsigned *row,unsigned *col,double *val,unsigned frow,unsigned fcol)
{
	unsigned i,j,r,c;
	double max,t;

	max = -BIG;
	for (i=frow;i<=n;i++)
		for (j=fcol;j<=p;j++)
			if ( (t=v(i,j)) > max ) {
				r = i;
				c = j;
				max = t;
			}
	(*row) = r;
	(*col) = c;
	(*val) = max;
}

void fwrksqm::minele(unsigned *row,unsigned *col,double *val,unsigned frow,unsigned fcol)
{
	unsigned i,j,r,c;
	double min,t;

	min = BIG;
	for (i=frow;i<=n;i++)
		for (j=fcol;j<=p;j++)
			if ( (t=v(i,j)) < min ) {
				r = i;
				c = j;
				min = t;
			}
	(*row) = r;
	(*col) = c;
	(*val) = min;
}

void fwrksqm::maxabs(unsigned *row,unsigned *col,double *val,unsigned frow,unsigned fcol)
{
	unsigned i,j,r,c;
	double max,t;

	max = 0.;
	for (i=frow;i<=n;i++)
		for (j=fcol;j<=p;j++)
			if ( (t=fabs(v(i,j))) > max ) {
				r = i;
				c = j;
				max = t;
			}
	(*row) = r;
	(*col) = c;
	(*val) = max;
}

void fwrksqm::minabs(unsigned *row,unsigned *col,double *val,unsigned frow,unsigned fcol)
{
	unsigned i,j,r,c;
	double min,t;

	min = BIG;
	for (i=frow;i<=n;i++)
		for (j=fcol;j<=p;j++)
			if ( (t=fabs(v(i,j))) < min ) {
				r = i;
				c = j;
				min = t;
			}
	(*row) = r;
	(*col) = c;
	(*val) = min;
}

ltmatrix::ltmatrix(unsigned dim)
	:  opsqmat(dim)
{
	data = new double[n*(n+1)/2];
}

ltmatrix::~ltmatrix()                     
{
	delete[] data;
}
		
double ltmatrix::v(unsigned i,unsigned j)
{
	assert(i >= j);
	assert(i <= n);
	return data[(unsigned)((i/2.)*(i-1))+j-1];
}

void ltmatrix::put(double x,unsigned i,unsigned j)
{
	assert(i >= j);
	assert(i <= n);
	data[(unsigned)((i/2.)*(i-1))+j-1] = x;
}

unsigned ltmatrix::getrow(unsigned rowi)
{
	assert(rowi <= n);
	return rowi;
}

unsigned ltmatrix::getcol(unsigned coli)
{
	assert(coli <= p);
	return coli;
}

unsigned ltmatrix::getrowp(unsigned orow)
{
	assert(orow <= n);
	return orow;
}

unsigned ltmatrix::getcolp(unsigned ocol)
{
	assert(ocol <= p);
	return ocol;
}

utmatrix::utmatrix(unsigned dim) 
  : ltmatrix(dim)  { }

utmatrix::~utmatrix()
{  }

double utmatrix::v(unsigned i,unsigned j)
{
	assert(i <= j);
	assert(j <= n);
	return ltmatrix::v(j,i);
}

void utmatrix::put(double x,unsigned i,unsigned j)
{
	assert(i <= j);
	assert(j <= n);
	ltmatrix::put(x,j,i);
}

symatrixb::symatrixb(unsigned dim)
	:  ltmatrix(dim)
{
	unsigned i;

	vi = new unsigned[n];
	iv = new unsigned[n];
	for (i=0;i<n;i++) vi[i] = iv[i] = i+1;
}

symatrixb::~symatrixb()
{ 
	delete[] vi; 
	delete[] iv;
}

void symatrixb::rewind(void)
{
	unsigned i;

	for (i=0;i<n;i++) vi[i] = iv[i] = i+1;
}

void symatrixb::swpvar(unsigned var1,unsigned var2)
{
	unsigned t1,t2;

	assert(var1 <= n && var2 <= n);
	iv[(t1=vi[var1-1])-1] = var2;
	iv[(t2=vi[var2-1])-1] = var1;
	vi[var1-1] = t2;
	vi[var2-1] = t1;
}

void symatrixb::asgvar(unsigned fvar,unsigned nvar,unsigned *var)
{
	unsigned i,*t;
	
	assert(fvar + nvar -1 <= n);
	t = new unsigned[nvar];
	for (i=0;i<nvar;i++)  iv[(t[i]=vi[var[i]-1])-1] = fvar+i;
	for (i=0;i<nvar;i++)  vi[fvar+i-1] = t[i];
	delete[] t;
}

unsigned symatrixb::getvar(unsigned vari)
{
	assert(vari <= n);
	return vi[vari-1];
}

unsigned symatrixb::getvarp(unsigned ovar)
{
	assert(ovar <= n);
	return iv[ovar-1];
}

double symatrixb::v(unsigned i,unsigned j)
{
	unsigned i0,j0;

	if ((i0=vi[i-1]) >= (j0=vi[j-1])) return ltmatrix::v(i0,j0);
	else return ltmatrix::v(j0,i0);
}

void symatrixb::put(double x,unsigned i,unsigned j)
{
	unsigned i0,j0;

	if ((i0=vi[i-1]) >= (j0=vi[j-1])) ltmatrix::put(x,i0,j0);
	else ltmatrix::put(x,j0,i0);
}

void symatrixb::cptowsqm(fwrksqm *fsm)
{
	unsigned i,j;

	fsm->det = det;
	for (i=1;i<=n;i++)
		for (j=1;j<=i;j++)  {
			fsm->put(v(i,j),i,j);
			fsm->put(v(i,j),j,i);
		}
}

void symatrixb::cpfrwsqm(fwrksqm *fsm)
{
	unsigned i,j;

	det = fsm->det;
	for (i=1;i<=n;i++)
		for (j=1;j<=i;j++)
			put(fsm->v(i,j),i,j);
}

void symatrixb::multcol(double sca,unsigned col,unsigned fele,unsigned lele)
{
	unsigned i;

	assert(col <= p && fele <= n && lele <= n);
	if (lele == 0) lele = n;
	for (i=fele;i<=lele;i++)  put(sca*v(i,col),i,col);
}

void symatrixb::coloper(unsigned col1,double sca,unsigned col2,unsigned frow,unsigned lrow)
{
	unsigned i;

	assert(col1 <= p && col2 <= p && frow <= n && lrow <= n);
	if (lrow == 0) lrow = n;
	for (i=frow;i<=lrow;i++)  put(v(i,col1)+sca*v(i,col2),i,col1);
}

void symatrixb::multrow(double sca,unsigned row,unsigned fele,unsigned lele)
{
	unsigned j;

	assert(row <= n && fele <= p && lele <= p);
	if (lele == 0) lele = p;
	for (j=fele;j<=lele;j++)  put(sca*v(row,j),row,j);
}

void symatrixb::rowoper(unsigned row1,double sca,unsigned row2,unsigned fcol,unsigned lcol)
{
	unsigned j;

	assert(row1 <= n && row2 <= n && fcol <= p && lcol <= p);
	if (lcol == 0) lcol = p;
	for (j=fcol;j<=lcol;j++)  put(v(row1,j)+sca*v(row2,j),row1,j);
}

void symatrixb::maxdele(unsigned *var,double *val,unsigned fvar)
{
	unsigned i,vi;
	double max,t;

	max = -BIG;
	for (i=fvar;i<=n;i++)
		if ( (t=v(i,i)) > max ) {
			vi = i;
			max = t;
		}
	(*var) = vi;
	(*val) = max;
}

void symatrixb::mindele(unsigned *var,double *val,unsigned fvar)
{
	unsigned i,vi;
	double min,t;

	min = BIG;
	for (i=fvar;i<=n;i++)
		if ( (t=v(i,i)) < min ) {
			vi = i;
			min = t;
		}
	(*var) = vi;
	(*val) = min;
}

void symatrixb::maxdabs(unsigned *var,double *val,unsigned fvar)
{
	unsigned i,vi;
	double max,t;

	max = 0.;
	for (i=fvar;i<=n;i++)
		if ( (t=fabs(v(i,i))) > max ) {
			vi = i;
			max = t;
		}
	(*var) = vi;
	(*val) = max;
}

void symatrixb::mindabs(unsigned *var,double *val,unsigned fvar)
{
	unsigned i,vi;
	double min,t;

	min = BIG;
	for (i=fvar;i<=n;i++)
		if ( (t=fabs(v(i,i))) < min ) {
			vi = i;
			min = t;
		}
	(*var) = vi;
	(*val) = min;
}

void symatrixb::maxabs(unsigned *row,unsigned *col,double *val,unsigned frow,unsigned fcol)
{
	unsigned i,j,r,c;
	double max,t;

	max = 0.;
	for (i=frow;i<=n;i++)
		for (j=fcol;j<=i;j++)
			if ( (t=fabs(v(i,j))) > max ) {
				r = i;
				c = j;
				max = t;
			}
	(*row) = r;
	(*col) = c;
	(*val) = max;
}

void symatrixb::minabs(unsigned *row,unsigned *col,double *val,unsigned frow,unsigned fcol)
{
	unsigned i,j,r,c;
	double min,t;

	min = BIG;
	for (i=frow;i<=n;i++)
		for (j=fcol;j<=i;j++)
			if ( (t=fabs(v(i,j))) < min ) {
				r = i;
				c = j;
				min = t;
			}
	(*row) = r;
	(*col) = c;
	(*val) = min;
}

fullmfact::fullmfact(unsigned dim)  
 :  factrz(dim)  
{
	L = new ltmatrix(dim);
	U = new utmatrix(dim);
}

fullmfact::~fullmfact()
{
	delete L; 
	delete U;
}

fullsqm::fullsqm(unsigned dim) 
  :  sqmatrix(dim) 
{ 
	sqmatrix::data = data = new fwrksqm(dim); 
}

fullsqm::~fullsqm()
{
	delete data;
}

void  cpmatdt(symatrixb *org,symatrixb *res)
{
	unsigned i,j;

	for (i=1;i<=org->n;i++)
		for (j=1;j<=i;j++)
			res->put(org->v(i,j),i,j);
}

void  addmat(symatrixb *m1,symatrixb *m2,symatrixb *res)
{
	unsigned i,j;

	for (i=1;i<=m1->n;i++)
		for (j=1;j<=i;j++)
			res->put(m1->v(i,j)+m2->v(i,j),i,j);
}

void  submat(symatrixb *m1,symatrixb *m2,symatrixb *res)
{
	unsigned i,j;

	for (i=1;i<=m1->n;i++)
		for (j=1;j<=i;j++)
			res->put(m1->v(i,j)-m2->v(i,j),i,j);
}

void  scamult(double sca,symatrixb *org,symatrixb *res)
{
	unsigned i,j;

	for (i=1;i<=org->n;i++)
		for (j=1;j<=i;j++)
			res->put(sca*org->v(i,j),i,j);
}

symatrix::symatrix(unsigned dim)
	:  sqmatrix(dim)
{
	sqmatrix::data = data = new symatrixb(dim);
	inv   = static_cast<symatrixb *>(NULL);
	egval = static_cast<vector *>(NULL);
	egvct = static_cast<vector **>(NULL);
}

symatrix::~symatrix(void)
{
	unsigned j;

	delete inv;
	delete egval;
	if (egvct != NULL) {
		for (j=0;j<rank;j++) delete egvct[j];
		delete egvct;
	}
}

double symatrix::cmpdet(void)
{
	fwrksqm m(n);
	data->cptowsqm(&m);
	return  det = m.trngrp();
}

char symatrix::invert(void)
{
	assert(det != 0.);
	if (inv == NULL) sqmatrix::inv = inv = new symatrixb(n);
	cpmatdt(data,inv);
	if (!inv->invertrp())  {
		det = 0.;
		return FALSE;
	}
	if (det == UNKNOWN) det = 1./inv->det;
	return TRUE;
}

void  cpmatdt(symatrix *org,symatrix *res)
{
	cpmatdt(org->data,res->data);
}

void  addmat(symatrix *m1,symatrix *m2,symatrix *res)
{
	addmat(m1->data,m2->data,res->data);
}

void  submat(symatrix *m1,symatrix *m2,symatrix *res)
{
	submat(m1->data,m2->data,res->data);
}

void  scamult(double sca,symatrix *org,symatrix *res)
{
	scamult(sca,org->data,res->data);
}

void fullmfact::getf(wrksqmat *fm)
{
	unsigned i,j;

	for (i=1;i<=n;i++)  {
		if (i>1) for (j=1;j<i;j++)  L->put(fm->v(i,j),i,j);
		L->put(1.,i,i);
		for (j=i;j<=n;j++) U->put(fm->v(i,j),i,j);
	}
}

