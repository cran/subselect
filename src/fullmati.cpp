#include <math.h>
#include <assert.h>
#include "fullmati.h"

const double    BIG     = 9999E99;

fullvct::fullvct(unsigned dim)    
  :  vector(dim)
{
	unsigned i;

	data = new double[n];
	ri = new unsigned[n];
	ir = new unsigned[n];
	for (i=0;i<n;i++) ri[i] = ir[i] = i+1;
}

fullvct::~fullvct()
{ 
	delete[] data; 
	delete[] ri; 
	delete[] ir; 
}

void fullvct::rewind()
{
	unsigned i;

	for (i=0;i<n;i++) ri[i] = ir[i] = i+1;
}

double fullvct::v(unsigned i)
{
	assert(i<=n);
	return data[ri[i-1]-1];
}

void fullvct::put(double x,unsigned i)
{
	assert(i<=n);
	data[ri[i-1]-1] = x;
}

void fullvct::swprows(unsigned row1,unsigned row2)
{
	unsigned t1,t2;

	assert(row1 <= n && row2 <= n);
	ir[(t1=ri[row1-1])-1] = row2;
	ir[(t2=ri[row2-1])-1] = row1;
	ri[row1-1] = t2;
	ri[row2-1] = t1;
}

void fullvct::asgrows(unsigned frow,unsigned nrows,unsigned *rows)
{
	unsigned i,*t;

	assert(frow + nrows - 1 <= n);
	t = new unsigned[nrows];
	for (i=0;i<nrows;i++)  ir[(t[i]=ri[rows[i]-1])-1] = frow+i;
	for (i=0;i<nrows;i++)  ri[frow+i-1] = t[i];
	delete[] t;
}

unsigned fullvct::getrow(unsigned rowi)
{
	assert(rowi<=n);
	return ri[rowi-1];
}

unsigned fullvct::getrowp(unsigned orow)
{
	assert(orow<=n);
	return ir[orow-1];
}

double fullvct::maxabs(unsigned fele)
{
	unsigned j;
	double t,max;

	assert(fele<=n);
	max = fabs(v(fele));
	for (j=fele+1;j<=n;j++) if ( (t=fabs(v(j))) > max) max = t;
	return max;
}

fullmat::fullmat(unsigned nobs,unsigned nvar)
	:  wrkmatrix(nobs,nvar)
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

fullmat::~fullmat(void)
{
	delete[] data;
	delete[] ri;
	delete[] ci;
	delete[] ir;
	delete[] ic;
}

void fullmat::rewind(void)
{
	unsigned i,j;

	for (i=0;i<n;i++) ri[i] = ir[i] = i+1;
	for (j=0;j<p;j++) ci[j] = ic[j] = j+1;
}

double fullmat::v(unsigned i,unsigned j)
{
	assert(i <= n && j <= n);
	return data[p*(ri[i-1]-1)+ci[j-1]-1];
}

void fullmat::put(double x,unsigned i,unsigned j)
{
	assert(i <= n && j <= p);
	data[p*(ri[i-1]-1)+ci[j-1]-1] = x;
}

void fullmat::swprows(unsigned row1,unsigned row2)
{
	unsigned t1,t2;

	assert(row1 <= n && row2 <= n);
	ir[(t1=ri[row1-1])-1] = row2;
	ir[(t2=ri[row2-1])-1] = row1;
	ri[row1-1] = t2;
	ri[row2-1] = t1;
}

void fullmat::swpcols(unsigned col1,unsigned col2)
{
	unsigned t1,t2;

	assert(col1 <= n && col2 <= n);
	ic[(t1=ci[col1-1])-1] = col2;
	ic[(t2=ci[col2-1])-1] = col1;
	ci[col1-1] = t2;
	ci[col2-1] = t1;
}

void fullmat::asgrows(unsigned frow,unsigned nrows,unsigned *rows)
{
	unsigned i,*t;

	assert(frow + nrows - 1 <= n);
	t = new unsigned[nrows];
	for (i=0;i<nrows;i++)  ir[(t[i]=ri[rows[i]-1])-1] = frow+i;
	for (i=0;i<nrows;i++)  ri[frow+i-1] = t[i];
	delete[] t;
}

void fullmat::asgcols(unsigned fcol,unsigned ncols,unsigned *cols)
{
	unsigned j,*t;

	assert(fcol + ncols -1 <= p);
	t = new unsigned[ncols];
	for (j=0;j<ncols;j++)  ic[(t[j]=ci[cols[j]-1])-1] = fcol+j;
	for (j=0;j<ncols;j++)  ci[fcol+j-1] = t[j];
	delete[] t;
}

unsigned fullmat::getrow(unsigned rowi)
{
	assert(rowi <= n);
	return ri[rowi-1];
}

unsigned fullmat::getcol(unsigned coli)
{
	assert(coli <= p);
	return ci[coli-1];
}

unsigned fullmat::getrowp(unsigned orow)
{
	assert(orow <= n);
	return ir[orow-1];
}

unsigned fullmat::getcolp(unsigned ocol)
{
	assert(ocol <= p);
	return ic[ocol-1];
}

void fullmat::multcol(double sca,unsigned col,unsigned fele)
{
	unsigned i;

	assert(col <= p && fele <= n);
	for (i=fele;i<=n;i++)  put(sca*v(i,col),i,col);
}

void fullmat::coloper(unsigned col1,double sca,unsigned col2,unsigned frow)
{
	unsigned i;

	assert(col1 <= p && col2 <= p && frow <= n);
	for (i=frow;i<=n;i++)  put(v(i,col1)+sca*v(i,col2),i,col1);
}

void fullmat::multrow(double sca,unsigned row,unsigned fele)
{
	unsigned j;

	assert(row <= n && fele <= p);
	for (j=fele;j<=p;j++)  put(sca*v(row,j),row,j);
}

void fullmat::rowoper(unsigned row1,double sca,unsigned row2,unsigned fcol)
{
	unsigned j;

	assert(row1 <= n && row2 <= n && fcol <= p);
	for (j=fcol;j<=p;j++)  put(v(row1,j)+sca*v(row2,j),row1,j);
}

void fullmat::maxele(unsigned *row,unsigned *col,double *val,unsigned frow,unsigned fcol)
{
	unsigned i,j,r,c;
	double max,t;

	max = -BIG;
	for (i=frow;i<=n;i++)
		for (j=fcol;j<=n;j++)
			if ( (t=v(i,j)) > max ) {
				r = i;
				c = j;
				max = t;
			}
	(*row) = r;
	(*col) = c;
	(*val) = max;
}

void fullmat::minele(unsigned *row,unsigned *col,double *val,unsigned frow,unsigned fcol)
{
	unsigned i,j,r,c;
	double min,t;

	min = BIG;
	for (i=frow;i<=n;i++)
		for (j=fcol;j<=n;j++)
			if ( (t=v(i,j)) < min ) {
				r = i;
				c = j;
				min = t;
			}
	(*row) = r;
	(*col) = c;
	(*val) = min;
}

void fullmat::maxabs(unsigned *row,unsigned *col,double *val,unsigned frow,unsigned fcol)
{
	unsigned i,j,r,c;
	double max,t;

	max = 0.;
	for (i=frow;i<=n;i++)
		for (j=fcol;j<=n;j++)
			if ( (t=fabs(v(i,j))) > max ) {
				r = i;
				c = j;
				max = t;
			}
	(*row) = r;
	(*col) = c;
	(*val) = max;
}

void fullmat::minabs(unsigned *row,unsigned *col,double *val,unsigned frow,unsigned fcol)
{
	unsigned i,j,r,c;
	double min,t;

	min = BIG;
	for (i=frow;i<=n;i++)
		for (j=fcol;j<=n;j++)
			if ( (t=fabs(v(i,j))) < min ) {
				r = i;
				c = j;
				min = t;
			}
	(*row) = r;
	(*col) = c;
	(*val) = min;
}

void  cpvectdt(vector *org,vector *res)
{
	unsigned j;

	assert(org && res);
	assert(org->n == res->n);
	for (j=1;j<=org->n;j++)   res->put(org->v(j),j);
	return;
}

void  scamult(double sca,vector *org,vector *res)
{
	unsigned j;

	assert(org && res);
	assert(org->n == res->n);
	for (j=1;j<=org->n;j++)	res->put(sca*org->v(j),j);
	return;
}

void  addvect(vector *v1,vector *v2,vector *res)
{
	unsigned j;

	assert(v1 && v2 && res);
	assert(v1->n == v2->n);
	assert(v1->n == res->n);
	for (j=1;j<=v1->n;j++)	res->put(v1->v(j)+v2->v(j),j);
	return;
}

void  subvect(vector *v1,vector *v2,vector *res)
{
	unsigned j;

	assert(v1 && v2 && res);
	assert(v1->n == v2->n);
	assert(v1->n == res->n);
	for (j=1;j<=v1->n;j++)	res->put(v1->v(j)-v2->v(j),j);
	return;
}

void  leftmult(vector *v,matrix *m,vector *res)
{
	unsigned j,k;

	assert(v != res);
	assert(v && m && res);
	assert(v->n == m->n);
	assert(m->p == res->n);
	for (j=1;j<=m->p;j++)  {
		res->put(v->v(1)*m->v(1,j),j);
		if (v->n > 1) for (k=2;k<=v->n;k++)
			res->put(res->v(j)+v->v(k)*m->v(k,j),j);
	}
	return;
}

void  rightmult(matrix *m,vector *v,vector *res)
{
	unsigned i,j;

	assert(v != res);
	assert(v && m && res);
	assert(v->n == m->p);
	assert(m->n == res->n);
	for (i=1;i<=m->n;i++)  {
		res->put(m->v(i,1)*v->v(1),i);
		if (m->p > 1) for (j=2;j<=m->p;j++)
			res->put(res->v(i)+m->v(i,j)*v->v(j),i);
	}
	return;
}

double  vectprod(vector *v1,vector *v2)
{
	unsigned j;
	double res;

	assert(v1 && v2);
	assert(v1->n == v2->n);
	res = v1->v(1) * v2->v(1);
	if (v1->n > 1)
		for (j=2;j<=v1->n;j++) res += v1->v(j)*v2->v(j);
	return res;
}


void  cpmatdt(matrix *org,matrix *res)
{
	unsigned i,j;

	assert(org && res);
	assert(org->n == res->n && org->p == res->p);
	for (i=1;i<=org->n;i++)
		for (j=1;j<=org->p;j++)
			res->put(org->v(i,j),i,j);
	return;
}

void  scamult(double sca,matrix *org,matrix *res)
{
	unsigned i,j;

	assert(org && res);
	assert(org->n == res->n && org->p == res->p);
	for (i=1;i<=org->n;i++)
		for (j=1;j<=org->p;j++)
			res->put(sca*org->v(i,j),i,j);
	return;
}

void  transpose(matrix *m,matrix *res)
{
	unsigned i,j;

	assert(m && res);
	assert(m->n == res->p && m->p == res->n);
	for (i=1;i<=m->n;i++)
		for (j=1;j<=m->p;j++)  res->put(m->v(i,j),j,i);
	return;
}

void  addmat(matrix *m1,matrix *m2,matrix *res)
{
	unsigned i,j;

	assert(m1 && m2 && res);
	assert(m1->n == m2->n && m1->p == m2->p);
	assert(m1->n == res->n && m1->p == res->p);
	for (i=1;i<=m1->n;i++)
		for (j=1;j<=m1->p;j++)
			res->put(m1->v(i,j)+m2->v(i,j),i,j);
	return;
}

void  addmtm(matrix *mt,matrix *m,matrix *res)
{
	unsigned i,j;

	assert(mt && m && res);
	assert(mt->p == m->n && mt->n == m->p);
	assert(mt->p == res->n && mt->n == res->p);
	for (i=1;i<=mt->p;i++)
		for (j=1;j<=mt->n;j++)
			res->put(mt->v(j,i)+m->v(i,j),i,j);
	return;
}

void  addmmt(matrix *m,matrix *mt,matrix *res)
{
	unsigned i,j;

	assert(m && mt && res);
	assert(m->n == mt->p && m->p == mt->n);
	assert(m->n == res->n && m->p == res->p);
	for (i=1;i<=m->n;i++)
		for (j=1;j<=m->p;j++)
			res->put(m->v(i,j)+mt->v(j,i),i,j);
	return;
}

void  addmtmt(matrix *mt1,matrix *mt2,matrix *res)
{
	unsigned i,j;

	assert(mt1 && mt2 && res);
	assert(mt1->p == mt2->p && mt1->n == mt2->n);
	assert(mt1->p == res->n && mt1->n == res->p);
	for (i=1;i<=mt1->p;i++)
		for (j=1;j<=mt1->n;j++)
			res->put(mt1->v(j,i)+mt2->v(j,i),i,j);
	return;
}

void  submat(matrix *m1,matrix *m2,matrix *res)
{
	unsigned i,j;

	assert(m1 && m2 && res);
	assert(m1->n == m2->n && m1->p == m2->p);
	assert(m1->n == res->n && m1->p == res->p);
	for (i=1;i<=m1->n;i++)
		for (j=1;j<=m1->p;j++)
			res->put(m1->v(i,j)-m2->v(i,j),i,j);
	return;
}

void  submtm(matrix *mt,matrix *m,matrix *res)
{
	unsigned i,j;

	assert(mt && mt && res);
	assert(mt->p == m->n && mt->n == m->p);
	assert(mt->p == res->n && mt->n == res->p);
	for (i=1;i<=mt->p;i++)
		for (j=1;j<=mt->n;j++)
			res->put(mt->v(j,i)-m->v(i,j),i,j);
	return;
}

void  submmt(matrix *m,matrix *mt,matrix *res)
{
	unsigned i,j;

	assert(m && mt && res);
	assert(m->n == mt->p && m->p == mt->n);
	assert(m->n == res->n && m->p == res->p);
	for (i=1;i<=m->n;i++)
		for (j=1;j<=m->p;j++)
			res->put(m->v(i,j)-mt->v(j,i),i,j);
	return;
}

void  submtmt(matrix *mt1,matrix *mt2,matrix *res)
{
	unsigned i,j;

	assert(mt1 && mt2 && res);
	assert(mt1->p == mt2->p && mt1->n == mt2->n);
	assert(mt1->p == res->n && mt1->n == res->p);
	for (i=1;i<=mt1->p;i++)
		for (j=1;j<=mt1->n;j++)
			res->put(mt1->v(j,i)-mt2->v(j,i),i,j);
	return;
}

void  multmat(matrix *m1,matrix *m2,matrix *res)
{
	unsigned i,j,k;

	assert(m1 && m2 && res);
	assert(m1 != res && m2 != res);
	assert(m1->p == m2->n);
	assert(m1->n == res->n && m2->p == res->p);
	for (i=1;i<=m1->n;i++)
		for (j=1;j<=m2->p;j++) {
				res->put(m1->v(i,1)*m2->v(1,j),i,j);
				if (m1->p > 1) for (k=2;k<=m1->p;k++)
					res->put(res->v(i,j)+m1->v(i,k)*m2->v(k,j),i,j);
		}
	return;
};

void  multmtm(matrix *mt,matrix *m,matrix *res)
{
	unsigned i,j,k;

	assert(mt && mt && res);
	assert(mt != res && m != res);
	assert(mt->n == m->n);
	assert(mt->p == res->n && m->p == res->p);
	for (i=1;i<=mt->p;i++)
		for (j=1;j<=m->p;j++) {
				res->put(mt->v(1,i)*m->v(1,j),i,j);
				if (mt->n > 1) for (k=2;k<=mt->n;k++)
					res->put(res->v(i,j)+mt->v(k,i)*m->v(k,j),i,j);
		}
	return;
};

void  multmmt(matrix *m,matrix *mt,matrix *res)
{
	unsigned i,j,k;

	assert(m && mt && res);
	assert(m != res && mt != res);
	assert(m->p == mt->p);
	assert(m->n == res->n && mt->n == res->p);
	for (i=1;i<=m->n;i++)
		for (j=1;j<=mt->n;j++) {
				res->put(m->v(i,1)*mt->v(j,1),i,j);
				if (m->p > 1) for (k=2;k<=m->p;k++)
					res->put(res->v(i,j)+m->v(i,k)*mt->v(j,k),i,j);
		}
	return;
};

void  multmtmt(matrix *mt1,matrix *mt2,matrix *res)
{
	unsigned i,j,k;

	assert(mt1 && mt2 && res);
	assert(mt1 != res && mt2 != res);
	assert(mt1->n == mt2->p);
	assert(mt1->p == res->n && mt2->n == res->p);
	for (i=1;i<=mt1->p;i++)
		for (j=1;j<=mt2->n;j++) {
				res->put(mt1->v(1,i)*mt2->v(j,1),i,j);
				if (mt1->n > 1) for (k=2;k<=mt1->n;k++)
					res->put(res->v(i,j)+mt1->v(k,i)*mt2->v(j,k),i,j);
		}
	return;
};


