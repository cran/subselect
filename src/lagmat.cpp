#include <math.h>
#include "fullsqmi.h"
#include "fullmati.h"
#include "lagmvct.h"

#define EPSLON 1E-8


lagvct::lagvct(unsigned k,unsigned nvar) 
  :  fullvct(k), lag_(nvar-k)   {  }

double lagvct::v(unsigned  i)
{
	return fullvct::v(i-lag_);
}

void lagvct::put(double x,unsigned i)
{
	fullvct::put(x,i-lag_);
}

kmatrix::kmatrix(unsigned k,unsigned nvar) 
  :  symatrixb(k), lag_(nvar-k), tmpvi(static_cast<unsigned *>(NULL)), 
     old(static_cast<double *>(NULL)) 
{  }

kmatrix::~kmatrix()
{
	delete[] tmpvi;
	delete[] old;
}

double kmatrix::LUFrp(double *d)
{
	double t;
	unsigned j,k;

	d[0] = v(1,1);
	for (j=1;j<n;j++) {
		if ( (t=v(j,j)) < EPSLON)  
			errmsg("\nProgram terminated because of numeric problems.\n");
		for (k=j+1;k<=n;k++)  rowoper(k,-v(k,j)/t,j,j+1);
		multcol(-1./t,j,j+1);
		if (j>1) d[j-1] = d[j-2]*t;
	}
	if (n>2) for (j=n-1;j>1;j--)
		for (k=1;k<j;k++) coloper(k,v(j,k),j,j+1);
	return (det = d[n-1] = d[n-2]*v(n,n));
}

void kmatrix::multcol(double sca,unsigned col,unsigned fele,unsigned lele)
{
	unsigned i;

	if (lele == 0) lele = n;
	for (i=fele;i<=lag_+lele;i++)  put(sca*v(i,col),i,col);
}

void kmatrix::rowoper(unsigned row1,double sca,unsigned row2,unsigned fcol,unsigned lcol)
{
	unsigned j;

	if (lcol == 0) lcol = p;
	for (j=fcol;j<=lag_+lcol;j++)  put(v(row1,j)+sca*v(row2,j),row1,j);
}

void kmatrix::coloper(unsigned col1,double sca,unsigned col2,unsigned frow,unsigned lrow)
{
	unsigned i;

	if (lrow == 0) lrow = n;
	for (i=frow;i<=lag_+lrow;i++)  put(v(i,col1)+sca*v(i,col2),i,col1);
}

double kmatrix::v(unsigned i,unsigned j)
{
	return symatrixb::v(i-lag_,j-lag_);
}

void kmatrix::put(double x,unsigned i,unsigned j)
{
	symatrixb::put(x,i-lag_,j-lag_);
}

void kmatrix::siminvrp(void)
{
	double t,t1;
	unsigned j,k;

	if (tmpvi == NULL) tmpvi = new unsigned[n];
	if (old == NULL) old = new double[n];

	for (j=1;j<=n;j++) {
		maxdabs(tmpvi+j-1,&t,j);
		t = v(tmpvi[j-1],tmpvi[j-1]);
		if (tmpvi[j-1] != j )	swpvar(tmpvi[j-1],j);
		for (k=1;k<j;k++)  old[k-1] = v(k,j);
		for (k=j+1;k<=n;k++)  old[k-1] = v(j,k);
		t = -1 / t;
		if (j > 1) multcol(t,j,1,j-1);
		put(t,j,j);
		if (j < n) multcol(t,j,j+1,n);
		for (k=1;k<=n;k++)  {
			if (k == j) continue;
			t1 = old[k-1];
			if (k < j) {
				coloper(k,t1,j,k,j-1);
				coloper(k,t1,j,j+1,n);
			}
			else coloper(k,t1,j,k,n);
		}
		if (fabs(t) < EPSLON)  errmsg("Error: Trying to invert singular matrix");
		if (det != NOTDF) {
			if (j == 1) det = -t;
			else  det *= -t;
		}
	}
	for (j=n;j>0;j--)  if (tmpvi[j-1] != j)	swpvar(tmpvi[j-1],j);
	return;
}

