#include <math.h>
#include "fullmati.h"
#include "fullsqmi.h"


#define EPSLON 1E-8

/********
 *
 *      diag --  Finds the eigenvalues and eigenvectors
 *					  of a symmetric matrix using the power method
 *
 **********/

typedef  vector* pvector;

void symatrix::diag(unsigned r)
{

	unsigned  i,j,d;
	double lbd,lbd0;
	fullmat m0(n,n);
	fullvct wv(n);
	fullvct wv1(n);

	rank = r;
	egval = new fullvct(r);
	egvct = new pvector[r];
	for (d=0;d<r;d++) egvct[d] = new fullvct(p);
	cpmatdt(data,&m0);
	for (d=1;d<=r;d++)  {
		for (j=1;j<=n;j++) wv.put(1./sqrt(n),j);
		lbd = 1.;
		for (;;)  {
			lbd0 = lbd;
			rightmult(&m0,&wv,&wv1);
			lbd = vectprod(&wv,&wv1);
			scamult(1./sqrt(vectprod(&wv1,&wv1)),&wv1,&wv);
			if  ( fabs(lbd-lbd0) <  EPSLON) break;
		}
		cpvectdt(&wv,egvct[d-1]);
		egval->put(lbd,d);
		for (i=1;i<=n;i++)
			for (j=1;j<=n;j++)  m0.put(m0.v(i,j)-lbd*wv.v(i)*wv.v(j),i,j);
	}
}

