#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream.h>
#include "fullsqmi.h"

const double EPSLON = 1E-8;

/********
 *
 *      symatrixb::invertrp - Replaces a symmetric matrix by its inverse
 *                 			   Uses an implementation of the Gauss-Jordan
 *                            algorithm based on symmetric sweep operators
 *
 **********/


char symatrixb::invertrp(void)
{

		  double d,t,t1,*old;
		  unsigned j,k,l,*vi;

		  /* initialize  auxiliary variables  */

		  vi = new unsigned[n];
		  old = new double[n];

		  /* loop through the matrix n times, performing the
			  Gauss-Jordan elimination calculations */

		  for (j=1;j<=n;j++) {

			 /* test for singularity, adjust the determinant and then
				 permute matrix so that element j,j is the diagonal
				 element with maximum absolute value                      */

				maxdabs(vi+j-1,&t,j);
				if (t < EPSLON) {
					det = 0.;
					delete[] vi;
					delete[] old;
					return FALSE;
				}
				t = v(vi[j-1],vi[j-1]);
				if (det == UNKNOWN)
					if (j == 1) d = t;
					else  d *= t;
				if (vi[j-1] != j )	swpvar(vi[j-1],j);

			/* now pivot on the element j,j */

				for (k=1;k<j;k++)  old[k-1] = v(k,j);
				for (k=j+1;k<=n;k++)  old[k-1] = v(j,k);
				if (j < n) t = -1 / t;
				else t = 1/t;
				if (j > 1) multcol(t,j,1,j-1);
				put(t,j,j);
				if (j < n) multcol(t,j,j+1,n);
				for (k=1;k<=n;k++)  {
					if (k == j) continue;
					t1 = old[k-1];
					if (j < n) {
						if (k < j) {
							coloper(k,t1,j,k,j-1);
							coloper(k,t1,j,j+1,n);
						}
						else coloper(k,t1,j,k,n);
					}
					else  for (l=k;l<n;l++)  put(-v(l,k)+t1*v(l,j),l,k);
				}
		  }

		  /*  permute matrix back, adjust determinant and deallocate memory */

		  for (j=n;j>0;j--)  if (vi[j-1] != j)	swpvar(vi[j-1],j);
		  if (det == UNKNOWN) det = 1./d;
		  else det = 1./det;
		  delete[] vi;
		  delete[] old;

		  return TRUE;
}

