#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sqmatb.h"
#include <iostream.h>


#define EPSLON 1E-8


/********
 *
 *      wrksqmat::invertrp - Replaces a matrix by its inverse
 *                 			  Uses the Gauss-Jordan algorithm
 *
 **********/

char wrksqmat::invertrp(void)
{

		  double t,t1,d;
		  unsigned j,k,*r,*c;
		  int sgn;

		  /* Prepare auxiliary variables */

			r = new unsigned[n];
			c = new unsigned[n];
			if (det == UNKNOWN) sgn = 1;

		  /* loop through the matrix n times, performing the
			  Gauss-Jordan elimination calculations */

		  for (j=1;j<=n;j++) {

			 /* test for singularity and then adjust the determinant */

				maxabs(r+j-1,c+j-1,&t,j,j);
				if (t < EPSLON) {
					det = 0.;
					delete[] r;
					delete[] c;
					return FALSE;
				}
				t = v(r[j-1],c[j-1]);
				if (det == UNKNOWN)
					if (j == 1) d = t;
					else  d *= t;

			/*  permute matrix so that element j,j has the maximum abs value  */

				if (r[j-1] != j )	{
					swprows(r[j-1],j);
					sgn *= -1;
				}
				if (c[j-1] != j) 	{
					swpcols(c[j-1],j);
					sgn *= -1;
				}

			/* now pivot on the element j,j */

				t = 1 / t;
				multcol(-t,j);
				put(t,j,j);
				for (k=1;k<=n;k++) {
					if (k == j) continue;
					t1 = v(j,k);
					coloper(k,t1,j);
					put(t1*t,j,k);
				}
		  }

		  /*  permute matrix back  */

		  for (j=n;j>0;j--)  {
				if (c[j-1] != j)	swprows(c[j-1],j);
				if (r[j-1] != j) 	swpcols(r[j-1],j);
		  }

		  /*  Adjust determinant and deallocate memory */

		  if (det == UNKNOWN) det = sgn*(1./d);
		  else det = 1./det;
		  delete[] r;
		  delete[] c;

		  return TRUE;
}

