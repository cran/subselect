#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sqmatb.h"

const double EPSLON = 1E-8;


/********
 *
 *      trngrp - Triangulizes a matrix by the Gaussian algorithm
 *
 **********/

double wrksqmat::trngrp(void)
{

		  double t,d;
		  unsigned j,k,r,c;
		  int sgn;

		  /* initialize the determinant */

		  if (det == UNKNOWN) {
				d = 1.0;
				sgn = 1;
		  }

		  /* loop through the matrix n times, performing the
			  Gaussian elimination calculations */

		  for (j=1;j<=n;j++) {

			 /* test for singularity and then adjust the determinant */

				maxabs(&r,&c,&t,j,j);
				if (t <  EPSLON) return (det = 0.);
				t = v(r,c);
				if (det == UNKNOWN)  d *= t;

			/*  permute matrix so that element j,j has the maximum abs value  */

				if (r != j )	{
					swprows(r,j);
					sgn *= -1;
				}
				if (c != j) 	{
					swpcols(c,j);
					sgn *= -1;
				}

				if (j<n)  {

			/* now pivot on the element j,j */

					for (k=j+1;k<=n;k++) {
						rowoper(k,-v(k,j)/t,j,j+1);
						put(0,k,j);
					}
				}
		  }

		  /*  Adjust determinant  */

		  if (det == UNKNOWN) det = sgn*d;

		  return det;
}