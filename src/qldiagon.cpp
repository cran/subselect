#include <stdlib.h>
#include <math.h>
#include "fullmati.h"
#include "fullsqmi.h"

#define EPSLON  1E-12

/********
 *
 *      diag --  Finds the eigenvalues and eigenvectors of a symmetric matrix
 *               using a tridiagonal reduction followed by the QL algorithm
 **********/

typedef  vector* pvector;
int cmptmp(const void *,const void *);
void triad(fullmat *,unsigned,fullvct *, fullvct *);
void ql(vector *, fullvct *,unsigned,fullmat *);
void errmsg(char *);
inline double sign(double a,double b)  { return b > 0 ? fabs(a) : -fabs(a); }

fullvct	  *tmp;

void symatrix::diag(unsigned r)
{

	unsigned  i,j,*dmyv;

	fullmat m0(n,n);						                             // Duplicate of original matrix 
    fullvct wv(n);                                                       // Storage allocation for 'intermediate' vectors 
	tmp = new fullvct(n);						                         
	egval = new fullvct(r);						                         // Storage allocation for eigenvalues
	egvct = new pvector[r];                                              // Storage allocation for eigenvectors
	for (i=0;i<r;i++) egvct[i] = new fullvct(p);	        

	rank = r;
	cpmatdt(data,&m0);
	triad(&m0,n,tmp,&wv);                                                // Triangular decomposition 
	ql(tmp,&wv,n,&m0);                                                   // Reduction of sym. trid. matrix 
	dmyv = new unsigned[n];												 // Sort in descending order of eigenvalues
	for (i=1;i<=n;i++)  dmyv[i-1] = i;
	qsort((void *)dmyv,n,sizeof(*dmyv),cmptmp);
	for (i=1;i<=r;i++)  {
		egval->put(tmp->v(dmyv[i-1]),i);                                 // Save eigenvalues
		for (j=1;j<=n;j++)  egvct[i-1]->put(m0.v(j,dmyv[i-1]),j);        // Save eigenvectors
	}

	delete[] dmyv;                                                       // Deallocate no longer needed sotorage space
	delete tmp;
}

int cmptmp(const void *a,const void *b)
{
	unsigned ai,bi;

	ai = *(unsigned *)a;
	bi = *(unsigned *)b;
	if (tmp->v(ai) < tmp->v(bi)) return 1;
		else if (tmp->v(ai) > tmp->v(bi)) return -1;
			else return 0;
}



/********
 *
 *      triad --  Householder Reduction a real, symmetric matrix to a symmetric, tridiagonal form 
 *                Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
 *                Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
 *                Springer-Verlag, 1976, pp. 489-494.
 *                W H Press et al., Numerical Recipes in C, Cambridge U P,
 *                1988, pp. 373-374.  
********/

void triad(fullmat *a, unsigned n, fullvct *d, fullvct *e)
{
	unsigned l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n;i>=2;i--)  {
		l = i - 1;
		h = scale = 0.0;
		if (l > 1)  {
			for (k=1;k<=l;k++)  scale += fabs(a->v(i,k));
			if (scale < EPSLON)  e->put(a->v(i,l),i);
			else  {
				for (k=1;k<=l;k++)  {
					a->put(a->v(i,k)/scale,i,k);
					h += a->v(i,k) * a->v(i,k);
				}
				f = a->v(i,l);
				g = (f>=0.0 ? -sqrt(h) : sqrt(h));
				e->put(scale*g,i);
				h -= f*g;
				a->put(f-g,i,l);
				f = 0.0;
				for (j=1;j<=l;j++)  {
					a->put(a->v(i,j)/h,j,i);
					g = 0.0;
					for (k=1;k<=j;k++)    g += a->v(j,k)*a->v(i,k);
					for (k=j+1;k<=l;k++)  g += a->v(k,j)*a->v(i,k);
					e->put(g/h,j);
					f += e->v(j)*a->v(i,j);
				}
				hh = f / (h + h);
				for (j=1;j<=l;j++)  {
					f = a->v(i,j);
					e->put(g = e->v(j)-hh*f,j);
					for (k=1;k<=j;k++)  a->put(a->v(j,k)-(f*e->v(k)+g*a->v(i,k)),j,k);
				}
			}
		}
		else  e->put(a->v(i,l),i);
		d->put(h,i);
    }
	d->put(0.0,1);
	e->put(0.0,1);
	for (i=1;i<=n;i++)  {
		l = i - 1;
		if (d->v(i))  {
			for (j=1;j<=l;j++)  {
				g = 0.0;
				for (k=1;k<=l;k++)  g += a->v(i,k)*a->v(k,j);
				for (k=1;k<=l;k++)  a->put(a->v(k,j)-g*a->v(k,i),k,j);
			}
		}
		d->put(a->v(i,i),i);
		a->put(1.0,i,i);
		for (j=1;j<=l;j++)  {
			a->put(0.0,i,j);
			a->put(0.0,j,i);
		}
	}
}

/********
 *
 *      ql --     Finds the eigenvalues and eigenvectors of a real symmetric, tridiagonal
 *                matrix by the QL algorithm
********/


void ql(vector *d, fullvct *e,unsigned n,fullmat *z)
{
	const unsigned MAXITER = 30;
	unsigned m,l,iter,i,k;
	double s, r, p,g,f,dd,c,b;

	for (i=2;i<=n;i++)  e->put(e->v(i),i-1);
	e->put(0.0,n);
	for (l=1;l<=n;l++)  {
		iter = 0;
		do
		{
			for (m=l;m<=n-1;m++)  {
				dd = fabs(d->v(m)) + fabs(d->v(m+1));
				if (fabs(e->v(m)) < EPSLON) break;
			}
			if (m != l)  {
				if (iter++ == MAXITER) errmsg("No convergence in QL");
				g = (d->v(l+1) - d->v(l)) / (2.0 * e->v(l));
				r = sqrt((g * g) + 1.0);
				g = d->v(m) - d->v(l) + e->v(l) / (g + sign(r, g));
				s = c = 1.0;
				p = 0.0;
				for (i=m-1;i>=l;i--)  {
					f = s * e->v(i);
					b = c * e->v(i);
					if (fabs(f) >= fabs(g))  {
						c = g / f;
						r = sqrt((c * c) + 1.0);
						e->put(f*r,i+1);
						c *= (s = 1.0/r);
					}
					else  {
						s = f / g;
						r = sqrt((s * s) + 1.0);
						e->put(g*r,i+1);
						s *= (c = 1.0/r);
                    }
					g = d->v(i+1) - p;
					r = (d->v(i) - g) * s + 2.0 * c * b;
					p = s * r;
					d->put(g+p,i+1);
					g = c * r - b;
					for (k=1;k<=n;k++)  {
						f = z->v(k,i+1);
						z->put(s*z->v(k,i)+c*f,k,i+1);
						z->put(c*z->v(k,i)-s*f,k,i);
                    }
                 }
                 d->put(d->v(l)-p,l);
                 e->put(g,l);
                 e->put(0.0,m);
             }
        }  while (m != l);
    }
 }


