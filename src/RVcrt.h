#ifndef RVDATA
#define RVDATA

#include <deque> 
#include "SpecialArrays.h"

typedef vector< vector<real> >  twodarray; 
using std::deque;

class rvgdata : public globaldata {
	public:
		rvgdata(vind);
		virtual		~rvgdata(void);
		real*		gettmpv(void)			{ return &tmpv[0];  }
		real*		getcndv(void)			{ return &cndv[0]; }
		real		trs2(void)	   const	{ return trs2_; }
		void		settrs2(real ts2)		{ trs2_ = ts2;  }
		twodarray&	getm1t(void)  			{ return m1t; }
		void		sets2(vind i,vind j,real val)	{ (*s2)(i,j) = val; }
		real		gets2(vind i,vind j) const	{ return (*s2)(i,j); }
	private:
		vind		p;  
		vector<real>	tmpv;
		vector<real>	cndv;
		symtwodarray*	s2;
		real		trs2_;
		twodarray	m1t;
};

class rvdata :  public subsetdata {
	public:
		rvdata(vind,vind,vind,rvgdata *,const deque<bool>&,vind *,real criterion);
		virtual ~rvdata(void);
		virtual real criterion(void)	const	{ return crt;  }
		virtual real indice(void)	const	{ return sqrt(crt/gdt->trs2()); } 
		virtual real updatecrt(vind *,vind *,vind,vind)	const;
		virtual void pivot(vind,vind,vind,vind,real,vind *,vind *,subsetdata *,bool);  
/*
Note: subsetdata pointer must point to rvgdata class or unpredictable behaviour will result (general subsetdata class was used in order to garantee upward compability)
*/
		virtual subsetdata *crcopy(vind lastvariab,vind partialnv) const
			{  return new rvdata(lastvariab,partialnv,p,gdt,varin,orgvar,crt);	}
		virtual void setorgvarl(vind* list)		{ orgvar = list; }
		virtual const real*	getbnds(void)	const	{ return 0; }
		void setcoefmatel(vind i,vind j,real val)	{ (*e)(i,j) = val; }
		void setcrt(real val)				{ crt = val; }
		rvgdata *getgdata(void) const			{ return gdt; }
		void  sets2m1(vind i,vind j,real val)		{ s2m1[i][j] = val; }
		real  gets2m1(vind i,vind j) const		{ return s2m1[i][j]; }
		void  cmpts2sm1(vind *,vind *,twodarray &,bool *,vind *,vind,vind,vind) const;  
/* Computes the S2*S^1 matrix product for the sub-matrices defined by the boolean list */
		real  frobenius(twodarray&,bool *) const;   
/* Computes the Frobenius norm for the sub-matrix defined by the boolean list  */
	private:
		vind			lastv;
		vind			p;
		vind			k;
		real			crt;
		deque<bool>		varin; 
		vind*			orgvar;
		symtwodarray*		e;
		vector<matvectarray *>	ivct;
		twodarray		s2m1;
		rvgdata*		gdt;
};

#endif
