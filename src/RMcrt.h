#ifndef RMDATA
#define RMDATA

#include <deque> 
#include "SpecialArrays.h"

using std::deque;

class rmgdata : public globaldata {
	public:
		rmgdata(vind);
		real*   gettmpv(void)	{ return &tmpv[0]; }
		real	trs(void)	{ return trs_; }
		void	settrs(real ts)	{ trs_ = ts;  }
	private:
		vind		p;  
		real		trs_;
		vector<real>	tmpv;	
};

class rmdata :  public subsetdata {
	public:
		rmdata(vind,vind,vind,rmgdata *,const deque<bool> &,real criterion=0.);
		virtual ~rmdata(void);
		virtual real criterion(void)	const	{ return crt;  }
		virtual real indice(void)	const	{ return sqrt(1.-crt/gdt->trs()); } 
		virtual real updatecrt(vind *,vind *,vind,vind) const;
		virtual void pivot(vind,vind,vind,vind,real,vind *,vind *,subsetdata *,bool); 
/*
Note: subsetdata pointer must point to rmgdata class or unpredictable behaviour will result (general subsetdata class was used in order to garantee upward compability)
*/
		virtual subsetdata *crcopy(vind lastvariab,vind partialnv) const
			{  return new rmdata(lastvariab,partialnv,p,gdt,varin,crt);  }
		virtual void setorgvarl(vind *) {  }
		virtual const real*	getbnds(void)	const	{ return 0; }	
		void setcoefmatel(vind i,vind j,real val)     	{ (*e)(i,j) = val;  }
		void setcrt(real val)			  	{ crt = val; }
		rmgdata*	getgdata(void)	const		{ return gdt;  }
	private:
		vind			lastv;
		vind			p;
		vind			k;
		real			crt;
		deque<bool>		varin; 
		symtwodarray*		e;
		vector<matvectarray *>	ovct;
		rmgdata*		gdt;
};

#endif
