#ifndef QFDATA
#define QFDATA

#include "SpecialArrays.h"

using std::vector;


class qfauxmem : public auxmemory {
	public:
		explicit qfauxmem(vind);
		virtual ~qfauxmem(void);
		real*   gettmpv(void)      { return &tmpv[0]; }
	protected:
		vind		r;  	
		vector<real>	tmpv;  	
};


class qfdata :  public subsetdata {
	public:
		qfdata(vind,vind,vind,qfauxmem *,real criterion);
		virtual ~qfdata(void);
		virtual real criterion(void)  const	{ return crt;  }
		virtual real updatecrt(vind *,vind *,vind,vind) const ;
		virtual void pivot(vind,vind,vind,vind,real,vind *,vind *,subsetdata *,bool);  
/* 
Note: subsetdata pointer must point to qfdata class or unpredictable behaviour will result  (general subsetdata class was used in order to garantee upward compability)
*/
		virtual subsetdata *crcopy(vind totalnv,vind partialnv)  const
			{  return new qfdata(totalnv,partialnv,r,auxmem,crt);  }
		virtual void setorgvarl(vind *) {  }
		virtual const real*	getbnds(void)    	const   { return 0; }	 
		void setvectel(vind i,vind j,real val)			{ ve[i][j] = val; }
		void setcoefmatel(vind i,vind j,real val)		{ (*e)(i,j) = val;  }
		void setcrt(real val)					{ crt = val; }
		qfauxmem*	getauxmem(void)			  const	{ return auxmem;  }
	protected:
		vind			p;
		vind			k;
		vind			r;
		real			crt;
		vector< vector<real> >	ve;
		symtwodarray*		e;
		qfauxmem*		auxmem;
};

#endif
