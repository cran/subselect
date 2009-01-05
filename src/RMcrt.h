#ifndef RMDATA
#define RMDATA

#include <deque> 
#include "SpecialArrays.h"

using std::deque;

namespace extendedleaps {

class rmgdata : public globaldata {
	public:
		rmgdata(vind nvar):   	p(nvar) { }
		real	trs(void)	{ return trs_; }
		void	settrs(real ts)	{ trs_ = ts;  }
	private:
		vind	p;  
		real	trs_;
};

class partialrmdata :  public partialdata {                 /* Data used in criterion RM updates   */
	public:
		explicit		partialrmdata(vind);
		virtual			~partialrmdata(void)    {  }	
		real*			gettmpv(void)		{ return &tmpv[0]; }
		virtual const real	getcrt(void) const	{ return crt; }
		const real		getpivotval(void) const	{ return pivotval; }
		void 			setcrt(real c) 		{ crt = c; }
		void 			setpivotval(real pv)	{ pivotval = pv; }
	protected:
		vind		p;  
		real		crt;
		real		pivotval;
		vector<real>	tmpv;  	
	friend class rmdata;
};

class rmdata :  public subsetdata {
	public:
		rmdata(vind lastvariab,vind nvtopiv,vind tnv,rmgdata *data,const deque<bool>& active,real criterion=0.);
		virtual ~rmdata(void);
		virtual const real criterion(void)	const	{ return crt;  }
		virtual void setcriterion(real c)		{ crt = c; }
		virtual const real indice(void)		const	{ return sqrt(1.-crt/gdt->trs()); } 
		virtual real updatecrt(direction d,mindices& mmind,vind var,partialdata* pdt) const;
		virtual void pivot(direction d,mindices& mmind,vind vp,vind t,partialdata* pdt,subsetdata* fdt,bool last);
/*
	Note: subsetdata pointer must point to rmgdata class or unpredictable behaviour will result 
	(general subsetdata class was used in order to garantee upward compability)
*/
		virtual subsetdata *crcopy(vind lastvariab,vind partialnv) const
			{  return new rmdata(lastvariab,partialnv,p,gdt,varin,crt);  }
		virtual void setorgvarl(vind *) {  }
		virtual const real*	getbnds(void)	const	{ return 0; }	
		void setcoefmatel(vind i,vind j,real val)	{ (*e)(i,j) = val;  }
		void setcrt(real val)			  	{ crt = val; }
		rmgdata*	getgdata(void)	const		{ return gdt;  }
	private:
		real updatecrt(direction d,itindex<d>& fmmind,vind var,vind varind,partialdata* newdtpnt) const;   
		real updatecrt(direction d,itindex<i>& fmmind,vind var,vind varind,partialdata* newdtpnt) const;   
		void pivot(direction d,lagindex<d>& prtmmit,itindex<d>& fmmind,vind vp,vind t,partialdata* newpdtpnt,subsetdata* newfdtpnt,bool last);
		void pivot(direction d,lagindex<i>& prtmmit,itindex<i>& fmmind,vind vp,vind t,partialdata* newpdtpnt,subsetdata* newfdtpnt,bool last);
		vind			lastv;
		vind			p;
		vind			k;
		real			crt;
		deque<bool>		varin; 
		symtwodarray*		e;
		vector<matvectarray *>	ovct;
		rmgdata*		gdt;
};

}

#endif
