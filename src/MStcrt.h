#ifndef MSTATDATA
#define MSTATDATA

namespace extendedleaps {

class symtwodarray;  /* forward declaration     */

class partialwilksdata :  public partialdata {     /* Data used in Wilks statistic updates    */
	public:
		partialwilksdata(vind nv,real w)  : nvar(nv), wilksst(w)	{  }
		virtual ~partialwilksdata(void)					{  }
		const real			getepivot(void) const		{ return epivot; }
		void 				setepivot(real pv)		{ epivot = pv; }
		const real			gettpivot(void) const		{ return tpivot; }
		void 				settpivot(real pv)		{ tpivot = pv; }
		virtual const real	getcrt(void)   const			{ return wilksst; }
		virtual void setcrt(real w)					{ wilksst = w; }	
	private:
		vind			nvar;
		real			epivot;
		real			tpivot;
		real			wilksst;
	friend class wilksdata;
};

class wilksdata :  public subsetdata {
	public:
		wilksdata(vind nv,vind tnv,vind nvtopiv,vind hr,real wst);
		virtual ~wilksdata(void);
		virtual const real criterion(void) const { return wilksst; }
		virtual void setcriterion(real w)        { wilksst = w; }	
		virtual const real indice(void)	const; 
		virtual void  getpdata(partialdata *pd);  
		virtual real updatecrt(direction dir,mindices& mmind,vind var,partialdata* pdt) const;
		virtual void pivot(direction,mindices& mmind,vind vp,vind t,partialdata* pdt,subsetdata* fdt,bool last);
/* 
	Note: partialdata and subsetdata pointer must point to partialwilksdata and wilksdata classes
		  or unpredictable behaviour will result  
		  (general partialdata and subsetdata classes were used in order to garantee upward compability)
*/
		virtual subsetdata *crcopy(vind totalnv,vind partialnv)  const
			{  return new wilksdata(nvar,totalnv,partialnv,hrank,wilksst);  }
		virtual const real*	getbnds(void)	const	{ return 0; }	
		void setematcoef(vind i,vind j,real val)   { (*emat)(i,j) = val;  }
		void settmatcoef(vind i,vind j,real val)   { (*tmat)(i,j) = val;  }
		virtual void setorgvarl(vind *) {  }
	private:
		real updatecrt(direction dir,vind varind,partialdata* newdtpnt) const;   
		void pivot(lagindex<d>& prtmmit,vind vp,vind t,partialdata* newpdtpnt,subsetdata* newfdtpnt,bool last);
		void pivot(lagindex<i>& prtmmit,vind vp,vind t,partialdata* newpdtpnt,subsetdata* newfdtpnt,bool last);
		vind		nvar;
		vind		p;
		vind		k;
		vind		hrank;
		real		wilksst;
		symtwodarray*	emat;
		symtwodarray*	tmat;
};

class partialtracedata :  public partialdata {     /* Data used in trace statistic updates	*/
	public:
		partialtracedata(vind nvars,vind hrank);
		virtual ~partialtracedata(void)			{ delete pqf;  }
		virtual const real	getcrt(void)	const;
		partialsqfdata*  getpqfdata(void)	const	{ return pqf; }
	protected:
		vind			nvar;
		partialsqfdata*		pqf;
	friend class tracedata;
};

class tracedata :  public subsetdata {
	public:
		tracedata(vind nv,vind tnv,vind nvtopiv,vind hr,real crt);
		virtual ~tracedata(void) { delete sqf; }
		sqfdata*  getqfdata(void) const	{ return sqf; };
		virtual const real criterion(void) const;
		virtual void setcriterion(real c);	
		virtual void setorgvarl(vind *) {  }
		virtual void  getpdata(partialdata *);  
		virtual real updatecrt(direction dir,mindices& mmind,vind var,partialdata* pdt) const;
		virtual void pivot(direction dir,mindices& mmind,vind vp,vind t,partialdata* pdt,subsetdata* fdt,bool last);
/* 
	Note: partialdata and subsetdata pointer must point to partialtracedata and tracedata classes
		  or unpredictable behaviour will result  
		  (general partialdata and subsetdata classes were used in order to garantee upward compability)
*/
		virtual const real*	getbnds(void)	const	{ return 0; }	
	protected:
		vind		hrank;
		vind		nvar;
		sqfdata*	sqf;
};

class bartpistdata : public tracedata {
	public:
		bartpistdata(vind nv,vind tnv,vind nvtopiv,vind hr,real crt) 
			:  tracedata(nv,tnv,nvtopiv,hr,crt) {  }
		virtual ~bartpistdata(void) { }
		virtual const real indice(void)	const;
		virtual subsetdata *crcopy(vind totalnv,vind partialnv)  const
			{  return new bartpistdata(nvar,totalnv,partialnv,hrank,criterion());  }
};

class lawlhotstdata : public tracedata {
	public:
		lawlhotstdata(vind nv,vind tnv,vind nvtopiv,vind hr,real crt) 
			:  tracedata(nv,tnv,nvtopiv,hr,crt) {  }
		virtual ~lawlhotstdata(void) { }
		virtual const real indice(void)	const;
		virtual subsetdata *crcopy(vind totalnv,vind partialnv)  const
			{  return new lawlhotstdata(nvar,totalnv,partialnv,hrank,criterion());  }
};

}

#endif
