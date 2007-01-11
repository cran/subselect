#ifndef VSQFDATA
#define VSQFDATA

using std::vector;

namespace extendedleaps {

class partialvsqfdata : public partialsqfdata    {
/* Data used in updates of sums of quadratic forms with a variable number of parcels  */
	public:
		explicit partialvsqfdata(vind nparcels,real vc0=0.);
		real*   gettmpvc(void)	{ return &tmpvc[0]; }
		virtual	~partialvsqfdata(void)  {  }	
	protected:
		vector<real>	tmpvc;
	friend class vsqfdata;
};

class vsqfdata :  public sqfdata {
/* Subset data for sums of quadratic forms with a variable number of parcels  */
	public:
		vsqfdata(vind tnv,vind nvtopiv,vind nparcels,real vc0,real sum);
		vsqfdata(vind tnv,vind nvtopiv,vind nparcels,const vector<real>& ovc,real sum);
		virtual ~vsqfdata(void);
		real*   getvc(void)	{ return &vc[0]; }
		void  setvc(real* x,vind nparcels);  
		void  setvc(real* x)   { setvc(x,r);  }  
		virtual real updatesum(direction d,mindices& mmind,vind var,vind dim,partialvsqfdata *pdt) const;
		virtual void pivot(direction d,mindices& mmind,vind vp,vind t,vind dim,partialvsqfdata* pdt,vsqfdata* fdt,bool last);
	private:
		real updatesum(direction d,vind varind,vind dim,partialvsqfdata* newdata) const;   
		template<accesstp tp> 
			void pivot(direction d,lagindex<tp>& prtmmit,vind vp,vind t,vind dim,partialvsqfdata* newpdata,vsqfdata* newfdata,bool last);
		vector<real>	vc;
};

}

#endif
