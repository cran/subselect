#ifndef VQFDATA
#define VQFDATA

using std::vector;
using leapsnbnds::vind;
using leapsnbnds::real;

class qfgdata;
class qfdata;

enum dir {forward,backward};

class vqfgdata : public globaldata {
	public:
		explicit vqfgdata(dir d) :  direction(d) { }
		dir getdirection(void)	    { return direction;  }
	private:
		dir			direction;
};

class vqfauxmem : public qfauxmem {
	public:
		explicit vqfauxmem(vind r) : qfauxmem(r) { tmpvc.reserve(r); }
		real*   gettmpvc(void)      { return &tmpvc[0]; }
	protected:
		vector<real>	tmpvc;  	
};


class vqfdata :  public qfdata {
	public:
		vqfdata(vind,vind,vind,vqfauxmem *,vqfgdata *,vind sbdimension,real criterion);
		virtual ~vqfdata(void);
		virtual real updatecrt(vind *,vind *,vind,vind) const;
		virtual void pivot(vind,vind,vind,vind,real,vind *,vind *,subsetdata *,bool);  
//		virtual void pivot(vind,vind,vind,vind,vind *,vind *,subsetdata *,bool);  
/*
Note: subsetdata pointer must point to vqfdata class or unpredictable behaviour will result (general subsetdata class was used in order to garantee upward compability)
 */
		virtual subsetdata *crcopy(vind totalnv,vind partialnv) const
			{  return new vqfdata(totalnv,partialnv,r,auxmem,gdt,dim,crt);  }
		virtual void setorgvarl(vind *) {  }
		virtual const real* getbnds(void)	const	{ return &vc[0]; }  
		vqfauxmem*	getauxmem(void) 	const	{ return auxmem; }
		vqfgdata*	getgdata(void) 		const	{ return gdt; }
	private:
		vind		dim; 
		vector<real>	vc;
		vqfauxmem*	auxmem;
		vqfgdata*       gdt;
};

#endif
