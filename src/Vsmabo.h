#ifndef VSDAOB
#define VSDAOB

using leapsnbnds::vind;
using leapsnbnds::real;


class globaldata {
	public:
		virtual ~globaldata(void)   {  }
};

class auxmemory {
	public:
		virtual ~auxmemory(void)   {  }
};

class subsetdata {
	public:
		virtual ~subsetdata(void)  {  }
		virtual real criterion(void) const = 0;
		virtual real indice(void)    const        { return criterion(); }
		virtual real updatecrt(vind *,vind *,vind,vind) const = 0;
		virtual void pivot(vind,vind,vind,vind,real,vind *,vind *,subsetdata *,bool) = 0;
		virtual subsetdata* crcopy(vind,vind) const = 0;
		virtual void setorgvarl(vind *) = 0;
		virtual const real*	getbnds(void) const = 0;	 
};

class subset {
	public:
		subset(vind,vind,subsetdata *,bool,vind);
		subset(vind * const,vind,vind,subsetdata *,bool,vind);
		~subset(void);
		const vind getithvar(vind i)        		{ return orgvarind[i];	}
		void setithvar(vind ele,vind val)		{ orgvarind[ele] = val; }
		void reorder(vind *);
		void asgvar(vind fvar,vind nv,vind *lagv);
		subsetdata& getdata(void)			{ return *data;  }
		subsetdata *getdatap(void)			{ return data;  }
		vind getp(void)                     		{ return p; }
		void copyvar(subset &);
		vind *getvar(void)				{ return var; }
		void setvar(vind ele,vind val)			{ var[ele-1] = val;  }
		vind getvar(vind ele)				{ return var[ele-1]; }
		void setvarp(vind ele,vind val)			{ orgvarpos[ele] = val; }
		vind getvarp(vind ele)				{ return orgvarpos[ele];}
		void setnvar (vind n)				{ k = n; }
		vind getnvar(void)				{ return k; }
		void pivot(vind vp,vind v1,vind vl,subset *newsp,bool last);
		void sort(vind,vind);
		vind getpmemorypos(vind i)          		{ return pmemorypos[i]; }
		vind getfmemorypos(vind i)          		{ return fmemorypos[i]; }
	private:
		vind		p;		//  Total number of variables 
		vind		k;              //  Number of variables in this subset
		vind*		var;            //  Array with current indices of variables in subset
		vind		frstvarpm;      //  Current index of first variable in partial memory 
		vind		lstvarpm;       //  Current index of last variable in partial memory
		vind*		orgvarind;	//  Array of original indices of the variables ordered by current indices
		vind*		orgvarpos;      //  Array of current indices of the variables ordered by original  indices
		vind*		pmemorypos;     //  Array of partial memory positions of variables ordered by current indices  
		vind*		fmemorypos;     //  Array of full memory positions of variables ordered by current indices 
		subsetdata*	data;		//  Pointer to data  
		bool		privatedata;	//  True if data should be created and destroyed by subset constructores and 
						//  destructores. False otherwise.
		void assgnmem(void);		//  Auxiliary member function for memory allocation
};

typedef  subset* pkspc;

class wrkspace  {
	public:
		wrkspace(bool,vind,vind,subsetdata *);
		~wrkspace(void);
		subset&    subsetat(vind i)           { return *(wrklst[i-1]); }  
		void pivot(vind,vind,vind,vind,vind);
	private:
		vind   p;
		vind   nwl;
		bool	   full;
		pkspc      *wrklst;
		void frontlsts(vind *,vind *,vind,vind,vind *);
};

extern vind *dmyv,*Flp;  
extern real *Fl;
extern wrkspace   *SW,*IW;

/*
            Declaration of main sscma routine 

     ( makes use of two pointers to subsetdata classes containing respectivelly
	   the data relative to the null and the full variable subsets                  )
*/

bool sscma(bool,subsetdata *,subsetdata *);


void symatpivot(real** const,real** const,const vind,const vind,const vind);
void fullvpivot(real* const ,real* const,real** const,const real,const vind,const vind,const vind);
void symatpivot(vind * const,real** const,real** const,const vind,const vind,const vind);
void fullvpivot(vind * const,real* const ,real* const,real** const,const real,const vind,const vind,const vind);

#endif
