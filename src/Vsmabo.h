#ifndef VSDAOB
#define VSDAOB

#include "ItIndices.h"

namespace extendedleaps {

enum direction {forward,backward};

class mindices  {
/*
                        Class with index pointers implementing the mapping from the current variable processing  sequence 
                        to their memory positions.  It includes four indices associated with the combinations of two 
                        memory indexing schemes and two acess methods.

                        Memory schemes:
	                   Full variable set - accessed through idfm and iifm 
	                   Pivoted variable set - accessed through idpm and iipm 
			Acess modes:
	                   Direct - idfm and idpm indices
			   Indirect - iifm and iipm indices
*/
	public:
		mindices(vind szfm,vind szpm,vind pmemlag);     
	/*
  		Constructor - direct acess in both memory indexing schemes
		Arguments:
			szfm - number of variables in full variable set
			szpm - number of variables in pivoted variable set
			pmemlag - current index of first variable in pivoted set
	*/
		mindices(vind szfm,vind szpm,vind pmemlag,vind* fmmlst);
	/*
  		Constructor - direct acess in pivoted set and indirect acess in full set
		Arguments:
			szfm - number of variables in full variable set
			szpm - number of variables in pivoted variable set
			pmemlag - current index of first variable in pivoted set
			fmlst - pointer to list with mapping of current variable sequence into full variable set
	*/
		mindices(vind szfm,vind szpm,vind pmemlag,vind* fmmlst,vind* pmmlst);
	/*
  		Constructor - direct acess in both memory indexing schemes
		Arguments:
			szfm - number of variables in full variable set
			szpm - number of variables in pivoted variable set
			pmemlag - current index of first variable in pivoted set
			fmlst - pointer to list with mapping of current variable sequence into full variable set
			pmlst - pointer to list with mapping of current variable sequence into pivoted variable set
	*/
		~mindices(void);    					/*  Destructor   */
		itindex<d>*	idfm(void)	{ return idfm_; }	
		lagindex<d>*	idpm(void)	{ return idpm_; }	
		itindex<i>*	iifm(void)	{ return iifm_; }
		lagindex<i>*	iipm(void)	{ return iipm_; }	
		void		asgnfmmiid(itindex<i>* i)	{ iifm_ = i; }
                  /* Assign list for indirect access  in full variable scheme    */
		void		asgnpmmiid(lagindex<i>* i)	{ iipm_ = i; } 
                  /* Assign list for i acess in pivoted variable scheme          */
		bool		direct(void)	{ return (iipm_ == 0); }
                   /* True if direct acces to pivoted variable set. False otherwise   */ 
	private:
		itindex<d>*	idfm_;		
		lagindex<d>*	idpm_;		
		itindex<i>*	iifm_;		
		lagindex<i>*	iipm_;		
};

class globaldata {
	public:
		virtual ~globaldata(void)   {  }
};

class partialdata {
	public:
		virtual ~partialdata(void)   {  }
		virtual const real	getcrt(void) const = 0;
};

class subsetdata {
/*
   General virtual class to store all necessary data related to a subset and update the comparison criterion
   Different versions of the algorithm are based on different specializations of this class
*/
	public:
		virtual ~subsetdata(void)  {  }               	/* Virtual destructor    */
		virtual const real criterion(void) const = 0; 	/* Returns comparison criterion   */
		virtual void setcriterion(real)	= 0;          	/* Sets comparison criterion      */
		virtual void  getpdata(partialdata* pd)  { setcriterion(pd->getcrt()); }
		virtual const real indice(void)    const { return criterion(); }
		   /* Returns comparison indice - a monotone function of the comparison
		      criterion, used for reporting the results                          */
		virtual real updatecrt(direction dir,mindices& mmind,vind var,partialdata* pdt) const = 0;
        /* 
		Updates and returns the comparison criterion 
		Arguments:
			d  -- direction of the update, set to: 
			  		"forward" when adding a new variable 
					"backward" when removing a current variable
                 	mmind -- mmindices object used to map the current variable sequence into their memory positions - see mindices definition
				Note: Two versions of updatecrt must be provided one, when mmind.direct() evalutes to true, using direct access 
				through mmind.idfm and mmin.idpm and the other, when mmind.direct() evalutes to false, using indirect access through 
				mmind.iifm and iipm. To avoid code duplication it is often convinient to employ function templates.
                	var -- current index of the pivot variable
			pdt -- pointer to a partialdata object where updated data, usefull to evaluate further subsets, can be stored
	*/
		virtual real updatecrt(direction dir,mindices& mmind,vind var,partialdata* pdt,real rqbound) const 
			{ return updatecrt(dir,mmind,var,pdt); }
	/* 
		Alternative version of updatecrt with additional argument:
			
			rqbond -- Required bound that new criterion has to satisfy for new subset to be worth considering. Usefull, 
				  when it is possible to find appropriate bounds, faster than the new criterion value.
				
		Should be redefined when, and only when, usebounds is redefined so that it evaluates to true.
	*/
		virtual bool usebounds(void)  { return false; }
		virtual void pivot(direction dir,mindices& mmind,vind vp,vind t,partialdata* pdt,subsetdata* fdt,bool last) = 0;
        /* 
		Updates the data necessary to compute the comparison criterion for future subsets
		Arguments:
			d  -- direction of the update, set to: 
			  		"forward" when adding a new variable 
					"backward" when removing a current variable
                 	mmind -- mmindices object used to map the current variable sequence into their memory positions - see mindices definition
				Note: Two versions of updatecrt must be provided one, when mmind.direct() evalutes to true, using direct access 
				through mmind.idfm and mmin.idpm and the other, when mmind.direct() evalutes to false, using indirect access through 
				mmind.iifm and iipm. To avoid code duplication it is often convinient to employ function templates.
			vp -- current index of the pivot variable
                 	t --  number of variables to be updated
			pdt -- pointer to a partialdata object containing the data already updated by	the updatecrt member function		          
			fdt -- pointer to a subsetdata object where further data, (not updated by by updatecrt), can be stored 
			last -- true if this subset data will not be used to evaluate further subsets. False otherwise. 
        */  
		virtual subsetdata *crcopy(vind totalnv,vind partialnv)  const = 0;
	/*    
		Creates a copy of current subset.
		Arguments:
			totalnv - Total number of variables
			partialnv - Number of variables that will serve as pivots in future updates derived from this subset
	*/  
		virtual void setorgvarl(vind* list) = 0;	/* Sets list of original variable indices to list   */
		virtual const real*	getbnds(void) const = 0;	
                   /* Returns pointer to list of criteria bounds for each subset size    */
};

class subset {
	public:
		subset(vind nvp,vind pnv,subsetdata *id,bool pdt,vind tnv); 
	/* 
		Cronstructor. Arguments:
			nvp - Number of variables pivoted by this subset parent search tree (in the original F&W algorithm nvp initially equals the total 
			      number of variables, and decreases by one after the first sorting rearrengment)
			pnv - Number of variables that can be pivoted from this subset
			id  - Pointer to data.
			pdt - True if data should be created and destroyed by subset constructores and destructores. False otherwise.
			tnv - Total number of variables
	*/
		subset(vind * const ivar,vind nvp,vind pnv,subsetdata *id,bool pdt,vind tnv);
	/* 
		Cronstructor. Aditional argument:
			ivar - Pointer of list of variable incides mapping their original memory positions into a specific processing order
        */
		~subset(void); 		/*  Destructor     */
		const vind getithvar(vind i)        	{ return orgvarind[i];	}
		void setithvar(vind ele,vind val)	{ orgvarind[ele] = val; }
		void reorder(vind *l);	/*  Assign processing order according to the list pointed by l   */
		void asgvar(vind fvar,vind nv,vind *lagv);
		subsetdata& getdata(void)		{ return *data;  }
		subsetdata *getdatap(void)		{ return data;  }
		vind getp(void)                     	{ return p; }
		void copyvar(subset &);
		vind *getvar(void)			{ return var; }
		void setvar(vind ele,vind val)		{ var[ele-1] = val;  }
		vind getvar(vind ele)			{ return var[ele-1]; }
		void setvarp(vind ele,vind val)		{ orgvarpos[ele] = val; }
		vind getvarp(vind ele)			{ return orgvarpos[ele];}
		void setnvar (vind n)			{ k = n; }
		vind getnvar(void)			{ return k; }
		void pivot(direction dir,vind vp,vind t,subset *newsp,bool last,bool usebnd,real acpbound);
		void sort(vind,vind);
		vind getpmemorypos(vind i)          { return pmemorypos[i]; }
		vind getfmemorypos(vind i)          { return fmemorypos[i]; }
	private:
		vind	p;		/*  Total number of variables    */
		vind	t;		/*  Number of variables that can be pivoted from this subset  */
		vind	k;		/*  Number of variables already pivoted in this subset        */
		vind*	var;		/*  Array with current indices of variables in subset         */
		vind	frstvarpm;	/*  Current index of first variable in pivoted set            */
		vind*	orgvarind;      /*  Array of original indices of the variables ordered by current indices  */
		vind*	orgvarpos;	/*  Array of current indices of the variables ordered by original  indices */
		vind*	pmemorypos;	/*  Array of pivoted variable memory positions ordered by current indices  */
		vind*	fmemorypos;	/*  Array of full varirable memory positions ordered by current indices    */
		mindices*  memii;
                     /*  Pointer to class of indice iteratores used to access variable memory positions */
		subsetdata*	data;		/*  Pointer to data    */
		bool		privatedata;
    /*  True if data should be created and destroyed by subset constructores and destructores. False otherwise. */
		void assgnmem(void);		/*  Auxiliary member function for memory allocation    */
};

typedef  subset* pkspc;

class wrkspace  {  /* General memory working space    */
	public:
		void initwrkspace(vind nv,subsetdata *data0,vind lstind,vind nvattop,vind nvatbot,vind* vattop,vind* vatbot);
		virtual ~wrkspace(void);
		subset&    subsetat(vind i)           { return *(wrklst[i-1]); } 
		bool	   usebounds(void)            { return usebounds_;}
		virtual void pivot(vind vp,vind t,vind li,vind lo,bool usebnd,real acpbound) = 0;
	private:
		vind	p;		/* Total number of variables                 */
		vind	nwl;    	/* Number of different memory positions      */
		pkspc   *wrklst;	/* Pointer to memory positions               */
		void	frontlsts(vind *l1,vind *l2,vind nl1,vind nl2,vind *ol);
	/*
		Switchs the order of the variable indices pointed by ol so that the indices pointed by l1 appear first, followed
		by the indices pointed by l2, and the remainning indices.
		Further arguments:
			nl1 - number of indices in the first list (pointed by l1)
			nl2 - number of indices in the second list (pointed by l2)
	*/
	protected:
		void pivotinit(subset *&,subset *&,vind,vind,vind);
	private:
		bool usebounds_;
};

class SRCwrkspace : public wrkspace  {	 /* Memory working space for forward searches    */
	public:
		SRCwrkspace(vind tp,vind nv,subsetdata *data0,vind* ivlst,vind* ovlst);
		virtual ~SRCwrkspace(void)     {   }
		virtual void pivot(vind vp,vind t,vind li,vind lo,bool usebnd,real acpbound);
};

class INVwrkspace : public wrkspace  {	/* Memory working space for backward searches   */
	public:
		INVwrkspace(vind tp,vind nv,subsetdata *data0,vind* ivlst,vind* ovlst);
		virtual ~INVwrkspace(void)     {   }
		virtual void pivot(vind vp,vind t,vind li,vind lo,bool usebnd,real acpbound);
};

extern vind *dmyv,*Flp;  
extern real *Fl;
extern SRCwrkspace* SW;
extern INVwrkspace* IW;

/*
            Declaration of main sscma routine 

( makes use of two pointers to subsetdata classes containing respectivelly the data relative to the null and the full variable subsets )

*/

bool sscma(bool,subsetdata *,subsetdata *);

}

#endif
