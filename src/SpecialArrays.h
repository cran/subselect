#ifndef ARRAYS
#define ARRAYS

#include <vector> 

using std::vector;

namespace extendedleaps {

class symtwodarray   {
/*  Symmetric two dimensional array class. Stores a symmetric matrix of real variables in compact form  */
	public:
		explicit symtwodarray(const vind dim);          	/*  Constructor  */
		symtwodarray(const symtwodarray&);			/*  Copy Constructor  */
		symtwodarray& operator=(const symtwodarray&);		/*  Assignment Operator  */
		~symtwodarray(void);                            	/*  Destructor  */
		real&  operator() (const vind i,const vind j)
/*  Subscripting: reading and writing.  Uses (,) sintaxe starting at (0,0)  */
		{ if (j<=i) return data[i][j];	else return data[j][i]; }
		const real  operator() (const vind i,const vind j) const
/*  Subscripting: read only. Uses (,) sintaxe starting at (0,0)  */
		{ if (j<=i) return data[i][j];	else return data[j][i]; }
	private:
		vind	dimension;  
		vector< vector<real> >	data;
};

template <class I>
void symatpivot(I&,const real,const symtwodarray&,symtwodarray&,const vind vp,const vind t);
template <class I>
void vectorpivot(I&,const vector<real>&,vector<real>&,const symtwodarray&,const real,const vind,const vind);

class matvectarray  {
/* One dimensional array whose elements can be stored in local memory OR refer to the row of a symmetric matrix  */
	public:
		matvectarray(const vind,symtwodarray*,vind const); 
		void setvalue(const vind j,const real val);		/* Writing  own data.  */
		const real  operator[] (const vind j) const;	
/* Subscripting: read only. Uses [] sintaxe starting at [0]  */
		void switchtoowndata(void)     {mat = 0; }
	private:
		vind		dimension;  
		symtwodarray*	mat;
		vind		matrowind;
		vector<real>	owndata;	
		matvectarray(const matvectarray&);		/*  Forbid copy construction  */
		matvectarray& operator=(const matvectarray&);	/*  Forbid direct assignment  */
};

template <class I>
void vectorpivot(I&,const matvectarray&,matvectarray&,const symtwodarray&,const real,const vind vp,const vind t);

template <class I>
void symatpivot(I& rowind,const real pivotvalue,const symtwodarray& im,symtwodarray& om,const vind vp,const vind t)
{
	I colind(rowind);
	vind pivotind=rowind[vp-1];
	real t0=1./pivotvalue,t1;

	rowind.reset(vp);
	for (vind i=0;i<t;rowind++,i++)  {
		t1 = im(rowind(),pivotind) * t0;
		colind.reset(vp);
		for (vind j=0;j<=i;colind++,j++) 
			om(i,j) = im(rowind(),colind()) - t1 * im(pivotind,colind());
	}
	#ifdef COUNTING 
	fpcnt += t*(t+3)/2;
	#endif
}

template <class I>
void vectorpivot(I& colind,const std::vector<real>& iv,std::vector<real>& ov,const symtwodarray& im,const real t1,const vind vp,const vind t)
{
	vind pivotind = colind[vp-1];      

	colind.reset(vp);
	for (vind j=0;j<t;colind++,j++)  
		ov[j] = iv[colind()] - t1 * im(pivotind,colind()); 
	#ifdef COUNTING 
	fpcnt += t;
	#endif
}

template <class I>
void vectorpivot(I& colind,const matvectarray& iv,matvectarray& ov,const symtwodarray& im,const real t1,const vind vp,const vind t)
{
	vind pivotind = colind[vp-1];

	colind.reset(vp);
	for (vind j=0;j<t;colind++,j++)  
		ov.setvalue(j,iv[colind()] - t1 * im(pivotind,colind())); 
	#ifdef COUNTING 
	fpcnt += t;
	#endif
}

}

#endif
