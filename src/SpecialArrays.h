#ifndef ARRAYS
#define ARRAYS

#include <vector> 

using std::vector;
using leapsnbnds::vind;
using leapsnbnds::real;

class symtwodarray   {		//  Symmetric two dimensional array class
				//  Stores a symmetric matrix of doubles in compact form
	public:
		explicit symtwodarray(const vind dim);                  //  Constructor
		symtwodarray(const symtwodarray&);			//  Copy Constructor
		symtwodarray& operator=(const symtwodarray&);		//  Assignment Operator
		~symtwodarray(void);                                    //  Destructor 
		real&  operator() (const vind i,const vind j)		//  Subscripting: reading and writing  
			{ if (j<=i) return data[i][j];			//  Uses (,) sintaxe starting at (0,0) 
		      else return data[j][i];         }                     
		const real  operator() (const vind i,const vind j) const	   
			{ if (j<=i) return data[i][j];			//  Subscripting: read only  
		      else return data[j][i];         }                 //  Uses (,) sintaxe starting at (0,0)
	private:
		vind	dimension;  
		vector< vector<real> >	data;
};

void symatpivot(vind* const,const real,const symtwodarray&,symtwodarray&,const vind,const vind,const vind);
void vectorpivot(vind* const,const vector<real>&,vector<real>&,const symtwodarray&,const real,const vind,const vind,const vind);

class matvectarray  {		// One dimensional array whose elements can be stored in local memory
				// OR refer to the row of a symmetric matrix
	public:
		matvectarray(const vind,symtwodarray*,vind const); 
		void setvalue(const vind j,const real val);	// Writing  own data.
		const real  operator[] (const vind j) const;	// Subscripting: read only. Uses [] sintaxe starting at [0]  
		void switchtoowndata(void)     {mat = 0; }
	private:
		vind		dimension;  
		symtwodarray*	mat;
		vind		matrowind;
		vector<real>	owndata;			
		matvectarray(const matvectarray&);		//  Forbid copy construction
		matvectarray& operator=(const matvectarray&);	//  Forbid direct assignment
};

void vectorpivot(vind* const,const matvectarray&,matvectarray&,const symtwodarray&,const real,const vind,const vind,const vind);


#endif
