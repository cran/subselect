#include <vector>
#include "Sscma.h"
#include "SpecialArrays.h"


symtwodarray::symtwodarray(const vind dim) 
  :	dimension(dim)
{	
	data.assign(dim,vector<real>());
	for (vind i=0;i<dim;i++) data[i].reserve(i+1);
}

symtwodarray::~symtwodarray()   {  }

symtwodarray::symtwodarray(const symtwodarray& org ) 
:	dimension(org.dimension), data(org.data)   {  }

symtwodarray& symtwodarray::operator=(const symtwodarray& org ) 
{	
	if (this != &org)  {
		dimension = org.dimension;
		data = org.data;
	}
	return *this;
}

matvectarray::matvectarray(const vind dim,symtwodarray* m,vind const mr) 
: mat(m), matrowind(mr), dimension(dim)
{
	owndata.reserve(dimension);
}


void matvectarray::setvalue(const vind j,const real val)	        
{ 
	owndata[j] = val; 
}										 

const real matvectarray::operator[] (const vind j) const	   
{ 
	if (mat) return (*mat)(matrowind,j);                 
	else return owndata[j];                             
}

void symatpivot(vind* const order,const real pivotvalue,const symtwodarray& im,symtwodarray& om,
				const vind vp,const vind v1,const vind vl)
{
	vind pivotind,inrowi,incoli;
	real t,t1;

	if (order) pivotind = order[vp-1];
	else pivotind = vp-1;
	t = 1./pivotvalue;
	for (vind i=0;i<vl-v1+1;i++)  {
		if (order) inrowi = order[v1+i-1];
		else inrowi = v1+i-1;
		t1 = im(inrowi,pivotind) * t;
		for (vind j=0;j<=i;j++) {
			if (order) incoli = order[v1+j-1];
			else incoli = v1+j-1;
			om(i,j) = im(inrowi,incoli) - t1 * im(pivotind,incoli);
		}
	}
	#ifdef COUNTING 
	fpcnt += (vl-v1+1)*(vl-v1+4)/2;
	#endif
}

void vectorpivot(vind* const order,const std::vector<real>& iv,std::vector<real>& ov,const symtwodarray& im,
				 const real t1,const vind vp,const vind v1,const vind vl)
{
	vind pivotind,incoli;

	if (order) pivotind = order[vp-1];
	else pivotind = vp-1;
	for (vind j=0;j<vl-v1+1;j++) {
		if (order) incoli = order[v1+j-1];
		else incoli = v1+j-1;
		ov[j] = iv[incoli] - t1 * im(pivotind,incoli); 
	}
	#ifdef COUNTING 
	fpcnt += vl-v1+1;
	#endif
}

void vectorpivot(vind* const order,const matvectarray& iv,matvectarray& ov,const symtwodarray& im,const real t1,
				 const vind vp,const vind v1,const vind vl)
{
	vind pivotind,incoli;

	if (order) pivotind = order[vp-1];
	else pivotind = vp-1;
	for (vind j=0;j<vl-v1+1;j++) {
		if (order) incoli = order[v1+j-1];
		else incoli = v1+j-1;
		ov.setvalue(j,iv[incoli] - t1 * im(pivotind,incoli)); 
	}
	#ifdef COUNTING 
	fpcnt += vl-v1+1;
	#endif
}
