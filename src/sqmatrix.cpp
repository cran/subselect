#include "sqmatb.h"

sqmatrixb::sqmatrixb(unsigned dim)
 : matrixb(dim,dim), det(UNKNOWN)  {  }

sqmatrixb::~sqmatrixb()
{  }

opsqmat::opsqmat(unsigned dim)  
 :  matrix(dim,dim), det(UNKNOWN)  {  }
		
opsqmat::~opsqmat()  
{  }

wrksqmat::wrksqmat(unsigned dim) 
  :   opsqmat(dim) { }
		
wrksqmat::~wrksqmat()
{ }

factrz::factrz(unsigned dim)  
  : n(dim)       {  }
		
factrz::~factrz()        
{ }

sqmatrix::sqmatrix(unsigned dim) 
  :  sqmatrixb(dim)  
{ 
	inv = static_cast<opsqmat *>(NULL); 
	LUF = static_cast<factrz *>(NULL); 
}
		
sqmatrix::~sqmatrix()	  
{ 
}
		
void sqmatrix::cmpdet()         
{ 
	wrksqmat *m = cwrksqm(n); 
	cpmatdt(data,m);
	det = m->trngrp(); 
	delete m;
}

char sqmatrix::invert()
{
	if (det == 0) {
		errmsg("\nError: trying to invert singular matrix.\n");
		return FALSE;
	}
	wrksqmat *m = cwrksqm(n);
	cpmatdt(data,m);
	m->det = det;
	if (!m->invertrp()) {
		det = 0.;
		delete m;
		return FALSE;
	}
	inv = csqmdat(n);
	cpmatdt(m,inv);
	inv->det = m->det;
	if (det == UNKNOWN) det = 1./inv->det;
	return TRUE;
}

void sqmatrix::LUfact()        
{
	wrksqmat *m=cwrksqm(n); 
	cpmatdt(data,m);
	m->det=det;	
	det=m->LUFrp(); 
	LUF=cfactrz(n);
	LUF->getf(m); 
	delete m;
}

void  leftmult(vector *v,sqmatrix *m,vector *res)
{
	leftmult(v,m->data,res);
}

void  rightmult(sqmatrix *m,vector *v,vector *res)
{
	rightmult(m->data,v,res);
}

void  cpmatdt(sqmatrix *org,sqmatrix *res)
{
	cpmatdt(org->data,res->data);
}

void  cpmatdt(sqmatrix *org,matrix *res)
{
	cpmatdt(org->data,res);
}

void  cpmatdt(matrix *org,sqmatrix *res)
{
	cpmatdt(org,res->data);
}

void  scamult(double sca,sqmatrix *org,sqmatrix *res)
{
	scamult(sca,org->data,res->data);
}

void  addmat(sqmatrix *m1,sqmatrix *m2,sqmatrix *res)
{
	addmat(m1->data,m2->data,res->data);
}

void  submat(sqmatrix *m1,sqmatrix *m2,sqmatrix *res)
{
	submat(m1->data,m2->data,res->data);
}

void  multmat(sqmatrix *m1,sqmatrix *m2,sqmatrix *res)
{
	multmat(m1->data,m2->data,res->data);
}

void  multmat(matrix *m1,matrix *m2,sqmatrix *res)
{
	multmat(m1,m2,res->data);
}

void  addmtm(sqmatrix *mt,sqmatrix *m,sqmatrix *res)
{
	addmtm(mt->data,m->data,res->data);
}

void  submtm(sqmatrix *mt,sqmatrix *m,sqmatrix *res)
{
	submtm(mt->data,m->data,res->data);
}

void  multmtm(sqmatrix *mt,sqmatrix *m,sqmatrix *res)
{
	multmtm(mt->data,m->data,res->data);
}

void  multmtm(matrix *mt,matrix *m,sqmatrix *res)
{
	multmtm(mt,m,res->data);
}

void  addmmt(sqmatrix *m,sqmatrix *mt,sqmatrix *res)
{
	addmmt(m->data,mt->data,res->data);
}

void  submmt(sqmatrix *m,sqmatrix *mt,sqmatrix *res)
{
	submmt(m->data,mt->data,res->data);
}

void  multmmt(sqmatrix *m,sqmatrix *mt,sqmatrix *res)
{
	multmmt(m->data,mt->data,res->data);
}

void  multmmt(matrix *m,matrix *mt,sqmatrix *res)
{
	multmmt(m,mt,res->data);
}

void  addmtmt(sqmatrix *mt1,sqmatrix *mt2,sqmatrix *res)
{
	addmtmt(mt1->data,mt2->data,res->data);
}

void  submtmt(sqmatrix *mt1,sqmatrix *mt2,sqmatrix *res)
{
	submtmt(mt1->data,mt2->data,res->data);
}

void  multmtmt(sqmatrix *mt1,sqmatrix *mt2,sqmatrix *res)
{
	multmtmt(mt1->data,mt2->data,res->data);
}

void  multmtmt(matrix *mt1,matrix *mt2,sqmatrix *res)
{
	multmtmt(mt1,mt2,res->data);
}

