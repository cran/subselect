#ifndef SQMATRIXB
#define SQMATRIXB

#include "matrixb.h"

class  sqmatrixb  : public matrixb  {
	public:
		double   det;
		sqmatrixb(unsigned);
		virtual ~sqmatrixb(void);
};

class  opsqmat :  public matrix  {
	public:
		double det;
		opsqmat(unsigned);
		virtual ~opsqmat(void);
};

class  wrksqmat:     public opsqmat {
	public:
		double det;
		wrksqmat(unsigned);
		virtual ~wrksqmat(void);
		virtual  void multcol(double,unsigned,unsigned fele=1,unsigned lele=0)               = 0;
		virtual  void coloper(unsigned,double,unsigned,unsigned frow=1,unsigned lrow=0)      = 0;
		virtual  void multrow(double,unsigned,unsigned fele=1,unsigned lele=0)               = 0;
		virtual  void rowoper(unsigned,double,unsigned,unsigned fcol=1,unsigned lcol=0)      = 0;
		virtual  void maxele(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1) = 0;
		virtual  void minele(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1) = 0;
		virtual  void maxabs(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1) = 0;
		virtual  void minabs(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1) = 0;
		char invertrp(void);
		double trngrp(void);
		double LUFrp(void);
};

class  factrz  {
	public:
		unsigned n;
		factrz(unsigned);
		virtual ~factrz(void);
		virtual void getf(wrksqmat *)         = 0;
		virtual double Lv(unsigned,unsigned)  = 0;
		virtual double Uv(unsigned,unsigned)  = 0;
};

class  sqmatrix  :  public sqmatrixb {
	public:
		opsqmat  *data;
		opsqmat  *inv;
		factrz   *LUF;
		sqmatrix(unsigned);
		virtual ~sqmatrix(void);
		void cmpdet(void);
		char invert(void);
		void LUfact(void);
		virtual  opsqmat     *csqmdat(unsigned) = 0;
		virtual  wrksqmat    *cwrksqm(unsigned) = 0;
		virtual  factrz      *cfactrz(unsigned) = 0;
};

void leftmult(vector *v,sqmatrix *m,vector *res);
void rightmult(sqmatrix *m,vector *v,vector *res);
void cpmatdt(sqmatrix *org,sqmatrix *res);
void cpmatdt(sqmatrix *org,matrix *res);
void cpmatdt(matrix *org,sqmatrix *res);
void scamult(double sca,sqmatrix *org,sqmatrix *res);
void addmat(sqmatrix *m1,sqmatrix *m2,sqmatrix *res);
void addmtm(sqmatrix *mt,sqmatrix *m,sqmatrix *res);
void addmmt(sqmatrix *m,sqmatrix *mt,sqmatrix *res);
void addmtmt(sqmatrix *mt1,sqmatrix *mt2,sqmatrix *res);
void submat(sqmatrix *m1,sqmatrix *m2,sqmatrix *res);
void submtm(sqmatrix *mt,sqmatrix *m,sqmatrix *res);
void submmt(sqmatrix *m,sqmatrix *mt,sqmatrix *res);
void submtmt(sqmatrix *mt1,sqmatrix *mt2,sqmatrix *res);
void multmat(sqmatrix *m1,sqmatrix *m2,sqmatrix *res);
void multmtm(sqmatrix *mt,sqmatrix *m,sqmatrix *res);
void multmmt(sqmatrix *m,sqmatrix *mt,sqmatrix *res);
void multmtmt(sqmatrix *mt1,sqmatrix *mt2,sqmatrix *res);
void multmat(matrix *m1,matrix *m2,sqmatrix *res);
void multmtm(matrix *mt,matrix *m,sqmatrix *res);
void multmmt(matrix *m,matrix *mt,sqmatrix *res);
void multmtmt(matrix *mt1,matrix *mt2,sqmatrix *res);

#endif
