#ifndef FSMATRIXI
#define FSMATRIXI

#include "sqmatb.h"

class  fwrksqm:     public wrksqmat  {
	protected:
		double   *data;
	public:
		unsigned *ri;
		unsigned *ir;
		unsigned *ci;
		unsigned *ic;
		fwrksqm(unsigned);
		virtual ~fwrksqm(void);
		void rewind(void);
		virtual double v(unsigned,unsigned);
		virtual void put(double,unsigned,unsigned);
		virtual void swprows(unsigned,unsigned);
		virtual void swpcols(unsigned,unsigned);
		virtual void asgrows(unsigned,unsigned,unsigned *);
		virtual void asgcols(unsigned,unsigned,unsigned *);
		virtual unsigned getrow(unsigned);
		virtual unsigned getcol(unsigned);
		virtual unsigned getrowp(unsigned);
		virtual unsigned getcolp(unsigned);
		virtual void multcol(double,unsigned,unsigned fele=1,unsigned lele=0);
		virtual void coloper(unsigned,double,unsigned,unsigned frow=1,unsigned lrow=0);
		virtual void multrow(double,unsigned,unsigned fele=1,unsigned lele=0);
		virtual void rowoper(unsigned,double,unsigned,unsigned fcol=1,unsigned lcol=0);
		virtual void maxele(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1);
		virtual void minele(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1);
		virtual void maxabs(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1);
		virtual void minabs(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1);
};

class ltmatrix :   public opsqmat  {
	private:
		double   *data;
	public:
		ltmatrix(unsigned);
		virtual ~ltmatrix(void);
		virtual double v(unsigned,unsigned);
		virtual void put(double,unsigned,unsigned);
		virtual unsigned getrow(unsigned);
		virtual unsigned getcol(unsigned);
		virtual unsigned getrowp(unsigned);
		virtual unsigned getcolp(unsigned);
	private:
		virtual void swprows(unsigned,unsigned)   { }
		virtual void swpcols(unsigned,unsigned)   { }
};

class utmatrix :  public ltmatrix  {
	public:
		utmatrix(unsigned);
		virtual ~utmatrix(void);
		virtual double v(unsigned,unsigned);
		virtual void put(double,unsigned,unsigned);
};

class symatrixb : public ltmatrix {
	public:
		unsigned    *vi;
		unsigned    *iv;
		symatrixb(unsigned);
		virtual ~symatrixb(void);  
		void rewind(void);
		virtual void swpvar(unsigned,unsigned);
		virtual void asgvar(unsigned,unsigned,unsigned *);
		virtual unsigned getvar(unsigned);
		virtual unsigned getvarp(unsigned);
		void cptowsqm(fwrksqm *);
		void cpfrwsqm(fwrksqm *);
		virtual double v(unsigned,unsigned);
		virtual void put(double,unsigned,unsigned);
		virtual void multcol(double,unsigned,unsigned fele=1,unsigned lele=0);
		virtual void coloper(unsigned,double,unsigned,unsigned frow=1,unsigned lrow=0);
		virtual void multrow(double,unsigned,unsigned fele=1,unsigned lele=0);
		virtual void rowoper(unsigned,double,unsigned,unsigned fcol=1,unsigned lcol=0);
		virtual void maxdele(unsigned *,double *,unsigned fvar=1);
		virtual void mindele(unsigned *,double *,unsigned fvar=1);
		virtual void maxdabs(unsigned *,double *,unsigned fvar=1);
		virtual void mindabs(unsigned *,double *,unsigned fvar=1);
		virtual void maxabs(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1);
		virtual void minabs(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1);
		char invertrp(void);
};

class fullmfact  :  public factrz  {
	public:
		ltmatrix  *L;
		utmatrix  *U;
		fullmfact(unsigned);
		virtual  ~fullmfact(void);
		virtual void getf(wrksqmat *);
		virtual double Lv(unsigned i,unsigned j)   {return L->v(i,j); }
		virtual double Uv(unsigned i,unsigned j)   {return U->v(i,j); }
};

class fullsqm :   public sqmatrix  {
	public:
		fullsqm(unsigned);
		virtual ~fullsqm(void);
		virtual opsqmat   *csqmdat(unsigned dim)    {return new fwrksqm(dim); }
		virtual wrksqmat  *cwrksqm(unsigned dim)    {return new fwrksqm(dim); }
		virtual factrz    *cfactrz(unsigned dim)    {return new fullmfact(dim); }
		double  v(unsigned i,unsigned j)            {return data->v(i,j); }
		void    put(double x,unsigned i,unsigned j) {data->put(x,i,j); }
};

class symatrix:   public sqmatrix {
	public:
		unsigned     rank;
		symatrixb    *data;
		symatrixb    *inv;
		vector       *egval;
		vector       **egvct;
		symatrix(unsigned);
		virtual ~symatrix(void);
		char invert(void);
		double cmpdet(void);
		virtual void diag(unsigned);
		virtual  opsqmat   *csqmdat(unsigned dim)   {return (data = new symatrixb(dim)); }
		virtual  wrksqmat  *cwrksqm(unsigned dim)   {return new fwrksqm(dim); }
		virtual  factrz    *cfactrz(unsigned dim)   {return new fullmfact(dim); }
		double  v(unsigned i,unsigned j)            {return data->v(i,j); }
		void    put(double x,unsigned i,unsigned j) {data->put(x,i,j); }
};

void cpmatdt(symatrixb *,symatrixb *);
void addmat(symatrixb *,symatrixb *,symatrixb *);
void submat(symatrixb *m1,symatrixb *m2,symatrixb *res);
void scamult(double sca,symatrixb *org,symatrixb *res);
void cpmatdt(symatrix *,symatrix *);
void addmat(symatrix *,symatrix *,symatrix *);
void submat(symatrix *m1,symatrix *m2,symatrix *res);
void scamult(double sca,symatrix *org,symatrix *res);


#endif
