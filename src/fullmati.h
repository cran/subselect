#ifndef FMATRIXI
#define FMATRIXI

#include "matrixb.h"

class fullvct :  public vector {
	public:
		double  *data;
		unsigned *ri;
		unsigned *ir;
		fullvct(unsigned);
		virtual ~fullvct(void);            
		void rewind(void);
		virtual double v(unsigned);
		virtual void put(double,unsigned);
		virtual void swprows(unsigned,unsigned);
		virtual void asgrows(unsigned,unsigned,unsigned *);
		virtual unsigned getrow(unsigned);
		virtual unsigned getrowp(unsigned);
		virtual double maxabs(unsigned fele=1);
};

class fullmat :  virtual public wrkmatrix  {
	public:
		double   *data;
		unsigned *ri;
		unsigned *ir;
		unsigned *ci;
		unsigned *ic;
		fullmat(unsigned,unsigned);
		virtual ~fullmat(void);
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
		virtual void multcol(double,unsigned,unsigned fele=1);
		virtual void coloper(unsigned,double,unsigned,unsigned frow=1);
		virtual void multrow(double,unsigned,unsigned fele=1);
		virtual void rowoper(unsigned,double,unsigned,unsigned fcol=1);
		virtual void maxele(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1);
		virtual void minele(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1);
		virtual void maxabs(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1);
		virtual void minabs(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1);
};

#endif
