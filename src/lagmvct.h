#ifndef LAGMVCT
#define LAGMVCT

#include "fullmati.h"

const double BIG   = 1E+8;
const double NOTDF = 9E+99;
const double UNKNW = 9E+98;

class lagvct : public fullvct  {
	public:
		lagvct(unsigned,unsigned);
		virtual ~lagvct(void)               {   }
		unsigned lag(void)                  { return lag_; }
		double v(unsigned);
		void   put(double,unsigned);
	private:
		unsigned lag_;
};

class kmatrix : public symatrixb {
	public:
		kmatrix(unsigned,unsigned);
		~kmatrix(void);
		unsigned lag(void)                  { return lag_; }
		double   v(unsigned,unsigned);
		void     put(double,unsigned,unsigned);
		void 	 multcol(double,unsigned,unsigned fele=1,unsigned lele=0);
		void 	 rowoper(unsigned,double,unsigned,unsigned fcol=1,unsigned lcol=0);
		void 	 coloper(unsigned,double,unsigned,unsigned frow=1,unsigned lrow=0);
		void     siminvrp(void);
		double	 LUFrp(double *);
	private:
		unsigned lag_;
		unsigned *tmpvi;
		double   *old;
};

#endif
