#ifndef VSDAO
#define VSDAO

#include "lagmvct.h"

const extern double EPS;
const extern double NOTDF;
const extern double UNKNW;

const unsigned V      = 1;
const unsigned PTR    = 2;
const unsigned NS2SM1 = 3;

const unsigned VIN    =  0;
const unsigned VOUT   =  1;

extern double *Fl;
extern unsigned *Flp;

class cndvct : public lagvct  {
	public:
		kmatrix   *mt;
		unsigned   mr;
		cndvct(unsigned k,unsigned p,kmatrix *m,unsigned r): lagvct(k,p) {mt=m; mr=r;}
		double v(unsigned);
};

class kspace {
	public:
		unsigned   ind;
		unsigned   p;
		unsigned   k;
		unsigned   r;           
		unsigned   p0;
		unsigned   nvar;
		unsigned   *var;
		unsigned   *srtv;
		unsigned   *varp;
		unsigned   *vst;
		kmatrix    *e;
		fwrksqm    *s2m1;
		lagvct     **ve;
		cndvct     **ovct;
		cndvct     **ivct;
		double     v;
		double     c;
		kspace(unsigned,unsigned,unsigned,unsigned,unsigned);
		kspace(unsigned *,unsigned,unsigned,unsigned,unsigned,unsigned);
		~kspace(void);
		void mvtofront(unsigned *,unsigned);
		void sort(unsigned *,unsigned,unsigned);
		double Ftolv(unsigned i);
		void asgvar(unsigned,unsigned,unsigned,unsigned *,unsigned *);
	private:
		void assgnmem(void);
};

typedef  kspace* pkspc;

class wrkspace  {
	public:
		unsigned   ind;
		unsigned   p;
		unsigned   r;
		unsigned   nwl;
		double     *ve;
		double     *cndv;
		double     *tmpv;
		symatrix   *E;
		symatrix   *T;
		fullvct    **hv;
		pkspc      *wrklst;
		wrkspace(unsigned,unsigned,unsigned,unsigned,symatrix *,symatrix *,fullvct **);
		~wrkspace(void);
		void pivot(unsigned,unsigned,unsigned,unsigned,unsigned);
	private:
		void frontlsts(unsigned *,unsigned *,unsigned,unsigned,unsigned *);
		void pivot(fullvct *,fullvct *,kmatrix *,double,unsigned,unsigned,unsigned);
		void pivot(symatrixb *,symatrixb *,unsigned,unsigned,unsigned);
		void cpivot(symatrixb *,symatrixb *,unsigned,unsigned,unsigned,unsigned,unsigned *);
};

extern wrkspace   *SW,*IW;

#endif
