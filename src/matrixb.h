#ifndef MATRIXB
#define MATRIXB

const double    UNKNOWN =    9E99;
const unsigned  NWD =          11;


#ifndef TRUEFALSE
#define TRUE   1
#define FALSE  0
#define TRUEFALSE
#endif

#ifndef NULL

#define NULL ((void *) 0)

// class  {
//	public:
//		template<class T>
//		operator T*() const  { return 0;  }
//	private:
//		void operator&() const;
// }    NULL;

#define NULLPOINT

#endif  

#ifndef interfacef
void prmtend(const char *);
void errmsg(char *);
#define interfacef
#endif

class  vector  {
	public:
		unsigned n;
		vector(unsigned dim) : n(dim)  {  }
		virtual ~vector(void)          {  }  
		virtual  double v(unsigned)            = 0;
		virtual  void put(double,unsigned)     = 0;
};

class  matrixb  {
	public:
		unsigned n;
		unsigned p;
		matrixb(unsigned nrows,unsigned ncols) 
			: n(nrows), p(ncols) {  }
		virtual ~matrixb(void)   {  }
};

class  matrix : public matrixb {
	public:
		matrix(unsigned nrows,unsigned ncols) : matrixb(nrows,ncols)   {  }
		virtual  ~matrix(void)                                         {  } 
		virtual  double v(unsigned,unsigned)                             = 0;
		virtual  void put(double,unsigned,unsigned)                      = 0;
		virtual  void swprows(unsigned,unsigned)                         = 0;
		virtual  void swpcols(unsigned,unsigned)                         = 0;
		virtual  unsigned getrow(unsigned)                               = 0;
		virtual  unsigned getcol(unsigned)                               = 0;
		virtual  unsigned getrowp(unsigned)                              = 0;
		virtual  unsigned getcolp(unsigned)                              = 0;
};

void cpvectdt(vector *org,vector *res);
void scamult(double sca,vector *org,vector *res);
void addvect(vector *v1,vector *v2,vector *res);
void subvect(vector *v1,vector *v2,vector *res);
void leftmult(vector *v,matrix *m,vector *res);
void rightmult(matrix *m,vector *v,vector *res);
double vectprod(vector *v1,vector *v2);
void cpmatdt(matrix *org,matrix *res);
void scamult(double sca,matrix *org,matrix *res);
void transpose(matrix *m,matrix *res);
void addmat(matrix *m1,matrix *m2,matrix *res);
void submat(matrix *m1,matrix *m2,matrix *res);
void multmat(matrix *m1,matrix *m2,matrix *res);
void addmtm(matrix *mt,matrix *m,matrix *res);
void submtm(matrix *mt,matrix *m,matrix *res);
void multmtm(matrix *mt,matrix *m,matrix *res);
void addmmt(matrix *m,matrix *mt,matrix *res);
void submmt(matrix *m,matrix *mt,matrix *res);
void multmmt(matrix *m,matrix *mt,matrix *res);
void addmtmt(matrix *mt1,matrix *mt2,matrix *res);
void submtmt(matrix *mt1,matrix *mt2,matrix *res);
void multmtmt(matrix *mt1,matrix *mt2,matrix *res);

class  wrkmatrix : public matrix {
	public:
		wrkmatrix(unsigned nrows,unsigned ncols) : matrix(nrows,ncols) { }
		virtual ~wrkmatrix(void)                                       { }
		virtual  void multcol(double,unsigned,unsigned fele=1)           = 0;
		virtual  void coloper(unsigned,double,unsigned,unsigned frow=1)  = 0;
		virtual  void multrow(double,unsigned,unsigned fele=1)           = 0;
		virtual  void rowoper(unsigned,double,unsigned,unsigned fcol=1)  = 0;
		virtual  void maxele(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1) = 0;
		virtual  void minele(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1) = 0;
		virtual  void maxabs(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1) = 0;
		virtual  void minabs(unsigned *,unsigned *,double *,unsigned frow=1,unsigned fcol=1) = 0;
};

#endif
