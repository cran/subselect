#ifndef VSDA
#define VSDA

template<class T> 
inline const T  SQR(const T x)  {return x*x;}
template<class T> 
inline const T  CUBE(const T x)  {return x*x*x;}
template<class T> 
inline const T  MIN_XY(const T x,const T y)  {return x < y ? x : y;}


const unsigned VARPL = 10;
const double   ZPRC  = 0.001;
const double   PPRC  = 0.001;
const double   YES   = 1E98;
const double   NO    = 1E99;
const unsigned SRC   = 1;
const unsigned INV   = 0;
const unsigned MIN   = 0;
const unsigned MAX   = 1;
const unsigned MS    = 5;
const unsigned C     = 1;
const unsigned RMAX  = 32767;
const double   DNO   = 1E99;


#ifndef TRUEFALSE
const unsigned TRUE  = 1;
const unsigned FALSE = 0;
#define TRUEFALSE
#endif


typedef  char*    pchar;


class sbset  {
	public:
		sbset(unsigned,unsigned);
		virtual ~sbset(void);
		virtual void save(int *,double *,unsigned);
		char*	   crtn(void)    { return crtn_;  }
		char*	   indn(void)    { return indn_;  }
		unsigned*  actvar(void)  { return actvar_; }
		unsigned   nvar(void)    { return nvar_;  }
		unsigned   r(void)       { return r_;     }
		double	   crt(void)     { return crt_;   }
	private:
		char*	   crtn_;
		char*      indn_;
		char       newset;
		unsigned   pos;
		unsigned   nvar_;
		unsigned   r_;
		unsigned*  actvar_;
		double	   crt_;
		sbset(const sbset&);
		sbset& operator=(const sbset&);
	friend sbset *csbset(unsigned,unsigned,unsigned *,char *,char *,double);
	friend void dsbset(sbset *);
};

class sbsetqe {
	public:
		unsigned  pos(void)				  { return pos_; }
		unsigned  assgnpos(unsigned p)    { return pos_ = p; }
		char      elegl(void)             { return elegl_; } 
		char      assgnelegl(char elg)    { return elegl_ = elg; }
		sbset*    set(void)               { return set_; }
		sbset*    assgnset(sbset *s)      { return set_ = s; }
		sbsetqe*  nxtset(void)            { return nxtset_; }
		sbsetqe*  assgnnset(sbsetqe *sp)  { return nxtset_ = sp; }
		sbsetqe*  prvset(void)            { return prvset_; }
		sbsetqe*  assgnpset(sbsetqe *sp)  { return prvset_ = sp; }
		sbsetqe(unsigned,sbset *);
		sbsetqe(unsigned);
	private:
		unsigned  pos_;
		char      elegl_;
		sbset*    set_;
		sbsetqe*  nxtset_;
		sbsetqe*  prvset_;
};

class subsetq {
	public:
		int      dir(void)	       { return dir_;  }
		unsigned maxele(void)      { return dir_;  }   
		unsigned nele(void)        { return nele_; }
		char     full(void)        { return full_; }
		char     all(void)         { return all_;  }
		sbsetqe* pmin(void)        { return pmin_;  }
		sbsetqe* pmax(void)        { return pmax_;  }
		subsetq(int,unsigned);
		~subsetq(void)             { delalle();  }
		double inserte(sbsetqe *);
		void   delel(sbsetqe *);
		void   delalle(void);
		void show(const char *,unsigned,int);
		virtual void save(int *,double *,unsigned,unsigned,int);
	private:
		const int      dir_;
		const unsigned maxele_;
		unsigned	   nele_;
		char		   full_;
		char		   all_;
		sbsetqe*       pmin_;
		sbsetqe*       pmax_;
		subsetq& operator=(const subsetq&);
};

typedef  subsetq*   psbsq;
typedef  sbsetqe*   psbstqe;
typedef  sbset*     psbst;

sbsetqe *csbsetqe(sbset *);
void dsbsetqe (sbsetqe *);
void asgmemory(void);
void crtwrksp(void);
void savfull(void);
void savfrst(void);
void fsort(void);
int  cmp(const void *,const void *);
int prcksp(unsigned,unsigned,unsigned,unsigned,unsigned);
void actvcnv(unsigned,unsigned,unsigned *,unsigned *);
char leap(unsigned,double,unsigned,unsigned);
void updtlst(void);
void freport(void);
void cleanup(void);

#ifndef interfacef
void prmtend(const char *);
void errmsg(char *);
#define interfacef
#endif


#endif
