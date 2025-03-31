#ifndef ITIND
#define ITIND

namespace extendedleaps {

enum accesstp {d,i};
template<accesstp tp> class itindex { };
template<accesstp tp> class lagindex { };


class indexbase {
	public:
		indexbase(vind n) : nele(n)			{ }
		virtual	~indexbase(void)			{  }
		virtual void	reset(vind i=0)			{ cur_ = i; }
		virtual const bool	more(void) const	{ return (cur_ < nele); }
		virtual void operator++(int)			{ cur_++; }
		virtual void operator+=(vind inc)		{ cur_ += inc; }
	protected:
		virtual const	vind	cur(void) const		{ return cur_; }
		vind	cur_;
		vind	nele;

};


class itindexd :  public indexbase {           /* Trivial d index type   */
	public:
		itindexd(vind n) : indexbase(n)			{ }
		virtual const vind operator()(void) const		{ return cur(); }
		virtual const vind	operator[](vind i) const	{ return i; }
};

class itindexi :  public indexbase {           /* Indirect index type    */
	public:
		itindexi(vind n,vind* il) : indexbase(n), indlist(il) { }
		itindexi(vind n,vector<vind>& il) : indexbase(n), indlist(&il[0]) { }
		virtual const vind	operator()(void) const		{ return indlist[cur()]; }
		virtual const vind	operator[](vind i) const	{ return indlist[i]; }
		virtual void		asglst(vind *lst)		{ indlist = lst; }
	protected:
		vind* indlist; 
};

class lagindexd : public itindexd  {  /* Lagged d index - implements an index offset   */
	public:
		lagindexd(vind n,vind lag) : itindexd(n)		{ lag_ = lag; }
		void setlag(vind lag)					{ lag_ = lag; }
		virtual void	reset(void)				{ cur_ = 0; }
		virtual void	reset(vind i)				{ cur_ = i-lag_; }
		virtual const vind	operator[](vind i) const	{ return i-lag_; }
	protected:
		vind lag_;
};

class lagindexi : public itindexi  {  /* Lagged i index - implements an index offset   */
	public:
		lagindexi(vind n,vind lag,vind* il) : itindexi(n,il)  { lag_ = lag; }	
		lagindexi(const vind n,const vind lag,vector<vind>& il) : itindexi(n,il) { lag_ = lag; }
		void setlag(vind lag)					                 { lag_ = lag; }
		virtual void	reset(void)				                 { cur_ = 0; }
		virtual void	reset(vind i)				                 { cur_ = i-lag_; }
		virtual const vind	operator[](vind i) const	                 { return indlist[i-lag_]; }
	protected:
		vind lag_;
};

}

#endif
