#ifndef EMREALS
#define EMREALS

#include <limits>
#include <cmath>

using namespace std;

namespace  ErrMReals {

template <class T>
class errmonitreal {  // a real value of type T, with a first order estimate of its rounding error
	public:
		errmonitreal(const T v,const T e,const bool dpec) : val(v) , er(e) {  }	// Constructor	
		errmonitreal(const T v,const T e) : val(v) , er(e) {  }	// Constructor	
		errmonitreal(const T v) : val(v), er(0.)	{  }	// Constructor	
		errmonitreal(void)	:	val(0.), er(0.)	{  }	// Constructor
		errmonitreal(const errmonitreal<T>& org)
		:  val(org.value()) , er(org.err())  {   }	// Copy Constructor	
		errmonitreal& operator=(const errmonitreal& org);	// Assignment Operator
		~errmonitreal(void)  {   }	// Destructor 
		errmonitreal& operator+=(const errmonitreal& rhs);	// Add Assign
		errmonitreal& operator-=(const errmonitreal& rhs);	// Subtract Assign
		errmonitreal& operator*=(const errmonitreal& rhs);	// Multiply Assign
		errmonitreal& operator/=(const errmonitreal& rhs);	// Divide Assign
		operator T() const	{ return val; }	// Convert to basic type
		const	T	value(void)	const	{ return val; }	// Return value
		const	T	err(void)	const	{ return er; }	// Return estimate of roundig error		
		const static T RNDERR;  // Machine precision for basic type
		static bool dropec;  // Drop error checking when set to true
	private:
		T  val;	// computed value
		T  er;	// estimate of the relative rounding error
};

template <class T>
bool errmonitreal<T>::dropec = false; // Error checking active by default 

template<class T> const T errmonitreal<T>::RNDERR = std::numeric_limits<T>::epsilon(); 

template<class T> const errmonitreal<T> operator+(const errmonitreal<T>& lhs,const errmonitreal<T>& rhs);
template<class T> const errmonitreal<T> operator-(const errmonitreal<T>& lhs,const errmonitreal<T>& rhs);
template<class T> const errmonitreal<T> operator*(const errmonitreal<T>& lhs,const errmonitreal<T>& rhs);
template<class T> inline const errmonitreal<T> operator/(const errmonitreal<T>& lhs,const errmonitreal<T>& rhs);

template<class T> bool errcheck(const errmonitreal<T> *realp,const double tol);
template<class T> bool errcheck(errmonitreal<T> **realpl,const double tol,const int nreals);


template<class T> inline bool errcheck(const errmonitreal<T> *realp,const double tol)
{ 
	if (!errmonitreal<T>::dropec) return (realp->err() <= tol);
	else return true;
};

template<class T> inline bool errcheck(errmonitreal<T> **realpl,const double tol,const int nreals)
{ 
	if (!errmonitreal<T>::dropec) for (int i=0;i<nreals;i++) if (realpl[i]->err() > tol) return false;
	return true;
};

inline bool errcheck(double *realp,const double tol)
{ 
	return true;
};

inline bool errcheck(double **realpl,const double tol,const int nreals)
{ 
	return true;
};

template<class T> 
errmonitreal<T>& errmonitreal<T>::operator=(const errmonitreal<T>& org)
{	
	if (this != &org)  {
		val = org.value();
		if (!dropec) er = org.err();
	}
	return *this;
}

template<class T>
errmonitreal<T>& errmonitreal<T>::operator+=(const errmonitreal<T>& rhs)
{	
	T res = val + rhs.value();
	if (!dropec) {
		T absres = std::fabs(res);
		if (absres < std::numeric_limits<T>::min()) er = std::numeric_limits<T>::max();
		else er = (fabs(val)*er+fabs(rhs.value())*rhs.err())/absres+RNDERR;
	}
	val = res;
	return *this;
}

template<class T>
errmonitreal<T>& errmonitreal<T>::operator-=(const errmonitreal<T>& rhs)
{	
	T res = val - rhs.value();;
	if (!dropec) {
		T absres = std::fabs(res);
		if (absres < std::numeric_limits<T>::min()) er = std::numeric_limits<T>::max();
		else er = (fabs(val)*er+fabs(rhs.value())*rhs.err())/absres+RNDERR;
	}
	val = res;
	return *this;
}

template<class T>
errmonitreal<T>& errmonitreal<T>::operator*=(const errmonitreal<T>& rhs)
{	
	val *= rhs.value();
	if (!dropec) er += rhs.err() + RNDERR;
	return *this;
}

template<class T>
errmonitreal<T>& errmonitreal<T>::operator/=(const errmonitreal<T>& rhs)
{	
	val /= rhs.value();
	if (!dropec) er += rhs.err() + RNDERR;
	return *this;
}

template<class T>
const errmonitreal<T> operator+(const errmonitreal<T>& lhs,const errmonitreal<T>& rhs)
{	
	T res = lhs.value() + rhs.value();
	if (!errmonitreal<T>::dropec) {
		T er,absres=std::fabs(res);
		if (absres < std::numeric_limits<T>::min()) er = std::numeric_limits<T>::max();
		else er = (std::fabs(lhs.value())*lhs.err()+std::fabs(rhs.value())*rhs.err())/absres + errmonitreal<T>::RNDERR;
		return errmonitreal<T>(res,er);
	}
	else return errmonitreal<T>(res);
}

template<class T>
const errmonitreal<T> operator-(const errmonitreal<T>& lhs,const errmonitreal<T>& rhs)
{	
	T res = lhs.value()-rhs.value();
	if (!errmonitreal<T>::dropec) {
		T er,absres=std::fabs(res);
		if (absres < std::numeric_limits<T>::min()) er = std::numeric_limits<T>::max();
		else er = (std::fabs(lhs.value())*lhs.err()+std::fabs(rhs.value())*rhs.err())/absres + errmonitreal<T>::RNDERR;
		return errmonitreal<T>(res,er);
	}
	else return errmonitreal<T>(res);
}

template<class T>
const errmonitreal<T> operator*(const errmonitreal<T>& lhs,const errmonitreal<T>& rhs)
{	
	if (!errmonitreal<T>::dropec) return errmonitreal<T>(lhs.value()*rhs.value(),lhs.err()+rhs.err()+errmonitreal<T>::RNDERR);
	else return errmonitreal<T>(lhs.value()*rhs.value());
}

template<class T>
const errmonitreal<T> operator/(const errmonitreal<T>& lhs,const errmonitreal<T>& rhs)
{	
	if (!errmonitreal<T>::dropec) return errmonitreal<T>(lhs.value()/rhs.value(),lhs.err()+rhs.err()+errmonitreal<T>::RNDERR);
	else return errmonitreal<T>(lhs.value()/rhs.value());
}

}

#endif
