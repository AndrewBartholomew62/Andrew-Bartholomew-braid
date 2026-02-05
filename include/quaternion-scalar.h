/************************************************************************
			 Quaternian class
			 
		  A. Bartholomew December, 2005
		  
A concrete class for quaternion<scalar>

This file was renamed from Quaternion.h to quaternion-scalar.h to avoid a conflict with 
another header file quaternion.h experienced in some operating systems, 
such as Windows and MAC OS X.

We cannot simply typedef quaternion<scalar> to Quaternion since we need to 
re-define some of the operators

************************************************************************/
//#include <scalar.h>
#include <quaternion.h>


class Quaternion
{
		quaternion<scalar> n;
	
public:

	
	Quaternion(scalar t1=scalar(0), scalar t2=scalar(0), scalar t3= scalar(0), scalar t4= scalar(0)): n(t1,t2,t3,t4) {}
	Quaternion(quaternion<scalar> q):	n(q) {}
    explicit Quaternion (int i) : n(scalar(i),scalar(0),scalar(0),scalar(0)) {}

	scalar& Qr() {return n.Qr();}
	scalar& Qi() {return n.Qi();}
	scalar& Qj() {return n.Qj();}
	scalar& Qk() {return n.Qk();}

	scalar QrC()const {return n.QrC();}
	scalar QiC()const {return n.QiC();}
	scalar QjC()const {return n.QjC();}
	scalar QkC()const {return n.QkC();}

	bool real() const {return QiC() == scalar(0) && QjC() == scalar(0) && QkC() == scalar(0);}
	
	//default copy assignment used
    Quaternion operator += (const Quaternion q) {n += q.n; return *this;}
    Quaternion operator -= (const Quaternion q) {n -= q.n; return *this;}
	Quaternion operator *= (const Quaternion q) {n *= q.n; return *this;}
    Quaternion operator /= (const Quaternion q);

    bool operator == (const Quaternion q) const {return n == q.n;}
	
	/* NOTE: THE COMPARISON OPERATORS ARE FOR THE POLYNOMIAL OUTPUT OPERATOR << ONLY
	*/
    bool operator > (const Quaternion) const;
    bool operator < (const Quaternion) const;

	static void Quaternion_error()	{ cout << "\n\nError in Quaternions!"; exit(0);}
    friend istream& operator >> (istream& s, Quaternion& q) {s >> q.n; return s;}
    friend ostream& operator << (ostream& s, const Quaternion q) {s << q.n; return s;}
    friend int num_len (const Quaternion q) {return num_len(q.n);}
    friend Quaternion conj (const Quaternion q) {Quaternion r(q); r.n = conj(r.n); return r;}
	friend scalar norm_sqrd (const Quaternion q) {return norm_sqrd(q.n);}

	/* These are friends to allow implicit conversion of their arguments, helper functions do not support 
	   implicit conversion.  The / operator calls /= since /= is defined below, the others manipulate the 
	   result directly since it involves one les function call.
	*/
	friend Quaternion operator + (const Quaternion a, const Quaternion b) {Quaternion r(a); r.n += b.n; return r;}
	friend Quaternion operator - (const Quaternion a, const Quaternion b) {Quaternion r(a); r.n -= b.n; return r;}
	friend Quaternion operator * (const Quaternion a, const Quaternion b) {Quaternion r(a); r.n *= b.n; return r;}
	friend Quaternion operator / (const Quaternion a, const Quaternion b) {Quaternion r(a); r /= b; return r;}
	friend bool operator != (const Quaternion a, const Quaternion b) {return a.n != b.n;}

	/* NOTE: THIS % OPERATOR RETURNS O FOR ALL OPERANDS WITHUT CHECKING THAT NORM_SQURD !=0 FOR THE 2ND OPERAND */
	friend Quaternion operator % (const Quaternion, const Quaternion) {Quaternion r; return r;}
};

/* This operator needs defining again, since the quaternion operator / is a template that 
   uses T(1) in place of the explicit scalar(1) here. */
inline Quaternion Quaternion::operator /= (const Quaternion q)
{
	Quaternion invq = conj(q);
	invq.n *= scalar(1)/norm_sqrd(q);
	*this *= invq;
	return *this;
}

inline bool Quaternion::operator > (const Quaternion q) const
{
	if (n.QiC()==scalar(0) && n.QjC()==scalar(0) && n.QkC()==scalar(0) && 
	    q.n.QiC()==scalar(0) && q.n.QjC()==scalar(0) && q.n.QkC()==scalar(0)) 
		return n.QrC() > q.n.QrC();
	else 
		return true;
}

inline bool Quaternion::operator < (const Quaternion q) const
{
	if (n.QiC()==scalar(0) && n.QjC()==scalar(0) && n.QkC()==scalar(0) && 
	    q.n.QiC()==scalar(0) && q.n.QjC()==scalar(0) && q.n.QkC()==scalar(0)) 
		return n.QrC() < q.n.QrC();
	else 
		return false;
}
