
/* Scalar is an abstract class providing the interface to our scalar types */
class Scalar
{

public:
	virtual Scalar* plus_eq (Scalar*) = 0;
	virtual Scalar* minus_eq (Scalar*) = 0;
	virtual Scalar* times_eq (Scalar*) = 0;
	virtual Scalar* divide_eq (Scalar*) = 0;
	virtual Scalar* remainder_eq (Scalar*) = 0;
	virtual bool eq (Scalar*) const = 0;
	virtual bool ne (Scalar*) const = 0;
	virtual bool gt (Scalar*) const = 0;
	virtual bool lt (Scalar*) const = 0;
	virtual bool ge (Scalar*) const = 0;
	virtual bool le (Scalar*) const = 0;

	virtual Scalar* abs_val() = 0;
	virtual void increment() = 0;
	virtual void decrement() = 0;

	virtual void print(ostream&) = 0;
	virtual void read(istream&) = 0;
	virtual int nl() = 0; // num_len
	virtual void dump(ostream&) const = 0;
	virtual void sanitize() = 0;

	virtual Scalar* clone() = 0; // virtual form of copy constructor	
	virtual ~Scalar() {}
};

/* scalar is a concrete class that hides the polymorphic, 
   pointer-oriented operation of the underlying Scalar class.
*/
class scalar
{
	Scalar* C;
	static int variant;
public:
	enum scalar_variant {MOD_P, BIGINT, BIGRATIONAL};
	static void set_variant(scalar_variant t);
	static void show_variant(ostream& s);
	
	scalar operator = (const scalar& c) {if (this != &c){if (C) delete C; C = c.C->clone();} return *this;} 
	scalar operator ++ () {C->increment(); return *this;} // prefix
	scalar operator ++ (int) {scalar t = *this; C->increment(); return t;} // postfix
	scalar operator -- () {C->decrement(); return *this;}
	scalar operator -- (int) {scalar t = *this; C->decrement(); return t;}
	scalar operator += (const scalar& c) {C->plus_eq(c.C); return *this;} 
	scalar operator -= (const scalar& c) {C->minus_eq(c.C); return *this;} 
	scalar operator *= (const scalar& c) {C->times_eq(c.C); return *this;} 
	scalar operator /= (const scalar& c) {C->divide_eq(c.C); return *this;} 
	scalar operator %= (const scalar& c) {C->remainder_eq(c.C); return *this;} 
	bool operator == (const scalar& c) const {return C->eq(c.C);} 
	bool operator != (const scalar& c) const {return C->ne(c.C);} 
	bool operator > (const scalar& c) const {return C->gt(c.C);} 
	bool operator < (const scalar& c) const {return C->lt(c.C);} 
	bool operator >= (const scalar& c) const {return C->ge(c.C);} 
	bool operator <= (const scalar& c) const {return C->le(c.C);} 

	friend scalar operator + (const scalar& a, const scalar& b) {scalar r(a); r.C = r.C->plus_eq(b.C); return r;} 
	friend scalar operator - (const scalar& a, const scalar& b) {scalar r(a); r.C = r.C->minus_eq(b.C); return r;}
	friend scalar operator * (const scalar& a, const scalar& b) {scalar r(a); r.C = r.C->times_eq(b.C); return r;}
	friend scalar operator / (const scalar& a, const scalar& b) {scalar r(a); r.C = r.C->divide_eq(b.C); return r;}
	friend scalar operator % (const scalar& a, const scalar& b) {scalar r(a); r.C = r.C->remainder_eq(b.C); return r;}

	void dump(ostream& s) const {C->dump(s);}
	void sanitize() {C->sanitize();}
	friend scalar abs(const scalar& c) {scalar result; delete result.C; result.C = c.C->abs_val(); return result;}
	friend ostream& operator << (ostream& s, const scalar& c) {(c.C)->print(s); return s;}
	friend istream& operator >> (istream& s, const scalar& c) {(c.C)->read(s); return s;}
	friend int num_len (const scalar c) {return c.C->nl();}
	
	scalar();
	scalar(const scalar& c) {C = c.C->clone();}
	scalar(int);
	~scalar() {delete C; C=0;}
};

struct scalar_error {
	scalar_error (char* message) {cout << "\nscalar error!" << message << endl;}
};


/* scalar construction will be done via function pointers 
   that will point to a function returning a pointer to
   an object of the appropriate type
   
*/
extern Scalar* (*make_scalar) ();
extern Scalar* (*make_scalar_int) (int);

/* so now the scalar constructor can call make_scalar 
   to get the right type of object
*/
inline scalar::scalar()
{
	C=make_scalar();
}

inline scalar::scalar(int a)
{
	C=make_scalar_int(a);
}


/* finally, include the scalar classes that are based on Scalar */
#include<mod-p-scalar.h>
#include<bigint-scalar.h>
#include<big-rational-scalar.h>

