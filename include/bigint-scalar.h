/***********************************************************************
			   bigint-scalar

	       A. Bartholomew 4th December 2005

	This is a version of bigint that can be used with the dynamic
	coefficient scalar mechanism rather than stand-alone as a concrete 
	class.

***********************************************************************/

/* 
This class assumes a 32-bit architecture.  If a 64 bit architecture is used
in the future, the bigint constructors from const int and const unsigned long should be changed
or BASE redefined to accommodate the larger architecture.  Currently these constructors assume that an
unsigned long cannot be larger than BASE^2 -1
*/

#define BASE 65536
#define MAX_DIGIT 65535
#define BIGINT_MSB 0x8000
#define BIGINT_BITS 0xFFFF

#include <vector>

class bigint : public Scalar
{
    vector<unsigned int> n;
    bool negative;
	
public:
	
	/* these are the Scalar virtuals */
	bigint* plus_eq (Scalar* a) {bigint* ptr = dynamic_cast<bigint*>(a); *this += *ptr; return this;}
	bigint* minus_eq (Scalar* a) {bigint* ptr = dynamic_cast<bigint*>(a); *this -= *ptr; return this;}
	bigint* times_eq (Scalar* a) {bigint* ptr = dynamic_cast<bigint*>(a); *this *= *ptr; return this;}
	bigint* divide_eq (Scalar* a) {bigint* ptr = dynamic_cast<bigint*>(a); *this /= *ptr; return this;}
	bigint* remainder_eq (Scalar* a) {bigint* ptr = dynamic_cast<bigint*>(a); *this %= *ptr; return this;}
	bool eq (Scalar* a) const {bigint* ptr = dynamic_cast<bigint*>(a); return *this == *ptr;}
	bool ne (Scalar* a) const {bigint* ptr = dynamic_cast<bigint*>(a); return *this != *ptr;}
	bool gt (Scalar* a) const {bigint* ptr = dynamic_cast<bigint*>(a); return *this > *ptr;}
	bool lt (Scalar* a) const {bigint* ptr = dynamic_cast<bigint*>(a); return *this < *ptr;}
	bool ge (Scalar* a) const {bigint* ptr = dynamic_cast<bigint*>(a); return *this >= *ptr;}
	bool le (Scalar* a) const {bigint* ptr = dynamic_cast<bigint*>(a); return *this <= *ptr;}

	bigint* abs_val() { bigint* result = new bigint(*this); result->negative = false; return result;}
	void increment() {(*this)++;}
	void decrement() {(*this)--;}

	void print(ostream& s ) {s << *this;}
	void read(istream& s) {s >> *this;}
	int nl() {return num_len(*this);}
	bigint* clone() {return new bigint(*this);} // virtual form of copy constructor	
	void dump(ostream& s) const;
	void sanitize();

	/* These are the native bigint members */
    bigint():n(1) {negative = false;}
    bigint(const int num);
    bigint(const unsigned long num);
	
	operator int(); // conversion operator to an int
	operator bool(); // conversion operator to a bool

	bool even() const {return !(n[0]%2);}

    bigint& operator ++ (); //prefix
    bigint operator ++ (int); //postfix
    bigint& operator -- ();
    bigint operator -- (int);
	
	friend void sum(bigint& result,const bigint& a, const bigint& b);
	friend void diff(bigint& result,const bigint& a, const bigint& b);

    friend bigint operator += (bigint& a, const bigint& b);
    friend bigint operator -= (bigint& a, const bigint& b);
    friend bigint operator *= (bigint& a, const bigint& b);
    friend bigint operator /= (bigint& a, const bigint& b);
    friend bigint operator %= (bigint& a, const bigint& b);
	friend bigint operator >>= (bigint& a, int n);
	friend bigint operator <<= (bigint& a, int n);
    friend bigint operator + (const bigint& a, const bigint& b);
    friend bigint operator - (const bigint& a, const bigint& b);
    friend bigint operator * (const bigint& a, const bigint& b);
    friend bigint operator / (const bigint& a, const bigint& b);
    friend bigint operator % (const bigint& a, const bigint& b);
    friend bool operator == (const bigint& a, const bigint& b);
    friend bool operator != (const bigint& a, const bigint& b);
    friend bool operator > (const bigint& a, const bigint& b);
    friend bool operator < (const bigint& a, const bigint& b);
    friend bool operator <= (const bigint& a, const bigint& b);
    friend bool operator >= (const bigint& a, const bigint& b);
    friend ostream& operator << (ostream& os, const bigint& a);
    friend istream& operator >> (istream& is, bigint& a);
    friend bigint abs (const bigint& a);
    friend int num_len(const bigint& a);
	friend bigint gcd (const bigint& u, const bigint& v);
};

struct bigint_control
{
	static unsigned int DEBUG; //bitmap
	enum parameters
	{
		general = 	0x00000001, // unused
		sanitize =  0x00000002, 
		add = 		0x00000004,
		subtract = 	0x00000008,
		multiply = 	0x00000010,
		divide = 	0x00000020,
		remainder = 0x00000020,
		equal =  	0x00000040,
		greater =  	0x00000080,
	    output =  	0x00000100,
		input =  	0x00000200,
		sum =  		0x00000400,
		diff =  	0x00000800,
		num_len =  	0x00001000,
		r_shift =  	0x00002000,
		l_shift =  	0x00004000,
		bool_conv =	0x00008000,
		gcd =  		0x00010000,
		// 'all' also supported
	};
};

struct bigint_error {
	bigint_error (string message) {cout << "\nbigint error!" << message << endl;}
};

