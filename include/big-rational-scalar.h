#include<rational.h>

class big_rational : public Scalar
{
	rational<bigint> r;
	
	public:
	
		big_rational():r(){}
		big_rational(int a):r(bigint(a)){}
		
		big_rational* plus_eq (Scalar* a) {big_rational* ptr = dynamic_cast<big_rational*>(a); r+=ptr->r; return this;}
		big_rational* minus_eq (Scalar* a) {big_rational* ptr = dynamic_cast<big_rational*>(a); r-=ptr->r; return this;}
		big_rational* times_eq (Scalar* a) {big_rational* ptr = dynamic_cast<big_rational*>(a); r*=ptr->r; return this;}
		big_rational* divide_eq (Scalar* a) {big_rational* ptr = dynamic_cast<big_rational*>(a); r/=ptr->r; return this;}
		big_rational* remainder_eq (Scalar* a) {r=bigint(0); return this;}
		bool eq (Scalar* a) const {big_rational* ptr = dynamic_cast<big_rational*>(a); return r == ptr->r;}
		bool ne (Scalar* a) const {big_rational* ptr = dynamic_cast<big_rational*>(a); return r != ptr->r;}
		bool gt (Scalar* a) const {big_rational* ptr = dynamic_cast<big_rational*>(a); return r > ptr->r;}
		bool lt (Scalar* a) const {big_rational* ptr = dynamic_cast<big_rational*>(a); return r < ptr->r;}
		bool ge (Scalar* a) const {big_rational* ptr = dynamic_cast<big_rational*>(a); return r >= ptr->r;}
		bool le (Scalar* a) const {big_rational* ptr = dynamic_cast<big_rational*>(a); return r <= ptr->r;}

		big_rational* abs_val() { big_rational* result = new big_rational(*this); result->r = abs(result->r); return result;}
		void increment() {}
		void decrement() {}
		void read(istream& s) {s >> r;}
		void print (ostream& s) {s << r;}
		int nl() {return num_len(r);}
		void dump(ostream& s) const;
		void sanitize() {}
		
		big_rational* clone() {return new big_rational(*this);}	
};

inline void big_rational::dump(ostream& s) const
{
	big_rational loc = *this;
	s << "numerator = "; 
	bigint& bn = loc.r.getn();	
	bn.dump(s); 
	s << " denominator = "; 
	bigint& bd = loc.r.getd(); 
	bd.dump(s);
	s << endl;
}

