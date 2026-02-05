
class rational_scalar : public Scalar
{
	rational<int> r;
	
	public:
	
		rational_scalar():r(){}
		rational_scalar(int a):r(bigint(a)){}
		
		rational_scalar* plus_eq (Scalar* a) {rational_scalar* ptr = dynamic_cast<rational_scalar*>(a); r+=ptr->r; return this;}
		rational_scalar* minus_eq (Scalar* a) {rational_scalar* ptr = dynamic_cast<rational_scalar*>(a); r-=ptr->r; return this;}
		rational_scalar* times_eq (Scalar* a) {rational_scalar* ptr = dynamic_cast<rational_scalar*>(a); r*=ptr->r; return this;}
		rational_scalar* divide_eq (Scalar* a) {rational_scalar* ptr = dynamic_cast<rational_scalar*>(a); r/=ptr->r; return this;}
		rational_scalar* remainder_eq (Scalar* a) {r=0; return this;}
		bool eq (Scalar* a) const {rational_scalar* ptr = dynamic_cast<rational_scalar*>(a); return r == ptr->r;}
		bool ne (Scalar* a) const {rational_scalar* ptr = dynamic_cast<rational_scalar*>(a); return r != ptr->r;}
		bool gt (Scalar* a) const {rational_scalar* ptr = dynamic_cast<rational_scalar*>(a); return r > ptr->r;}
		bool lt (Scalar* a) const {rational_scalar* ptr = dynamic_cast<rational_scalar*>(a); return r < ptr->r;}
		bool ge (Scalar* a) const {rational_scalar* ptr = dynamic_cast<rational_scalar*>(a); return r >= ptr->r;}
		bool le (Scalar* a) const {rational_scalar* ptr = dynamic_cast<rational_scalar*>(a); return r <= ptr->r;}

		rational_scalar* abs_val() { rational_scalar* result = new rational_scalar(*this); result->r = abs(result->r); return result;}
		void increment() {}
		void decrement() {}
		void read(istream& s) {s >> r;}
		void print (ostream& s) {s << r;}
		int nl() {return num_len(r);}
		void dump(ostream& s) const;
		void sanitize() {}
		
		rational_scalar* clone() {return new rational_scalar(*this);}	
};

inline void rational_scalar::dump(ostream& s) const
{
	rational_scalar loc = *this;
	s << "numerator = " << loc.r.getn() << " denominator = " << loc.r.getd() << endl;
}

