
class bigint_scalar : public Scalar
{
	bigint n;
	
	public:
	
		bigint_scalar():n(){}
		bigint_scalar(int a):n(bigint(a)){}
		
		bigint_scalar* plus_eq (Scalar* a) {bigint_scalar* ptr = dynamic_cast<bigint_scalar*>(a); n+=ptr->n; return this;}
		bigint_scalar* minus_eq (Scalar* a) {bigint_scalar* ptr = dynamic_cast<bigint_scalar*>(a); n-=ptr->n; return this;}
		bigint_scalar* times_eq (Scalar* a) {bigint_scalar* ptr = dynamic_cast<bigint_scalar*>(a); n*=ptr->n; return this;}
		bigint_scalar* divide_eq (Scalar* a) {bigint_scalar* ptr = dynamic_cast<bigint_scalar*>(a); n/=ptr->n; return this;}
		bigint_scalar* remainder_eq (Scalar* a) {bigint_scalar* ptr = dynamic_cast<bigint_scalar*>(a); n%=ptr->n; return this;}
		bool eq (Scalar* a) const {bigint_scalar* ptr = dynamic_cast<bigint_scalar*>(a); return n == ptr->n;}
		bool ne (Scalar* a) const {bigint_scalar* ptr = dynamic_cast<bigint_scalar*>(a); return n != ptr->n;}
		bool gt (Scalar* a) const {bigint_scalar* ptr = dynamic_cast<bigint_scalar*>(a); return n > ptr->n;}
		bool lt (Scalar* a) const {bigint_scalar* ptr = dynamic_cast<bigint_scalar*>(a); return n < ptr->n;}
		bool ge (Scalar* a) const {bigint_scalar* ptr = dynamic_cast<bigint_scalar*>(a); return n >= ptr->n;}
		bool le (Scalar* a) const {bigint_scalar* ptr = dynamic_cast<bigint_scalar*>(a); return n <= ptr->n;}

		bigint_scalar* abs_val() {bigint_scalar* result = new bigint_scalar(*this); result->n = abs(result->n); return result;}
		void increment() {n++;}
		void decrement() {n--;}
		void read(istream& s) {s >> n;}
		void print (ostream& s) {s << n;}
		int nl() {return num_len(n);}
		void dump(ostream& s) const {n.dump(s);}
		void sanitize() {n.sanitize();}
		
		bigint_scalar* clone() {return new bigint_scalar(*this);}	
};

