
class int_scalar : public Scalar
{
	int n; 
	
	public:
		int_scalar(){n=0;}
		int_scalar(int a){n=a;}
		
		int_scalar* plus_eq (Scalar* a) {int_scalar* ptr = dynamic_cast<int_scalar*>(a); n+=ptr->n; return this;}
		int_scalar* minus_eq (Scalar* a) {int_scalar* ptr = dynamic_cast<int_scalar*>(a); n=n-ptr->n; return this;}
		int_scalar* times_eq (Scalar* a) {int_scalar* ptr = dynamic_cast<int_scalar*>(a); n*=ptr->n; return this;}
		int_scalar* divide_eq (Scalar* a) {int_scalar* ptr = dynamic_cast<int_scalar*>(a); n/=ptr->n; return this;}
		int_scalar* remainder_eq (Scalar* a) {int_scalar* ptr = dynamic_cast<int_scalar*>(a); n%=ptr->n; return this;}
		bool eq (Scalar* a) const {int_scalar* ptr = dynamic_cast<int_scalar*>(a); return n == ptr->n;}
		bool ne (Scalar* a) const {int_scalar* ptr = dynamic_cast<int_scalar*>(a); return n != ptr->n;}
		bool gt (Scalar* a) const {int_scalar* ptr = dynamic_cast<int_scalar*>(a); return n > ptr->n;}
		bool lt (Scalar* a) const {int_scalar* ptr = dynamic_cast<int_scalar*>(a); return n < ptr->n;}
		bool ge (Scalar* a) const {int_scalar* ptr = dynamic_cast<int_scalar*>(a); return n >= ptr->n;}
		bool le (Scalar* a) const {int_scalar* ptr = dynamic_cast<int_scalar*>(a); return n <= ptr->n;}

		int_scalar* abs_val() { int_scalar* result = new int_scalar(*this); result->n = std::abs(result->n); return result;}
		void increment() {n++;}
		void decrement() {n--;}
		void read(istream& s) {s >> n;}
		void print (ostream& s) {s << n;}
		int nl() {return num_len(n);} 
		void dump(ostream& s) const {s << n;}
		void sanitize() {}
		
		int_scalar* clone() {return new int_scalar(*this);}	
};

