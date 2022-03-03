#include<vector>

class mod_p : public Scalar
{
	static int p;

	int n; // modulo p
	
	public:
		static void set_p (int);
		static vector<int> inv;
		static int  get_p () {return p;}
	
		mod_p(){n=0;}
		mod_p(int a){n=(a%p+p)%p;}
		
		mod_p* plus_eq (Scalar* a) {mod_p* ptr = dynamic_cast<mod_p*>(a); n+=ptr->n; n %= p; return this;}
		mod_p* minus_eq (Scalar* a) {mod_p* ptr = dynamic_cast<mod_p*>(a); n=n+p-ptr->n; n %= p; return this;}
		mod_p* times_eq (Scalar* a) {mod_p* ptr = dynamic_cast<mod_p*>(a); n*=ptr->n; n %= p; return this;}
		mod_p* divide_eq (Scalar* a) {mod_p* ptr = dynamic_cast<mod_p*>(a); n*=inv[ptr->n]; n%=p; return this;}
		mod_p* remainder_eq (Scalar* a) {n=0; return this;}
		bool eq (Scalar* a) const {mod_p* ptr = dynamic_cast<mod_p*>(a); return n == ptr->n;}
		bool ne (Scalar* a) const {mod_p* ptr = dynamic_cast<mod_p*>(a); return n != ptr->n;}
		bool gt (Scalar* a) const {mod_p* ptr = dynamic_cast<mod_p*>(a); return n > ptr->n;}
		bool lt (Scalar* a) const {mod_p* ptr = dynamic_cast<mod_p*>(a); return n < ptr->n;}
		bool ge (Scalar* a) const {mod_p* ptr = dynamic_cast<mod_p*>(a); return n >= ptr->n;}
		bool le (Scalar* a) const {mod_p* ptr = dynamic_cast<mod_p*>(a); return n <= ptr->n;}

		mod_p* abs_val() { mod_p* result = new mod_p(*this); result->n = std::abs(result->n); return result;}
		void increment() { n++; n%= p;}
		void decrement() { n = n+p-1; n%= p;}
		void read(istream& s) {s >> n; n %= p; n = (n+p)%p;} // n = (n+p)%p deals with negative input
		void print (ostream& s) {s << n;}
		int nl(); // num_len
		void dump(ostream& s) const {s << n << "(mod " << p << ")";}
		void sanitize() {}
		
		mod_p* clone() {return new mod_p(*this);}	
};


inline int mod_p::nl ()
{
	int num = 1;
	int loc = n;
	
	while (loc /= 10)
		num++;
		
	return num;
}

