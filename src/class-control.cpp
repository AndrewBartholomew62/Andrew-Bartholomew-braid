/****************************************************************************
 class-control contains the definition of all the statics required in the 
 various control structures used with the classes, and contains the scalar
 control functions, including setting the value of p for the mod_p scalar
 ****************************************************************************/

#include <string>
#include <sstream>

#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <iomanip>

using namespace std;

extern ofstream debug;

#include <util.h>

#include <initialize.h>

#ifdef INITIALIZE_SCALAR
#include <scalar.h>
#else
#ifdef INITIALIZE_BIGINT
#include <bigint.h>
#endif
#ifdef INITIALIZE_RATIONAL
#include <rational.h>
#endif
#endif

#ifdef INITIALIZE_SCALAR
	Scalar* make_bigint() {return new bigint;}
	Scalar* make_bigint(int a) {return new bigint(a);}
	Scalar* make_mod_p() {return new mod_p;}
	Scalar* make_mod_p(int a) {return new mod_p(a);}
	Scalar* make_big_rational() {return new big_rational;}
	Scalar* make_big_rational(int a) {return new big_rational(a);}
	
	/* here are the function pointers that will point to 
	   one of the above, set to the default variant
	*/
	Scalar* (*make_scalar) () = make_big_rational;
	Scalar* (*make_scalar_int) (int) = make_big_rational;
	
	
	/* The variant member of scalar is for debugging purposes only */
	int scalar::variant = scalar::BIGRATIONAL;
	
	/* We only need the if-else logic when we set  
	   the variant of the scalar, not when we create one.
	*/
	void scalar::set_variant(scalar_variant t)
	{
		if (t == BIGINT)
		{
			make_scalar = make_bigint;
			make_scalar_int = make_bigint;
			variant = BIGINT;
		}
		else if (t == BIGRATIONAL)
		{
			make_scalar = make_big_rational;
			make_scalar_int = make_big_rational;
			variant = BIGRATIONAL;
		}
		else
		{
			make_scalar = make_mod_p;
			make_scalar_int = make_mod_p;
			variant = MOD_P;
		}
	}
	
	
	void scalar::show_variant (ostream& s)
	{
		if (variant == BIGRATIONAL)
			s << "\nscalar::variant == BIGRATIONAL" << endl;
		else if (variant == BIGINT)
			s << "\nscalar::variant == BIGINT" << endl;
		else if (variant == MOD_P)
			s << "\nscalar::variant == MOD_P" << endl;
	}
#endif 

#ifdef INITIALIZE_BIGINT
	unsigned int bigint_control::DEBUG = 0;
#endif 

#ifdef INITIALIZE_RATIONAL
	bool rational_control::DEBUG = false;
#endif	


#ifdef INITIALIZE_MATRIX	
#include <matrix.h>
	unsigned int matrix_control::DEBUG = 0;
	bool matrix_control::COMFORT_DOTS = true;
	bool matrix_control::WAIT_INFO = false;
	bool matrix_control::SINGLE_LINE_OUTPUT = false;
	
	/* the number of rows into a determinant that are reported in debug.
	   if matrix_control::DET_DEBUG_LIMIT = m, then determinant debug 
	   is produced for sub-determinants of size >= m
	*/
	int matrix_control::DET_DEBUG_LIMIT = 0;
	int matrix_control::wait_threshold = 5;
	int matrix_control::wait_count = 0;
	int matrix_control::reset_count = 0;
#endif

#ifdef INITIALIZE_POLYNOMIAL	
#include <polynomial.h>
	unsigned int polynomial_control::DEBUG = 0;
	bool polynomial_control::WAIT_INFO = false;
	bool polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = true;
	bool polynomial_control::OUTPUT_PROXY_VARIABLES_ONLY = false;
	bool polynomial_control::WRITE_PARITY_PEER_CODES = false;
	bool polynomial_control::MOD_P = false;
	bool polynomial_control::TeX = false;
	int polynomial_control::wait_threshold = 10;
#endif


#ifdef INITIALIZE_QUATERNION
#include <quaternion.h>
	bool quaternion_control::DEBUG = false;
#endif


/* These are static definitions from each scalar type */
#ifdef INITIALIZE_MOD_P
	int mod_p::p = 5;
	vector<int> mod_p::inv;
	
	int mod_p_inverse(int n)
	{
		int p = mod_p::get_p();
		int x1 = 0;
		int x2 = 1;
		int y1 = 1;
		int y2 = 0;
		int a = n;
		int b = p;
		int x;
	
		while (b > 0)
		{
			int q = a/b;
			int r = a - q*b;
			x = x2 - q*x1;
			int y = y2 - q*y1;
			
			a = b;
			b = r;
			x2 = x1;
			x1 = x;
			y2 = y1;
			y1 = y;
		};
		
		
		return (x2%p + p)%p; // x2 may be large and negative
		
		return 1;
	}
	
	void mod_p::set_p (int i)
	{
		p = i;
		
		/* set up the inverses */
		inv.resize(i);
		inv[0] = 0; // just for safety
		inv[1] = 1;
		
		for (int j=2; j<i-1; j++)
			inv[j] = mod_p_inverse(j);
			
		inv[i-1] = i-1;
		
	}
#endif	
	
