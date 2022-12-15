/**************************************************
                  util.h

This header file is intended to hold common utilities not
related to any particular template class

**************************************************/
#include <sstream>

template <class T> T gcd (T i, T j)
{
    if (i == T(0))
		return (j);
    else
		if (j == T(0))
		    return (i);
		else
	    	return ( gcd<T>(j, i%j));
}

/* ext_gcd implements the extended Euclidean algorithm to evaluate not only the gcd 
   but the linear combination of the original integers that yields the gcd.  Thus 
   given a and b, put
   
       r1 = a = x1 a + y1 b; where x1=1 and y1=0
       r2 = b = x2 a + y2 b; where x2=0 and y2=1
   
   We start evaluating the gcd of a and b using Euclid's algorithm, so
   
       r1 = q3 r2 + r3
   so
       r3 = r1 - q3 r2
          = (x1 a + y1 b) - q3(x2 a +y2 b)
          = (x1 - q3 x2)a + (y1 - q3 y2)b
   
   and so on.  We arrive at x,y such that xa+yb=gcd(a,b)
*/
template <class T> T ext_gcd (T a, T b, T& x, T& y)
{
	T x1 = T(1);
	T y1 = T(0);
//	T x = T(0);
//	T y = T(1);
	x = T(0);
	y = T(1);
	
	T q = a/b;
	T r = a % b;
	
	for (;;)
	{
		if (r == T(0))
		{
			return b;
		}
		else
		{
			a = b;
			b = r;
			
			T t = x1;
			x1 = x;
			x = t - q*x;
			
			t = y1;
			y1 = y;
			y = t -q*y;
		}
	}
}


template <class T> inline void get_number (T& num, char* cptr)
{
	istringstream ss(cptr);
	ss >> num;
}

template <class T> inline void get_number (T& num, string s, string::size_type pos)
{
	istringstream ss(s.substr(pos));
	ss >> num;
}

/* Other utility functions defined in util.cpp */
int num_len (long n);
char* c_string(const string& s);
void print_state(istream& is);
bool fileExists(const std::string& filename);
void substitute_string_variable(string& str, char ch1, char ch2);
char next_variable (string variables);

// This function is used only for debugging .h files
void breakpoint();

