/**************************************************
                  util.h

This header file is intended to hold common utilities not
related to any particular template class

**************************************************/

#include <sstream>
#include <vector>
#include <list>

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
   
       a1 = a = x1 a + y1 b; where x1=1 and y1=0
       b1 = b = x a + y b; where x=0 and y=1
   
   We start evaluating the gcd of a and b using Euclid's algorithm, so
   
       a1 = q1 b1 + r1
   so
       r1 = a1 - q1 b1
          = (x1 a + y1 b) - q1(x a +y b)
          = (x1 - q1 x)a + (y1 - q1 y)b

	   a2 = b1 = x a + y b
	   b2 = r1 = (x1 - q1 x)a + (y1 - q1 y)b

	   and so on.  If rn == 0, gcd is bn

       r3 = r1 - q3 r2
          = (x1 a + y1 b) - q3(x2 a +y2 b)
          = (x1 - q3 x2)a + (y1 - q3 y2)b
   
   and so on.  We arrive at x,y such that xa+yb=gcd(a,b)
*/
template <class T> T ext_gcd (T a, T b, T& x, T& y)
{
	T x1 = T(1);
	T y1 = T(0);
	x = T(0);
	y = T(1);
	
	for (;;)
	{
		T r = (a % b + b)%b; // we always want r to be positive
		T q = (a-r)/b; // calculates the quotient correctly for our choice of r

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
int least_common_multiple (list<int>& l);
int least_common_multiple (vector<int>& v);

// This function is used only for debugging .h files
void breakpoint();
