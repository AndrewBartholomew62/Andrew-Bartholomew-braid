/*********************************************************************** 
             bigint A. Bartholomew 2nd October 2005

    bigint is an arbitrary precision integer class with an interface designed 
	to meet the requirements of the braid programme.  It is based on earlier 
	int128 and int64 code but uses the division algorithm in Knuth The Art 
	of Computer Programing Vol 2 (third edition), page 272 (Algorithm D).  

	The class uses radix 2^16 storing the numerals in reverse order in a vector,
	so that the ith place stores the numeral for the ith power of the radix.  The 
	sign is stored separately.

***********************************************************************/

#include <string>
#include <sstream>

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>

using namespace std;

extern ofstream debug;

#include <util.h>
#include <bigint.h>

unsigned long decimal[10];

/********************* Function prototypes ********************/
void sum(bigint& result,const bigint& a, const bigint& b);
void diff(bigint& result,const bigint& a, const bigint& b);
int ten_power (int exponent);
int num_len (unsigned int n);
ostream& operator << (ostream& os, const bigint& a);

/*
bigint::bigint()
{
    n.resize(1);
	n[0] = 0;
    negative = false;
	cout << "\nmade a bigint with vector at " << &n << endl;
}

bigint::bigint() : n(1)
{
    negative = false;
}
*/


bigint::bigint(const int nn)
{

	int num = nn;
    if (num <0)
    {
		negative = true;
		num *= -1;
    }
    else
		negative = false;

	int digits = (num > BASE ? 2 : 1);

    n.resize(digits);
	n[0] = num%BASE;
	
	if (digits > 1)
		n[1] = num/BASE;
}

bigint::bigint(const unsigned long num)
{
	negative = false;

	int digits = (num > BASE ? 2 : 1);
    n.resize(digits);
	n[0] = num%BASE;
	
	if (digits > 1)
		n[1] = num/BASE;
}

bigint::operator bool() // conversion operator to a bool
{

if (bigint_control::DEBUG & bigint_control::bool_conv)
{
	debug << "bigint::operator bool() : ";
	(*this).dump(debug);
	debug << endl;
}

	if (*this != bigint(0))
		return true;
	else
		return false;
}

bigint::operator int() // conversion operator to an int
{
	if (negative)
		return -1 * n[0];
	else
		return n[0];
}

bigint& bigint::operator ++ () //prefix
{
    bigint local = 1;
    *this += local;
    return *this;
}

bigint bigint::operator ++ (int) //postfix
{
    bigint tmp = *this;
    ++*this;
    return tmp;
}

bigint& bigint::operator -- () //prefix
{
    bigint local = 1;
    *this -= local;
    return *this;
}

bigint bigint::operator -- (int) //postfix
{
    bigint tmp = *this;
    --*this;
    return tmp;
}

void bigint::dump(ostream& s) const
{
	
	s << "bigint(";
	s << n.size() << ',';
	
	if (negative)
		s << '-';
	else 
		s << '+';
		
	for (vector<unsigned int>::const_iterator nptr = n.begin(); nptr != n.end(); nptr++)
		s << ',' << *nptr;
	s << ')' << flush;
}

/* sanitize removes any leading zeros from the list of digits */

void bigint::sanitize()
{

if (bigint_control::DEBUG & bigint_control::sanitize)
{
	debug << "bigint::sanitize : ";
	(*this).dump(debug);
	debug << endl;
}
	size_t size = n.size();

if (bigint_control::DEBUG & bigint_control::sanitize)
{
	debug << "bigint::sanitize : initial size = " << size << endl;
}
	
	for (vector<unsigned int>::reverse_iterator n_ptr=n.rbegin(); n_ptr != n.rend(); n_ptr++)
	{
		if (*n_ptr == 0)
			size--;
		else
			break;	
	}

if (bigint_control::DEBUG & bigint_control::sanitize)
{
	debug << "bigint::sanitize : final size = " << size << endl;
}
	if(size == 0)
		size=1;
		
	n.resize(size); // discards high end zeros.

if (bigint_control::DEBUG & bigint_control::sanitize)
{
	debug << "bigint::sanitize : after sanitizing (*this) = ";
	(*this).dump(debug);
	debug << endl;
}

}

bigint operator + (const bigint& a, const bigint& b)
{
if (bigint_control::DEBUG & bigint_control::add)
{
	debug << "bigint::operator + : a = " << a << 
	         "\n                     b = " << b << endl;
}
    bigint result;

    if (a.negative && b.negative)
    {
if (bigint_control::DEBUG & bigint_control::add)
	debug << "bigint::operator + : both negative" << endl;
		result.negative = true;
		sum(result,a,b);
    }
    else if (a.negative)
    {
		if ( abs(a) >= abs(b) )
		{
if (bigint_control::DEBUG & bigint_control::add)
	debug << "bigint::operator + : a negative b positive result negative" << endl;
	    	result.negative = true;
	    	diff (result,a,b);
		}
		else
		{
if (bigint_control::DEBUG & bigint_control::add)
	debug << "bigint::operator + : a negative b positive result positive" << endl;
		    result.negative = false;
		    diff (result,b,a);
		}
    }
    else if (b.negative)
    {
		if ( abs(a) >= abs(b) )
		{
if (bigint_control::DEBUG & bigint_control::add)
	debug << "bigint::operator + : a positive b negative result positive" << endl;
	    	result.negative = false;
	    	diff (result,a,b);
		}
		else
		{
if (bigint_control::DEBUG & bigint_control::add)
	debug << "bigint::operator + : a positive b negative result negative" << endl;
	    	result.negative = true;
	    	diff (result,b,a);
		}
    }
    else
    {
if (bigint_control::DEBUG & bigint_control::add)
	debug << "bigint::operator + : both positive" << endl;
		result.negative = false;
		sum (result,a,b);
    }
    return result;
}

bigint operator - (const bigint& a, const bigint& b)
{
if (bigint_control::DEBUG & bigint_control::subtract)
	debug << "bigint::operator - :" << endl;
	
	bigint minusb = b*-1;
	return a+minusb;
}

bigint operator * (const bigint& a, const bigint& b)
{
	unsigned long carry;
    bigint result;
	
if (bigint_control::DEBUG & bigint_control::multiply)
{
	debug << "bigint::operator * : a = " << a << 
             "\n                     b = " << b << endl;
	debug << "bigint::operator * : a = ";
	a.dump(debug);
	debug << endl;
	debug << "bigint::operator * : b = ";
	b.dump(debug);
	debug << endl;
}

    /* first check for zero product */
    if (a == 0 || b == 0)
    {

if (bigint_control::DEBUG & bigint_control::multiply)
	debug << "bigint::operator * : zero operand, returning zero" << endl;	

		return 0;
    }

    if (a == 1 )
    {

if (bigint_control::DEBUG & bigint_control::multiply)
	debug << "bigint::operator * : a = 1, returning b" << endl;	

		return b;
    }

    if (b == 1 )
    {

if (bigint_control::DEBUG & bigint_control::multiply)
	debug << "bigint::operator * : b = 1, returning a" << endl;	

		return a;
    }

    if ((a.negative || b.negative) && !(a.negative && b.negative))
        result.negative = true;
    else
        result.negative = false;
        
	result.n.resize(a.n.size()+b.n.size());
if (bigint_control::DEBUG & bigint_control::multiply)
{
	debug << "bigint::operator * : a.n.size() = " << a.n.size() << " b.n.size() = " << b.n.size()  << " result.n.size() = " << result.n.size() << endl;
}
	
    /* now multiply b through by a */
    for (size_t apower=0; apower < a.n.size(); apower++)
    {
		unsigned long loc_res;
		
		if (a.n[apower] == 0)
	    	continue;

if (bigint_control::DEBUG & bigint_control::multiply)
{
	debug << "bigint::operator * : multiply b by " << a.n[apower] << endl;
}
		size_t rpower = apower;
		carry = 0;

		for (size_t bpower = 0; bpower < b.n.size(); bpower++)
		{
		    loc_res = (unsigned long)a.n[apower]*(unsigned long)b.n[bpower]+carry;

if (bigint_control::DEBUG & bigint_control::multiply)
{
	debug << "bigint::operator * : rpower = " << rpower << " carry = " << carry << endl;
	debug << "bigint::operator * : initial result = " << loc_res <<endl;
}			
	    	carry = 0;
	    	if ( loc_res > MAX_DIGIT )
	    	{
				carry = loc_res / BASE;
				loc_res %= BASE;

if (bigint_control::DEBUG & bigint_control::multiply)
{
	debug << "bigint::operator * :     carry = " << carry << endl;
	debug << "bigint::operator * :     interim result = " << loc_res <<endl;
}			
	    	}

			loc_res += (unsigned long)result.n[rpower];
				
if (bigint_control::DEBUG & bigint_control::multiply)
{
	debug << "bigint::operator * :   add in " << result.n[rpower] << endl;
	debug << "bigint::operator * :   interim result = " << loc_res <<endl;
}			
	    	if ( loc_res > MAX_DIGIT )
	    	{
				loc_res %= BASE;
				carry++;

if (bigint_control::DEBUG & bigint_control::multiply)
{
	debug << "bigint::operator * :     carry = " << carry << endl;
	debug << "bigint::operator * :     interim result = " << loc_res <<endl;
}			
	    	}
			
			result.n[rpower++] = (unsigned int)loc_res;

if (bigint_control::DEBUG & bigint_control::multiply)
{
	debug << "bigint::operator * :   result.n[" << rpower << "] = " << loc_res <<endl;
}			
		}

	    if (carry)
		{
			result.n[rpower] = carry;

if (bigint_control::DEBUG & bigint_control::multiply)
{
	debug << "bigint::operator * :   carry left over" << endl;
	debug << "bigint::operator * :     result.n[" << rpower << "] = " << loc_res <<endl;
}			
		}
    }

	if (result.n[result.n.size()-1] == 0)
		result.n.pop_back(); // removes the zero in the most significant digit

if (bigint_control::DEBUG & bigint_control::multiply)
{
	debug << "bigint::operator * : result = " << result << endl;
	debug << "bigint::operator * : result = ";
	result.dump(debug);
	debug << endl;
}
   return result;
}

/* the division operator is based on the algorithm D in Knuth The Art of Computer Programing 
   Vol 2 (third edition), page 272 (Algorithm D).  We use u,v, n, m i and j as in Knuth.
	   */
bigint operator / (const bigint& a, const bigint& b)
{
	bigint u = a;
	bigint v = b;
   	bigint q; //quotient

if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / : u = " << u << 
	         "\n                     v = " << v << endl;
	debug << "bigint::operator / : u = "; 
	u.dump(debug);
	debug << endl;
	debug << "bigint::operator / : v = "; 
	v.dump(debug);
	debug << endl;
}

   	if (v == bigint(0))
   	{
    	cout << "bigint::operator / : zero divide error" << endl;
       	exit(0);
   	}

   	if (v == bigint(1))
   	{
    	q = u;
		
if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / : v = 1, returning u" << endl;
}
       	return q;
   	}

   	if (u == v)
   	{
if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / : u = v, returning 1" << endl;
}
    	q = 1;
       	return q;
   	}

   	if ( abs(v) > abs(u) )
   	{

if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / : abs(v) > abs(u), returning 0" << endl;
}
       	q = 0;
       	return q;
   	}

   	if (u.n.size() == 1  && v.n.size() == 1)
   	{

if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / : both u and v fit in unsigned long" << endl;
}
       	unsigned long t = u.n[0]/v.n[0];
		q = t;

       	if ((u.negative || v.negative) && !(u.negative && v.negative))
	   		q.negative = true;

if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / : returning q = " << q << endl;
}
       	return q;
   	}
	else if (v.n.size() == 1)
	{

		unsigned long r = 0;
		size_t n = u.n.size();
		q.n.resize(n);

if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / : using single digit denominator algorithm" << endl;
	debug << "bigint::operator / : n = " << n << endl;
}

		for (size_t j=n; j> 0; j--)
		{
		
if (bigint_control::DEBUG & bigint_control::divide)
	debug <<"bigint::operator / :   j = " << j-1 << endl;			

			/* set q_{j} = (r*b + u_{j})/v_{0}
			           r = (r*b + u_{j}) - q_{j} * v_{0}
			   loop is for j=n-1 j>=0 j--, but j is unsigned, hence translation
			*/
			r = r * BASE + u.n[j-1];
if (bigint_control::DEBUG & bigint_control::divide)
	debug <<"bigint::operator / :     r = " << r << endl;			
			q.n[j-1] = r / v.n[0];
if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / :     q.n[j] = " << q.n[j-1] << endl;			
	debug << "bigint::operator / :     product = " << q.n[j-1] * v.n[0] << endl;
}
			r -= q.n[j-1] * v.n[0];
if (bigint_control::DEBUG & bigint_control::divide)
	debug << "bigint::operator / :     r = " << r << endl;			
		}

	   	/* now set the sign of the quotient */
	   	if ((u.negative || v.negative) && !(u.negative && v.negative))
    	    q.negative = true;

if (bigint_control::DEBUG & bigint_control::divide)
	debug <<"bigint::operator / :   q = " << q << endl;			
		
		q.sanitize();
		
if (bigint_control::DEBUG & bigint_control::divide)
	debug <<"bigint::operator / :   after sanitizing q = " << q << endl;			

		return q;
	}

if (bigint_control::DEBUG & bigint_control::divide)
	debug << "bigint::operator / : full division algorithm required" << endl;

   	/* so we've go some work to do, with abs(v) < abs(u), 
	*/
	size_t n = v.n.size();
	size_t m = u.n.size() - n;
	q.n.resize(m+1);
	
	/* normalize */
	unsigned long d = BASE/(v.n[n-1] +1);

if (bigint_control::DEBUG & bigint_control::divide)
	debug << "bigint::operator / : normalizing factor d = " << d << endl;

	if (d != 1)
	{
		u *= bigint(d);
		v *= bigint(d);
	}
	
if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / : scaled u = " << u << endl;
	debug << "bigint::operator / : scaled v = " << v << endl;
}

	/* provide the introduction of the new digit position from step D1
	   if it has not already been done by the normalizing */
	u.n.resize(m+n+1); 

if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / : resized u = ";
	u.dump(debug);
	debug << endl;
	debug << "bigint::operator / : v = ";
	v.dump(debug);
	debug << endl;}
	
	/* this is the loop through steps D2 to D7
	   again, loop is really j=m; j>=0; j-- but j is again unsigned.  As the complexity is greater here,
	   we use J as the loop control variable and adjust to j for use in the code.
	*/
	for (size_t J=m+1; J>0; J--)
	{
		size_t j= J-1;  // so we loop j=m; j>=0; j--

if (bigint_control::DEBUG & bigint_control::divide)
	debug << "bigint::operator / : j = " << j << endl;
		
		/* each digit is <= BASE-1, so u[j+n]*BASE +u[j+n-1] < (BASE-1)*BASE + BASE-1 = BASE^2 -1
		   and this is less than BASE^2, so it fits in an unsigned long.
		*/
		unsigned long r_hat = u.n[j+n]*BASE + u.n[j+n-1];
		unsigned long q_hat = r_hat/v.n[n-1];
		r_hat -= q_hat * v.n[n-1];
		
if (bigint_control::DEBUG & bigint_control::divide)
	debug << "bigint::operator / :   q_hat = " << q_hat << " r_hat = " << r_hat << endl;
		
		if (q_hat == BASE || q_hat * v.n[n-2] > BASE * r_hat + u.n[j+n-2])
		{
			q_hat--;
			r_hat += v.n[n-1];

if (bigint_control::DEBUG & bigint_control::divide)
	debug << "bigint::operator / :   after first test q_hat = " << q_hat << " r_hat = " << r_hat << endl;

			if (r_hat < BASE && (q_hat == BASE || q_hat * v.n[n-2] > BASE * r_hat + u.n[j+n-2]))
			{
if (bigint_control::DEBUG & bigint_control::divide)
	debug << "bigint::operator / :   after second test q_hat = " << q_hat << endl;
				q_hat--;
			}
		}
		
		/* multiply and subtract */
		unsigned long borrow = 0;
		bigint p = bigint(q_hat) * v; //product;		
		p.n.resize(n+1);

if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / :   q_hat * v = " << p << endl;
	debug << "bigint::operator / :   q_hat * v = ";
	p.dump(debug);
	debug << endl;
}
		
		for (size_t i = 0; i < n+1; i++)
		{

if (bigint_control::DEBUG & bigint_control::divide)
	debug << "bigint::operator / :   subtraction: i = " << i << endl;
		
			/* Here we would naturally write 
               u.n[j+i] -= p.n[i] + borrow;
			   if (u.n[j+i] < 0)
			   	...
			   
			   but the vector is unsigned integeters, so we have to 
			   carry out the check for borrow first.
			*/
if (bigint_control::DEBUG & bigint_control::divide)
	debug << "bigint::operator / :     p.n[" << i << "] = " << p.n[i] << " borrow = " << borrow << 
	" u.n[" << j+i << "] = " << u.n[j+i] << endl;

			if (p.n[i] + borrow > u.n[j+i])
			{
				u.n[j+i] += BASE;
if (bigint_control::DEBUG & bigint_control::divide)
	debug << "bigint::operator / :     add BASE to u.n[" << j+i << "] to give " << u.n[j+i] << endl;
				u.n[j+i] -= p.n[i] + borrow;
				borrow = 1;
			}
			else
			{
				u.n[j+i] -= p.n[i] + borrow;
				borrow = 0;
			}
if (bigint_control::DEBUG & bigint_control::divide)
	debug << "bigint::operator / :     u.n[" << j+i << "] = " << u.n[j+i] << " borrow = " << borrow << endl;
		}
			
		/* Step D6 */
		if (borrow)
		{
if (bigint_control::DEBUG & bigint_control::divide)
	debug << "bigint::operator / :     result negative, add back " << endl;
			q_hat--;
			
			unsigned long carry = 0;
			
			for (size_t i = 0; i < n+1; i++)
			{
				u.n[j+i] += v.n[i] + carry;
				if (u.n[j+i] >= BASE) 
				{
					u.n[j+i] -= BASE;
					carry = 1;
				}
				else
					carry = 0;
			}
			
if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / :   q_hat = " << q_hat << endl;
	debug << "bigint::operator / :   u = ";
	u.dump(debug);
	debug << endl;
}
			// ignore any carry left over
			
		}

		q.n[j] = q_hat;
		
	}

   	/* finally, set the sign of the quotient */
   	if ((u.negative || v.negative) && !(u.negative && v.negative))
        q.negative = true;

if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigint::operator / : result = " << q << endl;
	debug << "bigint::operator / : result = ";
	q.dump(debug);
	debug << endl;
}
	
	q.sanitize();
   	return q;
}

bigint operator % (const bigint& a, const bigint& b)
{
if (bigint_control::DEBUG & bigint_control::remainder)
	debug << "bigint::operator % : " << endl;
    return a - a/b * b;
}

bigint operator += (bigint& a, const bigint& b)
{
    return a = a+b;
}

bigint operator -= (bigint& a, const bigint& b)
{
    return a = a-b;
}

bigint operator *= (bigint& a, const bigint& b)
{
    return a = a*b;
}

bigint operator /= (bigint& a, const bigint& b)
{
    return a = a/b;
}

bigint operator %= (bigint& a, const bigint& b)
{
    return a = a%b;
}

bigint operator >>= (bigint& a, int n)
{
	vector<unsigned int>::iterator aptr;

if (bigint_control::DEBUG & bigint_control::r_shift)
{
	debug << "bigint::operator >>= : a = " << a  << ", n = " << n << endl;
	debug << "bigint::operator >>= : a = ";
	a.dump(debug);
	debug << endl;	
}
	 

//cout << "\nin shift operator" << endl;	

	for (int i=0; i < n; i++)
	{
		aptr = a.n.begin();	
		
//cout << "\n\niteration " << i+1 << " of " << n	<< endl;

		if (a.n.size() == 1)
		{
//cout << "\nsize = 1, *aptr = " << *aptr << endl;			
			*aptr >>= 1;
//cout << "\nsize = 1, after right shift *aptr = " << *aptr << endl;			
		}
		else
		{
			// r.n.size() >= 2
			for ( ;aptr != a.n.end()-1; aptr++)
			{
//cout << "\n*aptr = " << *aptr << endl;			
				*aptr >>= 1;
//cout << "\nafter right shift *aptr = " << *aptr << endl;			
				if (*(aptr+1) & 1)
					*aptr |= BIGINT_MSB;
			}
//cout << "\nafter loop *aptr = " << *aptr << endl;			
			*aptr >>= 1;
//cout << "\nafter right shift *aptr = " << *aptr << endl;			
		}
	}

	a.sanitize();
	
if (bigint_control::DEBUG & bigint_control::r_shift)
{
	debug << "bigint::operator >>= : returning a = " << a << endl;
	debug << "bigint::operator >>= : a = ";
	a.dump(debug);
	debug << endl;
}
	return a;
}

bigint operator <<= (bigint& a, int n)
{
	vector<unsigned int>::reverse_iterator aptr;
	
if (bigint_control::DEBUG & bigint_control::l_shift)
{
	debug << "bigint::operator <<= : a = " << a  << ", n = " << n << endl;
	debug << "bigint::operator <<= : a = ";
	a.dump(debug);
	debug << endl;	
}

	for (int i=0; i < n; i++)
	{
		aptr = a.n.rbegin();
		
 		if (a.n.size() == 1)
		{
//cout << "\naptr = " << *aptr << endl;		
			if (*aptr & BIGINT_MSB)
			{
				a.n.resize(a.n.size()+1);
				*a.n.rbegin() = 1;
//cout << "\nafter resizing aptr = " << *aptr << endl;		
			}
			*aptr <<= 1;
			*aptr &= BIGINT_BITS;
		}
		else
		{
//cout << "\naptr = " << *aptr << endl;		
			if (*aptr & BIGINT_MSB)
			{
				a.n.resize(a.n.size()+1);
				*a.n.rbegin() = 1;
//cout << "\nafter resizing aptr = " << *aptr << endl;		
			}
			// r.n.size() >= 2
			for ( ;aptr != a.n.rend()-1; aptr++)
			{
				*aptr <<= 1;
				*aptr &= BIGINT_BITS;
				if (*(aptr+1) & BIGINT_MSB)
					*aptr |= 1;
			}
			*aptr <<= 1;
			*aptr &= BIGINT_BITS;
		}
	}

if (bigint_control::DEBUG & bigint_control::l_shift)
{
	debug << "bigint::operator <<= : returning a = " << a << endl;
	debug << "bigint::operator <<= : a = ";
	a.dump(debug);
	debug << endl;
}
	return a;
}

bool operator == (const bigint& aa, const bigint& bb)
{
	bigint a = aa;
	bigint b = bb;
	a.sanitize();
	b.sanitize();
/* necessary? */

if (bigint_control::DEBUG & bigint_control::equal)
{
	debug << "bigint::operator == : a = " << a << 
	         "\n                      b = " << b << endl;
}

    if (a.n.size() != b.n.size())
	{
if (bigint_control::DEBUG & bigint_control::equal)
	debug << "bigint::operator == : different sizes, return false" << endl; 
		return false;
	}

    if (a.negative != b.negative)
	{
if (bigint_control::DEBUG & bigint_control::equal)
	debug << "bigint::operator == : different signs, return false" << endl; 
		return false;
	}
	
    for (size_t i = 0; i < a.n.size(); i++)
    {
		if (a.n[i] != b.n[i])
		{
if (bigint_control::DEBUG & bigint_control::equal)
	debug << "bigint::operator == : a.n[" << i << "] != b.n[" << i << "], return false" << endl; 
	    	return false;
		}
    }
if (bigint_control::DEBUG & bigint_control::equal)
	debug << "bigint::operator == : returning true" << endl; 
    return true;
}

bool operator != (const bigint& a, const bigint& b)
{
	return (!(a == b));
}

bool operator > (const bigint& a, const bigint& b)
{

if (bigint_control::DEBUG & bigint_control::greater)
{
	debug << "bigint::operator > : a = " << a << 
	         "\n                     b = " << b << endl;
	debug << "bigint::operator > : a = ";
	a.dump(debug);
	debug << endl;
	debug << "bigint::operator > : b = ";
	b.dump(debug);
	debug << endl;
}
	if  (a == b)
	{
if (bigint_control::DEBUG & bigint_control::greater)
	debug << "bigint::operator > : a = b, return false" << endl; 
		return false;
	}
		
	if (a.negative && !b.negative)
	{
if (bigint_control::DEBUG & bigint_control::greater)
	debug << "bigint::operator > : a negative b positive, return false" << endl; 
		return false;
	}
	
	if (!a.negative && b.negative)
	{
if (bigint_control::DEBUG & bigint_control::greater)
	debug << "bigint::operator > : a positive b negative, return true" << endl; 
		return true;
	}
	
	/* take local copies and sanitize to remove leading zeros */
	bigint aa = a;
	bigint bb = b;
	aa.sanitize();
	bb.sanitize();

	const bigint* A;
	const bigint* B;
	
	if (a.negative) // b.negative also, check |b| > |a|
	{
		A = &bb;
		B = &aa;
	}
	else // !a.negative && !b.negative, check |a| > |b|
	{
		A = &aa;
		B = &bb;
	}

if (bigint_control::DEBUG & bigint_control::greater)
{
	debug << "bigint::operator > : checking A > B where A = " << *A << 
	       "\n                                          B = " << *B << endl; 
}
	
	/* so here we return true iff |A| > |B| */
	
	if (A->n.size() > B->n.size())
	{
if (bigint_control::DEBUG & bigint_control::greater)
	debug << "bigint::operator > : A size > B size, return true" << endl; 
		return true;
	}
	else if (A->n.size() < B->n.size())
	{
if (bigint_control::DEBUG & bigint_control::greater)
	debug << "bigint::operator > : A size < B size, return false" << endl; 
		return false;
	}
	
	bool result=false;
    for (size_t i = A->n.size(); i > 0; i--)
    {
		if (A->n[i-1] > B->n[i-1])
		{
if (bigint_control::DEBUG & bigint_control::greater)
	debug << "bigint::operator > : A.n[" << i-1 << "] > B.n[" << i-1 <<"] return true" << endl; 
	    	result = true;
			break;
		}
		else if (A->n[i-1] < B->n[i-1])
		{
if (bigint_control::DEBUG & bigint_control::greater)
	debug << "bigint::operator > : A.n[" << i-1 << "] < B.n[" << i-1 <<"] return false" << endl; 
			break; // result already false
		}
    }

	return result;
}

bool operator < (const bigint& a, const bigint& b)
{
  if (a == b || a > b)
	  return false;
  else
  	return true;
}

bool operator >= (const bigint& a, const bigint& b)
{
	return (!(a<b));
}

bool operator <= (const bigint& a, const bigint& b)
{
	return (!(a>b));
}

/*(This version dumps the output, for debugging the real << operator
ostream& operator << (ostream& os, const bigint& a)
{
	a.dump(os);
	return os;
}
*/

/* This is the real << operator */

ostream& operator << (ostream& os, const bigint& a)
//ostream& output(ostream& os, const bigint& a)
{	
	unsigned int loc_debug = bigint_control::DEBUG;
	bigint_control::DEBUG = 0;

if (bigint_control::DEBUG & bigint_control::output)
	debug << "bigint::operator << : a = " << a << endl;

    if (a == bigint(0))
	{
	
if (bigint_control::DEBUG & bigint_control::output)
{
	debug << "bigint::operator << : a == 0 ";
	a.dump(debug);
	debug << endl;
}
	
		os << 0;
	}
	else
	{
		bigint loc = a;
		bigint ten(10);
		ostringstream oss;

		/* we handle the sign from a  later */
		loc.negative = false;

		do
		{

if (bigint_control::DEBUG & bigint_control::output)
{		
	debug << "bigint::operator << : loc = ";
	loc.dump(debug);
	debug << endl;
}
			int digit = loc % ten;

if (bigint_control::DEBUG & bigint_control::output)
{		
	debug << "bigint::operator << : digit = " << digit << endl;
}
			oss << digit;
								
			bigint t1 = loc/ten;

if (bigint_control::DEBUG & bigint_control::output)
{		
	debug << "bigint::operator << : t1 = ";
	t1.dump(debug);
	debug << endl;
}
			loc = t1;
			
//			loc /= ten;
		}while (loc);
		
		/* now the sign, since oss holds digits in reverse order */
		if (a.negative)
			oss << '-';
	
		string str = oss.str();
				
if (bigint_control::DEBUG & bigint_control::output)
	debug << "bigint::operator << : oss = " << oss.str() << endl;
		
		string::reverse_iterator sptr = str.rbegin();
		while (sptr != str.rend())
			os << *sptr++;	
	}
	
	bigint_control::DEBUG = loc_debug;
	return os;	
}

istream& operator >> (istream& is, bigint& a)
{
    string input;
	char ch = 0;
	bool negative = false;
	bigint ten(10);

    /* set a to zero */
    a = 0; // sets nagative = false

	is >> ch; // skips leading whitespace and reads first char

if (bigint_control::DEBUG & bigint_control::input)
	debug << "bigint::operator >> : initial character from istream: " << ch << endl;
	
	if (is.fail()) // || ch == 0) no non-white characters have been read
	{
if (bigint_control::DEBUG & bigint_control::input)
	debug << "bigint::operator >> : initial read failed, returning with failbit set on istream at " << &is << endl;	
		return is;
	}
	
	if (ch == '-' || ch == '+')
	{
		char next = 0;
		next = is.peek();
if (bigint_control::DEBUG & bigint_control::input)
	debug << "bigint::operator >> : read sign, peeked at next character and found '" << next << "'" << endl;
		if (!isdigit(next)) //includes case that there is no next character and peek failed or next
		{
			is.putback(ch);
			is.setstate(ios_base::failbit);
if (bigint_control::DEBUG & bigint_control::input)
	debug << "bigint::operator >> : spurious sign returned to istream, setting failbit on istream at " << &is << endl;
			return is;
		}
		else if (ch == '-')
		{
	   		negative = true;
if (bigint_control::DEBUG & bigint_control::input)
	debug << "bigint::operator >> : detected good minus sign" << endl;	
		}
		else
		{
if (bigint_control::DEBUG & bigint_control::input)
	debug << "bigint::operator >> : detected good plus sign" << endl;	
		}
		
	}
	else if (isdigit(ch))
	{
		is.putback(ch);
if (bigint_control::DEBUG & bigint_control::input)
	debug << "bigint::operator >> : putting back digit " << ch << " for main loop execution" << endl;
	}
	else
	{
		is.putback(ch);
		is.setstate(ios_base::failbit);
if (bigint_control::DEBUG & bigint_control::input)
	debug << "bigint::operator >> : no number to read, putting back the " << ch << " and setting failbit on istream at " << &is << endl;
		return is;
	}
	
	/* here we know there is at least one digit to read from the istream */
	do 
	{
		is.get(ch);

		if (is.fail())
		{
			// some characters have been read but we have come to the end of an istringstream or ifstream
			is.clear();  // the fail is local to bigint if we've come to the end of an istream

if (bigint_control::DEBUG & bigint_control::input)
	debug << "bigint::operator >> : read characters to end of istream, clearing state flags on istream at " << &is << endl;

			return is; // the fail will be visible to the calling code if we've just read '+' or '-'
		}		
		else if (isdigit(ch))
		{

if (bigint_control::DEBUG & bigint_control::input)
	debug << "bigint::operator >> : next character from istream: " << ch << endl;

			if (a)
				a *= ten;
				
			switch (ch)
			{
				case '1': a += bigint(1); break;
				case '2': a += (bigint)2; break;
				case '3': a += (bigint)3; break;
				case '4': a += (bigint)4; break;
				case '5': a += (bigint)5; break;
				case '6': a += (bigint)6; break;
				case '7': a += (bigint)7; break;
				case '8': a += (bigint)8; break;
				case '9': a += (bigint)9; break;
				default:;
			}

if (bigint_control::DEBUG & bigint_control::input)
	debug << "bigint::operator >> :   a = " << a << endl;		

		}
		else
		{
			if (!isspace(ch))
				is.putback(ch);			
			break; // reached end of digits
		}
	} while (true);

if (negative)	
	a.negative = true;

if (bigint_control::DEBUG & bigint_control::input)
	debug << "bigint::operator >> : final value of a = " << a << endl;		

    return is;
}

bigint abs (const bigint& a)
{
   bigint result = a;
   result.negative = false;
   return result;
}

/* sum calculates the sum of the bigints at a and b and stores it at
   result.  The signs of a and b must have been taken care of elsewhere, 
   sum asumes the negative flag in result is correct.
*/
void sum(bigint& result,const bigint& aa, const bigint& bb)
{

if (bigint_control::DEBUG & bigint_control::sum)
{
	debug << "bigint::sum : a = " << aa << 
	         "\n              b = " << bb << endl;
	debug << "bigint::sum : a = ";
	aa.dump(debug);
	debug << endl;
	debug << "bigint::sum : b = ";
	bb.dump(debug);
	debug << endl;
}

    /* check for a == 0 */
    if (aa == 0)
    {
		result = bb;
		return;
    }

    if (bb == 0)
	{
		result = aa;
		return;
	}

	bigint a = aa;
	bigint b = bb;

    /* set result to a and add b to it */
    bool hold_negative = result.negative;
    result = a;
    result.negative = hold_negative;

	size_t result_size = (a.n.size() > b.n.size() ? a.n.size() : b.n.size());
	result.n.resize(result_size);
	a.n.resize(result_size);
	b.n.resize(result_size);

if (bigint_control::DEBUG & bigint_control::sum)
{
	debug << "bigint::sum : original a.n.size() = " << aa.n.size() << " original b.n.size() = " << bb.n.size()  << endl;
	debug << "bigint::sum : a.n.size() = " << a.n.size() << " b.n.size() = " << b.n.size()  << " result.n.size() = " << result.n.size() << endl;
	debug << "bigint::sum : a now = " << a << " b now = " << b << endl;
}
	
    int carry = 0;
    for (size_t i=0; i<b.n.size(); i++)
    {
		long loc_res = (long)a.n[i] + (long)b.n[i] + (long)carry;

if (bigint_control::DEBUG & bigint_control::sum)
{
	debug << "bigint::sum :   adding digits " << i << " carry = " << carry << endl;
	debug << "bigint::sum :     interim result = " << loc_res;
}
		carry = 0;


		if ( loc_res > MAX_DIGIT )
		{
		    carry = 1;
	    	loc_res %= BASE;

if (bigint_control::DEBUG & bigint_control::sum)
{
	debug << "bigint::sum :       carry = " << carry << " interim result = " << loc_res;
}
		}
		result.n[i] = (unsigned int)loc_res;

if (bigint_control::DEBUG & bigint_control::sum)
{
	debug << "bigint::sum :     result.n[" << i << "] = " << loc_res <<endl;
}			

    }

    if (carry)
	{
		result.n.resize (result_size+1);
		result.n[result_size] = 1;

if (bigint_control::DEBUG & bigint_control::sum)
{
	debug << "bigint::sum :   carry left over, result resized to " << result.n.size() << endl;
}

	}
	
if (bigint_control::DEBUG & bigint_control::sum)
{
	debug << "bigint::sum : result = " << result << endl;
	debug << "bigint::sum : result = ";
	result.dump(debug);
	debug << endl;
}
	result.sanitize();
}

/* diff calculates the difference between the bigints at a and b and stores it at
   result.  The signs must have been taken care of elsewhere and
   a must be greater than or equal to b.

   diff clears the negative flag in result if the difference is 0
*/
void diff(bigint& result,const bigint& a, const bigint& bb)
{

if (bigint_control::DEBUG & bigint_control::diff)
{
	debug << "bigint::diff : a = " << a << 
	         "\n               b = " << bb << endl;
	debug << "bigint::diff : a = ";
	a.dump(debug);
	debug << endl;
	debug << "bigint::diff : b = ";
	bb.dump(debug);
	debug << endl;
}

    /* copy the numeric part of a into result; note a >= b */
	result.n = a.n;
	
	bigint b = bb;
	b.n.resize(a.n.size());
	

if (bigint_control::DEBUG & bigint_control::diff)
{
	debug << "bigint::diff : a.n.size() = " << a.n.size() << " bb.n.size() = " << bb.n.size()  << 
  	         " b.n.size() = " << b.n.size()  << " result.n.size() = " << result.n.size() << endl;
}


    /* now subtract b from a */
    int borrow = 0;
	size_t i;
    for (i = 0; i < b.n.size(); i++)
    {
		long loc_res = (long)result.n[i] - (long)b.n[i] - (long)borrow;

if (bigint_control::DEBUG & bigint_control::diff)
{
	debug << "bigint::diff :   subtracting digits " << i << " borrow = " << borrow << endl;
	debug << "bigint::diff :     interim result = " << loc_res << endl;
}

		borrow = 0;
		if ( loc_res < 0 )
		{
		    borrow = 1;
	    	loc_res += BASE;

if (bigint_control::DEBUG & bigint_control::diff)
{
	debug << "bigint::diff :       borrow = " << borrow << " interim result = " << loc_res << endl;
}

		}
		result.n[i] = (unsigned int)loc_res;

if (bigint_control::DEBUG & bigint_control::diff)
{
	debug << "bigint::diff :     result.n[" << i << "] = " << loc_res <<endl;
}			

    }
	
	if (borrow)
		throw bigint_error("diff attempting to subtract b from a where b>a");

    /* check for the zero result */
    bool loc_sign = result.negative;
    result.negative = false;
    if (!(result == bigint(0)))
		result.negative = loc_sign;

if (bigint_control::DEBUG & bigint_control::diff)
{
	debug << "bigint::diff : result = " << result << endl;
	debug << "bigint::diff : result = ";
	result.dump(debug);
	debug << endl;
}

	result.sanitize();

}

int num_len (unsigned int n)
{
    int i = 0;

    do
    {
        i++;
        n /= 10;
    } while(n != 0);

    return i;
}

int num_len(const bigint& a)
{
if (bigint_control::DEBUG & bigint_control::num_len)
	debug << "bigint::num_len : a = " << a << endl;
	
//	string str;
	ostringstream oss;//(str);
	oss << a;
//output(oss,a);
	
if (bigint_control::DEBUG & bigint_control::num_len)
	debug << "bigint::num_len : str = " << oss.str() << endl;
	
	return oss.str().length();
}

/* The following binary gcd algorithm is Algorithm B from Knuth, The Art of
   Computer Programming vol 2, section 4.5.2.  The notation used here is
   chosen to be consistent with Knuth's.
*/
//int hcf_count;
bigint gcd (const bigint& a, const bigint& b)
{
	int k = 0;
	bigint t;
	bigint u(a);
	bigint v(b);
	
	/* could use abs but faster to set value directly 
	   and avoid the function call
	*/
	u.negative = false;
	v.negative = false;

	if (u == bigint(0))
		return v;
	
	if (v == bigint(0))
		return u;

//	hcf_count = 0;
//  cout << endl;
	
if (bigint_control::DEBUG & bigint_control::gcd)
	debug << "bigint::gcd : u = " << u << " v = " << v << endl;

	while (u.even() && v.even()) //B.1
	{
		k++;
		u >>= 1;
		v >>= 1;
if (bigint_control::DEBUG & bigint_control::gcd)
{
	debug  << "bigint::gcd : B.1 both u and v even, k incremented to " << k << endl;
	debug << "bigint::gcd :   u >>= " << u << " v >>= " << v << endl;
}
	}
	
	if (u.even())               // B.2
	{

if (bigint_control::DEBUG & bigint_control::gcd)
	debug << "bigint::gcd : u remains even" << endl;
		t = u;
	}
	else
		t = v * bigint(-1);
	
if (bigint_control::DEBUG & bigint_control::gcd)
	debug << "bigint::gcd : B.2 t = " << t << endl;

	do
	{

if (bigint_control::DEBUG & bigint_control::gcd)
{
	debug << "bigint::gcd : t remains non-zero" << t << endl;
//	debug << "\n" << "\t" << t << flush;
}

		while (t.even())       // B.4
		{
			t >>= 1;           // B.3

if (bigint_control::DEBUG & bigint_control::gcd)
	debug << "bigint::gcd : B.4   t remains even t >> = " << t << endl;
	
		}
			
		if (t > bigint(0))     // B.5
		{
			u = t;

if (bigint_control::DEBUG & bigint_control::gcd)
	debug << "bigint::gcd : B.5   u set to " << u << endl;
	
		}
		else
		{
			v = t * bigint(-1);
if (bigint_control::DEBUG & bigint_control::gcd)
	debug << "bigint::gcd : B.5   v set to " << v << endl;
	
		}

		t = u-v;              // B.6

if (bigint_control::DEBUG & bigint_control::gcd)
	debug << "bigint::gcd : B.6   t set to " << t << endl;

	} while (t != bigint(0));
	
	u <<= k;

if (bigint_control::DEBUG & bigint_control::gcd)
	debug << "bigint::gcd : u right shifted " << k << " to give " << u << endl;
	
	return u;
}

