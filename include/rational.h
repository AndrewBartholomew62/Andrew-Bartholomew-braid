/************************************************************************
			 Rational definitions

		  A. Bartholomew  27th February, 2002

	   Modified to accept control.h mechanisms 22nd April 2002

************************************************************************/

struct rational_control {
	static bool DEBUG;
};


template <class T> class rational;

/* first define the friend functions */

template <class T> istream& operator >> (istream& is, rational<T>& Q)
{
	bool negative = false;
	char ch;

	is >> ch; // skips leading whitespace and reads first char

if (rational_control::DEBUG)
	debug << "rational::operator >> : initial character from istream: " << ch << endl;
	
	if (is.fail()) // || ch == 0) no non-white characters have been read
	{
if (rational_control::DEBUG)
	debug << "rational::operator >> : initial read failed, returning with failbit set on istream at " << &is << endl;	
		return is;
	}
	
	if (ch == '-' || ch == '+')
	{
		char next = 0;
		next = is.peek();
if (rational_control::DEBUG)
	debug << "rational::operator >> : read sign, peeked at next character and found '" << next << "'" << endl;
		if (!isdigit(next)) //includes case that there is no next character and peek failed or next
		{
			is.putback(ch);
			is.setstate(ios_base::failbit);
if (rational_control::DEBUG)
	debug << "rational::operator >> : spurious sign returned to istream, setting failbit on istream at " << &is << endl;
			return is;
		}
		else if (ch == '-')
		{
	   		negative = true;
if (rational_control::DEBUG)
	debug << "rational::operator >> : detected good minus sign" << endl;	
		}
		else
		{
if (rational_control::DEBUG)
	debug << "rational::operator >> : detected good plus sign" << endl;	
		}
		
	}
	else if (isdigit(ch))
	{
		is.putback(ch);
if (rational_control::DEBUG)
	debug << "rational::operator >> : putting back digit " << ch << " for reading numerator" << endl;
	}
	else
	{
		is.putback(ch);
		is.setstate(ios_base::failbit);
if (rational_control::DEBUG)
	debug << "rational::operator >> : no number to read, putting back the " << ch << " and setting failbit on istream at " << &is << endl;
		return is;
	}
	
	/* here we know there is at least one digit to read from the istream */
	is >> Q.n;

if (rational_control::DEBUG)
	debug << "rational::operator >> : absolute value of numerator = " << Q.n  << endl;
	
	if (negative)
		 Q.n *= T(-1);

if (rational_control::DEBUG)
	debug << "rational::operator >> : final value of numerator = " << Q.n  << endl;

	ch = is.get();
	
	if (is.fail())
	{
		is.clear();  // the fail is local to bigint if we've come to the end of an istream
		Q.d = T(1);

if (rational_control::DEBUG)
{
	debug << "rational::operator >> : read characters to end of istream, clearing state flags on istream at " << &is << endl;
	debug << "rational::operator >> : read characters to end of istream, setting denominator = 1" << endl;
}
		return is;
		
	}
	 
	if (ch != '/')
	{		

if (rational_control::DEBUG)
	debug << "rational::operator >> : read character other than '/', returning it to istream and setting denominator = 1" << endl;

		is.putback(ch);
		Q.d = T(1);
	}
	else
	{
if (rational_control::DEBUG)
	debug << "rational::operator >> : read '/' , proceeding to read denominator" << endl;
		is >> Q.d;
if (rational_control::DEBUG)
	debug << "rational::operator >> : denominator = " << Q.d  << endl;
	}
	
	return is;
}

template <class T> inline ostream& operator << (ostream& s, rational<T> Q)
{
    if (Q.d == T(1))
		s << Q.n;
    else
    {
	    s << Q.n << '/' << Q.d;
    }
    return s;
}

template <class T> int num_len (const rational<T>& Q)
{
    int i = 1;

    if (Q.d == 1)
		return num_len(Q.n);
    else
    {
		i += num_len(Q.n);
		i += num_len(Q.d);
		return i;
    }
}

template <class T> inline rational<T> scaled_rational(rational<T>& Q)
{
	T g;
	
    if ( Q.d != T(1) )
    {
		g = gcd ( abs(Q.n), Q.d );
		Q.n /= g;
		Q.d /= g;
    }
	
	return Q;
}

template<class T> inline rational<T> operator + (const rational<T> a, const rational<T> b)
{
    rational<T> result = a;
    return result += b;
}

template <class T> inline rational<T> operator - (const rational<T> a, const rational<T> b)
{
    rational<T> result = a;
    return result -= b;
}


template <class T> inline rational<T> operator * (const rational<T> a, const rational<T> b)
{
    rational<T> result = a;
    return result *= b;
}

template <class T> inline rational<T> operator / (const rational<T> a, const rational<T> b)
{
    rational<T> result = a;
    return result /= b;
}

template <class T> inline rational<T> operator % (const rational<T> a, const rational<T> b)
{
    return rational<T>(); //i.e zero
}

template <class T> inline bool operator != (const rational<T> a, const rational<T> b)
{
    return (!(a==b));
}

template<class T> inline bool operator >= (const rational<T> a, const rational<T> b)
{
	return !(a<b);
}

template<class T> inline bool operator < (const rational<T> a, const rational<T> b)
{
	if (a == b || a > b )
		return false;
	else
		return true;
}

template<class T> inline bool operator <= (const rational<T> a, const rational<T> b)
{
	return !(a>b);
}

template<class T> inline rational<T> abs (const rational<T>& a)
{
	rational<T> loc = a;
	loc.n = abs(loc.n);
	return loc;
}


#include<rational.int.h>


template <class T> inline rational<T> rational<T>::operator += (const rational<T> a)
{
	if (d == a.d)
	{
		n += a.n;
	}
	/* don't optimize here, it becomes ambiguous for some T 
	else if (d%a.d == T(0))
	{
		n += d/a.d * a.n;
	}
	else if (a.d%d == T(0))
	{
		n *= a.d/d;
		n += a.n;
		d = a.d;	
	}
	*/
	else
	{
	    n *= a.d;
	    n += d*a.n;
	    d *= a.d;
	}

    return scaled_rational(*this);
}

template <class T> inline rational<T> rational<T>::operator -= (const rational<T> a)
{
	if (d == a.d)
	{
		n -= a.n;
	}
	/* don't optimize here, it becomes ambiguous for some T 
	else if (d/a.d != T(0) && d%a.d == T(0))
	{
		n -= d/a.d * a.n;
	}
	else if (a.d/d != T(0) && a.d%d == T(0))
	{
		n *= a.d/d;
		n -= a.n;
		d = a.d;	
	}
	*/
	else
	{
	    n *= a.d;
	    n -= d*a.n;
	    d *= a.d;
	}

    return scaled_rational(*this);
}
                           
template <class T> inline rational<T> rational<T>::operator *= (const rational<T> a)
{
	n *= a.n;
	d *= a.d;
    return scaled_rational(*this);
}

template <class T> inline rational<T> rational<T>::operator /= (const rational<T> a)
{
	n *= a.d;
	d *= a.n;
    return scaled_rational(*this);
}

template <class T> inline rational<T> rational<T>::operator %= (const rational<T> a)
{
	n = T(0);
	d = T(1);
    return *this;
}

template <class T> inline bool rational<T>::operator == (const rational<T> Q) const
{
    if (n == Q.n && d == Q.d)
		return true;
    else
        return false;
}

template <class T> inline bool rational<T>::operator > (const rational<T> Q) const
{
    if (n*Q.d > Q.n*d)
		return true;
    else
        return false;
}
