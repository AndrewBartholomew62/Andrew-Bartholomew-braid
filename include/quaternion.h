/************************************************************************
			 Quaternian template definitions

		  A. Bartholomew September, 2004
		  
Derived from non-template definitions by  A. Bartholomew  27th February, 2002


************************************************************************/

struct quaternion_control {
	static bool DEBUG;
};

template <class T> class quaternion;


/* First define the friend functions */

template <class T> istream& operator >> (istream& s, quaternion<T>& Q)
{
	char c;

	s >> c;


if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : initial character from istream: " << c << endl;
	
	if (s.fail()) // no non-white characters have been read
	{
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : initial read failed, returning with failbit set on istream at " << &s << endl;	
		return s;
	}
	
	if (c != '(')
	{
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : no opening '(' detected" << endl;
	
		if (isdigit(c))
		{
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : digit found, putting " << c << " back into istream ready to read real part" << endl;
			s.putback(c);
			s >> Q.r;
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : setting all other quaternionic parts to zero" << endl;
			Q.i = T(0);
			Q.j = T(0);
			Q.k = T(0);
			goto quaternion_input_done;
		}
		else
		{
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : no number to read, putting back the " << c << " and setting failbit on istream at " << &s << endl;
			Q = T(1);
			s.putback(c);
			s.setstate(ios_base::failbit);
		}
		goto quaternion_input_done;
	}
	
	s >> Q.r;
	s.get(c);
	if (c != ',')
	{		
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : no comma found after reading real part setting all other quaternionic parts to zero" << endl;
		s.putback(c);
		Q.i = T(0);
		Q.j = T(0);
		Q.k = T(0);
		goto quaternion_input_done;
	}
	else
	{
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : read comma after r-part, proceeding to read i-part" << endl;
	}

	s >> Q.i;
	s.get(c);
	if (c != ',')
	{		
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : no comma found after reading i-part setting all other quaternionic parts to zero" << endl;
		s.putback(c);
		Q.j = T(0);
		Q.k = T(0);
		goto quaternion_input_done;
	}
	else
	{
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : read comma after i-part, proceeding to read j-part" << endl;
	}

	s >> Q.j;
	s.get(c);
	if (c != ',')
	{		
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : no comma found after reading j-part setting k-part parts to zero" << endl;
		s.putback(c);
		Q.k = T(0);
		goto quaternion_input_done;
	}
	else
	{
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : read comma after j-part, proceeding to read k-part" << endl;
	}

	s >> Q.k;
	s.get(c);
	if (c != ')')
	{
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : failed to read ')' after k-part, read " << c << 
	         "instead, returning it to istream and continuing" << endl;
		s.putback(c);
	}
	else
	{
if (quaternion_control::DEBUG)
	debug << "quaternion::operator >> : read ')' after k-part" << endl;
	}

	
	quaternion_input_done:
	return s;
}

template <class T> ostream& operator << (ostream& s, const quaternion<T> Q)
{
    if (Q.i == T(0) && Q.j == T(0) && Q.k == T(0))
        s << Q.r;
    else
	        s << '(' << Q.r << ',' << Q.i << ',' << Q.j << ',' << Q.k << ')';
    return s;
}


template <class T> int num_len(const quaternion<T> Q)
{
    int length;

    if ( Q.i == T(0) && Q.j == T(0) && Q.k == T(0))
        return num_len(Q.r);
    else
    {
        length = 5;

		length += num_len(Q.r);
		length += num_len(Q.i);
		length += num_len(Q.j);
		length += num_len(Q.k);
		return length;
    }
}

template <class T> quaternion<T> conj(const quaternion<T> Q)
{
    quaternion<T> result = Q;
    result.i = Q.i * T(-1);
    result.j = Q.j * T(-1);
    result.k = Q.k * T(-1);
    return result;
}

template <class T> T norm_sqrd (const quaternion<T> Q)
{
    return  Q.r*Q.r + Q.i*Q.i + Q.j*Q.j + Q.k*Q.k;
}

template <class T> inline quaternion<T> operator + (const quaternion<T> a, const quaternion<T> b)
{
	quaternion<T> result = a;
	return result+= b;
}

template <class T> inline quaternion<T> operator - (const quaternion<T> a, const quaternion<T> b)
{
	quaternion<T> result = a;
	return result-= b;
}

template <class T> inline quaternion<T> operator * (const quaternion<T> a, const quaternion<T> b)
{
	quaternion<T> result = a;
	return result*= b;
}

template <class T> inline quaternion<T> operator / (const quaternion<T> a, const quaternion<T> b)
{
	quaternion<T> result = a;
	return result/= b;
}


template <class T> inline bool operator != (const quaternion<T> a, const quaternion<T> b)
{
        return !(a==b);
}


#include<quaternion.int.h>

template <class T> inline quaternion<T> quaternion<T>::operator += (const quaternion a)
{
    r += a.r;
    i += a.i;
    j += a.j;
    k += a.k;
    return *this;
}

template <class T> inline quaternion<T> quaternion<T>::operator -= (const quaternion a)
{
    r -= a.r;
    i -= a.i;
    j -= a.j;
    k -= a.k;
    return *this;
}

template <class T> quaternion<T> quaternion<T>::operator *= (const quaternion a)
{
	quaternion<T> result;
    result.r = r*a.r - i*a.i - j*a.j - k*a.k;
    result.i = r*a.i + i*a.r + j*a.k - k*a.j;
    result.j = r*a.j - i*a.k + j*a.r + k*a.i;
    result.k = r*a.k + i*a.j - j*a.i + k*a.r;
	
	return *this = result;
}

/* this division operator calculates t/a = (\conj(a)*t)/norm_squared(a)
   it is therefore not going to return the correct answer for integers
*/
template <class T> inline quaternion<T> quaternion<T>::operator /= (const quaternion<T> a)
{

	quaternion<T> temp = (T(1)/norm_sqrd(a));
	return *this *= conj(a) * temp;
}


template <class T> inline bool quaternion<T>::operator == (const quaternion Q) const
{
    if (r == Q.r && i == Q.i && j == Q.j && k == Q.k)
        return true;
    else
        return false;
}



