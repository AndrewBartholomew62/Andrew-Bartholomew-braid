/**********************************************************************

This header file defines local programme types used to control objects
used by the programme 
**********************************************************************/

#include <hash-defs.h>
#include <braid-control.h>
#include <braid-util.h>
#include <debug-control.h>
#include <generic-code.h>


template <class T> class Rational : public rational<T>
{
public:
	Rational<T> (T t1=T(0), T t2=T(1)): rational<T> (t1, t2) {}
	Rational<T> (string str, T t2=T(1)): rational<T> (T(str), t2) {}
	Rational<T>& operator = (rational<T> r) {this->setn(r.getn()); this->setd(r.getd()); return *this;}
	Rational<T> (const rational<T> r): rational<T>(r) {}
};

/* a Quaternion is a concrete class for quaternion<scalar> */

typedef matrix<Quaternion,Quaternion> Hmatrix;
typedef polynomial<Quaternion, char, int> Hpolynomial;
typedef matrix<Hpolynomial,Quaternion> Hpmatrix;

typedef Rational<polynomial<scalar,char,bigint> > Qpolynomial; 
typedef matrix<Qpolynomial, scalar> Qpmatrix;
typedef polynomial<int> hpolynomial; // Homfly polynomial

/* rational polynomial specializations */

template <typename T, typename V, typename U> ostream& operator << (ostream& os, const rational<polynomial<T,V,U> >& Q)
{
	rational<polynomial<T,V,U> > loc = Q;
	
    if (loc.getd() == T(1))
		os << loc.getn();
    else
    {
		if (loc.getn().get_pterms()->next)
			os << '(';
			
	    os << loc.getn();
		
		if (loc.getn().get_pterms()->next)
			os << ')';
		
		os << '/';
		
		if (loc.getd().get_pterms()->next)
			os << '(';

		os << loc.getd();

		if (loc.getd().get_pterms()->next)
			os << ')';		
    }
	
	return os;
}

template <typename T, typename V, typename U> rational<polynomial<T,V,U> > scaled_rational(rational<polynomial<T,V,U> >& Q)
{
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::scaled_rational: scaling " << Q.getn() << "/" << Q.getd() << endl;

//	if (Q.getn() == T(0))
	if (Q.getn().iszero())
		Q.setd(T(1));
	else if (Q.getn() == Q.getd())
	{
		Q.setn(T(1));
		Q.setd(T(1));
	}
	else
	{
		/* call simplify poly to cancel coefficients */
		Q = simplify_rational(Q);
		
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::scaled_rational: simplified to " << Q << endl;

		/* We only carry out these tests if there are non-unit terms in both polys */
		
		polynomial<T,V,U>& numerator = Q.getn();
		polynomial<T,V,U>& denominator = Q.getd();
		pterm<T,U>* n_ptr = numerator.get_pterms();
		pterm<T,U>* d_ptr = denominator.get_pterms();
		
		if ( (n_ptr->pl != U(0) || n_ptr->next != NULL) && (d_ptr->pl != U(0) || d_ptr->next != NULL))
		{

			/* First try to cancel as many variables as possible.  We take the minimum
			   degree of each variable in the terms of the numerator and denominator.
			   We have to look out for the two polys having different variables.
		    */
		
			vector<int> min_num_deg = numerator.min_degrees();
			vector<int> min_den_deg = denominator.min_degrees();
			
			/* reduce those variables in the numerator by those in the denominator
			   if necessary, and set to 1 any variable from the numerator not in 
			   the denominator.
			*/
			char* d_vars = denominator.getvars();
			for (int i = 0; i < numerator.numvars(); i++)
			{
				char var = numerator.getvar(i);
				if (char* location=strchr(d_vars,var))
				{
					int index = location - d_vars;
					if (min_num_deg[i] > min_den_deg[index])
						min_num_deg[i] = min_den_deg[index];
				}
				else
					min_num_deg[i] = 0;
			}

			/* now do the same for the denominator.
			*/
			char* n_vars = numerator.getvars();
			for (int i = 0; i < denominator.numvars(); i++)
			{
				char var = denominator.getvar(i);
				if (char* location=strchr(n_vars,var))
				{
					int index = location - n_vars;
					if (min_den_deg[i] > min_num_deg[index])
						min_den_deg[i] = min_num_deg[index];
				}
				else
					min_den_deg[i] = 0;
			}
if (rational_control::DEBUG)
{
	debug << "rational<polynomial<T,V,U> >::scaled_rational: min_num_deg = ";
	for (unsigned int i=0; i< min_num_deg.size(); i++)
		debug << min_num_deg[i] << ' ';
	debug << endl;
	debug << "rational<polynomial<T,V,U> >::scaled_rational: min_den_deg = ";
	for (unsigned int i=0; i< min_den_deg.size(); i++)
		debug << min_den_deg[i] << ' ';
	debug << endl;
}

			/* create a poly from min_num_deg and min_den_deg and divide the
		   	numerator and denominator respectively
			*/
			ostringstream num_oss;
			for (unsigned int i=0; i< min_num_deg.size(); i++)
			{
				if (min_num_deg[i] != 0)
					num_oss << numerator.getvar(i);
				if (min_num_deg[i] > 1)
					num_oss << '^' << min_num_deg[i];
			}

			if (num_oss.str().size())
			{
				polynomial<T,V,U> p_factor = num_oss.str();
			
 if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::scaled_rational: cancelling " << p_factor << " from numerator" << endl;

				numerator /= p_factor;

			}


			ostringstream den_oss;
			for (unsigned int i=0; i< min_den_deg.size(); i++)
			{
				if (min_den_deg[i] != 0)
					den_oss << denominator.getvar(i);
				if (min_den_deg[i] > 1)
					den_oss << '^' << min_den_deg[i];
			}

			if (den_oss.str().size())
			{
				polynomial<T,V,U> p_factor = den_oss.str();
			
 if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::scaled_rational: cancelling " << p_factor << " from denominator" << endl;

				denominator /= p_factor;

			}
		
			/* The next test substitutes each variable X with 1 to see if X=1 is a root
			   of both the numerator and denominator, if so we cancel 1-X from both.  
			   The test is applied repeatedly to catch (1-X)^n
			*/

			if (denominator.numvars() != 0)
			{

if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::scaled_rational: root test" << endl;

				bool found;
				do
				{
					d_vars = denominator.getvars();  // the pointer may have changed so set it again
					found = false;
					for (int i=0; i< numerator.numvars(); i++)
					{
						char var = numerator.getvar(i);
if (rational_control::DEBUG)
{
	debug << "rational<polynomial<T,V,U> >::scaled_rational:   numerator variable " << var << " in denominator = ";
	debug << (strchr(d_vars, var)? "true" : "false" ) << endl;
}
					
						if (strchr(d_vars, var) && numerator.has_root_1(var))
						{

if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::scaled_rational:     variable is root of numerator" << endl;
	
							if (denominator.has_root_1(var))
							{
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::scaled_rational:     variable is root of denominator" << endl;
								found = true;
								polynomial<T,V,U> p_factor(string("1-")+var);
								numerator /= p_factor;
								denominator /= p_factor;
							}
						}
					}
				} while (found && denominator.numvars() != 0);
			}
		}
		
		/* check for the denominator being -1 */
		polynomial<T,V,U> minus_one("-1");	
		if (denominator == minus_one)
		{
			numerator *= minus_one;
			denominator *= minus_one;			
		}	
	}
	
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::scaled_rational: returning " << Q << endl;
	return Q;
}

template <typename T, typename V, typename U> rational< polynomial<T,V,U> > simplify_rational (rational< polynomial<T,V,U> >& Poly)
{
	rational<polynomial<T,V,U> > result = Poly;

	if (Poly.getn() != polynomial<T,V,U>(T(0)) && Poly.getd() != polynomial<T,V,U>(T(1)))
	{
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::simplify_rational: simplifying " << Poly << endl;

		/* take the gcd of the numerator coefficients */
		const polynomial<T,V,U>& numerator = Poly.getn();
		pterm<T,U>* pterm_ptr = numerator.get_pterms();
		T g = abs(pterm_ptr->n);
		pterm_ptr = pterm_ptr->next;
		while (pterm_ptr != NULL)
		{		
			g = gcd(g, abs(pterm_ptr->n));
			pterm_ptr = pterm_ptr->next;
		}
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::simplify_rational:   gcd of numerator terms = " << g << endl;

		if ( g != T(1))
		{
				
			const polynomial<T,V,U>& denominator = Poly.getd();
			pterm_ptr = denominator.get_pterms();
			do
			{		
				g = gcd(g, abs(pterm_ptr->n));
				pterm_ptr = pterm_ptr->next;
			} while (pterm_ptr != NULL);

if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::simplify_rational:   gcd of numerator and denominator terms = " << g << endl;
			if (g != T(1))
			{	
				/* divide both numerator and denominator by g */
				pterm_ptr = result.getn().get_pterms();
				do
				{
					pterm_ptr->n /= g;
					pterm_ptr = pterm_ptr->next;
				} while (pterm_ptr != NULL);
				
				pterm_ptr = result.getd().get_pterms();
				do
				{
					pterm_ptr->n /= g;
					pterm_ptr = pterm_ptr->next;
				} while (pterm_ptr != NULL);
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::simplify_rational:   able to cancel gcd, result is " << result << endl;
			}
		}
	}
	return result;
}

struct alpha_char 
{
    bool operator()(char c) {return isalpha(c);}
};
