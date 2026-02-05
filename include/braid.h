/**********************************************************************
This header file defines local programme types used to control objects
used by the programme 

 **********************************************************************/
#include <hash-defs.h>
#include <braid-control.h>
#include <debug-control.h>
//#include <braid-util.h>
//#include <generic-code.h>

//#include <seifert.h>


struct alpha_char 
{
    bool operator()(char c) {return isalpha(c);}
};

//typedef polynomial<scalar, string, rational<int> > Cpolynomial;
typedef polynomial<int, string, scalar > Cpolynomial;

class colouring_data
{
public:
	int num_fixed_points;
	int num_cocycle_invariants;
//	vector<polynomial<int> > cocycle_invariant;
	vector<Cpolynomial> cocycle_invariant;

	colouring_data():num_fixed_points(0), num_cocycle_invariants(0){}
};

class generic_switch_data
{
public:
	string title;
	string definition;
	int size; // records the size of finite biracks
	int num_chain_generators; // the number of chain group generators in chain map
	bool cocycles_calculated;
	bool biquandle;  //used by finite switches and the BIRACK_POLYNOMIAL (for biquandles, period = 1 and even for doubled biracks we don't need to check W)
	
	/* for finite biracks we record the  up and down actions, so that when evaluating the cocycle
	   invariant we can pass the doubled birack to fixed_point_invariant along with the generic_switch_data
	   of the underlying birack and retain the up and down actions for cohomology generator evaluation
	*/
	matrix<int> Su;
	matrix<int> Sd;
	list<string> cocycle_string; 
	list<vector<scalar> > cocycle_scalar; 
	vector<int> chain_map;  // maps n-tuples to non-degenerate n-tuples for use with the cocycle invariant.
	
	generic_switch_data():title(""), definition(""),size(0),num_chain_generators(0),cocycles_calculated(false),biquandle(false){}
	generic_switch_data(string t,string s): title(t), definition(s),size(0),num_chain_generators(0),cocycles_calculated(false),biquandle(false){}
	
};

void print_switch_data(ostream& s, generic_switch_data& switch_data, string prefix="");

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
typedef polynomial<Quaternion, char> Hpolynomial;
typedef matrix<Hpolynomial,Quaternion> Hpmatrix;

typedef Rational<polynomial<scalar,char> > Qpolynomial; 
typedef matrix<Qpolynomial, scalar> Qpmatrix;
typedef polynomial<int> hpolynomial; // Homfly polynomial

/* rational polynomial specializations */

template <typename T, typename V, typename E> ostream& operator << (ostream& os, const rational<polynomial<T,V,E> >& Q)
{
	rational<polynomial<T,V,E> > loc = Q;
	
    if (loc.getd() == T(1))
		os << loc.getn();
    else
    {
		if (loc.getn().numterms()>1)
			os << '(';
			
	    os << loc.getn();
		
		if (loc.getn().numterms()>1)
			os << ')';
		
		os << '/';
		
		if (loc.getd().numterms()>1)
			os << '(';

		os << loc.getd();

		if (loc.getd().numterms()>1)
			os << ')';		
    }
	
	return os;
}

template <typename T, typename V, typename E> rational<polynomial<T,V,E> > scaled_rational(rational<polynomial<T,V,E> >& Q)
{
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,E> >::scaled_rational: scaling " << Q.getn() << "/" << Q.getd() << endl;

	if (Q.getn().is_zero())
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

		/* We only carry out these tests if the numerator and denominator both have integral
		   coefficients and there are non-unit terms in both polys 
		*/
		
		polynomial<T,V,E>& numerator = Q.getn();
		polynomial<T,V,E>& denominator = Q.getd();
		
		if (!numerator.is_unit() && !denominator.is_unit())
		{
			/* First try to cancel as many variables as possible.  We take the minimum degree of each variable 
			   in the terms of the numerator and denominator.  We have to look out for the two polys having different variables.
		    */
		
			vector<E> min_num_deg = numerator.min_degrees();
			vector<E> min_den_deg = denominator.min_degrees();
			
			/* reduce those variables in the numerator by those in the denominator if necessary, and set to 1 any 
			   variable from the numerator not in the denominator.
			*/
//			const vector<char>& d_vars = denominator.getvars();
			for (int i = 0; i < numerator.nv; i++)
			{
				char var = numerator.getvar(i);
				vector<char>::const_iterator d_pos = find(denominator.vc.begin(),denominator.vc.end(),var);
				if (d_pos != denominator.vc.end())
				{
					int index = d_pos - denominator.vc.begin();
					if (min_num_deg[i] > min_den_deg[index])
						min_num_deg[i] = min_den_deg[index];
				}
				else
					min_num_deg[i] = 0;
			}

			/* now do the same for the denominator.
			*/
//			const vector<char> n_vars = numerator.getvars();
			for (int i = 0; i < denominator.nv; i++)
			{
				char var = denominator.getvar(i);
				vector<char>::const_iterator n_pos = find(numerator.vc.begin(),numerator.vc.end(),var);
				if (n_pos != numerator.vc.end())
				{
					int index = n_pos - numerator.vc.begin();
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
				if (min_num_deg[i] != E(0))
					num_oss << numerator.getvar(i);
				if (min_num_deg[i] > E(1))
					num_oss << '^' << min_num_deg[i];
			}

			if (num_oss.str().size())
			{
				polynomial<T,V,E> p_factor = num_oss.str();
			
 if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::scaled_rational: cancelling " << p_factor << " from numerator" << endl;

				numerator /= p_factor;

			}


			ostringstream den_oss;
			for (unsigned int i=0; i< min_den_deg.size(); i++)
			{
				if (min_den_deg[i] != E(0))
					den_oss << denominator.getvar(i);
				if (min_den_deg[i] > E(1))
					den_oss << '^' << min_den_deg[i];
			}

			if (den_oss.str().size())
			{
				polynomial<T,V,E> p_factor = den_oss.str();
			
 if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::scaled_rational: cancelling " << p_factor << " from denominator" << endl;

				denominator /= p_factor;

			}
		
			/* The next test substitutes each variable X with 1 to see if X=1 is a root
			   of both the numerator and denominator, if so we cancel 1-X from both.  
			   The test is applied repeatedly to catch (1-X)^n
			*/

			if (denominator.nv != 0)
			{

if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,E> >::scaled_rational: root test" << endl;

				bool found;
				do
				{
//					const vector<char>& d_vars = denominator.getvars();  // the variables may have changed 
					found = false;
					for (int i=0; i< numerator.nv; i++)
					{
						char var = numerator.getvar(i);
						vector<char>::const_iterator d_pos = find(denominator.vc.begin(),denominator.vc.end(),var);

if (rational_control::DEBUG)
{
	debug << "rational<polynomial<T,V,E> >::scaled_rational:   numerator variable " << var << " in denominator = ";
	debug << (d_pos != denominator.vc.end() ? "true" : "false" ) << endl;
}
					
						if (d_pos != denominator.vc.end() && numerator.has_root_1(var))
						{

if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::scaled_rational:     variable is root of numerator" << endl;
	
							if (denominator.has_root_1(var))
							{
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::scaled_rational:     variable is root of denominator" << endl;
								found = true;
								polynomial<T,V,E> p_factor(string("1-")+var);
								numerator /= p_factor;
								denominator /= p_factor;
							}
						}
					}
				} while (found && denominator.nv != 0);
			}
		}
		
		/* check for the denominator being -1 */
		polynomial<T,V,E> minus_one("-1");	
		if (denominator == minus_one)
		{
			numerator *= minus_one;
			denominator *= minus_one;			
		}	
	}
	
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,E> >::scaled_rational: returning " << Q << endl;
	return Q;
}

template <typename T, typename V, typename E> rational< polynomial<T,V,E> > simplify_rational (rational< polynomial<T,V,E> >& Poly)
{
	rational<polynomial<T,V,E> > result = Poly;

	if (Poly.getn() != polynomial<T,V,E>(T(0)) && Poly.getd() != polynomial<T,V,E>(T(1)))
	{
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::simplify_rational: simplifying " << Poly << endl;

		/* take the gcd of the numerator coefficients */
		polynomial<T,V,E>& numerator = Poly.getn();
		typename list<pterm<T,E> >::iterator pterm_ptr = numerator.pt.begin();
		T g = abs(pterm_ptr->n);
		pterm_ptr++;
		while (pterm_ptr != numerator.pt.end())
		{
			g = gcd(g, abs(pterm_ptr->n));
			pterm_ptr++;
		}
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::simplify_rational:   gcd of numerator terms = " << g << endl;

		if ( g != T(1))
		{
				
			polynomial<T,V,E>& denominator = Poly.getd();
			pterm_ptr = denominator.pt.begin();
			g = gcd(g, abs(pterm_ptr->n));
			pterm_ptr++;
			while (pterm_ptr != denominator.pt.end())
			{
				g = gcd(g, abs(pterm_ptr->n));
				pterm_ptr++;
			}

if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::simplify_rational:   gcd of numerator and denominator terms = " << g << endl;
			if (g != T(1))
			{	
				polynomial<T,V,E> g_poly(g);
				result.getn() /= g_poly;
				result.getd() /= g_poly;
				
				/* divide both numerator and denominator by g 
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
				*/
if (rational_control::DEBUG)
	debug << "rational<polynomial<T,V,U> >::simplify_rational:   able to cancel gcd, result is " << result << endl;
			}
		}
	}
	return result;
}

