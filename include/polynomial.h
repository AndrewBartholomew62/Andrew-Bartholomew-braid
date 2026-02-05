/*********************************************************************************************************
                     Polynomial class October 2025
                     
The new polynomial class, started in October 2025 is a complete update of the May 2020 version using 
standard library constructs.

The class supports mapped variables to enable calculation of the parity bracket polynomial and to 
facilitate a better representation of the arrow polynomial variables L_1, L_2,... and K_1, K_2,....

The class supports laurent polynomials with rational or integral exponents.

The class uses lexicographic ordering of terms rather than the old place mechanism that was only supported 
with integral exponents.

The zero polynomial is represented by any polynomial structure with an empty list of pterms.

Unit polynomials have no variables but a single pterm whose exponent map is empty.

This file includes pterm and polynomial definitions, member functions, friends and the following non-friend templates
  template <typename K, typename V> typename map<K,V>::const_iterator find_value (const map<K,V>& m, V v)
  template void Laurent_scale <> (polynomial<T,V,E>& poly);	
  template polynomial<T,V,E> substitute_signed_polynomial_variable <> (const polynomial<T,V,E>& poly,  char c1, char c2, bool change_sign);
   
**********************************************************************************************************/

/*********************************************************************************************************
                                 structure class and friend templates
**********************************************************************************************************/

template <typename T, typename V, typename E> struct poly_div_res;
template <typename T, typename V, typename E> class polynomial;

template <typename T, typename V, typename E> inline polynomial<T,V,E> operator + (const polynomial<T,V,E>& a, const polynomial<T,V,E>& b);
template <typename T, typename V, typename E> inline polynomial<T,V,E> operator - (const polynomial<T,V,E>& a, const polynomial<T,V,E>& b);
template <typename T, typename V, typename E> inline polynomial<T,V,E> operator * (const polynomial<T,V,E>& a, const polynomial<T,V,E>& b);
template <typename T, typename V, typename E> inline polynomial<T,V,E> operator / (const polynomial<T,V,E>& a, const polynomial<T,V,E>& b);
template <typename T, typename V, typename E> inline polynomial<T,V,E> operator % (const polynomial<T,V,E>& a, const polynomial<T,V,E>& b);
template <typename T, typename V, typename E> inline bool operator != (const polynomial<T,V,E>& a, const polynomial<T,V,E>& b);

template <typename T, typename V, typename E> istream& operator >> (istream& s, polynomial<T,V,E>& poly);
template <typename T, typename V, typename E> ostream& operator << (ostream& os, const polynomial<T,V,E>& poly);
template <typename T, typename V, typename E> poly_div_res<T,V,E> divide_by (polynomial<T,V,E> numerator, polynomial<T,V,E> denominator);
template <typename T, typename V, typename E> polynomial<T,V,E> gcd (const polynomial<T,V,E>& i,  const polynomial<T,V,E>& j);
template <typename T, typename V, typename E> void dump (ostream& s, const polynomial<T,V,E>& poly,string prefix="");
template <typename T, typename V, typename E> void sanitize (polynomial<T,V,E>& poly);
template <typename T, typename V, typename E> void set_to_one (polynomial<T,V,E>& poly, char ch);
template <typename T, typename V, typename E> void substitute_polynomial_variable (polynomial<T,V,E>& poly, char c1, char c2);
template <typename T, typename V, typename E> polynomial<T,V,E> Laurent_factor (polynomial<T,V,E>& poly);
template <typename T, typename V, typename E> pair<string,map<char,V> > merge_variables(const polynomial<T,V,E>& a, polynomial<T,V,E>& b);
template <typename T, typename V, typename E> polynomial<T,V,E> lub_Laurent_factor (polynomial<T,V,E> a, polynomial<T,V,E> b);



/*********************************************************************************************************
                                 structures and classes
**********************************************************************************************************/

template <typename T,typename E=int> class pterm
{
public:	
    T 			      n; /* coefficient */
    map<char,E>       e; /* exponents */
    
    pterm<T,E> ():n(T(0)){}
    pterm<T,E> (T t):n(t){}
    
    int compare_variables(const pterm<T,E>& a, const vector<char>& vc) const;
};

template <typename T, typename V=string, typename E=int> class polynomial
{
public:
    int         		nv;     /* number of variables */
    vector<char>    	vc;     /* string of variable characters */
	map<char,V> 		vm; 	/* variable map */
    list<pterm<T,E> >   pt;     /* polynomial terms */

    /* constructors */
    polynomial<T,V,E>(): nv(0) {}
	polynomial<T,V,E>(T t): nv(0) {if (t != T(0)) pt.push_back(pterm<T,E>(t));}
//	polynomial<T,V,E>(char,V);  // monomial constructor, degree=coefficient=1
	polynomial<T,V,E>(string);
	
    static void add_pterm (polynomial<T,V,E>&, const pterm<T,E>&);
//	int printlen() const;
//	int numvars() const {return nv;}
	int numterms() const {return pt.size();}
	char getvar(int i) const {return vc[i];}
//	const vector<char>& getvars() const {return vc;}
//	list<pterm<T,E> >& get_pterms() {return pt;}
//	const list<pterm<T,E> >& const_get_pterms() {return pt;}
//	pterm<T,E> last_term() {return *pt.rbegin();}
	void set_varmap(char c, V v);// set variable mapping
	V get_varmap(char c) const {return vm.at(c);} // get variable mapping
	bool laurent () const; // test for laurent polynomial
	bool is_one () const; // test for polynomial == polynomial<T,V,E>(T(1))
	bool is_minus_one () const; // test for polynomial == polynomial<T,V,E>(T(-1))
	bool is_negative () const; // test for first term of polynomial having negative coefficient
	bool is_unit () const; // test for easily identifiable unit polynomials
	bool is_zero() const {return (pt.size() == 0? true:false);}
	bool non_zero() const {return !is_zero();}
//	int getdeg() const {return pt.begin()->e.begin()->second;}  // return the maximum degree of the lexicographically greatest variable
	E getdeg() const {return pt.begin()->e.begin()->second;}  // return the maximum degree of the lexicographically greatest variable
	vector<E> min_degrees() const;
	bool has_root_1 (char c) const {polynomial<T,V,E> loc = *this; set_to_one(loc,c); return (loc.pt.size() == 0);}

    /* and the operators */
    polynomial<T,V,E> operator += (const polynomial<T,V,E>&);
    polynomial<T,V,E> operator -= (const polynomial<T,V,E>&);
    polynomial<T,V,E> operator *= (const polynomial<T,V,E>&);
    polynomial<T,V,E> operator /= (const polynomial<T,V,E>&);
    polynomial<T,V,E> operator %= (const polynomial<T,V,E>&);
    bool operator == (const polynomial<T,V,E>&) const;

	/* These are friends to allow implicit conversion of their arguments,
       helper functions do not support implicit conversion 
	*/
	friend polynomial<T,V,E> operator + <> (const polynomial<T,V,E>&, const polynomial<T,V,E>&);
	friend polynomial<T,V,E> operator - <> (const polynomial<T,V,E>&, const polynomial<T,V,E>&);
	friend polynomial<T,V,E> operator * <> (const polynomial<T,V,E>&, const polynomial<T,V,E>&);
	friend polynomial<T,V,E> operator / <> (const polynomial<T,V,E>&, const polynomial<T,V,E>&);
	friend polynomial<T,V,E> operator % <> (const polynomial<T,V,E>&, const polynomial<T,V,E>&);
	friend bool operator != <> ( const polynomial<T,V,E>&, const polynomial<T,V,E>&);
	friend istream& operator >> <> (istream& s, polynomial<T,V,E>& poly);
    friend ostream& operator << <> (ostream& os, const polynomial<T,V,E>&);
    friend poly_div_res<T,V,E> divide_by<> (polynomial<T,V,E>, polynomial<T,V,E>);
    
	friend polynomial<T,V,E> gcd <> (const polynomial<T,V,E>&,  const polynomial<T,V,E>&);
	friend void dump <> (ostream& s, const polynomial<T,V,E>& poly, string);
	friend void sanitize <> (polynomial<T,V,E>& poly);
	friend void set_to_one <> (polynomial<T,V,E>& poly, char ch);	
	friend void substitute_polynomial_variable <> (polynomial<T,V,E>& poly, char c1, char c2);	
	friend polynomial<T,V,E> Laurent_factor <> (polynomial<T,V,E>& poly);	
	friend polynomial<T,V,E> lub_Laurent_factor <> (polynomial<T,V,E> a, polynomial<T,V,E> b);
	friend pair<string,map<char,V> > merge_variables<> (const polynomial<T,V,E>& a, polynomial<T,V,E>& b);
	
};

struct polynomial_control 
{

	static unsigned int DEBUG; //bitmap
	enum parameters
	{
		general=1, 
		sanitize=2, 
		add=4, 
		multiply=8, 
		divide=16, 
		gcd=32,
		input=64,
		//'all' also supported
	};

	static bool WAIT_INFO;
	static bool MOD_P; // used in the output function to control the generation of '-' signs
	static bool OUTPUT_PROXY_VARIABLES_ONLY;  // used for inverting variables
	static bool SUBSTITUTE_MAPPED_VARIABLES;  // used for polynomial output
	static bool WRITE_PARITY_PEER_CODES; // in addition to unoriented left preferred Gauss code
	static bool TeX; // output polynomials in TeX format.
	static int wait_threshold;
};

struct polynomial_error 
{
	polynomial_error (char* message) {cout << "\nPolynomial error! " << message << endl;}
};

template <typename T, typename V, typename E> struct poly_div_res
{
   	polynomial<T,V,E> q;
   	polynomial<T,V,E> r;
};


/*********************************************************************************************************
                                 pterm friends, operators and functions
**********************************************************************************************************/

/* compare_variables compares the variables of this pterm with those of pterm a based on the lexicographical
   order determined by the order of variables in vc.  Since polynomial variables are sorted, this lexicographical 
   order is determined by the order of variables produced by the C++ sort algorithm.
   
   The function returns:
   
     1 if this pterm has variables strictly greater than those of a, 
    -1 if this pterm has variables strictly less than those of a
     0 if this pterm has variables identical to those of a
*/
template <typename T, typename E> int pterm<T,E>::compare_variables(const pterm<T,E>& a, const vector<char>& vc) const
{
if (polynomial_control::DEBUG & polynomial_control::general) 
{
	debug << "polynomial::compare_variables: comparing pterm " << *this << " with pterm " << a << endl;
	debug << "polynomial::compare_variables: vc = ";
	for (size_t i=0; i< vc.size(); i++)
		debug << vc[i] << ' ';
	debug << endl;
}

	for (size_t i=0; i< vc.size(); i++)
	{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::compare_variables:   variable " << vc[i] << endl;

		bool variable_in_this = false;
		typename map<char,E>::const_iterator this_ptr=this->e.find(vc[i]);
		if (this_ptr != this->e.end())
			variable_in_this = true;

		bool variable_in_a = false;
		typename map<char,E>::const_iterator a_ptr=a.e.find(vc[i]);
		if (a_ptr != a.e.end())
			variable_in_a = true;

if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::compare_variables:     variable_in_this " << variable_in_this << " variable_in_a  = " << variable_in_a << endl;
			
		if (variable_in_this && !variable_in_a)
		{		
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::compare_variables:     variable not in a, this is bigger " << endl;
			return 1;
		}
		else if (!variable_in_this && variable_in_a)
		{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::compare_variables:     variable not in this, this is smaller " << endl;
			return -1;
		}
		else if (variable_in_this && variable_in_a)
		{
			if (this_ptr->second > a_ptr->second)
			{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::compare_variables:     variable in both, this exponent bigger" << endl;
				return 1;
			}
			else if (this_ptr->second < a_ptr->second)
			{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::compare_variables:     variable in both, this exponent smaller" << endl;
				return -1;
			}
		}
	}
	
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::compare_variables: pterm variables identical" << endl;
	
	return 0;
}

template <typename T, typename E> pterm<T,E> operator * (const pterm<T,E>& a, const pterm<T,E>& b)
{
	pterm<T,E> result;
	result.n = a.n * b.n;
	
	if (result.n != T(0))
	{
		result.e = a.e;
		
		typename map<char,E>::const_iterator mptr = b.e.begin();
		typename map<char,E>::iterator pos;
		
		while (mptr != b.e.end())
		{
			pos = result.e.find(mptr->first);
			if (pos != result.e.end())
			{
				pos->second += mptr->second;
				if (pos->second == E(0))
					result.e.erase(pos);
			}
			else
			{
				result.e[mptr->first] = mptr->second;
			}
			
			mptr++;
		}
	}
	
	return result;
}

template <typename T, typename E> ostream& operator << (ostream& os, const pterm<T,E>& pt)
{
	if (pt.n != T(1))
		os << pt.n;
		
	typename map<char,E>::const_iterator mptr = pt.e.begin();
	while (mptr != pt.e.end())
	{
		os << mptr->first;
		
		if (mptr->second != E(1))
			os << '^' << mptr->second;
		mptr++;
	}
	return os;
}

template <typename T, typename E> void dump (ostream& s, const pterm<T,E>& p, string prefix="")
{
	s << prefix << "n = " << p.n << endl;
	s << prefix << "e = ";
	typename map<char,E>::const_iterator mptr = p.e.begin();
	while (mptr != p.e.end())
	{
		s << '(' << mptr->first << ',' << mptr->second << ") ";
		mptr++;
	}
	s << endl;
}


/*********************************************************************************************************
                                 polynomial constructors
**********************************************************************************************************/


template <typename T, typename V, typename E> polynomial<T,V,E>::polynomial(string str)
{
	istringstream ss(str);
	polynomial<T,V,E> loc;
	ss >> loc;

	nv = loc.nv;
	vc = loc.vc;
	vm = loc.vm;
	pt = loc.pt;

}


/*********************************************************************************************************
                                 utility templates
**********************************************************************************************************/

template <typename K, typename V> typename map<K,V>::const_iterator find_value (const map<K,V>& m, V v)
{
	typename map<K,V>::const_iterator mptr=m.begin();
	while (mptr != m.end())
	{
		if (mptr->second == v)
			break;
		mptr++;
	}
	return mptr;
}


/*********************************************************************************************************
                                 polynomial member functions
**********************************************************************************************************/

/* add_pterm adds new_pterm into the list of pterms of poly so that the list remains in decending lexicographical order of variables */
template <typename T, typename V, typename E> void polynomial<T,V,E>::add_pterm (polynomial<T,V,E>& poly, const pterm<T,E>& new_pterm)
{
if (polynomial_control::DEBUG & polynomial_control::general) 
{
	debug << "polynomial::add_pterm: adding pterm " << new_pterm << " to " << poly << endl;
	dump(debug,poly,"polynomial::add_pterm: ");
}

	if (poly.pt.size() == 0)
	{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::add_pterm: no pterms in poly, pushback new_pterm" << endl;
		poly.pt.push_back(new_pterm);
	}
	else
	{
		bool added = false;
		typename list<pterm<T,E> >::iterator pptr = poly.pt.begin();
		while (pptr != poly.pt.end())
		{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::add_pterm: comparing " << new_pterm << " with " << *pptr << " based on " << poly.vc.size() << " variables" << endl;

			/* check whether the new_pterm variables are greater than those of *pptr */
			int comparison = new_pterm.compare_variables(*pptr,poly.vc);
			
			if (comparison > 0)
			{
if (polynomial_control::DEBUG & polynomial_control::general) 
{
	debug << "polynomial::add_pterm: pterm comparison > 0" << endl;
	dump(debug,poly,"polynomial::add_pterm: ");
}
	
				poly.pt.insert(pptr,new_pterm); // add the new pterm in front of pptr
				added = true;
				break;
			}
			else if (comparison == 0)
			{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::add_pterm: pterm comparison = 0" << endl;
				pptr->n += new_pterm.n;  //add coefficients
				
				if (pptr->n == T(0))
				{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::add_pterm:   resultant coefficient 0, removing pterm" << endl;
					poly.pt.erase(pptr);
				}
					
				added = true;
				break;
			}
			pptr++;
		}
		if (!added)
			poly.pt.push_back(new_pterm);
	}	
}

template <typename T, typename V, typename E> void polynomial<T,V,E>::set_varmap(char c, V v)
{
	if (nv && find(vc.begin(),vc.end(),c) != vc.end())
	{
		/* remove any existing mapping to v */
		if (vm.size() != 0)
		{
			typename map<char,V>::const_iterator mptr = vm.begin();
			bool found=false;
			while (mptr != vm.end())
			{
				if (mptr->second == v)
				{
					found = true;
					break;
				}
				mptr++;
			}
			if (found)
				vm.erase(mptr);
		}
		
		vm[c]=v; // set variable mapping
	}
}


template <typename T, typename V, typename E> bool polynomial<T,V,E>::laurent() const
{
	/* work along the pterms and determine whether the term involves negative exponents	*/
	typename list<pterm<T,E> >::const_iterator pterm_ptr = pt.begin();
	
	while (pterm_ptr != pt.end())
	{
		typename map<char,E>::const_iterator mptr = pterm_ptr->e.begin();
		while(mptr != pterm_ptr->e.end())
		{
			if (mptr->second < E(0))
			{
				return true;
			}
			mptr++;
		}
		pterm_ptr++;
	}
	return false;
}

template <typename T, typename V, typename E> bool polynomial<T,V,E>::is_one() const
{
	if (nv == 0 && pt.size() == 1 && pt.begin()->n == T(1))
		return true;
	else
		return false;
}

template <typename T, typename V, typename E> bool polynomial<T,V,E>::is_minus_one() const
{
	if (nv == 0 && pt.size() == 1 && pt.begin()->n == T(-1))
		return true;
	else
		return false;
}

template <typename T, typename V, typename E> bool polynomial<T,V,E>::is_negative() const
{
	if (pt.begin()->n < T(0))
		return true;
	else
		return false;
}

/* This unit test returns true if the polynomial is +/- 1 or is a monomial with coefficient +/- 1, 
   ie, has the form  +/- x_0^{+/- i_0}...x_n^{+/- i_n}
*/
template <typename T, typename V, typename E> bool polynomial<T,V,E>::is_unit() const
{	
	if (pt.size() == 1 && (pt.begin()->n == T(1) || pt.begin()->n == T(-1)))
		return true;
	else
		return false;
}

/* min_degrees evaluates the smallest variable degrees amongst all terms of a polynomial 
*/
template <typename T, typename V, typename E> vector<E> polynomial<T,V,E>::min_degrees() const
{
	if (nv == 0)
		return vector<E>();
	
	vector<E> min_exponent(nv);
	vector<bool> min_set(nv,false);
	
if (polynomial_control::DEBUG & polynomial_control::general)
{
	debug << "polynomial::min_degrees: *this" << endl;
	dump (debug, *this, "polynomial::min_degrees:");
}

	typename list<pterm<T,E> >::const_iterator pterm_ptr = pt.begin();
	while(pterm_ptr != pt.end())
	{		
		typename map<char,E>::const_iterator mptr = pterm_ptr->e.begin();
		while (mptr != pterm_ptr->e.end())
	    {	
			int index = find(vc.begin(),vc.end(),mptr->first)-vc.begin();
			
			if (!min_set[index] || mptr->second < min_exponent[index])
				min_exponent[index] = mptr->second;
        	
        	mptr++;
    	}
	
		pterm_ptr++;
	}
	
if (polynomial_control::DEBUG & polynomial_control::general)
{
	debug << "polynomial::min_degrees: min_degrees ";
	for (int i=0; i< nv; i++)
		debug << min_exponent[i] << ' ';
	debug << " calculated for " << *this << endl;
}

	return min_exponent;
}


/*********************************************************************************************************
                                 polynomial member operators
**********************************************************************************************************/

/* 
   The operator deals with the case that the variables of a and b may be different and the fact that the 
   variable degrees may be different in a and b, even if the variables overlap.  
   
   In the case of mapped variables, it is the mapped object that determines the variable uniqueness.  Thus,
   the character keys of the two polynomials need to be aligned by the operator in the event that the same 
   mapped object v is used with different character keys by the two polynomials.
*/
template <typename T, typename V, typename E> polynomial<T,V,E> polynomial<T,V,E>::operator += (const polynomial<T,V,E>& poly)
{
    if ( poly.is_zero() )
        return *this ;

    if (this->is_zero())
	{
		*this = poly;
		return *this;
	}
	
	polynomial<T,V,E>& a = *this;	
	polynomial<T,V,E> b = poly;	
	polynomial<T,V,E> result;

if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::operator += : merging variables" << endl;

	pair <string,map<char,V> > merged_data = merge_variables(a,b);
	string variables = merged_data.first;
	result.nv = variables.length();
	
	result.vc = vector<char>(result.nv);
	for (int i=0; i< result.nv;i++)
		result.vc[i] = variables[i];

	result.vm = merged_data.second;
	
    /* add in the pterms of a and b to the result */
	typename list<pterm<T,E> >::const_iterator pptr = a.pt.begin();
    while (pptr != a.pt.end())
    {
		polynomial<T,V,E>::add_pterm(result,*pptr);			
		pptr++;
    }

	pptr = b.pt.begin();
    while (pptr != b.pt.end())
    {
		polynomial<T,V,E>::add_pterm(result,*pptr);			
		pptr++;
    }

    /* finally sanitize the result before returning */
    sanitize (result);
	
	*this = result;
	return *this;
}

template <typename T, typename V, typename E> polynomial<T,V,E> polynomial<T,V,E>::operator -= (const polynomial<T,V,E>& poly)
{
	polynomial<T,V,E> loc; 
	loc = poly;
	loc *= 	polynomial<T,V,E>(T(-1));
	
	return *this = *this + loc;
}

template <typename T, typename V, typename E> polynomial<T,V,E> polynomial<T,V,E>::operator *= (const polynomial<T,V,E>& poly)
{
	polynomial<T,V,E> result;
	polynomial<T,V,E>& a = *this;
	polynomial<T,V,E> b = poly;

if (polynomial_control::DEBUG & polynomial_control::multiply) 
	debug << "polynomial::operator *= : a = " << a << ", b = " << b << endl;

if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::operator *= : merging variables" << endl;

	pair <string,map<char,V> > merged_data = merge_variables(a,b);
	string variables = merged_data.first;
	result.nv = variables.length();
	
	result.vc = vector<char>(result.nv);
	for (int i=0; i< result.nv;i++)
		result.vc[i] = variables[i];

	result.vm = merged_data.second;

	/* multiply through the list of b pterms by the list of a pterms */
	typename list<pterm<T,E> >::iterator a_pterm_ptr = a.pt.begin();
	while (a_pterm_ptr != a.pt.end())
	{
if (polynomial_control::DEBUG & polynomial_control::multiply) 
	debug << "polynomial::operator *= : a pterm " << *a_pterm_ptr << endl;
		
		typename list<pterm<T,E> >::iterator b_pterm_ptr = b.pt.begin();
		
		while (b_pterm_ptr != b.pt.end())
		{
if (polynomial_control::DEBUG & polynomial_control::multiply) 
	debug << "polynomial::operator *= : b pterm " << *b_pterm_ptr;
			
			pterm<T,E> pterm_product = (*a_pterm_ptr)*(*b_pterm_ptr);

if (polynomial_control::DEBUG & polynomial_control::multiply) 
	debug << " product = " << pterm_product << endl;
			
			if (pterm_product.n != T(0))
				add_pterm(result,pterm_product);
				
			b_pterm_ptr++;
		}
		
		a_pterm_ptr++;
	}

if (polynomial_control::DEBUG & polynomial_control::multiply) 
	debug << "polynomial::operator *= : result = " << result << endl;
	
    /* finally sanitize the result before returning */
    sanitize (result);

if (polynomial_control::DEBUG & polynomial_control::multiply) 
	debug << "polynomial::operator *= : complete" << endl;

	*this = result;
	return *this;
}

template <typename T, typename V, typename E> inline polynomial<T,V,E> polynomial<T,V,E>::operator /= (const polynomial& poly)
{
	poly_div_res<T,V,E> result = divide_by(*this, poly);
	
	return *this = result.q;
}

template <typename T, typename V, typename E> inline polynomial<T,V,E> polynomial<T,V,E>::operator %= (const polynomial& poly)
{
	poly_div_res<T,V,E> result = divide_by(*this, poly);
	
	return *this = result.r;
}

/* The == operator uses mapped variable values rather than the particular mapped character when determining equality */
template <typename T, typename V, typename E> bool polynomial<T,V,E>::operator == (const polynomial<T,V,E>& a) const
{
    /* first check for the zero poly */
    if (pt.size() == 0 && a.pt.size() == 0)
		return true;
	else if (pt.size() == 0 || a.pt.size() == 0) // one is zero but not both
		return false;

    /* next check for the unit poly "1" */
    if (this->is_one() && a.is_one())
		return true;
	else if (this->is_one() || a.is_one())
		return false;

	if (nv != a.nv)
		return false;
	
	if (vm.size() != a.vm.size())
		return false;

	if (pt.size() != a.pt.size())
		return false;
	
	/* Check the variables are the same. Each un-mapped variable in a must appear in this and each 
	   mapped variable in a must have a corresponding variable in this mapped to the same V.  Note
	   that at most one variable in a polynomial may be mapped to given V and we know at this point 
	   of the == operator that we have the same number of variables and the same size of variable maps
	*/
	if (vm.size() == 0)
	{
		for (int i = 0; i < nv; i++)
		{
			if (vc[i] != a.vc[i])
				return false;
		}
	}
	else
	{
		for (int i = 0; i < nv; i++)
		{
			if (a.vm.count(a.vc[i]) == 0)
			{
				if (find(vc.begin(),vc.end(),a.vc[i]) == vc.end())
					return false;
			}
			else
			{
				V av = a.get_varmap(a.vc[i]);
				if (find_value(vm,av) == vm.end())
					return false;
			}
		}
	}
		
	if (vm.size() !=0)
	{	
		typename list<pterm<T,E> >::const_iterator pterm_ptr = pt.begin();
		typename list<pterm<T,E> >::const_iterator a_pterm_ptr = a.pt.begin();

		while (pterm_ptr != pt.end())
		{
			if (pterm_ptr->n != a_pterm_ptr->n || pterm_ptr->e != a_pterm_ptr->e)
				return false;
			
			pterm_ptr++;
			a_pterm_ptr++;
		}
	}
	else
	{
		polynomial<T,V,E> a_loc = a;

if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::operator == : merging variables" << endl;
		
		merge_variables(*this,a_loc);
		typename list<pterm<T,E> >::const_iterator pterm_ptr = pt.begin();
		typename list<pterm<T,E> >::const_iterator a_pterm_ptr = a_loc.pt.begin();

		while (pterm_ptr != pt.end())
		{
			if (pterm_ptr->n != a_pterm_ptr->n || pterm_ptr->e != a_pterm_ptr->e)
				return false;
			
			pterm_ptr++;
			a_pterm_ptr++;
		}
	}
	
	return true;
}


/*********************************************************************************************************
                                 polynomial friend operators
**********************************************************************************************************/

template <typename T, typename V, typename E> inline polynomial<T,V,E> operator + (const polynomial<T,V,E>& a, const polynomial<T,V,E>& b)
{
    polynomial<T,V,E> result = a;
    return result += b;
}

template <typename T, typename V, typename E> inline polynomial<T,V,E> operator - (const polynomial<T,V,E>& a, const polynomial<T,V,E>& b)
{
    polynomial<T,V,E> result = a;
    return result -= b;
}

template <typename T, typename V, typename E> inline polynomial<T,V,E> operator * (const polynomial<T,V,E>& a, const polynomial<T,V,E>& b)
{
    polynomial<T,V,E> result = a;
    return result *= b;
}

template <typename T, typename V, typename E> inline polynomial<T,V,E> operator / (const polynomial<T,V,E>& a, const polynomial<T,V,E>& b)
{
    polynomial<T,V,E> result = a;
    return result /= b;
}

template <typename T, typename V, typename E> inline polynomial<T,V,E> operator % (const polynomial<T,V,E>& a, const polynomial<T,V,E>& b)
{
    polynomial<T,V,E> result = a;
    return result %= b;
}

template <typename T, typename V, typename E> inline bool operator != (const polynomial<T,V,E>& a, const polynomial<T,V,E>& b)
{
	return !(a==b);
}

template <typename T, typename V, typename E> istream& operator >> (istream& s, polynomial<T,V,E>& poly)
{
	string inbuf;
	s >> inbuf;

if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>: presented with string: " << inbuf << endl;

    if (inbuf == "0")
    {
        /* set poly_ptr to zero polynomial */
		poly = polynomial<T,V,E>(); //construct a T from the int and assign to poly

if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>: setting zero polynomial" << endl;

    }
    else if (inbuf == "1")
    {
        /* set poly_ptr to a polynomial with one term*/
		poly = polynomial<T,V,E>(T(1));

if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>: setting polynomial to 1" << endl;

    }
    else
    {
        /* 
		   First, we check whether the polynomial is enclosed in parentheses and 
		   remove them if it is.
        */
		
		char* c_inbuf = c_string(inbuf);
    	char* cptr = c_inbuf;
			
		if (*cptr == '(')
		{

if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>: polynomial begins with '(', checking for parentheses" << endl;

			/* parentheses are detected by an open bracket by the time
			   we reach the first alpha: ((...)x...
			*/
			int count = 1;
			do
			{
				cptr++;
				if (*cptr == '(')
					count++;
				else if (*cptr == ')')
					count--;
			} while (!isalpha(*cptr));
			
			if (count)
			{
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>: open parentheses detected" << endl;
				/* remove the parentheses from c_inbuf */
				strcpy(c_inbuf,&c_inbuf[1]);
				if (c_inbuf[strlen(c_inbuf)-1] != ')')
				{
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>: no matching ')', for opening parenthesis, ignoring and continuing to read polynomial" << endl;
				}
				else
				{
					/* check that the terminating ')' is not the end of a quaternion
					   we use the same technique as above, but working from the end 
					   of the string looking for an alpha.  We expect z^2+(1,1,1,1))
					   or similar, so can apply the same logic as before.  Note that this
					   approach does not support the degenerate case where the polynomial is actually 
					   a single coefficient (4), or (3/4) or ((1,1,1,1))
					*/
					cptr = &c_inbuf[strlen(c_inbuf)-1];
					int count = 1;
					do
					{
						cptr--;
						if (*cptr == ')')
							count++;
						else if (*cptr == '(')
							count--;
					} while (!isalpha(*cptr));
			
					if (count)
					{				
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>: close parentheses detected" << endl;
						c_inbuf[strlen(c_inbuf)-1] = '\0';
					}
					else
					{
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>: closing ')', does not match opening parenthesis, continuing to read polynomial" << endl;
					}
				}
				
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>: input buffer after removal of parentheses: " << c_inbuf << endl;
				
			}
		}
		
		/* scan the input to evaluate the number of terms in the polynomial by 
		   looking for + and - signs; at the same time, determine the number of 
		   variables and the variable characters.  
		*/		
		poly.nv = 0;
		string variables;
    	int num_terms = 1;
    	bool term_level = true;
    	cptr = c_inbuf;
		do
		{
	    	/* if cptr is an alpha it points at a variable of this
		       polynomial; if so, add it to the vc list if its not already
		       there
		    */
		    if (isalpha(*cptr) && variables.find(*cptr) == string::npos)
	    	{
				poly.nv++;
				variables += *cptr;
	    	}
	    	else if (*cptr == '(')
				term_level = false; /* entering a quaternion */
	    	else if (*cptr == ')')
				term_level = true; /* leaving a quaternion */
	    	else if (cptr != c_inbuf && term_level &&
			  	(*cptr == '+' || *cptr == '-') && *(cptr-1) != '^')
				num_terms++;
		} while (*++cptr != '\0');

if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>: number of variables = " << poly.nv << endl;

        /* sort the polynomial variables */
		if (poly.nv)
		{
			sort(variables.begin(), variables.end());
			poly.vc = vector<char>(variables.length());
			for (int i=0; i< poly.nv;i++)
				poly.vc[i] = variables[i];

if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>: variable characters: " << variables << endl;
			
		}

		/* set inbuf to the current value of c_inbuf, which may have had parentheses
		   removed, before creating an istringstream for reading the polynomial */
		inbuf = string(c_inbuf);   
		istringstream iss(inbuf);

        for (int i = 0; i< num_terms; i++)
        {
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>: reading term " << i+1 << ": " << endl;

	        pterm<T,E> new_pterm;

			/* input the coefficient, start by looking for a sign:
			   we have to deal with the cases -x -2x+2 x-2y x-y x+y
			   and the >> operator fails to read -1 in the case -x
			   so we handle the sign explicitly
			*/
			char c;
			int sign;
			
			iss >> c;

if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>:   initial character read from istream = " << c << endl;
			
			if (c == '-')
			{
				sign = -1;
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>:   explicit sign = " << sign << endl;
			}
			else if (c == '+')
			{
				sign = 1;
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>:   explicit sign = " << sign << endl;
			}
			else
			{
				sign = 1;
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>:   implicit sign = " << sign << endl;
				iss.putback(c);  // this should only be the first term with no leading sign
			}
				
			
			iss >> new_pterm.n;
			if (iss.fail())
			{
				new_pterm.n = T(1);

if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>:   implicit absolute coefficient = 1" << endl;
	
				iss.clear();

if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>:   clearing state flags on istream at " << &iss << endl;
			}
			else
			{
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>:   absolute coefficient = " << new_pterm.n << endl;
			}
			
			new_pterm.n *= T(sign);
			
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>:   coefficient value = " << new_pterm.n << endl;
			
			/* read the variables and their exponents */
			do
			{
				iss >> c;
				if (iss.fail()) // the last term may have no variable
				{
					break;
				}
				else
				{
					if (isalpha(c))
					{
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>:   read variable " << c << endl;
						char carrat = 0;
						iss >> carrat;
						E exp;
						if (carrat == '^')
						{
							iss >> exp;
							
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>:   read exponent = " << exp << endl;
						}
						else // includes case where the input to carrat fails due to the polynomial ending in a variable: e.g x^2+x
						{
							exp = 1;  
							iss.putback(carrat);
if (polynomial_control::DEBUG & polynomial_control::input)
	debug << "polynomial<T,V,E>::>>:   set exponent = 1" << endl;
						}
						
						new_pterm.e[c] = exp;
					}
					else if (c == '+' || c == '-')
					{
						iss.putback(c);
						break;
					}
				}
			} while (true);

			
            /* Now add the pterm to the polynomial */
	    	polynomial<T,V,E>::add_pterm(poly,new_pterm);
		}
    }
	
if (polynomial_control::DEBUG & polynomial_control::input)
{
	debug << "polynomial<T,V,E>::>>: final polynomial structure: " << endl;
	dump(debug,poly,"polynomial dump: ");
}
	sanitize(poly);
	
	return s;
}

template <typename T, typename V, typename E> ostream& operator << (ostream& os, const polynomial<T,V,E>& poly)
{
    bool first_term = true;
    
    /* set the width of the outstream back to 0 */
    os.width(0);

    /* check for the zero polynomial */
    if (poly.pt.size() == 0)
    {
        os << '0';
        return os;
    }

	typename list<pterm<T,E> >::const_iterator pterm_ptr = poly.pt.begin();
	
	while (pterm_ptr != poly.pt.end())
    {
        /* put out the coefficient.  If this is not the unit term
           and the absolute value of the coefficient is 1 we may
           dispense with the coefficient.  Note that we cannot really talk 
		   of absolute value here, since we may have T=quaternion<X> for 
		   some X and the absolute value of a quaternion is its modulus.
		   Therefore we have to check that the coefficient is neither 1 or -1.
		   
		   We have to pay particular attention to the situation where we are working 
		   mod p.  Here the value -1 is written p-1 and we *do* need to write the
		   coefficient, provided p != 2 when p-1=1.  Thus, we write the coefficient if
		   
		  a) The place is non-zero (not the unit term)
		  b) We're working mod p and the coefficient is -1 (but not also == 1, which is the p=2 case)
		  c) We're not working mod p and the coefficient is not 1 or -1.
		  
        */

        if (!first_term && T(pterm_ptr->n) > T(0))  // we should never have a pterm with a zero coefficient
            os << '+';
        else if (!polynomial_control::MOD_P && pterm_ptr->e.size() && pterm_ptr->n == T(-1)) // pterm_ptr->e.size() zero iff unit pterm
           os << '-';

        if (pterm_ptr->e.size() == 0 ||
		    (polynomial_control::MOD_P && pterm_ptr->n == T(-1) && pterm_ptr->n != T(1)) || 
			(pterm_ptr->n != T(1) && pterm_ptr->n != T(-1))
		   ) 
		{
            os << pterm_ptr->n;
		}
			
        /* now the variables */
        
		typename map<char,E>::const_iterator mptr = pterm_ptr->e.begin();
		while (mptr != pterm_ptr->e.end())
		{
			if (poly.vm.count(mptr->first) && polynomial_control::SUBSTITUTE_MAPPED_VARIABLES && !polynomial_control::OUTPUT_PROXY_VARIABLES_ONLY)
				os << poly.get_varmap(mptr->first);
			else
				os << mptr->first;

			if (mptr->second != E(1))
			{
				if (polynomial_control::TeX)
					os << "^{";
				else
					os << '^';
				
				os << mptr->second;
				
				if (polynomial_control::TeX)
					os << '}';
			}

			mptr++;
		}
        
        first_term = false;
		pterm_ptr++;
    } 


	if (!polynomial_control::SUBSTITUTE_MAPPED_VARIABLES && poly.vm.size() && !polynomial_control::OUTPUT_PROXY_VARIABLES_ONLY)
	{
//		os << " #";
		os << " : ";
		typename map<char,V>::const_iterator mptr = poly.vm.begin();
		while (mptr != poly.vm.end())
		{
//			os << '(' << mptr->first << ',' << mptr->second << ") ";
			os << mptr->first << '=' << mptr->second << " ";
			mptr++;
		}
	}
	
	return os;
}

template <typename T, typename V, typename E> poly_div_res<T,V,E> divide_by (polynomial<T,V,E> numerator, polynomial<T,V,E> denominator)
{
    poly_div_res<T,V,E> result;

	
if (polynomial_control::DEBUG & polynomial_control::divide)
{
    debug << "polynomial::divide_by: numerator = " << numerator << endl;
    debug << "polynomial::divide_by: denominator = " << denominator << endl;
}

    if (denominator.pt.size() == 0)
    {
		/* division by zero return error */
		cout << "\n\nError! Polynomial division by zero";
		exit(0);
    }

    if (numerator.pt.size() == 0)
    {
		/* return the zero polynomial */
		result.q = polynomial<T,V,E>(T(0));
		result.r = polynomial<T,V,E>(T(0));
		return result;
    }

    /* compare variables, the variables in the denominator must be contained in the numerator if a quotient is to exist.  
       We use merge_variables to align the variables of the denominator to the numerator, in case we have mapped variables
    */
    bool variables_ok = true;
	
    if (denominator.nv > numerator.nv)
		variables_ok = false;
    else
    {

if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "polynomial::divide_by : merging variables" << endl;
		
		merge_variables(numerator,denominator);
		for (int i=0; i< denominator.nv; i++)
		{
	    	if (find(numerator.vc.begin(), numerator.vc.end(), denominator.vc[i]) == numerator.vc.end())
	    	{
				variables_ok = false;
if (polynomial_control::DEBUG & polynomial_control::divide)
	debug << "polynomial::divide_by: variables not suitable for division" << endl;
				break;
	    	}
		}
    }

    /* If the denominator contains a variable not present in the numerator the quotient is zero and the 
       remainder is equal to the numerator
    */
    if ( !variables_ok)
    {
		result.q = polynomial<T,V,E>(T(0));
		result.r = numerator;
		return result;
    }

    /* so, here a quotient really does exist */

    /* we take the variables of the quotient to be those of the numerator, we can tidy up later. */
    result.q.nv = numerator.nv;
	result.q.vc = numerator.vc;
	result.q.vm = numerator.vm;

    /* form the quotient, we make a working copy of the numerator */
    polynomial<T,V,E> num = numerator;

    /* Next, prepare quot_term to accept the terms of the quotient.  
    polynomial<T,V,E> quot_term;
    quot_term.nv = numerator.nv;
	quot_term.vc = numerator.vc;
	quot_term.vm = numerator.vm;
	*/

    /* now start the division */
    while (num.pt.size() > 0)
    {
		/* evaluate the next quotient term, if it exists. */

		/* During the course of the enclosing loop num may be sanitized so that there is a variable in 
		   the denominator that does not appear in num any longer.
		*/
		bool division_possible = true;
		if (denominator.nv > num.nv)
	    	division_possible = false;
		else
		{
			// XX CHECK - DO WE NEED TO LOOK ONLY FOR THE VARIABLES IN THE FIRST TERM OF THE DENOMINATOR?
		    for (int i=0; i< denominator.nv; i++)
	    	{
		    	if (find(num.vc.begin(), num.vc.end(), denominator.vc[i]) == num.vc.end())
				{
		    		division_possible = false;
if (polynomial_control::DEBUG & polynomial_control::divide)
	debug << "polynomial::divide_by: numerator now contains variables not suitable for division" << endl;
		    		break;
				}
	    	}
		}

		if (!division_possible)
		    break;

		/* the term in the numerator with the highest degree is at the start of the pterm list, start by checking that
		   the coefficients divide exactly so we can create the next quotient term.
		*/
		pterm<T,E> quotient_pterm;

	    if (num.pt.begin()->n % denominator.pt.begin()->n != T(0))
    	{
			break;
    	}
    	else
    	{
			quotient_pterm.n = num.pt.begin()->n / denominator.pt.begin()->n;
		}

		/* The new term, if it exists, has exponents equal to the difference between the num and denominator exponent values.  
		   If a denominator exponent value is greater then the corresponding num value, the division is not possible.
		*/
		
		for (int i=0; i<numerator.nv; i++)
		{			
			if (num.pt.begin()->e.find(numerator.vc[i]) == num.pt.begin()->e.end())
				continue;  // variable does not appear in the first term of num
				
if (polynomial_control::DEBUG & polynomial_control::divide)
	debug << "polynomial::divide_by:   variable " << numerator.vc[i] << endl;

			/* if this variable character does not exist in the denominator, the exponent of numerator.vc[i] in the quotient pterm is 
			   its exponent in the greatest pterm of num.  If it does exist in the denominator, it is the difference between 
			   its exponent values in the greatest pterms of num and the denominator 
			*/
			quotient_pterm.e[numerator.vc[i]] = num.pt.begin()->e.find(numerator.vc[i])->second;
			
			typename map<char,E>::iterator dpos = denominator.pt.begin()->e.find(numerator.vc[i]);
			if (dpos != denominator.pt.begin()->e.end())
			{
if (polynomial_control::DEBUG & polynomial_control::divide)
	debug << "polynomial::divide_by:     variable exists in first term of denominator with exponent " << dpos->second << endl;
				quotient_pterm.e[numerator.vc[i]] -= dpos->second;														
			}
			
	    	if (quotient_pterm.e[numerator.vc[i]] < E(0))
	    	{
				division_possible = false;
if (polynomial_control::DEBUG & polynomial_control::divide)
	debug << "polynomial::divide_by: degrees not suitable to create quotient term" << endl;
				break;
		    }
			else if (quotient_pterm.e[numerator.vc[i]] == E(0))
			{
				/* variable cancels in the quotient, so remove it from the exponent map */
				typename map<char,E>::iterator qpos = quotient_pterm.e.find(numerator.vc[i]);

if (polynomial_control::DEBUG & polynomial_control::divide)
	debug << "polynomial::divide_by:     variable cancel in the quotient term" << endl;

				quotient_pterm.e.erase(qpos);
			}		    
		}

		if (!division_possible)
		    break; /* jump out asap */

if (polynomial_control::DEBUG & polynomial_control::divide)
    debug << "polynomial::divide_by: quotient term = " << quotient_pterm << endl;

		/* evaluate the product of quotient_pterm and denominator then subtract the product of the quotient term
		   and the denominator from the numerator.
		*/

		polynomial<T,V,E> product;
		polynomial<T,V,E>::add_pterm(product, quotient_pterm);		
		product *= denominator;

if (polynomial_control::DEBUG & polynomial_control::divide)
    debug << "polynomial::divide_by: product = " << product << endl;

		num -= product;

if(polynomial_control::DEBUG & polynomial_control::divide)
    debug << "polynomial::divide_by: new numerator = " << num << endl;

		/* finally add the quotient_pterm to the quotient */
		polynomial<T,V,E>::add_pterm(result.q, quotient_pterm);
    }

    /* if we have broken out of the above loop there is a non zero remainder, equal to num. */
    if (num.pt.size() > 0)
    {
		result.r = num;
	    sanitize(result.r);
    }
    else
    {
		/* remainder is zero */
		result.r = polynomial<T,V,E>(T(0));
    }
	
	sanitize(result.q);
	sanitize(result.r);


if (polynomial_control::DEBUG & polynomial_control::divide)
{
	debug << "polynomial::divide_by: result: " << endl;
	debug << "polynomial::divide_by:   numerator = " << numerator << endl;
	debug << "polynomial::divide_by:   denominator = " << denominator << endl;
	debug << "polynomial::divide_by:   quotient = " << result.q << endl;
	debug << "polynomial::divide_by:   remainder = " << result.r << endl;
}
    return result;
}



/*********************************************************************************************************
                                 polynomial friend functions
**********************************************************************************************************/

/* this gcd template is a specialization to handle
	- multivariate polynomials
    - polynomials with differente variables
    - laurent polynomials
   returning 1 in each case.
*/
template <typename T, typename V, typename E> polynomial<T,V,E> gcd (const polynomial<T,V,E>& i,  const polynomial<T,V,E>& j)
{
if (polynomial_control::DEBUG & polynomial_control::gcd)
	debug << "polynomial<T,V,U>::gcd: i = " << i << " j = " << j << endl;
	if ( i.nv > 1 || j.nv > 1 || (i.nv && j.nv && i.vc != j.vc)  || i.laurent() || j.laurent() )
	{
if (polynomial_control::DEBUG & polynomial_control::gcd)
	debug << "polynomial<T,V,U>::gcd: returning 1" << endl;
		return polynomial<T,V,E>(T(1));
	}
	else
	{

	    if (i.non_zero())
			if (j.non_zero())
			{
				if ((i/j).is_zero() && (j/i).is_zero()) // the i%j = i and j%i = j loop
				{
if (polynomial_control::DEBUG & polynomial_control::gcd)
	debug << "polynomial<T,V,U>::gcd: remainder loop detected, returning 1" << endl;
					return polynomial<T,V,E>(T(1));
				}
				
if (polynomial_control::DEBUG & polynomial_control::gcd)
	debug << "polynomial<T,V,U>::gcd: recursing with " << j << " and " << i%j << endl;
		    	return ( gcd<T,V,E>(j, i%j));
			}
			else
			{
if (polynomial_control::DEBUG & polynomial_control::gcd)
	debug << "polynomial<T,V,U>::gcd: returning " << i << endl;
			    return (i);
			}
	    else
		{
if (polynomial_control::DEBUG & polynomial_control::gcd)
	debug << "polynomial<T,V,U>::gcd: returning " << j << endl;
			return (j);
		}
	}
}

template <typename T, typename V, typename E> void dump (ostream& s, const polynomial<T,V,E>& poly,string prefix)
{
//    s << prefix << "polynomial located at    " << &poly << endl;
    s << prefix << "number of variables      = " << poly.nv << endl;
    s << prefix << "variable characters      = ";
    for (int i=0; i< poly.nv; i++)
		s << poly.vc[i];
	s << endl;
 
    s << prefix << "variable mapping         = ";
	if (poly.vm.size())
	{
		s << "size = " << poly.vm.size() << ": ";
		typename map<char,V>::const_iterator mptr = poly.vm.begin();
		while (mptr != poly.vm.end())
		{
			s << '(' << mptr->first << ',' << mptr->second << ") ";
			mptr++;
		}
	}
	else
	{
		s << "Empty";
	}
	s << endl;

	s << prefix << "polynomial pterms        = size " << poly.pt.size() << endl;

	if (poly.pt.size())
	{
		int index = 0;
		typename list<pterm<T,E> >::const_iterator pptr = poly.pt.begin();
		while (pptr != poly.pt.end())
		{	
			s << prefix << index << '.' << endl;
			dump (s, *pptr , prefix+"  ");
			index++;
			pptr++;
		}
	}
	else
	{
		s << prefix << "zero polynomial" << endl;
	}
}

/* sanitize checks that all of the variable characters of a polynomial are included 
   amongst the pterms.  If not, they are removed from the polynomial
*/
template <typename T, typename V, typename E> void sanitize (polynomial<T,V,E>& poly)
{

if (polynomial_control::DEBUG & polynomial_control::sanitize)
{
    debug << "polynomial::sanitize: sanitizing " << poly << endl;
    dump(debug,poly,"polynomial::sanitize:   ");
}

	vector<char>::iterator cptr = poly.vc.begin();
	while (cptr != poly.vc.end())
	{
		bool included = false;
		
		typename list<pterm<T,E> >::iterator pptr = poly.pt.begin();
		while (pptr != poly.pt.end())
		{	
			if (pptr->e.find(*cptr) != pptr->e.end())
			{
				included = true;
				break;
			}
			pptr++;
		}
		
		if (!included)
		{
if (polynomial_control::DEBUG & polynomial_control::sanitize)
    debug << "polynomial::sanitize: variable " << *cptr << " not included in pterm list " << endl;
			
			/* check if *cptr is a mapped variable */
			typename map<char,V>::iterator mptr = poly.vm.find(*cptr);
			if (mptr != poly.vm.end())
				poly.vm.erase(mptr);
				
			poly.vc.erase(cptr);
			poly.nv--;
		}
		else
		{
			cptr++;
		}
	}

if (polynomial_control::DEBUG & polynomial_control::sanitize)
{
    debug << "polynomial::sanitize: produced " << poly << endl;
    dump(debug,poly,"polynomial::sanitize:   ");
}
}

/* The function set_to_one sets the variable ch in poly to one */
template <typename T, typename V, typename E> void set_to_one (polynomial<T,V,E>& poly, char ch)
{

if (polynomial_control::DEBUG & polynomial_control::general)
    debug << "polynomial::set_to_one: setting variable " << ch << " to 1 in poly " << poly << endl;

    /* if this is the zero polynomial or if there are no variables, or if 
	   ch is not a variable of poly there's nothing to do
    */
    vector<char>::iterator ch_ptr = find(poly.vc.begin(),poly.vc.end(),ch);
    
    if (poly.pt.size() == 0 || poly.nv == 0 || ch_ptr == poly.vc.end())
	{
if (polynomial_control::DEBUG & polynomial_control::general)
    debug << "polynomial::set_to_one: nothing to do" << endl;
		return;
	}

	polynomial<T,V,E> new_poly;
	new_poly.nv = poly.nv -1;
	new_poly.vc = poly.vc;
	new_poly.vc.erase(find(new_poly.vc.begin(),new_poly.vc.end(),ch));
	new_poly.vm = poly.vm;

	/* if ch is mapped in poly, remove it from the variable map of new_poly */
	typename map<char,V>::iterator vm_ptr = poly.vm.find(ch);
	if (vm_ptr != poly.vm.end())
		new_poly.vm.erase (new_poly.vm.find(ch));
	
	int index = ch_ptr - poly.vc.begin();
	
if (polynomial_control::DEBUG & polynomial_control::general)
{
    dump(debug, poly,"polynomial::set_to_one: ");
	debug << endl;
    debug << "polynomial::set_to_one: variable index = " << index << endl;
}

	typename list<pterm<T,E> >::iterator pterm_ptr = poly.pt.begin();
	while (pterm_ptr != poly.pt.end())
	{
		/* remove ch from exponent map of this pterm, if it exists and add what is left to new_poly */
		typename map<char,E>::iterator mptr = pterm_ptr->e.find(ch);
		if (mptr != pterm_ptr->e.end())
			pterm_ptr->e.erase(mptr);
			
		polynomial<T,V,E>::add_pterm(new_poly,*pterm_ptr);
		
		pterm_ptr++;
	}

	poly = new_poly;
	
if (polynomial_control::DEBUG & polynomial_control::general)
{
    debug << "polynomial::set_to_one: resultant polynomial = " << poly;
	dump (debug, poly);
	debug << endl;
}

}

/* The function substitute_polynomial_variable substitutes character c1 for character c2 in poly. */
template <typename T, typename V, typename E> void substitute_polynomial_variable (polynomial<T,V,E>& poly, char c1, char c2)
{
	/* if c1 doesn't exist in poly, or if c2 is already used in poly, do nothing */
	if (find(poly.vc.begin(),poly.vc.end(),c1) == poly.vc.end() || find(poly.vc.begin(),poly.vc.end(),c2) != poly.vc.end())
//	if (!strchr(poly.vc,c1) || strchr(poly.vc,c2))
		return;

if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "substitute_polynomial_variable: substituting " << c1 << " for " << c2 << " in " << poly <<  endl;
		
	ostringstream oss;
	bool substitute_hold = polynomial_control::SUBSTITUTE_MAPPED_VARIABLES;
	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = false;
	oss << poly;
	
	string str_poly = oss.str();
	substitute_string_variable(str_poly,c1,c2);
	
	map<char,V> poly_map = poly.vm;
	if (poly.vm.count(c1))
	{
		V c1_V = poly_map[c1];
		poly_map.erase(c1);
		poly_map[c2] = c1_V;
	}
	
	poly = polynomial<T,V,E>(str_poly);
	poly.vm = poly_map;
	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = substitute_hold;
	
/*if (polynomial_control::DEBUG & polynomial_control::general) 
{
	debug << "substitute_polynomial_variable: polynomial changed to ";
	dump(debug,poly);
}
*/
}

/* Laurent_factor returns the smallest monomial that when multiplied by poly yields 
   a polynomial with no negative degrees 
*/
template <typename T, typename V, typename E> polynomial<T,V,E> Laurent_factor (polynomial<T,V,E>& poly)
{
if (polynomial_control::DEBUG & polynomial_control::general)
{
	debug << "polynomial::Laurent_factor: Polynomial = " << poly << endl;
	dump(debug,poly,"polynomial::Laurent_factor:   ");
}
	
	pterm<T,E> pterm_factor(T(1));
	vector<E> exponents(poly.nv);
	typename list<pterm<T,E> >::iterator pterm_ptr = poly.pt.begin();
	
	int index=0;
	
	while (pterm_ptr != poly.pt.end())
	{
if (polynomial_control::DEBUG & polynomial_control::general)
	debug << "polynomial::Laurent_factor: pterm " << index << ": " << *pterm_ptr << endl;
	
		typename map<char,E>::iterator mptr = pterm_ptr->e.begin();
		while (mptr != pterm_ptr->e.end())
		{
			int vc_index = find(poly.vc.begin(),poly.vc.end(),mptr->first) - poly.vc.begin();
if (polynomial_control::DEBUG & polynomial_control::general)
	debug << "polynomial::Laurent_factor:   mptr " << mptr->first << " vc_index = " << vc_index << endl;

			if (mptr->second < exponents[vc_index])
			{
if (polynomial_control::DEBUG & polynomial_control::general)
	debug << "polynomial::Laurent_factor:   replace exponents[ " << vc_index << "] with " << mptr->second << endl;
				exponents[vc_index] = mptr->second;
			}
				
			mptr++;
		}
		
		pterm_ptr++;
	}

if (polynomial_control::DEBUG & polynomial_control::general)
{
	debug << "polynomial::Laurent_factor:   final exponents ";
	for (int i=0; i< poly.nv; i++)
		debug << exponents[i] << ' ';
	debug << endl;
}
	
	for (int i=0; i< poly.nv; i++)
	{
		if (exponents[i] < E(0))
			pterm_factor.e[poly.vc[i]] = abs(exponents[i]);
	}
	
	polynomial<T,V,E> result;
	result.nv = poly.nv;
	result.vc = poly.vc;
	polynomial<T,V,E>::add_pterm(result,pterm_factor);
	sanitize(result); // tidy up the required variables
	
	
if (polynomial_control::DEBUG & polynomial_control::general)
{
	debug << "polynomial::Laurent_factor: Laurent scaling factor = " << result << endl;
}
	return result;
}

//  
/* merge_variables merges the variables and their mappings between two polynomials a and b, returning the merged
   variables and their mappings as a pair<string,map<char,V> >.

   For mapped variables, the mappings in a take precedence; that is, if v_a is mapped to V in a and v_b is mapped
   to V in b, in the merged variables, we use v_a to represent V, changing b accordingly.  Similarly, if v_b is 
   used as an unmapped variable in a, we choose a new v_b' for the mapping in b.

There are multiple cases depending upon the overlap of variables

   variable v_b in b does not appear in a
   --------------------------------------
	if v_b mapped in b // case 1   
		if there is a variable v_a in a mapped to the same V. // case 1.1
			if v_a appears in b // case 1.1.1
				substitute v_a in b for a new variable unused in a or b
				substitute v_b in b for v_a 
			else // case 1.1.2
				substitute v_b in b for v_a
		else // case 1.2
			add v_b and its mapping to a
    else // case 2
	  add v_b to the variable characters of a.  

	  
   variable v in b also appears in a
   ---------------------------------

	if v mapped in a // case 3
		if v mapped in b // case 3.1
			if v maps to the same V in both a and b // case 3.1.1
				nothing to do
			else // case 3.1.2
				if there is another variable v' in a that maps to the same V as v does in b // case 3.1.2.1
					if v' appears in b // case 3.1.2.1.1
						substitute v' in b for a new variable unused in a and b
						substitute v for v' in b
					else // case 3.1.2.1.2
						substitute v for v' in b
				else // case 3.1.2.2
					substitute v in b for a new variable unused in a or b 
					add the new variable and its mapping from b to a		   
	   else // case 3.2
			substitute v in b for a new variable unused in a or b
			add the new variable to a
	else // case 4
		if v mapped in b // case 4.1
				if there is another variable v' in a that maps to the same V as v does in b // case 4.1.1
					if v' appears in b // case 4.1.1.1
						substitute v' in b for a new variable unused in a and b
						substitute v for v' in b
					else // case 4.1.1.2
						substitute v for v' in b
				else // case 4.1.2
					substitute v in b for a new variable unused in a or b 
					add the new variable and its mapping from b to a		   
		else // case 4.2
			nothing to do   
*/
template <typename T, typename V, typename E> pair<string,map<char,V> > merge_variables(const polynomial<T,V,E>& a, polynomial<T,V,E>& b)
{	
	string variables;
	map<char,V> new_vm = a.vm;
	
	bool hold = polynomial_control::SUBSTITUTE_MAPPED_VARIABLES;
	
if (polynomial_control::DEBUG & polynomial_control::general) 
{
    polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = false;
	debug << "merge_variables: given polynomials a = " << a << " b = " << b <<  endl;
	debug << "merge_variables: dump of  a:" << endl;
	dump(debug,a,"merge_variables:   ");
	debug << "merge_variables: dump of  b:" << endl;
	dump(debug,b,"merge_variables:   ");
}	

	if (a.nv)
		variables = string(a.vc.begin(),a.vc.end());
	
	/* work with a copy of the variables in b, since substituting variables can
	   re-order the variables stored in the polynomial structure.
	*/
	string b_variables;
	if (b.nv !=0)
		b_variables = string(b.vc.begin(),b.vc.end());

	string combined_variables = variables+b_variables;

if (polynomial_control::DEBUG & polynomial_control::general) 
{
	debug << "merge_variables: b_variables " << b_variables <<  endl;
	debug << "merge_variables: combined_variables " << combined_variables <<  endl;
}
		
	for (int i=0; i<b.nv; i++)
    {
		char bv = b_variables[i];
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables: processing variable " << bv <<  endl;
	
		/* does the b variable appear in a? */
        if (a.vc.size() == 0 || find(a.vc.begin(),a.vc.end(),bv) == a.vc.end())
//        if (!a.vc || !strchr(a.vc,bv))
        {					
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   variable " << bv << " does not appear in a" <<  endl;
	
			
			/* is the b variable mapped */
			typename map<char,V>::const_iterator bm_ptr = b.vm.find(bv);
			if (bm_ptr != b.vm.end()) // case 1
			{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 1: variable " << bv << " is mapped in b" <<  endl;
	
				/* is there a variable v_a in a mapped to the same V. */
				typename map<char,V>::const_iterator am_ptr = find_value(a.vm,bm_ptr->second);
				if (am_ptr != a.vm.end()) // case 1.1
				{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 1.1: variable " << am_ptr->first << " in a mapped to the same V" <<  endl;
				
					/* check whether v_a appears in b */
					if (find(b.vc.begin(),b.vc.end(),am_ptr->first) != b.vc.end()) // case 1.1.1					
//					if (strchr(b.vc,am_ptr->first)) // case 1.1.1					
					{
						char new_variable = next_variable(combined_variables);
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 1.1.1: variable " << am_ptr->first << " appears in b, new_variable = " << new_variable <<  endl;
						substitute_polynomial_variable(b,am_ptr->first,new_variable);
						substitute_string_variable(b_variables,am_ptr->first,new_variable);

						substitute_polynomial_variable(b,bv,am_ptr->first);
						substitute_string_variable(b_variables,bv,am_ptr->first);
							
						combined_variables += new_variable;
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   b changed to " << b << ", b_variables = " << b_variables << endl;
					}
					else // case 1.1.2
					{
						substitute_polynomial_variable(b,bv,am_ptr->first);
						substitute_string_variable(b_variables,bv,am_ptr->first);
							
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 1.1.2: variable " << am_ptr->first << " does not appear in b, b changed to " << b << ", b_variables = " << b_variables << endl;
					}
				}
				else // case 1.2
				{
					variables += bv;				
//					a.vm[bm_ptr->first] = bm_ptr->second;
					new_vm[bm_ptr->first] = bm_ptr->second;
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 1.2: no variable in a mapped to the same V" <<  endl;
				}
			}		
			else // case 2
			{
				variables += bv;
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   cae 2: variable " << bv << " is not mapped in b" <<  endl;
			}
        }
		else
		{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   variable " << bv << " appears in a" <<  endl;
	
			/* is the variable mapped in a? */
			typename map<char,V>::const_iterator am_ptr = a.vm.find(bv);
			
			/* is the variable mapped in b? */
			typename map<char,V>::const_iterator bm_ptr = b.vm.find(bv);
			
			if (am_ptr != a.vm.end()) // case 3
			{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 3: variable " << bv << " is mapped in a to " << am_ptr->second <<  endl;
				if (bm_ptr != b.vm.end()) // case 3.1
				{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 3.1: variable " << bv << " is mapped in b to " << bm_ptr->second <<  endl;
	
					/* check whether v maps to the same V in both a and b */
					if (am_ptr->second != bm_ptr->second) // case 3.1.2
					{
						/*check whether there is another variable v' in a that maps to the same V as v does in b */
						typename map<char,V>::const_iterator am_ptr = find_value(a.vm,bm_ptr->second);
						if (am_ptr != a.vm.end()) // case 3.1.2.1
						{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 3.1.2.1: variable " << am_ptr->first << " in a is mapped to the same V as variable " << bv << " in b" <<  endl;
						
							/* check whether v' appears in b */
							if (find(b.vc.begin(),b.vc.end(),am_ptr->first) != b.vc.end()) // case 3.1.2.1.1				
//							if (strchr(b.vc,am_ptr->first)) // case 3.1.2.1.1				
							{
								char new_variable = next_variable(combined_variables);
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 3.1.2.1.1: new_variable for " << am_ptr->first << " in b is " << new_variable <<  endl;
								substitute_polynomial_variable(b,am_ptr->first,new_variable);
								substitute_string_variable(b_variables,am_ptr->first,new_variable);

								substitute_polynomial_variable(b,bv,am_ptr->first);
								substitute_string_variable(b_variables,bv,am_ptr->first);
							
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   b changed to " << b << ", b_variables = " << b_variables << endl;
							}
							else // case 3.1.2.1.2
							{
								substitute_polynomial_variable(b,bv,am_ptr->first);
								substitute_string_variable(b_variables,bv,am_ptr->first);
							
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 3.1.2.1.2: variable " << am_ptr->first << " does not appear in b, b changed to " << b << ", b_variables = " << b_variables << endl;
							}
						}												
						else // case 3.1.2.2
						{
							char new_variable = next_variable(combined_variables);					
							V b_map = bm_ptr->second;
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 3.1.2.2: new variable for " << b_map << " in b is " << new_variable << endl;
							substitute_polynomial_variable(b,bv,new_variable);
							substitute_string_variable(b_variables,bv,new_variable);
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   b updated to " << b << endl;
						
							combined_variables += new_variable;
							variables += new_variable;				
//							a.vm[new_variable] = b_map;					
							new_vm[new_variable] = b_map;					
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   new_vm[new_variable] = " << new_vm[new_variable] <<  endl;
						}
					}
					else // case 3.1.1 nothing to do					
					{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 3.1.1: variable " << bv << " maps to the same V in both a and b, nothing to do" <<  endl;
					}

				}
				else // case 3.2
				{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 3.2: variable " << bv << " is not mapped in b" <<  endl;
					char new_variable = next_variable(combined_variables);
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   new variable for " << bv << " in b is " << new_variable << endl;
					substitute_polynomial_variable(b,bv,new_variable);
					substitute_string_variable(b_variables,bv,new_variable);
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   b updated to " << b << endl;
					combined_variables += new_variable;
					variables += new_variable;				
				}
			}
			else // case 4
			{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 4: variable " << bv << " is not mapped in a" <<  endl;
	
				if (bm_ptr != b.vm.end()) // case 4.1
				{

					typename map<char,V>::const_iterator am_ptr = find_value(a.vm,bm_ptr->second);
					if (am_ptr != a.vm.end()) // case 4.1.1
					{								
						/* check whether v' appears in b */
						if (find(b.vc.begin(),b.vc.end(),am_ptr->first) != b.vc.end()) // case 4.1.1.1
//						if (strchr(b.vc,am_ptr->first)) // case 4.1.1.1
						{
							char new_variable = next_variable(combined_variables);
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   cae 4.1.1.1: new_variable for " << am_ptr->first << " in b is " << new_variable <<  endl;
							substitute_polynomial_variable(b,am_ptr->first,new_variable);
							substitute_string_variable(b_variables,am_ptr->first,new_variable);

							substitute_polynomial_variable(b,bv,am_ptr->first);
							substitute_string_variable(b_variables,bv,am_ptr->first);
							
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   b changed to " << b << ", b_variables = " << b_variables << endl;
						}
						else // case 4.1.1.2
						{
							substitute_polynomial_variable(b,bv,am_ptr->first);
							substitute_string_variable(b_variables,bv,am_ptr->first);
							
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 4.1.1.2: variable " << am_ptr->first << " does not appear in b, b changed to " << b << ", b_variables = " << b_variables << endl;
						}
					}
					else // case 4.1.2
					{
						char new_variable = next_variable(combined_variables);
						V b_map = bm_ptr->second;
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 4.1.2: new variable for " << b_map << " in b is " << new_variable << endl;
					
						substitute_polynomial_variable(b,bv,new_variable);
						substitute_string_variable(b_variables,bv,new_variable);
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   b updated to " << b << endl;
					
						combined_variables += new_variable;
						variables += new_variable;				
//						a.vm[new_variable] = b_map;					
						new_vm[new_variable] = b_map;					
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   new_vm[new_variable] = " << new_vm[new_variable] <<  endl;
					}
				}				
				else // case 4.2, nothing to do
				{
if (polynomial_control::DEBUG & polynomial_control::general) 
	debug << "merge_variables:   case 4.2: variable " << bv << " is not mapped in b, nothing to do" <<  endl;
				}
			}			
		}

//if (polynomial_control::DEBUG & polynomial_control::general) 
//	debug << "merge_variables: finished procesing variable " << bv << " i = " << i << "b.nv = " << b.nv <<  endl;
		
    }

	sort(variables.begin(), variables.end());
    
if (polynomial_control::DEBUG & polynomial_control::general) 
{
	debug << "merge_variables: produced polynomials a = " << a << " b = " << b <<  endl;
//	debug << "merge_variables: dump of  a:";
//	dump(debug,a,"merge_variables:   ");
//	debug << "merge_variables: dump of  b:";
//	dump(debug,b,"merge_variables:   ");
    polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = hold;
}
	
	return pair<string,map<char,V> >(variables,new_vm);
}

/* lub_Laurent_factor (least upper bound) takes two Laurent factors for polys a and b and returns 
   the smallest polynomial that is a Laurent factor for both a and b   
*/
template <typename T, typename V, typename E> polynomial<T,V,E> lub_Laurent_factor (polynomial<T,V,E> a, polynomial<T,V,E> b)
{

	if ( a.is_one())
		return b;
	else if ( b.is_one())
		return a;
	else
	{
		polynomial<T,V,E> a_factor = Laurent_factor(a);
		polynomial<T,V,E> b_factor = Laurent_factor(b);
		
		pair<string,map<char,V> > merge_data = merge_variables(a_factor,b_factor);
		
		polynomial<T,V,E> lub_factor;
		lub_factor.nv = merge_data.first.length();
		lub_factor.vc = vector<char>(lub_factor.nv);
		for (int i=0; i<lub_factor.nv; i++)
			lub_factor.vc[i] = merge_data.first[i];
		lub_factor.vm = merge_data.second;
		
		pterm<T,E> new_pterm(T(1));
		map<char,E>& a_map = a_factor.pt.begin()->e;
		map<char,E>& b_map = b_factor.pt.begin()->e;
		
		/*  a_factor and b_factor are sanitized polynomials, so all of their variable characters appear
		    in their (single) pterm.  Each of the merged variables therefore appears in at least one of
		    a_map or b_map and we want to create a pterm that maps to the greater of the exponent values.
		*/
		for (int i=0; i<lub_factor.nv; i++)
		{
			typename map<char,E>::const_iterator am_ptr = a_map.find(lub_factor.vc[i]);
			typename map<char,E>::const_iterator bm_ptr = b_map.find(lub_factor.vc[i]);
			
			if (am_ptr == a_map.end())
				new_pterm.e[lub_factor.vc[i]] = bm_ptr->second;
			else if (bm_ptr == b_map.end())
				new_pterm.e[lub_factor.vc[i]] = am_ptr->second;
			else
				new_pterm.e[lub_factor.vc[i]] = max(am_ptr->second,bm_ptr->second);
		}
		
		polynomial<T,V,E>::add_pterm(lub_factor,new_pterm);
		return lub_factor;				
	}
}

/*********************************************************************************************************
                                 non-friend function templates
**********************************************************************************************************/

template <typename T, typename V, typename E> void Laurent_scale (polynomial<T,V,E>& poly)
{
	poly *= Laurent_factor(poly);
	sanitize(poly);
}

/* substitute_signed_polynomial_variable substitutes c1 with c2 or, if change_sign == true, with -c2 */
template <typename T, typename V, typename E> polynomial<T,V,E> substitute_signed_polynomial_variable (const polynomial<T,V,E>& poly,  char c1, char c2, bool change_sign)
{
    /* check for the zero polynomial */
    if (poly.pt.size() == 0)
        return poly;

if (polynomial_control::DEBUG & polynomial_control::general)
{
	debug << "substitute_signed_polynomial_variable: presented with polynomial" << endl;
	dump(debug,poly,"substitute_signed_polynomial_variable: ");
}

	polynomial<T,V,E> result;
	result.nv = poly.nv;
	result.vc = poly.vc;
	
	for (int i=0; i< result.nv; i++)
	{
		if (result.vc[i] == c1)
		{
			result.vc[i] = c2;
			break;
		}
	}
	
	/* if c1 is mapped in poly, replace it with c2 in result.vm */
	result.vm = poly.vm;
	typename map<char,V>::const_iterator vm_ptr = poly.vm.find(c1);
	if (vm_ptr != poly.vm.end())
	{
		result.vm.erase(result.vm.find(c1));
		result.vm[c2] = vm_ptr->second;
	}

	typename list<pterm<T,E> >::const_iterator pterm_ptr = poly.pt.begin();
	
    while (pterm_ptr != poly.pt.end())
    {
if (polynomial_control::DEBUG & polynomial_control::general)
	debug << "substitute_signed_polynomial_variable: term "  << *pterm_ptr;

		pterm<T,E> new_pterm(*pterm_ptr);

		/* check whether this term includes c1 and, if so, substitute c1 with c2 in new_pterm and check 
		   whether the sign of the term will change.
		*/
		typename map<char,E>::const_iterator mptr = pterm_ptr->e.find(c1);
		if(mptr != pterm_ptr->e.end())
		{
			new_pterm.e.erase(new_pterm.e.find(c1));
			new_pterm.e[c2] = mptr->second;
			
			/* the sign changes iff the exponent of c1 is an odd integer */
//			if (change_sign && mptr->second.n % 2 && mptr->second.d == 1)
			if (change_sign && new_pterm.e[c2] % 2)
				new_pterm.n *= T(-1);
		}

if (polynomial_control::DEBUG & polynomial_control::general)
	debug << " becomes "  << new_pterm << endl;
		
		polynomial<T,V,E>::add_pterm(result,new_pterm);
		
        pterm_ptr++;
    }

if (polynomial_control::DEBUG & polynomial_control::general)
{
	debug << "substitute_signed_polynomial_variable: returning polynomial " << result << endl;
	dump(debug,result,"substitute_signed_polynomial_variable:   ");
}

    return result;
}

