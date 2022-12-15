/*********************************************************************************************
                     New polynomial class May 2020
                     
The new polynomial class, started in May 2020 is intended to support mapped variables to enable
calculation of the parity bracket polynomial and to facilitate a better representation of the arrow 
polynomial variables L_1, L_2,... and K_1, K_2,....

The requirements for the new polynomial class were:

   Calculations and processsing of un-mapped variable polynommials should be as fast as before 
   Complications in calculation arising from variable mappings should only be faced when necessary
    - suggests "if variable mapping then... else..."
   Require no changes to the existing input of normal polynomials
    - we should not need to map a char to a char
   Mapped polynomial variables will be established on an existing polynomial object, not as that object is created
    - with the exception that we need to create monomials (see below)
    - not all variables will have a mapping, in parity bracket we may have a character variable A and a graph variable n
      we cannot map A to a graph in the same was as with the arrow polnomial, where we could map A to "A" and, e.g. n to "K_1"
   However tempting, tidying up the pointer based variable characters, degrees and pterm list, is not required
    - unlikley we will ever need more than five or ten mapped variables.
   We will need to create a monomial of degree 1 and coefficient 1 (e.g. the analogue of P(x)=x )from an object of type V
   
*********************************************************************************************/
template<typename T, typename V=string, typename U=int> class polynomial
{
    int         nv;     /* number of variables */
    char*       vc;     /* string of variable characters */
    int*        vd;     /* array of (maximum absolute) variable degrees */
	map<char,V> vm; 	/* variable map */
    U         	l;      /* length of virtual array of pterms */
    pterm<T,U>*   p;    /* pointer to first term of polynomial */
    static polynomial<T,V,U> add_pterm (polynomial<T,V,U>&, pterm<T,U>*);

public:
	int printlen() const;
	int numvars() const {return nv;}
	int numterms() const;
	U getlength() {return l;}
	char getvar(int i) const {return vc[i];}
	char* getvars() const {return vc;}
	pterm<T,U>* get_pterms() const {return p;}
	void set_varmap(char c, V v);// set variable mapping
	V get_varmap(char c) const {return vm.at(c);} // get variable mapping
	bool laurent () const; // test for laurent polynomial
	bool isone () const; // test for polynomial == polynomial<T,V,U>(T(1))
	bool isunit () const; // test for easily identifiable unit polynomials
	bool nonzero() const {return p;}
	bool iszero() const {return !nonzero();}
	int* getdeg() const {return vd;}
	vector<int> min_degrees() const;
	bool has_root_1 (char c) const {polynomial<T,V,U> loc = *this; set_to_one(loc,c); return !loc.p;}

    /* declare the constructor and destructor functions */
    polynomial<T,V,U>();
	polynomial<T,V,U>(T);
	polynomial<T,V,U>(char,V);  // monomial constructor, degree=coefficient=1
	polynomial<T,V,U>(string);
	
    polynomial<T,V,U>(const polynomial<T,V,U>&); // copy constructor
    ~polynomial<T,V,U>();

    /* and the operators */
    polynomial<T,V,U> operator = (const polynomial<T,V,U>&); // copy assignment
    polynomial<T,V,U> operator += (const polynomial<T,V,U>&);
    polynomial<T,V,U> operator -= (const polynomial<T,V,U>&);
    polynomial<T,V,U> operator *= (const polynomial<T,V,U>&);
    polynomial<T,V,U> operator /= (const polynomial<T,V,U>&);
    polynomial<T,V,U> operator %= (const polynomial<T,V,U>&);
    bool operator == (const polynomial<T,V,U>&) const;
	friend polynomial<T,V,U> gcd <> (const polynomial<T,V,U>&,  const polynomial<T,V,U>&);
    friend ostream& operator << <> (ostream& os, const polynomial<T,V,U>&);
	friend istream& operator >> <> (istream& s, polynomial<T,V,U>& poly);
    friend poly_div_res<T,V,U> divide_by<> (polynomial<T,V,U>, polynomial<T,V,U>);
	friend void dump <> (ostream& s, const polynomial<T,V,U>& poly, string);
	friend void sanitize <> (polynomial<T,V,U>*);
	friend void set_to_one <> (polynomial<T,V,U>& poly, char ch);	
	friend void substitute_polynomial_variable <> (polynomial<T,V,U>& poly, char c1, char c2);	
	friend polynomial<T,V,U> substitute_signed_polynomial_variable <> (const polynomial<T,V,U>& poly,  char c1, char c2, bool change_sign);
	friend polynomial<T,V,U> Laurent_factor <> (polynomial<T,V,U>& poly);	
	friend void Laurent_scale <> (polynomial<T,V,U>& poly);	
	friend polynomial<T,V,U> lub_Laurent_factor <> (polynomial<T,V,U> a, polynomial<T,V,U> b);
	friend string merge_variables<> (polynomial<T,V,U>& a, polynomial<T,V,U>& b);
	
	/* These are friends to allow implicit conversion of their arguments,
       helper functions do not support implicit conversion 
	*/
	friend polynomial<T,V,U> operator + <> (const polynomial<T,V,U>&, const polynomial<T,V,U>&);
	friend polynomial<T,V,U> operator - <> (const polynomial<T,V,U>&, const polynomial<T,V,U>&);
	friend polynomial<T,V,U> operator * <> (const polynomial<T,V,U>&, const polynomial<T,V,U>&);
	friend polynomial<T,V,U> operator / <> (const polynomial<T,V,U>&, const polynomial<T,V,U>&);
	friend polynomial<T,V,U> operator % <> (const polynomial<T,V,U>&, const polynomial<T,V,U>&);
	friend bool operator != <> ( const polynomial<T,V,U>&, const polynomial<T,V,U>&);
};
