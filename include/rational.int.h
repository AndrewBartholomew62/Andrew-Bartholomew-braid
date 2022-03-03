template <class T=int> class rational
{
    T n;
    T d;

public:

	T& getn() {return n;}
	T& getd() {return d;}
	void setn(const T& t) {n = t;}
	void setd(const T& t) {d = t;}

	rational<T> (T t1=T(0), T t2=T(1)): n(t1), d(t2) {}
//	explicit rational<T> (string str, T t2=T(1)): n(str), d(t2) {}

	// default copy assignment used
    rational<T> operator += (const rational<T>);
    rational<T> operator -= (const rational<T>);
	rational<T> operator *= (const rational<T>);
    rational<T> operator /= (const rational<T>);
    rational<T> operator %= (const rational<T>);
    bool operator == (const rational<T>) const;
    bool operator > (const rational<T>) const;

    friend istream& operator >> <> (istream&, rational<T>&);
    friend ostream& operator << <> (ostream& , rational<T>);
    friend int num_len <> (const rational<T>& );
	friend rational<T> scaled_rational<>(rational<T>&);
	friend rational<T> abs <> (const rational<T>&);

	/* These are friends to allow implicit conversion of their arguments,
       helper functions do not support implicit conversion */
	friend rational<T> operator + <> (const rational<T>, const rational<T>);
	friend rational<T> operator - <> (const rational<T>, const rational<T>);
	friend rational<T> operator * <> (const rational<T>, const rational<T>);
	friend rational<T> operator / <> (const rational<T>, const rational<T>);
	friend rational<T> operator % <> (const rational<T> a, const rational<T> b);
	friend bool operator != <> (const rational<T>, const rational<T>);
	friend bool operator >= <> (const rational<T>, const rational<T>);
	friend bool operator < <> (const rational<T>, const rational<T>);
	friend bool operator <= <> (const rational<T>, const rational<T>);

};

