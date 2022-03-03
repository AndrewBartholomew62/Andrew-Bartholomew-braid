template <class T=int> class quaternion
{
    T	     r;
    T	     i;
    T	     j;
    T	     k;

public:
	quaternion(T t1=T(0), T t2=T(0), T t3=T(0), T t4=T(0)):
	r(t1), i(t2), j(t3), k(t4) {}
	
	//default copy assignment used
	
	T& Qr() {return r;}
	T& Qi() {return i;}
	T& Qj() {return j;}
	T& Qk() {return k;}
	T QrC() const {return r;}
	T QiC() const {return i;}
	T QjC() const {return j;}
	T QkC() const {return k;}
	
	
    quaternion<T> operator += (const quaternion<T>);
    quaternion<T> operator -= (const quaternion<T>);
	quaternion<T> operator *= (const quaternion<T>);
    quaternion<T> operator /= (const quaternion<T>);
    quaternion<T> operator %= (const quaternion<T>);

    bool operator == (const quaternion<T>) const;

	static void quaternion_error()	{ cout << "\n\nError in quaternions!"; exit(0);}
    friend istream& operator >> <> (istream&, quaternion<T>&);
    friend ostream& operator << <> (ostream& , const quaternion<T>);
    friend int num_len <> (const quaternion<T> );
    friend quaternion<T> conj <> (const quaternion<T>);
	friend T norm_sqrd <> (const quaternion<T>);

	/* These are friends to allow implicit conversion of their arguments,
       helper functions do not support implicit conversion */
	friend quaternion<T> operator +<> (const quaternion<T>, const quaternion<T>);
	friend quaternion<T> operator -<> (const quaternion<T>, const quaternion<T>);
	friend quaternion<T> operator *<> (const quaternion<T>, const quaternion<T>);
	friend quaternion<T> operator /<> (const quaternion<T>, const quaternion<T>);
	friend bool operator !=<> (const quaternion<T>, const quaternion<T>);
};
