/**********************************************************************
  matrix.h contains a matrix template definition based on the Matrix
  class described by Strustrup.

  NOTE: This representation of matrices stores them row by row in a valarray
        rather than column by column in a valarray.  This is because I think
		a row by row approach is more natural 

  this matrix.h was derived from the array-of-pointers template version
            A. Bartholomew, 3rd January, 2005
			
***********************************************************************/
#include <valarray>
#include <Slice_iter.h>

using namespace Slice_iterators;

struct matrix_control 
{
	static unsigned int DEBUG; //bitmap
	enum parameters
	{
		general=1,
		immanant=2, 
		multiply=4, 
		inverse=8, 
		echelon=16, 
		// 'all' also supported
	};

	static bool WAIT_INFO;
	static bool COMFORT_DOTS;
	static bool SINGLE_LINE_OUTPUT;
	static int DET_DEBUG_LIMIT;
	static int wait_threshold;
	static int wait_count;
	static int reset_count; // number of times wait_count has reached wait_threshold
};

struct matrix_error {
	matrix_error (string message) {cout << "\nMatrix error!" << message << endl;}
};


/* the Scalar type St is the type of the scalar necessary to initialize a T with
   zero, so we can write e.g. T det = T(St(0));  This is necessary to overcome the
   restriction of one implicit conversion between user defined types.
*/
template<class T, class St> class matrix;

/* First declare the friend functions */

template<class T, class St> inline matrix<T,St> operator + (const matrix<T,St>& a, const matrix<T,St>& b)
{
	matrix<T,St> result = a;
	return result += b;
}

template<class T, class St> inline matrix<T,St> operator - (const matrix<T,St>& a, const matrix<T,St>& b)
{
	matrix<T,St> result = a;
	return result -= b;
}

template<class T, class St> inline matrix<T,St> operator * (const matrix<T,St>& a, const matrix<T,St>& b)
{
	matrix<T,St> result = a;
	return result *= b;
}

template<class T, class St> inline matrix<T,St> operator * (const T& t, const matrix<T,St>& b)
{
	matrix<T,St> result = b;
	return result *= t;
}

template<class T, class St> inline bool operator != (const matrix<T,St>& a, const matrix<T,St>& b)
{
	return !(a==b);
}

template<class T, class St> ostream& operator << (ostream& outstream, const matrix<T,St>& m)
{
    if (!matrix_control::SINGLE_LINE_OUTPUT)
		outstream << "\n";
		
    for (size_t i = 0; i< m.numrows(); i++ )
    {
        for (size_t j = 0; j < m.numcols(); j++)
        {
			outstream << m[i][j];
			if (i != m.numrows()-1 || j != m.numcols()-1)
				outstream << ' ';
		}

	    if (!matrix_control::SINGLE_LINE_OUTPUT)
		    outstream << endl;
    }
    return outstream;
}

/* subtract this from n x I */
template<class T, class St> matrix<T,St> operator - (const int& num, const matrix<T,St>& M)
{
	matrix<T,St> result(M.numrows(),M.numcols());
	St minus_one = -1;
	for (size_t i = 0; i< M.numrows(); i++)
	{
		for (size_t j=0; j< M.numcols(); j++)
		{
			result[i][j] = M[i][j]*T(minus_one);
		}
		
		St Mnum = num;
		result[i][i] += T(Mnum);
	}
	return result;
}

template<class T, class St> bool operator == (const matrix<T,St>& M, const int num )
{
	for (int i = 0; i< M.numrows(); i++)
	{
		for (int j=0; j< i; j++)
		{
			if (M[i][j] != T(St(0)))
				return false;
		}
		if (M[i][i] != T(St(num)))
			return false;
		for (int j=i+1; j< M.numcols(); j++)
		{
			if (M[i][j] != T(St(0)))
				return false;
		}
	}
	return true;
}

template<class T, class St> inline bool operator != (const matrix<T,St>& M, const int num )
{
	return !(M == num);
}

template<class T, class St> inline bool operator == (const int num, const matrix<T,St>& M)
{
	return M == num;
}

template<class T, class St> inline bool operator != (const int num, const matrix<T,St>& M)
{
	return !(M == num);
}

template<class T, class St> inline T trace (const matrix<T,St> M)
{
	T result;
	
	for (int i = 0; i < M.numrows(); i++)
		result += M[i][i];
	
	return result;
}

template <class T, class St> 
int non_zero_count(matrix<T,St>& mat, int row)
{
   int count = 0;
   for ( unsigned int i = 0; i < mat.numcols(); i++)
       if ( mat[row][i] != T(0))
	   count += 1;
   return (count);
}

template <class T, class St> T gcdn (const matrix<T,St>& mat, int row)
{
    T g = abs(mat[row][0]);

    for (unsigned int i = 1; i<mat.numcols(); i++)
       g = gcd(g, abs(mat[row][i]));

    return (g);
}

/* if reduced form row reduce above the diagonal as well as below */
template <class T, class St> 
void echelon (matrix<T,St>& matrixref, bool field_coefficients, bool reduced_form = false, matrix<T,St>* P_ptr = 0)
{
	matrix<T,St> mat(matrixref);
    int m = mat.numrows();
	int n = mat.numcols();

if (matrix_control::DEBUG & matrix_control::echelon)
{
	debug << "echelon: presented with mat " << endl;
	print (mat,debug,3,"echelon:   ");
	debug << "echelon: field_coefficients = " << field_coefficients << " reduced_form = " << reduced_form << endl;
}

	/* P is the product of the elementary row matrices taking matrixref to mat */
	matrix<T,St> P(m,m);
	for (int i=0; i< m; i++)
	{
	    for (int j=0; j<i; j++)
			P[i][j]=0;
		
		P[i][i]=1;

	    for (int j=i+1; j<m; j++)
			P[i][j]=0;
	}

if (matrix_control::DEBUG & matrix_control::echelon)
	debug << "echelon: num rows m = " << m << " num cols n = " << n << endl;

    /* We use a row permutation vector, so that the i-th row of mat is stored in
       the perm[i] slice of the matrix valarray.
    */
    vector<int> perm(m); 

    for (int i = 0 ; i< m; i++)
		perm[i] = i;

    int r = 0; // r and c will be the position of the next pivot
    int c = 0;
    bool complete = false;

	/* If we're not using field coefficients, we shall divide down by GCDs later, since we cannot guarantee to be able to scale P 
	   by the same amount.  All we know is that at some point we have a mat with P*mat = matrixref and mat contains a row with a 
	   common factor but that does not mean that the matrix P contains a row with the same common factor 
	*/
    if (false && field_coefficients && non_zero_count(mat,0) > 1 )
    {
		T g = gcdn (mat,0);
		
		if (g !=T(1))
		{
			for (int i = 0; i < n ; i++)
			   mat[0][i] /= g;

			for (int i = 0; i < m ; i++)
			   P[0][i] /= g;

if (matrix_control::DEBUG & matrix_control::echelon)
	debug << "echelon: row 0 scaled by " << g << endl;
		}
    }


    do
    {
		int p = r;
		bool found = false;

		/* look for the next pivot */
		do
		{
	    	do // look down the current column for a pivot
	    	{
				if (mat[perm[p]][c] != T(0))
				{						
if (matrix_control::DEBUG & matrix_control::echelon)
	debug << "echelon: found non-zero pivot entry " << mat[perm[p]][c] << " at row " << p << ", column " << c << endl;
	
					if (p !=r)
					{
						/* interchange rows r and p */
						swap(perm[r],perm[p]);

if (matrix_control::DEBUG & matrix_control::echelon)
	debug << "echelon: interchanged rows " << r << " and " << p << endl;	
					}
					
		    		found = true;
	
				}
				else
		    		p += 1;
	    	} while ( !found && p < m);

	    	if (!found) // move to the next column
	    	{
				p = r;
				c += 1;
	    	}
		} while (!found && c < n );

		if (!found)
	    	complete = true;
		else
		{
			/* having interchanged perm[r] and perm[p], the next pivot is now element (perm[r],c) */

if (matrix_control::DEBUG & matrix_control::echelon)
{
	debug << "echelon: pivot (" << r << ',' << c << ')' << endl;
	debug << "echelon: pivot row          : ";
	for (int k=0; k<n; k++)
		debug << setw(3) << mat[perm[r]][k];
	debug << endl;
}
	    	if (field_coefficients &&  non_zero_count(mat,perm[r]) > 1)
	    	{
				// divide the row by the pivot so the pivot becomes 1				
				T g = mat[perm[r]][c];
				
				if (g != T(1))
				{
					for (int i= 0; i< n ; i++)
				    	mat[perm[r]][i] /= g;

					for (int i= 0; i< m ; i++)
				    	P[perm[r]][i] /= g;

if (matrix_control::DEBUG & matrix_control::echelon)
{
	debug << "echelon: row " << r << " scaled by " << g << endl;
	debug << "echelon: updated pivot row          : ";
	for (int k=0; k<n; k++)
		debug << setw(3) << mat[perm[r]][k];
	debug << endl;
}

if (matrix_control::DEBUG & matrix_control::echelon)
{
//	debug << "echelon: matrix row reduced to" << endl;
//	perm_print (mat,debug,3,perm,"echelon:   ");
	
	matrix<T,St> mat_test(m,n);
	for (int i=0; i< m; i++)
	for (int j=0; j< n; j++)
		mat_test[i][j] = mat[perm[i]][j];

		
	matrix<T,St> P_test(m,m);
	for (int i=0; i< m; i++)
	for (int j=0; j< m; j++)
		P_test[i][j] = P[perm[i]][j];

//	debug << "echelon: elementary matrix" << endl;
//	print (P_test,debug,3,"echelon:   ");
		
	matrix<T,St> PM = P_test*matrixref;
	
	if (PM != mat_test)
	{
		debug << "echelon: PM does not match mat " << endl;
		print (PM,debug,3,"echelon:   ");

		cout << "elementary matrix product does not match!" << endl;
		exit(0);
	}

	T det = determinant(P_test);
	debug << "echelon:   det(P_test) = " << det << endl;
		
}			

				}
	    	}
			
	    	for (int i = (reduced_form? 0:r+1) ; i < m; i++)
			{
				if (i == r)
				{
if (matrix_control::DEBUG & matrix_control::echelon)
	debug << "echelon:   omitting row " << i << endl;
					continue; //added for reduced echelon form
				}
					
				if (mat[perm[i]][c] != T(0))
				{
if (matrix_control::DEBUG & matrix_control::echelon)
	debug << "echelon:   procesing row " << i << endl;

		    		if(false && field_coefficients && non_zero_count(mat,perm[i]) > 1)
		    		{
						T g = mat[perm[i]][c];
						
						if (g != T(1))
						{
							for (int j = 0 ;j < n ; j++)
								mat[perm[i]][j] /= g;

							for (int j= 0; j< m ; j++)
						    	P[perm[i]][j] /= g;
								
if (matrix_control::DEBUG & matrix_control::echelon)
	debug << "echelon:   row " << i << " scaled by " << g << endl;
						}
		    		}

		    		T a = mat[perm[r]][c];
		    		T b = mat[perm[i]][c];

if (matrix_control::DEBUG & matrix_control::echelon)
{
	debug << "echelon:   replace row " << i << "      : ";
	for (int k=0; k<n; k++)
		debug << setw(3) << mat[perm[i]][k];
	debug << endl;
	debug << "echelon:   with a = mat(" << r <<',' << c << ") = " << a << " times row " << i << " : ";
	for (int k=0; k<n; k++)
		debug << setw(3) << a*mat[perm[i]][k];
	debug << endl;
	debug << "echelon:   minus b = mat(" << i <<',' << c << ") = " << b << " times row " << r << ": ";
	for (int k=0; k<n; k++)
		debug << setw(3) << b*mat[perm[r]][k];
	debug << endl;	
}

					/* We have to adjust the whole row of mat because we have added the reduced form, so 
					   cannot just adjust from column c as we did before: for (int j = c ; j< n; j++)
					*/			    		
		    		for (int j = 0 ; j< n; j++)
		    		   mat[perm[i]][j] = a * mat[perm[i]][j] - b * mat[perm[r]][j];

		    		for (int j = 0 ; j< m; j++)
		    		   P[perm[i]][j] = a * P[perm[i]][j] - b * P[perm[r]][j];
		    		 
if (matrix_control::DEBUG & matrix_control::echelon)
{
	debug << "echelon:   giving             : ";
	for (int k=0; k<n; k++)
		debug << setw(3) << mat[perm[i]][k];
	debug << endl;
}

if (matrix_control::DEBUG & matrix_control::echelon)
{
//	debug << "echelon: after single row operation matrix row reduced to" << endl;
//	perm_print (mat,debug,3,perm,"echelon:   ");
	
	matrix<T,St> mat_test(m,n);
	for (int i=0; i< m; i++)
	for (int j=0; j< n; j++)
		mat_test[i][j] = mat[perm[i]][j];
		
	matrix<T,St> P_test(m,m);
	for (int i=0; i< m; i++)
	for (int j=0; j< m; j++)
		P_test[i][j] = P[perm[i]][j];

//	debug << "echelon: elementary matrix" << endl;
//		print (P_test,debug,3,"echelon:   ");
		
	matrix<T,St> PM = P_test*matrixref;
	
	if (PM != mat_test)
	{
		debug << "echelon: PM does not match mat " << endl;
		print (PM,debug,3,"echelon:   ");
		cout << "elementary matrix product does not match!" << endl;
		exit(0);
	}	

	T det = determinant(P_test);
	debug << "echelon:   det(P_test) = " << det << endl; //<< ", size " << sizeof(det) << endl;
	
}


				}
			}
		} // end of processing the pivot

if (matrix_control::DEBUG & matrix_control::echelon)
{
	debug << "echelon: matrix after processing pivot (" << r << ',' << c << ')' << endl;
	matrix<T,St> mat_test(m,n);
	for (int i=0; i< m; i++)
	for (int j=0; j< n; j++)
		mat_test[i][j] = mat[perm[i]][j];
	print (mat_test,debug,3,"echelon:   ");
}

		if (c == n-1)
	    	complete = true;
		else  // move down the diagonal
		{
	    	c += 1;
	    	r += 1;
		}
    } while (!complete && r < m);

	/* final test */
	matrix<T,St> mat_test(m,n);
	for (int i=0; i< m; i++)
	for (int j=0; j< n; j++)
		mat_test[i][j] = mat[perm[i]][j];

	matrix<T,St> P_test(m,m);
	for (int i=0; i< m; i++)
	for (int j=0; j< m; j++)
		P_test[i][j] = P[perm[i]][j];

	matrix<T,St> PM = P_test*matrixref;

if (matrix_control::DEBUG & matrix_control::echelon)
{
//	debug << "echelon: final P " << endl;
//	print (P_test,debug,3,"echelon:   ");

	debug << "echelon: echelon form " << endl;
	print (mat_test,debug,3,"echelon:   ");
}
		
	if (PM != mat_test)
	{

if (matrix_control::DEBUG & matrix_control::echelon)
{
	debug << "echelon: PM does not match mat " << endl;
	print (PM,debug,3,"echelon:   ");
}	
		cout << "elementary matrix product does not match!" << endl;
		exit(0);
	}
	else
	{
if (matrix_control::DEBUG & matrix_control::general)
{
	debug << "echelon: elementary matrix product check matches calculated echelon form" << endl;
}	
	}		

	/* remove common factors from mat before undoing the row permutation into matrixref */
	bool scaled_result = false;
	
	for (int i=0 ; i< m ; i++)
	{	
    	if (!field_coefficients && non_zero_count(mat,perm[i]) > 1)
    	{
			T g = gcdn (mat,perm[i]);
			
			if (g !=T(1))
			{
				for (int j=i; j< n ; j++) // mat is in echelon form so mat[perm[i]][j] = 0 for j=0,...,i-1
			    	mat[perm[i]][j] /= g;
			
				scaled_result = true;
				
if (matrix_control::DEBUG & matrix_control::general)
	debug << "echelon: before returning, row " << i << " scaled by " << g << endl;
			}
    	}
				
	    for (int j = 0; j<n; j++)
			matrixref[i][j] = mat[perm[i]][j];
	}

	if (!scaled_result)
	{
if (matrix_control::DEBUG & matrix_control::general)
	debug << "echelon: no row scaling before returning" << endl;
	}

	if (P_ptr != 0)
		*P_ptr = P_test;	
}

/* this inverse function is based on echelon reduction */
template <class T, class St> 
matrix<T,St> inverse (const matrix<T,St>& M, bool field_coefficients = false)
{
	int n=M.numrows();


	
	matrix<T,St> MI(n,2*n);
	
	for (int i=0; i<n; i++)
	for (int j=0; j<n; j++)
		MI[i][j] = M[i][j];


	matrix<T,St> I(n,n);
	
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<i; j++)
			I[i][j] = 0;

		I[i][i] = T(1);

		for (int j=i+1; j<n; j++)
			I[i][j] = 0;
	}

	for (int i=0; i<n; i++)
	for (int j=0; j<n; j++)
		MI[i][n+j] = I[i][j];

if (matrix_control::DEBUG & matrix_control::inverse)
{
	debug << "inverse: initial MI:" << endl;
	print(MI, debug, 3, "inverse: ");
}
	echelon(MI,field_coefficients,true);

if (matrix_control::DEBUG & matrix_control::inverse)
{
	debug << "inverse: ecelon returned MI:" << endl;
	print(MI, debug, 3, "inverse: ");
}

	/* invert any rows that have -1 in the leading diagonal */
	for (int i=0; i<n; i++)
	{
		if (MI[i][i] == T(-1))
		{			
			for (int j=i; j<2*n; j++)
				MI[i][j] *= T(-1);
		}
	}
if (matrix_control::DEBUG & matrix_control::inverse)
{
	debug << "inverse: final MI:" << endl;
	print(MI, debug, 3, "inverse: ");
}
	

	matrix<T,St> Minv(n,n);


	for (int i=0; i<n; i++)
	for (int j=0; j<n; j++)
		Minv[i][j] = MI[i][n+j];

if (matrix_control::DEBUG & matrix_control::inverse)
{
	debug << "inverse: Minv:" << endl;
	print(Minv, debug, 3, "inverse: ");
}

	
	matrix<T,St> check = M*Minv;

	if (check != I)
	{

if (matrix_control::DEBUG & matrix_control::inverse)
{
	debug << "inverse: inverse check..." << endl;
	print(check, debug, 3, "inverse: ");
}
		cout << "\nincorrect inverse!" << endl;
		exit(0);
	}

if (matrix_control::DEBUG & matrix_control::inverse)
	debug << "OK" << endl;
	
if (matrix_control::DEBUG & matrix_control::general)
	debug << "inverse: inverse check OK" << endl;

	return Minv;
}

/* determinant uses a row and column permutation rperm and cperm, and the number of rows/columns n, 
   to evaluate the determinant of any sub-matrix of the square matrix M.  The recursion_level is used
   only for debugging purposes.

   The approach is to look along the top row and down the left column of the submatrix determined by rperm 
   and cperm, and to evaluate the determinant based on the smallest number of non-zero entries.  The function 
   then recurses to evaluate the sub-determinants; only when we get to n = 2 do we have a 2x2 submatrix 
   whereupon we calculate the determinant directly.
*/

template <class T, class St> 
T immanant (const matrix<T,St>& M, string title, int n, vector<int> rperm, vector<int> cperm, bool permanent, int recursion_level=0)
{
    T zero = T(St(0));
    T det = zero;
    T temp = zero;

	/* We are starting a new immanant calculation if recursion_level == 0 */
	if (!recursion_level)
	{
		matrix_control::reset_count = 0;

if (matrix_control::DEBUG & matrix_control::immanant)
{
	debug << "matrix::immanant: underlying matrix M = \n";
	print (M,debug,0,"matrix::immanant: ");
	debug << "matrix::immanant: n = " << n << endl;
	debug << "matrix::immanant: matrix_control::WAIT_INFO = " << (matrix_control::WAIT_INFO?"true":"false") << endl;
	debug << "matrix::immanant: matrix_control::COMFORT_DOTS = " << (matrix_control::COMFORT_DOTS?"true":"false") << endl;
	debug << "matrix::immanant: matrix_control::wait_threshold = " << matrix_control::wait_threshold << endl;
}
	}
	
	if (rperm.size() == 0)
	{
	
if (matrix_control::DEBUG & matrix_control::immanant)
	debug << "matrix::immanant: default parameters provided, creating permutations" << endl;
	
		n = M.numcols();
		rperm = vector<int>(M.numrows());
		for (size_t i=0; i< M.numrows(); i++)
			rperm[i] = i;
		cperm = vector<int>(n);
		for (int i=0; i< n; i++)
			cperm[i] = i;
	}


if (matrix_control::DEBUG & matrix_control::immanant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::immanant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "n = " << n << ", recursion_level = " << recursion_level << endl;

    debug << "matrix::immanant: ";
	for (int k=0;k<recursion_level;k++)
	    debug << "    ";
	debug << "rperm: ";
	for ( int k = 0 ; k< n; k++)
    	debug << rperm[k] << " ";
	debug << endl;

    debug << "matrix::immanant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "cperm: ";
    for (int k = 0 ; k< n; k++)
		debug << cperm[k] << " ";
	debug << endl;
}

	if (n==1)
	{
if (matrix_control::DEBUG & matrix_control::immanant && n >= matrix_control::DET_DEBUG_LIMIT)
	debug << "matrix::immanant: returning single element " << M[rperm[0]][cperm[0]] << endl;

		return M[rperm[0]][cperm[0]];
	}
    else if (n == 2)
    {
if (matrix_control::DEBUG & matrix_control::immanant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::immanant: ";
	for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "positive term: [" << rperm[0] << "][" << cperm[0] << "]*["
          << rperm[1] << "][" << cperm[1] << "] = " 
	      << M[rperm[0]][cperm[0]] * M[rperm[1]][cperm[1]];
	debug << endl;
    debug << "matrix::immanant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "negative term: [" << rperm[0] << "][" << cperm[1] << "]*["
          << rperm[1] << "][" << cperm[0] << "] = "
	      << M[rperm[0]][cperm[1]] * M[rperm[1]][cperm[0]];
	debug << endl;	
}

		if (permanent)
		{
			det =   M[rperm[0]][cperm[0]] * M[rperm[1]][cperm[1]]
			      + M[rperm[0]][cperm[1]] * M[rperm[1]][cperm[0]];
		}
		else
		{
			det =   M[rperm[0]][cperm[0]] * M[rperm[1]][cperm[1]]
			      - M[rperm[0]][cperm[1]] * M[rperm[1]][cperm[0]];
		}

    }
    else
    {
		vector<int> sub_r_perm(n-1);
		vector<int> sub_c_perm(n-1);
		bool evaluate_along_row = true;
		
		/* Look to see whether there are more zeros along the top row or left column
		   we evaluate the immanant based on the largest number of zero values
		*/
		int num_zeros = 0;
		for (int i=0; i < n; i++)
		{
			if (M[rperm[0]][cperm[i]] == zero)
				num_zeros++;
				
			if (M[rperm[i]][cperm[0]] == zero)
				num_zeros--;
		}
				
		if (num_zeros < 0)
			evaluate_along_row = false;

if (matrix_control::DEBUG & matrix_control::immanant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::immanant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "evaluating " << n << "-sub-immanant " << (evaluate_along_row? "along row" : "down column") 
	      << ", row-column num_zeros = " << num_zeros << endl;
}

		if (evaluate_along_row) // row permutation for all recursive calls at this level is the same
		{
			for (int i=1; i<n;i++)
	    		sub_r_perm[i-1] = rperm[i];

if (matrix_control::DEBUG & matrix_control::immanant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::immanant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "sub_r_perm (fixed for iterative calls at this level): ";
	for ( int k = 0 ; k< n-1; k++)
	   	debug << sub_r_perm[k] << " ";
	debug << endl;
}			
		}
		else // column permutation for all recursive calls at this level is the same
		{
			for (int i=1; i<n;i++)
	    		sub_c_perm[i-1] = cperm[i];
			
if (matrix_control::DEBUG & matrix_control::immanant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::immanant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "sub_c_perm (fixed for iterative calls at this level): ";
	for ( int k = 0 ; k< n-1; k++)
	   	debug << sub_c_perm[k] << " ";
	debug << endl;
}			
		}
				
		for (int i=0; i<n; i++)
		{
		
		    if ((evaluate_along_row && M[rperm[0]][cperm[i]] != T(St(0))) || (!evaluate_along_row && M[rperm[i]][cperm[0]] != T(St(0))))
	    	{
			
				if (evaluate_along_row)
				{
					for (int j=0; j<i;j++)
			    		sub_c_perm[j] = cperm[j];
					for (int j=i+1;j<n;j++)
			    		sub_c_perm[j-1] = cperm[j];

if (matrix_control::DEBUG & matrix_control::immanant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::immanant: ";
	for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << n-1 << "-multiplier =" << (i%2? " -1 * " :" ") << M[rperm[0]][cperm[i]] << endl;
}

					temp = M[rperm[0]][cperm[i]] * immanant(M,title, n-1, sub_r_perm, sub_c_perm, permanent, recursion_level+1);

if (matrix_control::DEBUG & matrix_control::immanant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::immanant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "additive term = " << temp << endl;
}

				}
				else
				{
					for (int j=0; j<i;j++)
			    		sub_r_perm[j] = rperm[j];
					for (int j=i+1;j<n;j++)
			    		sub_r_perm[j-1] = rperm[j];

if (matrix_control::DEBUG & matrix_control::immanant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::immanant: ";
	for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << n-1 << "-multiplier =" << (i%2? " -1 * " :" ") << M[rperm[i]][cperm[0]] << endl;
}

					temp = M[rperm[i]][cperm[0]] * immanant(M,title, n-1, sub_r_perm, sub_c_perm, permanent, recursion_level+1);

if (matrix_control::DEBUG & matrix_control::immanant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::immanant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << "additive term = " << temp << endl;
}
				}
				

//				if (i%2)
				if (!permanent && i%2)
					det -= temp;
				else
					det += temp;
	
if (matrix_control::DEBUG & matrix_control::immanant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::immanant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << n << "-sub-immanant: " << det << endl;
}

	    	}
		}
    }

	if (matrix_control::WAIT_INFO)
	{
	    if (n == matrix_control::wait_threshold)
	    {
			if (matrix_control::COMFORT_DOTS)
			{
				cout << ".";
   				cout.flush();
			}
			
			if (matrix_control::wait_count > 1000)
			{
				matrix_control::reset_count++;
	    		cout << "\nworking on " << (permanent? "permanent": "determinant") << " for " << title << " (" << matrix_control::reset_count << "), please wait\n";
	    		matrix_control::wait_count = 0;
			}
			else
	    		matrix_control::wait_count++;
    	}
	}

if (matrix_control::DEBUG & matrix_control::immanant && n >= matrix_control::DET_DEBUG_LIMIT)
{
    debug << "matrix::immanant: ";
    for (int k=0;k<recursion_level;k++)
		debug << "    ";
    debug << n << "-immanant: " << det << endl;
}

    return det;
}

template <class T, class St> 
T determinant (const matrix<T,St>& M, string title="untitled", int n=0, vector<int> rperm=vector<int>(0), vector<int> cperm=vector<int>(0))
{
	return immanant (M, title, n, rperm, cperm, false); // permanent = false
}

template <class T, class St> 
T permanent (const matrix<T,St>& M, string title="untitled", int n=0, vector<int> rperm=vector<int>(0), vector<int> cperm=vector<int>(0))
{
	return immanant (M, title, n, rperm, cperm, true); // permanent = true
}

#include<matrix.int.h>


template<class T, class St> matrix<T,St>::matrix (size_t r, size_t c)
{
    v = new valarray<T>(T(St(0)),r*c);
    rows = r;
    cols = c;
}

template<class T, class St> matrix<T,St>::matrix (size_t r, size_t c, T t)
{
    v = new valarray<T>(t,r*c);
    rows = r;
    cols = c;
}

/* the copy constructor */
template<class T, class St> matrix<T,St>::matrix (const matrix<T,St>& M)
{
	cols = M.cols;
	rows = M.rows;
	
	if (M.v)
	{
	    v = new valarray<T>(T(St(0)),rows*cols);
		*v = *M.v;
	}
	else
		v = 0;
}

/* row and column definitions - stored row by row in the representation */
template <class T, class St> inline Slice_iter<T> matrix<T,St>::row(size_t i)
{
	return Slice_iter<T>(v,slice(i*cols, cols, 1));
}

template <class T, class St> inline Cslice_iter<T> matrix<T,St>::row(size_t i) const
{
	return Cslice_iter<T>(v,slice(i*cols, cols, 1));
}

template <class T, class St> inline Slice_iter<T> matrix<T,St>::column(size_t i)
{
	return Slice_iter<T>(v,slice(i,rows,cols));
}

template <class T, class St> inline Cslice_iter<T> matrix<T,St>::column(size_t i) const
{
	return Cslice_iter<T>(v,slice(i,rows,cols));
}

/* copy assignment */
template<class T, class St> inline matrix<T,St>& matrix<T,St>::operator = (const matrix<T,St>& M)
{
	if (this != &M)
	{	
		cols = M.cols;
		rows = M.rows;
		
	    delete v;
	    if (M.v)
	    {
			v = new valarray<T>(T(St(0)),rows*cols);
			*v = *M.v;
		}
		else
			v=0;
		
	}
	return *this;
}


/* Now the operators */
template<class T, class St> matrix<T,St> matrix<T,St>::operator += (const matrix<T,St>& M)
{
	matrix<T,St>& loc = *this;
	
	if (loc.rows != M.rows || loc.cols != M.cols)
		throw(matrix_error ("incompatible sizes in += operator"));
	
	for (size_t i=0; i< loc.numrows(); i++)
	{
		for (size_t j=0; j< loc.numcols(); j++)
		{
			loc[i][j] += M[i][j];
		}
	}
	
	return loc;
}

template<class T, class St> matrix<T,St> matrix<T,St>::operator -= (const matrix<T,St>& M)
{
	matrix<T,St>& loc = *this;
	
	if (loc.rows != M.rows || loc.cols != M.cols)
		throw(matrix_error ("incompatible sizes in -= operator"));
	
	for (size_t i=0; i< loc.numrows(); i++)
	{
		for (size_t j=0; j< loc.numcols(); j++)
		{
			loc[i][j] -= M[i][j];
		}
	}
	
	return loc;
}

template<class T, class St> matrix<T,St> matrix<T,St>::operator *= (const matrix<T,St>& M)
{
	matrix<T,St>& loc = *this;

	if (loc.cols != M.rows )
	{
if (matrix_control::DEBUG & matrix_control::multiply)
	debug << "matrix::operator *= : loc.cols = " << loc.cols <<  "M.rows = " << M.rows << endl;
	
		throw(matrix_error ("incompatible sizes in *= operator"));
	}

	matrix<T,St> result(loc.rows, M.cols);

if (matrix_control::DEBUG & matrix_control::multiply)
{
	debug << "matrix::operator *= : \n\t(*this) = " << *this <<  "\n\t      M = " << M << endl;
	debug << "matrix::operator *= : multiplication produces " << loc.numrows() << " by " << M.numcols() << " result" << endl;
}


	for (size_t i=0; i< loc.rows; i++)
	{
		for (size_t j=0; j< M.cols; j++)
		{
if (matrix_control::DEBUG & matrix_control::multiply)
	debug << "matrix::operator *= : element " << i << "," << j <<" contributions" << endl;	
			for (size_t k=0; k< loc.cols; k++)
			{

if (matrix_control::DEBUG & matrix_control::multiply)
	debug << "matrix::operator *= : result[" << i << "][" << j <<"] stands at " << result[i][j] << endl;	
				result[i][j] += loc[i][k]*M[k][j];
if (matrix_control::DEBUG & matrix_control::multiply)
{
	T temp = loc[i][k]*M[k][j];
	debug << "matrix::operator *= :   " << temp;
}
			}
if (matrix_control::DEBUG & matrix_control::multiply)
	debug << endl;
		}
	}
	
	loc = result;
	return loc;
}


template<class T, class St> matrix<T,St> matrix<T,St>::operator *= (const T& t)
{
	matrix<T,St>& loc = *this;
	
	for (unsigned int i=0; i< loc.numrows(); i++)
	for (unsigned int j=0; j< loc.numcols(); j++)
		loc[i][j] *= t;

	return loc;
}


/* subtract n x I from this */
template<class T, class St> matrix<T,St> matrix<T,St>::operator -= (const int num)
{
	matrix<T,St>& loc = *this;
	matrix<T,St> result = loc;
	for (int i=0; i< loc.numrows(); i++)
	{
		result[i][i] -= T(1);
	}
	return result;
}

/* Currently this inverse function assumes that the determinant is non-zero */
template<class T, class St> matrix<T,St> matrix<T,St>::inverse(bool adjunct_only) const
{

	const matrix<T,St>& M = *this;
	int n = M.numrows();

cout << "calculating inverse of " << n << " by " << n << " matrix" << endl;
	
	matrix<T,St> invM(n,n);
	
	
	T det = determinant (M);

if (matrix_control::DEBUG & matrix_control::inverse)
	debug << "\nmatrix::inverse: det(M) = " << det << endl;

	vector<int> rperm(n);
	vector<int> cperm(n);
	for (int i=0; i<n; i++)
		rperm[i]=cperm[i]=i;

	if (det == T(St(0)))
		throw(matrix_error ("attempt to take inverse of singular matrix"));
	
	/* calculate the inverse of M using the adjoint method */
	for (int i = 0; i < n; i++)
	{
		/* take out row i */
		for (int k=0; k < i; k++)
			rperm[k] = k;
		for (int k=i+1; k<n; k++)
			rperm[k-1] = k;
	
		for (int j=0; j<n; j++)
		{
//cout << '#' << flush;			
			/* take out column j */
			for (int k=0; k < j; k++)
				cperm[k] = k;
			for (int k=j+1; k<n; k++)
				cperm[k-1] = k;
			
			/* element (j,i) of the adjoint is the signed n-1 x n-1 determinant of the 
			   matrix determined by rperm and cperm */
			invM[j][i] = determinant(M,"",n-1,rperm,cperm);

if (matrix_control::DEBUG & matrix_control::inverse)
	debug << "\nmatrix::inverse: initial invM[" << j << "][" << i << "] = " << invM[j][i] << endl;

			/* the sign is given by {-1}^{i+j} */
			St minus_one = St(-1);
			if ((i+j) % 2)
				invM[j][i] *= T(minus_one);

if (matrix_control::DEBUG & matrix_control::inverse)
	debug << "\nmatrix::inverse: after sign invM[" << j << "][" << i << "] = " << invM[j][i] << endl;
	
		}
	}
	
	/* divide down by the determinant to get the inverse unless we're
	   calculating the adjunct only
	*/
	if (!adjunct_only)
	{
		for (int i=0;i<n;i++)
		for (int j=0;j<n;j++)
			invM[i][j] /= det;
	}
	
	return invM;
}

template<class T, class St> bool matrix<T,St>::operator == (const matrix<T,St>& M) const
{
	const matrix<T,St>& loc = *this;
	
	if (loc.rows != M.rows || loc.cols != M.cols)
		return false;
	
	for (size_t i=0; i< loc.rows; i++)
	for (size_t j=0; j< loc.cols; j++)
	{
		if (loc[i][j] != M[i][j])
		{
//cout << "matrix elements " << i << ',' << j << " differ:" << " loc[i][j]  = " << loc[i][j] << ", M[i][j] = " << M[i][j] << endl;		
			return false;
		}
	}
	
	return true;
}

template<class T, class St> void matrix<T,St>::dump(ostream& os) const
{
	os << "\nvalarray at " << v << " rows = " << rows << " cols = " << cols << endl;
	for (size_t r = 0; r < rows ; r++)
	{
		os << "\nrow " << r;
		for (size_t c = 0; c < cols; c++)
		{
			os << "\n\t[" << c << "] at " << &((*this)[r][c]) << " value: " << (*this)[r][c];
		}
	}
	os << endl;
}

/* The next two functions are designed for r*N x c*N matrices that may be regarded as r x c matrices having 
   N x N matrix elements.  The functions set an N x N element in a destination matrix from an N x N element
   in a source matrix, and decrement an N x N element by 1, that is subtract the N x N identity matrix
*/
template <class T, class St> 
void set_matrix_N_element(matrix<T,St>& d_matrix, int d_N_row, int d_N_col, 
                          const matrix<T,St>& s_matrix, int s_N_row, int s_N_col, int N)
{
	int d_r_base = N*d_N_row;
	int d_c_base = N*d_N_col;
	int s_r_base = N*s_N_row;
	int s_c_base = N*s_N_col;
	
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
		d_matrix[d_r_base+i][d_c_base+j] = s_matrix[s_r_base+i][s_c_base+j];
}


template <class T, class St> 
void decrement_matrix_N_element(matrix<T,St>& d_matrix, int d_N_row, int d_N_col, int N)
{
	int d_r_base = N*d_N_row;
	int d_c_base = N*d_N_col;

	for (int i = 0; i < N; i++)
		d_matrix[d_r_base+i][d_c_base+i] += T("-1");
}

template <class T, class St> 
void print(const matrix<T,St>& m, ostream& s, int n, string prefix)
{
	for (unsigned int i=0; i< m.numrows(); i++)
		print (m[i],s,n,prefix);
}
