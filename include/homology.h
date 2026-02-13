/*************************************************************************************************
                  Generic homology calculation templates

                  A. Bartholomew 2rd June, 2025
                  
The original version of the Smith normal form function appeared in braid version 11.0 in 2006 but was subsequently 
removed as it was unused at the time.

The current version is restricted to integer matrices and does not support the row and column 
permutation provided by the original.  However, this implementation supports non-square matrices

template <class T, class St> void test_matrices (matrix<T,St>& D, matrix<T,St>& P, matrix<T,St>& Q, const matrix<T,St> A, string test_case)
template <class T, class St> int Smith_normal_form (const matrix<T,St> A, matrix<T,St>& D, matrix<T,St>& P, matrix<T,St>& Q,  bool field_coefficients = false)
template <class T, class St> matrix<T,St> SNF_inverse (const matrix<T,St>& M, bool field_coefficients = false)
template <class T, class St> homology_generators_return<T> homology_generators(matrix<T,St>& Delta_k, matrix<T,St>& Delta_kp1, bool cohomology, bool field_coefficients)	
*************************************************************************************************/

class homology_control
{
public:
	bool field_coefficients;
	bool silent_operation;
	bool test_Smith_normal_form;

	homology_control(): field_coefficients(false), silent_operation(false),test_Smith_normal_form(false) {}
};


template <class T> class homology_generators_return
{
public:
	int num_generators;
	int num_torsion_generators;
	list<vector<T> > generators;
	vector<T> torsion;

	homology_generators_return(): num_generators(0), num_torsion_generators(0) {}
};

template <class T, class St> void test_matrices (matrix<T,St>& D, matrix<T,St>& P, matrix<T,St>& Q, const matrix<T,St> A, string test_case)
{


if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "test_matrices: " << test_case;


if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "\ntest_matrices: P" << endl;
	print(P, debug, 3, "test_matrices:   ");
	debug << "test_matrices: Q" << endl;
	print(Q, debug, 3, "test_matrices:   ");
	debug << "test_matrices: D" << endl;
	print(D, debug, 3, "test_matrices:   ");
	debug << "test_matrices: " << test_case;
}
	
	// Now work out PAQ
	matrix<T,St> PAQ;
	PAQ = (P*A);

/*
if (false && debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "test_matrices: PA = " << endl;
	print(PAQ, debug, 3, "test_matrices:   ");
}
*/
	
	PAQ *= Q;

	// compare PAQ with D
	if (PAQ != D)
	{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "found PAQ !=D" << endl;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "\ntest_matrices: PAQ" << endl;
	print(PAQ, debug, 3, "test_matrices:   ");
}
	    cout << test_case << endl;
		cout << "Smith normal form found PAQ !=D" << endl;
		debug.flush();
		exit(0);
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "OK" << endl;
}


/* Smith_normal_form reduces the matrix A to Smith normal form, assigning to P and Q matrices such that
   PAQ = D, where D is diagonal.  These matrices need to be of the correct size at the call to this function.
   
   The function returns the rank of D.   
*/
template <class T, class St> int Smith_normal_form (const matrix<T,St> A, matrix<T,St>& D, matrix<T,St>& P, matrix<T,St>& Q,  homology_control h_control)														
{
	
	int n=A.numrows();
	int m=A.numcols();

	bool field_coefficients = h_control.field_coefficients;
	bool silent_operation=h_control.silent_operation;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "Smith_normal_form: calculating Smith normal form " << n << " rows	" << m << " columns, field_coefficients = " << field_coefficients << endl;

	if (!silent_operation)
		cout << "Calculating Smith normal form " << n << " rows " << m << " columns" << endl;

	for (int i=0; i<n; i++)
	{
	    for (int j=0; j<i; j++)
			P[i][j]=T(0);
		
		P[i][i]=T(1);

	    for (int j=i+1; j<n; j++)
			P[i][j]=T(0);
	}

	for (int i=0; i<m; i++)
	{
	    for (int j=0; j<i; j++)
			Q[i][j]=T(0);
		
		Q[i][i]=T(1);

	    for (int j=i+1; j<m; j++)
			Q[i][j]=T(0);
	}
		
    int rank = 0;
    	
    int pivot_r=0; // pivot_r and pivot_c will be the location of the next pivot.
    int pivot_c=0;
    
    bool complete = false;
	do
	{

		if (!silent_operation)
			cout << '.' << flush;

		if (D[pivot_r][pivot_c] == T(0) )
		{
			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: looking for non-zero pivot for (" << pivot_r << " " << pivot_c << ")" << endl;

			/* Try to move a nonzero element into position (pivot_r pivot_c), starting from r=pivot_r and c=pivot_c, we look
			   along row r, then down column c for a non-zero entry and if none is found, we increment
			   both r and c.  If we find a non-zero entry, we swap column i with column c and then swap
			   row i with row r.
			*/
			
			bool found_pivot = false;
			int step_r=0; // the number of row steps from r if nothing is found on the current row
			int step_c=0; // the number of column steps from c if nothing is found on the current column
			int r=pivot_r;
			int c;
			
			do
			{
				/* Look along row pivot_r from column pivot_c+step_c */
				for (c=pivot_c+step_c+1; c< m; c++)
				{
					if (D[r][c] != T(0))
					{							
						found_pivot = true; 
						break;
					}
				}
			
				if (!found_pivot)
				{										
					/* Look down column pivot_c+step_c */
					c = pivot_c+step_c;
					for (r=pivot_r+step_r+1; r< n; r++)
					{
						if (D[r][c] != T(0))
						{
							found_pivot = true;
							break;
						}
					}
				}
				
				if (!found_pivot)
				{
					step_c++;						
					step_r++;

					r=pivot_r+step_r;
					c=pivot_c+step_c;
					
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: moving down the diagonal to row " << r << " column " << c << endl;
	
					if (r == n || c == m)
					{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: no more pivots" << endl;
						break;
					}
					else
					{											
						if (D[r][c] != T(0))
						{
							found_pivot = true;
							break;
						}
					}					
				}				
			} while (!found_pivot);
			
			if (found_pivot)
			{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: found non-zero element in row " << r << " column " << c << endl;

				/* swap rows and columns to move pivot found at (r,c) into position (pivot_r,pivot_c) */				
				if (c != pivot_c)
				{
					/* interchange column pivot_c and column c in D and Q */

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: swapping columns " << pivot_c << " and " << c << endl;

					for (int j=0; j<n; j++)
						swap(D[j][pivot_c],D[j][c]);

					for (int j=0; j<m; j++)
						swap(Q[j][pivot_c],Q[j][c]);
				}
				
				if (r != pivot_r)
				{
					/* interchange row pivot_r and row r in D and P */
					
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: swapping rows " << pivot_r << " and " << r << endl;
	
							for (int j=0; j<m; j++)
								swap(D[pivot_r][j],D[r][j]);

							for (int j=0; j<n; j++)
								swap(P[pivot_r][j],P[r][j]);
				}
				
/*####################################################################*/
//test_matrices(D, P, Q, A,"Testing row and column interchange...");
/*####################################################################*/
			}
			else
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: couldn't find a non-zero element for (" << pivot_r << " " << pivot_c << ")" << endl;
				break;
			}									
		}
		
		
		rank++;
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: looking for non-zero entries around pivot (" << pivot_r << " " << pivot_c << ") = " << D[pivot_r][pivot_c] << endl;
	print (D, debug, 3, "Smith_normal_form: ");
}

		/* Now (pivot_r pivot_c) is non-zero.  If we're using field coefficients, divide down the row so that the pivot becomes 1. */
		if (field_coefficients)
		{		
			T factor = D[pivot_r][pivot_c];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: field coefficients, dividing row for pivot (" << pivot_r << " " << pivot_c << ") = " << factor << endl;		
			
			if (factor != T(1))
			{
				for (int j=pivot_c; j< m; j++)
					D[pivot_r][j] /= factor;
	
				/* collect the contribution to P, we are multiplying P on the left by the elementary matrix that divides row pivot_r by factor	*/
				for (int k=0; k<n; k++)
					P[pivot_r][k] /= factor;
									
/*####################################################################*/
//test_matrices(D, P, Q, A,"Testing row division by pivot value...");
/*####################################################################*/
			}
			else
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: row division not required" << endl;				
			}
		}
		
		
		
		/* Reduce along the row and down the column so the remainder of the ith row and ith column is zero.	
		   
		   The row operations carried out as a result of looking down column pivot_c may reduce the value of D[pivot_r][pivot_c], in which case we
		   process the row and column again, as per Cohen's algorithm.
		*/
		bool again = false;  // first time round do row and column loop
		do
		{
			/* look along the pivot_r row */

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: processing row " << pivot_r << " around pivot (" << pivot_r << " " << pivot_c << ") = " << D[pivot_r][pivot_c] << endl;

			bool found = false;
			for (int j=pivot_r+1; j< m; j++)
			{
				if (D[pivot_r][j] != T(0))
				{
					/* process */
					found = true;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: found non-zero element in row at (" << pivot_r << " " << j << ")" << endl;


					/* If D[pivot_r][pivot_r] divides D[pivot_r][j], we subtract the appropriate 
					    multiple of column [pivot_r] from column  [j].
					*/
					T quotient = D[pivot_r][j]/D[pivot_r][pivot_c];
					T remainder = D[pivot_r][j]%D[pivot_r][pivot_c];
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: test division quotient:  " << quotient << endl;
	debug << "Smith_normal_form: test division remainder: " << remainder << endl;
}
					if (remainder == T(0))
					{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: pivot (" << pivot_r << " " << pivot_c << ") = " << D[pivot_r][pivot_c] << " divides row element (" << pivot_r << " " << j << ") = " << D[pivot_r][j]<< endl;
		
						/* subtract quotient times column pivot_c from column j in D and Q */
						for (int k=pivot_r; k<n; k++)
							D[k][j] -= quotient * D[k][pivot_c];
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: matrix after column operation:" << endl;
	print (D, debug, 3, "Smith_normal_form: ");
}
						/* collect the contribution to Q, we are multiplying Q on the right by the
						   elementary matrix that subtracts quotient times column i from column j
						*/
						for (int k=0; k<m; k++)
							Q[k][j] -= quotient * Q[k][pivot_c];
	
					
/*####################################################################*/
//test_matrices(D, P, Q, A,"Testing column reduction by multiple...");
/*####################################################################*/
	
					}
					else
					{
if (false && debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: DEV: D matrix before right multiplication: " << endl;
	print (D, debug, 3, "Smith_normal_form: ");
	debug << "Smith_normal_form: DEV: Q matrix before right multiplication: " << endl;
	print (Q, debug, 3, "Smith_normal_form: ");
	debug << "Smith_normal_form: DEV: pivot_r = " << pivot_r << " j = " << j << endl;
}
							
						/* We want to replace D[pivot_r][pivot_c] and D[pivot_r][j] with their gcd and zero  respectively */
						T u = D[pivot_r][pivot_c];
						T v = D[pivot_r][j] ;
						T c = T(0);
						T d = T(0);
						
						T t = ext_gcd(u,v,d,c);
						
						/* ext_gcd returns d,c such that ud+vc=1, Cohen's algorithm requires ud-vc=1 */
						c*=-1;
							
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: gcd t = " << t << ", d = " << d << "c = " << c << endl;
	debug << "Smith_normal_form: check u = " << u << " v = " << v << " u*d-v*c = " << u*d-v*c << endl;
}
	
						/* In terms of Cohen's description of Smith normal form,  we now evaluate the multiples a and b	*/
						T a = u/t;
						T b = v/t;
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: check a*d-b*c = " << a*d-b*c << endl;
}
						
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: a = " << a << ", b = " << b << endl;
						
						/* Multiply D on the right by the matrix
						
						   d_i,i      -b_i,j
						  -c_j,i       a_j,j
						  
						   and 1s in the other diagonal places.  This replaces (i,i) and (i,j) with t and 0 respectively, so we can 
						   just set these values before calculating the rest of the entries.  For each row below row i we replace
						   
						   (k,i) = d*(k,i)-c*(k,j)
						   (k,j) = -b*(k,i)+a*(k,j)						   					   
						*/
						
						D[pivot_r][pivot_c] = t;
						D[pivot_r][j] = T(0);

						for (int k=pivot_r+1; k<n; k++)
						{
							T temp = D[k][pivot_c];
							D[k][pivot_c] = d*D[k][pivot_c]-c*D[k][j];
							D[k][j] = a*D[k][j]-b*temp;							
						}
						
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: matrix after right multiplication: " << endl;
	print (D, debug, 3, "Smith_normal_form: ");
}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: DEV: Q matrix after right multiplication: " << endl;
	print (Q, debug, 3, "Smith_normal_form: ");
}
	
						/* collect the contribution to Q, we are multiplying Q on the right by the matrix in the comment above */
						for (int k=0; k<m; k++)
						{
							T temp = Q[k][pivot_c];
							Q[k][pivot_c] = d*Q[k][pivot_c]-c*Q[k][j];
							Q[k][j] = a*Q[k][j]-b*temp;
						}
	
/*####################################################################*/
//test_matrices(D, P, Q, A,"Testing column hcf...");
/*####################################################################*/
					}
				}
			}
						
			if (!found && again)
			{	
				/* We've already processed the column below (pivot_r pivot_c) at least once and now can't find a non-zero element in the row after (pivot_r pivot_c), so we're done. */
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: no more non-zero entries around (" << pivot_r << " " << pivot_c << ")" << endl;
				break;
			}
			else // (found && again) || !again
			{
				/* Either we've just found a non-zero entry in the row after (pivot_r pivot_c) having already processed the column at least once (i.e again==true), 				
				   or we're going around the loop for the first time (again==false).  We know, however, that when we're here the row 
				   after (pivot_r pivot_c) is all zero.
				*/
				again = false;
				/* look down the column */
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: processing column " << pivot_c << " around pivot (" << pivot_r << " " << pivot_c << ") = " << D[pivot_r][pivot_c] << endl;

				for (int j=pivot_r+1; j< n; j++)
				{						
					if (D[j][pivot_c] != T(0))
					{
						/* process */
						again = true;
						
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: found non-zero in column element at (" << j << " " << pivot_c << ")" << endl;


						/* If D[pivot_r][pivot_c] divides D[j][pivot_c], we subtract the  appropriate 
						   multiple of row [pivot_r] from row [j]
						*/
						T quotient = D[j][pivot_c]/D[pivot_r][pivot_c];
						T remainder = D[j][pivot_c]%D[pivot_r][pivot_c];
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: test division quotient:  " << quotient << endl;
	debug << "Smith_normal_form: test division remainder: " << remainder << endl;
}
						if (remainder == T(0))
						{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: pivot (" << pivot_r << " " << pivot_c << ") = " << D[pivot_r][pivot_c] << " divides column element (" << j << " " << pivot_c << ") = " << D[j][pivot_c]<< endl;

							/* subtract quotient times row pivot_r from row j */
							for (int k=0; k<m; k++)
								D[j][k] -= quotient * D[pivot_r][k];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: matrix after row operation:" << endl;
	print (D, debug, 3, "Smith_normal_form: ");
	debug.flush();
}

							/* collect the contribution to P, we are multiplying P on the left by the 
							   elementary matrix that subtracts *div_res.q times row i from row j
							*/
							for (int k=0; k<n; k++)
								P[j][k] -= quotient * P[pivot_r][k];
									
/*####################################################################*/
//test_matrices(D, P, Q, A,"Testing row reduction by multiple...");
/*####################################################################*/
						}	
						else
						{

							/* This time we want to replace D[pivot_r][pivot_c] and D[j][pivot_c] 
							   with their gcd and zero respectively
							*/
							T u = D[pivot_r][pivot_c];
							T v = D[j][pivot_c];
							T c = T(0);
							T d = T(0);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: check u = " << u << " v = " << v << endl;
							
							T t = ext_gcd(u,v,d,c);

							/* ext_gcd returns d,c such that ud+vc=1, Chhen's algorithm requires ud-vc=1 */
							c*=-1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: gcd t = " << t << ", d = " << d << "c = " << c << endl;
	debug << "Smith_normal_form: check u = " << u << " v = " << v << " u*d-v*c = " << u*d-v*c << endl;
}

							/* In terms of Cohen's description of Smith normal form,  we now evaluate the multiples a and b	*/
							T a = u/t;
							T b = v/t;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Smith_normal_form: a = " << a << ", b = " << b << endl;
					
							/* Multiply D on the left by the matrix
						
							   d_i,i      -c_i,j
							  -b_j,i       a_j,j
						  
							   and 1s in the other diagonal places.  This replaces (i,i) and (j,i) with t and 0 respectively, 
							   so we can just set these values before calculating the rest of the entries.  For each column 
							   after column i we replace
						   
							   (i,k) = d*(i,k)-c*(j,k)
							   (j,k) = -b*(i,k)+a*(j,k)
							*/
							

							for (int k=0; k<m; k++)
							{
								T temp = D[pivot_r][k];
								D[pivot_r][k] = d*D[pivot_r][k]-c*D[j][k];
								D[j][k] = a*D[j][k]-b*temp;
							}				
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: matrix after left multiplication: " << endl;
	print (D, debug, 3, "Smith_normal_form: ");
	debug.flush();
}
							/* collect the contribution to P, we are multiplying P on the left by the matrix in the comment above */
							for (int k=0; k<n; k++)
							{
								T temp = P[pivot_r][k];
								P[pivot_r][k] = d*P[pivot_r][k]-c*P[j][k];
								P[j][k] = a*P[j][k]-b*temp;
							}
/*####################################################################*/
//test_matrices(D, P, Q, A,"Testing row hcf...");
/*####################################################################*/
						}
					}
				}				
			}				
		} while (again);
		
		pivot_r++;
		pivot_c++;
		
	} while (!complete && pivot_r < n && pivot_c < m);

	if (!silent_operation)
		cout << endl;


	/* invert any rows of D that have negative values in D[i][i] and 
	   collect the inversion in P
	*/
	for (int r=0; r< rank; r++)
	{
		if (D[r][r] < T(0))
		{
			for (int c=0; c< m; c++)
				D[r][c] *= T(-1);

			for (int c=0; c< n; c++)
				P[r][c] *= T(-1);

		}
	}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Smith_normal_form: reached Smith normal form" << endl;
	print (D, debug, 3, "Smith_normal_form: ");
}
	

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
   debug << "Smith_normal_form: P ended up as " << endl;	
	print (P, debug, 3, "Smith_normal_form: ");
    debug << "Smith_normal_form: Q ended up as " << endl;	
	print (Q, debug, 3, "Smith_normal_form: ");
	debug.flush();
}

/*####################################################################*/
if (h_control.test_Smith_normal_form)
	test_matrices(D, P, Q, A,"Testing Smith normal form...");
/*####################################################################*/
	
	return rank;
}

template <class T, class St> 
matrix<T,St> SNF_inverse (const matrix<T,St>& M, homology_control h_control)
{
	int n = M.numcols();
	
	matrix<T,St> P(n,n);
	matrix<T,St> Q(n,n);
	matrix<T,St> D = M;
	
	Smith_normal_form (M,D,P,Q,h_control);													
	
	matrix<T,St> Minv= Q*P;

	matrix<T,St> I(n,n);
	
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<i; j++)
			I[i][j] = 0;

		I[i][i] = T(1);

		for (int j=i+1; j<n; j++)
			I[i][j] = 0;
	}

	matrix<T,St> check = M*Minv;

	if (check != I)
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "SNF_inverse: inverse check..." << endl;
	print(check, debug, 3, "inverse: ");
}
		cout << "\nincorrect SNF_inverse!" << endl;
		exit(0);
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "SNF_OK" << endl;
	
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "SNF_inverse: inverse check OK" << endl;

	return Minv;	
}


/* homology_generators calculates the k-th homology or cohomology generators from the boundary matrices Delta_k : C_k --> C_{k-1} and Delta_kp1 : : C_{k+1} --> C_k.  
   The function works with coeficients only and has no understanding of the basis elements for the chain groups C_n.  It returns the number of generators, the number of 
   torsion generators, a vector of length the number of torsion generators containing the torsion of the corresponding generator and a list generators as vectors of 
   coefficients (type T), the first of which are the torsion generators (if any exist), followed by the free generators.
   
   For situations where we have an exact sequence of 0 --> C_0 --> C_1 --> ...  --> C_k --> 0, to indicate we are calculating H_0 or H^0, the matrix Delta_kp1 should be empty 
   and to indicate that we are calculating H_k or H^k, the matrix Delta_k should be empty.  These two conditions are reflected in the booleans zero_Delta_kp1 and 
   zero_Delta_k.  
   
   If zero_Delta_kp1 is true, for homology there are no boundary generators and so H_0 is generated by the basis for the kernel Z_0.  For cohomology, it means that all the
   co-chains in C^0 are cocycles, since they vanish on zero (the image of the map 0 --> C_0) and so Z^0 is all of C^0.

   If zero_Delta_k is true, for homology the whole of C_0 lies in the kernel of the map C_k --> 0 and for cohomology it means there are no coboundary cochains in C^k, meaning
   that H^k is generated by the basis for Z^k.
*/
template <class T, class St> homology_generators_return<T> homology_generators(matrix<T,St>& Delta_k, matrix<T,St>& Delta_kp1, bool cohomology, homology_control h_control)	
{
	int num_km1_generators = Delta_k.numrows(); // k minus 1
	int num_k_generators = (Delta_k.numcols() != 0?Delta_k.numcols():Delta_kp1.numrows());
	int num_kp1_generators = Delta_kp1.numcols(); // k plus 1
	
	bool field_coefficients = h_control.field_coefficients;
	bool silent_operation=h_control.silent_operation;
	
	homology_generators_return<T> return_data;
	
	bool zero_Delta_kp1 = (Delta_kp1.numcols()==0);
	bool zero_Delta_k = (Delta_k.numcols()==0);

/*	
if (debug_control::DEBUG >= debug_control::BASIC)	
{
	debug << "homology_generators: Delta_k:" << endl;
	print(Delta_k, debug,3,"homology_generators: ");
	debug << "homology_generators: Delta_kp1:" << endl;
	print(Delta_kp1, debug,3,"homology_generators: ");
}
*/
	matrix<T,St> Delta_k_P(num_km1_generators,num_km1_generators);
	matrix<T,St> Delta_k_Q(num_k_generators,num_k_generators);
	matrix<T,St> Delta_k_D = Delta_k;
	int Delta_k_rank = 0;
	
	if(!zero_Delta_k)
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "homology_generators: calculating Smith normal form of Delta_k" << endl;

//		Delta_k_rank = Smith_normal_form (Delta_k,Delta_k_D,Delta_k_P,Delta_k_Q,field_coefficients,silent_operation);													
		Delta_k_rank = Smith_normal_form (Delta_k,Delta_k_D,Delta_k_P,Delta_k_Q,h_control);													

		if (!silent_operation)
			cout << "calculated Smith normal form of Delta_k" << endl;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "homology_generators: Smith normal form for Delta_k" << endl;
	print (Delta_k_D, debug, 3, "homology_generators: ");

	debug << "homology_generators: Delta_k_P" << endl;
	print (Delta_k_P, debug, 3, "homology_generators: ");

	debug << "homology_generators: Delta_k_Q" << endl;
	print (Delta_k_Q, debug, 3, "homology_generators: ");

}
	}
	
	matrix<T,St> Delta_kp1_P(num_k_generators,num_k_generators);
	matrix<T,St> Delta_kp1_Q(num_kp1_generators,num_kp1_generators);
	matrix<T,St> Delta_kp1_D = Delta_kp1;
	int Delta_kp1_rank = 0;
	
	if (!zero_Delta_kp1)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "homology_generators: calculating Smith normal form of Delta_kp1" << endl;

//		Delta_kp1_rank = Smith_normal_form (Delta_kp1,Delta_kp1_D,Delta_kp1_P,Delta_kp1_Q,field_coefficients,silent_operation);													
		Delta_kp1_rank = Smith_normal_form (Delta_kp1,Delta_kp1_D,Delta_kp1_P,Delta_kp1_Q,h_control);													

	if (!silent_operation)
		cout << "calculated Smith normal form of Delta_kp1" << endl;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "homology_generators: Smith normal form for Delta_kp1" << endl;
	print (Delta_kp1_D, debug, 3, "homology_generators: ");

	debug << "homology_generators: Delta_kp1_P" << endl;
	print (Delta_kp1_P, debug, 3, "homology_generators: ");

	debug << "homology_generators: Delta_kp1_Q" << endl;
	print (Delta_kp1_Q, debug, 3, "homology_generators: ");

}
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "homology_generators: Delta_k_rank = " << Delta_k_rank <<endl;
	debug << "homology_generators: Delta_kp1_rank = " << Delta_kp1_rank <<endl;
	
	debug << "homology_generators: Delta_k is " << Delta_k.numrows() << " by " << Delta_k.numcols() << endl;
	debug << "homology_generators: Delta_kp1 is " << Delta_kp1.numrows() << " by " << Delta_kp1.numcols() << endl;	
}

	int num_cycle_generators = (cohomology? num_k_generators-Delta_kp1_rank : num_k_generators-Delta_k_rank);
	int num_boundary_generators = (cohomology? Delta_k_rank : Delta_kp1_rank);

if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "homology_generators: num_cycle_generators = " << num_cycle_generators << ", num_boundary_generators = " <<  num_boundary_generators << endl;

	matrix<T,St> Gk(num_k_generators,num_cycle_generators+num_boundary_generators);		

	if (cohomology)	
	{
		/*  H^k = Z^k/B^k = ker(\delta^k)/im(\delta^{k-1}) = ker(Delta_kp1^T)/im(Delta_k^T)
				
			Delta_kp1^T is n_{k+1} x n_k and ker(Delta_kp1^T) is the last n_k-r_{k+1} columns of P_{k+1}^T, i.e the last rows of P_{k+1}.  
			Since P_{k+1} has n_k rows, these last rows start at index n_k-(n_k-r_{k+1}) = r_{k+1}.

			Delta_k^T is n_k x n_{k-1} and im(Delta_k^T) is the first r_k columns of (Q_k^T)^{-1}D_k^T, i.e the first rows of D_kQ_k^{-1}
			
			If Delta_k_rank is zero, then H^k is the free group generated by the basis for ker(Delta_kp1^T).
			
			Otherwise, we create the matrix Gk = [Z^k,B^k] and reduce it to echelon form to determine the coefficient matrix N, which is a 
			presentation matrix for H^k. The columns of N are relators in the generators z_i of the free group Z^k that generate the 
			kernel, B^k, of the unique hommomorphism from that free group to H^k with the given generator sets for Z^k and B^k			
		*/	
		
		matrix<T,St> Delta_k_Q_inv;
		matrix<T,St> Delta_k_T_image;
		
		if (!zero_Delta_k)
		{
			if (!silent_operation)
				cout << "Starting inverse calculation for Delta_k_Q (" << Delta_k_Q.numrows() << "x" << Delta_k_Q.numcols() << ")..." << flush;
			
			if (field_coefficients)
				Delta_k_Q_inv = inverse(Delta_k_Q, field_coefficients);
			else
				Delta_k_Q_inv = SNF_inverse(Delta_k_Q, h_control);

			if (!silent_operation)
				cout << "done" << endl;
			

if (debug_control::DEBUG >= debug_control::BASIC)	
{
	debug << "homology_generators: Delta_k_Q_inv:" << endl;
	print(Delta_k_Q_inv, debug,3,"homology_generators: ");
}	
			Delta_k_T_image = Delta_k_D*Delta_k_Q_inv;

if (debug_control::DEBUG >= debug_control::BASIC)	
{
	debug << "homology_generators: Delta_k transpose column space Delta_k_D*Delta_k_Q_inv generated by the first Delta_k_rank rows of:" << endl;
	print(Delta_k_T_image, debug,3,"homology_generators: ");
}	
		}
			
		for (int i=0; i<num_cycle_generators; i++)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "homology_generators: Gk column " << i << " is Delta_kp1_P row " << Delta_kp1_rank+i << endl;
		
			if (zero_Delta_kp1)
			{
				/* The kernel is the whole of C_0 so set the first num_k_generators columns of G_K to the identity matrix */
				for (int j=0; j< i; j++)
					Gk[j][i] = 0;

				Gk[i][i] = 1;

				for (int j=i+1; j< num_k_generators; j++)
					Gk[j][i] = 0;
				
			}
			else
			{
				for (int j=0; j< num_k_generators; j++)
					Gk[j][i] = Delta_kp1_P[Delta_kp1_rank+i][j];
			}
		}

		for (int i=num_cycle_generators; i<num_cycle_generators+num_boundary_generators; i++)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "homology_generators: Gk column " << i << " is Delta_k_T_image row " << i-num_cycle_generators << endl;
			
			for (int j=0; j< num_k_generators; j++)
				Gk[j][i] = Delta_k_T_image[i-num_cycle_generators][j];
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "homology_generators: the cohomology matrix Gk=[Z^k|B^k] is n_k = " << num_k_generators << " by (n_k-r_{k+1})=" <<  num_k_generators-Delta_kp1_rank << " + r_k=" << Delta_k_rank << " i.e. " << Delta_k_rank+num_k_generators-Delta_kp1_rank << " columns" << endl;
	debug << "homology_generators: Gk:" << endl;
	print(Gk, debug,3,"homology_generators: ");
}
	}
	else // homology
	{
		/*  H_k = Z_k/B_k = ker(Delta_k)/im(Delta_kp1)
				
			if Delta_k is non-empty,
			Delta_k is n_{k-1} x n_k and ker(Delta_k) is the last n_k - r_k columns of Q_k
			
			if Delta_k is empty we are dealing with a zero_Delta_k case, so Z_k is the whole of C_k

			Delta_kp1 is n_k x n_{k+1} and im(Delta_kp1) is the  first r_{k+1} columns of (P_{k+1})^-1 D_{k+1}
				
			If Delta_kp1_rank is zero, then H_k is the free group generated by the basis for ker(Delta_k).
			
			Otherwise, we create the matrix Gk = [Z_k,B_k] and reduce it to echelon form to determine the coefficient matrix N, which is a 
			presentation matrix for H_k. The columns of N are relators in the generators z_i of the free group Z_k that generate the 
			kernel, B_k, of the unique hommomorphism from that free group to H_k with the given generator sets for Z_k and B_k
		*/	

		matrix<T,St> Delta_kp1_P_inv;
		matrix<T,St> Delta_kp1_image;
		
		if (!zero_Delta_kp1)
		{		
			if (!silent_operation)
				cout << "Starting inverse calculation for Delta_kp1_P (" << Delta_kp1_P.numrows() << "x" << Delta_kp1_P.numcols() << ")..." << endl;

			if (field_coefficients)
				Delta_kp1_P_inv = inverse(Delta_kp1_P, field_coefficients);
			else
				Delta_kp1_P_inv = SNF_inverse(Delta_kp1_P, h_control);

			if (!silent_operation)
				cout << "done" << endl;

if (debug_control::DEBUG >= debug_control::BASIC)	
{
	debug << "homology_generators: Delta_kp1_P_inv:" << endl;
	print(Delta_kp1_P_inv, debug,3,"homology_generators: ");

}	
			Delta_kp1_image = Delta_kp1_P_inv * Delta_kp1_D;

if (debug_control::DEBUG >= debug_control::BASIC)	
{
	debug << "homology_generators: Delta_kp1 column space Delta_kp1_P_inv * Delta_kp1_D generated by the first " << Delta_kp1_rank << " columns of:" << endl;
	print(Delta_kp1_image, debug,3,"homology_generators: ");
}
		}
/*
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
{
	debug << "homology_generators: Delta_kp1 column space Delta_kp1_P_inv * Delta_kp1_D generated by columns:" << endl;
	for (int c=0; c< Delta_kp1_rank; c++)
	{
		vector<T> k_chain(num_k_generators);
		for (int r=0; r< num_k_generators; r++)
			k_chain[r] = Delta_kp1_image[r][c];
		
		debug << "homology_generators: ";
		for (int r=0; r< num_k_generators; r++)
			debug << setw(2) << k_chain[r];
		
		debug << ' ';
		print_k_chain(debug, k_chain, switch_data,k);

		debug << endl;
	}
}
*/
		for (int i=0; i<num_cycle_generators; i++)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "homology_generators: Gk column " << i << " is Delta_k_Q column " << Delta_k_rank+i << endl;
		
			if (zero_Delta_k)
			{
				/* The kernel is the whole of C_k so set the first num_k_generators columns of G_K to the identity matrix */
				for (int j=0; j< i; j++)
					Gk[j][i] = 0;

				Gk[i][i] = 1;

				for (int j=i+1; j< num_k_generators; j++)
					Gk[j][i] = 0;
				
			}
			else
			{
				for (int j=0; j< num_k_generators; j++)
					Gk[j][i] = Delta_k_Q[j][Delta_k_rank+i];
			}
		}

		for (int i=num_cycle_generators; i<num_cycle_generators+num_boundary_generators; i++)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "homology_generators: Gk column " << i << " is Delta_kp1_image column" << i-num_cycle_generators << endl;
			
			for (int j=0; j< num_k_generators; j++)
				Gk[j][i] = Delta_kp1_image[j][i-num_cycle_generators];
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "homology_generators: the homology matrix Gk=[Z_k|B_k] is n_k = " << num_k_generators << " by (n_k-r_k)=" <<  num_cycle_generators << " + r_{k+1}=" << Delta_kp1_rank << " i.e. " << num_cycle_generators+num_boundary_generators << " columns" << endl;	
	debug << "homology_generators: Gk:" << endl;
	print(Gk, debug,3,"homology_generators: ");
}

	}

/*
if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "homology_generators: first " << num_cycle_generators << " columns of Gk written as the cocycle generators of Z^k:" << endl;
	for (int i = 0; i< num_cycle_generators; i++)
	{ 
		vector<T> generator(num_k_generators);
		for (int j=0; j< num_k_generators; j++)
			generator[j] = Gk[j][i];
		
		debug << "homology_generators:   " << i+1 << ": ";
		print_characteristic_function(debug, generator, reverse_k_chain_map, cohomology, n, k);			
		debug << endl;
	}
}
*/		
	matrix<T,St> Gk_copy = Gk;
	
	/* need to watch for the case that num_boundary_generators = num_cycle_generators = 0 */
	if (Gk.numcols() !=0)
	{
		if (!silent_operation)
			cout << "Starting echelon calculation for Gk (" << Gk.numrows() << "x" << Gk.numcols() << ")..." << flush;
			
		echelon(Gk_copy,field_coefficients,true);

		if (!silent_operation)
			cout << "done" << endl;

	}
		
	/* invert any rows with a pivot equal to -1, since we want the columns of N to be the
	   positive linear combinations of the first (num_k_generators-Delta_kp1_rank) columns
	*/
	for (int i=0; i< num_cycle_generators; i++)
	{
		if (Gk_copy[i][i] == T(-1))
		{
			for (int j=i; j< num_cycle_generators+num_boundary_generators; j++)
				Gk_copy[i][j] *= T(-1);
		}
	}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "homology_generators: reduced echelon form of Gk:" << endl;
	print(Gk_copy, debug,3,"homology_generators: ");
}
	

	int num_generators = 0;
	int num_torsion_generators = 0;
	vector<T>	torsion;
	
	if (num_boundary_generators !=0)
	{		
		/* the matrix N is the first num_cycle_generators rows of the last 
		   num_boundary_generators columns of the echelon form
		*/
		matrix<T,St> N(num_cycle_generators,num_boundary_generators);
		
		for (int r=0; r< num_cycle_generators; r++)
		for (int c=0; c< num_boundary_generators; c++)
			N[r][c] = Gk_copy[r][c+num_cycle_generators];
		
if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "homology_generators: N (" << num_cycle_generators << " by " << num_boundary_generators << "): " << endl;
	print(N, debug,3,"homology_generators: ");
	debug << flush;
}

		/* Check we have calculated N correctly: the matrix B_k should be Z_k * N, there are num_boundary_generators 
		   columns of N and (obviously) the same number of boundary generators 
		*/
		for (int r=0; r< num_k_generators; r++)
		for (int c=0; c< num_boundary_generators; c++)
		{
			T B_element = 0;
			for (int i=0; i< num_cycle_generators; i++)
				B_element += Gk[r][i] * N[i][c];
			
			if (B_element != Gk[r][c+num_cycle_generators])
			{
				cout << "Boundary element (" << r << ',' << c << ") is not the linear combination of cycles specified by the corresponding column of N" << endl;
				exit(0);
			}
		}
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "homology_generators: B_k = Z_k * N check OK:" << endl;

		/* evaluate the Smith normal form of N */
		matrix<T,St> N_P(N.numrows(),N.numrows());
		matrix<T,St> N_Q(N.numcols(),N.numcols());
		matrix<T,St> N_D = N;

if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "homology_generators: calculating Smith normal form of N" << endl;
			
//		int N_rank = Smith_normal_form (N,N_D,N_P,N_Q,field_coefficients,silent_operation); 
		int N_rank = Smith_normal_form (N,N_D,N_P,N_Q,h_control); 

	if (!silent_operation)
		cout << "calculated Smith normal form of N" << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "homology_generators: Smith normal form of N, rank " << N_rank << endl;
	print(N_D, debug,3,"homology_generators: ");

	debug << "homology_generators: N_P:" << endl;
	print(N_P, debug,3,"homology_generators: ");
}

		if (!silent_operation)
			cout << "Starting inverse calculation for N_P (" << N_P.numrows() << "x" << N_P.numcols() << ")..." << flush;
			
		matrix<T,St> N_P_inv(N.numrows(),N.numrows());	
		if (field_coefficients)
			N_P_inv = inverse(N_P, field_coefficients);
		else
			N_P_inv = SNF_inverse(N_P, h_control);

		if (!silent_operation)
			cout << "done" << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "homology_generators: N_P_inv:" << endl;
	print(N_P_inv, debug,3,"homology_generators: ");
}	
		/* look for torsion generators */
		for (int i = 0; i< N_rank; i++)
		{
			if (abs(N_D[i][i]) != T(1))
			{
				num_generators++;
				
				/* N_D[i][i] times the i-th new generator element z'_i for the free group on {z_i} is a boundary
				   and in terms of the original generators z'_i is the linear combination of the z_i determined
				   by the i-th column of N_P_inv.  Thus N_D[i][i]*(n_{1,i}z_1 + ... + n_{p,i}z_p), where
				   p = num_k_generators-Delta_kp1_rank, is a boundary and so (n_{1,i}z_1 + ... + n_{p,i}z_p) is a 
				   torsion generator.  				   
				*/				
				vector<T> generator(num_k_generators);
				
				for (int r=0; r< num_k_generators; r++)
				{				
					/* multiply along the first num_cycle_generators columns in the r-th row of Gk and down the i-th column of N */
					for (int c=0; c< num_cycle_generators; c++)
						generator[r] += Gk[r][c] * N_P_inv[c][i]; 				
				}
	
				return_data.generators.push_back(generator);				
//				ostringstream oss;
//				print_characteristic_function(oss, generator, reverse_k_chain_map, cohomology, n, k);			
//				cocycles.push_back(oss.str());
			}
		}
	
		/* record the torsion element of H^3*/
		num_torsion_generators = num_generators;		
		torsion = vector<T>(num_torsion_generators);
		int index=0;
		for (int i = 0; i< N_rank; i++)
		{
			if (abs(N_D[i][i]) != T(1))
			{
				torsion[index] = abs(N_D[i][i]);
				index++;
			}
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "homology_generators: num_torsion_generators = " << num_torsion_generators << endl;
	debug << "homology_generators: torsion: ";
	for (int i=0; i< num_torsion_generators; i++)
		debug << torsion[i] << ' ';
	debug << endl;
}				
		/* Look for infinite cyclic generators, the last (num_k_generators-Delta_kp1_rank)-N_rank 
		   columns of N_P_inv determine infinite cyclic generators
		*/

if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "homology_generators: " << num_cycle_generators - N_rank << " infinite cyclic generators: " << endl;
		
		for (int i = N_rank;  i< num_cycle_generators; i++)
		{
			num_generators++;

			vector<T> generator(num_k_generators);

			for (int r=0; r< num_k_generators; r++)
			{				
				for (int c=0; c< num_cycle_generators; c++)
					generator[r] += Gk[r][c] * N_P_inv[c][i];					
			}

if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "homology_generators:   ";
	for (int r=0; r< num_k_generators; r++)
		debug << generator[r] << ' ';
	debug << endl;
}				

			return_data.generators.push_back(generator);				

//			ostringstream oss;
//			print_characteristic_function(oss, generator, reverse_k_chain_map, cohomology, n, k);			
//			cocycles.push_back(oss.str());			
		}		
	} // end of clause if (Delta_k_rank !=0) 		
	else
	{
		for (int c = 0;  c< num_cycle_generators; c++)
		{
			num_generators++;

			vector<T> generator(num_k_generators);
			for (int r=0; r< num_k_generators; r++)
				generator[r] = Gk[r][c];

			return_data.generators.push_back(generator);				

//			ostringstream oss;
//			print_characteristic_function(oss, generator, reverse_k_chain_map, cohomology, n, k);			
//			cocycles.push_back(oss.str());			
		}
	}

	return_data.num_generators = num_generators;
	return_data.num_torsion_generators = num_torsion_generators;
	return_data.num_torsion_generators = num_torsion_generators;
	return_data.torsion = torsion;

if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "homology_generators:   found " << num_generators << " generators, " << num_torsion_generators << " torsion generators" << endl;;
	
	return return_data;
}
