/************************************************************************************************************

                                Switch invariant calculations and supporting functions

NOT USED IN CURRENT VERSION
                                

Hmatrix H22inverse(Hmatrix& M)

template <class T, class St> 
void poly_invariant(const matrix<T,St>& switch_matrix, const matrix<T,St>& switch_matrix_inverse, 
                    string input_string, string title)

template <class T, class St> 
T study_determinant (const matrix<T,St>& H_matrix, string title, int n=0, int* rperm=0, int* cperm=0)

matrix<polynomial<scalar,int>,scalar> C_study_matrix(const matrix<polynomial<scalar,int>,scalar>& H_matrix, int n=0, int m=0, int* rperm=0, int* cperm=0) 

template <class T, class U, class St> 
matrix<polynomial<T,U>,St> C_study_matrix(const matrix<polynomial<T,U>,St>& H_matrix, int n=0, int m=0, int* rperm=0, int* cperm=0)

template <class T, class St> T delta_1(const matrix<T,St>& M, int codimension, string title, T delta_0 = T("0"))

template <class T, class St> bool minor_determinant(T& hcf, const matrix<T,St>& M, int* rperm, int* cperm,
                       int codimension, int r1, int r2, int c1, int c2, int c3, string title)

Qpmatrix adj2(Qpmatrix& M)
void simplify_matrix (Qpmatrix& M)

template <class T1, class U, class T2, class St>  void normalize(polynomial<T1,U>& poly, const matrix<T2,St>& switch_matrix) {}

void normalize(Hpolynomial& poly, const Hpmatrix& switch_matrix)

template <class T> void equality_test (T x, T y, bool A)

template <class T, class St> 
bool braid_rep(matrix<T,St>*& Matrix_rep, const matrix<T,St>& switch_matrix, const matrix<T,St>& switch_matrix_inverse, string input_string,
                int num_terms, int num_strings)

template <class T, class St> 
void R_module_rep(matrix<T,St>*& Matrix_rep, const matrix<T,St>& switch_matrix, const matrix<T,St>& switch_matrix_inverse, 
               string input_string)

bool rat_minor_determinant (Qpmatrix* Matrix_rep, int N, int Nr1, int Nc1, int Nc2, int* rperm, int* cperm, string title, polynomial<scalar,bigint>& hcf)

void rat_poly_invariant(const Qpmatrix& switch_matrix, const Qpmatrix& switch_matrix_inverse,
                    string input_string, string title)

void commutative_automorphism_invariant(const Qpmatrix& phi, const Qpmatrix& psi, string input_string, string title)

************************************************************************************************************/
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iomanip>

using namespace std;

extern ifstream     input;
extern ofstream 	output;
extern ofstream     debug;

extern bool SIDEWAYS_SEARCH;

#include <util.h>
#include <quaternion-scalar.h>
#include <polynomial.h>
#include <matrix.h>
#include <ctype.h>
#include <braid.h>

string parse_long_knot_input_string (string input_string);
void read_immersion_code (generic_code_data& code_data, string input_string);
void read_peer_code (generic_code_data& code_data, string input_string);

bool valid_braid_input (string input_string, int& num_terms, int& num_strings);


Hmatrix H22inverse(Hmatrix& M)
{
	Quaternion& A = M[0][0];
	Quaternion& B = M[0][1];
	Quaternion& C = M[1][0];
	Quaternion& D = M[1][1];
	
	Hmatrix invM(2,2);
	
    /* check existence of inverse: i.e. C^-1D != A^-1B */
    if(conj(C)*D*Quaternion(norm_sqrd(A)) == conj(A)*B*Quaternion(norm_sqrd(C)))
    {
		throw(matrix_error("attempt to calculate inverse of a non-invertible 2x2 quaternionic matrix"));
	}
	else
	{	
		Quaternion C_inv = conj(C) * (Quaternion(1)/Quaternion(norm_sqrd(C)));
		Quaternion A_inv = conj(A) * (Quaternion(1)/Quaternion(norm_sqrd(A)));
		Quaternion delta = C_inv * D - A_inv * B;
		Quaternion delta_inv = conj(delta) * (Quaternion(1)/Quaternion(norm_sqrd(delta)));

		invM[0][0] = C_inv * D * delta_inv * A_inv;
		invM[0][1] = Quaternion(-1) * A_inv * B * delta_inv * C_inv;
		invM[1][0] = Quaternion(-1) * delta_inv * A_inv;
		invM[1][1] = delta_inv * C_inv;
		
	    return invM;
    }
}


/* poly_invariant calculates the Delta_0 and Delta_1 determined by the supplied switch matrix and its inverse for the 
   braid or immersion code given in the input string.  It is used only where the elements of the switch matrix are 
   of type polynomial<T,U> and not rational<polynomial<T,U> > ,i.e. Burau or quaternionic switches currently.  Note that 
   this means that the switch matrix N-factor is one in this function, which consequently is optimized for this case.
*/
template <class T, class St> 
void poly_invariant(const matrix<T,St>& switch_matrix, const matrix<T,St>& switch_matrix_inverse, 
                    string input_string, string title)
{
	int num_terms;
	int num_strings;
	
	bool calculate_polys = true;
	bool virtual_crossings_present = false;
	
	matrix<T,St>*		Matrix_rep;
	
	if (title.length())
   	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\n\n" << title << endl;
		if (!braid_control::RAW_OUTPUT)
			output << "\n\n" << title;
   	}		

	if (!braid_control::RAW_OUTPUT)
	{
		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
		output << input_string;		
	}
	
	/* note the N-factor of the switch matrix, this is done explicitly as a debug check */
	int switch_matrix_N_factor = switch_matrix.numcols()/2;
	
if(braid_control::DEBUG >= braid_control::DETAIL)
	debug << "poly_invariant: switch_matrix_N_factor = " << switch_matrix_N_factor << endl;

	/* The two cases of the next conditional clause are in-line rather than function calls only because
	   the code appear only once
	*/
//	if (input_string.find('(') != string::npos)
	if (input_string.find('(') != string::npos || input_string.find('[') != string::npos)
	{
		/* record whether this is a virtual knot or a classical one */
		if (input_string.find('*') != string::npos)
			virtual_crossings_present = true;

		/* check for a long knot */
		
		if (input_string.find('L') != string::npos)
		{
			braid_control::LONG_KNOT = true;
			
			input_string = parse_long_knot_input_string(input_string);
			
			
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "poly_invariant: long knot: " << input_string << endl;
	
			if (input_string == "link")
			{				
				if (!braid_control::SILENT_OPERATION)
					cout << "\nError! long knot indicator provided for the peer code of a link" << endl;
			
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Error! long knot indicator provided for the peer code of a link" << endl;
				}
					
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "poly_invariant: Error! long knot indicator provided for the peer code of a link" << endl;

				return;
			}
		}
		else
			braid_control::LONG_KNOT = false;
	
		R_module_rep(Matrix_rep, switch_matrix, switch_matrix_inverse, input_string);
	}
	else if (valid_braid_input(input_string, num_terms, num_strings))
	{ 

		/* remove any unwanted qualifiers from the input string */
		string::size_type pos = input_string.find('{');
		if (pos != string::npos)
			input_string = input_string.substr(0,pos);


		/* record whether this is a virtual knot or a classical one */
		if (input_string.find('t') != string::npos || input_string.find('T') != string::npos)
			virtual_crossings_present = true;
			

		/* When using a braid representation, we do not calculate polynomials if the representation is the identity
		   matrix.  This check is done within braid_rep */	   
		calculate_polys = braid_rep(Matrix_rep, switch_matrix, switch_matrix_inverse,input_string, num_terms, num_strings);
	} 

	/* Under normal circumstances Matrix_rep is created by R_module_rep or braid_rep.  However, if R_module_rep
	   determines that an immersion code contains only virtual crossings it returns a zero Matrix_rep pointer. */
	if (Matrix_rep == 0)
		return;
		
	/* Note the size of the matrix representation, this is either square or in the case of long knots
	   has one more column than row, so we note the number of rows.  For elements of VB(n) this is the 
	   number of strings, but for R_module representations it is twice the number of classical crossings.
	*/
	int matrix_rows = Matrix_rep->numrows();
	int matrix_cols = Matrix_rep->numcols();

    if (calculate_polys)
    {	
		/* First determine the codimension 0 invariant.  In the case of a long knot we have to remove successive columns  
		   and calculate the gcd of the (Study) determinants of remianing submatrices.  In the case of a closed knot
		   the codimension 0 invariant is just the (Study) determinant of the matrix representation.
		   
		   We then determine whether to calculate delta_1 based on the codimension zero value and the global booleans
		*/
		bool calc_delta_1;

		if (braid_control::WAIT_SWITCH)
		{
		    matrix_control::WAIT_INFO = true;
		    matrix_control::wait_threshold = braid_control::wait_threshold;
		    matrix_control::wait_count = 0;

		    /* the wait_threshold comes from the command line */
			if (!braid_control::SILENT_OPERATION)
				cout << "\nCalculating polynomials, please wait\n";
		}
		else
		    matrix_control::WAIT_INFO = false;

		/* We use delta_0 for the codimension 0 polynomial invariant in both the closed
		   knot and long knot cases.  This allows us to use common code for delta_1, the
		   codimension 1 polynomial invariant.
		*/ 		  
		T delta_0 = T("0");

		if (braid_control::LONG_KNOT)
		{	
			/* The codimension 0 polynomial invariant is formed by taking out each column of the 
			   matrix representation in turn and taking the determinant  of the rest, keeping a 
			   record of the gcd of the determinants calculated.
			*/
			int rperm[matrix_rows];
			int cperm[matrix_cols];
			bool unit_detected = false;

			for (int i = 0 ; i< matrix_rows; i++)
		    	rperm[i] = i;

			T op_0 = T("0"); // \def\op{\buildrel o \over p}
			T np_0 = T("0"); // \def\np{\buildrel n \over p}

			for (int i = 0 ; i< matrix_cols; i++)
			{
//				if (WAIT_SWITCH)
//				    matrix_control::wait_count = 0;

				/* take out column i */
				for (int j=0; j<i; j++)
    				cperm[j] = j;
				for (int j = i+1; j<= matrix_cols; j++)
					cperm[j-1] = j;

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "poly_invariant: removing column " << i << endl;

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "poly_invariant: cperm: ";
   	for (int j = 0 ; j< matrix_cols-1; j++)
		debug << cperm[j] << " ";
	debug << endl;
}
				T det = T("0");

				if (braid_control::QUATERNION)
					det =  study_determinant(*Matrix_rep, title, matrix_rows, rperm, cperm);
				else
					det =  determinant(*Matrix_rep, title, matrix_rows, rperm, cperm);

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "poly_invariant: p_0 generator det " << i << " = " << det << endl;

    			Laurent_scale (det);

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "poly_invariant: Laurent scaled p_0 generator det " << i << " = " << det << endl;


				if (i == 0)
				{
					op_0 = det;

					if (braid_control::WAIT_SWITCH)
					    matrix_control::wait_count = 0;

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "poly_invariant: op_0  set to det " << i << " = " << op_0 << endl;

				}
				else if (i == matrix_cols-1)
				{
					np_0 = det;

					if (braid_control::WAIT_SWITCH)
					    matrix_control::wait_count = 0;

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "poly_invariant: np_0  set to det " << i << " = " << np_0 << endl;

					if (unit_detected)
						break;  //we've jumped to here, so we're now done.

				}

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "poly_invariant:  p_0 gcd stands at " << delta_0 << endl;

				if (braid_control::DELTA_1_UNIT_CHECK && det.isunit())
				{
					unit_detected = true; 
					if (i != matrix_cols-1)
						i = matrix_cols-2; // jump to last column to determine np_0, i will be incremented by loop

					delta_0 = T("1");

					if (braid_control::EXTRA_OUTPUT)
					{
						if (!braid_control::SILENT_OPERATION)
							cout << "\nunit generator detected, terminating calculation" << endl;
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "unit generator detected, terminating calculation" << endl;
						}
					}

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "poly_invariant: unit generator detected, terminating calculation" << endl;
				}
				else
				{
			    	delta_0 = gcd(delta_0,det);					

					if (braid_control::DELTA_1_UNIT_CHECK && delta_0.isunit())
					{
						unit_detected = true; 
						if (i != matrix_cols-1)
							i = matrix_cols-2; // jump to last column to determine np_0, i will be incremented by loop

						if (braid_control::EXTRA_OUTPUT)
						{
							if (!braid_control::SILENT_OPERATION)
								cout << "\nunit gcd detected, terminating calculation" << endl;
							if (!braid_control::RAW_OUTPUT)
							{
								output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
								output << "unit gcd detected, terminating calculation" << endl;
							}
						}

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "poly_invariant: unit gcd detected, terminating calculation" << endl;
					}
				}
			}

			if (!braid_control::SILENT_OPERATION)
			{
	   			cout << "\np^{(0)} = " << delta_0 << endl;
	   			cout << "\n\\op^{(0)} = " << op_0 << endl;
	   			cout << "\n\\np^{(0)} = " << np_0 << endl;
			}
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "p^{(0)} = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << delta_0 << endl;
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "\\op^{(0)} = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << op_0 << endl;
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "\\np^{(0)} = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << np_0 << endl;

if (braid_control::DEBUG >= braid_control::SUMMARY)
{
   	debug << "poly_invariant: p_0  = p^{(0)} = " << delta_0 << endl;
   	debug << "poly_invariant: \\op^{(0)} = " << op_0 << endl;
   	debug << "poly_invariant: \\np^{(0)} = " << np_0 << endl;
}

			if (delta_0 == T("0"))		
				calc_delta_1 = true;
    		else
				calc_delta_1 = false;

		}
		else // closed knot
		{
			/* Delta_0 = 0 in the classical case, so only calculate it for virtuals, the codimension 0 
			   polynomial invariant is just the determinant of the presentation matrix, since it is a square matrix.
			*/
			if (virtual_crossings_present || braid_control::VERIFY_DELTA_0)
			{
		    	if (!braid_control::QUATERNION)
	    		{
					/* use deteminant function to calculate det(M-I) */
					delta_0 = determinant(*Matrix_rep, title);

					if (braid_control::ALEXANDER)
    				{
						/* set s=1 in delta if it involves s */
						if (delta_0.numvars() && strchr(delta_0.getvars(),'s'))
							set_to_one(delta_0,'s');
    				}

	    			Laurent_scale (delta_0);

					if (!braid_control::DISPLAY_DELTA_1_ONLY)
					{
						if (!braid_control::SILENT_OPERATION)
							cout << "\nDelta_0 = " << delta_0 << endl;
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "Delta_0 = ";
							if (braid_control::OUTPUT_AS_INPUT)
								output << "\n";
						}
						output << delta_0 << endl;
					}

if (braid_control::DEBUG >= braid_control::SUMMARY)
{
   	debug << "\npoly_invariant: Delta_0 = " << delta_0 << "\n";
}
	    		}
	    		else
	    		{
					/* we have to use the Study determinant */
					matrix<T,St>& matrix_rep = *Matrix_rep;
					delta_0 = study_determinant(matrix_rep, title);

					if (braid_control::NORMALIZING)
						normalize(delta_0, switch_matrix);

					if (!braid_control::DISPLAY_DELTA_1_ONLY)
					{
						if (!braid_control::SILENT_OPERATION)
							cout << "\nDelta_0 = " << delta_0 << endl;
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "Delta_0 = ";
							if (braid_control::OUTPUT_AS_INPUT)
								output << "\n";
						}
						output << delta_0 << endl;

					}

if (braid_control::DEBUG >= braid_control::SUMMARY)
{
   	debug << "\npoly_invariant: Study determinant = " << delta_0 << "\n";
}

	    		}

				if (!braid_control::SILENT_OPERATION && matrix_control::WAIT_INFO)
					cout << "\n";

				if (delta_0 == T("0"))		
					calc_delta_1 = true;
	    		else
					calc_delta_1 = false;

			}
			else
			{
	    		calc_delta_1 = true;

				if (!braid_control::DISPLAY_DELTA_1_ONLY)
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "\nInput contains no virtual crossings, so Delta_0 = 0" << endl;
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Input contains no virtual crossings, so Delta_0 = ";
						if (braid_control::OUTPUT_AS_INPUT)
							output << "\n";
					}
					output << "0" << endl;

				}
			}
		}
		
		if (!braid_control::CALCULATE_DELTA_0_ONLY && (braid_control::ALWAYS_CALCULATE_DELTA_1 || calc_delta_1))
		{
	    	if (braid_control::EXTRA_OUTPUT)
		    {
				/* if we're doing quaternions we don't need wait state
				   information if we're going to see the adjunct elements
				*/
				if ((braid_control::QUATERNION && matrix_cols > 4) ||
		    		(!braid_control::QUATERNION && matrix_cols > 5))
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "\n";
		    		matrix_control::WAIT_INFO = false;
				}
	    	}

//			if (braid_control::COMPLEX_STUDY_DELTA_1 || braid_control::REAL_STUDY_DELTA_1)
			if (braid_control::COMPLEX_STUDY_DELTA_1)
			{
				
				/* Here we map the (quaternionic) matrix representation into M_{2n,2m}(C[t,t^{-1}]) and  
				   determine the gcd of the codimension 1 determinants of the image.  If this gcd is 0 we 
				   determine the gcd of the codimension 2 determinants of the image.  This is an alternative 
				   approach to the normal procedure of taking the codimension 1 submatrices of M_n,m(H[t,t^{-1}), 
				   mapping each one into M_{2n,2m}(C[t,t^{-1}]), taking determinants and evaluating their gcd.
				
				   First, set up the Complex Study matrix and evaluate delta_1 from it. This call to delta_1 will  
				   remove rows and columns of the matrix in M_{2n,2m}(C).  
				*/
				int rperm[matrix_rows];
				int cperm[matrix_cols];

				for (int i = 0 ; i< matrix_rows; i++)
		    		rperm[i] = i;

				for (int i = 0 ; i< matrix_cols; i++)
		    		cperm[i] = i;

				matrix<T,St> C_study_matrix_rep = C_study_matrix(*Matrix_rep, matrix_rows, matrix_cols, rperm, cperm);
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "\npoly_invariant: complex Study matrix:\n";
	debug << C_study_matrix_rep;
	debug << "\n";
}	

				bool display_delta;
				if (braid_control::COMPLEX_STUDY_DELTA_1)
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "Evaluating the complex Study matrix" << endl;
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "\n\npoly_invariant: evaluating the complex Study matrix\n";	

					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Evaluating the complex Study matrix";
					}

					if (braid_control::EXTRA_OUTPUT)
					{
						if (!braid_control::SILENT_OPERATION)
							cout << "\nGenerators of Delta_1 calculated from the complex Study matrix" << endl;
						output << "\nGenerators of Delta_1 calculated from the complex Study matrix" << endl;
					}

					display_delta = false;
					/* It is mandatory to include delta_0 as a generator of delta_1 in this case, since
					   we are dealing with quaternions. */
					T delta = delta_1(C_study_matrix_rep, CODIMENSION_1, title, delta_0);

					int codimension;
					if (delta != T("0") && !braid_control::ALWAYS_CALCULATE_CODIMENSION_2_DELTA_1)
					{
						display_delta = true;
						codimension = CODIMENSION_1;
					}
					else
					{
						if (!braid_control::SILENT_OPERATION && matrix_control::WAIT_INFO)
						{
							cout << "\n";
							cout << "\nDelta_1^C=1" << endl;
						}
						
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "\n\npoly_invariant: Delta_1^C=1";
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "Delta_1^C=1";
						}

						/* take out two rows and columns from C_study_matrixptr at a time */
						delta = delta_1(C_study_matrix_rep, CODIMENSION_2, title, delta_0);

						if (delta == T("1"))
						{
							if (!braid_control::SILENT_OPERATION)
								cout << "\nDelta_2^C=1";

//							cout << "\nTrying real Study matrix" << endl;
if (braid_control::DEBUG >= braid_control::SUMMARY)
{	
//	debug << "\n\npoly_invariant: Delta_2^C=1, trying the real Study matrix" << endl;
	debug << "\n\npoly_invariant: Delta_2^C=1" << endl;
}	

							if (!braid_control::RAW_OUTPUT)
							{
								output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
								output << "Delta_2^C=1" << endl;
//								output << "\nTrying real Study matrix" << endl;
							}

//							REAL_STUDY_DELTA_1 = true;
							
						}
						else
						{
							display_delta = true;
							codimension = CODIMENSION_2;
						}
					}

					if (display_delta)
					{
						if (!braid_control::SILENT_OPERATION)
						{
							cout << "\nDelta_";
							cout << (codimension == 1? "1" : "2");
							cout << "^C = ";
						}
						
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "\nDelta_";
							output << (codimension == 1? "1" : "2");
							output << "^C = ";
						}

						if (braid_control::OUTPUT_AS_INPUT)
							output << "\n";
						if (!braid_control::SILENT_OPERATION)
							cout << delta << "\n";
		   				output << delta << "\n";

if (braid_control::DEBUG >= braid_control::SUMMARY)
{
	debug << "\npoly_invariant: Delta_";
	debug << (codimension == 1? "1" : "2");
	debug << "^C = " << delta << "\n";	
}
					}
				}


/******************************************************************
				if (REAL_STUDY_DELTA_1)
				{
					display_delta = false;
					cout << "\nEvaluating the real Study matrix";
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "\n\npoly_invariant: evaluating the real Study matrix\n";	

					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Evaluating the real Study matrix";
					}

					  calculate the Real Study matrix in M_{4n}(R) and
					  use R_study_matrixptr to point to it.  We do this via
					  a function call just in order to keep the code here simple
					
					matrix<T,St> R_study_matrix_rep = R_study_matrix(C_study_matrix_rep);

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "\npoly_invariant: real Study matrix:\n";
	debug << R_study_matrix_rep;
	debug << endl;
}					

					if (braid_control::EXTRA_OUTPUT)
					{
						cout << "\n\nGenerators of Delta_1 calculated from the real Study matrix";
						output << "\n\nGenerators of Delta_1 calculated from the real Study matrix";
					}

					// set wait info true, the Real Study matrix is big! 
					matrix_control::WAIT_INFO = true;
	    			matrix_control::wait_count = 0;

				    * if the wait switch was given, the wait_threshold
	    			   comes from the command line
	    			if (!WAIT_SWITCH)
						wait_threshold = (QUATERNION? 2* DEFAULT_WAIT_THRESHOLD: DEFAULT_WAIT_THRESHOLD);
					*

				    cout << "\nCalculating polynomials, please wait\n";

					T delta = delta_1(R_study_matrix_rep, CODIMENSION_1, title);

					int number;
					if (delta.nonzero())
					{
						number = CODIMENSION_1;
						display_delta = true;
					}
					else 
					{
						cout << "\nDelta_1^R=0";
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "\n\npoly_invariant: Delta_1^R=0";

						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "Delta_1^R=0";
						}

						* take out two rows and columns from R_study_matrixptr at a time *
						delta = delta_1(R_study_matrix_rep, CODIMENSION_2, title);

						if (delta == T("1"))
						{
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "\n\npoly_invariant: Delta_2^R=1\n";	

							cout << "\nDelta_2^R=1";
							if (!braid_control::RAW_OUTPUT)
							{
								output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
								output << "Delta_2^R=1";
							}

							cout << "\nOut of ideas, jumping off nearest cliff...";
						}
						else
						{
							number = CODIMENSION_2;
							display_delta = true;
						}
					}

					if (display_delta)
					{
						cout << "\nDelta_";
						cout << (number == 1? "1" : "2");
						cout << "^R = ";
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "\nDelta_";
							output << (number == 1? "1" : "2");
							output << "^R = ";
						}

						if (braid_control::OUTPUT_AS_INPUT)
							output << "\n";
		   				cout << delta << "\n";
		   				output << delta << "\n";

if (braid_control::DEBUG >= braid_control::SUMMARY)
{
	debug << "\npoly_invariant: Delta_";
	debug << (number == 1? "1" : "2");
	debug << "^R = " << delta << "\n";	
}
					}
				}
*****************************************************************/
			}
			else
			{
if (braid_control::DEBUG >= braid_control::SUMMARY && braid_control::QUATERNION)
	debug << "\npoly_invariant: evaluating codimension 1 determinant in M_{n,m}(H) using Study determinant\n";

				/* Use the normal procedure of taking the codimension 1 submatrices of the matrix 
				   representation, taking (Study) determinants and evaluating their gcd.
				   
					We include delta_0 as a generator for E_1 in all cases, though strictly it is 
					only required in the quaternionic (non-commutative) case but it keeps the code 
					simple and the overhead of this appraoach is small.
				*/	
				T delta = delta_1(*Matrix_rep, CODIMENSION_1,title,delta_0);
				
				if (braid_control::EXTRA_OUTPUT)
		    	{
					if (!braid_control::SILENT_OPERATION)
						cout << "\n";
					output << "\n";
		    	}
		    	else if (!braid_control::SILENT_OPERATION && matrix_control::WAIT_INFO)
		    	{
					/* we need a carridge return because of the wait info */
					cout << "\n";
		    	}

		    	if (braid_control::QUATERNION)
		    	{
					if (braid_control::NORMALIZING)
						normalize(delta,switch_matrix);

					if (!braid_control::SILENT_OPERATION)
						cout << (braid_control::LONG_KNOT? "\np^{(1)} = ":"\nDelta_1^H = ");
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << (braid_control::LONG_KNOT? "p^{(1)} = ":"Delta_1^H = ");
					}
				}
				else
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "\nAlexander polynomial = ";
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Alexander polynomial = ";
					}
				}

				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
				if (!braid_control::SILENT_OPERATION)
					cout << delta << "\n";
		    	output << delta << "\n";

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "\n\npoly_invariant: " << (braid_control::LONG_KNOT? "p^{(1)} = ":"Delta_1 = ") << delta << "\n";	

			}

		}
		else
		{
		    if (braid_control::DISPLAY_DELTA_1_ONLY)
		    {
				/* We need a place marker in the output to identify why we're not displaying anything */
				if (!braid_control::SILENT_OPERATION)
					cout << "\nDelta_0 non zero, Delta_1 not calculated." << endl;
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Delta_0 non zero, Delta_1 not calculated." << endl;
				}
		    }
		}
    }
	
	delete Matrix_rep;
}

/* only a template because the compiler doesn't know this function is not going to be called by 
   matrices with non quaternionic coefficients
*/
template <class T, class St> 
T study_determinant (const matrix<T,St>& H_matrix, string title, int n=0, int* rperm=0, int* cperm=0)
{
	/* call C_study_matrix to create the image in M_{2n}(C) under \psi of the
	   n rows and columns of matrixptr indicated by rperm and cperm.  The image
	   will be stored under C_study_matrixptr.
	*/
	matrix<T,St> C_matrix = C_study_matrix(H_matrix, n, n, rperm, cperm);

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
    debug << "\nstudy_determinant: Complex Study matrix:\n";
    debug << C_matrix;
    debug << "\n";
}

	return determinant(C_matrix,title);

}

/* C_study_matrix calculates the image under a mapping from M_{n,m}(H) to M_{2n,2m}(C)
   of the n rows and m columns of H_matptr indicated by rperm and cperm and stores
   it in a new Hpmatrix created at C_matptr.
   
   The mapping is either the mu-mapping, the default, or the rho-mapping, selected via the
   boolean variable STUDY_RHO_MAPPING.  The mu-mapping is based on the fact that any quaternion 
   w + xi +yj +zk may be written as Z + Wi, where Z = w + xi and W = y + zi are complex numbers.
   The rho-mapping is based on the fact that w + xi +yj +zk may also be written as Z + Wj, where 
   Z = w + zk and W = y + xk are complex numbers.
   
*/
matrix<polynomial<scalar,int>,scalar> C_study_matrix(const matrix<polynomial<scalar,int>,scalar>& H_matrix, int n=0, int m=0, int* rperm=0, int* cperm=0) 
{
	return H_matrix;
}

template <class T, class U, class St> 
matrix<polynomial<T,U>,St> C_study_matrix(const matrix<polynomial<T,U>,St>& H_matrix, int n=0, int m=0, int* rperm=0, int* cperm=0)
{

if (braid_control::DEBUG >= braid_control::BASIC)
{
	if (braid_control::STUDY_RHO_MAPPING)
	    debug << "\nC_study_matrix: using rho-mapping to create complex Study matrix:" << endl;
	else
	    debug << "\nC_study_matrix: using mu-mapping to create complex Study matrix:" << endl;
}
	bool clean_up_params = false;
	
	typedef int It;
		
	if (n == 0)
	{
		n = H_matrix.numrows();
		m = H_matrix.numcols();
		rperm = new int[n];
		for (int i=0; i< n; i++)
			rperm[i] = i;
		cperm = new int[m];
		for (int i=0; i< m; i++)
			cperm[i] = i;
		clean_up_params = true;
	}
	
    matrix<polynomial<T,U>,St> C_matrix(2*n,2*m);

    for (int i=0; i<n; i++)
    for (int j=0; j<m; j++)
    {
		/* take the quaternionic polynomial at H_matrix[rperm[i]][cperm[j]]
		   and create the 2x2 matrix      Z       W
								      -\bar{W} \bar{Z}
		   at C_matptr[2*i][2*j], where Z and W are given by either the mu-mapping
		   or the rho-mapping, dependent on the value of STUDY_RHO_MAPPING.

		   The \bar function is complex conjugation.
		   
		   The mu-mapping takes (r,i,j,k) to   (r,i)   (j,k)
		                                       (-j,k)   (r,-i)
									
		   The rho-mapping takes (r,i,j,k) to   (r,k)   (j,i)
		                                       (-j,i)   (r,-k)
											  
		   That is, i and k are interchanged.
		*/
		C_matrix[2*i][2*j] = H_matrix[rperm[i]][cperm[j]];
		
		/* Set the 00 position in relation to the  base, 2*i,2*j.  Moving along the pterms,
		   for the mu-mapping set the j and k  components of the coefficients to zero
		   for the rho-mapping first copy the k component to the i component, then set j and k to zero.
		*/
		pterm<T,U>* pterm_ptr = C_matrix[2*i][2*j].get_pterms();
		
    	while(pterm_ptr)
    	{
			if (braid_control::STUDY_RHO_MAPPING)
				pterm_ptr->n.Qi() = pterm_ptr->n.Qk();

			pterm_ptr->n.Qj() = It(0);
			pterm_ptr->n.Qk() = It(0);
			pterm_ptr = pterm_ptr->next;
    	}
    	sanitize(&C_matrix[2*i][2*j]);

		/* now copy the Z we have just created to the
		   11 relative position and create the conjugate
		*/
		C_matrix[2*i+1][2*j+1] = C_matrix[2*i][2*j];

		pterm_ptr = C_matrix[2*i+1][2*j+1].get_pterms();

	    while(pterm_ptr)
    	{
			pterm_ptr->n.Qi() *= It(-1);
			pterm_ptr = pterm_ptr->next;
    	}

		/* next create W at the relative 01 position */
		C_matrix[2*i][2*j+1] = H_matrix[rperm[i]][cperm[j]];

		/* For the mu-mapping copy j and k to r and i, then set j and k to zero.
		   For the rho-mapping just copy j to r and then set j and k to zero (i is correctly set in this case)
		*/
		pterm_ptr = C_matrix[2*i][2*j+1].get_pterms();
	    while(pterm_ptr)
	    {
			pterm_ptr->n.Qr() = pterm_ptr->n.Qj();
			
			if (!braid_control::STUDY_RHO_MAPPING)
				pterm_ptr->n.Qi() = pterm_ptr->n.Qk();
				
			pterm_ptr->n.Qj() = It(0);
			pterm_ptr->n.Qk() = It(0);
			pterm_ptr = pterm_ptr->next;
	    }
	    sanitize(&C_matrix[2*i][2*j+1]);

		/* finally copy the W we have just created to the
		   10 relative position and create the negative
		   conjugate - this ammounts to negating the r term
		*/
		C_matrix[2*i+1][2*j] = C_matrix[2*i][2*j+1];
		pterm_ptr = C_matrix[2*i+1][2*j].get_pterms();
	    while(pterm_ptr)
    	{
			pterm_ptr->n.Qr() *= It(-1);
			pterm_ptr = pterm_ptr->next;
    	}
    }

	if (clean_up_params)
	{
		delete[] rperm;
		delete[] cperm;
	}
	
	return C_matrix;
	
}


/* R_study_matrix calculates the image of C_matrix under the mapping from M_n(C)
   to M_{2n}(R).
*/
matrix<polynomial<scalar,int>,scalar> R_study_matrix(const matrix<polynomial<scalar,int>,scalar>& C_matrix) {return C_matrix;}

template <class T, class U, class St> 
matrix<polynomial<T,U>,St> R_study_matrix(const matrix<polynomial<T,U>,St>& C_matrix)
{
	int n = C_matrix.numcols();

	typedef int It;

    matrix<polynomial<T,U>, St> R_matrix = matrix<polynomial<T,U>,St>(2*n,2*n);

    for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
    {
		/* take the complex polynomial at C_matrix[i][j]
		   and create the 2x2 matrix   X   Y
								      -Y   X
		   at R_matptr[2*i][2*j].
		
		   Here X is the polynomial formed by setting all complex terms in 
		   coeficients of C_matrix[i][j] to zero, and Y is formed by
		   setting all the real components to the respective values of the 
		   complex terms, then setting the complex terms to zero.
		*/
		R_matrix[2*i][2*j] = C_matrix[i][j];
		
		/* move along the pterms setting the i components of the coefficients to
		   zero, the j and k components are already zero.  This is the 00
		   position in relation to the base, 2*i,2*j
		*/
		pterm<T,U>* pterm_ptr = R_matrix[2*i][2*j].get_pterms();
    	while(pterm_ptr)
    	{
			pterm_ptr->n.Qi() = It(0);
			pterm_ptr = pterm_ptr->next;
    	}
    	sanitize(&R_matrix[2*i][2*j]);

		/* now copy the X we have just created to the
		   11 relative position
		*/
		R_matrix[2*i+1][2*j+1] = R_matrix[2*i][2*j];

		/* next create Y at the relative 01 position */
		R_matrix[2*i][2*j+1] = C_matrix[i][j];

		/* copy i to r then delete i */
		pterm_ptr = R_matrix[2*i][2*j+1].get_pterms();
	    while(pterm_ptr)
	    {
			pterm_ptr->n.Qr() = pterm_ptr->n.Qi();
			pterm_ptr->n.Qi() = It(0);
			pterm_ptr = pterm_ptr->next;
	    }
	    sanitize(&R_matrix[2*i][2*j+1]);
		
		/* finally copy the Y we have just created to the
		   10 relative position and negate the r term
		*/
		R_matrix[2*i+1][2*j] = R_matrix[2*i][2*j+1];
		pterm_ptr = R_matrix[2*i+1][2*j].get_pterms();
	    while(pterm_ptr)
    	{
			pterm_ptr->n.Qr() *= It(-1);
			pterm_ptr = pterm_ptr->next;
    	}
    }
	
	return R_matrix;
}

/* delta_1 evaluates the gcd of the generators of the first elementary ideal, E_1,
   determined by the matrix M.  It is able to work with the intermediate
   ideals determined by mapping from M_{n,m}(H) into M_{2n,2m}(C) and evaluating codimension
   1 or 2 determinants of the image in M_{2n,2m}(C), providing M is given as 
   this image at the time of the call.  
*/
template <class T, class St> T delta_1(const matrix<T,St>& M, int codimension, string title, T delta_0 = T("0"))
{
	T hcf=delta_0; // needed for the non-commutative case where delta_0 is required as a generator of E_1
    Laurent_scale (hcf);

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "delta_1: hcf initialized to " << hcf << endl;
	
	int matrix_rows = M.numrows();
	int matrix_cols = M.numcols();
	
    int rperm[matrix_rows];
    int cperm[matrix_cols];
	
	if (codimension == 1)
	{

if (braid_control::DEBUG >= braid_control::SUMMARY)
    debug << "delta_1: codimension 1 generators of 1st elementary ideal, E_1:" << endl;

		for (int r1=0;r1<matrix_rows;r1++)
    	{
			if (!braid_control::SILENT_OPERATION && matrix_control::WAIT_INFO && matrix_control::COMFORT_DOTS)
			{
				cout << ".";
				cout.flush();
			}

			/* take out the r1-th row */
			for (int i=0; i<r1; i++)
			    rperm[i] = i;
			for (int i=r1+1; i<matrix_rows; i++)
			    rperm[i-1] = i;

			for (int c1=0;c1<matrix_cols;c1++)
			{
			    /* take out the c1-th column */
		    	for (int i=0; i<c1; i++)
				cperm[i] = i;
					
				if (matrix_cols > matrix_rows)
				{
					/* We're dealing with a long knot and there is one extra column than row */
					for (int c2= c1+1; c2< matrix_cols; c2++)
					{
						/* take out the c2-nd column */
						for (int i=c1+1; i<c2; i++)
						    cperm[i-1] = i;
					
			    		for (int i=c2+1; i<matrix_cols; i++)
							cperm[i-2] = i;
						
						if (minor_determinant(hcf, M, rperm, cperm, codimension, r1, 0, c1, c2, 0, title))
							goto delta_1_calculation_complete;
						
					}
				}
				else
				{				
			    	for (int i=c1+1; i<matrix_cols; i++)
					cperm[i-1] = i;

					if (minor_determinant(hcf, M, rperm, cperm, codimension, r1, 0, c1, 0, 0, title))
						goto delta_1_calculation_complete;
				}
			}
		}
	}
	else if (codimension == 2)
	{

if (braid_control::DEBUG >= braid_control::SUMMARY)
    debug << "delta_1: codimension 2 generators of intermediate 1st elementary ideal, E_1:" << endl;
	
		for (int r1=0;r1<matrix_rows;r1++)
    	{

			if (!braid_control::SILENT_OPERATION && matrix_control::WAIT_INFO && matrix_control::COMFORT_DOTS)
			{
				cout << ".";
				cout.flush();
			}

			/* take out the r1-th row */
			for (int i=0; i<r1; i++)
			    rperm[i] = i;
			
			for (int r2= r1+1; r2< matrix_rows; r2++)
			{
				/* take out the r2-nd row */
				for (int i=r1+1; i<r2; i++)
				    rperm[i-1] = i;
				
				for (int i=r2+1; i<matrix_rows; i++)
			    	rperm[i-2] = i;

				for (int c1=0;c1<matrix_cols;c1++)
				{
				    /* take out the c1-th column */
			    	for (int i=0; i<c1; i++)
						cperm[i] = i;
					
					for (int c2= c1+1; c2< matrix_cols; c2++)
					{
						/* take out the c2-nd column */
						for (int i=c1+1; i<c2; i++)
						    cperm[i-1] = i;
					
						if (matrix_cols > matrix_rows)
						{
							/* We're dealing with a long knot and there is one extra column than row */
							for (int c3= c2+1; c3< matrix_cols; c3++)
							{
								/* take out the c3-rd column */
								for (int i=c2+1; i<c3; i++)
						    		cperm[i-2] = i;
					
			    				for (int i=c3+1; i<matrix_cols; i++)
									cperm[i-3] = i;
						
								if (minor_determinant(hcf, M, rperm, cperm, codimension, r1, r2, c1, c2, c3, title))
									goto delta_1_calculation_complete;
							}
						}
						else
						{
				    		for (int i=c2+1; i<matrix_cols; i++)
								cperm[i-2] = i;
						
							if (minor_determinant(hcf, M, rperm, cperm, codimension, r1, r2, c1, c2, 0, title))
								goto delta_1_calculation_complete;
						}

					}
				}
			}
		}
    }
	else
	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\nError! Unsupported codimension parameter presented to function delta_1" << endl;
if (braid_control::DEBUG >= braid_control::SUMMARY)
    debug << "delta_1: Error! Unsupported codimension parameter" << endl;
		exit(0);
	}

delta_1_calculation_complete:
	return hcf;
}

/* minor determinant is a macro-like function that is only used because delta_1
   needs to cope with codimension 1 and 2 cases for both square and non-square matrices, 
   so this function provides common code for all cases.  The parameters at the call
   are therefore essentially the same as those for delta_1.
   
   The function was extended in version 12.0 to return a boolean indicating whether or not the
   minor determinant it calculates is a unit.  If so then the function returns true, otherwise
   it returns false.
*/
template <class T, class St> bool minor_determinant(T& hcf, const matrix<T,St>& M, int* rperm, int* cperm,
                       int codimension, int r1, int r2, int c1, int c2, int c3, string title)
{
	/* rperm and cperm will always indicate a square submatrix of M but there may be more columns than rows
	   in M, so we use the number of rows as the matrix_size to pass to the determinant functions*/	
    int matrix_size = M.numrows();
	
	T delta;

if (braid_control::DEBUG >= braid_control::DETAIL)
{
    debug << "minor_determinant: matrix_size = " << matrix_size << ", codimension = " << codimension
          << ", r1=" << r1 << ", r2=" << r2 << ", c1=" << c1 << ", c2=" << c2 << ", c3=" << c3 << endl;;
    debug << "minor_determinant: rperm: ";
	for ( int k = 0 ; k< matrix_size-codimension; k++)
    	debug << rperm[k] << " ";
	debug << endl;
    debug << "minor_determinant: cperm: ";
	for ( int k = 0 ; k< matrix_size-codimension; k++)
    	debug << cperm[k] << " ";
	debug << endl;

}
	

   	if (!braid_control::QUATERNION && matrix_size == 2)
		delta = M[rperm[0]][cperm[0]];
   	else
   	{
//		if (!braid_control::QUATERNION || braid_control::COMPLEX_STUDY_DELTA_1 || braid_control::REAL_STUDY_DELTA_1)
		if (!braid_control::QUATERNION || braid_control::COMPLEX_STUDY_DELTA_1)
		{
			/* we want the ordinary determinant of the given matrix,
			   either because we are working with the Alexander rep or
			   because the matrix contains complex or real entries
			   determined in the course of evaluating the Study versions of Delta_1
			*/	
    			delta = determinant(M,title,matrix_size-codimension,rperm,cperm);
			
		}
		else
		{
			/* QUATERNION && !COMPLEX_STUDY_DELTA_1 && !REAL_STUDY_DELTA_1,
			   so take study_determinant, M will be quaternionic matrix rep in this case*/
    		delta = study_determinant(M,title,matrix_size-codimension,rperm,cperm);
		}
   	}

    if (braid_control::EXTRA_OUTPUT)
   	{
		if (codimension == 1)
		{
			if (!braid_control::SILENT_OPERATION)
			{
				cout << "\nGenerator (" << r1 << "," << c1;
				if (braid_control::LONG_KNOT)
					cout << "&" << c2;
				cout << ") = " << delta;
			}

			output << "\nGenerator (" << r1 << "," << c1;
			if (braid_control::LONG_KNOT)
				output << "&" << c2;
			output << ") = " << delta;
		}
		else
		{
			if (!braid_control::SILENT_OPERATION)
			{
				cout << "\nGenerator (" << r1 << "&" << r2 << "," << c1 << "&" << c2;
				if (braid_control::LONG_KNOT)
					cout << "&" << c3;
				cout << ") = " << delta;
			}

			output << "\nGenerator (" << r1 << "&" << r2 << "," << c1 << "&" << c2;
			if (braid_control::LONG_KNOT)
				output << "&" << c3;
			output << ") = " << delta;
		}
   	}

if (braid_control::DEBUG >= braid_control::SUMMARY)
{
	if (codimension == 2)
		debug << "minor_determinant: generator (" << r1 << " & " << r2 << "," << c1 << " & " << c2 << ") = " << delta << endl;
	else
    	debug << "minor_determinant: generator (" << r1 << "," << c1 << ") = " << delta << endl;
}
    if (!braid_control::QUATERNION)
   	{
		/* set s=1 in delta if it involves s */
		if (delta.numvars() && strchr(delta.getvars(),'s'))
		    set_to_one(delta,'s');

		if (braid_control::EXTRA_OUTPUT)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\nwith s=1: " << delta;
		    output << "\nwith s=1: " << delta;
		}
   	}

    /* Now delta is a Laurent polynomial in one variable,
       scale it so that it becomes an ordinary polynomial
    */
    Laurent_scale (delta);

   	if (braid_control::EXTRA_OUTPUT)
    {
		if (!braid_control::SILENT_OPERATION)
			cout << "  .=. " << delta;
		output << "  .=. " << delta;
   	}

    sanitize(&delta);
	
	if (braid_control::DELTA_1_UNIT_CHECK && delta.isunit())
	{
		hcf = T("1");
		
	   	if (braid_control::EXTRA_OUTPUT)
    	{
			if (!braid_control::SILENT_OPERATION)
				cout << "\nunit generator detected, terminating calculation" << endl;
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "unit generator detected, terminating calculation" << endl;
			}
   		}
		
if (braid_control::DEBUG >= braid_control::SUMMARY)
    debug << "minor_determinant: unit generator detected, terminating calculation" << endl;

		return true; // unit generator encountered
	}
	else if (delta != T("0"))
	{
	    hcf = gcd(hcf,delta);

if (braid_control::DEBUG >= braid_control::SUMMARY)
    debug << "minor_determinant: gcd stands at " << hcf << endl;

		if (braid_control::DELTA_1_UNIT_CHECK && hcf.isunit())
		{
		
		   	if (braid_control::EXTRA_OUTPUT)
	    	{
				if (!braid_control::SILENT_OPERATION)
					cout << "\nunit gcd detected, terminating calculation" << endl;
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "unit gcd detected, terminating calculation" << endl;
				}
	   		}
		
if (braid_control::DEBUG >= braid_control::SUMMARY)
    debug << "minor_determinant: unit gcd detected, terminating calculation" << endl;

			return true; // unit generator encountered
		}

	}
	else
	{
if (braid_control::DEBUG >= braid_control::SUMMARY)
    debug << "minor_determinant: generator is zero, no change in gcd" << endl;
	}
	
	return false;
}

Qpmatrix adj2(Qpmatrix& M)
{
	Qpmatrix result(2,2);
	result[0][0] = M[1][1];
	result[0][1] = M[0][1]*Qpolynomial("-1");
	result[1][0] = M[1][0]*Qpolynomial("-1");
	result[1][1] = M[0][0];
	return result;
}

void simplify_matrix (Qpmatrix& M)
{
	for (size_t i = 0; i < M.numrows(); i++)
	for (size_t j = 0; j < M.numcols(); j++)
    	M[i][j] = simplify_rational(M[i][j]);
}




/* normalize takes a (Delta_i^H) Hpolynomial, Laurent scales it, ensures that the
   constant term is positive, then multiplies or divides by |A|^2, |B|^2, |C|^2, 
   |D|^2, or |Delta|^2 so that the constant term is as clost to 1 a possible.
   This function needs interactive help from the user!!!!!
*/
template <class T1, class U, class T2, class St>  void normalize(polynomial<T1,U>& poly, const matrix<T2,St>& switch_matrix) {}

void normalize(Hpolynomial& poly, const Hpmatrix& switch_matrix)
{
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "\nnormalize: normalizing " << poly << endl;

	if (poly.get_pterms() == NULL)
	{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "\nnormalize: doing nothing, poly is zero\n" << endl;
		return;
	}
	
	Laurent_scale(poly);
	pterm<Quaternion,int>* pterm_ptr = poly.get_pterms();

	if (pterm_ptr->n.Qr() < 0)
		poly *= Hpolynomial("-1");

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "\nnormalize: scaled poly = " << poly << endl;

	/* if the constant term is 1 already, there's nothing to do */
	if (pterm_ptr->n == Quaternion(1))
	{	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "\nnormalize: constant term is already 1, nothing to do" << endl;			
		return;
	}
	
	Quaternion A = switch_matrix[0][0].get_pterms()->n;
	Quaternion B = switch_matrix[0][1].get_pterms()->n;
	Quaternion C = switch_matrix[1][0].get_pterms()->n;
	Quaternion D = switch_matrix[1][1].get_pterms()->n;

	Hpolynomial factor[5];
	factor[0] = Hpolynomial(norm_sqrd(A));
	factor[1] = Hpolynomial(norm_sqrd(B));
	factor[2] = Hpolynomial(norm_sqrd(C));
	factor[3] = Hpolynomial(norm_sqrd(D));
	Quaternion delta = conj(B)*A*Quaternion(norm_sqrd(D)) - conj(D)*C*Quaternion(norm_sqrd(B));
	factor[4] = Hpolynomial(norm_sqrd(delta));
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "normalize: initial factors : ";
	for (int i=0; i<5; i++)
		debug << factor[i] << ' ';
	debug << endl;
}

	/* we're not interested in factors of 1, nor repeats */
	int num_factors = 5;
	for (int i=0; i<5; i++)
	{
		if (factor[i] == Hpolynomial("1"))
		{
			factor[i] = Hpolynomial("0");
			num_factors--;
		}
	}
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "normalize: zeroed factors : ";
	for (int i=0; i<5; i++)
		debug << factor[i] << ' ';
	debug << endl;
}

	if (num_factors == 0)
	{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "normalize: all factors are 1, we're done." << endl;
		return;
	}
	else
	{
		/* move factors left so that non repeated factor are flushed left in factor */
		for (int i=0; i<num_factors; i++)
		{
			/* find the next non-zero factor */
			for (int j=i; j<5; j++)
			{
				if (factor[j] != Hpolynomial("0"))
				{
					/* check it's not a duplicate */
					for (int k=0; k<i; k++)
					{
						if (factor[j] == factor[k])
						{
							factor[j] = Hpolynomial("0");
							num_factors--;
							break;
						}
					}
					
					/* if it's still non-zero, move it to the next location */
					if (factor[j] != Hpolynomial("0"))
					{
						factor[i] = factor[j];
					
						/* set factor[j] to zero unless we've not moved it at all! */
						if (j != i)
							factor[j] = Hpolynomial("0");
						
						break;
					}
				}
			}
		}		
	}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{	
	debug << "normalize: available factors ";
	for (int i=0; i<num_factors; i++)
		debug << i+1 << ": " << factor[i] << "    ";
	debug << endl;
}

	Hpolynomial cfactor = Hpolynomial("1"); // used for debugging
	int rfactor;
	do
	{
		cout << "\nPolynomial stands at " << poly << endl;	
		cout << "Available factors ";
		for (int i=0; i<num_factors; i++)
			cout << i+1 << ": " << factor[i] << "    ";
		cout << endl;
		
		cout << "\nEnter factor number you want to apply, zero if you're done";
		cout << "\nfactor: ";
		cin >> rfactor;
			
		if (rfactor)
		{
			Hpolynomial& pfactor = factor[rfactor-1];
			
			cout << "\nMultiply or divide? [m/d] ";
		
			char multiply_or_divide;
			cin >> multiply_or_divide;
			
			if (multiply_or_divide == 'd')
			{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "\nnormalize: dividing by ";
				cfactor /= pfactor;
				poly /= pfactor;
			}
			else
			{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "\nnormalize: multiplying by ";
				cfactor *= pfactor;
				poly *= pfactor;			
			}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << pfactor;
			
			pterm_ptr = poly.get_pterms(); // the pterms may have been changed
			if (pterm_ptr->n == Quaternion(1))
			{	
				cout << "\nConstant term now 1" << endl;
				rfactor = 0;
			}
		}
	} while (rfactor);

if(braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "\nnormalize: total multiplicative factor applied = " << cfactor;

	cout << "\nNormalized poly = " << poly << endl;	
}

template <class T> void equality_test (T x, T y, bool A)
{
	if (x == y)
	{
		if (!braid_control::SILENT_OPERATION)
		{
			if (A)
				cout << "\nSwitch elements A and D are equal" << endl;
			else
				cout << "Switch elements B and C are equal" << endl;
		}
		
		if (!braid_control::RAW_OUTPUT)
		{
			if (A)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Switch elements A and D are equal" << endl;
			}
			else
			{
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Switch elements B and C are equal" << endl;
			}
		}

if (braid_control::DEBUG >= braid_control::BASIC)
{
	if (A)
		debug << "braid::main: switch elements A and D are equal" << endl;
	else
		debug << "braid::main: switch elements B and C are equal" << endl;
}
	}
	else
	{
		if (!braid_control::SILENT_OPERATION)
		{
			if (A)
				cout << "\nSwitch elements A and D are not equal" << endl;
			else
				cout << "Switch elements B and C are not equal" << endl;
		}
		
		if (!braid_control::RAW_OUTPUT)
		{
			if (A)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Switch elements A and D are not equal" << endl;
			}
			else
			{
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Switch elements B and C are not equal" << endl;
			}
		}

if (braid_control::DEBUG >= braid_control::BASIC)
{
	if (A)
		debug << "braid::main: switch elements A and D are not equal" << endl;
	else
		debug << "braid::main: switch elements B and C are not equal" << endl;	
}
	}
}

template <class T, class St> 
bool braid_rep(matrix<T,St>*& Matrix_rep, const matrix<T,St>& switch_matrix, const matrix<T,St>& switch_matrix_inverse, string input_string,
                int num_terms, int num_strings)
{
	bool calculate_polys;
	
	matrix<T,St>*		Twist;
	int switch_matrix_N_factor = switch_matrix.numcols()/2;
	char* inbuf = c_string(input_string);
	
	/* create the matrix corresponding to virtual crossings */
	if (switch_matrix_N_factor == 1)
	{
		Twist = new matrix<T,St>(2,2);
		matrix<T,St>& twist = *Twist;
		
		twist[0][0] = T("0");  //T is a polynomial type, so has a string constructor
		twist[0][1] = T("1");
		twist[1][0] = T("1");
		twist[1][1] = T("0");
	}
	else
	{
		Twist = new matrix<T,St>(2*switch_matrix_N_factor,2*switch_matrix_N_factor);
		matrix<T,St>& twist = *Twist;
	
		for (int i=0; i< switch_matrix_N_factor; i++)
		for (int j=0; j< switch_matrix_N_factor; j++)
			twist[i][j] = T("0");
		
		for (int i=0; i< switch_matrix_N_factor; i++)
		{
			for (int j=0; j<i; j++)
				twist[i][j+switch_matrix_N_factor] = T("0");
				
			twist[i][i+switch_matrix_N_factor] = T("1");
			
			for (int j=i+1; j<switch_matrix_N_factor; j++)
				twist[i][j+switch_matrix_N_factor] = T("0");
		}

		for (int i=0; i< switch_matrix_N_factor; i++)
		{
			for (int j=0; j<i; j++)
				twist[i+switch_matrix_N_factor][j] = T("0");
				
			twist[i+switch_matrix_N_factor][i] = T("1");
			
			for (int j=i+1; j<switch_matrix_N_factor; j++)
				twist[i+switch_matrix_N_factor][j] = T("0");
		}

		for (int i=0; i< switch_matrix_N_factor; i++)
		for (int j=0; j< switch_matrix_N_factor; j++)
			twist[i+switch_matrix_N_factor][j+switch_matrix_N_factor] = T("0");
	}

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "\nbraid_rep: twist matrix:\n" << *Twist << endl;
}

	/* record whether this is a virtual knot or a classical one 
	if (strchr(inbuf, 't') || strchr(inbuf,'T'))
		virtual_crossings_present = true;
	*/
	
	int& N = switch_matrix_N_factor;
	
	matrix<T,St> identity(N*num_strings,N*num_strings);	
	for (int i=0; i<N*num_strings; i++)
	{
		for (int j=0; j<i; j++)
			identity[i][j] = T("0");
	
		identity[i][i] = T("1");
	
		for (int j=i+1; j<N*num_strings; j++)
			identity[i][j] = T("0");
	}
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "\nbraid_rep: identity:\n" << identity;


    /* Create the matrix and initialize to the identity */
    Matrix_rep = new matrix<T,St>(N*num_strings,N*num_strings);
	matrix<T,St>& matrix_rep = *Matrix_rep;
	matrix_rep = identity;		

    /* Now evaluate the product of the matrices determined by
       the word in inbuf.
    */
	bool S;
	bool inverse;
    char* cptr = inbuf;
    for ( int i = 0; i< num_terms; i++)
    {
		if (*cptr == '-')
		{
		    inverse = true;
	    	cptr++;
		}
		else
		    inverse = false;

		if (*cptr == 's' || *cptr == 'S')
		    S = true;
		else
		    S = false;
		cptr++;

		char* mark = cptr; /* mark where we start the number */

		/* look for the end of the number */
		while (isdigit(*cptr))
		    cptr++;

		int number;
		get_number(number, mark);

		matrix<T,St> term = identity;

		const matrix<T,St>* Mptr;
		
		if (S && inverse)
		{
			Mptr= &switch_matrix_inverse;
		}
		else if (S)
		{
			Mptr = &switch_matrix;
		}
		else /* use T */
		{
			Mptr = Twist;
		}
				
		/* determine which form of Sn we are using Sn_matrix_top_down == true iff
           Sn = I^{n-1} x S x I^{k-n-1}, otherwise Sn = I^{k-n-1} x S x I^{n-1}
		*/
		size_t base;
		if (braid_control::Sn_matrix_top_down)
			base = (number-1)*N; // offsets from 0, braid strings from 1
		else
			base = term.numcols() - switch_matrix.numcols() - (number-1)*N;  // 
			                                                                      
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "\nbraid_rep: base for S in term matrix is " << base << endl;																					  																					   
			
		size_t switch_size = switch_matrix.numcols();		
		for (size_t j=0; j < switch_size; j++)
		for (size_t k=0; k < switch_size; k++)
		{
			term[base+j][base+k] = (*Mptr)[j][k];
		}

		if (!braid_control::SILENT_OPERATION && braid_control::WAIT_SWITCH)
			cout << "term " << i+1 << endl;
			
		matrix_rep *= term;

if (braid_control::DEBUG >= braid_control::DETAIL)
{
    debug << "\nbraid_rep: term for right multiplication by ";
    if (inverse)
        debug << "-";

    if (S)
        debug << "s";
    else
        debug << "t";
    debug << number << ":\n";
    debug << term << endl;

}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
    debug << "\nbraid_rep: matrix after right multiplication by ";
    if (inverse)
        debug << "-";

    if (S)
        debug << "s";
    else
        debug << "t";
    debug << number << ":\n";
    debug << matrix_rep << endl;
}

	}

	if (braid_control::EXTRA_OUTPUT)
	{
		if (!braid_control::SILENT_OPERATION)
		{
			cout << "\nMatrix representation of braid word, M:\n";
			cout << matrix_rep << endl;
		}
		output << "\nMatrix representation of braid word:\n";
		output << matrix_rep << endl;
	}

if (braid_control::DEBUG >= braid_control::BASIC)
{
    debug << "\nbraid_rep: matrix representation of braid word, M:\n";
    debug << matrix_rep << endl;
}

	if (matrix_rep == identity)
	{
    	calculate_polys = false;
		
		if (!braid_control::SILENT_OPERATION)
			cout << "\nmatrix representation of braid word is the identity.\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << "matrix representation of braid word is the identity.\n";
		}
	}
	else
	{
		matrix_rep -= identity;
    	calculate_polys = true;

if (braid_control::DEBUG >= braid_control::DETAIL)
{
    debug << "\nbraid_rep: matrix representation minus the identity, M-I:\n";
    debug << matrix_rep << endl;
}
	}

	delete Twist;	
	delete[] inbuf;

	return calculate_polys;
}

template <class T, class St> 
void R_module_rep(matrix<T,St>*& Matrix_rep, const matrix<T,St>& switch_matrix, const matrix<T,St>& switch_matrix_inverse, 
               string input_string)
{
	bool OLD_SWITCH_ACTION = false;

	generic_code_data code_data;
	
	if (input_string.find('[') != string::npos)
		read_peer_code(code_data, input_string);
	else
		read_immersion_code(code_data, input_string);
	
	matrix<int>& code_table = code_data.code_table; 
	int num_crossings = code_data.num_crossings;
	int num_components = code_data.num_components;
	int num_classical = num_crossings;

	for (int i=0; i<num_crossings; i++)
	{
		if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
			num_classical--;
	}

	/* If num_classical == 0 we have been given a code containing only virtual crossings.
	   This yields a 0x0 matrix representation and a degenerate R-module.  This situation 
	   is indicated by setting the Matrix_rep pointer to zero.
	*/
	if (num_classical == 0)
	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\nError! Labelled peer and immersion codes must contain at least one real crossing" << endl;
		Matrix_rep = 0;
		return;
	}

	/* Determine the size of the R_module representation.  In the case of a knot or link, this is a 
	   2*N*num_classical square matrix where the switch_matrix elements are NxN matrices (for the 
	   Alexander and Quaternionic reps, N=1).  In the case of a long knot, the R_module representation
	   is a 2*N*num_classical x (2*N*num_classical + N) matrix, allowing for the extra variable.  
	*/
	int N = switch_matrix.numrows()/2; // N is used later, it is what is referred to elsewhere as the switch_matrix_N_factor
	int matrix_N_size = 2*num_classical;  // number of N-rows
	int matrix_size = N*matrix_N_size;  // number of actual rows

if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "R_module_rep: num_components = " << num_components << endl;
	debug << "R_module_rep: num_crossings = " << num_crossings << endl;
	debug << "R_module_rep: num_classical = " << num_classical << endl;
	debug << "R_module_rep: N (switch_matrix_N_factor) = " << N << endl;
	debug << "R_module_rep: matrix_N_size = " << matrix_N_size << endl;
	debug << "R_module_rep: matrix_size = " << matrix_size << endl;
}

	/* The matrix representation changed in version 12 of the braid programme to a more natural, but essentially equivalent,
	   formulation.  The code supports both the old and new format via the boolean variable OLD_SWITCH_ACTION to invoke the 
	   original action of a switch at a real crossing.
	   
	   The original switch action was determined from the equations
	
		A  B  x_2i    = x_2j         A  B  x_{2j-1}  = x_{2i+1}   
		C  D  x_{2j-1}  x_{2i+1}     C  D  x_2i        x_2j 
	
	   for type I and type II crossings respectively, where the matrix A B
																	   C D
	
	   is either S, S' or T, depending upon the label of the crossing.

	   This means that the action of the switch is such that for a positive crossing with labels

                    b d
				    a c

	   we have c=Aa+Bb d=Ca+Db, which is why the above equations have the  column vector (x_2i, x_{2j-1}) the "wrong way around" 
	   for the corresponding crossing Type.

	   Thus for the old switch action we have equations
	
	   A * x_2i + B * x_{2j-1} - x_2j = 0 and 
	   C * x_2i + D * x_{2j-1} - x_{2i+1} = 0
	
	   for type I crossings, and for type II crossings

	   A * x_{2j-1} + B * x_2i - x_{2i+1} = 0 and 
	   C * x_{2j-1} + D * x_2i - x_2j = 0
	   

	   The new switch action was determined from the equations
	
		A  B  x_{2j-1}  = x_{2i+1}         A  B  x_2i      = x_2j   
		C  D  x_2i        x_2j             C  D  x_{2j-1}    x_{2i+1} 
	
	   for type I and type II crossings respectively, using the same notation as above.  Here the switch action is more natural
	   in the sense that the action is determined by 'standard' matrix multiplication of a vector.  
	   
	   Thus for the new switch action we have equations
	
	   A * x_{2j-1} + B * x_2i - x_{2i+1} = 0 and 
	   C * x_{2j-1} + D * x_2i - x_2j = 0
	
	   for type I crossings, and for type II crossings

	   A * x_2i + B * x_{2j-1} - x_2j = 0 and 
	   C * x_2i + D * x_{2j-1} - x_{2i+1} = 0

	   Note that in effect the change in the switch action is just inverting the definition of a type I and type II crossing.  This is 
	   how the code handles the difference, governed by the OLD_SWITCH_ACTION variable, which may then be controlled by the command
	   line at a later date, if required.


	   This reconciles with the theory in the paper by noting that the theory ignores virtual crossings, so here the above equations
	   state that x_2i = x_{2i+1} and x_{2j-1} = x_2j at a virtual crossing.  To reduce the size of the determinant we need to deal with, we remove
	   any such variable identities expressed by virtual crossings are handled by using a set of trivial variable flags that indicate 
	   when a variable is identical to another, we then use a variable permutation to map to a smaller set of variables.
   
	   NB. In the general case the switch matrix may have N x N matrix entries A B C and D (where N is the switch_matrix_N_factor), 
	   in which case the variables described here are N-tuples.
	*/

	valarray<int> trivial_variable_flags(2*num_crossings); 
	valarray<int> variable(2*num_crossings);               
	
	
	for (int i=0; i<2*num_crossings; i++)
		trivial_variable_flags[i] = 0;
	
	for (int i=0; i<num_crossings; i++)
	{
		if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
		{
			/* The equations for this crossing state x_2i = x_{2i+1} and x_{2j-1} = x_2j, where the terminating
			   edges at the crossing are 2i and 2j-1.  Thus, the variables corresponding to the incoming edges
			   are trivial.
			*/
//			int twice_pi_minus1 = ((2*code_table[PERM][i] -1)+(2*num_crossings))%(2*num_crossings);
//if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
//	debug << "R_module_rep: crossing " << i << " trivial variables = " << 2*i << " and " << twice_pi_minus1 << endl;
//			trivial_variable_flags[2*i] = 1;
//			trivial_variable_flags[twice_pi_minus1] = 1;
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "R_module_rep: crossing " << i << " trivial variables = " << 2*i << " and " << code_table[OPEER][i] << endl;
			trivial_variable_flags[2*i] = 1;
			trivial_variable_flags[code_table[OPEER][i]] = 1;
		}
	}
	
if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "R_module_rep: trivial_variable_flags: ";
	for (int i=0; i<2*num_crossings; i++)
		debug << trivial_variable_flags[i] << ' ';
	debug << endl;
}	

	/* now we can determine the variable permutation between the variables
	   assigned to the semi-arcs between all crossings, and those assigned
	   to semi-arcs between classical crossings, ignoring virtual crossings.
	*/
	int count = 0;
	vector<int>& first_edge_on_component = code_data.first_edge_on_component;
	vector<int>& num_component_edges = code_data.num_component_edges;
	
	for (int i=0; i< num_components; i++)
	{
		int first_edge = first_edge_on_component[i];
		for (int j=0; j < num_component_edges[i]; j++)
		{
			variable[first_edge+j] = first_edge+j-count;
			count += trivial_variable_flags[first_edge+j];
		}
		
		/* We have to deal with the case that the last k crossings on the 
		   component are virtual crossing, in which case the last k variables 
		   on the component are the same as the variable corresponding to the first edge.
		   In the degenerate case where all the crossings on a component are virtual
		   this will result in all the edges of the component being associated with 
		   the same variable, being the variable belonging to the first edge of the 
		   next component.  This won't matter as we'll not be creating any relations
		   for these virtual crossings.
		*/
		for (int j=num_component_edges[i]-1; j>=0 && trivial_variable_flags[first_edge+j] == 1; j--)
			variable[first_edge+j] = variable[first_edge];
		
	}
	
//	for (int i=0; i<2*num_crossings; i++)
//	{
//		variable[i] = i-count;
//		count += trivial_variable_flags[i];
//	}
	
//if (braid_control::DEBUG >= braid_control::DETAIL)
//{
//	debug << "R_module_rep: intermediate variables: ";
//	for (int i=0; i<2*num_crossings; i++)
//		debug << variable[i] << ' ';
//	debug << endl;
//}	

//	if (trivial_variable_flags[2*num_crossings-1])
//	{
		/* x_{2num_crossings-1} = x_{2num_crossings} = x_0.  Since we may have also had x_{2num_crossings-2} = 
		   x_{2num_crossings-1} etc.  We may avoid looking back around the knot for the last real crossing and 
		   simply reduce all variable values modulo 2*num_classical.
		*/
//		for (int i=0; i<2*num_crossings; i++)
//		{
//			variable[i] %= 2*num_classical;
//		}
//	}
	
if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "R_module_rep: variables: ";
	for (int i=0; i<2*num_crossings; i++)
		debug << variable[i] << ' ';
	debug << endl;
}	
    /* Now create the matrix and evaluate the representation of the knot based upon the reduced variable set.
	   Each real crossing determines two rows of the matrix representation, representing one of the pairs of 
	   equations above, dependent on the crossing type.  We work through the real crossings adding a pair of 
	   rows to the matrix representation for each one.  (Note we are interested in the R-module generated by
	   this matrix, so the order in which we construct the rows is irrelevant.)
    */

	if (braid_control::LONG_KNOT)
	    Matrix_rep = new matrix<T,St>(matrix_size,matrix_size+N);
	else
	    Matrix_rep = new matrix<T,St>(matrix_size,matrix_size);
		
	matrix<T,St>& matrix_rep = *Matrix_rep;
	
	for (unsigned int i=0; i<matrix_rep.numrows(); i++)
	for (unsigned int j=0; j<matrix_rep.numcols(); j++)
		matrix_rep[i][j] = T("0");
	

	int row = 0;
	
	for (int i=0; i<num_crossings; i++)
	{
		/* only consider classical crossings */
		if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
			continue;
		
//		/* first evaluate 2p(i)-1 = 2j-1 = 2*code_table[PERM][i]-1 modulo 2*num_crossings */
//		int twice_pi_minus1 = ((2*code_table[PERM][i] -1)+(2*num_crossings))%(2*num_crossings);

		/* identify the odd_numbered_peer of the naming edge */
//		int odd_numbered_peer = code_table[OPEER][i];
	
		/* now determine the braid crossing type, POSITIVE, NEGATIVE or VIRTUAL */
		int crossing_type;		
		if (  (code_table[TYPE][i] == generic_code_data::TYPE1 && code_table[LABEL][i] == generic_code_data::POSITIVE)
				||(code_table[TYPE][i] == generic_code_data::TYPE2 && code_table[LABEL][i] == generic_code_data::NEGATIVE)
   					)
			crossing_type = braid_crossing_type::NEGATIVE;
		else
			crossing_type = braid_crossing_type::POSITIVE;

if (braid_control::DEBUG >= braid_control::DETAIL)
{	
	debug << "R_module_rep: crossing " << i;
	if ( crossing_type == braid_crossing_type::NEGATIVE )
		debug << " negative";
	else
		debug << " positive";

	debug <<  ", TYPE " << (code_table[TYPE][i] == generic_code_data::TYPE1? "I" : "II") << endl;
	debug << "R_module_rep: ic-variables\tp[i] = ";// << code_table[PERM][i];
//	debug << "\t2i = " << 2*i << "\t\t2p[i]-1 = " << twice_pi_minus1;
	debug << "\t2i = " << 2*i << "\t\t2j-1 = " << code_table[ODD_TERMINATING][i];
	debug << "\t2i+1 = " << 2*i+1 << "\t2p[i] = " << code_table[EVEN_ORIGINATING][i] << endl;
	
	debug << "R_module_rep: mapped to\t\t\t";//v(p[i]) = " << variable[code_table[PERM][i]];
	debug << "\tv(2i) = " << variable[2*i] << "\tv(2p[i]-1) = " << variable[code_table[ODD_TERMINATING][i]];
	debug << "\tv(2i+1) = " << variable[2*i+1] << "\tv(2p[i]) = " << variable[code_table[EVEN_ORIGINATING][i]] << endl;

	debug << "R_module_rep: row " << 2*row << "\nR_module_rep: \t\t";
}		
		/* first do row 2*row.  2*num_classical is the number of N-cols in the matrix_rep for a closed knot, 
		   there is an additional N-col in the case of a long knot but that is handled explicity below.
		*/
		for (int j=0; j<2*num_classical; j++)  
    	{
//			if (j == (code_table[TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[2*i] : variable[twice_pi_minus1]))
			if (j == (code_table[TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[code_table[EVEN_TERMINATING][i]] : variable[code_table[ODD_TERMINATING][i]]))
			{
				/* Record the A entry at j*/
				if (crossing_type == braid_crossing_type::POSITIVE)
				{
if (braid_control::DEBUG >= braid_control::DETAIL) debug << j << ":S00    ";
					set_matrix_N_element(matrix_rep,2*row,j,switch_matrix,0,0,N);
				}
				else
				{
if (braid_control::DEBUG >= braid_control::DETAIL) debug << j << ":iS00   ";
					set_matrix_N_element(matrix_rep,2*row,j,switch_matrix_inverse,0,0,N);				
				}
			}
//			else if (j == (code_table[TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[twice_pi_minus1] : variable[2*i]))
			else if (j == (code_table[TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[code_table[ODD_TERMINATING][i]] : variable[code_table[EVEN_TERMINATING][i]]))
			{
				/* Record the B entry at j */
				if (crossing_type == braid_crossing_type::POSITIVE)
				{
if (braid_control::DEBUG >= braid_control::DETAIL) debug << j << ":S01    ";
					set_matrix_N_element(matrix_rep,2*row,j,switch_matrix,0,1,N);				
				}
				else
				{
if (braid_control::DEBUG >= braid_control::DETAIL) debug << j << ":iS01   ";
					set_matrix_N_element(matrix_rep,2*row,j,switch_matrix_inverse,0,1,N);				
				}
			}
			else 
			{
				/* there is a possibility that we have to handle the case where
				   
					A  B  x_p  = x_p     A  B  x_p  = x_s   
					C  D  x_q    x_r     C  D  x_q    x_q 
				
				   due to a Reidemeister Type 1 portion of the knot encountering
				   two or more virtual crossings.  Such a case is theoretically
				   redundant but cannot be detected algorithmically until we reach 
				   here.
				   
				   We deal with this case by adding the -1 entry to the row
				   rather than setting an entry of the matrix to -1, hence the 
				   following code...
				*/
				
//				if (j == (code_table[TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[2*code_table[PERM][i]] : variable[2*i+1]))
				if (j == (code_table[TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[code_table[EVEN_ORIGINATING][i]] : variable[code_table[ODD_ORIGINATING][i]]))
				{
					/* The theoretical edge variables ignore virtual crossings and have been represented by variable[i] 
					   so that all immersion-semi-arcs after the last real crossing have variable[i]=0.  If this is a
					   long knot then infinity sits on this arc and at the end of the long knot we encounter 
					   variable[i]=0 uniquely at the last real crossing.  Depending on how many virtual crossings there 
					   are left to cross before we head off to infinity, we could find it on the 2*row or 2*row+1 of
					   either a type 1 or type 2 crossing but if it's on 2*row then it will be in this conditional
					   clause when j==0.
					*/
					if (braid_control::LONG_KNOT && j == 0) 
					{
						/* Add -1 to the polynomial at (*matrixptr)[2*row][matrix_N_size] */
if (braid_control::DEBUG >= braid_control::DETAIL) debug << matrix_N_size << ":-1     ";
						decrement_matrix_N_element(matrix_rep,2*row,matrix_N_size,N);
					}
					else
					{
						/* Add -1 to the polynomial at (*matrixptr)[2*row][j] */
if (braid_control::DEBUG >= braid_control::DETAIL) debug << j << ":-1     ";
						decrement_matrix_N_element(matrix_rep,2*row,j,N);
					}
				}
				else
				{
					/* We've already set the location to zero, so record it
					   in the debug file if necessary */
if (braid_control::DEBUG >= braid_control::DETAIL) debug << j << ":0      ";
				}
			}
		}
   
if (braid_control::DEBUG >= braid_control::DETAIL) 
	debug << "\nR_module_rep: row " << 2*row+1 << "\nR_module_rep: \t\t";

		/* now row 2*row+1 */	
		for (int j=0; j<2*num_classical; j++)  
    	{
//			if (j == (code_table[TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[2*i] : variable[twice_pi_minus1]))
			if (j == (code_table[TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[code_table[EVEN_TERMINATING][i]] : variable[code_table[ODD_TERMINATING][i]]))
			{
				/* Record the C entry at j */
				if (crossing_type == braid_crossing_type::POSITIVE)
				{
if (braid_control::DEBUG >= braid_control::DETAIL) debug << j << ":S10    ";
					set_matrix_N_element(matrix_rep,2*row+1,j,switch_matrix,1,0,N);				
				}
				else 
				{
if (braid_control::DEBUG >= braid_control::DETAIL) debug << j << ":iS10   ";
					set_matrix_N_element(matrix_rep,2*row+1,j,switch_matrix_inverse,1,0,N);				
				}
			}
//			else if (j == (code_table[TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[twice_pi_minus1] : variable[2*i]))
			else if (j == (code_table[TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[code_table[ODD_TERMINATING][i]] : variable[code_table[EVEN_TERMINATING][i]]))
			{
				/* Record the D entry at j */
				if (crossing_type == braid_crossing_type::POSITIVE)
				{
if (braid_control::DEBUG >= braid_control::DETAIL) debug << j << ":S11    ";
					set_matrix_N_element(matrix_rep,2*row+1,j,switch_matrix,1,1,N);				
				}
				else 
				{
if (braid_control::DEBUG >= braid_control::DETAIL) debug << j << ":iS11   ";
					set_matrix_N_element(matrix_rep,2*row+1,j,switch_matrix_inverse,1,1,N);				
				}
			}
			else 
			{				
//				if (j == (code_table[TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[2*i+1] : variable[2*code_table[PERM][i]]))
				if (j == (code_table[TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[code_table[ODD_ORIGINATING][i]] : variable[code_table[EVEN_ORIGINATING][i]]))
				{
					/* see comment for corresponding clause above */
					if (braid_control::LONG_KNOT && j == 0)
					{
						/* Add -1 to the polynomial at (*matrixptr)[2*row+1][matrix_N_size] */
if (braid_control::DEBUG >= braid_control::DETAIL) debug << matrix_N_size << ":-1     ";
						decrement_matrix_N_element(matrix_rep,2*row+1,matrix_N_size,N);
					}
					else
					{
						/* Add -1 to the polynomial at (*matrixptr)[2*row+1][j] */
if (braid_control::DEBUG >= braid_control::DETAIL) debug << j << ":-1     ";
						decrement_matrix_N_element(matrix_rep,2*row+1,j,N);
					}
				}
				else
				{
					/* We've already set the location to zero, so record it
					   in the debug file if necessary */
if (braid_control::DEBUG >= braid_control::DETAIL) debug << j << ":0      ";
				}
			}
		}

if (braid_control::DEBUG >= braid_control::DETAIL) 
	debug << endl;
		
		row++;
	}
	

	if (braid_control::EXTRA_OUTPUT)
	{
		if (!braid_control::SILENT_OPERATION)
		{
			cout << "\nMatrix representation of R-module determined by immersion code, M:\n";
			cout << matrix_rep << endl;
		}
		output << "\nMatrix representation of R-module determined by immersion code, M:\n";
		output << matrix_rep << endl;
	}

if (braid_control::DEBUG >= braid_control::BASIC)
{
    debug << "R_module_rep: matrix representation derived from generic code data:" << endl;
    debug << matrix_rep;
    debug << endl;
}

}

/* rat_minor_determinant was introduced for rat_poly_invariant to deal with both long knots and closed knots.  
   It evaluates a minor determinant and if necessary updates the numerator gcd, returning true if it detects a
   unit generator or gcd, and false otherwise.
*/
bool rat_minor_determinant (Qpmatrix* Matrix_rep, int N, int Nr1, int Nc1, int Nc2, int* rperm, int* cperm, string title, polynomial<scalar,bigint>& hcf)
{
	int matrix_rows = Matrix_rep->numrows();
	Qpolynomial det = determinant(*Matrix_rep, title,  matrix_rows - N, rperm, cperm);

if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "rat_minor_determinant: " << (braid_control::LONG_KNOT? "p^{(1)}":"Delta_1") << " generator " << N << "-row/column (" << Nr1 << "," << Nc1;
	if (braid_control::LONG_KNOT)
		debug << "&" << Nc2;
	debug << ")" << endl;
}
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "rat_minor_determinant: rperm: ";
    for (int k = 0 ; k< matrix_rows-N; k++)
		debug << rperm[k] << " ";
	debug << endl;
	debug << "rat_minor_determinant: cperm: ";
    for (int k = 0 ; k< matrix_rows-N; k++)
		debug << cperm[k] << " ";
	debug << endl;
}
					
	if (!braid_control::SILENT_OPERATION)
	{
		cout << (braid_control::LONG_KNOT? "\np^{(1)}":"\nDelta_1") << " generator (" << Nr1 << "," << Nc1;
		if (braid_control::LONG_KNOT)
			cout << "&" << Nc2;
		cout << ") = " << det;
	}
					
	if (!braid_control::RAW_OUTPUT)
	{
		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");

		output << (braid_control::LONG_KNOT? "p^{(1)}":"Delta_1") << " generator (" << Nr1 << "," << Nc1;
		if (braid_control::LONG_KNOT)
			output << "&" << Nc2;
		output << ") = " << det;
	}

if (braid_control::DEBUG >= braid_control::SUMMARY)
{
	debug << "rat_minor_determinant: ";
	debug << (braid_control::LONG_KNOT? "p^{(1)}":"Delta_1") << " generator det = " << det << endl;
}

	polynomial<scalar,bigint> det_numerator = det.getn();
   	sanitize(&det_numerator);
	polynomial<scalar,bigint> det_denominator = det.getd();
   	sanitize(&det_denominator);

	if (braid_control::DELTA_1_UNIT_CHECK && det_denominator.isunit() && det_numerator.isunit())
	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nunit generator detected, terminating calculation" << endl;
		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << "unit generator detected, terminating calculation" << endl;
		}
		
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "rat_minor_determinant: unit generator detected, terminating calculation" << endl;

		if (braid_control::NUMERATOR_GCD)
			hcf = polynomial<scalar,bigint>("1");					
	
		return true; 
	}
	else if (braid_control::NUMERATOR_GCD)
	{
		/******** Take hcf code for single variable examples *********/
    	hcf = gcd(hcf,det_numerator);

		if (!braid_control::SILENT_OPERATION)
			cout << "\n numerator gcd stands at " << hcf << endl;
		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
	    	output << "numerator gcd stands at " << hcf << endl;
		}

if (braid_control::DEBUG >= braid_control::SUMMARY)
   	debug << "rat_minor_determinant:  numerator gcd stands at " << hcf << endl;
	
		if (braid_control::DELTA_1_UNIT_CHECK && hcf.isunit())
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\nunit gcd detected, terminating calculation" << endl;
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "unit gcd detected, terminating calculation" << endl;
			}

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "rat_minor_determinant: unit gcd detected, terminating calculation" << endl;

			return true;; 
		}
	}
	
	return false;
	
}



/* rat_poly_invariant calculates the codimension 0 and codimension 1 invariants determined by the supplied 
   switch matrix and its inverse for the braid or immersion code given in the input string.  
   It is used only where the elements of the switch matrix are of type rational<polynomial<T,U> > and not 
   polynomial<T,U>, i.e. Matrix or Weyl switches. Note that this means that the switch matrix N-factor 
   is greater than one in this function.
*/
void rat_poly_invariant(const Qpmatrix& switch_matrix, const Qpmatrix& switch_matrix_inverse,
                    string input_string, string title)
{
	int 	num_terms;
	int 	num_strings;
	bool    calculate_polys = true;
	bool	virtual_crossings_present = false;
	
	Qpmatrix*		Matrix_rep;
	
	if (title.length())
   	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\n\n" << title << endl;
		if (!braid_control::RAW_OUTPUT)
			output << "\n\n" << title;
   	}		

	if (!braid_control::RAW_OUTPUT)
	{
		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
		output << input_string;		
	}

	/* note the N-factor of the switch matrix */
	int switch_matrix_N_factor = switch_matrix.numcols()/2;
	
if(braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "rat_poly_invariant: switch_matrix_N_factor = " << switch_matrix_N_factor << endl;

	/* The two cases of the next conditional clause are in-line rather than function calls only because
	   the code appear only once
	*/
//	if (input_string.find('(') != string::npos)
	if (input_string.find('(') != string::npos || input_string.find('[') != string::npos)
	{
		/* record whether this is a virtual knot or a classical one */
		if (input_string.find('*') != string::npos)
			virtual_crossings_present = true;

		/* check for a long knot */
		if (input_string.find('L') != string::npos)
		{
			braid_control::LONG_KNOT = true;

			input_string = parse_long_knot_input_string(input_string);
			
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "rat_poly_invariant: long knot: " << input_string << endl;
	
			if (input_string == "link")
			{				
				if (!braid_control::SILENT_OPERATION)
					cout << "\nError! long knot indicator provided for the peer code of a link" << endl;
			
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Error! long knot indicator provided for the peer code of a link" << endl;
				}
					
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "rat_poly_invariant: Error! long knot indicator provided for the peer code of a link" << endl;
	
				return;
			}
		}
		else
			braid_control::LONG_KNOT = false;
	
		R_module_rep(Matrix_rep, switch_matrix, switch_matrix_inverse, input_string);
	}
	else if (valid_braid_input(input_string, num_terms, num_strings))
	{ 

		/* remove any unwanted qualifiers from the input string */
		string::size_type pos = input_string.find('{');
		if (pos != string::npos)
			input_string = input_string.substr(0,pos);

		/* record whether this is a virtual knot or a classical one */
		if (input_string.find('t') != string::npos || input_string.find('T') != string::npos)
			virtual_crossings_present = true;
			
		/* When using a braid representation, we do not calculate polynomials if the representation is the identity
		   matrix.  This check is done within braid_rep */
		calculate_polys = braid_rep(Matrix_rep, switch_matrix, switch_matrix_inverse,input_string, num_terms, num_strings);
	} 

	/* Under normal circumstances Matrix_rep is created by R_module_rep or braid_rep.  However, if R_module_rep
	   determines that an immersion code contains only virtual crossings it returns a zero Matrix_rep pointer. */
	if (Matrix_rep == 0)
		return;
	
	/* Note the size of the matrix representation, this is either square or in the case of long knots
	   has one more N-column than N-row, so we note the number of N-rows.  For elements of VB(n) this is the 
	   N * number of strings, but for R_module representations it is 2 * N * number of classical crossings, 
	   where N is the switch matrix N-factor.  Here, we also need to  know the number of N-rows in the 
	   matrix representation so we can remove N-rows and N-columns.
	*/
	int N = switch_matrix_N_factor;
	int matrix_rows = Matrix_rep->numrows();
	int matrix_cols = Matrix_rep->numcols();
	int matrix_N_rows = matrix_rows/N;
	int matrix_N_cols = matrix_cols/N;
	
	if (calculate_polys)
	{
		/* First determine the codimension 0 invariant.  In the case of a long knot we have to remove successive columns  
		   and take the determinant of remianing submatrices, producing a set of generators for the first elementary ideal.
		   If we have been told, via braid_control::NUMERATOR_GCD, that we may calculate the gcd of the numerators we do so.
		   
		   In the case of a closed knot the codimension 0 invariant is just the determinant of the matrix representation.
		   
		   We then determine whether to calculate delta_1 based on the codimension zero value and the global booleans
		*/
		bool calc_delta_1;
	
		if (braid_control::WAIT_SWITCH)
		{
		    matrix_control::WAIT_INFO = true;
		    matrix_control::wait_count = 0;

		    /* the wait_threshold comes from the command line */
			if (!braid_control::SILENT_OPERATION)
				cout << "\nCalculating polynomials, please wait\n";
		}
		else
		    matrix_control::WAIT_INFO = false;

		Qpolynomial delta_0 = polynomial<scalar, bigint> ("0"); // required in 'else' of this clause and in the following clause
		
		if (braid_control::LONG_KNOT)
		{	
			/* take out each N-column of the matrix representation in turn and take the determinant
			   of the rest, keeping a record of the gcd of the determinants calculated.
			*/
			int rperm[matrix_rows];
			int cperm[matrix_cols];
			bool unit_detected = false;

			for (int i = 0 ; i< matrix_rows; i++)
		    	rperm[i] = i;

			/* hcf is used only in single variable examples */
			polynomial<scalar,bigint> hcf = polynomial<scalar,bigint>("0");

			Qpolynomial op_0 = polynomial<scalar,bigint>("0"); // \def\op{\buildrel o \over p}
			Qpolynomial np_0 = polynomial<scalar,bigint>("0"); // \def\np{\buildrel n \over p}

			for (int i = 0 ; i< matrix_N_cols; i++)
			{
				/* take out N-column i */
				for (int j=0; j<N*i; j++)
   					cperm[j] = j;
				for (int j=N*i+N; j<matrix_cols; j++)
  						cperm[j-N] = j;

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "rat_poly_invariant: removing N-column " << i << endl;

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "rat_poly_invariant: cperm: ";
   	for (int j = 0 ; j< matrix_cols-N; j++)
		debug << cperm[j] << " ";
	debug << endl;
}
				Qpolynomial det =  determinant(*Matrix_rep, title, matrix_rows, rperm, cperm);			
							
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "p^{(0)} generator " << i << " = " << det;
				}

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "rat_poly_invariant: p_0 generator det " << i << " = " << det << endl;

				if (i == 0)
				{
					op_0 = det;

					if (braid_control::WAIT_SWITCH)
					    matrix_control::wait_count = 0;

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "rat_poly_invariant: op_0  set to det " << i << " = " << op_0 << endl;

				}
				else if (i == matrix_N_cols-1)
				{
					np_0 = det;

					if (braid_control::WAIT_SWITCH)
					    matrix_control::wait_count = 0;

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "rat_poly_invariant: np_0  set to det " << i << " = " << np_0 << endl;
	
					if (unit_detected)
						break;  //we've jumped to here, so we're now done.

				}

				if (!braid_control::SILENT_OPERATION)
					cout << "\np^{(0)} generator " << i << " = " << det;

				polynomial<scalar,bigint> det_numerator = det.getn();
		    	sanitize(&det_numerator);
				polynomial<scalar,bigint> det_denominator = det.getd();
		    	sanitize(&det_denominator);

				if (braid_control::DELTA_1_UNIT_CHECK && det_denominator.isunit() && det_numerator.isunit())
				{
					unit_detected = true; 
					if (i != matrix_N_cols-1)
						i = matrix_N_cols-2; // jump to last column to determine np_0, i will be incremented by loop

					hcf = polynomial<scalar,bigint>("1");

					if (!braid_control::SILENT_OPERATION)
						cout << "\nunit generator detected, terminating calculation" << endl;
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "unit generator detected, terminating calculation" << endl;
					}

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "rat_poly_invariant: unit generator detected, terminating calculation" << endl;
				}
				else if (braid_control::NUMERATOR_GCD)
				{
					/* take hcf code for single variable examples */
					polynomial<scalar,bigint> generator = det.getn();
		    		sanitize(&generator);

		    		hcf = gcd(hcf,generator);
					
					if (!braid_control::SILENT_OPERATION)
						cout << "\n numerator gcd stands at " << hcf << endl;
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
		    			output << " numerator gcd stands at " << hcf << endl;
					}

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "rat_poly_invariant:  numerator gcd stands at " << hcf << endl;

					if (braid_control::DELTA_1_UNIT_CHECK && hcf.isunit())
					{
						if (!braid_control::SILENT_OPERATION)
							cout << "\nunit gcd detected, terminating calculation" << endl;
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "unit gcd detected, terminating calculation" << endl;
						}

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "rat_poly_invariant: unit gcd detected, terminating calculation" << endl;

						unit_detected = true; 
						if (i != matrix_N_cols-1)
							i = matrix_N_cols-2; // jump to last column to determine np_0, i will be incremented by loop
					}
				}
			}

			if (!braid_control::SILENT_OPERATION)
				cout << endl;

			if (braid_control::NUMERATOR_GCD)
			{
				if (!braid_control::SILENT_OPERATION)
					cout << (braid_control::LONG_KNOT? "p^{(0)} = ":"Delta_0 = ") << hcf << endl;	
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "numerator gcd  = ";
					if (braid_control::OUTPUT_AS_INPUT)
						output << "\n";
				}
				output << hcf << endl;
			}

			if (!braid_control::SILENT_OPERATION)
			{
				cout << "\n\\op^{(0)} = " << op_0 << endl;		
				cout << "\n\\np^{(0)} = " << np_0 << endl;		
   			}
   			
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
	   			output << "\\op^{(0)} = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << op_0 << endl;		
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
	   			output << "\\np^{(0)} = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << np_0 << endl;		

if (braid_control::DEBUG >= braid_control::SUMMARY)
{
   	debug << "rat_poly_invariant: p_0  = p^{(0)} = " << delta_0 << endl;
   	debug << "rat_poly_invariant: \\op^{(0)} = " << op_0 << endl;
   	debug << "rat_poly_invariant: \\np^{(0)} = " << np_0 << endl;
}

    		if (braid_control::NUMERATOR_GCD)
			{
				if (hcf == polynomial<scalar,bigint>("0"))
					calc_delta_1 = true;
				else
					calc_delta_1 = false;
			}
			else
				calc_delta_1 = true; // we may have only been able to produce generators for p^{(0)}, so calculate p^{(1)} regardless
			
		}
		else //closed knot
		{
			if (virtual_crossings_present || braid_control::VERIFY_DELTA_0)
			{		
				delta_0 = determinant (*Matrix_rep, title);

				if (!braid_control::SILENT_OPERATION)
					cout << "\nDelta_0 = " << delta_0 << endl;
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Delta_0 = ";
					if (braid_control::OUTPUT_AS_INPUT) 
						output << "\n";
				}
				output << delta_0 << endl;

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "rat_poly_invariant: Delta_0 = " << delta_0 << endl;

				if (delta_0 == Qpolynomial("0"))
					calc_delta_1 = true;
				else
					calc_delta_1 = false;
			}
			else
			{
	    		calc_delta_1 = true;

				if (!braid_control::DISPLAY_DELTA_1_ONLY)
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "\nInput contains no virtual crossings, so Delta_0 = 0" << endl;
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Input contains no virtual crossings, so Delta_0 = ";
						if (braid_control::OUTPUT_AS_INPUT)
							output << "\n";
					}
					output << "0" << endl;

				}
			}
		}

		if (!braid_control::CALCULATE_DELTA_0_ONLY && (braid_control::ALWAYS_CALCULATE_DELTA_1 || calc_delta_1))
		{
			/* hcf is used only in single variable examples */
			polynomial<scalar,bigint> hcf = delta_0.getn();

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "\nrat_poly_invariant: hcf initialized to " << hcf << endl;


			int rperm[matrix_rows];
			int cperm[matrix_cols];

			/* take out a (N) row and (N) column from the matrix representation 
			   and calculate the determinant of what is left
			*/
			for (int Nr1=0;Nr1<matrix_N_rows;Nr1++)
			{
				/* take out the N rows of the matrix rep starting from row Nr1*N row */
				for (int i=0; i<Nr1*N; i++)
	    			rperm[i] = i;
				for (int i = Nr1*N+N; i< matrix_rows; i++)
					rperm[i-N] = i;


				int Nc2; // declared here so in scope for debug output below;
				for (int Nc1= 0; Nc1< matrix_N_cols; Nc1++)
				{
				
					Qpolynomial det = polynomial<scalar, bigint> ("0");
					
					/* take out the N columns of the matrix rep starting from column N*j */
					for (int i=0; i<Nc1*N; i++)
		   				cperm[i] = i;

					if (braid_control::LONG_KNOT)
					{
						/* There is one extra N-column than N-row */
						for (Nc2= Nc1+1; Nc2< matrix_N_cols; Nc2++)
						{
							/* take out Nc2-nd N column*/
							for (int i=Nc1*N+N; i<Nc2*N; i++)
		   						cperm[i-N] = i;

							for (int i=Nc2*N+N; i<matrix_cols; i++)
		   						cperm[i-2*N] = i;
								
							if (rat_minor_determinant(Matrix_rep, N, Nr1, Nc1, Nc2, rperm, cperm, title, hcf))
								goto rat_poly_delta_1_calculation_complete;
						}
					}
					else
					{
						for (int i=Nc1*N+N; i<matrix_cols; i++)
	   						cperm[i-N] = i;
						
						if (rat_minor_determinant(Matrix_rep, N, Nr1, Nc1, Nc2, rperm, cperm, title, hcf))
							goto rat_poly_delta_1_calculation_complete;
					}
				}			
			}

rat_poly_delta_1_calculation_complete:

			if (!braid_control::SILENT_OPERATION)
				cout << endl;
			
			if (braid_control::NUMERATOR_GCD)
			{
				if (!braid_control::SILENT_OPERATION)
					cout << (braid_control::LONG_KNOT? "p^{(1)} = ":"Delta_1 = ") << hcf << endl;	
				
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "numerator gcd  = ";
					if (braid_control::OUTPUT_AS_INPUT)
						output << "\n";
				}
				output << hcf << endl;

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "rat_poly_invariant: " << (braid_control::LONG_KNOT? "p^{(1)} = ":"Delta_1 = ") << hcf << endl;	

			}
		}
	}
			
	delete Matrix_rep;
}


void commutative_automorphism_invariant(const Qpmatrix& phi, const Qpmatrix& psi, string input_string, string title)
{

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "commutative_automorphism_invariant: Calculating commutative automorphism invariant for " << title << " = " << input_string << endl;

	if (title.length())
   	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\n\n" << title << endl;
		if (!braid_control::RAW_OUTPUT)
			output << "\n\n" << title;
   	}		
   		
if (braid_control::DEBUG >= braid_control::SUMMARY)
{
	debug << "commutative_automorphism_invariant: phi = " << endl;
	print(phi,debug,3,"commutative_automorphism_invariant: ");
	debug << "commutative_automorphism_invariant: psi = " << endl;
	print(psi,debug,3,"commutative_automorphism_invariant: ");
}

	cout << "\nCalculating commutative automorphism invariant for " << title << " = " << input_string << endl;
	cout << "phi = " << endl;
	print(phi,cout,3,"commutative_automorphism_invariant: ");
	cout << "psi = " << endl;
	print(psi,cout,3,"commutative_automorphism_invariant: ");
	
	
}

