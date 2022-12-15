/************************************************************************
                  Support functions for braid tools

                  A. Bartholomew 25th November, 2001

string invstr(string& str)
void print_prog_params (ostream& os, string level, string prefix)
string long_knot_concatenation (string ic1, string ic2)
string set_long_knot_infinity_point (string input_string)
string parse_long_knot_input_string (string input_string)
bool XxX_permutation (matrix<int>& U, matrix<int>& D)
bool rows_are_permutations (matrix<int>&M)
bool Yang_Baxter_satisfied (matrix<int>& U, matrix<int>& D, bool rack_condition = false)
bool one_dominates(matrix<int>& Ua, matrix<int>& Da, matrix<int>& Ub, matrix<int>& Db)
bool two_dominates(matrix<int>& Ua, matrix<int>& Da, matrix<int>& Ub, matrix<int>& Db)
bool power_2(matrix<int>& U, matrix<int>& D)
int switch_order(matrix<int>& U, matrix<int>& D)
void display_fixed_point_switch(matrix<int>& M, ostream& os, bool number_from_zero)
void flip_braid(string& braid)
void invert_braid(string& braid)
void plane_reflect_braid(string& braid)
void line_reflect_braid(string& braid)
void Kamada_double_covering(string& braid)
void commutative_automorphism_invariant(const Qpmatrix& phi, const Qpmatrix& psi, string input_string, string title)

**************************************************************************/
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

//bool valid_braid_input (string input_string, int& num_terms, int& num_strings, bool raw_output, bool silent_output, bool output_as_input);
string braid_to_generic_code (string braid, int num_terms, int code_type);
//void parse_braid(string word, int num_terms, vector<int>& braid_num, vector<int>& type);
bool Yang_Baxter_satisfied (matrix<int>& U, matrix<int>& D, bool rack_condition = false);
bool rat_minor_determinant (Qpmatrix* Matrix_rep, int N, int Nr1, int Nc1, int Nc2, int* rperm, int* cperm, string title, polynomial<scalar,char,bigint>& hcf);

void print_prog_params (ostream& os, string level, string prefix)
{
	os << prefix << "braid_control::AFFINE_INDEX = " << braid_control::AFFINE_INDEX << endl;
	os << prefix << "braid_control::ALEXANDER = " << braid_control::ALEXANDER << endl;
	os << prefix << "braid_control::ARROW_POLYNOMIAL = " << braid_control::ARROW_POLYNOMIAL << endl;
	os << prefix << "braid_control::BURAU = " << braid_control::BURAU << endl;
	os << prefix << "braid_control::COMMUTATIVE_AUTOMORPHISM = " << braid_control::COMMUTATIVE_AUTOMORPHISM << endl;
	os << prefix << "braid_control::DOODLE_ALEXANDER = " << braid_control::DOODLE_ALEXANDER << endl;
	os << prefix << "braid_control::DOWKER_CODE = " << braid_control::DOWKER_CODE << endl;
	os << prefix << "braid_control::DYNNIKOV_TEST = " << braid_control::DYNNIKOV_TEST << endl;
	os << prefix << "braid_control::FIXED_POINT_INVARIANT = " << braid_control::FIXED_POINT_INVARIANT << endl;
	os << prefix << "braid_control::GAUSS_CODE = " << braid_control::GAUSS_CODE << endl;
	os << prefix << "braid_control::HOMFLY = " << braid_control::HOMFLY << endl;
	os << prefix << "braid_control::IMMERSION_CODE = " << braid_control::IMMERSION_CODE << endl;
	os << prefix << "braid_control::JONES_POLYNOMIAL = " << braid_control::JONES_POLYNOMIAL << endl;
	os << prefix << "braid_control::KAUFFMAN_BRACKET = " << braid_control::KAUFFMAN_BRACKET << endl;
	os << prefix << "braid_control::KNOTOID_BRACKET = " << braid_control::KNOTOID_BRACKET << endl;
	os << prefix << "braid_control::MATRIX = " << braid_control::MATRIX << endl;
	os << prefix << "braid_control::MOCK_ALEXANDER = " << braid_control::MOCK_ALEXANDER << endl;
	os << prefix << "braid_control::PARITY_ARROW = " << braid_control::PARITY_ARROW << endl;
	os << prefix << "braid_control::PARITY_BRACKET = " << braid_control::PARITY_BRACKET << endl;
	os << prefix << "braid_control::PEER_CODE = " << braid_control::PEER_CODE << endl;
	os << prefix << "braid_control::QUATERNION = " << braid_control::QUATERNION << endl;
	os << prefix << "braid_control::SAWOLLEK = " << braid_control::SAWOLLEK << endl;
	os << prefix << "braid_control::SATELLITE = " << braid_control::SATELLITE << endl;
	os << prefix << "braid_control::SWITCH_POLYNOMIAL_INVARIANT = " << braid_control::SWITCH_POLYNOMIAL_INVARIANT << endl;
	os << prefix << "braid_control::VOGEL_ALGORITHM = " << braid_control::VOGEL_ALGORITHM << endl;
	os << prefix << "braid_control::WELDED_BRAID = " << braid_control::WELDED_BRAID << endl;
	os << prefix << "braid_control::WEYL = " << braid_control::WEYL << endl;

	if (level == "detail" )
	{
		os << prefix << "braid_control::ALWAYS_CALCULATE_CODIMENSION_2_DELTA_1 = " << braid_control::ALWAYS_CALCULATE_CODIMENSION_2_DELTA_1 << endl;
		os << prefix << "braid_control::ALWAYS_CALCULATE_DELTA_1 = " << braid_control::ALWAYS_CALCULATE_DELTA_1 << endl;
		os << prefix << "braid_control::BIGELOW_KNOT_SEARCH = " << braid_control::BIGELOW_KNOT_SEARCH << endl;
		os << prefix << "braid_control::BYPASS_FUNDAMENTAL_EQUATION_CHECK = " << braid_control::BYPASS_FUNDAMENTAL_EQUATION_CHECK << endl;
		os << prefix << "braid_control::CALCULATE_DELTA_0_ONLY = " << braid_control::CALCULATE_DELTA_0_ONLY << endl;
		os << prefix << "braid_control::CALCULATE_MOD_P = " << braid_control::CALCULATE_MOD_P << endl;
		os << prefix << "braid_control::COMPLEX_STUDY_DELTA_1 = " << braid_control::COMPLEX_STUDY_DELTA_1 << endl;
		os << prefix << "braid_control::CUSTOM_WEYL = " << braid_control::CUSTOM_WEYL << endl;
		os << prefix << "braid_control::DELTA_1_UNIT_CHECK = " << braid_control::DELTA_1_UNIT_CHECK << endl;
		os << prefix << "braid_control::DISPLAY_DELTA_1_ONLY = " << braid_control::DISPLAY_DELTA_1_ONLY << endl;
		os << prefix << "braid_control::EQUALITY_TEST = " << braid_control::EQUALITY_TEST << endl;
		os << prefix << "braid_control::EXPANDED_BRACKET_POLYNOMIAL = " << braid_control::EXPANDED_BRACKET_POLYNOMIAL << endl;
		os << prefix << "braid_control::EXTRA_OUTPUT = " << braid_control::EXTRA_OUTPUT << endl;
		os << prefix << "braid_control::FLAT_VOGEL_MOVES = " << braid_control::FLAT_VOGEL_MOVES << endl;
		os << prefix << "braid_control::FLIP_BRAID = " << braid_control::FLIP_BRAID << endl;
		os << prefix << "braid_control::GCD_BACK_SUBSTITUTING = " << braid_control::GCD_BACK_SUBSTITUTING << endl;
		os << prefix << "braid_control::INVERT_BRAID = " << braid_control::INVERT_BRAID << endl;
		os << prefix << "braid_control::KAMADA_DOUBLE_COVERING = " << braid_control::KAMADA_DOUBLE_COVERING << endl;
		os << prefix << "braid_control::LINE_REFLECT_BRAID = " << braid_control::LINE_REFLECT_BRAID << endl;
		os << prefix << "braid_control::LONG_KNOT = " << braid_control::LONG_KNOT << endl;
		os << prefix << "braid_control::NORMALIZE_BRACKET = " << braid_control::NORMALIZE_BRACKET << endl;
		os << prefix << "braid_control::NORMALIZING_Q_POLYNOMIALS = " << braid_control::NORMALIZING_Q_POLYNOMIALS << endl;
		os << prefix << "braid_control::NUMERATOR_GCD = " << braid_control::NUMERATOR_GCD << endl;
		os << prefix << "braid_control::OUTPUT_AS_INPUT = " << braid_control::OUTPUT_AS_INPUT << endl;
		os << prefix << "braid_control::PLANE_REFLECT_BRAID = " << braid_control::PLANE_REFLECT_BRAID << endl;
		os << prefix << "braid_control::PRIME_WEYL = " << braid_control::PRIME_WEYL << endl;
		os << prefix << "braid_control::QUANTUM_WEYL = " << braid_control::QUANTUM_WEYL << endl;
		os << prefix << "braid_control::RACK_POLYNOMIAL = " << braid_control::RACK_POLYNOMIAL << endl;
		os << prefix << "braid_control::RAW_OUTPUT = " << braid_control::RAW_OUTPUT << endl;
		os << prefix << "braid_control::RELAXED_PARITY = " << braid_control::RELAXED_PARITY << endl;
		os << prefix << "braid_control::REMOVE_PEER_CODE_COMPONENT = " << braid_control::REMOVE_PEER_CODE_COMPONENT << endl;
		os << prefix << "braid_control::REMOVE_REIDEMEISTER_II_MOVES = " << braid_control::REMOVE_REIDEMEISTER_II_MOVES << endl;
		os << prefix << "braid_control::SILENT_OPERATION = " << braid_control::SILENT_OPERATION << endl;
		os << prefix << "braid_control::STATUS_INFORMATION = " << braid_control::STATUS_INFORMATION << endl;
		os << prefix << "braid_control::STUDY_RHO_MAPPING = " << braid_control::STUDY_RHO_MAPPING << endl;
		os << prefix << "braid_control::SWITCH_POWER = " << braid_control::SWITCH_POWER << endl;
		os << prefix << "braid_control::T_VARIABLE = " << braid_control::T_VARIABLE << endl;
		os << prefix << "braid_control::TMP_DIRECTORY = " << braid_control::TMP_DIRECTORY << endl;
		os << prefix << "braid_control::TRUNCATED_WEYL = " << braid_control::TRUNCATED_WEYL << endl;
		os << prefix << "braid_control::VERIFY_DELTA_0 = " << braid_control::VERIFY_DELTA_0 << endl;
		os << prefix << "braid_control::WAIT_SWITCH = " << braid_control::WAIT_SWITCH << endl;
	}
}

/* long_knot_concatenation yields a peer_code string representation of the concatination of the long knots
   determined by the input codes code1 and code2. The input codes can be either peer codes or immersion codes.

   To concatenate the long knots we create a combined generic_code_data structure and write it to a stringstream.
   Since writing only requires the head, num_components, OPEER, TYPE, LABEL and COMPONENT data, we calculate only
   these.  Note however, that since we are only dealing with knots, there is only one component to consider.
   
   We create a combined code table by juxtaposing code_table_1 and code_table_2 and then adjusting the OPEER row
   by incrementing those corresponding to code_data_2 by the number of edges corresponding to code_data_1.  
     
   The type and labels are unchanged by the concatenation and the COMPONENT row will correctly be populated with zeros.
   
*/
string long_knot_concatenation (string code1, string code2)
{

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "long_knot_concatenation: presented with strings :" << endl;
	debug << "long_knot_concatenation:   code1 = " << code1 << endl;
	debug << "long_knot_concatenation:   code2 = " << code2 << endl;
}
	generic_code_data code_data_1;
	generic_code_data code_data_2;

	read_code_data(code_data_1, code1);
	matrix<int>& code_table_1 = code_data_1.code_table;
	int n1 = code_data_1.num_crossings;

/*
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "long_knot_concatenation: code data produced from input code code1:" << endl;
	print (code_data_1,debug,"long_knot_concatenation: ");
}
*/
	read_code_data(code_data_2, code2);
	matrix<int>& code_table_2 = code_data_2.code_table;
	int n2 = code_data_2.num_crossings;
/*	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "long_knot_concatenation: code table produced from input code code2:" << endl;
	print (code_data_2,debug,"long_knot_concatenation: "); 
}
*/
	matrix<int> code_table(CODE_TABLE_SIZE,n1 + n2);
	
	for (int i = 0; i<CODE_TABLE_SIZE; i++)
	for (int j = 0; j<n1; j++)
		code_table[i][j] = code_table_1[i][j];
		
	for (int i = 0; i<CODE_TABLE_SIZE; i++)
	for (int j = 0; j<n2; j++)
		code_table[i][n1+j] = code_table_2[i][j];
		
	int num_edges_1 = 2* code_data_1.num_crossings;
		
	for (int i = 0 ; i<n2; i++)
		code_table[OPEER][n1+i] += num_edges_1;
	
	/* write the new peer code to a stringstream */
	ostringstream oss;
	code_data_1.type = generic_code_data::peer_code;
//	code_data_1.head = -1;
	code_data_1.num_crossings += code_data_2.num_crossings;
//	code_data_1.num_components += code_data_2.num_components - 1 ;
	code_data_1.code_table = code_table;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "long_knot_concatenation: concatenated code table:" << endl;
	print_code_data (debug,code_data_1,"long_knot_concatenation: ");
}

	write_peer_code (oss, code_data_1);
	return oss.str();
}


string set_long_knot_infinity_point (string input_string)
{
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "set_long_knot_infinity_point: presented with input_string: " << input_string << endl;

	generic_code_data code_data;
	read_code_data(code_data, input_string);

	if (code_data.num_components > 1)
	{
		return "link";
	}
	
	int shift = 0;
	
	char* inbuf = c_string(input_string);
	char* cptr = strchr(inbuf, 'L'); // We know there's an L in the string, otherwise this function wouldn't be called.
	cptr++; // step over L
	cptr++; // step over :
	if (*cptr == '+' || *cptr == '-')
	{
		/* a shift of infinity has been requested for the long knot */
		get_number(shift, cptr);
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "set_long_knot_infinity_point: long knot shift = " << shift << endl;

		/* move over number so it can be removed from the input string */
		cptr++;
		while (isdigit(*cptr))
			cptr++;
	}
			
//	input_string.erase(0,cptr-inbuf);

//if (debug_control::DEBUG >= debug_control::DETAIL)
//	debug << "set_long_knot_infinity_point: input string after L and shift removal " << input_string << endl;

	if (code_data.type == generic_code_data::peer_code && shift)
	{			
		/* The function renumber_peer_code moves the point at infinity forwards with respect
		   to the orientation, so if the requested shift is negtive, we have to adjust it accordingly.
		*/
		int num_crossings = code_data.num_crossings;

		// we do this repeatedly in case we've been given a huge negative shift.
		while (shift < 0)
			shift += 2*num_crossings;
		
		shift %= 2*num_crossings;
		
		vector<int> shift_vector(1); 
		shift_vector[0]=abs(shift);
		renumber_peer_code (code_data,shift_vector);

		/* write the new immersion code to a stringstream */
		ostringstream oss;
//		code_data.head = -1;
//		code_data.num_crossings = num_crossings;
//		code_data.code_table = code_table;

		write_code_data(oss,code_data);

		input_string = oss.str();
	}


if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "set_long_knot_infinity_point: returning immersion code " << input_string << endl;

	delete[] inbuf;

	return input_string;
}

string parse_long_knot_input_string (string input_string)
{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "parse_long_knot_input_string: presented with input string " << input_string << endl;

	/* Check to see if this is a concatenation of multiple long knots */
	string::size_type tilde = input_string.find('~');

	string result = set_long_knot_infinity_point(input_string.substr(0,tilde));

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "parse_long_knot_input_string: initial component: " << result << endl;

	while( result != "link" && tilde != string::npos) // result != "link" catches the case where the first component is a link
	{
		string::size_type next_tilde = input_string.find('~',tilde+1);		
		string component = set_long_knot_infinity_point(input_string.substr(tilde+1,next_tilde-tilde-1));

		if (component =="link")
		{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "parse_long_knot_input_string: next component: " << component << endl;

			result = component;
			break;
		}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "parse_long_knot_input_string: next component: " << component << endl;

		result = long_knot_concatenation (result, component);
		tilde = next_tilde;
	}
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "parse_long_knot_input_string: input string after parsing: " << result << endl;

	return result;

}

bool XxX_permutation (matrix<int>& U, matrix<int>& D)
{
	int n = U.numrows();
	
	for (int k=0; k< n; k++)
	{
		/* for each U[i][j]=k, D[j][i] must be a permutation */
		valarray<int> perm(0,n);
		for (int i=0; i< n; i++)
		for (int j=0; j< n; j++)
		{
			if (U[i][j] == k)
			{
				if (perm[D[j][i]])
					return false; // not a permutation
				else
					perm[D[j][i]] = 1;
			}
		}
	}				
	return true;
}

bool rows_are_permutations (matrix<int>&M)
{
	int n = M.numrows();
	
	for (int i=0; i<n; i++)
	{
		valarray<int> check_perm(0,n);
		for (int j=0; j< n; j++)
		{
			if (check_perm[M[i][j]])
				return false; // not a permutation
			else
				check_perm[M[i][j]] = 1;
		}				
	}
	return true;
}


/* If the boolean rack_condition is false (the default) the Yang_Baxter_satisfied function checks 
   the up down and rule of five conditions for a birack.  If rack_condition is true, it only
   tests the up condition for a rack.
*/

bool Yang_Baxter_satisfied (matrix<int>& U, matrix<int>& D,bool rack_condition)
{
	int n = U.numcols();

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "  S-Yang-Baxter rules: " << endl;

	vector<int> rof(3); //rule of five
	for (int i =0; i< n; i++)
	{
		rof[0]=i;
		for (int j =0; j< n; j++)
		{
			rof[1]=j;
			for (int k =0; k < n; k++)
			{
				rof[2]=k;
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "    rof: " << rof[0] << " " << rof[1] << " " << rof[2] << endl;

				int& a = rof[0];
				int& b = rof[1];
				int& c = rof[2];
					
				/* Up interchange */
				int& a_up_b = U[b][a];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      a_up_b = " << a_up_b <<endl;

				int& b_up_c = U[c][b];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      b_up_c = " << b_up_c <<endl;

				int& a_up_b__up_c = U[c][a_up_b];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      a_up_b__up_c = " << a_up_b__up_c <<endl;

				if (rack_condition)
				{
					/* a^bc = a^{cb^c} */
					int& a_up_c = U[c][a];						
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      a_up_c = " << a_up_c <<endl;

					if (a_up_b__up_c != U[b_up_c][a_up_c])
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      rack up interchange condition fails" << endl;
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "    fails" << endl;
	
						return false;
					}
					else
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      a_up_b__up_c == U[b_up_c][a_up_c]\n" << endl;
					}
				}
				else 
				{
					int& c_down_b = D[b][c];					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      c_down_b = " << c_down_b <<endl;

					int& a_up__c_down_b = U[c_down_b][a];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      a_up__c_down_b = " << a_up__c_down_b <<endl;
					
					if (a_up_b__up_c != U[b_up_c][a_up__c_down_b])
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      up interchange condition fails" << endl;
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "    fails" << endl;
	
						return false;
					}
					else
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      a_up_b__up_c == U[b_up_c][a_up__c_down_b]\n" << endl;
					}
				}

				if (!rack_condition)
				{
					/* Down interchange */
					int& c_up_b = U[b][c];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      c_up_b = " << c_up_b <<endl;

					int& b_down_c = D[c][b];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      b_down_c = " << b_down_c <<endl;

					int& a_down_b = D[b][a];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      a_down_b = " << a_down_b <<endl;

					int& a_down_b__down_c = D[c][a_down_b];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      a_down_b__down_c = " << a_down_b__down_c <<endl;

					int& a_down__c_up_b = D[c_up_b][a];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      a_down__c_up_b = " << a_down__c_up_b <<endl;
					
					if (a_down_b__down_c != D[b_down_c][a_down__c_up_b])
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      down interchange condition fails" << endl;
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "    fails" << endl;
	
						return false;
					}
					else
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      a_down_b__down_c == D[b_down_c][a_down__c_up_b]\n" << endl;
					}
					
					/* Rule of five */
					int& b_up_a = U[a][b];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      b_up_a = " << b_up_a <<endl;
	
					int& a_up_c = U[c][a];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      a_up_c = " << a_up_c <<endl;
	
					int& c_down_a = D[a][c];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      c_down_a = " << c_down_a <<endl;
	
					int& b_up__c_down_a = U[c_down_a][b];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      b_up__c_down_a = " << b_up__c_down_a <<endl;
		
					int& c_down__b_up_a = D[b_up_a][c];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      c_down__b_up_a = " << c_down__b_up_a <<endl;
	
					int& a_down_b__up___c_down__b_up_a = U[c_down__b_up_a][a_down_b];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      a_down_b__up___c_down__b_up_a = " << a_down_b__up___c_down__b_up_a <<endl;
	
					
					if (a_down_b__up___c_down__b_up_a != D[b_up__c_down_a][a_up_c])
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "      rule of five condition fails" << endl;
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "  fails" << endl;
	
						return false;
					}	
				}
			}
		}
	}
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "    OK" << endl;

	return true;
}


/* one_dominates tests whether the switch A, determined by up and down matrices Ua and Da
   (a candidate for S), 1-dominates the switch B, determined by up and down matrices Ub and Db
   (a candidate for T).
   
   The swith A 1-dominates B if 
   
   (A x I) (I x A) (B x 1) = (I x B) (A x I) (I x A)   (*)  
   
   Each case is checked by evaluating the matrix representing the composition on
   the lhs and rhs.

   Writing A(x,y) for the effect of the switch A, we have
   (A x I)(i1,i2,i3)=(A(i1,i2),i3) and (I x A)(i1,i2,i3)=(i1,A(i2,i3)), etc..
   
   Note that we are associating the labels i1, i2 and i3 to the braid strands in a top down 
   manner here, consistent with our convention of writing S(x,y) = (y^x,x_y).
   
   If we write S1 = Id x S, and S2 = S x Id, and similarly for T, this function checks the
   relationship
   
		S2S1T2 = T1S2S1
		
   which corresponds to the conditions imposed on S and T by the allowable weld move in a braid.   
   	   
   The effect of a switch S with up and down matrices U and D on the pair (x,y) is given by
   S(x,y) = (U[x][y],D[y][x]), since the (x,y) entry of U and D are y^x and y_x respectively
   and S(x,y) = (y^x,x_y).

   The function compares the lhs and rhs compositions in (*), the composition for the lhs gives 
   (i1,i2,i3) -> (l1,l2,l3), where (i1,i2,i3)->(j1,j2,j3) under A1, (j1,j2,j3)->(k1,k2,k3) 
   under A2 and (k1,k2,k3)->(l1,l2,l3) under B1.  Similarly for the rhs.
   
   We consider each triple (i1,i2,i3) and evaluate the corresponding (l1,l2,l3) for the lhs
   and (r1,r2,r3) for the rhs and then compare them, if they differ, then the test fails.
*/
bool one_dominates(matrix<int>& Ua, matrix<int>& Da, matrix<int>& Ub, matrix<int>& Db)
{
	int n = Ua.numcols();

	for (int i=0; i < n; i++)
	for (int j=0; j < n; j++)
	for (int k=0; k < n; k++)
	{
		int l1=i,l2=j,l3=k,temp;
		
		/* evaluate the lhs */
		
		/* (aptr x I)(x,y,z) = (Ua[x][y],Da[y][x],z) */
		temp = Ua[l1][l2];
		l2 = Da[l2][l1];
		l1 = temp;

		/* (I x aptr)(x,y,z) = (x, Ua[y][z],Da[z][y]) */
		temp = Ua[l2][l3];
		l3 = Da[l3][l2];
		l2 = temp;
		
		/* (bptr x I)(x,y,z) = (Ub[x][y],Db[y][x],z) */
		temp = Ub[l1][l2];
		l2 = Db[l2][l1];
		l1 = temp;

		int r1=i,r2=j,r3=k;

		/* evaluate the rhs */

		/* (I x bptr)(x,y,z) = (x, Ub[y][z],Db[z][y]) */
		temp = Ub[r2][r3];
		r3 = Db[r3][r2];
		r2 = temp;

		/* (aptr x I)(x,y,z) = (Ua[x][y],Da[y][x],z) */
		temp = Ua[r1][r2];
		r2 = Da[r2][r1];
		r1 = temp;
		
		/* (I x aptr)(x,y,z) = (x, Ua[y][z],Da[z][y]) */
		temp = Ua[r2][r3];
		r3 = Da[r3][r2];
		r2 = temp;
		
		if (r1 != l1 || r2 != l2 || r3 != l3)
			return false;
	}
	
	return true;
}

/* two_dominates tests whether the switch A, determined by up and down matrices Ua and Da
   (a candidate for S), 2-dominates the switch B, determined by up and down matrices Ub and Db
   (a candidate for T).
   
   The swith A 2-dominates B if 
   
   (I x A) (A x I) (I x B) = (B x I) (I x A) (A x I)   (*)  
    
   Each case is checked by evaluating the matrix representing the composition on
   the lhs and rhs.

   Writing A(x,y) for the effect of the switch A, we have
   (A x I)(i1,i2,i3)=(A(i1,i2),i3) and (I x A)(i1,i2,i3)=(i1,A(i2,i3)), etc..
   
   Note that we are associating the labels i1, i2 and i3 to the braid strands in a top down 
   manner here, consistent with our convention of writing S(x,y) = (y^x,x_y).
   
   If we write S1 = Id x S, and S2 = S x Id, and similarly for T, this function checks the
   relationship
   
		S1S2T1 = T2S1S2
		
   which corresponds to conditions imposed on S and T by the forbidden weld move in a braid.

   The effect of a switch S with up and down matrices U and D on the pair (x,y) is given by
   S(x,y) = (U[x][y],D[y][x]), since the (x,y) entry of U and D are y^x and y_x respectively
   and S(x,y) = (y^x,x_y).

   The function compares the lhs and rhs compositions in (*), the composition for the lhs gives 
   (i1,i2,i3) -> (l1,l2,l3), where (i1,i2,i3)->(j1,j2,j3) under A2, (j1,j2,j3)->(k1,k2,k3) 
   under A1 and (k1,k2,k3)->(l1,l2,l3) under B2.  Similarly for the rhs.
   
   We consider each triple (i1,i2,i3) and evaluate the corresponding (l1,l2,l3) for the lhs
   and (r1,r2,r3) for the rhs and then compare them, if they differ, then the test fails.
*/
bool two_dominates(matrix<int>& Ua, matrix<int>& Da, matrix<int>& Ub, matrix<int>& Db)
{
	int n = Ua.numcols();

	for (int i=0; i < n; i++)
	for (int j=0; j < n; j++)
	for (int k=0; k < n; k++)
	{
		int l1=i,l2=j,l3=k,temp;
		
		/* evaluate the lhs */
		
		/* (I x aptr)(x,y,z) = (x, Ua[y][z],Da[z][y]) */
		temp = Ua[l2][l3];
		l3 = Da[l3][l2];
		l2 = temp;

		/* (aptr x I)(x,y,z) = (Ua[x][y],Da[y][x],z) */
		temp = Ua[l1][l2];
		l2 = Da[l2][l1];
		l1 = temp;

		/* (I x bptr)(x,y,z) = (x, Ub[y][z],Db[z][y]) */
		temp = Ub[l2][l3];
		l3 = Db[l3][l2];
		l2 = temp;

		int r1=i,r2=j,r3=k;

		/* evaluate the rhs */
		
		/* (bptr x I)(x,y,z) = (Ub[x][y],Db[y][x],z) */
		temp = Ub[r1][r2];
		r2 = Db[r2][r1];
		r1 = temp;
		
		/* (I x aptr)(x,y,z) = (x, Ua[y][z],Da[z][y]) */
		temp = Ua[r2][r3];
		r3 = Da[r3][r2];
		r2 = temp;
		
		/* (aptr x I)(x,y,z) = (Ua[x][y],Da[y][x],z) */
		temp = Ua[r1][r2];
		r2 = Da[r2][r1];
		r1 = temp;

		if (r1 != l1 || r2 != l2 || r3 != l3)
			return false;
	}
	
	return true;
}

/* power_2 returns true if the switch indicated by the parameter has the property S^2=1
 
   The effect of a switch S with up and down matrices U and D on the pair (x,y) is given by
   S(x,y) = (U[x][y],D[y][x]), since the (x,y) entry of D and U are y_x and y^x respectively
   and S(x,y) = (y^x,x_y).
   
   Therefore power_2 looks at all combinations of (x,y) and checks that S(y^x,x_y) = (x,y)
*/ 
bool power_2(matrix<int>& U, matrix<int>& D)
{
	int n = U.numrows();
	
	for (int i=0; i < n; i++)
	for (int j=0; j < n; j++)
	{		
		if (U[ U[i][j] ][ D[j][i] ] != i)
			return false;

		if (D[ D[j][i] ][ U[i][j] ] != j)
			return false;
	}
	return true;
}

/* switch_order calculates the least k such that S^k(x,y) = (x,y).
   We shall calculate up and down matrices U^k and D^k for S^k so that
   S^k(x,y) = (Uk[x][y],Dk[y][x]) and will have
   S^k(x,y)=(x,y) if the (x,y) entry of U^k is x and the (y,x) entry of D^k is y,
   i.e the row i of U^k and D^k are both constant with entries equal to i.
   
   We know that the order of S divides n^2 since S is a permutation of XxX
*/
int switch_order(matrix<int>& U, matrix<int>& D)
{
	matrix<int> Uk(U),Dk(D);
	int n = U.numrows();		
	
	for (int k=1; k <= n*n; k++)
	{

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	bool saved_bool = matrix_control::SINGLE_LINE_OUTPUT;
	matrix_control::SINGLE_LINE_OUTPUT = true;
	debug << "braid::switch_order: k = " << k << " Uk = " << Uk << " Dk = " << Dk << endl;
	matrix_control::SINGLE_LINE_OUTPUT = saved_bool;
}

		/* test Uk and Dk for the identity */
		bool identity = true;
		
		for (int i=0; i<n && identity==true; i++)
		for (int j=0; j<n; j++)
		{
			if (Uk[i][j] != i)
			{
				identity = false;
				break;
			}
		}
		
		if (identity)
		{
			for (int i=0; i<n && identity==true; i++)
			for (int j=0; j<n; j++)
			{
				if (Dk[i][j] != i)
				{
					identity = false;
					break;
				}
			}
		}
		
		if (identity)
		{

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::switch_order: returning " << k << endl;

			return k;
		}
		else
		{
			/* calculate Uk and Dk for the next value of k.  We iterate using
			   S(x,y)=(y^x,x_y) = (U[x][y], D[y][x]) by writing
			   S^k(x,y) = (Uk[x][y],Dk[y][x]) and then we have
			   S^{k+1}(x,y) = S(S^k(x,y)) 
			                = S(Uk[x][y],Dk[y][x])
							= ((Dk[y][x])^(Uk[x][y]),(Uk[x][y])_(Dk[y][x]))
							= ( U[Uk[x][y]][Dk[y][x]], D[Dk[y][x]][Uk[x][y]] ) 
									   
			   Thus, 			   
			   
			   U{k+1}[x][y] = U[Uk[x][y]][Dk[y][x]]
			   D{k+1}[y][x] = D[Dk[y][x]][Uk[x][y]]
			   
			*/
			matrix<int> locU(n,n),locD(n,n);

			for (int i=0; i<n; i++)
			for (int j=0; j<n; j++)
			{
				locU[i][j] = U[Uk[i][j]][Dk[j][i]];
				locD[j][i] = D[Dk[j][i]][Uk[i][j]];
			}
			
			Uk = locU;
			Dk = locD;
		}
	}
	
	return -1; // error condition, we shouldn't ever get here
}

/* display_fixed_point_switch is used for reporting fixed-point switches since the internal
   representation numbers from zero and we want to report in the form
   presented by the user.
*/
void display_fixed_point_switch(matrix<int>& M, ostream& os, bool number_from_zero)
{
	int n = M.numcols();
	
	for (int i=0; i< n; i++)
	{
		for (int j=0; j< n; j++)
		{
			if (number_from_zero)			
				os << M[i][j];
			else
				os << M[i][j]+1;
			
			if (j < n-1)
				os << " ";
		}
			
		if (i < n-1)
			os << "  ";
	}
}

/* The Kamada double covering of a diagram is defined in "Colourings and doubled colourings of virtual doodles" for flat 
   virtuals and doodles, here implemented for braids.  Consider a braid reflected below within the plane, so that the 
   original strands 1,...,n and the reflected strands 1',...,n' are aligned as in the left side of the following left hand 
   diagram.  On the right side of the left hand diagram are the corresponding strand numbers of the double covering of the 
   braid.  Note that a crossing s_i in the original braid corresponds to the two crossings s_{n+i} and s_{n-i} in the 
   double cover.
   
                                                                                          s_i  
                   n         2n                   i+1 ---  --- n+i+1               i --    *    -- n+i
                                                         \/                            \       /
                   i         n+i                         /\                      i-1 ---o-----o--- n+i-1    
                                                  i --*--  --*-- n+i                     o   o
                   1         n+1                                                   1 -----o-o----- n+1   
                 -----------------                   B_i    B_i                            o
                   1'        n                                                     1 -----o-o----- n
                                                  i --*--  --*-- n-i+1                   o   o
                   i'        n-i+1                       \/                      i-1 ---o-----o--- n-i+2
                                                         /\                            /       \
                   n'        1                    i+1 ---  --- n-i                 i --    *    -- n-i+1
       
   
   The middle diagram shows the replacement of the cut system (the four points indicated by asterisks) around a flat crossing
   by two instances of the braid B_i.  That is, the replacement of crossing s_i in the Kamada double covering is 
   B_i s_{n+i} s_{n-i} B_i.

   The braid B_i is shown in the right diagram: we add the virtual crossings in the following order
    - top left quadrant, top down
    - bottom left quadrant, bottom up
    - centre crossing (t_n)
    - top right quadrant, bottom up
    - bottom right quadrant, top down
    
   The replacement of a virtual crossing t_i in the original braid is given by t_{n+i} t{n-i}.
*/
void Kamada_double_covering(string& braid)
{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::Kamada_double_covering: evaluating the double covering of braid " << braid << endl;

	int num_terms;
	int num_strings;

	if (valid_braid_input(braid, num_terms, num_strings, braid_control::SILENT_OPERATION, braid_control::RAW_OUTPUT, braid_control::OUTPUT_AS_INPUT))
	{
		vector<int> braid_num(num_terms);
	    vector<int> type(num_terms);

		parse_braid(braid, num_terms, braid_num, type);

		ostringstream oss;
		
		for (int i=0; i< num_terms; i++)
		{
			if (type[i] == braid_crossing_type::VIRTUAL)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "braid::Kamada_double_covering:   crossing " << i << " is virtual, replacing t" << braid_num[i] << " with "
	      << "t" << num_strings + braid_num[i] << " t" << num_strings - braid_num[i] << endl;
}

				oss << "t" << num_strings + braid_num[i] << "t" << num_strings - braid_num[i];
			}
			else 
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "braid::Kamada_double_covering:   crossing " << i << " is flat, replacing s" << braid_num[i] << " with "
	      << "B_i s" << num_strings + braid_num[i] << " s" << num_strings - braid_num[i] << " B_i" << endl;
}

				ostringstream B_i;
				
				for (int j=braid_num[i]-1; j >=1; j--)
				    B_i << "t" << num_strings + j;
				    
				for (int j=braid_num[i]-1; j >=1; j--)
				    B_i << "t" << num_strings - j; 

				B_i << "t" << num_strings;

				for (int j=1; j <=braid_num[i]-1; j++)
				    B_i << "t" << num_strings + j;
				    
				for (int j=1; j <=braid_num[i]-1; j++)
				    B_i << "t" << num_strings - j; 


if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "braid::Kamada_double_covering:     where B_i = " << B_i.str() << endl;
}				
				oss << B_i.str() << "s" << num_strings + braid_num[i] << "s" << num_strings - braid_num[i] << B_i.str();
			}
		}
		
		braid = oss.str();
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::Kamada_double_covering: to give braid " << braid << endl;

	}
	else
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::Kamada_double_covering: not a valid braid, doing nothing" << endl;
	}
}

/*  The commutative automorphism invariant determines the codimension 0 invariant (Delta_0) and, where necessary, the codimension 1
    invariant (Delta_1) of the R-module representation of the given labelled peer code determined by the F-biquandle defined from
    a pair of commuatative automorphism, \phi and \psi, of an abelian group.  That is F(b,a) = (a_b, b^a), where 
    a^b = \phi(a)-\phi(b)+\psi(b) and a_b = \psi(a).
    
    Note that although there are many similarities to R_module_rep in the code here, R_module_rep is using equations determined 
    by a switch and here we use equations determined by the sideways map of a switch.
*/
void commutative_automorphism_invariant(const Qpmatrix& phi, const Qpmatrix& psi, string input_string, string title)
{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "commutative_automorphism_invariant: Calculating commutative automorphism invariant for " << title << " = " << input_string << endl;

	if (title.length())
   	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\n\n" << title << "\n" << input_string << endl;
		if (!braid_control::RAW_OUTPUT)
			output << "\n\n" << title << "\n" << input_string << endl;
   	}		
   		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "commutative_automorphism_invariant: phi = " << endl;
	print(phi,debug,3,"commutative_automorphism_invariant: ");
	debug << "commutative_automorphism_invariant: psi = " << endl;
	print(psi,debug,3,"commutative_automorphism_invariant: ");
}

	generic_code_data code_data;
	
	if (input_string.find('(') != string::npos || input_string.find('[') != string::npos)
	{
		read_code_data (code_data,input_string);
	}
	else
	{
		int num_terms;
		int num_strings;
		if (valid_braid_input(input_string, num_terms, num_strings, braid_control::SILENT_OPERATION, braid_control::RAW_OUTPUT, braid_control::OUTPUT_AS_INPUT))
		{
			input_string = braid_to_generic_code(input_string, num_terms, generic_code_data::peer_code);
			read_peer_code(code_data, input_string);
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "commutative_automorphism_invariant: valid braid input converted to peer code " << input_string << endl;

		}
		else
		{
			return; // invalid braid.
		}
	}
	
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
		return;
	}

	/* Determine the size of the R_module representation.  In the case of a knot or link, this is a 
	   2*N*num_classical square matrix where the switch_matrix elements are NxN matrices (for the 
	   Alexander and Quaternionic reps, N=1).  In the case of a long knot, the R_module representation
	   is a 2*N*num_classical x (2*N*num_classical + N) matrix, allowing for the extra variable.  
	*/
	int N = phi.numrows(); // N is used later, it is what is referred to elsewhere as the switch_matrix_N_factor
	int matrix_N_size = 2*num_classical;  // number of N-rows
	int matrix_size = N*matrix_N_size;  // number of actual rows

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "commutative_automorphism_invariant: num_components = " << num_components << endl;
	debug << "commutative_automorphism_invariant: num_crossings = " << num_crossings << endl;
	debug << "commutative_automorphism_invariant: num_classical = " << num_classical << endl;
	debug << "commutative_automorphism_invariant: N (switch_matrix_N_factor) = " << N << endl;
	debug << "commutative_automorphism_invariant: matrix_N_size (number of N-rows) = " << matrix_N_size << endl;
	debug << "commutative_automorphism_invariant: matrix_size (number of actual rows) = " << matrix_size << endl;
}

	/* The equations for virtual crossings state x_2i = x_{2i+1} and x_{2j-1} = x_2j, where the terminating
	   edges at the crossing are 2i and 2j-1.  Thus, the variables corresponding to the incoming edges
	   are trivial.
	*/
	valarray<int> trivial_variable_flags(2*num_crossings); 
	valarray<int> variable(2*num_crossings);               
	
	
	for (int i=0; i<2*num_crossings; i++)
		trivial_variable_flags[i] = 0;
	
	for (int i=0; i<num_crossings; i++)
	{
		if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
		{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "commutative_automorphism_invariant: crossing " << i << " trivial variables = " << 2*i << " and " << code_table[OPEER][i] << endl;
			trivial_variable_flags[2*i] = 1;
			trivial_variable_flags[code_table[OPEER][i]] = 1;
		}
	}
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "commutative_automorphism_invariant: trivial_variable_flags: ";
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
	
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "commutative_automorphism_invariant: variables: ";
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
	
	Qpmatrix matrix_rep(matrix_size,matrix_size);
	
	for (unsigned int i=0; i<matrix_rep.numrows(); i++)
	for (unsigned int j=0; j<matrix_rep.numcols(); j++)
		matrix_rep[i][j] = Qpolynomial("0");

	/*  The matrix_rep equations are determined from the fact that a commuatative automorphism switch is the
	    sideways map of an F-biquandle F(b,a) = (a_b, b^a), where a^b = \phi(a)-\phi(b)+\psi(b) and 
	    a_b = \psi(a).
	    
                     a_b   b^a          a_b   b^a            b^a   a_b          b^a   a_b
                  2j-1 \   / 2i+1      2i \   / 2j        2j-1 \   / 2i+1       2i\   / 2j
                        \ /                \ /	                \ /                \ /
                         \                  \                    /                  /
                        / \                / \                  / \                / \
                    2i /   \ 2j      2j-1 /   \	2i+1        2i /   \ 2j       2j-1/   \ 2i+1
                      b     a            b     a              a     b            a     b
	    
	    By considering a positive type I and type II crossing, both of which respect the above map F we see
	    that the corresponding equations are
	    
	    
	                                \psi(x_{2j}=x_{2j-1}                                    \psi(x_{2i+1}=x_{2i}
	    \phi(x_{2i})-\phi(x_{2j})+\psi(x_{2j}))=x_{2i+1}    \phi(x_{2j-1})-\phi(x_{2i+1})+\psi(x_{2i+1}))=x_{2j}
	
	    A negative type I and type II crossing respects map F^*(a,b) = (b^a,a_b) so the corresponding equations 
	    are
	    
	    
	    \phi(x_{2j})-\phi(x_{2i})+\psi(x_{2i}))=x_{2j-1}    \phi(x_{2i+1})-\phi(x_{2j-1})+\psi(x_{2j-1}))=x_{2i}
	                                \psi(x_{2i}=x_{2i+1}                                    \psi(x_{2j-1}=x_{2j}
	
	   
	   Note that for the up operation we have (\psi-\phi)(x_p) for some p.
	   
	   There is a possibility that we have to handle the case where F(x_p,x_q) = (x_p,x_r), F(x_p,x_q) = (x_s,x_q),
	   or similarly for F^*, due to a Reidemeister I detour.  We therefore add the (negative) terms corresponding 
	   to the right hand side of the above equations to the matrix_rep entry, rather than setting that entry.
	   	
	*/

	Qpmatrix psi_minus_phi = psi - phi;
	
	
	int row = 0;
	
	for (int i=0; i<num_crossings; i++)
	{
		/* only consider classical crossings */
		if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
			continue;
		
		/* determine the braid crossing type, POSITIVE, NEGATIVE */
		int crossing_type;		
		if (  (code_table[TYPE][i] == generic_code_data::TYPE1 && code_table[LABEL][i] == generic_code_data::POSITIVE)
				||(code_table[TYPE][i] == generic_code_data::TYPE2 && code_table[LABEL][i] == generic_code_data::NEGATIVE)
   					)
			crossing_type = braid_crossing_type::NEGATIVE;
		else
			crossing_type = braid_crossing_type::POSITIVE;

if (debug_control::DEBUG >= debug_control::DETAIL)
{	
	debug << "commutative_automorphism_invariant: crossing " << i;
	if ( crossing_type == braid_crossing_type::NEGATIVE )
		debug << " negative";
	else
		debug << " positive";

	debug <<  ", TYPE " << (code_table[TYPE][i] == generic_code_data::TYPE1? "I" : "II") << endl;
	debug << "commutative_automorphism_invariant: peer code variables p[i] = ";// << code_table[PERM][i];
	debug << "\t2i = " << 2*i << "\t\t2j-1 = " << code_table[ODD_TERMINATING][i];
	debug << "\t2i+1 = " << 2*i+1 << "\t2p[i] = " << code_table[EVEN_ORIGINATING][i] << endl;
	
	debug << "commutative_automorphism_invariant: mapped to\t\t";//v(p[i]) = " << variable[code_table[PERM][i]];
	debug << "\tv(2i) = " << variable[2*i] << "\tv(2p[i]-1) = " << variable[code_table[ODD_TERMINATING][i]];
	debug << "\tv(2i+1) = " << variable[2*i+1] << "\tv(2p[i]) = " << variable[code_table[EVEN_ORIGINATING][i]] << endl;

	debug << "commutative_automorphism_invariant: row " << 2*row << " ";
}		
		/* first do the up action in row 2*row. */
		if (crossing_type == braid_crossing_type::POSITIVE)
		{
			if (code_table[TYPE][i] == generic_code_data::TYPE1)
			{
				set_matrix_N_element(matrix_rep,2*row,code_table[EVEN_TERMINATING][i],phi,0,0,N);
				set_matrix_N_element(matrix_rep,2*row,code_table[EVEN_ORIGINATING][i],psi_minus_phi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row,code_table[ODD_ORIGINATING][i],N);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << code_table[EVEN_TERMINATING][i] << ":phi    ";
	debug << code_table[EVEN_ORIGINATING][i] << ":psi_minus_phi    ";
	debug << code_table[ODD_ORIGINATING][i] << ":-1     " << endl;
}
			}
			else
			{
				set_matrix_N_element(matrix_rep,2*row,code_table[ODD_TERMINATING][i],phi,0,0,N);
				set_matrix_N_element(matrix_rep,2*row,code_table[ODD_ORIGINATING][i],psi_minus_phi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row,code_table[EVEN_ORIGINATING][i],N);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << code_table[ODD_TERMINATING][i] << ":phi    ";
	debug << code_table[ODD_ORIGINATING][i] << ":psi_minus_phi    ";
	debug << code_table[EVEN_ORIGINATING][i] << ":-1     " << endl;
}
			}
		}
		else
		{
			if (code_table[TYPE][i] == generic_code_data::TYPE1)
			{
				set_matrix_N_element(matrix_rep,2*row,code_table[EVEN_ORIGINATING][i],phi,0,0,N);
				set_matrix_N_element(matrix_rep,2*row,code_table[EVEN_TERMINATING][i],psi_minus_phi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row,code_table[ODD_TERMINATING][i],N);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << code_table[EVEN_ORIGINATING][i] << ":phi    ";
	debug << code_table[EVEN_TERMINATING][i] << ":psi_minus_phi    ";
	debug << code_table[ODD_TERMINATING][i] << ":-1     " << endl;
}
			}
			else
			{
				set_matrix_N_element(matrix_rep,2*row,code_table[ODD_ORIGINATING][i],phi,0,0,N);
				set_matrix_N_element(matrix_rep,2*row,code_table[ODD_TERMINATING][i],psi_minus_phi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row,code_table[EVEN_TERMINATING][i],N);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << code_table[ODD_ORIGINATING][i] << ":phi    ";
	debug << code_table[ODD_TERMINATING][i] << ":psi_minus_phi    ";
	debug << code_table[EVEN_TERMINATING][i] << ":-1     " << endl;
}
			}
		}
		
   
if (debug_control::DEBUG >= debug_control::DETAIL) 
	debug << "commutative_automorphism_invariant: row " << 2*row+1 << " ";

		/* now the down action in row 2*row+1 */	
		if (crossing_type == braid_crossing_type::POSITIVE)
		{
			if (code_table[TYPE][i] == generic_code_data::TYPE1)
			{
				set_matrix_N_element(matrix_rep,2*row+1,code_table[EVEN_ORIGINATING][i],psi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row+1,code_table[ODD_TERMINATING][i],N);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << code_table[EVEN_ORIGINATING][i] << ":psi    ";
	debug << code_table[ODD_TERMINATING][i] << ":-1     " << endl;
}
			}
			else
			{
				set_matrix_N_element(matrix_rep,2*row+1,code_table[ODD_ORIGINATING][i],psi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row+1,code_table[EVEN_TERMINATING][i],N);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << code_table[ODD_ORIGINATING][i] << ":psi    ";
	debug << code_table[EVEN_TERMINATING][i] << ":-1     " << endl;
}
			}
		}
		else
		{
			if (code_table[TYPE][i] == generic_code_data::TYPE1)
			{
				set_matrix_N_element(matrix_rep,2*row+1,code_table[EVEN_TERMINATING][i],psi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row+1,code_table[ODD_ORIGINATING][i],N);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << code_table[EVEN_TERMINATING][i] << ":psi    ";
	debug << code_table[ODD_ORIGINATING][i] << ":-1     " << endl;
}
			}
			else
			{
				set_matrix_N_element(matrix_rep,2*row+1,code_table[ODD_TERMINATING][i],psi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row+1,code_table[EVEN_ORIGINATING][i],N);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << code_table[ODD_TERMINATING][i] << ":psi    ";
	debug << code_table[EVEN_ORIGINATING][i] << ":-1     " << endl;
}
			}
		}
	
		row++;
	}
	

	if (braid_control::EXTRA_OUTPUT)
	{
		if (!braid_control::SILENT_OPERATION)
		{
			cout << "\n" << matrix_size << "-square matrix representation of R-module determined by peer code, M:\n";
			cout << matrix_rep << endl;
		}
		
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\nMatrix representation of R-module determined by peer code, M:\n";
			output << matrix_rep << endl;
		}
	}

if (debug_control::DEBUG >= debug_control::BASIC)
{
    debug << "commutative_automorphism_invariant: matrix representation derived from generic code data:" << endl;
    debug << matrix_rep;
    debug << endl;
}


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
		    
	Qpolynomial delta_0 = determinant (matrix_rep, title);

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

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "commutative_automorphism_invariant:Delta_0 = " << delta_0 << endl;

	if (!braid_control::CALCULATE_DELTA_0_ONLY && (braid_control::ALWAYS_CALCULATE_DELTA_1 || delta_0 == Qpolynomial("0")))
	{
		/* hcf is used only in single variable examples */
		polynomial<scalar,char,bigint> hcf = delta_0.getn();

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "\ncommutative_automorphism_invariant:hcf initialized to " << hcf << endl;

		int rperm[matrix_size];
		int cperm[matrix_size];

		/* take out a (N) row and (N) column from the matrix representation 
		   and calculate the determinant of what is left
		*/
		for (int Nr1=0;Nr1<matrix_N_size;Nr1++)
		{
			/* take out the N rows of the matrix rep starting from row Nr1*N row */
			for (int i=0; i<Nr1*N; i++)
    			rperm[i] = i;
			for (int i = Nr1*N+N; i< matrix_size; i++)
				rperm[i-N] = i;


			int Nc2=0; // needed for rat_minor_determinant - relates to long knots;
			for (int Nc1= 0; Nc1< matrix_N_size; Nc1++)
			{
			
				Qpolynomial det = polynomial<scalar,char,bigint> ("0");
				
				/* take out the N columns of the matrix rep starting from column N*j */
				for (int i=0; i<Nc1*N; i++)
	   				cperm[i] = i;

				for (int i=Nc1*N+N; i<matrix_size; i++)
					cperm[i-N] = i;
				
				if (rat_minor_determinant(&matrix_rep, N, Nr1, Nc1, Nc2, rperm, cperm, title, hcf))
					goto rat_poly_delta_1_calculation_complete;
			}			
		}

rat_poly_delta_1_calculation_complete:

		if (!braid_control::SILENT_OPERATION)
			cout << endl;
		
		if (braid_control::NUMERATOR_GCD)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "Delta_1 = " << hcf << endl;	
			
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "numerator gcd  = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << hcf << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "commutative_automorphism_invariant:" << "Delta_1 = " << hcf << endl;	

		}
	}
}
