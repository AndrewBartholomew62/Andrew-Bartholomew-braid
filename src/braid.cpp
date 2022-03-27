/**************************************************************************
	  This module contains the braid brokering function and those
	  functions it calls

void braid(string input_string, string title)
void dynnikov(string input_string, int num_strings, int num_terms)
int braid_to_dowker (string braid, int num_terms, int*& code)
polynomial<int> sawollek(string braid, int num_terms)
int number_of_components (char* braid, int num_terms)
string braid_to_gauss_code (string braid, int num_terms)
string braid_to_generic_code (string braid, int num_terms, int code_type)
hpolynomial homfly(vector<int> braid_num, vector<int> type, vector<int> component_record, 
                                                                       int basepoint, int level, bool virtual_braid)
hpolynomial virtual_homfly(vector<int> braid_num, vector<int> type, vector<int> master_component_record, int level)
int num_braid_terms(string word)
void set_component_record (int action, vector<int>& braid_num, vector<int>& type, vector<int>& component_record,
                           int base_crossing, int basepoint, int datum)
void write_braid (ostream& s, vector<int>& braid_num, vector<int>& type)
void fixed_point_invariant(matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, 
						   matrix<int>& Tu, matrix<int>& Td, ST_pair_type pair_type, string input_string, string title)
int num_fixed_points(matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, 
						   matrix<int>& Tu, matrix<int>& Td, string input_string)
void rack_poly_invariant(matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, matrix<int>& Tu, 
						 matrix<int>& Td, ST_pair_type pair_type, string input_string, string title)

	
 **************************************************************************/
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <valarray>
#include <iomanip>
#include <ctype.h>
#include <stdio.h>

using namespace std;


extern ofstream     debug;
extern ofstream     output;
extern ifstream     input;

#include <util.h>
#include <quaternion-scalar.h>
#include <polynomial.h>
#include <matrix.h>
#include <braid.h>
#include <gauss-orientation.h>
#include <reidemeister.h>

/********************* Function prototypes ***********************/
void help_info();
void tostring (char*& cptr, int n);
polynomial<int> sawollek(string braid, int num_terms);
//bool valid_braid_input (string input_string, int& num_terms, int& num_strings, bool raw_output, bool silent_output, bool output_as_input);
//void parse_braid(string word, int num_terms, vector<int>& braid_num, vector<int>& type);
void dynnikov(string input_string, int num_strings, int num_terms);
int braid_to_dowker (string braid, int num_terms, int*& code);
Hpolynomial study_determinant (const Qpmatrix& H_matrix, int n=0, int* rperm=0, int* cperm=0);
void braid_reduce (vector<int>& braid_num, vector<int>& type, int& basepoint, vector<int>& component_record);
int number_of_components (string braid, int num_terms);
string braid_to_gauss_code (string input, int num_terms);
string braid_to_generic_code (string braid, int num_terms, int code_type);
Qpmatrix C_study_matrix(const Qpmatrix& H_matrix, int n, int* rperm, int* cperm);
Qpmatrix R_study_matrix(const Qpmatrix& C_matrix);
void normalize(Hpolynomial& poly);
hpolynomial homfly(vector<int> braid_num, vector<int> type, vector<int> component_record, 
                                                                 int basepoint, int level, bool virtual_braid);
hpolynomial virtual_homfly(vector<int> braid_num, vector<int> type, vector<int> master_component_record, int level);
void set_component_record (int action, vector<int>& braid_num, vector<int>& type, vector<int>& component_record,
                           int base_crossing=0, int basepoint=0, int datum=0);
void write_braid (ostream& s, vector<int>& braid_num, vector<int>& type);
int num_fixed_points(matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, 
						   matrix<int>& Tu, matrix<int>& Td, string input_string);
bool distinguish_modified_string (matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, 
						   matrix<int>& Tu, matrix<int>& Td, string input_string, int num_terms);
//void flip_braid(string& braid);
//void invert_braid(string& braid);
//void line_reflect_braid(string& braid);
//void plane_reflect_braid(string& braid);
void Kamada_double_covering(string& braid);
void generic_code(string input_string, string title);
//int remove_Reidemeister_II(generic_code_data& code_data, vector<int>& component_flags);


void braid(string input_string, string title)
{
	int num_terms;
	int num_strings;
	bool flip_braid_qualifier = false;
	bool invert_braid_qualifier = false;
	bool plane_reflect_braid_qualifier = false;
	bool line_reflect_braid_qualifier = false;
	
	/* first read and then remove braid qualifiers from the input string */
	string::size_type pos = input_string.find('{');

   	if (pos != string::npos)
   	{
   		if (input_string.substr(pos).find("flip") != string::npos)
   			flip_braid_qualifier = true;

   		if (input_string.substr(pos).find("invert") != string::npos)
   			invert_braid_qualifier = true;

   		if (input_string.substr(pos).find("plane-reflect") != string::npos)
   			plane_reflect_braid_qualifier = true;

   		if (input_string.substr(pos).find("line-reflect") != string::npos)
   			line_reflect_braid_qualifier = true;
    			
   		input_string = input_string.substr(0,pos);

if (debug_control::DEBUG >= debug_control::DETAIL)
  	debug << "braid: after removing braid qualifiers, input_string =  " << input_string << endl;

   	}

	if (braid_control::FLIP_BRAID || flip_braid_qualifier)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "braid: flipping braid " << input_string << endl;
    
				flip_braid(input_string);    
	}

	if (braid_control::INVERT_BRAID || invert_braid_qualifier)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "braid: inverting braid " << input_string << endl;
    
				invert_braid(input_string);    
	}

	if (braid_control::PLANE_REFLECT_BRAID || plane_reflect_braid_qualifier)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "braid: plane-reflecting braid " << input_string << endl;
    
				plane_reflect_braid(input_string);    
	}

	if (braid_control::LINE_REFLECT_BRAID || line_reflect_braid_qualifier)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "braid: line-reflecting braid " << input_string << endl;
    
				line_reflect_braid(input_string);    
	}

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
					
	if (valid_braid_input(input_string, num_terms, num_strings, braid_control::SILENT_OPERATION, braid_control::RAW_OUTPUT, braid_control::OUTPUT_AS_INPUT))
	{
		if (braid_control::DOWKER_CODE)
		{
			int*		dowker_code;
		    /* we use dowker_code to hold the Dowker code once calculated */
	    	num_terms = braid_to_dowker(input_string,num_terms,dowker_code);
	    	if (num_terms == -1)
			{
				if (!braid_control::SILENT_OPERATION)
					cout << "\nDowker code does not support virtual braids.\n";
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Dowker code does not support virtual braids.\n";
				}
	    	}
			else if (num_terms == 0)
			{
				if (!braid_control::SILENT_OPERATION)
					cout << "\nClosure of braid is a link.\n";
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Closure of braid is a link.\n";
				}
			}
	    	else
	    	{
				if (!braid_control::SILENT_OPERATION)
				{
					cout << "\nDowker code = ";
					for (int i=0;i<num_terms;i++)
						cout << dowker_code[i] << " ";
				}
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Dowker code = ";
				}
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
				for (int i=0;i<num_terms;i++)
		    		output << dowker_code[i] << " ";
				if (!braid_control::SILENT_OPERATION)
					cout << "\n";
				output << "\n";
				delete dowker_code;
	    	}
		}
		else if (braid_control::GAUSS_CODE)
		{
			
			string pcode = braid_to_generic_code(input_string, num_terms, generic_code_data::peer_code);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid: pcode from braid word is " << pcode << endl;

			generic_code_data code_data;
			read_peer_code (code_data, pcode);
			
			ostringstream oss;
			
			if (braid_control::LPGD)
			{
				gauss_orientation_data gauss_data(code_data);
				gauss_orientation_data lp_gauss_data = left_preferred(gauss_data); //unoriented = false;
				write_gauss_data(lp_gauss_data,oss); 
			}
			else if (braid_control::ULPGD)
			{
				gauss_orientation_data gauss_data(code_data);
				gauss_orientation_data lp_gauss_data = left_preferred(gauss_data,true, code_data.immersion); //unoriented = true;
				write_gauss_data(lp_gauss_data,oss); 
			}
			else		
				write_gauss_code(oss,code_data);
			
//	    	string gcode = braid_to_gauss_code(input_string, num_terms);

			if (!braid_control::SILENT_OPERATION)
				cout << "\nGauss code = " << oss.str() << endl;
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Gauss code = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << '\n';
			}
			output << oss.str() << endl;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "\nbraid: Gauss code = " << oss.str() << endl;

		}
		else if (braid_control::HOMFLY)
		{
				vector<int> braid_num(num_terms);
			    vector<int> type(num_terms);

				parse_braid(input_string, num_terms, braid_num, type);

				/* For virtual braids we consider all orderings of the braid components.  There is a natural
				   numbering of braid components obtained by considering the lowest strand number of the braid
				   belonging to each component and then considering the components to be ordered according to the 
				   ascending order of their lowest strands.  We set up master_component_record[i] to record the 
				   lowest strand in the braid belonging to component i where i is the natural order of components 
				   in the input braid (in the current scope). 
				*/
				vector<int> master_component_record;
				set_component_record(CR_CREATE, braid_num, type, master_component_record);

				/* call braid_reduce to capture cases where there is an unnecessary singular crossing at the top or at
				   the bottom of the braid.  This avoids incorrect smoothing of braids like "s1" by the homfly function.
				*/
				int dummy_basepoint=1; // only required for the call to braid_reduce
				braid_reduce(braid_num, type, dummy_basepoint, master_component_record);

				int num_cpts = master_component_record.size();

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "braid: master_component_record = ";
	for (int i=0; i< num_cpts; i++)
		debug << master_component_record[i] << " ";
	debug << "\nbraid: num_cpts = " << num_cpts << endl;
}	
				
				hpolynomial hpoly; 
				
				if (input_string.find('t') == string::npos)
				{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "braid: classical braid, using master component_record for call to homfly function" << endl;

					hpoly = homfly(braid_num,type,master_component_record,1,0,false); //basepoint on strand 1, level=0, virtual=false

					if (!braid_control::SILENT_OPERATION)
						cout << "\n\nHOMFLY polynomial = " << hpoly << endl;

					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "HOMFLY polynomial = ";
						if (braid_control::OUTPUT_AS_INPUT)
							output << '\n';
					}
					output << hpoly << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid: HOMFLY polynomial = " << hpoly << endl;

				}
				else
				{
//if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
//	debug << "braid: virtual braid, using master_component_perm for call to virtual_homfly function" << endl;
	
//					hpoly = virtual_homfly(braid_num,type,master_component_record,0); // level=0	

					if (!braid_control::SILENT_OPERATION)
						cout << "\nThe HOMFLY polynomial remains undefined for virtual knots and links" << endl;
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "The HOMFLY polynomial remains undefined for virtual knots and links";
					}
					output << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid: The HOMFLY polynomial remains undefined for virtual knots and links" << endl;					
				}
		}
		else if (braid_control::IMMERSION_CODE)
		{
			string icode = braid_to_generic_code(input_string, num_terms, generic_code_data::immersion_code);
			if (icode.length())
			{
				if (!braid_control::SILENT_OPERATION)
					cout << "\nImmersion code = " << icode << endl;
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Immersion code = ";
					if (braid_control::OUTPUT_AS_INPUT)
						output << '\n';
				}
				output << icode << endl;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "\nbraid: \tImmersion code = " << icode << endl;
			}
			else
			{
				if (!braid_control::SILENT_OPERATION)
					cout << "\nNo immersion code defined, closure of braid is a link." << endl;
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "No immersion code defined, closure of braid is a link." << endl;
				}
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "\n\nbraid: no immersion code defined, closure of braid is a link." << endl;
			}
		}
		else if (braid_control::PEER_CODE)
		{
			string pcode = braid_to_generic_code(input_string, num_terms, generic_code_data::peer_code);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid: \tpcode from braid word is " << pcode << endl;

			generic_code_data code_data;
			read_peer_code (code_data, pcode);

			if (braid_control::REMOVE_REIDEMEISTER_II_MOVES)
			{
				vector<int> dummy_flags; // not tracking components
				remove_Reidemeister_II (code_data,dummy_flags);
			}

			if (!braid_control::SILENT_OPERATION)
			{
				cout << "\nPeer code = ";
				write_peer_code(cout,code_data);
				cout << endl;
			}
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Peer code = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << '\n';
			}
			write_peer_code(output,code_data);
			output << endl;
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "\nbraid: \tPeer code = ";
	write_peer_code(debug,code_data);
	debug << endl;
}
		}
		else if (braid_control::SAWOLLEK)
		{
		    /* call the function sawollek to evaluate the polynomial */
			polynomial<int>	delta = sawollek(input_string, num_terms);
			if (!braid_control::SILENT_OPERATION)
				cout << "\nSawollek's normalized Z-polynomial = " << delta << "\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Sawollek's normalized Z-polynomial = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << delta << "\n";
		}
		else if (braid_control::DYNNIKOV_TEST)
		{
			dynnikov(input_string, num_strings, num_terms);
		}
		else if (braid_control::AFFINE_INDEX)
		{
			/* convert to a labelled peer code and call generic_code, 
			   which will evaluate the affine index polynomial
			*/
			string pcode = braid_to_generic_code(input_string, num_terms, generic_code_data::peer_code);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "\nbraid: peer code from braid word = " << pcode << endl;

			generic_code(pcode, title);

		}
		else if (braid_control::STATUS_INFORMATION)
		{
			int num_cpts = number_of_components(input_string,num_terms);
			
			if (!braid_control::SILENT_OPERATION)
			{
				cout << num_terms << " terms" << "\n";
				cout << num_strings << " strands" << "\n";
				cout << num_cpts << " components" << "\n";
			}
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << num_terms << " terms";
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << num_strings << " strands";
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << num_cpts << " components";
			}
			
		}
		else
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\nUnknown or missing task for braid.\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Unknown or missing task for braid.\n";
			}
			exit(0);
		}
		
	}
} // End of function braid

/* The Dynnikov test determines whether the braid in inbuf is trivial or not */
void dynnikov(string input_string, int num_strings, int num_terms)
{
	char*		mark;
	int			number;
	bigint		b1plus,b1minus,b2plus,b2minus;
	bigint		a1dash,a2dash,b1dash,b2dash;
	bigint		temp;
	bigint		zero(0);
	bool     S;
	bool     inverse;


	/* assign alpha = (0,1,...,0,1) in Z^2n */
    bigint alpha[2*num_strings];

	for (int i=0; i<num_strings; i++)
    {
		alpha[2*i] = 0;
		alpha[2*i+1] = 1;
	}

if (debug_control::DEBUG >= debug_control::BASIC)
{
    debug << "\n\nDynnikov test: initial value of alpha = ";
    for (int i=0;i<2*num_strings; i++)
		debug << alpha[i] << " ";
}
	char* inbuf = c_string(input_string);
    char* cptr = inbuf;
    for (int i = 0; i< num_terms; i++)
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

		mark = cptr; /* mark where we start the number */

		/* look for the end of the number */
		while (isdigit(*cptr))
		    cptr++;

		get_number(number, mark);

if (debug_control::DEBUG >= debug_control::BASIC)
{
    debug << "\nDynnikov test: ";
    if (S)
    {
		if (inverse)
	    	debug << "-";
		else
	    	debug << " ";
		debug << "s";
    }
    else
		debug << " t";
    debug << number << ": ";
}

		if (!S)
		{
		    /* the twist swaps both a_i and a_{i+1}
		       and b_i and b_{i+1} do the a_i first
		    */
		    temp = alpha[2*(number-1)];
		    alpha[2*(number-1)] = alpha[2*number];
		    alpha[2*number] = temp;

		    temp = alpha[2*(number-1)+1];
		    alpha[2*(number-1)+1] = alpha[2*number+1];
		    alpha[2*number+1] = temp;
		}
		else
		{
		    /* calculate b_i^{+} and b_i^{-}
		       and b_{i+1}^{+} and b_{i+1}^{-}
		    */
		    b1plus = b1minus = alpha[2*(number-1)+1];
		    if (b1plus < zero)
				b1plus = zero;
		    if (b1minus > zero)
				b1minus = zero;
		    b2plus = b2minus = alpha[2*number+1];
		    if (b2plus < zero)
				b2plus = zero;
		    if (b2minus > zero)
				b2minus = zero;

		    if (inverse)
		    {
				/* set temp = (a_i - a_{i+1} + b1minus)^{+} */
				temp = alpha[2*(number-1)] - alpha[2*number] + b1minus;
				if (temp < zero)
			    	temp = zero;

				/* a1' = a_i - b1plus - temp */
				a1dash = alpha[2*(number-1)] - b1plus - temp;

				/* set temp = (a_{i+1} - a_i + b2plus)^{-} */
				temp = alpha[2*number] - alpha[2*(number-1)] + b2plus;
				if (temp > zero)
			    	temp = zero;

				/* a2' = a_{i+1} - b2minus - temp */
				a2dash = alpha[2*number] - b2minus - temp;

				/* set temp = (a_i - a_{i+1} + b1minus - b2plus)^{-} */
				temp = alpha[2*(number-1)] - alpha[2*number] + b1minus - b2plus;
				if (temp > zero)
			    	temp = zero;

				/* b1' = b_{i+1} + temp */
				b1dash = alpha[2*number+1] + temp;
				/* b2' = b_i - temp */
				b2dash = alpha[2*(number-1)+1] - temp;
		    }
		    else
		    {
				/* set temp = (a_{i+1} - a_i + b1minus)^{+} */
				temp = alpha[2*number] - alpha[2*(number-1)] + b1minus;
				if (temp < zero)
				    temp = zero;

				/* a1' = a_i + b1plus + temp */
				a1dash = alpha[2*(number-1)] + b1plus + temp;

				/* set temp = (a_i - a_{i+1} + b2plus)^{-} */
				temp = alpha[2*(number-1)] - alpha[2*number] + b2plus;
				if (temp > zero)
			    	temp = zero;

				/* a2' = a_{i+1} + b2minus + temp */
				a2dash = alpha[2*number] + b2minus + temp;

				/* set temp = (a_i - a_{i+1} - b1minus + b2plus)^{+} */
				temp = alpha[2*(number-1)] - alpha[2*number] - b1minus + b2plus;
				if (temp < zero)
			    	temp = zero;

				/* b1' = b_{i+1} - temp */
				b1dash = alpha[2*number+1] - temp;
				/* b2' = b_i + temp */
				b2dash = alpha[2*(number-1)+1] + temp;
		    }

		    /* now update alpha with a1dash, a2dash,
		       b1dash and b2dash
		    */
		    alpha[2*(number-1)] = a1dash;
		    alpha[2*(number-1)+1] = b1dash;
		    alpha[2*number] = a2dash;
		    alpha[2*number+1] = b2dash;
		}
if (debug_control::DEBUG >= debug_control::BASIC)
{
    debug << "Dynnikov test:  yields alpha = ";

    for (int j=0;j<2*num_strings; j++)
		debug << alpha[j] << " ";
}
	}

	/* now determine whether alpha has been changed */
    bool identity = true;
    for (int i=0; i<num_strings; i++)
		alpha[2*i+1] -= 1;

    for (int i=0; i<2*num_strings; i++)
    {
		if (alpha[i] != zero)
		{
		    identity = false;
		    break;
		}
    }

    if (identity)
    {
		if (!braid_control::SILENT_OPERATION)
			cout << "\nis the trivial braid\n";
		output << "\nis the trivial braid\n";
	}
	else
    {
		if (!braid_control::SILENT_OPERATION)
			cout << "\nis not the trivial braid\n";
		output << "\nis not the trivial braid\n";
    }

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    if (identity)
		debug << "\nDynnikov test: is the trivial braid\n";
    else
		debug << "\nDynnikov test: is not the trivial braid\n";
}
	delete[] inbuf;
}

/* the function braid_to_dowker calculates the dowker code for word in braid,
   storing the result at code.  The length of the dowker code is returned
   to the call, this is equal to the number of terms in the braid (the
   number of crossings).  If the braid is a link the function does not
   return anything at code, which is set to NULL, and the function returns
   zero.  If and error is detected in the braid word, or if the braid is
   virtual the function returns -1 and sets code = NULL.

   The dowker code assumes that the knot is alternating, and indicates by
   a negative label that the crossing is an over(under)-crossing when the
   alternating assumption expects an under(over)-crossing.
   
   The following is taken from the help screen of knotscape:
   
	    The notation used for knots is the Dowker notation, which is adapted from a
		notation proposed by Tait in the nineteenth century.  It is also equivalent
		to an enhanced version of the "Gauss" notation.  The idea is to "walk
		round" the knot curve, recording crossover data on the way.  Probably
		the quickest route to understanding the notation is to study the 
		Dowker Code Example, which illustrates the code for a 4-crossing diagram
		of the right-handed trefoil knot.
		
		Suppose that we are given an n-crossing alternating diagram D of a knot K.
		Choose arbitrarily a starting point and direction on K, and then
		travel around the knot K labelling in order from 1 to 2n the points of K
		which project to crossovers of the diagram D.  If we consider that two
		labelled points on K form a pair if they project to the same crossover
		of D, this pairing gives rise to a chart such as that given in the
		Dowker Code Example, where the top row consists of the odd numbers 
		from 1 to 2n-1 in order, and the bottom row is the sequence of corresponding
		even numbers.  Then the Dowker code for the diagram D together with the
		choice of starting point and direction, is this sequence of even numbers. 
		
		We may extend this code to non-alternating diagrams by taking the
		convention that prefixing a number by a minus sign corresponds to
		switching the corresponding crossover of the diagram D.
		
		The code depends in general on the choice of starting point and direction,
		but we can regard the standard code for the diagram as being the code
		which is lexicographically least.  Similarly the standard code for a
		knot type is the lexicographically smallest code, taken over all
		possible diagrams representing the knot type.

*/
int braid_to_dowker (string braid, int num_terms, int*& code)
{
    int i,j;
    int number;
    int string;
    int crossing;
    int next_label;
    int type;
    char* cptr;
    char* mark;
    bool string_found;
   
	/* record the signed value of i for each \sigma_i in braid_num */
    int braid_num[num_terms];
	
	char* inbuf = c_string(braid);
	cptr = inbuf;
	for ( i = 0; i< num_terms; i++)
	{
	    if (*cptr == '-')
	    {
			braid_num[i] = -1;
			cptr++;
	    }
	    else
			braid_num[i] = 1;

        cptr++;
        mark = cptr; /* mark where we start the number */

        /* look for the end of the number */
        while (isdigit(*cptr))
            cptr++;

	    get_number(number, mark);
	    braid_num[i] *= number;
	}
	
	delete[] inbuf;
	
    /* currently only consider real braids */
    if (braid.find('t') != string::npos)
    {
		code = NULL;
		return -1;
    }
    else
    {
		
		/* set out the array in which we will store the code, even
	   	   labels are stored at even locations, odd labels at odd
		   locations
		*/
    	int label[2*num_terms];
		
		for (i=0;i<2*num_terms;i++)
	    	label[i] = 0;

		/* we are going to start by looking for string 1, the next label
	   	   being 1.  The orientation of the crossings is currently unknown,
		   so set type = 0 (it will oscillate between 1 and -1 for over and
		   under crossings).
		*/
		string = 1;
		next_label = 1;
		type = 0;
		crossing = 0;

		do
		{
		    /* search for the involvement of 'string' in the next crossing,
		       starting from 'crossing'
	    	*/
		    string_found = false;
	    	for (i=crossing; i<num_terms;i++)
	    	{
				if (abs(braid_num[i]) == string || abs(braid_num[i]) == string-1 )
				{
				    string_found = true;
				    crossing = i;
				    break;
				}
		    }

		    if (!string_found)
		    {
				/* start looking at the start of the braid */
				for (i=0; i<crossing; i++)
				{
		    		if (abs(braid_num[i]) == string || abs(braid_num[i]) == string-1 )
		    		{
						string_found = true;
						crossing = i;
						break;
		    		}
				}
	    	}

		    /* so now braid_num[crossing] involves string, odd labels
		       are stored at odd locations in label and even labels
		       at even locations, so each crossing has two consecutive
		       places in label.
		    */

	    	if (label[2*crossing + next_label%2])
		    {
				/* we've gone round in a loop - i.e the braid
				   is a link, return zero and set code to NULL.
				*/
				code = NULL;
				return 0;
		    }
	    	else
			{
				label[2*crossing + next_label%2] = next_label;
				next_label++;
			}

	    	/* check to see if the crossing is of the correct type */
		    if(!type)
	    	{
				/* this is the first crossing, set the orientation */
				if (braid_num[crossing] < 0)
			    	type = 1;
				else
		    		type = -1;
	    	}
	    	else
	    	{
				/* we've set the crossing orientation at an earlier,
				   crossing, so multiply type by -1 to indicate the
				   type of crossing we expect at this stage
				*/
				type *= -1;

				/* now assess whether the thread we are following meets
				   a crossing of the wrong type.  If so, and the label we just
				   added was even, multiply that label by -1.  Note
				   that next_label was incremented at the assignment, so if
			   	   the last label allocated was even, it is now odd.

				   The sign enables us to evaluate codes for non-alternating
				   knots.

				   If we arrive on sting i at \sigma_i, then the sign of
				   sigma and type must be different for the type to be as we
				   expect; if we arrive on string i+1, their sign must be
				   the same.  If this is not the case, we have an unexpected
				   crossing.
				*/
				if (  (abs(braid_num[crossing]) == string && braid_num[crossing] * type > 0)
		   			|| (abs(braid_num[crossing]) != string && braid_num[crossing] * type < 0)
				   )
				{
				    /* the crossing is the wrong type */
				    if (next_label%2)
					label[2*crossing] *= -1;
				}
	    	}

		    /* change string to record the effect of the crossing on the
		       thread we are following
		    */
	    	if ( abs(braid_num[crossing]) == string )
	    	   	string++;
	    	else
				string--;

		    /* now look for the new string starting at the next crossing */
		    crossing++;
		} while (next_label <= 2*num_terms);

		/* now create the dowker code array under code */
    	code = new int[num_terms];
		
		for (i=0;i<num_terms;i++)
		{
		    /* find label 2i+1 in label */
	    	for (j=0;j<num_terms;j++)
	    	{
				if (label[2*j+1] == 2*i+1)
				{
		    		/* the corresponding even label
				       is stored one place back
		    		*/
				    code[i] = label[2*j];
				    break;
				}
	    	}
		}
		return num_terms;
    }
}

/* the function sawollek evaluates Sawollek's normalized Conway
   polynomial for braid.
*/
polynomial<int> sawollek(string braid, int num_terms)
{
    char* mark;
    int w_power = 0;
    int string_num;
    int number;
    int crossing;
    int real_crossing;
    int num_real_crossings = 0;
    int* height;
    int* crossing_num;
    bool string_found;

	height = new int[num_terms];
    
    crossing_num = new int[num_terms];
	
	typedef polynomial<int> Polynomial;
	
    /* record the signed height of each crossing in height,
       and the number of the crossing in crossing_num - zero in
       crossing_num will indicate a virtual crossing.  The sign will be
       used to initialize the Sawollek matrix M and to calculate the value
       w_power, the difference between the number of positive and negative
       crossings. After this it the sign will be set to positive.
    */
    char* inbuf = c_string(braid);
    char *cptr = inbuf;
    for (int i = 0; i< num_terms; i++)
    {
		if (*cptr == '-')
		{
	    	height[i] = -1;
	    	cptr++;
		}
		else
		    height[i] = 1;

		/* if this is a virtual crossing set crossing_num[i] to zero,
		   otherwise set it to the post-increment value of num_real_crossings
		*/
		if (*cptr++ == 't')
	    	crossing_num[i] = 0;
		else
		    crossing_num[i] = ++num_real_crossings;

		/* now determine the height */
		mark = cptr; /* mark where we start the number */

		/* look for the end of the number */
		while (isdigit(*cptr))
	    	cptr++;

		get_number(number, mark);
		height[i] *= number;
    }

    /* evaluate w_power the number of positive crossings
       minus the number of negative crossings
    */
    for (int i=0; i<num_terms;i++)
    {
		if (crossing_num[i])
		{
	    	if (height[i] > 0)
				w_power++;
	    	else
				w_power--;
		}
    }

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "\nsawolleck: braid crossing heights\n";
    for (int i=0; i<num_terms; i++)
		debug <<  height[i] << " ";
    debug << "\nsawolleck: braid crossing numbers\n";
    for (int i=0; i<num_terms; i++)
		debug <<  crossing_num[i] << " ";
    debug << "\nsawolleck: difference between positive and negative crossings = " << w_power;
}

    /* assign space for the Sawollek matrix M */
    matrix<Polynomial> Smatrix(2*num_real_crossings,2*num_real_crossings);

    for (int i=0;i<2*num_real_crossings;i++)
    for (int j=0;j<2*num_real_crossings;j++)
    {
		Smatrix[i][j] = Polynomial("0");
    }

    /* Initialize the Sawollek matrix to be the 2Nx2n diagonal
       matrix diag(M1,M2,...,MN) where N is the number of real
       crossings and Mi is either M+ or M- depending upon the
       sign of real crossing i, and M+ and M- are given by

       M+ =   1-x    -y        M- = 0    -x^-1y
		     -xy^-1   0            -y^-1  1-x^-1

    */
    real_crossing = 0;
    for (int i=0; i< num_terms; i++)
    {
		if (crossing_num[i])
		{
	    	/* this is a real crossing */
	    	if (height[i] <0)
	    	{
				/* (*matrixptr)[2*real_crossing][2*real_crossing]
		    		zero from initialization
				*/
				Smatrix[2*real_crossing][2*real_crossing+1] = Polynomial("-x^-1y");
				Smatrix[2*real_crossing+1][2*real_crossing] = Polynomial("-y^-1");
				Smatrix[2*real_crossing+1][2*real_crossing+1] = Polynomial("1-x^-1");

				height[i] *= -1;
	    	}
	    	else
	    	{
				Smatrix[2*real_crossing][2*real_crossing] = Polynomial("1-x");
				Smatrix[2*real_crossing][2*real_crossing+1] = Polynomial("-y");
				Smatrix[2*real_crossing+1][2*real_crossing] = Polynomial("-xy^-1");
				/* (*matrixptr)[2*real_crossing+1][2*real_crossing+1]
		   			zero from initialization
				*/
	    	}
	    	real_crossing++;
		}
    }

    /* Now determine the permutation P and at the same time evaluate
       the matrix M-P
    */
    for (int i=0; i<2*num_real_crossings; i++)
    {
		/* determine the crossing in the braid where we start
		   and set the string number;
		*/
		real_crossing = i/2+1;
		for (int j=0; j<num_terms; j++)
		{
	    	if (crossing_num[j] == real_crossing)
	    	{
				crossing = j;
				break;
		    }
		}

		if (i%2)
		    string_num = height[crossing] + 1;  // upper output string
		else
	    	string_num = height[crossing];  // lower output string

		/* look along braid to find the next real crossing involving
		   string.  We need the do loop because we may switch strands
		   through virtual crossings and need to follow the braid around
		   more than once.
		*/
		string_found = false;
		do
		{
	    	for (int j=crossing+1; j<num_terms; j++)
	    	{
				if (height[j] == string_num)
				{
		    		if (crossing_num[j] == 0)
						string_num++; // follow strand through virtual crossing
		    		else
		    		{
						string_found = true;
						/* subtract non-zero permutation element from M */
						Smatrix[2*(crossing_num[j]-1)][i] -= 1;
						break;
		    		}
				}
				else if (height[j] == string_num-1)
				{
				    if (crossing_num[j] == 0)
						string_num--; // follow strand through virtual crossing
		    		else
		    		{
						string_found = true;
						/* subtract non-zero permutation element from M */
						Smatrix[2*(crossing_num[j]-1)+1][i] -= 1;
						break;
		    		}
				}
	    	}

		    if (!string_found)
	    	{
				/* start looking at the start of the braid */
				for (int j=0; j<=crossing; j++)
				{
				    if (height[j] == string_num)
		    		{
						if (crossing_num[j] == 0)
			    			string_num++; // follow strand through virtual crossing
						else
						{
			    			string_found = true;
			    			/* subtract non-zero permutation element from M */
			    			Smatrix[2*(crossing_num[j]-1)][i] -= 1;
			    			break;
						}
		    		}
		    		else if (height[j] == string_num-1)
		    		{
						if (crossing_num[j] == 0)
			    			string_num--; // follow strand through virtual crossing
						else
						{
			    			string_found = true;
			    			/* subtract non-zero permutation element from M */
			    			Smatrix[2*(crossing_num[j]-1)+1][i] -= 1;
			    			break;
						}
		    		}
				}
	    	}
		} while (!string_found);
    }

if (debug_control::DEBUG >= debug_control::BASIC)
{
    debug << "\n\nsawolleck: Sawollek matrix (M-P)\n";
    debug << Smatrix;
}

    /* now calculate det(Smatrix) we want to include all the
       rows and columns, so set cperm[i] = rperm[i] = i
    */
	if (braid_control::WAIT_SWITCH)
	{
		matrix_control::WAIT_INFO = true;
		matrix_control::wait_threshold = braid_control::wait_threshold;
	}
    Polynomial poly = determinant(Smatrix);
	matrix_control::WAIT_INFO = false;
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "\nsawolleck: det(M-P) = " << poly;
}
    if (poly.nonzero())
    {
		/* now calculate the normalizing factor and include the w_power
		   to create the normalized form of the poly.  The normalizing
		   factor is minus the least power of x in poly, and x will be the
		   first variable, so we can just evaluate this exponent from
		   each place.
		*/
		int c_index = 0;
		if (poly.numvars() && strchr(poly.getvars(),'x'))
		{
	    	pterm<int,int>* pterm_ptr = poly.get_pterms();

	    	int factor = poly.getlength()/(2*(poly.getdeg()[0]+1));
	    	do
	    	{
				int power = pterm_ptr->pl/factor/2;
				if ((pterm_ptr->pl/factor)%2)
		    		power *= -1;

				if (pterm_ptr == poly.get_pterms())
		    		c_index = power;
				else if (power < c_index)
		    		c_index = power;

				pterm_ptr = pterm_ptr->next;

		    } while (pterm_ptr);
		}

		/* create the normalizing polynomial in word, then get_polynomial
		   into norm_factor.  Thus word is (-1^w_power)x^(-c_index), check
		   its not 1 or -1 first;
		*/
		if (c_index == 0)
		{
	    	if (w_power%2)
				poly *= -1;

if (debug_control::DEBUG >= debug_control::BASIC)
{
    debug << "\nsawolleck: normalizing poly= " << (w_power%2? -1 : 1);
}

		}
		else
		{
			ostringstream oss;
			
	    	if (w_power%2)
				oss << '-';
	    	oss <<  'x';
	    	if ( c_index != -1)
	    	{
				oss << '^';
				oss << -c_index;
	    	}

	    	poly *= Polynomial(oss.str());

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
    debug << "\nsawolleck: normalizing poly here= " << Polynomial(oss.str());
}

        }
    }

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "\nsawolleck: Sawollek's normalized Z-polynomial = " << poly << "\n";
}

	delete[] inbuf;
    delete[] height;
    delete[] crossing_num;
	
	return poly;
}

/* number_of_components was created from braid_to_gauss_code so there is code duplication
   here.  Difficult to see how this could have been avoided though....
*/
int number_of_components (string braid, int num_terms)
{
	valarray<int> braid_num(num_terms);
	
	char* inbuf = c_string(braid);
	char* cptr = inbuf; 
	for ( int i = 0; i< num_terms; i++)
	{
	    if (*cptr == '-')
			cptr++;
		
        cptr++;
        char* mark = cptr; /* mark where we start the number */

        /* look for the end of the number */
        while (isdigit(*cptr))
            cptr++;

		int number;
	    get_number(number, mark);
	    braid_num[i] = number;
	}

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "\nnumber_of_components: braid_num: ";
	for (int i=0; i<num_terms; i++)
		debug << braid_num[i] << ' ';
	debug << endl;
}    

	/* The number of components of the braid's closure is bounded by
	   the number of strands, work out the maximum in num_cpts,
	   we'll set it more accurately later.
	*/
	int num_cpts = 0;
	for (int i=0; i< num_terms; i++)
	   if (braid_num[i] > num_cpts)
		   num_cpts = braid_num[i];

	/* the number of strands is one greater than the highest
	   crossing number
	*/
	num_cpts++;

	/* We shall be tracing around the braid working and need to deal
	   with the case that the braid is a link.  We shall set a flag for
	   each lower and upper strand at each crossing that we traverse, so we
	   can detect the end of each link component and find the start of
	   another.  The flags will be a simple array with lower flags stored at
	   even locations and upper flags at odd locations.  We use lower and
	   upper because they are valid for virtual crossings and help in
	   identifying the next string.
	*/
	valarray<int> flags(0,2*num_terms);

	num_cpts = 0;
	bool new_component_found;
	do
	{
		/* look for the next starting  place in the braid, we want an lower
		   or upper strand that we've not yet visited.
		*/
		new_component_found = false;
		int crossing;
		int string_num;
		for (int i=0; i< 2*num_terms; i++)
		{
			if (flags[i] == 0)
			{
				new_component_found = true;

				/* work out which string and at which crossing we have found the
				   start of the next component.
				*/
				crossing = i/2;

				/* the string we want is either the corresponding braid_num
				   or the string above it, depending upon the flag
				*/
				string_num = braid_num[crossing];
				if (i%2)
					string_num++;

				break;
			}
		}

		if (!new_component_found)
			continue;


		/* now trace around the braid starting at the crossing we've
		   just identified
		*/
		bool component_not_finished = true;
		do
		{
			/* search for the involvement of 'string' in the next crossing,
			   starting from 'crossing', the first time around we'll find
			   it immediately.
			*/
			bool string_found = false;
			for (int i=crossing; i<num_terms;i++)
			{
				if (braid_num[i] == string_num || braid_num[i] == string_num-1 )
				{
					string_found = true;
					crossing = i;
					break;
				}
			}

			if (!string_found)
			{
				/* start looking at the start of the braid */
				for (int i=0; i<crossing; i++)
				{
					if (braid_num[i] == string_num || braid_num[i] == string_num-1 )
					{
						string_found = true;
						crossing = i;
						break;
					}
				}
			}

			/* Now braid_num[crossing] involves string; check that we've
			   not visited this crossing on this strand before.  
			*/
			int flag_location = 2*crossing;
			if (braid_num[crossing] != string_num)
				flag_location++;

			if (flags[flag_location])
			{
				/* we've been here before, so we've just completed a
				   component
				*/
				num_cpts++;
				component_not_finished = false;
			}
			else
			{
				/* set the flag to record the fact that we've been here
				   and increase the count of the length of this component
				*/
				flags[flag_location] = 1;

				/* change string to record the effect of the crossing on the
				   thread we are following
				*/
				if ( braid_num[crossing] == string_num )
					string_num++;
				else
					string_num--;

				/* now look for the new string starting at the next crossing */
				crossing++;
			}
		} while (component_not_finished);
	} while (new_component_found);

	delete[] inbuf;
	return num_cpts;
}

/* braid_to_gauss_code determines the Gauss code for virtual links.*/
string braid_to_gauss_code (string braid, int num_terms)
{
	/* record the value of i for each \sigma_i or \tau_i in braid_num
	   and the type of each crossing in type.  The crossing_perm is used 
	   to map between braid crossings, which might include virtual crossings
	   and real crossings.
	*/
	valarray<int> braid_num(num_terms);
    valarray<int> type(num_terms);
	valarray<int> crossing_perm(num_terms);
	int real_crossing_num = 0; // the Gauss code is written with crossings numbered from 1
	
	char* inbuf = c_string(braid);
	char* cptr = inbuf; 
	for ( int i = 0; i< num_terms; i++)
	{
	    if (*cptr == '-')
	    {
			type[i] = braid_crossing_type::NEGATIVE;
			cptr++;
	    }
	    else
			type[i] = braid_crossing_type::POSITIVE;

		if (*cptr == 't' || *cptr == 'T')
			type[i] = braid_crossing_type::VIRTUAL;
		else
			real_crossing_num++;  // so not incremented for virtual crossings
		
		crossing_perm[i] = real_crossing_num;
		
        cptr++;
        char* mark = cptr; /* mark where we start the number */

        /* look for the end of the number */
        while (isdigit(*cptr))
            cptr++;

		int number;
	    get_number(number, mark);
	    braid_num[i] = number;
	}

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "\n\ngauss_code: braid_num: ";
	for (int i=0; i<num_terms; i++)
		debug << braid_num[i] << ' ';
	debug << "\ngauss_code: type: ";
	for (int i=0; i<num_terms; i++)
		debug << type[i] << ' ';
	debug << "\ngauss_code: crossing_perm: ";
	for (int i=0; i<num_terms; i++)
		debug << crossing_perm[i] << ' ';
	debug << endl;
}    

	/* The number of components of the braid's closure is bounded by
	   the number of strands, work out the maximum in num_components,
	   we'll set it more accurately later.
	*/
	int num_cpts = 0;
	for (int i=0; i< num_terms; i++)
	   if (braid_num[i] > num_cpts)
		   num_cpts = braid_num[i];

	/* the number of strands is one greater than the highest
	   crossing number
	*/
	num_cpts++;

	valarray<int> cum_cpt_len(0,num_cpts);

	/* We shall be tracing around the braid working and need to deal
	   with the case that the braid is a link.  We shall set a flag for
	   each lower and upper strand at each crossing that we traverse, so we
	   can detect the end of each link component and find the start of
	   another.  The flags will be a simple array with lower flags stored at
	   even locations and upper flags at odd locations.  We use lower and
	   upper because they are valid for virtual crossings and help in
	   identifying the next string.
	*/
	valarray<int> flags(0,2*num_terms);

	/* set out the array in which we will store the code */
    valarray<int> gcode(0,2*num_terms);
	int place = 0; //used as an index into gcode

	num_cpts = 0;
	bool new_component_found;
	do
	{
		/* look for the next starting  place in the braid, we want an lower
		   or upper strand that we've not yet visited.
		*/
		new_component_found = false;
		int crossing;
		int string_num;
		for (int i=0; i< 2*num_terms; i++)
		{
			if (flags[i] == 0)
			{
				new_component_found = true;

				/* work out which string and at which crossing we have found the
				   start of the next component.
				*/
				crossing = i/2;

				/* the string we want is either the corresponding braid_num
				   or the string above it, depending upon the flag
				*/
				string_num = braid_num[crossing];
				if (i%2)
					string_num++;

				break;
			}
		}

		if (!new_component_found)
			continue;


		/* now trace around the braid starting at the crossing we've
		   just identified
		*/
		bool component_not_finished = true;
		int this_cpt_len = 0;
		do
		{
			/* search for the involvement of 'string' in the next crossing,
			   starting from 'crossing', the first time around we'll find
			   it immediately.
			*/
			bool string_found = false;
			for (int i=crossing; i<num_terms;i++)
			{
				if (braid_num[i] == string_num || braid_num[i] == string_num-1 )
				{
					string_found = true;
					crossing = i;
					break;
				}
			}

			if (!string_found)
			{
				/* start looking at the start of the braid */
				for (int i=0; i<crossing; i++)
				{
					if (braid_num[i] == string_num || braid_num[i] == string_num-1 )
					{
						string_found = true;
						crossing = i;
						break;
					}
				}
			}

			/* Now braid_num[crossing] involves string; check that we've
			   not visited this crossing on this strand before.  If not, and if 
			   this is not a virtual crossing, write the next place in gcode.  
			   
			   If this is a virtual crossing, we do not write anything to the Gauss 
			   code, but do set the flag to record the fact that we've been here.
			   
			   If we have been here before we're at the end of the current 
			   component, so record the length of the component we've just 
			   finished and look for another.
			*/
			int flag_location = 2*crossing;
			if (braid_num[crossing] != string_num)
				flag_location++;

			if (flags[flag_location])
			{
				/* we've been here before, so we've just completed a
				   component
				*/
				if (num_cpts)
					cum_cpt_len[num_cpts] =	cum_cpt_len[num_cpts-1]+this_cpt_len;
				else
					cum_cpt_len[num_cpts] = this_cpt_len;

				num_cpts++;

				component_not_finished = false;
			}
			else
			{
			    if (type[crossing] != braid_crossing_type::VIRTUAL)  
				{
					gcode[place] = crossing_perm[crossing];
					if ( (type[crossing] == braid_crossing_type::POSITIVE && braid_num[crossing] == string_num)
					|| (type[crossing] == braid_crossing_type::NEGATIVE && braid_num[crossing] == string_num - 1) )

						gcode[place]*= -1;

					place++;
					this_cpt_len++;
				}

				/* set the flag to record the fact that we've been here
				   and increase the count of the length of this component
				*/
				flags[flag_location] = 1;

				/* change string to record the effect of the crossing on the
				   thread we are following
				*/
				if ( braid_num[crossing] == string_num )
					string_num++;
				else
					string_num--;

				/* now look for the new string starting at the next crossing */
				crossing++;
			}
		} while (component_not_finished);
	} while (new_component_found);
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "\ngauss_code: gcode: ";
	for (int i=0; i<2*num_terms; i++)
		debug << gcode[i] << ' ';
	debug << "\ngauss_code: number of components = " << num_cpts;
	debug << "\ngauss_code: cumulative component lengths: ";
	for (int i=0; i< num_cpts; i++)
		debug << cum_cpt_len[i] << ' ';
	debug << endl;
}
    
	ostringstream oss;

	for (int i=0; i< num_cpts; i++)
	{
		for (int j=(i?cum_cpt_len[i-1]:0); j<cum_cpt_len[i]; j++)
		{
			oss << gcode[j];
			if (j < cum_cpt_len[i]-1)
				oss << ' ';
		}
		if (i< num_cpts-1)
		{
			oss <<  ',';
			oss << ' ';
		}
	}
	
	/* now write out the crossing types */
	oss << ' ';
	oss << '/';

	for (int i=0; i< num_terms; i++)
	{
		if (type[i] == 1)
			oss << " +";
		else if (type[i] == -1)
			oss << " -";
	}

	delete[] inbuf;
	return oss.str();
}

/* braid_to_generic_code evaluates the labelled code corresponding to a
   braid word.  If the closure of the braid is a link and the requested generic code type is a
   labelled immersion code, the function returns the empty string.
*/
string braid_to_generic_code (string braid, int num_terms, int code_type)
{   
	/* record the value of i for each \sigma_i or \tau_i in braid_num and
	   the type of each crossing, from a braid perspective, in braid_type.
	*/
    vector<int> braid_num(num_terms);
	
    vector<int> braid_type(num_terms);
	
	char* inbuf = c_string(braid);
	char* cptr = inbuf;
	for ( int i = 0; i< num_terms; i++)
	{
	    if (*cptr == '-')
	    {
			braid_type[i] = braid_crossing_type::NEGATIVE;
			cptr++;
	    }
	    else
			braid_type[i] = braid_crossing_type::POSITIVE;

		if (*cptr == 't' || *cptr == 'T')
			braid_type[i] = braid_crossing_type::VIRTUAL;

		cptr++;
		char* mark = cptr; /* mark where we start the number */

        /* look for the end of the number */
        while (isdigit(*cptr))
            cptr++;

	    get_number(braid_num[i], mark);
	}

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "braid_to_generic_code: braid_num: ";
	for (int i=0; i<num_terms; i++)
		debug << braid_num[i] << ' ';
	debug << "\nbraid_to_generic_code: braid_type: ";
	for (int i=0; i<num_terms; i++)
		debug << braid_type[i] << ' ';
	debug << endl;
}    

     matrix<int> code_table(CODE_TABLE_SIZE,num_terms);

	/* trace through the braid writing the odd and even peers into the code table
	   and noting the component for each naming edge.  We start the numbering of
	   the first component from the string 1 at the start of the braid, the first edge
	   being 0, to which 'edge' was initialized.
	   
	   Since we have to convert the numbering of the crosings from braid crosings to
	   peer code crossings, we begin by writing the edge labels into the EVEN_TERMINATING
	   and ODD_TERMINATING rows of the code table using braid crossing numbers, then we
	   sort the odd and even numbered peers into the OPEER and EPEER rows based on peer
	   code crossing numbers.  We record the braid crossing number of the ith peer-code 
	   crossing in braid_crossing_num.
	   
	   We initialize the EVEN_TERMINATING and ODD_TERMINATING rows of the code table to -1
	   to determine the components of the braid.  Since each component of an immersed link contains 
	   an even number of edges, when we complete a component we shall attempt
	   to write an even edge to a location of the code_table that has already had an edge assigned
	   (so the value will be positive).  We then need to find the start of another component, 
	   if one exists.  We look along the EVEN_TERMINATING row of the code_table to see if there is
	   another crossing for which we have not assigned an even terminating peer (i.e one for which
	   the entry is negative.  If there is one, we choose the first one for which we have already
	   assigned an odd terminating peer (i.e whose ODD_TRMINATING entry is positive).  We are 
	   guaranteed to be able to find such a crossing if there is another component to number. 
	   
	   To continue the numbering we need to understand on which strand the new component starts.
	   To do this we record, as a negative number, the peer strand number in code_table[EVEN_TERMINATING][i]
	   when we write a value to code_table[ODD_TERMINATING][i], if the even numbered peer is not yet known.
	*/
    vector<int> braid_crossing_num(num_terms);
	int edge=0;
	int string_num = 1;
	int crossing = 0;
	int component = 0;
	
	for (int i=0; i< num_terms; i++)
		code_table[EVEN_TERMINATING][i] = code_table[ODD_TERMINATING][i] = -1;

	bool complete = false;
	do
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid_to_generic_code: immersion edge " << edge << endl;
	
	    /* search for the next crossing involving 'string_num' in the ,
	       starting from 'crossing'
    	*/
	    bool string_found = false;
    	for (int i=crossing; i<num_terms;i++)
    	{
			if (braid_num[i] == string_num || braid_num[i] == string_num-1 )
			{
			    string_found = true;
			    crossing = i;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid_to_generic_code:  string " << string_num << " found at braid crossing " << i << endl;
			    break;
			}
	    }

	    if (!string_found)
	    {
			/* start looking at the start of the braid */
			for (int i=0; i<crossing; i++)
			{
				if (braid_num[i] == string_num || braid_num[i] == string_num-1 )
				{
					string_found = true;
					crossing = i;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid_to_generic_code:  string " << string_num << " found at braid crossing " << i << endl;
					break;
	    		}
			}
    	}

		/* so now braid_num[crossing] involves string_num. */
		int row = (edge%2? ODD_TERMINATING: EVEN_TERMINATING);
		

		if (code_table[row][crossing] >= 0)
		{
			/* we've reached the end of a component (and edge will therefore be even).  We
		       look along the EVEN_TERMINATING row of the code_table to see if there is
		       another crossing for which we have not assigned an even terminating peer.
		       If there is one, we choose the first one for which we have already assigned
		       an odd terminating peer.  We are guaranteed to be able to find such a crossing
		       if there is another component to number. 
		    */
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid_to_generic_code:  end of component detected" << endl;
		    complete = true;
			for (int i=0; i< num_terms; i++)
		    {
				if (code_table[EVEN_TERMINATING][i] < 0 && code_table[ODD_TERMINATING][i] >= 0)
				{
					component++;
					crossing = i;
					string_num = abs(code_table[EVEN_TERMINATING][i]);
					complete = false;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "braid_to_generic_code:  start of component " << component << " detected at crossing " << crossing << endl;
	debug << "braid_to_generic_code:  starting string number = " << string_num <<  endl;
}
					break;
				}
			}
		}
		else
		{
			/* continue with the current component */
			code_table[row][crossing] = edge;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid_to_generic_code:    setting code_table[" << (edge%2? "ODD_TERMINATING": "EVEN_TERMINATING") << "][" << crossing << "] = " << edge << endl;
			
			if (edge%2 == 0)
			{
				int peer_code_crossing = edge/2;
				braid_crossing_num[peer_code_crossing] = crossing;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid_to_generic_code:    setting braid_crossing_num[" << peer_code_crossing << "] = " << crossing << endl;
	
				/* record the component of the peer-code naming edge */
				code_table[COMPONENT][peer_code_crossing] = component;
			}
			else
			{
				/* write the string number on which the even numbered peer appears into 
				   the EVEN_TERMINATING row if it has not yet been assigned an edge number
				*/
				if (code_table[EVEN_TERMINATING][crossing] < 0)
				{
					code_table[EVEN_TERMINATING][crossing] = (braid_num[crossing] == string_num? string_num+1: string_num-1);
					code_table[EVEN_TERMINATING][crossing] *= -1;
				}
			}
	
			/* change string to record the effect of the crossing on the
			   thread we are following and record the type (TYPE1 or TYPE2)
			   and the label (POSITIVE, NEGATIVE, VIRTUAL) for the crossing
			   for the immersion code at the same time
			*/
			if ( braid_num[crossing] == string_num )
			{
				string_num++;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid_to_generic_code:    change to string " << string_num << endl;
				if (edge%2 == 0)
				{
					int peer_code_crossing = edge/2;
					code_table[TYPE][peer_code_crossing] = generic_code_data::TYPE1; 
	
					if (braid_type[crossing] == braid_crossing_type::VIRTUAL)
						code_table[LABEL][peer_code_crossing] = generic_code_data::VIRTUAL;
					else if (braid_type[crossing] == braid_crossing_type::NEGATIVE)
						code_table[LABEL][peer_code_crossing] = generic_code_data::POSITIVE; 
					else
						code_table[LABEL][peer_code_crossing] = generic_code_data::NEGATIVE; 
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "braid_to_generic_code:    crossing TYPE = TYPE1" << endl;
	debug << "braid_to_generic_code:    crossing LABEL = " << code_table[LABEL][peer_code_crossing] << endl;
}
				}
			}
			else
			{
				string_num--;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid_to_generic_code:    change to string " << string_num << endl;
				if (edge%2 == 0)
				{
					int peer_code_crossing = edge/2;
					code_table[TYPE][peer_code_crossing] = generic_code_data::TYPE2;
	
					if (braid_type[crossing] == braid_crossing_type::VIRTUAL)
						code_table[LABEL][peer_code_crossing] = braid_crossing_type::VIRTUAL;
					else if (braid_type[crossing] == braid_crossing_type::POSITIVE)
						code_table[LABEL][peer_code_crossing] = generic_code_data::POSITIVE;
					else
						code_table[LABEL][peer_code_crossing] = generic_code_data::NEGATIVE;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "braid_to_generic_code:    crossing TYPE = TYPE1" << endl;
	debug << "braid_to_generic_code:    crossing LABEL = " << code_table[LABEL][peer_code_crossing] << endl;
}
				}
			}
	
	
		    /* now look for the new string starting at the next crossing, which
			   will be our next edge
			*/
			edge++;
		    crossing++;
		}
	} while (!complete);

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "braid_to_generic_code: braid_crossing_num (conversion from peer code numbers to braid numbers): ";
	for (int i=0; i<num_terms; i++)
		debug << braid_crossing_num[i] << ' ';
	debug << "\nbraid_to_generic_code: EVEN_TERMINATING (braid order): ";
	for (int i=0; i<num_terms; i++)
		debug << code_table[EVEN_TERMINATING][i] << ' ';
	debug << "\nbraid_to_generic_code: ODD_TERMINATING  (braid order): ";
	for (int i=0; i<num_terms; i++)
		debug << code_table[ODD_TERMINATING][i] << ' ';
	debug << "\nbraid_to_generic_code: crossing types (peer code order): ";
	for (int i=0; i<num_terms; i++)
		debug << code_table[TYPE][i] << ' ';
	debug << "\nbraid_to_generic_code: labels (peer code order): ";
	for (int i=0; i<num_terms; i++)
		debug << code_table[LABEL][i] << ' ';
	debug << "\nbraid_to_generic_code: component (peer code order): ";
	for (int i=0; i<num_terms; i++)
		debug << code_table[COMPONENT][i] << ' ';
	debug << endl;
}

	/* write the odd and even peers into the OPEER and EPEER rows */
	for (int i=0; i< num_terms; i++)
	{
		code_table[OPEER][i] = code_table[ODD_TERMINATING][braid_crossing_num[i]];
		code_table[EPEER][(code_table[OPEER][i]-1)/2] = 2*i;
	}
	
	generic_code_data code_data;
	
	if (code_type == generic_code_data::peer_code)
		code_data.type = generic_code_data::peer_code;
	else if (code_type == generic_code_data::immersion_code)
		code_data.type = generic_code_data::immersion_code;
		
	code_data.code_table = code_table;
	code_data.num_crossings = num_terms;
	code_data.num_components = component+1;
	
//if (debug_control::DEBUG >= debug_control::DETAIL)
//	print (code_data, debug, "braid_to_generic_code: ");

	string result;
	if (code_type == generic_code_data::immersion_code && code_data.num_components > 1)
		result = "";  //actuall it initializes to this of course but it makes the logic explicit
	else
	{
		ostringstream oss;
		write_code_data(oss,code_data);
		result = oss.str();
	}
	
	delete[] inbuf;
	return result;
}

/* The homfly function uses the skein relation
		
		a^-1 P(L+) - a P(L-) = z P(L0)
		P(unknot) = 1

   where L+ is a positive crossing, L- a negative crossing and L0 a (Seifert) smoothed crossing
   to determine the homfly polynomial of a supplied braid word, using the number of bad crossings
   as the induction variable.  This choice of skein relationship was made to be consistent with the
   polynomials given by Chuck Livingston's calculator at http://www.indiana.edu/~knotinfo/.
   
   A crossing is bad if it is first encountered on an under-arc, and is good
   if it is first encountered on an over-arc.  If a classical braid has no bad crossings, following the 
   braid according to its orientation means one is always decending, i.e. it is the unlink.
   
   Consideration of the above skein relationship with the crossings and smoothing obtained from the immersion of 
   a figure 8 in the plane shows that P(O O) = D = z^-1(a^-1 - a)  and that P(U_n) = D^{n-1}, where U_n is the unlink
   comprised of n unknotted components.
   
   Note that the smoothing operation used by the homfly function does not change the component record, which means 
   it does not smooth singular crossings (Reidemeister I loops) at the top or bottom of a braid correctly.  Such braids 
   (e.g. s1, s2-s1s3s2-s3, s1-s2s3s1-s2 etc.) need to be braid reduced before being processed by the homfly function.
   This will always be the case for recursive calls but means that the code that first calls the homfly function 
   should make sure that the first bad crossing is not a singular crossing.

   Note also that the homfly function will not be passed a braid with only one term by a recursive call, due to the 
   braid_reduce function, which will identify it as the unknot.  Therefore when the homfly function passes a smoothed
   braid to the braid_reduce function it still contains at least one crossing.

*/
hpolynomial homfly(vector<int> braid_num, vector<int> type, vector<int> component_record, 
                                                                       int basepoint, int level, bool virtual_braid)
{
	if (!braid_control::SILENT_OPERATION && braid_control::WAIT_SWITCH && level <= braid_control::wait_threshold)
		cout << "." << flush;

	int num_cpts=component_record.size();
	int num_terms=braid_num.size();

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "presented with input braid ";
	write_braid(debug,braid_num,type);
	debug << " at level " << level << endl;
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "input braid has " << num_terms << " terms" << endl;

	if (debug_control::DEBUG >= debug_control::DETAIL)
	{
		debug << "homfly:(" << level << ") ";
		for (int i=0; i< level; i++)
			debug << "  ";
		debug << "braid_num = ";
		for ( int i=0; i< num_terms; i++)
			debug << braid_num[i] << " ";
		debug << "\nhomfly:(" << level << ") ";
		for (int i=0; i< level; i++)
			debug << "  ";
		debug << "     type = ";
		for ( int i=0; i< num_terms; i++)
			debug << type[i] << " ";
		debug << endl;
	}

	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "number of components = " << num_cpts << endl;
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "component_record = ";
	for (unsigned int i=0; i< component_record.size(); i++)
		debug << component_record[i] << " ";
	debug << "\nhomfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "initial basepoint on strand " << basepoint << endl;
}

	/* determine the number of strings in the braid */
	int num_strings = 0;
	if (num_terms)
	{
		for (int i=0; i< num_terms; i++)
		{
			if (braid_num[i] > num_strings)
				num_strings = braid_num[i];
		}
		num_strings++;
	}

	/* There may be other strings above the upper crossing in braid_num that are recorded in component_record.  
	   In this case the number of strings is the largest value stored in the component_record
    */
	for (int i=0; i<num_cpts; i++)
	{
		if (component_record[i] > num_strings)
			num_strings = component_record[i];
	}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "num_strings = " << num_strings << endl;
}
	/* Look for a bad crossing.  Potentially we have to trace each strand of the braid, although the order in which this is done is 
	   unknown, since we will start at the string indicated by the basepoint and follow the braid circuit.  Then, if we return to 
	   the start of a circuit we use the component record to choose the next strand.  To control this transit we maintain a record 
	   of the processed strands.  We also record the crossings visited, to handle the case that the braid is a link.  Link require
	   special treatment: as we descend the skein tree all crossings on a component may become good and we have to move the 
	   basepoint on to the next component (the order being determined by the component record).  In this case we have to mark all
	   crossings on components appearing earlier in the crossing record to the one containing the basepoint as good.  Moreover, we
	   cannot assume that the basepoint strand is the lowest strand belonging to a component, in the case of knot 9_35 the 
	   braid s1s1s2-s1s2s2s3-s2-s2s4-s3s2s4s3 is reduced to -s1-s1s3-s2s1s3s2 during the skein descent, which has just one component 
	   but the basepoint lies on strand 2.
	   
	   To handle this situation we maintain a second set of flags temp_crossing_visited and temp_processed_strand that record which 
	   crossings and strands we've visited as we trace around a component.  Then, if we haven't encountered the basepoint when the 
	   component is completed, we transcribe these temporary records into crossing_visited and processed_strand.
	   
	*/
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "starting search for a bad crossing" << endl;
}
	int bad_crossing=0;
	valarray<int> processed_strand(0,num_strings);
	valarray<int> temp_processed_strand(0,num_strings);
	valarray<int> crossing_visited(0,num_terms);
	valarray<int> temp_crossing_visited(0,num_terms);
	
	int cpt_index = 0; // component index
	int current_cpt_strands=0; // counts strands on a component
	
	bool passed_basepoint = false;

	int string_num = component_record[cpt_index];
	if (string_num == basepoint)
	{			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "  passing basepoint on strand " << string_num << endl;
}
		passed_basepoint = true;
	}


	for (int i=0; i< num_strings; i++)
	{
		if (processed_strand[string_num-1] || temp_processed_strand[string_num-1])
		{
		
			/* If we've not passed the basepoint, transcribe temp_crossing_visited and temp_processed_strand
			   into crossing_visited and and processed_strand clear the temporary records
			*/
			if (!passed_basepoint)
			{
				for (int j=0; j< num_terms; j++)
				{
					crossing_visited[j] += temp_crossing_visited[j];
					temp_crossing_visited[j] = 0;					
				}

				for (int j=0; j< num_strings; j++)
				{
					processed_strand[j] += temp_processed_strand[j];
					temp_processed_strand[j] = 0;					
				}
			}
			
			/* move to the next component in the order determined by component_record */
			cpt_index++;
			string_num = component_record[cpt_index];
			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "  completed component after " << current_cpt_strands << " strands" << endl;
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "  starting new component with strand " << string_num << endl;
}
			current_cpt_strands = 0;

		} 
		
		if (string_num == basepoint)
		{			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "  passing basepoint on strand " << string_num << endl;
}
			passed_basepoint = true;
				
			/* We may be in the middle of a component when we find the basepoint, in which case the
			   values recorded in component_record for the current component will be discarded, since 
			   we have to trace the component from the basepoint.   We also reduce i by current_cpt_strands
			   to allow the strands of the current component already processed to be processed again, this 
			   time 'from the basepoint', rather than from the the strand in component_record (which may be 
			   different).
			   
			   Note we have to clear the temporary records due to the conditional clause above and below that 
			   look at both the temporary and full records.
			*/
			for (int j=0; j< num_terms; j++)
				temp_crossing_visited[j] = 0;
			for (int j=0; j< num_strings; j++)
				temp_processed_strand[j] = 0;

			i -= current_cpt_strands;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "  basepoint detected midway through component, allowing " << current_cpt_strands << " additional strands to be processed" << endl;
}
		}

		/* record the fact that we've processed this strand */
		if (passed_basepoint)
			processed_strand[string_num-1] = 1;
		else
			temp_processed_strand[string_num-1] = 1;
		
		current_cpt_strands++;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "  starting processing of strand " << string_num << endl;
}

		/* now trace along the braid strand */
		for (int j=0; j< num_terms; j++)
		{
			/* does this crossing involve 'string_num' ?	*/
			if (braid_num[j] == string_num || braid_num[j] == string_num-1 )
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "    next crossing involving strand " << string_num << " is " << j+1 << endl;
}
				if (!crossing_visited[j] && !temp_crossing_visited[j])
				{		
					if ( type[j] != braid_crossing_type::VIRTUAL  && // we're at a real crossing
					     passed_basepoint && // otherwise we know it's good
			    	     ( (braid_num[j] == string_num && type[j] == braid_crossing_type::POSITIVE) ||
					      (braid_num[j] != string_num && type[j] == braid_crossing_type::NEGATIVE) // we've just arrived on an under-arc
                         )
					   )
					{
						/* this is a bad crossing, note that cpt_index will indicate the component 
						   on which this bad crossing was encountered.
						*/
						bad_crossing = j+1;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "bad crossing detected at crossing " << bad_crossing << endl;
}
						break;
					}			

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	if (passed_basepoint)
	{
		debug << "homfly:(" << level << ") ";
		for (int i=0; i< level; i++)
			debug << "  ";
		debug << "    crossing is virtual or good" << endl;
	}
	else
	{
		debug << "homfly:(" << level << ") ";
		for (int i=0; i< level; i++)
			debug << "  ";
		debug << "    crossing is virtual or potentially good, pending encounter with basepoint" << endl;
	}
}

					/* record the fact that we've visited this crossing and move on */
					if (passed_basepoint)
						crossing_visited[j] = 1;
					else
						temp_crossing_visited[j] = 1;
				}

				/* change string to record the effect of the crossing on the strand we are following */
				if ( braid_num[j] == string_num )
					string_num++;
				else
					string_num--;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "    new strand number is " << string_num << endl;
}

			}
		}

		if (bad_crossing)
			break;
	}
			
	if (!bad_crossing)
	{
//		hpolynomial delta = hpolynomial("D");
//		hpolynomial delta = hpolynomial("az^-1-a^-1z^-1");
		hpolynomial delta = hpolynomial("a^-1z^-1-az^-1");
		hpolynomial result = hpolynomial("1");
			
		for (int i=0; i<num_cpts-1; i++)
			result *= delta;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "no bad crossings detected, number of components = " << num_cpts << endl;
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "returning " << result << endl;
}
			
		return result;
	}
	else
	{
		bad_crossing--; // adjust so it refers to the correct offset in braid_num and type

		bool bad_crossing_positive = (type[bad_crossing] == braid_crossing_type::POSITIVE);

		/* change the type of the bad crossing for the good braid */
		if (bad_crossing_positive)
			type[bad_crossing] = braid_crossing_type::NEGATIVE;
		else
			type[bad_crossing] = braid_crossing_type::POSITIVE;
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "bad_crossing_positive =  " << (bad_crossing_positive? "true" : "false") << endl;
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "good braid before reduction: ";
	write_braid(debug,braid_num,type);
	debug << endl;

	if (debug_control::DEBUG >= debug_control::DETAIL)
	{
		debug << "homfly:(" << level << ") ";
		for (int i=0; i< level; i++)
			debug << "  ";
		debug << "braid_num = ";
		for ( unsigned int i=0; i< braid_num.size(); i++)
			debug << braid_num[i] << " ";
		debug << "\nhomfly:(" << level << ") ";
		for (int i=0; i< level; i++)
			debug << "  ";
		debug << "     type = ";
		for ( unsigned int i=0; i< type.size(); i++)
			debug << type[i] << " ";
		debug << endl;
	}
}
		/* now take out redundant terms etc. but save the current braid data 
		   before we do in order to calculate the smoothed braid
		*/
		vector<int> smoothed_braid_num = braid_num;
		vector<int> smoothed_type = type;
		vector<int> smoothed_component_record = component_record;
		
		int good_basepoint = basepoint;
		braid_reduce(braid_num, type, good_basepoint, component_record);
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "good braid after reduction: ";
	write_braid(debug,braid_num,type);
	debug << endl;


	if (debug_control::DEBUG >= debug_control::DETAIL)
	{
		debug << "homfly:(" << level << ") ";
		for (int i=0; i< level; i++)
			debug << "  ";
		debug << "braid_num = ";
		for ( unsigned int i=0; i< braid_num.size(); i++)
			debug << braid_num[i] << " ";
		debug << "\nhomfly:(" << level << ") ";
		for (int i=0; i< level; i++)
			debug << "  ";
		debug << "     type = ";
		for ( unsigned int i=0; i< type.size(); i++)
			debug << type[i] << " ";
		debug << endl;
	}

	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "good basepoint = " << good_basepoint << endl;
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "good component_record = ";
	for (unsigned int i=0; i< component_record.size(); i++)
		debug << component_record[i] << " ";
	debug << "\nhomfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "starting evaluation of smoothed braid" << endl;

}
		/* To calculate the smoothed braid we have to remove the bad crossing completely.
		   Note!  This smoothing operation does not change the component record, which means 
		   it does not smooth singular crossings at the top or bottom of a braid correctly.
		*/
	
		smoothed_braid_num.erase(smoothed_braid_num.begin()+bad_crossing);
		smoothed_type.erase(smoothed_type.begin()+bad_crossing);
		
		/* Determine the component record for the smoothed braid.  Smoothing a crossing can fuse two components or can
		break a component into two, but the number of strands is not affected and therefore the original coponent_record 
		remains accurate for the start of the smoothed braid.  In the case that smoothing splits a component into two, we
		need to identify the component on which the bad crossing was first encountered, which is recorded by cpt_index.
		*/
		int smoothed_basepoint = basepoint;
		set_component_record(CR_SMOOTHED_CROSSING, smoothed_braid_num, smoothed_type, smoothed_component_record,
		                     0,smoothed_basepoint, cpt_index);
		int num_smoothed_cpts = smoothed_component_record.size();

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "smoothed braid before reduction: ";
	write_braid(debug,smoothed_braid_num,smoothed_type);
	debug << endl;
	
	if (debug_control::DEBUG >= debug_control::DETAIL)
	{
		debug << "homfly:(" << level << ") ";
		for (int i=0; i< level; i++)
			debug << "  ";
		debug << "braid_num = ";
		for ( unsigned int i=0; i< smoothed_braid_num.size(); i++)
			debug << smoothed_braid_num[i] << " ";
		debug << "\nhomfly:(" << level << ") ";
		for (int i=0; i< level; i++)
			debug << "  ";
		debug << "     type = ";
		for ( unsigned int i=0; i< smoothed_type.size(); i++)
			debug << smoothed_type[i] << " ";
		debug << endl;
	}
	
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "number of smoothed components = " << num_smoothed_cpts << endl;
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "smoothed_component_record = ";
	for (unsigned int i=0; i< smoothed_component_record.size(); i++)
		debug << smoothed_component_record[i] << " ";
	debug << endl;						
}

		/* take out redundant terms etc.*/
		braid_reduce(smoothed_braid_num, smoothed_type, smoothed_basepoint, smoothed_component_record);
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "smoothed braid after reduction: ";
	write_braid(debug,smoothed_braid_num,smoothed_type);
	debug << endl;

	if (debug_control::DEBUG >= debug_control::DETAIL)
	{
		debug << "homfly:(" << level << ") ";
		for (int i=0; i< level; i++)
			debug << "  ";
		debug << "braid_num = ";
		for ( unsigned int i=0; i< smoothed_braid_num.size(); i++)
			debug << smoothed_braid_num[i] << " ";
		debug << "\nhomfly:(" << level << ") ";
		for (int i=0; i< level; i++)
			debug << "  ";
		debug << "     type = ";
		for ( unsigned int i=0; i< smoothed_type.size(); i++)
			debug << smoothed_type[i] << " ";
		debug << endl;
	}
	
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "smoothed basepoint = " << smoothed_basepoint << endl;
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "smoothed component_record = ";
	for (unsigned int i=0; i< smoothed_component_record.size(); i++)
		debug << smoothed_component_record[i] << " ";
	debug << endl;						
}


if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "recursing with good braid from level " << level << endl;
//	write_braid(debug,braid_num,type);
//	debug << " from level " << level << endl;
}

		hpolynomial hgood = (virtual_braid? 
								virtual_homfly(braid_num,type,component_record,level+1) :
								homfly(braid_num,type,component_record,good_basepoint,level+1,false)
							);

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "HOMFLY polynomial of good braid ";
	write_braid(debug,braid_num,type);
	debug << " at level " << level << " = " << hgood << endl;
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "recursing with smoothed braid from level " << level << endl;
//	write_braid(debug,smoothed_braid_num,smoothed_type);
//	debug << " from level " << level << endl;
}
		hpolynomial hsmoothed = (virtual_braid?
							virtual_homfly(smoothed_braid_num,smoothed_type,smoothed_component_record,level+1) :
							homfly(smoothed_braid_num,smoothed_type,smoothed_component_record,smoothed_basepoint,level+1,false)
								);
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "HOMFLY polynomial of smoothed braid ";
	write_braid(debug,smoothed_braid_num,smoothed_type);
	debug << " at level " << level << " = " << hsmoothed << endl;
}

		hpolynomial result;
		
		if (bad_crossing_positive)
		{
		    result = hpolynomial("z")*hsmoothed;
			result += hpolynomial("a")*hgood;
			result *= hpolynomial("a");

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "HOMFLY polynomial at level " << level << " = a * (a * hgood + z * hsmoothed)" << endl;
}
		}
		else
		{
		    result = hpolynomial("-z")*hsmoothed;
			result += hpolynomial("a^-1")*hgood;
			result *= hpolynomial("a^-1");

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "HOMFLY polynomial at level " << level << " = a^-1 * (a^-1 * hgood - z * hsmoothed)" << endl;
}
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "level " << level << " returning " << result << endl;
}
		return result;				
	}
}

/* virtual_homfly is an attempt to extend the definition of the HOMFLY polynomial to the virtual case,
   to date this has proved unsuccessful.
   
   For virtual braids we consider all orderings of the braid components, starting from the order
   presented at the call in the master_component_record.  To consider all orderings we maintain a 
   permutation component_perm that records a permutation of this natural order.
*/
hpolynomial virtual_homfly(vector<int> braid_num, vector<int> type, vector<int> master_component_record, int level)
{
	int num_cpts = master_component_record.size();
	vector<int> component_perm(num_cpts);
	for (int i=0; i< num_cpts; i++)
		component_perm[i] = i;

	hpolynomial vhpoly = hpolynomial("1");
	do
	{
		/* set the component_record to match the component_perm by permuting the master_component_record*/
		vector<int> component_record(num_cpts);
		for (int i=0; i< num_cpts; i++)
			component_record[i] = master_component_record[component_perm[i]];

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "virtual_homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "component_perm = ";
	for (int i=0; i< num_cpts; i++)
		debug << component_perm[i] << " ";
	debug << "\nvirtual_homfly:(" << level << ") ";
	for (int i=0; i< level; i++)
		debug << "  ";
	debug << "component_record = ";
	for (int i=0; i< num_cpts; i++)
		debug << component_record[i] << " ";
	debug << endl;	
}  
		/* when we call homfly, look for a bad_crossing starting from the lowest strand of the first component */
		vhpoly *= homfly(braid_num,type,component_record,component_record[0],level,true); //virtual=true
	} while (next_permutation(component_perm.begin(), component_perm.end()));
	
	return vhpoly;
}


int num_braid_terms(string word)
{
	int num_terms;
	int num_strings;
	
	if (word.length() == 0)
		return 0;
		
	if (!valid_braid_input(word, num_terms, num_strings, braid_control::SILENT_OPERATION, braid_control::RAW_OUTPUT, braid_control::OUTPUT_AS_INPUT))
		exit(0);
	
	return num_terms;
}

/* There is a natural numbering of braid components obtained by considering the lowest strand number of the braid
   belonging to each component and then considering the components to be ordered according to the ascending order 
   of their lowest strands.  Thus component_record[i] records the lowest strand in the braid belonging to component 
   i of this natural order.

   The operation of set_component_record is determined by the value of the action parameter, which takes one of the following 
   values:
   
   		CR_CREATE
		CR_SMOOTHED_CROSSING
		CR_RESET_TO_CROSSING
		CR_UNTWIST_CROSSING
		CR_INVERSE_PAIR
		  
   If the action is CR_CREATE, the component_record size should be zero at the call, and the function will create a 
   component record from braid_num and type.  In other cases the component_record  must be the original component_record 
   for the braid and the function will re-set the component_record to reflect the new braid structure otherwise.  The 
   function assumes braid crossings are numbered from 1.

   If the action is CR_SMOOTHED_CROSSING, the function resets the component record following the smoothing of a crossing.
   It traces around each component in turn from the original lowest strand recording which strands have been visited in 
   processed_strand and noting when we encounter the basepoint.  If we start a component from the original component record 
   and find that we have already visited that strand, we have fused two components by smoothing a crossing.  If when we have
   completed the original component_record we have unprocessed strands, we have broken a component into two; here there are 
   two cases to consider.  If we have already encountered the basepoint then the unprocessed strands represent the new component 
   and we can add the corresponding new entry to the end of the component record.  If we have not encountered the basepoint then 
   it must lie on one of the unprocessed strands and we have already processed the new component.  In this latter case we have to 
   move the new component to the end of the component record and insert the value for basepoint-component's component record in 
   the approriate location.  This location is passed as datum (it indicates the component on which the smoothed crossing was 
   originally encountered when it was identified as a bad crossing, since this is the component in the braid that is split 
   into two by the smoothing).
   
   If the action is CR_RESET_TO_CROSSING the function simply re-evaluates the component_record at the crossing specified by 
   datum in the braid.  That is, it determines the lowest strand of each component in front of the given crossing.
   
   If the action is CR_UNTWIST_CROSSING, the function resets the component_record after un-twisting the crossing given in datum, 
   a singlular lower crossing, and effectively removing the lower strand from the braid. Note that there may be other unlinked 
   unknotted strands below the lower crossing.  In this case the function traces around the component that starts on strand 
   braid_num[crossing] (we have to look for which component that is, since it might not be the one with index 0 in component_record)
   and resets the corresponding component_record entry to the second lowest strand, since the strand is effectively removed by the 
   untwisting.  It then reduces each entry in the component record by one to reflect the removal of the strand.
   
   If the action is CR_INVERSE_PAIR, the function resets the component record for the case where the base crossing (and 
   therefore the current component_record) lies between an inverse pair of crossings.  Here the value of the datum parameter 
   is the value i where s_i/-s_i or t_i/t_i are involved in the inverse pair.  If the two strands are on different components 
   and if a component_record value indicates we start a component on either of these two strands then it starts on the 
   other one after the inverse pair is removed.  If the two strands are on the same component, the component_record 
   remains unchanged after the pair is removed.
void set_component_record (int action, vector<int>& braid_num, vector<int>& type, vector<int>* component_record,
                           int base_crossing=0, int basepoint=0, int datum=0);
*/
void set_component_record (int action, vector<int>& braid_num, vector<int>& type, vector<int>& component_record,
                           int base_crossing, int basepoint, int datum)
{
	int num_terms = braid_num.size();
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "set_component_record: action = ";
	if (action == CR_CREATE)
		debug << "CR_CREATE" << endl;
	else if (action == CR_SMOOTHED_CROSSING)
		debug << "CR_SMOOTHED_CROSSING" << endl;
	else if (action == CR_RESET_TO_CROSSING)
		debug << "CR_RESET_TO_CROSSING" << endl;
	else if (action == CR_UNTWIST_CROSSING)
		debug << "CR_UNTWIST_CROSSING" << endl;
	else if (action == CR_INVERSE_PAIR)
		debug << "CR_INVERSE_PAIR" << endl;
	debug << "set_component_record: braid_num = ";
	for ( unsigned int i=0; i< braid_num.size(); i++)
		debug << braid_num[i] << " ";
	debug << "\nset_component_record:      type = ";
	for ( unsigned int i=0; i< type.size(); i++)
		debug << type[i] << " ";
	debug << "\nset_component_record: num_terms = " << num_terms << endl;
	debug << "set_component_record: current_component_record: ";
	for ( unsigned int i=0; i< component_record.size(); i++)
		debug << component_record[i] << " ";
	debug << "\nset_component_record: base_crossing = " << base_crossing << endl;
	debug << "set_component_record: basepoint = " << basepoint << endl;
	debug << "set_component_record: datum = " << datum << endl;
}
	
	if (num_terms == 0)
		return;

	int fused_cpt_index=0;	
	int first_inverse_pair_strand=0;
	int second_inverse_pair_strand=0;
	int num_cpts = component_record.size();
	
	int num_strings = 0;
	for (int i=0; i< num_terms; i++)
	{
		if (braid_num[i] > num_strings)
			num_strings = braid_num[i];
	}
	num_strings++;
	
	/* There may be other strings above the upper crossing in braid_num that are recorded in component_record.  
	   In this case the number of strings is the largest value stored in the component_record
    */
	for (int i=0; i<num_cpts; i++)
	{
		if (component_record[i] > num_strings)
			num_strings = component_record[i];
	}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: num_strings = " << num_strings << endl;

	valarray<int> processed_strand(0,num_strings);
	valarray<int> crossing_visited(0,num_terms);
	vector<int> new_component_record = component_record;

	int string_num;
	int cpt_index=0;
	
	if (action == CR_CREATE)
	{
		/* we're creating a component record for the braid, record the lowest strand of the first component */
		string_num = 1;
		component_record.push_back(1);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: creating component record, starting first component on strand 1" << endl;

	}
	else if (action == CR_UNTWIST_CROSSING)
	{
		string_num = braid_num[datum-1];
		
		/* set cpt_index to string_num's component in the component record and set the record to num_strings.  
		   Since we know the current starting strand for this component we don't need it again and can re-set 
		   it to the second lowest strand for the comonent in the exitsing component_record.
		*/
		for (; cpt_index < num_cpts; cpt_index++)
		{
			if (component_record[cpt_index] == string_num)
			{
				component_record[cpt_index] = num_strings;
				break;
			}
		}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "set_component_record: untwisting lower crossing " << datum << endl;
	debug << "set_component_record: braid_number of lower crossing and initial string_num = " << string_num << endl;
	debug << "set_component_record: cpt_index = " << cpt_index << ", component record for this component initialized to " << num_strings << endl;
}
	}
	else
	{
		string_num = component_record[cpt_index];

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "set_component_record: re-setting component_record ";
	if (action == CR_RESET_TO_CROSSING)
		debug << "for crossing " << datum << endl;
	else if (action == CR_SMOOTHED_CROSSING)
		debug << "after smoothing a crossing" << endl;
	else if (action == CR_INVERSE_PAIR)
		debug << "for the middle of an inverse pair" << endl;
}
	}
	
	if (action == CR_RESET_TO_CROSSING) 
	{
		for (int i=0; i< num_cpts; i++)
			new_component_record[i] = num_strings;
			
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "set_component_record: new_component_record initialized to: ";
	for (unsigned int i=0; i< new_component_record.size(); i++)
		debug << new_component_record[i] << " ";
	debug << endl;
}
			
	}

	bool passed_basepoint = false;
	if (string_num == basepoint)
	{

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: passing basepoint on strand " << string_num << endl;

		passed_basepoint = true;
	}
		
	for (int i=0; i< num_strings; i++) // this only controls the number of iterations, not the order of strand processing
	{
		if (processed_strand[string_num-1])
		{
			if (action == CR_UNTWIST_CROSSING)
			{
				/* we have finished tracing around the lower crossing's component, we reduce each component_record value
				   above the untwisted strand by one outside the for loop, since there might be only one component,
				   in which case we don't reach this point
				*/
				break;
			}
			else
			{
				if (action == CR_CREATE)
				{
					/* look for the next unprocessed strand, note that we are guaranteed to find one because the 
					   fact that we're in the for loop means there's another string to process.  When we find it, we push 
					   the lowest strand number onto the back of the component_record
					*/
					for (int j=0; j< num_strings; j++)
					{
						if (processed_strand[j] == 0)
						{
							string_num = j+1;
							component_record.push_back(string_num);
							cpt_index++; // used only to set processed_strand correctly

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: new component found starting with strand " << string_num << endl;
							break;
						}
					}
				}
				else
				{
					if (action == CR_SMOOTHED_CROSSING)
					{
						/* we're resetting after smoothing a crossing, which may have fused or split off a component.  If
						   it splits off a new component, we may have already reached the end of the component record.
						   If we have already have checked all the known components, we have a new one (because we're still 
						   in the for loop and so have unprocessed strands).  In this case we have to search for the new 
						   component, otherwise we can determine the next component's starting strand directly.		
						*/
						cpt_index++;

						if (cpt_index == num_cpts)
						{
							/* we have split component datum+1 into two, first find the lowest unprocessed strand */

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: component " << datum+1 << " has been split by smoothing the crossing" << endl;

							for (int j=0; j< num_strings; j++)
							{
								if (processed_strand[j] == 0)
								{
									string_num = j+1;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: new component found starting with strand " << string_num << endl;
									break;
								}
							}
							
							if (passed_basepoint)
							{						
								/* the new component goes at the end of the component_record */
								component_record.push_back(string_num);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: basepoint already encountered, adding " << string_num 
	      << " to the end of the component record" << endl;

							}
							else
							{
								/* datum2 is the offset of the new component, move it to the end and replace it */
								component_record.push_back(component_record[datum]);
								component_record[datum] = string_num;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: basepoint not yet encountered, moving " << component_record[datum] 
	      << " to the end of the component record" << endl;
							}
							cpt_index++; // used only to set processed_strand correctly
							
						}
						else
						{
							string_num = component_record[cpt_index];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: checking for new component on strand " << string_num << endl;

							if (processed_strand[string_num-1])
							{
								/* the current component has been fused with an earlier one, so we record the current
								   component number.  We remove string_num from the component_record outside of the for 
								   loop, since it's possible that the last component is fused, in which case we won't 
								   capture it here because the for loop will have already terminated.
								*/
								fused_cpt_index = cpt_index;

								/* There must be another component to process, since we've not finished processing all of the 
								   strands (we're still in the for loop), and there can be only one fused component.
								*/
								string_num = component_record[++cpt_index];
																
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "set_component_record: component " << fused_cpt_index << " has been fused with another component" << endl;
	debug << "set_component_record: moving to next component starting on strand " << string_num << endl;
}
							}
							else
							{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record:   valid start of next component" << endl;
							}
						}
					}
					else
					{
						/* move to the next component in the component_record. */					
						string_num = component_record[++cpt_index];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: moving to next component starting on strand " << string_num << endl;
					}
				}
			}
		}
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: starting processing of strand " << string_num
	      << ", cpt_index = " << cpt_index << ", from crossing " << base_crossing+1 << endl;

		if (string_num == basepoint)
		{

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: passing basepoint on strand " << string_num << endl;

			passed_basepoint = true;
		}

		/* record the fact that we've processed this strand we set the value to the current component number*/
		processed_strand[string_num-1] = cpt_index+1;
		
		if (action == CR_UNTWIST_CROSSING && string_num < component_record[cpt_index] && string_num != braid_num[datum-1])
		{
			/* check the component_record for the component to see if it is smaller 
			   than the current value but not the original lowest value
			*/
			component_record[cpt_index] = string_num;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: reducing record for cpt_index " << cpt_index << " to " << string_num << endl;

		}
		else if (action == CR_INVERSE_PAIR)
		{
			/* check whether this strand is involved in the inverse pair of crossings 
			   and note the component number either for the first or second time */
			if (string_num == datum || string_num == datum+1)
			{
				if (first_inverse_pair_strand == 0)
				{
					first_inverse_pair_strand = string_num;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: first_inverse_pair_strand = " << first_inverse_pair_strand << endl;
				}
				else
				{
					second_inverse_pair_strand = string_num;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: second_inverse_pair_strand = " << second_inverse_pair_strand << endl;
				}
			}
		}
		
		/* now trace around the braid strand from the base_crossing */
		for (int j=0; j< num_terms; j++)
		{
			int k = (j+base_crossing)%num_terms;
			
			/* set the new_component_record, if required */
			if (action == CR_RESET_TO_CROSSING && k==datum-1)
			{					
				/* check whether the component record should be reduced */
				if (string_num < new_component_record[cpt_index])
				{
					new_component_record[cpt_index] = string_num;					

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record:   reducing new component record of component " << cpt_index+1 << " to " << string_num << endl;
				}
			}

			/* does this crossing involve 'string_num' ? */
			if (braid_num[k] == string_num || braid_num[k] == string_num-1 ) 
			{

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record:   next crossing involving strand " << string_num << " is " << k+1 << endl;
		
				/* change string to record the effect of the crossing on the strand we are following */
				if ( braid_num[k] == string_num )
					string_num++;
				else
					string_num--;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record:   new strand number is " << string_num << endl;
			}
		}
	}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "set_component_record: strand to component mapping after tracing braid: ";
	for (unsigned int i=0; i< processed_strand.size(); i++)
		debug << processed_strand[i] << " ";
	debug << endl;
}


	if (action == CR_UNTWIST_CROSSING)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "set_component_record: after untwisting, component_record changed to: ";
	for (unsigned int i=0; i< component_record.size(); i++)
		debug << component_record[i] << " ";
	debug << endl;
	debug << "set_component_record: decrementing all component_record values above " << braid_num[datum-1] << endl;
}
		/* we reduce each component_record value above the untwisted strand by one 
		   to reflect the fact that a strand has been removed and we're done 
		*/
		for ( int i=0; i< num_cpts; i++)
		{
			if (component_record[i] > braid_num[datum-1])
				component_record[i]--;
		}
	}
	else if (action == CR_SMOOTHED_CROSSING && (fused_cpt_index != 0 || cpt_index < num_cpts-1 ))
	{
		/* we are resetting after smoothing a crossing and have fused two components.  If the fused_cpt_index has been
		   set, we have detected a fused crosing as we've processed the strands, if we've not got to the end of the 
		   component_record but have processed all the strands then we've fused the last component.
		*/
		if (fused_cpt_index)
			cpt_index = fused_cpt_index;
		else
			cpt_index++; // i.e cpt_index = num_cpts-1
			
		int other_cpt = processed_strand[component_record[cpt_index]-1];
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: smoothing crossing has fused component " << cpt_index+1 << " with " << other_cpt <<  endl;
	
		if (component_record[other_cpt-1] > component_record[cpt_index])
		{
			component_record[other_cpt-1] = component_record[cpt_index];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: component " << cpt_index+1 << " had a lower component_record, adjusting component " 
	      << other_cpt << endl;
		}

		component_record.erase(component_record.begin()+cpt_index);
	}
	else if (action == CR_RESET_TO_CROSSING)
	{
		component_record = new_component_record;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_component_record: setting component_record to new_component_record" << endl;
	}
	else if (action == CR_INVERSE_PAIR)
	{	
		if (processed_strand[first_inverse_pair_strand-1] == processed_strand[second_inverse_pair_strand-1])
		{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "set_component_record: inverse pair involves strands on the same component, no adjustment required" << endl;
		}
		else
		{
		
			for (unsigned int j=0; j<component_record.size(); j++)
			{
				if (component_record[j] == first_inverse_pair_strand)
					component_record[j] = second_inverse_pair_strand;
				else if (component_record[j] == second_inverse_pair_strand)
					component_record[j] = first_inverse_pair_strand;
			}		

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "set_component_record: checked component_record for inverse pair strands " 
	      << first_inverse_pair_strand << " and " << second_inverse_pair_strand << endl;
	
		}
	}


if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "set_component_record: component_record changed to: ";
	for (unsigned int i=0; i< component_record.size(); i++)
		debug << component_record[i] << " ";
	debug << endl;	
}
	
}

void write_braid (ostream& s, vector<int>& braid_num, vector<int>& type)
{
	if (braid_num.size())
	{
		for (unsigned int i=0; i< braid_num.size(); i++)
		{
			if (type[i] == braid_crossing_type::NEGATIVE)
				s << "-";

			if (type[i] == braid_crossing_type::VIRTUAL)
				s << "t";
			else
				s << "s";

			s << braid_num[i];
		}
	}
	else
		s <<"\"\"";
}

void fixed_point_invariant(matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, 
						   matrix<int>& Tu, matrix<int>& Td, ST_pair_type pair_type, string input_string, string title)
{		
	int num_strings;
	int num_terms;
	bool flip_braid_qualifier = false;
	bool invert_braid_qualifier = false;
	bool plane_reflect_braid_qualifier = false;
	bool line_reflect_braid_qualifier = false;
	bool double_braid_qualifier = false;

	if (valid_braid_input(input_string, num_terms, num_strings, braid_control::SILENT_OPERATION, braid_control::RAW_OUTPUT, braid_control::OUTPUT_AS_INPUT))
	{	

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
   
		/* If the braid qualifier "welded" is present, it is to be considered as a welded braid, 
		   in which case we need an essential welded pair.
    	   
    	   If the braid qualifier "doodle" is present, it is to be considered as a doodle braid, 
		   in which case we need an essential doodle pair.

    	   If the braid qualifier "flat" is present, it is to be considered as a flat braid, 
		   in which case we need a flat essential doodle pair.
		   
		   If no qualifiers are present the function requires an essential virtual pair.
    	   
    	   If the braid qualifier "flip" is present we'll flip the braid before calculating any
    	   fixed point invariants.  This braid qualifier allows us to flip individal braids
    	   without having to use the global FLIP_BRAID flag.
    	   
    	   If the braid qualifier "double" is present, in addition to "doodle" or "flat", we calculate 
    	   the Kamada double covering of the braid before calculating the number of fixed points.
    	*/
    	bool welded_braid = false;
    	bool doodle_braid = false;
    	bool flat_braid = false;
    	
    	string::size_type pos = input_string.find('{');
    	if (pos != string::npos)
    	{
    		if (input_string.substr(pos).find("welded") != string::npos)
    			welded_braid = true;

    		if (input_string.substr(pos).find("doodle") != string::npos)
    			doodle_braid = true;

    		if (input_string.substr(pos).find("flat") != string::npos)
    			flat_braid = true;

    		if (input_string.substr(pos).find("flip") != string::npos)
    			flip_braid_qualifier = true;

    		if (input_string.substr(pos).find("invert") != string::npos)
    			invert_braid_qualifier = true;

    		if (input_string.substr(pos).find("line-reflect") != string::npos)
    			line_reflect_braid_qualifier = true;

    		if (input_string.substr(pos).find("plane-reflect") != string::npos)
    			plane_reflect_braid_qualifier = true;

    		if (input_string.substr(pos).find("double") != string::npos)
				double_braid_qualifier = true;
    			
    		input_string = input_string.substr(0,pos);

if (debug_control::DEBUG >= debug_control::DETAIL)
  	debug << "fixed_point_invariant: after removing braid qualifiers, input_string =  " << input_string << endl;

    	}
    
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
   	debug << "fixed_point_invariant: welded_braid =  " << (welded_braid?"true":"false") << endl;
   	debug << "fixed_point_invariant: doodle_braid =  " << (doodle_braid?"true":"false") << endl;
   	debug << "fixed_point_invariant: flat_braid =  " << (flat_braid?"true":"false") << endl;
}

		/* check that the S and T we've been given are appropriate for the input braid */
		if (welded_braid && pair_type != ST_pair_type::ESSENTIAL_WELDED)
    	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
   	debug << "fixed_point_invariant: terminating since S and T are not an essential welded pair, as required for welded braids" << endl;
    
			if (!braid_control::SILENT_OPERATION)
				cout << "\nWelded braid requires essential weldeed pair, skipping" << endl;
    
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Welded braid requires essential welded pair, skipping" << endl;
			}
    
			return;
    	}
		else if (doodle_braid && pair_type != ST_pair_type::ESSENTIAL_DOODLE)
    	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
   	debug << "fixed_point_invariant: terminating since S and T are not an essential doodle pair, as required for doodle braids" << endl;
    
			if (!braid_control::SILENT_OPERATION)
				cout << "\nDoodle braid requires essential doodle pair, skipping" << endl;
    
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Doodle braid requires essential doodle pair, skipping" << endl;
			}
    
			return;
    	}
		else if (flat_braid && pair_type != ST_pair_type::FLAT_ESSENTIAL_VIRTUAL)
    	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
   	debug << "fixed_point_invariant: terminating since S and T are not a flat essential virtual pair, as required for flat braids" << endl;
    
			if (!braid_control::SILENT_OPERATION)
				cout << "\nFlat braid requires flat essential virtual pair, skipping" << endl;
    
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Flat braid requires flat essential virtual pair, skipping" << endl;
			}
    
			return;
    	}
		else if (!welded_braid && !doodle_braid && !flat_braid && pair_type != ST_pair_type::ESSENTIAL_VIRTUAL && pair_type != ST_pair_type::FLAT_ESSENTIAL_VIRTUAL)
    	{
			/* virtual braids may use either essential virtual pairs or flat essential virtual pairs */
if (debug_control::DEBUG >= debug_control::SUMMARY)
   	debug << "fixed_point_invariant: terminating since S and T are not an essential virtual pair, as required by the input braid" << endl;
    
			if (!braid_control::SILENT_OPERATION)
				cout << "\ninput braid requires an essential virtual pair, skipping" << endl;
    
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "input braid requires an essential virtual pair, skipping" << endl;
			}
    
			return;
    	}

		if (braid_control::BIGELOW_KNOT_SEARCH)
		{
			/* The Bigelow knot search was used to check braids of the form b=b_ix_1...x_k where b_i is one of the 
			   Bigelow braids (thought clearly any braid could be used), x_i is either s_i, -s_i or t_i and k is 
			   the number of strings in b_i.  For each braid we check whether b, bb, or bbb is a knot and if so 
			   whether it can be distinguished from the unknot via a fixed point invariant.
			*/
		
			/* set up an array of num_strings integers to control the additional crossing types */
			int crossing_type[num_strings];
			for (int i=0; i<num_strings; i++)
				crossing_type[i] = braid_crossing_type::VIRTUAL;
			
			bool found;	
			do
			{
				/* test the combination determined by crossing_type */
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "fixed_point_invariant: crossing_types: ";
	for (int i=0; i< num_strings; i++)
	{
		if (crossing_type[i] != braid_crossing_type::NEGATIVE)
			debug << ' ';
		debug << crossing_type[i] << ' ';
	}
	debug << endl;
}		
		
				/* create the additional word to append to the braid */
				ostringstream oss;
				for (int i=0; i< num_strings; i++)
				{
					if (crossing_type[i] == braid_crossing_type::NEGATIVE)
						oss << '-';
						
					if (crossing_type[i] == braid_crossing_type::VIRTUAL)
						oss << 't';
					else
						oss << 's';
						
					oss << i+1;
				}
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "fixed_point_invariant:   additional word: " << oss.str() << endl;
			
				string modified_input_string = input_string + oss.str();
				bool distinguished = false;
				int count = 1;
				
		        /* test modified_input_string */
		        distinguished = distinguish_modified_string (Su,Sd,invSu,invSd,Tu,Td,modified_input_string,num_terms+num_strings);
		        
				if (!distinguished)
				{
					/* test double modified_input_string */
					distinguished = distinguish_modified_string (Su,Sd,invSu,invSd,Tu,Td,
					                                             modified_input_string+modified_input_string,
					                                             2*(num_terms+num_strings));
				    if (distinguished)
				    {
						count = 2;
					}
					else
					{
						/* test treble modified_input_string */
						distinguished = distinguish_modified_string (Su,Sd,invSu,invSd,Tu,Td,
						                                             modified_input_string+modified_input_string+modified_input_string,
															         3*(num_terms+num_strings));
						if (distinguished)
						{
							count = 3;
						}
					}		
				}
				
				if (distinguished)
				{
					if (!braid_control::SILENT_OPERATION)
					{
						cout << "\nadding " << oss.str();
						if (count == 2)
							cout << " and doubling";
						else if (count == 3)
							cout << " and trebling";
					
					
						cout << " produces a knot that is distinguished from the unknot" << endl;
//						cout <<"\nNumber of fixed-points = " << fixed_points << endl;
					}
					
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "adding " << oss.str();
						if (count == 2)
							output << " and doubling";
						else if (count == 3)
							output << " and trebling";
							
						output << " produces a knot that is distinguished from the unknot";
		//				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
		//				output << "Number of fixed points = ";
					}
		//			output << fixed_points << endl;
					output << endl;
		    
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "fixed_point_invariant:   distinguishable from unknot" << endl;
		    
					return;
				}
		
				/* check whether there is another combination to try,
				   increment crossing_type as if it were an integer, 
				   assuming VIRTUAL < POSITIVE < NEGATIVE, but starting 
				   from the left
				*/
				found = false;
				int place = 0;
				do
				{
					if (crossing_type[place] == braid_crossing_type::VIRTUAL)
					{
						crossing_type[place] = braid_crossing_type::POSITIVE;
						found = true;
						break;
					}
					else if (crossing_type[place] == braid_crossing_type::POSITIVE)
					{
						crossing_type[place] = braid_crossing_type::NEGATIVE;
						found = true;
						break;
					}
					else
					{
						crossing_type[place] = braid_crossing_type::VIRTUAL;
						place++;
					}
					
				} while (place < num_strings);
				
			}while (found);
		}
		else
		{ 
			
			if (braid_control::FLIP_BRAID || flip_braid_qualifier)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "fixed_point_invariant: flipping braid " << input_string << endl;
    
				flip_braid(input_string);
    
			}

			if (braid_control::INVERT_BRAID || invert_braid_qualifier)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "fixed_point_invariant: inverting braid " << input_string << endl;
    
				invert_braid(input_string);
    
			}
			
			if (braid_control::LINE_REFLECT_BRAID || line_reflect_braid_qualifier)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "fixed_point_invariant: line-reflecting braid " << input_string << endl;
    
				line_reflect_braid(input_string);
    
			}


			if (braid_control::PLANE_REFLECT_BRAID || plane_reflect_braid_qualifier)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "fixed_point_invariant: plane-reflecting braid " << input_string << endl;
    
				plane_reflect_braid(input_string);
    
			}

			if ( (braid_control::KAMADA_DOUBLE_COVERING || double_braid_qualifier) && (doodle_braid || flat_braid) )
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "fixed_point_invariant: evaluating the Kamada double covering of braid " << input_string << endl;
    
				Kamada_double_covering(input_string);
    
			}
			
			int fixed_points = num_fixed_points(Su,Sd,invSu,invSd,Tu,Td,input_string);
	    	
			if (!braid_control::SILENT_OPERATION)
				cout <<"\nNumber of fixed-points = " << fixed_points << endl;
	    
			if (!braid_control::RAW_OUTPUT)
	    	{
	    		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
	    		output << "Number of fixed points = ";
	    	}
	    	output << fixed_points << endl;
		    
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "fixed_point_invariant: total number of fixed points = " << fixed_points << endl;
    
    
			/* test code 26-10-13 used to search for virtual knots with non-invariant 
			   fixed-point invariant under the flip isotopy
			
			flip_braid(input_string);
	    	fixed_points = num_fixed_points(Su,Sd,invSu,invSd,Tu,Td,input_string);
			output << fixed_points << endl;
			 end of test code 26-10-13 */
			    
		}
    }
}
					

/* num_fixed_points takes a braid as input string and evaluates the number of fixed points for
   that braid determined by S and T.  All checking of the appropriate use of S and T for input_string
   is done by fixed_point_invariant.
*/
int num_fixed_points(matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, 
						   matrix<int>& Tu, matrix<int>& Td, string input_string)
{
	int num_strings;
	int num_terms;
	int n = Su.numcols();
	
	if (valid_braid_input(input_string, num_terms, num_strings, braid_control::SILENT_OPERATION, braid_control::RAW_OUTPUT, braid_control::OUTPUT_AS_INPUT))
	{	
		char* str = c_string(input_string);
		
		/* work out the braid numbers for each term and identify whether the crossig is virtual or not */
		int braid_num[num_terms];
		bool virtual_crossing[num_terms];

		bool inverse;
    	char* cptr = str;
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
			    virtual_crossing[i] = false;
			else
			    virtual_crossing[i] = true;
				
			cptr++;

			char* mark = cptr; /* mark where we start the number */

			/* look for the end of the number */
			while (isdigit(*cptr))
		    	cptr++;

			get_number(braid_num[i], mark);
			if (inverse)
				braid_num[i] *= -1;
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
  	debug << "num_fixed_points: provided with input string " << input_string << endl;
	debug << "num_fixed_points: braid_num: ";
	for (int i = 0; i< num_terms; i++)
		debug << braid_num[i] << ' ';
	debug << "\nnum_fixed_points: virtual_crossing: ";
	for (int i = 0; i< num_terms; i++)
		debug << virtual_crossing[i] << ' ';
	debug << endl;
}


		/* The fixed point invariant calculates the number of fixed points Mv=v where v ranges through all
		   vectors of length num_strings over X_n and M is the cumulative effect on v of the S_i and T_i
		   determined by the braid word given as an input string.
		*/
		int v[num_strings];
		for (int i=0;i<num_strings; i++)
			v[i] = 0;

		int fixed_points = 0;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "num_fixed_points: calculating fixed points" << endl;

//output << "\n\n";

		bool not_finished = true;
		do
		{
//			for (int i=0;i<num_strings; i++)
//				cout << v[i] << ' ';
//			cout << endl;
	
			/* Evaluate the effect of the input braid on v to determine Mv.  Determine the braid 
			   number bn for the crossing and then apply the appropriate switch (S invS or T) to 
			   strands bn-1 and bn.  Thus for S we have S(x,y) = (y^x,x_y) so 
			   
			   S(Mv[bn], Mv[bn-1]) = (Mv[bn-1]^Mv[bn], Mv[bn]_Mv[bn-1])
			   					   = (Su[Mv[bn]][Mv[bn-1]],Sd[Mv[bn-1]][Mv[bn]])
								   
			   but for invS we have invS(x,y) = (y_\bar{x},x^\bar{y}), so
			   			   
			   invS(Mv[bn], Mv[bn-1]) = (Mv[bn-1]_\bar{Mv[bn]}, Mv[bn]^\bar{Mv[bn-1]})
			   					      = (invSd[Mv[bn]][Mv[bn-1]],invSu[Mv[bn-1]][Mv[bn]])
								   
			   and T is similar to S.
			*/
			int Mv[num_strings];			
			for (int i=0; i< num_strings; i++)
				Mv[i] = v[i];
		
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "num_fixed_points:  initial Mv = ";
	for (int i=0; i< num_strings; i++)
		debug << Mv[i] << " ";
	debug << endl;
}
//output << "v = ";
//for (int i=0; i< num_strings; i++)
//	output << Mv[i] << " ";


			for (int i=0; i< num_terms; i++)
			{

				if (virtual_crossing[i])
				{
					int bn = abs(braid_num[i]);
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "num_fixed_points:    braid_number = " << bn << " Mv[" << bn << "] = " << Mv[bn] 
	      << " Mv[" << bn-1 << "] = " << Mv[bn-1] << endl;
}
					int temp = Td[Mv[bn-1]][Mv[bn]];
					Mv[bn] = Tu[Mv[bn]][Mv[bn-1]];
					Mv[bn-1] = temp;
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "num_fixed_points:    after applying T_" << bn << " Mv[" << bn << "] = " << Mv[bn] 
	      << " Mv[" << bn-1 << "] = " << Mv[bn-1] << endl;
}
				}
				else if (braid_num[i] > 0)
				{
					int bn = braid_num[i];
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "num_fixed_points:    braid_number = " << bn << " Mv[" << bn << "] = " << Mv[bn] 
	      << " Mv[" << bn-1 << "] = " << Mv[bn-1] << endl;
}
					int temp = Sd[Mv[bn-1]][Mv[bn]];
					Mv[bn] = Su[Mv[bn]][Mv[bn-1]];
					Mv[bn-1] = temp;
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "num_fixed_points:    after applying S_" << bn << " Mv[" << bn << "] = " << Mv[bn] 
	      << " Mv[" << bn-1 << "] = " << Mv[bn-1] << endl;
}
				}
				else // (braid_num[i] < 0)
				{
					int bn = abs(braid_num[i]);
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "num_fixed_points:    braid_number = " << bn << " Mv[" << bn << "] = " << Mv[bn] 
	      << " Mv[" << bn-1 << "] = " << Mv[bn-1] << endl;
}
					int temp = invSu[Mv[bn-1]][Mv[bn]];
					Mv[bn] = invSd[Mv[bn]][Mv[bn-1]];
					Mv[bn-1] = temp;
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "num_fixed_points:    after applying invS_" << bn << " Mv[" << bn << "] = " << Mv[bn] 
	      << " Mv[" << bn-1 << "] = " << Mv[bn-1] << endl;
}
				}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "num_fixed_points:  i=" << i << " Mv = ";
	for (int j=0; j< num_strings; j++)
		debug << Mv[j] << " ";
	debug << endl;
}
			}

//output << " Mv = ";
//for (int i=0; i< num_strings; i++)
//	output << Mv[i] << " ";
			
			/* test if it's a fixed point */
			bool fixed_point = true;
			
			for (int i=0; i< num_strings; i++)
			{
				if (Mv[i] != v[i])
				{
					fixed_point = false;
					break;
				}
			}

			if (fixed_point)
			{					
				fixed_points++;
				
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "num_fixed_points:  fixed point" << endl;
	
//output << "invariant" << endl;

			}
			else
			{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "num_fixed_points:  not fixed point" << endl;
	
//output << endl;
			}
			
			/* increment v as a num-stings digit number base n */
			int digit;
			for (digit = num_strings-1; digit >=0; digit--)
			{
				if (++v[digit] == n)
					v[digit] = 0;
				else
					break;
			}
			if (digit == -1)
				not_finished = false;

		} while (not_finished);		
		
		delete [] str;		
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "num_fixed_points: number of fixed points = " << fixed_points << endl;
	
		return fixed_points;
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "num_fixed_points:  Error! returning -1" << endl;
	return -1; // error
}

void rack_poly_invariant(matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, matrix<int>& Tu, 
						 matrix<int>& Td, ST_pair_type pair_type, string input_string, string title)
{
	int num_strings;
	int num_terms;

	if (valid_braid_input(input_string, num_terms, num_strings, braid_control::SILENT_OPERATION, braid_control::RAW_OUTPUT, braid_control::OUTPUT_AS_INPUT))
	{
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

		/* first check that the S and T we've been given are appropriate for the input braid.  The braid
    	   qualifier will be {welded} if it is to be considered as a welded braid, in which case we need
    	   an essential welded pair.  The qualifier will be {doodle} if it is to be considered a doodle braid, 
    	   in which case an essential doodle pair is required.  If neither qualifier is present then an
    	   essential virtual pair is needed.
    	*/
    	bool welded_braid = false;
    	bool doodle_braid = false;
    	bool flat_braid = false;
    	
    	string::size_type pos = input_string.find('{');
    	if (pos != string::npos)
    	{
    		if (input_string.substr(pos).find("welded") != string::npos)
    			welded_braid = true;

    		if (input_string.substr(pos).find("doodle") != string::npos)
    			doodle_braid = true;

    		if (input_string.substr(pos).find("flat") != string::npos)
    			flat_braid = true;
    	}
    
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
   	debug << "rack_poly_invariant: welded_braid =  " << (welded_braid?"true":"false") << endl;
   	debug << "rack_poly_invariant: doodle_braid =  " << (doodle_braid?"true":"false") << endl;
}    

    
		if (doodle_braid)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
   	debug << "rack_poly_invariant: terminating, doodles not supported by rack_poly_invariant" << endl;
   	
			return;
		}
		else if (welded_braid && pair_type != ST_pair_type::ESSENTIAL_WELDED)
    	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
   	debug << "rack_poly_invariant: terminating since S and T are not an essential welded pair, as required for welded braids" << endl;
    
			if (!braid_control::SILENT_OPERATION)
				cout << "\nWelded braid requires essential welded pair, skipping" << endl;
    
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Welded braid requires essential welded pair, skipping" << endl;
			}
    
			return;
    	}
		else if (flat_braid && pair_type != ST_pair_type::FLAT_ESSENTIAL_VIRTUAL)
    	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
   	debug << "rack_poly_invariant: terminating since S and T are not a flat essential virtual pair, as required for flat braids" << endl;
    
			if (!braid_control::SILENT_OPERATION)
				cout << "\nFlat braids require a flat essential virtual pair, skipping" << endl;
    
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Flat braids require a flat essential virtual pair, skipping" << endl;
			}
    
			return;
    	}
		else if (!welded_braid && !flat_braid && pair_type != ST_pair_type::ESSENTIAL_VIRTUAL && pair_type != ST_pair_type::FLAT_ESSENTIAL_VIRTUAL)
    	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
   	debug << "rack_poly_invariant: terminating since S and T are not an essential virtual pair, as required for virtual braids" << endl;
    
			if (!braid_control::SILENT_OPERATION)
				cout << "\ninappropriate essential pair provided for virtual braid, skipping" << endl;
    
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "inappropriate essential pair provided for virtual braid, skipping" << endl;
			}
    
			return;
    	}

		/* Check to see if there's a rack-terms braid qualifier specified 
		   the default is braid_control::RACK_TERMS, which, if specified, overrides any
		   braid qualifier.
		*/
		int rack_terms = 0;
	
		if (pos != string::npos)
		{
	
			string::size_type rpos = input_string.substr(pos).find("rack-terms");
			if (rpos != string::npos)
			{
				get_number(rack_terms, input_string.substr(pos), rpos+11); //expected format is "rack-terms=n"

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "rack_poly_invariant: rack_terms read from braid qualifier, rack-terms =  " << rack_terms << endl;
			}
    			
			input_string = input_string.substr(0,pos);

if (debug_control::DEBUG >= debug_control::DETAIL)
  	debug << "rack_poly_invariant: after removing braid qualifiers, input_string =  " << input_string << endl;

		}
	
		/* if braid_control::RACK_TERMS is non-zero, it's been set as the default and
		   overrides any braid qualifier.  If there's no braid qualifier
	       and braid_control::RACK_TERMS has not been set, we apply a default value of n,
	       the size of the rack.
		*/
		int n = Su.numrows();
		if (braid_control::RACK_TERMS)
		{
			if (!braid_control::SILENT_OPERATION && rack_terms)
				cout << "\noverriding rack-terms in braid qualifier with the value " << rack_terms << " from the command line" << endl;

			rack_terms = braid_control::RACK_TERMS;

if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "rack_poly_invariant: using rack_terms provided from the command line, rack_terms =  " << rack_terms << endl;
		}
		else if (!rack_terms)
		{
			rack_terms = n;

if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "rack_poly_invariant: using default rack_terms value of n = " << rack_terms << endl;
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "rack_poly_invariant: applying " << rack_terms << " rack terms" << endl;
	
		/* The writhe is the number of positive terms in the braid 
		   minus the number of negative terms
		*/
		int writhe = num_terms;
		for (unsigned int i=0; i< input_string.length(); i++)
		{
			if (input_string[i] == '-')
				writhe--;
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "rack_poly_invariant: winding number of input braid =  " << writhe << endl;

		vector<int> rack_poly_coefficients(rack_terms+1);
		ostringstream oss;
    	int fixed_points = num_fixed_points(Su,Sd,invSu,invSd,Tu,Td,input_string);		
    	rack_poly_coefficients[0] = fixed_points;
		oss << fixed_points << "t^" << writhe;
		polynomial<int> rack_poly(oss.str());

if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "rack_poly_invariant: initial term of rack polynomial =  " << rack_poly << endl;
		
		string positive_string = input_string;
		string negative_string = input_string;
		bigint num_vectors = 1;
		for (int i=0; i< num_strings; i++)
			num_vectors *= bigint(n);
			
		for (int i=num_strings; i< num_strings+rack_terms; i++)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "term " << i - num_strings << ", " << num_vectors << " vectors to consider" << endl;
			num_vectors *= bigint(n); // ready for next loop
			
			ostringstream i_num;
			i_num << i;
			positive_string += "s";
			positive_string += i_num.str();
			negative_string += "-s";
			negative_string += i_num.str();

if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "rack_poly_invariant: adding crossing s" << i_num.str() << endl;

			fixed_points = num_fixed_points(Su,Sd,invSu,invSd,Tu,Td,positive_string);
			rack_poly_coefficients[i-num_strings+1] = fixed_points;

if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "rack_poly_invariant:   number of fixed points =  " << fixed_points << endl;

			ostringstream oss_p;
			oss_p << fixed_points << "t^" << writhe + i - num_strings + 1;
			polynomial<int> p_term(oss_p.str());

if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "rack_poly_invariant:   new polynomial term: " << p_term << endl;

			rack_poly += p_term;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
  	debug << "rack_poly_invariant:   rack polynomial updated to : " << rack_poly << endl;
//  	debug << "rack_poly_invariant: adding crossing -s" << i_num.str() << endl;
}

/* take out the calculation for the negative write

			fixed_points = num_fixed_points(Su,Sd,invSu,invSd,Tu,Td,negative_string);

if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "rack_poly_invariant:   number of fixed points =  " << fixed_points << endl;

			ostringstream oss_n;
			oss_n << fixed_points << "t^" << writhe - i + num_strings - 1;
			polynomial<int> n_term(oss_n.str());

if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "rack_poly_invariant:   new polynomial term: " << n_term << endl;
			
			rack_poly += n_term;

if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "rack_poly_invariant:   rack polynomial updated to : " << rack_poly << endl;	
*/	
		}
		
    	if (braid_control::EXTRA_OUTPUT)
    	{
			if (!braid_control::SILENT_OPERATION)
				cout << "\nRack polynomial = " << rack_poly << endl;
			output << "\nRack polynomial = " << rack_poly << endl;
    	}    	
    	else
    	{
			if (!braid_control::SILENT_OPERATION)
			{
				cout <<"\nRack polynomial coefficients = ";
				for (int i=0; i< rack_terms+1; i++)
					cout << rack_poly_coefficients[i] << ' ';
				cout << endl;
			}
		
			if (!braid_control::RAW_OUTPUT)
				output <<"\nRack polynomial coefficients = ";
			
			for (int i=0; i< rack_terms+1; i++)
				output << rack_poly_coefficients[i] << ' ';
			output << endl;			
		}
		
    
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "rack_poly_invariant: final rack_poly = " << rack_poly << endl;
	debug <<"\nRack polynomial = ";
	for (int i=0; i< rack_terms+1; i++)
		debug << rack_poly_coefficients[i] << ' ';
	debug << endl;
}
	}
}

bool distinguish_modified_string (matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, 
						   matrix<int>& Tu, matrix<int>& Td, string input_string, int num_terms)
{
	bool distinguished = false;
	int n = Su.numrows();	
	char* inbuf = c_string(input_string);
	int num_cpts=number_of_components(inbuf,num_terms);
	delete[] inbuf;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "distinguish_modified_string: input_string = " << input_string << endl;
	debug << "distinguish_modified_string: num_terms = " << num_terms << endl;
	debug << "distinguish_modified_string: number of components of resultant braid closure = " << num_cpts << endl;
}

	if ( num_cpts == 1)
	{
		int fixed_points = num_fixed_points(Su,Sd,invSu,invSd,Tu,Td,input_string);
    	   
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "distinguish_modified_string:   number of fixed points = " << fixed_points << endl;
    
		if (fixed_points != n)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "distinguish_modified_string:   distinguished from unknot" << endl;
			distinguished = true;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "distinguish_modified_string:   indistinguishable from unknot" << endl;
		}
	}
	else
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "distinguish_modified_string:   braid closure is a link" << endl;
	}		
	
	return distinguished;
}

