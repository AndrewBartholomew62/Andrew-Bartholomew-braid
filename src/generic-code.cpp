/**************************************************************************
	  This module contains the immersion, peer and gauss code read and 
	  write functions, the generic_code code brokering function and those 
	  functions it calls

void generic_code(string input_string, string title)
void read_immersion_code (generic_code_data& code_data, string input_string)
void write_immersion_code(ostream& s, generic_code_data& code_data)
void print_code_data(generic_code_data& code_data, ostream& s, string prefix)
void read_peer_code (generic_code_data& code_data, string input_string)
void write_peer_code(ostream& s, const generic_code_data& code_data, bool zig_zags, bool labelled)
void cycle_gauss_code_labels(generic_code_data& gauss_code_data, int edge)
string convert_gauss_code(string OU_gauss_code);
void read_gauss_code (generic_code_data& code_data, string input_string)
void write_gauss_code(ostream& s, generic_code_data& code_data)
void read_code_data (generic_code_data& code_data, string input_string)
void write_code_data(ostream& s, generic_code_data& code_data)
bool immersion_to_dowker (matrix<int>& code_table, vector<int>& code)
void affine_index_polynomial(string input_string, generic_code_data& code_data)
void renumber_peer_code(generic_code_data& code_data, vector<int> shift)
bool satellite_code_data(generic_code_data& knot_data, generic_code_data& satellite_data, int strands)
void assign_satellite_code(generic_code_data& satellite_data, int s_edge, int s_prime_edge, int s_component, int s_prime_component, int type, int crossing)
int remove_peer_code_component(generic_code_data& code_data, int component, vector<int>& component_flags)
vector<int> gauss_parity(generic_code_data& code_data)
int remove_virtual_components(generic_code_data& code_data, vector<int>& component_flags)
string minimal_peer_code_representation (generic_code_data peer_code_data)
generic_code_data partition_peer_code(generic_code_data& code_data, vector<int>& component_flags)
bool gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data)
void add_virtual_crossing(int location, vector<int>& fringe, vector<int>& num_virtual_crossings_on_gauss_arc, list<gc_pc_xlabels>& virtual_crossings) 
void assign_gauss_arc_direction(int arc_1, int arc_2, matrix<int>& gauss_code_table, int crossing, vector<int>& gauss_arc_direction)
string direction_to_string(int edge, vector<int>& gauss_arc_direction)
void remove_fringe_edge(int edge, vector<int>& current_fringe)
**************************************************************************/
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <valarray>
#include <list>
#include <iomanip>
#include <ctype.h>
#include <stdio.h>
#include <algorithm>

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


/********************* Function prototypes ***********************/
void bracket_polynomial(generic_code_data& code_data, int variant);
bool immersion_to_dowker (matrix<int>& code_table, vector<int>& code);
bool valid_knotoid_input(generic_code_data& code_data);
void print_code_data(generic_code_data& code_data, ostream& s, string prefix="");
void write_peer_code(ostream& s, const generic_code_data& code_data, bool zig_zags=false, bool labelled=true);
bool gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data);
void write_gauss_code(ostream& s, generic_code_data& code_data);
void read_code_data (generic_code_data& code_data, string input_string);
void write_code_data(ostream& s, generic_code_data& code_data);
void standard_rep (string& word);
string vogel (generic_code_data code_data);
void renumber_peer_code(generic_code_data& code_data, vector<int> shift);
int remove_edge_flags_from_peer_code(generic_code_data& code_data, vector<int>& edge_flag, vector<int>& component_flags);
string affine_index_polynomial(string input_string, generic_code_data& code_data);
bool satellite_code_data(generic_code_data& knot_data, generic_code_data& satellite_data, int strands);
void assign_satellite_code(generic_code_data& satellite_data, int s_edge, int s_prime_edge, int s_component, int s_prime_component, int type, int crossing);
int remove_peer_code_component(generic_code_data& code_data, int component, vector<int>& component_flags);
string parse_long_knot_input_string (string input_string);
vector<int> gauss_parity(generic_code_data& code_data);
list<int> vertex_span (generic_code_data& code_data, int initial_crossing, vector<int>* exclude=0);
void clear_component_flag(int component, vector<int>& component_flags);
void add_virtual_crossing(int location, vector<int>& fringe, vector<int>& num_virtual_crossings_on_gauss_arc, list<gc_pc_xlabels>& virtual_crossings);
void assign_gauss_arc_direction(int arc_1, int arc_2, matrix<int>& gauss_code_table, int crossing, vector<int>& gauss_arc_direction);
string direction_to_string(int edge, vector<int>& gauss_arc_direction);
void remove_fringe_edge(int edge, vector<int>& current_fringe);


/* generic_code is the broker function for an input_string in the form of a labelled code, that is
   either a labelled immersion code or a labelled peer code.  It will do one of the following:
   
    - evaluate link satellites, if required by the input code
    - call the Vogel function
    - call the Dowker code evaluation function (labelled immersion code input only)
    - call the peer code evaluation function
    - call the Gauss code evaluation function
    - call the bracket_polynomial calculation function for the appropriate variant
    - call the affine_index_polynomial function
    - calls renumber_peer_code to renumber or re-orient a diagram described by a labelled peer code
*/
void generic_code(string input_string, string title)
{

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "generic_code: provided with input string: " << input_string << endl;

	/* first isolate any qualifiers from the input string */
	string qualifier;
	unsigned int satellite_code_qualifier = 0;
	
	string::size_type pos = input_string.find('{');
	if (pos != string::npos)
	{
		qualifier = input_string.substr(pos,string::npos);
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "generic_code: isolated qualifier = " << qualifier << endl;
		input_string = input_string.substr(0,pos);
	}

	if (qualifier.find("satellite") != string::npos)
	{
		string::size_type pos = qualifier.find("satellite");
		string::size_type next = qualifier.find(',',pos); //next == string::npos if there's no next qualifier
		string satellite_string = qualifier.substr(pos,next);

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "generic_code: satellite_string: " << satellite_string << endl;

		pos = satellite_string.find('=');


		if (pos != string::npos)
		{
			char* c_satellite_string = c_string(satellite_string);			
			char* cptr = c_satellite_string;	
			while (*cptr != '=')
				cptr++;
				
			get_number(satellite_code_qualifier,++cptr);			
			delete c_satellite_string;
		}
		else
			satellite_code_qualifier = 2;


if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "generic_code: satellite qualifier provided, number of strands = " << satellite_code_qualifier << endl;

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
					
	generic_code_data code_data;
	read_code_data (code_data, input_string);
	
	if (code_data.immersion == generic_code_data::character::LONG_KNOT)
	{	
		input_string = parse_long_knot_input_string(input_string);
		
			
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "generic_code: long knot: " << input_string << endl;
	
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
	debug << "generic_code: Error! long knot indicator provided for the peer code of a link" << endl;

			return;
		}
	}

	if (satellite_code_qualifier)
	{
		if (braid_control::SATELLITE)
		{
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "generic_code: ignoring satellite code qualifier due to braid_control::SATELLITE being specified" << endl;
		}
		else
		{	
			generic_code_data satellite_data;
			
			if (satellite_code_data(code_data, satellite_data, satellite_code_qualifier))
			{
				code_data = satellite_data;

if (braid_control::DEBUG >= braid_control::SUMMARY)
{
	debug << "generic_code: replacing code data with satellite code data" << endl;
	print_code_data(code_data,debug,"generic_code: ");	

}
			}
	    	else
	    	{
				if (!braid_control::SILENT_OPERATION)
					cout << "\nThe satellite function does not support links.\n";
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "The satellite function does not support links.\n";
				}
	    	}
	    }	
	}
	
	if (braid_control::VOGEL_ALGORITHM)
	{
		if(code_data.type == generic_code_data::peer_code)
		{
		    /* call the function vogel to evaluate the braid word.  We check first whether all the
		       crossings are specified as flat crossings, in which case we require the Vogel
		       algorithm to use flat Vogel crossings.  If there are any non-flat crossings in the
		       input code, we require the --flat option to be specified if the intent is for the Vogel
		       algorithm to use flat crossings, otherwise it will use classical Vogel crossings.
		    */
		    bool flat_crossings_only = true;
		    int num_crossings = code_data.num_crossings;
		    matrix<int> code_table = code_data.code_table;
		    
		    for (int i=0; i<num_crossings; i++)
		    {
				if (code_table[LABEL][i] != generic_code_data::FLAT)
				{
					flat_crossings_only = false;
					break;
				}
			}
		    
		    if (flat_crossings_only)
		    {
				braid_control::FLAT_VOGEL_MOVES = true;
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "generic_code: input code contains only flat crossings, setting programme option FLAT_VOGEL_MOVES" << endl;
				
			}
		    
			string braid_word = vogel(code_data);
					
			if (!braid_control::SILENT_OPERATION)
				cout << "\nBraid word = " << braid_word << endl;
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Braid word = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << '\n';
			}
			output << braid_word << endl;
		
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "generic_code: braid word = " << braid_word << endl;
	
			if (braid_word.find("unlink") == string::npos &&
			    braid_word.find('t') == string::npos && 
				braid_word.find('T') == string::npos)
			{
				standard_rep(braid_word);
				if (!braid_control::SILENT_OPERATION)
					cout << "           = " << braid_word << endl;
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? ";" : "           ");
					output << "= " << braid_word << endl;
				}		
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "generic_code:            = " << braid_word << endl;
			}
		}
		else
		{
			cout << "\nError! The Vogel algorithm requires a labelled peer code as input.\n";
			exit(0);
		}
	}
	else if (braid_control::DOWKER_CODE)
	{
		if(code_data.type == generic_code_data::immersion_code)
		{
			vector<int>	dowker_code;
		    /* we use dowker_code to hold the Dowker code once calculated */
	
	    	if (immersion_to_dowker(code_data.code_table,dowker_code))
			{
				int num_terms = dowker_code.size();
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
	    	}
	    	else
	    	{
				if (!braid_control::SILENT_OPERATION)
					cout << "\nDowker code does not support virtual knots.\n";
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Dowker code does not support virtual knots.\n";
				}
	    	}
		}
		else
		{
			cout << "\nError! generic_code presented with a peer code for a call to immersion_to_dowker.\n";
			exit(0);
		}
	}
	else if (braid_control::PEER_CODE)
	{
		if(code_data.type == generic_code_data::gauss_code)
		{
			generic_code_data peer_code_data;
		    if( gauss_to_peer_code(code_data, peer_code_data) )
		    {
				if (!braid_control::SILENT_OPERATION)
				{
					cout << "\npeer code = ";
					write_peer_code(cout,peer_code_data);
				}
				
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "peer code = ";
				}
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
				write_peer_code(output,peer_code_data);
				if (!braid_control::SILENT_OPERATION)
					cout << "\n";
				output << "\n";
			}
		}
		else if(code_data.type == generic_code_data::immersion_code)
		{
			code_data.type = generic_code_data::peer_code;	
			if (!braid_control::SILENT_OPERATION)
			{
				cout << "\npeer code = ";
				write_peer_code(cout,code_data);
			}
			
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "peer code = ";
			}
			if (braid_control::OUTPUT_AS_INPUT)
				output << "\n";
			write_peer_code(output,code_data);
			if (!braid_control::SILENT_OPERATION)
				cout << "\n";
			output << "\n";
		}
		else if (code_data.type == generic_code_data::peer_code && qualifier.find("shift") != string::npos)
		{
			vector<int> shift(code_data.num_components);
			string shift_qualifier = qualifier.substr(qualifier.find("shift"), string::npos);

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "generic_code: shift_qualifier string = " << shift_qualifier << endl;

			char* qualifier_c_string = c_string(shift_qualifier);
			
			char* cptr = qualifier_c_string;
			
			while (*cptr != '[')
				cptr++;
				
			
			for (int i=0; i< code_data.num_components; i++)
			{
				/* step over the initial '[' or separating ',' */
				cptr++;
				
				/* need to accommodate the -0 syntax for component reversal with no shift */
				bool negative = false;
				if (*cptr == '-')
				{
					negative = true;
					cptr++;
				}	
				
				get_number(shift[i],cptr);
				
				if (!shift[i] && negative)
				{
					shift[i] = -1 * code_data.num_component_edges[i];
				}
				else if (negative)
				{
					shift[i] *=-1;
				}
				
				while (*cptr != ',' && *cptr != ']')
					cptr++;
			}
			
			delete[] qualifier_c_string;
			
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "generic_code: renumbering shift vector = ";
	for (int i=0; i< code_data.num_components; i++)
		debug << shift[i] << ' ';
	debug << endl;
}


			/* check the shift vector is valid */
			bool valid_shift_vector = true;
			int root_parity = (shift[0]%2? -1 : 1);
			if (shift[0] < 0)
				root_parity *= -1;
			
			for (int i=1; i< code_data.num_components; i++)
			{
				int parity = (shift[i]%2? -1 : 1);
				if (shift[i] < 0)
					parity *= -1;
				
				if (parity != root_parity)
				{
					valid_shift_vector = false;
					break;
				}
			}

			if (valid_shift_vector)
			{
				renumber_peer_code(code_data, shift);
				if (!braid_control::SILENT_OPERATION)
				{
					cout << "\npeer code = ";
					write_peer_code(cout,code_data);
				}
				
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "peer code = ";
				}
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
				write_peer_code(output,code_data);
				if (!braid_control::SILENT_OPERATION)
					cout << "\n";
				output << "\n";
			}
			else
			{
				cout << "\nError! Invalid shift vector provided in input with PEER_CODE option.\n";
				exit(0);
			}
		}
		else if (code_data.type == generic_code_data::peer_code && braid_control::REMOVE_PEER_CODE_COMPONENT)
		{
			if ((0 <= braid_control::REMOVE_COMPONENT) && (braid_control::REMOVE_COMPONENT < code_data.num_components))
			{
				vector<int> dummy_flags;  // not tracking any components being removed
				remove_peer_code_component(code_data, braid_control::REMOVE_COMPONENT,dummy_flags);

				if (!braid_control::SILENT_OPERATION)
				{
					cout << "\npeer code = ";
					write_peer_code(cout,code_data);
				}
				
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "peer code = ";
				}
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
				write_peer_code(output,code_data);
				if (!braid_control::SILENT_OPERATION)
					cout << "\n";
				output << "\n";
			}
			else
			{
				cout << "\nError! Invalid component number provided in input with REMOVE_PEER_CODE_COMPONENT option.\n";
				exit(0);
			}
		}
		else
		{
			cout << "\nError! generic_code requires a Gauss code, an immersion code. or a shifted peer code when PEER_CODE is true.\n";
			exit(0);
		}
	}
	else if (braid_control::GAUSS_CODE)
	{
		if (code_data.immersion == generic_code_data::character::PURE_KNOTOID && code_data.head != -1)
		{
			if (!valid_knotoid_input(code_data))
			{
				cout << "\nError! Gauss code task presented with an invalid knotoid.\n";
				exit(0);
			}
			
		}
		
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
			write_gauss_code(oss, code_data);
		
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

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "\nbraid::generic_code: Gauss code = " << oss.str() << endl;
		
	}
	else if (braid_control::KAUFFMAN_BRACKET)
	{
		bracket_polynomial(code_data,KAUFFMAN_VARIANT);
	}
	else if (braid_control::JONES_POLYNOMIAL)
	{
		bracket_polynomial(code_data,JONES_VARIANT);
	}
	else if (braid_control::KNOTOID_BRACKET)
	{
		bracket_polynomial(code_data,TURAEV_VARIANT);
	}
	else if (braid_control::ARROW_POLYNOMIAL)
	{
		bracket_polynomial(code_data,ARROW_VARIANT);
	}
	else if (braid_control::PARITY_BRACKET)
	{
		bracket_polynomial(code_data,PARITY_VARIANT);
	}
	else if (braid_control::PARITY_ARROW)
	{
		bracket_polynomial(code_data,PARITY_ARROW_VARIANT);
	}
	else if (braid_control::AFFINE_INDEX)
	{
		string poly = affine_index_polynomial(input_string, code_data);
				
		if (!braid_control::SILENT_OPERATION)
			cout << "\nAffine index polynomial = " << poly << endl;
		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << "Affine index polynomial = ";
			if (braid_control::OUTPUT_AS_INPUT)
				output << '\n';
		}
		output << poly << endl;
		
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "generic_code: Affine index polynomial = " << poly << endl;
	}
	else if (braid_control::SATELLITE)
	{
		generic_code_data satellite_data;
		
		if (satellite_code_data(code_data, satellite_data, braid_control::SATELLITE))
		{
			if (!braid_control::SILENT_OPERATION)
			{
				cout << "\npeer code = ";
				write_peer_code(cout,satellite_data);
			}
			
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "peer code = ";
			}
			if (braid_control::OUTPUT_AS_INPUT)
				output << "\n";
			write_peer_code(output,satellite_data);
			if (!braid_control::SILENT_OPERATION)
				cout << "\n";
			output << "\n";
		}
    	else
    	{
			if (!braid_control::SILENT_OPERATION)
				cout << "\nThe satellite function does not support links.\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "The satellite function does not support links.\n";
			}
    	}	
	}
	else
	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\nUnknown or missing task for label peer code.\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << "Unknown or missing task for label peer code.\n";
		}
		exit(0);
	}
	
} // End of function generic_code

void read_immersion_code (generic_code_data& code_data, string input_string)
{
	code_data.type = generic_code_data::immersion_code;
	code_data.head = -1;

	char* inbuf = c_string(input_string);
	int num_crossings = 0;
	char* cptr = strchr(inbuf,'/');
	cptr++;
	while (*cptr != '\0')
	{
		if (*cptr == '+' || *cptr == '-' || *cptr == '*')
			num_crossings++;

		cptr++;
	}	

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "read_immersion_code: num_crossings determined from labelled immersion code = " << num_crossings << endl;

	int num_edges = 2*num_crossings;
	code_data.num_crossings = num_crossings;
	code_data.num_components = 1;
	vector<int> num_component_edges(1);
	vector<int> first_edge_on_component(1);
	first_edge_on_component[0]=0;
	num_component_edges[0] = num_edges;
	
	matrix<int> code_table(CODE_TABLE_SIZE,num_crossings);
	
	for (int i=0; i< num_crossings; i++)
		code_table[COMPONENT][i] = 0;

	vector<int> perm(num_crossings);
	vector<int> type(num_crossings);
	
	cptr = strchr(inbuf,'(');
    while(*cptr != '/')
	{
		int count = 0;

		/* We write the crossing numbers in the order they appear in the immersion code 
		   into perm together with their corresponding type.  Then we unravel 
		   the permuatation writing the PERM and TYPE rows from INVERSE and LABEL.
		*/
		
		for (int i=0; i<num_crossings; i++)
			type[i] = generic_code_data::TYPE2;

		while (*cptr != ')')
		{		
			if (*cptr == '-')
			{
				type[count] = generic_code_data::TYPE1;
				cptr++;
			}
			else if (isdigit(*cptr))
			{
				char* mark = cptr;
				while (isdigit(*cptr)) cptr++;

				get_number (perm[count],mark);
				
				if (*cptr == '^')
				{
					/* we've found the knotoid head, note the 
					   corresponding crossing number */
					code_data.head = perm[count];
					code_data.immersion = generic_code_data::character::PURE_KNOTOID;
					cptr++;
				}
				
				count++;
			}
			else
				cptr++;
		}

		
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "read_immersion_code: number of terms in this part of perm = " << count << endl;
	debug << "read_immersion_code: read absolute terms: ";
	for (int i=0; i< count; i++)
		debug << perm[i] << ' ';
	debug << endl;		
	debug << "read_immersion_code: signs: ";
	for (int i=0; i< count; i++)
		debug << type[i] << ' ';
	debug << endl;
}
		
		/* work out this part of the permutation */
		for (int i=0; i<count-1; i++)
		{
			code_table[OPEER][perm[i]] = (2*perm[i+1]-1+num_edges)%num_edges;
			code_table[EPEER][(code_table[OPEER][perm[i]]-1)/2] = 2*perm[i];
			code_table[TYPE][perm[i]] = type[i];
		}
		code_table[OPEER][perm[count-1]] = (2*perm[0]-1+num_edges)%num_edges;
		code_table[EPEER][(code_table[OPEER][perm[count-1]]-1)/2] = 2*perm[count-1];
		code_table[TYPE][perm[count-1]] = type[count-1];

		cptr++; // move over ')'
		while (*cptr == ' ') cptr++;

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "read_immersion_code: interim code table produced from immersion code:" << endl;
	debug << "read_immersion_code:   type:    ";	
	for (int i=0; i< num_crossings; i++)
		debug << code_table[TYPE][i] << ' ';
	debug << endl;
	debug << "read_immersion_code:   OPEER: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_table[OPEER][i] << ' ';
	debug << endl;
	debug << "read_immersion_code:   EPEER: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_table[EPEER][i] << ' ';
	debug << endl;
	debug << "read_immersion_code:   label: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_table[LABEL][i] << ' ';
	debug << endl;
}

	}

	/* now do the labels */
	int count = 0;
	while (*cptr != '\0')
	{
		if (*cptr == '+')
			code_table[LABEL][count++] = generic_code_data::POSITIVE;
		else if (*cptr == '-')
			code_table[LABEL][count++] = generic_code_data::NEGATIVE;
		else if (*cptr == '*')
			code_table[LABEL][count++] = generic_code_data::VIRTUAL;
		cptr++;	
	}
	
	/* write the inverse permutation */
//	for (int i=0; i< num_crossings; i++)
//		code_table[INVERSE][code_table[PERM][i]] = i;

	/* now we can write the originating and terminating vertices and edges */
	vector<int> term_crossing(num_edges);
	vector<int> orig_crossing(num_edges);

/*	
	for (int i=0; i< num_crossings; i++)
	{
		term_crossing[2*i] = i;
		orig_crossing[2*i+1] = i;
		term_crossing[(2*code_table[PERM][i]-1 + num_edges) % num_edges] = i;
		orig_crossing[(2*code_table[PERM][i]) % num_edges] = i;
		
		code_table[EVEN_TERMINATING][i] = 2*i;
		code_table[ODD_ORIGINATING][i] = 2*i+1;
		code_table[ODD_TERMINATING][i] = (2*code_table[PERM][i]-1 + num_edges) % num_edges;
		code_table[EVEN_ORIGINATING][i] = (2*code_table[PERM][i]) % num_edges;
	}
*/
	for (int i=0; i< num_crossings; i++)
	{
		term_crossing[2*i] = i;
		orig_crossing[2*i+1] = i;
		term_crossing[code_table[OPEER][i]] = i;
		
		/* we need to identify the edge following the naming edge's peer */
		int component = code_table[COMPONENT][(code_table[OPEER][i]-1)/2];
		int peer_successor = (code_table[OPEER][i]+1 - first_edge_on_component[component])%
		                     num_component_edges[component] + first_edge_on_component[component];

		orig_crossing[peer_successor] = i;
		
		code_table[EVEN_TERMINATING][i] = 2*i;
		code_table[ODD_ORIGINATING][i] = 2*i+1;
		code_table[ODD_TERMINATING][i] = code_table[OPEER][i];
		code_table[EVEN_ORIGINATING][i] = peer_successor;
	}
	
	code_data.code_table = code_table;
	code_data.num_component_edges = num_component_edges;
	code_data.first_edge_on_component = first_edge_on_component;
	code_data.term_crossing = term_crossing;
	code_data.orig_crossing = orig_crossing;
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "read_immersion_code: code data produced from immersion code:" << endl;
	print_code_data(code_data,debug,"read_immersion_code: ");	
}

	delete[] inbuf;
}

/* This function only uses head, num_crossings and the OPEER, TYPE, and LABEL rows 
   of the code_table in the generic_code_data structure.
*/
void write_immersion_code(ostream& s, generic_code_data& code_data)
{
	matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;
	
	vector<int> flag(num_crossings); // used to track which crossings have been written
	for (int i=0; i<num_crossings; i++)
		flag[i] = 1;
		
	/* work out the permutation from the generic code data */
	vector<int> perm(num_crossings);
	for (int i=0; i< num_crossings; i++)
		perm[i] = ((code_table[OPEER][i]+1)/2)%num_crossings;
	
	bool found;
	do
	{
		found = false;
		int crossing;
		int start;
		
		/* look for another starting place in rptr */
		for (int i=0; i<num_crossings; i++)
		{
			if (flag[i])
			{
				crossing = start = i;
				flag[i] = 0;
				found = true;
				break;
			}
		}

		if (found)
		{
			s << "(";
			if (code_table[TYPE][crossing] == generic_code_data::TYPE1)
				s << '-';
			s << crossing;

			int next_crossing;
			do
			{
				next_crossing = perm[crossing];
				s << " ";
				if (code_table[TYPE][next_crossing] == generic_code_data::TYPE1)
					s << '-';
				s << next_crossing;

				if (next_crossing == code_data.head)
					s << '^';
					
				flag[next_crossing] = 0;
				crossing = next_crossing;
			} while (perm[crossing] != start);
			s << ")";
		}
	} while (found);

	s << " / ";
							
	for (int i=0; i< num_crossings; i++)
	{
		if (code_table[LABEL][i] == generic_code_data::POSITIVE)
			s << "+ ";
		if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
			s << "- ";
		if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
			s << "* ";
	}	
}

void print_code_data(generic_code_data& code_data, ostream& s, string prefix)
{
	int num_crossings = code_data.num_crossings;
	
	s << prefix << "code_type = ";
	switch(code_data.type)
	{
		case (generic_code_data::immersion_code): s << "immersion_code"; break;
		case (generic_code_data::peer_code): s << "peer_code"; break;
		case (generic_code_data::gauss_code): s << "gauss_code"; break;
		default: s << "unknown";
	};	
	s << endl;
	
	s << prefix << "immersion (character) = ";
	switch (code_data.immersion)
	{
		case (generic_code_data::character::CLOSED): s << "CLOSED"; break;
		case (generic_code_data::character::LONG_KNOT): s << "LONG_KNOT"; break;
		case (generic_code_data::character::PURE_KNOTOID): s << "PURE_KNOTOID"; break;
		case (generic_code_data::character::KNOTOID): s << "KNOTOID"; break;		
		default: s << "unknown";
	};	
	s << endl;
	
	s << prefix << "head = " << code_data.head << endl;
	s << prefix << "head_zig_zag_count = " << code_data.head_zig_zag_count << endl;
	s << prefix << "num_crossings = " << num_crossings << endl;
	s << prefix << "num_components = " << code_data.num_components << endl;
	s << prefix << "type:      ";	
	matrix<int>& code_table = code_data.code_table;
	
	for (int i=0; i< num_crossings; i++)
		s << code_table[TYPE][i] << ' ';
	s << endl;
	
/*	if (code_data.type == generic_code_data::immersion_code)
	{
		s << prefix << "perm:      ";
		for (int i=0; i< num_crossings; i++)
			s << code_table[PERM][i] << ' ';
		s << endl;
		s << prefix << "inverse:   ";
		for (int i=0; i< num_crossings; i++)
			s << code_table[INVERSE][i] << ' ';
		s << endl;
	}
	else
*/	
	{
		s << prefix << "odd peer: ";
		for (int i=0; i< num_crossings; i++)
			s << code_table[OPEER][i] << ' ';
		s << endl;
		s << prefix << "even peer:  ";
		for (int i=0; i< num_crossings; i++)
			s << code_table[EPEER][i] << ' ';
		s << endl;
	}
	
	s << prefix << "even term: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[EVEN_TERMINATING][i] << ' ';
	s << endl;
	s << prefix << "odd term:  ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[ODD_TERMINATING][i] << ' ';
	s << endl;
	s << prefix << "even orig: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[EVEN_ORIGINATING][i] << ' ';
	s << endl;
	s << prefix << "odd orig:  ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[ODD_ORIGINATING][i] << ' ';
	s << endl;
	s << prefix << "label:     ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[LABEL][i] << ' ';
	s << endl;	
	s << prefix << "component: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[COMPONENT][i] << ' ';
	s << endl;	

//	if (code_data.type == generic_code_data::peer_code)
	{
		s << prefix << "num_component_edges: ";
		for (int i=0; i< code_data.num_components; i++)
			s << code_data.num_component_edges[i] << ' ';
		s << endl;
		s << prefix << "first_edge_on_component: ";
		for (int i=0; i< code_data.num_components; i++)
			s << code_data.first_edge_on_component[i] << ' ';
		s << endl;
	}
	
	s << prefix << "term_crossing: ";
	for (int i=0; i< 2*num_crossings; i++)
		s << code_data.term_crossing[i] << ' ';
	s << endl;
	s << prefix << "orig_crossing: ";
	for (int i=0; i< 2*num_crossings; i++)
		s << code_data.orig_crossing[i] << ' ';
	s << endl;
	s << prefix << "shortcut_crossing: ";
	for (unsigned int i=0; i< code_data.shortcut_crossing.size(); i++)
		s << code_data.shortcut_crossing[i] << ' ';
	s << endl;	
	
	if (code_data.zig_zag_count.numrows() != 0)
	{
		s << prefix << "zig_zag_count: " << endl;
		print (code_data.zig_zag_count, s, 3, prefix);
	}
}

/* This function reads a peer code from the input_string into a generic_code_data object.  Within that
   object it creates a code_table that is a matrix<int>(9,num_crossings).  The function records in this
   matrix the crossing type, the peers of the even and odd edges 
   that terminate at the crossing, the component to which the naming edge belongs and the crossing label. 
   For each crossing it also records the odd and even terminating and originating vertices.
   It uses the #define TYPE, OPEER, EPEER, COMPONENT, LABEL, EVEN_TERMINATING, ODD_TERMINATING, 
   EVEN_ORIGINATING, ODD_ORIGINATING to index the rows of the matrix.
   If non-zero, the int* head is used to identify the first shortcut crossing in a knotoid (identified by 
   a '^' character after the crossing number).  This in turn identifies where in the peer code the head of 
   the knotoid is located. 
*/
void read_peer_code (generic_code_data& code_data, string input_string)
{
if (braid_control::DEBUG >= braid_control::EXHAUSTIVE)
	debug << "read_peer_code: provided with input_string " << input_string << endl;
	
	code_data.type = generic_code_data::peer_code;
	code_data.head = -1;

	if (input_string.find("K:") != string::npos)
		code_data.immersion = generic_code_data::character::KNOTOID;
	else if (input_string.find("L:") != string::npos)
		code_data.immersion = generic_code_data::character::LONG_KNOT;
	else
		code_data.immersion = generic_code_data::character::CLOSED;

	
	char* inbuf = c_string(input_string);
	int num_crossings = 0;
	char* cptr = strchr(inbuf,'/');
	cptr++;
	while (*cptr != '\0')
	{
		if (*cptr == '+' || *cptr == '-' || *cptr == '*' || *cptr == '#' )
			num_crossings++;

		cptr++;
	}	

if (braid_control::DEBUG >= braid_control::EXHAUSTIVE)
	debug << "read_peer_code: num_crossings determined from labelled peer code = " << num_crossings << endl;
	
	code_data.num_crossings = num_crossings;
	
	matrix<int> code_table(CODE_TABLE_SIZE,num_crossings);

	int num_edges = 2*num_crossings;
	int component = 0;

	/* assume all crossings are generic_code_data::TYPE2 and set those 
	   indicated as generic_code_data::TYPE1 accordingly.
	*/
	for (int i=0; i<num_crossings; i++)
		code_table[TYPE][i] = generic_code_data::TYPE2;
		
	cptr = strchr(inbuf,'[');
	cptr++;
    for (int i=0; i< num_crossings; i++)
	{
		/* We write the odd-numbered peer of edge 2i in code_table[OPEER][i] and the 
		   even-numbered peer of edge 2i+1 in code_table[EPEER][i] for i=0,...,num_crossings-1
		*/
		bool crossing_complete = false;
		do
		{
			if (*cptr == ']')
			{
				cout << "\nError! Not enough peers specified in peer code" << endl;

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "read_peer_code: Error! Not enough peers specified in peer code" << endl;
					
				exit(0);
			}
			if (*cptr == '-')
			{
				code_table[TYPE][i] = generic_code_data::TYPE1;
				cptr++;
			}
			else if (isdigit(*cptr))
			{
				char* mark = cptr;
				while (isdigit(*cptr)) cptr++;

				get_number (code_table[OPEER][i],mark);
				
				if (code_table[OPEER][i]%2 == 0)
				{
					cout << "\nError! Read even number " << code_table[OPEER][i] << " in peer code" << endl;

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "read_peer_code: Error!  Read even number " << code_table[OPEER][i] << " in a peer code" << endl;
					
					exit(0);
				}
				
				code_table[EPEER][(code_table[OPEER][i]-1)/2] = 2*i;
				code_table[COMPONENT][i] = component;
				
				if (*cptr == '^')
				{
					/* we've found the knotoid head, note the 
					corresponding crossing number */
					code_data.head = i;
					code_data.immersion = generic_code_data::character::PURE_KNOTOID;
					cptr++;
				}			
				crossing_complete = true;
			}
			else if (*cptr == ',')
			{
				/* we've come to the end of a component and are about to 
				   read the first crossing of the new component
				*/
				component++;
				cptr++;
			}
			else
				cptr++;		
		} while (!crossing_complete);
	}
	
	/* now do the labels */
	int count = 0;
	while (*cptr != '\0')
	{
		if (*cptr == '+')
			code_table[LABEL][count++] = generic_code_data::POSITIVE;
		else if (*cptr == '-')
			code_table[LABEL][count++] = generic_code_data::NEGATIVE;
		else if (*cptr == '*')
			code_table[LABEL][count++] = generic_code_data::VIRTUAL;
		else if (*cptr == '#')
			code_table[LABEL][count++] = generic_code_data::FLAT;
		cptr++;	
	}

	/* write the component data into code_data */
	int num_components = component+1;
	code_data.num_components = num_components;

	vector<int> num_component_edges(num_components);
	vector<int> first_edge_on_component(num_components);
	first_edge_on_component[0]=0;
	component=0;
	
	for (int i=1; i< num_crossings; i++)
	{
		if (code_table[COMPONENT][i] != code_table[COMPONENT][i-1])
		{
			num_component_edges[component] = 2*i-first_edge_on_component[component];
			component++;
			first_edge_on_component[component] = 2*i;
		}
	}
	
	num_component_edges[component] = num_edges - first_edge_on_component[component];
	
	/* now we can write the originating and terminating vertices and edges */
	vector<int> term_crossing(num_edges);
	vector<int> orig_crossing(num_edges);
	
	for (int i=0; i< num_crossings; i++)
	{
		term_crossing[2*i] = i;
		orig_crossing[2*i+1] = i;
		term_crossing[code_table[OPEER][i]] = i;
		
		/* we need to identify the edge following the naming edge's peer */
		int component = code_table[COMPONENT][(code_table[OPEER][i]-1)/2];
		int peer_successor = (code_table[OPEER][i]+1 - first_edge_on_component[component])%
		                     num_component_edges[component] + first_edge_on_component[component];

if (braid_control::DEBUG >= braid_control::EXHAUSTIVE)
	debug << "read_peer_code: crossing " << i << " peer_successor = " << peer_successor << endl;

		orig_crossing[peer_successor] = i;
		
		code_table[EVEN_TERMINATING][i] = 2*i;
		code_table[ODD_ORIGINATING][i] = 2*i+1;
		code_table[ODD_TERMINATING][i] = code_table[OPEER][i];
		code_table[EVEN_ORIGINATING][i] = peer_successor;
	}
	
	code_data.code_table = code_table;
	code_data.num_component_edges = num_component_edges;
	code_data.first_edge_on_component = first_edge_on_component;
	code_data.term_crossing = term_crossing;
	code_data.orig_crossing = orig_crossing;
	
if (braid_control::DEBUG >= braid_control::EXHAUSTIVE)
{
	debug << "read_peer_code: code data produced from peer code:" << endl;
	print_code_data(code_data,debug,"read_peer_code: ");	
}

	delete[] inbuf;
}

/* This function uses head, num_crossings and the OPEER, TYPE, LABEL 
   and COMPONENT rows of the code_table in the generic_code_data structure.
*/
void write_peer_code(ostream& s, const generic_code_data& code_data, bool zig_zags, bool labelled)
{
	const matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;
	
	if (code_data.immersion == generic_code_data::character::KNOTOID)
		s << "K:";
	else if (code_data.immersion == generic_code_data::character::LONG_KNOT)
		s << "L:";
	
	s << '[';
	
	if (code_table[TYPE][0] == generic_code_data::TYPE1)
		s << '-';
	if (code_data.num_crossings >0)
		s << code_table[OPEER][0];
		
	if (code_data.head == 0)
		s << "^ ";
	else if (num_crossings > 1 && code_table[COMPONENT][0] == code_table[COMPONENT][1])
		s << ' ';
		
	for (int i=1; i<num_crossings; i++)
	{
		if (code_table[COMPONENT][i] != code_table[COMPONENT][i-1])
			s << ", ";
			
		if (code_table[TYPE][i] == generic_code_data::TYPE1)
			s << '-';
		s << code_table[OPEER][i];
		if (code_data.head == i)
			s << '^';
			
		if ( i < num_crossings-1 && code_table[COMPONENT][i] == code_table[COMPONENT][i+1])
			s << ' ';
	}
	s << "]";
	if (labelled)
	{
		s << "/";
		for (int i=0; i< num_crossings; i++)
		{
			if (code_table[LABEL][i] == generic_code_data::POSITIVE)
				s << "+";
			else if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
				s << "-";
			else if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
				s << "*";
			else // (code_table[LABEL][i] == generic_code_data::FLAT)
				s << "#";
			
			if (i< num_crossings-1)
				s << " ";
		}	
	}
	
	if (zig_zags && code_data.zig_zag_count.numrows() != 0)
	{
		if (code_data.immersion == generic_code_data::character::KNOTOID || code_data.immersion == generic_code_data::character::LONG_KNOT)
		{
			s << " (";
			
			const matrix<int>& m = code_data.zig_zag_count;
			
			s << m[0][0] << ',' << code_data.head_zig_zag_count << ' ';
			
	        for (size_t j = 1; j < m.numcols(); j++)
	        {
				s << m[0][j] << ' ';
			}
			
		    for (size_t i = 1; i< m.numrows(); i++ )
		    {
		        for (size_t j = 0; j < m.numcols(); j++)
		        {
					s << m[i][j];
					if (i != m.numrows()-1 || j != m.numcols()-1)
						s << ' ';
				}
		    }
			s << ")";
		}
		else
		{
			bool save = matrix_control::SINGLE_LINE_OUTPUT;
			matrix_control::SINGLE_LINE_OUTPUT = true;
			s << " (" << code_data.zig_zag_count << ")";
			matrix_control::SINGLE_LINE_OUTPUT = save;
		}
	}
}

/* cycle_gauss_code_labels shifts the start of the numbering of the component containing edge back by one, 
   so that the number of each edge in the component is incremented (around the component)   
*/
void cycle_gauss_code_labels(generic_code_data& gauss_code_data, int edge)
{
	int num_crossings = gauss_code_data.num_crossings;
	int num_components = gauss_code_data.num_components;
	matrix<int>& gauss_table = gauss_code_data.code_table;
	vector<int>& num_component_edges = gauss_code_data.num_component_edges;
	vector<int>& first_edge_on_component = gauss_code_data.first_edge_on_component;

if (braid_control::DEBUG >= braid_control::EXHAUSTIVE)
{
	debug << "cycle_gauss_code_labels: initial gauss_code_data :" << endl;
	print_code_data(gauss_code_data,debug,"cycle_gauss_code_labels: ");	
	debug << "cycle_gauss_code_labels: cycling component containing edge " << edge << endl;
}

	/* identify the component containing edge */
	int component = 0;
	for (int i=1; i< num_components; i++)
	{
		if (first_edge_on_component[i] > edge)
			break;
		else
			component++;
	}

	int first = first_edge_on_component[component];
	int last = first_edge_on_component[component] + num_component_edges[component]-1;

if (braid_control::DEBUG >= braid_control::EXHAUSTIVE)
	debug << "cycle_gauss_code_labels: edge lies on component " << component << ", first edge on component = " << first << ", last edge = " << last << endl;

    /* we record the new labels in opeer and epeer */
    vector<int> opeer(num_crossings);
    vector<int> epeer(num_crossings);
    
    for (int i=0; i< num_crossings; i++)
    {
		opeer[i] = gauss_table[OPEER][i];
		epeer[i] = gauss_table[EPEER][i];
	}

	/* cycle the labels in opeer and epeer belonging to the identified component */
	for (int i=first; i <= last; i++)
	{
		/* find edge i in the code table */
		for (int j=0; j< num_crossings; j++)
		{
			if (gauss_table[EPEER][j] == i)
			{
				if (i == last)
					epeer[j] = first;
				else
					epeer[j]++;
					
				gauss_table[ODD_ORIGINATING][j] = (epeer[j]+1-first)%num_component_edges[component] + first;
				
				break;
			}
			else if (gauss_table[OPEER][j] == i)
			{
				if (i == last)
					opeer[j] = first;
				else
					opeer[j]++;
					
				gauss_table[EVEN_ORIGINATING][j] = (epeer[j]+1-first)%num_component_edges[component] + first;

				break;
			}
		}
	}

    for (int i=0; i< num_crossings; i++)
    {
		gauss_table[OPEER][i] = opeer[i];
		gauss_table[ODD_TERMINATING][i] = opeer[i];
		gauss_table[EPEER][i] = epeer[i];
		gauss_table[EVEN_TERMINATING][i] = epeer[i];
	}


if (braid_control::DEBUG >= braid_control::EXHAUSTIVE)
{
	debug << "cycle_gauss_code_labels: resultant gauss_code_data :" << endl;
	print_code_data(gauss_code_data,debug,"cycle_gauss_code_labels: ");	
}
	
}

string convert_gauss_code(string OU_gauss_code)
{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "convert_gauss_code: presented with OU_gauss_code: " << OU_gauss_code << endl;			

	ostringstream oss;
	
	/* count the number of crossings */
	int num_crossings = count(OU_gauss_code.begin(),OU_gauss_code.end(),'O');
	vector<int> sign(num_crossings);
	
	char* inbuf = c_string(OU_gauss_code);
	char* cptr = inbuf;
	bool over_arc;
	
	for (int i=0; i< 2*num_crossings; i++)
	{
		while (*cptr != 'O' && *cptr != 'U')
		{
			if (*cptr == ',')
				oss << ',';
				
			cptr++;
		}
		
		if (*cptr == 'O')
			over_arc = true;
		else
			over_arc = false;

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "convert_gauss_code: crossing " << "cptr = " << *cptr << " over_arc = " << (over_arc?"true":"false");			
		
		cptr++;
		
		int crossing_num;
		get_number(crossing_num, cptr);
		
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << ", crossing_num = " << crossing_num;			
		
		while (isdigit(*cptr))
			cptr++;
		
		if (*cptr == '-')
			sign[crossing_num-1] = -1;
		else
			sign[crossing_num-1] = 1;

		cptr++;
		
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << ", sign = " << sign[crossing_num-1] << endl;			

		if (!over_arc)
			oss << '-';
		oss << crossing_num << ' ';
	}
	
	oss << "/";
	for (int i=0; i< num_crossings; i++)
	{
		if (sign[i] == 1)
			oss << '+';
		else
			oss << '-';
	}
			
	delete[] inbuf;
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "convert_gauss_code: returning gauss_code: " << oss.str() << endl;			
	
	return oss.str();	
}


/* A Gauss code does not describe an immersion with crossing labels in the same way that a labelled peer code
   or a labelled immersion code does, since it does not include any information at all about virtual crossings.
   For this reason, if we number the semi-arcs that are described by a Gauss code then in the case of a virtual
   knot we do not always have an odd and even edge arriving at each crossing.  An example of this is the virtual 
   trefoil.  Therefore our use of generic code data for Gauss codes is different to that for peer codes and 
   immersion codes, and the conventions used here have been designed so that the bracket polynomial may be calculated
   for a virtual gauss code without changing the code developed for peer and immersion codes.
   
   A Gauss code for a classical link takes the form "1 -2 5 -4 3 -5, -1 2 -3 4 / ++--+", a Gauss code for a flat link,
   or a doodle has the form "L1 L2 R1 R2 L3 L4 R3 R4 / # # # #" (we distinguish this from the "Gauss orientation data"
   "L1 L2 R1 R2 L3 L4 R3 R4").  Hybrid Gauss codes such as "L1 2 R1 -2 3 -4 -3 4/#+--" are supported but are not used
   currently.  Virtual crossings are explicitly prohibited in Gauss codes, so it is not possible to describe the
   immersion of a virtual(-anything) diagram using the syntax of a Gauss code.  Labelled peer codes should be used
   for such diagrams.
   
   We trace the diagram (well, the classical crossings) described by the Gauss code, numbering semi-arcs without
   regard for any virtual crossings that may be required to realize the code.  We record the edge on which we arrive at 
   a crossing for the first time as the odd peer of the crossing and the edge on which we arrive at the crossing for 
   the second time as the even peer.
   
   The type of crossing is assigned based on the relative orientation of the first and second visits as follows:
           
    1st \ /              2nd \ /         
         X  =  Type 1         X  = Type 2    
    2nd / \              1st / \
   
   Where the orientation of all the semi-arcs is left to right.
   
   We convert the usual Gauss code label convention of identifying positive or negative crossings into the same
   convention used for immersion and peer codes, so that a crossing is labelled '+' if the "even" edge, that is the 
   edge on which we arrive at the crossing for the second time is the over-arc of the crossing, and '-' otherwise.  This
   allows the determination of crossing sign for generic code data to be unaffected.
      
   The odd and even originating and terminating edges are similarly assigned according to the first or second visit to the
   crossing, using the following convention:

    1st OT \ / OO             2nd ET \ / EO         
            X     =  Type 1           X  = Type 2    
    2nd ET / \ EO             1st OT / \ OO
   
   This allows component tracing as used by the bracket polynomial calculation to be unaffected.

   The component recorded for each crossing is the component associated with the edge on which we first arrive
   at the crossing.

*/
void read_gauss_code (generic_code_data& code_data, string input_string)
{
	if (input_string.find('O') != string::npos)
	{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "read_gauss_code: O detected in input string, converting Gauss code format." << endl;
		input_string = convert_gauss_code(input_string);
	}
	
	if (input_string.find('*') != string::npos)
	{
	    cout << "Error! Gauss codes may not describe crossings as virtual" << endl;
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "read_gauss_code: error in input string, '*' character found, indicating virtual crossing." << endl;
	
	    exit(0);
	}
	
	if (input_string.find("K:") != string::npos)
		code_data.immersion = generic_code_data::character::KNOTOID;
	else if (input_string.find("L:") != string::npos)
		code_data.immersion = generic_code_data::character::LONG_KNOT;
	else
		code_data.immersion = generic_code_data::character::CLOSED;
	
	code_data.type = generic_code_data::gauss_code;
	code_data.head = -1;
//	code_data.immersion = generic_code_data::character::CLOSED;
	
	int num_components = (int) count(input_string.begin(),input_string.end(),',') + 1;

	char* inbuf = c_string(input_string);
	int num_crossings = 0;
	char* cptr = strchr(inbuf,'/');
	cptr++;
	while (*cptr != '\0')
	{
		if (*cptr == '+' || *cptr == '-' || *cptr == '#' )
			num_crossings++;

		cptr++;
	}	

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "read_gauss_code: num_crossings determined from gauss code " << input_string << " is " << num_crossings << endl;
	
	code_data.num_crossings = num_crossings;
	
	matrix<int> code_table(CODE_TABLE_SIZE,num_crossings);
	
	int num_edges = 2*num_crossings;
	int component = 0;

	/* 	write the OPEER row of the code table as -1 to act as flags	*/
	for (int i=0; i<num_crossings; i++)
		code_table[OPEER][i] = -1;
		
	cptr = inbuf;
	int edge = 0;
	vector<int> gauss_code_crossing(num_edges);
	vector<int> term_crossing(num_edges);
	vector<int> orig_crossing(num_edges);
	vector<int> num_component_edges(num_components);
	vector<int> first_edge_on_component(num_components);
	
	first_edge_on_component[0] = 0;
	
    for (int i=0; i< 2 * num_crossings; i++)
	{
		/* We write the edge on which we arrive at crossing i for the first time in code_table[OPEER][i] and
		   the edge on which we arrive at the crossing for the second time in code_table[EPEER][i].
		*/
		bool crossing_complete = false;
		int arc_type = generic_code_data::VIRTUAL; // initialized to this as we never use this value with Gauss codes
		do
		{
			if (*cptr == '-')
			{
				arc_type = generic_code_data::NEGATIVE;
				cptr++;
			}
			else if (*cptr == '+')
			{
				arc_type = generic_code_data::POSITIVE;
				cptr++;
			}
			else if (*cptr == 'L')
			{
				arc_type = generic_code_data::LEFT;
				cptr++;
			}
			else if (*cptr == 'R')
			{
				arc_type = generic_code_data::RIGHT;
				cptr++;
			}
			else if (isdigit(*cptr))
			{
				if (arc_type == generic_code_data::VIRTUAL) // i.e. not yet set
				    arc_type = generic_code_data::POSITIVE;
			
				char* mark = cptr;
				while (isdigit(*cptr)) cptr++;
				int crossing;
				get_number (crossing,mark);
				
				/* Gauss codes number crossings from one, we number them from zero */
				crossing--;
				
				gauss_code_crossing[i] = crossing;
				
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "read_gauss_code: crossings = " << crossing << ", edge = " << edge << endl;
	
				if (code_table[OPEER][crossing] == -1)
				{
					/* this is our first visit to this crossing */
					
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "read_gauss_code:   first visit, arc_type = ";
	switch (arc_type)
	{
		case generic_code_data::POSITIVE: 
			debug << "POSITIVE"; 
			break;
		case generic_code_data::NEGATIVE : 
			debug << "NEGATIVE"; 
			break;
		case generic_code_data::LEFT : 
			debug << "LEFT"; 
			break;
		case generic_code_data::RIGHT : 
			debug << "RIGHT";
	} 
	debug << endl;
}
					
					code_table[OPEER][crossing] = edge;						
					code_table[COMPONENT][crossing] = component;
					
					/* If we are at a classical crossing we can assign the crossing label as we now know whether 
					   the second visit (our "even" edge) will arrive on the over-arc or under-arc.  If the crossing 
					   is flat we can use the arc_tytpe to set the TYPE of the crossing as well as the label.
					*/
					if (arc_type == generic_code_data::POSITIVE)
						code_table[LABEL][crossing] = generic_code_data::NEGATIVE;
					else if (arc_type == generic_code_data::NEGATIVE)
						code_table[LABEL][crossing] = generic_code_data::POSITIVE;
					else
					{
						code_table[LABEL][crossing] = generic_code_data::FLAT;
						
						if (arc_type == generic_code_data::LEFT)
							code_table[TYPE][crossing] = generic_code_data::TYPE2;
						else
							code_table[TYPE][crossing] = generic_code_data::TYPE1;
							
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "read_gauss_code:   set type = " << (code_table[TYPE][crossing] == generic_code_data::TYPE1? "TYPE1":"TYPE2") << endl;
					}
				}
				else
				{
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "read_gauss_code:   second visit, arc_type = ";
	switch (arc_type)
	{
		case generic_code_data::POSITIVE: 
			debug << "POSITIVE"; 
			break;
		case generic_code_data::NEGATIVE : 
			debug << "NEGATIVE"; 
			break;
		case generic_code_data::LEFT : 
			debug << "LEFT"; 
			break;
		case generic_code_data::RIGHT : 
			debug << "RIGHT";
	} 
	debug << endl;
}
					/* this is the second visit to this crossing */
					code_table[EPEER][crossing] = edge;
				}
								
				edge++;
				crossing_complete = true;
			}
			else if (*cptr == ',')
			{
				/* we've come to the end of a component and are about to 
				   read the first crossing of the new component
				*/
				num_component_edges[component] = edge - first_edge_on_component[component];
				component++;
				first_edge_on_component[component] = edge;
				cptr++;
			}
			else
				cptr++;		
		} while (!crossing_complete);
	}
	
	/* finish off num_component_edges for the last component */
	num_component_edges[component] = edge - first_edge_on_component[component];

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "read_gauss_code: gauss_code_crossing: ";
	for (int i=0; i< num_edges; i++)
		debug << gauss_code_crossing[i] << ' ';
	debug << endl;
}

	/* now we can identify the TYPE of each clasical crossing from the sign of the 
	   crossing in the Gauss code and the LABEL assignment.  Flat crossing have
	   had their TYPE set already.
	*/
	int count = 0;
	while (*cptr != '\0')
	{
		if (*cptr == '+')
		{
			if (code_table[LABEL][count] == generic_code_data::POSITIVE)
				code_table[TYPE][count] = generic_code_data::TYPE2;
			else
				code_table[TYPE][count] = generic_code_data::TYPE1;
				
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "read_gauss_code:   set type for crossing " << count << ": positive crossing, second visit " << (code_table[LABEL][count] == generic_code_data::POSITIVE? "over": "under")
	      << " type is " << (code_table[TYPE][count] == generic_code_data::TYPE1? "TYPE1":"TYPE2") << endl;
}
			count++;
		}
		else if (*cptr == '-')
		{
			if (code_table[LABEL][count] == generic_code_data::POSITIVE)
				code_table[TYPE][count] = generic_code_data::TYPE1;
			else
				code_table[TYPE][count] = generic_code_data::TYPE2;

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "read_gauss_code:   set type for crossing " << count << ": negative crossing, second visit " << (code_table[LABEL][count] == generic_code_data::POSITIVE? "over": "under")
	      << " type is " << (code_table[TYPE][count] == generic_code_data::TYPE1? "TYPE1":"TYPE2") << endl;
}

			count++;
		}
		else if (*cptr == '#')
		{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "read_gauss_code:   type for flat crossing " << count << " already set" << endl;
		    count++;
		}
		cptr++;	
	}
	
	/* identify the originating and terminating edges at each crossing.  We cannot do this
	   during our initial scan of the Gauss code as we need to know when we have reached the
	   end of each component and we don't know that until after we've processed the last crossing
	   on the component.
	*/
	vector<bool> first_visit_to_crossing(num_crossings);
	for (int i=0; i< num_crossings; i++)
		first_visit_to_crossing[i] = true;
	
	
	edge = 0;
	for (int i=0; i< num_components; i++)
	{
		for (int j=0; j< num_component_edges[i]; j++)
		{
			int crossing = gauss_code_crossing[edge++];
			if (first_visit_to_crossing[crossing])
			{
				code_table[ODD_TERMINATING][crossing] = first_edge_on_component[i] + j;
				code_table[EVEN_ORIGINATING][crossing] = first_edge_on_component[i] + (j+1)%num_component_edges[i];
				first_visit_to_crossing[crossing] = false;
			}
			else
			{				
				code_table[EVEN_TERMINATING][crossing] = first_edge_on_component[i] + j;
				code_table[ODD_ORIGINATING][crossing] = first_edge_on_component[i] + (j+1)%num_component_edges[i];
			}
			term_crossing[first_edge_on_component[i] + j] = crossing;
			orig_crossing[first_edge_on_component[i] + (j+1)%num_component_edges[i]] = crossing;
		}
	}
	
	/* write the data into code_data */
	code_data.num_components = num_components;	
	code_data.code_table = code_table;
	code_data.num_component_edges = num_component_edges;
	code_data.first_edge_on_component = first_edge_on_component;
	code_data.term_crossing = term_crossing;
	code_data.orig_crossing = orig_crossing;
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "read_gauss_code: code data produced from gauss code:" << endl;
	print_code_data(code_data,debug,"read_gauss_code: ");	
}

	delete[] inbuf;
}

/* write_gauss_code is capable of writing a gauss code from any kind of generic code data, 
   so it can convert from peer codes or immersion codes as well as output "native" Gauss codes
*/
void write_gauss_code(ostream& s, generic_code_data& code_data)
{
	matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;
	int num_components = code_data.first_edge_on_component.size();

	if (code_data.type == generic_code_data::gauss_code)
	{
		/* we are writing out a "native" Gauss code, so the code_data follows Gauss code conventions */
		
		/* we have to follow the edges around the Gauss code noting each crossing
		   we encounter as we do so.  Since Gauss codes record the comonent of the first
		   visit in the COMPONENT row of the code table, we cannot use this row to output
		   the ',' between components accurately, since there may be two or more additional 
		   components in the diagram that are not recorded in that row.  Instead, we have to use
		   the edges recorded in first_edge_on_component to determine when we reach the end of a
		   component.
	        */	  
		for (int i=0; i < 2*num_crossings; i++)
		{
			bool last_edge_on_component = false;
			for (int j=0; j< num_components; j++)
			{
				if (i+1 == code_data.first_edge_on_component[j])
				{
					last_edge_on_component = true;
					break;
				}
			}

			/* locate the edge in the code table */
			for (int j=0; j< num_crossings; j++)
			{
				if (code_table[OPEER][j] == i)
				{
					s << ' ';
	
					if (code_table[LABEL][j] == generic_code_data::POSITIVE)
						s << '-'; // we're going under the crossing
					s << j+1;
					break;
				}
				else if (code_table[EPEER][j] == i)
				{					
					s << ' ';
	
					if (code_table[LABEL][j] == generic_code_data::NEGATIVE)
						s << '-'; // we're going under the crossing
					s << j+1;
					break;
				}

			}
				if (last_edge_on_component)
					s << ',';
				
		}
	
		s << " / ";
		for (int i=0; i< num_crossings; i++)
		{
			if ( (code_table[LABEL][i] == generic_code_data::POSITIVE && code_table[TYPE][i] == generic_code_data::TYPE2) ||
			     (code_table[LABEL][i] == generic_code_data::NEGATIVE && code_table[TYPE][i] == generic_code_data::TYPE1)
			   )
				s << "+ ";
			else if ( (code_table[LABEL][i] == generic_code_data::POSITIVE && code_table[TYPE][i] == generic_code_data::TYPE1) ||
			          (code_table[LABEL][i] == generic_code_data::NEGATIVE && code_table[TYPE][i] == generic_code_data::TYPE2)
					)
				s << "- ";
			else
				s << "! "; // shouldn't get here
		}	
	}
	else
	{
		/* write Gauss code from immersion data in code_data */
		vector<int>& term_crossing = code_data.term_crossing;
		int num_edges = 2*num_crossings;
		vector<int> edge_flag(num_edges); // initialized to zero

		/* we need to re-number the crossings if we have virtual crossings in the immersion, 
		   since these are ignored by the Gauss code.
		*/
		int num_classical_crossings = num_crossings;
		bool pure_knotoid_code_data = false;
		vector<int>& shortcut_crossing = code_data.shortcut_crossing;
		
		for (int i=0; i< num_crossings; i++)
		{
			if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
				num_classical_crossings--;
		}
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "write_gauss_code: num_classical_crossings = " << num_classical_crossings << endl;

		if (code_data.immersion == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && shortcut_crossing.size())
		{
			pure_knotoid_code_data = true;
			for (unsigned int i=0; i< shortcut_crossing.size(); i++)
			{
				if (shortcut_crossing[i] != 0)
					num_classical_crossings--;
			}
if (braid_control::DEBUG >= braid_control::EXHAUSTIVE)
	debug << "write_gauss_code: knotoid: num_classical_crossings = " << num_classical_crossings << endl;
		}
	
		int num_classical_crossings_visited = 0;
		
		/* classical_crossing will be used to renumber the immersion crossings,
		   so that if classical_crossing[i] = j, then crossing j of the immersion
		   is renumbered as crossing i in the Gauss code.
		   
		   crossing_visited will be a flag to indicate whether the ith immersion crossing
		   has been recorded in classical_crossing or not.
		*/
		vector<int> classical_crossing(num_classical_crossings);
		vector<int> crossing_visited(num_crossings);

		int start=0;
		int edge=0;
		bool complete = false;


		if (code_data.immersion == generic_code_data::character::KNOTOID)
			s << "K:";
		else if (code_data.immersion == generic_code_data::character::LONG_KNOT)
			s << "L:";
		
		do 
		{
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "write_gauss_code: component_start = " << start << endl;
			/*	trace this component */
			do
			{	

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "write_gauss_code: edge = " << edge;
				edge_flag[edge] = 1;
				int next_crossing = term_crossing[edge];

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << ", next_crossing = " << next_crossing;
	
				if (code_table[LABEL][next_crossing] == generic_code_data::VIRTUAL)
				{
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << ", is virtual" << endl;
				}
				else if (pure_knotoid_code_data && shortcut_crossing[next_crossing])
				{
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << ", is a shortcut crossing" << endl;
				}
				else if ((edge%2 != 0 && code_table[LABEL][next_crossing] == generic_code_data::POSITIVE) ||
				    (edge%2 == 0 && code_table[LABEL][next_crossing] == generic_code_data::NEGATIVE))
				{
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << ", going under" << endl;
					s << '-';
				}
				else
				{
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << ", going over" << endl;
				}
				
				if (code_table[LABEL][next_crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[next_crossing]))
				{
					if(crossing_visited[next_crossing])
					{
						for (int i=0; i< num_classical_crossings_visited; i++)
						{
							if (classical_crossing[i] == next_crossing)
							{
								s << i+1; // Gauss crossings are numbered from 1 not zero

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "write_gauss_code:   second visit to Gauss crossing " << i+1<< endl;
							}
						}
					}
					else
					{

						classical_crossing[num_classical_crossings_visited] = next_crossing;
						crossing_visited[next_crossing] = 1;
						num_classical_crossings_visited++;
						
						s << num_classical_crossings_visited;

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "write_gauss_code:   first visit, becomes Gauss crossing " << num_classical_crossings_visited << endl;
					}

					if (edge%2)
						edge = code_table[EVEN_ORIGINATING][next_crossing];
					else
						edge = code_table[ODD_ORIGINATING][next_crossing];
						
					if (edge != start)
						s << ' ';

				}
				else
				{
					/* just move on around the component */				
					if (edge%2)
						edge = code_table[EVEN_ORIGINATING][next_crossing];
					else
						edge = code_table[ODD_ORIGINATING][next_crossing];

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "write_gauss_code:   doing nothing" << endl;
				}				
			} while (edge != start);

			/* look for the start of another component */
			complete = true;
			for (int i=0; i< num_edges; i++)
			{
				if (edge_flag[i] == 0)
				{
					complete = false;
					start = i;
					edge = start;
					break;
				}
			}
			
			if (!complete)
				s << ',';
				
		} while (!complete);
		
		s << '/';
		
		for (int i=0; i< num_classical_crossings; i++)
		{
			int crossing = classical_crossing[i];
			
			if (code_table[LABEL][crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[crossing])) // it should never be!
			{
				if ((code_table[TYPE][crossing] == generic_code_data::TYPE1 && code_table[LABEL][crossing] == generic_code_data::NEGATIVE) ||
				    (code_table[TYPE][crossing] == generic_code_data::TYPE2 && code_table[LABEL][crossing] == generic_code_data::POSITIVE))
					s << '+';
				else
					s << '-';
			}
			
			if (i< num_classical_crossings-1)
				s << " ";
			
		}
	}
}

void write_code_data(ostream& s, generic_code_data& code_data)
{
	if (code_data.type == generic_code_data::peer_code)
		write_peer_code(s, code_data);
	else if (code_data.type == generic_code_data::immersion_code)
		write_immersion_code(s,code_data);
	else
		write_gauss_code(s,code_data);
}

void read_code_data (generic_code_data& code_data, string input_string)
{
	if (input_string.find('[') != string::npos)
		read_peer_code(code_data, input_string);
	else if (input_string.find('(') != string::npos)
		read_immersion_code(code_data, input_string);
	else
		read_gauss_code(code_data, input_string);
}

/* the function immersion_to_dowker calculates the dowker code for the
   labelled immersion code in code_table, storing the result at code.  
   
   The function returns a boolean indicating whether a code has been
   written to 'code'.  If the code includes a virtual label the function  
   returns false and leaves code unchanged.  Otherwise the code is set
   to a vector containing the dowker code and the function returns true.

   Given that we have a code table most of the work for the 
   dowker code is done except that the dowker code numbers semi-arcs from 1 rather 
   than zero.  (See comments to the function braid_to_dowker for a 
   description of the dowker code).  Therefore the OPEER row of the 
   table gives the Dowker code once the odd-numbered peers have been incremented
   by one.   
   
   Note also that if we are given generic code data for an 
   alternating knot then all the labels will be '+'.  Those crossings for
   which the label is set to '-' are precisely those that require the 
   corresponding dowker code element to be negated.
*/
bool immersion_to_dowker (matrix<int>& code_table, vector<int>& code)
{
	int num_crossings=code_table.numcols();
	bool classical = true;
	
	vector<int> dowker_code(num_crossings);
	
	for (int i=0; i<num_crossings; i++)
		dowker_code[i] = code_table[OPEER][i]+1;
	
	for (int i=0; i<num_crossings; i++)
	{
		if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
		{
			dowker_code[i] *= -1;
		}
		else if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
		{
			classical = false;
			break;
		}
	}
	
	if (classical)
		code = dowker_code;
		
	return classical;
}

/* The affine_index_polynomial calculated here is that described in Kauffman's 
   paper "An Affine Index Polynomial Invariant of Virtual Knots" (arXiv:1211.1601v1).
   As well as being a non-trivial invariant of virtual knots it is useful to distinguish
   virtual and classical knots, since the affine index polynomial takes the value zero for
   a classical knot.

	The affine index polynomial does not support code data specifying flat crossings.

*/
string affine_index_polynomial(string input_string, generic_code_data& code_data)
{
	if (code_data.num_components != 1)
	{

if (braid_control::DEBUG >= braid_control::SUMMARY)
{
	debug << "affine_index_polynomial: code_data.num_components != 1, links are not supported by this function.";
}
		return "Error!";
	}
		
	int num_crossings = code_data.num_crossings;
	matrix<int>& code_table = code_data.code_table;

	/* check we've not been given flat crossings */
	for (int i=0; i< num_crossings; i++)
	{
		if (code_table[LABEL][i] == generic_code_data::FLAT)
		{

			if (!braid_control::SILENT_OPERATION)
				cout << "\n\naffine_index_polynomial does not support flat crossings, skipping";
	
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "\n\naffine_index_polynomial does not support flat crossings, skipping";
				if (braid_control::OUTPUT_AS_INPUT)
					output << '\n';
			}

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "affine_index_polynomial: input includes flat crossings, doing nothing" << endl;

			return "Error!";
		}
	}
	
	bool pure_knotoid_code_data = false;
	
	if (code_data.immersion == generic_code_data::character::PURE_KNOTOID && code_data.head != -1)
	{
		if (valid_knotoid_input(code_data))
		{
			pure_knotoid_code_data = true;
			
if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "affine_index_polynomial: shortcut_crossings: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_data.shortcut_crossing[i] << ' ';
	debug << endl;
}
		}
	}

	/* evaluate the Cheng colouring of the diagram described by code_data, starting at zero
	   by tracing around the diagram.  Note that we have to deal with peer, immersion and 
	   Gauss codes.
	*/
	int colour=0;
	vector<int> a_colour(num_crossings);
	vector<int> b_colour(num_crossings);
	int edge=0;
	
	for (int i=0; i< 2*num_crossings; i++)
	{
		int next_crossing = code_data.term_crossing[edge];
		
		if (code_table[LABEL][next_crossing] == generic_code_data::VIRTUAL || (pure_knotoid_code_data && code_data.shortcut_crossing[next_crossing]))
		{
			/* do nothing, just move around the diagram */
			edge++;
			continue;
		}				
		else if (code_table[LABEL][next_crossing] == generic_code_data::POSITIVE || code_table[LABEL][next_crossing] == generic_code_data::NEGATIVE)
		{
			bool a_arc = false;
			if ( (edge == code_table[ODD_TERMINATING][next_crossing] && code_table[TYPE][next_crossing] == generic_code_data::TYPE1)
			   ||(edge == code_table[EVEN_TERMINATING][next_crossing] && code_table[TYPE][next_crossing] == generic_code_data::TYPE2)
			   )
			{
				a_arc = true;
			}
			
			if (a_arc)
			{
				a_colour[next_crossing] = colour;
				colour--;
			}
			else
			{
				b_colour[next_crossing] = colour;
				colour++;
			}
			
			edge++;
		}
		else
		{
			cout << "Error! Unknown label type in affine_index_polynomial." << endl;
			exit(0);
		}
	}
	
if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "affine_index_polynomial: a_colour = ";
	for (int i=0; i< num_crossings; i++)
		debug << a_colour[i] << ' ';
	debug << endl;

	debug << "affine_index_polynomial: b_colour = ";
	for (int i=0; i< num_crossings; i++)
		debug << b_colour[i] << ' ';
	debug << endl;
}

	/* Now evaluate the weight W_plus for each classical crossing, we shall use the 
	   fact that W_minus = -W_plus when calulating the polynomil itself.
	*/
	vector<int> W_plus(num_crossings);
	
	for (int i=0; i< num_crossings; i++)
	{
		if (code_table[LABEL][i] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && code_data.shortcut_crossing[i]))
			W_plus[i] = a_colour[i] - b_colour[i] -1;
	}

if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "affine_index_polynomial: W_plus = ";
	for (int i=0; i< num_crossings; i++)
		debug << W_plus[i] << ' ';
	debug << endl;
}

	/* Note the sign of each crossing (positive or negative crossing, or virtual crossing).  
	*/
	vector<int> sign(num_crossings);

	for (int i=0; i< num_crossings; i++)
	{
		if (pure_knotoid_code_data && code_data.shortcut_crossing[i])
		{
			/* regard shortcut crossings as virtual for the purposes of sign, so we ignore them */
			sign[i] = braid_crossing_type::VIRTUAL;
		}
		else if (   (code_table[TYPE][i] == generic_code_data::TYPE1 && code_table[LABEL][i] == generic_code_data::NEGATIVE)
		    || (code_table[TYPE][i] == generic_code_data::TYPE2 && code_table[LABEL][i] == generic_code_data::POSITIVE)
		   )
		{
			/* positive crossing */
			sign[i]  = braid_crossing_type::POSITIVE;
	    }
		else if (   (code_table[TYPE][i] == generic_code_data::TYPE1 && code_table[LABEL][i] == generic_code_data::POSITIVE)
		          || (code_table[TYPE][i] == generic_code_data::TYPE2 && code_table[LABEL][i] == generic_code_data::NEGATIVE)
		        )
		{
			/* negative crossing */
			sign[i] = braid_crossing_type::NEGATIVE;
	    }
	    else
		{
			/* virtual crossing */
			sign[i] = braid_crossing_type::VIRTUAL;
	    }			    		    
	}
	
if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "affine_index_polynomial: crossing sign (positive = 1, negative = -1, virtual = 0): ";
	for (int i=0; i< num_crossings; i++)
		debug << sign[i] << ' ';
	debug << endl;
}

	/* Evaluate the writhe */
	
	int writhe = 0;
	for (int i=0; i< num_crossings; i++)
	{
		if (sign[i] == braid_crossing_type::POSITIVE)
			writhe ++;
		else if (sign[i] == braid_crossing_type::NEGATIVE)
			writhe --;
	}

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "affine_index_polynomial: writhe = " << writhe << endl;
	
	/* write the polynomial to a ostringstream, then read it back again to put it into standard format */
	ostringstream oss;
	
	for (int i=0; i< num_crossings; i++)
	{
//		if (code_table[LABEL][i] != generic_code_data::VIRTUAL)
		if (sign[i] != braid_crossing_type::VIRTUAL)
		{
			if (sign[i] == braid_crossing_type::POSITIVE)
				oss << '+';
			else // (sign[i] == braid_crossing_type::NEGATIVE)
				oss << '-';
				
			oss << "t^";
			
			if (sign[i] == braid_crossing_type::POSITIVE)
				oss << W_plus[i];
			else 
				oss << -W_plus[i];
		}
	}
	
	if (writhe < 0)
		oss << '+' << abs(writhe);
	else if (writhe > 0)
		oss << '-' << writhe;

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "affine_index_polynomial: initial_polynomial = " << oss.str() << endl;
	
	polynomial<int> initial_polynomial(oss.str());
	
	oss.str("");
	oss << initial_polynomial;

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "affine_index_polynomial: returning polynomial " << oss.str() << endl;

	return oss.str();
}

/* renumber_peer_code moves the starting point for the numbering of each component
   forwards with respect to the orientation by the number of semi-arcs given by 
   value of the shift vector corresponding to that component. 
   
   If the shift value is positive it indicates a forwards shift of the starting point for
   the numbering of the corresponding component.  If it is negative, it represents a forwards shift
   of the absolute value together with an orientation reversal.  
   
   If the code represents a link and any component is shifted by an odd number of edges, 
   we must have all intersecting components shifted by an odd number of edges, otherwise 
   we will violate the requirement to have an odd and even edge terminating at every crossing.
   Since the peer code has to be connected to be realizable, this means that every component 
   has to be shifted by an odd number of edges.  Thus we have to shift every component by an odd
   number of edges or every component by an even number of edges.  
  
   If a component's orientation is reversed we also violate the requirement to have an odd and 
   even edge terminating at every crossing.  Therefore an orientation reversal must be accompanied by
   an odd shift of the component being reversed or of those other components it meets at a crossing.	
   By considering orientation reversal an "odd" operation, a particular shift vector is valid if 
   every entry is odd or every entry is even.  Thus -2 is odd and -1 is even. It is the responsibility 
   of the calling code to ensure that the shift vector is valid.    
*/
void renumber_peer_code(generic_code_data& code_data, vector<int> shift)
{
    int head_semi_arc;
    if (code_data.head !=-1)
    {	
		if (code_data.code_table[LABEL][code_data.head] == generic_code_data::POSITIVE)
			head_semi_arc = code_data.code_table[OPEER][code_data.head];
		else if (code_data.code_table[LABEL][code_data.head] == generic_code_data::NEGATIVE)
			head_semi_arc = 2*code_data.head;
		else
		{
			cout << "\nError! Function renumber_peer_code presented with a knotoid code whose first shortcut crossing is virtual." << endl;
			exit(0);
		}
		
		if (shift[0] != 0 && shift[0] != -head_semi_arc)	
		{
if (braid_control::DEBUG >= braid_control::SUMMARY)
		debug << "renumber_peer_code: asked to renumber a knotoid with a shift vector that changes the shortcut, doing nothing" << endl;
		
		   return;
		}
    }
	
	int num_crossings = code_data.num_crossings;
	int num_edges = 2* num_crossings;
	int num_components = code_data.num_components;
	matrix<int>& code_table = code_data.code_table;
	

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
//	debug << "renumber_peer_code: given code data:" << endl;
//	print_code_data(code_data,debug,"renumber_peer_code: ");	
	debug << "renumber_peer_code: shift = ";
	for (int i=0; i< num_components; i++)
		debug << shift[i] << ' ';
	debug << endl;
}

	/* create an edge map that records in edge_map[i] the new number of edge i
	   after the shift
	*/
	vector<int> initial_edge_map(num_edges);
	vector<int>& first_edge_on_component = code_data.first_edge_on_component;
	vector<int>& num_component_edges = code_data.num_component_edges;
	
	for (int i=0; i< num_crossings; i++)
	{
		int component = code_table[COMPONENT][i];
		initial_edge_map[2*i] = (2*i - first_edge_on_component[component] - abs(shift[component]) + num_component_edges[component])%
		                num_component_edges[component] + first_edge_on_component[component];

		initial_edge_map[2*i+1] = (initial_edge_map[2*i] + 1 - first_edge_on_component[component])%
		                num_component_edges[component] + first_edge_on_component[component];		   
	}
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code: initial_edge_map ";
	for (int i=0; i< num_edges; i++)
		debug << setw(3) << i;
	debug << endl;
	debug << "renumber_peer_code:                  ";
	for (int i=0; i< num_edges; i++)
		debug << setw(3) << initial_edge_map[i];
	debug << endl;
}

	/* If any of the shift vector values are negative we have to reverse the numbering
	   of that component
	*/
	vector<int> edge_map(initial_edge_map);
	for (int i=0; i< num_components; i++)
	{
		if (shift[i] < 0)
		{
			for (int j=1; j < num_component_edges[i]; j++)
			{
				int edge = first_edge_on_component[i]+j;
				
				/* find edge in the initial edge_map */
				int position;
				for (int k=0; k< num_edges; k++)
				{
					if (initial_edge_map[k] == edge)
					{
						position = k;
						break;
					}
				}
				
				edge_map[position] = first_edge_on_component[i] + num_component_edges[i] - j;
			}
		}
	}
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code: edge_map ";
	for (int i=0; i< num_edges; i++)
		debug << setw(3) << i;
	debug << endl;
	debug << "renumber_peer_code:          ";
	for (int i=0; i< num_edges; i++)
		debug << setw(3) << edge_map[i];
	debug << endl;
}

	/* re-write the even and odd peers in the code_table, we must determine
	   the new peers for each crossing from the edge_map and then write them
	   into the correct location in the code_table.  Note that the crossing
	   numbering will change as a result of renumbering the edges.
	*/
	int even_edge;
	int odd_edge;
	matrix<int> new_code_table(code_table);
		
	for (int i=0; i< num_crossings; i++)
	{
		int component = code_table[COMPONENT][i];
		int component_edges = num_component_edges[component];
		int first_edge = first_edge_on_component[component];
		int peer_component = code_table[COMPONENT][(code_table[OPEER][i]-1)/2];
		int peer_component_edges = num_component_edges[peer_component];
		int peer_first_edge = first_edge_on_component[peer_component];
		
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code: old crossing " << i << ", component " << component 
	      << ", peer_component " << peer_component << endl; 
}

		if (edge_map[2*i]%2)
		{

			if (shift[component] < 0 && shift[peer_component] < 0)
			{

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " now odd" << endl;
	debug << "renumber_peer_code:   both components at old crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   type and label unchanged for new crossing" << endl;
}

				even_edge = (edge_map[2*i]-first_edge-1+component_edges)%component_edges + first_edge;
				odd_edge = (edge_map[code_table[OPEER][i]]-peer_first_edge-1+peer_component_edges)%peer_component_edges + peer_first_edge;

				/* the new naming edge is on the same strand as the old naming edge */
				new_code_table[LABEL][even_edge/2] = code_table[LABEL][i];

				/* since the old naming edge is now odd, reversing both strands leaves the crossing type unchanged */
				new_code_table[TYPE][even_edge/2] = code_table[TYPE][i];


			}
			else if (shift[component] < 0)
			{

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " now odd" << endl;
	debug << "renumber_peer_code:   naming edge component at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   type reversed but label unchanged for new crossing" << endl;
}

				even_edge = (edge_map[2*i]-first_edge-1+component_edges)%component_edges + first_edge;
				odd_edge = edge_map[code_table[OPEER][i]];

				/* the new naming edge is on the same strand as the old naming edge */
				new_code_table[LABEL][even_edge/2] = code_table[LABEL][i];

				/* since the old naming edge is now odd, reversing the old naming strand reverses crossing type */
				if (code_table[TYPE][i] == generic_code_data::TYPE1)
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE2;
				else
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE1;
			}
			else if (shift[peer_component] < 0)
			{

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " now odd" << endl;
	debug << "renumber_peer_code:   peer to naming edge component at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   type unchanged but label reversed for new crossing" << endl;
}

				odd_edge = edge_map[2*i];
				even_edge = (edge_map[code_table[OPEER][i]]-peer_first_edge-1+peer_component_edges)%peer_component_edges + peer_first_edge;

				/* the new naming edge is on the opposite strand to the old naming edge */
				if (code_table[LABEL][i] == generic_code_data::POSITIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::NEGATIVE;
				else if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::POSITIVE;
				else if (code_table[LABEL][i] == generic_code_data::FLAT)
					new_code_table[LABEL][even_edge/2] = generic_code_data::FLAT;
				else
					new_code_table[LABEL][even_edge/2] = generic_code_data::VIRTUAL;

				/* since the old naming edge is now odd, reversing the peer of the old naming 
				   strand leaves the crossing type unchanged*/
				new_code_table[TYPE][even_edge/2] = code_table[TYPE][i];
			}
			else
			{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " now odd" << endl;
	debug << "renumber_peer_code:   neither component at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   both type and label reversed for new crossing" << endl;
}

				odd_edge = edge_map[2*i];
				even_edge = edge_map[code_table[OPEER][i]];

				/* the new naming edge is on the opposite strand to the old naming edge */
				if (code_table[LABEL][i] == generic_code_data::POSITIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::NEGATIVE;
				else if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::POSITIVE;
				else if (code_table[LABEL][i] == generic_code_data::FLAT)
					new_code_table[LABEL][even_edge/2] = generic_code_data::FLAT;
				else
					new_code_table[LABEL][even_edge/2] = generic_code_data::VIRTUAL;

				/* since the old naming edge is now odd the crossing type is reversed */
				if (code_table[TYPE][i] == generic_code_data::TYPE1)
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE2;
				else
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE1;
			}			
		}
		else
		{
			if (shift[component] < 0 && shift[peer_component] < 0)
			{

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code: current code data:" << endl;
	print_code_data(code_data,debug,"renumber_peer_code: ");	
	debug << "renumber_peer_code: current code_table:" << endl;
	print(code_table,debug, 3, "renumber_peer_code: ");	
	debug << "renumber_peer_code: code_table[LABEL][" << i << "] = " << code_table[LABEL][i] << endl;
	
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " still even" << endl;
	debug << "renumber_peer_code:   both components at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   both type and label reversed for new crossing" << endl;
}

				odd_edge = (edge_map[2*i]-first_edge-1+component_edges)%component_edges + first_edge;
				even_edge = (edge_map[code_table[OPEER][i]]-peer_first_edge-1+peer_component_edges)%peer_component_edges + peer_first_edge;

				/* the new naming edge is on the opposite strand to the old naming edge */
				if (code_table[LABEL][i] == generic_code_data::POSITIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::NEGATIVE;
				else if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::POSITIVE;
				else if (code_table[LABEL][i] == generic_code_data::FLAT)
					new_code_table[LABEL][even_edge/2] = generic_code_data::FLAT;
				else
					new_code_table[LABEL][even_edge/2] = generic_code_data::VIRTUAL;

				/* since the old naming edge is still even, reversing both strands reverses the crossing type */
				if (code_table[TYPE][i] == generic_code_data::TYPE1)
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE2;
				else
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE1;
			}
			else if (shift[component] < 0)
			{

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " still even" << endl;
	debug << "renumber_peer_code:   naming edge component at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   type unchanged but label reversed for new crossing" << endl;
}

				odd_edge = (edge_map[2*i]-first_edge-1+component_edges)%component_edges + first_edge;
				even_edge = edge_map[code_table[OPEER][i]];

				/* the new naming edge is on the opposite strand to the old naming edge */
				if (code_table[LABEL][i] == generic_code_data::POSITIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::NEGATIVE;
				else if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::POSITIVE;
				else if (code_table[LABEL][i] == generic_code_data::FLAT)
					new_code_table[LABEL][even_edge/2] = generic_code_data::FLAT;
				else
					new_code_table[LABEL][even_edge/2] = generic_code_data::VIRTUAL;

				/* since the old naming edge is still even, reversing the old naming 
				   strand leaves the crossing type unchanged */
				new_code_table[TYPE][even_edge/2] = code_table[TYPE][i];
			}
			else if (shift[peer_component] < 0)
			{

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " still even" << endl;
	debug << "renumber_peer_code:   peer to naming edge component at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   type reversed but label unchanged for new crossing" << endl;
}

				even_edge = edge_map[2*i];
				odd_edge = (edge_map[code_table[OPEER][i]]-peer_first_edge-1+peer_component_edges)%peer_component_edges + peer_first_edge;

				/* the new naming edge is on the same strand as the old naming edge */
				new_code_table[LABEL][even_edge/2] = code_table[LABEL][i];

				/* since the old naming edge is still even, reversing the peer of the old naming 
				   strand reverses the crossing type */
				if (code_table[TYPE][i] == generic_code_data::TYPE1)
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE2;
				else
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE1;
			}
			else
			{

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " still even" << endl;
	debug << "renumber_peer_code:   neither component at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   both type and label unchanged for new crossing" << endl;
}

				even_edge = edge_map[2*i];
				odd_edge = edge_map[code_table[OPEER][i]];

				/* the new naming edge is on the same strand as the old naming edge */
				new_code_table[LABEL][even_edge/2] = code_table[LABEL][i];

				/* since the old naming edge is still even the crossing type is unchanged */
				new_code_table[TYPE][even_edge/2] = code_table[TYPE][i];
			}			
		}
		
		new_code_table[OPEER][even_edge/2] = odd_edge;
		new_code_table[EPEER][(odd_edge-1)/2] = even_edge;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code: renumbered as odd_edge = " << odd_edge << ", even_edge = " << even_edge << endl;
	debug << "renumber_peer_code: new type =  ";
	if (new_code_table[TYPE][even_edge/2] == generic_code_data::TYPE1)
		debug << "generic_code_data::TYPE1";
	else
		debug << "generic_code_data::TYPE2";
	
	debug << ", new label = ";
	if (new_code_table[LABEL][even_edge/2] == generic_code_data::POSITIVE)
		debug << "generic_code_data::POSITIVE";
	else if (new_code_table[LABEL][even_edge/2] == generic_code_data::NEGATIVE)
		debug << "generic_code_data::NEGATIVE";
	else if (new_code_table[LABEL][even_edge/2] == generic_code_data::VIRTUAL)
		debug << "generic_code_data::VIRTUAL";
	else
		debug << "generic_code_data::FLAT";
	debug << endl;	
}
	}
	
	/* now we can write the originating and terminating vertices and edges */
	vector<int> term_crossing(num_edges);
	vector<int> orig_crossing(num_edges);
	
	for (int i=0; i< num_crossings; i++)
	{
		term_crossing[2*i] = i;
		orig_crossing[2*i+1] = i;
		term_crossing[new_code_table[OPEER][i]] = i;
		
		/* we need to identify the edge following the naming edge's peer */
		int component = new_code_table[COMPONENT][(new_code_table[OPEER][i]-1)/2];
		int peer_successor = (new_code_table[OPEER][i]+1 - first_edge_on_component[component])%
		                     num_component_edges[component] + first_edge_on_component[component];

		orig_crossing[peer_successor] = i;
		
		new_code_table[EVEN_TERMINATING][i] = 2*i;
		new_code_table[ODD_ORIGINATING][i] = 2*i+1;
		new_code_table[ODD_TERMINATING][i] = new_code_table[OPEER][i];
		new_code_table[EVEN_ORIGINATING][i] = peer_successor;
	}
	
	code_data.code_table = new_code_table;
	code_data.term_crossing = term_crossing;
	code_data.orig_crossing = orig_crossing;

	/* re-set the head, if the code data is for a knotoid.  Note that if the head_semi_arc is even then 
	   the head is already set correctly, even if the orientation of the first componenet is reversed
   */
	if (code_data.head !=-1 && head_semi_arc %2)
    {	
		code_data.head = code_data.code_table[EPEER][(head_semi_arc-1)/2]/2;
    }

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "renumber_peer_code: renumbered code data:" << endl;
	print_code_data(code_data,debug,"renumber_peer_code: ");	
}

}

/* satellite_code_data calculates the code data of the satellite link, L, writing it to satellite_data, 
   formed by taking a number of concentric copies, S, given by strands, of the unlink in a solid torus 
   and a companion knot, K, described by knot_data.  For the purposes of calculating the satellite_data,
   these additional strands are added to the left of K, following the orientation of K determined by
   knot_data.
   
   If the writhe (sum of signs) of K is w, the solid torus is twisted -2w times, which produces a satellite 
   whose writhe is w*S.  Without the twists, the writhe of the satellite would be w*S^2, and each pair of 
   twists introduces S*(S-1) crossings of the same sign, so w*S^2-w*S(S-1) = w*S.  (As an aside, not relevent 
   to this function, note that adding a Reidemeister I move to the knot would require one more, or fewer, pair 
   of twists to the torus to maintain this write condition.)
   
   The function requires that knot_data be a knot, if it is a link it sets satellite_data to knot_data and 
   returns false, otherwise it calculates satellite data and returns true.
*/
bool satellite_code_data(generic_code_data& knot_data, generic_code_data& satellite_data, int strands)
{
	if (knot_data.num_components != 1)
	{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "satellite_code_data: knot data does not describe a knot" << endl;
		satellite_data = knot_data;
		return false;
	}
	else
	{
		int num_kcrossings = knot_data.num_crossings;		
		matrix<int>& kcode_table = knot_data.code_table;		

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "satellite_code_data: presented with knot data" << endl;
	print_code_data(knot_data, debug, "satellite_code_data:   ");
	debug << "satellite_code_data: number of strands required = " << strands << endl;
}

		int writhe = 0;
		
		for (int i=0; i< num_kcrossings; i++)
		{
			if (kcode_table[LABEL][i] == generic_code_data::POSITIVE || kcode_table[LABEL][i] == generic_code_data::NEGATIVE)
				writhe += kcode_table[LABEL][i];
		}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "satellite_code_data: writhe = " << writhe << endl;
		
		/* Evaluate the number of crossings in the satellite link, L: each crossing in the knot contributes S^2 crossings
		   to L and in the writhe twists there are S(S-1)/2 crossings in every twist, and there are two twists per writhe.
		*/
		int num_scrossings = num_kcrossings*strands*strands + abs(writhe)*strands*(strands-1);

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "satellite_code_data: num_scrossings = " << num_scrossings << endl;
		
		satellite_data = generic_code_data(num_scrossings,strands);
		matrix<int>& scode_table = satellite_data.code_table;		
		satellite_data.type = generic_code_data::peer_code;
		
		/* Where the strand of K crosses one other semi-arc, in the satellite, L, it crosses S semi-arcs, thus where there
		   were 2*num_kcrossings edges in the (one) component of K, each component of L (without the writhe twists) has 
		   2*num_kcrossings*S crossings.  Considering the writhe twists, each component encounters (S-1) crossings in each 
		   twist and there are two twists per writhe.
		*/
		for (int i=0; i< strands; i++)
			satellite_data.first_edge_on_component[i] = i*(2*abs(writhe)*(strands-1)+2*num_kcrossings*strands);

		for (int i=0; i< strands; i++)
			satellite_data.num_component_edges[i] = satellite_data.first_edge_on_component[1];

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "satellite_code_data: first edge on components of satellite: ";
	for (int i=0; i< strands; i++)
		debug << satellite_data.first_edge_on_component[i] << " ";
	debug << endl;
}

        /*
                                                             3
                                                  s'[s]  \  /
                                                       \  x    2
                                    (2j-1)*S [2i*S]  \  x  \  /                 \e
               2j-1 [2i]\   /                      \  x  \  x    1               x
                         \ /                        x  \  x  \  /              o/ \o
                          x                        / \  x  \  x   0            3   x    
                         / \                          x  \  x  \ /               e/ \e 
               2i [2j-1]/   \                        / \  x  \  x                2   x
                                                s[s']   x  \  x  \                 o/ \o
                                                       / \  x  \  3                1   x
                                                          x  \  2                    e/ \e
                                                         / \  1                      0
                                          2i*S [(2j-1)*S]   0



        If we consider the S^2 crossings in L corresponding to either a type I or type II crossing in K with terminating 
        edges 2i and 2j-1, we see that in strand s=0,...S-1, there are S terminating edges, "corresponding" to the edge 2i, 
        given by 
        
                      f_s + s%2 + 2i*S + k, for k=0,...,S-1,  where f_s is the first edge on strand s
                      
        We need to include s%2 because, as can be seen from the right hand diagram above, whilst all strands contain 
        2*num_kcrossings*S + 2*|w|*(S-1) edges, an even number, odd numbered strands have to be numbered with a shift of one
        (which we choose to be backwards) in order to ensure we always have a odd and even terminating edge at each crossing.
        
        In the case of a type I crossing, these edges terminate at crossings on strand s' = 0,...S-1 respectively.  In the 
        case of a type II crossing they terminate on strands s' = S-1,...,0 respectively.
        
        Each of the stands s' contains S edges, "corresponding" to edge 2j-1, given by
        
                      f_s' + s'%2 + (2j-1)*S + k, for k=0,...,S-1
        
        which themselves meet strands s = S-1,...,0 respectively in a type I crossing and strands s=0,...,S-1 in a type II
        crossing.
        
        For k=0,...,S-1, the enumeration S-1,...0 is just (S-1-k)

        */

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "satellite_code_data: evaluate satellite non-writhe crossings" << endl;

        for (int i=0; i< num_kcrossings; i++)
        {
	
			int odd_peer = kcode_table[OPEER][i];
			bool type_1_crossing = (kcode_table[TYPE][i] == generic_code_data::TYPE1? true: false);
			
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "satellite_code_data:   crossing " << i;
	debug << ", odd peer = " << odd_peer << ", type " << (type_1_crossing? "1":"2") << endl;
}			
			for (int s=0; s < strands; s++)
			{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "satellite_code_data:     strand " << s<< endl;
				
				for (int k=0; k < strands; k++)
				{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "satellite_code_data:       k= " << k << ": ";
					int s_edge = satellite_data.first_edge_on_component[s] + s%2 + 2*i*strands + k;
					int s_prime = (type_1_crossing? k: strands-1-k);
					int s_prime_edge = satellite_data.first_edge_on_component[s_prime] + s_prime%2 + odd_peer*strands + (type_1_crossing? strands-1-s: s);

					if ( s_prime_edge == satellite_data.first_edge_on_component[s_prime] + satellite_data.num_component_edges[s_prime])
							s_prime_edge = satellite_data.first_edge_on_component[s_prime];					

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "s_edge = " << s_edge << ", s' = " << s_prime << ", s_prime_edge = " << s_prime_edge << endl;


					int required_type;
					if ((type_1_crossing && s_prime_edge%2) || (!type_1_crossing && s_edge%2))
						required_type = generic_code_data::TYPE1;
					else
						required_type = generic_code_data::TYPE2;
						
					int required_crossing;
					if ((kcode_table[LABEL][i] == generic_code_data::POSITIVE && kcode_table[TYPE][i] == generic_code_data::TYPE2) ||
					     (kcode_table[LABEL][i] == generic_code_data::NEGATIVE && kcode_table[TYPE][i] == generic_code_data::TYPE1))
					    required_crossing = generic_code_data::POSITIVE;
					else if ((kcode_table[LABEL][i] == generic_code_data::POSITIVE && kcode_table[TYPE][i] == generic_code_data::TYPE1) ||
					     (kcode_table[LABEL][i] == generic_code_data::NEGATIVE && kcode_table[TYPE][i] == generic_code_data::TYPE2))
					    required_crossing = generic_code_data::NEGATIVE;
					else
						required_crossing = kcode_table[LABEL][i];

					assign_satellite_code(satellite_data, s_edge, s_prime_edge, s, s_prime, required_type, required_crossing);									
				}
			}
		}
		
		/* consider the writhe twists, set the base edge to the first edge on each strand involved in the
		   writhe twists
		*/
		vector<int> base_edge(strands);
		for (int i=0; i< strands; i++)
			base_edge[i] = satellite_data.first_edge_on_component[i]+i%2+2*num_kcrossings*strands;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "satellite_code_data: writhe twists base_edge: ";
	for (int i=0; i< strands; i++)
		debug << base_edge[i] << " ";
	debug << endl;
}

		
		for (int w=0; w< abs(writhe); w++)
		{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "satellite_code_data: writhe twist " << w+1 << endl;
	
			for (int s=0; s< strands-1; s++)
			{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "satellite_code_data:   strand " << s << endl;
				for (int k=0; k<=strands-2-s;k++)
				{
					int s_prime=strands-1-k;
					int s_edge = base_edge[s]+k;
					int s_prime_edge = base_edge[s_prime]+strands-2-s;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "satellite_code_data:   offset " << k << ", s_prime = " << s_prime << endl;
	debug << "satellite_code_data:   s_edge = " << s_edge << ", s_prime_edge = " << s_prime_edge << endl;
}
					assign_satellite_code(satellite_data, s_edge, s_prime_edge, s, s_prime, 
										  (s_prime_edge%2? generic_code_data::TYPE1 : generic_code_data::TYPE2),
										  (writhe < 0? generic_code_data::POSITIVE: generic_code_data::NEGATIVE));


					/* now the second meeting */
					s_edge=base_edge[s]+strands-1+s_prime-1;
						
					if ( s_edge == satellite_data.first_edge_on_component[s] + satellite_data.num_component_edges[s])
							s_edge = satellite_data.first_edge_on_component[s];					
					
					s_prime_edge = base_edge[s_prime]+strands-1+s;

					if ( s_prime_edge == satellite_data.first_edge_on_component[s_prime] + satellite_data.num_component_edges[s_prime])
							s_prime_edge = satellite_data.first_edge_on_component[s_prime];					

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "satellite_code_data:   second meeting s_edge = " << s_edge << ", s_prime_edge = " << s_prime_edge << endl;

					assign_satellite_code(satellite_data, s_edge, s_prime_edge, s, s_prime, 
										  (s_prime_edge%2? generic_code_data::TYPE2 : generic_code_data::TYPE1),
										  (writhe < 0? generic_code_data::POSITIVE: generic_code_data::NEGATIVE));
				}				
			}
			
			/* adjust base_edge for the next write twists */
			for (int i=0; i< strands; i++)
				base_edge[i] += 2*(strands-1);
			
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "satellite_code_data: updated writhe twists base_edge: ";
	for (int i=0; i< strands; i++)
		debug << base_edge[i] << " ";
	debug << endl;
}
		}

		/* complete the generic code data of the satellite by evaluating the originating and terminating crossings
		   for each edge
		*/
		for (int i =0; i < satellite_data.num_crossings; i++)
		{
			satellite_data.term_crossing[scode_table[EVEN_TERMINATING][i]] = i;
			satellite_data.term_crossing[scode_table[ODD_TERMINATING][i]] = i;
			satellite_data.orig_crossing[scode_table[EVEN_ORIGINATING][i]] = i;
			satellite_data.orig_crossing[scode_table[ODD_ORIGINATING][i]] = i;
		}

//		print_code_data(satellite_data, cout, "    ");
		
//		cout << "\nsatellite_code_data: number of strands = " << strands << " writhe = " << writhe << endl;
		return true;
	}
}

void assign_satellite_code(generic_code_data& satellite_data, int s_edge, int s_prime_edge, int s_component, int s_prime_component, int type, int crossing)
{
	matrix<int>& scode_table = satellite_data.code_table;
	
	int s_even;
	int s_odd;
	
	if (s_edge %2)
	{
		s_odd = s_edge;
		s_even = s_prime_edge;
	}
	else
	{
		s_odd = s_prime_edge;
		s_even = s_edge;
	}
					
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "assign_satellite_code_data:       s_even = " << s_even << ", s_odd = " << s_odd << endl;
	
	int s_crossing = s_even/2;

	scode_table[OPEER][s_crossing] = s_odd;
	scode_table[EPEER][(s_odd-1)/2] = s_even;

	if (type == generic_code_data::TYPE1)
		scode_table[TYPE][s_crossing] = generic_code_data::TYPE1;
	else
		scode_table[TYPE][s_crossing] = generic_code_data::TYPE2;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "assign_satellite_code_data:       crossing is " << (scode_table[TYPE][s_crossing] == generic_code_data::TYPE1? "Type I" : "Type II") << endl;
					
	/* crossing indicates whether we have a positive, negative crossing in the standard sense, so we adjust
	   for the peer code based on type.
	*/
	if (crossing == generic_code_data::POSITIVE)
	{
		if (scode_table[TYPE][s_crossing] == generic_code_data::TYPE2)
			scode_table[LABEL][s_crossing] = generic_code_data::POSITIVE;
		else
			scode_table[LABEL][s_crossing] = generic_code_data::NEGATIVE;
	}
	else if (crossing == generic_code_data::NEGATIVE)
	{
		if (scode_table[TYPE][s_crossing] == generic_code_data::TYPE1)
			scode_table[LABEL][s_crossing] = generic_code_data::POSITIVE;
		else
			scode_table[LABEL][s_crossing] = generic_code_data::NEGATIVE;
	}
	else
	{
		scode_table[LABEL][s_crossing] = crossing; // virtual, or flat
	}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "assign_satellite_code_data:       label is ";
	if (scode_table[LABEL][s_crossing] == generic_code_data::POSITIVE)
		debug << "positive" << endl;
	else if (scode_table[LABEL][s_crossing] == generic_code_data::NEGATIVE)
		debug << "negative" << endl;
	else if (scode_table[LABEL][s_crossing] == generic_code_data::VIRTUAL)
		debug << "virtual" << endl;
	else if (scode_table[LABEL][s_crossing] == generic_code_data::FLAT)
		debug << "flat" << endl;
	else
	    debug << scode_table[LABEL][s_crossing] << endl;
}
					
	scode_table[EVEN_TERMINATING][s_crossing] = s_even;
	scode_table[ODD_TERMINATING][s_crossing] = s_odd;				


	if (s_odd == s_edge)
	{
		if ( s_odd == satellite_data.first_edge_on_component[s_component] + satellite_data.num_component_edges[s_component] - 1)
			scode_table[EVEN_ORIGINATING][s_crossing] = satellite_data.first_edge_on_component[s_component];
		else
			scode_table[EVEN_ORIGINATING][s_crossing] = s_odd+1;

		scode_table[COMPONENT][s_crossing] = s_prime_component;
	}
	else
	{
		if ( s_odd == satellite_data.first_edge_on_component[s_prime_component] + satellite_data.num_component_edges[s_prime_component] - 1)
			scode_table[EVEN_ORIGINATING][s_crossing] = satellite_data.first_edge_on_component[s_prime_component];
		else
			scode_table[EVEN_ORIGINATING][s_crossing] = s_odd+1;					

		scode_table[COMPONENT][s_crossing] = s_component;
	}
	
	scode_table[ODD_ORIGINATING][s_crossing] = s_even+1;
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "assign_satellite_code_data:       even terminating edge is " << scode_table[EVEN_TERMINATING][s_crossing] << endl;
	debug << "assign_satellite_code_data:       odd terminating edge is " << scode_table[ODD_TERMINATING][s_crossing] << endl;
	debug << "assign_satellite_code_data:       even originating edge is " << scode_table[EVEN_ORIGINATING][s_crossing] << endl;
	debug << "assign_satellite_code_data:       odd originating edge is " << scode_table[ODD_ORIGINATING][s_crossing] << endl;
	debug << "assign_satellite_code_data:       component " << s_component << endl;
}

}

/* remove_peer_code_component identifies those edges belonging to the specified component 
   together with their peer edges and then calls remove_edge_flags_from_peer_code to remove them.  
   
   the function returns the number of components removed, which may be more than 1 if the component
   to be removed intesects other components that only meet the diagram in the component we are removing.
   An example of this is removing component zero of [-15 -7 9, 11 -13, -1^ 3, -5]/ * - - + - + + +
*/
int remove_peer_code_component(generic_code_data& code_data, int component, vector<int>& component_flags)
{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "remove_peer_code_component: removing component " << component << endl;

	int num_components_removed = 0;  //the return variable
	
	if (code_data.num_components > 1)
	{
		vector<int> edge_flag(2*code_data.num_crossings);
		
		for (int i= code_data.first_edge_on_component[component]; i < code_data.first_edge_on_component[component]+code_data.num_component_edges[component]; i++ )
		{
			int crossing = code_data.term_crossing[i];
			edge_flag[code_data.code_table[ODD_TERMINATING][crossing]] = 1;
			edge_flag[code_data.code_table[EVEN_TERMINATING][crossing]] = 1;
		}
		
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "remove_peer_code_component: edge_flag = ";
	for (int i=0; i < 2*code_data.num_crossings; i++)
		debug << edge_flag[i] << ' ';
	debug << endl;
}
		num_components_removed = remove_edge_flags_from_peer_code(code_data, edge_flag,component_flags);
	}
	else  // set the number of crossings to be zero to indicate the last component has been removed.
	{
		code_data.num_crossings = 0;
		num_components_removed = 1;
	}
	
	return num_components_removed;
}

/* gauss_parity returns a vector of length the number of immersion crossings in code data, whose entries
   corresponding to non-virtual and non-shortcut crossings are set to gauss_orientation_data::parity::EVEN
   or gauss_orientation_data::parity::ODD dependent on the number of intervening terms in the Gauss data.
*/
vector<int> gauss_parity(generic_code_data& code_data)
{
	int num_crossings = code_data.num_crossings;
	gauss_orientation_data gauss_data(code_data);
	
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_parity: presented with generic code ";
	write_code_data(debug, code_data);
	debug << endl;
	debug << "gauss_parity: gauss data ";
//	write_gauss_data(gauss_data, debug);
//	debug << endl;
	print_gauss_data(gauss_data, debug,"gauss_parity: ");
}
	int num_terms = gauss_data.num_terms;
	int num_gauss_crossings = num_terms/2;
	matrix<int>& orientation_matrix = gauss_data.orientation_matrix;
	
	vector<int> parity(num_crossings); // initializes to zero, i.e gauss_orientation_data::parity::NONE
	
	for (int i=0; i< num_gauss_crossings; i++)
	{
		int immersion_crossing = gauss_data.immersion_crossing[i];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_parity: gauss_crossing = " << i+1 << " immersion_crossing = " << immersion_crossing << endl;
	
		int start=0;
		for (int j=0; j< num_terms; j++)
		{
			if (orientation_matrix[1][j] == i+1)
			{
				start=j;
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_parity:   first ocurrence at index " << start << endl;
				break;
			}
		}
		
		for (int k=1; k< num_terms; k++)
		{
			if (orientation_matrix[1][start+k] == i+1)
			{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_parity:   second ocurrence at index " << start+k << endl;
				if (k%2)
					parity[immersion_crossing] = gauss_orientation_data::parity::EVEN;
				else
					parity[immersion_crossing] = gauss_orientation_data::parity::ODD;
				
				break;
			}
		}
	}

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_parity: crossing number ";
	for (int i=0; i< num_crossings; i++)
		debug << i << " ";
	debug << endl;
	debug << "gauss_parity: parity          ";
	for (int i=0; i< num_crossings; i++)
	{
		if (parity[i] == gauss_orientation_data::parity::ODD)		
			debug << "O ";
		else if (parity[i] == gauss_orientation_data::parity::EVEN)		
			debug << "E ";
		else
			debug << "N ";
	}
	debug << endl;
}
	
	return parity;
}

/* remove_virtual_components removes any component in the peer code that contains only virtual crossings 
   the function returns the number of crossings it has removed.
*/
int remove_virtual_components(generic_code_data& code_data, vector<int>& component_flags)
{
	int num_virtual_components = 0;
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "remove_virtual_components: processing peer code ";
	write_peer_code(debug,code_data);
	debug << endl;
}

	for (int i=0; i< code_data.num_components; i++)
	{
		bool virtual_component = true;
		for (int j=code_data.first_edge_on_component[i]; j< code_data.first_edge_on_component[i]+code_data.num_component_edges[i]; j++)
		{
			if (code_data.code_table[LABEL][code_data.term_crossing[j]] != generic_code_data::VIRTUAL)
			{
				virtual_component = false;
				break;
			}
		}
				
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "remove_virtual_components: component " << i << (virtual_component? " is " : " is not ") << "a virtual component" << endl;
	
		if (virtual_component && code_data.num_components > 1)
		{
			num_virtual_components += remove_peer_code_component(code_data, i,component_flags);			
		}
		else if (virtual_component)
		{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "remove_virtual_components: no other components, setting num_crossings = 0" << endl;
	
			num_virtual_components = 1;
			code_data.num_crossings = 0;
			for (unsigned int i=0; i< component_flags.size(); i++)
				component_flags[i] = 0;			
		}
	}
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "remove_virtual_components: removed " << num_virtual_components << " components" << endl;
	
	return num_virtual_components;
}

/* un-oriented minimal representation of peer codes only */
string minimal_peer_code_representation (generic_code_data peer_code_data)
{
	if (peer_code_data.type != generic_code_data::code_type::peer_code	)
		return "Error! minimal_peer_code_representation supports only peer codes";

	int num_crossings = peer_code_data.num_crossings;
	int num_components = peer_code_data.num_components;
	vector<int>& num_component_edges = peer_code_data.num_component_edges;
	
	generic_code_data minimum_peer_code_data(peer_code_data);


	/* If the code represents a link and any component is shifted by an odd number of edges, 
	   we must have all intersecting components shifted by an odd number of edges, otherwise 
	   we will violate the requirement to have an odd and even edge terminating at every crossing.
	   Since the peer code has to be connected to be realizable, this means that every component 
	   has to be shifted by an odd number of edges.  Thus we have to shift every component by an odd
	   number of edges or every component by an even number of edges.  
	  
	   If a component's orientation is reversed we also violate the requirement to have an odd and 
	   even edge terminating at every crossing.  Therefore an orientation reversal must be accompanied by
	   an odd shift of the component being reversed or of those other components it meets at a crossing.
		
	   We create a vector shift whose ith element lies in the range -num_component_edges[i] to
	   num_component_edges[i] -1.  A negative value for shift[i] indicates an orientation reversal after
	   a shift of abs(shift[i]) edges.  We include -num_component_edges[i] in the range because it 
	   corresponds to orientation reversal of a component without any shift (and we don't have -0 at 
	   our disposal of course). By considering orientation reversal an "odd" operation (in the 
	   light of the above coments), a particular shift vector is valid if every entry is odd or every entry
	   is even.  Thus -2 is odd and -1 is even.  Note no attempt is made to avoid checking the zero vector 
	   in the following code.
	*/
	vector<int> shift(num_components); 
	for (int j=0; j< num_components; j++)
		shift[j] = -1 * num_component_edges[j];
	
	bool next_shift_vector_found;
	do
	{
					   
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
   	debug << "minimal_peer_code_representation:     shift vector: ";
   	for (int j=0; j< num_components; j++)
		debug << shift[j] << ' ';
//	debug << endl;
}

		/* check the shift vector is valid */
		bool valid_shift_vector = true;
		int root_parity = (shift[0]%2? -1 : 1);
		if (shift[0] < 0)
			root_parity *= -1;
		
		for (int j=1; j< num_components; j++)
		{
			int parity = (shift[j]%2? -1 : 1);
			if (shift[j] < 0)
				parity *= -1;
			
			if (parity != root_parity)
			{
				valid_shift_vector = false;
				break;
			}
		}

		if (valid_shift_vector)
		{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
  	debug << " is valid, testing resultant renumbering" << endl;

			generic_code_data renumbered_code_data(peer_code_data);						
			renumber_peer_code(renumbered_code_data, shift);
			
			matrix<int>& mcode_table = minimum_peer_code_data.code_table;
			matrix<int>& rcode_table = renumbered_code_data.code_table;
							
							
			bool new_minimal_rep_found = false;
			/* check the lexicographic ordering of the peer code, ignoring components, that is
			   compare the product of the TYPE and OPEER
			*/
			for (int i=0; i< num_crossings; i++)
			{
				if (rcode_table[TYPE][i]*rcode_table[OPEER][i] < mcode_table[TYPE][i]*mcode_table[OPEER][i])
				{
					new_minimal_rep_found = true;
					break;
				}
				else if (rcode_table[TYPE][i]*rcode_table[OPEER][i] > mcode_table[TYPE][i]*mcode_table[OPEER][i])
				{
					break;
				}
			}
			
			if (new_minimal_rep_found)
			{
				minimum_peer_code_data = renumbered_code_data;
			}
		}
   	
		/* is there another valid shift vector? */
		next_shift_vector_found = false;
		for (int j=num_components-1; j >= 0; j--)
		{
			if (shift[j] < num_component_edges[j]-1)
			{
				shift[j] ++;
				next_shift_vector_found = true;
				break;
			}
			else
			{
				shift[j] = -1 * num_component_edges[j];
			}
		}
	} while(next_shift_vector_found);				

	ostringstream oss;
	write_peer_code(oss,minimum_peer_code_data);
	return oss.str();
}

/* partition_peer_code partitions code_data into the first connected component and the subsequent component(s).
   It returns generic_code_data describing the subsequent components(s) and adjusts code_data to reflect only the
   first connected component.
      
   If the code_data character is not CLOSED, code_data.head may not be -1, in which case the segment component must 
   be part of the first connected component, since the leg semi_arc is numbered zero.  That means the head and 
   immersion character is set correctly in code_data and that the default initialization of the head and immersion 
   character in subsequent_code_data is also correct.
*/					
generic_code_data partition_peer_code(generic_code_data& code_data, vector<int>& component_flags)
{
if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "partition_peer_code: presented with code data ";
	write_code_data(debug,code_data);	
	debug << endl;
	debug << "partition_peer_code: component_flags: ";
	for (unsigned int j=0; j< component_flags.size(); j++)
		debug << component_flags[j] << ' ';
	debug << endl;
}
	generic_code_data first_code_data;
	generic_code_data subsequent_code_data;
	subsequent_code_data.type = generic_code_data::code_type::peer_code;

	bool pure_knotoid_code_data = false;
	int head_semi_arc = -1;
	if (code_data.immersion == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && code_data.shortcut_crossing.size())
	{
		pure_knotoid_code_data = true;
		
		if (code_data.code_table[LABEL][code_data.head] == generic_code_data::POSITIVE)
			head_semi_arc = code_data.code_table[OPEER][code_data.head];
		else if (code_data.code_table[LABEL][code_data.head] == generic_code_data::NEGATIVE)
			head_semi_arc = 2*code_data.head;
		else
		{
			cout << "\nError! Function partition_peer_code presented with a knotoid code whose first shortcut crossing is virtual." << endl;
			exit(0);
		}
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "partition_peer_code: head_semi_arc = " << head_semi_arc << endl;			
	}

	bool track_zig_zag_counts;
	if (code_data.zig_zag_count.numrows() !=0)
		track_zig_zag_counts = true;
	else
		track_zig_zag_counts = false;

	list<int> span_list = vertex_span (code_data, 0);
	
	int num_f_crossings = span_list.size();
	
if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "partition_peer_code: vertex span reaches " << num_f_crossings << " crossings" << endl;
	debug << "partition_peer_code: span_list: ";
	list<int>::iterator lptr=span_list.begin();
	while (lptr != span_list.end())
	{
		debug << *lptr << ' ';
		lptr++;
	}
	debug << endl;
	
}
	if (num_f_crossings != code_data.num_crossings)
	{
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "partition_peer_code: code_data not connected, first connected component contains " << num_f_crossings << " crossings" << endl;

		/* code_data is not connected and num_f_crossings identifies the number of crossings belonging to the first connected component, 
		   which must comprise all of the crossings on the unicursal components making up the first connected component.
		*/
		vector<int> f_components(code_data.num_components);
		list<int>::iterator lptr=span_list.begin();
		while (lptr != span_list.end())
		{
			f_components[code_data.code_table[COMPONENT][*lptr]] = 1;
			lptr++;
		}
		
if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "partition_peer_code: f_components: " ;
	for (int i=0; i<code_data.num_components; i++)
		debug << f_components[i] << ' ';
	debug << endl;
}
		/* update component_flags to reflect the components included in the first part */
		for (int i=0; i< code_data.num_components; i++)
		{
			if (f_components[i] == 0)
				clear_component_flag(i,component_flags);
		}
		
		/* evaluate the peer code for each partition by creating a generic code object with enough information to write the 
		   partition's peer code to an ostringstream and then read it back again.
		   
		   We need the head, num_crossings and the OPEER, TYPE, LABEL and COMPONENT rows of the code_table in the generic_code_data 
		   structure to do this.		
		*/				
		int num_s_crossings = code_data.num_crossings - num_f_crossings;
		
		for (int i=0; i< 2; i++)
		{
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "partition_peer_code: determine the peer code of partition " << i << endl;
	
			int num_partition_crossings = (i==0? num_f_crossings: num_s_crossings);
			matrix<int> partition_code_table(CODE_TABLE_SIZE,num_partition_crossings,-1);
			matrix<int> partition_zig_zag_count(2,num_partition_crossings);
			
			/* write new edge labels into the even (row 0)  and odd (row 1) rows of new_partition_labels, 
			   in the location corresponding to the old edge labels in the EVEN_TERMINATING and ODD_TERMINATING
			   rows of code_data.code_table.  In row 2 of new_partition_labels we record the component number
			   of the unicursal component in the partition.
			*/
			matrix<int> new_partition_labels(3,code_data.num_crossings,-1);
			int edge = 0;
			int component = -1;
			int new_head = -1;
			
			for (int j=0; j< code_data.num_components; j++)
			{
				/* f_components is a vector of flags, with f_components[j] == 1 iff the jth component is part of the first partition,
				   therefore for the first partition we want flags with value 1 and for the subsequent partition, flags with value 0
				*/
				if (f_components[j] == i)
				{
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "partition_peer_code:   unicursal component " << j << " is not part of partition " << i << endl;
					
					continue;					
				}
				
				component++;
				
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "partition_peer_code:   unicursal component " << j << " is part of partition " << i << endl;
					
				int num_component_edges = code_data.num_component_edges[j];
				int first_edge = code_data.first_edge_on_component[j];
				
				for (int k=0; k < num_component_edges; k++)
				{
					int old_edge = first_edge + k;
					int row;
					int column;
										
					if (old_edge % 2 == 0)
					{
						row = 0; // even
						column = old_edge/2;
						
						new_partition_labels[2][column] = component;
					}
					else
					{
						row = 1; // odd
						column = code_data.code_table[EPEER][(old_edge-1)/2]/2;
					}
					
if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "partition_peer_code:     old_edge " << old_edge << " is (row) " << (row == 0? "EVEN_TERMINATING": "ODD_TERMINATING") 
	      <<  " edge at old crossing (column) " << column << endl;
	debug << "partition_peer_code:     new edge = " << edge;
	if (old_edge % 2 == 0)
		debug << ", component = " << component << endl;
	else
		debug << endl;
}					
					new_partition_labels[row][column] = edge;
					edge++;					
				}
						
				/* write the new edge labels into the correct location in the OPEER row of a partition code table, 
				   and set the TYPE LABEL and COMPONENT rows 
				*/
				for (int j=0; j < code_data.num_crossings; j++)
				{
					if (new_partition_labels[0][j] != -1)
					{
						int crossing = new_partition_labels[0][j]/2;
						partition_code_table[OPEER][crossing] = new_partition_labels[1][j];
						partition_code_table[COMPONENT][crossing] = new_partition_labels[2][j];
						partition_code_table[TYPE][crossing] = code_data.code_table[TYPE][j];
						partition_code_table[LABEL][crossing] = code_data.code_table[LABEL][j];
						
						if (track_zig_zag_counts)
						{
							partition_zig_zag_count[0][crossing] = code_data.zig_zag_count[0][j];
							partition_zig_zag_count[1][crossing] = code_data.zig_zag_count[1][j];
						}
					}
				}
				
			}
			
if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "partition_peer_code: new_partition_labels even terminating: ";
	for (int k=0; k< code_data.num_crossings; k++)
		debug << setw(3) << new_partition_labels[0][k];
	debug << "\npartition_peer_code: new_partition_labels odd terminating:  ";
	for (int k=0; k< code_data.num_crossings; k++)
		debug << setw(3) << new_partition_labels[1][k];
	debug << endl;
	debug << "partition_peer_code: partition_code_table" << endl;
	debug << "partition_peer_code:  odd peer: ";
	for (int j=0; j< num_partition_crossings; j++)
		debug << partition_code_table[OPEER][j] << ' ';
	debug << endl;
	debug << "partition_peer_code:  component: ";
	for (int j=0; j< num_partition_crossings; j++)
		debug << partition_code_table[COMPONENT][j] << ' ';
	debug << endl;
	debug << "partition_peer_code:  type: ";
	for (int j=0; j< num_partition_crossings; j++)
		debug << partition_code_table[TYPE][j] << ' ';
	debug << endl;
	debug << "partition_peer_code:  label: ";
	for (int j=0; j< num_partition_crossings; j++)
		debug << partition_code_table[LABEL][j] << ' ';
	debug << endl;
	
	if (track_zig_zag_counts)
	{
		debug << "partition_peer_code: partition_zig_zag_count: " << endl;
		print(partition_zig_zag_count,debug, 3,"partition_peer_code: ");	
	}
}
			generic_code_data partition_code_data;
			partition_code_data.type = generic_code_data::code_type::peer_code;
			partition_code_data.num_crossings = num_partition_crossings;
			partition_code_data.code_table = partition_code_table;
						
			ostringstream oss;
			write_peer_code(oss, partition_code_data);
			
			if (i==0)
			{
				read_peer_code(first_code_data, oss.str());

				first_code_data.immersion = code_data.immersion;
				first_code_data.head_zig_zag_count = code_data.head_zig_zag_count;				
				
				if (pure_knotoid_code_data)
				{
					first_code_data.head = new_head;
					
					if (head_semi_arc % 2)
						first_code_data.head = first_code_data.code_table[EPEER][(head_semi_arc-1)/2]/2;
					else
						first_code_data.head = head_semi_arc/2;
						
					valid_knotoid_input(first_code_data);  // sets the shortcut crossings
				}
					
			}
			else
			{
				read_peer_code(subsequent_code_data, oss.str());			
			}	

			if (track_zig_zag_counts)
			{
				/* check that this partition has a non-zero zig-zag_count */
				bool zig_zags_present = false;
				for (int j=0; j<num_partition_crossings; j++)
				{
					if (partition_zig_zag_count[0][j] != 0 || partition_zig_zag_count[1][j] !=0)
					{
						zig_zags_present = true;
						break;
					}
				}
				
				if (zig_zags_present || (i==0 && first_code_data.head_zig_zag_count != 0))
				{
					if (i==0)
						first_code_data.zig_zag_count = partition_zig_zag_count;
					else
						subsequent_code_data.zig_zag_count = partition_zig_zag_count;
				}					
			}
			
if (braid_control::DEBUG >= braid_control::BASIC)
{
	if (i==0)
	{
		debug << "partition_peer_code:   first partition peer code ";
		write_peer_code(debug,first_code_data);
		debug << endl;
		debug << "partition_peer_code:   first partition code_data" << endl;
		print_code_data(first_code_data, debug, "partition_peer_code:   ");
		debug << "partition_peer_code:   component_flags updated to: ";
		for (unsigned int j=0; j< component_flags.size(); j++)
			debug << component_flags[j] << ' ';
		debug << endl;
	}
	else
	{
		debug << "partition_peer_code:   subsequent partition peer code ";
		write_peer_code(debug,subsequent_code_data);
		debug << endl;
		debug << "partition_peer_code:   subsequent partition code_data" << endl;
		print_code_data(subsequent_code_data, debug, "partition_peer_code:   ");
	}
}							
		}
		
		code_data = first_code_data;														
	}
	return subsequent_code_data;
}


/* cycle_gauss_code_labels shifts the start of the numbering of the component containing edge back by one, 
   so that the number of each edge in the component is incremented (around the component)   
*/
void cycle_gauss_code_labels(matrix<int>& crossing_peers, vector<int>& num_component_edges, vector<int>& first_edge_on_component, int num_crossings,int num_components, int edge)
{

if (braid_control::DEBUG >= braid_control::EXHAUSTIVE)
{
	debug << "cycle_gauss_code_labels: initial crossing_peers:" << endl;
	print(crossing_peers,debug,4,"cycle_gauss_code_labels: ");	
	debug << "cycle_gauss_code_labels: cycling component containing edge " << edge << endl;
}

	/* identify the component containing edge */
	int component = 0;
	for (int i=1; i< num_components; i++)
	{
		if (first_edge_on_component[i] > edge)
			break;
		else
			component++;
	}

	int first = first_edge_on_component[component];
	int last = first_edge_on_component[component] + num_component_edges[component]-1;

if (braid_control::DEBUG >= braid_control::EXHAUSTIVE)
	debug << "cycle_gauss_code_labels: edge lies on component " << component << ", first edge on component = " << first << ", last edge = " << last << endl;

	matrix<int> new_crossing_peers = crossing_peers;
	
	/* cycle the labels in new_crossing_peers belonging to the identified component */
	for (int i=first; i <= last; i++)
	{
		/* find edge i in the code table */
		for (int j=0; j< num_crossings; j++)
		{
			if (crossing_peers[j][0] == i)
			{
				if (i == last)
					new_crossing_peers[j][0] = first;
				else
					new_crossing_peers[j][0]++;
									
				break;
			}
			else if (crossing_peers[j][1] == i)
			{
				if (i == last)
					new_crossing_peers[j][1] = first;
				else
					new_crossing_peers[j][1]++;
									
				break;
			}
		}
	}

	crossing_peers = new_crossing_peers;
	
if (braid_control::DEBUG >= braid_control::EXHAUSTIVE)
{
	debug << "cycle_gauss_code_labels: resultant crossing_peers:" << endl;
	print(crossing_peers,debug,4,"cycle_gauss_code_labels: ");	
}
	
}



/* gauss_to_peer_code converts the gauss code of a classical link diagram to a peer code.  It is up
   to the calling environment to ensure that the gauss_code_data supplied to this function is indeed that of 
   a clasical link rather than a virtual link, since Gauss codes contain no indication of the presence of
   virtual crossings.
   
   The first parameter is passed by value so that we can modify the code_table of the Gauss code if we need
   to shift edges.
   
   
The layout of a virtual diagram is constructed in layers top-down on the page.  Starting from the first crossing in the Gauss code, the crossing
is set out with the incoming edges on the left and the four edges are extended in downwards facing arcs ("combed"), as shown in the left diagram below
  
            -----------     -----------
           |           \   /           |                                                          ---------------------------------
           |            \ /            |                  +--------------------------+           |  +--------------------------+   |
           |             x             |                  |         diagram          |           |  |                          |   |
           |            / \            |                  +--------------------------+           |  +--------------------------+   |
           |           /   \           |                     |  |   . . .      |  |              |     |  |   . . .      |  |      |
           |        ---     ---        |                                                                                     ------
           |       |           |       |
           

The four edges are called the "fringe" of the diagram.  The algorithm proceeds by adding a crossing that connects to one or more edges of the fringe, introducing virtual crossings
(described below) as necessary, and those edges that do not connect to the frings are combed downwards, without introducing any further crossings to create a new fringe.  Thus at
each intermediate stage of the algorithm, we have a diagram with a fringe below it, as shown in the middle diagram above.  The number of edges in the fringe depends upon the number of 
matching edges in the last fringe and the last crossing that was added to the diagram.  At the end of the algorithm, the fringe will contain four edges, to which the last crossing is 
added.
           
In order to attach the next crossing to a fringe, we require that the edges in the fringe are adjacent to each other and appear in the same order as they do around the crossing. 
Moving from left to right along the fringe the edges must appear in the corresponding order when moving clockwise around the crossing.
If there
is just one matching edge, this is trivially satisfied, if there are two or more matching edges that are adjacent to each other but in the wrong order we can add virtual crossings to
permute the order of matching crossings in a new fringe.

Reidemeister I loops
====================

If the gauss code presented to the function contains a Reidemeister I loop then after the last Gauss crossing has been added, the fringe contains a pair of edges having the same 
label for each Reidemeister I loop present in the diagram.  Note that these pairs may be nested, if the diagram contains a Reidemeister I loop within another Reidemeister I loop.
   
   
*/
bool gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data)
{

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code: initial gauss_code_data :" << endl;
	print_code_data(gauss_code_data,debug,"gauss_to_peer_code: ");	
}

	int num_components = gauss_code_data.num_components;
	int num_gauss_crossings = gauss_code_data.num_crossings;
	int num_gauss_arcs = 2*num_gauss_crossings;
	matrix<int>& gauss_code_table = gauss_code_data.code_table;

	for (int i=0; i< num_gauss_crossings; i++)
	{
	    if (gauss_code_table[LABEL][i] != generic_code_data::POSITIVE  && gauss_code_table[LABEL][i] != generic_code_data::NEGATIVE && gauss_code_table[LABEL][i] != generic_code_data::FLAT)
		{
if (braid_control::DEBUG >= braid_control::SUMMARY)
{
	debug << "gauss_to_peer_code: gauss_code_data includes an unsupported crossing label" << endl;
	print_code_data(gauss_code_data,debug,"gauss_to_peer_code: ");	
}
		    return false;
		}
	}
	
	/* Evaluate the planar diagram data for each crossing. Starting at the incoming under-arc, or the second visit to the
	   crossing in the case of flat crossings, list the Gauss arcs at the crossing as we move anti-clockwise around the 
	   crossing.  Note that the convention for odd and even terminating edges in Gauss code data is slightly counter-
	   intuitive, since we retain the peer code approach of the labels associated with the arcs through a crossing being 
	   stored consecutively as ET--OO and OT--EO (in peer codes the odd and even parity changes through a crossing) whereas
	   in Gauss codes "odd" and "even" refer to the first and second arrivals, so one might expect labels to be stored 
	   consecutively as ET--EO and OT--OO but this is not the case.     
	   
	   In a gauss code, the type of crossing is assigned based on the relative orientation of the first and second visits as 
	   follows:
		   
			1st \ /              2nd \ /         
				 X  =  Type 1         X  = Type 2    
			2nd / \              1st / \
		
       The odd and even originating and terminating edges are similarly assigned according to the first or second visit to the
       crossing, using the following convention:
		
			1st OT \ / OO             2nd ET \ / EO         
					X     =  Type 1           X  = Type 2    
			2nd ET / \ EO             1st OT / \ OO
	*/
	matrix<int> PD_data(num_gauss_crossings,4); // planar diagram data
	
	for (int i=0; i< num_gauss_crossings; i++)
	{
		if (gauss_code_table[LABEL][i] == generic_code_data::POSITIVE)
		{
			PD_data[i][0] = gauss_code_table[ODD_TERMINATING][i];

			if (gauss_code_table[TYPE][i] == generic_code_data::TYPE1)
			{
				PD_data[i][1] = gauss_code_table[EVEN_TERMINATING][i];
				PD_data[i][2] = gauss_code_table[EVEN_ORIGINATING][i];
				PD_data[i][3] = gauss_code_table[ODD_ORIGINATING][i];
			}
			else
			{
				PD_data[i][1] = gauss_code_table[ODD_ORIGINATING][i];
				PD_data[i][2] = gauss_code_table[EVEN_ORIGINATING][i];
				PD_data[i][3] = gauss_code_table[EVEN_TERMINATING][i];
			}			
		}
		else // generic_code_data::NEGAITIVE or generic_code_data::FLAT
		{
			PD_data[i][0] = gauss_code_table[EVEN_TERMINATING][i];

			if (gauss_code_table[TYPE][i] == generic_code_data::TYPE1)
			{
				PD_data[i][1] = gauss_code_table[EVEN_ORIGINATING][i];
				PD_data[i][2] = gauss_code_table[ODD_ORIGINATING][i];
				PD_data[i][3] = gauss_code_table[ODD_TERMINATING][i];
			}
			else
			{
				PD_data[i][1] = gauss_code_table[ODD_TERMINATING][i];
				PD_data[i][2] = gauss_code_table[ODD_ORIGINATING][i];
				PD_data[i][3] = gauss_code_table[EVEN_ORIGINATING][i];
			}			
		}
	}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code: PD_data:" << endl;
	print (PD_data,debug,3,"gauss_to_peer_code: ");
}



	vector<int> gauss_arc_direction;
	list<gc_pc_xlabels> virtual_crossings;
	list<int> diametric_virtual_crossings;
	vector<int> num_virtual_crossings_on_gauss_arc;
	int optimal_initial_crossing;
	int min_num_virtual_crossings = -1;
	bool too_many_virtual_crossings = false;
	
	/* initial_crossing loop starts here */	
	for (int initial_crossing=0; initial_crossing < num_gauss_crossings; initial_crossing++)
	{


if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: consider candidate initial_crossing = " << initial_crossing << endl;
	
		vector<int> candidate_num_virtual_crossings_on_gauss_arc (num_gauss_arcs); // counts the number of virtual crossings on each Gauss arc
	
		/* We associate a direction to each Gauss arc by tracing the diagram in the direction determined by the arc numbering.  For each arc
		   there is a unique crossing in the diagram at which the arc originates.  The direction assigned to an arc a is determined based on 
		   whether the originating crossing for arc a+1 is placed before after the orginating crossing for arc a, in the course of the algorithm.  
		   If the originating crossing for arc (a+1) is placed earlier, arc a is assigned the direction gc_pc_xlabels::direction::UPWARDS, if it 
		   placed afterwards, arc a is assigned the direction gc_pc_xlabels::direction::DOWNWARDS.  If an arc represents a Reidemeiseter I loop
		   then we assign it the direction gc_pc_xlabels::direction::SIDEWAYS.
		*/
		vector<int> candidate_gauss_arc_direction (num_gauss_arcs,gc_pc_xlabels::direction::UNKNOWN);
		
		list<gc_pc_xlabels> candidate_virtual_crossings;
		list<int> candidate_diametric_virtual_crossings;
		vector<int> fringe (2*num_gauss_crossings);
		
		
		
		for (int i=0; i< 4; i++)
			fringe[i] = PD_data[initial_crossing][i];
			
		int num_fringe_edges = 4;
	
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: num_fringe_edges = " << num_fringe_edges << " fringe-check: ";
	for (int f=0; f< num_fringe_edges; f++)
	{
		int temp = fringe[f]-1;
		if (temp <0)
			temp += num_gauss_arcs;
		debug << temp << ' ';
	}
	debug << endl;
}		


		/* The initialization of the fringe places the incoming under-arc top-left of the first crossing.
		   
		   If the first crossing includes a Reidemeister I loop then we shall then correct the direction of the 
		   loop arc to SIDEWAYS.  Similarly, the first crossing might involve a linked component that contains 
		   just one classical crossing, as in [-7 5 -9, -1, -3]/- + - * -.  In this latter case, the crossing
		   will involve a Gauss arc at locations 0 and 2, or 1 and 3 in the fringe.  This combination of locations
		   cannot occur with a Reidemeister I loop and so allows us to distinguish the two cases.  In the case of a
		   linked component, we add a virtual crossing to bring the two instances of the Gauss arc together.
		*/
		
		for (int i=0; i< 2; i++)
			assign_gauss_arc_direction(fringe[0+i],fringe[2+i],gauss_code_table,initial_crossing,candidate_gauss_arc_direction);
	
			
		/* check whether the first crossing includes a Reidemeister I loop or a linked component containing
	       one classical crossing.
		*/
		bool crossing_loop_detected = false;
		for (int i=0; i< 2 && !crossing_loop_detected; i++)
		{
			for (int j=i+1; j < 4; j++)
			{
				if (fringe[i] == fringe[j])
				{	
					candidate_gauss_arc_direction[fringe[i]] = gc_pc_xlabels::direction::SIDEWAYS;
					
					if ((i-j+4)%4 == 2)
					{
	                    add_virtual_crossing(i, fringe, candidate_num_virtual_crossings_on_gauss_arc, candidate_virtual_crossings);
						remove_fringe_edge(fringe[j],fringe); // use j because we've just added the virtual crossing at i.
					
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: linked component containing a single classical crossing detected at first crossing" << endl;
					}
					else
					{
						remove_fringe_edge(fringe[i],fringe);
					
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: Reidemeister I loop detected at first crossing" << endl;
					}
	
					crossing_loop_detected = true;				
					num_fringe_edges = 2;				
					break;
				}
			}
		}
		
		/* check VC count  - continue around initial_crossing loop */
		if (min_num_virtual_crossings != -1 && static_cast<int>(candidate_virtual_crossings.size()) > min_num_virtual_crossings)
		{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
		
			continue;
		}
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code: num_fringe_edges = " << num_fringe_edges << " current fringe: ";
	for (int f=0; f< num_fringe_edges; f++)
		debug << fringe[f] << ' ';
	debug << endl;
}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code: initial candidate_gauss_arc_direction: " ;
	for (int a=0; a< num_gauss_arcs; a++)
		debug << direction_to_string(a,candidate_gauss_arc_direction);
	debug << endl;
}		

		vector<bool> considered_crossing(num_gauss_crossings);
		considered_crossing[initial_crossing] = true;
		bool complete = false;
		
		do
		{
			/* look for another crossing */
			complete = true;
			int best_crossing;
			int best_match_count;
			int best_precedence = 5;  // bigger than the maximum prececence we use
			int last_fringe_index;
			int best_last_fringe_index;
			int first_fringe_index;
			int first_PD_index;
			int last_PD_index;
			int crossing_edge_c; // used for diametrically_connecting_edge_pair
			int crossing_edge_d;
			bool diametrically_connecting_edge_pair = false;
		
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: considered_crossing : ";
	for (int j=0; j< num_gauss_crossings; j++)
		debug << considered_crossing[j] << ' ';
	debug << endl;
}
		
			for (int i=0; i< num_gauss_crossings; i++)
			{			
				if (considered_crossing[i] == false)
				{
					
					complete = false;
					
					int match_count=0;
					last_fringe_index = -1;
					last_PD_index = 0;
					
					for (int j=0; j< 4; j++)
					{
						for (int k=0; k < num_fringe_edges; k++)
						{
							if (PD_data[i][j] == fringe[k])
							{
								match_count++;
								
								if (k > last_fringe_index)
								{
									last_fringe_index = k;
									last_PD_index = j;
								}
								
								break;
							}
						}
					}

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: crossing " << i << " is unused and matches " << match_count << " edges of the fringe" << endl;
	debug << "gauss_to_peer_code: last_fringe_index = " << last_fringe_index << ", last_PD_index = " << last_PD_index << endl;
}						
					if (match_count == 0)
						continue;  // there must be another crossing to consider in this case
						
					/* We assign a precedence to each crossing that has gauss edges in common with the fringe and select the next
					   crossing to be the one with the lowest precedence.  The precedence is assigned as follows:
					
					         0 shared vertices are contiguous in the fringe and in the correct order
					         1 shared vertices are contiguous and in the correct order but wrap around the start and end of the fringe
					         2 shared vertices are contiguous but virtual crossing(s) need to be added to obtain the correct order
					         3 shared edges are not contiguous, virtual crossings and wrapping of edges around the top of the 
							   diagram are required				        
					*/
					
					int precedence=0;
								    
					/*  Are the matching edges contiguous in the fringe and in the correct order?  If they are, then moving from right 
					    to left along the fringe the edges must appear in the corresponding order when moving anti-clockwise around 
					    the crossing.  That is, from the last_fringe_index moving backwards and the last_PD_index moving forwards.  
					    If this is not the case, increment the precedence to 1,
					*/				
		
					for (int j = 0; j < match_count; j++)
					{
						if (PD_data[i][(last_PD_index + j) % 4] != fringe[last_fringe_index - j])
						{
						
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   matching crossings not contiguous and in the correct order in the fringe, precedence incremented to 1" << endl;
	
							precedence = 1; 
						}
					}
									    
					if (precedence == 1) 
					{ 
						/* Check if there is a contiguous set of matching edges in the right order if we wrap edges from the start
						   of the fringe to the end.  This will only be the case if the last edge in the fringe is a matching
						   edge.
						*/ 
						if (last_fringe_index == num_fringe_edges - 1)
						{
						
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   last edge of the fringe is a matching edge, check for a wrapped contiguous set of matches in the correct order" << endl;
	
							/* identify the last matching edge at the start of the fringe, look for the first edge in the fringe 
							   that is not a matching edge and step back from that if it is not the very first edge.
							*/	
						
						    for (last_fringe_index = 0; last_fringe_index < num_fringe_edges; last_fringe_index++) 
						    {
								bool found = false;
								
							    for (int j = 0; j < 4; j++)
							    {
									if (PD_data[i][j] == fringe[last_fringe_index]) 
									{
										  found = true;
										  last_PD_index = j;
										  break;
									}
								}
								  
							    if (!found)
									break;
						    }
					    
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	if (last_fringe_index == num_fringe_edges)
		debug << "gauss_to_peer_code:   all fringe edges are matching edges" << endl;
	else
		debug << "gauss_to_peer_code:   first edge of the fringe that is not a matching edge is " << fringe[last_fringe_index] << " at index " <<  last_fringe_index << endl;
}
						    
						    
						    if (last_fringe_index != 0) // first matching current edge is not the first one, so we can step back
						    {
						        last_fringe_index--;
	
								/* look backwards along the fringe from the last edge we've just identified, wrapping around
								   to the end of the fringe, looking for a contiguous set of matching edges in the right order.
								   We have to look anti-clockwise around the crossing
								*/
						        for (int j = 0; j < match_count; j++)
						        {
									if (PD_data[i][(last_PD_index + j) % 4] != fringe[(last_fringe_index - j + num_fringe_edges) % num_fringe_edges])
									{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   wrapping edges does not produce a contiguous set of matching edges in the correct order, precedence incremented to 2" << endl;
	
										precedence = 2; // no wrap-around possible, so virtual crossing needed
									}
							    }
							    
							    if (precedence == 1)
							    {
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   contiguous matching edges in the correct order found wrapping to end of fringe" << endl;
								}					
						    } 
						    else 
						    {
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code:   first edge in the fringe is not a matching edge so wrapping edges doesn't help, precedence incremented to 2" << endl;
	debug << "gauss_to_peer_code:   reset last_fringe_index to the end of the fringe" << endl;
}
	
						        last_fringe_index = num_fringe_edges - 1;
						        precedence = 2;
						    }
						} 
						else
						{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   last edge in the fringe is not a matching edge, precedence incremented to 2" << endl;
	
							precedence = 2;
						}
					}
	
					if (precedence == 2) 
					{ 
						/* check whether the matching edges are at least contiguous, if not in the correct order */
						
						for (int j = 0; j < match_count; j++) 
						{						  
							bool found = false;
							
							for (int k = 0; k < 4; k++)
							{
								if (last_fringe_index - j >= 0 && PD_data[i][k] == fringe[last_fringe_index - j]) 
							    {
									found = true;
									break;
							    }
							}
							
							if (!found) 
							{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   matching edges not even contiguous in the fringe, precedence incremented to 3" << endl;
								
								precedence = 3; // other edges must be moved out of the way
							    break;
							}
						}
						
						if (precedence == 2)
						{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   matching edges are adjacent in the fringe, though not in the correct order" << endl;
						}
					}
				
					if (precedence < best_precedence || (precedence == best_precedence && match_count > best_match_count)) 
					{
					    best_crossing = i;
					    best_precedence = precedence;
					    best_match_count = match_count;
					    best_last_fringe_index = last_fringe_index;
					}								
				}
			}
			
			if (complete)
				break;
				
			
			last_fringe_index = best_last_fringe_index;
		
		
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: next best crossing is crossing " << best_crossing << ": ";
	for (int i=0; i< 4; i++)
		debug << PD_data[best_crossing][i] << ' ';
	debug << endl;
	debug << "gauss_to_peer_code: matching best_match_count = " << best_match_count << " edges, precedence = " << best_precedence << endl;
	debug << "gauss_to_peer_code: best_last_fringe_index = " << best_last_fringe_index << endl;
}

			/* If the best precedence is non-zero, adjust the fringe so that the matching edges are contiguous and in the correct 
			   order by adding virtual crossings and wrapping edges around the top of the diagram as necessary, ready to connect
			   the crossing.
			   
			   If best_precedence == 3, using a combination of addingvirtual crossings and wrapping edges from one end of the fringe
			   to the other, we adjust the matching edges so that they become contiguous, thereby reducing to precedence 0 or 2.  In 
			   the case we have reduced to precedence 0, we will have already adjusted the fringe accordingly so there is nothing more
			   to do.
			   
			   Therefore, we test for best_precedence == 1, then for best_precedence == 3 (possibly reducing to precedence 2) and finally 
			   check whether we have a precedence 2 situation to deal with.		
			   
			   In order to connect the crossing correctly, we need to know the first_PD_index around the crossing and the
			   first_fringe_index, which we can calculate from last_fringe_index, once the matching edges are contiguous in the fringe.
			   19/9/21 TO START WITH WE'LL JUST RE-CALCULATE THESE ONCE ALL THE FRINGE ADJUSTMENTS HAVE BEEN MADE.
			     
			*/
			
		    if (best_precedence == 1) 
		    {
			    /* We don't need to add any virtual crossings but do need to wrap edges from one side of the fringe to the other, 
			       so start by identifying which side has more matching edges 
			    */  
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: precedence == 1, just need to wrap matching edges" << endl;
			  
			    int left_count = 0;
			    int right_count = 0;
			      
			    for (int i = 0; i < num_fringe_edges; i++) 
			    {
					bool found = false;					
					for (int j = 0; j < 4; j++)
					{
						if (fringe[i] == PD_data[best_crossing][j]) 
						{
							left_count++;
						    found = true;
						    break;
						}
					}
						  
					if (!found)
						break;
			    }
	
			    for (int i = num_fringe_edges -1; i >= 0 ; i--) 
			    {
					bool found = false;					
					for (int j = 0; j < 4; j++)
					{
						if (fringe[i] == PD_data[best_crossing][j]) 
						{
							right_count++;
						    found = true;
						    break;
						}
					}
						  
					if (!found)
						break;
			    }
		      		      

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   " << left_count << " edges match on the left of the fringe " << right_count << " match on the right" << endl;
		      
				vector<int> new_fringe (2*num_gauss_crossings);
	
			    if (right_count < left_count) 
			    {	  			  
					for (int j = 0; j < right_count; j++)
						new_fringe[j] = fringe[num_fringe_edges - right_count + j];
						  
					for (int j=0; j < num_fringe_edges-right_count; j++)
						  new_fringe[right_count+j] = fringe[j];
			    } 
			    else 
			    {
					for (int j = 0; j < num_fringe_edges - left_count; j++)
						new_fringe[j] = fringe[j+left_count];
						  
					for (int j = 0; j < left_count; j++)
						new_fringe[num_fringe_edges - left_count + j] = fringe[j];
				}
				
				fringe = new_fringe;

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code:   new_fringe: ";
	for (int f=0; f< num_fringe_edges; f++)
		debug << fringe[f] << ' ';
	debug << endl;
}			
		    }
		    else if (best_precedence == 3) 
		    { 
				/* the matching edges in the fringe are are not grouped together, so we need to consider wrapping edges from 
				   one end of the fringe to the other and then add virtual crossings to make them contiguous, we can then proceed
				   as in the case of precedence 2
				*/
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: precedence == 3, matching edges need grouping together" << endl;
			
				/* identify the indices in edges of the matching edges */
			    vector<int> matching_fringe_index(best_match_count); // edgepos[be]; // which edges connect to the crossing
			    int index=0;
			    for (int i = 0; i < num_fringe_edges; i++)
			    {
					for (int j = 0; j < 4; j++)
					{
						if (fringe[i] == PD_data[best_crossing][j]) 
						{
							matching_fringe_index[index++] = i;
							break;
						}
					}
					if (index == best_match_count)
						break;
				}

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code:   matching_fringe_index: ";
	for (int i = 0; i < best_match_count; i++)
		debug << matching_fringe_index[i] << ' ';
	debug << endl;
}
			  		    
			    /* We look for the best "anchor" edge around which to build a contiguous set of matching edges in the fringe
			       by minimizing the number of virtual crossings required. For each matching edge i we calculate an array 
			       of flags 'wrap' that indicate for each matching edge j,
				   whether there are fewer edges between i and j within the fringe (directly between edges i and j), or
				   fewer edges moving from the first of i and j to the start of the fringe, and wrapping to the end.  If 
				   it is shorter to wrap, the flag is set.
				   
				   We also maintain, for each matching edge i, a cumulative count of the smallest number of intervening edges
				   within or wrapping the end of the fringe and select the edge with the smallest cumulative count as the
				   best edge, which we shall use as the anchor.
	
				   Note that, moving away from a matching edge i towards either end of the fringe, once we have encountered
				   a matching edge j, for which it is shorter to go the other way from i and wrap around the fringe, the
				   same will be true for all other matching edges encountered beyond edge j.  That means that, if the
				   vector wrap contains any set flags, they will always be flushed left and/or right within wrap, they 
				   will never take the form "1 0 1 0", or "0 1 0".
				   
				   Note also that from a given edge  we will not need to wrap "on both sides".  To see this, consider 
				   edges i, j and k appearing in the fringe at indices i,j,k so we have:
				   
				       0...j...i...k...(n-1)
	
				   where there are n edges in the fringe.  Then if n-(i-j) < (i-j), so that wrapping from i around the 
				   end of the fringe to get to j is shorter than going directly, we have 
				   n - (k-i) = (n-k) + (i-j) + j 
				             > (n-k) + (n- (i-j)) + j 
				             = (n-k) + ((n-k) + (k-i) +j) + j 
				             = 2(n-k) + 2j + (k-i)
				             > (k-i)
	
				    so the direct path from i to k is shorter.
				   
				   
			    */
			    int fringe_anchor_edge_index;
			    int min_virtual_count = -1;
			    vector<bool> matching_fringe_index_wrap_flag(best_match_count);
			    vector<bool> best_matching_fringe_index_wrap_flag(best_match_count);
			    		    
			    for (int i = 0; i < best_match_count; i++) 
			    {
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   evauate wrap_flags for matching_fringe_index " << i << ", edge " << fringe[matching_fringe_index[i]] << endl;;

					int count = 0;
					for (int j = 0; j < best_match_count; j++) 
					{
						int direct = abs(matching_fringe_index[j] - matching_fringe_index[i]);  // distance from i to j within the fringe
						int wrapped = num_fringe_edges - direct; // distance wrapping around the end of the fringe
						
						count += (direct < wrapped ? direct : wrapped);  // count the smallest number of edges that would need to be moved
						matching_fringe_index_wrap_flag[j] = (wrapped < direct); 

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code:     fringe index " << matching_fringe_index[j] << ", edge " << fringe[matching_fringe_index[j]]
	      << " direct count = " << direct << ", wrapped count = " << wrapped << ", cumuative count = " << count << endl;
}
					}
				
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code:     wrap_flags for matching fringe index " << matching_fringe_index[i] << ": ";
	for (int j=0; j< best_match_count; j++)
		debug << matching_fringe_index_wrap_flag[j] << ' ';
	debug << endl;
}						
					if (min_virtual_count < 0 || count < min_virtual_count) 
					{
						min_virtual_count = count;
						fringe_anchor_edge_index = i;
						best_matching_fringe_index_wrap_flag = matching_fringe_index_wrap_flag;

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     setting fringe_anchor_edge_index to " << fringe_anchor_edge_index << endl;
					}
			    }

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   final fringe_anchor_edge_index = " << fringe_anchor_edge_index << ", edge " << fringe[matching_fringe_index[fringe_anchor_edge_index]] << endl;
	
				/* Check the matching edges to the left of the anchor edge and see whether they require wrapping over the top of the 
				   diagram in order to minimize the number of virtual crossings required.  If the wrap flag is set for any of 
				   these matching edges, add virtual crossings to move them (sequentially) to the left of the fringe until the 
				   reach the left edge, or are adjacent to another matching edge.  This results in all of the matching edges
				   to the left of the anchor edge that need wrapping over the top of the diagram being flush left in the fringe.
				*/
				int left_count=0; // counts the number of matching edges we need to wrap on the left of the anchor edge
				
			    for (int i = 0; i < fringe_anchor_edge_index; i++) // first try on the left 
			    {
					if (best_matching_fringe_index_wrap_flag[i])
					{
						left_count++;
						
						while (matching_fringe_index[i] > 0 && 
						       (i > 0 ? matching_fringe_index[i] != matching_fringe_index[i-1]+1 : true)
						      ) 
						{
							add_virtual_crossing(matching_fringe_index[i]-1, fringe, candidate_num_virtual_crossings_on_gauss_arc, candidate_virtual_crossings);
							matching_fringe_index[i]--;
						}
					}
					else
						break;
				}
				
				/* check VC count  - break out of 'do' */
				if (min_num_virtual_crossings != -1 && static_cast<int>(candidate_virtual_crossings.size()) > min_num_virtual_crossings)
				{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
				
					too_many_virtual_crossings = true;
					break;
				}
		
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   number of matching edges that need wrapping to the left of the anchor edge = " << left_count << endl;
			
			    if (left_count > 0) 
			    {	  			 
					/* all the matching edges that need wrapping to the right are now flush-left in the fringe, 
					   so we need to wrap the first left_count edges from the left of the fringe to the right */
					vector<int> new_fringe (2*num_gauss_crossings);
														 
					for (int i = 0; i < num_fringe_edges - left_count; i++)
						new_fringe[i] = fringe[left_count+i];
						  
					for (int i = 0; i < left_count; i++)
						new_fringe[num_fringe_edges - left_count + i] = fringe[i];
						
					fringe = new_fringe;
					
					/* cycle the matching_fringe_index so that the first left_count indices appear at the end of the vector and
					   decrement each index by the left_count (modulo best_match_count)
					*/
					vector<int> new_matching_fringe_index(best_match_count);
					for (int i=0; i< best_match_count; i++)
					{
						new_matching_fringe_index[i] = matching_fringe_index[(left_count +i)%best_match_count];
						new_matching_fringe_index[i] = (new_matching_fringe_index[i] - left_count + num_fringe_edges)%num_fringe_edges;
					}
																	
					matching_fringe_index = new_matching_fringe_index;
					fringe_anchor_edge_index -= left_count;

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code:   new_fringe: ";
	for (int f=0; f< num_fringe_edges; f++)
		debug << fringe[f] << ' ';
	debug << endl;
	debug << "gauss_to_peer_code:   new_matching_fringe_index: ";
	for (int f=0; f< best_match_count; f++)
		debug << matching_fringe_index[f] << ' ';
	debug << endl;
	debug << "gauss_to_peer_code:   reset fringe_anchor_edge_index to " << fringe_anchor_edge_index << endl;	
}			

				}
			    else 
			    {	
					/* If we did not have to wrap edges to the left of bestedge, we may need to do so to the right.  We sequentially 
					   move matching edges that need wrapping to the right of the anchor edge by adding virtual crossings so that they 
					   end up flush right in the fringe.
					*/				
				
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   check to the right of the anchor edge = " << fringe[fringe_anchor_edge_index] << endl;
				
					int right_count=0; // counts the number of matching edges we need to wrap on the right of the anchor edge
					
				    for (int i = best_match_count -1; i > fringe_anchor_edge_index; i--)
				    {
						if (best_matching_fringe_index_wrap_flag[i])
						{
							right_count++;
						
/*if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code:     i = " << i << ", matching_fringe_index[i] = " << matching_fringe_index[i] << ", while test = " 
	      << (i < best_match_count - 1? matching_fringe_index[i] != matching_fringe_index[i+1]-1: true) << endl;
}
*/						
							while (matching_fringe_index[i] < num_fringe_edges -1 && 
							       (i < best_match_count - 1? matching_fringe_index[i] != matching_fringe_index[i+1]-1: true)
							      ) 
							{
								add_virtual_crossing(matching_fringe_index[i], fringe, candidate_num_virtual_crossings_on_gauss_arc, candidate_virtual_crossings);
								matching_fringe_index[i]++;
							}
						}
						else
							break;
					}
	
					/* check VC count  - break out of 'do' */
					if (min_num_virtual_crossings != -1 && static_cast<int>(candidate_virtual_crossings.size()) > min_num_virtual_crossings)
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
					
						too_many_virtual_crossings = true;
						break;
					}
							
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   number of matching edges that need wrapping to the right of the anchor edge = " << right_count << endl;
				
				    if (right_count > 0) 
				    {	  			 
						/* all the matching edges that need wrapping to the left are now flush-right in the fringe, 
						   so we need to wrap the last right_count edges from the right of the fringe to the left */

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code:   old_fringe: ";
	for (int f=0; f< num_fringe_edges; f++)
		debug << fringe[f] << ' ';
	debug << endl;
}					
						vector<int> new_fringe (2*num_gauss_crossings);
						
						for (int i = 0; i < right_count; i++)
							new_fringe[i] = fringe[num_fringe_edges - right_count + i];
							
						for (int i = 0; i < num_fringe_edges - right_count; i++)
							new_fringe[right_count+ i] = fringe[i];
							  						
						fringe = new_fringe;
						
						/* cycle the matching_fringe_index so that the first right_count indices appear at the front 
						   of the vector and increment each index by right_count modulo best_match_count
						*/
						vector<int> new_matching_fringe_index(best_match_count);
						for (int i=0; i< best_match_count; i++)
						{
							new_matching_fringe_index[(right_count +i)%best_match_count] = matching_fringe_index[i];
							new_matching_fringe_index[(right_count +i)%best_match_count] = (new_matching_fringe_index[(right_count +i)%best_match_count] + right_count)%num_fringe_edges;
						}
											
						matching_fringe_index = new_matching_fringe_index;
						fringe_anchor_edge_index = (fringe_anchor_edge_index+right_count)%best_match_count;

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code:   new_fringe: ";
	for (int f=0; f< num_fringe_edges; f++)
		debug << fringe[f] << ' ';
	debug << endl;
	debug << "gauss_to_peer_code:   new_matching_fringe_index: ";
	for (int f=0; f< best_match_count; f++)
		debug << matching_fringe_index[f] << ' ';
	debug << endl;
	debug << "gauss_to_peer_code:   reset fringe_anchor_edge_index to " << fringe_anchor_edge_index << endl;	
}		
					}				
			    }
			    
				/* Now all of the matching edges can be moved towards fringe[fringe_anchor_edge_index] within the fringe. */
	
			    for (int i = fringe_anchor_edge_index+1; i < best_match_count; i++)
			    {
					while (matching_fringe_index[i] != matching_fringe_index[i-1]+1)
					{
						add_virtual_crossing(matching_fringe_index[i]-1, fringe, candidate_num_virtual_crossings_on_gauss_arc, candidate_virtual_crossings);
						matching_fringe_index[i]--;
					}
				}
				
			    for (int i = fringe_anchor_edge_index-1; i >= 0; i--)
			    {
					while (matching_fringe_index[i] != matching_fringe_index[i+1]-1)
					{
						add_virtual_crossing(matching_fringe_index[i], fringe, candidate_num_virtual_crossings_on_gauss_arc, candidate_virtual_crossings);
						matching_fringe_index[i]++;
					}
				}			
				
				/* check VC count  - break out of 'do' */
				if (min_num_virtual_crossings != -1 && static_cast<int>(candidate_virtual_crossings.size()) > min_num_virtual_crossings)
				{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
				
					too_many_virtual_crossings = true;
					break;
				}
					
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: fringe after moving edges twowards anchor edge: ";
	for (int f=0; f< num_fringe_edges; f++)
		debug << fringe[f] << ' ';
	debug << endl;
}		
			    /* The matching edges are now adjacent but may be out of order, so check whether more virtual crossings are needed
				   edgepos[0] indicates the index in edges of the first matching edge.
				*/
				first_PD_index=-1;  //initialized to catch errors
	
			    for (int j = 0; j < 4; j++)
			    {
					if (PD_data[best_crossing][j] == fringe[matching_fringe_index[0]])
					{
						first_PD_index = j;
						break; 
					}
				}
	
				/* first_PD_index is the index around the crossing that matches the first shared edge on the fringe.  Therefore
				   we have to look backwards around the crossings as we move forwards along the fringe to check whether we have
				   a contiguous set of matching edges in the correct order.
				*/
			    for (int i = 0; i < best_match_count; i++)
			    {
					if (fringe[matching_fringe_index[i]] != PD_data[best_crossing][(first_PD_index -i + 4) % 4]) 
					{
						 best_precedence = 2; // matching edges out of order, more virtual crossings needed
						 last_fringe_index = matching_fringe_index[best_match_count-1];

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: new_fringe: matching edges out of order, more virtual crossings needed, set best_precedence to 2" << endl;
	debug << "gauss_to_peer_code: last_fringe_index = " << last_fringe_index << endl;
}

						break;
					}
				}
		    }
	
		    if (best_precedence == 2) 
			{
				/* the matching edges are contiguous but in the wrong order.  There must be either 2, 3 or 4 matching edges on the crossing, otherwise 
				   best_precedence would be 0.  If there are fewer than 4 matching edges, there is a unque order of the matching edges in the fringe 
				   that permits the non-matching edges of the crossing to become part of the new fringe without introducing additional virtual crossings.  
				   If best_match_count == 4, any order in the fringe would work but we choose the one minimizing the number of virtual crossings required.
				   
				   We identify the PD_data index of the first matching edge in the fringe in first_PD_index.  Since the PD_data records edges
				   anticlockwise around the crossing, this edge is the last of the matching edges around the crossing when moving anticlockwise.
				*/
	
				first_PD_index=-1;  //initialized to catch errors
				diametrically_connecting_edge_pair = false;
			
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: precedence == 2, the " << best_match_count << " matching edges in the fringe are contiguous but are in the wrong order" << endl;
		  
				if (best_match_count < 4)  // AB be = 2 or 3
				{
					/* identify the indices in PD_data of the edges that connect to the fringe as we move through the fringe from left to right */
					vector<int> PD_index_of_connecting_edge(best_match_count,-1); 				
				
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   last_fringe_index =  " << last_fringe_index << ", best_match_count = " << best_match_count << endl;
	
					for (int i = 0; i < best_match_count; i++)
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   look for fringe edge " << fringe[last_fringe_index - best_match_count + 1 + i] << endl;

						for (int j = 0; j < 4; j++)
						{
							if (fringe[last_fringe_index - best_match_count + 1 + i] == PD_data[best_crossing][j]) 
							{
								PD_index_of_connecting_edge[i] = j;
								break;
							}
						}
					}
				
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: PD_index_of_connecting_edge: ";
	for (int i = 0; i < best_match_count; i++)
		debug << PD_index_of_connecting_edge[i] << ' ';
	debug << endl;
}
					
					if (best_match_count == 2 && abs(PD_index_of_connecting_edge[1] - PD_index_of_connecting_edge[0]) == 2) 
					{ 
						/* special case where diametrically opposite edges connect to the fringe
						
						   If the fringe edges are a and b, in that order, then we can rotate the crossing so that edge b connects
						   directly to the fringe.  That leaves edge a requiring a virtual crossing, along with edge d, which is the 
						   edge around the crossing between b and a when moving anti-clockwise around the crossing.
						
						         a     b
								|   d   |
								 \ /-\ / b
								  X   /        <----- X is a virtual crossing
								 / \-/ \ c
								|   a   |					
								d       c
								
						   In this case, the fringe edges a and b are replaced by edges d and c respectively.
						*/
						first_PD_index = PD_index_of_connecting_edge[0];
						diametrically_connecting_edge_pair = true;
						candidate_diametric_virtual_crossings.push_back(candidate_virtual_crossings.size());
				  				  
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: matching edges are diametrically opposite on the crossing, first index of matching edges = " << first_PD_index << endl;
	
						crossing_edge_c = PD_data[best_crossing][(PD_index_of_connecting_edge[0]+1)%4];
						crossing_edge_d = PD_data[best_crossing][(PD_index_of_connecting_edge[1]+1)%4];  // as in above diagram
					
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   crossing_edge_c = " << crossing_edge_c << " crossing_edge_d = " << crossing_edge_d << endl;
					
					
						/* add the sideways arc in place zero to be consistent with other sideways scenarios */
					    gc_pc_xlabels xlabels;
					    xlabels.gauss.first = crossing_edge_d;
					    xlabels.gauss.second = fringe[last_fringe_index - 1];
					    
					    candidate_num_virtual_crossings_on_gauss_arc[xlabels.gauss.first]++;
					    candidate_num_virtual_crossings_on_gauss_arc[xlabels.gauss.second]++;
					    
					    candidate_virtual_crossings.push_back(xlabels);	
	
						/* check VC count  - break out of 'do' */
						if (min_num_virtual_crossings != -1 && static_cast<int>(candidate_virtual_crossings.size()) > min_num_virtual_crossings)
						{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
						
							too_many_virtual_crossings = true;
							break;
						}
								    
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   added virtual crossing " << crossing_edge_d << ' ' << fringe[last_fringe_index - 1] << endl;

					} 
					else 
					{
						/* look for an index in the PD_data corresponding to an edge that doesn't connect above.  The indices that
						   DO connect are recorded in PD_index_of_connecting_edge, so we look for an index that is not in that vector
						*/
						int non_matching_PD_index;
						
						for (int i = 0; i < 4; i++) 
						{
							bool found = false;
							
							for (int j = 0; j < best_match_count; j++)
							{
								if (PD_index_of_connecting_edge[j] == i)
								{
									found = true;
									break;
								}
							}
							
							if (!found)
							{
								non_matching_PD_index = i;
								break;
							}
						}
						  
						/* work clockwise around the crossing until we reach an edge which does connect, clockwise around the crossing 
						   corresponds to backwards along PD_data
						*/
						for (int i=0; i< 4; i++)
						{ 
							bool found = false;
							
							for (int j = 0; j < best_match_count; j++)
							{
								if (PD_index_of_connecting_edge[j] == (non_matching_PD_index - i + 4)%4 )
								{
									found = true;
									break;
								}
							}
	
							if (found)
							{
								first_PD_index = (non_matching_PD_index - i + 4)%4;
								break;
							}
	
						}
					
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: first matching edge in the fringe has PD_data index " << first_PD_index << endl;
					}
				} 
				else 
				{ 				
					/* Here best_match_count = 4, so we could, in principle, start connecting the crossing to the fringe starting at any of the edges in the 
					   crossing. We know that the matching edges in the fringe are contiguous and start at first_fringe_index, so we look to arrange the 
					   fringe edges from first_fringe_index onwards in the order corresponding to the crossing edges encountered as we move clockwise 
					   around the crossing.  We look for the best starting point on the crossing that introduces the fewest number of virtual crossings; i.e. 
					   the starting point that gives us the fewest out of order edges on the fringe, starting at first_fringe_index.
					*/			
				  	
					int best_starting_PD_index; 
					int least_number_of_virtual_crossings = 10;
					int	first_fringe_index = last_fringe_index - best_match_count + 1;
	
					
					for (int i = 0; i < 4; i++) 
					{ 
						/* consider starting at index i around the crossing and count the number of virtual crossings required */
						int virtual_crossing_count = 0;
						  
						for (int j = 0; j < 4; j++) // for each edge around the crossing
						{
							/* look forward from first_PD_index in the fringe - we know all four edges are contiguous, so only need to look forward 4 places. */
							for (int k = 0; k < 4; k++)  
							{
								if (PD_data[best_crossing][(i - j + 4) % 4] == fringe[first_fringe_index + k]) 
								{
									int num_virtual_crossings_required = abs(k-j);
									if (num_virtual_crossings_required > virtual_crossing_count)  
										virtual_crossing_count = num_virtual_crossings_required;
									break;
								}
							}
						}

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: starting index " << i << ", edge " << PD_data[best_crossing][i] << ", requires adding " << virtual_crossing_count << " virtual crossings" << endl;
					  
						if (virtual_crossing_count < least_number_of_virtual_crossings) 
						{
							best_starting_PD_index = i;
							least_number_of_virtual_crossings = virtual_crossing_count;
						}
					}
					
					first_PD_index = best_starting_PD_index; 

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: first matching edge in the fringe has PD_data index " << first_PD_index << endl;
								
				}
				
				if (!diametrically_connecting_edge_pair)
				{				   
					/* We want to add the crossing to the current fringe by considering the edges clockwise around the crossing from 
					   first_PD_index forwards along the fringe from the first fringe index.
					*/		
					first_fringe_index = last_fringe_index - best_match_count + 1;
				
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: first matching edge in the fringe has fringe index " << first_fringe_index << endl;
				
					/* work clockwise from first_PD_index on crossing, adding virtual crossings to the fringe where necessary */
					for (int i = 0; i < best_match_count; i++) 
					{
						for (int j = 0; j < best_match_count; j++)
						{
							/* look forwards along fringe for the next matching edge clockwise around the crossing */
							if (fringe[first_fringe_index + j] == PD_data[best_crossing][(first_PD_index - i + 4) % 4]) 
							{
								/*  If we're out of order, we had to look further along the fringe from the first_fringe_index than we have moved around 
								    the crossing from first_PD_index; that is,  first_fringe_index+j > first_fringe_index+i 
								*/
								for (int k = first_fringe_index + j - 1; k >= first_fringe_index + i; k--)  
									add_virtual_crossing(k, fringe, candidate_num_virtual_crossings_on_gauss_arc, candidate_virtual_crossings);
																
// AB: probably don't want to break here to handle Reidemeister I loops where the same edge appears twice around a crossing?
							}
						}
					}
	
					/* check VC count  - break out of 'do' */
					if (min_num_virtual_crossings != -1 && static_cast<int>(candidate_virtual_crossings.size()) > min_num_virtual_crossings)
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
					
						too_many_virtual_crossings = true;
						break;
					}
			
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: fringe after re-ordering: ";
	for (int f=0; f< num_fringe_edges; f++)
		debug << fringe[f] << ' ';
	debug << endl;
}									  
				}
		    }
	
		    /* add the next crossing to the fringe, first calculate the first_fringe_index and the first_PD_index */
		    first_fringe_index = num_fringe_edges; // NB this is the index *after* the end of the fringe
		    
		    for (int j=0; j< 4; j++)
			{
			
//if (braid_control::DEBUG >= braid_control::DETAIL)
//	debug << "gauss_to_peer_code: check for crossing index " << j << ", edge " << PD_data[best_crossing][j] << " in the fringe" << endl;   
			
				for (int k=0; k< num_fringe_edges; k++)
				{
					if (PD_data[best_crossing][j] == fringe[k])
					{
//if (braid_control::DEBUG >= braid_control::DETAIL)
//	debug << "gauss_to_peer_code:   found at fringe index " << k << endl;   
						if (k < first_fringe_index)
						{
							first_fringe_index = k;
							first_PD_index = j;
						}
						
						break;
					}
				}
			}
					
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: ready to add crossing to the fringe: first_fringe_index = " << first_fringe_index << " first_PD_index = " << first_PD_index << endl;   
	    
			int new_num_fringe_edges = num_fringe_edges + 4 - 2* best_match_count;
	    
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: when next crossing is added, new_num_fringe_edges  = " << new_num_fringe_edges << endl;
	    	    
		    if (diametrically_connecting_edge_pair)
		    {
			    fringe[first_fringe_index] = crossing_edge_d;
			    fringe[first_fringe_index+1] = crossing_edge_c;				    
			    
				if (crossing_edge_c == crossing_edge_d)
				{
					candidate_gauss_arc_direction[crossing_edge_c] = gc_pc_xlabels::direction::SIDEWAYS;			
					remove_fringe_edge(crossing_edge_c, fringe);
					num_fringe_edges -=2;
					
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: loop component containing a single classical crossing detected" << endl;
	
				}
				else
				{
					assign_gauss_arc_direction(crossing_edge_d,crossing_edge_c,gauss_code_table, best_crossing,candidate_gauss_arc_direction);		
				}
			
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: two new fringe edges: direction of edges " << crossing_edge_d << ',' << crossing_edge_c 
	      << " are " << direction_to_string(crossing_edge_d,candidate_gauss_arc_direction) << direction_to_string(crossing_edge_c,candidate_gauss_arc_direction) << endl;
}
			}
			else
			{	    
			    vector<int> new_fringe (new_num_fringe_edges);
			    int index = 0;
			    for (int i=0; i< first_fringe_index; i++)
					new_fringe[index++] = fringe[i];
					
			    for (int i=0; i< 4 - best_match_count; i++)
					new_fringe[index++] = PD_data[best_crossing][(first_PD_index+1+i)%4];
	
			    for (int i=0; i< new_num_fringe_edges-first_fringe_index-(4 - best_match_count); i++)
					new_fringe[index++] = fringe[first_fringe_index+best_match_count+i];
	
				
				/* assign the direction to the new edges in the fringe */
				if (best_match_count == 1)
				{
					int arc_1 = PD_data[best_crossing][(first_PD_index+1)%4];
					int arc_2 = PD_data[best_crossing][(first_PD_index+2)%4];
					int arc_3 = PD_data[best_crossing][(first_PD_index+3)%4];
	
					if (arc_2 == arc_1)
					{
						candidate_gauss_arc_direction[arc_2] = gc_pc_xlabels::direction::SIDEWAYS;
						candidate_gauss_arc_direction[arc_3] = candidate_gauss_arc_direction[fringe[first_fringe_index]];
						
						remove_fringe_edge(arc_2, new_fringe);
						new_num_fringe_edges -=2;
					
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: Reidemeister I loop detected at crossing" << endl;
	
					}
					else if (arc_2 == arc_3)
					{
						candidate_gauss_arc_direction[arc_1] = candidate_gauss_arc_direction[fringe[first_fringe_index]];
						candidate_gauss_arc_direction[arc_2] = gc_pc_xlabels::direction::SIDEWAYS;
						
						remove_fringe_edge(arc_2, new_fringe);
						new_num_fringe_edges -=2;

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: Reidemeister I loop detected at crossing" << endl;
	
					}
					else if (arc_1 == arc_3)
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: linked component containing a single classical crossing detected at crossing" << endl;
	
						candidate_gauss_arc_direction[arc_1] = gc_pc_xlabels::direction::SIDEWAYS;
						candidate_gauss_arc_direction[arc_2] = candidate_gauss_arc_direction[fringe[first_fringe_index]];
						
						/* We need to add a virtual crossing for the loop component but the virtual crossing does not involve fringe edges,
						   so we cannot call add_virtual_crossing.  We record the sideways arc in the first place and the fringe arc
						   in the second place.
					   */
						   
						gc_pc_xlabels xlabels;
						xlabels.gauss.first = arc_1;
						xlabels.gauss.second = arc_2;					
						candidate_num_virtual_crossings_on_gauss_arc[xlabels.gauss.first]++;
						candidate_num_virtual_crossings_on_gauss_arc[xlabels.gauss.second]++;					
						candidate_virtual_crossings.push_back(xlabels);					   					
	
						remove_fringe_edge(arc_1, new_fringe);
						new_num_fringe_edges -=2;	
						
						/* check VC count  - break out of 'do' */
						if (min_num_virtual_crossings != -1 && static_cast<int>(candidate_virtual_crossings.size()) > min_num_virtual_crossings)
						{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
						
							too_many_virtual_crossings = true;
							break;
						}
										
					}
					else
					{
						assign_gauss_arc_direction(arc_1,arc_3,gauss_code_table,best_crossing,candidate_gauss_arc_direction);						
						candidate_gauss_arc_direction[arc_2] = candidate_gauss_arc_direction[fringe[first_fringe_index]];
					}

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: three new fringe edges: direction of edges " << arc_1 << ',' << arc_2 << ',' << arc_3
	      << " are " << direction_to_string(arc_1,candidate_gauss_arc_direction) << direction_to_string(arc_2,candidate_gauss_arc_direction)
	      << direction_to_string(arc_3,candidate_gauss_arc_direction) << endl;
}
				    				
				}
				else if (best_match_count == 2)
				{
					int new_fringe_edge_1 = PD_data[best_crossing][(first_PD_index+1)%4];
					int new_fringe_edge_2 = PD_data[best_crossing][(first_PD_index+2)%4];
					
					int direction_1 = candidate_gauss_arc_direction[fringe[first_fringe_index]];
					int direction_2 = candidate_gauss_arc_direction[fringe[first_fringe_index+1]];
					
					if (direction_1 == gc_pc_xlabels::direction::UPWARDS)
					{
						if (direction_2 == gc_pc_xlabels::direction::UPWARDS)
						{
							candidate_gauss_arc_direction[new_fringe_edge_1] = gc_pc_xlabels::direction::UPWARDS;
							candidate_gauss_arc_direction[new_fringe_edge_2] = gc_pc_xlabels::direction::UPWARDS;												
						}
						else
						{
							if (new_fringe_edge_1 == new_fringe_edge_2)
							{
								candidate_gauss_arc_direction[new_fringe_edge_1] = gc_pc_xlabels::direction::SIDEWAYS;
								remove_fringe_edge(new_fringe_edge_1, new_fringe);
								new_num_fringe_edges -=2;

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: Reidemeister I loop detected at crossing" << endl;
							}
							else
							{
								candidate_gauss_arc_direction[new_fringe_edge_1] = gc_pc_xlabels::direction::DOWNWARDS;
								candidate_gauss_arc_direction[new_fringe_edge_2] = gc_pc_xlabels::direction::UPWARDS;
							}
						}
					}
					else
					{
						if (direction_2 == gc_pc_xlabels::direction::UPWARDS)
						{
							if (new_fringe_edge_1 == new_fringe_edge_2)
							{
								candidate_gauss_arc_direction[new_fringe_edge_1] = gc_pc_xlabels::direction::SIDEWAYS;
								remove_fringe_edge(new_fringe_edge_1, new_fringe);
								new_num_fringe_edges -=2;
							
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: Reidemeister I loop detected at crossing" << endl;
							}
							else
							{
								candidate_gauss_arc_direction[new_fringe_edge_1] = gc_pc_xlabels::direction::UPWARDS;
								candidate_gauss_arc_direction[new_fringe_edge_2] = gc_pc_xlabels::direction::DOWNWARDS;
							}
						}
						else
						{
							candidate_gauss_arc_direction[new_fringe_edge_1] = gc_pc_xlabels::direction::DOWNWARDS;
							candidate_gauss_arc_direction[new_fringe_edge_2] = gc_pc_xlabels::direction::DOWNWARDS;
						}
					}				
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: two new fringe edges: direction of edges " << new_fringe_edge_1 << ',' << new_fringe_edge_2 
	      << " are " << direction_to_string(new_fringe_edge_1,candidate_gauss_arc_direction) << direction_to_string(new_fringe_edge_2,candidate_gauss_arc_direction) << endl;
}
				}
				else if (best_match_count == 3)
				{
					int up_count = 0;
					int down_count = 0;
	
			        for (int i=0; i< 3; i++)
			        {
						if (candidate_gauss_arc_direction[fringe[first_fringe_index+i]] == gc_pc_xlabels::direction::UPWARDS)
							up_count++;
						else 
							down_count++;
					}
					
					int new_fringe_edge = PD_data[best_crossing][(first_PD_index+1)%4];
					
					if (up_count == 2)
						candidate_gauss_arc_direction[new_fringe_edge] = gc_pc_xlabels::direction::UPWARDS;
					else
						candidate_gauss_arc_direction[new_fringe_edge] = gc_pc_xlabels::direction::DOWNWARDS;

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: one new fringe edge: direction of edge " << new_fringe_edge << " is " << direction_to_string(new_fringe_edge,candidate_gauss_arc_direction) << endl;
}
				}
							
				fringe = new_fringe;
				num_fringe_edges = new_num_fringe_edges;
			}	

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code: num_fringe_edges = " << num_fringe_edges << " fringe-check: ";
	for (int f=0; f< num_fringe_edges; f++)
	{
		int temp = fringe[f]-1;
		if (temp <0)
			temp += num_gauss_arcs;
		debug << temp << ' ';
	}
	debug << endl;
}		

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: fringe after adding new crossing: ";
	for (int f=0; f< num_fringe_edges; f++)
		debug << fringe[f] << ' ';
	debug << endl;
	
	debug << "gauss_to_peer_code: num_fringe_edges adjusted to: " << num_fringe_edges << endl;
}									  
		
			considered_crossing[best_crossing] = true;
	
		} while (!complete);
	
		if (too_many_virtual_crossings)
			continue; // around the initial_crossing loop

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: total number of virtual crossings added = " << candidate_virtual_crossings.size() << endl;

	if (candidate_virtual_crossings.size() != 0)
	{
		debug << "gauss_to_peer_code: candidate_virtual_crossings Gauss labels: " << endl;
		list<gc_pc_xlabels>::iterator lptr = candidate_virtual_crossings.begin();
		
		while (lptr != candidate_virtual_crossings.end())
		{
			debug << "gauss_to_peer_code: " << lptr->gauss.first << ' ' << lptr->gauss.second << endl;
			lptr++;
		}
	}
	
	debug << "gauss_to_peer_code: candidate_num_virtual_crossings_on_gauss_arc: ";
	for (int i=0; i< num_gauss_arcs; i++)
		debug << candidate_num_virtual_crossings_on_gauss_arc[i] << ' ';
	debug << endl;

	debug << "gauss_to_peer_code: final candidate_gauss_arc_direction: " ;
	for (int a=0; a< num_gauss_arcs; a++)
		debug << direction_to_string(a,candidate_gauss_arc_direction);
	debug << endl;	
}									  




		if (min_num_virtual_crossings == -1 || static_cast<int>(candidate_virtual_crossings.size()) < min_num_virtual_crossings)
		{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: initial crossing " << initial_crossing << " requires minimal number of virtual crossings" << endl;

			optimal_initial_crossing = initial_crossing;
			min_num_virtual_crossings = candidate_virtual_crossings.size();
			gauss_arc_direction = candidate_gauss_arc_direction;
			virtual_crossings = candidate_virtual_crossings;
			diametric_virtual_crossings = candidate_diametric_virtual_crossings;
			num_virtual_crossings_on_gauss_arc = candidate_num_virtual_crossings_on_gauss_arc;
		}
	} /* initial_crossing loop ends here */

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: minimum number of virtual crossings is " << virtual_crossings.size() << ", using initial_crossing " << optimal_initial_crossing<< endl;

	if (virtual_crossings.size() != 0)
	{
		debug << "gauss_to_peer_code: virtual_crossings Gauss labels: " << endl;
		list<gc_pc_xlabels>::iterator lptr = virtual_crossings.begin();
		
		while (lptr != virtual_crossings.end())
		{
			debug << "gauss_to_peer_code: " << lptr->gauss.first << ' ' << lptr->gauss.second << endl;
			lptr++;
		}
	}

	if (diametric_virtual_crossings.size() != 0)
	{
		debug << "gauss_to_peer_code: diametric_virtual_crossings ";
		list<int>::iterator lptr = diametric_virtual_crossings.begin();
		
		while (lptr != diametric_virtual_crossings.end())
		{
			debug << "gauss_to_peer_code: " << *lptr << ' ';
			lptr++;
		}
		
		debug << endl;
	}
	
	debug << "gauss_to_peer_code: num_virtual_crossings_on_gauss_arc: ";
	for (int i=0; i< num_gauss_arcs; i++)
		debug << num_virtual_crossings_on_gauss_arc[i] << ' ';
	debug << endl;

	debug << "gauss_to_peer_code: gauss_arc_direction: " ;
	for (int a=0; a< num_gauss_arcs; a++)
		debug << direction_to_string(a,gauss_arc_direction);
	debug << endl;	
}									  

	

	int num_immersion_crossings = num_gauss_crossings+virtual_crossings.size();

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: num_immersion_crossings = " << num_immersion_crossings << endl;


	/* write the diametrically opposite virtual crossings into a flag vector */
	vector<bool> diametric_virtual_crossing(virtual_crossings.size());
	
	if (diametric_virtual_crossings.size() != 0)
	{
		list<int>::iterator lptr = diametric_virtual_crossings.begin();		
		while (lptr != diametric_virtual_crossings.end())
		{
			diametric_virtual_crossing[*lptr] = true;;
			lptr++;
		}
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: diametric_virtual_crossing: ";
	for (unsigned int i=0; i< virtual_crossings.size(); i++)
		debug << diametric_virtual_crossing[i] << ' ';
	debug << endl;
}		
	}

	int next_immersion_arc = 0;
	matrix <int> immersion_crossing_peers(num_immersion_crossings,2); // the first num_gauss_crossing rows will record the Gauss crossing peers
	vector<int> first_immersion_edge_on_component(num_components);
	int component = 0;
	first_immersion_edge_on_component[component] = next_immersion_arc;
	

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "gauss_to_peer_code: calculate immersion edge labels" << endl;
	
	for (int i=0; i< num_gauss_arcs; i++)
	{
		if (i == gauss_code_data.first_edge_on_component[component]+gauss_code_data.num_component_edges[component])
		{
			component++;
			first_immersion_edge_on_component[component] = next_immersion_arc;
			
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "gauss_to_peer_code:   found start of component " << component << ", first immersion edge =  " << first_immersion_edge_on_component[component] << endl;
	
		}
		
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "gauss_to_peer_code:   gauss_arc " << i << ", contains " << num_virtual_crossings_on_gauss_arc[i] << " virtual crossings" << endl;

		if (num_virtual_crossings_on_gauss_arc[i] !=0)
		{		
			if (gauss_arc_direction[i] == gc_pc_xlabels::direction::DOWNWARDS)
			{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "gauss_to_peer_code:     DOWNWARDS arc" << endl;
	
				list<gc_pc_xlabels>::iterator lptr = virtual_crossings.begin();
				
				while (lptr != virtual_crossings.end())
				{
					if (lptr->gauss.first == i)
					{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code:     found arc in first place at virtual crossing " << lptr->gauss.first << ',' << lptr->gauss.second
	      << ", allocating immersion edge label " << next_immersion_arc << endl;
}	
						lptr->immersion.first = next_immersion_arc++;
					}
					else if (lptr->gauss.second == i)
					{

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code:     found arc in second place at virtual crossing " << lptr->gauss.first << ',' << lptr->gauss.second
	      << ", allocating immersion edge label " << next_immersion_arc << endl;
}
						lptr->immersion.second = next_immersion_arc++;
					}
					
					lptr++;
				}
			}
			else
			{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "gauss_to_peer_code:     UPWARDS or SIDEWAYS arc" << endl;
	
				list<gc_pc_xlabels>::reverse_iterator lptr = virtual_crossings.rbegin();
				
				while (lptr != virtual_crossings.rend())
				{
					if (lptr->gauss.first == i)
					{						
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code:     found arc in first place at virtual crossing " << lptr->gauss.first << ',' << lptr->gauss.second
	      << ", allocating immersion edge label " << next_immersion_arc << endl;
}
						lptr->immersion.first = next_immersion_arc++;
					}
					else if (lptr->gauss.second == i)
					{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code:     found arc in second place at virtual crossing " << lptr->gauss.first << ',' << lptr->gauss.second
	      << ", allocating immersion edge label " << next_immersion_arc << endl;
}						
						lptr->immersion.second = next_immersion_arc++;
					}
					
					lptr++;
				}
			}
			
		}

		for (int j=0; j< num_gauss_crossings; j++)
		{
			if(gauss_code_table[EVEN_TERMINATING][j] == i)
			{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "gauss_to_peer_code:   found arc as EVEN_TERMINATING at Gauss crossing " << j << ", allocating immersion edge label " << next_immersion_arc << endl;
	
				immersion_crossing_peers[j][0] = next_immersion_arc++;
				break;
			}
			else if(gauss_code_table[ODD_TERMINATING][j] == i)
			{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "gauss_to_peer_code:   found arc as ODD_TERMINATING at Gauss crossing " << j << ", allocating immersion edge label " << next_immersion_arc << endl;
	
				immersion_crossing_peers[j][1] = next_immersion_arc++;
				break;
			}
		}
	}
	
	/* write the virtual crossing immersion labels into immersion_crossing_peers and the gauss labels into a similar matrix, 
	   virtual_gauss_peers, so we can index into virtual_gauss_peers for the corresponding Gauss peers when we come to evaluate the
	   type.
	*/
	matrix <int> virtual_gauss_peers(num_immersion_crossings,2,-1); // the first num_gauss_crossing rows will be unused


	int index = num_gauss_crossings;
	list<gc_pc_xlabels>::iterator lptr = virtual_crossings.begin();
		
	while (lptr != virtual_crossings.end())
	{
		virtual_gauss_peers[index][0] = lptr->gauss.first;
		virtual_gauss_peers[index][1] = lptr->gauss.second;
		immersion_crossing_peers[index][0] = lptr->immersion.first;
		immersion_crossing_peers[index][1] = lptr->immersion.second;
		index++;
		lptr++;
	}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code: immersion_crossing_peers: " << endl;
	print(immersion_crossing_peers, debug, 4, "gauss_to_peer_code: ");
}	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code: virtual_gauss_peers: " << endl;
	print(virtual_gauss_peers, debug, 4, "gauss_to_peer_code: ");
}	
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "gauss_to_peer_code: next_immersion_arc  = " << next_immersion_arc << " (should be twice the number of immersion crossings)" << endl;
	
	/* evaluate the number of immersion edges on each component */
	vector<int> num_immersion_component_edges(num_components);
	for (int i=0; i<num_components; i++)
	{
		int next_first_edge = (i< num_components - 1? first_immersion_edge_on_component[i+1]: next_immersion_arc);
		num_immersion_component_edges[i] = next_first_edge - first_immersion_edge_on_component[i];
	}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code: num_immersion_component_edges: ";
	for (int i=0; i< num_components; i++)
		debug << num_immersion_component_edges[i] << ' ';
	debug << endl;
	debug << "gauss_to_peer_code: first_immersion_edge_on_component: ";
	for (int i=0; i< num_components; i++)
		debug << first_immersion_edge_on_component[i] << ' ';
	debug << endl;	
}
	/* There is the possibility that we have not ended up with an odd and even terminating (immersion) edge at each crossing,
	   so we may need to cycle the edge labels for some components.
	   
	   To ensure that we don't enter an infinite loop of shifting component numbering, we shall consider each component in
	   the order in which we first encounter them, rather than the order determined by the numbering of the gauss code.   
	*/		
	if (num_components > 1)
	{
		list<int> next_component;
		next_component.push_back(0);
		vector <bool> component_considered(num_components);
		
		bool complete = false; // only needed for debugging
		
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: check components for aligned odd and even parity edges" << endl;
		do
		{
			if (next_component.size() == 0)
			{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   stopped processing components due to empty list" << endl;
				complete = true;
				break;
			}
				
			int component = next_component.front();
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   front of next_component = " << component << endl;
		
			next_component.pop_front();
			if (component_considered[component])
				continue;
			else
				component_considered[component] = true;
			
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   processing component " << component << endl;
			
			/* work around this component */
			for (int i=0; i< num_immersion_component_edges[component]; i++)
			{
				int edge = first_immersion_edge_on_component[component]+i;
				int peer_edge;
			
				/* find edge in immersion_crossing_peers */
				int peer_component;
				
				for (int j=0; j< num_immersion_crossings; j++)
				{
					if (immersion_crossing_peers[j][0] == edge)
					{
						peer_edge = immersion_crossing_peers[j][1];
						if (peer_edge % 2 == edge % 2)
						{
							/* both the first and second visits to this crossing have the same parity, so cycle the
							   edges in the component containing peer_edge to align the parity with edge */
	
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     immersion_crossing_peers shows the peer of " << edge << " is " << peer_edge 
	      << ", cycle component containing " << peer_edge << endl;
	                        cycle_gauss_code_labels(immersion_crossing_peers,num_immersion_component_edges,first_immersion_edge_on_component,num_immersion_crossings,num_components,peer_edge);
						}
						
						peer_component = 0;
						for (int j=1; j< num_components; j++)
						{
							if (peer_edge >= first_immersion_edge_on_component[j])
								peer_component++;
							else
								break;
						}
						next_component.push_back(peer_component);
						
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:       peer component " << peer_component << endl;
										
						break;
					}
					else if (immersion_crossing_peers[j][1] == edge)
					{
						peer_edge = immersion_crossing_peers[j][0];
						if (peer_edge % 2 == edge % 2)
						{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     immersion_crossing_peers shows the peer of " << edge << " is " << peer_edge
	      << ", cycle component containing " << peer_edge << endl;
		
	                        cycle_gauss_code_labels(immersion_crossing_peers,num_immersion_component_edges,first_immersion_edge_on_component,num_immersion_crossings,num_components,peer_edge);
						}
						
						peer_component = 0;
						for (int j=1; j< num_components; j++)
						{
							if (peer_edge >= first_immersion_edge_on_component[j])
								peer_component++;
							else
								break;
						}
						next_component.push_back(peer_component);
						
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:       peer component " << peer_component << endl;
						
						break;
					}
				}
			}    
		} while (!complete);
	}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code: final immersion_crossing_peers: " << endl;
	print(immersion_crossing_peers, debug, 4, "gauss_to_peer_code: ");
}	

	matrix<int> peer_code_table(CODE_TABLE_SIZE,num_immersion_crossings);

	/* write the COMPONENT data into a code table for the peer code */
	for (int i=0; i< num_immersion_crossings; i++)
	{
		int component = 0;
		for (int j=1; j< num_components; j++)
		{
			if (2*i >= first_immersion_edge_on_component[j])
				component++;
			else
				break;
		}
		peer_code_table[COMPONENT][i] = component;
	}
	
	/* write the ODD and EVEN peers to peer_code_table */
	int even_edge;
	int odd_edge;
	
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: write immersion labels:" << endl;
	
	for (int i=0; i< num_immersion_crossings; i++)
	{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   immersion_crossing_peers row " << i;
	
		if (immersion_crossing_peers[i][0] % 2 == 0)
		{
			even_edge = immersion_crossing_peers[i][0];
			odd_edge = immersion_crossing_peers[i][1];
		}	
		else
		{
			even_edge = immersion_crossing_peers[i][1];
			odd_edge = immersion_crossing_peers[i][0];
		}	
		
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << " even immersion edge = " << even_edge << ", odd immersion edge = " << odd_edge << endl;

		peer_code_table[OPEER][even_edge/2] = odd_edge;
		peer_code_table[EPEER][(odd_edge-1)/2] = even_edge;
	}

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code: imersion EPEERs: ";
	for (int i=0; i< num_immersion_crossings; i++)
		debug << peer_code_table[EPEER][i] << ' ';
	debug << endl;
	debug << "gauss_to_peer_code: imersion OPEERs: ";
	for (int i=0; i< num_immersion_crossings; i++)
		debug << peer_code_table[OPEER][i] << ' ';
	debug << endl;
}				

	/* write the crossing TYPE data to peer_code_table. In a gauss code, the type of crossing is assigned based on the relative 
	   orientation of the first and second visits as follows:
		   
		1st \ /              2nd \ /         
			 X  =  Type 1         X  = Type 2    
		2nd / \              1st / \

       The peer code type of a Gauss crossings may therefore be determined from the parity of the immersion label assigned to the second visit.
       If the 2nd visit immersion label is even, the peer code type matches the Gauss code type, otherwise it is reversed.  Note that the first 
       and second visits are recorded in the rows of immersion_crossing_peers in places 1 and 0 respectively.
      
       For virtual crossings that are not added as a result of diametrically opposite matching fringe edges, we recorded the Gauss and immersion 
       labels in the rows of immersion_crossing_peers and virtual_gauss_peers in places 0 and 1 in the order that the Gauss arcs appeared in the 
       fringe from left to right.  Note however, that the immersion labels respected the  direction of the Gauss arc, as shown below.

        1st   2nd                       2nd    1st      
	                   ^   ^        ^                ^
	       \ /          \ /          \ /          \ /
			X            X            X            X
	       / \          / \          / \          / \
          v   v      2nd   1st      v   1st    2nd   v
          D   D        U   U        D   U        U   D
          
       If the Gauss arcs are both downwards, the first visit immersion label is the lower terminating immersion label; if they are both upwards, the
       first visit immersion label is the lower terminating immersion label.  If the first Gauss arc in the fringe is upwards and the second downwards, 
       the first visit immersion is the upper terminating immersion label and finally, if the first Gauss arc in the fringe is downwards and the second 
       upwards, the first visit immersion is again the upper terminating immersion label.
       
       If the virtual crossing was added as a result of diametrically opposite matching fringe edges, or if we have a linked component that contains just 
       one classical crossing, then regardless of whether the classical crossing were the first crossing considered (whereupon add_virtual_crossing is called) 
       or not (whereupon the virtual crossing was added manually, since it does not involve fringe edges), the Gauss and immersion labels are recorded with 
       the arc dorresponding to crossing_edge_d, or the sideways loop in place zero.
	*/
	for (int i = 0; i< num_immersion_crossings; i++)
	{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   " << (i < num_gauss_crossings? "Gauss": "virtual") << " crossing " << (i < num_gauss_crossings? i: i-num_gauss_crossings);
	
		int even_peer;
		bool even_peer_in_place_zero;
		
		if (immersion_crossing_peers[i][0] % 2 == 0)
		{
			even_peer = immersion_crossing_peers[i][0];
			even_peer_in_place_zero = true;
		}
		else
		{
			even_peer = immersion_crossing_peers[i][1];
			even_peer_in_place_zero = false;
		}
			
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << ", even_peer = " << even_peer << ", even_peer_in_place_zero (of immersion_crossing_peers) = " << even_peer_in_place_zero << endl;

		if (i < num_gauss_crossings)
		{
			if (even_peer_in_place_zero)
			{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to second visit to Gauss crossing " << i;
	
				peer_code_table[TYPE][even_peer/2] = gauss_code_table[TYPE][i];
			}
			else // even label corresponds to first visit to Gauss crossing
			{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to first visit to Gauss crossing " << i;

				if (gauss_code_table[TYPE][i] == generic_code_data::type::TYPE1)
					peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE2;
				else 
					peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE1;
			}
			
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << " crossing is (peer) " << (peer_code_table[TYPE][even_peer/2] == generic_code_data::type::TYPE1? "TYPE1": "TYPE2") << endl;
		}
		else
		{
			int first_fringe_edge_direction = gauss_arc_direction[virtual_gauss_peers[i][0]];
			int second_fringe_edge_direction = gauss_arc_direction[virtual_gauss_peers[i][1]];

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code:     first_fringe_edge_direction " ;
	debug << direction_to_string(virtual_gauss_peers[i][0],gauss_arc_direction);
	debug << ", second_fringe_edge_direction " ;
	debug << direction_to_string(virtual_gauss_peers[i][1],gauss_arc_direction);
	debug << endl;
}					

			if (diametric_virtual_crossing[i-num_gauss_crossings])
			{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     virtual crossing " << i-num_gauss_crossings << " added as a result of a diametrically opposite matching Gauss crossing" << endl;
	
				/* If the virtual crossing was introduced as a result of added a Gauss crossing with diametrically opposite matching edges, then it involves the first
				   matching fringe edge of the Gauss crossing and crossing_edge_d, as int he above diagram.  However, when we stored the immersion_crossing_peers and
				   virtual_gauss_peers we stored the labels and arcs of crossing_edge_d in position 0 to be consistent with sideways loops.  Here both the first fringe 
				   edge and crossing_edge_d have been assigned UPWARDS or DOWNWARDS directions, so we can determine the type by considering which of the two edges is
				   assigned the even immersion label together with the two directions, as shown below.
				
                                                                  e                         e
                          ^           ^             |             |             |           |             ^             ^
				          |           |             |             |             |           |             |             |
				     e ---o--> U   ---o--> U   D <--o--- e   D <--o---     e ---o--> U   ---o--> U   D <--o--- e   D <--o---     <<== first_fringe_edge_direction
				          |           |             |             |             |           |             |             | 
				          |           |             |             |             |           |             |             |
                                      e             v             v             v           v                           e
 
   				         II           I            II             I             I           II            I             II
				
				   In the diagram the horizontal arrow corresponds to crossing_edge_d and thus the first_fringe_edge_direction as noted above.  The
				   vertical arrow is actually the first fringe edge but because of the convention of storing crossing_edge_d in position zero, it is
				   regarded as the second_fringe_edge_direction here.  Thus,       
				   
				   first_fringe_edge_direction is upwards and second_fringe_edge_direction is upwards; or
                   first_fringe_edge_direction is downwards and second_fringe_edge_direction is downwards
                   
						even_peer_in_place_zero => TYPE II   !even_peer_in_place_zero => TYPE I
				         
				   first_fringe_edge_direction is upwards and second_fringe_edge_direction is downwards; or
                   first_fringe_edge_direction is downwards and second_fringe_edge_direction is upwards
                   
						even_peer_in_place_zero => TYPE I   !even_peer_in_place_zero => TYPE II
				   				
				*/
				if ((first_fringe_edge_direction == gc_pc_xlabels::direction::UPWARDS && second_fringe_edge_direction == gc_pc_xlabels::direction::UPWARDS) ||
				    (first_fringe_edge_direction == gc_pc_xlabels::direction::DOWNWARDS && second_fringe_edge_direction == gc_pc_xlabels::direction::DOWNWARDS)
			       )
			    {
					
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     both Gauss arcs have the same direction at virtual crossing " << i-num_gauss_crossings << endl;
	
					if (even_peer_in_place_zero)
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to non-matching fringe Gauss arc";
						peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE2;				
					}
					else
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to first matching fringe Gauss arc,";
	
						peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE1;				
					}
				}
				else
				{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     Gauss arcs have opposite direction at virtual crossing " << i-num_gauss_crossings << endl;

					if (even_peer_in_place_zero)
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to non-matching fringe Gauss arc";
						peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE1;				
					}
					else
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to first matching fringe Gauss arc,";
						peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE2;				
					}
				}
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << " crossing is (peer) " << (peer_code_table[TYPE][even_peer/2] == generic_code_data::type::TYPE1? "TYPE1": "TYPE2") << endl;
			}
			else if (first_fringe_edge_direction == gc_pc_xlabels::direction::SIDEWAYS)
			{
				/* If the first Gauss arc involved in a virtual crossings is a sideways arc, then it is a loop component 
				   that contains a single Gauss crossing.  In such a case we cannot identify whether the virtual crossing 
				   should be TYPE1 or TYPE2 based on the first or second visit, since the orientation of the loop component 
				   is not uniquely determined by that information.
				   
				   Instead, we have to look at the Gauss crossing on the loop component to determine the orientation of the
				   loop component.  Once we have the orientation, even_peer_in_place_zero tells us the immersion type of the 
				   virtual crossing.  The loop component is the first fringe edge, and we can determine whether the sideways 
				   loop is the second visit, so we have
				   
				          ^         ^                ^         ^                |         |                |         |
				          |         |                |         |                |         |                |         |
				   2nd ---o-->   <--o--- 2nd  1st ---o-->   <--o--- 1st  2nd ---o-->   <--o--- 2nd  1st ---o-->   <--o--- 1st 
				          |         |                |         |                |         |                |         | 
				          |         |                |         |                v         v                v         v
				          
				         II         I                I         II               I         II              II         I
				         
				   where the 'o' indicates a Gauss crossing of some kind and the horizontal arcs loop around either upwards or 
				   downwards depending on the direction of the non-sideways arc and whether the Gauss crossing preceeeds the 
				   virtual crossing introduced by the loop or not.  Thus:
				                                                                          

				   sideways loop is second visit && second_fringe_edge_direction is UPWARDS && !gauss_crossing_preceeds_virtual_crossing; or
				   sideways loop is first visit && second_fringe_edge_direction is DOWNWARDS && gauss_crossing_preceeds_virtual_crossing; or
				   sideways loop is second visit && second_fringe_edge_direction is DOWNWARDS && !gauss_crossing_preceeds_virtual_crossing; or
				   sideways loop is first visit && second_fringe_edge_direction is UPWARDS && gauss_crossing_preceeds_virtual_crossing;
				   
				       TYPE1 => anticlockwise loop TYPE2 => clockwise loop
				       
				   sideways loop is second visit && second_fringe_edge_direction is DOWNWARDS && gauss_crossing_preceeds_virtual_crossing; or
				   sideways loop is first visit && second_fringe_edge_direction is UPWARDS && !gauss_crossing_preceeds_virtual_crossing; or
				   sideways loop is second visit && second_fringe_edge_direction is UPWARDS && gauss_crossing_preceeds_virtual_crossing; or
				   sideways loop is first visit && second_fringe_edge_direction is DOWNWARDS && !gauss_crossing_preceeds_virtual_crossing; 
				   				   
				       TYPE1 => clockwise loop TYPE2 => anticlockwise loop
				   				
				*/
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     virtual crossing involves sideways Gauss arc " << virtual_gauss_peers[i][0] << endl;
	         
				/* find the Gauss crossing containing the sideways arc */
				int gauss_crossing;
				bool found = false;
				
				for (int j=0; j < num_gauss_crossings && !found; j++)
				{
					for (int k=0; k< 4 && !found; k++)
					{
						if (PD_data[j][k] == virtual_gauss_peers[i][0])
						{
							gauss_crossing = j;
							found = true;
						}						
					}
				}
				
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     sideways_arc " << virtual_gauss_peers[i][0] << " meets gauss_crossing " << gauss_crossing << endl;
	
				/* identify whether the sideways arc was the first or second visit to gauss_crossing */
				bool sideways_second_visit = (virtual_gauss_peers[i][0] > virtual_gauss_peers[i][1]);
				
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     sideways_second_visit to gauss_crossing " << gauss_crossing << " " << sideways_second_visit << endl;
				
				int gauss_crossing_type = gauss_code_table[TYPE][gauss_crossing];

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     gauss_crossing " << gauss_crossing << " is (Gauss) " << (gauss_crossing_type == generic_code_data::TYPE1? "TYPE1": "TYPE2") << endl;
	
				/* identify the immersion labels on the non-sideways arc.  For the virtual crossing, the label on the sideways arc is always stored
				   in position zero.  For Gauss crossings, positions 0 and 1 store the labels corresponding to the second and first visit respectively.
				*/
				int gauss_crossing_non_sideways_immersion_label = (sideways_second_visit?immersion_crossing_peers[gauss_crossing][1]:immersion_crossing_peers[gauss_crossing][0]);				
				int virtual_crossing_non_sideways_immersion_label = immersion_crossing_peers[i][1];

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "gauss_to_peer_code:     gauss_crossing_non_sideways_immersion_label = " << gauss_crossing_non_sideways_immersion_label << endl;
	debug << "gauss_to_peer_code:     virtual_crossing_non_sideways_immersion_label = " << virtual_crossing_non_sideways_immersion_label << endl;
}

				/* identify the non-sideways component */
				int non_sideways_component_even_edge = gauss_crossing_non_sideways_immersion_label;
				if (non_sideways_component_even_edge %2)
					non_sideways_component_even_edge = virtual_crossing_non_sideways_immersion_label;
				int non_sideways_component = peer_code_table[COMPONENT][non_sideways_component_even_edge/2];

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     non_sideways_component = " << non_sideways_component << endl;

				bool gauss_crossing_preceeds_virtual_crossing = false;

				int last_component_edge = first_immersion_edge_on_component[non_sideways_component] + num_immersion_component_edges[non_sideways_component]-1;

				if (gauss_crossing_non_sideways_immersion_label == last_component_edge && 
				    virtual_crossing_non_sideways_immersion_label == first_immersion_edge_on_component[non_sideways_component])				    
				{
					gauss_crossing_preceeds_virtual_crossing = true;
				}
				else if (virtual_crossing_non_sideways_immersion_label == last_component_edge && 
				         gauss_crossing_non_sideways_immersion_label== first_immersion_edge_on_component[non_sideways_component])				    
				{
					gauss_crossing_preceeds_virtual_crossing = false;
				}
				else
				{
					gauss_crossing_preceeds_virtual_crossing = (gauss_crossing_non_sideways_immersion_label < virtual_crossing_non_sideways_immersion_label);
				}
				    
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     gauss_crossing_preceeds_virtual_crossing = " << gauss_crossing_preceeds_virtual_crossing << endl;

				bool sideways_loop_clockwise;

				if ((sideways_second_visit && second_fringe_edge_direction == gc_pc_xlabels::direction::UPWARDS && !gauss_crossing_preceeds_virtual_crossing) ||
			        (!sideways_second_visit && second_fringe_edge_direction == gc_pc_xlabels::direction::DOWNWARDS && gauss_crossing_preceeds_virtual_crossing) ||
			        (sideways_second_visit && second_fringe_edge_direction == gc_pc_xlabels::direction::DOWNWARDS && !gauss_crossing_preceeds_virtual_crossing) ||
			        (!sideways_second_visit && second_fringe_edge_direction == gc_pc_xlabels::direction::UPWARDS && gauss_crossing_preceeds_virtual_crossing)
			       )
			    {
					if (gauss_crossing_type == generic_code_data::TYPE1)
						sideways_loop_clockwise = false;
					else
						sideways_loop_clockwise = true;						
				}
				else
				{
					if (gauss_crossing_type == generic_code_data::TYPE1)
						sideways_loop_clockwise = true;
					else
						sideways_loop_clockwise = false;						
				}

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     sideways_loop_clockwise = " << sideways_loop_clockwise << endl;

				/*  Now we know the orientation of the sideways loop we can work out the immersion type of the virtual crossing:
				
				   sideways loop is clockwise && second_fringe_edge_direction is UPWARDS && !gauss_crossing_preceeds_virtual_crossing; or
				   sideways loop is anti-clockwise && second_fringe_edge_direction is DOWNWARDS && gauss_crossing_preceeds_virtual_crossing; or
				   sideways loop is clockwise && second_fringe_edge_direction is DOWNWARDS && !gauss_crossing_preceeds_virtual_crossing; or
				   sideways loop is anti-clockwise && second_fringe_edge_direction is UPWARDS && gauss_crossing_preceeds_virtual_crossing; 
				   
				       even_peer_in_place_zero => TYPE1     !even_peer_in_place_zero => TYPE2
				       
				   sideways loop is clockwise && second_fringe_edge_direction is DOWNWARDS && gauss_crossing_preceeds_virtual_crossing; or
				   sideways loop is anti-clockwise && second_fringe_edge_direction is UPWARDS && !gauss_crossing_preceeds_virtual_crossing; or
				   sideways loop is clockwise && second_fringe_edge_direction is UPWARDS && gauss_crossing_preceeds_virtual_crossing; or
				   sideways loop is anti-clockwise && second_fringe_edge_direction is DOWNWARDS && !gauss_crossing_preceeds_virtual_crossing;
				   
				       even_peer_in_place_zero => TYPE2     !even_peer_in_place_zero => TYPE1
				       				       
				*/	
				
				if ((sideways_loop_clockwise && second_fringe_edge_direction == gc_pc_xlabels::direction::UPWARDS && !gauss_crossing_preceeds_virtual_crossing) ||
			        (!sideways_loop_clockwise && second_fringe_edge_direction == gc_pc_xlabels::direction::DOWNWARDS && gauss_crossing_preceeds_virtual_crossing) ||
			        (sideways_loop_clockwise && second_fringe_edge_direction == gc_pc_xlabels::direction::DOWNWARDS && !gauss_crossing_preceeds_virtual_crossing) ||
			        (!sideways_loop_clockwise && second_fringe_edge_direction == gc_pc_xlabels::direction::UPWARDS && gauss_crossing_preceeds_virtual_crossing)
			       )
			    {
					if (even_peer_in_place_zero)
						peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE1;				
					else
						peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE2;				
				}
				else
			    {
					if (even_peer_in_place_zero)
						peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE2;				
					else
						peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE1;				
				}

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     immersion crossing " << even_peer/2 << " is " << (peer_code_table[TYPE][even_peer/2] == generic_code_data::type::TYPE1? "TYPE1": "TYPE2") << endl;
				
			}
			else if ((first_fringe_edge_direction == gc_pc_xlabels::direction::DOWNWARDS && second_fringe_edge_direction == gc_pc_xlabels::direction::DOWNWARDS) ||
			    (first_fringe_edge_direction == gc_pc_xlabels::direction::UPWARDS && second_fringe_edge_direction == gc_pc_xlabels::direction::UPWARDS)
			   )
			{				
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     both Gauss arcs have the same direction at virtual crossing " << i-num_gauss_crossings << endl;
				if (even_peer_in_place_zero)
				{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to left fringe Gauss arc of the virtual crossing " << i-num_gauss_crossings;
	
					peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE1;
				}
				else
				{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to right fringe Gauss arc of the virtual crossing " << i-num_gauss_crossings;
	
					peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE2;
				}
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << " crossing is (peer) " << (peer_code_table[TYPE][even_peer/2] == generic_code_data::type::TYPE1? "TYPE1": "TYPE2") << endl;
			}
			else
			{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     Gauss arcs have opposite direction at virtual crossing " << i-num_gauss_crossings << endl;
				if (even_peer_in_place_zero)
				{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to left fringe Gauss arc of the virtual crossing " << i-num_gauss_crossings;
	
					peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE2;
				}
				else
				{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to right fringe Gauss arc of the virtual crossing " << i-num_gauss_crossings;
	
					peer_code_table[TYPE][even_peer/2] = generic_code_data::type::TYPE1;
				}
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << " crossing is (peer) " << (peer_code_table[TYPE][even_peer/2] == generic_code_data::type::TYPE1? "TYPE1": "TYPE2") << endl;
			}
		}
	}
		
	/* write the crossing LABEL data to peer_code_table */
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: write crossing labels:" << endl;
	
	for (int i = 0; i< num_immersion_crossings; i++)
	{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code:   " << (i < num_gauss_crossings? "Gauss": "virtual") << " crossing " << (i < num_gauss_crossings? i: i-num_gauss_crossings);
	
		int even_peer;
		bool even_peer_in_place_zero;
		
		if (immersion_crossing_peers[i][0] %2 == 0)
		{
			even_peer = immersion_crossing_peers[i][0];
			even_peer_in_place_zero = true;
		}
		else
		{
			even_peer = immersion_crossing_peers[i][1];
			even_peer_in_place_zero = false;
		}
			
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << " even_peer = " << even_peer << endl;

		if (i < num_gauss_crossings)
		{
			if (even_peer_in_place_zero)  // even label corresponds to second visit to Gauss crossing
			{
				if (gauss_code_table[LABEL][i] == generic_code_data::label::POSITIVE)  //second visit to Gauss crossing goes over
					peer_code_table[LABEL][even_peer/2] = generic_code_data::label::POSITIVE;
				else if (gauss_code_table[LABEL][i] == generic_code_data::label::NEGATIVE)  //second visit to Gauss crossing goes under
					peer_code_table[LABEL][even_peer/2] = generic_code_data::label::NEGATIVE;
				else if (gauss_code_table[LABEL][i] == generic_code_data::label::FLAT)
					peer_code_table[LABEL][even_peer/2] = generic_code_data::label::FLAT;
			}
			else // even label corresponds to first visit to Gauss crossing
			{
				if (gauss_code_table[LABEL][i] == generic_code_data::label::POSITIVE)  //second visit to Gauss crossing goes over
					peer_code_table[LABEL][even_peer/2] = generic_code_data::label::NEGATIVE;
				else if (gauss_code_table[LABEL][i] == generic_code_data::label::NEGATIVE)  //second visit to Gauss crossing goes under
					peer_code_table[LABEL][even_peer/2] = generic_code_data::label::POSITIVE;
				else if (gauss_code_table[LABEL][i] == generic_code_data::label::FLAT)
					peer_code_table[LABEL][even_peer/2] = generic_code_data::label::FLAT;
			}
		}
		else
		{
			peer_code_table[LABEL][even_peer/2] = generic_code_data::label::VIRTUAL;
		}
	}
		
	/* write the originating and terminating vertices and edges */
	vector<int> term_crossing(2*num_immersion_crossings);
	vector<int> orig_crossing(2*num_immersion_crossings);
	
	for (int i=0; i< num_immersion_crossings; i++)
	{
		term_crossing[2*i] = i;
		orig_crossing[2*i+1] = i;
		term_crossing[peer_code_table[OPEER][i]] = i;
		
		int component = peer_code_table[COMPONENT][(peer_code_table[OPEER][i]-1)/2];
		int peer_successor = (peer_code_table[OPEER][i]+1 - first_immersion_edge_on_component[component])%
		                     num_immersion_component_edges[component] + first_immersion_edge_on_component[component];

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "gauss_to_peer_code: crossing " << i << " peer_successor = " << peer_successor << endl;

		orig_crossing[peer_successor] = i;
		
		peer_code_table[EVEN_TERMINATING][i] = 2*i;
		peer_code_table[ODD_ORIGINATING][i] = 2*i+1;
		peer_code_table[ODD_TERMINATING][i] = peer_code_table[OPEER][i];
		peer_code_table[EVEN_ORIGINATING][i] = peer_successor;
	}
	
	peer_code_data.type = generic_code_data::peer_code;
	peer_code_data.head = -1;
	peer_code_data.head_zig_zag_count = gauss_code_data.head_zig_zag_count;
	peer_code_data.immersion = gauss_code_data.immersion;
	peer_code_data.num_crossings = num_immersion_crossings;
	peer_code_data.num_components = num_components;
	peer_code_data.code_table = peer_code_table;
	peer_code_data.num_component_edges = num_immersion_component_edges;
	peer_code_data.first_edge_on_component = first_immersion_edge_on_component;
	peer_code_data.term_crossing = term_crossing;
	peer_code_data.orig_crossing = orig_crossing;

	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "gauss_to_peer_code: peer code data produced from gauss code:" << endl;
	print_code_data(peer_code_data,debug,"gauss_to_peer_code: ");	
}

    return true;	
}

void add_virtual_crossing(int location, vector<int>& fringe, vector<int>& num_virtual_crossings_on_gauss_arc, list<gc_pc_xlabels>& virtual_crossings) 
{	
    
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "add_virtual_crossing: adding virtual crossing " << virtual_crossings.size()+1 << ": location = " << location << ", gauss edges " << fringe[location] << ',' << fringe[location+1] << endl;
	
//if (braid_control::DEBUG >= braid_control::SUMMARY)
//	debug << location << ',' << (fringe[location]-1+num_gauss_arcs)%num_gauss_arcs << ',' << (fringe[location+1]-1+num_gauss_arcs)%num_gauss_arcs << endl;

    gc_pc_xlabels xlabels;
    xlabels.gauss.first = fringe[location];
    xlabels.gauss.second = fringe[location+1];
    
    num_virtual_crossings_on_gauss_arc[xlabels.gauss.first]++;
    num_virtual_crossings_on_gauss_arc[xlabels.gauss.second]++;
    
    virtual_crossings.push_back(xlabels);
    
    swap (fringe[location],fringe[location+1]);   
}

/* arc_1 and arc_2 are both incident at the supplied Gauss crossing, we set the direction of the terminating arc to UPWARDS and
   that of the originating arc DOWNWARDS
*/
void assign_gauss_arc_direction(int arc_1, int arc_2, matrix<int>& gauss_code_table, int crossing, vector<int>& gauss_arc_direction)
{
		if (gauss_code_table[ODD_TERMINATING][crossing] == arc_1 || gauss_code_table[EVEN_TERMINATING][crossing] == arc_1 )
		{		
			gauss_arc_direction[arc_1] = gc_pc_xlabels::direction::UPWARDS;
			gauss_arc_direction[arc_2] = gc_pc_xlabels::direction::DOWNWARDS;
		}
		else
		{
			gauss_arc_direction[arc_1] = gc_pc_xlabels::direction::DOWNWARDS;
			gauss_arc_direction[arc_2] = gc_pc_xlabels::direction::UPWARDS;
		}
}

string direction_to_string(int edge, vector<int>& gauss_arc_direction)
{
	switch(gauss_arc_direction[edge])
	{
		case (gc_pc_xlabels::direction::UPWARDS): return "U ";
		case (gc_pc_xlabels::direction::DOWNWARDS): return "D ";
		case (gc_pc_xlabels::direction::SIDEWAYS): return "S ";
		default: return "- ";
	}
}

void remove_fringe_edge(int edge, vector<int>& current_fringe)
{
	int current_fringe_size = current_fringe.size();
	vector<int> new_fringe(current_fringe_size - 2);
	
	int index = 0;
	for (int i=0; i<current_fringe_size; i++)
	{
		if (current_fringe[i] != edge)
			new_fringe[index++] = current_fringe[i];
	}
	
	current_fringe = new_fringe;
}
