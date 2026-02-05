/**************************************************************************
void generic_code(string input_string, string title)
bool immersion_to_dowker (matrix<int>& code_table, vector<int>& code)
void affine_index_polynomial(string input_string, generic_code_data& code_data)
bool satellite_code_data(generic_code_data& knot_data, generic_code_data& satellite_data, int strands)
void assign_satellite_code(generic_code_data& satellite_data, int s_edge, int s_prime_edge, int s_component, int s_prime_component, int type, int crossing)
vector<int> gauss_parity(generic_code_data& code_data)
string minimal_peer_code_representation (generic_code_data peer_code_data)
void manturov_alexander(generic_code_data& code_data)
void mock_alexander(generic_code_data& code_data)
polynomial<int> nabla_k(generic_code_data& code_data, int starred_edge)	
int find_cycle(matrix<int>& cycle, int num_cycles, int num_left_cycles, int edge_1, int edge_2, bool left_cycle)
void doodle_Q_polynomial(generic_code_data code_data)
generic_code_data isolate_component(generic_code_data& code_data, int component,vector<int>& crossing_map)
bool smooth_diagram(generic_code_data& code_data, int smoothed_crossing, generic_code_data& smoothed_code_data)
int linking_number (generic_code_data& code_data, int crossing_type, int component_1, int component_2)
list<vector<int> > calculate_colourings( generic_code_data& code_data, matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, matrix<int>& Tu, matrix<int>& Td)
vector<Cpolynomial> cocycle_invariant(generic_code_data& code_data, generic_switch_data& switch_data, list<vector<int> >& colourings)
void peer_code_colouring_invariant(matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, matrix<int>& Tu, matrix<int>& Td,
			braid_control::ST_pair_type pair_type, string input_string, string title, generic_switch_data& switch_data, int period)
**************************************************************************/
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstring>
#include <valarray>
#include <list>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <map>

using namespace std;

namespace util
{
	string itos(long n);
}


extern ofstream     debug;
extern ofstream     output;
extern ifstream     input;

#include <util.h>
#include <scalar.h> // includes bigint.h and rational.h
#include <quaternion-scalar.h>
#include <polynomial.h>
#include <matrix.h>
#include <braid.h>
#include <braid-util.h>
#include <generic-code.h>
#include <gauss-orientation.h>

//bool CYCLE_KNOT_TYPE_NABLA_K = true;

/********************* Function prototypes ***********************/
void bracket_polynomial(generic_code_data& code_data, int variant);
bool immersion_to_dowker (matrix<int>& code_table, vector<int>& code);
bool gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data, bool optimal=true, vector<int>* gauss_crossing_perm=0, bool evaluate_gauss_crossing_perm=false);
//bool classical_gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data);
void standard_rep (string& word);
string vogel (generic_code_data code_data, int* turning_number_ptr=0);
string affine_index_polynomial(string input_string, generic_code_data& code_data);
bool satellite_code_data(generic_code_data& knot_data, generic_code_data& satellite_data, int strands);
void assign_satellite_code(generic_code_data& satellite_data, int s_edge, int s_prime_edge, int s_component, int s_prime_component, int type, int crossing);
int remove_peer_code_component(generic_code_data& code_data, int component, vector<int>& component_flags);
string parse_long_knot_input_string (string input_string);
void add_virtual_crossing(int location, vector<int>& fringe, vector<int>& num_virtual_crossings_on_gauss_arc, list<gc_pc_xlabels>& virtual_crossings);
void assign_gauss_arc_direction(int arc_1, int arc_2, matrix<int>& gauss_code_table, int crossing, vector<int>& gauss_arc_direction);
string direction_to_string(int edge, vector<int>& gauss_arc_direction);
void remove_fringe_edge(int edge, vector<int>& current_fringe);

void manturov_alexander(generic_code_data& code_data);
void mock_alexander(generic_code_data& code_data);
int find_cycle(matrix<int>& cycle, int num_cycles, int num_left_cycles, int edge_1, int edge_2, bool left_cycle);
list<vector<int> > hamiltonian_circuit(generic_code_data& code_data, bool list_all_circuits, bool count_circuits_only, bool edge_circuit,int include_edge);

void doodle_Q_polynomial(generic_code_data code_data);
bool smooth_diagram(generic_code_data& code_data, int smoothed_crossing, generic_code_data& smoothed_code_data);
generic_code_data isolate_component(generic_code_data& code_data, int component, vector<int>& crossing_map);
int linking_number (generic_code_data& code_data, int crossing_type, int component_1, int component_2);
void print_k_chain(ostream& os, const vector<int> k_chain, generic_switch_data& switch_data, int k);
void birack_homology_generators(generic_switch_data& switch_data, int k, bool cohomology=false, matrix<scalar>* _B=0);
void determine_cohomology_generators(generic_switch_data& switch_data);
void test_cohomology_generators(generic_switch_data& switch_data);
vector<int> smallest_parent_birack(matrix<int>& Su, matrix<int>& Sd, vector<int> labels);
generic_code_data add_Reidemeister_1_loop (generic_code_data& code_data, bool positive_crossing, bool positive_turn = true);
void display_fixed_point_switch(matrix<int>& M, ostream& os, bool number_from_zero);

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

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "generic_code: provided with input string: " << input_string << endl;

	/* first isolate any qualifiers from the input string */
	string qualifier;
	unsigned int satellite_strands = braid_control::SATELLITE;
	
	string::size_type pos = input_string.find('{');
	if (pos != string::npos)
	{
		qualifier = input_string.substr(pos,string::npos);
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "generic_code: isolated qualifier = " << qualifier << endl;
		input_string = input_string.substr(0,pos);
	}

	/* if we have a satellite qualifier, override the global valuse set by braid_control::SATELLITE with the value determined by the qualifier */
	if (qualifier.find("satellite") != string::npos)
	{
		string::size_type pos = qualifier.find("satellite");
		string::size_type next = qualifier.find(',',pos); //next == string::npos if there's no next qualifier
		string satellite_string = qualifier.substr(pos,next);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "generic_code: satellite_string: " << satellite_string << endl;

		pos = satellite_string.find('=');


		if (pos != string::npos)
		{
			char* c_satellite_string = c_string(satellite_string);			
			char* cptr = c_satellite_string;	
			while (*cptr != '=')
				cptr++;
				
			get_number(satellite_strands,++cptr);			
			delete c_satellite_string;
		}
		else
			satellite_strands = 2;


if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "generic_code: satellite qualifier provided, number of strands = " << satellite_strands << endl;

	}

	/* if we have a cycle qualifier, override the global value set by braid_control::INFINITE_CYCLE with the value determined by the qualifier */
	braid_control::INFINITE_CYCLE = braid_control::cycle::UNSPECIFIED;
	if (qualifier.find("cycle") != string::npos)
	{
		string::size_type pos = qualifier.find("cycle");
		string::size_type next = qualifier.find(',',pos); //next == string::npos if there's no next qualifier
		string cycle_string = qualifier.substr(pos,next);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "generic_code: cycle_string: " << cycle_string << endl;

		pos = cycle_string.find('=');


		if (pos != string::npos)
		{
			char* c_cycle_string = c_string(cycle_string);			
			char* cptr = c_cycle_string;	
			while (*cptr != '=')
				cptr++;
				
			get_number(braid_control::INFINITE_CYCLE,++cptr);			
			delete c_cycle_string;
		}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "generic_code: cycle qualifier provided, infinite turning cycle = " << braid_control::INFINITE_CYCLE << endl;

	}

	bool plane_reflect_code_qualifier = false;
	if (qualifier.find("plane-reflect") != string::npos)
	{
		plane_reflect_code_qualifier = true;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "generic_code: plane-reflect qualifier provided" << endl;
	}
		
	bool reverse_code_orientation_qualifier = false;
	if (qualifier.find("reverse") != string::npos)
	{
		reverse_code_orientation_qualifier = true;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "generic_code: reverse qualifier provided" << endl;
	}
		
	if (title.length())
   	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\n\n" << title << endl;
		if (!braid_control::RAW_OUTPUT)
			output << "\n\n" << title;
   	}
   	else if (!braid_control::SILENT_OPERATION)
   	{
		cout << "\n\n" << input_string << endl;
	}
		

	if (!braid_control::RAW_OUTPUT)
	{
		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
	   	output << input_string;		
	}
						
	if (input_string.find("L:") != string::npos)
	{	
		input_string = parse_long_knot_input_string(input_string);
		
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
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
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "generic_code: Error! long knot indicator provided for the peer code of a link" << endl;

			return;
		}
	}

	generic_code_data code_data;
	read_code_data (code_data, input_string);

	if (braid_control::PLANE_REFLECT_INPUT || plane_reflect_code_qualifier)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "generic_code: replacing code data ";
	write_code_data(debug,code_data);
	debug << endl;
}
		/* plane reflect the braid, which means inverting the positive and negative labels in code_data */
		for (int i=0; i< code_data.num_crossings; i++)
		{
			if (code_data.code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::POSITIVE)
				code_data.code_table[generic_code_data::table::LABEL][i] = generic_code_data::label::NEGATIVE;
			else if (code_data.code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::NEGATIVE)
				code_data.code_table[generic_code_data::table::LABEL][i] = generic_code_data::label::POSITIVE;
		}
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "generic_code: plane-reflect code data ";
	write_code_data(debug,code_data);
	debug << endl;
}
	}
	
	if (braid_control::REVERSE_INPUT_ORIENTATION || reverse_code_orientation_qualifier)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "generic_code: replacing code data ";
	write_code_data(debug,code_data);
	debug << endl;
}
		vector<int> shift(code_data.num_components);
		for (int i=0; i< code_data.num_components; i++)
			shift[i] = -1*code_data.num_component_edges[i];
			
		renumber_peer_code(code_data,shift);
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "generic_code: reverse orientation code data ";
	write_code_data(debug,code_data);
	debug << endl;
}

	   	if (!braid_control::SILENT_OPERATION)
	   	{
			cout << "reversed input ";
			write_code_data(cout,code_data);
			cout << endl;
		}

		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? ";" : "");
			output << "reversed input ";
			write_code_data(output,code_data);
			output << endl;
		}
	}

	
	if (satellite_strands)
	{

/*		if (braid_control::SATELLITE)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "generic_code: ignoring satellite code qualifier due to braid_control::SATELLITE being specified" << endl;
		}
		else
*/		
		{	
			generic_code_data satellite_data;
			
			if (satellite_code_data(code_data, satellite_data, satellite_strands))
			{
				code_data = satellite_data;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "generic_code: replacing code data with satellite code data" << endl;
	print_code_data(debug,code_data,"generic_code: ");	
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
				if (code_table[generic_code_data::table::LABEL][i] != generic_code_data::FLAT)
				{
					flat_crossings_only = false;
					break;
				}
			}
		    
		    if (flat_crossings_only)
		    {
				braid_control::FLAT_CROSSINGS = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "generic_code: input code contains only flat crossings, setting programme option FLAT_CROSSINGS" << endl;
				
			}
		    
		    int turning_number;
			string braid_word = vogel(code_data,&turning_number);
					
			if (braid_control::VOGEL_TURNING_NUMBER)
			{
				if (!braid_control::SILENT_OPERATION)
					cout << "\nTurning number = " << turning_number << endl;

				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Turning number = ";
					if (braid_control::OUTPUT_AS_INPUT)
						output << '\n';
				}
				output << turning_number << endl;
			}


			{
				if (!braid_control::SILENT_OPERATION)
					cout << "Braid word = " << braid_word << endl;
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Braid word = ";
					if (braid_control::OUTPUT_AS_INPUT)
						output << '\n';
				}
				output << braid_word << endl;
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
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
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "generic_code:            = " << braid_word << endl;
				}
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
		if(code_data.type == generic_code_data::peer_code || code_data.type == generic_code_data::immersion_code)
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
				output << "DT:";
				
				for (int i=0;i<num_terms;i++)
				{
		    		output << dowker_code[i];
		    		
		    		if (i< num_terms-1)
						output << " ";
		    	}
		    		
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
			cout << "\nError! generic_code presented with a a code that is not a peer code or an immersion code for a call to immersion_to_dowker.\n";
			exit(0);
		}
	}
	else if (braid_control::PEER_CODE)
	{
		if(code_data.type == generic_code_data::gauss_code)
		{
			generic_code_data peer_code_data;
//			if ((braid_control::CLASSICAL_INPUT && classical_gauss_to_peer_code(code_data, peer_code_data)) || gauss_to_peer_code(code_data, peer_code_data) )		    
			if (gauss_to_peer_code(code_data, peer_code_data))
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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
			
if (debug_control::DEBUG >= debug_control::DETAIL)
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
		else if (code_data.type == generic_code_data::peer_code)
		{
			if (braid_control::REMOVE_PEER_CODE_COMPONENT)
			{
				if (0 <= braid_control::REMOVE_COMPONENT && braid_control::REMOVE_COMPONENT < code_data.num_components)
				{
					vector<int> dummy_flags;  // not tracking any components being removed
					remove_peer_code_component(code_data, braid_control::REMOVE_COMPONENT,dummy_flags);
				}
				else
				{
					cout << "\nError! Invalid component number provided in input with REMOVE_PEER_CODE_COMPONENT option.\n";
					exit(0);
				}
			}

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
			cout << "\nError! Unsupported combination code data when PEER_CODE is true.\n"; // shouldn't ever get here
			exit(0);
		}
	}
	else if (braid_control::GAUSS_CODE)
	{
		if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && code_data.head != -1)
		{
			if (!valid_knotoid_input(code_data))
			{
				cout << "\nError! Gauss code task presented with an invalid knotoid.\n";
				exit(0);
			}
			
		}
		else if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID)
		{
			if (!valid_multi_linkoid_input(code_data))
			{
				cout << "\nError! Gauss code task presented with an invalid multi-linkoid.\n";
				exit(0);
			}
			
		}
		
		
		ostringstream oss;
		if (braid_control::PD_FORMAT)
		{
			if (code_data.type != generic_code_data::code_type::gauss_code)
			{
				ostringstream gauss_oss;
				write_gauss_code(gauss_oss,code_data);
				generic_code_data gauss_code_data;
				read_gauss_code(gauss_code_data,gauss_oss.str());
				write_planar_diagram(oss,gauss_code_data);
			}
			else
			{
				write_planar_diagram(oss,code_data); 
			}
		}
		else if (braid_control::LPGD)
		{
			gauss_orientation_data gauss_data(code_data);
			gauss_orientation_data lp_gauss_data = left_preferred(gauss_data); //unoriented = false;
			write_gauss_data(lp_gauss_data,oss); 
		}
		else if (braid_control::ULPGD)
		{
			gauss_orientation_data gauss_data(code_data);
			gauss_orientation_data lp_gauss_data = left_preferred(gauss_data,true, code_data.immersion_character); //unoriented = true;
			write_gauss_data(lp_gauss_data,oss); 
		}
		else if (braid_control::OPGC)
		{
			oss << over_preferred_gauss_code(code_data, false);
		}
		else if (braid_control::UOPGC)
		{
			oss << over_preferred_gauss_code(code_data, true);
		}
		else
			write_gauss_code(oss, code_data,braid_control::OU_FORMAT);
		
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
	debug << "\ngeneric_code: Gauss code = " << oss.str() << endl;
		
	}
	else if (braid_control::HAMILTONIAN)
	{
		if (code_data.immersion_character != generic_code_data::character::CLOSED)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\n\nHamiltonian circuits only defined for closed immersions, skipping";
	
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "\n\nHamiltonian circuits only defined for closed immersions, skipping";
				if (braid_control::OUTPUT_AS_INPUT)
					output << '\n';
			}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "generic_code: Hamiltonian circuits only defined for closed immersions, doing nothing" << endl;
		}
		else
		{
			list<vector<int> > circuit_list = hamiltonian_circuit(code_data, braid_control::HC_LIST_ALL, braid_control::HC_COUNT, braid_control::HC_EDGES, braid_control::HC_INCLUDE_EDGE);
			
			if (circuit_list.size() != 0)
			{
				list<vector<int> >::iterator cptr = circuit_list.begin();
				while (cptr != circuit_list.end())
				{
					vector<int>& circuit = *cptr;
					
						if (!braid_control::HC_COUNT)
						{
							if (!braid_control::SILENT_OPERATION)
							{
								cout << "Hamiltonian circuit ";
								for (int i=0; i< code_data.num_crossings; i++)
									cout << circuit[i] << ' ';
								cout << endl;
							}
							
							if (!braid_control::RAW_OUTPUT)
							{
								output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
								output << "Hamiltonian circuit ";
							}
							if (braid_control::OUTPUT_AS_INPUT)
								output << "\n";
			
							for (int i=0; i< code_data.num_crossings; i++)
								output << circuit[i] << ' ';
							output << endl;
			
							if (!braid_control::SILENT_OPERATION)
								cout << "\n";
								
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
		debug << "generic_code: hamiltonian_circuit: " << endl;
		for (int i=0; i< code_data.num_crossings; i++)
			debug << circuit[i] << ' ';
		debug << endl;
}
						}
						
					cptr++;
				}
			}
			else
			{
				if (!braid_control::SILENT_OPERATION)
					cout << "No Hamiltonian circuit found" << endl;
				
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "No Hamiltonian circuit found" << endl;
				}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
		debug << "generic_code: hamiltonian_circuit: no Hamiltonian circuit found" << endl;
			}
			
			if (braid_control::HC_COUNT)
			{
				if (!braid_control::SILENT_OPERATION)
					cout << "Found " << circuit_list.size() << " Hamiltonian circuits" << endl;
				
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Found " << circuit_list.size() << " Hamiltonian circuits" << endl;
				}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
		debug << "generic_code: hamiltonian_circuit: found " << circuit_list.size() << " Hamiltonian circuits" << endl;
			}			
		}
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
		if (code_data.type == generic_code_data::code_type::gauss_code)
		{
			generic_code_data peer_code_data;
			if (gauss_to_peer_code(code_data, peer_code_data, true))
			{
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "generic_code: nusing immersion ";
	write_peer_code(debug,peer_code_data);
	debug << " to calaulate the knotoid bracket polynomial" << endl;
}
				bracket_polynomial(peer_code_data,TURAEV_VARIANT);
			}
			else
			{
				if (!braid_control::SILENT_OPERATION)
					cout << "\n\nunsuccessful conversion from Gauss code to peer code, skipping";
		
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "\n\nunsuccessful conversion from Gauss code to peer code, skipping";
					if (braid_control::OUTPUT_AS_INPUT)
						output << '\n';
				}
		
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "generic_code: unsuccessful conversion from Gauss code to peer code, doing nothing" << endl;
		
				return;
			}
		}
		else
		{
			bracket_polynomial(code_data,TURAEV_VARIANT);
		}
	}
	else if (braid_control::ARROW_POLYNOMIAL)
	{
		bracket_polynomial(code_data,ARROW_VARIANT);
	}
	else if (braid_control::PARITY_BRACKET)
	{
		if (code_data.type == generic_code_data::code_type::gauss_code)
		{
			generic_code_data peer_code_data;
			if (gauss_to_peer_code(code_data, peer_code_data, true))
			{
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "generic_code: nusing immersion ";
	write_peer_code(debug,peer_code_data);
	debug << " to calaulate the parity bracket polynomial" << endl;
}

				bracket_polynomial(peer_code_data,PARITY_VARIANT);
			}
			else
			{
				if (!braid_control::SILENT_OPERATION)
					cout << "\n\nunsuccessful conversion from Gauss code to peer code, skipping";
		
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "\n\nunsuccessful conversion from Gauss code to peer code, skipping";
					if (braid_control::OUTPUT_AS_INPUT)
						output << '\n';
				}
		
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "unsuccessful conversion from Gauss code to peer code, doing nothing" << endl;
		
				return;
			}
		}
		else
		{
			bracket_polynomial(code_data,PARITY_VARIANT);
		}
	}
	else if (braid_control::PARITY_ARROW)
	{
		if (code_data.type == generic_code_data::code_type::gauss_code)
		{
			generic_code_data peer_code_data;
			if (gauss_to_peer_code(code_data, peer_code_data, true))
			{
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "generic_code: nusing immersion ";
	write_peer_code(debug,peer_code_data);
	debug << " to calaulate the parity arrow polynomial" << endl;
}				
				bracket_polynomial(peer_code_data,PARITY_ARROW_VARIANT);
			}
			else
			{
				if (!braid_control::SILENT_OPERATION)
					cout << "\n\nunsuccessful conversion from Gauss code to peer code, skipping" << endl;
		
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "\n\nunsuccessful conversion from Gauss code to peer code, skipping";
					if (braid_control::OUTPUT_AS_INPUT)
						output << endl;
				}
		
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "generic_code: unsuccessful conversion from Gauss code to peer code, doing nothing" << endl;
		
				return;
			}
		}
		else
		{
			bracket_polynomial(code_data,PARITY_ARROW_VARIANT);
		}
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
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "generic_code: Affine index polynomial = " << poly << endl;
	}
	else if (braid_control::MANTUROV_ALEXANDER)
	{
		if(code_data.type == generic_code_data::gauss_code)
		{
			generic_code_data peer_code_data;
		    if(gauss_to_peer_code(code_data, peer_code_data))
		    {				
				manturov_alexander(peer_code_data);
			}
		}
		else
		{
			manturov_alexander(code_data);
		}
	}
	else if (braid_control::MOCK_ALEXANDER)
	{
		if(code_data.type == generic_code_data::gauss_code)
		{
			generic_code_data peer_code_data;
		    if(gauss_to_peer_code(code_data, peer_code_data))
		    {				
				mock_alexander(peer_code_data);
			}
		}
		else
		{
			mock_alexander(code_data);
		}
	}
	else if (braid_control::PRIME_TEST)
	{
		if(code_data.type != generic_code_data::peer_code)
		{			
			generic_code_data peer_code_data;
			if (gauss_to_peer_code(code_data, peer_code_data, true))
			{
				code_data = peer_code_data;
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "generic_code: using immersion ";
	write_peer_code(debug,code_data);
	debug << " to to determine whether Gauss code is prime" << endl;
}
			}
			else
			{
				if (!braid_control::SILENT_OPERATION)
					cout << "\n\nunsuccessful conversion from Gauss code to peer code, skipping" << endl;
		
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "\n\nunsuccessful conversion from Gauss code to peer code, skipping";
					if (braid_control::OUTPUT_AS_INPUT)
						output << endl;
				}
		
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "generic_code: unsuccessful conversion from Gauss code to peer code, doing nothing" << endl;
		
				return;
			}
		}

		bool flat_crossings = true;
		for (int i=0; i< code_data.num_crossings; i++)
		{
			if (code_data.code_table[generic_code_data::table::LABEL][i] != generic_code_data::label::FLAT &&
			    code_data.code_table[generic_code_data::table::LABEL][i] != generic_code_data::label::VIRTUAL)
				flat_crossings = false;
		}
		
		if (braid_control::FLAT_CROSSINGS)
			flat_crossings = true;
					
		/* check input and set shortcut crossings for multi-linkoids and knotoids */
		if ( (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID && !valid_multi_linkoid_input(code_data))  ||    
		     ((code_data.immersion_character == generic_code_data::character::KNOTOID || code_data.immersion_character == generic_code_data::character::PURE_KNOTOID) 
		       && !valid_knotoid_input(code_data)))
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\n\nInvalid knotoid or multi-linkoid input to prime test, skipping" << endl;
	
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "\n\ninvalid knotoid or multi-linkoid input to prime test, skipping";
				if (braid_control::OUTPUT_AS_INPUT)
					output << endl;
			}
		
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "generic_code: invalid knotoid or multi-linkoid input to prime test, doing nothing" << endl;
		
			return;
		}

		if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID ||    
		    code_data.immersion_character == generic_code_data::character::KNOTOID ||    
		    code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
		{
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "generic_code: shortcut_crossings: ";
	for (size_t i=0; i< code_data.shortcut_crossing.size(); i++)
		debug << code_data.shortcut_crossing[i] << ' ';
	debug << endl;
}
		}
			
		int prime_diagram = three_connected(code_data,flat_crossings);
		
		if (prime_diagram == -1)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\nUndefined" << endl;
			if (!braid_control::RAW_OUTPUT)
				output << "\nUndefined";
		}
		else
		{					
			if (!braid_control::SILENT_OPERATION)
				cout << "\nInput is " << (prime_diagram? "prime": "not prime")  << endl;
			if (!braid_control::RAW_OUTPUT)
			{
				output << "\nInput is ";
			}
			output <<  (prime_diagram? "prime": "not prime")  << endl;	
		}
	}
	else if (braid_control::DOODLE_Q_POLYNOMIAL)
	{
		doodle_Q_polynomial(code_data);
	}
/*
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
*/	
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

/* the function immersion_to_dowker calculates the dowker code for the labelled peer code or immersion code in code_table, storing the result at code.  
   
   The function returns a boolean indicating whether a code has been written to 'code'.  If the code includes a virtual label the function  
   returns false and leaves code unchanged.  Otherwise the code is set to a vector containing the dowker code and the function returns true.

   Given that we have a code table most of the work for the dowker code is done except that the dowker code numbers semi-arcs from 1 rather 
   than zero.  (See comments to the function braid_to_dowker for a description of the dowker code).  Therefore the generic_code_data::table::OPEER row of the 
   table gives the Dowker code once the odd-numbered peers have been incremented by one.   
   
   Note also that if we are given generic code data for an alternating knot then all the labels will be '+'.  Those crossings for which the label 
   is set to '-' are precisely those that require the corresponding dowker code element to be negated.
*/
bool immersion_to_dowker (matrix<int>& code_table, vector<int>& code)
{
	int num_crossings=code_table.numcols();
	bool classical = true;
	
	vector<int> dowker_code(num_crossings);
	
	for (int i=0; i<num_crossings; i++)
		dowker_code[i] = code_table[generic_code_data::table::OPEER][i]+1;
	
	for (int i=0; i<num_crossings; i++)
	{
		if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
		{
			dowker_code[i] *= -1;
		}
		else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::VIRTUAL)
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

if (debug_control::DEBUG >= debug_control::SUMMARY)
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
		if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::FLAT)
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

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "affine_index_polynomial: input includes flat crossings, doing nothing" << endl;

			return "Error!";
		}
	}
	
	bool pure_knotoid_code_data = false;
	
	if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && code_data.head != -1)
	{
		if (valid_knotoid_input(code_data))
		{
			pure_knotoid_code_data = true;
			
if (debug_control::DEBUG >= debug_control::BASIC)
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
		
		if (code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::VIRTUAL || (pure_knotoid_code_data && code_data.shortcut_crossing[next_crossing]))
		{
			/* do nothing, just move around the diagram */
			edge++;
			continue;
		}				
		else if (code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::POSITIVE || code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::NEGATIVE)
		{
			bool a_arc = false;
			if ( (edge == code_table[generic_code_data::table::ODD_TERMINATING][next_crossing] && code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1)
			   ||(edge == code_table[generic_code_data::table::EVEN_TERMINATING][next_crossing] && code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE2)
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
	
if (debug_control::DEBUG >= debug_control::BASIC)
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
		if (code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && code_data.shortcut_crossing[i]))
			W_plus[i] = a_colour[i] - b_colour[i] -1;
	}

if (debug_control::DEBUG >= debug_control::BASIC)
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
			sign[i] = generic_braid_data::crossing_type::VIRTUAL;
		}
		else if (   (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1 && code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
		    || (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE2 && code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE)
		   )
		{
			/* positive crossing */
			sign[i]  = generic_braid_data::crossing_type::POSITIVE;
	    }
		else if (   (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1 && code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE)
		          || (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE2 && code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
		        )
		{
			/* negative crossing */
			sign[i] = generic_braid_data::crossing_type::NEGATIVE;
	    }
	    else
		{
			/* virtual crossing */
			sign[i] = generic_braid_data::crossing_type::VIRTUAL;
	    }			    		    
	}
	
if (debug_control::DEBUG >= debug_control::BASIC)
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
		if (sign[i] == generic_braid_data::crossing_type::POSITIVE)
			writhe ++;
		else if (sign[i] == generic_braid_data::crossing_type::NEGATIVE)
			writhe --;
	}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "affine_index_polynomial: writhe = " << writhe << endl;
	
	/* write the polynomial to a ostringstream, then read it back again to put it into standard format */
	ostringstream oss;
	
	for (int i=0; i< num_crossings; i++)
	{
//		if (code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL)
		if (sign[i] != generic_braid_data::crossing_type::VIRTUAL)
		{
			if (sign[i] == generic_braid_data::crossing_type::POSITIVE)
				oss << '+';
			else // (sign[i] == generic_braid_data::crossing_type::NEGATIVE)
				oss << '-';
				
			oss << "t^";
			
			if (sign[i] == generic_braid_data::crossing_type::POSITIVE)
				oss << W_plus[i];
			else 
				oss << -W_plus[i];
		}
	}
	
	if (writhe < 0)
		oss << '+' << abs(writhe);
	else if (writhe > 0)
		oss << '-' << writhe;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "affine_index_polynomial: initial_polynomial = " << oss.str() << endl;
	
	polynomial<int> initial_polynomial(oss.str());
	
	oss.str("");
	oss << initial_polynomial;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "affine_index_polynomial: returning polynomial " << oss.str() << endl;

	return oss.str();
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
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "satellite_code_data: knot data does not describe a knot" << endl;
		satellite_data = knot_data;
		return false;
	}
	else
	{
		int num_kcrossings = knot_data.num_crossings;		
		matrix<int>& kcode_table = knot_data.code_table;		

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "satellite_code_data: presented with knot data" << endl;
	print_code_data(debug,knot_data,"satellite_code_data:   ");
	debug << "satellite_code_data: number of strands required = " << strands << endl;
}

		int writhe = 0;
		
		for (int i=0; i< num_kcrossings; i++)
		{
			if (kcode_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE || kcode_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
				writhe += kcode_table[generic_code_data::table::LABEL][i];
		}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "satellite_code_data: writhe = " << writhe << endl;
		
		/* Evaluate the number of crossings in the satellite link, L: each crossing in the knot contributes S^2 crossings
		   to L and in the writhe twists there are S(S-1)/2 crossings in every twist, and there are two twists per writhe.
		*/
		int num_scrossings = num_kcrossings*strands*strands + abs(writhe)*strands*(strands-1);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "satellite_code_data: evaluate satellite non-writhe crossings" << endl;

        for (int i=0; i< num_kcrossings; i++)
        {
	
			int odd_peer = kcode_table[generic_code_data::table::OPEER][i];
			bool type_1_crossing = (kcode_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1? true: false);
			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "satellite_code_data:   crossing " << i;
	debug << ", odd peer = " << odd_peer << ", type " << (type_1_crossing? "1":"2") << endl;
}			
			for (int s=0; s < strands; s++)
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "satellite_code_data:     strand " << s<< endl;
				
				for (int k=0; k < strands; k++)
				{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "satellite_code_data:       k= " << k << ": ";
					int s_edge = satellite_data.first_edge_on_component[s] + s%2 + 2*i*strands + k;
					int s_prime = (type_1_crossing? k: strands-1-k);
					int s_prime_edge = satellite_data.first_edge_on_component[s_prime] + s_prime%2 + odd_peer*strands + (type_1_crossing? strands-1-s: s);

					if ( s_prime_edge == satellite_data.first_edge_on_component[s_prime] + satellite_data.num_component_edges[s_prime])
							s_prime_edge = satellite_data.first_edge_on_component[s_prime];					

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "s_edge = " << s_edge << ", s' = " << s_prime << ", s_prime_edge = " << s_prime_edge << endl;


					int required_type;
					if ((type_1_crossing && s_prime_edge%2) || (!type_1_crossing && s_edge%2))
						required_type = generic_code_data::TYPE1;
					else
						required_type = generic_code_data::TYPE2;
						
					int required_crossing;
					if ((kcode_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE && kcode_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE2) ||
					     (kcode_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE && kcode_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1))
					    required_crossing = generic_code_data::POSITIVE;
					else if ((kcode_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE && kcode_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1) ||
					     (kcode_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE && kcode_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE2))
					    required_crossing = generic_code_data::NEGATIVE;
					else
						required_crossing = kcode_table[generic_code_data::table::LABEL][i];

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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "satellite_code_data: writhe twists base_edge: ";
	for (int i=0; i< strands; i++)
		debug << base_edge[i] << " ";
	debug << endl;
}

		
		for (int w=0; w< abs(writhe); w++)
		{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "satellite_code_data: writhe twist " << w+1 << endl;
	
			for (int s=0; s< strands-1; s++)
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "satellite_code_data:   strand " << s << endl;
				for (int k=0; k<=strands-2-s;k++)
				{
					int s_prime=strands-1-k;
					int s_edge = base_edge[s]+k;
					int s_prime_edge = base_edge[s_prime]+strands-2-s;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "satellite_code_data:   second meeting s_edge = " << s_edge << ", s_prime_edge = " << s_prime_edge << endl;

					assign_satellite_code(satellite_data, s_edge, s_prime_edge, s, s_prime, 
										  (s_prime_edge%2? generic_code_data::TYPE2 : generic_code_data::TYPE1),
										  (writhe < 0? generic_code_data::POSITIVE: generic_code_data::NEGATIVE));
				}				
			}
			
			/* adjust base_edge for the next write twists */
			for (int i=0; i< strands; i++)
				base_edge[i] += 2*(strands-1);
			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
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
			satellite_data.term_crossing[scode_table[generic_code_data::table::EVEN_TERMINATING][i]] = i;
			satellite_data.term_crossing[scode_table[generic_code_data::table::ODD_TERMINATING][i]] = i;
			satellite_data.orig_crossing[scode_table[generic_code_data::table::EVEN_ORIGINATING][i]] = i;
			satellite_data.orig_crossing[scode_table[generic_code_data::table::ODD_ORIGINATING][i]] = i;
		}

//		print_code_data(cout, satellite_data, "    ");
		
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
					
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "assign_satellite_code_data:       s_even = " << s_even << ", s_odd = " << s_odd << endl;
	
	int s_crossing = s_even/2;

	scode_table[generic_code_data::table::OPEER][s_crossing] = s_odd;
	scode_table[generic_code_data::table::EPEER][(s_odd-1)/2] = s_even;

	if (type == generic_code_data::TYPE1)
		scode_table[generic_code_data::table::TYPE][s_crossing] = generic_code_data::TYPE1;
	else
		scode_table[generic_code_data::table::TYPE][s_crossing] = generic_code_data::TYPE2;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "assign_satellite_code_data:       crossing is " << (scode_table[generic_code_data::table::TYPE][s_crossing] == generic_code_data::TYPE1? "Type I" : "Type II") << endl;
					
	/* crossing indicates whether we have a positive, negative crossing in the standard sense, so we adjust
	   for the peer code based on type.
	*/
	if (crossing == generic_code_data::POSITIVE)
	{
		if (scode_table[generic_code_data::table::TYPE][s_crossing] == generic_code_data::TYPE2)
			scode_table[generic_code_data::table::LABEL][s_crossing] = generic_code_data::POSITIVE;
		else
			scode_table[generic_code_data::table::LABEL][s_crossing] = generic_code_data::NEGATIVE;
	}
	else if (crossing == generic_code_data::NEGATIVE)
	{
		if (scode_table[generic_code_data::table::TYPE][s_crossing] == generic_code_data::TYPE1)
			scode_table[generic_code_data::table::LABEL][s_crossing] = generic_code_data::POSITIVE;
		else
			scode_table[generic_code_data::table::LABEL][s_crossing] = generic_code_data::NEGATIVE;
	}
	else
	{
		scode_table[generic_code_data::table::LABEL][s_crossing] = crossing; // virtual, or flat
	}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "assign_satellite_code_data:       label is ";
	if (scode_table[generic_code_data::table::LABEL][s_crossing] == generic_code_data::POSITIVE)
		debug << "positive" << endl;
	else if (scode_table[generic_code_data::table::LABEL][s_crossing] == generic_code_data::NEGATIVE)
		debug << "negative" << endl;
	else if (scode_table[generic_code_data::table::LABEL][s_crossing] == generic_code_data::VIRTUAL)
		debug << "virtual" << endl;
	else if (scode_table[generic_code_data::table::LABEL][s_crossing] == generic_code_data::FLAT)
		debug << "flat" << endl;
	else
	    debug << scode_table[generic_code_data::table::LABEL][s_crossing] << endl;
}
					
	scode_table[generic_code_data::table::EVEN_TERMINATING][s_crossing] = s_even;
	scode_table[generic_code_data::table::ODD_TERMINATING][s_crossing] = s_odd;				


	if (s_odd == s_edge)
	{
		if ( s_odd == satellite_data.first_edge_on_component[s_component] + satellite_data.num_component_edges[s_component] - 1)
			scode_table[generic_code_data::table::EVEN_ORIGINATING][s_crossing] = satellite_data.first_edge_on_component[s_component];
		else
			scode_table[generic_code_data::table::EVEN_ORIGINATING][s_crossing] = s_odd+1;

		scode_table[generic_code_data::table::COMPONENT][s_crossing] = s_prime_component;
	}
	else
	{
		if ( s_odd == satellite_data.first_edge_on_component[s_prime_component] + satellite_data.num_component_edges[s_prime_component] - 1)
			scode_table[generic_code_data::table::EVEN_ORIGINATING][s_crossing] = satellite_data.first_edge_on_component[s_prime_component];
		else
			scode_table[generic_code_data::table::EVEN_ORIGINATING][s_crossing] = s_odd+1;					

		scode_table[generic_code_data::table::COMPONENT][s_crossing] = s_component;
	}
	
	scode_table[generic_code_data::table::ODD_ORIGINATING][s_crossing] = s_even+1;
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "assign_satellite_code_data:       even terminating edge is " << scode_table[generic_code_data::table::EVEN_TERMINATING][s_crossing] << endl;
	debug << "assign_satellite_code_data:       odd terminating edge is " << scode_table[generic_code_data::table::ODD_TERMINATING][s_crossing] << endl;
	debug << "assign_satellite_code_data:       even originating edge is " << scode_table[generic_code_data::table::EVEN_ORIGINATING][s_crossing] << endl;
	debug << "assign_satellite_code_data:       odd originating edge is " << scode_table[generic_code_data::table::ODD_ORIGINATING][s_crossing] << endl;
	debug << "assign_satellite_code_data:       component " << s_component << endl;
}

}

/* gauss_parity returns a vector of length the number of immersion crossings in code data, whose entries
   corresponding to non-virtual and non-shortcut crossings are set to gauss_orientation_data::parity::EVEN
   or gauss_orientation_data::parity::ODD dependent on the number of intervening terms in the Gauss data.
*/
vector<int> gauss_parity(generic_code_data& code_data)
{
	int num_crossings = code_data.num_crossings;
	gauss_orientation_data gauss_data(code_data);
	
if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_parity: gauss_crossing = " << i+1 << " immersion_crossing = " << immersion_crossing << endl;
	
		int start=0;
		for (int j=0; j< num_terms; j++)
		{
			if (orientation_matrix[1][j] == i+1)
			{
				start=j;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_parity:   first ocurrence at index " << start << endl;
				break;
			}
		}
		
		for (int k=1; k< num_terms; k++)
		{
			if (orientation_matrix[1][start+k] == i+1)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_parity:   second ocurrence at index " << start+k << endl;
				if (k%2)
					parity[immersion_crossing] = gauss_orientation_data::parity::EVEN;
				else
					parity[immersion_crossing] = gauss_orientation_data::parity::ODD;
				
				break;
			}
		}
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
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

#ifdef TAKEOUT
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
					   
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
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
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
  	debug << " is valid, testing resultant renumbering" << endl;

			generic_code_data renumbered_code_data(peer_code_data);						
			renumber_peer_code(renumbered_code_data, shift);
			
			matrix<int>& mcode_table = minimum_peer_code_data.code_table;
			matrix<int>& rcode_table = renumbered_code_data.code_table;
							
							
			bool new_minimal_rep_found = false;
			/* check the lexicographic ordering of the peer code, ignoring components, that is
			   compare the product of the generic_code_data::table::TYPE and generic_code_data::table::OPEER
			*/
			for (int i=0; i< num_crossings; i++)
			{
				if (rcode_table[generic_code_data::table::TYPE][i]*rcode_table[generic_code_data::table::OPEER][i] < mcode_table[generic_code_data::table::TYPE][i]*mcode_table[generic_code_data::table::OPEER][i])
				{
					new_minimal_rep_found = true;
					break;
				}
				else if (rcode_table[generic_code_data::table::TYPE][i]*rcode_table[generic_code_data::table::OPEER][i] > mcode_table[generic_code_data::table::TYPE][i]*mcode_table[generic_code_data::table::OPEER][i])
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
#endif

/* manturov_alexander implements the Alexander-like polynomial for flat virtual knots described in 
   V.O.Manturov, I.M.Nikonov. Flat-virtual knot: introduction and some invariants 	arXiv:2406.12864

   The current version only supports knots, not links
*/
void manturov_alexander(generic_code_data& code_data)
{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "manturov_alexander: presented with code_data ";
	write_code_data(debug, code_data);
	debug << endl;
}

	matrix<int>& code_table = code_data.code_table;			
	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;

	int num_classical_crossings = num_crossings;
	for (int i=0; i< num_crossings; i++)
	{
		if ((code_table[generic_code_data::table::LABEL][i] != generic_code_data::label::POSITIVE) &&
		    (code_table[generic_code_data::table::LABEL][i] != generic_code_data::label::NEGATIVE))
			num_classical_crossings--;
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "manturov_alexander: num_classical_crossings = " << num_classical_crossings << endl;
	
	/* for each edge determine the weight according to 
		
		classical crossing egress edges weight 1
		virtual crossings change weight by s ^{+-1}
	*/
	vector<polynomial<int> > edge_weight(num_edges);
	/* find the first originating under-arc edge at a classical crossing, we shall be numbering the classical crossings from the
	   origin crossing of the start edge.
	*/
	vector<int> classical_crossing_number(num_crossings,-1);
	int classical_crossing_index = 0;
	int start = 0;
	do
	{
		int origin_crossing = code_data.orig_crossing[start];
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "manturov_alexander: start = " << start << ", origin_crossing = " << origin_crossing << endl;

		if ((start%2 == 0 && code_table[generic_code_data::table::LABEL][origin_crossing] == generic_code_data::label::POSITIVE) ||
		    (start%2 == 1 && code_table[generic_code_data::table::LABEL][origin_crossing] == generic_code_data::label::NEGATIVE))
		{
			/* start is an even under-arc egress edge at a classical crossing */			
			classical_crossing_number[origin_crossing] = classical_crossing_index++;
			break;
		}
		else
		{
			start++;
		}
	} while(true);
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "manturov_alexander: start edge = " << start << endl;

	
	polynomial<int> weight("1");
	
	vector<int> over_terminating_arc(num_crossings,-1);
	vector<int> under_terminating_arc(num_crossings,-1);
	
	int arc_index = 0;
	
	for (int i=0; i< num_edges; i++)
	{
		int edge = (start+i)%num_edges;
		edge_weight[edge] = weight;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "manturov_alexander: edge = " << edge << ", weight " << edge_weight[edge] << endl;
		
		/* evaluate the weight of the next edge */
		int next_crossing = code_data.term_crossing[edge];
		if ((edge%2 == 1 && code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::label::POSITIVE) ||
		    (edge%2 == 0 && code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::label::NEGATIVE))
		{
			/* arrived on the under_arc of a classical crossing */
			weight = polynomial<int>("1");
			under_terminating_arc[next_crossing] = arc_index;
			arc_index++;
			
			if (classical_crossing_number[next_crossing] == -1)
				classical_crossing_number[next_crossing] = classical_crossing_index++;
		}
		else if ((edge%2 == 0 && code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::label::POSITIVE) ||
		    (edge%2 == 1 && code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::label::NEGATIVE))
		{
			/* arrived on the over arc  of a classical crossing */
			over_terminating_arc[next_crossing] = arc_index;

			if (classical_crossing_number[next_crossing] == -1)
				classical_crossing_number[next_crossing] = classical_crossing_index++;

		}
		else if (code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::label::VIRTUAL)
		{
			if ((edge%2 == 0 && code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::type::TYPE1) ||
			    (edge%2 == 1 && code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::type::TYPE2))
			{
				weight *= polynomial<int>("v^-1");
			}
			else
			{
				weight *= polynomial<int>("v");
			}			
		}		
		else if (code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::label::FLAT)
		{
			if ((edge%2 == 0 && code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::type::TYPE1) ||
			    (edge%2 == 1 && code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::type::TYPE2))
			{
				weight *= polynomial<int>("u^-1");
			}
			else
			{
				weight *= polynomial<int>("u");
			}			
		}		
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "manturov_alexander: edge_weight: ";
	for (int i=0; i< num_edges; i++)
		debug << edge_weight[i] << "  ";
	debug << endl;
	
	debug << "manturov_alexander: over_terminating_arc: ";
	for (int i=0; i< num_crossings; i++)
		debug << over_terminating_arc[i] << ' ';
	debug << endl;
	
	debug << "manturov_alexander: under_terminating_arc: ";
	for (int i=0; i< num_crossings; i++)
		debug << under_terminating_arc[i] << ' ';
	debug << endl;

	debug << "manturov_alexander: classical_crossing_number: ";
	for (int i=0; i< num_crossings; i++)
		debug << classical_crossing_number[i] << ' ';
	debug << endl;
}
	
	
	matrix<polynomial<int> > M(num_classical_crossings,num_classical_crossings);
	int classical_crossing = 0;
	for (int i=0; i<num_crossings; i++)
	{
		if ((code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::POSITIVE) ||
		   (code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::NEGATIVE))
		{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "manturov_alexander: crossing " << i << " is classical" << endl;


			int ingress_over_edge;
			int ingress_under_edge;
			
			if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::POSITIVE)
			{
				ingress_over_edge = code_table[generic_code_data::table::EVEN_TERMINATING][i];
				ingress_under_edge = code_table[generic_code_data::table::ODD_TERMINATING][i];
			}
			else
			{
				ingress_under_edge = code_table[generic_code_data::table::EVEN_TERMINATING][i];
				ingress_over_edge = code_table[generic_code_data::table::ODD_TERMINATING][i];
			}
		
			bool positive_crossing = true;
			if ((code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::POSITIVE &&
			    code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE1) ||
			    (code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::NEGATIVE &&
			    code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE2))
			{
				positive_crossing = false;
			}

			int crossing = classical_crossing_number[classical_crossing];

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "manturov_alexander:   positive_crossing = " << positive_crossing << ",   classical_crossing number = " << crossing << endl;
			
			
			M[crossing][crossing] = polynomial<int>("-1");
			M[crossing][under_terminating_arc[i]] += edge_weight[ingress_under_edge]*(positive_crossing?polynomial<int>("t"):polynomial<int>("t^-1"));
			M[crossing][over_terminating_arc[i]] += edge_weight[ingress_over_edge]*(positive_crossing?polynomial<int>("1-t"):polynomial<int>("1-t^-1"));
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "manturov_alexander:   ingress_under_edge = " << ingress_under_edge << ", weight = " << edge_weight[ingress_under_edge] << endl;
	debug << "manturov_alexander:   ingress_over_edge = " << ingress_over_edge << ", weight = " << edge_weight[ingress_over_edge] << endl;
	debug << "manturov_alexander:   under_terminating_arc = " << under_terminating_arc[i] << ", over_terminating_arc = " << over_terminating_arc[i] << endl;
}

		
			classical_crossing++;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "manturov_alexander: crossing " << i << " is not classical" << endl;
		}
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "manturov_alexander: M = " << endl;
	print(M,debug,0,"manturov_alexander: ");
}	
	
	if (num_classical_crossings > 1)
	{
		polynomial<int> det = determinant(M);
		cout << "Manturov Alexander polynomial = " << det << endl;
	}
	else
	{
		cout << "Manturov Alexander polynomial = " << M[0][0] << endl;
	}
	
}

polynomial<int> nabla_k(generic_code_data& code_data, int starred_edge);
//void nabla_k(generic_code_data& code_data, int starred_edge);

void mock_alexander(generic_code_data& code_data)
{

	bool pure_knotoid = (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID);

	if (code_data.type != generic_code_data::code_type::peer_code || (code_data.immersion_character != generic_code_data::character::PURE_KNOTOID && code_data.immersion_character != generic_code_data::character::KNOTOID))
	{
		if (!braid_control::SILENT_OPERATION)
			cout << "mock_alexander requires the peer code of a knotoid, skipping" << endl;

		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << "mock_alexander requires the peer code of a knotoid, skipping" << endl;
		}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "mock_alexander: mock_alexander requires the peer code of a knotoid, doing nothing" << endl;

		return;
	}
	else
	{
		if(!pure_knotoid || valid_knotoid_input(code_data))  // only need to check knotoid data if it's a pure knotoid
		{	
/*				
			if (CYCLE_KNOT_TYPE_NABLA_K)
			{
				polynomial<int> nabla_k_main = nabla_k(code_data,0);

				for (int i=1; i< 2*code_data.num_crossings; i++)
				{
					polynomial<int> nabla_k_i = nabla_k(code_data,i);
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "mock_alexander: nabla_k_" << i << " = " << nabla_k_i << endl;
	
					if (nabla_k_i != nabla_k_main)
					{
						cout << "nabla_k_" << i << " = " << nabla_k_i << " nabla_k_main = " << nabla_k_main << endl; 
						exit(0);
					}
				}
			}
*/
			polynomial<int> nabla_k_main = nabla_k(code_data,0);

			if (!braid_control::SILENT_OPERATION)
				cout << "mock_alexander polynomial = " << nabla_k_main << endl;
	
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "mock_alexander polynomial = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << '\n';
			}
			output << nabla_k_main << endl;
		}
		else
		{
	
			if (!braid_control::SILENT_OPERATION)
				cout << "\nmock_alexander: invalid knotoid input, skipping" << endl;
	
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "mock_alexander: invalid knotoid input, skipping" << endl;
			}
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "nabla_k: invalid knotoid input" << endl;
	
			return;
		}
	}
}	
	
polynomial<int> nabla_k(generic_code_data& code_data, int starred_edge)	
{
	bool pure_knotoid = (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID);
	
//	if(!pure_knotoid || valid_knotoid_input(code_data))  // only need to check knotoid data if it's a pure knotoid
	{		
		matrix<int>& code_table = code_data.code_table;			
		int head = code_data.head;	
		int num_crossings = code_data.num_crossings;
		int num_edges = 2*num_crossings;
		int num_cycles;
		int num_left_cycles;

		int num_classical_crossings = num_crossings;
		if (pure_knotoid)
		{
			for (int i=0; i< num_crossings; i++)
			{
				if (code_data.shortcut_crossing[i] != 0)
					num_classical_crossings--;
			}
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "nabla_k: starred_edge = " << starred_edge << endl;
	debug << "nabla_k: num_classical_crossings = " << num_classical_crossings << endl;
}
	
		vector<int> crossing_sign(num_crossings);
		
		for (int i=0; i< num_crossings; i++)
		{
			if (code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL)
			{
				if ((code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1 && code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE) ||
				    (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE2 && code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE))
					crossing_sign[i] = 1;
				else
					crossing_sign[i] = -1;
			}
		}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "nabla_k: crossing_sign: ";
	for (int i=0; i< num_crossings; i++)
		debug << crossing_sign[i] << ' ';
	debug << endl;
}

		matrix<int> cycle(num_crossings+2, num_edges+1);
		calculate_turning_cycles(code_data, cycle, num_left_cycles, num_cycles);

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "nabla_k: num_cycles = " << num_cycles << endl;
    debug << "nabla_k: number of left turning cycles = " << num_left_cycles;
	for (int i=0; i<num_left_cycles; i++)
	{
		debug << "\nnabla_k: cycle " << i << " length = " << cycle[i][0] << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " ";
	}
	debug << endl;
    debug << "nabla_k: number of right turning cycles = " << num_cycles - num_left_cycles;
	for (int i=num_left_cycles; i<num_cycles; i++)
	{
		debug << "\nnabla_k: cycle " << i << " length = " << cycle[i][0] << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " "	;
	}
	debug << endl;
}

		/* identify the region of the knotoid complement corresponding to each turning cycle. turning cycles sharing a shortcut
		   edge belong to the same region.
		*/
		vector<int> region(num_cycles);		
		for (int i=0; i< num_cycles; i++)
			region[i] = i;
			
		int first_shortcut_edge = starred_edge;  // set for the case of knot-type knotoids
		
		if (pure_knotoid)
		{
			if (code_table[generic_code_data::table::LABEL][head] == generic_code_data::POSITIVE)
				first_shortcut_edge = code_table[generic_code_data::table::OPEER][head];
			else if (code_table[generic_code_data::table::LABEL][head] == generic_code_data::NEGATIVE)
				first_shortcut_edge = 2*head;
			else
			{
				/* shouldn't get here as valid_knotoid_input checks the first shortcut crossing cannot be virtual */
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "nabla_k: head " << 0 << " indicates first shortcut crosing is virtual" << endl;
			}
		}
	
		int head_cycle = -1; // the first turning cycle containing the first shortcut edge.
		int first_starred_cycle = -1; // the first turning cycle containing the starred edge.
		int second_starred_cycle = -1; // the first turning cycle containing the starred edge.
		
		for (int i = first_shortcut_edge; i<= (pure_knotoid? code_data.num_component_edges[0]: starred_edge); i++)
		{
			int shortcut_edge = i%code_data.num_component_edges[0]; // allows us to cycle to include edge zero, is equal to zero in the case of knot-type knotoids

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "nabla_k: shortcut_edge " << shortcut_edge;
			
			/* identify first occurrence of the shortcut edge and note the first turning cycle containing the first edge
			   of the shortcut, since this will determine the region containing the head of the knotoid.
			*/
			int first_cycle = -1;
			int second_cycle;
			bool finished = false;
			for (int j=0; j< num_cycles; j++)
			{
				for (int k=1; k<= cycle[j][0]; k++)
				{
					if (abs(cycle[j][k]) == shortcut_edge)
					{
						if (first_cycle == -1)
						{
							first_cycle = j;							
							if (i == first_shortcut_edge)
								head_cycle = first_cycle;

							if (i == starred_edge)
								first_starred_cycle = first_cycle;
						}
						else
						{
							second_cycle = j;
							
							if (i == starred_edge)
								second_starred_cycle = second_cycle;
								
							finished = true;
						}
						
						break;
					}
				}
				if(finished)
					break;
			}
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << " first_cycle " << first_cycle << " second_cycle = " << second_cycle << endl;

			region[second_cycle] = first_cycle;
		}	

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "nabla_k: head_cycle =  " << head_cycle << endl;

		int index = 0;
		for (int i=0; i< num_cycles; i++)
		{
			if(region[i] == i)
				region[i] = index++;
			else
				region[i] = region[region[i]];
		}
		
		int num_regions = index;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "nabla_k: num_regions = " << num_regions << endl;
	debug << "nabla_k: cycle to region mapping: ";
	for (int i=0; i< num_cycles; i++)
		debug << region[i] << ' ';
	debug << endl;
	debug << "nabla_k: first_starred_cycle = " << first_starred_cycle << endl;
	debug << "nabla_k: second_starred_cycle = " << second_starred_cycle << endl;
	debug << "nabla_k: region corresponding to starred_cycle = " << region[first_starred_cycle] << endl;
}

		
		matrix<polynomial<int> > M(num_classical_crossings, num_regions);
		
		index = -1;
		for (int i=0; i< num_crossings; i++)
		{
			if (!pure_knotoid || code_data.shortcut_crossing[i] == 0)
			{
				index++;
//if (debug_control::DEBUG >= debug_control::SUMMARY)
//	debug << "nabla_k: crossing " << crossing << " index = " << index << endl;

				int edge_1;
				int edge_2;
				int turning_cycle;
				int crossing_region;
				
				edge_1 = code_table[generic_code_data::table::EVEN_TERMINATING][i];
				edge_2 = code_table[generic_code_data::table::ODD_TERMINATING][i];					
				turning_cycle = find_cycle(cycle, num_cycles, num_left_cycles, edge_1, edge_2, (code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE1?true:false));
				crossing_region = region[turning_cycle];					
				M[index][crossing_region] += (crossing_sign[i] == 1? polynomial<int>("-B") : polynomial<int>("-W"));
				
				edge_2 = code_table[generic_code_data::table::EVEN_ORIGINATING][i];
				turning_cycle = find_cycle(cycle, num_cycles, num_left_cycles, edge_1, edge_2, (code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE1?false:true)); 
				crossing_region = region[turning_cycle];					
				M[index][crossing_region] += polynomial<int>("1");
				
				edge_1 = code_table[generic_code_data::table::ODD_ORIGINATING][i];
				turning_cycle = find_cycle(cycle, num_cycles, num_left_cycles, edge_1, edge_2, (code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE1?true:false)); 
				crossing_region = region[turning_cycle];					
				M[index][crossing_region] += (crossing_sign[i] == 1? polynomial<int>("W") : polynomial<int>("B"));

				edge_2 = code_table[generic_code_data::table::ODD_TERMINATING][i];
				turning_cycle = find_cycle(cycle, num_cycles, num_left_cycles, edge_1, edge_2, (code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE1?false:true)); 
				crossing_region = region[turning_cycle];					
				M[index][crossing_region] += polynomial<int>("1");					
			}
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "nabla_k: M = " << endl;
	for (int i=0; i< num_classical_crossings; i++)
	{
		debug << "nabla_k:     ";
		for (int j=0; j< num_regions; j++)
			debug << M[i][j] << ' ';
		debug << endl;
	}
		
	debug << endl;
}
		
		/* nabla_k is determined by placing a star adjacent to the leg of the knotoid.  In a pure knotoid, the leg is always edge zero, and edge zero always appears
		   in the first turning cycle, the star lies in region zero, so to evaluate nabla_ we omit the first column of M to calculate the permanent
		   
		   In a knot-type knotoid, we're being given a peer code of a knot and could consider the leg to be on any edge
		*/ 
		int r_perm[num_classical_crossings];
		for (int i=0; i< num_classical_crossings; i++)
			r_perm[i] = i;
			
		int c_perm[num_classical_crossings];
//		for (int i=0; i< num_classical_crossings; i++)
//			c_perm[i] = i+1;			

		index = 0;
		for (int i=0; i< num_regions; i++)
		{
			if (i == region[first_starred_cycle])
				continue;
				
			c_perm[index++] = i;			
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "nabla_k: main M = " << endl;
	for (int i=0; i< num_classical_crossings; i++)
	{
		debug << "nabla_k:     ";
		for (int j=0; j< num_classical_crossings; j++)
			debug << M[i][c_perm[j]] << ' ';
		debug << endl;
	}
		
	debug << endl;
}
		
		polynomial<int> nabla_k_main = permanent(M, "", num_classical_crossings, r_perm, c_perm);

/*
		if (!braid_control::SILENT_OPERATION)
			cout << "nabla_k polynomial = " << nabla_k_main << endl;

		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << "nabla_k polynomial = ";
			if (braid_control::OUTPUT_AS_INPUT)
				output << '\n';
		}
		output << nabla_k_main << endl;
*/

		/* nabla_k_prime is determined by placing a star adjacent to the head of the knotoid.  For a pure knotoid, the region needs to be looked up using 
		   head_cycle so we omit the corresponding column of M to calculate the permanent, for a knot-type knotoid, nabla_k_prime = nabla_k
		*/ 
		
		if (pure_knotoid)
		{
			index = 0;
			for (int i=0; i< num_regions; i++)
			{
				if (i == region[head_cycle])
					continue;
					
				c_perm[index++] = i;			
			}
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "nabla_k: prime M = " << endl;
	for (int i=0; i< num_classical_crossings; i++)
	{
		debug << "nabla_k:     ";
		for (int j=0; j< num_classical_crossings; j++)
			debug << M[i][c_perm[j]] << ' ';
		debug << endl;
	}
		
	debug << endl;
}

		if (braid_control::WAIT_SWITCH)
		{
		    matrix_control::WAIT_INFO = true;
		    matrix_control::wait_threshold = braid_control::wait_threshold;
		    matrix_control::wait_count = 0;
		}
		
		polynomial<int> nabla_k_prime = permanent(M, "", num_classical_crossings, r_perm, c_perm);
		
		polynomial<int> step_1 = substitute_signed_polynomial_variable(nabla_k_prime,'W','A',true);
		polynomial<int> step_2 = substitute_signed_polynomial_variable(step_1,'B','W',true);
		polynomial<int> step_3 = substitute_signed_polynomial_variable(step_2,'A','B',false);

			

		if (step_3 != nabla_k_main)
		{
			cout << "nabla_k = " << nabla_k_main << " nabla_k_prime = " << nabla_k_prime << " step_3 = " << step_3 << endl; 
			exit(0);
		}
		else
		{
			return nabla_k_main;
		}
		
	}
}

int find_cycle(matrix<int>& cycle, int num_cycles, int num_left_cycles, int edge_1, int edge_2, bool left_cycle)
{
	int turning_cycle = -1;
	
	int start = (left_cycle? 0: num_left_cycles);
	int end = (left_cycle? num_left_cycles: num_cycles);
	bool found = false;
	
	for (int i=start; i< end; i++)
	{
		for (int j=1; j<= cycle[i][0]; j++)
		{
			if (abs(cycle[i][j]) == edge_1)
			{
				int succeeding_edge = (j<cycle[i][0]? j+1 : 1);
				if (abs(cycle[i][succeeding_edge]) == edge_2)
				{
					turning_cycle = i;
					found = true;
					break;
				}
			}
		}
		if (found)
			break;
	}
	
	return turning_cycle;
}

void doodle_Q_polynomial(generic_code_data code_data)
{
	int num_crossings = code_data.num_crossings;
	int num_components = code_data.num_components;
	
	if (num_components > 1)
	{
		cout << "doodle_Q_polynomial: too many component in input" << endl;
		return;
	}

	bool virtual_diagram = false;
	bool doodle_diagram = true;
	
	for (int i=0; i< num_crossings; i++)
	{
		if (code_data.code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::VIRTUAL)
			virtual_diagram = true;
		else if (code_data.code_table[generic_code_data::table::LABEL][i] != generic_code_data::label::FLAT)
		{
			doodle_diagram = false;
			break;
		}
	}
	
	if (!doodle_diagram)
	{
		cout << "doodle_Q_polynomial: not a doodle diagram" << endl;
		return;
	}

	if (!virtual_diagram)
	{
		cout << "doodle_Q_polynomial: not a virtual diagram" << endl;
		return;
	}

	polynomial<int> QD_t;
	
	/* Smooth each flat crossing to give a pair of components, Dp */
	for (int p_crossing=0; p_crossing< num_crossings; p_crossing++)
	{
		if (code_data.code_table[generic_code_data::table::LABEL][p_crossing] != generic_code_data::label::VIRTUAL)
		{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "doodle_Q_polynomial: smooth diagram at p_crossing " << p_crossing << endl;

			generic_code_data p_smoothed_code_data;
			bool shifted_start = smooth_diagram(code_data,p_crossing,p_smoothed_code_data);

			/* evaluate \sigma_D(p), the linking number based on flat crossings of the two components in Dp */
			int sigma_Dp = linking_number(p_smoothed_code_data,generic_code_data::label::FLAT,0,1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "doodle_Q_polynomial:   p_smoothed_code_data = ";
	write_peer_code(debug, p_smoothed_code_data);
	debug << endl;
	debug << "doodle_Q_polynomial:   sigma_Dp = " << sigma_Dp << endl;	
}
	
			polynomial<int> p_term;
						
			/* evaluate the sign of each component self intersection q_crossing from p_smoothed_code_data component 
			zero of p_smoothed_code_data is the upper component in the original code_data at p_crossing
			*/
			vector<int> q_sign(p_smoothed_code_data.num_crossings);
			for (int i=0; i< 2; i++)
			{
				bool odd_component;
				if ( (i==0 && code_data.code_table[generic_code_data::table::TYPE][p_crossing] == generic_code_data::type::TYPE1) || 
				     (i==1 && code_data.code_table[generic_code_data::table::TYPE][p_crossing] == generic_code_data::type::TYPE2)
				   )
					odd_component = true;
				else
					odd_component = false;
				
				/* trace the component from the starting point to identify the edge on which we first encounter the q_crossing,
				   the start is the first edge of the component unless we're considering the odd component and we had to adjust 
				   the edge labels when calculating p_smoothed_code_data.
				*/
				int start = p_smoothed_code_data.first_edge_on_component[i];
				
				if (odd_component && shifted_start)
					start++;
					
				for (int j=0; j< p_smoothed_code_data.num_component_edges[i]; j++)
				{
					int edge = (start+j-p_smoothed_code_data.first_edge_on_component[i])%p_smoothed_code_data.num_component_edges[i] + p_smoothed_code_data.first_edge_on_component[i];
					int crossing = p_smoothed_code_data.term_crossing[edge];
					int peer = (edge%2? p_smoothed_code_data.code_table[generic_code_data::table::EVEN_TERMINATING][crossing]:p_smoothed_code_data.code_table[generic_code_data::table::ODD_TERMINATING][crossing]);
					int peer_component = (edge%2? p_smoothed_code_data.code_table[generic_code_data::table::COMPONENT][peer/2]:p_smoothed_code_data.code_table[generic_code_data::table::COMPONENT][(peer-1)/2]);
					
					if (p_smoothed_code_data.code_table[generic_code_data::table::LABEL][crossing] == generic_code_data::label::FLAT && peer_component == i && q_sign[crossing] == 0) // first visit to a component flat self intersection
					{
						if ( (edge%2==1 && p_smoothed_code_data.code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::type::TYPE1) || 
						     (edge%2==0 && p_smoothed_code_data.code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::type::TYPE2)
						   )
							q_sign[crossing] = -1; // second visit arrives from the right
						else
							q_sign[crossing] = 1;
					}
				}				
			}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "doodle_Q_polynomial:   q_sign: ";	
	for (int i=0; i< p_smoothed_code_data.num_crossings; i++)
		debug << q_sign[i] << ' ';
	debug << endl;
}
					
			/* Consider each flat crossing q in the two components of Dp */
			for (int i=0; i< 2; i++)
			{
				vector<int> crossing_map;				
				generic_code_data Dp_component = isolate_component(p_smoothed_code_data,i,crossing_map);

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "doodle_Q_polynomial:   Dp_component " << i << ": ";	
	write_code_data(debug,Dp_component);
	debug << endl;
}			
				int num_Dp_crossings = Dp_component.num_crossings;
				
				for (int q_crossing=0; q_crossing < num_Dp_crossings; q_crossing++)
				{
					if (Dp_component.code_table[generic_code_data::table::LABEL][q_crossing] == generic_code_data::label::FLAT)
					{												
						generic_code_data q_smoothed_code_data;
						smooth_diagram(Dp_component,q_crossing,q_smoothed_code_data);
						
						/* evaluate the weight of D_p(q) */
						int weight_Dpq;
						ostringstream oss;
						
						if (q_sign[crossing_map[q_crossing]] == 1)
						{
							weight_Dpq = linking_number(q_smoothed_code_data,generic_code_data::label::VIRTUAL,0,1);
							if (weight_Dpq == 1)
								oss << "t-1";
							else
								oss << "t^" << weight_Dpq << "-1";
						}
						else
						{
							weight_Dpq = linking_number(q_smoothed_code_data,generic_code_data::label::VIRTUAL,1,0);
							if (weight_Dpq == 1)
								oss << "1-t";
							else
								oss << "1-t^" << weight_Dpq;
						}
						
						polynomial<int> q_term(oss.str());

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "doodle_Q_polynomial:     q_sign = " << q_sign[crossing_map[q_crossing]] << ", weight_Dpq = " << weight_Dpq << ", q_term = " << q_term << endl;	
						
						p_term += q_term;
					}
				}
			}			
			
			p_term *= polynomial<int>(sigma_Dp);
			QD_t += p_term;
		}
	}
	
	cout << "doodle_Q polynomial = " << QD_t << endl;
}

/* isolate component creates a peer code for the specified unicursal component in code_data, considering self intersections of the component only.
   The function returns to crossing_map the mapping from the isolated component crossings to those of code_data.
*/
generic_code_data isolate_component(generic_code_data& code_data, int component,vector<int>& crossing_map)
{
	int first_old_edge_on_component = code_data.first_edge_on_component[component];
	int num_old_edges_on_component = code_data.num_component_edges[component];
	int num_crossings_on_component = 0;
	int new_edge = 0;

	/* write new_peers as above and record the new type of the crossing.  We are guaranteed to have an odd and even new_label at each crossing, so if the
	   polarity of the first occurrence of the new labels at a self intersection agrees with that of the old labels, the new_type of the crossing will be
	   as before, otherwise it will be reversed.
	*/
	matrix<int> new_peers(2,code_data.num_crossings,-1);
	vector<int> new_type(code_data.num_crossings);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "isolated_component: first_old_edge_on_component = " << first_old_edge_on_component << ", num_old_edges_on_component = " << num_old_edges_on_component << endl;	
	
	for (int edge=first_old_edge_on_component; edge< first_old_edge_on_component+num_old_edges_on_component; edge++)
	{
		int crossing = code_data.term_crossing[edge];
		int peer = (edge%2? code_data.code_table[generic_code_data::table::EVEN_TERMINATING][crossing]:code_data.code_table[generic_code_data::table::ODD_TERMINATING][crossing]);
		int peer_component = (edge%2? code_data.code_table[generic_code_data::table::COMPONENT][peer/2]:code_data.code_table[generic_code_data::table::COMPONENT][(peer-1)/2]);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "isolated_component: edge " << edge << ", crossing = " << crossing << ", peer = " << peer << ", peer_component = " << peer_component << endl;	

		
		if (peer_component == component)
		{
			if (new_peers[0][crossing] == -1)
			{
				num_crossings_on_component++;
				new_peers[0][crossing] = new_edge;
				
				if (new_edge%2 == edge%2)
					new_type[crossing] = code_data.code_table[generic_code_data::table::TYPE][crossing];
				else
					new_type[crossing] = (code_data.code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::type::TYPE1? 
					                        generic_code_data::type::TYPE2: generic_code_data::type::TYPE1);
			}
			else
			{
				new_peers[1][crossing] = new_edge;
			}

			new_edge++;
		}
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "isolated_component: num_crossings_on_component = " << num_crossings_on_component << endl;
	debug << "isolated_component: new_peers" << endl;
	print(new_peers,debug, 3,"isolated_component: ");
	debug << "isolated_component: new_type: ";
	for (int i=0; i< code_data.num_crossings; i++)
		debug << new_type[i] << ' ';
	debug << endl;
}

	generic_code_data isolated_component;
	isolated_component.type = generic_code_data::code_type::peer_code;
	isolated_component.num_components = 1;
	isolated_component.num_crossings = num_crossings_on_component;
	isolated_component.num_component_edges = {2*num_crossings_on_component};
	isolated_component.code_table = matrix<int>(generic_code_data::table::CODE_TABLE_SIZE,num_crossings_on_component);

	crossing_map = vector<int>(num_crossings_on_component);
	
	for (int i=0; i<code_data.num_crossings; i++)
	{
		if (new_peers[0][i] == -1)
			continue;
			
		int odd_edge = new_peers[0][i];	
		int even_edge = new_peers[1][i];	
		
		if (even_edge%2)
			swap(odd_edge,even_edge);
			
		int new_crossing = even_edge/2;
	
		crossing_map[new_crossing] = i;
		
		isolated_component.code_table[generic_code_data::table::ODD_TERMINATING][new_crossing] = odd_edge;
		isolated_component.code_table[generic_code_data::table::EVEN_TERMINATING][new_crossing] = even_edge;
		isolated_component.code_table[generic_code_data::table::TYPE][new_crossing] = new_type[i];
		isolated_component.code_table[generic_code_data::table::LABEL][new_crossing] = code_data.code_table[generic_code_data::table::LABEL][i];
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "isolated_component: new_code_table" << endl;
	print(isolated_component.code_table,debug, 3,"isolated_component: ");
}

	ostringstream oss;
	write_code_data(oss,isolated_component);
	read_peer_code(isolated_component,oss.str());
			
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "isolated_component: isolated_component " ;	
	write_code_data(debug,isolated_component);
	debug << endl;
	print_code_data(debug,isolated_component,"isolated_component: ");

	debug << "isolated_component: crossing_map: ";	
	for (int i=0; i< num_crossings_on_component; i++)
		debug << crossing_map[i] << ' ';
	debug << endl;
}
	return isolated_component;
}

/* smooth_diagram creates the labelled peer code for the diagram obtained from code_data by Seifert smoothing smoothed_crossing.

   The function starts the labelling of the smoothed diagram from the egress edges of the smoothed_crossing, starting from the
   upper component when viewing the crossing in the standarad left to right orientation.
   
   If the function needs to adjust the edge labels on one of the components to ensure odd and even terminating labels, it adjusts
   the component that starts on the odd egress label.  This will be the case if the two components formed by the smoothing intersect each other
   since one component will have reversed new-label parity compared with the original labelling.
   
   The funtion returns true if the new edge labels needed adjusting.
*/
bool smooth_diagram(generic_code_data& code_data, int smoothed_crossing, generic_code_data& smoothed_code_data)
{
	/* check the number of components, currently we support only single component code_data */
	if (code_data.num_components > 1)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "smooth-diagram: presented with code_data containing more than one component, doing nothing" << endl;
		
		smoothed_code_data = code_data;
		return false;
	}
	
	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;
//	int num_components = code_data.num_components;
	matrix<int>& code_table = code_data.code_table;
	vector<int>& term_crossing = code_data.term_crossing;	

			/* Determine the peer code of the link obtained by Seifert smoothing crossing smoothed_crossing.  We trace the
			   diagram from each egress edge, starting with the upper edge with repect to the standard left to right
			   orientation, until we reach the corresponding ingress edge, recording the new edge labels at the index 
			   corresponding to the crossing number in code_data, the first label is assigned to row zero, the second to row 1.  
			   
			    Type 1                        Type 2
			                
			    2j-1 \   / 2i+1          2i   \   /  2j
			          \ /                      \ /         
			           X                        X   
			          / \                      / \
			    2i   /   \ 2j            2j-1 /   \ 2i+1
			    
			    
			    The first component at a type 1 crossing starts with the odd edge 2i+1 and the second component at the even edge 2j.  
			    For a type 2 crossing the parity is the other way round (2j, then 2i+1).  Starting at an odd old-label with an even new
			    label reverses the polarity of the labels assigned to that component.  If the two components intersect, this means that 
			    some crossings have two even or two odd new ingress labels, in which case we adjust the component that started at the 
			    odd label so we end up with consistent crossing types between the old and new labels.
			    
   			*/
			matrix<int> new_peers(2,num_crossings,-1);
			int new_edge = 0;
			int start;
			int end;

			/* first the upper component */
			if (code_table[generic_code_data::table::TYPE][smoothed_crossing] == generic_code_data::type::TYPE1)
			{
				start = 2*smoothed_crossing+1;
				end = code_table[generic_code_data::table::OPEER][smoothed_crossing];
			}
			else
			{
				start = (code_table[generic_code_data::table::OPEER][smoothed_crossing]+1)%num_edges;				
				end = 2*smoothed_crossing;
			}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "smooth-diagram:   first component start = " << start << ", end = " << end << endl;
			
			for (int i=0; i< num_edges; i++)
			{
				int next_edge = (start+i)%num_edges;
				if ( next_edge == end)
					break;
					
				int crossing = term_crossing[next_edge];
				if (new_peers[0][crossing] == -1)
					new_peers[0][crossing] = new_edge;
				else
					new_peers[1][crossing] = new_edge;
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "smooth-diagram:     next_edge = " << next_edge << ", crossing = " << crossing << ", new_edge = " << new_edge << endl;
				
				new_edge++;
			}

			int num_edges_on_first_component = new_edge;
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "smooth-diagram: num_edges_on_first_component = " << num_edges_on_first_component << endl;

			/* now the lower component */
			if (code_table[generic_code_data::table::TYPE][smoothed_crossing] == generic_code_data::type::TYPE2)
			{
				start = 2*smoothed_crossing+1;
				end = code_table[generic_code_data::table::OPEER][smoothed_crossing];
			}
			else
			{
				start = (code_table[generic_code_data::table::OPEER][smoothed_crossing]+1)%num_edges;				
				end = 2*smoothed_crossing;
			}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "smooth-diagram:   second component start = " << start << ", end = " << end << endl;
			
			for (int i=0; i< num_edges; i++)
			{
				int next_edge = (start+i)%num_edges;
				if (next_edge == end)
					break;
					
				int crossing = term_crossing[next_edge];
				if (new_peers[0][crossing] == -1)
					new_peers[0][crossing] = new_edge;
				else
					new_peers[1][crossing] = new_edge;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "smooth-diagram:     next_edge = " << next_edge << ", crossing = " << crossing << ", new_edge = " << new_edge << endl;

				new_edge++;
			}
			
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "smooth-diagram: new_peers" << endl;
	print(new_peers,debug, 3,"smooth-diagram: ");
}
			bool adjust_component = false;
			for (int i =0; i< num_crossings; i++)
			{
				if (i==smoothed_crossing)
					continue;
				if (new_peers[0][i]%2 == new_peers[1][i]%2)
				{
					adjust_component = true;
					break;
				}
			}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "smooth-diagram: adjust_component = " << adjust_component << endl;
	
			if (adjust_component)
			{
				int start;
				int end;
				int num_edges_on_component;
				int first_component_edge;
				
				if (code_table[generic_code_data::table::TYPE][smoothed_crossing] == generic_code_data::type::TYPE1)
				{
					start = 0;
					end = num_edges_on_first_component-1;
					first_component_edge = 0;
					num_edges_on_component = num_edges_on_first_component;
				}
				else
				{
					start = num_edges_on_first_component;
					end = num_edges-3;
					first_component_edge = num_edges_on_first_component;
					num_edges_on_component = num_edges - 2 - num_edges_on_first_component;
				}
				
				for (int i=0; i<num_crossings; i++)
				{
					if (i==smoothed_crossing)
						continue;
					
					if (new_peers[0][i] >= start && new_peers[0][i] <= end)
						new_peers[0][i] = (new_peers[0][i]+1-first_component_edge)%num_edges_on_component+first_component_edge;
					if (new_peers[1][i] >= start && new_peers[1][i] <= end)
						new_peers[1][i] = (new_peers[1][i]+1-first_component_edge)%num_edges_on_component+first_component_edge;

				}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "smooth-diagram: adjusted new_peers" << endl;
	print(new_peers,debug, 3,"smooth-diagram: ");
}
			}

			smoothed_code_data = code_data;
			smoothed_code_data.num_crossings = num_crossings-1;
			smoothed_code_data.num_components = 2;
			smoothed_code_data.num_component_edges = {num_edges_on_first_component, num_edges - 2 - num_edges_on_first_component};
			matrix<int> new_code_table(generic_code_data::table::CODE_TABLE_SIZE,num_crossings-1);

			for (int i=0; i<num_crossings; i++)
			{
				if (i==smoothed_crossing)
					continue;
					
				int odd_edge = new_peers[0][i];	
				int even_edge = new_peers[1][i];	
				
				if (even_edge%2)
					swap(odd_edge,even_edge);
					
				int new_crossing = even_edge/2;
				
				/* we have recorded the new labels at the offset of the corresponding old crossing, 
				   so can use the index i to determine the TYPE and LABEL of the new crossing
				*/
				new_code_table[generic_code_data::table::ODD_TERMINATING][new_crossing] = odd_edge;
				new_code_table[generic_code_data::table::EVEN_TERMINATING][new_crossing] = even_edge;
				new_code_table[generic_code_data::table::TYPE][new_crossing] = code_table[generic_code_data::table::TYPE][i];
				new_code_table[generic_code_data::table::LABEL][new_crossing] = code_table[generic_code_data::table::LABEL][i];
				new_code_table[generic_code_data::table::COMPONENT][new_crossing] = (even_edge >= num_edges_on_first_component?1:0);
			}

			smoothed_code_data.code_table = new_code_table;
			
			ostringstream oss;
			write_code_data(oss,smoothed_code_data);
			read_peer_code(smoothed_code_data,oss.str());
			
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "smooth-diagram: new_code_table" << endl;
	print(new_code_table,debug, 3,"smooth-diagram: ");
	debug << "smooth-diagram: smoothed_code_data: ";	
	write_code_data(debug,smoothed_code_data);
	debug << endl;
	print_code_data(debug,smoothed_code_data,"smooth-diagram: ");
}

	return adjust_component;
}

int linking_number (generic_code_data& code_data, int crossing_type, int component_1, int component_2)
{

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "linking_number: crossing_type = " << crossing_type << ", component_1 = " << component_1 << ", component_2 = " << component_2 << endl;	

	int link_count = 0;
		
	for (int edge=code_data.first_edge_on_component[component_1]; edge < code_data.first_edge_on_component[component_1] + code_data.num_component_edges[component_1]; edge++)
	{
		int crossing = code_data.term_crossing[edge];
		
		if (code_data.code_table[generic_code_data::table::LABEL][crossing] != crossing_type)
			continue;
		
		int peer = (edge%2? code_data.code_table[generic_code_data::table::EVEN_TERMINATING][crossing]:code_data.code_table[generic_code_data::table::ODD_TERMINATING][crossing]);
		int peer_component = (edge%2? code_data.code_table[generic_code_data::table::COMPONENT][peer/2]:code_data.code_table[generic_code_data::table::COMPONENT][(peer-1)/2]);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "linking_number: edge " << edge << ", crossing = " << crossing << ", peer = " << peer << ", peer_component = " << peer_component << endl;	

		
		if (peer_component == component_2)
		{
			if ((edge%2 == 1 && code_data.code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::type::TYPE1) ||
			    (edge%2 == 0 && code_data.code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::type::TYPE2))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "linking_number:   increment link_count" << endl;	
				link_count++; // component_2 crosses from the right
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "linking_number:   decrement link_count" << endl;	
				link_count--;
			}
		}
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "linking_number: link_count = " << link_count << endl;	
	
	return link_count;
}

/* calculate_colourings is a list based search for colourings of a diagram described by a labelled peer code.  Colourings are represented as a vector of integer
   colours with -1 indicating that no colour has been assigned to the corresponding edge.  A list, candidate_schemes, of potential colourings is  initialized 
   with a candidate indicating each possible colouring of edge zero.  The function then takes the first colouring on the list and attempts to extend the 
   colouring of a diagram at a single crossing.  The function searches the crossings in order for one having a colour assigned to at least one of the ingress edges 
   but not having a colour assigned to one (or both) of the egress edges. 
   
   If such a crossing is found with colour assigned to just one ingress edge, the function considers all possible colourings for the other ingress edge and, 
   for each one, determines the egress edge colours.  If the crossing has colour assigned to both ingress edges, the function simply determines the egress 
   colour for the given edges.  The calculated egress colours are compared to any that have already been assigned.  If the colouring is consistent, the extended 
   colouring is pushed onto the back of the list.  Inconsitent colourings are discarded.

   If no such crossing is found, provided the diagram is connected, every edge will have been assigned a colour, however, it is not necessarily the case that 
   they will be globally consistent, since the colour calculated for an egress edge may not be consistent with colours already assigned at the crossing on which 
   that edge terminates.  An example may be found in the figure 8 knot [-5 7 -1 3]/++++, where the colour calculated for edge three may, in conjunction with the 
   edge calculated earlier for edge 6, conflict with conflict with the colours calculated for edges 4 or 7.  We therefore check each crossing for consistency
   before accepting the colouring as valid.
   
   It is also possible that we may determine the same colouring twice, so the list of completed colouring_schemes is checked for a duplicate before a new
   scheme is added.
*/

//list<vector<int> > calculate_colourings( generic_code_data& code_data, generic_switch_data& switch_data)
list<vector<int> > calculate_colourings(generic_code_data& code_data, matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, matrix<int>& Tu, matrix<int>& Td)
{

	if (braid_control::EXTRA_OUTPUT)
	{
		output << "H: ";
		write_code_data(output, code_data);
		output << endl;
	}
	
	matrix<int>& code_table = code_data.code_table;
	
	vector<int> ingress_top_label(code_data.num_crossings);
	vector<int> ingress_bottom_label(code_data.num_crossings);
	vector<int> egress_top_label(code_data.num_crossings);
	vector<int> egress_bottom_label(code_data.num_crossings);
	vector<int> crossing_type(code_data.num_crossings);
	
	for (int i=0; i< code_data.num_crossings; i++)
	{
		int even_terminating = 2*i;
		int odd_terminating = code_table[generic_code_data::table::ODD_TERMINATING][i];
		int even_originating = code_table[generic_code_data::table::EVEN_ORIGINATING][i];
		int odd_originating = code_table[generic_code_data::table::ODD_ORIGINATING][i];
				
		if (code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE1)
		{
			ingress_bottom_label[i] = even_terminating;
			ingress_top_label[i] = odd_terminating;
			egress_top_label[i] = odd_originating;
			egress_bottom_label[i] = even_originating;
		}
		else
		{
			ingress_bottom_label[i] = odd_terminating;
			ingress_top_label[i] = even_terminating;
			egress_top_label[i] = even_originating;
			egress_bottom_label[i] = odd_originating;
		}
		
		if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::VIRTUAL)
		{
			crossing_type[i] = generic_code_data::label::VIRTUAL;
		}
		else if ((code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::POSITIVE &&
			code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE2) ||
			(code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::NEGATIVE &&
			     code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE1)
		   )
		{
			crossing_type[i] = generic_code_data::label::POSITIVE;
		}
		else 
		{
			crossing_type[i] = generic_code_data::label::NEGATIVE;
		}

	}
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "calculate_colourings: ingress_top_label    ";
	for (int i=0; i< code_data.num_crossings; i++)
		debug << ingress_top_label[i] << ' ';
	debug << endl;
	debug << "calculate_colourings: ingress_bottom_label ";
	for (int i=0; i< code_data.num_crossings; i++)
		debug << ingress_bottom_label[i] << ' ';
	debug << endl;
	debug << "calculate_colourings: egress_top_label     ";
	for (int i=0; i< code_data.num_crossings; i++)
		debug << egress_top_label[i] << ' ';
	debug << endl;
	debug << "calculate_colourings: egress_bottom_label  ";
	for (int i=0; i< code_data.num_crossings; i++)
		debug << egress_bottom_label[i] << ' ';
	debug << endl;
	debug << "calculate_colourings: crossing_type ";
	for (int i=0; i< code_data.num_crossings; i++)
	{
		if (crossing_type[i] == generic_code_data::label::VIRTUAL)		
			debug << "VIRTUAL ";
		if (crossing_type[i] == generic_code_data::label::POSITIVE)		
			debug << "POSITIVE ";
		if (crossing_type[i] == generic_code_data::label::NEGATIVE)		
			debug << "NEGATIVE ";
	}
	debug << endl;
}
	
	list<vector<int> > candidate_schemes;	
	list<vector<int> > colouring_schemes;	
	
	for (size_t i=0; i< Su.numrows(); i++)
	{
		vector<int> colouring(2*code_data.num_crossings,-1);
		colouring[0] = i;
		candidate_schemes.push_back(colouring);
	}
	
	while (candidate_schemes.size() !=0)
	{
//cout << '.' << flush;
		
		vector<int>& colouring = *candidate_schemes.begin();

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "calculate_colourings: candidate colouring ";
	for (int i=0; i< 2*code_data.num_crossings; i++)
		debug << colouring[i] << ' ';
	debug << endl;
}
			
		/* look for a crossing with no egress edge colouring and at least one ingress edge colouring */
		bool found = false;
		
		for (int i=0; i< code_data.num_crossings; i++)
		{
			int even_terminating = 2*i;
			int odd_terminating = code_table[generic_code_data::table::ODD_TERMINATING][i];
			int even_originating = code_table[generic_code_data::table::EVEN_ORIGINATING][i];
			int odd_originating = code_table[generic_code_data::table::ODD_ORIGINATING][i];
			
			bool odd_ingress_edge_coloured = (colouring[odd_terminating] != -1);
			bool even_ingress_edge_coloured = (colouring[even_terminating] != -1);	
			
			/* one egress edge may be coloured if it is an ingress edge to crossing zero */
			if ((colouring[even_originating] == -1 || colouring[odd_originating] == -1) && (even_ingress_edge_coloured || odd_ingress_edge_coloured))
			{
				found = true;
	
				/* determine the number of peer colours we need to test, if there are colours on both ingress edges it is just one,
				   otherwise, we have to test all possible colourings for the uncoloured ingress edge.
				*/
				int num_peer_colours = 1;
				if (!even_ingress_edge_coloured || !odd_ingress_edge_coloured)
					num_peer_colours = Su.numrows();
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "calculate_colourings: extending colouring at crossing " << i << ", num_peer_colours = " << num_peer_colours << endl;
	debug << "calculate_colourings:   even_ingress_edge_coloured = " << even_ingress_edge_coloured << ", odd_ingress_edge_coloured = " << odd_ingress_edge_coloured << endl;
}
			
				for (int j=0; j< num_peer_colours; j++)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_colourings:   peer colour test " << j << endl;

					if (!even_ingress_edge_coloured)
						colouring[even_terminating] = j;
					else if (!odd_ingress_edge_coloured)
						colouring[odd_terminating] = j;				

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "calculate_colourings:   test colouring ";
	for (int i=0; i< 2*code_data.num_crossings; i++)
		debug << colouring[i] << ' ';
	debug << endl;
}
									
					int ingress_top_colour = colouring[ingress_top_label[i]];
					int ingress_bottom_colour = colouring[ingress_bottom_label[i]];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_colourings:     ingress_top_colour = " << ingress_top_colour << " ingress_bottom_colour = " << ingress_bottom_colour << endl;

					int egress_top_colour;
					int egress_bottom_colour;
					
					
					if (crossing_type[i] == generic_code_data::label::VIRTUAL)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_colourings:     virtual crossing" << endl;
						int temp = Td[ingress_bottom_colour][ingress_top_colour];
						egress_top_colour = Tu[ingress_top_colour][ingress_bottom_colour];
						egress_bottom_colour = temp;					
					}
					else if (crossing_type[i] == generic_code_data::label::POSITIVE)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_colourings:     positive crossing" << endl;
						int temp = Sd[ingress_bottom_colour][ingress_top_colour];
						egress_top_colour = Su[ingress_top_colour][ingress_bottom_colour];
						egress_bottom_colour = temp;
					}
					else // negative_crossing)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_colourings:     negative crossing" << endl;
						int temp = invSu[ingress_bottom_colour][ingress_top_colour];
						egress_top_colour = invSd[ingress_top_colour][ingress_bottom_colour];
						egress_bottom_colour = temp;
					}	
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_colourings:     egress_top_colour = " << egress_top_colour << " egress_bottom_colour = " << egress_bottom_colour << endl;
					
					bool consistent_colouring = true;			
					
					/* Chec that the colouring is consistently defined, where there is an egress colouring already defined, also
					   check tha for doubled biracks, the up and down actions are defined, so that the egress colours exits!
					*/
					if (egress_top_colour == -1 || egress_bottom_colour == -1)// Coded as a check, should always have egress_top_colour == -1 iff egress_bottom_colour == -1

						consistent_colouring = false;
					else if(colouring[egress_top_label[i]] != -1 && colouring[egress_top_label[i]] != egress_top_colour)
						consistent_colouring = false;
					else if (colouring[egress_bottom_label[i]] != -1 && colouring[egress_bottom_label[i]] != egress_bottom_colour)
						consistent_colouring = false;
					
					if (consistent_colouring)
					{
						/* assign egress colours and recurse */
						vector<int> extended_colouring = colouring;
						extended_colouring[egress_top_label[i]] = egress_top_colour;
						extended_colouring[egress_bottom_label[i]] = egress_bottom_colour;					
						candidate_schemes.push_back(extended_colouring);											 

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "calculate_colourings:     consistent colouring extended to ";
	for (int i=0; i< 2*code_data.num_crossings; i++)
		debug << extended_colouring[i] << ' ';
	debug << endl;
}

					}
					else
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_colourings:     inconsistent colouring, ignoring" << endl;							
					}
				}						
			}
		}
	
		if (!found)
		{					
			/* check each crossing to make sure the colouring is consistent Necessary since the egress edge consistency check does not take into account that the 
			   terminating crossings of the egress edges may already have colourings that are not consistent with those calculated for a preceeding edge.
			*/
			bool consistent_colouring = true;			

			for (int i=0; i< code_data.num_crossings; i++)
			{
				int egress_top_colour;
				int egress_bottom_colour;
				
				if (crossing_type[i] == generic_code_data::label::VIRTUAL)
				{
					int temp = Td[colouring[ingress_bottom_label[i]]][colouring[ingress_top_label[i]]];
					egress_top_colour = Tu[colouring[ingress_top_label[i]]][colouring[ingress_bottom_label[i]]];
					egress_bottom_colour = temp;					
				}
				else if (crossing_type[i] == generic_code_data::label::POSITIVE)
				{
					int temp = Sd[colouring[ingress_bottom_label[i]]][colouring[ingress_top_label[i]]];
					egress_top_colour = Su[colouring[ingress_top_label[i]]][colouring[ingress_bottom_label[i]]];
					egress_bottom_colour = temp;
				}
				else // negative_crossing)
				{
					int temp = invSu[colouring[ingress_bottom_label[i]]][colouring[ingress_top_label[i]]];
					egress_top_colour = invSd[colouring[ingress_top_label[i]]][colouring[ingress_bottom_label[i]]];
					egress_bottom_colour = temp;
				}									
				
				if(colouring[egress_top_label[i]] != egress_top_colour)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_colourings:     inconsistence at crossing " << i << ": egress_top_colour = " << egress_top_colour << " colouring[egress_top_label[i]] = " << colouring[egress_top_label[i]] << endl;
					consistent_colouring = false;
				}
				else if (colouring[egress_bottom_label[i]] != egress_bottom_colour)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_colourings:     inconsistence at crossing " << i << ": egress_bottom_colour = " << egress_bottom_colour << " colouring[egress_bottom_label[i]] = " << colouring[egress_bottom_label[i]] << endl;
					consistent_colouring = false;
				}
			}
		
			if (consistent_colouring)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_colourings: colouring complete, checking for duplicates...";
	
				if (find(colouring_schemes.begin(),colouring_schemes.end(),colouring) == colouring_schemes.end())
				{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "new scheme found, adding ";
	for (int i=0; i< 2*code_data.num_crossings; i++)
		debug << setw(3) << colouring[i] << ' ';
	debug << "to colouring_schemes" << endl;
}
					
					if (braid_control::EXTRA_OUTPUT)
					{
						output << "H: ";
						for (int i=0; i< 2*code_data.num_crossings; i++)
							output << colouring[i] << ' ';
						output << endl;
					}

					colouring_schemes.push_back(colouring);
					
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "duplicate scheme found, ignoring" << endl;
				}
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_colourings:     inconsistent colouring, ignoring" << endl;							
			}
		}
		
		candidate_schemes.pop_front();	
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "calculate_colourings: found a total of " << colouring_schemes.size() << " colourings" << endl;							

	if (braid_control::EXTRA_OUTPUT)
	{
		if (colouring_schemes.size() == 0)
			output << "H: no colourings" << endl;
	}
	
	return colouring_schemes;
}

/* cocycle_invariant determines a vector of cocycle invariant polynomials, one for each of the cocycles in the switch data, with the
   k-chains determined by the given colourings of the given diagram.
   
   if braid_control::DOUBLE_BIRACKS the k-chains are 3-chains, otherwise they are 2-chains.  The function requires that the switch_data
   contains the appropriate type of cocycle.
*/
//vector<polynomial<int> > cocycle_invariant(generic_code_data& code_data, generic_switch_data& switch_data, list<vector<int> >& colourings)
vector<Cpolynomial> cocycle_invariant(generic_code_data& code_data, generic_switch_data& switch_data, list<vector<int> >& colourings)
{
	int n = switch_data.size;				
	int num_cocycles = switch_data.cocycle_scalar.size();
	int num_chain_generators = switch_data.num_chain_generators;
	matrix<int>& code_table = code_data.code_table;

	/* determine whether the scalar coefficients are a field */
	bool field_coefficients = true;
	if (scalar::variant == scalar::scalar_variant::INT || scalar::variant == scalar::scalar_variant::BIGINT)
		field_coefficients = false;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "cocycle_invariant: presented with switch_data: " << endl;
	print_switch_data(debug, switch_data, "cocycle_invariant:  ");
	debug << "cocycle_invariant: num_cocycles = " << num_cocycles << endl;	
}

	/* each colouring determines a k_chain to which we apply each cocycle, so to avoid working through the diagram crossing more than once we
	   create a matrix k_chain, whose rows will store the three chain for each colouring.  We therefore work through the crossings of the diagram
	   and for each one, determine the contribution to the k_chain for each colouring.
	*/
	matrix<int> k_chain(colourings.size(), num_chain_generators);  
//	vector<polynomial<int> > invariant(num_cocycles);
	vector<Cpolynomial> invariant(num_cocycles);

	/* For testing whether all k-chains are boundaries */
	matrix<int> B;

	
	for (int i=0; i< code_data.num_crossings; i++)
	{
		if (code_table[generic_code_data::table::LABEL][i] != generic_code_data::label::VIRTUAL)
		{							
			/* identify the ingress and egress bottom labels of the crossing */
			int even_terminating = 2*i;
			int odd_terminating = code_table[generic_code_data::table::ODD_TERMINATING][i];
			int even_originating = code_table[generic_code_data::table::EVEN_ORIGINATING][i];
			int odd_originating = code_table[generic_code_data::table::ODD_ORIGINATING][i];
			int	ingress_bottom_label;
			int	egress_bottom_label;
					
			if (code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE1)
			{
				ingress_bottom_label = even_terminating;
				egress_bottom_label = even_originating;
			}
			else
			{
				ingress_bottom_label = odd_terminating;
				egress_bottom_label = odd_originating;
			}

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "cocycle_invariant: crossing i=" << i << " colouring " << " ingress_bottom_label = " << ingress_bottom_label
	      << " egress_bottom_label = " << egress_bottom_label << endl;
}				
						
			bool positive_classical = false;
			
			if ((code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::POSITIVE &&
				code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE2) ||
				(code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::NEGATIVE &&
				     code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE1)
			   )
			{
				positive_classical = true;
			}
			else if ((code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::NEGATIVE &&
				code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE2) ||
				(code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::POSITIVE &&
				     code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE1)
			   )
			{
				positive_classical = false;
			}

			
			int colouring_index = 0;
			list<vector<int> >::iterator lptr =  colourings.begin();
			while (lptr != colourings.end())
			{
				vector<int>& colouring = *lptr;					

				/* identify the colour assigned to the lower two labels of the crossing */
				int under_arc_colour = colouring[ingress_bottom_label];
				int over_arc_colour = colouring[egress_bottom_label];
				
				if (!positive_classical)
					swap(under_arc_colour,over_arc_colour);

				/* The 3-tuple to which the crossing maps is (x,y,z) where x is the right label in both of the pairs
				   on which the twitch acts, y is the left label of the under arc and z the left label of the over-arc.
				*/
				int under_left = under_arc_colour/n;
				int under_right = under_arc_colour%n;
				int over_left = over_arc_colour/n;
				int over_right = over_arc_colour%n;
				int sign = (positive_classical?1:-1);

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "cocycle_invariant: crossing i=" << i << " colouring " << colouring_index << " under_arc_colour = " << under_arc_colour
	      << " over_arc_colour = " << over_arc_colour << endl;
	debug << "cocycle_invariant:     under_left = " << under_left << ", under_right = " << under_right 
	      << ", over_left = " << over_left << ", over_right = " << over_right << endl;
	if (braid_control::DOUBLE_BIRACKS)
		debug << "cocycle_invariant:     tuple = (" << under_right << ',' << under_left << ',' << over_left << ") sign " << sign << endl;
	else
		debug << "cocycle_invariant:     tuple = (" << under_arc_colour << ',' << over_arc_colour << ") sign " << sign << endl;
	
}				
				int tuple;
				if (braid_control::DOUBLE_BIRACKS)
					tuple = n*n*under_right+n*under_left+over_left;
				else
					tuple = n*under_arc_colour+over_arc_colour;

				if (braid_control::BIRACK_HOMOLOGY)
				{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "cocycle_invariant: tuple = " << tuple << endl;				
	
					/* positive crossing map to negative triple points, negative crossings to positive triple points */
					k_chain[colouring_index][tuple] -= sign;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "cocycle_invariant:   updated k_chain[" << colouring_index << "][" << tuple << "] to " << k_chain[colouring_index][tuple] << endl;		
	
				}
				else
				{
	
					int non_degenerate_tuple = switch_data.chain_map[tuple];

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "cocycle_invariant:     tuple = " << tuple << ", non_degenerate_tuple = " << non_degenerate_tuple << endl;				
											
					if (non_degenerate_tuple != -1)
					{
						/* positive crossing map to negative triple points, negative crossings to positive triple points */
						k_chain[colouring_index][non_degenerate_tuple] -= sign;
									
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "cocycle_invariant:     updated k_chain[" << colouring_index << "][" << non_degenerate_tuple << "] to " << k_chain[colouring_index][non_degenerate_tuple] << endl;		
					}
				}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "cocycle_invariant:  3-chain stands at ";
	for (int i=0; i< num_chain_generators; i++)
		debug << k_chain[colouring_index][i] << ' ';
	debug << endl;
}		

				colouring_index++;
				lptr++;
			}
		}					
	}

if (debug_control::DEBUG >= debug_control::BASIC)	
{
	debug << "cocycle_invariant: k_chain, rows of :"<< endl;
	print(k_chain, debug,3,"cocycle_invariant: ");
}				
	/* boundary test for k-chain */
	if (B.numcols())
	{
		matrix<int> B_copy=B;
		matrix<int> B_P;
		echelon(B_copy,field_coefficients,true,&B_P);

if (debug_control::DEBUG >= debug_control::DETAIL)	
{
	debug << "cocycle_invariant: change of basis matrix B_P: " << endl;
	print(B_P, debug,3,"cocycle_invariant: ");
}				
					
		matrix<int> chains(num_chain_generators,B.numcols()+k_chain.numrows());

		for (size_t c=0; c< B.numcols(); c++)
		for (int i=0; i< num_chain_generators; i++)
			chains[i][c] = B[i][c];
		
		for (size_t r=0; r< k_chain.numrows(); r++)
		for (int i=0; i< num_chain_generators; i++)
			chains[i][B.numcols()+r] = k_chain[r][i];

if (debug_control::DEBUG >= debug_control::DETAIL)	
{
	debug << "cocycle_invariant: chains: " << num_chain_generators << " x  " << B.numcols()+k_chain.numrows() << endl;
	print(chains, debug,3,"cocycle_invariant: ");
}				

		matrix<int> chains_COB = B_P*chains;
				
if (debug_control::DEBUG >= debug_control::DETAIL)	
{
	debug << "cocycle_invariant: chains_COB: " << endl;
	print(chains_COB, debug,3,"cocycle_invariant: ");
}				
	}
	/* end of boundary test for k-chain */
	
	for (size_t c=0; c< colourings.size(); c++)
	{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "cocycle_invariant:  3-chain corresponding to colouring " ;
		for (int i=0; i< num_chain_generators; i++)
			debug << setw(3) << k_chain[c][i];
		debug << endl;
}

		if(!braid_control::RAW_OUTPUT && braid_control::EXTRA_OUTPUT)
		{
			output << "k-chain ";
			vector<int> chain_vector(num_chain_generators);
			
			for (int i=0; i< num_chain_generators; i++)
			{
				output << k_chain[c][i] << ' ';
				chain_vector[i] = k_chain[c][i];
			}
			print_k_chain(output, chain_vector, switch_data,(braid_control::DOUBLE_BIRACKS?3:2));
			output << endl;
		}					
						
	
		int number_of_k_tuples = 0;
		for (int i=0; i< num_chain_generators; i++)
			number_of_k_tuples += abs(k_chain[c][i]);
					
//if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
//	debug << "cocycle_invariant:  3-chain " << c << " contains " << number_of_k_tuples << " 3-tuples" << endl;
	
		list<vector<scalar> >::iterator lptr = switch_data.cocycle_scalar.begin();
		int cocycle_index = 0;
		char variable_char = 's'; //(braid_control::BIRACK_POLYNOMIAL? 's': 't');

		while (lptr != switch_data.cocycle_scalar.end())
		{

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "cocycle_invariant:  cocycle " << cocycle_index << "                          ";
		vector<scalar>& cocycle=*lptr;
		
			scalar exponent = 0; // the exponent lies in the coefficient group of the cocycle, so is a scalar
			
			for (int i=0; i< num_chain_generators; i++)
			{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << setw(3) << cocycle[i];
					exponent += k_chain[c][i]*cocycle[i];
			}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << ": exponent = " << exponent;
	
			ostringstream oss;
			if (exponent == scalar(0))
				oss << "1";
			else if (exponent == scalar(1))
				oss << variable_char;
			else
				oss << variable_char << "^" << exponent;
	
			Cpolynomial term(oss.str());
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << ", term = " << term;
	
			invariant[cocycle_index] += term;
					
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << ", invariant = " << invariant[cocycle_index] << endl;
					
			cocycle_index++;
			lptr++;
		}
	}
	
	return invariant;
}

/* peer_code_colouring_invariant evaluates invariants based on colourings of a diagram described by code_data.  It supports
	
	 BIRACK_POLYNOMIAL, based number of colourings or COCYCLE_INVARIANT coefficients of t^{writhe}
	 COCYCLE_INVARIANT of a given diagram
	 number of colourings 
	
   In the first case, if the cocycle variant is selected the function evaluates multiple invariants, one for each of the cocycles recorded
   in the switch_data supplied to the function.  If we are using the number of colourings as coeffients, we can use REFINE_RACK_POLYNOMIAL
   to describe the number of colourings as a polynomial \Sum x_i s^i whose coefficient x_i are the number of colourings of size i.
*/
void peer_code_colouring_invariant(matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, matrix<int>& Tu, matrix<int>& Td,
			braid_control::ST_pair_type pair_type, string input_string, string title, generic_switch_data& switch_data, int period)
{	
	int n = Su.numcols();  // can't use switch_data.size, since Su,Sd may be the double of the underlying switch in switch_data
	int writhe = 0;
//	int turning_number = 0;

	/* first isolate any qualifiers from the input string */
	string qualifier;
	
	string::size_type pos = input_string.find('{');
	if (pos != string::npos)
	{
		qualifier = input_string.substr(pos,string::npos);
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "peer_code_colouring_invariant: isolated qualifier = " << qualifier << endl;
		input_string = input_string.substr(0,pos);
	}

	bool plane_reflect_code_qualifier = false;
	if (qualifier.find("plane-reflect") != string::npos)
	{
		plane_reflect_code_qualifier = true;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "peer_code_colouring_invariant: plane-reflect qualifier provided" << endl;
	}

	bool reverse_code_orientation_qualifier = false;
	if (qualifier.find("reverse") != string::npos)
	{
		reverse_code_orientation_qualifier = true;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "peer_code_colouring_invariant: reverse qualifier provided" << endl;
	}

	generic_code_data code_data;
	read_peer_code(code_data, input_string);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "peer_code_colouring_invariant: code data ";
	print_code_data(debug,code_data,"peer_code_colouring_invariant:   ");
}

	if (braid_control::PLANE_REFLECT_INPUT || plane_reflect_code_qualifier)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "peer_code_colouring_invariant: replacing code data ";
	write_code_data(debug,code_data);
	debug << endl;
}
		/* plane reflect the braid, which means inverting the positive and negative labels in code_data */
		for (int i=0; i< code_data.num_crossings; i++)
		{
			if (code_data.code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::POSITIVE)
				code_data.code_table[generic_code_data::table::LABEL][i] = generic_code_data::label::NEGATIVE;
			else if (code_data.code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::NEGATIVE)
				code_data.code_table[generic_code_data::table::LABEL][i] = generic_code_data::label::POSITIVE;
		}
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "peer_code_colouring_invariant: plane-reflect code data ";
	write_code_data(debug,code_data);
	debug << endl;
}

	   	if (!braid_control::SILENT_OPERATION)
	   	{
			cout << "reflected input ";
			write_code_data(cout,code_data);
			cout << endl;
		}

		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? ";" : "");
			output << "reflected input ";
			write_code_data(output,code_data);
			output << endl;
		}
	}	
	
	if (braid_control::REVERSE_INPUT_ORIENTATION || reverse_code_orientation_qualifier)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "generic_code: replacing code data ";
	write_code_data(debug,code_data);
	debug << endl;
}
		vector<int> shift(code_data.num_components);
		for (int i=0; i< code_data.num_components; i++)
			shift[i] = -1*code_data.num_component_edges[i];
			
		renumber_peer_code(code_data,shift);
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "generic_code: reverse orientation code data ";
	write_code_data(debug,code_data);
	debug << endl;
}
	   	if (!braid_control::SILENT_OPERATION)
	   	{
			cout << "reversed input ";
			write_code_data(cout,code_data);
			cout << endl;
		}

		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? ";" : "");
			output << "reversed input ";
			write_code_data(output,code_data);
			output << endl;
		}
	}

	
	
	vector<generic_code_data> birack_poly_input;

	/* if we're doing COCYCLE_INVARIANT, make sure we've calculated the homology generators */
	if (braid_control::COCYCLE_INVARIANT)
	{						
		bool check_cocycle_conditions = false;

		if (switch_data.cocycles_calculated == false)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "switch_data does not contain cocycles" << endl;
				
			bool hold_silent_operation = braid_control::SILENT_OPERATION;
			
			if (braid_control::RAW_OUTPUT)
				braid_control::SILENT_OPERATION = true;
				
			birack_homology_generators(switch_data,(braid_control::DOUBLE_BIRACKS?3:2),true); //cohomology = true
			
			braid_control::SILENT_OPERATION = hold_silent_operation;
			
			switch_data.cocycles_calculated = true;
			check_cocycle_conditions = true;
		}
		else
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "switch_data contains cocycles" << endl;
		}
				
		/* birack_homology_generators has calculated all the cohomology generators, recording them in switch_data.cocycle_string,
		   or we have read cocycles from an input file.  If we have not yet computed, in switch_data.cocycle_scalar, a 
		   list of those generators that satisfy the cocycle condition, do so now.
		*/
		if (switch_data.cocycle_scalar.size() == 0)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "peer_code_colouring_invariant: no integer cocycles found in switch data" << endl;
			determine_cohomology_generators(switch_data);		
								
			if (check_cocycle_conditions)
				test_cohomology_generators(switch_data);
			
		}
		else
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "peer_code_colouring_invariant: switch_data contains integer cocycles" << endl;
		}
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "peer_code_colouring_invariant: switch_data:" << endl;
	print_switch_data(debug,switch_data,"peer_code_colouring_invariant:   ");
}
	}

	if (braid_control::BIRACK_POLYNOMIAL)
	{
		writhe = code_data_writhe(code_data);

		/* if we are using COCYCLE_INVARIANT coefficients, we have one invariant for each cocycle in the switch_data, 
		   so as we iterate through the number of invariant terms we constuct a matrix of coefficient polynomials, one row for each 
		   cocycle.  If we are using the refined number of colourings, we just have one invariant to calculate.
		*/
		int num_invariants = (braid_control::COCYCLE_INVARIANT? switch_data.cocycle_scalar.size(): 1);	
		int num_invariant_terms  = (braid_control::COCYCLE_INVARIANT? 2*braid_control::birack_poly_writhe_limit+1: period);	
	
//		polynomial<int> birack_poly;
		Cpolynomial birack_poly;
		vector<int> rack_poly_coefficients(num_invariant_terms);


if (debug_control::DEBUG >= debug_control::SUMMARY)
{
  	debug << "peer_code_colouring_invariant: num_invariants = " << num_invariants << endl;
  	debug << "peer_code_colouring_invariant: writhe = " << writhe << endl; //", turning number = " << turning_number << endl;
}

		if (braid_control::COCYCLE_INVARIANT)
		{
			birack_poly_input = vector<generic_code_data>(2*braid_control::birack_poly_writhe_limit+1);
			
			/* if writhe > braid_control::birack_poly_writhe_limit, add negative crossings until writhe == braid_control::birack_poly_writhe_limit */
			for (int i=braid_control::birack_poly_writhe_limit; i< writhe; i++)
				code_data = add_Reidemeister_1_loop(code_data,false); // positive_crossing = false

			/* if writhe < -braid_control::birack_poly_writhe_limit, add positive crossings until writhe == -1*braid_control::birack_poly_writhe_limit */
			for (int i=writhe; i< -1*braid_control::birack_poly_writhe_limit; i++)
				code_data = add_Reidemeister_1_loop(code_data,true); // positive_crossing = false

			writhe = code_data_writhe(code_data);
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "peer_code_colouring_invariant: after adjustment to within braid_control::birack_poly_writhe_limit, writhe = " << writhe << endl;

			generic_code_data new_code_data = code_data;			
			birack_poly_input[writhe+braid_control::birack_poly_writhe_limit] = code_data;
			
			for (int i=writhe-1; i >= -1*braid_control::birack_poly_writhe_limit; i--)
			{
				new_code_data = add_Reidemeister_1_loop(new_code_data,false); // positive_crossing = false
				birack_poly_input[i+braid_control::birack_poly_writhe_limit] = new_code_data;
			}

			new_code_data = code_data;
			for (int i=writhe+1; i <= braid_control::birack_poly_writhe_limit; i++)
			{
				new_code_data = add_Reidemeister_1_loop(new_code_data,true); // positive_crossing = false
				birack_poly_input[i+braid_control::birack_poly_writhe_limit] = new_code_data;
			}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
  	debug << "peer_code_colouring_invariant: birack_poly_input:"<< endl;
  	for (int i=0; i< num_invariant_terms; i++)
  	{
		debug << "peer_code_colouring_invariant:   ";
		write_peer_code(debug,birack_poly_input[i]);
		debug << endl;
	}
}
		}

//		matrix<Cpolynomial,scalar> coefficient_poly(num_invariants,num_invariant_terms);
		matrix<Cpolynomial,int> coefficient_poly(num_invariants,num_invariant_terms);
	
		if (braid_control::BIRACK_POLYNOMIAL && !braid_control::SILENT_OPERATION) // && braid_control::EXTRA_OUTPUT)
			cout << "term " ;

		for (int i=0; i< num_invariant_terms; i++)
		{
			if (braid_control::BIRACK_POLYNOMIAL && !braid_control::SILENT_OPERATION) // && braid_control::EXTRA_OUTPUT)
				cout << i << ' ' << flush;
			
			if (braid_control::COCYCLE_INVARIANT)
			{
				code_data = birack_poly_input[i];

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
  	debug << "peer_code_colouring_invariant: term " << i << " code_data: ";
  	 write_peer_code(debug,code_data);
  	 debug << endl;
}  	 
			}
			else if (i==0)
			{
if (braid_control::BIRACK_POLYNOMIAL && debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "peer_code_colouring_invariant: term " << i << endl;
			}
			else
			{	
				code_data = add_Reidemeister_1_loop(code_data,true); // positive_crossing = true

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
  	debug << "peer_code_colouring_invariant: term " << i << " adding crossing " << code_data.num_crossings-1 << ", new code: ";
  	 write_peer_code(debug,code_data);
  	 debug << endl;
}  	 
				
			}
			
			list<vector<int> > colourings = calculate_colourings(code_data,Su,Sd,invSu,invSd,Tu,Td);
			
			int num_colourings = colourings.size();
					
			/* we need exponent to be positive but the % operator can give negative results: n = -1 p = 5 n%p = -1 (n%p+p)%p = 4 */
			int exponent = (braid_control::COCYCLE_INVARIANT? i:((writhe+i)%period+period)%period); 
			rack_poly_coefficients[exponent] = num_colourings;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "peer_code_colouring_invariant:  invariant term " << i << ", exponent = " << exponent << ", num_colourings = " << num_colourings << endl;
										
			if (braid_control::REFINE_RACK_POLYNOMIAL)
			{
				vector<int> image_size_count(n);
				list<vector<int> >::iterator lptr = colourings.begin();
				
				while (lptr != colourings.end())
				{
					/* identify the labels used in the colouring lptr*/
					vector<int> labels(n);
					for (size_t i=0; i< (*lptr).size(); i++)
						labels[(*lptr)[i]] = 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "peer_code_colouring_invariant:  cumulative label record ";
	for (int j=0; j< n; j++)
		debug << labels[j] << ' ';	
	debug << endl;
}							
					int num_labels=0;
					for (int i=0; i< n; i++)
						num_labels += labels[i];
						
					vector<int> label_values(num_labels);
					int index=0;
					for (int i=0; i< n; i++)
					{
						if (labels[i] == 1)
							label_values[index++] = i;
					}
		
					vector<int> parent_birack = smallest_parent_birack(Su,Sd,label_values);
		
					int image_size=parent_birack.size();
				

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "peer_code_colouring_invariant:  image_size = " << image_size;
	debug << " labels: ";
	for (int j=0; j< n; j++)
		debug << labels[j] << ' ';	
	debug << endl;

	debug << "peer_code_colouring_invariant:  parent_birack: ";
	for (int j=0; j< image_size; j++)
		debug << parent_birack[j] << ' ';	
	
	for (int j=image_size; j< n; j++)
		debug << "  ";
		
	debug << "numbering from 1: ";
	for (int j=0; j< image_size; j++)
		debug << parent_birack[j]+1 << ' ';	
	debug << endl;
}
					image_size_count[image_size-1] ++;
					lptr++;
				}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
  	debug << "peer_code_colouring_invariant:   number of colourings =  " << num_colourings << endl;
  	debug << "peer_code_colouring_invariant:   fixed point image size count =  ";
  	for (int j=0; j< n; j++)
		debug << image_size_count[j] << ' ';
	debug << endl;
}
												
//				polynomial<int> coeff_poly;
				Cpolynomial coeff_poly;
				
				for (int j=0; j<n; j++)
				{
					if (image_size_count[j] !=0)
					{
						ostringstream oss_c;
						oss_c << image_size_count[j];
						if (j > 0)
							oss_c << "s^" << j;
						
//						polynomial<int> c_term(oss_c.str());
						Cpolynomial c_term(oss_c.str());
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
  	debug << "peer_code_colouring_invariant:   coefficient term: " << c_term << endl;
		  	
						coeff_poly += c_term;
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "peer_code_colouring_invariant:   coeff_poly updated to : " << coeff_poly << endl;
		
					}
				}
				
				coefficient_poly[0][exponent] = coeff_poly;						
			}
			else if (braid_control::COCYCLE_INVARIANT)
			{
//				vector<polynomial<int> > invariant = cocycle_invariant(code_data,switch_data,colourings);		
				vector<Cpolynomial> invariant = cocycle_invariant(code_data,switch_data,colourings);		
						
				/* set coefficient_poly[exponent] to the first cocycle invariant for the colourings of this writhe */

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
  	debug << "peer_code_colouring_invariant:   cocycle invariants: ";
  	for (int i=0; i< num_invariants; i++)
		debug << invariant[i] << ' '; 
  	debug << endl;
}
				for (int i=0; i< num_invariants; i++)
					coefficient_poly[i][exponent] = invariant[i];												
			}							
			else
			{
				ostringstream oss_p;			
				oss_p << num_colourings << "t^" << exponent;
	
//				polynomial<int> p_term(oss_p.str());
				Cpolynomial p_term(oss_p.str());
			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
  	debug << "peer_code_colouring_invariant:   writhe = : " << writhe+i << ", polynomial exponent " << exponent << endl;
  	debug << "peer_code_colouring_invariant:   polynomial term: " << p_term << endl;
}
				
				birack_poly += p_term;
						
if (debug_control::DEBUG >= debug_control::SUMMARY)
  	debug << "peer_code_colouring_invariant:   birack polynomial updated to : " << birack_poly << endl;
			}
		}


		if (!braid_control::SILENT_OPERATION)
		{
			if (braid_control::REFINE_RACK_POLYNOMIAL || braid_control::COCYCLE_INVARIANT)
			{
				
				if (braid_control::BIRACK_POLYNOMIAL)
					cout << endl;
				
				for (int i=0; i< num_invariants; i++)
				{
					bool first_term = true;
					
					if(braid_control::COCYCLE_INVARIANT)
						cout << "Cocycle " << i+1 << " birack polynomial = ";
					else
						cout << "Birack polynomial = ";
						
					for (int j=num_invariant_terms-1; j >=0; j--)
					{
						if (coefficient_poly[i][j].non_zero())
						{
							int exponent = (braid_control::COCYCLE_INVARIANT?j-braid_control::birack_poly_writhe_limit:j);

							bool parentheses;
							if (coefficient_poly[i][j].pt.size() >1  && exponent != 0)
								parentheses = true;
							else
								parentheses = false;

							if (!first_term)
							{
								if  (parentheses || !coefficient_poly[i][j].is_negative())
									cout << '+';
								else if (coefficient_poly[i][j].is_minus_one())
									cout << '-';
							}
							
							if (parentheses) 						cout << "(";
							if (!coefficient_poly[i][j].is_one()|| exponent == 0)	cout << coefficient_poly[i][j];
							if (parentheses)						cout << ")";
							
							if (exponent != 0)
							{
								cout << 't';
								if (exponent !=1)
								{
									cout << '^';
									cout << exponent;
								}
							}							
							first_term = false;
						}
					}
					cout << endl;
				}
			}
			else
			{
				cout << "Birack polynomial = " << birack_poly << endl;
			}
		}

		
		if (braid_control::REFINE_RACK_POLYNOMIAL || braid_control::COCYCLE_INVARIANT)
		{			
			for (int i=0; i< num_invariants; i++)
			{
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					if(braid_control::COCYCLE_INVARIANT)
						output << "Cocycle " << i+1 << " birack polynomial = ";
					else
						output << "Birack polynomial = ";
				}

				bool first_term=true;
				for (int j=num_invariant_terms-1; j >=0; j--)
				{
					sanitize(coefficient_poly[i][j]);
					if (coefficient_poly[i][j].non_zero())
					{
						int exponent = (braid_control::COCYCLE_INVARIANT?j-braid_control::birack_poly_writhe_limit:j);

						bool parentheses;
						if (coefficient_poly[i][j].pt.size() >1  && exponent != 0)
							parentheses = true;
						else
							parentheses = false;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
  	debug << "peer_code_colouring_invariant:   parentheses required = " << parentheses << endl;

						if (!first_term)
						{
							if  (parentheses || !coefficient_poly[i][j].is_negative())
								output << '+';
							else if (coefficient_poly[i][j].is_minus_one())
								output << '-';
						}
						
						if (parentheses) 						output << "(";
						if (!coefficient_poly[i][j].is_one()|| exponent == 0)	output << coefficient_poly[i][j];
						if (parentheses)						output << ")";	

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "peer_code_colouring_invariant: coefficient_poly[" << i << "][" << j << "] = " << coefficient_poly[i][j] << endl; // << " is_one = " << unit_coefficient_poly << endl;
	
	if (debug_control::DEBUG >= debug_control::DETAIL)
		dump(debug,coefficient_poly[i][j],"peer_code_colouring_invariant:");
}						
						if (exponent != 0)
						{
							output << 't';
							if (exponent != 1)
							{
								output << '^';
								output << exponent;
							}
						}							
						first_term = false;
					}
				}
				output << endl;
			}
		}
		else
		{
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Birack polynomial = ";
			}
			
			output << birack_poly << endl;
		}
		
    	if (braid_control::EXTRA_OUTPUT)
    	{
			if (!braid_control::REFINE_RACK_POLYNOMIAL && !braid_control::COCYCLE_INVARIANT)
			{
				cout << "Birack polynomial coefficients t^0 to t^" << period-1 << ": ";

				for (int i=0; i < period; i++)
					cout << rack_poly_coefficients[i] << ' ';
				cout << endl;

		    	if (!braid_control::RAW_OUTPUT)
				{
					output << "Birack polynomial " << birack_poly << endl;
					output << "Birack polynomial coefficients t^0 to t^" << period-1 << ": ";					
					for (int i=0; i < period; i++)
						output << rack_poly_coefficients[i] << ' ';
					output << endl;			
				}
			}
		}	

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	if (!braid_control::REFINE_RACK_POLYNOMIAL)
		debug << "peer_code_colouring_invariant: final birack_poly = " << birack_poly << endl;
	debug <<"birack polynomial coefficients = ";
    
    if (braid_control::REFINE_RACK_POLYNOMIAL)
    {
		for (int i=0; i < period; i++)
			debug << coefficient_poly[0][i] << ' ';
	}
	else
	{
		for (int i=0; i < period; i++)
			debug << rack_poly_coefficients[i] << ' ';
	}
	debug << endl;
}				
	}
	else if (braid_control::COCYCLE_INVARIANT)
	{			
		list<vector<int> > colourings = calculate_colourings(code_data,Su,Sd,invSu,invSd,Tu,Td);

		int num_cocycles = switch_data.cocycle_scalar.size();

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "peer_code_colouring_invariant: calculating cocycle_invariant with num_cocycles = " << num_cocycles << endl;
			
		if (num_cocycles != 0)
		{
//			vector<polynomial<int> > invariant = cocycle_invariant(code_data,switch_data,colourings);		
			vector<Cpolynomial> invariant = cocycle_invariant(code_data,switch_data,colourings);		

			list<string>::iterator lptr=switch_data.cocycle_string.begin();
			int cocycle_index=0;
			
			while (lptr != switch_data.cocycle_string.end())
			{
				
				if (!braid_control::SILENT_OPERATION)
				{
					cout << "\ncocycle " << *lptr << endl;
					cout << "invariant " << invariant[cocycle_index] << endl;
				}
		    
				if (!braid_control::RAW_OUTPUT)
		    	{
					
if (true && invariant[cocycle_index].nv >0) // added only to identify non trivial output easily for inclusion in a paper 
{					
	output << (braid_control::OUTPUT_AS_INPUT? "\n;XX " : "\nXX ");					
	write_peer_code(output,code_data);

	int size = switch_data.size;
	matrix<int> twitch_u(size,size,-1);
	matrix<int> twitch_d(size,size,-1);
	
	matrix<int> inv_D(size,size,-1);
	
	for (int i=0; i< size; i++)
	for (int j=0; j< size; j++)
		inv_D[i][switch_data.Sd[i][j]] = j; // inverts the map D_i given by the row Sd[i][]

	for (int a=0; a< size; a++)
	for (int b=0; b< size; b++)
	{
		int x = inv_D[a][b];  // D_a^{-1}(b)


		twitch_d[a][b] = x;
		twitch_u[b][a] = switch_data.Su[x][a];
	}
	bool saved_bool = matrix_control::SINGLE_LINE_OUTPUT;
	matrix_control::SINGLE_LINE_OUTPUT = true;
	output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");					
	output << "XX Su = ";
	display_fixed_point_switch(twitch_u, output, false);
//	display_fixed_point_switch(switch_data.Su, output, false);
	output << " Sd = ";
	display_fixed_point_switch(twitch_d, output, false);
//	display_fixed_point_switch(switch_data.Sd, output, false);
	matrix_control::SINGLE_LINE_OUTPUT = saved_bool;

	output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
	output << "XX cocycle " << cocycle_index+1 << ": " << *lptr << endl;
	output << "XX invariant ";
}
else
{				
		    		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "cocycle " << cocycle_index+1 << ": " << *lptr << endl;
					output << "invariant ";
}				
		    	}
		    	output << invariant[cocycle_index] << endl;
		    	
				lptr++;
				cocycle_index++;
			}			
		}
	}
	else
	{
		list<vector<int> > colourings = calculate_colourings(code_data,Su,Sd,invSu,invSd,Tu,Td);

		int fixed_points = colourings.size();
    	
		if (!braid_control::SILENT_OPERATION)
			cout <<"\nNumber of fixed-points = " << fixed_points << endl;
    
		if (!braid_control::RAW_OUTPUT)
    	{
    		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
    		output << "Number of fixed points = ";
    	}
    	output << fixed_points << endl;
		    
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "peer_code_colouring_invariant: total number of fixed points = " << fixed_points << endl;
	}	
}
