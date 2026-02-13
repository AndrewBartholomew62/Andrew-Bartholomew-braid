/**************************************************************************

void bracket_polynomial(generic_code_data& code_data, int variant)
int amalgamate_zig_zag_counts(int a, int b)

braid-code-1-9-20.tar contains code prior to the return of vectors from the reidemeister removal functions.

unicursal_component_zig_zag_count contains the correct counts for a relaxed parity arrow polynomial 

 **************************************************************************/
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <valarray>
#include <list>
#include <map>
#include <iomanip>
#include <ctype.h>
#include <stdio.h>
#include <algorithm>

using namespace std;

extern ofstream     debug;
extern ofstream     output;
extern ifstream     input;

#include <util.h>
#include <scalar.h>
#include <quaternion-scalar.h>
#include <polynomial.h>
#include <matrix.h>
#include <braid.h>
#include <braid-util.h>
#include <generic-code.h>
#include <gauss-orientation.h>
#include <reidemeister.h>

class bracket_variable
{
public:
	enum variable_type {NONE,STRING,GRAPH};
	
	int _type;
	generic_code_data code_data;
//	gauss_orientation_data gauss_data;
	string _string; // unoriented left preferred Gauss code

	bracket_variable ();
	bracket_variable (generic_code_data& c);
	bracket_variable (string s);

	bool operator == (const bracket_variable& bv) const {return (_string == bv._string);}		
	bool operator != (const bracket_variable& bv) const {return (_string != bv._string);}		
	friend ostream& operator << (ostream& os, const bracket_variable&);
};

bracket_variable::bracket_variable()
{
	_type = bracket_variable::variable_type::NONE;
	code_data.num_crossings = 0;
}

/* For parity arrow, the generic_code_data may represent a knotoid with a non-zero head_zig_zag_count.  The construction 
   of gauss_orientation_data doesn't consider this class member, nor does it consider the count recorded in c.zig_zag_count 
   corresponding to the head semi_arc, since this arc terminates at a shortcut crossings.  Thus, the gauss_orientation_data 
   includes only zig-zag counts for odd crossings.  
*/
bracket_variable::bracket_variable(generic_code_data& c): code_data(c)
{
	_type = bracket_variable::variable_type::GRAPH;
	gauss_orientation_data g(code_data);

	ostringstream oss;
	
	/* if the code_data is a knotoid then we cannot use the left_preferred Gauss data, since
	   we have to start at the leg and proceed towards the head.
	*/
	gauss_orientation_data gauss_data = left_preferred(g,true,c.immersion_character); //unoriented = true	
	write_gauss_data(gauss_data,oss,true,c.immersion_character,c.head_zig_zag_count); // zig_zags = true
	_string = oss.str();
	
	/* write the knotoid or long knot head semi-arc count to the end of _string, before the ')' character
	if (c.immersion != generic_code_data::character::CLOSED)	
	{
		size_t pos = _string.find(')');
		if (pos != string::npos)
			_string.insert(pos," "+to_string(abs(c.head_zig_zag_count)));
	}
	*/
}

bracket_variable::bracket_variable(string s)
{
	_type = bracket_variable::variable_type::STRING;
	_string = s;
}


ostream& operator << (ostream& os, const bracket_variable& bv)
{
	if (bv._type == bracket_variable::variable_type::GRAPH)
	{
		os << bv._string;
/*		if (bv.code_data.zig_zag_count.numrows()) // i.e. if there's anything there
		{
			bool save = matrix_control::SINGLE_LINE_OUTPUT;
			matrix_control::SINGLE_LINE_OUTPUT = true;
			os << " (" << bv.code_data.zig_zag_count << ")";
			matrix_control::SINGLE_LINE_OUTPUT = save;
		}
*/		
		if (polynomial_control::WRITE_PARITY_PEER_CODES)
		{
			os << " ";
			write_peer_code(os,bv.code_data, true); // zig_zags = true
		}
	}
	else
	{
		os << bv._string;
	}
	return os;
}

void print_bracket_variable(const bracket_variable bv, ostream& os, string prefix)
{
	os << prefix << "_type: ";
	switch (bv._type)
	{
		case bracket_variable::variable_type::NONE: os << "NONE"; break;
		case bracket_variable::variable_type::STRING: os << "STRING"; break;
		case bracket_variable::variable_type::GRAPH: os << "GRAPH"; break;
	}
	os << endl;
	os << prefix << "code_data: ";
	write_peer_code(os,bv.code_data);
	os << endl;
	os << prefix << "zig-zag_count: " << endl;
	print(bv.code_data.zig_zag_count,os,4,prefix);
	os << prefix << "_string: " << bv._string << endl;
}


vector<int> gauss_parity(generic_code_data& code_data);
int remove_virtual_Reidemeister_I(generic_code_data& code_data,vector<int>& component_flags);
matrix<int> create_incidence_matrix (gauss_orientation_data& gauss_data);
int remove_virtual_components(generic_code_data& code_data,vector<int>& component_flags);
generic_code_data partition_peer_code(generic_code_data& code_data, vector<int>& component_flags);

/* The function bracket_polynomial evaluates various forms of bracket polynomial depending on the 
   variant parameter supplied.  The function is capable of evaluating
   
    - Kauffman's bracket polynomial
    - The Jones polynomial
    - Turaev's extended bracket polynomial for knotoids
	- The arrow polynomial
	- The parity bracket polynomial
	- The parity arrow polynomial
     
   All these bracket polynomials are calculated using the state summation model, the definition of which
   uses the smoothing convention that a counter-clockwise rotation of the over-crossing arc sweeps two 
   regions labelled A and a smoothing that connects the A regions is known as an A-smoothing, a smoothing 
   that does not connect these regions is known as a B-smoothing.  Thus for a positive and negative 
   crossing we have:

     \ /          \_/          \   /
   A  \  A  =      _            | |    
     / \          / \          /   \
               A-smoothing    B-smoothing

      A
     \ /         \   /          \_/
      /     =     | |            _ 
     / \         /   \          / \
      A        A-smoothing    B-smoothing


   We label a crossing that has been A-smoothed with A and a crossing that has been B-smoothed with an A^-1.  
   These labellings are reflected in the calculation of sigma for each bracket state.
 
   In a labelled code the labels indicate whether the naming arc is part of the over 
   arc or under arc.  Thus, a labelled code describes an unoriented link, notwithstanding the 
   fact that the enumeration of the semi-arcs of the code can be used to determine an orientation, 
   if required.  Moreover, labelled codes allow us to determine the aspect of the over-arc at each 
   crossing from the type and the label associated with each crossing.  IF (emphasis intended) we 
   endowed the immersion with the orientation induced by the semi-arc numbering then we may use
   the resultant positive and negative crossings to calculate bracket polynomials.
   
   Given this orientation it is convenient to describe an A-smoothing of a "positive" crossing and a 
   B-smoothing for a "negative" crossing as a Seifert-smoothing, since the smoothing connects 
   semi-arcs 2i <-> 2j and 2j-1 <-> 2i+1.  We then talk of a smoothing connecting 2i <-> 2j-1 and 
   2j <-> 2i+1 as a non-Seifert smoothing.
   
   With this terminology we can determine the Kauffman bracket polynomial crossing labels by 
   considering the orientation induced by the semi-arc numbering as follows:
   
   crossing sign  Smoothing         Kauffman label
    "positive"    Seifert (A)              A
    "positive"    non-Seifert (B)          A^-1
    "negative"    non-Seifert (A)          A
    "negative"    Seifert (B)              A^-1
   
   The Turaev extended bracket polynomial for knotoids requires us to consider a shortcut established
   between the head and the leg of the knotoid.  

   Note: for the Turaev extended bracket polynomial we require that the leg of the knotoid lies in 
   the semi-arc of the immersion numbered with zero in the code_data.  This is so we can identify
   the segment component of a smoothed diagram and in particular the semi-arc containing the leg.
   Allowing arbitary numbering would require the identification of the leg as well as the head within 
   the peer code, since one could not guarantee that the segment component numbering could start with 
   the leg and preserve the requirement for odd and even terminating edges at each crossing.

   The bracket polynomial does not support code data specifying flat crossings.

   The Kauffman bracket, Jones polynomial knotoid bracket and arrow polynomial are all calculated for both knots 
   & links and knotoids & multi-knotoids.  The parity bracket polynomial implementation supports only knots and
   knotoids, not links or multi-knotoids.
   
   Gauss codes are not accepted for the Turaev extended bracket polynomial, or the parity bracket polynomial, 
   since these require us to track shortcut and virtual crossings.  

*/
void bracket_polynomial(generic_code_data& code_data, int variant)
{
	int head = code_data.head;
	vector<int>& first_edge_on_component = code_data.first_edge_on_component;
	vector<int>& num_component_edges = code_data.num_component_edges;
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial: provided with code data: ";
	write_code_data(debug, code_data);
	debug << "\nbracket_polynomial: head = " << head << endl;
	debug << "bracket_polynomial: variant requested is ";
	
	if (variant == KAUFFMAN_VARIANT)
		debug << "Kauffman's bracket polynomial" << endl;
	else if (variant == JONES_VARIANT)
		debug << "Jones Polynomial" << endl;
	else if (variant == TURAEV_VARIANT)
		debug << "Turaev's extended bracket polynomial for knotoids" << endl;
	else if (variant == ARROW_VARIANT)
		debug << "arrow polynomial" << endl;
	else if (variant == PARITY_VARIANT)
		debug << "parity bracket polynomial" << endl;
	else if (variant == PARITY_ARROW_VARIANT)
		debug << "parity arrow polynomial" << endl;
	else 
		debug << "Unknown!" << endl;				
		
	if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
		debug << "bracket_polynomial: pure knotoid detected" << endl;

	if (code_data.immersion_character == generic_code_data::character::LONG_KNOT)
		debug << "bracket_polynomial: long knot detected" << endl;
		
	if (code_data.immersion_character == generic_code_data::character::KNOTOID)
		debug << "bracket_polynomial: knotoid detected" << endl;
}
		
	matrix<int>& code_table = code_data.code_table;				
	vector<int>& term_crossing = code_data.term_crossing;
	vector<int>& orig_crossing = code_data.orig_crossing;

	
	/* check we've not been given a Gauss code when variant == PARITY_VARIANT || PARITY_ARROW_VARIANT || TURAEV_VARIANT */
	if ((variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT || variant == TURAEV_VARIANT) && code_data.type != generic_code_data::code_type::peer_code)
	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nbracket polynomial requires a peer code as input for this variant, skipping";

		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << "\n\nbracket polynomial requires a peer code as input for this variant, skipping";
			if (braid_control::OUTPUT_AS_INPUT)
				output << '\n';
		}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial: bracket_polynomial variant provided with a Gauss code rather than a peer code, doing nothing" << endl;

		return;
	}
	
	/* check we've not been given multiple components when variant == PARITY_VARIANT or variant == PARITY_ARROW_VARIANT*/
	if ((variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT) && code_data.num_components > 1)
	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nbracket polynomial does not support multiple components for this variant, skipping";

		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << "\n\nbracket polynomial does support not multiple components for this variant, skipping";
			if (braid_control::OUTPUT_AS_INPUT)
				output << '\n';
		}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial: bracket_polynomial provided with multiple components, doing nothing" << endl;

		return;
	}
	
	int num_crossings = code_data.num_crossings;

	int num_edges = 2*num_crossings;
	
	/* arrow polynomial variables */
	int arrow_state_zig_zag_count = 0;
//	int arrow_state_first_non_seifert_turn;
	int arrow_state_last_non_seifert_turn = generic_code_data::VOID;
	int non_seifert_turn = generic_code_data::VOID;
			
	/* set up a vector of flags to identify the shortcut crossings, the vector is initialized
	   to zeros which is the correct indication for links, which have no shortcut crossings.
	   In the case of knotoids the flags will be set by valid_knotoid_input.
	*/
	int num_non_shortcut_crossings = num_crossings;  //all crossings for links
	bool pure_knotoid_code_data = false;
	
	if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && head != -1)
	{
		if (valid_knotoid_input(code_data))
		{
			pure_knotoid_code_data = true;
			
			vector<int>& shortcut_crossing = code_data.shortcut_crossing;
			
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial: shortcut_crossings: ";
	for (int i=0; i< num_crossings; i++)
		debug << shortcut_crossing[i] << ' ';
	debug << endl;
}

			for (int i=0; i< num_crossings; i++)
				num_non_shortcut_crossings -= abs(shortcut_crossing[i]);

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial: num_non_shortcut_crossings = " << num_non_shortcut_crossings << endl;
	
			if (num_non_shortcut_crossings == 0)
			{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "bracket_polynomial: no knotoid crossings, normalized bracket polynomial set to 1";

				if (!braid_control::SILENT_OPERATION)
					cout << "\n\nno knotoid crossings, normalized bracket polynomial = 1";
		
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "normalized bracket polynomial = ";
					if (braid_control::OUTPUT_AS_INPUT)
						output << '\n';
				}
				output << 1 << endl;
		
				return;
			}
		}
		else
		{

			if (!braid_control::SILENT_OPERATION)
				cout << "\nbracket_polynomial: invalid knotoid input, skipping" << endl;
	
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "\nbracket_polynomial: invalid knotoid input, skipping";
				if (braid_control::OUTPUT_AS_INPUT)
					output << '\n';
			}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial: invalid knotoid input" << endl;

			return;
		}
	}
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial: pure_knotoid_code_data = " << pure_knotoid_code_data << endl;
	
	/* parity bracket variables, crossing_parity records the ODD or EVEN parity of non-virtual, 
	   non-shortcut crossings, we have to do this after checking for knotoid data, so the shortcut
	   crossings are set ready for gauss_parity.
	   
	   Count the number of odd crossings so we can track zig-zags per-odd crossing for the parity arrow
	   polynomial (see parity_arrow_gauss_zig_zags below).
	*/
	vector<int> crossing_parity;
	int num_odd_crossings = 0;
	
	if (variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT)
	{
		crossing_parity = gauss_parity(code_data);

		for (int i=0; i< num_crossings; i++)
		{
			if (crossing_parity[i] == gauss_orientation_data::parity::ODD)	
				num_odd_crossings++;
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "bracket_polynomial: crossing_parity: ";
	for (int i=0; i< num_crossings; i++)
	{
		if (crossing_parity[i] == gauss_orientation_data::parity::ODD)		
			debug << "O ";
		else if (crossing_parity[i] == gauss_orientation_data::parity::EVEN)		
			debug << "E ";
		else
			debug << "N ";
	}
	debug << endl;
	debug << "bracket_polynomial: number of odd crossings = " << num_odd_crossings << endl;
}		
	}
	
	
	/* Note the sign of each non-shortcut crossing (positive or negative crossing, or virtual crossing)
	   as determined by the orientation implied by the code data.  
	   
	   For knots and links shortcut_crossing[i] is zero for all i so we consider all the link's crossings.

	   The sign is used to calculate the writhe and also to determine the aspect of
	   each crossing during the calculation of the bracket polynomial.
	   
		NOTE: we use the generic_braid_data::crossing_type rather than generic_code_data::label for sign, since it
		relates to the standard definition of positive and negative crossings, rather than peer code labels.
	*/
	vector<int> sign(num_crossings,generic_braid_data::crossing_type::NONE);
	vector<int>& shortcut_crossing = code_data.shortcut_crossing;
	
	for (int i=0; i< num_crossings; i++)
	{		
		if (pure_knotoid_code_data && shortcut_crossing[i])
			continue;
		else
		{
			if (   (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1 && code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
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
		    else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::FLAT)
		    {
				sign[i] = generic_braid_data::crossing_type::FLAT;
			}
		    else
			{
				/* virtual crossing, we're skipping shortcut crossings, which retain a sign of generic_braid_data::crossing_type::NONE */
				sign[i] = generic_braid_data::crossing_type::VIRTUAL;
		    }			    
		}
	}
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial: non-shortcut crossing sign (positive = 1, negative = -1, virtual = 0, flat = 2, none = 9): ";
	for (int i=0; i< num_crossings; i++)
		debug << sign[i] << ' ';
	debug << endl;
}

	/* Evaluate the normalizing factor (-A^3)^{-writhe} = (-1)^{writhe}A^{-3writhe} 
	   if were evaluating the parity bracket polynomial and braid_control::EVEN_WRITHE is true, we conly consider 
	   even crossing, otherwise all non-virtual and non-shortcut crossings are included in the writhe calculation
	*/
	
	int writhe = 0;
	for (int i=0; i< num_crossings; i++)
	{
		if ((variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT) 
		     && braid_control::EVEN_WRITHE && crossing_parity[i] != gauss_orientation_data::parity::EVEN)
			continue;
		else
		{
			if (sign[i] == generic_braid_data::crossing_type::POSITIVE)
				writhe ++;
			else if (sign[i] == generic_braid_data::crossing_type::NEGATIVE)
				writhe --;
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "bracket_polynomial: " << ((variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT) && braid_control::EVEN_WRITHE? "even writhe = ": "(full) writhe = ") << writhe << endl;
	
	ostringstream oss;
	
	if (writhe)
	{
		if (writhe%2)
			oss << "-";

		writhe *= -3;
		oss << "A^" << writhe;
	}
	else
		oss << "1";

    /* in the case of Turaev's knotoid polynomial, calculate the algebraic 
	   number of intersections of the knotoid and the shortcut and complete the normalizing factor   
	*/
	int algebraic_crossing_number=0;
	
    if (variant == TURAEV_VARIANT)
	{
		for (int i=0; i<num_crossings; i++)
			algebraic_crossing_number += shortcut_crossing[i];

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "bracket_polynomial: algebraic number of intersections of knotoid and shortcut = " 
	      << algebraic_crossing_number << endl;
}

		oss << "u^" << -1 * algebraic_crossing_number;
	}

	polynomial<int,bracket_variable> normalizing_factor(oss.str());
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "bracket_polynomial: normalizing factor = " << normalizing_factor << endl;
		
	int head_semi_arc;
	
	/* If we have a knotoid, determine the head_semi_arc */
	if (pure_knotoid_code_data)
	{
		if (code_table[generic_code_data::table::LABEL][head] == generic_code_data::POSITIVE)
			head_semi_arc = code_table[generic_code_data::table::OPEER][head];
		else if (code_table[generic_code_data::table::LABEL][head] == generic_code_data::NEGATIVE)
			head_semi_arc = 2*head;
	}
	else
	{
		head_semi_arc = num_edges;  //head_semi_arc is used for loop control below so needs to be set correctly in all cases.
	}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial: head_semi_arc = " << head_semi_arc << endl;

	/* determine the number of classical non-shortcut crossings that we want to smooth and a mapping between the 
	   non-shortcut-crossings and those crossings.  The mapping is required for the proper calculation of sigma 
	   when doing parity bracket but we only calculate crossing_parity when variant == PARITY_VARIANT or 
	   variant == PARITY_ARROW_VARIANT
	*/
	int num_state_smoothed_crossings = num_crossings;
	
	/* note the number of classical non-shortcut crossings for the case where we're doing parity
	   bracket and have only odd crossings to consider
	*/
	int num_classical_non_shortcut_crossings = num_crossings;
	
	for (int i=0; i< num_crossings; i++)
	{
		if (sign[i] != generic_braid_data::crossing_type::POSITIVE && sign[i] != generic_braid_data::crossing_type::NEGATIVE)
		{
			num_state_smoothed_crossings--;
			num_classical_non_shortcut_crossings--;
		}
		else if (variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT)
		{
			/* for the parity bracket polynomial and parity arrow polynomial we smooth only the even crossings, so adjust 
			   num_state_smoothed_crossings accordingly in this case
			*/
			// enum parity {NONE = 0, ODD = 1, EVEN = 2}
			if (crossing_parity[i] == gauss_orientation_data::parity::ODD)
				num_state_smoothed_crossings--;
		}
	}
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial: num_classical_non_shortcut_crossings = " << num_classical_non_shortcut_crossings << endl;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "bracket_polynomial: num_state_smoothed_crossings = " << num_state_smoothed_crossings << endl;

	
	/* set up the bracket_polynomial (which initializes to zero) and the state vector for classical 
	   non-shortcut crossings.  The state vector takes values -1 or 1.
	   
	   For variants other than PARITY_VARIANT or PARITY_ARROW_VARIANT, if there are no crossings to smooth, the bracket polynomial is 1.
	   In the case of the parity bracket, we may have no smoothed (even) crossing but do have graphical nodes
	   (odd crossings) to consider.
	*/
	polynomial<int,bracket_variable> bracket_poly;

	if (num_state_smoothed_crossings > 0 || ((variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT )&& num_classical_non_shortcut_crossings > 0))
	{
	
		bool finished=false;
	
		vector<int> state(num_state_smoothed_crossings);
		for (int i=0; i< num_state_smoothed_crossings; i++)
			state[i] = -1;
	
		int state_count=0;  // used for comfort dots
		do
		{
			/* determine the contribution to the bracket polynomial from the current state.  We need to count the 
			   number of components for this state and evaluate sigma, the "state sum".
			   
			   We count the components by tracing around the diagram following the smoothing determined by the state.

			   In the case of the Turaev variant for knotoids, we must also evaluate the algebraic crossing number 
			   of the segment component with the shortcut.			   
			*/
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "bracket_polynomial: state ";
	for (int i=0; i< num_state_smoothed_crossings; i++)
		debug << state[i] << ' ';
	debug << endl;
}
	
			/* map the state onto a record of the all crossings in order to accommodate 
			   virtual crossings and shortcut crossings, as these are never smoothed.
			   
			   We shall perform Seifert smoothing, where 2i <-> 2j and 2j-1 <-> 2i+1, on classical 
			   non-shortcut crossings whose state is 1 and non-Seifert smoothing, where 2i <-> 2j-1 
			   and 2j <-> 2i+1 on those whose state is -1.
			   
			   This is not the usual interpretation of the state, which is normally taken to identify 
			   A-smoothing or B-smoothings.  However, our interpretation produces an equivalent result
			   and we allow for our different interpretation when evaluating the state's contribution 
			   to the bracket polynomial.
			   
			   We maintain a boolean to identify whether we are following the underlying immersion 
			   orientation or not.		   
			*/
			
			vector<int> crossing_record(num_crossings);
			int place=0;
			for (int i=0; i< num_crossings; i++)
			{
				if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::VIRTUAL)
				{
					crossing_record[i] = VIRTUAL_CROSSING;
				}
				else if (pure_knotoid_code_data && shortcut_crossing[i])
				{
					crossing_record[i] = SHORTCUT_CROSSING;
				}
				else if ((variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT) && crossing_parity[i] == gauss_orientation_data::parity::ODD)
				{
					crossing_record[i] = ODD_CROSSING;
				}
				else // classical non-shortcut crossing
				{
					if (state[place] == -1)
						crossing_record[i] = NON_SEIFERT_SMOOTHED;
					else
						crossing_record[i] = SEIFERT_SMOOTHED;
					
					place++;
				}
			}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "bracket_polynomial: crossing_record (SEIFERT_SMOOTHED = 1, NON_SEIFERT_SMOOTHED = -1, VIRTUAL = 2, SHORTCUT_CROSSING = 3, ODD = 4)" << endl;
	debug << "bracket_polynomial: crossing_record ";
	for (int i=0; i< num_crossings; i++)
		debug << crossing_record[i] << ' ';
	debug << endl;
}
			/* component_flag records the component to which each edge belongs.  It's a flag to indicate whether 
			   the edge has been considered when tracing the state components, but is recorded as the component, 
			   counted from one, to aid debugging. 
			*/
			vector<int> component_flag(num_edges);	
			
			
			if (pure_knotoid_code_data)
			{
				/* set the component_flag for the shortcut edges to -1, so we don't attempt to start another component
				   on a shortcut edge, when tracing smoothed diagram component.  Note that valid_braid_input has already 
				   checked that the knotoid leg and conseqently the head lie on the first component of the knotoid diagram.
				*/
				int semi_arc=-1;  // the first immersion semi-arc of the shortcut
				
				if (code_table[generic_code_data::table::LABEL][head] == generic_code_data::POSITIVE)
					semi_arc = code_table[generic_code_data::table::OPEER][head];
				else if (code_table[generic_code_data::table::LABEL][head] == generic_code_data::NEGATIVE)
					semi_arc = 2*head;
				else
				{
					/* shouldn't get here as valid_knotoid_input checks the first shortcut crossing cannot be virtual */
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial: head " << 0 << " indicates first shortcut crosing is virtual" << endl;
				}
			
				for (int i = semi_arc; i< code_data.num_component_edges[0]; i++)
					component_flag[i] = -1;
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   initial component_flag with shortcut edges: ";
	for (int i=0; i< num_edges; i++)
		debug << component_flag[i] << ' ';
	debug << endl;
}
	
			}
			
			algebraic_crossing_number = 0;
			polynomial<int,bracket_variable> arrow_factor(1);	
			
			/* For the parity bracket we record the oriented edges belonging to each unicursal component in the rows of 
			   unicursal_component: if the state component traverses the edge forwards, it is recorded as a positive integer
			   and negative if the edge is traversed backwards.  The number of edges in each component is stored in column zero
			   and we count the number of unicursal components as we trace the state.
			   
			   We also maintain a companion matrix immersion_crossing_map that records the odd and virtual immersion crossings 
			   we encounter as we trace the unicursal components.  If traversing a unicursal state component edge e (forwards or 
			   backwards) arrives at odd or virtual immersion crossing c, then c is recorded  in immersion_crossing_map in the
			   location corresponding to the location in unicursal_component of e.
			*/
			int num_unicursal_components = 0;  // the total number of unicursal components in the state
			matrix<int> unicursal_component(2*num_crossings,2*num_crossings+1); 
			matrix<int> immersion_crossing_map(2*num_crossings,2*num_crossings+1,-1); 
			int component_index;  // used to assign unicursal_component
			
			/* For the arrow polynomial, we maintain a single zig-zag count per state component (see arrow_state_zig_zag_count below) 
			   and count the number of state components with irreducible cusps in num_zig_zag_components.
			*/
			int num_zig_zag_components = 0;

			/* For the parity arrow polynomial, we will eventually record the irreducible cusps on the edges between odd crossings in each graphical component
			   in a matrix whose columns correspond to immersion crossings, and store the counts for the edges terminating at an odd crossing in the row and
			   column of this matrix corresponding to the corresponding terminating edges.  That is, if edge e terminates at an odd crossing, e appears as 
			   either an odd or even terminating edge at crossing k, so the zig-zag count is stored in the location corresponding to the position of e in the 
			   generic_code_data::table::ODD_TERMINATING or generic_code_data::table::EVEN_TERMINATING rows of the code table.  Note, this approach means that if there are virtual or shortcut crossings between
			   odd crossings in a graphical component, then the zig-zag count between those two crossings is associated with the last edge between them in the 
			   immersion.
			   
			   When we start tracing unicursal state components, we could start anywhere on a component, so cannot immediately store the counts as desired.  
			   Instead, since we are counting the irreducible cusps on edges between odd crossings of graphical components, these edges correspond to the 
			   terms of the odd-crossing Gauss data of the state, although at this stage we do not know how many odd crossings there are on each unicursal 
			   component.  So, as we trace the state components, we note the position of the first Gauss data term on each component in first_odd_crossing_term 
			   and the position corresponding to the next odd (immersion) crossing we will encounter in next_odd_crossing_term, then count the number of zig-zags 
			   until we reach the next odd crossing, storing the counts in parity_arrow_gauss_zig_zags .  When we get to the end of a component, we need to 
			   amalgamate the count for the current next_odd_crossing_term with that of first_odd_crossing_term.
			   			   
			   Since the order in which we consider unicursal components when creating graphical components may be different to the order in which we trace
			   unicursal components initially, we note, in odd_immersion_crossing, the order in which we encounter odd immersion crossings when doing the 
			   initial trace.
			   
			   Some unicursal components may become disconnected simple closed curves as a result of the reduction of graphical components so for parity
			   arrow we also maintain the arrow_state_zig_zag_count and record the zig-zag count for each state component in unicursal_component_zig_zag_count.
			*/
			vector<int> parity_arrow_gauss_zig_zags(2*num_odd_crossings+1);  // +1 to allow for the end of the last component and knotoids
			vector<int> odd_immersion_crossing(2*num_odd_crossings);
			int first_odd_crossing_term = 0;
			int next_odd_crossing_term = 0;
			
			vector<int> unicursal_component_zig_zag_count(2*num_crossings);
			
			int parity_arrow_last_non_seifert_turn = generic_code_data::VOID;
			int knotoid_head_zig_zag_count = 0;

			bool components_complete;
			do
			{
				int start;
				components_complete = true;
				
				/* look for the start of another component, the first component will always be found with start == 0 */
	
				for (int i=0; i< num_edges; i++)
				{
					if (component_flag[i] == 0)
					{
						components_complete = false;
						num_unicursal_components++;
						start = i;
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial: new component found starting at edge " << start << endl;

if (variant == ARROW_VARIANT && debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial: arrow_factor currently " << arrow_factor << endl;
	
						break;
					}
				}
				
				if (!components_complete)
				{					
					/* trace the current component */
					int edge=start;
					bool forwards = true;				
					bool component_traced = false;				
					int next_crossing;

					
					/* with the arrow polynomial variant, we maintain a single zig-zag count per state component */
					arrow_state_zig_zag_count = 0;
					component_index = 1;
					
					bool encountered_odd_crossing_on_component = false;
					
					do
					{
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:   edge " << edge;
						
						/* mark the flag with the current component number, counted from 1 (to aid debugging) */
						component_flag[edge] = num_unicursal_components;
						unicursal_component[num_unicursal_components-1][component_index] = edge;
						
						if (forwards)
						{
							next_crossing = term_crossing[edge];
						}
						else
						{
							next_crossing = code_data.orig_crossing[edge];								
							unicursal_component[num_unicursal_components-1][component_index] *= -1;
						}
						component_index++;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << ", next crossing = " << next_crossing << endl;
													
						/* Identify the next edge after this crossing. In the case of immersion or peer
						   codes we use the odd or even polarity of the edge to identify the next edge.
						   In the case of Gauss codes if we are moving forwards (backwards) with respect
						   to the edge numbering we have to identify whether we are arriving on the odd
						   terminating (originating) edge.  We use a boolean, odd_visit, to indicate 
						   whether the edge polarity is odd, or whether we are dealing with the odd
						   originating (terminating) visit.
						*/
						bool odd_visit;
						if (forwards)
						{
							if (code_data.type == generic_code_data::gauss_code)
								odd_visit = (edge == code_table[generic_code_data::table::ODD_TERMINATING][next_crossing]);
							else
								odd_visit = (edge%2);
						}
						else
						{
							if (code_data.type == generic_code_data::gauss_code)
								odd_visit = (edge == code_table[generic_code_data::table::ODD_ORIGINATING][next_crossing]);
							else
								odd_visit = (edge%2);
						}
						
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:     edge " << edge << " is an " << (odd_visit? "odd" : "even") 
		  << " visit to crossing " << next_crossing << " with crossing_record = " << crossing_record[next_crossing] << endl;
}
						if (crossing_record[next_crossing] == VIRTUAL_CROSSING)
						{
							immersion_crossing_map[num_unicursal_components-1][component_index-1]= next_crossing;							
							
							/* just move on around the immersion.  We have to evaluate the component of the underlying
							   diagram on which the current edge lies.  This is not the same as the component of the smoothed
							   diagram recorded in components.
							   
							   In the case of a virtual knot or knotoid, the generic_code_data::table::COMPONENT row will correctly be set to zero for all crossings
							*/
							int component;
							if (edge%2)
								component = code_table[generic_code_data::table::COMPONENT][(edge-1)/2];
							else
								component = code_table[generic_code_data::table::COMPONENT][edge/2];
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     virtual crossing, edge " << edge << " lies on component " << component <<  endl;

							if (forwards)
								edge = (edge + 1 - first_edge_on_component[component])% num_component_edges[component] + first_edge_on_component[component];		   
							else
								edge = (edge + num_component_edges[component] - 1 - first_edge_on_component[component])%
								num_component_edges[component] + first_edge_on_component[component];		   
								
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     virtual crossing, move onto edge " << edge << (forwards? " forwards" : " backwards") << endl;
						}
						else if ((variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT) && crossing_parity[next_crossing] == gauss_orientation_data::parity::ODD)
						{								
							encountered_odd_crossing_on_component = true;
							immersion_crossing_map[num_unicursal_components-1][component_index-1] = next_crossing;
							
							/* Move on around the component */
							int component;
							if (edge%2)
								component = code_table[generic_code_data::table::COMPONENT][(edge-1)/2];
							else
								component = code_table[generic_code_data::table::COMPONENT][edge/2];
								
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     odd crossing, edge " << edge << " lies on component " << component <<  endl;

							if (forwards)
								edge = (edge + 1 - first_edge_on_component[component])% num_component_edges[component] + first_edge_on_component[component];		   
							else
								edge = (edge + num_component_edges[component] - 1 - first_edge_on_component[component])%
								num_component_edges[component] + first_edge_on_component[component];		   
																
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     odd crossing, move onto edge " << edge << (forwards? " forwards" : " backwards") << endl;	
				
							/* for the parity arrow polynomial variant, reset the last non-Seifert turn records 
							   and increment the next odd crossing term marker
							*/
							parity_arrow_last_non_seifert_turn = generic_code_data::VOID;
							odd_immersion_crossing[next_odd_crossing_term] = next_crossing;
							next_odd_crossing_term++;														
							
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     next_odd_crossing_term updated to " << next_odd_crossing_term << endl;	
						}
						else if (crossing_record[next_crossing] == SEIFERT_SMOOTHED)
						{
							/* Seifert-smoothed moves onto the adjacent semi-arc 
							   without changing direction, that is 
							   2i <-> 2j and 2j-1 <-> 2i+1
							*/
							if (forwards)
							{
								if (odd_visit)
									edge = code_table[generic_code_data::table::ODD_ORIGINATING][next_crossing];
								else
									edge = code_table[generic_code_data::table::EVEN_ORIGINATING][next_crossing];
							}
							else
							{
								if (odd_visit)
									edge = code_table[generic_code_data::table::ODD_TERMINATING][next_crossing];
								else
									edge = code_table[generic_code_data::table::EVEN_TERMINATING][next_crossing];

							}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     Seifert smoothed, move onto edge " << edge << (forwards? " forwards" : " backwards") << endl;
						}
						else if (crossing_record[next_crossing] == NON_SEIFERT_SMOOTHED)
						{
							/* non-Seifert-smoothed moves onto the adjacent semi-arc 
							   that changes direction, that is 2i <-> 2j-1 and 2j <-> 2i+1
							   
							   determine whether we are turning left or right for the evaluation of the arrow polynomial
							*/
							if (forwards)
							{
								if (odd_visit)
								{
									edge = code_table[generic_code_data::table::EVEN_TERMINATING][next_crossing];
									non_seifert_turn = (code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1? generic_code_data::RIGHT: generic_code_data::LEFT);
								}
								else
								{
									edge = code_table[generic_code_data::table::ODD_TERMINATING][next_crossing];
									non_seifert_turn = (code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1? generic_code_data::LEFT: generic_code_data::RIGHT);
									
								}

								forwards = false;
							}
							else
							{
								if (odd_visit)
								{
									edge = code_table[generic_code_data::table::EVEN_ORIGINATING][next_crossing];
									non_seifert_turn = (code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1? generic_code_data::LEFT: generic_code_data::RIGHT);
								}
								else
								{
									edge = code_table[generic_code_data::table::ODD_ORIGINATING][next_crossing];
									non_seifert_turn = (code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1? generic_code_data::RIGHT: generic_code_data::LEFT);
								}
								
								forwards = true;
							}
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:     non-Seifert smoothed, move onto edge " << edge << (forwards? " forwards" : " backwards") << endl;
	debug << "bracket_polynomial:     non-seifert_turn = " << (non_seifert_turn == generic_code_data::LEFT? "LEFT" : "RIGHT") << endl;
}
							/* for the (parity) arrow polynomial, consecutive non-Seifert smoothed left or right turns cancel */
							if (variant == ARROW_VARIANT || variant == PARITY_ARROW_VARIANT)
							{
								if (arrow_state_zig_zag_count == 0)
								{
									arrow_state_zig_zag_count++;
									arrow_state_last_non_seifert_turn = non_seifert_turn;
								}
								else if (non_seifert_turn == arrow_state_last_non_seifert_turn)
								{
									arrow_state_zig_zag_count--;
									
									/* the last non-Seifert turn needs to be inverted, since we have just cancelled a consecutive pair */
									if (arrow_state_zig_zag_count == 0)
										arrow_state_last_non_seifert_turn = generic_code_data::VOID;
									else if (arrow_state_last_non_seifert_turn == generic_code_data::LEFT)
										arrow_state_last_non_seifert_turn = generic_code_data::RIGHT;
									else
										arrow_state_last_non_seifert_turn = generic_code_data::LEFT;
								}
								else
								{
									arrow_state_zig_zag_count++;									
									arrow_state_last_non_seifert_turn = non_seifert_turn;
								}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
{
	debug << "bracket_polynomial:     arrow_state_zig_zag_count = " << arrow_state_zig_zag_count << endl;
	debug << "bracket_polynomial:     arrow_state_last_non_seifert_turn = ";
	switch(arrow_state_last_non_seifert_turn)
	{
		case generic_code_data::LEFT: debug << "LEFT"; break;
		case generic_code_data::RIGHT: debug << "RIGHT"; break;
		case generic_code_data::VOID: debug << "VOID"; break;
	}
	debug << endl;
}
							}
							
							if (variant == PARITY_ARROW_VARIANT)
							{																									
								if (parity_arrow_gauss_zig_zags[next_odd_crossing_term] == 0)
								{
									parity_arrow_gauss_zig_zags[next_odd_crossing_term]++;
									
									if (non_seifert_turn == generic_code_data::LEFT)
										parity_arrow_gauss_zig_zags[next_odd_crossing_term] *= -1;
										
									parity_arrow_last_non_seifert_turn = non_seifert_turn;
								}
								else if (non_seifert_turn == parity_arrow_last_non_seifert_turn)
								{
									if (parity_arrow_gauss_zig_zags[next_odd_crossing_term] > 0)
										parity_arrow_gauss_zig_zags[next_odd_crossing_term]--;
									else
										parity_arrow_gauss_zig_zags[next_odd_crossing_term]++;
									
									/* the last non-Seifert turn needs to be inverted, since we have just cancelled a consecutive pair */
									if (parity_arrow_last_non_seifert_turn == 0)
										parity_arrow_last_non_seifert_turn = generic_code_data::VOID;																			
									else if (parity_arrow_last_non_seifert_turn == generic_code_data::LEFT)
										parity_arrow_last_non_seifert_turn = generic_code_data::RIGHT;
									else
										parity_arrow_last_non_seifert_turn = generic_code_data::LEFT;									
								}
								else
								{
									if (parity_arrow_gauss_zig_zags[next_odd_crossing_term] > 0)
										parity_arrow_gauss_zig_zags[next_odd_crossing_term]++;
									else
										parity_arrow_gauss_zig_zags[next_odd_crossing_term]--;

									parity_arrow_last_non_seifert_turn = non_seifert_turn;
								}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
{
	debug << "bracket_polynomial:     parity_arrow_gauss_zig_zags[" << next_odd_crossing_term << "] = " << parity_arrow_gauss_zig_zags[next_odd_crossing_term] << endl;
	debug << "bracket_polynomial:     parity_arrow_gauss_zig_zags: ";
	for (int i=0; i< 2*num_odd_crossings+1; i++)
		debug << parity_arrow_gauss_zig_zags[i] << ' ';
	debug << endl;
	debug << "bracket_polynomial:     parity_arrow_last_non_seifert_turn = ";
	switch(parity_arrow_last_non_seifert_turn)
	{
		case generic_code_data::LEFT: debug << "LEFT"; break;
		case generic_code_data::RIGHT: debug << "RIGHT"; break;
		case generic_code_data::VOID: debug << "VOID"; break;
	}
	debug << endl;
}
							}
						}
						else // shortcut crossing
						{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     shortcut crossing" <<  endl;

							/* In this case we must have pure_knotoid_code_data and we've arrived at a shortcut 
							   crossing so if we're evaluating Turaev's knotoid extended bracket and we're on 
							   the segment component (i.e. if start == 0) we first check whether we've reached 
							   the head, in which case we stop tracing the component.  Otherwise we evaluate 
							   the algebraic number of intersections of segment with the shortcut.  This is the 
							   number of times the segment crosses the shortcut from right to left minus the 
							   number of times it crosses from left to right. 
							   
							   We are traversing the segment rather than the shortcut, so we'll always 
							   approach a shortcut crossing on the over-arc, but we may be going forwards
							   or backwards.  Also, the shortcut is oriented from the leg to the head 
							   for the purpuses of the algebraic number of intersections, rather than
							   the implicit head to leg orientation given by the semi-arc numbering.  
							   
							   Thus, if we arrive at a type 1 crossing on an odd edge going forwards, it is a 
							   negative crossing as far as the peer code is concerned and the sortcut 
							   should be considered oriented from 2i+1 -> 2i, so the segment crosses from 
							   the right.  If on the other hand we are going backwards we are entering 
							   a positive crossing (because we're on the over-arc) and the shortcut is 
							   oriented 2j -> 2j-1, so the segment again crosses from the right.  
							   
							   By a similar reasoning we have the following:
							   
								   edge    direction     crossing    PC sign  segment crosses from
								   odd     forwards       type 1      -ve             right
								   odd     backwards      type 1      +ve             right
								   odd     forwards       type 2      -ve             left
								   odd     backwards      type 2      +ve             left
								   even    forwards       type 1      +ve             left
								   even    backwards      type 1      -ve             left
								   even    forwards       type 2      +ve             right
								   even    backwards      type 2      -ve             right
							*/

							if (variant == TURAEV_VARIANT)
							{
								if (edge == head_semi_arc)  // for links head_semi_arc is set to num_edges, so we'll never reach here
								{
									component_traced = true;
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     head semi-arc reached, stop tracing segment component" << endl;
									break;
								}
								else if (start == 0)
								{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:     on segment component, evaluate contribution to algebraic_crossing_number" << endl;
									if (odd_visit)
									{	
										if (code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1)				
										{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:     segment crosses shortcut from the right" << endl;
											algebraic_crossing_number += 1;
										}
										else									
										{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:     segment crosses shortcut from the left" << endl;
											algebraic_crossing_number -= 1;
										}
									}
									else
									{
										if (code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1)				
										{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:     segment crosses shortcut from the left" << endl;
											algebraic_crossing_number -= 1;
										}
										else									
										{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:     segment crosses shortcut from the right" << endl;
											algebraic_crossing_number += 1;
										}
									}
								}
								else
								{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:     not on segment component, do nothing" << endl;
								}
							}	
	
							immersion_crossing_map[num_unicursal_components-1][component_index-1]= next_crossing;							
							
							/* just move on around the immersion.  We have to evaluate the component of the underlying
							   diagram on which the current edge lies.  This is not the same as the component of the smoothed
							   diagram recorded in components.
							*/
							int component;
							if (edge%2)
								component = code_table[generic_code_data::table::COMPONENT][(edge-1)/2];
							else
								component = code_table[generic_code_data::table::COMPONENT][edge/2];
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     shortcut crossing, edge " << edge << " lies on component " << component <<  endl;

							if (forwards)
								edge = (edge + 1 - first_edge_on_component[component])% num_component_edges[component] + first_edge_on_component[component];		   
							else
								edge = (edge + num_component_edges[component] - 1 - first_edge_on_component[component])%
						num_component_edges[component] + first_edge_on_component[component];		   
								
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     move onto edge " << edge << (forwards? " forwards" : " backwards") << endl;
						}
				
						if (edge == start)
						{
							unicursal_component[num_unicursal_components-1][0] = component_index-1;
							component_traced = true;							
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     finished tracing the component" << endl;
						}
						else
						{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     edge = " << edge << " start = " << start << " continuing" << endl;
						}
					} while (!component_traced);

					if (variant == ARROW_VARIANT)
					{
						/* evaluate the contribution to the arrow factor from this component
						
						   The number of cusps on a loop component must be even, since passing through a cusp reverses
						   our forwards or backwards motion with respect to the underlying orientation, so in order to return 
						   to the starting edge we must always go through an even number of such reversals.
						   
						   That means that wherever we start tracing a component, by the time we return to the start, we cannot have the
						   first and the last non-Seifert turns having the same parity, since we have removed consecutive pairs of like turns
						   and so are left with an even number of turns of alternating parity.
						*/				
if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "bracket_polynomial:     arrow_state_zig_zag_count = " << arrow_state_zig_zag_count << endl;
//						if ((pure_knotoid_code_data || code_data.immersion_character == generic_code_data::character::LONG_KNOT || code_data.immersion_character == generic_code_data::character::KNOT_TYPE_KNOTOID) && start == 0)
						if (code_data.immersion_character != generic_code_data::character::CLOSED && start == 0)
						{
							/* The segment component in a knotoid contributes an L variable */
							if (arrow_state_zig_zag_count)
							{
								num_zig_zag_components++;
								
								/* deferring setting the variable mapping and using characters a,b,c, etc.  means that the resulting polynomial is 
								   presented  with the variables in order K_1K_2, rather than the K_2K_1 that could occur if we map the variables
								   here and let the polynomial code choose the next variable character
								*/
								   
								char z = 'a'+arrow_state_zig_zag_count/2-1;
								arrow_factor *= polynomial<int,bracket_variable>(string(1,z));
							}
							else if ((code_data.immersion_character == generic_code_data::character::LONG_KNOT || code_data.immersion_character == generic_code_data::character::KNOTOID) && !braid_control::EXPANDED_BRACKET_POLYNOMIAL)
							{
								num_zig_zag_components++;  // not a zig-zag but we do want to use L rather than D for this component
								char z = 'L'; // long segment with no zig-zags
								arrow_factor *= polynomial<int,bracket_variable>(string(1,z));
							}
						}
						else
						{
							/* Non-segment components contribute a K variable */
							if (arrow_state_zig_zag_count)
							{
								num_zig_zag_components++;
								
								char z = 'n'+arrow_state_zig_zag_count/2-1;
								arrow_factor *= polynomial<int,bracket_variable>(string(1,z));
							}
						}
					}
					else if (variant == PARITY_ARROW_VARIANT)
					{		
						/* Record the arrow_state_zig_zag_count for this unicursal component in unicursal_component_zig_zag_count.
						   We can handle those unicursal components that do not contain odd crossings immediately but we also need 
						   this information for unicursal components that will get disconnected by the reduction of graphical components.
						*/
						unicursal_component_zig_zag_count[num_unicursal_components-1] = arrow_state_zig_zag_count;
						
						if (encountered_odd_crossing_on_component)
						{
//							if ((pure_knotoid_code_data || code_data.immersion_character == generic_code_data::character::LONG_KNOT || code_data.immersion_character == generic_code_data::character::KNOT_TYPE_KNOTOID) && start == 0)
							if (code_data.immersion_character != generic_code_data::character::CLOSED && start == 0)
							{
								/* note the number of zig-zags at the knotoid head, so we can assign it to the component_zig_zag_count
								   once we have evaluated the rest of the component peer code.
								*/
								knotoid_head_zig_zag_count = parity_arrow_gauss_zig_zags[next_odd_crossing_term];
								
if (debug_control::DEBUG >= debug_control::BASIC) 
	debug << "bracket_polynomial:     knotoid_head_zig_zag_count = " << knotoid_head_zig_zag_count << endl;
							}
							else
							{
								/* amalgamate the current zig-zag count for the next_odd_crossing_term with that for the first_odd_crossing_term,
								   clear the zig-zag count for next_odd_crossing_term ready and update first_odd_crossing_term, for any subsequent component
								*/
if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "bracket_polynomial:   amalgamate parity arrow zig-zag counts " << next_odd_crossing_term << " and " << first_odd_crossing_term << endl;
	
								int amalgamated_first_zig_zag_count = amalgamate_zig_zag_counts(parity_arrow_gauss_zig_zags[next_odd_crossing_term],parity_arrow_gauss_zig_zags[first_odd_crossing_term]);
								
if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "bracket_polynomial:   amalgamated_first_zig_zag_count =  " << amalgamated_first_zig_zag_count << endl;
										
								parity_arrow_gauss_zig_zags[first_odd_crossing_term] = amalgamated_first_zig_zag_count;
								
if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "bracket_polynomial:   parity_arrow_gauss_zig_zags[" << first_odd_crossing_term << "] updated to " << parity_arrow_gauss_zig_zags[first_odd_crossing_term] << endl;						
							}
							
							first_odd_crossing_term = next_odd_crossing_term;
						}
						else
						{
							/* we have a unicursal component containing no odd cossings, so treat this component as in the 
							   arrow polynomial, based on the arrow_state_zig_zag_count
							*/
														
							if (code_data.immersion_character != generic_code_data::character::CLOSED && start == 0)
							{
								/* The segment component in a knotoid contributes an L variable */
								if (arrow_state_zig_zag_count)
								{
									num_zig_zag_components++;
									
									char z = 'a'+arrow_state_zig_zag_count/2-1;
									arrow_factor *= polynomial<int,bracket_variable>(string(1,z));
								}
								else if (!braid_control::EXPANDED_BRACKET_POLYNOMIAL)
								{
									num_zig_zag_components++;  // not a zig-zag but we do want to use L rather than D for this component
									char z = 'L'; // long segment with no zig-zags
									arrow_factor *= polynomial<int,bracket_variable>(string(1,z));
								}								
							}
							else
							{
								/* Non-segment components contribute a K variable */
								if (arrow_state_zig_zag_count)
								{
									num_zig_zag_components++;
									
									char z = 'n'+arrow_state_zig_zag_count/2-1;
									arrow_factor *= polynomial<int,bracket_variable>(string(1,z));									
								}
							}
							
if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "bracket_polynomial:   no odd crossings in component, evaluate arrow_factor based on arrow_state_zig_zag_count = " << arrow_state_zig_zag_count << endl;						
						}
						
						/* we have to reset the next_odd_crossing_term counter whether or not there was an odd crossing on the last component */
						parity_arrow_gauss_zig_zags[next_odd_crossing_term] = 0;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
{
	debug << "bracket_polynomial:   parity_arrow_gauss_zig_zags: ";
	for (int i=0; i< 2*num_odd_crossings+1; i++)
		debug << parity_arrow_gauss_zig_zags[i] << ' ';
	debug << endl;
}						
					}
				}					
			} while (!components_complete);
			
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial: component_flag: ";
	for (int i=0; i< num_edges; i++)
		debug << component_flag[i] << ' ';
	debug << endl;
	debug << "bracket_polynomial: number of unicursal components = " << num_unicursal_components << endl;

	if (variant == ARROW_VARIANT || variant == PARITY_ARROW_VARIANT)
	{
		debug << "bracket_polynomial: num_zig_zag_components = " << num_zig_zag_components << endl;
		debug << "bracket_polynomial: arrow_factor after tracing components = " << arrow_factor << endl;
	}
	
	if (variant == TURAEV_VARIANT)
		debug << "bracket_polynomial: algebraic_crossing_number = " << algebraic_crossing_number << endl;
	
	if (variant == PARITY_ARROW_VARIANT)
	{
		debug << "bracket_polynomial: parity_arrow_gauss_zig_zags: ";
		for (int i=0; i< 2*num_odd_crossings+1; i++)
			debug << parity_arrow_gauss_zig_zags[i] << ' ';
		debug << endl;
		debug << "bracket_polynomial: odd_immersion_crossing: ";
		for (int i=0; i< 2*num_odd_crossings; i++)
			debug << odd_immersion_crossing[i] << ' ';
		debug << endl;
		debug << "bracket_polynomial: unicursal_component_zig_zag_count: ";
		for (int i=0; i< num_unicursal_components; i++)
			debug << unicursal_component_zig_zag_count[i] << ' ';
		debug << endl;
	}
}
			/* For the parity arrow polynomial, arrange the parity_arrow_gauss_zig_zags counts in a 2 x num_crossings matrix, 
			   with the count for the Gauss arc terminating at the first visit to the corresponding odd immersion crossing is in 
			   row 1 (first is odd), and that for the arc following the second visit is in row 0 (second is even).
			*/
			matrix<int> parity_arrow_state_zig_zags(2,num_crossings);
			vector<bool> first_parity_arrow_state_zig_zag_set(num_crossings);
			if (variant == PARITY_ARROW_VARIANT)
			{
				for (int i=0; i< 2*num_odd_crossings; i++)
				{
					if (!first_parity_arrow_state_zig_zag_set[odd_immersion_crossing[i]])
					{
						parity_arrow_state_zig_zags[1][odd_immersion_crossing[i]] = parity_arrow_gauss_zig_zags[i];
						first_parity_arrow_state_zig_zag_set[odd_immersion_crossing[i]] = true;
					}
					else
					{
						parity_arrow_state_zig_zags[0][odd_immersion_crossing[i]] = parity_arrow_gauss_zig_zags[i];
					}
				}
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial: parity_arrow_state_zig_zags (first visit row 1, second visit row 0):" << endl;
    print(parity_arrow_state_zig_zags, debug, 3, "bracket_polynomial:   ");
}
			}
	
		    /* Determine the value of sigma, which is based on the A-smoothing or B-smoothing
		       performed at each crossing, adding 1 for each A-smoothed crossing and -1 for
		       every B-smoothed crossing.  In the above we interpreted the state as Seifert
		       and non-Seifert smoothed, so from the table shown in the comment introducing the
		       function we see that 
		       
		       for state 1 (Seifert smoothing) for positive crossings we add  1 to sigma
		       for state 1 (Seifert smoothing) for negative crossings we add -1 to sigma
		       for state -1 (non-Seifert smoothing) for positive crossings we add -1 to sigma
		       for state -1 (non-Seifert smoothing) for negative crossings we add 1 to sigma
		       
		       Thus, the ith crossing contributes the product of the state (Seifert or non-Seifert) and 
		       the crossing sign (positive or negative) to sigma.
		       
		       We skip over any virtual crossings
		    */	
			int sigma=0;
			place=0;
			for (int i=0; i< num_crossings; i++)
			{
				if ((variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT) && crossing_parity[i] != gauss_orientation_data::parity::EVEN)
					continue;
					
				if (sign[i] == generic_braid_data::crossing_type::POSITIVE)
				{
					sigma += state[place];
					place++;
				}
				else if (sign[i] == generic_braid_data::crossing_type::NEGATIVE)
				{
					sigma += state[place] * -1;
					place++;
				}
			}
			
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial: sigma = " << sigma << endl;

			int num_non_graphical_u_cpts = 0;
			int num_virtual_components = 0;
			polynomial<int,bracket_variable> parity_factor(1); // holds graphical component variables
			polynomial<int,bracket_variable> relaxed_parity_factor(1); // holds M and J variables
			
			if (variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT)
			{
if (debug_control::DEBUG >= debug_control::BASIC)
{
	for (int i=0; i< num_unicursal_components; i++)
	{
		debug << "bracket_polynomial: unicursal_component " << i << " length " << unicursal_component[i][0] << ": ";		
		for (int j=1; j<= unicursal_component[i][0]; j++)
			debug << unicursal_component[i][j] << " ";
		debug << endl;
		debug << "bracket_polynomial: immersion_crossing_map " << i << " length " << unicursal_component[i][0] << ": ";		
		for (int j=1; j<= unicursal_component[i][0]; j++)
			debug << immersion_crossing_map[i][j] << " ";
		debug << endl;
		
	}
}
				/* determine the number of graphical components and the correspondance between unicursal components 
				   and grahical components.  We trace each unicursal state component and if we encounter an odd 
				   crossing, we note the terminating peer at that crossing and then look for that edge in the state 
				   components to see if it belongs to a different unicursal component that is part of the same graphical 
				   compnent as the one we're tracing.
				   
				   As we trace the unicursal state components, we record the number of graphical component crossings (odd
				   virtual or shortcut crossings) that we encounter in num_unicursal_cpt_g_crossings.
				*/
				vector<int> unicursal_cpt_to_g_component_map(num_unicursal_components, -1);
				vector<int> num_unicursal_cpt_g_crossings(num_unicursal_components);
				int num_graphical_components = 0;
				
				for (int i=0; i< num_unicursal_components; i++)
				{
					int num_unicursal_component_g_crossings = 0;
					
					for (int j=1; j<= unicursal_component[i][0]; j++)
					{
						int edge = unicursal_component[i][j];
						int crossing;
						
						if (edge<0)
							crossing = orig_crossing[abs(edge)];						
						else
							crossing = term_crossing[edge];						
						
						if (crossing_parity[crossing] == gauss_orientation_data::parity::ODD)
						{							
							num_unicursal_component_g_crossings++;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "bracket_polynomial: edge " << edge << ", crossing " << crossing 
	      << " is ODD, num_unicursal_component_g_crossings = " << num_unicursal_component_g_crossings << endl;
}							
							int peer;
							int peer_unicursal_component=0;
							
							if (edge>=0)
								peer = code_table[(edge%2?generic_code_data::table::EVEN_TERMINATING:generic_code_data::table::ODD_TERMINATING)][crossing];
							else
								peer = code_table[(edge%2?generic_code_data::table::ODD_TERMINATING:generic_code_data::table::EVEN_TERMINATING)][crossing];
																
							bool found = false;
							for (int k=0; k< num_unicursal_components; k++)
							{
								for (int l=1; l<= unicursal_component[k][0]; l++)
								{
									if (abs(unicursal_component[k][l]) == peer)
									{
										peer_unicursal_component = k;
										found = true;
										break;
									}
								}
								if (found)
									break;
							}
							
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:   peer edge to " << edge << " is " << peer << ", peer_unicursal_component = " << peer_unicursal_component << endl;
							
							if (unicursal_cpt_to_g_component_map[i] == -1)
							{
								if (unicursal_cpt_to_g_component_map[peer_unicursal_component] == -1)
								{
									/* new graphical component */
									unicursal_cpt_to_g_component_map[i] = num_graphical_components;								
									num_graphical_components++;
								}
								else
								{
									unicursal_cpt_to_g_component_map[i] = unicursal_cpt_to_g_component_map[peer_unicursal_component];
								}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:   unicursal component " << i << " is part of graphical component " << unicursal_cpt_to_g_component_map[i] << endl;
							}
							
							if (unicursal_cpt_to_g_component_map[peer_unicursal_component] == -1)
							{
								unicursal_cpt_to_g_component_map[peer_unicursal_component] = unicursal_cpt_to_g_component_map[i];
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:   unicursal component " << peer_unicursal_component << " is part of graphical component " << unicursal_cpt_to_g_component_map[peer_unicursal_component] << endl;
							}
							
						}
						else if (code_table[generic_code_data::table::LABEL][crossing] == generic_code_data::VIRTUAL)
						{
							num_unicursal_component_g_crossings++;
							
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "bracket_polynomial: edge " << edge << ", crossing " << crossing 
	      << " is VIRTUAL, num_unicursal_component_g_crossings = " << num_unicursal_component_g_crossings << endl;
}
						}
						else if (pure_knotoid_code_data && shortcut_crossing[crossing])
						{
							num_unicursal_component_g_crossings++;
							
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "bracket_polynomial: edge " << edge << ", crossing " << crossing 
	      << " is a shortcut crossing, num_unicursal_component_g_crossings = " << num_unicursal_component_g_crossings << endl;
}
						}
						else
						{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "bracket_polynomial: edge " << edge << ", crossing " << crossing 
	      << " is an even crossing, num_unicursal_component_g_crossings = " << num_unicursal_component_g_crossings << endl;
}
						}
					}
					num_unicursal_cpt_g_crossings[i] = num_unicursal_component_g_crossings;
				}
				
				/* count the number of unicursal state components corresponding to each graphical component */
				vector<int> graphical_cpt_u_component_count(num_graphical_components);
				for (int i=0; i< num_graphical_components; i++)
				{
					for (int j=0; j< num_unicursal_components; j++)
					{
						if (unicursal_cpt_to_g_component_map[j] == i)
							graphical_cpt_u_component_count[i]++;
					}
				}
				
				/* identify the number of unicursal state components that are not part of a graphical component */
				num_non_graphical_u_cpts = num_unicursal_components;
				for (int i=0; i< num_graphical_components; i++)
					num_non_graphical_u_cpts -= graphical_cpt_u_component_count[i];
				
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial: num_unicursal_cpt_g_crossings: ";
	for (int i=0; i< num_unicursal_components; i++)
		debug << num_unicursal_cpt_g_crossings[i] << " ";
	debug << endl;
	debug << "bracket_polynomial: num_graphical_components = " << num_graphical_components << endl;
	debug << "bracket_polynomial: unicursal_cpt_to_g_component_map: ";
	for (int i=0; i< num_unicursal_components; i++)
		debug << unicursal_cpt_to_g_component_map[i] << " ";
	debug << endl;
	debug << "bracket_polynomial: graphical_cpt_u_component_count: ";
	for (int i=0; i< num_graphical_components; i++)
		debug << graphical_cpt_u_component_count[i] << " ";
	debug << endl;
	debug << "bracket_polynomial: num_non_graphical_u_cpts = " << num_non_graphical_u_cpts << endl;
}
				/* count the number of (odd, virtual or shortcut) crossings on each graphical component
				
				   Note that a graphical component may intersect a unicursal component that is not part of
				   that graphical component (the intersection being comprised of virtual or shortcut crossings),
				   so we have to look for the number of self intersections amongst the unicursal components that
				   make up the graphical component.  We count the occurrences of each immersion crossing in the
				   unicursal components making up the graphical component and then count the number of crossings
				   with a crossing_count of 2.
				   
				   We clear the crossings whose crossing_count is not 2 in immersion_crossing_map by setting
				   the corresponding entry to -1, since they do not represents crossings in the graphical component
				*/
				vector<int> num_graphical_cpt_crossings(num_graphical_components);
				for (int i=0; i< num_graphical_components; i++)
				{
					vector<int> crossing_count(num_crossings);

					for (int j=0; j< num_unicursal_components; j++)
					{
						if (unicursal_cpt_to_g_component_map[j] != i)
							continue;
							
						for (int k=1; k <= unicursal_component[j][0]; k++)
						{
							if (immersion_crossing_map[j][k] >=0)
								crossing_count[immersion_crossing_map[j][k]]++;
						}
					}
					
					for (int j=0; j< num_crossings; j++)
					{
						if (crossing_count[j] == 2)
							num_graphical_cpt_crossings[i]++;
					}					
					
					for (int j=0; j< num_unicursal_components; j++)
					{
						if (unicursal_cpt_to_g_component_map[j] != i)
							continue;
							
						for (int k=1; k <= unicursal_component[j][0]; k++)
						{
							if (immersion_crossing_map[j][k] >=0 && crossing_count[immersion_crossing_map[j][k]] != 2)
								immersion_crossing_map[j][k] = -1;
						}
					}
				}
					
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial: num_graphical_cpt_crossings: ";
	for (int i=0; i< num_graphical_components; i++)
		debug << num_graphical_cpt_crossings[i] << " ";
	debug << endl;
	
	for (int i=0; i< num_unicursal_components; i++)
	{
		debug << "bracket_polynomial: graphical component self-crossing immersion_crossing_map " << i << " length " << unicursal_component[i][0] << ": ";		
		for (int j=1; j<= unicursal_component[i][0]; j++)
			debug << immersion_crossing_map[i][j] << " ";
		debug << endl;		
	}
}				
				/* now we can evaluate a peer code for each graphical component by creating a component generic code object 
				   with enough information to write the component's peer code to an ostringstream and then read it back again.
				   We need the head, num_crossings and the generic_code_data::table::OPEER, generic_code_data::table::TYPE, generic_code_data::table::LABEL and generic_code_data::table::COMPONENT rows of the code_table in the 
				   generic_code_data structure to do this.
				*/
				for (int i=0; i< num_graphical_components; i++)
				{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial: graphical component " << i << endl;
	
					/* We start by writing the terminating edges into a component_code_table matrix indexed by the order in which we 
					   encounter the odd, virtual and shortcut crossings in the unicursal state components that correspond to this 
					   graphical component.  The first visit is written in the ODD-TERMINATING row and the second visit in the 
					   EVEN-TERMINATING row.  We set up a map object whose key is an immersion crossing number and whose value is the
					   index in the component_code_table.  We refer to these immersion crossings as "markers"
					*/
					map<int, int> marker_crossing_map;
					
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "bracket_polynomial:   calculating marker_crossing_map" << endl;
					int index=0;					
					for (int j=0; j< num_unicursal_components; j++)
					{
						if (unicursal_cpt_to_g_component_map[j] != i)
							continue;
							
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "bracket_polynomial:     unicursal component " << j << endl;

							
						for (int k=1; k <= unicursal_component[j][0]; k++)
						{
							if (immersion_crossing_map[j][k] >= 0 && marker_crossing_map.find(immersion_crossing_map[j][k]) == marker_crossing_map.end())
							{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "bracket_polynomial:       add crossing " << immersion_crossing_map[j][k] << " with index " << index << endl;
								marker_crossing_map[immersion_crossing_map[j][k]] = index;
								index++;
							}
						}						
					}
					
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   marker_crossing_map: ";
	map<int, int>::iterator mptr = marker_crossing_map.begin();
	while (mptr != marker_crossing_map.end())
	{
		debug << "(" << mptr->first << "," << mptr->second << ") ";
		mptr++;
	}
	debug << endl;
}				
					int num_cpt_crossings = num_graphical_cpt_crossings[i];
					matrix<int> component_code_table(generic_code_data::table::CODE_TABLE_SIZE,num_cpt_crossings,-1);

					int edge=0;
					int component = -1;
					int component_head_semi_arc = -1;
					vector<int> first_edge_on_g_component(graphical_cpt_u_component_count[i]);
					
					/* Record the zig-zag counts for each Gauss arc in this graphical component; we store the counts in a matrix whose entries 
					   correspond to the EVEN-TERMINATING and ODD-TERMINATING edges in the component_code_table, that is, with the count for the 
					   first visit in row 1 (first is odd, so store in odd row), and that for the second visit in row 0.  By storing the counts 
					   this way, if we need to re-align the edge labels we do not need to update the zig-zag counts.
					*/
					matrix<int> initial_component_zig_zag_count(2,num_cpt_crossings);

					/* consider each unicursal state component that corresponds to this graphical component 
					   and note the first graphical component edge on each of its components.
					*/
					bool component_includes_unicursal_cpt_zero = false;
					for (int j=0; j< num_unicursal_components; j++)
					{
						if (unicursal_cpt_to_g_component_map[j] != i)
						{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:   unicursal component " << j << " not part of grahical component " << i << endl;
							continue;
						}
						
						component++;
						first_edge_on_g_component[component] = edge;								
						
						if (j==0)  // we have not continued, so unicursal component zero is part of this graphical component
							component_includes_unicursal_cpt_zero = true;
							
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:   unicursal component " << j << " forms component " << component << " of grahical component " << i << endl;

						/* write the unicursal state component's graphical crossings into the component_code_table */
						for (int k=1; k <= unicursal_component[j][0]; k++)
						{
							/* If the code_data we're working with is that of a pure knotoid, and we're tracing unicursal component zero, 
							   look for a head semi-arc in the graphical component. It is possible (as in the case of state 1 -1 -1 of
							   [11 -7 9 -1 -13 15 3 -5^]/+ * + + + + + -) that a smoothed state produces a knot-type knotoid graphical 
							   component, so we look for a self crossing of the graphical component within the shortcut.  The first such 
							   self crossing will be the component_head_semi_arc and its existence will indicate whether the graphical 
							   component is a pure or knot-type knotoid. 
							   
							   Note that the old head semi-arc may not correspond to the head semi-ard of the graphical component, as in 
							   the case of state -1 -1 of [7 -11 -13 -9 15 -1 -5 -3^]/+ + + + + * * -
							*/
							if (pure_knotoid_code_data && j==0)
							{
							    if (immersion_crossing_map[j][k] >=0 && abs(unicursal_component[j][k]) >= head_semi_arc)
							    {
									if (component_head_semi_arc == -1) // not set yet
										component_head_semi_arc = edge;
										
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:     knotoid: shortcut edge " << abs(unicursal_component[j][k]) << " is component_head_semi_arc = " 
	      << component_head_semi_arc << endl;
}
								}	
						    }
								
							if (immersion_crossing_map[j][k] >=0)
							{
								int graphical_cpt_crossing = marker_crossing_map[immersion_crossing_map[j][k]];
								if (component_code_table[generic_code_data::table::ODD_TERMINATING][graphical_cpt_crossing] < 0)
								{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:     first visit to graphical component crossing " << graphical_cpt_crossing << " on component edge " << edge << endl;
	
									component_code_table[generic_code_data::table::ODD_TERMINATING][graphical_cpt_crossing] = edge++;
									
									if (crossing_parity[immersion_crossing_map[j][k]] == gauss_orientation_data::parity::ODD)
										initial_component_zig_zag_count[1][graphical_cpt_crossing] = parity_arrow_state_zig_zags[1][immersion_crossing_map[j][k]];
								}
								else
								{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:     second visit to graphical component crossing " << graphical_cpt_crossing << " on component edge " << edge << endl;
	
									component_code_table[generic_code_data::table::EVEN_TERMINATING][graphical_cpt_crossing] = edge++;

									if (crossing_parity[immersion_crossing_map[j][k]] == gauss_orientation_data::parity::ODD)
										initial_component_zig_zag_count[0][graphical_cpt_crossing] = parity_arrow_state_zig_zags[0][immersion_crossing_map[j][k]];									
								}																
							}
						}						
					}

if (debug_control::DEBUG >= debug_control::BASIC)
{	
	debug << "bracket_polynomial:   first_edge_on_g_component: ";
	for (int l=0; l < graphical_cpt_u_component_count[i]; l++)
		debug << first_edge_on_g_component[l] << " ";
	debug << endl;
	debug << "bracket_polynomial:   initial first visit to graphical component crossings: ";
	for (int l=0; l < num_cpt_crossings; l++)
		debug << component_code_table[generic_code_data::table::ODD_TERMINATING][l] << " ";
	debug << endl;
	debug << "bracket_polynomial:   initial second visit to graphical component crossings: ";
	for (int l=0; l < num_cpt_crossings; l++)
		debug << component_code_table[generic_code_data::table::EVEN_TERMINATING][l] << " ";
	debug << endl;
	debug << "bracket_polynomial:   initial_component_zig_zag_count: : " << endl;
	print(initial_component_zig_zag_count,debug, 4,"bracket_polynomial:     ");
}						

					/* align the numbering of the components to ensure there is an odd and even terminating edge label at each crossing.
					   We take the front component from the component_tree and trace that component.  If we encounter a crossing that is 
					   mis-aligned, we cycle the edge labels on the peer component, and push that component number onto the back of the 
					   component tree.  In this manner we conduct a depth-first search of the component tree.
					*/
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:   aligning component numbering"<< endl;
					vector<int> aligned_labels(graphical_cpt_u_component_count[i]);
					list<int> component_tree;
					component_tree.push_back(0);
					aligned_labels[0] = 1;
					list<int>::iterator tree_ptr = component_tree.begin();
					while (tree_ptr != component_tree.end())
					{
						int component = *tree_ptr;
						int first = first_edge_on_g_component[component];
						int last;
						if (component == graphical_cpt_u_component_count[i]-1)
							last = 2*num_graphical_cpt_crossings[i]-1;
						else
							last = first_edge_on_g_component[component+1]-1;
							
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:     component " << component << ", first edge = " << first << ", last edge = " << last << endl;
	
						for (int edge=first; edge<=last; edge++)
						{
							int crossing;
							int peer;
							for (int k=0; k< num_cpt_crossings; k++)
							{
								if (component_code_table[generic_code_data::table::ODD_TERMINATING][k] == edge)
								{
									peer = component_code_table[generic_code_data::table::EVEN_TERMINATING][k];
									crossing = k;
									break;
								}
								else if (component_code_table[generic_code_data::table::EVEN_TERMINATING][k] == edge)
								{
									peer = component_code_table[generic_code_data::table::ODD_TERMINATING][k];
									crossing = k;
									break;
								}
							}

							int peer_component = 0;
							for (int k=1; k< graphical_cpt_u_component_count[i]; k++)
							{
								if (peer >= first_edge_on_g_component[k])
									peer_component++;
								else
									break;
							}
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:       edge " << edge << ", peer = " << peer << " lies on component " << peer_component << endl;

							if (edge%2 == peer%2)
							{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:       graphical component crossing " << crossing << " has incompatible terminating edge labels " << endl;
	
								int peer_first = first_edge_on_g_component[peer_component];
								int peer_last;
								if (peer_component == graphical_cpt_u_component_count[i]-1)
									peer_last = 2*num_graphical_cpt_crossings[i]-1;
								else
									peer_last = first_edge_on_g_component[peer_component+1]-1;
								
								int num_peer_component_edges = peer_last-peer_first+1;
								
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:       cycle the " << num_peer_component_edges << " edge labels on peer_component " << peer_component << ", peer_first = " << peer_first << " peer_last = " << peer_last << endl;
	
								for (int k=0; k< num_cpt_crossings; k++)
								{
									if (component_code_table[generic_code_data::table::ODD_TERMINATING][k] >= peer_first && component_code_table[generic_code_data::table::ODD_TERMINATING][k] <= peer_last)
									{
										component_code_table[generic_code_data::table::ODD_TERMINATING][k] = (component_code_table[generic_code_data::table::ODD_TERMINATING][k]-peer_first+1)%num_peer_component_edges+peer_first;
									}
									
									if (component_code_table[generic_code_data::table::EVEN_TERMINATING][k] >= peer_first && component_code_table[generic_code_data::table::EVEN_TERMINATING][k] <= peer_last)
									{
										component_code_table[generic_code_data::table::EVEN_TERMINATING][k] = (component_code_table[generic_code_data::table::EVEN_TERMINATING][k]-peer_first+1)%num_peer_component_edges+peer_first;
									}
								}
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{	
	debug << "bracket_polynomial:   cycled first visit to graphical component crossings: ";
	for (int l=0; l < num_cpt_crossings; l++)
		debug << component_code_table[generic_code_data::table::ODD_TERMINATING][l] << " ";
	debug << endl;
	debug << "bracket_polynomial:   cycled second visit to graphical component crossings: ";
	for (int l=0; l < num_cpt_crossings; l++)
		debug << component_code_table[generic_code_data::table::EVEN_TERMINATING][l] << " ";
	debug << endl;
}						
							}							
							
							if (aligned_labels[peer_component] == 0)
							{
								component_tree.push_back(peer_component);
								aligned_labels[peer_component] = 1;
							}
							
						}
						
						tree_ptr++;
					}

					
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   first visit to graphical component crossings: ";
	for (int l=0; l < num_cpt_crossings; l++)
		debug << component_code_table[generic_code_data::table::ODD_TERMINATING][l] << " ";
	debug << endl;
	debug << "bracket_polynomial:   second visit to graphical component crossings: ";
	for (int l=0; l < num_cpt_crossings; l++)
		debug << component_code_table[generic_code_data::table::EVEN_TERMINATING][l] << " ";
	debug << endl;
	debug << "bracket_polynomial:   write generic_code_data::table::OPEER, generic_code_data::table::COMPONENT and generic_code_data::table::TYPE: " << endl;
}						

					/* re-write the edges in rows generic_code_data::table::ODD_TERMINATING and EVEN_TRMINATING into the correct place in generic_code_data::table::OPEER and set generic_code_data::table::COMPONENT and generic_code_data::table::TYPE, 
					   re-write initial_component_zig_zag_count to component_zig_zag_count to align the zig-zag counts with the generic_code_data::table::ODD_TERMINATING
					   and generic_code_data::table::EVEN_TERMINATING edges in the correct column as determined by the naming edges.  We record the zig-zag count for the
					   edges corresponding to the ODD-TERMINATING or EVEN-TERMINATING edges in rows 1 and 0 respectively.
					*/					
					matrix<int> component_zig_zag_count(2,num_cpt_crossings,-1);
					
					for (int j=0; j< num_cpt_crossings; j++)
					{
						/* identify the marker crossing in immersion_crossing_map for graphical component crossing j */
						int marker = find_value(marker_crossing_map,j)->first;
							
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     graphical component crossing " << j << " immersion crossing marker = " << marker << endl;

						/* isolate the graphical component peers at the odd (first) and even (second) visits and 
						   note the corresponding zig-zag counts
						*/  
						int odd_c_peer = component_code_table[generic_code_data::table::ODD_TERMINATING][j];
						int even_c_peer = component_code_table[generic_code_data::table::EVEN_TERMINATING][j];
						int odd_c_zig_zag_count = initial_component_zig_zag_count[1][j];  
						int even_c_zig_zag_count = initial_component_zig_zag_count[0][j];

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:     found odd_c_peer = " << odd_c_peer << " even_c_peer = " << even_c_peer << endl;
	debug << "bracket_polynomial:     odd_c_zig_zag_count = " << odd_c_zig_zag_count << " even_c_zig_zag_count = " << even_c_zig_zag_count << endl;
}							
						/* identify the immersion peers corresponding to the graphical component peers */
						bool odd_assigned = false;
						bool done = false;
						
						/* Note, the odd_i_peer and even_i_peer are the immersion edge labels corresponding to the immersion edge incident
						   with the crossing at which odd_c_peer and even_c_peer arrive.  odd_i_peer (odd_i_peer) may be either odd or even labels, 
						   and may be positive or negative, since they are identified from unicursal_component.
						*/
						int odd_i_peer;  
						int even_i_peer;
							
//if (debug_control::DEBUG >= debug_control::BASIC)
//	debug << "bracket_polynomial:       look for marker " << marker << endl;
							
						for (int l=0; l < num_unicursal_components; l++)
						{
//if (debug_control::DEBUG >= debug_control::BASIC)
//	debug << "bracket_polynomial:         unicursal component " << l << endl;
							for (int m=1; m <= unicursal_component[l][0]; m++)
							{
//if (debug_control::DEBUG >= debug_control::BASIC)
//	debug << "bracket_polynomial:           crossing " << immersion_crossing_map[l][m];
								if (immersion_crossing_map[l][m] == marker)
								{
									if (!odd_assigned)
									{
//if (debug_control::DEBUG >= debug_control::BASIC)
//	debug << " assigning odd ocurrence" << endl;
										odd_i_peer = unicursal_component[l][m];
										odd_assigned = true;
									}
									else
									{
//if (debug_control::DEBUG >= debug_control::BASIC)
//	debug << " assigning even ocurrence" << endl;
										even_i_peer = unicursal_component[l][m];
										done = true;
										break;
									}
								}
								else
								{
//if (debug_control::DEBUG >= debug_control::BASIC)
//	debug << " no match" << endl;
								}
							}
							if (done)
								break;
						}

						/* now change over the odd and even peers to reflect parity rather than the first and second visit */
						if (even_c_peer%2) 
						{
							swap(odd_c_peer,even_c_peer);
							swap(odd_i_peer,even_i_peer);
							swap(odd_c_zig_zag_count,even_c_zig_zag_count);
						}
						
						/* now set generic_code_data::table::OPEER, generic_code_data::table::COMPONENT and generic_code_data::table::LABEL */
						component_code_table[generic_code_data::table::OPEER][even_c_peer/2] = odd_c_peer;
						component_zig_zag_count[1][even_c_peer/2] = odd_c_zig_zag_count;
						component_zig_zag_count[0][even_c_peer/2] = even_c_zig_zag_count;
						
						component_code_table[generic_code_data::table::COMPONENT][even_c_peer/2] = 0;
						for (int l=1; l < graphical_cpt_u_component_count[i]; l++)
						{
							if (even_c_peer >= first_edge_on_g_component[l])
								component_code_table[generic_code_data::table::COMPONENT][even_c_peer/2]++;
							else 
								break;
						}											
						
						/* for the generic_code_data::table::TYPE of the component peer code crossing, we have odd_c_peer arriving at the crossing on the edge
						   corresponding to odd_i_peer (which may be either odd or even, and either positive or negative), similarly for
						   even_c_peer.  
						   
						   2j-1 \ D / 2i+1        2i \ D / 2j
								 \ /                  \ /         
							   A  X  C              A  X  C    
								 / \                  / \
							 2i / B \ 2j        2j-1 / B \ 2i+1
							 
						   The combinations of odd_i_peer and even_i_peer correspond to the four regions around the crossing:
						   if both are positive, the c-peers are arriving at region A for both immersion type I and type II crossings.  Thus:
							- if the odd_i_peer is odd
								* if the immersion crossing is generic_code_data::table::TYPE 1, the component crossing must be generic_code_data::table::TYPE 1.
								* if the immersion crossing is generic_code_data::table::TYPE 2, the component crossing must be generic_code_data::table::TYPE 2.
							- if the odd_i_peer is even
								* if the immersion crossing is generic_code_data::table::TYPE 1, the component crossing must be generic_code_data::table::TYPE 2.
								* if the immersion crossing is generic_code_data::table::TYPE 2, the component crossing must be generic_code_data::table::TYPE 1.
							
						   if both are negative, the c-peers are arriving at region C for both immersion type I and type II crossings. Thus:
							- if the odd_i_peer is odd
								* if the immersion crossing is generic_code_data::table::TYPE 1, the component crossing must be generic_code_data::table::TYPE 2.
								* if the immersion crossing is generic_code_data::table::TYPE 2, the component crossing must be generic_code_data::table::TYPE 1.
							- if the odd_i_peer is even
								* if the immersion crossing is generic_code_data::table::TYPE 1, the component crossing must be generic_code_data::table::TYPE 1.
								* if the immersion crossing is generic_code_data::table::TYPE 2, the component crossing must be generic_code_data::table::TYPE 2.
						   
						   if the odd_i_peer is positive and the even_i_peer is negative, then the odd_i_peer (and consequently the odd_c_peer) 
						   must be on the left side of the above diagrams and the even_i_peer (and consequently the even_c_peer) must be on the 
						   right side).  Thus:
							- if the odd_i_peer is odd
								* if the immersion crossing is generic_code_data::table::TYPE 1, the c_peers arrive at region D, and the component crossing must be generic_code_data::table::TYPE 2.
								* if the immersion crossing is generic_code_data::table::TYPE 2, the c_peers arrive at region B, and the component crossing must be generic_code_data::table::TYPE 1.
							- If the odd_i_peer is even
								* if the immersion crossing is generic_code_data::table::TYPE 1, the c_peers arrive at region B and the component crossing must be generic_code_data::table::TYPE 1
								* if the immersion crossing is generic_code_data::table::TYPE 2, the c_peers arrive at region D and the component crossing must be generic_code_data::table::TYPE 2

						   if the even_i_peer is positive and the odd_i_peer is negative, then the even_i_peer (and consequently the even_c_peer) 
						   must be on the left side of the above diagrams and the odd_i_peer (and consequently the odd_c_peer) must be on the 
						   right side).  Thus:
							- if the even_i_peer is odd
								* if the immersion crossing is generic_code_data::table::TYPE 1, the c_peers arrive at region D, and the component crossing must be generic_code_data::table::TYPE 1.
								* if the immersion crossing is generic_code_data::table::TYPE 2, the c_peers arrive at region B, and the component crossing must be generic_code_data::table::TYPE 2.
							- If the even_i_peer is even
								* if the immersion crossing is generic_code_data::table::TYPE 1, the c_peers arrive at region B and the component crossing must be generic_code_data::table::TYPE 2
								* if the immersion crossing is generic_code_data::table::TYPE 2, the c_peers arrive at region D and the component crossing must be generic_code_data::table::TYPE 1
						*/
						if (odd_i_peer >= 0 && even_i_peer >= 0)
						{
							if ((odd_i_peer % 2 == 1 && code_table[generic_code_data::table::TYPE][marker] == generic_code_data::TYPE1) ||
								(odd_i_peer % 2 == 0 && code_table[generic_code_data::table::TYPE][marker] == generic_code_data::TYPE2))
								component_code_table[generic_code_data::table::TYPE][even_c_peer/2] = generic_code_data::TYPE1;
							else
								component_code_table[generic_code_data::table::TYPE][even_c_peer/2] = generic_code_data::TYPE2;
						}
						else if (odd_i_peer < 0 && even_i_peer < 0)
						{
							if ((abs(odd_i_peer) % 2 == 1 && code_table[generic_code_data::table::TYPE][marker] == generic_code_data::TYPE1) ||
								(abs(odd_i_peer) % 2 == 0 && code_table[generic_code_data::table::TYPE][marker] == generic_code_data::TYPE2))
								component_code_table[generic_code_data::table::TYPE][even_c_peer/2] = generic_code_data::TYPE2;
							else
								component_code_table[generic_code_data::table::TYPE][even_c_peer/2] = generic_code_data::TYPE1;
						}
						else if (odd_i_peer >= 0 && even_i_peer < 0)
						{
							if ((odd_i_peer % 2 == 1 && code_table[generic_code_data::table::TYPE][marker] == generic_code_data::TYPE1) ||
								(odd_i_peer % 2 == 0 && code_table[generic_code_data::table::TYPE][marker] == generic_code_data::TYPE2))
								component_code_table[generic_code_data::table::TYPE][even_c_peer/2] = generic_code_data::TYPE2;
							else
								component_code_table[generic_code_data::table::TYPE][even_c_peer/2] = generic_code_data::TYPE1;
						}														
						else if (odd_i_peer < 0 && even_i_peer >= 0)
						{
							if ((even_i_peer % 2 == 1 && code_table[generic_code_data::table::TYPE][marker] == generic_code_data::TYPE1) ||
								(even_i_peer % 2 == 0 && code_table[generic_code_data::table::TYPE][marker] == generic_code_data::TYPE2))
								component_code_table[generic_code_data::table::TYPE][even_c_peer/2] = generic_code_data::TYPE1;
							else
								component_code_table[generic_code_data::table::TYPE][even_c_peer/2] = generic_code_data::TYPE2;
						}
						
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:       odd_c_peer = " << odd_c_peer << ", even_c_peer = " << even_c_peer 
	      << ", odd_i_peer = " << odd_i_peer << ", even_i_peer = " << even_i_peer << endl;
	debug << "bracket_polynomial:       component crossing = " << even_c_peer/2 << endl;
	debug << "bracket_polynomial:         generic_code_data::table::OPEER = " << component_code_table[generic_code_data::table::OPEER][even_c_peer/2] << endl;
	debug << "bracket_polynomial:         generic_code_data::table::COMPONENT = " << component_code_table[generic_code_data::table::COMPONENT][even_c_peer/2] << endl;
	debug << "bracket_polynomial:         generic_code_data::table::TYPE = " << component_code_table[generic_code_data::table::TYPE][even_c_peer/2] << endl;
}							
					}
				
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:   write generic_code_data::table::LABEL: " << endl;				
	
					/* now we have the components assigned, we can set the label*/
					for (int j=0; j< num_cpt_crossings; j++)
					{
						/* identify the marker crossing in immersion_crossing_map for graphical component crossing k */
						int marker = find_value(marker_crossing_map,j)->first;
						int odd_c_peer = component_code_table[generic_code_data::table::ODD_TERMINATING][j];
						int even_c_peer = component_code_table[generic_code_data::table::EVEN_TERMINATING][j];
						if (even_c_peer%2)
							swap(odd_c_peer,even_c_peer);
							
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:     graphical component crossing " << even_c_peer/2 << " immersion crossing marker = " << marker << endl;
						if (code_table[generic_code_data::table::LABEL][marker] == generic_code_data::VIRTUAL)
						{
							component_code_table[generic_code_data::table::LABEL][even_c_peer/2] = generic_code_data::VIRTUAL;
						}
						else if (component_head_semi_arc != -1)  // this graphical component is a flat knotoid
						{
							if (odd_c_peer >= component_head_semi_arc && component_code_table[generic_code_data::table::COMPONENT][(odd_c_peer-1)/2] == 0)
								component_code_table[generic_code_data::table::LABEL][even_c_peer/2] = generic_code_data::POSITIVE;
							else if (even_c_peer >= component_head_semi_arc && component_code_table[generic_code_data::table::COMPONENT][even_c_peer/2] == 0)
								component_code_table[generic_code_data::table::LABEL][even_c_peer/2] = generic_code_data::NEGATIVE;
							else
								component_code_table[generic_code_data::table::LABEL][even_c_peer/2] = generic_code_data::FLAT;
						}								
						else 
						{
							component_code_table[generic_code_data::table::LABEL][even_c_peer/2] = generic_code_data::FLAT;
						}
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:       odd_c_peer = " << odd_c_peer << ", even_c_peer = " << even_c_peer << endl;
	debug << "bracket_polynomial:       generic_code_data::table::LABEL = " << component_code_table[generic_code_data::table::LABEL][even_c_peer/2] << endl;
}	
					}
					
					/* Now we can create a generic_code_data object for the component, write it to an ostringstream and read it back,
					   start by re-writing the final ODD_TERMINATING and EVEN_TERMINATING crosing data into the component_code_table 
					*/										
					for (int j=0; j< num_cpt_crossings; j++)
					{
						component_code_table[generic_code_data::table::ODD_TERMINATING][j] = component_code_table[generic_code_data::table::OPEER][j];
						component_code_table[generic_code_data::table::EVEN_TERMINATING][j] = 2*j;
					}

					generic_code_data component_peer_code;
					component_peer_code.type = generic_code_data::peer_code;
					component_peer_code.num_crossings = num_cpt_crossings;
					component_peer_code.num_components = graphical_cpt_u_component_count[i];
					component_peer_code.code_table = component_code_table;
					component_peer_code.num_component_edges = vector<int>(component_peer_code.num_components);						
					
					for (int j=0; j < component_peer_code.num_components-1; j++)
						component_peer_code.num_component_edges[j] = first_edge_on_g_component[j+1] - first_edge_on_g_component[j];						
					
					component_peer_code.num_component_edges[component_peer_code.num_components-1] = 2*num_cpt_crossings - first_edge_on_g_component[component_peer_code.num_components-1];						
					
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   graphical component " << i << " initial code_table" << endl;
    print(component_code_table, debug, 4, "bracket_polynomial:   ");
}				

					ostringstream oss;
					write_peer_code(oss,component_peer_code);

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial: component code_data ";
	write_peer_code(debug,component_peer_code,true);
	debug << endl;
    print_code_data(debug,component_peer_code,"bracket_polynomial:   ");
}				
					
					read_peer_code (component_peer_code, oss.str());

					
					/* We can't set the head until we have the generic_code_data::table::EPEER row of the code table, which is completed by read_peer_code,
					   so have to set both the head and the immersion character manually here.
					*/
					if (component_head_semi_arc != -1)
					{
						if (component_head_semi_arc %2 == 0)
							component_peer_code.head = component_head_semi_arc/2;
						else
							component_peer_code.head = component_peer_code.code_table[generic_code_data::table::EPEER][(component_head_semi_arc-1)/2]/2;							
							
						component_peer_code.immersion_character = generic_code_data::character::PURE_KNOTOID;
						
						/* call valid_knotoid_data to set shortcut_crossing */
						valid_knotoid_input(component_peer_code);
					}
					else if (component_includes_unicursal_cpt_zero)
					{
						/* if we're dealing with a pure knotoid but have not set component_head_semi_arc then there is no edge   
						   within the shortcut terminating at a self-crossing of graphical component zero, meaning that  the 
						   graphical component is a knot-type knotoid rather than a pure knotoid.
						*/
						if (pure_knotoid_code_data)
						{
							component_peer_code.immersion_character = generic_code_data::character::KNOTOID;
							
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:   knotoid: graphical component 0 is a knot-type knotoid not a pure knotoid" << endl;
						}
						else
						{
							component_peer_code.immersion_character = code_data.immersion_character; // needed for knot-type knotoids and long knot
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:   knotoid: graphical component 0 has the same immersion character as original code data" << endl;
						}
					}

					if (variant == PARITY_ARROW_VARIANT)
					{
						/* if the knotoid_head_zig_zag_count is non-zero, write it into the corresponding location of 
						component_zig_zag_count
						*/
						if (knotoid_head_zig_zag_count != 0)
						{
							component_peer_code.head_zig_zag_count = knotoid_head_zig_zag_count;
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:   set component_peer_code.head_zig_zag_count = " << component_peer_code.head_zig_zag_count << endl;
							
							if (component_peer_code.immersion_character == generic_code_data::character::PURE_KNOTOID)
							{
								component_zig_zag_count[(component_head_semi_arc %2?1:0)][component_peer_code.head] = knotoid_head_zig_zag_count;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   set " << (component_head_semi_arc %2?"odd":"even") << " terminating zig-zag count for component_head_semi_arc "
          << component_head_semi_arc << " at crossing " << component_peer_code.head << " to " << knotoid_head_zig_zag_count << endl;
}
							}
						}
							
						bool zig_zags_present = false;
						for (int j=0; j<num_cpt_crossings; j++)
						{
							if (component_zig_zag_count[0][j] != 0 || component_zig_zag_count[1][j] !=0)
							{
								zig_zags_present = true;
								break;
							}
						}
						
						if (zig_zags_present || knotoid_head_zig_zag_count != 0)
						{
							component_peer_code.zig_zag_count = component_zig_zag_count;
						}					
					}
					
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   graphical component " << i << " component_code_data" << endl;
    print_code_data(debug, component_peer_code, "bracket_polynomial:   ");
	debug << "bracket_polynomial:   graphical component " << i << " peer code: ";
	write_peer_code(debug,component_peer_code, true); // zig_zags = true
	debug << endl;
}				
											
					/* reduce the graphical component by removing virtual Reidemeiseter I and both graphical (flat) and virtual 
					   Reidemeister II configurations.  Reducing the component_peer_code may remove linked unknoted components, 
					   so we add them to num_non_graphical_u_cpts and track those components that have been removed in 
					   unicursal_component_flags, which we initialize from unicursal_cpt_to_g_component_map and then bleach to
					   identify the components in the i-th graphical component.
					*/
					vector<int> unicursal_component_flags = unicursal_cpt_to_g_component_map;
					for (int j=0; j < num_unicursal_components; j++)
					{
						if (unicursal_component_flags[j] == i)
							unicursal_component_flags[j] = 1;
						else
							unicursal_component_flags[j] = 0;
					}

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   initial unicursal_component_flags: ";
	for (int j=0; j< num_unicursal_components; j++)
		debug << unicursal_component_flags[j] << ' ';
	debug << endl;
}											
					num_non_graphical_u_cpts += remove_virtual_Reidemeister_I(component_peer_code, unicursal_component_flags);

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   removing virtual Reidemeister I configurations updates num_non_graphical_u_cpts to " << num_non_graphical_u_cpts << endl;
	debug << "bracket_polynomial:   unicursal_component_flags: ";
	for (int j=0; j< num_unicursal_components; j++)
		debug << unicursal_component_flags[j] << ' ';
	debug << endl;
}						
					
					if (component_peer_code.num_crossings > 1)
					{
						num_non_graphical_u_cpts += remove_Reidemeister_II(component_peer_code, unicursal_component_flags);
						
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   removing Reidemeister II configurations updates num_non_graphical_u_cpts to " << num_non_graphical_u_cpts << endl;
	debug << "bracket_polynomial:   unicursal_component_flags: ";
	for (int j=0; j< num_unicursal_components; j++)
		debug << unicursal_component_flags[j] << ' ';
	debug << endl;
}						
					}
					
					if(component_peer_code.num_crossings > 0)
					{
						num_non_graphical_u_cpts +=  remove_virtual_components(component_peer_code, unicursal_component_flags);
						
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   removing floating virtual components updates num_non_graphical_u_cpts to " << num_non_graphical_u_cpts << endl;
	debug << "bracket_polynomial:   unicursal_component_flags: ";
	for (int j=0; j< num_unicursal_components; j++)
		debug << unicursal_component_flags[j] << ' ';
	debug << endl;
}						
					}

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   final unicursal_component_flags: ";
	for (int j=0; j< num_unicursal_components; j++)
		debug << unicursal_component_flags[j] << ' ';
	debug << endl;
}											
					/* If we are evaluating the parity arrow polynomial and some components have been removed, we need to determine
					   whether any of those components contain irreducible cusps and update the arrow_factor accordingly.
					   we count the number of components containing irreducible cusps in num_zig_zag_components, so we can add the 
					   appropriate number of delta terms for this state.
					   
					*/
					if (variant == PARITY_ARROW_VARIANT)
					{
						for (int j=0; j< num_unicursal_components; j++)					
						{
							if (unicursal_cpt_to_g_component_map[j] == i && unicursal_component_flags[j] == 0)
							{
								/* the jth unicursal component originally belonged to the ith graphical component but has been removed */								
								int zig_zag_count = unicursal_component_zig_zag_count[j];
								
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:   removed component " << j << ", zig-zag count = " << zig_zag_count << endl;
	
								if (component_peer_code.immersion_character != generic_code_data::character::CLOSED && j == 0)
								{
									/* The segment component in a knotoid contributes an L variable */
									if (zig_zag_count)
									{
										num_zig_zag_components++;  
										char z = 'a'+zig_zag_count/2-1;
										arrow_factor *= polynomial<int,bracket_variable>(string(1,z));
									}
									else if (!braid_control::EXPANDED_BRACKET_POLYNOMIAL)
									{
										num_zig_zag_components++;  // not a zig-zag but we do want to use L rather than D for this component
										char z = 'L'; // long segment with no zig-zags
										arrow_factor *= polynomial<int,bracket_variable>(string(1,z));
									}
								}
								else
								{
									/* Non-segment components contribute a K variable */
									if (zig_zag_count)
									{
										num_zig_zag_components++;  
										char z = 'n'+zig_zag_count/2-1;
										arrow_factor *= polynomial<int,bracket_variable>(string(1,z));
									}
								}		
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:   arrow_factor updated to " << arrow_factor << endl;
								
							}
						}
					}
					
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:   num_zig_zag_components updated to " << num_zig_zag_components << endl;
	
					if (component_peer_code.num_crossings > 0)
					{
						/* The component_peer_code might represent a disconnected diagram, in which case it represents multiple 
						   irreducible graphical components.  We successively isolate the first connected component by writing the 
						   subsequent components code data to subsequent_code_data, adjusting the edge numbering accordingly.
						   
						   If the component_peer_code head is not -1, the component graphs is a flat knotoid, in which case the 
						   segment component must be part of the first connected component, since the leg semi_arc is numbered 
						   zero.  That means the segment component is the first unicursal component and belongs to graphical
						   component zero.  That means the head is correctly in component_peer_code and that the default initialization
						   of the head in subsequent_code_data is also correct.
						*/					
						bool subsequent_connected_component;
						do
						{
							/* take a copy of the current unicursal_component_flags to pass to partition_component_flags, so that 
							   we determine which unicursal components are included in the first connected graphical component
							*/
							vector<int> partition_component_flags = unicursal_component_flags;
							generic_code_data subsequent_code_data = partition_peer_code(component_peer_code,partition_component_flags);
								
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial: create bracket_variable from component_peer_code ";
	write_peer_code(debug,component_peer_code);
	debug << endl;
	debug << "bracket_polynomial: component_peer_code code_data" << endl;
	print_code_data(debug,component_peer_code,"bracket_polynomial:   ");
	debug << "bracket_polynomial:  partition_component_flags: ";
	for (int j=0; j< num_unicursal_components; j++)
		debug << partition_component_flags[j] << ' ';
	debug << endl;
}																	

							if (variant == PARITY_ARROW_VARIANT && (pure_knotoid_code_data || braid_control::RELAXED_PARITY))
							{
								/* clear the component_peer_code zig-zag counts and head_zig_zag_count and add the relaxed
								   unicursal component variables to parity_factor
								*/
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "bracket_polynomial:   evaluating relaxed parity arrow polynomial: unicursal_component_zig_zag_count: ";
	for (int i=0; i< num_unicursal_components; i++)
		debug << unicursal_component_zig_zag_count[i] << ' ';
	debug << endl;
}
								component_peer_code.zig_zag_count.clear();
								component_peer_code.head_zig_zag_count = 0;
								
								for (int j=0; j< num_unicursal_components; j++)					
								{
									if (partition_component_flags[j] == 1 && unicursal_component_zig_zag_count[j] != 0)
									{
										int zig_zag_count = unicursal_component_zig_zag_count[j];
								
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:    unicursal component " << j << " is part of component_peer_code, zig-zag count = " << zig_zag_count << endl;
	
										if (component_peer_code.immersion_character != generic_code_data::character::CLOSED && j == 0)
										{
											/* segment components contribute a M variable, we can re-use variable characters a 
											   onwards and n onwards here because we're adding to parity_factor not arrow_factor
											*/

/*											polynomial<int,bracket_variable> component_term("a");
											bracket_variable bv("M_"+to_string(zig_zag_count/2));
											component_term.set_varmap('a',bv);
											parity_factor *= component_term;
*/		
											char z = 'a'+zig_zag_count/2-1;
											relaxed_parity_factor *= polynomial<int,bracket_variable>(string(1,z));
										}
										else
										{
											/* Non-segment components contribute a J variable */
											
/*											polynomial<int,bracket_variable> component_term("a");
											bracket_variable bv("J_"+to_string(zig_zag_count/2));
											component_term.set_varmap('a',bv);
											parity_factor *= component_term;
*/											
											char z = 'n'+zig_zag_count/2-1;
											relaxed_parity_factor *= polynomial<int,bracket_variable>(string(1,z));
										}		
				
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:    relaxed_parity_factor updated to " << relaxed_parity_factor << endl;
								
									}
								}
							}

							/* add the reduced graphical component as a mapped polynomial variable to the parity_factor for this state */
							bracket_variable component_graph(component_peer_code);
													
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   component_graph " << endl;
	print_bracket_variable(component_graph,debug,"bracket_polynomial:     ");	
}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "bracket_polynomial:     component_graph code_data" << endl;
	print_code_data(debug,component_graph.code_data,"bracket_polynomial:       ");	
}
							
							polynomial<int,bracket_variable> component_poly("a");
							component_poly.set_varmap('a',component_graph);
														
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:   component_poly " << component_poly << endl;
							
							parity_factor *= component_poly;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "bracket_polynomial:   parity_factor updated to " << parity_factor <<endl;
								
							if (subsequent_code_data.num_crossings != 0)
							{
								component_peer_code = subsequent_code_data;
								subsequent_connected_component = true;
								
								/* clear the unicursal_component_flags for those components we've just dealt with 
								   in the first partition
								*/
								for (int j=0; j< num_unicursal_components; j++)
								{
									if (partition_component_flags[j] == 1)
										unicursal_component_flags[j] = 0;
								}
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   more connected irreducible graphical components to consider, component_peer_code set to ";
	write_peer_code(debug,component_peer_code);
	debug << endl;
}
							}
							else
							{
								subsequent_connected_component = false;
							}
						} while (subsequent_connected_component);						
					}
				}
			}
	
			/* build up the next term of the polynomial, first clear the ostringstream */		
			oss.clear();
			oss.str("");
			
			if (sigma !=0)
				oss << "A^" << sigma;
			else
				oss << "1";
	
			if (variant == TURAEV_VARIANT)
				oss << "u^" << algebraic_crossing_number;
	
			polynomial<int,bracket_variable> term(oss.str());			
			
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "bracket_polynomial:   term initially set to " << term << endl;
						
			if (variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT)
			{				
				int num_delta_terms = num_non_graphical_u_cpts+num_virtual_components;
				
				if (variant == PARITY_ARROW_VARIANT && !braid_control::ZIG_ZAG_DELTA) // the default
					num_delta_terms -= num_zig_zag_components;
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "bracket_polynomial:   num_delta_terms = " << num_delta_terms << endl;
				
				for (int i=0; i< num_delta_terms; i++)
				{
					if (braid_control::EXPANDED_BRACKET_POLYNOMIAL)
						term *= polynomial<int,bracket_variable>("-A^2-A^-2");	
					else
						term *= polynomial<int,bracket_variable>("D");	
				}
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "bracket_polynomial:   term = " << term << endl;
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "bracket_polynomial:   num_non_graphical_u_cpts = " << num_non_graphical_u_cpts << endl;
	debug << "bracket_polynomial:   num_virtual_components = " << num_virtual_components << endl;
}
			}
			else
			{
				/* For the Kauffman bracket and the Turaev extended bracket, the number of delta terms is one less 
				   than the number of unicursal components.  The same is true for the arrow polynomial if we are
				   following the old convention of including delta with zig-zag components (K_i and L_i variables).
				   However, the default for the arrow polynomial is for ZIG_ZAG_DELTA to be false, in which case
				   it is the number of non-zig-zag components.
				*/
				int num_delta_terms = num_unicursal_components-1;
				
				if (variant == ARROW_VARIANT && !braid_control::ZIG_ZAG_DELTA) // the default
					num_delta_terms -= num_zig_zag_components-1;
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "bracket_polynomial:   num_delta_terms = " << num_delta_terms << endl;
				
				for (int i=0; i< num_delta_terms; i++)
				{
					if (braid_control::EXPANDED_BRACKET_POLYNOMIAL)
						term *= polynomial<int, bracket_variable>("-A^2-A^-2");	
					else
						term *= polynomial<int, bracket_variable>("D");	
				}
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "bracket_polynomial:   term = " << term << endl;
			}
	
			if (variant == ARROW_VARIANT)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "bracket_polynomial:   term before arrow factor = " << term;
	
				term *= arrow_factor;			

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	bool map_variables = polynomial_control::SUBSTITUTE_MAPPED_VARIABLES;
	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = true;			
	debug << ", arrow_factor = " << arrow_factor << endl;
	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = map_variables;			
	debug << "bracket_polynomial:   final term = " << term << endl;
}
			}
			else if (variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "bracket_polynomial:   term before arrow_factor and parity_factor = " << term << endl;
				for (int i=0; i< arrow_factor.nv; i++)
				{
					char ch = arrow_factor.getvar(i);
					int index=0;
					if (ch != 'A' && ch != 'D' && ch < 'n')
					{
						index = ch-'a'+1;					
						bracket_variable bv("L_"+to_string(index));
						arrow_factor.set_varmap(ch,bv);
					}
					else if (ch != 'A' && ch != 'D')
					{
						index = ch-'n'+1;	
						bracket_variable bv("K_"+to_string(index));					
						arrow_factor.set_varmap(ch,bv);
					}
				}

				term *= arrow_factor; // == 1 in the case of PARITY_VARIANT
				term *= parity_factor;
				
				for (int i=0; i< relaxed_parity_factor.nv; i++)
				{
					char ch = relaxed_parity_factor.getvar(i);
					int index=0;
					if (ch < 'n')
					{
						index = ch-'a'+1;					
						bracket_variable bv("M_"+to_string(index));
						relaxed_parity_factor.set_varmap(ch,bv);
					}
					else 
					{
						index = ch-'n'+1;	
						bracket_variable bv("J_"+to_string(index));					
						relaxed_parity_factor.set_varmap(ch,bv);
					}
				}
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	bool map_variables = polynomial_control::SUBSTITUTE_MAPPED_VARIABLES;
	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = true;			
	
	debug << "bracket_polynomial:   arrow_factor = " << arrow_factor << ", parity_factor = " << parity_factor << ", relaxed_parity_factor = " << relaxed_parity_factor << endl;
	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = map_variables;			
}
				
				term *= relaxed_parity_factor; // == 1 unless we're doing the relaxed version of PARITY_ARROW

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	bool map_variables = polynomial_control::SUBSTITUTE_MAPPED_VARIABLES;
	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = false;			
	debug << "bracket_polynomial:   final term = " << term << endl;
	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = map_variables;			
}
			}

			bracket_poly += term;
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "bracket_polynomial:   updated bracket_poly = ";
	bool map_variables = polynomial_control::SUBSTITUTE_MAPPED_VARIABLES;
	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = false;
	debug << bracket_poly << endl;
	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = map_variables;
}

			
			/* increment state */
			if (num_state_smoothed_crossings > 0)
			{
				place = num_state_smoothed_crossings-1;
				do
				{
					if (state[place] == -1)
					{
						state[place] = 1;
						for (int i = place+1; i< num_state_smoothed_crossings; i++)
							state[i] = -1;
						break;
					}
					else if (place == 0)
						finished = true;
					else
						place--;
				}while (!finished);
			}
			else
			{
				finished = true;
			}
	
			if (!braid_control::SILENT_OPERATION && braid_control::WAIT_SWITCH)
			{
			    if (++state_count == braid_control::wait_threshold)
			    {
					cout << ".";
	   				cout.flush();
					
					if (braid_control::wait_count > 1000)
					{
						braid_control::reset_count++;
			    		cout << "\nworking bracket polynomial (" << braid_control::reset_count << ")\n";
			    		braid_control::wait_count = 0;
					}
					else
			    		braid_control::wait_count++;
			    		
					state_count = 0;
		    	}
			}
			
		} while (!finished);
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "bracket_polynomial: bracket_poly  before normalizing = " ;
	bool map_variables = polynomial_control::SUBSTITUTE_MAPPED_VARIABLES;
	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = false;
	debug << bracket_poly << endl;
	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = map_variables;
}		
		if (braid_control::NORMALIZE_BRACKET)
			bracket_poly *= normalizing_factor;
	}
	else
	{		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:  no classical non-shortcut crossings in code data to smooth, setting bracket_poly to 1" << endl;
		
		bracket_poly = polynomial<int,bracket_variable>(1);
	}
		
	if (braid_control::JONES_POLYNOMIAL)
	{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:  substituting A = t^{-1/4} in " << bracket_poly << endl;

		ostringstream oss;
		
		oss << bracket_poly;
		
		string kstring = oss.str();
		
		ostringstream jss;
			
		/* in the following loop, p is incremented within the body of the loop as appropriate */
		for (string::size_type p=0;p<kstring.length();) 
		{
			if (kstring[p] == 'A')
			{
				jss << "t";
				p++;
				
				if (kstring[p] != '^')
				{
					if (braid_control::TeX_POLYNOMIAL_OUTPUT)
						jss << "^{-1/4}";
					else				
						jss << "^-1/4";
				}
				
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:  incremented p to " << p << " kstring[p] = " << kstring[p] << endl;
			}
			else if (kstring[p] == '^')
			{
				p++;
				int exp;
				get_number(exp,kstring,p);
				
				if (exp != -4) // in which case the new exponent is 1
				{
					jss << "^";
					
					if (braid_control::TeX_POLYNOMIAL_OUTPUT)
						jss << '{';
					
					if (exp > 0)
						jss << '-';
	
					if (exp % 4 == 0)
						jss << abs(exp)/4;
					else if (exp % 2 == 0)
						jss << abs(exp)/2 << "/2";
					else
						jss << abs(exp) << "/4";
						
					if (braid_control::TeX_POLYNOMIAL_OUTPUT)
						jss << '}';

				}
					
				if (kstring[p] == '-')
					p++;
					
				while (isdigit(kstring[p]))
					p++;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:  idntified exp = " << exp << ", incremented p to " << p << " kstring[p] = " << kstring[p] << endl;

			}
			else
			{
				jss << kstring[p];
				p++;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "bracket_polynomial:  incremented p to " << p << " kstring[p] = " << kstring[p] << endl;
			}

		}


if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "bracket_polynomial: Jones polynomial = " << jss.str() << endl;

		if (!braid_control::SILENT_OPERATION)
			cout<< "\n\nJones polynomial = " << jss.str() << endl;
	
		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
	
			output << "Jones polynomial = ";
			if (braid_control::OUTPUT_AS_INPUT)
				output << '\n';
		}
		output << jss.str() << endl;
	}
	else
	{
		if (variant == ARROW_VARIANT)
		{
			for (int i=0; i< bracket_poly.nv; i++)
			{
				char ch = bracket_poly.getvar(i);
				int index=0;
				if (ch != 'A' && ch != 'D' &&ch < 'n')
				{
					index = ch-'a'+1;					
					bracket_variable bv("L_"+to_string(index));
					bracket_poly.set_varmap(ch,bv);
				}
				else if (ch != 'A' && ch != 'D')
				{
					index = ch-'n'+1;	
					bracket_variable bv("K_"+to_string(index));					
					bracket_poly.set_varmap(ch,bv);
				}
			}
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	if (variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT)
	{
		bool map_variables = polynomial_control::SUBSTITUTE_MAPPED_VARIABLES;
		polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = false;
		debug << "bracket_polynomial: bracket_polynomial = " << bracket_poly << endl;
		polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = map_variables;
	}
	else
		debug << "bracket_polynomial: bracket_polynomial = " << bracket_poly << endl;
}

		bool TeX = polynomial_control::TeX;
		
		if (braid_control::TeX_POLYNOMIAL_OUTPUT)
			polynomial_control::TeX = true;

		if (!braid_control::SILENT_OPERATION)
		{
			if (braid_control::NORMALIZE_BRACKET)
				cout << "\n\nnormalized ";
			
			switch (variant)
			{
				case KAUFFMAN_VARIANT: cout << "Kauffman bracket"; break;
				case JONES_VARIANT: cout << "Jones"; break;
				case TURAEV_VARIANT: cout << "Turaev extended bracket"; break;
				case ARROW_VARIANT: cout << "arrow"; break;
				case PARITY_VARIANT: cout << "parity bracket"; break;
				case PARITY_ARROW_VARIANT: cout << "parity arrow"; break;
				default: cout << "unknown";
			}	
				
			if (variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT)
			{
				bool map_variables = polynomial_control::SUBSTITUTE_MAPPED_VARIABLES;
				polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = false;			
				cout << " polynomial = " << bracket_poly << endl;
				polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = map_variables;
			}
			else
				cout<< " polynomial = " << bracket_poly << endl;			
		}
	
		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			
			if (braid_control::NORMALIZE_BRACKET)
				output << "normalized ";
	
			switch (variant)
			{
				case KAUFFMAN_VARIANT: output << "Kauffman bracket"; break;
				case JONES_VARIANT: output << "Jones"; break;
				case TURAEV_VARIANT: output << "Turaev extended bracket"; break;
				case ARROW_VARIANT: output << "arrow"; break;
				case PARITY_VARIANT: output << "parity bracket"; break;
				case PARITY_ARROW_VARIANT: output << "parity arrow"; break;
				default: output << "unknown";
			}	
	
			output << " polynomial = ";
			if (braid_control::OUTPUT_AS_INPUT)
				output << '\n';								
		}
		
		if (variant == PARITY_VARIANT || variant == PARITY_ARROW_VARIANT)
		{
			bool map_variables = polynomial_control::SUBSTITUTE_MAPPED_VARIABLES;
			polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = false;
			output << bracket_poly << endl;
			polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = map_variables;
		}
		else
			output << bracket_poly << endl;
			
		polynomial_control::TeX = TeX;
		
	}	
}

/* add_zig_zag_polynomial_factor takes a vector indicating the unicursal components removed by remove_virtual_Reidemeister_I, remove_Reidemeister_II or
   remove_virtual_component and checks zig_zag_counts to see if there are any irreducible cusps associated with those components.  If there
   are, it updates the polynomial "factor" with additional multiplicative terms related to those components and returns the number of components
   for which it has added a term (i.e the number of removed components containing irreducible cusps).  
   
   Those components not containing irreducible cusps contribute a delta term, which is handled separately.
*/
int add_zig_zag_polynomial_factor(polynomial<int,bracket_variable>& factor, bool closed_immersion, vector<int>& removed_components, vector<int>& zig_zag_counts)
{
	int component_count = 0;
	
	for (unsigned int i=0; i< removed_components.size(); i++)
	{
		int component = removed_components[i];
		int zig_zag_count = zig_zag_counts[component];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "add_zig_zag_polynomial_factor: removed component " << component << ", zig-zag count = " << zig_zag_count << endl;
	
		if (!closed_immersion && component == 0)
		{
			/* The segment component in a knotoid contributes an L variable */
			if (zig_zag_count)
			{
				component_count++;  
				char z = 'a'+zig_zag_count/2-1;
				factor *= polynomial<int,bracket_variable>(string(1,z));
			}
//			else if (long_knot && !braid_control::EXPANDED_BRACKET_POLYNOMIAL)
			else if (!braid_control::EXPANDED_BRACKET_POLYNOMIAL)
			{
				component_count++;  // not a zig-zag but we do want to use L rather than D for this component
				char z = 'L'; // long segment with no zig-zags
				factor *= polynomial<int,bracket_variable>(string(1,z));
			}
		}
		else
		{
			/* Non-segment components contribute a K variable */
			if (zig_zag_count)
			{
				component_count++;  
				char z = 'n'+zig_zag_count/2-1;
				factor *= polynomial<int,bracket_variable>(string(1,z));
			}
		}		
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "add_zig_zag_polynomial_factor: factor updated to " << factor << endl;
	
	}
	
	return component_count;
}

					
