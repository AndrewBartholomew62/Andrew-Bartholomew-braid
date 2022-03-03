/*************************************************************************
                      gauss-orientation.h

***************************************************************************/

/* The internal representation of gauss_orientation_data numbers crossings from zero.
   However, when writing and reading gauss_orientation data as gauss codes, the crossings are
   numbered from 1.

   By default Gauss orientation data ignores both virtual and knotoid shortcut crossings.  If the extended version 
   is being used, as indicated by the boolean extended in the constructor below, these crossings are included 
   but their presence is indicated by setting their crossing number in the orientation_matrix to be negative.
   
   For the parity arrow polynomial, gauss_orientation_data includes the number of irreducible cusps on each arc
   in zig_zag_count.  This matrix has two rows and num_terms/2 columns (corresponding to the crossings in the code),
   row 0 (even) stores the number of cusps on the arc preceeding the second (even) visit to the crossing and row 1 (odd)
   stores the number of cusps on the arc preceeding the first (odd) visit.
   
   Note that when we write a gaus_orientation_data object containing zig_zag_count data to an ostream, we write the
   zig-zag counts in the order corresponding to the terms of the orientation matrix and not indexed by the crossing
   numbers.
   
   In the current implementation, zig_zag_count data is only added when constructing gauss_orientation_data from a 
   labelled peer code.  The gauss_orientation_data records the absolute value of the count, not the sign recorded in
   generic_code_data.
*/    
class gauss_orientation_data
{
public:	
	int 			num_terms; // twice the number of classical crossings
	int 			num_components;
	bool 			extended_version; // identifies whether virtual and knotoid shortcut crossings are included
	vector<int>		num_terms_in_component; 
	vector<int>		start_of_component; //offset into orientation matrix for the start of each component
	vector<int>		immersion_crossing; // maps Gauss data crossings to immersion crossings
	matrix<int> 	zig_zag_count;  // used by parity arrow polynomial 2 x num_terms/2 matrix
	matrix<int>		orientation_matrix; // columns of row[0] = LEFT|RIGHT row[1] = <crossing-number>
	
	enum {LEFT = 3, RIGHT = 4};  // only need LEFT < RIGHT; values align with generic_code_data, just because...
	enum parity {NONE = 0, ODD = 1, EVEN = 2};  

	gauss_orientation_data (): num_terms(0), num_components(0), extended_version(false) {}
	gauss_orientation_data (generic_code_data& code_data,bool extended=false);
	int terms() {return(num_terms);}
	bool operator == (gauss_orientation_data&) const;	
	void reverse_all ();
	void reverse_components( vector<bool>&);
};

void print_gauss_data(gauss_orientation_data g, ostream& os, string prefix);
void write_gauss_data(gauss_orientation_data g, ostream& os, bool zig_zags=false, int immersion_character = generic_code_data::character::CLOSED, int head_zig_zag_count = 0);
gauss_orientation_data left_preferred(const gauss_orientation_data& g, bool unoriented=false, int immersion_character=generic_code_data::character::CLOSED);
