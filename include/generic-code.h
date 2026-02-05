/********************************************************************************************************************************************

   generic_code_data stores information relating to a peer code or a Gauss code.  It can handle also handle the deprecated immersion code and 
   supports Gauss codes described as Planar Diagram data.  The code is capable of describing knots and links, and knotoids.
   
   The original version of generic_code_data supported only knotoids and multiknotoids involving a single open segment but in April 2023 support 
   was introduced for multi-linkoids that supportmultiple open segments.  As in the original version, open segments are stored along with a 
   shortcut that passes everywhere under the diagram.  With only one open component the requirement was for the segment component to be numbered 
   first and the head crossing was indicated with a ^ character.  The head crossing was set to -1 when the generic code had no head crossing and 
   head was consequently used as a flag to indicate pure knotoids, or as a check that the code data was consistent with the immersion character
   had been set correctly.
   
   In the updated version, the extended peer code syntax  ^: open component with leg the first edge  @: knot-like open component $: open component 
   with leg the last edge we introduced num_open_components, initialized to zero, which replaced head for the above checks and recorded the head
   crossing for each non-knot-like open component in vector<int> head.  The first component should be open and numbering should start from the leg, 
   so that the new code is as close to the legacy format as possible.
   
   
   TO DO
   =====
   
   int head_zig_zag_count wil need migrating to vector<int> head_zig_zag_count if parity arrow supports multi-linkoids.
   
   
   
   code table rows:
   
   OPEER: for immersion and peer codes, at position i, the (odd) peer of the even edge 2i, for Gauss codes the 
          terminating edge on the first visit to the crossing
          
   EPEER: for immersion and peer codes, at position i, the (even) peer of the odd edge 2i+1, for Gauss codes the 
          terminating edge on the second visit to the crossing
   
   TYPE: at position i, the type of crossing i; type 1 crossings have the even edge (or second visit) as the lower
         terminating edge in the standard (braid) view of a crossing.
		 
   LABEL: at position i, whether the naming edge (or second visit) passes over, under, vitually or flat through the crossing
   
   EVEN_TERMINATING: at position i, for immersion and peer codes, the even terminating edge at crossing i, for 
                     Gauss codes the terminating edge on the second visit to the crossing
   
   ODD_TERMINATING: at position i, for immersion and peer codes, the odd terminating edge at crossing i, for 
                     Gauss codes the terminating edge on the first visit to the crossing
   
   EVEN_ORIGINATING: at position i, the edge following the ODD_TERMINATING edge at the crossing

   ODD_ORIGINATING: at position i, the edge following the EVEN_TERMINATING edge at the crossing
    
   COMPONENT: at position i, for immersion and peer codes, the component to which the naming edge at crossing i belongs, 
              for Gauss codes, the component to which the first visit to the crossing belongs
              
   For the parity arrow polynomial, generic_code_data includes the number of irreducible cusps on each edge    in zig_zag_count.  This matrix has two 
   rows and num_crossings columns (corresponding to the crossings in the immersion), row 0 (even) stores the number of cusps on the even terminating 
   edge of the crossing and row 1 (odd) stores the number of cusps on the odd terminating edge.  If the first cusp on the edge is a left turn the count 
   is recorded as a negative number and if the first cusp is a right turn it is recorded as a positive number.

***********************************************************************************************************************************************/

#include <list>
#include <vector>

class component_character
{
public:

	enum
	{
		CLOSED,
		PURE_START_LEG, // leg on even edge
		PURE_END_LEG,   // leg on odd edge
		KNOT_TYPE_START_LEG, // leg and hed on even edge
		KNOT_TYPE_END_LEG // leg and hed on odd edge
	};

	int type;
	int head_semi_arc;
	
	component_character():type(CLOSED),head_semi_arc(-1){}
};

class generic_code_data
{

public:

	enum code_type
	{
		unknown,
		immersion_code,
		peer_code,
		gauss_code,
	};
	
	enum label
	{
		NEGATIVE = -1, // naming arc is the under-arc
		VIRTUAL = 0,
		POSITIVE = 1, // naming arc is the over-arc
		FLAT = 2,
		LEFT = 3, // LEFT and RIGHT used for labelling flat crossings in Gauss codes
		RIGHT = 4,
		ODD = 5,
		SINGULAR = 6,
		VOID = 9 // used as a marker when removing virtual Reidemeister II detours and as a generic void value
	};

	enum gauss_mode 
	{
		UNDER = -1,
		OVER = 1
	};
	
	enum type 
	{
		UNKNOWN = 0,
		TYPE1 = -1,
		TYPE2 = 1
	};
	
    enum character 
    {
		CLOSED,
		LONG_KNOT,
		PURE_KNOTOID, // indicates a pure knotoid peer code
		KNOTOID,  // indicates knot-type knotoid peer codes and all knotoid Gauss codes
		MULTI_LINKOID  // with peer codes, requires use of the extended syntax @: knot-like open component $: open component with leg the last edge not the first
	};
	
	enum table // uses default emum values but included to be explicit and for protection
	{
		EPEER =				0,	
		OPEER =				1,	
		TYPE =  			2,
		LABEL =				3,
		EVEN_TERMINATING =	4,
		ODD_TERMINATING =	5,
		EVEN_ORIGINATING =	6,
		ODD_ORIGINATING =	7,
		COMPONENT =			8,	// the component to which the naming (even) edge belongs
		CODE_TABLE_SIZE =	9
	};
	
	int type; // code_type
	int head; // used for knotoids
	int num_open_components; 
	
	/* component_type is used for multi-linkoids.  If num_open_components > 0, component_type is a vector of num_components integers pairs that records the 
	   nature of the component and, if necessary, the first undercrossing in the shortcut corresponding to that component.
	   - which components are open and which are closed
	   - whether the component is knot-type open or pure
	   - whether the leg of the open component is even or odd 
	*/
	vector<component_character> component_type; 
	
	int head_zig_zag_count; // used with knotoids by parity arrow
	
    /* immersion records the character of the immersion described by the generic code data and is set to one 
	   the character enum's above.
	   
	   the immersion character LONG_KNOT is not used by switch polynomial invariant template functions, since 
	   they don't have access to generic code date.  Instead, they use braid_control::LONG_KNOT but for other 
	   functions, it it more natural to include the indication of long status in the generic code data.
    */
   	
	int immersion_character;  // a character enum value
	int num_crossings;
	int num_components;

	matrix<int> code_table;
	vector<int> num_component_edges;
	vector<int> first_edge_on_component;
	vector<int> term_crossing;
	vector<int> orig_crossing;
	vector<int> shortcut_crossing;  // used for knotoids
	matrix<int> zig_zag_count;  // used by parity arrow polynomial
	
	bool multi_virtual;
	vector<int> virtual_index; // used with multiple virtual knots


	generic_code_data(): type(unknown),num_open_components(0),head_zig_zag_count(0),immersion_character(character::CLOSED),num_crossings(0),num_components(0),multi_virtual(false) {}
	generic_code_data(int n, int c): type(unknown),head(-1),head_zig_zag_count(0),immersion_character(character::CLOSED),num_crossings(n),num_components(c),
	                                 code_table(matrix<int>(CODE_TABLE_SIZE,n)),num_component_edges(vector<int>(c)), 
									 first_edge_on_component(vector<int>(c)),term_crossing(vector<int>(2*n)), orig_crossing(vector<int>(2*n)), 
									 shortcut_crossing(vector<int>(n)), multi_virtual(false), virtual_index(vector<int>(n)){}
};

/* next the class structure for tracking the crossing labels when converting Gauss codes to labelled peer codes */

class gc_pc_xlabels
{
	public:
	
	pair<int,int> gauss;
	pair<int,int> immersion;
	
	enum move {X01, X02, X11, X12, X21, X22, X31, X32, X41, X42, X2OV1, X2OV2, X2UV1, X2UV2, V, I, WL, WR, C, END};
	enum direction {UNKNOWN,UPWARDS,DOWNWARDS,SIDEWAYS};
	
	gc_pc_xlabels(): gauss({-1,-1}), immersion({-1,-1}) {}
	
};

void read_immersion_code (generic_code_data& code_data, string input_string);
void write_immersion_code(ostream& s, generic_code_data& code_data);
void print_code_data(ostream& s, generic_code_data& code_data, string prefix="");
void read_peer_code (generic_code_data& code_data, string input_string);
void write_peer_code(ostream& s, const generic_code_data& code_data, bool zig_zags=false, bool labelled=true);
void read_gauss_code (generic_code_data& code_data, string input_string);
void write_gauss_code(ostream& s, generic_code_data& code_data, bool OU_FORMAT = false);
void read_code_data (generic_code_data& code_data, string input_string);
void write_code_data(ostream& s, generic_code_data& code_data);
void read_planar_diagram (generic_code_data& code_data, string input_string);
void write_planar_diagram(ostream& s, generic_code_data& code_data);
string convert_gauss_code(string OU_gauss_code);
void renumber_peer_code(generic_code_data& code_data, vector<int> shift);
int remove_peer_code_component(generic_code_data& code_data, int component, vector<int>& component_flags);
int remove_virtual_components(generic_code_data& code_data, vector<int>& component_flags);
generic_code_data partition_peer_code(generic_code_data& code_data, vector<int>& component_flags);
void cycle_gauss_code_labels(generic_code_data& gauss_code_data, int edge);
void cycle_gauss_code_labels(matrix<int>& crossing_peers, vector<int>& num_component_edges, vector<int>& first_edge_on_component, int num_crossings,int num_components, int edge);
list<int> vertex_span (generic_code_data& code_data, int initial_crossing, vector<int>* exclude=0);
bool calculate_turning_cycles(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles);
bool realizable_code_data(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles);
bool valid_knotoid_input(generic_code_data& code_data);
string over_preferred_gauss_code(generic_code_data& code_data, bool unoriented);
int amalgamate_zig_zag_counts(int a, int b);
string read_dowker_code (string input_string);
int three_connected(generic_code_data& code_data, bool flat_crossings);
bool valid_multi_linkoid_input(generic_code_data& code_data);
int code_data_writhe(generic_code_data& code_data);
