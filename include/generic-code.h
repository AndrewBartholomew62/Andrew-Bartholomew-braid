/* code table rows
   
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
              
   For the parity arrow polynomial, generic_code_data includes the number of irreducible cusps on each edge
   in zig_zag_count.  This matrix has two rows and num_crossings columns (corresponding to the crossings in the immersion),
   row 0 (even) stores the number of cusps on the even terminating edge of the crossing and row 1 (odd) stores the number of 
   cusps on the odd terminating edge.  If the first cusp on the edge is a left turn the count is recorded as a negative number
   and if the first cusp is a right turn it is recorded as a positive number.
              
*/

#define EPEER 				0	
#define OPEER 				1	
#define TYPE    			2
#define LABEL				3
#define EVEN_TERMINATING 	4
#define ODD_TERMINATING 	5
#define EVEN_ORIGINATING 	6
#define ODD_ORIGINATING 	7
#define COMPONENT 			8	// the component to which the naming (even) edge belongs
#define CODE_TABLE_SIZE 	9

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
		VOID = 9 // used as a marker when removing virtual Reidemeister II detours and as a generic void value
	};
	
	enum type 
	{
		TYPE1 = -1,
		TYPE2 = 1
	};
	
    enum character 
    {
		CLOSED,
		LONG_KNOT,
		PURE_KNOTOID, // indicates a pure knotoid peer code
		KNOTOID  // indicates knot-type knotoid peer codes and all knotoid Gauss codes
	};

	int type;
	int head; // used for knotoids
	int head_zig_zag_count; // used with knotoids by parity arrow
	
    /* immersion records the character of the immersion described by the generic code data and is set to one 
	   the character enum's above.
	   
	   the immersion character LONG_KNOT is not used by switch polynomial invariant template functions, since 
	   they don't have access to generic code date.  Instead, they use braid_control::LONG_KNOT but for other 
	   functions, it it more natural to include the indication of long status in the generic code data.
    */
   	
	int immersion;
	int num_crossings;
	int num_components;

	matrix<int> code_table;
	vector<int> num_component_edges;
	vector<int> first_edge_on_component;
	vector<int> term_crossing;
	vector<int> orig_crossing;
	vector<int> shortcut_crossing;  // used for knotoids
	matrix<int> zig_zag_count;  // used by parity arrow polynomial


	generic_code_data(): type(unknown),head(-1),head_zig_zag_count(0),immersion(character::CLOSED),num_crossings(0),num_components(0) {}
	generic_code_data(int n, int c): type(unknown),head(-1),head_zig_zag_count(0),immersion(character::CLOSED),num_crossings(n),num_components(c),
	                                 code_table(matrix<int>(CODE_TABLE_SIZE,n)),num_component_edges(vector<int>(c)), 
									 first_edge_on_component(vector<int>(c)),term_crossing(vector<int>(2*n)), orig_crossing(vector<int>(2*n)), 
									 shortcut_crossing(vector<int>(n)){}
};
