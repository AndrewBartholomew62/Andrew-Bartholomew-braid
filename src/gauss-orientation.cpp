/************************************************************************

Gauss orientation data records the orientation of non-virtual crossings in a diagram in the manner used
used by Gauss codes of doodles or flat knots or links (e.g "L1 L2 R1 R2 L3 L4 R3 R4 / # # # #").  The Gauss 
orientation data of a knotoid is determined by ignoring shortcut crossings as well as virtual crossings.

Each component of the diagram has an orientation associated with it determined by the semi-arc labels and 
we trace the diagram respecting this orientation.  As we traverse a crossing the Gauss orientation data 
records whether the other strand crosses our path from left to right or from right to left, according to
the orientation.

The orientation data is stored as an orientation_matrix comprising two rows and having two columns for each 
non-virtual crossing in the diagram.  The columns of the orientation matrix are of the form 

  row[0] = LEFT|RIGHT row[1] = <crossing-number>
 
and are numbered according to the order in which the non-virtual crossings are encountered as the diagram 
is traced.

gauss_orientation_data::gauss_orientation_data (generic_code_data& code_data)
gauss_orientation_data::gauss_orientation_data (k_dets& k)
void gauss_orientation_data::reverse_all()
void gauss_orientation_data::reverse_components( vector<bool>&)
void print_gauss_data(gauss_orientation_data g, ostream& os, string prefix)
void write_gauss_data(gauss_orientation_data g, ostream& os, bool zig_zags)
gauss_orientation_data left_preferred(gauss_orientation_data g, bool unoriented=false, int immersion_character = generic_code_data::character::CLOSED)
bool gauss_orientation_data::operator == (gauss_orientation_data& b) const

**************************************************************************/
using namespace std;
#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>
#include <iomanip>
#include <algorithm>


/********************* External variables ************************/
extern ofstream     debug;


#include <util.h>
#include <quaternion-scalar.h>
#include <polynomial.h>
#include <matrix.h>
#include <generic-code.h>
#include <debug-control.h>
#include <gauss-orientation.h>

bool operator < (const matrix<int>& a, const matrix<int>& b)
{
	if (a.numrows() != b.numrows() || a.numcols() != b.numcols())
		return false;
		
	for(unsigned int i= 0; i < a.numrows(); i++)
	for(unsigned int j= 0; j < a.numcols(); j++)
	{
		if (a[i][j] < b[i][j])
			return true;
		else if (a[i][j] > b[i][j])
			return false;
	}
	
	return false;
}

/* Gauss orientation data records the left or right polarity of classical crossings in a diagram.  By default it ignores both
   virtual and knotoid shortcut crossings.  If the extended_version is being used, as indicated by the boolean extended, these
   crossings are included but their presence is indicated by setting their crossing number in the Gauss data to be negative.
*/
gauss_orientation_data::gauss_orientation_data (generic_code_data& code_data, bool extended)
{
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_orientation_data::constructor(generic_code_data): code_data: " << endl;
	print_code_data(debug,code_data,"gauss_orientation_data::constructor(generic_code_data):   ");
}
	vector<int>& term_crossing = code_data.term_crossing;

	if (code_data.type == generic_code_data::gauss_code) 
	{
		num_terms = 2 * code_data.num_crossings;
		num_components = code_data.num_components;
		extended_version = false;
		
		num_terms_in_component = vector<int>(num_components);
		for (int i=0; i< num_components; i++)
			num_terms_in_component[i] = code_data.num_component_edges[i];
			
		start_of_component = vector<int>(num_components);
		for (int i=1; i< num_components; i++)
			start_of_component[i] = start_of_component[i-1] + num_terms_in_component[i-1];
		
		immersion_crossing = vector<int>(code_data.num_crossings);
		for (int i=0; i< code_data.num_crossings; i++)
			immersion_crossing[i] = i;

		orientation_matrix = matrix<int>(2,num_terms);
		classical_gauss_data = vector<int>(num_terms);
		classical_crossing_sign = vector<int>(code_data.num_crossings);
		
		matrix<int>& code_table = code_data.code_table;

		int index = 0;
		
		for (int i=0; i <  num_components; i++)
		{
			int start = code_data.first_edge_on_component[i];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_orientation_data::constructor(generic_code_data): start_of_component = " << start << endl;

			int edge=start;
				
			/*	trace this component */
			do
			{		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_orientation_data::constructor(generic_code_data): edge = " << edge;

				int next_crossing = term_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", next_crossing = " << next_crossing;

				bool first_visit = (code_table[OPEER][next_crossing] == edge? true: false);
				int next_crossing_orientation_polarity;
				
	
				if ((first_visit && code_table[TYPE][next_crossing] == generic_code_data::TYPE1) ||
					(!first_visit && code_table[TYPE][next_crossing] == generic_code_data::TYPE2))						
					next_crossing_orientation_polarity = gauss_orientation_data::RIGHT;
				else
					next_crossing_orientation_polarity = gauss_orientation_data::LEFT;	
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", polarity = " << (next_crossing_orientation_polarity == gauss_orientation_data::LEFT? "LEFT" : "RIGHT") << endl;
				
				orientation_matrix[0][index] = next_crossing_orientation_polarity;
				orientation_matrix[1][index] = next_crossing+1;

				if ((first_visit && code_table[LABEL][next_crossing] == generic_code_data::POSITIVE) ||
					(!first_visit && code_table[LABEL][next_crossing] == generic_code_data::NEGATIVE))						
					classical_gauss_data[index] = (next_crossing+1)*-1;
				else
					classical_gauss_data[index] = next_crossing+1;	
															
				index++;

				if (first_visit)
					edge = code_table[EVEN_ORIGINATING][next_crossing];
				else
					edge = code_table[ODD_ORIGINATING][next_crossing];
			} while (edge != start);				
		}
		
		for (int i=0; i< code_data.num_crossings; i++)
		{
			if ( (code_table[LABEL][i] == generic_code_data::POSITIVE && code_table[TYPE][i] == generic_code_data::TYPE2) ||
				 (code_table[LABEL][i] == generic_code_data::NEGATIVE && code_table[TYPE][i] == generic_code_data::TYPE1)
			   )
				classical_crossing_sign[i] = 1;
			else if ( (code_table[LABEL][i] == generic_code_data::POSITIVE && code_table[TYPE][i] == generic_code_data::TYPE1) ||
					  (code_table[LABEL][i] == generic_code_data::NEGATIVE && code_table[TYPE][i] == generic_code_data::TYPE2)
					)
				classical_crossing_sign[i] = -1;
		}	
			
	}
	else
	{
		extended_version = extended;
		matrix<int>& code_table = code_data.code_table;
	
		/* write Gauss code from immersion data in code_data */
		int num_edges = 2*code_data.num_crossings;
		vector<int> edge_flag(num_edges); // initialized to zero

		matrix<int>& peer_code_zig_zag_count = code_data.zig_zag_count;
		
		bool track_zig_zag_counts;
		if (code_data.zig_zag_count.numrows() !=0)
			track_zig_zag_counts = true;
		else
			track_zig_zag_counts = false;		
	
		/* we need to understand the mapping betwen immersion crossings and gauss crossings when we create the gauss code, 
		   since virtal crossings are ignored by the Gauss code.  We also want to record in the gauss_orientation_data the reverse
		   mapping between classical crossings and the immersion_crossing.
		   
		   If we are creating an extended version of Gauss orientation data, we consider all crossings to be "classical".
		*/
		int num_classical_crossings = code_data.num_crossings;
		bool pure_knotoid_code_data = false;
		vector<int>& shortcut_crossing = code_data.shortcut_crossing;
		
		if (extended_version)
		{
			if (code_data.immersion == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && shortcut_crossing.size())
				pure_knotoid_code_data = true;
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_orientation_data::constructor(generic_code_data): extended version: num_classical_crossings = " << num_classical_crossings << endl;
		}
		else
		{
			for (int i=0; i< code_data.num_crossings; i++)
			{
				if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
					num_classical_crossings--;
			}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_orientation_data::constructor(generic_code_data): num_classical_crossings = " << num_classical_crossings << endl;
		
			if (code_data.immersion == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && shortcut_crossing.size())
			{
				pure_knotoid_code_data = true;
				for (unsigned int i=0; i< shortcut_crossing.size(); i++)
				{
					if (shortcut_crossing[i] != 0)
						num_classical_crossings--;
				}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_orientation_data::constructor(generic_code_data): knotoid: num_classical_crossings = " << num_classical_crossings << endl;
			}
		}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_orientation_data::constructor(generic_code_data): pure_knotoid_code_data = " << pure_knotoid_code_data << endl;

		/* create the orientation matrix and write the number of crossings */
		num_components = code_data.num_components;
		num_terms_in_component = vector<int>(num_components);
		start_of_component = vector<int>(num_components);
		num_terms = 2*num_classical_crossings;
		immersion_crossing = vector<int>(num_classical_crossings);
		orientation_matrix = matrix<int>(2,num_terms);
		classical_gauss_data = vector<int>(num_terms);
		classical_crossing_sign = vector<int>(num_classical_crossings);
		
		if (track_zig_zag_counts)
			zig_zag_count = matrix<int>(2,num_classical_crossings);
		
	
		/* trace around the diagram writing the Gauss orientation data into the orientation matrix */
		int index = 0;  // used to write into the orientation matrix	
		int num_classical_crossings_visited = 0;
		
		/* classical_crossing will be used to renumber the immersion crossings,
		   so that if classical_crossing[i] = j, then crossing j of the immersion
		   is renumbered as crossing i in the Gauss code.
		   
		   crossing_visited will be a flag to indicate whether the ith immersion crossing
		   has been recorded in classical_crossing or not.
		*/
		vector<int> classical_crossing(num_classical_crossings);
		vector<int> crossing_visited(code_data.num_crossings);
	
		int start=0;
		int edge=0;
		int component=0;
		bool complete = false;
		
		do 
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_orientation_data::constructor(generic_code_data): start_of_component = " << start << endl;
	
			/* set the offset for the start of this component */
			if (component !=0)
			{
				start_of_component[component] = start_of_component[component-1]+num_terms_in_component[component-1];
			}
			else
			{
				start_of_component[component] = 0;
			}
			
			/*	trace this component */
			int component_terms = 0;
			do
			{	
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_orientation_data::constructor(generic_code_data): edge = " << edge;
				edge_flag[edge] = 1;
				int next_crossing = term_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", next_crossing = " << next_crossing;

				int next_crossing_orientation_polarity;
				int next_crossing_gauss_polarity;
				
				if (!extended_version && code_table[LABEL][next_crossing] == generic_code_data::VIRTUAL)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", is virtual" << endl;
				}
				else if (!extended_version && pure_knotoid_code_data && shortcut_crossing[next_crossing])
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", is a shortcut crossing" << endl;
				}
				else 
				{
					component_terms++;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", is a Gauss crossing" << endl;
	
				    if ((edge%2 == 0 && code_table[TYPE][next_crossing] == generic_code_data::TYPE1) ||
				        (edge%2 != 0 && code_table[TYPE][next_crossing] == generic_code_data::TYPE2))						
						next_crossing_orientation_polarity = gauss_orientation_data::LEFT;
					else
						next_crossing_orientation_polarity = gauss_orientation_data::RIGHT;

					if ((edge%2 == 1 && code_table[LABEL][next_crossing] == generic_code_data::POSITIVE) ||
						(edge%2 == 0 && code_table[LABEL][next_crossing] == generic_code_data::NEGATIVE))						
						next_crossing_gauss_polarity = generic_code_data::UNDER;
					else
						next_crossing_gauss_polarity = generic_code_data::OVER;	
				}
				
				if (!extended_version && (code_table[LABEL][next_crossing] == generic_code_data::VIRTUAL || (pure_knotoid_code_data && shortcut_crossing[next_crossing])))
				{
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_orientation_data::constructor(generic_code_data):   doing nothing" << endl;
					
					/* just move on around the component */				
					if (edge%2)
						edge = code_table[EVEN_ORIGINATING][next_crossing];
					else
						edge = code_table[ODD_ORIGINATING][next_crossing];
	
				}				
				else
				{
					if(crossing_visited[next_crossing])
					{
						for (int i=0; i< num_classical_crossings_visited; i++)
						{
							if (classical_crossing[i] == next_crossing)
							{
								orientation_matrix[0][index] = next_crossing_orientation_polarity;
								orientation_matrix[1][index] = i+1;
								
								if (code_table[LABEL][next_crossing] == generic_code_data::VIRTUAL || (pure_knotoid_code_data && shortcut_crossing[next_crossing]))
									orientation_matrix[1][index] *= -1;

								classical_gauss_data[index] = (i+1)*next_crossing_gauss_polarity;	


								if (track_zig_zag_counts)
								{
									zig_zag_count[0][i] = abs(peer_code_zig_zag_count[(edge%2?1:0)][next_crossing]);
								}
								
								index++;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_orientation_data::constructor(generic_code_data):   second visit to Gauss crossing " << i << endl;
								break;
							}
						}
					}
					else
					{
	
						classical_crossing[num_classical_crossings_visited] = next_crossing;
						crossing_visited[next_crossing] = 1;
						
						orientation_matrix[0][index] = next_crossing_orientation_polarity;
						orientation_matrix[1][index] = num_classical_crossings_visited+1;
	
						if (code_table[LABEL][next_crossing] == generic_code_data::VIRTUAL || (pure_knotoid_code_data && shortcut_crossing[next_crossing]))
							orientation_matrix[1][index] *= -1;

						classical_gauss_data[index] = (num_classical_crossings_visited+1)*next_crossing_gauss_polarity;	
	
						if (track_zig_zag_counts)
						{
							zig_zag_count[1][num_classical_crossings_visited] = abs(peer_code_zig_zag_count[(edge%2?1:0)][next_crossing]);
						}
						
						index++;
	
						immersion_crossing[num_classical_crossings_visited] = next_crossing;
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_orientation_data::constructor(generic_code_data):   first visit, becomes Gauss crossing " << num_classical_crossings_visited << endl;

						num_classical_crossings_visited++;
					}
	
					if (edge%2)
						edge = code_table[EVEN_ORIGINATING][next_crossing];
					else
						edge = code_table[ODD_ORIGINATING][next_crossing];
				}
			} while (edge != start);
	
			num_terms_in_component[component] = component_terms;
			component++;
			
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
		} while (!complete);

		for (int i=0; i< num_classical_crossings; i++)
		{
			int crossing = classical_crossing[i];
			
			if (code_table[LABEL][crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[crossing])) // it should never be!
			{
				if ((code_table[TYPE][crossing] == generic_code_data::TYPE1 && code_table[LABEL][crossing] == generic_code_data::NEGATIVE) ||
				    (code_table[TYPE][crossing] == generic_code_data::TYPE2 && code_table[LABEL][crossing] == generic_code_data::POSITIVE))
					classical_crossing_sign[i] = 1;
				else
					classical_crossing_sign[i] = -1;
			}
		}

	}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_orientation_data::constructor(generic_code_data): generic_orientation_data" << endl;
	print_gauss_data(*this,debug,"gauss_orientation_data::constructor(generic_code_data):   ");
}
}

void gauss_orientation_data::reverse_all()
{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	print_gauss_data(*this, debug, "gauss_orientation_data::reverse_all: forward ");	

		
	matrix<int> rev_m(2,num_terms);
	
	for (int i=0; i< num_terms; i++)
	{
		rev_m[0][num_terms-1-i] = orientation_matrix[0][i];
		rev_m[1][num_terms-1-i] = orientation_matrix[1][i];
	}
	
	orientation_matrix = rev_m;
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	print_gauss_data(*this, debug, "gauss_orientation_data::reverse_all: reversed ");	
}

/* reverse the orientation of those components whose reverse_flag is set */
void gauss_orientation_data::reverse_components( vector<bool>& reverse_flag)
{
	
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "gauss_orientation_data::reverse_components: reverse_flag ";	
	for (int i=0; i< num_components; i++)
		debug << reverse_flag[i] << " ";
	debug << endl;
}

	matrix<int> new_orientation_matrix = orientation_matrix;
	
	for (int i=0; i< num_components; i++)
	{
		if (reverse_flag[i] == false)
			continue;
			
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "gauss_orientation_data::reverse_components: reversing component" << i << endl;	
			
		int first = start_of_component[i];
		int last = start_of_component[i]+num_terms_in_component[i]-1;

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "gauss_orientation_data::reverse_components:   first component offset = " << first << " last component offset = " << last << endl;	
		
		for (int j=0; j< num_terms_in_component[i]; j++)
		{
			new_orientation_matrix[1][first+j] = orientation_matrix[1][last-j];
			
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "gauss_orientation_data::reverse_components:   term " << j << " new crossing " << new_orientation_matrix[1][first+j] << endl;	
			
			/* look in orientation matrix for the other place of this crossing, since we are still in the middle of re-writing 
			   new_orientation_matrix, so we may have three instances of the crossing in new_orientation_matrix at this point */
			int component = 0;
			int place=0;
			for (int k=0; k< num_terms; k++)
			{
				if (k==start_of_component[component]+num_terms_in_component[component])
					component++;
					
				if (k != last-j && orientation_matrix[1][k] == orientation_matrix[1][last-j])
				{
					place = k;
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "gauss_orientation_data::reverse_components:   found crossing at place " << place << " in component " << component << endl;	
					break;
				}
			}
			
			if (reverse_flag[component] == false)
			{
				new_orientation_matrix[0][first+j] = orientation_matrix[0][place];
				new_orientation_matrix[0][place] = orientation_matrix[0][last-j];
			}
			else
				new_orientation_matrix[0][first+j] = orientation_matrix[0][last-j]; 
				
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "gauss_orientation_data::reverse_components:   polarity of crossing at offset " << place
	      << (new_orientation_matrix[0][place] == gauss_orientation_data::LEFT? " LEFT": " RIGHT") << endl;	
	debug << "gauss_orientation_data::reverse_components:   polarity of crossing at offset " << first+j
	      << (new_orientation_matrix[0][first+j] == gauss_orientation_data::LEFT? " LEFT": " RIGHT") << endl;	
}
			
		}
	}
	orientation_matrix = new_orientation_matrix;
}

void print_gauss_data(gauss_orientation_data g, ostream& os, string prefix)
{
	os << prefix << "number of terms = " << g.num_terms << endl;
	os << prefix << "number of components = " << g.num_components << endl;
	os << prefix << "extended_version = " << (g.extended_version? "true":"false") << endl;
	os << prefix << "number of terms in component: ";
	for (int i=0; i< g.num_components; i++)
		os << g.num_terms_in_component[i] << ' ';
	os << endl;
	os << prefix << "start of component: ";
	for (int i=0; i< g.num_components; i++)
		os << g.start_of_component[i] << ' ';
	os << endl;
	os << prefix << "immersion_crossing: ";
	for (int i=0; i< g.num_terms/2; i++)
		os << g.immersion_crossing[i] << ' ';
	os << endl;
	os << prefix << "orientation matrix: " << endl;
	os << prefix;
	for (int i=0; i< g.num_terms; i++)
		os << (g.orientation_matrix[0][i] == gauss_orientation_data::LEFT? "  L" : "  R");
	os << endl << prefix;
	for (int i=0; i< g.num_terms; i++)
		os << setw(3) << g.orientation_matrix[1][i];
	os << endl;
	os << prefix << "classical_gauss_data: ";	
	for (unsigned int i=0; i< g.classical_gauss_data.size(); i++)
		os << g.classical_gauss_data[i] << ' ';
	os << endl;
	os << prefix << "classical_crossing_sign: ";
	for (unsigned int i=0; i< g.classical_crossing_sign.size(); i++)
		os << g.classical_crossing_sign[i] << ' ';
	os << endl;
	
	
	if (g.zig_zag_count.numrows() != 0)
	{
		os << prefix << "zig_zag_count: " << endl;
		print (g.zig_zag_count, os, 3, prefix);
	}	
}

/* Note that when we write a gauss_orientation_data object containing zig_zag_count data to an ostream, we write the
   zig-zag counts in the order corresponding to the terms of the orientation matrix and not indexed by the crossing
   numbers.
*/
void write_gauss_data(gauss_orientation_data g, ostream& os, bool zig_zags, int immersion_character, int head_zig_zag_count)
{
	int index = 0;
	for (int k=0; k < g.num_components; k++)
	{
		for (int l=0; l < g.num_terms_in_component[k]; l++)
		{
			os << (g.orientation_matrix[0][index] == gauss_orientation_data::LEFT? "L" : "R");
			os << g.orientation_matrix[1][index];
			
			if (l < g.num_terms_in_component[k]-1)
				os << ' ';
			index++;
		}
		if (k < g.num_components-1)
			os << ", ";
	}
		
	if (zig_zags && g.zig_zag_count.numrows() != 0)
	{		
		index=0;
		os << " (";
		vector<bool> first_occurrence_written(g.num_terms/2);
		for (int k=0; k < g.num_components; k++)
		{
			for (int l=0; l < g.num_terms_in_component[k]; l++)
			{			
				int crossing = g.orientation_matrix[1][index]-1;
				if(!first_occurrence_written[crossing])
				{
					os << g.zig_zag_count[1][crossing];
					first_occurrence_written[crossing] = true;
				}
				else
				{
					os << g.zig_zag_count[0][crossing];
				}
				
				if (l == g.num_terms_in_component[k]-1 && k ==0 && immersion_character != generic_code_data::character::CLOSED)
					os << ' ' << abs(head_zig_zag_count);
					
				if (l == g.num_terms_in_component[k]-1 && k < g.num_components-1)
					os << ',';

				if (index < g.num_terms-1)
					os << ' ';
				
				index++;
			}
		}
		
		os << ")";
	}
}

/* The left preferred (Gauss) representation of an immersion diagram is a Gauss code for the diagram having L1,...,Ln as a 
   subsequence that is minimal, lexicographically (where L < R), amongst all such Gauss codes.

   To evaluate the left preferred representation we consider every possible re-numbering of the crossings and (optionally) all 
   possible orientations of the diagram's unicursal components.   If the diagram is that of a knotoid or a long knot, the function
   assumes that the first component of the parameter gauss_data is the segment component, always considers that component first, 
   and always considers that component by starting at the crossing relating to the first term of gauss_data.
   
   If the diagram is not that of a knotoid or a long knot, the function considers all permutations of the order of the components 
   of the diagram and, for each component order, cycles each component so we start considering that component by starting at each 
   crossing.  Each such component order and component cycle determines a crossing permutation that results in the 
   'L' instances appearing in the order L1,...Ln and we can apply that permutation to the 'R' instances to obtain the corresponding 
   overall sequence.  We use the lexicographic vector comparison operator < to determine the component order and component cycle 
   combination that yields the minimal sequence.  
   
   In the case where the the diagram is not that of a knotioid or a long knot, this approach will always produce a left_preferred
   Gauss code that starts with the term L1.  If the diagram is a knotoid or a long knot, this may not be the case.
   
   Note: we cannot simply cycle the components starting with 'L' instances, since virtual diagrams may not have any 'L' instances 
   in some components.  For example: [9 -7 1, -3 5]/# * * # * has Gauss data R1 R2, L2 L1.
   
*/
gauss_orientation_data left_preferred(const gauss_orientation_data& gauss_data, bool unoriented, int immersion_character)
{

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "left_preferred: initial Gauss data ";	
	write_gauss_data(gauss_data,debug);
	debug << endl;
	print_gauss_data(gauss_data, debug, "left_preferred: ");
}
	
	vector<pair<int,int> > min_perm_weight;  // initializes to zero size
	vector<int> min_perm;
	vector<int> min_component_perm;
	vector<int> min_component_cycle;
	vector<int> min_component_start;
	vector<bool> min_r_flag;
	matrix<int> min_zig_zag_count;

	const matrix<int>& zig_zag_count = gauss_data.zig_zag_count;
	bool track_zig_zag_counts;
	
	if (zig_zag_count.numrows() !=0)
		track_zig_zag_counts = true;
	else
		track_zig_zag_counts = false;

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred: track_zig_zag_counts = " << track_zig_zag_counts << endl;			
	
	int num_components = gauss_data.num_components;
	int num_terms=gauss_data.num_terms;
	int min_component = (immersion_character == generic_code_data::character::CLOSED? 0: 1);

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred: min_component = " << min_component << endl;			
	
	vector<bool> r_flag(num_components); // initializes to false
	bool finished = false;
	
	gauss_orientation_data g = gauss_data;  // will eventually be the return value
	
	do
	{
		/* consider the orientation determined by the r_flags */
		
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred: reverse flags: ";
	for (int i=0; i< num_components; i++)
		debug << r_flag[i] << " ";
	debug << endl;
}
	
		matrix<int>& m = g.orientation_matrix;
	
		
		/* consider the components in every posible order */
		vector<int> component_perm(g.num_components);
		for (int i=0; i< g.num_components; i++)
			component_perm[i]=i;
		vector<int>::iterator component_perm_start = component_perm.begin();
		
		if (min_component !=0)
			component_perm_start++;
	
		do
		{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:   component_perm = ";	
	for (int i=0; i< g.num_components; i++)
		debug << component_perm[i] << ' ';
	debug << endl;
}	
			/* component_cycle will record the cyclic rotation of each component, enumarated
			   in the order given by g, by recording which term of the component we are going 
			   to start at.  We will then need to find this instance and evaluate the starting 
			   offset in g.orientation_matrix
			*/
			vector<int> component_cycle(g.num_components); // initializes to zero
	
			bool complete = false;		
			do
			{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:     component_cycle = ";	
	for (int i=0; i< g.num_components; i++)
		debug << component_cycle[i] << ' ';
	debug << endl;
}	
				/* begin by identifying the cyclic starting location for ordering each component 
				   (enumerated in the order given by g) according to component_cycle.
				*/
				vector<int> component_start(g.num_components);
				for (int i=0; i< g.num_components; i++)
					component_start[i]=g.start_of_component[i]+component_cycle[i];
				
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:       component_start: ";
	for (int k=0; k< g.num_components; k++)
		debug << component_start[k] << " ";
	debug << endl;
}
	
				/* evaluate the crossing permutation: perm[i] indicates the new crossing number of the old ith crossing */
				vector<int> perm(num_terms/2);
				int new_crossing = 0;
				
				for (int k=0; k< g.num_components; k++)
				{
					int component = component_perm[k];
					int num_component_terms = g.num_terms_in_component[component];
					int first_term = g.start_of_component[component];
					int offset = component_start[component]-first_term;

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:       component = " << component << ", offset " << setw(2) << offset << ", num_component_terms = " 
		  << num_component_terms << ", first_term = "  << first_term << endl;
	debug << "left_preferred:       component indices:";				
}	
					for (int i=0; i < num_component_terms; i++)
					{
						int index = first_term+(offset+i)%num_component_terms;
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "  " << index;
						if (m[0][index] == gauss_orientation_data::LEFT)
						{
							perm[abs(m[1][index])-1] = new_crossing++;
						}
					}
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << endl;			
				}

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:       perm = ";	
	for (int i=0; i< num_terms/2; i++)
		debug << perm[i] << ' ';
	debug << endl;
}

				/* Evaluate the lexicographic weight of the permutation. */
				vector<pair<int,int> > perm_weight(num_terms);
				int place = 0;
				for (int k=0; k< g.num_components; k++)
				{
					int component = component_perm[k];
					int num_component_terms = g.num_terms_in_component[component];
					int first_term = g.start_of_component[component];
					int offset = component_start[component]-first_term;
	
					for (int i=0; i < num_component_terms; i++)
					{
						int index = first_term+(offset+i)%num_component_terms;
	
						perm_weight[place] = {perm[abs(m[1][index])-1],m[0][index]};
						place++;
					}
				}
	
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:       perm_weight = ";
	for (int i=0; i< num_terms; i++)
		debug << '(' << perm_weight[i].first << "," << perm_weight[i].second << ")" << ' ';
	debug << endl;
}
				matrix<int> new_zig_zag_count(2,num_terms/2);
				if (track_zig_zag_counts)
				{
					/* To adjust the zig-zag counts to reflect the left preferred Gauss code, we need to consider the effect of perm, 
					   component_perm, r_flag, and component_start on the order in which we encounter the crossings.  The parity 
					   of the crossing occurrences is not relevant.
					   
					   perm simply permutes the crossing numbers, we evaluate the effect of this after the other factors have been taken 
					   into account.
					   
					   If the first and second occurences of a crossing occur in different components then the relative position of the first 
					   and second occurrences in the modified code is dependent only on the component_perm and not whether any particular
					   component has been cycled or reversed.  Similarly, if the first and second occurrences are in the same component, then the
					   relative position in the modified code is dependent only on whether that component has been cycled or reversed and not on
					   component_perm.
					   
					   Moreover, if we reverse the orientation of a component, then for each crossing on that component, the irreducible 
					   cusps encountered before the crossing when going in the reverse direction is the number of cusps encountered after the
					   crossing when going in the original direction.
					*/
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred: adjusting gauss_data.zig_zag_count" << endl;
	
					matrix<int> initial_new_zig_zag_count(2,num_terms/2);
					for (int i = 0; i< num_terms/2; i++)
					{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred: crossing " << i+1 << endl;

						/* look for the component and places of the two occurrences of crossing i+1 in the original gauss_data */
						int first_component = -1;
						int second_component = -1;
						int first_place;
						int second_place;
						
						int component = 0;
						for (int j=0; j< num_terms; j++)
						{			
							if (j >= gauss_data.start_of_component[component]+gauss_data.num_terms_in_component[component])
								component++;
								
							if (gauss_data.orientation_matrix[1][j] == i+1)
							{
								if (first_component == -1)
								{
									first_component = component;
									first_place = j;
								}
								else
								{
									second_component = component;
									second_place = j;
									break;
								}
							}
						}

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:   first_component = " << first_component << " second_component = " << second_component << endl;
	debug << "left_preferred:   first_place = " << first_place << " second_place = " << second_place << endl;
}
			
		
						/* if we have reversed the direction of the component then we have to adjust for the fact that we're going in 
						   the opposite direction now.  We have to look at the term in the original Gauss code following the first(second) 
						   place and use the counts from that term.  Note that we need to identify whether the following term is the 
						   first or second occurrence of that crossing in the original Gauss code.					   
						*/
						int first_subsequent_place = (first_place+1-gauss_data.start_of_component[first_component])%gauss_data.num_terms_in_component[first_component] + gauss_data.start_of_component[first_component];
						int second_subsequent_place = (second_place+1-gauss_data.start_of_component[second_component])%gauss_data.num_terms_in_component[second_component] + gauss_data.start_of_component[second_component];
						int first_subsequent_term = gauss_data.orientation_matrix[1][first_subsequent_place];
						int second_subsequent_term = gauss_data.orientation_matrix[1][second_subsequent_place];
						
						bool first_subsequent_occurrence;
						for (int j=0; j< num_terms; j++)
						{
							if (gauss_data.orientation_matrix[1][j] == first_subsequent_term)
							{
								if (j == first_subsequent_place)
									first_subsequent_occurrence = true;
								else
									first_subsequent_occurrence = false;
								
								break;
							}
						}
						
						bool second_subsequent_occurrence;
						for (int j=0; j< num_terms; j++)
						{
							if (gauss_data.orientation_matrix[1][j] == second_subsequent_term)
							{
								if (j == second_subsequent_place)
									second_subsequent_occurrence = true;
								else
									second_subsequent_occurrence = false;
								
								break;
							}
						}

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:   first_subsequent_place = " << first_subsequent_place << " second_subsequent_place = " << second_subsequent_place << endl;
	debug << "left_preferred:   first_subsequent_term = " << first_subsequent_term << " second_subsequent_term = " << second_subsequent_term << endl;
	debug << "left_preferred:   first_subsequent_occurrence = " << first_subsequent_occurrence << " second_subsequent_occurrence = " << second_subsequent_occurrence << endl;
}
			
						if (first_component == second_component)
						{
							/* we need to consider the relative position of the first and second occurrences (1&2) relative to the offset (O)
							   of the component_start within the component (i.e the component_cycle).  
							   
							   The occurrences are reversed (Y/N) as follows:
							   
							   component reversed:       O..1..2 Y   1..2..O Y   1..O..2 N    O=1..2 N    1..2=O Y
							   component not reversed:   O..1..2 N   1..2..O N   1..O..2 Y    O=1..2 N    1..2=O Y	

							   Note, however, that if the component has been reversed, we have to adjust the calculation of the offset to 
							   accomodate the fact that component_start relates to g (which has reversed components), not gauss_data;
							   
							*/
							int num_component_terms = gauss_data.num_terms_in_component[first_component];
//							int first_term = gauss_data.start_of_component[first_component];
//							int offset = component_start[first_component]-first_term;
							int offset;
							if (r_flag[first_component]) 
								offset = gauss_data.start_of_component[first_component]+num_component_terms-1-component_cycle[first_component];
//								offset = num_component_terms-component_cycle[first_component];
							else
								offset = component_start[first_component];
//								offset = component_cycle[first_component];
								
							bool reversed_occurrences = false;
							
							if ( offset == second_place ||
								(r_flag[first_component] && (first_place > offset || second_place < offset)) ||
								(!r_flag[first_component] && (first_place < offset && second_place > offset))
								)
							{
								 reversed_occurrences = true;
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred:   same component, offset = " << offset << " reversed occurrences" << endl;
							}
							else
							{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred:   same component, offset = " << offset << " not reversed" << endl;
							}
							
							if (reversed_occurrences)
							{
								
								if (r_flag[first_component])
								{
									/* Since we are reversing the first and second occurrences of crosing i, we base the count in row zero (the
									   second occurrence counter) on the subsequent term to the first occurrence in the original code.  
									*/
									initial_new_zig_zag_count[0][i] = zig_zag_count[(first_subsequent_occurrence?1:0)][first_subsequent_term-1];
									initial_new_zig_zag_count[1][i] = zig_zag_count[(second_subsequent_occurrence?1:0)][second_subsequent_term-1];

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:     component reversed, new first count for crossing " << i+1 << " is old " 
	      << (second_subsequent_occurrence? "first": "second") << " count for crossing " << second_subsequent_term << endl;
	debug << "left_preferred:                         new second count for crossing " << i+1 << " is old " 
	      << (first_subsequent_occurrence? "first": "second") << " count for crossing " << first_subsequent_term << endl;
}
									
								}
								else
								{
									initial_new_zig_zag_count[0][i] = zig_zag_count[1][i];
									initial_new_zig_zag_count[1][i] = zig_zag_count[0][i];

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred:     component not reversed, new first(second) count for crossing " << i+1 << " is old second(first) count for crossing " << i+1 << endl; 

								}
							}
							else
							{
								if (r_flag[first_component])
								{
									/* we still have to adjust for the fact that we're going in the opposite direction now but since
									   we are not reversing the first and second occurrences of crosing i, we base the count in row one
									   (the first occurrence counter) on the subsequent term to the first occurrence in the original code.  
									*/
									initial_new_zig_zag_count[1][i] = zig_zag_count[(first_subsequent_occurrence?1:0)][first_subsequent_term-1];
									initial_new_zig_zag_count[0][i] = zig_zag_count[(second_subsequent_occurrence?1:0)][second_subsequent_term-1];
									
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:     component reversed, new first count for crossing " << i+1 << " is old " 
	      << (first_subsequent_occurrence? "first": "second") << " count for crossing " << first_subsequent_term << endl;
	debug << "left_preferred:                         new second count for crossing " << i+1 << " is old " 
	      << (second_subsequent_occurrence? "first": "second") << " count for crossing " << second_subsequent_term << endl;
}
								}
								else
								{
									initial_new_zig_zag_count[0][i] = zig_zag_count[0][i];
									initial_new_zig_zag_count[1][i] = zig_zag_count[1][i];

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred:     component not reversed, new first(second) count for crossing " << i+1 << " is old first(second) count for crossing " << i+1 << endl; 
								}				
							}				
						}
						else
						{
							/* the first and second occurrences are reversed if second_component appears before first_component */
							bool reversed_occurrences;
							for (int j=0; j < num_components; j++)
							{
								if (component_perm[j] == first_component)
								{
									reversed_occurrences = false;
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred:   different components, not reversed" << endl;
									break;
								}
								else if (component_perm[j] == second_component)
								{
									reversed_occurrences = true;
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred:   different components, reversed occurrences" << endl;
									break;
								}
							}
							
							
							if (reversed_occurrences)
							{
								if (r_flag[second_component])
								{
									initial_new_zig_zag_count[1][i] = zig_zag_count[(second_subsequent_occurrence?1:0)][second_subsequent_term-1];
									
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:     second component reversed, new first count for crossing " << i+1 << " is old " 
	      << (second_subsequent_occurrence? "first": "second") << " count for crossing " << second_subsequent_term << endl;
}
								}
								else
								{
									initial_new_zig_zag_count[1][i] = zig_zag_count[0][i];
									
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred:     second component not reversed, new first count for crossing " << i+1 << " is old second count for crossing " << i+1 << endl; 
									
								}

								if (r_flag[first_component])
								{
									initial_new_zig_zag_count[0][i] = zig_zag_count[(first_subsequent_occurrence?1:0)][first_subsequent_term-1];
									
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:     first component reversed, new second count for crossing " << i+1 << " is old " 
	      << (first_subsequent_occurrence? "first": "second") << " count for crossing " << first_subsequent_term << endl;
}
								}
								else
								{
									initial_new_zig_zag_count[0][i] = zig_zag_count[1][i];
									
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred:     first component not reversed, new second count for crossing " << i+1 << " is old first count for crossing " << i+1 << endl; 
									
								}

//								initial_new_zig_zag_count[0][i] = zig_zag_count[1][i];
//								initial_new_zig_zag_count[1][i] = zig_zag_count[0][i];
							}
							else
							{
								if (r_flag[first_component])
								{
									initial_new_zig_zag_count[1][i] = zig_zag_count[(first_subsequent_occurrence?1:0)][first_subsequent_term-1];
									
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:     first component reversed, new first count for crossing " << i+1 << " is old " 
	      << (first_subsequent_occurrence? "first": "second") << " count for crossing " << first_subsequent_term << endl;
}
								}
								else
								{
									initial_new_zig_zag_count[1][i] = zig_zag_count[1][i];

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred:     first component not reversed, new first count for crossing " << i+1 << " is old first count for crossing " << i+1 << endl; 
									
								}

								if (r_flag[second_component])
								{
									initial_new_zig_zag_count[0][i] = zig_zag_count[(second_subsequent_occurrence?1:0)][second_subsequent_term-1];
									
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred:     second component reversed, new second count for crossing " << i+1 << " is old " 
	      << (second_subsequent_occurrence? "first": "second") << " count for crossing " << second_subsequent_term << endl;
}
								}
								else
								{
									initial_new_zig_zag_count[0][i] = zig_zag_count[0][i];
									
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "left_preferred:     second component not reversed, new second count for crossing " << i+1 << " is old second count for crossing " << i+1 << endl; 
									
								}
							}
						}
					}
		
					/* now adjust to reflect the crossing perm and clear any minus signs */
					for (int i = 0; i< num_terms/2; i++)
					{
						new_zig_zag_count[0][perm[i]] = initial_new_zig_zag_count[0][i];
						new_zig_zag_count[1][perm[i]] = initial_new_zig_zag_count[1][i];
					}
		
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred: adjusted initial_new_zig_zag_count by perm:" << endl;
	print(new_zig_zag_count, debug,3,"left_preferred: ");
}
				}

				
				if (min_perm_weight.size() == 0)
				{
					min_perm_weight = perm_weight;
					min_perm = perm;
					min_component_perm = component_perm;
					min_component_cycle = component_cycle;
					min_component_start = component_start;
					min_r_flag = r_flag;
					
					if (track_zig_zag_counts)
						min_zig_zag_count = new_zig_zag_count;
				}
				else
				{
					if (perm_weight < min_perm_weight || (perm_weight == min_perm_weight && new_zig_zag_count < min_zig_zag_count))
					{
						min_perm_weight = perm_weight;
						min_perm = perm;
						min_component_perm = component_perm;
						min_component_cycle = component_cycle;
						min_component_start = component_start;
						min_r_flag = r_flag;
						
						if (track_zig_zag_counts)
							min_zig_zag_count = new_zig_zag_count;
					}
				}
			
				if (g.num_components == 1 && min_component == 1)
				{
					break;
				}
				else
				{
					/* increment component_cycle lexicographically up to the number of terms in each component.	*/
					for (int k= g.num_components-1; k >=min_component; k--)
					{
						component_cycle[k]++;
						
						if (component_cycle[k] >= g.num_terms_in_component[k])
						{
							component_cycle[k] = 0;
							if (k==min_component)
							{
								complete = true;
								break;
							}
						}
						else
							break;
					}
				}
			} while (!complete);
		} while(next_permutation(component_perm_start,component_perm.end()));

		/* if we're writing an unoriented Gauss code, consider other orientations of the components
		   if the immersion character is generic_code_data::character::CLOSED, we can reverse the 
		   orientation of all components, if it is not then we cannot reverse the orientation of 
		   component zero
		*/
		if (unoriented)
		{			
			finished = true;
			for (int i=num_components-1; i>=min_component; i--)
			{
				if (r_flag[i])
					r_flag[i] = false;
				else
				{
					r_flag[i] = true;
					finished = false;
					break;
				}
			}
			
			if (!finished)
			{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred: advanced r_flag to: ";
	for (int i=0; i< num_components; i++)
		debug << r_flag[i] << " ";
	debug << endl;
}
				g = gauss_data;
				g.reverse_components(r_flag);

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred: updated Gauss data ";	
	write_gauss_data(g,debug);
	debug << endl;
}
			}
		}
		else
		{
			break;
		}
		
	} while (!finished);
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "left_preferred: min_r_flag = ";
	for (int i=0; i< num_components; i++)
		debug << min_r_flag[i] << " ";
	debug << endl;
	debug << "left_preferred: min_component_perm = ";	
	for (int i=0; i< num_components; i++)
			debug << min_component_perm[i] << ' ';
	debug << endl;
	debug << "left_preferred: min_component_cycle = ";	
	for (int i=0; i< num_components; i++)
			debug << min_component_cycle[i] << ' ';
	debug << endl;
	debug << "left_preferred: min_component_start = ";	
	for (int i=0; i< num_components; i++)
			debug << min_component_start[i] << ' ';
	debug << endl;
	debug << "left_preferred: min_perm = ";	
	for (int i=0; i< num_terms/2; i++)
		debug << min_perm[i] << ' ';
	debug << endl;
	debug << "left_preferred: min_perm_weight = ";
	for (int i=0; i< num_terms; i++)
		debug << '(' << min_perm_weight[i].first << "," << min_perm_weight[i].second << ")" << ' ';
	debug << endl;
}

	/* write the new orientation matrix from the crossing permutation of minimum weight */
	g = gauss_data;
	g.reverse_components(min_r_flag);
	matrix<int>& m = g.orientation_matrix;
	
	vector<int>	new_num_terms_in_component(num_components); 
	vector<int>	new_start_of_component(num_components); 
	vector<int>	new_immersion_crossing(num_terms/2); 
	
	matrix<int> new_m(2,num_terms);
/*	for (int i=0; i< num_terms; i++)
	{
		new_m[0][i] = m[0][i];
		new_m[1][i] = min_perm[abs(m[1][i])-1]+1;
		
		if (m[1][i] < 0)
			new_m[1][i] *= -1;
	}
*/
	int new_index = 0;
	for (int k=0; k< g.num_components; k++)
	{
		int component = min_component_perm[k];
		int num_component_terms = g.num_terms_in_component[component];
		int first_term = g.start_of_component[component];
		int offset = min_component_start[component]-first_term;

		new_num_terms_in_component[k] = num_component_terms;
		
		for (int i=0; i < num_component_terms; i++)
		{
			int old_index = first_term+(offset+i)%num_component_terms;
			{
				new_m[0][new_index] = m[0][old_index];
				new_m[1][new_index] = min_perm[abs(m[1][old_index])-1]+1;

				if (m[1][old_index] < 0)
					new_m[1][new_index] *= -1;		
			}
			new_index++;
		}			
	}

	for (int i=1; i< num_components; i++)
		new_start_of_component[i] = new_start_of_component[i-1]+new_num_terms_in_component[i-1];
		
	for (int i=0; i< num_terms/2; i++)
		new_immersion_crossing[min_perm[i]] = gauss_data.immersion_crossing[i];
	
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "left_preferred: new_m" << endl;	
	debug << "left_preferred:   ";	
	for (int i=0; i< num_terms; i++)
		debug << (new_m[0][i] == gauss_orientation_data::LEFT? " L" : " R");
	debug << "\nleft_preferred:   ";	
	for (int i=0; i< num_terms; i++)
		debug << setw(2) << new_m[1][i];
	debug << endl;

}	

	g.num_terms_in_component = new_num_terms_in_component;
	g.start_of_component = new_start_of_component;
	g.immersion_crossing = new_immersion_crossing;
	g.orientation_matrix = new_m;
	
	if (track_zig_zag_counts)
		g.zig_zag_count = min_zig_zag_count;


if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "left_preferred: returning gauss data: " << endl;	
	print_gauss_data(g, debug, "left_preferred: ");
}	
	
	return g;
}


/* if two sets of gauss_orientation_data are equal they represent the same diagram, in which case they 
   will yield the same left-preferred representation of the corresponding Gauss code.
   
   Note: this operator ignores the zig-zag counts within the gauss_orientation_data
*/
//bool operator == const (gauss_orientation_data& a, gauss_orientation_data& b)
bool gauss_orientation_data::operator == (gauss_orientation_data& b) const
{

//if (debug_control::DEBUG >= debug_control::DETAIL)
if (debug_control::DEBUG >= debug_control::BASIC)
{
	print_gauss_data(*this, debug, "gauss_orientation_data::operator == : a: ");	
	print_gauss_data(b, debug, "gauss_orientation_data::operator == : b: ");	
}
	gauss_orientation_data a_left_preferred = left_preferred(*this);
	gauss_orientation_data b_left_preferred = left_preferred(b);
	ostringstream asstr;
	write_gauss_data(a_left_preferred,asstr);
	ostringstream bsstr;
	write_gauss_data(b_left_preferred,bsstr);

//if (debug_control::DEBUG >= debug_control::DETAIL)
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "gauss_orientation_data::operator == : left_preferred a: " << asstr.str() << endl;	
	debug << "gauss_orientation_data::operator == : left_preferred b: " << bsstr.str() << endl;	
}
	bool equal = (asstr.str() == bsstr.str());

//if (debug_control::DEBUG >= debug_control::DETAIL)
if (debug_control::DEBUG >= debug_control::BASIC)
		debug << "gauss_orientation_data::operator == : equal = " << equal << endl;	

	return equal;	
}

