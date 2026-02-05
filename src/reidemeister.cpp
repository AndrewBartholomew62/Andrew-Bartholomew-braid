/**********************************************************************************************************
	  This module contains support functions relating to Reidemeister moves, as used by the braid programme 
	  
bool virtual_Reidemeister_I_present(generic_code_data code_data)
int remove_virtual_Reidemeister_I(generic_code_data& code_data, vector<int>& component_flags)
int remove_edge_flags_from_peer_code(generic_code_data& code_data, vector<int>& edge_flag, vector<int>& component_flags);
bool Reidemeister_II_labels (const int label1, const int label2)
bool simple_Reidemeister_II_labels (const int label1, const int label2)
bool Gauss_Reidemeister_II_labels(int start_edge, int end_bigon_edge_a, int end_bigon_edge_b, generic_code_data& code_data)
simple_Reidemeister_II_return simple_Reidemeister_II_present(generic_code_data& code_data, bool create_flags)
Gauss_Reidemeister_II_return Gauss_Reidemeister_II_present(generic_code_data& code_data, gauss_orientation_data& gauss_data)
virtual_Reidemeister_II_plus_return virtual_Reidemeister_II_plus(generic_code_data& code_data, matrix<int>& cycle, int num_cycles, int num_left_cycles, bool mark_special=false)
int remove_Reidemeister_II(generic_code_data& code_data, vector<int>& component_flags)	
void align_label_numbers (matrix<int>& edge_labels, int num_components,vector<int> new_first_component_edge,vector<int> new_last_component_edge)
void clear_component_flag(int component, vector<int>& component_flags)
Reidemeister_III_return Reidemeister_III_present (generic_code_data code_data)
**************************************************************************/
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <valarray>
#include <list>
#include <iomanip>
#include <ctype.h>
#include <stdio.h>
#include <algorithm>
#include <tuple>

using namespace std;

extern ofstream     debug;
extern ofstream     output;
extern ifstream     input;

#include <util.h>
#include <matrix.h>
#include <generic-code.h>
#include <debug-control.h>
#include <gauss-orientation.h>
#include <reidemeister.h>

//bool gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data, bool optimal = true);
bool gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data, bool optimal=true, vector<int>* gauss_crossing_perm=0, bool evaluate_gauss_crossing_perm=false);

/* virtual_Reidemeister_I_present detects virtual Reidemeister I monogons and also virtual Reidemeister I 
   detours that, in the case of knotoids, may also encounter the shortcut.
    
   The function creates extended gauss_orientation_data from the code_data and then looks for a 1-reducible
   configuration where consecutive terms in the Gauss data lie in the same component and involve the same 
   crossing.  For the virtual case, this means looking for a sequence Li,...,Ri or Ri,...,Li in which 
   crossing i and any intervening crossings are virtual or shortcut crossings.

*/
bool virtual_Reidemeister_I_present(generic_code_data code_data)
{

	if (code_data.type == generic_code_data::gauss_code)
	{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "virtual_Reidemeister_I_present: presented with a Gauss code, returning false;" << endl;
		
		return false;
	}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "virtual_Reidemeister_I_present: presented with code data ";
	write_code_data(debug,code_data);	
	debug << endl;
}
	gauss_orientation_data gauss_data(code_data, true); // true gives us the extended version
	int num_terms = gauss_data.num_terms;
	matrix<int>& orientation_matrix = gauss_data.orientation_matrix;
	vector<int> start_of_component = gauss_data.start_of_component;
	vector<int> num_terms_in_component = gauss_data.num_terms_in_component;

	vector<int>& shortcut_crossing = code_data.shortcut_crossing;
	bool ignore_shortcut = false;
	if (code_data.head != -1 && shortcut_crossing.size())
		ignore_shortcut = true;
	
	for (int i=0; i< num_terms; i++)
	{
		int immersion_crossing = gauss_data.immersion_crossing[abs(orientation_matrix[1][i])-1];
		
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "virtual_Reidemeister_I_present: term " << i << " immersion_crossing = " << immersion_crossing << ", label = " 
	      << code_data.code_table[generic_code_data::table::LABEL][immersion_crossing] << endl;
}	
		/* we're only interested in terms related to virtual and shortcut crossings */
		if (code_data.code_table[generic_code_data::table::LABEL][immersion_crossing] != generic_code_data::VIRTUAL && !(ignore_shortcut && shortcut_crossing[immersion_crossing]))
		{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "virtual_Reidemeister_I_present: term " << i << " does not correspond to a virtual crossing" << endl;
	
			continue;
		}
			
		int ith_term_component=0;
		for (int j=1; j<gauss_data.num_components; j++)
		{
			if (i >= start_of_component[j])
				ith_term_component++;
			else
				break;
		}
		
		/* we require only virtual or shortcut crossings between this term and the term 
		   corresponding to the other visit to this crossing.  If a virtual Reidemeister I configuation
		   exists, it may lie between the first and second occurrence in the code, or between the second and 
		   first (wrapping around the code).
		*/
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "virtual_Reidemeister_I_present: checking term " << i << ", component " << ith_term_component << endl;

		for (int j=1; j< num_terms-1; j++)
		{
			if (orientation_matrix[1][(i+j)%num_terms] > 0)
			{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "virtual_Reidemeister_I_present:   term " << (i+j)%num_terms << " not virtual, no virtual Reidemeister I at term  " << i << endl;
				break;
			}
			else if (orientation_matrix[1][(i+j)%num_terms] == orientation_matrix[1][i])
			{
				int succeeding_term = (i+j)%num_terms;			
				int succeeding_term_component = 0;
		
				for (int k=1; k<gauss_data.num_components; k++)
				{
					if (succeeding_term >= start_of_component[k])
						succeeding_term_component++;
					else
						break;
				}
		
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "virtual_Reidemeister_I_present:   term " << (i+j)%num_terms << " is the second encounter corresponding to term " << i 
	      << ", succeeding_term_component = " << succeeding_term_component << endl;
}

				if (ith_term_component == succeeding_term_component)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_I_present: found virtual Reidemeister I move in gauss_data term " << i << endl;
					return true;
				}
			}
			else
			{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "virtual_Reidemeister_I_present:   term " << (i+j)%num_terms << " is virtual, or a shortcut crossing" << endl;
			}
		}
	}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "virtual_Reidemeister_I_present: no virtual Reidemeister I move found" << endl;
	
	return false;
}

/* remove_virtual_Reidemeister_I removes virtual Reidemeister I monogons from code_data and also virtual Reidemeister I 
   detours.  That is, consider a virtual Reidemeister I move and replace the Reidemeister I edge (the edge that leaves 
   and returns to the Reidemeister I crossing) with a virtual detour.    
   
   The function calls remove_edge_flags_from_peer_code to remove the Reidemeister I configurations and returns the number
   of components removed by remove_edge_flags_from_peer_code.
*/

int remove_virtual_Reidemeister_I(generic_code_data& code_data, vector<int>& component_flags)
{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "remove_virtual_Reidemeister_I: presented with code data ";
	write_code_data(debug,code_data);	
	debug << endl;
}
	matrix<int>& code_table = code_data.code_table;
	vector<int>& first_edge_on_component = code_data.first_edge_on_component;
	vector<int>& num_component_edges = code_data.num_component_edges;
	vector<int>& term_crossing = code_data.term_crossing;
	int num_crossings = code_data.num_crossings;

	vector<int>& shortcut_crossing = code_data.shortcut_crossing;
	bool ignore_shortcut = false;
	int head_semi_arc = -1;
	
	int num_components_removed = 0;  //the return variable
	
	if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && shortcut_crossing.size())
	{
		ignore_shortcut = true;
		
		if (code_table[generic_code_data::table::LABEL][code_data.head] == generic_code_data::POSITIVE)
			head_semi_arc = code_table[generic_code_data::table::OPEER][code_data.head];
		else if (code_table[generic_code_data::table::LABEL][code_data.head] == generic_code_data::NEGATIVE)
			head_semi_arc = 2*code_data.head;
		else
		{
			cout << "\nError! Function remove_virtual_Reidemeister_I presented with a knotoid code whose first shortcut crossing is virtual." << endl;
			exit(0);
		}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I: head_semi_arc = " << head_semi_arc << endl;			
	}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I: processing code data with " << num_crossings << " crossings" << endl;

	/* We record which edges terminate at a crossing involved in a Reidemeister I move */
	vector<int> RI_edge_flag(2*num_crossings); 

	int Reidemeister_I_count = 0;

	/* we may arrive at a Reidemeister I crossing for the first time on an odd edge or an even edge */
	for (int start=0; start<2*num_crossings; start++)
	{	
		int crossing = term_crossing[start];
		int start_label = code_table[generic_code_data::table::LABEL][crossing];
		bool head_semi_arc_in_Reidemeister_I_loop = false;
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I: starting edge " << start << ", crossing " << crossing << ", label = " << start_label << endl;
	
		/* Check added for parity bracket as part of braid programme to force detection of virtual 
		   or shortcut Reidemeister I configurations only */
		if (start_label != generic_code_data::VIRTUAL && !(ignore_shortcut && code_data.shortcut_crossing[crossing]))
			continue;

		/* if the start edge lies within a Reidemeister I configuration that we have already considered, e.g edge 13 in 
		   [-11^ -13 7 9 5 3 -1]/+ + * - - - -, then we shall not look for further loops from this edge
		*/
		if (RI_edge_flag[start] == 1)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I: starting edge already identified as part of a Reidemeister I loop" << endl;
			continue;
		}
		
		int peer;
		int component;
		int peer_component;
		
		if (start%2) // arriving on odd edge
		{
			peer = code_table[generic_code_data::table::EPEER][(start-1)/2];
			component = code_table[generic_code_data::table::COMPONENT][(start-1)/2];
			peer_component = code_table[generic_code_data::table::COMPONENT][peer/2]; 
		}
		else
		{
			peer = code_table[generic_code_data::table::OPEER][start/2];
			component = code_table[generic_code_data::table::COMPONENT][start/2];
			peer_component = code_table[generic_code_data::table::COMPONENT][(peer-1)/2]; 
		}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_virtual_Reidemeister_I:   start = " << start << ", peer = " << peer << ", start component = " << component 
	      << ", peer_component = " << peer_component << endl;	
}		
		/* look around the component from edge i to see if there's a Reidemeister I detour */					
		if (component == peer_component)
		{
			vector<int> virtual_crossing_flag(num_crossings);
			
			for (int i=1; i< num_component_edges[component]; i++)
			{
				int edge = (start+i - first_edge_on_component[component])%num_component_edges[component] + first_edge_on_component[component];

				if (ignore_shortcut && edge == head_semi_arc)
				{
					head_semi_arc_in_Reidemeister_I_loop = true;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I:   knotoid head_semi_arc " << head_semi_arc << " encountered" << endl;						
				}
				
				/* we might have a component that is (modulo detours) a simple closed curve with either just one non-virtual 
				   crossing and one virtual crossing (involving a second component) , or with two virtual crossings.  In these 
				   cases we will trace from start and arrive back on the same edge, which is not a Reidemeister I condition.
				*/
				if (term_crossing[edge] == crossing && edge != start)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I:   edge " << edge << " arrives back at crossing " << crossing << endl;

					/* we may have been given the peer code of a knot-like (non-pure) knotoid where the Reidemeister I 
					   crossing is not virtual.  In this case we cannot remove the Reidemeister I configuration.
					   If the Reidemeister I crossing is either virtual or a shortcut crossing, and the Reidemeister I loop
					   includes the leg or the head of the knotoid, then a non-shortcut arc involved in crossing may
					   be retracted to remove the crosing, since any other crossings it encounters in the loop are themselves
					   either virtual or shortcut crossings.
					*/
					if (ignore_shortcut && start_label != generic_code_data::VIRTUAL && !code_data.shortcut_crossing[crossing] &&
					    component == 0 && (head_semi_arc_in_Reidemeister_I_loop || edge < start)) // latter implies leg is in the Reidemeister						
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I:   knotoid leg or head_semi_arc in Reidemeister I loop, not a valid Reidemeister I configuration" << endl;											
					}
					else
					{
						/* this case includes the simple Reidemeister I move, where virtual_crossing_flag is empty */
						Reidemeister_I_count++;
					
						RI_edge_flag[code_table[generic_code_data::table::ODD_TERMINATING][crossing]] = 1;
						RI_edge_flag[code_table[generic_code_data::table::EVEN_TERMINATING][crossing]] = 1;
						
						for (int i=0; i< num_crossings; i++)
						{
							if (virtual_crossing_flag[i] == 1)
							{
								RI_edge_flag[code_table[generic_code_data::table::ODD_TERMINATING][i]] = 1;
								RI_edge_flag[code_table[generic_code_data::table::EVEN_TERMINATING][i]] = 1;
							}
						}
					
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_virtual_Reidemeister_I:   Reidemeister_I configuration found, Reidemeister_I_count = "
	      << Reidemeister_I_count << endl;
}

						if (edge > start) // we might have cycled around the component 
						{
							start = edge;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I:   start adjusted to " << start << endl;

						}
					}
					
					break;
				}
				else if (code_table[generic_code_data::table::LABEL][term_crossing[edge]] == generic_code_data::VIRTUAL)
				{
					virtual_crossing_flag[term_crossing[edge]] = 1;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I:   edge " << edge << " terminates at a virtual crossing" << endl;
				}
				else if (ignore_shortcut && shortcut_crossing[term_crossing[edge]])
				{
					virtual_crossing_flag[term_crossing[edge]] = 1;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I:   edge " << edge << " terminates at a shortcut crossing" << endl;
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I:   edge " << edge << " terminates at a different, non-virtual, crossing" << endl;
	
					/* note that we cannot advance start to this non-virtual crossing, since there may be virtual Reidemeister I
					   detours within the virtual crossings we've just passed through 
					*/
					break;
				}
			}
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I:   crossing involves two components, moving on"<< endl;
		}		
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I: found " << Reidemeister_I_count << " Reidemeister I moves" << endl;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_virtual_Reidemeister_I: RI_edge_flags ";
	for (int i=0; i< 2*num_crossings; i++)
		debug << RI_edge_flag[i] << ' ';
	debug << endl;
}

	if (Reidemeister_I_count > 0)
	{
		num_components_removed = remove_edge_flags_from_peer_code(code_data, RI_edge_flag,component_flags);
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I: removing_edge_flags_from_peer_code removes " << num_components_removed << " components" << endl;
			
		/* recurse in case there are multiple twists to be removed */
		if (code_data.num_crossings > 0)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I: recursing..." << endl;
			num_components_removed += remove_virtual_Reidemeister_I(code_data,component_flags);
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I: reduced peer code to a set of disjoint simple closed curves" << endl;
		}
	}
	else
	{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "remove_virtual_Reidemeister_I: no Reidemeister I moves detected." << endl;	
	}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_virtual_Reidemeister_I: returning num_components_removed = " << num_components_removed << endl;	
	return num_components_removed;	
}


/* remove_edge_flags_from_peer_code renumbers a peer code omitting those edges identified by
   edge_flag.  It was written for use by remove_Reidemeister_1 and remove_Reidemeister_II  
   but is also used when removing a component from a peer code.
   
   Note that this function can produce the peer code of a disconnected diagram.  For example, 
   the peer code of two virtually interlinked trefoils is [5 11 7 1, 13 3 15 9]/+ * + + + * + +
   and if edges 2,3,10, and 11 are removed, this function produces the code 
   [3 5 1, 9 11 7]/+ + + + + +.
   
   The function can also remove from the diagram simple closed curve components that are connected 
   to the other components only by crossings that are being removed (a simple example is a trefoil
   and an unknot linked via a Reidemeister II move).  
   
   The function returns the difference between the number of components supplied to the function and the
   number returned.
*/
int remove_edge_flags_from_peer_code(generic_code_data& code_data, vector<int>& edge_flag, vector<int>& component_flags)
{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_edge_flags_from_peer_code: presented with code data ";
	write_code_data(debug,code_data);	
	debug << endl;
	debug << "remove_edge_flags_from_peer_code: component_flags: ";
	for (unsigned int j=0; j< component_flags.size(); j++)
		debug << component_flags[j] << ' ';
	debug << endl;
}

	int num_crossings = code_data.num_crossings;
	int num_components = code_data.num_components;
	matrix<int>& code_table = code_data.code_table;
	vector<int>& num_component_edges = code_data.num_component_edges;
	vector<int>& first_edge_on_component = code_data.first_edge_on_component;

	vector<int>& shortcut_crossing = code_data.shortcut_crossing;
		
	bool pure_knotoid_code_data = false;
	bool track_zig_zag_counts;
	if (code_data.zig_zag_count.numrows() !=0)
		track_zig_zag_counts = true;
	else
		track_zig_zag_counts = false;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: track_zig_zag_counts = " << track_zig_zag_counts << endl;			
		
	
	int head_semi_arc;
    if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && code_data.head !=-1 && shortcut_crossing.size())
    {	
		pure_knotoid_code_data = true;
		
		if (code_table[generic_code_data::table::LABEL][code_data.head] == generic_code_data::POSITIVE)
			head_semi_arc = code_table[generic_code_data::table::OPEER][code_data.head];
		else if (code_table[generic_code_data::table::LABEL][code_data.head] == generic_code_data::NEGATIVE)
			head_semi_arc = 2*code_data.head;
		else
		{
			cout << "\nError! Function remove_edge_flags_from_peer_code presented with a knotoid code whose first shortcut crossing is virtual." << endl;
			exit(0);
		}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: head_semi_arc = " << head_semi_arc << endl;			
	}
	
	/* We shall re-number the diagram missing out those edges that terminate at crossings involved
	   in the Reidemeister I/II moves we've identified in edge_flag.  We write the new edge labels into a 
	   2 x num_crossings matrix, new_edge_labels, that mimics the terminating edge data in the code_table 
	   in such a way that, if old edge e is renumbered e' we write e' in the same position in new_edge_labels
	   as e appears in the terminating edge rows of the code data.  This enables us to read the new peer 
	   information from new_edge_labels.
	   
	   When tracking zig-zag counts, this approach means that the component_zig_zag_count entries are correctly
	   aligned with the new edge labels.  Note, however, that removing (non-virtual) Reidemeister II configurations 
	   may result in zig-zag counts from different Gauss-arcs being amalgamated.
	*/
	matrix<int> new_edge_labels(2,num_crossings);
	for (int i = 0; i< num_crossings; i++)
		new_edge_labels[0][i] = new_edge_labels[1][i] = -1;

	int new_edge = 0;
	int new_edge_component = 0;
	vector<int> new_first_component_edge(num_components);  // correctly initialzes new_first_component_edge[0]=0
	vector<int> new_last_component_edge(num_components);
	
	matrix<int> initial_new_zig_zag_count(2,num_crossings);
	int amalgamated_zig_zag_count = 0;
	
	/* any new head semi-arc will be non-zero */
	int new_head_semi_arc = 0;  
	
	for (int old_edge = 0; old_edge < 2*num_crossings; old_edge++)
	{
		if (pure_knotoid_code_data && old_edge == head_semi_arc)
			new_head_semi_arc = new_edge;

		int old_edge_component;
		if (old_edge%2)
			old_edge_component = code_table[generic_code_data::table::COMPONENT][(old_edge-1)/2];
		else
			old_edge_component = code_table[generic_code_data::table::COMPONENT][(old_edge)/2];
			
		if (old_edge_component != new_edge_component)
		{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "remove_edge_flags_from_peer_code:   reached end of component " << new_edge_component << endl;				 
			new_last_component_edge[new_edge_component] = new_edge-1;

			
			/*  We may have reached the end of a component and have a non-zero amalgamated_zig_zag_count since we last encountered an edge that 
				is retained.  If this is a segment component, and we're working with a pure knotoid, then we will have already included the head
				semi-arc zig-zag count in amalgamated_zig_zag_count, so need only set the head_zig_zag_count.  If we're working with a knot-type
				knotoid or a long knot, we have to amalgamate the current count with the head_zig_zag_count.
				
				In the case of a segment component we need also to update the head semi-arc zig-zag count to keep it aligned with head_zig-zag-count.
				
				If this is a loop segment, we need to amalgamate the amalgamated_zig_zag_count with that of the first edge of the component that 
				terminates at an odd crossing.  Note that in loop components, apart form a pure knotoid's head semi-arc, zig-zag counts are only 
				recorded for odd crossings; that is, not on edges terminating at a virtual or shortcut crossing.

				
			*/
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: amalgamated_zig_zag_count = " << amalgamated_zig_zag_count << " at the end of component " << new_edge_component << endl;
	
			if (code_data.immersion_character != generic_code_data::character::CLOSED && new_edge_component == 0)
			{
				/* set, or amalgamate with, the head_zig_zag_count  - see above comment */
				if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
				{
					code_data.head_zig_zag_count = amalgamated_zig_zag_count;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   pure knotoid segment component, code_data.head_zig_zag_count set to " << code_data.head_zig_zag_count << endl;
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   knot-type knotoid or long knot segment component, amalgamated with code_data.head_zig_zag_count " << code_data.head_zig_zag_count;
						code_data.head_zig_zag_count = amalgamate_zig_zag_counts(amalgamated_zig_zag_count, code_data.head_zig_zag_count);
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " to give " << code_data.head_zig_zag_count << endl;
				}
			
				/* update head semi-arc zig-zag count */
				if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
				{
					if (code_data.head == -1)
					{
						cout << "\nError! Function remove_edge_flags_from_peer_code cannot update head semi-arc zig-zag count, code_data.head = -1" << endl;
						exit(0);
					}
					
					initial_new_zig_zag_count[(head_semi_arc %2?1:0)][code_data.head] = code_data.head_zig_zag_count;
				}
			}
			else
			{				
				/* amalgamate with the count associated with the first odd crossing on the component, note that this may be a component
				   that is being removed, in which case we'll not find any edge assigned a new label and the component will be identified
				   as being removed later, when we consider all components.
				*/
				bool edge_found = false;
				
				for (int j=0; j< num_component_edges[new_edge_component]; j++)
				{
					int edge = first_edge_on_component[new_edge_component]+j;
					int crossing = code_data.term_crossing[edge];
					int row = (edge%2 == 0?0:1);
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   old edge " << edge << " terminates at crossing " << crossing << " row = " << row << endl;
					
					if (code_table[generic_code_data::table::LABEL][crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[crossing]))
					{
						if (new_edge_labels[row][crossing] != -1)
						{								
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_edge_flags_from_peer_code:   loop component, amalgamating with count " << initial_new_zig_zag_count[row][crossing] << " for old edge " 
	      << edge << ": row = " << row << " crossing = " << crossing << endl;
}
							initial_new_zig_zag_count[row][crossing] = amalgamate_zig_zag_counts(amalgamated_zig_zag_count, initial_new_zig_zag_count[row][crossing]);
										
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   count updated to " << initial_new_zig_zag_count[row][crossing] << endl;
	
							edge_found = true;
							break;
						}
					}				
				}
				
				if (!edge_found)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   loop component, no edge found, component being removed" << endl;
				}					
			}
			amalgamated_zig_zag_count = 0;
			new_edge_component++;
			new_first_component_edge[new_edge_component] = new_edge;
		}

		if (track_zig_zag_counts)
		{
			/* amalgamate any zig-zags from old_edge with the current amalgamated_zig_zag_count */
			int crossing = code_data.term_crossing[old_edge];
			int old_edge_zig_zag_count = (old_edge %2 == 1? code_data.zig_zag_count[1][crossing]:code_data.zig_zag_count[0][crossing]);
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: old_edge " << old_edge << " terminates at crossing " << crossing << ", old_edge_zig_zag_count = " << old_edge_zig_zag_count << endl;				 
			
			amalgamated_zig_zag_count = amalgamate_zig_zag_counts(amalgamated_zig_zag_count, old_edge_zig_zag_count);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   amalgamated_zig_zag_count adjusted to " << amalgamated_zig_zag_count << endl;				 
		}
		
		if (edge_flag[old_edge] == 0)
		{
			/* write new_edge into new_edge_labels at the location corresponding to old_edge's place in the terminating edge data 
			   contained in the code table.  We know which row to look in because of the polarity of old_edge;				
			   
			   if the edge terminates at an odd crossing (i.e one that is not virtual or a shortcut crossing) write the amalgamated_zig_zag_count 
			   to the corresponding location in initial_new_zig_zag_count and reset  amalgamated_zig_zag_count to zero;
			*/
			int row = (old_edge%2 ? generic_code_data::table::ODD_TERMINATING: generic_code_data::table::EVEN_TERMINATING);
			int col=0;
			for (int i=0; i< num_crossings; i++)
			{
				if (code_table[row][i] == old_edge)
				{
					col = i;
					break;
				}
			}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   old_edge " << old_edge << " lies in position " << (row==generic_code_data::table::ODD_TERMINATING? "generic_code_data::table::ODD_TERMINATING" : "generic_code_data::table::EVEN_TERMINATING") << ", " << col << ", new_edge = " << new_edge << endl;				 
				
			/* we'll write the new edges corresponding to the old even edges into the first row of new_edge_labels 
			   and the new edges corresponding to the old odd edges into the second row */
			if (row == generic_code_data::table::ODD_TERMINATING)
				row = 1;
			else
				row = 0;
				
			new_edge_labels[row][col] = new_edge;				
			
			if (track_zig_zag_counts)
			{
				int crossing = code_data.term_crossing[old_edge];
				
				if (code_table[generic_code_data::table::LABEL][crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[crossing]))
				{					
					initial_new_zig_zag_count[row][col] = amalgamated_zig_zag_count;
					amalgamated_zig_zag_count = 0;				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   new_edge " << new_edge << " initial_new_zig_zag_count[" << row << "][" << col << "] = " << initial_new_zig_zag_count[row][col] << endl;
				}
			}
			new_edge++;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   old_edge " << old_edge << " terminates at a crossing to be removed" << endl;				 	
		}	
	}

	/*  If we have reached the end of the last component and have a non-zero amalgamated_zig_zag_count since we last encountered an edge that 
	    is retained, then we need to amalgamate this count with either the head_zig_zag_count or the first edge of the component otherwise.  
		See comments above.
	*/
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: amalgamated_zig_zag_count = " << amalgamated_zig_zag_count << " at the end of component " << new_edge_component << endl;
	
	if (code_data.immersion_character != generic_code_data::character::CLOSED && new_edge_component == 0)
	{
		/* set, or amalgamate with, the head_zig_zag_count */
		if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
		{
			code_data.head_zig_zag_count = amalgamated_zig_zag_count;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   pure knotoid segment component, code_data.head_zig_zag_count set to " << code_data.head_zig_zag_count << endl;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_edge_flags_from_peer_code:   knot-type knotoid or long knot segment component, amalgamated with code_data.head_zig_zag_count " << code_data.head_zig_zag_count;
}
			code_data.head_zig_zag_count = amalgamate_zig_zag_counts(amalgamated_zig_zag_count, code_data.head_zig_zag_count);
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " to give " << code_data.head_zig_zag_count << endl;
		}
					
		if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
			initial_new_zig_zag_count[(head_semi_arc %2?1:0)][code_data.head] = code_data.head_zig_zag_count;
	}
	else
	{				
		bool edge_found = false;
		
		for (int j=0; j< num_component_edges[new_edge_component]; j++)
		{
			int edge = first_edge_on_component[new_edge_component]+j;
			int crossing = code_data.term_crossing[edge];
			int row = (edge%2 == 0?0:1);
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   old edge " << edge << " terminates at crossing " << crossing << " row = " << row << endl;
					
			if (code_table[generic_code_data::table::LABEL][crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[crossing]))
			{
				if (new_edge_labels[row][crossing] != -1)
				{	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_edge_flags_from_peer_code:   loop component, amalgamating with count " << initial_new_zig_zag_count[row][crossing] << " for old edge " 
	      << edge << ": row = " << row << " crossing = " << crossing << endl;
}
					initial_new_zig_zag_count[row][crossing] = amalgamate_zig_zag_counts(amalgamated_zig_zag_count, initial_new_zig_zag_count[row][crossing]);
										
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   count updated to " << initial_new_zig_zag_count[row][crossing] << endl;
	
					edge_found = true;
					break;
				}
			}				
		}
		
		if (!edge_found)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   loop component, no edge found, component being removed" << endl;
		}				
	}
	amalgamated_zig_zag_count = 0;
	
	/* complete the assignment of new_last_component_edge, we can do this even though we 
	   haven't added in the new edges corresponding to moving the peer detour, since
	   we may calculate it directly.
	*/
	new_last_component_edge[new_edge_component] = new_edge-1;
					
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_edge_flags_from_peer_code: new_first_component_edge ";
	for (int i=0; i< num_components; i++)
		debug << new_first_component_edge[i] << ' ';
	debug << "\nremove_edge_flags_from_peer_code: new_last_component_edge ";
	for (int i=0; i< num_components; i++)
		debug << new_last_component_edge[i] << ' ';
	debug << "\nremove_edge_flags_from_peer_code: old edge_labels: " << "\nremove_edge_flags_from_peer_code: ";
	for (int i=0; i< num_crossings; i++)
		debug << setw(3) << code_table[generic_code_data::table::EVEN_TERMINATING][i];
	debug << "\nremove_edge_flags_from_peer_code: ";
	for (int i=0; i< num_crossings; i++)
		debug << setw(3) << code_table[generic_code_data::table::ODD_TERMINATING][i];
	debug << endl;
	debug << "remove_edge_flags_from_peer_code: new_edge_labels: " << endl;
	print(new_edge_labels, debug,3,"remove_edge_flags_from_peer_code: ");
	debug << "remove_edge_flags_from_peer_code: new_head_semi_arc = " << new_head_semi_arc << endl;
	
	if (track_zig_zag_counts)
	{
		debug << "remove_edge_flags_from_peer_code: initial_new_zig_zag_count: " << endl;
		print(initial_new_zig_zag_count,debug, 3,"remove_edge_flags_from_peer_code: ");	
	}	
}

	/* new_edge now indicates the number of edges in the reduced diagram, so we can note the new number of crossings */
	int new_num_crossings = new_edge/2;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: new_num_crossings = " << new_num_crossings << endl;		

	if (new_num_crossings > 0)
	{
		/* now we know the new_num_crossings, adjust the new_head_semi_arc, in case it was set at the end of the tracing process, 
		   as is the case with [-11^ -13 7 9 5 3 -1]/+ + * - - - -
		   This may result (e.g. in the above case) in new_head_semi_arc being set to zero.  However, we may need to renumber the
		   components so need not do anything in this case here.
		   
		   check first that we have not removed the fist component, as in [-11 -13, -1 -9 -5 -7 -3^]/# * # * # # +, whereupon 
		   we are left with a non-knotoid diagram.  We can detect this condition from new_last_component_edge[0], which will be
		   set to -1 if the first component contains no new edges.
		   
		*/
		if (pure_knotoid_code_data)
		{
			if (new_last_component_edge[0] == -1)
			{
				code_data.head = -1;
				pure_knotoid_code_data = false;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: removed first component completely, cleared code_data.head and pure_knotoid_code_data." << endl;			
			}
			else
			{
				new_head_semi_arc %= (new_last_component_edge[0]+1); 
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: reset new_head_semi_arc = " << new_head_semi_arc << endl;			
			}
		}

		/* If we're dealing with a knotoid, we need to ensure that the leg of the segment component is labelled with zero
		   but it may be the case that the old edge zero is a transverse skip edge or a peer detour skip edge, in which case
		   the new label corresponding to edge zero will have changed.  If this is the case, we renumber the first component
		   so that the correct numbering of the leg semi-arc is retained.
		*/
		
		if (pure_knotoid_code_data && new_edge_labels[0][0] > 0)  // i.e not removed ( == -1) but non-zero
		{
			int component_zero_shift = new_edge_labels[0][0];
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: knotoid: need to shift numbering of component zero by " << component_zero_shift << endl;	
	
			int new_num_edges = num_component_edges[0];
			for (int i = 0; i < num_component_edges[0]; i++)
				new_num_edges -= edge_flag[i];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: knotoid: component zero new_num_edges = " << new_num_edges << endl;
					
			for (int j=0; j< num_crossings; j++)
			{
				if (new_edge_labels[0][j] < new_num_edges)
					new_edge_labels[0][j] = (new_edge_labels[0][j] - component_zero_shift + new_num_edges)%new_num_edges;
					
				if (new_edge_labels[1][j] < new_num_edges)
						new_edge_labels[1][j] = (new_edge_labels[1][j] - component_zero_shift + new_num_edges)%new_num_edges;
			}					
			
			new_head_semi_arc = (new_head_semi_arc - component_zero_shift + new_num_edges)%new_num_edges;
		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_edge_flags_from_peer_code: knotoid: after shifting component zero new_edge_labels: " << endl;
	print(new_edge_labels, debug,3,"remove_edge_flags_from_peer_code: ");
	if (new_head_semi_arc)
		debug << "remove_edge_flags_from_peer_code: knotoid: after shifting component zero new_head_semi_arc = " << new_head_semi_arc << endl;
}
		}
		   
		/* For some peer codes of links, we may now have two even or two odd new edges terminating at the same crossing.
		   
		   If this is the case, we cycle the labels on individual components so that we have even and odd edges arriving 
		   at each crossing.  	   	   	   
		*/
		align_label_numbers (new_edge_labels, num_components,new_first_component_edge,new_last_component_edge);
		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_edge_flags_from_peer_code: final new_edge_labels: " << endl;
	print(new_edge_labels, debug,3,"remove_edge_flags_from_peer_code: ");
}

		/* identify any components removed from the peer code by removing the edges.  If there is such a component, 
		   none of its edges will have been assigned a new edge label
		*/
		int new_num_components = num_components;
		for (int i=0; i< num_components; i++)
		{
			bool component_removed = true;
			for (int j=0; j< num_component_edges[i]; j++)
			{
				int edge = first_edge_on_component[i]+j;
				int crossing = code_data.term_crossing[edge];
				int row = (edge%2 == 0?0:1);
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: old edge " << edge << " terminates at crossing " << crossing << " row = " << row << endl;
	
				if (new_edge_labels[row][crossing] != -1)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   new edge label " << new_edge_labels[row][crossing] << endl;
	
					component_removed = false;
					break;
				}				
			}
			
			if (component_removed)
			{
				new_num_components--;
				
				/* we adjust the immersion character below, if necessary */
				clear_component_flag(i,component_flags);
				
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_edge_flags_from_peer_code: removing edge flags disconnects component " << i << endl;
	debug << "remove_edge_flags_from_peer_code:   component_flags updated to: ";
	for (unsigned int j=0; j< component_flags.size(); j++)
		debug << component_flags[j] << ' ';
	debug << endl;
}
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code:   removing edge flags retains component " << i << endl;
			}
		}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: new_num_components = " << new_num_components << endl;

		/* to update the generic code data we create a new code_table from the existing one,
		   write the resultant generic code data to an ostringstream and then read it back in.
		   This avoids having to re-calculate the originating and terminating edges etc..  
			
		   We must recompute the component assignment from the number of new edge labels on each
		   component, since the naming edge for a crossing may now be on a different component than 
		   before.  To do this we evaluate the new first edge on each component
		*/
		matrix<int> new_code_table (generic_code_data::table::CODE_TABLE_SIZE,new_num_crossings);
		matrix<int> new_zig_zag_count (2, new_num_crossings);
		
		for (int i=0; i< num_crossings; i++)
		{
			if (new_edge_labels[0][i] != -1) 
			{
				int odd_edge = new_edge_labels[1][i];
				int even_edge = new_edge_labels[0][i];
				int odd_zig_zag_count = 0;
				int even_zig_zag_count = 0;
				if (track_zig_zag_counts)
				{
					odd_zig_zag_count = initial_new_zig_zag_count[1][i];
					even_zig_zag_count = initial_new_zig_zag_count[0][i];
				}
				
				if (even_edge % 2 == 1)
				{
					swap(odd_edge,even_edge);
					swap(odd_zig_zag_count,even_zig_zag_count);
				}
				
				int col = even_edge/2;
				
				new_code_table[generic_code_data::table::OPEER][col] = odd_edge;
				new_code_table[generic_code_data::table::EPEER][(odd_edge-1)/2] = even_edge;
				
				if (odd_edge == new_edge_labels[0][i])
				{
					/* in the new numbering, what was an even-numbered edge is now odd-numbered, therefore
					   the type and any classical label will be reversed 
					*/
					new_code_table[generic_code_data::table::TYPE][col] = (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1?generic_code_data::TYPE2:generic_code_data::TYPE1);
					
					if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE)
						new_code_table[generic_code_data::table::LABEL][col] = generic_code_data::NEGATIVE;
					else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
						new_code_table[generic_code_data::table::LABEL][col] = generic_code_data::POSITIVE;
					else
						new_code_table[generic_code_data::table::LABEL][col] = code_table[generic_code_data::table::LABEL][i];
				}
				else
				{
					new_code_table[generic_code_data::table::LABEL][col] = code_table[generic_code_data::table::LABEL][i];
					new_code_table[generic_code_data::table::TYPE][col] = code_table[generic_code_data::table::TYPE][i];
				}
				
				new_code_table[generic_code_data::table::ODD_TERMINATING][col] = odd_edge;
				new_code_table[generic_code_data::table::EVEN_TERMINATING][col] = even_edge;
				
				if (track_zig_zag_counts)
				{
					new_zig_zag_count[1][col] = odd_zig_zag_count;
					new_zig_zag_count[0][col] = even_zig_zag_count;
				}
				
				int component = 0;
				for (int j=1; j < num_components; j++)
				{
					if (even_edge >= new_first_component_edge[j])
						component++;
					else
						break;
				}
				new_code_table[generic_code_data::table::COMPONENT][col] = component;
			}
		}
		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_edge_flags_from_peer_code: new_code_table" << endl;
	debug << "remove_edge_flags_from_peer_code:  type: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::TYPE][i] << ' ';
	debug << endl;
	debug << "remove_edge_flags_from_peer_code:  odd peer: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::OPEER][i] << ' ';
	debug << endl;
	debug << "remove_edge_flags_from_peer_code:  even peer: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::EPEER][i] << ' ';
	debug << endl;
	debug << "remove_edge_flags_from_peer_code:  odd term: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::ODD_TERMINATING][i] << ' ';
	debug << endl;
	debug << "remove_edge_flags_from_peer_code:  even term: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::EVEN_TERMINATING][i] << ' ';
	debug << endl;
	debug << "remove_edge_flags_from_peer_code:  label: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::LABEL][i] << ' ';
	debug << endl;
	debug << "remove_edge_flags_from_peer_code:  component: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::COMPONENT][i] << ' ';
	debug << endl;
	
	if (track_zig_zag_counts)
	{
		debug << "remove_edge_flags_from_peer_code: new_zig_zag_count: " << endl;
		print(new_zig_zag_count,debug, 4,"remove_edge_flags_from_peer_code: ");	
	}

}

		/* copy the new code_table to code_data and set the num_crossings.  Set the new head for knotoids.
		   It may be that removing the edge labels has resulted in new_head_semi_arc being set to zero, in 
		   which case, we shall clear the head (i.e set it to -1), to reflect the fact that we have reduced 
		   the diagram to a knot-like knotoid.  We also set code_data.immersion_character to be 
		   generic_code_data::character::KNOTOID, so that we do not perform invalid Reidemeiseter II moves 
		   involving edge zero if we are evaluating the parity bracket or parity arrow polynomial.
		*/
		code_data.code_table = new_code_table;
		code_data.num_crossings = new_num_crossings;
		code_data.num_components = new_num_components;
		
		/* Determine the new number of component edges.  We have vectors new_first_component_edge and new_last_component_edge 
		  of length the old number of components, such that  new_last_component_edge[i] == new_first_component_edge[i]-1 if 
		  the i-th component has been removed.
		*/
		vector<int> new_num_component_edges(new_num_components);
		int component = 0;
		for (int i=0; i< num_components; i++)
		{
			if (new_last_component_edge[i] == new_first_component_edge[i]-1)
				continue; // component i has been removed
				
			new_num_component_edges[component] = new_last_component_edge[i]-new_first_component_edge[i]+1;
			component++;
		}
		
		code_data.num_component_edges = new_num_component_edges;
		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_edge_flags_from_peer_code: new_num_component_edges: ";
	for (int i=0; i< new_num_components; i++)
		debug << new_num_component_edges[i] << ' ';
	debug << endl;
}		

		if (track_zig_zag_counts)
		{
			bool zig_zags_present = false;
			for (int i=0; i<new_num_crossings; i++)
			{
				if (new_zig_zag_count[0][i] != 0 || new_zig_zag_count[1][i] !=0)
				{
					zig_zags_present = true;
					break;
				}
			}
			
			if (zig_zags_present)
				code_data.zig_zag_count = new_zig_zag_count;
			else
				code_data.zig_zag_count.clear();
		}
		
		if (pure_knotoid_code_data)
		{
			if (new_head_semi_arc == 0)
			{
				code_data.head = -1;
				pure_knotoid_code_data = false;
				

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: new_head_semi_arc is zero, clearing knotoid status of code data" << endl;
	
				code_data.immersion_character = generic_code_data::character::KNOTOID;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: set code_data.immersion_character = generic_code_data::character::KNOTOID to indicate knot-type knotoid" << endl;
	
			}
			else if (new_head_semi_arc % 2)
				code_data.head = code_data.code_table[generic_code_data::table::EPEER][(new_head_semi_arc-1)/2]/2;
			else
				code_data.head = new_head_semi_arc/2;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: new peer code head = " << code_data.head << endl;
			
		}

		/* write the modified code to a string */
		ostringstream oss;
		write_code_data (oss, code_data);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_edge_flags_from_peer_code: new peer code = " << oss.str() << endl;
	
		/* read back the peer code into code_data */
		read_peer_code(code_data, oss.str());
		
		/* if we're dealing with a knotoid, we need to reset the shortcut crossings, do this by using
		   valid_knotoid_input, which will assign a new vector to code_data.shortcut crossing, thereby
		   deleting what is already there.
		*/
		if (pure_knotoid_code_data)
		{
			if (!valid_knotoid_input(code_data))
			{
				cout << "\nError!  Function remove_edge_flags_from_peer_code has produced invalid knotoid code data" << endl;
				exit(0);
			}
		}
	}
	else
	{
		for (unsigned int i=0; i< component_flags.size(); i++)
			component_flags[i] = 0;
			
		code_data.num_crossings = 0;
		code_data.num_components = 0;
		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_edge_flags_from_peer_code: reduced peer code to a set of disjoint simple closed curves" << endl;
	debug << "remove_edge_flags_from_peer_code: set code_data.num_components to zero, clearing component_flags" << endl;
}
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_edge_flags_from_peer_code: new code_data: " << endl;
	print_code_data (debug, code_data, "remove_edge_flags_from_peer_code: \t");
	debug << "remove_edge_flags_from_peer_code: removed " << num_components - code_data.num_components << " components" << endl;
}		

	return (num_components - code_data.num_components);
}

bool Reidemeister_II_labels (const int label1, const int label2)
{
	if (
		(label1 == generic_code_data::POSITIVE && label2 == generic_code_data::NEGATIVE) ||
		(label2 == generic_code_data::POSITIVE && label1 == generic_code_data::NEGATIVE) ||
		(label1 == generic_code_data::VIRTUAL && label2 == generic_code_data::VIRTUAL) ||
		(label1 == generic_code_data::FLAT && label2 == generic_code_data::FLAT)
	   )
	{
		return true;
	}
	else
	{
		return false;
	}
}


bool simple_Reidemeister_II_labels (const int label1, const int label2)
{
	if (
		(label1 == generic_code_data::POSITIVE && label2 == generic_code_data::NEGATIVE) ||
		(label1 == generic_code_data::NEGATIVE && label2 == generic_code_data::POSITIVE) ||
		(label1 == generic_code_data::VIRTUAL  && label2 == generic_code_data::VIRTUAL) ||
		(label1 == generic_code_data::FLAT     && label2 == generic_code_data::FLAT)
	   )
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Gauss_Reidemeister_II_labels(int start_edge, int end_bigon_edge_a, int end_bigon_edge_b, generic_code_data& code_data)
{
	matrix<int>& code_table = code_data.code_table;
	int component;
	
	int start_crossing = code_data.term_crossing[start_edge];

	component = code_table[generic_code_data::table::COMPONENT][(start_edge %2 ? (start_edge-1)/2: start_edge/2)];

	int end_edge;
	
	for (int i = 0; i < code_data.num_component_edges[component]; i++)
	{
		int edge = (start_edge + i - code_data.first_edge_on_component[component]) % code_data.num_component_edges[component] + code_data.first_edge_on_component[component];
		if (edge == end_bigon_edge_a || edge == end_bigon_edge_b)
		{
			end_edge = edge;
			break;
		}
	}

	int end_crossing = code_data.term_crossing[end_edge];


if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Gauss_Reidemeister_II_labels: start_edge = " << start_edge << ", component = " << component << ", end_edge = " << end_edge << endl;	
	debug << "Gauss_Reidemeister_II_labels: start_crossing = " << start_crossing << ", end_crossing = " << end_crossing << endl;	
}

	if (code_table[generic_code_data::table::LABEL][start_crossing] == generic_code_data::FLAT && code_table[generic_code_data::table::LABEL][end_crossing] == generic_code_data::FLAT)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_labels: start and end crossings both FLAT, compatible Reidemeister II labels" << endl;	
		
		return true;
	}
	else if ((code_table[generic_code_data::table::LABEL][start_crossing] == generic_code_data::POSITIVE || code_table[generic_code_data::table::LABEL][start_crossing] == generic_code_data::NEGATIVE)
	       && (code_table[generic_code_data::table::LABEL][end_crossing] == generic_code_data::POSITIVE || code_table[generic_code_data::table::LABEL][end_crossing] == generic_code_data::NEGATIVE))
	{	

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_labels: both start and end crossings are classical" << endl;	

		bool start_under_arc;
		if (start_edge % 2)
		{
			if (code_table[generic_code_data::table::LABEL][start_crossing] == generic_code_data::POSITIVE)
				start_under_arc = true;
			else
				start_under_arc = false;
		}
		else
		{
			if (code_table[generic_code_data::table::LABEL][start_crossing] == generic_code_data::POSITIVE)
				start_under_arc = false;
			else
				start_under_arc = true;
		}
	
		bool end_under_arc;
		if (end_edge % 2)
		{
			if (code_table[generic_code_data::table::LABEL][end_crossing] == generic_code_data::POSITIVE)
				end_under_arc = true;
			else
				end_under_arc = false;
		}
		else
		{
			if (code_table[generic_code_data::table::LABEL][end_crossing] == generic_code_data::POSITIVE)
				end_under_arc = false;
			else
				end_under_arc = true;
		}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_labels: start_under_arc = " << start_under_arc << ", end _under_arc = " << end_under_arc << endl;	
	
		bool Reidemeister_II = false;
		
		if ((start_under_arc && end_under_arc) || (!start_under_arc && !end_under_arc))
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_labels: compatible Reidemeister II labels "<< endl;	
			Reidemeister_II = true;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_labels: incompatible Reidemeister II labels "<< endl;	
		}
			
		return Reidemeister_II;
	}
	else 
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_labels: start and end crossings of mixed type, or both VIRTUAL, incompatible Reidemeister II labels" << endl;	
		
		return false;
	}
}



/* simple_Reidemeister_II_present looks for "ordinary" Reidemeister II moves between adjacent 
   crossings, including adjacent virtual crossings.
   
   simple_Reidemeister_II_present checks whether the diagram represented by the code data 
   it is given contains a Reidemeister Type 2 move.  This is detected 
   by having the peer of 2j or 2j-2 being 2i+1 or 2i-1, when the peer of 2i is 2j-1,
   and having the two (real) crossings with different labels, or having both as virtual or 
   flat crossings.  To see this, consider the Reidemeister II move with the two strands 
   oriented in the same direction and then with opposing directions.  For each configuration 
   of strand, consider the ith crossing (naming edge 2i) being the first crossing on one 
   strand and then the second crossing on that strand.  In all these cases the peer of 2i
   is 2j-1.
   
   2i    __2j__    2i+2              2i   _2j-1_    2i+2
      \ /      \ /  			       \ /      \ /
1)	   x        x 			2)		    x        x
	  / \_2i+1_/ \				       / \_2i+1_/ \
 2j-1             2j+1		         2j            2j-2
 
   peer of 2j is 2i+1               peer of 2j-2 is 2i+1
   i.e peer of peer_successor       i.e peer of peer_predecessor
   is successor                     is successor


  2i-1   _2j-1_   2i+1             2i-1   __2j__   2i+1
      \ /      \ /  			       \ /      \ /
3)	   x        x 			4)		    x        x
	  / \__2i__/ \				       / \__2i__/ \
 2j-2             2j		       2j+1            2j-1
 
   peer of 2j-2 is 2i-1             peer of 2j is 2i-1
   i.e peer of peer_predecessor    i.e peer of peer_successor
   is predecessor                  is predecessor

 
   The types of the two crossings will be the same, if they are different the immersion code is 
   not realizable.  Reversing both the type values in the two left hand configurations above 
   interchanges the two strands, reversing both the values in the two right hand configurations 
   reverses the orientation of the strands.

	In the case of knotoids, the function will ignore any non-virtual Reidemeister II configuation 
	involving either semi-arc zero or the semi-arc indicated by head_semi_arc.  For knotoids the 
	function detects Reidemeister II bigons where both crossings are shortcut crossings.
*/
simple_Reidemeister_II_return simple_Reidemeister_II_present(generic_code_data& code_data, bool create_flags)
{	
	matrix<int>& code_table = code_data.code_table;
	vector<int>& first_edge_on_component = code_data.first_edge_on_component;
	vector<int>& num_component_edges = code_data.num_component_edges;
	int num_crossings = code_data.num_crossings;

	simple_Reidemeister_II_return _return;
	if (create_flags)
		_return.RII_edge_flag = vector<int>(2*num_crossings);

	vector<int>& RII_edge_flag = _return.RII_edge_flag;
	bool& Reidemeister_II_found = _return.Reidemeister_II_found;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "simple_Reidemeister_II_present: presented with code data ";
	write_code_data(debug,code_data);	
	debug << endl;
}

	bool pure_knotoid_code_data = false;
	int head_semi_arc = -1;
	
	if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && code_data.shortcut_crossing.size())
	{
		if (code_table[generic_code_data::table::LABEL][code_data.head] == generic_code_data::POSITIVE)
			head_semi_arc = code_table[generic_code_data::table::OPEER][code_data.head];
		else if (code_table[generic_code_data::table::LABEL][code_data.head] == generic_code_data::NEGATIVE)
			head_semi_arc = 2*code_data.head;
		else
		{
			cout << "\nError! Function simple_Reidemeister_II_present presented with a knotoid code whose first shortcut crossing is virtual." << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "simple_Reidemeister_II_present: head_semi_arc = " << head_semi_arc << endl;		
	
		pure_knotoid_code_data = true;
	}
	


	/* We stop as soon as we've found one RII move to avoid situations where we have what looks like two RII moves that share a
	   middle crossing.  We only want to perform one of these moves.
	*/

	for (int i=0; i<num_crossings  && !Reidemeister_II_found ; i++)
	{	
		int peer = code_table[generic_code_data::table::OPEER][i];
		int component = code_table[generic_code_data::table::COMPONENT][i];
		int successor = (2*i+1 - first_edge_on_component[component])%num_component_edges[component] + first_edge_on_component[component];
		int predecessor = (2*i-1 - first_edge_on_component[component]+num_component_edges[component])%
		                     num_component_edges[component] + first_edge_on_component[component];	
		
		/* the first edge on a component is always even, so peer-1 is on the same component as peer */
		int peer_component = code_table[generic_code_data::table::COMPONENT][(peer-1)/2]; 

		int peer_successor = (peer+1 - first_edge_on_component[peer_component])%
		                     num_component_edges[peer_component] + first_edge_on_component[peer_component];
		int peer_predecessor = (peer-1 - first_edge_on_component[peer_component]+num_component_edges[peer_component])%
		                     num_component_edges[peer_component] + first_edge_on_component[peer_component];

		/* the j crossing is the crossing number at which the peer_successor terminates */
		int j = peer_successor/2;
	
		/* the j1 crossing is the crossing number at which the peer_predecessor terminates */
		int j1 = peer_predecessor/2;

		
		if (code_table[generic_code_data::table::OPEER][j] == successor && simple_Reidemeister_II_labels(code_table[generic_code_data::table::LABEL][i],code_table[generic_code_data::table::LABEL][j]))
	    {
			/* 1) peer of 2j is 2i+1; i.e peer of peer_successor is successor */
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "simple_Reidemeister_II_present:  i = " << i << ", naming edge = " << 2*i << ", peer = " << peer 
	      << ", peer of peer_sucessor=" << peer_successor << " is sucessor=" << successor << ", " << endl;
}
			if (pure_knotoid_code_data && code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL &&
			    (successor == 0 || successor == head_semi_arc || peer_successor == 0 || peer_successor == head_semi_arc))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "simple_Reidemeister_II_present: case 1, bigon bounded by successor and peer_successor includes knotoid leg (0) or head (" << head_semi_arc << ")" << endl;
			}
			else if (code_data.immersion_character == generic_code_data::character::LONG_KNOT && 
			         code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL && (successor == 0 || peer_successor == 0 ))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "simple_Reidemeister_II_present: case 1, bigon bounded by successor and peer_successor includes long knot ends" << endl;
			}
			else if (code_data.immersion_character == generic_code_data::character::KNOTOID && 
			         code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL && (successor == 0 || peer_successor == 0 ))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	if (code_data.type == generic_code_data::peer_code)
		debug << "simple_Reidemeister_II_present: case 1, bigon bounded by successor and peer_successor includes knot-type knotoid leg and head" << endl;
	else
		debug << "simple_Reidemeister_II_present: case 1, bigon bounded by successor and peer_successor includes knot-type knotoid leg or head" << endl;
}
			}
			else
			{
				Reidemeister_II_found = true;
	
				if (create_flags)
					RII_edge_flag[2*i] = RII_edge_flag[peer] = RII_edge_flag[peer_successor] = RII_edge_flag[successor] =1;
			}
		}
		else if (code_table[generic_code_data::table::OPEER][j] == predecessor && simple_Reidemeister_II_labels(code_table[generic_code_data::table::LABEL][i],code_table[generic_code_data::table::LABEL][j]))
	    {
			/* 4) peer of 2j is 2i-1 i.e peer of peer_successor is predecessor */
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "simple_Reidemeister_II_present:  i = " << i << ", naming edge = " << 2*i << ", peer = " << peer 
	      << ", peer of peer_successor=" << peer_successor << " is predecessor=" << predecessor << ", " << endl;
}
			if (pure_knotoid_code_data && code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL &&
			    (2*i == 0 || 2*i == head_semi_arc || peer_successor == 0 || peer_successor == head_semi_arc))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "simple_Reidemeister_II_present: case 4, bigon bounded by naming edge and peer_successor includes knotoid leg (0) or head (" << head_semi_arc << ")" << endl;
			}
			else if (code_data.immersion_character == generic_code_data::character::LONG_KNOT && 
			         code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL && (2*i == 0 || peer_successor == 0))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "simple_Reidemeister_II_present: case 4, bigon bounded by naming edge and peer_successor includes long knot ends" << endl;
			}
			else if (code_data.immersion_character == generic_code_data::character::KNOTOID && 
			         code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL && (2*i == 0 || peer_successor == 0))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	if (code_data.type == generic_code_data::peer_code)
		debug << "simple_Reidemeister_II_present: case 4, bigon bounded by naming edge and peer_successor includes knot-type knotoid leg and head" << endl;
	else
		debug << "simple_Reidemeister_II_present: case 4, bigon bounded by naming edge and peer_successor includes knot-type knotoid leg or head" << endl;
}
			}
			else
			{
				Reidemeister_II_found = true;
	
				if (create_flags)
					RII_edge_flag[2*i] = RII_edge_flag[peer] = RII_edge_flag[peer_successor] = RII_edge_flag[predecessor] =1;
			}
		}
		else if (code_table[generic_code_data::table::OPEER][j1] == successor && simple_Reidemeister_II_labels(code_table[generic_code_data::table::LABEL][i],code_table[generic_code_data::table::LABEL][j1]))
		{
			/* 2) peer of 2j-2 is 2i+1 i.e peer of peer_predecessor is successor */
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "simple_Reidemeister_II_present:  i = " << i << ", naming edge = " << 2*i << ", peer = " << peer 
	      << ", peer of peer_predecessor=" << peer_predecessor << " is successor=" << successor << ", " << endl;
}
			if (pure_knotoid_code_data && code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL &&
			    (peer == 0 || peer == head_semi_arc || successor == 0 || successor == head_semi_arc))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "simple_Reidemeister_II_present: case 2, bigon bounded by peer and successor includes knotoid  leg (0) or head (" << head_semi_arc << ")" << endl;
			}
			else if (code_data.immersion_character == generic_code_data::character::LONG_KNOT && 
			         code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL && (peer == 0 || successor == 0))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "simple_Reidemeister_II_present: case 2, bigon bounded by peer and successor includes long knot ends" << endl;
			}
			else if (code_data.immersion_character == generic_code_data::character::KNOTOID && 
			         code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL && (peer == 0 || successor == 0))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	if (code_data.type == generic_code_data::peer_code)
		debug << "simple_Reidemeister_II_present: case 2, bigon bounded by peer and successor includes knot-type knotoid leg and head" << endl;
	else
		debug << "simple_Reidemeister_II_present: case 2, bigon bounded by peer and successor includes knot-type knotoid leg or head" << endl;
}
			}
			else
			{
				Reidemeister_II_found = true;
	
				if (create_flags)
					RII_edge_flag[2*i] = RII_edge_flag[peer] = RII_edge_flag[peer_predecessor] = RII_edge_flag[successor] =1;
			}
		}
		else if (code_table[generic_code_data::table::OPEER][j1] == predecessor && simple_Reidemeister_II_labels(code_table[generic_code_data::table::LABEL][i],code_table[generic_code_data::table::LABEL][j1]))
		{
			/* 3) peer of 2j-2 is 2i-1 i.e peer of peer_predecessor is predecessor */
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "simple_Reidemeister_II_present:  i = " << i << ", naming edge = " << 2*i << ", peer = " << peer 
	      << ", peer of peer_predecessor=" << peer_predecessor << " is predecessor=" << predecessor << ", " << endl;
}

			if (pure_knotoid_code_data && code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL &&
			    (2*i == 0 || 2*i == head_semi_arc || peer == 0 || peer == head_semi_arc))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "simple_Reidemeister_II_present: case 3, bigon bounded by naming edge and peer includes knotoid leg (0) or head (" << head_semi_arc << ")" << endl;
			}
			else if (code_data.immersion_character == generic_code_data::character::LONG_KNOT && 
			         code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL && (2*i == 0 || peer == 0))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "simple_Reidemeister_II_present: case 3, bigon bounded by naming edge and peer includes long knot ends" << endl;
			}
			else if (code_data.immersion_character == generic_code_data::character::KNOTOID && 
			         code_table[generic_code_data::table::LABEL][i] != generic_code_data::VIRTUAL && (2*i == 0 || peer == 0))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	if (code_data.type == generic_code_data::peer_code)
		debug << "simple_Reidemeister_II_present: case 3, bigon bounded by naming edge and peer includes knot-type knotoid leg and head" << endl;
	else
		debug << "simple_Reidemeister_II_present: case 3, bigon bounded by naming edge and peer includes knot-type knotoid leg or head" << endl;
}
			}
			else
			{
				Reidemeister_II_found = true;
	
				if (create_flags)
					RII_edge_flag[2*i] = RII_edge_flag[peer] = RII_edge_flag[peer_predecessor] = RII_edge_flag[predecessor] =1;
			}
		}
	}
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "simple_Reidemeister_II_present: Reidemeister_II_found = " << Reidemeister_II_found << endl;
	
	return _return;
}

/* Gauss_Reidemeister_II_present looks for Reidemeister II configurations between non-virtual crossings that can be detected
   from Gauss orientation data.  This includes Reidemeister II bigons whose edges involve only virtual crossings or, in the 
   case of knotoids, shortcut crossings.  The edges involved in such crossings may be moved to leave a bigon that may be 
   removed by a Reidemeister II move.
   
   The function also checks for components containing no classical crossings, even though these are examples of virtual
   Reidemeister II configurations.  However, since such components are floating components, removing them always 
   simplifies the diagram and so we always want to detect them, so do so here as it is simple to do and essentially comes
   for free.
   
   The function checks for one of the 2-reducible configurations in a Gauss code:
             ...LiRj...RiLj..., or ...RiLj...LiRj... (forwards conditions), 
             ...LiRj...LjRi..., or ...RiLj...RjLi... (backwards conditions)
   if it finds one of these conditions the function sets start to be the terminating edge of the initial Li or 
   Ri in the 2-reducible configuration it finds, and sets forwards_from_peer to indicate whether it found one of
   the forwards conditions or backwards conditions.
   
   The search for each consecutive pair is cyclic within the corresponding component of code_data.

   Since non-extended gauss_orientation_data, as used by Gauss_Reidemeister_II_present, ignores knotoid shortcut crossings 
   the function does not detect Reidemeister II bigons in the underlying immersion where one or both of the crossings is
   a shortcut crossing.  In the case of one crossing being a shortcut crossing and the other a classical crossing, we 
   cannot undo the bigon, so not being able to detect this case is not an issue.  The case where both crossings are 
   shortcut crossings is dealt with by the function virtual_Reidemeister_II_plus.
   
   If the head and leg of a knotoid lie within the boundary of the bigon, so the shortcut is a subset of one of the detour
   paths then, again, we cannot undo the Reidemeister II bigon in the underlying immersion.  In this case we have to ignore 
   what looks like a Reidemeister II detour configuration if edge zero, and consequently the crossing corresponding to the 
   first term of the Gauss code, lies in the boundary of the bigon.  

   The start and rendezvous crossings are all one of the following cases, depending on whether the detour is forwards
   or backwards from the peer of the start edge.  Here the labels are on the terminating edges at the crossing.
   
         Forwards:                             Backwards:                        Backwards:
     1)            2)                       3)            4)              5)=3)           6) = 4)
   R 2i          R 2j-1                                               L 2i   R 2j-1    L 2j-1  R 2i
        \ /           \ /  			       \ /           \ /              \ /               \ /
	     x             x 					x             x                x                 x
 	    / \           / \  				   / \           / \              / \               / \
   L 2j-1        L 2i		          R 2j-1  L 2i   R 2i   L 2j-1
   
   At the start crossing of a forward detour we always want to check the odd and even originating edges
   At the end crossing of a forward detour we always want to check the odd and even terminating edges
   
   At the start of a backwards detour
    - if start_edge is RIGHT and crossing is TYPE1 we want to check even terminating and even originating
    - if start_edge is RIGHT and crossing is TYPE2 we want to check odd terminating and odd originating
    - if start_edge is LEFT and crossing is TYPE2 we want to check even terminating and even originating
    - if start_edge is LEFT and crossing is TYPE1 we want to check odd terminating and odd originating
   
   At the end of a backwards detour
    - if end_edge is RIGHT and crossing is TYPE2 we want to check even terminating and even originating
    - if end_edge is RIGHT and crossing is TYPE1 we want to check odd terminating and odd originating
    - if end_edge is LEFT and crossing is TYPE1 we want to check even terminating and even originating
    - if end_edge is LEFT and crossing is TYPE2 we want to check odd terminating and odd originating
*/
Gauss_Reidemeister_II_return Gauss_Reidemeister_II_present(generic_code_data& code_data, gauss_orientation_data& gauss_data)
{
	Gauss_Reidemeister_II_return _return;	
	int& start_edge = _return.start_edge;
	bool& forwards_from_peer = _return.forwards_from_peer;
	bool&Reidemeister_II_found = _return.Reidemeister_II_found; 
	
	matrix<int>& code_table = code_data.code_table;

	if (code_data.type == generic_code_data::gauss_code)
	{
		cout << "Error!  Gauss_Reidemeister_II_present called with a Gauss code not a Peer code" << endl;
		exit (0);
	}
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Gauss_Reidemeister_II_present: gauss_data = ";
	write_gauss_data (gauss_data, debug);
	debug << endl;
	print_gauss_data (gauss_data, debug, "Gauss_Reidemeister_II_present: ");
}	

	int num_components = gauss_data.num_components;
	int num_terms = gauss_data.num_terms;
	matrix<int>& orientation_matrix = gauss_data.orientation_matrix;
	vector<int> start_of_component = gauss_data.start_of_component;
	vector<int> num_terms_in_component = gauss_data.num_terms_in_component;
	
	
	/* Check first that there is not a component of code_data containing only virtual crossings.  This is
	   a degenerate case of a Reidemeister_II detour, equivalent to a link diagram with a floating unknotted
	   component.  It is detected in the gauss_data as a component involving no terms of the gauss_data.
	   
	   Note, this check is also performed by Gauss_virtual_detour_present, but it is duplicated here since
	   remove_Reidemeister_II only calls this function and we want to catch this case in that instance.
	*/
	for (int i=0; i< num_components; i++)
	{
		if (num_terms_in_component[i] == 0)
		{		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present: floating, unknotted component detected" << endl;
			_return.floating_component = i+1;
			return _return; // return with Reidemeister_II_found set false, since here we have a floating component.
		}
	}
		
	int start_peer;  // peer of start_edge
	int start_gauss_crossing;
	int start_immersion_crossing;
	int start_edge_parity; // L or R of start_edge at crossing i first LiRj or RiLj sequence
	int end_edge_parity;  // L or R of end_edge at crossing j first LiRj or RiLj sequence
	int rendezvous_gauss_crossing;
	int rendezvous_immersion_crossing;
	int rendezvous_start_term;  // used for debugging
	int start_bigon_edge_a;
	int start_bigon_edge_b;
	int end_bigon_edge_a;
	int end_bigon_edge_b;
	
	int first_component = 0;
	
	for (int i=0; i< num_terms-1; i++)
	{
		if (i == start_of_component[first_component]+num_terms_in_component[first_component])
			first_component++;
		
		int successor = (i+1-start_of_component[first_component])%num_terms_in_component[first_component] + start_of_component[first_component];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present: term " << i << ", first_component " << first_component << ", succeeding term " << successor << endl;
		
		if (orientation_matrix[0][i] != orientation_matrix[0][successor] && orientation_matrix[1][i] != orientation_matrix[1][successor]) // LiRj or RiLj sequence
		{
			start_edge_parity = orientation_matrix[0][i];
			end_edge_parity = orientation_matrix[0][successor];
			start_gauss_crossing = orientation_matrix[1][i]-1;
			start_immersion_crossing = gauss_data.immersion_crossing[start_gauss_crossing];
			rendezvous_gauss_crossing = orientation_matrix[1][successor]-1;
			rendezvous_immersion_crossing = gauss_data.immersion_crossing[rendezvous_gauss_crossing];
						
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Gauss_Reidemeister_II_present: potential RII sequence found at term " << i << endl;
	debug << "Gauss_Reidemeister_II_present:   start_gauss_crossing " << start_gauss_crossing+1 << ", start_immersion_crossing = " << start_immersion_crossing << endl;
	debug << "Gauss_Reidemeister_II_present:   rendezvous_gauss_crossing " << rendezvous_gauss_crossing+1 << ", rendezvous_immersion_crossing = " << rendezvous_immersion_crossing << endl;
	debug << "Gauss_Reidemeister_II_present:   start label = " << code_table[generic_code_data::table::LABEL][start_immersion_crossing] << ", rendezvous label = " << code_table[generic_code_data::table::LABEL][rendezvous_immersion_crossing] << endl;
}

			/* we identify the start_edge and peer_edge by the nature of the first term at the start of a potential RII sequence.
			   Note that the identification of the start_edge as odd or even is independent of whether a Reidemeister II configuration
			   is forwards or backwards.
			*/
			
			if ((orientation_matrix[0][i] == gauss_orientation_data::LEFT && code_table[generic_code_data::table::TYPE][start_immersion_crossing] == generic_code_data::TYPE1) ||
			    (orientation_matrix[0][i] == gauss_orientation_data::RIGHT && code_table[generic_code_data::table::TYPE][start_immersion_crossing] == generic_code_data::TYPE2))
			{
				start_edge = code_table[generic_code_data::table::EVEN_TERMINATING][start_immersion_crossing];
				start_peer = code_table[generic_code_data::table::ODD_TERMINATING][start_immersion_crossing];
			}
			else
			{
				start_edge = code_table[generic_code_data::table::ODD_TERMINATING][start_immersion_crossing];
				start_peer = code_table[generic_code_data::table::EVEN_TERMINATING][start_immersion_crossing];
			}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present:   start_edge = " << start_edge << ", start_peer = " << start_peer << endl;

			/* We cannot determine the compatibility of the crossings involved in a potential Reidemeister II configuration until 
			   we know the bigon end edges, since there may be an odd or even number of crossings in the immersion between the 
			   start and rendezvous Gauss crossings.
			*/
	
			/* note that we don't need to check the LEFT or RIGHT status of the second occurrence for each crossing, 
			   since they will always be the opposite of the first occurrence
			*/
			int second_component = first_component;
			
			/* if we have found a first pair of a potential Reidemeister II detour at the end and beginning of a
			   component, we start looking for the second pair at the first term of the next component.  Otherwise
			   we start looking at the term after the (consecutive) pair of terms.
			   
			*/
			int second_start_term = (successor == i+1? i+2: i+1);  // if i+1 or i+2 == num_terms the following loop doesn't execute
			for (int j=second_start_term; j<num_terms; j++)
			{
				if (j == start_of_component[second_component]+num_terms_in_component[second_component])
					second_component++;

				successor = (j+1-start_of_component[second_component])%num_terms_in_component[second_component] + start_of_component[second_component];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present:   term " << j << ", second_component " << second_component << ", succeeding term " << successor << endl;
					
				/* In orientation data such as R1 L2 L3 L1 R3 R2, it is possible that the successor of j(=5) is the starting point i(=0) but there is
				   no Reidemeister II detour present; that is, i, i-successor, j, and j-successor must be distinct.

				   is              i         j js					   
				   R1 R2 L2 R3 R4 L5, L3 L4 R5 L1

				   We need only check that j-successor is distinct from i, since i, i-successor and j will be distinct by definition 
				*/
				
				if (successor != i)
				{
					if (orientation_matrix[1][j]-1 == start_gauss_crossing)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present:     term " << j << " involves start_gauss_crossing" << endl;
	
						if (orientation_matrix[1][successor]-1 == rendezvous_gauss_crossing)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present:     succeeding term " << successor << " involves rendezvous_gauss_crossing" << endl;

							start_bigon_edge_a = code_table[generic_code_data::table::EVEN_ORIGINATING][start_immersion_crossing];
							start_bigon_edge_b = code_table[generic_code_data::table::ODD_ORIGINATING][start_immersion_crossing];
							end_bigon_edge_a = code_table[generic_code_data::table::EVEN_TERMINATING][rendezvous_immersion_crossing];
							end_bigon_edge_b = code_table[generic_code_data::table::ODD_TERMINATING][rendezvous_immersion_crossing];

							forwards_from_peer = true;
							rendezvous_start_term = j;

							if ((Reidemeister_II_found = Gauss_Reidemeister_II_labels(start_edge, end_bigon_edge_a, end_bigon_edge_b, code_data)))
							{

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Gauss_Reidemeister_II_present: Reidemeister_II configuration found, with rendezvous starting at term " << j << ", moving forwards from start_peer " << endl;
	debug << "Gauss_Reidemeister_II_present: start_bigon_edges " << start_bigon_edge_a << " and " << start_bigon_edge_b << endl;
	debug << "Gauss_Reidemeister_II_present: end_bigon_edges " << end_bigon_edge_a << " and " << end_bigon_edge_b << endl;
}
								break;
							}
						}
						else
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present:     succeeding term " << successor << " does not involve rendezvous_gauss_crossing" << endl;
						}
					}
					else if (orientation_matrix[1][j]-1 == rendezvous_gauss_crossing)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present:     term " << j << " involves rendezvous_gauss_crossing" << endl;
	
						if (orientation_matrix[1][successor]-1 == start_gauss_crossing)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present:     succeeding term " << successor << " involves start_gauss_crossing" << endl;

							if ((start_edge_parity == gauss_orientation_data::RIGHT && code_table[generic_code_data::table::TYPE][start_immersion_crossing] == generic_code_data::TYPE1) ||
								(start_edge_parity == gauss_orientation_data::LEFT && code_table[generic_code_data::table::TYPE][start_immersion_crossing] == generic_code_data::TYPE2))
							{
								start_bigon_edge_a = code_table[generic_code_data::table::EVEN_TERMINATING][start_immersion_crossing];
								start_bigon_edge_b = code_table[generic_code_data::table::EVEN_ORIGINATING][start_immersion_crossing];
							}
							else
							{
								start_bigon_edge_a = code_table[generic_code_data::table::ODD_TERMINATING][start_immersion_crossing];
								start_bigon_edge_b = code_table[generic_code_data::table::ODD_ORIGINATING][start_immersion_crossing];
							}
			
							if ((end_edge_parity == gauss_orientation_data::RIGHT && code_table[generic_code_data::table::TYPE][rendezvous_immersion_crossing] == generic_code_data::TYPE2) ||
								(end_edge_parity == gauss_orientation_data::LEFT && code_table[generic_code_data::table::TYPE][rendezvous_immersion_crossing] == generic_code_data::TYPE1))
							{
								end_bigon_edge_a = code_table[generic_code_data::table::EVEN_TERMINATING][rendezvous_immersion_crossing];
								end_bigon_edge_b = code_table[generic_code_data::table::EVEN_ORIGINATING][rendezvous_immersion_crossing];
							}
							else
							{
								end_bigon_edge_a = code_table[generic_code_data::table::ODD_TERMINATING][rendezvous_immersion_crossing];
								end_bigon_edge_b = code_table[generic_code_data::table::ODD_ORIGINATING][rendezvous_immersion_crossing];
							}

							forwards_from_peer = false;						    
							rendezvous_start_term = j;

							if ((Reidemeister_II_found = Gauss_Reidemeister_II_labels(start_edge, end_bigon_edge_a, end_bigon_edge_b, code_data)))
							{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Gauss_Reidemeister_II_present: Reidemeister_II configuration found, with rendezvous starting at term " << j << ", moving backwards from start_peer " << endl;
	debug << "Gauss_Reidemeister_II_present: start_bigon_edges " << start_bigon_edge_a << " and " << start_bigon_edge_b << endl;
	debug << "Gauss_Reidemeister_II_present: end_bigon_edges " << end_bigon_edge_a << " and " << end_bigon_edge_b << endl;
}
								break;
							}
						}
						else
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present:     succeeding term " << successor << " does not involve start_gauss_crossing" << endl;
						}
					}
					else
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present:     term " << j << " does not involve start_gauss_crossing or rendezvous_gauss_crossing" << endl;
					}
				}
			}
			
			/* If it looks like we have a Reidemeister II configuration, for standard Reidemeister II configrations in knotoids check 
			   that the shortcut does not lie within the bigon boundary.  If this is the case we cannot remove the bigon in the the 
			   underlying immersion and so we ignore this Reidemeister configuration.  Note that this is Gauss Reidemeister II detection, 
			   so the Reidemeister II crossings wil be classical, hence we cannot have just one of the leg or head
			   semi-arcs in the bigon boundary without the other.
			   
			   If the shortcut lies within the bigon bundary, then the leg of the knotoid, numbered 0, lies in the bigon boundary.
			   
			   We have a similar scenario if we are dealing with a knot-type knotoid or a long knot , and again, we may use the 
			   presence of edge zero in the bigon boundary to detect this case.
			   
			   Since we don't know which of end_bigon_edge_a or end_bigon_edge_b lies on the same path as the start_edge, we trace
			   forward from start_edge until we reach one of the end edges, checking for edge zero as we do so.  If we don't find
			   edge zero, we can then check the peer path by comparing the start_peer with the edge at the other end of
			   the peer path, since the numbering will have wrapped if edge zero lies on the peer path.
			*/
			if (Reidemeister_II_found && code_data.immersion_character != generic_code_data::character::CLOSED)
			{
				int rendezvous_end_edge;
				
                               /* note that the start_edge does not lie in the bigon boundary but may be edge zero, 
				  in which case the Reidemeister II move is possible
                               */   
				for (int j = 1; j < code_data.num_component_edges[first_component]; j++)
				{
					int edge = (start_edge + j - code_data.first_edge_on_component[first_component]) % code_data.num_component_edges[first_component] + code_data.first_edge_on_component[first_component];
					
					if (edge == 0)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present: knotoid check: edge zero found on the start path, cannot remove the bigon in the underlying immersion" << endl;
	
						Reidemeister_II_found = false;
						break;
					}
						
					if (edge == end_bigon_edge_a)
					{
						rendezvous_end_edge = end_bigon_edge_b;
						break;
					}
					else if (edge == end_bigon_edge_b)
					{
						rendezvous_end_edge = end_bigon_edge_a;
						break;
					}
				}

				if (Reidemeister_II_found)
				{

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present: knotoid check: start path does not contain edge zero, rendezvous_end_edge identified as " << rendezvous_end_edge << endl;
				
					if (second_component == 0)
					{
						if ((forwards_from_peer && start_peer > rendezvous_end_edge) || (!forwards_from_peer && rendezvous_end_edge > start_peer))
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present: knotoid check: edge zero found on the rendezvous path, cannot remove the bigon in the underlying immersion" << endl;
							Reidemeister_II_found = false;
						}
					}
					
					if (Reidemeister_II_found)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present: knotoid check: rendezvous path does not contain edge zero" << endl;
					}
				}
			}
			
			
			if (Reidemeister_II_found)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Gauss_Reidemeister_II_present: Gauss_Reidemeister_II found, with rendezvous starting at term " << rendezvous_start_term 
	      << ", moving " << (forwards_from_peer? "forwards": "backwards") << " from start_peer " << endl;
	debug << "Gauss_Reidemeister_II_present: starting immersion crossing = " << start_immersion_crossing << ", rendezvous immersion crossing = " 
	      << rendezvous_immersion_crossing  << endl;
}
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Gauss_Reidemeister_II_present: no corresponding pair of terms found, not a Gauss_Reidemeister_II" << endl;
			}			
		}
		
		if (Reidemeister_II_found)
			break;
	}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE && !Reidemeister_II_found)
	debug << "Gauss_Reidemeister_II_present: no Gauss Reidemeister II move found" << endl;
	
	return _return;
}


/* virtual_Reidemeister_II_plus looks for simple Reidemeister II configurations where both crossings are virtual or shortcut crossings, 
   or Reidemeister II-like detours involving either virtual or shortcut crossings within the detour.  If the boolean mark_special is 
   true, the function sets the labels of the start and end crossings of the virtual Reidemeister II move to be generic_code_data::VOID.  
   This is to assist the function remove_Reidemeister_II
   
   The function traces the diagram in the manner determined by the edge numbering and, for each start_edge that terminates at a virtual 
   or shortcut crossing, considers all other jth_edges in the same component that can be reached via a contiguous sequence of virtual or 
   shortcut crossings.   Each such jth_edge it finds terminates a potential start path to a rendezvous crossing, so the function serches 
   forwards around the diagram for the subsequent visits to the start_crossing and rendezvous_crossing. For each subsequent visit, the 
   function looks forwards again to see if the other crossing can be reached by another contiguous sequence of virtual or shortcut 
   crossings, i.e a via a valid peer path.  The function returns "true" as soon as it finds a Reidemeister II configuration that reduces
   the number of crossings when the Reidemeister II configuration is removed.
   
   In fact, the function detects some virtual detour scenaios that are not Reidemeister II configurations but which bear a strong 
   similiarity with Reidemeister II configurations and which can also reduce the number of crossings when removed.
   
   Considering the start and rendezvous crossings oriented left to right as shown below: 
      
      A _   __ C      E _   __ G 
         \ /             \ /  	
	      x               x 		
 	  B _/ \__ D      F _/ \__ H
   
   The function will detect the following combinations of start path and rendezvous paths:
   
       1. start-path B-C...E-H rendezvous paths A-D...F-G (forwards from peer) and F-G...A-D (backwards from peer)
       2. start-path A-D...F-G rendezvous paths B-C...E-H (forwards from peer) and E-H...B-C (backwards from peer)
       3. start-path B-C...F-G rendezvous paths A-D...E-H (forwards from peer) and E-H...A-D (backwards from peer)
       4. start-path A-D...E-H rendezvous paths B-C...F-G (forwards from peer) and F-G...B-C (backwards from peer)
 
   The scenarios in 1 and 2 are standard reidemeister II consifurations in the sense that the Gauss orientation has the form LiRj...RiLj or 
   RiLj...LiRj and in these cases both the start and rendezvous crossings may be removed if the start or peer path is moved to lie parallel 
   to the other in such a manner as to remove the start crossing.
   
   The scenarios in 3 and 4 have Gauss orientations of the form LiLj...RiRj or RiRj...LiLj.  In these cases, if we remove the Reidemeister
   configuration as in the function remove_Reidemeister_II, by moving the peer path so that it is parallel to the start path and removing 
   the virtual crossing at the initial (i.e the peer_detour_entry_edge end), then it is necessary to retain the virtual crosing at the other 
   (terminating) end of the peer path.  This is because we shall end up on the opposite side of the start path to the continuation of the peer 
   path component.
   
   Thus the scenarios in 3 and 4 are examples of non-Reidemeister II detours that may be detected by the function and may still result in a 
   simplification of a diagram if the number of virtual crossings is reduced by moving the start of peer path, hence the word "plus" in name of the function.
*/
virtual_Reidemeister_II_plus_return virtual_Reidemeister_II_plus(generic_code_data& code_data, matrix<int>& cycle, int num_cycles, int num_left_cycles, bool mark_special)
{
	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;
	int start_peer;  // peer of start_edge only used for debugging
	int start_bigon_edge_a;
	int start_bigon_edge_b;
	int end_bigon_edge_a;
	int end_bigon_edge_b;
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "virtual_Reidemeister_II_plus (full): presented with code data ";
	write_code_data(debug,code_data);
	debug << endl;
}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	print_code_data(debug,code_data,"virtual_Reidemeister_II_plus: ");
}	
	
	matrix<int>& code_table = code_data.code_table;
	vector<int>& first_edge_on_component = code_data.first_edge_on_component;
	vector<int>& num_component_edges = code_data.num_component_edges;

	vector<int>& shortcut_crossing = code_data.shortcut_crossing;

	virtual_Reidemeister_II_plus_return	_return;
	int& start_edge = _return.start_edge;
	bool& forwards_from_peer = _return.forwards_from_peer;
	bool& move_peer_path_allowed = _return.move_peer_path_allowed;
	bool& standard_Reidemeister_II_configuration = _return.standard_Reidemeister_II_configuration;
	
	bool ignore_shortcut = false;
	if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && shortcut_crossing.size())
		ignore_shortcut = true;
	
	for (int i=0; i< num_edges; i++)
	{
		start_edge = i;
		int start_crossing = code_data.term_crossing[i];
		int start_edge_parity;
		
		if ((start_edge%2 == 1 && code_table[generic_code_data::table::TYPE][start_crossing] == generic_code_data::TYPE1) ||
		    (start_edge%2 == 0 && code_table[generic_code_data::table::TYPE][start_crossing] == generic_code_data::TYPE2))
			start_edge_parity = gauss_orientation_data::RIGHT;
		else
			start_edge_parity = gauss_orientation_data::LEFT;
		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "virtual_Reidemeister_II_plus: edge " << i << " crossing = " << start_crossing << ", label = " 
	      << code_table[generic_code_data::table::LABEL][start_crossing] << ", start_edge_parity = " << (start_edge_parity == gauss_orientation_data::RIGHT? "RIGHT" : "LEFT") << endl;
}	

		/* we're only interested in starting at edges related to virtual or shortcut crossings */
		if (code_table[generic_code_data::table::LABEL][start_crossing] != generic_code_data::VIRTUAL && !(ignore_shortcut && shortcut_crossing[start_crossing]))
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   edge " << i << " does not terminate at a virtual crossing or a shortcut crossing" << endl;
	
			continue;
		}

		if (code_table[generic_code_data::table::LABEL][start_crossing] == generic_code_data::VIRTUAL)
		{		

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   edge " << i << " terminates at a virtual crossing, check all subsequent crossings reached via virtual or shortcut crossings" << endl;
		}
		else
		{	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   edge " << i << " terminates at a shortcut crossing, check all subsequent crossings reached via virtual or shortcut crossings" << endl;
		}

		/* we track the virtual crossings on the start and peer paths to assess whether any Reidemeister II 
		   configuration that we find warrants detecting.  Since we are successively searching from the ith edge, 
		   we build up start_virtual_crossing_flag as we proceed from the ith edge and then complete 
		   peer_virtual_crossing_flag when we look for subsequent paths.  We also look out for shortcut crossings
		   in the start path.
		*/
		vector<int> start_virtual_crossing_flag(num_crossings);
		start_virtual_crossing_flag[start_crossing] = 1;  
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   set start_virtual_crossing_flag[" << start_crossing << ']' << endl;

		start_edge = i;
		if (start_edge%2 == 1)
			start_peer = code_table[generic_code_data::table::EVEN_TERMINATING][start_crossing];
		else
			start_peer = code_table[generic_code_data::table::ODD_TERMINATING][start_crossing];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   start_edge = " << start_edge << ", start_peer = " << start_peer << endl;

		/* Look forward from the ith edge around the component to which start_edge belongs at each subsequent edge for which all 
		   intervening edges are virtual or shortcut crossings and see if there is a virtual Reidemeister II move between them.  
		   If there is, then the subsequent path between the same two crossings, around the component on which the subsequent path
		   resides, will also comprise only virtual or shortcut crossings.
		*/
		int ith_edge_component = code_table[generic_code_data::table::COMPONENT][(start_edge%2 == 0? start_edge/2: (start_edge-1)/2)];
		int num_initial_component_edges = code_data.num_component_edges[ith_edge_component];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   look around component " << ith_edge_component << ", which has " << num_initial_component_edges << " edges" << endl;

		
		/* look forwards from edge i around the initial component*/
		for (int j=1; j< num_initial_component_edges-1; j++)
		{
			int jth_edge = (i-code_data.first_edge_on_component[ith_edge_component]+j)%num_initial_component_edges + code_data.first_edge_on_component[ith_edge_component];
			
			int jth_crossing = code_data.term_crossing[jth_edge];		

			int jth_edge_parity;
			
			if ((jth_edge%2 == 1 && code_table[generic_code_data::table::TYPE][jth_crossing] == generic_code_data::TYPE1) ||
				(jth_edge%2 == 0 && code_table[generic_code_data::table::TYPE][jth_crossing] == generic_code_data::TYPE2))
				jth_edge_parity = gauss_orientation_data::RIGHT;
			else
				jth_edge_parity = gauss_orientation_data::LEFT;
				
			
			if (code_table[generic_code_data::table::LABEL][jth_crossing] != generic_code_data::VIRTUAL && !(ignore_shortcut && shortcut_crossing[jth_crossing]))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:     edge " << jth_edge << " does not terminate at a virtual or shortcut crossing, stop checking crossings reached from edge " << i << endl;
				break;
			}
			else if (jth_crossing == start_crossing)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:     edge " << jth_edge << " returns to start_crossing " << start_crossing << " stop checking crossings reached from edge " << i << endl;
				break;
			}
			else
			{
				/* the jth_crossing is a potential rendezvous crossing */
				int rendezvous_crossing = jth_crossing;			
				start_virtual_crossing_flag[rendezvous_crossing] = 1;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:     set start_virtual_crossing_flag[" << rendezvous_crossing << ']' << endl;			
				
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "virtual_Reidemeister_II_plus:     checking for a virtual Reidemeister II configuration between start_edge " << start_edge << " and jth_edge " << jth_edge << endl;				
	debug << "virtual_Reidemeister_II_plus:     start_crossing = " << start_crossing << ", rendezvous_crossing = " << rendezvous_crossing << endl;
	debug << "virtual_Reidemeister_II_plus:     jth_edge_parity = " << (jth_edge_parity == gauss_orientation_data::RIGHT? "RIGHT" : "LEFT") << endl;
}	
				if (start_edge_parity == jth_edge_parity)
					standard_Reidemeister_II_configuration = false;
				else
					standard_Reidemeister_II_configuration = true;
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:     standard_Reidemeister_II_configuration = " << standard_Reidemeister_II_configuration << endl;					
				
				int peer_path_start_edge;
				int subsequent_path_end_crossing;
				
				/* Look for the start of the subsequent path between the crossings at which start_edge and jth_edge terminate. 
				   For any Reidemeister II configuration, we will detect its existence when we reach the first of the two 
				   possible start_edge semi-arcs, when tracing the diagram according to the edge numbering.  From this first
				   possible start_edge the subsequent path might lie on a subsequent component but cannot lie on an earlier
				   component, otherwise we would not be at the first possible start-edge. Therefore, in order to look for a
				   subsequent path we look around the component containing the jth_edge until we get back to start_edge and
				   then look around any subsequent components.
				   
				   We have to check both of the subsequent encounters with these crossings, since the fist one we find may not 
				   be the one that reveals a virtual Reidemeister II.  An example of this is, 
				   
				   [5 -13^, 19 -15 21 23 -17, -3 -9 1, 7, 11]/+-**++*+*+++
				   			   
				   In this case, we identify a potential Reidemeister II starting at edge 12 (crossing 6) and edge 4 
				   (crossing 2).  However, searching forwards from edge 12 we first enounter crossing 6 again 
				   at edge 17, before we encounter crossing 2 at edge 19.  It is the latter occurrence cycling through
				   edges 14, 15, 16, and 17 that reveals the virtual Reidemeister II move.				   
				*/
				bool subsequent_start_found = false;
				int subsequent_encounters_found = 0;
				
				for (int component = ith_edge_component; component< code_data.num_components && subsequent_encounters_found < 2; component++)
				{
					/* we will iterate through the number of component edges but for the ith_edge_component we need to start after the jth_edge
					   and finish before the start_edge.  We will evaluate the component edge from the index t within the for loop, so just 
					   work out how many edges we need to consider first.
					*/
					
					int num_component_edges_to_consider;			
					if (component == ith_edge_component)
					{
					    /* if there a n edges in the ith component, the first edge is f, the start_edge is s and the jth edge is j, then the offsets
						   of s and j are (s-f) and (j-f) respectively.  However, the start of the component labelling might occur between
						   the start_edge and the jth_edge, so that (j-f) < (s-f).  That means we need to consider n - ((j-f)+n - (s-f))%n - 1 
						   edges, i.e. n - (j-f+n-s+f)%n -1 = n - (j+n-s)%n -1
						*/
						num_component_edges_to_consider = num_component_edges[component] - 
														  (jth_edge + num_component_edges[component] - start_edge)%num_component_edges[component] - 1;
					}
					else
						num_component_edges_to_consider = num_component_edges[component];

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "virtual_Reidemeister_II_plus:     component " << component << " num_component_edges_to_consider = " << num_component_edges_to_consider << endl;
	debug << "virtual_Reidemeister_II_plus:     first_edge_on_component " << first_edge_on_component[component] << " num_component_edges = " << num_component_edges[component] << endl;
}					

					for (int t=0; t< num_component_edges_to_consider && subsequent_encounters_found < 2; t++)									
					{	
						int component_edge;
						
						if (component == ith_edge_component)
							component_edge = first_edge_on_component[component] + (jth_edge+1-first_edge_on_component[component]+t)%num_component_edges[component];
						else
							component_edge = first_edge_on_component[component] + t;
						
						if (code_data.term_crossing[component_edge] == start_crossing)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "virtual_Reidemeister_II_plus:       edge " << component_edge << " is the second encounter with the start crossing " 
	      << start_crossing << endl;
}
							peer_path_start_edge = component_edge;
							subsequent_path_end_crossing = rendezvous_crossing;
							subsequent_start_found = true;
							subsequent_encounters_found++;
							forwards_from_peer = true;
							
							/* we need the start and end bigon edges to determine whether this Reidemeister II
							   configuration warrants detecting.
							*/						
							start_bigon_edge_a = code_table[generic_code_data::table::EVEN_ORIGINATING][start_crossing];
							start_bigon_edge_b = code_table[generic_code_data::table::ODD_ORIGINATING][start_crossing];
							end_bigon_edge_a = code_table[generic_code_data::table::EVEN_TERMINATING][rendezvous_crossing];
							end_bigon_edge_b = code_table[generic_code_data::table::ODD_TERMINATING][rendezvous_crossing];					
						}
						else if (code_data.term_crossing[component_edge] == rendezvous_crossing)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:       edge " << component_edge  << " is the second encounter with the rendezvous crossing " << rendezvous_crossing << endl;

							peer_path_start_edge = component_edge;
							subsequent_path_end_crossing = start_crossing;
							subsequent_start_found = true;
							subsequent_encounters_found++;
							forwards_from_peer = false;
	
							if (start_edge%2 == 1)
							{
								start_bigon_edge_a = code_table[generic_code_data::table::EVEN_ORIGINATING][start_crossing]; // start path
								start_bigon_edge_b = code_table[generic_code_data::table::EVEN_TERMINATING][start_crossing]; // peer path
							}
							else
							{
								start_bigon_edge_a = code_table[generic_code_data::table::ODD_ORIGINATING][start_crossing]; // start path
								start_bigon_edge_b = code_table[generic_code_data::table::ODD_TERMINATING][start_crossing]; // peer path
							}
	
							if (peer_path_start_edge%2 ==1)
							{
								end_bigon_edge_a = code_table[generic_code_data::table::EVEN_TERMINATING][rendezvous_crossing]; // start path
								end_bigon_edge_b = code_table[generic_code_data::table::EVEN_ORIGINATING][rendezvous_crossing]; // peer path
							}
							else
							{
								end_bigon_edge_a = code_table[generic_code_data::table::ODD_TERMINATING][rendezvous_crossing]; // start path
								end_bigon_edge_b = code_table[generic_code_data::table::ODD_ORIGINATING][rendezvous_crossing]; // peer path
							}						
						}
						else
						{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "virtual_Reidemeister_II_plus:       edge " << component_edge << " not in the subsequent path" << endl;
						}
	
						if (subsequent_start_found)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "virtual_Reidemeister_II_plus:         peer_path_start_edge = " << peer_path_start_edge << " subsequent_path_end_crossing = " 
	      << subsequent_path_end_crossing << ", forwards_from_peer = " << forwards_from_peer << endl;
}
							/* Look around the subsequent_component to see if we reach the subsequent_path_end_crossing via only
							   virtual or shrtcut crosings.  Note if we encounter any shortcut crossings as we do this.						   
							*/						
							vector<int> peer_virtual_crossing_flag(num_crossings);
							peer_virtual_crossing_flag[start_crossing] = 1;
							peer_virtual_crossing_flag[rendezvous_crossing] = 1;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "virtual_Reidemeister_II_plus:         set peer_virtual_crossing_flag[" << start_crossing << ']' 
	      << " and peer_virtual_crossing_flag[" << rendezvous_crossing << ']' << endl;
}						
							bool virtual_Reidemeister_II = true;
	
							int subsequent_component = code_table[generic_code_data::table::COMPONENT][(peer_path_start_edge%2 == 0? peer_path_start_edge/2: (peer_path_start_edge-1)/2)];
								
							int num_subsequent_component_edges = code_data.num_component_edges[subsequent_component];						
							bool peer_path_shortcut_crossing = false;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "virtual_Reidemeister_II_plus:         look around component " << subsequent_component << ", which has " 
	      << num_subsequent_component_edges << " edges" << endl;
}				
							for (int k=1; k<= num_subsequent_component_edges-1; k++)
							{
								int kth_edge = (peer_path_start_edge - code_data.first_edge_on_component[subsequent_component]+k)%num_subsequent_component_edges 
											   + code_data.first_edge_on_component[subsequent_component];							
								int kth_crossing = code_data.term_crossing[kth_edge];
								
								if ( kth_crossing == subsequent_path_end_crossing)
								{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:           edge " << kth_edge << " is the second encounter with crossing " << subsequent_path_end_crossing << endl;
	
									break;
								}
								else if (code_table[generic_code_data::table::LABEL][kth_crossing] != generic_code_data::VIRTUAL && !(ignore_shortcut && shortcut_crossing[kth_crossing]))
								{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:           edge " << kth_edge << " does not terminate at a virtual or shortcut crossing" << endl;
									virtual_Reidemeister_II = false;
									subsequent_start_found = false;
									break;
								}
								else
								{
									peer_virtual_crossing_flag[kth_crossing] = 1;	
		
									if ((kth_edge%2 == 0 && code_table[generic_code_data::table::LABEL][kth_crossing] == generic_code_data::NEGATIVE) ||
									     (kth_edge%2 == 1 && code_table[generic_code_data::table::LABEL][kth_crossing] == generic_code_data::POSITIVE))
									{
										peer_path_shortcut_crossing = true;
									}
										
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "virtual_Reidemeister_II_plus:           edge " << kth_edge << " terminates at a virtual or shortcut crossing, set peer_virtual_crossing_flag[" << kth_crossing << ']' << endl;
	debug << "virtual_Reidemeister_II_plus:           peer_path_shortcut_crossing = " << peer_path_shortcut_crossing << endl;
}
								}
							}
	
							if (virtual_Reidemeister_II)
							{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "virtual_Reidemeister_II_plus: found virtual Reidemeister II move between edges " << i << " and " << jth_edge  
	      << ", crossings " << start_crossing << " and " << rendezvous_crossing << endl;
}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "virtual_Reidemeister_II_plus: start_bigon_edge_a = " << start_bigon_edge_a << " start_bigon_edge_b = " << start_bigon_edge_b
		  << " end_bigon_edge_a = " << end_bigon_edge_a << " end_bigon_edge_b = " << end_bigon_edge_b << endl;
	debug << "virtual_Reidemeister_II_plus: start_virtual_crossing_flag ";
	for (int k=0; k< num_crossings; k++)
		debug << start_virtual_crossing_flag[k] << ' ';
	debug << endl;
	debug << "virtual_Reidemeister_II_plus: peer_virtual_crossing_flag  ";
	for (int k=0; k< num_crossings; k++)
		debug << peer_virtual_crossing_flag[k] << ' ';
	debug << endl;
}
								/* we have to decide whether this Reidemeister II configuration warrants detecting.  That depends
								   upon the nature of the start and end crossings (virtual or shortcut) and whether removing the 
								   bigon would reduce the number of immersion crossings or not.  
	
								   If the start and end crossings lie in the same turning cycle, we can simply remove all
								   the crossings from the start and peer path, since at least one of the paths is a virtual 
								   detour not contianing part of the shortcut and we will always reduce the number of crossings
								   by removing the bigon in this case.
								   
								   If the start and end crossings are not in the same turning cycle, we want to move one 
								   of the paths to be parallel with the other to remove the Reidemeister II bigon and this 
								   process may add more virtual crossings than it removes, so we have to check whether that 
								   is the case.  								   
								*/
								bool same_turning_cycle = false;
								for (int i=0; i< num_cycles; i++)
								{
									int found_count = 0;
									
									for (int j=1; j<= cycle[i][0]; j++)
									{	
										/* an edge might be a start edge and a rendezvous edge, if it is, we want to count it twice */				
										if (abs(cycle[i][j]) == start_bigon_edge_a)
											found_count++;
										else if (abs(cycle[i][j]) == start_bigon_edge_b)
											found_count++;
										
										if (abs(cycle[i][j]) == end_bigon_edge_a)
											found_count++;
										else if (abs(cycle[i][j]) == end_bigon_edge_b)
											found_count++;
									}
									
									if (found_count == 4)
									{
										same_turning_cycle = true;
										
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus: start and rendezvous crossings lie in the same turning cycle, keep Reidemeister II configuration" << endl;

										break;
									}
								}
								
								if (!same_turning_cycle)
								{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "virtual_Reidemeister_II_plus: start and rendezvous crossings not in the same turning cycle, check crossing counts" << endl;
	debug << "virtual_Reidemeister_II_plus:   initial num_crossings = " << num_crossings << endl;
}
									bool move_start_path_allowed = true;
									move_peer_path_allowed = true;
																																
									/* We need to count virtual crossings, since it's not necessarily the case that we
									   want to consider this Reidemeister II configuration.  Note that the start and peer
									   path may intersect.
									*/							
									int num_start_virtual_crossings = 0;
									for (int i=0; i< num_crossings; i++)
									{
										if (start_virtual_crossing_flag[i]  && !peer_virtual_crossing_flag[i]) // not on both paths
											num_start_virtual_crossings++;
									}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   number of non-peer-path start path virtual crossings = " << num_start_virtual_crossings << endl;

									/* we will not have counted the start and rendezvous crossings in the above calculation */
									int new_num_crossings = num_crossings+num_start_virtual_crossings;
									for (int i=0; i< num_crossings; i++)
									{
										if (peer_virtual_crossing_flag[i])
											new_num_crossings--;
									}
									
									if (!standard_Reidemeister_II_configuration)
										new_num_crossings++;  // we can't remove both the start and rendezvous virtual crossings in these cases
									
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   moving the peer path results in new_num_crossings = " << new_num_crossings << endl;

									if (new_num_crossings >= num_crossings)
									{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   moving the peer path does not reduce the number of crossings"  << endl;

										move_peer_path_allowed = false;
									}																		
									else
									{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   moving the peer path reduces the number of crossings"  << endl;
									}																											
																		
									int num_peer_virtual_crossings = 0;
									for (int i=0; i< num_crossings; i++)
									{
										if (peer_virtual_crossing_flag[i]  && !start_virtual_crossing_flag[i]) // not on both paths
											num_peer_virtual_crossings++;
									}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   number of non-start-path peer path virtual crossings = " << num_peer_virtual_crossings << endl;
									
									/* we will not have counted the start and rendezvous crossings in the above calculation */
									new_num_crossings = num_crossings+num_peer_virtual_crossings;
									for (int i=0; i< num_crossings; i++)
									{
										if (start_virtual_crossing_flag[i])
											new_num_crossings--;
									}

									if (!standard_Reidemeister_II_configuration)
										new_num_crossings++;  // we can't remove both the start and rendezvous virtual crossings in these cases
									
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   moving the start path results in new_num_crossings = " << new_num_crossings << endl;
										
									if (new_num_crossings >= num_crossings)
									{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   moving the start path does not reduce the number of crossings"  << endl;

										move_start_path_allowed = false;
									}
									else
									{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:   moving the start path reduces the number of crossings"  << endl;
									}									

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "virtual_Reidemeister_II_plus:   move_peer_path_allowed = " << move_peer_path_allowed << endl;
	debug << "virtual_Reidemeister_II_plus:   move_start_path_allowed = " << move_start_path_allowed << endl;
}				
																	
									if (!move_peer_path_allowed && !move_start_path_allowed)
									{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus: removing the Reidemeister II configuration does not reduce the number of crossings, discarded" << endl;
	
											virtual_Reidemeister_II = false;
									}
									else
									{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus: removing the Reidemeister II configuration reduces the number of crossings, retained" << endl;
									}
								}
								   	
								if (virtual_Reidemeister_II)
								{
									if (mark_special)
									{
										code_data.code_table[generic_code_data::table::LABEL][start_crossing] = generic_code_data::VOID;
										code_data.code_table[generic_code_data::table::LABEL][rendezvous_crossing] = generic_code_data::VOID;
									}
									
									_return.Reidemeister_II_found = true;
									return _return;
								}						
							}
							else
							{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:         no virtual Reidemeister II peer path from peer_path_start_edge = " << peer_path_start_edge << endl;
							}
						}
					}
				}	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "virtual_Reidemeister_II_plus:     no virtual Reidemeister II move between start_edge " << start_edge << " and jth_edge " << jth_edge << endl;				
			}			
		}
	}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "virtual_Reidemeister_II_plus: no virtual Reidemeister II move found" << endl;
	
	_return.Reidemeister_II_found = false;
	return _return;
}


/* remove_Reidemeister_II recursively looks for simple Reidemeister II moves and then Reidemeister II detours, adjusting the code_data
   accordingly at each step.  
   
   By looking for simple moves first, we ensure that, if we ever detect a Reidemeister II detour, the two paths 
   between the starting and rendezvous crossing that comprise the Reidemeister II detour may be categorized as combinations of
   the following configurations (referred to as a to d, left to right).  
   

                |  |                                                                                      |
              +-----+                   |  |  |                    |                         ----------   |     ---------- 
			--o ... o--              ---o--#--o---              ---#---     ------           |        |   |     |        |        
	   \   /  +-----+  \   /    \   /   |  |  |   \   /    \   /   |   \   /      \   /   ---o--------+---#-----o--------+---    
		\ /             \ /      \ /  +-+--#--+-+  \ /      \ /    |    \ /        \ /       |        |   |     |        |   
		 x               x        x   |    |    |   x        x     |     o          x        |        |   |     |        |  
		/ \             / \      / \  +-+--#--+-+  / \      / \    |    / \        / \       |        |   |     |      
	   /   \  +-----+  /   \    /   \   |  |  |   /   \    /   \   |   /   \      /   \      ---------o---#------   
			--o ... o--              ---o--#--o---              ---#---     ------                    |   |
		      +-----+                   |  |  |                    |                                  |   |
                |  |

   Actually, d is a special case of c, where the two detour paths intersect and one of them also self intersects.
   
   If the Reidemeister II bigon intersects a knotoid shortcut, then the shortcut will either loop in and out, as in case a, 
   or will pass from top to bottom in cases b c or d as shown using # to indicate an under-crossing.  
   This is because the shortcut does not self-intersect.
   
   
   Case a is most simply handled, since the start and rendezvous crossings of the Reidemester II move lie in the same turning cycle.
   Here the boxes may involve both virtual and non-virtual crossings, though the detour paths encounter only virtuals.  Any transverse edges
   enter and leave the boxes on the outside of the Reidemeister II bigon.  In this case, removing the Reidemeister II removes all the virtual 
   crossings on both detour paths as well as the two non-virtual Reidemeister II crossings.
   
   In case b (c and d) both detour paths include virtual or shortcut crossings and we move one of the detour paths 
   (the peer detour path involving the peer of the (start) from which the Reidemeister II
   configuration is first identified) so that it is parallel with the other, as follows for case b (cases c and d are similar)
                
                                           |   |          
	                            \   -      |   |       -   /
		                         \ / \   +-+---+-+    / \ /   
		                          x   \  |       |   /   x    
		                         / \   \ +-+---+-+  /   / \   
                                /   \   \  |   |   /   /   \
	                           /     \   --o---o---   /     \  
			                          -----o---o------       
		                                   |   |                                                      |
   
   From here, the transverse virtual crossings may be moved across one of the non-virtual crossings and the Reidemeister II bigon removed.  
   The result is that whilst we lose those virtual crossings from the peer detour, those on the start detour are retained and new virtual 
   crossings introduced alongside them, as the peer detour is moved.
   
   Cases c and d demonstrate that we need to consider instersections between the two detours in a variety of ways.  Note that case a could be handled 
   in the same was as case b but it is quicker to handle it separately.  Note also that the above procedure may introduce virtual Reidemeister I 
   monogons but these are removed by subsequent recursions of the function.
         
   The function returns the number of components removed from the peer code as a result of disconnecting it, or removing simple closed closed curve 
   components that are connected to the other components only by crossings that are being removed.
*/
int remove_Reidemeister_II(generic_code_data& code_data, vector<int>& component_flags)
{
	int num_components_removed = 0;  //the return variable
	int num_crossings = code_data.num_crossings;
	int num_components = code_data.num_components;
	int floating_component;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: presented with code data ";
	write_code_data(debug,code_data);
	debug  << endl;
	debug << "remove_Reidemeister_II: component_flags:";
	for (unsigned int j=0; j< component_flags.size(); j++)
		debug << component_flags[j] << ' ';
	debug << endl;
	debug << "remove_Reidemeister_II:  peer code contains " << num_crossings << " crossings" << endl;
}

	gauss_orientation_data gauss_data(code_data);
	
	matrix<int>& code_table = code_data.code_table;
	matrix<int> copy_code_table = code_table;  // used to reset VOID crossings
	vector<int>& first_edge_on_component = code_data.first_edge_on_component;
	vector<int>& num_component_edges = code_data.num_component_edges;
	vector<int>& term_crossing = code_data.term_crossing;
	bool move_peer_path_allowed = true;
	bool standard_Reidemeister_II_configuration = false; 

	vector<int>& shortcut_crossing = code_data.shortcut_crossing;

	bool pure_knotoid_code_data = false;
	int head_semi_arc = -1;
	
	if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && shortcut_crossing.size())
	{
		if (code_table[generic_code_data::table::LABEL][code_data.head] == generic_code_data::POSITIVE)
			head_semi_arc = code_table[generic_code_data::table::OPEER][code_data.head];
		else if (code_table[generic_code_data::table::LABEL][code_data.head] == generic_code_data::NEGATIVE)
			head_semi_arc = 2*code_data.head;
		else
		{
			cout << "\nError! Function remove_Reidemeister_II presented with a knotoid code whose first shortcut crossing is virtual." << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: head_semi_arc = " << head_semi_arc << endl;		
	
		pure_knotoid_code_data = true;
	}

	bool track_zig_zag_counts;
	if (code_data.zig_zag_count.numrows() !=0)
		track_zig_zag_counts = true;
	else
		track_zig_zag_counts = false;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: track_zig_zag_counts = " << track_zig_zag_counts << endl;			
		
	int num_left_cycles=0;
	int num_cycles=0;
	
	matrix<int> cycle(num_crossings+2, 2*num_crossings+1);
	calculate_turning_cycles(code_data, cycle, num_left_cycles, num_cycles);
	
	/* We record which edges terminate at a crossing involved in a Reidemeister II move */
	simple_Reidemeister_II_return simple_Reidemeister_II_data = simple_Reidemeister_II_present(code_data, true); // create_flags = true
	vector<int>& RII_edge_flag = simple_Reidemeister_II_data.RII_edge_flag; 

	bool Reidemeister_II_found = simple_Reidemeister_II_data.Reidemeister_II_found;

	if (Reidemeister_II_found)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: simple Reidemister II move found" << endl;
	
		/* First check that removing the simple Reidemeister II move does not reduce the diagram to the unknot */
		int num_edges_remaining = 2*num_crossings;
		for (int i=0; i< 2*num_crossings; i++)
			num_edges_remaining -= RII_edge_flag[i];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:  num_edges_remaining = "<< num_edges_remaining << endl;
			
		if (num_edges_remaining == 0)
		{	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:  removing simple Reidemeister II moves reduces the diagram to as set of disjoint simple closed curves" << endl;
	
			/* we clear the component flags below in this case */
			code_data.num_crossings = 0;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:  removing RII_edge_flags ";
	for (int i=0; i< 2*num_crossings; i++)
		debug << RII_edge_flag[i] << ' ';
	debug << endl;
}
			num_components_removed += remove_edge_flags_from_peer_code(code_data,RII_edge_flag,component_flags);
		}
	}
	else
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: no simple Reidemeister II move found, looking for a Reidemeister II detour" << endl;
		/* look for a Reidemeister II detour */
		Gauss_Reidemeister_II_return Gauss_Reidemeister_II_data = Gauss_Reidemeister_II_present(code_data,gauss_data);

		int start = Gauss_Reidemeister_II_data.start_edge;
		bool forwards_from_peer = Gauss_Reidemeister_II_data.forwards_from_peer;
		floating_component = Gauss_Reidemeister_II_data.floating_component;

		Reidemeister_II_found = Gauss_Reidemeister_II_data.Reidemeister_II_found;

		if (floating_component != 0)  // then we have detected that component floating_component-1 contains only virtual crossings
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: floating component " << floating_component-1 << " detected by Gauss_Reideister_II_present" << endl;
	debug << "remove_Reidemeister_II: Reidemeister_II_found is " << Reidemeister_II_found << " (should be false, since here we have a floating component)" << endl;
}
			if (false && pure_knotoid_code_data && floating_component == 1)
			{
				pure_knotoid_code_data = false;
				code_data.head = -1;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: knotoid: floating segment component clearing pure_knotoid_code_data and setting code_data.head = -1" << endl;
				
			}
			
			num_components_removed = remove_peer_code_component(code_data, floating_component-1, component_flags);
		}
		else if (Reidemeister_II_found)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: (standard) Gauss Reidemister II detour found" << endl;
			standard_Reidemeister_II_configuration = true;
		}
		else
		{			
			virtual_Reidemeister_II_plus_return virtual_Reidemeister_II_plus_data = virtual_Reidemeister_II_plus(code_data, cycle, num_cycles, num_left_cycles, true); // mark_special=true
			                                                     
			Reidemeister_II_found = virtual_Reidemeister_II_plus_data.Reidemeister_II_found;
			start = virtual_Reidemeister_II_plus_data.start_edge;
			forwards_from_peer = virtual_Reidemeister_II_plus_data.forwards_from_peer;
			move_peer_path_allowed = virtual_Reidemeister_II_plus_data.move_peer_path_allowed;
			standard_Reidemeister_II_configuration = virtual_Reidemeister_II_plus_data.standard_Reidemeister_II_configuration;

			if (Reidemeister_II_found)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: virtual Reidemister II detour found" << endl;
	print_code_data(debug, code_data, "remove_Reidemeister_II: ");
}

	
				if (!move_peer_path_allowed) 
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: virtual Reidemister II detour requires start path to be moved, not the peer path" << endl;					
				}
			}			
		}
		
		if (Reidemeister_II_found)
		{
			/* identify the virtual crossings and edge flags associated with the detour */
			int start_crossing = term_crossing[start];	
			int start_label = code_table[generic_code_data::table::LABEL][start_crossing];

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: Reidemister II detour found" << endl;
	debug << "remove_Reidemeister_II: start = " << start << ", start_crossing " << start_crossing << ", label = ";
	switch (start_label)
	{
		case generic_code_data::POSITIVE: debug << "POSITIVE"; break;
		case generic_code_data::NEGATIVE: debug << "NEGATIVE"; break;
		case generic_code_data::VIRTUAL: debug << "VIRTUAL"; break;
		case generic_code_data::FLAT: debug << "FLAT"; break;
		case generic_code_data::VOID: debug << "VOID"; break;
	}
	debug << endl;
}
			int peer;
			int start_component;
			int peer_component;
			
			if (start%2) // arriving on odd edge
			{
				peer = code_table[generic_code_data::table::EPEER][(start-1)/2];
				start_component = code_table[generic_code_data::table::COMPONENT][(start-1)/2];
				peer_component = code_table[generic_code_data::table::COMPONENT][peer/2]; 
			}
			else
			{
				peer = code_table[generic_code_data::table::OPEER][start/2];
				start_component = code_table[generic_code_data::table::COMPONENT][start/2];
				peer_component = code_table[generic_code_data::table::COMPONENT][(peer-1)/2]; 
			}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: peer = " << peer << endl;	
	debug << "remove_Reidemeister_II: start component = " << start_component << ", peer_component = " << peer_component << endl;	
}
			vector<int> start_virtual_crossing_flag(num_crossings);
			vector<int> peer_virtual_crossing_flag(num_crossings);
			int edge;
			int rendezvous_crossing;
			
			/* determine the virtual crossings moving forwards from start.
			
			   If we're dealing with a knotoid, note whether the head_semi_arc or leg_semi_arc lies within start path.
			*/			
			bool head_semi_arc_in_start_path = false;
			bool leg_semi_arc_in_start_path = false;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   moving forwards from start..." << endl;
			for (int i=1; i< num_component_edges[start_component]; i++)
			{
				edge = (start+i - first_edge_on_component[start_component])%num_component_edges[start_component] + first_edge_on_component[start_component];
				
				if (edge == 0)
					leg_semi_arc_in_start_path = true;
				else if (edge == head_semi_arc)
					head_semi_arc_in_start_path = true;
					
				if (code_table[generic_code_data::table::LABEL][term_crossing[edge]] == generic_code_data::VOID)
				{
					rendezvous_crossing = term_crossing[edge];
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   edge " << edge << " terminates at the rendezvous crossing, " << rendezvous_crossing << endl;
	
					break;
				}
				else if (code_table[generic_code_data::table::LABEL][term_crossing[edge]] == generic_code_data::VIRTUAL)
				{
					start_virtual_crossing_flag[term_crossing[edge]] = 1;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   edge " << edge << " terminates at a virtual crossing" << endl;
				}
				else if (pure_knotoid_code_data && shortcut_crossing[term_crossing[edge]])
				{
					start_virtual_crossing_flag[term_crossing[edge]] = 1;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   edge " << edge << " terminates at a shortcut crossing" << endl;
				}
				else
				{
					rendezvous_crossing = term_crossing[edge];
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   edge " << edge << " terminates at the rendezvous crossing, " << rendezvous_crossing << endl;
	
					break;
				}
			}

			/* determine the virtual crossings moving forwards or backwards from peer.  Also record the peer_detour_entry_edge,
			   this is the terminating edge on the peer detour at the start or rendezvous crossing when the moving forward or
			   backwards from the start respectively.  It is used to compare with a transverse skip edge later in order to determine
			   which occurrence of a skip edge we encounter first.
			   
			   If we're dealing with a knotoid, note whether the head_semi_arc or leg_semi_arc lies within peer path.
			*/
			int peer_detour_entry_edge = peer;
			bool head_semi_arc_in_peer_path = false;
			bool leg_semi_arc_in_peer_path = false;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	if (forwards_from_peer)
		debug << "remove_Reidemeister_II:   moving forwards from peer..." << endl;
	else
		debug << "remove_Reidemeister_II:   moving backwards from peer..." << endl;
}
						
			for (int i=1; i< num_component_edges[peer_component]; i++)
			{
				if (forwards_from_peer)
					edge = (peer+i - first_edge_on_component[peer_component])%num_component_edges[peer_component] + first_edge_on_component[peer_component];
				else
					edge = (peer-i - first_edge_on_component[peer_component]+num_component_edges[peer_component])%num_component_edges[peer_component] + first_edge_on_component[peer_component];

				if (edge == 0)
					leg_semi_arc_in_peer_path = true;
				else if (edge == head_semi_arc)
					head_semi_arc_in_peer_path = true;

				if (code_table[generic_code_data::table::LABEL][term_crossing[edge]] == generic_code_data::VOID)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   edge " << edge << " terminates at the rendezvous crossing, " << rendezvous_crossing << endl;
					if (!forwards_from_peer)
						peer_detour_entry_edge = edge;
					break;
				}
				if (code_table[generic_code_data::table::LABEL][term_crossing[edge]] == generic_code_data::VIRTUAL)
				{
					peer_virtual_crossing_flag[term_crossing[edge]] = 1;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   edge " << edge << " terminates at a virtual crossing" << endl;	
				}
				else if (pure_knotoid_code_data && shortcut_crossing[term_crossing[edge]])
				{
					peer_virtual_crossing_flag[term_crossing[edge]] = 1;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   edge " << edge << " terminates at a shortcut crossing" << endl;	
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   edge " << edge << " terminates at the rendezvous crossing, " << rendezvous_crossing << endl;
					if (!forwards_from_peer)
						peer_detour_entry_edge = edge;
					break;
				}
			}
			
			/* Now we know the start and rendezvous crossings we can reset any VOID labels */
			if (code_table[generic_code_data::table::LABEL][start_crossing] == generic_code_data::VOID)
			{
				code_table[generic_code_data::table::LABEL][start_crossing] = copy_code_table[generic_code_data::table::LABEL][start_crossing];
				code_table[generic_code_data::table::LABEL][rendezvous_crossing] = copy_code_table[generic_code_data::table::LABEL][rendezvous_crossing];
			}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   peer_detour_entry_edge = " << peer_detour_entry_edge << endl;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: start_virtual_crossing_flag ";
	for (int k=0; k< num_crossings; k++)
		debug << start_virtual_crossing_flag[k] << ' ';
	debug << endl;
	debug << "remove_Reidemeister_II: peer_virtual_crossing_flag  ";
	for (int k=0; k< num_crossings; k++)
		debug << peer_virtual_crossing_flag[k] << ' ';
	debug << endl;
	
	if (pure_knotoid_code_data)
	{
		debug << "remove_Reidemeister_II: head_semi_arc_in_start_path = " << head_semi_arc_in_start_path << endl;
		debug << "remove_Reidemeister_II: leg_semi_arc_in_start_path = " << leg_semi_arc_in_start_path << endl;
		debug << "remove_Reidemeister_II: head_semi_arc_in_peer_path = " << head_semi_arc_in_peer_path << endl;
		debug << "remove_Reidemeister_II: leg_semi_arc_in_peer_path = " << leg_semi_arc_in_peer_path << endl;
	}
}

			/* Determine whether the start and rendevous crossings lie in the same turning cycle.  This will tell
			   us whether we need to preserve the virtual crossings on the start path and move those on the peer path,
			   or whether we can simply remove them all.
			   
			   We know the successor to start and we have edge indicating the terminating edge at the rendezvous
			   crossing on the peer path.  We can therefore determine the start_path_rendezvous_edge, i.e. the other terminating 
			   edge at rendezvous crossing (which must be the one on the start path).  If these two start-path edges lie in 
			   the the same turning cycle as either the peer successor and edge, or peer and edge's successor (depending on 
			   whether we've traced forwards or backwards, then we're in the same turning cycle.			   
			*/
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: checking turning cycles" << endl;
			
			bool same_turning_cycle = false;
			
			int start_successor = (start+1 - first_edge_on_component[start_component])%
												num_component_edges[start_component] + first_edge_on_component[start_component];
			int start_path_rendezvous_edge = (edge %2 == 1? code_table[generic_code_data::table::EPEER][(edge-1)/2]: code_table[generic_code_data::table::OPEER][rendezvous_crossing]);
		
			int start_path_start_edge;
			int peer_path_start_edge;
			int peer_path_rendezvous_edge;
			
			if (forwards_from_peer)
			{
				start_path_start_edge = start_successor;
				peer_path_start_edge = (peer+1 - first_edge_on_component[peer_component])%
												num_component_edges[peer_component] + first_edge_on_component[peer_component];
				peer_path_rendezvous_edge = edge;
			}
			else
			{
				start_path_start_edge = start_successor;
				peer_path_start_edge = peer;
				peer_path_rendezvous_edge = (edge+1 - first_edge_on_component[peer_component])%
												num_component_edges[peer_component] + first_edge_on_component[peer_component];
			}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   starting edges " << start_path_start_edge << " and " << peer_path_start_edge << endl;
	debug << "remove_Reidemeister_II:   rendezvous edges " << start_path_rendezvous_edge << " and " << peer_path_rendezvous_edge << endl;
	debug << "remove_Reidemeister_II:   start_path_rendezvous_edge =" << start_path_rendezvous_edge << endl;
}
				
			for (int i=0; i< num_cycles; i++)
			{
				int found_count = 0;
				
				for (int j=1; j<= cycle[i][0]; j++)
				{	
					/* an edge might be a start edge and a rendezvous edge, if it is, we want to count it twice */				
					if (abs(cycle[i][j]) == start_path_start_edge)
						found_count++;
						
					if (abs(cycle[i][j]) == peer_path_start_edge)
						found_count++;
					
					if (abs(cycle[i][j]) == start_path_rendezvous_edge)
						found_count++;
						
					if (abs(cycle[i][j]) == peer_path_rendezvous_edge)
						found_count++;
				}
				
				if (found_count == 4)
				{
					same_turning_cycle = true;
					break;
				}
			}

			if (same_turning_cycle)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: start and rendezvous edges lie in the same turning cycle" << endl;

				/* we just need to remove all the virtual crossings we've identified plus the two non-virtual crossings */
				
				RII_edge_flag[code_table[generic_code_data::table::ODD_TERMINATING][start_crossing]] = 1;
				RII_edge_flag[code_table[generic_code_data::table::EVEN_TERMINATING][start_crossing]] = 1;
				RII_edge_flag[code_table[generic_code_data::table::ODD_TERMINATING][rendezvous_crossing]] = 1;
				RII_edge_flag[code_table[generic_code_data::table::EVEN_TERMINATING][rendezvous_crossing]] = 1;
				
				for (int i=0; i< num_crossings; i++)
				{
					if (start_virtual_crossing_flag[i] || peer_virtual_crossing_flag[i])
					{
						RII_edge_flag[code_table[generic_code_data::table::ODD_TERMINATING][i]] = 1;
						RII_edge_flag[code_table[generic_code_data::table::EVEN_TERMINATING][i]] = 1;
					}
				}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: removing RII_edge_flags ";
	for (int i=0; i< 2*num_crossings; i++)
		debug << RII_edge_flag[i] << ' ';
	debug << endl;
}
				/* check that removing the Reidemeister II detour does not reduce the diagram to the unknot */
				int num_edges_remaining = 2*num_crossings;
				for (int i=0; i< 2*num_crossings; i++)
					num_edges_remaining -= RII_edge_flag[i];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: num_edges_remaining = "<< num_edges_remaining << endl;

			
				if (num_edges_remaining == 0)
				{	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:  removing Reidemeister II move in the same turning cycle reduces the diagram to as set of disjoint simple closed curves" << endl;
	
					code_data.num_crossings = 0;
				}
				else
				{
					num_components_removed += remove_edge_flags_from_peer_code(code_data,RII_edge_flag,component_flags);
				}
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: start and rendezvous edges lie in different turning cycles" << endl;

				/* If we need to move the start path rather than the peer path, swap start_edge and peer_entry_edge,
				   start_component and peer_component, start_path_start_edge and peer_path_start_edge, rendezvous_edge_1 and peer_path_rendezvous_edge,
				   start_virtual_crossing_flag and peer_virtual_crossing_flag
				*/
				if (!move_peer_path_allowed)
				{					
					swap(start, peer_detour_entry_edge);
					swap(start_component, peer_component);
					swap(start_virtual_crossing_flag,peer_virtual_crossing_flag);
					swap(head_semi_arc_in_start_path,head_semi_arc_in_peer_path);
					swap(leg_semi_arc_in_start_path,leg_semi_arc_in_peer_path);
					
					if (forwards_from_peer)
					{
						swap(start_path_start_edge, peer_path_start_edge);
						swap(start_path_rendezvous_edge, peer_path_rendezvous_edge);
						peer = peer_detour_entry_edge; // was the start
					}
					else
					{
						swap(start_path_start_edge, peer_path_rendezvous_edge);
						swap(start_path_rendezvous_edge, peer_path_start_edge);
						peer = peer_path_start_edge; // was the start_path_rendezvous_edge
					}	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: moving the start path not the peer path" << endl;
	debug << "remove_Reidemeister_II:   start reset to " << start << ", peer_detour_entry_edge reset to " << peer_detour_entry_edge << endl;
	debug << "remove_Reidemeister_II:   peer reset to " << peer << endl;
	debug << "remove_Reidemeister_II:   start_component reset to " << start_component << ", peer_component reset to " << peer_component << endl;
	debug << "remove_Reidemeister_II:   start_path_start_edge reset to " << start_path_start_edge << ", peer_path_start_edge reset to " << peer_path_start_edge << endl;
	debug << "remove_Reidemeister_II:   start_path_rendezvous_edge " << start_path_rendezvous_edge << ", peer_path_rendezvous_edge reset to " << peer_path_rendezvous_edge << endl;
	debug << "remove_Reidemeister_II:   start_virtual_crossing_flag now ";
	for (int k=0; k< num_crossings; k++)
		debug << start_virtual_crossing_flag[k] << ' ';
	debug << endl;
	debug << "remove_Reidemeister_II:   peer_virtual_crossing_flag now  ";
	for (int k=0; k< num_crossings; k++)
		debug << peer_virtual_crossing_flag[k] << ' ';
	debug << endl;
	
	if (pure_knotoid_code_data)
	{
		debug << "remove_Reidemeister_II:   head_semi_arc_in_start_path reset to " << head_semi_arc_in_start_path << endl;
		debug << "remove_Reidemeister_II:   head_semi_arc_in_peer_path reset to " << head_semi_arc_in_peer_path << endl;
		debug << "remove_Reidemeister_II:   leg_semi_arc_in_start_path reset to " << leg_semi_arc_in_start_path << endl;
		debug << "remove_Reidemeister_II:   leg_semi_arc_in_peer_path reset to " << leg_semi_arc_in_peer_path << endl;
	}
	
}
					if (forwards_from_peer)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   forwards from peer: start crossing remains" << start_crossing 
	      << ", rendezvous_crossing remains " << rendezvous_crossing << endl;
}
					}
					else
					{
						swap(start_crossing,rendezvous_crossing);
						
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   backwards from peer: start crossing reset to " << start_crossing 
	      << ", rendezvous_crossing reset to " << rendezvous_crossing << endl;
}
					}

				}
				
				/* We are going to move the peer detour so that it lies alongside the start detour, then move those virtual
				   crossings (or shortcut crossings) across the Reidemeister II move, so we can remove it.  
				   
				   In the case of a knotoid, if we have the leg or head_semi_arc as a transverse edge incident with the
				   start detour, then the additional crossing added when moving the peer detour is a virtual or shortcut
				   crossing to match that on the start path.  See the diagram below.
				   
				   We count the number of additional virtual crossing when the peer detour is moved.  In some cases, e.g.
				   [-3 9 -15 -13 -11 1 -5 -7]/# # * # # * * *, the peer detour crosses the start detour, so that there is a
				   virtual crossing on both detours.  These crossings do not appear when we move the peer detour, since it
				   is moved to become parallel with the start detour.  However, we need to record how many there are initially,
				   since we have to set the generic_code_data::table::TYPE of the corresponding additional virtual crossings on the peer detour when it has
				   been moved.  Therefore we first count the initial number of start_virtual_crossings, then adjust them to 
				   remove those also lying on the peer detour.  Once this is done, the number of additional crossings added 
				   when the peer detour is moved is just the number of virtual crossings remaining on the start detour.  
				   Finally, we adjust for the two Reidemeister II crossings and the number of virtual crossings on the peer 
				   detour to calculate the new number of crossings.
				   
				   Note that when we are considering a non-standard (virtual) Reidemeister II configuration, we do not remove
				   the crossing at the terminating end of the peer path.  This is because the LiLj...RiRj or vice versa nature
				   of this configuration means that the moved peer path would be on the opposite side of the start path to the 
				   continuation of the peer path component if we tried to remove both crossings, so we have to retain one of them.  
				   Since we always wish to use the start crossing polarity to determine whether the peer path stays left or right 
				   of the start path, when the peer path moves backwards from the start, it is the start virtual crossing that has 
				   to be retained in this case.
				*/

				int num_initial_start_virtual_crossings = 0;
				for (int i=0; i< num_crossings; i++)
				{
					if (start_virtual_crossing_flag[i])
						num_initial_start_virtual_crossings++;
				}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: num_initial_start_virtual_crossings (includes peer path intersections) = " << num_initial_start_virtual_crossings << endl;

				for (int i=0; i< num_crossings; i++)
				{
					if (start_virtual_crossing_flag[i] && peer_virtual_crossing_flag[i])
					{
						start_virtual_crossing_flag[i] = 0;  // has the effect of removing this virtual crossing
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: crossing " << i << " lies on both the start and peer detours, crossing will be removed when peer detour is moved" << endl;
					}
				}
						
				int num_start_virtual_crossings = 0;
				for (int i=0; i< num_crossings; i++)
				{
					if (start_virtual_crossing_flag[i])
						num_start_virtual_crossings++;
				}

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "remove_Reidemeister_II: num_start_virtual_crossings (excluding peer path intersections) = " << num_start_virtual_crossings << endl;
	
				int new_num_crossings = num_crossings-2+num_start_virtual_crossings;
				
				if (!standard_Reidemeister_II_configuration)
				{
					new_num_crossings++;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: non-standard Reidemeister II condition, retaining rendezvous crossing" << endl;
				}
				
				for (int i=0; i< num_crossings; i++)
				{
					if (peer_virtual_crossing_flag[i])
						new_num_crossings--;
				}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: number of crossings remaining after moving the peer detour, new_num_crossings = " << new_num_crossings << endl;
			
				if (new_num_crossings == 0)
				{	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:  removing Reidemeister II move in different turning cycles reduces the diagram to as set of disjoint simple closed curves" << endl;
	
					/* we clear the component flags below in this case */
					code_data.num_crossings = 0;
				}				
				else
				{			
					/* identify the edges terminating at virtual crossings from the peer detour plus the two non-virtual crossings */
									
					if (standard_Reidemeister_II_configuration)
					{
						RII_edge_flag[code_table[generic_code_data::table::ODD_TERMINATING][start_crossing]] = 1;
						RII_edge_flag[code_table[generic_code_data::table::EVEN_TERMINATING][start_crossing]] = 1;
						RII_edge_flag[code_table[generic_code_data::table::ODD_TERMINATING][rendezvous_crossing]] = 1;
						RII_edge_flag[code_table[generic_code_data::table::EVEN_TERMINATING][rendezvous_crossing]] = 1;
					}
					else
					{
						if (forwards_from_peer)
						{
							RII_edge_flag[code_table[generic_code_data::table::ODD_TERMINATING][start_crossing]] = 1;
							RII_edge_flag[code_table[generic_code_data::table::EVEN_TERMINATING][start_crossing]] = 1;
							
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: non-standard Reidemeister II condition, forwards from peer, retaining rendezvous crossing terminating edges ";
	debug << code_table[generic_code_data::table::ODD_TERMINATING][rendezvous_crossing] << " and " << code_table[generic_code_data::table::EVEN_TERMINATING][rendezvous_crossing] << endl;
}
						}
						else
						{
							RII_edge_flag[code_table[generic_code_data::table::ODD_TERMINATING][rendezvous_crossing]] = 1;
							RII_edge_flag[code_table[generic_code_data::table::EVEN_TERMINATING][rendezvous_crossing]] = 1;
							
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: non-standard Reidemeister II condition, backwards from peer, retaining start crossing terminating edges ";
	debug << code_table[generic_code_data::table::ODD_TERMINATING][start_crossing] << " and " << code_table[generic_code_data::table::EVEN_TERMINATING][start_crossing] << endl;
}
						}
					}
					
					for (int i=0; i< num_crossings; i++)
					{
						if (peer_virtual_crossing_flag[i])
						{
							RII_edge_flag[code_table[generic_code_data::table::ODD_TERMINATING][i]] = 1;
							RII_edge_flag[code_table[generic_code_data::table::EVEN_TERMINATING][i]] = 1;
						}
					}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: removing RII_edge_flags ";
	for (int i=0; i< 2*num_crossings; i++)
		debug << RII_edge_flag[i] << ' ';
	debug << endl;
}
	
					
					/* identify whether the peer detour stays to the right or the left of the start detour 
					   (with respect to the start detour's orientation) after it has been moved 
					*/
					int start_polarity;				
					if ((start%2 == 0 && code_table[generic_code_data::table::TYPE][start_crossing] == generic_code_data::TYPE1) || 
					    (start%2 == 1 && code_table[generic_code_data::table::TYPE][start_crossing] == generic_code_data::TYPE2)
					   )
					{
						start_polarity = gauss_orientation_data::LEFT;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: start_crossing is a LEFT crossing from start edge perspective" << endl;
					}
					else
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: start_crossing is a RIGHT crossing from start edge perspective" << endl;
						start_polarity = gauss_orientation_data::RIGHT;
					}
					
					bool stay_left_of_start_detour;
					if (standard_Reidemeister_II_configuration)
					{
						if ((forwards_from_peer && start_polarity == gauss_orientation_data::LEFT) ||
							(!forwards_from_peer && start_polarity == gauss_orientation_data::RIGHT)
						   )
						{
							stay_left_of_start_detour = true;
						}
						else
						{
							stay_left_of_start_detour = false;
						}
					}
					else
					{
						/* for non-standard Reidemeister II configurations, the moved peer detour crosses the start 
						   detour at the start crossing but for it's passage past any transverse edges it stays left 
						   or right of the start detour dependent only on the polarity of the start crossing.   
						*/
						if (start_polarity == gauss_orientation_data::LEFT) 
						{
							stay_left_of_start_detour = true;
						}
						else
						{
							stay_left_of_start_detour = false;
						}
					}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   stay_left_of_start_detour = " << stay_left_of_start_detour << endl;
	
					/* Re-number the diagram missing out those edges that terminate at the virtual crossings we are removing
					   but allowing for the additional virtual crossings when the peer detour is moved. For those crossings 
					   from the original immersion that remain after the detour is moved, this renumbering provides the 
					   edge number of the terminating edges in the final state of the immersion.
					   
					   We write the new edge labels into a 2 x (num_crossings+num_start_virtual_crossings) matrix, new_edge_labels, 
					   that mimics the terminating edge data in the code_table in such a way that, if old edge e is renumbered e' 
					   we write e' in the same position in new_edge_labels as e appears in the terminating edge rows of the code data.				   
					   the terminating edges of the new crossings are added to the last num_start_virtual_crossings columns of
					   new_edge_labels later.
					   
					   As we write the new labels, we keep track of the new first and last edges on each component, which we shall use
					   for controlling any component renumbering that may be required later.  					   
					*/
					int num_new_edge_label_columns = num_crossings+num_start_virtual_crossings;
					vector<int> new_first_component_edge(num_components);  // correctly initialzes new_first_component_edge[0]=0
					vector<int> new_last_component_edge(num_components);

					matrix<int> initial_new_zig_zag_count(2,num_crossings);
					int amalgamated_zig_zag_count = 0;					
					
					matrix<int> new_edge_labels(2,num_new_edge_label_columns,-1);
					
					/* The new edge labels will be tracked by new_edge, which will be incremented for each edge that remains in the 
					   immersion, i.e. for those not identified by RII_edge_flag.  In order to adjust new_edge to accommodate (i.e skip 
					   over initially) those edges that will be created by moving the peer detour, we need to identify (i) an edge (that we 
					   refer to as the peer_detour_skip_edge) in the original immersion numbering at either the start or rendezvous 
					   crossing on the peer detour, and (ii) those edges transverse to the start detour. These are the edges before which 
					   additional crossings are added.
					   
					   The edge at the start or rendezvous crossings we are interested in depends on whether we are moving forwards or 
					   backwards from the peer edge when moving along the peer detour from the start crossing.  If we are moving forwards 
					   we want the edge following the rendezvous crossing on the peer detour.  If we are moving backwards we want the peer's
					   successor at the start crossing.  These edges are identified as edge A and B respectively below.
					   
					   
					                            ^             ^                                     ^     ^
					                      s\   /         \   /                        s\   /p        \   /
					                        \ /           \ /                           \ /           \ / 
					                         x             x                             x             x
					                        / \           / \                           / \           / \
					                      p/   \         /   \ edge A           edge B /   \         /   \
					                            v             v                       v     v
					                  
					   When we encounter edge A or B we increment new_edge by num_start_virtual_crossings to account for moving the peer detour.
					   
					   The transverse edges are dependent on the LEFT or RIGHT polarity of the virtual crossings on the start detour, which is 
					   recorded for each virtual crossing of the start detour in start_virtual_polarity.  The transverse edges we are interested
					   in are shown as edge C and edge D respectively below.
					   
	                                              				                                        
	                                              ^				                                        ^
	                                              |				                                        |
	                                              +				                                        +
	                                    |         | edge D		                              |         |
	                                    |         |                                           |         |
		                             ================= peer detour moves to here	   -------o---------o------->  start detour    
					                    |         |	                                          |         | edge D
					             edge C |         |				                              |         |
					             -------o---------o------->  start detour                  ================= peer detour moves to here
					                    |         |				                              |         |
					                    |         |				                       edge C |         |
					                    v          				                              +         
					                                				                          |          
	                                                                                          v
	
					   When we encounter an edge C or D we increment new_edge by one to account for the transverse path crossing the  
					   peer detour after it has been moved.				            

					   When we move the peer detour, if the transverse skip edge is either edge zeo , or the head_semi_arc, then the 
					   new crossing matches that on the start detour.  That is, if the transverse skip edge is zero, then case C on 
					   the lhs and D on the rhs introduce a virtual crossing, whereas case D on the lhs and C on the rhs introduce a 
					   shortcut crossing.  Similarly, if the transverse skip edge is the head_semi_arc, then then case C on the lhs 
					   and D on the rhs introduce a shortcut crossing, whereas case D on the lhs and C on the rhs introduce a virtual
					   crossing.
					*/
					int peer_detour_skip_edge;
					
					if (forwards_from_peer)
					{
						peer_detour_skip_edge = (peer_path_rendezvous_edge+1 - first_edge_on_component[peer_component])%
													num_component_edges[peer_component] + first_edge_on_component[peer_component];
					}
					else
					{
						peer_detour_skip_edge = (peer+1 - first_edge_on_component[peer_component])%
													num_component_edges[peer_component] + first_edge_on_component[peer_component];
					}

					if (!standard_Reidemeister_II_configuration)
					{
						/* the peer_detour_skip_edge needs moving back one edge within the peer component to account for the 
						   fact that we are retaining the crossing at the end of the peer path
						*/
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: non-standard Reidemeister II configuration, moving the initial peer_detour_skip_edge = " 
	      << peer_detour_skip_edge << " back one edge" << endl;
}
						peer_detour_skip_edge = (peer_detour_skip_edge-1+num_component_edges[peer_component]- first_edge_on_component[peer_component])%num_component_edges[peer_component] + first_edge_on_component[peer_component];
					}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: peer_detour_skip_edge = " << peer_detour_skip_edge << endl;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: evaluate start_virtual_polarity..." << endl;
		
					vector<int> start_virtual_polarity(num_start_virtual_crossings);			
					int start_virtual_index = 0;
					for (int i=0; i< num_initial_start_virtual_crossings; i++)
					{
						edge = (start+1+i - first_edge_on_component[start_component])%num_component_edges[start_component] + first_edge_on_component[start_component];

						/* check that this edge terminates at a crossing that remains on the start detour when the peer
						   detour has moved.  If the peer detour intersects the start detour, we may lose crossings from
						   the start detour.
						*/
						if (start_virtual_crossing_flag[term_crossing[edge]] == 1)
						{	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   start detour edge " << edge << " terminates at crossing " << term_crossing[edge] << ", which is " << 
	         (code_table[generic_code_data::table::TYPE][term_crossing[edge]] == generic_code_data::TYPE1? "TYPE1" : "TYPE2") << endl;
}					
							if ((edge%2 == 0 && code_table[generic_code_data::table::TYPE][term_crossing[edge]] == generic_code_data::TYPE1) ||
							    (edge%2 == 1 && code_table[generic_code_data::table::TYPE][term_crossing[edge]] == generic_code_data::TYPE2)
							   )
							{
								start_virtual_polarity[start_virtual_index] = gauss_orientation_data::LEFT;
							}
							else
							{
								start_virtual_polarity[start_virtual_index] = gauss_orientation_data::RIGHT;
							}
							
							start_virtual_index++;
						}
						else
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   edge " << edge << " terminates at a virtual crossing that lies on the peer detour, so has been removed" << endl; 
						}
						
					}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: start_virtual_polarity ";
	for (int i=0; i< num_start_virtual_crossings; i++)
		debug << (start_virtual_polarity[i] == gauss_orientation_data::LEFT? "L " : "R ");
	debug << endl;
}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: evaluate transverse_skip_edge..." << endl;
				
					vector<int> transverse_skip_edge(num_start_virtual_crossings);
					vector<int> transverse_incident_edge(num_start_virtual_crossings);

					bool edge_zero_a_transverse_skip_edge = false;
					start_virtual_index = 0;
					for (int i=0; i< num_initial_start_virtual_crossings; i++)
					{
						edge = (start+1+i - first_edge_on_component[start_component])%num_component_edges[start_component] + first_edge_on_component[start_component];
						
						/* check that this edge terminates at a crossing that remains on the start detour when the peer
						   detour has moved.  If the peer detour intersects the start detour, we may lose crossings from
						   the start detour.
						*/
						if (start_virtual_crossing_flag[term_crossing[edge]] == 1)
						{						
							int transverse_peer = (edge%2? code_table[generic_code_data::table::EPEER][(edge-1)/2]:code_table[generic_code_data::table::OPEER][edge/2]);
							int transverse_peer_component = (transverse_peer%2? code_table[generic_code_data::table::COMPONENT][(transverse_peer-1)/2]:code_table[generic_code_data::table::COMPONENT][(transverse_peer)/2]);
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   start detour edge " << edge << " has peer " << transverse_peer << " on component " << transverse_peer_component 
	      << ", start_virtual_index = " << start_virtual_index << endl; 
}	
							transverse_incident_edge[start_virtual_index] = transverse_peer;
					
							if ((stay_left_of_start_detour && start_virtual_polarity[start_virtual_index] == gauss_orientation_data::LEFT) ||
		                        (!stay_left_of_start_detour && start_virtual_polarity[start_virtual_index] == gauss_orientation_data::RIGHT)
		                       )
		                    {
								transverse_skip_edge[start_virtual_index] = transverse_peer;
							}
							else
							{
								transverse_skip_edge[start_virtual_index] = (transverse_peer+1 - first_edge_on_component[transverse_peer_component])%num_component_edges[transverse_peer_component] 
								                                             + first_edge_on_component[transverse_peer_component];																			 
							}
							
							if (transverse_skip_edge[start_virtual_index] == 0)
								edge_zero_a_transverse_skip_edge = true;
							
							start_virtual_index++;
						}
						else
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   edge " << edge << " terminates at a virtual crossing that lies on the peer detour, so has been removed" << endl; 
						}
					}
		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: transverse_skip_edge ";
	for (int i=0; i< num_start_virtual_crossings; i++)
		debug << transverse_skip_edge[i] << ' ';
	debug << endl;
	debug << "remove_Reidemeister_II: transverse_incident_edge ";
	for (int i=0; i< num_start_virtual_crossings; i++)
		debug << transverse_incident_edge[i] << ' ';
	debug << endl;
}
	
					/* It is possible, e.g. edge 10 in [-15 19 13 -11 -3,1 17 -9 5 -7]/# # # # # * * * * *, that the peer_detour_skip_edge is
					   also in transverse_skip_edge.  (In such cases, removing the Reidemeister II detour leaves a virtual Reidemeister I
					   move.)  Depending on the edge numbering, we may need to adjust the new edge labels first as a result of skipping over 
					   the moved peer detour and then skipping transversely back over it, or first skipping a transverse skip edge and then 
					   as a peer detour skip edge.  
					  					   
					   We identify which adjustment is needed first by comparing the peer_detour_skip_edge with peer_detour_entry_edge: 
					   these edges are on the same component and we will encounter the lower of the edge labels first.  We shall use a 
					   boolean, adjust_as_transverse_skip_first as a flag to indicate when we need to adjust new_edge as a transverse skip edge 
					   first, (since the check for peer_detour_skip_edge comes first in the loop clause below).
					   
					   If we adjust_as_transverse_skip_first then we have to adjust the new edge labels (to allow for the additional crossings when 
					   moving the peer detour) when we come to the end of the component containing the peer_detour_skip_edge.  This is to ensure that
					   the next component new edge label numbering starts with the right edge label.  This will only be the case if the edge label
					   numbering of the component containing the peer detour starts at the peer_detour_skip edge.					   
					*/
					bool adjust_as_transverse_skip_first = false;
					
					vector<int>::iterator peer_detour_skip_edge_ptr = find(transverse_skip_edge.begin(), transverse_skip_edge.end(), peer_detour_skip_edge);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	if (peer_detour_skip_edge_ptr != transverse_skip_edge.end())
		debug << "remove_Reidemeister_II: peer_detour_skip_edge found in transverse_skip_edge" << endl;				 
	else
		debug << "remove_Reidemeister_II: peer_detour_skip_edge not found in transverse_skip_edge" << endl;				 	
}

					if (peer_detour_skip_edge_ptr != transverse_skip_edge.end() && peer_detour_skip_edge < peer_detour_entry_edge)
					{
						adjust_as_transverse_skip_first = true;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: peer_detour_skip_edge is encountered as a transverse skip edge before we reach the peer detour entry edge" << endl;				 
					}
					else
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: peer_detour_skip_edge either not a transverse skip edge or is encountered first on the component as the peer detour skip edge" << endl;				
					}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: adjust_as_transverse_skip_first = " << adjust_as_transverse_skip_first << endl;
		
					/* We are now ready to write the new_edge_labels, we start with those that have been retained from the original 
					   immersion and then add in those new edges resulting from moving the peer path.
					
					   For knotoids, we also have to keep track of the head semi-arc.  This may lie:
					   
					       1. within the peer path and will be moved along with the peer detour
						   2. in a transverse skip edge; that is, in a strand incident with the start path 
						      that intersects the peer path after it has been moved
						   3. in a transverse strand incident with the peer path before it is moved
						   4. within the start path
						   5. not incident with (i.e not terminating at a crossing within) the start path or the peer path

					   In cases 3,4 and 5, the new_head_semi_arc will simply be the new_edge corresponding to the original 
					   head_semi_arc.  Note that in case 3, the peer path is moved to the other side of the start path to 
					   that containing the head semi-arc (even if the Reidemeister II condition is not a standard 
					   Reidemeister II), so if the second shortcut crossings lies in the start path, the ingress shortcut 
					   edge at this crossing will be retained and becomes the new_head_semi_arc.  Note that in case 4, 
					   the head_semi_arc may terminate at one of the Reidemeister crossings, in which case it is
					   an edge that will be removed, but this approach remains correct.
					   
					   In case 2, the new_head_semi_arc will be the value of new_edge before it is incremented as a result of 
					   encountering the transverse skip edge, so again is the new_edge corresponding to the original 
					   head_semi_arc.
					   
					   In case 1, any peer path virtual crossings may be removed from the segment component so that in the 
					   original diagram the head semi_arc is the bigon end edge.  Thus when the peer path is moved, the 
					   new head_semi_arc is the first_new_peer_detour_edge (see below), since we shall always approach 
					   the head following the orientation of the segment component.
					*/
					int new_edge = 0;
					int new_edge_component = 0;
					
					vector<int> new_transverse_peer_edge(num_start_virtual_crossings);
					
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: renumbering immersion, new_edge initialized to " << new_edge << endl;				 
							
					int new_head_semi_arc = -1;
					int first_new_peer_detour_edge = -1; // set to aid debugging

					for (int old_edge = 0; old_edge < 2*num_crossings; old_edge++)
					{
						if (pure_knotoid_code_data && old_edge == head_semi_arc)
						{
							new_head_semi_arc = new_edge;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   encountered head_semi_arc at old_edge " << old_edge << " , setting  new_head_semi_arc = " << new_head_semi_arc << endl;				 
						}
												
						int old_edge_component;
						if (old_edge%2)
							old_edge_component = code_table[generic_code_data::table::COMPONENT][(old_edge-1)/2];
						else
							old_edge_component = code_table[generic_code_data::table::COMPONENT][(old_edge)/2];
							
						/* Since the start and rendezvous crossings are in different turning cycles, there is at least one
						   virtual crossing on the start detour.  This means there will always be at least one crossing 
						   remaining on every component, after moving the peer detour.  This in turn means that we may always
						   set the new_first_edge to be new_edge when we detect a change in component via old_edge.
						*/
						if (old_edge_component != new_edge_component)
						{							
							if (adjust_as_transverse_skip_first && peer_component == new_edge_component)
							{
								/* If we need to adjust the peer_detour_skip edge as a transverse skip edge first, then we have to adjust 
								   new_edge when we come to the end of the component to account for the new start virtual crossings.
								   
								   In the case of a knotoid, if adjust_as_transverse_skip_first is true (i.e the peer_detour_skip_edge
								   is less than the peer_detour_entry edge), then the leg of the knotoid is encountered after the 
								   peer_detour_entry_edge as we cycle the first component.  In such a case, if the peer_detour_skip
								   edge were the head_semi_arc then we would have a floating first component, which would already have 
								   been removed.  Therefore, we do not need to check for any adjustment to the head_semi_arc here, as we do
								   when adjust_as_transverse_skip_first is false and we increment new_edge by num_start_virtual_crossings
								*/
																
								new_edge += num_start_virtual_crossings;
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   reached end of component containing peer detour skip edge without adjusting for  num_start_virtual_crossings, new_edge incremented to " << new_edge << endl;				 
}
							}
							
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "remove_Reidemeister_II:   reached end of component " << new_edge_component << endl;				 
							new_last_component_edge[new_edge_component] = new_edge-1;
							
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: amalgamated_zig_zag_count = " << amalgamated_zig_zag_count << " at the end of component " << new_edge_component << endl;
	
							if (code_data.immersion_character != generic_code_data::character::CLOSED && new_edge_component == 0)
							{
								/* set, or amalgamate with, the head_zig_zag_count */
								if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
								{
									code_data.head_zig_zag_count = amalgamated_zig_zag_count;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   pure knotoid segment component, code_data.head_zig_zag_count set to " << code_data.head_zig_zag_count << endl;
								}
								else
								{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   knot-type knotoid or long knot segment component, amalgamated with code_data.head_zig_zag_count " << code_data.head_zig_zag_count;
										code_data.head_zig_zag_count = amalgamate_zig_zag_counts(amalgamated_zig_zag_count, code_data.head_zig_zag_count);
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " to give " << code_data.head_zig_zag_count << endl;
								}
								
								/* update head semi-arc zig-zag count */
								if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
									initial_new_zig_zag_count[(head_semi_arc %2?1:0)][code_data.head] = code_data.head_zig_zag_count;
							}
							else
							{				
								/* amalgamate with the count associated with the first odd crossing on the component, note that this may be a component
								   that is being removed, in which case we'll not find any edge assigned a new label and the component will be identified
								   as being removed later, when we consider all components.
								*/
								bool edge_found = false;
								
								for (int j=0; j< num_component_edges[new_edge_component]; j++)
								{
									int edge = first_edge_on_component[new_edge_component]+j;
									int crossing = code_data.term_crossing[edge];
									int row = (edge%2 == 0?0:1);
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   old edge " << edge << " terminates at crossing " << crossing << " row = " << row << endl;
					
									if (code_table[generic_code_data::table::LABEL][crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[crossing]))
									{
										if (new_edge_labels[row][crossing] != -1)
										{		
											
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   loop component, amalgamating with count " << initial_new_zig_zag_count[row][crossing] << " for old edge " 
	      << edge << ": row = " << row << " crossing = " << crossing << endl;
}
											initial_new_zig_zag_count[row][crossing] = amalgamate_zig_zag_counts(amalgamated_zig_zag_count, initial_new_zig_zag_count[row][crossing]);
										
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   count updated to " << initial_new_zig_zag_count[row][crossing] << endl;
											
											edge_found = true;
											break;
										}
									}				
								}
								
								if (!edge_found)
								{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   loop component, no edge found, component being removed" << endl;
								}					
							}
							amalgamated_zig_zag_count = 0;
							
							new_edge_component++;
							new_first_component_edge[new_edge_component] = new_edge;
						}
						
						if (old_edge == peer_detour_skip_edge)
						{
							if (adjust_as_transverse_skip_first)
							{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   peer_detour_skip_edge " << old_edge << " is in transverse_skip_edge and is encountered as a transverse edge first" << endl;				 
							}
							else
							{
								/* In cases such as [-3 5 -7 1^]/# # + -, the head semi-arc is 6, which is also the peer detour 
								   skip edge when we come to moving the peer detour.  That means we need to increment the new_head_semi_arc
								   that we have already identified to accommodate the start virtual crossings.
								   
								   Note that we may have removed the old head_semi_arc and are about to reset the new_head_semi_arc to zero,
								   as in the case of [9 -7 1, -3 5^]/# # * * +, so we have to check that we're still on new component zero
								   before adjusting new_head_semi_arc.
								*/
								bool incremented_new_head_semi_arc = false;
								
								if (new_edge_component == 0 && new_edge == new_head_semi_arc)
								{
									new_head_semi_arc += num_start_virtual_crossings;
									incremented_new_head_semi_arc = true;
								}
								
								new_edge += num_start_virtual_crossings;
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   old_edge " << old_edge << " is peer_detour_skip_edge, new_edge incremented to " << new_edge << endl;				 
	if (incremented_new_head_semi_arc)
		debug << "remove_Reidemeister_II:   incremented new_head_semi_arc to " << new_head_semi_arc << endl;				 
}								
							}
						}
	
						vector<int>::iterator eptr = find(transverse_incident_edge.begin(), transverse_incident_edge.end(), old_edge);
						bool increment_new_edge_twice_after_assignment = false;
						
						if (eptr != transverse_incident_edge.end())
						{

							int offset = eptr - transverse_incident_edge.begin();
							
							if ((stay_left_of_start_detour && start_virtual_polarity[offset] == gauss_orientation_data::LEFT) ||
		                        (!stay_left_of_start_detour && start_virtual_polarity[offset] == gauss_orientation_data::RIGHT)
		                       )
							{
								new_transverse_peer_edge[offset] = new_edge;
								new_edge++; // before assignment of the new label for the incident edge
								increment_new_edge_twice_after_assignment = false;
							}
							else
							{
								new_transverse_peer_edge[offset] = new_edge+1;
								increment_new_edge_twice_after_assignment = true;
							}
																			
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   old_edge " << old_edge << " is in transverse_incident_edge, at offset " << offset << endl;
	debug << "remove_Reidemeister_II:   new_transverse_peer_edge = " << new_transverse_peer_edge[offset] 
	      <<  ", before new label assignment, new_edge stands at " << new_edge << endl;				 
}	
						}

						if (old_edge == peer_detour_entry_edge)
						{
							first_new_peer_detour_edge = new_edge;
							
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   first_new_peer_detour_edge = " << first_new_peer_detour_edge << endl;
						}
						
						if (track_zig_zag_counts)
						{
							/* amalgamate any zig-zags from old_edge with the current amalgamated_zig_zag_count */
							int crossing = term_crossing[old_edge];
							int old_edge_zig_zag_count = (old_edge %2 == 1? code_data.zig_zag_count[1][crossing]:code_data.zig_zag_count[0][crossing]);
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: old_edge " << old_edge << " terminates at crossing " << crossing << ", old_edge_zig_zag_count = " << old_edge_zig_zag_count << endl;				 
			
							amalgamated_zig_zag_count = amalgamate_zig_zag_counts(amalgamated_zig_zag_count, old_edge_zig_zag_count);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: amalgamated_zig_zag_count adjusted to " << amalgamated_zig_zag_count << endl;				 
						}
						
						if (RII_edge_flag[old_edge] == 0)
						{
							/* write new_edge into new_edge_labels at the location corresponding to 
							   old_edge's place in the terminating edge data contained in the code table.  We
							   know which row to look in because of the polarity of old_edge.
							*/
							int row = (old_edge%2 ? generic_code_data::table::ODD_TERMINATING: generic_code_data::table::EVEN_TERMINATING);
							int col=0;
							for (int i=0; i< num_crossings; i++)
							{
								if (code_table[row][i] == old_edge)
								{
									col = i;
									break;
								}
							}
							
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   old_edge " << old_edge << " lies in position " << (row==generic_code_data::table::ODD_TERMINATING? "generic_code_data::table::ODD_TERMINATING" : "generic_code_data::table::EVEN_TERMINATING") 
	      << ", " << col << "; new_edge = " << new_edge << endl;				 
}					
							/* we'll write the new edges corresponding to the old even edges into the first row of new_edge_labels 
							   and the new edges corresponding to the old odd edges into the second row */
							if (row == generic_code_data::table::ODD_TERMINATING)
								row = 1;
							else
								row = 0;
								
							new_edge_labels[row][col] = new_edge;																	
							new_edge++;

							if (increment_new_edge_twice_after_assignment)
							{
								new_edge++;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   new_edge incremented again after assignment to " << new_edge << endl;				 
							}
							
							if (track_zig_zag_counts)
							{
								int crossing = term_crossing[edge];
								
								if (code_table[generic_code_data::table::LABEL][crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[crossing]))
								{					
									initial_new_zig_zag_count[row][col] = amalgamated_zig_zag_count;
									amalgamated_zig_zag_count = 0;				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   new_edge " << new_edge << " initial_new_zig_zag_count[" << row << "][" << col << "] = " << initial_new_zig_zag_count[row][col] << endl;
								}
							}
						}
						else
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   old_edge " << old_edge << " terminates at a Reidemeister detour, or peer virtual, crossing, new_edge remains " << new_edge << endl;
						}	
					}									
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: amalgamated_zig_zag_count = " << amalgamated_zig_zag_count << " at the end of component " << new_edge_component << endl;
	
					if (code_data.immersion_character != generic_code_data::character::CLOSED && new_edge_component == 0)
					{
						/* set, or amalgamate with, the head_zig_zag_count */
						if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
						{
							code_data.head_zig_zag_count = amalgamated_zig_zag_count;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   pure knotoid segment component, code_data.head_zig_zag_count set to " << code_data.head_zig_zag_count << endl;
						}
						else
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   knot-type knotoid or long knot segment component, amalgamated with code_data.head_zig_zag_count " << code_data.head_zig_zag_count;
								code_data.head_zig_zag_count = amalgamate_zig_zag_counts(amalgamated_zig_zag_count, code_data.head_zig_zag_count);
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " to give " << code_data.head_zig_zag_count << endl;
						}
			
						/* update head semi-arc zig-zag count */
						if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
							initial_new_zig_zag_count[(head_semi_arc %2?1:0)][code_data.head] = code_data.head_zig_zag_count;
					}
					else
					{				
						/* amalgamate with the count associated with the first odd crossing on the component, note that this may be a component
						   that is being removed, in which case we'll not find any edge assigned a new label and the component will be identified
						   as being removed later, when we consider all components.
						*/
						bool edge_found = false;
						
						for (int j=0; j< num_component_edges[new_edge_component]; j++)
						{
							int edge = first_edge_on_component[new_edge_component]+j;
							int crossing = code_data.term_crossing[edge];
							int row = (edge%2 == 0?0:1);
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   old edge " << edge << " terminates at crossing " << crossing << " row = " << row << endl;
					
							if (code_table[generic_code_data::table::LABEL][crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[crossing]))
							{
								if (new_edge_labels[row][crossing] != -1)
								{		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   loop component, amalgamating with count " << initial_new_zig_zag_count[row][crossing] << " for old edge " 
	      << edge << ": row = " << row << " crossing = " << crossing << endl;
}
									initial_new_zig_zag_count[row][crossing] = amalgamate_zig_zag_counts(amalgamated_zig_zag_count, initial_new_zig_zag_count[row][crossing]);
										
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   count updated to " << initial_new_zig_zag_count[row][crossing] << endl;
	
									edge_found = true;
									break;
								}
							}				
						}
						
						if (!edge_found)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   loop component, no edge found, component being removed" << endl;
						}					
					}
					amalgamated_zig_zag_count = 0;
							
					/* complete the assignment of new_last_component_edge, we can do this even though we 
					   haven't added in the new edges corresponding to moving the peer detour, since
					   we may calculate it directly.
					*/
					new_last_component_edge[new_edge_component] = 2*new_num_crossings -1;
					
					/* If the peer_detour_entry_edge becomes the first edge in its component after the peer
					   detour is moved, as in [-5 -23 13 -1 -21^ 3, -11 -19 7 -9 -15 17]/ -*----+--+-*, then the
					   above assignment of the first_new_peer_detour_edge will be incorrect, since new_edge had
					   been incremented in anticipation of another edge on this component that does not exist.
					   We therefore adjust first_new_peer_detour_edge to new_first_component_edge[peer_component]
					   in this case.
					*/
					if (first_new_peer_detour_edge > new_last_component_edge[peer_component])
					{
						first_new_peer_detour_edge = new_first_component_edge[peer_component];
						
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: adjusted first_new_peer_detour_edge to " << first_new_peer_detour_edge << " to accommodate component edge label wrap " << endl;
					}
					
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: new_first_component_edge ";
	for (int i=0; i< num_components; i++)
		debug << new_first_component_edge[i] << ' ';
	debug << "\nremove_Reidemeister_II: new_last_component_edge ";
	for (int i=0; i< num_components; i++)
		debug << new_last_component_edge[i] << ' ';
	debug << "\nremove_Reidemeister_II: new_transverse_peer_edge ";
	for (int i=0; i< num_start_virtual_crossings; i++)
		debug << new_transverse_peer_edge[i] << ' ';
	debug << endl;
}	

					if (pure_knotoid_code_data)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: knotoid: after assigning new edge labels new_head_semi_arc = " << new_head_semi_arc << endl;
	
						new_head_semi_arc %= (new_last_component_edge[0]+1); 
						
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: knotoid: reset new_head_semi_arc = " << new_head_semi_arc << endl;			
					}
						
					/* now add in the new edges resulting from moving the peer detour.  They are recorded in the order
					   determined by increasing value of new_edge (e.g. if moving backwards from peer, for increasing i, 
					   the columns represent the crossings in oder moving from the rendezvous crossing to the start crossing
					*/
					new_edge = first_new_peer_detour_edge;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: renumbering new crossings resulting from moving peer detour, first_new_peer_detour_edge = " << new_edge << endl;			
					for (int i=0; i< num_start_virtual_crossings; i++)
					{
						int transverse_virtual_crossing;
						if (forwards_from_peer)
							transverse_virtual_crossing = i;
						else
							transverse_virtual_crossing = num_start_virtual_crossings-1-i;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: new edge " << new_edge << ", transverse_virtual_crossing = " << transverse_virtual_crossing 
	      << " peer edge " << new_transverse_peer_edge[transverse_virtual_crossing] << endl;				 
}
						
						if (new_edge%2)
						{
							new_edge_labels[0][num_crossings+i] = new_transverse_peer_edge[transverse_virtual_crossing];
							new_edge_labels[1][num_crossings+i] = new_edge;
						}
						else
						{
							new_edge_labels[0][num_crossings+i] = new_edge;
							new_edge_labels[1][num_crossings+i] = new_transverse_peer_edge[transverse_virtual_crossing];
						}
						new_edge++;
					}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: old edge_labels" << "\nremove_Reidemeister_II: ";
	for (int i=0; i< num_crossings; i++)
		debug << setw(3) << code_table[generic_code_data::table::EVEN_TERMINATING][i];
	debug << "\nremove_Reidemeister_II: ";
	for (int i=0; i< num_crossings; i++)
		debug << setw(3) << code_table[generic_code_data::table::ODD_TERMINATING][i];
	debug << endl;
	debug << "remove_Reidemeister_II: new_edge_labels: " << endl;
	print(new_edge_labels, debug, 3, "remove_Reidemeister_II: ");
	
	if (track_zig_zag_counts)
	{
		debug << "remove_Reidemeister_II: initial_new_zig_zag_count: " << endl;
		print(initial_new_zig_zag_count,debug, 4,"remove_Reidemeister_II: ");	
	}
	
}

					/* If we're dealing with an open diagram, we need to ensure that the leg of the segment component is labelled with zero but it may 
					   be the case that the old edge zero is removed, or moved (as part of the peer detour), or is assigned a non-zero new edge label.
					   
					   The value of new_edge is incremented before assignment:
						 a) if the old_edge is the peer_detour_skip_edge
						 b) if the old_edge is a transverse_skip edge incident with the start path
											   
					   If old edge zero is the peer_detour_skip edge, but is not also a transverse skip edge, then it will be assigned a non-zero new edge 
					   label if the start and rendezvous crossings are in different turning cycles, since there is at least one crossing in the start detour 
					   in that case, we may use the new edge label to determine the required shift.  If old edge zero is also a transverse skip edge, then 
					   it will correctly be assigned label zero (as may be seen in the horrendously complex example:
					   [-21 -39 19 -17 45 -43 41 -37, -31 33 -35 -15 13 -11 9 -7^ 5 -3 -23 -1 25 -27 29]/ * + * * - - - - - * * + + + + + - + - + + - +   

					   If edge zero is a transverse skip edge that terminates on the start path, then the additional crossing introduced by moving the peer is
					   virtual and old edge zero is correctly incremented to new edge label 1.  
										   
					   If old_edge zero is the peer_detour_entry_edge or lies within the peer detour, it is removed (or moved parallel to the start path) but, 
					   in these cases, we can extend the shortcut to the peer_detour_skip_edge, since there are only virtual or shortcut crossings on the bigon 
					   boundary.  Moerover, the Reidemeister II crossings are removed, so if old edge zero is the peer_detour entry edge we only meet virtual 
					   crossings having moved the peer detour before we get to the peer_skip_edge.  If edge zero lies within the peer detour, then we have a 
					   virtual_Reidemeister_II_plus configuration, since Gauss Reidemeister II configurations are discounted if the shortcut lies within the bigon 
					   boundary, so it is permitted to move the leg along with the peer path, since all crossings are virtual or shortcut crossings.  Note that 
					   in this cases that we may have one virtual and one shortcut Reidemeister crossing.
					   
					   However, whilst we can extend the shortcut to the the peer_detour_skip_edge, that edge will likely have been given a non-zero new edge 
					   label zero.  This is because there is usually at least one crossing in the peer detour, since the start and rendezvous crossings are in 
					   different turning cycles.  Only in some non-standard Reidemeister II configurations, such as that arising in 
					   [13 -5 -11 9, 1 -3 7]/# * * * * # * do we see no start virtual crossings and yet the start and end bigon edges in different turning cycles.  
					   We therefore look for the peer_detour_skip_edge in the code_table and use the corresponding new_edge_label to determine the renumbering 
					   required to have the leg numbered zero.
					   
					   If old edge zero lies in the start detour, then either the edge is not removed (if it meets a transverse edge), or (in the case that it 
					   meets a start virtual crossing that is also in the peer detour, or it meets the rendezvous crossing) the first edge that is not removed 
					   is correctly assigned new edge label zero.
					   
					   If edge zero is the start edge and is not also a transverse skip edge then it is removed but the leg of the segment component is still 
					   numbered zero, so we need do nothing.  
					   
					   When dealing with with knot-type knotoids or long knots, it is possible that edge zero is both the start edge and a transverse skip edge.  
					   For example,  in K[-7 5 9, -1 3]/# * # * *, there is a Gauss Reidemeiseter II move with start = 0, peer = 7, rendesvous edges 4 (on the 
					   start path) and 9 (on the peer path) we have a standard Reidemeister II move where the start path loops around the rendezvous crossing 
					   and intersects the peer path in a virtual crossing before reaching the rendezvous crossing.
					   
								 _____________________
								/     ______          |
						   \   /     |      \   /     |
							\ /      |       \ /      | 
							 x       |        x       |
							/ \      |       / \      |
						   /   \     |      /   \     |
						  s=0	-----o-----      |    |
									 |           |    |
									  -----------o-----
												 |					   
						
					   However, in such a case, edge zero does not terminate at a start virtual crossing (since these do not include the reidemeister II crossings), 
					   so it is not incremented before the assignment of the new label. 
					*/
										
					if (code_data.immersion_character != generic_code_data::character::CLOSED) 
					{					
						int component_zero_shift = 0;
						
						if (new_edge_labels[0][0] >= 0)  // >= just so we get the debug message even if edge zero remains labelled with zero
						{
														
							if (peer_detour_skip_edge == 0 && !edge_zero_a_transverse_skip_edge)
							{
								component_zero_shift = new_edge_labels[0][0];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: knotoid: old edge zero is the peer detour skip edge but not a transverse_skip_edge and is assigned new edge label = " << new_edge_labels[0][0] << endl;	
							}
							else
							{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: knotoid: old edge zero retained and assigned the correct new edge label " << new_edge_labels[0][0] << endl;	
							}
						}
						else if (new_edge_labels[0][0] == -1)
						{						
							if (peer_detour_entry_edge == 0 || leg_semi_arc_in_peer_path)
							{
								int row = (peer_detour_skip_edge%2 ? generic_code_data::table::ODD_TERMINATING: generic_code_data::table::EVEN_TERMINATING);
								int col=0;
								for (int i=0; i< num_crossings; i++)
								{
									if (code_table[row][i] == peer_detour_skip_edge)
									{
										col = i;
										break;
									}
								}						
								
								if (row == generic_code_data::table::ODD_TERMINATING)
									row = 1;
								else
									row = 0;
								
								component_zero_shift = new_edge_labels[row][col];
							
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: knotoid: edge zero not retained, peer_detour_skip_edge " << peer_detour_skip_edge << " found in code_table at row=" 
	      << row << " col=" << col << endl;	
}
							}
							else if (start == 0)
							{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: knotoid: old edge zero removed but edge zero still identifies the knotoid leg" << endl;	
							}							
							else
							{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: knotoid: edge zero not retained, but new edge label 0 correctly assigned to new leg" << endl;	
							}
						}
						
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   need to shift numbering of component zero by " << component_zero_shift << endl;	
						
						if (component_zero_shift)
						{
							int new_num_edges =new_last_component_edge[0]+1;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: knotoid: component zero new_num_edges = " << new_num_edges << endl;
									
							for (int j=0; j< num_new_edge_label_columns; j++)
							{
								if (new_edge_labels[0][j] != -1 && new_edge_labels[0][j] < new_num_edges)
									new_edge_labels[0][j] = (new_edge_labels[0][j] - component_zero_shift + new_num_edges)%new_num_edges;
									
								if (new_edge_labels[1][j] != -1 && new_edge_labels[1][j] < new_num_edges)
										new_edge_labels[1][j] = (new_edge_labels[1][j] - component_zero_shift + new_num_edges)%new_num_edges;
							}					


							new_head_semi_arc = (new_head_semi_arc - component_zero_shift + new_num_edges)%new_num_edges;
							if (peer_component == 0)
								first_new_peer_detour_edge = (first_new_peer_detour_edge - component_zero_shift + new_num_edges)%new_num_edges;
							
						
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: knotoid: after shifting component zero new_edge_labels: " << endl;
	print(new_edge_labels, debug,3,"remove_Reidemeister_II: ");
	debug << "remove_Reidemeister_II: knotoid: after shifting component zero new_head_semi_arc = " << new_head_semi_arc << endl;
	if (peer_component == 0)
		debug << "remove_Reidemeister_II: knotoid: after shifting component zero first_new_peer_detour_edge = " << first_new_peer_detour_edge << endl;
}
						}
					}
					
					/* Identify where in new_edge_labels the first_new_peer_detour_edge lies, so that we may adjust the new_head_semi_arc
					   when the head_semi_arc lies within the peer path (case 1 above)
					*/
					int first_new_peer_detour_edge_row;
					int first_new_peer_detour_edge_col;
					bool first_new_peer_detour_edge_found = false;
					for (int i=0; i< 2; i++)
					{
						for (int j=0; j < num_new_edge_label_columns; j++)
						{
							if (new_edge_labels[i][j] == first_new_peer_detour_edge)
							{
								first_new_peer_detour_edge_row = i;
								first_new_peer_detour_edge_col = j;
								first_new_peer_detour_edge_found = true;
								break;
							}
						}
						if (first_new_peer_detour_edge_found)
							break;
					}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: first_new_peer_detour_edge_row = " << first_new_peer_detour_edge_row 
	      << " first_new_peer_detour_edge_col = " << first_new_peer_detour_edge_col << endl;
}					

					/* For some peer codes of links, we may now have two even or two odd new edges terminating at the same crossing.
					   
					   If this is the case, we cycle the labels on individual components so that we have even and odd edges arriving 
					   at each crossing.  	   	   
					   
					   Note also that new_edge now indicates the number of edges in the reduced diagram, so we note first the new number of
					   crossings.
					*/
					align_label_numbers (new_edge_labels, num_components,new_first_component_edge,new_last_component_edge);
		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: final renumbered new_edge_labels: " << endl;
	print(new_edge_labels, debug,3,"remove_Reidemeister_II: ");
}

					/* identify any components removed from the peer code by removing the edges.  If there is such a component, 
					   none of its edges will have been assigned a new edge label
					*/
					for (int i=0; i< num_components; i++)
					{
						bool component_removed = true;
						for (int j=0; j< num_component_edges[i]; j++)
						{
							int edge = first_edge_on_component[i]+j;
							int crossing = code_data.term_crossing[edge];
							int row = (edge%2 == 0?0:1);
							
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: old edge " << edge << " terminates at crossing " << crossing << " row = " << row << endl;
				
							if (new_edge_labels[row][crossing] != -1)
							{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   new edge label " << new_edge_labels[row][crossing] << endl;
				
								component_removed = false;
								break;
							}				
						}
						
						if (component_removed)
						{
							clear_component_flag(i,component_flags);
							
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: removing edge flags disconnects component " << i << endl;
	debug << "remove_Reidemeister_II:   component_flags updated to: ";
	for (unsigned int j=0; j< component_flags.size(); j++)
		debug << component_flags[j] << ' ';
	debug << endl;
}
						}
						else
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   removing edge flags retains component " << i << endl;
						}
					}

					/* Adjust the new_head_semi_arc if the head_semi_arc lay within the peer path (case 1 above).  
					   When the peer path is moved, the new head_semi_arc is the first_new_peer_detour_edge, 
					   since we shall always approach the head following the orientation of the segment component.
					   
					   Since the components may have been renumbered, we have to use the position in new_edge_labels 
					   to identify the new_head_semi_arc.  Note that for non-standard Reidemeister II configurations,
					   we may have the start and rendezvous crossings in different turning cycles but still have
					   num_start_virtual_crossings = 0, as in the case of [11 -13^ -9 1 -3 15 -5 7]/- + * * - + - -.
					   That's why we need to track first_new_peer_detour_edge_row and first_new_peer_detour_edge_col. 
					*/
					if (pure_knotoid_code_data)
					{
						if (head_semi_arc_in_peer_path) // case 1
						{
							new_head_semi_arc = new_edge_labels[first_new_peer_detour_edge_row][first_new_peer_detour_edge_col];
							
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: knotoid: head semi_arc is in the peer path, setting head_semi_arc to new label on first_new_peer_detour_edge = " << first_new_peer_detour_edge << endl;
						}
						
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: knotoid: new_head_semi_arc = " << new_head_semi_arc << endl;
	
						/* If we have a knotoid that ends up with new_head_semi_arc = 0, as in the case:
						
						   [-15 -27 19 -29 -3 1 -23 25 -11 7 13^, -5 -21 -9 17]/ * + * - * * - + - * - - + * -
						
						   then the leg and the head are on the same semi-arc (zero), that is the knotoid ceases to 
						   be a pure knotoid and we shall therefore leave the head set to -1, identifying that the 
						   knotoid has reduced to a knot.  					   
						*/
						if (new_head_semi_arc == 0)
						{
							/* clear pure_knotoid_code_data so we assign labels correctly below, the head is set below */
							pure_knotoid_code_data = false;
							code_data.immersion_character = generic_code_data::character::KNOTOID;
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: new_head_semi_arc is zero, clearing knotoid status of code data" << endl;
	debug << "remove_Reidemeister_II: set code_data.immersion_character = generic_code_data::character::KNOTOID to indicate knot-type knotoid" << endl;
}
							
						}
						
					}
					
					/* to update the generic code data we create a new code_table from the existing one,
					   write the resultant generic code data to an ostringstream and then read it back in.
					   This avoids having to re-calculate the originating and terminating edges etc..  
					    
					   We must recompute the component assignment from the number of new edge labels on each
					   component, since the naming edge for a crossing may now be on a different component than 
					   before.  
					*/
					matrix<int> new_code_table (generic_code_data::table::CODE_TABLE_SIZE,new_num_crossings);
					matrix<int> new_zig_zag_count (2, new_num_crossings);

					for (int i=0; i< num_new_edge_label_columns; i++)
					{
						if (new_edge_labels[0][i] != -1) 
						{						
							int odd_edge = new_edge_labels[1][i];
							int even_edge = new_edge_labels[0][i];
							int odd_zig_zag_count = 0;
							int even_zig_zag_count = 0;
							
							/* the last entries in new_edge_labels are those relating to the moved peer path */
							if (track_zig_zag_counts && i < num_crossings) 
							{
								odd_zig_zag_count = initial_new_zig_zag_count[1][i];
								even_zig_zag_count = initial_new_zig_zag_count[0][i];
							}
							
							if (even_edge % 2 == 1)
							{
								swap(odd_edge,even_edge);
								swap(odd_zig_zag_count,even_zig_zag_count);
							}
							
							int col = even_edge/2;
							
							new_code_table[generic_code_data::table::OPEER][col] = odd_edge;
							new_code_table[generic_code_data::table::EPEER][(odd_edge-1)/2] = even_edge;
							
							/* We arbitrarily recorded the new_edge labels of the virtual crossings introduced by moving the peer
							   with the even edge in row zero, but they do not have a corresponding old crossing from which we can
							   take the type, so we defer handling the type and label for these crosings until we've assigned the 
							   rest of new_code_table.
                            */                        
                            if (i< num_crossings)
                            {    
								if (odd_edge == new_edge_labels[0][i])
								{
									/* in the new numbering, what was an even-numbered edge is now odd-numbered, therefore
									   the type and any classical label will be reversed 
									*/
									new_code_table[generic_code_data::table::TYPE][col] = (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1?generic_code_data::TYPE2:generic_code_data::TYPE1);
									
									if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE)
										new_code_table[generic_code_data::table::LABEL][col] = generic_code_data::NEGATIVE;
									else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
										new_code_table[generic_code_data::table::LABEL][col] = generic_code_data::POSITIVE;
									else
										new_code_table[generic_code_data::table::LABEL][col] = code_table[generic_code_data::table::LABEL][i];
								}
								else
								{
									new_code_table[generic_code_data::table::LABEL][col] = code_table[generic_code_data::table::LABEL][i];
									new_code_table[generic_code_data::table::TYPE][col] = code_table[generic_code_data::table::TYPE][i];
								}
							}

							new_code_table[generic_code_data::table::EVEN_TERMINATING][col] = even_edge;
							new_code_table[generic_code_data::table::ODD_TERMINATING][col] = odd_edge;

							if (track_zig_zag_counts)
							{
								new_zig_zag_count[1][col] = odd_zig_zag_count;
								new_zig_zag_count[0][col] = even_zig_zag_count;
							}
							
							int component = 0;
							for (int j=1; j < num_components; j++)
							{
								if (even_edge >= new_first_component_edge[j])
									component++;
								else
									break;
							}
							new_code_table[generic_code_data::table::COMPONENT][col] = component;
						}
					}
					
					/* The generic_code_data::table::TYPE of the virtual crossings and both the generic_code_data::table::TYPE and generic_code_data::table::LABEL of shortcut crossings
					   introduced by moving the peer detour may be determined from the generic_code_data::table::TYPE and generic_code_data::table::LABEL of 
					   the start virtual crossing on the same transverse path, since the start detour 
					   and moved peer detour are parallel.
					   
					   If moving forwards from the peer, the peer detour crossings have the opposite type
					   to those on the start crossing.  If moving backwards from the peer, the peer detour
                       crossings have the same type.  This is true regardless of whether the peer detour is
                       to the left or right of the start detour and the LEFT or RIGHT polarity of the transverse
                       path at the corresponding crossings.  Two examples are shown below.                                              				                                         
						
					   The label of shortcut crossings is always the opposite of the start crossing label
					   
						  o          e				                              ^         ^
							|         |       		                              |         |
							|I        |II                                         |II       |I
					  ------+---------+--------> peer/start               --------+---------+------>  start/peer      
		                e   |    o    |	                                      e   |    o    |       
		                   e|        o|				                             o|        e|
		                    |II       |I			                              |II       |I
		             -------+---------+------->  start/peer               <-------+---------+------ peer/start               
		                o   |    e    |				                              |    o    |   e
		                    |         |				                              |         |
		                    v         v				                            e          o	
		                    
		               We work forwards along the start detour, the number of edges to consider being determined by the initial
		               number of virtual crossings on the start detour (before any were removed as a result of also lying on
		               the peer detour).
					*/
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: set generic_code_data::table::TYPE and generic_code_data::table::LABEL of virtual crossings introduced by moving the peer detour" << endl;
	
					int index = 0;
					for (int i=0; i< num_initial_start_virtual_crossings; i++)
					{
						edge = (start+1+i - first_edge_on_component[start_component])%num_component_edges[start_component] + first_edge_on_component[start_component];
						
						/* check that this edge terminates at a crossing that remains on the start detour when the peer
						   detour has moved.  If the peer detour intersects the start detour, we may lose crossings from
						   the start detour.
						*/
						if (start_virtual_crossing_flag[term_crossing[edge]] == 1)
						{
							/* for knotoids, check whether the old start path edge lies within the shortcut, since in this case
							   the new crossing on the moved peer path will be virtual.
							*/
							bool start_path_edge_in_shortcut = false;
							if (pure_knotoid_code_data && edge >= head_semi_arc)
							{
								if (num_components == 1 || ( num_components > 1 && edge < first_edge_on_component[1]))
								{
									start_path_edge_in_shortcut = true;
									
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: knotoid: start path edge " << edge << " lies within the shortcut" << endl;
								}
							}
							
							int start_new_edge_label_column = term_crossing[edge];
							int row_zero_new_label = new_edge_labels[0][start_new_edge_label_column];
							int row_one_new_label = new_edge_labels[1][start_new_edge_label_column];
							
							int new_start_path_edge_label;
							if (edge%2 == 0)
								new_start_path_edge_label = row_zero_new_label;
							else
								new_start_path_edge_label = row_one_new_label;
							
							int start_new_even_edge = (row_zero_new_label %2? row_one_new_label: row_zero_new_label);
							
							int start_detour_crossing = start_new_even_edge/2;											
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   start path old edge label " << edge << ", start_new_edge_label_column " << start_new_edge_label_column << endl;
	debug << "remove_Reidemeister_II:   new_start_path_edge_label = " << new_start_path_edge_label << ", start_new_even_edge " << start_new_even_edge 
	      << ", start_detour_crossing " << start_detour_crossing << endl;
}

							int peer_new_edge_label_column;
							int peer_detour_crossing;
							int peer_new_even_edge;
							if (forwards_from_peer)
							{
								peer_new_edge_label_column = num_crossings+index;
								peer_new_even_edge = (new_edge_labels[0][peer_new_edge_label_column] % 2 ? new_edge_labels[1][peer_new_edge_label_column]: new_edge_labels[0][peer_new_edge_label_column]);								
								peer_detour_crossing = peer_new_even_edge/2;
								if (new_code_table[generic_code_data::table::TYPE][start_detour_crossing] == generic_code_data::TYPE1)
									new_code_table[generic_code_data::table::TYPE][peer_detour_crossing] = generic_code_data::TYPE2;
								else
									new_code_table[generic_code_data::table::TYPE][peer_detour_crossing] = generic_code_data::TYPE1;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   moving forwards from peer, peer_new_edge_label_column = " << peer_new_edge_label_column << ", peer_new_even_edge " << peer_new_even_edge 
		  << ", peer_detour_crossing = " << peer_detour_crossing << endl;
	debug << "remove_Reidemeister_II:   peer_detour_crossing type is " << (new_code_table[generic_code_data::table::TYPE][peer_detour_crossing] == generic_code_data::TYPE1? "TYPE1" : "TYPE2") << endl;
}
							}
							else
							{
								peer_new_edge_label_column = num_new_edge_label_columns-1-index;
								peer_new_even_edge = (new_edge_labels[0][peer_new_edge_label_column] % 2 ? new_edge_labels[1][peer_new_edge_label_column]: new_edge_labels[0][peer_new_edge_label_column]);								
								peer_detour_crossing = peer_new_even_edge/2;

								new_code_table[generic_code_data::table::TYPE][peer_detour_crossing] = new_code_table[generic_code_data::table::TYPE][start_detour_crossing];

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   moving backwards from peer, peer_new_edge_label_column = " << peer_new_edge_label_column << ", peer_new_even_edge " << peer_new_even_edge 
		  << ", peer_detour_crossing = " << peer_detour_crossing << endl;
	debug << "remove_Reidemeister_II:   peer_detour_crossing type is " << (new_code_table[generic_code_data::table::TYPE][peer_detour_crossing] == generic_code_data::TYPE1? "TYPE1" : "TYPE2") << endl;
}
							}
	
							/* Set the label of the peer_detour_crossing.  For knotoids, we have to determine whether the terminating edge on the 
							   peer path lies within the shortcut.  From the above diagrams it can be seen that the terminating edge on the peer 
							   path has the opposite polarity to the terminating edge on the start path (i.e new_start_path_edge_label), so we can determine which entry 
							   in column peer_new_edge_label_column of the new_edge_labels we need to consider.
							*/
							bool peer_path_edge_in_shortcut = false;
							int peer_path_terminating_edge;
							if (new_start_path_edge_label%2 ==1)
								peer_path_terminating_edge = peer_new_even_edge;
							else
								peer_path_terminating_edge = (new_edge_labels[0][peer_new_edge_label_column] % 2 ? new_edge_labels[0][peer_new_edge_label_column]: new_edge_labels[1][peer_new_edge_label_column]);								
								
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   peer_path_terminating_edge = " << peer_path_terminating_edge << endl;
								
							if (pure_knotoid_code_data && peer_path_terminating_edge >= new_head_semi_arc)
							{
								if (num_components == 1 || (num_components > 1 && peer_path_terminating_edge <= new_last_component_edge[0]))
								{
									peer_path_edge_in_shortcut = true;
									
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: knotoid: peer path edge " << peer_path_terminating_edge << " lies within the shortcut" << endl;
								}
							}
							
							if (!peer_path_edge_in_shortcut &&(new_code_table[generic_code_data::table::LABEL][start_detour_crossing] == generic_code_data::VIRTUAL || start_path_edge_in_shortcut))
							{
								new_code_table[generic_code_data::table::LABEL][peer_detour_crossing] = generic_code_data::VIRTUAL;
							}
							else if (peer_path_edge_in_shortcut)
							{
								if (peer_path_terminating_edge%2 == 0)
									new_code_table[generic_code_data::table::LABEL][peer_detour_crossing] = generic_code_data::NEGATIVE;
								else
									new_code_table[generic_code_data::table::LABEL][peer_detour_crossing] = generic_code_data::POSITIVE;
							}
							else
							{
								if (new_code_table[generic_code_data::table::LABEL][start_detour_crossing] == generic_code_data::POSITIVE)
									new_code_table[generic_code_data::table::LABEL][peer_detour_crossing] = generic_code_data::NEGATIVE;
								else
									new_code_table[generic_code_data::table::LABEL][peer_detour_crossing] = generic_code_data::POSITIVE;
							}
						
							index++;
							
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II:   peer_detour_crossing label is ";
	switch(new_code_table[generic_code_data::table::LABEL][peer_detour_crossing])
	{
		case generic_code_data::POSITIVE: debug << "POSITIVE"; break;
		case generic_code_data::NEGATIVE: debug << "NEGATIVE"; break;
		case generic_code_data::VIRTUAL: debug << "VIRTUAL"; break;
		default: debug << "ERROR!";
	}
	debug << endl;
}							
						}
						else
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   edge " << edge << " terminates at a virtual crossing that lies on the peer detour, so has been removed" << endl; 
						}
						
					}
						
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: new_code_table" << endl;
	debug << "remove_Reidemeister_II:  type: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::TYPE][i] << ' ';
	debug << endl;
	debug << "remove_Reidemeister_II:  odd peer: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::OPEER][i] << ' ';
	debug << endl;
	debug << "remove_Reidemeister_II:  even peer: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::EPEER][i] << ' ';
	debug << endl;
	debug << "remove_Reidemeister_II:  odd term: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::ODD_TERMINATING][i] << ' ';
	debug << endl;
	debug << "remove_Reidemeister_II:  even term: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::EVEN_TERMINATING][i] << ' ';
	debug << endl;
	debug << "remove_Reidemeister_II:  label: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::LABEL][i] << ' ';
	debug << endl;
	debug << "remove_Reidemeister_II:  component: ";
	for (int i=0; i< new_num_crossings; i++)
		debug << new_code_table[generic_code_data::table::COMPONENT][i] << ' ';
	debug << endl;
	
	if (track_zig_zag_counts)
	{
		debug << "remove_Reidemeister_II: new_zig_zag_count: " << endl;
		print(new_zig_zag_count,debug, 3,"remove_Reidemeister_II: ");	
	}
	
}
				
					/* copy the new code_table to code_data and set the num_crossings */
					code_data.code_table = new_code_table;
					code_data.num_crossings = new_num_crossings;				
					code_data.head = -1;  // for knotoids we'll set this correctly once the new peer code is read back in.

					if (track_zig_zag_counts)
					{
						bool zig_zags_present = false;
						for (int i=0; i<new_num_crossings; i++)
						{
							if (new_zig_zag_count[0][i] != 0 || new_zig_zag_count[1][i] !=0)
							{
								zig_zags_present = true;
								break;
							}
						}
						
						if (zig_zags_present)
							code_data.zig_zag_count = new_zig_zag_count;
						else
							code_data.zig_zag_count.clear();
					}

					for (int i=0; i< code_data.num_crossings; i++)
						code_data.num_component_edges[i] = new_last_component_edge[i] - new_first_component_edge[i]+1;
					
					/* write the modified code to a string */
					ostringstream oss;
					write_code_data (oss, code_data);
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: new peer code with head cleared = " << oss.str() << endl;
					
					/* read back the peer code into code_data */
					read_peer_code(code_data, oss.str());
				
					/* We can't set the head until we have the generic_code_data::table::EPEER row of the code table, which is completed by read_peer_code,
					   so have to set both the head and the immersion character manually here.
					
					   If new_head_semi_arc == 0, we have already cleared pure_knotoid_code_data 
					*/
					if (pure_knotoid_code_data)
					{
						if (new_head_semi_arc == 0)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: new_head_semi_arc is zero, already cleared knotoid status of code data" << endl;
						}
						else if (new_head_semi_arc % 2)
							code_data.head = code_data.code_table[generic_code_data::table::EPEER][(new_head_semi_arc-1)/2]/2;
						else
							code_data.head = new_head_semi_arc/2;

						code_data.immersion_character = generic_code_data::character::PURE_KNOTOID;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: new peer code head = " << code_data.head << endl;
							
					}
				
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: new code_data: " << endl;
	print_code_data (debug, code_data, "remove_Reidemeister_II:   ");
}		
				}
			}
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: no Reidemister II detour found" << endl;
		}
	}

	if (Reidemeister_II_found || floating_component != 0)
	{
		if (code_data.num_crossings > 0)
		{
			/* if we're dealing with a knotoid, we need to reset the shortcut crossings, do this by using
			   valid_knotoid_input, which will assign a new vector to code_data.shortcut crossing, thereby
			   deleting what is already there.
			   
			   We cannot rely on pure_knotoid_code_data to see whether we still have to check that we still have 
			   valid knotoid data, since remove_edge_flags_from_peer_code may have cleared the head if the
			   diagram reduces to a knot.
			*/
			if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && code_data.head != -1) 
			{
				if (!valid_knotoid_input(code_data))
				{
					cout << "\nError!  Function remove_Reidemeister_II has produced invalid knotoid code data" << endl;
					exit(0);
				}
			}
			
			/* Remove any Reidemeister I moves we may have created by the Reidemeister II move
			   and recurse in case we've got multiple "stacked" Reidemeister II moves, providing we've
			   still got at least three crossings.
			*/
			num_components_removed += remove_virtual_Reidemeister_I(code_data,component_flags);
			
			if (code_data.num_crossings > 1)
			{
			
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "remove_Reidemeister_II: after removing virtual Reidemeister_I moves peer code = ";
	write_code_data(debug,code_data);
	debug << endl;
	debug << "remove_Reidemeister_II: recursing..." << endl;
}
				num_components_removed += remove_Reidemeister_II(code_data,component_flags);
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II:   removing Reidemeister II configurations updates num_components_removed to " << num_components_removed << endl;
				
			}
		}
		else
		{
			num_components_removed = num_components;		
			for (unsigned int i=0; i< component_flags.size(); i++)
				component_flags[i] = 0;
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "remove_Reidemeister_II: no crossings remaining, clearing component_flags, num_components_removed = " << num_components_removed << endl;
		}
	}
	else
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "remove_Reidemeister_II: no Reidemeister II moves detected." << endl;	
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "remove_Reidemeister_II: returning num_components_removed = " << num_components_removed << endl;	
	
	return 	num_components_removed;	
}

/* align_label_numbers aligns the entries of edge_labels to ensure there is an odd and even terminating edge label at 
   each crossing.  The matrix edge_labels must have two rows and num_crossings columns.
   We take the front component from the component_tree and trace that component.  If we encounter a crossing that is 
   mis-aligned, we cycle the edge labels on the peer component, and push that component number onto the back of the 
   component tree.  In this manner we conduct a depth-first search of the component tree.
*/
void align_label_numbers (matrix<int>& edge_labels, int num_components,vector<int> new_first_component_edge,vector<int> new_last_component_edge)
{
	int num_crossings = edge_labels.numcols();
	vector<int> aligned_labels(num_components);
	bool start_component_found;
	
	/* at least one of the components has crossings but we may have removed floating components, so
	   we start the list with the first component that has a new_last_component_edge not equal to -1
	*/
	int start_component = 0;
	while (new_last_component_edge[start_component] == -1)
	{
		aligned_labels[start_component] = 1;
		start_component++;
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "align_label_numbers: first component with remaining crossings is " << start_component << endl;
	debug << "align_label_numbers: aligned_labels: ";
	for (int i=0; i< num_components; i++)
		debug << aligned_labels[i] << ' ';
	debug << endl;
}

	do
	{
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "align_label_numbers: start_component = " << start_component << endl;
	
		aligned_labels[start_component] = 1;

		list<int> component_tree;
		component_tree.push_back(start_component);
		list<int>::iterator tree_ptr = component_tree.begin();
		
		while (tree_ptr != component_tree.end())
		{
			int component = *tree_ptr;
			int first = new_first_component_edge[component];
			int last = new_last_component_edge[component];
							
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "align_label_numbers:   component " << component << ", first edge = " << first << ", last edge = " << last << endl;
	
			for (int edge=first; edge<=last; edge++)
			{
				int column;
				int peer;		
				
				for (int i=0; i< num_crossings; i++)
				{
					if (edge_labels[0][i] == edge)
					{
						peer = edge_labels[1][i];
						column = i;
						break;
					}
					else if (edge_labels[1][i] == edge)
					{
						peer = edge_labels[0][i];
						column = i;
						break;
					}
				}
				
				int peer_component = 0;
				for (int i=1; i< num_components; i++)
				{
					if (peer >= new_first_component_edge[i])
						peer_component++;
					else
						break;
				}

			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "align_label_numbers:     edge " << edge << ", peer = " << peer << " lies on component " << peer_component << endl;

				if (edge%2 == peer%2)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "align_label_numbers:     edge labels in column " << column << " has incompatible terminating edge labels " << endl;
	
					int peer_first_edge = new_first_component_edge[peer_component];				
					int peer_last_edge = new_last_component_edge[peer_component];
					int peer_num_edges = peer_last_edge - peer_first_edge + 1;
															
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "align_label_numbers:       cycle the " << peer_num_edges << " edge labels on peer_component " << peer_component << ", peer_first_edge = " << peer_first_edge << " peer_last_edge = " << peer_last_edge << endl;
	
					for (int j=0; j< num_crossings; j++)
					{
						if (edge_labels[0][j] <= peer_last_edge && edge_labels[0][j] >= peer_first_edge)
							edge_labels[0][j] = (edge_labels[0][j] - peer_first_edge +1)%peer_num_edges + peer_first_edge;
							
						if (edge_labels[1][j] <= peer_last_edge && edge_labels[1][j] >= peer_first_edge)
							edge_labels[1][j] = (edge_labels[1][j] - peer_first_edge +1)%peer_num_edges + peer_first_edge;
					}					

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "align_label_numbers:       renumbered new_edge_labels: " << endl;
	print(edge_labels, debug,3,"align_label_numbers:       ");
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

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "align_label_numbers: reached the end of the component_tree" << endl;
	debug << "align_label_numbers: aligned_labels: ";
	for (int i=0; i< num_components; i++)
		debug << aligned_labels[i] << ' ';
	debug << endl;
}
	
		/* look for another start_component */
		start_component_found = false;
		for (int i=0; i< num_components; i++)
		{
			if (aligned_labels[i] == 0)
			{
				start_component = i;
				start_component_found = true;
				break;
			}
		}
		
	} while (start_component_found);
}

void clear_component_flag(int component, vector<int>& component_flags)
{
	int count = -1;
	for (unsigned int i=0; i< component_flags.size(); i++)
	{
		if (component_flags[i] == 1)
		{
			count++;
			if (count == component)
			{
				component_flags[i] = 0;
				break;
			}
		}
	}
}

/* Reidemeister_III_present is based on the very clever enumeration of the eight possible R3 configurations that 
   appears in Jeremy Green's code.  Given a topmost arc of an R3 configuration going over consecutive crossings 
   a and b in the Gauss code, we seek a crossing, c, adjacent in the Gauss code to the two under-crossing terms 
   of a and b.  To reject situations where a virtual crossing prevents an R3 move from being performed, as in:

				               | a                     | b
				         -------------------------------------> 
				               |                       |
				         ------x----------             |
				               |         | c           |
				               ----------o--------------
				                         |
										 
   we shall require that a term relating to crossing c is not adjacent to both the under-crossing terms of a and b.
   
   The enumeration of the eight cases allows us to verify the combination of signs amongst the three crossings
   involved in the Reidemeister III move.
   The indices of the under-crossings of a and b are identified as under-a_index and under_b_index.  The indices and crossings
   to the left and right of these indices are identified as left_a_index, left_a_crossing, right_a_index, right_a_crossing and
   similarly for b.
   
   In the following the control variables, k, m2 and m3 have been chosen to correspond to those used by Jeremy.
   
   We then check for the following conditions
    - left_a_crossing = left_b_crossing    k=0  m2=1   case 1
	- left_a_crossing = right_b_crossing   k=0  m2=1   case 2
    - right_a_crossing = left_b_crossing   k=1  m2=-1  case 1
	- right_a_crossing = right_b_crossing  k=1  m2=-1  case 2

	k=0 m2=1, case1         k=0, m2=1, case 2       k=1, m2=-1, case 1      k=1, m2=-1, case 2
	   ^     ^                 ^                             ^
	  a|     |b               a|     |b               a|     |b               a|     |b
	------------->          ------------->          ------------->          ------------->
	 + |     | +             + |     | -             - |     | +             - |     | -
	    \   /                   \   /                   \   /                   \   /
        c / +                   c / -                   c / -                   c / +
		/   \                   /   \                   /   \                   /   \
	   |     |                 |     |                 |     |                 |     |
	                           v                             v                 v     v
							   
	Note that in the above diagrams the strand through crossing_a at the left_a_crossing or right_a_crossing 
	is the under-arc at crossing c.
	
	Note also that the sign of the crossings are as follows, when m3=1:
	
    case 1:  sign(a) = m3 * sign(c) and  sign(a) = m2 * sign(b)
    case 2: -sign(a) = m3 * sign(c) and -sign(a) = m2 * sign(b)
	   
   Next, we reverse crossings a and b and repeat the above checks with m3=-1;  This gives rise to the following diagrams
   
	k=2 m2=1, case1         k=2, m2=1, case 2       k=3, m2=-1, case 1      k=3, m2=-1, case 2
	   ^      ^                      ^                 ^      
	  b|     |a               b|     |a               b|     |a               b|     |a
	------------->          ------------->          ------------->          ------------->
	 + |     | +             - |     | +             + |     | -             - |     | -
	    \   /                   \   /                   \   /                   \   /
        c \ -                   c \ +                   c \ +                   c \ -
		/   \                   /   \                   /   \                   /   \
	   |     |                 |     |                 |     |                 |     |
	                                 v                 v                       v     v

   Note that in the above diagrams we still have that the strand through crossing_a at the left_a_crossing 
   or right_a_crossing is the under-arc at crossing c.
   
   The signs in these diagrams are then given by the same relations as above but with m3 = -1
   
   The fact that the strand through crossing_a is the under-arc at crossing_c is important, since non-Reidemeister III
   configurations do not have this property, even if all other conditions are matched.  An example of this is given by
   1 2 -1 3 -2 -3 / + + - , where the signs match for the over-arc at terms 1 and 2 checking for the right_a_index=3 
   but the term at index 3 is the over-arc not the under-arc.
   
   The function returns a set of integer vectors of length 6 <3 indices><3 crossing numbers>.  The indices are the indices 
   into the Gauss code (i.e the ortientation matrix or classical_gauss_data in gauss_orientation_data) of the first crossing 
   on each of the thre arcs involved in the configuration, consistent with the orientation determined by the code numbering.  
   The crossing numbers are the crossings of the over-arc followed by the dominated crossing. These crossing numbers are taken 
   from the classical_gauss_data and, as such, are numbered from 1.
   
   We return indices to avoid ambiguities such as in the case of -1 2 -3 4 1 3 -4 -2/+ - + +, where the over-arc 4 1 and 
   crossing 3 form a Reidemeister III configuration but it is not clear from the code which of the -2 or 2 crossing either 
   side of the under crossing -1 is part of the same arc.  The crossing numbers are for human convenience when comparing
   debug or programme output with diagrams.
   
   Note that, being indices, the return values are numbered from zero, so although the indices correspond to Gauss edges in a diagram, 
   the Gauss edges are generally numbered from 1.
   
*/
Reidemeister_III_return Reidemeister_III_present (generic_code_data code_data)
{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "Reidemeister_III_present: presented with code data ";
	write_code_data(debug,code_data);	
	debug << endl;
}

	Reidemeister_III_return _return;
	
	gauss_orientation_data gdata(code_data);
	
	
	list<vector<int> > Reidemeister_III_list;
	
	int num_terms = gdata.num_terms;
	int num_components = gdata.num_components;
	
	vector<int>& gauss_data = gdata.classical_gauss_data;
	vector<int>& crossing_sign = gdata.classical_crossing_sign;
	vector<int>& num_terms_in_component = gdata.num_terms_in_component;
	vector<int>& start_of_component = gdata.start_of_component;
	
	for (int c=0; c < num_components; c++)
	for (int i=0; i< num_terms_in_component[c]; i++)
	{
		int term = start_of_component[c]+i;
		int successor = start_of_component[c] + (i+1)%num_terms_in_component[c];
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Reidemeister_III_present: term = " << term << " successor = " << successor << endl;
		
		if (term != successor && gauss_data[term] > 0 && gauss_data[successor] > 0)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Reidemeister_III_present:   consecutive over crossings found at terms " << term << " and " << successor << endl;

			/* we have a pair of consecutive over-crossings, numbered from 1 */
			int crossing_a = gauss_data[term];
			int crossing_b = gauss_data[successor];
			int crossing_a_index = term;
			int crossing_b_index = successor;
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Reidemeister_III_present:   crossing_a = " << crossing_a << " crossing_b = " << crossing_b << endl;
			
			/* find the under arcs of crossing_a and crossing_b */
			int under_a_index;
			int under_b_index;
			for (int j=0; j < num_terms; j++)
			{
				if (gauss_data[j] == crossing_a * -1)
					under_a_index = j;					
				else if (gauss_data[j] == crossing_b * -1)
					under_b_index = j;					
			}
			
			int under_a_component=0;
			for (int j=1; j< num_components; j++)
			{
				if (under_a_index >= start_of_component[j])
					under_a_component++;
				else
					break;
			}

			int under_b_component=0;
			for (int j=1; j< num_components; j++)
			{
				if (under_b_index >= start_of_component[j])
					under_b_component++;
				else
					break;
			}
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Reidemeister_III_present:   under_a_component = " << under_a_component << " under_b_component = " << under_b_component << endl;
			
			
			int left_a_index = (under_a_index-start_of_component[under_a_component]-1+num_terms_in_component[under_a_component])%num_terms_in_component[under_a_component] + start_of_component[under_a_component];
			int right_a_index = (under_a_index-start_of_component[under_a_component]+1)%num_terms_in_component[under_a_component] + start_of_component[under_a_component];
			int left_b_index = (under_b_index-start_of_component[under_b_component]-1+num_terms_in_component[under_b_component])%num_terms_in_component[under_b_component] + start_of_component[under_b_component];
			int right_b_index = (under_b_index-start_of_component[under_b_component]+1)%num_terms_in_component[under_b_component] + start_of_component[under_b_component];
			int left_a_crossing = abs(gauss_data[left_a_index]);
			int right_a_crossing = abs(gauss_data[right_a_index]);
			int left_b_crossing = abs(gauss_data[left_b_index]);
			int right_b_crossing = abs(gauss_data[right_b_index]);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Reidemeister_III_present:   under_a_index = " << under_a_index << " under_b_index = " << under_b_index << endl;
	debug << "Reidemeister_III_present:   left_a_index = " << left_a_index << " left_a_crossing = " << left_a_crossing << endl;
	debug << "Reidemeister_III_present:   right_a_index = " << right_a_index << " right_a_crossing = " << right_a_crossing << endl;
	debug << "Reidemeister_III_present:   left_b_index = " << left_b_index << " left_b_crossing = " << left_b_crossing << endl;
	debug << "Reidemeister_III_present:   right_b_index = " << right_b_index << " right_b_crossing = " << right_b_crossing << endl;
}
			
			/* Start the enumeration described above */
			int m2 = 0;
			int m3 = 1;
			for (int k=0; k< 4; k++)
			{
				
				if (k==2)
				{
					/* invert m3 and swap the a and b references */
					m3 = -1;
					swap(crossing_a,crossing_b);
					swap(crossing_a_index,crossing_b_index);
					swap(under_a_index,under_b_index);
					swap(left_a_index,left_b_index);
					swap(right_a_index,right_b_index);
					swap(left_a_crossing,left_b_crossing);
					swap(right_a_crossing,right_b_crossing);
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Reidemeister_III_present:   swap crossings a and b " << endl;
				}
				
				int index_c;
				int crossing_c;
				if (k%2 == 0) // check for crossing_c as the left_a_crossing
				{
					index_c = left_a_index;
					crossing_c = left_a_crossing;
					m2 = 1;
				}
				else // check for crossing_c as the right_a_crossing
				{
					index_c = right_a_index;
					crossing_c = right_a_crossing;
					m2 = -1;
				}
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Reidemeister_III_present:   k = " << k << (k%2 == 0? " left": " right") << "_a check, index_c = " << index_c << " crossing_c = " 
	      << crossing_c << " m2 = " << m2 << " m3 = " << m3 << endl;
}
			
				/* check whether the crossings left and right of the under_a_index can also be found either left or right of the under_b_index 
				   We need to ensure there is more than one clasical crossing between under_a_index and under_b_index, since otherwise there would
				   be a virtual crossing preventing a Reidemeister III move.
				   
				   The condition gauss_data[index_c] < 0 ensures that the strand through crossing_a is the under-arc at crossing_c
				*/
//			if (t >= 0 && t != jju && t != j && t != jj) 
//			if (index_c != under_b_index && index_c != crossing_a_index && t != crossing_b_index)
				if (crossing_c != crossing_a && crossing_c != crossing_b && gauss_data[index_c] < 0)
				{					
					/* case 1 check the left_b_index */
					if (left_b_crossing == crossing_c && left_b_index != index_c && left_b_crossing != crossing_a && left_b_crossing != crossing_b)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Reidemeister_III_present:     case 1: left_b_crossing " << left_b_crossing << " == crossing_c " << crossing_c << endl;
						/* check the crossing signs: case 1:  sign(a) = m3 * sign(c) and  sign(a) = m2 * sign(b) */
						if (crossing_sign[crossing_a-1] == m3 * crossing_sign[crossing_c-1] && crossing_sign[crossing_a-1] == m2 * crossing_sign[crossing_b-1])
						{					
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Reidemeister_III_present:     found Reidemeister III " << crossing_a << ' ' << crossing_b << ' '<< crossing_c << endl;
	debug << "Reidemeister_III_present:     crossing_sign[crossing_a] = " << crossing_sign[crossing_a-1] << " m3 = " << m3 
	      << " crossing_sign[crossing_c] = " << crossing_sign[crossing_c-1] << " m2 = " << m2 << " crossing_sign[crossing_b] = " << crossing_sign[crossing_b-1] << endl;
}
							vector<int> R3_indices(6);
							switch(k)
							{
								case 0: R3_indices = {crossing_a_index,left_a_index,left_b_index, crossing_a, crossing_b, crossing_c};break;
								case 1: R3_indices = {crossing_a_index,under_a_index,left_b_index, crossing_a, crossing_b, crossing_c};break;
								case 2: R3_indices = {crossing_b_index,left_b_index,left_a_index, crossing_a, crossing_b, crossing_c};break;
								case 3: R3_indices = {crossing_b_index,left_b_index,under_a_index, crossing_a, crossing_b, crossing_c};break;
							}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Reidemeister_III_present:     Reidemeister III indices ";
	for (int i=0; i< 6; i++)
		debug << R3_indices[i] << ' ';
	debug << endl;
}

							Reidemeister_III_list.push_back(R3_indices);
							
							_return.Reidemeister_III_found = true;
						}
						else
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Reidemeister_III_present:     case 1 crossing signs do not match a == m3*c: " << (crossing_sign[crossing_a-1] == m3 * crossing_sign[crossing_c-1])
	      << " a = m2 * b: " << (crossing_sign[crossing_a-1] == m2 * crossing_sign[crossing_b-1]) << endl;
}
						}
					}

					/* case 2 check the right_b_index */
					if (right_b_crossing == crossing_c && right_b_index != index_c && right_b_crossing != crossing_a && right_b_crossing != crossing_b)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Reidemeister_III_present:     case 2: right_b_crossing " << right_b_crossing << " == crossing_c " << crossing_c << endl;
						/* check the crossing signs:  case 2: -sign(a) = m3 * sign(c) and -sign(a) = m2 * sign(b) */
						if (crossing_sign[crossing_a-1] * -1 == m3 * crossing_sign[crossing_c-1] && crossing_sign[crossing_a-1] * -1 == m2 * crossing_sign[crossing_b-1])
						{					
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Reidemeister_III_present:     found Reidemeister III " << crossing_a << ' ' << crossing_b << ' '<< crossing_c << endl;
	debug << "Reidemeister_III_present:     crossing_sign[crossing_a] = " << crossing_sign[crossing_a-1] << " m3 = " << m3 
	      << " crossing_sign[crossing_c] = " << crossing_sign[crossing_c-1] << " m2 = " << m2 << " crossing_sign[crossing_b] = " << crossing_sign[crossing_b-1] << endl;
}
							vector<int> R3_indices(6);
							switch(k)
							{
								case 0: R3_indices = {crossing_a_index,left_a_index,under_b_index, crossing_a, crossing_b, crossing_c};break;
								case 1: R3_indices = {crossing_a_index,under_a_index,under_b_index, crossing_a, crossing_b, crossing_c};break;
								case 2: R3_indices = {crossing_b_index,under_b_index,left_a_index, crossing_a, crossing_b, crossing_c};break;
								case 3: R3_indices = {crossing_b_index,under_b_index,under_a_index, crossing_a, crossing_b, crossing_c};break;
							}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Reidemeister_III_present:     Reidemeister III indices ";
	for (int i=0; i< 6; i++)
		debug << R3_indices[i] << ' ';
	debug << endl;
}
							Reidemeister_III_list.push_back(R3_indices);
							
							_return.Reidemeister_III_found = true;
						}
						else
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Reidemeister_III_present:     case 2 crossing signs do not match -a == m3*c: " << (crossing_sign[crossing_a-1] * -1 == m3 * crossing_sign[crossing_c-1])
	      << " -a = m2 * b: " << (crossing_sign[crossing_a-1] * -1 == m2 * crossing_sign[crossing_b-1]) << endl;
}
						}
					}
				}			
			}
		}
	}

	if (_return.Reidemeister_III_found)
	{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "Reidemeister_III_present: found " << Reidemeister_III_list.size() << " Reidemeister III configurations" << endl;
		matrix<int> R3_configuration(Reidemeister_III_list.size(),6);
		
		list<vector<int> >::iterator lptr = Reidemeister_III_list.begin();
		int index=0;
		
		while (lptr !=	Reidemeister_III_list.end())
		{
			R3_configuration[index][0] = (*lptr)[0];
			R3_configuration[index][1] = (*lptr)[1];
			R3_configuration[index][2] = (*lptr)[2];
			R3_configuration[index][3] = (*lptr)[3];
			R3_configuration[index][4] = (*lptr)[4];
			R3_configuration[index][5] = (*lptr)[5];
			index++;
			lptr++;
		}
		_return.R3_configuration = R3_configuration;
	}
	
	return _return;
}

/*  To perform the Reidemeister III move we number the two new crossings by interchanging the numbering of the 
    two over-arc crossings, as shown below (i.e interchange a and b) so the signs of the new a and b crossings 
    are the same as the signs of the old a and b.
    

	   ^     ^                 ^                             ^
	  a|     |b               a|     |b               a|     |b               a|     |b
	------------->          ------------->          ------------->          ------------->
	 + |     | +             + |     | -             - |     | +             - |     | -
	    \   /                   \   /                   \   /                   \   /
        c / +                   c / -                   c / -                   c / +
		/   \                   /   \                   /   \                   /   \
	   |     |                 |     |                 |     |                 |     |
	                           v                             v                 v     v

	   ^     ^                 ^                             ^
	   |     |                 |     |                 |     |                 |     | 
	    \   /                   \   /                   \   /                   \   /
        c / +                   c / -                   c / -                   c / +
		/   \                   /   \                   /   \                   /   \
	  b|     |a               b|     |a               b|     |a               b|     |a
	------------->          ------------->          ------------->          ------------->
	 + |     | +             - |     | +             + |     | -             - |     | -
	                           v                             v                 v     v
	                           
	we then interchage the under-arc terms of the old a annd b with the adjacent term corresponding to crossing c
	
	The over_index indicates the first term involved in the over-arc, the a_index the first term on the under_arc 
	through a and the b_index the first term on the under-arc through b.
*/
generic_code_data Reidemeister_III_move (generic_code_data code_data, int over_index, int a_index, int b_index)
{
	gauss_orientation_data gdata(code_data);

	int num_terms = gdata.num_terms;
	int num_components = gdata.num_components;
	
	vector<int>& gauss_data = gdata.classical_gauss_data;
	vector<int>& crossing_sign = gdata.classical_crossing_sign;
	vector<int>& num_terms_in_component = gdata.num_terms_in_component;
	vector<int>& start_of_component = gdata.start_of_component;

	int over_index_component=0;
	for (int j=1; j< num_components; j++)
	{
		if (over_index >= start_of_component[j])
			over_index_component++;
		else
			break;
	}

	int under_a_component=0;
	for (int j=1; j< num_components; j++)
	{
		if (a_index >= start_of_component[j])
			under_a_component++;
		else
			break;
	}

	int under_b_component=0;
	for (int j=1; j< num_components; j++)
	{
		if (b_index >= start_of_component[j])
			under_b_component++;
		else
			break;
	}

	int over_index_successor = (over_index-start_of_component[over_index_component]+1)%num_terms_in_component[over_index_component] + start_of_component[over_index_component];
	int a_index_successor = (a_index-start_of_component[under_a_component]+1)%num_terms_in_component[under_a_component] + start_of_component[under_a_component];
	int b_index_successor = (b_index-start_of_component[under_b_component]+1)%num_terms_in_component[under_b_component] + start_of_component[under_b_component];
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Reidemeister_III_move: performing move at indices " << over_index << ' ' << a_index << ' ' << b_index << endl;
	debug << "Reidemeister_III_move: components " << over_index_component << ' ' << under_a_component << " and " << under_b_component << endl;
	debug << "Reidemeister_III_move: arcs " << gauss_data[over_index] << ' ' << gauss_data[over_index_successor] << ", " 
	                                       << gauss_data[a_index] << ' ' << gauss_data[a_index_successor] << ", " 
	                                       << gauss_data[b_index] << ' ' << gauss_data[b_index_successor] << endl;
	debug << "Reidemeister_III_move: initial gauss_data: ";
	for (int i=0; i< num_terms; i++)
		debug << gauss_data[i] << ' ';
	debug << endl;
}

	vector<int> new_gauss_data = gauss_data;
	
	new_gauss_data[over_index] = gauss_data[over_index_successor];
	new_gauss_data[over_index_successor] = gauss_data[over_index];
	new_gauss_data[a_index] = gauss_data[a_index_successor];
	new_gauss_data[a_index_successor] = gauss_data[a_index];
	new_gauss_data[b_index] = gauss_data[b_index_successor];
	new_gauss_data[b_index_successor] = gauss_data[b_index];
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Reidemeister_III_move:     new gauss_data: ";
	for (int i=0; i< num_terms; i++)
		debug << new_gauss_data[i] << ' ';
	debug << endl;
}

	ostringstream oss;
	int index = 0;
	
	for (int i=0; i< num_components; i++)
	{
		for (int j=0; j<num_terms_in_component[i]; j++)
		{
			oss << new_gauss_data[index++];
			if (j < num_terms_in_component[i] -1)
				oss << ' ';
		}
		
		if (i < num_components-1)
			oss << ',';
	}
	
	oss << '/';
	
	for (int i=0; i< num_terms/2; i++)
	{
		if (crossing_sign[i] > 0)
			oss << "+ ";
		else
			oss << "- ";
	}	
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "Reidemeister_III_move: oss.str() =  " << oss.str() << endl;

	generic_code_data updated_gauss_code_data;
	read_gauss_code(updated_gauss_code_data,oss.str());
	
	if (code_data.type == generic_code_data::peer_code)
	{
		generic_code_data updated_peer_code_data;
		gauss_to_peer_code(updated_gauss_code_data, updated_peer_code_data);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Reidemeister_III_move: returning updated_peer_code_data: ";
	write_code_data(debug, updated_peer_code_data);
	debug << endl;
	print_code_data(debug,updated_peer_code_data,"Reidemeister_III_move: ");
}	

		return updated_peer_code_data;
	}
	else
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "Reidemeister_III_move: returning updated_gauss_code_data: ";
	write_code_data(debug, updated_gauss_code_data);
	debug << endl;
	print_code_data(debug,updated_gauss_code_data,"Reidemeister_III_move: ");
}	
		return updated_gauss_code_data;
	}
}
