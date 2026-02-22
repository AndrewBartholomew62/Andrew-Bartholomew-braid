/**************************************************************************
string convert_gauss_code(string OU_gauss_code);
void renumber_peer_code(generic_code_data& code_data, vector<int> shift)
int remove_peer_code_component(generic_code_data& code_data, int component, vector<int>& component_flags)
int remove_virtual_components(generic_code_data& code_data, vector<int>& component_flags)
generic_code_data partition_peer_code(generic_code_data& code_data, vector<int>& component_flags)
void cycle_gauss_code_labels(generic_code_data& gauss_code_data, int edge)
void cycle_gauss_code_labels(matrix<int>& crossing_peers, vector<int>& num_component_edges, vector<int>& first_edge_on_component, int num_crossings,int num_components, int edge)
list<int> vertex_span (generic_code_data& peer_code_data, int initial_crossing, vector<int>* exclude)
bool calculate_turning_cycles(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles)
bool realizable_code_data(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles)
bool valid_knotoid_input(generic_code_data& code_data)
string over_preferred_gauss_code(generic_code_data& code_data, bool unoriented)
vector<int> classical_gauss_data(generic_code_data& code_data)
int amalgamate_zig_zag_counts(int a, int b)
void trace_component(generic_code_data& code_data, int start, int end, int component, int base_component, vector<bool>& visited_component, vector<bool>& visited_crossing,
bool check_2_separating_edges(generic_code_data& code_data,int cycle_a,int cycle_b,int edge_1,int edge_2)
bool check_2_separating_crossings(generic_code_data& code_data,int cycle_a,int cycle_b,int crossing_1,int crossing_2)
bool three_connected(generic_code_data& code_data, bool flat_crossings)
bool valid_multi_linkoid_input(generic_code_data& code_data)
int  code_data_writhe(generic_code_data& code_data)

**************************************************************************/
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <valarray>
#include <iomanip>
#include <list>

using namespace std;

namespace util
{
	string itos(long n);
}

extern ofstream     debug;
extern ofstream     output;
extern ifstream     input;

#include <util.h>
#include <matrix.h>
#include <generic-code.h>
#include <debug-control.h>
#include <gauss-orientation.h>
#include <reidemeister.h>


/* convert_gauss_code converts from OU to standard Gauss code format i.e from
   O1-O2+U1-O3+O4+U2+U3+U4+ to 1 2 -1 3 4 -2 -3 -4/- + + +
*/
string convert_gauss_code(string OU_gauss_code)
{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "convert_gauss_code: crossing " << "cptr = " << *cptr << " over_arc = " << (over_arc?"true":"false");			
		
		cptr++;
		
		int crossing_num;
		get_number(crossing_num, cptr);
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", crossing_num = " << crossing_num;			
		
		while (isdigit(*cptr))
			cptr++;
		
		if (*cptr == '-')
			sign[crossing_num-1] = -1;
		else
			sign[crossing_num-1] = 1;

		cptr++;
		
if (debug_control::DEBUG >= debug_control::DETAIL)
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
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "convert_gauss_code: returning gauss_code: " << oss.str() << endl;			
	
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
		if (code_data.code_table[generic_code_data::table::LABEL][code_data.head] == generic_code_data::POSITIVE)
			head_semi_arc = code_data.code_table[generic_code_data::table::OPEER][code_data.head];
		else if (code_data.code_table[generic_code_data::table::LABEL][code_data.head] == generic_code_data::NEGATIVE)
			head_semi_arc = 2*code_data.head;
		else
		{
			cout << "\nError! Function renumber_peer_code presented with a knotoid code whose first shortcut crossing is virtual." << endl;
			exit(0);
		}
		
		if (shift[0] != 0 && shift[0] != -head_semi_arc)	
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
		debug << "renumber_peer_code: asked to renumber a knotoid with a shift vector that changes the shortcut, doing nothing" << endl;
		
		   return;
		}
    }
	
	int num_crossings = code_data.num_crossings;
	int num_edges = 2* num_crossings;
	int num_components = code_data.num_components;
	matrix<int>& code_table = code_data.code_table;
	

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
//	debug << "renumber_peer_code: given code data:" << endl;
//	print_code_data(debug,code_data,"renumber_peer_code: ");	
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
		int component = code_table[generic_code_data::table::COMPONENT][i];
		initial_edge_map[2*i] = (2*i - first_edge_on_component[component] - abs(shift[component]) + num_component_edges[component])%
		                num_component_edges[component] + first_edge_on_component[component];

		initial_edge_map[2*i+1] = (initial_edge_map[2*i] + 1 - first_edge_on_component[component])%
		                num_component_edges[component] + first_edge_on_component[component];		   
	}
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
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
				int position=-1;
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
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
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
		int component = code_table[generic_code_data::table::COMPONENT][i];
		int component_edges = num_component_edges[component];
		int first_edge = first_edge_on_component[component];
		int peer_component = code_table[generic_code_data::table::COMPONENT][(code_table[generic_code_data::table::OPEER][i]-1)/2];
		int peer_component_edges = num_component_edges[peer_component];
		int peer_first_edge = first_edge_on_component[peer_component];
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code: old crossing " << i << ", component " << component 
	      << ", peer_component " << peer_component << endl; 
}

		if (edge_map[2*i]%2)
		{

			if (shift[component] < 0 && shift[peer_component] < 0)
			{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " now odd" << endl;
	debug << "renumber_peer_code:   both components at old crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   type and label unchanged for new crossing" << endl;
}

				even_edge = (edge_map[2*i]-first_edge-1+component_edges)%component_edges + first_edge;
				odd_edge = (edge_map[code_table[generic_code_data::table::OPEER][i]]-peer_first_edge-1+peer_component_edges)%peer_component_edges + peer_first_edge;

				/* the new naming edge is on the same strand as the old naming edge */
				new_code_table[generic_code_data::table::LABEL][even_edge/2] = code_table[generic_code_data::table::LABEL][i];

				/* since the old naming edge is now odd, reversing both strands leaves the crossing type unchanged */
				new_code_table[generic_code_data::table::TYPE][even_edge/2] = code_table[generic_code_data::table::TYPE][i];


			}
			else if (shift[component] < 0)
			{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " now odd" << endl;
	debug << "renumber_peer_code:   naming edge component at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   type reversed but label unchanged for new crossing" << endl;
}

				even_edge = (edge_map[2*i]-first_edge-1+component_edges)%component_edges + first_edge;
				odd_edge = edge_map[code_table[generic_code_data::table::OPEER][i]];

				/* the new naming edge is on the same strand as the old naming edge */
				new_code_table[generic_code_data::table::LABEL][even_edge/2] = code_table[generic_code_data::table::LABEL][i];

				/* since the old naming edge is now odd, reversing the old naming strand reverses crossing type */
				if (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1)
					new_code_table[generic_code_data::table::TYPE][even_edge/2] = generic_code_data::TYPE2;
				else
					new_code_table[generic_code_data::table::TYPE][even_edge/2] = generic_code_data::TYPE1;
			}
			else if (shift[peer_component] < 0)
			{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " now odd" << endl;
	debug << "renumber_peer_code:   peer to naming edge component at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   type unchanged but label reversed for new crossing" << endl;
}

				odd_edge = edge_map[2*i];
				even_edge = (edge_map[code_table[generic_code_data::table::OPEER][i]]-peer_first_edge-1+peer_component_edges)%peer_component_edges + peer_first_edge;

				/* the new naming edge is on the opposite strand to the old naming edge */
				if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE)
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::NEGATIVE;
				else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::POSITIVE;
				else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::FLAT)
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::FLAT;
				else
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::VIRTUAL;

				/* since the old naming edge is now odd, reversing the peer of the old naming 
				   strand leaves the crossing type unchanged*/
				new_code_table[generic_code_data::table::TYPE][even_edge/2] = code_table[generic_code_data::table::TYPE][i];
			}
			else
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " now odd" << endl;
	debug << "renumber_peer_code:   neither component at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   both type and label reversed for new crossing" << endl;
}

				odd_edge = edge_map[2*i];
				even_edge = edge_map[code_table[generic_code_data::table::OPEER][i]];

				/* the new naming edge is on the opposite strand to the old naming edge */
				if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE)
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::NEGATIVE;
				else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::POSITIVE;
				else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::FLAT)
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::FLAT;
				else
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::VIRTUAL;

				/* since the old naming edge is now odd the crossing type is reversed */
				if (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1)
					new_code_table[generic_code_data::table::TYPE][even_edge/2] = generic_code_data::TYPE2;
				else
					new_code_table[generic_code_data::table::TYPE][even_edge/2] = generic_code_data::TYPE1;
			}			
		}
		else
		{
			if (shift[component] < 0 && shift[peer_component] < 0)
			{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code: current code data:" << endl;
	print_code_data(debug,code_data,"renumber_peer_code: ");	
	debug << "renumber_peer_code: current code_table:" << endl;
	print(code_table,debug, 3, "renumber_peer_code: ");	
	debug << "renumber_peer_code: code_table[generic_code_data::table::LABEL][" << i << "] = " << code_table[generic_code_data::table::LABEL][i] << endl;
	
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " still even" << endl;
	debug << "renumber_peer_code:   both components at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   both type and label reversed for new crossing" << endl;
}

				odd_edge = (edge_map[2*i]-first_edge-1+component_edges)%component_edges + first_edge;
				even_edge = (edge_map[code_table[generic_code_data::table::OPEER][i]]-peer_first_edge-1+peer_component_edges)%peer_component_edges + peer_first_edge;

				/* the new naming edge is on the opposite strand to the old naming edge */
				if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE)
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::NEGATIVE;
				else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::POSITIVE;
				else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::FLAT)
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::FLAT;
				else
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::VIRTUAL;

				/* since the old naming edge is still even, reversing both strands reverses the crossing type */
				if (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1)
					new_code_table[generic_code_data::table::TYPE][even_edge/2] = generic_code_data::TYPE2;
				else
					new_code_table[generic_code_data::table::TYPE][even_edge/2] = generic_code_data::TYPE1;
			}
			else if (shift[component] < 0)
			{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " still even" << endl;
	debug << "renumber_peer_code:   naming edge component at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   type unchanged but label reversed for new crossing" << endl;
}

				odd_edge = (edge_map[2*i]-first_edge-1+component_edges)%component_edges + first_edge;
				even_edge = edge_map[code_table[generic_code_data::table::OPEER][i]];

				/* the new naming edge is on the opposite strand to the old naming edge */
				if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE)
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::NEGATIVE;
				else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::POSITIVE;
				else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::FLAT)
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::FLAT;
				else
					new_code_table[generic_code_data::table::LABEL][even_edge/2] = generic_code_data::VIRTUAL;

				/* since the old naming edge is still even, reversing the old naming 
				   strand leaves the crossing type unchanged */
				new_code_table[generic_code_data::table::TYPE][even_edge/2] = code_table[generic_code_data::table::TYPE][i];
			}
			else if (shift[peer_component] < 0)
			{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " still even" << endl;
	debug << "renumber_peer_code:   peer to naming edge component at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   type reversed but label unchanged for new crossing" << endl;
}

				even_edge = edge_map[2*i];
				odd_edge = (edge_map[code_table[generic_code_data::table::OPEER][i]]-peer_first_edge-1+peer_component_edges)%peer_component_edges + peer_first_edge;

				/* the new naming edge is on the same strand as the old naming edge */
				new_code_table[generic_code_data::table::LABEL][even_edge/2] = code_table[generic_code_data::table::LABEL][i];

				/* since the old naming edge is still even, reversing the peer of the old naming 
				   strand reverses the crossing type */
				if (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1)
					new_code_table[generic_code_data::table::TYPE][even_edge/2] = generic_code_data::TYPE2;
				else
					new_code_table[generic_code_data::table::TYPE][even_edge/2] = generic_code_data::TYPE1;
			}
			else
			{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " still even" << endl;
	debug << "renumber_peer_code:   neither component at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   both type and label unchanged for new crossing" << endl;
}

				even_edge = edge_map[2*i];
				odd_edge = edge_map[code_table[generic_code_data::table::OPEER][i]];

				/* the new naming edge is on the same strand as the old naming edge */
				new_code_table[generic_code_data::table::LABEL][even_edge/2] = code_table[generic_code_data::table::LABEL][i];

				/* since the old naming edge is still even the crossing type is unchanged */
				new_code_table[generic_code_data::table::TYPE][even_edge/2] = code_table[generic_code_data::table::TYPE][i];
			}			
		}
		
		new_code_table[generic_code_data::table::OPEER][even_edge/2] = odd_edge;
		new_code_table[generic_code_data::table::EPEER][(odd_edge-1)/2] = even_edge;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code: renumbered as odd_edge = " << odd_edge << ", even_edge = " << even_edge << endl;
	debug << "renumber_peer_code: new type =  ";
	if (new_code_table[generic_code_data::table::TYPE][even_edge/2] == generic_code_data::TYPE1)
		debug << "generic_code_data::TYPE1";
	else
		debug << "generic_code_data::TYPE2";
	
	debug << ", new label = ";
	if (new_code_table[generic_code_data::table::LABEL][even_edge/2] == generic_code_data::POSITIVE)
		debug << "generic_code_data::POSITIVE";
	else if (new_code_table[generic_code_data::table::LABEL][even_edge/2] == generic_code_data::NEGATIVE)
		debug << "generic_code_data::NEGATIVE";
	else if (new_code_table[generic_code_data::table::LABEL][even_edge/2] == generic_code_data::VIRTUAL)
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
		term_crossing[new_code_table[generic_code_data::table::OPEER][i]] = i;
		
		/* we need to identify the edge following the naming edge's peer */
		int component = new_code_table[generic_code_data::table::COMPONENT][(new_code_table[generic_code_data::table::OPEER][i]-1)/2];
		int peer_successor = (new_code_table[generic_code_data::table::OPEER][i]+1 - first_edge_on_component[component])%
		                     num_component_edges[component] + first_edge_on_component[component];

		orig_crossing[peer_successor] = i;
		
		new_code_table[generic_code_data::table::EVEN_TERMINATING][i] = 2*i;
		new_code_table[generic_code_data::table::ODD_ORIGINATING][i] = 2*i+1;
		new_code_table[generic_code_data::table::ODD_TERMINATING][i] = new_code_table[generic_code_data::table::OPEER][i];
		new_code_table[generic_code_data::table::EVEN_ORIGINATING][i] = peer_successor;
	}
	
	code_data.code_table = new_code_table;
	code_data.term_crossing = term_crossing;
	code_data.orig_crossing = orig_crossing;

	/* re-set the head, if the code data is for a knotoid.  Note that if the head_semi_arc is even then 
	   the head is already set correctly, even if the orientation of the first componenet is reversed
   */
	if (code_data.head !=-1 && head_semi_arc %2)
    {	
		code_data.head = code_data.code_table[generic_code_data::table::EPEER][(head_semi_arc-1)/2]/2;
    }

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code: renumbered code data: ";
	write_peer_code(debug,code_data);
	debug << endl;	
	print_code_data(debug,code_data,"renumber_peer_code: ");	
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
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "remove_peer_code_component: removing component " << component << endl;

	int num_components_removed = 0;  //the return variable
	
	if (code_data.num_components > 1)
	{
		vector<int> edge_flag(2*code_data.num_crossings);
		
		for (int i= code_data.first_edge_on_component[component]; i < code_data.first_edge_on_component[component]+code_data.num_component_edges[component]; i++ )
		{
			int crossing = code_data.term_crossing[i];
			edge_flag[code_data.code_table[generic_code_data::table::ODD_TERMINATING][crossing]] = 1;
			edge_flag[code_data.code_table[generic_code_data::table::EVEN_TERMINATING][crossing]] = 1;
		}
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
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

/* remove_virtual_components removes any component in the peer code that contains only virtual crossings 
   the function returns the number of crossings it has removed.
*/
int remove_virtual_components(generic_code_data& code_data, vector<int>& component_flags)
{
	int num_virtual_components = 0;
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
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
			if (code_data.code_table[generic_code_data::table::LABEL][code_data.term_crossing[j]] != generic_code_data::VIRTUAL)
			{
				virtual_component = false;
				break;
			}
		}
				
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "remove_virtual_components: component " << i << (virtual_component? " is " : " is not ") << "a virtual component" << endl;
	
		if (virtual_component && code_data.num_components > 1)
		{
			num_virtual_components += remove_peer_code_component(code_data, i,component_flags);			
		}
		else if (virtual_component)
		{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "remove_virtual_components: no other components, setting num_crossings = 0" << endl;
	
			num_virtual_components = 1;
			code_data.num_crossings = 0;
			for (unsigned int i=0; i< component_flags.size(); i++)
				component_flags[i] = 0;			
		}
	}
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "remove_virtual_components: removed " << num_virtual_components << " components" << endl;
	
	return num_virtual_components;
}

/* partition_peer_code partitions code_data into the first connected component and the subsequent component(s).
   It returns generic_code_data describing the subsequent components(s) and adjusts code_data to reflect only the
   first connected component.
      
   If the code_data character is not CLOSED, code_data.head may not be -1, in which case the segment component must 
   be part of the first connected component, since the leg semi_arc is numbered zero.  That means the head and 
   immersion character is set correctly in code_data and that the default initialization of the head and immersion 
   character in subsequent_code_data is also correct.
   
   The vector component_flags is used by the calling code to track the unicursal components included in the first 
   connected component, partition_peer_code will clear the i-th set component_flag (which may not be at index i)
   if the i-th unicursal component of code_data does not belong to the first connected component.
   
   Both peer codes and Gauss codes are supported.
*/					
generic_code_data partition_peer_code(generic_code_data& code_data, vector<int>& component_flags)
{
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "partition_peer_code: presented with code data ";
	write_code_data(debug,code_data);	
	debug << endl;
	print_code_data(debug,code_data);
	debug << endl;
	debug << "partition_peer_code: component_flags: ";
	for (unsigned int j=0; j< component_flags.size(); j++)
		debug << component_flags[j] << ' ';
	debug << endl;
}
	generic_code_data first_code_data;
	generic_code_data subsequent_code_data;

	bool pure_knotoid_code_data = false;
	int head_semi_arc = -1;
	if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && code_data.shortcut_crossing.size())
	{
		pure_knotoid_code_data = true;
		
		if (code_data.code_table[generic_code_data::table::LABEL][code_data.head] == generic_code_data::POSITIVE)
			head_semi_arc = code_data.code_table[generic_code_data::table::OPEER][code_data.head];
		else if (code_data.code_table[generic_code_data::table::LABEL][code_data.head] == generic_code_data::NEGATIVE)
			head_semi_arc = 2*code_data.head;
		else
		{
			cout << "\nError! Function partition_peer_code presented with a knotoid code whose first shortcut crossing is virtual." << endl;
			exit(0);
		}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "partition_peer_code: head_semi_arc = " << head_semi_arc << endl;			
	}

	bool track_zig_zag_counts;
	if (code_data.zig_zag_count.numrows() !=0)
		track_zig_zag_counts = true;
	else
		track_zig_zag_counts = false;

	list<int> span_list = vertex_span (code_data, 0);
	
	int num_f_crossings = span_list.size();
	
if (debug_control::DEBUG >= debug_control::BASIC)
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
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "partition_peer_code: code_data not connected, first connected component contains " << num_f_crossings << " crossings" << endl;

		/* code_data is not connected and num_f_crossings identifies the number of crossings belonging to the first connected component, 
		   which must comprise all of the crossings on the unicursal components making up the first connected component.
		*/
		vector<int> f_components(code_data.num_components);
		list<int>::iterator lptr=span_list.begin();
		while (lptr != span_list.end())
		{
			f_components[code_data.code_table[generic_code_data::table::COMPONENT][*lptr]] = 1;
			lptr++;
		}
		
		int num_f_components=0;
		for (int i=0; i<code_data.num_components; i++)
			num_f_components += f_components[i];
		
if (debug_control::DEBUG >= debug_control::BASIC)
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
		   
		   We need the head, num_crossings and the generic_code_data::table::OPEER, generic_code_data::table::TYPE, generic_code_data::table::LABEL and generic_code_data::table::COMPONENT rows of the code_table in the generic_code_data 
		   structure to do this.		
		*/				
		int num_s_crossings = code_data.num_crossings - num_f_crossings;
		int num_s_components = code_data.num_components - num_f_components;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "partition_peer_code: num_f_components = " << num_f_components << " num_f_crossings = " << num_f_crossings << " num_s_components = " << num_s_components << " num_s_crossings = " << num_s_crossings << endl;
		
		for (int i=0; i< 2; i++)
		{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "partition_peer_code: determine the peer code of partition " << i << endl;
	
			int num_partition_crossings = (i==0? num_f_crossings: num_s_crossings);
			int num_partition_components = (i==0? num_f_components: num_s_components);
			matrix<int> partition_code_table(generic_code_data::table::CODE_TABLE_SIZE,num_partition_crossings,-1);
			matrix<int> partition_zig_zag_count(2,num_partition_crossings);
			vector<int> first_edge_on_partition_cpt(num_partition_components); // used for writing Gauss codes

			
			/* write new edge labels into the even (row 0)  and odd (row 1) rows of new_partition_labels, 
			   in the location corresponding to the old edge labels in the generic_code_data::table::EVEN_TERMINATING and generic_code_data::table::ODD_TERMINATING
			   rows of code_data.code_table.  In row 2 of new_partition_labels we record the component number
			   of the unicursal component in the partition.
			   
			   If we are working with a Gauss code, the third row of new_partition_labels records the crossing
			   number as determined by the new edge labels.
			*/
			matrix<int> new_partition_labels(4,code_data.num_crossings,-1);
			int edge = 0;
			int partition_crossing = 0;
			int component = -1;
			
			for (int j=0; j< code_data.num_components; j++)
			{
				/* f_components is a vector of flags, with f_components[j] == 1 iff the jth component is part of the first partition,
				   therefore for the first partition we want flags with value 1 and for the subsequent partition, flags with value 0
				*/
				if (f_components[j] == i)
				{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "partition_peer_code:   unicursal component " << j << " is not part of partition " << i << endl;
					
					continue;					
				}
				
				component++;
				
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "partition_peer_code:   unicursal component " << j << " is part of partition " << i << endl;
					
				first_edge_on_partition_cpt[component] = edge;
				
				int num_component_edges = code_data.num_component_edges[j];
				int first_edge = code_data.first_edge_on_component[j];
				
				for (int k=0; k < num_component_edges; k++)
				{
					int old_edge = first_edge + k;
					int row;
					int column;
										
					if (code_data.type == generic_code_data::gauss_code)
					{
						column = code_data.term_crossing[old_edge];
						if (code_data.code_table[generic_code_data::table::EVEN_TERMINATING][column] == old_edge)
						{
							row = 0;
						}
						else
						{
							row = 1;
							new_partition_labels[2][column] = component;
							new_partition_labels[3][column] = partition_crossing++;
						}
							
					}
					else
					{
						if (old_edge % 2 == 0)
						{
							row = 0; // even
							column = old_edge/2;
							
							new_partition_labels[2][column] = component;
						}
						else
						{
							row = 1; // odd
							column = code_data.code_table[generic_code_data::table::EPEER][(old_edge-1)/2]/2;
						}
					}
					
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "partition_peer_code:     old_edge " << old_edge << " is (row) " << (row == 0? "generic_code_data::table::EVEN_TERMINATING": "generic_code_data::table::ODD_TERMINATING") 
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
			}

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "partition_peer_code: first_edge_on_partition_cpt " ;
	for (int j=0; j< num_partition_components; j++)
		debug << first_edge_on_partition_cpt[j] << ' ';
	debug << endl;
}					
			
			/* write the new edge labels into the correct location in the generic_code_data::table::EPEER and generic_code_data::table::OPEER rows of the partition code table, 
			   and set the generic_code_data::table::TYPE generic_code_data::table::LABEL and generic_code_data::table::COMPONENT rows 
			*/
			for (int j=0; j < code_data.num_crossings; j++)
			{
				if (new_partition_labels[0][j] != -1)
				{
					int crossing;
					
					if (code_data.type == generic_code_data::gauss_code)
						crossing = new_partition_labels[3][j];
					else
						crossing = new_partition_labels[0][j]/2;
						
					partition_code_table[generic_code_data::table::EPEER][crossing] = new_partition_labels[0][j];
					partition_code_table[generic_code_data::table::OPEER][crossing] = new_partition_labels[1][j];
					partition_code_table[generic_code_data::table::COMPONENT][crossing] = new_partition_labels[2][j];
					partition_code_table[generic_code_data::table::TYPE][crossing] = code_data.code_table[generic_code_data::table::TYPE][j];
					partition_code_table[generic_code_data::table::LABEL][crossing] = code_data.code_table[generic_code_data::table::LABEL][j];
					
					if (track_zig_zag_counts)
					{
						partition_zig_zag_count[0][crossing] = code_data.zig_zag_count[0][j];
						partition_zig_zag_count[1][crossing] = code_data.zig_zag_count[1][j];
					}
				}
			}
				
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "partition_peer_code: new_partition_labels even terminating: ";
	for (int k=0; k< code_data.num_crossings; k++)
		debug << setw(3) << new_partition_labels[0][k];
	debug << "\npartition_peer_code: new_partition_labels odd terminating:  ";
	for (int k=0; k< code_data.num_crossings; k++)
		debug << setw(3) << new_partition_labels[1][k];
	debug << "\npartition_peer_code: new_partition_labels component:  ";
	for (int k=0; k< code_data.num_crossings; k++)
		debug << setw(3) << new_partition_labels[2][k];
	if (code_data.type == generic_code_data::gauss_code)
	{
		debug << "\npartition_peer_code: new_partition_labels partition crossing:  ";
		for (int k=0; k< code_data.num_crossings; k++)
			debug << setw(3) << new_partition_labels[3][k];
	}
	debug << endl;
	debug << "partition_peer_code: partition_code_table" << endl;
	debug << "partition_peer_code:  odd peer: ";
	for (int j=0; j< num_partition_crossings; j++)
		debug << partition_code_table[generic_code_data::table::OPEER][j] << ' ';
	debug << endl;
	debug << "partition_peer_code:  component: ";
	for (int j=0; j< num_partition_crossings; j++)
		debug << partition_code_table[generic_code_data::table::COMPONENT][j] << ' ';
	debug << endl;
	debug << "partition_peer_code:  type: ";
	for (int j=0; j< num_partition_crossings; j++)
		debug << partition_code_table[generic_code_data::table::TYPE][j] << ' ';
	debug << endl;
	debug << "partition_peer_code:  label: ";
	for (int j=0; j< num_partition_crossings; j++)
		debug << partition_code_table[generic_code_data::table::LABEL][j] << ' ';
	debug << endl;
	
	if (track_zig_zag_counts)
	{
		debug << "partition_peer_code: partition_zig_zag_count: " << endl;
		print(partition_zig_zag_count,debug, 3,"partition_peer_code: ");	
	}
}
			generic_code_data partition_code_data;
			partition_code_data.type = code_data.type;
			partition_code_data.immersion_character = code_data.immersion_character;
			partition_code_data.num_components = num_partition_components;
			partition_code_data.num_crossings = num_partition_crossings;
			partition_code_data.code_table = partition_code_table;
			partition_code_data.first_edge_on_component = first_edge_on_partition_cpt;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "partition_peer_code: partition_code_data: " << endl;
	print_code_data(debug,partition_code_data,"partition_peer_code:  ");
}
						
			ostringstream oss;
			write_code_data(oss, partition_code_data);
			
			if (i==0)
			{
				read_code_data(first_code_data, oss.str());

				first_code_data.immersion_character = code_data.immersion_character;
				first_code_data.head_zig_zag_count = code_data.head_zig_zag_count;				
				
				if (pure_knotoid_code_data)
				{					
					if (head_semi_arc % 2)
						first_code_data.head = first_code_data.code_table[generic_code_data::table::EPEER][(head_semi_arc-1)/2]/2;
					else
						first_code_data.head = head_semi_arc/2;
						
					valid_knotoid_input(first_code_data);  // sets the shortcut crossings
				}
					
			}
			else
			{
				read_code_data(subsequent_code_data, oss.str());			
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
			
if (debug_control::DEBUG >= debug_control::BASIC)
{
	if (i==0)
	{
		debug << "partition_peer_code:   first partition peer code ";
		write_code_data(debug,first_code_data);
		debug << endl;
		debug << "partition_peer_code:   first partition code_data" << endl;
		print_code_data(debug, first_code_data, "partition_peer_code:   ");
		debug << "partition_peer_code:   component_flags updated to: ";
		for (unsigned int j=0; j< component_flags.size(); j++)
			debug << component_flags[j] << ' ';
		debug << endl;
	}
	else
	{
		debug << "partition_peer_code:   subsequent partition peer code ";
		write_code_data(debug,subsequent_code_data);
		debug << endl;
		debug << "partition_peer_code:   subsequent partition code_data" << endl;
		print_code_data(debug, subsequent_code_data, "partition_peer_code:   ");
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
void cycle_gauss_code_labels(generic_code_data& gauss_code_data, int edge)
{
	int num_crossings = gauss_code_data.num_crossings;
	int num_components = gauss_code_data.num_components;
	matrix<int>& gauss_table = gauss_code_data.code_table;
	vector<int>& num_component_edges = gauss_code_data.num_component_edges;
	vector<int>& first_edge_on_component = gauss_code_data.first_edge_on_component;

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "cycle_gauss_code_labels: initial gauss_code_data :" << endl;
	print_code_data(debug,gauss_code_data,"cycle_gauss_code_labels: ");	
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

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "cycle_gauss_code_labels: edge lies on component " << component << ", first edge on component = " << first << ", last edge = " << last << endl;

    /* we record the new labels in opeer and epeer */
    vector<int> opeer(num_crossings);
    vector<int> epeer(num_crossings);
    
    for (int i=0; i< num_crossings; i++)
    {
		opeer[i] = gauss_table[generic_code_data::table::OPEER][i];
		epeer[i] = gauss_table[generic_code_data::table::EPEER][i];
	}

	/* cycle the labels in opeer and epeer belonging to the identified component */
	for (int i=first; i <= last; i++)
	{
		/* find edge i in the code table */
		for (int j=0; j< num_crossings; j++)
		{
			if (gauss_table[generic_code_data::table::EPEER][j] == i)
			{
				if (i == last)
					epeer[j] = first;
				else
					epeer[j]++;
					
				gauss_table[generic_code_data::table::ODD_ORIGINATING][j] = (epeer[j]+1-first)%num_component_edges[component] + first;
				
				break;
			}
			else if (gauss_table[generic_code_data::table::OPEER][j] == i)
			{
				if (i == last)
					opeer[j] = first;
				else
					opeer[j]++;
					
				gauss_table[generic_code_data::table::EVEN_ORIGINATING][j] = (epeer[j]+1-first)%num_component_edges[component] + first;

				break;
			}
		}
	}

    for (int i=0; i< num_crossings; i++)
    {
		gauss_table[generic_code_data::table::OPEER][i] = opeer[i];
		gauss_table[generic_code_data::table::ODD_TERMINATING][i] = opeer[i];
		gauss_table[generic_code_data::table::EPEER][i] = epeer[i];
		gauss_table[generic_code_data::table::EVEN_TERMINATING][i] = epeer[i];
	}


if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "cycle_gauss_code_labels: resultant gauss_code_data :" << endl;
	print_code_data(debug,gauss_code_data,"cycle_gauss_code_labels: ");	
}
	
}

/* cycle_gauss_code_labels shifts the start of the numbering of the component containing edge back by one, 
   so that the number of each edge in the component is incremented (around the component)   
*/
void cycle_gauss_code_labels(matrix<int>& crossing_peers, vector<int>& num_component_edges, vector<int>& first_edge_on_component, int num_crossings,int num_components, int edge)
{

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
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

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
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
	
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "cycle_gauss_code_labels: resultant crossing_peers:" << endl;
	print(crossing_peers,debug,4,"cycle_gauss_code_labels: ");	
}
	
}


/* identify_gauss_crossings returns a vector, gauss_crossing_map that maps Gauss crossing i to  gauss_crossing_map[i] */
void identify_gauss_maps(generic_code_data& code_data,vector<int>& gauss_crossing_map, vector<int>& gauss_arc_map)
{

	matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;
	int num_components = code_data.num_components;
	
	vector<int>& term_crossing = code_data.term_crossing;
//	int num_edges = 2*num_crossings;
//	vector<int> edge_flag(num_edges); // initialized to zero
	vector<bool> component_traced(num_components); // initialized to zero

		
	int num_classical_crossings = num_crossings;
	bool pure_knotoid_code_data = false;
	vector<int>& shortcut_crossing = code_data.shortcut_crossing;
	
	for (int i=0; i< num_crossings; i++)
	{
		if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::VIRTUAL)
			num_classical_crossings--;
	}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_maps: num_classical_crossings = " << num_classical_crossings << endl;

	if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && shortcut_crossing.size())
	{
		pure_knotoid_code_data = true;
		for (unsigned int i=0; i< shortcut_crossing.size(); i++)
		{
			if (shortcut_crossing[i] != 0)
				num_classical_crossings--;
		}
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "identify_gauss_maps: knotoid: num_classical_crossings = " << num_classical_crossings << endl;
	}

	int num_classical_crossings_visited = 0;
		
	/* gauss_crossing_map records the immersion crossings corresponding to Gauss crossings
	   so that gauss_crossing_map[i] is the immersion crossing corresponding to the ith Gauss crossing.
	   
	   crossing_visited will be a flag to indicate whether the ith immersion crossing
	   has been recorded in gauss_crossing_map or not.
	*/
	gauss_crossing_map = vector<int>(num_classical_crossings);
	vector<int> crossing_visited(num_crossings);
	gauss_arc_map = vector<int>(num_edges);
	
	vector<int> first_gauss_arc_on_component(num_components);
	vector<int> num_gauss_arcs_on_component(num_components);

	int component = 0;
	int component_count=1;
	int start=0;
	int edge=0;
	int gauss_arc = 0;
	bool complete = false;

	component_traced[0] = true;
		
	do 
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_maps: component_start = " << start << endl;
			
		first_gauss_arc_on_component[component] = gauss_arc;
				
		/*	trace this component */
		do
		{			
			gauss_arc_map[edge] = gauss_arc;
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_maps: edge = " << edge << ", gauss_arc " << gauss_arc;
//			edge_flag[edge] = 1;
			int next_crossing = term_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", next_crossing = " << next_crossing;
		
			if (code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::VIRTUAL)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", is virtual" << endl;
			}
			else if (pure_knotoid_code_data && shortcut_crossing[next_crossing])
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", is a shortcut crossing" << endl;
			}
			else if ((edge%2 != 0 && code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::POSITIVE) ||
			    (edge%2 == 0 && code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::NEGATIVE))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", going under" << endl;	
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", going over" << endl;
			}
				
			if (code_table[generic_code_data::table::LABEL][next_crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[next_crossing]))
			{
				if(crossing_visited[next_crossing])
				{
					for (int i=0; i< num_classical_crossings_visited; i++)
					{
						if (gauss_crossing_map[i] == next_crossing)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_maps:   second visit to Gauss crossing " << i+1<< endl;
	
							gauss_arc++;							
							break;
						}
					}
				}
				else
				{

					gauss_crossing_map[num_classical_crossings_visited] = next_crossing;
					crossing_visited[next_crossing] = 1;
					num_classical_crossings_visited++;
					gauss_arc++;
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_maps:   first visit, becomes Gauss crossing " << num_classical_crossings_visited << endl;
				}

				if (edge%2)
					edge = code_table[generic_code_data::table::EVEN_ORIGINATING][next_crossing];
				else
					edge = code_table[generic_code_data::table::ODD_ORIGINATING][next_crossing];					
			}
			else
			{
				/* just move on around the component */				
				if (edge%2)
					edge = code_table[generic_code_data::table::EVEN_ORIGINATING][next_crossing];
				else
					edge = code_table[generic_code_data::table::ODD_ORIGINATING][next_crossing];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_maps:   doing nothing" << endl;
			}				
		} while (edge != start);
		
		/* at the end of the component, gauss_arc records the first edge of the next component, since we have returned to the start 
		   so we can determine the number of Gauss arcs on this component and reset the gauss_arc_map where it set to gauss_arc.
		*/
		num_gauss_arcs_on_component[component] = gauss_arc - first_gauss_arc_on_component[component];
		for (int i=0; i< num_edges; i++)
		{
			if (gauss_arc_map[i] == gauss_arc)
				gauss_arc_map[i] = first_gauss_arc_on_component[component];
		}
		

		/* look for the start of another component */
		complete = true;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_maps: look for start of next component" << endl;
	
		for (int i=0; i< num_components; i++)
		{
			if (component_traced[i] == false)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_maps:   component " << i << " not traced";

				if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID && component_count < code_data.num_open_components)
				{
					if (code_data.component_type[i].type == component_character::CLOSED)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", CLOSED component but component_count = " << component_count << " so skipping" << endl;
						continue;  // conside all open components first
					}
				}
		
				component = i;
				component_traced[i] = true;
				component_count++;
				complete = false;

				if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID && 
				    (code_data.component_type[i].type == component_character::PURE_END_LEG || code_data.component_type[i].type == component_character::KNOT_TYPE_END_LEG))					
					start = code_data.first_edge_on_component[i]+code_data.num_component_edges[i]-1;
				else
					start = code_data.first_edge_on_component[i];

				edge = start;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", acceptable, component_count = " << component_count << " starting edge on component = " << start << endl;

				break;
			}
		}

	} while (!complete);
}

/* The function vertex_span counts the number of vertices that can be reached from the initial crossing in the given generic code data.  
   This can be used to determine whether the generic code data can be realized by a connected immersion.  

   If the boolean gauss_arcs is true, the function uses Gauss arcs rather than immersion arcs to move between vertices and builds up a list of
   reachable Gauss crossings.
   
   The function takes no account of the labels attached to the crosings, nor the character of the immersion, so would include shortcut
   edges in the case of a knotoid.  Such immersions may be handled using the optional exclude parameter.
   
   Optionally, the exclude vector may specify a number of edges to exclude from the search, this may be use to determine whether 
   the complement of the excluded edges is connected.
*/
list<int> vertex_span (generic_code_data& code_data, int initial_crossing, vector<int>* exclude)
{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "vertex_span: checking code ";
	write_code_data (debug, code_data);
	debug << endl;
	if (exclude != 0 && exclude->size())
	{
		debug << "vertex_span: exclude: ";
		for (unsigned int i=0; i< exclude->size(); i++)
			debug << (*exclude)[i] << ' ';
		debug << endl;
	}
	
	debug << "vertex_span: starting from crossing " << initial_crossing << endl;
}
	
	matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;
	list<int> vertex_list;
	vertex_list.push_back(initial_crossing);
	list<int>::iterator lptr = vertex_list.begin();
	vector<int> crossing_flag(num_crossings,0);
	crossing_flag[initial_crossing] = 1;
	vector<int> edge_flag(2*num_crossings,0);
	vector<int>& orig_crossing = code_data.orig_crossing;
	vector<int>& term_crossing = code_data.term_crossing;


	if (exclude != 0 && exclude->size())
	{
		
		for (unsigned int i=0; i< exclude->size(); i++)
			edge_flag[(*exclude)[i]] = 1;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "vertex_span: exclude edges ";
    for (int i=0; i< 2*num_crossings; i++)
		debug << edge_flag[i] << ' ';
    debug << endl;
}
	}
	else
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span: no excluded edges " << endl;
	}
	
	while (lptr != vertex_list.end())
	{
		/* Look along any edges not already considered to see if we can extend the list of vertices.  
		   We look back along terminating edges and forwards along originating edges.  If
		   edge_flag[i] is 1 then we have either looked along this edge already or have excluded it 
		   from the search; crossing_flag[i] indicates whether we already have crossing i on the list.
	    */
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span:   from crossing " << *lptr << endl;

		int edge = code_table[generic_code_data::table::EVEN_TERMINATING][*lptr];
		int crossing;
		
		if (edge_flag[edge] == 0)
		{
			crossing = orig_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span:     going back along " << edge << " takes us to " << crossing;

		    if(crossing_flag[crossing] == 0)
		    {
				vertex_list.push_back(crossing);
				crossing_flag[crossing] = 1;

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", adding " << crossing << " to vertex_list" << endl;
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", which is already on the list " << endl;
			}
			edge_flag[edge] = 1;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span:     don't need to consider edge " << edge << endl;
		}
		
		edge = code_table[generic_code_data::table::ODD_TERMINATING][*lptr];
		
		if (edge_flag[edge] == 0)
		{
			crossing = orig_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span:     going back along " << code_table[generic_code_data::table::ODD_TERMINATING][*lptr] << " takes us to " << crossing;

		    if(crossing_flag[crossing] == 0)
		    {
				vertex_list.push_back(crossing);
				crossing_flag[crossing] = 1;

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", adding " << crossing << " to vertex_list" << endl;
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", which is already on the list " << endl;
			}
			edge_flag[edge] = 1;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span:     don't need to consider edge " << edge << endl;
		}
		
		edge = code_table[generic_code_data::table::EVEN_ORIGINATING][*lptr];

		if (edge_flag[edge] == 0)
		{
			crossing = term_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span:     going forwards along " << code_table[generic_code_data::table::EVEN_ORIGINATING][*lptr] << " takes us to " << crossing;

		    if(crossing_flag[crossing] == 0)
		    {
				vertex_list.push_back(crossing);
				crossing_flag[crossing] = 1;

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", adding " << crossing << " to vertex_list" << endl;
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", which is already on the list " << endl;
			}
			edge_flag[edge] = 1;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span:     don't need to consider edge " << edge << endl;
		}
		
		edge = code_table[generic_code_data::table::ODD_ORIGINATING][*lptr];
		if (edge_flag[edge] == 0)
		{
			crossing = term_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span:     going forwards along " << code_table[generic_code_data::table::ODD_ORIGINATING][*lptr] << " takes us to " << crossing;

		    if(crossing_flag[crossing] == 0)
		    {
				vertex_list.push_back(crossing);
				crossing_flag[crossing] = 1;

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", adding " << crossing << " to vertex_list" << endl;
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", which is already on the list " << endl;
			}
			edge_flag[edge] = 1;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span:     don't need to consider edge " << edge << endl;
		}
		
		lptr++;
	}
	return vertex_list;
}

/* calculate_turning_cycles looks for left and right turning cycles in the given code_data and returns
   false if it is unable to calculate a set of cycles that may be realized.  From Eulers formula a
   realizable diagram cannot have more than num_crossings+2 turning cycles.  If the search therefore 
   exceeds this number of cycles the code_data cannot be realizable and the function returns false.
   If the search completes without reaching this limit the function returns true but this does not indicate
   that the code_data is necessarily realizable.  Note that the num_crossings+2 limit can only be breached
   when looking for right turning cycles, given that we look for left turning cycles first.
*/
bool calculate_turning_cycles(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles)
{

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "calculate_turning_cycles: code data: ";
	write_peer_code (debug, code_data);
	debug << endl;
	print_code_data(debug,code_data,"calculate_turning_cycles: ");	
}

	matrix<int>& code_table = code_data.code_table;
	vector<int>& term_crossing = code_data.term_crossing;
	vector<int>& orig_crossing = code_data.orig_crossing;
	
	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;

	num_cycles = 0;
	
	/* First look for left turning cycles */
	for (int i=0; i<2*num_crossings; i++)
	{
		
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "calculate_turning_cycles: edge = " << i;

    	/* does edge i already appear in a cycle ? */
		bool found = false;
		for (int j=0; j<num_cycles; j++)
		{
			for (int k=1; k<= cycle[j][0]; k++)
			{
				if (abs(cycle[j][k]) == i)
				{

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << " found in left turning cycle " << j << endl;

					found = true;
					break;
				}
			}
		}

		if (!found)
		{

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << " not found in current left turning cycles "<< endl;

			/* start a new cycle */
			int column = 1;

			/* we always traverse odd edges backwards */
			int edge = (i%2? -i: i);
			cycle [num_cycles][column++] = edge;
			bool complete = false;

			/* a cycle cannot be longer than num_edges so we check that we do indeed
			   cycle within that limit.  This is used to check that the component map
			   within the code_table is valid, since an unrealizable component map can
			   result in infinite loops being calculated in the terminating and 
			   originating edges of the crossings where we never return to the start
			   of a cycle.
			*/
			for (int j=0; !complete && j< num_edges; j++)
			{
				/* determine which vertex the current edge takes us to and 
	   			the next edge we turn onto */
				if (edge % 2) // edge is odd (and negative)
				{
					int vertex = orig_crossing[-edge];
					if (code_table[generic_code_data::table::TYPE][vertex] == generic_code_data::TYPE1)
		    			edge = code_table[generic_code_data::table::EVEN_ORIGINATING][vertex];
					else
						edge = -code_table[generic_code_data::table::ODD_TERMINATING][vertex];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;

				}
				else // edge is even (and positive)
				{
					int vertex = term_crossing[edge];
					if (code_table[generic_code_data::table::TYPE][vertex] == generic_code_data::TYPE1)
						edge = -code_table[generic_code_data::table::ODD_TERMINATING][vertex];
					else
		    			edge = code_table[generic_code_data::table::EVEN_ORIGINATING][vertex];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;
				}

    
				if (edge == cycle[num_cycles][1])
				{
					complete = true;
					cycle[num_cycles][0] = column-1;
					num_cycles++;
				}
				else
				{
					cycle[num_cycles][column++] = edge;
				}				
			}
			
			if (!complete)
			{
				/* we've encounterd an infinte loop */
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "calculate_turning_cycles: exceeded the maximum possible length of a turning cycle in a realizable code" << endl;
				return false;
			}
		}
	}	

	/* record the number of left cycles */
	num_left_cycles = num_cycles;
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
    debug << "calculate_turning_cycles: number of left turning cycles = " << num_left_cycles;
	for (int i=0; i<num_left_cycles; i++)
	{
		debug << "\ncalculate_turning_cycles: cycle " << i << " length = " << cycle[i][0] << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " ";
	}
	debug << endl;
}
		
	/* Now look for right turning cycles */

	for (int i=0; i<2*num_crossings; i++)
	{

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "calculate_turning_cycles: edge = " << i;

		/* does edge i already appear in a right cycle ? */
		bool found = false;
		for (int j=num_left_cycles; j<num_cycles; j++)
		{
			for (int k=1; k<= cycle[j][0]; k++)
			{
				if (abs(cycle[j][k]) == i)
				{
					
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << " found in right turning cycle " << j << endl;
    
    					found = true;
					break;
				}
			}
		}

		if (!found)
		{

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << " not found in current right turning cycles "<< endl;

			/* check we've not exceeded the maximum number of turning cycles */
			if (num_cycles == num_crossings+2)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "calculate_turning_cycles: exceeded the maximum possible number of turning cycles in a realizable code" << endl;
				return false;
			}
			
			/* start a new cycle */
			int column = 1;
			
			/* we always traverse odd edges backwards */
			int edge = (i%2? -i: i);
			cycle [num_cycles][column++] = edge;
			bool complete = false;

			do
			{
				/* determine which vertex the current edge takes us to and 
	   			the next edge we turn onto */
				if (edge % 2) // edge is odd (and negative)
				{
					int vertex = orig_crossing[-edge];
					if (code_table[generic_code_data::table::TYPE][vertex] == generic_code_data::TYPE1)
						edge = -code_table[generic_code_data::table::ODD_TERMINATING][vertex];
					else
		    			edge = code_table[generic_code_data::table::EVEN_ORIGINATING][vertex];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;

				}
				else // edge is even (and positive)
				{
					int vertex = term_crossing[edge];
					if (code_table[generic_code_data::table::TYPE][vertex] == generic_code_data::TYPE1)
		    			edge = code_table[generic_code_data::table::EVEN_ORIGINATING][vertex];
					else
						edge = -code_table[generic_code_data::table::ODD_TERMINATING][vertex];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;

				}

				if (edge == cycle[num_cycles][1])
				{
					complete = true;
					cycle[num_cycles][0] = column-1;
					num_cycles++;
				}
				else
				{
					cycle[num_cycles][column++] = edge;
				}				
			} while(!complete);			
		}
	}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
    debug << "calculate_turning_cycles: number of right turning cycles = " << num_cycles - num_left_cycles;
	for (int i=num_left_cycles; i<num_cycles; i++)
	{
		debug << "\ncalculate_turning_cycles: cycle " << i << " length = " << cycle[i][0] << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " "	;
	}
	debug << endl;
}
	return true;
}


bool realizable_code_data(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles)
{
	int num_crossings = code_data.num_crossings;

	if (!calculate_turning_cycles(code_data, cycle, num_left_cycles, num_cycles))
		return false;

	if (num_cycles != num_crossings+2)
		return false;
	
	return true;
}

/* The generic code data code table and crossing number identify a valid knotoid if the
   specified crossing is the start of a shortcut.  A shortcut from the head
   to the leg of a knotoid must pass under all the strands of the knotoid; thus
   if the label associated with the head crossing is '+' then the head lies on the
   odd-numbered semi-arc of the diagram, and if the label is '-' it lies
   on the even-numbered semi-arc.  From this arc onwards all the crossings we 
   encounter must be approached on the under-arc if this is a valid knotoid 
   combination of code_table and head.
   
   A shortcut is an embeded arc, so it is not permitted to contain self-intersections.

   Furthermore, we require that the semi-arc containing the leg of a knotoid is labelled
   zero.  This is required so that in the case of a multi-knotoid we can identify the segment 
   component of a smoothed diagram when calculating the Turaev extended bracket polynomial.

   Allowing arbitary numbering would require the identification of the leg as well as the head within 
   the peer code, since one could not guarantee that the segment component numbering could start with 
   the leg and preserve the requirement for odd and even terminating edges at each crossing.
       
   As we check the validity of the input we note the shortcut crossings.  We set shortcut_crossing[i] 
   to be non-zero if the ith crossing is a shortcut crossing.  We set it to 1 if the knotoid crosses the 
   shortcut from right to left and -1 if it crosses from left to right, orienting the shortcut from leg to head.
   Note, however, that the shortcut is oriented from the leg to the head and we will be traversing it in 
   the opposite direction.  
*/
bool valid_knotoid_input(generic_code_data& code_data)
{
	matrix<int>& code_table = code_data.code_table;
	vector<int> shortcut_crossing(code_data.num_crossings);
	int head = code_data.head;
	
	/* we must have been given a crossing indicating the position of the head */
	if (head == -1)
	{
		
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "valid_knotoid_input: no indication of knotoid head provided, returning false" << endl;
	
		return false;
	}
		
	
	int semi_arc;	
	int peer;


if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "valid_knotoid_input: code_table " << endl;
	print(code_table,debug,3,"valid_knotoid_input: ");	
	debug << "valid_knotoid_input: head =" << head << endl;
}

	int component;
	if (code_table[generic_code_data::table::LABEL][head] == generic_code_data::POSITIVE)
	{
		semi_arc = code_table[generic_code_data::table::OPEER][head];
		peer = 2*head;
		component = code_table[generic_code_data::table::COMPONENT][(semi_arc-1)/2];
	}
	else if (code_table[generic_code_data::table::LABEL][head] == generic_code_data::NEGATIVE)
	{
		semi_arc = 2*head;
		peer = code_table[generic_code_data::table::OPEER][head];
		component = code_table[generic_code_data::table::COMPONENT][head];
	}
	else
	{
		/* the first shortcut crossing cannot be virtual */
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input: head " << 0 << " indicates first shortcut crossing is virtual" << endl;
		return false;
	}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "valid_knotoid_input: head " << head << " lies on semi-arc " << semi_arc << " in component " << component << " with peer edge " << peer << endl;
	
	if (component != 0)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input: invalid knotoid code, the knotoid leg must lie on the semi-arc numbered 0" << endl;
		return false;
	}

	if (semi_arc < peer && peer < code_data.num_component_edges[0])
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input: invalid knotoid code, peer edge also lies in the shortcut, which may not self-intersect" << endl;
		return false;
	}

    /* set the shortcut crossing flag for the head crossing */
    if (semi_arc % 2)
    {
		if (code_table[generic_code_data::table::TYPE][head] == generic_code_data::TYPE1)				
		{
			shortcut_crossing[head] = -1; 
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   knotoid crosses shortcut at first shortcut crossing from the left" << endl;
		}
		else
		{
			shortcut_crossing[head] = 1;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   knotoid crosses shortcut at first shortcut crossing from the right" << endl;
		}
	}
	else
	{
		if (code_table[generic_code_data::table::TYPE][head] == generic_code_data::TYPE1)				
		{
			shortcut_crossing[head] = 1; 
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   knotoid crosses shortcut at first shortcut crossing from the right" << endl;
		}
		else
		{
			shortcut_crossing[head] = -1;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   knotoid crosses shortcut at first shortcut crossing from the left" << endl;
		}
	}
		
	bool valid = true;

	for (int i = semi_arc+1; i< code_data.num_component_edges[0]; i++)
	{

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input: semi-arc " << i << endl;
	
		/* check that semi-arc i is the under-arc at the next crossing */
		if (i%2)
		{				
			int crossing = code_data.term_crossing[i];
			peer = code_table[generic_code_data::table::EVEN_TERMINATING][crossing];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   crossing " << crossing << ", peer " << peer << endl;
			
			/* peer edge must not lie in the shortcut and the label for this crossing must be '+' if odd arc is under-arc */
			if (semi_arc < peer && peer < code_data.num_component_edges[0])
			{
				valid = false;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input: invalid knotoid code, peer edge also lies in the shortcut, which may not self-intersect" << endl;
				break;
			}
			else if (code_table[generic_code_data::table::LABEL][crossing] != generic_code_data::POSITIVE)
			{
				valid = false;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   not positive, as required for valid knotoid" << endl;
				break;
			}
			else
			{
			    /* The sortcut is oriented from the leg to the head, so we are tracing it in the
				   reverse direction.  Thus, we must view the knotoid from the perspective of 
				   looking back at the crossing from the originating shortcut edge
				*/
				if (code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::TYPE1)				
				{
					shortcut_crossing[crossing] = -1; 
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   OK, knotoid crosses from the left" << endl;
				}
				else
				{
					shortcut_crossing[crossing] = 1;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   OK, knotoid crosses from the right" << endl;
				}
					
			}
		}
		else
		{
			int crossing = i/2 ;
			peer = code_table[generic_code_data::table::ODD_TERMINATING][crossing];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   crossing " << crossing << ", peer " << peer << endl;

			/* peer edge must not lie in the shortcut and the label for this crossing must be '-' if even arc is under-arc */
			if (semi_arc < peer && peer < code_data.num_component_edges[0])
			{
				valid = false;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input: invalid knotoid code, peer edge also lies in the shortcut, which may not self-intersect" << endl;
				break;
			}
			else if (code_table[generic_code_data::table::LABEL][crossing] != generic_code_data::NEGATIVE)
			{
				valid = false;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   not negative, as required for valid knotoid" << endl;
				break;
			}
			else
			{
			    /* see comment above for the odd case */
				if (code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::TYPE1)				
				{
					shortcut_crossing[crossing] = 1; 
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   OK, knotoid crosses from the right" << endl;
				}
				else
				{
					shortcut_crossing[crossing] = -1;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   OK, knotoid crosses from the left" << endl;
				}
			}
		}
	}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input: returning " << valid << endl;
	
	code_data.shortcut_crossing = shortcut_crossing;
	return valid;
}

/* The over_preferred_gauss_code is a unique Gauss code representation of a virtual diagram, suitable for
   virtual knots and links, welded knots and links, virtual knotoids and multiknotoids.  It mimics the
   left preferred Gauss code defined for doodles with preference given to over-crossings rather than under-
   crossings but adds a comparison of crossing signs, with positive crossings preferred to negative crossings.
   
   The over_preferred_gauss_code for a diagram is a Gauss code describing the diagram that contains
   O1, ..., On as a subsequence and is lexicographically minimal amongst all such codes.  If the 
   boolean unoriented is true, we determine the over_preferred_gauss_code by considering both orientations 
   of each unicursal component in all combinations.
   
   For each given orientation of the unicursal components, we consider all possible orderings of the components
   and all possible starting positions around the component relative to the starting location implicit within the 
   code_data parameter supplied to the function.  For each combination of component ordering and starting position, 
   we construct a vector<pair<int,int> > representation of the Gauss code of the form {<new_crossing_number>,<O|U>}
   where the new_crossing number is determined by the order in which we encounter the crossings.  We then append the 
   crossing signs obtained by that choice of orientation and starting position and compare the concatenation of the vector
   and signs lexicographcally.  That is, we use the lexicographical vector comparison operator to compare the vectors and
   if they are equal, compare the signs manually (preferring + over -) to select the minimal representation.
   
   For links, we have to consider the case that a diagram admits two different Gauss codes having the same sequence
   {<crossing>,<O|U>} and signs.  For example, the unoriented immersion [13 -15 17 -1 3 -5, 11 7 -9] / * + - - + * + * - admits the
   two Gauss code descriptions 
   
           1 2 -3 -5 -1 3 4 -6,-2 -4 5 6/+ - + + - +  and
           1 2 -3 -5,-1 3 4 -6 -2 -4 5 6/+ - + + - +
   
   which differ only in the position of the component separator.  To handle this we consider the number of terms in each 
   component and prefer the smallest vector lexicographically.  Thus, the first of the above codes has component lengths
   6 and 4, whereas the second has component lengths 4 and 6, so we prefer the second form.
   
   Every Gauss code has a unique over-preferered representation, thus there is a unique over-preferred description of each 
   diagram.  However, reflecting a diagram in a disjoint line of the plane results in a diagram that may differ from the 
   original but has the same Gauss code, so the over-preferred Gauss code does not determine a unique diagram
*/
string over_preferred_gauss_code(generic_code_data& code_data, bool unoriented)
{

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "over_preferred_gauss_code: initial code data ";	
	write_code_data(debug,code_data);
	debug << endl;
	print_code_data(debug, code_data, "over_preferred_gauss_code: ");
}

	generic_code_data g; // Gauss code data
	
	if (code_data.type == generic_code_data::code_type::gauss_code)
	{
		g = code_data;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "over_preferred_gauss_code: gauss_code_data is the same as code_data";	

	}
	else
	{
		ostringstream oss;
		write_gauss_code(oss, code_data);
		read_gauss_code(g, oss.str());
		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "over_preferred_gauss_code: gauss_code_data ";	
	write_code_data(debug,g);
	debug << endl;
	print_code_data(debug, g, "over_preferred_gauss_code: ");
}
		
	}
	
	matrix<int>& code_table = g.code_table;
	vector<pair<int,int> > min_perm_weight;  // initializes to zero size
	vector<int> min_perm;
	vector<char> min_perm_sign;
	vector<int> min_component_perm;
	vector<int> min_component_cycle;
	vector<int> min_component_length; // the number of terms in each component when considered in the order determined by min_component_perm
	vector<int> min_term_crossing;
	vector<bool> min_r_flag;
	
	int num_crossings = g.num_crossings;
	int num_components = g.num_components;
	int num_terms=2*num_crossings;
	int min_component = (g.immersion_character == generic_code_data::character::CLOSED? 0: 1);

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "over_preferred_gauss_code: min_component = " << min_component << endl;			
	
	vector<bool> r_flag(num_components); // initializes to false
	bool finished = false;
	
	vector<int> initial_term_crossing = g.term_crossing;
	vector<int> term_crossing = initial_term_crossing;

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "over_preferred_gauss_code: initial_term_crossing: ";	
	for (int i=0; i< num_terms; i++)
		debug << initial_term_crossing[i] << ' ';
	debug << endl;
}

	/* identify the components associated with the first and second visits to crossing i in the original numbering */
	vector<int> first_component(num_crossings);
	vector<int> second_component(num_crossings);
	
	for (int i=0; i< num_crossings; i++)
	{
		int component = 0;
		for (int j=1; j< num_components; j++)
		{
			if (code_table[generic_code_data::table::OPEER][i] >= g.first_edge_on_component[j])
				component++;
			else
				break;
		}
		first_component[i] = component;
	
		component = 0;
		for (int j=1; j< num_components; j++)
		{
			if (code_table[generic_code_data::table::EPEER][i] >= g.first_edge_on_component[j])
				component++;
			else
				break;
		}
		second_component[i] = component;
		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "over_preferred_gauss_code: original crossing " << i << " first visit on component " << first_component[i] 
	      << " second visit on component " << second_component[i] << endl;	
}
	}
	
	
	do
	{
		/* consider the orientation determined by the r_flags */
		
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "over_preferred_gauss_code: reverse flags: ";
	for (int i=0; i< num_components; i++)
		debug << r_flag[i] << " ";
	debug << endl;
}
			
		/* consider the components in every posible order */
		vector<int> component_perm(num_components);
		for (int i=0; i< num_components; i++)
			component_perm[i]=i;
		vector<int>::iterator component_perm_start = component_perm.begin();
		
		if (min_component !=0)
			component_perm_start++;
	
		do
		{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "over_preferred_gauss_code:   component_perm = ";	
	for (int i=0; i< num_components; i++)
		debug << component_perm[i] << ' ';
	debug << endl;
}	

			vector<int> component_length(num_components);
			for (int i=0; i< num_components; i++)
				component_length[i] = g.num_component_edges[component_perm[i]];
			
			/* component_cycle will record the cyclic rotation of each component, enumerated
			   in the order given by g, by recording which term of the component we are going 
			   to start at.  We will then need to find this instance and evaluate the starting 
			   offset in g.orientation_matrix
			*/
			vector<int> component_cycle(num_components); // initializes to zero
	
			bool complete = false;		
			do
			{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "over_preferred_gauss_code:     component_cycle = ";	
	for (int i=0; i< num_components; i++)
		debug << component_cycle[i] << ' ';
	debug << endl;
}	
	
				/* evaluate the crossing permutation: perm[i] indicates the new crossing number of the old ith crossing */
				vector<int> perm(num_crossings);
				int new_crossing = 0;
				int new_edge = -1;
				for (int k=0; k< num_components; k++)
				{
					int component = component_perm[k];
					int num_component_terms = g.num_component_edges[component];

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "over_preferred_gauss_code:       component = " << component << ", num_component_terms = " 
		  << num_component_terms << endl; //" <edge>,<terminating_crossing><O|U>: " << endl;
}	
					/* r_flag[component] tells us whether this component has been reversed.  If it has, the entries in term_crossing corresponding to 
					   this component have been reversed (based on the order in initial_term_crossing).  Thus, term_crossing[new_edge] identifies the
					   crossing at which we arrive on new_edge if we renumbered the diagram using the current r_flags component_perm, using 
					   component_cycle to determine the starting point on each component.  
					   
					   If the component has not been reversed, the index into term_crossing identifies the old edge number (i.e. as used in the 
					   initial code data numbering) that arrives at that crossing.  We need to know this old edge number in order to determine 
					   whether we are arriving at the crossing on an over-arc or an under-arc.
					   
					   If the component has been reversed we have the following situation in term_crossing:
					   
					     first_edge_on_component
						            |
						            |      component_cycle
					                |=======>|
					                v        v
					         -------+---------------------------+-----------
							        |<---num_component_terms--->|
					         -------+---------------------------+-----------								
							                     forwards
							         ------->o------------------=======> new edge numbering
									                             wrap																 
							                     backwards
							         <-------o<-----------------======= old edge numbering
									                             wrap
																 
					   Thus, the effective component_cycle for the old edge numbering is (num_component_terms-component_cycle)
					   and we move backwards around the old numbering.
					*/
					for (int i=0; i < num_component_terms; i++)
					{
//						int new_edge = g.first_edge_on_component[component] + i;
						new_edge++;
						int index = (component_cycle[component]+i)%num_component_terms + g.first_edge_on_component[component];
						
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "over_preferred_gauss_code:         new edge " << new_edge << " index = " << index;
									
						int terminating_crossing = term_crossing[index];
						
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << ", terminating_crossing  = "<< terminating_crossing << endl;
						
						bool over_arc = true;
						int old_edge = index;
						
						/* If we have reversed the orientation of this component, we need to look for the succeeding edge around the 
						   component in the code table to determine whether we are on the over_arc or under arc.
						*/
						if (r_flag[component])
						{
							int reverse_component_cycle = (num_component_terms-1)-component_cycle[component];
							int old_offset = reverse_component_cycle-i;
							int old_cpt_edge = (old_offset+num_component_terms)%num_component_terms;
							old_edge =  old_cpt_edge + g.first_edge_on_component[component];
							
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "over_preferred_gauss_code:         reversed component, reverse_component_cycle = " << reverse_component_cycle << " old_offset = " << old_offset << " old_cpt_edge = " << old_cpt_edge << " old_edge = " << old_edge << endl;
						}
						
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "over_preferred_gauss_code:         old_edge = " << old_edge << " code_table[generic_code_data::table::OPEER][terminating_crossing] = " << code_table[generic_code_data::table::OPEER][terminating_crossing] << " code_table[generic_code_data::table::EPEER][terminating_crossing] = " << code_table[generic_code_data::table::EPEER][terminating_crossing] << endl;
}

						if (code_table[generic_code_data::table::OPEER][terminating_crossing] == old_edge)
						{	
							if (code_table[generic_code_data::table::LABEL][terminating_crossing] == generic_code_data::POSITIVE)
								over_arc = false;
						}
						else if (code_table[generic_code_data::table::EPEER][terminating_crossing] == old_edge)
						{					
							if (code_table[generic_code_data::table::LABEL][terminating_crossing] == generic_code_data::NEGATIVE)
								over_arc = false;
						}

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "over_preferred_gauss_code:           over_arc = " << over_arc << endl;
						
						if (over_arc) // if we're going over
						{
							perm[terminating_crossing] = new_crossing++;
						}
					}
				}

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "over_preferred_gauss_code:       perm = ";	
	for (int i=0; i< num_crossings; i++)
		debug << perm[i] << ' ';
	debug << endl;
}

				/* Evaluate the lexicographic weight of the permutation. The perm weight is a vector of pairs
				   (<permuted crossing number>,<O|U>) calculated from the start of each component, we
				   record O and U as 0 and 1 respectively so that O is preferred (minimal)
				*/
				vector<pair<int,int> > perm_weight(num_terms);
				int place = 0;
				for (int k=0; k< num_components; k++)
				{
					int component = component_perm[k];
					int num_component_terms = g.num_component_edges[component];
	
					for (int i=0; i < num_component_terms; i++)
					{
	//					int new_edge = g.first_edge_on_component[component] + i;
						int index = (component_cycle[component]+i)%num_component_terms + g.first_edge_on_component[component];								
						int terminating_crossing = term_crossing[index];				
						bool over_arc = true;
						int old_edge = index;
						
						if (r_flag[component])
						{
							int reverse_component_cycle = (num_component_terms-1)-component_cycle[component];
							int old_offset = reverse_component_cycle-i;
							int old_cpt_edge = (old_offset+num_component_terms)%num_component_terms;
							old_edge =  old_cpt_edge + g.first_edge_on_component[component];
						}

						if (code_table[generic_code_data::table::OPEER][terminating_crossing] == old_edge)
						{	
							if (code_table[generic_code_data::table::LABEL][terminating_crossing] == generic_code_data::POSITIVE)
								over_arc = false;
						}
						else if (code_table[generic_code_data::table::EPEER][terminating_crossing] == old_edge)
						{					
							if (code_table[generic_code_data::table::LABEL][terminating_crossing] == generic_code_data::NEGATIVE)
								over_arc = false;
						}
	
						perm_weight[place] = {perm[terminating_crossing],(over_arc? 0:1)};
						place++;
					}
				}
	
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "over_preferred_gauss_code:       perm_weight = ";
	for (int i=0; i< num_terms; i++)
		debug << '(' << perm_weight[i].first << "," << (perm_weight[i].second?'U':'O') << ")" << ' ';
	debug << endl;
}

				/* to determine the crossing signs we have to identify the components involved at each crossing and see whether they have been reversed.  
				   If exactly one of the components is reversed, the sign changes from the original numbering, otherwise it remains the same.  We need to 
				   determine the inverse permutation to min_perm, so that as we consider the crossing in the new order, we can identify the old crossing 
				   number for the purpose of identifying the components and the original sign.
				*/
			
				vector<char> perm_sign(num_crossings,'!');
				vector<int> inv_perm(num_crossings);
				for (int i=0; i< num_crossings; i++)
					inv_perm[perm[i]] = i;
					
				for (int i=0; i< num_crossings; i++)
				{
					bool reverse_sign = false;
					if ((r_flag[first_component[inv_perm[i]]] && !r_flag[second_component[inv_perm[i]]]) || 
					    (r_flag[second_component[inv_perm[i]]] && !r_flag[first_component[inv_perm[i]]])
					   )
						reverse_sign = true;
						
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "over_preferred_gauss_code:         reverse sign at crossing " << inv_perm[i] << ' ' << reverse_sign << endl;	
									
					if ( (code_table[generic_code_data::table::LABEL][inv_perm[i]] == generic_code_data::POSITIVE && code_table[generic_code_data::table::TYPE][inv_perm[i]] == generic_code_data::TYPE2) ||
					     (code_table[generic_code_data::table::LABEL][inv_perm[i]] == generic_code_data::NEGATIVE && code_table[generic_code_data::table::TYPE][inv_perm[i]] == generic_code_data::TYPE1)
					   )
					{
						if (reverse_sign)
							perm_sign[i] = '-';
						else
							perm_sign[i] = '+';
					}
					else if ( (code_table[generic_code_data::table::LABEL][inv_perm[i]] == generic_code_data::POSITIVE && code_table[generic_code_data::table::TYPE][inv_perm[i]] == generic_code_data::TYPE1) ||
					          (code_table[generic_code_data::table::LABEL][inv_perm[i]] == generic_code_data::NEGATIVE && code_table[generic_code_data::table::TYPE][inv_perm[i]] == generic_code_data::TYPE2)
							)
					{
						if (reverse_sign)
							perm_sign[i] = '+';
						else
							perm_sign[i] = '-';
					}
				}	

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "over_preferred_gauss_code:       perm_sign = ";
	for (int i=0; i< num_crossings; i++)
		debug << perm_sign[i] << ' ';
	debug << endl;
}
				
				bool minimum_weight = false;
				if (min_perm_weight.size() == 0 || (perm_weight < min_perm_weight))
				{
					minimum_weight = true;
				}				
				else if (perm_weight == min_perm_weight)
				{
					/* compare signs */
					bool equal_signs = true;
					
					for (int i=0; i< num_crossings; i++)
					{
						if (perm_sign[i] == '+' && min_perm_sign[i] == '-')
						{
							equal_signs = false;
							minimum_weight = true;
							break;
						}
						else if (perm_sign[i] == '-' && min_perm_sign[i] == '+')
						{
							equal_signs = false;
							break; // perm_sign lexicographically bigger than min_perm_sign
						}
					}
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "over_preferred_gauss_code:       perm_weight equal to min_perm_weight, perm_sign" << (minimum_weight? "" : " not") << " better than min_perm_sign ";
	for (int i=0; i< num_crossings; i++)
		debug << min_perm_sign[i] << ' ';
	debug << endl;
}
					if (equal_signs && num_components >1)
					{
						/* compare component lengths */
						if (component_length < min_component_length)
							minimum_weight = true;
							
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "over_preferred_gauss_code:       perm_weight equal to min_perm_weight and perm_sign equials min_perm_sign; component_length" << (minimum_weight? "" : " not") << " better than min_component_length ";
	for (int i=0; i< num_components; i++)
		debug << min_component_length[i] << ' ';
	debug << endl;
}
							
					}
				}

				if (minimum_weight)
				{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	if (min_perm_weight.size() != 0)
	{
		debug << "over_preferred_gauss_code:       better than current min_perm_weight = ";
		for (int i=0; i< num_terms; i++)
			debug << '(' << min_perm_weight[i].first << "," << (min_perm_weight[i].second?'U':'O') << ")" << ' ';
		debug << endl;
	}
}
					min_perm_weight = perm_weight;
					min_perm_sign = perm_sign;
					min_perm = perm;
					min_component_perm = component_perm;
					min_component_cycle = component_cycle;
					min_component_length = component_length;
					min_term_crossing = term_crossing;
					min_r_flag = r_flag;					
				}
				else
				{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	if (min_perm_weight.size() != 0)
	{
		debug << "over_preferred_gauss_code:       not better than current min_perm_weight = ";
		for (int i=0; i< num_terms; i++)
			debug << '(' << min_perm_weight[i].first << "," << (min_perm_weight[i].second?'U':'O') << ")" << ' ';
		debug << endl;
	}
}
				}
			
				if (num_components == 1 && min_component == 1)
				{
					break;
				}
				else
				{
					/* increment component_cycle lexicographically up to the number of terms in each component.	*/
					for (int k= num_components-1; k >=min_component; k--)
					{
						component_cycle[k]++;
						
						if (component_cycle[k] >= g.num_component_edges[k])
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

		/* if we're dealing with an unoriented Gauss code, consider other orientations of the components
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
	debug << "over_preferred_gauss_code: advanced r_flag to: ";
	for (int i=0; i< num_components; i++)
		debug << r_flag[i] << " ";
	debug << endl;
}
				vector<int> new_term_crossing = initial_term_crossing;
				
				for (int i=0; i< num_components; i++)
				{
					if (r_flag[i] == false)
						continue;
						
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "over_preferred_gauss_code: reversing component " << i << endl;	
						
					int first = g.first_edge_on_component[i];
					int last = g.first_edge_on_component[i]+g.num_component_edges[i]-1;
			
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "over_preferred_gauss_code:   first component edge = " << first << " last component edge = " << last << endl;	
						
					for (int j=0; j< g.num_component_edges[i]; j++)
						new_term_crossing[first+j] = initial_term_crossing[last-j];
				}
									
				term_crossing = new_term_crossing;
				
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "over_preferred_gauss_code: updated term_crossing: ";	
	for (int i=0; i< num_terms; i++)
		debug << term_crossing[i] << ' ';
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
	debug << "over_preferred_gauss_code: min_r_flag = ";
	for (int i=0; i< num_components; i++)
		debug << min_r_flag[i] << " ";
	debug << endl;
	debug << "over_preferred_gauss_code: min_component_perm = ";	
	for (int i=0; i< num_components; i++)
			debug << min_component_perm[i] << ' ';
	debug << endl;
	debug << "over_preferred_gauss_code: min_component_cycle = ";	
	for (int i=0; i< num_components; i++)
			debug << min_component_cycle[i] << ' ';
	debug << endl;
	debug << "over_preferred_gauss_code: min_perm = ";	
	for (int i=0; i< num_crossings; i++)
		debug << min_perm[i] << ' ';
	debug << endl;
	debug << "over_preferred_gauss_code: min_perm_sign = ";	
	for (int i=0; i< num_crossings; i++)
		debug << min_perm_sign[i] << ' ';
	debug << endl;
	debug << "over_preferred_gauss_code: min_perm_weight = ";
	for (int i=0; i< num_terms; i++)
		debug << '(' << min_perm_weight[i].first << "," << (min_perm_weight[i].second?'U':'O') << ")" << ' ';
	debug << endl;
	debug << "over_preferred_gauss_code: min_term_crossing = ";
	for (int i=0; i< num_terms; i++)
		debug << min_term_crossing[i] << ' ';
	debug << endl;
}

	ostringstream oss;
//	int place = 0;
	for (int k=0; k< num_components; k++)
	{
		int component = min_component_perm[k];
		int num_component_terms = g.num_component_edges[component];

		for (int i=0; i < num_component_terms; i++)
		{
			int index = (min_component_cycle[component]+i)%num_component_terms + g.first_edge_on_component[component];								
			int terminating_crossing = min_term_crossing[index];				
			bool over_arc = true;
			int old_edge = index;
			
			if (min_r_flag[component])
			{
				int reverse_component_cycle = (num_component_terms-1)-min_component_cycle[component];
				int old_offset = reverse_component_cycle-i;
				int old_cpt_edge = (old_offset+num_component_terms)%num_component_terms;
				old_edge =  old_cpt_edge + g.first_edge_on_component[component];
			}

			if (code_table[generic_code_data::table::OPEER][terminating_crossing] == old_edge)
			{	
				if (code_table[generic_code_data::table::LABEL][terminating_crossing] == generic_code_data::POSITIVE)
					over_arc = false;
			}
			else if (code_table[generic_code_data::table::EPEER][terminating_crossing] == old_edge)
			{					
				if (code_table[generic_code_data::table::LABEL][terminating_crossing] == generic_code_data::NEGATIVE)
					over_arc = false;
			}
			
			if (!over_arc)
				oss << '-';
			
			oss << min_perm[terminating_crossing]+1;
			
			if (i < num_component_terms-1)
				oss << ' ';
				
//			place++;
		}
		if (k < num_components-1)
			oss << ',';
	}
	oss << '/';

	for (int i=0; i< num_crossings; i++)
	{
		oss << min_perm_sign[i];
		if (i < num_crossings - 1)
			oss << ' ';
	}
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "over_preferred_gauss_code: over preferred Gauss code string = " << oss.str() << endl;

	return oss.str();
}

int amalgamate_zig_zag_counts(int a, int b)
{
	if (a == 0)
		return b;
	
	if (b == 0)
		return a;
		
	int result;
	int a_first_parity;
	int b_first_parity;
	int a_last_parity;
	int b_last_parity;
	
	if (a<0)
	{ 
		a_first_parity = generic_code_data::label::LEFT;
		a = abs(a);
		a_last_parity = (a%2 == 1? generic_code_data::label::LEFT : generic_code_data::label::RIGHT);
	}
	else
	{
		a_first_parity = generic_code_data::label::RIGHT;
		a_last_parity = (a%2 == 1? generic_code_data::label::RIGHT : generic_code_data::label::LEFT);
	}

	if (b<0)
	{ 
		b_first_parity = generic_code_data::label::LEFT;
		b = abs(b);
		b_last_parity = (b%2 == 1? generic_code_data::label::LEFT : generic_code_data::label::RIGHT);
	}
	else
	{
		b_first_parity = generic_code_data::label::RIGHT;
		b_last_parity = (b%2 == 1? generic_code_data::label::RIGHT : generic_code_data::label::LEFT);
	}

//cout << "a_first_parity = " << (a_first_parity == generic_code_data::label::LEFT? "LEFT":"RIGHT") << endl;
//cout << "a_last_parity = " << (a_last_parity == generic_code_data::label::LEFT? "LEFT":"RIGHT") << endl;
//cout << "b_first_parity = " << (b_first_parity == generic_code_data::label::LEFT? "LEFT":"RIGHT") << endl;
//cout << "b_last_parity = " << (b_last_parity == generic_code_data::label::LEFT? "LEFT":"RIGHT") << endl;

	if (a_last_parity == b_first_parity)
	{
		result = a-b;
//cout << "cancelling parity, difference = " << result << endl;
		if (result < 0)
		{
			if ((result % 2 == 0 && b_last_parity == generic_code_data::label::LEFT) ||
			    (result % 2 != 0 && b_last_parity == generic_code_data::label::RIGHT))
			{
//cout << "reverse negative result, result%2=" << result%2 << endl;
				result *=-1; // working backwards, so first parity is now right
			}
		}
		else
		{
			if (a_first_parity == generic_code_data::label::LEFT)
			{
				result *=-1;
//cout << "reverse positive result" << endl;
			}
		}
	}
	else
	{
		result = a+b;
//cout << "additive parity, sum = " << result << endl;
		if (a_first_parity == generic_code_data::label::LEFT)
			result *=-1;
	}
	
	return result;
}

/* An immersion diagram D of a classical knot or link that does not contain any Reidemeister I loops, is non-prime if there exist two edges 
   e1 and e2 that are both members of the same two turning cycles c1 and c2.  To see this, let v1 and v2 be points in the in the interior of
   e1 and e2 respectively and let R1 and R2 be the regions of the immersion's complement bounded by c1 and c2.  Choose the point at infinity 
   so that R1 is the infinite region.  We may join v1 and v2 by paths p1 and p2 in the regions R1 and R2 so that 
   p1 \cap D = p2 \cap D = v1 \cup v2.  Thus p1 \cup p2 is a simple closed curve that separates the crossings of D into two non-empty subsets, 
   so that e1 and e2 realize D as a connected sum.   

   A Reidemeister I loop in a semi-arc creates two edges either side of the Reidemeister I crossing in the same two turning cycles.  Hence the 
   requirement that there are no Reidemeister I loops.
   
   The shadow of a non-prime knot or link is 2-separable and hence not 3-connected.  We can use the same technique to detect doodles or flat
   knots that are 2-connected but not 3-connected, but in these cases we can also look for "crossing connected sums", where there is a simple 
   closed curve intersecting the diagram in two crossings, such that the curve meets opposite regions at each crossing.  (Note that there is 
   a simple closed curve that meets the diagram in the two separating crossings if we have an "edge connected sum" but that curve doesn't 
   meet opposite regions of the separating crossings.)
   
   Turning cycles move forwards along even edges and backwards along odd edges, so opposite regions at a crossing both belong to left or to 
   right turning cycles, not one of each, as can be seen from the diagram below

     Type 1                Type 2
        
    o \ R / o             e \ L / e
       \ /                   \ /         
    L   X  L              R   X   R
       / \                   / \
    e / R \ e             o / L \ o
   
   If flat_crossings is true three_connected looks for both edge connected sums, as above, and crossing connected sums.  For the latter, we 
   construct a matrix, cycle crossing, that records the crossings encountered by each turning cycle.  We then look for a pair of crossings 
   in the same two left or right turning cycles  that are not adjacent to each other in both turning cycles.  If they were adjacent in both 
   cycles, we would have found a bigon, so removing the two end vertices does not separate the remaining vertices and consequently doesn't 
   indicate a 2-separable diagram.   
       
   In the case of virtual knots and links, or virtual doodles,there are cases such as [-5 -7 -11 -1 -13 -3 -9]/ * # # # * # * demonstrating
   that one cannot apply the same techniques when the diagram has virtual crossings.  Instead, we test whether removing two classical 
   crossings disconnects the diagram by checking the gauss_vertex_span of their complement.  That is the span of vertices reachable via Gauss arcs.
   
   We can handle knotoids and multi-linkoids in a similar manner to virtuals but using vertex_span, based on immersion edges rather than Gauss
   arcs, since since there are no virtual crossing other than shortcut crossings that are virtual closure virtuals (for multi-linkoids) or shortcut
   crossings that are treated much like virtuals (classical knotoids).

   Virtual knotoids and multi-linkoids, that contain virtual crossings other than within shortcuts are handled using gauss_vertex_span, as in the virtual
   knot or link cases.
   
   Long knots are currently not supported by this function.
   
   Function returns
   -1 Undefined
    1 three_connected
    0 not three_connected
*/
int three_connected(generic_code_data& peer_code_data, bool flat_crossings)
{	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "three_connected: code = ";
	write_peer_code (debug, peer_code_data);
	debug << endl;
	debug << "three_connected: peer_code_data:" << endl;
	print_code_data(debug, peer_code_data, "three_connected: ");
	debug << "three_connected: flat_crossings = " << flat_crossings << endl;
}
	int num_crossings = peer_code_data.num_crossings;
	int num_edges = 2*num_crossings;
	int num_left_cycles=0;
	int num_cycles=0;

	matrix<int>& code_table = peer_code_data.code_table;

	bool shortcut_in_diagram = false;
	if (peer_code_data.immersion_character != generic_code_data::character::CLOSED && peer_code_data.shortcut_crossing.size())
		shortcut_in_diagram = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "three_connected: shortcut_in_diagram = " << shortcut_in_diagram << endl;

	int num_real_crossings = 0;
	bool virtual_diagram = false;
	
	
	for (int i=0; i< num_crossings; i++)
	{
//if (debug_control::DEBUG >= debug_control::SUMMARY)
//    debug << "three_connected: crossing label = " << peer_code_data.code_table[generic_code_data::table::LABEL][i] << endl;

		if (peer_code_data.code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::VIRTUAL && !(shortcut_in_diagram && peer_code_data.shortcut_crossing[i]))
			virtual_diagram = true; // virtual and not a shortcut (virtual closure) crossing
			
		if (peer_code_data.code_table[generic_code_data::table::LABEL][i] != generic_code_data::label::VIRTUAL && !(shortcut_in_diagram && peer_code_data.shortcut_crossing[i]))
			num_real_crossings++;
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "three_connected: virtual_diagram (i.e.  non-shortcut virtuals present) = " << virtual_diagram << ", num_real_crossings = " << num_real_crossings << endl;


	if (num_real_crossings < 4)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "three_connected: found only " << num_real_crossings << " classical and non-shortcut, crossings, returning false" << endl;
    
		return -1;
	}
	
	if (peer_code_data.immersion_character == generic_code_data::character::LONG_KNOT)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "three_connected: long knots are not suported, returning false" << endl;
    
		return -1;
	}

	/* record the real crossings */
	vector<int> real_crossing(num_real_crossings);
	int index=0;
	for (int i=0; i< num_crossings; i++)
	{				
		if (peer_code_data.code_table[generic_code_data::table::LABEL][i] != generic_code_data::label::VIRTUAL && (!shortcut_in_diagram || !peer_code_data.shortcut_crossing[i]))
			real_crossing[index++] = i;
	}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "three_connected: classical crossings:";
	for (int i=0; i< num_real_crossings; i++)
		debug << real_crossing[i] << ' ';
	debug << endl;
}
	bool real_crossings_flat = true;
	for (int i=0; i< num_real_crossings; i++)
	{
		if (peer_code_data.code_table[generic_code_data::table::LABEL][i] != generic_code_data::label::FLAT)
		{
			real_crossings_flat = false;
			break;
		}
	}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected: real_crossings_flat = " << real_crossings_flat << endl;

	int prime = 1; // i.e. 3-connected

	if (virtual_diagram == false && peer_code_data.immersion_character == generic_code_data::character::CLOSED) // clause deals with both flat and non-flat crossings
	{

		matrix<int> cycle (num_crossings+2, num_edges+1);
		calculate_turning_cycles(peer_code_data, cycle, num_left_cycles, num_cycles);
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
    debug << "three_connected: number of turning cycles = " << num_cycles;
	for (int i=0; i<num_cycles; i++)
	{
		debug << "\nthree_connected: cycle " << i << " length = " << cycle[i][0] << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " "	;
	}
	debug << endl;
}
	
		for (int i=0; i< num_cycles && prime; i++)
		{
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected: i-cycle " << i << " length = " << cycle[i][0] << endl;	
			
			vector<int> common_edge_flag(cycle[i][0]);
				
			for (int j=0; j< num_cycles && prime; j++)
			{					
				if (j==i)
					continue;
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected:   j-cycle " << j;
									
	//			for (int k=0; k< cycle[i][0]; k++)
	//				common_edge_flag[k] = 0;
	
				/* look for an edge of cycle i in cycle j */
				for (int k1=1; k1<= cycle[i][0] && prime; k1++)
				{				
					for (int l1=1; l1<= cycle[j][0] && prime; l1++)
					{
						if (cycle[j][l1] == cycle[i][k1])
						{
	//						common_edge_flag[k-1] = 1;
	//						break;
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << ", found common edge " << cycle[i][k1];
	
							/* look for a subsequent edge common to both cycles */
							for (int k2=k1+1; k2<= cycle[i][0] && prime; k2++)
							{
								for (int l2=1; l2<= cycle[j][0] && prime; l2++)
								{
									if (cycle[j][l2] == cycle[i][k2])
									{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " and common edge " << cycle[i][k2] << endl;
	
										/* we have a pair of edges in the same two cycles */
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "three_connected: found an edge connected sum between regions " << i << " and " << j << " at edges " 
	      << cycle[i][k1] << " and " << cycle[i][k2] << endl;
}
										prime = 0; // not three_connected
									}
								}
							}
						}					
					}
				}
				
				if (prime)
				{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << ", no pair of edges in common" << endl;				
				}	
			}			
		}
	
		if (prime)
		{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected: no region realizing an edge connected sum found" << endl;
		}
		
		if (prime && (real_crossings_flat || flat_crossings))
		{
	
			/* record the crossings belonging to each turning cycle */
			matrix<int> cycle_crossing (num_crossings+2, num_edges+1);
	
			for (int i=0; i< num_cycles; i++)
			{
				cycle_crossing[i][0] = cycle[i][0];
				for (int j=1; j<= cycle[i][0]; j++)
				{
					int crossing;
					if (cycle[i][j] < 0)
						crossing = peer_code_data.orig_crossing[abs(cycle[i][j])];
					else
						crossing = peer_code_data.term_crossing[cycle[i][j]];
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected: cycle edge " << cycle[i][j] << " arrives at crossing = " << crossing << endl;
					
					cycle_crossing[i][j] = crossing;
				}
			}
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
    debug << "three_connected: cycle_crossing:";
	for (int i=0; i<num_cycles; i++)
	{
		debug << "\nthree_connected: cycle " << i << " length = " << cycle_crossing[i][0] << ": ";
		for (int j=1; j<=cycle_crossing[i][0]; j++)
			debug << cycle_crossing[i][j] << " "	;
	}
	debug << endl;
}
		
			/* Look for pairs of crossings in different left, or different right turning cycles */
			for (int LR=0; LR < 2; LR++)
			{
				int cycle_search_start;
				int cycle_search_end;
				if (LR == 0)
				{
					cycle_search_start = 0;
					cycle_search_end = num_left_cycles;
				}
				else
				{
					cycle_search_start = num_left_cycles;
					cycle_search_end = num_cycles;
				}
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected: cycle_search_start " << cycle_search_start << ", cycle_search_end " << cycle_search_end << endl;
				
				for (int i=cycle_search_start; i< cycle_search_end-1 && prime; i++)
				{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected: checking cycle " << i << " for crossing-connected sum" << endl;
		
					/* we are looking for pairs of crossings in the same two cycles, so only need check up to the 
					   penultimate crossing of cycle i for the first crossing, since the last crossing will be found
					   as a second common crossing, if there are common crossings between cycle i and cycle k
					*/
					for (int j1=1; j1<= cycle_crossing[i][0]-1 && prime; j1++) // j1 will index the first common crossing we find in cycle i
					{
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected:   i-cycle first crossing " << cycle_crossing[i][j1];
		
						/* look for another left turning cycle containing cycle_crossing[i][j1] */
						bool first_found = false;
						for (int k=i+1; k< cycle_search_end && prime; k++)
						{
							for (int l1=1; l1 <= cycle_crossing[k][0] && prime; l1++) // l1 will index the first common crossing we find in cycle k
							{
								if (cycle_crossing[k][l1] == cycle_crossing[i][j1])
								{
									first_found = true;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " found in cycle " << k << endl;
		
									/* We have one common crossing between cycles i and k, so now look for a second */
									for (int j2=j1+1; j2 <= cycle_crossing[i][0] && prime; j2++) // j2 will index the second common crossing we find in cycle i
									{
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected:   i-cycle second crossing " << cycle_crossing[i][j2];
			
										bool second_found = false;								
										for (int l2=1; l2 <= cycle_crossing[k][0] && prime; l2++) // l2 will index the second common crossing we find in cycle k
										{
											if (cycle_crossing[k][l2] == cycle_crossing[i][j2])
											{
												second_found = true;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " also found in cycle " << k << endl;
												
												/* we have two common crossings in cycles i and k, so check whether they are adjacent crossings in both
												   cycles.  If they are, then these crossings are part of a bigon and do not separate the diagram, so 
												   do not indicate the diagram is non-prime
												*/
												bool adjacent_i = false;
												if (abs(j2-j1) == 1 || (j2 == cycle[i][0] && j1 == 1))
													adjacent_i = true;
		
												bool adjacent_k = false;
												if (abs(l2-l1) == 1 || (l2 == cycle[k][0] && l1 == 1))
													adjacent_k = true;
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected:   adjacent_i = " << adjacent_i << ", adjacent_k = " << adjacent_k << endl;
													
												if (!adjacent_i || !adjacent_k)
												{
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "three_connected: found a crossing connected sum between regions " << i << " and " << k << " at crossings " 
	      << cycle_crossing[i][j1] << " and " << cycle_crossing[i][j2] << endl;
}
													
													prime = 0; // not three_connected
												}
												else
												{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected:   crossings do not recognize a connected sum" << endl;
												}
											}
										}
										
										if (!second_found)
										{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " not found in cycle " << k << endl;
										}
									}
								}
							}
						}
						
						if (!first_found)
						{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " not found in any other left cycle " << endl;
						}
					}
				}
			} // end of LR loop
				
			if (prime)
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected: no region realizing a crossing connected sum found" << endl;
			}
		}

	}
	else if (real_crossings_flat || flat_crossings)
	{
		/* This case includes classical knotoids and multi-knotoids. We consider all combinations to check whether removing 
		   two real (classical, non-shortcut) crossings separates the diagram.
		*/
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected: pairwise crossing test with vertex_span required" << endl;


		/* if we're dealing with a virtual diagram, prepare the gauss_crossing_map and gauss_arc_map */
		generic_code_data gauss_code_data;
		vector<int> gauss_crossing_map; 
		vector<int> gauss_arc_map; 
		if (virtual_diagram)
		{
			ostringstream oss;
			write_gauss_code(oss,peer_code_data);
			read_gauss_code(gauss_code_data,oss.str());
			
			identify_gauss_maps(peer_code_data, gauss_crossing_map, gauss_arc_map);  // maps Gauss crossings to immersion crossings

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
    debug << "three_connected: gauss_crossing_map: ";
    for (int i=0; i< gauss_code_data.num_crossings; i++)
		debug << gauss_crossing_map[i] << ' ';
	debug << endl;
    debug << "three_connected: gauss_arc_map: ";
    for (int i=0; i< 2*peer_code_data.num_crossings; i++)
		debug << gauss_arc_map[i] << ' ';
	debug << endl;
}		
		}

		/* identify and count the shortcut edges in the diagram */
		vector<int> shortcut_edge(num_edges);
		int num_shortcut_edges = 0;
		
		if (peer_code_data.immersion_character != generic_code_data::character::CLOSED)
		{
			if (peer_code_data.immersion_character == generic_code_data::character::KNOTOID)
			{
				shortcut_edge[0] = 1;
				num_shortcut_edges = 1;
			}
			else if (peer_code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
			{
				int semi_arc=-1;  // the first immersion semi-arc of the shortcut
				int head = peer_code_data.head;
								
				if (code_table[generic_code_data::table::LABEL][head] == generic_code_data::POSITIVE)
					semi_arc = code_table[generic_code_data::table::OPEER][head];
				else if (code_table[generic_code_data::table::LABEL][head] == generic_code_data::NEGATIVE)
					semi_arc = 2*head;
				
				for (int i=semi_arc; i< peer_code_data.num_component_edges[0]; i++)
				{
					shortcut_edge[i] = 1;
					num_shortcut_edges++;
				}
				shortcut_edge[0] = 1; // edge 0 contains the leg of the knotoid
				num_shortcut_edges++;				
			}
			else // (peer_code_data.immersion_character == generic_code_data::character::MULTI_LINKOID)
			{
				for (size_t i=0; i< peer_code_data.component_type.size(); i++)  // safer than int i=0 to num_components
				{
					if (peer_code_data.component_type[i].type == component_character::PURE_START_LEG || 
					    peer_code_data.component_type[i].type == component_character::PURE_END_LEG)
					{
						for (int j=peer_code_data.component_type[i].head_semi_arc; j< peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]; j++)
						{
							shortcut_edge[j] = 1;
							num_shortcut_edges++;
						}
						
						if (peer_code_data.component_type[i].type == component_character::PURE_START_LEG) // leg on first, even edge
						shortcut_edge[peer_code_data.first_edge_on_component[i]] = 1; // first edge contains the leg of this open component
						num_shortcut_edges++;				

					}
					else if (peer_code_data.component_type[i].type == component_character::KNOT_TYPE_START_LEG) // leg and hed on even edge
					{
						shortcut_edge[peer_code_data.first_edge_on_component[i]] = 1;
						num_shortcut_edges += 1;
					}
					else if (peer_code_data.component_type[i].type == component_character::KNOT_TYPE_END_LEG) // leg and hed on odd edge
					{
						shortcut_edge[peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-1] = 1; // last edge on component
						num_shortcut_edges += 1;
					}
					
				}
			}
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "three_connected: shortcut_edge: ";
	for (int i=0; i< num_edges; i++)
		debug << shortcut_edge[i] << ' ';
	debug << "\nthree_connected: num_shortcut_edges = " << num_shortcut_edges << endl;

}
		}

		vector<int> exclude(num_shortcut_edges+8);

		/* add the shortcut edges to the end of exclude; if the diagram has virtual crossings other than in the virtual closure of a
		   multi-linkoid, use gauss_arc_map to map the edges to Gauss arcs 
		*/
		int index = 8;
		for (int i=0; i< num_edges; i++)
		{
			if (shortcut_edge[i] == 1)
				exclude[index++] = (virtual_diagram?gauss_arc_map[i]:i);
		}
				
		/* for each pair of classical crossings remove the edges incident with those crossings and see if the remainder of the diagram is disconnected */
		for (int i=0; i< num_real_crossings-1 && prime; i++) // there are at least 4
		{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected: testing classical crossing " << real_crossing[i] << endl;
	
			exclude[0] = code_table[generic_code_data::table::EVEN_TERMINATING][real_crossing[i]];
			exclude[1] = code_table[generic_code_data::table::ODD_ORIGINATING][real_crossing[i]];
			exclude[2] = code_table[generic_code_data::table::ODD_TERMINATING][real_crossing[i]];
			exclude[3] = code_table[generic_code_data::table::EVEN_ORIGINATING][real_crossing[i]];
			
			if (virtual_diagram)
			{
				exclude[0] = gauss_arc_map[exclude[0]];
				exclude[1] = gauss_arc_map[exclude[1]];
				exclude[2] = gauss_arc_map[exclude[2]];
				exclude[3] = gauss_arc_map[exclude[3]];
			}
	
			for (int j=i+1; j< num_real_crossings && prime; j++)
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected:   crossing " << real_crossing[j] << endl;
	
				exclude[4] = code_table[generic_code_data::table::EVEN_TERMINATING][real_crossing[j]];
				exclude[5] = code_table[generic_code_data::table::ODD_ORIGINATING][real_crossing[j]];
				exclude[6] = code_table[generic_code_data::table::ODD_TERMINATING][real_crossing[j]];
				exclude[7] = code_table[generic_code_data::table::EVEN_ORIGINATING][real_crossing[j]];

				if (virtual_diagram)
				{
					exclude[4] = gauss_arc_map[exclude[4]];
					exclude[5] = gauss_arc_map[exclude[5]];
					exclude[6] = gauss_arc_map[exclude[6]];
					exclude[7] = gauss_arc_map[exclude[7]];
				}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "three_connected:   exclude " << (virtual_diagram? "Gauss arcs: ": "edges: ");
	for (int k=0; k< 8; k++)
		debug << exclude[k] << ' ';
	debug << endl;
}
				list<int> span;
				for (int k=0; k< num_real_crossings; k++)
				{
					if (k!=i && k !=j)
					{						
						if (virtual_diagram)
						{
							int gauss_initial_crossing = -1;
							for (int i=0; i < gauss_code_data.num_crossings; i++)
							{
								if (gauss_crossing_map[i] == real_crossing[k])
								{
									gauss_initial_crossing = i;
									break;
								}
							}
							
							if (gauss_initial_crossing == -1)
							{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "three_connected: Error! three_connected could not map inital_crossing onto a Gauss crossing" << endl;
				    
								cerr << "Error! three_connected could not map inital_crossing onto a Gauss crossing" << endl;
								exit(0);
							}
							
							span = vertex_span (gauss_code_data, gauss_initial_crossing, &exclude);

						}
						else
						{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected:   start vertex span test at crossing " << real_crossing[k] << endl;
	
							span = vertex_span (peer_code_data, real_crossing[k], &exclude);
						}
						break;
					}
				}
						
				list<int>::iterator lptr = span.begin();
				int real_span = 0;
						
				while (lptr != span.end())
				{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected:     vertex " << (virtual_diagram?gauss_crossing_map[*lptr]:*lptr);
							
					if (find(real_crossing.begin(),real_crossing.end(),(virtual_diagram?gauss_crossing_map[*lptr]:*lptr)) != real_crossing.end())
					{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " is a real crossing" << endl;
						real_span++;
					}
					else
					{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " is not a real crossing" << endl;
					}
					lptr++;							
				}
						
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected:   real_span = " << real_span << endl;
						
				if (real_span != num_real_crossings-2)
				{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected: classical crossings " << real_crossing[i] << " and " << real_crossing[j] << " 2-separate" << endl;
					prime = 0;  // diagram is 2-separable
					break;
				}
			}
			
			if (prime)
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "three_connected: classical crossings " << real_crossing[i] << " does not 2-separate with any other classical crossing" << endl;
			}	
		}	
	}
	else
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "three_connected: input has crossings that are not flat and is not a CLOSED classical diagram, returning undefined" << endl;
    
		prime = -1;
	}
	
	return prime;
}

/* valid_multi_linkoid_input checks that all the crossings encountered on pure knotoid open components after the head_semi_arc
   are virtual and that the virtual closure contains no self intersections.
   
   As the check is carried out, the function creates and sets the shortcut crossings, by analogy with valid_knotoid_input.  
   Here, we simply set shortcut_crossing[i] to be 1 if the ith crossing is a shortcut crossing. 
   the opposite direction.  
*/
bool valid_multi_linkoid_input(generic_code_data& peer_code_data)
{
	if (peer_code_data.immersion_character != generic_code_data::character::MULTI_LINKOID)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "valid_multi_linkoid_input: presented with code data that is not a multi-linkoid, returning false" << endl;
    
		return false;
	}
	
	vector<int> shortcut_crossing(peer_code_data.num_crossings);
	
	for (int component = 0; component < peer_code_data.num_components; component++)
	{
		if (peer_code_data.component_type[component].type == component_character::PURE_START_LEG || peer_code_data.component_type[component].type == component_character::PURE_END_LEG)
		{
			int start = peer_code_data.component_type[component].head_semi_arc;
			int end = peer_code_data.first_edge_on_component[component]+peer_code_data.num_component_edges[component]-1;
			
			if (peer_code_data.component_type[component].type == component_character::PURE_END_LEG)
				end--;  // last edge of the component is the leg of the knotoid
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_multi_linkoid_input: component " << component << " shortcut start = " << start << ", end = " << end << endl;

			for (int edge=start; edge <= end; edge++)
			{
				int crossing = peer_code_data.term_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_multi_linkoid_input: edge " << edge << ", crossing " << crossing;
				
				int peer;
				if (edge%2)
					peer = peer_code_data.code_table[generic_code_data::table::EVEN_TERMINATING][crossing];
				else
					peer = peer_code_data.code_table[generic_code_data::table::ODD_TERMINATING][crossing];
				
				if (edge < peer && peer <= end)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", peer << " << peer << " also lies in the shortcut, invalid multi-linkoid" << endl;
					return false;
				}

				if (peer_code_data.code_table[generic_code_data::table::LABEL][crossing] == generic_code_data::label::VIRTUAL)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", crossing is virtual, setting shortcut_crossing[" << crossing <<"] = 1" << endl;
					shortcut_crossing[crossing] = 1;
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", crossing is not virtual, invalid multi-linkoid" << endl;
					return false;
				}
			}
		}
	}

	peer_code_data.shortcut_crossing = shortcut_crossing;
	return true;
}

int  code_data_writhe(generic_code_data& code_data)
{
	int writhe=0;
	
	for (int i=0; i< code_data.num_crossings; i++)
	{
		int label = code_data.code_table[generic_code_data::table::LABEL][i];
		int type = code_data.code_table[generic_code_data::table::TYPE][i];
		
		if ((label == generic_code_data::label::POSITIVE && type == generic_code_data::type::TYPE2) || (label == generic_code_data::label::NEGATIVE && type == generic_code_data::type::TYPE1))
		{
			writhe++;
		}
		else if ((label == generic_code_data::label::POSITIVE && type == generic_code_data::type::TYPE1) || (label == generic_code_data::label::NEGATIVE && type == generic_code_data::type::TYPE2))
		{
			writhe--;
		}
	}
	return writhe;
}

/* add_Reidemeister_1_loop adds a Reidemeister loop of the specified crossing type in the first edge of the last component of the diagram 
   described by code_data.  This edge is always even, so the new crossing is then num_crossings and the new odd peer 2*num_crossings+1.
   
   By default, the function adds a positive turn; that is anticlockwise in the plane, so , so the type of the crossing is type 1.
   Then for a positive crossing the label is '-' and for a negative crossing the label is '+'.
*/   
generic_code_data add_Reidemeister_1_loop (generic_code_data& code_data, bool positive_crossing, bool positive_turn = true)
{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "add_Reidemeister_1_loop: presented with code ";
	write_code_data(debug,code_data);
	debug << endl;
}

	if (code_data.type != generic_code_data::code_type::peer_code)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "add_Reidemeister_1_loop: presented with code_data that is not a peer-code, returning code_data" << endl;
	
		return code_data;
	}
	else
	{
		ostringstream oss;
		write_peer_code(oss,code_data);
		string old_code = oss.str();
		string::size_type pos = old_code.find(']');
		int new_peer = 2*code_data.num_crossings+1;
		
		string new_code = old_code.substr(0,pos);
		
		if (positive_turn)  // anticlockwise turn, so type 1 crossing
		{
			new_code += " -" + util::itos(new_peer);
			if (positive_crossing)
				new_code += old_code.substr(pos,string::npos) + " -";
			else
				new_code += old_code.substr(pos,string::npos) + " +";
		}
		else // clockwise turn, so type 2 crossing
		{
			new_code += " " + util::itos(new_peer);
			
			if (positive_crossing)
				new_code += old_code.substr(pos,string::npos) + " +";
			else
				new_code += old_code.substr(pos,string::npos) + " -";
		}
				
		generic_code_data new_code_data;
		read_peer_code(new_code_data,new_code);
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "add_Reidemeister_1_loop: produced code ";
	write_code_data(debug,new_code_data);
	debug << endl;
}
		return new_code_data;
	}
}

