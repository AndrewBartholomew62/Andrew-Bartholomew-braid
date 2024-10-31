/**************************************************************************

void hamiltonian_circuit(generic_code_data& code_data)

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
#include <matrix.h>
#include <debug-control.h>
#include <generic-code.h>
#include <gauss-orientation.h>

bool gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data, bool optimal=true, vector<int>* gauss_crossing_perm=0, bool evaluate_gauss_crossing_perm=false);

list<vector<int> > hamiltonian_circuit(generic_code_data& code_data, bool list_all_circuits, bool count_circuits_only, bool edge_circuit, int include_edge)
{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "hamiltonian_circuit: provided with code data ";
	write_code_data(debug,code_data);
	debug << endl;
	print_code_data(debug,code_data,"hamiltonian_circuit: ");
}

	list<vector<int> > circuit_list;

	generic_code_data peer_code_data;
	if(code_data.type == generic_code_data::gauss_code)
	{
		if (gauss_to_peer_code(code_data, peer_code_data))
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "hamiltonian_circuit: converted Gauss code to peer code ";
	write_code_data(debug,peer_code_data);
	debug << endl;
	print_code_data(debug,peer_code_data,"hamiltonian_circuit: ");
}
		}
		else
		{		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "generic_code: Error converting Gauss code to peer code, doing nothing" << endl;
	
			return circuit_list;
		}
	}
	else
	{
		peer_code_data = code_data;
	}
		
//	bool hamiltonian_found = false;	
	int num_crossings = peer_code_data.num_crossings;
	int num_components = code_data.num_components;
	int num_edges = 2*num_crossings;
	matrix<int>& code_table = peer_code_data.code_table;
	vector<int>& first_edge_on_component = peer_code_data.first_edge_on_component;
	vector<int>& num_component_edges = peer_code_data.num_component_edges;

	bool virtual_crossing = false;
	for (int i=0; i< num_crossings; i++)
	{
		if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::VIRTUAL)
		{
			virtual_crossing = true;
			break;
		}
	}
	
	if (virtual_crossing)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "hamiltonian_circuit: Hamiltonian circuits not defined for diagrams including virtual crossings, doing nothing" << endl;
	
		return circuit_list;
	}

	int num_cycles;
	int num_left_cycles;

	matrix<int> cycle(num_crossings+2, num_edges+1);
	calculate_turning_cycles(peer_code_data, cycle, num_left_cycles, num_cycles);

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "hamiltonian_circuit: num_cycles = " << num_cycles << endl;
    debug << "hamiltonian_circuit: number of left turning cycles = " << num_left_cycles;
	for (int i=0; i<num_left_cycles; i++)
	{
		debug << "\nhamiltonian_circuit: cycle " << i << " length = " << cycle[i][0] << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " ";
	}
	debug << endl;
    debug << "hamiltonian_circuit: number of right turning cycles = " << num_cycles - num_left_cycles;
	for (int i=num_left_cycles; i<num_cycles; i++)
	{
		debug << "\nhamiltonian_circuit: cycle " << i << " length = " << cycle[i][0] << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " "	;
	}
	debug << endl;
}

	/* create a crossing_matrix that records the regions, indexed by the row in cycle of the corresponding tuning cycle, that appear
	   to the left when moving anti-clockwise around the crossing from the naming edge
	*/
	matrix<int> crossing_region(num_crossings,4);
	
	for (int i=0; i< num_crossings; i++)
	{
		/* look for the naming edge in the left cycles */
		int region=0;
		bool found = false;
		for (int j=0; j<num_left_cycles && !found; j++)
		{
			for (int k=1; k<=cycle[j][0] && !found; k++) 
			{
				if (abs(cycle[j][k]) == 2*i)
				{
					region = j;
					found = true;
				}
			}
		}
		
		crossing_region[i][0] = region;
		
		/* look for the naming edge in the right cycles */
		found = false;
		for (int j=num_left_cycles; j<num_cycles && !found; j++)
		{
			for (int k=1; k<=cycle[j][0] && !found; k++) 
			{
				if (abs(cycle[j][k]) == 2*i)
				{
					region = j;
					found = true;
				}
			}
		}
		
		crossing_region[i][1] = region;

		/* look for the naming edge successor in the left cycles */
		found = false;
		for (int j=0; j<num_left_cycles && !found; j++)
		{
			for (int k=1; k<=cycle[j][0] && !found; k++) 
			{
				if (abs(cycle[j][k]) == 2*i+1)
				{
					region = j;
					found = true;
				}
			}
		}
		
		crossing_region[i][2] = region;
		
		/* look for the naming edge successor in the right cycles */
		found = false;
		for (int j=num_left_cycles; j<num_cycles && !found; j++)
		{
			for (int k=1; k<=cycle[j][0] && !found; k++) 
			{
				if (abs(cycle[j][k]) == 2*i+1)
				{
					region = j;
					found = true;
				}
			}
		}
		
		crossing_region[i][3] = region;
		
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "hamiltonian_circuit: crossing_region: " << endl;
	print(crossing_region, debug, 3, "hamiltonian_circuit: ");
}

	/* note the max and min region at each crossing */
	vector<int> max_region(num_crossings);
	vector<int> min_region(num_crossings);

	for (int i=0; i< num_crossings; i++)
	{
		max_region[i] = min_region[i] = crossing_region[i][0];
		for (int j=1; j<4; j++)
		{
			if (crossing_region[i][j] > max_region[i])
				max_region[i] = crossing_region[i][j];			
			else if (crossing_region[i][j] < min_region[i])
				min_region[i] = crossing_region[i][j];			
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{	
	debug << "hamiltonian_circuit: max_region: ";
	for (int i=0; i< num_crossings; i++)
		debug << max_region[i] << ' ';
	debug << endl;

	debug << "hamiltonian_circuit: min_region: ";
	for (int i=0; i< num_crossings; i++)
		debug << min_region[i] << ' ';
	debug << endl;
}		

	/* create a place_map matrix of integers that holds the (r,c) indices of all instances of a region within crossing_region
	   as consecutive entries in each row.
	   num_boundary_vertices will hold the number of valid entries: that is, the number of crossings in the boundary
	   of each region
	*/
	matrix<int> place_map(num_cycles,2*num_crossings);
	vector<int> num_boundary_vertices(num_cycles);
	
	for (int i=0; i< num_crossings; i++)
	{
		for (int j=0; j< 4; j++)
		{
			int region = crossing_region[i][j];
			place_map[region][2*num_boundary_vertices[region]] = i;
			place_map[region][2*num_boundary_vertices[region]+1] = j;
			num_boundary_vertices[region]++;
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{	
	debug << "hamiltonian_circuit: place_map: " << endl;
	for (int i=0; i< num_cycles; i++)
	{
		debug << "hamiltonian_circuit:   " << num_boundary_vertices[i] << ": ";
		for (int j=0; j< num_boundary_vertices[i]; j++)
			debug << '(' << place_map[i][2*j] << ',' << place_map[i][2*j+1] << ") ";
		debug << endl;
	}
}
	
	/* We record the colour associated with each region in a vector that we increment as a binary number RED = 0, BLUE = 1 */
	int RED = 0;
	int BLUE= 1;

	vector<int> colour_map(num_cycles,RED);
			
	bool not_finished = true;
	do
	{
		bool contains_include_edge = true; 
		
		if (include_edge != -1)
			contains_include_edge = false;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{	
	debug << "hamiltonian_circuit: colour_map: ";
	for (int i=0; i< num_cycles; i++)
		debug << (colour_map[i] == RED? "R ": "B ");
	debug << endl;
}		
		/* record the colour of the regions around each crossing when coloured with colour_map */
		matrix<int> crossing_colour(num_crossings,4);
		for (int i=0; i< num_cycles; i++)
		{
			for (int j=0; j< num_boundary_vertices[i]; j++)
				crossing_colour[place_map[i][2*j]][place_map[i][2*j+1]] = colour_map[i];				
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{	
	debug << "hamiltonian_circuit: crossing_colour: " << endl;
	for (int i=0; i< num_crossings; i++)
	{
		debug << "hamiltonian_circuit:   ";
		for (int j=0; j< 4; j++)
			debug << (crossing_colour[i][j] == RED? "R ": "B ");
		debug << endl;
	}
}

		/* reject any colour_map that doesn't result in at least one region of each colour at each crossing 
		   if we sum the region colours at each crossing we should not get a total of zero or four.
		*/
		bool valid_colour_map = true;
		
		for (int i=0; i< num_crossings; i++)
		{
			int colour_sum = 0;
			for (int j=0; j< 4; j++)
				colour_sum += crossing_colour[i][j];
				
			if (colour_sum == 0) // all RED
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "hamiltonian_circuit: all regions around crossing " << i << " coloured RED." << endl;

				valid_colour_map = false;				
				colour_map[max_region[i]] = BLUE;
				
				for (int j=max_region[i]+1; j< num_cycles; j++)
					colour_map[j] = RED;
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
{	
	debug << "hamiltonian_circuit: incremented colour_map: ";
	for (int i=0; i< num_cycles; i++)
		debug << (colour_map[i] == RED? "R ": "B ");
	debug << endl;
}
				break;
			}
			else if (colour_sum == 4) // all BLUE
			{				
if (debug_control::DEBUG >= debug_control::SUMMARY)
{	
	debug << "hamiltonian_circuit: all regions around crossing " << i << " coloured BLUE" << endl;
}
				valid_colour_map = false;
				not_finished = false;
//				for (int j= min_region[i]; j>=0; j--)
				for (int j= max_region[i]; j>=0; j--)
				{
					if (colour_map[j] == RED)
					{
						colour_map[j] = BLUE;

						for (int k=j+1; k< num_cycles; k++)
							colour_map[k] = RED;
						not_finished = true;
						break;
					}
				}
if (debug_control::DEBUG >= debug_control::SUMMARY)
{	
	if (not_finished)
	{
		debug << "hamiltonian_circuit: incremented colour_map: ";
		for (int i=0; i< num_cycles; i++)
			debug << (colour_map[i] == RED? "R ": "B ");
		debug << endl;
	}
	else
	{
		debug << "hamiltonian_circuit: cannot increment colour_map further, finished" << endl;
	}
}
				break;
			}
		}
		
		if (valid_colour_map)
		{
			/* check we don't have opposite regions of the same colour at any crossing */
			for (int i=0; i< num_crossings && valid_colour_map; i++)
			{
				if (crossing_colour[i][0] == crossing_colour[i][2] && crossing_colour[i][1] == crossing_colour[i][3])
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "hamiltonian_circuit: opposite regions at crossing " << i << " have the same colour" << endl;

					valid_colour_map = false;
					
					/* increment the region_colour, if possible */
					not_finished = false;
					for (int i = num_cycles-1; i>=0; i--)
					{
						if (colour_map[i] == RED)
						{
							colour_map[i] = BLUE;
							
							for (int j=i+1; j< num_cycles; j++)
								colour_map[j] = RED;
							not_finished = true;
							break;
						}
					}

					if (not_finished)
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{	
	debug << "hamiltonian_circuit: incremented colour_map: ";
	for (int i=0; i< num_cycles; i++)
		debug << (colour_map[i] == RED? "R ": "B ");
	debug << endl;
}
					}
					else
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
		debug << "hamiltonian_circuit: cannot increment colour_map further, finished" << endl;
					}
					
					break;
				}
			}
		}
		
		if (valid_colour_map)
		{
			vector<int> circuit(num_crossings);
			circuit[0] = 0;
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
		debug << "hamiltonian_circuit: trace boundary of colour_map starting at crossing zero" << endl;
		
			/* find an edge incident with crossing zero that lies in the colour boundary and determine whether we are going forwards
			   or backwards with respect to the label orientation when we leave the crossing on the selected edge.
			*/
			int edge = 0;
			bool forwards;

			int peer = code_table[generic_code_data::table::OPEER][0];
			int peer_component = 0;
			for (int i=1; i< num_components; i++)
			{
				if (peer >= first_edge_on_component[i])
					peer_component++;
				else
					break;
			}			
			int peer_successor = (peer+1 - first_edge_on_component[peer_component])%
		                     num_component_edges[peer_component] + first_edge_on_component[peer_component];

			for (int i=0; i<4; i++)
			{
				if (crossing_colour[0][i] != crossing_colour[0][(i+1)%4])
				{
					if (code_table[generic_code_data::table::TYPE][0] == generic_code_data::TYPE1)
					{
						switch (i)
						{
							case 0: edge = 0; forwards = false; break;
							case 1: edge = peer_successor; forwards = true; break;
							case 2: edge = 1; forwards = true; break;
							case 3: edge = peer; forwards = false; break;
						}
					}
					else
					{
						switch (i)
						{
							case 0: edge = 0; forwards = false; break;
							case 1: edge = peer; forwards = false; break;
							case 2: edge = 1; forwards = true; break;
							case 3: edge = peer_successor; forwards = true; break;
						}
					}
					break;
				}
			}
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
		debug << "hamiltonian_circuit:   found boundary edge " << edge << " at crossing 0, forwards = " << forwards << endl;

			if (edge_circuit)
				circuit[0] = edge;
			
			if (include_edge != -1 && edge == include_edge)
				contains_include_edge = true;
			
			int edge_count = 0;						
			int next_crossing = 0;
			do
			{				
				edge_count++;
				
				/* identify the crossing we reach going along edge */
				next_crossing = (forwards? peer_code_data.term_crossing[edge] : peer_code_data.orig_crossing[edge]);
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
		debug << "hamiltonian_circuit:   next_crossing = " << next_crossing;
		debug << (code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1? " TYPE1" : " TYPE2") << endl;
}

				if (next_crossing !=0)
				{			

					if (!edge_circuit)
						circuit[edge_count] = next_crossing;					

					peer = code_table[generic_code_data::table::OPEER][next_crossing];
					peer_component = 0;
					for (int i=1; i< num_components; i++)
					{
						if (peer >= first_edge_on_component[i])
							peer_component++;
						else
							break;
					}			
					peer_successor = (peer+1 - first_edge_on_component[peer_component])%
				                     num_component_edges[peer_component] + first_edge_on_component[peer_component];
					
					/* identify the unique other boundary edge at next_crossing */
					int boundary_edge;
					for (int i=0; i<4; i++)
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
		debug << "hamiltonian_circuit:     edge " << i << " boundary colours "
		      << (crossing_colour[next_crossing][i] == RED? 'R': 'B') << " and " << (crossing_colour[next_crossing][(i+1)%4] == RED? 'R': 'B') << endl;
}

						if (crossing_colour[next_crossing][i] != crossing_colour[next_crossing][(i+1)%4])
						{
							if (code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1)
							{
								switch (i)
								{
									case 0: boundary_edge = 2*next_crossing; forwards = false; break;
									case 1: boundary_edge = peer_successor; forwards = true; break;
									case 2: boundary_edge = 2*next_crossing+1; forwards = true; break;
									case 3: boundary_edge = peer; forwards = false; break;
								}
							}
							else
							{
								switch (i)
								{
									case 0: boundary_edge = 2*next_crossing; forwards = false; break;
									case 1: boundary_edge = peer; forwards = false; break;
									case 2: boundary_edge = 2*next_crossing+1; forwards = true; break;
									case 3: boundary_edge = peer_successor; forwards = true; break;
								}
							}

if (debug_control::DEBUG >= debug_control::SUMMARY)
		debug << "hamiltonian_circuit:   found boundary edge " << boundary_edge << endl;
							
							if (boundary_edge != edge)
							{
								edge = boundary_edge;
if (debug_control::DEBUG >= debug_control::SUMMARY)
		debug << "hamiltonian_circuit:   next boundary edge " << edge << " at crossing " << next_crossing << ", forwards = " << forwards << endl;
								break;
							}
						}
					}
					
					if (edge_circuit)
						circuit[edge_count] = edge;

					if (include_edge != -1 && edge == include_edge)
						contains_include_edge = true;
					
				}
			} while (next_crossing !=0); // not back where we started
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "hamiltonian_circuit: returned to crossing zero with edge_count " << edge_count << endl;
		
			if (edge_count == num_crossings && contains_include_edge) 
			{
				circuit_list.push_back(circuit);				
				if (!list_all_circuits && !count_circuits_only)
					break;				
			}
		}
		
		/* increment the region_colour, if possible */
		if (valid_colour_map)
		{
			not_finished = false;
			for (int i = num_cycles-1; i>=0; i--)
			{
				if (colour_map[i] == RED)
				{
					colour_map[i] = BLUE;
					
					for (int j=i+1; j< num_cycles; j++)
						colour_map[j] = RED;
					not_finished = true;
					break;
				}
			}
		}
		
	} while (not_finished);
	
	/* remove the duplicates from circuit_list, we shall find circuits twice by interchanging RED and BLUE
	   but we shall always traverse circuits in the same direction, regardless of colour polarity.
	*/
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "hamiltonian_circuit: initial circuit_list contains " << circuit_list.size() << " entries" << endl;
	
	if (circuit_list.size() > 0)
	{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "hamiltonian_circuit: remove circuit_list duplicates" << endl;

		list<vector<int> >::iterator aptr = circuit_list.begin();
		
		while (aptr != circuit_list.end())
		{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "hamiltonian_circuit:  aptr: " ;
	for (size_t i=0; i< aptr->size(); i++)
		debug << (*aptr)[i] << ' ';
	debug << endl;
}

			list<vector<int> >::iterator bptr = aptr;
			bptr++;
			
			while (bptr != circuit_list.end())
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "hamiltonian_circuit:  bptr: " ;
	for (size_t i=0; i< bptr->size(); i++)
		debug << (*bptr)[i] << ' ';
}
				if (*bptr == *aptr)
				{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << ":  duplicate found" << endl;
	
					bptr = circuit_list.erase(bptr);
				}
				else	
				{			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << endl;

					bptr++;
				}
			}
			
			aptr++;
		}
	}
	return circuit_list;
}
