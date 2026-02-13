/********************************************************************************************************

The algorithm as implemented by JG and here is not optimal in all cases.  For example, in the diagram
[7 11 9 13 3 5 1]/- * * - + + + one can perform a Reidemeister III move by moving edge 13 over crossing 0
thereby creating a Reidemeister I loop at edge zero, producing immersion [-13 11 9 -1 3 5 -7]/- * * + + + + 
The Gauss code of this diagram is -1 -2 -3 -4 2 -5 3 4 5 1/+ - + + - but gauss_to_peer produces 
[-15 -9 -11 -1 -13 -5 -3 -7]/- - - + * * * + with one more virtual crossing than necessary.  The reason for the
additional crosing is that the algorithm has to choose between equal best choices: the next crossing, or (in the above case)
which anchor edge to use when rearranging the fringe ready for adding a crossing and without looking forward it is 
impossible to know whether one or other of the immediate choices is better.  The only way to produce an optimal diagram 
would be to consider all such equal best choices at each step, which is likely to extend the running time unreasonably.
                                                                         
void add_virtual_crossing(int location, vector<int>& fringe, vector<int>& num_virtual_crossings_on_gauss_arc, list<gc_pc_xlabels>& virtual_crossings) 
void assign_gauss_arc_direction(int arc_1, int arc_2, matrix<int>& gauss_code_table, int crossing, vector<int>& gauss_arc_direction)
string direction_to_string(int edge, vector<int>& gauss_arc_direction)
void remove_fringe_edge(int edge, vector<int>& current_fringe)
bool gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data)
bool classical_gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data)
********************************************************************************************************/
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <valarray>
#include <list>
#include <iomanip>
#include <ctype.h>

using namespace std;

extern ofstream     debug;
extern ofstream     output;
extern ifstream     input;

#include <util.h>
#include <matrix.h>
#include <generic-code.h>
#include <debug-control.h>
#include <gauss-orientation.h>


/********************* Function prototypes ***********************/

void add_virtual_crossing(int location, vector<int>& fringe, vector<int>& num_virtual_crossings_on_gauss_arc, list<gc_pc_xlabels>& virtual_crossings) 
{	
    
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "add_virtual_crossing: adding virtual crossing " << virtual_crossings.size()+1 << ": location = " << location << ", gauss edges " << fringe[location] << ',' << fringe[location+1] << endl;
	
//if (debug_control::DEBUG >= debug_control::SUMMARY)
//	debug << location << ',' << (fringe[location]-1+num_gauss_arcs)%num_gauss_arcs << ',' << (fringe[location+1]-1+num_gauss_arcs)%num_gauss_arcs << endl;

    gc_pc_xlabels xlabels;
    xlabels.gauss.first = fringe[location];
    xlabels.gauss.second = fringe[location+1];
    
    num_virtual_crossings_on_gauss_arc[xlabels.gauss.first]++;
    num_virtual_crossings_on_gauss_arc[xlabels.gauss.second]++;
    
    virtual_crossings.push_back(xlabels);
    
    swap (fringe[location],fringe[location+1]);   
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

/* arc_1 and arc_2 are both incident at the supplied Gauss crossing, we set the direction of the terminating arc to UPWARDS and
   that of the originating arc DOWNWARDS
*/
void assign_gauss_arc_direction(int arc_1, int arc_2, matrix<int>& gauss_code_table, int crossing, vector<int>& gauss_arc_direction)
{
	
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "assign_gauss_arc_direction: crossing = " << crossing << ", arc_1 = " << arc_1 << ", arc_2 = " << arc_2 << endl;
	debug << "assign_gauss_arc_direction: gauss_code_table[generic_code_data::table::ODD_TERMINATING][crossing] = " << gauss_code_table[generic_code_data::table::ODD_TERMINATING][crossing] 
	      << " gauss_code_table[generic_code_data::table::EVEN_TERMINATING][crossing] = " << gauss_code_table[generic_code_data::table::EVEN_TERMINATING][crossing] << endl;
}

		if ( (gauss_code_table[generic_code_data::table::ODD_TERMINATING][crossing] == arc_1 && gauss_code_table[generic_code_data::table::EVEN_ORIGINATING][crossing] == arc_2) 
		    || (gauss_code_table[generic_code_data::table::EVEN_TERMINATING][crossing] == arc_1 && gauss_code_table[generic_code_data::table::ODD_ORIGINATING][crossing] == arc_2))
		{		
			gauss_arc_direction[arc_1] = gc_pc_xlabels::direction::UPWARDS;
			gauss_arc_direction[arc_2] = gc_pc_xlabels::direction::DOWNWARDS;
		}
		else
		{
			gauss_arc_direction[arc_1] = gc_pc_xlabels::direction::DOWNWARDS;
			gauss_arc_direction[arc_2] = gc_pc_xlabels::direction::UPWARDS;
		}
		
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "assign_gauss_arc_direction: gauss_arc_direction[arc_1] = " << direction_to_string(arc_1,gauss_arc_direction)
	      << " gauss_arc_direction[arc_2] = " << direction_to_string(arc_2,gauss_arc_direction) << endl;
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
bool gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data, bool optimal, vector<int>* gauss_crossing_map, bool evaluate_gauss_crossing_perm)
{

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code: initial gauss_code_data :" << endl;
	print_code_data(debug,gauss_code_data,"gauss_to_peer_code: ");	
}

	int num_components = gauss_code_data.num_components;
	int num_gauss_crossings = gauss_code_data.num_crossings;
	int num_gauss_arcs = 2*num_gauss_crossings;
	matrix<int>& gauss_code_table = gauss_code_data.code_table;

	for (int i=0; i< num_gauss_crossings; i++)
	{
	    if (gauss_code_table[generic_code_data::table::LABEL][i] != generic_code_data::POSITIVE  && gauss_code_table[generic_code_data::table::LABEL][i] != generic_code_data::NEGATIVE && gauss_code_table[generic_code_data::table::LABEL][i] != generic_code_data::FLAT)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "gauss_to_peer_code: gauss_code_data includes an unsupported crossing label" << endl;
	print_code_data(debug,gauss_code_data,"gauss_to_peer_code: ");	
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
		if (gauss_code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE)
		{
			PD_data[i][0] = gauss_code_table[generic_code_data::table::ODD_TERMINATING][i];

			if (gauss_code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1)
			{
				PD_data[i][1] = gauss_code_table[generic_code_data::table::EVEN_TERMINATING][i];
				PD_data[i][2] = gauss_code_table[generic_code_data::table::EVEN_ORIGINATING][i];
				PD_data[i][3] = gauss_code_table[generic_code_data::table::ODD_ORIGINATING][i];
			}
			else
			{
				PD_data[i][1] = gauss_code_table[generic_code_data::table::ODD_ORIGINATING][i];
				PD_data[i][2] = gauss_code_table[generic_code_data::table::EVEN_ORIGINATING][i];
				PD_data[i][3] = gauss_code_table[generic_code_data::table::EVEN_TERMINATING][i];
			}			
		}
		else // generic_code_data::NEGAITIVE or generic_code_data::FLAT
		{
			PD_data[i][0] = gauss_code_table[generic_code_data::table::EVEN_TERMINATING][i];

			if (gauss_code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1)
			{
				PD_data[i][1] = gauss_code_table[generic_code_data::table::EVEN_ORIGINATING][i];
				PD_data[i][2] = gauss_code_table[generic_code_data::table::ODD_ORIGINATING][i];
				PD_data[i][3] = gauss_code_table[generic_code_data::table::ODD_TERMINATING][i];
			}
			else
			{
				PD_data[i][1] = gauss_code_table[generic_code_data::table::ODD_TERMINATING][i];
				PD_data[i][2] = gauss_code_table[generic_code_data::table::ODD_ORIGINATING][i];
				PD_data[i][3] = gauss_code_table[generic_code_data::table::EVEN_ORIGINATING][i];
			}			
		}
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
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


if (debug_control::DEBUG >= debug_control::DETAIL)
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
		vector<int> fringe (max(4,2*num_gauss_crossings)); // need at least four fringe edges, even for codes like 1,-1/+
		
		
		
		for (int i=0; i< 4; i++)
			fringe[i] = PD_data[initial_crossing][i];
			
		int num_fringe_edges = 4;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
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
		{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "gauss_to_peer_code: assigning directions: i = " << i << ", fringe[0+i] = " << fringe[0+i] << " fringe[2+i] = " << fringe[2+i] << endl;
	
			assign_gauss_arc_direction(fringe[0+i],fringe[2+i],gauss_code_table,initial_crossing,candidate_gauss_arc_direction);
		}
			
		/* check whether the first crossing includes a Reidemeister I loop or a linked component containing
	       one classical crossing.
		*/
		bool crossing_loop_detected = false;
		for (int i=0; i < 3 && !crossing_loop_detected; i++)
		{
			for (int j=i+1; j < 4; j++)
			{
				
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "gauss_to_peer_code: checking fringe i = " << i << " fringe[i] = " << fringe[i] << " j = " << j << " fringe[j] = " << fringe[j] << endl;

				if (fringe[i] == fringe[j])
				{	
					candidate_gauss_arc_direction[fringe[i]] = gc_pc_xlabels::direction::SIDEWAYS;
					
					if ((i-j+4)%4 == 2)
					{
	                    add_virtual_crossing(i, fringe, candidate_num_virtual_crossings_on_gauss_arc, candidate_virtual_crossings);
						remove_fringe_edge(fringe[j],fringe); // use j because we've just added the virtual crossing at i.
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: linked component containing a single classical crossing detected at first crossing" << endl;
					}
					else
					{
						remove_fringe_edge(fringe[i],fringe);
					
if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
		
			continue;
		}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code: num_fringe_edges = " << num_fringe_edges << " current fringe: ";
	for (int f=0; f< num_fringe_edges; f++)
		debug << fringe[f] << ' ';
	debug << endl;
}

if (debug_control::DEBUG >= debug_control::DETAIL)
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
			int best_match_count = 0;
			int best_precedence = 5;  // bigger than the maximum prececence we use
			int last_fringe_index;
			int best_last_fringe_index;
			int first_fringe_index;
			int first_PD_index;
			int last_PD_index;
			int crossing_edge_c; // used for diametrically_connecting_edge_pair
			int crossing_edge_d;
			bool diametrically_connecting_edge_pair = false;
		
if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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

/*							
if (match_count == 1 && PD_data[i][last_PD_index] == 0 && gauss_code_data.immersion_character == generic_code_data::character::KNOTOID)
{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   don't want to match directly onto leg of knotoid, precedence incremented to 4" << endl;
	
	precedence = 4;
}
*/
					/*  Are the matching edges contiguous in the fringe and in the correct order?  If they are, then moving from right 
					    to left along the fringe the edges must appear in the corresponding order when moving anti-clockwise around 
					    the crossing.  That is, from the last_fringe_index moving backwards and the last_PD_index moving forwards.  
					    If this is not the case, increment the precedence to 1,
					*/												
		
					for (int j = 0; j < match_count; j++)
					{
						if (PD_data[i][(last_PD_index + j) % 4] != fringe[last_fringe_index - j])
						{
						
if (debug_control::DEBUG >= debug_control::DETAIL)
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
						
if (debug_control::DEBUG >= debug_control::DETAIL)
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
					    
if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   wrapping edges does not produce a contiguous set of matching edges in the correct order, precedence incremented to 2" << endl;
	
										precedence = 2; // no wrap-around possible, so virtual crossing needed
									}
							    }
							    
							    if (precedence == 1)
							    {
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   contiguous matching edges in the correct order found wrapping to end of fringe" << endl;
								}					
						    } 
						    else 
						    {
if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   matching edges not even contiguous in the fringe, precedence incremented to 3" << endl;
								
								precedence = 3; // other edges must be moved out of the way
							    break;
							}
						}
						
						if (precedence == 2)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   matching edges are adjacent in the fringe, though not in the correct order" << endl;
						}
					}
				
					if (precedence < best_precedence || (precedence == best_precedence && match_count > best_match_count)) 
//					if (match_count > best_match_count || (match_count == best_match_count && precedence < best_precedence)) 
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
		
		
if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
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
		      		      

if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   evauate wrap_flags for matching_fringe_index " << i << ", edge " << fringe[matching_fringe_index[i]] << endl;;

					int count = 0;
					for (int j = 0; j < best_match_count; j++) 
					{
						int direct = abs(matching_fringe_index[j] - matching_fringe_index[i]);  // distance from i to j within the fringe
						int wrapped = num_fringe_edges - direct; // distance wrapping around the end of the fringe
						
						count += (direct < wrapped ? direct : wrapped);  // count the smallest number of edges that would need to be moved

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code:     fringe index " << matching_fringe_index[j] << ", edge " << fringe[matching_fringe_index[j]]
	      << " direct count = " << direct << ", wrapped count = " << wrapped << ", cumuative count = " << count << endl;
}

						if (wrapped == direct)
						{
							/* if the matching fringe edges are already in the correct order, we prefer the direct approach, if not, wrapping
							   will introduce one fewer virtual crossing since we will not need to re-order the two edges */
							   
							int PD_index_i;
							int PD_index_j;
							for (int k = 0; k < 4; k++)
							{
								if (fringe[matching_fringe_index[i]] == PD_data[best_crossing][k]) 
									PD_index_i = k;
								else if (fringe[matching_fringe_index[j]] == PD_data[best_crossing][k]) 
									PD_index_j = k;								
							}
							
							/* We only consider the cases where the two matching fringe edges are adjacent, since if they are diametrically opposite we have 
							   the special case, which is handled separately.
							   
							   to determine whether we prefer the direct approach:
							   if matching_fringe_index[i] is less than matching_fringe_index[j] we look for fringe[j] adjacent clockwise from fringe[i]
							   if matching_fringe_index[i] is greater than matching_fringe_index[j] we look for fringe[j] adjacent anti-clockwise from fringe[i]
							*/
							if ((matching_fringe_index[i] < matching_fringe_index[j] && PD_index_i == (PD_index_j+1)%4) ||
							    (matching_fringe_index[j] < matching_fringe_index[i] && PD_index_j == (PD_index_i+1)%4) )
							{
								matching_fringe_index_wrap_flag[j] = false; // prefer the direct route
								
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code:     prefer the direct approach: fringe index " << matching_fringe_index[i] << " is " 
	      << (matching_fringe_index[i] < matching_fringe_index[j]? "smaller":"greater") 
	      << " than fringe index " << matching_fringe_index[j]  << endl;
	debug << "gauss_to_peer_code:     fringe[" << fringe[matching_fringe_index[i]] << "] lies at PD index " << PD_index_i 
	      << ", fringe[" << fringe[matching_fringe_index[j]] << "] lies at PD index " << PD_index_j  << endl;
}

							}
							else
							{
								matching_fringe_index_wrap_flag[j] = true;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code:     prefer the wrapped approach: fringe index " << matching_fringe_index[i] << " is " 
	      << (matching_fringe_index[i] < matching_fringe_index[j]? "smaller":"greater") 
	      << " than fringe index " << matching_fringe_index[j]  << endl;
	debug << "gauss_to_peer_code:     fringe edge " << fringe[matching_fringe_index[i]] << " lies at PD index " << PD_index_i 
	      << ", fringe edge " << fringe[matching_fringe_index[j]] << " lies at PD index " << PD_index_j  << endl;
}
							}
						}
						else
						{
							matching_fringe_index_wrap_flag[j] = (wrapped < direct); 
						}
					}
				
if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     setting fringe_anchor_edge_index to " << fringe_anchor_edge_index << endl;
					}
			    }

if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
				
					too_many_virtual_crossings = true;
					break;
				}
		
if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   check to the right of the anchor edge = " << fringe[fringe_anchor_edge_index] << endl;
				
					int right_count=0; // counts the number of matching edges we need to wrap on the right of the anchor edge
					
				    for (int i = best_match_count -1; i > fringe_anchor_edge_index; i--)
				    {
						if (best_matching_fringe_index_wrap_flag[i])
						{
							right_count++;
						
/*if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
					
						too_many_virtual_crossings = true;
						break;
					}
							
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   number of matching edges that need wrapping to the right of the anchor edge = " << right_count << endl;
				
				    if (right_count > 0) 
				    {	  			 
						/* all the matching edges that need wrapping to the left are now flush-right in the fringe, 
						   so we need to wrap the last right_count edges from the right of the fringe to the left */

if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
				
					too_many_virtual_crossings = true;
					break;
				}
					
if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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
				   best_precedence would be 0.  If there are fewer than 4 matching edges, there is a unique order of the matching edges in the fringe 
				   that permits the non-matching edges of the crossing to become part of the new fringe without introducing additional virtual crossings.  
				   If best_match_count == 4, any order in the fringe would work but we choose the one minimizing the number of virtual crossings required.
				   
				   We identify the PD_data index of the first matching edge in the fringe in first_PD_index.  Since the PD_data records edges
				   anticlockwise around the crossing, this edge is the last of the matching edges around the crossing when moving anticlockwise.
				*/
	
				first_PD_index=-1;  //initialized to catch errors
				diametrically_connecting_edge_pair = false;
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: precedence == 2, the " << best_match_count << " matching edges in the fringe are contiguous but are in the wrong order" << endl;
		  
				if (best_match_count < 4)  // AB be = 2 or 3
				{
					/* identify the indices in PD_data of the edges that connect to the fringe as we move through the fringe from left to right */
					vector<int> PD_index_of_connecting_edge(best_match_count,-1); 				
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   last_fringe_index =  " << last_fringe_index << ", best_match_count = " << best_match_count << endl;
	
					for (int i = 0; i < best_match_count; i++)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
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
				
if (debug_control::DEBUG >= debug_control::DETAIL)
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
				  				  
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: matching edges are diametrically opposite on the crossing, first index of matching edges = " << first_PD_index << endl;
	
						crossing_edge_c = PD_data[best_crossing][(PD_index_of_connecting_edge[0]+1)%4];
						crossing_edge_d = PD_data[best_crossing][(PD_index_of_connecting_edge[1]+1)%4];  // as in above diagram
					
if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
						
							too_many_virtual_crossings = true;
							break;
						}
								    
if (debug_control::DEBUG >= debug_control::DETAIL)
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
					
if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: starting index " << i << ", edge " << PD_data[best_crossing][i] << ", requires adding " << virtual_crossing_count << " virtual crossings" << endl;
					  
						if (virtual_crossing_count < least_number_of_virtual_crossings) 
						{
							best_starting_PD_index = i;
							least_number_of_virtual_crossings = virtual_crossing_count;
						}
					}
					
					first_PD_index = best_starting_PD_index; 

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: first matching edge in the fringe has PD_data index " << first_PD_index << endl;
								
				}
				
				if (!diametrically_connecting_edge_pair)
				{				   
					/* We want to add the crossing to the current fringe by considering the edges clockwise around the crossing from 
					   first_PD_index forwards along the fringe from the first fringe index.
					*/		
					first_fringe_index = last_fringe_index - best_match_count + 1;
				
if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: choice of initial_crossing " << initial_crossing << " introduces more than " << min_num_virtual_crossings << " virtual crossings, jump to next initial_crossing" << endl;
					
						too_many_virtual_crossings = true;
						break;
					}
			
if (debug_control::DEBUG >= debug_control::DETAIL)
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
			
//if (debug_control::DEBUG >= debug_control::DETAIL)
//	debug << "gauss_to_peer_code: check for crossing index " << j << ", edge " << PD_data[best_crossing][j] << " in the fringe" << endl;   
			
				for (int k=0; k< num_fringe_edges; k++)
				{
					if (PD_data[best_crossing][j] == fringe[k])
					{
//if (debug_control::DEBUG >= debug_control::DETAIL)
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
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: ready to add crossing to the fringe: first_fringe_index = " << first_fringe_index << " first_PD_index = " << first_PD_index << endl;   
	    
			int new_num_fringe_edges = num_fringe_edges + 4 - 2* best_match_count;
	    
if (debug_control::DEBUG >= debug_control::DETAIL)
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
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: loop component containing a single classical crossing detected" << endl;
	
				}
				else
				{
					assign_gauss_arc_direction(crossing_edge_d,crossing_edge_c,gauss_code_table, best_crossing,candidate_gauss_arc_direction);		
				}
			
if (debug_control::DEBUG >= debug_control::DETAIL)
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
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: Reidemeister I loop detected at crossing" << endl;
	
					}
					else if (arc_2 == arc_3)
					{
						candidate_gauss_arc_direction[arc_1] = candidate_gauss_arc_direction[fringe[first_fringe_index]];
						candidate_gauss_arc_direction[arc_2] = gc_pc_xlabels::direction::SIDEWAYS;
						
						remove_fringe_edge(arc_2, new_fringe);
						new_num_fringe_edges -=2;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: Reidemeister I loop detected at crossing" << endl;
	
					}
					else if (arc_1 == arc_3)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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
							
if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code: two new fringe edges: direction of edges " << new_fringe_edge_1 << ',' << new_fringe_edge_2 
	      << " are " << direction_to_string(new_fringe_edge_1,candidate_gauss_arc_direction) << direction_to_string(new_fringe_edge_2,candidate_gauss_arc_direction) << endl;
}
				}
				else if (best_match_count == 3)
				{
					int up_count = 0;
//					int down_count = 0;
	
			        for (int i=0; i< 3; i++)
			        {
						if (candidate_gauss_arc_direction[fringe[first_fringe_index+i]] == gc_pc_xlabels::direction::UPWARDS)
							up_count++;
//						else 
//							down_count++;
					}
					
					int new_fringe_edge = PD_data[best_crossing][(first_PD_index+1)%4];
					
					if (up_count == 2)
						candidate_gauss_arc_direction[new_fringe_edge] = gc_pc_xlabels::direction::UPWARDS;
					else
						candidate_gauss_arc_direction[new_fringe_edge] = gc_pc_xlabels::direction::DOWNWARDS;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code: one new fringe edge: direction of edge " << new_fringe_edge << " is " << direction_to_string(new_fringe_edge,candidate_gauss_arc_direction) << endl;
}
				}
							
				fringe = new_fringe;
				num_fringe_edges = new_num_fringe_edges;
			}	

if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: initial crossing " << initial_crossing << " requires minimal number of virtual crossings" << endl;

			optimal_initial_crossing = initial_crossing;
			min_num_virtual_crossings = candidate_virtual_crossings.size();
			gauss_arc_direction = candidate_gauss_arc_direction;
			virtual_crossings = candidate_virtual_crossings;
			diametric_virtual_crossings = candidate_diametric_virtual_crossings;
			num_virtual_crossings_on_gauss_arc = candidate_num_virtual_crossings_on_gauss_arc;
		}
		
		if (!optimal)
			break;
			
	} /* initial_crossing loop ends here */

if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
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
	

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: calculate immersion edge labels" << endl;
	
	for (int i=0; i< num_gauss_arcs; i++)
	{
		if (i == gauss_code_data.first_edge_on_component[component]+gauss_code_data.num_component_edges[component])
		{
			component++;
			first_immersion_edge_on_component[component] = next_immersion_arc;
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   found start of component " << component << ", first immersion edge =  " << first_immersion_edge_on_component[component] << endl;
	
		}
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   gauss_arc " << i << ", contains " << num_virtual_crossings_on_gauss_arc[i] << " virtual crossings" << endl;

		if (num_virtual_crossings_on_gauss_arc[i] !=0)
		{		
			if (gauss_arc_direction[i] == gc_pc_xlabels::direction::DOWNWARDS)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     DOWNWARDS arc" << endl;
	
				list<gc_pc_xlabels>::iterator lptr = virtual_crossings.begin();
				
				while (lptr != virtual_crossings.end())
				{
					if (lptr->gauss.first == i)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code:     found arc in first place at virtual crossing " << lptr->gauss.first << ',' << lptr->gauss.second
	      << ", allocating immersion edge label " << next_immersion_arc << endl;
}	
						lptr->immersion.first = next_immersion_arc++;
					}
					else if (lptr->gauss.second == i)
					{

if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     UPWARDS or SIDEWAYS arc" << endl;
	
				list<gc_pc_xlabels>::reverse_iterator lptr = virtual_crossings.rbegin();
				
				while (lptr != virtual_crossings.rend())
				{
					if (lptr->gauss.first == i)
					{						
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code:     found arc in first place at virtual crossing " << lptr->gauss.first << ',' << lptr->gauss.second
	      << ", allocating immersion edge label " << next_immersion_arc << endl;
}
						lptr->immersion.first = next_immersion_arc++;
					}
					else if (lptr->gauss.second == i)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
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
			if(gauss_code_table[generic_code_data::table::EVEN_TERMINATING][j] == i)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   found arc as generic_code_data::table::EVEN_TERMINATING at Gauss crossing " << j << ", allocating immersion edge label " << next_immersion_arc << endl;
	
				immersion_crossing_peers[j][0] = next_immersion_arc++;
				break;
			}
			else if(gauss_code_table[generic_code_data::table::ODD_TERMINATING][j] == i)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   found arc as generic_code_data::table::ODD_TERMINATING at Gauss crossing " << j << ", allocating immersion edge label " << next_immersion_arc << endl;
	
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

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code: immersion_crossing_peers: " << endl;
	print(immersion_crossing_peers, debug, 4, "gauss_to_peer_code: ");
}	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code: virtual_gauss_peers: " << endl;
	print(virtual_gauss_peers, debug, 4, "gauss_to_peer_code: ");
}	
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: next_immersion_arc  = " << next_immersion_arc << " (should be twice the number of immersion crossings)" << endl;
	
	/* evaluate the number of immersion edges on each component */
	vector<int> num_immersion_component_edges(num_components);
	for (int i=0; i<num_components; i++)
	{
		int next_first_edge = (i< num_components - 1? first_immersion_edge_on_component[i+1]: next_immersion_arc);
		num_immersion_component_edges[i] = next_first_edge - first_immersion_edge_on_component[i];
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
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
	vector<bool> component_cycled(num_components); // initializes to false
	
	if (num_components > 1)
	{
		list<int> next_component;
		next_component.push_back(0);
		vector <bool> component_considered(num_components);
		
		bool complete = false; // only needed for debugging
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: check components for aligned odd and even parity edges" << endl;
		do
		{
			if (next_component.size() == 0)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   stopped processing components due to empty list" << endl;
				complete = true;
				break;
			}
				
			int component = next_component.front();
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   front of next_component = " << component << endl;
		
			next_component.pop_front();
			if (component_considered[component])
				continue;
			else
				component_considered[component] = true;
			
if (debug_control::DEBUG >= debug_control::DETAIL)
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
						peer_component = 0;
						for (int j=1; j< num_components; j++)
						{
							if (peer_edge >= first_immersion_edge_on_component[j])
								peer_component++;
							else
								break;
						}
						next_component.push_back(peer_component);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:       peer component " << peer_component << endl;

						if (peer_edge % 2 == edge % 2)
						{
							/* both the first and second visits to this crossing have the same parity, so cycle the
							   edges in the component containing peer_edge to align the parity with edge */
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     immersion_crossing_peers shows the peer of " << edge << " is " << peer_edge 
	      << ", cycle component containing " << peer_edge << endl;
	                        cycle_gauss_code_labels(immersion_crossing_peers,num_immersion_component_edges,first_immersion_edge_on_component,num_immersion_crossings,num_components,peer_edge);
	                        component_cycled[peer_component] = true;
						}
																						
						break;
					}
					else if (immersion_crossing_peers[j][1] == edge)
					{
						peer_edge = immersion_crossing_peers[j][0];
						peer_component = 0;
						for (int j=1; j< num_components; j++)
						{
							if (peer_edge >= first_immersion_edge_on_component[j])
								peer_component++;
							else
								break;
						}
						next_component.push_back(peer_component);
						
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:       peer component " << peer_component << endl;

						if (peer_edge % 2 == edge % 2)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     immersion_crossing_peers shows the peer of " << edge << " is " << peer_edge
	      << ", cycle component containing " << peer_edge << endl;
		
	                        cycle_gauss_code_labels(immersion_crossing_peers,num_immersion_component_edges,first_immersion_edge_on_component,num_immersion_crossings,num_components,peer_edge);
	                        component_cycled[peer_component] = true;
						}
												
						break;
					}
				}
			}    
		} while (!complete);
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code: final immersion_crossing_peers: " << endl;
	print(immersion_crossing_peers, debug, 4, "gauss_to_peer_code: ");
}	

	matrix<int> peer_code_table(generic_code_data::table::CODE_TABLE_SIZE,num_immersion_crossings);

	/* write the generic_code_data::table::COMPONENT data into a code table for the peer code */
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
		peer_code_table[generic_code_data::table::COMPONENT][i] = component;
	}
	
	/* write the ODD and EVEN peers to peer_code_table */
	int even_edge;
	int odd_edge;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: write immersion labels:" << endl;
	
	for (int i=0; i< num_immersion_crossings; i++)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
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
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " even immersion edge = " << even_edge << ", odd immersion edge = " << odd_edge << endl;

		peer_code_table[generic_code_data::table::OPEER][even_edge/2] = odd_edge;
		peer_code_table[generic_code_data::table::EPEER][(odd_edge-1)/2] = even_edge;
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code: imersion EPEERs: ";
	for (int i=0; i< num_immersion_crossings; i++)
		debug << peer_code_table[generic_code_data::table::EPEER][i] << ' ';
	debug << endl;
	debug << "gauss_to_peer_code: imersion OPEERs: ";
	for (int i=0; i< num_immersion_crossings; i++)
		debug << peer_code_table[generic_code_data::table::OPEER][i] << ' ';
	debug << endl;
}				

	/* write the crossing generic_code_data::table::TYPE data to peer_code_table. In a gauss code, the type of crossing is assigned based on the relative 
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
if (debug_control::DEBUG >= debug_control::DETAIL)
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
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", even_peer = " << even_peer << ", even_peer_in_place_zero (of immersion_crossing_peers) = " << even_peer_in_place_zero << endl;

		if (i < num_gauss_crossings)
		{
			if (even_peer_in_place_zero)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to second visit to Gauss crossing " << i;
	
				peer_code_table[generic_code_data::table::TYPE][even_peer/2] = gauss_code_table[generic_code_data::table::TYPE][i];
			}
			else // even label corresponds to first visit to Gauss crossing
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to first visit to Gauss crossing " << i;

				if (gauss_code_table[generic_code_data::table::TYPE][i] == generic_code_data::type::TYPE1)
					peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE2;
				else 
					peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE1;
			}
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " crossing is (peer) " << (peer_code_table[generic_code_data::table::TYPE][even_peer/2] == generic_code_data::type::TYPE1? "TYPE1": "TYPE2") << endl;
		}
		else
		{
			int first_fringe_edge_direction = gauss_arc_direction[virtual_gauss_peers[i][0]];
			int second_fringe_edge_direction = gauss_arc_direction[virtual_gauss_peers[i][1]];

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code:     first_fringe_edge_direction " ;
	debug << direction_to_string(virtual_gauss_peers[i][0],gauss_arc_direction);
	debug << ", second_fringe_edge_direction " ;
	debug << direction_to_string(virtual_gauss_peers[i][1],gauss_arc_direction);
	debug << endl;
}					

			if (first_fringe_edge_direction == gc_pc_xlabels::direction::SIDEWAYS)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
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
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     sideways_arc " << virtual_gauss_peers[i][0] << " meets gauss_crossing " << gauss_crossing << endl;
	
				/* identify whether the sideways arc was the first or second visit to gauss_crossing */
				bool sideways_second_visit = (virtual_gauss_peers[i][0] > virtual_gauss_peers[i][1]);
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     sideways_second_visit to gauss_crossing " << gauss_crossing << " " << sideways_second_visit << endl;
				
				int gauss_crossing_type = gauss_code_table[generic_code_data::table::TYPE][gauss_crossing];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     gauss_crossing " << gauss_crossing << " is (Gauss) " << (gauss_crossing_type == generic_code_data::TYPE1? "TYPE1": "TYPE2") << endl;
	
				/* identify the immersion labels on the non-sideways arc.  For the virtual crossing, the label on the sideways arc is always stored
				   in position zero.  For Gauss crossings, positions 0 and 1 store the labels corresponding to the second and first visit respectively.
				*/
				int gauss_crossing_non_sideways_immersion_label = (sideways_second_visit?immersion_crossing_peers[gauss_crossing][1]:immersion_crossing_peers[gauss_crossing][0]);				
				int virtual_crossing_non_sideways_immersion_label = immersion_crossing_peers[i][1];

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code:     gauss_crossing_non_sideways_immersion_label = " << gauss_crossing_non_sideways_immersion_label << endl;
	debug << "gauss_to_peer_code:     virtual_crossing_non_sideways_immersion_label = " << virtual_crossing_non_sideways_immersion_label << endl;
}

				/* identify the non-sideways component */
				int non_sideways_component_even_edge = gauss_crossing_non_sideways_immersion_label;
				if (non_sideways_component_even_edge %2)
					non_sideways_component_even_edge = virtual_crossing_non_sideways_immersion_label;
				int non_sideways_component = peer_code_table[generic_code_data::table::COMPONENT][non_sideways_component_even_edge/2];

if (debug_control::DEBUG >= debug_control::DETAIL)
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
				    
if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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
						peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE1;				
					else
						peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE2;				
				}
				else
			    {
					if (even_peer_in_place_zero)
						peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE2;				
					else
						peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE1;				
				}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     immersion crossing " << even_peer/2 << " is " << (peer_code_table[generic_code_data::table::TYPE][even_peer/2] == generic_code_data::type::TYPE1? "TYPE1": "TYPE2") << endl;
				
			}
			else if (diametric_virtual_crossing[i-num_gauss_crossings])
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     virtual crossing " << i-num_gauss_crossings << " added as a result of a diametrically opposite matching Gauss crossing" << endl;
	
				/* If the virtual crossing was introduced as a result of added a Gauss crossing with diametrically opposite matching edges, then it involves the first
				   matching fringe edge of the Gauss crossing and crossing_edge_d, as in the above diagram.  However, when we stored the immersion_crossing_peers and
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
                   
						even_peer_in_place_zero => generic_code_data::table::TYPE II   !even_peer_in_place_zero => generic_code_data::table::TYPE I
				         
				   first_fringe_edge_direction is upwards and second_fringe_edge_direction is downwards; or
                   first_fringe_edge_direction is downwards and second_fringe_edge_direction is upwards
                   
						even_peer_in_place_zero => generic_code_data::table::TYPE I   !even_peer_in_place_zero => generic_code_data::table::TYPE II
				   				
				*/
				if ((first_fringe_edge_direction == gc_pc_xlabels::direction::UPWARDS && second_fringe_edge_direction == gc_pc_xlabels::direction::UPWARDS) ||
				    (first_fringe_edge_direction == gc_pc_xlabels::direction::DOWNWARDS && second_fringe_edge_direction == gc_pc_xlabels::direction::DOWNWARDS)
			       )
			    {
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     both Gauss arcs have the same direction at virtual crossing " << i-num_gauss_crossings << endl;
	
					if (even_peer_in_place_zero)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to first matching fringe Gauss arc,";
						peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE2;				
					}
					else
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to second matching fringe Gauss arc,";
	
						peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE1;				
					}
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     Gauss arcs have opposite direction at virtual crossing " << i-num_gauss_crossings << endl;

					if (even_peer_in_place_zero)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to first matching fringe Gauss arc,";
						peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE1;				
					}
					else
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to second matching fringe Gauss arc,";
						peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE2;				
					}
				}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " crossing is (peer) " << (peer_code_table[generic_code_data::table::TYPE][even_peer/2] == generic_code_data::type::TYPE1? "TYPE1": "TYPE2") << endl;
			}
			else if ((first_fringe_edge_direction == gc_pc_xlabels::direction::DOWNWARDS && second_fringe_edge_direction == gc_pc_xlabels::direction::DOWNWARDS) ||
			    (first_fringe_edge_direction == gc_pc_xlabels::direction::UPWARDS && second_fringe_edge_direction == gc_pc_xlabels::direction::UPWARDS)
			   )
			{				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     both Gauss arcs have the same direction at virtual crossing " << i-num_gauss_crossings << endl;
				if (even_peer_in_place_zero)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to left fringe Gauss arc of the virtual crossing " << i-num_gauss_crossings;
	
					peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE1;
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to right fringe Gauss arc of the virtual crossing " << i-num_gauss_crossings;
	
					peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE2;
				}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " crossing is (peer) " << (peer_code_table[generic_code_data::table::TYPE][even_peer/2] == generic_code_data::type::TYPE1? "TYPE1": "TYPE2") << endl;
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     Gauss arcs have opposite direction at virtual crossing " << i-num_gauss_crossings << endl;
				if (even_peer_in_place_zero)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to left fringe Gauss arc of the virtual crossing " << i-num_gauss_crossings;
	
					peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE2;
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:     even_peer assigned to right fringe Gauss arc of the virtual crossing " << i-num_gauss_crossings;
	
					peer_code_table[generic_code_data::table::TYPE][even_peer/2] = generic_code_data::type::TYPE1;
				}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " crossing is (peer) " << (peer_code_table[generic_code_data::table::TYPE][even_peer/2] == generic_code_data::type::TYPE1? "TYPE1": "TYPE2") << endl;
			}
		}
	}
		
	/* write the crossing generic_code_data::table::LABEL data to peer_code_table */
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: write crossing labels:" << endl;
	
	for (int i = 0; i< num_immersion_crossings; i++)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
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
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " even_peer = " << even_peer << endl;

		if (i < num_gauss_crossings)
		{
			if (even_peer_in_place_zero)  // even label corresponds to second visit to Gauss crossing
			{
				if (gauss_code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::POSITIVE)  //second visit to Gauss crossing goes over
					peer_code_table[generic_code_data::table::LABEL][even_peer/2] = generic_code_data::label::POSITIVE;
				else if (gauss_code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::NEGATIVE)  //second visit to Gauss crossing goes under
					peer_code_table[generic_code_data::table::LABEL][even_peer/2] = generic_code_data::label::NEGATIVE;
				else if (gauss_code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::FLAT)
					peer_code_table[generic_code_data::table::LABEL][even_peer/2] = generic_code_data::label::FLAT;
			}
			else // even label corresponds to first visit to Gauss crossing
			{
				if (gauss_code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::POSITIVE)  //second visit to Gauss crossing goes over
					peer_code_table[generic_code_data::table::LABEL][even_peer/2] = generic_code_data::label::NEGATIVE;
				else if (gauss_code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::NEGATIVE)  //second visit to Gauss crossing goes under
					peer_code_table[generic_code_data::table::LABEL][even_peer/2] = generic_code_data::label::POSITIVE;
				else if (gauss_code_table[generic_code_data::table::LABEL][i] == generic_code_data::label::FLAT)
					peer_code_table[generic_code_data::table::LABEL][even_peer/2] = generic_code_data::label::FLAT;
			}
		}
		else
		{
			peer_code_table[generic_code_data::table::LABEL][even_peer/2] = generic_code_data::label::VIRTUAL;
		}
	}
		
	/* write the originating and terminating vertices and edges */
	vector<int> term_crossing(2*num_immersion_crossings);
	vector<int> orig_crossing(2*num_immersion_crossings);
	
	for (int i=0; i< num_immersion_crossings; i++)
	{
		term_crossing[2*i] = i;
		orig_crossing[2*i+1] = i;
		term_crossing[peer_code_table[generic_code_data::table::OPEER][i]] = i;
		
		int component = peer_code_table[generic_code_data::table::COMPONENT][(peer_code_table[generic_code_data::table::OPEER][i]-1)/2];
		int peer_successor = (peer_code_table[generic_code_data::table::OPEER][i]+1 - first_immersion_edge_on_component[component])%
		                     num_immersion_component_edges[component] + first_immersion_edge_on_component[component];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: crossing " << i << " peer_successor = " << peer_successor << endl;

		orig_crossing[peer_successor] = i;
		
		peer_code_table[generic_code_data::table::EVEN_TERMINATING][i] = 2*i;
		peer_code_table[generic_code_data::table::ODD_ORIGINATING][i] = 2*i+1;
		peer_code_table[generic_code_data::table::ODD_TERMINATING][i] = peer_code_table[generic_code_data::table::OPEER][i];
		peer_code_table[generic_code_data::table::EVEN_ORIGINATING][i] = peer_successor;
	}
	
	peer_code_data.type = generic_code_data::peer_code;
	peer_code_data.num_open_components = gauss_code_data.num_open_components;
	peer_code_data.component_type = vector<component_character>(num_components);
	peer_code_data.head = -1;
	peer_code_data.head_zig_zag_count = gauss_code_data.head_zig_zag_count;
	peer_code_data.immersion_character = gauss_code_data.immersion_character;
	peer_code_data.num_crossings = num_immersion_crossings;
	peer_code_data.num_components = num_components;
	peer_code_data.code_table = peer_code_table;
	peer_code_data.num_component_edges = num_immersion_component_edges;
	peer_code_data.first_edge_on_component = first_immersion_edge_on_component;
	peer_code_data.term_crossing = term_crossing;
	peer_code_data.orig_crossing = orig_crossing;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code: peer code data produced from gauss code: ";
	write_peer_code(debug,peer_code_data);
	debug << endl;
	print_code_data(debug,peer_code_data,"gauss_to_peer_code: ");	
}

	bool component_zero_odd_shift = false;
	
	if ((peer_code_data.immersion_character == generic_code_data::character::KNOTOID || peer_code_data.immersion_character == generic_code_data::character::MULTI_LINKOID) && num_virtual_crossings_on_gauss_arc[0] != 0)
	{
		/* If the Gauss code is that of a KNOTOID or a MULTI_LINKOID and we have added virtual crossings to Gauss arc zero, then we have 
		   produced the peer code of the virtual closure of the knotoid.  In this case we renumber the peer code so that 
		   edge zero is the leg of the knotoid and, if we are dealing with a KNOTOID, adjust the immersion character to be PURE_KNOTOID.  
		   
		   Since the leg is always on component zero, and we adjust the components other than zero to be consistent 
		   with the numbering of component zero.  The shift for the renumbering required is the number of virtual 
		   crossings on Gauss arc zero.
		   
		   If the Gauss code is that of a KNOTOID but no virtual crossings have been added to Gauss arc zero, then we have
		   the peer code of a knot-type knotoid and so do not need to adjust the peer code. 
		*/
		vector<int> shift_vector(num_components);
		shift_vector[0] = num_virtual_crossings_on_gauss_arc[0];
		if (num_virtual_crossings_on_gauss_arc[0]%2 == 1)
		{
			component_zero_odd_shift = true;
			for (int i=1; i< num_components; i++)
				shift_vector[i] = 1;
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "gauss_to_peer_code: renumber peer code to accommodate knotoid gauss code, shift vector: ";
	for (int i=0; i< num_components; i++)
		debug << shift_vector[i] << ' ';
	debug << endl;
}
		renumber_peer_code(peer_code_data,shift_vector);

		int terminating_edge_at_head = peer_code_data.num_component_edges[0] - num_virtual_crossings_on_gauss_arc[0];
		int head_naming_edge = (terminating_edge_at_head%2?peer_code_data.code_table[generic_code_data::table::EPEER][(terminating_edge_at_head-1)/2]:terminating_edge_at_head);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: terminating_edge_at_head of component zero = " << terminating_edge_at_head << ", head_naming_edge = " << head_naming_edge << endl;

		/* set the head crossing */
		if (peer_code_data.immersion_character == generic_code_data::character::KNOTOID)
		{
			peer_code_data.head = head_naming_edge/2;
			peer_code_data.immersion_character = generic_code_data::character::PURE_KNOTOID;
		}
		else if (peer_code_data.immersion_character == generic_code_data::character::MULTI_LINKOID)
		{
			peer_code_data.component_type[0].type = component_character::PURE_START_LEG;
//			peer_code_data.component_type[0].head = head_naming_edge/2;
			peer_code_data.component_type[0].head_semi_arc = terminating_edge_at_head;
		}
	
		/* multi-linkoids are represented as the virtual closure of the multi-linkoid, so the virtual crossings we have already is correct,
		   knotoids are represented with the shortcut passing everywhere under, so the crossing labels need adjusting.
		*/
		if (peer_code_data.immersion_character != generic_code_data::character::MULTI_LINKOID)
		{
			for (int i = terminating_edge_at_head; i < peer_code_data.num_component_edges[0]; i++)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: shortcut edge " << i;
	
				if (i%2 == 0)
				{
					peer_code_data.code_table[generic_code_data::table::LABEL][i/2] = generic_code_data::NEGATIVE;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " crossing " << i/2 << " becomes NEGATIVE" << endl;
				}
				else
				{
					peer_code_data.code_table[generic_code_data::table::LABEL][peer_code_data.code_table[generic_code_data::table::EPEER][(i-1)/2]/2] = generic_code_data::POSITIVE;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " crossing " << peer_code_data.code_table[generic_code_data::table::EPEER][(i-1)/2]/2 << " becomes POSITIVE" << endl;
				}
			}
		}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code: adjusted peer code data for knotoid gauss code:" << endl;
	print_code_data(debug,peer_code_data,"gauss_to_peer_code: ");	
}

		/* update the immersion_crossing_peers for the Gauss crossings so we determine the correct gauss_crossing_map */
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "gauss_to_peer_code: update immersion_crossing_peers to accommodate shift vector:" << endl;
	
		for (int i=0; i< num_gauss_crossings; i++)
		{
			for (int j=0; j < 2; j++)
			{
				int old_label = immersion_crossing_peers[i][j];

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "gauss_to_peer_code:   Gauss crossing " << i << " label " << j << " = " << old_label << endl;

				int component = 0;
				while (old_label >= peer_code_data.first_edge_on_component[component]+peer_code_data.num_component_edges[component])
					component++;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "gauss_to_peer_code:   component " << component << endl;
					
				int new_label = (old_label - peer_code_data.first_edge_on_component[component] - shift_vector[component] + peer_code_data.num_component_edges[component])%
		                peer_code_data.num_component_edges[component] + peer_code_data.first_edge_on_component[component];

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "gauss_to_peer_code:   new label = " << new_label << endl;

				immersion_crossing_peers[i][j] = new_label;
			}
		}	
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "gauss_to_peer_code: updated immersion_crossing_peers: " << endl;
	print(immersion_crossing_peers, debug, 4, "gauss_to_peer_code: ");
}	
	}
	
	/* If we are dealing with a multi-linkoid, we may have other open components to handle.  The immersion labels have been assigned to components starting
	   with any virtual crossings added to the first Gauss arc.  Each component's immersion labels begin with an even label but to ensure that each crossing
	   has a odd and an even terminating edge, components other than the zero-th component may have been cycled by moving the start of the numbering backwards 
	   by one edge.  The boolean vector component_cycled records which components have been cycled.
	   
	   In the peer code we are creating, open components should have the leg on either the first (even) immersion edge or the last (odd) immersion edge, so we
	   have to shift subsequent components to ensure this is the case.  Since conponent zero has already been considered, and possibly shifted, we may already 
	   have adjusted the other components accordingly.  That means we only want to shift subsequent components by an even number of edges in order to retain the
	   odd and even terminating edges at a crossing.  The boolean component_zero_odd_shift tells us whether component zero has had an odd shift, in which case 
	   the subsequent components have all been shifted forwards by one edge.  Thus, some components have been shifted back when the labels were fist allocated 
	   and then forwards if component_zero_odd_shift is true.  This results in multiple cases to consider when shifting the subsequent open components to ensure
	   that the leg is correctly positioned.  We record the character of the open components in the component_type vector of the generic_code_data.	   
	*/
	if (peer_code_data.immersion_character == generic_code_data::character::MULTI_LINKOID)
	{
		for (int i=1; i< peer_code_data.num_open_components; i++)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code: multi-knotoid open component " << i << ':' << endl;
			
			vector<int> shift_vector(num_components);
			int num_virtual_crossings_on_first_gauss_arc = num_virtual_crossings_on_gauss_arc[gauss_code_data.first_edge_on_component[i]];
			int terminating_edge_at_head=-1;
			int last_shortcut_edge=-1;
			
			if (component_cycled[i]) // component start moved back by one
			{
				if (component_zero_odd_shift) // component start moved forwards by one
				{
					if (num_virtual_crossings_on_first_gauss_arc == 0)
					{
						/* knot-like component, leg is the first edge no further shift required */
						peer_code_data.component_type[i].type = component_character::KNOT_TYPE_START_LEG;
					}
					else if (num_virtual_crossings_on_first_gauss_arc % 2 == 1)
					{ 
						/* shift numbering start forwards by num_virtual_crossings_on_first_gauss_arc + 1, so that leg lies on last edge */
						shift_vector[i] = num_virtual_crossings_on_first_gauss_arc+1;
						peer_code_data.component_type[i].type = component_character::PURE_END_LEG;

						terminating_edge_at_head = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-num_virtual_crossings_on_first_gauss_arc-1;
						last_shortcut_edge = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-1;
					}
					else
					{						
						/* shift numbering start forwards by num_virtual_crossings_on_first_gauss_arc, so that leg lies on first edge */
						shift_vector[i] = num_virtual_crossings_on_first_gauss_arc;
						peer_code_data.component_type[i].type = component_character::PURE_START_LEG;

						terminating_edge_at_head = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-num_virtual_crossings_on_first_gauss_arc;
						last_shortcut_edge = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i];
					}
				}
				else
				{
					if (num_virtual_crossings_on_first_gauss_arc == 0)
					{
						/* knot-like component, shift forwards by 2 so the leg is the last edge*/
						shift_vector[i] = 2;
						peer_code_data.component_type[i].type = component_character::KNOT_TYPE_END_LEG;
					}
					else if (num_virtual_crossings_on_first_gauss_arc % 2 == 1)
					{ 
						/* shift numbering start forwards by num_virtual_crossings_on_first_gauss_arc + 1, so that leg lies on first edge */
						shift_vector[i] = num_virtual_crossings_on_first_gauss_arc+1;
						peer_code_data.component_type[i].type = component_character::PURE_START_LEG;

						terminating_edge_at_head = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-num_virtual_crossings_on_first_gauss_arc;
						last_shortcut_edge = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i];
					}
					else
					{						
						/* shift numbering start forwards by num_virtual_crossings_on_first_gauss_arc + 2, so that leg lies on last edge */
						shift_vector[i] = num_virtual_crossings_on_first_gauss_arc+2;
						peer_code_data.component_type[i].type = component_character::PURE_END_LEG;

						terminating_edge_at_head = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-num_virtual_crossings_on_first_gauss_arc-1;
						last_shortcut_edge = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-1;
					}
				}
			}
			else
			{
				if (component_zero_odd_shift) // component start moved forwards by one
				{
					if (num_virtual_crossings_on_first_gauss_arc == 0)
					{
						/* knot-like component, leg is the last edge no further shift required */
						peer_code_data.component_type[i].type = component_character::KNOT_TYPE_END_LEG;
					}
					else if (num_virtual_crossings_on_first_gauss_arc % 2 == 1)
					{ 
						/* shift numbering start forwards by num_virtual_crossings_on_first_gauss_arc - 1, so that leg lies on first edge */
						shift_vector[i] = num_virtual_crossings_on_first_gauss_arc-1;
						peer_code_data.component_type[i].type = component_character::PURE_START_LEG;

						terminating_edge_at_head = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-num_virtual_crossings_on_first_gauss_arc;
						last_shortcut_edge = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i];
					}
					else
					{						
						/* shift numbering start forwards by num_virtual_crossings_on_first_gauss_arc, so that leg lies on last edge */
						shift_vector[i] = num_virtual_crossings_on_first_gauss_arc;
						peer_code_data.component_type[i].type = component_character::PURE_END_LEG;

						terminating_edge_at_head = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-num_virtual_crossings_on_first_gauss_arc-1;
						last_shortcut_edge = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-1;
					}
				}
				else
				{
					if (num_virtual_crossings_on_first_gauss_arc == 0)
					{
						/* knot-like component, leg is the first edge no further shift required */
						peer_code_data.component_type[i].type = component_character::KNOT_TYPE_START_LEG;
					}
					else if (num_virtual_crossings_on_first_gauss_arc % 2 == 1)
					{ 
						/* shift numbering start forwards by num_virtual_crossings_on_first_gauss_arc + 1, so that leg lies on last edge */
						shift_vector[i] = num_virtual_crossings_on_first_gauss_arc+1;
						peer_code_data.component_type[i].type = component_character::PURE_END_LEG;

						terminating_edge_at_head = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-num_virtual_crossings_on_first_gauss_arc-1;
						last_shortcut_edge = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-1;
					}
					else
					{						
						/* shift numbering start forwards by num_virtual_crossings_on_first_gauss_arc, so that leg lies on first edge */
						shift_vector[i] = num_virtual_crossings_on_first_gauss_arc;
						peer_code_data.component_type[i].type = component_character::PURE_START_LEG;

						terminating_edge_at_head = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i]-num_virtual_crossings_on_first_gauss_arc;
						last_shortcut_edge = peer_code_data.first_edge_on_component[i]+peer_code_data.num_component_edges[i];
					}
				}
			}
			
			/* set the component head for pure knotoid components */
			if (terminating_edge_at_head != -1)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   terminating_edge_at_head = " << terminating_edge_at_head << ", last_shortcut_edge = " << last_shortcut_edge << endl;
	
				int head_naming_edge = (terminating_edge_at_head%2?peer_code_data.code_table[generic_code_data::table::EPEER][(terminating_edge_at_head-1)/2]:terminating_edge_at_head);
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "gauss_to_peer_code:   head_naming_edge = " << head_naming_edge << endl;

//				peer_code_data.component_type[i].head = head_naming_edge/2;
				peer_code_data.component_type[i].head_semi_arc = terminating_edge_at_head;
			}

			renumber_peer_code(peer_code_data,shift_vector);

			/* change any shortcut virtual crossings to be under-crossings
			
			October 2023: multi-knotoids are represented as their virtual closure, so no change is required
				
			for (int i = terminating_edge_at_head; i < last_shortcut_edge; i++)
			{
//if (debug_control::DEBUG >= debug_control::DETAIL)
//	debug << "gauss_to_peer_code: shortcut edge " << i;
	
				if (i%2 == 0)
				{
					peer_code_data.code_table[generic_code_data::table::LABEL][i/2] = generic_code_data::NEGATIVE;
//if (debug_control::DEBUG >= debug_control::DETAIL)
//	debug << " crossing " << i/2;
				}
				else
				{
					peer_code_data.code_table[generic_code_data::table::LABEL][peer_code_data.code_table[generic_code_data::table::EPEER][(i-1)/2]/2] = generic_code_data::POSITIVE;
				}
			}
			*/
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "gauss_to_peer_code:   adjusted peer code data for knotoid gauss code:" << endl;
	print_code_data(debug,peer_code_data,"gauss_to_peer_code:   ");	
}
		
			/* update the immersion_crossing_peers for the Gauss crossings so we determine the correct gauss_crossing_map */
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "gauss_to_peer_code:   update immersion_crossing_peers to accommodate shift vector:" << endl;
		
			for (int i=0; i< num_gauss_crossings; i++)
			{
				for (int j=0; j < 2; j++)
				{
					int old_label = immersion_crossing_peers[i][j];
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "gauss_to_peer_code:     Gauss crossing " << i << " label " << j << " = " << old_label << endl;
	
					int component = 0;
					while (old_label >= peer_code_data.first_edge_on_component[component]+peer_code_data.num_component_edges[component])
						component++;
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "gauss_to_peer_code:     component " << component << endl;
						
					int new_label = (old_label - peer_code_data.first_edge_on_component[component] - shift_vector[component] + peer_code_data.num_component_edges[component])%
			                peer_code_data.num_component_edges[component] + peer_code_data.first_edge_on_component[component];
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "gauss_to_peer_code:     new label = " << new_label << endl;
	
					immersion_crossing_peers[i][j] = new_label;
				}
			}	
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "gauss_to_peer_code: updated immersion_crossing_peers:   " << endl;
	print(immersion_crossing_peers, debug, 4, "gauss_to_peer_code:   ");
}	
		}		
	}

	if (evaluate_gauss_crossing_perm)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "gauss_to_peer_code: evaluate gauss_crossing_map:" << endl;

		for (int i=0; i< num_gauss_crossings; i++)
		{	
			int even_peer = (immersion_crossing_peers[i][0]%2 == 0 ? immersion_crossing_peers[i][0] : immersion_crossing_peers[i][1]);			
			(*gauss_crossing_map)[i] = even_peer/2;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "gauss_to_peer_code:   Gauss crossing " << i << ", even_peer = " << even_peer << endl;

		}		
	}

    return true;	
}
