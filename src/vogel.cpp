/***********************************************************************************************************
			   Vogel's algorithm

		  A. Bartholomew, started 20th October, 2002
		  modified for virtuals 27th October, 2003
		  updated for braid 10 re-write 1st May 2005
		  updated for labelled peer codes January 2011
		  updated to handle flat crossings January 2015
		  updated to calculate turning numbers August 2024

This module contains an implementation of the Yamada-Vogel algorithm to determine a braid representation of a 
link given by a labelled peer code.  The labelled peer code is read into a generic_code_data object by the 
generic_code brokering function, which is presented to the vogel function contained herein.  The code table 
in the generic code data is used to determine Seifert circles, the Seifert graph and from there, whether a Vogel 
move is required.  If we carry out a Vogel move, then we re-evaluate enough of the code table to write a new 
labelled peer code which can be re-read to a new generic_code_data object and we go around the loop again.  
If we determine from the Seifert graph that no Vogel move is required then we use the Seifert circles and 
crossing types to evaluate the braid representation of the knot.  Full details of the algorithm may be found 
at www.layer8.co.uk/maths

In August 2024 the function was updated to calculate turning numbers.  Since the algorithm involves only Reidemeister 
2 moves it does not change the turning number of the initial diagram.  The algorithm terminates with one or two sets of 
concentric Seifert circles and we calculate the turning number by keeping track of the point at infinity as given by the 
initial peer code.  For this purpose we use the concept of the infinite turning cycle, as used by the draw programme.  We
allow the user to specify the infinite turning cycle with a qualifier {cycle=n}, again following the draw programme.

The first edge of the infinite turning cycle is the orientation_reference_edge, the polarity of this edge, together with
an understanding of whether the infinite turning cycle is a left or right turning cycle tells us whether the point at infinity
is to the right or left of this edge as we follow its natural orientation determined by the edge numbering (note we traverse
odd edges backwards in turning cycles)

Specifically, the point at infinity lies the following side of the reference edge

orientation_                    Infinite turning cycle                       
reference_edge             left cycle              right cycle        
polarity                   (clockwise)           (anti-clockwise)

           odd                right                  left
           even                left                  right

Once we have a chain of concentric Seifert circles we can look for the circle containing the orientation reference edge, call this 
the reference circle.  We may move Seifert circles around the 2-sphere as necessary to reposition the referenced point at infinity 
to the non-compact component of the plane, so that the reference circle is the outermost circle of a set of concentric Seifert circles.  
The clockwise or anti-clockwise nature of the reference circle is determined by that of the orientation reference edge, which inherits
its sense from the from the infinite turning cycle.  We therefore know how many braid strands are closed clockwise or anti-clockwise, 
since all those circles in the compact region of the reference circle's complement have the same sense as the reference circle and the 
other Seifert circles have the opposite sense.

**************************************************************************/
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <vector>
#include <iomanip>

using namespace std;

/********************* External variables ***********************/
extern bool		OUTPUT_AS_INPUT;
extern ofstream		output;
extern ofstream		debug;



#include <util.h>
#include <matrix.h>
#include <hash-defs.h>
#include <braid-control.h>
#include <braid-util.h>
#include <generic-code.h>

#define OVER      0
#define UNDER     1
#define OVERMARK  2
#define UNDERMARK 3

/* used for calculating turning numbers */
#define ANTICLOCKWISE  1
#define CLOCKWISE     -1

/********************* Function prototypes ***********************/
void vogel_error (string errstring);
string braid_reduce (string word);
string virtual_vogel (string inbuf);
int on_circle(matrix<int>& seifert_circle, int circle, int crossing);

string vogel (generic_code_data code_data, int* turning_number_ptr)
{
	char* 		mark;
	int			num_left_cycles=0;
	int			num_cycles=0;
	int			crossing;
	int			vertex;
	int			edge;
	int			graph_edge_type;
	int			circle;
	int			c1;
	int			c2;
	int 		e1;
	int 		e2;
	bool		found;
	bool		closed_braid;
	
	/* variables used for calculating turning numbers */
	bool		initial_code_data = true;
	int			orientation_reference_edge; // the first edge in the infinite turning cycle
	int 		orientation_reference_polarity; // CLOCKWISE or ANTICLOCKWISE according to the infinite turning cycle
	bool		infinity_right_of_reference_edge;

	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;
	
	ostringstream dev_ss;
	write_peer_code(dev_ss, code_data);
	string dev_str = dev_ss.str();

if (braid_control::VOGEL_DEBUG)
{
	debug << "vogel: presented with code_data = ";
	write_peer_code(debug, code_data);
	debug << "\nvogel: code_data:" << endl;
	print_code_data(debug,code_data,"vogel: ");
}

	if (braid_control::VOGEL_TURNING_NUMBER && turning_number_ptr == 0)
	{
if (braid_control::VOGEL_DEBUG)
	debug << "vogel: Error! turning_number_ptr is zero when braid_control::VOGEL_TURNING_NUMBER is true" << endl;
	
		return "";
	}
	
	/* Count the number of Seifert circles, this will not change even if we
	   apply Vogel moves. We will be storing Seifert circles as lists of edges
	   and as lists of crossings, in both cases in a matrix<int>, whose first
	   column will record the length of that Seifert circle (the same for both 
	   edges and crossings).  We shall use num_edges as an upper bound on the 
	   length of a Seifert circle.
	*/
	vector<int> edge_flag(num_edges);

	int num_circles = 0;
	int start;
	bool complete = false;
	
	do
	{
		/* look for an edge in edge_flag we have not yet encountered */
		found = false;
		for (int i=0; i<num_edges; i++)
		{
			if (edge_flag[i] == 0)
			{
				start = i;
				edge = start;
				found = true;
				break;
			}
		}
	
		if (found)
		{
			num_circles++;

			/* trace around the Seifert circle from start, counting edges */
			do
			{	
				edge_flag[edge] = 1;
				
				/* identify the crossing at which the next edge terminates */
				crossing = code_data.term_crossing[edge];
				
				/* determine the next edge in the Seifert circle */
				edge = (edge % 2? code_data.code_table[generic_code_data::table::ODD_ORIGINATING][crossing]: code_data.code_table[generic_code_data::table::EVEN_ORIGINATING][crossing]);

			} while (edge != start);
			
		}
		else
			complete = true;
	} while (!complete);

if (braid_control::VOGEL_DEBUG)
    debug << "vogel: number of Seifert circles = " << num_circles << endl;

	/* reset the edge_flags ready for the first time round the main loop 
	for (int i=0; i< num_edges; i++)
		edge_flag[i] = 0;
	*/

	/* We declare seifert_graph here, since it will never
	   need to change in size.
	*/
	matrix<int> seifert_graph(num_circles,num_circles);	
	
	/* This is the start of the main loop clause that evaluates the Seifert
	   circles and graph and determines whether a Vogel move is required.
	   The while clause that closes the loop is dependent on a bool
	   closed_braid that is set during the loop itself when we test the
	   Seifert graph.
	*/
	matrix<int> seifert_circle(num_circles,num_edges+1);
	matrix<int> seifert_e_circle(num_circles,num_edges+1);
	do
	{	
		/* we need to define code_data within the scope of the "do" loop
		   since we will be overwriting the code_data object after performing
		   a Vogel move 
		*/
		matrix<int>& code_table = code_data.code_table;

		/* First, identify the turning cycles from the code table. so we
		   have a record of the regions of the knot's compliment.  There
		   cannot be more than num_crossings + 2 cycles, and they cannot 
		   be more than num_edges in length, we add 1 to the length of a 
   		   cycle as we store the actual length in column 0.  For
		   a description of turning cycles see the description of labelled
		   peer codes at www.layer8.co.uk/maths
		*/
		matrix<int> cycle(num_crossings+2, num_edges+1);
		calculate_turning_cycles(code_data, cycle, num_left_cycles, num_cycles);
		
		/* in the first iteration of this loop, set the orientation_reference_edge and record whether it lies
		   in a left or right turning cycle */
		if (braid_control::VOGEL_TURNING_NUMBER && initial_code_data)
		{
			if (!braid_control::BIRACK_POLYNOMIAL)
			{
				/* write the initial code data turning cycles to the output file */
			    output << "\n\nturning_cycles:" << endl;
				for (int i=0; i<num_cycles; i++)
				{
					output << "  cycle " << i << " length = " << cycle[i][0] << ": ";
					for (int j=1; j<=cycle[i][0]; j++)
						output << cycle[i][j] << " ";
					output << endl;
				}
			}
			
			if (braid_control::INFINITE_CYCLE == braid_control::cycle::UNSPECIFIED)
			{
				/* determine the infinite region as the longest turning cycle */
				braid_control::INFINITE_CYCLE = 0;

				for (int i=1; i< num_cycles; i++)
				{
					if (cycle[i][0] > cycle[braid_control::INFINITE_CYCLE][0])
						braid_control::INFINITE_CYCLE = i;
				}
			}
			
if (braid_control::VOGEL_DEBUG)
    debug << "vogel: infinite turning cycle = " << braid_control::INFINITE_CYCLE << endl;			

			/* identify the first edge in the infinite turning cycle and record the cycle polarity */
			orientation_reference_edge = abs(cycle[braid_control::INFINITE_CYCLE][1]);
			
if (braid_control::VOGEL_DEBUG)
    debug << "vogel: orientation_reference_edge = " << orientation_reference_edge << endl;			
			
			if (braid_control::INFINITE_CYCLE < num_left_cycles )
			{
				orientation_reference_polarity = CLOCKWISE;
if (braid_control::VOGEL_DEBUG)
    debug << "vogel: infinite turning cycle oriented clockwise" << endl;			
			}
			else
			{
				orientation_reference_polarity = ANTICLOCKWISE;
if (braid_control::VOGEL_DEBUG)
    debug << "vogel: infinite turning cycle oriented anti-clockwise" << endl;			
			}

			if ( (orientation_reference_edge%2==0 && orientation_reference_polarity ==  ANTICLOCKWISE) ||
			     (orientation_reference_edge%2==1 && orientation_reference_polarity ==  CLOCKWISE) )
			{
				infinity_right_of_reference_edge = true;
			}
			else
			{
				infinity_right_of_reference_edge = false;
			}
			
if (braid_control::VOGEL_DEBUG)
    debug << "vogel: infinity_right_of_reference_edge = " << (infinity_right_of_reference_edge?"true":"false") << endl;			
			
			initial_code_data = false;
		}
		
		
		/* Now calculate the seifert circles from the code_table, we use Seifert_Circle for crossings, and seifert_e_circle for the edges.
		
		   We use a local object Seifert_Circle and Seifert_E_Circle due to an error encountered with peer code [-1 -3 -5 -7]/+ + + + when we use the object
		   seifert_circle declared outside the scope of this "do" clause to record Seifert circles.  Although Vogel moves do not change the 
		   number of Seifert circles and the code below re-writes the length of each Seifert circle as it is calculated, with a single object 
		   declared at a wider scope, the crossing information for a Seifert circle was being overridden when the next Seifert circle length 
		   was written.  The above peer code had three Vogel moves applied to give peer code [-7 9 13 17 -11 1 -15 3 -19 5]/+ + + + + - + - + -, 
		   whereupon writing the length of the third Seifert circle overwrote part of the second Seifert circle.  This did not happen if the 
		   vogel function were presented with the latter peer code from the outset.  Change history has been left here in case this is a 
		   compiler error.
		*/
		matrix<int> Seifert_E_Circle(num_circles,num_edges+1);
		matrix<int> Seifert_Circle(num_circles,num_edges+1);
		
		/* We shall use num_circles as an index into seifert_circle
		   but need to count rows from zero, so this assignment keeps the
		   code below the same as that above
		*/
		num_circles = -1;

		/* we define a new, local, edge_flag vector here so it picks up the
		   new value of num_edges in each iteration of the loop
		*/
		vector<int> edge_flag(num_edges);


if (braid_control::VOGEL_DEBUG)
    debug << "vogel: check, num_edges = " << num_edges << ", edge_flag size = " << edge_flag.size() << endl;

		complete = false;
		do
		{
			/* look for an edge in edge_flag we have not yet encountered */
			found = false;
			for (int i=0; i<num_edges; i++)
			{
				if (edge_flag[i] == 0)
				{
					start = i;
					edge = start;

//if (braid_control::VOGEL_DEBUG)
//    debug << "vogelXX:   start of seifert circle at edge " << edge << "(start = " << start << ")" << endl;

					found = true;
					break;
				}
			}
		
			if (found)
			{
				num_circles++;
				int offset = 0;

				/* trace around the Seifert circle from start, counting edges */
				do
				{	
					edge_flag[edge] = 1;
					offset++;
					
					Seifert_E_Circle[num_circles][offset] = edge;
					
					/* identify the crossing at which the next edge terminates */
					crossing = code_data.term_crossing[edge];

//if (braid_control::VOGEL_DEBUG)
//    debug << "vogelXX:   edge " <<  edge << " terminates at crossing " << crossing << ", offset = " << offset << endl;

					
//					seifert_circle[num_circles][offset] = crossing;
					Seifert_Circle[num_circles][offset] = crossing;
					
					/* determine the next edge in the Seifer circle */
					edge = (edge % 2? code_table[generic_code_data::table::ODD_ORIGINATING][crossing]: code_table[generic_code_data::table::EVEN_ORIGINATING][crossing]);

//if (braid_control::VOGEL_DEBUG)
//    debug << "vogelXX:   next edge " << edge << endl;
				
				} while (edge != start);	
	
/*			
if (braid_control::VOGEL_DEBUG)
{
    debug << "vogelXX:   Seifert circle crossings from seifert_circle before:" << endl;
    for (int i=0;i<num_circles;i++)
    {
//		debug << "vogelXX:   " << i << ". vertex = " << i+1 << ", length = "  << seifert_circle[i][0] << ": ";
//		for (int j=1;j<=seifert_circle[i][0];j++)
//		    debug << seifert_circle[i][j] << " ";
		
		debug << "vogelXX:   " << i << ". vertex = " << i+1 << ", length = "  << Seifert_Circle[i][0] << ": ";
		for (int j=1;j<=Seifert_Circle[i][0];j++)
		    debug << Seifert_Circle[i][j] << " ";
		debug << endl;
    }
    debug << "vogelXX:   Seifert circle edges from seifert_circle before:" << endl;
    for (int i=0;i<num_circles;i++)
    {
		debug << "vogelXX:   " << i << ". vertex = " << i+1 << ", length = "  << Seifert_E_Circle[i][0] << ": ";
		for (int j=1;j<=Seifert_E_Circle[i][0];j++)
		    debug << Seifert_E_Circle[i][j] << " ";
		debug << endl;
    }
}
*/

				/* write the length of the Seifert circle into the first
	   			   column of seifert_e_circle and seifert_circle */
	 			Seifert_E_Circle[num_circles][0] = offset;
//	 			seifert_circle[num_circles][0] = offset;
	 			Seifert_Circle[num_circles][0] = offset;

/*
if (braid_control::VOGEL_DEBUG)
{
    debug << "vogelXX:   Seifert circle crossings from seifert_circle after:" << endl;
    for (int i=0;i<=num_circles;i++)
    {
//		debug << "vogelXX:   " << i << ". vertex = " << i+1 << ", length = "  << seifert_circle[i][0] << ": ";
//		for (int j=1;j<=seifert_circle[i][0];j++)
//		    debug << seifert_circle[i][j] << " ";

		debug << "vogelXX:   " << i << ". vertex = " << i+1 << ", length = "  << Seifert_Circle[i][0] << ": ";
		for (int j=1;j<=Seifert_Circle[i][0];j++)
		    debug << Seifert_Circle[i][j] << " ";
		debug << endl;
    }
    debug << "vogelXX:   Seifert circle edges from seifert_circle after:" << endl;
    for (int i=0;i<num_circles;i++)
    {
		debug << "vogelXX:   " << i << ". vertex = " << i+1 << ", length = "  << Seifert_E_Circle[i][0] << ": ";
		for (int j=1;j<=Seifert_E_Circle[i][0];j++)
		    debug << Seifert_E_Circle[i][j] << " ";
		debug << endl;
    }
}
*/

			}
			else
			{
				/* We copy-assign Seifert_Circle to seifert_circle here since subsequent use of seifert_circle is read-only */
				seifert_circle = Seifert_Circle;				
				seifert_e_circle = Seifert_E_Circle;				
				complete = true;
			}
				
				
		} while (!complete);
			
		/* adjust num_circles to its correct value */
		num_circles++;
		
if (braid_control::VOGEL_DEBUG)
{
    debug << "vogel: Seifert circle edges from seifert_e_circle:" << endl;
    for (int i=0;i<num_circles;i++)
    {
		debug << "vogel: " << i << ". vertex = " << i+1 << ", length = "  << seifert_e_circle[i][0] << ": ";
		for (int j=1;j<=seifert_e_circle[i][0];j++)
		    debug << seifert_e_circle[i][j] << " ";
		debug << endl;
    }
    debug << "vogel: Seifert circle crossings from seifert_circle:" << endl;
    for (int i=0;i<num_circles;i++)
    {
		debug << "vogel: " << i << ". vertex = " << i+1 << ", length = "  << seifert_circle[i][0] << ": ";
		for (int j=1;j<=seifert_circle[i][0];j++)
		    debug << seifert_circle[i][j] << " ";
		debug << endl;
    }
}

		
	    /* From the set of Seifert circle crossings we may determine the edges
	       of the reduced Seifert graph of the link diagram.  We store this
	       information in seifert_graph, where row i records the vertices
	       adjacent to vertex i.  The valency of vertex i is stored in
	       column 0 of the row.

		   If the peer-code Type of the crossings between two Seifert circles
		   is TYPE1 then the corresponding edge in seifert_graph is negative,
		   otherwise (TYPE2) it is positive.
    	*/

    	for (int i=0; i< num_circles; i++)
    	{
			edge = 0;
			for (int j=1;j<=seifert_circle[i][0];j++)
			{
		    	crossing = seifert_circle[i][j];

		    	/* look for another Seifert circle containing crossing */
	    		found = false;
	    		for (int k=0; k<num_circles;k++)
	    		{
					if (k==i)
		    			continue;

					for (int l=1; l<=seifert_circle[k][0];l++)
					{
				    	if (seifert_circle[k][l] == crossing)
		    			{
							circle = k;
							found = true;
		    			}
					}
					if (found) 
						break;
	    		}

		    	/* If we haven't yet recorded an edge to the vertex represented
		       	   by circle, do so now
		    	*/
	    		found = false;
	    		for (int k=1; k<= edge; k++)
	    		{
					if (abs(seifert_graph[i][k]) == circle+1)
					{
			    		found = true;
		    			break;
					}
	    		}

		    	if (!found)
	    		{
					edge++;
					seifert_graph[i][edge] = circle+1;
					
					/* look for the Type of the crossing in the code table as an
					   example of the type of crossing between these two Seifert
					   circles, and set the edge in the graph to be negative if they
					   are joined by Type 1 crossings.
					*/
					if (code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::TYPE1)
						seifert_graph[i][edge] *= -1;
	    		}
			}

			/* record the number of edges we have found */
			seifert_graph[i][0] = edge;
    	}

if (braid_control::VOGEL_DEBUG)
{
   	debug << "vogel: reduced Seifert graph:" << endl;
   	for (int i=0;i<num_circles;i++)
   	{
		debug << "vogel: " << i << ". vertex = " << i+1 << ", length = "  << seifert_graph[i][0] << ": ";
		for (int j=1;j<=seifert_graph[i][0];j++)
    		debug << seifert_graph[i][j] << " ";
		debug << endl;
   	}
}
		
		/* Look through the Seifert graph to see if a Vogel move is required.
		   We need a Vogel move if there is a vertex in the boundary of two
		   edges of the same type.
		*/
		closed_braid = true;

		for (int i=0; i< num_circles; i++)
		{
			if ( seifert_graph[i][0] > 2 || 
				(seifert_graph[i][0] == 2 && seifert_graph[i][1]*seifert_graph[i][2] > 0)
			   )
			{
				closed_braid = false;
				circle = i;
				break;
			}
		}

		if (braid_control::VOGEL_HEIGHT_ONLY)
		{		
			if (closed_braid)
			    return "height is zero";
			else
				return "non-zero height";
		}
				
		if (!closed_braid)
		{
			/* We can perform a Vogel move between two of the Seifert circles
			   reachable from circle, search for two that both have an edge in
			   the same turning cycle, since we can perform the move between
				those two edges.  We are guaranteed to find such a pair.
			*/
			found = false;
			int k;
			for (int i=1; i< seifert_graph[circle][0]; i++)
			{
				c1 = abs(seifert_graph[circle][i]) - 1;

				/* save the type of the edges for later */
				if (seifert_graph[circle][i] < 0)
					graph_edge_type = generic_code_data::TYPE1;
				else
					graph_edge_type = generic_code_data::TYPE2;
					
				for (int j=i+1; j<= seifert_graph[circle][0]; j++)
				{
					/* Only test c1 with this circle if the corresponding edges
					   in the seifert_graph are of the same type */
					if (seifert_graph[circle][i]*seifert_graph[circle][j] < 0)
						continue;
						
					c2 = abs(seifert_graph[circle][j]) - 1;

					
if (braid_control::VOGEL_DEBUG)
{
	debug << "vogel: testing Seifert circles " << c1+1 << " and " << c2+1 << " for a Vogel move" << endl;
}	
					/* see if circles c1 and c2 have an edge in the same turning cycle */
					for (k=0; k < num_cycles; k++)
					{
						/* Does c1 have an edge in cycle k ? */
						found = false;
						for (int l=1; l<= cycle[k][0]; l++)
						{
							edge = abs(cycle[k][l]);
							
							/* look for edge in c1 */
							for (int m=1; m<= seifert_e_circle[c1][0]; m++)
							{
								if (seifert_e_circle[c1][m] == edge)
								{
									e1 = edge;
									found = true;	
									break;
								}
							}
							if (found)
								break; // don't need to look for another edge in cycle k
						}	
						
						if (found)
						{
							/* Does c2 have an edge in cycle k ? */
							found = false;
							for (int l=1; l<= cycle[k][0]; l++)
							{
								edge = abs(cycle[k][l]);
							
								/* look for edge in c2 */
								for (int m=1; m<= seifert_e_circle[c2][0]; m++)
								{
									if (seifert_e_circle[c2][m] == edge)
									{
										e2 = edge;
										found = true;
										break;
									}
								}
								if (found)
									break; // don't need to look for another edge in cycle k
							}	
						}
						
						if (found)
							break; // don't need to check other cycles
					}
					if (found)
					break;  // don't need to check any more c2 candidates					
				}
				if (found)
					break;  // don't need to check any more c1 candidates					
			}
if (braid_control::VOGEL_DEBUG)
{
	debug << "vogel: may apply a Vogel move between edges " << e1 << " and " << e2;
	debug << " across cycle ";
	for (int i=1; i<= cycle[k][0]; i++)
		debug << cycle[k][i] << ' ';
	debug << endl;
}

			
			/* Apply the Vogel move by determining the code table for the peer code of the modified diagram	*/
			vector<int> new_edge(num_edges);
		
			/* set e1 to be the low edge */
			if (e2 < e1)
				swap(e1,e2);

			/* write the new edge numbers */
			for (int i=0; i< e1; i++)
				new_edge[i] = i;
			for (int i=e1; i<e2; i++)
				new_edge[i] = i+2;
			for (int i=e2; i<num_edges; i++)
				new_edge[i] = i+4; 
		
if (braid_control::VOGEL_DEBUG)
{
	debug << "vogel: new edge numbers: ";
	for (int i=0; i<num_edges; i++)
		debug << new_edge[i] << ' ';
	debug << endl;
}

			if (braid_control::VOGEL_TURNING_NUMBER)
			{
				/* update the orientation_reference_edge with it's new edge number */
				if (orientation_reference_edge > e1)
				{
					orientation_reference_edge = new_edge[orientation_reference_edge];
				
if (braid_control::VOGEL_DEBUG)
	debug << "vogel: orientation_reference_edge updated to " << orientation_reference_edge << endl;
				}
				
			}
			
			/* determine the crossing numbers in the new code for the two additional crossings introduced by the Vogel move.
			   also determine the component on which each edge resides.
			*/
			int component_e1;
			int component_e2;
			
			if (e1%2)
			{
				/* both edges are odd */
				c1 = (e1+1)/2;
				c2 = (e2+3)/2; // allowing for the fact that new_edge[e2] == e2+2

				component_e1 = code_table[generic_code_data::table::COMPONENT][(e1-1)/2];
				component_e2 = code_table[generic_code_data::table::COMPONENT][(e2-1)/2];
			}
			else
			{
				c1 = e1/2;
				c2 = (e2+2)/2; // allowing for the fact that new_edge[e2] == e2+2
				
				component_e1 = code_table[generic_code_data::table::COMPONENT][e1/2];
				component_e2 = code_table[generic_code_data::table::COMPONENT][e2/2];
			}

if (braid_control::VOGEL_DEBUG)
{
	debug << "vogel: new crossings will have numbers " << c1 << " and " << c2 << endl;
	debug << "vogel: edge " << e1 << " lies on component " << component_e1 << ", edge " << e2 << " lies on component " << component_e2 << endl;
}	
		
			vector<int> new_crossing(num_crossings);
		
			/* write the new crossing numbers */
			for (int i=0; i< c1; i++)
				new_crossing[i] = i;
			for (int i=c1; i<c2-1; i++)
				new_crossing[i] = i+1;
			for (int i=c2-1; i<num_crossings; i++)
				new_crossing[i] = i+2;
		
if (braid_control::VOGEL_DEBUG)
{
	debug << "vogel: new crossing numbers: ";
	for (int i=0; i<num_crossings; i++)
		debug << new_crossing[i] << ' ';
	debug << endl;
}

			generic_code_data new_code_data;
			new_code_data.head = -1;
			new_code_data.type = generic_code_data::peer_code;
			new_code_data.num_crossings = num_crossings+2;
			new_code_data.num_components = code_data.num_components;
					
			new_code_data.num_component_edges = code_data.num_component_edges;
			new_code_data.num_component_edges[component_e1] += 2;
			new_code_data.num_component_edges[component_e2] += 2;
			
			matrix<int> new_code_table(generic_code_data::table::CODE_TABLE_SIZE,num_crossings+2);
			
			for (int i=0; i< num_crossings; i++)
			{
				/* evaluate the (old) odd edge incoming at this crossing */
				edge = code_table[generic_code_data::table::OPEER][i];
				
				/* now determine the the new odd edge incoming at the 
				   same crossing using the new numbering
				*/
				new_code_table[generic_code_data::table::OPEER][new_crossing[i]] = new_edge[edge];
			
				/* The type of the crossing doesn't change, nor does the label or component*/
				new_code_table[generic_code_data::table::TYPE][new_crossing[i]] = code_table[generic_code_data::table::TYPE][i];
				new_code_table[generic_code_data::table::LABEL][new_crossing[i]] = code_table[generic_code_data::table::LABEL][i];
				new_code_table[generic_code_data::table::COMPONENT][new_crossing[i]] = code_table[generic_code_data::table::COMPONENT][i];
			}
			
			/* Now fill in the peer edges and components for the two new crossings c1 and c2 */
			if (e1%2)
			{
				new_code_table[generic_code_data::table::OPEER][c1] = e2+2;
				new_code_table[generic_code_data::table::OPEER][c2] = e1;
				
				/* the naming edge for c1 is e1+1 which belongs to the same
				   component as the edge preceeding e1
				*/
				int preceeding_edge = code_table[generic_code_data::table::EVEN_TERMINATING][code_data.orig_crossing[e1]];
				new_code_table[generic_code_data::table::COMPONENT][c1] = code_table[generic_code_data::table::COMPONENT][preceeding_edge/2];
				
				/* the naming edge for c2 is e2+3 which belongs to the same 
				   component as the edge preceeding e2
				*/
				preceeding_edge = code_table[generic_code_data::table::EVEN_TERMINATING][code_data.orig_crossing[e2]];
				new_code_table[generic_code_data::table::COMPONENT][c2] = code_table[generic_code_data::table::COMPONENT][preceeding_edge/2];
			}
			else
			{
				new_code_table[generic_code_data::table::OPEER][c1] = e2+3;
				new_code_table[generic_code_data::table::OPEER][c2] = e1+1;

				/* the naming edge for c1 is e1 which belongs to the same
				   component as the old edge e1
				*/
				new_code_table[generic_code_data::table::COMPONENT][c1] = code_table[generic_code_data::table::COMPONENT][e1/2];
				
				/* the naming edge for c2 is e2+2 which belongs to the same 
				   component as the old edge e2
				*/
				new_code_table[generic_code_data::table::COMPONENT][c2] = code_table[generic_code_data::table::COMPONENT][e2/2];
			}
		
			/* both crossings c1 and c2 have the same type, it is the opposite
			   of the type noted above when we identified the two Seifert circles
			   involved in the move.
			*/
			if (graph_edge_type == generic_code_data::TYPE1)
			{
				new_code_table[generic_code_data::table::TYPE][c1] = generic_code_data::TYPE2;
				new_code_table[generic_code_data::table::TYPE][c2] = generic_code_data::TYPE2;
			}
			else
			{
				new_code_table[generic_code_data::table::TYPE][c1] = generic_code_data::TYPE1;
				new_code_table[generic_code_data::table::TYPE][c2] = generic_code_data::TYPE1;
			}
		
			if (braid_control::FLAT_CROSSINGS)
			{
				new_code_table[generic_code_data::table::LABEL][c1] = generic_code_data::FLAT;
				new_code_table[generic_code_data::table::LABEL][c2] = generic_code_data::FLAT;
			}
			else
			{
				/* we are free to choose which of c1 and c2 has label + and which - */
				new_code_table[generic_code_data::table::LABEL][c1] = generic_code_data::POSITIVE;
				new_code_table[generic_code_data::table::LABEL][c2] = generic_code_data::NEGATIVE;
			}
			
			/* Added 6-8-24  to accommodate the updated write_peer_code function supporting multilinkoids.  This needs ODD_TERMINATING and EVEN_TERMINATING data */
			for (int i=0; i< new_code_data.num_crossings; i++)
			{
				new_code_table[generic_code_data::table::ODD_TERMINATING][i] = new_code_table[generic_code_data::table::OPEER][i];
				new_code_table[generic_code_data::table::EVEN_TERMINATING][i] = 2*i;
			}


			new_code_data.code_table = new_code_table;
			
			/* write the peer code to an ostringstream and read it back into 
			   code_data ready for the next loop
			*/
if (braid_control::VOGEL_DEBUG)
{
	debug << "vogel: new_code_data:" << endl;
	print_code_data(debug,new_code_data,"vogel: ");
	debug << "\nvogel: new_code_data = ";
	write_peer_code(debug, new_code_data);
}
			ostringstream oss;
			write_peer_code(oss, new_code_data);
			read_peer_code(code_data, oss.str());

dev_str = oss.str();
		
			/* increase num_crossings and reset num_edges*/
			num_crossings += 2;
			num_edges = 2*num_crossings;
			
if (braid_control::VOGEL_DEBUG)
{
	debug << "vogel: new number of crossings = " << num_crossings << endl;	
	debug << "vogel: new code table after Vogel move:" << endl;
	debug << "vogel: peer of even edges: ";
	for (int i=0; i< num_crossings; i++)
		debug << new_code_table[generic_code_data::table::OPEER][i] << ' ';
	debug << endl;
	debug << "vogel: type: ";	
	for (int i=0; i< num_crossings; i++)
		debug << new_code_table[generic_code_data::table::TYPE][i] << ' ';
	debug << endl;
	debug << "vogel: label: ";
	for (int i=0; i< num_crossings; i++)
		debug << new_code_table[generic_code_data::table::LABEL][i] << ' ';
	debug << endl;
	debug << "vogel: component: ";
	for (int i=0; i< num_crossings; i++)
		debug << new_code_table[generic_code_data::table::COMPONENT][i] << ' ';
	debug << endl;

	debug << "vogel: peer code written from new_code_table: " << oss.str() << endl;
	debug << "vogel: new code data after reading new peer_code: " << endl;
	print_code_data(debug,code_data,"vogel: ");
}
			if (braid_control::EXTRA_OUTPUT)
			{
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Vogel move " << e1 << " and " << e2 << ": ";
					write_peer_code(output,code_data);
					output << endl;

				}
			}
		}	
	} while (!closed_braid);


	/* Now evaluate the braid crossing types, positive, negative or virtual, 
	   that will be used when determining the braid word
	*/
	vector<int> crossing_type(num_crossings);
	matrix<int>& code_table = code_data.code_table;
	for (int i=0; i< num_crossings; i++)
	{
		if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::VIRTUAL)
			crossing_type[i] = generic_braid_data::crossing_type::VIRTUAL;
		else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::FLAT)
			crossing_type[i] = generic_braid_data::crossing_type::FLAT;
		else
		{
			if ((code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE && code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1)
			  ||(code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE && code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE2))
				crossing_type[i] = generic_braid_data::crossing_type::POSITIVE;
			else
				crossing_type[i] = generic_braid_data::crossing_type::NEGATIVE;		
		}
	}

if (braid_control::VOGEL_DEBUG)
{
    debug << "vogel: braid crossing types = ";
    for (int i=0; i<num_crossings; i++)
		debug << crossing_type[i] << " ";
	debug << endl;
}
	
	/* To determine the braid word we set up a chain vector that lists the Seifert circles
	   corresponding to the vertices of the Seifert graph in the order determined by the 
	   graph, numbered from 1.
	   	   
	   We also establish a strand vector identifies the position in the chain for each Seifert 
	   circle, thereby numbering each Seifert circle as a braid strand.  The strands are 
	   also numbered from 1, so that the strand number for the first Seifert circle 
	   in the chain is 1.
	*/

	vector<int> chain(num_circles);

	/* start chain from the first univalent vertex found in seifert_graph */
	for (int i=0; i< num_circles; i++)
	{
		if (seifert_graph[i][0] == 1)
		{
			chain[0] = i+1;
			chain[1] = abs(seifert_graph[i][1]);
			break;
		}
	}
	
	/* now set up the rest of the chain */
	for (int i=2; i< num_circles; i++)
	{
		/* identify the last vertex and its preceeding circle in the chain */
		vertex = chain[i-1]-1;
		circle = chain[i-2];
		
		if (abs(seifert_graph[vertex][1]) == circle)
			chain[i] = abs(seifert_graph[vertex][2]);
		else
			chain[i] = abs(seifert_graph[vertex][1]);
	}
	

if (braid_control::VOGEL_DEBUG)
{
	debug << "vogel: chain: ";
   	for (int j=0;j<num_circles;j++)
		debug << chain[j] << " ";
	debug << endl;
}

	if (braid_control::VOGEL_TURNING_NUMBER)
	{
		/* Identify where in the chain of Seifert circles the orientation reference edge lies */	
		int turning_number;	
		bool found = false;
		int location = -1;
		
		for (int i=0; i< num_circles; i++)
		{
			/* check whether the orientation_reference_edge lies in the ith Seifert circle in chain.*/
			for (int j=1; j<= seifert_e_circle[chain[i]-1][0]; j++)
			{
				if (seifert_e_circle[chain[i]-1][j] == orientation_reference_edge)
				{
					location = i;
					found = true;
				}
			}
			if (found)
				break;
		}
		
if (braid_control::VOGEL_DEBUG)
	debug << "vogel: orientation_reference_edge " << orientation_reference_edge << " located in Seifert circle at chain offset " << location << endl;
		
		if (location == 0 || location == num_circles-1)
		{
if (braid_control::VOGEL_DEBUG)
	debug << "vogel: Seifert circles form a single set of concentric circles " << endl;
	
			turning_number = num_circles*orientation_reference_polarity;
		}
		else
		{
if (braid_control::VOGEL_DEBUG)
	debug << "vogel: Seifert circles form two sets of concentric circles " << endl;

			/* The Seifert circle containing the orientation_reference_edge is deemed to be oriented according to the orientation_reference_polarity
			   but we do not yet know from which end of the chain we were looking (the "clockwise" or "anti-clockwise" Seifert circles end).
			   We check the Seifert circle at location-1 in the chain and determine the TYPE of the edge in the reduced Seifert graph
			   between the two Seifert circles.  This will tell us, together with the odd or even nature of the orientation_reference_edge, whether the 
			   Seifert circle chain[0] is to the left or the right of the orientation_reference_edge.  Comparing that information with the location of 
			   infinity determined by the orientation_reference_polarity allows us to calculate the turning number.
			   
			   Start by locating the edge from chain[location] to chain[location-1] in the seifert graph
			*/
			int preceeding_location_type;
			for (int i=1; i<= seifert_graph[chain[location]-1][0]; i++)
			{
				if (abs(seifert_graph[chain[location]-1][i]) == chain[location-1])
				{
					if(seifert_graph[chain[location]-1][i] < 0)
						preceeding_location_type = generic_code_data::TYPE1;
					else
						preceeding_location_type = generic_code_data::TYPE2;

if (braid_control::VOGEL_DEBUG)
{
	debug << "vogel: preceeding Seifert circle connected to that of orientation_reference_edge by crossings of type " <<
    (preceeding_location_type ==  generic_code_data::TYPE1?"TYPE1":"TYPE2") << endl;				
}	
					break;
					
				}
			}
			
			bool start_of_chain_to_left = false;
			
			if ( (orientation_reference_edge%2==0 && preceeding_location_type ==  generic_code_data::TYPE1) ||
			     (orientation_reference_edge%2==1 && preceeding_location_type ==  generic_code_data::TYPE2) )
			{
				start_of_chain_to_left = true;
			}

if (braid_control::VOGEL_DEBUG)
{
	debug << "vogel: start of chain lies to left of Seifert circle containing edge " << orientation_reference_edge << " is " <<
    (start_of_chain_to_left?"true":"false") << endl;		    		
    debug << "vogel: infinity_right_of_reference_edge = " << (infinity_right_of_reference_edge?"true":"false") << endl;			
}
			if ((start_of_chain_to_left && infinity_right_of_reference_edge) || (!start_of_chain_to_left && !infinity_right_of_reference_edge))
			{
if (braid_control::VOGEL_DEBUG)
	debug << "vogel: " << location+1 << " Seifert circles with the same and " << (num_circles-1-location) << " with the opposite orientation as the reference edge" << endl;

				/* The first location+1 Seifert circles have the same orientation as the reference edge, the remainder have the opposite orientation */
				turning_number = (location+1)*orientation_reference_polarity + (num_circles-1-location)*orientation_reference_polarity*-1;
			}
			else
			{
if (braid_control::VOGEL_DEBUG)
	debug << "vogel: " << location << "Seifert circles with the same and " << (num_circles-location) << " with the opposite orientation as the reference edge" << endl;

				/* The first location Seifert circles have the same orientation as the reference edge, the remainder have the opposite orientation */
				turning_number = location*orientation_reference_polarity + (num_circles-location)*orientation_reference_polarity*-1;		
			}
		}

if (braid_control::VOGEL_DEBUG)
	debug << "vogel: turning_number = " << turning_number << endl;
			
			(*turning_number_ptr) = turning_number;
	}


	/* next set up the strand numbers */	
	vector<int> strand(num_circles);
	
	for (int i=0; i<num_circles; i++)
    	strand[chain[i]-1] = i+1;

if (braid_control::VOGEL_DEBUG)
{	
	debug << "vogel: strand numbers: ";
    for (int i=0;i<num_circles;i++)
		debug << strand[i] << " ";
    debug << endl;
}


/* Now we are in a position to start reading off the braid word.
   We work along the chain successively evaluating the contribution
   of each child into its parent.
   
   The contributions to the braid word will be stored as strings
   attached to word.  word[i] is a pointer to the string that corresponds
   to circle i, NOT strand i.
*/
	
	char** word;
	if (!(word = new char*[num_circles]))
	{
	    vogel_error("couldn't allocate memory for word.\n");
	}
	for (int i = 0; i< num_circles; i++)
	    word[i] = NULL;

	/* evaluate an upper bound for the word lengths by assuming every crossing is negative
	   and involves the highest crossing number, we add to this the maximum length of Seifert
	   circles (as crossings) to allow for the '#' and '%' characters we will be
	   including as we evaluate the braid word
	*/
	int word_length = 0;
	for (int i=0; i< num_circles; i++)
	{
		if (seifert_circle[i][0] > word_length)
			word_length = seifert_circle[i][0];
	}
	word_length += num_crossings * 2; // to allow for "-s"

	ostringstream oss;
	oss << num_circles-1; // the highest crossing number
	word_length += num_crossings * oss.str().length();
	
	
if (braid_control::VOGEL_DEBUG)
    debug << "vogel: word length upper bound:" << word_length << endl;

	for (int i=0; i<num_circles; i++)
	{
		if (!(word[i] = new char[word_length+1]))
		{
    		vogel_error("couldn't allocate memory for word[i].\n");
		}

		int j;
		for (j=0; j<seifert_circle[i][0]; j++)
	    	*(word[i]+j) = '#';
		*(word[i]+j) = '\0';
	}

if (braid_control::VOGEL_DEBUG)
{
    debug << "vogel: initial circle words:\n";
    for (int i=0; i<num_circles; i++)
    {
		if (word[i])
		    debug << "vogel: " << i << ". vertex = " << i+1 << ": " << word[i] << endl;
    }
}

    /********* Start of word creation ********/

    char* cptr; //circle pointer
    char* pptr; //parent pointer
    char* xptr; //general pointer

    for (int i=0; i< num_circles-1; i++)
    {	
		int child = chain[i]-1;
	    int parent = chain[i+1]-1;    
	    int p_strand = strand[parent];

	    /* determine the number of characters we are going to insert into the parent word,
	       scan along the parent row in seifert_circle looking for crossings on the child 
	    */
	    for (int j=1; j<=seifert_circle[parent][0]; j++)
	    {
			int k=on_circle(seifert_circle,child, seifert_circle[parent][j]);
			if (k)
			{
				/* crossing seifert_circle[parent][j] is the kth crossing on the child */
				int num_chars = (crossing_type[seifert_circle[parent][j]] == generic_braid_data::crossing_type::NEGATIVE? 2 : 1);
				ostringstream oss;
				oss << p_strand-1;
				num_chars += oss.str().length();

		    	if (i>0)
		    	{
					/* there may be further word contributions from child's word, resulting from putting
					   the grandchild circle into child.  We look in the child word at the kth crossing, 
					   setting cptr to the character following the kth # or %.  It will be a # since the 
					   kth crossing on the child is on the parent.
					*/
					cptr = word[child];
					int count = 0;
					do
					{
			    		if (*cptr == '#' || *cptr == '%')
							count++;
			    		cptr++;
					} while (count<k);

					/* count from cptr until we reach another '#'.  Wrap at the end of the child word and 
					   skip any '%' characters.  Leave cptr where it is.  This gives the numnber of 
					   characters we have to copy into the parent from the child
					*/
					xptr = cptr;
					while (*xptr != '#')
					{
			    		if (*xptr == '\0')
			    		{
							xptr = word[child];
							if (*xptr == '#')
				    			break;
			    		}
			    		if (*xptr != '%')
							num_chars++;
			    		xptr++;
					}
		    	}

		    	/* Next make space in the parent word for the additional characaters.
		    	   Look for the j-th crossing in the parent
		    	*/
		    	xptr = word[parent];
		    	int count = 0;
		    	do
		    	{
					if (*xptr == '#' || *xptr == '%')
			    		count++;
					xptr++;
		    	} while (count<j);

		    	/* xptr has been advanced beyond the last #, so move back one place, then move 
		    	   the parent word to the right by num_chars, from the character following
			       the '#' at xptr.
			    */
		    	xptr--;

		    	/* set pptr to the NULL at the end of the word */
		    	pptr = word[parent]+strlen(word[parent]);
		    	int l=pptr-xptr; // the number of chars to move
		    	xptr = pptr+num_chars; // the new end point
		    	for (; l>0; l--)
					*xptr-- = *pptr--;

		    	/* now we can add the charcters to the parent word, we change the '#' to a '%' to indicate 
		    	   that the next characters came from a child child.  pptr currently indicates the '#'
			    */
			    *pptr++ = '%';

			    /* first add the crossing */
		    	if (crossing_type[seifert_circle[parent][j]] == generic_braid_data::crossing_type::NEGATIVE)
					*pptr++ = '-';
//		    	if (crossing_type[seifert_circle[parent][j]] == generic_braid_data::crossing_type::VIRTUAL)
		    	if (crossing_type[seifert_circle[parent][j]] == generic_braid_data::crossing_type::VIRTUAL || crossing_type[seifert_circle[parent][j]] == generic_braid_data::crossing_type::FLAT)
					*pptr++ = 't';
		    	else // generic_braid_data::crossing_type::POSITIVE or generic_braid_data::crossing_type::FLAT
					*pptr++ = 's';

				pptr += oss.str().copy(pptr,string::npos);
	               
			    if (i>0)
		    	{
					/* copy any other word contributions from child's word.  cptr is correctly placed to
					   check so just copy from  cptr into parent word until we reach a '#'.  As before we 
					   wrap at the child word  end and skip any '%' characters.
					*/
					while (*cptr != '#')
					{
			    		if (*cptr == '\0')
			    		{
							cptr = word[child];
							if (*cptr == '#')
							    break;
			    		}
			    		if (*cptr != '%')
							*pptr++ = *cptr;
			    		cptr++;
					}
		    	}
			}
	    }
	
if (braid_control::VOGEL_DEBUG)
{
    debug << "vogel::create_word: word " << parent+1 << " after adding " << child+1 << ": " << word[parent] << endl;
}
    }

    /********* End of word creation ********/

if (braid_control::VOGEL_DEBUG)
{
    debug << "vogel: words after creation" << endl;
    for (int i=0; i<num_circles; i++)
    {
		if (word[i])
	    	debug << "vogel: " << i << ". vertex = " << i+1 << ": " << word[i] << endl;
    }
}

	/* take out all of the '%' characters from the word and reduce before returning it */

	// clear the oss ready for the new string
	oss.str("");

	int tree_root = chain[num_circles-1];
	mark = word[tree_root-1];
	for(unsigned int i=0; i< strlen(word[tree_root-1]); i++)
	{
	    if (*mark != '%')
			oss << *mark;
	    mark++;
	}

if (braid_control::VOGEL_DEBUG)
	debug << "vogel: braid word before reduction: " << oss.str() << endl;

	string braid_word = oss.str();
	
	if (braid_control::REDUCE_BRAIDS)
		braid_word = braid_reduce(oss.str());

if (braid_control::VOGEL_DEBUG)
	debug << "vogel: braid word after reduction: " << braid_word << endl;

	/* delete the word strings */
	for (int i=0; i<num_circles; i++)
	{
	    if (word[i])
			delete[] word[i];
	}

	delete[] word;
	
	return braid_word;
}
