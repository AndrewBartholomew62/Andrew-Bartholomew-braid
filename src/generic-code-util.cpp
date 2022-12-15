/**************************************************************************
void read_immersion_code (generic_code_data& code_data, string input_string)
void write_immersion_code(ostream& s, generic_code_data& code_data)
void print_code_data(ostream& s, generic_code_data& code_data, string prefix)
void read_peer_code (generic_code_data& code_data, string input_string)
void write_peer_code(ostream& s, const generic_code_data& code_data, bool zig_zags, bool labelled)
void read_gauss_code (generic_code_data& code_data, string input_string)
void write_gauss_code(ostream& s, generic_code_data& code_data, bool OU_FORMAT)
void read_planar_diagram (generic_code_data& code_data, string input_string)
void write_planar_diagram(ostream& s, generic_code_data& code_data)
void read_code_data (generic_code_data& code_data, string input_string)
void write_code_data(ostream& s, generic_code_data& code_data)
string convert_gauss_code(string OU_gauss_code);
void renumber_peer_code(generic_code_data& code_data, vector<int> shift)
int remove_peer_code_component(generic_code_data& code_data, int component, vector<int>& component_flags)
int remove_virtual_components(generic_code_data& code_data, vector<int>& component_flags)
generic_code_data partition_peer_code(generic_code_data& code_data, vector<int>& component_flags)
void cycle_gauss_code_labels(generic_code_data& gauss_code_data, int edge)
void cycle_gauss_code_labels(matrix<int>& crossing_peers, vector<int>& num_component_edges, vector<int>& first_edge_on_component, int num_crossings,int num_components, int edge)
list<int> vertex_span (generic_code_data& code_data, int initial_crossing, vector<int>* exclude)
bool calculate_turning_cycles(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles)
bool realizable_code_data(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles)
bool valid_knotoid_input(generic_code_data& code_data)
string over_preferred_gauss_code(generic_code_data& code_data, bool unoriented)
vector<int> classical_gauss_data(generic_code_data& code_data)
int amalgamate_zig_zag_counts(int a, int b)
string read_dowker_code (string input_string)
**************************************************************************/
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <valarray>
#include <iomanip>
#include <list>

using namespace std;

extern ofstream     debug;
extern ofstream     output;
extern ifstream     input;

#include <util.h>
#include <quaternion-scalar.h>
#include <polynomial.h>
#include <matrix.h>
#include <generic-code.h>
#include <debug-control.h>
#include <gauss-orientation.h>
#include <reidemeister.h>


void read_immersion_code (generic_code_data& code_data, string input_string)
{
	code_data.type = generic_code_data::immersion_code;
	code_data.head = -1;

	if (input_string.find("K:") != string::npos)
		code_data.immersion = generic_code_data::character::KNOTOID;
	else if (input_string.find("L:") != string::npos)
		code_data.immersion = generic_code_data::character::LONG_KNOT;
	else
		code_data.immersion = generic_code_data::character::CLOSED;

	char* inbuf = c_string(input_string);
	int num_crossings = 0;
	char* cptr = strchr(inbuf,'/');
	cptr++;
	while (*cptr != '\0')
	{
		if (*cptr == '+' || *cptr == '-' || *cptr == '*')
			num_crossings++;

		cptr++;
	}	

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "read_immersion_code: num_crossings determined from labelled immersion code = " << num_crossings << endl;

	int num_edges = 2*num_crossings;
	code_data.num_crossings = num_crossings;
	code_data.num_components = 1;
	vector<int> num_component_edges(1);
	vector<int> first_edge_on_component(1);
	first_edge_on_component[0]=0;
	num_component_edges[0] = num_edges;
	
	matrix<int> code_table(CODE_TABLE_SIZE,num_crossings);
	
	for (int i=0; i< num_crossings; i++)
		code_table[COMPONENT][i] = 0;

	vector<int> perm(num_crossings);
	vector<int> type(num_crossings);
	
	cptr = strchr(inbuf,'(');
    while(*cptr != '/')
	{
		int count = 0;

		/* We write the crossing numbers in the order they appear in the immersion code 
		   into perm together with their corresponding type.  Then we unravel 
		   the permuatation writing the PERM and TYPE rows from INVERSE and LABEL.
		*/
		
		for (int i=0; i<num_crossings; i++)
			type[i] = generic_code_data::TYPE2;

		while (*cptr != ')')
		{		
			if (*cptr == '-')
			{
				type[count] = generic_code_data::TYPE1;
				cptr++;
			}
			else if (isdigit(*cptr))
			{
				char* mark = cptr;
				while (isdigit(*cptr)) cptr++;

				get_number (perm[count],mark);
				
				if (*cptr == '^')
				{
					/* we've found the knotoid head, note the 
					   corresponding crossing number */
					code_data.head = perm[count];
					code_data.immersion = generic_code_data::character::PURE_KNOTOID;
					cptr++;
				}
				
				count++;
			}
			else
				cptr++;
		}

		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "read_immersion_code: number of terms in this part of perm = " << count << endl;
	debug << "read_immersion_code: read absolute terms: ";
	for (int i=0; i< count; i++)
		debug << perm[i] << ' ';
	debug << endl;		
	debug << "read_immersion_code: signs: ";
	for (int i=0; i< count; i++)
		debug << type[i] << ' ';
	debug << endl;
}
		
		/* work out this part of the permutation */
		for (int i=0; i<count-1; i++)
		{
			code_table[OPEER][perm[i]] = (2*perm[i+1]-1+num_edges)%num_edges;
			code_table[EPEER][(code_table[OPEER][perm[i]]-1)/2] = 2*perm[i];
			code_table[TYPE][perm[i]] = type[i];
		}
		code_table[OPEER][perm[count-1]] = (2*perm[0]-1+num_edges)%num_edges;
		code_table[EPEER][(code_table[OPEER][perm[count-1]]-1)/2] = 2*perm[count-1];
		code_table[TYPE][perm[count-1]] = type[count-1];

		cptr++; // move over ')'
		while (*cptr == ' ') cptr++;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "read_immersion_code: interim code table produced from immersion code:" << endl;
	debug << "read_immersion_code:   type:    ";	
	for (int i=0; i< num_crossings; i++)
		debug << code_table[TYPE][i] << ' ';
	debug << endl;
	debug << "read_immersion_code:   OPEER: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_table[OPEER][i] << ' ';
	debug << endl;
	debug << "read_immersion_code:   EPEER: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_table[EPEER][i] << ' ';
	debug << endl;
	debug << "read_immersion_code:   label: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_table[LABEL][i] << ' ';
	debug << endl;
}

	}

	/* now do the labels */
	int count = 0;
	while (*cptr != '\0')
	{
		if (*cptr == '+')
			code_table[LABEL][count++] = generic_code_data::POSITIVE;
		else if (*cptr == '-')
			code_table[LABEL][count++] = generic_code_data::NEGATIVE;
		else if (*cptr == '*')
			code_table[LABEL][count++] = generic_code_data::VIRTUAL;
		cptr++;	
	}
	
	/* write the inverse permutation */
//	for (int i=0; i< num_crossings; i++)
//		code_table[INVERSE][code_table[PERM][i]] = i;

	/* now we can write the originating and terminating vertices and edges */
	vector<int> term_crossing(num_edges);
	vector<int> orig_crossing(num_edges);

/*	
	for (int i=0; i< num_crossings; i++)
	{
		term_crossing[2*i] = i;
		orig_crossing[2*i+1] = i;
		term_crossing[(2*code_table[PERM][i]-1 + num_edges) % num_edges] = i;
		orig_crossing[(2*code_table[PERM][i]) % num_edges] = i;
		
		code_table[EVEN_TERMINATING][i] = 2*i;
		code_table[ODD_ORIGINATING][i] = 2*i+1;
		code_table[ODD_TERMINATING][i] = (2*code_table[PERM][i]-1 + num_edges) % num_edges;
		code_table[EVEN_ORIGINATING][i] = (2*code_table[PERM][i]) % num_edges;
	}
*/
	for (int i=0; i< num_crossings; i++)
	{
		term_crossing[2*i] = i;
		orig_crossing[2*i+1] = i;
		term_crossing[code_table[OPEER][i]] = i;
		
		/* we need to identify the edge following the naming edge's peer */
		int component = code_table[COMPONENT][(code_table[OPEER][i]-1)/2];
		int peer_successor = (code_table[OPEER][i]+1 - first_edge_on_component[component])%
		                     num_component_edges[component] + first_edge_on_component[component];

		orig_crossing[peer_successor] = i;
		
		code_table[EVEN_TERMINATING][i] = 2*i;
		code_table[ODD_ORIGINATING][i] = 2*i+1;
		code_table[ODD_TERMINATING][i] = code_table[OPEER][i];
		code_table[EVEN_ORIGINATING][i] = peer_successor;
	}
	
	code_data.code_table = code_table;
	code_data.num_component_edges = num_component_edges;
	code_data.first_edge_on_component = first_edge_on_component;
	code_data.term_crossing = term_crossing;
	code_data.orig_crossing = orig_crossing;
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "read_immersion_code: code data produced from immersion code:" << endl;
	print_code_data(debug,code_data,"read_immersion_code: ");	
}

	delete[] inbuf;
}

/* This function only uses head, num_crossings and the OPEER, TYPE, and LABEL rows 
   of the code_table in the generic_code_data structure.
*/
void write_immersion_code(ostream& s, generic_code_data& code_data)
{
	matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;

	if (code_data.immersion == generic_code_data::character::KNOTOID)
		s << "K:";
	else if (code_data.immersion == generic_code_data::character::LONG_KNOT)
		s << "L:";
	
	vector<int> flag(num_crossings); // used to track which crossings have been written
	for (int i=0; i<num_crossings; i++)
		flag[i] = 1;
		
	/* work out the permutation from the generic code data */
	vector<int> perm(num_crossings);
	for (int i=0; i< num_crossings; i++)
		perm[i] = ((code_table[OPEER][i]+1)/2)%num_crossings;
	
	bool found;
	do
	{
		found = false;
		int crossing;
		int start;
		
		/* look for another starting place in rptr */
		for (int i=0; i<num_crossings; i++)
		{
			if (flag[i])
			{
				crossing = start = i;
				flag[i] = 0;
				found = true;
				break;
			}
		}

		if (found)
		{
			s << "(";
			if (code_table[TYPE][crossing] == generic_code_data::TYPE1)
				s << '-';
			s << crossing;

			int next_crossing;
			do
			{
				next_crossing = perm[crossing];
				s << " ";
				if (code_table[TYPE][next_crossing] == generic_code_data::TYPE1)
					s << '-';
				s << next_crossing;

				if (next_crossing == code_data.head)
					s << '^';
					
				flag[next_crossing] = 0;
				crossing = next_crossing;
			} while (perm[crossing] != start);
			s << ")";
		}
	} while (found);

	s << " / ";
							
	for (int i=0; i< num_crossings; i++)
	{
		if (code_table[LABEL][i] == generic_code_data::POSITIVE)
			s << "+ ";
		if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
			s << "- ";
		if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
			s << "* ";
	}	
}

void print_code_data(ostream& s, generic_code_data& code_data, string prefix)
{
	int num_crossings = code_data.num_crossings;
	
	s << prefix << "code_type = ";
	switch(code_data.type)
	{
		case (generic_code_data::immersion_code): s << "immersion_code"; break;
		case (generic_code_data::peer_code): s << "peer_code"; break;
		case (generic_code_data::gauss_code): s << "gauss_code"; break;
		default: s << "unknown";
	};	
	s << endl;
	
	s << prefix << "immersion (character) = ";
	switch (code_data.immersion)
	{
		case (generic_code_data::character::CLOSED): s << "CLOSED"; break;
		case (generic_code_data::character::LONG_KNOT): s << "LONG_KNOT"; break;
		case (generic_code_data::character::PURE_KNOTOID): s << "PURE_KNOTOID"; break;
		case (generic_code_data::character::KNOTOID): s << "KNOTOID"; break;		
		default: s << "unknown";
	};	
	s << endl;
	
	s << prefix << "head = " << code_data.head << endl;
	s << prefix << "head_zig_zag_count = " << code_data.head_zig_zag_count << endl;
	s << prefix << "num_crossings = " << num_crossings << endl;
	s << prefix << "num_components = " << code_data.num_components << endl;
	s << prefix << "type:      ";	
	matrix<int>& code_table = code_data.code_table;
	
	for (int i=0; i< num_crossings; i++)
		s << code_table[TYPE][i] << ' ';
	s << endl;
	
/*	if (code_data.type == generic_code_data::immersion_code)
	{
		s << prefix << "perm:      ";
		for (int i=0; i< num_crossings; i++)
			s << code_table[PERM][i] << ' ';
		s << endl;
		s << prefix << "inverse:   ";
		for (int i=0; i< num_crossings; i++)
			s << code_table[INVERSE][i] << ' ';
		s << endl;
	}
	else
*/	
	{
		s << prefix << "odd peer: ";
		for (int i=0; i< num_crossings; i++)
			s << code_table[OPEER][i] << ' ';
		s << endl;
		s << prefix << "even peer:  ";
		for (int i=0; i< num_crossings; i++)
			s << code_table[EPEER][i] << ' ';
		s << endl;
	}
	
	s << prefix << "even term: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[EVEN_TERMINATING][i] << ' ';
	s << endl;
	s << prefix << "odd term:  ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[ODD_TERMINATING][i] << ' ';
	s << endl;
	s << prefix << "even orig: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[EVEN_ORIGINATING][i] << ' ';
	s << endl;
	s << prefix << "odd orig:  ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[ODD_ORIGINATING][i] << ' ';
	s << endl;
	s << prefix << "label:     ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[LABEL][i] << ' ';
	s << endl;	
	s << prefix << "component: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[COMPONENT][i] << ' ';
	s << endl;	

//	if (code_data.type == generic_code_data::peer_code)
	{
		s << prefix << "num_component_edges: ";
		for (unsigned int i=0; i< code_data.num_component_edges.size(); i++)
			s << code_data.num_component_edges[i] << ' ';
		s << endl;
		s << prefix << "first_edge_on_component: ";
		for (unsigned int i=0; i< code_data.first_edge_on_component.size(); i++)
			s << code_data.first_edge_on_component[i] << ' ';
		s << endl;
	}
	
	s << prefix << "term_crossing: ";
	for (unsigned int i=0; i< code_data.term_crossing.size(); i++)
		s << code_data.term_crossing[i] << ' ';
	s << endl;
	s << prefix << "orig_crossing: ";
	for (unsigned int i=0; i< code_data.orig_crossing.size(); i++)
		s << code_data.orig_crossing[i] << ' ';
	s << endl;
	s << prefix << "shortcut_crossing: ";
	for (unsigned int i=0; i< code_data.shortcut_crossing.size(); i++)
		s << code_data.shortcut_crossing[i] << ' ';
	s << endl;	
	
	if (code_data.zig_zag_count.numrows() != 0)
	{
		s << prefix << "zig_zag_count: " << endl;
		print (code_data.zig_zag_count, s, 3, prefix);
	}
}

/* This function reads a peer code from the input_string into a generic_code_data object.  Within that
   object it creates a code_table that is a matrix<int>(9,num_crossings).  The function records in this
   matrix the crossing type, the peers of the even and odd edges 
   that terminate at the crossing, the component to which the naming edge belongs and the crossing label. 
   For each crossing it also records the odd and even terminating and originating vertices.
   It uses the #define TYPE, OPEER, EPEER, COMPONENT, LABEL, EVEN_TERMINATING, ODD_TERMINATING, 
   EVEN_ORIGINATING, ODD_ORIGINATING to index the rows of the matrix.
   If non-zero, the int* head is used to identify the first shortcut crossing in a knotoid (identified by 
   a '^' character after the crossing number).  This in turn identifies where in the peer code the head of 
   the knotoid is located. 
*/
void read_peer_code (generic_code_data& code_data, string input_string)
{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code: provided with input_string " << input_string << endl;
	
	code_data.type = generic_code_data::peer_code;
	code_data.head = -1;

	if (input_string.find("K:") != string::npos)
		code_data.immersion = generic_code_data::character::KNOTOID;
	else if (input_string.find("L:") != string::npos)
		code_data.immersion = generic_code_data::character::LONG_KNOT;
	else
		code_data.immersion = generic_code_data::character::CLOSED;

	
	char* inbuf = c_string(input_string);
	int num_crossings = 0;
	char* cptr = strchr(inbuf,'/');
	cptr++;
	while (*cptr != '\0')
	{
		if (*cptr == '+' || *cptr == '-' || *cptr == '*' || *cptr == '#' )
			num_crossings++;

		cptr++;
	}	

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code: num_crossings determined from labelled peer code = " << num_crossings << endl;
	
	code_data.num_crossings = num_crossings;
	
	matrix<int> code_table(CODE_TABLE_SIZE,num_crossings);

	int num_edges = 2*num_crossings;
	int component = 0;

	/* assume all crossings are generic_code_data::TYPE2 and set those 
	   indicated as generic_code_data::TYPE1 accordingly.
	*/
	for (int i=0; i<num_crossings; i++)
		code_table[TYPE][i] = generic_code_data::TYPE2;
		
	cptr = strchr(inbuf,'[');
	cptr++;
    for (int i=0; i< num_crossings; i++)
	{
		/* We write the odd-numbered peer of edge 2i in code_table[OPEER][i] and the 
		   even-numbered peer of edge 2i+1 in code_table[EPEER][i] for i=0,...,num_crossings-1
		*/
		bool crossing_complete = false;
		do
		{
			if (*cptr == ']')
			{
				cout << "\nError! Not enough peers specified in peer code" << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_peer_code: Error! Not enough peers specified in peer code" << endl;
					
				exit(0);
			}
			if (*cptr == '-')
			{
				code_table[TYPE][i] = generic_code_data::TYPE1;
				cptr++;
			}
			else if (isdigit(*cptr))
			{
				char* mark = cptr;
				while (isdigit(*cptr)) cptr++;

				get_number (code_table[OPEER][i],mark);
				
				if (code_table[OPEER][i]%2 == 0)
				{
					cout << "\nError! Read even number " << code_table[OPEER][i] << " in peer code" << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_peer_code: Error!  Read even number " << code_table[OPEER][i] << " in a peer code" << endl;
					
					exit(0);
				}
				
				code_table[EPEER][(code_table[OPEER][i]-1)/2] = 2*i;
				code_table[COMPONENT][i] = component;
				
				if (*cptr == '^')
				{
					/* we've found the knotoid head, note the 
					corresponding crossing number */
					code_data.head = i;
					code_data.immersion = generic_code_data::character::PURE_KNOTOID;
					cptr++;
				}			
				crossing_complete = true;
			}
			else if (*cptr == ',')
			{
				/* we've come to the end of a component and are about to 
				   read the first crossing of the new component
				*/
				component++;
				cptr++;
			}
			else
				cptr++;		
		} while (!crossing_complete);
	}
	
	/* now do the labels */
	int count = 0;
	while (*cptr != '\0')
	{
		if (*cptr == '+')
			code_table[LABEL][count++] = generic_code_data::POSITIVE;
		else if (*cptr == '-')
			code_table[LABEL][count++] = generic_code_data::NEGATIVE;
		else if (*cptr == '*')
			code_table[LABEL][count++] = generic_code_data::VIRTUAL;
		else if (*cptr == '#')
			code_table[LABEL][count++] = generic_code_data::FLAT;
		cptr++;	
	}

	/* write the component data into code_data */
	int num_components = component+1;
	code_data.num_components = num_components;

	vector<int> num_component_edges(num_components);
	vector<int> first_edge_on_component(num_components);
	first_edge_on_component[0]=0;
	component=0;
	
	for (int i=1; i< num_crossings; i++)
	{
		if (code_table[COMPONENT][i] != code_table[COMPONENT][i-1])
		{
			num_component_edges[component] = 2*i-first_edge_on_component[component];
			component++;
			first_edge_on_component[component] = 2*i;
		}
	}
	
	num_component_edges[component] = num_edges - first_edge_on_component[component];
	
	/* now we can write the originating and terminating vertices and edges */
	vector<int> term_crossing(num_edges);
	vector<int> orig_crossing(num_edges);
	
	for (int i=0; i< num_crossings; i++)
	{
		term_crossing[2*i] = i;
		orig_crossing[2*i+1] = i;
		term_crossing[code_table[OPEER][i]] = i;
		
		/* we need to identify the edge following the naming edge's peer */
		int component = code_table[COMPONENT][(code_table[OPEER][i]-1)/2];
		int peer_successor = (code_table[OPEER][i]+1 - first_edge_on_component[component])%
		                     num_component_edges[component] + first_edge_on_component[component];

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code: crossing " << i << " peer_successor = " << peer_successor << endl;

		orig_crossing[peer_successor] = i;
		
		code_table[EVEN_TERMINATING][i] = 2*i;
		code_table[ODD_ORIGINATING][i] = 2*i+1;
		code_table[ODD_TERMINATING][i] = code_table[OPEER][i];
		code_table[EVEN_ORIGINATING][i] = peer_successor;
	}
	
	code_data.code_table = code_table;
	code_data.num_component_edges = num_component_edges;
	code_data.first_edge_on_component = first_edge_on_component;
	code_data.term_crossing = term_crossing;
	code_data.orig_crossing = orig_crossing;
	
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "read_peer_code: code data produced from peer code:" << endl;
	print_code_data(debug,code_data,"read_peer_code: ");	
}

	delete[] inbuf;
}

/* This function uses head, num_crossings and the OPEER, TYPE, LABEL 
   and COMPONENT rows of the code_table in the generic_code_data structure.
*/
void write_peer_code(ostream& s, const generic_code_data& code_data, bool zig_zags, bool labelled)
{

	if (code_data.type != generic_code_data::peer_code)
	{
		cout << "Error!  write_peer_code called with generic_code_data that is not a peer code" << endl;
		exit (0);
	}

	const matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;
	
	if (code_data.immersion == generic_code_data::character::KNOTOID)
		s << "K:";
	else if (code_data.immersion == generic_code_data::character::LONG_KNOT)
		s << "L:";
	
	s << '[';
	
	if (code_table[TYPE][0] == generic_code_data::TYPE1)
		s << '-';
	if (code_data.num_crossings >0)
		s << code_table[OPEER][0];
		
	if (code_data.head == 0)
		s << "^ ";
	else if (num_crossings > 1 && code_table[COMPONENT][0] == code_table[COMPONENT][1])
		s << ' ';
		
	for (int i=1; i<num_crossings; i++)
	{
		if (code_table[COMPONENT][i] != code_table[COMPONENT][i-1])
			s << ", ";
			
		if (code_table[TYPE][i] == generic_code_data::TYPE1)
			s << '-';
		s << code_table[OPEER][i];
		if (code_data.head == i)
			s << '^';
			
		if ( i < num_crossings-1 && code_table[COMPONENT][i] == code_table[COMPONENT][i+1])
			s << ' ';
	}
	s << "]";
	if (labelled)
	{
		s << "/";
		for (int i=0; i< num_crossings; i++)
		{
			if (code_table[LABEL][i] == generic_code_data::POSITIVE)
				s << "+";
			else if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
				s << "-";
			else if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
				s << "*";
			else // (code_table[LABEL][i] == generic_code_data::FLAT)
				s << "#";
			
			if (i< num_crossings-1)
				s << " ";
		}	
	}
	
	if (zig_zags && code_data.zig_zag_count.numrows() != 0)
	{
		if (code_data.immersion == generic_code_data::character::KNOTOID || code_data.immersion == generic_code_data::character::LONG_KNOT)
		{
			s << " (";
			
			const matrix<int>& m = code_data.zig_zag_count;
			
			s << m[0][0] << ',' << code_data.head_zig_zag_count << ' ';
			
	        for (size_t j = 1; j < m.numcols(); j++)
	        {
				s << m[0][j] << ' ';
			}
			
		    for (size_t i = 1; i< m.numrows(); i++ )
		    {
		        for (size_t j = 0; j < m.numcols(); j++)
		        {
					s << m[i][j];
					if (i != m.numrows()-1 || j != m.numcols()-1)
						s << ' ';
				}
		    }
			s << ")";
		}
		else
		{
			bool save = matrix_control::SINGLE_LINE_OUTPUT;
			matrix_control::SINGLE_LINE_OUTPUT = true;
			s << " (" << code_data.zig_zag_count << ")";
			matrix_control::SINGLE_LINE_OUTPUT = save;
		}
	}
}

/* A Gauss code does not describe an immersion with crossing labels in the same way that a labelled peer code
   or a labelled immersion code does, since it does not include any information at all about virtual crossings.
   For this reason, if we number the semi-arcs that are described by a Gauss code then in the case of a virtual
   knot we do not always have an odd and even edge arriving at each crossing.  An example of this is the virtual 
   trefoil.  Therefore our use of generic code data for Gauss codes is different to that for peer codes and 
   immersion codes, and the conventions used here have been designed so that the bracket polynomial may be calculated
   for a virtual gauss code without changing the code developed for peer and immersion codes.
   
   A Gauss code for a classical link takes the form "1 -2 5 -4 3 -5, -1 2 -3 4 / ++--+", a Gauss code for a flat link,
   or a doodle has the form "L1 L2 R1 R2 L3 L4 R3 R4 / # # # #" (we distinguish this from the "Gauss orientation data"
   "L1 L2 R1 R2 L3 L4 R3 R4").  Hybrid Gauss codes such as "L1 2 R1 -2 3 -4 -3 4/#+--" are supported but are not used
   currently.  Virtual crossings are explicitly prohibited in Gauss codes, so it is not possible to describe the
   immersion of a virtual(-anything) diagram using the syntax of a Gauss code.  Labelled peer codes should be used
   for such diagrams.
   
   We trace the diagram (well, the classical crossings) described by the Gauss code, numbering semi-arcs without
   regard for any virtual crossings that may be required to realize the code.  We record the edge on which we arrive at 
   a crossing for the first time as the odd peer of the crossing and the edge on which we arrive at the crossing for 
   the second time as the even peer.
   
   The type of crossing is assigned based on the relative orientation of the first and second visits as follows:
           
    1st \ /              2nd \ /         
         X  =  Type 1         X  = Type 2    
    2nd / \              1st / \
   
   Where the orientation of all the semi-arcs is left to right.
   
   We convert the usual Gauss code label convention of identifying positive or negative crossings into the same
   convention used for immersion and peer codes, so that a crossing is labelled '+' if the "even" edge, that is the 
   edge on which we arrive at the crossing for the second time is the over-arc of the crossing, and '-' otherwise.  This
   allows the determination of crossing sign for generic code data to be unaffected.
      
   The odd and even originating and terminating edges are similarly assigned according to the first or second visit to the
   crossing, using the following convention:

    1st OT \ / OO             2nd ET \ / EO         
            X     =  Type 1           X  = Type 2    
    2nd ET / \ EO             1st OT / \ OO
   
   This allows component tracing as used by the bracket polynomial calculation to be unaffected.

   The component recorded for each crossing is the component associated with the edge on which we first arrive
   at the crossing.

*/
void read_gauss_code (generic_code_data& code_data, string input_string)
{
	if (input_string.find('O') != string::npos)
	{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "read_gauss_code: O detected in input string, converting Gauss code format." << endl;
		input_string = convert_gauss_code(input_string);
	}
	
	if (input_string.find('*') != string::npos)
	{
	    cout << "Error! Gauss codes may not describe crossings as virtual" << endl;
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "read_gauss_code: error in input string, '*' character found, indicating virtual crossing." << endl;
	
	    exit(0);
	}
	
	bool knotoid_or_long_indicator = false;
	
	if (input_string.find("K:") != string::npos)
	{
		code_data.immersion = generic_code_data::character::KNOTOID;
		knotoid_or_long_indicator = true;
	}
	else if (input_string.find("L:") != string::npos)
	{
		code_data.immersion = generic_code_data::character::LONG_KNOT;
		knotoid_or_long_indicator = true;
	}
	else
		code_data.immersion = generic_code_data::character::CLOSED;
	
	code_data.type = generic_code_data::gauss_code;
	code_data.head = -1;
//	code_data.immersion = generic_code_data::character::CLOSED;
	
	int num_components = (int) count(input_string.begin(),input_string.end(),',') + 1;

	char* inbuf = c_string(input_string);
	int num_crossings = 0;
	char* cptr = strchr(inbuf,'/');
	cptr++;
	while (*cptr != '\0')
	{
		if (*cptr == '+' || *cptr == '-' || *cptr == '#' )
			num_crossings++;

		cptr++;
	}	

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "read_gauss_code: num_crossings determined from gauss code " << input_string << " is " << num_crossings << endl;
	
	code_data.num_crossings = num_crossings;
	
	matrix<int> code_table(CODE_TABLE_SIZE,num_crossings);
	
	int num_edges = 2*num_crossings;
	int component = 0;

	/* 	write the OPEER row of the code table as -1 to act as flags	*/
	for (int i=0; i<num_crossings; i++)
		code_table[OPEER][i] = -1;
		
	cptr = inbuf;
	
	if (knotoid_or_long_indicator)
	{
		while (*cptr != ':')
			cptr++;
		cptr++; // step over the ':'
	}

	int edge = 0;
	vector<int> gauss_code_crossing(num_edges);
	vector<int> term_crossing(num_edges);
	vector<int> orig_crossing(num_edges);
	vector<int> num_component_edges(num_components);
	vector<int> first_edge_on_component(num_components);
	
	first_edge_on_component[0] = 0;
	
    for (int i=0; i< 2 * num_crossings; i++)
	{
		/* We write the edge on which we arrive at crossing i for the first time in code_table[OPEER][i] and
		   the edge on which we arrive at the crossing for the second time in code_table[EPEER][i].
		*/
		bool crossing_complete = false;
		int arc_type = generic_code_data::VIRTUAL; // initialized to this as we never use this value with Gauss codes
		do
		{
			if (*cptr == '-')
			{
				arc_type = generic_code_data::NEGATIVE;
				cptr++;
			}
			else if (*cptr == '+')
			{
				arc_type = generic_code_data::POSITIVE;
				cptr++;
			}
			else if (*cptr == 'L')
			{
				arc_type = generic_code_data::LEFT;
				cptr++;
			}
			else if (*cptr == 'R')
			{
				arc_type = generic_code_data::RIGHT;
				cptr++;
			}
			else if (isdigit(*cptr))
			{
				if (arc_type == generic_code_data::VIRTUAL) // i.e. not yet set
				    arc_type = generic_code_data::POSITIVE;
			
				char* mark = cptr;
				while (isdigit(*cptr)) cptr++;
				int crossing;
				get_number (crossing,mark);
				
				/* Gauss codes number crossings from one, we number them from zero */
				crossing--;
				
				gauss_code_crossing[i] = crossing;
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "read_gauss_code: crossings = " << crossing << ", edge = " << edge << endl;
	
				if (code_table[OPEER][crossing] == -1)
				{
					/* this is our first visit to this crossing */
					
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "read_gauss_code:   first visit, arc_type = ";
	switch (arc_type)
	{
		case generic_code_data::POSITIVE: 
			debug << "POSITIVE"; 
			break;
		case generic_code_data::NEGATIVE : 
			debug << "NEGATIVE"; 
			break;
		case generic_code_data::LEFT : 
			debug << "LEFT"; 
			break;
		case generic_code_data::RIGHT : 
			debug << "RIGHT";
	} 
	debug << endl;
}
					
					code_table[OPEER][crossing] = edge;						
					code_table[COMPONENT][crossing] = component;
					
					/* If we are at a classical crossing we can assign the crossing label as we now know whether 
					   the second visit (our "even" edge) will arrive on the over-arc or under-arc.  If the crossing 
					   is flat we can use the arc_type to set the TYPE of the crossing as well as the label.
					*/
					if (arc_type == generic_code_data::POSITIVE)
						code_table[LABEL][crossing] = generic_code_data::NEGATIVE;
					else if (arc_type == generic_code_data::NEGATIVE)
						code_table[LABEL][crossing] = generic_code_data::POSITIVE;
					else
					{
						code_table[LABEL][crossing] = generic_code_data::FLAT;
						
						if (arc_type == generic_code_data::LEFT)
							code_table[TYPE][crossing] = generic_code_data::TYPE2;
						else
							code_table[TYPE][crossing] = generic_code_data::TYPE1;
							
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "read_gauss_code:   set type = " << (code_table[TYPE][crossing] == generic_code_data::TYPE1? "TYPE1":"TYPE2") << endl;
					}
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "read_gauss_code:   second visit, arc_type = ";
	switch (arc_type)
	{
		case generic_code_data::POSITIVE: 
			debug << "POSITIVE"; 
			break;
		case generic_code_data::NEGATIVE : 
			debug << "NEGATIVE"; 
			break;
		case generic_code_data::LEFT : 
			debug << "LEFT"; 
			break;
		case generic_code_data::RIGHT : 
			debug << "RIGHT";
	} 
	debug << endl;
}
					/* this is the second visit to this crossing */
					code_table[EPEER][crossing] = edge;
				}
								
				edge++;
				crossing_complete = true;
			}
			else if (*cptr == ',')
			{
				/* we've come to the end of a component and are about to 
				   read the first crossing of the new component
				*/
				num_component_edges[component] = edge - first_edge_on_component[component];
				component++;
				first_edge_on_component[component] = edge;
				cptr++;
			}
			else
				cptr++;		
		} while (!crossing_complete);
	}
	
	/* finish off num_component_edges for the last component */
	num_component_edges[component] = edge - first_edge_on_component[component];

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "read_gauss_code: gauss_code_crossing: ";
	for (int i=0; i< num_edges; i++)
		debug << gauss_code_crossing[i] << ' ';
	debug << endl;
}

	/* now we can identify the TYPE of each clasical crossing from the sign of the 
	   crossing in the Gauss code and the LABEL assignment.  Flat crossing have
	   had their TYPE set already.
	*/
	int count = 0;
	while (*cptr != '\0')
	{
		if (*cptr == '+')
		{
			if (code_table[LABEL][count] == generic_code_data::POSITIVE)
				code_table[TYPE][count] = generic_code_data::TYPE2;
			else
				code_table[TYPE][count] = generic_code_data::TYPE1;
				
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "read_gauss_code:   set type for crossing " << count << ": positive crossing, second visit " << (code_table[LABEL][count] == generic_code_data::POSITIVE? "over": "under")
	      << " type is " << (code_table[TYPE][count] == generic_code_data::TYPE1? "TYPE1":"TYPE2") << endl;
}
			count++;
		}
		else if (*cptr == '-')
		{
			if (code_table[LABEL][count] == generic_code_data::POSITIVE)
				code_table[TYPE][count] = generic_code_data::TYPE1;
			else
				code_table[TYPE][count] = generic_code_data::TYPE2;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "read_gauss_code:   set type for crossing " << count << ": negative crossing, second visit " << (code_table[LABEL][count] == generic_code_data::POSITIVE? "over": "under")
	      << " type is " << (code_table[TYPE][count] == generic_code_data::TYPE1? "TYPE1":"TYPE2") << endl;
}

			count++;
		}
		else if (*cptr == '#')
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "read_gauss_code:   type for flat crossing " << count << " already set" << endl;
		    count++;
		}
		cptr++;	
	}
	
	/* identify the originating and terminating edges at each crossing.  We cannot do this
	   during our initial scan of the Gauss code as we need to know when we have reached the
	   end of each component and we don't know that until after we've processed the last crossing
	   on the component.
	*/
	vector<bool> first_visit_to_crossing(num_crossings);
	for (int i=0; i< num_crossings; i++)
		first_visit_to_crossing[i] = true;
	
	
	edge = 0;
	for (int i=0; i< num_components; i++)
	{
		for (int j=0; j< num_component_edges[i]; j++)
		{
			int crossing = gauss_code_crossing[edge++];
			if (first_visit_to_crossing[crossing])
			{
				code_table[ODD_TERMINATING][crossing] = first_edge_on_component[i] + j;
				code_table[EVEN_ORIGINATING][crossing] = first_edge_on_component[i] + (j+1)%num_component_edges[i];
				first_visit_to_crossing[crossing] = false;
			}
			else
			{				
				code_table[EVEN_TERMINATING][crossing] = first_edge_on_component[i] + j;
				code_table[ODD_ORIGINATING][crossing] = first_edge_on_component[i] + (j+1)%num_component_edges[i];
			}
			term_crossing[first_edge_on_component[i] + j] = crossing;
			orig_crossing[first_edge_on_component[i] + (j+1)%num_component_edges[i]] = crossing;
		}
	}
	
	/* write the data into code_data */
	code_data.num_components = num_components;	
	code_data.code_table = code_table;
	code_data.num_component_edges = num_component_edges;
	code_data.first_edge_on_component = first_edge_on_component;
	code_data.term_crossing = term_crossing;
	code_data.orig_crossing = orig_crossing;
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "read_gauss_code: code data produced from gauss code:" << endl;
	print_code_data(debug,code_data,"read_gauss_code: ");	
}

	delete[] inbuf;
}

/* write_gauss_code is capable of writing a gauss code from any kind of generic code data, 
   so it can convert from peer codes or immersion codes as well as output "native" Gauss codes
*/
void write_gauss_code(ostream& s, generic_code_data& code_data, bool OU_FORMAT)
{	
	int debug_save = debug_control::DEBUG;
	debug_control::DEBUG = debug_control::OFF;

	matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;
	int num_components = code_data.first_edge_on_component.size();

	if (code_data.type == generic_code_data::gauss_code)
	{
		/* we are writing out a "native" Gauss code, so the code_data follows Gauss code conventions */
		
		/* we have to follow the edges around the Gauss code noting each crossing
		   we encounter as we do so.  Since Gauss codes record the comonent of the first
		   visit in the COMPONENT row of the code table, we cannot use this row to output
		   the ',' between components accurately, since there may be two or more additional 
		   components in the diagram that are not recorded in that row.  Instead, we have to use
		   the edges recorded in first_edge_on_component to determine when we reach the end of a
		   component.
        */	  

		if (code_data.immersion == generic_code_data::character::KNOTOID)
			s << "K:";
		else if (code_data.immersion == generic_code_data::character::LONG_KNOT)
			s << "L:";

		bool add_space = false;

		for (int i=0; i < 2*num_crossings; i++)
		{
			bool last_edge_on_component = false;
			for (int j=0; j< num_components; j++)
			{
				if (i+1 == code_data.first_edge_on_component[j])
				{
					last_edge_on_component = true;
					break;
				}
			}

			/* locate the edge in the code table */
			for (int j=0; j< num_crossings; j++)
			{
				if (code_table[OPEER][j] == i)
				{
					if (!OU_FORMAT && add_space)
						s << ' ';
	
					if (code_table[LABEL][j] == generic_code_data::POSITIVE)
					{	
						if (OU_FORMAT)
							s << 'U'; // we're going under the crossing						
						else
							s << '-'; // we're going under the crossing
					}
					else if (OU_FORMAT)
					{
						s << 'O'; // we're going over the crossing						
					}

					s << j+1;
					
					if (OU_FORMAT)
					{
						if ( (code_table[LABEL][j] == generic_code_data::POSITIVE && code_table[TYPE][j] == generic_code_data::TYPE2) ||
						     (code_table[LABEL][j] == generic_code_data::NEGATIVE && code_table[TYPE][j] == generic_code_data::TYPE1)
						   )
							s << "+";
						else if ( (code_table[LABEL][j] == generic_code_data::POSITIVE && code_table[TYPE][j] == generic_code_data::TYPE1) ||
						          (code_table[LABEL][j] == generic_code_data::NEGATIVE && code_table[TYPE][j] == generic_code_data::TYPE2)
								)
							s << "-";
						else
							s << "!"; // shouldn't get here
					}
					break;
				}
				else if (code_table[EPEER][j] == i)
				{					
					if (!OU_FORMAT && add_space)
						s << ' ';
	
					if (code_table[LABEL][j] == generic_code_data::NEGATIVE)
					{	
						if (OU_FORMAT)
							s << 'U'; // we're going under the crossing						
						else
							s << '-'; // we're going under the crossing
					}
					else if (OU_FORMAT)
					{
						s << 'O'; // we're going over the crossing						
					}
					
					s << j+1;

					if (OU_FORMAT)
					{
						if ( (code_table[LABEL][j] == generic_code_data::POSITIVE && code_table[TYPE][j] == generic_code_data::TYPE2) ||
						     (code_table[LABEL][j] == generic_code_data::NEGATIVE && code_table[TYPE][j] == generic_code_data::TYPE1)
						   )
							s << "+";
						else if ( (code_table[LABEL][j] == generic_code_data::POSITIVE && code_table[TYPE][j] == generic_code_data::TYPE1) ||
						          (code_table[LABEL][j] == generic_code_data::NEGATIVE && code_table[TYPE][j] == generic_code_data::TYPE2)
								)
							s << "-";
						else
							s << "!"; // shouldn't get here
					}
					
					break;
				}
				add_space = true;
			}
			
			if (last_edge_on_component)
				s << ',';	
		}

		if (!OU_FORMAT)
		{
			s << '/';
			for (int i=0; i< num_crossings; i++)
			{
				if ( (code_table[LABEL][i] == generic_code_data::POSITIVE && code_table[TYPE][i] == generic_code_data::TYPE2) ||
				     (code_table[LABEL][i] == generic_code_data::NEGATIVE && code_table[TYPE][i] == generic_code_data::TYPE1)
				   )
					s << "+";
				else if ( (code_table[LABEL][i] == generic_code_data::POSITIVE && code_table[TYPE][i] == generic_code_data::TYPE1) ||
				          (code_table[LABEL][i] == generic_code_data::NEGATIVE && code_table[TYPE][i] == generic_code_data::TYPE2)
						)
					s << "-";
				else
					s << "! "; // shouldn't get here

				if (i< num_crossings-1)
					s << ' ';
					
			}	
		}
	}
	else
	{
		/* write Gauss code from immersion data in code_data */
		vector<int>& term_crossing = code_data.term_crossing;
		int num_edges = 2*num_crossings;
		vector<int> edge_flag(num_edges); // initialized to zero

		ostringstream oss;
		
		/* we need to re-number the crossings if we have virtual crossings in the immersion, 
		   since these are ignored by the Gauss code.
		*/
		int num_classical_crossings = num_crossings;
		bool pure_knotoid_code_data = false;
		vector<int>& shortcut_crossing = code_data.shortcut_crossing;
		
		for (int i=0; i< num_crossings; i++)
		{
			if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
				num_classical_crossings--;
		}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_gauss_code: num_classical_crossings = " << num_classical_crossings << endl;

		if (code_data.immersion == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && shortcut_crossing.size())
		{
			pure_knotoid_code_data = true;
			for (unsigned int i=0; i< shortcut_crossing.size(); i++)
			{
				if (shortcut_crossing[i] != 0)
					num_classical_crossings--;
			}
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "write_gauss_code: knotoid: num_classical_crossings = " << num_classical_crossings << endl;
		}
	
		int num_classical_crossings_visited = 0;
		
		/* classical_crossing will be used to renumber the immersion crossings,
		   so that if classical_crossing[i] = j, then crossing j of the immersion
		   is renumbered as crossing i in the Gauss code.
		   
		   crossing_visited will be a flag to indicate whether the ith immersion crossing
		   has been recorded in classical_crossing or not.
		*/
		vector<int> classical_crossing(num_classical_crossings);
		vector<int> crossing_visited(num_crossings);

		int start=0;
		int edge=0;
		bool complete = false;


		if (code_data.immersion == generic_code_data::character::KNOTOID || code_data.immersion == generic_code_data::character::PURE_KNOTOID)
			oss << "K:";
		else if (code_data.immersion == generic_code_data::character::LONG_KNOT)
			oss << "L:";
		
		bool add_space = false;
		
		do 
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_gauss_code: component_start = " << start << endl;
			
			/*	trace this component */
			do
			{	

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_gauss_code: edge = " << edge;
				edge_flag[edge] = 1;
				int next_crossing = term_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", next_crossing = " << next_crossing;
		
				if (code_table[LABEL][next_crossing] == generic_code_data::VIRTUAL)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", is virtual" << endl;
				}
				else if (pure_knotoid_code_data && shortcut_crossing[next_crossing])
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", is a shortcut crossing" << endl;
				}
				else if ((edge%2 != 0 && code_table[LABEL][next_crossing] == generic_code_data::POSITIVE) ||
				    (edge%2 == 0 && code_table[LABEL][next_crossing] == generic_code_data::NEGATIVE))
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", going under" << endl;
	
					if (!OU_FORMAT && add_space)
						oss << ' ';
									
					if (OU_FORMAT)
						oss << 'U'; // we're going under the crossing						
					else
						oss << '-'; // we're going under the crossing					
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", going over" << endl;

					if (OU_FORMAT)
						oss << 'O'; // we're going over the crossing						

					if (!OU_FORMAT && add_space)
						oss << ' ';					
				}
				
				if (code_table[LABEL][next_crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[next_crossing]))
				{
					if(crossing_visited[next_crossing])
					{
						for (int i=0; i< num_classical_crossings_visited; i++)
						{
							if (classical_crossing[i] == next_crossing)
							{
								oss << i+1; // Gauss crossings are numbered from 1 not zero

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_gauss_code:   second visit to Gauss crossing " << i+1<< endl;
							}
						}
					}
					else
					{

						classical_crossing[num_classical_crossings_visited] = next_crossing;
						crossing_visited[next_crossing] = 1;
						num_classical_crossings_visited++;
						
						oss << num_classical_crossings_visited;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_gauss_code:   first visit, becomes Gauss crossing " << num_classical_crossings_visited << endl;
					}

					if (OU_FORMAT)
					{
						if ( (code_table[LABEL][next_crossing] == generic_code_data::POSITIVE && code_table[TYPE][next_crossing] == generic_code_data::TYPE2) ||
						     (code_table[LABEL][next_crossing] == generic_code_data::NEGATIVE && code_table[TYPE][next_crossing] == generic_code_data::TYPE1)
						   )
							oss << "+";
						else if ( (code_table[LABEL][next_crossing] == generic_code_data::POSITIVE && code_table[TYPE][next_crossing] == generic_code_data::TYPE1) ||
						          (code_table[LABEL][next_crossing] == generic_code_data::NEGATIVE && code_table[TYPE][next_crossing] == generic_code_data::TYPE2)
								)
							oss << "-";
						else
							oss << "!"; // shouldn't get here
					}



					if (edge%2)
						edge = code_table[EVEN_ORIGINATING][next_crossing];
					else
						edge = code_table[ODD_ORIGINATING][next_crossing];
	
					add_space = true;
					
//					if (edge != start)
//						oss << ' ';

				}
				else
				{
					/* just move on around the component */				
					if (edge%2)
						edge = code_table[EVEN_ORIGINATING][next_crossing];
					else
						edge = code_table[ODD_ORIGINATING][next_crossing];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_gauss_code:   doing nothing" << endl;
				}				
			} while (edge != start);

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
			
			if (!complete)
				oss << ',';
				
		} while (!complete);


		if (!OU_FORMAT)
		{		
			oss << '/';
			
			for (int i=0; i< num_classical_crossings; i++)
			{
				int crossing = classical_crossing[i];
				
				if (code_table[LABEL][crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[crossing])) // it should never be!
				{
					if ((code_table[TYPE][crossing] == generic_code_data::TYPE1 && code_table[LABEL][crossing] == generic_code_data::NEGATIVE) ||
					    (code_table[TYPE][crossing] == generic_code_data::TYPE2 && code_table[LABEL][crossing] == generic_code_data::POSITIVE))
						oss << '+';
					else
						oss << '-';
				}
				
				if (i< num_classical_crossings-1)
					oss << ' ';
				
			}
		}
		
		s << oss.str();
	}
	debug_control::DEBUG = debug_save;
}

/* read_planar_diagram converts the PD data in input_string to a Gauss code and reads the resultant Gauss code into code_data */
void read_planar_diagram (generic_code_data& code_data, string input_string)
{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_planar_diagram: provided with input string: " << input_string << endl;
	
	int num_crossings = count(input_string.begin(),input_string.end(),'[');
	
	matrix<int> PD_data(num_crossings,4);

	istringstream iss(input_string);
	char c = 0;

	for (int i=0; i< num_crossings; i++)
	{
		while (c != '[')
			iss >> c;
		
		for (int j=0; j< 4; j++)
		{
			iss >> PD_data[i][j];
			iss >> c;
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_planar_diagram: PD_data: " << endl;
	print (PD_data,debug,3,"read_planar_diagram: ");
}
	
	/* planar diagram crossings are described starting at the ingress under_arc and working anti-clockwise around the crossing.
	   Thus, columns 0 and 2 of PD_data are the terminating and originating under-arcs respectively and columns 1 and 3 are the
	   over-arcs but we need to compare the values to determine which is the terminating and which the originating arc.
	*/
	vector<bool> edge_flag(2*num_crossings);  //used to record when an edge has been encountered
	vector<int> code_term(2*num_crossings);  // records the Gauss code terms in the order we determine them
	vector<int> num_component_edges;
	vector<int> start_of_component;
	vector<int> ingress_place(num_crossings);
			
	int num_components=0;
	int term_index=0;
	bool complete = false;
	
	int code_start_term=0;
	int code_start_component = 0;
	
	while (!complete)
	{
		complete = true;
		/* find the next start of a component. look first for an ingress under-arc that has not been considered */
		int start = 0;
		int crossing = -1;
		int place = -1;
		for (int j=0; j< num_crossings; j++)
		{
			if (edge_flag[PD_data[j][0]-1] == false)
			{
				start=PD_data[j][0];  // we know this is an ingress edge
				crossing = j;
				place = 0;
				complete = false;
				break;
			}
		}
		
		/* If we have a component that passes over every Gauss crossing that it encounters we need to compare the labels
		   on the two over-arcs to determine the orientation of the component.  If the componant has less than three Gauss
		   crossings, the planar diagram data doesn not contain sufficient information to determine the orientation unambiguously,
		   so we are free to choose either over semi-arc as the ingress.  If there are three or more crossings and the two semi-arc
		   labels are contiguous, then the smaller is the ingress edge and if they are not contiguous it is the greater that is the
		   ingress edge (and we are at the last crossing on the component).
		*/
		if (start == 0)
		{
			for (int j=0; j< num_crossings; j++)
			{
				if (edge_flag[PD_data[j][1]-1] == false)
				{
					int edge_1 = PD_data[j][1];
					int edge_2 = PD_data[j][3];

					/* the crossing sign will already have been set in this case */
					if (edge_2 == edge_1+1)
					{
						start=edge_1; 
						place = 1;
					}
					else if (edge_1 == edge_2+1)
					{
						start=edge_2; 
						place = 3;
					}
					else if (edge_1 < edge_2)
					{
						start=edge_2; 
						place = 3;
					}
					else // includes the case where edge_1 == edge_2
					{
						start=edge_1; 
						place = 1;
					}
						
					crossing = j;
					complete = false;
					break;
				}
			}
		}
		
		if(complete)
			break;
		
		start_of_component.push_back(term_index);
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_planar_diagram: start of next component at edge " << start << endl;
						
		int next_edge = start;
		int component_edge_count = 0;
		do 
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_planar_diagram:   next_edge = " << next_edge << ", crossing = " << crossing << ", place = " << place << endl;
	
			if (next_edge == 1 && crossing == 0)
			{
				code_start_term = term_index;
				code_start_component = num_components;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_planar_diagram:   code_start_term = " << code_start_term << ", code_start_component = " << code_start_component << endl;
			}
				
			edge_flag[next_edge-1] = true;
			component_edge_count++;

			if (place == 1 || place == 3)
				ingress_place[crossing] = place;
		
			code_term[term_index] = crossing+1;
			if (place%2 == 0)
				code_term[term_index] *= -1;  // we've arrived on the under_arc

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_planar_diagram:   code_term = " << code_term[term_index] << endl;
				
			next_edge = PD_data[crossing][(place+2)%4];
			term_index++;
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_planar_diagram:   succeeding edge = " << next_edge << endl;
			
			
			bool found = false;
			for (int j=0; j< num_crossings && !found; j++)
			{
				for (int k=0; k< 4; k++)
				{
					if (PD_data[j][k] == next_edge && !(j==crossing && k==(place+2)%4))
					{
						found = true;
						crossing = j;
						place = k;
						break;
					}
				}
			}					
			
			if (found)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_planar_diagram:   found succeeding edge at crossing " << crossing << ", place " << place << endl;
			}
			else
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_planar_diagram:   Error! failed to find succeeding edge";
			}
			
		} while(next_edge != start);		

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_planar_diagram: code_term: ";
	for (int j=0; j< 2*num_crossings; j++)
		debug << code_term[j] << ' ';
	debug << endl;
	debug << "read_planar_diagram: component_edge_count = " << component_edge_count << endl;
}
		num_component_edges.push_back(component_edge_count);
		num_components++;
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_planar_diagram: num_component_edges: ";
	for (int i=0; i< num_components; i++)
		debug << num_component_edges[i] << ' ';
	debug << endl;
	debug << "read_planar_diagram: ingress_place: ";
	for (int i=0; i< num_crossings; i++)
		debug << ingress_place[i] << ' ';
	debug << endl;
}

	vector<int> crossing_sign(num_crossings);
	
	for (int i=0; i< num_crossings; i++)
	{
		if (ingress_place[i] == 1)
			crossing_sign[i] = generic_code_data::NEGATIVE;
		else 
			crossing_sign[i] = generic_code_data::POSITIVE;			
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_planar_diagram: crossing_sign: ";
	for (int i=0; i< num_crossings; i++)
		debug << crossing_sign[i] << ' ';
	debug << endl;
	debug << "read_planar_diagram: write Gauss code" << endl;
}

	ostringstream oss;

	if (input_string.find("K:") != string::npos)
		oss << "K:";
	else if (input_string.find("L:") != string::npos)
		oss << "L:";

	for (int i=0; i< num_components; i++)
	{
		int component = (code_start_component+i)%num_components;
		int start;
		
		if (i==0)
			start = code_start_term;
		else
			start = start_of_component[component];

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_planar_diagram:   component = " << component << ", start = " << start << endl;
			
		for (int j=0; j< num_component_edges[component]; j++)
		{
			int term = (start - start_of_component[component] + j)%num_component_edges[component] + start_of_component[component];
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_planar_diagram:     j = " << j << ", term = " << term << ", code_term[" << term << "] = " << code_term[term] <<  endl;

			oss << code_term[term];
			if (j <num_component_edges[component] - 1)
				oss << ' ';
		}
		
		if (i < num_components-1)
			oss << ',';
	}	

	oss << '/';
	for (int i=0; i< num_crossings; i++)
	{
		if (crossing_sign[i] == generic_code_data::POSITIVE)
			oss << "+ ";
		else
			oss << "- ";
	}

	read_gauss_code(code_data,oss.str());
}

void write_planar_diagram(ostream& s, generic_code_data& code_data)
{	
	if (code_data.type != generic_code_data::code_type::gauss_code)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_planar_diagram: Error! presented with generic code data that is not a Gauss code." << endl;
	
		cout << "Error! presented with generic code data that is not a Gauss code." << endl;
		exit(0);
	}
	else
	{
		int num_crossings = code_data.num_crossings;
		
		matrix<int> PD_data(num_crossings,4); // planar diagram data
	
		matrix<int>& gauss_code_table = code_data.code_table;
		
		for (int i=0; i< num_crossings; i++)
		{
			if (gauss_code_table[LABEL][i] == generic_code_data::POSITIVE)
			{
				PD_data[i][0] = gauss_code_table[ODD_TERMINATING][i];
	
				if (gauss_code_table[TYPE][i] == generic_code_data::TYPE1)
				{
					PD_data[i][1] = gauss_code_table[EVEN_TERMINATING][i];
					PD_data[i][2] = gauss_code_table[EVEN_ORIGINATING][i];
					PD_data[i][3] = gauss_code_table[ODD_ORIGINATING][i];
				}
				else
				{
					PD_data[i][1] = gauss_code_table[ODD_ORIGINATING][i];
					PD_data[i][2] = gauss_code_table[EVEN_ORIGINATING][i];
					PD_data[i][3] = gauss_code_table[EVEN_TERMINATING][i];
				}			
			}
			else // generic_code_data::NEGAITIVE or generic_code_data::FLAT
			{
				PD_data[i][0] = gauss_code_table[EVEN_TERMINATING][i];
	
				if (gauss_code_table[TYPE][i] == generic_code_data::TYPE1)
				{
					PD_data[i][1] = gauss_code_table[EVEN_ORIGINATING][i];
					PD_data[i][2] = gauss_code_table[ODD_ORIGINATING][i];
					PD_data[i][3] = gauss_code_table[ODD_TERMINATING][i];
				}
				else
				{
					PD_data[i][1] = gauss_code_table[ODD_TERMINATING][i];
					PD_data[i][2] = gauss_code_table[ODD_ORIGINATING][i];
					PD_data[i][3] = gauss_code_table[EVEN_ORIGINATING][i];
				}			
			}
		}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "write_planar_diagram: PD_data:" << endl;
	print (PD_data,debug,3,"write_planar_diagram: ");
}		
		for (int i=0; i< num_crossings; i++)
		{
			s << "X[" << PD_data[i][0]+1 << ',' << PD_data[i][1]+1 << ',' << PD_data[i][2]+1 << ',' << PD_data[i][3]+1 << ']';
			if (i< num_crossings-1)
				s << ' ';
		}		
	}
}

void write_code_data(ostream& s, generic_code_data& code_data)
{
	if (code_data.type == generic_code_data::peer_code)
		write_peer_code(s, code_data);
	else if (code_data.type == generic_code_data::immersion_code)
		write_immersion_code(s,code_data);
	else
		write_gauss_code(s,code_data);
}

void read_code_data (generic_code_data& code_data, string input_string)
{
	if (input_string.find('X') != string::npos)
		read_planar_diagram(code_data, input_string);
	else if (input_string.find('[') != string::npos)
		read_peer_code(code_data, input_string);
	else if (input_string.find('(') != string::npos)
		read_immersion_code(code_data, input_string);
	else if (input_string.find("DT:") != string::npos)
	{
		string peer_code = read_dowker_code(input_string.substr(input_string.find(':')+1));
		read_peer_code(code_data, peer_code);
	}
	else
		read_gauss_code(code_data, input_string);
}

string convert_gauss_code(string OU_gauss_code)
{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "convert_gauss_code: presented with OU_gauss_code: " << OU_gauss_code << endl;			

	ostringstream oss;

	if (OU_gauss_code.find("K:") != string::npos)
		oss << "K:";
	else if (OU_gauss_code.find("L:") != string::npos)
		oss << "L:";

	
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
		if (code_data.code_table[LABEL][code_data.head] == generic_code_data::POSITIVE)
			head_semi_arc = code_data.code_table[OPEER][code_data.head];
		else if (code_data.code_table[LABEL][code_data.head] == generic_code_data::NEGATIVE)
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
		int component = code_table[COMPONENT][i];
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
				int position;
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
		int component = code_table[COMPONENT][i];
		int component_edges = num_component_edges[component];
		int first_edge = first_edge_on_component[component];
		int peer_component = code_table[COMPONENT][(code_table[OPEER][i]-1)/2];
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
				odd_edge = (edge_map[code_table[OPEER][i]]-peer_first_edge-1+peer_component_edges)%peer_component_edges + peer_first_edge;

				/* the new naming edge is on the same strand as the old naming edge */
				new_code_table[LABEL][even_edge/2] = code_table[LABEL][i];

				/* since the old naming edge is now odd, reversing both strands leaves the crossing type unchanged */
				new_code_table[TYPE][even_edge/2] = code_table[TYPE][i];


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
				odd_edge = edge_map[code_table[OPEER][i]];

				/* the new naming edge is on the same strand as the old naming edge */
				new_code_table[LABEL][even_edge/2] = code_table[LABEL][i];

				/* since the old naming edge is now odd, reversing the old naming strand reverses crossing type */
				if (code_table[TYPE][i] == generic_code_data::TYPE1)
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE2;
				else
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE1;
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
				even_edge = (edge_map[code_table[OPEER][i]]-peer_first_edge-1+peer_component_edges)%peer_component_edges + peer_first_edge;

				/* the new naming edge is on the opposite strand to the old naming edge */
				if (code_table[LABEL][i] == generic_code_data::POSITIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::NEGATIVE;
				else if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::POSITIVE;
				else if (code_table[LABEL][i] == generic_code_data::FLAT)
					new_code_table[LABEL][even_edge/2] = generic_code_data::FLAT;
				else
					new_code_table[LABEL][even_edge/2] = generic_code_data::VIRTUAL;

				/* since the old naming edge is now odd, reversing the peer of the old naming 
				   strand leaves the crossing type unchanged*/
				new_code_table[TYPE][even_edge/2] = code_table[TYPE][i];
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
				even_edge = edge_map[code_table[OPEER][i]];

				/* the new naming edge is on the opposite strand to the old naming edge */
				if (code_table[LABEL][i] == generic_code_data::POSITIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::NEGATIVE;
				else if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::POSITIVE;
				else if (code_table[LABEL][i] == generic_code_data::FLAT)
					new_code_table[LABEL][even_edge/2] = generic_code_data::FLAT;
				else
					new_code_table[LABEL][even_edge/2] = generic_code_data::VIRTUAL;

				/* since the old naming edge is now odd the crossing type is reversed */
				if (code_table[TYPE][i] == generic_code_data::TYPE1)
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE2;
				else
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE1;
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
	debug << "renumber_peer_code: code_table[LABEL][" << i << "] = " << code_table[LABEL][i] << endl;
	
	debug << "renumber_peer_code:   naming edge of old crossing " << i << " still even" << endl;
	debug << "renumber_peer_code:   both components at crossing " << i << " reversed" << endl;
	debug << "renumber_peer_code:   both type and label reversed for new crossing" << endl;
}

				odd_edge = (edge_map[2*i]-first_edge-1+component_edges)%component_edges + first_edge;
				even_edge = (edge_map[code_table[OPEER][i]]-peer_first_edge-1+peer_component_edges)%peer_component_edges + peer_first_edge;

				/* the new naming edge is on the opposite strand to the old naming edge */
				if (code_table[LABEL][i] == generic_code_data::POSITIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::NEGATIVE;
				else if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::POSITIVE;
				else if (code_table[LABEL][i] == generic_code_data::FLAT)
					new_code_table[LABEL][even_edge/2] = generic_code_data::FLAT;
				else
					new_code_table[LABEL][even_edge/2] = generic_code_data::VIRTUAL;

				/* since the old naming edge is still even, reversing both strands reverses the crossing type */
				if (code_table[TYPE][i] == generic_code_data::TYPE1)
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE2;
				else
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE1;
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
				even_edge = edge_map[code_table[OPEER][i]];

				/* the new naming edge is on the opposite strand to the old naming edge */
				if (code_table[LABEL][i] == generic_code_data::POSITIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::NEGATIVE;
				else if (code_table[LABEL][i] == generic_code_data::NEGATIVE)
					new_code_table[LABEL][even_edge/2] = generic_code_data::POSITIVE;
				else if (code_table[LABEL][i] == generic_code_data::FLAT)
					new_code_table[LABEL][even_edge/2] = generic_code_data::FLAT;
				else
					new_code_table[LABEL][even_edge/2] = generic_code_data::VIRTUAL;

				/* since the old naming edge is still even, reversing the old naming 
				   strand leaves the crossing type unchanged */
				new_code_table[TYPE][even_edge/2] = code_table[TYPE][i];
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
				odd_edge = (edge_map[code_table[OPEER][i]]-peer_first_edge-1+peer_component_edges)%peer_component_edges + peer_first_edge;

				/* the new naming edge is on the same strand as the old naming edge */
				new_code_table[LABEL][even_edge/2] = code_table[LABEL][i];

				/* since the old naming edge is still even, reversing the peer of the old naming 
				   strand reverses the crossing type */
				if (code_table[TYPE][i] == generic_code_data::TYPE1)
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE2;
				else
					new_code_table[TYPE][even_edge/2] = generic_code_data::TYPE1;
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
				odd_edge = edge_map[code_table[OPEER][i]];

				/* the new naming edge is on the same strand as the old naming edge */
				new_code_table[LABEL][even_edge/2] = code_table[LABEL][i];

				/* since the old naming edge is still even the crossing type is unchanged */
				new_code_table[TYPE][even_edge/2] = code_table[TYPE][i];
			}			
		}
		
		new_code_table[OPEER][even_edge/2] = odd_edge;
		new_code_table[EPEER][(odd_edge-1)/2] = even_edge;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code: renumbered as odd_edge = " << odd_edge << ", even_edge = " << even_edge << endl;
	debug << "renumber_peer_code: new type =  ";
	if (new_code_table[TYPE][even_edge/2] == generic_code_data::TYPE1)
		debug << "generic_code_data::TYPE1";
	else
		debug << "generic_code_data::TYPE2";
	
	debug << ", new label = ";
	if (new_code_table[LABEL][even_edge/2] == generic_code_data::POSITIVE)
		debug << "generic_code_data::POSITIVE";
	else if (new_code_table[LABEL][even_edge/2] == generic_code_data::NEGATIVE)
		debug << "generic_code_data::NEGATIVE";
	else if (new_code_table[LABEL][even_edge/2] == generic_code_data::VIRTUAL)
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
		term_crossing[new_code_table[OPEER][i]] = i;
		
		/* we need to identify the edge following the naming edge's peer */
		int component = new_code_table[COMPONENT][(new_code_table[OPEER][i]-1)/2];
		int peer_successor = (new_code_table[OPEER][i]+1 - first_edge_on_component[component])%
		                     num_component_edges[component] + first_edge_on_component[component];

		orig_crossing[peer_successor] = i;
		
		new_code_table[EVEN_TERMINATING][i] = 2*i;
		new_code_table[ODD_ORIGINATING][i] = 2*i+1;
		new_code_table[ODD_TERMINATING][i] = new_code_table[OPEER][i];
		new_code_table[EVEN_ORIGINATING][i] = peer_successor;
	}
	
	code_data.code_table = new_code_table;
	code_data.term_crossing = term_crossing;
	code_data.orig_crossing = orig_crossing;

	/* re-set the head, if the code data is for a knotoid.  Note that if the head_semi_arc is even then 
	   the head is already set correctly, even if the orientation of the first componenet is reversed
   */
	if (code_data.head !=-1 && head_semi_arc %2)
    {	
		code_data.head = code_data.code_table[EPEER][(head_semi_arc-1)/2]/2;
    }

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "renumber_peer_code: renumbered code data:" << endl;
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
			edge_flag[code_data.code_table[ODD_TERMINATING][crossing]] = 1;
			edge_flag[code_data.code_table[EVEN_TERMINATING][crossing]] = 1;
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
			if (code_data.code_table[LABEL][code_data.term_crossing[j]] != generic_code_data::VIRTUAL)
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
	if (code_data.immersion == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && code_data.shortcut_crossing.size())
	{
		pure_knotoid_code_data = true;
		
		if (code_data.code_table[LABEL][code_data.head] == generic_code_data::POSITIVE)
			head_semi_arc = code_data.code_table[OPEER][code_data.head];
		else if (code_data.code_table[LABEL][code_data.head] == generic_code_data::NEGATIVE)
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
			f_components[code_data.code_table[COMPONENT][*lptr]] = 1;
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
		   
		   We need the head, num_crossings and the OPEER, TYPE, LABEL and COMPONENT rows of the code_table in the generic_code_data 
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
			matrix<int> partition_code_table(CODE_TABLE_SIZE,num_partition_crossings,-1);
			matrix<int> partition_zig_zag_count(2,num_partition_crossings);
			vector<int> first_edge_on_partition_cpt(num_partition_components); // used for writing Gauss codes

			
			/* write new edge labels into the even (row 0)  and odd (row 1) rows of new_partition_labels, 
			   in the location corresponding to the old edge labels in the EVEN_TERMINATING and ODD_TERMINATING
			   rows of code_data.code_table.  In row 2 of new_partition_labels we record the component number
			   of the unicursal component in the partition.
			   
			   If we are working with a Gauss code, the third row of new_partition_labels records the crossing
			   number as determined by the new edge labels.
			*/
			matrix<int> new_partition_labels(4,code_data.num_crossings,-1);
			int edge = 0;
			int partition_crossing = 0;
			int component = -1;
			int new_head = -1;
			
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
						if (code_data.code_table[EVEN_TERMINATING][column] == old_edge)
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
							column = code_data.code_table[EPEER][(old_edge-1)/2]/2;
						}
					}
					
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "partition_peer_code:     old_edge " << old_edge << " is (row) " << (row == 0? "EVEN_TERMINATING": "ODD_TERMINATING") 
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
			
			/* write the new edge labels into the correct location in the EPEER and OPEER rows of the partition code table, 
			   and set the TYPE LABEL and COMPONENT rows 
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
						
					partition_code_table[EPEER][crossing] = new_partition_labels[0][j];
					partition_code_table[OPEER][crossing] = new_partition_labels[1][j];
					partition_code_table[COMPONENT][crossing] = new_partition_labels[2][j];
					partition_code_table[TYPE][crossing] = code_data.code_table[TYPE][j];
					partition_code_table[LABEL][crossing] = code_data.code_table[LABEL][j];
					
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
		debug << partition_code_table[OPEER][j] << ' ';
	debug << endl;
	debug << "partition_peer_code:  component: ";
	for (int j=0; j< num_partition_crossings; j++)
		debug << partition_code_table[COMPONENT][j] << ' ';
	debug << endl;
	debug << "partition_peer_code:  type: ";
	for (int j=0; j< num_partition_crossings; j++)
		debug << partition_code_table[TYPE][j] << ' ';
	debug << endl;
	debug << "partition_peer_code:  label: ";
	for (int j=0; j< num_partition_crossings; j++)
		debug << partition_code_table[LABEL][j] << ' ';
	debug << endl;
	
	if (track_zig_zag_counts)
	{
		debug << "partition_peer_code: partition_zig_zag_count: " << endl;
		print(partition_zig_zag_count,debug, 3,"partition_peer_code: ");	
	}
}
			generic_code_data partition_code_data;
			partition_code_data.type = code_data.type;
			partition_code_data.immersion = code_data.immersion;
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

				first_code_data.immersion = code_data.immersion;
				first_code_data.head_zig_zag_count = code_data.head_zig_zag_count;				
				
				if (pure_knotoid_code_data)
				{
					first_code_data.head = new_head;
					
					if (head_semi_arc % 2)
						first_code_data.head = first_code_data.code_table[EPEER][(head_semi_arc-1)/2]/2;
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
		opeer[i] = gauss_table[OPEER][i];
		epeer[i] = gauss_table[EPEER][i];
	}

	/* cycle the labels in opeer and epeer belonging to the identified component */
	for (int i=first; i <= last; i++)
	{
		/* find edge i in the code table */
		for (int j=0; j< num_crossings; j++)
		{
			if (gauss_table[EPEER][j] == i)
			{
				if (i == last)
					epeer[j] = first;
				else
					epeer[j]++;
					
				gauss_table[ODD_ORIGINATING][j] = (epeer[j]+1-first)%num_component_edges[component] + first;
				
				break;
			}
			else if (gauss_table[OPEER][j] == i)
			{
				if (i == last)
					opeer[j] = first;
				else
					opeer[j]++;
					
				gauss_table[EVEN_ORIGINATING][j] = (epeer[j]+1-first)%num_component_edges[component] + first;

				break;
			}
		}
	}

    for (int i=0; i< num_crossings; i++)
    {
		gauss_table[OPEER][i] = opeer[i];
		gauss_table[ODD_TERMINATING][i] = opeer[i];
		gauss_table[EPEER][i] = epeer[i];
		gauss_table[EVEN_TERMINATING][i] = epeer[i];
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

/* The function vertex_span counts the number of vertices that can be reached from the 
   initial crossing in the given generic code data.  This can be used to determine whether the
   generic code data can be realized by a connected immersion.  Optionally, the exclude vector
   may specify a number of edges to exclude from the search, this may be use to determine whether 
   the complement of the excluded edges is connected.
*/
list<int> vertex_span (generic_code_data& code_data, int initial_crossing, vector<int>* exclude)
{
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
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "vertex_span: checking code ";
	write_code_data (debug, code_data);
	debug << "\nvertex_span: starting from crossing " << initial_crossing;
}

	if (exclude != 0 && exclude->size())
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << " but excluding edges ";
    for (unsigned int i=0; i< exclude->size(); i++)
		debug << (*exclude)[i] << ' ';
    debug << endl;
}
		for (unsigned int i=0; i< exclude->size(); i++)
			edge_flag[(*exclude)[i]] = 1;
	}
	else
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << " with no excluded edges " << endl;
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

		int edge = code_table[EVEN_TERMINATING][*lptr];
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
		
		edge = code_table[ODD_TERMINATING][*lptr];
		
		if (edge_flag[edge] == 0)
		{
			crossing = orig_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span:     going back along " << code_table[ODD_TERMINATING][*lptr] << " takes us to " << crossing;

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
		
		edge = code_table[EVEN_ORIGINATING][*lptr];

		if (edge_flag[edge] == 0)
		{
			crossing = term_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span:     going forwards along " << code_table[EVEN_ORIGINATING][*lptr] << " takes us to " << crossing;

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
		
		edge = code_table[ODD_ORIGINATING][*lptr];
		if (edge_flag[edge] == 0)
		{
			crossing = term_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "vertex_span:     going forwards along " << code_table[ODD_ORIGINATING][*lptr] << " takes us to " << crossing;

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
//				if (cycle[j][k] == i)
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
//			int edge = i;
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
//					int vertex = orig_crossing[edge];
					if (code_table[TYPE][vertex] == generic_code_data::TYPE1)
		    			edge = code_table[EVEN_ORIGINATING][vertex];
					else
						edge = -code_table[ODD_TERMINATING][vertex];
//						edge = code_table[ODD_TERMINATING][vertex];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;

				}
				else // edge is even (and positive)
				{
					int vertex = term_crossing[edge];
					if (code_table[TYPE][vertex] == generic_code_data::TYPE1)
						edge = -code_table[ODD_TERMINATING][vertex];
//						edge = code_table[ODD_TERMINATING][vertex];
					else
		    			edge = code_table[EVEN_ORIGINATING][vertex];

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
		
if (debug_control::DEBUG >= debug_control::DETAIL)
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
//				if (cycle[j][k] == i)
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
//			int edge = i;
			cycle [num_cycles][column++] = edge;
			bool complete = false;

			do
			{
				/* determine which vertex the current edge takes us to and 
	   			the next edge we turn onto */
				if (edge % 2) // edge is odd (and negative)
				{
					int vertex = orig_crossing[-edge];
//					int vertex = orig_crossing[edge];
					if (code_table[TYPE][vertex] == generic_code_data::TYPE1)
						edge = -code_table[ODD_TERMINATING][vertex];
//						edge = code_table[ODD_TERMINATING][vertex];
					else
		    			edge = code_table[EVEN_ORIGINATING][vertex];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;

				}
				else // edge is even (and positive)
				{
					int vertex = term_crossing[edge];
					if (code_table[TYPE][vertex] == generic_code_data::TYPE1)
		    			edge = code_table[EVEN_ORIGINATING][vertex];
					else
						edge = -code_table[ODD_TERMINATING][vertex];
//						edge = code_table[ODD_TERMINATING][vertex];

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

if (debug_control::DEBUG >= debug_control::DETAIL)
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

/* in the case of links it is possible that the code data we have been given is disconnected, as in 
   [-3 -5 -1, -9 -11 7]/# # # # # #, therefore we first check that the code is connected, then check
   the Euler characteristic.  In the process of this check we return a set of turning cycles to the call.
*/
bool realizable_code_data(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles)
{
	int num_crossings = code_data.num_crossings;
	int num_edges = 2* num_crossings;
	
	/* first check that the code is connected by calling vertex_span starting at crossing 0 
	   and not excluding any edges.  If we cannot reach every crossing, the code cannot be
	   realized by a connected immersion.
	*/
	int crossing_span = vertex_span(code_data,0).size();
	if (crossing_span != num_crossings)
	{		
if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "realizable_code_data: vertex span is only " << crossing_span
         << " crossings, code cannot be realized by a connected immersion" << endl;
}   
		return false;
	}
	else
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "realizable_code_data: vertex span includes all the crossings, code can be realized by a connected immersion" << endl;
	}
	
	/* now calculate the turning cycles, this includes additional checks on the realizable nature of the code */	
	if (!calculate_turning_cycles(code_data, cycle, num_left_cycles, num_cycles))
		return false;

	/* next check each edge appears once in a left handed
	   and once in a right handed cycle */
	   
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "realizable_code_data: check edges appear exactly once in left and right turning cycles" << endl;

	bool realizable = true;
	for (int i=0; i<num_edges; i++)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "realizable_code_data:   edge " << i;
    
		int count = 0;
		for (int j=0; j<num_left_cycles; j++)
		{
			for (int k=1; k<= cycle[j][0]; k++)
			{
				if (abs(cycle[j][k]) == i)
				{
					if (++count == 2)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << " appears twice in left turning cycle " << j << endl;
						realizable = false;
						break;
					}
				}
			}
			if (!realizable)
				break;
		}
		
		if (!count)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << " does not appear at all in left turning cycles";
			realizable = false;
		}

		if (!realizable)
		{
			break;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << " appears exactly once in left turning cycles";
		}

		/* we're only still here if count = 1, so now check the right
		   turning cycles */

		count = 0;
		for (int j=num_left_cycles; j<num_cycles; j++)
		{
			for (int k=1; k<= cycle[j][0]; k++)
			{
				if (abs(cycle[j][k]) == i)
				{
					if (++count == 2)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << " but twice in right turning cycle " << j << endl;
						realizable = false;
						break;
					}
				}
			}
			if (!realizable)
				break;
		}
		
		if (!count)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << " but not at all in right turning cycles";
			realizable = false;
		}

		if (!realizable)
		{
			break;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << " and in right turning cycles" << endl;
		}
	}

	/* Now check that the number of cycles is equal to num_crossings+2 */
	if (num_cycles != num_crossings+2)
	{
		realizable = false;

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "realizable_code_data: number of cycles != num_crossings+2" << endl;
	}
	
	if (realizable)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "realizable_code_data: code is able to be realized by a connected immersion" << endl;
	}
	else
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "realizable_code_data: code cannot be realized by a connected immersion" << endl;
	}
	
	return realizable;
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
      
   As we check the validity of the input we note the shortcut crossings and also
   the algebraic number of intersections of the knotoid and the shortcut.  This is 
   the number of times the knotoid crosses the shortcut from right to left minus 
   the number of times it crosses from left to right.  However, for the definition
   the shortcut is oriented from the leg to the head and we will be traversing it in 
   the opposite direction.  
   
   We set shortcut_crossing[i] to be non-zero if the ith crossing is a shortcut crossing.
   We set it to 1 if the knotoid crosses the shortcut from right to left and -1 if it crosses from 
   left to right, orienting the shortcut from leg to head.
   
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
	if (code_table[LABEL][head] == generic_code_data::POSITIVE)
	{
		semi_arc = code_table[OPEER][head];
		peer = 2*head;
		component = code_table[COMPONENT][(semi_arc-1)/2];
	}
	else if (code_table[LABEL][head] == generic_code_data::NEGATIVE)
	{
		semi_arc = 2*head;
		peer = code_table[OPEER][head];
		component = code_table[COMPONENT][head];
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
		if (code_table[TYPE][head] == generic_code_data::TYPE1)				
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
		if (code_table[TYPE][head] == generic_code_data::TYPE1)				
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
			peer = code_table[EVEN_TERMINATING][crossing];

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
			else if (code_table[LABEL][crossing] != generic_code_data::POSITIVE)
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
				if (code_table[TYPE][crossing] == generic_code_data::TYPE1)				
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
			peer = code_table[ODD_TERMINATING][crossing];

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
			else if (code_table[LABEL][crossing] != generic_code_data::NEGATIVE)
			{
				valid = false;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "valid_knotoid_input:   not negative, as required for valid knotoid" << endl;
				break;
			}
			else
			{
			    /* see comment above for the odd case */
				if (code_table[TYPE][crossing] == generic_code_data::TYPE1)				
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
	int min_component = (g.immersion == generic_code_data::character::CLOSED? 0: 1);

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
			if (code_table[OPEER][i] >= g.first_edge_on_component[j])
				component++;
			else
				break;
		}
		first_component[i] = component;
	
		component = 0;
		for (int j=1; j< num_components; j++)
		{
			if (code_table[EPEER][i] >= g.first_edge_on_component[j])
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
	debug << "over_preferred_gauss_code:         old_edge = " << old_edge << " code_table[OPEER][terminating_crossing] = " << code_table[OPEER][terminating_crossing] << " code_table[EPEER][terminating_crossing] = " << code_table[EPEER][terminating_crossing] << endl;
}

						if (code_table[OPEER][terminating_crossing] == old_edge)
						{	
							if (code_table[LABEL][terminating_crossing] == generic_code_data::POSITIVE)
								over_arc = false;
						}
						else if (code_table[EPEER][terminating_crossing] == old_edge)
						{					
							if (code_table[LABEL][terminating_crossing] == generic_code_data::NEGATIVE)
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

						if (code_table[OPEER][terminating_crossing] == old_edge)
						{	
							if (code_table[LABEL][terminating_crossing] == generic_code_data::POSITIVE)
								over_arc = false;
						}
						else if (code_table[EPEER][terminating_crossing] == old_edge)
						{					
							if (code_table[LABEL][terminating_crossing] == generic_code_data::NEGATIVE)
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
									
					if ( (code_table[LABEL][inv_perm[i]] == generic_code_data::POSITIVE && code_table[TYPE][inv_perm[i]] == generic_code_data::TYPE2) ||
					     (code_table[LABEL][inv_perm[i]] == generic_code_data::NEGATIVE && code_table[TYPE][inv_perm[i]] == generic_code_data::TYPE1)
					   )
					{
						if (reverse_sign)
							perm_sign[i] = '-';
						else
							perm_sign[i] = '+';
					}
					else if ( (code_table[LABEL][inv_perm[i]] == generic_code_data::POSITIVE && code_table[TYPE][inv_perm[i]] == generic_code_data::TYPE1) ||
					          (code_table[LABEL][inv_perm[i]] == generic_code_data::NEGATIVE && code_table[TYPE][inv_perm[i]] == generic_code_data::TYPE2)
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
	int place = 0;
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

			if (code_table[OPEER][terminating_crossing] == old_edge)
			{	
				if (code_table[LABEL][terminating_crossing] == generic_code_data::POSITIVE)
					over_arc = false;
			}
			else if (code_table[EPEER][terminating_crossing] == old_edge)
			{					
				if (code_table[LABEL][terminating_crossing] == generic_code_data::NEGATIVE)
					over_arc = false;
			}
			
			if (!over_arc)
				oss << '-';
			
			oss << min_perm[terminating_crossing]+1;
			
			if (i < num_component_terms-1)
				oss << ' ';
				
			place++;
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


/* read_dowker_code creates a labelled peer code string corresponding to the Dowker-Thistlethwaite code of a prime knot, based on the conventions used by knotscape, where 
   the first Dowker-Thistlethwaite label (label 1) is assigned to an over-arc. That means that for alternating knots all the terms of the Dowker-Thistlethwaite code are 
   positive, since the even labels are always allocated to under-arcs at a crossing.  Since Dowker-Thistlethwaite codes are invariant under reflection in a line in the 
   plane disjoint from the diagram, knotscape fixes the TYPE of the crossing involving label 1 and evaluates the other crossings' TYPE using an algorithm based on 
   chess-board colouring that is reproduced here.
   
   Dowker-Thistlethwaite codes as used by knotscape allocate odd and even labels to shadow crossings but for the purposes of evaluating regions (i.e turning cycles)
   the labels are considered to be allocated to the originating semi-arcs at the crossing, rather than at the terminating semi-arcs as in the case of a peer code.

   The Dowker-Thistlethwaite code specifies the even label corresponding to odd labels 1, 3, 5,... in a sequence of terms that we consider numbered from zero; that is the i-th
   term is the even label associated with odd label 2*1+1.  Thus, we "decrement and shift backwards by one"; that is, subtract one from each Dowker-Thistlethwaite label 
   corresponding to the i-th term of a Dowker-Thistlethwaite code and consider the label to number the preceeding terminating edge at that crossing.  Thus, we obtain valid 
   peer-code labels for crossing i as numbered by labelled peer codes.   
*/
string read_dowker_code (string input_string)
{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_dowker_code: presented with input_string " << input_string << endl;
	
	istringstream iss(input_string);
	vector<int> dowker_code;
	int term;
	int num_crossings = 0;
	
	while(iss >> term)
	{
		dowker_code.push_back(term);
		num_crossings++;
	}
	
	int num_edges = 2*num_crossings;
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_dowker_code: num_crossings = " << num_crossings << ", num_edges = " << num_edges << endl;
	debug << "read_dowker_code: dowker_code: ";
	for (int i=0; i< num_crossings; i++)
		debug << dowker_code[i] << ' ';
	debug << endl;
}	

	/* Knowing the dowker code we may identify the peer of each edge, which we store in a single vector, since we are not building generic code data */
	vector<int> peer(num_edges);
	for (int i=0; i< num_crossings; i++)
	{
		peer[2*i] = abs(dowker_code[i])-1;
		peer[peer[2*i]] = 2*i;
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_dowker_code: peer: ";
	for (int i=0; i< num_edges; i++)
		debug << peer[i] << ' ';
	debug << endl;
}	

	/* Dowker-Thistlethwaite codes number semi-arcs from 1 and record the even peers of odd edges, so if a dowker_code term is negative
	   it indicates that the even dowker-edge is an over-arc.  That in turn means that the corresponding odd peer-label is an over-arc and 
	   therefore that the even peer code naming-edge is an under-arc.
	*/
	vector<int> crossing_label(num_crossings);
	for (int i=0; i< num_crossings; i++)
	{
		if (dowker_code[i] < 0)
		{
			crossing_label[i] = generic_code_data::label::NEGATIVE;
			dowker_code[i] *= -1;
		}
		else
		{
			crossing_label[i] = generic_code_data::label::POSITIVE;
		}
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_dowker_code: crossing_label: ";
	for (int i=0; i< num_crossings; i++)
		debug << crossing_label[i] << ' ';
	debug << endl;
}	

	/* To determine the crossing types we follow the approach used by knotscape where we successively consider sub-paths in the underlying immersion that trace the 
	   diagram from a 'base point' crossing back to the same crossing.  We leave the base point crossing on the odd originating edge and return on the odd terminating 
	   edge.
	   	   
	   The sub-path corresponding to a base point crossing may be a simple closed curve or may have self intersections.  However, if we consider the sub-diagram 
	   comprised of the sub-path and its self intersections only, ignoring any crossings with edges not in the sub-path, we have a diagram that may be chess-board 
	   coloured black and white. 
	   
	   With respect to this chess-board colouring, tracing around the sub-path from the base point we create a record of the whether the region to the right of the 
	   edge is black or white, which we refer to as the 'polarity' of the edge.  For the 'remaining' edges of the original diagram we can then identify whether they 
	   lie within a black or white region of the coloured sub-diagram and, again, refer to this as the polarity of the edge.  This yields four possibilities as 
	   indicated below, where the vertical edges represent (part of) an edge of the sub-path and the horizontal edge represents a pair of edges in the original diagram 
	   that intersect the vertical edge in a crossing of the original diagram that is not a self-intersection of the sub-path.  From polarity shown at the arrow heads 
	   it can be seen that for the ingress edges at the central crossing, the vertical edges lie adjacent in a clockwise direction around the crossing from the horizontal 
	   edge when the ingress horizontal and vertical polarity are the same (cases b and d) and anti-clockwise if they differ (cases a and c).  
	        

            a)                         b)                        c)                       d)
	                +                          +                         -                        -
	            B   ^   W                   B  ^  W                   W  ^  B                  W  ^  B
	                |                          |                         |                        |
	                |                          |                         |                        |
	        - ----->+------> +         - <-----+<----- +         + ----->+------> -       + <-----+<----- -  
	                ^                          ^                         ^                        ^
	                |                          |                         |                        |
	            B   |   W                   B  |  W                   W  |  B                  W  |  B
	   
	                +                          +                         -                        -

       The clockwise or anti-clockwise relationship between the ingress edges is used to set the crossing TYPE relative to the TYPE of the base point crossing.  Since the 
       Dowker-Thistlethwaite code is invariant under reflection in a line of the plane, we may start by setting the TYPE of crossing zero arbitrarily, ase described in more 
       detail below.  Note that in the above diagrams, for a sub-path edge, + indicates a black region lies to the left and - indicates white lies to the left.  For the remaining 
       edges + indicates that the edge lies in a white region and - that it lies in a black region.   We could have considered the colours to be reversed (+ meaning white to 
       the left of the subpath etc.), we would still have different cases describing the configuration of the crossing relative to the initial conditions.  In practice the
       convention used is determined by the base point crossing, as described below.

	   We record the polarity for each edge of the original diagram, so that successive edges in the subpath between self-intersections of the sub-path have the same 
	   polarity.  That is, the polarity of successive edges in the sub-path changes only if we encounter a crossing that is a self intersection of the sub-path.  Similarly
	   the polarity of successive remaining edges does not change until we encounter a crossing whose transverse strand lies in the sub-path.
	   
	   On arriving at a crossing on an ingress edge, we identify the transverse strand that we encounter as being part of the sub-path by testing whether the peer ingress 
	   edge lies in the sub-path.  That means (as can be seen from the diagrams below) that when we arrive back at the base point crossing at the end of the sub-path, given 
	   that the peer ingress edge is not part of the sub-path the first remaining edge will have the same polarity as the last sub-path edge.  Moreover, since we set the 
	   polarity of the sub_path_start_edge to + to initialize the algorithm, we do not test the peer of the even terminating edge at the base-point crossing that preceeds the 
	   sub_path_start_edge.  That means that terminating edge and the sub_path_start edge are assigned the same polarity, though the meaning of that polarity is different for 
	   these edges, as described above.
	   	   
       As noted above, at the base-point crossing, we consider the sub-path starting on the odd originating edge and terminating on the odd terminating edge.  We assume the 
       quadrant bounded by the start and end edges is coloured C, as shown below and initialize the polarity of the sub-path start edge as + (the upper vertical edge in the 
       diagrams below), establishing the polarity convention for this sub-path.  Clearly, the end edge of the sub_path_end_edge will have the same + polarity (the ingress 
       horizontal edges).
                           
	                +                                          +          
	                ^ start                                    ^  start
	            C   |                                          |  C       
	                |                                          |          
	   2j-1 + ----->+------> +                         + <-----+<----- + 2j-1 
	        end     ^                                          ^       end    
	                |                                          |          
	                |                                          |          
                    +  2i                                  2i  +    

              
       The base-point diagrams above show a TYPE1 crossing on the left and a TYPE2 crossing on the right.  We set the type of crossing zero arbitrarily, thereby determining 
       which of the mirror images descibed by the Dowker-Thistlethwaite code the resulting peer code will describe, and start setting the other crossings' TYPE using the 
       subpath starting at edge 1.  Once the type of crossing zero is set, the initialization of the polarity from the base point crossing means that the four cases a) through 
       d) above determine the type of crossing relative to that of the base point crossings.       
	*/
	
	vector<int> crossing_type(num_crossings,generic_code_data::type::UNKNOWN);
	crossing_type[0] = generic_code_data::type::TYPE2;

	vector<int> candidate_base_point(num_crossings);
	
	int base_point = 0;
	
	bool finished = false;
	
	do
	{
		int sub_path_start_edge = 2*base_point+1;
		int sub_path_end_edge = peer[2*base_point]; 
		
		vector<int> sub_path(num_edges);
		sub_path[sub_path_start_edge] = 1;
		
		for (int i=0; i< num_edges; i++)
		{
			int edge = (sub_path_start_edge+i)%num_edges;
			sub_path[edge] = 1;
			
			if (edge == sub_path_end_edge)
				break;
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_dowker_code:   sub_path: ";
	for (int i=0; i< num_edges; i++)
		debug << sub_path[i] << ' ';
	debug << endl;
}	
						
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_dowker_code: next base_point crossing = " << base_point << ", sub_path_start_edge = " << sub_path_start_edge << ", sub_path_end_edge = " << sub_path_end_edge << endl;		
		
		vector<int> polarity(num_edges);
		polarity[sub_path_start_edge] = 1;
		
		/* identify the polarity of the edges in the subpath and remaining edges */
		for (int i=0; i < num_edges-1; i++)
		{
			int edge = (sub_path_start_edge+i)%num_edges;
			int peer_edge = peer[edge];
			int succeeding_edge = (edge+1)%num_edges;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_dowker_code:   i = " << i << ", edge = " << edge << ", peer_edge = " << peer_edge << ", succeeding_edge = " << succeeding_edge << endl;		
			
			if (sub_path[peer_edge] == 1)
			{
				polarity[succeeding_edge] = polarity[edge]*-1;
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_dowker_code:     peer_edge " << peer_edge << " lies on the sub-path from " << sub_path_start_edge << " to " << sub_path_end_edge 
	      << ", polarity[" << succeeding_edge << "] = " << polarity[succeeding_edge] << endl;		
}	      
			}
			else
			{
				polarity[succeeding_edge] = polarity[edge];
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_dowker_code:     peer_edge " << peer_edge << " does not lie on the sub-path from " << sub_path_start_edge << " to " << sub_path_end_edge 
	      << ", polarity[" << succeeding_edge << "] = " << polarity[succeeding_edge] << endl;		
}	      
			}
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_dowker_code:   polarity: ";
	for (int i=0; i< num_edges; i++)
		debug << polarity[i] << ' ';
	debug << endl;
}	
		
		/* we set as many crossing types as we can based on the sub-path by considering the non-sub-path edges in turn. */				
		for (int i=0; i < num_edges; i++)
		{
			if (sub_path[i] == 1)
				continue;
				
			int preceding_edge = (i-1+num_edges)%num_edges;
			int peer_edge = peer[i];
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_dowker_code:   edge " << i << ", preceding_edge = " << preceding_edge << ", peer_edge = " << peer_edge << endl;

			if (sub_path[peer_edge] == 1)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_dowker_code:     peer_edge " << peer_edge << " lies on the sub-path from " << sub_path_start_edge << " to " << sub_path_end_edge 
	      << ", polarity[" << i << "] = " << polarity[i] << ", polarity[" << peer_edge << "] = " << polarity[peer_edge] << endl;
}	      
				int crossing = (i%2==0? i/2: peer_edge/2);
				
				if (crossing_type[crossing] == 0)
				{
					int local_type_polarity =  polarity[i]*polarity[peer_edge];
					if ( (local_type_polarity == 1  && i%2 == 0) || (local_type_polarity == -1  && i%2 == 1) )
						crossing_type[crossing] = crossing_type[base_point];
					else
						crossing_type[crossing] = (crossing_type[base_point] == generic_code_data::type::TYPE1? generic_code_data::type::TYPE2 : generic_code_data::type::TYPE1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_dowker_code:     local_type_polarity = " << local_type_polarity << " crossing_type[" << crossing << "] set to " 
	      << (crossing_type[crossing] == generic_code_data::type::TYPE1? "TYPE1" : "TYPE2") << endl;
}	      

					/* If the peer of the preceeding edge is the preceeding edge of the peer, then we learn nothing new by considering the sub-path from crossing 
					   since, other than this edge and its peer, the subpath from crossing comprises the remaining edges for the subpath of the crossing involving 
					   the preceeding edge and its peer.  That means we have already set all of the crossing types that we could on the sub-path from crossing.
					*/
					if (abs(peer[preceding_edge] - peer_edge) != 1 && abs(peer[preceding_edge] - peer_edge) != num_edges-1)
					{
						candidate_base_point[crossing] = 1;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_dowker_code:     set candidate_base_point[" << crossing << "] to 1" << endl;
					}
				}
				else
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_dowker_code:     crossing_type[ " << crossing << "] already set to " 
	      << (crossing_type[crossing] == generic_code_data::type::TYPE1? "TYPE1" : "TYPE2") << endl;
}	      
				}
			}
			else
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_dowker_code:     peer_edge " << peer_edge << " does not lie on the sub-path from " << sub_path_start_edge << " to " << sub_path_end_edge 
	      << ", polarity[" << i << "] = " << polarity[i] << ", polarity[" << peer_edge << "] = " << polarity[peer_edge] << endl;
}	      
			}
		}
				
		candidate_base_point[base_point] = 0;
		
		finished = true;
		for (int i=0; i<num_crossings; i++)
		{
			if (candidate_base_point[i] !=0)
			{
				base_point = i;
				finished = false;
				break;
			}
		}		
	} while(!finished);

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "read_dowker_code: crossing_type: ";
	for (int i=0; i< num_crossings; i++)
		debug << crossing_type[i] << ' ';
	debug << endl;
}	
	
	ostringstream oss;
	oss << '[';
	for (int i=0; i< num_crossings; i++)
	{
		if (crossing_type[i] == generic_code_data::type::TYPE1)
			oss << '-';
		oss << dowker_code[i] -1 ;
		
		if (i < num_crossings-1)
			oss << ' ';				
	}

	oss << "]/";
	for (int i=0; i< num_crossings; i++)
	{
		if (crossing_label[i] == generic_code_data::label::POSITIVE)
			oss << '+';
		else
			oss << '-';
			
		if (i < num_crossings-1)
			oss << ' ';				
	}
	
	return oss.str();
}
