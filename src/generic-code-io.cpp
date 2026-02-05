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
		code_data.immersion_character = generic_code_data::character::KNOTOID;
	else if (input_string.find("L:") != string::npos)
		code_data.immersion_character = generic_code_data::character::LONG_KNOT;
	else
		code_data.immersion_character = generic_code_data::character::CLOSED;

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
	
	matrix<int> code_table(generic_code_data::table::CODE_TABLE_SIZE,num_crossings);
	
	for (int i=0; i< num_crossings; i++)
		code_table[generic_code_data::table::COMPONENT][i] = 0;

	vector<int> perm(num_crossings);
	vector<int> type(num_crossings);
	
	cptr = strchr(inbuf,'(');
    while(*cptr != '/')
	{
		int count = 0;

		/* We write the crossing numbers in the order they appear in the immersion code 
		   into perm together with their corresponding type.  Then we unravel 
		   the permuatation writing the PERM and generic_code_data::table::TYPE rows from INVERSE and generic_code_data::table::LABEL.
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
					code_data.immersion_character = generic_code_data::character::PURE_KNOTOID;
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
			code_table[generic_code_data::table::OPEER][perm[i]] = (2*perm[i+1]-1+num_edges)%num_edges;
			code_table[generic_code_data::table::EPEER][(code_table[generic_code_data::table::OPEER][perm[i]]-1)/2] = 2*perm[i];
			code_table[generic_code_data::table::TYPE][perm[i]] = type[i];
		}
		code_table[generic_code_data::table::OPEER][perm[count-1]] = (2*perm[0]-1+num_edges)%num_edges;
		code_table[generic_code_data::table::EPEER][(code_table[generic_code_data::table::OPEER][perm[count-1]]-1)/2] = 2*perm[count-1];
		code_table[generic_code_data::table::TYPE][perm[count-1]] = type[count-1];

		cptr++; // move over ')'
		while (*cptr == ' ') cptr++;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "read_immersion_code: interim code table produced from immersion code:" << endl;
	debug << "read_immersion_code:   type:    ";	
	for (int i=0; i< num_crossings; i++)
		debug << code_table[generic_code_data::table::TYPE][i] << ' ';
	debug << endl;
	debug << "read_immersion_code:   generic_code_data::table::OPEER: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_table[generic_code_data::table::OPEER][i] << ' ';
	debug << endl;
	debug << "read_immersion_code:   generic_code_data::table::EPEER: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_table[generic_code_data::table::EPEER][i] << ' ';
	debug << endl;
	debug << "read_immersion_code:   label: ";
	for (int i=0; i< num_crossings; i++)
		debug << code_table[generic_code_data::table::LABEL][i] << ' ';
	debug << endl;
}

	}

	/* now do the labels */
	int count = 0;
	while (*cptr != '\0')
	{
		if (*cptr == '+')
			code_table[generic_code_data::table::LABEL][count++] = generic_code_data::POSITIVE;
		else if (*cptr == '-')
			code_table[generic_code_data::table::LABEL][count++] = generic_code_data::NEGATIVE;
		else if (*cptr == '*')
			code_table[generic_code_data::table::LABEL][count++] = generic_code_data::VIRTUAL;
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
		
		code_table[generic_code_data::table::EVEN_TERMINATING][i] = 2*i;
		code_table[generic_code_data::table::ODD_ORIGINATING][i] = 2*i+1;
		code_table[generic_code_data::table::ODD_TERMINATING][i] = (2*code_table[PERM][i]-1 + num_edges) % num_edges;
		code_table[generic_code_data::table::EVEN_ORIGINATING][i] = (2*code_table[PERM][i]) % num_edges;
	}
*/
	for (int i=0; i< num_crossings; i++)
	{
		term_crossing[2*i] = i;
		orig_crossing[2*i+1] = i;
		term_crossing[code_table[generic_code_data::table::OPEER][i]] = i;
		
		/* we need to identify the edge following the naming edge's peer */
		int component = code_table[generic_code_data::table::COMPONENT][(code_table[generic_code_data::table::OPEER][i]-1)/2];
		int peer_successor = (code_table[generic_code_data::table::OPEER][i]+1 - first_edge_on_component[component])%
		                     num_component_edges[component] + first_edge_on_component[component];

		orig_crossing[peer_successor] = i;
		
		code_table[generic_code_data::table::EVEN_TERMINATING][i] = 2*i;
		code_table[generic_code_data::table::ODD_ORIGINATING][i] = 2*i+1;
		code_table[generic_code_data::table::ODD_TERMINATING][i] = code_table[generic_code_data::table::OPEER][i];
		code_table[generic_code_data::table::EVEN_ORIGINATING][i] = peer_successor;
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

/* This function only uses head, num_crossings and the generic_code_data::table::OPEER, generic_code_data::table::TYPE, and generic_code_data::table::LABEL rows 
   of the code_table in the generic_code_data structure.
*/
void write_immersion_code(ostream& s, generic_code_data& code_data)
{
	matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;

	if (code_data.immersion_character == generic_code_data::character::KNOTOID)
		s << "K:";
	else if (code_data.immersion_character == generic_code_data::character::LONG_KNOT)
		s << "L:";
	
	vector<int> flag(num_crossings); // used to track which crossings have been written
	for (int i=0; i<num_crossings; i++)
		flag[i] = 1;
		
	/* work out the permutation from the generic code data */
	vector<int> perm(num_crossings);
	for (int i=0; i< num_crossings; i++)
		perm[i] = ((code_table[generic_code_data::table::OPEER][i]+1)/2)%num_crossings;
	
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
			if (code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::TYPE1)
				s << '-';
			s << crossing;

			int next_crossing;
			do
			{
				next_crossing = perm[crossing];
				s << " ";
				if (code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1)
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
		if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE)
			s << "+ ";
		if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
			s << "- ";
		if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::VIRTUAL)
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
	switch (code_data.immersion_character)
	{
		case (generic_code_data::character::CLOSED): s << "CLOSED"; break;
		case (generic_code_data::character::LONG_KNOT): s << "LONG_KNOT"; break;
		case (generic_code_data::character::PURE_KNOTOID): s << "PURE_KNOTOID"; break;
		case (generic_code_data::character::KNOTOID): s << "KNOTOID"; break;		
		case (generic_code_data::character::MULTI_LINKOID): s << "MULTI_LINKOID"; break;		
		default: s << "unknown";
	};	
	s << endl;
	
	s << prefix << "multi_virtual = " << code_data.multi_virtual << endl;
	if (code_data.multi_virtual)
	{
		s << prefix << "virtual_index: ";
		for (int i=0; i< num_crossings; i++)
			s << code_data.virtual_index[i] << ' ';
		s << endl;
	}
	
	s << prefix << "head = " << code_data.head << endl;
	s << prefix << "num_open_components = " << code_data.num_open_components << endl;
	s << prefix << "head_zig_zag_count = " << code_data.head_zig_zag_count << endl;
	s << prefix << "num_crossings = " << num_crossings << endl;
	s << prefix << "num_components = " << code_data.num_components << endl;
	
	if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID)
	{
		s << prefix << "component_type = ";
		for (int i=0; i< code_data.num_components; i++)
		{
			s << i << ":(";
			switch (code_data.component_type[i].type)
			{
				case(component_character::CLOSED): s << "CLOSED, "; break;
				case(component_character::PURE_START_LEG): s << "PURE_START_LEG," << code_data.component_type[i].head_semi_arc; break;
				case(component_character::PURE_END_LEG): s << "PURE_END_LEG," << code_data.component_type[i].head_semi_arc; break;
				case(component_character::KNOT_TYPE_START_LEG): s << "KNOT_TYPE_START_LEG,"; break;
				case(component_character::KNOT_TYPE_END_LEG): s << "KNOT_TYPE_END_LEG,"; break;
			};
			s << ") ";
		}
		s << endl;
	}
	
	s << prefix << "type:      ";	
	matrix<int>& code_table = code_data.code_table;
	
	for (int i=0; i< num_crossings; i++)
		s << code_table[generic_code_data::table::TYPE][i] << ' ';
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
			s << code_table[generic_code_data::table::OPEER][i] << ' ';
		s << endl;
		s << prefix << "even peer:  ";
		for (int i=0; i< num_crossings; i++)
			s << code_table[generic_code_data::table::EPEER][i] << ' ';
		s << endl;
	}
	
	s << prefix << "even term: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[generic_code_data::table::EVEN_TERMINATING][i] << ' ';
	s << endl;
	s << prefix << "odd term:  ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[generic_code_data::table::ODD_TERMINATING][i] << ' ';
	s << endl;
	s << prefix << "even orig: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[generic_code_data::table::EVEN_ORIGINATING][i] << ' ';
	s << endl;
	s << prefix << "odd orig:  ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[generic_code_data::table::ODD_ORIGINATING][i] << ' ';
	s << endl;
	s << prefix << "label:     ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[generic_code_data::table::LABEL][i] << ' ';
	s << endl;	
	s << prefix << "component: ";
	for (int i=0; i< num_crossings; i++)
		s << code_table[generic_code_data::table::COMPONENT][i] << ' ';
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
   It uses the #define generic_code_data::table::TYPE, generic_code_data::table::OPEER, generic_code_data::table::EPEER, generic_code_data::table::COMPONENT, generic_code_data::table::LABEL, generic_code_data::table::EVEN_TERMINATING, generic_code_data::table::ODD_TERMINATING, 
   generic_code_data::table::EVEN_ORIGINATING, generic_code_data::table::ODD_ORIGINATING to index the rows of the matrix.
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
	{
		code_data.immersion_character = generic_code_data::character::KNOTOID;
		code_data.num_open_components = 1;
	}
	else if (input_string.find("L:") != string::npos)
		code_data.immersion_character = generic_code_data::character::LONG_KNOT;
	else if (input_string.find("$") != string::npos || input_string.find("%") != string::npos ||
	         count(input_string.begin(),input_string.end(),'^') > 1)
	{
		int num_components = count(input_string.begin(),input_string.end(),',')+1;
		
		code_data.immersion_character = generic_code_data::character::MULTI_LINKOID;
		code_data.component_type = vector<component_character>(num_components);

		code_data.num_open_components = count(input_string.begin(),input_string.end(),'^') + count(input_string.begin(),input_string.end(),'$')
		                                + count(input_string.begin(),input_string.end(),'%');
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_peer_code: multi-linkoid inialization determined " << num_components << " components " << code_data.num_open_components << " open components"  << endl;
	
	}
	else
		code_data.immersion_character = generic_code_data::character::CLOSED;

	
	char* inbuf = c_string(input_string);
	int num_crossings = 0;
	char* cptr = strchr(inbuf,'/');
	cptr++;

	while (*cptr != '\0')
	{
		if (*cptr == '+' || *cptr == '-' || *cptr == '*' || *cptr == '#' || *cptr == '@' )
			num_crossings++;

		cptr++;
	}	

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code: num_crossings determined from labelled peer code = " << num_crossings << endl;
	
	code_data.num_crossings = num_crossings;
	
	matrix<int> code_table(generic_code_data::table::CODE_TABLE_SIZE,num_crossings);

	int num_edges = 2*num_crossings;
	int component = 0;
//	bool component_start = true; // used for multi-linkoid knot type components 

	/* assume all crossings are generic_code_data::TYPE2 and set those 
	   indicated as generic_code_data::TYPE1 accordingly.
	*/
	for (int i=0; i<num_crossings; i++)
		code_table[generic_code_data::table::TYPE][i] = generic_code_data::TYPE2;
		
	cptr = strchr(inbuf,'[');
	cptr++;
	vector<int> odd_head_semi_arc(num_crossings,generic_code_data::character::CLOSED);  // odd_head_semi_arc[i] records the component type if edge label 2i+1 is a head semi_arc
	
    for (int i=0; i< num_crossings; i++)
	{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code:   crossing " << i << ": " << endl;

		/* We write the odd-numbered peer of edge 2i in code_table[generic_code_data::table::OPEER][i] and the 
		   even-numbered peer of edge 2i+1 in code_table[generic_code_data::table::EPEER][i] for i=0,...,num_crossings-1
		*/
		bool crossing_complete = false;
		int record_odd_head_semi_arc = generic_code_data::character::CLOSED;
		do
		{				
			if (*cptr == ']')
			{
				cout << "\nError! More labels than peers specified in peer code" << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_peer_code: Error! More labels than peers specified in peer code" << endl;
					
				exit(0);
			}
			if (*cptr == '-')
			{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code:     type I indicator " << "\'-\'" << endl;
	
				code_table[generic_code_data::table::TYPE][i] = generic_code_data::TYPE1;
				cptr++;
			}
			else if (isdigit(*cptr))
			{
				char* mark = cptr;
				while (isdigit(*cptr)) cptr++;

				get_number (code_table[generic_code_data::table::OPEER][i],mark);

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code:     read odd peer " << code_table[generic_code_data::table::OPEER][i] << endl;
				
				if (code_table[generic_code_data::table::OPEER][i]%2 == 0)
				{
					cout << "\nError! Read even number " << code_table[generic_code_data::table::OPEER][i] << " in peer code" << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_peer_code: Error!  Read even number " << code_table[generic_code_data::table::OPEER][i] << " in a peer code" << endl;
					
					exit(0);
				}
				
				code_table[generic_code_data::table::EPEER][(code_table[generic_code_data::table::OPEER][i]-1)/2] = 2*i;
				code_table[generic_code_data::table::COMPONENT][i] = component;

				if (record_odd_head_semi_arc != generic_code_data::character::CLOSED)					
				{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code:     recording odd_head_smi_arc at index " << (code_table[generic_code_data::table::OPEER][i]-1)/2 << endl;

					odd_head_semi_arc[(code_table[generic_code_data::table::OPEER][i]-1)/2] = record_odd_head_semi_arc; // this odd peer is a head semi-arc
					record_odd_head_semi_arc = generic_code_data::character::CLOSED;
				}
				
				if (*cptr == '^')
				{
					/* we've found the knotoid head, note the corresponding crossing number */
					if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID)
					{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code:     PURE_START_LEG even head indicator " << "\'^\'" << endl;

						code_data.component_type[component].type = component_character::PURE_START_LEG;
						code_data.component_type[component].head_semi_arc = 2*i; // the even peer terminating at this crossing
					}
					else
					{
						code_data.head = i;
						code_data.immersion_character = generic_code_data::character::PURE_KNOTOID;
						code_data.num_open_components = 1;
					}
					cptr++;
				}			
				else if (*cptr == '$' && code_data.immersion_character == generic_code_data::character::MULTI_LINKOID)
				{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code:     PURE_END_LEG even head indicator " << "\'$\'" << endl;

					code_data.component_type[component].type = component_character::PURE_END_LEG;
					code_data.component_type[component].head_semi_arc = 2*i; // the even peer terminating at this crossing
					cptr++;
				}			
				else if (*cptr == '%' && code_data.immersion_character == generic_code_data::character::MULTI_LINKOID) // % sign after peer
				{

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code:     KNOT_TYPE_END_LEG odd head indicator " << "\'%\'" << endl;

					code_data.component_type[component].type = component_character::KNOT_TYPE_END_LEG;
					cptr++;
				}			
				
				crossing_complete = true;
//				component_start = false;
				
			}
			else if (*cptr == '^')
			{
				/* we've found the head of a multi-linkoid component but we haven't read the odd peer yet */
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code:     PURE_START_LEG odd head indicator " << "\'^\'" << endl;
				record_odd_head_semi_arc = component_character::PURE_START_LEG;
				cptr++;
			}			
			else if (*cptr == '$' && code_data.immersion_character == generic_code_data::character::MULTI_LINKOID)
			{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code:     PURE_END_LEG odd head indicator " << "\'$\'" << endl;

				record_odd_head_semi_arc = component_character::PURE_END_LEG;
				cptr++;
			}			
			else if (*cptr == '%' && code_data.immersion_character == generic_code_data::character::MULTI_LINKOID) // % before peer
			{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code:     KNOT_TYPE_START_LEG even head indicator " << "\'%\'" << endl;
	
				code_data.component_type[component].type = component_character::KNOT_TYPE_START_LEG;
				cptr++;
			}			
			else if (*cptr == ',')
			{
				/* we've come to the end of a component and are about to 
				   read the first crossing of the new component
				*/
				component++;
				cptr++;
//				component_start = true;
			}
			else
				cptr++;		
		} while (!crossing_complete);
	}
	
	/* complete the head_semi_arc of any components with odd head semi-arcs */
	if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID)
	{
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "read_peer_code:  odd_head_semi_arc: ";
	for (int i=0; i< num_crossings; i++)
		debug << odd_head_semi_arc[i] << ' ';
	debug << endl;
}

		for (int i=0; i< num_crossings; i++)
		{	
			if (odd_head_semi_arc[i] != component_character::CLOSED)
			{
				int semi_arc = 2*i+1;
				int component = code_table[generic_code_data::table::COMPONENT][(semi_arc-1)/2];
				code_data.component_type[component].type = odd_head_semi_arc[i];
				code_data.component_type[component].head_semi_arc = semi_arc;
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code:    odd_head_semi_arc[" << i << "] = " << odd_head_semi_arc[i] << ", semi_arc = " << semi_arc << ", component = " << component << endl;
			}
		}
	}
	
	/* now do the labels */
	int count = 0;
	vector<int> virtual_index(num_crossings);
	while (*cptr != '\0')
	{
		if (isdigit(*cptr))
		{
			cout << "\nError! More peers than labels specified in peer code" << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_peer_code: Error! More peers than labels specified in peer code" << endl;
					
			exit(0);
		}

		
		if (*cptr == '+')
			code_table[generic_code_data::table::LABEL][count++] = generic_code_data::POSITIVE;
		else if (*cptr == '-')
			code_table[generic_code_data::table::LABEL][count++] = generic_code_data::NEGATIVE;
		else if (*cptr == '*')
		{
			code_table[generic_code_data::table::LABEL][count++] = generic_code_data::VIRTUAL;
			if (isdigit(*(cptr+1)))
			{
				char* mark = cptr+1;
				while (isdigit(*(cptr+1))) cptr++;

				get_number (virtual_index[count-1],mark); // count has already been incremented
				
				if (virtual_index[count-1] != 0 && virtual_index[count-1] != 1)
				{
					code_data.multi_virtual = true;
				}
			}
			else
			{
				virtual_index[count-1] = 1;
			}
		}
		else if (*cptr == '#')
			code_table[generic_code_data::table::LABEL][count++] = generic_code_data::FLAT;
		else if (*cptr == '@')
			code_table[generic_code_data::table::LABEL][count++] = generic_code_data::SINGULAR;
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
		if (code_table[generic_code_data::table::COMPONENT][i] != code_table[generic_code_data::table::COMPONENT][i-1])
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
		term_crossing[code_table[generic_code_data::table::OPEER][i]] = i;
		
		/* we need to identify the edge following the naming edge's peer */
		int component = code_table[generic_code_data::table::COMPONENT][(code_table[generic_code_data::table::OPEER][i]-1)/2];
		int peer_successor = (code_table[generic_code_data::table::OPEER][i]+1 - first_edge_on_component[component])%
		                     num_component_edges[component] + first_edge_on_component[component];

if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "read_peer_code: crossing " << i << " peer_successor = " << peer_successor << endl;

		orig_crossing[peer_successor] = i;
		
		code_table[generic_code_data::table::EVEN_TERMINATING][i] = 2*i;
		code_table[generic_code_data::table::ODD_ORIGINATING][i] = 2*i+1;
		code_table[generic_code_data::table::ODD_TERMINATING][i] = code_table[generic_code_data::table::OPEER][i];
		code_table[generic_code_data::table::EVEN_ORIGINATING][i] = peer_successor;
	}
	
	code_data.code_table = code_table;
	code_data.num_component_edges = num_component_edges;
	code_data.first_edge_on_component = first_edge_on_component;
	code_data.term_crossing = term_crossing;
	code_data.orig_crossing = orig_crossing;
	code_data.virtual_index = virtual_index;
	
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
{
	debug << "read_peer_code: code data produced from peer code:" << endl;
	print_code_data(debug,code_data,"read_peer_code: ");	
}

	delete[] inbuf;
}

/* This function uses:

	immersion_character, 
	num_crossings, 
	num_components, 
    num_component_edges
	generic_code_data::table::ODD_TERMINATING, 
	generic_code_data::table::EVEN_TERMINATING
    generic_code_data::table::TYPE, 
    generic_code_data::table::LABEL
    generic_code_data::table::COMPONENT 
   
   If the code to be written is a PURE_KNOTOID, then the component_type vector is required.
   If the code to be written is a MULTI_LINKOID, then the head is required.
   If the code to be written has multi_virtual = true, then virtual_index is required
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
	
	if (code_data.immersion_character == generic_code_data::character::KNOTOID)
		s << "K:";
	else if (code_data.immersion_character == generic_code_data::character::LONG_KNOT)
		s << "L:";
	
	s << '[';
	
	int crossing = 0;
	
	if (code_data.num_crossings > 0)
	{
		for (int i=0; i< code_data.num_components; i++)
		{

			if (crossing > 0 && code_table[generic_code_data::table::COMPONENT][crossing] != code_table[generic_code_data::table::COMPONENT][crossing-1])
				s << ", ";
			
			if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID && code_data.component_type[i].type == component_character::KNOT_TYPE_START_LEG)
				s << '%';
			
			for (int j = 0; j< code_data.num_component_edges[i]/2; j++)
			{
				int odd_component = code_table[generic_code_data::table::COMPONENT][(code_table[generic_code_data::table::ODD_TERMINATING][crossing] -1)/2];
				int even_component = code_table[generic_code_data::table::COMPONENT][code_table[generic_code_data::table::EVEN_TERMINATING][crossing]/2];
				
				if (code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::TYPE1)
					s << '-';
					
				if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID && code_data.component_type[odd_component].type == component_character::PURE_START_LEG && 
				     code_data.component_type[odd_component].head_semi_arc == code_table[generic_code_data::table::ODD_TERMINATING][crossing])				   
					s << '^';
				else if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID && code_data.component_type[odd_component].type == component_character::PURE_END_LEG 
				         && code_data.component_type[odd_component].head_semi_arc == code_table[generic_code_data::table::ODD_TERMINATING][crossing])
					s << '$';
					
//				s << code_table[generic_code_data::table::OPEER][crossing];
				s << code_table[generic_code_data::table::ODD_TERMINATING][crossing];
				
				if ((code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && code_data.head == crossing) || 
				    (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID && code_data.component_type[even_component].type == component_character::PURE_START_LEG && 
				     code_data.component_type[even_component].head_semi_arc == code_table[generic_code_data::table::EVEN_TERMINATING][crossing])
				   )
					s << '^';
				else if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID && code_data.component_type[even_component].type == component_character::PURE_END_LEG 
				         && code_data.component_type[even_component].head_semi_arc == code_table[generic_code_data::table::EVEN_TERMINATING][crossing])
					s << '$';
				
				if ( crossing < num_crossings-1 && code_table[generic_code_data::table::COMPONENT][crossing] == code_table[generic_code_data::table::COMPONENT][crossing+1])
					s << ' ';

				crossing++;
			}
	
			if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID && code_data.component_type[i].type == component_character::KNOT_TYPE_END_LEG)
				s << '%';
		}
	}
	
/*	
	if (code_table[generic_code_data::table::TYPE][0] == generic_code_data::TYPE1)
		s << '-';
	if (code_data.num_crossings >0)
		s << code_table[generic_code_data::table::OPEER][0];
		
	if (code_data.head == 0)
		s << "^ ";
	else if (num_crossings > 1 && code_table[generic_code_data::table::COMPONENT][0] == code_table[generic_code_data::table::COMPONENT][1])
		s << ' ';
		
	for (int i=1; i<num_crossings; i++)
	{
		if (code_table[generic_code_data::table::COMPONENT][i] != code_table[generic_code_data::table::COMPONENT][i-1])
			s << ", ";
			
		if (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1)
			s << '-';
		s << code_table[generic_code_data::table::OPEER][i];
		if (code_data.head == i)
			s << '^';
			
		if ( i < num_crossings-1 && code_table[generic_code_data::table::COMPONENT][i] == code_table[generic_code_data::table::COMPONENT][i+1])
			s << ' ';
	}
*/
	
	s << "]";
	if (labelled)
	{
		s << "/";
		for (int i=0; i< num_crossings; i++)
		{
			if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE)
				s << "+";
			else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
				s << "-";
			else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::VIRTUAL)
			{
				s << "*";
				
				if (code_data.multi_virtual && code_data.virtual_index[i] != 0 && code_data.virtual_index[i] != 1)
					s << code_data.virtual_index[i];
			}
			else if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::FLAT)
				s << "#";
			else // if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
				s << "@";
			
			if (i< num_crossings-1)
				s << " ";
		}	
	}
	
	if (zig_zags && code_data.zig_zag_count.numrows() != 0)
	{
		if (code_data.immersion_character == generic_code_data::character::KNOTOID || code_data.immersion_character == generic_code_data::character::LONG_KNOT)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "read_gauss_code: O detected in input string, converting Gauss code format." << endl;
		input_string = convert_gauss_code(input_string);
	}
	
	if (input_string.find('*') != string::npos)
	{
	    cout << "Error! Gauss codes may not describe crossings as virtual" << endl;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "read_gauss_code: error in input string, '*' character found, indicating virtual crossing." << endl;
	
	    exit(0);
	}
	
	bool knotoid_or_long_indicator = false;
	int num_components = (int) count(input_string.begin(),input_string.end(),',') + 1;
	
	if (input_string.find("K:") != string::npos)
	{
		code_data.immersion_character = generic_code_data::character::KNOTOID;
		knotoid_or_long_indicator = true;
	}
	else if (input_string.find("K(") != string::npos)
	{
		code_data.immersion_character = generic_code_data::character::MULTI_LINKOID;
		get_number(code_data.num_open_components,input_string,input_string.find('(')+1);
		knotoid_or_long_indicator = true;		
		code_data.component_type = vector<component_character>(num_components);
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "read_gauss_code: num_open_components determined from gauss code " << input_string << " is " << code_data.num_open_components << endl;
	}
	else if (input_string.find("L:") != string::npos)
	{
		code_data.immersion_character = generic_code_data::character::LONG_KNOT;
		knotoid_or_long_indicator = true;
	}
	else
		code_data.immersion_character = generic_code_data::character::CLOSED;
	
	code_data.type = generic_code_data::gauss_code;
	code_data.head = -1;
//	code_data.immersion_character = generic_code_data::character::CLOSED;
	
//	int num_components = (int) count(input_string.begin(),input_string.end(),',') + 1;

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

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "read_gauss_code: num_crossings determined from gauss code " << input_string << " is " << num_crossings << endl;
	
	code_data.num_crossings = num_crossings;
	
	matrix<int> code_table(generic_code_data::table::CODE_TABLE_SIZE,num_crossings);
	
	int num_edges = 2*num_crossings;
	int component = 0;

	/* 	write the generic_code_data::table::OPEER row of the code table as -1 to act as flags	*/
	for (int i=0; i<num_crossings; i++)
		code_table[generic_code_data::table::OPEER][i] = -1;
		
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
		/* We write the edge on which we arrive at crossing i for the first time in code_table[generic_code_data::table::OPEER][i] and
		   the edge on which we arrive at the crossing for the second time in code_table[generic_code_data::table::EPEER][i].
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
	
				if (code_table[generic_code_data::table::OPEER][crossing] == -1)
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
					
					code_table[generic_code_data::table::OPEER][crossing] = edge;						
					code_table[generic_code_data::table::COMPONENT][crossing] = component;
					
					/* If we are at a classical crossing we can assign the crossing label as we now know whether 
					   the second visit (our "even" edge) will arrive on the over-arc or under-arc.  If the crossing 
					   is flat we can use the arc_type to set the generic_code_data::table::TYPE of the crossing as well as the label.
					*/
					if (arc_type == generic_code_data::POSITIVE)
						code_table[generic_code_data::table::LABEL][crossing] = generic_code_data::NEGATIVE;
					else if (arc_type == generic_code_data::NEGATIVE)
						code_table[generic_code_data::table::LABEL][crossing] = generic_code_data::POSITIVE;
					else
					{
						code_table[generic_code_data::table::LABEL][crossing] = generic_code_data::FLAT;
						
						if (arc_type == generic_code_data::LEFT)
							code_table[generic_code_data::table::TYPE][crossing] = generic_code_data::TYPE2;
						else
							code_table[generic_code_data::table::TYPE][crossing] = generic_code_data::TYPE1;
							
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "read_gauss_code:   set type = " << (code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::TYPE1? "TYPE1":"TYPE2") << endl;
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
					code_table[generic_code_data::table::EPEER][crossing] = edge;
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

	/* now we can identify the generic_code_data::table::TYPE of each clasical crossing from the sign of the 
	   crossing in the Gauss code and the generic_code_data::table::LABEL assignment.  Flat crossing have
	   had their generic_code_data::table::TYPE set already.
	*/
	int count = 0;
	while (*cptr != '\0')
	{
		if (*cptr == '+')
		{
			if (code_table[generic_code_data::table::LABEL][count] == generic_code_data::POSITIVE)
				code_table[generic_code_data::table::TYPE][count] = generic_code_data::TYPE2;
			else
				code_table[generic_code_data::table::TYPE][count] = generic_code_data::TYPE1;
				
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "read_gauss_code:   set type for crossing " << count << ": positive crossing, second visit " << (code_table[generic_code_data::table::LABEL][count] == generic_code_data::POSITIVE? "over": "under")
	      << " type is " << (code_table[generic_code_data::table::TYPE][count] == generic_code_data::TYPE1? "TYPE1":"TYPE2") << endl;
}
			count++;
		}
		else if (*cptr == '-')
		{
			if (code_table[generic_code_data::table::LABEL][count] == generic_code_data::POSITIVE)
				code_table[generic_code_data::table::TYPE][count] = generic_code_data::TYPE1;
			else
				code_table[generic_code_data::table::TYPE][count] = generic_code_data::TYPE2;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "read_gauss_code:   set type for crossing " << count << ": negative crossing, second visit " << (code_table[generic_code_data::table::LABEL][count] == generic_code_data::POSITIVE? "over": "under")
	      << " type is " << (code_table[generic_code_data::table::TYPE][count] == generic_code_data::TYPE1? "TYPE1":"TYPE2") << endl;
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
				code_table[generic_code_data::table::ODD_TERMINATING][crossing] = first_edge_on_component[i] + j;
				code_table[generic_code_data::table::EVEN_ORIGINATING][crossing] = first_edge_on_component[i] + (j+1)%num_component_edges[i];
				first_visit_to_crossing[crossing] = false;
			}
			else
			{				
				code_table[generic_code_data::table::EVEN_TERMINATING][crossing] = first_edge_on_component[i] + j;
				code_table[generic_code_data::table::ODD_ORIGINATING][crossing] = first_edge_on_component[i] + (j+1)%num_component_edges[i];
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
	
if (debug_control::DEBUG >= debug_control::DETAIL)
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
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "write_gauss_code: OU_FORMAT = " << OU_FORMAT << endl;
	print_code_data(debug,code_data,"write_gauss_code: ");
}

	int debug_save = debug_control::DEBUG;
//	debug_control::DEBUG = debug_control::OFF;

	matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;
	int num_components = code_data.first_edge_on_component.size();

	if (code_data.type == generic_code_data::gauss_code)
	{
		/* we are writing out a "native" Gauss code, so the code_data follows Gauss code conventions */
		
		/* we have to follow the edges around the Gauss code noting each crossing
		   we encounter as we do so.  Since Gauss codes record the comonent of the first
		   visit in the generic_code_data::table::COMPONENT row of the code table, we cannot use this row to output
		   the ',' between components accurately, since there may be two or more additional 
		   components in the diagram that are not recorded in that row.  Instead, we have to use
		   the edges recorded in first_edge_on_component to determine when we reach the end of a
		   component.
		   
		   Note that in the case of knotoids and multi-linkoids, the open component(s) will have been recorded
		   first and so will be written first, as required.
        */	  

		if (code_data.immersion_character == generic_code_data::character::KNOTOID)
			s << "K:";
		else if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID)
			s << "K(" << code_data.num_open_components << "):";
		else if (code_data.immersion_character == generic_code_data::character::LONG_KNOT)
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
				if (code_table[generic_code_data::table::OPEER][j] == i)
				{
					if (!OU_FORMAT && add_space)
						s << ' ';
	
					if (code_table[generic_code_data::table::LABEL][j] == generic_code_data::FLAT)
					{
						if (code_table[generic_code_data::table::TYPE][j] == generic_code_data::TYPE1)
							s << "R";
						else
							s << "L";						
					}
					else if (code_table[generic_code_data::table::LABEL][j] == generic_code_data::POSITIVE)
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
						if ( (code_table[generic_code_data::table::LABEL][j] == generic_code_data::POSITIVE && code_table[generic_code_data::table::TYPE][j] == generic_code_data::TYPE2) ||
						     (code_table[generic_code_data::table::LABEL][j] == generic_code_data::NEGATIVE && code_table[generic_code_data::table::TYPE][j] == generic_code_data::TYPE1)
						   )
							s << "+";
						else if ( (code_table[generic_code_data::table::LABEL][j] == generic_code_data::POSITIVE && code_table[generic_code_data::table::TYPE][j] == generic_code_data::TYPE1) ||
						          (code_table[generic_code_data::table::LABEL][j] == generic_code_data::NEGATIVE && code_table[generic_code_data::table::TYPE][j] == generic_code_data::TYPE2)
								)
							s << "-";
					}
					break;
				}
				else if (code_table[generic_code_data::table::EPEER][j] == i)
				{					
					if (!OU_FORMAT && add_space)
						s << ' ';
	
					if (code_table[generic_code_data::table::LABEL][j] == generic_code_data::FLAT)
					{
						if (code_table[generic_code_data::table::TYPE][j] == generic_code_data::TYPE1)
							s << "L";						
						else
							s << "R";
					}
					else if (code_table[generic_code_data::table::LABEL][j] == generic_code_data::NEGATIVE)
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
						if ( (code_table[generic_code_data::table::LABEL][j] == generic_code_data::POSITIVE && code_table[generic_code_data::table::TYPE][j] == generic_code_data::TYPE2) ||
						     (code_table[generic_code_data::table::LABEL][j] == generic_code_data::NEGATIVE && code_table[generic_code_data::table::TYPE][j] == generic_code_data::TYPE1)
						   )
							s << "+";
						else if ( (code_table[generic_code_data::table::LABEL][j] == generic_code_data::POSITIVE && code_table[generic_code_data::table::TYPE][j] == generic_code_data::TYPE1) ||
						          (code_table[generic_code_data::table::LABEL][j] == generic_code_data::NEGATIVE && code_table[generic_code_data::table::TYPE][j] == generic_code_data::TYPE2)
								)
							s << "-";
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
				if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::FLAT)
					s << "#";
				else if ( (code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE && code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE2) ||
				     (code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE && code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1)
				   )
					s << "+";
				else if ( (code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE && code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1) ||
				          (code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE && code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE2)
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
		vector<bool> component_traced(num_components); // initialized to zero

		ostringstream oss;
		
		/* we need to re-number the crossings if we have virtual crossings in the immersion, 
		   since these are ignored by the Gauss code.
		*/
		int num_classical_crossings = num_crossings;
		bool pure_knotoid_code_data = false;
		vector<int>& shortcut_crossing = code_data.shortcut_crossing;
		
		for (int i=0; i< num_crossings; i++)
		{
			if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::VIRTUAL)
				num_classical_crossings--;
		}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_gauss_code: num_classical_crossings = " << num_classical_crossings << endl;

		if ((code_data.immersion_character == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && shortcut_crossing.size()) ||
		    (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID && shortcut_crossing.size()))
		{
			pure_knotoid_code_data = true;
			if (code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
			{
				for (unsigned int i=0; i< shortcut_crossing.size(); i++)
				{
					if (shortcut_crossing[i] != 0)
						num_classical_crossings--;
				}
			}
if (debug_control::DEBUG >= debug_control::DETAIL)
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
		
		int component_count=1;
		int start=0;
		int edge=0;
		bool complete = false;

		component_traced[0] = true;

		if (code_data.immersion_character == generic_code_data::character::KNOTOID || code_data.immersion_character == generic_code_data::character::PURE_KNOTOID)
			oss << "K:";
		else if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID)
			s << "K(" << code_data.num_open_components << "):";
		else if (code_data.immersion_character == generic_code_data::character::LONG_KNOT)
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
				else if (code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::FLAT)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", is flat" << endl;

					if (!OU_FORMAT && add_space)
						oss << ' ';					
	
					if ((edge%2 != 0 && code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1) ||
					    (edge%2 == 0 && code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE2))
						oss << "R";
					else
						oss << "L";						
				}
				else if ((edge%2 != 0 && code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::POSITIVE) ||
				    (edge%2 == 0 && code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::NEGATIVE))
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
				
				if (code_table[generic_code_data::table::LABEL][next_crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[next_crossing]))
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
						if ( (code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::POSITIVE && code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE2) ||
						     (code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::NEGATIVE && code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1)
						   )
							oss << "+";
						else if ( (code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::POSITIVE && code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE1) ||
						          (code_table[generic_code_data::table::LABEL][next_crossing] == generic_code_data::NEGATIVE && code_table[generic_code_data::table::TYPE][next_crossing] == generic_code_data::TYPE2)
								)
							oss << "-";
					}



					if (edge%2)
						edge = code_table[generic_code_data::table::EVEN_ORIGINATING][next_crossing];
					else
						edge = code_table[generic_code_data::table::ODD_ORIGINATING][next_crossing];
	
					add_space = true;
					

				}
				else
				{
					/* just move on around the component */				
					if (edge%2)
						edge = code_table[generic_code_data::table::EVEN_ORIGINATING][next_crossing];
					else
						edge = code_table[generic_code_data::table::ODD_ORIGINATING][next_crossing];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_gauss_code:   doing nothing" << endl;
				}				
			} while (edge != start);

			/* look for the start of another component */
			complete = true;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_gauss_code: look for start of next component" << endl;
	
			for (int i=0; i< num_components; i++)
			{
				if (component_traced[i] == false)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_gauss_code:   component " << i << " not traced";

					if (code_data.immersion_character == generic_code_data::character::MULTI_LINKOID && component_count < code_data.num_open_components)
					{
						if (code_data.component_type[i].type == component_character::CLOSED)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", CLOSED component but component_count = " << component_count << " so skipping" << endl;
							continue;  // conside all open components first
						}
					}
					
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
			
			if (!complete)
				oss << ',';
				
		} while (!complete);


		if (!OU_FORMAT)
		{		
			oss << '/';
			
			for (int i=0; i< num_classical_crossings; i++)
			{
				int crossing = classical_crossing[i];
				
				if (code_table[generic_code_data::table::LABEL][crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[crossing])) // it should never be!
				{
					if (code_table[generic_code_data::table::LABEL][crossing] == generic_code_data::FLAT)
						oss << '#';
					else if ((code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::TYPE1 && code_table[generic_code_data::table::LABEL][crossing] == generic_code_data::NEGATIVE) ||
					    (code_table[generic_code_data::table::TYPE][crossing] == generic_code_data::TYPE2 && code_table[generic_code_data::table::LABEL][crossing] == generic_code_data::POSITIVE))
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

/* read_planar_diagram converts PD data in input_string of the form X[3,1,4,2] X[4,2,5,3] X[7,6,8,5] X[6,1,7,8] to a Gauss code and reads the resultant 
   Gauss code into code_data.  If the PD data includes a zero as an edge label we increment all elements to produce a code that is numbered from 1.
*/
void read_planar_diagram (generic_code_data& code_data, string input_string)
{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "read_planar_diagram: provided with input string: " << input_string << endl;
	
	int num_crossings = count(input_string.begin(),input_string.end(),'[');
	
	matrix<int> PD_data(num_crossings,4);

	istringstream iss(input_string);
	char c = 0;
	bool zero_label_detected = false;

	for (int i=0; i< num_crossings; i++)
	{
		while (c != '[')
			iss >> c;
		
		for (int j=0; j< 4; j++)
		{
			iss >> PD_data[i][j];
			iss >> c; // the comma
			
			if (PD_data[i][j] == 0)
				zero_label_detected = true;
		}
	}
	
	if (zero_label_detected)
	{
		for (int i=0; i< num_crossings; i++)
		{	
			for (int j=0; j< 4; j++)
				PD_data[i][j]++;
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
	else if (input_string.find("K(") == string::npos && input_string.find('(') != string::npos)
		read_immersion_code(code_data, input_string);
	else if (input_string.find("DT:") != string::npos)
	{
		string peer_code = read_dowker_code(input_string.substr(input_string.find(':')+1));
		read_peer_code(code_data, peer_code);
	}
	else
		read_gauss_code(code_data, input_string);
}


/* read_dowker_code creates a labelled peer code string corresponding to the Dowker-Thistlethwaite code of a prime knot, based on the conventions used by knotscape, where 
   the first Dowker-Thistlethwaite label (label 1) is assigned to an over-arc. That means that for alternating knots all the terms of the Dowker-Thistlethwaite code are 
   positive, since the even labels are always allocated to under-arcs at a crossing.  Since Dowker-Thistlethwaite codes are invariant under reflection in a line in the 
   plane disjoint from the diagram, knotscape fixes the generic_code_data::table::TYPE of the crossing involving label 1 and evaluates the other crossings' generic_code_data::table::TYPE using an algorithm based on 
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

       The clockwise or anti-clockwise relationship between the ingress edges is used to set the crossing generic_code_data::table::TYPE relative to the generic_code_data::table::TYPE of the base point crossing.  Since the 
       Dowker-Thistlethwaite code is invariant under reflection in a line of the plane, we may start by setting the generic_code_data::table::TYPE of crossing zero arbitrarily, ase described in more 
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
       which of the mirror images descibed by the Dowker-Thistlethwaite code the resulting peer code will describe, and start setting the other crossings' generic_code_data::table::TYPE using the 
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
