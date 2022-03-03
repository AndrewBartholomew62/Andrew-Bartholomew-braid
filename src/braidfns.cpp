/************************************************************************
                  Support functions for braid tools

                  A. Bartholomew 25th November, 2001

string invstr(string& str)
void get_word (ifstream& input, string& buffer, string& title)
void print_prog_params (ostream& os, string level, string prefix)
void convert_rep (char*& word)
string long_knot_concatenation (string ic1, string ic2)
string set_long_knot_infinity_point (string input_string)
string parse_long_knot_input_string (string input_string)
bool XxX_permutation (matrix<int>& U, matrix<int>& D)
bool rows_are_permutations (matrix<int>&M)
bool Yang_Baxter_satisfied (matrix<int>& U, matrix<int>& D, bool rack_condition = false)
bool one_dominates(matrix<int>& Ua, matrix<int>& Da, matrix<int>& Ub, matrix<int>& Db)
bool two_dominates(matrix<int>& Ua, matrix<int>& Da, matrix<int>& Ub, matrix<int>& Db)
bool power_2(matrix<int>& U, matrix<int>& D)
int switch_order(matrix<int>& U, matrix<int>& D)
void display_fixed_point_switch(matrix<int>& M, ostream& os, bool number_from_zero)
bool valid_knotoid_input(generic_code_data& code_data, vector<int>& shortcut_crossing)
void flip_braid(string& braid)
void invert_braid(string& braid)
void plane_reflect_braid(string& braid)
void line_reflect_braid(string& braid)
void Kamada_double_covering(string& braid)
void commutative_automorphism_invariant(const Qpmatrix& phi, const Qpmatrix& psi, string input_string, string title)

**************************************************************************/
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iomanip>

using namespace std;

extern ifstream     input;
extern ofstream 	output;
extern ofstream     debug;

extern bool SIDEWAYS_SEARCH;

#include <util.h>
#include <quaternion-scalar.h>
#include <polynomial.h>
#include <matrix.h>
#include <ctype.h>
#include <braid.h>

void read_peer_code (generic_code_data& code_data, string input_string);
void write_peer_code(ostream& s, const generic_code_data& code_data, bool zig_zags=false, bool labelled=true);
void read_immersion_code (generic_code_data& code_data, string input_string);
void write_immersion_code(ostream& s, generic_code_data& code_data);
void print_code_data(generic_code_data& code_data, ostream& s, string prefix="");
void write_code_data(ostream& s, generic_code_data& code_data);
void read_code_data (generic_code_data& code_data, string input_string);
void renumber_peer_code(generic_code_data& code_data, vector<int> shift);
bool valid_braid_input (string input_string, int& num_terms, int& num_strings);
string braid_to_generic_code (string braid, int num_terms, int code_type);
void parse_braid(string word, int num_terms, vector<int>& braid_num, vector<int>& type);
bool Yang_Baxter_satisfied (matrix<int>& U, matrix<int>& D, bool rack_condition = false);
bool rat_minor_determinant (Qpmatrix* Matrix_rep, int N, int Nr1, int Nc1, int Nc2, int* rperm, int* cperm, string title, polynomial<scalar,char,bigint>& hcf);

/* invstr creates the inverse of src and returns it to the call.
   The string src must be composed of si, -si, and ti components only.
   The inverse of ti is again ti.
*/
string invstr(string& str)
{
    char*	src = c_string(str);
    char*	loc_buf = new char[2*str.size() + 1]; // a slight overestimate, but always enough
    char*	mark = src;
    bool	not_finished = true;

    /* start at the end of the string and work back */
    char* cptr = strchr(src,0);
    if (strlen(src))
	    cptr--;  /* if src is the null string, cptr = src */

    char* lptr = loc_buf;
    do
    {
		if (isdigit(*cptr))
	    	cptr--;
		else
		{
		    if (*cptr == 's')
	    	{
				if (cptr !=src && *(cptr-1) == '-')
			    	mark = cptr-1;
				else
				{
		    		*lptr++ = '-';
		    		mark = cptr;
				}
				*lptr++ = 's';
	    	}
	    	else if (*cptr == 't')
	    	{
				mark = cptr;
				*lptr++ = 't';
		    }

		    cptr++;
		    while (isdigit(*cptr))
				*lptr++ = *cptr++;
	    	if (mark == src)
				not_finished = false;
	    	else
				cptr = mark-1;
		}
    } while (not_finished);

    *lptr = '\0';
    return string(loc_buf);
}

/* get_word retrieves the next braid word, Gauss code labelled peer code or 
   labelled immersion code from the input file, carrying out any manipulation 
   required, places the resultant word in buffer and any title string in title.  
   The function does not remove the '--' leader from the title and sets the title 
   to the empty string if none is found.  Any qualifiers are appended to the 
   input string, the function does not remove the braces {} from around qualifiers.
*/
void get_word (ifstream& input, string& buffer, string& title)
{
    string 	next_line;
	string qualifiers;
    string 	loc_buf;
	string	w1,w2,w3,w4,w5,w6,w7,w8;
	string* target;
	char*	lptr;
    bool   	word_found = false;
    bool	accept_braid_word = false;
	bool	accept_immersion_code = false;
	bool	accept_peer_code = false;
	bool	accept_gauss_code = false;
	bool	accepted_braid_input = false;
	bool	accepted_non_braid_input = false;
    bool	not_finished;

    /* start by setting the buffer and title to the empty string */
    buffer.clear();
    title.clear();

	if (braid_control::FIXED_POINT_INVARIANT || braid_control::RACK_POLYNOMIAL)
	{
		accept_braid_word = true;

		if (braid_control::first_time)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\n\nReading braid words from input.\n\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << "\n\n";
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Reading braid words from input.\n\n";
			}
			braid_control::first_time = false;
		}
	}
	else if (braid_control::SWITCH_POLYNOMIAL_INVARIANT)
	{
		accept_immersion_code = true;
		accept_peer_code = true;
		accept_braid_word = true;

		if (braid_control::first_time)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\n\nReading braid words, labelled immersion codes and labelled peer codes from input.\n\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << "\n\n";
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Reading braid words, labelled immersion codes and labelled peer codes from input.\n\n";
			}
			braid_control::first_time = false;
		}
	}
	else if (braid_control::VOGEL_ALGORITHM || braid_control::SATELLITE)
	{
		accept_peer_code = true;

		if (braid_control::first_time)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\n\nReading labelled peer codes from input.\n\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << "\n\n";
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Reading labelled peer codes from input.\n\n";
			}
			braid_control::first_time = false;
		}
	}
	else if (braid_control::GAUSS_CODE)
	{
		accept_immersion_code = true;
		accept_peer_code = true;
		accept_gauss_code = true;
		accept_braid_word = true;

		if (braid_control::first_time)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\n\nReading braid words, labelled immersion codes and labelled peer codes from input.\n\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << "\n\n";
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Reading braid words, labelled immersion codes and labelled peer codes from input.\n\n";
			}
			braid_control::first_time = false;
		}
	}
	else if (braid_control::PEER_CODE)
	{
		accept_immersion_code = true;
		accept_peer_code = true;
		accept_gauss_code = true;
		accept_braid_word = true;

		if (braid_control::first_time)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\n\nReading braid words, labelled immersion codes and labelled peer codes from input.\n\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << "\n\n";
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Reading braid words, labelled immersion codes and labelled peer codes from input.\n\n";
			}
			braid_control::first_time = false;
		}
	}
	else if (braid_control::DOWKER_CODE)
	{
		accept_braid_word = true;
		accept_immersion_code = true;

		if (braid_control::first_time)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\n\nReading braid words and labelled immersion codes from input.\n\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << "\n\n";
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Reading braid words and labelled immersion codes from input.\n\n";
			}
			braid_control::first_time = false;
		}
	}
	else if (braid_control::KAUFFMAN_BRACKET || braid_control::JONES_POLYNOMIAL || braid_control::ARROW_POLYNOMIAL)
	{
		accept_immersion_code = true;
		accept_peer_code = true;	
		accept_gauss_code = true;  // only Gauss codes containing classical crossings are supported currently, not flat or doodle crossings

		if (braid_control::first_time)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\n\nReading labelled immersion codes, labelled peer codes and Gauss codes from input.\n\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << "\n\n";
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Reading labelled immersion codes, labelled peer codes and Gauss codes from input.\n\n";
			}
			braid_control::first_time = false;
		}
	}
	else if (braid_control::KNOTOID_BRACKET || braid_control::PARITY_BRACKET || braid_control::PARITY_ARROW)
	{
		accept_immersion_code = true;
		accept_peer_code = true;

		if (braid_control::first_time)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\n\nReading labelled immersion codes and labelled peer codes from input.\n\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << "\n\n";
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Reading labelled immersion codes and labelled peer codes from input.\n\n";
			}
			braid_control::first_time = false;
		}
	}
	else if (braid_control::AFFINE_INDEX)
	{
		accept_immersion_code = true;
		accept_peer_code = true;
		accept_gauss_code = true; // only Gauss codes for classical diagrams are supported currently, not flat or doodle diagrams
		accept_braid_word = true;

		if (braid_control::first_time)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\n\nReading braid words, labelled immersion codes, labelled peer codes and Gauss codes from input.\n\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << "\n\n";
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Reading braid words, labelled immersion codes, labelled peer codes and Gauss codes from input.\n\n";
			}
			braid_control::first_time = false;
		}
	}
	else
	{
		accept_braid_word = true;

		if (braid_control::first_time)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\n\nReading braid words from input.\n\n";
			if (!braid_control::RAW_OUTPUT)
			{
				output << "\n\n";
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Reading braid words from input.\n\n";
			}
			braid_control::first_time = false;
		}
	}


if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
    debug << "braid::get_word: accepting ";
    if (accept_braid_word)
		debug << "braid words ";
    if (accept_immersion_code)
		debug << "immersion codes ";
    if (accept_peer_code)
		debug << "peer codes ";
    if (accept_gauss_code)
		debug << "Gauss codes ";
    debug << endl;
}
	
	while (!word_found && getline(input, next_line))
    {		

if (braid_control::DEBUG >= braid_control::DETAIL)
    debug << "braid::get_word: read line: " << next_line << endl;

		/* kill the <cr> at the end of the line, if one exists */
		string::size_type pos = next_line.find("\r");
		if (pos != string::npos)
		     next_line[pos] = ' ';
		
		/* remove any comments from the end of the line */
		pos = next_line.find(';');
		if (pos != string::npos)
			next_line = next_line.substr(0,pos);

 		
		if (next_line.find("exit") != string::npos)
			break;
			
		/* move any qualifiers from the end of the line into the qualifiers string 
		   by assigning the qualifier sting here any misplaced qualifiers in the input file 
		   are ignored.  Only when qualifiers are added to the end of a braid statement or the
		   last line of a peer code or immersion code or Gauss code will they be acted upon.
		   
		   Also, removing the qualifiers from the current line at this stage avoids any 
		   confusion between qualifiers and braid statement
		*/
		
		pos = next_line.find('{');
		if (pos != string::npos && next_line.find("--") != 0) // titles may contain braces
		{
			qualifiers = next_line.substr(pos);

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "braid::get_word: read qualifiers " << qualifiers << " from line " << next_line << endl;

			next_line = next_line.substr(0,pos);
			
		}

		char* line_buf = c_string(next_line);
				
	    /* if this line contains a switch, or programme options ignore it */
	    if ( strchr(line_buf,'[') && !strchr(line_buf,'/') && !strchr(line_buf,'\\'))
		{
if (braid_control::DEBUG >= braid_control::DETAIL)
    debug << "braid::get_word: line contains programme options, ignoring line" << endl;
			goto done_with_line;
		}
		else if (  (strchr(line_buf,'s') || strchr(line_buf,'S')) && 
			        strchr(line_buf,'=') && 
				    !(strchr(line_buf,'w') && isdigit(*(strchr(line_buf,'w')+1))) //not a braid statement
			    )
		{
if (braid_control::DEBUG >= braid_control::DETAIL)
    debug << "braid::get_word: line contains a switch, ignoring line" << endl;
			goto done_with_line;
		}

		
	    /* is there any whitespace at the beginning of the line? */
	    lptr = line_buf;
	    while (isspace(*lptr))
			lptr++;

	    if (strlen(lptr) && *lptr != '\n')
	    {
			
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "braid::get_word: next word line = " << lptr << endl;

			/* check for a title line */
			char* cptr = lptr;
			if (*cptr == '-' && *(cptr+1) == '-')
			{
		    	title = lptr;

				goto done_with_line;
			}

			/* the target of the line parsing is either one of the word buffers 
			   w1-w8, or is 'buffer', if the input is a braid statement, a peer code, an 
			   immersion code, or a Gauss code, .  An immersion code is 
			   indicated by the presence of a '(' character, a peer code by a '[' character.
			   A Gauss code is indicated by the presence of a '/' or '\' character but no '(' or '['.			
			*/
			if (!accepted_non_braid_input && !accepted_braid_input)
			{
				if (strchr(lptr,'(') && accept_immersion_code)
				{
					accepted_non_braid_input = true;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "braid::get_word: detected acceptable start of immersion code" << endl;
				}
				else if (strchr(line_buf,'[') && accept_peer_code)
				{
					accepted_non_braid_input = true;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "braid::get_word: detected acceptable start of peer code" << endl;
				}
				else if ((strchr(lptr,'/') || strchr(lptr,'O') || strchr(lptr,'\\')) && accept_gauss_code)
				{
					accepted_non_braid_input = true;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "braid::get_word: detected acceptable start of Gauss code" << endl;
				}
				else if ((strchr(line_buf,'s') || strchr(line_buf,'S') || 
				          strchr(line_buf,'t') || strchr(line_buf,'T')) && accept_braid_word)
				{
					accepted_braid_input = true;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "braid::get_word: detected acceptable start of braid word" << endl;
				}
				else
					goto done_with_line;
			}
			
			/* Here we have the kind of input we're looking for on the next_line */
			if (accepted_non_braid_input && (accept_gauss_code || accept_immersion_code || accept_peer_code ))
			{
			    /* Take out line escapes and build up either the
				   labelled immersion code or the Gauss code in buffer
			    */
			    cptr = strchr(lptr,'\\');
			    if (cptr)
					*cptr = ' ';

			    /* copy the line into buffer and decide if there's more to come */
				buffer += string(lptr);

			    /* This test means we cannot break lines after
				   the '/' character in the input file
				*/
				if (strchr(lptr,'/') || strchr(lptr,'O'))
				{
					word_found = true;					
					
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "braid::get_word: adding qualifiers " << qualifiers << " to buffer " << buffer << endl;;

					buffer += qualifiers;
				}
			    
				goto done_with_line;
			}
			
			/* Here we're looking for a braid word, first check whether this is 
			   an assignment statement or a braid statement, if it's a braid
			   statement, we're done once we've parsed this line. 
			
				Note the boolean accepted_braid_input is not actually needed in the
				current code but has been included for completeness and readability.
			*/
			cptr = strchr(lptr, '=');
			
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "braid::get_word: looking for a braid word...";

			if (cptr) // assignment statement
			{
			    cptr = strchr(lptr,'w');
		    	switch (*++cptr)
		    	{
					case '1': target = &w1; break;
					case '2': target = &w2; break;
					case '3': target = &w3; break;
					case '4': target = &w4; break;
					case '5': target = &w5; break;
					case '6': target = &w6; break;
					case '7': target = &w7; break;
					case '8': target = &w8; break;
					default: target = &w1;
		    	}				
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "line is an assignment statement" << endl;
			}
			else // braid statement
			{
			    target = &buffer;
		    	word_found = true;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "line is a braid statement" << endl;;
			}

			/* now parse the line into target, first move cptr to the
			   start of the line contents, i.e. past any = sign.
			*/

			cptr = strchr(lptr, '=');
			if (cptr)
			    cptr++;
			else
			    cptr = lptr;

			/* move through whitespace */
			while (isspace(*cptr))
			    cptr++;

			(*target).clear();
			not_finished = true;
			do
			{
			    if (isspace(*cptr) || *cptr == '\0')
					not_finished = false;
		    	else if (*cptr == 'w')
		    	{
					switch (*++cptr)
					{
			    		case '1': *target += w1; break;
			    		case '2': *target += w2; break;
			    		case '3': *target += w3; break;
			    		case '4': *target += w4; break;
			    		case '5': *target += w5; break;
			    		case '6': *target += w6; break;
			    		case '7': *target += w7; break;
			    		case '8': *target += w8; break;
			    		default: *target += w1;
					}
					cptr++;
		    	}
		    	else if (*cptr == '-' && *(cptr+1) == 'w')
		    	{
					cptr++; /* moves cptr to the w character */
					switch (*++cptr)
					{
					    case '1': *target += invstr(w1); break;
					    case '2': *target += invstr(w2); break;
				    	case '3': *target += invstr(w3); break;
				    	case '4': *target += invstr(w4); break;
			    		case '5': *target += invstr(w5); break;
				    	case '6': *target += invstr(w6); break;
				    	case '7': *target += invstr(w7); break;
				    	case '8': *target += invstr(w8); break;
				    	default: *target += invstr(w1);
					}
					cptr++;
		    	}
		    	else
		    	{
					/* copy the characters up to the next whitespace or 'w' into 
					   a local buffer and appent to target.
				    */
					char* copy_buf = new char[strlen(cptr)+1];
					char* sptr = copy_buf;
				
					while (!isspace(*cptr) && *cptr != '\0' && *cptr != 'w')
			    		*sptr++ = *cptr++;
					if (*cptr == 'w' && *(cptr-1) == '-')
					{
			    		/* move back one */
			    		cptr--;
			    		sptr--;
					}
					*sptr = '\0';
					*target += string(copy_buf);

					delete[] copy_buf;
		    	}
			} while (not_finished);
			
			if (word_found)
			{
				/* we have just parsed a braid statement into the buffer so we
				   append any qualifiers provided with that braid statement
				*/
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "braid::get_word: adding qualifiers " << qualifiers << " to buffer " << buffer << endl;;

				buffer += qualifiers;
			}
	    }
done_with_line:
		delete[] line_buf;
    } //end of while (getline(input, next_line) && !word_found)

    if (!word_found)
    	buffer = "exit";
		
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "braid::get_word: returning buffer " << buffer << endl;;
		
}

void print_prog_params (ostream& os, string level, string prefix)
{
	os << prefix << "braid_control::AFFINE_INDEX = " << braid_control::AFFINE_INDEX << endl;
	os << prefix << "braid_control::ALEXANDER = " << braid_control::ALEXANDER << endl;
	os << prefix << "braid_control::ARROW_POLYNOMIAL = " << braid_control::ARROW_POLYNOMIAL << endl;
	os << prefix << "braid_control::BURAU = " << braid_control::BURAU << endl;
	os << prefix << "braid_control::COMMUTATIVE_AUTOMORPHISM = " << braid_control::COMMUTATIVE_AUTOMORPHISM << endl;
	os << prefix << "braid_control::DOODLE_ALEXANDER = " << braid_control::DOODLE_ALEXANDER << endl;
	os << prefix << "braid_control::DOWKER_CODE = " << braid_control::DOWKER_CODE << endl;
	os << prefix << "braid_control::DYNNIKOV_TEST = " << braid_control::DYNNIKOV_TEST << endl;
	os << prefix << "braid_control::FIXED_POINT_INVARIANT = " << braid_control::FIXED_POINT_INVARIANT << endl;
	os << prefix << "braid_control::GAUSS_CODE = " << braid_control::GAUSS_CODE << endl;
	os << prefix << "braid_control::HOMFLY = " << braid_control::HOMFLY << endl;
	os << prefix << "braid_control::IMMERSION_CODE = " << braid_control::IMMERSION_CODE << endl;
	os << prefix << "braid_control::JONES_POLYNOMIAL = " << braid_control::JONES_POLYNOMIAL << endl;
	os << prefix << "braid_control::KAUFFMAN_BRACKET = " << braid_control::KAUFFMAN_BRACKET << endl;
	os << prefix << "braid_control::KNOTOID_BRACKET = " << braid_control::KNOTOID_BRACKET << endl;
	os << prefix << "braid_control::MATRIX = " << braid_control::MATRIX << endl;
	os << prefix << "braid_control::PARITY_ARROW = " << braid_control::PARITY_ARROW << endl;
	os << prefix << "braid_control::PARITY_BRACKET = " << braid_control::PARITY_BRACKET << endl;
	os << prefix << "braid_control::PEER_CODE = " << braid_control::PEER_CODE << endl;
	os << prefix << "braid_control::QUATERNION = " << braid_control::QUATERNION << endl;
	os << prefix << "braid_control::SAWOLLEK = " << braid_control::SAWOLLEK << endl;
	os << prefix << "braid_control::SATELLITE = " << braid_control::SATELLITE << endl;
	os << prefix << "braid_control::SWITCH_POLYNOMIAL_INVARIANT = " << braid_control::SWITCH_POLYNOMIAL_INVARIANT << endl;
	os << prefix << "braid_control::VOGEL_ALGORITHM = " << braid_control::VOGEL_ALGORITHM << endl;
	os << prefix << "braid_control::WELDED_BRAID = " << braid_control::WELDED_BRAID << endl;
	os << prefix << "braid_control::WEYL = " << braid_control::WEYL << endl;

	if (level == "detail" )
	{
		os << prefix << "braid_control::ALWAYS_CALCULATE_CODIMENSION_2_DELTA_1 = " << braid_control::ALWAYS_CALCULATE_CODIMENSION_2_DELTA_1 << endl;
		os << prefix << "braid_control::ALWAYS_CALCULATE_DELTA_1 = " << braid_control::ALWAYS_CALCULATE_DELTA_1 << endl;
		os << prefix << "braid_control::BIGELOW_KNOT_SEARCH = " << braid_control::BIGELOW_KNOT_SEARCH << endl;
		os << prefix << "braid_control::BYPASS_FUNDAMENTAL_EQUATION_CHECK = " << braid_control::BYPASS_FUNDAMENTAL_EQUATION_CHECK << endl;
		os << prefix << "braid_control::CALCULATE_DELTA_0_ONLY = " << braid_control::CALCULATE_DELTA_0_ONLY << endl;
		os << prefix << "braid_control::CALCULATE_MOD_P = " << braid_control::CALCULATE_MOD_P << endl;
		os << prefix << "braid_control::COMPLEX_STUDY_DELTA_1 = " << braid_control::COMPLEX_STUDY_DELTA_1 << endl;
		os << prefix << "braid_control::CUSTOM_WEYL = " << braid_control::CUSTOM_WEYL << endl;
		os << prefix << "braid_control::DELTA_1_UNIT_CHECK = " << braid_control::DELTA_1_UNIT_CHECK << endl;
		os << prefix << "braid_control::DISPLAY_DELTA_1_ONLY = " << braid_control::DISPLAY_DELTA_1_ONLY << endl;
		os << prefix << "braid_control::EQUALITY_TEST = " << braid_control::EQUALITY_TEST << endl;
		os << prefix << "braid_control::EXPANDED_BRACKET_POLYNOMIAL = " << braid_control::EXPANDED_BRACKET_POLYNOMIAL << endl;
		os << prefix << "braid_control::EXTRA_OUTPUT = " << braid_control::EXTRA_OUTPUT << endl;
		os << prefix << "braid_control::FLAT_VOGEL_MOVES = " << braid_control::FLAT_VOGEL_MOVES << endl;
		os << prefix << "braid_control::FLIP_BRAID = " << braid_control::FLIP_BRAID << endl;
		os << prefix << "braid_control::GCD_BACK_SUBSTITUTING = " << braid_control::GCD_BACK_SUBSTITUTING << endl;
		os << prefix << "braid_control::INVERT_BRAID = " << braid_control::INVERT_BRAID << endl;
		os << prefix << "braid_control::KAMADA_DOUBLE_COVERING = " << braid_control::KAMADA_DOUBLE_COVERING << endl;
		os << prefix << "braid_control::LINE_REFLECT_BRAID = " << braid_control::LINE_REFLECT_BRAID << endl;
		os << prefix << "braid_control::LONG_KNOT = " << braid_control::LONG_KNOT << endl;
		os << prefix << "braid_control::NORMALIZE_BRACKET = " << braid_control::NORMALIZE_BRACKET << endl;
		os << prefix << "braid_control::NORMALIZING_Q_POLYNOMIALS = " << braid_control::NORMALIZING_Q_POLYNOMIALS << endl;
		os << prefix << "braid_control::NUMERATOR_GCD = " << braid_control::NUMERATOR_GCD << endl;
		os << prefix << "braid_control::OUTPUT_AS_INPUT = " << braid_control::OUTPUT_AS_INPUT << endl;
		os << prefix << "braid_control::PLANE_REFLECT_BRAID = " << braid_control::PLANE_REFLECT_BRAID << endl;
		os << prefix << "braid_control::PRIME_WEYL = " << braid_control::PRIME_WEYL << endl;
		os << prefix << "braid_control::QUANTUM_WEYL = " << braid_control::QUANTUM_WEYL << endl;
		os << prefix << "braid_control::RACK_POLYNOMIAL = " << braid_control::RACK_POLYNOMIAL << endl;
		os << prefix << "braid_control::RAW_OUTPUT = " << braid_control::RAW_OUTPUT << endl;
		os << prefix << "braid_control::RELAXED_PARITY = " << braid_control::RELAXED_PARITY << endl;
		os << prefix << "braid_control::REMOVE_PEER_CODE_COMPONENT = " << braid_control::REMOVE_PEER_CODE_COMPONENT << endl;
		os << prefix << "braid_control::REMOVE_REIDEMEISTER_II_MOVES = " << braid_control::REMOVE_REIDEMEISTER_II_MOVES << endl;
		os << prefix << "braid_control::SILENT_OPERATION = " << braid_control::SILENT_OPERATION << endl;
		os << prefix << "braid_control::STATUS_INFORMATION = " << braid_control::STATUS_INFORMATION << endl;
		os << prefix << "braid_control::STUDY_RHO_MAPPING = " << braid_control::STUDY_RHO_MAPPING << endl;
		os << prefix << "braid_control::SWITCH_POWER = " << braid_control::SWITCH_POWER << endl;
		os << prefix << "braid_control::T_VARIABLE = " << braid_control::T_VARIABLE << endl;
		os << prefix << "braid_control::TMP_DIRECTORY = " << braid_control::TMP_DIRECTORY << endl;
		os << prefix << "braid_control::TRUNCATED_WEYL = " << braid_control::TRUNCATED_WEYL << endl;
		os << prefix << "braid_control::VERIFY_DELTA_0 = " << braid_control::VERIFY_DELTA_0 << endl;
		os << prefix << "braid_control::WAIT_SWITCH = " << braid_control::WAIT_SWITCH << endl;
	}
}


/* convert_rep takes a braid word in the standard representation using a,b,c,A,B,C etc.
   and converts it into si, -si format.
*/
void convert_rep (char*& word)
{
    string local;
    char* wptr=word;
    bool not_complete = true;
    bool inverse;

	
	do
	{
		if (*wptr == '\0')
			not_complete = false;	
	    else if (!isalpha(*wptr))
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "Error in input: non-alphabetical character " << *wptr << " \n";
			strcpy(word,"");
			return;
		}
		else
		{
			if (*wptr >= 'a')
				inverse = false;
			else
				inverse = false;

			/* write this term into local */
			if (inverse)
				local += "-";
			local += "s";
			ostringstream oss(local);
			oss << *wptr - (inverse?'A':'a')+1;
		}
		wptr++;
	    
	} while (not_complete);

	/* tidy up local, copy to word */
	delete[] word;
	word = c_string(local);
}

/* long_knot_concatenation yields a peer_code string representation of the concatination of the long knots
   determined by the input codes code1 and code2. The input codes can be either peer codes or immersion codes.

   To concatenate the long knots we create a combined generic_code_data structure and write it to a stringstream.
   Since writing only requires the head, num_components, OPEER, TYPE, LABEL and COMPONENT data, we calculate only
   these.  Note however, that since we are only dealing with knots, there is only one component to consider.
   
   We create a combined code table by juxtaposing code_table_1 and code_table_2 and then adjusting the OPEER row
   by incrementing those corresponding to code_data_2 by the number of edges corresponding to code_data_1.  
     
   The type and labels are unchanged by the concatenation and the COMPONENT row will correctly be populated with zeros.
   
*/
string long_knot_concatenation (string code1, string code2)
{

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "long_knot_concatenation: presented with strings :" << endl;
	debug << "long_knot_concatenation:   code1 = " << code1 << endl;
	debug << "long_knot_concatenation:   code2 = " << code2 << endl;
}
	generic_code_data code_data_1;
	generic_code_data code_data_2;

	read_code_data(code_data_1, code1);
	matrix<int>& code_table_1 = code_data_1.code_table;
	int n1 = code_data_1.num_crossings;

/*
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "long_knot_concatenation: code data produced from input code code1:" << endl;
	print (code_data_1,debug,"long_knot_concatenation: ");
}
*/
	read_code_data(code_data_2, code2);
	matrix<int>& code_table_2 = code_data_2.code_table;
	int n2 = code_data_2.num_crossings;
/*	
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "long_knot_concatenation: code table produced from input code code2:" << endl;
	print (code_data_2,debug,"long_knot_concatenation: "); 
}
*/
	matrix<int> code_table(CODE_TABLE_SIZE,n1 + n2);
	
	for (int i = 0; i<CODE_TABLE_SIZE; i++)
	for (int j = 0; j<n1; j++)
		code_table[i][j] = code_table_1[i][j];
		
	for (int i = 0; i<CODE_TABLE_SIZE; i++)
	for (int j = 0; j<n2; j++)
		code_table[i][n1+j] = code_table_2[i][j];
		
	int num_edges_1 = 2* code_data_1.num_crossings;
		
	for (int i = 0 ; i<n2; i++)
		code_table[OPEER][n1+i] += num_edges_1;
	
	/* write the new peer code to a stringstream */
	ostringstream oss;
	code_data_1.type = generic_code_data::peer_code;
//	code_data_1.head = -1;
	code_data_1.num_crossings += code_data_2.num_crossings;
//	code_data_1.num_components += code_data_2.num_components - 1 ;
	code_data_1.code_table = code_table;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
{
	debug << "long_knot_concatenation: concatenated code table:" << endl;
	print_code_data (code_data_1, debug,"long_knot_concatenation: ");
}

	write_peer_code (oss, code_data_1);
	return oss.str();
}


string set_long_knot_infinity_point (string input_string)
{
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "set_long_knot_infinity_point: presented with input_string: " << input_string << endl;

	generic_code_data code_data;
	read_code_data(code_data, input_string);

	if (code_data.num_components > 1)
	{
		return "link";
	}
	
	int shift = 0;
	
	char* inbuf = c_string(input_string);
	char* cptr = strchr(inbuf, 'L'); // We know there's an L in the string, otherwise this function wouldn't be called.
	cptr++; // step over L
	cptr++; // step over :
	if (*cptr == '+' || *cptr == '-')
	{
		/* a shift of infinity has been requested for the long knot */
		get_number(shift, cptr);
		
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "set_long_knot_infinity_point: long knot shift = " << shift << endl;

		/* move over number so it can be removed from the input string */
		cptr++;
		while (isdigit(*cptr))
			cptr++;
	}
			
	input_string.erase(0,cptr-inbuf);

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "set_long_knot_infinity_point: input string after L and shift removal " << input_string << endl;

	if (shift)
	{			
		/* The function renumber_peer_code moves the point at infinity forwards with respect
		   to the orientation, so if the requested shift is negtive, we have to adjust it accordingly.
		*/
		int num_crossings = code_data.num_crossings;

		// we do this repeatedly in case we've been given a huge negative shift.
		while (shift < 0)
			shift += 2*num_crossings;
		
		shift %= 2*num_crossings;
		
		vector<int> shift_vector(1); 
		shift_vector[0]=abs(shift);
		renumber_peer_code (code_data,shift_vector);

		/* write the new immersion code to a stringstream */
		ostringstream oss;
//		code_data.head = -1;
//		code_data.num_crossings = num_crossings;
//		code_data.code_table = code_table;

		write_code_data(oss,code_data);

		input_string = oss.str();
	}


if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "set_long_knot_infinity_point: returning immersion code " << input_string << endl;

	delete[] inbuf;

	return input_string;
}

string parse_long_knot_input_string (string input_string)
{

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "parse_long_knot_input_string: presented with input string " << input_string << endl;

	/* Check to see if this is a concatenation of multiple long knots */
	string::size_type tilde = input_string.find('~');

	string result = set_long_knot_infinity_point(input_string.substr(0,tilde));

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "parse_long_knot_input_string: initial component: " << result << endl;

	while( result != "link" && tilde != string::npos) // result != "link" catches the case where the first component is a link
	{
		string::size_type next_tilde = input_string.find('~',tilde+1);		
		string component = set_long_knot_infinity_point(input_string.substr(tilde+1,next_tilde-tilde-1));

		if (component =="link")
		{

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "parse_long_knot_input_string: next component: " << component << endl;

			result = component;
			break;
		}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "parse_long_knot_input_string: next component: " << component << endl;

		result = long_knot_concatenation (result, component);
		tilde = next_tilde;
	}
		
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "parse_long_knot_input_string: input string after parsing: " << result << endl;

	return result;

}

bool XxX_permutation (matrix<int>& U, matrix<int>& D)
{
	int n = U.numrows();
	
	for (int k=0; k< n; k++)
	{
		/* for each U[i][j]=k, D[j][i] must be a permutation */
		valarray<int> perm(0,n);
		for (int i=0; i< n; i++)
		for (int j=0; j< n; j++)
		{
			if (U[i][j] == k)
			{
				if (perm[D[j][i]])
					return false; // not a permutation
				else
					perm[D[j][i]] = 1;
			}
		}
	}				
	return true;
}

bool rows_are_permutations (matrix<int>&M)
{
	int n = M.numrows();
	
	for (int i=0; i<n; i++)
	{
		valarray<int> check_perm(0,n);
		for (int j=0; j< n; j++)
		{
			if (check_perm[M[i][j]])
				return false; // not a permutation
			else
				check_perm[M[i][j]] = 1;
		}				
	}
	return true;
}


/* If the boolean rack_condition is false (the default) the Yang_Baxter_satisfied function checks 
   the up down and rule of five conditions for a birack.  If rack_condition is true, it only
   tests the up condition for a rack.
*/

bool Yang_Baxter_satisfied (matrix<int>& U, matrix<int>& D,bool rack_condition)
{
	int n = U.numcols();

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "  S-Yang-Baxter rules: " << endl;

	vector<int> rof(3); //rule of five
	for (int i =0; i< n; i++)
	{
		rof[0]=i;
		for (int j =0; j< n; j++)
		{
			rof[1]=j;
			for (int k =0; k < n; k++)
			{
				rof[2]=k;
				
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "    rof: " << rof[0] << " " << rof[1] << " " << rof[2] << endl;

				int& a = rof[0];
				int& b = rof[1];
				int& c = rof[2];
					
				/* Up interchange */
				int& a_up_b = U[b][a];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      a_up_b = " << a_up_b <<endl;

				int& b_up_c = U[c][b];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      b_up_c = " << b_up_c <<endl;

				int& a_up_b__up_c = U[c][a_up_b];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      a_up_b__up_c = " << a_up_b__up_c <<endl;

				if (rack_condition)
				{
					/* a^bc = a^{cb^c} */
					int& a_up_c = U[c][a];						
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      a_up_c = " << a_up_c <<endl;

					if (a_up_b__up_c != U[b_up_c][a_up_c])
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      rack up interchange condition fails" << endl;
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "    fails" << endl;
	
						return false;
					}
					else
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      a_up_b__up_c == U[b_up_c][a_up_c]\n" << endl;
					}
				}
				else 
				{
					int& c_down_b = D[b][c];					
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      c_down_b = " << c_down_b <<endl;

					int& a_up__c_down_b = U[c_down_b][a];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      a_up__c_down_b = " << a_up__c_down_b <<endl;
					
					if (a_up_b__up_c != U[b_up_c][a_up__c_down_b])
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      up interchange condition fails" << endl;
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "    fails" << endl;
	
						return false;
					}
					else
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      a_up_b__up_c == U[b_up_c][a_up__c_down_b]\n" << endl;
					}
				}

				if (!rack_condition)
				{
					/* Down interchange */
					int& c_up_b = U[b][c];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      c_up_b = " << c_up_b <<endl;

					int& b_down_c = D[c][b];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      b_down_c = " << b_down_c <<endl;

					int& a_down_b = D[b][a];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      a_down_b = " << a_down_b <<endl;

					int& a_down_b__down_c = D[c][a_down_b];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      a_down_b__down_c = " << a_down_b__down_c <<endl;

					int& a_down__c_up_b = D[c_up_b][a];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      a_down__c_up_b = " << a_down__c_up_b <<endl;
					
					if (a_down_b__down_c != D[b_down_c][a_down__c_up_b])
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      down interchange condition fails" << endl;
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "    fails" << endl;
	
						return false;
					}
					else
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      a_down_b__down_c == D[b_down_c][a_down__c_up_b]\n" << endl;
					}
					
					/* Rule of five */
					int& b_up_a = U[a][b];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      b_up_a = " << b_up_a <<endl;
	
					int& a_up_c = U[c][a];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      a_up_c = " << a_up_c <<endl;
	
					int& c_down_a = D[a][c];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      c_down_a = " << c_down_a <<endl;
	
					int& b_up__c_down_a = U[c_down_a][b];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      b_up__c_down_a = " << b_up__c_down_a <<endl;
		
					int& c_down__b_up_a = D[b_up_a][c];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      c_down__b_up_a = " << c_down__b_up_a <<endl;
	
					int& a_down_b__up___c_down__b_up_a = U[c_down__b_up_a][a_down_b];
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      a_down_b__up___c_down__b_up_a = " << a_down_b__up___c_down__b_up_a <<endl;
	
					
					if (a_down_b__up___c_down__b_up_a != D[b_up__c_down_a][a_up_c])
					{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "      rule of five condition fails" << endl;
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "  fails" << endl;
	
						return false;
					}	
				}
			}
		}
	}
	
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "    OK" << endl;

	return true;
}


/* one_dominates tests whether the switch A, determined by up and down matrices Ua and Da
   (a candidate for S), 1-dominates the switch B, determined by up and down matrices Ub and Db
   (a candidate for T).
   
   The swith A 1-dominates B if 
   
   (A x I) (I x A) (B x 1) = (I x B) (A x I) (I x A)   (*)  
   
   Each case is checked by evaluating the matrix representing the composition on
   the lhs and rhs.

   Writing A(x,y) for the effect of the switch A, we have
   (A x I)(i1,i2,i3)=(A(i1,i2),i3) and (I x A)(i1,i2,i3)=(i1,A(i2,i3)), etc..
   
   Note that we are associating the labels i1, i2 and i3 to the braid strands in a top down 
   manner here, consistent with our convention of writing S(x,y) = (y^x,x_y).
   
   If we write S1 = Id x S, and S2 = S x Id, and similarly for T, this function checks the
   relationship
   
		S2S1T2 = T1S2S1
		
   which corresponds to the conditions imposed on S and T by the allowable weld move in a braid.   
   	   
   The effect of a switch S with up and down matrices U and D on the pair (x,y) is given by
   S(x,y) = (U[x][y],D[y][x]), since the (x,y) entry of U and D are y^x and y_x respectively
   and S(x,y) = (y^x,x_y).

   The function compares the lhs and rhs compositions in (*), the composition for the lhs gives 
   (i1,i2,i3) -> (l1,l2,l3), where (i1,i2,i3)->(j1,j2,j3) under A1, (j1,j2,j3)->(k1,k2,k3) 
   under A2 and (k1,k2,k3)->(l1,l2,l3) under B1.  Similarly for the rhs.
   
   We consider each triple (i1,i2,i3) and evaluate the corresponding (l1,l2,l3) for the lhs
   and (r1,r2,r3) for the rhs and then compare them, if they differ, then the test fails.
*/
bool one_dominates(matrix<int>& Ua, matrix<int>& Da, matrix<int>& Ub, matrix<int>& Db)
{
	int n = Ua.numcols();

	for (int i=0; i < n; i++)
	for (int j=0; j < n; j++)
	for (int k=0; k < n; k++)
	{
		int l1=i,l2=j,l3=k,temp;
		
		/* evaluate the lhs */
		
		/* (aptr x I)(x,y,z) = (Ua[x][y],Da[y][x],z) */
		temp = Ua[l1][l2];
		l2 = Da[l2][l1];
		l1 = temp;

		/* (I x aptr)(x,y,z) = (x, Ua[y][z],Da[z][y]) */
		temp = Ua[l2][l3];
		l3 = Da[l3][l2];
		l2 = temp;
		
		/* (bptr x I)(x,y,z) = (Ub[x][y],Db[y][x],z) */
		temp = Ub[l1][l2];
		l2 = Db[l2][l1];
		l1 = temp;

		int r1=i,r2=j,r3=k;

		/* evaluate the rhs */

		/* (I x bptr)(x,y,z) = (x, Ub[y][z],Db[z][y]) */
		temp = Ub[r2][r3];
		r3 = Db[r3][r2];
		r2 = temp;

		/* (aptr x I)(x,y,z) = (Ua[x][y],Da[y][x],z) */
		temp = Ua[r1][r2];
		r2 = Da[r2][r1];
		r1 = temp;
		
		/* (I x aptr)(x,y,z) = (x, Ua[y][z],Da[z][y]) */
		temp = Ua[r2][r3];
		r3 = Da[r3][r2];
		r2 = temp;
		
		if (r1 != l1 || r2 != l2 || r3 != l3)
			return false;
	}
	
	return true;
}

/* two_dominates tests whether the switch A, determined by up and down matrices Ua and Da
   (a candidate for S), 2-dominates the switch B, determined by up and down matrices Ub and Db
   (a candidate for T).
   
   The swith A 2-dominates B if 
   
   (I x A) (A x I) (I x B) = (B x I) (I x A) (A x I)   (*)  
    
   Each case is checked by evaluating the matrix representing the composition on
   the lhs and rhs.

   Writing A(x,y) for the effect of the switch A, we have
   (A x I)(i1,i2,i3)=(A(i1,i2),i3) and (I x A)(i1,i2,i3)=(i1,A(i2,i3)), etc..
   
   Note that we are associating the labels i1, i2 and i3 to the braid strands in a top down 
   manner here, consistent with our convention of writing S(x,y) = (y^x,x_y).
   
   If we write S1 = Id x S, and S2 = S x Id, and similarly for T, this function checks the
   relationship
   
		S1S2T1 = T2S1S2
		
   which corresponds to conditions imposed on S and T by the forbidden weld move in a braid.

   The effect of a switch S with up and down matrices U and D on the pair (x,y) is given by
   S(x,y) = (U[x][y],D[y][x]), since the (x,y) entry of U and D are y^x and y_x respectively
   and S(x,y) = (y^x,x_y).

   The function compares the lhs and rhs compositions in (*), the composition for the lhs gives 
   (i1,i2,i3) -> (l1,l2,l3), where (i1,i2,i3)->(j1,j2,j3) under A2, (j1,j2,j3)->(k1,k2,k3) 
   under A1 and (k1,k2,k3)->(l1,l2,l3) under B2.  Similarly for the rhs.
   
   We consider each triple (i1,i2,i3) and evaluate the corresponding (l1,l2,l3) for the lhs
   and (r1,r2,r3) for the rhs and then compare them, if they differ, then the test fails.
*/
bool two_dominates(matrix<int>& Ua, matrix<int>& Da, matrix<int>& Ub, matrix<int>& Db)
{
	int n = Ua.numcols();

	for (int i=0; i < n; i++)
	for (int j=0; j < n; j++)
	for (int k=0; k < n; k++)
	{
		int l1=i,l2=j,l3=k,temp;
		
		/* evaluate the lhs */
		
		/* (I x aptr)(x,y,z) = (x, Ua[y][z],Da[z][y]) */
		temp = Ua[l2][l3];
		l3 = Da[l3][l2];
		l2 = temp;

		/* (aptr x I)(x,y,z) = (Ua[x][y],Da[y][x],z) */
		temp = Ua[l1][l2];
		l2 = Da[l2][l1];
		l1 = temp;

		/* (I x bptr)(x,y,z) = (x, Ub[y][z],Db[z][y]) */
		temp = Ub[l2][l3];
		l3 = Db[l3][l2];
		l2 = temp;

		int r1=i,r2=j,r3=k;

		/* evaluate the rhs */
		
		/* (bptr x I)(x,y,z) = (Ub[x][y],Db[y][x],z) */
		temp = Ub[r1][r2];
		r2 = Db[r2][r1];
		r1 = temp;
		
		/* (I x aptr)(x,y,z) = (x, Ua[y][z],Da[z][y]) */
		temp = Ua[r2][r3];
		r3 = Da[r3][r2];
		r2 = temp;
		
		/* (aptr x I)(x,y,z) = (Ua[x][y],Da[y][x],z) */
		temp = Ua[r1][r2];
		r2 = Da[r2][r1];
		r1 = temp;

		if (r1 != l1 || r2 != l2 || r3 != l3)
			return false;
	}
	
	return true;
}

/* power_2 returns true if the switch indicated by the parameter has the property S^2=1
 
   The effect of a switch S with up and down matrices U and D on the pair (x,y) is given by
   S(x,y) = (U[x][y],D[y][x]), since the (x,y) entry of D and U are y_x and y^x respectively
   and S(x,y) = (y^x,x_y).
   
   Therefore power_2 looks at all combinations of (x,y) and checks that S(y^x,x_y) = (x,y)
*/ 
bool power_2(matrix<int>& U, matrix<int>& D)
{
	int n = U.numrows();
	
	for (int i=0; i < n; i++)
	for (int j=0; j < n; j++)
	{		
		if (U[ U[i][j] ][ D[j][i] ] != i)
			return false;

		if (D[ D[j][i] ][ U[i][j] ] != j)
			return false;
	}
	return true;
}

/* switch_order calculates the least k such that S^k(x,y) = (x,y).
   We shall calculate up and down matrices U^k and D^k for S^k so that
   S^k(x,y) = (Uk[x][y],Dk[y][x]) and will have
   S^k(x,y)=(x,y) if the (x,y) entry of U^k is x and the (y,x) entry of D^k is y,
   i.e the row i of U^k and D^k are both constant with entries equal to i.
   
   We know that the order of S divides n^2 since S is a permutation of XxX
*/
int switch_order(matrix<int>& U, matrix<int>& D)
{
	matrix<int> Uk(U),Dk(D);
	int n = U.numrows();		
	
	for (int k=1; k <= n*n; k++)
	{

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	bool saved_bool = matrix_control::SINGLE_LINE_OUTPUT;
	matrix_control::SINGLE_LINE_OUTPUT = true;
	debug << "braid::switch_order: k = " << k << " Uk = " << Uk << " Dk = " << Dk << endl;
	matrix_control::SINGLE_LINE_OUTPUT = saved_bool;
}

		/* test Uk and Dk for the identity */
		bool identity = true;
		
		for (int i=0; i<n && identity==true; i++)
		for (int j=0; j<n; j++)
		{
			if (Uk[i][j] != i)
			{
				identity = false;
				break;
			}
		}
		
		if (identity)
		{
			for (int i=0; i<n && identity==true; i++)
			for (int j=0; j<n; j++)
			{
				if (Dk[i][j] != i)
				{
					identity = false;
					break;
				}
			}
		}
		
		if (identity)
		{

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::switch_order: returning " << k << endl;

			return k;
		}
		else
		{
			/* calculate Uk and Dk for the next value of k.  We iterate using
			   S(x,y)=(y^x,x_y) = (U[x][y], D[y][x]) by writing
			   S^k(x,y) = (Uk[x][y],Dk[y][x]) and then we have
			   S^{k+1}(x,y) = S(S^k(x,y)) 
			                = S(Uk[x][y],Dk[y][x])
							= ((Dk[y][x])^(Uk[x][y]),(Uk[x][y])_(Dk[y][x]))
							= ( U[Uk[x][y]][Dk[y][x]], D[Dk[y][x]][Uk[x][y]] ) 
									   
			   Thus, 			   
			   
			   U{k+1}[x][y] = U[Uk[x][y]][Dk[y][x]]
			   D{k+1}[y][x] = D[Dk[y][x]][Uk[x][y]]
			   
			*/
			matrix<int> locU(n,n),locD(n,n);

			for (int i=0; i<n; i++)
			for (int j=0; j<n; j++)
			{
				locU[i][j] = U[Uk[i][j]][Dk[j][i]];
				locD[j][i] = D[Dk[j][i]][Uk[i][j]];
			}
			
			Uk = locU;
			Dk = locD;
		}
	}
	
	return -1; // error condition, we shouldn't ever get here
}

/* display_fixed_point_switch is used for reporting fixed-point switches since the internal
   representation numbers from zero and we want to report in the form
   presented by the user.
*/
void display_fixed_point_switch(matrix<int>& M, ostream& os, bool number_from_zero)
{
	int n = M.numcols();
	
	for (int i=0; i< n; i++)
	{
		for (int j=0; j< n; j++)
		{
			if (number_from_zero)			
				os << M[i][j];
			else
				os << M[i][j]+1;
			
			if (j < n-1)
				os << " ";
		}
			
		if (i < n-1)
			os << "  ";
	}
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
		
if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "valid_knotoid_input: no indication of knotoid head provided, returning false" << endl;
	
		return false;
	}
		
	
	int semi_arc;	
	int peer;


if (braid_control::DEBUG >= braid_control::BASIC)
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
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input: head " << 0 << " indicates first shortcut crossing is virtual" << endl;
		return false;
	}

if (braid_control::DEBUG >= braid_control::BASIC)
	debug << "valid_knotoid_input: head " << head << " lies on semi-arc " << semi_arc << " in component " << component << " with peer edge " << peer << endl;
	
	if (component != 0)
	{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input: invalid knotoid code, the knotoid leg must lie on the semi-arc numbered 0" << endl;
		return false;
	}

	if (semi_arc < peer && peer < code_data.num_component_edges[0])
	{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input: invalid knotoid code, peer edge also lies in the shortcut, which may not self-intersect" << endl;
		return false;
	}

    /* set the shortcut crossing flag for the head crossing */
    if (semi_arc % 2)
    {
		if (code_table[TYPE][head] == generic_code_data::TYPE1)				
		{
			shortcut_crossing[head] = -1; 
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input:   knotoid crosses shortcut at first shortcut crossing from the left" << endl;
		}
		else
		{
			shortcut_crossing[head] = 1;
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input:   knotoid crosses shortcut at first shortcut crossing from the right" << endl;
		}
	}
	else
	{
		if (code_table[TYPE][head] == generic_code_data::TYPE1)				
		{
			shortcut_crossing[head] = 1; 
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input:   knotoid crosses shortcut at first shortcut crossing from the right" << endl;
		}
		else
		{
			shortcut_crossing[head] = -1;
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input:   knotoid crosses shortcut at first shortcut crossing from the left" << endl;
		}
	}
		
	bool valid = true;

	for (int i = semi_arc+1; i< code_data.num_component_edges[0]; i++)
	{

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input: semi-arc " << i << endl;
	
		/* check that semi-arc i is the under-arc at the next crossing */
		if (i%2)
		{				
			int crossing = code_data.term_crossing[i];
			peer = code_table[EVEN_TERMINATING][crossing];

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input:   crossing " << crossing << ", peer " << peer << endl;
			
			/* peer edge must not lie in the shortcut and the label for this crossing must be '+' if odd arc is under-arc */
			if (semi_arc < peer && peer < code_data.num_component_edges[0])
			{
				valid = false;
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input: invalid knotoid code, peer edge also lies in the shortcut, which may not self-intersect" << endl;
				break;
			}
			else if (code_table[LABEL][crossing] != generic_code_data::POSITIVE)
			{
				valid = false;
if (braid_control::DEBUG >= braid_control::DETAIL)
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
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input:   OK, knotoid crosses from the left" << endl;
				}
				else
				{
					shortcut_crossing[crossing] = 1;
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input:   OK, knotoid crosses from the right" << endl;
				}
					
			}
		}
		else
		{
			int crossing = i/2 ;
			peer = code_table[ODD_TERMINATING][crossing];

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input:   crossing " << crossing << ", peer " << peer << endl;

			/* peer edge must not lie in the shortcut and the label for this crossing must be '-' if even arc is under-arc */
			if (semi_arc < peer && peer < code_data.num_component_edges[0])
			{
				valid = false;
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input: invalid knotoid code, peer edge also lies in the shortcut, which may not self-intersect" << endl;
				break;
			}
			else if (code_table[LABEL][crossing] != generic_code_data::NEGATIVE)
			{
				valid = false;
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input:   not negative, as required for valid knotoid" << endl;
				break;
			}
			else
			{
			    /* see comment above for the odd case */
				if (code_table[TYPE][crossing] == generic_code_data::TYPE1)				
				{
					shortcut_crossing[crossing] = 1; 
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input:   OK, knotoid crosses from the right" << endl;
				}
				else
				{
					shortcut_crossing[crossing] = -1;
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input:   OK, knotoid crosses from the left" << endl;
				}
			}
		}
	}
	
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "valid_knotoid_input: returning " << valid << endl;
	
	code_data.shortcut_crossing = shortcut_crossing;
	return valid;
}



/* flip_braid reverses the strand numbering of a braid whilst leaving the crossing type
   unchanged.  Thus, for a braid on n strands, crossing +/-s_i is changed to +/-s_{n-i}
   and t_i is changed to t_{n-i}.  This has the same effect as an ambient isotopy of R^3 
   that turns over the braid, when considered to be laid on the plane R^2.
*/
void flip_braid(string& braid)
{

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::flip_braid: flipping braid " << braid << endl;

	int num_terms;
	int num_strings;

	if (valid_braid_input(braid, num_terms, num_strings))
	{
		vector<int> braid_num(num_terms);
	    vector<int> type(num_terms);

		parse_braid(braid, num_terms, braid_num, type);

		ostringstream oss;
		
		for (int i=0; i< num_terms; i++)
		{
			if (type[i] == braid_crossing_type::POSITIVE)
				oss << "s";
			else if (type[i] == braid_crossing_type::NEGATIVE)
				oss << "-s";
			else 
				oss << "t";
			
			oss << num_strings - braid_num[i];
		}
		
		braid = oss.str();
		
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::flip_braid: to give braid " << braid << endl;

	}
	else
	{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::flip_braid: not a valid braid, doing nothing" << endl;
	}
	
}

/* invert_braid reverses the order of the braid generators and their signs.
   This has the same effect as reflecting the braid in a vertical line.
*/
void invert_braid(string& braid)
{

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::invert_braid: inverting (vertical line reflecting) braid " << braid << endl;

	int num_terms;
	int num_strings;

	if (valid_braid_input(braid, num_terms, num_strings))
	{
		vector<int> braid_num(num_terms);
	    vector<int> type(num_terms);

		parse_braid(braid, num_terms, braid_num, type);

		ostringstream oss;
		
		for (int i=num_terms-1; i >= 0; i--)
		{
			if (type[i] == braid_crossing_type::POSITIVE)
				oss << "-s";
			else if (type[i] == braid_crossing_type::NEGATIVE)
				oss << "s";
			else 
				oss << "t";
			
			oss << braid_num[i];
		}
		
		braid = oss.str();
		
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::invert_braid: to give braid " << braid << endl;

	}
	else
	{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::invert_braid: not a valid braid, doing nothing" << endl;
	}
	
}

/* plane_reflect_braid interchanges the sign of classical crossings.
   This has the same effect as reflecting the braid in the plane of the diagram.
   This is the same as Jeremy green's "vertical mirror"
*/
void plane_reflect_braid(string& braid)
{

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::plane_reflect_braid: reflecting braid " << braid << endl;

	int num_terms;
	int num_strings;

	if (valid_braid_input(braid, num_terms, num_strings))
	{
		vector<int> braid_num(num_terms);
	    vector<int> type(num_terms);

		parse_braid(braid, num_terms, braid_num, type);

		ostringstream oss;
		
		for (int i=0; i< num_terms; i++)
		{
			if (type[i] == braid_crossing_type::POSITIVE)
				oss << "-s";
			else if (type[i] == braid_crossing_type::NEGATIVE)
				oss << "s";
			else 
				oss << "t";
			
			oss << braid_num[i];
		}
		
		braid = oss.str();
		
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::plane_reflect_braid: to give braid " << braid << endl;

	}
	else
	{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::plane_reflect_braid: not a valid braid, doing nothing" << endl;
	}
	
}

/* line_reflect_braid reflects the braid in a horizontal line south of the braid diagram.
   It reverses the strand numbering of a braid and toggles the classical crossing types.
   Thus, for a braid on n strands, crossing +/-s_i is changed to -/+s_{n-i}
   and t_i is changed to t_{n-i}.  
*/
void line_reflect_braid(string& braid)
{

if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::line_reflect_braid: (horizontal) line reflecting braid " << braid << endl;

	int num_terms;
	int num_strings;

	if (valid_braid_input(braid, num_terms, num_strings))
	{
		vector<int> braid_num(num_terms);
	    vector<int> type(num_terms);

		parse_braid(braid, num_terms, braid_num, type);

		ostringstream oss;
		
		for (int i=0; i< num_terms; i++)
		{
			if (type[i] == braid_crossing_type::POSITIVE)
				oss << "-s";
			else if (type[i] == braid_crossing_type::NEGATIVE)
				oss << "s";
			else 
				oss << "t";
			
			oss << num_strings - braid_num[i];
		}
		
		braid = oss.str();
		
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::line_reflect_braid: to give braid " << braid << endl;

	}
	else
	{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::line_reflect_braid: not a valid braid, doing nothing" << endl;
	}
	
}

/* The Kamada double covering of a diagram is defined in "Colourings and doubled colourings of virtual doodles" for flat 
   virtuals and doodles, here implemented for braids.  Consider a braid reflected below within the plane, so that the 
   original strands 1,...,n and the reflected strands 1',...,n' are aligned as in the left side of the following left hand 
   diagram.  On the right side of the left hand diagram are the corresponding strand numbers of the double covering of the 
   braid.  Note that a crossing s_i in the original braid corresponds to the two crossings s_{n+i} and s_{n-i} in the 
   double cover.
   
                                                                                          s_i  
                   n         2n                   i+1 ---  --- n+i+1               i --    *    -- n+i
                                                         \/                            \       /
                   i         n+i                         /\                      i-1 ---o-----o--- n+i-1    
                                                  i --*--  --*-- n+i                     o   o
                   1         n+1                                                   1 -----o-o----- n+1   
                 -----------------                   B_i    B_i                            o
                   1'        n                                                     1 -----o-o----- n
                                                  i --*--  --*-- n-i+1                   o   o
                   i'        n-i+1                       \/                      i-1 ---o-----o--- n-i+2
                                                         /\                            /       \
                   n'        1                    i+1 ---  --- n-i                 i --    *    -- n-i+1
       
   
   The middle diagram shows the replacement of the cut system (the four points indicated by asterisks) around a flat crossing
   by two instances of the braid B_i.  That is, the replacement of crossing s_i in the Kamada double covering is 
   B_i s_{n+i} s_{n-i} B_i.

   The braid B_i is shown in the right diagram: we add the virtual crossings in the following order
    - top left quadrant, top down
    - bottom left quadrant, bottom up
    - centre crossing (t_n)
    - top right quadrant, bottom up
    - bottom right quadrant, top down
    
   The replacement of a virtual crossing t_i in the original braid is given by t_{n+i} t{n-i}.
*/
void Kamada_double_covering(string& braid)
{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::Kamada_double_covering: evaluating the double covering of braid " << braid << endl;

	int num_terms;
	int num_strings;

	if (valid_braid_input(braid, num_terms, num_strings))
	{
		vector<int> braid_num(num_terms);
	    vector<int> type(num_terms);

		parse_braid(braid, num_terms, braid_num, type);

		ostringstream oss;
		
		for (int i=0; i< num_terms; i++)
		{
			if (type[i] == braid_crossing_type::VIRTUAL)
			{
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "braid::Kamada_double_covering:   crossing " << i << " is virtual, replacing t" << braid_num[i] << " with "
	      << "t" << num_strings + braid_num[i] << " t" << num_strings - braid_num[i] << endl;
}

				oss << "t" << num_strings + braid_num[i] << "t" << num_strings - braid_num[i];
			}
			else 
			{
if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "braid::Kamada_double_covering:   crossing " << i << " is flat, replacing s" << braid_num[i] << " with "
	      << "B_i s" << num_strings + braid_num[i] << " s" << num_strings - braid_num[i] << " B_i" << endl;
}

				ostringstream B_i;
				
				for (int j=braid_num[i]-1; j >=1; j--)
				    B_i << "t" << num_strings + j;
				    
				for (int j=braid_num[i]-1; j >=1; j--)
				    B_i << "t" << num_strings - j; 

				B_i << "t" << num_strings;

				for (int j=1; j <=braid_num[i]-1; j++)
				    B_i << "t" << num_strings + j;
				    
				for (int j=1; j <=braid_num[i]-1; j++)
				    B_i << "t" << num_strings - j; 


if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << "braid::Kamada_double_covering:     where B_i = " << B_i.str() << endl;
}				
				oss << B_i.str() << "s" << num_strings + braid_num[i] << "s" << num_strings - braid_num[i] << B_i.str();
			}
		}
		
		braid = oss.str();
		
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::Kamada_double_covering: to give braid " << braid << endl;

	}
	else
	{
if (braid_control::DEBUG >= braid_control::DETAIL)
	debug << "braid::Kamada_double_covering: not a valid braid, doing nothing" << endl;
	}
}

/*  The commutative automorphism invariant determines the codimension 0 invariant (Delta_0) and, where necessary, the codimension 1
    invariant (Delta_1) of the R-module representation of the given labelled peer code determined by the F-biquandle defined from
    a pair of commuatative automorphism, \phi and \psi, of an abelian group.  That is F(b,a) = (a_b, b^a), where 
    a^b = \phi(a)-\phi(b)+\psi(b) and a_b = \psi(a).
    
    Note that although there are many similarities to R_module_rep in the code here, R_module_rep is using equations determined 
    by a switch and here we use equations determined by the sideways map of a switch.
*/
void commutative_automorphism_invariant(const Qpmatrix& phi, const Qpmatrix& psi, string input_string, string title)
{

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "commutative_automorphism_invariant: Calculating commutative automorphism invariant for " << title << " = " << input_string << endl;

	if (title.length())
   	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\n\n" << title << "\n" << input_string << endl;
		if (!braid_control::RAW_OUTPUT)
			output << "\n\n" << title << "\n" << input_string << endl;
   	}		
   		
if (braid_control::DEBUG >= braid_control::SUMMARY)
{
	debug << "commutative_automorphism_invariant: phi = " << endl;
	print(phi,debug,3,"commutative_automorphism_invariant: ");
	debug << "commutative_automorphism_invariant: psi = " << endl;
	print(psi,debug,3,"commutative_automorphism_invariant: ");
}

	generic_code_data code_data;
	
	if (input_string.find('(') != string::npos || input_string.find('[') != string::npos)
	{
		read_code_data (code_data,input_string);
	}
	else
	{
		int num_terms;
		int num_strings;
		if (valid_braid_input(input_string, num_terms, num_strings))
		{
			input_string = braid_to_generic_code(input_string, num_terms, generic_code_data::peer_code);
			read_peer_code(code_data, input_string);
			
if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "commutative_automorphism_invariant: valid braid input converted to peer code " << input_string << endl;

		}
		else
		{
			return; // invalid braid.
		}
	}
	
	matrix<int>& code_table = code_data.code_table; 
	int num_crossings = code_data.num_crossings;
	int num_components = code_data.num_components;
	int num_classical = num_crossings;

	for (int i=0; i<num_crossings; i++)
	{
		if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
			num_classical--;
	}

	/* If num_classical == 0 we have been given a code containing only virtual crossings.
	   This yields a 0x0 matrix representation and a degenerate R-module.  This situation 
	   is indicated by setting the Matrix_rep pointer to zero.
	*/
	if (num_classical == 0)
	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\nError! Labelled peer and immersion codes must contain at least one real crossing" << endl;
		return;
	}

	/* Determine the size of the R_module representation.  In the case of a knot or link, this is a 
	   2*N*num_classical square matrix where the switch_matrix elements are NxN matrices (for the 
	   Alexander and Quaternionic reps, N=1).  In the case of a long knot, the R_module representation
	   is a 2*N*num_classical x (2*N*num_classical + N) matrix, allowing for the extra variable.  
	*/
	int N = phi.numrows(); // N is used later, it is what is referred to elsewhere as the switch_matrix_N_factor
	int matrix_N_size = 2*num_classical;  // number of N-rows
	int matrix_size = N*matrix_N_size;  // number of actual rows

if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "commutative_automorphism_invariant: num_components = " << num_components << endl;
	debug << "commutative_automorphism_invariant: num_crossings = " << num_crossings << endl;
	debug << "commutative_automorphism_invariant: num_classical = " << num_classical << endl;
	debug << "commutative_automorphism_invariant: N (switch_matrix_N_factor) = " << N << endl;
	debug << "commutative_automorphism_invariant: matrix_N_size (number of N-rows) = " << matrix_N_size << endl;
	debug << "commutative_automorphism_invariant: matrix_size (number of actual rows) = " << matrix_size << endl;
}

	/* The equations for virtual crossings state x_2i = x_{2i+1} and x_{2j-1} = x_2j, where the terminating
	   edges at the crossing are 2i and 2j-1.  Thus, the variables corresponding to the incoming edges
	   are trivial.
	*/
	valarray<int> trivial_variable_flags(2*num_crossings); 
	valarray<int> variable(2*num_crossings);               
	
	
	for (int i=0; i<2*num_crossings; i++)
		trivial_variable_flags[i] = 0;
	
	for (int i=0; i<num_crossings; i++)
	{
		if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
		{
if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
	debug << "commutative_automorphism_invariant: crossing " << i << " trivial variables = " << 2*i << " and " << code_table[OPEER][i] << endl;
			trivial_variable_flags[2*i] = 1;
			trivial_variable_flags[code_table[OPEER][i]] = 1;
		}
	}
	
if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "commutative_automorphism_invariant: trivial_variable_flags: ";
	for (int i=0; i<2*num_crossings; i++)
		debug << trivial_variable_flags[i] << ' ';
	debug << endl;
}	

	/* now we can determine the variable permutation between the variables
	   assigned to the semi-arcs between all crossings, and those assigned
	   to semi-arcs between classical crossings, ignoring virtual crossings.
	*/
	int count = 0;
	vector<int>& first_edge_on_component = code_data.first_edge_on_component;
	vector<int>& num_component_edges = code_data.num_component_edges;
	
	for (int i=0; i< num_components; i++)
	{
		int first_edge = first_edge_on_component[i];
		for (int j=0; j < num_component_edges[i]; j++)
		{
			variable[first_edge+j] = first_edge+j-count;
			count += trivial_variable_flags[first_edge+j];
		}
		
		/* We have to deal with the case that the last k crossings on the 
		   component are virtual crossing, in which case the last k variables 
		   on the component are the same as the variable corresponding to the first edge.
		   In the degenerate case where all the crossings on a component are virtual
		   this will result in all the edges of the component being associated with 
		   the same variable, being the variable belonging to the first edge of the 
		   next component.  This won't matter as we'll not be creating any relations
		   for these virtual crossings.
		*/
		for (int j=num_component_edges[i]-1; j>=0 && trivial_variable_flags[first_edge+j] == 1; j--)
			variable[first_edge+j] = variable[first_edge];
		
	}
	
	
if (braid_control::DEBUG >= braid_control::BASIC)
{
	debug << "commutative_automorphism_invariant: variables: ";
	for (int i=0; i<2*num_crossings; i++)
		debug << variable[i] << ' ';
	debug << endl;
}	
    /* Now create the matrix and evaluate the representation of the knot based upon the reduced variable set.
	   Each real crossing determines two rows of the matrix representation, representing one of the pairs of 
	   equations above, dependent on the crossing type.  We work through the real crossings adding a pair of 
	   rows to the matrix representation for each one.  (Note we are interested in the R-module generated by
	   this matrix, so the order in which we construct the rows is irrelevant.)
    */
	
	Qpmatrix matrix_rep(matrix_size,matrix_size);
	
	for (unsigned int i=0; i<matrix_rep.numrows(); i++)
	for (unsigned int j=0; j<matrix_rep.numcols(); j++)
		matrix_rep[i][j] = Qpolynomial("0");

	/*  The matrix_rep equations are determined from the fact that a commuatative automorphism switch is the
	    sideways map of an F-biquandle F(b,a) = (a_b, b^a), where a^b = \phi(a)-\phi(b)+\psi(b) and 
	    a_b = \psi(a).
	    
                     a_b   b^a          a_b   b^a            b^a   a_b          b^a   a_b
                  2j-1 \   / 2i+1      2i \   / 2j        2j-1 \   / 2i+1       2i\   / 2j
                        \ /                \ /	                \ /                \ /
                         \                  \                    /                  /
                        / \                / \                  / \                / \
                    2i /   \ 2j      2j-1 /   \	2i+1        2i /   \ 2j       2j-1/   \ 2i+1
                      b     a            b     a              a     b            a     b
	    
	    By considering a positive type I and type II crossing, both of which respect the above map F we see
	    that the corresponding equations are
	    
	    
	                                \psi(x_{2j}=x_{2j-1}                                    \psi(x_{2i+1}=x_{2i}
	    \phi(x_{2i})-\phi(x_{2j})+\psi(x_{2j}))=x_{2i+1}    \phi(x_{2j-1})-\phi(x_{2i+1})+\psi(x_{2i+1}))=x_{2j}
	
	    A negative type I and type II crossing respects map F^*(a,b) = (b^a,a_b) so the corresponding equations 
	    are
	    
	    
	    \phi(x_{2j})-\phi(x_{2i})+\psi(x_{2i}))=x_{2j-1}    \phi(x_{2i+1})-\phi(x_{2j-1})+\psi(x_{2j-1}))=x_{2i}
	                                \psi(x_{2i}=x_{2i+1}                                    \psi(x_{2j-1}=x_{2j}
	
	   
	   Note that for the up operation we have (\psi-\phi)(x_p) for some p.
	   
	   There is a possibility that we have to handle the case where F(x_p,x_q) = (x_p,x_r), F(x_p,x_q) = (x_s,x_q),
	   or similarly for F^*, due to a Reidemeister I detour.  We therefore add the (negative) terms corresponding 
	   to the right hand side of the above equations to the matrix_rep entry, rather than setting that entry.
	   	
	*/

	Qpmatrix psi_minus_phi = psi - phi;
	
	
	int row = 0;
	
	for (int i=0; i<num_crossings; i++)
	{
		/* only consider classical crossings */
		if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
			continue;
		
		/* determine the braid crossing type, POSITIVE, NEGATIVE */
		int crossing_type;		
		if (  (code_table[TYPE][i] == generic_code_data::TYPE1 && code_table[LABEL][i] == generic_code_data::POSITIVE)
				||(code_table[TYPE][i] == generic_code_data::TYPE2 && code_table[LABEL][i] == generic_code_data::NEGATIVE)
   					)
			crossing_type = braid_crossing_type::NEGATIVE;
		else
			crossing_type = braid_crossing_type::POSITIVE;

if (braid_control::DEBUG >= braid_control::DETAIL)
{	
	debug << "commutative_automorphism_invariant: crossing " << i;
	if ( crossing_type == braid_crossing_type::NEGATIVE )
		debug << " negative";
	else
		debug << " positive";

	debug <<  ", TYPE " << (code_table[TYPE][i] == generic_code_data::TYPE1? "I" : "II") << endl;
	debug << "commutative_automorphism_invariant: peer code variables p[i] = ";// << code_table[PERM][i];
	debug << "\t2i = " << 2*i << "\t\t2j-1 = " << code_table[ODD_TERMINATING][i];
	debug << "\t2i+1 = " << 2*i+1 << "\t2p[i] = " << code_table[EVEN_ORIGINATING][i] << endl;
	
	debug << "commutative_automorphism_invariant: mapped to\t\t";//v(p[i]) = " << variable[code_table[PERM][i]];
	debug << "\tv(2i) = " << variable[2*i] << "\tv(2p[i]-1) = " << variable[code_table[ODD_TERMINATING][i]];
	debug << "\tv(2i+1) = " << variable[2*i+1] << "\tv(2p[i]) = " << variable[code_table[EVEN_ORIGINATING][i]] << endl;

	debug << "commutative_automorphism_invariant: row " << 2*row << " ";
}		
		/* first do the up action in row 2*row. */
		if (crossing_type == braid_crossing_type::POSITIVE)
		{
			if (code_table[TYPE][i] == generic_code_data::TYPE1)
			{
				set_matrix_N_element(matrix_rep,2*row,code_table[EVEN_TERMINATING][i],phi,0,0,N);
				set_matrix_N_element(matrix_rep,2*row,code_table[EVEN_ORIGINATING][i],psi_minus_phi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row,code_table[ODD_ORIGINATING][i],N);

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << code_table[EVEN_TERMINATING][i] << ":phi    ";
	debug << code_table[EVEN_ORIGINATING][i] << ":psi_minus_phi    ";
	debug << code_table[ODD_ORIGINATING][i] << ":-1     " << endl;
}
			}
			else
			{
				set_matrix_N_element(matrix_rep,2*row,code_table[ODD_TERMINATING][i],phi,0,0,N);
				set_matrix_N_element(matrix_rep,2*row,code_table[ODD_ORIGINATING][i],psi_minus_phi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row,code_table[EVEN_ORIGINATING][i],N);

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << code_table[ODD_TERMINATING][i] << ":phi    ";
	debug << code_table[ODD_ORIGINATING][i] << ":psi_minus_phi    ";
	debug << code_table[EVEN_ORIGINATING][i] << ":-1     " << endl;
}
			}
		}
		else
		{
			if (code_table[TYPE][i] == generic_code_data::TYPE1)
			{
				set_matrix_N_element(matrix_rep,2*row,code_table[EVEN_ORIGINATING][i],phi,0,0,N);
				set_matrix_N_element(matrix_rep,2*row,code_table[EVEN_TERMINATING][i],psi_minus_phi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row,code_table[ODD_TERMINATING][i],N);

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << code_table[EVEN_ORIGINATING][i] << ":phi    ";
	debug << code_table[EVEN_TERMINATING][i] << ":psi_minus_phi    ";
	debug << code_table[ODD_TERMINATING][i] << ":-1     " << endl;
}
			}
			else
			{
				set_matrix_N_element(matrix_rep,2*row,code_table[ODD_ORIGINATING][i],phi,0,0,N);
				set_matrix_N_element(matrix_rep,2*row,code_table[ODD_TERMINATING][i],psi_minus_phi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row,code_table[EVEN_TERMINATING][i],N);

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << code_table[ODD_ORIGINATING][i] << ":phi    ";
	debug << code_table[ODD_TERMINATING][i] << ":psi_minus_phi    ";
	debug << code_table[EVEN_TERMINATING][i] << ":-1     " << endl;
}
			}
		}
		
   
if (braid_control::DEBUG >= braid_control::DETAIL) 
	debug << "commutative_automorphism_invariant: row " << 2*row+1 << " ";

		/* now the down action in row 2*row+1 */	
		if (crossing_type == braid_crossing_type::POSITIVE)
		{
			if (code_table[TYPE][i] == generic_code_data::TYPE1)
			{
				set_matrix_N_element(matrix_rep,2*row+1,code_table[EVEN_ORIGINATING][i],psi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row+1,code_table[ODD_TERMINATING][i],N);

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << code_table[EVEN_ORIGINATING][i] << ":psi    ";
	debug << code_table[ODD_TERMINATING][i] << ":-1     " << endl;
}
			}
			else
			{
				set_matrix_N_element(matrix_rep,2*row+1,code_table[ODD_ORIGINATING][i],psi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row+1,code_table[EVEN_TERMINATING][i],N);

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << code_table[ODD_ORIGINATING][i] << ":psi    ";
	debug << code_table[EVEN_TERMINATING][i] << ":-1     " << endl;
}
			}
		}
		else
		{
			if (code_table[TYPE][i] == generic_code_data::TYPE1)
			{
				set_matrix_N_element(matrix_rep,2*row+1,code_table[EVEN_TERMINATING][i],psi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row+1,code_table[ODD_ORIGINATING][i],N);

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << code_table[EVEN_TERMINATING][i] << ":psi    ";
	debug << code_table[ODD_ORIGINATING][i] << ":-1     " << endl;
}
			}
			else
			{
				set_matrix_N_element(matrix_rep,2*row+1,code_table[ODD_TERMINATING][i],psi,0,0,N);
				decrement_matrix_N_element(matrix_rep,2*row+1,code_table[EVEN_ORIGINATING][i],N);

if (braid_control::DEBUG >= braid_control::DETAIL)
{
	debug << code_table[ODD_TERMINATING][i] << ":psi    ";
	debug << code_table[EVEN_ORIGINATING][i] << ":-1     " << endl;
}
			}
		}
	
		row++;
	}
	

	if (braid_control::EXTRA_OUTPUT)
	{
		if (!braid_control::SILENT_OPERATION)
		{
			cout << "\n" << matrix_size << "-square matrix representation of R-module determined by peer code, M:\n";
			cout << matrix_rep << endl;
		}
		
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\nMatrix representation of R-module determined by peer code, M:\n";
			output << matrix_rep << endl;
		}
	}

if (braid_control::DEBUG >= braid_control::BASIC)
{
    debug << "commutative_automorphism_invariant: matrix representation derived from generic code data:" << endl;
    debug << matrix_rep;
    debug << endl;
}


	if (braid_control::WAIT_SWITCH)
	{
	    matrix_control::WAIT_INFO = true;
	    matrix_control::wait_threshold = braid_control::wait_threshold;
	    matrix_control::wait_count = 0;

	    /* the wait_threshold comes from the command line */
		if (!braid_control::SILENT_OPERATION)
			cout << "\nCalculating polynomials, please wait\n";
	}
	else
	    matrix_control::WAIT_INFO = false;
		    
	Qpolynomial delta_0 = determinant (matrix_rep, title);

	if (!braid_control::SILENT_OPERATION)
		cout << "\nDelta_0 = " << delta_0 << endl;
	if (!braid_control::RAW_OUTPUT)
	{
		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
		output << "Delta_0 = ";
		if (braid_control::OUTPUT_AS_INPUT) 
			output << "\n";
	}
	output << delta_0 << endl;

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "commutative_automorphism_invariant:Delta_0 = " << delta_0 << endl;

	if (!braid_control::CALCULATE_DELTA_0_ONLY && (braid_control::ALWAYS_CALCULATE_DELTA_1 || delta_0 == Qpolynomial("0")))
	{
		/* hcf is used only in single variable examples */
		polynomial<scalar,char,bigint> hcf = delta_0.getn();

if (braid_control::DEBUG >= braid_control::INTERMEDIATE)
    debug << "\ncommutative_automorphism_invariant:hcf initialized to " << hcf << endl;

		int rperm[matrix_size];
		int cperm[matrix_size];

		/* take out a (N) row and (N) column from the matrix representation 
		   and calculate the determinant of what is left
		*/
		for (int Nr1=0;Nr1<matrix_N_size;Nr1++)
		{
			/* take out the N rows of the matrix rep starting from row Nr1*N row */
			for (int i=0; i<Nr1*N; i++)
    			rperm[i] = i;
			for (int i = Nr1*N+N; i< matrix_size; i++)
				rperm[i-N] = i;


			int Nc2=0; // needed for rat_minor_determinant - relates to long knots;
			for (int Nc1= 0; Nc1< matrix_N_size; Nc1++)
			{
			
				Qpolynomial det = polynomial<scalar,char,bigint> ("0");
				
				/* take out the N columns of the matrix rep starting from column N*j */
				for (int i=0; i<Nc1*N; i++)
	   				cperm[i] = i;

				for (int i=Nc1*N+N; i<matrix_size; i++)
					cperm[i-N] = i;
				
				if (rat_minor_determinant(&matrix_rep, N, Nr1, Nc1, Nc2, rperm, cperm, title, hcf))
					goto rat_poly_delta_1_calculation_complete;
			}			
		}

rat_poly_delta_1_calculation_complete:

		if (!braid_control::SILENT_OPERATION)
			cout << endl;
		
		if (braid_control::NUMERATOR_GCD)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "Delta_1 = " << hcf << endl;	
			
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "numerator gcd  = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << hcf << endl;

if (braid_control::DEBUG >= braid_control::SUMMARY)
	debug << "commutative_automorphism_invariant:" << "Delta_1 = " << hcf << endl;	

		}
	}
}
