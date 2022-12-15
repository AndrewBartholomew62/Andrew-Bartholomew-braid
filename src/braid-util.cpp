/************************************************************************
string invstr(string& str)
void flip_braid(string& braid)
void invert_braid(string& braid)
void plane_reflect_braid(string& braid)
void line_reflect_braid(string& braid)
void parse_braid(string word, int num_terms, vector<int>& braid_num, vector<int>& type)
void convert_rep (char*& word, bool silent)
bool valid_braid_input (string& input_string, int& num_terms, int& num_strings, bool raw_output, bool silent_output, bool output_as_input)
void parse_braid(string word, int num_terms, vector<int>& braid_num, vector<int>& type);
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
#include <debug-control.h>
#include <braid-util.h>


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

/* flip_braid reverses the strand numbering of a braid whilst leaving the crossing type
   unchanged.  Thus, for a braid on n strands, crossing +/-s_i is changed to +/-s_{n-i}
   and t_i is changed to t_{n-i}.  This has the same effect as an ambient isotopy of R^3 
   that turns over the braid, when considered to be laid on the plane R^2.
*/
void flip_braid(string& braid)
{

if (debug_control::DEBUG >= debug_control::DETAIL)
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
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::flip_braid: to give braid " << braid << endl;

	}
	else
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::flip_braid: not a valid braid, doing nothing" << endl;
	}
	
}

/* invert_braid reverses the order of the braid generators and their signs.
   This has the same effect as reflecting the braid in a vertical line.
*/
void invert_braid(string& braid)
{

if (debug_control::DEBUG >= debug_control::DETAIL)
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
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::invert_braid: to give braid " << braid << endl;

	}
	else
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::invert_braid: not a valid braid, doing nothing" << endl;
	}
	
}

/* plane_reflect_braid interchanges the sign of classical crossings.
   This has the same effect as reflecting the braid in the plane of the diagram.
   This is the same as Jeremy green's "vertical mirror"
*/
void plane_reflect_braid(string& braid)
{

if (debug_control::DEBUG >= debug_control::DETAIL)
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
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::plane_reflect_braid: to give braid " << braid << endl;

	}
	else
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
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

if (debug_control::DEBUG >= debug_control::DETAIL)
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
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::line_reflect_braid: to give braid " << braid << endl;

	}
	else
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::line_reflect_braid: not a valid braid, doing nothing" << endl;
	}
	
}

void parse_braid(string word, int num_terms, vector<int>& braid_num, vector<int>& type)
{
	char* braid_word = c_string(word);
	char* cptr = braid_word;
	
	for ( int i = 0; i< num_terms; i++)
	{
	    if (*cptr == '-')
	    {
			type[i] = braid_crossing_type::NEGATIVE;
			cptr++;
	    }
	    else if (*cptr == 't' || *cptr == 'T')
			type[i] = braid_crossing_type::VIRTUAL;
        else
			type[i] = braid_crossing_type::POSITIVE;

        cptr++;
        char* mark = cptr; /* mark where we start the number */

        /* look for the end of the number */
        while (isdigit(*cptr))
            cptr++;

		int number;
	    get_number(number, mark);
	    braid_num[i] = number;
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "parse_braid: braid_num: ";
	for (int i=0; i<num_terms; i++)
		debug << braid_num[i] << ' ';
	debug << "\nparse_braid: type: ";
	for (int i=0; i<num_terms; i++)
		debug << type[i] << ' ';
	debug << endl;
}    

	delete[] braid_word;
}

/* convert_rep takes a braid word in the standard representation using a,b,c,A,B,C etc.
   and converts it into si, -si format.
*/
void convert_rep (string& word, bool silent)
{
    bool inverse;
	ostringstream oss;
	
	for (size_t i=0; i< word.length(); i++)
	{
	    if (!isalpha(word[i]))
		{
			if (!silent)
				cout << "Error in input: non-alphabetical character " << word[i] << " \n";
			word = "";
			return;
		}
		else
		{
			if (word[i] >= 'a')
				inverse = false;
			else
				inverse = true;

			/* write this term into local */
			if (inverse)
				oss << '-';
			oss << 's';
			oss << word[i] - (inverse?'A':'a')+1;
		}
	    
	}

	word = oss.str();
}

/* valid_braid_input returns a boolean indicating whether input is a valid braid, 
   num_terms and num_strings may or may not have been set if the input is invalid 
   
   If the braid includes the doodle qualifier, the function checks that only positive
   classical crossings have been specified.
   
   input_string needs to be a reference, since we might be converting the representation if we're presented with an abcABC format braid
*/
bool valid_braid_input (string& input_string, int& num_terms, int& num_strings, bool raw_output, bool silent_output, bool output_as_input)
{
	char*       mark;
	int         number;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::valid_braid_input: presented with input_string " << input_string << endl;
	
	/* first check and then remove any unwanted qualifiers from the input string */
	bool doodle_braid = false;
	string::size_type pos = input_string.find('{');
	
	if (pos != string::npos)
	{
		
   		if (input_string.substr(pos).find("doodle") != string::npos)
   			doodle_braid = true;

		input_string = input_string.substr(0,pos);
	}

	pos = input_string.find('-');
	if (pos != string::npos && doodle_braid)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::valid_braid_input: doodle braid contains '-' character, only positive crossings should be specified" << endl;
		
		return false;
	}

	if (input_string.find('s') == string::npos && input_string.find('S') == string::npos && input_string.find('t') == string::npos && input_string.find('T') == string::npos)
	{
		convert_rep(input_string,silent_output); //assume its a traditional representation of a classical braid
		if (input_string.length() ==0)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "braid::valid_braid_input: convert_rep found error in input " << endl;

			return false; // convert_rep has found an error
		}
		else
		{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "braid::valid_braid_input: converted to = " << input_string << endl;
		}
	}
	
	char* inbuf = c_string(input_string);
	
    /* evaluate the number of terms in the word by counting the number of
       s and t characters - we scan the input at this point
    */
    char* cptr = inbuf;
    num_terms = 0;
	do
    {
        if (isalpha(*cptr))
        {
            num_terms++;
            if (isalpha(*(cptr+1)) || *(cptr+1) == '\0')
            {
				if (!silent_output)
					cout << "Error in input: must specify suffix for " << *cptr << "\n";

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "braid::valid_braid_input: Error in input: must specify suffix for " << *cptr << endl;

                return false;
            }
        }
        else if (*cptr == '-')
        {
            if (!isalpha(*(cptr+1)) || *(cptr+1) == '\0')
            {
                if (*(cptr+1) == '\0')
				{
					if (!silent_output)
						cout << "Error in input: word ends in '-'!\n";

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "braid::valid_braid_input: Error in input: word ends in '-'!" << endl;
				}
                else
				{
					if (!silent_output)
						cout << "Error in input: must specify suffix for " << *(cptr-1) << "\n";

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "braid::valid_braid_input: Error in input: must specify suffix for " << *(cptr-1) << endl;
				}
                return false;
            }
            else if (*(cptr+1) == 't' || *(cptr+1) == 'T')
            {
				if (!silent_output)
					cout << "Error in input: t crossings cannot be negative\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "braid::valid_braid_input: Error in input: t crossings cannot be negative" << endl;
				return false;
            }
        }
    } while (*++cptr != '\0');

	/* Work out how many strings we have in the braid */
	num_strings = 0;
	cptr = inbuf;
	for ( int i = 0; i< num_terms; i++)
	{
		if (*cptr == '-')
			cptr++;

		cptr++;
		mark = cptr; /* mark where we start the number */

		/* look for the end of the number */
		while (isdigit(*cptr))
			cptr++;

    	get_number(number, mark);
			if (number > num_strings)
			num_strings = number;
	}
	num_strings++;
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
    debug << "braid::valid_braid_input: number of strands = " << num_strings << endl;
    debug << "braid::valid_braid_input: number of terms = " << num_terms << endl;
}
		
	delete[] inbuf;
	return true;
}

