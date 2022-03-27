/************************************************************************
void get_input_word (ifstream& input, string& buffer, string& title)
**************************************************************************/
#include <fstream>
#include <cstring>

using namespace std;

extern ifstream     input;
extern ofstream     debug;

#include <util.h>
#include <debug-control.h>
#include <input-control.h>

unsigned int input_control::ACCEPT_MAP = 0;

string invstr(string& str);

/* get_input_word retrieves the next braid word, Gauss code labelled peer code or 
   labelled immersion code from the input file, carrying out any manipulation 
   required, places the resultant word in buffer and any title string in title.  
   The function does not remove the '--' leader from the title and sets the title 
   to the empty string if none is found.  Any qualifiers are appended to the 
   input string, the function does not remove the braces {} from around qualifiers.
*/
void get_input_word (ifstream& input, string& buffer, string& title)
{
    string 	next_line;
	string qualifiers;
    string 	loc_buf;
	string	w1,w2,w3,w4,w5,w6,w7,w8;
	string* target;
	char*	lptr;
    bool   	word_found = false;
	bool	accepted_braid_input = false;
	bool	accepted_non_braid_input = false;
    bool	not_finished;

    /* start by setting the buffer and title to the empty string */
    buffer.clear();
    title.clear();

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
    debug << "get_input_word: accepting ";
    if (input_control::ACCEPT_MAP & input_control::braid_word)
		debug << "braid words ";
    if (input_control::ACCEPT_MAP & input_control::immersion_code)
		debug << "immersion codes ";
    if (input_control::ACCEPT_MAP & input_control::peer_code)
		debug << "peer codes ";
    if (input_control::ACCEPT_MAP & input_control::gauss_code)
		debug << "Gauss codes ";
    if (input_control::ACCEPT_MAP & input_control::planar_diagram)
		debug << "planar diagrams ";
    if (input_control::ACCEPT_MAP & input_control::lace_code)
		debug << "lace codes ";
    debug << endl;
}
	
	while (!word_found && getline(input, next_line))
    {		

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "get_input_word: read line: " << next_line << endl;

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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "get_input_word: read qualifiers " << qualifiers << " from line " << next_line << endl;

			next_line = next_line.substr(0,pos);
			
		}

		char* line_buf = c_string(next_line);
				
	    /* if this line contains a switch, or programme options ignore it */
	    if ( strchr(line_buf,'[') && !strchr(line_buf,'/') && !strchr(line_buf,'\\') && !strchr(line_buf,'X'))
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "get_input_word: line contains programme options, ignoring line" << endl;
			goto done_with_line;
		}
		else if (  (strchr(line_buf,'s') || strchr(line_buf,'S')) && 
			        strchr(line_buf,'=') && 
				    !(strchr(line_buf,'w') && isdigit(*(strchr(line_buf,'w')+1))) //not a braid statement
			    )
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "get_input_word: line contains a switch, ignoring line" << endl;
			goto done_with_line;
		}

		
	    /* is there any whitespace at the beginning of the line? */
	    lptr = line_buf;
	    while (isspace(*lptr))
			lptr++;

	    if (strlen(lptr) && *lptr != '\n')
	    {
			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "get_input_word: next word line = " << lptr << endl;

			/* check for a title line */
			char* cptr = lptr;
			if (*cptr == '-' && *(cptr+1) == '-')
			{
		    	title = lptr;

				goto done_with_line;
			}

			/* the target of the line parsing is either one of the word buffers w1-w8, or is 'buffer', if the input is a braid statement, 
			   a peer code, an immersion code, a Gauss code, a planar diagram or a lace, .  
			   
			    - an immersion code is indicated by the presence of a '(' character but no '|'
			    - a peer code is indicated by the presence of a '[' character
			    - a Gauss code is indicated by the presence of a '/', an 'O' or '\' character but no '(' or '['.	
			    - a braid word is indicated by the presents of an 's','S', 't' or 'T'
			    - a planar diagram is indicated by the presence of an 'X'
			    - a lace code is indicated by the presence of a '|'
			*/
			if (!accepted_non_braid_input && !accepted_braid_input)
			{
				if (strchr(lptr,'(') && !strchr(lptr,'|') && (input_control::ACCEPT_MAP & input_control::immersion_code))
				{
					accepted_non_braid_input = true;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "get_input_word: detected acceptable start of immersion code" << endl;
				}
				else if (strchr(line_buf,'[') && (input_control::ACCEPT_MAP & input_control::peer_code))
				{
					accepted_non_braid_input = true;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "get_input_word: detected acceptable start of peer code" << endl;
				}
				else if ((strchr(lptr,'/') || strchr(lptr,'O') || strchr(lptr,'\\')) && (input_control::ACCEPT_MAP & input_control::gauss_code))
				{
					accepted_non_braid_input = true;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "get_input_word: detected acceptable start of Gauss code" << endl;
				}
				else if ((strchr(line_buf,'s') || strchr(line_buf,'S') || 
				          strchr(line_buf,'t') || strchr(line_buf,'T')) && (input_control::ACCEPT_MAP & input_control::braid_word))
				{
					accepted_braid_input = true;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "get_input_word: detected acceptable start of braid word" << endl;
				}
				else if (strchr(line_buf,'X') && (input_control::ACCEPT_MAP & input_control::planar_diagram))
				{
					accepted_non_braid_input = true;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "get_input_word: detected acceptable start of planar diagram" << endl;
				}
				else if (strchr(lptr,'|') && (input_control::ACCEPT_MAP & input_control::lace_code))
				{
					accepted_non_braid_input = true;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "get_input_word: detected acceptable start of lace code" << endl;
				}
				else
					goto done_with_line;
			}
			
			/* Here we have the kind of input we're looking for on the next_line */
			if ( accepted_non_braid_input && 
			     ( (input_control::ACCEPT_MAP & input_control::gauss_code) || 
			       (input_control::ACCEPT_MAP & input_control::immersion_code) || 
			       (input_control::ACCEPT_MAP & input_control::peer_code) ||
			       (input_control::ACCEPT_MAP & input_control::planar_diagram) ||
			       (input_control::ACCEPT_MAP & input_control::lace_code)
			     ) 
			   )
			{
			    /* Take out line escapes and build up either the
				   labelled immersion code or the Gauss code in buffer
			    */
			    bool escape_present = false;
			    cptr = strchr(lptr,'\\');
			    if (cptr)
			    {
					escape_present = true;
					*cptr = ' ';
				}

			    /* copy the line into buffer and decide if there's more to come */
				buffer += string(lptr);

			    /* This test means we cannot break lines after
				   the '/' character in the input file
				*/
				if (strchr(lptr,'/') || strchr(lptr,'|') || strchr(lptr,'O') || (strchr(lptr,'X') && !escape_present) )
				{
					word_found = true;					
					
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "get_input_word: adding qualifiers " << qualifiers << " to buffer " << buffer << endl;;

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
			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "get_input_word: looking for a braid word...";

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
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "line is an assignment statement" << endl;
			}
			else // braid statement
			{
			    target = &buffer;
		    	word_found = true;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
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
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "get_input_word: adding qualifiers " << qualifiers << " to buffer " << buffer << endl;;

				buffer += qualifiers;
			}
	    }
done_with_line:
		delete[] line_buf;
    } //end of while (getline(input, next_line) && !word_found)

    if (!word_found)
    	buffer = "exit";
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "get_input_word: returning buffer " << buffer << endl;;
		
}
