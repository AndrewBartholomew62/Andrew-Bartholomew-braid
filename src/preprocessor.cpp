/**************************************************************************
  A simple pre-processor for braid programme input files, designed to support
  systematic, modular testing.  The function is defined in a separate file 
  to allow for the building of a standalone preprocessor for regression testing 
  with old versions.
  
  Supported directives:
    #common   options to be included for all test cases
    #define   macro definition for directory specification as in makefiles
			  #define source /home/bart/source/braid
			  #define textfiles ($source)/textfiles
    #end      stop including the current file at this point
    #include[{n}] <filename> specifies the name of a file to include at this point of the input
                  the file <filename> may reside in a subdirectory for example #include fred or 
                  #include subdir/fred.  If the optional {n} parameter is provided, only the first 
                  n lines of <filename> are included.
                  
	#result <name> <result-file> the <result-file> file name is expanded as for include files
	#script <test-script> the <test-script> file name is expanded as for include files
	#test <name> [option,option,...]
    
  Note that #end and exit are different.  Lines containing exit are included in the
  output and may cause the braid programme to terminate processing, if the line is
  not commented out.  The #end directive simply aborts the #include.
              
 **************************************************************************/
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <list>
#include <map>

using namespace std;
#include <debug-control.h>
#include <util.h>

extern ofstream     debug;
int recursion_level = 0;

string expand_string (string input_string,map<string,string>& var_map);

/* max_lines == -2 indicates an override to include a single line from each include file */
string preprocess_file (ostream& os, ifstream& input, list<pair<string,int> > history, int max_lines = -1, string testname="")
{
	bool first_line_sent_to_output = false;
	bool switch_line_sent_to_output = false;
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: max_lines = " << max_lines << ", testname = " << testname << endl;
	
	string return_string = "";
	
	string next_line;
	
	map<string,string> var_map;
	
	int line_number = 0;
	int output_lines = 0;
	
	while(getline(input,next_line))
	{
		line_number++;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: line_number " << line_number << ": " << next_line << endl;

		if (next_line[0] == '#')
		{
			if (next_line.find("#common") != string::npos && testname != "")
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: #common line= " << next_line << endl;

				size_t pos = 7;
				while (next_line[pos] == ' ' || next_line[pos] == '\t')
					pos++;
											
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: common options = " << next_line.substr(pos) << " found at line " << line_number << endl;
				
				os << next_line.substr(pos) << endl;
			}
			else if (next_line.find("#define") != string::npos)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: #define line= " << next_line << endl;

				string def_map = next_line.substr(7);
				istringstream iss(def_map);
				string field1, field2;
				iss >> field1;
				iss >> field2;
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: #define " << field1 << " to " << field2 << endl;

				var_map[field1] = field2;
			}
			else if (next_line.find("#end") != string::npos)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: #end encountered in " << history.rbegin()->first << ", stopping #include at line " << line_number << endl;
				break;
			}
			else if (next_line.find("#include") != string::npos)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: #include line= " << next_line << endl;

				int max_lines_to_include = -1;
				size_t pos = next_line.find('{');
				if (pos != string::npos)
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: specifies max_lines " << next_line.substr(pos) << endl;

					get_number(max_lines_to_include,next_line,pos+1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: max_lines_to_include = " << max_lines_to_include << endl;
					pos  = next_line.find('}')+1;
				}
				else
					pos = 8;

				while (next_line[pos] == ' ' || next_line[pos] == '\t')
					pos++;
					
				next_line = next_line.substr(pos);
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: include_string = " << next_line << endl;
	
				
				string include_file = expand_string(next_line,var_map);
				
				ifstream is(include_file);
				if (!is)
				{
					list<pair<string,int> >::iterator hptr = history.begin();
					cout << "\nIn root file " << hptr->first << endl;
					hptr++;
					while (hptr != history.end())
					{
						cout << "in included file " << hptr->first << " at line number " << hptr->second << endl;
						hptr++;
					}
					cout << "could not find included file " << include_file << " at line " << line_number << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	list<pair<string,int> >::iterator hptr = history.begin();
	debug << "preprocess_file: in root file " << hptr->first << endl;
	hptr++;
	while (hptr != history.end())
	{
		debug << "preprocess_file: in included file " << hptr->first << " at line number " << hptr->second << endl;
		hptr++;
	}
	debug << "preprocess_file: could not find included file " << include_file << " at line " << line_number << endl;
}
					exit(0);
				}
				else
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: including " << include_file << " at line " << line_number << endl;

					list<pair<string,int> > include_history = history;
					include_history.push_back({include_file,line_number});	

					recursion_level++;
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: recursing to level " << recursion_level << endl;											
					preprocess_file(os,is,include_history,(max_lines == -2? max_lines: max_lines_to_include));
				}
			}
			else if (next_line.find("#result") != string::npos || next_line.find("#test") != string::npos)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: #result line= " << next_line << endl;
				bool found_result = (next_line.find("#result") != string::npos?true:false);
				size_t pos = (found_result?7:5);
				while (next_line[pos] == ' ' || next_line[pos] == '\t')
					pos++;
				size_t end = pos+1;
				while (next_line[end] != ' ' && next_line[end] != '\t')
					end++;


				string name = next_line.substr(pos,end-pos);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: " << (found_result?"#result ":"#test ") << name << " found at line " << line_number << endl;

				pos = end;
				if (name == testname)
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: name matches testname" << endl;
					
					if (found_result)
					{
						while (next_line[pos] == ' ' || next_line[pos] == '\t')
							pos++;
						
						/* return the expanded result filename */
						return_string = expand_string(next_line.substr(pos),var_map);
					}
					else
					{
						os << next_line.substr(next_line.find('[',end)) << endl;
					}
					
				}				
			}
			else if (next_line.find("#script") != string::npos)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: #script line= " << next_line << endl;
				size_t pos = 7;
				while (next_line[pos] == ' ' || next_line[pos] == '\t')
					pos++;
					
				next_line = next_line.substr(pos);
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: script string = " << next_line << endl;
				
				string script_file = expand_string(next_line,var_map);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: script string = " << next_line << ", script_file = " << script_file << endl;
	
				os << script_file << endl;
			}
		}
		else
		{				
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "preprocess_file: first_line_sent_to_output = " << first_line_sent_to_output;
	debug << " switch_line_sent_to_output = " << switch_line_sent_to_output << endl;
}

			if (first_line_sent_to_output)
			{
				/* If we are in an included file, and max_lines == -2, then only put out a second line if it is a cocycle */
				if (recursion_level > 0 && max_lines == -2)
				{
					if(switch_line_sent_to_output && next_line.find('C') != string::npos)
						os << next_line << endl;
				}
				else
				{
					os << next_line << endl;
				}	
			}
			else
			{
				os << next_line << endl;
				if (next_line.find('S') != string::npos)
					switch_line_sent_to_output = true;
			}
			
			if (next_line.length() && next_line[0] != ';' && next_line.substr(0,2) != "--")
			{
				++output_lines;
				
				if (recursion_level > 0 && max_lines == -2 && first_line_sent_to_output)
					break;
				else
					first_line_sent_to_output = true;
			}

			if (output_lines == max_lines)
				break;
		}
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "preprocess_file: returning " << return_string << " at level " << recursion_level << endl;
	
	if (recursion_level>0)
		recursion_level--;

	return return_string;
}

string expand_string (string input_string,map<string,string>& var_map)
{	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "expand_string: presented with input_string " << input_string << endl;

	bool expansion_complete = false;
	do
	{
		ostringstream oss;
		size_t pos=0;
		while (pos < input_string.length())
		{
			size_t next = input_string.find("$(",pos);
			if (next == string::npos)
			{
				oss << input_string.substr(pos);
				if (pos == 0)
					expansion_complete = true;
				break;
			}
			else
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "expand_string: found variable after " << pos << " at " << next << endl;

				oss << input_string.substr(pos,next-pos);
				size_t end = input_string.find(')',next);
				string variable = input_string.substr(next+2,end-next-2);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "expand_string: variable = " << variable << endl;

				if(var_map.find(variable) != var_map.end())
					oss << var_map[variable];
				pos = end+1;
			}
		}					
		input_string = oss.str();
	} while (!expansion_complete);
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "expand_string: returning " << input_string << endl;

	return input_string;
}
