/***********************************************************************
			   preprocessor-main.cpp

		A. Barthomew November 2025
		
***********************************************************************/
#include <fstream>
#include <iostream>
//#include <cstdlib>
#include <cstring>
//#include <iomanip>
//#include <algorithm>
#include <list>

using namespace std;

#include <debug-control.h>
#include <util.h>

extern ofstream    debug;

string preprocess_file (ostream& os, ifstream& input, list<pair<string,int> > history, int max_lines = -1, string testname="");

int main(int argc, char* argv[])
{
	debug_control::DEBUG = debug_control::SUMMARY;
	
	string filename;

	if (argc == 1)
	{
		cout << "Usage: preprocess <filename1>" << endl;
	}
	else
	{
		for (int i=1; i< argc; i++)
		{
			if ( *(argv[i]) == '-')
			{		
				if (*(argv[i]+1) == '-')
				{
					string argument = argv[i];
//					if (argument.find("suite") != string::npos)
//						running_test_suite = true;
//					else if (argument.find("summary") != string::npos)
//						summary_test = true;
				}
				else
				{
					char* cptr = strchr(argv[i], '#');
					if (cptr)
					{
						debug_control::DEBUG = debug_control::SUMMARY;
						debug.open("preprocess.dbg");
		    			if (!debug)
		    			{
		        			cout << "\nError opening debug file\n";
		        			exit(0);
		    			}
						else
							debug << "Debug information from preprocess\n\n";				
					}		
				}
			}
			else if (filename.length() == 0)
			{
				filename = argv[i];
			}
/*			else if (testname.length() == 0)
			{
				testname = argv[i];
			}
*/		}
	}	

	string root(filename);
	ifstream input(root);
	if (!input)
	{
		cout << "Error! Could not open input file " << filename << endl;
		exit(0);
	}

	list<pair<string,int> > history;
	history.push_back({filename,0});
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "main: preprocessing_file " << filename << endl;
	
	string result_file = preprocess_file(cout, input, history);
	cout << "\npreprocess_file returns " << result_file << endl;
	
	if (debug_control::DEBUG)
		debug.close();
	return 0;
}


/* **************  Functions required for setting up programme specific debug ***************/
void set_main_debug_option(char* start, char* end) {}
void set_main_debug_option_parameter(char* pptr, string option) {}
void main_display_default_options() 
{
	cout << "\t\tnone" << endl;
}




