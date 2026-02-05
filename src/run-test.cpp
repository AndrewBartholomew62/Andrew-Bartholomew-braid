/*********************************************************************************************************
			   run-test.cpp		A. Barthomew November 2025

run-test provides an automated testing facility for the braid programme.  

Tests are described in test scripts that use the processor directive 

#test <test-name> [braid-options]
#result <test-name> <result-file>

to identify one or more tests that may be run with the same set of switches and input codes.  The input
data may appear in the test script or be #included by a preprocesor directive, the only restriction is
that all the tests in a script use the same input data.

run-test can run an individual test contained in a test script or all of them.  It scans the test script
for #test directives and either checks that the specified <test-name> exists or builds a list of al the
<test-name> strings in the file.

For each test to be run, run-test sequentially uses the preprocessor to produce a command_line_input_file 
called _tmpfile for a given test.  The _tmpfile produced by the preprocessor includes all #common programe
options and any that appear in the specified #test directive and returns the file containing the expected 
result, as given by the #result directive corresponding to the <test-name.

The _tmpfile is then passed to the braid proramme via a system call and the braid.out file compared to the
result file using diff.  The failure of a test does not halt operation, run-test just proceeds to the next 
test.

Using the command line option --suite, run-test can run a suite of test scripts listed in a file <suite-name>
This file is list of test scripts that are desired to be run together prefixed with the #script preprocessor 
directive.  The file may also contain comments and use the preprocessor variable definitions to construct 
filenames.  The programme collates the output from each test script into a single output file <suite-name>.out

		
*********************************************************************************************************/
#include <fstream>
#include <iostream>
#include <cstring>
#include <list>

using namespace std;

#include <debug-control.h>
#include <util.h>

extern ofstream    debug;

string preprocess_file (ostream& os, ifstream& input, list<pair<string,int> > history, int max_lines = -1, string testname="");

void run_test(ofstream& output, string filename, bool summary_test, string testname="");

bool LIST_TESTS = false;

int main(int argc, char* argv[])
{
	debug_control::DEBUG = debug_control::OFF;

	bool running_test_suite = false; 
	bool summary_test = false; // a summary test #includes just the first line of each file
	
	string filename;
	string testname;
	
	if (argc == 1)
	{
		cout << "Usage: 'run-test <filename>' runs every test in <filename>" << endl;
		cout << "       'run-test --suite <suite-name>' the suite of test files listed in <suite-name>" << endl;
		cout << "       'run-test --suite --summary <suite-name>' one item from each included file in the suite of tests" << endl;
		cout << "       'run-test <filename> <testname>' runs <testname> in <filename>" << endl;
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
					if (argument.find("suite") != string::npos)
						running_test_suite = true;
					else if (argument.find("summary") != string::npos)
						summary_test = true;
				}
				else
				{
					char* cptr = strchr(argv[i], '#');
					if (cptr)
					{
						debug_control::DEBUG = debug_control::SUMMARY;
						debug.open("run-test.dbg");
		    			if (!debug)
		    			{
		        			cout << "\nError opening debug file\n";
		        			exit(0);
		    			}
						else
							debug << "Debug information from run-test\n\n";				
					}
		
					cptr = strchr(argv[i], 'l');
					if (cptr)
						LIST_TESTS = true;
				}
			}
			else if (filename.length() == 0)
			{
				filename = argv[i];
			}
			else if (testname.length() == 0)
			{
				testname = argv[i];
			}
		}
	}	
	
	/* preprocess produces a _suitefile that has expanded any variables
	   used in #script directives.  We then call run_test with each script name
	*/
	if (running_test_suite)
	{
		ofstream 	output;		
		output.open(filename+".out");
		if (!output)
		{
			cout << "\nError opening output file " << filename+".out" << endl;
			exit(0);
		}
		
		ifstream input(filename);
		if (!input)
		{
			cout << "Error opening " << filename << endl;
			exit(0);
		}
		
		ofstream os("_suitefile");
		if (!os)
		{
			cout << "Error!  run_test could not open _suitefile for output" << endl;
			exit(0);
		}

		list<pair<string,int> > history;
		history.push_back({filename,0});
		preprocess_file(os,input,history); 
		
		input.close();
		input.open("_suitefile");
		if (!input)
		{
			cout << "Error opening _suitefile for input" << endl;
			exit(0);
		}
		
		string next_line;		
		while (getline(input,next_line))
		{
			if (next_line.length() && next_line[0] != ';') // allow comments
			{
				istringstream iss(next_line);
				string test_script;
				iss >> test_script;
				
				cout << "Run test script " << test_script << endl;
				
				run_test(output, test_script,summary_test);				
			}
		}

	}
	else
	{
		ofstream 	output;
		output.open("run-test.out");
		if (!output)
		{
			cout << "\nError opening output file\n";
			exit(0);
		}
				
		if (testname.length() != 0)
		{
			string testname(argv[2]);
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "run-test main: running test " << testname << " in " << filename << endl;
		
			run_test(output, filename,summary_test,testname);
		}
		else
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "run-test main: running all tests in " << filename << endl;
		
			run_test(output, filename,summary_test);
		}
		
		if (debug_control::DEBUG)
			debug.close();
			
		output.close();
	}
	
	return 0;
}

void run_test(ofstream& output, string filename, bool summary_test, string testname)
{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "run-test: presented with test " << testname << " in " << filename << endl;

	list<string> test_list;

	ifstream input (filename);
	
	if  (!input)
	{
		cout << "Error!  run_test could not open " << filename << endl;
		output << "Error!  run_test could not open " << filename << endl;
		return;
	}
	
	
	if (!LIST_TESTS)
	{
		if (testname == "")
		{
			cout << "running all tests in " << filename << endl;
			output << "running all tests in " << filename << endl;
		}
		else
		{
			cout << "running test " << testname << " in " << filename << endl;
			output << "running test " << testname << " in " << filename << endl;
		}
	}


	string next_line;
	while (getline(input,next_line))
	{		
		if (next_line[0] == '#' && next_line.find("#test") != string::npos)
		{
			size_t pos = 5;
			while (next_line[pos] == ' ' || next_line[pos] == '\t')
				pos++;
			size_t end = pos+1;
			while (next_line[end] != ' ' && next_line[end] != '\t')
				end++;

			if (testname == "" || testname == next_line.substr(pos,end-pos))
				test_list.push_back(next_line.substr(pos,end-pos));
		}	
	}
	

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "run-test: test_list size = " << test_list.size() << endl;
	
	if (LIST_TESTS)
	{
		list<string>::iterator tptr = test_list.begin();
		while (tptr != test_list.end())
		{
			cout << *tptr << endl;
			tptr++;
		}
		
		return;
	}
	
	if (test_list.size() == 0)
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "run-test: no matching test in " << filename << endl;
	
		if (testname == "")
		{
			cout << "no tests found in " << filename << endl;
			output << "no tests found in " << filename << endl;
		}
		else
		{
			cout << "test not found in " << filename << endl;
			output << "test not found in " << filename << endl;
		}
		
	}
	
	list<string>::iterator lptr = test_list.begin();
	while (lptr != test_list.end())
	{
		cout << "    test " << *lptr;
		output << "    test " << *lptr;

		list<pair<string,int> > history;
		history.push_back({filename,0});
		
		ofstream os("_tmpfile");
		if (!os)
		{
			cout << "Error!  run_test could not open _tmpfile for output" << endl;
			output << "Error!  run_test could not open _tmpfile for output" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "run-test: running preprocess_file for testname " << *lptr << endl;

		input.clear();
		input.seekg(0);

		string result_file = preprocess_file(os,input,history,(summary_test? -2:-1),*lptr); // max_lines = -1; i.e. all the file -2 summary tests

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "run-test: result_file = " << result_file << endl;

		system("rm braid.out");
		system("braid _tmpfile");
		
		ifstream rf(result_file);
		if (!rf)
		{
			cout << "\ncannot open result file " << result_file << endl;
			output << "\ncannot open result file " << result_file << endl;
		}
		else
		{		
			if (summary_test)
			{
				string result_line;
				getline(rf,result_line);
				
				ifstream bf("braid.out");
				string braid_line;
				getline(bf,braid_line);
				
				if (result_line != braid_line)
				{
					cout << "...failed" << endl;
					output << "...failed" << endl;
				}
				else
				{
					cout << "...OK" << endl;
					output << "...OK" << endl;
				}
			}
			else
			{
				rf.close();

				string command = "diff braid.out " + result_file + " > _diff_output";

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "run-test: bash command string = " << command << endl;

				system(command.c_str());
				ifstream diff("_diff_output");
				
				string next_line;
				if (getline(diff,next_line)) 
				{
					cout << "...failed" << endl;
					output << "...failed" << endl;
					output << "        " << next_line << endl;
					while (getline(diff,next_line)) 
						output << "        " << next_line << endl;				
				}
				else
				{
					cout << "...OK" << endl;
					output << "...OK" << endl;
				}
			}
		}

		lptr++;
	}
	
}

/* **************  Functions required for setting up programme specific debug ***************/
void set_main_debug_option(char* start, char* end) {}
void set_main_debug_option_parameter(char* pptr, string option) {}
void main_display_default_options() 
{
	cout << "\t\tnone" << endl;
}




