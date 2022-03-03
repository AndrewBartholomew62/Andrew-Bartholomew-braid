/**************************************************
                  util.cpp

This file is intended to hold common utilities 

**************************************************/

#include <iostream>
#include <ios>
#include <string>
#include <sys/stat.h>

using namespace std;

int num_len (long n)
{
    int i = 0;

    if (n<0)
    {
        i++;
        n *= -1;
    }

    do
    {
        i++;
        n /= 10;
    } while(n != 0);

    return i;
}

char* c_string(const string& s)
{
	char* p = new char [s.length()+1];
	p[s.copy(p,string::npos)] = 0; // add terminator
	return p;
}

void print_state(istream& is)
{
cout << "\nprint_state:" << endl;
	ios_base::iostate s = is.rdstate();
cout << "\niostate: " << s << endl;
	if (s&ios_base::badbit) cout << "\nbad bit set";
	if (s&ios_base::failbit) cout << "\nfail bit set";
	if (s&ios_base::eofbit) cout << "\neof bit set";
	if (s&ios_base::eofbit) cout << "\ngood bit set";
	cout << endl;
}

//This function is used only for debugging .h files
void breakpoint() {}


// Function: fileExists
/**
    Check if a file exists
@param[in] filename - the name of the file to check

@return    true if the file exists, else false

*/
bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}

void stringSubstitute(string& str, char ch1, char ch2) 
{
  for (unsigned int i = 0; i < str.length(); ++i) 
  {
    if (str[i] == ch1)
      str[i] = ch2;
  }
}

char next_variable (string variables)
{
	char ch;
	bool found = false;
	for (int i=0; i<26; i++)
	{
		ch='a'+i;
		if (variables.find(ch) == string::npos)
		{
			found = true;
			break;
		}
	}
	
	if (!found)
	{
		for (int i=0; i<26; i++)
		{
			ch='A'+i;
			if (variables.find(ch) == string::npos)
			{
				found = true;
				break;
			}
		}
	}
	
	if (found)
		return ch;
	else
	{
		cout << "\nError! next_variable cannot find another variable for " << variables << endl;
		exit(0);
	}
}
