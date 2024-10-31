/**************************************************
                  util.cpp

This file is intended to hold common utilities 

**************************************************/

#include <iostream>
#include <ios>
#include <string>
#include <sys/stat.h>

using namespace std;

#include <util.h>

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

void substitute_string_variable(string& str, char ch1, char ch2) 
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

namespace util {
string itos (long n)
{
	unsigned int negative = 0;
    if (n<0)
    {
		negative++;
		n *= -1;
    }

	/* count the digits in n */
	int num_digits = 0;
	long nn = n;
	do
	{
		num_digits++;
		nn /= 10;
	} while (nn != 0);

	/* size of the string is num_digits + negative, the use of the '-'
	   character removes the requirement to set it explicitly.
    */
	string res(num_digits+negative,'-');
	string::reverse_iterator ptr = res.rbegin();
	
    do
    {
		switch (n%10)
		{
	    	case 0: *ptr++ = '0'; break;
	    	case 1: *ptr++ = '1'; break;
	     	case 2: *ptr++ = '2'; break;
	     	case 3: *ptr++ = '3'; break;
	     	case 4: *ptr++ = '4'; break;
	     	case 5: *ptr++ = '5'; break;
	     	case 6: *ptr++ = '6'; break;
	     	case 7: *ptr++ = '7'; break;
	     	case 8: *ptr++ = '8'; break;
	     	case 9: *ptr++ = '9'; break;
		}
		n /= 10;
    } while(n != 0);
	
	return res;
}
}


int least_common_multiple (list<int>& l)
{
	if (l.size() == 0)
		return 0;
	
	list<int>::iterator lptr = l.begin();
	
	int lcm = *lptr;
	
	while (lptr != l.end())
	{
		if (*lptr != 1 && *lptr != lcm)
		{
			int g=gcd(lcm, *lptr);
			lcm *= *lptr;
			lcm /= g;
		}
		
		lptr++;
	}
	
	return lcm;
}

int least_common_multiple (vector<int>& v)
{
	int lcm = v[0];
	
	for (size_t i=1; i< v.size(); i++)
	{
		if (v[i] != 1 && v[i] != lcm)
		{
			int g=gcd(lcm, v[i]);
			lcm *= v[i];
			lcm /= g;
		}
	}
	
	return lcm;
}

