/***********************************************************************************************
 debug.cpp defines the debug ofstream and provides the functions to set debug booleans from the 
 command line.
 
 The debug subsystem works by having the main programme call debug_setup with a pointer to a 
 command line short option starting with -# and taking the form:
 
 -#[option:option{parameter:parameter}:option]

debug_setup returns true if [...] is present, and false otherwise, allowing default debug options 
to be set by the calling programme.

if debug options are provided, debug_setup calls set_debug_option for each 'option' or 
'option{parameter:parameter}' on the command line.  The function set_debug_option uses the #defines 
in initialize.h to determine what class-based debugging to activate.  It sets any general debug parameters 
for the activated classes and calls check_debug_option_parameters for each of their options that may have 
been specified.  For non-class-based options, set_debug_option calls set_main_debug_option, which has to be
provided by the main programme.

set_main_debug_option_parameter provides the main-specific equivalent to set_debug_option for options
relevant to the main programme.  This function also calls check_debug_option_parameters to set any debug
parameters associated with those options.

check_debug_option_parameters simply works through parameters, calling set_debug_option_parameter for each 
one.  The function set_debug_option_parameter is also "initialize-aware", setting the parameters only for 
those classes that have been initialized, based on the 'option' parameter it is given.  If this option is
main programme specific, set_debug_option_parameter calls set_main_debug_option_parameter, which also has 
to be provided by the main programme.

The function debug_help may be called by the main programme, if help is offered, via any means preferred by
the calling programme, usually using the # character in the help option.  Class based help is provided by a
call to set_debug_option with two zero parameters.  Help on the main debug options is prompted by a call to
set_main_debug_option with two zeros but it is up to that programme to determine what help is offered. 
The debug_help function also displays default debug options (as set by set_debug_option) and calls
main_display_default_options, again provided by the main programme, di display any default debug parameters 
for the programme specific options that may have been set by set_main_debug_option.
 
*********************************************************************************************/

#include <string>
#include <sstream>

#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <map>

using namespace std;

ofstream debug;

/* The preprocessor directives used here should be provided by the main programme to describe what debug
   capabilities the debug subsystem should provide.
   
   INITIALIZE_BIGINT
   INITIALIZE_BIGREAL
   INITIALIZE_MATRIX
   INITIALIZE_POLYNOMIAL
   INITIALIZE_QUATERNION
   INITIALIZE_RATIONAL
   INITIALIZE_SURD

   
*/
#include <debug-control.h>
#include <initialize.h>


int debug_control::DEBUG = debug_control::OFF; 

#include <util.h>
#include <algorithm>

//#ifdef INITIALIZE_SCALAR
//	#include <scalar.h>
//#endif


#ifdef INITIALIZE_BIGINT
	#include <bigint.h>
#endif

#ifdef INITIALIZE_BIGREAL
	#include <bigreal_control.h>
#endif


#ifdef INITIALIZE_MATRIX	
	#include <matrix.h>
#endif

#ifdef INITIALIZE_QUATERNION
	#include <quaternion.h>
#endif

#ifdef INITIALIZE_RATIONAL
	#include <rational.h>
#endif

#ifdef INITIALIZE_POLYNOMIAL
	#include <polynomial.h>
#endif

#ifdef INITIALIZE_SURD
	#include <surd_control.h>
#endif

//void set_debug_option(char* start, char* end);
void set_debug_option(string option);
void check_debug_option_parameters(size_t pos, string option);
void set_debug_option_parameter(string parameter, string option);

/* These functions have to be provided by the calling programme */
void set_main_debug_option(string option);
void set_main_debug_option_parameter(string parameter, string option);
void main_display_default_options();

/* The parameter to debug_setup is the argv[] pointer that contains the debug option
   The format of the debug option is #[option:option{parameter:parameter}:option] etc.
   and this function removes the [...] from the argument to allow the calling programme 
   to process other options independently of the debug options.
   
   The function returns true if [...] is present, and false otherwise, this allows the default 
   debug options to be set by the calling programme
*/
bool debug_setup(string argument)
{
	size_t start = argument.find('#');
	
	if (argument[++start] == '[')
	{
		/* make sure there's an end to the option definitions*/
		if (argument.find(']') == string::npos)
		{
			cout << "\nError in command line: debug options must be enclosed within []\n";
			exit(0);
		}

		do
		{
			start++; // step over '[', or ':'

			if (!isalpha(argument[start]))
			{
				cout << "\nError in command line: debug options must start with a keyword\n";
				exit(0);
			}

			size_t end = start;

			while (argument[end] != ':' && argument[end] != ']')
			{						
				if (argument[end] == '{')
				{					
					end++;
					/* move to the end of the option parameters */
					while (argument[end] != '}')
					{
						 if (argument[end] == ']')
						 {
							 cout << "\nError in command line: debug option parameters must be enclosed in {}\n";
							 exit(0);
						 }
						 else
							end++;									 
					}
					end++; // step over '}'
				}
				else
					end++;
			}

			if (end>start) // deals with the degenerate case -#[]
				set_debug_option(argument.substr(start, end-start));
			start = end;

		} while (argument[start] != ']');
		
		return true;
	}
	else
		return false;
}

void set_debug_option(string option)
{
	size_t pos = option.find('{');
	
//cout << "set_debug_option isolated option string " << option << endl;

	/* if we've been given the empty string display debug help information */
	if (option == "")
	{
#ifdef INITIALIZE_BIGINT
		cout << "\t\tbigint{san:+:-:*:/:%:rs:ls:bc:==:gt:out:in:sum:diff:num_len:gcd:all}, bitmap: no default" << endl;
#endif
#ifdef INITIALIZE_BIGREAL
		cout << "\t\tbigreal{+:-:*:/:==:gt:out:in:sum:diff:num_len:carry:all}, bitmap: no default" << endl;
#endif
#ifdef INITIALIZE_MATRIX
		cout << "\t\tmatrix{gen:det:inv:*:all}, bitmap: default general" << endl;
#endif
#ifdef INITIALIZE_POLYNOMIAL
		cout << "\t\tpoly{no-gen:san:+:*:/:in:gcd:all}, bitmap: general debug included unless no-gen specified" << endl;
#endif
#ifdef INITIALIZE_QUATERNION
		cout << "\t\tquaternion, boolean" << endl;
#endif
#ifdef INITIALIZE_RATIONAL
		cout << "\t\trational, boolean" << endl;
#endif
#ifdef INITIALIZE_SURD
		cout << "\t\tsurd{+:-:*:==:gt:in:newton:all}, bitmap: no default" << endl;
#endif
		return;
	}

	if (option.find("bigint") != string::npos)	
	{
#ifdef INITIALIZE_BIGINT
		if (pos != string::npos)
		{
			/* check for any parameters */
			check_debug_option_parameters(pos,option);
		}
#endif
	}
	else if (option.find("bigreal") != string::npos)	
	{
#ifdef INITIALIZE_BIGREAL
		if (pos != string::npos)
		{
			/* check for any parameters */
			check_debug_option_parameters(pos,option);
		}
#endif
	}
	else if (option.find("matrix") != string::npos)	
	{
#ifdef INITIALIZE_MATRIX
		matrix_control::DEBUG |= matrix_control::general; // set default debug option
		debug << "debug: setting debug option matrix_control::general\n";		
		if (pos != string::npos)
		{
			/* check for any parameters */
			check_debug_option_parameters(pos,option);
		}
#endif
	}
	else if (option.find("poly") != string::npos)	
	{
#ifdef INITIALIZE_POLYNOMIAL
		polynomial_control::DEBUG |= polynomial_control::general; // set default debug option
		debug << "debug: setting debug option polynomial_control::general\n";		

		if (pos != string::npos)
		{
			/* check for any parameters */
			check_debug_option_parameters(pos,option);
		}
#endif
	}
	else if (option.find("rational") != string::npos)	
	{
#ifdef INITIALIZE_RATIONAL
		rational_control::DEBUG = true;
		debug << "debug: setting debug option rational_control::DEBUG\n";		

#endif
	}
	else if (option.find("quaternion") != string::npos)	
	{
#ifdef INITIALIZE_QUATERNION
		quaternion_control::DEBUG = true;
		debug << "debug: setting debug option quaternion_control::DEBUG\n";		
#endif
	}
	else if (option.find("surd") != string::npos)	
	{
#ifdef INITIALIZE_SURD
		if (pos != string::npos)
		{
			/* check for any parameters */
			check_debug_option_parameters(pos,option);
		}
#endif
	}
	else
		set_main_debug_option(option); // function has to be provided by calling programme


	debug.flush();
}

/* work through the parameters from pos and call set_debug_parameter_option for each one 
*/
void check_debug_option_parameters(size_t pos, string option)
{
	string base_option = option.substr(0,pos);
//cout << "check_debug_option_parameters base option = " << base_option << endl;

	bool end_of_parameters = false;	
	size_t start = ++pos; // move past opening {
	do
	{
		if (option[pos] == ':' || option[pos] == '}')
		{
		
			if (option[pos] == '}')
				end_of_parameters = true;
	
//cout << "check_debug_option_parameters parameter " << option.substr(start,pos-start) << endl;
				
			set_debug_option_parameter(option.substr(start,pos-start), base_option);

			if (end_of_parameters)
				break;
			else
				start = ++pos;
		}
		else
			pos++;

	} while(true); // we break out of this loop
}

void set_debug_option_parameter(string parameter, string option)
{
	if (option == "bigint")
	{
#ifdef INITIALIZE_BIGINT
		if (parameter == "san")
		{
			bigint_control::DEBUG |= bigint_control::sanitize;
			debug << "debug: setting debug option bigint_control::sanitize\n";		
		}
		else if (parameter == "+")
		{
			bigint_control::DEBUG |= bigint_control::add;
			debug << "debug: setting debug option bigint_control::add\n";		
		}
		else if (parameter == "-")
		{
			bigint_control::DEBUG |= bigint_control::subtract;
			debug << "debug: setting debug option bigint_control::subtract\n";		
		}
		else if (parameter == "*")
		{
			bigint_control::DEBUG |= bigint_control::multiply;
			debug << "debug: setting debug option bigint_control::multiply\n";		
		}
		else if (parameter == "/")
		{
			bigint_control::DEBUG |= bigint_control::divide;
			debug << "debug: setting debug option bigint_control::divide\n";		
		}
		else if (parameter == "%")
		{
			bigint_control::DEBUG |= bigint_control::remainder;
			debug << "debug: setting debug option bigint_control::remainder\n";		
		}
		else if (parameter == "rs")
		{
			bigint_control::DEBUG |= bigint_control::r_shift;
			debug << "debug: setting debug option bigint_control::r_shift\n";		
		}
		else if (parameter == "ls")
		{
			bigint_control::DEBUG |= bigint_control::l_shift;
			debug << "debug: setting debug option bigint_control::l_shift\n";		
		}
		else if (parameter == "bc")
		{
			bigint_control::DEBUG |= bigint_control::bool_conv;
			debug << "debug: setting debug option bigint_control::bool_conv\n";		
		}
		else if (parameter == "==")
		{
			bigint_control::DEBUG |= bigint_control::equal;
			debug << "debug: setting debug option bigint_control::equal\n";		
		}
		else if (parameter == "gt")
		{
			bigint_control::DEBUG |= bigint_control::greater;
			debug << "debug: setting debug option bigint_control::greater\n";		
		}
		else if (parameter == "out")
		{
			bigint_control::DEBUG |= bigint_control::output;
			debug << "debug: setting debug option bigint_control::output\n";		
		}
		else if (parameter == "in")
		{
			bigint_control::DEBUG |= bigint_control::input;
			debug << "debug: setting debug option bigint_control::input\n";		
		}
		else if (parameter == "sum")
		{
			bigint_control::DEBUG |= bigint_control::sum;
			debug << "debug: setting debug option bigint_control::sum\n";		
		}
		else if (parameter == "diff")
		{
			bigint_control::DEBUG |= bigint_control::diff;
			debug << "debug: setting debug option bigint_control::diff\n";		
		}
		else if (parameter == "num_len")
		{
			bigint_control::DEBUG |= bigint_control::num_len;
			debug << "debug: setting debug option bigint_control::num_len\n";		
		}
		else if (parameter == "gcd")
		{
			bigint_control::DEBUG |= bigint_control::gcd;
			debug << "debug: setting debug option bigint_control::gcd\n";		
		}
		else if (parameter == "all")
		{
			bigint_control::DEBUG = 0xFFFFFFFF;
			debug << "debug: setting debug option bigint_control::all\n";		
		}
#endif
	}
	else if (option == "bigreal")
	{
#ifdef INITIALIZE_BIGREAL
//		if (parameter == "san")
//		{
//			bigreal_control::DEBUG |= bigreal_control::sanitize;
//			debug << "debug: setting debug option bigreal_control::sanitize\n";		
//		}
		if (parameter == "+")
		{
			bigreal_control::DEBUG |= bigreal_control::add;
			debug << "debug: setting debug option bigreal_control::add\n";		
		}
		else if (parameter == "-")
		{
			bigreal_control::DEBUG |= bigreal_control::subtract;
			debug << "debug: setting debug option bigreal_control::subtract\n";		
		}
		else if (parameter == "*")
		{
			bigreal_control::DEBUG |= bigreal_control::multiply;
			debug << "debug: setting debug option bigreal_control::multiply\n";		
		}
		else if (parameter == "/")
		{
			bigreal_control::DEBUG |= bigreal_control::divide;
			debug << "debug: setting debug option bigreal_control::divide\n";		
		}
		else if (parameter == "==")
		{
			bigreal_control::DEBUG |= bigreal_control::equal;
			debug << "debug: setting debug option bigreal_control::equal\n";		
		}
		else if (parameter == "gt")
		{
			bigreal_control::DEBUG |= bigreal_control::greater;
			debug << "debug: setting debug option bigreal_control::greater\n";		
		}
		else if (parameter == "out")
		{
			bigreal_control::DEBUG |= bigreal_control::output;
			debug << "debug: setting debug option bigreal_control::output\n";		
		}
		else if (parameter == "in")
		{
			bigreal_control::DEBUG |= bigreal_control::input;
			debug << "debug: setting debug option bigreal_control::input\n";		
		}
		else if (parameter == "sum")
		{
			bigreal_control::DEBUG |= bigreal_control::sum;
			debug << "debug: setting debug option bigreal_control::sum\n";		
		}
		else if (parameter == "diff")
		{
			bigreal_control::DEBUG |= bigreal_control::diff;
			debug << "debug: setting debug option bigreal_control::diff\n";		
		}
		else if (parameter == "num_len")
		{
			bigreal_control::DEBUG |= bigreal_control::num_len;
			debug << "debug: setting debug option bigreal_control::num_len\n";		
		}
		else if (parameter == "carry")
		{
			bigreal_control::DEBUG |= bigreal_control::carry;
			debug << "debug: setting debug option bigreal_control::carry\n";		
		}
		else if (parameter == "all")
		{
			bigreal_control::DEBUG = 0xFFFFFFFF;
			debug << "debug: setting debug option bigreal_control::all\n";		
		}
#endif
	}
	else if (option == "matrix")
	{
#ifdef INITIALIZE_MATRIX	
		if (parameter == "gen")
		{
			matrix_control::DEBUG |= matrix_control::general;
			debug << "debug: setting debug option matrix_control::general\n";		
		}
		if (parameter == "det")
		{
			matrix_control::DEBUG |= matrix_control::immanant;
			debug << "debug: setting debug option matrix_control::immanant\n";		
		}
		if (parameter == "ech")
		{
			matrix_control::DEBUG |= matrix_control::echelon;
			debug << "debug: setting debug option matrix_control::echelon\n";		
		}
		if (parameter == "*")
		{
			matrix_control::DEBUG |= matrix_control::multiply;
			debug << "debug: setting debug option matrix_control::multiply\n";		
		}
		if (parameter == "inv")
		{
			matrix_control::DEBUG |= matrix_control::inverse;
			debug << "debug: setting debug option matrix_control::inverse\n";		
		}
		else if (parameter == "all")
		{
			matrix_control::DEBUG = 0xFFFFFFFF;
			debug << "debug: setting debug option matrix_control::all\n";		
		}
#endif
	}
	else if (option == "poly")
	{
#ifdef INITIALIZE_POLYNOMIAL
		if (parameter == "no-gen")
		{
			polynomial_control::DEBUG >>= 1;
			polynomial_control::DEBUG <<= 1;
			debug << "debug: clearing debug option polynomial_control::general\n";		
		}
		else if (parameter == "san")
		{
			polynomial_control::DEBUG |= polynomial_control::sanitize;
			debug << "debug: setting debug option polynomial_control::sanitize\n";		
		}
		else if (parameter == "+")
		{
			polynomial_control::DEBUG |= polynomial_control::add;
			debug << "debug: setting debug option polynomial_control::add\n";		
		}
		else if (parameter == "*")
		{
			polynomial_control::DEBUG |= polynomial_control::multiply;
			debug << "debug: setting debug option polynomial_control::multiply\n";		
		}
		else if (parameter == "/")
		{
			polynomial_control::DEBUG |= polynomial_control::divide;
			debug << "debug: setting debug option polynomial_control::divide\n";		
		}
		else if (parameter == "gcd")
		{
			polynomial_control::DEBUG |= polynomial_control::gcd;
			debug << "debug: setting debug option polynomial_control::gcd\n";		
		}
		else if (parameter == "in")
		{
			polynomial_control::DEBUG |= polynomial_control::input;
			debug << "debug: setting debug option polynomial_control::input\n";		
		}
		else if (parameter == "all")
		{
			polynomial_control::DEBUG = 0xFFFFFFFF;
			debug << "debug: setting debug option polynomial_control::all\n";		
		}
#endif
	}	
	if (option == "surd")
	{
#ifdef INITIALIZE_SURD
		if (parameter == "+")
		{
			surd_control::DEBUG |= surd_control::add;
			debug << "debug: setting debug option surd_control::add\n";		
		}
		else if (parameter == "-")
		{
			surd_control::DEBUG |= surd_control::subtract;
			debug << "debug: setting debug option surd_control::subtract\n";		
		}
		else if (parameter == "*")
		{
			surd_control::DEBUG |= surd_control::times;
			debug << "debug: setting debug option surd_control::times\n";		
		}
		else if (parameter == "==")
		{
			surd_control::DEBUG |= surd_control::equal;
			debug << "debug: setting debug option surd_control::equal\n";		
		}
		else if (parameter == "gt")
		{
			surd_control::DEBUG |= surd_control::greater;
			debug << "debug: setting debug option surd_control::greater\n";		
		}
		else if (parameter == "in")
		{
			surd_control::DEBUG |= surd_control::input;
			debug << "debug: setting debug option surd_control::input\n";		
		}
		else if (parameter == "newton")
		{
			surd_control::DEBUG |= surd_control::newton;
			debug << "debug: setting debug option surd_control::newton\n";		
		}
		else if (parameter == "all")
		{
			surd_control::DEBUG = 0xFFFFFFFF;
			debug << "debug: setting debug option surd_control::all\n";		
		}
#endif
	}
	else
		set_main_debug_option_parameter(parameter, option); // has to be provided by the calling programme
}

void debug_help ()
{
	cout << "\n\tdebug options may be specified using the syntax -#[option:option{parameter:parameter}:option]" << endl;
	cout << "\tsupported options and parameters:" << endl;
	set_debug_option("");
	set_main_debug_option("");
	cout << "\tdefault options set by omitting [...]" << endl;
	cout << "\tdefault options:" << endl;
	main_display_default_options();
	cout << endl;
	exit(0);
}
