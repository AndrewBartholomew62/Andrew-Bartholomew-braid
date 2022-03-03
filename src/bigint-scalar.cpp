/*********************************************************************** 
			bigint A. Bartholomew 2nd October 2005

    bigint is an arbitrary precision integer class with an interface designed 
	to meet the requirements of the braid programme.  It is based on earlier 
	int128 and int64 code but uses the division algorithm in Knuth The Art 
	of Computer Programing Vol 2 (third edition), page 272 (Algorithm D).  

	The class uses radix 2^16 storing the numerals in reverse order in a vector,
	so that the ith place stores the numeral for the ith power of the radix.  The 
	sign is stored separately.

***********************************************************************/

#include <string>
#include <sstream>

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>

using namespace std;

extern ofstream debug;

#include <util.h>
#include <scalar.h>

/* This is the code itself, separated into a separate file to allow for 
   stand-alone bigints as well as this scalar variant
*/
#include "bigint-cpp.cpp"
