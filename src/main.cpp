/********************************************************************************************************************************************************
The braid programme was originally written in 2001 to calculate the Burau 
matrix representation of a virtual braid word.

Since then, many other tools have been added accessed through command line 
options, as detailed in the revision history below.

User instructions and a programme overview are available in the user
documentation, available at www.layer8.co.uk/maths/braids 
           
   Version 1:     Coding started 25th November 2001
   Version 2:     Added quarternion coefficients to polynommials, December 2001
   Version 2.1:   Re-wrote code to suport polynomial declarations rather than
    			  always manipulating pointers, this made the code a little
    			  more elegant, though it necessitated contructors and the
    			  assignment operator.  Coded January, 2002.
   Version 2.2:   Added input files with word algebra, February 2002.
   Version 3:     Added determinant and adjoint calculation March/April 2002
   Version 4:     Added Dynnikov test April 2002
   Version 4.1:   Support for Dynnikov test of virtual braid
   Version 5:     Added gcd calculation for Alexander polynomial
   Version 6:     Added calculation of Dowker code
   Version 6.1:   Added switch specification in input file June 2002
   Version 7:     Added Sawolleck polynomial calaulation, started August 2002
   Version 8:     Added Vogel algorithm, started October 2002, completed April 2003
   Version 9:	  Added Gauss code and immersion code calculation, alternative 
   				  Delta_1 computation, polynomial invariants from immersion codes
				  and normalizing quaternionic polynomial invariants, completed November 2003
   Version 10:    Big clean up of the code to improve the C++; introduced templates 
   			  	  removed the control.h concept used in previous versions.  
				  Arbitrary precision integer arithmetic added
				  Support for matrix switches added.
				  This version dropped explicit DOS support.				
				  Completed November 2005
   Version 11:    Added Weyl algebra switches for polynomial invariants.  Made optional the 
                  use of the t variable with quaternionic, matrix and Weyl switches.  Made the 
				  calculation of Delta_1 the default, even when Delta_0 is non-zero (with an 
				  option to override).  Removed the non-standard time class used previously and 
				  replaced it with ctime. Completed may 2006
	Version 11.1: Added debug subsystem, re-wrote polynomial and bigint input routines.
	Version 11.2: Added test for A=D and B=C.
	Version 12.0: Added long knot support
	Version 12.1: Bug fixes
	Version 12.2: Added concatenation product for long knots.
	Version 12.3: Bug fixes and added error some input string error checking.
	Version 13.0: Added the HOMFLY polynomial and the (erroneous) Nicholson polynomial for flat virtual knots. 
	              it also changed the CLI syntax for options, introducing the -- long option and - short option
				  syntax.
	Version 14.0: Added the finite biquandle fixed-point invariant for virtual and welded braids, the rack-polynomial
	              invariant.  Added switch titles and introduced the notion of braid qualifiers.
	Version 15.0: Added the bracket polynomial for knotoids March 2010, extended dowker code to accommodate
	              labelled immersion codes as input April 2010.  The latter development was done in aid of the
	              development of the draw programme.
	Version 16.0: Added support of labelled peer codes for the bracket polynomial (2010) and 
	              the Vogel algorithm (January 2011).  Added the silent option for batch processing.  Removed
	              the Nicholson polynomial.
	Version 16.1: Added the flip-braid option to turn over a braid before calculating fixed point invariants. (September 2011)
	Version 16.2: Added {flip} braid qualifier as an alternative to flip-braid option. (November 2012)
	Version 17.0: Added Kauffman bracket and Jones Polynomial to the function bracket_polynomial, added support
	              for Kauffman's affine index polynomial invariant (January 2013)
	Version 17.1: Extended peer code evaluation to remove any Reidemeister II moves unless explicitly told otherwise.
	Version 18.0: Added doodle fixed-point invariants (January 2015)
	Version 18.1: Added invert-braid, line-reflect-braid and plane-reflect-braid options and {invert, line-reflect, 
	              plane-reflect} braid qualifiers (March 2015)
	Version 19.0: Added satellite knot calculation for r-parallel cables of knots (March 2017).
	Version 19.1: Modified write_gauss_code and read_gauss_code to handle doodle Gauss codes. (July 2017)
	Version 20.0: Added support for commutative automorphism switches (February 2018)
	Version 20.1: Added Kamada double covering calculation for braids (April 2018)
	Version 20.2: Added option for the Vogel algorithm to report the height of diagram only (February 2019)
	Version 21.0: Added support for links and multi-knotoids in bracket_polynomial, added gauss code to peer code conversion for classical links (November 2019)
	Version 21.1: Modified the bracket polynomial so that the kauffman-bracket and jones-polynomial options evaluated the normalized
	              bracket polynomial and Jones polynomial for knotoids and multi-knotoids as well as knots and links.  The knotoid-bracket
				  option calculates Turaev's extended bracket polynomial for knotoids and multi-knotoids as before.
				  Added TeX output format option for bracket polynomials (April 2020)
	Version 22.0: Added support for the arrow polynomial for classical and virtual knots, links, knotoids and multi-knotoids 
	              Added no-expand-bracket and no-normalize-bracket options 
				  Added support for "An Alexander type invariant for doodles" (May 2020)
				  Added mapped polynomials (6/6/20) for improved arrow polynomial presentation and to prepare for the parity bracket polynomial
				  Added support for knotoids to affine_index_polynomial
				  Added support for knotoids to the Gauss code task
				  Added the lpgd and ulpgd options to the Gauss code task
				  Added the parity bracket polynomial for classical or virtual knots or knotoids
				  Added the parity arrow polynomial for classical or virtual knots, and the relaxed parity arrow polynomial  for  classical and virtual knotoids
	Version 23.0  Updated the syntax for knot-type knotoids and long knots to be "K:" and "L:"
	              Extended the ability to convert Gauss codes to labelled peer codes based on Jeremy Green's algorithm for drawing virtual knots from Gauss codes.
	              Added support for reading Gauss codes of the form sequence{(O|U)<crossing-num><crossing-sign>, e.g. O1-O2+U1-O3+O4+U2+U3+U4+ (October 2021)
	Version 23.1  Code overhaul to provide common code base for braid, draw and vlist programmes (November 2021)
	Version 24.0  Added the opgc and uopgc options to the Gauss code task, added OU format for Gauss code output (December 2021)	
	Version 25.0  Added planar diagram support as an alternative input format, added supprt for Gauss code and planar diagram
	              input to the knotoid bracket, parity bracket and parity arrow polynomial tasks (March 2022)
	Version 25.1  Updated the handling of input files containing switches specifying the t-variable.  Updated the behaviour of "satellite" to make it an option 
                  rather than a task.  This allowed the cabling of entries in an input file prior to carrying out the specified task.  Corrected an error in the 
                  Vogel algorithm (lines 412 and 434) that caused segmentation error in some cases by not using abs.
                  Updated the sawollek implementation of matrix P to be more clearly aligned with Sawollek's paper. (September 2022)
                  Added draft mock Alexander option for knotoids (November 2022)
                  Added the ability to read Dowker-Thistlethwaite codes for prime knots and convert them to labelled peer codes (December 2022).
    Version 25.2  Added support for multi-linkoids to the peer code and Gauss code input and output functions and modified gauss_to_peer_code to convert
                  multi-knotoid Gauss codes to peer codes.(April 2023)
    Version 26.0  Added Hamiltonian task to find Hamiltonian circuits in a diagram of a classical or flat knot or link.  Removed the start menu. (October 2023)
    Version 27.0  Added prime task to test whether (the shadow of) a diagram is prime; i.e 3-connected. Cleaned up the handling of multi-linkoids (October 2023)
    Version 28.0  Added turning number calculation(August 2024), changed RACK_POLYNOMIAL task to accept peer codes in addition to braid words.  Cleaned up the
                  validation of braid input and harmonized the checking of braid qualifiers across the various tasks that use them.  
                  Implemented the Q-polynomial invariant for single component doodles(October 2024)
	
********************************************************************************************************************************************************/
using namespace std;
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <stdio.h>
#include <valarray>
#include <ctime>
#include <csignal>
#include <list>
#include <iomanip>

extern ofstream    debug;
ofstream 	output;
ifstream    input;


#include <util.h>
#include <quaternion-scalar.h>
#include <polynomial.h>
#include <matrix.h>
#include <input-control.h>
#include <braid.h>
#include <gauss-orientation.h>

/******************** Control Variables **********************/
string braid_control::version  = "28.0";

int braid_control::gcd_initializing;

bool braid_control::AFFINE_INDEX = false;
bool braid_control::ALEXANDER = false;
bool braid_control::ALWAYS_CALCULATE_CODIMENSION_2_DELTA_1 = false;
bool braid_control::ALWAYS_CALCULATE_DELTA_1 = true;
bool braid_control::ARROW_POLYNOMIAL = false;
bool braid_control::BIGELOW_KNOT_SEARCH = false;
bool braid_control::BURAU = false;
bool braid_control::BYPASS_FUNDAMENTAL_EQUATION_CHECK = false;
bool braid_control::CALCULATE_DELTA_0_ONLY = false;
bool braid_control::CALCULATE_MOD_P = false;
bool braid_control::CLASSICAL_ONLY = false;
bool braid_control::COMMUTATIVE_AUTOMORPHISM = false;
bool braid_control::COMPLEX_STUDY_DELTA_1 = false;
bool braid_control::CUSTOM_WEYL = false;
bool braid_control::DELTA_1_UNIT_CHECK = true;
bool braid_control::DEVELOPMENT_MODE = false;
bool braid_control::DISPLAY_DELTA_1_ONLY = false;
bool braid_control::DOODLE_ALEXANDER = false;
bool braid_control::DOODLE_Q_POLYNOMIAL = false;
bool braid_control::DOWKER_CODE = false;
bool braid_control::DYNNIKOV_TEST = false;
bool braid_control::EQUALITY_TEST = false;
bool braid_control::EVEN_WRITHE = true;
bool braid_control::EXPANDED_BRACKET_POLYNOMIAL = true;
bool braid_control::EXTRA_OUTPUT = false;
bool braid_control::FIXED_POINT_INVARIANT = false;
bool braid_control::FLIP_BRAID = false;
bool braid_control::FLAT_CROSSINGS = false;
bool braid_control::GAUSS_CODE = false;
bool braid_control::GCD_BACK_SUBSTITUTING = false;
bool braid_control::HAMILTONIAN = false;
bool braid_control::HC_COUNT = false;
bool braid_control::HC_EDGES = false;
bool braid_control::HC_LIST_ALL = false;
bool braid_control::HOMFLY = false;
bool braid_control::IMMERSION_CODE = false;
bool braid_control::INVERT_BRAID = false;
bool braid_control::JONES_POLYNOMIAL = false;
bool braid_control::KAMADA_DOUBLE_COVERING = false;
bool braid_control::KAUFFMAN_BRACKET = false;
bool braid_control::KNOTOID_BRACKET = false;
bool braid_control::LINE_REFLECT_BRAID = false;
bool braid_control::LONG_KNOT = false;
bool braid_control::MATRIX = false;
bool braid_control::MOCK_ALEXANDER = false;
bool braid_control::NORMALIZING_Q_POLYNOMIALS = false;
bool braid_control::NORMALIZE_BRACKET = true;
bool braid_control::NUMERATOR_GCD = false;
bool braid_control::OPGC = false; // over preferred Gaus code
bool braid_control::OU_FORMAT = false;
bool braid_control::OUTPUT_AS_INPUT = false;
bool braid_control::PARITY_ARROW = false;
bool braid_control::PARITY_BRACKET = false;
bool braid_control::PD_FORMAT = false;
bool braid_control::PEER_CODE = false;
bool braid_control::PLANE_REFLECT_BRAID = false;
bool braid_control::PRIME_WEYL = false;
bool braid_control::PRIME_TEST = false;
bool braid_control::QUANTUM_WEYL = false;
bool braid_control::QUATERNION = false;
bool braid_control::RACK_POLYNOMIAL = false;
bool braid_control::RAW_OUTPUT = false;
bool braid_control::REDUCE_BRAIDS = true;
bool braid_control::RELAXED_PARITY = false;
bool braid_control::REMOVE_REIDEMEISTER_II_MOVES = true;
bool braid_control::REMOVE_PEER_CODE_COMPONENT = false;
//bool braid_control::REAL_STUDY_DELTA_1 = false;
bool braid_control::SAWOLLEK = false;
bool braid_control::SILENT_OPERATION = false;
bool braid_control::STATUS_INFORMATION = false;
bool braid_control::STUDY_RHO_MAPPING = false;
bool braid_control::SWITCH_POLYNOMIAL_INVARIANT = false;
bool braid_control::TEST_MODE = false;
bool braid_control::TeX_POLYNOMIAL_OUTPUT = false;
bool braid_control::TMP_DIRECTORY = false;
bool braid_control::TRUNCATED_WEYL = false;
bool braid_control::T_VARIABLE = true;
bool braid_control::LPGD = false; //  left preferred Gaus data
bool braid_control::ULPGD = false; // unoriented left preferred Gaus data
bool braid_control::UOPGC = false; // unoriented over preferred Gaus code
bool braid_control::VERIFY_DELTA_0 = false;
bool braid_control::VOGEL_ALGORITHM = false;
bool braid_control::VOGEL_HEIGHT_ONLY = false;
bool braid_control::VOGEL_TURNING_NUMBER = false;
bool braid_control::WAIT_SWITCH = false;
bool braid_control::WELDED_BRAID = false;
bool braid_control::WEYL = false;
bool braid_control::ZIG_ZAG_DELTA = false;


/* Sn_matrix_top_down is true iff Sn = I^{n-1} x S x I^{k-n-1}, so
   otherwise Sn = I^{k-n-1} x S x I^{n-1}. Our default setting of false
   is consistent with numbering braids bottom-up, so S_1 places S in the bottom 
   right hand corner of a matrix.
*/

bool braid_control::Sn_matrix_top_down = false; 
bool braid_control::first_time = true;

int braid_control::REMOVE_COMPONENT=0; // used to identify a component of a peer code to be removed
int	braid_control::SWITCH_POWER=0; // used to control whether powers of switches are calculated.

bool braid_control::VOGEL_DEBUG = false;

int braid_control::HC_INCLUDE_EDGE = -1;  // initialises to the "no specific edge to include" value

int	braid_control::INFINITE_CYCLE = cycle::UNSPECIFIED; // used by turning-number calculation initialises to "no specific turning cycle" value

/* SATELLITE acts as a flag and indicates the number of strands to be added, if the default of 2 is not required */
int braid_control::SATELLITE = 0;


int braid_control::wait_threshold = 1;
int braid_control::wait_count = 0;
int braid_control::reset_count = 0; 

/********************* Function prototypes ***********************/
void braid(string input_string, string title);
bool get_next_input_string(string input_file,string& input_string, string& title);
template <class T, class St> void report_switch_matrix (const matrix<T,St>& switch_matrix, 
														const matrix<T,St>& switch_matrix_inverse, string message);
template <class T, class St> void poly_invariant(const matrix<T,St>& switch_matrix, const matrix<T,St>& switch_matrix_inverse,
                    string input_string, string title);
void help_info();
void debug_help();
void set_programme_long_option(char* cptr, string source);
void set_programme_short_option(char* cptr);
bool debug_setup(char* argv_ptr);
void check_debug_option_parameters(char* pptr, string option);
void set_cli_suboptions(char* optr, string option);
void print_prog_params (ostream& os, string level, string prefix);
Hmatrix H22inverse(Hmatrix& M);
Hpmatrix Hpmatrix_from_q_switch(const Hmatrix& q_switch, bool t_variable);
Qpmatrix adj2(Qpmatrix& M);
void simplify_matrix (Qpmatrix& M);
string parse_long_knot_input_string (string input_string);
void display_fixed_point_switch(matrix<int>& M, ostream& os, bool number_from_zero);
string vogel (generic_code_data code_data, int* turning_number_ptr=0);
			
template <typename T, typename V, typename U, typename St> 
matrix<polynomial<T,V,U>,St> C_study_matrix(const matrix<polynomial<T,V,U>,St>& H_matrix, int n=0, int m=0, int* rperm=0, int* cperm=0);

void rat_poly_invariant(const Qpmatrix& switch_matrix, const Qpmatrix& switch_matrix_inverse,
                    string input_string, string title);

template <typename T, typename V, typename U, typename St>  void normalize(polynomial<T,V,U>& poly, const matrix<T,St>& switch_matrix);
void normalize(Hpolynomial& poly, const Hpmatrix& switch_matrix);


template <typename T, typename St> 
bool braid_rep(matrix<T,St>*& Matrix_rep, const matrix<T,St>& switch_matrix, const matrix<T,St>& switch_matrix_inverse, string input_string,
                int num_terms, int num_strings);
				
template <typename T> void equality_test (T x, T y, bool A);

bool XxX_permutation (matrix<int>& U, matrix<int>& D);
bool rows_are_permutations (matrix<int>&M);
bool normal_Yang_Baxter_rules (matrix<int>& U, matrix<int>& D, bool rack_condition = false);
bool Yang_Baxter_satisfied (matrix<int>& U, matrix<int>& D, bool rack_condition = false);
bool one_dominates(matrix<int>& Ua, matrix<int>& Da, matrix<int>& Ub, matrix<int>& Db);
bool two_dominates(matrix<int>& Ua, matrix<int>& Da, matrix<int>& Ub, matrix<int>& Db);
bool power_2(matrix<int>& U, matrix<int>& D);
int switch_order(matrix<int>& U, matrix<int>& D);
void fixed_point_invariant(matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, 
						   matrix<int>& Tu, matrix<int>& Td, braid_control::ST_pair_type pair_type, string input_string, string title);

void rack_poly_invariant(matrix<int>& Su, matrix<int>& Sd, matrix<int>& invSu, matrix<int>& invSd, matrix<int>& Tu, 
						 matrix<int>& Td, braid_control::ST_pair_type pair_type, string input_string, string title, int writhe, int turning_number);

void commutative_automorphism_invariant(const Qpmatrix& phi, const Qpmatrix& psi, string input_string, string title);			
void set_acceptable_input_type();

void sigfpe_handler (int sig) 
{
	cout << "Error, braid programme received SIGFPE!" << endl;
	exit(1);
}

void generic_code(string input_string, string title);

/******************* Main Function ************************/
int main (int argc, char* argv[])
{
    string 	input_file;
    string 	output_file;
	string 	input_string;
	string 	title;
    bool    input_file_provided = false;
    bool    output_file_provided = false;

	/* take a copy of the original command line in case we're debugging and remove
	   the debug details.  We will want to record the original command line for the debug file.
	*/
	
	vector<string> original_programme_arguments(argc); 
	for (int i = 0; i< argc; i++)
		original_programme_arguments[i] = string(argv[i]);

	vector<pair<string,string> > switch_cache; // pair=(title,switch)

	signal(SIGFPE, sigfpe_handler); // register the SIGFPE handler
	time_t start_time = time(0);

    /* Determine command line options and whether input is to
       come from the screen or a file */
    if (argc == 1)
	{

		cout << "\nThis is A. Bartholomew's braid programme, v" << braid_control::version;
		cout << "\nbraid --<task> [-<short_options>][--<long_option>][<infile>[<outfile>]]" << endl;
		cout << "  braid -h to view help file" << endl;
		cout << "  braid -H! for detailed usage information" << endl;
    	exit(0);

#ifdef TAKEOUT		

		/* Start menu deprecated as of version 26 */


		/* Put up a menu for the user to decide what to do */
		cout << "\nThis is A. Bartholomew's braid programme, v" << braid_control::version;	
		cout << "\n\nPlease choose from the following options:\n";
		cout << "\n 0. Quit the program without doing anything.";
		cout << "\n 1. Burau matrix representation of a braid.";
		cout << "\n 2. Alexander matrix representation of a braid.";
		cout << "\n 3. Quaternionic matrix representation of a braid.";
		cout << "\n 4. Matrix-switch representation of a braid.";
		cout << "\n 5. Weyl algebra switch representation of a braid.";
		cout << "\n 6. Commutative automorphism switch representation of a peer code.";
		cout << "\n 7. Sawollek's normalized Conway polynomial of a braid.";
		cout << "\n 8. HOMFLY polynomial for the closure of a braid.";
		cout << "\n 9. Finite biquandle fixed point invariant of a braid.";
		cout << "\n10. Rack-polynomial invariant of a braid.";
		cout << "\n11. Kauffman backet polynomial of a knot, link, or (multi-)knotoid.";
		cout << "\n12. Jones polynomial of a knot, link, or (multi-)knotoid.";
		cout << "\n13. Affine index polynomial of a knot or link.";
		cout << "\n14. Turaev extended bracket polynomial of a knotoid.";
		cout << "\n15. Arrow polynomial of a knot, link or (multi-)knotoid.";
		cout << "\n16. Parity bracket polynomial of a knot, link or (multi-)knotoid.";
		cout << "\n17. Parity arrow polynomial of a knot, link or (multi-)knotoid.";
		cout << "\n18. Dynnikov test for the trivial braid.";
		cout << "\n19. Dowker code for the closure of a braid or an immersion code.";
		cout << "\n20. Gauss code for the closure of a braid.";
		cout << "\n21. Immersion code for the closure of a braid.";
		cout << "\n22. Peer code for the closure of a braid or an immersion code.";
		cout << "\n23. Peer code for the r-parallel cable satellite of a knot's peer code.";
		cout << "\n24. Display information about a braid.";
		cout << "\n25. Vogel's algorithm to determine a braid word from a knot or link";
		cout << "\n h. View the help screens.";

		bool invalid_input;			
		do
		{		
			invalid_input = false;			
			cout << "\n\nbraid: ";
			string selection;
			getline(cin, selection);
		
			if (selection == "0" || selection == "q" || selection == "Q")
			    return 0;
			else if (selection == "1")
			{
			    braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
				braid_control::BURAU = true;
			}
			else if (selection == "2")
			{
			    braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
				braid_control::ALEXANDER = true;
			}
			else if (selection == "3")
			{		
				braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
				braid_control::QUATERNION = true;
			}
			else if (selection == "4")
			{		
				braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
				braid_control::MATRIX = true;
				polynomial_control::WAIT_INFO = true;
			}
			else if (selection == "5")
			{		
				braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
				braid_control::WEYL = true;
				polynomial_control::WAIT_INFO = true;
			}
			else if (selection == "6")
			{		
				braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
				braid_control::COMMUTATIVE_AUTOMORPHISM = true;
				polynomial_control::WAIT_INFO = true;
			}
			else if (selection == "7")
			    braid_control::SAWOLLEK = true;
			else if (selection == "8")
			    braid_control::HOMFLY = true;
			else if (selection == "9")
			{
			    braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
			    braid_control::FIXED_POINT_INVARIANT = true;
			}
			else if (selection == "10")
			{
			    braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
				braid_control::RACK_POLYNOMIAL = true;
			}
			else if (selection == "11")
				braid_control::KAUFFMAN_BRACKET = true;
			else if (selection == "12")
				braid_control::JONES_POLYNOMIAL = true;
			else if (selection == "13")
				braid_control::AFFINE_INDEX = true;
			else if (selection == "14")
				braid_control::KNOTOID_BRACKET = true;
			else if (selection == "15")
			    braid_control::ARROW_POLYNOMIAL = true;
			else if (selection == "16")
			    braid_control::PARITY_BRACKET = true;
			else if (selection == "17")
			    braid_control::PARITY_ARROW = true;
			else if (selection == "18")
			    braid_control::DYNNIKOV_TEST = true;
			else if (selection == "19")
			    braid_control::DOWKER_CODE = true;
			else if (selection == "20")
			    braid_control::GAUSS_CODE = true;
			else if (selection == "21")
				braid_control::IMMERSION_CODE = true;
			else if (selection == "22")
			    braid_control::PEER_CODE = true;
			else if (selection == "23")
			{
			    braid_control::SATELLITE = 2;
			    braid_control::PEER_CODE = true;
			}
			else if (selection == "24")
			    braid_control::STATUS_INFORMATION = true;
			else if (selection == "25")
			    braid_control::VOGEL_ALGORITHM = true;
			else if (selection == "h" || selection == "help")
			    help_info();
			else
			{
				cout << "\nInvalid input, please select again.";
				invalid_input = true;
			}
		} while (invalid_input);
#endif
	}
	else
    {
		for (int i=1; i < argc; i++)
		{
			if (*(argv[i]) == '-')
			{
				if (*(argv[i]+1) == '-')
					set_programme_long_option(argv[i], "command line");
				else 
					set_programme_short_option(argv[i]);
			}
			else if (!input_file.length()) 
			{
				input_file = argv[i];
	    		input_file_provided = true;
			}
			else if (!output_file.length()) 
			{
				output_file = argv[i];
		    	output_file_provided = true;
			}
			else
			{
		    	cout << "Usage braid --<task> [-<options>][<infile>[<outfile>]], type braid -h for help.\n";
		    	exit(0);
			}
		}	
    }

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "full command line: ";
	for (int i = 0; i< argc; i++)
		debug << original_programme_arguments[i] << ' ';
	debug << endl;
	debug << boolalpha;
}

	if (braid_control::DEVELOPMENT_MODE)
	{
		
/*
		list<int> test;
		cout << "terms: ";
		do 
		{
			int next;
			cin >> next;
			if (next < 0)
				break;
			else
				test.push_back(next);
		} while(true);
		
		cout << "lcm = " << least_common_multiple(test) << endl;
exit(0);		
*/
		
/*
vector<int> braid_permutation (string braid, int num_terms, int num_strings);
	

void renumber_peer_code(generic_code_data& code_data, vector<int> shift);
matrix<int> permutation_cycles (generic_code_data& code_data);

bool gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data, bool optimal=true, vector<int>* gauss_crossing_perm=0, bool evaluate_gauss_crossing_perm=false);		

//		input_control::ACCEPT_MAP |= input_control::planar_diagram;
//		input_control::ACCEPT_MAP |= input_control::gauss_code;

//if (debug_control::DEBUG >= debug_control::SUMMARY)
//	debug << "development: input_control::ACCEPT_MAP set to 0x" << hex << input_control::ACCEPT_MAP << dec << endl;
*/
		
		string buffer;
/*		
		string title;
    	input.open(input_file.c_str());
    	if(!input)
    	{
			cout << "\nError opening input file\n";
			exit(0);
    	}	
		
		while (getline(input,buffer))
		{
			if (buffer[0] == ';' || (buffer[0] == '-' && buffer[1] == '-'))
				continue;
				
			cout << buffer << endl;
*/
		do
		{
//			cout << "enter peer code: ";
			cout << "enter braid: ";
			getline(cin, buffer);
			
			
			if (buffer == "q")
				break;
	
//			generic_code_data code_data;
//			read_peer_code(code_data,buffer);
			
//			seifert_circle_data SCD(code_data);
		
//			seifert_circle_data SCD(buffer);
			
//			print_seifert_data(cout,SCD," ");

//if (debug_control::DEBUG >= debug_control::SUMMARY)
//	print_seifert_data(debug,SCD," ");
			
/*			vector<int> perm = braid_permutation(buffer,num_terms,num_strings);

			cout << "  braid_permutation: ";
			for (size_t i=0; i < perm.size(); i++)
				cout << perm[i] << ' ';
			cout << endl;
			
			int turning_number = 1;
			for (size_t i=0; i < perm.size(); i++)
			{
				for (size_t j=i+1; j < perm.size(); j++)
				{
					if (perm[j] < perm[i])
						turning_number++;
				}
			}

			cout << "  turning_number = " << turning_number << endl;
*/			
		} while(true);
		
//continue;				
exit(0);

#ifdef TAKEOUT
			generic_code_data code_data;
			read_code_data(code_data,buffer);
	
			int num_cycles;
			int num_left_cycles;
			matrix<int> cycle(code_data.num_crossings+2, 2*code_data.num_crossings+1);
			if (calculate_turning_cycles(code_data, cycle, num_left_cycles, num_cycles))
			{
			    int genus = (2 + code_data.num_crossings - num_cycles)/2;

			    cout << "genus = " << genus << endl;

                            cout << "turning_cycles:" << endl;
                            for (int i=0; i<num_cycles; i++)
		            {
		                 cout << "  cycle " << i << ": ";
	                         for (int j=1; j<=cycle[i][0]; j++)
	                            cout << cycle[i][j] << " ";

	                         cout << endl;
			    }
			}
			else
			{
			    cout << "failed to calculate turning cycles" << endl;
			}

			continue;

			for (int i=0; i<=2*code_data.num_crossings; i++)
			{
				generic_code_data shifted_code_data = code_data;
				vector<int> shift_vector = {i};
				renumber_peer_code(shifted_code_data,shift_vector);
				
				
				vector<int> perm(code_data.num_crossings);
				vector<int> inv_perm(code_data.num_crossings);
				
				for (int j=0; j< code_data.num_crossings; j++)
				{
					int odd_terminating = shifted_code_data.code_table[generic_code_data::table::EVEN_TERMINATING][j]+1;
					int even_terminating = shifted_code_data.code_table[generic_code_data::table::ODD_TERMINATING][j]+1;
					int crossing = even_terminating/2;
					perm[crossing-1] = (odd_terminating+1)/2;

					if (shifted_code_data.code_table[generic_code_data::table::TYPE][(odd_terminating-1)/2]==generic_code_data::type::TYPE2)
						perm[crossing-1]*=-1;
				}
					
				
				cout << "(";
				for (int j=0; j< code_data.num_crossings; j++)
				{
					cout << setw(3) << j+1;
//					if (j< code_data.num_crossings-1)
//						cout << ' ';
				}
				cout << ")";
				
				cout << "\t";
				write_code_data(cout,shifted_code_data);
				cout << "\n(";

				for (int j=0; j< code_data.num_crossings; j++)
				{
					cout << setw(3) << perm[j];
//					if (j< code_data.num_crossings-1)
//						cout << ' ';
				}
				cout << ")\n";
				
				cout << endl;				
				
				matrix<int> peer_perm_cycles = permutation_cycles(shifted_code_data);
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "development: peer_perm_cycles: " << endl;
	for (unsigned int i=0; i< peer_perm_cycles.numrows(); i++)
	{
		debug << "development:   cycle " << i << "(" << peer_perm_cycles[i][0] << "): ";
		for (int j=1; j<= peer_perm_cycles[i][0]; j++)
			debug << peer_perm_cycles[i][j] << ' ';
		debug << endl;
	}
}
				if (peer_perm_cycles.numrows() == 1 && i == 1)
				{
					cout << "single cycle ";
					write_code_data(cout,shifted_code_data);				
//					cout << " shift " << i << endl;
					cout << endl;
				}

			}
			continue;
		
/*			bool connected_sum = !three_connected(code_data,true);
			
			if (connected_sum)
				cout << "non_prime" << endl;
			else
				cout << "prime" << endl;
*/				

/*
			simple_Reidemeister_II_return simple_Reidemeister_II_data = simple_Reidemeister_II_present(code_data, false); // create_flags = false
			if (simple_Reidemeister_II_data.Reidemeister_II_found)
				cout << "non_minimal" << endl;
			else
				cout << "minimal" << endl;
*/			
				
			
			write_code_data(cout,code_data);
			cout << endl;
/*
			generic_code_data peer_code_data;
			gauss_to_peer_code(code_data, peer_code_data);			
			write_code_data(cout,peer_code_data);
//			print_code_data(cout,peer_code_data);
			cout << endl;
*/

continue;			
		//	int num_cycles;
		//	int num_left_cycles;
	//		matrix<int> cycle(code_data.num_crossings+2, 2*code_data.num_crossings+1);
			calculate_turning_cycles(code_data, cycle, num_left_cycles, num_cycles);
			vector<int> cycle_length_count(code_data.num_crossings+2);
			
			for (int i=0; i< num_cycles; i++)
				cycle_length_count[cycle[i][0]]++;
			
			write_code_data(cout,code_data);
			cout << " num_crossings = " << code_data.num_crossings << ": ";
		
			for (int i=3; i<code_data.num_crossings+2; i++)
			{
				cout << cycle_length_count[i];
				
				bool more = false;
				for (int j=i+1; j< code_data.num_crossings+2; j++)
				{
					if (cycle_length_count[j] !=0)
					{
						more = true;
						break;
					}
				}
				
				if (more)
					cout << ',';
				else
					break;
			}
			cout << endl;
		}
		
exit(0);

		string input_string;// = "[5 11 7 1, 13 3 15 9]/+ * + + * * # #";
		cout << "Code data: ";
		getline(cin,input_string);
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "development: input_string = " << input_string << endl;
	
		generic_code_data code_data;
		read_peer_code(code_data,input_string);
	
//void print_code_data(generic_code_data& code_data, ostream& s, string prefix);
//void write_peer_code(ostream& s, const generic_code_data& code_data, bool zig_zags=false, bool labelled=true);
bool remove_edge_flags_from_peer_code(generic_code_data& code_data, vector<int>& edge_flag);
int remove_Reidemeister_II(generic_code_data& code_data, vector<int>& component_flags);
//int vertex_span (generic_code_data& peer_code_data, int initial_crossing, vector<int>* exclude=0);
//bool valid_knotoid_input(generic_code_data& code_data);
generic_code_data partition_peer_code(generic_code_data& code_data);
/*
		matrix<int> input_zig_zag_count(2, code_data.num_crossings);
		
		cout << "zig-zag counts: ";
			
		for (int i=0; i< 2; i++)
		for (int j=0; j< code_data.num_crossings; j++)
			cin >> input_zig_zag_count[i][j];
			
		code_data.zig_zag_count = input_zig_zag_count;
		
		cout << "head zig-zag count: ";
		cin >> code_data.head_zig_zag_count;

		valid_knotoid_input(code_data);
*/	
/*		cout << "peer code: ";
		write_peer_code(cout,code_data);
		cout << endl;
		print_code_data(code_data, cout,"");
		cout << endl;
		
		gauss_orientation_data gauss_data(code_data);
		cout << "Gauss data: ";
		write_gauss_data(gauss_data,cout,true,code_data.immersion,code_data.head_zig_zag_count);
		cout << endl;
*/		
		int num_cpt_flags;
		cout << "num component flags: ";
		cin >> num_cpt_flags;
			
		vector<int> component_flags(num_cpt_flags);
		cout << "component flags: ";
		
		for (int j=0; j< num_cpt_flags; j++)
			cin >> component_flags[j];
		
		valid_knotoid_input(code_data);
		remove_Reidemeister_II(code_data,component_flags);
		
		cout << "after removing Reidemeister II, component flags: ";
		
		for (int j=0; j< num_cpt_flags; j++)
			cout  << component_flags[j] << ' ';
		cout << endl;
		
		
/*
		gauss_orientation_data reduced_gauss_data(code_data);
		cout << "reduced Gauss data: ";
		write_gauss_data(reduced_gauss_data,cout,true,code_data.immersion,code_data.head_zig_zag_count);
		cout << endl;
*/

/*
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "Initial code_data:" << endl;
	print_code_data(code_data, debug,"");
	debug << endl;
}		
		generic_code_data subsequent_code_data = partition_peer_code(code_data);
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
		debug << "partitioned code_data:" << endl;
		print_code_data(code_data, debug,"");
		debug << endl;
		debug << "subsequent code_data:" <<  endl;
		print_code_data(subsequent_code_data, debug,"");
		debug << endl;
}
*/
		exit(0);
#endif
		
	} /* END OF DEVELOPMENT MODE CLAUSE */
	
	
	if (input_file_provided)
	{
    	input.open(input_file.c_str());
    	if(!input)
    	{
			cout << "\nError opening input file\n";
			exit(0);
    	}	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::main: input file: " << input_file << "\n";

		/* Read any options included in the input file */
		string next_line;
		while(getline(input,next_line))
		{
			char* line_buf = c_string(next_line);
			char* cptr = strchr(line_buf,';');
			if (cptr)
	    		*cptr = '\0';

			if (strlen(line_buf))
			{
				cptr = strchr(line_buf,'[');
				if (cptr && !strchr(line_buf,'\\') && !strchr(line_buf,'/') && next_line.find("X[")== string::npos && next_line.find("S[") == string::npos && next_line.find("s[") == string::npos) 
				{		
					/* this is an options line and not a line containing a peer code, a planar diagram or a switch specifying the t-variable.  Make sure there's an end to the option definitions on this line 
					*/
					if (!strchr(line_buf,']'))
					{
						cout << "\nError in " << input_file << ": programme options must be enclosed within []\n";
						exit(0);
					}
					do
					{
						cptr++;
						/* skip whitespace */
						while (isspace(*cptr))
							cptr++;
						set_programme_long_option(cptr++, "input file");
						while (*cptr != ',' && *cptr != ']')
							cptr++;
					} while (*cptr != ']');
				}
			}
			
			delete[] line_buf;
    	}
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
	print_prog_params(debug, "detail", "print_prog_params: ");

    /* prepare output file */
    if (!output_file_provided && braid_control::TMP_DIRECTORY)
		output_file = "/tmp/braid.out";
	else if (!output_file_provided)
		output_file = "braid.out";
	
    output.open (output_file.c_str());

    if (!output)
    {
		cout << "\nError opening output file " << output_file << endl;
//		system("pwd");
//		system ("ps -ef | grep braid");
		exit(0);
    }
    else
    {
		if (!braid_control::RAW_OUTPUT)
		{
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Output from braid v" << braid_control::version << "\n";

			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Command line: ";
			for (int i = 0; i< argc; i++)
				output << original_programme_arguments[i] << ' ';
			output << "\n";

/*			if (input_file_provided)
			{
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Input file: " << input_file << "\n";
			}
*/
		}
    }


    if (!input_file_provided)
    {
	
		if (braid_control::QUATERNION && braid_control::MATRIX)
		{
			cout << "\nYou need to provide an input file if you want to specify both the quaternionic";
			cout << "\nand matrix-switch options together.  With both, the default switch cannot be detrmined.\n" << endl;
			exit(0);
		}

		if (argc > 1)
		{
			cout << "\nThis is A. Bartholomew's braid programme, v" << braid_control::version;
			cout << "\n\nYou have started the programme to carry out the task described below.\n";
			cout << "If you are unsure whether this is what you intended type the word help at the prompt.\n";
			cout << "Type exit at the prompt to leave the programme.\n";
		}
		
		cout << "\n\nThe programme will";

		if (braid_control::VOGEL_TURNING_NUMBER)
		{
		    cout << " determine the turning number of the planar diagram\n";
	    	cout << "determined by the labelled peer codes that you enter.\n";
		}
		else if (braid_control::VOGEL_ALGORITHM)
		{
		    cout << " determine a braid word representation for the links\n";
	    	cout << "determined by the labelled peer codes that you enter.\n";
		}
		else if (braid_control::SAWOLLEK)
		{
		    cout << " evaluate Sawollek's normalized Conway polynomial\n";
	    	cout << "for the links determined by the braid words that you enter.\n";
		}
		else if (braid_control::HOMFLY)
		{
		    cout << " evaluate the HOMFLY polynomial for the closure";
	    	cout << "\nof the classical braid words that you enter.\n";
		}
		else if (braid_control::HAMILTONIAN)
		{
		    cout << " evaluate a Hamiltonian circuit for the shadow of the diagram specified";

		}
		else if (braid_control::RACK_POLYNOMIAL)
		{
		    cout << " evaluate the rack-polynomial invariant for the";
	    	cout << "\nbraid words that you enter.\n";
		}
		else if (braid_control::JONES_POLYNOMIAL)
		{
		    cout << " evaluate the Jones polynomial for knot, link or (multi-)knotoid ";
	    	cout << "\nspecified by a labelled peer code, labelled immersion code";
	    	cout << "\nor Gauss code that you enter.\n";
		}
		else if (braid_control::KAUFFMAN_BRACKET)
		{
		    cout << " evaluate the Kauffman bracket polynomial for a knot, link or (multi-)knotoid ";
	    	cout << "\nspecified by a labelled peer code, labelled immersion code or (for knots and knotoids only)";
	    	cout << "\nGauss code that you enter.\n";
		}
		else if (braid_control::KNOTOID_BRACKET)
		{
		    cout << " evaluate Turaev's extended bracket polynomial for a knotoid specified";
	    	cout << "\nby a peer code or immersion code that you enter.\n";
		}
		else if (braid_control::ARROW_POLYNOMIAL)
		{
		    cout << " evaluate the arrow polynomial for a knot, link or (multi-)knotoid ";
	    	cout << "\nspecified by a labelled peer code, labelled immersion code or (for knots and knotoids";
	    	cout << "\n only) Gauss code that you enter.\n";
		}
		else if (braid_control::PARITY_BRACKET)
		{
		    cout << " evaluate the parity bracket polynomial for a knot, link or (multi-)knotoid ";
	    	cout << "\nspecified by a labelled peer code or labelled immersion code that you enter.\n";
		}
		else if (braid_control::PARITY_ARROW)
		{
		    cout << " evaluate the parity arrow polynomial for a knot, link or (multi-)knotoid ";
	    	cout << "\nspecified by a labelled peer code or labelled immersion code that you enter.\n";
		}
		else if (braid_control::FIXED_POINT_INVARIANT)
		{
		    cout << " evaluate the fixed point invariant for the closure of the braid words that you enter.\n";
		}
		else if (braid_control::DYNNIKOV_TEST)
		{
		    cout << " determine whether the braids described by the braid words that you enter are trivial.\n";
		}
		else if (braid_control::DOWKER_CODE)
		{
	    	cout << " determine the Dowker code of the closure of the braids that you enter, provided they are knots.\n";
		}
		else if (braid_control::GAUSS_CODE || braid_control::IMMERSION_CODE || braid_control::PEER_CODE)
		{
			cout << " determine the ";
			if (braid_control::GAUSS_CODE)
				cout << "Gauss";
			else if (braid_control::IMMERSION_CODE)
				cout << "immersion";
			else
				cout << "peer";
				
			cout<< " code of the closure";
			cout << "\nof the braids ";
			if (braid_control::PEER_CODE) 
				cout << "or shifted peer codes, Gauss codes, immersion codes, or planar diagrams ";
			cout << "that you enter";
			cout << (braid_control::IMMERSION_CODE?  ", provided they are knots.\n":".\n");
		}
		else if (braid_control::AFFINE_INDEX)
		{
		    cout << " evaluate the affine index polynomial for a classical";
	    	cout << "\nor virtual knot specified by a braid word, a labelled peer code,";
	    	cout << "\nlabelled immersion code or Gauss code that you enter.\n";
		}
		else if (braid_control::DISPLAY_DELTA_1_ONLY)  
		{                          
		    cout << " evaluate the";

	    	if (braid_control::QUATERNION || braid_control::MATRIX)
				cout << " 1st ideal ";
	    	else
				cout << " Alexander ";

    		cout << "polynomial of the links\ndetermined by the braid words that you enter.\n";
		}
		else if (braid_control::SATELLITE)
		{
	    	cout << " determine the labelled peer code of the r-parallel cable";
	    	cout << "\nsatellite of the peer codes that you enter, provided they are knots.\n";
		}
		else
		{
		    /* full SWITCH_POLYNOMIAL_INVARIANT, either Burau or quaternionic */
			cout << " evaluate ";
		    if (braid_control::QUATERNION)
				cout << "a quaternionic matrix ";
		    else if (braid_control::MATRIX)
				cout << "a matrix-switch ";
		    else if (braid_control::WEYL)
				cout << "a Weyl-algebra-switch ";
		    else if (braid_control::COMMUTATIVE_AUTOMORPHISM)
				cout << "a commutative automorphism switch ";
		    else
				cout << "the Burau or Alexander matrix ";
	    	cout << "representation, M, of\n";
			cout << "each braid word or peer code that you enter, then calculate the \n";
	    	if (braid_control::ALEXANDER || braid_control::BURAU) 
				 cout << "1st ideal polynomial, and the 0th ideal polynomial in the virtual case." << endl;
			else
			{
				cout << "0th";

		    	if (!braid_control::MATRIX && !braid_control::WEYL)
			    	cout << " and 1st ideal polynomials." << endl;
	    		else
		    		cout << " ideal polynomial and 1st ideal polynomial generators." << endl;
			}
		}
		if (argc == 1)
			cout <<"\nType help at the prompt to view the help screens.\n";

		if (braid_control::VOGEL_ALGORITHM)
		{
			cout << "\nEnter a labelled peer code\n";
			cout << "\nE.g. [9 19 -15, 1 -17, -7 -21 -5 13, -11 3]/- - + - - - - + + - -";
		}
		if (braid_control::SATELLITE)
		{
			cout << "\nEnter a labelled peer code for a knot\n";
			cout << "\nE.g. [-3 7 -11 -9 1 -5]/+ + + + + +";
		}
		else if (braid_control::DOWKER_CODE || braid_control::HOMFLY)
		{
			cout << "\nEnter the word in si";		
			cout << "\nE.g. s1-s2-s3-s3";
		}
		else if (braid_control::HAMILTONIAN)
		{
			cout << "\nEnter a braid, peer code, Gauss code, planar diagram or braid of a classical knot or link";		
		}
		else if (braid_control::KAUFFMAN_BRACKET || braid_control::JONES_POLYNOMIAL || braid_control::ARROW_POLYNOMIAL 
		         || braid_control::PARITY_BRACKET || braid_control::PARITY_ARROW)
		{
			cout << "\nEnter a labelled peer code, a labelled immersion code or a Gauss code\n";
			cout << "\nE.g. [-3 7 -11 -9 1 -5]/+ + + + + + or (-0 -2)(1 4)(-3 -5) / + + + + + +"
			     << "\nor -1 2 -3 1 -2 3 / + + +";
		}
		else if (braid_control::AFFINE_INDEX)
		{
			cout << "\nEnter a braid word, a labelled peer code, a labelled immersion code or a Gauss code\n";
			cout << "\nE.g. s1-s2-s3-s3, [-3 7 -11 -9 1 -5]/+ + + + + +," 
			     << "\n(-0 -2)(1 4)(-3 -5) / + + + + + + or -1 2 -3 1 -2 3 / + + +";
		}
		else if (braid_control::KNOTOID_BRACKET)
		{
			cout << "\nEnter a labelled peer code or a labelled immersion code for a knotoid\n";
			cout << "\nE.g. [-3 7 -11^ -9 1 -5]/+ + + + + + or (-0 -2^)(1 4)(-3 -5) / + + + + + +";
		}
		else
		{
			cout << "\nEnter the word in s and t characters";		
			
			if (braid_control::WEYL)
				cout << "\nE.g. s1s2t1s3s3";
			else
				cout << "\nE.g. s1-s2t1-s3-s3";
		}
    }
    
    set_acceptable_input_type();
	
try {

	if (braid_control::SWITCH_POLYNOMIAL_INVARIANT)
	{	
		bool		fundamental_equation_satisfied;
	
		/* If a file has been provided we first scan the input file to see if it contains 
		   any relevant switches, and it so we add them to the switch_cache.
	    */
	    if ((braid_control::QUATERNION || braid_control::MATRIX || braid_control::WEYL || braid_control::FIXED_POINT_INVARIANT || 
	         braid_control::RACK_POLYNOMIAL  || braid_control::COMMUTATIVE_AUTOMORPHISM ) && input_file_provided
	       )
	    {
			string next_line;
		
			/* reset input ready for read */
			input.clear(); // state flags
			input.seekg(0);
			while(getline(input,next_line))
			{
				if (next_line.find(';') != string::npos)
					next_line.erase(next_line.find(';'));

//if (debug_control::DEBUG >= debug_control::DETAIL)
//	debug << "braid::main: read " << next_line<< " from input file" << endl;
				
				if (next_line.length())
				{				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::main: processing input file line: " << next_line << endl;

					if (next_line[0] == '-' && next_line[1] == '-')
					{
						title = next_line;
					}
					else if ((next_line.find('S') != string::npos || next_line.find('s') != string::npos ) && 
					     next_line.find('=') != string::npos &&	next_line.find('w') == string::npos)
					{
						if (next_line.find('A') != string::npos) 
						{
							if (braid_control::COMMUTATIVE_AUTOMORPHISM)
							{
								switch_cache.push_back(pair<string,string>(title,next_line));
								title.clear();
							}
						}
						else if (next_line.find('F') != string::npos) 
						{
							/* We check F before W since welded switches (S = FW[T]) and doodle switches (S = FD[T]) are subtypes of F */
							if (braid_control::FIXED_POINT_INVARIANT || braid_control::RACK_POLYNOMIAL)
							{
								switch_cache.push_back(pair<string,string>(title,next_line));
								title.clear();
							}
						}
						else if (next_line.find('W') != string::npos)						
						{
							if (braid_control::WEYL)
							{
								switch_cache.push_back(pair<string,string>(title,next_line));
								title.clear();
							}
						}
						else
						{

							/* to distinguish between matrix and quaternion switches, we note that quaternion
							   switch specifications only contain commas within quaternions whlst matrix
							   switches do not contain commas within parentheses (if parentheses are used at all)
							   since with matrix switches parentheses are used only to delimit polynomials.
							*/
							char* c_next_line = c_string(next_line);
							char* cptr = c_next_line;
							int level = 0;

							do
							{
								if (*cptr == '(')
									level++;
								else if (*cptr == ')')
									level--;
								else if (!level && *cptr == ',')
									break;

								cptr++;
							} while (*cptr != '\0');

							/* when we break, if *cptr == ',' we have a matrix switch 
							   otherwise *cptr = 0 and we have a quaternionic switch
							*/

							if (braid_control::MATRIX && *cptr)
							{
								switch_cache.push_back(pair<string,string>(title,next_line));
								title.clear();

//if (debug_control::DEBUG >= debug_control::DETAIL)
//	debug << "braid::main: added matrix switch " << next_line<< " to switch cache" << endl;
							}
							else if (braid_control::QUATERNION && !(*cptr))
							{
								switch_cache.push_back(pair<string,string>(title,next_line));
								title.clear();
//if (debug_control::DEBUG >= debug_control::DETAIL)
//	debug << "braid::main: added quaternionic switch " << next_line<< " to switch cache" << endl;
							}
							delete[] c_next_line;
						}
					}
				}
	    	}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::main: inputfile contains " << switch_cache.size() << " switches" << endl;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	if (switch_cache.size())
	{
		debug << "braid::main: switch cache contains: " << endl;
		for (unsigned int i=0; i < switch_cache.size(); i++)
		{
			debug << "braid::main:   " << (switch_cache[i].first == "" ? "untitled" : switch_cache[i].first ) << endl;
			debug << "braid::main:   " << switch_cache[i].second << endl;
		}
	}
	else
		debug << "braid::main: switch cache is empty" << endl;	
}
	
		}

    	for (unsigned int s=0; s< (switch_cache.size()? switch_cache.size(): 1) ; s++)
    	{	
			string next_switch;
			string switch_title;
//			first_time = true; // with this switch
			
			if (switch_cache.size())
			{
				/* get the s-th switch from the switch_cache */	
				next_switch = switch_cache[s].second;
				switch_title = switch_cache[s].first;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::main: next switch from switch_cache: " << next_switch << endl;

			}
			else if (braid_control::QUATERNION)
			{
				// set next switch to the Budapest switch, T_VARIABLE is true by default
				next_switch = "(1,1,0,0) (0,0,-1,0)";
			}
			else if (braid_control::MATRIX)
			{
				// set next switch to the 3-variable polynomial switch
//				next_switch =  "1, -1, 1, 0, 1+x+x^2, 1+x, x, 1"; //primrose
//				next_switch =  "1, -1, 1, 0, 1+xy+x^2/y, y+x, x, y"; //rocky
				next_switch =  "2, 0, x, 2, y, z, xyz-2y^2-2/2z, xz-2y/2"; // 3-variable
				braid_control::T_VARIABLE = false;
			}
			else if (braid_control::WEYL)
			{
				// set next switch to the prime 3x3 Weyl switch
				next_switch = "WP,3,3";
				braid_control::T_VARIABLE = false;
			}
			else if (braid_control::COMMUTATIVE_AUTOMORPHISM)
			{
				next_switch = "A x,0,0,y,z,0,0,w";
			}
			else if (braid_control::FIXED_POINT_INVARIANT || braid_control::RACK_POLYNOMIAL)
			{
				/* The default essential virtual pair is S = BQ^3_{6} T = BQ^3_{7} */
				next_switch = "FT 1 0 2 0 2 1 2 1 0 1 2 0 1 2 0 1 2 0 1 2 0 1 2 0 1 2 0 2 0 1 2 0 1 2 0 1"; 
//				next_switch = "F 1 0 2 0 2 1 2 1 0 1 2 0 1 2 0 1 2 0";  // as above but with the default twist
			}


			if (input_file_provided)
			{
				/* reset input ready for calls to get_input_word */
				input.clear(); // state flags
				input.seekg(0);
			}

			if (braid_control::QUATERNION)
			{
				/* Set the scalar variant to rational<bigint> unless we are already doing mod-p (this is the
				   default variant of scalar at the time of writing, but setting it here future-proof the code.
				*/
				if (!braid_control::CALCULATE_MOD_P)
					scalar::set_variant(scalar::BIGRATIONAL);	
			
				istringstream ss(next_switch);
				char ch;
				do 
				{
					ss >> ch;
				}while (ch != '(');
				ss.putback(ch);
			
                /* If an old format input file is used, C and D specified may be specified 
				   in the file, but we ignore them here and calculate them ourselves
				*/
				
				Hmatrix q_switch(2,2);
				ss >> q_switch[0][0] >> q_switch[0][1];

				Quaternion& A = q_switch[0][0];
				Quaternion& B = q_switch[0][1];
				
				typedef Quaternion Q;

				/* Check that A, B, and 1-A are invertible */
				
				Quaternion A1 = Q(1)-A;
				
				Quaternion invA;
				Quaternion invB;
				
				if (A != Q(0) && B != Q(0) && A1 != Q(0))
				{
					invA = Q(1)/A;
					invB = Q(1)/B;

					Quaternion AB = invA * invB * A * B;
					Quaternion A1AB = A1 * AB;
				
					/* Check the fundamental equation  [B,(1-A)(A,B)] = 0
					   Here, (X,Y) = invX * invY * X * Y and [X,Y] = X*Y - Y*X
					   we therefore check B * (1-A)(A,B) == (1-A)(A,B) * B
					*/
					if (B*A1AB == A1AB*B)
						fundamental_equation_satisfied = true;
					else
			    		fundamental_equation_satisfied = false;
				}
				else
		    		fundamental_equation_satisfied = false;
				
				if (fundamental_equation_satisfied || braid_control::BYPASS_FUNDAMENTAL_EQUATION_CHECK)
				{
				    Quaternion C = invA * invB * A * A1;
					Quaternion D = Q(1) - invA * invB * A * B;

					if (braid_control::EQUALITY_TEST)
					{
						equality_test(A,D,true);
						equality_test(B,C,false);
					}

					q_switch[1][0] = C;
					q_switch[1][1] = D;
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "braid::main: fundamenatal equation satisfied:\n" << "A = " << A << " B = " << B << " C = " << C << " D = " << D << endl;

					Quaternion det = determinant(q_switch);
					if (det != Quaternion(0))
					{
				    	/* create the switch_matrix and its inverse from q_matrix and, if we're 
						   using the t-variable, multiply  B and invB by t, and C and invC by t^-1 */
						Hpmatrix switch_matrix = Hpmatrix_from_q_switch(q_switch, braid_control::T_VARIABLE);				

						Hmatrix q_switch_inverse = H22inverse(q_switch);
					
						Hpmatrix switch_matrix_inverse = Hpmatrix_from_q_switch(q_switch_inverse, braid_control::T_VARIABLE);


						if (switch_cache.size() && switch_cache[s].first != "")
						{
							if (!braid_control::SILENT_OPERATION)
								cout << "\n\n" << switch_cache[s].first << endl;
							if (!braid_control::RAW_OUTPUT)
								output << "\n\n" << switch_cache[s].first << endl;
						}
						else
							report_switch_matrix(switch_matrix, switch_matrix_inverse, "quaternionic");
						
						if (braid_control::SWITCH_POWER)
						{
							Hpmatrix power = switch_matrix;
							for (int i=1; i< braid_control::SWITCH_POWER; i++)
								power *= switch_matrix;
							
							if (!braid_control::SILENT_OPERATION)
								cout << "\nswitch matrix raised to the power " << braid_control::SWITCH_POWER << ":\n" << power << endl;
	
							if (!braid_control::RAW_OUTPUT)
							{
								output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
								output << "switch matrix raised to the power " << braid_control::SWITCH_POWER << ":";
								for (unsigned int i=0; i < power.numrows(); i++)
								{
									output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
									for (unsigned int j=0; j < power.numcols(); j++)
										output << power[i][j] << " ";
								}
							}
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::main: switch power = " << braid_control::SWITCH_POWER << "\n" << power << endl;
	
						}

						while (get_next_input_string(input_file,input_string, title))
							poly_invariant(switch_matrix, switch_matrix_inverse, input_string, title);			
					
					}
				}
				else
				{
					if (!braid_control::SILENT_OPERATION)
					{
						cout << "\nSwitch elements " << A << " " << B;
						cout << "\ndo not satisfy the unit and fundamental equation requirements .\n";
					}
	
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Switch elements" << A << " " << B;
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "\ndo not satisfy the unit and fundamental equation requirements .\n";
					}
		    		output.flush();

if (debug_control::DEBUG >= debug_control::BASIC)	
{
	debug << "braid::main: switch elements " << A << " " << B << endl;
	debug << "braid::main: do not satisfy the unit and fundamental equation requirements ." << endl;
}

				}
    		}
			else if (braid_control::MATRIX) 
			{
				braid_control::T_VARIABLE = false;
				
				/* Set the scalar variant to bigint unless we are already doing mod-p */
				if (!braid_control::CALCULATE_MOD_P)
					scalar::set_variant(scalar::BIGINT);	

				int switch_terms = count(next_switch.begin(),next_switch.end(),',')+1;				
				int matrix_size = static_cast<int>(sqrt(float(switch_terms/2)));
				
				/* test matrix_size: if we have missed a comma or some terms, matrix_size will be too small */
				int size_check = 2*matrix_size*matrix_size;
				if (size_check < switch_terms)
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "\nError! Terms or comma missing from switch: " << next_switch << endl;

					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Error! Terms or comma missing from switch: " << next_switch << endl;
					}
					
if (debug_control::DEBUG >= debug_control::SUMMARY)	
    debug << "\nError! Terms or comma missing from switch: " << next_switch << endl;

					continue; // with the next switch
				}

				/* Check for the t variable */
				if (next_switch.find("[t]") != string::npos)
				{
					braid_control::T_VARIABLE = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: t-variable read from input switch " << endl;
				}
				
				
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "braid::main: matrix switch size = " << matrix_size << endl;

				Qpmatrix A(matrix_size,matrix_size);
				Qpmatrix B(matrix_size,matrix_size);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "braid::main: reading matrix A from next_switch:" << endl;

				/* read in A */				
				string::size_type posA = next_switch.find('=')+1;
				for (int i=0; i< matrix_size; i++)
				for (int j=0; j< matrix_size; j++)
				{
					string::size_type posB = next_switch.find(',',posA);
					string term = next_switch.substr(posA, (posB==string::npos? posB: posB-posA));

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "braid::main:  next term = " << term << endl;

					A[i][j] = Qpolynomial(term.substr(0,term.find('/')));

					string::size_type slash_pos = term.find('/');
					if (slash_pos != string::npos)
					{
						A[i][j] /= Qpolynomial(term.substr(slash_pos+1));
					}
					
					posA = posB+1;
				}
				
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "braid::main: A = \n" << A << endl;
	debug << "braid::main: reading matrix B from next_switch:" << endl;
}

				/* now B, note that B is set to the value provided, possibly multiplied through by t */				
				for (int i=0; i< matrix_size; i++)
				for (int j=0; j< matrix_size; j++)
				{
					string::size_type posB = next_switch.find(',',posA);
					string term = next_switch.substr(posA, (posB==string::npos? posB: posB-posA));

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "braid::main:  next term = " << term << endl;

					B[i][j] = Qpolynomial(term.substr(0,term.find('/')));
					
					if (braid_control::T_VARIABLE)
						B[i][j] *= Qpolynomial("t");

					string::size_type slash_pos = term.find('/');
					if (slash_pos != string::npos)
					{
						Qpolynomial d(term.substr(slash_pos+1));
						B[i][j] /= d;
					}
					
					posA = posB+1;
				}
				
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "braid::main: B = \n" << B << endl;
				

				/* Now work out 1-A */
				Qpmatrix A1 = A;
				
				for (int i=0; i< matrix_size; i++)
				{
					for (int j=0; j<i; j++)
						A1[i][j] *= Qpolynomial("-1");
						
					A1[i][i] = Qpolynomial("1") - A1[i][i];
					
					for (int j=i+1; j<matrix_size; j++)
						A1[i][j] *= Qpolynomial("-1");
				}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "braid::main: 1-A = \n" << A1 << endl;

				/* this is cheating really, it saves calling numcols etc 
				   these inverse values will be changed if A and B are non-singular
				*/
				Qpmatrix invA = A;
				Qpmatrix invB = B;

				if (determinant(A) != Qpolynomial("0") && 
				    determinant(B) != Qpolynomial("0") && 
					determinant(A1) != Qpolynomial("0"))
				{
					invA = A.inverse();
					invB = B.inverse();
					
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "braid::main: invA = \n" << invA << endl;
	Qpmatrix t1 = invA * A;
	debug << "braid::main: check = \n" << t1 << endl;
	debug << "braid::main: invB = \n" << invB << endl;
	Qpmatrix t2 = invB * B;
	debug << "braid::main: check = \n" << t2 << endl;

}


					/* Calculate (1-A)(A,B), where (A,B) = invA * invB * A * B */
					Qpmatrix AB = invA * invB * A * B;
					Qpmatrix A1AB = A1 * AB;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "braid::main: AB = " << AB << endl;
	debug << "braid::main: A1AB = " << A1AB << endl;
	debug << "braid::main: B*A1AB = " << B*A1AB << endl;
	debug << "braid::main: A1AB*B = " << A1AB*B << endl;
}
				
					/* Check the fundamental equation  [B,(1-A)(A,B)] = 0
					   Here, (X,Y) = invX * invY * X * Y and [X,Y] = X*Y - Y*X
					   we therefore check B * (1-A)(A,B) == (1-A)(A,B) * B, and we
					   have calculated A1AB = (1-A)(AB) above.
					*/
					
					if (B*A1AB - A1AB*B == A-A)
			    		fundamental_equation_satisfied = true;
					else
			    		fundamental_equation_satisfied = false;
				}
				else
				{
		    		fundamental_equation_satisfied = false;
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "braid::main: determinant(A) = " << determinant(A) << endl;
	debug << "braid::main: determinant(B) = " << determinant(B) << endl;
	debug << "braid::main: determinant(1-A) = " << determinant(A1) << endl;
}
		    	}

				if (fundamental_equation_satisfied)
				{

					// C = A^{-1}B^{-1}A(1-A) and
					if (!braid_control::SILENT_OPERATION)
						cout << "\ncalculating C" << endl;
				    Qpmatrix C = invA * invB * A * A1;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "braid::main: fundamenatal equation satisfied:" << endl;
	debug << "braid::main: C = \n" << C << endl;
}	
	
					if (!braid_control::SILENT_OPERATION)
						cout << "calculating D" << endl;
					Qpmatrix D = 1 - invA * invB * A * B;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "braid::main: D = \n" << D << endl;


if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "braid::YB-test ";
	debug << "A = (" << A[0][0] << ", " << A[0][1] << ", " <<  A[1][0] << ", " <<  A[1][1] << "), ";
	debug << "B = (" << B[0][0] << ", " << B[0][1] << ", " <<  B[1][0] << ", " <<  B[1][1] << "), ";
	debug << "C = (" << C[0][0] << ", " << C[0][1] << ", " <<  C[1][0] << ", " <<  C[1][1] << "), ";
	debug << "D = (" << D[0][0] << ", " << D[0][1] << ", " <<  D[1][0] << ", " <<  D[1][1] << ')' << endl;
}


					if (braid_control::EQUALITY_TEST)
					{
						equality_test(A,D,true);
						equality_test(B,C,false);
					}

					Qpmatrix switch_matrix(2*matrix_size,2*matrix_size);
					
					for (int i=0; i< matrix_size; i++)
					for (int j=0; j< matrix_size; j++)
						switch_matrix[i][j] = A[i][j];
				
					for (int i=0; i< matrix_size; i++)
					for (int j=0; j< matrix_size; j++)
						switch_matrix[i][j+matrix_size] = B[i][j];

					for (int i=0; i< matrix_size; i++)
					for (int j=0; j< matrix_size; j++)
						switch_matrix[i+matrix_size][j] = C[i][j];

					for (int i=0; i< matrix_size; i++)
					for (int j=0; j< matrix_size; j++)
						switch_matrix[i+matrix_size][j+matrix_size] = D[i][j];
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "braid::main: switch_matrix prior to cancellation = \n" << switch_matrix << endl;

					simplify_matrix(switch_matrix);


if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "braid::main: switch_matrix = \n" << switch_matrix << endl;


					if (!braid_control::SILENT_OPERATION)
						cout << "calculating switch matrix inverse" << endl;
					matrix_control::wait_threshold = 3;
					Qpmatrix switch_matrix_inverse = switch_matrix.inverse();

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "braid::main: inverse of S prior to cancellation = \n" << switch_matrix_inverse << endl;


					simplify_matrix(switch_matrix_inverse);

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "braid::main: switch_matrix_inverse = \n" << switch_matrix_inverse << endl;
	Qpmatrix t3 = switch_matrix_inverse * switch_matrix;
	debug << "braid::main: check = \n" << t3 << endl;
}

					if (switch_cache.size() && switch_cache[s].first != "")
					{
						if (!braid_control::SILENT_OPERATION)
							cout << "\n\n" << switch_cache[s].first << endl;
						if (!braid_control::RAW_OUTPUT)
							output << "\n\n" << switch_cache[s].first << endl;
					}
					else
						report_switch_matrix(switch_matrix, switch_matrix_inverse, "matrix");

					if (braid_control::SWITCH_POWER)
					{
						Qpmatrix power = switch_matrix;
						for (int i=1; i< braid_control::SWITCH_POWER; i++)
							power *= switch_matrix;
							
						if (!braid_control::SILENT_OPERATION)
							cout << "\nswitch matrix raised to the power " << braid_control::SWITCH_POWER << ":\n" << power << endl;
	
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "switch matrix raised to the power " << braid_control::SWITCH_POWER << ":";
							for (unsigned int i=0; i < power.numrows(); i++)
							{
								output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
								for (unsigned int j=0; j < power.numcols(); j++)
									output << power[i][j] << " ";
							}
						}
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::main: switch power = " << braid_control::SWITCH_POWER << "\n" << power << endl;
	
					}

					while (get_next_input_string(input_file,input_string, title))
						rat_poly_invariant(switch_matrix, switch_matrix_inverse, input_string, title);			
				}
				else
				{
					if (!braid_control::SILENT_OPERATION)
					{
						cout << "\nSwitch elements A = \n" << A << "\nB= \n" << B;
						cout << "\ndo not satisfy the unit and fundamental equation requirements .\n";
					}
	
					if (!braid_control::RAW_OUTPUT)
					{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "Switch elements:";
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "A = ";
							for (size_t i=0; i< A.numrows(); i++)
							{
								output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
								for (size_t j=0; j< A.numcols(); j++)
								{
									output << A[i][j] << " ";
								}
							}

							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "B = ";
							for (size_t i=0; i< B.numrows(); i++)
							{
								output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
								for (size_t j=0; j< B.numcols(); j++)
								{
									output << B[i][j] << " ";
								}
							}
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "\ndo not satisfy the unit and fundamental equation requirements .\n";
					}
		    		output.flush();

if (debug_control::DEBUG >= debug_control::BASIC)	
{
	debug << "braid::main: switch elements A = " << A << "\nB= " << B << endl;
	debug << "braid::main: do not satisfy the unit and fundamental equation requirements .\n" << endl;
}
				
				}
			}
			else if (braid_control::WEYL)
			{
							
				// these are controls and we may have more than one
				// switch in the input file, hence we need to set them all
				braid_control::CUSTOM_WEYL = false;
				braid_control::PRIME_WEYL = false;
				braid_control::QUANTUM_WEYL = false;
				braid_control::TRUNCATED_WEYL = false;
				braid_control::T_VARIABLE = false;
				braid_control::NUMERATOR_GCD = false;

				if (next_switch.find('U') != string::npos)
					braid_control::CUSTOM_WEYL = true;				
					
				if (next_switch.find('T') != string::npos)
					braid_control::TRUNCATED_WEYL = true; 
				else if (next_switch.find('P') != string::npos)
					braid_control::PRIME_WEYL = true;
				else if (next_switch.find('Q') != string::npos)
					braid_control::QUANTUM_WEYL = true; 

				/* Check for the t variable */
				if (next_switch.find("[t]") != string::npos)
				{
					braid_control::T_VARIABLE = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: t-variable read from input switch " << endl;
				}

					
				/* Determine the values of n and, where necessary, p*/
				int n;
				int p;

				if (braid_control::CUSTOM_WEYL)
				{
					//  first determine the value p 
					char* str = c_string(next_switch);
					char* c1 = str;
					while (!isdigit(*c1) && *c1 != '\0')
						c1++;

					if (*c1 != '\0')
					{
						// c1 points at the start of the first digit 
						get_number(p, c1);
					}
					
					delete [] str;

					// now n, which we determine from the number of terms
					int switch_terms = -1;
					string::size_type pos = 0;
					while ((pos = next_switch.find(',',pos)) != string::npos)
					{
						switch_terms++;
						pos++;
					}

					n = static_cast<int>(sqrt(float(switch_terms/2)));
					

					/* test n: if we have missed a comma or some terms, n will be too small */
					int size_check = 2*n*n;
					if (size_check < switch_terms)
					{
						if (!braid_control::SILENT_OPERATION)
							cout << "\nError! Terms or comma missing from switch: " << next_switch << endl;

						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "Error! Terms or comma missing from switch: " << next_switch << endl;							
						}
					
if (debug_control::DEBUG >= debug_control::SUMMARY)	
    debug << "\nError! Terms or comma missing from switch: " << next_switch << endl;

						continue; // with the next switch
					}
					
				}
				else
				{	
					char* str = c_string(next_switch);
					char* c1 = str;
					while (!isdigit(*c1) && *c1 != '\0')
						c1++;
				
					if (*c1 != '\0')
					{
				
						/* c1 points at the start of the first digit */
						get_number(n, c1);
	
						while (*c1 != ',' && *c1 != '\0')
							c1++;

						if (*c1 != '\0')
						{
							c1++; //step over the comma
							while (!isdigit(*c1) && *c1 != '\0')
								c1++;						
						}

						if (*c1 != '\0')
						{
							get_number(p, c1);		
						}
						else if (braid_control::QUANTUM_WEYL)
						{
							p=0; // default is integers
						}
						else
						{
								// find the first prime dividing n
								p = 2;
								while (n%p)
									p++;
						}

					}
					else 
					{
						/* set default values if none are provided */
						if (braid_control::QUANTUM_WEYL)
						{
							n = 2;
							p = 0;
						}
						else
							n = p = 2;							
					}
				
					delete[] str;
				}
				
				/* set the scalar variant */				
				if (p) 
				{
					scalar::set_variant(scalar::MOD_P);	
					polynomial_control::MOD_P = true;
					mod_p::set_p(p);
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "braid::main: Weyl algebra switch requires MOD_P scalars with p = " << mod_p::get_p() << endl;
	debug << "braid::main: setting polynomial_control::MOD_P = true" << endl;
}
				}
				else if (!braid_control::CALCULATE_MOD_P)
				{
					scalar::set_variant(scalar::BIGINT); // here p==0 but we may have been given the global mod-p option
					polynomial_control::MOD_P = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "braid::main: Weyl algebra switch requires BIGINT scalars" << endl;
	debug << "braid::main: setting polynomial_control::MOD_P = false" << endl;
}
				}
				else
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::main: Weyl algebra switch using globally specified MOD_P scalars with p = " << mod_p::get_p() << endl;
				}

				typedef Qpolynomial Wpolynomial;
				typedef matrix<Wpolynomial,scalar> Wmatrix;
				Wmatrix B(n,n);
				Wmatrix V(n,n);

				bool valid_Weyl_algebra = true;  // check only needed for custom algebras
				
				if (braid_control::CUSTOM_WEYL)
				{
					if (next_switch.find('G') != string::npos)
						braid_control::NUMERATOR_GCD = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "braid::main: custom Weyl algebra switch, n = " << n << ", p = " << p << endl;
	if (braid_control::QUANTUM_WEYL)
		debug << "braid::main: custom version of quantum Weyl algebra switch" << endl;
	if (braid_control::NUMERATOR_GCD)
		debug << "braid::main: calculating gcd of Delta_1 generators' numerator" << endl;
}

					/* read in B (= U) note that B is set to the value provided, possibly multiplied through by t */				
					string::size_type posA = next_switch.find(','); // S = WUV,p,...
					posA = next_switch.find(',',posA+1)+1; // moves over p
					for (int i=0; i< n; i++)
					for (int j=0; j< n; j++)
					{
						string::size_type posB = next_switch.find(',',posA);
						string term = next_switch.substr(posA, (posB==string::npos? posB: posB-posA));

						B[i][j] = Wpolynomial(term.substr(0,term.find('/')));

						if (braid_control::T_VARIABLE)
							B[i][j] *= Wpolynomial("t");

						string::size_type slash_pos = term.find('/');
						if (slash_pos != string::npos)
						{
							B[i][j] /= Wpolynomial(term.substr(slash_pos+1));
						}

						posA = posB+1;
					}

					/* now V, */
					for (int i=0; i< n; i++)
					for (int j=0; j< n; j++)
					{
						string::size_type posB = next_switch.find(',',posA);
						string term = next_switch.substr(posA, (posB==string::npos? posB: posB-posA));

						V[i][j] = Wpolynomial(term.substr(0,term.find('/')));


						string::size_type slash_pos = term.find('/');
						if (slash_pos != string::npos)
						{
							V[i][j] /= Wpolynomial(term.substr(slash_pos+1));
						}

						posA = posB+1;
					}

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "braid::main: B (= U) = \n" << B << endl;
	debug << "braid::main: V = \n" << V << endl;
}					
					/* Check that U and V define a valid Weyl algebra.  For a custom version
					   of the quantum Weyl algebra we require that UV-qVU = I and for the others
					   that UV-VU = I
					*/
					Wmatrix	identity(n,n);
					for (int i=0; i< n; i++)
					{
						for (int j=0; j< i; j++)
							identity[i][j] = Wpolynomial("0");

						identity[i][i] = Wpolynomial("1");

						for (int j=i+1; j< n; j++)
							identity[i][j] = Wpolynomial("0");
					}

					Wmatrix weyl_check = V*B;
						
					if (braid_control::QUANTUM_WEYL)
					{
						Wpolynomial q("q");
						for (int i=0; i< n; i++)
						for (int j=0; j< n; j++)
							weyl_check[i][j] *= q;
								
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "braid::main: Weyl algebra check: B*V - qV*B = " << B*V - weyl_check << endl;
					}
					else
					{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "braid::main: Weyl algebra check: B*V - V*B = " << B*V - weyl_check << endl;
					}
						
					if (B*V - weyl_check != identity)
						valid_Weyl_algebra = false;
				}
				else // standard Weyl algebras
				{
					valid_Weyl_algebra = true;
					
					if (braid_control::PRIME_WEYL)
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::main: prime Weyl algebra switch with n = " << n << ", p = " << p << endl;

						/* Assign B (= U) */
						char v = 'a';
						for (int i = 0; i< n; i++)
						{
							B[i][i] = Wpolynomial("1");
							if (i > 0)
							{
								string s(1,v++);
								B[i-1][i] = Wpolynomial(s);
							}
						}

						/* Assign V */
						v = 'a';
						for (int i = 0; i< n; i++)
						{
							V[i][i] = Wpolynomial("1");
							if (i > 0)
							{
								if (i % p)
								{
									ostringstream oss;
									oss << i;
									V[i][i-1] = Wpolynomial(oss.str());
									V[i][i-1] /= Wpolynomial(string(1,v++));
								}
								else
									v++;
							}
						}
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "braid::main: B (= U) = \n" << B << endl;
	debug << "braid::main: V = \n" << V << endl;
}					

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::main: Weyl algebra check: B*V - V*B = " << B*V - V*B << endl;
					}
					else if (braid_control::TRUNCATED_WEYL)
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::main: truncated Weyl algebra switch with n = " << n << ", p = " << p << endl;

						char var = 'a';
						char Var = 'A';

						vector<Wpolynomial> I_coeff(n, Wpolynomial("0"));
						vector<Wpolynomial> J_coeff(n, Wpolynomial("0"));
						vector<Wpolynomial> K_coeff(n, Wpolynomial("0"));

						for (int i=0; i< n; i++)
						{
							string s(1,var++);
							string S(1,Var++);

							I_coeff[i] = Wpolynomial(s);
							J_coeff[i] = Wpolynomial(S);
						}
						
						
						/* Calculate K_coeff[r] for r = 1,..., n-1 using the difference equation based on
						
						   I'      = i_1 + 2i_{2}x + ... + (n-1)i_{n-1}x^{n-2}
						   I'^{-1) = k_0 + k_{1}x + ... + k_{n-1}x^{n-1}
						   
						   multiplying out and equating coefficients (since I' * I'^{-1} = 1).  Thus the coefficient of
						   X^0 is 1 and the coefficient of x^r=0 for r=1,...n-1.
						   
						   We have to do this in two stages, since the coefficient of x^{n-1} does not involve
						   as many terms as the other coefficients. This is because when collecting terms in x^{n-1} in the product 
						   I' * I'^{-1} we can can only go up to (n-1)i_{n-1}x^{n-2} * k_{1}x
						
						
						   Note also that when evaluating K_coeff[r] from the collected terms, in each case there is
						   a factor i_1 involved, since the first contribution to the coefficients of x^r is the term
						   i_1k_r.  Therefore we sum from the 2i_{2}x * k_{r-1}x^{r-1} and then multiply through by
						   -i_1^{-1}; we refer to this multiplicative term as the 'i1factor'.
						*/
						
						
						K_coeff[0] = Wpolynomial("1"); 
						K_coeff[0] /= Wpolynomial("b"); // i_1^{-1}
						Wpolynomial i1factor = Wpolynomial("-1") * K_coeff[0];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::main: i1 factor for K_coeffs = " << i1factor << endl;					

						for (int r=1; r < n-1; r++)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::main: r=" << r << endl;

							Wpolynomial sum = Wpolynomial("0");
							for (int k = 1; k <= r; k++)
							{
								ostringstream oss;
								oss << k+1;
								Wpolynomial term  = Wpolynomial(oss.str());
								term *= I_coeff[k+1] * K_coeff[r-k];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::main: \t" << k+1 << "i" << k+1 << "k" << r-k << endl;							
								sum += term;							
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::main: \t\tterm = " << term << "\tsum = " << sum << endl;
							}

							K_coeff[r] = i1factor * sum;
						}

						/* Now calculate K_coeff[n-1] */
						{ // new scope to make sum a local variable again
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::main: r=" << n-1 << endl;

							Wpolynomial sum = Wpolynomial("0");
							for (int k = 1; k <= n-2; k++)
							{
								ostringstream oss;
								oss << k+1;
								Wpolynomial term  = Wpolynomial(oss.str());
								term *= I_coeff[k+1] * K_coeff[n-1-k];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::main: \t" << k+1 << "i" << k+1 << "k" << n-1-k << endl;							
								sum += term;							
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::main: \t\tterm = " << term << "\tsum = " << sum << endl;
							}
							K_coeff[n-1] = i1factor * sum;
						}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "braid::main: I coefficients: ";
	for (int i=0; i < n; i++)
		debug << I_coeff[i] << ' ';
	debug << endl;
	debug << "braid::main: J coefficients: ";
	for (int i=0; i < n; i++)
		debug << J_coeff[i] << ' ';
	debug << endl;
	debug << "braid::main: K coefficients: ";
	for (int i=0; i < n; i++)
		debug << K_coeff[i] << ' ';
	debug << endl;
	
	/* If the K coefficients have been calculated correctly, multiplying out K by I
	   should give 1
	*/
	debug << "braid::main: check that I*K = 1; I*K = ";
	Wpolynomial x("x");
	Wpolynomial xpower("x");
	Wpolynomial Ipoly = I_coeff[1];
	Wpolynomial Kpoly = K_coeff[0];

	for (int i=1; i< n-1; i++)
	{
		ostringstream oss;
		oss << i+1;
		Wpolynomial t(oss.str());
		Ipoly += t*I_coeff[i+1]*xpower;
		Kpoly += K_coeff[i]*xpower;
		xpower *= x;
	}

	Kpoly += K_coeff[n-1]*xpower;
	debug << Ipoly * Kpoly << endl;
}

						/* The Fenn-Turaev paper, Weyl algebras and knots gives the matrix
						   representation of the truncated Weyl algebra with u and v interchanged,
						   so that uv-vu = -1.  Here we assign them correctly so that BV-VB = 1
						*/
						   
						/* Assign B (= U) */
						for (int c=0; c < n; c++)
						{					
							for (int r=0; r < n; r++)
							{
								if (c >=r)
								{
									B[r][c] = I_coeff[c-r];
								}
							}
						}

						/* Assign V */
						for (int r=0; r < n; r++)
						{					
							Wpolynomial t = Wpolynomial("0");
							if (r % p) // the k factor is non-zero mod p
							{
								ostringstream oss;
								oss << r;
								t = Wpolynomial(oss.str());
							}

							for (int c=0; c < n; c++)
							{
								/* the index of j is c-r when c >= r*/
								if (c >= r)
								{
									V[r][c] = J_coeff[c-r];
								}

								/* the index of k is c-r+1 when c >= r-1 
								   note the condition that r%p means we don't 
								   have to worry about row zero.
								*/
								if (r % p && c >= r-1)
								{
									V[r][c] += t * K_coeff[c-r+1];
								}
							}
						}
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "braid::main: B (= U) = \n" << B << endl;
	debug << "braid::main: V = \n" << V << endl;
}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::main: Weyl algebra check: B*V - V*B = " << B*V - V*B << endl;
					}
					else if (braid_control::QUANTUM_WEYL)
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::main: quantum Weyl algebra switch with n = " << n << ", p = " << p << endl;
						/* Assign B (= U) */
						for (int i = 0; i< n; i++)
						{
							ostringstream s1;
							if (i<n-2)
								s1 << "q^" << n-1-i;
							else if (i == n-2)
								s1 << "q";
						
							s1 << "a";				
							B[i][i] = Wpolynomial(s1.str());
											
							ostringstream s2;					
						
							if (i > 0)
							{
								if (i<n-2)
									s2 << "b^" << n-1-i;
								else if (i == n-2)
									s2 << "b";
	
								s2 <<"d";				

								B[i-1][i] = Wpolynomial(s2.str());
							}
						}

						/* Assign V */
						for (int i = 0; i< n; i++)
						{
							ostringstream s1;

							if (i==0)
								s1 << "1";
							else if (i==1)
								s1 << "q";
							else if (i>1)
								s1 << "q^" <<i;

							V[i][i] = Wpolynomial(s1.str());

							// divide V[i][i] by aq^{n-1}-aq^n)
							ostringstream s2;
							if (n == 2)
								s2 << "aq-aq^" << n;
							else
								s2 << "aq^" << n-1 << "-aq^" << n;				

							V[i][i] /= Wpolynomial(s2.str());

											
							ostringstream s3;
						
							if (i > 0)
							{
								if (i==2)
									s3 << "q";
								else if (i>2)
									s3 << "q^" << i-1;
	
								s3 << "e";								
							
								V[i-1][i] = Wpolynomial(s3.str());
							
								ostringstream s4;
							
								if (i==2)
									s4 << "b";
								else if (i>2)
									s4 << "b^" << i-1;
							
								if (i>1)
									V[i-1][i] /= Wpolynomial(s4.str());								
							}
						}

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "braid::main: B (= U) = \n" << B << endl;
	debug << "braid::main: V = \n" << V << endl;
}					

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	Wmatrix weyl_check = V*B;					
	Wpolynomial q("q");
	for (int i=0; i< n; i++)
	for (int j=0; j< n; j++)
		weyl_check[i][j] *= q;
								
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::main: Weyl algebra check: B*V - qV*B = " << B*V - weyl_check << endl;
}
					}
				}
				
				if (valid_Weyl_algebra)
				{
					Wmatrix invB = B.inverse();
					Wmatrix A = V.inverse() * invB;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "braid::main: V.inverse() = \n" << V.inverse() << endl;
	debug << "braid::main: check = " << V*V.inverse() << endl;
	debug << "braid::main: B.inverse() = \n" << invB << endl;
	debug << "braid::main: check = " << B*invB << endl;
	debug << "braid::main: A = \n" << A << endl;
}					

					/* Now work out 1-A */
					Wmatrix A1 = A;

					for (int i=0; i< n; i++)
					{
						for (int j=0; j<i; j++)
							A1[i][j] *= Wpolynomial("-1");

						A1[i][i] = Wpolynomial("1") - A1[i][i];

						for (int j=i+1; j<n; j++)
							A1[i][j] *= Wpolynomial("-1");
					}

					/* since V has been scaled, the inverse we calculate here is
					   also a scaled form of the inverse; the scale factor is
					   V_scale_factor^n
					*/
					Wmatrix invA=B*V;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "braid::main: 1-A = \n" << A1 << endl;
	debug << "braid::main: Ainv () = U*V = \n" << invA << endl;
	debug << "braid::main: check = \n"  << A*invA << endl;
}

					/* We know the fundamental equation is satisfied, 
					   so go ahead and calculate C and D
					*/

					// C = A^{-1}B^{-1}A(1-A)
			    	Wmatrix C = invA * invB * A * A1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "braid::main: C = \n" << C << endl;

					Wmatrix D = 1 - invA * invB * A * B;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "braid::main: D = \n" << D << endl;

					if (braid_control::EQUALITY_TEST)
					{
						equality_test(A,D,true);
						equality_test(B,C,false);
					}

					Wmatrix switch_matrix(2*n,2*n);

					for (int i=0; i< n; i++)
					for (int j=0; j< n; j++)
						switch_matrix[i][j] = A[i][j];

					for (int i=0; i< n; i++)
					for (int j=0; j< n; j++)
						switch_matrix[i][j+n] = B[i][j];

					for (int i=0; i< n; i++)
					for (int j=0; j< n; j++)
						switch_matrix[i+n][j] = C[i][j];

					for (int i=0; i< n; i++)
					for (int j=0; j< n; j++)
						switch_matrix[i+n][j+n] = D[i][j];

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "braid::main: switch_matrix  = \n" << switch_matrix << endl;


					/* If we're dealing with the quantum Weyl algebra qS' = S +(q-1), so we just take the 
					   inverse as usual.  For other cases, S^2 = I, so the inverse is just S */
					Wmatrix switch_matrix_inverse = switch_matrix;

					if (braid_control::QUANTUM_WEYL)
						switch_matrix_inverse = switch_matrix.inverse();

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "braid::main: Weyl switch matrix inverse:\n" << switch_matrix_inverse << endl;
	debug << "braid::main: Weyl switch matrix inverse check:\n" << switch_matrix * switch_matrix_inverse << endl;
	
	if (braid_control::QUANTUM_WEYL)
	{
		Wmatrix Ssquared = switch_matrix*switch_matrix;
		Wmatrix q1S = switch_matrix;
		Wpolynomial q1("q-1");
		for (int i=0; i< 2*n; i++)
		for (int j=0; j< 2*n; j++)
			q1S[i][j] *= q1;
		debug << "braid::main: custom quantum Weyl switch matrix check that S^2 + (q-1)S = q:  S^2 + (q-1)S = \n" 
		      << Ssquared + q1S << endl;
	}
}

					if (switch_cache.size() && switch_cache[s].first != "")
					{
						if (!braid_control::SILENT_OPERATION)
							cout << "\n\n" << switch_cache[s].first << endl;
						if (!braid_control::RAW_OUTPUT)
							output << "\n\n" << switch_cache[s].first << endl;
					}
					else
						report_switch_matrix(switch_matrix, switch_matrix_inverse, (braid_control::QUANTUM_WEYL? "quantum Weyl":"Weyl"));

					if (braid_control::SWITCH_POWER)
					{
						Wmatrix power = switch_matrix;
						for (int i=1; i< braid_control::SWITCH_POWER; i++)
							power *= switch_matrix;
							
						if (!braid_control::SILENT_OPERATION)
							cout << "\nswitch matrix raised to the power " << braid_control::SWITCH_POWER << ":\n" << power << endl;
	
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "switch matrix raised to the power " << braid_control::SWITCH_POWER << ":";
							for (unsigned int i=0; i < power.numrows(); i++)
							{
								output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
								for (unsigned int j=0; j < power.numcols(); j++)
									output << power[i][j] << " ";
							}
						}
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::main: switch power = " << braid_control::SWITCH_POWER << "\n" << power << endl;
	
					}

					while (get_next_input_string(input_file,input_string, title))
						rat_poly_invariant(switch_matrix, switch_matrix_inverse, input_string, title);			
				}
				else
				{
					if (!braid_control::SILENT_OPERATION)
					{
						cout << "\nMatrices\nU = \n" << B << "\nV= \n" << V;
						cout << "\ndo not define a Weyl algebra.\n";
					}

					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Matrices\n";
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "U = ";
						for (size_t i=0; i< B.numrows(); i++)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							for (size_t j=0; j< B.numcols(); j++)
							{
								output << B[i][j] << " ";
							}
						}

						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "V = ";
						for (size_t i=0; i< V.numrows(); i++)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							for (size_t j=0; j< V.numcols(); j++)
							{
								output << V[i][j] << " ";
							}
						}
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "\ndo not define a Weyl algebra.\n";
					}
    				output.flush();

if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "braid::main: matrices\nU = \n" << B << "\nV= \n" << V << endl;
	debug << "braid::main: do not define a Weyl algebra.\n";
}
				
				}
			}
			else if (braid_control::FIXED_POINT_INVARIANT || braid_control::RACK_POLYNOMIAL)
			{
				/*  determine the value of n, whether the switch is numbered from 
				    zero and how many terms we've been given
				*/
				int size=0;
				int matrix_element;
				int switch_terms=0;
				char* str = c_string(next_switch);
				char* cptr = str;
				
				/* we assume the switch is not numbered from zero and correct if
				   a zero matrix element is found
				*/

if (debug_control::DEBUG >= debug_control::DETAIL)	
	debug << "braid::main: fixed-point switch matrix elements: ";
	
				bool number_from_zero = false;
				while (*cptr != '\0')
				{
					/* look for start of next digit */
					while (!isdigit(*cptr) && *cptr != '\0')
						cptr++;
					
					if (*cptr != '\0')
					{
						get_number(matrix_element,cptr);
						switch_terms++;

if (debug_control::DEBUG >= debug_control::DETAIL)	
	debug << matrix_element << " ";
					
						if (matrix_element > size)
							size = matrix_element;
						
						if (matrix_element == 0)
							number_from_zero = true;
					
						/* skip over the number */
						while (isdigit(*cptr) && *cptr != '\0')
							cptr++;
					}
				}

if (debug_control::DEBUG >= debug_control::DETAIL)	
	debug << endl;
				
				if (number_from_zero)
					size++;
				
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
{
	debug << "braid::main: fixed-point switch size = " << size << endl;
	debug << "braid::main: fixed-point number of switch terms detected = " << switch_terms << endl;
}

				/* do we have a twist matrix specified */
				bool default_twist = (strchr(str,'T')?false:true);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "braid::main: fixed-point switch default_twist = " << (default_twist?"true":"false") << endl;
				
				/* Is this switch claimed to be an essential welded pair, suitable for use with welded braids? */
				bool essential_welded_pair = (strchr(str,'W')?true:false);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "braid::main: fixed-point switch essential_welded_pair = " << (essential_welded_pair?"true":"false") << endl;

				/* Is this switch claimed to be an essential doodle pair, suitable for use with doodles? */
				bool essential_doodle_pair = (strchr(str,'D')?true:false);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "braid::main: fixed-point switch essential_doodle_pair = " << (essential_doodle_pair?"true":"false") << endl;
	
				/* If we are to evaluate welded invariants, we have to have a non-default twist */
				if (essential_welded_pair && default_twist)
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
    debug << "\nError! Welded invariants requested with the default twist: " << next_switch << endl;
					continue;
				}
		
				/* test the number of switch_terms provided against the detected size of X_n */
				if (switch_terms != (default_twist? 2*size*size: 4*size*size))
				{
					if (!braid_control::SILENT_OPERATION)
					{
						cout << "\nError! Incorrect number of terms in switch: " << next_switch
							 << " found " << switch_terms << " terms, expected " <<
							(default_twist? 2*size*size: 4*size*size) << endl;
					}

					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Error! Incorrect number of terms in switch: " << next_switch << endl;
					}
					
if (debug_control::DEBUG >= debug_control::SUMMARY)	
    debug << "\nError! Incorrect number of terms in switch: " << next_switch << endl;

					continue; // with the next switch
				}

				matrix<int> Su(size,size), Sd(size,size), Tu(size,size), Td(size,size);
				cptr = str;
				
				/* read in the matrix values, we have e.g. S = FT U(0,0) U(0,1) ... 
				   we know there are the correct numnber of terms present, so we don't 
				   need to worry about hitting a terminating null prematurely
				*/
				cptr = str;
				
				/* read in Su and Sd */
				for (int i=0; i< size; i++)
				for (int j=0; j< size; j++)
				{
					while (!isdigit(*cptr))
						cptr++;
					
					get_number(Su[i][j], cptr);

					if (!number_from_zero)
						Su[i][j]--;

					while (isdigit(*cptr))
						cptr++;
					
				}
				
				for (int i=0; i< size; i++)
				for (int j=0; j< size; j++)
				{
					while (!isdigit(*cptr))
						cptr++;
					
					get_number(Sd[i][j], cptr);

					if (!number_from_zero)
						Sd[i][j]--;
					
					while (isdigit(*cptr))
						cptr++;
				}
				
				/* If there's a non-default twist, read it in, otherwise assign the default twist
				   T(x,y)=(y^x,x_y)=(y,x)  In the case of the default twist Tu has (x,y) entry 
				   y^x=y and Td has (y,x) entry x_y=x; thus Td has (x,y) entry y_x=y and so 
				   both Tu and Td have all their rows 0...(n-1).
				*/
				for (int i=0; i< size; i++)
				for (int j=0; j< size; j++)
				{
					if (!default_twist)
					{
						while (!isdigit(*cptr))
							cptr++;
					
						get_number(Tu[i][j], cptr);

						if (!number_from_zero)
							Tu[i][j]--;

						while (isdigit(*cptr))
							cptr++;
					}
					else
						Tu[i][j] = j;
					
				}
				
				for (int i=0; i< size; i++)
				for (int j=0; j< size; j++)
				{
					if (!default_twist)
					{
						while (!isdigit(*cptr))
							cptr++;
					
						get_number(Td[i][j], cptr);

						if (!number_from_zero)
							Td[i][j]--;

						while (isdigit(*cptr))
							cptr++;
					}
					else
						Td[i][j] = j;
					
				}

				delete [] str;
				
if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	bool loc_newline = matrix_control::SINGLE_LINE_OUTPUT;
	matrix_control::SINGLE_LINE_OUTPUT = true;
	
	debug << "braid::main: fixed-point switch Su = " << Su << endl;
	debug << "braid::main: fixed-point switch Sd = " << Sd << endl;
	debug << "braid::main: fixed-point switch Tu = " << Tu << endl;
	debug << "braid::main: fixed-point switch Td = " << Td << endl;
	
	matrix_control::SINGLE_LINE_OUTPUT = loc_newline;
}

				/* Check the validity of S and (if necessary) T: both S and T must be permutations 
				   of X x X, they must have rows that are permutations, T must have order 2 and satisfy 
				   S-Yang-Baxter.  The switch S may be required to satisfy, or not to satisfy S-Yang-Baxter,
				   dependent on whether we are working with doodle invariants, so this test is left for later.
				   Similarly, S may be required to have order 2 (to be an essential doodle pair or a flat essential
				   virtual pair)				   
				*/

				if (!XxX_permutation (Su, Sd))
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch S is not a permutation of XxX" << endl;
					continue;
				}	
			
				/* The rows of Su, Sd, Tu and Td must be permutations */
				if (!rows_are_permutations(Su))
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch matrix Su has rows that are not permutations" << endl;
					continue;
				}
				
				if (!rows_are_permutations(Sd))
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch matrix Sd has rows that are not permutations" << endl;
					continue;
				}
				
				/* Now do the same tests for T if it is not the default twist, and also check that T^2 = Id and
				   that T one-dominates S if we've been given a non-default twist (recall that if T has order 2 
				   T > S <=> T >> S).
				*/
				if (!default_twist)
				{
					if (!XxX_permutation (Tu, Td))
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch T is not a permutation of XxX" << endl;
						continue;
					}

					if (!rows_are_permutations(Tu))
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch matrix Tu has rows that are not permutations" << endl;
						continue;
					}
				
					if (!rows_are_permutations(Td))
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch matrix Td has rows that are not permutations" << endl;
						continue;
					}
				
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: checking the Yang-Baxter relation for fixed-point switch T" << endl;

					if (!Yang_Baxter_satisfied(Tu, Td))
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch T does not satisfy the Yang-Baxter relation" << endl;
						continue;
					}		

					if (switch_order(Tu,Td) == 2)
					{						
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch T has order 2" << endl;
	
						if (one_dominates(Tu,Td,Su,Sd))
						{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch T one-dominates S" << endl;
						}
						else
						{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch T does not one-dominate S" << endl;
							continue;
						}
					}
					else					
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch T does not have order 2" << endl;
						continue;				
					}   
				}

				/* Now verify that we have an essential virtual pair, essential welded pair or an essential doodle pair,
				   as specified in the switch definition. 
				   
				   If S and T are claimed to be an essential doodle pair then S is requried not to satisfy S-Yang-Baxter,
				   otherwise it must satisfy S-Yang-Baxter.  Moreover, if S and T are claimed to be an essential doodle 
				   pair or an essential virtual pair then S must not one-dominate, nor two-dominate T. If S and T are 
				   claimed to be an essential welded pair, then S must one-dominate T but fail to two-dominate T.
				   
				   A finite switch is explicitly indicated as an essential welded pair or an essential doodle pair, if neither of 
				   these indications are present it is considered to be an essential virtual pair.
				*/
				bool S_of_order_2;
				if (switch_order(Su, Sd) == 2)
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch S has order 2" << endl;
						S_of_order_2 = true;
				}
				else
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch S does not have order 2" << endl;
						S_of_order_2 = false;
				}

if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: checking the Yang-Baxter relation for fixed-point switch S" << endl;

				bool S_Yang_Baxter;
				if (Yang_Baxter_satisfied(Su, Sd))
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch S satisfies the Yang-Baxter relation" << endl;
						S_Yang_Baxter = true;
				}
				else
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch S does not satisfy the Yang-Baxter relation" << endl;
						S_Yang_Baxter = false;
				}
				
				bool S_one_dominates_T;
				if (one_dominates(Su,Sd,Tu,Td))
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch S one-dominates T" << endl;
					S_one_dominates_T = true;
				}
				else
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch S does not one-dominate T" << endl;
					S_one_dominates_T = false;

				}
				
				bool S_two_dominates_T;
				if (two_dominates(Su,Sd,Tu,Td))
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch S two-dominates T" << endl;
					S_two_dominates_T = true;
				}
				else
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switch S does not two-dominate T" << endl;
					S_two_dominates_T = false;

				}

				braid_control::ST_pair_type pair_type;
				if (essential_doodle_pair)
				{
					if (!S_of_order_2 || S_Yang_Baxter || S_one_dominates_T || S_two_dominates_T)
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switches S and T do not satisfy the conditions for an essential doodle pair" << endl;
						continue; // to the next switch
					}
					else
					{
						pair_type = braid_control::ST_pair_type::ESSENTIAL_DOODLE;
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switches S and T validated as an essential doodle pair" << endl;
					}
				}
				else if (essential_welded_pair)
				{
					if (!S_Yang_Baxter || !S_one_dominates_T || S_two_dominates_T)
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switches S and T do not satisfy the conditions for an essential welded pair" << endl;
						continue; // to the next switch
					}
					else
					{
						pair_type = braid_control::ST_pair_type::ESSENTIAL_WELDED;
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switches S and T validated as an essential welded pair" << endl;
					}
				}
				else // essential virtual pair
				{
//					if (!S_Yang_Baxter || S_one_dominates_T || S_two_dominates_T)
					if (!S_Yang_Baxter)
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switchs S does not satisfy Yang Baxter equations" << endl;
						continue; // to the next switch
					}
					else if (!braid_control::CLASSICAL_ONLY && (S_one_dominates_T || S_two_dominates_T))
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switches S and T do not satisfy the conditions for an essential virtual pair" << endl;
						continue; // to the next switch
					}
					else
					{
						if(S_of_order_2)
						{
							pair_type = braid_control::ST_pair_type::FLAT_ESSENTIAL_VIRTUAL;
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switches S and T validated as a flat essential virtual pair" << endl;
						}
						else
						{
							pair_type = braid_control::ST_pair_type::ESSENTIAL_VIRTUAL;
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "braid::main: fixed-point switches S and T validated as an essential virtual pair" << endl;
						}
					}
				}
				

				/* calculate the inverse of S, invS(x,y) = (y_{\bar{x}},x^{\bar{y}}); we create 
				   matrices invSu and invSd so that the (x,y) entry of invSu is y^{\bar x} and 
				   the (x,y) entry of invSd is y_{\bar x}.
	   
				   From the Reidemeister II move with a positive crossing followed by a negative crossing
				   we have that (b^a)^{\bar{a_b}}=b and (a_b)_{\bar{b^a}} = a, so that the (a_b, b^a) entry
				   of invSu is b and the (b^a, a_b) entry of invSd is a.
				*/
		
				matrix<int> invSu(size,size), invSd(size,size);
				
				for (int a=0; a< size; a++)
				for (int b=0; b< size; b++)
				{
					int b_up_a = Su[a][b];
					int a_down_b = Sd[b][a];
						
					invSu[a_down_b][b_up_a]=b;
					invSd[b_up_a][a_down_b]=a;		
				}


if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	bool saved_bool = matrix_control::SINGLE_LINE_OUTPUT;
	matrix_control::SINGLE_LINE_OUTPUT = true;
	debug << "braid::main: inverse of S: invSu = " << invSu << " invSd = " << invSd << endl;
	matrix_control::SINGLE_LINE_OUTPUT = saved_bool;
}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "braid::main: fixed-point switch S has order " << switch_order(Su,Sd) << endl;
	debug << "braid::main: fixed-point switch invS has order " << switch_order(invSu,invSd) << endl;

	matrix<int> checka(size,size), checkb(size,size);
	
	/* we should have invS(S(a,b)) = (a,b) that is invS(b^a,a_b) = (a,b) and we have 
	   invS(x,y) = (y_{\bar x}, x^{\bar y}), so invS(b^a,a_b) = ((a_b)_{\bar b^a}, (b^a)^{\bar a_b})
	   
	   Thus we set checka[a][b] to be (a_b)_{\bar b^a} and checkb[a][b] to be (b^a)^{\bar a_b}, 
	   so that checka[a][b] should be a and checkb[a][b] should be b, that is row i of checka 
	   should be all i and column j of checkb should all be j   
	*/

	for (int a=0; a< size; a++)
	for (int b=0; b< size; b++)
	{
		int  b_up_a = Su[a][b];
		int a_down_b = Sd[b][a];
		checka[a][b] = invSd[b_up_a][a_down_b];
		checkb[a][b] = invSu[a_down_b][b_up_a];
	}
	debug << "braid::main: fixed-point switch checka = " << checka; 
	debug << "braid::main: fixed-point switch checkb = " << checkb;
}

				/* report the switch, we do this explicitly as report_switch is designed for linear switch matrices */

				if (!braid_control::SILENT_OPERATION)
				{
					if (essential_welded_pair)
						cout << "\nEssential welded pair";
					else if (essential_doodle_pair)
						cout << "\nEssential doodle pair";
					else if (!braid_control::CLASSICAL_ONLY)
						cout << "\nEssential virtual pair";
				}

				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");

					if (essential_welded_pair)
						output << "Essential welded pair";
					else if (essential_doodle_pair)
						output << "Essential doodle pair";
					else if (!braid_control::CLASSICAL_ONLY)
						output << "Essential virtual pair";
				}
				
				if (switch_cache.size() && switch_cache[s].first != "")
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "\n" << switch_cache[s].first;
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << switch_cache[s].first;
					}
				}

				bool saved_bool = matrix_control::SINGLE_LINE_OUTPUT;
				matrix_control::SINGLE_LINE_OUTPUT = true;
				if (!braid_control::SILENT_OPERATION)
				{
					cout << "\nSu = ";
					display_fixed_point_switch(Su, cout, number_from_zero);
					cout << "  Sd = ";
					display_fixed_point_switch(Sd, cout, number_from_zero);
					cout << "\nTu = ";
					display_fixed_point_switch(Tu, cout, number_from_zero);
					cout << "  Td = ";
					display_fixed_point_switch(Td, cout, number_from_zero);
					cout << endl;
					cout << "\nInverse invSu = ";
					display_fixed_point_switch(invSu, cout, number_from_zero);
					cout << "  invSd = ";
					display_fixed_point_switch(invSd, cout, number_from_zero);
					cout << "\n" << endl;
				}
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");					
					output << "Su = ";
					display_fixed_point_switch(Su, output, number_from_zero);
					output << " Sd = ";
					display_fixed_point_switch(Sd, output, number_from_zero);
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output <<  "Tu = ";
					display_fixed_point_switch(Tu, output, number_from_zero);
					output << " Td = ";
					display_fixed_point_switch(Td, output, number_from_zero);
					output << endl;
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output <<  "Inverse invSu = ";
					display_fixed_point_switch(invSu, output, number_from_zero);
					output << " invSd = ";
					display_fixed_point_switch(invSd, output, number_from_zero);
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				}
				matrix_control::SINGLE_LINE_OUTPUT = saved_bool;

//if (debug_control::DEBUG >= debug_control::SUMMARY)
//{
//	debug << (essential_pair? "Essential": "Virtual");
//	debug << " pair\nSu = " << Su << " Sd = " << Sd << "\nTu = " << Tu << " Td = " << Td << endl;
//	debug << "Inverse invSu = " << invSu << " invSd = " << invSd << endl;
//}
				
				/* We could implement evaluating the switch power here but it's probably not very useful */					
//				if (braid_control::SWITCH_POWER)
//				{
//				}

				if (braid_control::RACK_POLYNOMIAL)
				{
					while (get_next_input_string(input_file,input_string, title))
					{

						/* The code here replicates that found in generic-code.cpp, since that function is only called
						   when we are not processing a braid_control::SWITCH_POLYNOMIAL_INVARIANT.  In Version 28.0
						   the RACK_POLYNOMIAL task changed to accept peer codes in addition to braids, so we can 
						   handle explicit turning numbers, using Vogel to calculate an equivalent braid, before 
						   calculating the rack_poly_invariant.  Note that the Vogel algorithm does not change the writhe 
						   or the turning number of the input.
						   
						   If we are given a braid, we assume all strands are closed anti-clockwise, so the turning 
						   nunber equals the number of braid strands.						   						   
						*/

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "braid::main: provided with input string: " << input_string << endl;

					    int turning_number;
					    int writhe;
					    string braid_word;

						if (input_string.find('(') != string::npos || input_string.find('[') != string::npos)
						{
							/* first isolate any qualifiers from the input string */
							string qualifier;						
							string::size_type pos = input_string.find('{');
							if (pos != string::npos)
							{
								qualifier = input_string.substr(pos,string::npos);

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "braid::main: isolated qualifier = " << qualifier << endl;
								input_string = input_string.substr(0,pos);
							}
	
							/* if we have a cycle qualifier, override the global value set by braid_control::INFINITE_CYCLE 
							   with the value determined by the qualifier 
							*/
							braid_control::INFINITE_CYCLE = braid_control::cycle::UNSPECIFIED;
							if (qualifier.find("cycle") != string::npos)
							{
								string::size_type pos = qualifier.find("cycle");
								string::size_type next = qualifier.find(',',pos); //next == string::npos if there's no next qualifier
								string cycle_string = qualifier.substr(pos,next);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "braid::main: cycle_string: " << cycle_string << endl;

							pos = cycle_string.find('=');
					
					
							if (pos != string::npos)
							{
								char* c_cycle_string = c_string(cycle_string);			
								char* cptr = c_cycle_string;	
								while (*cptr != '=')
									cptr++;
									
								get_number(braid_control::INFINITE_CYCLE,++cptr);			
								delete c_cycle_string;
							}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "braid::main: cycle qualifier provided, infinite turning cycle = " << braid_control::INFINITE_CYCLE << endl;

							}

							generic_code_data code_data;
							read_code_data (code_data, input_string);
								
						    /* call the function vogel to evaluate the braid word.  We check first whether all the
						       crossings are specified as flat crossings, in which case we require the Vogel
						       algorithm to use flat Vogel crossings.  If there are any non-flat crossings in the
						       input code, we require the --flat option to be specified if the intent is for the Vogel
						       algorithm to use flat crossings, otherwise it will use classical Vogel crossings.
						    */
						    bool flat_crossings_only = true;
						    int num_crossings = code_data.num_crossings;
						    matrix<int> code_table = code_data.code_table;
						    
						    for (int i=0; i<num_crossings; i++)
						    {
								if (code_table[generic_code_data::table::LABEL][i] != generic_code_data::FLAT)
								{
									flat_crossings_only = false;
									break;
								}
							}
						    
						    if (flat_crossings_only)
						    {
								braid_control::FLAT_CROSSINGS = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::main: input code contains only flat crossings, setting programme option FLAT_CROSSINGS" << endl;
							
							}
						    
//						    int turning_number;
//							string braid_word = vogel(code_data,&turning_number);
							braid_word = vogel(code_data,&turning_number);
							
//							int writhe = code_data_writhe(code_data);
							writhe = code_data_writhe(code_data);
						}
						else
						{
							int num_terms;
							int num_strings;
							if (valid_braid_input(input_string, num_terms, num_strings, braid_control::SILENT_OPERATION, braid_control::RAW_OUTPUT, braid_control::OUTPUT_AS_INPUT))
							{
								braid_word = input_string;
								turning_number = num_strings;
								writhe = count(input_string.begin(),input_string.end(),'s') - 2*count(input_string.begin(),input_string.end(),'-');
							}
							
						}
																	
//						rack_poly_invariant(Su, Sd, invSu, invSd, Tu, Td, pair_type, input_string, title);			
						rack_poly_invariant(Su, Sd, invSu, invSd, Tu, Td, pair_type, braid_word, title, writhe, turning_number);			

					}
				}
				else
				{
					while (get_next_input_string(input_file,input_string, title))
						fixed_point_invariant(Su, Sd, invSu, invSd, Tu, Td, pair_type, input_string, title);			
				}
			}
			else if (braid_control::COMMUTATIVE_AUTOMORPHISM)
			{
				int switch_terms = count(next_switch.begin(),next_switch.end(),',')+1;
				int matrix_size = static_cast<int>(sqrt(float(switch_terms/2)));
				
				/* test matrix_size: if we have missed a comma or some terms, matrix_size will be too small */
				int size_check = 2*matrix_size*matrix_size;
				if (size_check < switch_terms)
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "\nError! Terms or comma missing from switch: " << next_switch << endl;

					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Error! Terms or comma missing from switch: " << next_switch << endl;
					}
					
if (debug_control::DEBUG >= debug_control::SUMMARY)	
    debug << "\nError! Terms or comma missing from switch: " << next_switch << endl;

					continue; // with the next switch
				}

				
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "braid::main: commutative automorphism matrix size = " << matrix_size << endl;

				Qpmatrix phi(matrix_size,matrix_size);
				Qpmatrix psi(matrix_size,matrix_size);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "braid::main: reading phi from next_switch:" << endl;

				/* read in phi */				
				string::size_type posA = next_switch.find('A')+1;
				for (int i=0; i< matrix_size; i++)
				for (int j=0; j< matrix_size; j++)
				{
					string::size_type posB = next_switch.find(',',posA);
					string term = next_switch.substr(posA, (posB==string::npos? posB: posB-posA));

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "braid::main:  next term = " << term << endl;

					phi[i][j] = Qpolynomial(term.substr(0,term.find('/')));

					string::size_type slash_pos = term.find('/');
					if (slash_pos != string::npos)
					{
						phi[i][j] /= Qpolynomial(term.substr(slash_pos+1));
					}
					
					posA = posB+1;
				}
				
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "braid::main: phi = \n" << phi << endl;
	debug << "braid::main: reading psi from next_switch:" << endl;
}

				/* now B */				
				for (int i=0; i< matrix_size; i++)
				for (int j=0; j< matrix_size; j++)
				{
					string::size_type posB = next_switch.find(',',posA);
					string term = next_switch.substr(posA, (posB==string::npos? posB: posB-posA));

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "braid::main:  next term = " << term << endl;

					psi[i][j] = Qpolynomial(term.substr(0,term.find('/')));
					
					string::size_type slash_pos = term.find('/');
					if (slash_pos != string::npos)
					{
						psi[i][j] /= Qpolynomial(term.substr(slash_pos+1));
					}
					
					posA = posB+1;
				}
				
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "braid::main: psi = \n" << psi << endl;



				/* report the switch, we do this explicitly as report_switch is designed for linear switch matrices */

				if (!braid_control::SILENT_OPERATION)
					cout << "\nCommutative automprphism switch ";

				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "\nCommutative automprphism switch ";
				}
				
				if (switch_title != "")
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "\n" << switch_title;
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");					
						output << "\n" << switch_title;
					}
				}

				bool saved_bool = matrix_control::SINGLE_LINE_OUTPUT;
				matrix_control::SINGLE_LINE_OUTPUT = true;
				if (!braid_control::SILENT_OPERATION)
				{
					cout << "\n\\phi = " << endl;
					print(phi, cout, 3, "    ");
					cout << "\n\\psi = " << endl;
					print(psi, cout, 3, "    ");
					cout << endl;
				}
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");					
					output << "\\phi = " << endl;
					print(phi, output, 3, "    ");
					output << "\n\\psi = " << endl;
					print(psi, output, 3, "    ");
				}
				matrix_control::SINGLE_LINE_OUTPUT = saved_bool;


					
				while (get_next_input_string(input_file,input_string, title))
					commutative_automorphism_invariant(phi, psi, input_string, title);			
				
			}
			else
			{
				/* use the Burau switch, or in the case of braid input with the doodle qualifier, the doodle-Alexander representation.
				   In the latter case, we still set the switch matrix to the Burau switch, since we explicitly write the terms of the
				   doodle-Alexander representation in braid_rep, as required.
				   
				   This clause is also the default case where SWITCH_POLYNOMIAL_INVARIANT == false */
				typedef matrix<polynomial<scalar>,scalar> Bmatrix;
				Bmatrix switch_matrix(2,2);
				Bmatrix switch_matrix_inverse(2,2);		
				
				/* Set the scalar variant to bigint unless we are already doing mod-p */
				if (!braid_control::CALCULATE_MOD_P)
					scalar::set_variant(scalar::BIGINT);	
					
//				string str = "0 s t 1-st"; // old Burau definition
				string str = "1-st t s 0"; // current Burau definition
				istringstream ss(str);
				ss >> switch_matrix[0][0] >> switch_matrix[0][1] >> switch_matrix[1][0] >> switch_matrix[1][1];

//				string str_inv = "1-s^-1t^-1 t^-1 s^-1 0"; // old Burau inverse
				string str_inv = "0 s^-1 t^-1 1-s^-1t^-1"; // current Burau inverse
				istringstream ss_inv(str_inv);
				ss_inv >> switch_matrix_inverse[0][0] >> switch_matrix_inverse[0][1] 
				   >> switch_matrix_inverse[1][0] >> switch_matrix_inverse[1][1];
			
				report_switch_matrix(switch_matrix, switch_matrix_inverse, "Burau");
			
				if (braid_control::SWITCH_POWER)
				{
					Bmatrix power = switch_matrix;
					for (int i=1; i< braid_control::SWITCH_POWER; i++)
						power *= switch_matrix;
						
					if (!braid_control::SILENT_OPERATION)
						cout << "\nswitch matrix raised to the power " << braid_control::SWITCH_POWER << ":\n" << power << endl;

					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "switch matrix raised to the power " << braid_control::SWITCH_POWER << ":";
						for (unsigned int i=0; i < power.numrows(); i++)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							for (unsigned int j=0; j < power.numcols(); j++)
								output << power[i][j] << " ";
						}
					}
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "braid::main: switch power = " << braid_control::SWITCH_POWER << "\n" << power << endl;
	
				}
				
				/* The input string here may be a braid or a peer code */
				while (get_next_input_string(input_file,input_string, title))
				{

					string::size_type pos = input_string.find('{');
					
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "braid::input_string = " << input_string << endl;

					if (pos != string::npos && input_string.substr(pos).find("doodle") != string::npos)				
						braid_control::DOODLE_ALEXANDER = true;
					else
						braid_control::DOODLE_ALEXANDER = false;
						
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "braid::braid_control::DOODLE_ALEXANDER = " << braid_control::DOODLE_ALEXANDER << endl;
									
					poly_invariant(switch_matrix, switch_matrix_inverse, input_string, title);	
				}
				
			}
	    }
	}
	else
	{
	    if (input_file_provided)
	    {
			/* reset input ready for read */
			input.clear(); // state flags
			input.seekg(0);
		}

		while (get_next_input_string(input_file,input_string, title))
		{
			if (input_string.find('/') != string::npos || (input_string.find('O') != string::npos && input_string.find('U') != string::npos) || input_string.find('X') != string::npos || input_string.find("DT:") != string::npos)
			{
				generic_code(input_string, title);
			}
			else if (find_if(input_string.begin(), input_string.end(), alpha_char()) != input_string.end())  // if the input string contains an alphabetical character
			{
				braid(input_string, title);
			}
			else  
			{
				cout << "Error! Unknown code " << input_string << endl;
				exit(0);
//				string peer_code = dowker_to_peer_code (input_string);
//				cout << peer_code << endl;
			}
		}
	}
}
catch (const bad_alloc& __e) 
{
	cout << "ERROR! Out of memory, bad_alloc thown: " << __e.what() << endl;
}
catch (...)
{
	cout << "ERROR! Exception thrown that is not caught explicitly" << endl;
}


	/* Evaluate time taken */
	time_t end_time = time(0);
	int time_exp = static_cast<int>(difftime(end_time, start_time));
	int hh = time_exp/3600;
	int mm = (time_exp - hh * 3600)/60;
	int ss = time_exp - hh * 3600 - mm * 60;
	if (!braid_control::SILENT_OPERATION)
		cout << "\nRunning time: " << hh << " hours " << mm	<< " minutes " << ss << " seconds" << endl;

	if (!braid_control::RAW_OUTPUT)
	{
		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
		output << "Running time: " << hh << " hours " << mm	<< " minutes " << ss << " seconds" << endl;
	}

	if (input_file_provided)
		input.close();
	
    output.flush();
    output.close();


if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "\nRunning time: " << hh << " hours " << mm	<< " minutes " << ss << " seconds" << endl;
    debug.close();
}

	return 0;
} /* End of main */

Hpmatrix Hpmatrix_from_q_switch(const Hmatrix& q_switch, bool t_variable)
{
	Hpmatrix switch_matrix(2,2);
	string t = "t";
	string tinv = "t^-1";
	
	switch_matrix[0][0] = Hpolynomial(Quaternion(q_switch[0][0]));
	
	if (t_variable)
	{
		switch_matrix[0][1] = Hpolynomial(Quaternion(q_switch[0][1])) * Hpolynomial(t);
		switch_matrix[1][0] = Hpolynomial(Quaternion(q_switch[1][0])) * Hpolynomial(tinv);
	}
	else
	{
		switch_matrix[0][1] = Hpolynomial(Quaternion(q_switch[0][1]));
		switch_matrix[1][0] = Hpolynomial(Quaternion(q_switch[1][0]));
	}
	
	switch_matrix[1][1] = Hpolynomial(Quaternion(q_switch[1][1]));

	return switch_matrix;
}

void set_programme_long_option(char* cptr, string source)
{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "set_programme_long_option: passed " << cptr << " from " << source << endl;

	char  loc_buf[40];
	char* c1;
	char* c2;

	c1 = cptr;

	/* take out any leading space */
	while (*c1 == ' ')
		c1++;
	
	/* take out any leading -- */
	while (*c1 == '-')
		c1++;
	
	c2 = loc_buf;
	while (isalpha(*c1) || isdigit(*c1) || *c1 == '-' || *c1 == '_')
		*c2++ = *c1++;
    *c2 = '\0';
	
	if (!strcmp(loc_buf,"affine-index"))
	{
   		braid_control::AFFINE_INDEX = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: AFFINE_INDEX read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"alexander"))
	{
   		braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
		braid_control::ALEXANDER = true;
		braid_control::BURAU = false;
		braid_control::COMMUTATIVE_AUTOMORPHISM = false;
    	braid_control::FIXED_POINT_INVARIANT = false;
		braid_control::MATRIX = false;
		braid_control::QUATERNION = false;
   		braid_control::RACK_POLYNOMIAL = false;
		braid_control::WEYL = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: SWITCH_POLYNOMIAL_INVARIANT read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"arrow-polynomial"))
	{
   		braid_control::ARROW_POLYNOMIAL = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: ARROW_POLYNOMIAL read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"automorphism"))
	{
   		braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
		braid_control::ALEXANDER = false;
		braid_control::BURAU = false;
		braid_control::COMMUTATIVE_AUTOMORPHISM = true;
    	braid_control::FIXED_POINT_INVARIANT = false;
		braid_control::MATRIX = false;
		braid_control::QUATERNION = false;
   		braid_control::RACK_POLYNOMIAL = false;
		braid_control::WEYL = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: COMMUTATIVE_AUTOMORPHISM read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"burau"))
	{
   		braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
		braid_control::ALEXANDER = false;
		braid_control::BURAU = true;
		braid_control::COMMUTATIVE_AUTOMORPHISM = false;
    	braid_control::FIXED_POINT_INVARIANT = false;
		braid_control::MATRIX = false;
		braid_control::QUATERNION = false;
   		braid_control::RACK_POLYNOMIAL = false;
		braid_control::WEYL = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: BURAU read from " << source << endl;
	}	
	else if (!strcmp(loc_buf,"classical"))
	{
    	braid_control::CLASSICAL_ONLY = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: CLASSICAL_ONLY read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"complex-delta1"))
	{
    	braid_control::COMPLEX_STUDY_DELTA_1= true;
		braid_control::QUATERNION = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: COMPLEX_STUDY_DELTA_1 read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"delta1-only"))
	{
    	braid_control::DISPLAY_DELTA_1_ONLY = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: DISPLAY_DELTA_1_ONLY read from " << source << endl;
	}
/*    else if (!strcmp(loc_buf,"doodle"))
	{
		braid_control::DOODLE_CONDITIONS = true;
		braid_control::FLAT_CROSSINGS = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_options: DOODLE_CONDITIONS read from " << source << endl;
	} */
    else if (!strcmp(loc_buf,"development"))
	{
		braid_control::DEVELOPMENT_MODE = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: DEVELOPMENT_MODE read from " << source << endl;
	}
    else if (!strcmp(loc_buf,"double-braid"))
	{
		braid_control::KAMADA_DOUBLE_COVERING = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: KAMADA_DOUBLE_COVERING read from " << source << endl;
	}
    else if (!strcmp(loc_buf,"dowker"))
	{
		braid_control::DOWKER_CODE = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: DOWKER_CODE read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"dynnikov"))
	{
    	braid_control::DYNNIKOV_TEST = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: DYNNIKOV_TEST read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"equality"))
	{
    	braid_control::EQUALITY_TEST = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: EQUALITY_TEST read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"extra-output"))
	{
    	braid_control::EXTRA_OUTPUT = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: EXTRA_OUTPUT read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"fixed-point"))
	{
		braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
		braid_control::ALEXANDER = false;
		braid_control::BURAU = false;
		braid_control::COMMUTATIVE_AUTOMORPHISM = false;
    	braid_control::FIXED_POINT_INVARIANT = true;
		braid_control::MATRIX = false;
		braid_control::QUATERNION = false;
   		braid_control::RACK_POLYNOMIAL = false;
    	braid_control::WEYL = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: FIXED_POINT_INVARIANT read from " << source << endl;
	}
    else if (!strcmp(loc_buf,"flat-crossings"))
	{
		braid_control::FLAT_CROSSINGS = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_options: FLAT_CROSSINGS read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"flip-braid"))
	{
		braid_control::FLIP_BRAID = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: FLIP_BRAID read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"format"))
	{
    	braid_control::OUTPUT_AS_INPUT = true;
		braid_control::EXTRA_OUTPUT = false; 
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: braid_control::OUTPUT_AS_INPUT read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"gauss"))
	{
		braid_control::GAUSS_CODE = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: GAUSS_CODE read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"hamiltonian"))
	{
    	braid_control::HAMILTONIAN = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: HAMILTONIAN read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"HC-count"))
	{
    	braid_control::HC_COUNT = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: HC_COUNT read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"HC-edges"))
	{
    	braid_control::HC_EDGES = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: HC_EDGES read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"HC-include-edge"))
	{
		if (*c1 == '=')
		{
			get_number(braid_control::HC_INCLUDE_EDGE,++c1);
		}
		else
		{
			cout << "\nYou must specify an edge to include if you use the HC-include-edge option, e.g. HC-include-edge=12" << endl;
			exit(0);
		}
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: HC_INCLUDE_EDGE read from " << source << ", HC_INCLUDE_EDGE = " << braid_control::HC_INCLUDE_EDGE << endl;
	}
	else if (!strcmp(loc_buf,"HC-list-all"))
	{
    	braid_control::HC_LIST_ALL = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: HC_LIST_ALL read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"homfly"))
	{
    	braid_control::HOMFLY = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: HOMFLY read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"immersion"))
	{
		braid_control::IMMERSION_CODE = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: IMMERSION_CODE read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"info"))
	{
		braid_control::STATUS_INFORMATION = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: STATUS_INFORMATION read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"invert-braid"))
	{
		braid_control::INVERT_BRAID = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: INVERT_BRAID read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"doodle-Q-poly"))
	{
   		braid_control::DOODLE_Q_POLYNOMIAL = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: DOODLE_Q_POLYNOMIAL read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"jones-polynomial"))
	{
   		braid_control::JONES_POLYNOMIAL = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: JONES_POLYNOMIAL read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"kauffman-bracket"))
	{
   		braid_control::KAUFFMAN_BRACKET = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: KAUFFMAN_BRACKET read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"knotoid-bracket"))
	{
   		braid_control::KNOTOID_BRACKET = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: KNOTOID_BRACKET read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"line-reflect-braid"))
	{
		braid_control::LINE_REFLECT_BRAID = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: LINE_REFLECT_BRAID read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"lpgd"))
	{
    	braid_control::LPGD = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: braid_control::LPGD read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"matrix"))
	{
   		braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
		braid_control::ALEXANDER = false;
		braid_control::BURAU = false;
		braid_control::COMMUTATIVE_AUTOMORPHISM = false;
    	braid_control::FIXED_POINT_INVARIANT = false;
		braid_control::MATRIX = true;
		braid_control::QUATERNION = false;
   		braid_control::RACK_POLYNOMIAL = false;
		braid_control::WEYL = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: MATRIX read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"mock"))
	{
   		braid_control::MOCK_ALEXANDER = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: MOCK_ALEXANDER read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"mod-p"))
	{
		braid_control::CALCULATE_MOD_P = true;
		scalar::set_variant(scalar::MOD_P);	
		polynomial_control::MOD_P = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting polynomial_control::MOD_P = true as a result of " << source << " option" << endl;

		int p_value;

		if (*c1 == '=')
		{
			get_number(p_value,++c1);
			mod_p::set_p(p_value);
		}
		else
		{
			cout << "\nYou must specify a prime if you use the mod-p option, e.g. mod-p=7" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: MOD_P read from " << source << ", p = " << p_value << endl;
	}
	else if (!strcmp(loc_buf,"no-auto-delta1"))
	{
		braid_control::ALWAYS_CALCULATE_DELTA_1 = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: clear ALWAYS_CALCULATE_DELTA_1 read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"no-even-writhe"))
	{
    	braid_control::EVEN_WRITHE = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: clear EVEN_WRITHE read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"no-expanded-bracket"))
	{
		braid_control::EXPANDED_BRACKET_POLYNOMIAL = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: clear EXPANDED_BRACKET_POLYNOMIAL read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"normalize-quaternions"))
	{
		braid_control::NORMALIZING_Q_POLYNOMIALS = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: NORMALIZING_Q_POLYNOMIALS read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"no-normalize-bracket"))
	{
		braid_control::NORMALIZE_BRACKET = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: clear NORMALIZE_BRACKET read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"opgc"))
	{
    	braid_control::OPGC = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: braid_control::OPGC read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"OU-format"))
	{
		braid_control::OU_FORMAT = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: OU_FORMAT read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"parity-arrow"))
	{
   		braid_control::PARITY_ARROW = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: PARITY_ARROW read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"parity-bracket"))
	{
   		braid_control::PARITY_BRACKET = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: PARITY_BRACKET read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"PD-format"))
	{
    	braid_control::PD_FORMAT = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: braid_control::PD_FORMAT read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"peer"))
	{
		braid_control::PEER_CODE = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: PEER_CODE read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"plane-reflect-braid"))
	{
		braid_control::PLANE_REFLECT_BRAID = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: PLANE_REFLECT_BRAID read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"power"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting SWITCH_POWER as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(braid_control::SWITCH_POWER,++c1);
			if (braid_control::SWITCH_POWER < 2)
				braid_control::SWITCH_POWER = 2;
		}
		else
		{
			cout << "\nYou must specify an integer n >= 2 if you use the power option, e.g. power=4" << endl;
			exit(0);
		}
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: SWITCH_POWER=" << braid_control::SWITCH_POWER << " read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"prime"))
	{
		braid_control::PRIME_TEST = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: PRIME_TEST read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"quaternion"))
	{
   		braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
		braid_control::ALEXANDER = false;
		braid_control::BURAU = false;
		braid_control::COMMUTATIVE_AUTOMORPHISM = false;
    	braid_control::FIXED_POINT_INVARIANT = false;
		braid_control::MATRIX = false;
		braid_control::QUATERNION = true;
   		braid_control::RACK_POLYNOMIAL = false;
		braid_control::WEYL = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: QUATERNION read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"rack-polynomial"))
	{
		braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
		braid_control::ALEXANDER = false;
		braid_control::BURAU = false;
		braid_control::COMMUTATIVE_AUTOMORPHISM = false;
    	braid_control::FIXED_POINT_INVARIANT = false;
		braid_control::MATRIX = false;
		braid_control::QUATERNION = false;
   		braid_control::RACK_POLYNOMIAL = true;
    	braid_control::WEYL = false;

    	braid_control::VOGEL_TURNING_NUMBER = true;
		braid_control::REDUCE_BRAIDS = false;
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: RACK_POLYNOMIAL read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"raw-output"))
	{
		braid_control::RAW_OUTPUT = true;
		braid_control::EXTRA_OUTPUT = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: braid_control::RAW_OUTPUT read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"relaxed-parity"))
	{
    	braid_control::RELAXED_PARITY = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: braid_control::RELAXED_PARITY read from " << source << endl;
	}	
	else if (!strcmp(loc_buf,"remove"))
	{
		braid_control::REMOVE_PEER_CODE_COMPONENT = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting REMOVE_PEER_CODE_COMPONENT as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(braid_control::REMOVE_COMPONENT,++c1);
		}
		else
		{
			cout << "\nYou must specify an integer n >= 0 if you use the remove option, e.g. remove=2" << endl;
			exit(0);
		}
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: REMOVE_COMPONENT=" << braid_control::REMOVE_COMPONENT << " read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"rho"))
	{
		braid_control::STUDY_RHO_MAPPING = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: STUDY_RHO_MAPPING read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"satellite"))
	{
		if (*c1 == '=')
			get_number(braid_control::SATELLITE,++c1);
		else
			braid_control::SATELLITE = 2;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: SATELLITE  = " << braid_control::SATELLITE << " read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"sawollek"))
	{
    	braid_control::SAWOLLEK = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: SAWOLLEK read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"show-parity-peer-codes"))
	{
    		polynomial_control::WRITE_PARITY_PEER_CODES = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: polynomial_control::WRITE_PARITY_PEER_CODES set from " << source << endl;
	}	
	else if (!strcmp(loc_buf,"show-varmaps"))
	{
    	polynomial_control::SUBSTITUTE_MAPPED_VARIABLES = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: polynomial_control::SUBSTITUTE_MAPPED_VARIABLES cleared from " << source << endl;
	}	
	else if (!strcmp(loc_buf,"silent"))
	{
    	braid_control::SILENT_OPERATION = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: braid_control::SILENT_OPERATION read from " << source << endl;
	}	
	else if (!strcmp(loc_buf,"TeX-polynomials"))
	{
    	braid_control::TeX_POLYNOMIAL_OUTPUT = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: braid_control::TeX_POLYNOMIAL_OUTPUT read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"turning-number"))
	{
    	braid_control::VOGEL_ALGORITHM = true;
    	braid_control::VOGEL_TURNING_NUMBER = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: turning-number read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"ulpgd"))
	{
    	braid_control::ULPGD = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: braid_control::ULPGD read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"uopgc"))
	{
    	braid_control::UOPGC = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: braid_control::UOPGC read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"vogel"))
	{
    	braid_control::VOGEL_ALGORITHM = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: VOGEL_ALGORITHM read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"vogel-height"))
	{
    	braid_control::VOGEL_HEIGHT_ONLY = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: VOGEL_HEIGHT_ONLY read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"wait"))
	{
    	braid_control::WAIT_SWITCH = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: WAIT_SWITCH read from " << source << endl;
   		if (*c1 == '=')
   		{
			int threshold;
			get_number(threshold,++c1);
	    	braid_control::wait_threshold = max(threshold,2);
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "    wait_threshold " << braid_control::wait_threshold << " read from " << source << endl;
   		}
	}
	else if (!strcmp(loc_buf,"weyl"))
	{
   		braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
		braid_control::ALEXANDER = false;
		braid_control::BURAU = false;
		braid_control::COMMUTATIVE_AUTOMORPHISM = false;
    	braid_control::FIXED_POINT_INVARIANT = false;
		braid_control::MATRIX = false;
		braid_control::QUATERNION = false;
   		braid_control::RACK_POLYNOMIAL = false;
    	braid_control::WEYL = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: WEYL_SWITCH read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"zig-zag-delta"))
	{
		braid_control::ZIG_ZAG_DELTA = true;	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: ZIG_ZAG_DELTA read from " << source << endl;
	}
	else
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: invalid long option " << loc_buf << " read from " << source << endl;
	
		cout << "Invalid long option " << loc_buf << endl;
		exit(0);		
	}

    /* Make sure we've not overspecified the options */
	if (braid_control::SWITCH_POLYNOMIAL_INVARIANT + braid_control::DOWKER_CODE + braid_control::GAUSS_CODE + 
	    braid_control::IMMERSION_CODE + braid_control::DYNNIKOV_TEST + braid_control::SAWOLLEK + braid_control::VOGEL_ALGORITHM > 1)
	{
		cout << "\n\nError: incompatible options specified after reading " << source << endl;
		if (braid_control::ALEXANDER)
			cout << "Alexander version of ";
		if (braid_control::SWITCH_POLYNOMIAL_INVARIANT)
			cout << "Polynomial invariants selected\n";
		if (braid_control::FIXED_POINT_INVARIANT)
			cout << "Fixed point invariant selected\n";
		if (braid_control::DOWKER_CODE)
			cout << "Dowker code selected\n";
		if (braid_control::GAUSS_CODE)
			cout << "Gauss code selected\n";
		if (braid_control::IMMERSION_CODE)
			cout << "Immersion code selected\n";
		if (braid_control::DYNNIKOV_TEST)
			cout << "Dynnikov test selected\n";		
		if (braid_control::RACK_POLYNOMIAL)
			cout << "Rack polynomial invariants selected\n";
		if (braid_control::SAWOLLEK)
			cout << "Sawollek polynomial selected\n";
		if (braid_control::VOGEL_ALGORITHM)
			cout << "Vogel algorithm selected\n";
		exit(0);
	}		

}


void set_programme_short_option(char* cptr)
{
	if (strchr(cptr, 'H') && strchr(cptr, '!'))
	{
		if (strchr(cptr,'#'))
		{
			debug_help();
		}
		else
		{
			cout << "\nUsage braid --<task> [-<short_options>][--<long_option>][<infile>[<outfile>]]\n\n";

			cout << "<task> =\n";
			cout << "  affine-index: affine index polynomial invariant\n";
			cout << "  alexander: Alexander polynomial invariant\n";
			cout << "  arrow-polynomial: the arrow polynomial invariant of a classical or virtual knot, link, knotoid or multi-knotoid\n";
			cout << "  automorphism: evaluate a commutative automorphism switch polynomial invariant\n";
	    	cout << "  burau: Burau polynomial invariant\n";
	    	cout << "  doodle-Q-poly: The Q-polynomial for doodles with one component\n";
	    	cout << "  dowker: Dowker code\n";
			cout << "  dynnikov: Dynnikov test\n";
			cout << "  fixed-point: Fixed-point invariant of a braid\n";
	    	cout << "  gauss: Gauss code\n";
			cout << "  hamiltonian: determine a Hamiltonian circuit in the shadow of a knot or link\n";
			cout << "  homfly: HOMFLY polynomial of a braid\n";
			cout << "  immersion: immersion code\n";
			cout << "  info: display status information about the braid\n";
			cout << "  jones-polynomial: calculate the Jones polynomial\n";
			cout << "  kauffman-bracket: calculate the normalized Kauffman bracket polynomial\n";
			cout << "  knotoid-bracket: calculate the Turaev extended bracket polynomial of a knotoid\n";
			cout << "  matrix: matrix-switch polynomial invariant\n";
			cout << "  mock: mock Alexander polynomial invariants of knotoids\n";
			cout << "  parity-arrow: calculate the normalized parity arrow polynomial\n";
			cout << "  parity-bracket: calculate the normalized parity bracket polynomial\n";
			cout << "  peer: peer code\n";
			cout << "  prime: determine whether a given diagram is prime; i.e is not a connected sum or has a 3-connected shadow\n";
	    	cout << "  quaternion: quaternionic-switch polynomial invariant\n";
			cout << "  rack-polynomial: calculate the rack-polynomial invariant\n";
			cout << "  sawollek: calculate Sawollek's normalized Conway polynomial\n";
			cout << "  turning-number: evaluate the turning number of a given diagram\n";
			cout << "  vogel: Vogel algorithm\n";
			cout << "  weyl: Weyl-algebra-switch polynomial invariant\n\n";

			cout << "A <task> is just an example of a <long_option>. Any of the programme option keywords\n";
			cout << "that may appear as an input file programme option may be used as a <long_option>.  The other\n";
			cout << "<long_options> available are:\n\n";

			cout << "  complex-delta1         calculate Delta_1^C rather than Delta_1^H for quaternionic switches\n";
			cout << "  delta1-only            display polynomial output for Delta_1 only\n";
			cout << "  double-braid           calculate the Kamada double covering of all the braids in the input file\n";
			cout << "  equality               test for A=D or B=C in switch when calculating switch polynomial invariants\n";
			cout << "  extra-output           display additional polynomial invariant output\n";
			cout << "  fixed-point            finite biquandle fixed point invariant\n";
			cout << "  flat-crossings         create flat Reidemeister II moves when executing the Vogel algorithm\n";
			cout << "                         consider crossings to be flat when testing for prime knots, so include crossing test\n";
			cout << "  flip-braid             flip all the braids in the input file\n";
			cout << "  format                 format the output file so that it may be used as an input file later\n";
			cout << "  invert-braid           invert all the braids in the input file\n";
			cout << "  HC-count               count the number of Hamiltonian circuits in a diagram\n";
			cout << "  HC-edges               create Hamiltonian circuits from edges rather than crossings\n";
			cout << "  HC-list-all            find all of the Hamiltonian circuits in a diagram\n";
			cout << "  line-reflect           reflect all the braids in the file in a horizontal line drawn south of the braid\n";
			cout << "  lpgd                   calculate the left preferred Gauss code, rather than a standard gauss code\n";
			cout << "  mod-p=n                calculate mod p with p=n (only used for non-Weyl algebra switches)\n";
			cout << "  no-auto-delta1         only calculate Delta_1 if Delta_0 is zero\n";
			cout << "  no-even-writhe         normalize the parity bracket polynomial with the full writhe rather than the even writhe\n";
			cout << "  no-expanded-bracket    do not expand D=(-A^2-A^{-2}) in bracket polynomials\n";
			cout << "  no-normalize-bracket   do not normalize bracket polynomial invariants\n";
			cout << "  normalize-quaternions  normalize quaternionic polynomial invariants\n";
			cout << "  opgc                   calculate the over preferred Gauss code, rather than a standard Gauss code\n";
			cout << "  OU-format              write Gauss codes as a sequence (O|U)<crossing-num><crossing-sign>\n";
			cout << "  PD-format              write Gauss code as a planar diagram\n";
			cout << "  plane-reflect          reflect all the braids in the file in the plane of the page\n";
			cout << "  power=n                evaluate the nth power of the switch when calculating switch polynomial invariants\n";
			cout << "  raw-output             produce raw output, that is the result only without descriptive text\n";
			cout << "  relaxed-parity         evaluate the relaxed variant of the parity arrow polynomial\n";
			cout << "  remove=n               remove the n-th component from a peer code\n";
			cout << "  rho                    use the Study rho mapping for calculating Study determinants\n";
			cout << "  satellite[=n]:         determine the peer code of the n-parallel cable satellite of a knot's peer code before carrying out the required programme task\n";
			cout << "                         n: default value 2\n";
			cout << "  show-parity-peer-codes show peer codes in addition to unoriented left preferred Gauss codes in parity bracket polynomial output\n";
			cout << "  show-varmaps           show variable mappings instead of substituting mapped variables in polynomial output\n";
			cout << "  silent                 do not generate any output to the command line (stdout)\n";
			cout << "  TeX-polynomials        display output polynomials in TeX format\n";
			cout << "  ulpgd                  calculate the unoriented left preferred Gauss code, rather than a standard gauss code\n";
			cout << "  uopgc                  calculate the unoriented over preferred Gauss code, rather than a standard gauss code\n";
			cout << "  wait[=n]               display determinant wait information, (based on nxn minors, so larger n produces less frequent output)\n";
			cout << "  zig-zag-delta          include delta with K_i and Lambda_i variables when calculating the arrow polynomial\n\n";
		
			cout << "<short_option> =\n";
	    	cout << "  #: debug\n";
			cout << "  0: calculate Delta_0 only\n";
			cout << "  c[{2}]: complex Study Delta_1 (turns on q)\n";
			cout << "     {2}: always calculate codimension 2 determinant from complex Study Delta_1\n";
			cout << "  d: evaluate the Kamada double covering for braids (used only by the fixed-point task)\n";
			cout << "  D: verify Delta_0 = 0 in the classical case\n";
	    	cout << "  e: Test for the switch equality condition A=D and B=C\n";
			cout << "  F: bypass the fundamental equation check for switches\n";
	    	cout << "  h: help screen\n";
	    	cout << "  H!: this help screen\n";
	    	cout << "  #H!: display debug help screen\n";
			cout << "  I: format output as a valid input file\n";
			cout << "  M: do not remove Reidemeister II moves when calculating labelled peer codes from braid words\n";
			cout << "  N: normalize quaternionic polynomial invariants\n";
			cout << "  o: additional output\n";
			cout << "  O: raw output\n";
	    	cout << "  p=n: use coefficients mod p, where p=n\n";
	    	cout << "  P: display polynomial wait indicator\n";
//					cout << "r: real Study matrix Delta_1 (turns on q)\n";
			cout << "  R: use rho-mapping for Study determinants\n";
			cout << "  S: silent operation\n";
//			cout << "  t=n: set the number of additional +ve/-ve terms for the rack-polynomial invariant (default " << braid_control::RACK_TERMS << ")\n";
			cout << "  T: output in /tmp\n";
			cout << "  U: do not test for units in calculation of Delta_1\n";
			cout << "  V: Do not use the t-variable with quaternionic switches\n";
			cout << "  W[=n]: force wait information to be displayed, if n is supplied, set the wait threshold to n\n";
			cout << "  x: test mode, for develpment testing only\n";
			cout << "  z: Do not calculate Delta_1 when Delta_0 is non-zero\n";
	    	cout << "  Z: display Delta_1 polynomials only\n";
	    	exit(0);
		}
	}

	/* First we check where the output and debug files should be created */
	if (strchr(cptr, 'T'))
		braid_control::TMP_DIRECTORY = true;

	char* dptr = strchr(cptr, '#');
	if (dptr)
	{

    	/* establish a debug file */
		if (braid_control::TMP_DIRECTORY)
	    	debug.open ("/tmp/braid.dbg"); 
		else	
	    	debug.open ("braid.dbg"); 

    	if (!debug)
    	{
        	cout << "\nError opening debug file\n";
        	exit(0);
    	}
		else
			debug << "Debug information from braid version " << braid_control::version << "\n\n";

		if (!debug_setup(cptr))  // could probably be dptr, but the original code used cptr
		{
			debug_control::DEBUG = debug_control::SUMMARY;
			debug << "set_programme_short_option: default debug options set" << endl;
		}

	}

	if (strchr(cptr, '0'))
		braid_control::CALCULATE_DELTA_0_ONLY = true;


	char* c_option = strchr(cptr, 'c');
	if (c_option)
	{
		braid_control::COMPLEX_STUDY_DELTA_1 = true;
		braid_control::SWITCH_POLYNOMIAL_INVARIANT = true;
	    braid_control::QUATERNION = true;
		set_cli_suboptions(c_option, "c");
	}

	if (strchr(cptr, 'B'))
		braid_control::BIGELOW_KNOT_SEARCH = true;

	if (strchr(cptr, 'd'))
		braid_control::KAMADA_DOUBLE_COVERING = true;

	if (strchr(cptr, 'D'))
		braid_control::VERIFY_DELTA_0 = true;

	if (strchr(cptr, 'e'))
		braid_control::EQUALITY_TEST = true;

	if (strchr(cptr, 'F'))
		braid_control::BYPASS_FUNDAMENTAL_EQUATION_CHECK = true;

	if (strchr(cptr, 'h'))
		help_info();

	if (strchr(cptr, 'I'))
	{
		braid_control::OUTPUT_AS_INPUT = true;
		braid_control::EXTRA_OUTPUT = false; // turns off additional poly-invariant output
	}

	if (strchr(cptr, 'M'))
		braid_control::REMOVE_REIDEMEISTER_II_MOVES = false;

	if (strchr(cptr, 'N'))
		braid_control::NORMALIZING_Q_POLYNOMIALS = true;

	if (strchr(cptr, 'o'))
		braid_control::EXTRA_OUTPUT = true;

	if (strchr(cptr, 'O'))
	{
		braid_control::RAW_OUTPUT = true;
		braid_control::EXTRA_OUTPUT = false;
	}

	if (strchr(cptr, 'p'))
	{
		braid_control::CALCULATE_MOD_P = true;
		scalar::set_variant(scalar::MOD_P);	
		polynomial_control::MOD_P = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
debug << "set_programme_short_option: setting polynomial_control::MOD_P = true as a result of command line option" << endl;

	    char* c1 = strchr(cptr,'p')+1;
	    if (*c1 == '=')
	    {
			int modulus;
			get_number(modulus,++c1);
			mod_p::set_p(modulus);
	    }
		else
		{
			cout << "\nYou must specify a prime if you use the mod p option, e.g. 'braid -qp=7'" << endl;
			exit(0);
		}
	}

	if (strchr(cptr, 'P'))
		polynomial_control::WAIT_INFO = true;

//			if (strchr(cptr, 'r'))
//			{
//		    	REAL_STUDY_DELTA_1 = true;
//		    	SWITCH_POLYNOMIAL_INVARIANT = true;
//	    		QUATERNION = true;
//			}

	if (strchr(cptr, 'R'))
		braid_control::STUDY_RHO_MAPPING = true;

	if (strchr(cptr, 'S'))
		braid_control::SILENT_OPERATION = true;

/*	if (strchr(cptr, 't'))
	{
	    char* c1 = strchr(cptr,'t')+1;
	    if (*c1 == '=')
		{
			get_number(braid_control::RACK_TERMS,++c1);
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: setting RACK_TERMS = " << braid_control::RACK_TERMS << " as a result of command line option" << endl;
		}
		
		if (!braid_control::RACK_TERMS)
		{
			cout << "\nYou must specify the number of terms if you use the t option, e.g. 'braid -t=7'" << endl;
			exit(0);
		}
	}
*/		
	if (strchr(cptr, 'U'))
		braid_control::DELTA_1_UNIT_CHECK = false;

	if (strchr(cptr, 'V'))
		braid_control::T_VARIABLE = false;

	if (strchr(cptr, 'W'))
	{
		braid_control::WAIT_SWITCH = true;
	    char* c1 = strchr(cptr,'W')+1;
	    if (*c1 == '=')
	    {
			int threshold;
			get_number(threshold,++c1);
		    braid_control::wait_threshold = max(threshold,2);
	    }
	}

	if (strchr(cptr, 'x'))
	    braid_control::TEST_MODE = true;

	if (strchr(cptr, 'z'))
	    braid_control::ALWAYS_CALCULATE_DELTA_1 = false;

	if (strchr(cptr, 'Z'))
		braid_control::DISPLAY_DELTA_1_ONLY = true;
}


bool get_next_input_string(string input_file,string& input_string, string& title)
{
	bool success = true;

	if (input_file.length())
	{
		/* get next word from input */
		get_input_word(input, input_string, title);
	}
	else
	{
		do
		{
			cout << "\n\ninput: ";
			getline(cin, input_string); // has to be getline so there can be spaces in the string
		} while (input_string.length() == 0);
	}
	
	if (input_string == "exit" || input_string == "q" || input_string =="Q")
		success = false;
	else if (input_string == "help")
	{
	    help_info();
		exit(0);
	} 	
	else    	
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "braid::get_next_input_string: got string " << input_string << endl;
}
	}	
	
//	if (!braid_control::SILENT_OPERATION && success && title.length() == 0)
//		cout << "\n\n" << input_string << endl;

	return success;
}

void help_info()
{
	ifstream helpfile;

	/* check the helpfile is present by attempting to open it */
	helpfile.open ("helpinfo.txt");

	if (!helpfile)
	{
		cout << "\nI can't find helpinfo.txt, it is available online at";
		cout << "\nwww.layer8.co.uk/maths" << endl;
		exit(0);
	}
	else
	{
		helpfile.close();
		system("more helpinfo.txt");
	}
	exit(0);
}

Hmatrix H22inverse(Hmatrix& M)
{
	Quaternion& A = M[0][0];
	Quaternion& B = M[0][1];
	Quaternion& C = M[1][0];
	Quaternion& D = M[1][1];
	
	Hmatrix invM(2,2);
	
    /* check existence of inverse: i.e. C^-1D != A^-1B */
    if(conj(C)*D*Quaternion(norm_sqrd(A)) == conj(A)*B*Quaternion(norm_sqrd(C)))
    {
		throw(matrix_error("attempt to calculate inverse of a non-invertible 2x2 quaternionic matrix"));
	}
	else
	{	
		Quaternion C_inv = conj(C) * (Quaternion(1)/Quaternion(norm_sqrd(C)));
		Quaternion A_inv = conj(A) * (Quaternion(1)/Quaternion(norm_sqrd(A)));
		Quaternion delta = C_inv * D - A_inv * B;
		Quaternion delta_inv = conj(delta) * (Quaternion(1)/Quaternion(norm_sqrd(delta)));

		invM[0][0] = C_inv * D * delta_inv * A_inv;
		invM[0][1] = Quaternion(-1) * A_inv * B * delta_inv * C_inv;
		invM[1][0] = Quaternion(-1) * delta_inv * A_inv;
		invM[1][1] = delta_inv * C_inv;
		
	    return invM;
    }
}

template <typename T, typename St> void report_switch_matrix (const matrix<T,St>& switch_matrix, const matrix<T,St>&
switch_matrix_inverse, string message)
{

	if (!braid_control::SILENT_OPERATION)
	{
		if (braid_control::WEYL)
			cout << "\n" << message << " switch S, where\n";
		else
			cout << "\n" << message << " switch S, with inverse S', where\n";
			
		cout << "\nS = " << switch_matrix;
		
		if (!braid_control::WEYL)
			cout << "\nS' = " << switch_matrix_inverse;
	
		/* We use polynomial_control::MOD_P to determine mod-p calculation since CALCULATE_MOD_P is
		   only true if mod-p has been set globally, which is not the case when dealing with Weyl algebra
		   switches, which may specify mod-p on a per-switch basis. */
		if (polynomial_control::MOD_P)
			cout << "\n\nCalculating mod " << mod_p::get_p();
				
		cout << endl;
	}
					
	if (!braid_control::RAW_OUTPUT)
	{
		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
		
		if (braid_control::WEYL)
			output << "Using " << message << " switch S, where";
		else
			output << "Using " << message << " switch S, with inverse S', where";
			
		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
		output << "S = ";
		for (size_t i=0; i< switch_matrix.numrows(); i++)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			for (size_t j=0; j< switch_matrix.numcols(); j++)
			{
				output << switch_matrix[i][j] << " ";
			}
		}

		if (!braid_control::WEYL)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n\n;" : "\n\n");
			output << "S'= ";
			for (size_t i=0; i< switch_matrix_inverse.numrows(); i++)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				for (size_t j=0; j< switch_matrix_inverse.numcols(); j++)
				{
					output << switch_matrix_inverse[i][j] << " ";
				}
			}
		}
		
		if (polynomial_control::MOD_P)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n\n;" : "\n\n");
			output << "Calculating mod " << mod_p::get_p();
		}

		output << endl;
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	if (braid_control::WEYL)
		debug << "\nUsing " << message << " switch S, where\n";
	else
		debug << "\nUsing " << message << " switch S, with inverse S', where\n";
	
	debug << "\nS = " << switch_matrix;

	if (!braid_control::WEYL)
		debug << "\nS'= " << switch_matrix_inverse;
	
	if (polynomial_control::MOD_P)
		debug << "\n\nCalculating mod " << mod_p::get_p();

	debug  << endl;
}


}


/* poly_invariant calculates the Delta_0 and Delta_1 determined by the supplied switch matrix and its inverse for the 
   braid or immersion code given in the input string.  It is used only where the elements of the switch matrix are 
   of type polynomial<T,V,U> and not rational<polynomial<T,V,U> > ,i.e. Burau or quaternionic switches currently.  Note that 
   this means that the switch matrix N-factor is one in this function, which consequently is optimized for this case.
*/
template <typename T, typename St> 
void poly_invariant(const matrix<T,St>& switch_matrix, const matrix<T,St>& switch_matrix_inverse, 
                    string input_string, string title)
{
	int num_terms;
	int num_strings;
	
	bool calculate_polys = true;
	bool virtual_crossings_present = false;
	
	matrix<T,St>*		Matrix_rep=0;
	
	if (title.length())
   	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\n\n" << title << "\n" << input_string << endl;
		if (!braid_control::RAW_OUTPUT)
			output << "\n\n" << title << "\n" << input_string << endl;
   	}		

	if (!braid_control::RAW_OUTPUT)
	{
		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
		output << input_string;		
	}
	
	/* note the N-factor of the switch matrix, this is done explicitly as a debug check */
	int switch_matrix_N_factor = switch_matrix.numcols()/2;
	
if(debug_control::DEBUG >= debug_control::DETAIL)
	debug << "poly_invariant: switch_matrix_N_factor = " << switch_matrix_N_factor << endl;

	/* If we're calculating the doodle-Alexander polynomial and have been given a peer code,
	   convert the input string to a braid
	*/
	if (braid_control::DOODLE_ALEXANDER && (input_string.find('(') != string::npos || input_string.find('[') != string::npos))
	{
		/* remove any unwanted qualifiers from the input string */
		string::size_type pos = input_string.find('{');
		if (pos != string::npos)
			input_string = input_string.substr(0,pos);
			
if(debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "poly_invariant: converting peer or immersion code to a braid for DOODLE_ALEXANDER" << endl;
			
		generic_code_data code_data;
		read_code_data(code_data,input_string);
		input_string = vogel(code_data);

if(debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "poly_invariant: input_string set to " << input_string << endl;		
	}
	
	if (input_string.find('(') != string::npos || input_string.find('[') != string::npos)
	{
		/* record whether this is a virtual knot or a classical one */
		if (input_string.find('*') != string::npos)
			virtual_crossings_present = true;

		/* check for a long knot */
		
		if (input_string.find("L:") != string::npos)
		{
			braid_control::LONG_KNOT = true;
			
			input_string = parse_long_knot_input_string(input_string);
			
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "poly_invariant: long knot: " << input_string << endl;
	
			if (input_string == "link")
			{				
				if (!braid_control::SILENT_OPERATION)
					cout << "\nError! long knot indicator provided for the peer code of a link" << endl;
			
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Error! long knot indicator provided for the peer code of a link" << endl;
				}
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "poly_invariant: Error! long knot indicator provided for the peer code of a link" << endl;

				return;
			}
		}
		else
			braid_control::LONG_KNOT = false;
	
		R_module_rep(Matrix_rep, switch_matrix, switch_matrix_inverse, input_string);
	}
	else if (valid_braid_input(input_string, num_terms, num_strings, braid_control::SILENT_OPERATION, braid_control::RAW_OUTPUT, braid_control::OUTPUT_AS_INPUT))
	{ 

		/* remove any unwanted qualifiers from the input string */
		string::size_type pos = input_string.find('{');
		if (pos != string::npos)
			input_string = input_string.substr(0,pos);


		/* record whether this is a virtual knot or a classical one */
		if (input_string.find('t') != string::npos || input_string.find('T') != string::npos)
			virtual_crossings_present = true;
			

		/* When using a braid representation, we do not calculate polynomials if the representation is the identity
		   matrix.  This check is done within braid_rep */	   
		calculate_polys = braid_rep(Matrix_rep, switch_matrix, switch_matrix_inverse,input_string, num_terms, num_strings);
	} 

	/* Under normal circumstances Matrix_rep is created by R_module_rep or braid_rep.  However, if R_module_rep
	   determines that an immersion code contains only virtual crossings it returns a zero Matrix_rep pointer. */
	if (Matrix_rep == 0)
		return;
		
	/* Note the size of the matrix representation, this is either square or in the case of long knots
	   has one more column than row, so we note the number of rows.  For elements of VB(n) this is the 
	   number of strings, but for R_module representations it is twice the number of classical crossings.
	*/
	int matrix_rows = Matrix_rep->numrows();
	int matrix_cols = Matrix_rep->numcols();

    if (calculate_polys)
    {	
		/* First determine the codimension 0 invariant.  In the case of a long knot we have to remove successive columns  
		   and calculate the gcd of the (Study) determinants of remianing submatrices.  In the case of a closed knot
		   the codimension 0 invariant is just the (Study) determinant of the matrix representation.
		   
		   We then determine whether to calculate delta_1 based on the codimension zero value and the global booleans
		*/
		bool calc_delta_1;

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

		/* We use delta_0 for the codimension 0 polynomial invariant in both the closed
		   knot and long knot cases.  This allows us to use common code for delta_1, the
		   codimension 1 polynomial invariant.
		*/ 		  
		T delta_0 = T("0");

		/* If we're evaluating the DOODLE_ALEXANDER polynomial, since num_strings >= 2, we simply divide the determinant 
		   of the matrix representation by -x^{num-strings-1}U_{num_strings-1}(1/x), where U_n(z) is the Chebyshev polynomial:
		   U_0(z)=1, U_1(z)=2z, U_n(z) = 2zU_{n-1} -U_{n-2}
		   
		*/
		if (braid_control::DOODLE_ALEXANDER)
		{
			delta_0 = determinant(*Matrix_rep, title);
			calc_delta_1 = false;

			T Un_minus_2 = T("1");
			T Un_minus_1 = T("2x^-1");
			T Un = Un_minus_1;  // initializes to U_1
			
			/* if num_strings == 2, we want U_1, which is what Un_minus_1 was intialized with.  If num_strings > 2,
			   we execute this loop num_strings - 2 times, giving U_{num_strings-1}
			*/
			for (int i=0; i< num_strings-2; i++)
			{
				Un = T("2x^-1")*Un_minus_1 - Un_minus_2;
				Un_minus_2 = Un_minus_1;
				Un_minus_1 = Un;
			}
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "poly_invariant: U_(n-1) = " << Un << endl;
			
			for (int i=0; i< num_strings-1; i++)
				Un *= T("-x");

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "poly_invariant: (-x)^(n-1) scaled Un = " << Un << endl;
			
			delta_0 /= Un;
			
			if (!braid_control::SILENT_OPERATION)
				cout << "\nDoodle-Alexander polynomial = ";
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Doodle-Alexander polynomial = ";
			}

			if (braid_control::OUTPUT_AS_INPUT)
				output << "\n";
			if (!braid_control::SILENT_OPERATION)
				cout << delta_0 << "\n";
			output << delta_0 << "\n";
			
		}
		else if (braid_control::LONG_KNOT)
		{	
			/* The codimension 0 polynomial invariant is formed by taking out each column of the 
			   matrix representation in turn and taking the determinant  of the rest, keeping a 
			   record of the gcd of the determinants calculated.
			*/
			int rperm[matrix_rows];
			int cperm[matrix_cols];
			bool unit_detected = false;

			for (int i = 0 ; i< matrix_rows; i++)
		    	rperm[i] = i;

			T op_0 = T("0"); // \def\op{\buildrel o \over p}
			T np_0 = T("0"); // \def\np{\buildrel n \over p}

			for (int i = 0 ; i< matrix_cols; i++)
			{
//				if (WAIT_SWITCH)
//				    matrix_control::wait_count = 0;

				/* take out column i */
				for (int j=0; j<i; j++)
    				cperm[j] = j;
				for (int j = i+1; j<= matrix_cols; j++)
					cperm[j-1] = j;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "poly_invariant: removing column " << i << endl;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "poly_invariant: cperm: ";
   	for (int j = 0 ; j< matrix_cols-1; j++)
		debug << cperm[j] << " ";
	debug << endl;
}
				T det = T("0");

				if (braid_control::QUATERNION)
					det =  study_determinant(*Matrix_rep, title, matrix_rows, rperm, cperm);
				else
					det =  determinant(*Matrix_rep, title, matrix_rows, rperm, cperm);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "poly_invariant: p_0 generator det " << i << " = " << det << endl;

    			Laurent_scale (det);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "poly_invariant: Laurent scaled p_0 generator det " << i << " = " << det << endl;


				if (i == 0)
				{
					op_0 = det;

					if (braid_control::WAIT_SWITCH)
					    matrix_control::wait_count = 0;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "poly_invariant: op_0  set to det " << i << " = " << op_0 << endl;

				}
				else if (i == matrix_cols-1)
				{
					np_0 = det;

					if (braid_control::WAIT_SWITCH)
					    matrix_control::wait_count = 0;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "poly_invariant: np_0  set to det " << i << " = " << np_0 << endl;

					if (unit_detected)
						break;  //we've jumped to here, so we're now done.

				}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "poly_invariant:  p_0 gcd stands at " << delta_0 << endl;

				if (braid_control::DELTA_1_UNIT_CHECK && det.isunit())
				{
					unit_detected = true; 
					if (i != matrix_cols-1)
						i = matrix_cols-2; // jump to last column to determine np_0, i will be incremented by loop

					delta_0 = T("1");

					if (braid_control::EXTRA_OUTPUT)
					{
						if (!braid_control::SILENT_OPERATION)
							cout << "\nunit generator detected, terminating calculation" << endl;
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "unit generator detected, terminating calculation" << endl;
						}
					}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "poly_invariant: unit generator detected, terminating calculation" << endl;
				}
				else
				{
			    	delta_0 = gcd(delta_0,det);					

					if (braid_control::DELTA_1_UNIT_CHECK && delta_0.isunit())
					{
						unit_detected = true; 
						if (i != matrix_cols-1)
							i = matrix_cols-2; // jump to last column to determine np_0, i will be incremented by loop

						if (braid_control::EXTRA_OUTPUT)
						{
							if (!braid_control::SILENT_OPERATION)
								cout << "\nunit gcd detected, terminating calculation" << endl;
							if (!braid_control::RAW_OUTPUT)
							{
								output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
								output << "unit gcd detected, terminating calculation" << endl;
							}
						}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "poly_invariant: unit gcd detected, terminating calculation" << endl;
					}
				}
			}

			if (!braid_control::SILENT_OPERATION)
			{
	   			cout << "\np^{(0)} = " << delta_0 << endl;
	   			cout << "\n\\op^{(0)} = " << op_0 << endl;
	   			cout << "\n\\np^{(0)} = " << np_0 << endl;
			}
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "p^{(0)} = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << delta_0 << endl;
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "\\op^{(0)} = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << op_0 << endl;
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "\\np^{(0)} = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << np_0 << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
   	debug << "poly_invariant: p_0  = p^{(0)} = " << delta_0 << endl;
   	debug << "poly_invariant: \\op^{(0)} = " << op_0 << endl;
   	debug << "poly_invariant: \\np^{(0)} = " << np_0 << endl;
}

			if (delta_0 == T("0"))		
				calc_delta_1 = true;
    		else
				calc_delta_1 = false;

		}
		else // closed knot
		{
			/* Delta_0 = 0 in the classical case, so only calculate it for virtuals, the codimension 0 
			   polynomial invariant is just the determinant of the presentation matrix, since it is a square matrix.
			*/
			if (virtual_crossings_present || braid_control::VERIFY_DELTA_0)
			{
		    	if (!braid_control::QUATERNION)
	    		{
					/* use deteminant function to calculate det(M-I) */
					delta_0 = determinant(*Matrix_rep, title);

					if (braid_control::ALEXANDER)
    				{
						/* set s=1 in delta if it involves s */
						if (delta_0.numvars() && strchr(delta_0.getvars(),'s'))
							set_to_one(delta_0,'s');
    				}

					Laurent_scale (delta_0);

					if (!braid_control::DISPLAY_DELTA_1_ONLY)
					{
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
					}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
   	debug << "\npoly_invariant: Delta_0 = " << delta_0 << "\n";
}
	    		}
	    		else
	    		{
					/* we have to use the Study determinant */
					matrix<T,St>& matrix_rep = *Matrix_rep;
					delta_0 = study_determinant(matrix_rep, title);

					if (braid_control::NORMALIZING_Q_POLYNOMIALS)
						normalize(delta_0, switch_matrix);

					if (!braid_control::DISPLAY_DELTA_1_ONLY)
					{
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

					}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
   	debug << "\npoly_invariant: Study determinant = " << delta_0 << "\n";
}

	    		}

				if (!braid_control::SILENT_OPERATION && matrix_control::WAIT_INFO)
					cout << "\n";

				if (delta_0 == T("0"))		
					calc_delta_1 = true;
	    		else
					calc_delta_1 = false;

			}
			else
			{
	    		calc_delta_1 = true;

				if (!braid_control::DISPLAY_DELTA_1_ONLY)
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "\nInput contains no virtual crossings, so Delta_0 = 0" << endl;
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Input contains no virtual crossings, so Delta_0 = ";
						if (braid_control::OUTPUT_AS_INPUT)
							output << "\n";
					}
					output << "0" << endl;

				}
			}
		}
		
		if (!braid_control::DOODLE_ALEXANDER && !braid_control::CALCULATE_DELTA_0_ONLY && (braid_control::ALWAYS_CALCULATE_DELTA_1 || calc_delta_1))
		{
	    	if (braid_control::EXTRA_OUTPUT)
		    {
				/* if we're doing quaternions we don't need wait state
				   information if we're going to see the adjunct elements
				*/
				if ((braid_control::QUATERNION && matrix_cols > 4) ||
		    		(!braid_control::QUATERNION && matrix_cols > 5))
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "\n";
		    		matrix_control::WAIT_INFO = false;
				}
	    	}

//			if (braid_control::COMPLEX_STUDY_DELTA_1 || braid_control::REAL_STUDY_DELTA_1)
			if (braid_control::COMPLEX_STUDY_DELTA_1)
			{
				
				/* Here we map the (quaternionic) matrix representation into M_{2n,2m}(C[t,t^{-1}]) and  
				   determine the gcd of the codimension 1 determinants of the image.  If this gcd is 0 we 
				   determine the gcd of the codimension 2 determinants of the image.  This is an alternative 
				   approach to the normal procedure of taking the codimension 1 submatrices of M_n,m(H[t,t^{-1}), 
				   mapping each one into M_{2n,2m}(C[t,t^{-1}]), taking determinants and evaluating their gcd.
				
				   First, set up the Complex Study matrix and evaluate delta_1 from it. This call to delta_1 will  
				   remove rows and columns of the matrix in M_{2n,2m}(C).  
				*/
				int rperm[matrix_rows];
				int cperm[matrix_cols];

				for (int i = 0 ; i< matrix_rows; i++)
		    		rperm[i] = i;

				for (int i = 0 ; i< matrix_cols; i++)
		    		cperm[i] = i;

				matrix<T,St> C_study_matrix_rep = C_study_matrix(*Matrix_rep, matrix_rows, matrix_cols, rperm, cperm);
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "\npoly_invariant: complex Study matrix:\n";
	debug << C_study_matrix_rep;
	debug << "\n";
}	

				bool display_delta;
				if (braid_control::COMPLEX_STUDY_DELTA_1)
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "Evaluating the complex Study matrix" << endl;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "\n\npoly_invariant: evaluating the complex Study matrix\n";	

					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Evaluating the complex Study matrix";
					}

					if (braid_control::EXTRA_OUTPUT)
					{
						if (!braid_control::SILENT_OPERATION)
							cout << "\nGenerators of Delta_1 calculated from the complex Study matrix" << endl;
						output << "\nGenerators of Delta_1 calculated from the complex Study matrix" << endl;
					}

					display_delta = false;
					/* It is mandatory to include delta_0 as a generator of delta_1 in this case, since
					   we are dealing with quaternions. */
					T delta = delta_1(C_study_matrix_rep, CODIMENSION_1, title, delta_0);

					int codimension;
					if (delta != T("0") && !braid_control::ALWAYS_CALCULATE_CODIMENSION_2_DELTA_1)
					{
						display_delta = true;
						codimension = CODIMENSION_1;
					}
					else
					{
						if (!braid_control::SILENT_OPERATION && matrix_control::WAIT_INFO)
						{
							cout << "\n";
							cout << "\nDelta_1^C=1" << endl;
						}
						
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "\n\npoly_invariant: Delta_1^C=1";
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "Delta_1^C=1";
						}

						/* take out two rows and columns from C_study_matrixptr at a time */
						delta = delta_1(C_study_matrix_rep, CODIMENSION_2, title, delta_0);

						if (delta == T("1"))
						{
							if (!braid_control::SILENT_OPERATION)
								cout << "\nDelta_2^C=1";

//							cout << "\nTrying real Study matrix" << endl;
if (debug_control::DEBUG >= debug_control::SUMMARY)
{	
//	debug << "\n\npoly_invariant: Delta_2^C=1, trying the real Study matrix" << endl;
	debug << "\n\npoly_invariant: Delta_2^C=1" << endl;
}	

							if (!braid_control::RAW_OUTPUT)
							{
								output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
								output << "Delta_2^C=1" << endl;
//								output << "\nTrying real Study matrix" << endl;
							}

//							REAL_STUDY_DELTA_1 = true;
							
						}
						else
						{
							display_delta = true;
							codimension = CODIMENSION_2;
						}
					}

					if (display_delta)
					{
						if (!braid_control::SILENT_OPERATION)
						{
							cout << "\nDelta_";
							cout << (codimension == 1? "1" : "2");
							cout << "^C = ";
						}
						
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "\nDelta_";
							output << (codimension == 1? "1" : "2");
							output << "^C = ";
						}

						if (braid_control::OUTPUT_AS_INPUT)
							output << "\n";
						if (!braid_control::SILENT_OPERATION)
							cout << delta << "\n";
		   				output << delta << "\n";

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "\npoly_invariant: Delta_";
	debug << (codimension == 1? "1" : "2");
	debug << "^C = " << delta << "\n";	
}
					}
				}


/******************************************************************
				if (REAL_STUDY_DELTA_1)
				{
					display_delta = false;
					cout << "\nEvaluating the real Study matrix";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "\n\npoly_invariant: evaluating the real Study matrix\n";	

					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Evaluating the real Study matrix";
					}

					  calculate the Real Study matrix in M_{4n}(R) and
					  use R_study_matrixptr to point to it.  We do this via
					  a function call just in order to keep the code here simple
					
					matrix<T,St> R_study_matrix_rep = R_study_matrix(C_study_matrix_rep);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "\npoly_invariant: real Study matrix:\n";
	debug << R_study_matrix_rep;
	debug << endl;
}					

					if (braid_control::EXTRA_OUTPUT)
					{
						cout << "\n\nGenerators of Delta_1 calculated from the real Study matrix";
						output << "\n\nGenerators of Delta_1 calculated from the real Study matrix";
					}

					// set wait info true, the Real Study matrix is big! 
					matrix_control::WAIT_INFO = true;
	    			matrix_control::wait_count = 0;

				    * if the wait switch was given, the wait_threshold
	    			   comes from the command line
	    			if (!WAIT_SWITCH)
						wait_threshold = (QUATERNION? 2* DEFAULT_WAIT_THRESHOLD: DEFAULT_WAIT_THRESHOLD);
					*

				    cout << "\nCalculating polynomials, please wait\n";

					T delta = delta_1(R_study_matrix_rep, CODIMENSION_1, title);

					int number;
					if (delta.nonzero())
					{
						number = CODIMENSION_1;
						display_delta = true;
					}
					else 
					{
						cout << "\nDelta_1^R=0";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "\n\npoly_invariant: Delta_1^R=0";

						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "Delta_1^R=0";
						}

						* take out two rows and columns from R_study_matrixptr at a time *
						delta = delta_1(R_study_matrix_rep, CODIMENSION_2, title);

						if (delta == T("1"))
						{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "\n\npoly_invariant: Delta_2^R=1\n";	

							cout << "\nDelta_2^R=1";
							if (!braid_control::RAW_OUTPUT)
							{
								output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
								output << "Delta_2^R=1";
							}

							cout << "\nOut of ideas, jumping off nearest cliff...";
						}
						else
						{
							number = CODIMENSION_2;
							display_delta = true;
						}
					}

					if (display_delta)
					{
						cout << "\nDelta_";
						cout << (number == 1? "1" : "2");
						cout << "^R = ";
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "\nDelta_";
							output << (number == 1? "1" : "2");
							output << "^R = ";
						}

						if (braid_control::OUTPUT_AS_INPUT)
							output << "\n";
		   				cout << delta << "\n";
		   				output << delta << "\n";

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "\npoly_invariant: Delta_";
	debug << (number == 1? "1" : "2");
	debug << "^R = " << delta << "\n";	
}
					}
				}
*****************************************************************/
			}
			else
			{
if (debug_control::DEBUG >= debug_control::SUMMARY && braid_control::QUATERNION)
	debug << "\npoly_invariant: evaluating codimension 1 determinant in M_{n,m}(H) using Study determinant\n";

				/* Use the normal procedure of taking the codimension 1 submatrices of the matrix 
				   representation, taking (Study) determinants and evaluating their gcd.
				   
					We include delta_0 as a generator for E_1 in all cases, though strictly it is 
					only required in the quaternionic (non-commutative) case but it keeps the code 
					simple and the overhead of this appraoach is small.
				*/	
				T delta = delta_1(*Matrix_rep, CODIMENSION_1,title,delta_0);
				
				if (braid_control::EXTRA_OUTPUT)
		    	{
					if (!braid_control::SILENT_OPERATION)
						cout << "\n";
					output << "\n";
		    	}
		    	else if (!braid_control::SILENT_OPERATION && matrix_control::WAIT_INFO)
		    	{
					/* we need a carridge return because of the wait info */
					cout << "\n";
		    	}

		    	if (braid_control::QUATERNION)
		    	{
					if (braid_control::NORMALIZING_Q_POLYNOMIALS)
						normalize(delta,switch_matrix);

					if (!braid_control::SILENT_OPERATION)
						cout << (braid_control::LONG_KNOT? "\np^{(1)} = ":"\nDelta_1^H = ");
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << (braid_control::LONG_KNOT? "p^{(1)} = ":"Delta_1^H = ");
					}
				}
				else
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "\nAlexander polynomial = ";
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Alexander polynomial = ";
					}
				}

				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
				if (!braid_control::SILENT_OPERATION)
					cout << delta << "\n";
		    	output << delta << "\n";

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "\n\npoly_invariant: " << (braid_control::LONG_KNOT? "p^{(1)} = ":"Delta_1 = ") << delta << "\n";	

			}

		}
		else
		{
		    if (braid_control::DISPLAY_DELTA_1_ONLY)
		    {
				/* We need a place marker in the output to identify why we're not displaying anything */
				if (!braid_control::SILENT_OPERATION)
					cout << "\nDelta_0 non zero, Delta_1 not calculated." << endl;
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Delta_0 non zero, Delta_1 not calculated." << endl;
				}
		    }
		}
    }
	
	delete Matrix_rep;
}

/* only a template because the compiler doesn't know this function is not going to be called by 
   matrices with non quaternionic coefficients
*/
template <typename T, typename St> 
T study_determinant (const matrix<T,St>& H_matrix, string title, int n=0, int* rperm=0, int* cperm=0)
{
	/* call C_study_matrix to create the image in M_{2n}(C) under \psi of the
	   n rows and columns of matrixptr indicated by rperm and cperm.  The image
	   will be stored under C_study_matrixptr.
	*/
	matrix<T,St> C_matrix = C_study_matrix(H_matrix, n, n, rperm, cperm);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
    debug << "\nstudy_determinant: Complex Study matrix:\n";
    debug << C_matrix;
    debug << "\n";
}

	return determinant(C_matrix,title);

}

/* C_study_matrix calculates the image under a mapping from M_{n,m}(H) to M_{2n,2m}(C)
   of the n rows and m columns of H_matptr indicated by rperm and cperm and stores
   it in a new Hpmatrix created at C_matptr.
   
   The mapping is either the mu-mapping, the default, or the rho-mapping, selected via the
   boolean variable STUDY_RHO_MAPPING.  The mu-mapping is based on the fact that any quaternion 
   w + xi +yj +zk may be written as Z + Wi, where Z = w + xi and W = y + zi are complex numbers.
   The rho-mapping is based on the fact that w + xi +yj +zk may also be written as Z + Wj, where 
   Z = w + zk and W = y + xk are complex numbers.
   
*/
matrix<polynomial<scalar,string,int>,scalar> C_study_matrix(const matrix<polynomial<scalar,string,int>,scalar>& H_matrix, int n=0, int m=0, int* rperm=0, int* cperm=0) 
{
	return H_matrix;
}

template <typename T, typename V, typename U, typename St> 
matrix<polynomial<T,V,U>,St> C_study_matrix(const matrix<polynomial<T,V,U>,St>& H_matrix, int n, int m, int* rperm, int* cperm)
{

if (debug_control::DEBUG >= debug_control::BASIC)
{
	if (braid_control::STUDY_RHO_MAPPING)
	    debug << "\nC_study_matrix: using rho-mapping to create complex Study matrix:" << endl;
	else
	    debug << "\nC_study_matrix: using mu-mapping to create complex Study matrix:" << endl;
}
	bool clean_up_params = false;
	
	typedef int It;
		
	if (n == 0)
	{
		n = H_matrix.numrows();
		m = H_matrix.numcols();
		rperm = new int[n];
		for (int i=0; i< n; i++)
			rperm[i] = i;
		cperm = new int[m];
		for (int i=0; i< m; i++)
			cperm[i] = i;
		clean_up_params = true;
	}
	
    matrix<polynomial<T,V,U>,St> C_matrix(2*n,2*m);

    for (int i=0; i<n; i++)
    for (int j=0; j<m; j++)
    {
		/* take the quaternionic polynomial at H_matrix[rperm[i]][cperm[j]]
		   and create the 2x2 matrix      Z       W
								      -\bar{W} \bar{Z}
		   at C_matptr[2*i][2*j], where Z and W are given by either the mu-mapping
		   or the rho-mapping, dependent on the value of STUDY_RHO_MAPPING.

		   The \bar function is complex conjugation.
		   
		   The mu-mapping takes (r,i,j,k) to   (r,i)   (j,k)
		                                       (-j,k)   (r,-i)
									
		   The rho-mapping takes (r,i,j,k) to   (r,k)   (j,i)
		                                       (-j,i)   (r,-k)
											  
		   That is, i and k are interchanged.
		*/
		C_matrix[2*i][2*j] = H_matrix[rperm[i]][cperm[j]];
		
		/* Set the 00 position in relation to the  base, 2*i,2*j.  Moving along the pterms,
		   for the mu-mapping set the j and k  components of the coefficients to zero
		   for the rho-mapping first copy the k component to the i component, then set j and k to zero.
		*/
		pterm<T,U>* pterm_ptr = C_matrix[2*i][2*j].get_pterms();
		
    	while(pterm_ptr)
    	{
			if (braid_control::STUDY_RHO_MAPPING)
				pterm_ptr->n.Qi() = pterm_ptr->n.Qk();

			pterm_ptr->n.Qj() = It(0);
			pterm_ptr->n.Qk() = It(0);
			pterm_ptr = pterm_ptr->next;
    	}
    	sanitize(&C_matrix[2*i][2*j]);

		/* now copy the Z we have just created to the
		   11 relative position and create the conjugate
		*/
		C_matrix[2*i+1][2*j+1] = C_matrix[2*i][2*j];

		pterm_ptr = C_matrix[2*i+1][2*j+1].get_pterms();

	    while(pterm_ptr)
    	{
			pterm_ptr->n.Qi() *= It(-1);
			pterm_ptr = pterm_ptr->next;
    	}

		/* next create W at the relative 01 position */
		C_matrix[2*i][2*j+1] = H_matrix[rperm[i]][cperm[j]];

		/* For the mu-mapping copy j and k to r and i, then set j and k to zero.
		   For the rho-mapping just copy j to r and then set j and k to zero (i is correctly set in this case)
		*/
		pterm_ptr = C_matrix[2*i][2*j+1].get_pterms();
	    while(pterm_ptr)
	    {
			pterm_ptr->n.Qr() = pterm_ptr->n.Qj();
			
			if (!braid_control::STUDY_RHO_MAPPING)
				pterm_ptr->n.Qi() = pterm_ptr->n.Qk();
				
			pterm_ptr->n.Qj() = It(0);
			pterm_ptr->n.Qk() = It(0);
			pterm_ptr = pterm_ptr->next;
	    }
	    sanitize(&C_matrix[2*i][2*j+1]);

		/* finally copy the W we have just created to the
		   10 relative position and create the negative
		   conjugate - this ammounts to negating the r term
		*/
		C_matrix[2*i+1][2*j] = C_matrix[2*i][2*j+1];
		pterm_ptr = C_matrix[2*i+1][2*j].get_pterms();
	    while(pterm_ptr)
    	{
			pterm_ptr->n.Qr() *= It(-1);
			pterm_ptr = pterm_ptr->next;
    	}
    }

	if (clean_up_params)
	{
		delete[] rperm;
		delete[] cperm;
	}
	
	return C_matrix;
	
}


/* R_study_matrix calculates the image of C_matrix under the mapping from M_n(C)
   to M_{2n}(R).
*/
matrix<polynomial<scalar,int>,scalar> R_study_matrix(const matrix<polynomial<scalar,int>,scalar>& C_matrix) {return C_matrix;}

template <typename T, typename V, typename U, typename St> 
matrix<polynomial<T,V,U>,St> R_study_matrix(const matrix<polynomial<T,V,U>,St>& C_matrix)
{
	int n = C_matrix.numcols();

	typedef int It;

    matrix<polynomial<T,V,U>, St> R_matrix = matrix<polynomial<T,V,U>,St>(2*n,2*n);

    for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
    {
		/* take the complex polynomial at C_matrix[i][j]
		   and create the 2x2 matrix   X   Y
								      -Y   X
		   at R_matptr[2*i][2*j].
		
		   Here X is the polynomial formed by setting all complex terms in 
		   coeficients of C_matrix[i][j] to zero, and Y is formed by
		   setting all the real components to the respective values of the 
		   complex terms, then setting the complex terms to zero.
		*/
		R_matrix[2*i][2*j] = C_matrix[i][j];
		
		/* move along the pterms setting the i components of the coefficients to
		   zero, the j and k components are already zero.  This is the 00
		   position in relation to the base, 2*i,2*j
		*/
		pterm<T,U>* pterm_ptr = R_matrix[2*i][2*j].get_pterms();
    	while(pterm_ptr)
    	{
			pterm_ptr->n.Qi() = It(0);
			pterm_ptr = pterm_ptr->next;
    	}
    	sanitize(&R_matrix[2*i][2*j]);

		/* now copy the X we have just created to the
		   11 relative position
		*/
		R_matrix[2*i+1][2*j+1] = R_matrix[2*i][2*j];

		/* next create Y at the relative 01 position */
		R_matrix[2*i][2*j+1] = C_matrix[i][j];

		/* copy i to r then delete i */
		pterm_ptr = R_matrix[2*i][2*j+1].get_pterms();
	    while(pterm_ptr)
	    {
			pterm_ptr->n.Qr() = pterm_ptr->n.Qi();
			pterm_ptr->n.Qi() = It(0);
			pterm_ptr = pterm_ptr->next;
	    }
	    sanitize(&R_matrix[2*i][2*j+1]);
		
		/* finally copy the Y we have just created to the
		   10 relative position and negate the r term
		*/
		R_matrix[2*i+1][2*j] = R_matrix[2*i][2*j+1];
		pterm_ptr = R_matrix[2*i+1][2*j].get_pterms();
	    while(pterm_ptr)
    	{
			pterm_ptr->n.Qr() *= It(-1);
			pterm_ptr = pterm_ptr->next;
    	}
    }
	
	return R_matrix;
}

/* delta_1 evaluates the gcd of the generators of the first elementary ideal, E_1,
   determined by the matrix M.  It is able to work with the intermediate
   ideals determined by mapping from M_{n,m}(H) into M_{2n,2m}(C) and evaluating codimension
   1 or 2 determinants of the image in M_{2n,2m}(C), providing M is given as 
   this image at the time of the call.  
*/
template <typename T, typename St> T delta_1(const matrix<T,St>& M, int codimension, string title, T delta_0 = T("0"))
{
	T hcf=delta_0; // needed for the non-commutative case where delta_0 is required as a generator of E_1
    Laurent_scale (hcf);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "delta_1: hcf initialized to " << hcf << endl;
	
	int matrix_rows = M.numrows();
	int matrix_cols = M.numcols();
	
    int rperm[matrix_rows];
    int cperm[matrix_cols];
	
	if (codimension == 1)
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "delta_1: codimension 1 generators of 1st elementary ideal, E_1:" << endl;

		for (int r1=0;r1<matrix_rows;r1++)
    	{
			if (!braid_control::SILENT_OPERATION && matrix_control::WAIT_INFO && matrix_control::COMFORT_DOTS)
			{
				cout << ".";
				cout.flush();
			}

			/* take out the r1-th row */
			for (int i=0; i<r1; i++)
			    rperm[i] = i;
			for (int i=r1+1; i<matrix_rows; i++)
			    rperm[i-1] = i;

			for (int c1=0;c1<matrix_cols;c1++)
			{
			    /* take out the c1-th column */
		    	for (int i=0; i<c1; i++)
				cperm[i] = i;
					
				if (matrix_cols > matrix_rows)
				{
					/* We're dealing with a long knot and there is one extra column than row */
					for (int c2= c1+1; c2< matrix_cols; c2++)
					{
						/* take out the c2-nd column */
						for (int i=c1+1; i<c2; i++)
						    cperm[i-1] = i;
					
			    		for (int i=c2+1; i<matrix_cols; i++)
							cperm[i-2] = i;
						
						if (minor_determinant(hcf, M, rperm, cperm, codimension, r1, 0, c1, c2, 0, title))
							goto delta_1_calculation_complete;
						
					}
				}
				else
				{				
			    	for (int i=c1+1; i<matrix_cols; i++)
					cperm[i-1] = i;

					if (minor_determinant(hcf, M, rperm, cperm, codimension, r1, 0, c1, 0, 0, title))
						goto delta_1_calculation_complete;
				}
			}
		}
	}
	else if (codimension == 2)
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "delta_1: codimension 2 generators of intermediate 1st elementary ideal, E_1:" << endl;
	
		for (int r1=0;r1<matrix_rows;r1++)
    	{

			if (!braid_control::SILENT_OPERATION && matrix_control::WAIT_INFO && matrix_control::COMFORT_DOTS)
			{
				cout << ".";
				cout.flush();
			}

			/* take out the r1-th row */
			for (int i=0; i<r1; i++)
			    rperm[i] = i;
			
			for (int r2= r1+1; r2< matrix_rows; r2++)
			{
				/* take out the r2-nd row */
				for (int i=r1+1; i<r2; i++)
				    rperm[i-1] = i;
				
				for (int i=r2+1; i<matrix_rows; i++)
			    	rperm[i-2] = i;

				for (int c1=0;c1<matrix_cols;c1++)
				{
				    /* take out the c1-th column */
			    	for (int i=0; i<c1; i++)
						cperm[i] = i;
					
					for (int c2= c1+1; c2< matrix_cols; c2++)
					{
						/* take out the c2-nd column */
						for (int i=c1+1; i<c2; i++)
						    cperm[i-1] = i;
					
						if (matrix_cols > matrix_rows)
						{
							/* We're dealing with a long knot and there is one extra column than row */
							for (int c3= c2+1; c3< matrix_cols; c3++)
							{
								/* take out the c3-rd column */
								for (int i=c2+1; i<c3; i++)
						    		cperm[i-2] = i;
					
			    				for (int i=c3+1; i<matrix_cols; i++)
									cperm[i-3] = i;
						
								if (minor_determinant(hcf, M, rperm, cperm, codimension, r1, r2, c1, c2, c3, title))
									goto delta_1_calculation_complete;
							}
						}
						else
						{
				    		for (int i=c2+1; i<matrix_cols; i++)
								cperm[i-2] = i;
						
							if (minor_determinant(hcf, M, rperm, cperm, codimension, r1, r2, c1, c2, 0, title))
								goto delta_1_calculation_complete;
						}

					}
				}
			}
		}
    }
	else
	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\nError! Unsupported codimension parameter presented to function delta_1" << endl;
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "delta_1: Error! Unsupported codimension parameter" << endl;
		exit(0);
	}

delta_1_calculation_complete:
	return hcf;
}

/* minor determinant is a macro-like function that is only used because delta_1
   needs to cope with codimension 1 and 2 cases for both square and non-square matrices, 
   so this function provides common code for all cases.  The parameters at the call
   are therefore essentially the same as those for delta_1.
   
   The function was extended in version 12.0 to return a boolean indicating whether or not the
   minor determinant it calculates is a unit.  If so then the function returns true, otherwise
   it returns false.
*/
template <typename T, typename St> bool minor_determinant(T& hcf, const matrix<T,St>& M, int* rperm, int* cperm,
                       int codimension, int r1, int r2, int c1, int c2, int c3, string title)
{
	/* rperm and cperm will always indicate a square submatrix of M but there may be more columns than rows
	   in M, so we use the number of rows as the matrix_size to pass to the determinant functions*/	
    int matrix_size = M.numrows();
	
	T delta;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "minor_determinant: matrix_size = " << matrix_size << ", codimension = " << codimension
          << ", r1=" << r1 << ", r2=" << r2 << ", c1=" << c1 << ", c2=" << c2 << ", c3=" << c3 << endl;;
    debug << "minor_determinant: rperm: ";
	for ( int k = 0 ; k< matrix_size-codimension; k++)
    	debug << rperm[k] << " ";
	debug << endl;
    debug << "minor_determinant: cperm: ";
	for ( int k = 0 ; k< matrix_size-codimension; k++)
    	debug << cperm[k] << " ";
	debug << endl;

}
	

   	if (!braid_control::QUATERNION && matrix_size == 2)
		delta = M[rperm[0]][cperm[0]];
   	else
   	{
//		if (!braid_control::QUATERNION || braid_control::COMPLEX_STUDY_DELTA_1 || braid_control::REAL_STUDY_DELTA_1)
		if (!braid_control::QUATERNION || braid_control::COMPLEX_STUDY_DELTA_1)
		{
			/* we want the ordinary determinant of the given matrix,
			   either because we are working with the Alexander rep or
			   because the matrix contains complex or real entries
			   determined in the course of evaluating the Study versions of Delta_1
			*/	
    			delta = determinant(M,title,matrix_size-codimension,rperm,cperm);
			
		}
		else
		{
			/* QUATERNION && !COMPLEX_STUDY_DELTA_1 && !REAL_STUDY_DELTA_1,
			   so take study_determinant, M will be quaternionic matrix rep in this case*/
    		delta = study_determinant(M,title,matrix_size-codimension,rperm,cperm);
		}
   	}

    if (braid_control::EXTRA_OUTPUT)
   	{
		if (codimension == 1)
		{
			if (!braid_control::SILENT_OPERATION)
			{
				cout << "\nGenerator (" << r1 << "," << c1;
				if (braid_control::LONG_KNOT)
					cout << "&" << c2;
				cout << ") = " << delta;
			}

			output << "\nGenerator (" << r1 << "," << c1;
			if (braid_control::LONG_KNOT)
				output << "&" << c2;
			output << ") = " << delta;
		}
		else
		{
			if (!braid_control::SILENT_OPERATION)
			{
				cout << "\nGenerator (" << r1 << "&" << r2 << "," << c1 << "&" << c2;
				if (braid_control::LONG_KNOT)
					cout << "&" << c3;
				cout << ") = " << delta;
			}

			output << "\nGenerator (" << r1 << "&" << r2 << "," << c1 << "&" << c2;
			if (braid_control::LONG_KNOT)
				output << "&" << c3;
			output << ") = " << delta;
		}
   	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	if (codimension == 2)
		debug << "minor_determinant: generator (" << r1 << " & " << r2 << "," << c1 << " & " << c2 << ") = " << delta << endl;
	else
    	debug << "minor_determinant: generator (" << r1 << "," << c1 << ") = " << delta << endl;
}
    if (!braid_control::QUATERNION)
   	{
		/* set s=1 in delta if it involves s */
		if (delta.numvars() && strchr(delta.getvars(),'s'))
		    set_to_one(delta,'s');

		if (braid_control::EXTRA_OUTPUT)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\nwith s=1: " << delta;
		    output << "\nwith s=1: " << delta;
		}
   	}

    /* Now delta is a Laurent polynomial in one variable,
       scale it so that it becomes an ordinary polynomial
    */
    Laurent_scale (delta);

   	if (braid_control::EXTRA_OUTPUT)
    {
		if (!braid_control::SILENT_OPERATION)
			cout << "  .=. " << delta;
		output << "  .=. " << delta;
   	}

    sanitize(&delta);
	
	if (braid_control::DELTA_1_UNIT_CHECK && delta.isunit())
	{
		hcf = T("1");
		
	   	if (braid_control::EXTRA_OUTPUT)
    	{
			if (!braid_control::SILENT_OPERATION)
				cout << "\nunit generator detected, terminating calculation" << endl;
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "unit generator detected, terminating calculation" << endl;
			}
   		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "minor_determinant: unit generator detected, terminating calculation" << endl;

		return true; // unit generator encountered
	}
	else if (delta != T("0"))
	{
	    hcf = gcd(hcf,delta);

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "minor_determinant: gcd stands at " << hcf << endl;

		if (braid_control::DELTA_1_UNIT_CHECK && hcf.isunit())
		{
		
		   	if (braid_control::EXTRA_OUTPUT)
	    	{
				if (!braid_control::SILENT_OPERATION)
					cout << "\nunit gcd detected, terminating calculation" << endl;
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "unit gcd detected, terminating calculation" << endl;
				}
	   		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "minor_determinant: unit gcd detected, terminating calculation" << endl;

			return true; // unit generator encountered
		}

	}
	else
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "minor_determinant: generator is zero, no change in gcd" << endl;
	}
	
	return false;
}

Qpmatrix adj2(Qpmatrix& M)
{
	Qpmatrix result(2,2);
	result[0][0] = M[1][1];
	result[0][1] = M[0][1]*Qpolynomial("-1");
	result[1][0] = M[1][0]*Qpolynomial("-1");
	result[1][1] = M[0][0];
	return result;
}

void simplify_matrix (Qpmatrix& M)
{
	for (size_t i = 0; i < M.numrows(); i++)
	for (size_t j = 0; j < M.numcols(); j++)
    	M[i][j] = simplify_rational(M[i][j]);
}



template <typename T, typename St> 
bool braid_rep(matrix<T,St>*& Matrix_rep, const matrix<T,St>& switch_matrix, const matrix<T,St>& switch_matrix_inverse, string input_string,
                int num_terms, int num_strings)
{
	bool calculate_polys;
	
	matrix<T,St>*		Twist;
	int switch_matrix_N_factor = switch_matrix.numcols()/2;
	char* inbuf = c_string(input_string);
	
	/* create the matrix corresponding to virtual crossings */
	if (switch_matrix_N_factor == 1)
	{
		Twist = new matrix<T,St>(2,2);
		matrix<T,St>& twist = *Twist;
		
		twist[0][0] = T("0");  //T is a polynomial type, so has a string constructor
		twist[0][1] = T("1");
		twist[1][0] = T("1");
		twist[1][1] = T("0");
	}
	else
	{
		Twist = new matrix<T,St>(2*switch_matrix_N_factor,2*switch_matrix_N_factor);
		matrix<T,St>& twist = *Twist;
	
		for (int i=0; i< switch_matrix_N_factor; i++)
		for (int j=0; j< switch_matrix_N_factor; j++)
			twist[i][j] = T("0");
		
		for (int i=0; i< switch_matrix_N_factor; i++)
		{
			for (int j=0; j<i; j++)
				twist[i][j+switch_matrix_N_factor] = T("0");
				
			twist[i][i+switch_matrix_N_factor] = T("1");
			
			for (int j=i+1; j<switch_matrix_N_factor; j++)
				twist[i][j+switch_matrix_N_factor] = T("0");
		}

		for (int i=0; i< switch_matrix_N_factor; i++)
		{
			for (int j=0; j<i; j++)
				twist[i+switch_matrix_N_factor][j] = T("0");
				
			twist[i+switch_matrix_N_factor][i] = T("1");
			
			for (int j=i+1; j<switch_matrix_N_factor; j++)
				twist[i+switch_matrix_N_factor][j] = T("0");
		}

		for (int i=0; i< switch_matrix_N_factor; i++)
		for (int j=0; j< switch_matrix_N_factor; j++)
			twist[i+switch_matrix_N_factor][j+switch_matrix_N_factor] = T("0");
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "\nbraid_rep: twist matrix:\n" << *Twist << endl;
}

	/* record whether this is a virtual knot or a classical one 
	if (strchr(inbuf, 't') || strchr(inbuf,'T'))
		virtual_crossings_present = true;
	*/
	
	int& N = switch_matrix_N_factor;

	int identity_size;
	
	/*	The doodle Alexander polynomial uses a representation of the twin group based on matrices of size
		num_strings-1.  Note that num_strings >=2.
	*/
	if (braid_control::DOODLE_ALEXANDER)
		identity_size = num_strings-1;	
	else
		identity_size = N*num_strings;	
		
	matrix<T,St> identity(identity_size,identity_size);	
	
//		matrix<T,St> identity(N*num_strings,N*num_strings);	
	
//	for (int i=0; i<N*num_strings; i++)
	for (int i=0; i<identity_size; i++)
	{
		for (int j=0; j<i; j++)
			identity[i][j] = T("0");
	
		identity[i][i] = T("1");
	
		for (int j=i+1; j<identity_size; j++)
//		for (int j=i+1; j<N*num_strings; j++)
			identity[i][j] = T("0");
	}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "\nbraid_rep: identity:\n" << identity;


    /* Create the matrix and initialize to the identity */
	if (braid_control::DOODLE_ALEXANDER)
		Matrix_rep = new matrix<T,St>(num_strings-1,num_strings-1);
	else
		Matrix_rep = new matrix<T,St>(N*num_strings,N*num_strings);
	matrix<T,St>& matrix_rep = *Matrix_rep;
	matrix_rep = identity;		

    /* Now evaluate the product of the matrices determined by
       the word in inbuf.
    */
	bool S;
	bool inverse;
    char* cptr = inbuf;
    for ( int i = 0; i< num_terms; i++)
    {
		if (*cptr == '-')
		{
		    inverse = true;
	    	cptr++;
		}
		else
		    inverse = false;

		if (*cptr == 's' || *cptr == 'S')
		    S = true;
		else
		    S = false;
		cptr++;

		char* mark = cptr; /* mark where we start the number */

		/* look for the end of the number */
		while (isdigit(*cptr))
		    cptr++;

		int number;
		get_number(number, mark);

		matrix<T,St> term = identity;

		const matrix<T,St>* Mptr;
		
		if (S && inverse)
		{
			Mptr= &switch_matrix_inverse;
		}
		else if (S)
		{
			Mptr = &switch_matrix;
		}
		else /* use T */
		{
			Mptr = Twist;
		}
				
		/* determine which form of Sn we are using Sn_matrix_top_down == true iff
           Sn = I^{n-1} x S x I^{k-n-1}, otherwise Sn = I^{k-n-1} x S x I^{n-1}
		*/
		size_t base;
		if (braid_control::Sn_matrix_top_down)
			base = (number-1)*N; // offsets from 0, braid strings from 1
		else
			base = term.numcols() - switch_matrix.numcols() - (number-1)*N;  // 
			                                                                      
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "\nbraid_rep: base for S in term matrix is " << base << endl;					
	
		/* If we're calculating the doodle Alexander polynomial, set the term 
		   explicitly, otherwise set it from Mptr;
		*/
																	  																					   
	    if (braid_control::DOODLE_ALEXANDER)
		{
			if (number == 1)
			{
				term[0][0] = T("-1");
				term[1][0] = T("x");
			}
			else if (number == num_strings-1)
			{
				/* If num_strings == 2, number will always be 1.  Therefore here num_strings >= 3 and number >=2 
				   term is a (num_strings-1,num_strings-1) matrix, so indices up to num_strings-2
				*/
				term[num_strings-3][num_strings-2] = T("x");
				term[num_strings-2][num_strings-2] = T("-1");
			}
			else // 2 <= number < num_strings-1
			{
				term[number-2][number-1] = T("x");
				term[number-1][number-1] = T("-1");
				term[number][number-1] = T("x");
			}
		}
		else
		{
			size_t switch_size = switch_matrix.numcols();		
			for (size_t j=0; j < switch_size; j++)
			for (size_t k=0; k < switch_size; k++)
			{
				term[base+j][base+k] = (*Mptr)[j][k];
			}
		}
		
		if (!braid_control::SILENT_OPERATION && braid_control::WAIT_SWITCH)
			cout << "term " << i+1 << endl;
			
		matrix_rep *= term;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "\nbraid_rep: term for right multiplication by ";
    if (inverse)
        debug << "-";

    if (S)
        debug << "s";
    else
        debug << "t";
    debug << number << ":\n";
    debug << term << endl;

}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
    debug << "\nbraid_rep: matrix after right multiplication by ";
    if (inverse)
        debug << "-";

    if (S)
        debug << "s";
    else
        debug << "t";
    debug << number << ":\n";
    debug << matrix_rep << endl;
}

	}

	if (braid_control::EXTRA_OUTPUT)
	{
		if (!braid_control::SILENT_OPERATION)
		{
			cout << "\nMatrix representation of braid word, M:\n";
			cout << matrix_rep << endl;
		}
		output << "\nMatrix representation of braid word:\n";
		output << matrix_rep << endl;
	}

if (debug_control::DEBUG >= debug_control::BASIC)
{
    debug << "\nbraid_rep: matrix representation of braid word, M:\n";
    debug << matrix_rep << endl;
}

	if (matrix_rep == identity)
	{
    	calculate_polys = false;
		
		if (!braid_control::SILENT_OPERATION)
			cout << "\nmatrix representation of braid word is the identity.\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << "matrix representation of braid word is the identity.\n";
		}
	}
	else
	{
		matrix_rep -= identity;
    	calculate_polys = true;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "\nbraid_rep: matrix representation minus the identity, M-I:\n";
    debug << matrix_rep << endl;
}
	}

	delete Twist;	
	delete[] inbuf;

	return calculate_polys;
}

template <typename T, typename St> 
void R_module_rep(matrix<T,St>*& Matrix_rep, const matrix<T,St>& switch_matrix, const matrix<T,St>& switch_matrix_inverse, 
               string input_string)
{
	bool OLD_SWITCH_ACTION = false;

	generic_code_data code_data;
	
	if (input_string.find('[') != string::npos)
		read_peer_code(code_data, input_string);
	else
		read_immersion_code(code_data, input_string);
	
	matrix<int>& code_table = code_data.code_table; 
	int num_crossings = code_data.num_crossings;
	int num_components = code_data.num_components;
	int num_classical = num_crossings;

	for (int i=0; i<num_crossings; i++)
	{
		if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::VIRTUAL)
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
		Matrix_rep = 0;
		return;
	}

	/* Determine the size of the R_module representation.  In the case of a knot or link, this is a 
	   2*N*num_classical square matrix where the switch_matrix elements are NxN matrices (for the 
	   Alexander and Quaternionic reps, N=1).  In the case of a long knot, the R_module representation
	   is a 2*N*num_classical x (2*N*num_classical + N) matrix, allowing for the extra variable.  
	*/
	int N = switch_matrix.numrows()/2; // N is used later, it is what is referred to elsewhere as the switch_matrix_N_factor
	int matrix_N_size = 2*num_classical;  // number of N-rows
	int matrix_size = N*matrix_N_size;  // number of actual rows

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "R_module_rep: num_components = " << num_components << endl;
	debug << "R_module_rep: num_crossings = " << num_crossings << endl;
	debug << "R_module_rep: num_classical = " << num_classical << endl;
	debug << "R_module_rep: N (switch_matrix_N_factor) = " << N << endl;
	debug << "R_module_rep: matrix_N_size = " << matrix_N_size << endl;
	debug << "R_module_rep: matrix_size = " << matrix_size << endl;
}

	/* The matrix representation changed in version 12 of the braid programme to a more natural, but essentially equivalent,
	   formulation.  The code supports both the old and new format via the boolean variable OLD_SWITCH_ACTION to invoke the 
	   original action of a switch at a real crossing.
	   
	   The original switch action was determined from the equations
	
		A  B  x_2i    = x_2j         A  B  x_{2j-1}  = x_{2i+1}   
		C  D  x_{2j-1}  x_{2i+1}     C  D  x_2i        x_2j 
	
	   for type I and type II crossings respectively, where the matrix A B
																	   C D
	
	   is either S, S' or T, depending upon the label of the crossing.

	   This means that the action of the switch is such that for a positive crossing with labels

                    b d
				    a c

	   we have c=Aa+Bb d=Ca+Db, which is why the above equations have the  column vector (x_2i, x_{2j-1}) the "wrong way around" 
	   for the corresponding crossing Type.

	   Thus for the old switch action we have equations
	
	   A * x_2i + B * x_{2j-1} - x_2j = 0 and 
	   C * x_2i + D * x_{2j-1} - x_{2i+1} = 0
	
	   for type I crossings, and for type II crossings

	   A * x_{2j-1} + B * x_2i - x_{2i+1} = 0 and 
	   C * x_{2j-1} + D * x_2i - x_2j = 0
	   

	   The new switch action was determined from the equations
	
		A  B  x_{2j-1}  = x_{2i+1}         A  B  x_2i      = x_2j   
		C  D  x_2i        x_2j             C  D  x_{2j-1}    x_{2i+1} 
	
	   for type I and type II crossings respectively, using the same notation as above.  Here the switch action is more natural
	   in the sense that the action is determined by 'standard' matrix multiplication of a vector.  
	   
	   Thus for the new switch action we have equations
	
	   A * x_{2j-1} + B * x_2i - x_{2i+1} = 0 and 
	   C * x_{2j-1} + D * x_2i - x_2j = 0
	
	   for type I crossings, and for type II crossings

	   A * x_2i + B * x_{2j-1} - x_2j = 0 and 
	   C * x_2i + D * x_{2j-1} - x_{2i+1} = 0

	   Note that in effect the change in the switch action is just inverting the definition of a type I and type II crossing.  This is 
	   how the code handles the difference, governed by the OLD_SWITCH_ACTION variable, which may then be controlled by the command
	   line at a later date, if required.


	   This reconciles with the theory in the paper by noting that the theory ignores virtual crossings, so here the above equations
	   state that x_2i = x_{2i+1} and x_{2j-1} = x_2j at a virtual crossing.  To reduce the size of the determinant we need to deal with, we remove
	   any such variable identities expressed by virtual crossings are handled by using a set of trivial variable flags that indicate 
	   when a variable is identical to another, we then use a variable permutation to map to a smaller set of variables.
   
	   NB. In the general case the switch matrix may have N x N matrix entries A B C and D (where N is the switch_matrix_N_factor), 
	   in which case the variables described here are N-tuples.
	*/

	valarray<int> trivial_variable_flags(2*num_crossings); 
	valarray<int> variable(2*num_crossings);               
	
	
	for (int i=0; i<2*num_crossings; i++)
		trivial_variable_flags[i] = 0;
	
	for (int i=0; i<num_crossings; i++)
	{
		if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::VIRTUAL)
		{
			/* The equations for this crossing state x_2i = x_{2i+1} and x_{2j-1} = x_2j, where the terminating
			   edges at the crossing are 2i and 2j-1.  Thus, the variables corresponding to the incoming edges
			   are trivial.
			*/
//			int twice_pi_minus1 = ((2*code_table[PERM][i] -1)+(2*num_crossings))%(2*num_crossings);
//if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
//	debug << "R_module_rep: crossing " << i << " trivial variables = " << 2*i << " and " << twice_pi_minus1 << endl;
//			trivial_variable_flags[2*i] = 1;
//			trivial_variable_flags[twice_pi_minus1] = 1;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "R_module_rep: crossing " << i << " trivial variables = " << 2*i << " and " << code_table[generic_code_data::table::OPEER][i] << endl;
			trivial_variable_flags[2*i] = 1;
			trivial_variable_flags[code_table[generic_code_data::table::OPEER][i]] = 1;
		}
	}
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "R_module_rep: trivial_variable_flags: ";
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
	
//	for (int i=0; i<2*num_crossings; i++)
//	{
//		variable[i] = i-count;
//		count += trivial_variable_flags[i];
//	}
	
//if (debug_control::DEBUG >= debug_control::DETAIL)
//{
//	debug << "R_module_rep: intermediate variables: ";
//	for (int i=0; i<2*num_crossings; i++)
//		debug << variable[i] << ' ';
//	debug << endl;
//}	

//	if (trivial_variable_flags[2*num_crossings-1])
//	{
		/* x_{2num_crossings-1} = x_{2num_crossings} = x_0.  Since we may have also had x_{2num_crossings-2} = 
		   x_{2num_crossings-1} etc.  We may avoid looking back around the knot for the last real crossing and 
		   simply reduce all variable values modulo 2*num_classical.
		*/
//		for (int i=0; i<2*num_crossings; i++)
//		{
//			variable[i] %= 2*num_classical;
//		}
//	}
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "R_module_rep: variables: ";
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

	if (braid_control::LONG_KNOT)
	    Matrix_rep = new matrix<T,St>(matrix_size,matrix_size+N);
	else
	    Matrix_rep = new matrix<T,St>(matrix_size,matrix_size);
		
	matrix<T,St>& matrix_rep = *Matrix_rep;
	
	for (unsigned int i=0; i<matrix_rep.numrows(); i++)
	for (unsigned int j=0; j<matrix_rep.numcols(); j++)
		matrix_rep[i][j] = T("0");
	

	int row = 0;
	
	for (int i=0; i<num_crossings; i++)
	{
		/* only consider classical crossings */
		if (code_table[generic_code_data::table::LABEL][i] == generic_code_data::VIRTUAL)
			continue;
		
//		/* first evaluate 2p(i)-1 = 2j-1 = 2*code_table[PERM][i]-1 modulo 2*num_crossings */
//		int twice_pi_minus1 = ((2*code_table[PERM][i] -1)+(2*num_crossings))%(2*num_crossings);

		/* identify the odd_numbered_peer of the naming edge */
//		int odd_numbered_peer = code_table[generic_code_data::table::OPEER][i];
	
		/* now determine the braid crossing type, POSITIVE, NEGATIVE or VIRTUAL */
		int crossing_type;		
		if (  (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1 && code_table[generic_code_data::table::LABEL][i] == generic_code_data::POSITIVE)
				||(code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE2 && code_table[generic_code_data::table::LABEL][i] == generic_code_data::NEGATIVE)
   					)
			crossing_type = generic_braid_data::crossing_type::NEGATIVE;
		else
			crossing_type = generic_braid_data::crossing_type::POSITIVE;

if (debug_control::DEBUG >= debug_control::DETAIL)
{	
	debug << "R_module_rep: crossing " << i;
	if ( crossing_type == generic_braid_data::crossing_type::NEGATIVE )
		debug << " negative";
	else
		debug << " positive";

	debug <<  ", generic_code_data::table::TYPE " << (code_table[generic_code_data::table::TYPE][i] == generic_code_data::TYPE1? "I" : "II") << endl;
	debug << "R_module_rep: ic-variables\tp[i] = ";// << code_table[PERM][i];
//	debug << "\t2i = " << 2*i << "\t\t2p[i]-1 = " << twice_pi_minus1;
	debug << "\t2i = " << 2*i << "\t\t2j-1 = " << code_table[generic_code_data::table::ODD_TERMINATING][i];
	debug << "\t2i+1 = " << 2*i+1 << "\t2p[i] = " << code_table[generic_code_data::table::EVEN_ORIGINATING][i] << endl;
	
	debug << "R_module_rep: mapped to\t\t\t";//v(p[i]) = " << variable[code_table[PERM][i]];
	debug << "\tv(2i) = " << variable[2*i] << "\tv(2p[i]-1) = " << variable[code_table[generic_code_data::table::ODD_TERMINATING][i]];
	debug << "\tv(2i+1) = " << variable[2*i+1] << "\tv(2p[i]) = " << variable[code_table[generic_code_data::table::EVEN_ORIGINATING][i]] << endl;

	debug << "R_module_rep: row " << 2*row << "\nR_module_rep: \t\t";
}		
		/* first do row 2*row.  2*num_classical is the number of N-cols in the matrix_rep for a closed knot, 
		   there is an additional N-col in the case of a long knot but that is handled explicity below.
		*/
		for (int j=0; j<2*num_classical; j++)  
    	{
//			if (j == (code_table[generic_code_data::table::TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[2*i] : variable[twice_pi_minus1]))
			if (j == (code_table[generic_code_data::table::TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[code_table[generic_code_data::table::EVEN_TERMINATING][i]] : variable[code_table[generic_code_data::table::ODD_TERMINATING][i]]))
			{
				/* Record the A entry at j*/
				if (crossing_type == generic_braid_data::crossing_type::POSITIVE)
				{
if (debug_control::DEBUG >= debug_control::DETAIL) debug << j << ":S00    ";
					set_matrix_N_element(matrix_rep,2*row,j,switch_matrix,0,0,N);
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL) debug << j << ":iS00   ";
					set_matrix_N_element(matrix_rep,2*row,j,switch_matrix_inverse,0,0,N);				
				}
			}
//			else if (j == (code_table[generic_code_data::table::TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[twice_pi_minus1] : variable[2*i]))
			else if (j == (code_table[generic_code_data::table::TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[code_table[generic_code_data::table::ODD_TERMINATING][i]] : variable[code_table[generic_code_data::table::EVEN_TERMINATING][i]]))
			{
				/* Record the B entry at j */
				if (crossing_type == generic_braid_data::crossing_type::POSITIVE)
				{
if (debug_control::DEBUG >= debug_control::DETAIL) debug << j << ":S01    ";
					set_matrix_N_element(matrix_rep,2*row,j,switch_matrix,0,1,N);				
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL) debug << j << ":iS01   ";
					set_matrix_N_element(matrix_rep,2*row,j,switch_matrix_inverse,0,1,N);				
				}
			}
			else 
			{
				/* there is a possibility that we have to handle the case where
				   
					A  B  x_p  = x_p     A  B  x_p  = x_s   
					C  D  x_q    x_r     C  D  x_q    x_q 
				
				   due to a Reidemeister Type 1 portion of the knot encountering
				   two or more virtual crossings.  Such a case is theoretically
				   redundant but cannot be detected algorithmically until we reach 
				   here.
				   
				   We deal with this case by adding the -1 entry to the row
				   rather than setting an entry of the matrix to -1, hence the 
				   following code...
				*/
				
//				if (j == (code_table[generic_code_data::table::TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[2*code_table[PERM][i]] : variable[2*i+1]))
				if (j == (code_table[generic_code_data::table::TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[code_table[generic_code_data::table::EVEN_ORIGINATING][i]] : variable[code_table[generic_code_data::table::ODD_ORIGINATING][i]]))
				{
					/* The theoretical edge variables ignore virtual crossings and have been represented by variable[i] 
					   so that all immersion-semi-arcs after the last real crossing have variable[i]=0.  If this is a
					   long knot then infinity sits on this arc and at the end of the long knot we encounter 
					   variable[i]=0 uniquely at the last real crossing.  Depending on how many virtual crossings there 
					   are left to cross before we head off to infinity, we could find it on the 2*row or 2*row+1 of
					   either a type 1 or type 2 crossing but if it's on 2*row then it will be in this conditional
					   clause when j==0.
					*/
					if (braid_control::LONG_KNOT && j == 0) 
					{
						/* Add -1 to the polynomial at (*matrixptr)[2*row][matrix_N_size] */
if (debug_control::DEBUG >= debug_control::DETAIL) debug << matrix_N_size << ":-1     ";
						decrement_matrix_N_element(matrix_rep,2*row,matrix_N_size,N);
					}
					else
					{
						/* Add -1 to the polynomial at (*matrixptr)[2*row][j] */
if (debug_control::DEBUG >= debug_control::DETAIL) debug << j << ":-1     ";
						decrement_matrix_N_element(matrix_rep,2*row,j,N);
					}
				}
				else
				{
					/* We've already set the location to zero, so record it
					   in the debug file if necessary */
if (debug_control::DEBUG >= debug_control::DETAIL) debug << j << ":0      ";
				}
			}
		}
   
if (debug_control::DEBUG >= debug_control::DETAIL) 
	debug << "\nR_module_rep: row " << 2*row+1 << "\nR_module_rep: \t\t";

		/* now row 2*row+1 */	
		for (int j=0; j<2*num_classical; j++)  
    	{
//			if (j == (code_table[generic_code_data::table::TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[2*i] : variable[twice_pi_minus1]))
			if (j == (code_table[generic_code_data::table::TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[code_table[generic_code_data::table::EVEN_TERMINATING][i]] : variable[code_table[generic_code_data::table::ODD_TERMINATING][i]]))
			{
				/* Record the C entry at j */
				if (crossing_type == generic_braid_data::crossing_type::POSITIVE)
				{
if (debug_control::DEBUG >= debug_control::DETAIL) debug << j << ":S10    ";
					set_matrix_N_element(matrix_rep,2*row+1,j,switch_matrix,1,0,N);				
				}
				else 
				{
if (debug_control::DEBUG >= debug_control::DETAIL) debug << j << ":iS10   ";
					set_matrix_N_element(matrix_rep,2*row+1,j,switch_matrix_inverse,1,0,N);				
				}
			}
//			else if (j == (code_table[generic_code_data::table::TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[twice_pi_minus1] : variable[2*i]))
			else if (j == (code_table[generic_code_data::table::TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[code_table[generic_code_data::table::ODD_TERMINATING][i]] : variable[code_table[generic_code_data::table::EVEN_TERMINATING][i]]))
			{
				/* Record the D entry at j */
				if (crossing_type == generic_braid_data::crossing_type::POSITIVE)
				{
if (debug_control::DEBUG >= debug_control::DETAIL) debug << j << ":S11    ";
					set_matrix_N_element(matrix_rep,2*row+1,j,switch_matrix,1,1,N);				
				}
				else 
				{
if (debug_control::DEBUG >= debug_control::DETAIL) debug << j << ":iS11   ";
					set_matrix_N_element(matrix_rep,2*row+1,j,switch_matrix_inverse,1,1,N);				
				}
			}
			else 
			{				
//				if (j == (code_table[generic_code_data::table::TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[2*i+1] : variable[2*code_table[PERM][i]]))
				if (j == (code_table[generic_code_data::table::TYPE][i] == (OLD_SWITCH_ACTION?generic_code_data::TYPE1:generic_code_data::TYPE2) ? variable[code_table[generic_code_data::table::ODD_ORIGINATING][i]] : variable[code_table[generic_code_data::table::EVEN_ORIGINATING][i]]))
				{
					/* see comment for corresponding clause above */
					if (braid_control::LONG_KNOT && j == 0)
					{
						/* Add -1 to the polynomial at (*matrixptr)[2*row+1][matrix_N_size] */
if (debug_control::DEBUG >= debug_control::DETAIL) debug << matrix_N_size << ":-1     ";
						decrement_matrix_N_element(matrix_rep,2*row+1,matrix_N_size,N);
					}
					else
					{
						/* Add -1 to the polynomial at (*matrixptr)[2*row+1][j] */
if (debug_control::DEBUG >= debug_control::DETAIL) debug << j << ":-1     ";
						decrement_matrix_N_element(matrix_rep,2*row+1,j,N);
					}
				}
				else
				{
					/* We've already set the location to zero, so record it
					   in the debug file if necessary */
if (debug_control::DEBUG >= debug_control::DETAIL) debug << j << ":0      ";
				}
			}
		}

if (debug_control::DEBUG >= debug_control::DETAIL) 
	debug << endl;
		
		row++;
	}
	

	if (braid_control::EXTRA_OUTPUT)
	{
		if (!braid_control::SILENT_OPERATION)
		{
			cout << "\nMatrix representation of R-module determined by immersion code, M:\n";
			cout << matrix_rep << endl;
		}
		output << "\nMatrix representation of R-module determined by immersion code, M:\n";
		output << matrix_rep << endl;
	}

if (debug_control::DEBUG >= debug_control::BASIC)
{
    debug << "R_module_rep: matrix representation derived from generic code data:" << endl;
    debug << matrix_rep;
    debug << endl;
}

}

/* rat_minor_determinant was introduced for rat_poly_invariant to deal with both long knots and closed knots.  
   It evaluates a minor determinant and if necessary updates the numerator gcd, returning true if it detects a
   unit generator or gcd, and false otherwise.
*/
bool rat_minor_determinant (Qpmatrix* Matrix_rep, int N, int Nr1, int Nc1, int Nc2, int* rperm, int* cperm, string title, polynomial<scalar,char,bigint>& hcf)
{
	int matrix_rows = Matrix_rep->numrows();
	Qpolynomial det = determinant(*Matrix_rep, title,  matrix_rows - N, rperm, cperm);

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "rat_minor_determinant: " << (braid_control::LONG_KNOT? "p^{(1)}":"Delta_1") << " generator " << N << "-row/column (" << Nr1 << "," << Nc1;
	if (braid_control::LONG_KNOT)
		debug << "&" << Nc2;
	debug << ")" << endl;
}
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "rat_minor_determinant: rperm: ";
    for (int k = 0 ; k< matrix_rows-N; k++)
		debug << rperm[k] << " ";
	debug << endl;
	debug << "rat_minor_determinant: cperm: ";
    for (int k = 0 ; k< matrix_rows-N; k++)
		debug << cperm[k] << " ";
	debug << endl;
}
					
	if (!braid_control::SILENT_OPERATION)
	{
		cout << (braid_control::LONG_KNOT? "\np^{(1)}":"\nDelta_1") << " generator (" << Nr1 << "," << Nc1;
		if (braid_control::LONG_KNOT)
			cout << "&" << Nc2;
		cout << ") = " << det;
	}
					
	if (!braid_control::RAW_OUTPUT)
	{
		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");

		output << (braid_control::LONG_KNOT? "p^{(1)}":"Delta_1") << " generator (" << Nr1 << "," << Nc1;
		if (braid_control::LONG_KNOT)
			output << "&" << Nc2;
		output << ") = " << det;
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "rat_minor_determinant: ";
	debug << (braid_control::LONG_KNOT? "p^{(1)}":"Delta_1") << " generator det = " << det << endl;
}

	polynomial<scalar,char,bigint> det_numerator = det.getn();
   	sanitize(&det_numerator);
	polynomial<scalar,char,bigint> det_denominator = det.getd();
   	sanitize(&det_denominator);

	if (braid_control::DELTA_1_UNIT_CHECK && det_denominator.isunit() && det_numerator.isunit())
	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nunit generator detected, terminating calculation" << endl;
		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << "unit generator detected, terminating calculation" << endl;
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "rat_minor_determinant: unit generator detected, terminating calculation" << endl;

		if (braid_control::NUMERATOR_GCD)
			hcf = polynomial<scalar,char,bigint>("1");					
	
		return true; 
	}
	else if (braid_control::NUMERATOR_GCD)
	{
		/******** Take hcf code for single variable examples *********/
    	hcf = gcd(hcf,det_numerator);

		if (!braid_control::SILENT_OPERATION)
			cout << "\n numerator gcd stands at " << hcf << endl;
		if (!braid_control::RAW_OUTPUT)
		{
			output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
	    	output << "numerator gcd stands at " << hcf << endl;
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
   	debug << "rat_minor_determinant:  numerator gcd stands at " << hcf << endl;
	
		if (braid_control::DELTA_1_UNIT_CHECK && hcf.isunit())
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "\nunit gcd detected, terminating calculation" << endl;
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "unit gcd detected, terminating calculation" << endl;
			}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "rat_minor_determinant: unit gcd detected, terminating calculation" << endl;

			return true;; 
		}
	}
	
	return false;
	
}


/* rat_poly_invariant calculates the codimension 0 and codimension 1 invariants determined by the supplied 
   switch matrix and its inverse for the braid or immersion code given in the input string.  
   It is used only where the elements of the switch matrix are of type rational<polynomial<T,V,U> > and not 
   polynomial<T,V,U>, i.e. Matrix or Weyl switches. Note that this means that the switch matrix N-factor 
   is greater than one in this function.
*/
void rat_poly_invariant(const Qpmatrix& switch_matrix, const Qpmatrix& switch_matrix_inverse,
                    string input_string, string title)
{
	int 	num_terms;
	int 	num_strings;
	bool    calculate_polys = true;
	bool	virtual_crossings_present = false;
	
	Qpmatrix*		Matrix_rep;
	
	if (title.length())
   	{
		if (!braid_control::SILENT_OPERATION)
			cout << "\n\n" << title << "\n" << input_string << endl;
		if (!braid_control::RAW_OUTPUT)
			output << "\n\n" << title << "\n" << input_string << endl;
   	}		

	if (!braid_control::RAW_OUTPUT)
	{
		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
		output << input_string;		
	}

	/* note the N-factor of the switch matrix */
	int switch_matrix_N_factor = switch_matrix.numcols()/2;
	
if(debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "rat_poly_invariant: switch_matrix_N_factor = " << switch_matrix_N_factor << endl;

	/* The two cases of the next conditional clause are in-line rather than function calls only because
	   the code appear only once
	*/
//	if (input_string.find('(') != string::npos)
	if (input_string.find('(') != string::npos || input_string.find('[') != string::npos)
	{
		/* record whether this is a virtual knot or a classical one */
		if (input_string.find('*') != string::npos)
			virtual_crossings_present = true;

		/* check for a long knot */
		if (input_string.find("L:") != string::npos)
		{
			braid_control::LONG_KNOT = true;

			input_string = parse_long_knot_input_string(input_string);
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "rat_poly_invariant: long knot: " << input_string << endl;
	
			if (input_string == "link")
			{				
				if (!braid_control::SILENT_OPERATION)
					cout << "\nError! long knot indicator provided for the peer code of a link" << endl;
			
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "Error! long knot indicator provided for the peer code of a link" << endl;
				}
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "rat_poly_invariant: Error! long knot indicator provided for the peer code of a link" << endl;
	
				return;
			}
		}
		else
			braid_control::LONG_KNOT = false;
	
		R_module_rep(Matrix_rep, switch_matrix, switch_matrix_inverse, input_string);
	}
	else if (valid_braid_input(input_string, num_terms, num_strings, braid_control::SILENT_OPERATION, braid_control::RAW_OUTPUT, braid_control::OUTPUT_AS_INPUT))
	{ 

		/* remove any unwanted qualifiers from the input string */
		string::size_type pos = input_string.find('{');
		if (pos != string::npos)
			input_string = input_string.substr(0,pos);

		/* record whether this is a virtual knot or a classical one */
		if (input_string.find('t') != string::npos || input_string.find('T') != string::npos)
			virtual_crossings_present = true;
			
		/* When using a braid representation, we do not calculate polynomials if the representation is the identity
		   matrix.  This check is done within braid_rep */
		calculate_polys = braid_rep(Matrix_rep, switch_matrix, switch_matrix_inverse,input_string, num_terms, num_strings);
	} 

	/* Under normal circumstances Matrix_rep is created by R_module_rep or braid_rep.  However, if R_module_rep
	   determines that an immersion code contains only virtual crossings it returns a zero Matrix_rep pointer. */
	if (Matrix_rep == 0)
		return;
	
	/* Note the size of the matrix representation, this is either square or in the case of long knots
	   has one more N-column than N-row, so we note the number of N-rows.  For elements of VB(n) this is the 
	   N * number of strings, but for R_module representations it is 2 * N * number of classical crossings, 
	   where N is the switch matrix N-factor.  Here, we also need to  know the number of N-rows in the 
	   matrix representation so we can remove N-rows and N-columns.
	*/
	int N = switch_matrix_N_factor;
	int matrix_rows = Matrix_rep->numrows();
	int matrix_cols = Matrix_rep->numcols();
	int matrix_N_rows = matrix_rows/N;
	int matrix_N_cols = matrix_cols/N;
	
	if (calculate_polys)
	{
		/* First determine the codimension 0 invariant.  In the case of a long knot we have to remove successive columns  
		   and take the determinant of remianing submatrices, producing a set of generators for the first elementary ideal.
		   If we have been told, via braid_control::NUMERATOR_GCD, that we may calculate the gcd of the numerators we do so.
		   
		   In the case of a closed knot the codimension 0 invariant is just the determinant of the matrix representation.
		   
		   We then determine whether to calculate delta_1 based on the codimension zero value and the global booleans
		*/
		bool calc_delta_1;
	
		if (braid_control::WAIT_SWITCH)
		{
		    matrix_control::WAIT_INFO = true;
		    matrix_control::wait_count = 0;

		    /* the wait_threshold comes from the command line */
			if (!braid_control::SILENT_OPERATION)
				cout << "\nCalculating polynomials, please wait\n";
		}
		else
		    matrix_control::WAIT_INFO = false;

		Qpolynomial delta_0 = polynomial<scalar,char,bigint> ("0"); // required in 'else' of this clause and in the following clause
		
		if (braid_control::LONG_KNOT)
		{	
			/* take out each N-column of the matrix representation in turn and take the determinant
			   of the rest, keeping a record of the gcd of the determinants calculated.
			*/
			int rperm[matrix_rows];
			int cperm[matrix_cols];
			bool unit_detected = false;

			for (int i = 0 ; i< matrix_rows; i++)
		    	rperm[i] = i;

			/* hcf is used only in single variable examples */
			polynomial<scalar,char,bigint> hcf = polynomial<scalar,char,bigint>("0");

			Qpolynomial op_0 = polynomial<scalar,char,bigint>("0"); // \def\op{\buildrel o \over p}
			Qpolynomial np_0 = polynomial<scalar,char,bigint>("0"); // \def\np{\buildrel n \over p}

			for (int i = 0 ; i< matrix_N_cols; i++)
			{
				/* take out N-column i */
				for (int j=0; j<N*i; j++)
   					cperm[j] = j;
				for (int j=N*i+N; j<matrix_cols; j++)
  						cperm[j-N] = j;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "rat_poly_invariant: removing N-column " << i << endl;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "rat_poly_invariant: cperm: ";
   	for (int j = 0 ; j< matrix_cols-N; j++)
		debug << cperm[j] << " ";
	debug << endl;
}
				Qpolynomial det =  determinant(*Matrix_rep, title, matrix_rows, rperm, cperm);			
							
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "p^{(0)} generator " << i << " = " << det;
				}

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "rat_poly_invariant: p_0 generator det " << i << " = " << det << endl;

				if (i == 0)
				{
					op_0 = det;

					if (braid_control::WAIT_SWITCH)
					    matrix_control::wait_count = 0;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "rat_poly_invariant: op_0  set to det " << i << " = " << op_0 << endl;

				}
				else if (i == matrix_N_cols-1)
				{
					np_0 = det;

					if (braid_control::WAIT_SWITCH)
					    matrix_control::wait_count = 0;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "rat_poly_invariant: np_0  set to det " << i << " = " << np_0 << endl;
	
					if (unit_detected)
						break;  //we've jumped to here, so we're now done.

				}

				if (!braid_control::SILENT_OPERATION)
					cout << "\np^{(0)} generator " << i << " = " << det;

				polynomial<scalar,char,bigint> det_numerator = det.getn();
		    	sanitize(&det_numerator);
				polynomial<scalar,char,bigint> det_denominator = det.getd();
		    	sanitize(&det_denominator);

				if (braid_control::DELTA_1_UNIT_CHECK && det_denominator.isunit() && det_numerator.isunit())
				{
					unit_detected = true; 
					if (i != matrix_N_cols-1)
						i = matrix_N_cols-2; // jump to last column to determine np_0, i will be incremented by loop

					hcf = polynomial<scalar,char,bigint>("1");

					if (!braid_control::SILENT_OPERATION)
						cout << "\nunit generator detected, terminating calculation" << endl;
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "unit generator detected, terminating calculation" << endl;
					}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "rat_poly_invariant: unit generator detected, terminating calculation" << endl;
				}
				else if (braid_control::NUMERATOR_GCD)
				{
					/* take hcf code for single variable examples */
					polynomial<scalar,char,bigint> generator = det.getn();
		    		sanitize(&generator);

		    		hcf = gcd(hcf,generator);
					
					if (!braid_control::SILENT_OPERATION)
						cout << "\n numerator gcd stands at " << hcf << endl;
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
		    			output << " numerator gcd stands at " << hcf << endl;
					}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "rat_poly_invariant:  numerator gcd stands at " << hcf << endl;

					if (braid_control::DELTA_1_UNIT_CHECK && hcf.isunit())
					{
						if (!braid_control::SILENT_OPERATION)
							cout << "\nunit gcd detected, terminating calculation" << endl;
						if (!braid_control::RAW_OUTPUT)
						{
							output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
							output << "unit gcd detected, terminating calculation" << endl;
						}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "rat_poly_invariant: unit gcd detected, terminating calculation" << endl;

						unit_detected = true; 
						if (i != matrix_N_cols-1)
							i = matrix_N_cols-2; // jump to last column to determine np_0, i will be incremented by loop
					}
				}
			}

			if (!braid_control::SILENT_OPERATION)
				cout << endl;

			if (braid_control::NUMERATOR_GCD)
			{
				if (!braid_control::SILENT_OPERATION)
					cout << (braid_control::LONG_KNOT? "p^{(0)} = ":"Delta_0 = ") << hcf << endl;	
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "numerator gcd  = ";
					if (braid_control::OUTPUT_AS_INPUT)
						output << "\n";
				}
				output << hcf << endl;
			}

			if (!braid_control::SILENT_OPERATION)
			{
				cout << "\n\\op^{(0)} = " << op_0 << endl;		
				cout << "\n\\np^{(0)} = " << np_0 << endl;		
   			}
   			
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
	   			output << "\\op^{(0)} = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << op_0 << endl;		
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
	   			output << "\\np^{(0)} = ";
				if (braid_control::OUTPUT_AS_INPUT)
					output << "\n";
			}
			output << np_0 << endl;		

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
   	debug << "rat_poly_invariant: p_0  = p^{(0)} = " << delta_0 << endl;
   	debug << "rat_poly_invariant: \\op^{(0)} = " << op_0 << endl;
   	debug << "rat_poly_invariant: \\np^{(0)} = " << np_0 << endl;
}

    		if (braid_control::NUMERATOR_GCD)
			{
				if (hcf == polynomial<scalar,char,bigint>("0"))
					calc_delta_1 = true;
				else
					calc_delta_1 = false;
			}
			else
				calc_delta_1 = true; // we may have only been able to produce generators for p^{(0)}, so calculate p^{(1)} regardless
			
		}
		else //closed knot
		{
			if (virtual_crossings_present || braid_control::VERIFY_DELTA_0)
			{		
				delta_0 = determinant (*Matrix_rep, title);

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

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "rat_poly_invariant: Delta_0 = " << delta_0 << endl;

				if (delta_0 == Qpolynomial("0"))
					calc_delta_1 = true;
				else
					calc_delta_1 = false;
			}
			else
			{
	    		calc_delta_1 = true;

				if (!braid_control::DISPLAY_DELTA_1_ONLY)
				{
					if (!braid_control::SILENT_OPERATION)
						cout << "\nInput contains no virtual crossings, so Delta_0 = 0" << endl;
					if (!braid_control::RAW_OUTPUT)
					{
						output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
						output << "Input contains no virtual crossings, so Delta_0 = ";
						if (braid_control::OUTPUT_AS_INPUT)
							output << "\n";
					}
					output << "0" << endl;

				}
			}
		}

		if (!braid_control::CALCULATE_DELTA_0_ONLY && (braid_control::ALWAYS_CALCULATE_DELTA_1 || calc_delta_1))
		{
			/* hcf is used only in single variable examples */
			polynomial<scalar,char,bigint> hcf = delta_0.getn();

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
    debug << "\nrat_poly_invariant: hcf initialized to " << hcf << endl;


			int rperm[matrix_rows];
			int cperm[matrix_cols];

			/* take out a (N) row and (N) column from the matrix representation 
			   and calculate the determinant of what is left
			*/
			for (int Nr1=0;Nr1<matrix_N_rows;Nr1++)
			{
				/* take out the N rows of the matrix rep starting from row Nr1*N row */
				for (int i=0; i<Nr1*N; i++)
	    			rperm[i] = i;
				for (int i = Nr1*N+N; i< matrix_rows; i++)
					rperm[i-N] = i;


				int Nc2; // declared here so in scope for debug output below;
				for (int Nc1= 0; Nc1< matrix_N_cols; Nc1++)
				{
				
					Qpolynomial det = polynomial<scalar,char,bigint> ("0");
					
					/* take out the N columns of the matrix rep starting from column N*j */
					for (int i=0; i<Nc1*N; i++)
		   				cperm[i] = i;

					if (braid_control::LONG_KNOT)
					{
						/* There is one extra N-column than N-row */
						for (Nc2= Nc1+1; Nc2< matrix_N_cols; Nc2++)
						{
							/* take out Nc2-nd N column*/
							for (int i=Nc1*N+N; i<Nc2*N; i++)
		   						cperm[i-N] = i;

							for (int i=Nc2*N+N; i<matrix_cols; i++)
		   						cperm[i-2*N] = i;
								
							if (rat_minor_determinant(Matrix_rep, N, Nr1, Nc1, Nc2, rperm, cperm, title, hcf))
								goto rat_poly_delta_1_calculation_complete;
						}
					}
					else
					{
						for (int i=Nc1*N+N; i<matrix_cols; i++)
	   						cperm[i-N] = i;
						
						if (rat_minor_determinant(Matrix_rep, N, Nr1, Nc1, Nc2, rperm, cperm, title, hcf))
							goto rat_poly_delta_1_calculation_complete;
					}
				}			
			}

rat_poly_delta_1_calculation_complete:

			if (!braid_control::SILENT_OPERATION)
				cout << endl;
			
			if (braid_control::NUMERATOR_GCD)
			{
				if (!braid_control::SILENT_OPERATION)
					cout << (braid_control::LONG_KNOT? "p^{(1)} = ":"Delta_1 = ") << hcf << endl;	
				
				if (!braid_control::RAW_OUTPUT)
				{
					output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
					output << "numerator gcd  = ";
					if (braid_control::OUTPUT_AS_INPUT)
						output << "\n";
				}
				output << hcf << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "rat_poly_invariant: " << (braid_control::LONG_KNOT? "p^{(1)} = ":"Delta_1 = ") << hcf << endl;	

			}
		}
	}
			
	delete Matrix_rep;
}

/* normalize takes a (Delta_i^H) Hpolynomial, Laurent scales it, ensures that the
   constant term is positive, then multiplies or divides by |A|^2, |B|^2, |C|^2, 
   |D|^2, or |Delta|^2 so that the constant term is as clost to 1 a possible.
   This function needs interactive help from the user!!!!!
*/
template <typename T1, typename U, typename T2, typename St>  void normalize(polynomial<T1,U>& poly, const matrix<T2,St>& switch_matrix) {}

void normalize(Hpolynomial& poly, const Hpmatrix& switch_matrix)
{
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "\nnormalize: normalizing " << poly << endl;

	if (poly.get_pterms() == NULL)
	{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "\nnormalize: doing nothing, poly is zero\n" << endl;
		return;
	}
	
	Laurent_scale(poly);
	pterm<Quaternion,int>* pterm_ptr = poly.get_pterms();

	if (pterm_ptr->n.Qr() < 0)
		poly *= Hpolynomial("-1");

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "\nnormalize: scaled poly = " << poly << endl;

	/* if the constant term is 1 already, there's nothing to do */
	if (pterm_ptr->n == Quaternion(1))
	{	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "\nnormalize: constant term is already 1, nothing to do" << endl;			
		return;
	}
	
	Quaternion A = switch_matrix[0][0].get_pterms()->n;
	Quaternion B = switch_matrix[0][1].get_pterms()->n;
	Quaternion C = switch_matrix[1][0].get_pterms()->n;
	Quaternion D = switch_matrix[1][1].get_pterms()->n;

	Hpolynomial factor[5];
	factor[0] = Hpolynomial(norm_sqrd(A));
	factor[1] = Hpolynomial(norm_sqrd(B));
	factor[2] = Hpolynomial(norm_sqrd(C));
	factor[3] = Hpolynomial(norm_sqrd(D));
	Quaternion delta = conj(B)*A*Quaternion(norm_sqrd(D)) - conj(D)*C*Quaternion(norm_sqrd(B));
	factor[4] = Hpolynomial(norm_sqrd(delta));
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "normalize: initial factors : ";
	for (int i=0; i<5; i++)
		debug << factor[i] << ' ';
	debug << endl;
}

	/* we're not interested in factors of 1, nor repeats */
	int num_factors = 5;
	for (int i=0; i<5; i++)
	{
		if (factor[i] == Hpolynomial("1"))
		{
			factor[i] = Hpolynomial("0");
			num_factors--;
		}
	}
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "normalize: zeroed factors : ";
	for (int i=0; i<5; i++)
		debug << factor[i] << ' ';
	debug << endl;
}

	if (num_factors == 0)
	{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "normalize: all factors are 1, we're done." << endl;
		return;
	}
	else
	{
		/* move factors left so that non repeated factor are flushed left in factor */
		for (int i=0; i<num_factors; i++)
		{
			/* find the next non-zero factor */
			for (int j=i; j<5; j++)
			{
				if (factor[j] != Hpolynomial("0"))
				{
					/* check it's not a duplicate */
					for (int k=0; k<i; k++)
					{
						if (factor[j] == factor[k])
						{
							factor[j] = Hpolynomial("0");
							num_factors--;
							break;
						}
					}
					
					/* if it's still non-zero, move it to the next location */
					if (factor[j] != Hpolynomial("0"))
					{
						factor[i] = factor[j];
					
						/* set factor[j] to zero unless we've not moved it at all! */
						if (j != i)
							factor[j] = Hpolynomial("0");
						
						break;
					}
				}
			}
		}		
	}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{	
	debug << "normalize: available factors ";
	for (int i=0; i<num_factors; i++)
		debug << i+1 << ": " << factor[i] << "    ";
	debug << endl;
}

	Hpolynomial cfactor = Hpolynomial("1"); // used for debugging
	int rfactor;
	do
	{
		cout << "\nPolynomial stands at " << poly << endl;	
		cout << "Available factors ";
		for (int i=0; i<num_factors; i++)
			cout << i+1 << ": " << factor[i] << "    ";
		cout << endl;
		
		cout << "\nEnter factor number you want to apply, zero if you're done";
		cout << "\nfactor: ";
		cin >> rfactor;
			
		if (rfactor)
		{
			Hpolynomial& pfactor = factor[rfactor-1];
			
			cout << "\nMultiply or divide? [m/d] ";
		
			char multiply_or_divide;
			cin >> multiply_or_divide;
			
			if (multiply_or_divide == 'd')
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "\nnormalize: dividing by ";
				cfactor /= pfactor;
				poly /= pfactor;
			}
			else
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "\nnormalize: multiplying by ";
				cfactor *= pfactor;
				poly *= pfactor;			
			}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << pfactor;
			
			pterm_ptr = poly.get_pterms(); // the pterms may have been changed
			if (pterm_ptr->n == Quaternion(1))
			{	
				cout << "\nConstant term now 1" << endl;
				rfactor = 0;
			}
		}
	} while (rfactor);

if(debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "\nnormalize: total multiplicative factor applied = " << cfactor;

	cout << "\nNormalized poly = " << poly << endl;	
}

template <typename T> void equality_test (T x, T y, bool A)
{
	if (x == y)
	{
		if (!braid_control::SILENT_OPERATION)
		{
			if (A)
				cout << "\nSwitch elements A and D are equal" << endl;
			else
				cout << "Switch elements B and C are equal" << endl;
		}
		
		if (!braid_control::RAW_OUTPUT)
		{
			if (A)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Switch elements A and D are equal" << endl;
			}
			else
			{
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Switch elements B and C are equal" << endl;
			}
		}

if (debug_control::DEBUG >= debug_control::BASIC)
{
	if (A)
		debug << "braid::main: switch elements A and D are equal" << endl;
	else
		debug << "braid::main: switch elements B and C are equal" << endl;
}
	}
	else
	{
		if (!braid_control::SILENT_OPERATION)
		{
			if (A)
				cout << "\nSwitch elements A and D are not equal" << endl;
			else
				cout << "Switch elements B and C are not equal" << endl;
		}
		
		if (!braid_control::RAW_OUTPUT)
		{
			if (A)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << "Switch elements A and D are not equal" << endl;
			}
			else
			{
				if (braid_control::OUTPUT_AS_INPUT)
					output << ';';
				output << "Switch elements B and C are not equal" << endl;
			}
		}

if (debug_control::DEBUG >= debug_control::BASIC)
{
	if (A)
		debug << "braid::main: switch elements A and D are not equal" << endl;
	else
		debug << "braid::main: switch elements B and C are not equal" << endl;	
}
	}
}


/***************  Functions required for setting up braid specific debug **************/

void set_main_debug_option_parameter(char* pptr, string option);
void set_main_debug_option(char* start, char* end)
{
	char  loc_buf[end-start+2];
	char* c1 = start;
	char* c2 = loc_buf;

	/* if both start and end are zero, display debug help information.
	   this has been included here in this manner so that each time a debug option
	   is added, the help will be updated (hopefully!)
	*/
	if (start == 0 && end == 0)
	{
		set_main_debug_option_parameter(0,"braid");
		cout << "\t\tvogel, boolean" << endl;
		return;
	}

	do
	{
		*c2++ = *c1++;
	} while (c1 <= end);
	
	*c2 = '\0';

	char* pptr = strchr(loc_buf,'{');
	if (pptr)
	{
		*pptr++ = '\0';
	}

	/* now, even if there are parameters, the first part of loc_buf 
	   is a C-string that identifies the option */
	
	if (!strcmp(loc_buf,"braid"))
	{
		debug_control::DEBUG = debug_control::SUMMARY;
		debug << "main::set_main_debug_option: setting debug option debug_control::DEBUG = debug_control::SUMMARY\n";		
		
		if (pptr)
		{
			/* check for any parameters */
			check_debug_option_parameters(pptr, "braid");
		}
	}
	else if (!strcmp(loc_buf,"vogel"))
	{
		braid_control::VOGEL_DEBUG = true;
		debug << "main::set_main_debug_option: setting debug option braid_control::VOGEL_DEBUG\n";		
	}

	debug.flush();

}

void set_main_debug_option_parameter(char* pptr, string option)
{
	if (option == "braid")
	{
		if (!pptr)
		{
			cout << "\t\tbraid{summary|1:basic|2:intermediate|3:detail|4}, integer: default 0=off, no parameters sets summary" << endl;
		}
		else
		{
			if (!strcmp(pptr,"summary") || !strcmp(pptr,"1") )
			{
				debug_control::DEBUG = debug_control::SUMMARY;
				debug << "main::set_main_debug_option_parameter: setting debug option debug_control::DEBUG = debug_control::BASIC\n";		
			}
			if (!strcmp(pptr,"basic") || !strcmp(pptr,"2") )
			{
				debug_control::DEBUG = debug_control::BASIC;
				debug << "main::set_main_debug_option_parameter: setting debug option debug_control::DEBUG = debug_control::BASIC\n";		
			}
			else if (!strcmp(pptr,"intermediate") || !strcmp(pptr,"3"))
			{
				debug_control::DEBUG = debug_control::INTERMEDIATE;
				debug << "main::set_main_debug_option_parameter: setting debug option debug_control::DEBUG = debug_control::INTERMEDIATE\n";		
			}
			else if (!strcmp(pptr,"detail") || !strcmp(pptr,"4"))
			{
				debug_control::DEBUG = debug_control::DETAIL;
				debug << "main::set_main_debug_option_parameter: setting debug option debug_control::DEBUG = debug_control::DETAIL\n";		
			}
			else if (!strcmp(pptr,"exhaustive") || !strcmp(pptr,"5"))
			{
				debug_control::DEBUG = debug_control::EXHAUSTIVE;
				debug << "main::set_main_debug_option_parameter: setting debug option debug_control::DEBUG = debug_control::EXHAUSTIVE\n";		
			}
		}
	}
}

void main_display_default_options()
{
	cout << "\t\tbraid{summary}" << endl;
}

/* The option '2' for the complex Study Delta_1 option clashed with p=2 or W=2, 
   so the concept of cli suboptions similar to debug parsing was introduced.
*/
void set_cli_suboptions(char* optr, string option)
{
	char* start_of_option = optr;
	char* cptr = start_of_option+1;
	
	if (*cptr != '{')
		return;

	if (!strchr(cptr+1,'}'))
	{
		cout << "\nError! Unterminated command line suboption, '}' missing." << endl;
		exit(0);
	}
	else
		cptr++;

	if (option == "c")
	{
		do
		{
			if (*cptr == '2')
				braid_control::ALWAYS_CALCULATE_CODIMENSION_2_DELTA_1 = true;

			cptr++;			
		} while (*cptr != '}');
	}
	
	/* remove the {...} we've just processed */
	strcpy(start_of_option+1,cptr+1);
}

void set_acceptable_input_type()
{
//	if (braid_control::FIXED_POINT_INVARIANT || braid_control::RACK_POLYNOMIAL)
	if (braid_control::FIXED_POINT_INVARIANT)
	{
		input_control::ACCEPT_MAP |= input_control::braid_word;

		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nReading braid words from input.\n\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\n\n";
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Reading braid words from input.\n\n";
		}
	}
	else if (braid_control::RACK_POLYNOMIAL)
	{
		input_control::ACCEPT_MAP |= input_control::braid_word;
		input_control::ACCEPT_MAP |= input_control::peer_code;

		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nReading braid words and labelled peer codes from input.\n\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\n\n";
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Reading braid words  and labelled peer codes from input.\n\n";
		}
	}
	else if (braid_control::SWITCH_POLYNOMIAL_INVARIANT)
	{
		input_control::ACCEPT_MAP |= input_control::immersion_code;
		input_control::ACCEPT_MAP |= input_control::peer_code;
		input_control::ACCEPT_MAP |= input_control::braid_word;

		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nReading braid words, labelled immersion codes and labelled peer codes from input.\n\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\n\n";
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Reading braid words, labelled immersion codes and labelled peer codes from input.\n\n";
		}
	}
	else if (braid_control::VOGEL_ALGORITHM || braid_control::SATELLITE || braid_control::RACK_POLYNOMIAL || braid_control::DOODLE_Q_POLYNOMIAL)
	{
		input_control::ACCEPT_MAP |= input_control::peer_code;

		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nReading labelled peer codes from input.\n\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\n\n";
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Reading labelled peer codes from input.\n\n";
		}
	}
	else if (braid_control::GAUSS_CODE || braid_control::PEER_CODE || braid_control::HAMILTONIAN)
	{
		input_control::ACCEPT_MAP |= input_control::immersion_code;
		input_control::ACCEPT_MAP |= input_control::peer_code;
		input_control::ACCEPT_MAP |= input_control::gauss_code;
		input_control::ACCEPT_MAP |= input_control::braid_word;
		input_control::ACCEPT_MAP |= input_control::planar_diagram;
		input_control::ACCEPT_MAP |= input_control::dowker_code;

		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nReading braid words, labelled immersion codes, labelled peer codes and planar diagrams from input.\n\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\n\n";
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Reading braid words, labelled immersion codes, labelled peer codes and planar diagrams from input.\n\n";
		}
	}
	else if (braid_control::DOWKER_CODE)
	{
		input_control::ACCEPT_MAP |= input_control::braid_word;
		input_control::ACCEPT_MAP |= input_control::immersion_code;
		input_control::ACCEPT_MAP |= input_control::peer_code;

		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nReading braid words and labelled immersion codes from input.\n\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\n\n";
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Reading braid words and labelled immersion codes from input.\n\n";
		}
	}
	else if (braid_control::KAUFFMAN_BRACKET || braid_control::JONES_POLYNOMIAL || braid_control::ARROW_POLYNOMIAL)
	{
		input_control::ACCEPT_MAP |= input_control::immersion_code;
		input_control::ACCEPT_MAP |= input_control::peer_code;
		input_control::ACCEPT_MAP |= input_control::gauss_code;  // only Gauss codes containing classical crossings are supported currently, not flat or doodle crossings
		input_control::ACCEPT_MAP |= input_control::planar_diagram;


		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nReading labelled immersion codes, labelled peer codes and Gauss codes from input.\n\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\n\n";
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Reading labelled immersion codes, labelled peer codes and Gauss codes from input.\n\n";
		}
	}
	else if (braid_control::KNOTOID_BRACKET || braid_control::PARITY_BRACKET || braid_control::PARITY_ARROW)
	{
		input_control::ACCEPT_MAP |= input_control::immersion_code;
		input_control::ACCEPT_MAP |= input_control::peer_code;
		input_control::ACCEPT_MAP |= input_control::gauss_code;

		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nReading labelled peer codes and Gauss codes from input.\n\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\n\n";
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Reading labelled peer codes and Gauss codes from input.\n\n";
		}
	}
	else if (braid_control::AFFINE_INDEX)
	{
		input_control::ACCEPT_MAP |= input_control::immersion_code;
		input_control::ACCEPT_MAP |= input_control::peer_code;
		input_control::ACCEPT_MAP |= input_control::gauss_code; // only Gauss codes for classical diagrams are supported currently, not flat or doodle diagrams
		input_control::ACCEPT_MAP |= input_control::braid_word;
		input_control::ACCEPT_MAP |= input_control::planar_diagram;		

		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nReading braid words, labelled immersion codes, labelled peer codes and Gauss codes from input.\n\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\n\n";
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Reading braid words, labelled immersion codes, labelled peer codes and Gauss codes from input.\n\n";
		}
	}
	else if (braid_control::MOCK_ALEXANDER)
	{
		input_control::ACCEPT_MAP |= input_control::peer_code;
		input_control::ACCEPT_MAP |= input_control::gauss_code;
		input_control::ACCEPT_MAP |= input_control::planar_diagram;		

		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nReading labelled peer codes and Gauss codes from input.\n\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\n\n";
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Reading labelled peer codes and Gauss codes from input.\n\n";
		}
	}
	else if (braid_control::PRIME_TEST)
	{
		input_control::ACCEPT_MAP |= input_control::peer_code;
		input_control::ACCEPT_MAP |= input_control::gauss_code;
		input_control::ACCEPT_MAP |= input_control::planar_diagram;		

		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nReading labelled peer codes and Gauss codes from input.\n\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\n\n";
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Reading labelled peer codes and Gauss codes from input.\n\n";
		}
	}
	else
	{
		input_control::ACCEPT_MAP |= input_control::braid_word;

		if (!braid_control::SILENT_OPERATION)
			cout << "\n\nReading braid words from input.\n\n";
		if (!braid_control::RAW_OUTPUT)
		{
			output << "\n\n";
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';
			output << "Reading braid words from input.\n\n";
		}
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_acceptable_input_type: input_control::ACCEPT_MAP set to 0x" << hex << input_control::ACCEPT_MAP << dec << endl;
}
