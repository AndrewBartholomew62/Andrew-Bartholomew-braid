/*************************************************************************************************
                  Birack homology and cohomology calculation

                     A. Bartholomew 2rd June, 2025
                  
void perm_print(matrix<int>& M, ostream& os, int w, vector<int>& rperm, string prefix)
void print_characteristic_function(ostream& os, vector<scalar>& generator, vector<int>& reverse_map, bool cohomology, int n)
void birack_homology_generators(generic_switch_data& switch_data, int k, bool cohomology=false, matrix<scalar>* _B=0)
void calculate_boundary(int n, vector<int>& tuple, matrix<scalar>& Delta, int domain_index, vector<int>& codomain_map, matrix<int>& twitch_u, matrix<int>& twitch_d)
void test_cohomology_generators(generic_switch_data& switch_data)
*************************************************************************************************/
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <climits>
#include <map>

using namespace std;

extern ifstream     input;
extern ofstream 	output;
extern ofstream     debug;

extern bool SIDEWAYS_SEARCH;


#include <util.h>
#include <scalar.h> // includes bigint.h and rational.h
#include <quaternion-scalar.h>
#include <polynomial.h>
#include <matrix.h>
#include <ctype.h>
#include <braid.h>
#include <homology.h>

void print_k_chain(ostream& os, const vector<scalar> k_chain, generic_switch_data& switch_data, int k);
void calculate_boundary(int n, vector<int>& tuple, matrix<scalar>& Delta, int domain_index, vector<int>& codomain_map, matrix<int>& twitch_u, matrix<int>& twitch_d);


void perm_print(matrix<int>& M, ostream& os, int w, vector<int>& rperm, string prefix)
{
	size_t r=M.numrows();
	size_t c=M.numcols();
	for (size_t i=0; i< r; i++)
	{
		os << prefix;
		for (size_t j=0; j< c; j++)
			os << setw(w) << M[rperm[i]][j];
		os << endl;
	}
}

/* print_characteristic_function writes the given generator as a linear combination of the tuples of dimension k that make up the basis
   in which generator is a set of coefficients.  The k-tuples are enumerated as numbers base n.
   
   The value n is the size of the underlying birack
*/
void print_characteristic_function(ostream& os, vector<scalar>& generator, vector<int>& reverse_map, bool cohomology, int n, int k)
{
	bool first = true;
	int initial_factor = 1;
	for (int j=0; j< k-1; j++)
		initial_factor *= n;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "print_characteristic_function: generator: ";
	for (size_t i=0; i< generator.size(); i++)
		debug << generator[i] << ' ';
	debug << "n = " << n << ", k = " << k << ", cohomology = " << cohomology << endl;
	debug << "print_characteristic_function: initial_factor = " << initial_factor << endl;
}
			
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "print_characteristic_function: generator ";
	
	ostringstream oss;
	
	for (size_t i=0; i< generator.size(); i++)
	{
		scalar coefficient = generator[i];
		
		if (coefficient != scalar(0))
		{			
			int k_tuple = (braid_control::BIRACK_HOMOLOGY?i:reverse_map[i]);				

//if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
//	debug << "print_characteristic_function: i = " << i << " coefficient = " << coefficient << ", k_tuple = " << k_tuple << endl;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << coefficient << "(tuple=" << k_tuple;

			vector<int> p(k);
			int factor = initial_factor;
			for (int j=0; j< k-1; j++)
			{
				p[j] = k_tuple/factor;
				k_tuple %= factor;
				factor /=n;
			}
			p[k-1] = k_tuple;
			
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	for (int j=0; j< k; j++)
		debug << "(, p" << j+1 << '=' << p[j];
	debug << ") " << endl;
}
			/* We need to accommodate mod-p scalars but write the non-zero element of Z_2 as 1 not -1 */
			if (coefficient == scalar(-1) && !(braid_control::CALCULATE_MOD_P && mod_p::get_p() == 2))
				oss << '-';  
			else if (coefficient < scalar(0))
				oss << coefficient;
			else
			{
				if (!first)
					oss << '+';			
		
				if (coefficient != scalar(1))
					oss << coefficient;
			}
			
			if (cohomology)
				oss << 'X';

			oss << '(';
			for (int j=0; j< k-1; j++)
				oss << p[j] << ',';
			oss << p[k-1] << ')';
			
			first = false;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << 0 << ' ';
		}
	}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "\nprint_characteristic_function: generates " << oss.str() << endl;
	
	os << oss.str();
}

/* create_chain_map creates a map from all k-tuples in generators {0,..,n} to non-degenerate k-tuples in chain_map, returning the
   number of non-degenerate tuples found
*/
int create_chain_map(vector<int>& chain_map, int n, int k)
{
	vector<int> tuple(k);
	int num_generators = 0;	
	
	bool complete = false;
	do
	{
		/* check for degeneracy */
		bool degenerate = false;
		for (int i=0; i< k-1; i++)
		{
			if (tuple[i] == tuple[i+1])
			{
				degenerate=true;
				break;
			}
		}
		
		if (!degenerate)
		{
			int factor = 1;
			int place = 0;
			for (int i=k-1; i>=0; i--)
			{
				place += tuple[i]*factor;
				factor *= n;
			}
			chain_map[place] = num_generators++;
		}
		
		/* is there another tuple?*/
		complete = true;
		for (int i=k-1; i>=0; i--)
		{
			if (tuple[i] < n-1)
			{
				tuple[i]++;
				for (int j=i+1; j< k; j++)
					tuple[j] = 0;
				complete = false;
				break;
			}
		}			
	} while (!complete);
		
	return num_generators;
}

void create_boundary_map(matrix<scalar>& Delta, vector<int>& chain_map, matrix<int>& twitch_u,matrix<int>& twitch_d, int n, int k)
{
	vector<int> tuple(k);
	int index = 0;
	
	bool complete = false;
	do
	{
		/* check for degeneracy */
		bool degenerate = false;

		if (!braid_control::BIRACK_HOMOLOGY)
		{
			for (int i=0; i< k-1; i++)
			{
				if (tuple[i] == tuple[i+1])
				{
					degenerate=true;
					break;
				}
			}
		}
		
		if (!degenerate)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)	
{
	debug << "create_boundary_map: partial_ " << k << " of tuple = (";
	for (int i=0; i< k; i++)
		debug << tuple[i] << ',';
	debug << "), index = " << index << endl;
}
				calculate_boundary(n,tuple,Delta,index,chain_map,twitch_u,twitch_d);
				index++;
		}
		
		/* is there another tuple?*/
		complete = true;
		for (int i=k-1; i>=0; i--)
		{
			if (tuple[i] < n-1)
			{
				tuple[i]++;
				for (int j=i+1; j< k; j++)
					tuple[j] = 0;
				complete = false;
				break;
			}
		}			
	} while (!complete);
}

/* birack_homology_generators calculates generators for either the k-th homology or cohomology of the finite switch provided, dependent on the value of the
   cohomology boolean.  If braid_control::BIRACK_HOMOLOGY is false, the default, the function determines the biquandle variant of the homology
   or cohomology, setting degenerate tuples to zero.  Otherwise the birack variant is calculated, taking into account the degenerate tuples.
*/
void birack_homology_generators(generic_switch_data& switch_data, int k, bool cohomology=false, matrix<scalar>* _B=0)
{
	if (braid_control::BIRACK_HOMOLOGY)
		braid_control::USE_BIGINT = true;

	/* determine whether the scalar coefficients are a field */
	bool field_coefficients = true;
	if (scalar::variant == scalar::scalar_variant::INT || scalar::variant == scalar::scalar_variant::BIGINT)
		field_coefficients = false;
		
	/* evaluate the sideways twitch map from the switch map in the switch_data */

	int n = switch_data.size;
		
	matrix<int> twitch_u(n,n);
	matrix<int> twitch_d(n,n);
	
	matrix<int> inv_D(n,n);
	
	for (int i=0; i< n; i++)
	for (int j=0; j< n; j++)
		inv_D[i][switch_data.Sd[i][j]] = j; // inverts the map D_i given by the row switch_data.Sd[i][]
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "birack_homology_generators: inv_D:" << endl;
	print (inv_D,debug,3,"colouring_invariant: ");
}
	
	/* The sideways twitch map of a switch is T(a,b) = (t_a(b),t^b(a)), i.e S_{-}^{+} (a,b) = (b_{a^{-1}}, a^{b_{a^{-1}}}), 
	   where ^ and _ are switch actions.
	*/
	for (int a=0; a< n; a++)
	for (int b=0; b< n; b++)
	{
		int x = inv_D[a][b];  // D_a^{-1}(b)

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "birack_homology_generators:   a = " << a << " b = " << b << " x = " << x << endl;
	
		twitch_d[a][b] = x;
		twitch_u[b][a] = switch_data.Su[x][a];
	}


if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "birack_homology_generators: twitch_u:" << endl;
	print (twitch_u,debug,3,"birack_homology_generators: ");
	debug << "birack_homology_generators: twitch_d:" << endl;
	print (twitch_d,debug,3,"birack_homology_generators: ");
	
}		

	list<string>& cocycles = switch_data.cocycle_string;
	
	if (!braid_control::SILENT_OPERATION)
		cout << "calculating homology generators, cohomology flag = " << cohomology << endl;

	int n_km1_power = 1;
	for (int i=0; i< k-1; i++)
		n_km1_power *= n;
		
	int n_k_power = n_km1_power*n;
	int n_kp1_power = n_k_power*n;
	
	int num_km1_generators=0; // k minus 1
	int num_k_generators=0;
	int num_kp1_generators=0; // k plus 1
	
	/* we wish to enumerate the non-degenerate n-tuples for n=3,4 which we do by enumerating all n-tuples 
	   and mapping the non-degenerate ones to an enumeration that allows us to count them at the same time.
	   
	   For the 3-tuples, we record the reverse mapping so we can recover the 3-tuple from the enumeration of the
	   non-degenerate tuples.
	*/
	vector<int> km1_chain_map(n_km1_power,-1);
	vector<int> k_chain_map(n_k_power,-1);
	vector<int> reverse_k_chain_map(n_k_power,-1); // only the first num_k_generators places will have values
	vector<int> kp1_chain_map(n*n_k_power,-1);
	
	/* We enumarate the n-tuples x_1,...,x_n by ordering them in lexicographical order, if we are not doing
	   BIRACK_HOMOLOGY, we omit any that include a pair x_i = x_{i+1}
	*/
	if (braid_control::BIRACK_HOMOLOGY)
	{
		num_km1_generators = n*n;
		num_k_generators = n* num_km1_generators;
		num_kp1_generators = n* num_k_generators;
	}
	else
	{		
		num_km1_generators = create_chain_map(km1_chain_map,n,k-1);
		num_k_generators = create_chain_map(k_chain_map,n,k);		
		num_kp1_generators = create_chain_map(kp1_chain_map,n,k+1);
		
		for (int i=0; i< n_k_power; i++)
		{
			if (k_chain_map[i] != -1)
				reverse_k_chain_map[k_chain_map[i]] = i;
		}					
	}

	/* if k== 2 or 3 record the k_chain_map in the switch_data for invariant calculation */
	if (k==2 || k==3)
	{
		switch_data.num_chain_generators = num_k_generators;
		switch_data.chain_map = k_chain_map;
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)	
	debug << "birack_homology_generators:  k = " << k << ", num_km1_generators = " << num_km1_generators << " num_k_generators = " << num_k_generators << " num_kp1_generators = " << num_kp1_generators << endl;

if (debug_control::DEBUG >= debug_control::DETAIL)	
{
	if (!braid_control::BIRACK_HOMOLOGY)
	{
		debug << "birack_homology_generators:  km1_chain_map: " << endl;
		for (int i=0; i< n_km1_power; i++)
		{
			int place = i;
			debug << "birack_homology_generators: " << i << ": ";
			int factor = n_km1_power/n;
			for (int j=0; j<k-2; j++)
			{
				debug << place/factor << ' ';
				place %= factor;
				factor /=n;
			}
			debug << place << " = " << km1_chain_map[i] << endl;		
		}
	
		debug << "birack_homology_generators:  k_chain_map: " << endl;
		for (int i=0; i< n_k_power; i++)
		{
			int place = i;
			debug << "birack_homology_generators: "  << i << ": ";
			int factor = n_k_power/n;
			for (int j=0; j<k-1; j++)
			{
				debug << place/factor << ' ';
				place %= factor;
				factor /=n;
			}
			debug << place << " = " << k_chain_map[i] << endl;		
		}
		
		debug << "birack_homology_generators:  reverse k_chain_map: ";
		for (int i=0; i< num_k_generators; i++)
			debug << reverse_k_chain_map[i] << ' ';
		debug << endl;
		
	
		debug << "birack_homology_generators:  kp1_chain_map: " << endl;
		for (int i=0; i< n_kp1_power; i++)
		{
			int place = i;
			debug << "birack_homology_generators: "  << i << ": ";
			int factor = n_kp1_power/n;
			for (int j=0; j<k; j++)
			{
				debug << place/factor << ' ';
				place %= factor;
				factor /=n;
			}
			debug << place << " = " << kp1_chain_map[i] << endl;		
		}
	}
}	
	
	matrix<scalar> Delta_k(num_km1_generators,num_k_generators);
	matrix<scalar> Delta_kp1(num_k_generators,num_kp1_generators);
	
	create_boundary_map(Delta_k, km1_chain_map, twitch_u, twitch_d,n,k);
	create_boundary_map(Delta_kp1, k_chain_map, twitch_u, twitch_d,n,k+1);

if (debug_control::DEBUG >= debug_control::BASIC)	
{
	debug << "birack_homology_generators: Delta_k:" << endl;
	print(Delta_k, debug,3,"birack_homology_generators: ");
	debug << "birack_homology_generators: Delta_kp1:" << endl;
	print(Delta_kp1, debug,3,"birack_homology_generators: ");
}

if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "birack_homology_generators: Delta_k rows written as the coboundary generators of B^k; i.e the co-boundary of the characteristic functions X_i^{k-1}:" << endl;
	for (int i = 0; i< num_km1_generators; i++)
	{ 
		vector<scalar> coboundary(num_k_generators);
		for (int j=0; j< num_k_generators; j++)
			coboundary[j] = Delta_k[i][j];
		
		debug << "birack_homology_generators:   " << i+1 << ": ";
		print_characteristic_function(debug, coboundary, reverse_k_chain_map, cohomology, n, k);			
		debug << endl;
	}
}

	homology_control h_control;
	h_control.field_coefficients = field_coefficients;
	h_control.silent_operation = braid_control::SILENT_OPERATION;

//	homology_generators_return<scalar> homology_data = homology_generators(Delta_k, Delta_kp1, cohomology, field_coefficients, braid_control::SILENT_OPERATION);
	homology_generators_return<scalar> homology_data = homology_generators(Delta_k, Delta_kp1, cohomology, h_control);
	
	int num_generators = homology_data.num_generators;
	int num_torsion_generators = homology_data.num_torsion_generators;
	vector<scalar>& torsion = homology_data.torsion;
	
	list<vector<scalar> >::iterator gptr = homology_data.generators.begin();
	while (gptr !=homology_data.generators.end())
	{
		ostringstream oss;
		print_characteristic_function(oss, *gptr, reverse_k_chain_map, cohomology, n, k);			
		cocycles.push_back(oss.str());
		gptr++;
	}
			
if (debug_control::DEBUG >= debug_control::SUMMARY)	
{
	debug << "birack_homology_generators: found " << num_generators << (cohomology? " cohomology": " homology") << " generators" << endl;
	debug << "birack_homology_generators: switch_data: " << endl;
	print_switch_data(debug, switch_data, "birack_homology_generators:  ");
}

	if (!braid_control::SILENT_OPERATION)
	{
		cout << "Found " << num_generators << (cohomology? " cohomology": " homology") << " generators" << endl;
		
		if (num_generators !=0)
		{
			if (cohomology)
			{
				cout << (braid_control::BIRACK_HOMOLOGY?"H_B^":"H_{BQ}^");
				cout << k;
				if (scalar::variant == scalar::scalar_variant::MOD_P)
					cout << "(X,Z_" << mod_p::get_p() << ')';
				else if (scalar::variant == scalar::scalar_variant::RATIONAL || scalar::variant == scalar::scalar_variant::BIGRATIONAL)
					cout << "(X,Q)";
				else
					cout << "(X,Z)";
				
				cout << " is isomorphic to ";
			}
			else
			{
				cout << (braid_control::BIRACK_HOMOLOGY?"H^B^":"H^{BQ}^");
				cout << k;
				if (scalar::variant == scalar::scalar_variant::MOD_P)
					cout << "(X,Z_" << mod_p::get_p() << ')';
				else if (scalar::variant == scalar::scalar_variant::RATIONAL || scalar::variant == scalar::scalar_variant::BIGRATIONAL)
					cout << "(X,Q)";
				else
					cout << "(X,Z)";
				cout << " is isomorphic to ";
			}
			
			for (int i=0; i<num_torsion_generators; i++)
			{
				cout << "Z_" << torsion[i];
				if (i< num_generators-1)
					cout << " + ";
			}
			
			for (int i=num_torsion_generators; i<num_generators; i++)
			{
				if (scalar::variant == scalar::scalar_variant::MOD_P)
					cout << "Z_" << mod_p::get_p();
				else if (scalar::variant == scalar::scalar_variant::RATIONAL || scalar::variant == scalar::scalar_variant::BIGRATIONAL)
					cout << "Q";
				else
					cout << "Z";
				if (i< num_generators-1)
					cout << " + ";
			}
			cout << endl;
		}
	}
	
	if (!braid_control::RAW_OUTPUT)
	{
		output << "\n";
		if (braid_control::OUTPUT_AS_INPUT)
			output << ';';
		output << "Found " << num_generators << (cohomology? " cohomology": " homology") << " generators" << endl;

		if (num_generators !=0)
		{
			if (braid_control::OUTPUT_AS_INPUT)
				output << ';';

			if (cohomology)
			{
				output << (braid_control::BIRACK_HOMOLOGY?"H_B^":"H_{BQ}^");
				output << k;
				if (scalar::variant == scalar::scalar_variant::MOD_P)
					output << "(X,Z_" << mod_p::get_p() << ')';
				else if (scalar::variant == scalar::scalar_variant::RATIONAL || scalar::variant == scalar::scalar_variant::BIGRATIONAL)
					output << "(X,Q)";
				else
					output << "(X,Z)";
				output << " is isomorphic to ";
			}
			else
			{
				output << (braid_control::BIRACK_HOMOLOGY?"H^B^":"H^{BQ}^");
				output << k;
				if (scalar::variant == scalar::scalar_variant::MOD_P)
					output << "(X,Z_" << mod_p::get_p() << ')';
				else if (scalar::variant == scalar::scalar_variant::RATIONAL || scalar::variant == scalar::scalar_variant::BIGRATIONAL)
					output << "(X,Q)";
				else
					output << "(X,Z)";
				output << " is isomorphic to ";
			}

			for (int i=0; i<num_torsion_generators; i++)
			{
				output << "Z_" << torsion[i];
				if (i< num_generators-1)
					output << " + ";
			}
			
			for (int i=num_torsion_generators; i<num_generators; i++)
			{
				if (scalar::variant == scalar::scalar_variant::MOD_P)
					output << "Z_" << mod_p::get_p();
				else if (scalar::variant == scalar::scalar_variant::RATIONAL || scalar::variant == scalar::scalar_variant::BIGRATIONAL)
					output << "Q";
				else
					output << "Z";
				if (i< num_generators-1)
					output << " + ";
			}
			output << endl;
		}
	}

	/* if braid_control::HOMOLOGY or braid_control::COHOMOLOGY, write the list of generators to the output file */
	if (braid_control::HOMOLOGY || braid_control::COHOMOLOGY)
	{
		if (cocycles.size() == 0)
		{
			if (!braid_control::SILENT_OPERATION)
				cout << "no " << (cohomology? "cohomology":"homology") << " generators found" << endl;
				    
//			if (!braid_control::RAW_OUTPUT)
		 	{
				output << "no " << (cohomology? "cohomology":"homology") << " generators found" << endl;
		 	}	
		}
		else
		{
			list<string>::iterator lptr = cocycles.begin();		
			if (!braid_control::SILENT_OPERATION)
			{
				cout << cocycles.size() << (cohomology? " cohomology":" homology") << " generators:" << endl;
				lptr = cocycles.begin();
				while (lptr != cocycles.end())
				{
					cout << "C=" << *lptr << endl;
					lptr++;
				};
			}
	
			if (!braid_control::RAW_OUTPUT)
			{
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				output << cocycles.size() << (cohomology? " cohomology":" homology") << " generators:";
		   		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			}
		   	
			lptr = cocycles.begin();
			while (lptr != cocycles.end())
			{
				output << "C=" << *lptr;
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				lptr++;
			};
			output << flush;
		}
	}		
}

/* calculate_boundary evaluates the boundary of the given tuple that represents a tuple.size() cube in a finite birack of size n.  Note that n is *not* the size
   of the tuple.  The funtion stores the boundary in the column domain_index of the matrix Delta.  It is recorded as the coefficients in the tuples
   of the codomain, that is tuples of length tuple.size()-1.
   
   The codomain_map is used when braid_control::BIRACK_HOMOLOGY is false, since then the co-domain tuples need checking to ignore any that are degenerate.
*/
void calculate_boundary(int n, vector<int>& tuple, matrix<scalar>& Delta, int domain_index, vector<int>& codomain_map, matrix<int>& twitch_u, matrix<int>& twitch_d)
{
	
if (debug_control::DEBUG >= debug_control::DETAIL)	
{
	bool loc_newline = matrix_control::SINGLE_LINE_OUTPUT;
	matrix_control::SINGLE_LINE_OUTPUT = true;
	
	debug << "calculate_boundary: twitch_u = " << endl;
	print (twitch_u,debug,3,"calculate_boundary:   ");
	debug << "calculate_boundary: twitch_d = " << endl;
	print (twitch_d,debug,3,"calculate_boundary:   ");
	
	matrix_control::SINGLE_LINE_OUTPUT = loc_newline;
	
	debug << "calculate_boundary: codomain_map: ";
	for (size_t i=0; i< codomain_map.size(); i++)
		debug << codomain_map[i] << ' ';
	debug << endl;
}
	
	int t = tuple.size()-1; // the size of each boundary element tuple
	
	int max_factor = 1; // for calculating codomain_map hashes
	for (int i=0; i< t-1; i++)
		max_factor *= n;

if (debug_control::DEBUG >= debug_control::DETAIL)	
{
	debug << "calculate_boundary: partial_" << t+1 << " of t_tuple = ";
	for (int i=0; i<=t; i++)
		debug << tuple[i] << ' ';
	debug << endl;
	debug << "calculate_boundary: domain_index = " << domain_index << endl;	
	debug << "calculate_boundary: boundary element size t = " << t << " max_factor = " << max_factor << endl;
}

	for (int i=0; i< t+1; i++)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)	
	debug << "calculate_boundary:   i= " << i << " tuple element = " << tuple[i] << endl;

		vector<int> boundary_element(t);
		
		/* first boundary element */	
		for (int j=0; j<i; j++)
			boundary_element[j] = tuple[j];
		for (int j=i+1; j<=t; j++)
			boundary_element[j-1] = tuple[j];

if (debug_control::DEBUG >= debug_control::DETAIL)	
{
	debug << "calculate_boundary:     first boundary element ";
	for (int j=0; j< t; j++)
		debug << boundary_element[j] << ' ';
	debug << endl;
}
			
		scalar sign = scalar((i%2?1:-1));
		int hash=0;
		int factor = max_factor;
		for (int i=0; i< t; i++)
		{
			hash += boundary_element[i]*factor;
			factor/=n;
		}
		
		if (braid_control::BIRACK_HOMOLOGY)
		{
			/* the hash is the boundary element coefficient we are interested in */
			
if (debug_control::DEBUG >= debug_control::DETAIL)	
	debug << "calculate_boundary:       boundary element = " << hash << ", sign = " << sign << endl;

			Delta[hash][domain_index] += sign;

		}
		else
		{
			/* check the hash corresponds to a non-degenerate boundary element */

if (debug_control::DEBUG >= debug_control::DETAIL)	
	debug << "calculate_boundary:       boundary element hash = " << hash << ", codomain_map " << codomain_map[hash] << ", sign = " << sign << endl;

			if (codomain_map[hash] != -1)
				Delta[codomain_map[hash]][domain_index] += sign;
		}
		
		/* second boundary element */	
if (debug_control::DEBUG >= debug_control::DETAIL)	
	debug << "calculate_boundary:     second boundary element terms:" << endl;
	
		for (int j=0; j<i; j++)
		{
			boundary_element[j] = twitch_u[tuple[i]][tuple[j]];
if (debug_control::DEBUG >= debug_control::DETAIL)	
	debug << "calculate_boundary:       j = " << j << ": " << tuple[j] << " up " << tuple[i] << " = " << boundary_element[j] << endl;
		}
		
		for (int j=i+1; j<=t; j++)
		{
			boundary_element[j-1] = twitch_d[tuple[i]][tuple[j]];
if (debug_control::DEBUG >= debug_control::DETAIL)	
	debug << "calculate_boundary:       j = " << j << ": " << tuple[j] << " down " << tuple[i] << " = " << boundary_element[j-1] << endl;
		}			

if (debug_control::DEBUG >= debug_control::DETAIL)	
{
	debug << "calculate_boundary:     second boundary element ";
	for (int j=0; j< t; j++)
		debug << boundary_element[j] << ' ';
	debug << endl;
}
		sign = (i%2?-1:1);
		hash=0;
		factor = max_factor;
		for (int i=0; i< t; i++)
		{
			hash += boundary_element[i]*factor;
			factor/=n;
		}
		
		if (braid_control::BIRACK_HOMOLOGY)
		{
			/* the hash is the boundary element coefficient we are interested in */
			
if (debug_control::DEBUG >= debug_control::DETAIL)	
	debug << "calculate_boundary:       boundary element = " << hash << ", sign = " << sign << endl;

			Delta[hash][domain_index] += sign;

		}
		else
		{
			/* check the hash corresponds to a non-degenerate boundary element */
		
if (debug_control::DEBUG >= debug_control::DETAIL)	
	debug << "calculate_boundary:       boundary element hash = " << hash << ", codomain_map " << codomain_map[hash] << ", sign = " << sign << endl;

			if (codomain_map[hash] != -1)
				Delta[codomain_map[hash]][domain_index] += sign;	
		}
	}

if (debug_control::DEBUG >= debug_control::DETAIL)	
{
	debug << "calculate_boundary: produced boundary ";
	for (size_t i=0; i< Delta.numrows(); i++)
		debug << Delta[i][domain_index] << ' ';
	debug << endl;
}	
	
}

/* test_cohomology_generators checks the following cocycle condition for each of the scalar generators stored in switch_data

	For 3-cocycles:
			   LHS (a,b,c)+(a^c,b^c,d_c)+(a,c,d)+(b_a,c_a,d_a)
			   RHS (b,c,d)+(a^b,c_b,d_b)+(a,b,d)+(a^d,b^d,c^d)				  
	For 2-cocycles:
			   LHS (b,c)+(a^b,c_b)+(a,b)
			   RHS (b_a,c_a)+(a,c)+(a^c,b^c)				  
	
   it is used as a development test aid rather than providing calculation functionality
   
   9/1/26 Added the check for the "pseudo-rack" condition: \theta(b_a,c_a,d_a) = \theta(b,c,d) to look for biracks that are 
   not racks giving doubled colouring invariants
*/
void test_2_cohomology_generators(generic_switch_data& switch_data);
void test_3_cohomology_generators(generic_switch_data& switch_data);

void test_cohomology_generators(generic_switch_data& switch_data)
{
	if (braid_control::DOUBLE_BIRACKS)
		test_3_cohomology_generators(switch_data);
	else
		test_2_cohomology_generators(switch_data);
}
	
void test_2_cohomology_generators(generic_switch_data& switch_data)
{	
	if (!braid_control::SILENT_OPERATION)
		cout << "testing 2-cohomology generators" << endl;
	
	int n = switch_data.size;
	int num_chain_generators = switch_data.num_chain_generators;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "test_2_cohomology_generators: switch_data" << endl;
	print_switch_data(debug,switch_data,"test_2_cohomology_generators:   ");
}
	/* evaluate the sideways twitch map from the switch map in the switch_data */		
	matrix<int> twitch_u(n,n);
	matrix<int> twitch_d(n,n);
	
	matrix<int> inv_D(n,n);
	
	for (int i=0; i< n; i++)
	for (int j=0; j< n; j++)
		inv_D[i][switch_data.Sd[i][j]] = j; // inverts the map D_i given by the row switch_data.Sd[i][]
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "test_2_cohomology_generators: inv_D:" << endl;
	print (inv_D,debug,3,"test_2_cohomology_generators: ");
}
	
	/* The sideways twitch map of a switch is T(a,b) = (t_a(b),t^b(a)), i.e S_{-}^{+} (a,b) = (b_{a^{-1}}, a^{b_{a^{-1}}}), 
	   where ^ and _ are switch actions.
	*/
	for (int a=0; a< n; a++)
	for (int b=0; b< n; b++)
	{
		int x = inv_D[a][b];  // D_a^{-1}(b)

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "test_2_cohomology_generators:   a = " << a << " b = " << b << " x = " << x << endl;
	
		twitch_d[a][b] = x;
		twitch_u[b][a] = switch_data.Su[x][a];
	}


if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "test_2_cohomology_generators: twitch_u:" << endl;
	print (twitch_u,debug,3,"test_2_cohomology_generators: ");
	debug << "test_2_cohomology_generators: twitch_d:" << endl;
	print (twitch_d,debug,3,"test_2_cohomology_generators: ");
	
}		
	
	list<string>::iterator strptr = switch_data.cocycle_string.begin();
	list<vector<scalar> >::iterator lptr = switch_data.cocycle_scalar.begin();
	vector<int>& chain_map = switch_data.chain_map;
	
	while (lptr != switch_data.cocycle_scalar.end())
	{
		vector<scalar>& generator = *lptr;
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
{
	debug << "test_2_cohomology_generators: cohomology generator: " << *strptr << endl;
	debug << "test_2_cohomology_generators: cohomology generator: ";
	for (int i=0; i< num_chain_generators; i++)
		debug << generator[i] << ' ';
	debug << endl;
}
						
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "test_2_cohomology_generators: checking cocycle conditions: " << endl;
					
		for (int a=0; a<n; a++)
		for (int b=0; b<n; b++)
		for (int c=0; c<n; c++)
		{
			int non_degenerate_term;
				
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "test_2_cohomology_generators:   a = " << a << " b = " << b << " c = " << c;
	
			/* cocycle condition 
			   LHS (b,c)+(a^b,c_b)+(a,b)
			   RHS (b_a,c_a)+(a,c)+(a^c,b^c)				  
			*/
			vector<int> condition1_lhs(num_chain_generators);
			vector<int> condition1_rhs(num_chain_generators);


					
			int term = b*n+c;
			
			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_lhs[term] += 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L1 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
					condition1_lhs[non_degenerate_term] += 1;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L1 = " << term << "==" << non_degenerate_term;
			}
				
			int a_up_b = twitch_u[b][a];
			int c_down_b = twitch_d[b][c];
			term = a_up_b*n+c_down_b;

			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_lhs[term] += 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L2 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
					condition1_lhs[non_degenerate_term] += 1;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L2 = " << term << "==" << non_degenerate_term;
			}

			term = a*n+b;
			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_lhs[term] += 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L3 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
					condition1_lhs[non_degenerate_term] += 1;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L3 = " << term << "==" << non_degenerate_term;
			}

			int b_down_a = twitch_d[a][b];
			int c_down_a = twitch_d[a][c];
			term = b_down_a*n+c_down_a;

			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_rhs[term] += 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R1 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
				{
					condition1_rhs[non_degenerate_term] += 1;
				}
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R1 = " << term << "==" << non_degenerate_term;
			}

			term = a*n+c;
			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_rhs[term] += 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R2 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
					condition1_rhs[non_degenerate_term] += 1;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R2 = " << term << "==" << non_degenerate_term;
			}

			int a_up_c = twitch_u[c][a];
			int b_up_c = twitch_u[c][b];
			term = a_up_c*n+b_up_c;

			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_rhs[term] += 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R3 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
					condition1_rhs[non_degenerate_term] += 1;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R3 = " << term << "==" << non_degenerate_term;
			}
					
			scalar lhs_value=0;
			for (int i=0; i< num_chain_generators; i++)
				lhs_value += generator[i]*condition1_lhs[i];

			scalar rhs_value=0;
			for (int i=0; i< num_chain_generators; i++)
				rhs_value += generator[i]*condition1_rhs[i];
		
			if (lhs_value != rhs_value)
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
{
	debug << endl;
	debug << "test_2_cohomology_generators: condition1_lhs      : ";
	for (int i=0; i< num_chain_generators; i++)
		debug << setw(2) << condition1_lhs[i] << ' ';
	debug << endl;
	debug << "test_2_cohomology_generators: condition1_rhs      : ";
	for (int i=0; i< num_chain_generators; i++)
		debug << setw(2) << condition1_rhs[i] << ' ';
	debug << endl;
	debug << "test_2_cohomology_generators:   generator fails condition a=" << a << ", b=" << b << ", c=" << c
	      << " lhs_value = " << lhs_value << ", rhs_value = " << rhs_value  << endl;
}
				cout << "Error! Cohomology generator fails cocycle condition." << endl;
				exit(0);
			}
			else
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " lhs = rhs" << endl;
			}							
		}
				
		lptr++;
		strptr++;
	}

	/* write the final list of generators to the output file */
	if (switch_data.cocycle_string.size() == 0)
	{
		if (!braid_control::SILENT_OPERATION)
			cout << "no cohomology generators satisfying cocycle condition:" << endl;
			    
		if (!braid_control::RAW_OUTPUT)
	 	{
	   		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << "no cohomology generators satisfying cocycle condition:" << endl;
	 	}	
	}
	else
	{
		
		if (!braid_control::SILENT_OPERATION)
		{
			cout << switch_data.cocycle_string.size() << " cohomology generators satisfying cocycle condition:" << endl;
			strptr = switch_data.cocycle_string.begin();
			while (strptr != switch_data.cocycle_string.end())
			{
				cout << "C=" << *strptr << endl;
				strptr++;
			};
		}

		if (!braid_control::RAW_OUTPUT)
		{
	   		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << switch_data.cocycle_string.size() << " cohomology generators satisfying cocycle condition:";
	   		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			strptr = switch_data.cocycle_string.begin();
			while (strptr != switch_data.cocycle_string.end())
			{
				output << "C=" << *strptr;
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				strptr++;
			};
			output << flush;
		}

	}	
}

void test_3_cohomology_generators(generic_switch_data& switch_data)
{	
	if (!braid_control::SILENT_OPERATION)
		cout << "testing 3-cohomology generators" << endl;
	
	int n = switch_data.size;
	int nn = n*n;
	int num_chain_generators = switch_data.num_chain_generators;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "test_3_cohomology_generators: switch_data" << endl;
	print_switch_data(debug,switch_data,"test_3_cohomology_generators:   ");
}
	/* evaluate the sideways twitch map from the switch map in the switch_data */		
	matrix<int> twitch_u(n,n);
	matrix<int> twitch_d(n,n);
	
	matrix<int> inv_D(n,n);
	
	for (int i=0; i< n; i++)
	for (int j=0; j< n; j++)
		inv_D[i][switch_data.Sd[i][j]] = j; // inverts the map D_i given by the row switch_data.Sd[i][]
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "test_3_cohomology_generators: inv_D:" << endl;
	print (inv_D,debug,3,"test_3_cohomology_generators: ");
}
	
	/* The sideways twitch map of a switch is T(a,b) = (t_a(b),t^b(a)), i.e S_{-}^{+} (a,b) = (b_{a^{-1}}, a^{b_{a^{-1}}}), 
	   where ^ and _ are switch actions.
	*/
	for (int a=0; a< n; a++)
	for (int b=0; b< n; b++)
	{
		int x = inv_D[a][b];  // D_a^{-1}(b)

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "test_3_cohomology_generators:   a = " << a << " b = " << b << " x = " << x << endl;
	
		twitch_d[a][b] = x;
		twitch_u[b][a] = switch_data.Su[x][a];
	}


if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "test_3_cohomology_generators: twitch_u:" << endl;
	print (twitch_u,debug,3,"test_3_cohomology_generators: ");
	debug << "test_3_cohomology_generators: twitch_d:" << endl;
	print (twitch_d,debug,3,"test_3_cohomology_generators: ");
	
}		
	
	list<string>::iterator strptr = switch_data.cocycle_string.begin();
	list<vector<scalar> >::iterator lptr = switch_data.cocycle_scalar.begin();
	vector<int>& chain_map = switch_data.chain_map;
	
	while (lptr != switch_data.cocycle_scalar.end())
	{
		vector<scalar>& generator = *lptr;
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
{
	debug << "test_3_cohomology_generators: cohomology generator: " << *strptr << endl;
	debug << "test_3_cohomology_generators: cohomology generator: ";
	for (int i=0; i< num_chain_generators; i++)
		debug << generator[i] << ' ';
	debug << endl;
}
						
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "test_3_cohomology_generators: checking cocycle conditions: " << endl;
		
		bool pseudo_rack = true;
				
		for (int a=0; a<n; a++)
		for (int b=0; b<n; b++)
		for (int c=0; c<n; c++)
		for (int d=0; d<n; d++)
		{
			int non_degenerate_term;
				
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << "test_3_cohomology_generators:   a = " << a << " b = " << b << " c = " << c << " d = " << d;
	
			/* cocycle condition */
			vector<int> condition1_lhs(num_chain_generators);
			vector<int> condition1_rhs(num_chain_generators);
			
			/* pseudo-rack condition */
			vector<int> condition2_lhs(num_chain_generators);
			vector<int> condition2_rhs(num_chain_generators);
			
			/* set conditions: the bottom, middle, top ingress labels at each R3 move
			   involved in a tetrahedral move, determined by the co-originating generators
			   of the corresponding R3 cube.  See Carter, state sums paper, Figure 2.			   
			*/
			int term = a*nn+b*n+c;
			
			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_lhs[term] += 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L1 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
					condition1_lhs[non_degenerate_term] += 1;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L1 = " << term << "==" << non_degenerate_term;
			}
				
			int a_up_c = twitch_u[c][a];
			int b_up_c = twitch_u[c][b];
			int d_down_c = twitch_d[c][d];
			term = a_up_c*nn+b_up_c*n+d_down_c;

			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_lhs[term] += 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L2 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
					condition1_lhs[non_degenerate_term] += 1;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L2 = " << term << "==" << non_degenerate_term;
			}

			term = a*nn+c*n+d;
			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_lhs[term] += 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L3 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
					condition1_lhs[non_degenerate_term] += 1;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L3 = " << term << "==" << non_degenerate_term;
			}

			int b_down_a = twitch_d[a][b];
			int c_down_a = twitch_d[a][c];
			int d_down_a = twitch_d[a][d];
			term = b_down_a*nn+c_down_a*n+d_down_a;

			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_lhs[term] += 1;
				condition2_lhs[term] += 1; // pseudo-rack condition

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L4 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
				{
					condition1_lhs[non_degenerate_term] += 1;
					condition2_lhs[non_degenerate_term] += 1;
				}
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term L4 = " << term << "==" << non_degenerate_term;
			}

			term = b*nn+c*n+d;
			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_rhs[term] += 1;
				condition2_rhs[term] += 1; //pseudo-rack ondition

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R1 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
				{
					condition1_rhs[non_degenerate_term] += 1;
					condition2_rhs[non_degenerate_term] += 1;
				}
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R1 = " << term << "==" << non_degenerate_term;
			}

			int a_up_b = twitch_u[b][a];
			int c_down_b = twitch_d[b][c];
			int d_down_b = twitch_d[b][d];
			term = a_up_b*nn+c_down_b*n+d_down_b;
			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_rhs[term] += 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R2 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
					condition1_rhs[non_degenerate_term] += 1;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R2 = " << term << "==" << non_degenerate_term;
			}
			
			term = a*nn+b*n+d;
			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_rhs[term] += 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R3 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
					condition1_rhs[non_degenerate_term] += 1;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R3 = " << term << "==" << non_degenerate_term;
			}
			
			int a_up_d = twitch_u[d][a];
			int b_up_d = twitch_u[d][b];
			int c_up_d = twitch_u[d][c];
			term = a_up_d*nn+b_up_d*n+c_up_d;
			if (braid_control::BIRACK_HOMOLOGY)
			{
				condition1_rhs[term] += 1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R4 = " << term;
			}
			else
			{
				non_degenerate_term = chain_map[term];
				if (non_degenerate_term != -1)
					condition1_rhs[non_degenerate_term] += 1;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " term R4 = " << term << "==" << non_degenerate_term;
			}
			
			scalar lhs_value=0;
			for (int i=0; i< num_chain_generators; i++)
				lhs_value += generator[i]*condition1_lhs[i];

			scalar rhs_value=0;
			for (int i=0; i< num_chain_generators; i++)
				rhs_value += generator[i]*condition1_rhs[i];
		
			if (lhs_value != rhs_value)
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
{
	debug << endl;
	debug << "test_3_cohomology_generators: condition1_lhs      : ";
	for (int i=0; i< num_chain_generators; i++)
		debug << setw(2) << condition1_lhs[i] << ' ';
	debug << endl;
	debug << "test_3_cohomology_generators: condition1_rhs      : ";
	for (int i=0; i< num_chain_generators; i++)
		debug << setw(2) << condition1_rhs[i] << ' ';
	debug << endl;
	debug << "test_3_cohomology_generators:   generator fails condition a=" << a << ", b=" << b << ", c=" << c << ", d=" << d 
	      << " lhs_value = " << lhs_value << ", rhs_value = " << rhs_value  << endl;
}
				cout << "Error! Cohomology generator fails cocycle condition." << endl;
				exit(0);
			}
			else
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)	
	debug << " lhs = rhs" << endl;
			}				
			
			scalar lhs2_value=0;
			for (int i=0; i< num_chain_generators; i++)
				lhs2_value += generator[i]*condition2_lhs[i];

			scalar rhs2_value=0;
			for (int i=0; i< num_chain_generators; i++)
				rhs2_value += generator[i]*condition2_rhs[i];
		
			if (lhs2_value != rhs2_value)
				pseudo_rack = false; 
	
		}
		
		if (!pseudo_rack)
		{
			lptr = switch_data.cocycle_scalar.erase(lptr);
			strptr = switch_data.cocycle_string.erase(strptr);
		}
		else
		{	
			lptr++;
			strptr++;
		}
	}

	/* write the final list of generators to the output file */
	if (switch_data.cocycle_string.size() == 0)
	{
		if (!braid_control::SILENT_OPERATION)
			cout << "no cohomology generators satisfying cocycle condition:" << endl;
			    
		if (!braid_control::RAW_OUTPUT)
	 	{
	   		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << "no cohomology generators satisfying cocycle condition:" << endl;
	 	}	
	}
	else
	{
		
		if (!braid_control::SILENT_OPERATION)
		{
			cout << switch_data.cocycle_string.size() << " cohomology generators satisfying cocycle condition:" << endl;
			strptr = switch_data.cocycle_string.begin();
			while (strptr != switch_data.cocycle_string.end())
			{
				cout << "C=" << *strptr << endl;
				strptr++;
			};
		}

		if (!braid_control::RAW_OUTPUT)
		{
	   		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			output << switch_data.cocycle_string.size() << " cohomology generators satisfying cocycle condition:";
	   		output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
			strptr = switch_data.cocycle_string.begin();
			while (strptr != switch_data.cocycle_string.end())
			{
				output << "C=" << *strptr;
				output << (braid_control::OUTPUT_AS_INPUT? "\n;" : "\n");
				strptr++;
			};
			output << flush;
		}

	}	
}

