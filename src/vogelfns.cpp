/************************************************************************
		   Functions for Vogel's algorithm

		  A. Bartholomew 20th October, 2002

**************************************************************************/
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <vector>
#include <list>
#include <iomanip>

using namespace std;

/********************* External variables ***********************/
//extern int* 		strand;
//extern int*			word_length;
//extern int			num_chains;
//extern int			tree_root;
//extern int*			crossing_type;
extern ofstream 	output;
extern ofstream     debug;


#include <util.h>
#include <matrix.h>
#include <braid-control.h>
#include <hash-defs.h>
#include <generic-code.h>

/********************* Function prototypes ***********************/
bool valid_braid_input (string input_string, int& num_terms, int& num_strings);
void parse_braid(string word, int num_terms, vector<int>& braid_num, vector<int>& type);
void set_component_record (int action, vector<int>& braid_num, vector<int>& type, vector<int>& component_record,
                           int base_crossing=0, int basepoint=0, int datum=0);
void vogel_error (string errstring);
int num_braid_terms(string word);
string braid_reduce (string word);
void braid_reduce (vector<int>& braid_num, vector<int>& type, int& basepoint, vector<int>& component_record);
void write_braid (ostream& s, vector<int>& braid_num, vector<int>& type);
void write_peer_code(ostream& s, const generic_code_data& code_data, bool zig_zags=false, bool labelled=true);
void print_code_data(generic_code_data& code_data, ostream& s, string prefix="");


void vogel_error (string errstring)
{
    cout << "\nError in Vogel algorithm: ";
    cout << errstring;
    cout.flush();
    output << "\nError in Vogel Algorithm: ";
    output << errstring;
    output.flush();
    output.close();

if (braid_control::VOGEL_DEBUG)
{
    debug << "vogel: error in Vogel Algorithm: ";
    debug << errstring;
    debug.flush();
    debug.close();
}

    exit(0);
}

int find_crossing (matrix<int>& seifert_circle, int circle, int preceeding)
{
    bool found;
    int place;

    found = false;
    for (int i=1; i<= seifert_circle[circle][0]; i++)
    {
		for (int j=1; j<=seifert_circle[preceeding][0]; j++)
		{
	    	if (seifert_circle[preceeding][j] ==
				   seifert_circle[circle][i] )
	    	{
				found = true;
				place = i;
				break;
	    	}
		}
		if (found) break;
    }
    return place;
}

/* on_circle returns the index of crossing on circle, in the range 1...n if it exists
   and zero otherwise
*/
int on_circle(matrix<int>& seifert_circle, int circle, int crossing)
{
    int i;
    for (i=1; i<=seifert_circle[circle][0]; i++)
    {
		if (seifert_circle[circle][i] == crossing)
	    	return i;
    }
    return 0;
}


/* reset_basepoint is used to identify the correct base_crossing for a basepoint on a given base_strand.  When the 
   braid_reduce function is called by the homfly function, although the braid_strand is correct, the base_crossing may 
   not be (due to smoothing crossings in the homfly skein relationship).  Similarly, when braid_reduce removes braid terms 
   the braid_strand may be  identified correctly but the base_crossing is not able to be determined by calculation around 
   the removed terms.   Therefore reset_basepoint checks the base_strand at the call by seeing if it is involved in a crossing. 
   If not, the base_strand has become an unknotted unlinked component of the braid, so we look for an alternative basepoint 
   by moving through the braid components according to the component_record.
   
   As we move along a strand looking to reset the base_crossing, we may pass crossings on other strands that mean the component 
   record is incorrect at the new base crossing.  We therefore have to reset the component record if the base_crossing is changed.
    
   At the call, we do not know on which component the basepoint lies, however, if we need to change components because the initial
   base strand has become unknotted and unlinked, then we will find the initial base_strand in the component record.  Note that
   we cannot assume in general that the basepoint lies on the lowest strand (e.g. if reset_basepoint is called after reducing
   -s1s2-s1-s3s2-s3s1 by removing the first and last terms.)

   We only need to check components after the one containing the basepoint since all previous components are good.  If we 
   reach the end of the component_record and have not found a new base_crossing we leave the existing value unchanged, although 
   the base strand (component) will have been moved to the last one in the component_record.
*/
void reset_basepoint (int& base_strand, int& base_crossing, vector<int>& braid_num, vector<int>& type,
                      int num_terms, vector<int>& component_record)
{
	int original_base_crossing = base_crossing;

if (braid_control::DEBUG >= braid_control::DETAIL) 
{
	debug << "reset_basepoint: braid_num: ";
	for (int i=0; i< num_terms; i++)
		debug << braid_num[i] << " ";
	debug << "\nreset_basepoint:      type: ";
	for (int i=0; i< num_terms; i++)
		debug << type[i] << " ";
	debug << endl;
	debug << "reset_basepoint: initial base_strand = " << base_strand << " original base_crossing = " << 
	                                                                                       original_base_crossing+1 << endl;
	debug << "reset_basepoint: component_record = ";
	for (unsigned int i=0; i< component_record.size(); i++)
		debug << component_record[i] << " ";
	debug << endl;						
}
	
	unsigned int cpt_index = 0; 
	bool cpt_index_not_set = true;
	bool found = false;
		
	while (true) // we break out of this loop
	{
		/* look for a crossing on strand component_record[i] */
if (braid_control::DEBUG >= braid_control::DETAIL) 
	debug << "reset_basepoint: looking for basepoint on strand " << base_strand << endl;


		/* Look along base_strand from base crossing for a crossing involving base_strand, 
		   looping back to the beginning of the strand if one is not found */
		for (int j=0; j < num_terms; j++)
		{
			int k = (j+base_crossing)%num_terms;
			
			if (braid_num[k] == base_strand || braid_num[k] == base_strand -1)
			{
				base_crossing = k;
				found = true;

if (braid_control::DEBUG >= braid_control::DETAIL) 
	debug << "reset_basepoint: found strand " << base_strand << " at crossing " << k+1 << endl;
				break;
			}
		}

		if (found)
		{
			if (base_crossing != original_base_crossing)
			{
				/* re-evaluate the component_record for the new base_crossing, we then have to rotate 
				   the braid so that the base_crossing is at the beginning again.  Note that set_component_record
				   counts crossings from 1 so we have to adjust the base_crossing accordingly.
				*/

				set_component_record(CR_RESET_TO_CROSSING, braid_num, type, component_record, 
				                                                          original_base_crossing, base_strand, base_crossing+1);
				rotate(braid_num.begin(), braid_num.begin()+base_crossing, braid_num.end());
				rotate(type.begin(), type.begin()+base_crossing, type.end());
				base_crossing = 0;

if (braid_control::DEBUG >= braid_control::DETAIL) 
{
debug << "reset_basepoint: base_crossing has moved, rotated braid to make base_crossing zero" << endl;
debug << "reset_basepoint: braid_num = ";
for ( unsigned int i=0; i< braid_num.size(); i++)
	debug << braid_num[i] << " ";
debug << "\nreset_basepoint:      type = ";
for ( unsigned int i=0; i< type.size(); i++)
	debug << type[i] << " ";
debug << "\nreset_basepoint: new component_record = ";
for (unsigned int i=0; i< component_record.size(); i++)
	debug << component_record[i] << " ";
debug << endl;						
}
			}

			break;
		}
		else
		{
			/* reset base_strand to the first strand of the next component.  First we have
			   to check to see whether we have identified the place of the original base_strand 
			   in the component_record.
			*/
			if (cpt_index_not_set)
			{
				/* we have to look through the component_record to find the base_strand
				   in order to initialize cpt_rindex correctly.
				*/
				for (unsigned int i=0; i< component_record.size(); i++)
				{
					if (component_record[i] == base_strand)
						cpt_index = i;
				}
				
				cpt_index_not_set = false;
			}
			
			/* now we can move forwards in the component_record */
			if (++cpt_index == component_record.size())
				break; // there are no more crossings after the current base_crossing in component_record order.
			else
				base_strand = component_record[cpt_index];
		}
	}
}


/* This is a wrapper function for the braid_reduce function that follows it, taking a string and
   determining the braid_num and type prior to calling the main braid_reduce function 
*/
string braid_reduce (string word)
{

if (braid_control::DEBUG >= braid_control::DETAIL) 
	debug << "braid_reduce: presented with " << word << endl;

	int num_terms = num_braid_terms(word);
		
if (braid_control::DEBUG >= braid_control::DETAIL) 
	debug << "braid_reduce: num_terms = " << num_terms << endl;

	/* determine the crossing numbers, types and the number of strings in the initial braid */
	
	vector<int> braid_num(num_terms);
    vector<int> type(num_terms);
	
	parse_braid(word, num_terms, braid_num, type);
	
	vector<int> component_record;
	set_component_record(CR_CREATE, braid_num, type, component_record);
	int num_cpts = component_record.size();
	int basepoint = 1;
	
    braid_reduce(braid_num, type, basepoint, component_record);
	
	/* adjust num_terms to reflect the number of terms after reduction */
	num_terms = braid_num.size();

if (braid_control::DEBUG >= braid_control::DETAIL) 
	debug << "braid_reduce: adjusted num_terms after reduction = " << num_terms << endl;
			
	/* now write out the resulting braid starting at the basepoint and return it. */
	ostringstream oss;
	
	if (num_terms == 0)
		oss << "unlink: " << num_cpts << " unknotted components";
	else
	{	
		/* check if there are any strands above the upper crossing left in the braid */
		int upper_crossing=0;
		for (int i=0; i<num_terms; i++)
		{
			if (braid_num[i] > upper_crossing)
				upper_crossing = braid_num[i];
		}
		
if (braid_control::DEBUG >= braid_control::DETAIL) 
	debug << "braid_reduce: upper_crossing = " << upper_crossing << endl;

		int max_strand = 0;
		for (unsigned int i=0; i<component_record.size(); i++)
		{
			if (component_record[i] > max_strand)
				max_strand = component_record[i];
		}
		
if (braid_control::DEBUG >= braid_control::DETAIL) 
	debug << "braid_reduce: max_strand = " << max_strand << endl;

//		if (max_strand > upper_crossing) changed 28/2/11
		if (max_strand > upper_crossing+1)
		{
			for (int i=0; i< num_terms; i++)
				braid_num[i] += max_strand-upper_crossing;
		}

		write_braid (oss, braid_num, type);
		
	}

	return oss.str();
}


/* braid_reduce tidies up a braid word, removing all positive and  negative pairs of real 
   crossings, then looking for a single occurrence of +/- s_i or t_i where i is either 1 
   or num_strings-1.  If the parameter basepoint is provided at the call it identifies the 
   strand (numbered from 1) at the start of the braid that represents the basepoint of the 
   braid's orientation.  In this case the location of the basepoint is tracked by the 
   function and the returned braid is written so that the basepoint lies at the start of 
   the braid again.  This means that if there is a crossing involving the base_strand in the 
   returned braid, the braid is cycled so that the first crossing after the basepoint is the 
   first term in the returned braid.  If there is no crossing involving the base_strand in the 
   returned braid, the base_strand has been reduced to an unlinked, unknotted component, so the 
   basepoint is moved to the next linked or knotted component according to the ordering of 
   components given at the call.  The braid is always returned with the basepoint at the front 
   of the braid, as above and the resulting basepoint strand is written to the basepoint parameter 
   for the return.
   
   The reduce function cycles the braid during its operation and so recalculates the component_record
   before returning, if we are tracking basepoints.  Note that the component record is only correct for
   the braid at the base_crossing, as the braid is cycled by braid_reduce.
*/
void braid_reduce (vector<int>& braid_num, vector<int>& type, int& basepoint, vector<int>& component_record)
{
	/* record the basepoint strand and crossing for the braid.  Braid strands are numbered from one, but if
	   we end up with a basepoint on an unknotted unlinked component we shall set the base_strand parameter to
	   zero to indicate this.  Crossings (braid terms, not crossings between strands) are numbered from zero. 
	*/
	int base_strand = 1;
	int base_crossing = 0;
	int upper_crossing = 0;
	int upper_count;
	int lower_crossing;
	int lower_count;
	int num_terms = braid_num.size();
	
	base_strand = basepoint;	
		
if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
{
	debug << "braid_reduce: initial base_crossing = " << base_crossing + 1 << ", base_strand = " << base_strand << endl;
	debug << "braid_reduce: component record = ";
	for ( unsigned int i=0; i< component_record.size(); i++)
		debug << component_record[i] << " ";
	debug << endl;
}
	/* call reset_basepoint to adjust the base_crossing to it's correct location, since base_crossing currently 
	   indicates the first crossing encountered when following the braid orientation from the given basepoint.  
	   Note that braid_strand could change as a result of this adjustment.		
	*/		
	reset_basepoint(base_strand, base_crossing, braid_num, type, num_terms, component_record);

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
{
	debug << "braid_reduce: adjusted base_crossing = " << base_crossing + 1 << ", base_strand = " << base_strand << endl;
	debug << "braid_reduce: braid_num: ";
	for (int i=0; i< num_terms; i++)
		debug << braid_num[i] << " ";
	debug << "\nbraid_reduce:      type: ";
	for (int i=0; i< num_terms; i++)
		debug << type[i] << " ";
	debug << endl;
}

	int num_initial_components = component_record.size();

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "braid_reduce: num_initial_components = " << num_initial_components << endl;

	int passes_remaining = 2; // we want to go around the following loop at least twice because we'll cycle the braid
	do
	{
		if (num_terms == 0)
		{
			/* we've already reduced the braid to the empty word. */
			break;
		}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "braid_reduce: scanning crossings, passes_remaining = " << passes_remaining << endl;
	
		passes_remaining--; 

		bool reduced = false;
		/* look for pairs of adjacent inverse crossings and remove any we find.*/
		
		for (int i=0; i < num_terms-1; i++)
		{		
			/* compare the ith and i+1st crossing */

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "braid_reduce: i = " << i << ", comparing braid_num entries " << braid_num[i] << " and " << braid_num[i+1] 
	      << ", with types " << type[i] << " and " << type[i+1] << endl;

			if (	braid_num[i] == braid_num[i+1] && 
					( (type[i] == braid_crossing_type::POSITIVE && type[i+1] == braid_crossing_type::NEGATIVE) || 
					  (type[i] == braid_crossing_type::NEGATIVE && type[i+1] == braid_crossing_type::POSITIVE) || 
					  (type[i] == braid_crossing_type::VIRTUAL &&  type[i+1] == braid_crossing_type::VIRTUAL) 
					)
			   )
			{

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "braid_reduce: inverse pair found at braid_num value " << braid_num[i] <<endl;

				/* we've got an inverse pair, first adjust the basepoint and copmponent_record if necessary.  
				   If the basepoint lies in front of the inverse pair, it remains on the same strand after 
				   they have been removed, and remains at the same crossing offset in braid_num and type. 
				   If it lies in between the pair to be removed, it changes strand through the second crossing 
				   and it's offset is decreased by one.
				   
				   If the base_crossing lies in between the inverse pair the component record is affected
				*/
				if (i+1 == base_crossing)
				{
					/* Adjust the component_record.  The two strands between the inverse pair are braid_num[i]
					   and braid_num[i]+1, if a component_record value indicates we start a component on either of
					   these two strands then it starts on the other one after the inverse pair is removed.
					*/
					set_component_record(CR_INVERSE_PAIR,braid_num, type, component_record,base_crossing,basepoint,braid_num[i]);

					if ( braid_num[i+1] == base_strand )
						base_strand++;
					else
						base_strand--;		
					
					base_crossing--;
				}

				/* now adjust braid_num and type */
				braid_num.erase(braid_num.begin()+i);
				braid_num.erase(braid_num.begin()+i);
				type.erase(type.begin()+i);
				type.erase(type.begin()+i);
				num_terms -= 2;
				reduced = true;		

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
{
	debug << "braid_reduce: reduced braid to: " << endl;
	debug << "braid_reduce: braid_num: ";
	for (unsigned int i=0; i< braid_num.size(); i++)
		debug << braid_num[i] << " ";
	debug << "\nbraid_reduce:      type: ";
	for (int i=0; i< num_terms; i++)
		debug << type[i] << " ";
	debug << endl;
}
				/* re-adjust base_crossing offset for the new number of terms and move forwards through the 
				   braid until we reach a crossing involving base_strand 
				*/
				if (num_terms)
				{
					if (base_crossing > i+1)
					{
						base_crossing -= 2;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "braid_reduce: base_crossing reduced by 2 to " << base_crossing << endl;
					}
					else if (base_crossing > i)
					{
						base_crossing--;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "braid_reduce: base_crossing decremented to " << base_crossing << endl;
					}
	
					reset_basepoint(base_strand, base_crossing, braid_num, type, num_terms, component_record);
				}
										
if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "braid_reduce: base_crossing = " << base_crossing + 1 << ", base_strand = " << base_strand << endl;

				/* step i backwards ready for the next loop */
				i--;
			}
		
		}

		/* identify the upper_crossing and lower_crossing and the number of their ocurrences left in the braid 
		   the upper_crossing has already been calculated, so we only need count its occurrences
		*/
		upper_crossing = 0;
		upper_count = 0;

		for (int i=0; i< num_terms; i++)
		{
			if (braid_num[i] > upper_crossing)
			{
				upper_crossing = braid_num[i];
				upper_count = 1;
			}
			else if (braid_num[i] == upper_crossing)
				upper_count++;
		}
		
if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "braid_reduce: upper_crossing = " << upper_crossing << ", upper_count = " << upper_count << endl;

		lower_crossing = upper_crossing;
		lower_count = 0;

		for (int i=0; i< num_terms; i++)
		{
			if (braid_num[i] < lower_crossing)
			{
				lower_crossing = braid_num[i];
				lower_count = 1;
			}
			else if (braid_num[i] == lower_crossing)
				lower_count++;
		}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "braid_reduce: lower_crossing = " << lower_crossing << ", lower_count = " << lower_count << endl;

		/* remove any single upper or lower crossing */
		if (upper_count == 1)
		{
			for (int i=0; i< num_terms; i++)
			{
				if (braid_num[i] == upper_crossing)
				{

					/* If the basepoint lies on the top strand, it will move to the strand 
					   below, otherwise it will remain the same.  
					*/
					if (i == base_crossing)
					{
						if ( braid_num[i] == base_strand-1)
							base_strand--;		

					}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
{
	debug << "braid_reduce: removing upper_crossing with number = " << braid_num[i] << ", type = " << type[i] << endl;
	debug << "braid_reduce: base_crossing = " << base_crossing + 1 << ", base_strand = " << base_strand << endl;
}

					/* Reduce any component_record values that lie above the upper strand that
					   is being untwisted to reflect the removal of the strand.
					*/
					for (unsigned int j=0; j < component_record.size(); j++)
					{
						if (component_record[j] > upper_crossing+1)
							component_record[j]--;
					}

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
{
	debug << "braid_reduce: component record adjusted to: ";
	for ( unsigned int i=0; i< component_record.size(); i++)
		debug << component_record[i] << " ";
	debug << endl;
}
					braid_num.erase(braid_num.begin()+i);
					type.erase(type.begin()+i);
					num_terms--;
								
					/* re-adjust base_crossing offset for the new number of terms and move forwards through the 
					   braid until we reach a crossing involving base_strand   The base_crossing offset into braid_num
					   and type will be unchanged, unless the erased crossing is at the end.
					*/
					if (num_terms)
					{
						if (base_crossing > i)
						{
							base_crossing--;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "braid_reduce: base_crossing decremented to " << base_crossing << endl;
						}
							
						reset_basepoint (base_strand, base_crossing, braid_num, type, num_terms, component_record);
					}
						
					reduced = true;		
					break;
				}
			}
		}

		if (lower_count == 1 && num_terms)  //avoids words like "s2" where upper_count has just reduced the braid to nothing
		{
			for (int i=0; i< num_terms; i++)
			{
				if (braid_num[i] == lower_crossing)
				{

					/* If the basepoint lies in front of the lower crossing on the lower strand, it will end up on its 
					   current strand number, having moved up through the crossing but then drop down a strand as the 
					   crossing is untwisted.  Otherwise (regardless of at which crossing the basepoint lies) the 
					   basepoint will end up one strand lower once the lower crossing has been untwisted.
					   
					   If it is not the case that the basepoint lies in front of the lower crossing and is on the 
					   lower of the strands involved in the lower crossing, the basepoint will be decremented.
					*/ 
					if (!(i == base_crossing && lower_crossing == base_strand))
							base_strand--;		

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
{
	debug << "braid_reduce: removing lower_crossing with number = " << braid_num[i] << ", type = " << type[i] << endl;
	debug << "braid_reduce: base_crossing = " << base_crossing + 1 << ", base_strand = " << base_strand << endl;
}

					/* Note that in this case the number of strands has been reduced, so the component_record is inaccurate.
					   The function set_component_record will adjust the component record to take account of the fact that a 
					   strand is being removed.
					*/
					if (num_terms != 0)
						set_component_record (CR_UNTWIST_CROSSING, braid_num, type, component_record, base_crossing, basepoint, i+1);

					braid_num.erase(braid_num.begin()+i);
					type.erase(type.begin()+i);
					num_terms--;

					/* Now reduce the braid numbers by one to reflect the removal of the strand.
					   (All of them, since they are all above the lower crossing .)
					*/
					for (int j=0; j< num_terms; j++)
							braid_num[j]--;

					/* re-adjust base_crossing offset for the new number of terms and reset the basepoint. The offset into 
					   braid_num  and type will be unchanged, unless the crossing is at the end.*/					
					if (num_terms)
					{
						if (base_crossing > i)
						{
							base_crossing--;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "braid_reduce: base_crossing decremented to " << base_crossing << endl;

						}
						reset_basepoint (base_strand, base_crossing, braid_num, type, num_terms, component_record);
					}

					reduced = true;		
					break;
				}
			}
		}

		if (reduced)
			passes_remaining++; // scan braid an additional time
		
		/* cycle the braid */
			
		if (num_terms != 0)
		{
			rotate(braid_num.begin(), braid_num.begin()+1, braid_num.end());
			rotate(type.begin(), type.begin()+1, type.end());
			
			/* adjust the basepoint, it has cycled backwards */
			if (num_terms != 0)
				base_crossing = (base_crossing + num_terms - 1) % num_terms;

if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
{
	debug << "braid_reduce: after cycling, braid becomes: " << endl;
	debug << "braid_reduce: braid_num: ";
	for (int i=0; i< num_terms; i++)
		debug << braid_num[i] << " ";
	debug << "\nbraid_reduce:      type: ";
	for (int i=0; i< num_terms; i++)
		debug << type[i] << " ";
	debug << "\nbraid_reduce: base_crossing = " << base_crossing + 1 << ", base_strand = " << base_strand << endl;
}
		}
				
	} while (passes_remaining);


if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
	debug << "braid_reduce: final num_terms = " << num_terms << endl;

	/* write the base strand back to the basepoint parameter */
	basepoint = base_strand;

	/* rotate the braid_num and type so that the base_crossing appears at the front */
	rotate(braid_num.begin(), braid_num.begin()+base_crossing, braid_num.end());
	rotate(type.begin(), type.begin()+base_crossing, type.end());
	
if (braid_control::DEBUG >= braid_control::INTERMEDIATE) 
{
	debug << "braid_reduce: returning braid: ";
	write_braid(debug,braid_num,type);
	debug << endl;
}
	
}

/* standard_rep takes a braid word and converts it into the standard representation 
   using a,b,c,A,B,C etc.

   If the braid word includes t symbols it does nothing
*/
void standard_rep (string& braid_word)
{
    char* local;
    char* mark;
    char* lptr;
    int number;
    bool not_complete = true;
    bool inverse;

	if (braid_word.find('t') != string::npos || braid_word.find('T') != string::npos)
		return;
	
	char* word = c_string(braid_word);
    char* wptr=word;
	
    local = new char[strlen(word)+1];
	lptr = local;
	
	do
	{
	    if (*wptr == '-')
		{
			inverse = true;
			wptr++;
		}
		else
		    inverse = false;

		/* skip over s character */
		wptr++;
		mark = wptr; /* mark where we start the number */

   		/* look for the end of the number */
   		while (isdigit(*wptr))
			wptr++;

	    get_number(number, mark);

		/* write this term into local */
		*lptr++ = (inverse?'A':'a')+number-1;
		
	    if (*wptr == '\0')
			not_complete = false;
	} while (not_complete);

	/* tidy up local, copy to word */
	*lptr = '\0';
	
	braid_word = local;
	
	delete[] word;
	delete[] local;
}

/* calculate_turning_cycles looks for left and right turning cycles in the given code_data and returns
   false if it is unable to calculate a set of cycles that may be realized.  From Eulers formula a
   realizable diagram cannot have more than num_crossings+2 turning cycles.  If the search therefore 
   exceeds this number of cycles the code_data cannot be realizable and the function returns false.
   If the search completes without reaching this limit the function returns true but this does not indicate
   that the code_data is necessarily realizable.  Note that the num_crossings+2 limit can only be breached
   when looking for right turning cycles, given that we look for left turning cycles first.
*/
bool calculate_turning_cycles(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles)
{

if (braid_control::VOGEL_DEBUG)
{
	debug << "calculate_turning_cycles: code data: ";
	write_peer_code (debug, code_data);
	debug << endl;
	print_code_data(code_data,debug,"calculate_turning_cycles: ");	
}

	matrix<int>& code_table = code_data.code_table;
	vector<int> term_crossing = code_data.term_crossing;
	vector<int> orig_crossing = code_data.orig_crossing;
	
	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;

	num_cycles = 0;
	
	/* First look for left turning cycles */
	for (int i=0; i<2*num_crossings; i++)
	{
		
if (braid_control::VOGEL_DEBUG)
    debug << "calculate_turning_cycles: edge = " << i;

    	/* does edge i already appear in a cycle ? */
		bool found = false;
		for (int j=0; j<num_cycles; j++)
		{
			for (int k=1; k<= cycle[j][0]; k++)
			{
//				if (abs(cycle[j][k]) == i)
				if (cycle[j][k] == i)
				{

if (braid_control::VOGEL_DEBUG)
    debug << " found in left turning cycle " << j << endl;

					found = true;
					break;
				}
			}
		}

		if (!found)
		{

if (braid_control::VOGEL_DEBUG)
    debug << " not found in current left turning cycles "<< endl;

			/* start a new cycle */
			int column = 1;

			/* we always traverse odd edges backwards */
//			int edge = (i%2? -i: i);
			int edge = i;
			cycle [num_cycles][column++] = edge;
			bool complete = false;

			/* a cycle cannot be longer than num_edges so we check that we do indeed
			   cycle within that limit.  This is used to check that the component map
			   within the code_table is valid, since an unrealizable component map can
			   result in infinite loops being calculated in the terminating and 
			   originating edges of the crossings where we never return to the start
			   of a cycle.
			*/
			for (int j=0; !complete && j< num_edges; j++)
			{
				/* determine which vertex the current edge takes us to and 
	   			the next edge we turn onto */
				if (edge % 2) // edge is odd (and negative)
				{
//					int vertex = orig_crossing[-edge];
					int vertex = orig_crossing[edge];
					if (code_table[TYPE][vertex] == generic_code_data::TYPE1)
		    			edge = code_table[EVEN_ORIGINATING][vertex];
					else
//						edge = -code_table[ODD_TERMINATING][vertex];
						edge = code_table[ODD_TERMINATING][vertex];

if (braid_control::VOGEL_DEBUG)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;

				}
				else // edge is even (and positive)
				{
					int vertex = term_crossing[edge];
					if (code_table[TYPE][vertex] == generic_code_data::TYPE1)
//						edge = -code_table[ODD_TERMINATING][vertex];
						edge = code_table[ODD_TERMINATING][vertex];
					else
		    			edge = code_table[EVEN_ORIGINATING][vertex];

if (braid_control::VOGEL_DEBUG)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;
				}

    
				if (edge == cycle[num_cycles][1])
				{
					complete = true;
					cycle[num_cycles][0] = column-1;
					num_cycles++;
				}
				else
				{
					cycle[num_cycles][column++] = edge;
				}				
			}
			
			if (!complete)
			{
				/* we've encounterd an infinte loop */
if (braid_control::VOGEL_DEBUG)
    debug << "calculate_turning_cycles: exceeded the maximum possible length of a turning cycle in a realizable code" << endl;
				return false;
			}
		}
	}	

	/* record the number of left cycles */
	num_left_cycles = num_cycles;
		
if (braid_control::VOGEL_DEBUG)
{
    debug << "calculate_turning_cycles: number of left turning cycles = " << num_left_cycles;
	for (int i=0; i<num_left_cycles; i++)
	{
		debug << "\ncalculate_turning_cycles: cycle " << i << " length = " << cycle[i][0] << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " ";
	}
	debug << endl;
}
		
	/* Now look for right turning cycles */

	for (int i=0; i<2*num_crossings; i++)
	{

if (braid_control::VOGEL_DEBUG)
    debug << "calculate_turning_cycles: edge = " << i;

		/* does edge i already appear in a right cycle ? */
		bool found = false;
		for (int j=num_left_cycles; j<num_cycles; j++)
		{
			for (int k=1; k<= cycle[j][0]; k++)
			{
//				if (abs(cycle[j][k]) == i)
				if (cycle[j][k] == i)
				{
					
if (braid_control::VOGEL_DEBUG)
    debug << " found in right turning cycle " << j << endl;
    
    					found = true;
					break;
				}
			}
		}

		if (!found)
		{

if (braid_control::VOGEL_DEBUG)
    debug << " not found in current right turning cycles "<< endl;

			/* check we've not exceeded the maximum number of turning cycles */
			if (num_cycles == num_crossings+2)
			{
if (braid_control::VOGEL_DEBUG)
    debug << "calculate_turning_cycles: exceeded the maximum possible number of turning cycles in a realizable code" << endl;
				return false;
			}
			
			/* start a new cycle */
			int column = 1;
			
			/* we always traverse odd edges backwards */
//			int edge = (i%2? -i: i);
			int edge = i;
			cycle [num_cycles][column++] = edge;
			bool complete = false;

			do
			{
				/* determine which vertex the current edge takes us to and 
	   			the next edge we turn onto */
				if (edge % 2) // edge is odd (and negative)
				{
//					int vertex = orig_crossing[-edge];
					int vertex = orig_crossing[edge];
					if (code_table[TYPE][vertex] == generic_code_data::TYPE1)
//						edge = -code_table[ODD_TERMINATING][vertex];
						edge = code_table[ODD_TERMINATING][vertex];
					else
		    			edge = code_table[EVEN_ORIGINATING][vertex];

if (braid_control::VOGEL_DEBUG)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;

				}
				else // edge is even (and positive)
				{
					int vertex = term_crossing[edge];
					if (code_table[TYPE][vertex] == generic_code_data::TYPE1)
		    			edge = code_table[EVEN_ORIGINATING][vertex];
					else
//						edge = -code_table[ODD_TERMINATING][vertex];
						edge = code_table[ODD_TERMINATING][vertex];

if (braid_control::VOGEL_DEBUG)
    debug << "calculate_turning_cycles:   takes us to crossing " << vertex << " next edge = " << edge << endl;

				}

				if (edge == cycle[num_cycles][1])
				{
					complete = true;
					cycle[num_cycles][0] = column-1;
					num_cycles++;
				}
				else
				{
					cycle[num_cycles][column++] = edge;
				}				
			} while(!complete);			
		}
	}

if (braid_control::VOGEL_DEBUG)
{
    debug << "calculate_turning_cycles: number of right turning cycles = " << num_cycles - num_left_cycles;
	for (int i=num_left_cycles; i<num_cycles; i++)
	{
		debug << "\ncalculate_turning_cycles: cycle " << i << " length = " << cycle[i][0] << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " "	;
	}
	debug << endl;
}
	return true;
}

/* The function vertex_span counts the number of vertices that can be reached from the 
   initial crossing in the given generic code data.  This can be used to determine whether the
   generic code data can be realized by a connected immersion.  Optionally, the exclude vector
   may specify a number of edges to exclude from the search, this may be use to determine whether 
   the complement of the excluded edges is connected.
*/
list<int> vertex_span (generic_code_data& code_data, int initial_crossing, vector<int>* exclude=0)
{
	matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;
	list<int> vertex_list;
	vertex_list.push_back(initial_crossing);
	list<int>::iterator lptr = vertex_list.begin();
	vector<int> crossing_flag(num_crossings,0);
	crossing_flag[initial_crossing] = 1;
	vector<int> edge_flag(2*num_crossings,0);
	vector<int>& orig_crossing = code_data.orig_crossing;
	vector<int>& term_crossing = code_data.term_crossing;
	
if (braid_control::VOGEL_DEBUG)
{
    debug << "vertex_span: checking code ";
	write_peer_code (debug, code_data);
	debug << "\nvertex_span: starting from crossing " << initial_crossing;
}

	if (exclude != 0 && exclude->size())
	{
if (braid_control::VOGEL_DEBUG)
{
    debug << " but excluding edges ";
    for (unsigned int i=0; i< exclude->size(); i++)
		debug << (*exclude)[i] << ' ';
    debug << endl;
}
		for (unsigned int i=0; i< exclude->size(); i++)
			edge_flag[(*exclude)[i]] = 1;
	}
	else
	{
if (braid_control::VOGEL_DEBUG)
    debug << " with no excluded edges " << endl;
	}
	
	while (lptr != vertex_list.end())
	{
		/* Look along any edges not already considered to see if we can extend the list of vertices.  
		   We look back along terminating edges and forwards along originating edges.  If
		   edge_flag[i] is 1 then we have either looked along this edge already or have excluded it 
		   from the search; crossing_flag[i] indicates whether we already have crossing i on the list.
	    */
if (braid_control::VOGEL_DEBUG)
    debug << "vertex_span:   from crossing " << *lptr << endl;

		int edge = code_table[EVEN_TERMINATING][*lptr];
		int crossing;
		
		if (edge_flag[edge] == 0)
		{
			crossing = orig_crossing[edge];

if (braid_control::VOGEL_DEBUG)
    debug << "vertex_span:     going back along " << edge << " takes us to " << crossing;

		    if(crossing_flag[crossing] == 0)
		    {
				vertex_list.push_back(crossing);
				crossing_flag[crossing] = 1;

if (braid_control::VOGEL_DEBUG)
    debug << ", adding " << crossing << " to vertex_list" << endl;
			}
			else
			{
if (braid_control::VOGEL_DEBUG)
    debug << ", which is already on the list " << endl;
			}
			edge_flag[edge] = 1;
		}
		else
		{
if (braid_control::VOGEL_DEBUG)
    debug << "vertex_span:     don't need to consider edge " << edge << endl;
		}
		
		edge = code_table[ODD_TERMINATING][*lptr];
		
		if (edge_flag[edge] == 0)
		{
			crossing = orig_crossing[edge];

if (braid_control::VOGEL_DEBUG)
    debug << "vertex_span:     going back along " << code_table[ODD_TERMINATING][*lptr] << " takes us to " << crossing;

		    if(crossing_flag[crossing] == 0)
		    {
				vertex_list.push_back(crossing);
				crossing_flag[crossing] = 1;

if (braid_control::VOGEL_DEBUG)
    debug << ", adding " << crossing << " to vertex_list" << endl;
			}
			else
			{
if (braid_control::VOGEL_DEBUG)
    debug << ", which is already on the list " << endl;
			}
			edge_flag[edge] = 1;
		}
		else
		{
if (braid_control::VOGEL_DEBUG)
    debug << "vertex_span:     don't need to consider edge " << edge << endl;
		}
		
		edge = code_table[EVEN_ORIGINATING][*lptr];

		if (edge_flag[edge] == 0)
		{
			crossing = term_crossing[edge];

if (braid_control::VOGEL_DEBUG)
    debug << "vertex_span:     going forwards along " << code_table[EVEN_ORIGINATING][*lptr] << " takes us to " << crossing;

		    if(crossing_flag[crossing] == 0)
		    {
				vertex_list.push_back(crossing);
				crossing_flag[crossing] = 1;

if (braid_control::VOGEL_DEBUG)
    debug << ", adding " << crossing << " to vertex_list" << endl;
			}
			else
			{
if (braid_control::VOGEL_DEBUG)
    debug << ", which is already on the list " << endl;
			}
			edge_flag[edge] = 1;
		}
		else
		{
if (braid_control::VOGEL_DEBUG)
    debug << "vertex_span:     don't need to consider edge " << edge << endl;
		}
		
		edge = code_table[ODD_ORIGINATING][*lptr];
		if (edge_flag[edge] == 0)
		{
			crossing = term_crossing[edge];

if (braid_control::VOGEL_DEBUG)
    debug << "vertex_span:     going forwards along " << code_table[ODD_ORIGINATING][*lptr] << " takes us to " << crossing;

		    if(crossing_flag[crossing] == 0)
		    {
				vertex_list.push_back(crossing);
				crossing_flag[crossing] = 1;

if (braid_control::VOGEL_DEBUG)
    debug << ", adding " << crossing << " to vertex_list" << endl;
			}
			else
			{
if (braid_control::VOGEL_DEBUG)
    debug << ", which is already on the list " << endl;
			}
			edge_flag[edge] = 1;
		}
		else
		{
if (braid_control::VOGEL_DEBUG)
    debug << "vertex_span:     don't need to consider edge " << edge << endl;
		}
		
		lptr++;
	}
	return vertex_list;
}

/* in the case of links it is possible that the code data we have been given is disconnected, as in 
   [-3 -5 -1, -9 -11 7]/# # # # # #, therefore we first check that the code is connected, then check
   the Euler characteristic.  In the process of this check we return a set of turning cycles to the call.
*/
bool realizable_code_data(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles)
{
	int num_crossings = code_data.num_crossings;
	int num_edges = 2* num_crossings;
	
	/* first check that the code is connected by calling vertex_span starting at crossing 0 
	   and not excluding any edges.  If we cannot reach every crossing, the code cannot be
	   realized by a connected immersion.
	*/
	int crossing_span = vertex_span(code_data,0).size();
	if (crossing_span != num_crossings)
	{		
if (braid_control::VOGEL_DEBUG)
{
    debug << "realizable_code_data: vertex span is only " << crossing_span
         << " crossings, code cannot be realized by a connected immersion" << endl;
}   
		return false;
	}
	else
	{
if (braid_control::VOGEL_DEBUG)
    debug << "realizable_code_data: vertex span includes all the crossings, code can be realized by a connected immersion" << endl;
	}
	
	/* now calculate the turning cycles, this includes additional checks on the realizable nature of the code */	
	if (!calculate_turning_cycles(code_data, cycle, num_left_cycles, num_cycles))
		return false;

	/* next check each edge appears once in a left handed
	   and once in a right handed cycle */
	   
if (braid_control::VOGEL_DEBUG)
    debug << "realizable_code_data: check edges appear exactly once in left and right turning cycles" << endl;

	bool realizable = true;
	for (int i=0; i<num_edges; i++)
	{
if (braid_control::VOGEL_DEBUG)
    debug << "realizable_code_data:   edge " << i;
    
		int count = 0;
		for (int j=0; j<num_left_cycles; j++)
		{
			for (int k=1; k<= cycle[j][0]; k++)
			{
				if (abs(cycle[j][k]) == i)
				{
					if (++count == 2)
					{
if (braid_control::VOGEL_DEBUG)
    debug << " appears twice in left turning cycle " << j << endl;
						realizable = false;
						break;
					}
				}
			}
			if (!realizable)
				break;
		}
		
		if (!count)
		{
if (braid_control::VOGEL_DEBUG)
    debug << " does not appear at all in left turning cycles";
			realizable = false;
		}

		if (!realizable)
		{
			break;
		}
		else
		{
if (braid_control::VOGEL_DEBUG)
    debug << " appears exactly once in left turning cycles";
		}

		/* we're only still here if count = 1, so now check the right
		   turning cycles */

		count = 0;
		for (int j=num_left_cycles; j<num_cycles; j++)
		{
			for (int k=1; k<= cycle[j][0]; k++)
			{
				if (abs(cycle[j][k]) == i)
				{
					if (++count == 2)
					{
if (braid_control::VOGEL_DEBUG)
    debug << " but twice in right turning cycle " << j << endl;
						realizable = false;
						break;
					}
				}
			}
			if (!realizable)
				break;
		}
		
		if (!count)
		{
if (braid_control::VOGEL_DEBUG)
    debug << " but not at all in right turning cycles";
			realizable = false;
		}

		if (!realizable)
		{
			break;
		}
		else
		{
if (braid_control::VOGEL_DEBUG)
    debug << " and in right turning cycles" << endl;
		}
	}

	/* Now check that the number of cycles is equal to num_crossings+2 */
	if (num_cycles != num_crossings+2)
	{
		realizable = false;

if (braid_control::VOGEL_DEBUG)
    debug << "realizable_code_data: number of cycles != num_crossings+2" << endl;
	}
	
	if (realizable)
	{
if (braid_control::VOGEL_DEBUG)
    debug << "realizable_code_data: code is able to be realized by a connected immersion" << endl;
	}
	else
	{
if (braid_control::VOGEL_DEBUG)
    debug << "realizable_code_data: code cannot be realized by a connected immersion" << endl;
	}
	
	return realizable;
}

