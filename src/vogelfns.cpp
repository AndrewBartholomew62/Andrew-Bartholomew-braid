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
#include <braid-util.h>
#include <debug-control.h>
#include <hash-defs.h>
#include <generic-code.h>

/********************* Function prototypes ***********************/
void parse_braid(string word, int num_terms, vector<int>& braid_num, vector<int>& type);
void set_component_record (int action, vector<int>& braid_num, vector<int>& type, vector<int>& component_record,
                           int base_crossing=0, int basepoint=0, int datum=0);
void vogel_error (string errstring);
int num_braid_terms(string word);
string braid_reduce (string word);
void braid_reduce (vector<int>& braid_num, vector<int>& type, int& basepoint, vector<int>& component_record);
void write_braid (ostream& s, vector<int>& braid_num, vector<int>& type);

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

if (debug_control::DEBUG >= debug_control::DETAIL) 
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
if (debug_control::DEBUG >= debug_control::DETAIL) 
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

if (debug_control::DEBUG >= debug_control::DETAIL) 
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

if (debug_control::DEBUG >= debug_control::DETAIL) 
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

if (debug_control::DEBUG >= debug_control::DETAIL) 
	debug << "braid_reduce: presented with " << word << endl;

	int num_terms = num_braid_terms(word);
		
if (debug_control::DEBUG >= debug_control::DETAIL) 
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

if (debug_control::DEBUG >= debug_control::DETAIL) 
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
		
if (debug_control::DEBUG >= debug_control::DETAIL) 
	debug << "braid_reduce: upper_crossing = " << upper_crossing << endl;

		int max_strand = 0;
		for (unsigned int i=0; i<component_record.size(); i++)
		{
			if (component_record[i] > max_strand)
				max_strand = component_record[i];
		}
		
if (debug_control::DEBUG >= debug_control::DETAIL) 
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
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "braid_reduce: num_initial_components = " << num_initial_components << endl;

	int passes_remaining = 2; // we want to go around the following loop at least twice because we'll cycle the braid
	do
	{
		if (num_terms == 0)
		{
			/* we've already reduced the braid to the empty word. */
			break;
		}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "braid_reduce: scanning crossings, passes_remaining = " << passes_remaining << endl;
	
		passes_remaining--; 

		bool reduced = false;
		/* look for pairs of adjacent inverse crossings and remove any we find.*/
		
		for (int i=0; i < num_terms-1; i++)
		{		
			/* compare the ith and i+1st crossing */

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "braid_reduce: i = " << i << ", comparing braid_num entries " << braid_num[i] << " and " << braid_num[i+1] 
	      << ", with types " << type[i] << " and " << type[i+1] << endl;

			if (	braid_num[i] == braid_num[i+1] && 
					( (type[i] == braid_crossing_type::POSITIVE && type[i+1] == braid_crossing_type::NEGATIVE) || 
					  (type[i] == braid_crossing_type::NEGATIVE && type[i+1] == braid_crossing_type::POSITIVE) || 
					  (type[i] == braid_crossing_type::VIRTUAL &&  type[i+1] == braid_crossing_type::VIRTUAL) 
					)
			   )
			{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "braid_reduce: base_crossing reduced by 2 to " << base_crossing << endl;
					}
					else if (base_crossing > i)
					{
						base_crossing--;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "braid_reduce: base_crossing decremented to " << base_crossing << endl;
					}
	
					reset_basepoint(base_strand, base_crossing, braid_num, type, num_terms, component_record);
				}
										
if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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


if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
	debug << "braid_reduce: final num_terms = " << num_terms << endl;

	/* write the base strand back to the basepoint parameter */
	basepoint = base_strand;

	/* rotate the braid_num and type so that the base_crossing appears at the front */
	rotate(braid_num.begin(), braid_num.begin()+base_crossing, braid_num.end());
	rotate(type.begin(), type.begin()+base_crossing, type.end());
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE) 
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

