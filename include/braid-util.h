//class generic_braid_data
namespace generic_braid_data
{

//public:

	/* braid_crossing_type provides an enumeration of the braid crossing types */
	enum crossing_type 
	{ 
		NEGATIVE = -1, 
		VIRTUAL = 0, 
		POSITIVE = 1, 
		FLAT = 2, 
		NONE=9,
	};

	enum character
	{
		unknown = 0x00000000,
		doodle  = 0x00000001,
		welded  = 0x00000002,
		flat    = 0x00000004,
	};

/*	   	
	unsigned int braid_character;
	int num_terms;
	int num_strands;
	//int rack_terms;  // used by rack_poly_invariant
	int turning_number;
	string braid_word;
	string title;

	generic_braid_data(): braid_character(character::unknown),num_terms(0),num_strands(0),turning_number(0),braid_word(""),title("") {}
	generic_braid_data(string w, string t); // parses qualifiers
*/
	
};


string invstr(string& str);
void flip_braid(string& braid, int num_terms, int num_strings);
void invert_braid(string& braid, int num_terms);
void plane_reflect_braid(string& braid, int num_terms);
void line_reflect_braid(string& braid, int num_terms, int num_strings);
void parse_braid(string word, int num_terms, vector<int>& braid_num, vector<int>& type);
void convert_rep (char*& word, bool silent);
int num_braid_strands(string braid_string);
bool valid_braid_input (string& input_string, int& num_terms, int& num_strings, bool raw_output=true, bool silent_output=true, bool output_as_input=false);
void shift_braid(string& braid, int num_terms);
