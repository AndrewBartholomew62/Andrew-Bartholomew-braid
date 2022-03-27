/* braid_crossing_type provides an enumeration of the braid crossing types */
class braid_crossing_type 
{ 
	public:
	enum {NEGATIVE = -1, VIRTUAL = 0, POSITIVE = 1, FLAT = 2, NONE=9};
};

string invstr(string& str);
void flip_braid(string& braid);
void invert_braid(string& braid);
void plane_reflect_braid(string& braid);
void line_reflect_braid(string& braid);
void parse_braid(string word, int num_terms, vector<int>& braid_num, vector<int>& type);
void convert_rep (char*& word, bool silent);
bool valid_braid_input (string input_string, int& num_terms, int& num_strings, bool raw_output=true, bool silent_output=true, bool output_as_input=false);
