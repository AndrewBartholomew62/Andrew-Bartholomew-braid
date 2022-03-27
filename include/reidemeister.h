/*************************************************************************
                      reidemeister.h

***************************************************************************/

class simple_Reidemeister_II_return
{
public:
	bool Reidemeister_II_found;
	vector<int> RII_edge_flag;
	
	simple_Reidemeister_II_return(): Reidemeister_II_found(false) {}
};

class Reidemeister_III_return
{
public:
	bool Reidemeister_III_found;
	matrix<int> R3_configuration;
	
	Reidemeister_III_return(): Reidemeister_III_found(false) {}
};


class Gauss_Reidemeister_II_return
{
public:
	bool Reidemeister_II_found;
	int start_edge;
	bool forwards_from_peer;
	int floating_component;
	
	Gauss_Reidemeister_II_return(): Reidemeister_II_found(false),start_edge(0),forwards_from_peer(false),
	                                floating_component(0) {}
};

class virtual_Reidemeister_II_plus_return
{
public:
	bool Reidemeister_II_found;
	int start_edge;
	bool forwards_from_peer;
	bool move_peer_path_allowed;
	bool standard_Reidemeister_II_configuration;
	
	virtual_Reidemeister_II_plus_return(): Reidemeister_II_found(false),start_edge(0),forwards_from_peer(false),
	                                move_peer_path_allowed(false),standard_Reidemeister_II_configuration(false) {}
};


bool virtual_Reidemeister_I_present(generic_code_data code_data);
int remove_virtual_Reidemeister_I(generic_code_data& code_data, vector<int>& component_flags);
int remove_edge_flags_from_peer_code(generic_code_data& code_data, vector<int>& edge_flag, vector<int>& component_flags);
bool Reidemeister_II_labels (const int label1, const int label2);
bool simple_Reidemeister_II_labels (const int label1, const int label2);
bool Gauss_Reidemeister_II_labels(int start_edge, int end_bigon_edge_a, int end_bigon_edge_b, generic_code_data& code_data);
simple_Reidemeister_II_return simple_Reidemeister_II_present(generic_code_data& code_data, bool create_flags);
Gauss_Reidemeister_II_return Gauss_Reidemeister_II_present(generic_code_data& code_data, gauss_orientation_data& gauss_data);
virtual_Reidemeister_II_plus_return virtual_Reidemeister_II_plus(generic_code_data& code_data, matrix<int>& cycle, int num_cycles, int num_left_cycles, bool mark_special=false);
int remove_Reidemeister_II(generic_code_data& code_data, vector<int>& component_flags);
void align_label_numbers (matrix<int>& edge_labels, int num_components,vector<int> new_first_component_edge,vector<int> new_last_component_edge);
void clear_component_flag(int component, vector<int>& component_flags);
Reidemeister_III_return Reidemeister_III_present (generic_code_data code_data);
generic_code_data Reidemeister_III_move (generic_code_data code_data, int a, int b, int c);
