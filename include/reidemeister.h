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
