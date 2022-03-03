/**********************************************************************

This header file defines the control structure and enumeration types
used by the braid programme.

**********************************************************************/

struct braid_control 
{
	static string 		version;

	static int			gcd_initializing;

	static bool		AFFINE_INDEX;
	static bool		ALEXANDER;
	static bool		ALWAYS_CALCULATE_CODIMENSION_2_DELTA_1;
	static bool		ALWAYS_CALCULATE_DELTA_1;
	static bool		ARROW_POLYNOMIAL;
	static bool		BIGELOW_KNOT_SEARCH;
	static bool		BURAU;
	static bool		BYPASS_FUNDAMENTAL_EQUATION_CHECK;
	static bool		CALCULATE_DELTA_0_ONLY;
	static bool		CALCULATE_MOD_P;
	static bool		COMMUTATIVE_AUTOMORPHISM;
	static bool		COMPLEX_STUDY_DELTA_1;
	static bool    	CUSTOM_WEYL;
	static bool		DELTA_1_UNIT_CHECK;
	static bool     DEVELOPMENT_MODE;	
	static bool		DISPLAY_DELTA_1_ONLY;
	static bool		DOODLE_ALEXANDER;
	static bool		DOWKER_CODE;
	static bool		DYNNIKOV_TEST;
	static bool		EQUALITY_TEST;
	static bool		EVEN_WRITHE;
	static bool		EXPANDED_BRACKET_POLYNOMIAL;
	static bool     EXTRA_OUTPUT;
	static bool		FIXED_POINT_INVARIANT;
	static bool		FLAT_VOGEL_MOVES;
	static bool		FLIP_BRAID;
	static bool		GAUSS_CODE;
	static bool		GCD_BACK_SUBSTITUTING;
	static bool		HOMFLY;
	static bool		IMMERSION_CODE;
	static bool		INVERT_BRAID;
	static bool		JONES_POLYNOMIAL;
	static bool     KAMADA_DOUBLE_COVERING;
	static bool		KAUFFMAN_BRACKET;
	static bool		KNOTOID_BRACKET;
	static bool		LONG_KNOT;  // used by the switch polynomial invariant template functions that don't have access to generic code data
	static bool		LPGD;
	static bool		MATRIX;
	static bool     NORMALIZING_Q_POLYNOMIALS;
	static bool		NORMALIZE_BRACKET;
	static bool		NUMERATOR_GCD;
	static bool     OUTPUT_AS_INPUT;
	static bool		PARITY_ARROW;
	static bool		PARITY_BRACKET;
	static bool    	PEER_CODE;
	static bool    	PRIME_WEYL;
	static bool    	QUANTUM_WEYL;
	static bool     QUATERNION;
	static bool     RACK_POLYNOMIAL;
	static bool     RAW_OUTPUT;
	static bool		LINE_REFLECT_BRAID;
	static bool		PLANE_REFLECT_BRAID;
	static bool		RELAXED_PARITY;
	static bool		REMOVE_REIDEMEISTER_II_MOVES;
	static bool		REMOVE_PEER_CODE_COMPONENT;
//	static bool		REAL_STUDY_DELTA_1;
	static bool		SAWOLLEK;
	static bool		SILENT_OPERATION;
	static bool		STATUS_INFORMATION;
	static bool		STUDY_RHO_MAPPING;
	static bool		SWITCH_POLYNOMIAL_INVARIANT;
	static bool		TEST_MODE;
	static bool		TeX_POLYNOMIAL_OUTPUT;
	static bool		TMP_DIRECTORY;
	static bool    	TRUNCATED_WEYL;
	static bool		T_VARIABLE;
	static bool		ULPGD;
	static bool		VERIFY_DELTA_0;
	static bool		VOGEL_ALGORITHM;
	static bool		VOGEL_HEIGHT_ONLY;
	static bool    	WAIT_SWITCH;
	static bool		WELDED_BRAID;
	static bool    	WEYL;
	static bool 	ZIG_ZAG_DELTA;


	static bool		Sn_matrix_top_down; //true if Sn = Sn = I^{n-1} x S x I^{k-n-1}, false if Sn = I^{k-n-1} x S x I^{n-1}
	static bool     first_time;

	static int 		REMOVE_COMPONENT; // used to identify a component of a peer code to be removed
	static int		SWITCH_POWER; // used to control whether powers of switches are calculated.

	/* RACK_TERMS is the default for number of positive and negative terms added to
	a braid to calculate the rack polynomial.  The default is set in the function
	rack_poly_invariant
	*/
	static int		RACK_TERMS; 
	
	/* SATELLITE acts as a flag and indicates the number of strands to be added, if the default of 2 is not required */
	static int		SATELLITE;

	enum level {OFF, SUMMARY, BASIC, INTERMEDIATE, DETAIL, EXHAUSTIVE};
	static int DEBUG;
	static bool VOGEL_DEBUG;
	
	static int wait_threshold;
	static int wait_count;
	static int reset_count; // number of times wait_count has reached wait_threshold
	
};

/* braid_crossing_type provides an enumeration of the braid crossing types */
class braid_crossing_type 
{ 
	public:
	enum {NEGATIVE = -1, VIRTUAL = 0, POSITIVE = 1, FLAT = 2, NONE=9};
};
	
/* ST_pair_type is an enumeration of the type of pairs of finite switches S and T */
enum class ST_pair_type {FLAT_ESSENTIAL_VIRTUAL, ESSENTIAL_VIRTUAL, ESSENTIAL_WELDED, ESSENTIAL_DOODLE};
