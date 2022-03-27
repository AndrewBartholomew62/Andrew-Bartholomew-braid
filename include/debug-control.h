/**********************************************************************

This header file defines the control structure and enumeration types
used by the debug subsystem.

**********************************************************************/

struct debug_control 
{
	enum level {OFF, SUMMARY, BASIC, INTERMEDIATE, DETAIL, EXHAUSTIVE};
	static int DEBUG;	
};
