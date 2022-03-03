/**********************************************************************

This header file defines the #defines required by the programme
**********************************************************************/

/* component record values used by homfly polynomial calculation */
#define CR_CREATE				0
#define CR_SMOOTHED_CROSSING	1
#define CR_RESET_TO_CROSSING	2
#define CR_UNTWIST_CROSSING		3
#define CR_INVERSE_PAIR			4

#define CODIMENSION_1	1
#define CODIMENSION_2	2

/* bracket polynomial options */
#define KAUFFMAN_VARIANT	 1
#define JONES_VARIANT		 2
#define TURAEV_VARIANT		 3
#define ARROW_VARIANT		 4
#define PARITY_VARIANT		 5
#define PARITY_ARROW_VARIANT 6

/* smoothing options for bracket polynomial calculation */
#define SEIFERT_SMOOTHED 1
#define NON_SEIFERT_SMOOTHED -1
#define VIRTUAL_CROSSING 2
#define SHORTCUT_CROSSING	3
#define ODD_CROSSING	4
