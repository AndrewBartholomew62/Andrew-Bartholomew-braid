struct input_control
{
	static unsigned int ACCEPT_MAP; //bitmap
	enum input_type
	{
		general         = 0x00000001, // unused
		braid_word      = 0x00000002, 
		immersion_code  = 0x00000004,
		peer_code       = 0x00000008,
		gauss_code      = 0x00000010,
		planar_diagram  = 0x00000020,
		lace_code       = 0x00000040,
		dowker_code     = 0x00000080,
/*		remainder = 0x00000020,
		equal =  	0x00000040,
		greater =  	0x00000080,
	    output =  	0x00000100,
		input =  	0x00000200,
		sum =  		0x00000400,
		diff =  	0x00000800,
		num_len =  	0x00001000,
		r_shift =  	0x00002000,
		l_shift =  	0x00004000,
		bool_conv =	0x00008000,
		gcd =  		0x00010000,
		carry =		0x00020000,
		// 'all' also supported
*/		
	};
};

void get_input_word (ifstream& input, string& buffer, string& title);
