#ifndef _RANLUXGEN_INC_
#define _RANLUXGEN_INC_

// random number generator class - uses ranlux generator.
// this code refactored from implementation by Luescher

#define RANLUX_MAX ( (1<<24)-1 )
#define mask_lo 0x00ffffffUL  // 2^24 - 1
#define mask_hi ~0x00ffffffUL
#define two24 16777216        // 2^24

class RanluxGenerator
{
private:

	// internal generator structure
	//------------------------------
	typedef struct {
		unsigned int i;
		unsigned int j;
		unsigned int n;
		unsigned int skip;
		unsigned int carry;
		unsigned long int u[24];
	} ranlux_state_t;


	// internal generator state
	//--------------------------
	ranlux_state_t local_ranlux_state;


	// incrementation of the generator state
	//---------------------------------------
	inline unsigned long int increment_state()
	{
		unsigned int i = local_ranlux_state.i;
		unsigned int j = local_ranlux_state.j;
		long int delta = local_ranlux_state.u[j] - local_ranlux_state.u[i] 
				- local_ranlux_state.carry;

		if (delta & mask_hi)
		{
			local_ranlux_state.carry = 1;
			delta &= mask_lo;
		}
		else 
		{
			local_ranlux_state.carry = 0;
		}

		local_ranlux_state.u[i] = delta;

		if (i==0)
			i = 23;
		else
			i--;

		local_ranlux_state.i = i;

		if (j == 0)
			j = 23;
		else
			j--;

		local_ranlux_state.j = j;

		return delta;
	}


	// set generator state
	//---------------------
	void ranlux_set(unsigned long int s)
	{
		int i;
		long int seed;

		if (s==0)
			s = 314159265;      /* default seed is 314159265 */

		seed = s;

		/* This is the initialization algorithm of F. James, widely in use
			for RANLUX. */

		for (i=0;i<24;i++)
		{
			unsigned long int k = seed/53668;
			seed = 40014*(seed-k*53668)-k*12211;
			if (seed<0)
				seed += 2147483563;
			
			local_ranlux_state.u[i] = seed%two24;
		}

		local_ranlux_state.i = 23;
		local_ranlux_state.j = 9;
		local_ranlux_state.n = 0;
		local_ranlux_state.skip = 389-24; // 389 => best decorrelation

		if (local_ranlux_state.u[23]&mask_hi)
			local_ranlux_state.carry = 1;
		else 
			local_ranlux_state.carry = 0;
	}


	// generator initialization
	//--------------------------
	void ranlux_init( unsigned long int seed)
	{
		// seed the generator
		ranlux_set(seed);
	}


	// get random number
	//-------------------
	unsigned long int ranlux_get()
	{
		const unsigned int skip = local_ranlux_state.skip;
		unsigned long int r = increment_state();
  
		local_ranlux_state.n++;

		if (local_ranlux_state.n == 24)
		{
			unsigned int i;
			local_ranlux_state.n = 0;
			
			for (i = 0; i < skip; i++)
				increment_state();
		}

		return r;
	}

public:
	float frand( void )
	{
		return (float)ranlux_get() / (float)RANLUX_MAX;
	}

	int irand( int iMax )
	{
		return ((int)( (iMax+1)*frand() )%iMax);
	}

	void init( unsigned long int s = 0 )
	{
		ranlux_set( s );
	}



};

#endif
