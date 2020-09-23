
#include "rng.h"
//#include "params.h"
//#include "fips202.h"
//#include "randombytes.h"


int aigis_randombytes(unsigned char * r, unsigned long long r_byts)
{
	return aigis_rand_byts(r_byts, r);
}