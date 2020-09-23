#include <stdint.h>
#include <time.h>
#include "api.h"
#include "Alg.h"
//#include "Parameters.h"
//#include "Rounding.h"
#include "randombytes.h"
#define NTESTS 100
#define MLEN 59

int main(void)
{		
	uint32_t i;
	uint32_t ret;

	unsigned char m[MLEN];
	unsigned char sm[AIGIS_SIG_BYTES];
	unsigned char pk[AIGIS_SIG_PUBLICKEYBYTES];
	unsigned char sk[AIGIS_SIG_SECRETKEYBYTES];
	
	unsigned long long pk_bytes, sk_bytes, sm_bytes;

	for(i=0; i<NTESTS; i++) 
	{
		aigis_randombytes(m, MLEN);

		aigis_sig_keygen(pk, &pk_bytes, sk, &sk_bytes);
		aigis_sig_sign(sk, sk_bytes, m, MLEN, sm, &sm_bytes);
        
		ret = aigis_sig_verf(pk, pk_bytes, sm, sm_bytes, m, MLEN);
        
		if(ret!=1)
		{
			printf("ret = %d, Verification failed. \n", ret);
			getchar();
			return -1;
		}
 	}

	printf("Finished successfully. \n\n");
	return 0;
}

