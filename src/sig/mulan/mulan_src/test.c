#include <stdint.h>
#include <time.h>
#include "Alg.h"
#include "Parameters.h"
#include "Rounding.h"

#define NTESTS 10000
#define MLEN 59

int main(void)
{		
	uint32_t i;
	uint32_t ret;

	unsigned char m[MLEN];
	unsigned char sm[MLEN + MULAN_CRYPTO_BYTES];
	unsigned char pk[MULAN_CRYPTO_PUBLICKEYBYTES];
	unsigned char sk[MULAN_CRYPTO_SECRETKEYBYTES];
	
	unsigned long long pk_bytes, sk_bytes, sm_bytes;

	for(i=0; i<NTESTS; i++) 
	{
		mulan_randombytes(m, MLEN);

		mulan_sig_keygen(pk, &pk_bytes, sk, &sk_bytes);
		mulan_sig_sign(sk, sk_bytes, m, MLEN, sm, &sm_bytes);
		ret = mulan_sig_verf(pk, pk_bytes, sm, sm_bytes, m, MLEN);

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

