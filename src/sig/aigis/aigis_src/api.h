#ifndef AIGIS_API_H
#define AIGIS_API_H

#include "params.h"

#define AIGIS_SIG_SECRETKEYBYTES AIGIS_SK_SIZE_PACKED
#define AIGIS_SIG_PUBLICKEYBYTES AIGIS_PK_SIZE_PACKED
#define AIGIS_SIG_BYTES AIGIS_SIG_SIZE_PACKED

#define SIG_ALGNAME "Aigis-sig"

            
int msig_keygen(unsigned char *pk, unsigned char *sk);
                   
int msig_sign(unsigned char *sk, 
                   unsigned char *m, unsigned long long mlen, 
                   unsigned char *sm, unsigned long long *smlen);                     

int msig_verf(unsigned char *pk,
                   unsigned char *sm, unsigned long long smlen,
                   unsigned char *m, unsigned long long mlen);
#endif
