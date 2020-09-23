#ifndef PACKING_H
#define PACKING_H

#include "polyvec.h"

void aigis_pack_pk(unsigned char pk[AIGIS_PK_SIZE_PACKED],
             const unsigned char rho[SEEDBYTES], const aigis_polyveck *t1);
void aigis_pack_sk(unsigned char sk[AIGIS_SK_SIZE_PACKED],
             const unsigned char buf[2*SEEDBYTES + CRHBYTES],
             const aigis_polyvecl *s1,
             const aigis_polyveck *s2,
             const aigis_polyveck *t0);
void aigis_pack_sig(unsigned char sm[AIGIS_SIG_SIZE_PACKED],
              const aigis_polyvecl *z, const aigis_polyveck *h, const aigis_poly *c);

void aigis_unpack_pk(unsigned char rho[SEEDBYTES], aigis_polyveck *t1,
               const unsigned char sk[AIGIS_SK_SIZE_PACKED]);
void aigis_unpack_sk(unsigned char buf[2*SEEDBYTES + CRHBYTES],
               aigis_polyvecl *s1,
               aigis_polyveck *s2,
               aigis_polyveck *t0,
               const unsigned char sk[AIGIS_SK_SIZE_PACKED]);
void aigis_unpack_sig(aigis_polyvecl *z, aigis_polyveck *h, aigis_poly *c,
                const unsigned char sm[AIGIS_SIG_SIZE_PACKED]);

#endif
