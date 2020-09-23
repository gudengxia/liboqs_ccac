
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "packing.h"

/*************************************************
* pack the public key pk, 
* where pk = rho|t1
**************************************************/
void aigis_pack_pk(unsigned char pk[AIGIS_PK_SIZE_PACKED],
             const unsigned char rho[SEEDBYTES],
             const aigis_polyveck *t1)
{
  unsigned int i;

  for(i = 0; i < SEEDBYTES; ++i)
    pk[i] = rho[i];
  pk += SEEDBYTES;

  for(i = 0; i < PARAM_K; ++i)
    aigis_polyt1_pack(pk + i*POLT1_SIZE_PACKED, t1->vec+i);
}
void aigis_unpack_pk(unsigned char rho[SEEDBYTES],
               aigis_polyveck *t1,
               const unsigned char pk[AIGIS_PK_SIZE_PACKED])
{
  unsigned int i;

  for(i = 0; i < SEEDBYTES; ++i)
    rho[i] = pk[i];
  pk += SEEDBYTES;

  for(i = 0; i < PARAM_K; ++i)
    aigis_polyt1_unpack(t1->vec+i, pk + i*POLT1_SIZE_PACKED);
}

/*************************************************
* pack the secret key sk, 
* where sk = rho|key|hash(pk)|s1|s2|t0
**************************************************/
void aigis_pack_sk(unsigned char sk[AIGIS_SK_SIZE_PACKED],
             const unsigned char buf[2*SEEDBYTES + CRHBYTES],
             const aigis_polyvecl *s1,
             const aigis_polyveck *s2,
             const aigis_polyveck *t0)
{
  unsigned int i;

  for(i = 0; i < 2*SEEDBYTES + CRHBYTES; ++i)
    sk[i] = buf[i];
  sk += 2*SEEDBYTES + CRHBYTES;

  for(i = 0; i < PARAM_L; ++i)
    aigis_polyeta1_pack(sk + i*POLETA1_SIZE_PACKED, s1->vec+i);
  sk += PARAM_L*POLETA1_SIZE_PACKED;

  for(i = 0; i < PARAM_K; ++i)
    aigis_polyeta2_pack(sk + i*POLETA2_SIZE_PACKED, s2->vec+i);
  sk += PARAM_K*POLETA2_SIZE_PACKED;

  for(i = 0; i < PARAM_K; ++i)
    aigis_polyt0_pack(sk + i*POLT0_SIZE_PACKED, t0->vec+i);
}
void aigis_unpack_sk(unsigned char buf[2*SEEDBYTES + CRHBYTES],
               aigis_polyvecl *s1,
               aigis_polyveck *s2,
               aigis_polyveck *t0,
               const unsigned char sk[AIGIS_SK_SIZE_PACKED])
{
  unsigned int i;

  for(i = 0; i < 2*SEEDBYTES + CRHBYTES; ++i)
    buf[i] = sk[i];
  sk += 2*SEEDBYTES + CRHBYTES;

  for(i=0; i < PARAM_L; ++i)
    aigis_polyeta1_unpack(s1->vec+i, sk + i*POLETA1_SIZE_PACKED);
  sk += PARAM_L*POLETA1_SIZE_PACKED;
  
  for(i=0; i < PARAM_K; ++i)
    aigis_polyeta2_unpack(s2->vec+i, sk + i*POLETA2_SIZE_PACKED);
  sk += PARAM_K*POLETA2_SIZE_PACKED;

  for(i=0; i < PARAM_K; ++i)
    aigis_polyt0_unpack(t0->vec+i, sk + i*POLT0_SIZE_PACKED);
}

/*************************************************
* pack the signature sm, 
* where sm = z|h|c
**************************************************/
void aigis_pack_sig(unsigned char sm[AIGIS_SIG_SIZE_PACKED],
              const aigis_polyvecl *z,
              const aigis_polyveck *h,
              const aigis_poly *c)
{
  unsigned int i, j, k;
  uint64_t signs, mask;

  for(i = 0; i < PARAM_L; ++i)
    aigis_polyzpack(sm + i*POLZ_SIZE_PACKED, z->vec+i);
  sm += PARAM_L*POLZ_SIZE_PACKED;

  /* Encode h */
  k = 0;
  for(i = 0; i < PARAM_K; ++i) {
    for(j = 0; j < PARAM_N; ++j)
      if(h->vec[i].coeffs[j] == 1)
        sm[k++] = j;

    sm[OMEGA + i] = k;
  }
  while(k < OMEGA) sm[k++] = 0;
  sm += OMEGA + PARAM_K;
  
  /* Encode c */
  signs = 0;
  mask = 1;
  for(i = 0; i < PARAM_N/8; ++i) {
    sm[i] = 0;
    for(j = 0; j < 8; ++j) {
      if(c->coeffs[8*i+j] != 0) {
        sm[i] |= (1 << j);
        if(c->coeffs[8*i+j] == (PARAM_Q - 1)) signs |= mask;
        mask <<= 1;
      }
    }
  }
  sm += PARAM_N/8;
  for(i = 0; i < 8; ++i)
    sm[i] = signs >> 8*i;
}
void aigis_unpack_sig(aigis_polyvecl *z,
                aigis_polyveck *h,
                aigis_poly *c,
                const unsigned char sm[AIGIS_SIG_SIZE_PACKED])
{
  unsigned int i, j, k;
  uint64_t signs, mask;

  for(i = 0; i < PARAM_L; ++i)
    aigis_polyzunpack(z->vec+i, sm + i*POLZ_SIZE_PACKED);
  sm += PARAM_L*POLZ_SIZE_PACKED;

  /* Decode h */
  k = 0;
  for(i = 0; i < PARAM_K; ++i) {
    for(j = 0; j < PARAM_N; ++j)
      h->vec[i].coeffs[j] = 0;

    for(j = k; j < sm[OMEGA + i]; ++j)
      h->vec[i].coeffs[sm[j]] = 1;

    k = sm[OMEGA + i];
  }
  sm += OMEGA + PARAM_K;

  /* Decode c */
  for(i = 0; i < PARAM_N; ++i)
    c->coeffs[i] = 0;

  signs = 0;
  for(i = 0; i < 8; ++i)
    signs |= (uint64_t)sm[PARAM_N/8+i] << 8*i;

  mask = 1;
  for(i = 0; i < PARAM_N/8; ++i) {
    for(j = 0; j < 8; ++j) {
      if((sm[i] >> j) & 0x01) {
        c->coeffs[8*i+j] = (signs & mask) ? PARAM_Q - 1 : 1;
        mask <<= 1;
      }
    }
  }
}
