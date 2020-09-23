#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"
#include "fips202.h"

/*typedef __declspec(align(32)) struct {
  uint32_t coeffs[PARAM_N];
} aigis_poly;*/
typedef struct {
  uint32_t coeffs[PARAM_N];
} __attribute__((aligned(32))) aigis_poly; //fzhang

void aigis_poly_freeze2q(aigis_poly *a);
void aigis_poly_freeze4q(aigis_poly *a);
void aigis_poly_barrat_reduce(aigis_poly *a);
void aigis_poly_decompose(aigis_poly *r1,aigis_poly *r0, const aigis_poly *a);
void aigis_poly_power2round(aigis_poly *r1,aigis_poly *r0, const aigis_poly *a);
uint32_t aigis_poly_make_hint(aigis_poly *h, const aigis_poly *a, const aigis_poly *b);

void aigis_poly_add(aigis_poly *c, const aigis_poly *a, const aigis_poly *b);
void aigis_poly_sub(aigis_poly *c, const aigis_poly *a, const aigis_poly *b);
void aigis_poly_neg(aigis_poly *a);
void aigis_poly_shiftl(aigis_poly *a, unsigned int k);

void aigis_poly_ntt(aigis_poly *a);
void aigis_poly_invntt_montgomery(aigis_poly *a);
void aigis_poly_pointwise_invmontgomery(aigis_poly *c, const aigis_poly *a, const aigis_poly *b);

int  aigis_poly_chknorm(const aigis_poly *a, uint32_t B);
void aigis_poly_uniform(aigis_poly *a, unsigned char *buf);
void aigis_poly_uniform_eta1(aigis_poly *a,
                      const unsigned char seed[SEEDBYTES],
                      unsigned char nonce);
void aigis_poly_uniform_eta1_3x(aigis_poly *a0,
                          aigis_poly *a1,
                          aigis_poly *a2,
                          const unsigned char seed[SEEDBYTES],
                          unsigned char nonce0,
                          unsigned char nonce1,
                          unsigned char nonce2);
void aigis_poly_uniform_eta1_4x(aigis_poly *a0,
                         aigis_poly *a1,
                         aigis_poly *a2,
                         aigis_poly *a3,
                         const unsigned char seed[SEEDBYTES], 
                         unsigned char nonce0,
                         unsigned char nonce1,
                         unsigned char nonce2,
                         unsigned char nonce3);

void aigis_poly_uniform_eta2(aigis_poly *a,
                      const unsigned char seed[SEEDBYTES],
                      unsigned char nonce);
void aigis_poly_uniform_eta2_4x(aigis_poly *a0,
                         aigis_poly *a1,
                         aigis_poly *a2,
                         aigis_poly *a3,
                         const unsigned char seed[SEEDBYTES], 
                         unsigned char nonce0,
                         unsigned char nonce1,
                         unsigned char nonce2,
                         unsigned char nonce3);
void aigis_poly_uniform_eta2_2x(aigis_poly *a0,
	aigis_poly *a1,
	const unsigned char seed[SEEDBYTES],
	unsigned char nonce0,
	unsigned char nonce1);
void aigis_poly_uniform_gamma1m1(aigis_poly *a,
                           const unsigned char seed[SEEDBYTES + CRHBYTES],
                           uint16_t nonce);
void aigis_poly_uniform_gamma1m1_3x(aigis_poly *a0,
                              aigis_poly *a1,
                              aigis_poly *a2,
                              const unsigned char seed[SEEDBYTES + CRHBYTES],
                              uint16_t nonce0,
                              uint16_t nonce1,
                              uint16_t nonce2);
void aigis_poly_uniform_gamma1m1_4x(aigis_poly *a0,
                              aigis_poly *a1,
                              aigis_poly *a2,
                              aigis_poly *a3,
                              const unsigned char seed[SEEDBYTES + CRHBYTES],
                              uint16_t nonce0,
                              uint16_t nonce1,
                              uint16_t nonce2,
                              uint16_t nonce3);                        
void aigis_polyeta1_pack(unsigned char *r, const aigis_poly *a);
void aigis_polyeta1_unpack(aigis_poly *r, const unsigned char *a);
void aigis_polyeta2_pack(unsigned char *r, const aigis_poly *a);
void aigis_polyeta2_unpack(aigis_poly *r, const unsigned char *a);

void aigis_polyt1_pack(unsigned char *r, const aigis_poly *a);
void aigis_polyt1_unpack(aigis_poly *r, const unsigned char *a);

void aigis_polyt0_pack(unsigned char *r, const aigis_poly *a);
void aigis_polyt0_unpack(aigis_poly *r, const unsigned char *a);

void aigis_polyzpack(unsigned char *r, const aigis_poly *a);
void aigis_polyzunpack(aigis_poly *r, const unsigned char *a);

void aigis_polyw1_pack(unsigned char *r, const aigis_poly *a);
#endif
