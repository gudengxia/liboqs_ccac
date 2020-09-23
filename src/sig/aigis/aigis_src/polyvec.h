#ifndef POLYVEC_H
#define POLYVEC_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

/* Vectors of polynomials of length L */
typedef struct {
  aigis_poly vec[PARAM_L];
} aigis_polyvecl;

void aigis_polyvecl_freeze(aigis_polyvecl *v);
void aigis_polyvecl_freeze2q(aigis_polyvecl *v);
void aigis_polyvecl_freeze4q(aigis_polyvecl *v);

void aigis_polyvecl_add(aigis_polyvecl *w, const aigis_polyvecl *u, const aigis_polyvecl *v);

void aigis_polyvecl_ntt(aigis_polyvecl *v);
void aigis_polyvecl_pointwise_acc_invmontgomery(aigis_poly *w,
                                          const aigis_polyvecl *u,
                                          const aigis_polyvecl *v);

int aigis_polyvecl_chknorm(const aigis_polyvecl *v, uint32_t B);



/* Vectors of polynomials of length K */
typedef struct {
  aigis_poly vec[PARAM_K];
} aigis_polyveck;

void aigis_polyveck_freeze(aigis_polyveck *v);
void aigis_polyveck_freeze2q(aigis_polyveck *v);
void aigis_polyveck_freeze4q(aigis_polyveck *v);

void aigis_polyveck_add(aigis_polyveck *w, const aigis_polyveck *u, const aigis_polyveck *v);
void aigis_polyveck_sub(aigis_polyveck *w, const aigis_polyveck *u, const aigis_polyveck *v);
void aigis_polyveck_neg(aigis_polyveck *v);
void aigis_polyveck_shiftl(aigis_polyveck *v, unsigned int k);

void aigis_polyveck_ntt(aigis_polyveck *v);
void aigis_polyveck_invntt_montgomery(aigis_polyveck *v);

int aigis_polyveck_chknorm(const aigis_polyveck *v, uint32_t B);

void aigis_polyveck_power2round(aigis_polyveck *v1, aigis_polyveck *v0, const aigis_polyveck *v);
void aigis_polyveck_decompose(aigis_polyveck *v1, aigis_polyveck *v0, const aigis_polyveck *v);
unsigned int aigis_polyveck_make_hint(aigis_polyveck *h,
                                const aigis_polyveck *u,
                                const aigis_polyveck *v);
void aigis_polyveck_use_hint(aigis_polyveck *w, const aigis_polyveck *v, const aigis_polyveck *h);

#endif
