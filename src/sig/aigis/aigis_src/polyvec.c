#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "ntt.h"


void aigis_polyvecl_freeze2q(aigis_polyvecl *v) {
    unsigned int i;
    
    for(i = 0; i < PARAM_L; ++i)
        aigis_poly_freeze2q(&v->vec[i]);
}

void aigis_polyvecl_freeze4q(aigis_polyvecl *v) {
    unsigned int i;  
    for(i = 0; i < PARAM_L; ++i)
        aigis_poly_freeze4q(&v->vec[i]);
}

void aigis_polyvecl_add(aigis_polyvecl *w, const aigis_polyvecl *u, const aigis_polyvecl *v) {
  unsigned int i;

  for(i = 0; i < PARAM_L; ++i)
    aigis_poly_add(w->vec+i, u->vec+i, v->vec+i);
}

void aigis_polyvecl_ntt(aigis_polyvecl *v) {
  unsigned int i;

  for(i = 0; i < PARAM_L; ++i)
    aigis_poly_ntt(v->vec+i);
}

void aigis_polyvecl_pointwise_acc_invmontgomery(aigis_poly *w,
                                          const aigis_polyvecl *u,
                                          const aigis_polyvecl *v) 
{
  unsigned int i;
  aigis_poly t;

  aigis_poly_pointwise_invmontgomery(w, u->vec+0, v->vec+0);

  for(i = 1; i < PARAM_L; ++i) {
    aigis_poly_pointwise_invmontgomery(&t, u->vec+i, v->vec+i);
    aigis_poly_add(w, w, &t);
  }
  aigis_poly_barrat_reduce(w);
}

int aigis_polyvecl_chknorm(const aigis_polyvecl *v, uint32_t bound)  {
  unsigned int i;
  int ret = 0;

  for(i = 0; i < PARAM_L; ++i)
    ret |= aigis_poly_chknorm(v->vec+i, bound);

  return ret;
}

void aigis_polyveck_freeze2q(aigis_polyveck *v)  {
    unsigned int i;
    for(i = 0; i < PARAM_K; ++i)
        aigis_poly_freeze2q(&v->vec[i]);
}

void aigis_polyveck_freeze4q(aigis_polyveck *v)  {
    unsigned int i;
    for(i = 0; i < PARAM_K; ++i)
        aigis_poly_freeze4q(&v->vec[i]);
}

void aigis_polyveck_add(aigis_polyveck *w, const aigis_polyveck *u, const aigis_polyveck *v) {
  unsigned int i;

  for(i = 0; i < PARAM_K; ++i)
    aigis_poly_add(w->vec+i, u->vec+i, v->vec+i);
}

void aigis_polyveck_sub(aigis_polyveck *w, const aigis_polyveck *u, const aigis_polyveck *v) {
  unsigned int i;

  for(i = 0; i < PARAM_K; ++i)
    aigis_poly_sub(w->vec+i, u->vec+i, v->vec+i);
}

void aigis_polyveck_neg(aigis_polyveck *v) { 
  unsigned int i;

  for(i = 0; i < PARAM_K; ++i)
    aigis_poly_neg(v->vec+i);
}

void aigis_polyveck_shiftl(aigis_polyveck *v, unsigned int k) { 
  unsigned int i;

  for(i = 0; i < PARAM_K; ++i)
    aigis_poly_shiftl(v->vec+i, k);
}

void aigis_polyveck_ntt(aigis_polyveck *v) {
  unsigned int i;

  for(i = 0; i < PARAM_K; ++i)
    aigis_poly_ntt(v->vec+i);
}

void aigis_polyveck_invntt_montgomery(aigis_polyveck *v) {
  unsigned int i;

  for(i = 0; i < PARAM_K; ++i)
    aigis_poly_invntt_montgomery(v->vec+i);
}

int aigis_polyveck_chknorm(const aigis_polyveck *v, uint32_t bound) {
  unsigned int i;
  int ret = 0;

  for(i = 0; i < PARAM_K; ++i)
    ret |= aigis_poly_chknorm(v->vec+i, bound);

  return ret;
}

void aigis_polyveck_power2round(aigis_polyveck *v1, aigis_polyveck *v0, const aigis_polyveck *v) {
  unsigned int i;
  for(i = 0; i < PARAM_K; ++i)
  	aigis_poly_power2round(&v1->vec[i],&v0->vec[i], &v->vec[i]);
}

void aigis_polyveck_decompose(aigis_polyveck *v1, aigis_polyveck *v0, const aigis_polyveck *v) {
  unsigned int i;
  for(i = 0; i < PARAM_K; ++i)
    aigis_poly_decompose(&v1->vec[i],&v0->vec[i],&v->vec[i]);
}

unsigned int aigis_polyveck_make_hint(aigis_polyveck *h,
                                const aigis_polyveck *u,
                                const aigis_polyveck *v)
{
  unsigned int i,s = 0;
  for(i = 0; i < PARAM_K; ++i)
      s+=aigis_poly_make_hint(&h->vec[i],&u->vec[i],&v->vec[i]);
  
  return s;
}

void aigis_polyveck_use_hint(aigis_polyveck *w, const aigis_polyveck *u, const aigis_polyveck *h) {
  unsigned int i, j;
  aigis_poly v1,v0;
  for(i = 0; i < PARAM_K; ++i)
  {
   aigis_poly_decompose(&v1,&v0,&u->vec[i]);
    for(j = 0; j < PARAM_N; ++j)
    {
     if(h->vec[i].coeffs[j] == 0)
     w->vec[i].coeffs[j]=v1.coeffs[j];
  else if(v0.coeffs[j] > PARAM_Q)
    w->vec[i].coeffs[j] = (v1.coeffs[j] == (PARAM_Q - 1)/ALPHA - 1) ? 0 : v1.coeffs[j] + 1;
  else
    w->vec[i].coeffs[j] = (v1.coeffs[j] == 0) ? (PARAM_Q - 1)/ALPHA - 1 : v1.coeffs[j] - 1;
    }
  }
}
