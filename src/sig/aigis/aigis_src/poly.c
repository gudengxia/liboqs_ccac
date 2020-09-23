#include <stdint.h>
#include <immintrin.h>
#include "fips202.h"
#include "fips202x4.h"
#include "params.h"
#include "ntt.h"
#include "poly.h"
#include  "aes256ctr.h"



void aigis_poly_freeze2q(aigis_poly *a) 
{
    int i;
    __m256i *pa = (__m256i *) a->coeffs;
	__m256i q8x = _mm256_set1_epi32(PARAM_Q);
    __m256i temp;
	
    for(i = 0; i < PARAM_N/8; ++i)
    {
      pa[i] = _mm256_sub_epi32(pa[i],q8x);
      temp  = _mm256_srai_epi32(pa[i],31);
      temp  = _mm256_and_si256(temp,q8x);
      pa[i] = _mm256_add_epi32(pa[i],temp);      
    }
}
void aigis_poly_freeze4q(aigis_poly *a) 
{
    int i;
    __m256i *pa  = (__m256i *) a->coeffs;
    __m256i q8x  = _mm256_set1_epi32(PARAM_Q);
    __m256i dq8x = _mm256_set1_epi32(2*PARAM_Q);
    __m256i temp;
  
    for(i = 0; i < PARAM_N/8; ++i)
    {
      pa[i] = _mm256_sub_epi32(pa[i],dq8x);
      temp  = _mm256_srai_epi32(pa[i],31);
      temp  = _mm256_and_si256(temp,dq8x);
      pa[i] = _mm256_add_epi32(pa[i],temp);
      pa[i] = _mm256_sub_epi32(pa[i],q8x);
      temp  = _mm256_srai_epi32(pa[i],31);
      temp  = _mm256_and_si256(temp,q8x);
      pa[i] = _mm256_add_epi32(pa[i],temp);      
    }
}
void aigis_poly_barrat_reduce(aigis_poly *a)
{
	int i;
    __m256i *pa  = (__m256i *) a->coeffs;
#ifndef _WIN64
    __m256i q4x  = _mm256_set1_epi32(PARAM_Q);
#else
	__m256i q4x  = _mm256_set1_epi64x(PARAM_Q);
#endif
    __m256i temp;
    
    for(i = 0; i < PARAM_N/8; ++i)
    {
#if PARAM_Q == 2021377    
      temp  = _mm256_srai_epi32(pa[i],21);
#elif PARAM_Q == 3870721
      temp  = _mm256_srai_epi32(pa[i],22);
#endif
      temp  = _mm256_mullo_epi32(temp,q4x);
      pa[i] = _mm256_sub_epi32(pa[i],temp);    
    }
}
void aigis_poly_decompose(aigis_poly *r1,aigis_poly *r0, const aigis_poly *a)
{
  int i;
  __m256i * pa = (__m256i *) a->coeffs;
  __m256i * pr0 = (__m256i *) r0->coeffs;
  __m256i * pr1 = (__m256i *) r1->coeffs;
  
  __m256i q8x = _mm256_set1_epi32(PARAM_Q);
  __m256i alpha8x = _mm256_set1_epi32(ALPHA);
  __m256i barrat_const8x = _mm256_set1_epi32(3);
  __m256i mask = _mm256_set1_epi32(1);
  __m256i gamma2p18x = _mm256_set1_epi32(GAMMA2+1);
  __m256i gamma2m18x = _mm256_set1_epi32(GAMMA2-1);
  __m256i  u,t,t1,t2;
 
  for(i = 0; i < PARAM_N/8; ++i)
  {
    
    u = _mm256_mullo_epi32(pa[i],barrat_const8x);
#if ALPHA == 336896
    u = _mm256_srli_epi32(u,20);
#elif ALPHA == 645120
    u = _mm256_srli_epi32(u,21);
#endif
    u  = _mm256_add_epi32(u,mask);
    t  = _mm256_mullo_epi32(u,alpha8x);
    t  = _mm256_sub_epi32(pa[i],t);
    
    t1 = _mm256_srai_epi32(t,31);
    t2 = _mm256_and_si256(t1,mask);
    u  = _mm256_sub_epi32(u,t2);
    t1 = _mm256_and_si256(t1,alpha8x);
    t  = _mm256_add_epi32(t,t1);
    
    
    t  = _mm256_sub_epi32(t,gamma2p18x);
    t1 = _mm256_srai_epi32(t,31);
    t1 = _mm256_and_si256(t1,alpha8x);
    t  = _mm256_add_epi32(t,t1);
    t  = _mm256_sub_epi32(t,gamma2m18x);
    
    t1 = _mm256_srai_epi32(t,31);
    t1 = _mm256_and_si256(t1,mask);
    u  = _mm256_add_epi32(u,t1);
    
    t1 = _mm256_srai_epi32(u,1);
    t2 = _mm256_srai_epi32(u,2);
    t1 = _mm256_and_si256(t1,t2);
    t1 = _mm256_and_si256(t1,mask);
    
    t2 = _mm256_add_epi32(t,q8x);
    pr0[i] = _mm256_sub_epi32(t2,t1);
    t2 = _mm256_sub_epi32(mask,t1);
    pr1[i] = _mm256_mullo_epi32(u,t2); 
  }	
}

void aigis_poly_power2round(aigis_poly *r1,aigis_poly *r0, const aigis_poly *a)
{
  int i;
  __m256i * pa = (__m256i *) a->coeffs;
  __m256i * pr0 = (__m256i *) r0->coeffs;
  __m256i * pr1 = (__m256i *) r1->coeffs;
  __m256i q8x = _mm256_set1_epi32(PARAM_Q);
  __m256i t1,t2,t3;
  __m256i mask,powerdm1p1,powerdm1m1; 
  t1 = _mm256_set1_epi32((1<<PARAM_D));
  mask = _mm256_set1_epi32((1<<PARAM_D)-1);
  powerdm1p1 = _mm256_set1_epi32((1<<(PARAM_D-1))+1);
  powerdm1m1 = _mm256_set1_epi32((1<<(PARAM_D-1))-1);
  
  for(i = 0; i < PARAM_N/8; ++i)
  {
  	t2 = _mm256_and_si256(pa[i],mask);
  	t2 = _mm256_sub_epi32(t2,powerdm1p1);
  	t3 = _mm256_srai_epi32(t2,31);
  	t3 = _mm256_and_si256(t3,t1);
  	t2 = _mm256_add_epi32(t2,t3);
  	t2 = _mm256_sub_epi32(t2,powerdm1m1);
  	pr0[i] = _mm256_add_epi32(t2,q8x);
  	t2 = _mm256_sub_epi32(pa[i],t2);
  	pr1[i] = _mm256_srai_epi32(t2,PARAM_D);
  }
}
uint32_t aigis_poly_make_hint(aigis_poly *h, const aigis_poly*a, const aigis_poly*b)
{
	int i;
	uint32_t s = 0;
	for (i = 0; i < PARAM_N; ++i)
		if (a->coeffs[i] <= GAMMA2 || a->coeffs[i]> PARAM_Q - GAMMA2 || (a->coeffs[i] == PARAM_Q - GAMMA2 && b->coeffs[i] == 0))
			h->coeffs[i] = 0;
		else
		{
			h->coeffs[i] = 1;
			s++;
		}
	return s;
}
void aigis_poly_add(aigis_poly *c, const aigis_poly *a, const aigis_poly *b)  {
  unsigned int i;
  
  __m256i * pa = (__m256i *) a->coeffs;
  __m256i * pb = (__m256i *) b->coeffs;
  __m256i * pc = (__m256i *) c->coeffs;
  	
  for(i = 0; i < PARAM_N/8; ++i)
    pc[i] = _mm256_add_epi32(pa[i],pb[i]);
     
}

void aigis_poly_sub(aigis_poly *c, const aigis_poly *a, const aigis_poly *b) {
  unsigned int i;

  __m256i * pa = (__m256i *) a->coeffs;
  __m256i * pb = (__m256i *) b->coeffs;
  __m256i * pc = (__m256i *) c->coeffs;
  __m256i dq8x = _mm256_set1_epi32(2*PARAM_Q);
  __m256i temp;
  		
  for(i = 0; i < PARAM_N/8; ++i)
  {
    temp  = _mm256_add_epi32(pa[i],dq8x);
    pc[i] = _mm256_sub_epi32(temp,pb[i]);
  }
    
}

void aigis_poly_neg(aigis_poly *a) {
  unsigned int i;
    
  __m256i * pa = (__m256i *) a->coeffs;
  __m256i dq8x = _mm256_set1_epi32(2*PARAM_Q);

  for(i = 0; i < PARAM_N/8; ++i)
    pa[i] = _mm256_sub_epi32(dq8x,pa[i]);
}

void aigis_poly_shiftl(aigis_poly *a, unsigned int k) {
  unsigned int i;
  __m256i * pa = (__m256i *) a->coeffs;
  for(i = 0; i < PARAM_N/8; ++i)
    pa[i] = _mm256_slli_epi32(pa[i],k);
}

void aigis_poly_ntt(aigis_poly *a) {
  aigis_ntt(a->coeffs);
}

void aigis_poly_invntt_montgomery(aigis_poly *a) 
{
	aigis_invntt(a->coeffs); 
}
void aigis_poly_pointwise_invmontgomery(aigis_poly *c, const aigis_poly *a, const aigis_poly *b) 
{
  int i;
  __m256i * pa = (__m256i *) a->coeffs;
  __m256i * pb = (__m256i *) b->coeffs;
  __m256i * pc = (__m256i *) c->coeffs;
  __m256i q8x = _mm256_set1_epi32(PARAM_Q);
  __m256i qinv8x = _mm256_set1_epi32(QINV);
  
  __m256i  t1,t2,t3,t4;
  
  for(i = 0; i < PARAM_N/8; ++i)
  {
      //mul
      t1 = _mm256_srli_epi64(pa[i],32);
      t2 = _mm256_srli_epi64(pb[i],32);
      
      t2 = _mm256_mul_epu32(t1,t2);
      t1 = _mm256_mul_epu32(pa[i],pb[i]);
      
      //reduce 
      t3 = _mm256_mul_epu32(t1,qinv8x);
      t4 = _mm256_mul_epu32(t2,qinv8x);
      
      t3 = _mm256_mul_epu32(t3,q8x);
      t4 = _mm256_mul_epu32(t4,q8x);
      
      t1 = _mm256_add_epi64(t1,t3);
      t2 = _mm256_add_epi64(t2,t4);
      
      t1 = _mm256_srli_epi64(t1,32);
      
      pc[i] = _mm256_blend_epi32(t1,t2,0xAA);
  }
}

/*
int aigis_poly_chknorm(const aigis_poly *a, uint32_t B) {
  unsigned int i;
  int32_t t;

  for(i = 0; i < PARAM_N; ++i) {
    t = (PARAM_Q-1)/2 - a->coeffs[i];
    t ^= (t >> 31);
    t = (PARAM_Q-1)/2 - t;

    if((uint32_t)t >= B)
      return 1;
  }

  return 0;
}
*/


int  aigis_poly_chknorm(const aigis_poly *a, uint32_t B)
{
	unsigned int i;
	__m256i * pa = (__m256i *) a->coeffs;
	__m256i hq8x = _mm256_set1_epi32((PARAM_Q - 1) / 2);
	__m256i B8x = _mm256_set1_epi32(B);
	__m256i t0, t1, tmask, r = _mm256_set1_epi32(0);


	for (i = 0; i < PARAM_N / 8; ++i) {
		t0 = _mm256_sub_epi32(hq8x, pa[i]);
		tmask = _mm256_srai_epi32(t0, 31);
		t0 = _mm256_xor_si256(t0, tmask);
		t0 = _mm256_sub_epi32(hq8x, t0);
		t1 = _mm256_max_epi32(t0, B8x);
		t1 = _mm256_cmpeq_epi32(t1, t0);
		r = _mm256_or_si256(r, t1);
	}

	return _mm256_movemask_epi8(r);//the value is non-zero if some a[i] > B
}
/*void aigis_poly_uniform(aigis_poly *a, unsigned char *buf) {
  unsigned int ctr, pos;
  uint32_t t;

  ctr = pos = 0;
  while(ctr < PARAM_N) {
    t  = buf[pos++];
    t |= (uint32_t)buf[pos++] << 8;
    t |= (uint32_t)buf[pos++] << 16;
#if QBITS == 21
    t &= 0x1FFFFF;
#elif QBITS == 22
    t &= 0x3FFFFF;
#endif

    if(t < PARAM_Q)
      a->coeffs[ctr++] = t;
  }
}*/

static const uint64_t idx[256][4] = {
{0x800000008,0x800000008,0x800000008,0x800000008},
{0x800000000,0x800000008,0x800000008,0x800000008},
{0x800000001,0x800000008,0x800000008,0x800000008},
{0x100000000,0x800000008,0x800000008,0x800000008},
{0x800000002,0x800000008,0x800000008,0x800000008},
{0x200000000,0x800000008,0x800000008,0x800000008},
{0x200000001,0x800000008,0x800000008,0x800000008},
{0x100000000,0x800000002,0x800000008,0x800000008},
{0x800000003,0x800000008,0x800000008,0x800000008},
{0x300000000,0x800000008,0x800000008,0x800000008},
{0x300000001,0x800000008,0x800000008,0x800000008},
{0x100000000,0x800000003,0x800000008,0x800000008},
{0x300000002,0x800000008,0x800000008,0x800000008},
{0x200000000,0x800000003,0x800000008,0x800000008},
{0x200000001,0x800000003,0x800000008,0x800000008},
{0x100000000,0x300000002,0x800000008,0x800000008},
{0x800000004,0x800000008,0x800000008,0x800000008},
{0x400000000,0x800000008,0x800000008,0x800000008},
{0x400000001,0x800000008,0x800000008,0x800000008},
{0x100000000,0x800000004,0x800000008,0x800000008},
{0x400000002,0x800000008,0x800000008,0x800000008},
{0x200000000,0x800000004,0x800000008,0x800000008},
{0x200000001,0x800000004,0x800000008,0x800000008},
{0x100000000,0x400000002,0x800000008,0x800000008},
{0x400000003,0x800000008,0x800000008,0x800000008},
{0x300000000,0x800000004,0x800000008,0x800000008},
{0x300000001,0x800000004,0x800000008,0x800000008},
{0x100000000,0x400000003,0x800000008,0x800000008},
{0x300000002,0x800000004,0x800000008,0x800000008},
{0x200000000,0x400000003,0x800000008,0x800000008},
{0x200000001,0x400000003,0x800000008,0x800000008},
{0x100000000,0x300000002,0x800000004,0x800000008},
{0x800000005,0x800000008,0x800000008,0x800000008},
{0x500000000,0x800000008,0x800000008,0x800000008},
{0x500000001,0x800000008,0x800000008,0x800000008},
{0x100000000,0x800000005,0x800000008,0x800000008},
{0x500000002,0x800000008,0x800000008,0x800000008},
{0x200000000,0x800000005,0x800000008,0x800000008},
{0x200000001,0x800000005,0x800000008,0x800000008},
{0x100000000,0x500000002,0x800000008,0x800000008},
{0x500000003,0x800000008,0x800000008,0x800000008},
{0x300000000,0x800000005,0x800000008,0x800000008},
{0x300000001,0x800000005,0x800000008,0x800000008},
{0x100000000,0x500000003,0x800000008,0x800000008},
{0x300000002,0x800000005,0x800000008,0x800000008},
{0x200000000,0x500000003,0x800000008,0x800000008},
{0x200000001,0x500000003,0x800000008,0x800000008},
{0x100000000,0x300000002,0x800000005,0x800000008},
{0x500000004,0x800000008,0x800000008,0x800000008},
{0x400000000,0x800000005,0x800000008,0x800000008},
{0x400000001,0x800000005,0x800000008,0x800000008},
{0x100000000,0x500000004,0x800000008,0x800000008},
{0x400000002,0x800000005,0x800000008,0x800000008},
{0x200000000,0x500000004,0x800000008,0x800000008},
{0x200000001,0x500000004,0x800000008,0x800000008},
{0x100000000,0x400000002,0x800000005,0x800000008},
{0x400000003,0x800000005,0x800000008,0x800000008},
{0x300000000,0x500000004,0x800000008,0x800000008},
{0x300000001,0x500000004,0x800000008,0x800000008},
{0x100000000,0x400000003,0x800000005,0x800000008},
{0x300000002,0x500000004,0x800000008,0x800000008},
{0x200000000,0x400000003,0x800000005,0x800000008},
{0x200000001,0x400000003,0x800000005,0x800000008},
{0x100000000,0x300000002,0x500000004,0x800000008},
{0x800000006,0x800000008,0x800000008,0x800000008},
{0x600000000,0x800000008,0x800000008,0x800000008},
{0x600000001,0x800000008,0x800000008,0x800000008},
{0x100000000,0x800000006,0x800000008,0x800000008},
{0x600000002,0x800000008,0x800000008,0x800000008},
{0x200000000,0x800000006,0x800000008,0x800000008},
{0x200000001,0x800000006,0x800000008,0x800000008},
{0x100000000,0x600000002,0x800000008,0x800000008},
{0x600000003,0x800000008,0x800000008,0x800000008},
{0x300000000,0x800000006,0x800000008,0x800000008},
{0x300000001,0x800000006,0x800000008,0x800000008},
{0x100000000,0x600000003,0x800000008,0x800000008},
{0x300000002,0x800000006,0x800000008,0x800000008},
{0x200000000,0x600000003,0x800000008,0x800000008},
{0x200000001,0x600000003,0x800000008,0x800000008},
{0x100000000,0x300000002,0x800000006,0x800000008},
{0x600000004,0x800000008,0x800000008,0x800000008},
{0x400000000,0x800000006,0x800000008,0x800000008},
{0x400000001,0x800000006,0x800000008,0x800000008},
{0x100000000,0x600000004,0x800000008,0x800000008},
{0x400000002,0x800000006,0x800000008,0x800000008},
{0x200000000,0x600000004,0x800000008,0x800000008},
{0x200000001,0x600000004,0x800000008,0x800000008},
{0x100000000,0x400000002,0x800000006,0x800000008},
{0x400000003,0x800000006,0x800000008,0x800000008},
{0x300000000,0x600000004,0x800000008,0x800000008},
{0x300000001,0x600000004,0x800000008,0x800000008},
{0x100000000,0x400000003,0x800000006,0x800000008},
{0x300000002,0x600000004,0x800000008,0x800000008},
{0x200000000,0x400000003,0x800000006,0x800000008},
{0x200000001,0x400000003,0x800000006,0x800000008},
{0x100000000,0x300000002,0x600000004,0x800000008},
{0x600000005,0x800000008,0x800000008,0x800000008},
{0x500000000,0x800000006,0x800000008,0x800000008},
{0x500000001,0x800000006,0x800000008,0x800000008},
{0x100000000,0x600000005,0x800000008,0x800000008},
{0x500000002,0x800000006,0x800000008,0x800000008},
{0x200000000,0x600000005,0x800000008,0x800000008},
{0x200000001,0x600000005,0x800000008,0x800000008},
{0x100000000,0x500000002,0x800000006,0x800000008},
{0x500000003,0x800000006,0x800000008,0x800000008},
{0x300000000,0x600000005,0x800000008,0x800000008},
{0x300000001,0x600000005,0x800000008,0x800000008},
{0x100000000,0x500000003,0x800000006,0x800000008},
{0x300000002,0x600000005,0x800000008,0x800000008},
{0x200000000,0x500000003,0x800000006,0x800000008},
{0x200000001,0x500000003,0x800000006,0x800000008},
{0x100000000,0x300000002,0x600000005,0x800000008},
{0x500000004,0x800000006,0x800000008,0x800000008},
{0x400000000,0x600000005,0x800000008,0x800000008},
{0x400000001,0x600000005,0x800000008,0x800000008},
{0x100000000,0x500000004,0x800000006,0x800000008},
{0x400000002,0x600000005,0x800000008,0x800000008},
{0x200000000,0x500000004,0x800000006,0x800000008},
{0x200000001,0x500000004,0x800000006,0x800000008},
{0x100000000,0x400000002,0x600000005,0x800000008},
{0x400000003,0x600000005,0x800000008,0x800000008},
{0x300000000,0x500000004,0x800000006,0x800000008},
{0x300000001,0x500000004,0x800000006,0x800000008},
{0x100000000,0x400000003,0x600000005,0x800000008},
{0x300000002,0x500000004,0x800000006,0x800000008},
{0x200000000,0x400000003,0x600000005,0x800000008},
{0x200000001,0x400000003,0x600000005,0x800000008},
{0x100000000,0x300000002,0x500000004,0x800000006},
{0x800000007,0x800000008,0x800000008,0x800000008},
{0x700000000,0x800000008,0x800000008,0x800000008},
{0x700000001,0x800000008,0x800000008,0x800000008},
{0x100000000,0x800000007,0x800000008,0x800000008},
{0x700000002,0x800000008,0x800000008,0x800000008},
{0x200000000,0x800000007,0x800000008,0x800000008},
{0x200000001,0x800000007,0x800000008,0x800000008},
{0x100000000,0x700000002,0x800000008,0x800000008},
{0x700000003,0x800000008,0x800000008,0x800000008},
{0x300000000,0x800000007,0x800000008,0x800000008},
{0x300000001,0x800000007,0x800000008,0x800000008},
{0x100000000,0x700000003,0x800000008,0x800000008},
{0x300000002,0x800000007,0x800000008,0x800000008},
{0x200000000,0x700000003,0x800000008,0x800000008},
{0x200000001,0x700000003,0x800000008,0x800000008},
{0x100000000,0x300000002,0x800000007,0x800000008},
{0x700000004,0x800000008,0x800000008,0x800000008},
{0x400000000,0x800000007,0x800000008,0x800000008},
{0x400000001,0x800000007,0x800000008,0x800000008},
{0x100000000,0x700000004,0x800000008,0x800000008},
{0x400000002,0x800000007,0x800000008,0x800000008},
{0x200000000,0x700000004,0x800000008,0x800000008},
{0x200000001,0x700000004,0x800000008,0x800000008},
{0x100000000,0x400000002,0x800000007,0x800000008},
{0x400000003,0x800000007,0x800000008,0x800000008},
{0x300000000,0x700000004,0x800000008,0x800000008},
{0x300000001,0x700000004,0x800000008,0x800000008},
{0x100000000,0x400000003,0x800000007,0x800000008},
{0x300000002,0x700000004,0x800000008,0x800000008},
{0x200000000,0x400000003,0x800000007,0x800000008},
{0x200000001,0x400000003,0x800000007,0x800000008},
{0x100000000,0x300000002,0x700000004,0x800000008},
{0x700000005,0x800000008,0x800000008,0x800000008},
{0x500000000,0x800000007,0x800000008,0x800000008},
{0x500000001,0x800000007,0x800000008,0x800000008},
{0x100000000,0x700000005,0x800000008,0x800000008},
{0x500000002,0x800000007,0x800000008,0x800000008},
{0x200000000,0x700000005,0x800000008,0x800000008},
{0x200000001,0x700000005,0x800000008,0x800000008},
{0x100000000,0x500000002,0x800000007,0x800000008},
{0x500000003,0x800000007,0x800000008,0x800000008},
{0x300000000,0x700000005,0x800000008,0x800000008},
{0x300000001,0x700000005,0x800000008,0x800000008},
{0x100000000,0x500000003,0x800000007,0x800000008},
{0x300000002,0x700000005,0x800000008,0x800000008},
{0x200000000,0x500000003,0x800000007,0x800000008},
{0x200000001,0x500000003,0x800000007,0x800000008},
{0x100000000,0x300000002,0x700000005,0x800000008},
{0x500000004,0x800000007,0x800000008,0x800000008},
{0x400000000,0x700000005,0x800000008,0x800000008},
{0x400000001,0x700000005,0x800000008,0x800000008},
{0x100000000,0x500000004,0x800000007,0x800000008},
{0x400000002,0x700000005,0x800000008,0x800000008},
{0x200000000,0x500000004,0x800000007,0x800000008},
{0x200000001,0x500000004,0x800000007,0x800000008},
{0x100000000,0x400000002,0x700000005,0x800000008},
{0x400000003,0x700000005,0x800000008,0x800000008},
{0x300000000,0x500000004,0x800000007,0x800000008},
{0x300000001,0x500000004,0x800000007,0x800000008},
{0x100000000,0x400000003,0x700000005,0x800000008},
{0x300000002,0x500000004,0x800000007,0x800000008},
{0x200000000,0x400000003,0x700000005,0x800000008},
{0x200000001,0x400000003,0x700000005,0x800000008},
{0x100000000,0x300000002,0x500000004,0x800000007},
{0x700000006,0x800000008,0x800000008,0x800000008},
{0x600000000,0x800000007,0x800000008,0x800000008},
{0x600000001,0x800000007,0x800000008,0x800000008},
{0x100000000,0x700000006,0x800000008,0x800000008},
{0x600000002,0x800000007,0x800000008,0x800000008},
{0x200000000,0x700000006,0x800000008,0x800000008},
{0x200000001,0x700000006,0x800000008,0x800000008},
{0x100000000,0x600000002,0x800000007,0x800000008},
{0x600000003,0x800000007,0x800000008,0x800000008},
{0x300000000,0x700000006,0x800000008,0x800000008},
{0x300000001,0x700000006,0x800000008,0x800000008},
{0x100000000,0x600000003,0x800000007,0x800000008},
{0x300000002,0x700000006,0x800000008,0x800000008},
{0x200000000,0x600000003,0x800000007,0x800000008},
{0x200000001,0x600000003,0x800000007,0x800000008},
{0x100000000,0x300000002,0x700000006,0x800000008},
{0x600000004,0x800000007,0x800000008,0x800000008},
{0x400000000,0x700000006,0x800000008,0x800000008},
{0x400000001,0x700000006,0x800000008,0x800000008},
{0x100000000,0x600000004,0x800000007,0x800000008},
{0x400000002,0x700000006,0x800000008,0x800000008},
{0x200000000,0x600000004,0x800000007,0x800000008},
{0x200000001,0x600000004,0x800000007,0x800000008},
{0x100000000,0x400000002,0x700000006,0x800000008},
{0x400000003,0x700000006,0x800000008,0x800000008},
{0x300000000,0x600000004,0x800000007,0x800000008},
{0x300000001,0x600000004,0x800000007,0x800000008},
{0x100000000,0x400000003,0x700000006,0x800000008},
{0x300000002,0x600000004,0x800000007,0x800000008},
{0x200000000,0x400000003,0x700000006,0x800000008},
{0x200000001,0x400000003,0x700000006,0x800000008},
{0x100000000,0x300000002,0x600000004,0x800000007},
{0x600000005,0x800000007,0x800000008,0x800000008},
{0x500000000,0x700000006,0x800000008,0x800000008},
{0x500000001,0x700000006,0x800000008,0x800000008},
{0x100000000,0x600000005,0x800000007,0x800000008},
{0x500000002,0x700000006,0x800000008,0x800000008},
{0x200000000,0x600000005,0x800000007,0x800000008},
{0x200000001,0x600000005,0x800000007,0x800000008},
{0x100000000,0x500000002,0x700000006,0x800000008},
{0x500000003,0x700000006,0x800000008,0x800000008},
{0x300000000,0x600000005,0x800000007,0x800000008},
{0x300000001,0x600000005,0x800000007,0x800000008},
{0x100000000,0x500000003,0x700000006,0x800000008},
{0x300000002,0x600000005,0x800000007,0x800000008},
{0x200000000,0x500000003,0x700000006,0x800000008},
{0x200000001,0x500000003,0x700000006,0x800000008},
{0x100000000,0x300000002,0x600000005,0x800000007},
{0x500000004,0x700000006,0x800000008,0x800000008},
{0x400000000,0x600000005,0x800000007,0x800000008},
{0x400000001,0x600000005,0x800000007,0x800000008},
{0x100000000,0x500000004,0x700000006,0x800000008},
{0x400000002,0x600000005,0x800000007,0x800000008},
{0x200000000,0x500000004,0x700000006,0x800000008},
{0x200000001,0x500000004,0x700000006,0x800000008},
{0x100000000,0x400000002,0x600000005,0x800000007},
{0x400000003,0x600000005,0x800000007,0x800000008},
{0x300000000,0x500000004,0x700000006,0x800000008},
{0x300000001,0x500000004,0x700000006,0x800000008},
{0x100000000,0x400000003,0x600000005,0x800000007},
{0x300000002,0x500000004,0x700000006,0x800000008},
{0x200000000,0x400000003,0x600000005,0x800000007},
{0x200000001,0x400000003,0x600000005,0x800000007},
{0x100000000,0x300000002,0x500000004,0x700000006},
};
const unsigned int popcount[256] = {
0,1,1,2,1,2,2,3,
1,2,2,3,2,3,3,4,
1,2,2,3,2,3,3,4,
2,3,3,4,3,4,4,5,
1,2,2,3,2,3,3,4,
2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,
3,4,4,5,4,5,5,6,
1,2,2,3,2,3,3,4,
2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,
3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,
3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,
4,5,5,6,5,6,6,7,
1,2,2,3,2,3,3,4,
2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,
3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,
3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,
4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,
3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,
4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,
4,5,5,6,5,6,6,7,
4,5,5,6,5,6,6,7,
5,6,6,7,6,7,7,8
};
void aigis_poly_uniform(aigis_poly *a, unsigned char *buf)
{
	unsigned int ctr, pos, offset;
	uint32_t t;
	__m256i tmp, index;
	__m128i tmp0, index0;
	__m256i q8x = _mm256_set1_epi32(PARAM_Q);
#if QBITS == 21
	__m256i mask = _mm256_set1_epi32(0x1FFFFF);
#elif QBITS == 22
	__m256i mask = _mm256_set1_epi32(0x3FFFFF);
#endif
	ctr = pos = 0;

	while (ctr + 8 < PARAM_N)
	{
		tmp = _mm256_set_epi32(*(uint32_t*)&buf[pos + 21], *(uint32_t*)&buf[pos + 18], *(uint32_t*)&buf[pos + 15], *(uint32_t*)&buf[pos + 12], *(uint32_t*)&buf[pos + 9], *(uint32_t*)&buf[pos + 6], *(uint32_t*)&buf[pos + 3], *(uint32_t*)&buf[pos]);

		tmp = _mm256_and_si256(tmp, mask);

		index = _mm256_cmpgt_epi32(q8x, tmp);

		t = _mm256_movemask_ps(_mm256_castsi256_ps(index));

		index = _mm256_loadu_si256((__m256i *)idx[t]);
		tmp = _mm256_permutevar8x32_epi32(tmp, index); // permute good values to the bottom of tmp2

		_mm256_storeu_si256((__m256i *)&a->coeffs[ctr], tmp);

		offset = popcount[t];

		ctr += offset;

		pos += 24;
	}

	while (ctr + 4 < PARAM_N)
	{
		tmp0 = _mm_set_epi32(*(uint32_t*)&buf[pos + 9], *(uint32_t*)&buf[pos + 6], *(uint32_t*)&buf[pos + 3], *(uint32_t*)&buf[pos]);

		tmp0 = _mm_and_si128(tmp0, _mm256_castsi256_si128(mask));

		index0 = _mm_cmpgt_epi32(_mm256_castsi256_si128(q8x), tmp0);

		t = _mm_movemask_ps(_mm_castsi128_ps(index0));

		index0 = _mm_loadu_si128((__m128i *)idx[t]);
		tmp0 = _mm_castps_si128(_mm_permutevar_ps(_mm_castsi128_ps(tmp0), index0)); // permute good values to the bottom of tmp2

		_mm_storeu_si128((__m128i *)&a->coeffs[ctr], tmp0);

		offset = popcount[t];

		ctr += offset;

		pos += 12;
	}

	while (ctr < PARAM_N) {
		t = buf[pos++];
		t |= (uint32_t)buf[pos++] << 8;
		t |= (uint32_t)buf[pos++] << 16;
#if QBITS == 21
		t &= 0x1FFFFF;
#elif QBITS == 22
		t &= 0x3FFFFF;
#endif

		if (t < PARAM_Q)
			a->coeffs[ctr++] = t;
	}
}


static void rej_eta1(uint32_t *a,const unsigned char *buf)
{
#if ETA1 > 3
#error "rej_eta1() assumes ETA1 <= 3"
#endif
    unsigned int ctr, pos;
    unsigned char t[8];
    
    ctr = pos = 0;
#if ETA1 == 1
	do {
		t[0] = buf[pos] & 0x03;
		t[1] = (buf[pos] >> 2) & 0x03;
		t[2] = (buf[pos] >> 4) & 0x03;
		t[3] = (buf[pos++] >> 6) & 0x03;

		if (t[0] <= 2 * ETA1)
			a[ctr++] = PARAM_Q + ETA1 - t[0];
		if (t[1] <= 2 * ETA1)
			a[ctr++] = PARAM_Q + ETA1 - t[1];
		if (t[2] <= 2 * ETA1)
			a[ctr++] = PARAM_Q + ETA1 - t[2];
		if (t[3] <= 2 * ETA1)
			a[ctr++] = PARAM_Q + ETA1 - t[3];
	} while (ctr < PARAM_N - 4);

	do {
		t[0] = buf[pos] & 0x03;
		t[1] = (buf[pos] >> 2) & 0x03;
		t[2] = (buf[pos] >> 4) & 0x03;
		t[3] = (buf[pos++] >> 6) & 0x03;

		if (t[0] <= 2 * ETA1)
			a[ctr++] = PARAM_Q + ETA1 - t[0];
		if (t[1] <= 2 * ETA1 && ctr < PARAM_N)
			a[ctr++] = PARAM_Q + ETA1 - t[1];
		if (t[2] <= 2 * ETA1 && ctr < PARAM_N)
			a[ctr++] = PARAM_Q + ETA1 - t[2];
		if (t[3] <= 2 * ETA1 && ctr < PARAM_N)
			a[ctr++] = PARAM_Q + ETA1 - t[3];
	} while (ctr < PARAM_N);
#else
    do {
        
        t[0] = buf[pos] & 0x07;
        t[1] = (buf[pos] >> 3) & 0x07;
        t[2] = (buf[pos] >> 6) | ((buf[pos+1]& 0x1)<<2);
        t[3] = (buf[++pos] >> 1) & 0x07;
        t[4] = (buf[pos] >> 4) & 0x07;
        t[5] = (buf[pos] >> 7) | ((buf[pos+1]& 0x3)<<1);
        t[6] = (buf[++pos] >> 2) & 0x07;
        t[7] = buf[pos++] >> 5;
        
        
        if(t[0] <= 2*ETA1)
            a[ctr++] = PARAM_Q + ETA1 - t[0];
        if(t[1] <= 2*ETA1)
            a[ctr++] = PARAM_Q + ETA1 - t[1];
        if(t[2] <= 2*ETA1)
            a[ctr++] = PARAM_Q + ETA1 - t[2];
        if(t[3] <= 2*ETA1)
            a[ctr++] = PARAM_Q + ETA1 - t[3];
        if(t[4] <= 2*ETA1)
            a[ctr++] = PARAM_Q + ETA1 - t[4];
        if(t[5] <= 2*ETA1)
            a[ctr++] = PARAM_Q + ETA1 - t[5];
        if(t[6] <= 2*ETA1)
            a[ctr++] = PARAM_Q + ETA1 - t[6];
        if(t[7] <= 2*ETA1)
            a[ctr++] = PARAM_Q + ETA1 - t[7];
        
    }while(ctr < PARAM_N-8);
    
    
    do {
        
        t[0] = buf[pos] & 0x07;
        t[1] = (buf[pos] >> 3) & 0x07;
        t[2] = (buf[pos] >> 6) | ((buf[pos+1]& 0x1)<<2);
        t[3] = (buf[++pos] >> 1) & 0x07;
        t[4] = (buf[pos] >> 4) & 0x07;
        t[5] = (buf[pos] >> 7) | ((buf[pos+1]& 0x3)<<1);
        t[6] = (buf[++pos] >> 2) & 0x07;
        t[7] = buf[pos++] >> 5;
        
        
        if(t[0] <= 2*ETA1)
            a[ctr++] = PARAM_Q + ETA1 - t[0];
        if(t[1] <= 2*ETA1 && ctr < PARAM_N)
            a[ctr++] = PARAM_Q + ETA1 - t[1];
        if(t[2] <= 2*ETA1 && ctr < PARAM_N)
            a[ctr++] = PARAM_Q + ETA1 - t[2];
        if(t[3] <= 2*ETA1 && ctr < PARAM_N)
            a[ctr++] = PARAM_Q + ETA1 - t[3];
        if(t[4] <= 2*ETA1 && ctr < PARAM_N)
            a[ctr++] = PARAM_Q + ETA1 - t[4];
        if(t[5] <= 2*ETA1 && ctr < PARAM_N)
            a[ctr++] = PARAM_Q + ETA1 - t[5];
        if(t[6] <= 2*ETA1 && ctr < PARAM_N)
            a[ctr++] = PARAM_Q + ETA1 - t[6];
        if(t[7] <= 2*ETA1 && ctr < PARAM_N)
            a[ctr++] = PARAM_Q + ETA1 - t[7];
        
    }while(ctr < PARAM_N);
#endif
}

static unsigned int rej_eta2(uint32_t *a,unsigned int len,const unsigned char *buf)
{
#if ETA2 >7 || ETA2 <3
#error "rej_eta2() assumes 3 <= ETA2 <=7"
#endif
    
    unsigned int ctr=0, pos=0;
    unsigned char t0, t1;
    
    do {
#if ETA2 == 3
        t0 = buf[pos] & 0x07;
        t1 = buf[pos++] >> 5;
#else
        t0 = buf[pos] & 0x0F;
        t1 = buf[pos++] >> 4;
#endif
        if(t0 <= 2*ETA2)
            a[ctr++] = PARAM_Q + ETA2 - t0;
        if(t1 <= 2*ETA2)
            a[ctr++] = PARAM_Q + ETA2 - t1;
    }while(ctr < len-2);
    
    
    do {
#if ETA2 == 3
        t0 = buf[pos] & 0x07;
        t1 = buf[pos++] >> 5;
#else
        t0 = buf[pos] & 0x0F;
        t1 = buf[pos++] >> 4;
#endif
        
        if(t0 <= 2*ETA2)
            a[ctr++] = PARAM_Q + ETA2 - t0;
        if(t1 <= 2*ETA2 && ctr<len)
            a[ctr++] = PARAM_Q + ETA2 - t1;
    }while(ctr < len);
    
    return pos;
}



void aigis_poly_uniform_eta1(aigis_poly *a,
                      const unsigned char seed[SEEDBYTES], 
                      unsigned char nonce)
{
#if ETA1>3 
#error "rej_eta1() assumes ETA1 <=3"
#endif
 
#ifdef USE_AES
  aigis_aes256ctr_ctx state;
  unsigned char outbuf[256]; 

  aigis_aes256ctr_init(&state, seed, nonce);
  aigis_aes256ctr_squeezeblocks(outbuf, 2, &state); /* Probability we need more than 2 blocks: < 2^{-131}.*/
#else
  unsigned char inbuf[SEEDBYTES + 1];
  unsigned char outbuf[2 * SHAKE256_RATE];

  memcpy(inbuf, seed, SEEDBYTES);
  inbuf[SEEDBYTES] = nonce;

  shake256(outbuf, sizeof(outbuf),inbuf, SEEDBYTES + 1);
#endif
  rej_eta1(a->coeffs,outbuf);
}

void aigis_poly_uniform_eta1_3x(aigis_poly *a0,
                          aigis_poly *a1,
                          aigis_poly *a2,
                          const unsigned char seed[SEEDBYTES],
                          unsigned char nonce0,
                          unsigned char nonce1,
                          unsigned char nonce2)
{
#if ETA1>3 
#error "rej_eta1() assumes ETA1 <=3"
#endif
 
  unsigned int i;
  unsigned char inbuf[4][SEEDBYTES + 1];
  unsigned char outbuf[4][2*SHAKE256_RATE];

  for(i = 0; i < SEEDBYTES; ++i) {
        inbuf[0][i] = seed[i];
        inbuf[1][i] = seed[i];
        inbuf[2][i] = seed[i];
        inbuf[3][i] = seed[i];
    }
   inbuf[0][SEEDBYTES] = nonce0;
   inbuf[1][SEEDBYTES] = nonce1;
   inbuf[2][SEEDBYTES] = nonce2;
   inbuf[3][SEEDBYTES] = nonce2+1;

  shake256_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], 2*SHAKE256_RATE, inbuf[0], inbuf[1], inbuf[2], inbuf[3],SEEDBYTES + 1);

  rej_eta1(a0->coeffs,outbuf[0]);
  rej_eta1(a1->coeffs,outbuf[1]);
  rej_eta1(a2->coeffs,outbuf[2]);
}

void aigis_poly_uniform_eta1_4x(aigis_poly *a0,
                          aigis_poly *a1,
                          aigis_poly *a2,
                          aigis_poly *a3,
                          const unsigned char seed[SEEDBYTES],
                          unsigned char nonce0,
                          unsigned char nonce1,
                          unsigned char nonce2,
                          unsigned char nonce3)
{
    unsigned int i;
    unsigned char inbuf[4][SEEDBYTES + 1];
    unsigned char outbuf[4][2*SHAKE256_RATE];
    
    for(i = 0; i < SEEDBYTES; ++i) {
        inbuf[0][i] = seed[i];
        inbuf[1][i] = seed[i];
        inbuf[2][i] = seed[i];
        inbuf[3][i] = seed[i];
    }
    inbuf[0][SEEDBYTES] = nonce0;
    inbuf[1][SEEDBYTES] = nonce1;
    inbuf[2][SEEDBYTES] = nonce2;
    inbuf[3][SEEDBYTES] = nonce3;
    
    shake256_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], 2*SHAKE256_RATE, inbuf[0], inbuf[1], inbuf[2], inbuf[3],SEEDBYTES + 1);
    
    rej_eta1(a0->coeffs,outbuf[0]);
    rej_eta1(a1->coeffs,outbuf[1]);
    rej_eta1(a2->coeffs,outbuf[2]);
    rej_eta1(a3->coeffs,outbuf[3]);
}


void aigis_poly_uniform_eta2_4x(aigis_poly *a0,
                          aigis_poly *a1,
                          aigis_poly *a2,
                          aigis_poly *a3,
                          const unsigned char seed[SEEDBYTES],
                          unsigned char nonce0,
                          unsigned char nonce1,
                          unsigned char nonce2,
                          unsigned char nonce3)
{
    unsigned int i,len;
    unsigned char inbuf[4][SEEDBYTES + 1];
#if ETA2==3
    unsigned char outbuf[4][2*SHAKE256_RATE];
    len = 2*SHAKE256_RATE;
#elif ETA2==5 || ETA2 == 4
	unsigned char outbuf[4][3 * SHAKE256_RATE];
	len = 3 * SHAKE256_RATE;
#else
	unsigned char outbuf[4][2 * SHAKE256_RATE];
	len = 2 * SHAKE256_RATE;
#endif
    
    for(i = 0; i < SEEDBYTES; ++i) {
        inbuf[0][i] = seed[i];
        inbuf[1][i] = seed[i];
        inbuf[2][i] = seed[i];
        inbuf[3][i] = seed[i];
    }
    inbuf[0][SEEDBYTES] = nonce0;
    inbuf[1][SEEDBYTES] = nonce1;
    inbuf[2][SEEDBYTES] = nonce2;
    inbuf[3][SEEDBYTES] = nonce3;
    
    shake256_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], len, inbuf[0], inbuf[1], inbuf[2], inbuf[3],SEEDBYTES + 1);
    
    rej_eta2(a0->coeffs,PARAM_N,outbuf[0]);
    rej_eta2(a1->coeffs,PARAM_N,outbuf[1]);
    rej_eta2(a2->coeffs,PARAM_N,outbuf[2]);
    rej_eta2(a3->coeffs,PARAM_N,outbuf[3]);
}
void aigis_poly_uniform_eta2_2x(aigis_poly *a0,
	aigis_poly *a1,
	const unsigned char seed[SEEDBYTES],
	unsigned char nonce0,
	unsigned char nonce1)
{
	unsigned int i, len;
	unsigned char inbuf[4][SEEDBYTES + 1];
#if ETA2==3
	unsigned char outbuf[4][2 * SHAKE256_RATE];
	len = 2 * SHAKE256_RATE;
#elif ETA2==5 || ETA2 == 4
	unsigned char outbuf[4][3 * SHAKE256_RATE];
	len = 3 * SHAKE256_RATE;
#else
	unsigned char outbuf[4][2 * SHAKE256_RATE];
	len = 2 * SHAKE256_RATE;
#endif

	for (i = 0; i < SEEDBYTES; ++i) {
		inbuf[0][i] = seed[i];
		inbuf[1][i] = seed[i];
		inbuf[2][i] = seed[i];
		inbuf[3][i] = seed[i];
	}
	inbuf[0][SEEDBYTES] = nonce0;
	inbuf[1][SEEDBYTES] = nonce1;
	inbuf[2][SEEDBYTES] = nonce1+1;
	inbuf[3][SEEDBYTES] = nonce1+2;

	shake256_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], len, inbuf[0], inbuf[1], inbuf[2], inbuf[3], SEEDBYTES + 1);

	rej_eta2(a0->coeffs, PARAM_N, outbuf[0]);
	rej_eta2(a1->coeffs, PARAM_N, outbuf[1]);
}

void aigis_poly_uniform_eta2(aigis_poly *a,
                      const unsigned char seed[SEEDBYTES], 
                      unsigned char nonce)
{
  //unsigned int i;
  //unsigned char inbuf[SEEDBYTES + 1];

#ifdef USE_AES
  aigis_aes256ctr_ctx state;
#if ETA2==3
  unsigned char outbuf[256];
#elif ETA2==4
  unsigned char outbuf[512];
#elif ETA2 ==5
  unsigned char outbuf[384];
#else 
  unsigned char outbuf[256];
#endif
  aigis_aes256ctr_init(&state, seed, nonce);
  aigis_aes256ctr_squeezeblocks(outbuf, sizeof(outbuf) / 128, &state);
  rej_eta2(a->coeffs, PARAM_N, outbuf);
#else
  for (i = 0; i < SEEDBYTES; ++i)
	  inbuf[i] = seed[i];
  inbuf[SEEDBYTES] = nonce;

#if ETA2==3
  unsigned char outbuf[2 * SHAKE256_RATE];
  shake256(outbuf, sizeof(outbuf), inbuf, SEEDBYTES + 1);
  rej_eta2(a->coeffs, PARAM_N, outbuf);
#elif ETA2==5
  unsigned int pos;
  uint64_t state[25];
  unsigned char outbuf[3 * SHAKE256_RATE];

  shake256_absorb_aigis_mulan(state, inbuf, SEEDBYTES + 1); //fzhang
  shake256_squeezeblocks_aigis_mulan(outbuf, 2, state); //fzhang
  pos = rej_eta2(a->coeffs, 223, outbuf);

  /* Probability we need more than 2 blocks to generate 223 elements: < 2^{-133}.*/
  /* Probability we need more than 85 bytes to generate 33 elements: < 2^{-133}.*/

  if (2 * SHAKE256_RATE - pos < 85)
	  shake256_squeezeblocks_aigis_mulan(&outbuf[2 * SHAKE256_RATE], 1, state); //fzhang

  rej_eta2(&a->coeffs[223], 33, &outbuf[pos]);
#else
#if ETA2 == 4
  unsigned char outbuf[3 * SHAKE256_RATE];
#else 
  unsigned char outbuf[2 * SHAKE256_RATE];
#endif
  shake256(outbuf, sizeof(outbuf), inbuf, SEEDBYTES + 1);
  rej_eta2(a->coeffs, PARAM_N, outbuf);
#endif
#endif

}

void aigis_poly_uniform_gamma1m1(aigis_poly *a,
                           const unsigned char seed[SEEDBYTES + CRHBYTES],
                           uint16_t nonce)
{
#if GAMMA1 != 131072
#error "aigis_poly_uniform_gamma1m1() assumes GAMMA1 == 131072"
#endif

  unsigned int ctr=0, pos=0;
  uint32_t t0, t1;

#ifdef USE_AES
  aigis_aes256ctr_ctx state;
  unsigned char outbuf[640];
  aigis_aes256ctr_init(&state, seed, nonce);
  aigis_aes256ctr_squeezeblocks(outbuf, sizeof(outbuf) / 128, &state);
#else
  unsigned char outbuf[5 * SHAKE256_RATE];
  unsigned char inbuf[SEEDBYTES + CRHBYTES + 2];

  memcpy(inbuf, seed, SEEDBYTES + CRHBYTES);

  inbuf[SEEDBYTES + CRHBYTES] = nonce & 0xFF;
  inbuf[SEEDBYTES + CRHBYTES + 1] = nonce >> 8;

  shake256(outbuf, sizeof(outbuf), inbuf, SEEDBYTES + CRHBYTES + 2);
#endif
  
  do{
    t0  = outbuf[pos];
    t0 |= (uint32_t)outbuf[pos + 1] << 8;
    t0 |= (uint32_t)outbuf[pos + 2] << 16;

	t1  = outbuf[pos + 2] >> 4;
    t1 |= (uint32_t)outbuf[pos + 3] << 4;
    t1 |= (uint32_t)outbuf[pos + 4] << 12;

	t0 &= 0x3FFFF;
	t1 &= 0x3FFFF;

    if(t0 <= 2*GAMMA1)
    		a->coeffs[ctr++] = PARAM_Q + GAMMA1 - 1  - t0;
    if(t1 <= 2*GAMMA1)
		a->coeffs[ctr++] = PARAM_Q + GAMMA1 - 1  - t1;
     
    pos += 5;
  }while(ctr < PARAM_N-2);
  
  
  do{
    t0  = outbuf[pos];
    t0 |= (uint32_t)outbuf[pos + 1] << 8;
    t0 |= (uint32_t)outbuf[pos + 2] << 16;

	t1  = outbuf[pos + 2] >> 4;
    t1 |= (uint32_t)outbuf[pos + 3] << 4;
    t1 |= (uint32_t)outbuf[pos + 4] << 12;

	t0 &= 0x3FFFF;
	t1 &= 0x3FFFF;

    if(t0 <= 2*GAMMA1)
    		a->coeffs[ctr++] = PARAM_Q + GAMMA1 - 1 - t0;
    if(t1 <= 2*GAMMA1 && ctr< PARAM_N)
		a->coeffs[ctr++] = PARAM_Q + GAMMA1 - 1 - t1;
      
    pos += 5;
  }while(ctr < PARAM_N);

}

static void uniform_gamma1m1(aigis_poly *a,unsigned char *buf)
{
#if GAMMA1 != 131072
#error "aigis_poly_uniform_gamma1m1() assumes GAMMA1 == 131072"
#endif

  unsigned int  ctr=0, pos=0;
  uint32_t t0, t1;
  do{
    t0  = buf[pos];
    t0 |= (uint32_t)buf[pos + 1] << 8;
    t0 |= (uint32_t)buf[pos + 2] << 16;


	t1  = buf[pos + 2] >> 4;
    t1 |= (uint32_t)buf[pos + 3] << 4;
    t1 |= (uint32_t)buf[pos + 4] << 12;

	t0 &= 0x3FFFF;
	t1 &= 0x3FFFF;

    if(t0 <= 2*GAMMA1)
    		a->coeffs[ctr++] = PARAM_Q + GAMMA1 - 1 - t0;
    if(t1 <= 2*GAMMA1)
		a->coeffs[ctr++] = PARAM_Q + GAMMA1 - 1 - t1;

    pos += 5;
  }while(ctr < PARAM_N-2);
  
  do{
    t0  = buf[pos];
    t0 |= (uint32_t)buf[pos + 1] << 8;
    t0 |= (uint32_t)buf[pos + 2] << 16;


	t1  = buf[pos + 2] >> 4;
    t1 |= (uint32_t)buf[pos + 3] << 4;
    t1 |= (uint32_t)buf[pos + 4] << 12;

	t0 &= 0x3FFFF;
	t1 &= 0x3FFFF;

    if(t0 <= 2*GAMMA1)
    		a->coeffs[ctr++] = PARAM_Q + GAMMA1 - 1 - t0;
    if(t1 <= 2*GAMMA1 && ctr< PARAM_N)
		a->coeffs[ctr++] = PARAM_Q + GAMMA1 - 1 - t1;

    pos += 5;
  }while(ctr < PARAM_N);
}
void aigis_poly_uniform_gamma1m1_3x(aigis_poly *a0,
                              aigis_poly *a1,
                              aigis_poly *a2,
                              const unsigned char seed[SEEDBYTES + CRHBYTES],
                              uint16_t nonce0,
                              uint16_t nonce1,
                              uint16_t nonce2)
{
  unsigned int i;
  unsigned char inbuf[4][SEEDBYTES + CRHBYTES + 2];
  unsigned char outbuf[4][640];

  for(i = 0; i < SEEDBYTES + CRHBYTES; ++i) {
    inbuf[0][i] = seed[i];
    inbuf[1][i] = seed[i];
    inbuf[2][i] = seed[i];
    inbuf[3][i] = seed[i];
  }
  inbuf[0][SEEDBYTES + CRHBYTES] = nonce0 & 0xFF;
  inbuf[0][SEEDBYTES + CRHBYTES + 1] = nonce0 >> 8;
  inbuf[1][SEEDBYTES + CRHBYTES] = nonce1 & 0xFF;
  inbuf[1][SEEDBYTES + CRHBYTES + 1] = nonce1 >> 8;
  inbuf[2][SEEDBYTES + CRHBYTES] = nonce2 & 0xFF;
  inbuf[2][SEEDBYTES + CRHBYTES + 1] = nonce2 >> 8;
  inbuf[3][SEEDBYTES + CRHBYTES] = (nonce2+1) & 0xFF;
  inbuf[3][SEEDBYTES + CRHBYTES + 1] = (nonce2+1) >> 8;

  shake256_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], 640, inbuf[0], inbuf[1], inbuf[2], inbuf[3],SEEDBYTES + CRHBYTES + 2);
  
  uniform_gamma1m1(a0,outbuf[0]);
  uniform_gamma1m1(a1,outbuf[1]);
  uniform_gamma1m1(a2,outbuf[2]);
}
void aigis_poly_uniform_gamma1m1_4x(aigis_poly *a0,
                              aigis_poly *a1,
                              aigis_poly *a2,
                              aigis_poly *a3,
                              const unsigned char seed[SEEDBYTES + CRHBYTES],
                              uint16_t nonce0,
                              uint16_t nonce1,
                              uint16_t nonce2,
                              uint16_t nonce3)
{
  unsigned int i;
  unsigned char inbuf[4][SEEDBYTES + CRHBYTES + 2];
  unsigned char outbuf[4][640];

  for(i = 0; i < SEEDBYTES + CRHBYTES; ++i) {
    inbuf[0][i] = seed[i];
    inbuf[1][i] = seed[i];
    inbuf[2][i] = seed[i];
    inbuf[3][i] = seed[i];
  }
  inbuf[0][SEEDBYTES + CRHBYTES] = nonce0 & 0xFF;
  inbuf[0][SEEDBYTES + CRHBYTES + 1] = nonce0 >> 8;
  inbuf[1][SEEDBYTES + CRHBYTES] = nonce1 & 0xFF;
  inbuf[1][SEEDBYTES + CRHBYTES + 1] = nonce1 >> 8;
  inbuf[2][SEEDBYTES + CRHBYTES] = nonce2 & 0xFF;
  inbuf[2][SEEDBYTES + CRHBYTES + 1] = nonce2 >> 8;
  inbuf[3][SEEDBYTES + CRHBYTES] = nonce3 & 0xFF;
  inbuf[3][SEEDBYTES + CRHBYTES + 1] = nonce3 >> 8;

  shake256_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], 640, inbuf[0], inbuf[1], inbuf[2], inbuf[3],SEEDBYTES + CRHBYTES + 2);
  
  uniform_gamma1m1(a0,outbuf[0]);
  uniform_gamma1m1(a1,outbuf[1]);
  uniform_gamma1m1(a2,outbuf[2]);
  uniform_gamma1m1(a3,outbuf[3]);
  
}

void aigis_polyeta1_pack(unsigned char *r, const aigis_poly *a) {
#if ETA1 > 3
#error "aigis_polyeta1_pack() assumes ETA1 <= 3"
#endif
    unsigned int i;
    unsigned char t[8];
#if ETA1 == 1
	for (i = 0; i < PARAM_N / 4; ++i) {
		t[0] = PARAM_Q + ETA1 - a->coeffs[4 * i + 0];
		t[1] = PARAM_Q + ETA1 - a->coeffs[4 * i + 1];
		t[2] = PARAM_Q + ETA1 - a->coeffs[4 * i + 2];
		t[3] = PARAM_Q + ETA1 - a->coeffs[4 * i + 3];
		r[i] = t[0] | (t[1] << 2) | (t[2] << 4) | (t[3] << 6);
	}
#else    
    for(i = 0; i < PARAM_N/8; ++i) {
        t[0] = PARAM_Q + ETA1 - a->coeffs[8*i+0];
        t[1] = PARAM_Q + ETA1 - a->coeffs[8*i+1];
        t[2] = PARAM_Q + ETA1 - a->coeffs[8*i+2];
        t[3] = PARAM_Q + ETA1 - a->coeffs[8*i+3];
        t[4] = PARAM_Q + ETA1 - a->coeffs[8*i+4];
        t[5] = PARAM_Q + ETA1 - a->coeffs[8*i+5];
        t[6] = PARAM_Q + ETA1 - a->coeffs[8*i+6];
        t[7] = PARAM_Q + ETA1 - a->coeffs[8*i+7];
        
        r[3*i+0]  = t[0];
        r[3*i+0] |= t[1] << 3;
        r[3*i+0] |= t[2] << 6;
        r[3*i+1]  = t[2] >> 2;
        r[3*i+1] |= t[3] << 1;
        r[3*i+1] |= t[4] << 4;
        r[3*i+1] |= t[5] << 7;
        r[3*i+2]  = t[5] >> 1;
        r[3*i+2] |= t[6] << 2;
        r[3*i+2] |= t[7] << 5;
    }
#endif
}
void aigis_polyeta2_pack(unsigned char *r, const aigis_poly *a) {
#if ETA2 > 7
#error "aigis_polyeta2_pack() assumes ETA2 <= 7"
#endif
    unsigned int i;
    unsigned char t[8];
    
#if ETA2 <= 3
    for(i = 0; i < PARAM_N/8; ++i) {
        t[0] = PARAM_Q + ETA2 - a->coeffs[8*i+0];
        t[1] = PARAM_Q + ETA2 - a->coeffs[8*i+1];
        t[2] = PARAM_Q + ETA2 - a->coeffs[8*i+2];
        t[3] = PARAM_Q + ETA2 - a->coeffs[8*i+3];
        t[4] = PARAM_Q + ETA2 - a->coeffs[8*i+4];
        t[5] = PARAM_Q + ETA2 - a->coeffs[8*i+5];
        t[6] = PARAM_Q + ETA2 - a->coeffs[8*i+6];
        t[7] = PARAM_Q + ETA2 - a->coeffs[8*i+7];
        
        r[3*i+0]  = t[0];
        r[3*i+0] |= t[1] << 3;
        r[3*i+0] |= t[2] << 6;
        r[3*i+1]  = t[2] >> 2;
        r[3*i+1] |= t[3] << 1;
        r[3*i+1] |= t[4] << 4;
        r[3*i+1] |= t[5] << 7;
        r[3*i+2]  = t[5] >> 1;
        r[3*i+2] |= t[6] << 2;
        r[3*i+2] |= t[7] << 5;
    }
#else
    for(i = 0; i < PARAM_N/2; ++i) {
        t[0] = PARAM_Q + ETA2 - a->coeffs[2*i+0];
        t[1] = PARAM_Q + ETA2 - a->coeffs[2*i+1];
        r[i] = t[0] | (t[1] << 4);
    }
#endif
}

void aigis_polyeta1_unpack(aigis_poly *r, const unsigned char *a)
{
#if ETA1 > 3
#error "aigis_polyeta1_unpack() assumes ETA1 <= 3"
#endif
    
    unsigned int i;
#if ETA1 == 1
	for (i = 0; i < PARAM_N / 4; ++i) {
		r->coeffs[4 * i + 0] = a[i] & 0x03;
		r->coeffs[4 * i + 1] = (a[i] >> 2) & 0x03;
		r->coeffs[4 * i + 2] = (a[i] >> 4) & 0x03;
		r->coeffs[4 * i + 3] = (a[i] >> 6) & 0x03;

		r->coeffs[4 * i + 0] = PARAM_Q + ETA1 - r->coeffs[4 * i + 0];
		r->coeffs[4 * i + 1] = PARAM_Q + ETA1 - r->coeffs[4 * i + 1];
		r->coeffs[4 * i + 2] = PARAM_Q + ETA1 - r->coeffs[4 * i + 2];
		r->coeffs[4 * i + 3] = PARAM_Q + ETA1 - r->coeffs[4 * i + 3];
	}
#else    
    for(i = 0; i < PARAM_N/8; ++i) {
        r->coeffs[8*i+0] = a[3*i+0] & 0x07;
        r->coeffs[8*i+1] = (a[3*i+0] >> 3) & 0x07;
        r->coeffs[8*i+2] = (a[3*i+0] >> 6) | ((a[3*i+1] & 0x01) << 2);
        r->coeffs[8*i+3] = (a[3*i+1] >> 1) & 0x07;
        r->coeffs[8*i+4] = (a[3*i+1] >> 4) & 0x07;
        r->coeffs[8*i+5] = (a[3*i+1] >> 7) | ((a[3*i+2] & 0x03) << 1);
        r->coeffs[8*i+6] = (a[3*i+2] >> 2) & 0x07;
        r->coeffs[8*i+7] = (a[3*i+2] >> 5);
        
        r->coeffs[8*i+0] = PARAM_Q + ETA1 - r->coeffs[8*i+0];
        r->coeffs[8*i+1] = PARAM_Q + ETA1 - r->coeffs[8*i+1];
        r->coeffs[8*i+2] = PARAM_Q + ETA1 - r->coeffs[8*i+2];
        r->coeffs[8*i+3] = PARAM_Q + ETA1 - r->coeffs[8*i+3];
        r->coeffs[8*i+4] = PARAM_Q + ETA1 - r->coeffs[8*i+4];
        r->coeffs[8*i+5] = PARAM_Q + ETA1 - r->coeffs[8*i+5];
        r->coeffs[8*i+6] = PARAM_Q + ETA1 - r->coeffs[8*i+6];
        r->coeffs[8*i+7] = PARAM_Q + ETA1 - r->coeffs[8*i+7];
    }
#endif
}
void aigis_polyeta2_unpack(aigis_poly *r, const unsigned char *a)
{
#if ETA2 > 7
#error "aigis_polyeta2_unpack() assumes ETA2 <= 7"
#endif
    
    unsigned int i;
#if ETA2 <= 3
    for(i = 0; i < PARAM_N/8; ++i) {
        r->coeffs[8*i+0] = a[3*i+0] & 0x07;
        r->coeffs[8*i+1] = (a[3*i+0] >> 3) & 0x07;
        r->coeffs[8*i+2] = (a[3*i+0] >> 6) | ((a[3*i+1] & 0x01) << 2);
        r->coeffs[8*i+3] = (a[3*i+1] >> 1) & 0x07;
        r->coeffs[8*i+4] = (a[3*i+1] >> 4) & 0x07;
        r->coeffs[8*i+5] = (a[3*i+1] >> 7) | ((a[3*i+2] & 0x03) << 1);
        r->coeffs[8*i+6] = (a[3*i+2] >> 2) & 0x07;
        r->coeffs[8*i+7] = (a[3*i+2] >> 5);
        
        r->coeffs[8*i+0] = PARAM_Q + ETA2 - r->coeffs[8*i+0];
        r->coeffs[8*i+1] = PARAM_Q + ETA2 - r->coeffs[8*i+1];
        r->coeffs[8*i+2] = PARAM_Q + ETA2 - r->coeffs[8*i+2];
        r->coeffs[8*i+3] = PARAM_Q + ETA2 - r->coeffs[8*i+3];
        r->coeffs[8*i+4] = PARAM_Q + ETA2 - r->coeffs[8*i+4];
        r->coeffs[8*i+5] = PARAM_Q + ETA2 - r->coeffs[8*i+5];
        r->coeffs[8*i+6] = PARAM_Q + ETA2 - r->coeffs[8*i+6];
        r->coeffs[8*i+7] = PARAM_Q + ETA2 - r->coeffs[8*i+7];
    }
#else
    for(i = 0; i < PARAM_N/2; ++i) {
        r->coeffs[2*i+0] = a[i] & 0x0F;
        r->coeffs[2*i+1] = a[i] >> 4;
        r->coeffs[2*i+0] = PARAM_Q + ETA2 - r->coeffs[2*i+0];
        r->coeffs[2*i+1] = PARAM_Q + ETA2 - r->coeffs[2*i+1];
    }
#endif
}

void aigis_polyt1_pack(unsigned char *r, const aigis_poly *a)
{
#if QBITS - PARAM_D != 8
#error "aigis_polyt1_pack() assumes QBITS - PARAM_D == 8"
#endif
    
    unsigned int i;
    for(i = 0; i < PARAM_N; ++i)
        r[i]  =  a->coeffs[i];
}

void aigis_polyt1_unpack(aigis_poly *r, const unsigned char *a)
{
#if QBITS - PARAM_D != 8
#error "aigis_polyt1_unpack() assumes QBITS - PARAM_D == 8"
#endif
    
    unsigned int i;
    for(i = 0; i < PARAM_N; ++i)
        r->coeffs[i]  =  a[i];
}

void aigis_polyt0_pack(unsigned char *r, const aigis_poly *a) {
    
#if PARAM_D!=13 && PARAM_D!=14
#error "aigis_polyt0_unpack() assumes PARAM_D== 13 or 14"
#endif
    
    unsigned int i;
#if PARAM_D == 13
    uint32_t t[8];
    for(i = 0; i < PARAM_N/8; ++i) {
        t[0] = PARAM_Q + (1 << (PARAM_D-1)) - a->coeffs[8*i+0];
        t[1] = PARAM_Q + (1 << (PARAM_D-1)) - a->coeffs[8*i+1];
        t[2] = PARAM_Q + (1 << (PARAM_D-1)) - a->coeffs[8*i+2];
        t[3] = PARAM_Q + (1 << (PARAM_D-1)) - a->coeffs[8*i+3];
        t[4] = PARAM_Q + (1 << (PARAM_D-1)) - a->coeffs[8*i+4];
        t[5] = PARAM_Q + (1 << (PARAM_D-1)) - a->coeffs[8*i+5];
        t[6] = PARAM_Q + (1 << (PARAM_D-1)) - a->coeffs[8*i+6];
        t[7] = PARAM_Q + (1 << (PARAM_D-1)) - a->coeffs[8*i+7];
        
        r[13*i+0]   =  t[0];
        r[13*i+1]   =  t[0] >> 8;
        r[13*i+1]  |=  t[1] << 5;
        r[13*i+2]   =  t[1] >> 3;
        r[13*i+3]   =  t[1] >> 11;
        r[13*i+3]  |=  t[2] << 2;
        r[13*i+4]   =  t[2] >> 6;
        r[13*i+4]  |=  t[3] << 7;
        r[13*i+5]   =  t[3] >> 1;
        r[13*i+6]   =  t[3] >> 9;
        r[13*i+6]  |=  t[4] << 4;
        r[13*i+7]   =  t[4] >> 4;
        r[13*i+8]   =  t[4] >> 12;
        r[13*i+8]  |=  t[5] << 1;
        r[13*i+9]   =  t[5] >> 7;
        r[13*i+9]  |=  t[6] << 6;
        r[13*i+10]  =  t[6] >> 2;
        r[13*i+11]  =  t[6] >> 10;
        r[13*i+11] |=  t[7] << 3;
        r[13*i+12]  =  t[7] >> 5;
    }
#elif PARAM_D==14
    uint32_t t[4];
    for(i = 0; i < PARAM_N/4; ++i) {
        t[0] = PARAM_Q + (1 << (PARAM_D-1)) - a->coeffs[4*i+0];
        t[1] = PARAM_Q + (1 << (PARAM_D-1)) - a->coeffs[4*i+1];
        t[2] = PARAM_Q + (1 << (PARAM_D-1)) - a->coeffs[4*i+2];
        t[3] = PARAM_Q + (1 << (PARAM_D-1)) - a->coeffs[4*i+3];
        
        r[7*i+0]  =  t[0];
        r[7*i+1]  =  t[0] >> 8;
        r[7*i+1] |=  t[1] << 6;
        r[7*i+2]  =  t[1] >> 2;
        r[7*i+3]  =  t[1] >> 10;
        r[7*i+3] |=  t[2] << 4;
        r[7*i+4]  =  t[2] >> 4;
        r[7*i+5]  =  t[2] >> 12;
        r[7*i+5] |=  t[3] << 2;
        r[7*i+6]  =  t[3] >> 6;
    }
#endif
}

void aigis_polyt0_unpack(aigis_poly *r, const unsigned char *a)
{
#if PARAM_D!=13 && PARAM_D!=14
#error "aigis_polyt0_unpack() assumes PARAM_D== 13 or 14"
#endif
    
    unsigned int i;
#if PARAM_D==13
    for(i = 0; i < PARAM_N/8; ++i) {
        
        r->coeffs[8*i+0]  = a[13*i+0];
        r->coeffs[8*i+0] |= (uint32_t)(a[13*i+1] & 0x1F)<<8;
        
        r->coeffs[8*i+1]  = a[13*i+1]>>5;
        r->coeffs[8*i+1] |= (uint32_t)a[13*i+2]<< 3;
        r->coeffs[8*i+1] |= (uint32_t)(a[13*i+3] & 0x3)<< 11;
        
        r->coeffs[8*i+2]  = a[13*i+3]>>2;
        r->coeffs[8*i+2] |= (uint32_t)(a[13*i+4] & 0x7F)<< 6;
        
        r->coeffs[8*i+3]  = a[13*i+4]>>7;
        r->coeffs[8*i+3] |= (uint32_t)a[13*i+5]<< 1;
        r->coeffs[8*i+3] |= (uint32_t)(a[13*i+6] & 0x0F)<< 9;
        
        r->coeffs[8*i+4]  = a[13*i+ 6]>>4;
        r->coeffs[8*i+4] |= (uint32_t)a[13*i+ 7]<< 4;
        r->coeffs[8*i+4] |= (uint32_t)(a[13*i+8] & 0x01)<< 12;
        
        r->coeffs[8*i+5]  = a[13*i+8]>>1;
        r->coeffs[8*i+5] |= (uint32_t)(a[13*i+9] & 0x3F)<< 7;
        
        r->coeffs[8*i+6]  = a[13*i+9]>>6;
        r->coeffs[8*i+6] |= (uint32_t)a[13*i+10]<< 2;
        r->coeffs[8*i+6] |= (uint32_t)(a[13*i+11] & 0x07)<< 10;
        
        r->coeffs[8*i+7]  = a[13*i+11]>>3;
        r->coeffs[8*i+7] |= (uint32_t)a[13*i+12]<< 5;
        
        
        r->coeffs[8*i+0] = PARAM_Q + (1 << (PARAM_D-1)) - r->coeffs[8*i+0];
        r->coeffs[8*i+1] = PARAM_Q + (1 << (PARAM_D-1)) - r->coeffs[8*i+1];
        r->coeffs[8*i+2] = PARAM_Q + (1 << (PARAM_D-1)) - r->coeffs[8*i+2];
        r->coeffs[8*i+3] = PARAM_Q + (1 << (PARAM_D-1)) - r->coeffs[8*i+3];
        r->coeffs[8*i+4] = PARAM_Q + (1 << (PARAM_D-1)) - r->coeffs[8*i+4];
        r->coeffs[8*i+5] = PARAM_Q + (1 << (PARAM_D-1)) - r->coeffs[8*i+5];
        r->coeffs[8*i+6] = PARAM_Q + (1 << (PARAM_D-1)) - r->coeffs[8*i+6];
        r->coeffs[8*i+7] = PARAM_Q + (1 << (PARAM_D-1)) - r->coeffs[8*i+7];
        
    }
#elif PARAM_D==14
    for(i = 0; i < PARAM_N/4; ++i) {
        r->coeffs[4*i+0]  = a[7*i+0];
        r->coeffs[4*i+0] |= (uint32_t)(a[7*i+1] & 0x3F) << 8;
        
        r->coeffs[4*i+1]  = a[7*i+1] >> 6;
        r->coeffs[4*i+1] |= (uint32_t)a[7*i+2] << 2;
        r->coeffs[4*i+1] |= (uint32_t)(a[7*i+3] & 0x0F) << 10;
        
        r->coeffs[4*i+2]  = a[7*i+3] >> 4;
        r->coeffs[4*i+2] |= (uint32_t)a[7*i+4] << 4;
        r->coeffs[4*i+2] |= (uint32_t)(a[7*i+5] & 0x03) << 12;
        
        r->coeffs[4*i+3]  = a[7*i+5] >> 2;
        r->coeffs[4*i+3] |= (uint32_t)a[7*i+6] << 6;
        
        r->coeffs[4*i+0] = PARAM_Q + (1 << (PARAM_D-1)) - r->coeffs[4*i+0];
        r->coeffs[4*i+1] = PARAM_Q + (1 << (PARAM_D-1)) - r->coeffs[4*i+1];
        r->coeffs[4*i+2] = PARAM_Q + (1 << (PARAM_D-1)) - r->coeffs[4*i+2];
        r->coeffs[4*i+3] = PARAM_Q + (1 << (PARAM_D-1)) - r->coeffs[4*i+3];
    }
#endif
}

void aigis_polyzpack(unsigned char *r, const aigis_poly *a) {
#if GAMMA1 - BETA1 > (1 << 17)
#error "aigis_polyzpack() assumes GAMMA1 - BETA1 <= 2^{17}"
#endif
    unsigned int i;
    
    uint32_t t[4];
    for(i = 0; i < PARAM_N/4; ++i) {
        /* Map to {0,...,2*GAMMA1 - 2} */ // 18-bit
        t[0] = GAMMA1 - 1 - a->coeffs[4*i+0];
        t[0] += ((int32_t)t[0] >> 31) & PARAM_Q;
        t[1] = GAMMA1 - 1 - a->coeffs[4*i+1];
        t[1] += ((int32_t)t[1] >> 31) & PARAM_Q;
        t[2] = GAMMA1 - 1 - a->coeffs[4*i+2];
        t[2] += ((int32_t)t[2] >> 31) & PARAM_Q;
        t[3] = GAMMA1 - 1 - a->coeffs[4*i+3];
        t[3] += ((int32_t)t[3] >> 31) & PARAM_Q;
        
        r[9*i+0]  = t[0];
        r[9*i+1]  = t[0] >> 8;
        r[9*i+2]  = t[0] >> 16;
        r[9*i+2] |= t[1] << 2;
        r[9*i+3]  = t[1] >> 6;
        r[9*i+4]  = t[1] >> 14;
        r[9*i+4] |= t[2] << 4;
        r[9*i+5]  = t[2] >> 4;
        r[9*i+6]  = t[2] >> 12;
        r[9*i+6] |= t[3] << 6;
        r[9*i+7]  = t[3] >> 2;
        r[9*i+8]  = t[3] >> 10;
        
    }
}

void aigis_polyzunpack(aigis_poly *r, const unsigned char *a) {
    
#if GAMMA1 - BETA1 > (1 << 17)
#error "aigis_polyzunpack() assumes GAMMA1 - BETA1 <= 2^{17}"
#endif
    
    unsigned int i;
    for(i = 0; i < PARAM_N/4; ++i) {
        r->coeffs[4*i+0]  = a[9*i+0];
        r->coeffs[4*i+0] |= (uint32_t)a[9*i+1] << 8;
        r->coeffs[4*i+0] |= (uint32_t)(a[9*i+2] & 0x03) << 16;
        r->coeffs[4*i+0] = GAMMA1 - 1 - r->coeffs[4*i+0];
        r->coeffs[4*i+0] += ((int32_t)r->coeffs[4*i+0] >> 31) & PARAM_Q;
        
        r->coeffs[4*i+1]  = a[9*i+2] >> 2;
        r->coeffs[4*i+1] |= (uint32_t)a[9*i+3] << 6;
        r->coeffs[4*i+1] |= (uint32_t)(a[9*i+4] & 0x0F) << 14;
        r->coeffs[4*i+1] = GAMMA1 - 1 - r->coeffs[4*i+1];
        r->coeffs[4*i+1] += ((int32_t)r->coeffs[4*i+1] >> 31) & PARAM_Q;
        
        r->coeffs[4*i+2]  = a[9*i+4] >> 4;
        r->coeffs[4*i+2] |= (uint32_t)a[9*i+5] << 4;
        r->coeffs[4*i+2] |= (uint32_t)(a[9*i+6] & 0x3F) << 12;
        r->coeffs[4*i+2] = GAMMA1 - 1 - r->coeffs[4*i+2];
        r->coeffs[4*i+2] += ((int32_t)r->coeffs[4*i+2] >> 31) & PARAM_Q;
        
        
        r->coeffs[4*i+3]  = a[9*i+6] >> 6;
        r->coeffs[4*i+3] |= (uint32_t)a[9*i+7] << 2;
        r->coeffs[4*i+3] |= (uint32_t)a[9*i+8] << 10;
        r->coeffs[4*i+3] = GAMMA1 - 1 - r->coeffs[4*i+3];
        r->coeffs[4*i+3] += ((int32_t)r->coeffs[4*i+3] >> 31) & PARAM_Q;
    }
}

void aigis_polyw1_pack(unsigned char *r, const aigis_poly *a) 
{
#if PARAM_Q/ALPHA > 8
#error "aigis_polyw1_pack() assumes PARAM_Q/ALPHA -1 <= 7"
#endif
    unsigned int i;
    
    for(i = 0; i < PARAM_N/8; ++i)
    {
        r[3*i+0] = a->coeffs[8*i+0]      | (a->coeffs[8*i+1] << 3) | (a->coeffs[8*i+ 2] << 6);
        r[3*i+1] = (a->coeffs[8*i+2]>>2) | (a->coeffs[8*i+3] << 1) | (a->coeffs[8*i+ 4] << 4) | (a->coeffs[8*i+ 5] << 7);
        r[3*i+2] = (a->coeffs[8*i+5]>>1) | (a->coeffs[8*i+6] << 2) | (a->coeffs[8*i+ 7] << 5);
    }
}
