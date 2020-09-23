#include <stdint.h>
#include<stdio.h>
#include "api.h"
#include "params.h"
#include "sign.h"
#include "fips202.h"
#include "fips202x4.h"
#include "poly.h"
#include "polyvec.h"
#include "packing.h"
#include "randombytes.h"
#include "aes256ctr.h"


#ifdef USE_AES
void expand_mat(aigis_polyvecl mat[PARAM_K], const unsigned char rho[SEEDBYTES]) {
	unsigned int i, j;
#if PARAM_L == 3 && PARAM_K == 4
	unsigned int nblocks = 7;
	unsigned char buf[896];
#elif PARAM_L == 4 && PARAM_K == 5
	unsigned int nblocks = 8;
	unsigned char buf[1024];
#elif PARAM_L == 5 && PARAM_K == 6
	unsigned int nblocks = 9;
	unsigned char buf[1152];
#endif
	aigis_aes256ctr_ctx state;

	for (i = 0; i < PARAM_K; ++i)
		for (j = 0; j < PARAM_L; ++j)
		{
			aigis_aes256ctr_init(&state, rho, (i << 8) + j);
			aigis_aes256ctr_squeezeblocks(buf, nblocks, &state);
			aigis_poly_uniform(&mat[i].vec[j], buf);
		}
}
#elif PARAM_L == 3 && PARAM_K == 4
void expand_mat(aigis_polyvecl mat[4], const unsigned char rho[SEEDBYTES])
{
  unsigned int i;
  unsigned char inbuf[4][SEEDBYTES + 1];
  unsigned char outbuf[4][6*SHAKE128_RATE];

  for(i = 0; i < SEEDBYTES; ++i) {
    inbuf[0][i] = rho[i];
    inbuf[1][i] = rho[i];
    inbuf[2][i] = rho[i];
    inbuf[3][i] = rho[i];
  }

  for(i = 0; i < 3; ++i) {
    inbuf[0][SEEDBYTES] = 0 + (i << 4);
    inbuf[1][SEEDBYTES] = 1 + (i << 4);
    inbuf[2][SEEDBYTES] = 2 + (i << 4);
    inbuf[3][SEEDBYTES] = 3 + (i << 4);

    shake128_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], 6*SHAKE128_RATE,
                inbuf[0], inbuf[1], inbuf[2], inbuf[3], SEEDBYTES + 1);

    aigis_poly_uniform(&mat[0].vec[i], outbuf[0]);
    aigis_poly_uniform(&mat[1].vec[i], outbuf[1]);
    aigis_poly_uniform(&mat[2].vec[i], outbuf[2]);
    aigis_poly_uniform(&mat[3].vec[i], outbuf[3]);
  }
}

#elif PARAM_L == 4 && PARAM_K == 5

void expand_mat(aigis_polyvecl *mat, const unsigned char rho[SEEDBYTES])
{
  unsigned int i;
  unsigned char inbuf[4][SEEDBYTES + 1];
  unsigned char outbuf[4][7*SHAKE128_RATE];

  for(i = 0; i < SEEDBYTES; ++i) {
    inbuf[0][i] = rho[i];
    inbuf[1][i] = rho[i];
    inbuf[2][i] = rho[i];
    inbuf[3][i] = rho[i];
  }

  for(i = 0; i < 5; ++i) {
    inbuf[0][SEEDBYTES] = i + (0 << 4);
    inbuf[1][SEEDBYTES] = i + (1 << 4);
    inbuf[2][SEEDBYTES] = i + (2 << 4);
    inbuf[3][SEEDBYTES] = i + (3 << 4);

    shake128_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], 7*SHAKE128_RATE, inbuf[0], inbuf[1], inbuf[2], inbuf[3], SEEDBYTES + 1);

    aigis_poly_uniform(&mat[i].vec[0], outbuf[0]);
    aigis_poly_uniform(&mat[i].vec[1], outbuf[1]);
    aigis_poly_uniform(&mat[i].vec[2], outbuf[2]);
    aigis_poly_uniform(&mat[i].vec[3], outbuf[3]);
  }
}
#elif PARAM_L == 5 && PARAM_K == 6

void expand_mat(aigis_polyvecl mat[6], const unsigned char rho[SEEDBYTES])
{
	unsigned int i;
	unsigned char inbuf[4][SEEDBYTES + 1];
	unsigned char outbuf[4][8 * SHAKE128_RATE];

	for (i = 0; i < SEEDBYTES; ++i) {
		inbuf[0][i] = rho[i];
		inbuf[1][i] = rho[i];
		inbuf[2][i] = rho[i];
		inbuf[3][i] = rho[i];
	}

	for (i = 0; i < 4; i += 2) {
		inbuf[0][SEEDBYTES] = 0 + (i << 4);
		inbuf[1][SEEDBYTES] = 1 + (i << 4);
		inbuf[2][SEEDBYTES] = 2 + (i << 4);
		inbuf[3][SEEDBYTES] = 3 + (i << 4);

		shake128_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], 8 * SHAKE128_RATE,
			inbuf[0], inbuf[1], inbuf[2], inbuf[3], SEEDBYTES + 1);

		aigis_poly_uniform(&mat[0].vec[i], outbuf[0]);
		aigis_poly_uniform(&mat[1].vec[i], outbuf[1]);
		aigis_poly_uniform(&mat[2].vec[i], outbuf[2]);
		aigis_poly_uniform(&mat[3].vec[i], outbuf[3]);

		inbuf[0][SEEDBYTES] = 4 + (i << 4);
		inbuf[1][SEEDBYTES] = 5 + (i << 4);
		inbuf[2][SEEDBYTES] = 0 + ((i + 1) << 4);
		inbuf[3][SEEDBYTES] = 1 + ((i + 1) << 4);

		shake128_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], 8 * SHAKE128_RATE,
			inbuf[0], inbuf[1], inbuf[2], inbuf[3], SEEDBYTES + 1);

		aigis_poly_uniform(&mat[4].vec[i], outbuf[0]);
		aigis_poly_uniform(&mat[5].vec[i], outbuf[1]);
		aigis_poly_uniform(&mat[0].vec[i + 1], outbuf[2]);
		aigis_poly_uniform(&mat[1].vec[i + 1], outbuf[3]);

		inbuf[0][SEEDBYTES] = 2 + ((i + 1) << 4);
		inbuf[1][SEEDBYTES] = 3 + ((i + 1) << 4);
		inbuf[2][SEEDBYTES] = 4 + ((i + 1) << 4);
		inbuf[3][SEEDBYTES] = 5 + ((i + 1) << 4);

		shake128_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], 8 * SHAKE128_RATE,
			inbuf[0], inbuf[1], inbuf[2], inbuf[3], SEEDBYTES + 1);

		aigis_poly_uniform(&mat[2].vec[i + 1], outbuf[0]);
		aigis_poly_uniform(&mat[3].vec[i + 1], outbuf[1]);
		aigis_poly_uniform(&mat[4].vec[i + 1], outbuf[2]);
		aigis_poly_uniform(&mat[5].vec[i + 1], outbuf[3]);
	}

	inbuf[0][SEEDBYTES] = 0 + (4 << 4);
	inbuf[1][SEEDBYTES] = 1 + (4 << 4);
	inbuf[2][SEEDBYTES] = 2 + (4 << 4);
	inbuf[3][SEEDBYTES] = 3 + (4 << 4);

	shake128_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], 8 * SHAKE128_RATE,
		inbuf[0], inbuf[1], inbuf[2], inbuf[3], SEEDBYTES + 1);

	aigis_poly_uniform(&mat[0].vec[4], outbuf[0]);
	aigis_poly_uniform(&mat[1].vec[4], outbuf[1]);
	aigis_poly_uniform(&mat[2].vec[4], outbuf[2]);
	aigis_poly_uniform(&mat[3].vec[4], outbuf[3]);

	inbuf[0][SEEDBYTES] = 4 + (4 << 4);
	inbuf[1][SEEDBYTES] = 5 + (4 << 4);

	shake128_4x(outbuf[0], outbuf[1], outbuf[2], outbuf[3], 8 * SHAKE128_RATE,
		inbuf[0], inbuf[1], inbuf[2], inbuf[3], SEEDBYTES + 1);

	aigis_poly_uniform(&mat[4].vec[4], outbuf[0]);
	aigis_poly_uniform(&mat[5].vec[4], outbuf[1]);
}
#endif

void challenge(aigis_poly *c,
               const unsigned char mu[CRHBYTES],
               const aigis_polyveck *w1) 
{
  unsigned int i, b, pos;
  unsigned char inbuf[CRHBYTES + PARAM_K*POLW1_SIZE_PACKED];
  unsigned char outbuf[SHAKE256_RATE];
  uint64_t state[25], signs, mask;

  for(i = 0; i < CRHBYTES; ++i)
    inbuf[i] = mu[i];
  for(i = 0; i < PARAM_K; ++i)
    aigis_polyw1_pack(inbuf + CRHBYTES + i*POLW1_SIZE_PACKED, w1->vec+i);

  shake256_absorb_aigis_mulan(state, inbuf, sizeof(inbuf)); //fzhang
  shake256_squeezeblocks_aigis_mulan(outbuf, 1, state);  //fzhang

  signs = 0;
  for(i = 0; i < 8; ++i)
    signs |= (uint64_t)outbuf[i] << 8*i;

  pos = 8;
  mask = 1;

  for(i = 0; i < PARAM_N; ++i)
    c->coeffs[i] = 0;

  for(i = 196; i < 256; ++i) {
    do {
      if(pos >= SHAKE256_RATE) {
        shake256_squeezeblocks_aigis_mulan(outbuf, 1, state); //fzhang
        pos = 0;
      }

      b = outbuf[pos++];
    } while(b > i);

    c->coeffs[i] = c->coeffs[b];
    c->coeffs[b] = (signs & mask) ? PARAM_Q - 1 : 1;
    mask <<= 1;
  }
}

/*************************************************
* generate a pair of public key pk and secret key sk, 
* where pk = rho|t1
*       sk = rho|key|hash(pk)|s1|s2|t0
**************************************************/
int msig_keygen(unsigned char *pk, unsigned char *sk) 
{
  unsigned int i;
  unsigned char buf[3*SEEDBYTES + CRHBYTES]; //buf = r|rho|key|hash(pk)
  uint16_t nonce = 0;
  aigis_polyvecl mat[PARAM_K];
  aigis_polyvecl s1, s1hat;
  aigis_polyveck s2, t, t1, t0;
  
  aigis_randombytes(buf, SEEDBYTES);
  shake256(buf, 3*SEEDBYTES, buf, SEEDBYTES); 
  

  expand_mat(mat, &buf[SEEDBYTES]);

#ifdef USE_AES
  for (i = 0; i < PARAM_L; ++i)
	  aigis_poly_uniform_eta1(&s1.vec[i], buf, nonce++);
#elif PARAM_L==3
  aigis_poly_uniform_eta1_3x(&s1.vec[0],&s1.vec[1],&s1.vec[2],buf, nonce,nonce+1,nonce+2);
  nonce += 3;
#elif PARAM_L==4
  aigis_poly_uniform_eta1_4x(&s1.vec[0],
                         &s1.vec[1],
                         &s1.vec[2],
                         &s1.vec[3],buf,nonce,nonce+1,nonce+2,nonce+3);
  nonce += 4;
#elif PARAM_L==5
  aigis_poly_uniform_eta1_4x(&s1.vec[0],
	  &s1.vec[1],
	  &s1.vec[2],
	  &s1.vec[3],buf,nonce,nonce+1,nonce+2,nonce+3);
  nonce += 4;
  aigis_poly_uniform_eta1(&s1.vec[4], buf, nonce++);
#endif
    
#ifdef USE_AES
  for (i = 0; i < PARAM_K; ++i)
	  aigis_poly_uniform_eta2(&s2.vec[i], buf, nonce++);
#elif PARAM_K == 4
  aigis_poly_uniform_eta2_4x(&s2.vec[0],
                         &s2.vec[1],
                         &s2.vec[2],
                         &s2.vec[3],buf,nonce,nonce+1,nonce+2,nonce+3);
#elif PARAM_K == 5
  aigis_poly_uniform_eta2_4x(&s2.vec[0],
                         &s2.vec[1],
                         &s2.vec[2],
                         &s2.vec[3],buf,nonce,nonce+1,nonce+2,nonce+3);
  aigis_poly_uniform_eta2(&s2.vec[4], buf, nonce+4);
#elif PARAM_K == 6
  aigis_poly_uniform_eta2_4x(&s2.vec[0],
	  &s2.vec[1],
	  &s2.vec[2],
	  &s2.vec[3], buf, nonce, nonce + 1, nonce + 2, nonce + 3);
  aigis_poly_uniform_eta2(&s2.vec[4], buf, nonce + 4);
  aigis_poly_uniform_eta2(&s2.vec[5], buf, nonce + 5);
#endif

  s1hat = s1;
  aigis_polyvecl_ntt(&s1hat);
  for(i = 0; i < PARAM_K; ++i) {

    aigis_polyvecl_pointwise_acc_invmontgomery(&t.vec[i], mat+i, &s1hat);//output coefficient  <=2*Q
    aigis_poly_invntt_montgomery(t.vec+i);//output coefficient  <=2*Q
  }


  aigis_polyveck_add(&t, &t, &s2);//output coefficient <=4*Q
  aigis_polyveck_freeze4q(&t);
  aigis_polyveck_power2round(&t1, &t0, &t);
  aigis_pack_pk(pk, &buf[SEEDBYTES], &t1);
  

  shake256(&buf[3*SEEDBYTES], CRHBYTES, pk, AIGIS_SIG_PUBLICKEYBYTES);
  aigis_pack_sk(sk,&buf[SEEDBYTES], &s1, &s2, &t0);

  return 0;

}

/*************************************************
* create a signature sm on message m, where 
* sm = z|h|c
**************************************************/
int msig_sign(unsigned char *sk,
	unsigned char *m, unsigned long long mlen,
	unsigned char *sm, unsigned long long *smlen)
{
	unsigned long long i;
	unsigned int n;
	unsigned char *buf;
	uint16_t nonce = 0;
	aigis_poly     c, chat;
	aigis_polyvecl mat[PARAM_K], s1, y, yhat, z;
	aigis_polyveck s2, t0, w, w1;
	aigis_polyveck h, wcs2, ct0, tmp;

	buf = (unsigned char *)malloc(2 * SEEDBYTES + CRHBYTES + mlen);
	if (buf == NULL)
		return -1;
	aigis_unpack_sk(buf, &s1, &s2, &t0, sk);

	for (i = 0; i<mlen; i++)
		buf[2 * SEEDBYTES + CRHBYTES + i] = m[i];
	shake256(&buf[2 * SEEDBYTES], CRHBYTES, &buf[2 * SEEDBYTES], CRHBYTES + mlen);

	expand_mat(mat, buf);
	aigis_polyvecl_ntt(&s1);
	aigis_polyveck_ntt(&s2);
	aigis_polyveck_ntt(&t0);

#ifdef USE_AES
	shake256(&buf[SEEDBYTES], SEEDBYTES, &buf[SEEDBYTES], SEEDBYTES + CRHBYTES);
#endif

rej:
#ifdef USE_AES
	for (i = 0; i < PARAM_L; ++i)
		aigis_poly_uniform_gamma1m1(&y.vec[i], &buf[SEEDBYTES], nonce++);
#elif PARAM_L == 3
	aigis_poly_uniform_gamma1m1_3x(&y.vec[0], &y.vec[1], &y.vec[2], &buf[SEEDBYTES],
		nonce, nonce + 1, nonce + 2);
	nonce += 3;
#elif PARAM_L == 4
	aigis_poly_uniform_gamma1m1_4x(&y.vec[0], &y.vec[1], &y.vec[2], &y.vec[3], &buf[SEEDBYTES],
		nonce, nonce + 1, nonce + 2, nonce + 3);
	nonce += 4;
#elif PARAM_L == 5
	aigis_poly_uniform_gamma1m1_4x(&y.vec[0], &y.vec[1], &y.vec[2], &y.vec[3], &buf[SEEDBYTES],
		nonce, nonce + 1, nonce + 2, nonce + 3);
	aigis_poly_uniform_gamma1m1(&y.vec[4], &buf[SEEDBYTES], nonce+4);
	nonce += 5;
#endif 

	yhat = y;
	aigis_polyvecl_ntt(&yhat);
	for (i = 0; i < PARAM_K; ++i) {
		aigis_polyvecl_pointwise_acc_invmontgomery(w.vec + i, mat + i, &yhat);//output coefficient  <=2*Q
		aigis_poly_invntt_montgomery(w.vec + i);//output coefficient  <= 2*Q
	}

	aigis_polyveck_freeze2q(&w);
	aigis_polyveck_decompose(&w1, &tmp, &w);
	challenge(&c, &buf[2 * SEEDBYTES], &w1);

	chat = c;
	aigis_poly_ntt(&chat);

	

	//check z
	for (i = 0; i < PARAM_L; ++i) {
		aigis_poly_pointwise_invmontgomery(z.vec + i, &chat, s1.vec + i);//output coefficient  <=2*Q
		aigis_poly_invntt_montgomery(z.vec + i);//output coefficient  <=2*Q
	}
	aigis_polyvecl_add(&z, &z, &y);//output coefficient  <=4*Q
	aigis_polyvecl_freeze4q(&z);
	if (aigis_polyvecl_chknorm(&z, GAMMA1 - BETA1))
		goto rej;

	//check r0 and r1=w1?
	for (i = 0; i < PARAM_K; ++i) {
		aigis_poly_pointwise_invmontgomery(wcs2.vec + i, &chat, s2.vec + i);//output coefficient  <=2*Q
		aigis_poly_invntt_montgomery(wcs2.vec + i);//output coefficient  <=2*Q
	}
	aigis_polyveck_sub(&tmp, &tmp, &wcs2);//output coefficient  <=4*Q
	aigis_polyveck_freeze4q(&tmp);
	if (aigis_polyveck_chknorm(&tmp, GAMMA2 - BETA2))
		goto rej;

	//check ct0
	for (i = 0; i < PARAM_K; ++i) {
		aigis_poly_pointwise_invmontgomery(ct0.vec + i, &chat, t0.vec + i);//output coefficient  <=2*Q
		aigis_poly_invntt_montgomery(ct0.vec + i);//output coefficient  <=2*Q
	}
	aigis_polyveck_freeze2q(&ct0);
	if (aigis_polyveck_chknorm(&ct0, GAMMA2))
		goto rej;

	//check hw(h)
	aigis_polyveck_add(&tmp, &tmp, &ct0);//output coefficient  <=2*Q
	aigis_polyveck_freeze2q(&tmp);
	n = aigis_polyveck_make_hint(&h, &tmp, &w1);
	if (n > OMEGA)
		goto rej;

	aigis_pack_sig(sm, &z, &h, &c);

	*smlen = AIGIS_SIG_BYTES;

	free(buf);
	return 0;
}

int msig_verf(unsigned char *pk,
	unsigned char *sm, unsigned long long smlen,
	unsigned char *m, unsigned long long mlen)
{
	unsigned long long i;
	unsigned char rho[SEEDBYTES];
	unsigned char *buf;
	aigis_poly     c, chat, cp;
	aigis_polyvecl mat[PARAM_K], z;
	aigis_polyveck t1, w1, h, tmp1, tmp2;

	if (smlen < AIGIS_SIG_BYTES)
		return 0;

	aigis_unpack_pk(rho, &t1, pk);
	aigis_unpack_sig(&z, &h, &c, sm);
	if (aigis_polyvecl_chknorm(&z, GAMMA1 - BETA1))
		return 0;

	buf = (unsigned char*)malloc(CRHBYTES + mlen);
	if (buf == NULL)
		return -1;
	for (i = 0; i<mlen; i++)
		buf[CRHBYTES + i] = m[i];

	shake256(buf, CRHBYTES, pk, AIGIS_SIG_PUBLICKEYBYTES);
	shake256(buf, CRHBYTES, buf, CRHBYTES + mlen);

	expand_mat(mat, rho);

	aigis_polyvecl_ntt(&z);
	for (i = 0; i < PARAM_K; ++i)
		aigis_polyvecl_pointwise_acc_invmontgomery(tmp1.vec + i, mat + i, &z);//output coefficient  <=2*Q

	chat = c;
	aigis_poly_ntt(&chat);
	aigis_polyveck_shiftl(&t1, PARAM_D);
	aigis_polyveck_ntt(&t1);
	for (i = 0; i < PARAM_K; ++i)
		aigis_poly_pointwise_invmontgomery(tmp2.vec + i, &chat, t1.vec + i);//output coefficient  <=2*Q

	aigis_polyveck_sub(&tmp1, &tmp1, &tmp2);//output coefficient  <=4*Q
	aigis_polyveck_invntt_montgomery(&tmp1);//output coefficient  <= 2*Q

	aigis_polyveck_freeze2q(&tmp1);
	aigis_polyveck_use_hint(&w1, &tmp1, &h);

	challenge(&cp, buf, &w1);
	for (i = 0; i < PARAM_N; ++i)
		if (c.coeffs[i] != cp.coeffs[i])
			return 0;

	free(buf);
	return 1;
}