#include <stdint.h>
#include "Parameters.h"
#include "io.h"
#include "Alg.h"
#include "Polynomial.h"

int mulan_sig_keygen(
		unsigned char * pk, unsigned long long * pk_bytes, 
		unsigned char * sk, unsigned long long * sk_bytes) 
{
	uint32_t i, j;
	unsigned char seedbuf[3*SEEDBYTES];
	unsigned char * rho_se = seedbuf;			// for (s,e)
	unsigned char * rho_A = seedbuf + SEEDBYTES;	// for the matrix A
	unsigned char * key = seedbuf + 2*SEEDBYTES;	// for Key. 

	Mulan_Polynomial A[ROWS_A][COLUMNS_A];
	Mulan_Polynomial t[DIM_t], t1[DIM_t], t0[DIM_t]; 
	Mulan_Polynomial s[DIM_s];
	Mulan_Polynomial e[DIM_e]; 
	Mulan_Polynomial temp;
	 
	// sk1 = (rho_A, key, s, e)
	// sk2 = (tr, t0)
	unsigned char * sk1 = sk + 0;
	unsigned char * sk2 = sk + SIZE_sk1;
	unsigned char * tr = sk2 + 0;

	// check invalid argument
	if(NULL == pk || NULL == sk)
		return -2;

	// step 1: generate random bytes. 
	mulan_randombytes(seedbuf,SEEDBYTES);
	shake256(seedbuf, 3*SEEDBYTES, seedbuf, SEEDBYTES);
	
	// step 2: generate matrix A. 
	Mulan_Generate_A(A, rho_A);

	// step 3: generate s, e. 
	for(i = 0; i < DIM_s; i++)
		Mulan_Poly_uniform_eta(s+i, rho_se, i);	
	for(i = 0; i < DIM_e; i++)
		Mulan_Poly_uniform_eta(e+i, rho_se, DIM_s+i);
	
	// step 4: pack sk1 = (seed_A, key, s, e)
	Mulan_pack_sk1(sk1, rho_A, key, s, e);

	//step 5: compute t = As+e.  
 	Mulan_PolyVec_forward_ntt(s, DIM_s);
	for(i=0; i<ROWS_K; i++)
	{	
		Mulan_Poly_pointwise_invmontgomery(t+i, A[i]+0, s+0);
		for(j=1; j<COLUMNS_A; j++)
		{
			Mulan_Poly_pointwise_invmontgomery(&temp, A[i]+j, s+j);
			Mulan_Poly_add(t+i, t+i, &temp);
		}
		Mulan_Poly_invntt_frominvmont(t+i);
	}
	Mulan_PolyVec_add(t, t, e);
	
	// step 6: compute (t1, t0) <- (t1, t0) := Mulan_Power2Round(t). 
	Mulan_PolyVec_Power2Round(t1, t0, t);
	
	// step 7: pack pk = (rho_A, t1).
	Mulan_pack_pk(pk, rho_A, t1);

	// step 8: tr = CRH(pk) = CRH(rho_A, t1)
	shake256(tr, CRHBYTES, pk, SIZE_pk);
	
	// step 9: pack sk2 = (tr, t0);
 	Mulan_pack_sk2(sk2, tr, t0);
	
	// step 10: return
	*pk_bytes = MULAN_SIG_PUBLICKEYBYTES;
	*sk_bytes = MULAN_SIG_SECRETKEYBYTES;
	return 0;
}

int mulan_sig_sign(
	unsigned char * sk, unsigned long long sk_bytes,
	unsigned char * m, unsigned long long mlen,
	unsigned char * sm, unsigned long long * smlen)
{
	uint32_t i, j, counter, flag;
	uint32_t nonce = 0;

	Mulan_Polynomial A[ROWS_A][COLUMNS_A];
	Mulan_Polynomial s[DIM_s], e[DIM_e];
	Mulan_Polynomial t0[DIM_t0];
	
	Mulan_Polynomial u[ROWS_K], v[ROWS_K];
	
	Mulan_Polynomial c, c_backup; 
	Mulan_Polynomial y[DIM_y];
	Mulan_Polynomial z[DIM_z];
	Mulan_Polynomial w[DIM_w], w1[DIM_w1], w0[DIM_w1];
	
	Mulan_Polynomial r0[ROWS_K], r1[ROWS_K];
	Mulan_Polynomial h[ROWS_K];

	unsigned char * sk1 = sk;
	unsigned char * sk2 = sk + SIZE_sk1;
	
	unsigned char buffer[CRHBYTES+MULAN_MESSAGE_LENGTH_MAX];
	unsigned char * tr = buffer; 
	
	// seebuf = (seed_A, seed_key, mu) 
	unsigned char seedbuf[2*SEEDBYTES + CRHBYTES];
	unsigned char * rho_A = seedbuf + 0;
	unsigned char * key = seedbuf + SEEDBYTES; 
	unsigned char * mu = seedbuf + 2*SEEDBYTES;

	// add by lmwen to check invalid argument
	if(NULL ==sk ||NULL == m || NULL == sm || NULL == smlen || sk_bytes!=MULAN_SIG_SECRETKEYBYTES)
		return -1;
	
	if(mlen>MULAN_MESSAGE_LENGTH_MAX || mlen==0)
	{
		// printf("For ease of demonstration, the signing algorithm only handles message of length:  0< length<=%d. \n", MULAN_MESSAGE_LENGTH_MAX);
		return -2;
	}
		
	// step1: unpack sk1, sk2, and generate A. 
	Mulan_unpack_sk1(A, rho_A, key, s, e, sk1);
	Mulan_unpack_sk2(tr, t0, sk2);
	
	// step 2: mu = CRH(tr, message); 
	memcpy(buffer+CRHBYTES, m, mlen);
	shake256(mu, CRHBYTES, buffer, CRHBYTES + mlen);
	
	
	Mulan_PolyVec_forward_ntt(s, DIM_s);
	Mulan_PolyVec_forward_ntt(e, DIM_e);
	Mulan_PolyVec_forward_ntt(t0, DIM_t0);
	
	while(1)
	{
		// step 3: sample y
		for(i=0; i < DIM_y; i++)
		{
			Mulan_poly_uniform_gamma1m1(y+i, key, nonce++);
			Mulan_Poly_assignment(z+i, y+i);
		}

		// step 4: compute w = Ay
		Mulan_PolyVec_forward_ntt(y, DIM_y);
		for(i=0; i<ROWS_K; i++)
		{	
			Mulan_Poly_pointwise_invmontgomery(w+i, A[i]+0, y+0);
			for(j=1; j<COLUMNS_A; j++)
			{
				Mulan_Poly_pointwise_invmontgomery(r0+i, A[i]+j, y+j);
				Mulan_Poly_add(w+i, w+i, r0+i);
			}
			Mulan_Poly_invntt_frominvmont(w+i);
		}
		// Mulan_PolyVec_freeze(w, DIM_w);
				
		// step 5: compute w1 
		Mulan_PolyVec_decompose(w1, w0, w, DIM_w); 

		// step 6: sample c
		Mulan_challenge(&c, mu, w1);
		Mulan_Poly_assignment(&c_backup, &c);
		
		// step 7: z := y + c*s.
		// step 8: flag1 = Norm{z} vs. bound.
		Mulan_PolyVec_forward_ntt(&c, DIM_c);
		for(j=0; j<DIM_s; j++)
		{
			Mulan_Poly_pointwise_invmontgomery(v+j, &c, s+j);
			Mulan_Poly_invntt_frominvmont(v+j);
			Mulan_Poly_add(z+j, z+j, v+j);
			flag = Mulan_Poly_check_norm(z+j, NORMBOUND_z);
			if(flag) break;
		}
		if(j<DIM_s) continue;	
		

		for(j=0; j<ROWS_K; j++)
		{
			// step 9: u := w-ce, 
			Mulan_Poly_pointwise_invmontgomery(u+j, &c, e+j);
			Mulan_Poly_invntt_frominvmont(u+j);
			Mulan_Poly_sub(u+j, w+j, u+j);

			// step 10: [r1, r0] := Mulan_Decompose(u[j]), where u:=w-ce. 
			Mulan_Poly_decompose(r1+j, r0+j, u+j);

			// step 11: flag2 = (Norm{r0} >= bound) 
			flag = Mulan_Poly_check_norm(r0+j, NORMBOUND_r0);
			if(flag) break;

			// step 12: flag3 = (r1<>w1)
			flag = Mulan_Poly_compare(r1+j, w1+j);
			if(flag) break;
		}
		if(j<ROWS_K) continue;
	
		for(j=0; j<ROWS_K; j++)
		{
			//step 13: v := ct0. 
			Mulan_Poly_pointwise_invmontgomery(v+j, &c, t0+j);
			Mulan_Poly_invntt_frominvmont(v+j);
			
			//step 14: flag4 := (Norm(v[i]) >= bound)
			flag = Mulan_Poly_check_norm(v+j, NORMBOUND_ct0);
			if(flag)  break; 
		}
		if(j<ROWS_K) continue;

		// step 15: flag5 = #1's in h vs. OMEGA.
		Mulan_PolyVec_add(v, u, v);
		counter = Mulan_PolyVec_make_hint(h, v, u);
		if(counter>OMEGA) continue; 
		
		break;
	}	

	// step 16: pack signature=(z,h,c). 
	Mulan_PolyVec_freeze(z, DIM_z);
	Mulan_pack_signature(sm, z, h, &c_backup); 
	*smlen = MULAN_SIG_BYTES;

	// step 17: return. 
	return 0;
}

int mulan_sig_verf(
	unsigned char * pk, unsigned long long pk_bytes,
	unsigned char * sm, unsigned long long smlen,
	unsigned char * m, unsigned long long mlen)
{
	Mulan_Polynomial A[ROWS_A][COLUMNS_A], t1[DIM_t1];
	Mulan_Polynomial z[DIM_z], h[DIM_h], c; 
	Mulan_Polynomial c_backup; 
	
	Mulan_Polynomial u[ROWS_K];
	Mulan_Polynomial v[ROWS_K];
	Mulan_Polynomial temp;

	unsigned char rho_A[SEEDBYTES]; 
	unsigned char mu[CRHBYTES];
	unsigned char buffer[CRHBYTES+MULAN_MESSAGE_LENGTH_MAX];

	uint32_t flag; 
	uint32_t i, j; 

	// check invalid argument
	if(NULL ==pk ||NULL == m || NULL == sm || smlen != MULAN_SIG_BYTES || pk_bytes!=MULAN_SIG_PUBLICKEYBYTES)
		return -1;
	
	if(mlen>MULAN_MESSAGE_LENGTH_MAX || mlen==0)
	{
		printf("For ease of demonstration, the signing algorithm only handles message of length:  0< length<=%d. \n", MULAN_MESSAGE_LENGTH_MAX);
		return -2;
	}

	// step 1: unpack pk, and generate A. 
	Mulan_unpack_pk(A, rho_A, t1, pk);

	// step 2: unpack (z, h, c) from sm.
	Mulan_unpack_signature(z, h, &c, sm);

	// step 3: flag1: norm(z) vs. NORMBOUND_z.
	flag = Mulan_PolyVec_check_norm(z, DIM_z, NORMBOUND_z);
	if(flag) // indicatin failure. 
		return 1;
	
	// step 4: mu := CRH(tr, message)
	shake256(buffer, CRHBYTES, pk, SIZE_pk);
	memcpy(buffer+CRHBYTES, m, mlen);
	shake256(mu, CRHBYTES,buffer, CRHBYTES + mlen);
	
 	// step 5: u := Az-ct1, v := ct1. 
	Mulan_Poly_assignment(&c_backup, &c);
	Mulan_PolyVec_forward_ntt(&c, DIM_c);
	Mulan_PolyVec_forward_ntt(z, DIM_z);
	for(i=0; i<ROWS_A; i++)
	{
		Mulan_Poly_pointwise_invmontgomery(u+i, A[i]+0, z+0);
		for(j=1; j<COLUMNS_A; j++)
		{
			Mulan_Poly_pointwise_invmontgomery(&temp, A[i]+j, z+j);
			Mulan_Poly_add(u+i, u+i, &temp);
		}

		Mulan_Poly_forward_ntt(t1+i);	
		Mulan_Poly_pointwise_invmontgomery(v+i, &c, t1+i);
		Mulan_Poly_sub(u+i, u+i, v+i);
		Mulan_Poly_invntt_frominvmont(u+i);
	}

	// step 6: w1' := UseHint(h, Az-ct1).
	Mulan_PolyVec_use_hint(v, h, u);

	// step 7: temp := H(mu, w1')
	Mulan_challenge(&temp, mu, v);

	// step 8: flag3 = Mulan_Poly_compare(&temp, &backup_c);
	flag = Mulan_Poly_compare(&temp, &c_backup);
	
	return (flag==0);
}

