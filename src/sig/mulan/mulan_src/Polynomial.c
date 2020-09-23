#include "Polynomial.h"
#include "fips202.h"
const uint32_t uniform_eta_table[] = {
		Q-2,Q-1,Q,Q+1,Q+2,
		Q-2,Q-1,Q,Q+1,Q+2,
		Q-2,Q-1,Q,Q+1,Q+2};
		 
uint32_t Mulan_montgomery_reduce(uint64_t a) {
	uint64_t t;

	t = a * (uint32_t)(2226768383);
	t &= (1ULL << 32) - 1;
	t *= Q;
	t = a + t;
	t >>= 32;
	return (uint32_t)t; // warning is generated here. 
}


void Mulan_Poly_forward_ntt(Mulan_Polynomial *a)
{
	uint32_t *p = a->coefficients;
	unsigned int len, start, j, k;
	uint32_t zeta, t;

	k = 1;
	for (len = 128; len > 0; len >>= 1) {
		for (start = 0; start < N; start = j + len) {
			zeta = mulan_zetas[k++];
			for (j = start; j < start + len; ++j) {
				t = Mulan_montgomery_reduce((uint64_t)zeta * p[j + len]);
				p[j + len] = p[j] + 2 * Q - t;
				p[j] = p[j] + t;
			}
		}
	}
}		
void Mulan_PolyVec_forward_ntt(Mulan_Polynomial * p, const uint32_t dim)
{
	for(uint32_t i=0; i<dim; i++)
		Mulan_Poly_forward_ntt(p+i);
}


// compute the uniform Mulan_Polynomial Y based on "seed+nonce".
void Mulan_poly_uniform_gamma1m1(
		Mulan_Polynomial * a,
		const unsigned char seed[SEEDBYTES + CRHBYTES],
		uint16_t nonce)
{
	unsigned int ctr;
	unsigned char inbuf[SEEDBYTES + CRHBYTES + 2];
	unsigned char outbuf[5*SHAKE256_RATE];
	uint64_t state[25];

	memcpy(inbuf, seed, SEEDBYTES + CRHBYTES);

	inbuf[SEEDBYTES + CRHBYTES] = (nonce & 0xFF);
	inbuf[SEEDBYTES + CRHBYTES + 1] = (nonce >> 8);

	shake256_absorb_aigis_mulan(state, inbuf, SEEDBYTES + CRHBYTES + 2); //fzhang
	shake256_squeezeblocks_aigis_mulan(outbuf, 5, state);  				//fzhang

	ctr = Mulan_rej_gamma1m1(a->coefficients, N, outbuf, 5*SHAKE256_RATE); 

	if(ctr < N)
	{
		shake256_squeezeblocks_aigis_mulan(outbuf, 1, state);  //fzhang
		Mulan_rej_gamma1m1(a->coefficients + ctr, N - ctr, outbuf, SHAKE256_RATE);
	}
}


unsigned int Mulan_rej_gamma1m1(
		uint32_t * a,
		unsigned int len,
		const unsigned char *buf,
		unsigned int buflen)
{
	unsigned int ctr, pos;
	uint32_t t0, t1;

	// [0, 2r-2], [0, 4r-3]. 
	const uint32_t beta = 2*GAMMA1-1;
	const uint32_t beta2 = 4*GAMMA1-3;
	ctr = pos = 0;
	while(ctr < len) 
	{
		t0  = buf[pos];
		t0 |= (uint32_t)buf[pos + 1] << 8;
		t0 |= (uint32_t)buf[pos + 2] << 16;
		t0 &= 0xFFFFF;

		t1  = buf[pos + 2] >> 4;
		t1 |= (uint32_t)buf[pos + 3] << 4;
		t1 |= (uint32_t)buf[pos + 4] << 12;
		
		pos += 5;

		if(t0 <= beta2)
		{
			t0 -= (beta & (((int32_t)(beta-1-t0))>>31));	
			a[ctr++] = Q + (GAMMA1-1) - t0;
		}

		if(t1 <= beta2 && ctr < len)
		{
			t1 -= (beta & (((int32_t)(beta-1-t1))>>31));
			a[ctr++] = Q + (GAMMA1-1) - t1;
		}

		if(pos > buflen - 5)
		  break;
	}

	return ctr;
}


// [t1, t0] := Mulan_Power2Round(t)
void Mulan_PolyVec_Power2Round(
		Mulan_Polynomial * t1, 
		Mulan_Polynomial * t0, 
		const Mulan_Polynomial * t)
{	
	uint32_t i, j;
	for(i=0; i<ROWS_A; i++)
		for(j=0; j<N; j++)
			Mulan_Power2Round(
					t1[i].coefficients+j, 
					t0[i].coefficients+j, 
					t[i].coefficients[j]);
}



void Mulan_Poly_add(Mulan_Polynomial *c, const Mulan_Polynomial *a, const Mulan_Polynomial *b)
{
	for(uint32_t k=0; k<N; k++)
		c->coefficients[k] = a->coefficients[k] + b->coefficients[k];
}

// c[] := a[] + b[];
void Mulan_PolyVec_add(
		Mulan_Polynomial * c, 
		const Mulan_Polynomial * a, const Mulan_Polynomial * b)
{
	uint32_t i;
	
	for(i=0; i<ROWS_K; i++)
		Mulan_Poly_add(c+i, a+i, b+i);

}

uint32_t Mulan_Poly_compare(const Mulan_Polynomial *a, const Mulan_Polynomial *b)
{
	for(uint32_t k=0; k<N; k++)
		if(a->coefficients[k]!=b->coefficients[k])
			return 1;
	return 0;
}


void Mulan_Poly_assignment(Mulan_Polynomial *b, const Mulan_Polynomial *a)
{
	for(uint32_t k=0; k<N; k++)
		b->coefficients[k] = a->coefficients[k];
}

void Mulan_Poly_sub(
		Mulan_Polynomial *c, 
		const Mulan_Polynomial *a, 
		const Mulan_Polynomial *b)
{
	for(uint32_t k =0; k<N; k++)
		c->coefficients[k] = (a->coefficients[k]+(Q<<1)-b->coefficients[k]);
}

// compute the small poly based on "seed+nonce"
void Mulan_Poly_uniform_eta(
		Mulan_Polynomial * poly,
		const unsigned char * seed, 
		const unsigned char nonce)
{
  	unsigned int ctr;
  	unsigned char inbuf[SEEDBYTES + 1];

  	unsigned char outbuf[2*SHAKE256_RATE];
  	uint64_t state[25]; 
	
	memcpy(inbuf, seed, SEEDBYTES);
  	inbuf[SEEDBYTES] = nonce;
	
  	shake256_absorb_aigis_mulan(state, inbuf, SEEDBYTES + 1); //fzhang
  	shake256_squeezeblocks_aigis_mulan(outbuf, 2, state);	  //fzhang
	
  	ctr = Mulan_Rej_eta(poly->coefficients, N, outbuf, 2*SHAKE256_RATE);
  	if(ctr < N) 
	{
    	shake256_squeezeblocks_aigis_mulan(outbuf, 1, state); //fzhang
    	Mulan_Rej_eta(poly->coefficients + ctr, N - ctr, outbuf, SHAKE256_RATE);
  	}
}


// do the rejection sampling. 
unsigned int Mulan_Rej_eta(
		uint32_t * ptr,
		unsigned int len,
		const unsigned char *buf,
		unsigned int buflen)
{
	unsigned int ctr, pos;
	unsigned char t0, t1;
	
	ctr = pos = 0;
	while(ctr<len) 
	{
		t0 = buf[pos] & 0x0F;
		t1 = buf[pos++] >> 4;

		if(t0!=0x0F)
			ptr[ctr++] = uniform_eta_table[t0];

		if(t1!=0x0F && ctr < len)
			ptr[ctr++] = uniform_eta_table[t1];
		
    	if(pos >= buflen)
      		break;
	}

	return ctr;
}

// generate the matrix A based on the seed rho[]. 

void Mulan_Generate_A(
		Mulan_Polynomial matrix[ROWS_A][COLUMNS_A], 
		const unsigned char rho[SEEDBYTES])
{
	uint32_t i, j, k, pos, ctr;
	unsigned char inbuf[SEEDBYTES + 1];

	unsigned char outbuf[5*SHAKE128_RATE]; //  >=910
	uint32_t temp;
	unsigned char collector[280];
	uint32_t c_pos = 0;
	
	memcpy(inbuf, rho, SEEDBYTES);

	for(i=0; i<ROWS_A; i++) 
		for(j=0; j<COLUMNS_A; j++) 
		{
			ctr = pos = 0;
			inbuf[SEEDBYTES] = i + (j << 4);
			shake128(outbuf, sizeof(outbuf), inbuf, SEEDBYTES + 1);
			
			for(k=c_pos=0; k<240; k++)
				collector[k] = 0;

			while(ctr < N && pos+3<5*SHAKE128_RATE)
			{
				temp  = outbuf[pos++]; 
				temp |= (uint32_t)outbuf[pos++] << 8;
				temp |= (uint32_t)outbuf[pos] << 16;
				temp &= 0x1FFFFF;
				
				collector[c_pos] = (outbuf[pos] & 0xE0)>>5;
				pos++; c_pos++;

				if(temp<Q)
					matrix[i][j].coefficients[ctr++] = temp;
			}
			if(ctr==N)
				continue;

			for(k=0; k<238 && ctr<N; k+=7)
			{
				temp = collector[k];
				temp |= collector[k+1]<<3;
				temp |= collector[k+2]<<6;
				temp |= collector[k+3]<<9;
				temp |= collector[k+4]<<12;
				temp |= collector[k+5]<<15;
				temp |= collector[k+6]<<18;
				
				if(temp<Q)
					matrix[i][j].coefficients[ctr++] = temp;
			}
		}
}

void Mulan_Poly_invntt_frominvmont(Mulan_Polynomial *a) 
{
	uint32_t *p = a->coefficients;
	unsigned int start, len, j, k;
	uint32_t t, zeta;
	const uint32_t f = 1478235U;

	k = 0;
	for (len = 1; len < N; len <<= 1) {
		for (start = 0; start < N; start = j + len) {
			zeta = mulan_zetas_inv[k++];
			for (j = start; j < start + len; ++j) {
				t = p[j];
				p[j] = t + p[j + len];
				p[j + len] = t + 256 * Q - p[j + len];
				p[j + len] = Mulan_montgomery_reduce((uint64_t)zeta * p[j + len]);
			}
		}
	}

	for (j = 0; j < N; ++j) {
		p[j] = Mulan_montgomery_reduce((uint64_t)f * p[j]);
	}
}

void Mulan_Poly_pointwise_invmontgomery(Mulan_Polynomial *c, const Mulan_Polynomial *a, const Mulan_Polynomial *b) {
	unsigned int i;

	for (i = 0; i < N; ++i)
		c->coefficients[i] = Mulan_montgomery_reduce((uint64_t)a->coefficients[i] * b->coefficients[i]);
}

// [v1, v0] := Mulan_Decompose(v)
void Mulan_PolyVec_decompose(
		Mulan_Polynomial * v1, 
		Mulan_Polynomial * v0, 
		const Mulan_Polynomial * v, 
		const uint32_t dim)
{
	for(uint32_t j=0; j<dim; j++)
		for(uint32_t k=0; k<N; k++)
			v1[j].coefficients[k] = Mulan_Decompose(v0[j].coefficients+k, v[j].coefficients[k]);
}

// generate c <- B_{60}. 
void Mulan_challenge(
		Mulan_Polynomial * c, 
		const unsigned char mu[CRHBYTES], 
		const Mulan_Polynomial w1[ROWS_K])
{
	unsigned int i, b, pos;
	unsigned char inbuf[CRHBYTES + SIZE_w1];
	unsigned char outbuf[SHAKE256_RATE];
	uint64_t state[25], signs, mask;

	memcpy(inbuf, mu, CRHBYTES);
	
	Mulan_pack_w1(inbuf+CRHBYTES, w1);

	shake256_absorb_aigis_mulan(state, inbuf, sizeof(inbuf)); //fzhang
	shake256_squeezeblocks_aigis_mulan(outbuf, 1, state);	  //fzhang

	for(signs=i=0; i < 8; i++)
		signs |= (uint64_t)outbuf[i] << 8*i;

	pos = 8;
	mask = 1;

	for(i = 0; i < N; i++)
		c->coefficients[i] = 0;

	for(i = 196; i < 256; i++) 
	{
		while(1)
		{
	  		if(pos >= SHAKE256_RATE) 
			{
				shake256_squeezeblocks_aigis_mulan(outbuf, 1, state); //fzhang
				pos = 0;
	  		}
			b = outbuf[pos++];
			if(b<=i)
				break;
		};

		c->coefficients[i] = c->coefficients[b];
		c->coefficients[b] = (signs & mask) ? Q - 1 : 1;
		mask <<= 1;
	}
}

// compute the standard form of every coefficient of ptr[]. 
void Mulan_PolyVec_freeze(Mulan_Polynomial * ptr, uint32_t dim)
{
	uint32_t i, j; 
	uint32_t temp;
	for(i=0; i<dim; i++)
		for(j=0; j<N; j++)
		{
			temp = ptr[i].coefficients[j];
			temp -= ((temp>>21)*Q);
			ptr[i].coefficients[j] = temp - (Q & (((int32_t)(Q-1-temp))>>31));
		}
}



// return 0 iff the absolute value of every coefficient in p[] is < bound. 
int Mulan_PolyVec_check_norm(
		const Mulan_Polynomial * p,
		const uint32_t dim,
		const uint32_t bound)
{
	unsigned int i, j;
  	uint32_t temp; 

	for(i=0; i<dim; i++)
		for(j=0; j<N; j++)
		{
			temp = p[i].coefficients[j];
			temp = (temp - (temp>>21)*Q);
			temp -= (Q & (((int32_t)(Q-1-temp))>>31));
			if(temp>Q/2)
				temp = Q-temp;

			if(temp>= bound)
			  return 1;			
		}

	return 0;
}


// Mulan_HighBits(a) vs. Mulan_HighBits(b)
uint32_t Mulan_PolyVec_make_hint(
		Mulan_Polynomial h[ROWS_K],
		const Mulan_Polynomial a[ROWS_K], 
		const Mulan_Polynomial b[ROWS_K])
{
	uint32_t i, j, counter = 0;

	for(i=0; i<ROWS_K; i++)
		for(j=0; j<N; j++)
		{
			h[i].coefficients[j] = Mulan_Make_hint(a[i].coefficients[j], b[i].coefficients[j]);
			counter += h[i].coefficients[j];
		}

	return counter;
}

// dest := Mulan_Use_hint(r, h)
void Mulan_PolyVec_use_hint(
		Mulan_Polynomial dest[ROWS_K], 
		const Mulan_Polynomial h[ROWS_K], 
		const Mulan_Polynomial r[ROWS_K])
{	
	uint32_t i, j;
	
	for(i=0; i<ROWS_K; i++)
		for(j=0; j<N; j++)
			dest[i].coefficients[j] = Mulan_Use_hint(h[i].coefficients[j], r[i].coefficients[j]);
}

// pack w1 into the string p
void Mulan_pack_w1(
		unsigned char * p,
		const Mulan_Polynomial w1[DIM_w1])
{
	uint32_t j, k;

	for(j=0;j<DIM_w1; j++, p += BYTES_w1)
		for(k=0; k<N/8; k++) 
		{
			p[3*k+0]  = w1[j].coefficients[8*k+0];
			p[3*k+0] |= w1[j].coefficients[8*k+1] << 3;
			p[3*k+0] |= w1[j].coefficients[8*k+2] << 6;
			
			p[3*k+1]  = w1[j].coefficients[8*k+2] >> 2;
			p[3*k+1] |= w1[j].coefficients[8*k+3] << 1;
			p[3*k+1] |= w1[j].coefficients[8*k+4] << 4;
			p[3*k+1] |= w1[j].coefficients[8*k+5] << 7;
			
			p[3*k+2]  = w1[j].coefficients[8*k+5] >> 1;
			p[3*k+2] |= w1[j].coefficients[8*k+6] << 2;
			p[3*k+2] |= w1[j].coefficients[8*k+7] << 5;
		}
}



