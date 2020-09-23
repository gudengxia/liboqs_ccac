#include "io.h"


// pack sk1 := (seed_A, seed_key, s, e)
void Mulan_pack_sk1(
		unsigned char * sk1, 
		const unsigned char * rho_A, 
		const unsigned char * key,
		const Mulan_Polynomial s[DIM_s],
		const Mulan_Polynomial e[DIM_e])
{
	unsigned char * p = sk1;
	uint32_t i = 0;
	
	memcpy(p, rho_A, SEEDBYTES); 	
	memcpy(p+SEEDBYTES, key, SEEDBYTES);
	p += (SEEDBYTES*2);
	
	for(i=0; i<DIM_s; i++, p += BYTES_s)
		Mulan_pack_polyeta(p, s+i);
	
	for(i=0; i<DIM_e; i++, p += BYTES_e)
		Mulan_pack_polyeta(p, e+i);
}

// unpack sk1 := (rho_A, key, s, e)
void Mulan_unpack_sk1(
		Mulan_Polynomial A[ROWS_A][COLUMNS_A], unsigned char * rho_A, 
		unsigned char * key,
		Mulan_Polynomial s[DIM_s], 
		Mulan_Polynomial e[DIM_e], 
		const unsigned char * sk1)
{
	uint32_t i;
	unsigned char * p = (unsigned char *)sk1;

	for(i=0; i<SEEDBYTES; i++)
		rho_A[i] = p[i];
	for(i=0; i<SEEDBYTES; i++)
		key[i] = p[SEEDBYTES+i];
	p += 2*SEEDBYTES;
	
	Mulan_Generate_A(A, rho_A);
	
	for(i=0; i<DIM_s; i++, p += BYTES_s)
		Mulan_unpack_polyeta(s+i, p);
	
	for(i=0; i<DIM_e; i++, p += BYTES_e)
		Mulan_unpack_polyeta(e+i, p);
}

void Mulan_pack_polyeta(unsigned char *r, const Mulan_Polynomial *a)
{ 
	uint32_t i;
	unsigned char t[8];
	
  	for(i=0; i<N/8; i++) 
	{
		t[0] = Q + ETA - a->coefficients[8*i+0];
		t[1] = Q + ETA - a->coefficients[8*i+1];
		t[2] = Q + ETA - a->coefficients[8*i+2];
		t[3] = Q + ETA - a->coefficients[8*i+3];
		t[4] = Q + ETA - a->coefficients[8*i+4];
		t[5] = Q + ETA - a->coefficients[8*i+5];
		t[6] = Q + ETA - a->coefficients[8*i+6];
		t[7] = Q + ETA - a->coefficients[8*i+7];

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
}

void Mulan_unpack_polyeta(Mulan_Polynomial *r, const unsigned char *a)
{
	uint32_t * value = r->coefficients; 

	for(uint32_t i = 0; i < N/8; ++i) 
	{
		value[8*i+0] = a[3*i+0] & 0x07;
		value[8*i+1] = (a[3*i+0] >> 3) & 0x07;
		value[8*i+2] = (a[3*i+0] >> 6) | ((a[3*i+1] & 0x01) << 2);
		value[8*i+3] = (a[3*i+1] >> 1) & 0x07;
		value[8*i+4] = (a[3*i+1] >> 4) & 0x07;
		value[8*i+5] = (a[3*i+1] >> 7) | ((a[3*i+2] & 0x03) << 1);
		value[8*i+6] = (a[3*i+2] >> 2) & 0x07;
		value[8*i+7] = (a[3*i+2] >> 5);

		value[8*i+0] = Q + ETA - value[8*i+0];
		value[8*i+1] = Q + ETA - value[8*i+1];
		value[8*i+2] = Q + ETA - value[8*i+2];
		value[8*i+3] = Q + ETA - value[8*i+3];
		value[8*i+4] = Q + ETA - value[8*i+4];
		value[8*i+5] = Q + ETA - value[8*i+5];
		value[8*i+6] = Q + ETA - value[8*i+6];
		value[8*i+7] = Q + ETA - value[8*i+7];
	}

}


// pack pk := (seed_A, t1)
void Mulan_pack_pk(
		unsigned char * pk, 
		const unsigned char * rho_A, 
		const Mulan_Polynomial t1[DIM_t1])
{
	unsigned char * p = pk;
	uint32_t j, k;
	
	memcpy(p, rho_A, SEEDBYTES);
	p+=SEEDBYTES;
	
	for(j=0; j<DIM_t1; j++)
		for(k=0; k<N; k++)
			*p++ = (unsigned char)((t1[j].coefficients[k]) & 0xFF);
}

// unpack pk := (seed_A, t1)
void Mulan_unpack_pk(
		Mulan_Polynomial A[ROWS_A][COLUMNS_A], unsigned char * rho_A, 
		Mulan_Polynomial t1[DIM_t1], 
		const unsigned char * pk)
{
	uint32_t j, k;
	unsigned char * p = (unsigned char*)pk;
	
	memcpy(rho_A, p, SEEDBYTES);
	p += SEEDBYTES;
	Mulan_Generate_A(A, rho_A);

	for(j=0; j<DIM_t1; j++)
		for(k=0; k<N; k++)
			t1[j].coefficients[k] = (uint32_t)(*p++)<<D;
}


// pack sk2 := (tr, t0)
void Mulan_pack_sk2(
		unsigned char * sk2, 
		unsigned char * tr, 
		const Mulan_Polynomial t0[DIM_t0])
{
	unsigned char * p = sk2; 
	
	if(p!=tr)
		memcpy(p, tr, CRHBYTES);
	p += CRHBYTES;
	
	Mulan_pack_t0(p, t0);
}


// unpack sk2 := (tr, t0)
void Mulan_unpack_sk2(
		unsigned char * tr, 
		Mulan_Polynomial t0[DIM_t0], 
		const unsigned char * sk2)
{
	unsigned char * p = (unsigned char *)sk2; 

	if(tr!=p)
		memcpy(tr, p, CRHBYTES);
	p += CRHBYTES;

	Mulan_unpack_t0(t0, p);
}


void Mulan_pack_t0(unsigned char * p, const Mulan_Polynomial t0[DIM_t0])
{
	uint32_t j, k;
	uint32_t a[8];
	
	for(j=0; j<DIM_t0; j++)
		for(k=0; k<N; k+=8)
		{
			// p = ptr+j*BYTES_t0+(k>>3)*BITS_t0+;
			
			a[0] = (Q+4096) - t0[j].coefficients[k];
			a[1] = (Q+4096) - t0[j].coefficients[k+1];
			a[2] = (Q+4096) - t0[j].coefficients[k+2];
			a[3] = (Q+4096) - t0[j].coefficients[k+3];
			a[4] = (Q+4096) - t0[j].coefficients[k+4];
			a[5] = (Q+4096) - t0[j].coefficients[k+5];
			a[6] = (Q+4096) - t0[j].coefficients[k+6];
			a[7] = (Q+4096) - t0[j].coefficients[k+7];
			
			p[j*BYTES_t0+(k>>3)*BITS_t0+0] = (a[0] & 0xFF);
			p[j*BYTES_t0+(k>>3)*BITS_t0+1] = ((a[0] & 0x1F00)>>8);
			p[j*BYTES_t0+(k>>3)*BITS_t0+1] |= ((a[1]&0x7)<<5);
			p[j*BYTES_t0+(k>>3)*BITS_t0+2] = (a[1] & 0x7F8)>>3;
			p[j*BYTES_t0+(k>>3)*BITS_t0+3] = (a[1] & 0x1800)>>11;
			p[j*BYTES_t0+(k>>3)*BITS_t0+3] |= (a[2] & 0x3F)<<2;
			p[j*BYTES_t0+(k>>3)*BITS_t0+4] = (a[2] & 0x1FC0)>>6;
			p[j*BYTES_t0+(k>>3)*BITS_t0+4] |= (a[3] & 0x01)<<7;
			p[j*BYTES_t0+(k>>3)*BITS_t0+5] = (a[3] & 0x01FE)>>1;
			p[j*BYTES_t0+(k>>3)*BITS_t0+6] = (a[3] & 0x1E00)>>9;
			p[j*BYTES_t0+(k>>3)*BITS_t0+6] |= (a[4] & 0x0F)<<4;
			p[j*BYTES_t0+(k>>3)*BITS_t0+7] = (a[4] & 0x0FF0)>>4;
			p[j*BYTES_t0+(k>>3)*BITS_t0+8] = (a[4] & 0x1000)>>12;
			p[j*BYTES_t0+(k>>3)*BITS_t0+8] |= (a[5] & 0x7F)<<1;
			p[j*BYTES_t0+(k>>3)*BITS_t0+9] = (a[5] & 0x1F80)>>7;
			p[j*BYTES_t0+(k>>3)*BITS_t0+9] |= (a[6] & 0x03)<<6;
			p[j*BYTES_t0+(k>>3)*BITS_t0+10] = (a[6] & 0x03FC)>>2;
			p[j*BYTES_t0+(k>>3)*BITS_t0+11] = (a[6] & 0x1C00)>>10;
			p[j*BYTES_t0+(k>>3)*BITS_t0+11] |= (a[7] & 0x1F)<<3;
			p[j*BYTES_t0+(k>>3)*BITS_t0+12] = (a[7] & 0x1FE0)>>5;
		}
	
}

void Mulan_unpack_t0(Mulan_Polynomial t0[DIM_t0], const unsigned char *p)
{
	uint32_t j, k;
 	uint32_t b[8];
	
	for(j=0; j<DIM_t0; j++)
		for(k=0; k<N; k+=8)
		{
			b[0] = (p[j*BYTES_t0+(k>>3)*BITS_t0+0] & 0xFF);
			b[0] |= (p[j*BYTES_t0+(k>>3)*BITS_t0+1] & 0x1F)<<8;
			b[1] = (p[j*BYTES_t0+(k>>3)*BITS_t0+1] & 0xE0)>>5;
			b[1] |= (p[j*BYTES_t0+(k>>3)*BITS_t0+2]<<3);
			b[1] |= (p[j*BYTES_t0+(k>>3)*BITS_t0+3] & 0x03)<<11;
			b[2] = (p[j*BYTES_t0+(k>>3)*BITS_t0+3] & 0x0FC)>>2;
			b[2] |= (p[j*BYTES_t0+(k>>3)*BITS_t0+4] & 0x07F)<<6;
			b[3] = (p[j*BYTES_t0+(k>>3)*BITS_t0+4] & 0x80)>>7;
			b[3] |= (p[j*BYTES_t0+(k>>3)*BITS_t0+5]<<1);
			b[3] |= (p[j*BYTES_t0+(k>>3)*BITS_t0+6] & 0x0F)<<9;
			b[4] = (p[j*BYTES_t0+(k>>3)*BITS_t0+6] & 0xF0)>>4;
			b[4] |= (p[j*BYTES_t0+(k>>3)*BITS_t0+7]<<4);
			b[4] |= (p[j*BYTES_t0+(k>>3)*BITS_t0+8] & 0x01)<<12;
			b[5] = (p[j*BYTES_t0+(k>>3)*BITS_t0+8] & 0xFE)>>1;
			b[5] |= (p[j*BYTES_t0+(k>>3)*BITS_t0+9] & 0x3F)<<7;
			b[6] = (p[j*BYTES_t0+(k>>3)*BITS_t0+9] & 0xC0)>>6;
			b[6] |= (p[j*BYTES_t0+(k>>3)*BITS_t0+10] & 0xFF)<<2;
			b[6] |= (p[j*BYTES_t0+(k>>3)*BITS_t0+11] & 0x07)<<10;
			b[7] = (p[j*BYTES_t0+(k>>3)*BITS_t0+11] & 0xF8)>>3;
			b[7] |= (p[j*BYTES_t0+(k>>3)*BITS_t0+12] & 0xFF)<<5;

			t0[j].coefficients[k+0] = (Q+4096) - b[0];
			t0[j].coefficients[k+1] = (Q+4096) - b[1];
			t0[j].coefficients[k+2] = (Q+4096) - b[2];
			t0[j].coefficients[k+3] = (Q+4096) - b[3];
			t0[j].coefficients[k+4] = (Q+4096) - b[4];
			t0[j].coefficients[k+5] = (Q+4096) - b[5];
			t0[j].coefficients[k+6] = (Q+4096) - b[6];
			t0[j].coefficients[k+7] = (Q+4096) - b[7];
		}
}

void Mulan_pack_z(
		unsigned char * p, 
		const Mulan_Polynomial z[DIM_z])
{
	uint32_t j, k;
	uint32_t a[8], b[8];
	
	for(j=0; j<DIM_z; j++)
		for(k=0; k<N; k+=8)
		{			
			b[0] = z[j].coefficients[k+0];
			b[0] = (b[0]+(Q & (((int32_t)(b[0]-Q/2))>>31)));
			a[0] = Q+(Q/K-U_s) - b[0];
			
			b[1] = z[j].coefficients[k+1];
			b[1] = (b[1]+(Q & (((int32_t)(b[1]-Q/2))>>31)));
			a[1] = Q+(Q/K-U_s) - b[1];

			b[2] = z[j].coefficients[k+2];
			b[2] = (b[2]+(Q & (((int32_t)(b[2]-Q/2))>>31)));
			a[2] = Q+(Q/K-U_s) - b[2];

			b[3] = z[j].coefficients[k+3];
			b[3] = (b[3]+(Q & (((int32_t)(b[3]-Q/2))>>31)));
			a[3] = Q+(Q/K-U_s) - b[3];
			
			b[4] = z[j].coefficients[k+4];
			b[4] = (b[4]+(Q & (((int32_t)(b[4]-Q/4))>>31)));
			a[4] = Q+(Q/K-U_s) - b[4];
			
			b[5] = z[j].coefficients[k+5];
			b[5] = (b[5]+(Q & (((int32_t)(b[5]-Q/2))>>31)));
			a[5] = Q+(Q/K-U_s) - b[5];

			b[6] = z[j].coefficients[k+6];
			b[6] = (b[6]+(Q & (((int32_t)(b[6]-Q/2))>>31)));
			a[6] = Q+(Q/K-U_s) - b[6];

			b[7] = z[j].coefficients[k+7];
			b[7] = (b[7]+(Q & (((int32_t)(b[7]-Q/2))>>31)));
			a[7] = Q+(Q/K-U_s) - b[7];
			
			p[j*BYTES_z+(k>>3)*BITS_z+0] = (a[0] >> 11);
			p[j*BYTES_z+(k>>3)*BITS_z+1] = (a[0] >> 3);
			p[j*BYTES_z+(k>>3)*BITS_z+2] = (a[1] >> 11);
			p[j*BYTES_z+(k>>3)*BITS_z+3] = (a[1] >> 3);
			p[j*BYTES_z+(k>>3)*BITS_z+4] = (a[2] >> 11);
			p[j*BYTES_z+(k>>3)*BITS_z+5] = (a[2] >> 3);
			p[j*BYTES_z+(k>>3)*BITS_z+6] = (a[3] >> 11);
			p[j*BYTES_z+(k>>3)*BITS_z+7] = (a[3] >> 3);
			p[j*BYTES_z+(k>>3)*BITS_z+8] = (a[4] >> 11);
			p[j*BYTES_z+(k>>3)*BITS_z+9] = (a[4] >> 3);
			p[j*BYTES_z+(k>>3)*BITS_z+10] = (a[5] >> 11);
			p[j*BYTES_z+(k>>3)*BITS_z+11] = (a[5] >> 3);
			p[j*BYTES_z+(k>>3)*BITS_z+12] = (a[6] >> 11);
			p[j*BYTES_z+(k>>3)*BITS_z+13] = (a[6] >> 3);
			p[j*BYTES_z+(k>>3)*BITS_z+14] = (a[7] >> 11);
			p[j*BYTES_z+(k>>3)*BITS_z+15] = (a[7] >> 3);
			p[j*BYTES_z+(k>>3)*BITS_z+16] = ((a[0] & 7)<<5) | ((a[1] & 7)<<2) | ((a[2] & 6)>>1);
			p[j*BYTES_z+(k>>3)*BITS_z+17] =  ((a[2] & 1)<<7)|(( a[3] & 7)<<4)|(( a[4] & 7)<<1)|((a[5] & 4)>>2);
			p[j*BYTES_z+(k>>3)*BITS_z+18] = ((a[5] & 3)<<6) | ((a[6] & 7)<<3)|(a[7] & 7);
		}
}


void Mulan_unpack_z(
		Mulan_Polynomial z[DIM_z], 
		const unsigned char * p)
{
	uint32_t j, k;
	uint32_t b[8];
	
	for(j=0; j<DIM_z; j++)
		for(k=0; k<N; k+=8)
		{	
			b[0] = (p[j*BYTES_z+(k>>3)*BITS_z+0] << 11) | (p[j*BYTES_z+(k>>3)*BITS_z+1] << 3) | ((p[j*BYTES_z+(k>>3)*BITS_z+16] >> 5) & 7);
			b[1] = (p[j*BYTES_z+(k>>3)*BITS_z+2] << 11) | (p[j*BYTES_z+(k>>3)*BITS_z+3] << 3) | ((p[j*BYTES_z+(k>>3)*BITS_z+16] & 0x1C) >> 2);
			b[2] = (p[j*BYTES_z+(k>>3)*BITS_z+4] << 11) | (p[j*BYTES_z+(k>>3)*BITS_z+5] << 3) | ((p[j*BYTES_z+(k>>3)*BITS_z+16] & 3) << 1) | (p[j*BYTES_z+(k>>3)*BITS_z+17] >> 7);
			b[3] = (p[j*BYTES_z+(k>>3)*BITS_z+6] << 11) | (p[j*BYTES_z+(k>>3)*BITS_z+7] << 3) | ((p[j*BYTES_z+(k>>3)*BITS_z+17] & 0x70) >> 4);
			b[4] = (p[j*BYTES_z+(k>>3)*BITS_z+8] << 11) | (p[j*BYTES_z+(k>>3)*BITS_z+9] << 3) | ((p[j*BYTES_z+(k>>3)*BITS_z+17] & 0xE) >> 1);
			b[5] = (p[j*BYTES_z+(k>>3)*BITS_z+10] << 11) | (p[j*BYTES_z+(k>>3)*BITS_z+11] << 3) | ((p[j*BYTES_z+(k>>3)*BITS_z+17] & 1) << 2 )| ((p[j*BYTES_z+(k>>3)*BITS_z+18] & 0xc0) >> 6);
			b[6] = (p[j*BYTES_z+(k>>3)*BITS_z+12] << 11) | (p[j*BYTES_z+(k>>3)*BITS_z+13] << 3) | ((p[j*BYTES_z+(k>>3)*BITS_z+18] & 0x38) >> 3);
			b[7] = (p[j*BYTES_z+(k>>3)*BITS_z+14] << 11) | (p[j*BYTES_z+(k>>3)*BITS_z+15] << 3) | (p[j*BYTES_z+(k>>3)*BITS_z+18] & 7);
			
			z[j].coefficients[k+0] = Q+(Q/K-U_s) - b[0];
			z[j].coefficients[k+1] = Q+(Q/K-U_s) - b[1];
			z[j].coefficients[k+2] = Q+(Q/K-U_s) - b[2];
			z[j].coefficients[k+3] = Q+(Q/K-U_s) - b[3];
			z[j].coefficients[k+4] = Q+(Q/K-U_s) - b[4];
			z[j].coefficients[k+5] = Q+(Q/K-U_s) - b[5];
			z[j].coefficients[k+6] = Q+(Q/K-U_s) - b[6];
			z[j].coefficients[k+7] = Q+(Q/K-U_s) - b[7];
		}
}


// pack h into the string p.
void Mulan_pack_h(
		unsigned char * p, 
		const Mulan_Polynomial h[DIM_h])
{
	uint32_t i, j, k;
	
	for(k=i=0; i<DIM_h; i++) 
	{
		for(j=0; j<N; j++)
			if(h[i].coefficients[j] == 1)
				p[k++] = j;

		p[OMEGA+i] = k;
	}

	while(k<OMEGA) p[k++] = 0;
}


// unpack h from the string p.
void Mulan_unpack_h(
		Mulan_Polynomial h[DIM_h], 
		const unsigned char * p)
{
	uint32_t i, j, k; 

	for(k=i=0; i<DIM_h; i++) 
	{
		for(j=0; j<N; h[i].coefficients[j++] = 0);
		for(j=k; j<p[OMEGA+i]; h[i].coefficients[p[j++]] = 1);
		k = p[OMEGA+i];
	}
}


// pack c into the string p. 
void Mulan_pack_c(unsigned char * p, const Mulan_Polynomial *c)
{
	uint64_t signs, mask;
	uint32_t i, j; 
	
	for(signs=i=0, mask=1; i<(N>>3); i++) 
		for(p[i]=j=0; j<8; j++) 
			if(c->coefficients[(i<<3)+j] != 0) 
			{
				p[i] |= (1<<j);
				if(c->coefficients[(i<<3)+j] == (Q - 1)) 
					signs |= mask;
				mask <<= 1;
			}

	for(p += (N>>3), i=0; i<8; i++)
		p[i] = (signs>>(i<<3));
}


// unpack c from the string p. 
void Mulan_unpack_c(Mulan_Polynomial * c, const unsigned char * p)
{
	uint64_t signs, mask;
	uint32_t i, j; 

	for(i=0; i<N; c->coefficients[i++] = 0);

	for(signs=i=0; i<8; i++)
		signs |= (uint64_t)p[(N>>3)+i] << (i<<3);

	for(mask=1, i=0; i<(N>>3); i++) 
		for(j=0; j<8; j++) 
	  		if((p[i] >> j) & 0x01) 
			{
				c->coefficients[(i<<3)+j] = (signs & mask) ? (Q-1) : 1;
				mask <<= 1;
	  		}
}


// pack signature := (z,h,c)
void Mulan_pack_signature(
		unsigned char * signature, 
		const Mulan_Polynomial z[DIM_z], 
		const Mulan_Polynomial h[DIM_h], 
		const Mulan_Polynomial *c)
{
	unsigned char * p = signature; 
	
	Mulan_pack_z(p, z);
	p += SIZE_z; Mulan_pack_h(p, h);
	p += SIZE_h; Mulan_pack_c(p, c);
}


// unpack signature := (z,h,c)
void Mulan_unpack_signature(
		Mulan_Polynomial z[DIM_z], 
		Mulan_Polynomial h[DIM_h], 
		Mulan_Polynomial *c, 
		const unsigned char * signature)
{
	unsigned char * p = (unsigned char*)signature; 

	Mulan_unpack_z(z, p);
	p += SIZE_z; Mulan_unpack_h(h, p);
	p += SIZE_h; Mulan_unpack_c(c, p);
}
