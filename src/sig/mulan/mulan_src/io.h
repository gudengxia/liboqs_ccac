#ifndef MULAN_IO_H
#define MULAN_IO_H

#include "Polynomial.h"
				
void Mulan_pack_sk1(
		unsigned char * sk1, 
		const unsigned char * rho_A, 
		const unsigned char * key,
		const Mulan_Polynomial s[DIM_s],
		const Mulan_Polynomial e[DIM_e]);
void Mulan_unpack_sk1(
		Mulan_Polynomial A[ROWS_A][COLUMNS_A], unsigned char * rho_A, 
		unsigned char * key,
		Mulan_Polynomial s[DIM_s], 
		Mulan_Polynomial e[DIM_e], 
		const unsigned char * sk1);
		
void Mulan_pack_polyeta(unsigned char *r, const Mulan_Polynomial *a);
void Mulan_unpack_polyeta(Mulan_Polynomial *r, const unsigned char *a);

void Mulan_pack_sk2(
		unsigned char * sk2, 
		unsigned char * tr, 
		const Mulan_Polynomial t0[DIM_t0]);
void Mulan_unpack_sk2(
		unsigned char * tr, 
		Mulan_Polynomial t0[DIM_t0], 
		const unsigned char * sk2);

void Mulan_pack_pk(
		unsigned char * pk, 
		const unsigned char * rho_A, 
		const Mulan_Polynomial t1[DIM_t1]);
void Mulan_unpack_pk(
		Mulan_Polynomial A[ROWS_A][COLUMNS_A], 
		unsigned char * rho_A, 
		Mulan_Polynomial t1[DIM_t1], 
		const unsigned char * pk);		
		
void Mulan_pack_c(
		unsigned char * p, 
		const Mulan_Polynomial * c);
void Mulan_unpack_c(
		Mulan_Polynomial * c, 
		const unsigned char * p);

void Mulan_pack_h(
		unsigned char * p, 
		const Mulan_Polynomial h[DIM_h]);
void Mulan_unpack_h(
		Mulan_Polynomial h[DIM_h], 
		const unsigned char * p);
				
void Mulan_pack_z(
		unsigned char * p, 
		const Mulan_Polynomial z[DIM_z]);
void Mulan_unpack_z(
		Mulan_Polynomial z[DIM_z], 
		const unsigned char * p);
		
void Mulan_pack_signature(
		unsigned char * signature, 
		const Mulan_Polynomial z[DIM_z], 
		const Mulan_Polynomial h[DIM_h], 
		const Mulan_Polynomial *c);
void Mulan_unpack_signature(
		Mulan_Polynomial z[DIM_z], 
		Mulan_Polynomial h[DIM_h], 
		Mulan_Polynomial *c, 
		const unsigned char * signature);


void Mulan_unpack_t0(Mulan_Polynomial t0[DIM_t0], const unsigned char *p);
void Mulan_pack_t0(unsigned char * p, const Mulan_Polynomial t0[DIM_t0]);

void Mulan_pack_polyeta(unsigned char *r, const Mulan_Polynomial *a);
void Mulan_unpack_polyeta(Mulan_Polynomial *r, const unsigned char *a);


#endif