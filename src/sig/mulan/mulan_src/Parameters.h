#ifndef MULAN_PARAMETERS_H
#define MULAN_PARAMETERS_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "fips202.h" 
#include "randombytes.h"

#define N			256U
#define	Q		1952257U

#define ROWS_K 		5U
#define COLUMNS_ELL 4U

#define	ROWS_A		(ROWS_K)
#define	COLUMNS_A 	(COLUMNS_ELL)

#define SEEDBYTES	32U	
#define CRHBYTES 	48U


#define DIM_s		(COLUMNS_ELL)
#define	DIM_e		(ROWS_K)
#define	DIM_t		(ROWS_K)
#define	DIM_t1		(ROWS_K)
#define	DIM_t0		(ROWS_K)
#define	DIM_c		(1U)
#define	DIM_y		(COLUMNS_ELL)
#define	DIM_z		(COLUMNS_ELL)
#define DIM_w 		(ROWS_K)
#define DIM_w1 		(ROWS_K)
#define DIM_h		(ROWS_K)

#define D			13U
#define K			8U
#define ETA_s		2U
#define ETA_e		2U
#define	ETA			2U
#define	U_s			120U
#define	U_e			120U

#define BITS_s		(3U)
#define BITS_e		(3U)
#define BITS_t1 	(8U)
#define	BITS_t0		(13U)
#define BITS_z		(19U)
#define BITS_w1		(3U)


#define BYTES_s 	(BITS_s<<5)
#define	BYTES_e		(BITS_e<<5)
#define BYTES_t1	(BITS_t1<<5)
#define	BYTES_t0	(BITS_t0<<5)
#define	BYTES_z		(BITS_z<<5)
#define BYTES_w1	(BITS_w1<<5)


#define OMEGA 		96U
#define	GAMMA1 		(Q/K)

#define SIZE_s		(DIM_s * BYTES_s)
#define SIZE_e		(DIM_e * BYTES_e)	
#define SIZE_t1 	(DIM_t1 * BYTES_t1)
#define SIZE_t0 	(DIM_t0 * BYTES_t0)	
#define SIZE_w1 	(DIM_w1 * BYTES_w1)
#define SIZE_z		(DIM_z * BYTES_z)
#define SIZE_c		((N>>3)+8)
#define SIZE_h		(OMEGA + DIM_h)

#define	SIZE_sk1 	(SEEDBYTES + SEEDBYTES + SIZE_s + SIZE_e)
#define SIZE_sk2 	(CRHBYTES + SIZE_t0)
#define SIZE_sk 	(SIZE_sk1 + SIZE_sk2)
#define SIZE_pk 	(SEEDBYTES + SIZE_t1)


#define MULAN_CRYPTO_BYTES 			(SIZE_z+SIZE_h+SIZE_c)
#define MULAN_CRYPTO_PUBLICKEYBYTES 	SIZE_pk
#define	MULAN_CRYPTO_SECRETKEYBYTES 	SIZE_sk


#define NORMBOUND_z		(Q/K-U_s)
#define NORMBOUND_r0	(Q/2-K*U_e)
#define	NORMBOUND_ct0	(Q/(2*K))

#endif