#ifndef MULAN_SIGNATURE_H
#define MULAN_SIGNATURE_H

#include "Parameters.h"
#include "Polynomial.h"
#include "io.h"

#define MULAN_SIG_SECRETKEYBYTES 		MULAN_CRYPTO_SECRETKEYBYTES
#define MULAN_SIG_PUBLICKEYBYTES 		MULAN_CRYPTO_PUBLICKEYBYTES
#define MULAN_SIG_BYTES 				MULAN_CRYPTO_BYTES
#define SIG_ALGNAME 			"Mulan"

#define MULAN_MESSAGE_LENGTH_MAX		3300
#define	MULAN_MESSAGE_LENGTH_IN_TEST	59

int mulan_sig_keygen(
		unsigned char * pk, unsigned long long * pk_bytes, 
		unsigned char * sk, unsigned long long * sk_bytes); 

int mulan_sig_sign(
	unsigned char * sk, unsigned long long sk_bytes,
	unsigned char * m, unsigned long long mlen,
	unsigned char * sm, unsigned long long * smlen);

int mulan_sig_verf(
	unsigned char * pk, unsigned long long pk_bytes,
	unsigned char * sm, unsigned long long smlen,
	unsigned char * m, unsigned long long mlen);

#endif
