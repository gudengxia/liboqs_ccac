// SPDX-License-Identifier: MIT

#ifndef FIPS202_H
#define FIPS202_H

#include <oqs/sha3.h>

#define SHAKE128_RATE OQS_SHA3_SHAKE128_RATE
#define shake128 OQS_SHA3_shake128
#define shake128ctx OQS_SHA3_shake128_ctx
#define shake128_absorb OQS_SHA3_shake128_absorb
#define shake128_squeezeblocks OQS_SHA3_shake128_squeezeblocks
#define shake128_ctx_release OQS_SHA3_shake128_ctx_release
#define shake128_ctx_clone OQS_SHA3_shake128_ctx_clone

#define SHAKE256_RATE OQS_SHA3_SHAKE256_RATE
#define shake256 OQS_SHA3_shake256
#define shake256ctx OQS_SHA3_shake256_ctx
#define shake256_absorb OQS_SHA3_shake256_absorb
#define shake256_squeezeblocks OQS_SHA3_shake256_squeezeblocks
#define shake256_ctx_release OQS_SHA3_shake256_ctx_release
#define shake256_ctx_clone OQS_SHA3_shake256_ctx_clone

#define SHA3_256_RATE OQS_SHA3_SHA3_256_RATE
#define sha3_256 OQS_SHA3_sha3_256
#define sha3_256_inc_init OQS_SHA3_sha3_256_inc_init
#define sha3_256_inc_absorb OQS_SHA3_sha3_256_inc_absorb
#define sha3_256_inc_finalize OQS_SHA3_sha3_256_inc_finalize
#define sha3_256_inc_ctx_clone OQS_SHA3_sha3_256_inc_ctx_clone
#define sha3_256_inc_ctx_release OQS_SHA3_sha3_256_inc_ctx_release

#define SHA3_384_RATE OQS_SHA3_SHA3_384_RATE
#define sha3_384 OQS_SHA3_sha3_384
#define sha3_384_inc_init OQS_SHA3_sha3_384_inc_init
#define sha3_384_inc_absorb OQS_SHA3_sha3_384_inc_absorb
#define sha3_384_inc_finalize OQS_SHA3_sha3_384_inc_finalize
#define sha3_384_inc_ctx_clone OQS_SHA3_sha3_384_inc_ctx_clone
#define sha3_384_inc_ctx_release OQS_SHA3_sha3_384_inc_ctx_release

#define SHA3_512_RATE OQS_SHA3_SHA3_512_RATE
#define sha3_512 OQS_SHA3_sha3_512
#define sha3_512_inc_init OQS_SHA3_sha3_512_inc_init
#define sha3_512_inc_absorb OQS_SHA3_sha3_512_inc_absorb
#define sha3_512_inc_finalize OQS_SHA3_sha3_512_inc_finalize
#define sha3_512_inc_ctx_clone OQS_SHA3_sha3_512_inc_ctx_clone
#define sha3_512_inc_ctx_release OQS_SHA3_sha3_512_inc_ctx_release

#define shake128incctx OQS_SHA3_shake128_inc_ctx
#define shake128_inc_init OQS_SHA3_shake128_inc_init
#define shake128_inc_absorb OQS_SHA3_shake128_inc_absorb
#define shake128_inc_finalize OQS_SHA3_shake128_inc_finalize
#define shake128_inc_squeeze OQS_SHA3_shake128_inc_squeeze
#define shake128_inc_ctx_release OQS_SHA3_shake128_inc_ctx_release
#define shake128_inc_ctx_clone OQS_SHA3_shake128_inc_ctx_clone

#define shake256incctx OQS_SHA3_shake256_inc_ctx
#define shake256_inc_init OQS_SHA3_shake256_inc_init
#define shake256_inc_absorb OQS_SHA3_shake256_inc_absorb
#define shake256_inc_finalize OQS_SHA3_shake256_inc_finalize
#define shake256_inc_squeeze OQS_SHA3_shake256_inc_squeeze
#define shake256_inc_ctx_release OQS_SHA3_shake256_inc_ctx_release
#define shake256_inc_ctx_clone OQS_SHA3_shake256_inc_ctx_clone

/**************************************fzhang start********************************************/
void shake256_absorb_aigis_mulan(uint64_t *s,const unsigned char *input, unsigned long long inlen); 
void shake256_squeezeblocks_aigis_mulan(unsigned char *output, unsigned long nblocks, uint64_t *s);
/**************************************fzhang end**********************************************/
#endif