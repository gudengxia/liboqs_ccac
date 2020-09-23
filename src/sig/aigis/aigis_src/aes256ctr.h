#ifndef AIGIS_AES256CTR_H
#define AIGIS_AES256CTR_H

#include <stdint.h>
#include <immintrin.h>

typedef struct {
  __m128i rkeys[16];
  __m128i n;
} aigis_aes256ctr_ctx;

void aigis_aes256ctr_init(aigis_aes256ctr_ctx *state,
                    const unsigned char *key,
                    uint16_t nonce);
void aigis_aes256ctr_select(aigis_aes256ctr_ctx *state, uint16_t nonce);
void aigis_aes256ctr_squeezeblocks(unsigned char *out,
                             unsigned long long nblocks,
                             aigis_aes256ctr_ctx *state);

void aigis_aes256ctr_prf(unsigned char *out,
                   unsigned long long outlen,
                   const unsigned char *seed,
                   unsigned char nonce);

#endif
