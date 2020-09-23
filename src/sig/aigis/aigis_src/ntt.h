#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"



void aigis_ntt(uint32_t p[PARAM_N]);
void aigis_invntt(uint32_t p[PARAM_N]);
#endif
