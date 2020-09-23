#ifndef ROUNDING_H
#define ROUNDING_H

#include "Parameters.h"

void Mulan_Power2Round(uint32_t *high, uint32_t *low, const uint32_t value);

uint32_t Mulan_Decompose(uint32_t *r0_prime, const uint32_t r);
uint32_t Mulan_HighBits(uint32_t r);

uint32_t Mulan_Make_hint(uint32_t a, uint32_t b); 
uint32_t Mulan_Use_hint(uint32_t h, uint32_t r);


#endif