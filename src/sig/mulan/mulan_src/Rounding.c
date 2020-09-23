#include "Rounding.h"


// [a1, a0] := Mulan_Power2Round(a)
void Mulan_Power2Round(uint32_t *a1, uint32_t *a0, const uint32_t a) 
{	
	uint32_t temp, t;

	temp = (a - (a>>21)*Q);
	temp -= (Q & (((int32_t)(Q-temp))>>31));
	t = (temp & 0x1FFF); 
	*a0 = Q + t - (0x2000 & ((((int32_t)(0x1000-t))>>31)));
	*a1 = ((temp+Q)-(*a0))>>D;
}

// [r1, r0_prime] <- Mulan_Decompose(r)
// [R1, r0_prime] <- Con(r)
uint32_t Mulan_Decompose(uint32_t *r0_prime, const uint32_t r)
{
	uint32_t b;
	uint32_t r1;
	uint32_t temp;

	b = r - (r>>21)*Q;
	b -= Q & (((int32_t)(Q-1-b))>>31);
	b<<=3;
	temp = b - (b>>21)*Q;
	temp -= Q & (((int32_t)(Q-1-temp))>>31);

	*r0_prime = temp + (Q & (((int32_t)(temp-1-Q/2))>>31));

	temp = Q + b - *r0_prime;
	r1 = (temp>>21);
	r1 += ((uint32_t)((1+r1)*Q-1-temp))>>31;	

	return (r1 & 0x07);
}

// compute the high-bits of r. 
uint32_t Mulan_HighBits(uint32_t r)
{
	uint32_t r0;
	return  Mulan_Decompose(&r0, r);
}


// return 0 iff Mulan_HighBits(a) != Mulan_HighBits(b)
uint32_t Mulan_Make_hint(uint32_t a, uint32_t b)
{
	return (Mulan_HighBits(a) != Mulan_HighBits(b));
}


// return Mulan_Use_hint(r, h)
uint32_t Mulan_Use_hint(uint32_t h, uint32_t r)
{
	uint32_t r1, r0;
	r1 = Mulan_Decompose(&r0, r);

	if (h == 0)
		return r1;

	if (r0>Q)
		return ((r1+1) & 0x07);
	else
		return ((r1+K-1)& 0x07);
}
