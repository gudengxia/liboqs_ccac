#ifndef PARAMS_H
#define PARAMS_H


#define PARAMS 2

#define USE_AES //use AES-NI

#define SEEDBYTES 32
#define CRHBYTES 48
#define PARAM_N 256

#if PARAMS == 1
#define PARAM_Q 2021377
#define MONT 1562548U // 2^32 % Q
#define QINV 2849953791U // -q^(-1) mod 2^32
#define QBITS 21
#define ROOT_OF_UNITY 79
#define PARAM_D 13
#define GAMMA1 131072
#define GAMMA2 168448
#define ALPHA (2*GAMMA2)


#define PARAM_K 4
#define PARAM_L 3
#define ETA1 2
#define ETA2 3
#define SETA1BITS 3
#define SETA2BITS 3

#define BETA1 120
#define BETA2 175
#define OMEGA 80

#define POLZ_SIZE_PACKED (PARAM_N*18/8)
#define POLW1_SIZE_PACKED ((PARAM_N*3)/8)

#elif PARAMS == 2
#define PARAM_Q 3870721
#define MONT 2337707U // 2^32 % Q
#define QINV 2671448063U // -q^(-1) mod 2^32
#define QBITS 22
#define ROOT_OF_UNITY 19602
#define PARAM_D 14
#define GAMMA1 131072
#define GAMMA2 322560
#define ALPHA (2*GAMMA2)


#define PARAM_K 5
#define PARAM_L 4
#define ETA1 2
#define ETA2 5
#define SETA1BITS 3
#define SETA2BITS 4
#define BETA1 120
#define BETA2 275
#define OMEGA 96

#define POLZ_SIZE_PACKED (PARAM_N*18/8)
#define POLW1_SIZE_PACKED ((PARAM_N*3)/8)


#elif PARAMS == 22
#define PARAM_Q 3870721
#define MONT 2337707U // 2^32 % Q
#define QINV 2671448063U // -q^(-1) mod 2^32
#define QBITS 22
#define ROOT_OF_UNITY 19602
#define PARAM_D 14
#define GAMMA1 131072
#define GAMMA2 322560
#define ALPHA (2*GAMMA2)


#define PARAM_K 5
#define PARAM_L 4
#define ETA1 3
#define ETA2 5
#define SETA1BITS 3
#define SETA2BITS 4
#define BETA1 170
#define BETA2 275
#define OMEGA 96

#define POLZ_SIZE_PACKED (PARAM_N*18/8)
#define POLW1_SIZE_PACKED ((PARAM_N*3)/8)

#elif PARAMS == 3
#define PARAM_Q 3870721
#define MONT 2337707U // 2^32 % Q
#define QINV 2671448063U // -q^(-1) mod 2^32
#define QBITS 22
#define ROOT_OF_UNITY 19602
#define PARAM_D 14
#define GAMMA1 131072
#define GAMMA2 322560
#define ALPHA (2*GAMMA2)

#define PARAM_K 6
#define PARAM_L 5
#define ETA1 1
#define ETA2 5
#define SETA1BITS 2
#define SETA2BITS 4
#define BETA1 60
#define BETA2 275
#define OMEGA 120

#define POLZ_SIZE_PACKED (PARAM_N*18/8)
#define POLW1_SIZE_PACKED ((PARAM_N*3)/8)
#endif


#define POL_SIZE_PACKED ((PARAM_N*QBITS)/8)
#define POLT1_SIZE_PACKED ((PARAM_N*(QBITS - PARAM_D))/8)
#define POLT0_SIZE_PACKED ((PARAM_N*PARAM_D)/8)
#define POLETA1_SIZE_PACKED ((PARAM_N*SETA1BITS)/8)
#define POLETA2_SIZE_PACKED ((PARAM_N*SETA2BITS)/8)

#define POLVECK_SIZE_PACKED (PARAM_K*POL_SIZE_PACKED)
#define POLVECL_SIZE_PACKED (PARAM_L*POL_SIZE_PACKED)
#define AIGIS_PK_SIZE_PACKED (SEEDBYTES + PARAM_K*POLT1_SIZE_PACKED)
#define AIGIS_SK_SIZE_PACKED (2*SEEDBYTES + PARAM_L*POLETA1_SIZE_PACKED + PARAM_K*POLETA2_SIZE_PACKED + CRHBYTES + PARAM_K*POLT0_SIZE_PACKED)
#define AIGIS_SIG_SIZE_PACKED (PARAM_L*POLZ_SIZE_PACKED + (OMEGA + PARAM_K) + (PARAM_N/8 + 8))
#endif
