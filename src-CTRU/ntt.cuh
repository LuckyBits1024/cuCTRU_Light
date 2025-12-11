#ifndef CTRU_NTT_H
#define NTT_H

#include <stdint.h>

#if (CTRU_N == 512)
#define ROOT_DIMENSION 8
#elif (CTRU_N == 1024)
#define ROOT_DIMENSION 16
#endif

// #define LEVEL_ONE_FACTOR 154
#define LEVEL_ONE_FACTOR 62   // 5^32 mod 769

extern __device__ int16_t zetas[64];
extern __device__ int16_t zetas_inv[64];

__device__ void ntt_big(int16_t *a, const int16_t *in);
__device__ void ntt_small(int16_t *a, const int16_t *in);
//__device__ void invntt(int16_t *a, const int16_t *in);
__device__ void invntt_without_mont(int16_t *a, const int16_t *in);
//__device__ void basemul(int16_t *c, const int16_t *a, const int16_t *b, const int16_t zeta);

__device__ void ntt(int16_t regs[8], int16_t *s_ntt);
__device__ void inv_ntt(int16_t regs[8], int16_t *s_ntt);
__device__ void basemul(int16_t a[8], int16_t b[8], int16_t c[8], int16_t zeta);
__device__ void base_mul(int16_t a[8], int16_t b[8], int16_t zeta);
#endif