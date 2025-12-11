#ifndef NTT_H
#define NTT_H

#include <stdint.h>

#if (CTRU_N == 512)
#define ROOT_DIMENSION 8
#elif (CTRU_N == 1024)
#define ROOT_DIMENSION 16
#endif

// #define LEVEL_ONE_FACTOR 154
#define LEVEL_ONE_FACTOR 62   // 5^32 mod 769

extern int16_t zetas[64];
extern int16_t zetas_inv[64];

void ntt(int16_t a[]);
void ntt_big(int16_t *a, const int16_t *in);
void ntt_small(int16_t *a, const int16_t *in);
// void invntt(int16_t *a, const int16_t *in);
void invntt(int16_t a[]);
void invntt_without_mont(int16_t *a, const int16_t *in);
void basemul(int16_t *c, const int16_t *a, const int16_t *b, const int16_t zeta);

#endif