#ifndef NTT_H
#define NTT_H

#include <stdint.h>

#define ROOT_DIMENSION 8

//#define LEVEL_ONE_FACTOR 154
#define LEVEL_ONE_FACTOR 62   // 5^32 mod 769

extern int16_t zetas[64];
extern int16_t zetas_inv[64];
extern int16_t zetas_withoutmont[64];

void ntt(int16_t *a);
void invntt(int16_t *a);
void basemul(int16_t *c, const int16_t *a, const int16_t *b, const int16_t zeta);

#endif