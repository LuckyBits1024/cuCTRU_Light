#ifndef CTRU_POLY_CUH
#define CTRU_POLY_CUH

#include <stdint.h>
#include "params.h"

typedef struct{
    int16_t coeffs[CTRU_N];
} poly;

__device__ void poly_reduce(poly &a);
__device__ void poly_freeze(poly &a);
__device__ void poly_fqcsubq(poly &a);
__device__ void poly_add(poly &c, const poly &a, const poly &b);
__device__ void poly_double(poly &b, const poly &a);
__device__ void poly_multi_p(poly &b, const poly &a);
__device__ void poly_ntt(int16_t regs[8], int16_t *s_ntt);
__device__ void poly_invntt(int16_t regs[8], int16_t *s_ntt);
__device__ void poly_basemul(poly &c, const poly &a, const poly &b);

#endif