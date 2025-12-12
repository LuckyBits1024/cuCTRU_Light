#ifndef CTRU_PACK_CUH
#define CTRU_PACK_CUH

#include <stdint.h>
#include "params.h"
#include "poly.cuh"

__device__ void pack_pk(unsigned char *r, const poly &a);
__device__ void unpack_pk(int16_t regs[8], const unsigned char *a);
__device__ void pack_sk(unsigned char *r, const poly &a);
__device__ void unpack_sk(poly &r, const unsigned char *a);
__device__ void pack_ct(unsigned char *r, const int16_t regs[8]);
__device__ void unpack_ct(int16_t r[CTRU_N], const unsigned char *a);
__device__ void pack_sk_f(unsigned char *r, const poly &a);
__device__ void unpack_sk_f(poly &r, const unsigned char *a);

#endif