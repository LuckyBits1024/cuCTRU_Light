#include <cstddef>
#include <cstdint>
#include <stdint.h>
#include <vector_types.h>
#include "params.h"
#include "reduce.cuh"
#include "ntt.cuh"
#include "poly.cuh"
#include "coding.cuh"
#include "cbd.cuh"
#include "inv.cuh"

__device__ void poly_reduce(poly &a){
    a.coeffs[8 * threadIdx.x] = barrett_reduce(a.coeffs[8 * threadIdx.x]);
    a.coeffs[8 * threadIdx.x + 1] = barrett_reduce(a.coeffs[8 * threadIdx.x + 1]);
    a.coeffs[8 * threadIdx.x + 2] = barrett_reduce(a.coeffs[8 * threadIdx.x + 2]);
    a.coeffs[8 * threadIdx.x + 3] = barrett_reduce(a.coeffs[8 * threadIdx.x + 3]);
    a.coeffs[8 * threadIdx.x + 4] = barrett_reduce(a.coeffs[8 * threadIdx.x + 4]);
    a.coeffs[8 * threadIdx.x + 5] = barrett_reduce(a.coeffs[8 * threadIdx.x + 5]);
    a.coeffs[8 * threadIdx.x + 6] = barrett_reduce(a.coeffs[8 * threadIdx.x + 6]);
    a.coeffs[8 * threadIdx.x + 7] = barrett_reduce(a.coeffs[8 * threadIdx.x + 7]);
}

__device__ void poly_freeze(poly &a){
    poly_reduce(a);
    a.coeffs[8 * threadIdx.x] = fqcsubq(a.coeffs[8 * threadIdx.x]);
    a.coeffs[8 * threadIdx.x + 1] = fqcsubq(a.coeffs[8 * threadIdx.x + 1]);
    a.coeffs[8 * threadIdx.x + 2] = fqcsubq(a.coeffs[8 * threadIdx.x + 2]);
    a.coeffs[8 * threadIdx.x + 3] = fqcsubq(a.coeffs[8 * threadIdx.x + 3]);
    a.coeffs[8 * threadIdx.x + 4] = fqcsubq(a.coeffs[8 * threadIdx.x + 4]);
    a.coeffs[8 * threadIdx.x + 5] = fqcsubq(a.coeffs[8 * threadIdx.x + 5]);
    a.coeffs[8 * threadIdx.x + 6] = fqcsubq(a.coeffs[8 * threadIdx.x + 6]);
    a.coeffs[8 * threadIdx.x + 7] = fqcsubq(a.coeffs[8 * threadIdx.x + 7]);
}

__device__ void poly_add(poly &c, const poly &a, const poly &b){
    c.coeffs[8 * threadIdx.x] = a.coeffs[8 * threadIdx.x] + b.coeffs[8 * threadIdx.x];
    c.coeffs[8 * threadIdx.x + 1] = a.coeffs[8 * threadIdx.x + 1] + b.coeffs[8 * threadIdx.x + 1];
    c.coeffs[8 * threadIdx.x + 2] = a.coeffs[8 * threadIdx.x + 2] + b.coeffs[8 * threadIdx.x + 2];
    c.coeffs[8 * threadIdx.x + 3] = a.coeffs[8 * threadIdx.x + 3] + b.coeffs[8 * threadIdx.x + 3];
    c.coeffs[8 * threadIdx.x + 4] = a.coeffs[8 * threadIdx.x + 4] + b.coeffs[8 * threadIdx.x + 4];
    c.coeffs[8 * threadIdx.x + 5] = a.coeffs[8 * threadIdx.x + 5] + b.coeffs[8 * threadIdx.x + 5];
    c.coeffs[8 * threadIdx.x + 6] = a.coeffs[8 * threadIdx.x + 6] + b.coeffs[8 * threadIdx.x + 6];
    c.coeffs[8 * threadIdx.x + 7] = a.coeffs[8 * threadIdx.x + 7] + b.coeffs[8 * threadIdx.x + 7];
}

__device__ void poly_double(poly &b, const poly &a){
    b.coeffs[8 * threadIdx.x] = 2 * a.coeffs[8 * threadIdx.x];
    b.coeffs[8 * threadIdx.x + 1] = 2 * a.coeffs[8 * threadIdx.x + 1];
    b.coeffs[8 * threadIdx.x + 2] = 2 * a.coeffs[8 * threadIdx.x + 2];
    b.coeffs[8 * threadIdx.x + 3] = 2 * a.coeffs[8 * threadIdx.x + 3];
    b.coeffs[8 * threadIdx.x + 4] = 2 * a.coeffs[8 * threadIdx.x + 4];
    b.coeffs[8 * threadIdx.x + 5] = 2 * a.coeffs[8 * threadIdx.x + 5];
    b.coeffs[8 * threadIdx.x + 6] = 2 * a.coeffs[8 * threadIdx.x + 6];
    b.coeffs[8 * threadIdx.x + 7] = 2 * a.coeffs[8 * threadIdx.x + 7];
}

__device__ void poly_tomont(poly &a){
    const int16_t t = 867;
    a.coeffs[8 * threadIdx.x] = fqmul(a.coeffs[8 * threadIdx.x], t);
    a.coeffs[8 * threadIdx.x + 1] = fqmul(a.coeffs[8 * threadIdx.x + 1], t);
    a.coeffs[8 * threadIdx.x + 2] = fqmul(a.coeffs[8 * threadIdx.x + 2], t);
    a.coeffs[8 * threadIdx.x + 3] = fqmul(a.coeffs[8 * threadIdx.x + 3], t);
    a.coeffs[8 * threadIdx.x + 4] = fqmul(a.coeffs[8 * threadIdx.x + 4], t);
    a.coeffs[8 * threadIdx.x + 5] = fqmul(a.coeffs[8 * threadIdx.x + 5], t);
    a.coeffs[8 * threadIdx.x + 6] = fqmul(a.coeffs[8 * threadIdx.x + 6], t);
    a.coeffs[8 * threadIdx.x + 7] = fqmul(a.coeffs[8 * threadIdx.x + 7], t);
}

__device__ void poly_frommont(poly &a){
    a.coeffs[8 * threadIdx.x] = fqmul(a.coeffs[8 * threadIdx.x], 1);
    a.coeffs[8 * threadIdx.x + 1] = fqmul(a.coeffs[8 * threadIdx.x + 1], 1);
    a.coeffs[8 * threadIdx.x + 2] = fqmul(a.coeffs[8 * threadIdx.x + 2], 1);
    a.coeffs[8 * threadIdx.x + 3] = fqmul(a.coeffs[8 * threadIdx.x + 3], 1);
    a.coeffs[8 * threadIdx.x + 4] = fqmul(a.coeffs[8 * threadIdx.x + 4], 1);
    a.coeffs[8 * threadIdx.x + 5] = fqmul(a.coeffs[8 * threadIdx.x + 5], 1);
    a.coeffs[8 * threadIdx.x + 6] = fqmul(a.coeffs[8 * threadIdx.x + 6], 1);
    a.coeffs[8 * threadIdx.x + 7] = fqmul(a.coeffs[8 * threadIdx.x + 7], 1);
}

__device__ void poly_ntt(int16_t regs[8], int16_t *s_ntt){
    ntt(regs, s_ntt);
}

__device__ void poly_invntt(int16_t regs[8], int16_t *s_ntt){
    inv_ntt(regs, s_ntt);
}

__device__ void poly_sample_keygen(poly *a, const unsigned char *buf)
{
  cbd1(a, buf);
}

__device__ void poly_sample_enc(poly *a, const unsigned char *buf)
{
  cbd1(a, buf);
}

__device__ void poly_basemul(poly *c, const poly *a, const poly *b){
    if (threadIdx.x < CTRU_N / 16) {
        int16_t z = zetas_base[32 + threadIdx.x];
        int16_t regsA[8], regsB[8], regsC[8];
        for (int k = 0; k < 8; ++k) {
            regsA[k] = a->coeffs[threadIdx.x * 16 + k];
            regsB[k] = b->coeffs[threadIdx.x * 16 + k];
        }
        basemul(regsA, regsB, regsC, z);
        for (int k = 0; k < 8; ++k) {
            c->coeffs[threadIdx.x * 16 + k] = regsC[k];
        }
    } else if (threadIdx.x < 2 * (CTRU_N / 16)) {
        int16_t z = -zetas_base[32 + (threadIdx.x - (CTRU_N / 16))];
        int16_t regsA[8], regsB[8], regsC[8];
        for (int k = 0; k < 8; ++k) {
            regsA[k] = a->coeffs[(threadIdx.x - (CTRU_N / 16)) * 16 + 8 + k];
            regsB[k] = b->coeffs[(threadIdx.x - (CTRU_N / 16)) * 16 + 8 + k];
        }
        basemul(regsA, regsB, regsC, z);
        for (int k = 0; k < 8; ++k) {
            c->coeffs[(threadIdx.x - (CTRU_N / 16)) * 16 + 8 + k] = regsC[k];
        }
    }
}
