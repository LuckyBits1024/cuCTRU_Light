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

__device__ void pack_pk(unsigned char *r, const poly &a){
    r[10 * threadIdx.x + 0] = (a.coeffs[8 * threadIdx.x] >> 0);
    r[10 * threadIdx.x + 1] = (a.coeffs[8 * threadIdx.x] >> 8) | (a.coeffs[8 * threadIdx.x + 1] << 2);
    r[10 * threadIdx.x + 2] = (a.coeffs[8 * threadIdx.x + 1] >> 6) | (a.coeffs[8 * threadIdx.x + 2] << 4);
    r[10 * threadIdx.x + 3] = (a.coeffs[8 * threadIdx.x + 2] >> 4) | (a.coeffs[8 * threadIdx.x + 3] << 6);
    r[10 * threadIdx.x + 4] = (a.coeffs[8 * threadIdx.x + 3] >> 2);

    r[10 * threadIdx.x + 5] = (a.coeffs[8 * threadIdx.x + 4] >> 0);
    r[10 * threadIdx.x + 6] = (a.coeffs[8 * threadIdx.x + 4] >> 8) | (a.coeffs[8 * threadIdx.x + 5] << 2);
    r[10 * threadIdx.x + 7] = (a.coeffs[8 * threadIdx.x + 5] >> 6) | (a.coeffs[8 * threadIdx.x + 6] << 4);
    r[10 * threadIdx.x + 8] = (a.coeffs[8 * threadIdx.x + 6] >> 4) | (a.coeffs[8 * threadIdx.x + 7] << 6);
    r[10 * threadIdx.x + 9] = (a.coeffs[8 * threadIdx.x + 7] >> 2);
}

__device__ void unpack_pk(int16_t regs[8], const unsigned char *a){
    regs[0] = ((a[10 * threadIdx.x + 0] >> 0) | ((uint16_t)a[10 * threadIdx.x + 1] << 8)) & 0x3FF;
    regs[1] = ((a[10 * threadIdx.x + 1] >> 2) | ((uint16_t)a[10 * threadIdx.x + 2] << 6)) & 0x3FF;
    regs[2] = ((a[10 * threadIdx.x + 2] >> 4) | ((uint16_t)a[10 * threadIdx.x + 3] << 4)) & 0x3FF;
    regs[3] = ((a[10 * threadIdx.x + 3] >> 6) | ((uint16_t)a[10 * threadIdx.x + 4] << 2)) & 0x3FF;

    regs[4] = ((a[10 * threadIdx.x + 4] >> 0) | ((uint16_t)a[10 * threadIdx.x + 5] << 8)) & 0x3FF;
    regs[5] = ((a[10 * threadIdx.x + 5] >> 2) | ((uint16_t)a[10 * threadIdx.x + 6] << 6)) & 0x3FF;
    regs[6] = ((a[10 * threadIdx.x + 6] >> 4) | ((uint16_t)a[10 * threadIdx.x + 7] << 4)) & 0x3FF;
    regs[7] = ((a[10 * threadIdx.x + 7] >> 6) | ((uint16_t)a[10 * threadIdx.x + 8] << 2)) & 0x3FF;
}
__device__ void pack_sk(unsigned char *r, const poly &a){
    uint8_t t[8];
    t[0] = CTRU_BOUND - a.coeffs[8 * threadIdx.x + 0];
    t[1] = CTRU_BOUND - a.coeffs[8 * threadIdx.x + 1];
    t[2] = CTRU_BOUND - a.coeffs[8 * threadIdx.x + 2];
    t[3] = CTRU_BOUND - a.coeffs[8 * threadIdx.x + 3];
    t[4] = CTRU_BOUND - a.coeffs[8 * threadIdx.x + 4];
    t[5] = CTRU_BOUND - a.coeffs[8 * threadIdx.x + 5];
    t[6] = CTRU_BOUND - a.coeffs[8 * threadIdx.x + 6];
    t[7] = CTRU_BOUND - a.coeffs[8 * threadIdx.x + 7];

    r[3 * threadIdx.x + 0] = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
    r[3 * threadIdx.x + 1] = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
    r[3 * threadIdx.x + 2] = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
}

__device__ void unpack_sk(poly &r, const unsigned char *a){
    r.coeffs[8 * threadIdx.x + 0] = CTRU_BOUND - ((a[3 * threadIdx.x + 0] >> 0) & 0x7);
    r.coeffs[8 * threadIdx.x + 1] = CTRU_BOUND - ((a[3 * threadIdx.x + 0] >> 3) & 0x7);
    r.coeffs[8 * threadIdx.x + 2] = CTRU_BOUND - (((a[3 * threadIdx.x + 0] >> 6) | (a[3 * threadIdx.x + 1] << 2)) & 0x7);
    r.coeffs[8 * threadIdx.x + 3] = CTRU_BOUND - ((a[3 * threadIdx.x + 1] >> 1) & 0x7);
    r.coeffs[8 * threadIdx.x + 4] = CTRU_BOUND - ((a[3 * threadIdx.x + 1] >> 4) & 0x7);
    r.coeffs[8 * threadIdx.x + 5] = CTRU_BOUND - (((a[3 * threadIdx.x + 1] >> 7) | (a[3 * threadIdx.x + 2] << 1)) & 0x7);
    r.coeffs[8 * threadIdx.x + 6] = CTRU_BOUND - ((a[3 * threadIdx.x + 2] >> 2) & 0x7);
    r.coeffs[8 * threadIdx.x + 7] = CTRU_BOUND - ((a[3 * threadIdx.x + 2] >> 5) & 0x7);
}

__device__ void pack_ct(unsigned char *r, const int16_t regs[8]){

    r[8 * threadIdx.x + 0] = (uint8_t)regs[0];
    r[8 * threadIdx.x + 1] = (uint8_t)regs[1];
    r[8 * threadIdx.x + 2] = (uint8_t)regs[2];
    r[8 * threadIdx.x + 3] = (uint8_t)regs[3];
    r[8 * threadIdx.x + 4] = (uint8_t)regs[4];
    r[8 * threadIdx.x + 5] = (uint8_t)regs[5];
    r[8 * threadIdx.x + 6] = (uint8_t)regs[6];
    r[8 * threadIdx.x + 7] = (uint8_t)regs[7];
}

__device__ void unpack_ct(int16_t r[CTRU_N], const unsigned char *a){
    r[8 * threadIdx.x + 0] = (int16_t) a[8 * threadIdx.x + 0];
    r[8 * threadIdx.x + 1] = (int16_t) a[8 * threadIdx.x + 1];
    r[8 * threadIdx.x + 2] = (int16_t) a[8 * threadIdx.x + 2];
    r[8 * threadIdx.x + 3] = (int16_t) a[8 * threadIdx.x + 3];
    r[8 * threadIdx.x + 4] = (int16_t) a[8 * threadIdx.x + 4];
    r[8 * threadIdx.x + 5] = (int16_t) a[8 * threadIdx.x + 5];
    r[8 * threadIdx.x + 6] = (int16_t) a[8 * threadIdx.x + 6];
    r[8 * threadIdx.x + 7] = (int16_t) a[8 * threadIdx.x + 7];
}

//正确性待验证
__device__ void pack_sk_f(unsigned char *r, const poly &a){
    for(int i=threadIdx.x; i<CTRU_N/5; i+=64){
        unsigned char c = 0;
        c = (a.coeffs[5*i+4] + CTRU_ETA) & 255;
        c = (3*c + a.coeffs[5*i+3] + CTRU_ETA) & 255;
        c = (3*c + a.coeffs[5*i+2] + CTRU_ETA) & 255;
        c = (3*c + a.coeffs[5*i+1] + CTRU_ETA) & 255;
        c = (3*c + a.coeffs[5*i+0] + CTRU_ETA) & 255;
        r[i] = c;
    }
#if CTRU_N > (CTRU_N / 5) * 5
    if(threadIdx.x == 0){
        int i = CTRU_N / 5;
        unsigned char c = 0;
        for(int j = CTRU_N - (5*i) - 1; j>=0; j--) c = (3*c + a.coeffs[5*i+j] + CTRU_ETA) & 255;
        r[i] = c;
    }
#endif
}

__device__ void unpack_sk_f(poly &r, const unsigned char *a){
    for(int i=threadIdx.x; i<CTRU_N/5; i+=64){
        unsigned char c = a[i];
        r.coeffs[5*i+0] = (c%3 - CTRU_ETA) << 1; c/=3;
        r.coeffs[5*i+1] = (c%3 - CTRU_ETA) << 1; c/=3;
        r.coeffs[5*i+2] = (c%3 - CTRU_ETA) << 1; c/=3;
        r.coeffs[5*i+3] = (c%3 - CTRU_ETA) << 1; c/=3;
        r.coeffs[5*i+4] = (c%3 - CTRU_ETA) << 1;
    }
#if CTRU_N > (CTRU_N / 5) * 5
    if(threadIdx.x == 0){
        int i = CTRU_N/5; unsigned char c = a[i]; int j=0;
        while((5*i+j)<CTRU_N){ r.coeffs[5*i+j] = (c%3 - CTRU_ETA) << 1; c/=3; j++; }
        r.coeffs[0] += 1;
    }
#else
    if(threadIdx.x == 0){ r.coeffs[0] += 1; }
#endif
}