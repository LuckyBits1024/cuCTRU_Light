#ifndef CTRU_REDUCE_H
#define CTRU_REDUCE_H

#include <stdint.h>
#include "params.h"

extern __device__ int16_t fq_inverse_table[CTRU_Q];


#define MONT 171   // 2^16 mod 769
#define QINV 64769   // 769^(-1) mod 2^16

// #define MONT 154    // 2^16 mod 641
// #define QINV 15745    //641^(-1) mod 2^16

// #define BARRETT_V 26174
#define BARRETT_V 21817

__device__ int16_t montgomery_reduce(int32_t a);
__device__ int16_t barrett_reduce(int16_t a);
__device__ int16_t fqcsubq(int16_t a);
__device__ int16_t fqmul(int16_t a, int16_t b);
__device__ int16_t fqinv(int16_t a);
__device__ int16_t fquniform();

#endif
