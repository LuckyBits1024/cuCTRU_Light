#ifndef CTRU_H
#define CTRU_H

#include "params.h"
#include <cstddef>

__global__ void ctru_keygen(uint8_t *pk, uint8_t *sk, uint8_t *coins, size_t keygen_mem_pool_pitch);
__global__ void ctru_enc(uint8_t *ct, const uint8_t *pk, const uint8_t *m, const uint8_t *coins, size_t encaps_mem_pool_pitch);
__global__ void ctru_dec(uint8_t *m, const uint8_t *ct, const uint8_t *sk, size_t decaps_mem_pool_pitch);
__global__ void ctru_ct2(const uint8_t *d_ct, const uint8_t *d_sk, uint8_t *d_ct2, int32_t *d_fail, size_t decaps_mem_pool_pitch);
__global__ void ctru_getk2(uint8_t *d_k, const int32_t *d_fail, const uint8_t *d_buf, const uint8_t *d_buf2, size_t decaps_mem_pool_pitch);

#endif //CTRU_CUH