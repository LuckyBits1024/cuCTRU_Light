// #include "kem.cuh"
#include "poly.cuh"
#include "pack.cuh"
#include "reduce.cuh"
#include "cbd.cuh"
#include "ntt.cuh"
#include "coding.cuh"
#include "inv.cuh"
#include "fips202.cuh"
#include "params.h"

__global__ void ctru_keygen(uint8_t *g_pk, uint8_t *g_sk, uint8_t *g_coin, size_t keygen_mem_pool_pitch){
    int r = 0;
    int16_t regs[8];
    __shared__ int16_t s_ntt[CTRU_N];
    __shared__ int16_t t_ntt[CTRU_N];

    uint8_t *pk = g_pk + blockIdx.x * keygen_mem_pool_pitch;
    uint8_t *sk = g_sk + blockIdx.x * keygen_mem_pool_pitch;
    uint8_t *coins = g_coin + blockIdx.x * keygen_mem_pool_pitch;

    restart:

    uint64_t ctr = threadIdx.x;
    for(size_t j = threadIdx.x; j < (CTRU_SEEDBYTES + 7) / 8; j += blockDim.x){
        for(size_t i = 0; i < 8; ++i){
            coins[j * 8 + i] = ctr >> 8 * i;
        }
        ctr++;
    }
    gpu_keccak keccak;
    size_t out_blocks = (CTRU_N + (SHAKE128_RATE - 1)) / SHAKE128_RATE;
    keccak.notmp_shake<SHAKE128_RATE>(coins, out_blocks, coins, CTRU_SEEDBYTES);

    // poly_sample_keygen(&f, coins);
    uint8_t b0 = coins[2 * threadIdx.x + 0];
    uint8_t b1 = coins[2 * threadIdx.x + 1];
    regs[0] = 2 * (((b0 >> 0) & 1) - ((b0 >> 1) & 1));
    regs[1] = 2 * (((b0 >> 2) & 1) - ((b0 >> 3) & 1));
    regs[2] = 2 * (((b0 >> 4) & 1) - ((b0 >> 5) & 1));
    regs[3] = 2 * (((b0 >> 6) & 1) - ((b0 >> 7) & 1));
    regs[4] = 2 * (((b1 >> 0) & 1) - ((b1 >> 1) & 1));
    regs[5] = 2 * (((b1 >> 2) & 1) - ((b1 >> 3) & 1));
    regs[6] = 2 * (((b1 >> 4) & 1) - ((b1 >> 5) & 1));
    regs[7] = 2 * (((b1 >> 6) & 1) - ((b1 >> 7) & 1));

    // f.coeffs[0] += 1;
    if(threadIdx.x == 0){
        regs[0] += 1;
    }
    __syncthreads();

    // pack_sk(sk, &f); 另外一个pack_f是更加压缩的形式
    
    sk[3 * threadIdx.x + 0] = (regs[0] >> 0) | (regs[1] << 3) | (regs[2] << 6);
    sk[3 * threadIdx.x + 1] = (regs[2] >> 2) | (regs[3] << 1) | (regs[4] << 4) | (regs[5] << 7);
    sk[3 * threadIdx.x + 2] = (regs[5] >> 1) | (regs[6] << 2) | (regs[7] << 5);
    
    // regs2SM
    for(size_t i = 0; i < 8; ++i){
        s_ntt[8 * threadIdx.x + i] = regs[i];
    }
    // poly_ntt(&f);
    for(size_t i = 0; i < 8; ++i){
        regs[i] = s_ntt[64 * i + threadIdx.x];
    }
    poly_ntt(regs, s_ntt);
    for(size_t i = 0; i < 8; ++i){
        s_ntt[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8] = regs[i];
    }
    __syncthreads();
    // SM2regs
    for(size_t i = 0; i < 8; ++i){
        regs[i] = s_ntt[8 * threadIdx.x + i];
    }
    // r = poly_baseinv_opt(&finv, &f);
    int i = threadIdx.x >> 1;             // i = b/2，0..31
    int second = threadIdx.x & 1;         // 偶数块=0（前8个），奇数块=1（后8个）
    int16_t z = zetas_withoutmont[32 + i];
    if (second)
        z = -z;
    r += rq_inverse_opt(t_ntt + 8 * threadIdx.x, regs, z);

    if(r != 0){
        goto restart;
    }
    //之后再处理g

    // poly_sample_keygen(&g, coins + CTRU_COINBYTES_KEYGEN / 2);
    b0 = coins[CTRU_COINBYTES_KEYGEN / 2 + 2 * threadIdx.x + 0];
    b1 = coins[CTRU_COINBYTES_KEYGEN / 2 + 2 * threadIdx.x + 1];
    regs[0] = 2 * (((b0 >> 0) & 1) - ((b0 >> 1) & 1));
    regs[1] = 2 * (((b0 >> 2) & 1) - ((b0 >> 3) & 1));
    regs[2] = 2 * (((b0 >> 4) & 1) - ((b0 >> 5) & 1));
    regs[3] = 2 * (((b0 >> 6) & 1) - ((b0 >> 7) & 1));
    regs[4] = 2 * (((b1 >> 0) & 1) - ((b1 >> 1) & 1));
    regs[5] = 2 * (((b1 >> 2) & 1) - ((b1 >> 3) & 1));
    regs[6] = 2 * (((b1 >> 4) & 1) - ((b1 >> 5) & 1));
    regs[7] = 2 * (((b1 >> 6) & 1) - ((b1 >> 7) & 1));

    // regs2SM
    for(size_t i = 0; i < 8; ++i){
        s_ntt[8 * threadIdx.x + i] = regs[i];
    }
    // poly_ntt(&g);
    for(size_t i = 0; i < 8; ++i){
        regs[i] = s_ntt[64 * i + threadIdx.x];
    }
    poly_ntt(regs, s_ntt);
    for(size_t i = 0; i < 8; ++i){
        s_ntt[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8] = regs[i];
    }
    __syncthreads();

    // poly_basemul(&hhat, &g, &finv);
    base_mul(regs, t_ntt + 8 * threadIdx.x, zetas_withoutmont[threadIdx.x]);

    // regs2SM
    for(size_t i = 0; i < 8; ++i){
        s_ntt[8 * threadIdx.x + i] = regs[i];
    }
    // poly_invntt(&hhat);
    for(size_t i = 0; i < 8; ++i)
        regs[i] = s_ntt[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8];

    poly_invntt(regs, s_ntt);

    for(size_t i = 0; i < 8; ++i){
        s_ntt[64 * i + threadIdx.x] = regs[i];
    }
    
    __syncthreads();

    // poly_freeze(&hhat);
    for(size_t i = 0; i < 8; ++i){
        regs[i] = barrett_reduce(s_ntt[8 * threadIdx.x + i]);
        regs[i] = fqcsubq(regs[i]);
    }
    
    // pack_pk
    pk[10 * threadIdx.x + 0] = (regs[0] >> 0);
    pk[10 * threadIdx.x + 1] = (regs[0] >> 8) | (regs[1] << 2);
    pk[10 * threadIdx.x + 2] = (regs[1] >> 6) | (regs[2] << 4);
    pk[10 * threadIdx.x + 3] = (regs[2] >> 4) | (regs[3] << 6);
    pk[10 * threadIdx.x + 4] = (regs[3] >> 2);

    pk[10 * threadIdx.x + 5] = (regs[4] >> 0);
    pk[10 * threadIdx.x + 6] = (regs[4] >> 8) | (regs[5] << 2);
    pk[10 * threadIdx.x + 7] = (regs[5] >> 6) | (regs[6] << 4);
    pk[10 * threadIdx.x + 8] = (regs[6] >> 4) | (regs[7] << 6);
    pk[10 * threadIdx.x + 9] = (regs[7] >> 2);

    // return r;

}
__global__ void ctru_enc(uint8_t *d_ct, const uint8_t *d_pk, const uint8_t *d_m, const uint8_t *d_coins, size_t encaps_mem_pool_pitch){
    __shared__ int16_t s_ntt[CTRU_N];
    __shared__ int16_t r_tmp[CTRU_N];// 暂存r
    int16_t regs[8];

    uint8_t *ct = d_ct + blockIdx.x * encaps_mem_pool_pitch;
    auto pk = d_pk + blockIdx.x * encaps_mem_pool_pitch;
    auto m = d_m + blockIdx.x * encaps_mem_pool_pitch;
    auto coins = d_coins + blockIdx.x * encaps_mem_pool_pitch;

    // poly_sample_enc(&r, coins);
    uint8_t b0 = coins[2 * threadIdx.x + 0];
    uint8_t b1 = coins[2 * threadIdx.x + 1];

    regs[0] = 2 * (((b0 >> 0) & 1) - ((b0 >> 1) & 1));
    regs[1] = 2 * (((b0 >> 2) & 1) - ((b0 >> 3) & 1));
    regs[2] = 2 * (((b0 >> 4) & 1) - ((b0 >> 5) & 1));
    regs[3] = 2 * (((b0 >> 6) & 1) - ((b0 >> 7) & 1));
    regs[4] = 2 * (((b1 >> 0) & 1) - ((b1 >> 1) & 1));
    regs[5] = 2 * (((b1 >> 2) & 1) - ((b1 >> 3) & 1));
    regs[6] = 2 * (((b1 >> 4) & 1) - ((b1 >> 5) & 1));
    regs[7] = 2 * (((b1 >> 6) & 1) - ((b1 >> 7) & 1));

    // regs2SM
    for(size_t i = 0; i < 8; ++i){
        r_tmp[8 * threadIdx.x + i] = regs[i];
    }
    // poly_ntt(&r);
    for(size_t i = 0; i < 8; ++i){
        regs[i] = r_tmp[64 * i + threadIdx.x];
    }
    poly_ntt(regs, s_ntt);
    for(size_t i = 0; i < 8; ++i){
        r_tmp[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8] = regs[i];
    }
    __syncthreads();

    // unpack_pk(&hhat, pk);
    unpack_pk(regs, pk);

    // regs2SM
    for(size_t i = 0; i < 8; ++i){
        s_ntt[8 * threadIdx.x + i] = regs[i];
    }
    // poly_ntt(&hhat);
    for(size_t i = 0; i < 8; ++i){
        regs[i] = s_ntt[64 * i + threadIdx.x];
    }
    poly_ntt(regs, s_ntt);
    for(size_t i = 0; i < 8; ++i){
        s_ntt[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8] = regs[i];
    }
    __syncthreads();

    // regs2SM
    for(size_t i = 0; i < 8; ++i){
        s_ntt[8 * threadIdx.x + i] = regs[i];
    }

    // poly_basemul(&sigma, &hhat, &r);
    base_mul(regs, r_tmp + 8 * threadIdx.x, zetas1024_base[threadIdx.x]);

    // regs2SM
    for(size_t i = 0; i < 8; ++i){
        s_ntt[8 * threadIdx.x + i] = regs[i];
    }
    // poly_invntt(&sigma);
    for(size_t i = 0; i < 8; ++i)
        regs[i] = s_ntt[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8];

    poly_invntt(regs, s_ntt);

    for(size_t i = 0; i < 8; ++i){
        s_ntt[64 * i + threadIdx.x] = regs[i];
    }
    __syncthreads();

    // poly_invntt(&sigma);
    for(size_t i = 0; i < 8; ++i){
        regs[i] = barrett_reduce(s_ntt[8 * threadIdx.x + i]);
        regs[i] = fqcsubq(regs[i]);
    }

    // poly_encode_compress(&c, &sigma, m);
    poly_encode_compress(regs, m);

    // pack_ct(ct, &c)
    pack_ct(ct, regs);
}