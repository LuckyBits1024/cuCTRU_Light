#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "params.h"
#include "reduce.cuh"
#include "ntt.cuh"

//test
// 原始域 没有乘R
__device__ int16_t zetas_base[64] = {
    154, 640, 318, 256, 625, 100, 250, 40, 631, 383, 637, 25, 359, 160, 400, 64,
    609, 200, 500, 80, 636, 512, 639, 333, 77, 320, 159, 128, 633, 50, 125, 20,
    29, 620, 268, 248, 305, 177, 122, 199, 431, 351, 557, 525, 488, 155, 67, 62,
    610, 354, 244, 398, 536, 496, 599, 583, 335, 310, 134, 124, 473, 409, 61, 420};

__device__ int16_t zetas_inv_base[64] = {
    -420, -61, -409, -473, -124, -134, -310, -335, -583, -599, -496, -536, -398, -244, -354, -610,
    -62, -67, -155, -488, -525, -557, -351, -431, -199, -122, -177, -305, -248, -268, -620, -29,
    -20, -125, -50, -633, -128, -159, -320, -77, -333, -639, -512, -636, -80, -500, -200, -609,
    -64, -400, -160, -359, -25, -637, -383, -631, -40, -250, -100, -625, -256, -318, -640, 10};

__device__ int16_t zetas[64] = {
    171, 605, 81, 408, 766, 583, 120, 519, 753, 546, 640, 461, 432, 638, 407, 626,
    430, 514, 487, 203, 694, 733, 693, 671, 369, 577, 620, 759, 34, 570, 178, 270,
    86, 718, 405, 502, 754, 608, 600, 288, 689, 423, 124, 767, 622, 114, 497, 54,
    612, 263, 128, 246, 394, 589, 389, 279, 307, 578, 24, 719, 170, 543, 121, 581};

__device__ int16_t zetas_inv[64] = {
    -581, -121, -543, -170, -719, -24, -578, -307, -279, -389, -589, -394, -246, -128, -263, -612,
    -54, -497, -114, -622, -767, -124, -423, -689, -288, -600, -608, -754, -502, -405, -718, -86,
    -270, -178, -570, -34, -759, -620, -577, -369, -671, -693, -733, -694, -203, -487, -514, -430,
    -626, -407, -638, -432, -461, -640, -546, -753, -519, -120, -583, -766, -408, -81, -605, 541};  // 64^(-1) mod 769

// CT

__device__ void ntt_butt(int16_t &a, int16_t &b, int16_t zeta){
    int16_t t = fqmul(zeta, b);
    b = barrett_reduce(a - t);
    a = barrett_reduce(a + t);
}

__device__ void inv_ntt_butt(int16_t &a, int16_t &b, int16_t zeta){
    int16_t t = a;
    a = barrett_reduce(t + b);
    b = barrett_reduce(t - b);
    b = fqmul(zeta, b);
}

__device__ void ntt(int16_t regs[8], int16_t *s_ntt){

    // level 1 取下标系数为256的两个数进行蝴蝶变换
    ntt_butt(regs[0], regs[4], zetas[1]);
    ntt_butt(regs[1], regs[5], zetas[1]);
    ntt_butt(regs[2], regs[6], zetas[1]);
    ntt_butt(regs[3], regs[7], zetas[1]);

    // level 2  128
    ntt_butt(regs[0], regs[2], zetas[2]);
    ntt_butt(regs[1], regs[3], zetas[2]);
    ntt_butt(regs[4], regs[6], zetas[3]);
    ntt_butt(regs[5], regs[7], zetas[3]);

    // level 3   64
    ntt_butt(regs[0], regs[1], zetas[4]);
    ntt_butt(regs[2], regs[3], zetas[5]);
    ntt_butt(regs[4], regs[5], zetas[6]);
    ntt_butt(regs[6], regs[7], zetas[7]);

    // SMEM exchange
    #pragma unroll
    for (size_t i = 0; i < 8; i++)
        s_ntt[i * 64 + threadIdx.x] = regs[i]; //以64为间隔进行摆放
    __syncwarp();

    #pragma unroll
    for (size_t i = 0; i < 8; i++)
    //每8个线程对应64个系数，8个线程组成一个线程快。因此threadIdx.x/8用于确定位于哪一个数据块中
    // threadIdx.x & 7定位该线程在线程快内的偏移量，i*4定位是该线程中的哪一个数
        regs[i] = s_ntt[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8]; 
    __syncwarp();

    // level 4
    ntt_butt(regs[0], regs[4], zetas[8 + threadIdx.x / 8]);
    ntt_butt(regs[1], regs[5], zetas[8 + threadIdx.x / 8]);
    ntt_butt(regs[2], regs[6], zetas[8 + threadIdx.x / 8]);
    ntt_butt(regs[3], regs[7], zetas[8 + threadIdx.x / 8]);

    // level 5
    ntt_butt(regs[0], regs[2], zetas[16 + (threadIdx.x / 8) * 2]);
    ntt_butt(regs[1], regs[3], zetas[16 + (threadIdx.x / 8) * 2]);
    ntt_butt(regs[4], regs[6], zetas[17 + (threadIdx.x / 8) * 2]);
    ntt_butt(regs[5], regs[7], zetas[17 + (threadIdx.x / 8) * 2]);

    // level 6
    ntt_butt(regs[0], regs[1], zetas[32 + (threadIdx.x / 8) * 4]);
    ntt_butt(regs[2], regs[3], zetas[33 + (threadIdx.x / 8) * 4]);
    ntt_butt(regs[4], regs[5], zetas[34 + (threadIdx.x / 8) * 4]);
    ntt_butt(regs[6], regs[7], zetas[35 + (threadIdx.x / 8) * 4]);

// #pragma unroll
//     for (size_t i = 0; i < 8; i++)
//         //怎么读出来的，就怎么写回去
//         s_ntt[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8] = regs[i];
//     __syncwarp();
// #pragma unroll
//     for (size_t i = 0; i < 8; i++)
//         regs[i] = s_ntt[threadIdx.x * 8 + i]; //每个线程取8个系数即可

}


//invntt
__device__ void inv_ntt(int16_t regs[8], int16_t *s_ntt){
    // level 6
    inv_ntt_butt(regs[0], regs[1], zetas_inv[32 + (threadIdx.x / 8) * 4]);
    inv_ntt_butt(regs[2], regs[3], zetas_inv[33 + (threadIdx.x / 8) * 4]);
    inv_ntt_butt(regs[4], regs[5], zetas_inv[34 + (threadIdx.x / 8) * 4]);
    inv_ntt_butt(regs[6], regs[7], zetas_inv[35 + (threadIdx.x / 8) * 4]);

    // level 5
    inv_ntt_butt(regs[0], regs[2], zetas_inv[16 + (threadIdx.x / 8) * 2]);
    inv_ntt_butt(regs[1], regs[3], zetas_inv[16 + (threadIdx.x / 8) * 2]);
    inv_ntt_butt(regs[4], regs[6], zetas_inv[17 + (threadIdx.x / 8) * 2]);
    inv_ntt_butt(regs[5], regs[7], zetas_inv[17 + (threadIdx.x / 8) * 2]);

    // level 4
    inv_ntt_butt(regs[0], regs[4], zetas_inv[8 + threadIdx.x / 8]);
    inv_ntt_butt(regs[1], regs[5], zetas_inv[8 + threadIdx.x / 8]);
    inv_ntt_butt(regs[2], regs[6], zetas_inv[8 + threadIdx.x / 8]);
    inv_ntt_butt(regs[3], regs[7], zetas_inv[8 + threadIdx.x / 8]);
    #pragma unroll
    for (size_t i = 0; i < 8; i++)
    //每8个线程对应64个系数，8个线程组成一个线程快。因此threadIdx.x/8用于确定位于哪一个数据块中
    // threadIdx.x & 7定位该线程在线程快内的偏移量，i*4定位是该线程中的哪一个数
        s_ntt[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8] = regs[i]; 
    __syncwarp();

    // SMEM exchange
    #pragma unroll
    for (size_t i = 0; i < 8; i++)
         regs[i] = s_ntt[i * 64 + threadIdx.x]; //以64为间隔进行摆放
    __syncwarp();

    // level 3   64
    inv_ntt_butt(regs[0], regs[1], zetas_inv[4]);
    inv_ntt_butt(regs[2], regs[3], zetas_inv[5]);
    inv_ntt_butt(regs[4], regs[5], zetas_inv[6]);
    inv_ntt_butt(regs[6], regs[7], zetas_inv[7]);

    // level 2  128
    inv_ntt_butt(regs[0], regs[2], zetas_inv[2]);
    inv_ntt_butt(regs[1], regs[3], zetas_inv[2]);
    inv_ntt_butt(regs[4], regs[6], zetas_inv[3]);
    inv_ntt_butt(regs[5], regs[7], zetas_inv[3]);

    // level 1 取下标系数为256的两个数进行蝴蝶变换
    inv_ntt_butt(regs[0], regs[4], zetas_inv[1]);
    inv_ntt_butt(regs[1], regs[5], zetas_inv[1]);
    inv_ntt_butt(regs[2], regs[6], zetas_inv[1]);
    inv_ntt_butt(regs[3], regs[7], zetas_inv[1]);

    // 只是测试的时候用 把逆ntt的结果映射正常域
    // #pragma unroll
    // for (size_t i = 0; i < 8; i++)
    //     regs[i] = fqmul(regs[i], 255);
    // #pragma unroll
    // for (size_t i = 0; i < 8; i++)
    //     regs[i] = fqcsubq(regs[i]);

}
__device__ int16_t CALC_D(int16_t a[8], int16_t b[8], int16_t x, int16_t y, int16_t d[8]){
    return fqmul((a[x] + a[y]), (b[x] + b[y])) - d[x] - d[y];
}

__device__ void basemul(int16_t a[8], int16_t b[8], int16_t c[8], int16_t zeta){
    int16_t d[8];
    for(size_t i = 0; i < 8; ++i) d[i] = fqmul(a[i], b[i]);
    c[0] = d[0] + fqmul((CALC_D(a, b, 1, 7, d) + CALC_D(a, b, 2, 6, d) + CALC_D(a, b, 3, 5, d) + d[4]), zeta);
    c[1] = barrett_reduce(CALC_D(a, b, 0, 1, d) + fqmul((CALC_D(a, b, 2, 7, d) + CALC_D(a, b, 3, 6, d) + CALC_D(a, b, 4, 5, d)), zeta));
    c[2] = barrett_reduce(CALC_D(a, b, 0, 2, d) + d[1] + fqmul((CALC_D(a, b, 3, 7, d) + CALC_D(a, b, 4, 6, d) + d[5]), zeta));
    c[3] = barrett_reduce(CALC_D(a, b, 0, 3, d) + CALC_D(a, b, 1, 2, d) + fqmul((CALC_D(a, b, 4, 7, d) + CALC_D(a, b, 5, 6, d)), zeta));
    c[4] = barrett_reduce(CALC_D(a, b, 0, 4, d) + CALC_D(a, b, 1, 3, d) + d[2] + fqmul((CALC_D(a, b, 5, 7, d) + d[6]), zeta));  
    c[5] = barrett_reduce(CALC_D(a, b, 0, 5, d) + CALC_D(a, b, 1, 4, d) + CALC_D(a, b, 2, 3, d)) + fqmul(CALC_D(a, b, 6, 7, d), zeta);
    c[6] = barrett_reduce(CALC_D(a, b, 0, 6, d) + CALC_D(a, b, 1, 5, d) + CALC_D(a, b, 2, 4, d)) + d[3] + fqmul(d[7], zeta);
    c[7] = barrett_reduce(CALC_D(a, b, 0, 7, d) + CALC_D(a, b, 1, 6, d)) + barrett_reduce(CALC_D(a, b, 2, 5, d) + CALC_D(a, b, 3, 4, d));
}

__device__ void base_mul(int16_t a[8], int16_t b[8], int16_t zeta){
    int16_t c[8], d[8];
    for(size_t i = 0; i < 8; ++i) d[i] = fqmul(a[i], b[i]);
    c[0] = d[0] + fqmul((CALC_D(a, b, 1, 7, d) + CALC_D(a, b, 2, 6, d) + CALC_D(a, b, 3, 5, d) + d[4]), zeta);
    c[1] = barrett_reduce(CALC_D(a, b, 0, 1, d) + fqmul((CALC_D(a, b, 2, 7, d) + CALC_D(a, b, 3, 6, d) + CALC_D(a, b, 4, 5, d)), zeta));
    c[2] = barrett_reduce(CALC_D(a, b, 0, 2, d) + d[1] + fqmul((CALC_D(a, b, 3, 7, d) + CALC_D(a, b, 4, 6, d) + d[5]), zeta));
    c[3] = barrett_reduce(CALC_D(a, b, 0, 3, d) + CALC_D(a, b, 1, 2, d) + fqmul((CALC_D(a, b, 4, 7, d) + CALC_D(a, b, 5, 6, d)), zeta));
    c[4] = barrett_reduce(CALC_D(a, b, 0, 4, d) + CALC_D(a, b, 1, 3, d) + d[2] + fqmul((CALC_D(a, b, 5, 7, d) + d[6]), zeta));  
    c[5] = barrett_reduce(CALC_D(a, b, 0, 5, d) + CALC_D(a, b, 1, 4, d) + CALC_D(a, b, 2, 3, d)) + fqmul(CALC_D(a, b, 6, 7, d), zeta);
    c[6] = barrett_reduce(CALC_D(a, b, 0, 6, d) + CALC_D(a, b, 1, 5, d) + CALC_D(a, b, 2, 4, d)) + d[3] + fqmul(d[7], zeta);
    c[7] = barrett_reduce(CALC_D(a, b, 0, 7, d) + CALC_D(a, b, 1, 6, d)) + barrett_reduce(CALC_D(a, b, 2, 5, d) + CALC_D(a, b, 3, 4, d));
    for(size_t i = 0; i < 8; ++i) a[i] = c[i];
}


