#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "reduce.cuh"
#include "randombytes.cuh"
#include "ntt.cuh"

__global__ void test_mul(int16_t *a, int16_t *b, int16_t *c){
    __shared__ int16_t sa_ntt[CTRU_N];
    __shared__ int16_t sb_ntt[CTRU_N];

    int16_t regsA[8];
    int16_t regsB[8];
    int16_t regsC[8];

    for(size_t i = 0; i < 8; ++i){
        regsA[i] = a[128 * i + threadIdx.x];
        regsB[i] = b[128 * i + threadIdx.x];
    }
    ntt(regsA, sa_ntt);
    ntt(regsB, sb_ntt);
    basemul(regsA, regsB, regsC, zetas1024_base[threadIdx.x]);
    inv_ntt(regsC, sa_ntt);
    for(size_t i = 0; i < 8; ++i){
        c[i * 128 + threadIdx.x] = regsC[i];
    }
}

int main(){
    int16_t *h_polyvec, *a_polyvec, *b_polyvec, *c_polyvec;
    size_t size = CTRU_N * sizeof(int16_t);
    h_polyvec = (int16_t *)malloc(size);
    cudaMalloc((void**)&a_polyvec, size);
    cudaMalloc((void**)&b_polyvec, size);
    cudaMalloc((void**)&c_polyvec, size);

    srand(time(NULL));
    // 指定范围
    int min = 0;
    int max = CTRU_Q;
    for(int i = 0; i < CTRU_N; ++i){
        h_polyvec[i] = min + rand() % (max - min + 1);
    }
    printf("生成的随机数组A为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d, ", h_polyvec[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }
    cudaMemcpy(a_polyvec, h_polyvec, size, cudaMemcpyHostToDevice);

    for(int i = 0; i < CTRU_N; ++i){
        h_polyvec[i] = min + rand() % (max - min + 1);
    }
    printf("生成的随机数组B为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d, ", h_polyvec[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }
    cudaMemcpy(b_polyvec, h_polyvec, size, cudaMemcpyHostToDevice);
    
    test_mul<<<1, 128>>>(a_polyvec, b_polyvec, c_polyvec);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_polyvec, c_polyvec, size, cudaMemcpyDeviceToHost);
    printf("A*B为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d ", h_polyvec[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }

    // 清理资源
    free(h_polyvec);
    cudaFree(a_polyvec);
    cudaFree(b_polyvec);
    cudaFree(c_polyvec);
    return 0;
}