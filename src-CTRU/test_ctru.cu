#include <cstdint>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include "params.h"
#include "ctru.cuh"
#include "test_vector.cuh"
// nvcc -std=c++17 -O2 -rdc=true   -I./src-CTRU   ./src-CTRU/test_ctru.cu   ./src-CTRU/ctru.cu   ./src-CTRU/ntt.cu   ./src-CTRU/reduce.cu   ./src-CTRU/randombytes.cu   ./src-CTRU/pack.cu   ./src-CTRU/poly.cu   ./src-CTRU/inv.cu   ./src-CTRU/coding.cu   ./src-CTRU/cbd.cu   ./src-CTRU/fips202.cu   -o ./test_ctru_gpu && ./test_ctru_gpu
// keygen
int main(){
    uint8_t *pk, *sk, *coins;
    uint8_t *d_pk, *d_sk, *d_coins;
    pk = (uint8_t *)malloc(CTRU_KEM_PUBLICKEYBYTES * sizeof(uint8_t));
    sk = (uint8_t *)malloc(CTRU_KEM_SECRETKEYBYTES * sizeof(uint8_t));
    coins = (uint8_t *)malloc(CTRU_N * sizeof(uint8_t));
    cudaMalloc((void**)&d_pk, CTRU_KEM_PUBLICKEYBYTES * sizeof(uint8_t));
    cudaMalloc((void**)&d_sk, CTRU_KEM_SECRETKEYBYTES * sizeof(uint8_t));
    cudaMalloc((void**)&d_coins, CTRU_N * sizeof(uint8_t));

    srand(time(NULL));
    // 指定范围
    int min = 0, max = 256;
    for(int i = 0; i < CTRU_N; ++i){
        coins[i] = kFixedNTTInput[i];
    }
    // for(int i = 0; i < CTRU_N; ++i){
    //     coins[i] = min + rand() % (max - min + 1);
    // }
    // printf("生成的随机数组coins为\n");
    // for(int i = 0; i < CTRU_N; ++i){
    //     printf("%d, ", coins[i]);
    //     if((i + 1) % 16 == 0) printf("\n");
    // }
    cudaMemcpy(d_coins, coins, CTRU_N * sizeof(uint8_t), cudaMemcpyHostToDevice);
    
    ctru_keygen<<<1, 64>>>(d_pk, d_sk, d_coins, 0);
    cudaDeviceSynchronize();
    
    cudaMemcpy(pk, d_pk, CTRU_KEM_PUBLICKEYBYTES * sizeof(uint8_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(sk, d_sk, CTRU_KEM_SECRETKEYBYTES * sizeof(uint8_t), cudaMemcpyDeviceToHost);
    printf("生成的私钥为\n");
    for(int i = 0; i < 512 + 1536 + 32; ++i){
        printf("%d ", sk[i]);
        if((i + 1) % 32 == 0) printf("\n");
    }
    printf("生成的公钥为\n");
    for(int i = 0; i < 1536; ++i){
        printf("%d ", pk[i]);
        if((i + 1) % 32 == 0) printf("\n");
    }

    // 清理资源
    free(pk);
    free(sk);
    free(coins);
    cudaFree(d_pk);
    cudaFree(d_sk);
    cudaFree(d_coins);
    return 0;
}