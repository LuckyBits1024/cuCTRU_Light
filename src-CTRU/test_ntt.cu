// nvcc -std=c++17 -O2 -rdc=true -I./src-CTRU ./src-CTRU/test_ntt.cu ./src-CTRU/ntt.cu ./src-CTRU/reduce.cu ./src-CTRU/randombytes.cu -lm -o ./test_ntt_gpu && ./test_ntt_gpu
#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "reduce.cuh"
#include "randombytes.cuh"
#include "ntt.cuh"
#include "test_vector.cuh"

__global__ void test_ntt(int16_t *a_polyvec, int16_t *b_polyvec){
    __shared__ int16_t sa_ntt[CTRU_N];
    __shared__ int16_t sb_ntt[CTRU_N];

    int16_t regsA[8];
    int16_t regsB[8];

    for (size_t i = 0; i < 8; ++i){
        regsA[i] = a_polyvec[64 * i + threadIdx.x];
        regsB[i] = b_polyvec[64 * i + threadIdx.x];
        // printf("线程 %d  寄存器的值为 %d \n",threadIdx.x,regs[i]);
    }

    //写回到global mem 
    ntt(regsA, sa_ntt);
    ntt(regsB, sb_ntt);
    for (size_t i = 0; i < 8; ++i){
        a_polyvec[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8] = regsA[i];
        b_polyvec[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8] = regsB[i];
    }

}

__global__ void test_inv_ntt(int16_t *g_polyvec){
    __shared__ int16_t s_ntt[CTRU_N];
    int16_t regs[8];

    for(size_t i = 0; i < 8; ++i)
        regs[i] = g_polyvec[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8];

    inv_ntt(regs, s_ntt);

    for(size_t i = 0; i < 8; ++i){
        g_polyvec[64 * i + threadIdx.x] = regs[i];
    }
}
__global__ void test_base_mul(int16_t *a, int16_t *b, int16_t *c, int16_t zeta){
    basemul(a, b, c, zeta);
}


//test_ntt/test_ntt7681
int main1(){
    int16_t *h_polyvec, *a_polyvec, *b_polyvec;
    size_t size = CTRU_N * sizeof(int16_t);
    h_polyvec = (int16_t *)malloc(size);
    cudaMalloc((void**)&a_polyvec, size);
    cudaMalloc((void**)&b_polyvec, size);

    printf("使用固定数组A\n");
    for(int i = 0; i < CTRU_N; ++i){
        h_polyvec[i] = kFixedNTTInput[i];
        printf("%d, ", h_polyvec[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }
    cudaMemcpy(a_polyvec, h_polyvec, size, cudaMemcpyHostToDevice);

    printf("使用固定数组B（与A相同）\n");
    for(int i = 0; i < CTRU_N; ++i){
        h_polyvec[i] = kFixedNTTInput[i];
        printf("%d, ", h_polyvec[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }
    cudaMemcpy(b_polyvec, h_polyvec, size, cudaMemcpyHostToDevice);
    
    test_ntt<<<1, 64>>>(a_polyvec, b_polyvec);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_polyvec, a_polyvec, size, cudaMemcpyDeviceToHost);
    printf("NTT之后的数组A为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d ", h_polyvec[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }

    cudaMemcpy(h_polyvec, b_polyvec, size, cudaMemcpyDeviceToHost);
    printf("NTT之后的数组B为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d ", h_polyvec[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }

    // 清理资源
    free(h_polyvec);
    cudaFree(a_polyvec);
    cudaFree(b_polyvec);
    return 0;
}

//test_inv_ntt/test_inv_ntt7681
int main2(){
    int16_t *h_polyvec, *a_polyvec;
    size_t size = CTRU_N * sizeof(int16_t);
    h_polyvec = (int16_t *)malloc(size);
    cudaMalloc((void**)&a_polyvec, size);

    printf("使用固定数组作为逆NTT输入\n");
    for(int i = 0; i < CTRU_N; ++i){
        h_polyvec[i] = kFixedInvNTTInput[i];
        printf("%d, ", h_polyvec[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }
    cudaMemcpy(a_polyvec, h_polyvec, size, cudaMemcpyHostToDevice);
    
    test_inv_ntt<<<1, 64>>>(a_polyvec);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_polyvec, a_polyvec, size, cudaMemcpyDeviceToHost);
    printf("INV_NTT之后的数组A为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d ", h_polyvec[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }

    // 清理资源
    free(h_polyvec);
    cudaFree(a_polyvec);
    return 0;
}

//test_base_mul
int main3(){
    int16_t n = 8;
    int16_t *h_polyvec, *a_polyvec, *b_polyvec, *c_polyvec;
    size_t size = n * sizeof(int16_t);
    h_polyvec = (int16_t *)malloc(size);
    cudaMalloc((void**)&a_polyvec, size);
    cudaMalloc((void**)&b_polyvec, size);
    cudaMalloc((void**)&c_polyvec, size);

    srand(time(NULL));
    // 指定范围
    int min = 0;
    int max = CTRU_Q;
    for(int i = 0; i < n; ++i) h_polyvec[i] = min + rand() % (max - min + 1);
    printf("生成的随机数组A为\n");
    for(int i = 0; i < n; ++i){
        printf("%d, ", h_polyvec[i]);
    }
    printf("\n");
    cudaMemcpy(a_polyvec, h_polyvec, size, cudaMemcpyHostToDevice);

    for(int i = 0; i < n; ++i) h_polyvec[i] = min + rand() % (max - min + 1);
    printf("生成的随机数组B为\n");
    for(int i = 0; i < n; ++i){
        printf("%d, ", h_polyvec[i]);
    }
    printf("\n");
    cudaMemcpy(b_polyvec, h_polyvec, size, cudaMemcpyHostToDevice);

    test_base_mul<<<1, 32>>>(a_polyvec, b_polyvec, c_polyvec, 1117);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_polyvec, c_polyvec, size, cudaMemcpyDeviceToHost);
    printf("base_mul之后的数组为\n");
    for(int i = 0; i < n; ++i){
        printf("%d ", h_polyvec[i]);
    }
    printf("\n");

    return 0;
}

int main(){
    // main1();
    // main2();
    main3();
}