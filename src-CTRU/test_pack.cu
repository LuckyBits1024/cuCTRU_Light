#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include "params.h"
#include "pack.cuh"
#include "ntt.cuh"
#include "poly.cuh"
#include "test_vector.cuh"

//nvcc -std=c++17 -O2 -rdc=true -I. ./test_pack.cu ./pack.cu ./poly.cu ./ntt.cu ./reduce.cu ./randombytes.cu -lm -o ./test_pack_gpu && ./test_pack_gpu
__global__ void test_pack_pk(unsigned char *r, const poly &a){
    pack_pk(r, a);
}
__global__ void test_unpack_pk(int16_t *out, const unsigned char *a){
    int16_t regs[8];
    unpack_pk(regs, a);
    for(int i = 0; i < 8; ++i){
        out[8 * threadIdx.x + i] = regs[i];
    }
}
__global__ void test_pack_sk(unsigned char *r, const poly &a){
    pack_sk(r, a);
}
__global__ void test_unpack_sk(poly *r, const unsigned char *a){
    unpack_sk(*r, a);
}
__global__ void test_pack_ct(unsigned char *r, const poly &a){
    int16_t regs[8];
    for(int i = 0; i < 8; ++i){
        regs[i] = a.coeffs[8 * threadIdx.x + i];
    }
    pack_ct(r, regs);
}
__global__ void test_unpack_ct(int16_t *out, const unsigned char *a){
    unpack_ct(out, a);
}

int main(){
    unsigned char *h_pk, *h_sk, *h_ct;
    unsigned char *r_pk, *r_sk, *r_ct;
    int16_t *h_regs_pk, *d_regs_pk, *h_ct_out, *d_ct_out;
    poly *h_polyvec, *a_polyvec, *h_usk_out, *usk_out;
    size_t poly_size = sizeof(poly);
    h_pk = (unsigned char *)malloc(CTRU_PKE_PUBLICKEYBYTES * sizeof(unsigned char));
    h_sk = (unsigned char *)malloc(CTRU_PKE_SECRETKEYBYTES_1 * sizeof(unsigned char));
    h_ct = (unsigned char *)malloc(CTRU_PKE_CIPHERTEXTBYTES * sizeof(unsigned char));
    h_regs_pk = (int16_t *)malloc(CTRU_N * sizeof(int16_t));
    h_ct_out = (int16_t *)malloc(CTRU_N * sizeof(int16_t));
    h_polyvec = (poly *)malloc(poly_size);
    h_usk_out = (poly *)malloc(poly_size);
    cudaMalloc((void**)&r_pk, CTRU_PKE_PUBLICKEYBYTES * sizeof(unsigned char));
    cudaMalloc((void**)&r_sk, CTRU_PKE_SECRETKEYBYTES_1 * sizeof(unsigned char));
    cudaMalloc((void**)&r_ct, CTRU_PKE_CIPHERTEXTBYTES * sizeof(unsigned char));
    cudaMalloc((void**)&d_regs_pk, CTRU_N * sizeof(int16_t));
    cudaMalloc((void**)&d_ct_out, CTRU_N * sizeof(int16_t));
    cudaMalloc((void**)&a_polyvec, poly_size);
    cudaMalloc((void**)&usk_out, poly_size);

    srand(time(NULL));
    // 指定范围
    int min = 0;
    int max = CTRU_Q;
    // 使用固定数组
    for(int i = 0; i < CTRU_N; ++i){
        h_polyvec->coeffs[i] = kFixedNTTInput[i];
    }

    cudaMemcpy(a_polyvec, h_polyvec, poly_size, cudaMemcpyHostToDevice);
    
    test_pack_pk<<<1, 64>>>(r_pk, *a_polyvec);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_pk, r_pk, CTRU_PKE_PUBLICKEYBYTES * sizeof(unsigned char), cudaMemcpyDeviceToHost);
    printf("pack_pk之后的数组为\n");
    for(int i = 0; i < CTRU_PKE_PUBLICKEYBYTES; ++i){
        printf("%d ", h_pk[i]);
        if((i + 1) % 32 == 0) printf("\n");
    }

    test_unpack_pk<<<1, 64>>>(d_regs_pk, r_pk);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_regs_pk, d_regs_pk, CTRU_N * sizeof(int16_t), cudaMemcpyDeviceToHost);
    printf("unpack_pk之后的寄存器为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d ", h_regs_pk[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }

    test_pack_sk<<<1, 64>>>(r_sk, *a_polyvec);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_sk, r_sk, CTRU_PKE_SECRETKEYBYTES_1 * sizeof(unsigned char), cudaMemcpyDeviceToHost);
    printf("pack_sk之后的数组为\n");
    for(int i = 0; i < CTRU_PKE_SECRETKEYBYTES_1; ++i){
        printf("%d ", h_sk[i]);
        if((i + 1) % 32 == 0) printf("\n");
    }

    test_unpack_sk<<<1, 64>>>(usk_out, r_sk);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_usk_out, usk_out, poly_size, cudaMemcpyDeviceToHost);
    printf("unpack_sk之后的多项式为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d ", h_usk_out->coeffs[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }

    test_pack_ct<<<1, 64>>>(r_ct, *a_polyvec);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_ct, r_ct, CTRU_PKE_CIPHERTEXTBYTES * sizeof(unsigned char), cudaMemcpyDeviceToHost);
    printf("pack_ct之后的数组为\n");
    for(int i = 0; i < CTRU_PKE_CIPHERTEXTBYTES; ++i){
        printf("%d ", h_ct[i]);
        if((i + 1) % 32 == 0) printf("\n");
    }

    test_unpack_ct<<<1, 64>>>(d_ct_out, r_ct);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_ct_out, d_ct_out, CTRU_N * sizeof(int16_t), cudaMemcpyDeviceToHost);
    printf("unpack_ct之后的多项式为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d ", h_ct_out[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }

    free(h_pk);
    free(h_sk);
    free(h_ct);
    free(h_regs_pk);
    free(h_ct_out);
    free(h_polyvec);
    free(h_usk_out);
    cudaFree(a_polyvec);
    cudaFree(usk_out);
    cudaFree(r_pk);
    cudaFree(r_sk);
    cudaFree(r_ct);
    cudaFree(d_regs_pk);
    cudaFree(d_ct_out);
    return 0;
}