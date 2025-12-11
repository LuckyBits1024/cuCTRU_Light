#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include "params.h"
#include "pack.cuh"

__global__ void test_pack_pk(unsigned char *r, const poly &a){
    pack_pk(r, a);
}

int main(){
    unsigned char *h, *r;
    poly *h_polyvec, *a_polyvec;
    size_t size = sizeof(poly);
    h = (unsigned char *)malloc(1536 * sizeof(unsigned char));
    h_polyvec = (poly *)malloc(size);
    cudaMalloc((void**)&r, 1536 * sizeof(unsigned char));
    cudaMalloc((void**)&a_polyvec, size);

    srand(time(NULL));
    // 指定范围
    int min = 0;
    int max = CTRU_Q;
    for(int i = 0; i < CTRU_N; ++i){
        h_polyvec -> coeffs[i] = min + rand() % (max - min + 1);
    }
    printf("生成的随机数组A为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d, ", h_polyvec -> coeffs[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }
    cudaMemcpy(a_polyvec, h_polyvec, size, cudaMemcpyHostToDevice);
    
    test_pack_pk<<<1, 64>>>(r, *a_polyvec);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h, r, 1536 * sizeof(unsigned char), cudaMemcpyDeviceToHost);
    printf("pack_pk之后的数组为\n");
    for(int i = 0; i < CTRU_N * 12 / 8; ++i){
        printf("%d ", h[i]);
        if((i + 1) % 32 == 0) printf("\n");
    }

    // 清理资源
    free(h);
    free(h_polyvec);
    cudaFree(a_polyvec);
    cudaFree(r);
    return 0;
}