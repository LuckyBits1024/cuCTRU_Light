#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "params.h"
#include "pack.h"
#include "poly.h"
#include "test_vector.h"
// gcc -Wall -O3 -march=native -mtune=native -fomit-frame-pointer -Wno-unknown-pragmas -I./src ./src/test_pack.c ./src/pack.c -o ./test_pack_cpu && ./test_pack_cpu
int main(){
    unsigned char *h_pk, *h_sk, *h_ct;
    poly a, r_pk, r_sk, r_ct;
    h_pk = (unsigned char *)malloc(CTRU_PKE_PUBLICKEYBYTES * sizeof(unsigned char));
    h_sk = (unsigned char *)malloc(CTRU_PKE_SECRETKEYBYTES_1 * sizeof(unsigned char));
    h_ct = (unsigned char *)malloc(CTRU_PKE_CIPHERTEXTBYTES * sizeof(unsigned char));
    for(int i = 0; i < CTRU_N; ++i){
        a.coeffs[i] = kFixedNTTInput[i];
    }
    // printf("生成的随机数组A为\n");
    // for(int i = 0; i < CTRU_N; ++i){
    //     printf("%d, ", a.coeffs[i]);
    //     if((i + 1) % 16 == 0) printf("\n");
    // }
    pack_pk(h_pk, &a);
    printf("pack_pk之后的数组为\n");
    for(int i = 0; i < CTRU_PKE_PUBLICKEYBYTES; ++i){
        printf("%d ", h_pk[i]);
        if((i + 1) % 32 == 0) printf("\n");
    }
    unpack_pk(&r_pk, h_pk);
    printf("unpack_pk之后的多项式为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d ", r_pk.coeffs[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }
    pack_sk(h_sk, &a);
    printf("pack_sk之后的数组为\n");
    for(int i = 0; i < CTRU_PKE_SECRETKEYBYTES_1; ++i){
        printf("%d ", h_sk[i]);
        if((i + 1) % 32 == 0) printf("\n");
    }
    unpack_sk(&r_sk, h_sk);
    printf("unpack_sk之后的多项式为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d ", r_sk.coeffs[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }
    pack_ct(h_ct, &a);
    printf("pack_ct之后的数组为\n");
    for(int i = 0; i < CTRU_PKE_CIPHERTEXTBYTES; ++i){
        printf("%d ", h_ct[i]);
        if((i + 1) % 32 == 0) printf("\n");
    }
    unpack_ct(&r_ct, h_ct);
    printf("unpack_ct之后的多项式为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d ", r_ct.coeffs[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }
    free(h_pk);
    free(h_sk);
    free(h_ct);
    return 0;
}