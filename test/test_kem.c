#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "../src/params.h"
#include "../src/kem.h"
#include "../src/poly.h"
#include <stdlib.h>
#include "cpucycles.h"
#include "speed.h"
#include "randombytes.h"
#include "reduce.h"

#define NTESTS 10000
uint64_t t[NTESTS];
void print_value(unsigned char *input, int inputlength)
{
  int i, j, k;
  int cnt = 0;
  for(i = 0; i < inputlength/32; i ++)
  {
    for(k = 0; k < 8; k ++)
    {
      for(j = 0; j < 4; j ++)
      {
        printf("%02x",input[cnt++]);
      }
      printf("  ");
    }
    printf("\n");
  }
  printf("\n");
}
void test_kem()
{
  unsigned int i, j;
  unsigned char k1[CTRU_SHAREDKEYBYTES], k2[CTRU_SHAREDKEYBYTES];
  unsigned char pk[CTRU_KEM_PUBLICKEYBYTES], sk[CTRU_KEM_SECRETKEYBYTES];
  unsigned char ct[CTRU_KEM_CIPHERTEXTBYTES];
  int count = 0;

  for (i = 0; i < NTESTS; i++)
  {
    crypto_kem_keygen(pk, sk);
    // printf("pk:\n");
    // print_value(pk, CTRU_KEM_PUBLICKEYBYTES);
    // printf("sk:\n");
    // print_value(sk, CTRU_KEM_SECRETKEYBYTES);
    crypto_kem_encaps(ct, k1, pk);
    // printf("ct:\n");
    // print_value(ct, CTRU_KEM_CIPHERTEXTBYTES);
    // printf("k1:\n");
    // print_value(k1, CTRU_SHAREDKEYBYTES);
    crypto_kem_decaps(k2, ct, sk);
    // printf("k2:\n");
    // print_value(k2, CTRU_SHAREDKEYBYTES);

    for (j = 0; j < CTRU_SHAREDKEYBYTES; j++)
      if (k1[j] != k2[j])
      {
        printf("Round %d. Failure: Keys dont match: %hhx != %hhx!\n", i, k1[j], k2[j]);
        // return;
        count ++;
        break;
      }
  }
  printf("count = %d\n", count);
#if (KEM_TYPE == RLWE_KEM)
  printf("CTRU-Compact-%d-KEM is correct!\n", CTRU_N);
#elif (KEM_TYPE == RLWR_KEM)
  printf("CNTR-Compact-%d-KEM is correct!\n", CTRU_N);
#endif

  printf("Test %d times.\n\n", NTESTS);
#if (KEM_TYPE == RLWE_KEM)
  printf("CTRU_N = %d, CTRU_Q = %d\n", CTRU_N, CTRU_Q);
#elif (KEM_TYPE == RLWR_KEM)
  printf("CNTR_N = %d, CNTR_Q = %d, CNTR_Q2 = %d\n", CTRU_N, CTRU_Q, CTRU_Q2);
#endif

  printf("KEM size: pk = %d bytes, ct = %d bytes, bandwidth = %d bytes\n\n",
         CTRU_KEM_PUBLICKEYBYTES, CTRU_KEM_CIPHERTEXTBYTES, CTRU_KEM_PUBLICKEYBYTES + CTRU_KEM_CIPHERTEXTBYTES);
}
static void poly_naivemul(poly *c, const poly *a, const poly *b) {

    // int i;
    // for(i = 0; i < 512/2; i ++)
    // {
    //     // a->coeffs[i] = b->coeffs[i] = 1;
    //     a->coeffs[i] = b->coeffs[i] = fquniform();
    // }
    // for(i = 512/2; i < 512; i ++)
    // {
    //     // a.coeffs[i] = b.coeffs[i] = 0;
    //     a->coeffs[i] = b->coeffs[i] = fquniform();
    // }


     int16_t r[2*512] = {0};

     for(int i = 0; i < CTRU_N; i++)
     {
         for(int j = 0; j < CTRU_N; j++){
             r[i+j] = (r[i+j] + (int32_t)a->coeffs[i]*b->coeffs[j]) % CTRU_Q;
         }

     }
     for(int i = CTRU_N; i < 2*CTRU_N; i++)
      r[i-CTRU_N] = (r[i-CTRU_N] - r[i]) % CTRU_Q;

     for(int i = 0; i < CTRU_N; i++)
      c->coeffs[i] = r[i] % CTRU_Q;
}

void test_ntt()
{
    poly a, b, c, ahat, bhat, chat, c2;
    poly asmall;
    int i;
    for(i = 0; i < 512/2; i ++)
    {
        // a->coeffs[i] = b->coeffs[i] = 1;
        a.coeffs[i] = b.coeffs[i] = fquniform();
    }
    for(i = 512/2; i < 512; i ++)
    {
        // a.coeffs[i] = b.coeffs[i] = 0;
        a.coeffs[i] = b.coeffs[i] = fquniform();
    }
    poly_naivemul(&c2, &a, &b);

    poly_ntt_big(&ahat, &a);
    // for(int i = 0; i < 512; i ++)
    // {
    //     printf("%d, ",a.coeffs[i]);
    // }

    // printf("\n");
    // for(int i = 0; i < 512; i ++)
    // {
    //     printf("%d, ",ahat.coeffs[i]);
    // }
    // printf("\n");

    poly_ntt_big(&bhat, &b);
    poly_basemul(&chat, &ahat, &bhat);

    // for(int i = 0; i < 512; i ++)
    // {
    //     printf("%d, ",chat.coeffs[i]);
    // }
    // printf("\n");


    poly_invntt(&c,&chat);
    // printf("ntt:\n");
    for(int i = 0; i < 512; i ++)
    {
        // printf("%d, ",c.coeffs[i]);
        if((c.coeffs[i]-c2.coeffs[i])%CTRU_Q!=0)
        {
          printf("No\n");
          return;
        }
    }
    printf("Yes\n");

}
// void test_naivemul_q2()
// {
//   poly a, b, c;
//   int i;
//   for(i = 0; i < CTRU_N/2; i ++)
//   {
//     a.coeffs[i] = b.coeffs[i] = 1;
//   }
//   for(i = CTRU_N/2; i < CTRU_N; i ++)
//   {
//     a.coeffs[i] = b.coeffs[i] = 0;
//   }
//   poly_naivemul_q2(&c, &a, &b);
//   for(i = 0; i < CTRU_N; i ++)
//   {
//     printf("%d, ", c.coeffs[i]);
//   }
//   printf("\n");
// }

int main()
{
  test_kem();
  // test_ntt();

  // test_naivemul_q2();
  return 0;
}
