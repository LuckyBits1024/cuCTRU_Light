#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "params.h"
#include "api.h"
#include "poly.h"
#include "randombytes.h"
#define NTESTS 1

void test_kem()
{
  unsigned int i, j;
  unsigned char k1[CTRU_SHAREDKEYBYTES], k2[CTRU_SHAREDKEYBYTES];
  unsigned char pk[CTRU_KEM_PUBLICKEYBYTES], sk[CTRU_KEM_SECRETKEYBYTES];
  unsigned char ct[CTRU_KEM_CIPHERTEXTBYTES];
  unsigned long long ss_byts1, ss_byts2, ct_byts, pk_byts, sk_byts;

  for (i = 0; i < NTESTS; i++)
  {
    kem_keygen(pk, &pk_byts, sk, &sk_byts);
    kem_enc(pk, pk_byts, k1, &ss_byts1, ct, &ct_byts);
    kem_dec(sk, sk_byts, ct, ct_byts, k2, &ss_byts2);

    for (j = 0; j < CTRU_SHAREDKEYBYTES; j++)
      if (k1[j] != k2[j])
      {
        printf("Round %d. Failure: Keys dont match: %hhx != %hhx!\n", i, k1[j], k2[j]);
        return;
      }
  }

  printf("CTRU-%d-KEM is correct!\n", CTRU_N);

  printf("Test %d times.\n\n", NTESTS);
  printf("CTRU_N = %d, CTRU_Q = %d, CTRU_Q2 = %d\n", CTRU_N, CTRU_Q, CTRU_Q2);
  printf("KEM size: sk = %d bytes, pk = %d bytes, ct = %d bytes\n\n",
        CTRU_KEM_SECRETKEYBYTES, CTRU_KEM_PUBLICKEYBYTES, CTRU_KEM_CIPHERTEXTBYTES);
}

int main()
{
  test_kem();
  return 0;
}
