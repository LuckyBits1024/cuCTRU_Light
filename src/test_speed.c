#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "api.h"
#include "params.h"
#include "ctru.h"
#include "poly.h"
#include "cpucycles.h"
#include "speed.h"
#include "randombytes.h"

#define NTESTS 100000

uint64_t t[NTESTS];

void test_speed_kem()
{
  printf("\n");

  printf("CTRU-%d-%d-KEM\n\n", CTRU_N, CTRU_Q);
  unsigned int i;
  unsigned char k1[CTRU_SHAREDKEYBYTES], k2[CTRU_SHAREDKEYBYTES];
  unsigned char pk[CTRU_KEM_PUBLICKEYBYTES], sk[CTRU_KEM_SECRETKEYBYTES];
  unsigned char ct[CTRU_KEM_CIPHERTEXTBYTES];
  unsigned long long ss_byts1, ss_byts2, ct_byts, pk_byts, sk_byts;

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    kem_keygen(pk, &pk_byts, sk, &sk_byts);
  }
  print_results("ctru_kem_keygen: ", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    kem_enc(pk, pk_byts, k1, &ss_byts1, ct, &ct_byts);
  }
  print_results("ctru_kem_encaps: ", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    kem_dec(sk, sk_byts, ct, ct_byts, k2, &ss_byts2);
  }
  print_results("ctru_kem_decaps: ", t, NTESTS);

}

int main()
{
  test_speed_kem();
  return 0;
}
