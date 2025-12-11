#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "randombytes.h"
#include "ctru.h"
#include "params.h"
#include "api.h"
#include "cpucycles.h"
#include "speed.h"

#define NTESTS 10000
uint64_t t[NTESTS];

void test_pke()
{
  unsigned int i;
  unsigned char ct[CTRU_PKE_CIPHERTEXTBYTES];
  unsigned char m[CTRU_MSGBYTES],  m2[CTRU_MSGBYTES];
  unsigned char coins[CTRU_COINBYTES_KEYGEN];
  unsigned char pk[CTRU_PKE_PUBLICKEYBYTES], sk[CTRU_PKE_SECRETKEYBYTES];
  randombytes(coins, CTRU_COINBYTES_KEYGEN);
  randombytes(m, CTRU_MSGBYTES);


    for (i = 0; i < NTESTS; i++)
    {
        t[i] = cpucycles();
        pke_keygen(pk, sk, coins);
    }
    print_results("ctru_pke_keygen: ", t, NTESTS);

    for (i = 0; i < NTESTS; i++)
    {
        t[i] = cpucycles();
        pke_enc(ct, pk, m, coins);
    }
    print_results("ctru_pke_enc: ", t, NTESTS);

    for (i = 0; i < NTESTS; i++)
    {
        t[i] = cpucycles();
        pke_dec(m2, ct, sk);
    }
    print_results("ctru_pke_dec: ", t, NTESTS);


 for (i = 0; i < NTESTS; i++)
 {
   pke_keygen(pk, sk, coins);
   pke_enc(ct, pk, m, coins);
   pke_dec(m2, ct, sk);
   for(int j = 0; j < CTRU_MSGBYTES; j ++)
   {
       if(m[j]!=m2[j])
       {
           printf("error key\n");
           return;
       }
   }
 }
}

int main()
{
  test_pke();
  return 0;
}
