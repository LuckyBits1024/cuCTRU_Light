#include <stddef.h>
#include "randombytes.h"
#include "symmetric_crypto.h"
#include "params.h"
#include "ctru.h"
#include "api.h"

unsigned char seed[CTRU_SEEDBYTES] = {};
unsigned long long rand_get_sd_byts()
{
  return CTRU_SEEDBYTES;
}
int rand_init(unsigned char * s, unsigned long long s_byts)
{
  randombytes(s, s_byts);
  for(int i = 0; i < s_byts; i ++)
  {
    seed[i] = s[i];
  }
  return 0;
}
int rand_byts(unsigned long long r_byts, unsigned char * r)
{
  crypto_hash_shake256(r, r_byts, seed, CTRU_SEEDBYTES);
  return 0;
}
unsigned long long kem_get_pk_byts()
{
  return CTRU_KEM_PUBLICKEYBYTES;
}
unsigned long long kem_get_sk_byts()
{
  return CTRU_KEM_SECRETKEYBYTES;
}
unsigned long long kem_get_ss_byts()
{
  return CTRU_SHAREDKEYBYTES;
}
unsigned long long kem_get_ct_byts()
{
  return CTRU_KEM_CIPHERTEXTBYTES;
}


// int kem_keygen(unsigned char * pk, unsigned long long * pk_byts, unsigned char * sk, unsigned long long * sk_byts)
// {
//   unsigned int i;
//   unsigned char coins[CTRU_COINBYTES_KEYGEN];

//   do
//   {
//     rand_init(seed, CTRU_SEEDBYTES);
//     rand_byts(CTRU_COINBYTES_KEYGEN, coins);
//   } while (pke_keygen(pk, sk, coins));

//   for (i = 0; i < CTRU_PKE_PUBLICKEYBYTES; ++i)
//     sk[i + CTRU_PKE_SECRETKEYBYTES] = pk[i];

//   randombytes(sk + CTRU_PKE_SECRETKEYBYTES + CTRU_PKE_PUBLICKEYBYTES, CTRU_Z_BYTES);
//   *pk_byts = CTRU_KEM_PUBLICKEYBYTES;
//   *sk_byts = CTRU_KEM_SECRETKEYBYTES;
//   return 0;
// }

int kem_keygen(unsigned char * pk, unsigned long long * pk_byts, unsigned char * sk, unsigned long long * sk_byts)
{
  unsigned int i;
  unsigned char coins[CTRU_COINBYTES_KEYGEN];

  do
  {
    rand_init(seed, CTRU_SEEDBYTES);
    rand_byts(CTRU_COINBYTES_KEYGEN, coins);
  } while (pke_keygen(pk, sk, coins));

  for (i = 0; i < CTRU_PKE_PUBLICKEYBYTES; ++i)
    sk[i + CTRU_PKE_SECRETKEYBYTES] = pk[i];

  randombytes(sk + CTRU_PKE_SECRETKEYBYTES + CTRU_PKE_PUBLICKEYBYTES, CTRU_Z_BYTES);
  *pk_byts = CTRU_KEM_PUBLICKEYBYTES;
  *sk_byts = CTRU_KEM_SECRETKEYBYTES;
  return 0;
}

int kem_enc(unsigned char * pk, unsigned long long pk_byts,
unsigned char * ss, unsigned long long * ss_byts,
unsigned char * ct, unsigned long long * ct_byts)
{
  unsigned int i;
  unsigned char buf[CTRU_SHAREDKEYBYTES + CTRU_COINBYTES_ENC], m[CTRU_MSGBYTES];
  //buf2存ID(pk)||M
  unsigned char buf2[CTRU_PREFIXHASHBYTES + CTRU_MSGBYTES];
  randombytes(m, CTRU_MSGBYTES);

  for (i = 0; i < CTRU_PREFIXHASHBYTES; ++i)
    buf2[i] = pk[i];
  for (i = 0; i < CTRU_MSGBYTES; ++i)
    buf2[i+CTRU_PREFIXHASHBYTES] = m[i];

  crypto_hash_sha3512(buf, buf2, CTRU_PREFIXHASHBYTES + CTRU_MSGBYTES);
  crypto_hash_shake256(buf + CTRU_SHAREDKEYBYTES, CTRU_COINBYTES_ENC, buf + CTRU_SHAREDKEYBYTES, CTRU_SHAREDKEYBYTES);

  pke_enc(ct, pk, m, buf + CTRU_SHAREDKEYBYTES);

  for (i = 0; i < CTRU_SHAREDKEYBYTES; ++i)
    ss[i] = buf[i];
  *ss_byts = CTRU_SHAREDKEYBYTES;
  *ct_byts = CTRU_KEM_CIPHERTEXTBYTES;
  return 0;
}

// int kem_dec(
// unsigned char * sk, unsigned long long sk_byts,
// unsigned char * ct, unsigned long long ct_byts,
// unsigned char * ss, unsigned long long * ss_byts)
// {
//   unsigned int i;
//   //buf存放导出密钥及采样需要的coin，buf2存~K，由于是SHA3_512的输出，需要等于64字节
//   unsigned char buf[CTRU_SHAREDKEYBYTES + CTRU_COINBYTES_ENC], buf2[SHA3512_OUTPUTBYTES], m[CTRU_PREFIXHASHBYTES + CTRU_MSGBYTES];
//   unsigned char ct2[CTRU_PKE_CIPHERTEXTBYTES + CTRU_Z_BYTES + CTRU_PREFIXHASHBYTES];
//   int16_t t;
//   int32_t fail;

//   pke_dec(m + CTRU_PREFIXHASHBYTES, ct, sk);

//   for (i = 0; i < CTRU_PREFIXHASHBYTES; ++i)
//     m[i] = sk[i + CTRU_PKE_SECRETKEYBYTES];

//   crypto_hash_sha3512(buf, m, CTRU_PREFIXHASHBYTES + CTRU_MSGBYTES);
//   crypto_hash_shake256(buf + CTRU_SHAREDKEYBYTES, CTRU_COINBYTES_ENC, buf + CTRU_SHAREDKEYBYTES, CTRU_SHAREDKEYBYTES);

//   pke_enc(ct2, sk + CTRU_PKE_SECRETKEYBYTES, m + CTRU_PREFIXHASHBYTES, buf + CTRU_SHAREDKEYBYTES);

//   t = 0;
//   for (i = 0; i < CTRU_PKE_CIPHERTEXTBYTES; ++i)
//     t |= ct[i] ^ ct2[i];

//   fail = (uint16_t)t;
//   fail = (-fail) >> 31;


//   //concatenate c
//   for(i = 0; i < CTRU_KEM_CIPHERTEXTBYTES; ++i)
//   {
//     ct2[i] = ct[i];
//   }
//   //concatenate z
//   for(i = 0; i < CTRU_Z_BYTES; ++i)
//   {
//     ct2[i+CTRU_KEM_CIPHERTEXTBYTES] = sk[i + CTRU_PKE_SECRETKEYBYTES + CTRU_PKE_PUBLICKEYBYTES];
//   }
//   //H(c,z)
//   crypto_hash_sha3512(buf2, ct2, CTRU_KEM_CIPHERTEXTBYTES + CTRU_Z_BYTES);
 

//   for (i = 0; i < CTRU_SHAREDKEYBYTES; ++i)
//     ss[i] = buf[i] ^ ((-fail) & (buf[i] ^ buf2[i]));
//   *ss_byts = CTRU_SHAREDKEYBYTES;
//   return fail;
// }

int kem_dec(
unsigned char * sk, unsigned long long sk_byts,
unsigned char * ct, unsigned long long ct_byts,
unsigned char * ss, unsigned long long * ss_byts)
{
  unsigned int i;
  //buf存放导出密钥及采样需要的coin，buf2存~K，由于是SHA3_512的输出，需要等于64字节
  unsigned char buf[CTRU_SHAREDKEYBYTES + CTRU_COINBYTES_ENC], buf2[SHA3512_OUTPUTBYTES], m[CTRU_PREFIXHASHBYTES + CTRU_MSGBYTES];
  unsigned char ct2[CTRU_PKE_CIPHERTEXTBYTES + CTRU_Z_BYTES + CTRU_PREFIXHASHBYTES];
  int16_t t;
  int32_t fail;

  pke_dec(m + CTRU_PREFIXHASHBYTES, ct, sk);

  for (i = 0; i < CTRU_PREFIXHASHBYTES; ++i)
    m[i] = sk[i + CTRU_PKE_SECRETKEYBYTES];

  crypto_hash_sha3512(buf, m, CTRU_PREFIXHASHBYTES + CTRU_MSGBYTES);
  crypto_hash_shake256(buf + CTRU_SHAREDKEYBYTES, CTRU_COINBYTES_ENC, buf + CTRU_SHAREDKEYBYTES, CTRU_SHAREDKEYBYTES);

  pke_enc(ct2, sk + CTRU_PKE_SECRETKEYBYTES, m + CTRU_PREFIXHASHBYTES, buf + CTRU_SHAREDKEYBYTES);

  t = 0;
  for (i = 0; i < CTRU_PKE_CIPHERTEXTBYTES; ++i)
    t |= ct[i] ^ ct2[i];

  fail = (uint16_t)t;
  fail = (-fail) >> 31;


  //concatenate c
  for(i = 0; i < CTRU_KEM_CIPHERTEXTBYTES; ++i)
  {
    ct2[i] = ct[i];
  }
  //concatenate z
  for(i = 0; i < CTRU_Z_BYTES; ++i)
  {
    ct2[i+CTRU_KEM_CIPHERTEXTBYTES] = sk[i + CTRU_PKE_SECRETKEYBYTES + CTRU_PKE_PUBLICKEYBYTES];
  }
  //H(c,z)
  crypto_hash_sha3512(buf2, ct2, CTRU_KEM_CIPHERTEXTBYTES + CTRU_Z_BYTES);
 

  for (i = 0; i < CTRU_SHAREDKEYBYTES; ++i)
    ss[i] = buf[i] ^ ((-fail) & (buf[i] ^ buf2[i]));
  *ss_byts = CTRU_SHAREDKEYBYTES;
  return fail;
}

