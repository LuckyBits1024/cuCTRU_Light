#include <stddef.h>
#include "randombytes.h"
#include "symmetric_crypto.h"
#include "params.h"
#include "pke.h"
#include "poly.h"
#include "pack.h"

#if (USING_NEV_FO == 0) 

#if (USING_SHA2_VARIANT == 0)

int crypto_kem_keygen(unsigned char *pk,
                      unsigned char *sk)
{
  unsigned int i;
  unsigned char coins[CTRU_KEYGEN_COIN_BYTES];

  do
  {
    randombytes(coins, CTRU_SEEDBYTES);
    crypto_hash_shake256(coins, CTRU_KEYGEN_COIN_BYTES, coins, CTRU_SEEDBYTES);

  } while (crypto_pke_keygen(pk, sk, coins));

  for (i = 0; i < CTRU_PKE_PUBLICKEYBYTES; ++i)
    sk[i + CTRU_PKE_SECRETKEYBYTES] = pk[i];

  randombytes(sk + CTRU_PKE_SECRETKEYBYTES + CTRU_PKE_PUBLICKEYBYTES, CTRU_SEEDBYTES);

  return 0;
}

int crypto_kem_encaps(unsigned char *ct,
                      unsigned char *k,
                      const unsigned char *pk)
{
  unsigned int i;
  unsigned char buf[CTRU_SHAREDKEYBYTES + CTRU_ENC_COIN_BYTES], m[CTRU_PREFIXHASHBYTES + CTRU_MSGBYTES];

  randombytes(buf, CTRU_SEEDBYTES);
  crypto_hash_shake256(m + CTRU_PREFIXHASHBYTES, CTRU_MSGBYTES, buf, CTRU_SEEDBYTES);

  for (i = 0; i < CTRU_PREFIXHASHBYTES; ++i)
    m[i] = pk[i];

  crypto_hash_sha3_512(buf, m, CTRU_PREFIXHASHBYTES + CTRU_MSGBYTES);
  crypto_hash_shake256(buf + CTRU_SHAREDKEYBYTES, CTRU_ENC_COIN_BYTES, buf + CTRU_SHAREDKEYBYTES, CTRU_SEEDBYTES);

  crypto_pke_enc(ct, pk, m + CTRU_PREFIXHASHBYTES, buf + CTRU_SHAREDKEYBYTES);

  for (i = 0; i < CTRU_SHAREDKEYBYTES; ++i)
    k[i] = buf[i];

  return 0;
}

int crypto_kem_decaps(unsigned char *k,
                      const unsigned char *ct,
                      const unsigned char *sk)
{
  unsigned int i;
  unsigned char buf[CTRU_SHAREDKEYBYTES + CTRU_ENC_COIN_BYTES], buf2[CTRU_SHAREDKEYBYTES * 2], m[CTRU_PREFIXHASHBYTES + CTRU_MSGBYTES];
  unsigned char ct2[CTRU_PKE_CIPHERTEXTBYTES + CTRU_SEEDBYTES + CTRU_PREFIXHASHBYTES];
  int16_t t;
  int32_t fail;

  crypto_pke_dec(m + CTRU_PREFIXHASHBYTES, ct, sk);

  for (i = 0; i < CTRU_PREFIXHASHBYTES; ++i)
    m[i] = sk[i + CTRU_PKE_SECRETKEYBYTES];

  crypto_hash_sha3_512(buf, m, CTRU_PREFIXHASHBYTES + CTRU_MSGBYTES);
  crypto_hash_shake256(buf + CTRU_SHAREDKEYBYTES, CTRU_ENC_COIN_BYTES, buf + CTRU_SHAREDKEYBYTES, CTRU_SEEDBYTES);

  crypto_pke_enc(ct2, sk + CTRU_PKE_SECRETKEYBYTES, m + CTRU_PREFIXHASHBYTES, buf + CTRU_SHAREDKEYBYTES);

  t = 0;
  for (i = 0; i < CTRU_PKE_CIPHERTEXTBYTES; ++i)
    t |= ct[i] ^ ct2[i];

  fail = (uint16_t)t;
  fail = (-fail) >> 31;

  // for (i = 0; i < CTRU_PREFIXHASHBYTES; ++i)
  //   ct2[i] = sk[i + CTRU_PKE_SECRETKEYBYTES];
  // for (i = 0; i < CTRU_SEEDBYTES; ++i)
  //   ct2[i + CTRU_PREFIXHASHBYTES] = sk[i + CTRU_PKE_SECRETKEYBYTES + CTRU_PKE_PUBLICKEYBYTES];
  // for (i = 0; i < CTRU_PKE_CIPHERTEXTBYTES; ++i)
  //   ct2[i + CTRU_PREFIXHASHBYTES + CTRU_SEEDBYTES] = ct[i];
  // crypto_hash_sha3_512(buf2, ct2, CTRU_PKE_CIPHERTEXTBYTES + CTRU_SEEDBYTES + CTRU_PREFIXHASHBYTES);
   //concatenate c
  for(i = 0; i < CTRU_KEM_CIPHERTEXTBYTES; ++i)
  {
    ct2[i] = ct[i];
  }
  //concatenate z
  for(i = 0; i < CTRU_SEEDBYTES; ++i)
  {
    ct2[i+CTRU_KEM_CIPHERTEXTBYTES] = sk[i + CTRU_PKE_SECRETKEYBYTES + CTRU_PKE_PUBLICKEYBYTES];
  }
  //H(c,z)
  crypto_hash_sha3_512(buf2, ct2, CTRU_KEM_CIPHERTEXTBYTES + CTRU_SEEDBYTES);
  for (i = 0; i < CTRU_SHAREDKEYBYTES; ++i)
    k[i] = buf[i] ^ ((-fail) & (buf[i] ^ buf2[i]));

  return fail;
}

#elif (USING_SHA2_VARIANT == 1) 
static const unsigned char nonce[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

int crypto_kem_keygen(unsigned char *pk,
                      unsigned char *sk)
{
  unsigned int i;
  unsigned char coins[CTRU_KEYGEN_COIN_BYTES];

  do
  {
    randombytes(coins, CTRU_SEEDBYTES);
    crypto_stream(coins, CTRU_KEYGEN_COIN_BYTES, nonce, coins);
  } while (crypto_pke_keygen(pk, sk, coins));

  for (i = 0; i < CTRU_PKE_PUBLICKEYBYTES; ++i)
    sk[i + CTRU_PKE_SECRETKEYBYTES] = pk[i];

  return 0;
}

int crypto_kem_encaps(unsigned char *ct,
                      unsigned char *k,
                      const unsigned char *pk)
{
  unsigned int i;
  unsigned char buf[CTRU_SHAREDKEYBYTES + CTRU_ENC_COIN_BYTES], m[CTRU_MSGBYTES];

  randombytes(buf, CTRU_SEEDBYTES);
  crypto_stream(m, CTRU_MSGBYTES, nonce, buf);

  crypto_hash_sha_512(buf, m, CTRU_MSGBYTES);
  crypto_stream(buf + CTRU_SHAREDKEYBYTES, CTRU_ENC_COIN_BYTES, nonce, buf + CTRU_SHAREDKEYBYTES);

  crypto_pke_enc(ct, pk, m, buf + CTRU_SHAREDKEYBYTES);

  for (i = 0; i < CTRU_SHAREDKEYBYTES; ++i)
    k[i] = buf[i];

  return 0;
}

int crypto_kem_decaps(unsigned char *k,
                      const unsigned char *ct,
                      const unsigned char *sk)
{
  unsigned int i;
  unsigned char buf[CTRU_SHAREDKEYBYTES + CTRU_ENC_COIN_BYTES], m[CTRU_MSGBYTES];
  unsigned char ct2[CTRU_PKE_CIPHERTEXTBYTES];
  int16_t t;
  int32_t fail;

  crypto_pke_dec(m, ct, sk);

  crypto_hash_sha_512(buf, m, CTRU_MSGBYTES);
  crypto_stream(buf + CTRU_SHAREDKEYBYTES, CTRU_ENC_COIN_BYTES, nonce, buf + CTRU_SHAREDKEYBYTES);

  crypto_pke_enc(ct2, sk + CTRU_PKE_SECRETKEYBYTES, m, buf + CTRU_SHAREDKEYBYTES);

  t = 0;
  for (i = 0; i < CTRU_PKE_CIPHERTEXTBYTES; ++i)
    t |= ct[i] ^ ct2[i];

  fail = (uint16_t)t;
  fail = (-fail) >> 31;
  for (i = 0; i < CTRU_SHAREDKEYBYTES; ++i)
    k[i] = buf[i] & ~(-fail);

  return fail;
}

#endif

#elif (USING_NEV_FO == 1) 

int crypto_kem_keygen(unsigned char *pk,
                      unsigned char *sk)
{
  unsigned int i;
  unsigned char coins[CTRU_KEYGEN_COIN_BYTES];

  do
  {
    randombytes(coins, CTRU_SEEDBYTES);
    crypto_hash_shake256(coins, CTRU_KEYGEN_COIN_BYTES, coins, CTRU_SEEDBYTES);

  } while (crypto_pke_keygen(pk, sk, coins));

  for (i = 0; i < CTRU_PKE_PUBLICKEYBYTES; ++i)
    sk[i + CTRU_PKE_SECRETKEYBYTES] = pk[i];

  crypto_hash_sha3_256(sk + CTRU_PKE_SECRETKEYBYTES + CTRU_PKE_PUBLICKEYBYTES, pk, CTRU_PKE_PUBLICKEYBYTES);

  randombytes(sk + CTRU_PKE_SECRETKEYBYTES + CTRU_PKE_PUBLICKEYBYTES + CTRU_HASH_PUBLICKEYBYTES, CTRU_SEEDBYTES);

  return 0;
}

int crypto_kem_encaps(unsigned char *ct,
                      unsigned char *k,
                      const unsigned char *pk)
{
  unsigned int i;
  unsigned char buf[CTRU_SHAREDKEYBYTES + CTRU_ENC_COIN_BYTES], m[CTRU_HASH_PUBLICKEYBYTES + CTRU_MSGBYTES];
  unsigned char buf2[CTRU_SHAREDKEYBYTES + CTRU_PKE_CIPHERTEXTBYTES];

  randombytes(buf, CTRU_SEEDBYTES);
  crypto_hash_sha3_256(m, buf, CTRU_SEEDBYTES);
  crypto_hash_sha3_256(m + CTRU_MSGBYTES, pk, CTRU_PKE_PUBLICKEYBYTES);

  crypto_hash_sha3_512(buf, m, CTRU_HASH_PUBLICKEYBYTES + CTRU_MSGBYTES);
  crypto_hash_shake256(buf + CTRU_SHAREDKEYBYTES, CTRU_ENC_COIN_BYTES, buf + CTRU_SHAREDKEYBYTES, CTRU_SEEDBYTES);

  crypto_pke_enc(ct, pk, m, buf + CTRU_SHAREDKEYBYTES);

  for (i = 0; i < CTRU_SHAREDKEYBYTES; ++i)
    buf2[i] = buf[i];

  for (i = 0; i < CTRU_PKE_CIPHERTEXTBYTES; ++i)
    buf2[i + CTRU_SHAREDKEYBYTES] = ct[i];

  crypto_hash_sha3_256(k, buf2, CTRU_SHAREDKEYBYTES + CTRU_PKE_CIPHERTEXTBYTES);

  return 0;
}

int crypto_kem_decaps(unsigned char *k,
                      const unsigned char *ct,
                      const unsigned char *sk)
{
  unsigned int i;
  unsigned char buf[CTRU_SHAREDKEYBYTES + CTRU_ENC_COIN_BYTES], m[CTRU_HASH_PUBLICKEYBYTES + CTRU_MSGBYTES];
  unsigned char ct2[CTRU_PKE_CIPHERTEXTBYTES + CTRU_SHAREDKEYBYTES];
  int16_t t;
  int32_t fail;

  crypto_pke_dec(m, ct, sk);

  for (i = 0; i < CTRU_HASH_PUBLICKEYBYTES; ++i)
    m[i + CTRU_MSGBYTES] = sk[i + CTRU_PKE_SECRETKEYBYTES + CTRU_PKE_PUBLICKEYBYTES];

  crypto_hash_sha3_512(buf, m, CTRU_HASH_PUBLICKEYBYTES + CTRU_MSGBYTES);
  crypto_hash_shake256(buf + CTRU_SHAREDKEYBYTES, CTRU_ENC_COIN_BYTES, buf + CTRU_SHAREDKEYBYTES, CTRU_SEEDBYTES);

  crypto_pke_enc(ct2, sk + CTRU_PKE_SECRETKEYBYTES, m, buf + CTRU_SHAREDKEYBYTES);

  t = 0;
  for (i = 0; i < CTRU_PKE_CIPHERTEXTBYTES; ++i)
    t |= ct[i] ^ ct2[i];

  fail = (uint16_t)t;
  fail = (-fail) >> 31;

  for (i = 0; i < CTRU_SHAREDKEYBYTES; ++i)
    ct2[i] = buf[i] ^ ((-fail) & (buf[i] ^ sk[i + CTRU_PKE_SECRETKEYBYTES + CTRU_PKE_PUBLICKEYBYTES + CTRU_HASH_PUBLICKEYBYTES]));

  for (i = 0; i < CTRU_PKE_CIPHERTEXTBYTES; ++i)
    ct2[i + CTRU_SHAREDKEYBYTES] = ct[i];

  crypto_hash_sha3_256(k, ct2, CTRU_SHAREDKEYBYTES + CTRU_PKE_CIPHERTEXTBYTES);

  return fail;
}

#endif