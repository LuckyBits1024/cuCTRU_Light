#include "params.h"
#include "poly.h"
#include "pack.h"
#include <stdio.h>
// int pke_keygen(unsigned char pk[CTRU_PKE_PUBLICKEYBYTES],
//                       unsigned char sk[CTRU_PKE_SECRETKEYBYTES],
//                       const unsigned char coins[CTRU_COINBYTES_KEYGEN])
// {
//   int r;
//   poly f, g, hhat, finv;

//   poly_sample_keygen(&f, coins);
//   poly_sample_keygen(&g, coins + CTRU_COINBYTES_KEYGEN / 2);

//   poly_double(&f, &f);
//   f.coeffs[0] += 1;
//   //先pack，因为后面值会变
//   pack_sk(sk, &f);
//   poly_ntt(&f);
//   poly_ntt(&g);

//   r = poly_baseinv(&finv, &f);
//   poly_basemul(&hhat, &g, &finv);
//   poly_invntt(&hhat);
//   poly_freeze(&hhat);
//   pack_pk(pk, &hhat);

//   return r;
// }

int pke_keygen(unsigned char pk[CTRU_PKE_PUBLICKEYBYTES],
                      unsigned char sk[CTRU_PKE_SECRETKEYBYTES],
                      const unsigned char coins[CTRU_COINBYTES_KEYGEN])
{
  int r;
  poly f, g, hhat, finv;

  poly_sample_keygen(&f, coins);
  poly_sample_keygen(&g, coins + CTRU_COINBYTES_KEYGEN / 2);

  poly_double(&f, &f);
  f.coeffs[0] += 1;
  //先pack，因为后面值会变
  pack_sk(sk, &f);
  poly_ntt(&f);
  poly_ntt(&g);

  r = poly_baseinv_opt(&finv, &f);
  poly_basemul(&hhat, &g, &finv);
  poly_invntt(&hhat);
  poly_freeze(&hhat);
  pack_pk(pk, &hhat);

  return r;
}

void pke_enc(unsigned char ct[CTRU_PKE_CIPHERTEXTBYTES],
                    const unsigned char pk[CTRU_PKE_PUBLICKEYBYTES],
                    const unsigned char m[CTRU_MSGBYTES],
                    const unsigned char coins[CTRU_COINBYTES_ENC])
{
  poly r, sigma, c, hhat;

  unpack_pk(&hhat, pk);
  poly_sample_enc(&r, coins);
  poly_ntt(&hhat);
  poly_ntt(&r);
  poly_basemul(&sigma, &hhat, &r);
  poly_invntt(&sigma);
  poly_freeze(&sigma);
  poly_encode_compress(&c, &sigma, m);
  pack_ct(ct, &c);
}

void pke_dec(unsigned char m[CTRU_MSGBYTES],
                    const unsigned char ct[CTRU_PKE_CIPHERTEXTBYTES],
                    const unsigned char sk[CTRU_PKE_SECRETKEYBYTES])
{
  poly c, fhat;

  unpack_ct(&c, ct);
  unpack_sk(&fhat, sk);
  
  poly_decode(m, &c, &fhat);
}