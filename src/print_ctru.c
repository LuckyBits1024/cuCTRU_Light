#include "params.h"
#include "poly.h"
#include "pack.h"
#include <stdio.h>


static void print_poly(poly input)
{
  int i, k;
  int cnt = 0;
  for(i = 0; i < CTRU_N/16; i ++)
  {
    for(k = 0; k < 8; k ++)
    {
      printf("%04hX",input.coeffs[cnt++]);
      printf(" ");
      printf("%04hX",input.coeffs[cnt++]);
      if(k != 7)
        printf(" ");
    }
    printf("\n");
  }
  printf("\n");
}

int pke_keygen(unsigned char pk[CTRU_PKE_PUBLICKEYBYTES],
                      unsigned char sk[CTRU_PKE_SECRETKEYBYTES],
                      const unsigned char coins[CTRU_COINBYTES_KEYGEN])
{
  int r;
  poly f, g, hhat, finv;

  poly_sample_keygen(&f, coins);
  poly_sample_keygen(&g, coins + CTRU_COINBYTES_KEYGEN / 2);

  printf("多项式f'\n");
  print_poly(f);
  printf("多项式g\n");
  print_poly(g);

  pack_sk_f(sk, &f);

  printf("CTRU.PKE.KeyGen 私钥sk\n");
  print_poly(f);

  poly_double(&f, &f);
  f.coeffs[0] += 1;
  
  poly_ntt(&f);
  poly_ntt(&g);

  r = poly_baseinv_opt(&finv, &f);
  poly_basemul(&hhat, &g, &finv);
  poly_invntt(&hhat);
  poly_freeze(&hhat);
  printf("CTRU.PKE.KeyGen 公钥pk\n");
  print_poly(hhat);
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
  printf("CTRU.PKE.Enc 密文ct\n");
  print_poly(c);
  pack_ct(ct, &c);
}

void pke_dec(unsigned char m[CTRU_MSGBYTES],
                    const unsigned char ct[CTRU_PKE_CIPHERTEXTBYTES],
                    const unsigned char sk[CTRU_PKE_SECRETKEYBYTES])
{
  poly c, fhat;

  unpack_ct(&c, ct);
  unpack_sk_f(&fhat, sk);
  
  poly_decode_opt(m, &c, &fhat);
}