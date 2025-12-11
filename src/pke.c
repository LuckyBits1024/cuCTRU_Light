#include "params.h"
#include "poly.h"
#include "pack.h"
#include <stdio.h>


void poly_naivemul(poly *c, const poly *a, const poly *b) {
    unsigned int i,j;
    int16_t r[2*CTRU_N] = {0};

    for(i = 0; i < CTRU_N; i++)
        for(j = 0; j < CTRU_N; j++)
            r[i+j] = (r[i+j] + a->coeffs[i]*b->coeffs[j]) % CTRU_Q;

    for(i = CTRU_N; i < 2*CTRU_N; i++)
        r[i-CTRU_N] = (r[i-CTRU_N] - r[i]) % CTRU_Q;

    for(i = 0; i < CTRU_N; i++)
        c->coeffs[i] = r[i] % CTRU_Q;

}

int crypto_pke_keygen(unsigned char pk[CTRU_PKE_PUBLICKEYBYTES],
                      unsigned char sk[CTRU_PKE_SECRETKEYBYTES],
                      const unsigned char coins[CTRU_KEYGEN_COIN_BYTES])
{
  int r;
  poly f, fhat, g, hhat, finv;

  poly_sample_keygen(&f, coins);
  poly_sample_keygen(&g, coins + CTRU_KEYGEN_COIN_BYTES / 2);

  poly_multi_p(&f, &f);
  f.coeffs[0] += 1;

  poly_ntt_small(&fhat, &f);
  poly_ntt_small(&g, &g);

  r = poly_baseinv(&finv, &fhat);
  poly_basemul(&hhat, &g, &finv);
  poly_doublemont(&hhat);

  poly_freeze(&hhat);
  poly_freeze(&fhat);

  //  printf("\nhhat before pack\n");
  // for(int i = 0; i < CTRU_N; i ++)
  // {
  //   printf("%d, ",hhat.coeffs[i]);
  // }
  // printf("\n");
  pack_pk(pk, &hhat);
  pack_sk(sk, &fhat);

  return r;
}

#if (KEM_TYPE == RLWE_KEM)
void crypto_pke_enc(unsigned char ct[CTRU_PKE_CIPHERTEXTBYTES],
                    const unsigned char pk[CTRU_PKE_PUBLICKEYBYTES],
                    const unsigned char m[CTRU_MSGBYTES],
                    const unsigned char coins[CTRU_ENC_COIN_BYTES])
{
  poly r, sigmahat, chat, hhat, e, mp;

  unpack_pk(&hhat, pk);

  poly_sample_enc(&r, coins);
  poly_sample_enc(&e, coins + CTRU_ENC_COIN_BYTES / 2);

  poly_ntt_small(&r, &r);
  poly_basemul(&sigmahat, &hhat, &r);

  poly_encode_e8(&mp, m);
  poly_add(&mp, &e, &mp);
  poly_ntt_big(&mp, &mp);
  poly_add(&chat, &sigmahat, &mp);

  poly_freeze(&chat);
  pack_ct(ct, &chat);
}

void crypto_pke_dec(unsigned char m[CTRU_MSGBYTES],
                    const unsigned char ct[CTRU_PKE_CIPHERTEXTBYTES],
                    const unsigned char sk[CTRU_PKE_SECRETKEYBYTES])
{
  poly chat, fhat, mp;

  unpack_ct(&chat, ct);
  unpack_sk(&fhat, sk);

  poly_basemul(&mp, &chat, &fhat);
  poly_invntt(&mp, &mp);
  poly_fqcsubq(&mp);

  poly_decode_e8(m, &mp);
}

#elif (KEM_TYPE == RLWR_KEM)
void crypto_pke_enc(unsigned char ct[CTRU_PKE_CIPHERTEXTBYTES],
                    const unsigned char pk[CTRU_PKE_PUBLICKEYBYTES],
                    const unsigned char m[CTRU_MSGBYTES],
                    const unsigned char coins[CTRU_ENC_COIN_BYTES])
{
  poly r, sigmahat, c, hhat, mp;

  unpack_pk(&hhat, pk);

  
  poly_sample_enc(&r, coins);
  poly_ntt_small(&r, &r);
  poly_basemul(&sigmahat, &hhat, &r);
  poly_invntt_without_mont(&sigmahat, &sigmahat);
  poly_freeze(&sigmahat);
  poly_encode_e8(&mp, m);
  poly_add(&c, &sigmahat, &mp);
  // poly_freeze(&c);
  poly_fqcsubq(&c);
  poly_compress(&c, &c);
  

  pack_ct(ct, &c);
}

void crypto_pke_dec(unsigned char m[CTRU_MSGBYTES],
                    const unsigned char ct[CTRU_PKE_CIPHERTEXTBYTES],
                    const unsigned char sk[CTRU_PKE_SECRETKEYBYTES])
{
    poly mp,c,fhat;

    unpack_ct(&c, ct);
    unpack_sk(&fhat, sk);
    poly_decompress(&c, &c);

    poly_invntt_without_mont(&fhat,&fhat);

    poly_naivemul(&mp, &c, &fhat);
    poly_fqcsubq(&mp);
    poly_decode_e8(m, &mp);
}

//void crypto_pke_dec(unsigned char m[CTRU_MSGBYTES],
//                    const unsigned char ct[CTRU_PKE_CIPHERTEXTBYTES],
//                    const unsigned char sk[CTRU_PKE_SECRETKEYBYTES])
//{
//  poly c, fhat, mp;
//
//  unpack_ct(&c, ct);
//
//  unpack_sk(&fhat, sk);
//
//
//  poly_decompress(&c, &c);
//  poly_ntt_big(&c, &c);
//  poly_basemul(&mp, &c, &fhat);
//  poly_invntt(&mp, &mp);
//  poly_fqcsubq(&mp);
//  // poly_freeze(&mp);
//  poly_decode_e8(m, &mp);
//}
#endif
