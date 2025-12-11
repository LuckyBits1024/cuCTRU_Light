#include <stdint.h>
#include "cpucycles.h"
#include "params.h"
#include "reduce.h"
#include "ntt.h"
#include "poly.h"
#include "coding.h"
#include "cbd.h"
#include "inverse.h"

void poly_reduce(poly *a)
{
  for (int i = 0; i < CTRU_N; ++i)
    a->coeffs[i] = barrett_reduce(a->coeffs[i]);
}

void poly_freeze(poly *a)
{
  poly_reduce(a);
  for (int i = 0; i < CTRU_N; ++i)
    a->coeffs[i] = fqcsubq(a->coeffs[i]);
}

void poly_fqcsubq(poly *a)
{
  for (int i = 0; i < CTRU_N; ++i)
    a->coeffs[i] = fqcsubq(a->coeffs[i]);
}

void poly_add(poly *c, const poly *a, const poly *b)
{
  for (int i = 0; i < CTRU_N; ++i)
    c->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}

void poly_double(poly *b, const poly *a)
{
  for (int i = 0; i < CTRU_N; ++i)
    b->coeffs[i] = 2 * a->coeffs[i];
}

void poly_multi_p(poly *b, const poly *a)
{
#if (CTRU_N == 512)
  for (int i = 0; i < CTRU_N; ++i)
    b->coeffs[i] = a->coeffs[i] + a->coeffs[i];
#elif (CTRU_N == 1024)
  for (int i = 0; i < 512; ++i)
  {
    int16_t tmp = a->coeffs[i];
    b->coeffs[i] = a->coeffs[512 + i] + tmp;
    b->coeffs[512 + i] = a->coeffs[512 + i] - tmp;
  }
#endif
}

void poly_tomont(poly *a)
{
  const int16_t t = (MONT * MONT) % CTRU_Q;
  for (int i = 0; i < CTRU_N; ++i)
    a->coeffs[i] = fqmul(a->coeffs[i], t);
}

void poly_doublemont(poly *a)
{
  const int16_t t = (((MONT * MONT) % CTRU_Q) * MONT) % CTRU_Q;
  for (int i = 0; i < CTRU_N; ++i)
    a->coeffs[i] = fqmul(a->coeffs[i], t);
}

void poly_frommont(poly *a)
{
  for (int i = 0; i < CTRU_N; ++i)
    a->coeffs[i] = fqmul(a->coeffs[i], 1);
}

#if (KEM_TYPE == RLWE_KEM)
void poly_sample_keygen(poly *a, const unsigned char *buf)
{
#if (CTRU_N == 512)
  cbd1(a, buf);
#elif (CTRU_N == 1024)
  ternary_1of6(a, buf);
#endif
}

void poly_sample_enc(poly *a, const unsigned char *buf)
{
#if (CTRU_N == 512)
  cbd1(a, buf);
#elif (CTRU_N == 1024)
  ternary_1of6(a, buf);
#endif
}
#elif (KEM_TYPE == RLWR_KEM)
void poly_sample_keygen(poly *a, const unsigned char *buf)
{
#if (CTRU_N == 512)
  cbd1(a, buf);
#elif (CTRU_N == 1024)
  ternary_1of5(a, buf);
#endif
}

void poly_sample_enc(poly *a, const unsigned char *buf)
{
#if (CTRU_N == 512)
  cbd1(a, buf);
#elif (CTRU_N == 1024)
  cbd1(a, buf);
#endif
}

#endif

void poly_ntt_big(poly *b, const poly *a)
{
  ntt_big(b->coeffs, a->coeffs);
}

void poly_ntt_small(poly *b, const poly *a)
{
  ntt_small(b->coeffs, a->coeffs);
}

void poly_invntt(poly *b, const poly *a)
{
  invntt(b->coeffs, a->coeffs);
}

void poly_invntt_without_mont(poly *b, const poly *a)
{
  invntt_without_mont(b->coeffs, a->coeffs);
}

void poly_basemul(poly *c, const poly *a, const poly *b)
{
#if (CTRU_N == 512)
  for (int i = 0; i < CTRU_N / 16; ++i)
  {
    basemul(c->coeffs + 16 * i,
            a->coeffs + 16 * i,
            b->coeffs + 16 * i,
            zetas[32 + i]);
    basemul(c->coeffs + 16 * i + 8,
            a->coeffs + 16 * i + 8,
            b->coeffs + 16 * i + 8,
            -zetas[32 + i]);
  }
#elif (CTRU_N == 1024)
  for (int i = 0; i < CTRU_N / 32; ++i)
  {
    basemul(c->coeffs + 32 * i,
            a->coeffs + 32 * i,
            b->coeffs + 32 * i,
            zetas[32 + i]);
    basemul(c->coeffs + 32 * i + 16,
            a->coeffs + 32 * i + 16,
            b->coeffs + 32 * i + 16,
            -zetas[32 + i]);
  }
#endif
}

int poly_baseinv(poly *b, const poly *a)
{
  int r = 0;
#if (CTRU_N == 512)
  for (int i = 0; i < CTRU_N / 16; ++i)
  {
    r += rq_inverse(b->coeffs + 16 * i,
                    a->coeffs + 16 * i,
                    zetas[32 + i]);
    r += rq_inverse(b->coeffs + 16 * i + 8,
                    a->coeffs + 16 * i + 8,
                    -zetas[32 + i]);
  }
#elif (CTRU_N == 1024)
  for (int i = 0; i < CTRU_N / 32; ++i)
  {
    r += rq_inverse(b->coeffs + 32 * i,
                    a->coeffs + 32 * i,
                    zetas[32 + i]);
    r += rq_inverse(b->coeffs + 32 * i + 16,
                    a->coeffs + 32 * i + 16,
                    -zetas[32 + i]);
  }
#endif
  return r;
}

void poly_encode_e8(poly *mp, const unsigned char *msg)
{
  unsigned int i, j;
  int16_t mask;
  uint8_t mh[CTRU_MSGBYTES * 2];
  uint8_t tmp;

  for (i = 0; i < CTRU_MSGBYTES; i++)
  {
    tmp = msg[i] & 0xF;
    mh[2 * i] = encode_e8(tmp);

    tmp = (msg[i] >> 4) & 0xF;
    mh[2 * i + 1] = encode_e8(tmp);
  }

  for (i = 0; i < 512 / 8; i++)
  {
    for (j = 0; j < 8; j++)
    {
      mask = -(int16_t)((mh[i] >> j) & 1);
      mp->coeffs[8 * i + j] = (mask & ((CTRU_Q) >> 1));
    }
  }

#if (CTRU_N == 1024)
  for (i = 0; i < 512; i++)
  {
    mp->coeffs[512 + i] = mp->coeffs[i];
  }
#endif
  // poly_freeze(mp);
}

void poly_decode_e8(unsigned char *msg, const poly *mp)
{
  unsigned int i, j;
  uint32_t tmp_mp[8];

  for (i = 0; i < CTRU_MSGBYTES; i++)
  {
    msg[i] = 0;
  }

  for (i = 0; i < 512 / 8; i++)
  {
    for (j = 0; j < 8; j++)
    {
      tmp_mp[j] = (uint32_t)mp->coeffs[8 * i + j];
    }
    msg[i >> 1] |= decode_e8(tmp_mp) << ((i & 1) << 2);
  }
}

#if (KEM_TYPE == RLWR_KEM)
void poly_compress(poly *b, const poly *a)
{
  int16_t t;
  for (int i = 0; i < CTRU_N; i++)
  {
    t = (int32_t)((a->coeffs[i] << CTRU_LOGQ2) + (CTRU_Q >> 1)) / CTRU_Q;
    b->coeffs[i] = t & (CTRU_Q2 - 1);
  }
}

void poly_decompress(poly *b, const poly *a)
{
  for (int i = 0; i < CTRU_N; i++)
  {
    b->coeffs[i] = (a->coeffs[i] * CTRU_Q + (CTRU_LOGQ2 >> 1)) >> CTRU_LOGQ2;
  }
}
#endif
