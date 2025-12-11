#include <stdint.h>
#include "cpucycles.h"
#include "params.h"
#include "reduce.h"
#include "poly.h"
#include "coding.h"
#include "cbd.h"
#include "inverse.h"
#include "ntt.h"
void poly_naivemul_q(poly *c, const poly *a, const poly *b, int16_t q) {

     int16_t r[2*CTRU_N] = {0};

     for(int i = 0; i < CTRU_N; i++)
     {
         for(int j = 0; j < CTRU_N; j++){
             r[i+j] = (r[i+j] + (int32_t)a->coeffs[i]*b->coeffs[j]) % q;
         }

     }
     for(int i = CTRU_N; i < 2*CTRU_N; i++)
      r[i-CTRU_N] = (r[i-CTRU_N] - r[i]);

     for(int i = 0; i < CTRU_N; i++)
      c->coeffs[i] = r[i] % q;
}

void poly_naivemul_q2(poly *c, const poly *a, const poly *b, int16_t Q) {

     int16_t r[2*CTRU_N] = {0};

     for(int i = 0; i < CTRU_N; i++)
     {
         for(int j = 0; j < CTRU_N; j++){
             r[i+j] = (r[i+j] + (int16_t)a->coeffs[i]*b->coeffs[j]);
         }

     }
     for(int i = CTRU_N; i < 2*CTRU_N; i++)
      r[i-CTRU_N] = (r[i-CTRU_N] - r[i]);

     for(int i = 0; i < CTRU_N; i++)
      c->coeffs[i] = r[i] & (Q-1);
}

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
  for (int i = 0; i < CTRU_N; ++i)
    b->coeffs[i] = a->coeffs[i] + a->coeffs[i];

}


void poly_sample_keygen(poly *a, const unsigned char *buf)
{
  cbd1(a, buf);
}

void poly_sample_enc(poly *a, const unsigned char *buf)
{
  cbd1(a, buf);
}



int poly_baseinv(poly *b, const poly *a)
{
  int r = 0;
  for (int i = 0; i < CTRU_N / 16; ++i)
  {
    r += rq_inverse(b->coeffs + 16 * i,
                    a->coeffs + 16 * i,
                    zetas[32 + i]);
    r += rq_inverse(b->coeffs + 16 * i + 8,
                    a->coeffs + 16 * i + 8,
                    -zetas[32 + i]);
  }

  return r;
}
int poly_baseinv_opt(poly *b, const poly *a)
{
  int r = 0;
  for (int i = 0; i < CTRU_N / 16; ++i)
  {
    if(r != 0)   return r;
    r += rq_inverse_opt(b->coeffs + 16 * i,
                    a->coeffs + 16 * i,
                    zetas_withoutmont[32 + i]);
    r += rq_inverse_opt(b->coeffs + 16 * i + 8,
                    a->coeffs + 16 * i + 8,
                    -zetas_withoutmont[32 + i]);
  }

  return r;

}






void poly_encode_compress(poly *c,
                          const poly *sigma,
                          const unsigned char *msg)
{
  unsigned int i, j;
  int16_t mask;
  uint8_t mh[CTRU_N / 8];
  uint8_t tmp;
  // int16_t t;
  int32_t t;
  for (i = 0; i < CTRU_MSGBYTES; i++)
  {
    tmp = msg[i] & 0xF;
    mh[2 * i] = encode_e8(tmp);

    tmp = (msg[i] >> 4) & 0xF;
    mh[2 * i + 1] = encode_e8(tmp);
  }

  for (i = 0; i < CTRU_N / 8; i++)
  {
    for (j = 0; j < 8; j++)
    {
      mask = -(int16_t)((mh[i] >> j) & 1);
      t = (int32_t)((sigma->coeffs[8 * i + j] << CTRU_LOGQ2) + (CTRU_Q >> 1)) / CTRU_Q;
      t = t + (mask & (CTRU_Q2 >> 1));
      c->coeffs[8 * i + j] = t & (CTRU_Q2 - 1);
    }
  }
}

void poly_decode(unsigned char *msg,
                 const poly *c,
                 const poly *f)
{
  unsigned int i, j;
  poly mp;
  uint32_t tmp_mp[8];

  poly_naivemul_q2(&mp, c, f, CTRU_Q2);
  

  for (i = 0; i < CTRU_MSGBYTES; i++)
  {
    msg[i] = 0;
  }

  for (i = 0; i < CTRU_N / 8; i++)
  {
    for (j = 0; j < 8; j++)
    {
      tmp_mp[j] = (uint32_t)mp.coeffs[8 * i + j];
    }
    msg[i >> 1] |= decode_e8(tmp_mp) << ((i & 1) << 2);
  }
}
void poly_decode_opt(unsigned char *msg,
                 const poly *c,
                 const poly *f)
{
  unsigned int i, j;
  poly mp;
  uint32_t tmp_mp[8];

  poly_naivemul_q2(&mp, c, f, CTRU_Q2);
  

  for (i = 0; i < CTRU_MSGBYTES; i++)
  {
    msg[i] = 0;
  }

  for (i = 0; i < CTRU_N / 8; i++)
  {
    for (j = 0; j < 8; j++)
    {
      tmp_mp[j] = (uint32_t)mp.coeffs[8 * i + j];
    }
    msg[i >> 1] |= decode_e8(tmp_mp) << ((i & 1) << 2);
  }
}
void poly_ntt(poly *b)
{
  ntt(b->coeffs);
}

void poly_invntt(poly *b)
{
  invntt(b->coeffs);
}

void poly_basemul(poly *c, const poly *a, const poly *b)
{
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
}


