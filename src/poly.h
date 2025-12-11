#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

typedef struct
{
  int16_t coeffs[CTRU_N];
} poly;

void poly_reduce(poly *a);
void poly_freeze(poly *a);
void poly_fqcsubq(poly *a);

void poly_add(poly *c, const poly *a, const poly *b);
void poly_double(poly *b, const poly *a);
void poly_multi_p(poly *b, const poly *a);

void poly_tomont(poly *a);
void poly_frommont(poly *a);
void poly_doublemont(poly *a);

void poly_ntt_big(poly *b, const poly *a);
void poly_ntt_small(poly *b, const poly *a);
void poly_invntt(poly *b, const poly *a);
void poly_invntt_without_mont(poly *b, const poly *a);
void poly_basemul(poly *c, const poly *a, const poly *b);
int poly_baseinv(poly *b, const poly *a);

void poly_sample_keygen(poly *a, const unsigned char *buf);
void poly_sample_enc(poly *a, const unsigned char *buf);

void poly_encode_compress(poly *c, const poly *sigma, const unsigned char *m);
void poly_decode(unsigned char *m, const poly *c, const poly *fhat);

void poly_encode_e8(poly *mp, const unsigned char *msg);
void poly_decode_e8(unsigned char *msg, const poly *mp);

void poly_compress(poly *b, const poly *a);
void poly_decompress(poly *b, const poly *a);

#endif
