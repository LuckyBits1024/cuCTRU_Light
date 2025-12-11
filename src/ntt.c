#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "reduce.h"
#include "ntt.h"
#include "inverse.h"

int16_t zetas_withoutmont[64] = {-1, -62, 40, 173, 27, 136, -311, -57, 144, -300, -377, -304, -43, -359, 182, -251, -25, -12, 231, -289, -94, 324, -85, 113, -245, 190, -197, 90, -306, 253, -64, -123, -5, -310, 200, 96, 
135, -89, -17, -285, -49, 38, -347, 18, -215, -257, 141, 283, -125, -60, -383, 93, 299, 82, 344, -204, 313, 181, -216, -319, 8, -273, -320, 154};

int16_t zetas[64] = {
    171, 605, 81, 408, 766, 583, 120, 519, 753, 546, 640, 461, 432, 638, 407, 626,
    430, 514, 487, 203, 694, 733, 693, 671, 369, 577, 620, 759, 34, 570, 178, 270,
    86, 718, 405, 502, 754, 608, 600, 288, 689, 423, 124, 767, 622, 114, 497, 54,
    612, 263, 128, 246, 394, 589, 389, 279, 307, 578, 24, 719, 170, 543, 121, 581};

int16_t zetas_inv[64] = {
    -581, -121, -543, -170, -719, -24, -578, -307, -279, -389, -589, -394, -246, -128, -263, -612,
    -54, -497, -114, -622, -767, -124, -423, -689, -288, -600, -608, -754, -502, -405, -718, -86,
    -270, -178, -570, -34, -759, -620, -577, -369, -671, -693, -733, -694, -203, -487, -514, -430,
    -626, -407, -638, -432, -461, -640, -546, -753, -519, -120, -583, -766, -408, -81, -605, 541};  // 64^(-1) mod 769

void ntt(int16_t *a)
{
  unsigned int len, start, j, k = 1;
  int16_t t, zeta;



  for (len = CTRU_N / 2; len >= ROOT_DIMENSION; len >>= 1)
  {
    for (start = 0; start < CTRU_N; start = j + len)
    {
      zeta = zetas[k++];
      for (j = start; j < start + len; ++j)
      {
        t = fqmul(zeta, a[j + len]);
        a[j + len] = barrett_reduce(a[j] - t);
        a[j] = barrett_reduce(a[j] + t);
      }
    }
  }
}



void invntt(int16_t *a)
{
  unsigned int start, len, j, k;
  int16_t t, zeta;


  k = 0;
  for (len = ROOT_DIMENSION; len <= CTRU_N / 2; len <<= 1)
  {
    for (start = 0; start < CTRU_N; start = j + len)
    {
      zeta = zetas_inv[k++];
      for (j = start; j < start + len; ++j)
      {
        t = a[j];
        a[j] = barrett_reduce(t + a[j + len]);
        a[j + len] = barrett_reduce(t - a[j + len]);
        a[j + len] = fqmul(zeta, a[j + len]);
      }
    }
  }

  for (j = 0; j < CTRU_N; ++j)
  {
    a[j] = fqmul(a[j], zetas_inv[63]);
  }
}


#define CALC_D(a, b, x, y, d) (fqmul((a[x] + a[y]), (b[x] + b[y])) - d[x] - d[y])

void basemul(int16_t *c, const int16_t *a, const int16_t *b, const int16_t zeta)
{
#if (CTRU_N == 512)
  int16_t d[8];
  for (int i = 0; i < 8; i++)
    d[i] = fqmul(a[i], b[i]);

  c[0] = d[0] + fqmul((CALC_D(a, b, 1, 7, d) + CALC_D(a, b, 2, 6, d) + CALC_D(a, b, 3, 5, d) + d[4]), zeta);
  c[1] = CALC_D(a, b, 0, 1, d) + fqmul((CALC_D(a, b, 2, 7, d) + CALC_D(a, b, 3, 6, d) + CALC_D(a, b, 4, 5, d)), zeta);
  c[2] = barrett_reduce(CALC_D(a, b, 0, 2, d) + d[1] + fqmul((CALC_D(a, b, 3, 7, d) + CALC_D(a, b, 4, 6, d) + d[5]), zeta));
  c[3] = barrett_reduce(CALC_D(a, b, 0, 3, d) + CALC_D(a, b, 1, 2, d) + fqmul((CALC_D(a, b, 4, 7, d) + CALC_D(a, b, 5, 6, d)), zeta));
  c[4] = barrett_reduce(CALC_D(a, b, 0, 4, d) + CALC_D(a, b, 1, 3, d) + d[2] + fqmul((CALC_D(a, b, 5, 7, d) + d[6]), zeta));
  c[5] = barrett_reduce(CALC_D(a, b, 0, 5, d) + CALC_D(a, b, 1, 4, d) + CALC_D(a, b, 2, 3, d)) + fqmul(CALC_D(a, b, 6, 7, d), zeta);
  c[6] = barrett_reduce(CALC_D(a, b, 0, 6, d) + CALC_D(a, b, 1, 5, d) + CALC_D(a, b, 2, 4, d)) + d[3] + fqmul(d[7], zeta);
  c[7] = barrett_reduce(CALC_D(a, b, 0, 7, d) + CALC_D(a, b, 1, 6, d) + CALC_D(a, b, 2, 5, d) + CALC_D(a, b, 3, 4, d));

#elif (CTRU_N == 1024)
  int16_t d[16];
  for (int i = 0; i < 16; i++)
    d[i] = fqmul(a[i], b[i]);

  c[0] = d[0] + fqmul((CALC_D(a, b, 1, 15, d) + CALC_D(a, b, 2, 14, d) + CALC_D(a, b, 3, 13, d) + CALC_D(a, b, 4, 12, d) + CALC_D(a, b, 5, 11, d) + CALC_D(a, b, 6, 10, d) + CALC_D(a, b, 7, 9, d) + d[8]), zeta);
  c[1] = CALC_D(a, b, 0, 1, d) + fqmul((CALC_D(a, b, 2, 15, d) + CALC_D(a, b, 3, 14, d) + CALC_D(a, b, 4, 13, d) + CALC_D(a, b, 5, 12, d) + CALC_D(a, b, 6, 11, d) + CALC_D(a, b, 7, 10, d) + CALC_D(a, b, 8, 9, d)), zeta);
  c[2] = barrett_reduce(CALC_D(a, b, 0, 2, d) + d[1] + fqmul((CALC_D(a, b, 3, 15, d) + CALC_D(a, b, 4, 14, d) + CALC_D(a, b, 5, 13, d) + CALC_D(a, b, 6, 12, d) + CALC_D(a, b, 7, 11, d) + CALC_D(a, b, 8, 10, d) + d[9]), zeta));
  c[3] = barrett_reduce(CALC_D(a, b, 0, 3, d) + CALC_D(a, b, 1, 2, d) + fqmul((CALC_D(a, b, 4, 15, d) + CALC_D(a, b, 5, 14, d) + CALC_D(a, b, 6, 13, d) + CALC_D(a, b, 7, 12, d) + CALC_D(a, b, 8, 11, d) + CALC_D(a, b, 9, 10, d)), zeta));
  c[4] = barrett_reduce(CALC_D(a, b, 0, 4, d) + CALC_D(a, b, 1, 3, d) + d[2] + fqmul((CALC_D(a, b, 5, 15, d) + CALC_D(a, b, 6, 14, d) + CALC_D(a, b, 7, 13, d) + CALC_D(a, b, 8, 12, d) + CALC_D(a, b, 9, 11, d) + d[10]), zeta));
  c[5] = barrett_reduce(CALC_D(a, b, 0, 5, d) + CALC_D(a, b, 1, 4, d) + CALC_D(a, b, 2, 3, d) + fqmul((CALC_D(a, b, 6, 15, d) + CALC_D(a, b, 7, 14, d) + CALC_D(a, b, 8, 13, d) + CALC_D(a, b, 9, 12, d) + CALC_D(a, b, 10, 11, d)), zeta));
  c[6] = barrett_reduce(CALC_D(a, b, 0, 6, d) + CALC_D(a, b, 1, 5, d) + CALC_D(a, b, 2, 4, d) + d[3] + fqmul((CALC_D(a, b, 7, 15, d) + CALC_D(a, b, 8, 14, d) + CALC_D(a, b, 9, 13, d) + CALC_D(a, b, 10, 12, d) + d[11]), zeta));
  c[7] = barrett_reduce(CALC_D(a, b, 0, 7, d) + CALC_D(a, b, 1, 6, d) + CALC_D(a, b, 2, 5, d) + CALC_D(a, b, 3, 4, d) + fqmul((CALC_D(a, b, 8, 15, d) + CALC_D(a, b, 9, 14, d) + CALC_D(a, b, 10, 13, d) + CALC_D(a, b, 11, 12, d)), zeta));
  c[8] = barrett_reduce(CALC_D(a, b, 0, 8, d) + CALC_D(a, b, 1, 7, d) + CALC_D(a, b, 2, 6, d) + CALC_D(a, b, 3, 5, d) + d[4] + fqmul((CALC_D(a, b, 9, 15, d) + CALC_D(a, b, 10, 14, d) + CALC_D(a, b, 11, 13, d) + d[12]), zeta));
  c[9] = barrett_reduce(CALC_D(a, b, 0, 9, d) + CALC_D(a, b, 1, 8, d) + CALC_D(a, b, 2, 7, d) + CALC_D(a, b, 3, 6, d) + CALC_D(a, b, 4, 5, d) + fqmul((CALC_D(a, b, 10, 15, d) + CALC_D(a, b, 11, 14, d) + CALC_D(a, b, 12, 13, d)), zeta));
  c[10] = barrett_reduce(CALC_D(a, b, 0, 10, d) + CALC_D(a, b, 1, 9, d) + CALC_D(a, b, 2, 8, d) + CALC_D(a, b, 3, 7, d) + CALC_D(a, b, 4, 6, d) + d[5] + fqmul((CALC_D(a, b, 11, 15, d) + CALC_D(a, b, 12, 14, d) + d[13]), zeta));
  c[11] = barrett_reduce(CALC_D(a, b, 0, 11, d) + CALC_D(a, b, 1, 10, d) + CALC_D(a, b, 2, 9, d) + CALC_D(a, b, 3, 8, d) + CALC_D(a, b, 4, 7, d) + CALC_D(a, b, 5, 6, d) + fqmul((CALC_D(a, b, 12, 15, d) + CALC_D(a, b, 13, 14, d)), zeta));
  c[12] = barrett_reduce(CALC_D(a, b, 0, 12, d) + CALC_D(a, b, 1, 11, d) + CALC_D(a, b, 2, 10, d) + CALC_D(a, b, 3, 9, d) + CALC_D(a, b, 4, 8, d) + CALC_D(a, b, 5, 7, d) + d[6] + fqmul((CALC_D(a, b, 13, 15, d) + d[14]), zeta));
  c[13] = barrett_reduce(CALC_D(a, b, 0, 13, d) + CALC_D(a, b, 1, 12, d) + CALC_D(a, b, 2, 11, d) + CALC_D(a, b, 3, 10, d) + CALC_D(a, b, 4, 9, d) + CALC_D(a, b, 5, 8, d) + CALC_D(a, b, 6, 7, d) + fqmul((CALC_D(a, b, 14, 15, d)), zeta));
  c[14] = barrett_reduce(CALC_D(a, b, 0, 14, d) + CALC_D(a, b, 1, 13, d) + CALC_D(a, b, 2, 12, d) + CALC_D(a, b, 3, 11, d) + CALC_D(a, b, 4, 10, d) + CALC_D(a, b, 5, 9, d) + CALC_D(a, b, 6, 8, d) + d[7] + fqmul((d[15]), zeta));
  c[15] = barrett_reduce(CALC_D(a, b, 0, 15, d) + CALC_D(a, b, 1, 14, d) + CALC_D(a, b, 2, 13, d) + CALC_D(a, b, 3, 12, d) + CALC_D(a, b, 4, 11, d) + CALC_D(a, b, 5, 10, d) + CALC_D(a, b, 6, 9, d) + CALC_D(a, b, 7, 8, d));
#endif
}
