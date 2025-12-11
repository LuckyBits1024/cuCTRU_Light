#include <stdint.h>
#include <stdio.h>
#include "poly.h"
#include "params.h"
#include "randombytes.h"
#include "reduce.h"
int16_t gen_zetas_withoutmont[64] = {-1, -62, 40, 173, 27, 136, -311, -57, 144, -300, -377, -304, -43, -359, 182, -251, -25, -12, 231, -289, -94, 324, -85, 113, -245, 190, -197, 90, -306, 253, -64, -123, -5, -310, 200, 96, 
135, -89, -17, -285, -49, 38, -347, 18, -215, -257, 141, 283, -125, -60, -383, 93, 299, 82, 344, -204, 313, 181, -216, -319, 8, -273, -320, 154};
int16_t gen_zetas[64] = {
    171, 605, 81, 408, 766, 583, 120, 519, 753, 546, 640, 461, 432, 638, 407, 626,
    430, 514, 487, 203, 694, 733, 693, 671, 369, 577, 620, 759, 34, 570, 178, 270,
    86, 718, 405, 502, 754, 608, 600, 288, 689, 423, 124, 767, 622, 114, 497, 54,
    612, 263, 128, 246, 394, 589, 389, 279, 307, 578, 24, 719, 170, 543, 121, 581};
void generate_zetas()
{
    for(int i = 0; i < 64; i ++)
    {
        int16_t zeta = gen_zetas[i];
        gen_zetas_withoutmont[i] = fqmul(-zeta, 1);
        printf("%d, ",gen_zetas_withoutmont[i]);
    }
    printf("\n");
}
void test_ntt()
{
    poly a, b, c, c2;
    unsigned char coin[CTRU_COINBYTES_KEYGEN];
    randombytes(coin, CTRU_COINBYTES_KEYGEN);
    poly_sample_enc(&a, coin);
    poly_sample_enc(&b, coin + CTRU_COINBYTES_KEYGEN / 2);
    poly_naivemul_q(&c2, &a, &b, CTRU_Q);

    poly_ntt(&a);
    poly_ntt(&b);
    poly_basemul(&c, &a, &b);
    poly_invntt(&c);
    for(int i = 0; i < CTRU_N; i ++)
    {
        if(c.coeffs[i] != c2.coeffs[i])
        {
          printf("error!\n");
          return;
        }
    }
}
void test_inverse()
{
    poly a, b, c;
    unsigned char coin[CTRU_COINBYTES_KEYGEN];
    randombytes(coin, CTRU_COINBYTES_KEYGEN);
    poly_sample_enc(&a, coin);
    
    

    poly_ntt(&a);
    poly_ntt(&b);
    poly_baseinv(&b, &a);
    poly_basemul(&c, &a, &b);
    poly_invntt(&c);
    for(int i = 0; i < CTRU_N; i ++)
    {
        printf("%d, ", c.coeffs[i]);
    }
    printf("\n");
}
int main()
{
    // test_ntt();
    // test_inverse();
    generate_zetas();
}
