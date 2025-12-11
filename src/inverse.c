#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "params.h"
#include "inverse.h"
#include "reduce.h"
#include "poly.h"
#include "ntt.h"

/* Based on the reference implementation of NTRU Prime (NIST 3rd round submission)
 * by Daniel J. Bernstein, Chitchanok Chuengsatiansup, Tanja Lange, Christine van Vredendaal.
 * It can be used for q = 641.
 * */

void uint32_divmod_uint14(uint32_t *y, uint16_t *r, uint32_t x, uint16_t m)
{
    uint32_t w = 0x80000000;
    uint32_t qpart;
    uint32_t mask;

    w /= m;

    *y = 0;
    qpart = (x * (uint64_t)w) >> 31;
    x -= qpart * m;
    *y += qpart;

    qpart = (x * (uint64_t)w) >> 31;
    x -= qpart * m;
    *y += qpart;

    x -= m;
    *y += 1;
    mask = -(x >> 31);
    x += mask & (uint32_t)m;
    *y += mask;

    *r = x;
}

void int32_divmod_uint14(int32_t *y, uint16_t *r, int32_t x, uint16_t m)
{
    uint32_t uq, uq2;
    uint16_t ur, ur2;
    uint32_t mask;

    uint32_divmod_uint14(&uq, &ur, 0x80000000 + (uint32_t)x, m);
    uint32_divmod_uint14(&uq2, &ur2, 0x80000000, m);

    ur -= ur2;
    uq -= uq2;

    mask = -(uint32_t)(ur >> 15);
    ur += mask & m;
    uq += mask;
    *r = ur;
    *y = uq;
}

uint16_t int32_mod_uint14(int32_t x, uint16_t m)
{
    int32_t y;
    uint16_t r;
    int32_divmod_uint14(&y, &r, x, m);
    return r;
}

int16_t fq_freeze(int32_t x)
{

    const int16_t half_q = (CTRU_Q - 1) >> 1;
    return int32_mod_uint14(x + half_q, CTRU_Q) - half_q;
}

int int16_nonzero_mask(int16_t x)
{
    uint16_t u = x;
    uint32_t w = u;
    w = -w;
    w >>= 31;
    return -w;
}

int int16_negative_mask(int16_t x)
{
    uint16_t u = x;
    u >>= 15;
    return -(int)u;
}
int rq_inverse_opt(int16_t finv[ROOT_DIMENSION], const int16_t f[ROOT_DIMENSION], const int16_t zeta)
{
    int16_t Phi[ROOT_DIMENSION + 1] = {1, 0, 0, 0, 0, 0, 0, 0, zeta};
    int16_t V[ROOT_DIMENSION + 1] = {0};
    int16_t S[ROOT_DIMENSION + 1] = {1, 0, 0, 0, 0, 0, 0, 0, 0};
    int16_t F[ROOT_DIMENSION + 1] = {f[7], f[6], f[5], f[4], f[3], f[2], f[1], f[0], 0};
    int i, loop, swap, t;
    int Delta = 1;
    int32_t Phi0, F0;
    int16_t scale;

    // for (i = 0; i < ROOT_DIMENSION; ++i) {
    //     Phi[i] = 0;
    //     F[ROOT_DIMENSION - 1 - i] = f[i];
    //     V[i] = 0;
    //     S[i] = 0;
    // }
    // Phi[0] = 1;
    // Phi[ROOT_DIMENSION] = zeta;
    // F[ROOT_DIMENSION] = 0;

    // Delta = 1;

    // V[ROOT_DIMENSION] = 0;
    // S[ROOT_DIMENSION] = 0;
    // S[0] = 1;



//    for (i = 0; i < ROOT_DIMENSION; ++i)
//        Phi[i] = 0;
//    Phi[0] = 1;
//
//    Phi[ROOT_DIMENSION] = zeta;
//
//    for (i = 0; i < ROOT_DIMENSION; ++i)
//        F[ROOT_DIMENSION - 1 - i] = f[i];
//    F[ROOT_DIMENSION] = 0;
//
//    Delta = 1;
//
//    for (i = 0; i < ROOT_DIMENSION + 1; ++i)
//        V[i] = 0;
//
//    for (i = 0; i < ROOT_DIMENSION + 1; ++i)
//        S[i] = 0;
//    S[0] = 1;

    for (loop = 0; loop < 2 * ROOT_DIMENSION - 1; ++loop)
    {
        for (i = ROOT_DIMENSION; i > 0; --i)
            V[i] = V[i - 1];
        V[0] = 0;

        swap = int16_negative_mask(-Delta) & int16_nonzero_mask(F[0]);

        for (i = 0; i < ROOT_DIMENSION + 1; ++i)
        {
            t = swap & (Phi[i] ^ F[i]);
            Phi[i] ^= t;
            F[i] ^= t;
            t = swap & (V[i] ^ S[i]);
            V[i] ^= t;
            S[i] ^= t;
        }

        Delta ^= swap & (Delta ^ -Delta);
        Delta++;

        Phi0 = Phi[0];
        F0 = F[0];

        for (i = 0; i < ROOT_DIMENSION + 1; ++i)
            F[i] = fq_freeze(Phi0 * F[i] - F0 * Phi[i]);
        for (i = 0; i < ROOT_DIMENSION; ++i)
            F[i] = F[i + 1];
        F[ROOT_DIMENSION] = 0;

        for (i = 0; i < ROOT_DIMENSION + 1; ++i)
            S[i] = fq_freeze(Phi0 * S[i] - F0 * V[i]);
    }

    scale = Phi[0];
    scale += (scale >> 15) & CTRU_Q;
    scale = fq_inverse_table[scale];
    for (i = 0; i < ROOT_DIMENSION; ++i)
        finv[i] = fq_freeze(scale * (int32_t)V[ROOT_DIMENSION - 1 - i]);

    return int16_nonzero_mask(Delta);
}
int rq_inverse(int16_t finv[ROOT_DIMENSION], const int16_t f[ROOT_DIMENSION], const int16_t zeta)
{
    int16_t Phi[ROOT_DIMENSION + 1], F[ROOT_DIMENSION + 1], V[ROOT_DIMENSION + 1], S[ROOT_DIMENSION + 1];
    int i, loop, Delta, swap, t;
    int32_t Phi0, F0;
    int16_t scale;

    for (i = 0; i < ROOT_DIMENSION; ++i)
        Phi[i] = 0;
    Phi[0] = 1;
    Phi[ROOT_DIMENSION] = fqmul(-zeta, 1);


    for (i = 0; i < ROOT_DIMENSION; ++i)
        F[ROOT_DIMENSION - 1 - i] = f[i];
    F[ROOT_DIMENSION] = 0;

    Delta = 1;

    for (i = 0; i < ROOT_DIMENSION + 1; ++i)
        V[i] = 0;

    for (i = 0; i < ROOT_DIMENSION + 1; ++i)
        S[i] = 0;
    S[0] = 1;
    
    for (loop = 0; loop < 2 * ROOT_DIMENSION - 1; ++loop)
    {
        for (i = ROOT_DIMENSION; i > 0; --i)
            V[i] = V[i - 1];
        V[0] = 0;

        swap = int16_negative_mask(-Delta) & int16_nonzero_mask(F[0]);

        for (i = 0; i < ROOT_DIMENSION + 1; ++i)
        {
            t = swap & (Phi[i] ^ F[i]);
            Phi[i] ^= t;
            F[i] ^= t;
            t = swap & (V[i] ^ S[i]);
            V[i] ^= t;
            S[i] ^= t;
        }

        Delta ^= swap & (Delta ^ -Delta);
        Delta++;

        Phi0 = Phi[0];
        F0 = F[0];

        for (i = 0; i < ROOT_DIMENSION + 1; ++i)
            F[i] = fq_freeze(Phi0 * F[i] - F0 * Phi[i]);
        for (i = 0; i < ROOT_DIMENSION; ++i)
            F[i] = F[i + 1];
        F[ROOT_DIMENSION] = 0;

        for (i = 0; i < ROOT_DIMENSION + 1; ++i)
            S[i] = fq_freeze(Phi0 * S[i] - F0 * V[i]);
    }

    scale = Phi[0];
    scale += (scale >> 15) & CTRU_Q;
    scale = fq_inverse_table[scale];
    for (i = 0; i < ROOT_DIMENSION; ++i)
        finv[i] = fq_freeze(scale * (int32_t)V[ROOT_DIMENSION - 1 - i]);

    return int16_nonzero_mask(Delta);
}