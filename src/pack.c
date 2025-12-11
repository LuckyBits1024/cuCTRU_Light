#include "pack.h"
#include "inverse.h"

/* Based on the reference implementation of NTRU Prime (NIST 3rd round submission)
 * by Daniel J. Bernstein, Chitchanok Chuengsatiansup, Tanja Lange, Christine van Vredendaal.
 * It can be used for q = 641.
 * */

uint16_t uint32_mod_uint14(uint32_t x, uint16_t m)
{
    uint32_t q;
    uint16_t r;
    uint32_divmod_uint14(&q, &r, x, m);
    return r;
}

void Encode(unsigned char *out, const uint16_t *R, const uint16_t *M, const long long len)
{
    if (len == 1)
    {
        uint16_t r = R[0];
        uint16_t m = M[0];
        while (m > 1)
        {
            *out++ = r;
            r >>= 8;
            m = (m + 255) >> 8;
        }
    }

    if (len > 1)
    {
        uint16_t R2[(len + 1) / 2];
        uint16_t M2[(len + 1) / 2];
        long long i;

        for (i = 0; i < len - 1; i += 2)
        {
            uint32_t m0 = M[i];
            uint32_t r = R[i] + R[i + 1] * m0;
            uint32_t m = M[i + 1] * m0;

            while (m >= 1024)
            {
                *out++ = r;
                r >>= 8;
                m = (m + 255) >> 8;
            }
            R2[i / 2] = r;
            M2[i / 2] = m;
        }

        if (i < len)
        {
            R2[i / 2] = R[i];
            M2[i / 2] = M[i];
        }

        Encode(out, R2, M2, (len + 1) / 2);
    }
}

void Decode(uint16_t *out, const unsigned char *S, const uint16_t *M, const long long len)
{
    if (len == 1)
    {
        if (M[0] == 1)
            *out = 0;
        else if (M[0] <= 256)
            *out = uint32_mod_uint14(S[0], M[0]);
        else
            *out = uint32_mod_uint14(S[0] + (((uint16_t)S[1]) << 8), M[0]);
    }

    if (len > 1)
    {
        uint16_t R2[(len + 1) / 2];
        uint16_t M2[(len + 1) / 2];
        uint16_t bottomr[len / 2];
        uint32_t bottomt[len / 2];
        long long i;
        for (i = 0; i < len - 1; i += 2)
        {
            uint32_t m = M[i] * (uint32_t)M[i + 1];
            if (m > 256 * 1023)
            {
                bottomt[i / 2] = 256 * 256;
                bottomr[i / 2] = S[0] + 256 * S[1];
                S += 2;
                M2[i / 2] = (((m + 255) >> 8) + 255) >> 8;
            }
            else if (m >= 1024)
            {
                bottomt[i / 2] = 256;
                bottomr[i / 2] = S[0];
                S += 1;
                M2[i / 2] = (m + 255) >> 8;
            }
            else
            {
                bottomt[i / 2] = 1;
                bottomr[i / 2] = 0;
                M2[i / 2] = m;
            }
        }
        if (i < len)
            M2[i / 2] = M[i];

        Decode(R2, S, M2, (len + 1) / 2);

        for (i = 0; i < len - 1; i += 2)
        {
            uint32_t r = bottomr[i / 2];
            uint32_t r1;
            uint16_t r0;
            r += bottomt[i / 2] * R2[i / 2];
            uint32_divmod_uint14(&r1, &r0, r, M[i]);
            r1 = uint32_mod_uint14(r1, M[i + 1]);
            *out++ = r0;
            *out++ = r1;
        }
        if (i < len)
            *out++ = R2[i / 2];
    }
}

 void pack_pk(unsigned char *r, const poly *a)
 {
     for(int i = 0; i < CTRU_N/4; i ++)
     {
         r[5*i+0] = a->coeffs[4*i+0] & 0xFF;
         r[5*i+1] = (a->coeffs[4*i+0] >> 8) | (a->coeffs[4*i+1] << 2);
         r[5*i+2] = (a->coeffs[4*i+1] >> 6) | (a->coeffs[4*i+2] << 4);
         r[5*i+3] = (a->coeffs[4*i+2] >> 4) | (a->coeffs[4*i+3] << 6);
         r[5*i+4] = (a->coeffs[4*i+3] >> 2) & 0xFF;
     }
 }
 void unpack_pk(poly *r, const unsigned char *a)
 {
     for(int i = 0; i < CTRU_N/4; i ++)
     {
         r->coeffs[4*i+0] =((a[5*i+0] >> 0) | (a[5*i+1] << 8)) & 0x3FF;
         r->coeffs[4*i+1] =((a[5*i+1] >> 2) | (a[5*i+2] << 6)) & 0x3FF;
         r->coeffs[4*i+2] =((a[5*i+2] >> 4) | (a[5*i+3] << 4)) & 0x3FF;
         r->coeffs[4*i+3] =((a[5*i+3] >> 6) | (a[5*i+4] << 2)) & 0x3FF;
     }
 }
//void pack_pk(unsigned char *r, const poly *a)
//{
//    uint16_t R[CTRU_N], M[CTRU_N];
//
//    for (int i = 0; i < CTRU_N; ++i)
//        R[i] = (uint16_t)a->coeffs[i];
//
//    for (int i = 0; i < CTRU_N; ++i)
//        M[i] = CTRU_Q;
//
//    Encode(r, R, M, CTRU_N);
//}
//
//void unpack_pk(poly *r, const unsigned char *a)
//{
//    uint16_t R[CTRU_N], M[CTRU_N];
//
//    for (int i = 0; i < CTRU_N; ++i)
//        M[i] = CTRU_Q;
//
//    Decode(R, a, M, CTRU_N);
//
//    for (int i = 0; i < CTRU_N; ++i)
//        r->coeffs[i] = ((int16_t)R[i]);
//}

void pack_sk(unsigned char *r, const poly *a)
{
    unsigned int i;
    for (i = 0; i < CTRU_N / 4; ++i)
    {
        r[5 * i + 0] = (a->coeffs[4 * i + 0] >> 0);
        r[5 * i + 1] = (a->coeffs[4 * i + 0] >> 8) | ((int16_t)a->coeffs[4 * i + 1] << 2);
        r[5 * i + 2] = (a->coeffs[4 * i + 1] >> 6) | ((int16_t)a->coeffs[4 * i + 2] << 4);
        r[5 * i + 3] = (a->coeffs[4 * i + 2] >> 4) | ((int16_t)a->coeffs[4 * i + 3] << 6);
        r[5 * i + 4] = (a->coeffs[4 * i + 3] >> 2);
    }
}

void unpack_sk(poly *r, const unsigned char *a)
{
    unsigned int i;
    for (i = 0; i < CTRU_N / 4; ++i)
    {
        r->coeffs[4 * i + 0] = ((a[5 * i + 0] >> 0) | ((uint16_t)a[5 * i + 1] << 8)) & 0x3FF;
        r->coeffs[4 * i + 1] = ((a[5 * i + 1] >> 2) | ((uint16_t)a[5 * i + 2] << 6)) & 0x3FF;
        r->coeffs[4 * i + 2] = ((a[5 * i + 2] >> 4) | ((uint16_t)a[5 * i + 3] << 4)) & 0x3FF;
        r->coeffs[4 * i + 3] = ((a[5 * i + 3] >> 6) | ((uint16_t)a[5 * i + 4] << 2)) & 0x3FF;
    }
}

#if (KEM_TYPE == RLWE_KEM)
void pack_ct(unsigned char *r, const poly *a)
{
    uint16_t R[CTRU_N], M[CTRU_N];

    for (int i = 0; i < CTRU_N; ++i)
        R[i] = (uint16_t)a->coeffs[i];

    for (int i = 0; i < CTRU_N; ++i)
        M[i] = CTRU_Q;

    Encode(r, R, M, CTRU_N);
}

void unpack_ct(poly *r, const unsigned char *a)
{
    uint16_t R[CTRU_N], M[CTRU_N];

    for (int i = 0; i < CTRU_N; ++i)
        M[i] = CTRU_Q;

    Decode(R, a, M, CTRU_N);

    for (int i = 0; i < CTRU_N; ++i)
        r->coeffs[i] = ((int16_t)R[i]);
}

#elif (KEM_TYPE == RLWR_KEM)
void pack_ct(unsigned char *r, const poly *a)
{
    unsigned int i;
#if (CTRU_Q2 == 256)
    for (i = 0; i < CTRU_N; ++i)
    {
        r[i] = (uint8_t)a->coeffs[i];
    }
#elif (CTRU_Q2 == 512)
    for (i = 0; i < CTRU_N / 8; ++i)
    {
        r[9 * i + 0] = a->coeffs[8 * i + 0] >> 0;
        r[9 * i + 1] = a->coeffs[8 * i + 0] >> 8 | ((int16_t)a->coeffs[8 * i + 1] << 1);
        r[9 * i + 2] = a->coeffs[8 * i + 1] >> 7 | ((int16_t)a->coeffs[8 * i + 2] << 2);
        r[9 * i + 3] = a->coeffs[8 * i + 2] >> 6 | ((int16_t)a->coeffs[8 * i + 3] << 3);
        r[9 * i + 4] = a->coeffs[8 * i + 3] >> 5 | ((int16_t)a->coeffs[8 * i + 4] << 4);
        r[9 * i + 5] = a->coeffs[8 * i + 4] >> 4 | ((int16_t)a->coeffs[8 * i + 5] << 5);
        r[9 * i + 6] = a->coeffs[8 * i + 5] >> 3 | ((int16_t)a->coeffs[8 * i + 6] << 6);
        r[9 * i + 7] = a->coeffs[8 * i + 6] >> 2 | ((int16_t)a->coeffs[8 * i + 7] << 7);
        r[9 * i + 8] = a->coeffs[8 * i + 7] >> 1;
    }
#endif
}

void unpack_ct(poly *r, const unsigned char *a)
{
    unsigned int i;
#if (CTRU_Q2 == 256)
    for (i = 0; i < CTRU_N; ++i)
    {
        r->coeffs[i] = (int16_t)a[i];
    }
#elif (CTRU_Q2 == 512)
    for (i = 0; i < CTRU_N / 8; ++i)
    {
        r->coeffs[8 * i + 0] = ((a[9 * i + 0] >> 0) | ((uint16_t)a[9 * i + 1] << 8)) & 0x1FF;
        r->coeffs[8 * i + 1] = ((a[9 * i + 1] >> 1) | ((uint16_t)a[9 * i + 2] << 7)) & 0x1FF;
        r->coeffs[8 * i + 2] = ((a[9 * i + 2] >> 2) | ((uint16_t)a[9 * i + 3] << 6)) & 0x1FF;
        r->coeffs[8 * i + 3] = ((a[9 * i + 3] >> 3) | ((uint16_t)a[9 * i + 4] << 5)) & 0x1FF;
        r->coeffs[8 * i + 4] = ((a[9 * i + 4] >> 4) | ((uint16_t)a[9 * i + 5] << 4)) & 0x1FF;
        r->coeffs[8 * i + 5] = ((a[9 * i + 5] >> 5) | ((uint16_t)a[9 * i + 6] << 3)) & 0x1FF;
        r->coeffs[8 * i + 6] = ((a[9 * i + 6] >> 6) | ((uint16_t)a[9 * i + 7] << 2)) & 0x1FF;
        r->coeffs[8 * i + 7] = ((a[9 * i + 7] >> 7) | ((uint16_t)a[9 * i + 8] << 1)) & 0x1FF;
    }
#endif
}
#endif
