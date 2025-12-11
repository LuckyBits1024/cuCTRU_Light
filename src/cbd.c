#include <stdint.h>
#include "params.h"
#include "cbd.h"

void cbd1(poly *r, const uint8_t buf[])
{
  int i;
  for (i = 0; i < CTRU_N / 4; i++)
  {
    r->coeffs[4 * i + 0] = ((buf[i] >> 0) & 1) - ((buf[i] >> 1) & 1);
    r->coeffs[4 * i + 1] = ((buf[i] >> 2) & 1) - ((buf[i] >> 3) & 1);
    r->coeffs[4 * i + 2] = ((buf[i] >> 4) & 1) - ((buf[i] >> 5) & 1);
    r->coeffs[4 * i + 3] = ((buf[i] >> 6) & 1) - ((buf[i] >> 7) & 1);
  }
}

void ternary_1of5(poly *r, const uint8_t buf[])
{
  int count = 0, offset = 0;

  while (count < CTRU_N)
  {
    uint16_t tmp = buf[offset] & 0xF;
    if (tmp >= 15)
    {
    }
    else if (tmp <= 2)
    {
      r->coeffs[count++] = -1;
    }
    else if (tmp >= 12)
    {
      r->coeffs[count++] = 1;
    }
    else
    {
      r->coeffs[count++] = 0;
    }

    tmp = (buf[offset++] >> 4) & 0xF;
    if (tmp >= 15)
    {
    }
    else if (tmp <= 2)
    {
      r->coeffs[count++] = -1;
    }
    else if (tmp >= 12)
    {
      r->coeffs[count++] = 1;
    }
    else
    {
      r->coeffs[count++] = 0;
    }
  }
}

void ternary_1of6(poly *r, const uint8_t buf[])
{
  int count = 0, offset = 0;

  while (count < CTRU_N)
  {
    uint16_t tmp = buf[offset] & 0xF;
    if (tmp >= 12)
    {
    }
    else if (tmp <= 1)
    {
      r->coeffs[count++] = -1;
    }
    else if (tmp >= 10)
    {
      r->coeffs[count++] = 1;
    }
    else
    {
      r->coeffs[count++] = 0;
    }

    tmp = (buf[offset++] >> 4) & 0xF;
    if (tmp >= 12)
    {
    }
    else if (tmp <= 1)
    {
      r->coeffs[count++] = -1;
    }
    else if (tmp >= 10)
    {
      r->coeffs[count++] = 1;
    }
    else
    {
      r->coeffs[count++] = 0;
    }
  }
}
