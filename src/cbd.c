#include <stdint.h>
#include "params.h"
#include "cbd.h"

void cbd1(poly *r, const uint8_t buf[CTRU_CBD1_BYTES])
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
