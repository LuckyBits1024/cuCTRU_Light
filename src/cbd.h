#ifndef CBD_H
#define CBD_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

void cbd1(poly *r, const uint8_t buf[]);
void ternary_1of5(poly *r, const uint8_t buf[]);
void ternary_1of6(poly *r, const uint8_t buf[]);

#endif
