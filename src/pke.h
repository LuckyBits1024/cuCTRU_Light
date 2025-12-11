#ifndef NTRU_H
#define NTRU_H

#include "poly.h"

int crypto_pke_keygen(unsigned char pk[CTRU_PKE_PUBLICKEYBYTES],
                      unsigned char sk[CTRU_PKE_SECRETKEYBYTES],
                      const unsigned char coins[CTRU_N]);

void crypto_pke_enc(unsigned char ct[CTRU_PKE_CIPHERTEXTBYTES],
                    const unsigned char pk[CTRU_PKE_PUBLICKEYBYTES],
                    const unsigned char m[CTRU_MSGBYTES],
                    const unsigned char coins[CTRU_ENC_COIN_BYTES]);

void crypto_pke_dec(unsigned char m[CTRU_MSGBYTES],
                    const unsigned char ct[CTRU_PKE_CIPHERTEXTBYTES],
                    const unsigned char sk[CTRU_PKE_SECRETKEYBYTES]);

#endif
