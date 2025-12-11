#ifndef PARAMS_H
#define PARAMS_H

#define RLWE_KEM 1
#define RLWR_KEM 2

#ifndef KEM_TYPE          /* which type: RLWE or RLWR*/
//#define KEM_TYPE RLWE_KEM /* It is RLWE_KEM */
#define KEM_TYPE RLWR_KEM /* It is RLWR_KEM */
#endif

#define CTRU_N 512 /* Change this for different security strengths */
// #define CTRU_N 1024

#ifndef USING_SHA2_VARIANT
#define USING_SHA2_VARIANT 0 /* if using sha3 */
#endif

#ifndef USING_NEV_FO
#define USING_NEV_FO 0 /* if using CTRU's FO */
#endif

// #define CTRU_Q 641
#define CTRU_Q 769//这里是选不同的环
#define CTRU_LOGQ 10

#if (KEM_TYPE == RLWE_KEM)
#if (CTRU_N == 512)
#define CTRU_KEYGEN_COIN_BYTES (CTRU_N / 2)
#define CTRU_ENC_COIN_BYTES (CTRU_N / 2)
#define CTRU_PKE_PUBLICKEYBYTES (640)
#define CTRU_PKE_CIPHERTEXTBYTES (603)
#elif (CTRU_N == 1024)
#define CTRU_KEYGEN_COIN_BYTES ((4 * CTRU_N / 3) + 100)
#define CTRU_ENC_COIN_BYTES ((4 * CTRU_N / 3) + 100)
#define CTRU_PKE_PUBLICKEYBYTES (1206)
#define CTRU_PKE_CIPHERTEXTBYTES (1206)
#endif

#elif (KEM_TYPE == RLWR_KEM)
#if (CTRU_N == 512)
#define CTRU_Q2 256  /* Change this for the ciphertext modulus */
#define CTRU_LOGQ2 8 /* Change this for the ciphertext modulus */
#define CTRU_KEYGEN_COIN_BYTES (CTRU_N / 2)
#define CTRU_ENC_COIN_BYTES (CTRU_N / 4)
#define CTRU_PKE_PUBLICKEYBYTES 512*10/8
#define CTRU_PKE_CIPHERTEXTBYTES (CTRU_N)
#elif (CTRU_N == 1024)
#define CTRU_Q2 512
#define CTRU_LOGQ2 9
#define CTRU_KEYGEN_COIN_BYTES ((16 * CTRU_N / 15) + 100)
#define CTRU_ENC_COIN_BYTES (CTRU_N / 4)
#define CTRU_PKE_PUBLICKEYBYTES (1206)
#define CTRU_PKE_CIPHERTEXTBYTES (9 * CTRU_N / 8)
#endif
#endif

#define CTRU_SEEDBYTES 32
#define CTRU_SHAREDKEYBYTES 32
#define CTRU_MSGBYTES 32
#define CTRU_PREFIXHASHBYTES 33

#define CTRU_PKE_SECRETKEYBYTES (10 * CTRU_N / 8)

#if (USING_NEV_FO == 0)
#define CTRU_KEM_PUBLICKEYBYTES CTRU_PKE_PUBLICKEYBYTES
#define CTRU_KEM_SECRETKEYBYTES (CTRU_PKE_SECRETKEYBYTES + CTRU_PKE_PUBLICKEYBYTES + CTRU_SEEDBYTES)
#define CTRU_KEM_CIPHERTEXTBYTES CTRU_PKE_CIPHERTEXTBYTES
#endif

#endif