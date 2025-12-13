#include <stdint.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include "params.h"
#include "ctru.h"
#include "test_vector.h"

// CPU 版 keygen 测试：逻辑上对应 GPU 的 test_ctru.cu，
// 使用固定的 kFixedNTTInput 作为 coins，调用 pke_keygen，然后打印 sk 和 pk。
int main(void)
{
    unsigned char pk[CTRU_PKE_PUBLICKEYBYTES];
    unsigned char sk[CTRU_PKE_SECRETKEYBYTES];
    unsigned char coins[CTRU_COINBYTES_KEYGEN];

    // 用固定测试向量初始化 coins（按字节截断）
    // GPU 版本里虽然把 coins 复制到设备，但 kernel 内部会自己生成随机种子；
    // 这里我们显式用固定 coins 来让 CPU 版 keygen 可重复。
    for (int i = 0; i < CTRU_COINBYTES_KEYGEN; i++) {
        coins[i] = (unsigned char)(kFixedNTTInput[i % CTRU_N] & 0xFF);
    }

    // 调用 CPU 版 PKE KeyGen
    pke_keygen(pk, sk, coins);

    printf("生成的私钥为\n");
    for (int i = 0; i < CTRU_PKE_SECRETKEYBYTES; i++) {
        printf("%d ", sk[i]);
        if ((i + 1) % 32 == 0)
            printf("\n");
    }
    printf("\n");

    printf("生成的公钥为\n");
    for (int i = 0; i < CTRU_PKE_PUBLICKEYBYTES; i++) {
        printf("%d ", pk[i]);
        if ((i + 1) % 32 == 0)
            printf("\n");
    }
    printf("\n");

    return 0;
}