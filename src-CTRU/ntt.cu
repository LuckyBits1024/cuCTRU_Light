#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "params.h"
#include "reduce.cuh"
#include "ntt.cuh"

//test
// int16_t zetas[64] = {
//     154, 640, 318, 256, 625, 100, 250, 40, 631, 383, 637, 25, 359, 160, 400, 64,
//     609, 200, 500, 80, 636, 512, 639, 333, 77, 320, 159, 128, 633, 50, 125, 20,
//     29, 620, 268, 248, 305, 177, 122, 199, 431, 351, 557, 525, 488, 155, 67, 62,
//     610, 354, 244, 398, 536, 496, 599, 583, 335, 310, 134, 124, 473, 409, 61, 420};

// int16_t zetas_inv[64] = {
//     -420, -61, -409, -473, -124, -134, -310, -335, -583, -599, -496, -536, -398, -244, -354, -610,
//     -62, -67, -155, -488, -525, -557, -351, -431, -199, -122, -177, -305, -248, -268, -620, -29,
//     -20, -125, -50, -633, -128, -159, -320, -77, -333, -639, -512, -636, -80, -500, -200, -609,
//     -64, -400, -160, -359, -25, -637, -383, -631, -40, -250, -100, -625, -256, -318, -640, 10};
__device__ int16_t zetas[64] = {
    171, 605, 81, 408, 766, 583, 120, 519, 753, 546, 640, 461, 432, 638, 407, 626,
    430, 514, 487, 203, 694, 733, 693, 671, 369, 577, 620, 759, 34, 570, 178, 270,
    86, 718, 405, 502, 754, 608, 600, 288, 689, 423, 124, 767, 622, 114, 497, 54,
    612, 263, 128, 246, 394, 589, 389, 279, 307, 578, 24, 719, 170, 543, 121, 581};

__device__ int16_t zetas_inv[64] = {
    -581, -121, -543, -170, -719, -24, -578, -307, -279, -389, -589, -394, -246, -128, -263, -612,
    -54, -497, -114, -622, -767, -124, -423, -689, -288, -600, -608, -754, -502, -405, -718, -86,
    -270, -178, -570, -34, -759, -620, -577, -369, -671, -693, -733, -694, -203, -487, -514, -430,
    -626, -407, -638, -432, -461, -640, -546, -753, -519, -120, -583, -766, -408, -81, -605, 541};  // 64^(-1) mod 769

// CT

__device__  void ntt_butt(int16_t &a, int16_t &b, const int16_t zeta) {
    int16_t t = fqmul(zeta, b);
    b = a - t;
    a = a + t;
}
//s_ntt 共享内存
__device__ void ntt(int16_t regs[8], int16_t *s_ntt){

    // level 1 取下标系数为256的两个数进行蝴蝶变换
    ntt_butt(regs[0], regs[4], zetas[1]);
    ntt_butt(regs[1], regs[5], zetas[1]);
    ntt_butt(regs[2], regs[6], zetas[1]);
    ntt_butt(regs[3], regs[7], zetas[1]);

    // level 2  128
    ntt_butt(regs[0], regs[2], zetas[2]);
    ntt_butt(regs[1], regs[3], zetas[2]);
    ntt_butt(regs[4], regs[6], zetas[3]);
    ntt_butt(regs[5], regs[7], zetas[3]);

    // level 3   64
    ntt_butt(regs[0], regs[1], zetas[4]);
    ntt_butt(regs[2], regs[3], zetas[5]);
    ntt_butt(regs[4], regs[5], zetas[6]);
    ntt_butt(regs[6], regs[7], zetas[7]);

    // SMEM exchange
    #pragma unroll
    for (size_t i = 0; i < 8; i++)
        s_ntt[i * 64 + threadIdx.x] = regs[i]; //以64为间隔进行摆放
    __syncwarp();

    #pragma unroll
    for (size_t i = 0; i < 8; i++)
    //每8个线程对应64个系数，8个线程组成一个线程快。因此threadIdx.x/8用于确定位于哪一个数据块中
    // threadIdx.x & 7定位该线程在线程快内的偏移量，i*4定位是该线程中的哪一个数
        regs[i] = s_ntt[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8]; 
    __syncwarp();

    // level 4
    ntt_butt(regs[0], regs[4], zetas[8 + threadIdx.x / 8]);
    ntt_butt(regs[1], regs[5], zetas[8 + threadIdx.x / 8]);
    ntt_butt(regs[2], regs[6], zetas[8 + threadIdx.x / 8]);
    ntt_butt(regs[3], regs[7], zetas[8 + threadIdx.x / 8]);

    // level 5
    ntt_butt(regs[0], regs[2], zetas[16 + (threadIdx.x / 8) * 2]);
    ntt_butt(regs[1], regs[3], zetas[16 + (threadIdx.x / 8) * 2]);
    ntt_butt(regs[4], regs[6], zetas[17 + (threadIdx.x / 8) * 2]);
    ntt_butt(regs[5], regs[7], zetas[17 + (threadIdx.x / 8) * 2]);

    // level 6
    ntt_butt(regs[0], regs[1], zetas[32 + (threadIdx.x / 8) * 4]);
    ntt_butt(regs[2], regs[3], zetas[33 + (threadIdx.x / 8) * 4]);
    ntt_butt(regs[4], regs[5], zetas[34 + (threadIdx.x / 8) * 4]);
    ntt_butt(regs[6], regs[7], zetas[35 + (threadIdx.x / 8) * 4]);

// #pragma unroll
//     for (size_t i = 0; i < 8; i++)
//         //怎么读出来的，就怎么写回去
//         s_ntt[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8] = regs[i];
//     __syncwarp();
// #pragma unroll
//     for (size_t i = 0; i < 8; i++)
//         regs[i] = s_ntt[threadIdx.x * 8 + i]; //每个线程取8个系数即可

}
__device__ void inv_ntt_butt(int16_t &a, int16_t &b, int16_t zeta){
    int16_t t = a;
    a = barrett_reduce(t + b);
    b = barrett_reduce(t - b);
    b = fqmul(zeta, b);
}

//invntt
__device__ void invntt(int16_t regs[8], int16_t *s_ntt){
    // level 6
    inv_ntt_butt(regs[0], regs[1], zetas_inv[32 + (threadIdx.x / 8) * 4]);
    inv_ntt_butt(regs[2], regs[3], zetas_inv[33 + (threadIdx.x / 8) * 4]);
    inv_ntt_butt(regs[4], regs[5], zetas_inv[34 + (threadIdx.x / 8) * 4]);
    inv_ntt_butt(regs[6], regs[7], zetas_inv[35 + (threadIdx.x / 8) * 4]);

    // level 5
    inv_ntt_butt(regs[0], regs[2], zetas_inv[16 + (threadIdx.x / 8) * 2]);
    inv_ntt_butt(regs[1], regs[3], zetas_inv[16 + (threadIdx.x / 8) * 2]);
    inv_ntt_butt(regs[4], regs[6], zetas_inv[17 + (threadIdx.x / 8) * 2]);
    inv_ntt_butt(regs[5], regs[7], zetas_inv[17 + (threadIdx.x / 8) * 2]);

    // level 4
    inv_ntt_butt(regs[0], regs[4], zetas_inv[8 + threadIdx.x / 8]);
    inv_ntt_butt(regs[1], regs[5], zetas_inv[8 + threadIdx.x / 8]);
    inv_ntt_butt(regs[2], regs[6], zetas_inv[8 + threadIdx.x / 8]);
    inv_ntt_butt(regs[3], regs[7], zetas_inv[8 + threadIdx.x / 8]);
    #pragma unroll
    for (size_t i = 0; i < 8; i++)
    //每8个线程对应64个系数，8个线程组成一个线程快。因此threadIdx.x/8用于确定位于哪一个数据块中
    // threadIdx.x & 7定位该线程在线程快内的偏移量，i*4定位是该线程中的哪一个数
        s_ntt[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8] = regs[i]; 
    __syncwarp();

    // SMEM exchange
    #pragma unroll
    for (size_t i = 0; i < 8; i++)
         regs[i] = s_ntt[i * 64 + threadIdx.x]; //以64为间隔进行摆放
    __syncwarp();

    // level 3   64
    inv_ntt_butt(regs[0], regs[1], zetas_inv[4]);
    inv_ntt_butt(regs[2], regs[3], zetas_inv[5]);
    inv_ntt_butt(regs[4], regs[5], zetas_inv[6]);
    inv_ntt_butt(regs[6], regs[7], zetas_inv[7]);

    // level 2  128
    inv_ntt_butt(regs[0], regs[2], zetas_inv[2]);
    inv_ntt_butt(regs[1], regs[3], zetas_inv[2]);
    inv_ntt_butt(regs[4], regs[6], zetas_inv[3]);
    inv_ntt_butt(regs[5], regs[7], zetas_inv[3]);

    // level 1 取下标系数为256的两个数进行蝴蝶变换
    inv_ntt_butt(regs[0], regs[4], zetas_inv[1]);
    inv_ntt_butt(regs[1], regs[5], zetas_inv[1]);
    inv_ntt_butt(regs[2], regs[6], zetas_inv[1]);
    inv_ntt_butt(regs[3], regs[7], zetas_inv[1]);

}

__global__ void test_ntt(int16_t *g_polyvec){
    __shared__ int16_t s_ntt[CTRU_N];
    int16_t regs[8];

    for (size_t i = 0; i < 8; ++i){
        regs[i] = g_polyvec[64 * i + threadIdx.x];
        // printf("线程 %d  寄存器的值为 %d \n",threadIdx.x,regs[i]);
    }

    //写回到global mem 
    ntt(regs, s_ntt);
    for (size_t i = 0; i < 8; ++i){
        g_polyvec[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8] = regs[i];
    }

}

__global__ void test_inv_ntt(int16_t *g_polyvec){
    __shared__ int16_t s_ntt[CTRU_N];
    int16_t regs[8];

    for(size_t i = 0; i < 8; ++i)
        regs[i] = g_polyvec[(threadIdx.x / 8) * 64 + (threadIdx.x & 7) + i * 8];

    invntt(regs, s_ntt);

    for(size_t i = 0; i < 8; ++i){
        g_polyvec[64 * i + threadIdx.x] = regs[i];
    }
}

int main(){
    // int16_t *h_polyvec, *g_polyvec;
    int16_t  *g_polyvec , *inv_g_polyvec;
    size_t size = CTRU_N * sizeof(int16_t);
    // h_polyvec = (int16_t *)malloc(size);
    cudaMalloc((void**)&g_polyvec, size);
    srand(time(NULL));
    // 指定范围
    int min = 0;
    int max = CTRU_Q;

    // for(int i = 0; i < CTRU_N; ++i){
    //     h_polyvec[i] = min + rand() % (max - min + 1);
    // }
    // printf("生成的随机数组为\n");
    // for(int i = 0; i < CTRU_N; ++i){
    //     printf("%d, ", h_polyvec[i]);
    //     if((i + 1) % 16 == 0) printf("\n");
    // }

    int16_t h_polyvec[CTRU_N] = {
        494, 675, 181, 490, 588, 147, 140, 252, 481, 713, 318, 146, 579, 731, 287, 481, 
        635, 155, 93, 611, 587, 343, 570, 676, 765, 202, 264, 28, 727, 680, 468, 451, 
        585, 649, 323, 555, 179, 464, 37, 42, 559, 355, 188, 520, 316, 476, 231, 181, 
        13, 476, 174, 752, 49, 127, 658, 44, 329, 152, 72, 286, 214, 692, 119, 181, 
        572, 443, 736, 133, 289, 3, 175, 78, 511, 363, 750, 209, 221, 212, 543, 386, 
        688, 99, 369, 738, 226, 409, 164, 555, 714, 237, 223, 310, 311, 343, 492, 113, 
        168, 610, 246, 457, 766, 421, 687, 507, 167, 667, 98, 540, 109, 641, 157, 180, 
        741, 678, 148, 197, 317, 312, 135, 413, 701, 510, 724, 243, 235, 598, 508, 403, 
        438, 137, 242, 586, 558, 159, 475, 107, 209, 574, 648, 470, 445, 35, 650, 568, 
        713, 180, 148, 412, 493, 283, 208, 576, 175, 314, 201, 411, 142, 710, 196, 732, 
        77, 439, 549, 17, 598, 254, 277, 37, 210, 155, 660, 656, 342, 540, 606, 437, 
        721, 754, 231, 596, 419, 439, 402, 595, 753, 756, 236, 277, 696, 432, 240, 155, 
        253, 171, 324, 82, 425, 753, 271, 18, 138, 161, 674, 632, 84, 662, 299, 35, 
        647, 683, 631, 296, 504, 415, 121, 640, 401, 509, 147, 479, 324, 539, 16, 577, 
        710, 493, 41, 518, 476, 465, 536, 767, 626, 592, 11, 92, 484, 463, 127, 361, 
        376, 140, 40, 110, 708, 161, 750, 491, 53, 280, 353, 529, 49, 369, 488, 760, 
        244, 530, 508, 103, 225, 426, 100, 233, 248, 111, 326, 732, 574, 605, 476, 332, 
        128, 516, 595, 218, 59, 575, 709, 264, 85, 292, 23, 287, 44, 512, 277, 440, 
        424, 167, 543, 649, 745, 25, 264, 223, 289, 742, 337, 245, 578, 195, 578, 88, 
        711, 403, 306, 153, 360, 245, 417, 598, 690, 593, 115, 116, 487, 544, 556, 141, 
        93, 482, 20, 68, 507, 436, 443, 26, 409, 162, 424, 217, 358, 384, 457, 451, 
        17, 763, 604, 529, 390, 404, 509, 310, 227, 624, 578, 96, 550, 365, 237, 25, 
        77, 257, 245, 736, 75, 688, 145, 484, 233, 569, 701, 591, 183, 388, 424, 352, 
        533, 411, 263, 154, 197, 3, 616, 576, 9, 425, 672, 560, 172, 139, 585, 401, 
        548, 61, 519, 623, 131, 664, 338, 364, 615, 421, 337, 180, 192, 144, 532, 725, 
        555, 26, 109, 752, 181, 108, 558, 190, 533, 460, 750, 87, 599, 718, 488, 377, 
        161, 237, 382, 292, 284, 720, 39, 129, 372, 528, 310, 564, 672, 224, 671, 457, 
        402, 11, 591, 583, 271, 379, 4, 186, 221, 136, 273, 202, 236, 143, 731, 397, 
        532, 344, 72, 46, 294, 263, 328, 48, 21, 20, 764, 694, 396, 666, 533, 29, 
        59, 507, 764, 330, 116, 150, 516, 490, 287, 171, 692, 523, 466, 654, 303, 228, 
        228, 375, 427, 674, 638, 137, 723, 41, 157, 717, 117, 553, 613, 651, 734, 54, 
        540, 729, 384, 38, 109, 282, 528, 548, 453, 451, 454, 149, 487, 757, 530, 97, 
    };

    cudaMemcpy(g_polyvec, h_polyvec, size, cudaMemcpyHostToDevice);

    // 计算NTT的执行时间
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);
    // for(int i = 0;i < 100000;i++){
    //     test_ntt<<<1, 64>>>(g_polyvec);
    // }
    test_ntt<<<500, 64>>>(g_polyvec);
    cudaDeviceSynchronize();

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("NTT执行时间为: %f ms\n", elapsedTime);

    cudaDeviceSynchronize();
    //将device中的值拷回到主机端
    cudaMemcpy(h_polyvec, g_polyvec, size, cudaMemcpyDeviceToHost);

    printf("NTT之后的数组为\n");
    for(int i = 0; i < CTRU_N; ++i){
        printf("%d ", h_polyvec[i]);
        if((i + 1) % 16 == 0) printf("\n");
    }

    int16_t inv_h_polyvec[CTRU_N] = {
        86 ,690 ,461 ,583 ,350 ,230 ,108 ,243 ,19 ,299 ,336 ,359 ,625 ,723 ,570 ,225 ,
        545 ,688 ,232 ,692 ,702 ,677 ,11 ,668 ,746 ,486 ,184 ,682 ,653 ,240 ,643 ,656 ,
        262 ,493 ,436 ,240 ,185 ,761 ,341 ,457 ,723 ,194 ,108 ,564 ,75 ,685 ,306 ,524 ,
        362 ,421 ,252 ,392 ,208 ,84 ,746 ,187 ,153 ,521 ,67 ,684 ,396 ,127 ,616 ,331 ,
        337 ,440 ,282 ,485 ,612 ,482 ,224 ,95 ,88 ,429 ,233 ,89 ,716 ,333 ,74 ,274 ,
        412 ,511 ,763 ,550 ,756 ,605 ,557 ,750 ,167 ,616 ,153 ,445 ,346 ,699 ,400 ,455 ,
        507 ,629 ,694 ,21 ,720 ,142 ,765 ,693 ,643 ,701 ,635 ,228 ,515 ,55 ,334 ,355 ,
        3 ,678 ,203 ,685 ,65 ,501 ,544 ,494 ,326 ,10 ,763 ,137 ,423 ,739 ,658 ,167 ,
        589 ,665 ,616 ,674 ,474 ,339 ,690 ,283 ,290 ,374 ,159 ,479 ,198 ,256 ,320 ,93 ,
        516 ,271 ,220 ,403 ,558 ,684 ,645 ,503 ,314 ,498 ,213 ,61 ,358 ,743 ,178 ,599 ,
        186 ,600 ,47 ,755 ,723 ,497 ,428 ,410 ,195 ,413 ,46 ,767 ,120 ,29 ,678 ,99 ,
        624 ,708 ,453 ,711 ,418 ,580 ,122 ,106 ,419 ,200 ,119 ,572 ,159 ,24 ,222 ,296 ,
        268 ,543 ,729 ,93 ,274 ,183 ,688 ,265 ,379 ,726 ,533 ,431 ,213 ,69 ,107 ,68 ,
        425 ,18 ,521 ,225 ,405 ,436 ,719 ,637 ,701 ,124 ,15 ,276 ,40 ,406 ,195 ,587 ,
        215 ,631 ,711 ,556 ,290 ,251 ,577 ,699 ,430 ,48 ,83 ,234 ,649 ,477 ,517 ,638 ,
        545 ,561 ,563 ,83 ,421 ,432 ,755 ,242 ,114 ,21 ,58 ,15 ,397 ,727 ,737 ,682 ,
        536 ,629 ,330 ,667 ,218 ,215 ,357 ,596 ,611 ,86 ,573 ,651 ,8 ,30 ,85 ,527 ,
        450 ,92 ,242 ,376 ,268 ,400 ,58 ,737 ,373 ,504 ,538 ,745 ,518 ,81 ,361 ,241 ,
        36 ,694 ,721 ,10 ,556 ,118 ,764 ,207 ,581 ,705 ,226 ,467 ,309 ,84 ,509 ,753 ,
        172 ,240 ,119 ,531 ,484 ,266 ,198 ,508 ,548 ,136 ,506 ,169 ,735 ,245 ,735 ,570 ,
        676 ,248 ,734 ,282 ,216 ,341 ,162 ,701 ,580 ,122 ,198 ,461 ,226 ,55 ,671 ,169 ,
        179 ,73 ,515 ,542 ,82 ,200 ,300 ,200 ,58 ,713 ,689 ,474 ,258 ,479 ,202 ,527 ,
        470 ,69 ,137 ,421 ,539 ,622 ,82 ,742 ,713 ,411 ,200 ,371 ,113 ,44 ,249 ,731 ,
        582 ,94 ,302 ,719 ,430 ,298 ,472 ,483 ,261 ,148 ,461 ,450 ,2 ,54 ,71 ,454 ,
        76 ,467 ,333 ,382 ,271 ,233 ,219 ,289 ,113 ,559 ,64 ,695 ,723 ,42 ,213 ,535 ,
        568 ,530 ,614 ,752 ,627 ,651 ,220 ,397 ,389 ,116 ,456 ,45 ,444 ,634 ,267 ,299 ,
        81 ,190 ,44 ,188 ,194 ,331 ,439 ,611 ,404 ,444 ,491 ,657 ,98 ,228 ,22 ,756 ,
        323 ,563 ,498 ,44 ,722 ,327 ,732 ,313 ,180 ,747 ,764 ,691 ,446 ,163 ,177 ,493 ,
        267 ,12 ,573 ,459 ,609 ,721 ,642 ,546 ,90 ,493 ,593 ,487 ,312 ,101 ,32 ,595 ,
        42 ,112 ,227 ,136 ,714 ,505 ,191 ,48 ,120 ,739 ,506 ,442 ,540 ,508 ,215 ,458 ,
        208 ,465 ,387 ,755 ,634 ,117 ,187 ,557 ,602 ,226 ,316 ,256 ,557 ,579 ,524 ,119 ,
        460 ,605 ,26 ,10 ,596 ,756 ,669 ,687 ,50 ,74 ,617 ,240 ,604 ,606 ,609 ,33 
    };

    cudaMemcpy(inv_g_polyvec, inv_h_polyvec, size, cudaMemcpyHostToDevice);
    
   // 计算逆NTT的执行时间
    cudaEventRecord(start, 0);

    // for(int i = 0;i < 100000;i++){
    //     test_inv_ntt<<<1, 64>>>(inv_g_polyvec);
    // }
    test_inv_ntt<<<1, 64>>>(inv_g_polyvec);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("逆NTT执行时间为: %f ms\n", elapsedTime);

    cudaDeviceSynchronize();

    //将device中的值拷回到主机端
    cudaMemcpy(inv_h_polyvec, inv_g_polyvec, size, cudaMemcpyDeviceToHost);

    // printf("逆NTT之后的数组为\n");
    // for(int i = 0; i < CTRU_N; ++i){
    //     printf("%d ", inv_h_polyvec[i]);
    //     if((i + 1) % 16 == 0) printf("\n");
    // }
    // // 清理资源
    // free(h_polyvec);
    // cudaFree(g_polyvec);
    return 0;
}
