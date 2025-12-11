// #include <cuda_runtime.h>
// #include <stdint.h>
// #include <curand_kernel.h>
// #include "params.h"
// #include "../common/randombytes.h"
// #include "reduce.h"
// #include <stdio.h>

// // 定义蒙哥马利约简内核
// __global__ void montgomery_reduce_kernel(int32_t* a, int16_t* result, size_t n) {
//     unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     if (idx < n) {
//         int32_t t;
//         int16_t u;

//         u = a[idx] * QINV;
//         t = (int32_t)u * CTRU_Q;
//         t = a[idx] - t;
//         t >>= 16;
//         result[idx] = t;
//     }
// }

// // 定义巴雷特约简内核
// __global__ void barrett_reduce_kernel(int16_t* a, int16_t* result, size_t n) {
//     unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     if (idx < n) {
//         int32_t t;

//         t = BARRETT_V * (int32_t)a[idx];
//         t >>= 24;
//         t *= CTRU_Q;
//         result[idx] = a[idx] - t;
//     }
// }

// // 定义FQC减法内核
// __global__ void fqcsubq_kernel(int16_t* a, int16_t* result, size_t n) {
//     unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     if (idx < n) {
//         int16_t t = a[idx];

//         t += (t >> 15) & CTRU_Q;
//         t -= CTRU_Q;
//         t += (t >> 15) & CTRU_Q;
//         result[idx] = t;
//     }
// }

// // 定义FQ逆元内核
// __global__ void fqinv_kernel(int16_t* a, int16_t* result, size_t n) {
//     unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     if (idx < n) {
//         unsigned int i;
//         int16_t t = 1;
//         for (i = 1; i <= CTRU_Q - 2; ++i) {
//             t = (int32_t)t * (int32_t)a[idx] % CTRU_Q;
//         }
//         result[idx] = t;
//     }
// }
// // 初始化CURAND状态的内核
// __global__ void setup_curand_states(curandState *states, unsigned long seed) {
//     int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     if (idx < gridDim.x * blockDim.x) {
//         curand_init(seed, idx, 0, &states[idx]);
//     }
// }


// // 定义生成均匀随机数的内核
// __global__ void fquniform_kernel(curandState *state, int16_t* result, size_t n) {
//     unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     if (idx < n) {
//         int16_t r;
//         do {
//             r = curand(state + idx) & 0x3FF;
//         } while (r >= CTRU_Q);
//         result[idx] = r;
//     }
// }

// // 主机端代码
// void run_montgomery_reduce(int32_t* host_a, int16_t* host_result, size_t n) {
//     int32_t *dev_a;
//     int16_t *dev_result;
//     size_t bytes = n * sizeof(int32_t);

//     // 分配设备内存
//     cudaMalloc((void**)&dev_a, bytes);
//     cudaMalloc((void**)&dev_result, n * sizeof(int16_t));

//     // 将数据从主机复制到设备
//     cudaMemcpy(dev_a, host_a, bytes, cudaMemcpyHostToDevice);

//     // 定义线程块大小
//     int threadsPerBlock = 1024;
//     int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;

//     // 启动内核
//     montgomery_reduce_kernel<<<blocksPerGrid, threadsPerBlock>>>(dev_a, dev_result, n);

//     // 检查错误
//     cudaError_t err = cudaGetLastError();
//     if (err != cudaSuccess)
//         printf("Error: %s\n", cudaGetErrorString(err));

//     // 将结果从设备复制回主机
//     cudaMemcpy(host_result, dev_result, n * sizeof(int16_t), cudaMemcpyDeviceToHost);

//     // 释放设备内存
//     cudaFree(dev_a);
//     cudaFree(dev_result);
// }

// void run_barrett_reduce(int16_t* host_a, int16_t* host_result, size_t n) {
//     int16_t *dev_a, *dev_result;
//     size_t bytes = n * sizeof(int16_t);

//     // 分配设备内存
//     cudaMalloc((void**)&dev_a, bytes);
//     cudaMalloc((void**)&dev_result, bytes);

//     // 将数据从主机复制到设备
//     cudaMemcpy(dev_a, host_a, bytes, cudaMemcpyHostToDevice);

//     // 定义线程块大小
//     int threadsPerBlock = 256;
//     int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;

//     // 启动内核
//     barrett_reduce_kernel<<<blocksPerGrid, threadsPerBlock>>>(dev_a, dev_result, n);

//     // 检查错误
//     cudaError_t err = cudaGetLastError();
//     if (err != cudaSuccess)
//         printf("Error: %s\n", cudaGetErrorString(err));

//     // 将结果从设备复制回主机
//     cudaMemcpy(host_result, dev_result, bytes, cudaMemcpyDeviceToHost);

//     // 释放设备内存
//     cudaFree(dev_a);
//     cudaFree(dev_result);
// }

// void run_fqcsubq(int16_t* host_a, int16_t* host_result, size_t n) {
//     int16_t *dev_a, *dev_result;
//     size_t bytes = n * sizeof(int16_t);

//     // 分配设备内存
//     cudaMalloc((void**)&dev_a, bytes);
//     cudaMalloc((void**)&dev_result, bytes);

//     // 将数据从主机复制到设备
//     cudaMemcpy(dev_a, host_a, bytes, cudaMemcpyHostToDevice);

//     // 定义线程块大小
//     int threadsPerBlock = 256;
//     int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;

//     // 启动内核
//     fqcsubq_kernel<<<blocksPerGrid, threadsPerBlock>>>(dev_a, dev_result, n);

//     // 检查错误
//     cudaError_t err = cudaGetLastError();
//     if (err != cudaSuccess)
//         printf("Error: %s\n", cudaGetErrorString(err));

//     // 将结果从设备复制回主机
//     cudaMemcpy(host_result, dev_result, bytes, cudaMemcpyDeviceToHost);

//     // 释放设备内存
//     cudaFree(dev_a);
//     cudaFree(dev_result);
// }

// void run_fqinv(int16_t* host_a, int16_t* host_result, size_t n) {
//     int16_t *dev_a, *dev_result;
//     size_t bytes = n * sizeof(int16_t);

//     // 分配设备内存
//     cudaMalloc((void**)&dev_a, bytes);
//     cudaMalloc((void**)&dev_result, bytes);

//     // 将数据从主机复制到设备
//     cudaMemcpy(dev_a, host_a, bytes, cudaMemcpyHostToDevice);

//     // 定义线程块大小
//     int threadsPerBlock = 256;
//     int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;

//     // 启动内核
//     fqinv_kernel<<<blocksPerGrid, threadsPerBlock>>>(dev_a, dev_result, n);

//     // 检查错误
//     cudaError_t err = cudaGetLastError();
//     if (err != cudaSuccess)
//         printf("Error: %s\n", cudaGetErrorString(err));

//     // 将结果从设备复制回主机
//     cudaMemcpy(host_result, dev_result, bytes, cudaMemcpyDeviceToHost);

//     // 释放设备内存
//     cudaFree(dev_a);
//     cudaFree(dev_result);
// }

// void run_fquniform(int16_t* host_result, size_t n) {
//     int16_t *dev_result;
//     size_t bytes = n * sizeof(int16_t);
//     curandState *dev_states;

//     // 分配设备内存
//     cudaMalloc((void**)&dev_result, bytes);
//     cudaMalloc((void**)&dev_states, n * sizeof(curandState));

//     // 初始化随机数状态
//     setup_curand_states<<<(n + 255) / 256, 256>>>(dev_states, time(NULL));

//     // 定义线程块大小
//     int threadsPerBlock = 256;
//     int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;

//     // 启动内核
//     fquniform_kernel<<<blocksPerGrid, threadsPerBlock>>>(dev_states, dev_result, n);

//     // 检查错误
//     cudaError_t err = cudaGetLastError();
//     if (err != cudaSuccess)
//         printf("Error: %s\n", cudaGetErrorString(err));

//     // 将结果从设备复制回主机
//     cudaMemcpy(host_result, dev_result, bytes, cudaMemcpyDeviceToHost);

//     // 释放设备内存
//     cudaFree(dev_result);
//     cudaFree(dev_states);
// }



// int main() {
//     const size_t n = 10000000; // 处理的数据数量
//     int32_t *host_a = (int32_t*)malloc(n * sizeof(int32_t));
//     int16_t *host_result = (int16_t*)malloc(n * sizeof(int16_t));
//     // 初始化输入数据
//     for (size_t i = 0; i < n; ++i) {
//         host_a[i] = (int32_t)(rand() % 10000);
//     }
//     // for (size_t i = 0; i < n; ++i) {
//     //     printf("%d ,",host_a[i]);
//     //     if(i % 16 == 0){
//     //         printf("\n");
//     //     }
//     // }
//     cudaEvent_t start, stop;
//     float milliseconds = 0;

//     // 初始化CUDA事件
//     cudaEventCreate(&start);
//     cudaEventCreate(&stop);


//     // 运行蒙哥马利约简
//     cudaEventRecord(start, 0);

//     run_montgomery_reduce(host_a, host_result, n);

//     cudaEventRecord(stop, 0);
//     cudaEventSynchronize(stop);
//     cudaEventElapsedTime(&milliseconds, start, stop);
//     printf("Montgomery reduce time: %f ms\n", milliseconds);

//     printf("Montgomery reduce results:\n");
//     for (size_t i = 0; i < 10; ++i) { // 打印前10个结果
//         printf("%d ", host_result[i]);
//     }
//     printf("\n");

//     // // 运行巴雷特约简
//     // run_barrett_reduce(host_result, host_result, n);
//     // printf("Barrett reduce results:\n");
//     // for (size_t i = 0; i < 10; ++i) { // 打印前10个结果
//     //     printf("%d ", host_result[i]);
//     // }
//     // printf("\n");

//     // // 运行FQC减法
//     // run_fqcsubq(host_result, host_result, n);
//     // printf("FQC sub Q results:\n");
//     // for (size_t i = 0; i < 10; ++i) { // 打印前10个结果
//     //     printf("%d ", host_result[i]);
//     // }
//     // printf("\n");

//     // // 运行FQ逆元
//     // run_fqinv(host_result, host_result, n);
//     // printf("FQ inverse results:\n");
//     // for (size_t i = 0; i < 10; ++i) { // 打印前10个结果
//     //     printf("%d ", host_result[i]);
//     // }
//     // printf("\n");

//     // // 生成均匀随机数
//     // run_fquniform(host_result, n);
//     // printf("Uniform random numbers:\n");
//     // for (size_t i = 0; i < 10; ++i) { // 打印前10个结果
//     //     printf("%d ", host_result[i]);
//     // }
//     // printf("\n");

//     // 释放CUDA事件
//     cudaEventDestroy(start);
//     cudaEventDestroy(stop);

//     // 释放主机内存
//     free(host_a);
//     free(host_result);

//     return 0;
// }
