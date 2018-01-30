#ifndef _CUDA_CUH
#define _CUDA_CUH 1

#include <cuda.h>
#include <cuda_runtime.h>

void cuda_setup(float G, float gama, float omega, float dt, float h);

__global__ void cuda_iterate(complex* l, complex* r);
__device__ complex cuda_f(complex c);

#endif //_CUDA_CUH
