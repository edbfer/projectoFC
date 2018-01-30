#ifndef _CUDA_CUH
#define _CUDA_CUH 1

#include <cuda.h>
#include <cuda_runtime.h>

#include "matriz.h"

void cuda_setup(float G, float gama, float omega, float dt, float h);

matriz cuda_doround(matriz& l);
float cuda_norm(matriz& l);
__global__ void cuda_psi1(complex* l, complex* r, int n, int m);
__global__ void cuda_psi2(complex* l, complex* psi1, complex* r, int n, int m);
__global__ void cuda_psin(complex* l, complex* psi2, complex* r, int n, int m);
__global__ void cuda_reduction(int n, int m, complex* l, float* out);
__device__ complex cuda_f(int i, int j, complex c, complex c1, complex c2, complex c3, complex c4);


#endif //_CUDA_CUH
