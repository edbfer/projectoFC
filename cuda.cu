#include "cuda.cuh"

__device__ float G;
__device__ float gama;
__device__ float omega;
__device__ float dt;
__device__ float h;

void cuda_setup(float G, float gama, float omega, float dt, float h)
{
  cudaMemcpyToSymbol("G",)
}
