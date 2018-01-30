#include "cuda.cuh"

#include <iostream>
#include <cstdio>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <sstream>
#include <fstream>

__constant__ float cG;
__constant__ float cgama;
__constant__ float comega;
__constant__ float cdt;
__constant__ float ch;
extern float h;
extern float dt;

using namespace std;

__global__ void cuda_test()
{
  printf("device: G = %lf, gamma = %lf, omega = %lf, dt = %lf, h = %lf\n", cG, cgama, comega, cdt, ch);
}

void cuda_setup(float G, float gama, float omega, float dt, float h)
{
  cudaError_t err;
  err = cudaMemcpyToSymbol(cG, &G, sizeof(float));
  if(!(err == cudaSuccess)) {cout << "Failed CUDA Operation, in " << __FILE__ << " at line " << __LINE__ << ": " << cudaGetErrorString(err) << endl;}
  err = cudaMemcpyToSymbol(cgama, &gama, sizeof(float), 0, cudaMemcpyHostToDevice);
  if(!(err == cudaSuccess)) {cout << "Failed CUDA Operation, in " << __FILE__ << " at line " << __LINE__ << ": " << cudaGetErrorString(err) << endl;}
  err = cudaMemcpyToSymbol(comega, &omega, sizeof(float), 0, cudaMemcpyHostToDevice);
  if(!(err == cudaSuccess)) {cout << "Failed CUDA Operation, in " << __FILE__ << " at line " << __LINE__ << ": " << cudaGetErrorString(err) << endl;}
  err = cudaMemcpyToSymbol(cdt, &dt, sizeof(float), 0, cudaMemcpyHostToDevice);
  if(!(err == cudaSuccess)) {cout << "Failed CUDA Operation, in " << __FILE__ << " at line " << __LINE__ << ": " << cudaGetErrorString(err) << endl;}
  err = cudaMemcpyToSymbol(ch, &h, sizeof(float), 0, cudaMemcpyHostToDevice);
  if(!(err == cudaSuccess)) {cout << "Failed CUDA Operation, in " << __FILE__ << " at line " << __LINE__ << ": " << cudaGetErrorString(err) << endl;}

  cudaThreadSynchronize();

  cuda_test<<<1, 1>>>();
}

float cuda_norm(complex* l, float* linhas, int n, int m)
{
  dim3 tpb(8);
  dim3 nb(n / tpb.x);
  cuda_reduction<<<nb, tpb>>>(n, m, l , linhas);
  //cuda_reduction<<<nb, tpb>>>(n/8, 8, linhas, linhas);
  //cuda_reduction<<<nb, tpb>>>(1, n/8, linhas, linhas);
  float* host = new float[n];
  cudaMemcpy(host, linhas, sizeof(float)*n, cudaMemcpyDeviceToHost);
  float v = 0;
  for(int i = 0; i<n; i++)
  {
    v = v + host[i];
  }
  return v;
}

float cuda_norm(matriz& l)
{
  complex* val;
  float* linhas;
  cudaMalloc(&val, sizeof(complex)*l.n*l.m);
  cudaMalloc(&linhas, sizeof(float)*l.n);
  cudaThreadSynchronize();

  cudaMemcpy(val, l.mat, sizeof(complex)*l.n*l.m, cudaMemcpyHostToDevice);

  dim3 tpb(8);
  dim3 nb(l.n / tpb.x);
  cuda_reduction<<<nb, tpb>>>(l.n, l.m, val, linhas);

  float* host = new float[l.n];
  cudaMemcpy(host, linhas, sizeof(float)*l.n, cudaMemcpyDeviceToHost);

  float v = 0;
  for(int i = 0; i<l.n; i++)
  {
    v = v + host[i];
  }

  cudaFree(val);
  cudaFree(linhas);
  delete[] host;

  return v*h*h;

}

__global__ void cuda_reduction(int n, int m, complex* val, float* linhas)
{
  int x = (blockIdx.x * blockDim.x) + threadIdx.x;

  if((x > n) || (x < 0))
    return;

  float v = 0;
  for(int i = 0; i < m; i++)
  {
    float mod = val[x * m + i].mod();
    v = v + mod*mod*ch*ch;
  }

  linhas[x] = v;
}

matriz cuda_doround(matriz& l)
{
  complex* psi1, *lmat, *psi2, *r;

  int n = l.n, m = l.m;
  cudaMalloc(&lmat, sizeof(complex)*n*m);
  cudaMalloc(&psi1, sizeof(complex)*n*m);
  cudaMalloc(&psi2, sizeof(complex)*n*m);
  cudaMalloc(&r, sizeof(complex)*n*m);
  cudaThreadSynchronize();

  cudaMemcpy(lmat, l.mat, sizeof(complex)*n*m, cudaMemcpyHostToDevice);

  dim3 tpb(8, 4);
  dim3 nb(n / tpb.x, m / tpb.y);
  cuda_psi1<<<nb, tpb>>>(lmat, psi1, n, m);
  cudaThreadSynchronize();
  cuda_psi2<<<nb, tpb>>>(lmat, psi1, psi2, n, m);
  cudaThreadSynchronize();
  cuda_psin<<<nb, tpb>>>(lmat, psi2, r, n, m);
  cudaThreadSynchronize();

  matriz res(128, 128);

  cudaMemcpy(res.mat, r, sizeof(complex)*n*m, cudaMemcpyDeviceToHost);

  float norm = cuda_norm(res);
  res = res * complex(1/norm, 0.0f);

  cout << "Norma:" << norm << endl;

  cudaFree(lmat);
  cudaFree(psi1);
  cudaFree(psi2);
  cudaFree(r);
  return res;

}

__global__ void cuda_divideandmod(int n, int m, float val, complex* l, complex* r)
{
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;
  int j = (blockIdx.y * blockDim.y) + threadIdx.y;

  complex c = l[i * m + j] / val;
  float modd = c.mod();
  r[i * m + j] = modd*modd;
}

void cuda_simulate(matriz& l, float tmax)
{
    matriz res(l.n, l.m);
    complex *mat, *psi1, *psi2, *psi3, *toHost;
    float* linhas;

    cudaMalloc(&linhas, sizeof(float)*l.n);
    cudaMalloc(&mat, sizeof(complex)*l.n*l.m);
    cudaMalloc(&psi1, sizeof(complex)*l.n*l.m);
    cudaMalloc(&psi2, sizeof(complex)*l.n*l.m);
    cudaMalloc(&psi3, sizeof(complex)*l.n*l.m);
    cudaMalloc(&toHost, sizeof(complex)*l.n*l.m);
    cudaThreadSynchronize();

    cudaMemcpy(mat, l.mat, sizeof(complex)*l.n*l.m, cudaMemcpyHostToDevice);

    dim3 tpb(8, 4);
    dim3 nb(l.n / tpb.x, l.m / tpb.y);

    int iter = 0;
    for(float t = 0.0f; t<tmax; t += dt)
    {
      cuda_psi1<<<nb, tpb>>>(mat, psi1, l.n, l.m);
      cuda_psi2<<<nb, tpb>>>(mat, psi1, psi2, l.n, l.m);
      cuda_psin<<<nb, tpb>>>(mat, psi2, psi3, l.n, l.m);

      float f = cuda_norm(psi3, linhas, l.n, l.m);
      f = sqrt(f);
      cuda_divideandmod<<<nb, tpb>>>(l.n, l.m, f, psi3, mat);
      if(iter % 10 == 0)
      {
        cudaMemcpy(toHost, mat, sizeof(complex)*l.n*l.m, cudaMemcpyDeviceToDevice);
        pid_t pid = fork();
        if(pid == 0)
        {
          cudaMemcpyAsync(res.mat, toHost, sizeof(complex)*l.n*l.m, cudaMemcpyDeviceToHost);
          stringstream s;
          s << "dados/t" << iter << ".txt";
          ofstream o(s.str().c_str());
          res.printCoord(o);
          o.close();
          kill(getpid(), SIGTERM);
        }
      }
      iter++;
    }

    cudaFree(linhas); cudaFree(psi1); cudaFree(psi2); cudaFree(psi3); cudaFree(mat);
    return;
}

__device__ complex cuda_f(int i, int j, complex c, complex c1, complex c2, complex c3, complex c4)
{
  float x = -10.0f + (i-1)*ch;
  float y = -10.0f + (j-1)*ch;

  /*complex lapx = (c2 - 2.0f*c + c1)/(ch*ch);
  complex lapy = (c4 - 2.0f*c + c3)/(ch*ch);*/
  complex lap = ((c2 - 2.0f*c + c1)/(ch*ch) + (c4 - 2.0f*c + c3)/(ch*ch))*(-0.5f);

  complex p2 = c * ((x*x) + (y*y))/(2.0f);
  float n = c.mod();
  n *= n;
  complex p3 = c * (cG*n);

  complex dx = ((c2 - c1)/(2.0f*ch))*y;
  complex dy = ((c4 - c3)/(2.0f*ch))*x;
  /*dx = dx * y;
  dy = dy * x;*/
  complex p4 = (dy - dx) * complex(0, comega);

  complex r = (lap + p2 + p3 - p4)/(complex(-cgama, 1.0f));
  return r
  ;
}

__global__ void cuda_psi1(complex* l, complex* r, int n, int m)
{
  //RungeKuttadevMat
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;
  int j = (blockIdx.y * blockDim.y) + threadIdx.y;

  if((i <= 0) || (i >= 127) || (j <= 0) || (j >= 127))
    return;

  complex c = l[i * m + j];
  complex c1 = l[(i-1) * m + j];
  complex c2 = l[(i+1) * m + j];
  complex c3 = l[i * m + (j-1)];
  complex c4 = l[i * m + (j+1)];

  complex ft = cuda_f(i, j, c, c1, c2, c3, c4);
  r[i * m + j] = c + (ft*cdt);
}

__global__ void cuda_psi2(complex* l, complex* psi1, complex* r, int n, int m)
{
  //RungeKutta
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;
  int j = (blockIdx.y * blockDim.y) + threadIdx.y;

  if((i <= 0) || (i >= 127) || (j <= 0) || (j >= 127))
    return;

  complex c = psi1[i * m + j];
  complex la = l[i * m + j];
  complex c1 = psi1[(i-1) * m + j];
  complex c2 = psi1[(i+1) * m + j];
  complex c3 = psi1[i * m + (j-1)];
  complex c4 = psi1[i * m + (j+1)];

  complex ft = cuda_f(i, j, c, c1, c2, c3, c4);
  ft = ft*0.25f*cdt;
  complex ant = c*0.25f;
  complex lat = la * 0.75f;

  r[i * m + j] = ft + ant + lat;
}

__global__ void cuda_psin(complex* l, complex* psi2, complex* r, int n, int m)
{
  //RungeKutta
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;
  int j = (blockIdx.y * blockDim.y) + threadIdx.y;

  if((i <= 0) || (i >= 127) || (j <= 0) || (j >= 127))
    return;

  complex c = psi2[i * m + j];
  complex la = l[i * m + j];
  complex c1 = psi2[(i-1) * m + j];
  complex c2 = psi2[(i+1) * m + j];
  complex c3 = psi2[i * m + (j-1)];
  complex c4 = psi2[i * m + (j+1)];

  complex ft = cuda_f(i, j, c, c1, c2, c3, c4);
  ft = ft * (cdt * (2.0f/3.0f));
  complex ant = c * (2.0f/3.0f);
  complex lat = la * (1.0f/3.0f);

  r[i * m + j] = ft + ant + lat;
}
