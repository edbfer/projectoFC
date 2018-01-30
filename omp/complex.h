#pragma once

#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>

using namespace std;

class complex
{

public:
	float real;
	float im;

	__host__ __device__ complex();
	__host__ __device__ complex(int real);
	__host__ __device__ complex(float real);
	__host__ __device__ complex(double real);
	__host__ __device__ explicit complex(float real, float im);
	__host__ __device__ complex(const complex& c);

	__host__ __device__ friend complex operator+(const complex& c1, const complex& c2);
	__host__ __device__ friend complex operator-(const complex& c1, const complex& c2);
	__host__ __device__ complex operator-() const;
	__host__ __device__ complex operator~() const;
	__host__ __device__ friend complex operator*(const complex& c1, const complex& c2);
	__host__ __device__ friend complex operator*(const float a, const complex& c1);
	__host__ __device__ friend complex operator*(const complex& c1, const float a);
	__host__ __device__ friend complex operator/(const complex& c1, const complex& c2);
	__host__ __device__ friend complex operator/(const complex& c1, const float a);
	__host__ __device__ friend complex operator/(const float a, const complex& c1);

	__host__ __device__ friend int operator==(const complex& c1, const complex& c2);

	__host__ __device__ operator float()
	{
		return real;
	}

	friend ostream& operator<<(ostream& out, complex& c1);
	friend istream& operator>>(istream& in, complex& c1);

	__host__ __device__ float mod() const;
};
