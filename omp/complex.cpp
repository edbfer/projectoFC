#include "complex.h"
#include <cmath>

using namespace std;


__host__ __device__ complex::complex() : real(0.), im(0.)
{}

__host__ __device__ complex::complex(int real) : real(real), im(0.)
{}

__host__ __device__ complex::complex(float real): real(real), im(0.)
{}

__host__ __device__ complex::complex(double real): real(real), im(0.)
{}

__host__ __device__ complex::complex(float real, float im) : real(real), im(im)
{}

__host__ __device__ complex::complex(const complex& c): real(c.real), im(c.im)
{}

__host__ __device__ complex operator+(const complex& c1, const complex& c2)
{
	complex res(c1.real + c2.real, c1.im + c2.im);
	return res;
}

__host__ __device__ complex operator-(const complex& c1, const complex& c2)
{
	complex res(c1.real - c2.real, c1.im - c2.im);
	return res;
}

/*complex complex::operator-()
{
	complex res(-real, -im);
	return res;
}*/

__host__ __device__ complex complex::operator~() const
{
	complex res(real, -im);
	return res;
}

__host__ __device__ complex operator*(const complex& c1, const complex& c2)
{
	complex res(c1.real*c2.real - c1.im*c2.im, c1.real*c2.im + c2.real*c1.im);
	return res;
}

__host__ __device__ complex operator*(const float a, const complex& c1)
{
	complex res(c1.real * a, c1.im * a);
	return res;
}

__host__ __device__ complex operator*(const complex& c1, const float a)
{
	complex res(c1.real * a, c1.im * a);
	return res;
}

__host__ __device__ complex operator/(const complex& c1, const complex& c2)
{
	/*complex res(0, 0);
	float m = c2.real * c2.real + c2.im * c2.im;
	res.real = (c1.real * c2.real + c1.im * c2.im)/m;
	res.im	= (c1.im * c2.real - c1.real * c2.im)/m;
	return res;*/
	/*complex res = (c1*~c2)/(~c2*c2).real;
	return res;*/
	complex res, conj = ~c2;
	float norma;
	norma = c2.mod();
	norma = norma*norma;
	res.real = (conj*c1).real/norma;
	res.im = (conj*c1).im/norma;
	return res;

}

__host__ __device__ complex operator/(const complex& c1, const float a)
{
	complex res(c1.real / a, c1.im / a);
	return res;
}

__host__ __device__ complex operator/(const float a, const complex& c1)
{
	complex res(c1.real / a, c1.im / a);
	return res;
}

__host__ __device__ int operator==(const complex & c1, const complex & c2)
{
	if ((c1.real == c2.real) && (c1.im == c2.im))
	return 1;
	return 0;
}

ostream& operator<<(ostream& out, complex& c1)
{
	const char* sign = (c1.im > 0) ? "+" : ((c1.im < 0) ? "-" : "+");
	out << c1.real << sign << fabs(c1.im) << "i";
	//out << c1.real << " ";
	return out;
}

istream& operator>>(istream& in, complex& c1)
{
	in >> c1.real >> c1.im;
	return in;
}

__host__ __device__ float complex::mod() const
{
	return sqrt(real*real + im*im);
}
