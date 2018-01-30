#include <cstdlib>
#include <cmath>
#include <iostream>
#include "matriz.h"
#include "complex.h"
#include <vector>
#include <cstring>
#include <omp.h>

using namespace std;

matriz::matriz(const int n, const int m): n(n), m(m)
{
	mat = new complex[n*m];
	dirty = 1;
	det = 0.;
}

matriz::matriz(const matriz& m1): n(m1.n), m(m1.m)
{
	mat = new complex[n*m];
	memcpy(mat, m1.mat, sizeof(complex)*n*m);
	det = m1.det;
	dirty = m1.dirty;
}

matriz::~matriz()
{
	delete[] mat;
}

matriz matriz::id(int n)
{
	matriz res(n, n);
	for(int i = 0; i<n; i++)
	{
		res(i, i) = ((complex) 1.);
	}
	return res;
}

matriz permutation(int n, int l1, int l2)
{
	matriz res = matriz::id(n);
	res.swapLine(l1, l2);
	return res;
}

matriz matriz::random(int n, int m)
{
	matriz res(n, m);
	for(int i = 0; i<n; i++)
	{
		for(int j = 0; j<m; j++)
		{
			res(i, j) = (complex) ((double)rand()/(double) RAND_MAX)*10;
		}
	}
	return res;
}

matriz matriz::hermitian(int n)
{
	matriz res(n, n);
	for(int i = 0; i<n; i++){
		for(int j = 0; j<i+1; j++)
		{
			res(i, j) = (complex) ((double)rand()/(double) RAND_MAX)*10;
			res(j, i) = res(i, j);
		}
	}
	return res;
}

matriz matriz::tridiagonal(complex d0, complex dp, complex d1, int x)
{
	/*matriz res(x, x);
	for (int i = 1; i<x-1-1-1; i++)
	{
		res(i, i) = dp;
		int j2 = i + 1;
		res(i, j2) = d1;
		res(j2, i) = d0;
	}
	int j2 = x-1-1;
	res(j2, j2) = dp;
	res()
	return res;*/
}

complex& matriz::operator()(const int& x, const int& y)
{
	return mat[x * m + y];
}

complex& matriz::operator()(const int& x, const int& y) const
{
	return mat[x * m + y];
}

complex* matriz::operator[](const size_t i) const
{
	return (mat + i*m);
}

matriz& matriz::operator=(const matriz& m1)
{
	delete[] mat;

	n = m1.n;
	m = m1.m;

	mat = new complex[n*m];
	for(int i = 0; i<n; i++)
	{
		for(int j = 0; j<m; j++)
		{
			(*this)(i, j) = m1(i, j);
		}
	}
	//memcpy(mat, m1.mat, sizeof(complex)*n*m);

	dirty = m1.dirty;
	det = m1.det;

	return *this;
}

/*matriz matriz::operator-()
{
	for(int i = 0; i<n; i++)
	{
		for(int j = 0; j<m; j++)
		{
			(*this)(i, j) = -(*this)(i, j);
		}
	}
	return *this;
}*/

matriz matriz::operator++(int)
{
	for(int i = 0; i<n; i++)
	{
		for(int j = 0; j<m; j++)
		{
			(*this)(i, j) = (*this)(i, j) + (complex) 1.;
		}
	}
	return *this;
}

matriz transpose(matriz& m1)
{
	matriz r(m1.m, m1.n);
	for(int i = 0; i<m1.n; i++)
	{
		for(int j = 0; j<m1.m; j++)
		{
			r(j, i) = m1(i, j);
		}
	}
	return r;
}


matriz extend(matriz& m1, matriz& m2)
{
	matriz t(m1.n, m1.m+m2.m);
	for(int i = 0; i<t.n; i++)
	{
		for(int j = 0; j<m1.m; j++)
		{
			t(i, j) = m1(i, j);
		}

		for(int j = 0; j<m2.m; j++)
		{
			t(i, j+m1.m) = m2(i, j);
		}
	}

	return t;
}

matriz extract(matriz& m1, int x0, int y0, int x1, int y1)
{
	matriz r(x1 - x0 + 1, y1 - y0 + 1);
	for (int i = x0; i <=x1; i++)
	{
		for (int j = y0; j <= y1; j++)
		{
			r(i - x0, j - y0) = m1(i, j);
		}
	}
	return r;
}

void matriz::swapLine(int l1, int l2)
{
	complex* r = new complex[m];
	memcpy(r, (*this)[l1], sizeof(complex)*m);
	memcpy((*this)[l1], (*this)[l2], sizeof(complex)*m);
	memcpy((*this)[l2], r, sizeof(complex)*m);
}

matriz matriz::swapColForVector(int col, matriz& vec)
{
	for (int i = 0; i < n; i++)
	{
		(*this)(i, col) = vec(i, 0);
	}
	return *this;
}

matriz operator^(const matriz& m1, const int e)
{
	matriz r = m1;
	for(int i = 0; i<e; i++)
	{
		r = r * m1;
	}
	return r;
}

matriz operator+(const matriz& m1, const matriz& m2)
{
	matriz res(m1.n, m1.m);
	for(int i = 0; i<m1.n; i++)
	{
		for(int j = 0; j<m1.m; j++)
		{
			res(i, j) = m1(i, j) + m2(i, j);
		}
	}
	return res;
}
// matriz sum_cuda(matriz& m1, matriz& m2);
matriz operator-(const matriz& m1, const matriz& m2)
{
	matriz res(m1.n, m1.m);
	for(int i = 0; i<m1.n; i++)
	{
		for(int j = 0; j<m1.m; j++)
		{
			res(i, j) = m1(i, j) - m2(i, j);
		}
	}
	return res;
}

matriz operator*(const matriz& m1, const matriz& m2)
{
	matriz newm(m1.n, m2.m);
	for (int i = 0; i<newm.m; i++)
	{
		for(int j = 0; j<newm.n; j++)
		{
			complex v = 0.;
			for(int k = 0; k<m1.m; k++)
			{
				complex v1 = m2(k, i);
				complex v2 = m1(j, k);
				v = v + v1 * v2;
			}
			newm(j, i) = v;
		}
	}

	return newm;
}

matriz matriz::operator~()
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			~(*this)(i, j);
		}
	}
	return *this;
}
// matriz multiply_cuda(matriz& m1, matriz& m2);
matriz operator*(const matriz& m1, complex v)
{
	matriz res = m1;
	int i = 0, j = 0;
	for(; i<m1.n; i++)
	{
		for(j = 0; j<m1.m; j++)
		{
			res(i, j) = m1(i, j) * v;
		}
	}

	return res;
}

complex dot(matriz& m1, matriz& m2)
{
	matriz temp = transpose(m1)*(~m2);
	return temp(0, 0);
}

ostream& operator<<(ostream& out, matriz& m1)
{
	for (int i = 0; i<m1.n; i++)
	{
		for(int j = 0; j<m1.m; j++)
		{
			out << m1(i, j) << " ";
		}

		out << endl;
	}
	return out;
}

istream& operator>> (istream& in, matriz& m1)
{
	for (int i = 0; i<m1.n; i++)
	{
		for(int j = 0; j<m1.m; j++)
		{
			in >> m1(i, j);
		}
	}
	return in;
}

void matriz::fill(istream& in)
{
	while(!in.eof())
	{
		int x, y;
		in >> x >> y;
		in >> (*this)(x, y);
	}
}

void matriz::fill(complex val)
{
	#pragma omp parallel for
	for(int i = 0; i<n; i++)
	{
		#pragma omp simd
		for(int j = 0; j<m; j++)
		{
			(*this)(i, j) = val;
		}
	}
}

void matriz::printCoord(ostream& o)
{
	for(int i = 0; i<n; i++)
	{
		for(int j = 0; j<m; j++)
		{
			o << i << "\t" << j << "\t" << (*this)(i, j) << endl;
		}
	}
}

matriz matriz::transform(complex (*func) (int i, int j, complex vij))
{
	matriz res(n, m);

	for(int i = 0; i<n; i++)
	{
		for(int j = 0; j<m; j++)
		{
			res(i, j) = func(i, j, (*this)(i, j));
		}
	}
	return res;
}
