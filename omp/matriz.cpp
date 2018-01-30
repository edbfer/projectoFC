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

	dirty = m1.dirty;
	det = m1.det;

	return *this;
}

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
