#ifndef _MATRIZ_H
#define _MATRIZ_H 1

#include <cmath>
#include <iostream>
#include "complex.h"
#include <vector>

using namespace std;

struct swapOut
{
	complex pivot;
	int swapped;
	int l1;
	int l2;
};

struct gaussInfo
{
	int ispermutation;
	int p1;
	int p2;
	complex v;
};

class matriz
{
public:
	//complex ** mat;
	complex * mat;
	int n, m;
	mutable complex det;
	mutable int dirty;
	mutable int sigma;

	matriz(const int n = 3, const int m = 3);
	matriz(const matriz& m1);
	~matriz();

	void fill(istream& in);
	void fill(complex c);
	matriz transform(complex (*f) (int i, int j, complex f));

	complex& operator()(const int& n, const int& m) const;
	complex& operator()(const int& n, const int& m);

	complex* operator[](const size_t i) const;
	matriz& operator=(const matriz& m1);

	matriz operator-();
	matriz operator++(int);
	matriz operator~();

	friend matriz transpose(matriz& m1);
	static matriz inverte(matriz& m1);
	friend matriz extend(matriz& m1, matriz& m2);
	friend matriz extract(matriz& m1, int x0, int y0, int x1, int y1);

	friend matriz operator^(const matriz& m1, const int e);
	friend matriz operator+(const matriz& m1, const matriz& m2);
	//friend matriz sum_cuda(matriz& m1, matriz& m2);
	friend matriz operator-(const matriz& m1, const matriz& m2);
	friend matriz operator*(const matriz& m1, const matriz& m2);
	//matriz operator*=(complex a);
	friend matriz operator~(const matriz& m1);
	//friend matriz multiply_cuda(matriz& m1, matriz& m2);
	friend matriz operator*(const matriz& m1, complex v);

	friend complex dot(matriz& m1, matriz& m2);

	friend ostream& operator<<(ostream& out, matriz& m1);
	friend istream& operator>> (istream& in, matriz& m1);

	void printCoord(ostream& o);
};

#endif //_MATRIX_H
