#pragma once

#include <iostream>

using namespace std;

class complex
{

public:
	float real;
	float im;

	complex();
	complex(int real);
	complex(float real);
	complex(double real);
	explicit complex(float real, float im);
	complex(const complex& c);

	friend complex operator+(const complex& c1, const complex& c2);
	friend complex operator-(const complex& c1, const complex& c2);
	complex operator-() const;
	complex operator~() const;
	friend complex operator*(const complex& c1, const complex& c2);
	friend complex operator*(const float a, const complex& c1);
	friend complex operator*(const complex& c1, const float a);
	friend complex operator/(const complex& c1, const complex& c2);
	friend complex operator/(const complex& c1, const float a);
	friend complex operator/(const float a, const complex& c1);

	friend int operator==(const complex& c1, const complex& c2);

	operator float()
	{
		return real;
	}

	friend ostream& operator<<(ostream& out, complex& c1);
	friend istream& operator>>(istream& in, complex& c1);

	float mod() const;
};
