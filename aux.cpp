#include "aux.h"

extern float h;
extern float G;

complex aux::x2y2(int i, int j, complex c)
{
  float x = -10.0f + i*h;
  float y = -10.0f + j*h;
  complex m = ((x*x) + (y*y))/2.0f;
  //cout << "("<< x << "," << y << "," << "," << c << "," << c*m << ")" << endl;
  return c*m;

}

complex aux::my(int i, int j, complex c)
{
  float x = -10.0f + i*h;
  return c * x;
}


complex aux::gnorm(int i, int j, complex c)
{
  float n = c.mod();
  return c * (G*n*n);
}

complex aux::mx(int i, int j, complex c)
{
  float x = -10.0f + i*h;
  return c * x;
}

complex aux::conjugate(int i, int j, complex c)
{
  complex a =  ~c;
  //cout << "~c: " << a << " " << "c: " << c << endl;
  return c * (~c);
}

complex aux::norma(matriz& psi)
{
  matriz r = psi.transform(aux::conjugate);
  complex v = 0;
  for(int i = 1; i<psi.n-1; i++)
  {
    for(int j = 1; j<psi.m-1; j++)
    {
      v = v + r(i, j)*h*h;
    }
  }
  return v;
}

matriz aux::dx(matriz& l)
{
  matriz res(l.n, l.m);
  /*matriz d = matriz::tridiagonal(-1, 0, 1, l.n);
  res = d * l;*/
  for(int i = 1; i<l.n-1; i++)
  {
    for(int j = 1; j<l.m-1; j++)
    {
      res(i, j) = l(i + 1, j) - l(i - 1, j);
      res(i, j) = res(i, j) / (2.0f*h);
    }
  }
  return res;
}

matriz aux::dy(matriz& l)
{
  matriz res(l.n, l.m);
  /*matriz d = matriz::tridiagonal(-1, 0, 1, l.n);
  matriz t = transpose(l);
  res = d * t;*/
  for(int i = 1; i<l.n-1; i++)
  {
    for(int j = 1; j<l.m-1; j++)
    {
      res(i, j) = l(i, j + 1) + l(i, j - 1);
      res(i, j) = res(i, j) / (2.0f*h);
    }
  }
  return res;
}

matriz aux::d2x(matriz& l)
{
  matriz res(l.n, l.m);
  /*matriz d = matriz::tridiagonal(1, -2, 1, l.n);
  cout << d << endl;
  res = d * l;*/
  for(int i = 1; i<l.n-1; i++)
  {
    for(int j = 1; j<l.m-1; j++)
    {
      res(i, j) = l(i + 1, j) - l(i, j)*2.0f + l(i - 1, j);
      res(i, j) = res(i, j) / (h*h);
    }
  }
  //complex m = (1.0f)/(h*h);
  //res = res * m;
  return res;
}

matriz aux::d2y(matriz& l)
{
  matriz res(l.n, l.m);
  /*matriz d = matriz::tridiagonal(1, -2, 1, l.n);
  matriz t = transpose(l);
  res = d * t;*/
  for(int i = 1; i<l.n-1; i++)
  {
    for(int j = 1; j<l.m-1; j++)
    {
      res(i, j) = l(i, j + 1) - l(i, j)*2.0f + l(i, j - 1);
      res(i, j) = res(i, j) / (h*h);
    }
  }
  return res;
}

void aux::printCoord(ostream &o, matriz &l)
{
  for(int i = 0; i<l.n; i++)
  {
    for(int j = 0; j<l.m; j++)
    {
      /*float x = -10.0f + i*h;
      float y = -10.0f + j*h;*/
      o << i << "\t" << j << "\t" << l(i, j) << endl;
    }
  }
}
