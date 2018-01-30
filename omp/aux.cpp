#include "aux.h"

extern float h;
extern float G;

complex aux::x2y2(int i, int j, complex c)
{
  float x = -10.0f + i*h;
  float y = -10.0f + j*h;
  complex m = ((x*x) + (y*y))/2.0f;
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
  return c * (~c);
}

matriz aux::dx(matriz& l)
{
  matriz res(l.n, l.m);
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
  for(int i = 1; i<l.n-1; i++)
  {
    for(int j = 1; j<l.m-1; j++)
    {
      res(i, j) = l(i + 1, j) - l(i, j)*2.0f + l(i - 1, j);
      res(i, j) = res(i, j) / (h*h);
    }
  }
  return res;
}

matriz aux::d2y(matriz& l)
{
  matriz res(l.n, l.m);
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
      o << i << "\t" << j << "\t" << l(i, j) << endl;
    }
  }
}

void aux::norma(matriz& l)
{
  float n = 0;
  for(int i = 0; i<l.n; i++)
  {
    for(int j = 0; j<l.m; j++)
    {
      float modu = l(i, j).mod();
      n = n + modu*modu*h*h;
    }
  }
  cout << "Norma: " << n << endl;
  float multi = sqrt(n);

  for(int i = 0; i<l.n; i++)
  {
    for(int j = 0; j<l.m; j++)
    {
      l(i, j) = l(i, j) / multi;
    }
  }
}
