#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>

#include "complex.h"
#include "matriz.h"
#include "aux.h"

using namespace std;

float h;
float ix;
float iy;
float G;
float g;
float omega;
float dt;

matriz f(matriz& l)
{
  matriz lapx = aux::d2x(l);
  matriz lapy = aux::d2y(l);
  matriz laplaciano = lapx + lapy;
  laplaciano = laplaciano * (complex)(-0.5f);

  matriz par2 = l.transform(aux::x2y2);
  matriz par3 = l.transform(aux::gnorm);

  matriz dx = aux::dx(l);
  matriz dy = aux::dy(l);
  dx = dx.transform(aux::my);
  dy = dy.transform(aux::mx);
  matriz par4 = (dy - dx) * (complex(0, omega));

  matriz res = laplaciano + par2 + par3 - par4;
  res = res * (complex)((complex(1, 0)/complex(-g, 1)));
  return res;
}

complex f(int i, int j, complex c, complex c1, complex c2, complex c3, complex c4)
{
    if((i == 0) || (i == 127) || (j == 0) || (j == 127))
      return 0;

      float x = -10.0f + (i-1)*h;
      float y = -10.0f + (j-1)*h;

      complex r;
      complex lapx = (c2 - 2.0f*c + c1)/(h*h);
      complex lapy = (c4 - 2.0f*c + c3)/(h*h);
      complex lap = lapx + lapy;
      lap = lap * (-0.5f);

      complex p2 = c * ((x*x)/2.0f + (y*y)/2.0f);
      float n = c.mod();
      n *= n;
      complex p3 = c * (G*n);

      complex dx = (c2 - c1)/(2.0f*h);
      complex dy = (c4 - c3)/(2.0f*h);
      dx = dx * y;
      dy = dy * x;
      complex p4 = (dy - dx) * complex(0.0f, omega);

      r = lap + p2 + p3 - p4;

      complex m = complex(-g, 1.0f);
      complex resultado = r/m;
      return  resultado;
}


matriz runge_kutta(matriz& l)
{
  matriz r(128, 128);
  matriz psi(128, 128);
  matriz psi2(128, 128);

  #pragma omp parallel for //colapse()
  for(int i = 1; i<l.n-1; i++)
  {
    for(int j = 1; j<l.m-1; j++)
    {
        complex func = f(i, j, l(i,j), l(i-1,j), l(i+1,j), l(i,j-1), l(i,j+1));
        func = func * dt;
        psi(i, j) = l(i,j) + func;
    }
  }
  #pragma omp parallel for
  for(int i = 1; i<l.n-1; i++)
  {
    for(int j = 1; j<l.m-1; j++)
    {
        complex func = f(i, j, psi(i, j), psi(i-1,j), psi(i+1,j), psi(i,j-1), psi(i,j+1));
        func = func * (dt*0.25f);
        complex antigo = psi(i, j);
        antigo = antigo * 0.25f;
        complex lat = l(i, j)*0.75f;
        psi2(i, j) = lat + antigo + func;
    }
  }
  #pragma omp parallel for
  for(int i = 1; i<l.n-1; i++)
  {
    for(int j = 1; j<l.m-1; j++)
    {
      complex func =  f(i, j, psi2(i, j), psi2(i-1,j), psi2(i+1,j), psi2(i,j-1), psi2(i,j+1));
      func = func * (dt*(2.0f/3.0f));
      complex antigo = psi2(i, j);
      antigo = antigo * (2.0f/3.0f);
      complex lat = l(i, j);
      lat = lat*(1.0f/3.0f);
      r(i, j) = lat + antigo + func;
    }
  }
  return r;
}


int main(int argc, char const *argv[]) {

  ifstream cond("params.txt");
  cond >> g >> omega >> G;
  cond.close();

  matriz psi(128, 128);
  h = 20.0f/128.0f;
  ix = -10.0f;
  iy = -10.0f;
  dt = 0.0061f;

  psi.fill(1.0f);

  ifstream c("cond.txt");
  psi.fill(c);
  c.close();

  aux::norma(psi);

  float tmax;
  cout << "Tempo máximo: " << endl;
  cin >> tmax;

  matriz temp(128, 128);
  int iter = 0;

  for(float t = 0.0f; t<tmax; t += dt)
  {
    temp = runge_kutta(psi);
    aux::norma(temp);

    cout << "Iteração: " << iter <<  " <=> T: " << t << endl;

    stringstream f, op;
    f << "dados/t" << iter++ << ".txt";
    ofstream o(f.str().c_str());
    for(int i = 0; i<temp.n; i++)
    {
      for(int j = 0; j<temp.m; j++)
      {
        float x = -10.0f + i*h;
        float y = -10.0f + j*h;
        float mod = temp(i, j).mod();
        o << x << "\t" << y << "\t" << mod*mod << endl;
      }
      o << "\n";
    }

    psi = temp;
  }
  return 0;
}
