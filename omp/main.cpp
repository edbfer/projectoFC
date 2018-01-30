#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>

#include "complex.h"
#include "matriz.h"
#include "aux.h"
#include "cuda.cuh"

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
  //return dx;
}

complex f(int i, int j, complex c, complex c1, complex c2, complex c3, complex c4)
{
  /*matriz res(l.n, l.m);
  #pragma omp parallel for
  for(int i = 1; i<l.n-1; i++)iterate
  {
    for(int j = 1; j<l.m-1; j++)
    {*/

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

  //return res;
}

matriz runge_kutta(matriz& l)
{
  matriz fl = f(l);
  fl = fl * (complex)dt;
  matriz psi1 = l + fl;
  matriz fpsi1 = f(psi1);
  fpsi1 = fpsi1 * complex(0.25f*dt, 0);
  matriz psi2 = (l*complex(0.75f, 0));
  psi2 = psi2 + (psi1*complex(0.25f, 0));
  psi2 = psi2 + fpsi1;// + (psi1*(complex)0.25f) + (fpsi1*(complex)(0.25f*dt));
  matriz fpsi2 = f(psi2);
  fpsi2 = fpsi2 *(complex((2.0f/3.0f)*dt, 0));
  matriz res = (l*(complex((1.0f/3.0f), 0)));
  res = res + (psi2*(complex((2.0f/3.0f), 0)));
  res = res + fpsi2;
  return res;
  //return fl;
}

matriz engolindo_sapos(matriz& l)
{
  matriz r(128, 128);
  matriz psi(128, 128);
  matriz psi2(128, 128);
  //#pragma omp parallel for
  //ofstream p1("psi1.txt");
  #pragma omp parallel for
  for(int i = 1; i<l.n-1; i++)
  {
    for(int j = 1; j<l.m-1; j++)
    {
        complex func = f(i, j, l(i,j), l(i-1,j), l(i+1,j), l(i,j-1), l(i,j+1));
        func = func * dt;
        psi(i, j) = l(i,j) + func;
    }
  }
  //p1 << psi;
  //ofstream p2("psi2.txt");
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
  //p2 << r;
  //ofstream p3("psi3.txt");
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
  //p3 << r;
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

  //cuda_setup(G, g, omega, dt, h);

  psi.fill(1.0f);

  //psi.fill(complex(1.0f, 1.0f));
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
    temp = engolindo_sapos(psi);
    aux::norma(temp);
    //temp = runge_kutta(psi);
    //temp = cuda_doround(psi);
    /*ofstream ouch("lloo.txt");
    aux::printCoord(ouch, temp);
      int lel; cin >> lel; cout << "running" << endl;*/

    /*float no = 0.;// = aux::norma(temp);
    for(int i = 0; i<temp.n; i++)
    {
      for(int j = 0; j<temp.m; j++)
      {
        float mod = temp(i, j).mod();
        no = no + mod*mod*h*h;
      }
    }
    float multi = sqrt(no);
    //#pragma omp parallel for
    for(int i = 0; i<temp.n; i++)
    {
      for(int j = 0; j<temp.m; j++)
      {
        temp(i, j) = temp(i, j) / multi;
      }
    }*/

    cout << "Iteração: " << iter << /*" <=> Norma: " << no <<*/ " <=> T: " << t << endl;

    stringstream f, op;
    f << "dados/t" << iter++ << ".txt";
    //op << "dados/gnuplot" << iter << ".txt";
    ofstream o(f.str().c_str());
    //ofstream os(op.str().c_str());
    for(int i = 0; i<temp.n; i++)
    {
      for(int j = 0; j<temp.m; j++)
      {
        float x = -10.0f + i*h;
        float y = -10.0f + j*h;
        float mod = temp(i, j).mod();
        //os << i << "\t" << j << "\t" << temp(i, j) << endl;
        o << x << "\t" << y << "\t" << mod*mod << endl;
      }
      o << "\n";
    }

    psi = temp;
  }

  /*
  cout << res;
  int lel; cin >> lel; cout << "running" << endl;
  */
  return 0;
}
