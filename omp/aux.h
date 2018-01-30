#ifndef _AUX_H
#define _AUX_H 1

#include "matriz.h"
#include "complex.h"

namespace aux
{

  complex conjugate(int i, int j, complex c);
  complex gnorm(int i, int j, complex c);
  complex x2y2(int i, int j, complex c);
  complex mx(int i, int j, complex c);
  complex my(int i, int j, complex c);

  matriz dx(matriz& l);
  matriz dy(matriz& l);
  matriz d2x(matriz& l);
  matriz d2y(matriz& l);

  void printCoord(ostream& o, matriz& l);
  void norma(matriz& l);
}

#endif
