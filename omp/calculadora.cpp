#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

//#include "cuda.cuh"

#include "complex.h"

using namespace std;


int main ()
{
  double d,c;

  complex <double> c1;
  complex <double> c2;
  complex <double> result;

  cout << "Número complexo 1 (introduza parte real e de seguida imaginária)" << endl;
  cin >> c1;
  cout << endl;


  cout << "Selecione a operação que deseja executar:" << endl;
  cout << "1) Soma" << endl;
  cout << "2) Subtração" << endl;
  cout << "3) Multiplicação (entre complexos)" << endl;
  cout << "4) Multiplicação (por escalar)" << endl;
  cout << "5) Divisão (entre complexos)" << endl;
  cout << "6) Divisão (por escalar)" << endl;

  cin >> d;

  if(d!=1 && d!=2 && d!=3 && d!=4 && d!=5 && d!=6)
    {
      cerr << "Inválido!" << endl;
      exit(1);
    }

  if(d==1 || d==2 || d==3 || d==5)
    {
      cout << "Número complexo 2" << endl;
      cin >> c2;
      cout << endl;

      if(d==1)
	result=c1+c2;

      if(d==2)
	result=c1-c2;

      if(d==3)
	result=c1*c2;

      if(d==5)
	result=c1/c2;
    }

  if(d==4 || d==6)
    {
      cout << "Escalar" << endl;;
      cin >> c;
      cout << endl;

      if(d==4)
	result=c*c1;

      if(d==6)
	result=c1/c;
    }

  cout << "Resultado" << endl;
  cout << result << endl;

  return 0;

}
