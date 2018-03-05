#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>


int main(int argc, char ** argv)
{
  double x [] = {1,0,0};
  double y [] = {0,1,1};
  double z [3];
  double a = 2;
  int size = sizeof(x)/sizeof(*x);

  VecSum(x,y,z,&a,size);

  for(int i = 0; i < size; i++){
    std::cout << z[i] << " ";
  }
  std::cout << std::endl;

  DotProd(x,y,&a,size);
  std::cout << "a = " << a << std::endl;

  double TMP [size][size] = {{2,0,0},{0,3,0},{0,0,4}};
  double **A;
  A = new double *[size];
  for(int i = 0; i < size; i++)
  {
    A[i] = TMP[i];
  }

  PrintMat(A,size);

  VecMatSum(x,y,z,A,size);

  PrintVec(z,size);

  return 0;
}
