#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include <iostream>


void VecSum(double *x, double *y, double *z, double *a, int size)
{
  for(int i = 0; i < size; i++){
    z[i] = *a*x[i] + y[i];
  }
};

void VecMatSum(double *x, double *y, double *z, double **A, int size)
{
  for(int i = 0; i < size; i++){
    double tmp = x[i];
    for (int j = 0; j < size; j++){
      tmp += A[i][j]*y[j];
    }
    z[i] = tmp;
  }
};

void DotProd(double *x, double *y, double *a, int size)
{
  *a = 0;
  for(int i = 0; i < size; i++){
    *a += x[i]*y[i];
  }
}

void PrintMat(double **A,int size)
{
  std::cout << "Print a Matrix: " << std::endl;
  for(int i = 0; i < size; i++)
  {
    for(int j = 0; j < size; j++)
    {
      std::cout << A[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void PrintVec(double *x,int size)
{
    for(int i = 0; i < size; i++){
      std::cout << x[i] << " ";
    }
    std::cout << std::endl;
}
