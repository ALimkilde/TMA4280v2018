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

void DotProd(double *x, double *y, double *a, int size)
{
  *a = 0;
  for(int i = 0; i < size; i++){
    *a += x[i]*y[i];
  }
}
