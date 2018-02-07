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

  VecSum(x,y,z,&a,sizeof(x)/sizeof(*x));

  for(int i = 0; i < sizeof(x)/sizeof(*x); i++){
    std::cout << z[i] << " ";
  }
  std::cout << std::endl;

  DotProd(x,y,&a,sizeof(x)/sizeof(*x));
  std::cout << "a = " << a << std::endl;


  return 0;
}
