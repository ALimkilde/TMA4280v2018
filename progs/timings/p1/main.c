#include "../../VecOp/utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>


int main(int argc, char ** argv)
{
  int n = atoi(argv[1]);

  double a[n];
  double b[n];
  double z[n];

  srand(time(0));
  double gamma = rand()%10;

  for(int i=0; i<n; i++)
  {
    a[i] = (rand()%100);
    b[i] = (rand()%100);
  }

  double tstart = clock();
  VecSum(a,b,z,&gamma,n);
  double t = clock();
  t -= tstart;
  std::cout << "Time spend: " << (t)/CLOCKS_PER_SEC <<  "s" << std::endl;
  return 0;
}
