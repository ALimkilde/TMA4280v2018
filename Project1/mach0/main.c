#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

double myarctan(double x,int n){
  double Sn = 0;
  double tmp;
  for (int i = 1; i <= n; i++)
  {
    tmp = 2*i-1;
    Sn += pow(-1.0,i-1) * pow(x,tmp) / tmp;
  }
  return Sn;
}

double utest(){
  std::cout << "Perform Unit Test" << std::endl;
  double x = 0.5;
  double goal = x - pow(x,3)/3 + pow(x,5)/5;
  double val = myarctan(x,3);
  std::cout << "To be zero: " << goal - val << std::endl;
}

double mypi(int n){
  double Sn = 0;
  Sn += myarctan(1.0/5.0,n);
  Sn *= 4;
  Sn -= myarctan(1.0/239.0,n);
  return Sn*4;
}

double vtest(){
  std::cout << "Perform Verification Test" << std::endl;

  std::ofstream file;
  file.open("vtest.txt");

  for(int k = 0; k < 5; k++){
    int n = pow(2,k);
    double pin = mypi(n);
    file << n << ", " << std::abs(pin - M_PI) << "\n";
  }
  file.close();
}

int main(int argc, char ** argv){
  int n = atoi(argv[1]);
  if (n == -1){utest(); return 0;}
  if (n == -2){vtest(); return 0;}
  std::cout << "mypi = " << mypi(n) << std::endl;
  std::cout << "Pi = " << M_PI << std::endl;
  return 0;
}
