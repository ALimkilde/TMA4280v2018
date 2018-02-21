#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

double zetafun(int n){
  double Sn = 0;
  for (int i = 1; i <= n; i++)
  {
    Sn += 1.0/(i*i);
  }
  return Sn;
}

double mypi(int n){
  double z = zetafun(n);
  return sqrt(z*6);
}

void utest(){
  std::cout << "Perform Unit Test" << std::endl;
  double goal = 1.0 + 1.0/4.0 + 1.0/9.0;
  double val = zetafun(3);
  std::cout << "To be zero: " << goal - val << std::endl;
}

double vtest(){
  std::cout << "Perform Verification Test" << std::endl;

  std::ofstream file;
  file.open("vtest.txt");

  for(int k = 1; k < 25; k++){
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
  double Sn = zetafun(n);
  std::cout << "Sn = " << Sn << std::endl;
  std::cout << "Pi^2/6 = " << M_PI * M_PI / 6<< std::endl;
  return 0;
}
