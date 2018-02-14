#include "hello_utils.h"
#include "goodbye_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>


int main(int argc, char ** argv)
{
  int ntimes = atoi(argv[1]);
  int in;
  printf("Type 1 for Hello 2 for Goodbye\n");
  std::cin >> in;
  for(int i = 0; i < ntimes;i++){
    switch(in){
      case 1:
        print_hello(); break;
      case 2:
        print_goodbye(); break;
    }
  }
  return 0;
}
