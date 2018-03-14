#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <mpi.h>
#include <cstring>
#include <chrono>

double zetafun(int n){
  double Sn = 0;
  for (int i = 1; i <= n; i++)
  {
    Sn += 1.0/(i*i);
  }
  return Sn;
}

void genzeta(int n,double* vec){

  for (int i = 1; i <= n; i++)
  {
    vec[i-1] = 1.0/(i*i);
  }
  return;
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

double mypiMPI(int n, int argc, char ** argv, double * elapsed){
  int nprocs,rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if ((nprocs & (nprocs - 1)) != 0) {
    if (rank == 0){std::cout << "N.o. processes should be a power of two!" << std::endl;}
    MPI_Finalize();
    return 0;
  }
  // Measure time at rank 0
  double t1;
  if (rank == 0){
    t1 = MPI_Wtime();
  }

  // Sum recieved vector
  double sum = 0;
  for (int i = rank+1; i < n; i += nprocs){
      sum += 1.0/(i*i);
  }
  
  // Reduce all partial sums (Here we simply send to root process)
  double Sn = 0;
  MPI_Reduce(&sum, &Sn, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    *elapsed = MPI_Wtime() - t1;
    std::cout <<  sqrt(Sn*6) << std::endl;
    return sqrt(Sn*6);
  }

    return 0;
  }

void timemypiMPI(int argc, char **  argv){
  MPI_Init(&argc, &argv);
  int nprocs,rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::ofstream efile, tfile;
  if (rank == 0){
    std::string name = "res/Errors" + std::to_string(nprocs) + ".txt";
    efile.open(name);

    name = "res/Times" + std::to_string(nprocs) + ".txt";
    tfile.open(name);
  }
  for(int i = 4; i < 25; i++){
    int n = pow(2,i);
    double elapsed;

    double mpi = mypiMPI(n,argc, argv,&elapsed);
    if (rank == 0) {
      tfile << n << " " << elapsed << std::endl;
      efile << n << " " << M_PI-mpi << std::endl;
    }
  }
  if (rank == 0){
    efile.close();
    tfile.close();
  }
  MPI_Finalize();
  return;
}

int main(int argc, char ** argv){
  int n = atoi(argv[1]);
  if (n == -1){utest(); return 0;}
  if (n == -2){vtest(); return 0;}
  if (n == -3){timemypiMPI(argc,argv); return 0;}

  double elapsed;
  MPI_Init(&argc, &argv);
  int nprocs,rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double mpi = mypiMPI(n,argc, argv,&elapsed);
  if (rank == 0) {
    std::cout << mpi << std::endl;
  }

  MPI_Finalize();
  return 0;
}
