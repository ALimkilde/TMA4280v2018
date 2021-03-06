#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <mpi.h>
#include <cstring>
#include <chrono>
#include <omp.h>

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

void genmach(int n,double x, double* vec){

  #pragma omp parallel for
  for (int i = 1; i <= n; i++)
  {
    double tmp = 2*i-1;
    vec[i-1] = pow(-1.0,i-1) * pow(x,tmp) / tmp;
  }
  return;
}

void getvec(int n, double* vec){
  double* vec1 = new double[n];
  genmach(n,1.0/5.0,vec1);

  double* vec2 = new double[n];
  genmach(n,1.0/239.0,vec2);

  #pragma omp parallel for
  for (int i = 0; i < n; i++){
    vec[i] = vec1[i]*4 - vec2[i];
  }
  delete [] vec1;
  delete [] vec2;

}

double mypi(int n){
  double Sn = 0;
  Sn += myarctan(1.0/5.0,n);
  Sn *= 4;
  Sn -= myarctan(1.0/239.0,n);
  return Sn*4;
}

double utest(){
  std::cout << "Perform Unit Test" << std::endl;
  double x = 0.5;
  double goal = x - pow(x,3)/3 + pow(x,5)/5;
  double val = myarctan(x,3);
  std::cout << "To be zero: " << goal - val << std::endl;
}

double vtest(){
  std::cout << "Perform Verification Test" << std::endl;

  std::ofstream file;
  file.open("vtest.txt");

  for(int k = 1; k < 5; k++){
    int n = pow(2,k);
    double pin = mypi(n);
    file << n << ", " << std::abs(pin - M_PI) << "\n";
  }
  file.close();
}

double mypiMPI(int n, int argc, char ** argv){
  int nprocs,rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if ((nprocs & (nprocs - 1)) != 0) {
    if (rank == 0){std::cout << "N.o. processes should be a power of two!" << std::endl;}
    MPI_Finalize();
    return 0;
  }
  // Distribute a vector with elements to other processes
  if (rank == 0){
    // Generate vector elements
    double* vec = new double[n];
    getvec(n,vec);

    // Send to other ranks (nprocs -1 in total)
    int nm1 = nprocs-1;
    int iter = 0;
    for (int i = 1; i < nprocs; i++){
      int size = n/nm1;
      if (i <= n%nm1){size++;}

      // Copy vector
      double* sendvec = new double[size];
      std::memcpy(sendvec,vec+iter,size*sizeof(double));
      iter += size;
      //std::cout << "Rank " << i << " was sent " << size << " elements" << std::endl;

      // Send size and vector
      MPI_Send(&size,1,MPI_INT, i, 0,MPI_COMM_WORLD);
      MPI_Send(sendvec,size,MPI_DOUBLE,i,1,MPI_COMM_WORLD);
      delete [] sendvec;
    }
    delete [] vec;
  }

  // For processes that are not the root
    if (rank > 0){
      // Recieve size and vector
      int size;
      MPI_Recv(&size,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      //std::cout << "... Rank " << rank << " recieved " << size << " elements" << std::endl;

      double* recvec = new double[size+10];
      MPI_Recv(recvec,size,MPI_DOUBLE,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

      // Sum recieved vector
      double sum = 0;
      #pragma omp parallel
      {
        double loc_res = 0;
      #pragma omp for
        for (int i = 0; i < size; i++){
          loc_res += recvec[i];
        }
      #pragma omp critical
        sum += loc_res;
      }
      delete [] recvec;
      // Reduce all partial sums (Here we simply send to root process)
      double tmp;
      MPI_Reduce(&sum,&tmp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    }

  // At root process; reduce partial sums
    if (rank == 0){
      double Sn;
      double tmp = 0;
      MPI_Reduce(&tmp, &Sn, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (rank == 0){
        //std::cout << "Sn = " << Sn << std::endl;
        //std::cout << "Pi^2/6 = " << M_PI * M_PI / 6<< std::endl;
        return Sn*4;
      }
    }
    return 0;
  }

void timemypiMPI(int argc, char **  argv){
  MPI_Init(&argc, &argv);
  int nprocs,rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int m;
  if (argc > 1){
    m = atoi(argv[2]);
    omp_set_dynamic(0);     // Disable dynamic teams
    omp_set_num_threads(m); // Use m threads
    }

  std::ofstream efile, tfile;
  if (rank == 0){
    std::cout << "Using " << m << " omp threads" << std::endl;
    std::string name = "res2/Errors_np" + std::to_string(nprocs) + "_m" + std::to_string(m) + ".txt";
    efile.open(name);

    name = "res2/Times_np" + std::to_string(nprocs) + "_m" + std::to_string(m) + ".txt";
    tfile.open(name);
  }
  for(int i = 4; i < 20; i++){
    int n = pow(2,i);
    auto start = std::chrono::high_resolution_clock::now();

    if (rank == 0) {
      start = std::chrono::high_resolution_clock::now();
    }
    double mpi = mypiMPI(n,argc, argv);
    if (rank == 0) {
      auto stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = stop - start;
      tfile << n << " " << elapsed.count() << std::endl;
      efile << n << " " << mpi - M_PI << std::endl;
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

  MPI_Init(&argc, &argv);
  int nprocs,rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double mpi = mypiMPI(n,argc, argv);
  if (rank == 0) {
    std::cout << mpi << std::endl;
  }

  MPI_Finalize();
  return 0;
}
