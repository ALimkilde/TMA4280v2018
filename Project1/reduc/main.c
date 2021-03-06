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

double mypiMPI(int n, int argc, char ** argv, double* elapsed){
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
    genzeta(n,vec);

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

  double sum = 0;
  double Sn;
  // For processes that are not the root
    if (rank > 0){
      // Recieve size and vector
      int size;
      MPI_Recv(&size,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      //std::cout << "... Rank " << rank << " recieved " << size << " elements" << std::endl;

      double* recvec = new double[size+10];
      MPI_Recv(recvec,size,MPI_DOUBLE,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

      // Sum recieved vector
      for (int i = 0; i < size; i++){
        sum += recvec[i];
      }
      delete [] recvec;
      // Reduce all partial sums (Here we simply send to root process)
    }

    // Measure time
    double t1;
    if (rank == 0){
      t1 = MPI_Wtime();
    }

    //MPI_Allreduce(&sum,&Sn,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for (int k = nprocs/2; k >= 1; k /= 2){
      //std::cout << "k = "<< k << std::endl;
      int r = (rank + k)%(k*2);
      //std::cout << "r = "<< r << std::endl;
      std::cout << "Rank " << rank << " trying to send to " << r << std::endl;
      double in;
      MPI_Sendrecv(&sum,1,MPI_DOUBLE,r,1,&in,1,MPI_DOUBLE,r,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      sum += in;
    }

    if (rank == 0){
      double t2 = MPI_Wtime();
      *elapsed = t2-t1;
    }

  // At root process; reduce partial sums
    if (rank == 0){
      //double Sn;
      //double tmp = 0;
      //MPI_Reduce(&tmp, &Sn, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (rank == 0){
        //std::cout << "Sn = " << Sn << std::endl;
        //std::cout << "Pi^2/6 = " << M_PI * M_PI / 6<< std::endl;
        return sqrt(sum*6);
      }
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

  MPI_Init(&argc, &argv);
  int nprocs,rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double elapsed;
  double mpi = mypiMPI(n,argc, argv, &elapsed);
  if (rank == 0) {
    std::cout << mpi << std::endl;
    std::cout << "Elapsed time " << elapsed << std::endl;
  }

  MPI_Finalize();
  return 0;
}
