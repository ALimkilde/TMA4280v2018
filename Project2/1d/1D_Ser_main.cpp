#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <mpi.h>
#include <cstring>
#include <tma/types.h>
#include <tma/mesh/mesh.h>
#include <tma/cell/interval.h>
#include <tma/elements/lagrange.h>
#include <tma/mapping.h>
#include <tma/quadrature.h>
// namespace tma

void printvec(double* vec, int size){
    for(int i = 0; i < size; i++){
      std::cout << vec[i] << std::endl;
    }
}

class Vector{
private:
  double* data;
  int* ghosts;
public:
  int size,offset,n_ghosts;
  Vector(int a, int b){
    size = a;
    offset = b;
    data = new double[size];
  }
  double& operator[] (int x) {
          return data[x];
  }
  void set_ghosts(int n, int* arr){
    n_ghosts = n;
    ghosts = new int[n];
    for(int i = 0; i < n; i++){
      ghosts[i] = arr[i];
    }
  }
  void set_ghosts(int i){
    n_ghosts = 1;
    ghosts = new int[1];
    ghosts[0] = i;
  }
  int get_n_ghosts(){
    return n_ghosts;
  }
  int get_ghosts(int* out){
    for(int i = 0; i < n_ghosts; i++){
      out[i] = ghosts[i];
    }
  }
};

class CRS{
public:
  int* O;
  int* C;
  double* V;

  CRS(){

  }
  CRS(int m, int z){
    set_size_and_offset(m,z);
  }
  void set_size_and_offset(int m, int z){
    O = new int[m];
    C = new int[z];
    V = new double[z];
  }
  void insert(int k, int l, double val){

  }
};

class FullSquareMatrix{
public:
  int size_;
  double* vals_;

  FullSquareMatrix(int size){
    vals_ = new double[size*size];
    size_ = size;
    for(int i = 0; i < size; i++){
      for(int j = 0; j < size; j++){
        vals_[i + size*j] = 0;
      }
    }
  }

  void insert_add(int k, int l, double val){
    vals_[k+size_*l] += val;
  }

  void dump(){
    for(int i = 0; i < size_; i++){
      std::cout << vals_[i];
      for(int j = 1; j < size_; j++){
        std::cout << "\t" << vals_[i + size_*j];
      }
      std::cout << std::endl;
    }
  }

  void mult(double* b, double* x){
    for(int i = 0; i < size_; i++){
      for(int j = 0; j < size_; j++){
        x[i] += vals_[i + size_*j]*b[j];
        // std::cout << "mvprod at " << i << ", " << j << ":" << vals_[i + size_*j]*b[j] << std::endl;
      }
    }
  }

};

class Matrix{
public:
  CRS Diagonal;
  CRS OffDiagonal;

  int offset,n_rows;
  Matrix(int a, int b){
    n_rows = a;
    offset = b;
  }
};

class Function{
public:
  double eval(double x){
    return cos(2*M_PI*x);
  }
};

struct UnitIntervalMesh{
  tma::mesh<tma::interval,1> M;
  int rank_,nprocs_;
  uint ncells_,nverts_;

  UnitIntervalMesh(uint n, int rank, int nprocs):
  rank_(rank), nprocs_(nprocs),
  ncells_(n/nprocs),
  nverts_(n/nprocs+1),
  M(n/nprocs,n/nprocs+1)
  {
    int m = M.topo().nverts();
    // Set the topology...
    for(int i = 0; i < ncells_; i++){
      *M.topo()(i) = i;
      *(M.topo()(i)+1) = i+1;
    }

    // Set the geometry
    for(int i = 0; i < nverts_; i++){
      *M.geom()(i) = double((ncells_*rank) + i)/n;
    }
  }

  uint local_to_global(uint i){
    return uint(i + ncells_*rank_);
  }
  int owner_global(uint j){
    return j/nverts_;
  }
  int owner_local(uint i){
    return owner_global(local_to_global(i));
  }

  void assemble_stiffness_matrix(FullSquareMatrix mat){

    for(int cell = 0; cell < ncells_; cell++){
      int k = *M.topo()(cell);
      int l = *(M.topo()(cell)+1);

      // Gauss quadrature
      tma::QR qr;
      uint num = qr.interval_.num_from_deg(0);
      double * abscissas = new double[num];
      double * weights = new double[num];
      qr.interval_.get_quadrature(0,abscissas,weights);
      for(int i = 0; i < num; i++){

        double** dphis = new double*[1];
        dphis[0] = new double[2];

        double* x = new double[1];
        x[0] = abscissas[i];
        double w = weights[i];

        double** x2 = new double*[2];
        x2[0] = new double[1];
        x2[1] = new double[1];

        x2[0][0] = *M.geom()(k);
        x2[1][0] = *M.geom()(l);

        tma::P1 p1;
        p1.interval_.d(x,dphis);
        tma::Jacobian<1> J;
        J(x2);
        double jac = J.K[0][0];
        double val = w*jac*dphis[0][0]*dphis[0][1]*jac*J.det;
        mat.insert_add(k,l,val);
        mat.insert_add(l,k,val);

        double val2 = w*jac*dphis[0][0]*dphis[0][0]*jac*J.det;
        mat.insert_add(k,k,val2);

        double val3 = w*jac*dphis[0][1]*dphis[0][1]*jac*J.det;
        mat.insert_add(l,l,val3);
      }

    }
  }

  void assemble_mass_matrix(FullSquareMatrix mat){

    for(int cell = 0; cell < ncells_; cell++){
      int k = *M.topo()(cell);
      int l = *(M.topo()(cell)+1);

      // Gauss quadrature
      tma::QR qr;
      uint num = qr.interval_.num_from_deg(2);
      double * abscissas = new double[num];
      double * weights = new double[num];
      qr.interval_.get_quadrature(2,abscissas,weights);
      for(int i = 0; i < num; i++){

        double* phis = new double[2];

        double* x = new double[1];
        x[0] = abscissas[i];
        double w = weights[i];

        double** x2 = new double*[2];
        x2[0] = new double[1];
        x2[1] = new double[1];

        x2[0][0] = *M.geom()(k);
        x2[1][0] = *M.geom()(l);

        tma::P1 p1;
        p1.interval_(x,phis);
        tma::Jacobian<1> J;
        J(x2);
        double val = w*phis[0]*phis[1]*J.det;
        mat.insert_add(k,l,val);
        mat.insert_add(l,k,val);

        double val2 = w*phis[0]*phis[0]*J.det;
        mat.insert_add(k,k,val2);

        double val3 = w*phis[1]*phis[1]*J.det;
        mat.insert_add(l,l,val3);
      }

    }
  }

  void interpolate(Function fun, double * vec){
    for(int i = 0; i < nverts_; i++){
      vec[i] = fun.eval(*M.geom()(i));
    }
  }
};

int main(int argc, char ** argv){
  int n = atoi(argv[1]);

  if (n%2 != 0){
    std::cout << "N should be a power of 2..." << std::endl;
    return 0;
  }

  MPI_Init(&argc, &argv);
  int nprocs,rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  UnitIntervalMesh UM = UnitIntervalMesh(n,rank,nprocs);

  if (rank == 0){
    FullSquareMatrix K = FullSquareMatrix(n+1);
    // K.dump();
    UM.assemble_stiffness_matrix(K);
    std::cout << "Stiffness Matrix:" << std::endl;
    K.dump();

    FullSquareMatrix M = FullSquareMatrix(n+1);
    UM.assemble_mass_matrix(M);
    std::cout << "Mass Matrix:" << std::endl;
    M.dump();

    Function cosfun;
    double * vec = new double[n+1];
    UM.interpolate(cosfun,vec);
    std::cout << "vec: " << std::endl;
    printvec(vec,n+1);
    double * Kvec = new double[n+1];
    K.mult(vec,Kvec);
    std::cout << "\n Kvec : " << std::endl;
    printvec(Kvec,n+1);

    double * Mvec = new double[n+1];
    M.mult(vec,Mvec);
    std::cout << "\n Mvec: " << std::endl;
    printvec(Mvec,n+1);

    double * r = new double[n+1];
    double rnorm = 0;
    for(int i = 0; i < n+1; i++){
      r[i] = Kvec[i] - 4*M_PI*M_PI*Mvec[i];
      rnorm += r[i]*r[i];
    }

    std::cout << "\n\n residual: " << std::endl;
    printvec(r,n+1);

    std::cout << "\n\n .... rnorm = " << sqrt(rnorm) << std::endl;





  }
	// MeshPartition mesh(rank,nprocs,n);
  // //mesh.print_range();
  // mesh.generate_geometry();
  // mesh.print_geometry();

  if (rank == 0){
    // UM.M.topo().dump();

    //int i = 2;
    //std::cout << "Owner of " << i << " is " << mesh.get_owner(i) << std::endl;
    std::cout << "=== Done ===\n" << std::endl;
  }

  MPI_Finalize();
  return 0;
}
