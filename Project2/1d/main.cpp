#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
// #include <math.h>
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

int MAX = 2e7;

void printvec(double* vec, int size){
    for(int i = 0; i < size; i++){
      std::cout << vec[i] << std::endl;
    }
}


class Vector{
public:
  double* data;
  int* ghosts,*s_with;
  int size,offset,n_ghosts;
  Vector(int a, int b){
    size = a;
    offset = b;
    data = new double[size];
    s_with = new int[size];
    for (int i = 0; i < size; i++){data[i] = 0;}
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
  void set_shared_with(int rank, int nprocs){
    for(int i = 0; i < size; i++){
      s_with[i] = rank;
    }
    if (rank > 0){
      s_with[0]--;
    }
    if (rank < nprocs-1){
      s_with[size-1]++;
    }
  }
  int shared_with(int i){
    return s_with[i];
  }
  void insert_add(int i, double val){
    data[i] += val;
  }
  void dump(){
    printvec(data,size);
  }
  void insert_from_global(int i_g,double val){
    int i_l = i_g - offset;
    // std::cout << i_g << " was converted to " << i_l << std::endl;
    data[i_l] = val;
  }
  void insert_add_from_global(int i_g,double val){
    int i_l = i_g - offset;
    // std::cout << i_g << " was converted to " << i_l << std::endl;
    data[i_l] += val;
  }
  double norm(){
    double local_sum = 0;
    for(int i = 0; i < size; i++){
      local_sum += data[i]*data[i];
    }
    double global_sum;
    MPI_Allreduce(&local_sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return sqrt(global_sum);
  }
  double inner_prod(Vector vec){
    double local_sum = 0;
    for(int i = 0; i < size; i++){
      local_sum += data[i]*vec.data[i];
    }
    double global_sum;
    MPI_Allreduce(&local_sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return global_sum;
  }
  void add_vec(double alpha,Vector vec){
    for(int i = 0; i < size; i++){
      data[i] += alpha*vec[i];
    }
  }
  void add_vec(double alpha,Vector vec, Vector out){
    for(int i = 0; i < size; i++){
      out.data[i] = alpha*vec[i];
    }
  }
  void copy(Vector vec){
    for (int i = 0; i < size; i++){
      data[i] = vec.data[i];
    }
  }
  void mult_with_const(double c){
    for (int i = 0; i < size; i++){
      data[i] = c*data[i];
    }
  }
  void set_to_const(double c){
    for (int i = 0; i < size; i++){
      data[i] = c;
    }
  }
  double max_norm(){
    double local_max = 0;
    for(int i = 0; i < size; i++){
      if (local_max < std::abs(data[i])){
        local_max = std::abs(data[i]);
      }
    }
    std::cout << "local max " << local_max << std::endl;
    double global_max;
    MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    return global_max;
  }
};


class CRS{
public:
  int* O;
  int* C;
  double* V;

  int m_,z_;

  CRS(){

  }
  CRS(int m, int z){
    m_ = m;
    z_ = z;
    set_size_and_offset(m,z);
  }
  void set_size_and_offset(int m, int z){
    O = new int[m];
    C = new int[z];
    V = new double[z];
    O[0] = 0;
    for (int i = 1; i < m; i++){
      O[i] = 0;
    }
    for (int i = 0; i < z; i++){
      C[i] = MAX;
      V[i] = 0;
    }
  }
  void insert(int r, int c, double val){
    // Insert value and column index and push the rest...
    int i = O[r];

    int imax;
    if (r+1 < m_){
      imax = O[r+1];
    } else {
      imax = z_;
    }


    while(C[i] < c && i < imax){ i++; }

    double tmpval = V[i];
    int tmpcol = C[i];
    V[i] = val;
    C[i] = c;

    for(int k = i+1; k < z_; k++){
      double tmpval2 = V[k];
      int tmpcol2 = C[k];
      V[k] = tmpval;
      C[k] = tmpcol;
      tmpval = tmpval2;
      tmpcol = tmpcol2;
    }

    // Update offsets
    for(int k = r+1; k < m_; k++){ O[k]++; }

  }
  int exists(int i, int j){
    int c = O[i];
    int cn = z_;

    if (i+1 < m_){
      cn = O[i+1];
    }

    while(c < cn){
      if (C[c] == j){
        return c;
      }
      c++;
    }
    return -1;
  }
  double operator()(int i, int j){
    int ex = exists(i,j);
    if (ex != -1){
      return V[ex];
    } else {
      return 0;
    }
  }
  template<class vec>
  void mult(vec x, vec y){
    for(int r = 0; r < m_; r++){
      int c = O[r];
      int cn = z_;
      y[r] = 0;
      if(r+1 < m_){cn = O[r+1];}
      while(c < cn && C[c] < MAX){
        y[r] += V[c]*x[C[c]];
        c++;
      }
    }
  }
  void insert_add(int i, int j, double val){
    int ex = exists(i,j);
    if (ex == -1){
      insert(i,j,val);
    } else {
      V[ex] += val;
    }
  }
  void add_matrix(int alpha, CRS M){
    for (int i = 0; i < z_; i++){
      V[i] += alpha*M.V[i];
    }
  }
  void dump_full(){
    // Get max col num
    int cmax = 0;
    for(int i = 0; i < z_; i++){
      if (C[i] > cmax && C[i] < MAX){cmax = C[i];}
    }

    int iter = 0;
    for(int r = 0; r < m_; r++){
      int Orp1;
      if (r+1 < m_){
        Orp1 = O[r+1];
      } else {
        Orp1 = z_;
      }
      for(int c = 0; c <= cmax; c++){
        if(c == C[iter] && iter < Orp1){
          std::cout << V[iter] << "\t";
          iter++;
        } else {
          std::cout << "- \t";
        }
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  void dump(){
    std::cout << "C: \t";
    for(int i = 0; i < z_; i++){
      if(C[i] == MAX){ std::cout << "- "; break;}
      std::cout << C[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "V: \t";
    for(int i = 0; i < z_; i++){
      std::cout << V[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "O: \t";
    for(int i = 0; i < m_; i++){
      std::cout << O[i] << " ";
      if(O[i] == MAX){ std::cout << "- "; break;}
    }
    std::cout << std::endl;
    std::cout << std::endl;
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
  int rank_, nprocs_;
  Matrix(int a, int b,int rank, int nprocs){
    n_rows = a;
    offset = b;
    rank_ = rank;
    nprocs_ = nprocs;
  }
  void initialize_1d(int nverts_global,int nverts_local){
    int m = nverts_local;
    int zd = nverts_local + 2*(nverts_local-1);
    int zm = 3*nverts_local - zd;

    Diagonal = CRS(m,zd);
    OffDiagonal = CRS(m,zm);
  }
  void dump(){
    std::cout << "D:"<< std::endl;
    Diagonal.dump_full();
    std::cout << "M:" << std::endl;
    OffDiagonal.dump_full();
  }
  void insert_add_from_global(int i_g, int j_g, double val){
    int i_l = i_g-offset;
    int j_l;
    // std::cout << "the offset here is " << offset << std::endl;
    // std::cout << "\t " << i_g << ", " << "is converted to row" << i_l << std::endl;
    if (j_g < offset){
      j_l = j_g;
      OffDiagonal.insert_add(i_l,j_l,val);
    } else if (j_g >= offset + n_rows){
      j_l = j_g - n_rows;
      OffDiagonal.insert_add(i_l,j_l,val);
    } else {
      j_l = j_g - offset;
      Diagonal.insert_add(i_l,j_l,val);
    }
    // std::cout << "\t " << j_g << ", " << "is converted to col" << j_l << std::endl;
    return;
  }
  void mult(Vector x, Vector y){
    // Multiply with Diagonal Matrix.
    Diagonal.mult(x,y);    // Send contributions
    for(int i = 0; i < x.size; i++){
        int rank = x.shared_with(i);
        // std::cout << "at rank " << rank_ << ", "<< i << " is shared with " << rank << std::endl;
        int tag = 0;
        if (i > 0){tag = 1;}
        if (rank != rank_){
          MPI_Send(&x[i],1,MPI_DOUBLE,rank,tag,MPI_COMM_WORLD);
          // std::cout << "rank " << rank_ << " sends " << tag << " to " << rank << std::endl;
        }
    }

    // Multiply with OffDiagonal Matrix.
    for(int row = 0; row < n_rows; row++){
      int c1 = OffDiagonal.O[row];
      int c2 = OffDiagonal.z_;
      if (row+1 < n_rows){c2 = OffDiagonal.O[row+1];}
      for(int o = c1; o < c2; o++){
        int rank = rank_;
        double col = OffDiagonal.C[o];
        if (col >= MAX || col < 0 ){ break; }
        // std::cout << "col = " << col << std::endl;
        int tag;
        if (col >= offset){rank++; tag = 0;}
        if (col < offset){rank--; tag = 1;}
        double val;
        // std::cout << "rank " << rank_ << " tries to recv " << tag << " from " << rank << std::endl;
        MPI_Recv(&val,1,MPI_DOUBLE,rank,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        // std::cout << "rank " << rank_ << " recieved " << tag << " from " << rank << std::endl;
        y[row] += val*OffDiagonal.V[o];
      }
    }
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
  int rank_,nprocs_,n_;
  uint ncells_,nverts_;
  int * local_to_DoF;
  int owned_verts;

  UnitIntervalMesh(uint n, int rank, int nprocs):
  rank_(rank), nprocs_(nprocs), n_(n),
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

    // Count owned verticies
    owned_verts = 0;
    for(int i = 0; i < nverts_; i++){
      if (is_owner_local_vert(i)){
        owned_verts++;
      }
    }

    local_to_DoF = new int[owned_verts];
    int iter = 0;
    for(int i = 0; i < nverts_; i++){
      if (is_owner_local_vert(i)){
        local_to_DoF[i] = iter;
        iter++;
      }
    }

  }

  uint local_to_global_cell(uint i){
    return uint(i + ncells_*rank_);
  }
  uint local_to_global_vert(uint i){
    return uint(i + (nverts_-1)*rank_);
  }
  int cell_to_vert_global(int cell, int i){
    int c_g = local_to_global_cell(cell);
    if (i == 0){
      return c_g;
    } else if (i == 1){
      return c_g + 1;
    }
  }
  int owner_global_cell(uint j){
    return j/ncells_;
  }
  int owner_local_cell(uint i){
    return owner_global_cell(local_to_global_cell(i));
  }
  bool is_owner_local_cell(uint i){
    return (owner_local_cell(i) == rank_);
  }
  bool is_owner_local_vert(uint i){
    return (owner_local_vert(i) == rank_);
  }
  int owner_local_vert(uint i){
    if (i == 0 && rank_ > 0){
      return rank_-1;
    } else {
      return rank_;
    }
  }
  bool is_shared_local(int i){
    if (rank_ == 0){
      return i == nverts_-1;
    } else if (rank_ == nprocs_-1){
      return i == 0;
    } else {
      return i == nverts_-1 || i == 0;
    }
  }
  int shared_with(int i){
    if (i == 0){
      return rank_-1;
    } else if (i == nverts_-1){
      return rank_+1;
    }
  }
  template<class MatClass>
  void assemble_stiffness_matrix(MatClass mat){

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

  template<class MatClass>
  void assemble_mass_matrix(MatClass mat){

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

  double get_stiffness_element(int cell, int m, int n){
    int k = *M.topo()(cell);
    int l = *(M.topo()(cell)+1);

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
      if (m != n){
        return w*jac*dphis[0][0]*dphis[0][1]*jac*J.det;
      } else if (m == 0){
        return w*jac*dphis[0][0]*dphis[0][0]*jac*J.det;
      } else if (m == 1){
        return w*jac*dphis[0][1]*dphis[0][1]*jac*J.det;
      }

    }
  }

  double get_mass_element(int cell, int m, int n){
    int k = *M.topo()(cell);
    int l = *(M.topo()(cell)+1);

    tma::QR qr;
    uint num = qr.interval_.num_from_deg(2);
    double * abscissas = new double[num];
    double * weights = new double[num];
    qr.interval_.get_quadrature(2,abscissas,weights);
    double out = 0;
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
      double jac = J.K[0][0];
      if (m != n){
        out += w*phis[0]*phis[1]*J.det;
      } else if (m == 0){
        out += w*phis[0]*phis[0]*J.det;
      } else if (m == 1){
        out += w*phis[1]*phis[1]*J.det;
      }

    }
    return out;
  }
  double get_rhs_element(int cell, int m, Function fun){
    int k = *M.topo()(cell);
    int l = *(M.topo()(cell)+1);

    tma::QR qr;
    uint num = qr.interval_.num_from_deg(9);
    double * abscissas = new double[num];
    double * weights = new double[num];
    qr.interval_.get_quadrature(9,abscissas,weights);
    double out = 0;
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
      double jac = J.K[0][0];

      double xtilde = x[0]*x2[1][0] + (1-x[0])*x2[0][0];

      double funval = fun.eval(xtilde);

      if (m == 0){
        out += w*phis[0]*funval*J.det;
      } else if (m == 1){
        out += w*phis[1]*funval*J.det;
      }

    }
    return out;
  }

  void send_contribution(int i, double val, int rank, int tag){
    // std::cout << "rank " << rank_ << " tries to send " << i << " to rank " << rank << std::endl;
    // MPI_Request req;
    MPI_Send(&i,1,MPI_INT,rank,tag+n_*0,MPI_COMM_WORLD);
    MPI_Send(&val,1,MPI_DOUBLE,rank,tag+n_*1,MPI_COMM_WORLD);
  }
  void recieve_contribution(int i, int rank, int tag, double * val, int * ind){
    // std::cout << "rank " << rank_ << " tries to recieve " << i+n_*tag << ", " << " from rank " << rank <<std::endl;
    // MPI_Request req;
    MPI_Recv(ind,1,MPI_INT,rank,tag+n_*0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(val,1,MPI_DOUBLE,rank,tag+n_*1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }

  void get_nrows_and_offset(int * nrowsin, int *offsetin){
    int offset = MAX;
    int nverts_owned = 0;
    for(int i = 0; i < nverts_; i++){
      if (owner_local_vert(i) == rank_){
        nverts_owned++;
        // std::cout << rank_ << ": " << i << " is mapped to " << local_to_global_vert(i) << std::endl;
        if(local_to_global_vert(i) < offset){
          offset = local_to_global_vert(i);
        }
      }
    }
    *nrowsin = nverts_owned;
    *offsetin = offset;
  }
  Vector initialize_1d_vector(){
    int offset,nverts_owned;
    get_nrows_and_offset(&nverts_owned,&offset);
    Vector b = Vector(nverts_owned,offset);
    b.set_shared_with(rank_,nprocs_);
    return b;
  }
  Matrix initialize_1d_matrix(){
    int offset,nverts_owned;
    get_nrows_and_offset(&nverts_owned,&offset);
    // std::cout << "rank " << rank_ << " has offset " << offset << std::endl;
    // std::cout << "rank " << rank_ << " has DoF " << nverts_owned << std::endl;
    Matrix K = Matrix(nverts_owned,offset,rank_,nprocs_);
    K.initialize_1d(n_,nverts_owned);
    return K;
  }
  Matrix assemble_stiffness_matrix_parallel(){
    Matrix K = initialize_1d_matrix();
    for(int cell = 0; cell < ncells_; cell++){
      int k = *M.topo()(cell);
      int l = *(M.topo()(cell)+1);

      int m = local_to_DoF[k];
      int n = local_to_DoF[l];

      int m_g = cell_to_vert_global(cell,0);
      int n_g = cell_to_vert_global(cell,1);

      if(is_owner_local_vert(k) && is_owner_local_vert(l)){
        // std::cout << "rank "<< rank_ << " owns " << k << ", " << l << std::endl;
        // std::cout << "with local DoF incicies " << m << ", " << n << std::endl;
        K.Diagonal.insert_add(m,n,get_stiffness_element(cell,0,1));
        K.Diagonal.insert_add(n,m,get_stiffness_element(cell,0,1));
        K.Diagonal.insert_add(m,m,get_stiffness_element(cell,0,0));
        K.Diagonal.insert_add(n,n,get_stiffness_element(cell,1,1));
      } else if (is_owner_local_vert(k)){
        K.Diagonal.insert_add(m,m,get_stiffness_element(cell,0,0));
        K.insert_add_from_global(m_g,n_g,get_stiffness_element(cell,0,1));
        send_contribution(m_g,get_stiffness_element(cell,1,0),owner_local_vert(l),0);
        send_contribution(m_g,get_stiffness_element(cell,1,1),owner_local_vert(l),1);
      } else if (is_owner_local_vert(l)){
        // std::cout << "rank " << rank_ << " is haunted ... " << m_g << ", ";
        // std::cout << n_g << std::endl;
        // std::cout << "Tries to send " << n_g << " to rank " << owner_local_vert(k) << std::endl;
        K.Diagonal.insert_add(n,n,get_stiffness_element(cell,1,1));
        // std::cout << get_stiffness_element(cell,1,0) << " is about to be inserted at " << n_g << ", " << m_g << std::endl;
        K.insert_add_from_global(n_g,m_g,get_stiffness_element(cell,1,0));
        send_contribution(n_g,get_stiffness_element(cell,0,1), owner_local_vert(k),0);
        send_contribution(n_g,get_stiffness_element(cell,0,0), owner_local_vert(k),1);
      }

      if (is_shared_local(k) && is_owner_local_vert(k)){
        int r = shared_with(k);
        double val; int ind;
        recieve_contribution(m_g,r,0,&val,&ind);
        K.insert_add_from_global(m_g,ind,val);
        recieve_contribution(m_g,r,1,&val,&ind);
        K.insert_add_from_global(m_g,m_g,val);
      }
      if (is_shared_local(l) && is_owner_local_vert(l)){
        int r = shared_with(l);
        double val; int ind;
        recieve_contribution(n_g,r,0,&val,&ind);
        K.insert_add_from_global(n_g,ind,val);
        recieve_contribution(n_g,r,1,&val,&ind);
        K.insert_add_from_global(n_g,n_g,val);
      }
    }
    return K;
  }
  Matrix assemble_mass_matrix_parallel(){
    Matrix K = initialize_1d_matrix();
    for(int cell = 0; cell < ncells_; cell++){
      int k = *M.topo()(cell);
      int l = *(M.topo()(cell)+1);

      int m = local_to_DoF[k];
      int n = local_to_DoF[l];

      int m_g = cell_to_vert_global(cell,0);
      int n_g = cell_to_vert_global(cell,1);

      if(is_owner_local_vert(k) && is_owner_local_vert(l)){
        // std::cout << "rank "<< rank_ << " owns " << k << ", " << l << std::endl;
        // std::cout << "with local DoF incicies " << m << ", " << n << std::endl;
        K.Diagonal.insert_add(m,n,get_mass_element(cell,0,1));
        K.Diagonal.insert_add(n,m,get_mass_element(cell,0,1));
        K.Diagonal.insert_add(m,m,get_mass_element(cell,0,0));
        K.Diagonal.insert_add(n,n,get_mass_element(cell,1,1));
      } else if (is_owner_local_vert(k)){
        K.Diagonal.insert_add(m,m,get_mass_element(cell,0,0));
        K.insert_add_from_global(m_g,n_g,get_mass_element(cell,0,1));
        send_contribution(m_g,get_mass_element(cell,1,0),owner_local_vert(l),0);
        send_contribution(m_g,get_mass_element(cell,1,1),owner_local_vert(l),1);
      } else if (is_owner_local_vert(l)){
        // std::cout << "rank " << rank_ << " is haunted ... " << m_g << ", ";
        // std::cout << n_g << std::endl;
        // std::cout << "Tries to send " << n_g << " to rank " << owner_local_vert(k) << std::endl;
        K.Diagonal.insert_add(n,n,get_mass_element(cell,1,1));
        // std::cout << get_stiffness_element(cell,1,0) << " is about to be inserted at " << n_g << ", " << m_g << std::endl;
        K.insert_add_from_global(n_g,m_g,get_mass_element(cell,1,0));
        send_contribution(n_g,get_mass_element(cell,0,1), owner_local_vert(k),0);
        send_contribution(n_g,get_mass_element(cell,0,0), owner_local_vert(k),1);
      }

      if (is_shared_local(k) && is_owner_local_vert(k)){
        int r = shared_with(k);
        double val; int ind;
        recieve_contribution(m_g,r,0,&val,&ind);
        K.insert_add_from_global(m_g,ind,val);
        recieve_contribution(m_g,r,1,&val,&ind);
        K.insert_add_from_global(m_g,m_g,val);
      }
      if (is_shared_local(l) && is_owner_local_vert(l)){
        int r = shared_with(l);
        double val; int ind;
        recieve_contribution(n_g,r,0,&val,&ind);
        K.insert_add_from_global(n_g,ind,val);
        recieve_contribution(n_g,r,1,&val,&ind);
        K.insert_add_from_global(n_g,n_g,val);
      }
    }
    return K;
  }

  Vector assemble_rhs_vector_parallel(Function fun){
    Vector Vec = initialize_1d_vector();
    for(int cell = 0; cell < ncells_; cell++){
      int k = *M.topo()(cell);
      int l = *(M.topo()(cell)+1);

      int m = local_to_DoF[k];
      int n = local_to_DoF[l];

      int m_g = cell_to_vert_global(cell,0);
      int n_g = cell_to_vert_global(cell,1);

      if(is_owner_local_vert(k) && is_owner_local_vert(l)){
        // std::cout << "rank "<< rank_ << " owns " << k << ", " << l << std::endl;
        // std::cout << "with local DoF incicies " << m << ", " << n << std::endl;
        Vec.insert_add_from_global(m_g,get_rhs_element(cell,0,fun));
        Vec.insert_add_from_global(n_g,get_rhs_element(cell,1,fun));
      } else if (is_owner_local_vert(k)){
        Vec.insert_add_from_global(m_g,get_rhs_element(cell,0,fun));
        send_contribution(m_g,get_rhs_element(cell,1,fun),owner_local_vert(l),0);
      } else if (is_owner_local_vert(l)){
        // std::cout << "rank " << rank_ << " is haunted ... " << m_g << ", ";
        // std::cout << n_g << std::endl;
        // std::cout << "Tries to send " << n_g << " to rank " << owner_local_vert(k) << std::endl;
        Vec.insert_add_from_global(n_g,get_rhs_element(cell,1,fun));
        // std::cout << get_stiffness_element(cell,1,0) << " is about to be inserted at " << n_g << ", " << m_g << std::endl;
        send_contribution(n_g,get_rhs_element(cell,0,fun), owner_local_vert(k),1);
      }

      if (is_shared_local(k) && is_owner_local_vert(k)){
        int r = shared_with(k);
        double val; int ind;
        recieve_contribution(m_g,r,0,&val,&ind);
        Vec.insert_add_from_global(m_g,val);
      }
      if (is_shared_local(l) && is_owner_local_vert(l)){
        int r = shared_with(l);
        double val; int ind;
        recieve_contribution(n_g,r,1,&val,&ind);
        Vec.insert_add_from_global(n_g,val);
      }
    }
    return Vec;
  }
  Vector interpolate_parallel(Function fun){
    Vector vec = initialize_1d_vector();
    for(int i = 0; i < nverts_; i++){
      // std::cout << "At rank " << rank_ << " : " << i << " is owned by " << owner_local_vert(i) << std::endl;
      if (is_owner_local_vert(i)){
        double val = fun.eval(*M.geom()(i));
        int i_g = local_to_global_vert(i);
        vec.insert_from_global(i_g,val);
        // std::cout << "insert " << val << " at " << i_g << std::endl;
      }
    }
    return vec;

  }
  void interpolate(Function fun, double * vec){
    for(int i = 0; i < nverts_; i++){
      vec[i] = fun.eval(*M.geom()(i));
    }
  }
};

void CG(Matrix A, Vector x, Vector b, Vector r, Vector p, Vector Ap){
  A.mult(x,Ap);
  r.copy(b);
  r.add_vec(-1,Ap);
  p.copy(r);
  // std::cout << "iter " << -1 << ", rnorm " << r.norm() << std::endl;

  int maxiter = 20000;
  double tol = 1e-8;

  double alpha,beta,rinner;
  for(int k = 0; k < maxiter; k++){
    A.mult(p,Ap);
    rinner = r.inner_prod(r);
    alpha = rinner/p.inner_prod(Ap);
    x.add_vec(alpha,p);
    r.add_vec(-alpha,Ap);
    // std::cout << "iter " << k << ", rnorm " << r.norm() << std::endl;
    if (r.norm() < tol){return;}
    beta = r.inner_prod(r)/rinner;
    p.copy(r);
    p.add_vec(beta,p);
  }
}

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
  Matrix K = UM.assemble_stiffness_matrix_parallel();
  Matrix M = UM.assemble_mass_matrix_parallel();
  Vector f = UM.initialize_1d_vector();
  Function fun = Function();
  Vector fint = UM.interpolate_parallel(fun);
  Vector Mfint = UM.initialize_1d_vector();
  M.mult(fint,Mfint);

  // Vector res = UM.initialize_1d_vector();
  Vector Ku = UM.initialize_1d_vector();

  // M.mult(fint,f);
  Vector u = UM.initialize_1d_vector();
  Vector r = UM.initialize_1d_vector();
  Vector p = UM.initialize_1d_vector();
  Vector Ap = UM.initialize_1d_vector();

  f = UM.assemble_rhs_vector_parallel(fun);
  // M.mult(fint,f);
  f.mult_with_const(4*M_PI*M_PI);

  u.set_to_const(0);
  // u.copy(fint);

  Vector res = UM.initialize_1d_vector();
  CG(K,u,f,r,p,Ap);
  K.mult(u,Ku);
  res.copy(Ku);
  res.add_vec(-1,f);

  Vector err = UM.initialize_1d_vector();
  Vector Mu = UM.initialize_1d_vector();
  err.copy(u);
  err.add_vec(-1,fint);
  double resnorm = res.norm();
  double errnorm = err.max_norm();
  Mfint.mult_with_const(4*M_PI*M_PI);

  if (rank == 0){
    K.dump();
    // M.dump();
    // M.OffDiagonal.dump();
    std::cout << "u: " << std::endl;
    u.dump();
    std::cout << "fint: " << std::endl;
    fint.dump();
    // std::cout << "Mfint: " << std::endl;
    // Mfint.dump();
    // std::cout << "f: " << std::endl;
    // f.dump();
    // std::cout << "res: " << std::endl;
    // res.dump();
    std::cout << "norm of res " << resnorm << std::endl;
    std::cout << "norm of err " << errnorm << std::endl;

  }
  // if (rank == 0){
  //   CRS K = CRS(3,4);
  //   K.insert(0,0,1);
  //   K.insert(1,1,1);
  //   //K.insert(2,1,1);
  //   K.insert(2,2,1);
  //   K.dump();
  //   K.dump_full();
  //
  //   double * x = new double[3];
  //   double * y = new double[3];
  //   x[0] = 1;
  //   x[1] = 2;
  //   x[2] = 3;
  //
  //   K.mult(x,y);
  //   printvec(y,3);
  // }
  // if (rank == 0){
  //   // FullSquareMatrix K = FullSquareMatrix(n+1);
  //   CRS K = CRS(n+1,3*(n-1) + 4);
  //   // K.dump();
  //   UM.assemble_stiffness_matrix(K);
  //   K.dump();
  //   K.dump_full();
  //   // std::cout << "Stiffness Matrix:" << std::endl;
  //   // K.dump();
  //
  //   // FullSquareMatrix M = FullSquareMatrix(n+1);
  //   CRS M = CRS(n+1,3*(n-1) + 4);
  //   UM.assemble_mass_matrix(M);
  //   M.dump_full();
  //   // std::cout << "Mass Matrix:" << std::endl;
  //   // M.dump();
  //
  //   Function cosfun;
  //   double * vec = new double[n+1];
  //   UM.interpolate(cosfun,vec);
  //   // std::cout << "vec: " << std::endl;
  //   // printvec(vec,n+1);
  //   double * Kvec = new double[n+1];
  //   K.mult(vec,Kvec);
  //   // std::cout << "\n Kvec : " << std::endl;
  //   // printvec(Kvec,n+1);
  //
  //   double * Mvec = new double[n+1];
  //   M.mult(vec,Mvec);
  //   // std::cout << "\n Mvec: " << std::endl;
  //   // printvec(Mvec,n+1);
  //
  //   double * r = new double[n+1];
  //   double rnorm = 0;
  //   for(int i = 0; i < n+1; i++){
  //     r[i] = Kvec[i] - 4*M_PI*M_PI*Mvec[i];
  //     rnorm += r[i]*r[i];
  //   }

    // std::cout << "\n\n residual: " << std::endl;
    // printvec(r,n+1);

    // std::cout << "\n\n .... rnorm = " << sqrt(rnorm) << std::endl;


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
