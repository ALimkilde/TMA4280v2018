void NotImpl(){
  throw std::logic_error{"Function not yet implemented."};
  return;
}

class Mesh{
public:
  int* topology;
  double* geometry;
  int n_cells,n_verts;
};

class MeshPartition: public Mesh{
public:
  int shared[];
  int ghosts[];
  int rank,nprocs,size,offset,n_verts_local;
  MeshPartition(int rankin, int nprocsin, int Nin){
    rank = rankin;
    nprocs = nprocsin;
    n_cells = Nin;            // Only holds for 1d... FIXIT
    n_verts =  n_cells+1;      // Only holds for 1d... FIXIT
    size = n_cells/nprocs;
    n_verts_local = size + 1;  // Only holds for 1d... FIXIT
    offset = rank*size;
    generate_geometry();
    generate_topology();

    return;
  }
  void print_range(){
    std::cout << "rank " << rank << " has range (" << size << "," << offset << ")" << std::endl;
  }
  void print_geometry(){
    std::cout << "Geometry at rank " << rank << " is: [";
    for(int i = 0; i < n_verts_local; i++){
      std::cout << geometry[i] << ", ";
    }
    std::cout << "]" << std::endl;
    return;
  }
  int get_owner(int i){
    return i/(n_verts_local);
    return 0;
  }
  int local_to_global(int i){
    return offset + i;
    return 0;
  }
  void generate_topology(){  // Only holds for 1d... FIXIT
    topology = new int[size*2];

    for(int i = 0; i < size; i++){
     topology[2*i] = i;
     topology[2*i+1] = i+1;
    }
    //return;
  }
  void generate_geometry(){ // Only holds for 1d... FIXIT
    geometry = new double[size + 1];

    for(int i = 0; i < n_verts_local; i++){
      geometry[i] = double(offset + i)/n_cells;
    }
  }
};
