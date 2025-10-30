#ifndef LEVEL_HEADER
#define LEVEL_HEADER

class Level {

public:
  Level() = default; // default constructor
  ~Level();          // default destructor

  // NOTE: we = delete anything unused until actually needed, then it gets implemented
  // copy constructor and assigment operator
  Level(const Level &) = delete;
  Level &operator=(const Level &) = delete;

  // move constructor and move assigment operator
  Level(const Level &&) = delete;
  Level &operator=(const Level &&) = delete;

  void allocateMemory();
  void initialize();

  void defineInterpolation(vector_double *V);
  void setupCoarseGridCorrection();
  void coarseOperatorAllocateMemory();
  void setupCoarseOperator(vector_double *V);

  // distributed: non-idling processes of previos level
  // gathered: non-idling processes of current level
  operator_struct<double> op_double;    // distributed
  operator_struct<double> oe_op_double; // odd_even
  Schwarz s_double;
  interpolation_double_struct is_double; // interpolation / aggregation
  gathering_double_struct gs_double;     // gathering parameters and buffers
  gmres_double_struct dummy_p_double;    // dummy gmres struct
  profiling_double_struct prof_double;   // profiling

  Gmres gmres;
  Gmres gmresSmoother;
  Gmres gmresCoarseSolve;
  Gmres gmresKCycle;

  // communication
  MPI_Request *reqs;
  MPI_Comm comm;
  int parent_rank, idle = 0, neighbor_rank[8], num_processes, num_processes_dir[4];
  // lattice
  int *global_lattice;
  int *local_lattice;
  int *block_lattice;
  int num_eig_vect;
  int num_parent_eig_vect;
  int coarsening[4];
  int global_splitting[4];
  int periodic_bc[4];
  int comm_offset[4];
  // degrees of freedom on a site
  // 12 on fine lattice (i.e., complex d.o.f.)
  // 2*num_eig_vect on coarser lattices
  int num_lattice_site_var;
  int currentLevel;
  int depth = 0;
  // number of sites in local volume + ghost shell (either fw or bw)
  int num_lattice_sites;
  // number of sites in local volume
  int num_inner_lattice_sites;
  int num_boundary_sites[4];
  // complex d.o.f. in local volume + ghost shell = num_lattice_sites * num_lattice_site_var
  long int vector_size;
  // complex d.o.f. in local volume = num_inner_lattice_sites * num_lattice_site_var
  long int inner_vector_size;
  long int schwarz_vector_size;
  int D_size;
  int clover_size;
  int block_size;
  // buffer vectors
  vector_double vbuf_double[9], sbuf_double[2];
  // storage + daggered-operator bufferes
  vector_double x;
  // local solver parameters
  double tol, relax_fac;
  int n_cy, post_smooth_iter, block_iter, setup_iter;

#if defined(GCRODR) || defined(POLYPREC)
  // 'bool', if on H will be copied
  int dup_H;
#endif

  // neighboring levels
  Level *next_level = nullptr;
  Level *previous_level = nullptr;

  int type = _FINE;
  int odd_even;
  double csw;
  Global *global;
};

#endif
