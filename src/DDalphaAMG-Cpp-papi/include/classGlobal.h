#ifndef GLOBAL_HEADER
#define GLOBAL_HEADER

class Global {

public:
  Global() = default;
  ~Global();

  // global class should be unique, so no copy constructor or assignment operator allowed here
  Global(const Global &) = delete;
  Global &operator=(const Global &) = delete;

  // similarly no moving is required, so no move constructor or move assignment operator
  Global(const Global &&) = delete;
  Global &operator=(const Global &&) = delete;

  // post-instantiation data allocation
  void allocateDataAfterReadingInputFile();

  // read input files
  void readGlobalInfo();
  void readNoDefaultInfo();
  void readSolverInfo();
  void readGeometryInfo();
  void readTestvectorIOData();
  void readEvaluationParameters();
  void readKcycleData();
  void validateParameters();

  // setup of global/fine operator
  void initialize();
  void setupFineDiracOperator();

  // for debugging
  void printParams();

  FILE *inputfile, *logfile;
  void *save_pt;

  Gmres gmres;

  //   gmres_double_struct p;
  //   gmres_MP_struct p_MP;
  operator_struct<double> op_double;

  // communication
  MPI_Comm comm_cart;
  MPI_Group global_comm_group;
  MPI_Request sreqs[8], rreqs[8];
  int num_processes, my_rank, my_coords[4], tv_io_single_file, num_openmp_processes;
  // string buffers
  char in[STRINGLENGTH], in_clov[STRINGLENGTH], source_list[STRINGLENGTH],
      tv_io_file_name[STRINGLENGTH];
  // geometry, method parameters
  int num_levels, num_desired_levels, process_grid[4], in_format, **global_lattice, **local_lattice,
      **block_lattice, *post_smooth_iter, *block_iter, *setup_iter, *ncycle, method, odd_even, rhs,
      propagator_coords[4], interpolation, randomize, *num_eig_vect, num_coarse_eig_vect, kcycle,
      mixed_precision, restart, max_restart, kcycle_restart, kcycle_max_restart, coarse_iter,
      coarse_restart;
  double tol, coarse_tol, kcycle_tol, csw, rho, *relax_fac;
#ifdef GCRODR
  int gcrodr_k;
  int gcrodr_upd_itrs;
#endif

#ifdef POLYPREC
  int polyprec_d;
#endif

#ifdef BLOCK_JACOBI
  int local_polyprec_d;
#endif

  // profiling, analysis, output
  int coarse_iter_count, iter_count, iterator, print, conf_flag, setup_flag, in_setup;
  double coarse_time, prec_time, *output_table[8], cur_storage, max_storage, total_time, plaq_hopp,
      plaq_clov, norm_res, plaq, bicgstab_tol, twisted_bc[4], test;

  double coarsest_time;

  double m0, setup_m0;

#ifdef HAVE_TM
  // twisted mass parameters
  int downprop;
  double mu, setup_mu, mu_odd_shift, mu_even_shift, *mu_factor = nullptr;
#endif

#ifdef HAVE_TM1p1
  int n_flavours = 1;
  double epsbar, epsbar_ig5_odd_shift, epsbar_ig5_even_shift, *epsbar_factor = nullptr;
#endif

  // index functions for external usage
  int (*conf_index_fct)(), (*vector_index_fct)();
  int *odd_even_table;
  int (*Cart_rank)(MPI_Comm comm, const int coords[], int *rank);
  int (*Cart_coords)(MPI_Comm comm, int rank, int maxdims, int coords[]);

  // bc: 0 dirichlet, 1 periodic, 2 anti-periodic
  int bc;

  complex_double gamma[4][16];
  var_table vt;

  // mostly useful, as of now, for PP and GCRODR
  int on_solve;

  int low_level_meas;
  double avg_b1;
  double avg_b2;
  double avg_crst;

#ifdef PERS_COMMS
  // use persistent communications at the coarsest level
  int use_pers_comms1, use_pers_comms2;
  int pers_comms_id1, pers_comms_id2;

  int pers_comms_nrZs, pers_comms_nrZas, pers_comms_nrZxs;

  // for plus hopping (plus doesn't mean +mu here)
  MPI_Request *pers_comms_recvrs_plus[8];
  MPI_Request *pers_comms_sendrs_plus[8];
  // for minus hopping
  MPI_Request *pers_comms_recvrs_minus[8];
  MPI_Request *pers_comms_sendrs_minus[8];
#endif

  // stuff added to fully initialize Global without info of Level
  int periodic_bc[4];
  int num_processes_dir[4];
  int comm_offset[4];
  int neighbor_rank[8];
  int num_lattice_site_var;
  int num_lattice_sites;
  int num_inner_lattice_sites;
  int inner_vector_size;
  int vector_size;
  int schwarz_vector_size;
  vector_double vbuf_double[9];
};

#endif
