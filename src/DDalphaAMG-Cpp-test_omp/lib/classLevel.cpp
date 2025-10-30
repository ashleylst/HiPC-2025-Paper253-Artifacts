
#include "main.h"

void Level::allocateMemory() {

  int onePlusOneFactor = 1;
#ifdef HAVE_TM1p1
  onePlusOneFactor *= 2;
#endif

  for (int i = 0; i < 9; i++)
    vbuf_double[i] = new complex_double[onePlusOneFactor * vector_size];

  for (int i = 0; i < 2; i++)
    sbuf_double[i] = new complex_double[onePlusOneFactor * schwarz_vector_size];
}

Level::~Level() {

  for (int i = 0; i < 9; i++)
    delete[] vbuf_double[i];
  for (int i = 0; i < 2; i++)
    delete[] sbuf_double[i];

  if (type != _FINE) {
    operator_double_free(&(op_double), _ORDINARY, this);
    //std::cout << &(op_double.D) << std::endl;
    gathering_double_free(&gs_double, this);
  }
  if (type != _COARSEST) {
    coarsening_index_table_double_free(&(is_double), this);
    interpolation_double_free(this);
  }
}

void Level::initialize() {

  // initialize class members
  operator_double_init(&(op_double));
  operator_double_init(&(oe_op_double));
  interpolation_double_struct_init(&(is_double));

  // copy data from global class NOTE: Will hopefully be removed at some point
  currentLevel = global->num_levels - depth;
  global_lattice = global->global_lattice[depth];
  local_lattice = global->local_lattice[depth];
  block_lattice = global->block_lattice[depth];

  post_smooth_iter = global->post_smooth_iter[depth];
  n_cy = global->ncycle[depth];
  relax_fac = global->relax_fac[depth];
  block_iter = global->block_iter[depth];
  setup_iter = global->setup_iter[depth];

  tol = global->tol;
  odd_even = global->odd_even;
  csw = global->csw;

  num_eig_vect = global->num_eig_vect[depth];
  num_parent_eig_vect = (type == _FINE) ? 6 : previous_level->num_eig_vect;
  num_lattice_site_var = 2 * num_parent_eig_vect;

  for (int i = 0; i < 8; i++) {
    neighbor_rank[i] = global->neighbor_rank[i];
  }
  num_processes = 1;
  for (int mu = 0; mu < 4; mu++) {
    num_processes_dir[mu] =
        global_lattice[mu] / local_lattice[mu]; // NOTE: same variable as global splitting?
    num_processes *= num_processes_dir[mu];
    comm_offset[mu] = (type == _FINE)
                          ? 1
                          : (previous_level->num_processes_dir[mu] / num_processes_dir[mu]) *
                                previous_level->comm_offset[mu];
    coarsening[mu] = (type == _COARSEST || global->method == 0)
                         ? local_lattice[mu]
                         : global_lattice[mu] / global->global_lattice[depth + 1][mu];
    global_splitting[mu] = global_lattice[mu] / local_lattice[mu];

    periodic_bc[mu] = 1;
  }

  // compute vector sizes (former data_layout_init)
  num_inner_lattice_sites = 1;
  for (int i = 0; i < 4; i++)
    num_inner_lattice_sites *= local_lattice[i];
  num_lattice_sites = num_inner_lattice_sites;
  inner_vector_size = num_inner_lattice_sites * num_lattice_site_var;

  int j = num_lattice_sites; // j = num_inner_lattice_sites?! NOTE: double check, then get rid of j
  for (int i = 0; i < 4; i++)
    num_lattice_sites += j / local_lattice[i];

  vector_size = num_lattice_sites * num_lattice_site_var;
  schwarz_vector_size = 2 * vector_size - inner_vector_size;

  // sanity checks //TODO move to more sane location together with all other asserts, probably right
  // after reading input file
  for (int mu = 0; mu < 4; mu++)
    ASSERT(IMPLIES(global->rhs == 3,
                   ASCENDING(0, global->propagator_coords[mu], global_lattice[mu] - 1)));
  ASSERT(IMPLIES(global->method > 0, n_cy > 0));

  Level::allocateMemory();

  if (type != _FINE)
    //     cart_define(MPI_COMM_WORLD, this);
    //   else
    neighbor_define(this);

  // setup gathering struct (also determines idle processes on coarse levels)
  if (type != _FINE) {
    gathering_double_init(&gs_double, this);
    gathering_double_setup_new(&gs_double, this);
  }

  // setup Dirac operator on each level
  if (!idle) {
    if (type == _FINE) {
      op_double = global->op_double; // NOTE: this is not a real copy! Just matching pointers
                                     // (this->op_double not properly alloc'd)
    } else {
      Level::coarseOperatorAllocateMemory();
      Level::setupCoarseOperator(previous_level->is_double.interpolation);
    }

    // initialize fine solver
    if (type == _FINE) {
      switch (global->method) {
      case 0:
        gmres.initialize(&global->op_double, global->restart, global->max_restart,
                         inner_vector_size, tol, _GLOBAL_FGMRES, _NOTHING, 0, 1, vector_size,
                         nullptr, d_plus_clover_double, this);
        break; // pure GMRES on D
      case 1:  // methods 1-4 all initialize a right preconditioned FMGRES
      case 2:  // on D with different smoothers, so cases 1-3 can fall-through
      case 3:  // to case 4
      case 4:
        gmres.initialize(&global->op_double, global->restart, global->max_restart,
                         inner_vector_size, tol, _GLOBAL_FGMRES, _RIGHT, 0, 1, vector_size,
                         preconditioner, d_plus_clover_double, this);
        break;
      case 6:
        gmres.initialize(&global->op_double, global->restart, global->max_restart,
                         inner_vector_size, tol, _GLOBAL_FGMRES, _RIGHT, 0, 1, vector_size,
                         preconditioner, g5D_plus_clover_double, this);
        break; // GMRES + MG on Q (with GMRES smoothing)
      default:
        error0("Invalid outer solver method (global->method = %d) \n", global->method);
        break;
      }
    }

    // initialize smoother
    void (*matVec)(vector_double, vector_double, operator_struct<double> *, Level *) =
        d_plus_clover_double;
    s_double.initialize(matVec, this); // Initializes Schwarz smoother on all levels
    if (type != _COARSEST) {
      switch (global->method) {
      case 0:
        break; // no MG -> no smoother
      case 1:
        break; // Schwarz already handled above
      case 2:
        break; // Schwarz already handled above
      case 3:
        break; // Schwarz already handled above
      case 4:
        if (type == _FINE)
          matVec = (odd_even) ? apply_schur_complement_double : d_plus_clover_double;
        else
          matVec = (odd_even) ? coarse_apply_schur_complement_double : apply_coarse_operator_double;
        gmresSmoother.initialize((odd_even) ? &oe_op_double : &s_double.op, block_iter, 2,
                                 inner_vector_size, tol, _COARSE_GMRES, _NOTHING, 0, 0, vector_size,
                                 nullptr, matVec, this);
        break;
      case 6:
        if (type == _FINE)
          matVec = (odd_even) ? g5D_apply_schur_complement_double : g5D_plus_clover_double;
        else
          matVec = (odd_even) ? g5D_coarse_apply_schur_complement_double
                              : g5D_apply_coarse_operator_double;
        gmresSmoother.initialize((odd_even) ? &oe_op_double : &s_double.op, block_iter, 2,
                                 inner_vector_size, tol, _COARSE_GMRES, _NOTHING, 0, 0, vector_size,
                                 nullptr, matVec, this);
        break;
      default:
        error0("Invalid smoother (global->method = %d) \n", global->method);
        break;
      }
    }

    // setup odd even operator
    if (global->method >= 4 && odd_even) {
      switch (type) {
      case _FINE:
        oddeven_setup_double(&(op_double), this);
        break;
      case _INTERMEDIATE: {
        coarse_oddeven_alloc_double(this);
        coarse_oddeven_setup_double(&(s_double.op), _REORDER, this);
        break;
      }
      case _COARSEST: {
        coarse_oddeven_alloc_double(this);
        coarse_oddeven_setup_double(&(s_double.op), _NO_REORDERING, this);
        break;
      }
      default:
        error0("Invalid level type for odd even setup\n");
        break;
      }
    }

    if (type != _COARSEST) {
      coarsening_index_table_double_alloc(&(is_double), this);
      coarsening_index_table_double_define(&(is_double), &(s_double), this);
    }

    // setup k cycle gmres on intermediate levels
    if (type == _INTERMEDIATE && global->kcycle) {
      gmresKCycle.initialize(
          &s_double.op, global->kcycle_restart, global->kcycle_max_restart, vector_size,
          global->kcycle_tol, _K_CYCLE, _RIGHT, 0, 0, vector_size, vcycle_double,
          global->method == 6 ? g5D_apply_coarse_operator_double : apply_coarse_operator_double,
          this);
    }

    // initialize coarsest level solver TODO: init other solvers
    if (type == _COARSEST) {
      if (odd_even)
        gmresCoarseSolve.initialize(&oe_op_double, global->coarse_iter, global->coarse_restart,
                                    inner_vector_size, global->coarse_tol, _COARSE_GMRES, _NOTHING,
                                    0, 0, vector_size, nullptr,
                                    global->method == 6 ? g5D_coarse_apply_schur_complement_double
                                                        : coarse_apply_schur_complement_double,
                                    this);
      else
        gmresCoarseSolve.initialize(
            &op_double, global->coarse_iter, global->coarse_restart, inner_vector_size,
            global->coarse_tol, _COARSE_GMRES, _NOTHING, 0, 0, vector_size, nullptr,
            global->method == 6 ? g5D_apply_coarse_operator_double : apply_coarse_operator_double,
            this);
    }

    // obtain test vectors and interpolation operator
    if (type != _COARSEST) {
      interpolation_double_alloc(this);
      defineInterpolation(nullptr); // no initial testvectors available on fine level
    }
    // setup coarse grid correction
    if (type == _COARSEST)
      setupCoarseGridCorrection();
  }
}

void Level::defineInterpolation(vector_double *V) {

  double norm;

  if (V == nullptr) {
    if (global->print > 0)
      printf0("initial definition --- depth: %d\n", depth);

    for (int i = 0; i < num_eig_vect; i++) {
      //       if ( type == _FINE ) {
      //
      //             vector_double_define_random(vbuf_double[0], 0, inner_vector_size, this);
      // TODO: This ensures a deterministic outcome after setup, for better performance use above
      // line (random initial guesses)
      for (int kk = 0; kk < inner_vector_size; kk++)
        vbuf_double[0][kk] = kk % (i + 2) == 0 ? i : 1;
      //
      //       }

      smoother_double(is_double.test_vector[i], nullptr, vbuf_double[0], 1, _NO_RES, this);
      //       smoother_double(vbuf_double[0], nullptr, is_double.test_vector[i], global->method >=
      //       4 ? 1 : 2,
      //                       _NO_RES, this);
      //       smoother_double(is_double.test_vector[i], nullptr, vbuf_double[0], global->method >=
      //       4 ? 1 : 3,
      //                       _NO_RES, this);
    }

    for (int i = 0; i < num_eig_vect; i++) {
      norm = global_norm_double(is_double.test_vector[i], 0, inner_vector_size, this);
      printf0("norm TV[%d] = %.15lf\n", i, norm);
      vector_double_real_scale(is_double.test_vector[i], is_double.test_vector[i], 1.0 / norm, 0,
                               inner_vector_size, this);
    }
  } else {
    for (int i = 0; i < num_eig_vect; i++) {
      trans_double(is_double.test_vector[i], V[i], s_double.op.translation_table, this);
    }
  }

  //   error0("");

  for (int i = 0; i < num_eig_vect; i++) {
    vector_double_copy(is_double.interpolation[i], is_double.test_vector[i], 0, inner_vector_size,
                       this);
  }
  gram_schmidt_on_aggregates_double(is_double.interpolation, num_eig_vect, this);
  define_interpolation_double_operator(is_double.interpolation, this);
}

void Level::setupCoarseGridCorrection() {

  if (!idle) {
    // free memory in case of re-allocation
    operator_double_free(&(s_double.op), _ORDINARY, this);

    operator_double_init(&(s_double.op));
    operator_double_alloc(&(s_double.op), _ORDINARY, this);
    operator_double_define(&(s_double.op), this);
    conf_double_gather(&(s_double.op), &(op_double), this);

    if (type != _COARSEST) {
      schwarz_double_boundary_update(&(s_double), this);
      if (global->method >= 4 && odd_even) {
        coarse_oddeven_alloc_double(this);
        coarse_oddeven_setup_double(&(s_double.op), _REORDER, this);
      }
      coarse_operator_double_set_couplings(&(s_double.op), this);
    }
    if (type == _COARSEST) {
      if (odd_even) {
        coarse_oddeven_alloc_double(this);
        coarse_oddeven_setup_double(&(s_double.op), _NO_REORDERING, this);
      } else {
        coarse_operator_double_set_couplings(&(s_double.op), this);
      }
    }
  }
}

void Level::coarseOperatorAllocateMemory() {

  int nd = num_inner_lattice_sites, k = num_parent_eig_vect * 2;
  D_size = k * k * 4 * nd;
  clover_size = ((k * (k + 1)) / 2) * nd;
  block_size = ((k / 2 * (k / 2 + 1))) * nd;

  operator_double_alloc(&(op_double), _ORDINARY, this);
}

void Level::setupCoarseOperator(vector_double *V) {

  int n = previous_level->num_eig_vect;
  double t = -MPI_Wtime();
  vector_double buffer1 = previous_level->vbuf_double[4], buffer2 = previous_level->vbuf_double[5];

  void (*aggregate_self_coupling)(vector_double, vector_double, vector_double, Schwarz *, Level *) =
      (previous_level->type == _FINE) ? d_plus_clover_aggregate_double
                                      : coarse_aggregate_self_couplings_double;
  void (*aggregate_neighbor_coupling)(vector_double, vector_double, vector_double, const int,
                                      Schwarz *, Level *) =
      (previous_level->type == _FINE) ? d_neighbor_aggregate_double
                                      : coarse_aggregate_neighbor_couplings_double;
  void (*aggregate_block)(vector_double, vector_double, vector_double, config_double, Level *) =
      (previous_level->type == _FINE) ? diagonal_aggregate_double
                                      : coarse_aggregate_block_diagonal_double;

  operator_double_define(&(op_double), this);

  for (int i = 0; i < D_size; i++)
    op_double.D[i] = 0;
  for (int i = 0; i < clover_size; i++)
    op_double.clover[i] = 0;
  for (int i = 0; i < block_size; i++)
    op_double.odd_proj[i] = 0;

  // for all test vectors V[i]:
  for (int i = 0; i < n; i++) {
    for (int mu = 0; mu < 4; mu++) {
      // update ghost cells of V[i]
      negative_sendrecv_double(V[i], mu, &(previous_level->s_double.op.c), previous_level);
    }
    // apply self coupling of block-and-2spin-restricted dirac operator for each aggregate
    aggregate_self_coupling(buffer1, buffer2, V[i], &(previous_level->s_double), previous_level);
    // calculate selfcoupling entries of the coarse grid operator
    set_coarse_self_coupling_double(buffer1, buffer2, V, i, previous_level);
    // odd_proj
    aggregate_block(buffer1, buffer2, V[i], previous_level->s_double.op.odd_proj, previous_level);
    set_block_diagonal_double(buffer1, buffer2, V, i, op_double.odd_proj, previous_level);
    for (int mu = 0; mu < 4; mu++) {
      // finish updating ghostcells of V[i]
      negative_wait_double(mu, &(previous_level->s_double.op.c), previous_level);
      // apply 2spin-restricted dirac operator for direction mu for all aggregates
      aggregate_neighbor_coupling(buffer1, buffer2, V[i], mu, &(previous_level->s_double),
                                  previous_level);
      set_coarse_neighbor_coupling_double(buffer1, buffer2, V, mu, i, previous_level);
    }
  }

  coarse_operator_double_setup_finalize(previous_level);

  t += MPI_Wtime();
  if (global->print > 0)
    printf0("depth: %d, time spent for setting up next coarser operator: %lf seconds\n", depth, t);
}
