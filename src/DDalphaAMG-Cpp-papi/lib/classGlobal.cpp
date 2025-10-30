
#include "main.h"

void Global::allocateDataAfterReadingInputFile() {

//   global_lattice = new int*[num_levels];
//   local_lattice = new int*[num_levels];
//   block_lattice = new int*[num_levels];
//   for( int i=0;i<num_levels; i++) {
//     global_lattice[i] = new int[4];
//     local_lattice[i] = new int[4];
//     block_lattice[i] = new int[4];
//   }
NEW2D(global_lattice, int, 4, num_levels);
NEW2D(local_lattice, int, 4, num_levels);
NEW2D(block_lattice, int, 4, num_levels);

post_smooth_iter = new int[num_levels];
ncycle = new int[num_levels];
relax_fac = new double[num_levels];
#ifdef HAVE_TM
mu_factor = new double[num_levels];
#endif
#ifdef HAVE_TM1p1
epsbar_factor = new double[num_levels];
#endif
block_iter = new int[num_levels];
setup_iter = new int[num_levels];
num_eig_vect = new int[num_levels];
}

Global::~Global() {
//   for( int i=0;i<num_levels; i++) {
//     delete[] global_lattice[i];
//     delete[] local_lattice[i];
//     delete[] block_lattice[i];
//   }
//   delete[] global_lattice;
//   delete[] local_lattice;
//   delete[] block_lattice;

DELETE2D(global_lattice);
DELETE2D(local_lattice);
DELETE2D(block_lattice);

delete[] post_smooth_iter;
delete[] ncycle;
delete[] relax_fac;
#ifdef HAVE_TM
delete[] mu_factor;
#endif
#ifdef HAVE_TM1p1
delete[] epsbar_factor;
#endif
delete[] block_iter;
delete[] setup_iter;
delete[] num_eig_vect;

delete[] odd_even_table;
}

void Global::readGlobalInfo() {

  // Note: There is actually no default set for the three following values
  // Though, when using the code as a library, no configuration paths are required.
  save_pt = &(in);
  in[0] = '\0';
  read_parameter(&save_pt, "configuration:", "%s", 1, inputfile, _NO_DEFAULT_SET);

  save_pt = &(in_format);
  in_format = _STANDARD;
  read_parameter(&save_pt, "format:", "%d", 1, inputfile, _DEFAULT_SET);

  // right hand side
  save_pt = &(rhs);
  rhs = 1;
  read_parameter(&save_pt, "right hand side:", "%d", 1, inputfile, _DEFAULT_SET);
  if (rhs == 4) { // TODO: enum for rhs
    save_pt = &(source_list);
    read_parameter(&save_pt, "source list:", "%s", 1, inputfile, _NO_DEFAULT_SET);
  } else if (rhs == 3) {
    save_pt = propagator_coords;
    read_parameter(&save_pt, "propagator coordinates:", "%d", 4, inputfile, _DEFAULT_SET);
  }

  save_pt = &(num_levels);
  //   num_levels = 2;
  read_parameter(&save_pt, "number of levels:", "%d", 1, inputfile, _DEFAULT_SET);
  num_levels = MAX(1, num_levels);
  num_desired_levels = num_levels;

  save_pt = &(bc);
  bc = _PERIODIC;
  read_parameter(&save_pt, "boundary conditions:", "%d", 1, inputfile, _DEFAULT_SET);

  if (bc == _TWISTED) {
    save_pt = twisted_bc;
    for (int i = 0; i < 4; i++)
      twisted_bc[i] = 0;
    read_parameter(&save_pt, "twisted boundary conditions:", "%d", 4, inputfile, _DEFAULT_SET);
    for (int i = 0; i < 4; i++)
      twisted_bc[i] *= M_PI;
  }
}

void Global::readNoDefaultInfo() {

  // global lattice
  save_pt = global_lattice[0];
  read_parameter(&save_pt, "d0 global lattice:", "%d", 4, inputfile, _NO_DEFAULT_SET);

  // local lattice
  save_pt = local_lattice[0];
  read_parameter(&save_pt, "d0 local lattice:", "%d", 4, inputfile, _NO_DEFAULT_SET);

  // block lattice
  save_pt = block_lattice[0];
  read_parameter(&save_pt, "d0 block lattice:", "%d", 4, inputfile, _NO_DEFAULT_SET);

  // Wilson mass
  save_pt = &(m0);
  m0 = 0;
  read_parameter(&save_pt, "m0:", "%lf", 1, inputfile, _DEFAULT_SET);
  if (m0 == 0) {
    double kappa = 0;
    save_pt = &(kappa);
    read_parameter(&save_pt, "kappa:", "%lf", 1, inputfile, _DEFAULT_SET);
    ASSERT(kappa != 0);
    m0 = 1. / (2. * kappa) - 4.; // setting m0 from kappa
  }
  save_pt = &(csw);
  read_parameter(&save_pt, "csw:", "%lf", 1, inputfile, _NO_DEFAULT_SET);

#ifdef HAVE_TM
  save_pt = &(mu);
  mu = 0;
  read_parameter(&save_pt, "mu:", "%lf", 1, inputfile, _DEFAULT_SET);
  if (mu == 0) {
    read_parameter(&save_pt, "2KappaMu:", "%lf", 1, inputfile, _DEFAULT_SET);
    mu = mu * (4. + m0);
  }
#endif

#ifdef HAVE_TM1p1
  save_pt = &(epsbar);
  epsbar = 0;
  read_parameter(&save_pt, "epsbar:", "%lf", 1, inputfile, _DEFAULT_SET);
#endif
}

void Global::readSolverInfo() {

  save_pt = &(mixed_precision);
  mixed_precision = 2;
  read_parameter(&save_pt, "mixed precision:", "%d", 1, inputfile, _DEFAULT_SET);
  if (num_levels == 1)
    interpolation = 0;
  else {
    save_pt = &(interpolation);
    interpolation = 2;
    read_parameter(&save_pt, "interpolation:", "%d", 1, inputfile, _DEFAULT_SET);
  }

  save_pt = &(randomize);
  randomize = 0;
  read_parameter(&save_pt, "randomize test vectors:", "%d", 1, inputfile, _DEFAULT_SET);
  save_pt = &(coarse_iter);
  coarse_iter = 200;
  read_parameter(&save_pt, "coarse grid iterations:", "%d", 1, inputfile, _DEFAULT_SET);
  save_pt = &(coarse_restart);
  coarse_restart = 10;
  read_parameter(&save_pt, "coarse grid restarts:", "%d", 1, inputfile, _DEFAULT_SET);
  save_pt = &(coarse_tol);
  coarse_tol = 1E-1;
  read_parameter(&save_pt, "coarse grid tolerance:", "%le", 1, inputfile, _DEFAULT_SET);
  save_pt = &(odd_even);
  odd_even = 1;
  read_parameter(&save_pt, "odd even preconditioning:", "%d", 1, inputfile, _DEFAULT_SET);
#ifdef GCRODR
  save_pt = &(gcrodr_k);
  gcrodr_k = 15;
  read_parameter(&save_pt, "coarse grid gcrodr_k:", "%d", 1, inputfile, _DEFAULT_SET);
  save_pt = &(gcrodr_upd_itrs);
  gcrodr_upd_itrs = 5;
  read_parameter(&save_pt, "coarse grid gcrodr_upd_itrs:", "%d", 1, inputfile, _DEFAULT_SET);
#endif

#ifdef POLYPREC
  save_pt = &(polyprec_d);
  polyprec_d = 5;
  read_parameter(&save_pt, "coarse grid polyprec_d:", "%d", 1, inputfile, _DEFAULT_SET);
  polyprec_d++;
#endif

#ifdef BLOCK_JACOBI
  save_pt = &(local_polyprec_d);
  local_polyprec_d = 5;
  read_parameter(&save_pt, "coarse grid local_polyprec_d:", "%d", 1, inputfile, _DEFAULT_SET);
  local_polyprec_d++;
#endif

  save_pt = &(low_level_meas);
  low_level_meas = 0;
  read_parameter(&save_pt, "low level meas:", "%d", 1, inputfile, _DEFAULT_SET);

  save_pt = &(setup_m0);
  setup_m0 = m0;
  read_parameter(&save_pt, "setup m0:", "%lf", 1, inputfile, _DEFAULT_SET);
#ifdef HAVE_TM
  save_pt = &(mu_odd_shift);
  mu_odd_shift = 0;
  read_parameter(&save_pt, "mu odd shift:", "%lf", 1, inputfile, _DEFAULT_SET);
  save_pt = &(mu_even_shift);
  mu_even_shift = 0;
  read_parameter(&save_pt, "mu even shift:", "%lf", 1, inputfile, _DEFAULT_SET);
  save_pt = &(setup_mu);
  setup_mu = mu;
  read_parameter(&save_pt, "setup mu:", "%lf", 1, inputfile, _DEFAULT_SET);
#endif

#ifdef HAVE_TM1p1
  save_pt = &(epsbar_ig5_odd_shift);
  epsbar_ig5_odd_shift = 0;
  read_parameter(&save_pt, "epsbar odd shift:", "%lf", 1, inputfile, _DEFAULT_SET);
  save_pt = &(epsbar_ig5_even_shift);
  epsbar_ig5_even_shift = 0;
  read_parameter(&save_pt, "epsbar even shift:", "%lf", 1, inputfile, _DEFAULT_SET);
#endif

  save_pt = &(method);
  method = 2;
  read_parameter(&save_pt, "method:", "%d", 1, inputfile, _DEFAULT_SET);
  save_pt = &(restart);
  restart = 30;
  read_parameter(&save_pt, "iterations between restarts:", "%d", 1, inputfile, _DEFAULT_SET);
  save_pt = &(max_restart);
  max_restart = 20;
  read_parameter(&save_pt, "maximum of restarts:", "%d", 1, inputfile, _DEFAULT_SET);
  save_pt = &(tol);
  tol = 1E-10;
  read_parameter(&save_pt, "tolerance for relative residual:", "%le", 1, inputfile, _DEFAULT_SET);
  save_pt = &(print);
  print = 0;
  read_parameter(&save_pt, "print mode:", "%d", 1, inputfile, _DEFAULT_SET);
#ifdef HAVE_TM
  save_pt = &(downprop);
  downprop = 1;
  read_parameter(&save_pt, "addDownPropagator:", "%d", 1, inputfile, _DEFAULT_SET);
#endif

  // set seed for random number generator (rng)
  int rng_seed = 1234 * my_rank;
  rng_seed += randomize ? time(0) : 0;
  srand(rng_seed);
}

void Global::readGeometryInfo() {

  char inputstr[STRINGLENGTH];
  int mu, nb, nls, nlls, flag;

  for (int i = 0; i < num_levels; i++) {

    if (i > 0) {
      // global lattice
      sprintf(inputstr, "d%d global lattice:", i);
      save_pt = global_lattice[i];

      if (!read_parameter(&save_pt, inputstr, "%d", 4, inputfile, _DEFAULT_SET)) {
        nls = 1;
        for (mu = 0; mu < 4; mu++) {
          global_lattice[i][mu] = global_lattice[i - 1][mu] / block_lattice[i - 1][mu];
          nls *= global_lattice[i][mu];
        }
        if (odd_even && nls < 2) {
          warning0(
              "lattice dimensions not valid for a %d-level method, choosing a %d-level method\n",
              num_levels, i);
          num_levels = i;
          break;
        }
      }

      // local lattice
      sprintf(inputstr, "d%d local lattice:", i);
      save_pt = local_lattice[i];

      if (!read_parameter(&save_pt, inputstr, "%d", 4, inputfile, _DEFAULT_SET)) {
        nls = 1;
        nlls = 1;
        for (mu = 0; mu < 4; mu++) {
          local_lattice[i][mu] = local_lattice[i - 1][mu] / block_lattice[i - 1][mu];
          nlls *= local_lattice[i][mu];
          nls *= global_lattice[i][mu];
        }
        if (odd_even && nlls < 2) {
          if (nls / nlls > 1) {
            mu = shortest_dir(local_lattice[i]);
            if (global_lattice[i][mu] > local_lattice[i][mu]) {
              local_lattice[i][mu] *=
                  lcm(local_lattice[i][mu], global_lattice[i][mu] / local_lattice[i][mu]);
            }
          }
        }
      }

      // block lattice
      for (mu = 0; mu < 4; mu++)
        block_lattice[i][mu] = 1;
      if (i < num_levels - 1) {
        sprintf(inputstr, "d%d block lattice:", i);
        save_pt = block_lattice[i];
        if (!read_parameter(&save_pt, inputstr, "%d", 4, inputfile, _DEFAULT_SET)) {
          nls = 1;
          nb = 1;
          flag = 1;
          for (mu = 0; mu < 4; mu++) {
            if (DIVIDES(2, global_lattice[i][mu])) {
              block_lattice[i][mu] = 2;
            } else if (DIVIDES(3, global_lattice[i][mu])) {
              block_lattice[i][mu] = 3;
            } else {
              warning0("lattice dimensions not valid for a %d-level method, choosing a %d-level "
                        "method\n",
                        num_levels, i + 1);
              num_levels = i + 1;
              block_lattice[i][mu] = 1;
              flag = 0;
              break;
            }
            nb *= local_lattice[i][mu] / block_lattice[i][mu];

            if (local_lattice[i][mu] < block_lattice[i][mu]) {
              local_lattice[i][mu] *= block_lattice[i][mu];
              if (!DIVIDES(local_lattice[i][mu], global_lattice[i][mu])) {
                local_lattice[i][mu] /= block_lattice[i][mu];
              }
              warning0("lattice dimensions not valid for a %d-level method, choosing a %d-level "
                        "method\n",
                        num_levels, i + 1);
              num_levels = i + 1;
              block_lattice[i][mu] = 1;
              flag = 0;
              break;
            }
          }

          if (flag == 1 && method == 2 && nb == 1) {
            mu = shortest_dir(local_lattice[i]);
            if (global_lattice[i][mu] > local_lattice[i][mu]) {
              local_lattice[i][mu] *=
                  lcm(local_lattice[i][mu], global_lattice[i][mu] / local_lattice[i][mu]);
            }
          }
        }
      }
    }
#ifdef DEBUG
    printf0("level: %d, gl: %3d %3d %3d %3d\n", i, global_lattice[i][0], global_lattice[i][1],
            global_lattice[i][2], global_lattice[i][3]);

    printf0("level: %d, ll: %3d %3d %3d %3d\n", i, local_lattice[i][0], local_lattice[i][1],
            local_lattice[i][2], local_lattice[i][3]);

    printf0("level: %d, bl: %3d %3d %3d %3d\n\n", i, block_lattice[i][0], block_lattice[i][1],
            block_lattice[i][2], block_lattice[i][3]);
#endif

    sprintf(inputstr, "d%d post smooth iter:", i);
    save_pt = &(post_smooth_iter[i]);
    post_smooth_iter[i] = 4;
    read_parameter(&save_pt, inputstr, "%d", 1, inputfile, _DEFAULT_SET);

    sprintf(inputstr, "d%d preconditioner cycles:", i);
    save_pt = &(ncycle[i]);
    ncycle[i] = 1;
    read_parameter(&save_pt, inputstr, "%d", 1, inputfile, _DEFAULT_SET);

    sprintf(inputstr, "d%d relaxation factor:", i);
    save_pt = &(relax_fac[i]);
    relax_fac[i] = 1.0;
    read_parameter(&save_pt, inputstr, "%lf", 1, inputfile, _DEFAULT_SET);

    sprintf(inputstr, "d%d block iter:", i);
    save_pt = &(block_iter[i]);
    block_iter[i] = 4;
    read_parameter(&save_pt, inputstr, "%d", 1, inputfile, _DEFAULT_SET);

    sprintf(inputstr, "d%d setup iter:", i);
    save_pt = &(setup_iter[i]);
    if (i == 0)
      setup_iter[i] = 5;
    else if (i == 1)
      setup_iter[i] = 3;
    else if (i > 1)
      setup_iter[i] = 2;
    read_parameter(&save_pt, inputstr, "%d", 1, inputfile, _DEFAULT_SET);

#ifdef HAVE_TM
    sprintf(inputstr, "d%d mu factor:", i);
    save_pt = &(mu_factor[i]);
    mu_factor[i] = 1;
    read_parameter(&save_pt, inputstr, "%lf", 1, inputfile, _DEFAULT_SET);
#endif

#ifdef HAVE_TM1p1
    sprintf(inputstr, "d%d epsbar factor:", i);
    save_pt = &(epsbar_factor[i]);
    epsbar_factor[i] = 1;
    read_parameter(&save_pt, inputstr, "%lf", 1, inputfile, _DEFAULT_SET);
#endif

    sprintf(inputstr, "d%d test vectors:", i);
    save_pt = &(num_eig_vect[i]);
    if (i == 0)
      num_eig_vect[i] = 24;
    else
      num_eig_vect[i] = 28;
    read_parameter(&save_pt, inputstr, "%d", 1, inputfile, _DEFAULT_SET);
  }
}

void Global::readTestvectorIOData() {

  if (interpolation == 2) {
    save_pt = &(tv_io_single_file);
    tv_io_single_file = 1;
    read_parameter(&save_pt, "test vector io from single file:", "%d", 1, inputfile,
                    _DEFAULT_SET);
    save_pt = &(tv_io_file_name);
    read_parameter(&save_pt, "test vector io file name:", "%s", 1, inputfile, _NO_DEFAULT_SET);
  }
}

void Global::readEvaluationParameters() {

  save_pt = &(vt.evaluation);
  vt.evaluation = 0;
  read_parameter(&save_pt, "evaluation:", "%d", 1, inputfile, _DEFAULT_SET);
  if (vt.evaluation) {
    save_pt = &(vt.scan_var);
    read_parameter(&save_pt, "scan variable:", "%s", 1, inputfile, _NO_DEFAULT_SET);
    save_pt = &(vt.start_val);
    read_parameter(&save_pt, "start value:", "%lf", 1, inputfile, _NO_DEFAULT_SET);
    save_pt = &(vt.end_val);
    read_parameter(&save_pt, "end value:", "%lf", 1, inputfile, _NO_DEFAULT_SET);
    save_pt = &(vt.step_size);
    read_parameter(&save_pt, "step size:", "%lf", 1, inputfile, _NO_DEFAULT_SET);
    save_pt = &(vt.multiplicative);
    read_parameter(&save_pt, "multiplicative:", "%d", 1, inputfile, _NO_DEFAULT_SET);
    save_pt = &(vt.shift_update);
    read_parameter(&save_pt, "shift update:", "%d", 1, inputfile, _NO_DEFAULT_SET);
    save_pt = &(vt.re_setup);
    read_parameter(&save_pt, "setup update:", "%d", 1, inputfile, _NO_DEFAULT_SET);
    save_pt = &(vt.track_error);
    vt.track_error = 0;
    read_parameter(&save_pt, "track error:", "%d", 1, inputfile, _NO_DEFAULT_SET);
    save_pt = &(vt.track_cgn_error);
    vt.track_cgn_error = 0;
    read_parameter(&save_pt, "compare with CGN error:", "%d", 1, inputfile, _DEFAULT_SET);
    save_pt = &(vt.average_over);
    vt.average_over = 1;
    read_parameter(&save_pt, "average over:", "%d", 1, inputfile, _DEFAULT_SET);
  }
}

void Global::readKcycleData() {

  save_pt = &(kcycle);
  kcycle = 1;
  read_parameter(&save_pt, "kcycle:", "%d", 1, inputfile, _DEFAULT_SET);
  save_pt = &(kcycle_restart);
  kcycle_restart = 5;
  read_parameter(&save_pt, "kcycle length:", "%d", 1, inputfile, _DEFAULT_SET);
  save_pt = &(kcycle_max_restart);
  kcycle_max_restart = 2;
  read_parameter(&save_pt, "kcycle restarts:", "%d", 1, inputfile, _DEFAULT_SET);
  save_pt = &(kcycle_tol);
  kcycle_tol = 1E-1;
  read_parameter(&save_pt, "kcycle tolerance:", "%le", 1, inputfile, _DEFAULT_SET);
}

void Global::validateParameters() {

  int i;
  int mu;

#ifdef SSE
  if (!odd_even)
    warning0("The SSE implementation is based on the odd-even preconditioned code.\
  \n         Switch on odd-even preconditioning in the input file.\n");
  ASSERT(odd_even);
#endif

  ASSERT(ASCENDING(0, rhs, 2));
  ASSERT(ASCENDING(-1, method, 5));

  ASSERT(IMPLIES(vt.evaluation, rhs <= 2));
#ifdef _20TV
  ASSERT(num_eig_vect == 20);
#endif

  for (i = 0; i < num_levels; i++)
    for (mu = 0; mu < 4; mu++)
      ASSERT(DIVIDES(local_lattice[i][mu], global_lattice[i][mu]));

  for (i = 0; i < num_levels - 1; i++)
    for (mu = 0; mu < 4; mu++) {
      ASSERT(DIVIDES(global_lattice[i + 1][mu], global_lattice[i][mu]));
      ASSERT(DIVIDES(block_lattice[i][mu], local_lattice[i][mu]));
      ASSERT(DIVIDES(global_lattice[i][mu] / global_lattice[i + 1][mu], local_lattice[i][mu]));
      ASSERT(DIVIDES(block_lattice[i][mu], global_lattice[i][mu] / global_lattice[i + 1][mu]));
#ifdef SSE
      if (block_lattice[i][mu] != global_lattice[i][mu] / global_lattice[i + 1][mu])
        warning0("when using SSE, Schwarz block size and aggregate size have to match.\n");
      ASSERT(block_lattice[i][mu] == global_lattice[i][mu] / global_lattice[i + 1][mu]);
      // it works everywhere but we have some problem with the vector size.
      // TODO: check all vectors allocated with size level->inner_vector_size
      ASSERT(num_eig_vect[i] % SIMD_LENGTH_double == 0);
#endif
    }

  if (odd_even) {
    int coarse_sites_per_core = 1;
    for (mu = 0; mu < 4; mu++) {
      ASSERT(DIVIDES(2, global_lattice[num_levels - 1][mu]));
      coarse_sites_per_core *= local_lattice[num_levels - 1][mu];
    }
    ASSERT(DIVIDES(2, coarse_sites_per_core));
  }

  if (method == 2) {
    for (i = 0; i < num_levels - 1; i++) {
      int num_blocks = 1;
      for (mu = 0; mu < 4; mu++) {
        num_blocks *= (local_lattice[i][mu] / block_lattice[i][mu]);
      }
      ASSERT(num_blocks >= 2);
    }
  }

  for (mu = 0; mu < 4; mu++)
    ASSERT(IMPLIES(rhs == 3, ASCENDING(0, propagator_coords[mu], global_lattice[0][mu] - 1)));

  ASSERT(IMPLIES(method > 0 && interpolation > 0, coarse_iter > 0));
  ASSERT(IMPLIES(method > 0 && interpolation > 0, coarse_restart > 0));
  ASSERT(IMPLIES(method > 0 && interpolation > 0, 0 < coarse_tol && coarse_tol < 1));
  //   ASSERT(IMPLIES(method > 0, level->n_cy > 0));
  ASSERT(max_restart > 0);
  ASSERT(0 < tol && tol < 1);
  ASSERT(ASCENDING(0, kcycle, 1));
  ASSERT(IMPLIES(kcycle && method > 0, kcycle_restart > 0));
  ASSERT(IMPLIES(kcycle && method > 0, kcycle_max_restart > 0));
  ASSERT(IMPLIES(kcycle && method > 0, 0 < kcycle_tol && kcycle_tol < 1));

  // LIST OF CASES WHICH SHOULD WORK, BUT DO NOT (TODO)

#ifdef SSE
  ASSERT(mixed_precision);
#endif

  // TODO: Could work without, but you need to fix the setup phase.
  for (i = 0; i < num_levels - 2; i++)
    ASSERT(num_eig_vect[i] <= num_eig_vect[i + 1]);

    // TODO: for some reason mixed_precision=0 do not work with num_levels>2
    //   if (num_levels > 2 && interpolation)
    //     ASSERT(mixed_precision);

#ifdef HAVE_TM1p1
  // TODO: method = 6 not supported with HAVE_TM1p1. To fix all the g5D functions
  ASSERT(method != 6);
#endif
}

void Global::initialize() {

  // compute number of MPI processes
  num_processes = 1;
  for (int i = 0; i < 4; i++) {
    process_grid[i] = global_lattice[0][i] / local_lattice[0][i];
    num_processes *= process_grid[i];
    periodic_bc[i] = 1;
    comm_offset[i] = 1;
  }

  // setup communicators
  Cart_rank = MPI_Cart_rank;
  Cart_coords = MPI_Cart_coords;
  cart_define(MPI_COMM_WORLD, this);

  // compute vector sizes and number of lattice sites for finest lattice
  num_lattice_site_var = 12;
  num_inner_lattice_sites = 1;
  for (int i = 0; i < 4; i++)
    num_inner_lattice_sites *= local_lattice[0][i];
  num_lattice_sites = num_inner_lattice_sites;
  inner_vector_size = num_inner_lattice_sites * num_lattice_site_var;

  for (int i = 0; i < 4; i++)
    num_lattice_sites += num_inner_lattice_sites / local_lattice[0][i];

  vector_size = num_lattice_sites * num_lattice_site_var;
  schwarz_vector_size = 2 * vector_size - inner_vector_size;

  // setup operator for finest level
  setup_m0 = m0;
  operator_double_init(&(op_double));
  operator_double_alloc(&(op_double), _ORDINARY, this);
  operator_double_define(&(op_double), this);

  // define odd even odd even table
  int t, z, y, x, k = 0, oe_offset = 0;
  for (int i = 0; i < 4; i++)
    oe_offset += (local_lattice[0][i] * (my_coords[i] / comm_offset[i])) % 2;
  oe_offset = oe_offset % 2;

  odd_even_table = new int[num_inner_lattice_sites];
  for (t = 0; t < local_lattice[0][_T]; t++)
    for (z = 0; z < local_lattice[0][_Z]; z++)
      for (y = 0; y < local_lattice[0][_Y]; y++)
        for (x = 0; x < local_lattice[0][_X]; x++) {
          odd_even_table[k] = ((t + z + y + x + oe_offset) % 2 == 1) ? _ODD : _EVEN;
          k++;
        }
}

void Global::setupFineDiracOperator() {

  double t0, t1;
  int i, j, t, z, y, x, mu;
  int *ll = local_lattice[0], onb[4];

  SU3_storage U = nullptr;
  complex_double phase[4];

  // read configuration
  config_double hopp = new complex_double[3 * inner_vector_size];
  if (in_format == _LIME)
    lime_read_conf((double *)(hopp), in, &(plaq_hopp), this);
  else
    read_conf((double *)(hopp), in, &(plaq_hopp), this);

  for (i = 0; i < 4; i++)
    onb[i] = (my_coords[i] == process_grid[i] - 1) ? 1 : 0;

  if (print > 0)
    printf0("%s\n", CLIFFORD_BASIS);
  if (bc == _ANTIPERIODIC)
    printf0("antiperiodic in time");
  else if (bc == _TWISTED)
    printf0("twisted (%.2f, %.2f, %.2f, %.2f)", twisted_bc[0], twisted_bc[1], twisted_bc[2],
            twisted_bc[3]);
  else
    printf0("periodic in time");
  printf0(" boundary conditions\n");

  t0 = MPI_Wtime();

  // read and store configuration
  SU3_storage_alloc(&U, local_lattice[0]);

  if (bc == _ANTIPERIODIC && onb[_T]) {
    phase[_Z] = 1;
    phase[_Y] = 1;
    phase[_X] = 1;
    for (t = 1, i = 0; t < ll[_T] + 1; t++) {
      if (t < ll[_T])
        phase[_T] = 1;
      else
        phase[_T] = -1;
      for (z = 1; z < ll[_Z] + 1; z++)
        for (y = 1; y < ll[_Y] + 1; y++)
          for (x = 1; x < ll[_X] + 1; x++)
            for (mu = 0; mu < 4; mu++)
              for (j = 0; j < 9; j++, i++) {
                op_double.D[i] = 0.5 * phase[mu] * hopp[i];
                U[t][z][y][x][mu][j] = phase[mu] * hopp[i];
              }
    }
  } else if (bc == _TWISTED && (onb[_T] || onb[_Z] || onb[_Y] || onb[_X])) {
    // TODO
    warning0("Twisted boundary conditions not supported outside the library.\n");
    for (t = 1, i = 0; t < ll[_T] + 1; t++) {
      if (!onb[_T] || t < ll[_T] || twisted_bc[_T] == 0)
        phase[_T] = 1;
      else
        phase[_T] = exp(I * twisted_bc[_T]);
      for (z = 1; z < ll[_Z] + 1; z++) {
        if (!onb[_Z] || z < ll[_Z] || twisted_bc[_Z] == 0)
          phase[_Z] = 1;
        else
          phase[_Z] = exp(I * twisted_bc[_Z]);
        for (y = 1; y < ll[_Y] + 1; y++) {
          if (!onb[_Y] || y < ll[_Y] || twisted_bc[_Y] == 0)
            phase[_Y] = 1;
          else
            phase[_Y] = exp(-I * twisted_bc[_Y]);
          for (x = 1; x < ll[_X] + 1; x++) {
            if (!onb[_X] || x < ll[_X] || twisted_bc[_X] == 0)
              phase[_X] = 1;
            else
              phase[_X] = exp(-I * twisted_bc[_X]);
            for (mu = 0; mu < 4; mu++)
              for (j = 0; j < 9; j++, i++) {
                op_double.D[i] = 0.5 * phase[mu] * hopp[i];
                U[t][z][y][x][mu][j] = phase[mu] * hopp[i];
              }
          }
        }
      }
    }
  }

  else
    for (t = 1, i = 0; t < ll[_T] + 1; t++)
      for (z = 1; z < ll[_Z] + 1; z++)
        for (y = 1; y < ll[_Y] + 1; y++)
          for (x = 1; x < ll[_X] + 1; x++)
            for (mu = 0; mu < 4; mu++)
              for (j = 0; j < 9; j++, i++) {
                op_double.D[i] = 0.5 * hopp[i];
                U[t][z][y][x][mu][j] = hopp[i];
              }

  delete[] hopp;

  SU3_ghost_update(&U, this);
  if (print > 0)
    printf0("Configuration stored...\n");

  // define gamma matrices and compute clover term
  for (i = 0; i < 4; i++)
    for (j = 0; j < 16; j++)
      gamma[i][j] = 0;

  gamma[_T][GAMMA_T_SPIN0_CO] = GAMMA_T_SPIN0_VAL;
  gamma[_T][4 + GAMMA_T_SPIN1_CO] = GAMMA_T_SPIN1_VAL;
  gamma[_T][8 + GAMMA_T_SPIN2_CO] = GAMMA_T_SPIN2_VAL;
  gamma[_T][12 + GAMMA_T_SPIN3_CO] = GAMMA_T_SPIN3_VAL;

  gamma[_Z][GAMMA_Z_SPIN0_CO] = GAMMA_Z_SPIN0_VAL;
  gamma[_Z][4 + GAMMA_Z_SPIN1_CO] = GAMMA_Z_SPIN1_VAL;
  gamma[_Z][8 + GAMMA_Z_SPIN2_CO] = GAMMA_Z_SPIN2_VAL;
  gamma[_Z][12 + GAMMA_Z_SPIN3_CO] = GAMMA_Z_SPIN3_VAL;

  gamma[_Y][GAMMA_Y_SPIN0_CO] = GAMMA_Y_SPIN0_VAL;
  gamma[_Y][4 + GAMMA_Y_SPIN1_CO] = GAMMA_Y_SPIN1_VAL;
  gamma[_Y][8 + GAMMA_Y_SPIN2_CO] = GAMMA_Y_SPIN2_VAL;
  gamma[_Y][12 + GAMMA_Y_SPIN3_CO] = GAMMA_Y_SPIN3_VAL;

  gamma[_X][GAMMA_X_SPIN0_CO] = GAMMA_X_SPIN0_VAL;
  gamma[_X][4 + GAMMA_X_SPIN1_CO] = GAMMA_X_SPIN1_VAL;
  gamma[_X][8 + GAMMA_X_SPIN2_CO] = GAMMA_X_SPIN2_VAL;
  gamma[_X][12 + GAMMA_X_SPIN3_CO] = GAMMA_X_SPIN3_VAL;

  compute_clover_term(U, this);

  // calculate the plaquette
  plaq_clov = calc_plaq(U, this);
  if (print > 0)
    printf0("average plaquette: %.13lf\n", plaq_clov);

  SU3_storage_free(&U, local_lattice[0]);

  t1 = MPI_Wtime();
  if (print > 0) {
    printf0("\n+----------------------------------------------------------+\n");
    printf0("| read in and set up the parallel dirac operator           |\n");
    printf0("| elapsed wall clock time: %-12g seconds            |\n", t1 - t0);
    printf0("+----------------------------------------------------------+\n");
  }
}

// void Global::printParams() {
//   printf0("num_processes: %d\n", num_processes);
//   for (int mu = 0; mu < 4; mu++)
//     printf("my_rank: %d my_coords[%d]: %d\n", my_rank, mu, my_coords[mu]);
//   printf0("tv_io_single_file: %d\n", tv_io_single_file);
//   printf0("num_openmp_processes: %d\n", num_openmp_processes);
//   printf0("num_levels: %d\n", num_levels);
//   printf0("num_desired_levels: %d\n", num_desired_levels);
//   for (int mu = 0; mu < 4; mu++)
//     printf0("process_grid[%d]: %d\n", mu, process_grid[mu]);
//   for (int i = 0; i < num_levels; i++)
//     for (int mu = 0; mu < 4; mu++)
//       printf0("global_lattice[%d][%d]: %d\n", i, mu, global_lattice[i][mu]);
//   for (int i = 0; i < num_levels; i++)
//     for (int mu = 0; mu < 4; mu++)
//       printf("my_rank: %d local_lattice[%d][%d]: %d\n", my_rank, i, mu, local_lattice[i][mu]);
//   for (int i = 0; i < num_levels; i++)
//     for (int mu = 0; mu < 4; mu++)
//       printf0("block_lattice[%d][%d]: %d\n", i, mu, block_lattice[i][mu]);
//   for (int i = 0; i < num_levels; i++)
//     printf0("post_smooth_iter[%d]: %d\n", i, post_smooth_iter[i]);
//   for (int i = 0; i < num_levels; i++)
//     printf0("block_iter[%d]: %d\n", i, block_iter[i]);
//   for (int i = 0; i < num_levels; i++)
//     printf0("setup_iter[%d]: %d\n", i, setup_iter[i]);
//   for (int i = 0; i < num_levels; i++)
//     printf0("ncycle[%d]: %d\n", i, ncycle[i]);
//   printf0("method: %d\n", method);
//   printf0("odd_even: %d\n", odd_even);
//   printf0("rhs: %d\n", rhs);
//   for (int mu = 0; mu < 4; mu++)
//     printf0("propagator_coords[%d]: %d\n", mu, propagator_coords[mu]);
//   printf0("interpolation: %d\n\n", interpolation);
//   printf0("randomize: %d\n", randomize);
//   for (int i = 0; i < num_levels; i++)
//     printf0("num_eig_vect[%d]: %d\n", i, num_eig_vect[i]);
//   printf0("num_coarse_eig_vect: %d\n", num_coarse_eig_vect);
//   printf0("kcycle: %d\n", kcycle);
//   printf0("mixed_precision: %d\n", mixed_precision);
//   printf0("restart: %d\n", restart);
//   printf0("max_restart: %d\n", max_restart);
//   printf0("kcycle_restart: %d\n", kcycle_restart);
//   printf0("kcycle_max_restart: %d\n", kcycle_max_restart);
//   printf0("coarse_iter: %d\n", coarse_iter);
//   printf0("coarse_restart: %d\n", coarse_restart);
//   printf0("tol: %lf\n", tol);
//   printf0("coarse_tol: %lf\n", coarse_tol);
//   printf0("kcycle_tol: %lf\n", kcycle_tol);
//   printf0("csw: %lf\n", csw);
//   printf0("rho: %lf\n", rho);
//   for (int i = 0; i < num_levels; i++)
//     printf0("relax_fac[%d]: %lf\n", i, relax_fac[i]);
//   printf0("m0: %lf\n", m0);
//   printf0("setup_m0: %lf\n", setup_m0);
//   printf0("bc: %d\n", bc);
//   for (int mu = 0; mu < 4; mu++)
//     printf0("periodic_bc[%d]: %d\n", mu, periodic_bc[mu]);
//   for (int mu = 0; mu < 4; mu++)
//     printf0("num_processes_dir[%d]: %d\n", mu, num_processes_dir[mu]);
//   for (int mu = 0; mu < 4; mu++)
//     printf0("comm_offset[%d]: %d\n", mu, comm_offset[mu]);
//   for (int mu = 0; mu < 8; mu++)
//     printf0("neighbor_rank[%d]: %d\n", mu, neighbor_rank[mu]);
//   for (int mu = 0; mu < 4; mu++)
//     printf0("process_grid[%d]: %d\n", mu, process_grid[mu]);
//   printf0("num_lattice_site_var: %d\n", num_lattice_site_var);
//   printf0("num_lattice_sites: %d\n", num_lattice_sites);
//   printf0("num_inner_lattice_sites: %d\n", num_inner_lattice_sites);
//   printf0("inner_vector_size: %d\n", inner_vector_size);
//   printf0("vector_size: %d\n", vector_size);
//   printf0("schwarz_vector_size: %d\n", schwarz_vector_size);
// }
