
#include "main.h"


/**
 * @brief Allocating memory for Solver class.
 *
 */
void Solver::allocateMemory() {
  x = new complex_double[vectorLength];
  b = new complex_double[vectorLength];
  r = new complex_double[vectorLength];
  w = new complex_double[vectorLength];

  NEW2D(V, complex_double, vectorLength, restartLength + 1); // TODO: Get rid of the +1's?
  NEW2D(Z, complex_double, vectorLength, restartLength + 1);
  //NEW2D(H, complex_double, restartLength + 1, restartLength + 1);
}

Solver::Solver()
    : Solver(nullptr, 0, 0, 0, 0.0, 0, 0, 0, nullptr, nullptr, nullptr) {
} // delegating default constructor

Solver::Solver(operator_struct<double> *matrix, int resLen, int resNum, int probSize, double tol,
               int t, int print, int vecLen,
               void (*precond)(vector_double, vector_double, vector_double, int, Level *),
               void (*evalOp)(vector_double, vector_double, operator_struct<double> *, Level *),
               Level *lev)
    : preconditioner(precond), applyOperator(evalOp), level(lev), matrix(matrix), tolerance(tol),
      problemSize(probSize), restartNumber(resNum), restartLength(resLen), timing(t),
      verbose(print), vectorLength(vecLen) {
  Solver::allocateMemory();
}

Solver::~Solver() {
  delete[] x;
  delete[] b;
  delete[] r;
  delete[] w;
  DELETE2D(V);
  DELETE2D(Z);
  //DELETE2D(H);
}

void Gmres::allocateMemory() {
  sinus = new complex_double[restartLength + 1]; // TODO: Get rid of the +1's?
  cosinus = new complex_double[restartLength + 1];
  gamma = new complex_double[restartLength + 1];
  y = new complex_double[restartLength + 2];
}

void Gmres::initialize(operator_struct<double> *op, int resLen, int resNum, int probSize, double tol,
                       int t, int precKind, int tim, int print, int vecLen,
                       void (*precond)(vector_double, vector_double, vector_double, int, Level *),
                       void (*evalOp)(vector_double, vector_double, operator_struct<double> *,
                                      Level *),
                       Level *lev) {

  matrix = op;
  restartLength = resLen;
  restartNumber = resNum;
  problemSize = probSize;
  tolerance = tol;
  type = t;
  preconditionerKind = precKind;
  timing = tim;
  verbose = print;
  vectorLength = vecLen;
  vectorEnd = problemSize;
  preconditioner = precond;
  applyOperator = evalOp;
  level = lev;

  // free pointers in case of re-allocation
  delete[] sinus;
  delete[] cosinus;
  delete[] gamma;
  delete[] y;

  delete[] x;
  delete[] b;
  delete[] r;
  delete[] w;
  DELETE2D(V);
  DELETE2D(Z);
  //DELETE2D(H);

  Solver::allocateMemory();
  Gmres::allocateMemory();
}

Gmres::Gmres()
    : Gmres(nullptr, 0, 0, 0, 0.0, 0, 0, 0, 0, 0, nullptr, nullptr, nullptr) {
} // delegating default constructor

Gmres::Gmres(operator_struct<double> *matrix, int resLen, int resNum, int probSize, double tol,
             int t, int precKind, int tim, int print, int vecLen,
             void (*precond)(vector_double, vector_double, vector_double, int, Level *),
             void (*evalOp)(vector_double, vector_double, operator_struct<double> *, Level *),
             Level *lev)
    : Solver(matrix, resLen, resNum, probSize, tol, tim, print, vecLen, precond, evalOp, lev),
      type(t), preconditionerKind(precKind) {
  Gmres::allocateMemory();
}

Gmres::~Gmres() {
  delete[] sinus;
  delete[] cosinus;
  delete[] gamma;
  delete[] y;
}

static inline int idx(int i, int j){
  return (i*i + 3*i)/2+j;
}

/**
 * Implementation of Arnoldi method
 * @return
 */
int Gmres::arnoldiStep(complex_double *H) {
  // TODO: Explain why the arnoldi step always return 1, the c code returns 0 in line 793 linsolve_generic.c

  Global *global = level->global;

  if (preconditioner != nullptr) {
    if (preconditionerKind == _LEFT) {
      /// Z[0] = A*V[j]
      applyOperator(Z[0], V[j], matrix, level);
      /// w = M^-1 * Z[0] = M^-1 * A * V[j]
      preconditioner(w, nullptr, Z[0], _NO_RES, level);
    } else {
      //TODO: note that mixed_precision var isn't used in the program, consider deleting it
      if (global->mixed_precision == 2 && (global->method == 1 || global->method == 2)) {
        preconditioner(Z[j], w, V[j], _NO_RES, level);
      } else {
        /// Z[j] = M^-1 * V[j]
        preconditioner(Z[j], nullptr, V[j], _NO_RES, level);
        /// w = A*Z[j]
        applyOperator(w, Z[j], matrix, level);
      }
    }
  } else {
    /// w = A*V[j]
    applyOperator(w, V[j], matrix, level);
  }

  /// orthogonalization
  auto *tmp = new complex_double[j + 2];
  /// tmp[c] = conj(V[c][i]) * w[i], where c in 0...j and i in vectorStart..vectorEnd
  process_multi_inner_product_double(j + 1, tmp, V, w, vectorStart, vectorEnd, level);

  /// y[i] = <V[i], w>, i in 0...j
  for (int i = 0; i <= j; i++) {
    y[i] = tmp[i];
  }

  /// H[j][i] = y[i] , where i in 0...j
  if (global->num_processes > 1) {
    MPI_Allreduce(y, &H[idx(j,0)], j + 1, MPI_COMPLEX_double, MPI_SUM,
                  (level->depth == 0) ? global->comm_cart : level->gs_double.level_comm);
  } else {
    for (int i = 0; i <= j; i++)
      H[idx(j,i)] = y[i];
  }

  /// w = w - H[j][i] * V[i]
  for (int i = 0; i <= j; i++)
    vector_double_saxpy(w, w, V[i], -H[idx(j,i)], vectorStart, vectorEnd, level);
#ifdef REORTH
// TODO Redo Reorth
#endif

  /// H[j][j+1] = ||w||_2
  H[idx(j,j+1)] = global_norm_double(w, vectorStart, vectorEnd, level);

  /// V[j+1] = w / H[j][j+1]
  if (abs(H[idx(j,j+1)]) > 1e-15)
    vector_double_scale(V[j + 1], w, real(1.0 / H[idx(j,j+1)]), vectorStart, vectorEnd, level);

  delete[] tmp;

  return 1;
}

/**
 * Update QR factorization and apply previous Givens rotation.
 */
void Gmres::qrUpdate(complex_double *H) {
  /// apply Givens rotation matrix [-s c |c_bar s_bar] to H up to j
  for (int i = 0; i < j; i++) {
    beta = (-sinus[i]) * H[idx(j,i)] + (cosinus[i]) * H[idx(j, i+1)];
    H[idx(j,i)] = conj(cosinus[i]) * H[idx(j, i)] + conj(sinus[i]) * H[idx(j, i+1)];
    H[idx(j,i+1)] = beta;
  }
  /// compute current Givens rotation
  beta = (complex_double)sqrt(NORM_SQUARE_double(H[idx(j,j)]) + NORM_SQUARE_double(H[idx(j,j+1)]));
  sinus[j] = H[idx(j,j+1)] / beta;
  cosinus[j] = H[idx(j,j)] / beta;
  /// update right column
  gamma[j + 1] = (-sinus[j]) * gamma[j];
  gamma[j] = conj(cosinus[j]) * gamma[j];
  /// apply current Givens rotation
  H[idx(j,j)] = beta;
  H[idx(j, j+1)] = 0;
}

void Gmres::computeSolution(int outerLoop, complex_double *H) {
  /// if it is GMRES with right preconditioning, x is updated with M^-1 * V
  vector_double *basis = (preconditioner != nullptr && preconditionerKind == _RIGHT) ? Z : V;
  /// backward substitution
  for (int i = j; i >= 0; i--) {
    y[i] = gamma[i];
    for (int k = i + 1; k <= j; k++) {
      y[i] -= H[idx(k, i)] * y[k];
    }
    y[i] /= H[idx(i, i)];
  }

  /// x = x + V*y (or x = x+Z*y for FGMRES)
  if (outerLoop) {
    for (int i = 0; i <= j; i++) {
      vector_double_saxpy(x, x, basis[i], y[i], vectorStart, vectorEnd, level);
    }
  } else {
    vector_double_scale(x, basis[0], y[0], vectorStart, vectorEnd, level);
    for (int i = 1; i <= j; i++) {
      vector_double_saxpy(x, x, basis[i], y[i], vectorStart, vectorEnd, level);
    }
  }
}

vector_double Gmres::solve() {

  Global *global = level->global;
  auto *H = new complex_double[((restartLength+1)*(restartLength+1) + 3*(restartLength+1))/2];

  if (level->type == _FINE)
    elapsedTime = -MPI_Wtime();

  if (verbose && global->print > 0)
    printf0("\n+----------------------------------------------------------+\n");
  for (int outerLoop = 0; outerLoop < restartNumber && finished == 0; outerLoop++) {
    if (outerLoop == 0 && initialGuessIsZero) {
      /// set residual flag to no residual
      res = _NO_RES;
      /// set r = b
      vector_double_copy(r, b, vectorStart, vectorEnd, level);
    } else {
      res = _RES;
      if (preconditioner != nullptr && preconditionerKind == _LEFT) {
        /// Z[0] = Ax
        applyOperator(Z[0], x, matrix, level);
        /// w = M^-1 * Z[0] = M^-1 * Ax
        preconditioner(w, nullptr, Z[0], _NO_RES, level);
      } else {
        /// w = Ax
        applyOperator(w, x, matrix, level);
      }
      /// r = b - w
      vector_double_minus(r, b, w, vectorStart, vectorEnd, level);
    }
    /// norm = ||r||_2 Euclidean norm of r
    norm = global_norm_double(r, vectorStart, vectorEnd, level);
    gamma[0] = norm; // gamma_0 = norm(r)

    if (outerLoop == 0) {

      if (level->type == _FINE && !initialGuessIsZero) {
        normR0 = global_norm_double(b, vectorStart, vectorEnd, level);
      } else {
        normR0 = norm;
      }
    }
    /// V[0] = 1/gamma[0] * r
    vector_double_scale(V[0], r, 1.0 / gamma[0], vectorStart, vectorEnd,
                        level); // v_0 = r / gamma_0

    for (int innerLoop = 0; innerLoop < restartLength && finished == 0; innerLoop++) {
      j = innerLoop;
      iter++;

      /// simply do one step of Arnoldi
      if (!arnoldiStep(H)) {
        printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j + 1, j);
        break;
      }

      if (abs(H[idx(j,j+1)]) > tolerance / 10) {
        qrUpdate(H);
        gammaJPlusOne = abs(gamma[j + 1]);

#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
        if (iter % 10 == 0 || preconditioner != nullptr || level->depth > 0) {

          if (verbose && global->print > 0)
            printf0("| approx. rel. res. after  %-6d iterations: %e |\n", iter,
                    gammaJPlusOne / normR0);
        }
#endif

        if (gammaJPlusOne / normR0 < tolerance ||
            gammaJPlusOne / normR0 > 1E+5) { // if satisfied ... stop
          finished = 1;
          if (gammaJPlusOne / normR0 > 1E+5)
            printf0("Divergence of fgmres_double, iter = %d, depth=%d\n", iter, level->depth);
        }
      } else {
        finished = 1;
        break;
      }
    } // end of a single restart
    computeSolution((res == _NO_RES) ? outerLoop : 1, H);

  } // end of fgmres

  if (level->type == _FINE) {
    elapsedTime += MPI_Wtime();
    global->total_time = elapsedTime;
    global->iter_count = iter;
    global->norm_res = gammaJPlusOne / normR0;
  }

  if (verbose) {
#ifdef FGMRES_RESTEST
    applyOperator(w, x, matrix, level);
    vector_double_minus(r, b, w, vectorStart, vectorEnd, level);
    norm = global_norm_double(r, vectorStart, vectorEnd, level);
#else
    norm = gammaJPlusOne;
#endif
    global->norm_res = norm / normR0;
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    if (global->print > 0)

      printf0("+----------------------------------------------------------+\n\n");
#endif
    printf0("+----------------------------------------------------------+\n");
    printf0("|       FGMRES iterations: %-6d coarse average: 0        |\n", iter); // TODO Fix coarse average
    printf0("| exact relative residual: ||r||/||b|| = %e      |\n", norm / normR0);
    printf0("| elapsed wall clock time: %-8.4lf seconds                |\n", elapsedTime);
    if (global->coarse_time > 0)
      printf0("|        coarse grid time: %-8.4lf seconds (%04.1lf%%)        |\n",
              global->coarse_time, 100 * (global->coarse_time / elapsedTime));
    printf0("|      coarsest grid time: %-8.4lf seconds (%04.1lf%%)        |\n",
            global->coarsest_time, 100 * (global->coarsest_time / elapsedTime));
    printf0("|  consumed core minutes: %-8.2le (solve only)            |\n",
            (elapsedTime * global->num_processes / 60.0));
    printf0("|    max used mem/MPIproc: %-8.2le GB                     |\n",
            0.0); // TODO: Fix mem usage(?)
    printf0("+----------------------------------------------------------+\n");
  }

  if (level->type == _COARSEST) {
    global->coarse_iter_count += iter;
  }
  if (level->type == _FINE && (timing || verbose) && global->vt.p_end == nullptr) {
    if (global->method != 6)
      prof_print(level);
  }

  /// clean up / reset relevant class members
  finished = 0;
  iter = 0;
  delete[] H;

  return x;
}

void Schwarz::allocateMemory() {

  Global *global = level->global;

  // large vectors
  Dphi = new complex_double[vectorLength];
  latestIter = new complex_double[vectorLength];
  latestIter2 = new complex_double[vectorLength];

  // block index fields
  index_table[_T] = new int[MAX(1, direction_length[_T] + direction_length[_Z] +
                                       direction_length[_Y] + direction_length[_X])];
  index_table[_Z] = index_table[_T] + direction_length[_T];
  index_table[_Y] = index_table[_Z] + direction_length[_Z];
  index_table[_X] = index_table[_Y] + direction_length[_Y];
  if (level->type == _FINE && level->odd_even) {
    odd_even_index_table[_T] = new int[MAX(1, direction_length[_T] + direction_length[_Z] +
                                                  direction_length[_Y] + direction_length[_X])];
    odd_even_index_table[_Z] = odd_even_index_table[_T] + direction_length[_T];
    odd_even_index_table[_Y] = odd_even_index_table[_Z] + direction_length[_Z];
    odd_even_index_table[_X] = odd_even_index_table[_Y] + direction_length[_Y];
  }

  // create list of Schwarz blocks
  if (global->method == _SIXTEEN_COLOR) {
    block_list = new int *[16];
    block_list[0] = new int[number_of_blocks];
    int j = number_of_blocks / 16;
    for (int i = 1; i < 16; i++)
      block_list[i] = block_list[0] + i * j;
  } else if (global->method == _TWO_COLOR) {
    block_list_length = new int[8];
    block_list = new int *[8];
    for (int i = 0; i < 8; i++) {
      block_list[i] = new int[number_of_blocks];
    }
  }

  operator_double_init(&op);
  if (level->type != _COARSEST)
    operator_double_alloc(&op, _SCHWARZ, level);
  else
    operator_double_alloc(&op, _ORDINARY, level);

  block = new block_struct[number_of_blocks];
  int n = 0;
  for (int mu = 0; mu < 4; mu++) {
    int i = 1;
    for (int nu = 0; nu < 4; nu++) {
      if (mu != nu) {
        i *= level->block_lattice[nu];
      }
    }
    block_boundary_length[2 * mu] = n;
    block_boundary_length[2 * mu + 1] = n + 2 * i;
    n += 4 * i;
  }
  block_boundary_length[8] = n;

  for (int i = 0; i < number_of_blocks; i++) {
    block[i].bt = new int[n];
  }

  // buffer vectors
  int svs = level->schwarz_vector_size,
      vs = (level->type == _FINE) ? level->inner_vector_size : level->vector_size;

#ifdef HAVE_TM1p1
  svs *= 2;
  vs *= 2;
#endif

  if (level->type == _FINE) {
    NEW2D(odd_even_buffer, complex_double, vs, 4);
  }
  NEW2D(local_minres_buffer, complex_double, svs, 3);

  // TODO: new hugepages macro
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
  if (level->type == _FINE) {
    MALLOC_HUGEPAGES(op.D_vectorized, double,
                     2 * 4 * (2 * level->vector_size - level->inner_vector_size),
                     4 * SIMD_LENGTH_double);
    MALLOC_HUGEPAGES(op.D_transformed_vectorized, double,
                     2 * 4 * (2 * level->vector_size - level->inner_vector_size),
                     4 * SIMD_LENGTH_double);
  }
#endif
#ifdef OPTIMIZED_SELF_COUPLING_double
  if (level->type == _FINE) {
    MALLOC_HUGEPAGES(op.clover_vectorized, double, 2 * 6 * level->inner_vector_size,
                     4 * SIMD_LENGTH_double);
#ifdef HAVE_TM1p1
    MALLOC_HUGEPAGES(op.clover_doublet_vectorized, double, 4 * 2 * 6 * level->inner_vector_size,
                     4 * SIMD_LENGTH_double);
#endif
  }
#endif
}

void Schwarz::initialize(void (*evalOp)(vector_double, vector_double, operator_struct<double> *,
                                        Level *),
                         Level *level_in) {

  Global *global = level_in->global;

  // set input parameters
  level = level_in;
  type = global->method;
  problemSize = level->inner_vector_size;
  vectorLength = level->schwarz_vector_size;
  number_of_cycles = level->n_cy;

  int *block_lattice = level->block_lattice;

  // compute block parameters
  number_of_blocks = 1;
  number_of_block_sites = 1;
  for (int mu = 0; mu < 4; mu++) {
    number_of_block_sites *= block_lattice[mu];
    number_of_blocks *= level->local_lattice[mu] / block_lattice[mu];
  }
  block_vector_size = number_of_block_sites * level->num_lattice_site_var;

  direction_length[_T] =
      (block_lattice[_T] - 1) * block_lattice[_Z] * block_lattice[_Y] * block_lattice[_X];
  direction_length[_Z] =
      block_lattice[_T] * (block_lattice[_Z] - 1) * block_lattice[_Y] * block_lattice[_X];
  direction_length[_Y] =
      block_lattice[_T] * block_lattice[_Z] * (block_lattice[_Y] - 1) * block_lattice[_X];
  direction_length[_X] =
      block_lattice[_T] * block_lattice[_Z] * block_lattice[_Y] * (block_lattice[_X] - 1);

  // free memory in case of re-allocation
  delete[] x;
  delete[] b;
  delete[] r;
  delete[] w;
  DELETE2D(V);
  DELETE2D(Z);
  //DELETE2D(H);

  // allocate memory
  Solver::allocateMemory();
  Schwarz::allocateMemory();

  // set function pointers
  applyOperator = evalOp;
  blockOperator =
      (level->type == _FINE) ? block_d_plus_clover_double : coarse_block_operator_double;
  boundaryOperator = (level->type == _FINE) ? &Schwarz::blockBoundaryOperator
                                            : &Schwarz::coarseBlockBoundaryOperator;
  blockSolve = (level->type == _FINE && level->odd_even) ? &Schwarz::blockSolveOddeven
                                                         : &Schwarz::localMinres;

  // define block layout
  Schwarz::defineLayout();

  // setup operators
  switch (level->type) {
  case _FINE:
    Schwarz::setupFineOperator(&global->op_double);
    break;
  case _INTERMEDIATE:
    Schwarz::setupIntermediateOperator();
    break;
  case _COARSEST:
    Schwarz::setupCoarseOperator();
    break;
  default:
    error0("invalid level type!\n");
    break;
  }

  // setup odd even operators
  if (level->odd_even)
    setupOddEvenOperator();

  matrix = &op;
}

Schwarz::~Schwarz() {
  delete[] Dphi;
  delete[] latestIter;
  delete[] latestIter2;
  for (int i = 0; i < number_of_blocks; i++) {
    delete[] block[i].bt;
  }
  delete[] block;
  delete[] index_table[_T];
  if (level->type == _FINE && level->odd_even) {
    delete[] odd_even_index_table[_T];
  }
  if (level->global->method == _SIXTEEN_COLOR) {
    delete[] block_list[0];
    delete[] block_list;
  } else if (level->global->method == _TWO_COLOR) {
    delete[] block_list_length;
    for (int i = 0; i < 8; i++) {
      delete[] block_list[i];
    }
    delete[] block_list;
  }

  if (level->type != _COARSEST) {
    //std::cout << op.D << std::endl;
    operator_double_free(&op, _SCHWARZ, level);
  }
  else {
    //std::cout << op.D << std::endl;
    operator_double_free(&op, _ORDINARY, level);
  }
  if(level->type == _FINE) {
    DELETE2D(odd_even_buffer);
  }
  DELETE2D(local_minres_buffer);
}

void Schwarz::solve() {
  switch (type) {
  case _ADDITIVE:
    Schwarz::additiveSolve();
    break;
  case _TWO_COLOR:
    Schwarz::twoColorSolve();
    break;
  case _SIXTEEN_COLOR:
    Schwarz::sixteenColorSolve();
    break;
  default:
    error0("Invalid Schwarz solver\n");
    break;
  }
}

void Schwarz::additiveSolve() {

  if (has_residual == _NO_RES) {
    vector_double_copy(r, b, 0, number_of_blocks * block_vector_size, level);
    vector_double_define(x, 0, 0, level->schwarz_vector_size);
  } else {
    vector_double_copy(latestIter, x, 0, number_of_blocks * block_vector_size, level);
  }

  for (int k = 0; k < number_of_cycles; k++) {
    if (has_residual == _RES) {
      for (int mu = 0; mu < 4; mu++) {
        ghost_update_double(latestIter, mu, +1, &(op.c), level);
        ghost_update_double(latestIter, mu, -1, &(op.c), level);
      }
    }
    for (int i = 0; i < number_of_blocks; i++) {
      // for all blocks of current color NOT involved in communication
      if (block[i].no_comm) {
        // calculate block residual
        if (has_residual == _RES) {
          if (k == 0) {
            blockOperator(Dphi, latestIter, block[i].start * level->num_lattice_site_var, this,
                          level);
            (this->*boundaryOperator)(Dphi, latestIter, _POSITIVE, i);
            vector_double_minus(r, b, Dphi, block[i].start * level->num_lattice_site_var,
                                block[i].start * level->num_lattice_site_var + block_vector_size,
                                level);
          } else {
            (this->*boundaryOperator)(r, latestIter, _NEGATIVE, i);
          }
        }
        // local minres updates x, r and latest iter
        (this->*blockSolve)(x, r, latestIter2, block[i].start * level->num_lattice_site_var);
      }
    }

    if (has_residual == _RES) {
      for (int mu = 0; mu < 4; mu++) {
        ghost_update_wait_double(latestIter, mu, +1, &(op.c), level);
        ghost_update_wait_double(latestIter, mu, -1, &(op.c), level);
      }
    }

    for (int i = 0; i < number_of_blocks; i++) {
      // for all blocks of current color involved in communication
      if (!block[i].no_comm) {
        // calculate block residual
        if (has_residual == _RES) {
          if (k == 0) {
            blockOperator(Dphi, latestIter, block[i].start * level->num_lattice_site_var, this,
                          level);
            (this->*boundaryOperator)(Dphi, latestIter, _POSITIVE, i);
            vector_double_minus(r, b, Dphi, block[i].start * level->num_lattice_site_var,
                                block[i].start * level->num_lattice_site_var + block_vector_size,
                                level);
          } else {
            (this->*boundaryOperator)(r, latestIter, _NEGATIVE, i);
          }
        }
        // local minres updates x, r and latest iter
        (this->*blockSolve)(x, r, latestIter2, block[i].start * level->num_lattice_site_var);
      }
    }
    has_residual = _RES;
    std::swap(latestIter, latestIter2);
  }

  for (int i = 0; i < number_of_blocks; i++) {
    if (level->relax_fac != 1.0)
      vector_double_scale(x, x, level->relax_fac, block[i].start * level->num_lattice_site_var,
                          block[i].start * level->num_lattice_site_var + block_vector_size, level);
  }

  // TODO: this should be D_phi
  //calculate D * phi with help of the almost computed residual
  /*if (Dphi != nullptr) {
    for (int mu = 0; mu < 4; mu++) {
      ghost_update_double(latestIter, mu, +1, &(op.c), level);
      ghost_update_double(latestIter, mu, -1, &(op.c), level);
    }

    for (int i = 0; i < number_of_blocks; i++) {
      if (block[i].no_comm) {
        (this->*boundaryOperator)(r, latestIter, _NEGATIVE, i);
        vector_double_minus(Dphi, b, r, block[i].start * level->num_lattice_site_var,
                            block[i].start * level->num_lattice_site_var + block_vector_size,
                            level);
        if (level->relax_fac != 1.0)
          vector_double_scale(
              Dphi, Dphi, level->relax_fac, block[i].start * level->num_lattice_site_var,
              block[i].start * level->num_lattice_site_var + block_vector_size, level);
      }
    }

    for (int mu = 0; mu < 4; mu++) {
      ghost_update_wait_double(latestIter, mu, +1, &(op.c), level);
      ghost_update_wait_double(latestIter, mu, -1, &(op.c), level);
    }

    for (int i = 0; i < number_of_blocks; i++) {
      if (!block[i].no_comm) {
        (this->*boundaryOperator)(r, latestIter, _NEGATIVE, i);
        vector_double_minus(Dphi, b, r, block[i].start * level->num_lattice_site_var,
                            block[i].start * level->num_lattice_site_var + block_vector_size,
                            level);
        if (level->relax_fac != 1.0)
          vector_double_scale(
              Dphi, Dphi, level->relax_fac, block[i].start * level->num_lattice_site_var,
              block[i].start * level->num_lattice_site_var + block_vector_size, level);
      }
    }
  }*/

#ifdef SCHWARZ_RES
  /*if (Dphi == nullptr) {
    for (int mu = 0; mu < 4; mu++) {
      ghost_update_double(latestIter, mu, +1, &(op.c), level);
      ghost_update_double(latestIter, mu, -1, &(op.c), level);
    }

    for (int i = 0; i < number_of_blocks; i++) {
      if (block[i].no_comm) {
        (this->*boundaryOperator)(r, latestIter, _NEGATIVE, i);
      }
    }

    for (int mu = 0; mu < 4; mu++) {
      ghost_update_wait_double(latestIter, mu, +1, &(op.c), level);
      ghost_update_wait_double(latestIter, mu, -1, &(op.c), level);
    }

    for (int i = 0; i < number_of_blocks; i++) {
      if (!block[i].no_comm) {
        (this->*boundaryOperator)(r, latestIter, _NEGATIVE, i);
      }
    }
  }*/
  double r_norm = global_norm_double(r, 0, level->inner_vector_size, level);
  char number[3];
  sprintf(number, "%2d", 31 + level->depth);
  printf0("\033[1;%2sm|", number);
  printf0(" ---- depth: %d, c: %d, schwarz iter %2d, norm: %11.6le |", level->depth,
          number_of_colors, k, r_norm);
  printf0("\033[0m\n");
  fflush(0);
#endif
}

void Schwarz::twoColorSolve() {

  initialResidual = has_residual;
  residualCommunication = has_residual;

  communicate[0] = ghost_update_wait_double;
  communicate[1] = ghost_update_double;
  int communicationDirection[8] = {
      +1, -1, -1, +1, -1,
      +1, +1, -1}; // TODO: Move declaration to .h file somehow without needing to set each entry
                   // individually (maybe not possible)

  vector_double D_phi = nullptr; // TODO: check later if/how this is still useful

  start = 0;
  end = level->inner_vector_size;

  if (has_residual == _NO_RES) {
    vector_double_copy(r, b, start, end, level);
    vector_double_define(x, 0, start, end);

    vector_double_define(x, 0, level->inner_vector_size, level->schwarz_vector_size);

  } else {
    for (int mu = 0; mu < 4; mu++)
      ghost_update_double(x, mu, +1, &(op.c), level);
    for (int mu = 0; mu < 4; mu++)
      ghost_update_double(x, mu, -1, &(op.c), level);
  }

  // perform the Schwarz iteration, solve the block systems
  for (k = 0; k < number_of_cycles; k++) {
    for (int step = 0; step < 8; step++) {
      for (int i = 0; i < block_list_length[step]; i++) {
        int index_table = block_list[step][i];
        if (has_residual == _RES) {
          if (k == 0 && initialResidual == _RES) {
            blockOperator(Dphi, x, block[index_table].start * level->num_lattice_site_var, this,
                          level);
            (this->*boundaryOperator)(Dphi, x, _POSITIVE, index_table);
            vector_double_minus(
                r, b, Dphi, block[index_table].start * level->num_lattice_site_var,
                block[index_table].start * level->num_lattice_site_var + block_vector_size, level);
          } else {
            (this->*boundaryOperator)(r, latestIter, _NEGATIVE, index_table);
          }
        }
        // local minres updates x, r and latest iter
        (this->*blockSolve)(x, r, latestIter,
                            block[index_table].start * level->num_lattice_site_var);
      }

      if (residualCommunication == _RES &&
          !(k == number_of_cycles - 1 && (step == 6 || step == 7) && D_phi == nullptr)) {
        for (int mu = 0; mu < 4; mu++) {
          communicate[(step % 4) / 2]((k == 0 && step < 6 && initialResidual == _RES) ? x
                                                                                      : latestIter,

                                      mu, communicationDirection[step], &(op.c), level);
        }
      }

      if (k == 0 && step == 5)
        has_residual = _RES;
      if (k == 0 && step == 1)
        residualCommunication = _RES;
    }
  }

  // copy phi = x
  if (level->relax_fac != 1.0)
    vector_double_scale(x, x, level->relax_fac, start, end, level);

  // calculate D * phi from r
  if (D_phi != nullptr) {
    for (int step = 4; step < 8; step++) {
      for (int i = 0; i < block_list_length[step]; i++) {
        int index_table = block_list[step][i];
        vector_double_minus(
            D_phi, b, r, block[index_table].start * level->num_lattice_site_var,
            block[index_table].start * level->num_lattice_site_var + block_vector_size, level);
        if (level->relax_fac != 1.0) {
          vector_double_scale(
              D_phi, D_phi, level->relax_fac,
              block[index_table].start * level->num_lattice_site_var,
              block[index_table].start * level->num_lattice_site_var + block_vector_size, level);
        }
      }
    }

    for (int step = 0; step < 4; step++) {
      for (int i = 0; i < block_list_length[step]; i++) {
        int index_table = block_list[step][i];
        (this->*boundaryOperator)(r, latestIter, _NEGATIVE, index_table);
        vector_double_minus(
            D_phi, b, r, block[index_table].start * level->num_lattice_site_var,
            block[index_table].start * level->num_lattice_site_var + block_vector_size, level);
        if (level->relax_fac != 1.0) {
          vector_double_scale(
              D_phi, D_phi, level->relax_fac,
              block[index_table].start * level->num_lattice_site_var,
              block[index_table].start * level->num_lattice_site_var + block_vector_size, level);
        }
      }
      if (step == 0 || step == 1) {
        for (int mu = 0; mu < 4; mu++) {
          communicate[0](latestIter, mu, communicationDirection[step], &(op.c), level);
        }
      }
    }
  }

#ifdef SCHWARZ_RES
  int nb = number_of_blocks;
  if (D_phi == nullptr) {
    for (int mu = 0; mu < 4; mu++) {
      ghost_update_double(latestIter, mu, +1, &(op.c), level);
      ghost_update_double(latestIter, mu, -1, &(op.c), level);
    }

    for (int i = 0; i < nb; i++)
      if (block[i].no_comm)
        (this->*boundaryOperator)(r, latestIter, _NEGATIVE, i);

    for (int mu = 0; mu < 4; mu++) {
      ghost_update_wait_double(latestIter, mu, +1, &(op.c), level);
      ghost_update_wait_double(latestIter, mu, -1, &(op.c), level);
    }

    for (int i = 0; i < nb; i++)
      if (!block[i].no_comm)
        (this->*boundaryOperator)(r, latestIter, _NEGATIVE, i);
  }
  double r_norm = global_norm_double(r, 0, level->inner_vector_size, level);
  char number[3];
  sprintf(number, "%2d", 31 + level->depth);
  printf0("\033[1;%2sm|", number);
  printf0(" ---- depth: %d, c: %d, schwarz iter %2d, norm: %11.6le |", level->depth,
          number_of_colors, k, r_norm);
  printf0("\033[0m\n");
  fflush(0);
#endif
}

void Schwarz::sixteenColorSolve() { error0("sixteenColorSolve not yet implemented\n"); }

#include "schwarzUtils.h"
// #include "schwarzUtils.inc"//TODO: make this extension work
