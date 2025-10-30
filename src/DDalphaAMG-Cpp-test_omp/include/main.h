/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern
 * Leder.
 *
 * This file is part of the DDalphaAMG solver library.
 *
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 *
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 *
 */

#include <complex>
#include <random>
#include <valarray>
#include <mpi.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include "gtest/gtest.h"

#ifndef MAIN_HEADER
#define MAIN_HEADER

#define STRINGLENGTH 500
#define _FILE_OFFSET_BITS 64
#define EPS_double 1E-14

// #define HAVE_TM    // flag for enable twisted mass
// #define HAVE_TM1p1 // flag for enable doublet for twisted mass
#undef INIT_ONE_PREC // flag undef for enabling additional features in the lib

#define FOR2(e)                                                                                    \
  { e e }
#define FOR3(e)                                                                                    \
  { e e e }
#define FOR4(e)                                                                                    \
  { e e e e }
#define FOR6(e)                                                                                    \
  { e e e e e e }
#define FOR10(e)                                                                                   \
  { e e e e e e e e e e }
#define FOR12(e)                                                                                   \
  { e e e e e e e e e e e e }
#define FOR20(e)                                                                                   \
  { FOR10(e) FOR10(e) }
#define FOR24(e)                                                                                   \
  { FOR12(e) FOR12(e) }
#define FOR36(e)                                                                                   \
  { FOR12(e) FOR12(e) FOR12(e) }
#define FOR40(e)                                                                                   \
  { FOR20(e) FOR20(e) }
#define FOR42(e)                                                                                   \
  { FOR36(e) FOR6(e) }

#define SQUARE(e) (e) * (e)
#define NORM_SQUARE_double(e) SQUARE(real(e)) + SQUARE(imag(e))
#define CSPLIT(e) real(e), imag(e)

#define MPI_double MPI_DOUBLE
#define MPI_COMPLEX_double MPI_DOUBLE_COMPLEX
#define MPI_COMPLEX_float MPI_COMPLEX
#define I std::complex<double>(0.0, 1.0)

#if defined(GCRODR) || defined(POLYPREC)
#define geev_double LAPACKE_zgeev
#define ggev_double LAPACKE_zggev
#define geqr2_double LAPACKE_zgeqr2
#define ungqr_double LAPACKE_zungqr
#define trtri_double LAPACKE_ztrtri
#define gesv_double LAPACKE_zgesv
#define gels_double LAPACKE_zgels
#endif

#define NEW2D(variable, kind, rows, columns)                                                       \
  do {                                                                                             \
    variable = new kind *[(columns)];                                                              \
    variable[0] = new kind[(rows) * (columns)];                                                    \
    for (int i = 1; i < (columns); i++)                                                            \
      variable[i] = variable[i - 1] + (rows);                                                      \
  } while (0)

#define DELETE2D(variable)                                                                         \
  do {                                                                                             \
    delete[] variable[0];                                                                          \
    delete[] variable;                                                                             \
  } while (0)

#ifdef SSE
#define MALLOC(variable, kind, length)                                                             \
  do {                                                                                             \
    if (variable != NULL) {                                                                        \
      printf0("malloc of \"%s\" failed: pointer is not NULL (%s:%d).\n", #variable, __FILE__,      \
              __LINE__);                                                                           \
    }                                                                                              \
    if ((length) > 0) {                                                                            \
      variable = (kind *)memalign(64, sizeof(kind) * (length));                                    \
    }                                                                                              \
    if (variable == NULL && (length) > 0) {                                                        \
      error0("malloc of \"%s\" failed: no memory allocated (%s:%d).\n", #variable, __FILE__,       \
             __LINE__);                                                                            \
    }                                                                                              \
  } while (0)
#else
#define MALLOC(variable, kind, length)                                                             \
  do {                                                                                             \
    if (variable != NULL) {                                                                        \
      printf0("malloc of \"%s\" failed: pointer is not NULL (%s:%d).\n", #variable, __FILE__,      \
              __LINE__);                                                                           \
    }                                                                                              \
    if ((length) > 0) {                                                                            \
      variable = (kind *)malloc(sizeof(kind) * (length));                                          \
    }                                                                                              \
    if (variable == NULL && (length) > 0) {                                                        \
      error0("malloc of \"%s\" failed: no memory allocated (%s:%d).\n", #variable, __FILE__,       \
             __LINE__);                                                                            \
    }                                                                                              \
  } while (0)
#endif

#define FREE(variable, kind, length)                                                               \
  do {                                                                                             \
    if (variable != NULL) {                                                                        \
      free(variable);                                                                              \
      variable = NULL;                                                                             \
    } else {                                                                                       \
      printf0("multiple free of \"%s\"? pointer is already NULL (%s:%d).\n", #variable, __FILE__,  \
              __LINE__);                                                                           \
    }                                                                                              \
  } while (0)

// allocate and deallocate macros (hugepages, aligned)
#include <fcntl.h>
#include <sys/mman.h>
#define HUGE_PAGE_SIZE (2 * 1024 * 1024)
#define ROUND_UP_TO_FULL_PAGE(x) (((x) + HUGE_PAGE_SIZE - 1) / HUGE_PAGE_SIZE * HUGE_PAGE_SIZE)

#define MALLOC_HUGEPAGES(variable, kind, length, alignment)                                        \
  do {                                                                                             \
    if (variable != NULL) {                                                                        \
      printf0("malloc of \"%s\" failed: pointer is not NULL (%s:%d).\n", #variable, __FILE__,      \
              __LINE__);                                                                           \
    }                                                                                              \
    if ((length) > 0) {                                                                            \
      variable = (kind *)memalign(alignment, sizeof(kind) * ((size_t)length));                     \
    }                                                                                              \
    if (variable == NULL && (length) > 0) {                                                        \
      error0("malloc of \"%s\" failed: no memory allocated (%s:%d).\n", #variable, __FILE__,       \
             __LINE__);                                                                            \
    }                                                                                              \
  } while (0)

#define FREE_HUGEPAGES(addr, kind, length) FREE(addr, kind, length)

#ifdef DEBUG
#define DPRINTF0 printf0
#else
#define DPRINTF0(ARGS, ...)
#endif

#define ASSERT(expression)                                                                         \
  do {                                                                                             \
    if (!(expression)) {                                                                           \
      error0("assertion \"%s\" failed (%s:%d)\n       bad choice of input parameters (please "     \
             "read the user manual in /doc).\n",                                                   \
             #expression, __FILE__, __LINE__);                                                     \
    }                                                                                              \
  } while (0)

#define IMPLIES(A, B) !(A) || (B)
#define XOR(A, B) ((A) && !(B)) || (!(A) && (B))
#define NAND(A, B) !((A) && (B))
#define DIVIDES(A, B) A == 0 || ((double)(B) / (double)(A) - (double)((int)(B) / (int)(A))) == 0
#define ASCENDING(A, B, C) ((A) <= (B)) && ((B) <= (C))
#define MAX(A, B) ((A > B) ? A : B)
#define MIN(A, B) ((A < B) ? A : B)

#ifdef DEBUG
#define DEBUGOUTPUT_ARRAY(A, FORMAT, INDEX)                                                        \
  do {                                                                                             \
    char TMPSTRING[100];                                                                           \
    sprintf(TMPSTRING, "%s[%d] = %s\n", #A, INDEX, FORMAT);                                        \
    printf0(TMPSTRING, A[INDEX]);                                                                  \
  } while (0)
#else
#define DEBUGOUTPUT_ARRAY(A, FORMAT, INDEX)
#endif

#ifdef DEBUG
#define DEBUGOUTPUT(A, FORMAT)                                                                     \
  do {                                                                                             \
    char TMPSTRING[100];                                                                           \
    sprintf(TMPSTRING, "%s = %s\n", #A, FORMAT);                                                   \
    printf0(TMPSTRING, A);                                                                         \
  } while (0)
#else
#define DEBUGOUTPUT(A, FORMAT)
#endif

using namespace std;

class Level;
class Global;

typedef struct block_struct {
  int start, color, no_comm, *bt;
} block_struct;

// #include "threading.h"
#include "vectorization_control.h"
#include "enums.h"
#include "main_pre_def_generic.h"
#include "classQR.h"
#include "classSolver.h"

typedef struct plot_table_line {

  double values[_NUM_OPTB];
  struct plot_table_line *next;

} plot_table_line;

typedef struct var_table_entry {

  void *pt;
  char name[STRINGLENGTH];
  char datatype[20];
  struct var_table_entry *next;

} var_table_entry;

typedef struct var_table {

  int evaluation, multiplicative, shift_update, re_setup, track_error, track_cgn_error,
      average_over;
  char scan_var[STRINGLENGTH];
  double start_val, end_val, step_size, *output_table[6];
  var_table_entry *entry, *iterator;
  plot_table_line *p, *p_end;

} var_table;

typedef struct confbuffer_struct {

  double *data;
  struct confbuffer_struct *next;

} confbuffer_struct;

typedef struct {

  gmres_double_struct dp;

} gmres_MP_struct;

#include "classLevel.h"
#include "classGlobal.h"
#include "classVector.h"

extern int my_rank;
// NOTE: Do these functions really go here?

static inline void printf0(const char *format, ...) {
  if (my_rank == 0) {
    va_list argpt;
    va_start(argpt, format);
    vprintf(format, argpt);
    va_end(argpt);
    fflush(0);
  }
}

static inline void printf1(const char *format, ...) {
  if (my_rank == 1) {
    va_list argpt;
    va_start(argpt, format);
    vprintf(format, argpt);
    va_end(argpt);
    fflush(0);
  }
}

static inline void warning(const char *format, ...) {
  printf("\x1b[31mwarning, rank %d: ", my_rank);
  va_list argpt;
  va_start(argpt, format);
  vprintf(format, argpt);
  va_end(argpt);
  printf("\x1b[0m");
  fflush(0);
}

static inline void warning0(const char *format, ...) {
  if (my_rank == 0) {
    printf("\x1b[31mwarning: ");
    va_list argpt;
    va_start(argpt, format);
    vprintf(format, argpt);
    va_end(argpt);
    printf("\x1b[0m");
    fflush(0);
  }
}

static inline void error0(const char *format, ...) {
  if (my_rank == 0) {
    printf("\x1b[31merror: ");
    va_list argpt;
    va_start(argpt, format);
    vprintf(format, argpt);
    va_end(argpt);
    printf("\x1b[0m");
    fflush(0);
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
}

static inline int shortest_dir(int *data) {
  int min = 0, mu;
  for (mu = 1; mu < 4; mu++)
    if (data[mu] < data[min])
      min = mu;
  return min;
}

static inline int gcd(int a, int b) {
  if (b == 0)
    return a;
  return gcd(b, a % b);
}

static inline int lcm(int a, int b) { return (a * b / gcd(a, b)); }

#endif

// functions
#include "clifford.h"

#ifdef SSE
#include "blas_vectorized.h"
#include "sse_blas_vectorized.h"
#include "sse_coarse_operator_generic.h"
#include "sse_complex_double_intrinsic.h"
#include "sse_interpolation_generic.h"
#include "sse_linalg_generic.h"
#include "vectorization_dirac_generic.h"
#else
// no intrinsics
#include "interpolation_generic.h"
#endif

#include "coarse_oddeven_generic.h"
#include "coarse_operator_generic.h"
#include "coarsening_generic.h"
#include "data_layout.h"
#include "dirac.h"
#include "dirac_generic.h"
#include "dirac_block_vpsite.h"
#include "dirac_block_vpentry.h"
#include "dirac_templated_generic.h"
#include "gathering_generic.h"
#include "ghost.h"
#include "ghost_generic.h"
#include "init.h"
#include "init_generic.h"
#include "io.h"
#include "linalg.h"
#include "linalg_generic.h"
#include "linsolve.h"
#include "linsolve_generic.h"
#include "main_post_def_generic.h"
#include "oddeven_generic.h"
#include "oddeven_generic_template.h"
#include "oddeven_block.h"
#include "operator_generic.h"
#include "preconditioner.h"
#include "schwarz_generic.h"
#include "setup_generic.h"
#include "solver_analysis.h"
#include "classBlockSolver.h"
#include "block_gmres_solve.h"
#include "block.h"
// #include "top_level.h"
#include "var_table.h"
#include "vcycle_generic.h"
#ifdef HAVE_LIME
#include <lime.h>
#include <lime_config.h>
//#include <dcap-overload.h>
#include <lime_defs.h>
#include <lime_header.h>
#include <lime_reader.h>
#include <lime_writer.h>
#endif
#include "lime_io.h"

#ifdef GCRODR
#include "gcrodr_generic.h"
#endif

#ifdef BLOCK_JACOBI
#include "block_jacobi_generic.h"
#include "local_polyprec_generic.h"
#endif

#if defined(GCRODR) || defined(POLYPREC)
#include <lapacke.h>
#ifdef GCRODR
//#include <mkl_scalapack.h>
//#include <mkl_blacs.h>
//#include <mkl_pblas.h>
#endif
#include "lapackwrap_generic.h"
#endif

#ifdef POLYPREC
#include "polyprec_generic.h"
#endif

// for unit-Testing:
#include "unit_testing.h"
extern int test_argc;
extern char **test_argv;
