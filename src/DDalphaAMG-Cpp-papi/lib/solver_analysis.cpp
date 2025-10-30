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

#include "main.h"

void test_routine(Level *l) {
//   Global *global = l->global;
//   if (global->method >= 0) {
//
//     global->test = 0;
//     if (l->depth == 0) {
// #ifdef HAVE_TM1p1
//       if (global->n_flavours == 2)
//         printf0("\nRunning tests with D = TM doublet operator:\n");
//       else
// #endif
// #ifdef HAVE_TM
//         printf0("\nRunning tests with D = TM Wilson operator:\n");
// #else
//       printf0("\nRunning tests with D = Wilson operator:\n");
// #endif
//     }
//
//     if (global->mixed_precision) {
//       operator_double_test_routine(&(l->s_double.op), l);
//       if (global->method > 0 && global->method < 4)
//         schwarz_double_mvm_testfun(&(l->s_double), l);
//       if (global->method > 0 && global->method < 4 && global->odd_even)
//         block_oddeven_double_test(l);
//     } else {
//       operator_double_test_routine(&(l->s_double.op), l);
//       if (global->method > 0 && global->method < 4)
//         schwarz_double_mvm_testfun(&(l->s_double), l);
//       if (global->method > 0 && global->method < 4 && global->odd_even)
//         block_oddeven_double_test(l);
//     }
//
//     if (global->interpolation && global->method > 0) {
//       if (global->mixed_precision)
//         coarse_operator_double_test_routine(l);
//       else
//         coarse_operator_double_test_routine(l);
//     }
//
//     if (global->test < 1e-5)
//       printf0("TESTS passed, highest error %e < 1e-5\n", global->test);
//     else
//       warning0("some TESTS not passed, highest error %e > 1e-5\n", global->test);
//     printf0("\n");
//   }
//
// #ifdef HAVE_TM1p1
//   if (global->n_flavours == 1 && (global->epsbar != 0 || global->epsbar_ig5_odd_shift != 0 ||
//                                   global->epsbar_ig5_odd_shift != 0)) {
//
//     if (global->method >= 0) {
//
//       global->test = 0;
//       printf0("Running tests with D = TM doublet operator:\n");
//
//       data_layout_n_flavours(2, l);
//
//       if (global->mixed_precision)
//         two_flavours_test_double(&(l->s_double.op), l);
//       else
//         two_flavours_test_double(&(l->s_double.op), l);
//
//       if (global->mixed_precision) {
//         operator_double_test_routine(&(l->s_double.op), l);
//         if (global->method > 0 && global->method < 4)
//           schwarz_double_mvm_testfun(&(l->s_double), l);
//         if (global->method > 0 && global->method < 4 && global->odd_even)
//           block_oddeven_double_test(l);
//       } else {
//         operator_double_test_routine(&(l->s_double.op), l);
//         if (global->method > 0 && global->method < 4)
//           schwarz_double_mvm_testfun(&(l->s_double), l);
//         if (global->method > 0 && global->method < 4 && global->odd_even)
//           block_oddeven_double_test(l);
//       }
//
//       if (global->interpolation && global->method > 0) {
//         if (global->mixed_precision)
//           coarse_operator_double_test_routine(l);
//         else
//           coarse_operator_double_test_routine(l);
//       }
//
//       if (global->test < 1e-5)
//         printf0("TESTS passed, highest error %e < 1e-5\n", global->test);
//       else
//         warning0("some TESTS not passed, highest error %e > 1e-5\n", global->test);
//       printf0("\n");
//
//       data_layout_n_flavours(1, l);
//     }
//   }
// #endif
}

void prof_init(Level *l) {
  Global *global = l->global;
  if (l->depth == 0) {
    global->coarse_time = 0;
    global->coarse_iter_count = 0;
    global->coarsest_time = 0;
  }
  prof_double_init(l);
  prof_double_init(l);
  if (l->next_level != nullptr)
    prof_init(l->next_level);
}

double prof_print(Level *l) {
  double flop = 0;

#ifdef PROFILING
  Global *global = l->global;
  if (l != nullptr && global->print > 0) {
    if (l->depth == 0)
      printf0("\n+----------------------------------------------------------+\n");
    if (l->depth == 0)
      printf0("| solver profiling                                         |\n");
    printf0("+----------------------------------------------------------+\n");
    printf0("| depth: %3d / level: %3d                time    ( count ) |\n", l->depth,
            l->currentLevel);
    printf0("+----------------------------------------------------------+\n");
    flop += prof_double_print(l);
    flop += prof_double_print(l);
    flop += prof_print(l->next_level);
    if (l->depth == 0) {
      int *ll = l->local_lattice;
      printf0("+----------------------------------------------------------+\n");
      printf0("| flop/lattice site: %9.2le                             |\n", flop);
      printf0("| flop/s/MPIprocess: %9.2le                             |\n",
              (flop / global->total_time) * ll[_T] * ll[_Z] * ll[_Y] * ll[_X]);
      printf0("+----------------------------------------------------------+\n\n");
    }
  }
#endif
  return flop;
}
