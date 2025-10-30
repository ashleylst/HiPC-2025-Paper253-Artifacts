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

#ifndef OPTIMIZED_LINALG_double
void process_multi_inner_product_MP(int count, complex_double *results, vector_double *phi,
                                    vector_double psi, int start, int end, Level *l) {

  PROF_double_START(_PIP);
  int i;
  for (int c = 0; c < count; c++)
    results[c] = 0.0;

  for (int c = 0; c < count; c++) {
    for (i = start; i < end;) {
      FOR12(results[c] += (complex_double)conj(phi[c][i]) * psi[i]; i++;)
    }
  }

  PROF_double_STOP(_PIP, (double)(end - start) / (double)l->inner_vector_size);
}
#endif

double global_norm_MP(vector_double x, int start, int end, Level *l) {

  PROF_double_START(_GIP);

  int i;
  double local_alpha = 0, global_alpha = 0;

  for (i = start; i < end;)
    //     FOR12(local_alpha += (complex_double)NORM_SQUARE_double(x[i]); i++;) // Why was a norm
    //     cast to complex?
    FOR12(local_alpha += NORM_SQUARE_double(x[i]); i++;)

  if (l->global->num_processes > 1) {
    PROF_double_START(_ALLR);
    MPI_Allreduce(&local_alpha, &global_alpha, 1, MPI_double, MPI_SUM,
                  (l->depth == 0) ? l->global->comm_cart : l->gs_double.level_comm);
    PROF_double_STOP(_ALLR, 1);
    PROF_double_STOP(_GIP, (double)(end - start) / (double)l->inner_vector_size);
    return sqrt((double)global_alpha);
  } else {
    // all threads need the result of the norm
    PROF_double_STOP(_GIP, (double)(end - start) / (double)l->inner_vector_size);
    return sqrt((double)local_alpha);
  }

  return 0; //NOTE: Should not be reachable, its here to please the compiler gods. (Also at several other locations)
}
