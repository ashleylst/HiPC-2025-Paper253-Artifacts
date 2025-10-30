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

#ifndef SCHWARZ_double_HEADER
#define SCHWARZ_double_HEADER

void smoother_double_def(Level *l);
void smoother_double_free(Level *l);

void schwarz_double_init(Schwarz *s, Level *l);
void schwarz_double_alloc(Schwarz *s, Level *l);
void schwarz_layout_double_define(Schwarz *s, Level *l);
void schwarz_double_boundary_update(Schwarz *s, Level *l);

void schwarz_double_setup(Schwarz *s, operator_struct<double> *op_in, Level *l);

void schwarz_double_def(Schwarz *s, operator_struct<double> *op, Level *l);
void schwarz_double_free(Schwarz *s, Level *l);

void trans_double(vector_double out, vector_double in, int *tt, Level *l);
void trans_back_double(vector_double out, vector_double in, int *tt, Level *l);

void schwarz_double_mvm_testfun(Schwarz *s, Level *l);

static inline int connect_link_double(int t, int z, int y, int x, int mu, int dir, int *dt, int *it,
                                      Schwarz *s, Level *l) {

  int coord[4];
  coord[_T] = t;
  coord[_Z] = z;
  coord[_Y] = y;
  coord[_X] = x;
  coord[mu] += dir;
  if (l->global_splitting[mu] > 1) {
    return site_mod_index(coord[_T], coord[_Z], coord[_Y], coord[_X], dt, it);
  } else {
    coord[mu] = (coord[mu] + l->local_lattice[mu]) % l->local_lattice[mu];
    return site_index(coord[_T], coord[_Z], coord[_Y], coord[_X], dt, it);
  }
}

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
static inline void set_double_D_vectorized(double *out1, double *out2, complex_double *in) {
  // out1: column major, out2: row major
  for (int i = 0; i < 3; i++) {   // column
    for (int j = 0; j < 3; j++) { // row
      out1[8 * i + j] = real(in[3 * j + i]);
      out1[8 * i + 4 + j] = imag(in[3 * j + i]);
      out2[8 * i + j] = real(in[j + 3 * i]);
      out2[8 * i + 4 + j] = imag(in[j + 3 * i]);
    }
    out1[8 * i + 3] = 0.0;
    out1[8 * i + 7] = 0.0;
    out2[8 * i + 3] = 0.0;
    out2[8 * i + 7] = 0.0;
  }
}
#endif

#endif
