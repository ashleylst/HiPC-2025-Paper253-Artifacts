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

#if (!defined(SSE) || !defined(INTERPOLATION_OPERATOR_LAYOUT_OPTIMIZED_double))

void interpolation_double_alloc(Level *l) {

  int n = l->num_eig_vect;

  l->is_double.op = new complex_double[n * l->inner_vector_size];
  l->is_double.eigenvalues = new complex_double[n];
  l->is_double.test_vector = new complex_double *[n];
  l->is_double.test_vector[0] = new complex_double[n * l->inner_vector_size];
  l->is_double.interpolation = new complex_double *[n];
  l->is_double.interpolation[0] = new complex_double[n * l->vector_size];

  for (int k = 1; k < n; k++) {
    l->is_double.interpolation[k] = l->is_double.interpolation[0] + k * l->vector_size;
    l->is_double.test_vector[k] = l->is_double.test_vector[0] + k * l->inner_vector_size;
  }
}

void interpolation_double_dummy_alloc(Level *l) {

  l->is_double.test_vector = new complex_double *[l->num_eig_vect];
  l->is_double.interpolation = new complex_double *[l->num_eig_vect];
}

void interpolation_double_dummy_free(Level *l) {

  delete[] l->is_double.test_vector;
  delete[] l->is_double.interpolation;
}

void interpolation_double_free(Level *l) {

  delete[] l->is_double.test_vector[0];
  delete[] l->is_double.eigenvalues;
  delete[] l->is_double.test_vector;
  delete[] l->is_double.interpolation[0];
  delete[] l->is_double.interpolation;
  delete[] l->is_double.op;
}

void define_interpolation_double_operator(complex_double **interpolation, Level *l) {

  int j, num_eig_vect = l->num_eig_vect;
  complex_double *op = l->is_double.op;

  int start = 0;
  int end = l->inner_vector_size;

  op += start * num_eig_vect;
  for (int i = start; i < end; i++)
    for (j = 0; j < num_eig_vect; j++) {
      *op = interpolation[j][i];
      op++;
    }
}

void interpolate_double(vector_double phi, vector_double phi_c, Level *l) {

  PROF_double_START(_PR);
  int i, j, k, k1, k2, num_aggregates = l->is_double.num_agg, num_eig_vect = l->num_eig_vect,
                       sign = 1, num_parent_eig_vect = l->num_parent_eig_vect,
                       aggregate_sites = l->num_inner_lattice_sites / num_aggregates;
  complex_double *op = l->is_double.op, *phi_pt = phi,
                 *phi_c_pt = l->next_level->gs_double.transfer_buffer;

  vector_double_distribute(phi_c_pt, phi_c, l->next_level);

#ifdef HAVE_TM1p1
  if (g.n_flavours == 2)
    for (i = 0; i < num_aggregates; i++) {
      phi_pt = phi + i * 2 * 2 * num_parent_eig_vect * aggregate_sites;
      phi_c_pt = l->next_level->gs_double.transfer_buffer + i * 2 * 2 * num_eig_vect;
      op = l->is_double.op + i * 2 * num_eig_vect * num_parent_eig_vect * aggregate_sites;
      for (k = 0; k < aggregate_sites; k++) {
        for (k1 = 0; k1 < 2; k1++) {
          for (k2 = 0; k2 < num_parent_eig_vect; k2++) {
            for (j = 0; j < num_eig_vect; j++) {
              *phi_pt += phi_c_pt[j] * (*op);
              op++;
            }
            phi_pt++;
          }
          op -= num_parent_eig_vect * num_eig_vect;
          phi_c_pt += num_eig_vect;
          for (k2 = 0; k2 < num_parent_eig_vect; k2++) {
            for (j = 0; j < num_eig_vect; j++) {
              *phi_pt += phi_c_pt[j] * (*op);
              op++;
            }
            phi_pt++;
          }
          phi_c_pt -= num_eig_vect;
          phi_c_pt += sign * 2 * num_eig_vect;
          sign *= -1;
        }
      }
    }
  else
#endif
    for (i = 0; i < num_aggregates; i++) {
      phi_pt = phi + i * 2 * num_parent_eig_vect * aggregate_sites;
      phi_c_pt = l->next_level->gs_double.transfer_buffer + i * 2 * num_eig_vect;
      op = l->is_double.op + i * 2 * num_eig_vect * num_parent_eig_vect * aggregate_sites;
      for (k = 0; k < aggregate_sites; k++) {
        for (k1 = 0; k1 < 2; k1++) {
          for (k2 = 0; k2 < num_parent_eig_vect; k2++) {
            for (j = 0; j < num_eig_vect; j++) {
              *phi_pt += phi_c_pt[j] * (*op);
              op++;
            }
            phi_pt++;
          }
          phi_c_pt += sign * num_eig_vect;
          sign *= -1;
        }
      }
    }

  PROF_double_STOP(_PR, 1);
}

void interpolate3_double(vector_double phi, vector_double phi_c, Level *l) {

  PROF_double_START(_PR);
  int i, j, k, k1, k2, num_aggregates = l->is_double.num_agg, num_eig_vect = l->num_eig_vect,
                       num_parent_eig_vect = l->num_parent_eig_vect,
                       aggregate_sites = l->num_inner_lattice_sites / num_aggregates;
  complex_double *op = l->is_double.op, *phi_pt = phi,
                 *phi_c_pt = l->next_level->gs_double.transfer_buffer;

  vector_double_distribute(phi_c_pt, phi_c, l->next_level);

#ifdef HAVE_TM1p1
  if (g.n_flavours == 2)
    for (i = 0; i < num_aggregates; i++) {
      phi_pt = phi + i * 2 * 2 * num_parent_eig_vect * aggregate_sites;
      phi_c_pt = l->next_level->gs_double.transfer_buffer + i * 2 * 2 * num_eig_vect;
      int sign = 1;
      op = l->is_double.op + i * 2 * num_eig_vect * num_parent_eig_vect * aggregate_sites;
      for (k = 0; k < aggregate_sites; k++) {
        for (k1 = 0; k1 < 2; k1++) {
          for (k2 = 0; k2 < num_parent_eig_vect; k2++) {
            *phi_pt = phi_c_pt[0] * (*op);
            op++;
            for (j = 1; j < num_eig_vect; j++) {
              *phi_pt += phi_c_pt[j] * (*op);
              op++;
            }
            phi_pt++;
          }
          op -= num_parent_eig_vect * num_eig_vect;
          phi_c_pt += num_eig_vect;
          for (k2 = 0; k2 < num_parent_eig_vect; k2++) {
            *phi_pt = phi_c_pt[0] * (*op);
            op++;
            for (j = 1; j < num_eig_vect; j++) {
              *phi_pt += phi_c_pt[j] * (*op);
              op++;
            }
            phi_pt++;
          }
          phi_c_pt -= num_eig_vect;
          phi_c_pt += sign * 2 * num_eig_vect;
          sign *= -1;
        }
      }
    }
  else
#endif
    for (i = 0; i < num_aggregates; i++) {
      phi_pt = phi + i * 2 * num_parent_eig_vect * aggregate_sites;
      phi_c_pt = l->next_level->gs_double.transfer_buffer + i * 2 * num_eig_vect;
      int sign = 1;
      op = l->is_double.op + i * 2 * num_eig_vect * num_parent_eig_vect * aggregate_sites;
      for (k = 0; k < aggregate_sites; k++) {
        for (k1 = 0; k1 < 2; k1++) {
          for (k2 = 0; k2 < num_parent_eig_vect; k2++) {
            *phi_pt = phi_c_pt[0] * (*op);
            op++;
            for (j = 1; j < num_eig_vect; j++) {
              *phi_pt += phi_c_pt[j] * (*op);
              op++;
            }
            phi_pt++;
          }
          phi_c_pt += sign * num_eig_vect;
          sign *= -1;
        }
      }
    }
  PROF_double_STOP(_PR, 1);
}

void restrict_double(vector_double phi_c, vector_double phi, Level *l) {

  PROF_double_START(_PR);
  int i, j, k, k1, k2, num_aggregates = l->is_double.num_agg, num_eig_vect = l->num_eig_vect,
                       sign = 1, num_parent_eig_vect = l->num_parent_eig_vect,
                       aggregate_sites = l->num_inner_lattice_sites / num_aggregates;
  complex_double *op = l->is_double.op, *phi_pt = phi,
                 *phi_c_pt = l->next_level->gs_double.transfer_buffer;

#ifdef HAVE_TM1p1
  if (g.n_flavours == 2)
    for (i = 0; i < num_aggregates; i++) {
      phi_pt = phi + i * 2 * 2 * num_parent_eig_vect * aggregate_sites;
      phi_c_pt = l->next_level->gs_double.transfer_buffer + i * 2 * 2 * num_eig_vect;
      op = l->is_double.op + i * 2 * num_eig_vect * num_parent_eig_vect * aggregate_sites;

      for (j = 0; j < 2 * 2 * num_eig_vect; j++)
        phi_c_pt[j] = 0;

      for (k = 0; k < aggregate_sites; k++) {
        for (k1 = 0; k1 < 2; k1++) {
          for (k2 = 0; k2 < num_parent_eig_vect; k2++) {
            for (j = 0; j < num_eig_vect; j++) {
              phi_c_pt[j] += *phi_pt * conj(*op);
              op++;
            }
            phi_pt++;
          }
          op -= num_parent_eig_vect * num_eig_vect;
          phi_c_pt += num_eig_vect;
          for (k2 = 0; k2 < num_parent_eig_vect; k2++) {
            for (j = 0; j < num_eig_vect; j++) {
              phi_c_pt[j] += *phi_pt * conj(*op);
              op++;
            }
            phi_pt++;
          }
          phi_c_pt -= num_eig_vect;
          phi_c_pt += sign * 2 * num_eig_vect;
          sign *= -1;
        }
      }
    }
  else
#endif
    for (i = 0; i < num_aggregates; i++) {
      phi_pt = phi + i * 2 * num_parent_eig_vect * aggregate_sites;
      phi_c_pt = l->next_level->gs_double.transfer_buffer + i * 2 * num_eig_vect;
      op = l->is_double.op + i * 2 * num_eig_vect * num_parent_eig_vect * aggregate_sites;

      for (j = 0; j < 2 * num_eig_vect; j++)
        phi_c_pt[j] = 0;

      for (k = 0; k < aggregate_sites; k++) {
        for (k1 = 0; k1 < 2; k1++) {
          for (k2 = 0; k2 < num_parent_eig_vect; k2++) {
            for (j = 0; j < num_eig_vect; j++) {
              phi_c_pt[j] += *phi_pt * conj(*op);
              op++;
            }
            phi_pt++;
          }
          phi_c_pt += sign * num_eig_vect;
          sign *= -1;
        }
      }
    }

  vector_double_gather(phi_c, l->next_level->gs_double.transfer_buffer, l->next_level);
  PROF_double_STOP(_PR, 1);
}

#endif
