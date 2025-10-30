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

void interpolation_double_struct_init(interpolation_double_struct *is) {

  is->agg_index[_T] = nullptr;
  is->agg_boundary_index[_T] = nullptr;
  is->agg_boundary_neighbor[_T] = nullptr;
  is->op = nullptr;
  is->test_vector = nullptr;
  is->interpolation = nullptr;
  is->eigenvalues = nullptr;
  is->tmp = nullptr;
  is->bootstrap_vector = nullptr;
  is->bootstrap_eigenvalues = nullptr;
}

void coarsening_index_table_double_alloc(interpolation_double_struct *is, Level *l) {

  int i, j, mu, t, z, y, x, a0, b0, c0, d0, agg_split[4], agg_size[4], *count[4];

  is->num_agg = 1;
  for (mu = 0; mu < 4; mu++) {
    agg_split[mu] = l->local_lattice[mu] / l->coarsening[mu];
    agg_size[mu] = l->coarsening[mu];
    is->num_agg *= agg_split[mu];
  }

  count[_T] = &t;
  count[_Z] = &z;
  count[_Y] = &y;
  count[_X] = &x;
  for (mu = 0; mu < 4; mu++) {
    i = 0;
    j = 0;
    for (d0 = 0; d0 < agg_split[_T]; d0++)
      for (c0 = 0; c0 < agg_split[_Z]; c0++)
        for (b0 = 0; b0 < agg_split[_Y]; b0++)
          for (a0 = 0; a0 < agg_split[_X]; a0++)

            for (t = d0 * agg_size[_T]; t < (d0 + 1) * agg_size[_T]; t++)
              for (z = c0 * agg_size[_Z]; z < (c0 + 1) * agg_size[_Z]; z++)
                for (y = b0 * agg_size[_Y]; y < (b0 + 1) * agg_size[_Y]; y++)
                  for (x = a0 * agg_size[_X]; x < (a0 + 1) * agg_size[_X]; x++)

                    if ((*(count[mu]) + 1) % agg_size[mu] != 0)
                      i++;
                    else
                      j++;

    // number of lattice sites in local volume
    // that have a neighbor site in mu-direction which belongs to the same aggregate
    is->agg_length[mu] = i;
    // number of lattice sites in local volume
    // that have a neighbor site in mu-direction which belongs to a different aggregate
    is->agg_boundary_length[mu] = j;
  }

  // index table for contributions to the self couplings of the coarse operator
  is->agg_index[_T] =
      new int[is->agg_length[_T] + is->agg_length[_Z] + is->agg_length[_Y] + is->agg_length[_X]];
  // index table for contributions to neighbor couplings of the coarse operator
  is->agg_boundary_index[_T] = new int[is->agg_boundary_length[_T] + is->agg_boundary_length[_Z] +
                                       is->agg_boundary_length[_Y] + is->agg_boundary_length[_X]];
  // corresponging neighbors of the sites in agg_boundary_index
  is->agg_boundary_neighbor[_T] =
      new int[is->agg_boundary_length[_T] + is->agg_boundary_length[_Z] +
              is->agg_boundary_length[_Y] + is->agg_boundary_length[_X]];

  // offsets for the directions
  for (mu = 1; mu < 4; mu++) {
    is->agg_index[mu] = is->agg_index[mu - 1] + is->agg_length[mu - 1];
    is->agg_boundary_index[mu] = is->agg_boundary_index[mu - 1] + is->agg_boundary_length[mu - 1];
    is->agg_boundary_neighbor[mu] =
        is->agg_boundary_neighbor[mu - 1] + is->agg_boundary_length[mu - 1];
  }
}

void coarsening_index_table_double_free(interpolation_double_struct *is, Level *l) {

  int mu;

  delete[] is->agg_index[_T];
  delete[] is->agg_boundary_index[_T];
  delete[] is->agg_boundary_neighbor[_T];

  for (mu = 1; mu < 4; mu++) {
    is->agg_index[mu] = nullptr;
    is->agg_boundary_index[mu] = nullptr;
    is->agg_boundary_neighbor[mu] = nullptr;
  }
}

void coarsening_index_table_double_define(interpolation_double_struct *is, Schwarz *s, Level *l) {

  int i, j, k, mu, t, z, y, x, a0, b0, c0, d0, a1, b1, c1, d1, stride, offset, agg_split[4],
      block_split[4], block_size[4], agg_size[4],
      *index_table = s->op.index_table, *table_dim = s->op.table_dim, *count[4], *index_dir,
      *boundary_index_dir, *boundary_neighbor_index_dir, *neighbor = s->op.neighbor_table;

  for (mu = 0; mu < 4; mu++) {
    agg_split[mu] = l->local_lattice[mu] / l->coarsening[mu];
    agg_size[mu] = l->coarsening[mu];
    block_split[mu] = l->coarsening[mu] / l->block_lattice[mu];
    block_size[mu] = l->block_lattice[mu];
  }

  stride = (l->depth == 0) ? 4 : 5;
  offset = (l->depth == 0) ? 0 : 1;
  count[_T] = &t;
  count[_Z] = &z;
  count[_Y] = &y;
  count[_X] = &x;
  // filling index tables according to the schwarz operator layout
  for (mu = 0; mu < 4; mu++) {

    i = 0;
    j = 0;
    index_dir = is->agg_index[mu];
    boundary_index_dir = is->agg_boundary_index[mu];
    boundary_neighbor_index_dir = is->agg_boundary_neighbor[mu];

    for (d0 = 0; d0 < agg_split[_T]; d0++)
      for (c0 = 0; c0 < agg_split[_Z]; c0++)
        for (b0 = 0; b0 < agg_split[_Y]; b0++)
          for (a0 = 0; a0 < agg_split[_X]; a0++)

            for (d1 = d0 * block_split[_T]; d1 < (d0 + 1) * block_split[_T]; d1++)
              for (c1 = c0 * block_split[_Z]; c1 < (c0 + 1) * block_split[_Z]; c1++)
                for (b1 = b0 * block_split[_Y]; b1 < (b0 + 1) * block_split[_Y]; b1++)
                  for (a1 = a0 * block_split[_X]; a1 < (a0 + 1) * block_split[_X]; a1++)

                    for (t = d1 * block_size[_T]; t < (d1 + 1) * block_size[_T]; t++)
                      for (z = c1 * block_size[_Z]; z < (c1 + 1) * block_size[_Z]; z++)
                        for (y = b1 * block_size[_Y]; y < (b1 + 1) * block_size[_Y]; y++)
                          for (x = a1 * block_size[_X]; x < (a1 + 1) * block_size[_X]; x++)

                            if ((*(count[mu]) + 1) % agg_size[mu] != 0) {
                              index_dir[i] = site_index(t, z, y, x, table_dim, index_table);
                              i++;
                            } else {
                              k = site_index(t, z, y, x, table_dim, index_table);
                              boundary_index_dir[j] = k;
                              boundary_neighbor_index_dir[j] = neighbor[stride * k + mu + offset];
                              j++;
                            }
  }
}
