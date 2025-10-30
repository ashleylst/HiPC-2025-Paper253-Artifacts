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

void smoother_double_def(Level *l) {

  //   Global *global = l->global;
  //
  //   if (global->method >= 0)
  //     schwarz_double_def(&(l->s_double), &(global->op_double), l);
  //
  //   l->p_double.op = &(l->s_double.op);
  //   l->p_double.v_start = 0;
  //   l->p_double.v_end = l->inner_vector_size;
  // #ifdef BLOCK_JACOBI
  //   if (l->currentLevel == 0) {
  //     l->p_double.block_jacobi_double.local_p.op = &(l->s_double.op);
  //     l->p_double.block_jacobi_double.local_p.v_start = 0;
  //     l->p_double.block_jacobi_double.local_p.v_end = l->inner_vector_size;
  //   }
  // #endif
  //   if (global->method == 6) {
  //     l->p_double.eval_operator =
  //         (l->depth > 0) ? g5D_apply_coarse_operator_double : g5D_plus_clover_double;
  //   } else {
  //     l->p_double.eval_operator =
  //         (l->depth > 0) ? apply_coarse_operator_double : d_plus_clover_double;
  //   }
}

void smoother_double_free(Level *l) {

  if (l->global->method >= 0)
    schwarz_double_free(&(l->s_double), l);
}

void schwarz_double_init(Schwarz *s, Level *l) {

  operator_double_init(&(s->op));

  s->index_table[_T] = nullptr;
  s->odd_even_index_table[_T] = nullptr;
  s->block = nullptr;
  s->odd_even_buffer = nullptr;
  s->local_minres_buffer = nullptr;
  s->block_list = nullptr;
  s->block_list_length = nullptr;
  s->number_of_colors = 0;
}

void schwarz_double_alloc(Schwarz *s, Level *l) {
  /*
    int i, j, n, mu, nu, *bl = l->block_lattice;
    Global *global = l->global;
    if (global->method == 4) {
      fgmres_double_struct_alloc(
          l->block_iter, 1, (l->depth == 0) ? l->inner_vector_size : l->vector_size, EPS_double,
          _COARSE_GMRES, _NOTHING, nullptr,
          (l->depth == 0) ? (global->odd_even ? apply_schur_complement_double :
  d_plus_clover_double) : (global->odd_even ? coarse_apply_schur_complement_double :
  apply_coarse_operator_double),
          &(l->sp_double), l);
    } else if (global->method == 5) {
      fgmres_double_struct_alloc(
          5, 1, (l->depth == 0) ? l->inner_vector_size : l->vector_size, EPS_double, _COARSE_GMRES,
          _NOTHING, nullptr,
          (l->depth == 0) ? (global->odd_even ? apply_schur_complement_double :
  d_plus_clover_double) : (global->odd_even ? coarse_apply_schur_complement_double :
  apply_coarse_operator_double),
          &(l->sp_double), l);
    } else if (global->method == 6) {
      fgmres_double_struct_alloc(
          l->block_iter, 1, (l->depth == 0) ? l->inner_vector_size : l->vector_size, EPS_double,
          _COARSE_GMRES, _NOTHING, nullptr,
          (l->depth == 0)
              ? (global->odd_even ? g5D_apply_schur_complement_double : g5D_plus_clover_double)
              : (global->odd_even ? g5D_coarse_apply_schur_complement_double
                                  : g5D_apply_coarse_operator_double),
          &(l->sp_double), l);
    }

    operator_double_alloc(&(s->op), _SCHWARZ, l);
    if (l->currentLevel > 0 && l->depth > 0)
      l->p_double.op = &(s->op);

    s->direction_length[_T] = (bl[_T] - 1) * bl[_Z] * bl[_Y] * bl[_X];
    s->direction_length[_Z] = bl[_T] * (bl[_Z] - 1) * bl[_Y] * bl[_X];
    s->direction_length[_Y] = bl[_T] * bl[_Z] * (bl[_Y] - 1) * bl[_X];
    s->direction_length[_X] = bl[_T] * bl[_Z] * bl[_Y] * (bl[_X] - 1);

    MALLOC(s->index_table[_T], int,
           MAX(1, s->direction_length[_T] + s->direction_length[_Z] + s->direction_length[_Y] +
                      s->direction_length[_X]));
    s->index_table[_Z] = s->index_table[_T] + s->direction_length[_T];
    s->index_table[_Y] = s->index_table[_Z] + s->direction_length[_Z];
    s->index_table[_X] = s->index_table[_Y] + s->direction_length[_Y];

    if (l->depth == 0 && global->odd_even) {
      MALLOC(s->odd_even_index_table[_T], int,
             MAX(1, s->direction_length[_T] + s->direction_length[_Z] + s->direction_length[_Y] +
                        s->direction_length[_X]));
      s->odd_even_index_table[_Z] = s->odd_even_index_table[_T] + s->direction_length[_T];
      s->odd_even_index_table[_Y] = s->odd_even_index_table[_Z] + s->direction_length[_Z];
      s->odd_even_index_table[_X] = s->odd_even_index_table[_Y] + s->direction_length[_Y];
    }

    s->number_of_blocks = 1;
    s->number_of_block_sites = 1;
    for (mu = 0; mu < 4; mu++) {
      s->number_of_block_sites *= bl[mu];
      s->number_of_blocks *= l->local_lattice[mu] / bl[mu];
    }
    s->block_vector_size = s->number_of_block_sites * l->num_lattice_site_var;

    if (global->method == 3) {
      MALLOC(s->block_list, int *, 16);
      s->block_list[0] = nullptr;
      MALLOC(s->block_list[0], int, s->number_of_blocks);
      j = s->number_of_blocks / 16;
      for (i = 1; i < 16; i++)
        s->block_list[i] = s->block_list[0] + i * j;
    } else if (global->method == 2) {
      MALLOC(s->block_list_length, int, 8);
      MALLOC(s->block_list, int *, 8);
      for (i = 0; i < 8; i++) {
        s->block_list[i] = nullptr;
        MALLOC(s->block_list[i], int, s->number_of_blocks);
      }
    }

    MALLOC(s->block, block_struct, s->number_of_blocks);

    int svs = l->schwarz_vector_size, vs = (l->depth == 0) ? l->inner_vector_size : l->vector_size;

  #ifdef HAVE_TM1p1
    svs *= 2;
    vs *= 2;
  #endif

    if (l->depth == 0) {
      NEW2D(s->odd_even_buffer, complex_double, vs, 4);
    }

    n = 0;
    for (mu = 0; mu < 4; mu++) {
      i = 1;
      for (nu = 0; nu < 4; nu++) {
        if (mu != nu) {
          i *= bl[nu];
        }
      }
      s->block_boundary_length[2 * mu] = n;
      s->block_boundary_length[2 * mu + 1] = n + 2 * i;
      n += 4 * i;
    }
    s->block_boundary_length[8] = n;

    for (i = 0; i < s->number_of_blocks; i++) {
      s->block[i].bt = nullptr;
      MALLOC(s->block[i].bt, int, n);
    }

    MALLOC(l->sbuf_double[0], complex_double, 2 * vs);
    l->sbuf_double[1] = l->sbuf_double[0] + vs;

    // these buffers are introduced to make local_minres_double thread-safe
    NEW2D(s->local_minres_buffer, complex_double, 3, svs);

  #ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    if (l->depth == 0) {
      MALLOC_HUGEPAGES(s->op.D_vectorized, double,
                       2 * 4 * (2 * l->vector_size - l->inner_vector_size), 4 * SIMD_LENGTH_double);
      MALLOC_HUGEPAGES(s->op.D_transformed_vectorized, double,
                       2 * 4 * (2 * l->vector_size - l->inner_vector_size), 4 * SIMD_LENGTH_double);
    }
  #endif
  #ifdef OPTIMIZED_SELF_COUPLING_double
    if (l->depth == 0) {
      MALLOC_HUGEPAGES(s->op.clover_vectorized, double, 2 * 6 * l->inner_vector_size,
                       4 * SIMD_LENGTH_double);
  #ifdef HAVE_TM1p1
      MALLOC_HUGEPAGES(s->op.clover_doublet_vectorized, double, 4 * 2 * 6 * l->inner_vector_size,
                       4 * SIMD_LENGTH_double);
  #endif
    }
  #endif*/
}

void schwarz_double_free(Schwarz *s, Level *l) {

  int i, n, mu, nu, *bl = l->block_lattice;
  Global *global = l->global;
  if (global->method == 4 || global->method == 5 || global->method == 6)
    //     fgmres_double_struct_free(&(l->sp_double), l);

    FREE(s->index_table[_T], int,
         MAX(1, s->direction_length[_T] + s->direction_length[_Z] + s->direction_length[_Y] +
                    s->direction_length[_X]));
  s->index_table[_Z] = nullptr;
  s->index_table[_Y] = nullptr;
  s->index_table[_X] = nullptr;

  if (l->depth == 0 && global->odd_even) {
    FREE(s->odd_even_index_table[_T], int,
         MAX(1, s->direction_length[_T] + s->direction_length[_Z] + s->direction_length[_Y] +
                    s->direction_length[_X]));
    s->odd_even_index_table[_Z] = nullptr;
    s->odd_even_index_table[_Y] = nullptr;
    s->odd_even_index_table[_X] = nullptr;
  }

  n = 0;
  for (mu = 0; mu < 4; mu++) {
    i = 1;
    for (nu = 0; nu < 4; nu++) {
      if (mu != nu) {
        i *= bl[nu];
      }
    }
    n += 4 * i;
  }

  for (i = 0; i < s->number_of_blocks; i++) {
    FREE(s->block[i].bt, int, n);
  }

  if (global->method == 3) {
    FREE(s->block_list[0], int, s->number_of_blocks);
    FREE(s->block_list, int *, 16);
  } else if (global->method == 2) {
    FREE(s->block_list_length, int, 8);
    for (i = 0; i < 8; i++)
      FREE(s->block_list[i], int, s->number_of_blocks);
    FREE(s->block_list, int *, 8);
  }

  FREE(s->block, block_struct, s->number_of_blocks);

  operator_double_free(&(s->op), _SCHWARZ, l);

  FREE(l->sbuf_double[0], complex_double, 2 * vs);
  l->sbuf_double[1] = nullptr;

  FREE(s->local_minres_buffer[0], complex_double, 0);
  FREE(s->local_minres_buffer[1], complex_double, 0);
  FREE(s->local_minres_buffer[2], complex_double, 0);
  s->local_minres_buffer[0] = nullptr;
  s->local_minres_buffer[1] = nullptr;
  s->local_minres_buffer[2] = nullptr;

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
  if (l->depth == 0) {
    FREE_HUGEPAGES(s->op.D_vectorized, double, 2 * 4 * (2 * l->vector_size - l->inner_vector_size));
    FREE_HUGEPAGES(s->op.D_transformed_vectorized, double,
                   2 * 4 * (2 * l->vector_size - l->inner_vector_size));
  }
#endif
#ifdef OPTIMIZED_SELF_COUPLING_double
  if (l->depth == 0) {
    FREE_HUGEPAGES(s->op.clover_vectorized, double, 2 * 6 * l->inner_vector_size);
#ifdef HAVE_TM1p1
    FREE_HUGEPAGES(s->op.clover_doublet_vectorized, double, 4 * 2 * 6 * l->inner_vector_size);
#endif
  }
#endif
}

void schwarz_layout_double_define(Schwarz *s, Level *l) {

  int a0, b0, c0, d0, a1, b1, c1, d1, block_split[4], block_size[4], agg_split[4], i, j, k, mu,
      index, x, y, z, t, ls[4], le[4], l_st[4], l_en[4],
      *dt = s->op.table_dim, *dt_mod = s->op.table_mod_dim, *it = s->op.index_table, *count[4];
  Global *global = l->global;
  // Define coloring
  if (global->method == 1) // Additive
    s->number_of_colors = 1;
  else if (global->method == 2) // Red-Black
    s->number_of_colors = 2;
  else if (global->method == 3) { // 16 Color
    int flag = 0;
    for (mu = 0; mu < 4; mu++) {
      if ((l->local_lattice[mu] / l->block_lattice[mu]) % 2 == 1)
        flag = 1;
    }
    if (flag == 0)
      s->number_of_colors = 16;
    else {
      s->number_of_colors = 2;
      printf0("depth: %d, switching to red black schwarz as smoother\n", l->depth);
    }
  }

  const int sigma[16] = {0, 1, 3, 2, 6, 4, 5, 7, 15, 14, 12, 13, 9, 11, 10, 8};
  const int color_to_comm[16][2] = {{_T, -1}, {_X, +1}, {_Y, +1}, {_X, -1}, {_Z, +1}, {_Y, -1},
                                    {_X, +1}, {_Y, +1}, {_T, +1}, {_X, -1}, {_Y, -1}, {_X, +1},
                                    {_Z, -1}, {_Y, +1}, {_X, -1}, {_Y, -1}};
  int color_counter[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  s->number_of_block_sites = 1;
  s->block_odd_even_offset = 0;
  s->number_of_aggregates = 1;

  if (global->method == 2) {
    for (i = 0; i < 8; i++)
      s->block_list_length[i] = 0;
  }

  for (mu = 0; mu < 4; mu++) {
    s->number_of_block_sites *= l->block_lattice[mu];
    s->block_odd_even_offset += ((l->local_lattice[mu] / l->block_lattice[mu]) *
                                 (global->my_coords[mu] / l->comm_offset[mu])) %
                                2;
    ls[mu] = 0;
    le[mu] = ls[mu] + l->local_lattice[mu];
    dt[mu] = l->local_lattice[mu] + 2;
    dt_mod[mu] = l->local_lattice[mu] + 2;
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
    agg_split[mu] = l->local_lattice[mu] / l->coarsening[mu];
    block_split[mu] = l->coarsening[mu] / l->block_lattice[mu];
    block_size[mu] = l->block_lattice[mu];
    s->number_of_aggregates *= agg_split[mu];
  }
  s->block_odd_even_offset = s->block_odd_even_offset % 2;
  s->block_vector_size = s->number_of_block_sites * l->num_lattice_site_var;

  i = 0;
  j = 0;
  // inner hyper cuboid
  count[_T] = &d1;
  count[_Z] = &c1;
  count[_Y] = &b1;
  count[_X] = &a1;
  for (d0 = 0; d0 < agg_split[_T]; d0++)
    for (c0 = 0; c0 < agg_split[_Z]; c0++)
      for (b0 = 0; b0 < agg_split[_Y]; b0++)
        for (a0 = 0; a0 < agg_split[_X]; a0++) {

          for (d1 = d0 * block_split[_T]; d1 < (d0 + 1) * block_split[_T]; d1++)
            for (c1 = c0 * block_split[_Z]; c1 < (c0 + 1) * block_split[_Z]; c1++)
              for (b1 = b0 * block_split[_Y]; b1 < (b0 + 1) * block_split[_Y]; b1++)
                for (a1 = a0 * block_split[_X]; a1 < (a0 + 1) * block_split[_X]; a1++) {

                  s->block[j].start = i;
                  s->block[j].no_comm = 1;
                  if (s->number_of_colors == 1) {
                    s->block[j].color = 0;
                  } else if (s->number_of_colors == 2) {
                    s->block[j].color = (d1 + c1 + b1 + a1 + s->block_odd_even_offset) % 2;
                  } else if (s->number_of_colors == 16) {
                    for (k = 0; k < 16; k++)
                      if (sigma[k] == 8 * (d1 % 2) + 4 * (c1 % 2) + 2 * (b1 % 2) + 1 * (a1 % 2)) {
                        s->block[j].color = k;
                        s->block_list[k][color_counter[k]] = j;
                        color_counter[k]++;
                        break;
                      }
                  }

                  if (s->number_of_colors == 1 || s->number_of_colors == 2) {
                    for (mu = 0; mu < 4; mu++) {
                      if (((*count[mu]) == 0) || ((*count[mu] + 1) == le[mu] / block_size[mu]))
                        s->block[j].no_comm = 0;
                    }

                    if (s->number_of_colors == 2) {
                      // calculate boundary correspondence of the block
                      int count_plus = 0, count_minus = 0, count_inner = 0, index;
                      for (mu = 0; mu < 4; mu++) {
                        if ((*count[mu]) == 0)
                          count_minus++;
                        if ((*count[mu] + 1) == le[mu] / block_size[mu])
                          count_plus++;
                        if ((*count[mu]) != 0 && (*count[mu] + 1) != le[mu] / block_size[mu])
                          count_inner++;
                      }

                      if (count_inner == 4) {
                        index = 4 * s->block[j].color;
                      } else if (count_minus == 0) {
                        if (s->block[j].color == 0)
                          index = 1;
                        else
                          index = 7;
                      } else if (count_plus == 0) {
                        if (s->block[j].color == 0)
                          index = 3;
                        else
                          index = 5;
                      } else {
                        index = 2 + 4 * s->block[j].color;
                      }

                      s->block_list[index][s->block_list_length[index]] = j;
                      s->block_list_length[index]++;
                    }

                  } else if (s->number_of_colors == 16) {
                    k = s->block[j].color;
                    if (k == 0) {
                      for (mu = 0; mu < 4; mu++) {
                        if ((*count[mu]) == 0)
                          s->block[j].no_comm = 0;
                      }
                    } else {
                      mu = color_to_comm[k][0];
                      if ((color_to_comm[k][1] == +1 &&
                           (*count[mu] + 1) == le[mu] / block_size[mu]) ||
                          (color_to_comm[k][1] == -1 && (*count[mu]) == 0))
                        s->block[j].no_comm = 0;
                    }
                  }

                  j++;

                  // set up index table
                  if (l->depth == 0 && global->odd_even) {
                    // odd even on the blocks
                    // even sites
                    for (t = d1 * block_size[_T]; t < (d1 + 1) * block_size[_T]; t++)
                      for (z = c1 * block_size[_Z]; z < (c1 + 1) * block_size[_Z]; z++)
                        for (y = b1 * block_size[_Y]; y < (b1 + 1) * block_size[_Y]; y++)
                          for (x = a1 * block_size[_X]; x < (a1 + 1) * block_size[_X]; x++) {
                            if (((t - d1 * block_size[_T]) + (z - c1 * block_size[_Z]) +
                                 (y - b1 * block_size[_Y]) + (x - a1 * block_size[_X])) %
                                    2 ==
                                0) {
                              index = lex_index(t, z, y, x, dt);
                              it[index] = i;
                              i++;
                            }
                          }
                    // odd sites
                    for (t = d1 * block_size[_T]; t < (d1 + 1) * block_size[_T]; t++)
                      for (z = c1 * block_size[_Z]; z < (c1 + 1) * block_size[_Z]; z++)
                        for (y = b1 * block_size[_Y]; y < (b1 + 1) * block_size[_Y]; y++)
                          for (x = a1 * block_size[_X]; x < (a1 + 1) * block_size[_X]; x++) {
                            if (((t - d1 * block_size[_T]) + (z - c1 * block_size[_Z]) +
                                 (y - b1 * block_size[_Y]) + (x - a1 * block_size[_X])) %
                                    2 ==
                                1) {
                              index = lex_index(t, z, y, x, dt);
                              it[index] = i;
                              i++;
                            }
                          }
                  } else {
                    // no odd even
                    for (t = d1 * block_size[_T]; t < (d1 + 1) * block_size[_T]; t++)
                      for (z = c1 * block_size[_Z]; z < (c1 + 1) * block_size[_Z]; z++)
                        for (y = b1 * block_size[_Y]; y < (b1 + 1) * block_size[_Y]; y++)
                          for (x = a1 * block_size[_X]; x < (a1 + 1) * block_size[_X]; x++) {
                            index = lex_index(t, z, y, x, dt);
                            it[index] = i;
                            i++;
                          }
                  }
                }
        }

  // boundaries
  for (mu = 0; mu < 4; mu++) {
    l_st[mu] = le[mu];
    l_en[mu] = le[mu] + 1;
    for (t = l_st[_T]; t < l_en[_T]; t++)
      for (z = l_st[_Z]; z < l_en[_Z]; z++)
        for (y = l_st[_Y]; y < l_en[_Y]; y++)
          for (x = l_st[_X]; x < l_en[_X]; x++) {
            index = lex_index(t, z, y, x, dt);
            it[index] = i;
            i++;
          }
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }

  // define negative boundaries
  for (mu = 0; mu < 4; mu++) {
    l_st[mu] = ls[mu] - 1;
    l_en[mu] = ls[mu];
    for (t = l_st[_T]; t < l_en[_T]; t++)
      for (z = l_st[_Z]; z < l_en[_Z]; z++)
        for (y = l_st[_Y]; y < l_en[_Y]; y++)
          for (x = l_st[_X]; x < l_en[_X]; x++) {
            index = lex_mod_index(t, z, y, x, dt);
            it[index] = i;
            i++;
          }
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }

  i = 0;
  j = 0;
  // block boundary table
  for (d0 = 0; d0 < agg_split[_T]; d0++)
    for (c0 = 0; c0 < agg_split[_Z]; c0++)
      for (b0 = 0; b0 < agg_split[_Y]; b0++)
        for (a0 = 0; a0 < agg_split[_X]; a0++)

          for (d1 = d0 * block_split[_T]; d1 < (d0 + 1) * block_split[_T]; d1++)
            for (c1 = c0 * block_split[_Z]; c1 < (c0 + 1) * block_split[_Z]; c1++)
              for (b1 = b0 * block_split[_Y]; b1 < (b0 + 1) * block_split[_Y]; b1++)
                for (a1 = a0 * block_split[_X]; a1 < (a0 + 1) * block_split[_X]; a1++) {
                  // for all blocks
                  i = 0;
                  int block_start[4], block_end[4], tmp;

                  block_start[_T] = d1 * block_size[_T];
                  block_start[_Z] = c1 * block_size[_Z];
                  block_start[_Y] = b1 * block_size[_Y];
                  block_start[_X] = a1 * block_size[_X];
                  block_end[_T] = (d1 + 1) * block_size[_T];
                  block_end[_Z] = (c1 + 1) * block_size[_Z];
                  block_end[_Y] = (b1 + 1) * block_size[_Y];
                  block_end[_X] = (a1 + 1) * block_size[_X];

                  for (mu = 0; mu < 4; mu++) {
                    tmp = block_start[mu];

                    // minus dir
                    block_start[mu] = block_end[mu] - 1;
                    for (t = block_start[_T]; t < block_end[_T]; t++)
                      for (z = block_start[_Z]; z < block_end[_Z]; z++)
                        for (y = block_start[_Y]; y < block_end[_Y]; y++)
                          for (x = block_start[_X]; x < block_end[_X]; x++) {
                            s->block[j].bt[i] = site_index(t, z, y, x, dt, it);
                            i++;
                            s->block[j].bt[i] =
                                connect_link_double(t, z, y, x, mu, +1, dt, it, s, l);
                            i++;
                          }

                    block_start[mu] = tmp;
                    tmp = block_end[mu];

                    // plus dir
                    block_end[mu] = block_start[mu] + 1;
                    for (t = block_start[_T]; t < block_end[_T]; t++)
                      for (z = block_start[_Z]; z < block_end[_Z]; z++)
                        for (y = block_start[_Y]; y < block_end[_Y]; y++)
                          for (x = block_start[_X]; x < block_end[_X]; x++) {
                            s->block[j].bt[i] = site_index(t, z, y, x, dt, it);
                            i++;
                            s->block[j].bt[i] =
                                connect_link_double(t, z, y, x, mu, -1, dt, it, s, l);
                            i++;
                          }
                    block_end[mu] = tmp;
                  }
                  j++;
                }

  // index table for block dirac operator
  if (l->depth == 0 && global->odd_even) {
    count[_T] = &t;
    count[_Z] = &z;
    count[_Y] = &y;
    count[_X] = &x;
    i = 0;
    j = 0;
    for (t = 0; t < block_size[_T]; t++)
      for (z = 0; z < block_size[_Z]; z++)
        for (y = 0; y < block_size[_Y]; y++)
          for (x = 0; x < block_size[_X]; x++) {
            if ((t + z + y + x) % 2 == 0) {
              i++;
            } else {
              j++;
            }
          }
    s->number_of_even_block_sites = i;
    s->number_of_odd_block_sites = j;

    for (mu = 0; mu < 4; mu++) {
      // even sites, plus dir ( = odd sites, minus dir )
      i = 0;
      j = 0;
      for (t = 0; t < block_size[_T]; t++)
        for (z = 0; z < block_size[_Z]; z++)
          for (y = 0; y < block_size[_Y]; y++)
            for (x = 0; x < block_size[_X]; x++) {
              if ((t + z + y + x) % 2 == 0) {
                if (*(count[mu]) < block_size[mu] - 1) {
                  s->odd_even_index_table[mu][j] = i;
                  j++;
                }
                i++;
              }
            }
      s->direction_length_even[mu] = j;
      // odd sites, plus dir ( = even sites, minus dir )
      j = 0;
      for (t = 0; t < block_size[_T]; t++)
        for (z = 0; z < block_size[_Z]; z++)
          for (y = 0; y < block_size[_Y]; y++)
            for (x = 0; x < block_size[_X]; x++) {
              if ((t + z + y + x) % 2 == 1) {
                if (*(count[mu]) < block_size[mu] - 1) {
                  s->odd_even_index_table[mu][s->direction_length_even[mu] + j] = i;
                  j++;
                }
                i++;
              }
            }
      s->direction_length_odd[mu] = j;
    }
  }

  count[_T] = &t;
  count[_Z] = &z;
  count[_Y] = &y;
  count[_X] = &x;
  for (mu = 0; mu < 4; mu++) {
    j = 0;
    for (t = 0; t < block_size[_T]; t++)
      for (z = 0; z < block_size[_Z]; z++)
        for (y = 0; y < block_size[_Y]; y++)
          for (x = 0; x < block_size[_X]; x++) {
            if (*(count[mu]) < block_size[mu] - 1) {
              s->index_table[mu][j] = site_index(t, z, y, x, dt, it);
              j++;
            }
          }
  }

  // define neighbor table (for the application of the entire operator),
  // negative inner boundary table (for communication),
  // translation table (for translation to lexicographical site ordnering)
  define_nt_bt_tt(s->op.neighbor_table, s->op.backward_neighbor_table, s->op.c.boundary_table,
                  s->op.translation_table, it, dt, l);
}

void schwarz_double_boundary_update(Schwarz *s, Level *l) {

  /*********************************************************************************
   * Updates the current level hopping term in "s->op.D" on the process boundaries
   * in all negative directions. This is necessary for enabling Schwarz to perform
   * local block residual updates on demand.
   *********************************************************************************/

  int i, t, z, y, x, mu, nu, index, *it = s->op.index_table, *dt = s->op.table_dim, ls[4], le[4],
                                    buf_length[4], link_size;
  vector_double buf[4] = {nullptr, nullptr, nullptr, nullptr},
                rbuf[4] = {nullptr, nullptr, nullptr, nullptr};
  config_double D = s->op.D;
  Global *global = l->global;
  for (mu = 0; mu < 4; mu++) {
    ls[mu] = 0;
    le[mu] = l->local_lattice[mu];
    buf_length[mu] = 0;
  }

  if (l->depth == 0)
    link_size = 4 * 9;
  else
    link_size = 4 * SQUARE(l->num_lattice_site_var);

  // allocate buffers
  for (mu = 0; mu < 4; mu++) {
    if (l->global_splitting[mu] > 1) {
      buf_length[mu] = link_size;
      for (nu = 0; nu < 4; nu++) {
        if (nu != mu)
          buf_length[mu] *= le[nu];
      }
      MALLOC(buf[mu], complex_double, buf_length[mu]);
      MALLOC(rbuf[mu], complex_double, buf_length[mu]);
    }
  }

  // post recv for desired directions
  for (mu = 0; mu < 4; mu++) {
    if (l->global_splitting[mu] > 1) {
      MPI_Irecv(rbuf[mu], buf_length[mu], MPI_COMPLEX_double, l->neighbor_rank[2 * mu + 1],
                2 * mu + 1, global->comm_cart, &(s->op.c.rreqs[2 * mu + 1]));
    }
  }

  // buffer data for send and send it
  for (mu = 0; mu < 4; mu++) {
    if (l->global_splitting[mu] > 1) {
      ls[mu] = l->local_lattice[mu] - 1;
      i = 0;
      for (t = ls[_T]; t < le[_T]; t++)
        for (z = ls[_Z]; z < le[_Z]; z++)
          for (y = ls[_Y]; y < le[_Y]; y++)
            for (x = ls[_X]; x < le[_X]; x++) {
              index = site_index(t, z, y, x, dt, it);
              vector_double_copy(buf[mu] + i * link_size, D + index * link_size, 0, link_size, l);
              i++;
            }
      MPI_Isend(buf[mu], buf_length[mu], MPI_COMPLEX_double, l->neighbor_rank[2 * mu], 2 * mu + 1,
                global->comm_cart, &(s->op.c.sreqs[2 * mu + 1]));
      ls[mu] = 0;
    }
  }

  // store links in desired ordering after recv
  for (mu = 0; mu < 4; mu++) {
    if (l->global_splitting[mu] > 1) {
      MPI_Wait(&(s->op.c.rreqs[2 * mu + 1]), MPI_STATUS_IGNORE);
      ls[mu] = -1;
      le[mu] = 0;
      i = 0;
      for (t = ls[_T]; t < le[_T]; t++)
        for (z = ls[_Z]; z < le[_Z]; z++)
          for (y = ls[_Y]; y < le[_Y]; y++)
            for (x = ls[_X]; x < le[_X]; x++) {
              index = site_mod_index(t, z, y, x, dt, it);
              vector_double_copy(D + index * link_size, rbuf[mu] + i * link_size, 0, link_size, l);
              i++;
            }
      ls[mu] = 0;
      le[mu] = l->local_lattice[mu];
    }
  }

  // free buffers
  for (mu = 0; mu < 4; mu++) {
    if (l->global_splitting[mu] > 1) {
      MPI_Wait(&(s->op.c.sreqs[2 * mu + 1]), MPI_STATUS_IGNORE);
      FREE(buf[mu], complex_double, buf_length[mu]);
      FREE(rbuf[mu], complex_double, buf_length[mu]);
    }
  }
}

void schwarz_double_setup(Schwarz *s, operator_struct<double> *op_in, Level *l) {

  /*********************************************************************************
   * Copies the Dirac operator and the clover term from op_in into the Schwarz
   * struct (this function is depth 0 only).
   * - operator_struct<double> *op_in: Input operator.
   *********************************************************************************/

  int i, index, n = l->num_inner_lattice_sites, *tt = s->op.translation_table;
  config_double D_out_pt, clover_out_pt, odd_proj_out_pt;
  config_double D_in_pt = op_in->D, clover_in_pt = op_in->clover, odd_proj_in_pt = op_in->odd_proj;
  Global *global = l->global;
  s->op.m0 = op_in->m0;

  for (i = 0; i < n; i++) {
    index = tt[i];
    D_out_pt = s->op.D + 36 * index;
    FOR36(*D_out_pt = (complex_double)*D_in_pt; D_out_pt++; D_in_pt++;);
  }

  if (global->csw != 0) {
    for (i = 0; i < n; i++) {
      index = tt[i];
      clover_out_pt = s->op.clover + 42 * index;
      FOR42(*clover_out_pt = (complex_double)*clover_in_pt; clover_out_pt++; clover_in_pt++;);
    }
  } else {
    for (i = 0; i < n; i++) {
      index = tt[i];
      clover_out_pt = s->op.clover + 12 * index;
      FOR12(*clover_out_pt = (complex_double)*clover_in_pt; clover_out_pt++; clover_in_pt++;);
    }
  }

  for (i = 0; i < n; i++) {
    index = tt[i];
    odd_proj_out_pt = s->op.odd_proj + 12 * index;
    FOR12(*odd_proj_out_pt = (complex_double)*odd_proj_in_pt; odd_proj_out_pt++; odd_proj_in_pt++;);
  }

#ifdef HAVE_TM
  tm_term_double_setup((double)(global->mu_factor[l->depth] * op_in->mu),
                       (double)(global->mu_factor[l->depth] * op_in->mu_even_shift),
                       (double)(global->mu_factor[l->depth] * op_in->mu_odd_shift), &(s->op), l,
                       no_threading);
#endif

#ifdef HAVE_TM1p1
  epsbar_term_double_setup((double)(global->epsbar_factor[l->depth] * op_in->epsbar),
                           (double)(global->epsbar_factor[l->depth] * op_in->epsbar_ig5_even_shift),
                           (double)(global->epsbar_factor[l->depth] * op_in->epsbar_ig5_odd_shift),
                           &(s->op), l, no_threading);
#endif

  schwarz_double_boundary_update(s, l);

  operator_double_set_couplings(&(s->op), l);

  if (global->odd_even)
    schwarz_double_oddeven_setup(s, l);
}

void trans_double(vector_double out, vector_double in, int *tt, Level *l) {

  int i, index;
  vector_double out_pt = out;
  vector_double in_pt = in;
  //   int start = threading->start_site[l->depth];
  //   int end = threading->end_site[l->depth];
  int start = 0;
  int end = l->inner_vector_size / 12;

  // this function seems to do some data reordering, barriers ensure that everything is in sync

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2)
    for (i = start; i < end; i++) {
      index = tt[i];
      out_pt = out + 24 * index;
      in_pt = in + 24 * i;
      FOR24(*out_pt = (complex_double)*in_pt; out_pt++; in_pt++;)
    }
  else
#endif
    for (i = start; i < end; i++) {
      index = tt[i];
      out_pt = out + 12 * index;
      in_pt = in + 12 * i;
      FOR12(*out_pt = (complex_double)*in_pt; out_pt++; in_pt++;)
    }
}

void trans_back_double(vector_double out, vector_double in, int *tt, Level *l) {

  int i, index;
  vector_double out_pt = out;
  vector_double in_pt = in;
  //   int start = threading->start_site[l->depth];
  //   int end = threading->end_site[l->depth];
  int start = 0;
  int end = l->inner_vector_size / 12;

  // this function seems to do some data reordering, barriers ensure that everything is in sync

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2)
    for (i = start; i < end; i++) {
      index = tt[i];
      in_pt = in + 24 * index;
      out_pt = out + 24 * i;
      FOR24(*out_pt = (complex_double)*in_pt; out_pt++; in_pt++;)
    }
  else
#endif
    for (i = start; i < end; i++) {
      index = tt[i];
      in_pt = in + 12 * index;
      out_pt = out + 12 * i;
      FOR12(*out_pt = (complex_double)*in_pt; out_pt++; in_pt++;)
    }
}

void schwarz_double_def(Schwarz *s, operator_struct<double> *op, Level *l) {

  schwarz_double_alloc(s, l);
  schwarz_layout_double_define(s, l);
  schwarz_double_setup(s, op, l);
}

void schwarz_double_mvm_testfun(Schwarz *s, Level *l) {
  /*


    int mu, i, nb = s->number_of_blocks, svs = l->schwarz_vector_size, ivs = l->inner_vector_size,
    vs = l->vector_size; double diff; vector_double v1 = nullptr, v2 = nullptr, v3 = nullptr; void
    (*block_op)(vector_double eta, vector_double phi, int start, Schwarz *s, Level *l,
    struct Thread *threading) = (l->depth == 0) ? block_d_plus_clover_double :
    coarse_block_operator_double; void (*boundary_op)(vector_double eta, vector_double phi, int k,
    Schwarz *s, Level *l) = (l->depth == 0) ? block_double_boundary_op :
    coarse_block_double_boundary_op; void (*op)(vector_double eta, vector_double phi,
    operator_struct<double> *op, Level *l = (l->depth == 0) ?
    d_plus_clover_double : apply_coarse_operator_double;

    MALLOC(v1, complex_double, svs);
    MALLOC(v2, complex_double, vs);
    MALLOC(v3, complex_double, vs);

    vector_double_define_random(v1, 0, ivs, l);

    op(v3, v1, &(s->op), l, no_threading);

    for (mu = 0; mu < 4; mu++) {
      ghost_update_double(v1, mu, +1, &(s->op.c), l);
      ghost_update_double(v1, mu, -1, &(s->op.c), l);
    }

    for (mu = 0; mu < 4; mu++) {
      ghost_update_wait_double(v1, mu, +1, &(s->op.c), l);
      ghost_update_wait_double(v1, mu, -1, &(s->op.c), l);
    }

    for (i = 0; i < nb; i++) {
      block_op(v2, v1, s->block[i].start * l->num_lattice_site_var, s, l, no_threading);
      boundary_op(v2, v1, i, s, l);
    }

    vector_double_minus(v3, v3, v2, 0, l->inner_vector_size, l);
    diff = global_norm_double(v3, 0, l->inner_vector_size, l, no_threading) /
           global_norm_double(v2, 0, l->inner_vector_size, l, no_threading);

    test0_double("depth: %d, correctness of local residual vector: %le\n", l->depth, diff);

    FREE(v1, complex_double, l->schwarz_vector_size);
    FREE(v2, complex_double, l->vector_size);
    FREE(v3, complex_double, l->vector_size);

    */
}
