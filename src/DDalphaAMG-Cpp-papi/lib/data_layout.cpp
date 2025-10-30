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

void data_layout_init(Level *l) {

  int i, j;

  l->num_inner_lattice_sites = 1;
  for (i = 0; i < 4; i++)
    l->num_inner_lattice_sites *= l->local_lattice[i];
  l->num_lattice_sites = l->num_inner_lattice_sites;
  l->inner_vector_size = l->num_inner_lattice_sites * l->num_lattice_site_var;

  j = l->num_lattice_sites;
  for (i = 0; i < 4; i++)
    l->num_lattice_sites += j / l->local_lattice[i];

  l->vector_size = l->num_lattice_sites * l->num_lattice_site_var;
  l->schwarz_vector_size = 2 * l->vector_size - l->inner_vector_size;
}

void data_layout_n_flavours(int nf, Level *l) {

  ASSERT(nf > 0);
  ASSERT(l->type == _FINE);

#ifdef HAVE_TM1p1
  ASSERT(nf <= 2);

  if (l->global->n_flavours == nf)
    return;
  else {

    l->global->n_flavours = nf;

    Level *l_tmp = l;

    while (1) {
      if (l_tmp->type == _FINE)
        l_tmp->num_lattice_site_var = nf * 12;
      else
        l_tmp->num_lattice_site_var = nf * 2 * l_tmp->num_parent_eig_vect;

      l_tmp->inner_vector_size = l_tmp->num_inner_lattice_sites * l_tmp->num_lattice_site_var;

      l_tmp->vector_size = l_tmp->num_lattice_sites * l_tmp->num_lattice_site_var;
      l_tmp->schwarz_vector_size = 2 * l_tmp->vector_size - l_tmp->inner_vector_size;

      if (l_tmp->type == _FINE) {
        l->global->p.v_end = l_tmp->inner_vector_size;
        l->global->p_MP.dp.v_end = l_tmp->inner_vector_size;
      }

      if (l->global->mixed_precision) {
        l_tmp->s_double.block_vector_size =
            l_tmp->s_double.num_block_sites * l_tmp->num_lattice_site_var;
        l_tmp->p_double.v_end = l_tmp->inner_vector_size;
#ifdef BLOCK_JACOBI
        if (l_tmp->currentLevel == 0)
          l_tmp->p_double.block_jacobi_double.local_p.v_end = l_tmp->inner_vector_size;
#endif
        l_tmp->sp_double.v_end = l_tmp->inner_vector_size;
        l_tmp->dummy_p_double.v_end = l_tmp->inner_vector_size;
        if ((l->global->method >= 4 && l->global->odd_even) ||
            (!l_tmp->idle && l_tmp->currentLevel == 0 && l->global->odd_even)) {
          if (l_tmp->currentLevel == 0) {
            l_tmp->p_double.v_end =
                l_tmp->oe_op_double.num_even_sites * l_tmp->num_lattice_site_var;
#ifdef BLOCK_JACOBI
            if (l_tmp->currentLevel == 0)
              l_tmp->p_double.block_jacobi_double.local_p.v_end =
                  l_tmp->oe_op_double.num_even_sites * l_tmp->num_lattice_site_var;
#endif
          } else {
            l_tmp->sp_double.v_end =
                l_tmp->oe_op_double.num_even_sites * l_tmp->num_lattice_site_var;
          }
        }

      } else {
        l_tmp->s_double.block_vector_size =
            l_tmp->s_double.num_block_sites * l_tmp->num_lattice_site_var;
        l_tmp->p_double.v_end = l_tmp->inner_vector_size;
#ifdef BLOCK_JACOBI
        if (l_tmp->currentLevel == 0)
          l_tmp->p_double.block_jacobi_double.local_p.v_end = l_tmp->inner_vector_size;
#endif
        l_tmp->sp_double.v_end = l_tmp->inner_vector_size;
        l_tmp->dummy_p_double.v_end = l_tmp->inner_vector_size;
        if ((l->global->method >= 4 && l->global->odd_even) ||
            (!l_tmp->idle && l_tmp->currentLevel == 0 && l->global->odd_even)) {
          if (l_tmp->currentLevel == 0) {
            l_tmp->p_double.v_end =
                l_tmp->oe_op_double.num_even_sites * l_tmp->num_lattice_site_var;
#ifdef BLOCK_JACOBI
            if (l_tmp->currentLevel == 0)
              l_tmp->p_double.block_jacobi_double.local_p.v_end =
                  l_tmp->oe_op_double.num_even_sites * l_tmp->num_lattice_site_var;
#endif
          } else {
            l_tmp->sp_double.v_end =
                l_tmp->oe_op_double.num_even_sites * l_tmp->num_lattice_site_var;
          }
        }
      }

      if (l->currentLevel == 0 || l_tmp->next_level == nullptr)
        break;

      l_tmp = l_tmp->next_level;
    }
  }

#else
  ASSERT(nf == 1);
#endif
}

void define_eot(int *eot, int *N, Level *l) {

  int i, mu, t, z, y, x, ls[4], le[4], oe_offset = 0;

  for (mu = 0; mu < 4; mu++)
    oe_offset += (l->local_lattice[mu] * (l->global->my_coords[mu] / l->comm_offset[mu])) % 2;
  oe_offset = oe_offset % 2;

  for (mu = 0; mu < 4; mu++) {
    ls[mu] = 0;
    le[mu] = ls[mu] + l->local_lattice[mu];
  }

  i = 0;
  for (t = 0; t < le[_T]; t++)
    for (z = 0; z < le[_Z]; z++)
      for (y = 0; y < le[_Y]; y++)
        for (x = 0; x < le[_X]; x++)
          if ((t + z + y + x + oe_offset) % 2 == 0) {
            eot[lex_index(t, z, y, x, N)] = i;
            i++;
          }

  for (t = 0; t < le[_T]; t++)
    for (z = 0; z < le[_Z]; z++)
      for (y = 0; y < le[_Y]; y++)
        for (x = 0; x < le[_X]; x++)
          if ((t + z + y + x + oe_offset) % 2 == 1) {
            eot[lex_index(t, z, y, x, N)] = i;
            i++;
          }

  for (mu = 0; mu < 4; mu++) {
    ls[mu] = le[mu];
    le[mu]++;

    for (t = ls[_T]; t < le[_T]; t++)
      for (z = ls[_Z]; z < le[_Z]; z++)
        for (y = ls[_Y]; y < le[_Y]; y++)
          for (x = ls[_X]; x < le[_X]; x++)
            if ((t + z + y + x + oe_offset) % 2 == 0) {
              eot[lex_index(t, z, y, x, N)] = i;
              i++;
            }

    for (t = ls[_T]; t < le[_T]; t++)
      for (z = ls[_Z]; z < le[_Z]; z++)
        for (y = ls[_Y]; y < le[_Y]; y++)
          for (x = ls[_X]; x < le[_X]; x++)
            if ((t + z + y + x + oe_offset) % 2 == 1) {
              eot[lex_index(t, z, y, x, N)] = i;
              i++;
            }

    ls[mu] = 0;
    le[mu]--;
  }
}

void define_eo_bt(int **bt, int *eot, int *n_ebs, int *n_obs, int *n_bs, int *N, Level *l) {

  int i, t, z, y, x, mu, nu, le[4], bs, oe_offset = 0, *bt_mu;

  for (mu = 0; mu < 4; mu++) {
    le[mu] = l->local_lattice[mu];
  }

  for (mu = 0; mu < 4; mu++)
    oe_offset += (l->local_lattice[mu] * (l->global->my_coords[mu] / l->comm_offset[mu])) % 2;
  oe_offset = oe_offset % 2;

  for (mu = 0; mu < 4; mu++) {
    bt_mu = bt[2 * mu];
    bs = 1;
    le[mu] = 1;
    for (nu = 0; nu < 4; nu++)
      bs *= le[nu];

    i = 0;
    for (t = 0; t < le[_T]; t++)
      for (z = 0; z < le[_Z]; z++)
        for (y = 0; y < le[_Y]; y++)
          for (x = 0; x < le[_X]; x++)
            if ((t + z + y + x + oe_offset) % 2 == 0) {
              bt_mu[i] = site_index(t, z, y, x, N, eot);
              i++;
            }
    n_ebs[2 * mu] = i;
    n_ebs[2 * mu + 1] = i;

    for (t = 0; t < le[_T]; t++)
      for (z = 0; z < le[_Z]; z++)
        for (y = 0; y < le[_Y]; y++)
          for (x = 0; x < le[_X]; x++)
            if ((t + z + y + x + oe_offset) % 2 == 1) {
              bt_mu[i] = site_index(t, z, y, x, N, eot);
              i++;
            }

    n_obs[2 * mu] = i - n_ebs[2 * mu];
    n_obs[2 * mu + 1] = i - n_ebs[2 * mu + 1];
    n_bs[2 * mu] = i;
    n_bs[2 * mu + 1] = i;
    le[mu] = l->local_lattice[mu];
  }
}

void define_nt_bt_tt(int *nt, int *backward_nt, int **bt, int *tt, int *it, int *dt, Level *l) {

  /*********************************************************************************
   * Defines neighbor table (for the application of the entire operator), negative
   * inner boundary table (for communication) and translation table (for translation
   * to lexicographical site ordnering).
   * - int *nt: neighbor table
   * - int **bt: boundary table
   * - int *tt: translation table
   * - int *it: index table
   * - int *dt: dimension table
   *********************************************************************************/

  ASSERT(dt != nullptr && it != nullptr);

  int i, mu, pos, t, z, y, x, ls[4], le[4], l_st[4], l_en[4], offset, stride, *bt_mu,
      *gs = l->global_splitting;

  for (mu = 0; mu < 4; mu++) {
    ls[mu] = 0;
    le[mu] = ls[mu] + l->local_lattice[mu];
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }

  for (int k = 0; k < (l->type == _FINE ? 4 : 5) * l->num_inner_lattice_sites; k++)
    nt[k] = -1;

  // define neighbor table
  stride = (l->type == _FINE) ? 4 : 5;
  offset = (l->type == _FINE) ? 0 : 1;
  for (t = ls[_T]; t < le[_T]; t++)
    for (z = ls[_Z]; z < le[_Z]; z++)
      for (y = ls[_Y]; y < le[_Y]; y++)
        for (x = ls[_X]; x < le[_X]; x++) {
          pos = site_index(t, z, y, x, dt, it);
          if (offset)
            nt[stride * pos] = pos;
          nt[stride * pos + offset + _T] =
              site_index((gs[_T] > 1) ? t + 1 : (t + 1) % le[_T], z, y, x, dt, it); // T dir
          nt[stride * pos + offset + _Z] =
              site_index(t, (gs[_Z] > 1) ? z + 1 : (z + 1) % le[_Z], y, x, dt, it); // Z dir
          nt[stride * pos + offset + _Y] =
              site_index(t, z, (gs[_Y] > 1) ? y + 1 : (y + 1) % le[_Y], x, dt, it); // Y dir
          nt[stride * pos + offset + _X] =
              site_index(t, z, y, (gs[_X] > 1) ? x + 1 : (x + 1) % le[_X], dt, it); // X dir
        }

  // define backward neighbor table
  stride = (l->type == _FINE) ? 4 : 5;
  offset = (l->type == _FINE) ? 0 : 1;
  for (t = ls[_T]; t < le[_T]; t++)
    for (z = ls[_Z]; z < le[_Z]; z++)
      for (y = ls[_Y]; y < le[_Y]; y++)
        for (x = ls[_X]; x < le[_X]; x++) {
          pos = site_index(t, z, y, x, dt, it);
          if (offset)
            backward_nt[stride * pos] = pos;
          backward_nt[stride * pos + offset + _T] =
              site_index((gs[_T] > 1) ? (t - 1 + dt[_T]) % dt[_T] : (t - 1 + le[_T]) % le[_T], z, y,
                         x, dt, it); // T dir
          backward_nt[stride * pos + offset + _Z] =
              site_index(t, (gs[_Z] > 1) ? (z - 1 + dt[_Z]) % dt[_Z] : (z - 1 + le[_Z]) % le[_Z], y,
                         x, dt, it); // Z dir
          backward_nt[stride * pos + offset + _Y] =
              site_index(t, z, (gs[_Y] > 1) ? (y - 1 + dt[_Y]) % dt[_Y] : (y - 1 + le[_Y]) % le[_Y],
                         x, dt, it); // Y dir
          backward_nt[stride * pos + offset + _X] = site_index(
              t, z, y, (gs[_X] > 1) ? (x - 1 + dt[_X]) % dt[_X] : (x - 1 + le[_X]) % le[_X], dt,
              it); // X dir
        }

  if (bt != nullptr) {
    for (mu = 0; mu < 4; mu++) {
      // define negative boundary table for communication
      l_en[mu] = l_st[mu] + 1;
      bt_mu = bt[2 * mu + 1];
      i = 0;
      for (t = l_st[_T]; t < l_en[_T]; t++)
        for (z = l_st[_Z]; z < l_en[_Z]; z++)
          for (y = l_st[_Y]; y < l_en[_Y]; y++)
            for (x = l_st[_X]; x < l_en[_X]; x++) {
              bt_mu[i] = site_index(t, z, y, x, dt, it);
              i++;
            }
      l_en[mu] = le[mu];

      // define positive boundary table for communication (if desired)
      if (bt[2 * mu] != bt[2 * mu + 1]) {
        l_st[mu] = le[mu] - 1;
        bt_mu = bt[2 * mu];
        i = 0;
        for (t = l_st[_T]; t < l_en[_T]; t++)
          for (z = l_st[_Z]; z < l_en[_Z]; z++)
            for (y = l_st[_Y]; y < l_en[_Y]; y++)
              for (x = l_st[_X]; x < l_en[_X]; x++) {
                bt_mu[i] = site_index(t, z, y, x, dt, it);
                i++;
              }
        l_st[mu] = ls[mu];
      }
    }
  }

  // define layout translation table
  // for translation to lexicographical site ordering
  //   printf0("\ntranslation table:\n");
  if (tt != nullptr) {
    i = 0;
    for (t = 0; t < l->local_lattice[_T]; t++)
      for (z = 0; z < l->local_lattice[_Z]; z++)
        for (y = 0; y < l->local_lattice[_Y]; y++)
          for (x = 0; x < l->local_lattice[_X]; x++) {
            tt[i] = site_index(t, z, y, x, dt, it);
            i++;
          }
  }
}

void define_nt_bt_tt(int *nt, int *backward_nt, int **bt, int *tt, int *it, int *dt, Global *g) {

  /*********************************************************************************
   * Defines neighbor table (for the application of the entire operator), negative
   * inner boundary table (for communication) and translation table (for translation
   * to lexicographical site ordnering).
   * - int *nt: neighbor table
   * - int **bt: boundary table
   * - int *tt: translation table
   * - int *it: index table
   * - int *dt: dimension table
   *********************************************************************************/

  ASSERT(dt != nullptr && it != nullptr);

  int i, mu, pos, t, z, y, x, ls[4], le[4], l_st[4], l_en[4], offset = 0, stride = 4, *bt_mu,
                                                              *gs = g->process_grid;

  for (mu = 0; mu < 4; mu++) {
    ls[mu] = 0;
    le[mu] = ls[mu] + g->local_lattice[0][mu];
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }

  // define neighbor table
  for (t = ls[_T]; t < le[_T]; t++)
    for (z = ls[_Z]; z < le[_Z]; z++)
      for (y = ls[_Y]; y < le[_Y]; y++)
        for (x = ls[_X]; x < le[_X]; x++) {
          pos = site_index(t, z, y, x, dt, it);
          if (offset)
            nt[stride * pos] = pos;
          nt[stride * pos + offset + _T] =
              site_index((gs[_T] > 1) ? t + 1 : (t + 1) % le[_T], z, y, x, dt, it); // T dir
          nt[stride * pos + offset + _Z] =
              site_index(t, (gs[_Z] > 1) ? z + 1 : (z + 1) % le[_Z], y, x, dt, it); // Z dir
          nt[stride * pos + offset + _Y] =
              site_index(t, z, (gs[_Y] > 1) ? y + 1 : (y + 1) % le[_Y], x, dt, it); // Y dir
          nt[stride * pos + offset + _X] =
              site_index(t, z, y, (gs[_X] > 1) ? x + 1 : (x + 1) % le[_X], dt, it); // X dir
        }

  // define backward neighbor table
  for (t = ls[_T]; t < le[_T]; t++)
    for (z = ls[_Z]; z < le[_Z]; z++)
      for (y = ls[_Y]; y < le[_Y]; y++)
        for (x = ls[_X]; x < le[_X]; x++) {
          pos = site_index(t, z, y, x, dt, it);
          if (offset)
            backward_nt[stride * pos] = pos;
          backward_nt[stride * pos + offset + _T] =
              site_index((gs[_T] > 1) ? (t - 1 + dt[_T]) % dt[_T] : (t - 1 + le[_T]) % le[_T], z, y,
                         x, dt, it); // T dir
          backward_nt[stride * pos + offset + _Z] =
              site_index(t, (gs[_Z] > 1) ? (z - 1 + dt[_Z]) % dt[_Z] : (z - 1 + le[_Z]) % le[_Z], y,
                         x, dt, it); // Z dir
          backward_nt[stride * pos + offset + _Y] =
              site_index(t, z, (gs[_Y] > 1) ? (y - 1 + dt[_Y]) % dt[_Y] : (y - 1 + le[_Y]) % le[_Y],
                         x, dt, it); // Y dir
          backward_nt[stride * pos + offset + _X] = site_index(
              t, z, y, (gs[_X] > 1) ? (x - 1 + dt[_X]) % dt[_X] : (x - 1 + le[_X]) % le[_X], dt,
              it); // X dir
        }

  if (bt != nullptr) {
    for (mu = 0; mu < 4; mu++) {
      // define negative boundary table for communication
      l_en[mu] = l_st[mu] + 1;
      bt_mu = bt[2 * mu + 1];
      i = 0;
      for (t = l_st[_T]; t < l_en[_T]; t++)
        for (z = l_st[_Z]; z < l_en[_Z]; z++)
          for (y = l_st[_Y]; y < l_en[_Y]; y++)
            for (x = l_st[_X]; x < l_en[_X]; x++) {
              bt_mu[i] = site_index(t, z, y, x, dt, it);
              i++;
            }
      l_en[mu] = le[mu];

      // define positive boundary table for communication (if desired)
      if (bt[2 * mu] != bt[2 * mu + 1]) {
        l_st[mu] = le[mu] - 1;
        bt_mu = bt[2 * mu];
        i = 0;
        for (t = l_st[_T]; t < l_en[_T]; t++)
          for (z = l_st[_Z]; z < l_en[_Z]; z++)
            for (y = l_st[_Y]; y < l_en[_Y]; y++)
              for (x = l_st[_X]; x < l_en[_X]; x++) {
                bt_mu[i] = site_index(t, z, y, x, dt, it);
                i++;
              }
        l_st[mu] = ls[mu];
      }
    }
  }

  // define layout translation table
  // for translation to lexicographical site ordering
  //   printf0("\ntranslation table:\n");
  if (tt != nullptr) {
    i = 0;
    for (t = 0; t < g->local_lattice[0][_T]; t++)
      for (z = 0; z < g->local_lattice[0][_Z]; z++)
        for (y = 0; y < g->local_lattice[0][_Y]; y++)
          for (x = 0; x < g->local_lattice[0][_X]; x++) {
            tt[i] = site_index(t, z, y, x, dt, it);
            i++;
          }
  }
}
