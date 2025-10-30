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

void operator_double_init(operator_struct<double> *op) {

  op->index_table = nullptr;
  op->neighbor_table = nullptr;
  op->backward_neighbor_table = nullptr;
  op->translation_table = nullptr;
  op->D = nullptr;
  op->D_vectorized = nullptr;
  op->D_transformed_vectorized = nullptr;
  op->clover = nullptr;
  op->clover_oo_inv = nullptr;
  op->clover_vectorized = nullptr;
  op->clover_oo_inv_vectorized = nullptr;
  op->m0 = 0;
  op->odd_proj = nullptr;
  op->oe_offset = 0;
#ifdef HAVE_TM
  op->mu = 0;
  op->mu_even_shift = 0;
  op->mu_odd_shift = 0;
  op->tm_term = nullptr;
#endif
#ifdef HAVE_TM1p1
  op->epsbar = 0;
  op->epsbar_ig5_even_shift = 0;
  op->epsbar_ig5_odd_shift = 0;
  op->epsbar_term = nullptr;
  op->clover_doublet_oo_inv = nullptr;
  op->clover_doublet_vectorized = nullptr;
  op->clover_doublet_oo_inv_vectorized = nullptr;
#endif

  for (int mu = 0; mu < 4; mu++)
    op->config_boundary_table[mu] = nullptr;

  for (int i = 0; i < 8; i++) {
    op->c.boundary_table[i] = nullptr;
    op->c.buffer[i] = nullptr;
    op->c.in_use[i] = 0;
  }
  op->c.comm = 1;
  op->buffer = nullptr;
}

void operator_double_alloc_projection_buffers(operator_struct<double> *op, Level *l) {

  // when used as preconditioner we usually do not need the projection buffers, unless
  // global->method >= 4: then oddeven_setup_double() is called in init.c, method_setup().
  if (l->depth == 0) {
    int its = (l->num_lattice_site_var / 2) * l->num_lattice_sites;
#ifdef HAVE_TM1p1
    its *= 2;
#endif
  }
}

void operator_double_alloc_projection_buffers(operator_struct<double> *op, Global *g) {

  // when used as preconditioner we usually do not need the projection buffers, unless
  // global->method >= 4: then oddeven_setup_double() is called in init.c, method_setup().
  int its = (g->num_lattice_site_var / 2) * g->num_lattice_sites;
#ifdef HAVE_TM1p1
  its *= 2;
#endif
}

void operator_double_alloc(operator_struct<double> *op, const int type, Level *l) {

  /*********************************************************************************
   * Allocates space for setting up an operator.
   * - operator_struct<double> *op: operator struct for which space is allocated.
   * - const int type: Defines the data layout type of the operator.
   * Possible values are: { _ORDINARY, _SCHWARZ }
   *********************************************************************************/

  int mu, nu, its = 1, its_boundary, nls, clover_site_size, coupling_site_size;
  Global *global = l->global;
  if (l->depth == 0) {
    clover_site_size = 42;
    coupling_site_size = 4 * 9;
  } else {
    clover_site_size = (l->num_lattice_site_var * (l->num_lattice_site_var + 1)) / 2;
    coupling_site_size = 4 * l->num_lattice_site_var * l->num_lattice_site_var;
  }

  if (type == _SCHWARZ) {
    its_boundary = 2;
  } else {
    its_boundary = 1;
  }
  for (mu = 0; mu < 4; mu++) {
    its *= (l->local_lattice[mu] + its_boundary);
  }

  nls = (type == _SCHWARZ) ? (2 * l->num_lattice_sites - l->num_inner_lattice_sites)
                           : l->num_inner_lattice_sites;

  op->D = new complex_double[coupling_site_size * nls];
  op->clover = new complex_double[clover_site_size * l->num_inner_lattice_sites];

  int block_site_size =
      (l->depth == 0) ? 12 : (l->num_lattice_site_var / 2 * (l->num_lattice_site_var / 2 + 1));
  op->odd_proj = new complex_double[block_site_size * l->num_inner_lattice_sites];
#ifdef HAVE_TM
  op->tm_term = new complex_double[block_site_size * l->num_inner_lattice_sites];
#endif
#ifdef HAVE_TM1p1
  op->epsbar_term = new complex_double[block_site_size * l->num_inner_lattice_sites];
#endif

  op->index_table = new int[its];
  if (type == _ODDEVEN) {
    op->neighbor_table = new int[5 * its];
    op->backward_neighbor_table = new int[5 * its];
  } else {
    op->neighbor_table = new int[(l->depth == 0 ? 4 : 5) * l->num_inner_lattice_sites];
    op->backward_neighbor_table = new int[(l->depth == 0 ? 4 : 5) * l->num_inner_lattice_sites];
  }
  op->translation_table = new int[l->num_inner_lattice_sites];

  for (int k = 0; k < its; k++)
    op->index_table[k] = 0;

  if (type == _SCHWARZ && l->depth == 0 && global->odd_even) {
#ifndef OPTIMIZED_SELF_COUPLING_double

    if (global->csw) {
#ifdef HAVE_TM // we use LU here
      op->clover_oo_inv = new complex_double[72 * (l->num_inner_lattice_sites / 2 + 1)];
#else
      op->clover_oo_inv =
          new complex_double[clover_site_size * (l->num_inner_lattice_sites / 2 + 1)];
#endif
    }
#ifdef HAVE_TM1p1
    op->clover_doublet_oo_inv =
        new complex_double[12 * 12 * 2 * (l->num_inner_lattice_sites / 2 + 1)];
#endif

#else
    if (global->csw)
      MALLOC_HUGEPAGES(op->clover_oo_inv_vectorized, double,
                       144 * (l->num_inner_lattice_sites / 2 + 1), 4 * SIMD_LENGTH_double);
#ifdef HAVE_TM1p1
    MALLOC_HUGEPAGES(op->clover_doublet_oo_inv_vectorized, double,
                     2 * 2 * 144 * (l->num_inner_lattice_sites / 2 + 1), 4 * SIMD_LENGTH_double);
#endif

#endif
  }

  if (type != _ODDEVEN)
    operator_double_alloc_projection_buffers(op, l);

  ghost_alloc_double(0, &(op->c), l);

  for (mu = 0; mu < 4; mu++) {
    its = 1;
    for (nu = 0; nu < 4; nu++) {
      if (mu != nu) {
        its *= l->local_lattice[nu];
      }
    }
    op->c.num_boundary_sites[2 * mu] = its;
    op->c.num_boundary_sites[2 * mu + 1] = its;
    op->c.boundary_table[2 * mu] = new int[its];
    if (type == _SCHWARZ) {
      op->c.boundary_table[2 * mu + 1] = new int[its];
      op->config_boundary_table[mu] = new int[its];
    } else {
      op->c.boundary_table[2 * mu + 1] = op->c.boundary_table[2 * mu];
    }
  }
}

void operator_double_alloc(operator_struct<double> *op, const int type, Global *g) {

  /*********************************************************************************
   * Allocates space for setting up an operator.
   * - operator_struct<double> *op: operator struct for which space is allocated.
   * - const int type: Defines the data layout type of the operator.
   * Possible values are: { _ORDINARY, _SCHWARZ, _ODDEVEN }
   *********************************************************************************/

  int mu, nu, its = 1, its_boundary, nls, clover_site_size = 42, coupling_site_size = 4 * 9,
              block_site_size = 12;

  if (type == _SCHWARZ) {
    its_boundary = 2;
  } else {
    its_boundary = 1;
  }
  for (mu = 0; mu < 4; mu++) {
    its *= (g->local_lattice[0][mu] + its_boundary);
  }

  nls = (type == _SCHWARZ) ? (2 * g->num_lattice_sites - g->num_inner_lattice_sites)
                           : g->num_inner_lattice_sites;

  op->D = new complex_double[coupling_site_size * nls];
  op->clover = new complex_double[clover_site_size * g->num_inner_lattice_sites];
  op->odd_proj = new complex_double[block_site_size * g->num_inner_lattice_sites];
#ifdef HAVE_TM
  op->tm_term = new complex_double[block_site_size * g->num_inner_lattice_sites];
#endif
#ifdef HAVE_TM1p1
  op->epsbar_term = new complex_double[block_site_size * g->num_inner_lattice_sites];
#endif

  op->index_table = new int[its];
  if (type == _ODDEVEN) {
    op->neighbor_table = new int[5 * its];
    op->backward_neighbor_table = new int[5 * its];
  } else {
    op->neighbor_table = new int[4 * g->num_inner_lattice_sites];
    op->backward_neighbor_table = new int[4 * g->num_inner_lattice_sites];
  }
  op->translation_table = new int[g->num_inner_lattice_sites];

  if (type == _SCHWARZ && g->odd_even) {
#ifndef OPTIMIZED_SELF_COUPLING_double
    if (g->csw) {
#ifdef HAVE_TM // we use LU here
      op->clover_oo_inv = new complex_double[72 * ((g->num_inner_lattice_sites / 2) + 1)];
#else
      op->clover_oo_inv =
          new complex_double[clover_site_size * ((g->num_inner_lattice_sites / 2) + 1)];
#endif
    }
#ifdef HAVE_TM1p1
    op->clover_doublet_oo_inv =
        new complex_double[12 * 12 * 2 * ((g->num_inner_lattice_sites / 2) + 1)];
#endif

#else
    if (g->csw)
      MALLOC_HUGEPAGES(op->clover_oo_inv_vectorized, double,
                       144 * (g->num_inner_lattice_sites / 2 + 1), 4 * SIMD_LENGTH_double);
#ifdef HAVE_TM1p1
    MALLOC_HUGEPAGES(op->clover_doublet_oo_inv_vectorized, double,
                     2 * 2 * 144 * (g->num_inner_lattice_sites / 2 + 1), 4 * SIMD_LENGTH_double);
#endif

#endif
  }

  if (type != _ODDEVEN)
    operator_double_alloc_projection_buffers(op, g);

  ghost_alloc_double(0, &(op->c), g);

  for (mu = 0; mu < 4; mu++) {
    its = 1;
    for (nu = 0; nu < 4; nu++) {
      if (mu != nu) {
        its *= g->local_lattice[0][nu];
      }
    }
    op->c.num_boundary_sites[2 * mu] = its;
    op->c.num_boundary_sites[2 * mu + 1] = its;
    op->c.boundary_table[2 * mu] = new int[its];
    if (type == _SCHWARZ) {
      op->c.boundary_table[2 * mu + 1] = new int[its];
      op->config_boundary_table[mu] = new int[its];
    } else {
      op->c.boundary_table[2 * mu + 1] = op->c.boundary_table[2 * mu];
    }
  }
}

void operator_double_free(operator_struct<double> *op, const int type, Level *l) {

  int mu;
  Global *global = l->global;

  delete[] op->D;
  delete[] op->clover;
  delete[] op->odd_proj;
#ifdef HAVE_TM
  delete[] op->tm_term;
#endif
  if (type == _SCHWARZ && l->type == _FINE && l->odd_even) {
#ifndef OPTIMIZED_SELF_COUPLING_double
    if (global->csw) {
      delete[] op->clover_oo_inv;
    }
#ifdef HAVE_TM1p1
    delete[] op->clover_doublet_oo_inv;
#endif

#else
    if (global->csw)
      FREE_HUGEPAGES(op->clover_oo_inv_vectorized, double,
                     144 * (l->num_inner_lattice_sites / 2 + 1));
#ifdef HAVE_TM1p1
    FREE_HUGEPAGES(op->clover_doublet_oo_inv_vectorized, double,
                   2 * 2 * 144 * (l->num_inner_lattice_sites / 2 + 1));
#endif

#endif
  }

#ifdef HAVE_TM1p1
  delete[] op->epsbar_term;
#endif
  delete[] op->index_table;
  if (type == _ODDEVEN) {
    delete[] op->neighbor_table;
    delete[] op->backward_neighbor_table;
  } else {
    delete[] op->neighbor_table;
    delete[] op->backward_neighbor_table;
  }
  delete[] op->translation_table;

  ghost_free_double(&(op->c), l);

  for (mu = 0; mu < 4; mu++) {
    delete[] op->c.boundary_table[2 * mu];
    if (type == _SCHWARZ) {
      delete[] op->c.boundary_table[2 * mu + 1];
      delete[] op->config_boundary_table[mu];
    } else {
      op->c.boundary_table[2 * mu + 1] = nullptr;
    }
  }
}

void operator_double_define(operator_struct<double> *op, Level *l) {

  int i, mu, t, z, y, x, *it = op->index_table, ls[4], le[4], l_st[4], l_en[4], *dt = op->table_dim;

  for (mu = 0; mu < 4; mu++) {
    dt[mu] = l->local_lattice[mu] + 1;
    ls[mu] = 0;
    le[mu] = ls[mu] + l->local_lattice[mu];
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }

  // define index table
  // lexicographic inner cuboid and
  // lexicographic +T,+Z,+Y,+X boundaries
  i = 0;
  // inner hyper cuboid
  for (t = ls[_T]; t < le[_T]; t++)
    for (z = ls[_Z]; z < le[_Z]; z++)
      for (y = ls[_Y]; y < le[_Y]; y++)
        for (x = ls[_X]; x < le[_X]; x++) {
          it[lex_index(t, z, y, x, dt)] = i;
          i++;
        }
  // boundaries (buffers)
  for (mu = 0; mu < 4; mu++) {
    l_st[mu] = le[mu];
    l_en[mu] = le[mu] + 1;

    for (t = l_st[_T]; t < l_en[_T]; t++)
      for (z = l_st[_Z]; z < l_en[_Z]; z++)
        for (y = l_st[_Y]; y < l_en[_Y]; y++)
          for (x = l_st[_X]; x < l_en[_X]; x++) {
            it[lex_index(t, z, y, x, dt)] = i;
            i++;
          }

    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }

  // define neighbor table (for the application of the entire operator),
  // negative inner boundary table (for communication),
  // translation table (for translation to lexicographical site ordnering)
  define_nt_bt_tt(op->neighbor_table, op->backward_neighbor_table, op->c.boundary_table,
                  op->translation_table, it, dt, l);
}

void operator_double_define(operator_struct<double> *op, Global *g) {

  int i, mu, t, z, y, x, *it = op->index_table, ls[4], le[4], l_st[4], l_en[4], *dt = op->table_dim;

  for (mu = 0; mu < 4; mu++) {
    dt[mu] = g->local_lattice[0][mu] + 1;
    ls[mu] = 0;
    le[mu] = ls[mu] + g->local_lattice[0][mu];
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }

  // define index table
  // lexicographic inner cuboid and
  // lexicographic +T,+Z,+Y,+X boundaries
  i = 0;
  // inner hyper cuboid
  for (t = ls[_T]; t < le[_T]; t++)
    for (z = ls[_Z]; z < le[_Z]; z++)
      for (y = ls[_Y]; y < le[_Y]; y++)
        for (x = ls[_X]; x < le[_X]; x++) {
          it[lex_index(t, z, y, x, dt)] = i;
          i++;
        }
  // boundaries (buffers)
  for (mu = 0; mu < 4; mu++) {
    l_st[mu] = le[mu];
    l_en[mu] = le[mu] + 1;

    for (t = l_st[_T]; t < l_en[_T]; t++)
      for (z = l_st[_Z]; z < l_en[_Z]; z++)
        for (y = l_st[_Y]; y < l_en[_Y]; y++)
          for (x = l_st[_X]; x < l_en[_X]; x++) {
            it[lex_index(t, z, y, x, dt)] = i;
            i++;
          }

    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }

  // define neighbor table (for the application of the entire operator),
  // negative inner boundary table (for communication),
  // translation table (for translation to lexicographical site ordnering)
  define_nt_bt_tt(op->neighbor_table, op->backward_neighbor_table, op->c.boundary_table,
                  op->translation_table, it, dt, g);
}

void operator_double_set_couplings(operator_struct<double> *op, Level *l) {

  operator_double_set_self_couplings(op, l);
  operator_double_set_neighbor_couplings(op, l);
}

void operator_double_set_neighbor_couplings(operator_struct<double> *op, Level *l) {

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
  int i, n = 2 * l->num_lattice_sites - l->num_inner_lattice_sites;

  for (i = 0; i < n; i++) {
    double *D_vectorized = op->D_vectorized + 96 * i;
    double *D_transformed_vectorized = op->D_transformed_vectorized + 96 * i;
    complex_double *D_pt = op->D + 36 * i;
    for (int mu = 0; mu < 4; mu++)
      set_double_D_vectorized(D_vectorized + 24 * mu, D_transformed_vectorized + 24 * mu,
                              D_pt + 9 * mu);
  }
#endif
}

void operator_double_set_self_couplings(operator_struct<double> *op, Level *l) {

#ifdef OPTIMIZED_SELF_COUPLING_double
  int i, n = l->num_inner_lattice_sites;
  Global *global = l->global;
  if (global->csw != 0)
    for (i = 0; i < n; i++) {
      double *clover_vectorized_pt = op->clover_vectorized + 144 * i;
      config_double clover_pt = op->clover + 42 * i;
      sse_set_clover_double(clover_vectorized_pt, clover_pt);
#ifdef HAVE_TM1p1
      double *clover_doublet_vectorized_pt = op->clover_doublet_vectorized + 288 * i;
      sse_set_clover_doublet_double(clover_doublet_vectorized_pt, clover_pt);
#endif
#ifdef HAVE_TM
      config_double tm_term_pt = op->tm_term + 12 * i;
      sse_add_diagonal_clover_double(clover_vectorized_pt, tm_term_pt);
#ifdef HAVE_TM1p1
      sse_add_diagonal_clover_doublet_double(clover_doublet_vectorized_pt, tm_term_pt);
#endif
#endif
    }
#endif
}
