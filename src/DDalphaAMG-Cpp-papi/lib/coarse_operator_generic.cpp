/*
 * copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern
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

void coarse_operator_double_alloc(Level *l) {

  int nd = l->next_level->num_inner_lattice_sites, k = l->next_level->num_parent_eig_vect * 2;
  l->next_level->D_size = k * k * 4 * nd;
  l->next_level->clover_size = ((k * (k + 1)) / 2) * nd;
  l->next_level->block_size = ((k / 2 * (k / 2 + 1))) * nd;

  operator_double_alloc(&(l->next_level->op_double), _ORDINARY, l->next_level);
}

void coarse_operator_double_free(Level *l) {

  operator_double_free(&(l->next_level->op_double), _ORDINARY, l->next_level);

  coarse_operator_double_free_vectorized(&(l->next_level->s_double.op), l->next_level);
}

void coarse_operator_double_free_vectorized(operator_struct<double> *op, Level *l) {

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
  if (op->D_vectorized != nullptr) {
    int n2 = (l->depth > 0 && l->type != _COARSEST)
                 ? (2 * l->num_lattice_sites - l->num_inner_lattice_sites)
                 : l->num_inner_lattice_sites;
    int column_offset = 2 * SIMD_LENGTH_double *
                        ((l->num_parent_eig_vect + SIMD_LENGTH_double - 1) / SIMD_LENGTH_double);
    // 2 is for complex, 4 is for 4 directions
    FREE_HUGEPAGES(op->D_vectorized, OPERATOR_TYPE_double,
                   2 * 4 * 2 * l->num_parent_eig_vect * column_offset * n2);
    FREE_HUGEPAGES(op->D_transformed_vectorized, OPERATOR_TYPE_double,
                   2 * 4 * 2 * l->num_parent_eig_vect * column_offset * n2);
  }
#endif

#ifdef OPTIMIZED_SELF_COUPLING_double
  if (op->clover_vectorized != nullptr) {
    int n = l->num_inner_lattice_sites;
    int column_offset =
        SIMD_LENGTH_double *
        ((2 * l->num_parent_eig_vect + SIMD_LENGTH_double - 1) / SIMD_LENGTH_double);
    FREE_HUGEPAGES(op->clover_vectorized, OPERATOR_TYPE_double,
                   2 * 2 * l->num_parent_eig_vect * column_offset * n);
#ifdef HAVE_TM1p1
    int column_doublet_offset =
        SIMD_LENGTH_double *
        ((4 * l->num_parent_eig_vect + SIMD_LENGTH_double - 1) / SIMD_LENGTH_double);
    FREE_HUGEPAGES(op->clover_doublet_vectorized, OPERATOR_TYPE_double,
                   2 * 4 * l->num_parent_eig_vect * column_doublet_offset * n);
#endif
  }
#endif
}

void coarse_operator_double_setup(vector_double *V, Level *l) {

  double t0, t1;
  t0 = MPI_Wtime();

  vector_double buffer1 = l->vbuf_double[4], buffer2 = l->vbuf_double[5];

  int mu, n = l->num_eig_vect, i, j, D_size = l->next_level->D_size,
          clover_size = l->next_level->clover_size, block_size = l->next_level->block_size;
  void (*aggregate_self_coupling)(
      vector_double eta1, vector_double eta2, vector_double phi, Schwarz * s,
      Level * l) =
      (l->depth == 0) ? d_plus_clover_aggregate_double : coarse_aggregate_self_couplings_double,
              (*aggregate_neighbor_coupling)(vector_double eta1, vector_double eta2,
                                             vector_double phi, const int mu, Schwarz *s,
                                             Level *l) =
                  (l->depth == 0) ? d_neighbor_aggregate_double
                                  : coarse_aggregate_neighbor_couplings_double;
  void (*aggregate_block)(vector_double eta1, vector_double eta2, vector_double phi,
                          config_double diag, Level * l) =
      (l->depth == 0) ? diagonal_aggregate_double : coarse_aggregate_block_diagonal_double;

  operator_double_define(&(l->next_level->op_double), l->next_level);

  for (j = 0; j < D_size; j++)
    l->next_level->op_double.D[j] = 0;
  for (j = 0; j < clover_size; j++)
    l->next_level->op_double.clover[j] = 0;
  for (j = 0; j < block_size; j++)
    l->next_level->op_double.odd_proj[j] = 0;

  // for all test vectors V[i]:
  for (i = 0; i < n; i++) {
    for (mu = 0; mu < 4; mu++) {
      // update ghost cells of V[i]
      negative_sendrecv_double(V[i], mu, &(l->s_double.op.c), l);
    }
    // apply self coupling of block-and-2spin-restricted dirac operator for each aggregate
    aggregate_self_coupling(buffer1, buffer2, V[i], &(l->s_double), l);
    // calculate selfcoupling entries of the coarse grid operator
    set_coarse_self_coupling_double(buffer1, buffer2, V, i, l);
    // odd_proj
    aggregate_block(buffer1, buffer2, V[i], l->s_double.op.odd_proj, l);
    set_block_diagonal_double(buffer1, buffer2, V, i, l->next_level->op_double.odd_proj, l);

    for (mu = 0; mu < 4; mu++) {
      // finish updating ghostcells of V[i]
      negative_wait_double(mu, &(l->s_double.op.c), l);
      // apply 2spin-restricted dirac operator for direction mu for all aggregates
      aggregate_neighbor_coupling(buffer1, buffer2, V[i], mu, &(l->s_double), l);
      set_coarse_neighbor_coupling_double(buffer1, buffer2, V, mu, i, l);
    }
  }

  coarse_operator_double_setup_finalize(l);

  t1 = MPI_Wtime();
  printf0("depth: %d, time spent for setting up next coarser operator: %lf seconds\n", l->depth,
          t1 - t0);
}

void coarse_operator_double_setup_finalize(Level *l) {

  l->next_level->op_double.m0 = l->s_double.op.m0;
#ifdef HAVE_TM
  // tm_term
  double mf = (l->global->mu_factor[l->depth])
                  ? l->global->mu_factor[l->next_level->depth] / l->global->mu_factor[l->depth]
                  : 0;
  if (mf * l->s_double.op.mu + mf * l->s_double.op.mu_even_shift == 0 &&
      mf * l->s_double.op.mu + mf * l->s_double.op.mu_odd_shift == 0)
    vector_double_define(l->next_level->op_double.tm_term, 0, 0, l->next_level->block_size);
  else
    tm_term_double_setup(mf * l->s_double.op.mu, mf * l->s_double.op.mu_even_shift,
                         mf * l->s_double.op.mu_odd_shift, &(l->next_level->op_double),
                         l->next_level);
#endif
#ifdef HAVE_TM1p1
  // eps_term
  double ef =
      (l->global->epsbar_factor[l->depth])
          ? l->global->epsbar_factor[l->next_level->depth] / l->global->epsbar_factor[l->depth]
          : 0;
  if (ef * l->s_double.op.epsbar == 0 && ef * l->s_double.op.epsbar_ig5_even_shift == 0 &&
      ef * l->s_double.op.epsbar_ig5_odd_shift == 0)
    vector_double_define(l->next_level->op_double.epsbar_term, 0, 0, l->next_level->block_size);
  else
    epsbar_term_double_setup(ef * l->s_double.op.epsbar, ef * l->s_double.op.epsbar_ig5_even_shift,
                             ef * l->s_double.op.epsbar_ig5_odd_shift, &(l->next_level->op_double),
                             l->next_level);
#endif
}

void set_block_diagonal_double(vector_double spin_0_1, vector_double spin_2_3, vector_double *V,
                               const int n, config_double block, Level *l) {

  // U(x) = [ A 0      , A=A*, D=D*
  //          0 D ]
  // storage order: upper triangle of A, upper triangle of D, columnwise
  // suitable for tm_term and odd_proj

  int i, j, k, m, k1, k2,
      num_aggregates = l->is_double.num_agg, num_eig_vect = l->next_level->num_parent_eig_vect,
      aggregate_size = l->num_inner_lattice_sites * l->num_parent_eig_vect * 2 / num_aggregates,
      offset = l->num_parent_eig_vect, block_site_size = (num_eig_vect * (num_eig_vect + 1));
  vector_double spin_0_1_pt, spin_2_3_pt, interpolation_data;
  config_double block_pt;

  for (k = 0; k <= n; k++) {
    k1 = (n * (n + 1)) / 2 + k;
    k2 = (n * (n + 1)) / 2 + k + block_site_size / 2;

    for (j = 0; j < num_aggregates; j++) {
      spin_0_1_pt = spin_0_1 + j * aggregate_size;
      spin_2_3_pt = spin_2_3 + j * aggregate_size;
      interpolation_data = V[k] + j * aggregate_size;
      block_pt = block + j * block_site_size;

      for (i = 0; i < aggregate_size;) {
        // A
        for (m = 0; m < offset; m++, i++)
          block_pt[k1] += conj(interpolation_data[i]) * spin_0_1_pt[i];
        // D
        for (m = 0; m < offset; m++, i++)
          block_pt[k2] += conj(interpolation_data[i]) * spin_2_3_pt[i];
      }
    }
  }
}

void set_coarse_self_coupling_double(vector_double spin_0_1, vector_double spin_2_3,
                                     vector_double *V, const int n, Level *l) {

  int i, j, k, m, k1, k2,
      num_aggregates = l->is_double.num_agg, num_eig_vect = l->next_level->num_parent_eig_vect,
      aggregate_size = l->num_inner_lattice_sites * l->num_parent_eig_vect * 2 / num_aggregates,
      offset = l->num_parent_eig_vect, clover_site_size = (num_eig_vect * (2 * num_eig_vect + 1));
  vector_double spin_0_1_pt, spin_2_3_pt, interpolation_data;
  config_double clover_pt, clover = l->next_level->op_double.clover;

  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  for (k = 0; k <= n; k++) {
    k1 = (n * (n + 1)) / 2 + k;
    k2 = (n * (n + 1)) / 2 + k + (num_eig_vect * (num_eig_vect + 1)) / 2;

    for (j = 0; j < num_aggregates; j++) {
      spin_0_1_pt = spin_0_1 + j * aggregate_size;
      spin_2_3_pt = spin_2_3 + j * aggregate_size;
      interpolation_data = V[k] + j * aggregate_size;
      clover_pt = clover + j * clover_site_size;

      for (i = 0; i < aggregate_size;) {
        // A
        for (m = 0; m < offset; m++, i++)
          clover_pt[k1] += conj(interpolation_data[i]) * spin_0_1_pt[i];
        // D
        for (m = 0; m < offset; m++, i++)
          clover_pt[k2] += conj(interpolation_data[i]) * spin_2_3_pt[i];
      }
    }
  }

  for (k = 0; k < num_eig_vect; k++) {
    k1 = num_eig_vect * (num_eig_vect + 1 + n) + k;

    for (j = 0; j < num_aggregates; j++) {
      spin_0_1_pt = spin_0_1 + j * aggregate_size;
      spin_2_3_pt = spin_2_3 + j * aggregate_size;
      interpolation_data = V[k] + j * aggregate_size;
      clover_pt = clover + j * clover_site_size;

      for (i = 0; i < aggregate_size;) {
        // B
        for (m = 0; m < offset; m++, i++)
          clover_pt[k1] += conj(interpolation_data[i]) * spin_2_3_pt[i];
        i += offset;
      }
    }
  }
}

void set_coarse_neighbor_coupling_double(vector_double spin_0_1, vector_double spin_2_3,
                                         vector_double *V, const int mu, const int n, Level *l) {

  int i, i1, j, k, k1, k2, m,
      num_aggregates = l->is_double.num_agg, num_eig_vect = l->next_level->num_parent_eig_vect,
      offset = l->num_parent_eig_vect, nlsv = l->num_parent_eig_vect * 2,
      D_link_size = num_eig_vect * num_eig_vect * 4,
      *index_dir = l->is_double.agg_boundary_index[mu],
      aggregate_boundary_sites = l->is_double.agg_boundary_length[mu] / num_aggregates;

  vector_double spin_0_1_pt, spin_2_3_pt, interpolation_data;
  config_double D_pt, D = l->next_level->op_double.D;

  // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
  //             C D ]                        -B*  D* ]
  // storage order: A, C, B, D, each column wise
  for (k = 0; k < num_eig_vect; k++) {
    k1 = n * num_eig_vect + k;
    k2 = (n + num_eig_vect) * num_eig_vect + k;
    i1 = 0;
    for (j = 0; j < num_aggregates; j++) {
      D_pt = D + (j * 4 + mu) * D_link_size;

      for (i = 0; i < aggregate_boundary_sites; i++) {
        spin_0_1_pt = spin_0_1 + nlsv * index_dir[i1];
        interpolation_data = V[k] + nlsv * index_dir[i1];
        i1++;
        // A
        for (m = 0; m < offset; m++)
          D_pt[k1] += conj(interpolation_data[m]) * spin_0_1_pt[m];
        // C
        for (; m < 2 * offset; m++)
          D_pt[k2] += conj(interpolation_data[m]) * spin_0_1_pt[m];
      }
    }

    k1 = (n + 2 * num_eig_vect) * num_eig_vect + k;
    k2 = (n + 3 * num_eig_vect) * num_eig_vect + k;
    i1 = 0;
    for (j = 0; j < num_aggregates; j++) {
      D_pt = D + (j * 4 + mu) * D_link_size;

      for (i = 0; i < aggregate_boundary_sites; i++) {
        spin_2_3_pt = spin_2_3 + nlsv * index_dir[i1];
        interpolation_data = V[k] + nlsv * index_dir[i1];
        i1++;
        // B
        for (m = 0; m < offset; m++)
          D_pt[k1] += conj(interpolation_data[m]) * spin_2_3_pt[m];
        // D
        for (; m < 2 * offset; m++)
          D_pt[k2] += conj(interpolation_data[m]) * spin_2_3_pt[m];
      }
    }
  }
}

void coarse_block_operator_double(vector_double eta, vector_double phi, int start, Schwarz *s,
                                  Level *l) {

  int n = s->number_of_block_sites, *length = s->direction_length, **index = s->index_table, *ind,
      *neighbor = s->op.neighbor_table, m = l->num_lattice_site_var,
      num_eig_vect = l->num_parent_eig_vect;
  vector_double lphi = phi + start, leta = eta + start;

  // site-wise self coupling
#ifndef OPTIMIZED_COARSE_SELF_COUPLING_double
  coarse_self_couplings_double(eta, phi, &(s->op), (start / m), (start / m) + n, l);
#else
  coarse_self_couplings_double_vectorized(eta, phi, &(s->op), (start / m), (start / m) + n, l);
#endif

  // inner block couplings
#ifndef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_double
  int hopp_size = 4 * SQUARE(num_eig_vect * 2);
  config_double D_pt, D = s->op.D + (start / m) * hopp_size;

  for (int mu = 0; mu < 4; mu++) {
    ind = index[mu]; // mu direction
    for (int i = 0; i < length[mu]; i++) {
      int k = ind[i];
      int j = neighbor[5 * k + mu + 1];
      D_pt = D + hopp_size * k + (hopp_size / 4) * mu;
      coarse_hopp_double(leta + m * k, lphi + m * j, D_pt, l);
      coarse_daggered_hopp_double(leta + m * j, lphi + m * k, D_pt, l);
    }
  }
#else
  int column_offset =
      2 * SIMD_LENGTH_double * ((num_eig_vect + SIMD_LENGTH_double - 1) / SIMD_LENGTH_double);
  int vectorized_link_offset = 2 * 2 * num_eig_vect * column_offset;
  for (int mu = 0; mu < 4; mu++) {
    OPERATOR_TYPE_double *Dplus =
        s->op.D_vectorized + (start / m) * 4 * vectorized_link_offset + mu * vectorized_link_offset;
    OPERATOR_TYPE_double *Dminus = s->op.D_transformed_vectorized +
                                   (start / m) * 4 * vectorized_link_offset +
                                   mu * vectorized_link_offset;
    ind = index[mu]; // mu direction
    for (int i = 0; i < length[mu]; i++) {
      int k = ind[i];
      int j = neighbor[5 * k + mu + 1];
      // hopp
      coarse_hopp_double_vectorized(leta + m * k, lphi + m * j,
                                    Dplus + 4 * vectorized_link_offset * k, l);
      // daggered hopp
      coarse_hopp_double_vectorized(leta + m * j, lphi + m * k,
                                    Dminus + 4 * vectorized_link_offset * k, l);
    }
  }
#endif
}

void coarse_aggregate_self_couplings_double(vector_double eta1, vector_double eta2,
                                            vector_double phi, Schwarz *s, Level *l) {

  int i, mu, index1, index2, length, *index_dir,
      *neighbor = s->op.neighbor_table, n = l->num_lattice_site_var, Dls = n * n, Dss = 4 * n * n;
  vector_double eta1_pt, eta2_pt, phi_pt;
  config_double D_pt, D = s->op.D;

  vector_double_define(eta1, 0, 0, l->vector_size);
  vector_double_define(eta2, 0, 0, l->vector_size);
  coarse_spinwise_self_couplings_double(eta1, eta2, phi, s->op.clover, l->inner_vector_size, l);

  for (mu = 0; mu < 4; mu++) { // direction mu
    length = l->is_double.agg_length[mu];
    index_dir = l->is_double.agg_index[mu];
    for (i = 0; i < length; i++) {
      index1 = index_dir[i];
      index2 = neighbor[5 * index1 + mu + 1];
      D_pt = D + Dss * index1 + Dls * mu;
      phi_pt = phi + n * index2;
      eta1_pt = eta1 + n * index1;
      eta2_pt = eta2 + n * index1;
      coarse_spinwise_n_hopp_double(eta1_pt, eta2_pt, phi_pt, D_pt, l);
      phi_pt = phi + n * index1;
      eta1_pt = eta1 + n * index2;
      eta2_pt = eta2 + n * index2;
      coarse_spinwise_n_daggered_hopp_double(eta1_pt, eta2_pt, phi_pt, D_pt, l);
    }
  }
}

void coarse_aggregate_neighbor_couplings_double(vector_double eta1, vector_double eta2,
                                                vector_double phi, const int mu, Schwarz *s,
                                                Level *l) {

  int i, index1, index2, length = l->is_double.agg_boundary_length[mu],
                         *index_dir = l->is_double.agg_boundary_index[mu],
                         *neighbor = l->is_double.agg_boundary_neighbor[mu],
                         n = l->num_lattice_site_var, Dls = n * n, Dss = 4 * n * n;
  vector_double eta1_pt, eta2_pt, phi_pt;
  config_double D_pt, D = s->op.D;

  vector_double_define(eta1, 0, 0, l->vector_size);
  vector_double_define(eta2, 0, 0, l->vector_size);

  // requires the positive boundaries of phi to be communicated befor
  for (i = 0; i < length; i++) {
    index1 = index_dir[i];
    index2 = neighbor[i];
    D_pt = D + Dss * index1 + Dls * mu;
    phi_pt = phi + n * index2;
    eta1_pt = eta1 + n * index1;
    eta2_pt = eta2 + n * index1;
    coarse_spinwise_hopp_double(eta1_pt, eta2_pt, phi_pt, D_pt, l);
  }
}

void coarse_self_couplings_double(vector_double eta, vector_double phi, operator_struct<double> *op,
                                  int start, int end, Level *l) {

  int num_eig_vect = l->num_parent_eig_vect, vector_size = l->num_lattice_site_var,
      clover_size = (2 * num_eig_vect * num_eig_vect + num_eig_vect);

  coarse_self_couplings_clover_double(eta + start * vector_size, phi + start * vector_size,
                                      op->clover + start * clover_size, (end - start) * vector_size,
                                      l);
#ifdef HAVE_TM // tm_term
  if (op->mu + op->mu_odd_shift != 0.0 || op->mu + op->mu_even_shift != 0.0)
    coarse_add_anti_block_diagonal_double(eta + start * vector_size, phi + start * vector_size,
                                          op->tm_term +
                                              start * (num_eig_vect * num_eig_vect + num_eig_vect),
                                          (end - start) * vector_size, l);
#endif
#ifdef HAVE_TM1p1 // eps_term
  if (l->global->n_flavours == 2 &&
      (op->epsbar != 0 || op->epsbar_ig5_odd_shift != 0 || op->epsbar_ig5_odd_shift != 0))
    coarse_add_doublet_coupling_double(eta + start * vector_size, phi + start * vector_size,
                                       op->epsbar_term +
                                           start * (num_eig_vect * num_eig_vect + num_eig_vect),
                                       (end - start) * vector_size, l);
#endif
}

void coarse_aggregate_block_diagonal_double(vector_double eta1, vector_double eta2,
                                            vector_double phi, config_double block, Level *l) {
  int length = l->inner_vector_size, num_eig_vect = l->num_parent_eig_vect,
      block_step_size = (num_eig_vect * (num_eig_vect + 1)) / 2;
  config_double block_pt = block;
  vector_double phi_pt = phi, eta1_pt = eta1, eta2_pt = eta2, phi_end_pt = phi + length;
  // U(x) = [ A 0      , A=A*, D=D*
  //          0 D ]
  // storage order: upper triangle of A, upper triangle of D, columnwise
  // diagonal coupling
  while (phi_pt < phi_end_pt) {
    // A
    mvp_double(eta1_pt, block_pt, phi_pt, num_eig_vect);
    vector_double_define(eta2_pt, 0, 0, num_eig_vect);
    block_pt += block_step_size;
    eta1_pt += num_eig_vect;
    eta2_pt += num_eig_vect;
    phi_pt += num_eig_vect;
    // D
    vector_double_define(eta1_pt, 0, 0, num_eig_vect);
    mvp_double(eta2_pt, block_pt, phi_pt, num_eig_vect);
    block_pt += block_step_size;
    eta1_pt += num_eig_vect;
    eta2_pt += num_eig_vect;
    phi_pt += num_eig_vect;
  }
}

void coarse_spinwise_self_couplings_double(vector_double eta1, vector_double eta2,
                                           vector_double phi, config_double clover, int length,
                                           Level *l) {

  int num_eig_vect = l->num_parent_eig_vect,
      clover_step_size1 = (num_eig_vect * (num_eig_vect + 1)) / 2,
      clover_step_size2 = SQUARE(num_eig_vect);
  config_double clover_pt = clover;
  vector_double phi_pt = phi, eta1_pt = eta1, eta2_pt = eta2 + num_eig_vect,
                phi_end_pt = phi + length;
  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  while (phi_pt < phi_end_pt) {
    // A
    mvp_double(eta1_pt, clover_pt, phi_pt, num_eig_vect);
    clover_pt += clover_step_size1;
    phi_pt += num_eig_vect;
    eta1_pt += num_eig_vect;
    // D
    mvp_double(eta2_pt, clover_pt, phi_pt, num_eig_vect);
    clover_pt += clover_step_size1;
    phi_pt -= num_eig_vect;
    eta2_pt -= num_eig_vect;
    // C = -B*
    nmvh_double(eta1_pt, clover_pt, phi_pt, num_eig_vect);
    phi_pt += num_eig_vect;
    eta1_pt += num_eig_vect;
    // B
    mv_double(eta2_pt, clover_pt, phi_pt, num_eig_vect);
    clover_pt += clover_step_size2;
    phi_pt += num_eig_vect;
    eta2_pt += 3 * num_eig_vect;
  }
}

void coarse_operator_double_set_couplings(operator_struct<double> *op, Level *l) {

  coarse_operator_double_set_neighbor_couplings(op, l);
  coarse_operator_double_set_self_couplings(op, l);
}

void coarse_operator_double_set_neighbor_couplings(operator_struct<double> *op, Level *l) {

#ifdef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_double
  int nc_size = SQUARE(l->num_parent_eig_vect * 2);
  int n1, n2;
  int column_offset = 2 * SIMD_LENGTH_double *
                      ((l->num_parent_eig_vect + SIMD_LENGTH_double - 1) / SIMD_LENGTH_double);
  int offset_v = 4 * l->num_parent_eig_vect * column_offset;

  if (l->depth > 0 && l->type != _COARSEST) {
    n1 = l->num_lattice_sites;
    n2 = 2 * l->num_lattice_sites - l->num_inner_lattice_sites;
  } else {
    n1 = l->num_inner_lattice_sites;
    n2 = l->num_inner_lattice_sites;
  }
  int start, end;
  compute_core_start_end_custom(0, n1, &start, &end, l, 1);
  int n_per_core = end - start;

  if (op->D_vectorized == nullptr) {
    // 2 is for complex, 4 is for 4 directions
    MALLOC_HUGEPAGES(op->D_vectorized, OPERATOR_TYPE_double, 4 * offset_v * n2, 64);
    MALLOC_HUGEPAGES(op->D_transformed_vectorized, OPERATOR_TYPE_double, 4 * offset_v * n2, 64);
  }

  copy_coarse_operator_to_vectorized_layout_double(op->D + 4 * start * nc_size,
                                                   op->D_vectorized + 4 * start * offset_v,
                                                   n_per_core, l->num_parent_eig_vect);
  copy_coarse_operator_to_transformed_vectorized_layout_double(
      op->D + 4 * start * nc_size, op->D_transformed_vectorized + 4 * start * offset_v, n_per_core,
      l->num_parent_eig_vect);
  // vectorize negative boundary
  if (n2 > n1) {
    compute_core_start_end_custom(n1, n2, &start, &end, l, 1);
    n_per_core = end - start;
    copy_coarse_operator_to_vectorized_layout_double(op->D + 4 * start * nc_size,
                                                     op->D_vectorized + 4 * start * offset_v,
                                                     n_per_core, l->num_parent_eig_vect);
    copy_coarse_operator_to_transformed_vectorized_layout_double(
        op->D + 4 * start * nc_size, op->D_transformed_vectorized + 4 * start * offset_v,
        n_per_core, l->num_parent_eig_vect);
  }

#endif
}

void coarse_operator_double_set_self_couplings(operator_struct<double> *op, Level *l) {

#ifdef OPTIMIZED_COARSE_SELF_COUPLING_double
  int n = l->num_inner_lattice_sites, nv = l->num_parent_eig_vect;
  int sc_size = (nv) * (nv * 2 + 1);
  int start, end;
  compute_core_start_end_custom(0, n, &start, &end, l, 1);
  int n_per_core = end - start;

  int column_offset = SIMD_LENGTH_double * ((2 * nv + SIMD_LENGTH_double - 1) / SIMD_LENGTH_double);
  int offset_v = 2 * 2 * nv * column_offset;
  if (op->clover_vectorized == nullptr) {

    MALLOC_HUGEPAGES(op->clover_vectorized, OPERATOR_TYPE_double, offset_v * n, 64);
  }
  copy_coarse_operator_clover_to_vectorized_layout_double(
      op->clover + start * sc_size, op->clover_vectorized + start * offset_v, n_per_core, nv);
#ifdef HAVE_TM
  int tm_size = (nv) * (nv + 1);
  if (op->mu + op->mu_odd_shift != 0.0 || op->mu + op->mu_even_shift != 0.0)
    add_tm_term_to_vectorized_layout_double(
        op->tm_term + start * tm_size, op->clover_vectorized + start * offset_v, n_per_core, nv);
#endif

#ifdef HAVE_TM1p1
  int column_doublet_offset =
      SIMD_LENGTH_double * ((4 * nv + SIMD_LENGTH_double - 1) / SIMD_LENGTH_double);
  int offset_doublet_v = 2 * 4 * nv * column_doublet_offset;
  int eps_size = (nv) * (nv + 1);
  if (op->clover_doublet_vectorized == nullptr) {

    MALLOC_HUGEPAGES(op->clover_doublet_vectorized, OPERATOR_TYPE_double, offset_doublet_v * n, 64);
  }
  copy_coarse_operator_clover_to_doublet_vectorized_layout_double(
      op->clover + start * sc_size, op->clover_doublet_vectorized + start * offset_doublet_v,
      n_per_core, nv);
  if (op->epsbar != 0 || op->epsbar_ig5_odd_shift != 0 || op->epsbar_ig5_odd_shift != 0)
    add_epsbar_term_to_doublet_vectorized_layout_double(
        op->epsbar_term + start * eps_size,
        op->clover_doublet_vectorized + start * offset_doublet_v, n_per_core, nv);
#ifdef HAVE_TM
  if (op->mu + op->mu_odd_shift != 0.0 || op->mu + op->mu_even_shift != 0.0)
    add_tm_term_to_doublet_vectorized_layout_double(
        op->tm_term + start * tm_size, op->clover_doublet_vectorized + start * offset_doublet_v,
        n_per_core, nv);
#endif
#endif

#endif
}

void coarse_gamma5_double(vector_double eta, vector_double phi, int start, int end, Level *l) {

  int j, k = l->num_lattice_site_var / 2;
  vector_double eta_end;

  eta_end = eta + end;
  phi += start;
  eta += start;

  if (eta != phi) {
    while (eta < eta_end) {
      for (j = 0; j < k; j++) {
        *eta = -(*phi);
        eta++;
        phi++;
      }
      for (j = 0; j < k; j++) {
        *eta = *phi;
        eta++;
        phi++;
      }
    }
  } else {
    while (eta < eta_end) {
      for (j = 0; j < k; j++) {
        *eta = -(*eta);
        eta++;
      }
      eta += k;
    }
  }
}

void coarse_tau1_gamma5_double(vector_double eta, vector_double phi, int start, int end, Level *l) {

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2) {
    int j, k = l->num_lattice_site_var / 4;
    vector_double eta_end;

    eta_end = eta + end;
    phi += start;
    eta += start;

    ASSERT(eta != phi);
    while (eta < eta_end) {
      phi += k;
      for (j = 0; j < k; j++) {
        *eta = -(*phi);
        eta++;
        phi++;
      }
      phi -= 2 * k;
      for (j = 0; j < k; j++) {
        *eta = -(*phi);
        eta++;
        phi++;
      }
      phi += 2 * k;
      for (j = 0; j < k; j++) {
        *eta = *phi;
        eta++;
        phi++;
      }
      phi -= 2 * k;
      for (j = 0; j < k; j++) {
        *eta = *phi;
        eta++;
        phi++;
      }
      phi += k;
    }
  } else
#endif
  {
    warning0("coarse_tau1_gamma5_double called with global->n_flavours != 2\n");
    coarse_gamma5_double(eta, phi, start, end, l);
  }
}

void apply_coarse_operator_double(vector_double eta, vector_double phi, operator_struct<double> *op,
                                  Level *l) {

  //   PROF_double_START(_SC);
  int start = 0;
  int end = l->num_inner_lattice_sites;

#ifndef OPTIMIZED_COARSE_SELF_COUPLING_double
  coarse_self_couplings_double(eta, phi, op, start, end, l);
#else
  coarse_self_couplings_double_vectorized(eta, phi, op, start, end, l);
#endif
  //   PROF_double_STOP(_SC, 1);

  //   PROF_double_START(_NC);

#ifndef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_double
  coarse_hopping_term_double(eta, phi, op, _FULL_SYSTEM, l);
#else
  coarse_hopping_term_double_vectorized(eta, phi, op, _FULL_SYSTEM, l);
#endif
  //   PROF_double_STOP(_NC, 1);
}

void g5D_apply_coarse_operator_double(vector_double eta, vector_double phi,
                                      operator_struct<double> *op, Level *l) {

  apply_coarse_operator_double(eta, phi, op, l);
  coarse_gamma5_double(eta, eta, 0, l->inner_vector_size, l);
}

void apply_coarse_operator_dagger_double(vector_double eta, vector_double phi,
                                         operator_struct<double> *op, Level *l) {

  // NOTE: vector start and vector end COULD be wrong, not sure though
  coarse_gamma5_double(l->vbuf_double[3], phi, 0, l->vector_size, l);
  apply_coarse_operator_double(eta, l->vbuf_double[3], op, l);
  coarse_gamma5_double(eta, eta, 0, l->vector_size, l);
}
