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

void clover_double(vector_double eta, vector_double phi, operator_struct<double> *op, int start,
                   int end, Level *l) {

  int nv = l->num_lattice_site_var;
  vector_double lphi = phi + start, leta = eta + start;
  vector_double leta_end = eta + end;

#ifdef HAVE_TM
  config_double tm_term = op->tm_term + (start / nv) * 12;
#endif

  if (l->global->csw == 0.0) {

    config_double clover = op->clover + (start / nv) * 12;
#ifdef HAVE_TM1p1
    if (l->global->n_flavours == 2) {
#ifdef HAVE_TM
      if (l->global->mu + l->global->mu_odd_shift != 0.0 ||
          l->global->mu + l->global->mu_even_shift != 0.0)
        while (leta < leta_end) {
          FOR6(*leta = (*lphi) * ((*clover) + (*tm_term)); leta++; lphi++; clover++; tm_term++;);
          clover -= 6;
          tm_term -= 6;
          FOR6(*leta = (*lphi) * ((*clover) - (*tm_term)); leta++; lphi++; clover++; tm_term++;);
          FOR6(*leta = (*lphi) * ((*clover) + (*tm_term)); leta++; lphi++; clover++; tm_term++;);
          clover -= 6;
          tm_term -= 6;
          FOR6(*leta = (*lphi) * ((*clover) - (*tm_term)); leta++; lphi++; clover++; tm_term++;);
        }
      else
#endif
        while (leta < leta_end) {
          FOR6(*leta = (*lphi) * (*clover); leta++; lphi++; clover++;);
          clover -= 6;
          FOR12(*leta = (*lphi) * (*clover); leta++; lphi++; clover++;);
          clover -= 6;
          FOR6(*leta = (*lphi) * (*clover); leta++; lphi++; clover++;);
        }
    } else {
#endif
#ifdef HAVE_TM
      if (l->global->mu + l->global->mu_odd_shift != 0.0 ||
          l->global->mu + l->global->mu_even_shift != 0.0) {
        while (leta < leta_end)
          FOR12(*leta = (*lphi) * ((*clover) + (*tm_term)); leta++; lphi++; clover++; tm_term++;);
      } else
#endif
        while (leta < leta_end)
          FOR12(*leta = (*lphi) * (*clover); leta++; lphi++; clover++;);
#ifdef HAVE_TM1p1
    }
#endif

  } else {

#ifndef OPTIMIZED_SELF_COUPLING_double

    config_double clover = op->clover + (start / nv) * 42;
#ifdef HAVE_TM1p1
    if (l->global->n_flavours == 2) {
#ifdef HAVE_TM
      if (l->global->mu + l->global->mu_odd_shift != 0.0 ||
          l->global->mu + l->global->mu_even_shift != 0.0)
        while (leta < leta_end) {
          doublet_site_clover_double(leta, lphi, clover);
          clover += 42;
          FOR6(*leta += (*lphi) * (*tm_term); leta++; lphi++; tm_term++;);
          tm_term -= 6;
          FOR6(*leta -= (*lphi) * (*tm_term); leta++; lphi++; tm_term++;);
          FOR6(*leta += (*lphi) * (*tm_term); leta++; lphi++; tm_term++;);
          tm_term -= 6;
          FOR6(*leta -= (*lphi) * (*tm_term); leta++; lphi++; tm_term++;);
        }
      else
#endif
        while (leta < leta_end) {
          doublet_site_clover_double(leta, lphi, clover);
          leta += 24;
          lphi += 24;
          clover += 42;
        }
    } else {
#endif
#ifdef HAVE_TM
      if (l->global->mu + l->global->mu_odd_shift != 0.0 ||
          l->global->mu + l->global->mu_even_shift != 0.0)
        while (leta < leta_end) {
          site_clover_double(leta, lphi, clover);
          FOR12(*leta += (*lphi) * (*tm_term); leta++; lphi++; tm_term++;);
          clover += 42;
        }
      else
#endif
        while (leta < leta_end) {
          site_clover_double(leta, lphi, clover);
          leta += 12;
          lphi += 12;
          clover += 42;
        }
#ifdef HAVE_TM1p1
    }
#endif

#else

#ifdef HAVE_TM1p1
    double *clover =
        (l->global->n_flavours == 2) ? op->clover_doublet_vectorized : op->clover_vectorized;
#else
    double *clover = op->clover_vectorized;
#endif
    clover += start * 12;
    while (leta < leta_end) { // tm_term included in the clover vectorized
      sse_site_clover_double((double *)leta, (double *)lphi, clover);
      leta += nv;
      lphi += nv;
      clover += 12 * nv;
    }

#endif
  }

#ifdef HAVE_TM1p1
  config_double eps_term = op->epsbar_term + (start / nv) * 12;
  lphi = phi + start, leta = eta + start;
  if (l->global->n_flavours == 2 &&
      (l->global->epsbar != 0 || l->global->epsbar_ig5_odd_shift != 0 ||
       l->global->epsbar_ig5_odd_shift != 0))
    while (leta < leta_end) {
      lphi += 6;
      FOR6(*leta += (*lphi) * (*eps_term); leta++; lphi++; eps_term++;)
      lphi -= 12;
      eps_term -= 6;
      FOR6(*leta += (*lphi) * (*eps_term); leta++; lphi++; eps_term++;)
      lphi += 6;
    }
#endif

#ifdef PROFILING
  PROF_double_STOP(_SC, 1);
#endif
}

static void spin0and1_clover_double(vector_double eta, vector_double phi, config_double clover,
                                    Level *l) {

  vector_double eta_end = eta + l->inner_vector_size;
  if (l->global->csw == 0.0) {
    while (eta < eta_end) {
      FOR6(*eta = (*phi) * (*clover); eta++; phi++; clover++;)
      FOR6(*eta = 0; eta++;)
      phi += 6;
      clover += 6;
    }
  } else {
    while (eta < eta_end) {
      spin0and1_site_clover_double(eta, phi, clover);
      eta += 12;
      phi += 12;
      clover += 42;
    }
  }
}

static void spin2and3_clover_double(vector_double eta, vector_double phi, config_double clover,
                                    Level *l) {

  vector_double eta_end = eta + l->inner_vector_size;
  if (l->global->csw == 0.0) {
    while (eta < eta_end) {
      phi += 6;
      clover += 6;
      FOR6(*eta = 0; eta++;)
      FOR6(*eta = (*phi) * (*clover); eta++; phi++; clover++;)
    }
  } else {
    while (eta < eta_end) {
      spin2and3_site_clover_double(eta, phi, clover);
      eta += 12;
      phi += 12;
      clover += 42;
    }
  }
}

void block_d_plus_clover_double(vector_double eta, vector_double phi, int start, Schwarz *s,
                                Level *l) {

  int n = s->number_of_block_sites, *length = s->direction_length, **index = s->index_table,
      *neighbor = s->op.neighbor_table, nv = l->num_lattice_site_var;
  vector_double lphi = phi + start, leta = eta + start;

  // clover term
  clover_double(eta, phi, &(s->op), start, start + nv * n, l);

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double // block operator vectorized just in the double
                                          // environment
  double *Dplus = s->op.D_vectorized + (start / nv) * 96;
  double *Dminus = s->op.D_transformed_vectorized + (start / nv) * 96;
  for (int mu = 0; mu < 4; mu++) {
    block_oddeven_plus_coupling_double((double *)leta, Dplus, (double *)lphi, mu, 0, length[mu],
                                       index[mu], neighbor);
    block_oddeven_minus_coupling_double((double *)leta, Dminus, (double *)lphi, mu, 0, length[mu],
                                        index[mu], neighbor);
  }
#else
  int i, j, k, *ind;
  config_double D_pt;
  config_double D = s->op.D + (start / nv) * 36;
#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2) {
    complex_double buf1[50] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                   *buf2 = buf1 + 12, *buf3 = buf2 + 12, *buf4 = buf3 + 12;
    // inner block couplings
    ind = index[_T]; // T direction
    for (i = 0; i < length[_T]; i++) {
      k = ind[i];
      j = neighbor[4 * k + _T];
      D_pt = D + 36 * k + 9 * _T;
      dprn_T_double(buf1, lphi + 24 * k);    // (1+gamma_T) phi(x) + projection
      dprp_T_double(buf2, lphi + 24 * j);    // (1-gamma_T) phi(x+hat{T}) + projection
      mvmh_double(buf3, D_pt, buf1);         // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvmh_double(buf3 + 3, D_pt, buf1 + 3); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvmh_double(buf3 + 6, D_pt, buf1 + 6); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvmh_double(buf3 + 9, D_pt, buf1 + 9); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvm_double(buf4, D_pt, buf2);          // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      mvm_double(buf4 + 3, D_pt, buf2 + 3);  // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      mvm_double(buf4 + 6, D_pt, buf2 + 6);  // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      mvm_double(buf4 + 9, D_pt, buf2 + 9);  // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      dpbn_su3_T_double(
          buf3, leta + 24 * j); // eta(x+hat{T}) -= U_T(x)^dagger(x) (1+gamma_T) phi(x) + lift back
      dpbp_su3_T_double(buf4,
                        leta + 24 * k); // eta(x) -= U_T(x) (1-gamma_T) phi(x+hat{T}) + lift back
    }
    ind = index[_Z]; // Z direction
    for (i = 0; i < length[_Z]; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Z];
      D_pt = D + 36 * k + 9 * _Z;
      dprn_Z_double(buf1, lphi + 24 * k);    // (1+gamma_Z) phi(x) + projection
      dprp_Z_double(buf2, lphi + 24 * j);    // (1-gamma_Z) phi(x+hat{Z}) + projection
      mvmh_double(buf3, D_pt, buf1);         // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvmh_double(buf3 + 3, D_pt, buf1 + 3); // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvmh_double(buf3 + 6, D_pt, buf1 + 6); // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvmh_double(buf3 + 9, D_pt, buf1 + 9); // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvm_double(buf4, D_pt, buf2);          // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      mvm_double(buf4 + 3, D_pt, buf2 + 3);  // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      mvm_double(buf4 + 6, D_pt, buf2 + 6);  // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      mvm_double(buf4 + 9, D_pt, buf2 + 9);  // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      dpbn_su3_Z_double(
          buf3, leta + 24 * j); // eta(x+hat{Z}) -= U_Z(x)^dagger(x) (1+gamma_Z) phi(x) + lift back
      dpbp_su3_Z_double(buf4,
                        leta + 24 * k); // eta(x) -= U_Z(x) (1-gamma_Z) phi(x+hat{Z}) + lift back
    }
    ind = index[_Y]; // Y direction
    for (i = 0; i < length[_Y]; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Y];
      D_pt = D + 36 * k + 9 * _Y;
      dprn_Y_double(buf1, lphi + 24 * k);    // (1+gamma_Y) phi(x) + projection
      dprp_Y_double(buf2, lphi + 24 * j);    // (1-gamma_Y) phi(x+hat{Y}) + projection
      mvmh_double(buf3, D_pt, buf1);         // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvmh_double(buf3 + 3, D_pt, buf1 + 3); // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvmh_double(buf3 + 6, D_pt, buf1 + 6); // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvmh_double(buf3 + 9, D_pt, buf1 + 9); // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvm_double(buf4, D_pt, buf2);          // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      mvm_double(buf4 + 3, D_pt, buf2 + 3);  // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      mvm_double(buf4 + 6, D_pt, buf2 + 6);  // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      mvm_double(buf4 + 9, D_pt, buf2 + 9);  // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      dpbn_su3_Y_double(
          buf3, leta + 24 * j); // eta(x+hat{Y}) -= U_Y(x)^dagger(x) (1+gamma_Y) phi(x) + lift back
      dpbp_su3_Y_double(buf4,
                        leta + 24 * k); // eta(x) -= U_Y(x) (1-gamma_Y) phi(x+hat{Y}) + lift back
    }
    ind = index[_X]; // X direction
    for (i = 0; i < length[_X]; i++) {
      k = ind[i];
      j = neighbor[4 * k + _X];
      D_pt = D + 36 * k + 9 * _X;
      dprn_X_double(buf1, lphi + 24 * k);    // (1+gamma_X) phi(x) + projection
      dprp_X_double(buf2, lphi + 24 * j);    // (1-gamma_X) phi(x+hat{X}) + projection
      mvmh_double(buf3, D_pt, buf1);         // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvmh_double(buf3 + 3, D_pt, buf1 + 3); // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvmh_double(buf3 + 6, D_pt, buf1 + 6); // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvmh_double(buf3 + 9, D_pt, buf1 + 9); // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvm_double(buf4, D_pt, buf2);          // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      mvm_double(buf4 + 3, D_pt, buf2 + 3);  // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      mvm_double(buf4 + 6, D_pt, buf2 + 6);  // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      mvm_double(buf4 + 9, D_pt, buf2 + 9);  // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      dpbn_su3_X_double(
          buf3, leta + 24 * j); // eta(x+hat{X}) -= U_X(x)^dagger(x) (1+gamma_X) phi(x) + lift back
      dpbp_su3_X_double(buf4,
                        leta + 24 * k); // eta(x) -= U_X(x) (1-gamma_X) phi(x+hat{X}) + lift back
    }
  } else {
#endif
    complex_double buf1[25] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                   *buf2 = buf1 + 6, *buf3 = buf2 + 6, *buf4 = buf3 + 6;
    // inner block couplings
    ind = index[_T]; // T direction
    for (i = 0; i < length[_T]; i++) {
      k = ind[i];
      j = neighbor[4 * k + _T];
      D_pt = D + 36 * k + 9 * _T;
      prn_T_double(buf1, lphi + 12 * k);     // (1+gamma_T) phi(x) + projection
      prp_T_double(buf2, lphi + 12 * j);     // (1-gamma_T) phi(x+hat{T}) + projection
      mvmh_double(buf3, D_pt, buf1);         // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvmh_double(buf3 + 3, D_pt, buf1 + 3); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvm_double(buf4, D_pt, buf2);          // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      mvm_double(buf4 + 3, D_pt, buf2 + 3);  // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      pbn_su3_T_double(
          buf3, leta + 12 * j); // eta(x+hat{T}) -= U_T(x)^dagger(x) (1+gamma_T) phi(x) + lift back
      pbp_su3_T_double(buf4,
                       leta + 12 * k); // eta(x) -= U_T(x) (1-gamma_T) phi(x+hat{T}) + lift back
    }
    ind = index[_Z]; // Z direction
    for (i = 0; i < length[_Z]; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Z];
      D_pt = D + 36 * k + 9 * _Z;
      prn_Z_double(buf1, lphi + 12 * k);     // (1+gamma_Z) phi(x) + projection
      prp_Z_double(buf2, lphi + 12 * j);     // (1-gamma_Z) phi(x+hat{Z}) + projection
      mvmh_double(buf3, D_pt, buf1);         // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvmh_double(buf3 + 3, D_pt, buf1 + 3); // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvm_double(buf4, D_pt, buf2);          // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      mvm_double(buf4 + 3, D_pt, buf2 + 3);  // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      pbn_su3_Z_double(
          buf3, leta + 12 * j); // eta(x+hat{Z}) -= U_Z(x)^dagger(x) (1+gamma_Z) phi(x) + lift back
      pbp_su3_Z_double(buf4,
                       leta + 12 * k); // eta(x) -= U_Z(x) (1-gamma_Z) phi(x+hat{Z}) + lift back
    }
    ind = index[_Y]; // Y direction
    for (i = 0; i < length[_Y]; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Y];
      D_pt = D + 36 * k + 9 * _Y;
      prn_Y_double(buf1, lphi + 12 * k);     // (1+gamma_Y) phi(x) + projection
      prp_Y_double(buf2, lphi + 12 * j);     // (1-gamma_Y) phi(x+hat{Y}) + projection
      mvmh_double(buf3, D_pt, buf1);         // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvmh_double(buf3 + 3, D_pt, buf1 + 3); // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvm_double(buf4, D_pt, buf2);          // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      mvm_double(buf4 + 3, D_pt, buf2 + 3);  // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      pbn_su3_Y_double(
          buf3, leta + 12 * j); // eta(x+hat{Y}) -= U_Y(x)^dagger(x) (1+gamma_Y) phi(x) + lift back
      pbp_su3_Y_double(buf4,
                       leta + 12 * k); // eta(x) -= U_Y(x) (1-gamma_Y) phi(x+hat{Y}) + lift back
    }
    ind = index[_X]; // X direction
    for (i = 0; i < length[_X]; i++) {
      k = ind[i];
      j = neighbor[4 * k + _X];
      D_pt = D + 36 * k + 9 * _X;
      prn_X_double(buf1, lphi + 12 * k);     // (1+gamma_X) phi(x) + projection
      prp_X_double(buf2, lphi + 12 * j);     // (1-gamma_X) phi(x+hat{X}) + projection
      mvmh_double(buf3, D_pt, buf1);         // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvmh_double(buf3 + 3, D_pt, buf1 + 3); // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvm_double(buf4, D_pt, buf2);          // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      mvm_double(buf4 + 3, D_pt, buf2 + 3);  // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      pbn_su3_X_double(
          buf3, leta + 12 * j); // eta(x+hat{X}) -= U_X(x)^dagger(x) (1+gamma_X) phi(x) + lift back
      pbp_su3_X_double(buf4,
                       leta + 12 * k); // eta(x) -= U_X(x) (1-gamma_X) phi(x+hat{X}) + lift back
    }
#ifdef HAVE_TM1p1
  }
#endif
#endif
}

void d_plus_clover_double(vector_double eta, vector_double phi, operator_struct<double> *op,
                          Level *l) {

  int n = l->num_inner_lattice_sites, *neighbor = op->neighbor_table, nv = l->num_lattice_site_var;
  int start = 0, end =  nv * n;

  auto *prnT = new complex_double[nv/2 * l->num_lattice_sites];
  auto *prnZ = new complex_double[nv/2 * l->num_lattice_sites];
  auto *prnY = new complex_double[nv/2 * l->num_lattice_sites];
  auto *prnX = new complex_double[nv/2 * l->num_lattice_sites];

  auto *prpT = new complex_double[nv/2 * l->num_lattice_sites];
  auto *prpZ = new complex_double[nv/2 * l->num_lattice_sites];
  auto *prpY = new complex_double[nv/2 * l->num_lattice_sites];
  auto *prpX = new complex_double[nv/2 * l->num_lattice_sites];

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
  complex_double *prn[4] = {op->prnT, op->prnZ, op->prnY, op->prnX};
  complex_double *prp[4] = {op->prpT, op->prpZ, op->prpY, op->prpX};
#else
  int i, j, *nb_pt;
  vector_double phi_pt, eta_pt, end_pt;
  config_double D_pt;
#endif

  clover_double(eta, phi, op, start, end, l);

  PROF_double_START(_NC);

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2) {
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    dprp_double(prn, phi, start, end);
#else
    complex_double pbuf[12];
    for (i = start / 2, phi_pt = phi + start; i < end / 2; i += nv / 2, phi_pt += nv) {
      dprp_T_double(op->prnT + i, phi_pt);
      dprp_Z_double(op->prnZ + i, phi_pt);
      dprp_Y_double(op->prnY + i, phi_pt);
      dprp_X_double(op->prnX + i, phi_pt);
    }
#endif
    // start communication in negative direction
    ghost_sendrecv_double(op->prnT, _T, -1, &(op->c), _FULL_SYSTEM, l);
    ghost_sendrecv_double(op->prnZ, _Z, -1, &(op->c), _FULL_SYSTEM, l);
    ghost_sendrecv_double(op->prnY, _Y, -1, &(op->c), _FULL_SYSTEM, l);
    ghost_sendrecv_double(op->prnX, _X, -1, &(op->c), _FULL_SYSTEM, l);

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    dprn_su3_double(prp, phi, op, neighbor, start, end);
#else
    // project plus dir and multiply with U dagger
    for (phi_pt = phi + start, end_pt = phi + end, D_pt = op->D + ((start / nv) * 36),
        nb_pt = neighbor + ((start / nv) * 4);
         phi_pt < end_pt; phi_pt += nv) {
      // T dir
      j = nv / 2 * (*nb_pt);
      nb_pt++;
      dprn_T_double(pbuf, phi_pt);
      mvmh_double(op->prpT + j, D_pt, pbuf);
      mvmh_double(op->prpT + j + 3, D_pt, pbuf + 3);
      mvmh_double(op->prpT + j + 6, D_pt, pbuf + 6);
      mvmh_double(op->prpT + j + 9, D_pt, pbuf + 9);
      D_pt += 9;
      // Z dir
      j = nv / 2 * (*nb_pt);
      nb_pt++;
      dprn_Z_double(pbuf, phi_pt);
      mvmh_double(op->prpZ + j, D_pt, pbuf);
      mvmh_double(op->prpZ + j + 3, D_pt, pbuf + 3);
      mvmh_double(op->prpZ + j + 6, D_pt, pbuf + 6);
      mvmh_double(op->prpZ + j + 9, D_pt, pbuf + 9);
      D_pt += 9;
      // Y dir
      j = nv / 2 * (*nb_pt);
      nb_pt++;
      dprn_Y_double(pbuf, phi_pt);
      mvmh_double(op->prpY + j, D_pt, pbuf);
      mvmh_double(op->prpY + j + 3, D_pt, pbuf + 3);
      mvmh_double(op->prpY + j + 6, D_pt, pbuf + 6);
      mvmh_double(op->prpY + j + 9, D_pt, pbuf + 9);
      D_pt += 9;
      // X dir
      j = nv / 2 * (*nb_pt);
      nb_pt++;
      dprn_X_double(pbuf, phi_pt);
      mvmh_double(op->prpX + j, D_pt, pbuf);
      mvmh_double(op->prpX + j + 3, D_pt, pbuf + 3);
      mvmh_double(op->prpX + j + 6, D_pt, pbuf + 6);
      mvmh_double(op->prpX + j + 9, D_pt, pbuf + 9);
      D_pt += 9;
    }
#endif

    // start communication in positive direction

    ghost_sendrecv_double(op->prpT, _T, +1, &(op->c), _FULL_SYSTEM, l);
    ghost_sendrecv_double(op->prpZ, _Z, +1, &(op->c), _FULL_SYSTEM, l);
    ghost_sendrecv_double(op->prpY, _Y, +1, &(op->c), _FULL_SYSTEM, l);
    ghost_sendrecv_double(op->prpX, _X, +1, &(op->c), _FULL_SYSTEM, l);
    // wait for communication in negative direction
    ghost_wait_double(op->prnT, _T, -1, &(op->c), _FULL_SYSTEM, l);
    ghost_wait_double(op->prnZ, _Z, -1, &(op->c), _FULL_SYSTEM, l);
    ghost_wait_double(op->prnY, _Y, -1, &(op->c), _FULL_SYSTEM, l);
    ghost_wait_double(op->prnX, _X, -1, &(op->c), _FULL_SYSTEM, l);

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    su3_dpbp_double(eta, prn, op, neighbor, start, end);
#else
    // multiply with U and lift up minus dir
    for (eta_pt = eta + start, end_pt = eta + end, D_pt = op->D + (start / nv) * 36,
        nb_pt = neighbor + (start / nv) * 4;
         eta_pt < end_pt; eta_pt += nv) {
      // T dir
      j = nv / 2 * (*nb_pt);
      nb_pt++;
      mvm_double(pbuf, D_pt, op->prnT + j);
      mvm_double(pbuf + 3, D_pt, op->prnT + j + 3);
      mvm_double(pbuf + 6, D_pt, op->prnT + j + 6);
      mvm_double(pbuf + 9, D_pt, op->prnT + j + 9);
      dpbp_su3_T_double(pbuf, eta_pt);
      D_pt += 9;
      // Z dir
      j = nv / 2 * (*nb_pt);
      nb_pt++;
      mvm_double(pbuf, D_pt, op->prnZ + j);
      mvm_double(pbuf + 3, D_pt, op->prnZ + j + 3);
      mvm_double(pbuf + 6, D_pt, op->prnZ + j + 6);
      mvm_double(pbuf + 9, D_pt, op->prnZ + j + 9);
      dpbp_su3_Z_double(pbuf, eta_pt);
      D_pt += 9;
      // Y dir
      j = nv / 2 * (*nb_pt);
      nb_pt++;
      mvm_double(pbuf, D_pt, op->prnY + j);
      mvm_double(pbuf + 3, D_pt, op->prnY + j + 3);
      mvm_double(pbuf + 6, D_pt, op->prnY + j + 6);
      mvm_double(pbuf + 9, D_pt, op->prnY + j + 9);
      dpbp_su3_Y_double(pbuf, eta_pt);
      D_pt += 9;
      // X dir
      j = nv / 2 * (*nb_pt);
      nb_pt++;
      mvm_double(pbuf, D_pt, op->prnX + j);
      mvm_double(pbuf + 3, D_pt, op->prnX + j + 3);
      mvm_double(pbuf + 6, D_pt, op->prnX + j + 6);
      mvm_double(pbuf + 9, D_pt, op->prnX + j + 9);
      dpbp_su3_X_double(pbuf, eta_pt);
      D_pt += 9;
    }
#endif

    // wait for communication in positive direction

    ghost_wait_double(op->prpT, _T, +1, &(op->c), _FULL_SYSTEM, l);
    ghost_wait_double(op->prpZ, _Z, +1, &(op->c), _FULL_SYSTEM, l);
    ghost_wait_double(op->prpY, _Y, +1, &(op->c), _FULL_SYSTEM, l);
    ghost_wait_double(op->prpX, _X, +1, &(op->c), _FULL_SYSTEM, l);

    // lift up plus dir
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    dpbn_double(eta, prp, start, end);
#else
    for (i = start / 2, eta_pt = eta + start; i < end / 2; i += 12, eta_pt += 24) {
      dpbn_su3_T_double(op->prpT + i, eta_pt);
      dpbn_su3_Z_double(op->prpZ + i, eta_pt);
      dpbn_su3_Y_double(op->prpY + i, eta_pt);
      dpbn_su3_X_double(op->prpX + i, eta_pt);
    }
#endif
  } else {
#endif

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    prp_double(prn, phi, start, end);
#else
  complex_double pbuf[6];
  /// compute (I4 - Gamma_mu) kronecker I_3 * psi(x)
  for (i = start / 2, phi_pt = phi + start; i < end / 2; i += 6, phi_pt += 12) {
    prp_T_double(prnT + i, phi_pt);
    prp_Z_double(prnZ + i, phi_pt);
    prp_Y_double(prnY + i, phi_pt);
    prp_X_double(prnX + i, phi_pt);
  }
#endif
    // start communication in negative direction
    /// send prn (pi_mu_minus kronecker psi(x+mu_hat)) to mu-1 neighbours,
    /// receive corresponding prn from mu+1 neighbours
    ghost_sendrecv_double(prnT, _T, -1, &(op->c), _FULL_SYSTEM, l);
    ghost_sendrecv_double(prnZ, _Z, -1, &(op->c), _FULL_SYSTEM, l);
    ghost_sendrecv_double(prnY, _Y, -1, &(op->c), _FULL_SYSTEM, l);
    ghost_sendrecv_double(prnX, _X, -1, &(op->c), _FULL_SYSTEM, l);

    // project plus dir and multiply with U dagger
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    prn_su3_double(prp, phi, op, neighbor, start, end);
#else
  for (phi_pt = phi + start, end_pt = phi + end, D_pt = op->D + (start * 3),
      nb_pt = neighbor + ((start / 12) * 4);
       phi_pt < end_pt; phi_pt += 12) {
    // T dir
    j = 6 * (*nb_pt);
    nb_pt++;
    prn_T_double(pbuf, phi_pt); /// (I4 + Gamma_mu) kronecker I3 * psi(x)
    mvmh_double(prpT + j, D_pt, pbuf); /// apply U_mu^H(x), first half
    mvmh_double(prpT + j + 3, D_pt, pbuf + 3); /// second half
    D_pt += 9;
    // Z dir
    j = 6 * (*nb_pt);
    nb_pt++;
    prn_Z_double(pbuf, phi_pt);
    mvmh_double(prpZ + j, D_pt, pbuf);
    mvmh_double(prpZ + j + 3, D_pt, pbuf + 3);
    D_pt += 9;
    // Y dir
    j = 6 * (*nb_pt);
    nb_pt++;
    prn_Y_double(pbuf, phi_pt);
    mvmh_double(prpY + j, D_pt, pbuf);
    mvmh_double(prpY + j + 3, D_pt, pbuf + 3);
    D_pt += 9;
    // X dir
    j = 6 * (*nb_pt);
    nb_pt++;
    prn_X_double(pbuf, phi_pt);
    mvmh_double(prpX + j, D_pt, pbuf);
    mvmh_double(prpX + j + 3, D_pt, pbuf + 3);
    D_pt += 9;
  }
#endif

    // start communication in positive direction
    /// send prp to mu+1 neighbour, receive corresponding prp from mu-1 neighbour
    ghost_sendrecv_double(prpT, _T, +1, &(op->c), _FULL_SYSTEM, l);
    ghost_sendrecv_double(prpZ, _Z, +1, &(op->c), _FULL_SYSTEM, l);
    ghost_sendrecv_double(prpY, _Y, +1, &(op->c), _FULL_SYSTEM, l);
    ghost_sendrecv_double(prpX, _X, +1, &(op->c), _FULL_SYSTEM, l);

    // wait for communication in negative direction
    ghost_wait_double(prnT, _T, -1, &(op->c), _FULL_SYSTEM, l);
    ghost_wait_double(prnZ, _Z, -1, &(op->c), _FULL_SYSTEM, l);
    ghost_wait_double(prnY, _Y, -1, &(op->c), _FULL_SYSTEM, l);
    ghost_wait_double(prnX, _X, -1, &(op->c), _FULL_SYSTEM, l);

    // multiply with U and lift up minus dir
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    su3_pbp_double(eta, prn, op, neighbor, start, end);
#else
  for (eta_pt = eta + start, end_pt = eta + end, D_pt = op->D + start * 3,
      nb_pt = neighbor + (start / 12) * 4;
       eta_pt < end_pt; eta_pt += 12) {
    // T dir
    j = 6 * (*nb_pt);
    nb_pt++;
    mvm_double(pbuf, D_pt, prnT + j); /// apply U_mu(x) to prn, first half of x
    mvm_double(pbuf + 3, D_pt, prnT + j + 3); /// second half of x
    /// each line of second half of eta can be represented by lines of first half of eta,
    /// here we apply the effect of pi_mu_minus kronecker U_mu(x) * psi(x+mu_hat)
    pbp_su3_T_double(pbuf, eta_pt);
    D_pt += 9;
    // Z dir
    j = 6 * (*nb_pt);
    nb_pt++;
    mvm_double(pbuf, D_pt, prnZ + j);
    mvm_double(pbuf + 3, D_pt, prnZ + j + 3);
    pbp_su3_Z_double(pbuf, eta_pt);
    D_pt += 9;
    // Y dir
    j = 6 * (*nb_pt);
    nb_pt++;
    mvm_double(pbuf, D_pt, prnY + j);
    mvm_double(pbuf + 3, D_pt, prnY + j + 3);
    pbp_su3_Y_double(pbuf, eta_pt);
    D_pt += 9;
    // X dir
    j = 6 * (*nb_pt);
    nb_pt++;
    mvm_double(pbuf, D_pt, prnX + j);
    mvm_double(pbuf + 3, D_pt, prnX + j + 3);
    pbp_su3_X_double(pbuf, eta_pt);
    D_pt += 9;
  }
#endif

    // wait for communication in positive direction
    ghost_wait_double(prpT, _T, +1, &(op->c), _FULL_SYSTEM, l);
    ghost_wait_double(prpZ, _Z, +1, &(op->c), _FULL_SYSTEM, l);
    ghost_wait_double(prpY, _Y, +1, &(op->c), _FULL_SYSTEM, l);
    ghost_wait_double(prpX, _X, +1, &(op->c), _FULL_SYSTEM, l);

    // lift up plus dir
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    pbn_double(eta, prp, start, end);
#else
  for (i = start / 2, eta_pt = eta + start; i < end / 2; i += 6, eta_pt += 12) {
    /// each line of second half of eta can be represented by lines of first half of eta,
    /// here we apply the effect of pi_mu_plus kronecker U_mu(x-mu_hat) * psi(x-mu_hat)
    pbn_su3_T_double(prpT + i, eta_pt);
    pbn_su3_Z_double(prpZ + i, eta_pt);
    pbn_su3_Y_double(prpY + i, eta_pt);
    pbn_su3_X_double(prpX + i, eta_pt);
  }
#endif
#ifdef HAVE_TM1p1
  }
#endif

  delete[] prnT;
  delete[] prnZ;
  delete[] prnY;
  delete[] prnX;

  delete[] prpT;
  delete[] prpZ;
  delete[] prpY;
  delete[] prpX;
  PROF_double_STOP(_NC, 1);
}

void gamma5_double(vector_double eta, vector_double phi, Level *l) {

  ASSERT(l->type == _FINE);

  vector_double eta_end = eta + l->inner_vector_size;
#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2) {
    while (eta < eta_end) {
      FOR12(*eta = -(*phi); phi++; eta++;)
      FOR12(*eta = (*phi); phi++; eta++;)
    }
  } else
#endif
    while (eta < eta_end) {
      FOR6(*eta = -(*phi); phi++; eta++;)
      FOR6(*eta = (*phi); phi++; eta++;)
    }
}

void tau1_gamma5_double(vector_double eta, vector_double phi, Level *l) {

  ASSERT(l->type == _FINE);

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2) {
    vector_double eta_end = eta + l->inner_vector_size;
    complex_double b[6];
    while (eta < eta_end) {
      int i = 0;
      FOR6(b[i] = (*phi); phi++; i++;);
      FOR6(*eta = -(*phi); phi++; eta++;);
      i = 0;
      FOR6(*eta = -b[i]; eta++; i++;);
      i = 0;
      FOR6(b[i] = (*phi); phi++; i++;);
      FOR6(*eta = (*phi); phi++; eta++;);
      i = 0;
      FOR6(*eta = b[i]; eta++; i++;);
    }
  } else
#endif
  {

    warning0("tau1_gamma5_double called with l->global->n_flavours != 2\n");

    gamma5_double(eta, phi, l);
  }
}

void set_even_to_zero_double(vector_double eta, vector_double phi, Level *l) {

  ASSERT(l->type == _FINE);

  int i = 0;
  vector_double eta_end = eta + l->inner_vector_size;

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2)
    while (eta < eta_end) {
      if (l->global->odd_even_table[i] == _ODD) {
        FOR24(*eta = (*phi); phi++; eta++;);
      } else if (l->global->odd_even_table[i] == _EVEN) {
        FOR24(*eta = 0; phi++; eta++;);
      }
      i++;
    }
  else
#endif
    while (eta < eta_end) {
      if (l->global->odd_even_table[i] == _ODD) {
        FOR12(*eta = (*phi); phi++; eta++;);
      } else if (l->global->odd_even_table[i] == _EVEN) {
        FOR12(*eta = 0; phi++; eta++;);
      }
      i++;
    }
}

void gamma5_set_even_to_zero_double(vector_double eta, vector_double phi, Level *l) {

  ASSERT(l->type == _FINE);

  int i = 0;
  vector_double eta_end = eta + l->inner_vector_size;

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2)
    while (eta < eta_end) {
      if (l->global->odd_even_table[i] == _ODD) {
        FOR12(*eta = -(*phi); phi++; eta++;);
        FOR12(*eta = (*phi); phi++; eta++;);
      } else if (l->global->odd_even_table[i] == _EVEN) {
        FOR24(*eta = 0; phi++; eta++;);
      }
      i++;
    }
  else
#endif
    while (eta < eta_end) {
      if (l->global->odd_even_table[i] == _ODD) {
        FOR6(*eta = -(*phi); phi++; eta++;);
        FOR6(*eta = (*phi); phi++; eta++;);
      } else if (l->global->odd_even_table[i] == _EVEN) {
        FOR12(*eta = 0; phi++; eta++;);
      }
      i++;
    }
}

void tau1_gamma5_set_even_to_zero_double(vector_double eta, vector_double phi, Level *l) {

  ASSERT(l->type == _FINE);

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2) {
    int i = 0;
    vector_double eta_end = eta + l->inner_vector_size;

    complex_double b[6];
    while (eta < eta_end) {
      if (l->global->odd_even_table[i] == _ODD) {
        int i = 0;
        FOR6(b[i] = (*phi); phi++; i++;);
        FOR6(*eta = -(*phi); phi++; eta++;);
        i = 0;
        FOR6(*eta = -b[i]; eta++; i++;);
        i = 0;
        FOR6(b[i] = (*phi); phi++; i++;);
        FOR6(*eta = (*phi); phi++; eta++;);
        i = 0;
        FOR6(*eta = b[i]; eta++; i++;);
      } else if (l->global->odd_even_table[i] == _EVEN) {
        FOR24(*eta = 0; phi++; eta++;);
      }
      i++;
    }
  } else
#endif
  {

    warning0("tau1_gamma5_set_even_to_zero_double called with l->global->n_flavours != 2\n");

    gamma5_set_even_to_zero_double(eta, phi, l);
  }
}

void set_odd_to_zero_double(vector_double eta, vector_double phi, Level *l) {

  int i = 0;
  vector_double eta_end = eta + l->inner_vector_size;

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2)
    while (eta < eta_end) {
      if (l->global->odd_even_table[i] == _EVEN) {
        FOR24(*eta = (*phi); phi++; eta++;);
      } else if (l->global->odd_even_table[i] == _ODD) {
        FOR24(*eta = 0; phi++; eta++;);
      }
      i++;
    }
  else
#endif
    while (eta < eta_end) {
      if (l->global->odd_even_table[i] == _EVEN) {
        FOR12(*eta = (*phi); phi++; eta++;);
      } else if (l->global->odd_even_table[i] == _ODD) {
        FOR12(*eta = 0; phi++; eta++;);
      }
      i++;
    }
}

void gamma5_set_odd_to_zero_double(vector_double eta, vector_double phi, Level *l) {

  int i = 0;
  vector_double eta_end = eta + l->inner_vector_size;

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2)
    while (eta < eta_end) {
      if (l->global->odd_even_table[i] == _EVEN) {
        FOR12(*eta = -(*phi); phi++; eta++;);
        FOR12(*eta = (*phi); phi++; eta++;);
      } else if (l->global->odd_even_table[i] == _ODD) {
        FOR24(*eta = 0; phi++; eta++;);
      }
      i++;
    }
  else
#endif
    while (eta < eta_end) {
      if (l->global->odd_even_table[i] == _EVEN) {
        FOR6(*eta = -(*phi); phi++; eta++;);
        FOR6(*eta = (*phi); phi++; eta++;);
      } else if (l->global->odd_even_table[i] == _ODD) {
        FOR12(*eta = 0; phi++; eta++;);
      }
      i++;
    }
}

void tau1_gamma5_set_odd_to_zero_double(vector_double eta, vector_double phi, Level *l) {

  ASSERT(l->type == _FINE);

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2) {
    int i = 0;
    vector_double eta_end = eta + l->inner_vector_size;

    complex_double b[6];
    while (eta < eta_end) {
      if (l->global->odd_even_table[i] == _EVEN) {
        int i = 0;
        FOR6(b[i] = (*phi); phi++; i++;);
        FOR6(*eta = -(*phi); phi++; eta++;);
        i = 0;
        FOR6(*eta = -b[i]; eta++; i++;);
        i = 0;
        FOR6(b[i] = (*phi); phi++; i++;);
        FOR6(*eta = (*phi); phi++; eta++;);
        i = 0;
        FOR6(*eta = b[i]; eta++; i++;);
      } else if (l->global->odd_even_table[i] == _ODD) {
        FOR24(*eta = 0; phi++; eta++;);
      }
      i++;
    }
  } else
#endif
  {

    warning0("tau1_gamma5_set_odd_to_zero_double called with l->global->n_flavours != 2\n");

    gamma5_set_odd_to_zero_double(eta, phi, l);
  }
}

void scale_even_odd_double(vector_double eta, vector_double phi, complex_double even,
                           complex_double odd, Level *l) {

  int i = 0;
  vector_double eta_end = eta + l->inner_vector_size;

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2)
    while (eta < eta_end) {
      if (l->global->odd_even_table[i] == _EVEN) {
        FOR24(*eta = even * (*phi); phi++; eta++;);
      } else if (l->global->odd_even_table[i] == _ODD) {
        FOR24(*eta = odd * (*phi); phi++; eta++;);
      }
      i++;
    }
  else
#endif
    while (eta < eta_end) {
      if (l->global->odd_even_table[i] == _EVEN) {
        FOR12(*eta = even * (*phi); phi++; eta++;);
      } else if (l->global->odd_even_table[i] == _ODD) {
        FOR12(*eta = odd * (*phi); phi++; eta++;);
      }
      i++;
    }
}

void two_flavours_to_serial_double(vector_double flav1, vector_double flav2, vector_double serial,
                                   Level *l) {

#ifdef HAVE_TM1p1

  /*
   * Order: spin0and1 of flav1
   *        spin0and1 of flav2
   *        spin2and3 of flav1
   *        spin2and3 of flav2
   */
  vector_double serial_end;

  if (l->global->n_flavours == 2) {
    serial_end = serial + l->inner_vector_size;
  } else {
    serial_end = serial + l->inner_vector_size * 2;
  }

  while (serial < serial_end) {
    FOR6(*serial = (*flav1); serial++; flav1++;)
    FOR6(*serial = (*flav2); serial++; flav2++;)
    FOR6(*serial = (*flav1); serial++; flav1++;)
    FOR6(*serial = (*flav2); serial++; flav2++;)
  }
#else

  warning0("two_flavours_to_serial_double called without HAVE_TM1p1 defined\n");

#endif
}

void serial_to_two_flavours_double(vector_double flav1, vector_double flav2, vector_double serial,
                                   Level *l) {

#ifdef HAVE_TM1p1
  vector_double serial_end;

  if (l->global->n_flavours == 2) {
    serial_end = serial + l->inner_vector_size;

  } else {
    serial_end = serial + l->inner_vector_size * 2;
  }

  while (serial < serial_end) {
    FOR6(*flav1 = (*serial); serial++; flav1++;)
    FOR6(*flav2 = (*serial); serial++; flav2++;)
    FOR6(*flav1 = (*serial); serial++; flav1++;)
    FOR6(*flav2 = (*serial); serial++; flav2++;)
  }
#else

  warning0("two_flavours_to_serial_double called without HAVE_TM1p1 defined\n");

#endif
}

void g5D_plus_clover_double(vector_double eta, vector_double phi, operator_struct<double> *op,
                            Level *l) {
  d_plus_clover_double(eta, phi, op, l);
  gamma5_double(eta, eta, l);
}

void diagonal_aggregate_double(vector_double eta1, vector_double eta2, vector_double phi,
                               config_double diag, Level *l) {

  vector_double eta_end = eta1 + l->inner_vector_size;

  while (eta1 < eta_end) {
    FOR6(*eta1 = (*phi) * (*diag); *eta2 = 0; eta1++; eta2++; phi++; diag++;);
    FOR6(*eta2 = (*phi) * (*diag); *eta1 = 0; eta1++; eta2++; phi++; diag++;);
  }
}

void d_plus_clover_aggregate_double(vector_double eta1, vector_double eta2, vector_double phi,
                                    Schwarz *s, Level *l) {

  int i, length, index1, index2, *index_dir, *neighbor = s->op.neighbor_table;
  vector_double eta1_pt, eta2_pt, phi_pt;
  complex_double buffer1[12], buffer2[12];
  config_double D_pt, D = s->op.D;

  // add clover term/shift
  spin0and1_clover_double(eta1, phi, s->op.clover, l);
  spin2and3_clover_double(eta2, phi, s->op.clover, l);

  // T dir
  length = l->is_double.agg_length[_T];
  index_dir = l->is_double.agg_index[_T];
  for (i = 0; i < length; i++) {
    index1 = index_dir[i];
    index2 = neighbor[4 * index1 + _T];
    phi_pt = phi + 12 * index2;
    D_pt = D + 36 * index1 + 9 * _T;
    mvm_double(buffer1, D_pt, phi_pt);
    mvm_double(buffer1 + 3, D_pt, phi_pt + 3);
    mvm_double(buffer1 + 6, D_pt, phi_pt + 6);
    mvm_double(buffer1 + 9, D_pt, phi_pt + 9);
    phi_pt = phi + 12 * index1;
    mvmh_double(buffer2, D_pt, phi_pt);
    mvmh_double(buffer2 + 3, D_pt, phi_pt + 3);
    mvmh_double(buffer2 + 6, D_pt, phi_pt + 6);
    mvmh_double(buffer2 + 9, D_pt, phi_pt + 9);
    eta1_pt = eta1 + 12 * index1;
    eta2_pt = eta2 + 12 * index1;
    twospin_p_T_double(eta1_pt, eta2_pt, buffer1);
    eta1_pt = eta1 + 12 * index2;
    eta2_pt = eta2 + 12 * index2;
    twospin_n_T_double(eta1_pt, eta2_pt, buffer2);
  }
  // Z dir
  length = l->is_double.agg_length[_Z];
  index_dir = l->is_double.agg_index[_Z];
  for (i = 0; i < length; i++) {
    index1 = index_dir[i];
    index2 = neighbor[4 * index1 + _Z];
    phi_pt = phi + 12 * index2;
    D_pt = D + 36 * index1 + 9 * _Z;
    mvm_double(buffer1, D_pt, phi_pt);
    mvm_double(buffer1 + 3, D_pt, phi_pt + 3);
    mvm_double(buffer1 + 6, D_pt, phi_pt + 6);
    mvm_double(buffer1 + 9, D_pt, phi_pt + 9);
    phi_pt = phi + 12 * index1;
    mvmh_double(buffer2, D_pt, phi_pt);
    mvmh_double(buffer2 + 3, D_pt, phi_pt + 3);
    mvmh_double(buffer2 + 6, D_pt, phi_pt + 6);
    mvmh_double(buffer2 + 9, D_pt, phi_pt + 9);
    eta1_pt = eta1 + 12 * index1;
    eta2_pt = eta2 + 12 * index1;
    twospin_p_Z_double(eta1_pt, eta2_pt, buffer1);
    eta1_pt = eta1 + 12 * index2;
    eta2_pt = eta2 + 12 * index2;
    twospin_n_Z_double(eta1_pt, eta2_pt, buffer2);
  }
  // Y dir
  length = l->is_double.agg_length[_Y];
  index_dir = l->is_double.agg_index[_Y];
  for (i = 0; i < length; i++) {
    index1 = index_dir[i];
    index2 = neighbor[4 * index1 + _Y];
    phi_pt = phi + 12 * index2;
    D_pt = D + 36 * index1 + 9 * _Y;
    mvm_double(buffer1, D_pt, phi_pt);
    mvm_double(buffer1 + 3, D_pt, phi_pt + 3);
    mvm_double(buffer1 + 6, D_pt, phi_pt + 6);
    mvm_double(buffer1 + 9, D_pt, phi_pt + 9);
    phi_pt = phi + 12 * index1;
    mvmh_double(buffer2, D_pt, phi_pt);
    mvmh_double(buffer2 + 3, D_pt, phi_pt + 3);
    mvmh_double(buffer2 + 6, D_pt, phi_pt + 6);
    mvmh_double(buffer2 + 9, D_pt, phi_pt + 9);
    eta1_pt = eta1 + 12 * index1;
    eta2_pt = eta2 + 12 * index1;
    twospin_p_Y_double(eta1_pt, eta2_pt, buffer1);
    eta1_pt = eta1 + 12 * index2;
    eta2_pt = eta2 + 12 * index2;
    twospin_n_Y_double(eta1_pt, eta2_pt, buffer2);
  }
  // X dir
  length = l->is_double.agg_length[_X];
  index_dir = l->is_double.agg_index[_X];
  for (i = 0; i < length; i++) {
    index1 = index_dir[i];
    index2 = neighbor[4 * index1 + _X];
    phi_pt = phi + 12 * index2;
    D_pt = D + 36 * index1 + 9 * _X;
    mvm_double(buffer1, D_pt, phi_pt);
    mvm_double(buffer1 + 3, D_pt, phi_pt + 3);
    mvm_double(buffer1 + 6, D_pt, phi_pt + 6);
    mvm_double(buffer1 + 9, D_pt, phi_pt + 9);
    phi_pt = phi + 12 * index1;
    mvmh_double(buffer2, D_pt, phi_pt);
    mvmh_double(buffer2 + 3, D_pt, phi_pt + 3);
    mvmh_double(buffer2 + 6, D_pt, phi_pt + 6);
    mvmh_double(buffer2 + 9, D_pt, phi_pt + 9);
    eta1_pt = eta1 + 12 * index1;
    eta2_pt = eta2 + 12 * index1;
    twospin_p_X_double(eta1_pt, eta2_pt, buffer1);
    eta1_pt = eta1 + 12 * index2;
    eta2_pt = eta2 + 12 * index2;
    twospin_n_X_double(eta1_pt, eta2_pt, buffer2);
  }
}

void d_neighbor_aggregate_double(vector_double eta1, vector_double eta2, vector_double phi,
                                 const int mu, Schwarz *s, Level *l) {

  int i, length, index1, index2, *index_dir, *neighbor;
  vector_double eta1_pt, eta2_pt, phi_pt;
  complex_double buffer1[12];
  config_double D_pt, D = s->op.D;

  length = l->is_double.agg_boundary_length[mu];
  index_dir = l->is_double.agg_boundary_index[mu];
  neighbor = l->is_double.agg_boundary_neighbor[mu];

  // requires the positive boundaries of phi to be communicated befor
  if (mu == _T) {
    // T dir
    for (i = 0; i < length; i++) {
      index1 = index_dir[i];
      index2 = neighbor[i];
      phi_pt = phi + 12 * index2;
      D_pt = D + 36 * index1 + 9 * _T;
      mvm_double(buffer1, D_pt, phi_pt);
      mvm_double(buffer1 + 3, D_pt, phi_pt + 3);
      mvm_double(buffer1 + 6, D_pt, phi_pt + 6);
      mvm_double(buffer1 + 9, D_pt, phi_pt + 9);
      eta1_pt = eta1 + 12 * index1;
      eta2_pt = eta2 + 12 * index1;
      twospin2_p_T_double(eta1_pt, eta2_pt, buffer1);
    }
  } else if (mu == _Z) {
    // Z dir
    for (i = 0; i < length; i++) {
      index1 = index_dir[i];
      index2 = neighbor[i];
      phi_pt = phi + 12 * index2;
      D_pt = D + 36 * index1 + 9 * _Z;
      mvm_double(buffer1, D_pt, phi_pt);
      mvm_double(buffer1 + 3, D_pt, phi_pt + 3);
      mvm_double(buffer1 + 6, D_pt, phi_pt + 6);
      mvm_double(buffer1 + 9, D_pt, phi_pt + 9);
      eta1_pt = eta1 + 12 * index1;
      eta2_pt = eta2 + 12 * index1;
      twospin2_p_Z_double(eta1_pt, eta2_pt, buffer1);
    }
  } else if (mu == _Y) {
    // Y dir
    for (i = 0; i < length; i++) {
      index1 = index_dir[i];
      index2 = neighbor[i];
      phi_pt = phi + 12 * index2;
      D_pt = D + 36 * index1 + 9 * _Y;
      mvm_double(buffer1, D_pt, phi_pt);
      mvm_double(buffer1 + 3, D_pt, phi_pt + 3);
      mvm_double(buffer1 + 6, D_pt, phi_pt + 6);
      mvm_double(buffer1 + 9, D_pt, phi_pt + 9);
      eta1_pt = eta1 + 12 * index1;
      eta2_pt = eta2 + 12 * index1;
      twospin2_p_Y_double(eta1_pt, eta2_pt, buffer1);
    }
  } else if (mu == _X) {
    // X dir
    for (i = 0; i < length; i++) {
      index1 = index_dir[i];
      index2 = neighbor[i];
      phi_pt = phi + 12 * index2;
      D_pt = D + 36 * index1 + 9 * _X;
      mvm_double(buffer1, D_pt, phi_pt);
      mvm_double(buffer1 + 3, D_pt, phi_pt + 3);
      mvm_double(buffer1 + 6, D_pt, phi_pt + 6);
      mvm_double(buffer1 + 9, D_pt, phi_pt + 9);
      eta1_pt = eta1 + 12 * index1;
      eta2_pt = eta2 + 12 * index1;
      twospin2_p_X_double(eta1_pt, eta2_pt, buffer1);
    }
  }
}

void apply_twisted_bc_to_vector_double(vector_double eta, vector_double phi, double *theta,
                                       Level *l) {
  int t, z, y, x, i;
  int *gl = l->global_lattice, sl[4];
  double phase[4];
  complex_double twisted_bc;
  for (i = 0; i < 4; i++)
    sl[i] = l->local_lattice[i] * l->global->my_coords[i];

  for (t = 0; t < l->local_lattice[0]; t++) {
    phase[_T] = theta[_T] * ((double)sl[_T] + t) / (double)gl[_T];
    for (z = 0; z < l->local_lattice[1]; z++) {
      phase[_Z] = phase[_T] + theta[_Z] * ((double)sl[_Z] + z) / (double)gl[_Z];
      for (y = 0; y < l->local_lattice[2]; y++) {
        phase[_Y] = phase[_Z] + theta[_Y] * ((double)sl[_Y] + y) / (double)gl[_Y];
        for (x = 0; x < l->local_lattice[3]; x++) {
          phase[_X] = phase[_Y] + theta[_X] * ((double)sl[_X] + x) / (double)gl[_X];
          twisted_bc = exp(I * phase[_X]);
#ifdef HAVE_TM1p1
          if (l->global->n_flavours == 2) {
            FOR24(*eta = (*phi) * twisted_bc; phi++; eta++;);
          } else
#endif
          {
            FOR12(*eta = (*phi) * twisted_bc; phi++; eta++;)
          }
        }
      }
    }
  }
}

void operator_updates_double(Level *l) {

  if (l->currentLevel > 0) {
    if (!l->idle) {
#ifdef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_double
      coarse_operator_double_setup_vectorized(l->is_double.operator, l);

#else

      coarse_operator_double_setup(l->is_double.interpolation, l);
#endif
      conf_double_gather(&(l->next_level->s_double.op), &(l->next_level->op_double), l->next_level);

      if (!l->next_level->idle && l->next_level->currentLevel > 0) {

        schwarz_double_boundary_update(&(l->next_level->s_double), l->next_level);

        if (l->global->method >= 4 && l->global->odd_even) {
          coarse_oddeven_setup_double(&(l->next_level->s_double.op), _REORDER, l->next_level);
        } else {
          coarse_operator_double_set_couplings(&(l->next_level->s_double.op), l->next_level);
        }
      }
      if (!l->next_level->idle && l->next_level->currentLevel == 0 && l->global->odd_even) {
        coarse_oddeven_setup_double(&(l->next_level->s_double.op), _NO_REORDERING, l->next_level);
      } else if (!l->next_level->idle && l->next_level->currentLevel == 0) {
        coarse_operator_double_set_couplings(&(l->next_level->s_double.op), l->next_level);
      }
      operator_updates_double(l->next_level);
    }
  }
}

void m0_update_double(double m0, operator_struct<double> *op, Level *l) {

  config_double clover = op->clover;

  if (clover != nullptr && op->m0 != m0) {
    int i, j;
    complex_double m0_diff = m0 - op->m0;

    op->m0 = m0;

    if (m0_diff != 0.0) {
      if (l->type == _FINE) {
        int start = 0;
        int n = l->num_inner_lattice_sites;
        clover += start * (l->global->csw ? 42 : 12);
        for (i = 0; i < n; i++) {
          for (j = 0; j < 12; j++) {
            clover[j] += m0_diff;
          }
          // clover term diag also stored as complex, so size is 2*15+2*6 = 42
          clover += (l->global->csw ? 42 : 12);
        }
      } else {
        int start = 0;
        int n = l->num_inner_lattice_sites;
        int k = l->num_parent_eig_vect;
        int sc_size = (l->num_parent_eig_vect) * (l->num_parent_eig_vect * 2 + 1);
        clover += start * sc_size;
        for (i = 0; i < n; i++) {
          for (j = 0; j < k; j++) {
            if (j > 0)
              clover += j + 1;
            *clover += m0_diff;
          }
          clover++;
          for (j = 0; j < k; j++) {
            if (j > 0)
              clover += j + 1;
            *clover += m0_diff;
          }
          clover += 1 + SQUARE(k);
        }
      }
    }
  }
}

void tm_term_double_setup(double mu, double even, double odd, operator_struct<double> *op,
                          Level *l) {

#ifdef HAVE_TM

  config_double tm_term = op->tm_term;
  if (tm_term != nullptr) {
    config_double odd_proj = op->odd_proj;
    complex_double shift = I * mu;
    complex_double even_shift = I * even;
    complex_double odd_shift = I * odd;

    op->mu = mu;
    op->mu_even_shift = even;
    op->mu_odd_shift = odd;

    int i, j;
    int n = l->num_inner_lattice_sites;

    if (l->type == _FINE) {
      complex_double tm_shift;

      for (i = 0; i < n; i++) {
        if (imag(even_shift) == 0. && imag(odd_shift) == 0.)
          tm_shift = shift;
        else
          tm_shift = shift + even_shift + odd_proj[0] * (odd_shift - even_shift);
        FOR6(*tm_term = -tm_shift; tm_term++;)
        FOR6(*tm_term = tm_shift; tm_term++;)
        odd_proj += 12;
      }
    } else {
      int k, m = l->num_parent_eig_vect;
      int tm_size = m * (m + 1);

      if (imag(even_shift) == 0. && imag(odd_shift) == 0.) {

        for (i = 0; i < n; i++) {
          for (j = 0; j < m; j++) {
            for (k = 0; k < j; k++)
              tm_term[k] = 0;
            tm_term += j;
            *tm_term = -1. * shift;
            tm_term++;
          }

          for (j = 0; j < m; j++) {
            for (k = 0; k < j; k++)
              tm_term[k] = 0;
            tm_term += j;
            *tm_term = shift;
            tm_term++;
          }
        }
      } else {
        complex_double odd_factor = odd_shift - even_shift;

        for (i = 0; i < n; i++) {
          for (j = 0; j < m; j++) {
            for (k = 0; k < j; k++)
              tm_term[k] = -1. * odd_factor * odd_proj[k];
            tm_term += j;
            odd_proj += j;
            *tm_term = -1. * (shift + even_shift + odd_factor * (*odd_proj));
            tm_term++;
            odd_proj++;
          }

          for (j = 0; j < m; j++) {
            for (k = 0; k < j; k++)
              tm_term[k] = odd_factor * odd_proj[k];
            tm_term += j;
            odd_proj += j;
            *tm_term = (shift + even_shift + odd_factor * (*odd_proj));
            tm_term++;
            odd_proj++;
          }
        }
      }
    }
  }
#endif
}

void epsbar_term_double_setup(double epsbar, double even, double odd, operator_struct<double> *op,
                              Level *l) {

#ifdef HAVE_TM1p1

  config_double eps_term = op->epsbar_term;
  if (eps_term != nullptr) {
    config_double odd_proj = op->odd_proj;
    complex_double shift = -epsbar;
    complex_double even_shift = I * even;
    complex_double odd_shift = I * odd;

    op->epsbar = epsbar;
    op->epsbar_ig5_even_shift = even;
    op->epsbar_ig5_odd_shift = odd;

    int i, j;
    int n = l->num_inner_lattice_sites;

    if (l->type == _FINE) {

      if (imag(even_shift) == 0. && imag(odd_shift) == 0.)
        for (i = 0; i < n; i++) {
          FOR12(*eps_term = shift; eps_term++;);
        }
      else
        for (i = 0; i < n; i++) {
          complex_double ig5_shift = even_shift + (*odd_proj) * (odd_shift - even_shift);
          FOR6(*eps_term = shift - ig5_shift; eps_term++;);
          FOR6(*eps_term = shift + ig5_shift; eps_term++;);
          odd_proj += 12;
        }
    } else {
      int k, m = l->num_parent_eig_vect;
      int eps_size = m * (m + 1);

      if (imag(even_shift) == 0. && imag(odd_shift) == 0.) {
        for (i = 0; i < 2 * n; i++) {
          for (j = 0; j < m; j++) {
            for (k = 0; k < j; k++)
              eps_term[k] = 0;
            eps_term += j;
            *eps_term = shift;
            eps_term++;
          }
        }
      } else {
        complex_double odd_factor = odd_shift - even_shift;

        for (i = 0; i < n; i++) {
          for (j = 0; j < m; j++) {
            for (k = 0; k < j; k++)
              eps_term[k] = -1. * odd_factor * odd_proj[k];
            eps_term += j;
            odd_proj += j;
            *eps_term = shift - (even_shift + odd_factor * (*odd_proj));
            eps_term++;
            odd_proj++;
          }

          for (j = 0; j < m; j++) {
            for (k = 0; k < j; k++)
              eps_term[k] = odd_factor * odd_proj[k];
            eps_term += j;
            odd_proj += j;
            *eps_term = shift + (even_shift + odd_factor * (*odd_proj));
            eps_term++;
            odd_proj++;
          }
        }
      }
    }
  }
#endif
}

void two_flavours_test_double(operator_struct<double> *op, Level *l) {

#ifdef HAVE_TM1p1
  double diff;

  vector_double vd1 = nullptr, vd2, vd3, vd4, vdd1, vdd2, vdd3, vdd4;
  vector_double vpp1 = nullptr, vpp2;

  ASSERT(l->global->n_flavours == 2);

  data_layout_n_flavours(1, l);

  int ivs = l->inner_vector_size;

  MALLOC(vd1, complex_double, 4 * ivs + 2 * 4 * ivs);
  MALLOC(vpp1, complex_double, 2 * 2 * ivs);

  vd2 = vd1 + ivs;
  vd3 = vd2 + ivs;
  vd4 = vd3 + ivs;
  vdd1 = vd4 + ivs;
  vdd2 = vdd1 + 2 * ivs;
  vdd3 = vdd2 + 2 * ivs;
  vdd4 = vdd3 + 2 * ivs;
  vpp2 = vpp1 + 2 * ivs;

  vector_double_define_random(vd1, 0, l->inner_vector_size, l);
  vector_double_define_random(vd2, 0, l->inner_vector_size, l);
  apply_operator_double(vd3, vd1, &(l->global->p), l);
#ifdef HAVE_TM
  vector_double_real_scale(l->global->op_double.tm_term, l->global->op_double.tm_term, -1, 0,
                           l->inner_vector_size, l);
#endif
  apply_operator_double(vd4, vd2, &(l->global->p), l);
#ifdef HAVE_TM
  vector_double_real_scale(l->global->op_double.tm_term, l->global->op_double.tm_term, -1, 0,
                           l->inner_vector_size, l);
#endif
  add_diagonal_double(vd3, vd2, l->global->op_double.epsbar_term, l->inner_vector_size);
  add_diagonal_double(vd4, vd1, l->global->op_double.epsbar_term, l->inner_vector_size);

  two_flavours_to_serial_double(vd1, vd2, vdd1, l);
  two_flavours_to_serial_double(vd3, vd4, vdd2, l);

  data_layout_n_flavours(2, l);

  trans_double(vpp1, vdd1, op->translation_table, l);
  apply_operator_double(vpp2, vpp1, &(l->p_double), l);
  trans_back_double(vdd3, vpp2, op->translation_table, l);

  vector_double_minus(vdd4, vdd3, vdd2, 0, l->inner_vector_size, l);
  diff = global_norm_double(vdd4, 0, l->inner_vector_size, l) /
         global_norm_double(vdd3, 0, l->inner_vector_size, l);

  test0_double("depth: %d, correctness of doublet Dirac operator double: %le\n", l->depth, diff);

  FREE(vd1, complex_double, 4 * ivs + 2 * 4 * ivs);
  FREE(vpp1, complex_double, 2 * 2 * ivs);

  if (l->global->method >= 4 && l->global->odd_even)
    oddeven_double_test(l);

#endif
}
