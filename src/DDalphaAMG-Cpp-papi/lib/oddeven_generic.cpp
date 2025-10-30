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

void selfcoupling_cholesky_decomposition_double(const config_double output, config_double input) {

  /*********************************************************************************
   * Performs a cholesky decomposition for a selfcoupling term.
   * Input = [ A 0 ]   , A=A*, B=B*
   *         [ 0 B ]
   * Output ordering: diag(A), diag(B), triu(A,1) row major, triu(B,1) row major
   *                  (matlab notation).
   *********************************************************************************/

  register int i, j, k;
  int n, offset[4] = {0, 12, 6, 27};
  config_double in_pt;
  config_double out_pt = output;
  complex_double L[6][6], s;

  for (n = 0; n < 2; n++) {
    // construct initial L = tril(A) for n=0, L = tril(B) for n=1, L row major
    in_pt = input + offset[2 * n];
    for (j = 0; j < 6; j++) {
      L[j][j] = (complex_double)*in_pt;
      in_pt++;
    }

    in_pt = input + offset[2 * n + 1];
    for (j = 0; j < 5; j++) {
      for (i = j + 1; i < 6; i++) {
        L[i][j] = (complex_double)conj(*in_pt);
        in_pt++;
      }
    }

    // calculate cholesky factor L
    for (i = 0; i < 6; i++) {
      for (j = 0; j <= i; j++) {
        s = L[i][j];
        for (k = 0; k < j; k++)
          s = s - L[i][k] * conj(L[j][k]);
        if (i > j)
          L[i][j] = s / L[j][j];
        else if (abs(L[i][i]) > EPS_double)
          L[i][i] = sqrt(s);
        else {
          L[i][i] = sqrt(s);
        }
      }
    }

    // output = tril(L) row major
    for (i = 0; i < 6; i++) {
      for (j = 0; j <= i; j++) {
        *out_pt = L[i][j];
        out_pt++;
      }
    }
  }
}

#ifdef HAVE_TM
void selfcoupling_LU_decomposition_double(const config_double output, config_double input) {

  /*********************************************************************************
   * Performs a LU decomposition for a selfcoupling term.
   * Input = [ A 0 ]   , A=A*, B=B* (diagonals excluded)
   *         [ 0 B ]
   * Input ordering: diag(A), diag(B), triu(A,1) row major, triu(B,1) row major
   *                  (matlab notation).
   * Output ordering: triu(L,1) + tril(U,0), i.e., output contains L and U without
   *                  the diagonal of L which is equal to 1
   *********************************************************************************/

  register int i, j, k;
  int n;
  config_double in_pt;
  config_double out_pt;

  int offset[4] = {0, 12, 6, 27};

  // construct initial L = A for n=0, L = B for n=1, L row major
  for (n = 0; n < 2; n++) {

    out_pt = output + n * 36;

    in_pt = input + offset[2 * n];
    for (j = 0; j < 6; j++) {
      out_pt[6 * j + j] = (complex_double)*in_pt;
      in_pt++;
    }

    in_pt = input + offset[2 * n + 1];
    for (j = 0; j < 5; j++) {
      for (i = j + 1; i < 6; i++) {
        out_pt[6 * j + i] = (complex_double)*in_pt;
        out_pt[6 * i + j] = (complex_double)conj(*in_pt);
        in_pt++;
      }
    }

    // calculate LU
    for (k = 0; k < 5; k++) {
      for (i = k + 1; i < 6; i++) {
        out_pt[6 * i + k] =
            out_pt[6 * i + k] / out_pt[6 * k + k]; // L: out(i,k) = out(i,k)/out(k,k)
        for (j = k + 1; j < 6; j++)
          out_pt[6 * i + j] =
              out_pt[6 * i + j] -
              out_pt[6 * i + k] * out_pt[6 * k + j]; // U: out(i,j) = out(i,j)-out(i,k)*out(k,j)
      }
    }
  }
}
#endif

#ifdef HAVE_TM1p1
void selfcoupling_LU_doublet_decomposition_double(const config_double output, config_double input) {

  /*********************************************************************************
   * Performs a LU decomposition for a selfcoupling term.
   * Input = [ A 0 ]   , A=A*, B=B* (diagonals excluded)
   *         [ 0 B ]
   * Input ordering: diag(A), diag(B), triu(A,1) row major, triu(B,1) row major
   *                  (matlab notation).
   * Output ordering: triu(L,1) + tril(U,0), i.e., output contains L and U without
   *                  the diagonal of L which is equal to 1
   *********************************************************************************/

  register int i, j, k;
  int n;
  config_double in_pt;
  config_double out_pt = output;

  int offset[8] = {0, 12, 24, 54, 6, 18, 39, 60};

  // construct initial L = A for n=0, L = B for n=1, L row major
  for (n = 0; n < 2; n++) {

    out_pt = output + n * 144;

    in_pt = input + offset[4 * n];
    for (j = 0; j < 6; j++) {
      out_pt[12 * j + j] = (complex_double)*in_pt;
      in_pt++;
    }
    in_pt = input + offset[4 * n + 1];
    for (j = 0; j < 6; j++) {
      out_pt[12 * (j + 6) + (j + 6)] = (complex_double)*in_pt;
      in_pt++;
    }

    in_pt = input + offset[4 * n + 2];
    for (j = 0; j < 5; j++) {
      for (i = j + 1; i < 6; i++) {
        out_pt[12 * (j + 6) + (i + 6)] = out_pt[12 * j + i] = (complex_double)*in_pt;
        out_pt[12 * (i + 6) + (j + 6)] = out_pt[12 * i + j] = (complex_double)conj(*in_pt);
        in_pt++;
      }
    }

    in_pt = input + offset[4 * n + 3];
    for (j = 0; j < 6; j++) {
      for (i = 0; i < 6; i++) {
        out_pt[12 * (j + 6) + i] = out_pt[12 * j + (i + 6)] = 0;
      }
    }
    for (j = 0; j < 6; j++) {
      out_pt[12 * (j + 6) + j] = out_pt[12 * j + (j + 6)] = (complex_double)*in_pt;
      in_pt++;
    }

    // calculate LU
    for (k = 0; k < 11; k++) {
      for (i = k + 1; i < 12; i++) {
        out_pt[12 * i + k] =
            out_pt[12 * i + k] / out_pt[12 * k + k]; // L: out(i,k) = out(i,k)/out(k,k)
        for (j = k + 1; j < 12; j++)
          out_pt[12 * i + j] =
              out_pt[12 * i + j] -
              out_pt[12 * i + k] * out_pt[12 * k + j]; // U: out(i,j) = out(i,j)-out(i,k)*out(k,j)
      }
    }
  }
}
#endif

static inline void LLH_perform_fwd_bwd_subs_double(vector_double x, vector_double b,
                                                   config_double L) {

  /*********************************************************************************
   * Solves L*(L^H)*x = b for x, i.e., the clover coupling for a single lattice
   * site.
   * - vector_double b: Right hand side.
   * - vector_double x: Solution.
   * - config_double L: Cholesky factor ( lower triangular matrix )
   *********************************************************************************/

  register int i, j;
  int n;

  for (n = 0; n < 2; n++) {
    // forward substitution with L
    for (i = 0; i < 6; i++) {
      x[i] = b[i];
      for (j = 0; j < i; j++) {
        x[i] = x[i] - *L * x[j];
        L++;
      }
      x[i] = x[i] / *L;
      L++;
    }
    L -= 21;
    // backward substitution with L^H
    for (i = 5; i >= 0; i--) {
      for (j = i + 1; j < 6; j++) {
        x[i] = x[i] - conj(L[(j * (j + 1)) / 2 + i]) * x[j];
      }
      x[i] = x[i] / conj(L[(i * (i + 1)) / 2 + i]);
    }
    x += 6;
    b += 6;
    L += 21;
  }
}

static inline void LU_perform_fwd_bwd_subs_double(vector_double x, vector_double b,
                                                  config_double LU) {

  /*********************************************************************************
   * Solves L*U*x = b for x, i.e., the clover coupling for a single lattice
   * site.
   * - vector_double b: Right hand side.
   * - vector_double x: Solution.
   * - config_double L: Lower matrix from modified LU decomposition
   * Note: U is given by u_{ii}=1, u_{ij}=l_{ji}* / l_{ii}
   *********************************************************************************/

  register int i, j, n;

#ifdef HAVE_TM1p1
  if (global->n_flavours == 2)
    for (n = 0; n < 2; n++) {
      // solve x = U^(-1) L^(-1) b
      // forward substitution with L
      for (i = 0; i < 12; i++) {
        x[i] = b[i];
        for (j = 0; j < i; j++) {
          x[i] = x[i] - LU[i * 12 + j] * x[j];
        }
      }
      // backward substitution with U
      for (i = 12 - 1; i >= 0; i--) {
        for (j = i + 1; j < 12; j++) {
          x[i] = x[i] - LU[i * 12 + j] * x[j];
        }
        x[i] = x[i] / LU[i * (12 + 1)];
      }
      x += 12;
      b += 12;
      LU += 12 * 12;
    }
  else
#endif
    for (n = 0; n < 2; n++) {
      // solve x = U^(-1) L^(-1) b
      // forward substitution with L
      for (i = 0; i < 6; i++) {
        x[i] = b[i];
        for (j = 0; j < i; j++) {
          x[i] = x[i] - LU[i * 6 + j] * x[j];
        }
      }
      // backward substitution with U
      for (i = 6 - 1; i >= 0; i--) {
        for (j = i + 1; j < 6; j++) {
          x[i] = x[i] - LU[i * 6 + j] * x[j];
        }
        x[i] = x[i] / LU[i * (6 + 1)];
      }
      x += 6;
      b += 6;
      LU += 6 * 6;
    }
}

static inline void LLH_multiply_double(vector_double y, vector_double x, config_double L) {

  /*********************************************************************************
   * Applies the clover coupling term to a vector, by multiplying L^H
   * and then L.
   * - vector_double x: Input vector.
   * - vector_double y: Output vector.
   * - config_double L: Cholesky factor ( lower triangular matrix )
   *********************************************************************************/

  register int i, j;
  int n;
  complex_double z[6];

  for (n = 0; n < 2; n++) {
    // z = L^H x
    for (j = 0; j < 6; j++) {   // columns
      for (i = 0; i < j; i++) { // rows
        z[i] += conj(*L) * x[j];
        L++;
      }
      z[j] = conj(*L) * x[j];
      L++;
    }
    L -= 21;
    // y = L*z;
    for (i = 0; i < 6; i++) { // rows
      y[i] = *L * z[0];
      L++;
      for (j = 1; j <= i; j++) { // columns
        y[i] += *L * z[j];
        L++;
      }
    }
    x += 6;
    y += 6;
  }
}

static inline void LU_multiply_double(vector_double y, vector_double x, config_double LU) {

  /*********************************************************************************
   * Applies the clover coupling term to a vector, by multiplying L^H
   * and then L.
   * - vector_double x: Input vector.
   * - vector_double y: Output vector.
   * - config_double LU: LU decomposition
   *********************************************************************************/

  register int i, j, n;

#ifdef HAVE_TM1p1
  if (global->n_flavours == 2)
    for (n = 0; n < 2; n++) {
      for (i = 0; i < 12; i++) {
        y[i] = LU[i * (12 + 1)] * x[i];
        for (j = i + 1; j < 12; j++)
          y[i] += LU[i * 12 + j] * x[j];
      }
      // multiplication with L
      for (i = 12 - 1; i > 0; i--)
        for (j = 0; j < i; j++)
          y[i] += LU[i * 12 + j] * y[j];

      x += 12;
      y += 12;
      LU += 12 * 12;
    }
  else
#endif
    for (n = 0; n < 2; n++) {
      for (i = 0; i < 6; i++) {
        y[i] = LU[i * (6 + 1)] * x[i];
        for (j = i + 1; j < 6; j++)
          y[i] += LU[i * 6 + j] * x[j];
      }
      // multiplication with L
      for (i = 6 - 1; i > 0; i--)
        for (j = 0; j < i; j++)
          y[i] += LU[i * 6 + j] * y[j];

      x += 6;
      y += 6;
      LU += 6 * 6;
    }
}

void diag_ee_double(vector_double y, vector_double x, operator_struct<double> *op, Level *l,
                    int start, int end) {

  /*********************************************************************************
   * Applies the even-even block of the odd even decomposition to a vector.
   * - vector_double x: Input vector.
   * - vector_double y: Output vector.
   *********************************************************************************/

  Global *global = l->global;

#ifdef HAVE_TM1p1
  if (global->n_flavours == 2) {
    x += start;
    y += start;
#ifdef OPTIMIZED_SELF_COUPLING_double
    double *sc_pt = op->clover_doublet_vectorized + (start / 24) * 288;
    double *x_pt = (double *)x;
    double *y_pt = (double *)y;
    for (int i = start; i < end; i += 24) {
      sse_site_clover_double(y_pt, x_pt, sc_pt);
      y_pt += 2 * 24;
      x_pt += 2 * 24;
      sc_pt += 288;
    }
    config_double epsbar_term = op->epsbar_term + (start / 24) * 12;
    if (global->n_flavours == 2 &&
        (op->epsbar != 0 || op->epsbar_ig5_odd_shift != 0 || op->epsbar_ig5_odd_shift != 0))
      apply_doublet_coupling_double(x, y, epsbar_term, end - start);
#else
    config_double sc = op->clover_doublet_oo_inv + (start / 24) * 288;
    // diagonal blocks applied to the even sites
    for (int i = start; i < end; i += 24) {
      LU_multiply_double(y, x, sc);
      y += 24;
      x += 24;
      sc += 288;
    }
#endif
  } else {
#endif
    x += start;
    y += start;
    if (global->csw) {
#ifdef OPTIMIZED_SELF_COUPLING_double
      double *sc_pt = op->clover_vectorized + (start / 12) * 144;
      double *x_pt = (double *)x;
      double *y_pt = (double *)y;
      for (int i = start; i < end; i += 12) {
        sse_site_clover_double(y_pt, x_pt, sc_pt);
        y_pt += 2 * 12;
        x_pt += 2 * 12;
        sc_pt += 144;
      }
#elif defined(HAVE_TM)
    config_double sc = op->clover + (start / 12) * 72;
    // diagonal blocks applied to the even sites
    for (int i = start; i < end; i += 12) {
      LU_multiply_double(y, x, sc);
      y += 12;
      x += 12;
      sc += 72;
    }
#else
    config_double sc = op->clover + (start / 12) * 42;
    // diagonal blocks applied to the even sites
    for (int i = start; i < end; i += 12) {
      LLH_multiply_double(y, x, sc);
      y += 12;
      x += 12;
      sc += 42;
    }
#endif
    } else {
      config_double sc = op->clover + start;
      for (int i = start; i < end; i += 12) {
        FOR12(*y = (*x) * (*sc); y++; x++; sc++;)
      }
    }
#ifdef HAVE_TM1p1
  }
#endif
}

// for debugging only
void diag_ee_inv_double(vector_double y, vector_double x, operator_struct<double> *op, Level *l) {

  Global *global = l->global;

#ifdef HAVE_TM1p1
  if (global->n_flavours == 2) {
    int i, n1 = op->num_even_sites;
    config_double sc = op->clover_doublet_oo_inv;
    // diagonal blocks applied to the even sites
    for (i = 0; i < n1; i++) {
      LU_perform_fwd_bwd_subs_double(y, x, sc);
      y += 24;
      x += 24;
      sc += 288;
    }
  } else {
#endif
    int i, n1 = op->num_even_sites;
    config_double sc = op->clover;
    if (global->csw) {
      // diagonal blocks applied to the even sites
      for (i = 0; i < n1; i++) {
#ifndef HAVE_TM
        LLH_perform_fwd_bwd_subs_double(y, x, sc);
        y += 12;
        x += 12;
        sc += 42;
#else
      LU_perform_fwd_bwd_subs_double(y, x, sc);
      y += 12;
      x += 12;
      sc += 72;
#endif
      }
    } else {
      for (i = 0; i < n1; i++) {
        FOR12(*y = (*x) / (*sc); y++; x++; sc++;)
      }
    }
#ifdef HAVE_TM1p1
  }
#endif
}

// for debugging only
void diag_oo_double(vector_double y, vector_double x, operator_struct<double> *op, Level *l) {

  /*********************************************************************************
   * Applies the odd-odd block of the odd even decomposition to a vector.
   * - vector_double x: Input vector.
   * - vector_double y: Output vector.
   *********************************************************************************/

  Global *global = l->global;

#ifdef HAVE_TM1p1
  if (global->n_flavours == 2) {
    int i, n1 = op->num_even_sites, n2 = op->num_odd_sites;
    config_double sc = op->clover_doublet_oo_inv + n1 * 288;
    x += n1 * 24;
    y += n1 * 24;
    // diagonal blocks applied to the even sites
    for (i = 0; i < n2; i++) {
      LU_multiply_double(y, x, sc);
      y += 24;
      x += 24;
      sc += 288;
    }
  } else {
#endif
    int i, n1 = op->num_even_sites, n2 = op->num_odd_sites;
    config_double sc = op->clover;
    x += n1 * 12;
    y += n1 * 12;
    // diagonal blocks applied to the odd sites
    if (global->csw) {
#ifndef HAVE_TM
      sc += n1 * 42;
      for (i = 0; i < n2; i++) {
        LLH_multiply_double(y, x, sc);
        y += 12;
        x += 12;
        sc += 42;
      }
#else
    sc += n1 * 72;
    for (i = 0; i < n2; i++) {
      LU_multiply_double(y, x, sc);
      y += 12;
      x += 12;
      sc += 72;
    }
#endif
    } else {
      sc += n1 * 12;
      for (i = 0; i < n2; i++) {
        FOR12(*y = (*x) * (*sc); y++; x++; sc++;)
      }
    }
#ifdef HAVE_TM1p1
  }
#endif
}

void diag_oo_inv_double(vector_double y, vector_double x, operator_struct<double> *op, Level *l,
                        int start, int end) {

  Global *global = l->global;

#ifdef HAVE_TM1p1
  if (global->n_flavours == 2) {
    x += start;
    y += start;
    // inverse diagonal blocks applied to the odd sites
#ifdef OPTIMIZED_SELF_COUPLING_double
    double *sc_pt = op->clover_doublet_oo_inv_vectorized + (start / 24) * 2 * 288;
    double *x_pt = (double *)x;
    double *y_pt = (double *)y;
    for (int i = start; i < end; i += 24) {
      sse_site_clover_doublet_double(y_pt, x_pt, sc_pt);
      y_pt += 2 * 24;
      x_pt += 2 * 24;
      sc_pt += 2 * 288;
    }
#else
    config_double sc = op->clover_doublet_oo_inv + (start / 24) * 288;
    for (int i = start; i < end; i += 24) {
      LU_perform_fwd_bwd_subs_double(y, x, sc);
      y += 24;
      x += 24;
      sc += 288;
    }
#endif
  } else {
#endif
    config_double sc = op->clover;
    x += start;
    y += start;
    // inverse diagonal blocks applied to the odd sites
    if (global->csw) {
#ifdef OPTIMIZED_SELF_COUPLING_double
      double *sc_pt = op->clover_vectorized + 2 * 2 * (3 * start);
      double *x_pt = (double *)x;
      double *y_pt = (double *)y;
      for (int i = start; i < end; i += 12) {
        sse_site_clover_double(y_pt, x_pt, sc_pt);
        y_pt += 2 * 12;
        x_pt += 2 * 12;
        sc_pt += 2 * 2 * 36;
      }
#elif defined(HAVE_TM)
    sc += (start / 12) * 72;
    for (int i = start; i < end; i += 12) {
      LU_perform_fwd_bwd_subs_double(y, x, sc);
      y += 12;
      x += 12;
      sc += 72;
    }
#else
    sc += (start / 12) * 42;
    for (int i = start; i < end; i += 12) {
      LLH_perform_fwd_bwd_subs_double(y, x, sc);
      y += 12;
      x += 12;
      sc += 42;
    }
#endif
    } else {
      sc += start;
      for (int i = start; i < end; i += 12) {
        FOR12(*y = (*x) / (*sc); y++; x++; sc++;)
      }
    }
#ifdef HAVE_TM1p1
  }
#endif
}

void oddeven_setup_double(operator_struct<double> *in, Level *l) {

  /*********************************************************************************
   * Reorder data layouts and index tables to allow for odd even preconditioning.
   *********************************************************************************/

  Global *global = l->global;
  int j, k, k_e, k_o, n = l->num_inner_lattice_sites, oe_offset = 0, mu, nu,
                      sc_size = global->csw ? 42 : 12, lu_dec_size = 42, bs, **bt = nullptr,
                      *eot = nullptr, *nt = nullptr, *tt = nullptr, t, z, y, x, le[4], N[4];
  config_double sc_in = in->clover, nc_in = in->D;
  config_double Aee = nullptr, Aoo = nullptr;
  operator_struct<double> *op = &(l->oe_op_double);

  op->m0 = in->m0;

#ifdef HAVE_TM
  op->mu = in->mu;
  op->mu_even_shift = in->mu_even_shift;
  op->mu_odd_shift = in->mu_odd_shift;

  lu_dec_size = 72;
  config_double tm_term_in = in->tm_term;
#endif

  for (mu = 0; mu < 4; mu++) {
    le[mu] = l->local_lattice[mu];
    N[mu] = le[mu] + 1;
    op->table_dim[mu] = N[mu];
  }

  for (mu = 0; mu < 4; mu++)
    oe_offset += (l->local_lattice[mu] * (global->my_coords[mu] / l->comm_offset[mu])) % 2;
  oe_offset = oe_offset % 2;

  // estimate site numbers
  op->num_even_sites = 0;
  op->num_odd_sites = 0;
  op->oe_offset = oe_offset;

  for (t = 0; t < le[_T]; t++)
    for (z = 0; z < le[_Z]; z++)
      for (y = 0; y < le[_Y]; y++)
        for (x = 0; x < le[_X]; x++) {
          if ((t + z + y + x + oe_offset) % 2 == 1) {
            op->num_odd_sites++;
          } else {
            op->num_even_sites++;
          }
        }

  // re-order clover term (i.e., self coupling)
  if (global->csw) {
    MALLOC(op->clover, complex_double, lu_dec_size * n);
    Aee = op->clover;
    Aoo = op->clover + op->num_even_sites * lu_dec_size;
    /* TODO: fix the vectorized part
#ifdef OPTIMIZED_SELF_COUPLING_double
MALLOC_HUGEPAGES( op->clover_vectorized, double, l->num_inner_lattice_sites*2*2*36,
4*SIMD_LENGTH_double ); double *Aee_vectorized = op->clover_vectorized; double *Aoo_vectorized =
op->clover_vectorized + op->num_even_sites*2*2*36; #endif
    */
    for (t = 0; t < le[_T]; t++)
      for (z = 0; z < le[_Z]; z++)
        for (y = 0; y < le[_Y]; y++)
          for (x = 0; x < le[_X]; x++) {
            if ((t + z + y + x + oe_offset) % 2 == 1) {
              // odd sites
              /* TODO: fix the vectorized part
#ifdef OPTIMIZED_SELF_COUPLING_double
              double tmp[144] __attribute__((aligned(64)));
              sse_set_clover_double( tmp, sc_in );
#ifdef HAVE_TM
              if (global->mu + global->mu_odd_shift != 0.0 || global->mu + global->mu_even_shift !=
0.0 ) sse_add_diagonal_clover_double( tmp, tm_term_in ); #endif sse_site_clover_invert_double( tmp,
Aoo_vectorized ); Aoo_vectorized += 2*2*36; #endif
              */
#ifndef HAVE_TM
              selfcoupling_cholesky_decomposition_double(Aoo, sc_in);
#else
              complex_double buffer[42];
              for (int i = 0; i < 42; i++)
                buffer[i] = sc_in[i];
              for (int i = 0; i < 12; i++)
                buffer[i] += tm_term_in[i];
              selfcoupling_LU_decomposition_double(Aoo, buffer);
#endif
              Aoo += lu_dec_size;
            } else {
              // even sites
              /*
#ifdef OPTIMIZED_SELF_COUPLING_double
              sse_set_clover_double( Aee_vectorized, sc_in );
#ifdef HAVE_TM
              sse_add_diagonal_clover_double( Aee_vectorized, tm_term_in );
#endif
              Aee_vectorized += 2*2*36;
#endif
              */
#ifndef HAVE_TM
              selfcoupling_cholesky_decomposition_double(Aee, sc_in);
#else
              complex_double buffer[42];
              for (int i = 0; i < 42; i++)
                buffer[i] = (complex_double)sc_in[i];
              for (int i = 0; i < 12; i++)
                buffer[i] += (complex_double)tm_term_in[i];
              selfcoupling_LU_decomposition_double(Aee, buffer);
#endif
              Aee += lu_dec_size;
            }
            sc_in += sc_size;
#ifdef HAVE_TM
            tm_term_in += 12;
#endif
          }
  } else {
    MALLOC(op->clover, complex_double, 12 * n);
    Aee = op->clover;
    Aoo = op->clover + op->num_even_sites * 12;

    for (t = 0; t < le[_T]; t++)
      for (z = 0; z < le[_Z]; z++)
        for (y = 0; y < le[_Y]; y++)
          for (x = 0; x < le[_X]; x++) {
            if ((t + z + y + x + oe_offset) % 2 == 1) {
              // odd sites
              FOR12(*Aoo = *sc_in; Aoo++; sc_in++;)
            } else {
              // even sites
              FOR12(*Aee = *sc_in; Aee++; sc_in++;)
            }
          }
  }

#ifdef HAVE_TM1p1

  int lu_doublet_dec_size = 288;
  config_double eps_term_in = in->epsbar_term;
  sc_in = in->clover;
#ifdef HAVE_TM
  tm_term_in = in->tm_term;
#endif
  op->epsbar = in->epsbar;
  op->epsbar_ig5_even_shift = in->epsbar_ig5_even_shift;
  op->epsbar_ig5_odd_shift = in->epsbar_ig5_odd_shift;

  // re-order clover term (i.e., self coupling)
  MALLOC(op->clover_doublet_oo_inv, complex_double, lu_doublet_dec_size * n);
  Aee = op->clover_doublet_oo_inv;
  Aoo = op->clover_doublet_oo_inv + op->num_even_sites * lu_doublet_dec_size;
  /*
#ifdef OPTIMIZED_SELF_COUPLING_double
  MALLOC_HUGEPAGES( op->clover_doublet_vectorized, double, l->num_inner_lattice_sites*2*4*36,
4*SIMD_LENGTH_double ); MALLOC_HUGEPAGES( op->clover_doublet_oo_inv_vectorized, double,
op->num_odd_sites*2*2*144, 4*SIMD_LENGTH_double ); double *Aee_vectorized =
op->clover_doublet_vectorized; double *Aoo_vectorized = op->clover_doublet_vectorized +
op->num_even_sites*288; double *Aoo_inverse_vectorized = op->clover_doublet_oo_inv_vectorized;
#endif
  */
  for (t = 0; t < le[_T]; t++)
    for (z = 0; z < le[_Z]; z++)
      for (y = 0; y < le[_Y]; y++)
        for (x = 0; x < le[_X]; x++) {
          if ((t + z + y + x + oe_offset) % 2 == 1) {
            // odd sites
            /*
#ifdef OPTIMIZED_SELF_COUPLING_double
            sse_set_clover_doublet_double( Aoo_vectorized, sc_in );
#ifdef HAVE_TM
            if (global->mu + global->mu_odd_shift != 0.0 || global->mu + global->mu_even_shift !=
0.0 ) sse_add_diagonal_clover_doublet_double( Aoo_vectorized, tm_term_in ); #endif complex_double
eps_term[12]; for(int i=0; i<12; i++) { eps_term[i] = eps_term_in[i];
            }
            sse_site_clover_doublet_invert_double( Aoo_vectorized, (config_double) eps_term,
Aoo_inverse_vectorized ); Aoo_vectorized += 288; Aoo_inverse_vectorized += 2*288; #endif
            */
            complex_double buffer[66];
            if (global->csw) {
              for (int i = 0; i < 12; i++) // 0-23
                buffer[i + 12] = buffer[i] = (complex_double)sc_in[i];
              for (int i = 12; i < 42; i++) // 24-53
                buffer[i + 12] = (complex_double)sc_in[i];
            } else {
              for (int i = 0; i < 12; i++) // 0-23
                buffer[i + 12] = buffer[i] = (complex_double)sc_in[i];
              for (int i = 12; i < 42; i++) // 24-53
                buffer[i + 12] = 0;
            }
            for (int i = 0; i < 12; i++) // 54-65
              buffer[i + 54] = (complex_double)eps_term_in[i];
#ifdef HAVE_TM
            if (global->mu + global->mu_odd_shift != 0.0 ||
                global->mu + global->mu_even_shift != 0.0)
              for (int i = 0; i < 12; i++) {
                buffer[i] += (complex_double)tm_term_in[i];
                buffer[i + 12] -= (complex_double)tm_term_in[i];
              }
#endif
            selfcoupling_LU_doublet_decomposition_double(Aoo, buffer);
            Aoo += lu_doublet_dec_size;
          } else {
            // even sites
            /*
#ifdef OPTIMIZED_SELF_COUPLING_double
            sse_set_clover_doublet_double( Aee_vectorized, sc_in );
#ifdef HAVE_TM
            sse_add_diagonal_clover_doublet_double( Aee_vectorized, tm_term_in );
#endif
            Aee_vectorized += 288;
#endif
            */
            complex_double buffer[66];
            if (global->csw) {
              for (int i = 0; i < 12; i++) // 0-23
                buffer[i + 12] = buffer[i] = (complex_double)sc_in[i];
              for (int i = 12; i < 42; i++) // 24-53
                buffer[i + 12] = (complex_double)sc_in[i];
            } else {
              for (int i = 0; i < 12; i++) // 0-23
                buffer[i + 12] = buffer[i] = (complex_double)sc_in[i];
              for (int i = 12; i < 42; i++) // 24-53
                buffer[i + 12] = 0;
            }
            for (int i = 0; i < 12; i++) // 54-65
              buffer[i + 54] = (complex_double)eps_term_in[i];
#ifdef HAVE_TM
            if (global->mu + global->mu_odd_shift != 0.0 ||
                global->mu + global->mu_even_shift != 0.0)
              for (int i = 0; i < 12; i++) {
                buffer[i] += (complex_double)tm_term_in[i];
                buffer[i + 12] -= (complex_double)tm_term_in[i];
              }
#endif
            selfcoupling_LU_doublet_decomposition_double(Aee, buffer);
            Aee += lu_doublet_dec_size;
          }
          sc_in += sc_size;
          eps_term_in += 12;
#ifdef HAVE_TM
          tm_term_in += 12;
#endif
        }
#endif

  // re-order hopping term (i.e., nearest neighbor coupling)
  MALLOC(op->D, complex_double, 36 * n);

  k = 0;
  k_e = 0;
  k_o = 0;
  for (t = 0; t < le[_T]; t++)
    for (z = 0; z < le[_Z]; z++)
      for (y = 0; y < le[_Y]; y++)
        for (x = 0; x < le[_X]; x++) {
          if ((t + z + y + x + oe_offset) % 2 == 1) {
            for (j = 0; j < 36; j++)
              (op->D + (k_e + op->num_even_sites) * 36)[j] = (complex_double)(nc_in + k * 36)[j];
            k_e++;
          } else {
            for (j = 0; j < 36; j++)
              (op->D + k_o * 36)[j] = (complex_double)(nc_in + k * 36)[j];
            k_o++;
          }
          k++;
        }

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double // D_vectorized just used in the double environment
  MALLOC_HUGEPAGES(op->D_vectorized, double, 2 * 4 * l->inner_vector_size, 4 * SIMD_LENGTH_double);
  MALLOC_HUGEPAGES(op->D_transformed_vectorized, double, 2 * 4 * l->inner_vector_size,
                   4 * SIMD_LENGTH_double);
  for (int i = 0; i < l->num_inner_lattice_sites; i++) {
    double *D_vectorized = op->D_vectorized + 96 * i;
    double *D_transformed_vectorized = op->D_transformed_vectorized + 96 * i;
    complex_double *D_out_pt = op->D + 36 * i;
    for (int mu = 0; mu < 4; mu++) {
      set_double_D_vectorized(D_vectorized + 24 * mu, D_transformed_vectorized + 24 * mu,
                              D_out_pt + 9 * mu);
    }
  }
#endif

  // define data layout
  MALLOC(op->index_table, int, N[_T] * N[_Z] * N[_Y] * N[_X]);
  eot = op->index_table;

  define_eot(eot, N, l);

  // neighbor table, translation table
  MALLOC(op->neighbor_table, int, 5 * N[_T] * N[_Z] * N[_Y] * N[_X]);
  MALLOC(op->backward_neighbor_table, int, 5 * N[_T] * N[_Z] * N[_Y] * N[_X]);
  MALLOC(op->translation_table, int, le[_T] * le[_Z] * le[_Y] * le[_X]);
  nt = op->neighbor_table;
  tt = op->translation_table;

  define_nt_bt_tt(nt, op->backward_neighbor_table, nullptr, tt, eot, N, l);

  // boundary table
  for (mu = 0; mu < 4; mu++) {
    bs = 1;
    le[mu] = 1;
    for (nu = 0; nu < 4; nu++)
      bs *= le[nu];

    MALLOC(op->c.boundary_table[2 * mu], int, bs);
    op->c.boundary_table[2 * mu + 1] = op->c.boundary_table[2 * mu];

    le[mu] = l->local_lattice[mu];
  }

  bt = op->c.boundary_table;
  define_eo_bt(bt, eot, op->c.num_even_boundary_sites, op->c.num_odd_boundary_sites,
               op->c.num_boundary_sites, N, l);

  j = (l->num_lattice_site_var / 2) * l->num_lattice_sites;
#ifdef HAVE_TM1p1
  j *= 2;
#endif

  MALLOC(op->buffer, complex_double *, 2);
  op->buffer[0] = nullptr;
#ifdef HAVE_TM1p1
  MALLOC(op->buffer[0], complex_double, 4 * l->vector_size);
  op->buffer[1] = op->buffer[0] + 2 * l->vector_size;
#else
  MALLOC(op->buffer[0], complex_double, 2 * l->vector_size);
  op->buffer[1] = op->buffer[0] + l->vector_size;
#endif
  ghost_alloc_double(0, &(op->c), l);
  ghost_sendrecv_init_double(_COARSE_GLOBAL, &(op->c), l);
  l->gmresSmoother.vectorEnd = op->num_even_sites * l->num_lattice_site_var;
}

void oddeven_free_double(Level *l) {

  Global *global = l->global;
  int mu, nu, *ll = l->local_lattice, bs;
#ifdef HAVE_TM
  lu_dec_size = 72;
#endif

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
  FREE_HUGEPAGES(l->oe_op_double.D_vectorized, double, 2 * 4 * l->inner_vector_size);
  FREE_HUGEPAGES(l->oe_op_double.D_transformed_vectorized, double, 2 * 4 * l->inner_vector_size);
#endif
#ifdef OPTIMIZED_SELF_COUPLING_double
  FREE_HUGEPAGES(l->oe_op_double.clover_vectorized, double,
                 l->num_inner_lattice_sites * 2 * 2 * 36);
#ifdef HAVE_TM1p1
  FREE_HUGEPAGES(l->oe_op_double.clover_doublet_vectorized, double,
                 l->num_inner_lattice_sites * 2 * 4 * 36);
  FREE_HUGEPAGES(l->oe_op_double.clover_doublet_oo_inv_vectorized, double,
                 l->num_inner_lattice_sites * 2 * 2 * 144);
#endif
#endif

  ghost_free_double(&(l->oe_op_double.c), l);
  FREE(l->oe_op_double.D, complex_double, 4 * 9 * l->num_inner_lattice_sites);
  if (global->csw)
    FREE(l->oe_op_double.clover, complex_double, 42 * l->num_inner_lattice_sites);
  else
    FREE(l->oe_op_double.clover, complex_double, 12 * l->num_inner_lattice_sites);
  FREE(l->oe_op_double.index_table, int, (ll[_T] + 1) * (ll[_Z] + 1) * (ll[_Y] + 1) * (ll[_X] + 1));
  FREE(l->oe_op_double.neighbor_table, int,
       5 * (ll[_T] + 1) * (ll[_Z] + 1) * (ll[_Y] + 1) * (ll[_X] + 1));
  FREE(l->oe_op_double.backward_neighbor_table, int,
       5 * (ll[_T] + 1) * (ll[_Z] + 1) * (ll[_Y] + 1) * (ll[_X] + 1));
  FREE(l->oe_op_double.translation_table, int, ll[_T] * ll[_Z] * ll[_Y] * ll[_X]);

  for (mu = 0; mu < 4; mu++) {
    bs = 1;
    for (nu = 0; nu < 4; nu++)
      if (mu != nu)
        bs *= ll[nu];

    FREE(l->oe_op_double.c.boundary_table[2 * mu], int, bs);
    l->oe_op_double.c.boundary_table[2 * mu + 1] = nullptr;
  }

#ifdef HAVE_TM1p1
  FREE(l->oe_op_double.buffer[0], complex_double, 4 * l->vector_size);
#else
  FREE(l->oe_op_double.buffer[0], complex_double, 2 * l->vector_size);
#endif
  FREE(l->oe_op_double.buffer, complex_double *, 2);
#ifdef HAVE_TM1p1
  FREE(l->oe_op_double.prnT, complex_double,
       2 * (l->num_lattice_site_var / 2) * l->num_lattice_sites * 8);
  FREE(l->oe_op_double.clover_doublet_oo_inv, complex_double, 288 * l->num_inner_lattice_sites);
#else

#endif
}

void oddeven_to_serial_double(vector_double out, vector_double in, Level *l) {

  /*********************************************************************************
   * Translates a vector from an odd even double precision layout to a serial
   * double precision layout.
   *********************************************************************************/

  int k, nsv = l->num_lattice_site_var;

  for (int i = 0; i < l->num_inner_lattice_sites; i++) {
    k = l->oe_op_double.translation_table[i];
    for (int j = 0; j < nsv; j++) {
      out[i * nsv + j] = (complex_double)in[k * nsv + j];
    }
  }
}

void serial_to_oddeven_double(vector_double out, vector_double in, Level *l) {

  /*********************************************************************************
   * Translates a vector from a serial double precision layout to an odd even
   * double precision layout.
   *********************************************************************************/

  int k, nsv = l->num_lattice_site_var;

  for (int i = 0; i < l->num_inner_lattice_sites; i++) {
    k = l->oe_op_double.translation_table[i];
    for (int j = 0; j < nsv; j++) {
      out[k * nsv + j] = (complex_double)in[i * nsv + j];
    }
  }
}

void oddeven_to_block_double(vector_double out, vector_double in, Level *l) {

  int k, m, nsv = l->num_lattice_site_var;

  for (int i = 0; i < l->num_inner_lattice_sites; i++) {
    k = l->oe_op_double.translation_table[i];
    m = l->s_double.op.translation_table[i];
    for (int j = 0; j < nsv; j++) {
      out[m * nsv + j] = in[k * nsv + j];
    }
  }
}

void block_to_oddeven_double(vector_double out, vector_double in, Level *l) {

  int k, m, nsv = l->num_lattice_site_var;

  for (int i = 0; i < l->num_inner_lattice_sites; i++) {
    k = l->oe_op_double.translation_table[i];
    m = l->s_double.op.translation_table[i];
    for (int j = 0; j < nsv; j++) {
      out[k * nsv + j] = in[m * nsv + j];
    }
  }
}

void hopping_term_double(vector_double eta, vector_double phi, operator_struct<double> *op,
                         const int amount, Level *l) {

  int start_even, end_even, start_odd, end_odd,
      n = l->num_inner_lattice_sites, *neighbor = op->neighbor_table, start = 0,
      plus_dir_param = _FULL_SYSTEM, minus_dir_param = _FULL_SYSTEM;

  auto *prnT = new complex_double[(l->num_lattice_site_var)/2 * l->num_lattice_sites];
  auto *prnZ = new complex_double[(l->num_lattice_site_var)/2 * l->num_lattice_sites];
  auto *prnY = new complex_double[(l->num_lattice_site_var)/2 * l->num_lattice_sites];
  auto *prnX = new complex_double[(l->num_lattice_site_var)/2 * l->num_lattice_sites];

  auto *prpT = new complex_double[(l->num_lattice_site_var)/2 * l->num_lattice_sites];
  auto *prpZ = new complex_double[(l->num_lattice_site_var)/2 * l->num_lattice_sites];
  auto *prpY = new complex_double[(l->num_lattice_site_var)/2 * l->num_lattice_sites];
  auto *prpX = new complex_double[(l->num_lattice_site_var)/2 * l->num_lattice_sites];

  if (amount == _EVEN_SITES || amount == _ODD_SITES) {
    start_even = 0;
    end_even = op->num_even_sites;
    start_odd = op->num_even_sites;
    end_odd = op->num_even_sites + op->num_odd_sites;
  }

  if (amount == _EVEN_SITES) {
    start = start_odd, n = end_odd;
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if (amount == _ODD_SITES) {
    start = start_even, n = end_even;
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
  complex_double *prn[4] = {op->prnT, op->prnZ, op->prnY, op->prnX};
  complex_double *prp[4] = {op->prpT, op->prpZ, op->prpY, op->prpX};
#else
  int i, *nb_pt;
  vector_double phi_pt, eta_pt, end_pt;
  config_double D_pt;
#endif

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2) {
    // project in negative directions
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    dprp_double(prn, phi, 24 * start, 24 * n);
#else
    complex_double pbuf[12];
    for (i = 12 * start, phi_pt = phi + 24 * start; i < 12 * n; i += 12, phi_pt += 24) {
      dprp_T_double(op->prnT + i, phi_pt);
      dprp_Z_double(op->prnZ + i, phi_pt);
      dprp_Y_double(op->prnY + i, phi_pt);
      dprp_X_double(op->prnX + i, phi_pt);
    }
#endif
    // start communication in negative direction
    ghost_sendrecv_double(op->prnT, _T, -1, &(op->c), minus_dir_param, l);
    ghost_sendrecv_double(op->prnZ, _Z, -1, &(op->c), minus_dir_param, l);
    ghost_sendrecv_double(op->prnY, _Y, -1, &(op->c), minus_dir_param, l);
    ghost_sendrecv_double(op->prnX, _X, -1, &(op->c), minus_dir_param, l);
    // project plus dir and multiply with U dagger
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    dprn_su3_double(prp, phi, op, neighbor, 24 * start, 24 * n);
#else
    for (phi_pt = phi + 24 * start, end_pt = phi + 24 * n, D_pt = op->D + 36 * start,
        nb_pt = neighbor + 4 * start;
         phi_pt < end_pt; phi_pt += 24) {
      // T dir
      i = 12 * (*nb_pt);
      nb_pt++;
      dprn_T_double(pbuf, phi_pt);
      mvmh_double(op->prpT + i, D_pt, pbuf);
      mvmh_double(op->prpT + i + 3, D_pt, pbuf + 3);
      mvmh_double(op->prpT + i + 6, D_pt, pbuf + 6);
      mvmh_double(op->prpT + i + 9, D_pt, pbuf + 9);
      D_pt += 9;
      // Z dir
      i = 12 * (*nb_pt);
      nb_pt++;
      dprn_Z_double(pbuf, phi_pt);
      mvmh_double(op->prpZ + i, D_pt, pbuf);
      mvmh_double(op->prpZ + i + 3, D_pt, pbuf + 3);
      mvmh_double(op->prpZ + i + 6, D_pt, pbuf + 6);
      mvmh_double(op->prpZ + i + 9, D_pt, pbuf + 9);
      D_pt += 9;
      // Y dir
      i = 12 * (*nb_pt);
      nb_pt++;
      dprn_Y_double(pbuf, phi_pt);
      mvmh_double(op->prpY + i, D_pt, pbuf);
      mvmh_double(op->prpY + i + 3, D_pt, pbuf + 3);
      mvmh_double(op->prpY + i + 6, D_pt, pbuf + 6);
      mvmh_double(op->prpY + i + 9, D_pt, pbuf + 9);
      D_pt += 9;
      // X dir
      i = 12 * (*nb_pt);
      nb_pt++;
      dprn_X_double(pbuf, phi_pt);
      mvmh_double(op->prpX + i, D_pt, pbuf);
      mvmh_double(op->prpX + i + 3, D_pt, pbuf + 3);
      mvmh_double(op->prpX + i + 6, D_pt, pbuf + 6);
      mvmh_double(op->prpX + i + 9, D_pt, pbuf + 9);
      D_pt += 9;
    }
#endif
    if (amount == _EVEN_SITES) {
      start = start_even, n = end_even;
    } else if (amount == _ODD_SITES) {
      start = start_odd, n = end_odd;
    }
    // start communication in positive direction
    ghost_sendrecv_double(op->prpT, _T, +1, &(op->c), plus_dir_param, l);
    ghost_sendrecv_double(op->prpZ, _Z, +1, &(op->c), plus_dir_param, l);
    ghost_sendrecv_double(op->prpY, _Y, +1, &(op->c), plus_dir_param, l);
    ghost_sendrecv_double(op->prpX, _X, +1, &(op->c), plus_dir_param, l);
    // wait for communication in negative direction
    ghost_wait_double(op->prnT, _T, -1, &(op->c), minus_dir_param, l);
    ghost_wait_double(op->prnZ, _Z, -1, &(op->c), minus_dir_param, l);
    ghost_wait_double(op->prnY, _Y, -1, &(op->c), minus_dir_param, l);
    ghost_wait_double(op->prnX, _X, -1, &(op->c), minus_dir_param, l);
    // multiply with U and lift up minus dir
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    su3_dpbp_double(eta, prn, op, neighbor, 24 * start, 24 * n);
#else
    for (eta_pt = eta + 24 * start, end_pt = eta + 24 * n, D_pt = op->D + 36 * start,
        nb_pt = neighbor + 4 * start;
         eta_pt < end_pt; eta_pt += 24) {
      // T dir
      i = 12 * (*nb_pt);
      nb_pt++;
      mvm_double(pbuf, D_pt, op->prnT + i);
      mvm_double(pbuf + 3, D_pt, op->prnT + i + 3);
      mvm_double(pbuf + 6, D_pt, op->prnT + i + 6);
      mvm_double(pbuf + 9, D_pt, op->prnT + i + 9);
      dpbp_su3_T_double(pbuf, eta_pt);
      D_pt += 9;
      // Z dir
      i = 12 * (*nb_pt);
      nb_pt++;
      mvm_double(pbuf, D_pt, op->prnZ + i);
      mvm_double(pbuf + 3, D_pt, op->prnZ + i + 3);
      mvm_double(pbuf + 6, D_pt, op->prnZ + i + 6);
      mvm_double(pbuf + 9, D_pt, op->prnZ + i + 9);
      dpbp_su3_Z_double(pbuf, eta_pt);
      D_pt += 9;
      // Y dir
      i = 12 * (*nb_pt);
      nb_pt++;
      mvm_double(pbuf, D_pt, op->prnY + i);
      mvm_double(pbuf + 3, D_pt, op->prnY + i + 3);
      mvm_double(pbuf + 6, D_pt, op->prnY + i + 6);
      mvm_double(pbuf + 9, D_pt, op->prnY + i + 9);
      dpbp_su3_Y_double(pbuf, eta_pt);
      D_pt += 9;
      // X dir
      i = 12 * (*nb_pt);
      nb_pt++;
      mvm_double(pbuf, D_pt, op->prnX + i);
      mvm_double(pbuf + 3, D_pt, op->prnX + i + 3);
      mvm_double(pbuf + 6, D_pt, op->prnX + i + 6);
      mvm_double(pbuf + 9, D_pt, op->prnX + i + 9);
      dpbp_su3_X_double(pbuf, eta_pt);
      D_pt += 9;
    }
#endif
    // wait for communication in positive direction
    ghost_wait_double(op->prpT, _T, +1, &(op->c), plus_dir_param, l);
    ghost_wait_double(op->prpZ, _Z, +1, &(op->c), plus_dir_param, l);
    ghost_wait_double(op->prpY, _Y, +1, &(op->c), plus_dir_param, l);
    ghost_wait_double(op->prpX, _X, +1, &(op->c), plus_dir_param, l);
    // lift up plus dir
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    dpbn_double(eta, prp, 24 * start, 24 * n);
#else
    for (i = 12 * start, eta_pt = eta + 24 * start; i < 12 * n; i += 12, eta_pt += 24) {
      dpbn_su3_T_double(op->prpT + i, eta_pt);
      dpbn_su3_Z_double(op->prpZ + i, eta_pt);
      dpbn_su3_Y_double(op->prpY + i, eta_pt);
      dpbn_su3_X_double(op->prpX + i, eta_pt);
    }
#endif
  } else {
#endif
    // project in negative directions
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    prp_double(prn, phi, 12 * start, 12 * n);
#else
  complex_double pbuf[6];
  for (i = 6 * start, phi_pt = phi + 12 * start; i < 6 * n; i += 6, phi_pt += 12) {
    prp_T_double(prnT + i, phi_pt);
    prp_Z_double(prnZ + i, phi_pt);
    prp_Y_double(prnY + i, phi_pt);
    prp_X_double(prnX + i, phi_pt);
  }
#endif
    // start communication in negative direction
    ghost_sendrecv_double(prnT, _T, -1, &(op->c), minus_dir_param, l);
    ghost_sendrecv_double(prnZ, _Z, -1, &(op->c), minus_dir_param, l);
    ghost_sendrecv_double(prnY, _Y, -1, &(op->c), minus_dir_param, l);
    ghost_sendrecv_double(prnX, _X, -1, &(op->c), minus_dir_param, l);
    // project plus dir and multiply with U dagger
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    prn_su3_double(prp, phi, op, neighbor, 12 * start, 12 * n);
#else
  for (phi_pt = phi + 12 * start, end_pt = phi + 12 * n, D_pt = op->D + 36 * start,
      nb_pt = neighbor + 4 * start;
       phi_pt < end_pt; phi_pt += 12) {
    // T dir
    i = 6 * (*nb_pt);
    nb_pt++;
    prn_T_double(pbuf, phi_pt);
    mvmh_double(prpT + i, D_pt, pbuf);
    mvmh_double(prpT + i + 3, D_pt, pbuf + 3);
    D_pt += 9;
    // Z dir
    i = 6 * (*nb_pt);
    nb_pt++;
    prn_Z_double(pbuf, phi_pt);
    mvmh_double(prpZ + i, D_pt, pbuf);
    mvmh_double(prpZ + i + 3, D_pt, pbuf + 3);
    D_pt += 9;
    // Y dir
    i = 6 * (*nb_pt);
    nb_pt++;
    prn_Y_double(pbuf, phi_pt);
    mvmh_double(prpY + i, D_pt, pbuf);
    mvmh_double(prpY + i + 3, D_pt, pbuf + 3);
    D_pt += 9;
    // X dir
    i = 6 * (*nb_pt);
    nb_pt++;
    prn_X_double(pbuf, phi_pt);
    mvmh_double(prpX + i, D_pt, pbuf);
    mvmh_double(prpX + i + 3, D_pt, pbuf + 3);
    D_pt += 9;
  }
#endif
    if (amount == _EVEN_SITES) {
      start = start_even, n = end_even;
    } else if (amount == _ODD_SITES) {
      start = start_odd, n = end_odd;
    }
    // start communication in positive direction
    ghost_sendrecv_double(prpT, _T, +1, &(op->c), plus_dir_param, l);
    ghost_sendrecv_double(prpZ, _Z, +1, &(op->c), plus_dir_param, l);
    ghost_sendrecv_double(prpY, _Y, +1, &(op->c), plus_dir_param, l);
    ghost_sendrecv_double(prpX, _X, +1, &(op->c), plus_dir_param, l);
    // wait for communication in negative direction
    ghost_wait_double(prnT, _T, -1, &(op->c), minus_dir_param, l);
    ghost_wait_double(prnZ, _Z, -1, &(op->c), minus_dir_param, l);
    ghost_wait_double(prnY, _Y, -1, &(op->c), minus_dir_param, l);
    ghost_wait_double(prnX, _X, -1, &(op->c), minus_dir_param, l);
    // multiply with U and lift up minus dir
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    su3_pbp_double(eta, prn, op, neighbor, 12 * start, 12 * n);
#else
  for (eta_pt = eta + 12 * start, end_pt = eta + 12 * n, D_pt = op->D + 36 * start,
      nb_pt = neighbor + 4 * start;
       eta_pt < end_pt; eta_pt += 12) {
    // T dir
    i = 6 * (*nb_pt);
    nb_pt++;
    mvm_double(pbuf, D_pt, prnT + i);
    mvm_double(pbuf + 3, D_pt, prnT + i + 3);
    pbp_su3_T_double(pbuf, eta_pt);
    D_pt += 9;
    // Z dir
    i = 6 * (*nb_pt);
    nb_pt++;
    mvm_double(pbuf, D_pt, prnZ + i);
    mvm_double(pbuf + 3, D_pt, prnZ + i + 3);
    pbp_su3_Z_double(pbuf, eta_pt);
    D_pt += 9;
    // Y dir
    i = 6 * (*nb_pt);
    nb_pt++;
    mvm_double(pbuf, D_pt, prnY + i);
    mvm_double(pbuf + 3, D_pt, prnY + i + 3);
    pbp_su3_Y_double(pbuf, eta_pt);
    D_pt += 9;
    // X dir
    i = 6 * (*nb_pt);
    nb_pt++;
    mvm_double(pbuf, D_pt, prnX + i);
    mvm_double(pbuf + 3, D_pt, prnX + i + 3);
    pbp_su3_X_double(pbuf, eta_pt);
    D_pt += 9;
  }
#endif
    // wait for communication in positive direction
    ghost_wait_double(prpT, _T, +1, &(op->c), plus_dir_param, l);
    ghost_wait_double(prpZ, _Z, +1, &(op->c), plus_dir_param, l);
    ghost_wait_double(prpY, _Y, +1, &(op->c), plus_dir_param, l);
    ghost_wait_double(prpX, _X, +1, &(op->c), plus_dir_param, l);
    // lift up plus dir
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
    pbn_double(eta, prp, 12 * start, 12 * n);
#else
  for (i = 6 * start, eta_pt = eta + 12 * start; i < 6 * n; i += 6, eta_pt += 12) {
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
}

void apply_schur_complement_double(vector_double out, vector_double in, operator_struct<double> *op,
                                   Level *l) {

  /*********************************************************************************
   * Applies the Schur complement to a vector.
   *********************************************************************************/

  // start and end indices for vector functions
  int start_even = 0;
  int end_even = op->num_even_sites * l->num_lattice_site_var;
  int start_odd = op->num_even_sites * l->num_lattice_site_var;
  int end_odd = l->inner_vector_size;

  vector_double *tmp = op->buffer;

  vector_double_define(tmp[0], 0, start_odd, end_odd);
  vector_double_define(tmp[0], 0, start_even, end_even);
  PROF_double_START(_NC);
  PROF_double_START(_SC);
  diag_ee_double(out, in, op, l, start_even, end_even);
  PROF_double_STOP(_SC, 1);
  hopping_term_double(tmp[0], in, op, _ODD_SITES, l);
  PROF_double_STOP(_NC, 0);

  PROF_double_START(_SC);
  diag_oo_inv_double(tmp[1], tmp[0], op, l, start_odd, end_odd);
  PROF_double_STOP(_SC, 0);
  PROF_double_START(_NC);
  hopping_term_double(tmp[0], tmp[1], op, _EVEN_SITES, l);
  PROF_double_STOP(_NC, 1);
  vector_double_minus(out, out, tmp[0], start_even, end_even, l);
}

void solve_oddeven_double(Gmres *p, operator_struct<double> *op, Level *l) {

  int start = op->num_even_sites * l->num_lattice_site_var;
  int end = l->inner_vector_size;
  Global *global = l->global;
  vector_double tmp = op->buffer[0];

  // odd to even
  PROF_double_START(_SC);
  diag_oo_inv_double(tmp, p->b, op, l, start, end);

  PROF_double_STOP(_SC, 0);
  vector_double_scale(tmp, tmp, -1, start, end, l);
  PROF_double_START(_NC);
  hopping_term_double(p->b, tmp, op, _EVEN_SITES, l);
  PROF_double_STOP(_NC, 0);

  if (global->method == 4)
    p->solve();
  diag_oo_inv_double(p->x, p->b, op, l, start, end);

  // even to odd
  vector_double_define(tmp, 0, start, end);
  PROF_double_START(_NC);
  hopping_term_double(tmp, p->x, op, _ODD_SITES, l);
  PROF_double_STOP(_NC, 1);
  PROF_double_START(_SC);
  diag_oo_inv_double(p->b, tmp, op, l, start, end);
  PROF_double_STOP(_SC, 1);
  vector_double_minus(p->x, p->x, p->b, start, end, l);
}

void g5_double(vector_double eta, vector_double phi, int start, int end, Level *l) {
  if (eta != phi) {
    vector_double eta_end = eta + end;
    eta += start;
    phi += start;
    while (eta < eta_end) {
      FOR6(*eta = -(*phi); phi++; eta++;)
      FOR6(*eta = (*phi); phi++; eta++;)
    }
  } else {
    vector_double eta_end = eta + end;
    eta += start;
    phi += start;
    while (eta < eta_end) {
      FOR6(*eta = -(*phi); phi++; eta++;)
      eta += 6;
      phi += 6;
    }
  }
}

void minus_g5_double(vector_double eta, vector_double phi, int start, int end, Level *l) {
  if (eta != phi) {
    vector_double eta_end = eta + end;
    eta += start;
    phi += start;
    while (eta < eta_end) {
      FOR6(*eta = (*phi); phi++; eta++;)
      FOR6(*eta = -(*phi); phi++; eta++;)
    }
  } else {
    vector_double eta_end = eta + end;
    eta += start;
    phi += start;
    while (eta < eta_end) {
      eta += 6;
      phi += 6;
      FOR6(*eta = -(*phi); phi++; eta++;)
    }
  }
}

void g5D_apply_schur_complement_double(vector_double out, vector_double in,
                                       operator_struct<double> *op, Level *l) {

  /*********************************************************************************
   * Applies the Schur complement to a vector.
   *********************************************************************************/

  // start and end indices for vector functions
  int start_even = 0;
  int end_even = op->num_even_sites * l->num_lattice_site_var;
  int start_odd = op->num_even_sites * l->num_lattice_site_var;
  int end_odd = l->inner_vector_size;

  vector_double *tmp = op->buffer;

  vector_double_define(tmp[0], 0, start_odd, end_odd);
  vector_double_define(tmp[0], 0, start_even, end_even);
  PROF_double_START(_NC);

  PROF_double_START(_SC);
  diag_ee_double(out, in, op, l, start_even, end_even);
  PROF_double_STOP(_SC, 1);
  hopping_term_double(tmp[0], in, op, _ODD_SITES, l);
  PROF_double_STOP(_NC, 0);

  PROF_double_START(_SC);
  diag_oo_inv_double(tmp[1], tmp[0], op, l, start_odd, end_odd);
  PROF_double_STOP(_SC, 0);
  PROF_double_START(_NC);
  hopping_term_double(tmp[0], tmp[1], op, _EVEN_SITES, l);
  PROF_double_STOP(_NC, 1);
  vector_double_minus(out, out, tmp[0], start_even, end_even, l);
  g5_double(out, out, start_even, end_even, l);
  //   g5_double( out, out, start_odd, end_odd, l );
}

void g5D_solve_oddeven_double(Gmres *p, operator_struct<double> *op, Level *l) {

  int start_even = 0;
  int end_even = op->num_even_sites * l->num_lattice_site_var;
  int start_odd = op->num_even_sites * l->num_lattice_site_var;
  int end_odd = l->inner_vector_size;
  Global *global = l->global;
  vector_double tmp = op->buffer[0];

  // odd to even
  PROF_double_START(_SC);
  diag_oo_inv_double(tmp, p->b, op, l, start_odd, end_odd);
  PROF_double_STOP(_SC, 0);
  //   g5_double( tmp, tmp, start_odd, end_odd, l );
  //   vector_double_scale( tmp, tmp, -1, start_odd, end_odd, l );
  minus_g5_double(tmp, tmp, start_odd, end_odd, l);
  PROF_double_START(_NC);
  vector_double_define(p->x, 0, start_even, end_even);
  hopping_term_double(p->x, tmp, op, _EVEN_SITES, l);
  PROF_double_STOP(_NC, 0);
  g5_double(p->x, p->x, start_even, end_even, l);
  vector_double_plus(p->b, p->b, p->x, start_even, end_even, l);

  ASSERT(global->method == 6);
  p->solve();
  diag_oo_inv_double(p->x, p->b, op, l, start_odd, end_odd);
  g5_double(p->x, p->x, start_odd, end_odd, l);

  // even to odd
  vector_double_define(tmp, 0, start_odd, end_odd);
  PROF_double_START(_NC);
  hopping_term_double(tmp, p->x, op, _ODD_SITES, l);
  PROF_double_STOP(_NC, 1);
  PROF_double_START(_SC);
  diag_oo_inv_double(p->b, tmp, op, l, start_odd, end_odd);
  PROF_double_STOP(_SC, 1);
  vector_double_minus(p->x, p->x, p->b, start_odd, end_odd, l);
}

// ----- block odd even -----------------------------------------------------------

void schwarz_double_oddeven_setup(Schwarz *s, Level *l) {

  int mu, i, d0, c0, b0, a0, d1, c1, b1, a1, t, z, y, x, agg_split[4], block_split[4],
      block_size[4];
  operator_struct<double> *op = &(s->op);
  int n1 = s->number_of_even_block_sites;
  Global *global = l->global;
#ifdef HAVE_TM
  config_double tm_term_pt = op->tm_term;
#endif

  for (mu = 0; mu < 4; mu++) {
    agg_split[mu] = l->local_lattice[mu] / l->coarsening[mu];
    block_split[mu] = l->coarsening[mu] / l->block_lattice[mu];
    block_size[mu] = l->block_lattice[mu];
  }

  if (global->csw) {
#ifndef OPTIMIZED_SELF_COUPLING_double
    config_double clover_pt = op->clover, clover_oo_inv_pt = op->clover_oo_inv;
    complex_double buffer[42];
    int cs = 42;
#else
    double *clover_pt = op->clover_vectorized, *clover_oo_inv_pt = op->clover_oo_inv_vectorized;
    int cs = 144;
#endif
    for (d0 = 0; d0 < agg_split[_T]; d0++)
      for (c0 = 0; c0 < agg_split[_Z]; c0++)
        for (b0 = 0; b0 < agg_split[_Y]; b0++)
          for (a0 = 0; a0 < agg_split[_X]; a0++)

            for (d1 = d0 * block_split[_T]; d1 < (d0 + 1) * block_split[_T]; d1++)
              for (c1 = c0 * block_split[_Z]; c1 < (c0 + 1) * block_split[_Z]; c1++)
                for (b1 = b0 * block_split[_Y]; b1 < (b0 + 1) * block_split[_Y]; b1++)
                  for (a1 = a0 * block_split[_X]; a1 < (a0 + 1) * block_split[_X]; a1++) {

                    // skipping even sites
                    clover_pt += n1 * cs;
#ifdef HAVE_TM
                    tm_term_pt += n1 * 12;
#endif
                    for (t = d1 * block_size[_T]; t < (d1 + 1) * block_size[_T]; t++)
                      for (z = c1 * block_size[_Z]; z < (c1 + 1) * block_size[_Z]; z++)
                        for (y = b1 * block_size[_Y]; y < (b1 + 1) * block_size[_Y]; y++)
                          for (x = a1 * block_size[_X]; x < (a1 + 1) * block_size[_X]; x++) {
                            if (((t - d1 * block_size[_T]) + (z - c1 * block_size[_Z]) +
                                 (y - b1 * block_size[_Y]) + (x - a1 * block_size[_X])) %
                                    2 ==
                                1) {
#ifndef OPTIMIZED_SELF_COUPLING_double

                              for (i = 0; i < 42; i++)
                                buffer[i] = (complex_double)clover_pt[i];
#ifdef HAVE_TM
                              if (global->mu + global->mu_odd_shift != 0.0 ||
                                  global->mu + global->mu_even_shift != 0.0)
                                for (i = 0; i < 12; i++)
                                  buffer[i] += (complex_double)tm_term_pt[i];
                              tm_term_pt += 12;
                              selfcoupling_LU_decomposition_double(clover_oo_inv_pt, buffer);
                              clover_oo_inv_pt += 72;
#else
                              selfcoupling_cholesky_decomposition_double(clover_oo_inv_pt, buffer);
                              clover_oo_inv_pt += 42;
#endif

#else
                              sse_site_clover_invert_double(clover_pt, clover_oo_inv_pt);
                              clover_oo_inv_pt += 144;
#endif
                              clover_pt += cs;
                            }
                          }
                  }
  }

#ifdef HAVE_TM1p1
#ifndef OPTIMIZED_SELF_COUPLING_double
  complex_double buffer[66];
  config_double clover_oo_inv_pt = op->clover_doublet_oo_inv, clover_pt = op->clover;
  int cs = global->csw ? 42 : 12;
#else
  double *clover_pt = global->csw ? op->clover_doublet_vectorized : (double *)op->clover,
         *clover_oo_inv_pt = op->clover_doublet_oo_inv_vectorized;
  int cs = global->csw ? 288 : 24;
#endif
  config_double eps_term_pt = op->epsbar_term;
#ifdef HAVE_TM
  tm_term_pt = op->tm_term;
#endif

  for (d0 = 0; d0 < agg_split[_T]; d0++)
    for (c0 = 0; c0 < agg_split[_Z]; c0++)
      for (b0 = 0; b0 < agg_split[_Y]; b0++)
        for (a0 = 0; a0 < agg_split[_X]; a0++)

          for (d1 = d0 * block_split[_T]; d1 < (d0 + 1) * block_split[_T]; d1++)
            for (c1 = c0 * block_split[_Z]; c1 < (c0 + 1) * block_split[_Z]; c1++)
              for (b1 = b0 * block_split[_Y]; b1 < (b0 + 1) * block_split[_Y]; b1++)
                for (a1 = a0 * block_split[_X]; a1 < (a0 + 1) * block_split[_X]; a1++) {

                  // skipping even sites
                  clover_pt += n1 * cs;
                  eps_term_pt += n1 * 12;
#ifdef HAVE_TM
                  tm_term_pt += n1 * 12;
#endif
                  for (t = d1 * block_size[_T]; t < (d1 + 1) * block_size[_T]; t++)
                    for (z = c1 * block_size[_Z]; z < (c1 + 1) * block_size[_Z]; z++)
                      for (y = b1 * block_size[_Y]; y < (b1 + 1) * block_size[_Y]; y++)
                        for (x = a1 * block_size[_X]; x < (a1 + 1) * block_size[_X]; x++) {
                          if (((t - d1 * block_size[_T]) + (z - c1 * block_size[_Z]) +
                               (y - b1 * block_size[_Y]) + (x - a1 * block_size[_X])) %
                                  2 ==
                              1) {

#ifndef OPTIMIZED_SELF_COUPLING_double
                            if (global->csw) {
                              for (i = 0; i < 12; i++) // 0-23
                                buffer[i + 12] = buffer[i] = (complex_double)clover_pt[i];
                              for (i = 12; i < 42; i++) // 24-53
                                buffer[i + 12] = (complex_double)clover_pt[i];
                            } else {
                              for (i = 0; i < 12; i++) // 0-23
                                buffer[i + 12] = buffer[i] = (complex_double)clover_pt[i];
                              for (i = 12; i < 42; i++) // 24-53
                                buffer[i + 12] = 0;
                            }
                            for (i = 0; i < 12; i++) // 54-65
                              buffer[i + 54] = (complex_double)eps_term_pt[i];
#ifdef HAVE_TM
                            if (global->mu + global->mu_odd_shift != 0.0 ||
                                global->mu + global->mu_even_shift != 0.0)
                              for (int i = 0; i < 12; i++) { // 0-23
                                buffer[i] += (complex_double)tm_term_pt[i];
                                buffer[i + 12] -= (complex_double)tm_term_pt[i];
                              }
                            tm_term_pt += 12;
#endif
                            eps_term_pt += 12;
                            clover_pt += cs;
                            selfcoupling_LU_doublet_decomposition_double(clover_oo_inv_pt, buffer);
                            clover_oo_inv_pt += 288;
#else
                            if (global->csw) {
                              sse_site_clover_doublet_invert_double(clover_pt, eps_term_pt,
                                                                    clover_oo_inv_pt);
                            } else {
#ifdef HAVE_TM
                              for (i = 0; i < 6; i++) { // we temporaly save in clover_oo_inv_pt
                                clover_oo_inv_pt[2 * i] = clover_pt[2 * i] + real(tm_term_pt[i]);
                                clover_oo_inv_pt[2 * i + 1] =
                                    clover_pt[2 * i + 1] + imag(tm_term_pt[i]);
                                clover_oo_inv_pt[2 * i + 12] =
                                    clover_pt[2 * i] - real(tm_term_pt[i]);
                                clover_oo_inv_pt[2 * i + 13] =
                                    clover_pt[2 * i + 1] - imag(tm_term_pt[i]);
                              }
                              for (i = 6; i < 12; i++) {
                                clover_oo_inv_pt[2 * i + 12] =
                                    clover_pt[2 * i] + real(tm_term_pt[i]);
                                clover_oo_inv_pt[2 * i + 13] =
                                    clover_pt[2 * i + 1] + imag(tm_term_pt[i]);
                                clover_oo_inv_pt[2 * i + 24] =
                                    clover_pt[2 * i] - real(tm_term_pt[i]);
                                clover_oo_inv_pt[2 * i + 25] =
                                    clover_pt[2 * i + 1] - imag(tm_term_pt[i]);
                              }
                              tm_term_pt += 12;
#else
                              for (i = 0; i < 6; i++) {
                                clover_oo_inv_pt[2 * i + 12] = clover_oo_inv_pt[2 * i] =
                                    clover_pt[2 * i];
                                clover_oo_inv_pt[2 * i + 13] = clover_oo_inv_pt[2 * i + 1] =
                                    clover_pt[2 * i + 1];
                              }
                              for (i = 6; i < 12; i++) {
                                clover_oo_inv_pt[2 * i + 24] = clover_oo_inv_pt[2 * i + 12] =
                                    clover_pt[2 * i];
                                clover_oo_inv_pt[2 * i + 25] = clover_oo_inv_pt[2 * i + 13] =
                                    clover_pt[2 * i + 1];
                              }
#endif
                              sse_site_clover_doublet_invert_double(clover_oo_inv_pt, eps_term_pt,
                                                                    clover_oo_inv_pt);
                            }

                            clover_pt += cs;
                            eps_term_pt += 12;
                            clover_oo_inv_pt += 2 * 288;
#endif
                          }
                        }
                }
#endif
}

void block_diag_ee_double(vector_double eta, vector_double phi, int start, Schwarz *s, Level *l) {

  int n1 = s->number_of_even_block_sites, nv = l->num_lattice_site_var;
  clover_double(eta, phi, &(s->op), start, start + nv * n1, l);
}

// diagonal blocks applied to the odd sites of a block
void block_diag_oo_double(vector_double eta, vector_double phi, int start, Schwarz *s, Level *l) {

  Global *global = l->global;

#ifdef OPTIMIZED_SELF_COUPLING_double
  // we don't have the LU decomposition here, for debugging only
  int n1 = s->number_of_even_block_sites, n2 = s->number_of_odd_block_sites,
      nv = l->num_lattice_site_var;
  clover_double(eta, phi, &(s->op), start + nv * n1, start + nv * (n1 + n2), l);

#else

  int i, n1 = s->number_of_even_block_sites, n2 = s->number_of_odd_block_sites;
#ifdef HAVE_TM1p1
  if (global->n_flavours == 2) {
    int block_num = start / 24 / (n1 + n2);
    //    config_double clover = s->op.clover_doublet_oo_inv+n1*288+(start/24)*288;
    config_double clover = s->op.clover_doublet_oo_inv + (start / 24 - block_num * n1) * 288;
    vector_double lphi = phi + n1 * 24 + start, leta = eta + n1 * 24 + start;
    for (i = 0; i < n2; i++) {
      LU_multiply_double(leta, lphi, clover);
      leta += 24;
      lphi += 24;
      clover += 288;
    }
  } else {
#endif
    vector_double lphi = phi + n1 * 12 + start, leta = eta + n1 * 12 + start;
    if (global->csw) {
      int block_num = start / 12 / (n1 + n2);
#ifndef HAVE_TM
      config_double clover = s->op.clover_oo_inv + (start / 12 - block_num * n1) * 42;
      for (i = 0; i < n2; i++) {
        LLH_multiply_double(leta, lphi, clover);
        leta += 12;
        lphi += 12;
        clover += 42;
      }
#else
    config_double clover = s->op.clover_oo_inv + (start / 12 - block_num * n1) * 72;
    for (i = 0; i < n2; i++) {
      LU_multiply_double(leta, lphi, clover);
      leta += 12;
      lphi += 12;
      clover += 72;
    }
#endif
    } else {
      config_double clover = s->op.clover + n1 * 12 + start;
#ifndef HAVE_TM
      for (i = 0; i < 12 * n2; i++)
        leta[i] = lphi[i] * (clover[i]);
#else
    config_double tm_term = s->op.tm_term + n1 * 12 + start;
    for (i = 0; i < 12 * n2; i++)
      leta[i] = lphi[i] * (clover[i] + tm_term[i]);
#endif
    }
#ifdef HAVE_TM1p1
  }
#endif

#endif
}

// inverted diagonal blocks applied to the odd sites of a block
void block_diag_oo_inv_double(vector_double eta, vector_double phi, int start, Schwarz *s,
                              Level *l) {

  int i, n1 = s->number_of_even_block_sites, n2 = s->number_of_odd_block_sites;
  Global *global = l->global;
#ifdef HAVE_TM1p1
  if (global->n_flavours == 2) {

    vector_double lphi = phi + n1 * 24 + start, leta = eta + n1 * 24 + start;
    int block_num = start / 24 / (n1 + n2);
#ifndef OPTIMIZED_SELF_COUPLING_double
    config_double clover = s->op.clover_doublet_oo_inv + (start / 24 - block_num * n1) * 288;
    for (i = 0; i < n2; i++) {
      LU_perform_fwd_bwd_subs_double(leta, lphi, clover);
      leta += 24;
      lphi += 24;
      clover += 288;
    }
#else
    double *clover_vectorized =
        s->op.clover_doublet_oo_inv_vectorized + (start / 24 - block_num * n1) * 2 * 288;
    for (i = 0; i < n2; i++) {
      sse_site_clover_doublet_double((double *)leta, (double *)lphi, clover_vectorized);
      leta += 24;
      lphi += 24;
      clover_vectorized += 2 * 288;
    }
#endif

  } else {
#endif

    vector_double lphi = phi + n1 * 12 + start, leta = eta + n1 * 12 + start;
    if (global->csw) {
      int block_num = start / 12 / (n1 + n2);
#ifndef OPTIMIZED_SELF_COUPLING_double
#ifndef HAVE_TM
      config_double clover = s->op.clover_oo_inv + (start / 12 - block_num * n1) * 42;
      for (i = 0; i < n2; i++) {
        LLH_perform_fwd_bwd_subs_double(leta, lphi, clover);
        leta += 12;
        lphi += 12;
        clover += 42;
      }
#else
      config_double clover = s->op.clover_oo_inv + (start / 12 - block_num * n1) * 72;
      for (i = 0; i < n2; i++) {
        LU_perform_fwd_bwd_subs_double(leta, lphi, clover);
        leta += 12;
        lphi += 12;
        clover += 72;
      }
#endif
#else
    double *clover_vectorized =
        s->op.clover_oo_inv_vectorized + (start / 12 - block_num * n1) * 144;
    for (i = 0; i < n2; i++) {
      sse_site_clover_double((double *)leta, (double *)lphi, clover_vectorized);
      leta += 12;
      lphi += 12;
      clover_vectorized += 144;
    }
#endif
    } else {
      config_double clover = s->op.clover + n1 * 12 + start;
#ifndef HAVE_TM
      for (i = 0; i < 12 * n2; i++)
        leta[i] = lphi[i] / (clover[i]);
#else
    config_double tm_term = s->op.tm_term + n1 * 12 + start;
    for (i = 0; i < 12 * n2; i++)
      leta[i] = lphi[i] / (clover[i] + tm_term[i]);
#endif
    }
#ifdef HAVE_TM1p1
  }
#endif
}

void block_hopping_term_double(vector_double eta, vector_double phi, int start, int amount,
                               Schwarz *s, Level *l) {

  int a1, a2, n1, n2, *length_even = s->direction_length_even,
                      *length_odd = s->direction_length_odd, **index = s->odd_even_index_table,
                      *neighbor = s->op.neighbor_table, nv = l->num_lattice_site_var;

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
  double *Dplus = s->op.D_vectorized + (start / nv) * 96;
  double *Dminus = s->op.D_transformed_vectorized + (start / nv) * 96;

  for (int mu = 0; mu < 4; mu++) {
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[mu];
      a2 = n1;
      n2 = a2 + length_odd[mu];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[mu];
      n1 = a1 + length_odd[mu];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[mu] + length_odd[mu];
      a2 = 0;
      n2 = n1;
    }
    block_oddeven_plus_coupling_double((double *)(eta + start), Dplus, (double *)(phi + start), mu,
                                       a1, n1, index[mu], neighbor);
    block_oddeven_minus_coupling_double((double *)(eta + start), Dminus, (double *)(phi + start),
                                        mu, a2, n2, index[mu], neighbor);
  }

#else
  config_double D = s->op.D + (start / nv) * 36;
  int i, j, k, *ind;
  config_double D_pt;
  vector_double lphi = phi + start, leta = eta + start;

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2) {
    complex_double buf1[25] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                   *buf2 = buf1 + 12;
    // T direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_T];
      a2 = n1;
      n2 = a2 + length_odd[_T];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_T];
      n1 = a1 + length_odd[_T];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_T] + length_odd[_T];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +T coupling
    ind = index[_T];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _T];
      D_pt = D + 36 * k + 9 * _T;
      dprp_T_double(buf1, lphi + 24 * j);
      mvm_double(buf2, D_pt, buf1);
      mvm_double(buf2 + 3, D_pt, buf1 + 3);
      mvm_double(buf2 + 6, D_pt, buf1 + 6);
      mvm_double(buf2 + 9, D_pt, buf1 + 9);
      dpbp_su3_T_double(buf2, leta + 24 * k);
    }
    // "amount" of a block, -T coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _T];
      D_pt = D + 36 * k + 9 * _T;
      dprn_T_double(buf1, lphi + 24 * k);
      mvmh_double(buf2, D_pt, buf1);
      mvmh_double(buf2 + 3, D_pt, buf1 + 3);
      mvmh_double(buf2 + 6, D_pt, buf1 + 6);
      mvmh_double(buf2 + 9, D_pt, buf1 + 9);
      dpbn_su3_T_double(buf2, leta + 24 * j);
    }

    // Z direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_Z];
      a2 = n1;
      n2 = a2 + length_odd[_Z];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_Z];
      n1 = a1 + length_odd[_Z];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_Z] + length_odd[_Z];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +Z coupling
    ind = index[_Z];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Z];
      D_pt = D + 36 * k + 9 * _Z;
      dprp_Z_double(buf1, lphi + 24 * j);
      mvm_double(buf2, D_pt, buf1);
      mvm_double(buf2 + 3, D_pt, buf1 + 3);
      mvm_double(buf2 + 6, D_pt, buf1 + 6);
      mvm_double(buf2 + 9, D_pt, buf1 + 9);
      dpbp_su3_Z_double(buf2, leta + 24 * k);
    }
    // "amount" of a block, -Z coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Z];
      D_pt = D + 36 * k + 9 * _Z;
      dprn_Z_double(buf1, lphi + 24 * k);
      mvmh_double(buf2, D_pt, buf1);
      mvmh_double(buf2 + 3, D_pt, buf1 + 3);
      mvmh_double(buf2 + 6, D_pt, buf1 + 6);
      mvmh_double(buf2 + 9, D_pt, buf1 + 9);
      dpbn_su3_Z_double(buf2, leta + 24 * j);
    }

    // Y direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_Y];
      a2 = n1;
      n2 = a2 + length_odd[_Y];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_Y];
      n1 = a1 + length_odd[_Y];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_Y] + length_odd[_Y];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +Y coupling
    ind = index[_Y];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Y];
      D_pt = D + 36 * k + 9 * _Y;
      dprp_Y_double(buf1, lphi + 24 * j);
      mvm_double(buf2, D_pt, buf1);
      mvm_double(buf2 + 3, D_pt, buf1 + 3);
      mvm_double(buf2 + 6, D_pt, buf1 + 6);
      mvm_double(buf2 + 9, D_pt, buf1 + 9);
      dpbp_su3_Y_double(buf2, leta + 24 * k);
    }
    // "amount" of a block, -Y coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Y];
      D_pt = D + 36 * k + 9 * _Y;
      dprn_Y_double(buf1, lphi + 24 * k);
      mvmh_double(buf2, D_pt, buf1);
      mvmh_double(buf2 + 3, D_pt, buf1 + 3);
      mvmh_double(buf2 + 6, D_pt, buf1 + 6);
      mvmh_double(buf2 + 9, D_pt, buf1 + 9);
      dpbn_su3_Y_double(buf2, leta + 24 * j);
    }

    // X direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_X];
      a2 = n1;
      n2 = a2 + length_odd[_X];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_X];
      n1 = a1 + length_odd[_X];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_X] + length_odd[_X];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +X coupling
    ind = index[_X];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _X];
      D_pt = D + 36 * k + 9 * _X;
      dprp_X_double(buf1, lphi + 24 * j);
      mvm_double(buf2, D_pt, buf1);
      mvm_double(buf2 + 3, D_pt, buf1 + 3);
      mvm_double(buf2 + 6, D_pt, buf1 + 6);
      mvm_double(buf2 + 9, D_pt, buf1 + 9);
      dpbp_su3_X_double(buf2, leta + 24 * k);
    }
    // "amount" of a block, -X coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _X];
      D_pt = D + 36 * k + 9 * _X;
      dprn_X_double(buf1, lphi + 24 * k);
      mvmh_double(buf2, D_pt, buf1);
      mvmh_double(buf2 + 3, D_pt, buf1 + 3);
      mvmh_double(buf2 + 6, D_pt, buf1 + 6);
      mvmh_double(buf2 + 9, D_pt, buf1 + 9);
      dpbn_su3_X_double(buf2, leta + 24 * j);
    }
  } else {
#endif
    complex_double buf1[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, *buf2 = buf1 + 6;
    // T direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_T];
      a2 = n1;
      n2 = a2 + length_odd[_T];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_T];
      n1 = a1 + length_odd[_T];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_T] + length_odd[_T];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +T coupling
    ind = index[_T];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _T];
      D_pt = D + 36 * k + 9 * _T;
      prp_T_double(buf1, lphi + 12 * j);
      mvm_double(buf2, D_pt, buf1);
      mvm_double(buf2 + 3, D_pt, buf1 + 3);
      pbp_su3_T_double(buf2, leta + 12 * k);
    }
    // "amount" of a block, -T coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _T];
      D_pt = D + 36 * k + 9 * _T;
      prn_T_double(buf1, lphi + 12 * k);
      mvmh_double(buf2, D_pt, buf1);
      mvmh_double(buf2 + 3, D_pt, buf1 + 3);
      pbn_su3_T_double(buf2, leta + 12 * j);
    }

    // Z direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_Z];
      a2 = n1;
      n2 = a2 + length_odd[_Z];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_Z];
      n1 = a1 + length_odd[_Z];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_Z] + length_odd[_Z];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +Z coupling
    ind = index[_Z];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Z];
      D_pt = D + 36 * k + 9 * _Z;
      prp_Z_double(buf1, lphi + 12 * j);
      mvm_double(buf2, D_pt, buf1);
      mvm_double(buf2 + 3, D_pt, buf1 + 3);
      pbp_su3_Z_double(buf2, leta + 12 * k);
    }
    // "amount" of a block, -Z coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Z];
      D_pt = D + 36 * k + 9 * _Z;
      prn_Z_double(buf1, lphi + 12 * k);
      mvmh_double(buf2, D_pt, buf1);
      mvmh_double(buf2 + 3, D_pt, buf1 + 3);
      pbn_su3_Z_double(buf2, leta + 12 * j);
    }

    // Y direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_Y];
      a2 = n1;
      n2 = a2 + length_odd[_Y];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_Y];
      n1 = a1 + length_odd[_Y];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_Y] + length_odd[_Y];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +Y coupling
    ind = index[_Y];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Y];
      D_pt = D + 36 * k + 9 * _Y;
      prp_Y_double(buf1, lphi + 12 * j);
      mvm_double(buf2, D_pt, buf1);
      mvm_double(buf2 + 3, D_pt, buf1 + 3);
      pbp_su3_Y_double(buf2, leta + 12 * k);
    }
    // "amount" of a block, -Y coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Y];
      D_pt = D + 36 * k + 9 * _Y;
      prn_Y_double(buf1, lphi + 12 * k);
      mvmh_double(buf2, D_pt, buf1);
      mvmh_double(buf2 + 3, D_pt, buf1 + 3);
      pbn_su3_Y_double(buf2, leta + 12 * j);
    }

    // X direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_X];
      a2 = n1;
      n2 = a2 + length_odd[_X];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_X];
      n1 = a1 + length_odd[_X];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_X] + length_odd[_X];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +X coupling
    ind = index[_X];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _X];
      D_pt = D + 36 * k + 9 * _X;
      prp_X_double(buf1, lphi + 12 * j);
      mvm_double(buf2, D_pt, buf1);
      mvm_double(buf2 + 3, D_pt, buf1 + 3);
      pbp_su3_X_double(buf2, leta + 12 * k);
    }
    // "amount" of a block, -X coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _X];
      D_pt = D + 36 * k + 9 * _X;
      prn_X_double(buf1, lphi + 12 * k);
      mvmh_double(buf2, D_pt, buf1);
      mvmh_double(buf2 + 3, D_pt, buf1 + 3);
      pbn_su3_X_double(buf2, leta + 12 * j);
    }
#ifdef HAVE_TM1p1
  }
#endif
#endif
}

void block_n_hopping_term_double(vector_double eta, vector_double phi, int start, int amount,
                                 Schwarz *s, Level *l) {

  int a1, a2, n1, n2, *length_even = s->direction_length_even,
                      *length_odd = s->direction_length_odd, **index = s->odd_even_index_table,
                      *neighbor = s->op.neighbor_table, nv = l->num_lattice_site_var;

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
  double *Dplus = s->op.D_vectorized + (start / nv) * 96;
  double *Dminus = s->op.D_transformed_vectorized + (start / nv) * 96;

  for (int mu = 0; mu < 4; mu++) {
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[mu];
      a2 = n1;
      n2 = a2 + length_odd[mu];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[mu];
      n1 = a1 + length_odd[mu];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[mu] + length_odd[mu];
      a2 = 0;
      n2 = n1;
    }
    block_oddeven_nplus_coupling_double((double *)(eta + start), Dplus, (double *)(phi + start), mu,
                                        a1, n1, index[mu], neighbor);
    block_oddeven_nminus_coupling_double((double *)(eta + start), Dminus, (double *)(phi + start),
                                         mu, a2, n2, index[mu], neighbor);
  }

#else
  int i, j, k, *ind;
  vector_double lphi = phi + start, leta = eta + start;
  config_double D_pt, D = s->op.D + (start / nv) * 36;

#ifdef HAVE_TM1p1
  if (l->global->n_flavours == 2) {
    complex_double buf1[24], *buf2 = buf1 + 12;
    // T direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_T];
      a2 = n1;
      n2 = a2 + length_odd[_T];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_T];
      n1 = a1 + length_odd[_T];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_T] + length_odd[_T];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +T coupling
    ind = index[_T];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _T];
      D_pt = D + 36 * k + 9 * _T;
      dprp_T_double(buf1, lphi + 24 * j);
      nmvm_double(buf2, D_pt, buf1);
      nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      nmvm_double(buf2 + 6, D_pt, buf1 + 6);
      nmvm_double(buf2 + 9, D_pt, buf1 + 9);
      dpbp_su3_T_double(buf2, leta + 24 * k);
    }
    // "amount" of a block, -T coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _T];
      D_pt = D + 36 * k + 9 * _T;
      dprn_T_double(buf1, lphi + 24 * k);
      nmvmh_double(buf2, D_pt, buf1);
      nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      nmvmh_double(buf2 + 6, D_pt, buf1 + 6);
      nmvmh_double(buf2 + 9, D_pt, buf1 + 9);
      dpbn_su3_T_double(buf2, leta + 24 * j);
    }

    // Z direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_Z];
      a2 = n1;
      n2 = a2 + length_odd[_Z];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_Z];
      n1 = a1 + length_odd[_Z];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_Z] + length_odd[_Z];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +Z coupling
    ind = index[_Z];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Z];
      D_pt = D + 36 * k + 9 * _Z;
      dprp_Z_double(buf1, lphi + 24 * j);
      nmvm_double(buf2, D_pt, buf1);
      nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      nmvm_double(buf2 + 6, D_pt, buf1 + 6);
      nmvm_double(buf2 + 9, D_pt, buf1 + 9);
      dpbp_su3_Z_double(buf2, leta + 24 * k);
    }
    // "amount" of a block, -Z coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Z];
      D_pt = D + 36 * k + 9 * _Z;
      dprn_Z_double(buf1, lphi + 24 * k);
      nmvmh_double(buf2, D_pt, buf1);
      nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      nmvmh_double(buf2 + 6, D_pt, buf1 + 6);
      nmvmh_double(buf2 + 9, D_pt, buf1 + 9);
      dpbn_su3_Z_double(buf2, leta + 24 * j);
    }

    // Y direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_Y];
      a2 = n1;
      n2 = a2 + length_odd[_Y];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_Y];
      n1 = a1 + length_odd[_Y];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_Y] + length_odd[_Y];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +Y coupling
    ind = index[_Y];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Y];
      D_pt = D + 36 * k + 9 * _Y;
      dprp_Y_double(buf1, lphi + 24 * j);
      nmvm_double(buf2, D_pt, buf1);
      nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      nmvm_double(buf2 + 6, D_pt, buf1 + 6);
      nmvm_double(buf2 + 9, D_pt, buf1 + 9);
      dpbp_su3_Y_double(buf2, leta + 24 * k);
    }
    // "amount" of a block, -Y coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Y];
      D_pt = D + 36 * k + 9 * _Y;
      dprn_Y_double(buf1, lphi + 24 * k);
      nmvmh_double(buf2, D_pt, buf1);
      nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      nmvmh_double(buf2 + 6, D_pt, buf1 + 6);
      nmvmh_double(buf2 + 9, D_pt, buf1 + 9);
      dpbn_su3_Y_double(buf2, leta + 24 * j);
    }

    // X direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_X];
      a2 = n1;
      n2 = a2 + length_odd[_X];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_X];
      n1 = a1 + length_odd[_X];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_X] + length_odd[_X];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +X coupling
    ind = index[_X];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _X];
      D_pt = D + 36 * k + 9 * _X;
      dprp_X_double(buf1, lphi + 24 * j);
      nmvm_double(buf2, D_pt, buf1);
      nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      nmvm_double(buf2 + 6, D_pt, buf1 + 6);
      nmvm_double(buf2 + 9, D_pt, buf1 + 9);
      dpbp_su3_X_double(buf2, leta + 24 * k);
    }
    // "amount" of a block, -X coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _X];
      D_pt = D + 36 * k + 9 * _X;
      dprn_X_double(buf1, lphi + 24 * k);
      nmvmh_double(buf2, D_pt, buf1);
      nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      nmvmh_double(buf2 + 6, D_pt, buf1 + 6);
      nmvmh_double(buf2 + 9, D_pt, buf1 + 9);
      dpbn_su3_X_double(buf2, leta + 24 * j);
    }
  } else {
#endif
    complex_double buf1[12], *buf2 = buf1 + 6;

    // T direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_T];
      a2 = n1;
      n2 = a2 + length_odd[_T];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_T];
      n1 = a1 + length_odd[_T];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_T] + length_odd[_T];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +T coupling
    ind = index[_T];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _T];
      D_pt = D + 36 * k + 9 * _T;
      prp_T_double(buf1, lphi + 12 * j);
      nmvm_double(buf2, D_pt, buf1);
      nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      pbp_su3_T_double(buf2, leta + 12 * k);
    }
    // "amount" of a block, -T coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _T];
      D_pt = D + 36 * k + 9 * _T;
      prn_T_double(buf1, lphi + 12 * k);
      nmvmh_double(buf2, D_pt, buf1);
      nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      pbn_su3_T_double(buf2, leta + 12 * j);
    }

    // Z direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_Z];
      a2 = n1;
      n2 = a2 + length_odd[_Z];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_Z];
      n1 = a1 + length_odd[_Z];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_Z] + length_odd[_Z];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +Z coupling
    ind = index[_Z];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Z];
      D_pt = D + 36 * k + 9 * _Z;
      prp_Z_double(buf1, lphi + 12 * j);
      nmvm_double(buf2, D_pt, buf1);
      nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      pbp_su3_Z_double(buf2, leta + 12 * k);
    }
    // "amount" of a block, -Z coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Z];
      D_pt = D + 36 * k + 9 * _Z;
      prn_Z_double(buf1, lphi + 12 * k);
      nmvmh_double(buf2, D_pt, buf1);
      nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      pbn_su3_Z_double(buf2, leta + 12 * j);
    }

    // Y direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_Y];
      a2 = n1;
      n2 = a2 + length_odd[_Y];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_Y];
      n1 = a1 + length_odd[_Y];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_Y] + length_odd[_Y];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +Y coupling
    ind = index[_Y];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Y];
      D_pt = D + 36 * k + 9 * _Y;
      prp_Y_double(buf1, lphi + 12 * j);
      nmvm_double(buf2, D_pt, buf1);
      nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      pbp_su3_Y_double(buf2, leta + 12 * k);
    }
    // "amount" of a block, -Y coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _Y];
      D_pt = D + 36 * k + 9 * _Y;
      prn_Y_double(buf1, lphi + 12 * k);
      nmvmh_double(buf2, D_pt, buf1);
      nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      pbn_su3_Y_double(buf2, leta + 12 * j);
    }

    // X direction
    if (amount == _EVEN_SITES) {
      a1 = 0;
      n1 = length_even[_X];
      a2 = n1;
      n2 = a2 + length_odd[_X];
    } else if (amount == _ODD_SITES) {
      a1 = length_even[_X];
      n1 = a1 + length_odd[_X];
      a2 = 0;
      n2 = a1;
    } else {
      a1 = 0;
      n1 = length_even[_X] + length_odd[_X];
      a2 = 0;
      n2 = n1;
    }
    // "amount" of a block, +X coupling
    ind = index[_X];
    for (i = a1; i < n1; i++) {
      k = ind[i];
      j = neighbor[4 * k + _X];
      D_pt = D + 36 * k + 9 * _X;
      prp_X_double(buf1, lphi + 12 * j);
      nmvm_double(buf2, D_pt, buf1);
      nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      pbp_su3_X_double(buf2, leta + 12 * k);
    }
    // "amount" of a block, -X coupling
    for (i = a2; i < n2; i++) {
      k = ind[i];
      j = neighbor[4 * k + _X];
      D_pt = D + 36 * k + 9 * _X;
      prn_X_double(buf1, lphi + 12 * k);
      nmvmh_double(buf2, D_pt, buf1);
      nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      pbn_su3_X_double(buf2, leta + 12 * j);
    }
#ifdef HAVE_TM1p1
  }
#endif
#endif
}

void apply_block_schur_complement_double(vector_double out, vector_double in, int start, Schwarz *s,
                                         Level *l) {

  block_diag_ee_double(out, in, start, s, l);
  vector_double_define(s->odd_even_buffer[0], 0,
                       start + l->num_lattice_site_var * s->number_of_even_block_sites,
                       start + s->block_vector_size);
  block_hopping_term_double(s->odd_even_buffer[0], in, start, _ODD_SITES, s, l);
  block_diag_oo_inv_double(s->odd_even_buffer[1], s->odd_even_buffer[0], start, s, l);
  block_n_hopping_term_double(out, s->odd_even_buffer[1], start, _EVEN_SITES, s, l);
}
