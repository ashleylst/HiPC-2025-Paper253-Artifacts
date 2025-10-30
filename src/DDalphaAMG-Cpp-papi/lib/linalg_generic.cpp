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

#include "sse_double_intrinsic.h"
#include "sse_linalg.h"
#include "sse_linalg_generic.h"

#ifndef OPTIMIZED_LINALG_double
complex_double global_inner_product_double(vector_double phi, vector_double psi, int start, int end,
                                           Level *l) {

  PROF_double_START(_GIP);
  complex_double local_alpha = 0, global_alpha = 0;

  VECTOR_FOR(int i = start, i < end, local_alpha += conj(phi[i]) * psi[i], i++, l);

  if (l->global->num_processes > 1) {

    PROF_double_START(_ALLR);
    MPI_Allreduce(&local_alpha, &global_alpha, 1, MPI_COMPLEX_double, MPI_SUM,
                  (l->depth == 0) ? l->global->comm_cart : l->gs_double.level_comm);
    PROF_double_STOP(_ALLR, 1);
    PROF_double_STOP(_GIP, (double)(end - start) / (double)l->inner_vector_size);
    return global_alpha;
  } else {
    PROF_double_STOP(_GIP, (double)(end - start) / (double)l->inner_vector_size);
    return local_alpha;
  }

  return 0;
}
#endif

complex_double process_inner_product_double(vector_double phi, vector_double psi, int start,
                                            int end, Level *l) {

  PROF_double_START(_PIP);
  complex_double local_alpha = 0;

  VECTOR_FOR(int i = start, i < end, local_alpha += conj(phi[i]) * psi[i], i++, l);

  PROF_double_STOP(_PIP, (double)(end - start) / (double)l->inner_vector_size);

  return local_alpha;
}

#if !defined(OPTIMIZED_LINALG_double)
void process_multi_inner_product_double(int count, complex_double *results, vector_double *phi,
                                        vector_double psi, int start, int end, Level *l) {

  PROF_double_START(_PIP);
  int i;
  for (int c = 0; c < count; c++)
    results[c] = 0.0;

  if (l->depth == 0) {
    for (int c = 0; c < count; c++)
      for (i = start; i < end;)
        FOR12(results[c] += conj(phi[c][i]) * psi[i]; i++;)
  } else {
#ifdef _M10TV
    for (int c = 0; c < count; c++)
      for (i = start; i < end;)
        FOR20(results[c] += conj(phi[c][i]) * psi[i]; i++;)
#else
    for (int c = 0; c < count; c++)
      for (i = start; i < end;)
        FOR2(results[c] += conj(phi[c][i]) * psi[i]; i++;)
#endif
  }

  PROF_double_STOP(_PIP, (double)(end - start) / (double)l->inner_vector_size);
}
#endif

complex_double local_xy_over_xx_double(vector_double phi, vector_double psi, int start, int end,
                                       Level *l) {

  complex_double numerator = 0.0;
  double denominator = 0.0;

  VECTOR_FOR(int i = start, i < end, numerator += conj(phi[i]) * psi[i];
             denominator += NORM_SQUARE_double(phi[i]), i++, l);

  if (abs(denominator) < EPS_double) {
    return 0.0;
  }

  return numerator / denominator;
}

#ifndef OPTIMIZED_LINALG_double
double global_norm_double(vector_double x, int start, int end, Level *l) {

  PROF_double_START(_GIP);

  double local_alpha = 0, global_alpha = 0;

  for(int i = start; i < end; i++){
    local_alpha += NORM_SQUARE_double(x[i]);
  }

  if (l->global->num_processes > 1) {

    PROF_double_START(_ALLR);
    MPI_Allreduce(&local_alpha, &global_alpha, 1, MPI_double, MPI_SUM,
                  (l->depth == 0) ? l->global->comm_cart : l->gs_double.level_comm);
    PROF_double_STOP(_ALLR, 1);
    PROF_double_STOP(_GIP, (double)(end - start) / (double)l->inner_vector_size);
    return (double)sqrt((double)global_alpha);
  } else {
    PROF_double_STOP(_GIP, (double)(end - start) / (double)l->inner_vector_size);
    return (double)sqrt((double)local_alpha);
  }

  return 0;
}
#endif

double process_norm_double(vector_double x, int start, int end, Level *l) {

  double local_alpha = 0;
  PROF_double_START(_PIP);

  VECTOR_FOR(int i = start, i < end, local_alpha += NORM_SQUARE_double(x[i]), i++, l);

  PROF_double_STOP(_PIP, (double)(end - start) / (double)l->inner_vector_size);

  return (double)sqrt((double)local_alpha);
}

void vector_double_define(vector_double phi, complex_double value, int start, int end) {

  ASSERT(phi != nullptr);
  for (int i = start; i < end; i++)
    phi[i] = value;
}

void vector_double_define_random(vector_double phi, int start, int end) {

  ASSERT(phi != nullptr);
  ASSERT(start != end);

  std::random_device seed;
  std::mt19937 randomInt(seed());
  std::uniform_real_distribution<> unit(1.0, 10.0);

  PROF_double_START(_SET);
  for (int i = start; i < end; i++)
    phi[i] = unit(randomInt) + unit(randomInt) * I;

  PROF_double_STOP(_SET, 1);
}

void vector_double_plus(vector_double z, vector_double x, vector_double y, int start, int end,
                        Level *l) {

  if (start != end)
    PROF_double_START(_LA2);

  VECTOR_FOR(int i = start, i < end, z[i] = x[i] + y[i], i++, l);

  if (start != end)
    PROF_double_STOP(_LA2, (double)(end - start) / (double)l->inner_vector_size);
}

void vector_double_minus(vector_double z, vector_double x, vector_double y, int start, int end,
                         Level *l) {

  if (start != end)
    PROF_double_START(_LA2);

  for(int i = start; i < end; i++){
    z[i] = x[i] - y[i];
  }

  if (start != end)
    PROF_double_STOP(_LA2, (double)(end - start) / (double)l->inner_vector_size);
}

#ifndef OPTIMIZED_LINALG_double
void vector_double_scale(vector_double z, vector_double x, complex_double alpha, int start, int end,
                         Level *l) {

  if (start != end)
    PROF_double_START(_LA6);

  VECTOR_FOR(int i = start, i < end, z[i] = alpha * x[i], i++, l);

  if (start != end)
    PROF_double_STOP(_LA6, (double)(end - start) / (double)l->inner_vector_size);
}
#endif

void vector_double_real_scale(vector_double z, vector_double x, complex_double alpha, int start,
                              int end, Level *l) {

  double *r_z = (double *)z, *r_x = (double *)x, r_alpha = real(alpha);
  int r_start = 2 * start, r_end = 2 * end;

  if (start != end)
    PROF_double_START(_LA2);

  REAL_VECTOR_FOR(int i = r_start, i < r_end, r_z[i] = r_alpha * r_x[i], i++, l);

  if (start != end)
    PROF_double_STOP(_LA2, (double)(end - start) / (double)l->inner_vector_size);
}

void vector_double_copy(vector_double z, vector_double x, int start, int end, Level *l) {

  if (start != end)
    PROF_double_START(_CPY);

  VECTOR_FOR(int i = start, i < end, z[i] = x[i], i++, l);

  if (start != end)
    PROF_double_STOP(_CPY, (double)(end - start) / (double)l->inner_vector_size);
}

void vector_double_copy_3(vector_double z, vector_double x, int start, int end, Level *l){
  for(int i = start; i < end; i++){
    z[i] = x[i-start];
  }
}

#ifndef OPTIMIZED_LINALG_double
void vector_double_saxpy(vector_double z, vector_double x, vector_double y, complex_double alpha,
                         int start, int end, Level *l) {

  if (start != end)
    PROF_double_START(_LA8);

  VECTOR_FOR(int i = start, i < end, z[i] = x[i] + alpha * y[i], i++, l);

  if (start != end)
    PROF_double_STOP(_LA8, (double)(end - start) / (double)l->inner_vector_size);
}
#endif

#ifndef OPTIMIZED_LINALG_double
void vector_double_multi_saxpy(vector_double z, vector_double *V, complex_double *alpha, int sign,
                               int count, int start, int end, Level *l) {

  if (start != end)
    PROF_double_START(_LA8);

  complex_double *alpha_signed = nullptr;
  MALLOC(alpha_signed, complex_double, count);

  for (int c = 0; c < count; c++) {
    alpha_signed[c] = (double)sign * alpha[c];
  }

  for (int c = 0; c < count; c++) {
    for (int i = start; i < end;) {
      FOR12(z[i] += V[c][i] * alpha_signed[c]; i++;)
    }
  }

  if (start != end)
    PROF_double_STOP(_LA8, (double)(count));
  FREE(alpha_signed, complex_double, count);
}
#endif

void vector_double_projection(vector_double out, complex_double *innerProducts, vector_double in,
                              int k, vector_double *V, int orthogonal, Level *l) {
  int start = 0, end = l->inner_vector_size;

  vector_double v_tmp = nullptr;
  complex_double *ip_buffer = nullptr;

  MALLOC(ip_buffer, complex_double, 2 * k);
  MALLOC(v_tmp, complex_double, l->inner_vector_size);
  vector_double_define(v_tmp, 0, start, end);

  process_multi_inner_product_double(k, innerProducts, V, in, start, end, l);

  for (int i = 0; i < k; i++) {
    ip_buffer[i] = innerProducts[i];
  }
  MPI_Allreduce(ip_buffer, ip_buffer + k, k, MPI_COMPLEX_double, MPI_SUM,
                (l->depth == 0) ? l->global->comm_cart : l->gs_double.level_comm);

  for (int i = 0; i < k; i++) {
    innerProducts[i] = (ip_buffer + k)[i];
    printf0("innerProducts[%d]=%lf%+lfi\n", i, CSPLIT(innerProducts[i]));
  }

  vector_double_multi_saxpy(v_tmp, V, ip_buffer + k, 1, k, start, end, l);

  if (orthogonal)
    vector_double_minus(out, in, v_tmp, start, end, l);
  else
    vector_double_copy(out, v_tmp, start, end, l);

  FREE(ip_buffer, complex_double, 2 * k);
  FREE(v_tmp, complex_double, l->inner_vector_size);
}

void gram_schmidt_on_aggregates_double(vector_double *V, const int num_vect, Level *l) {

  PROF_double_START(_GRAM_SCHMIDT_ON_AGGREGATES);

  int i, j, k, k1, k2, num_aggregates = l->s_double.number_of_aggregates,
                       aggregate_size = l->inner_vector_size / num_aggregates,
                       offset = l->num_lattice_site_var / 2;

  complex_double alpha1, alpha2;
  vector_double v_pt1, v_pt2;
  double norm1, norm2;

  for (j = 0; j < num_aggregates; j++) {
    for (k1 = 0; k1 < num_vect; k1++) {
      v_pt1 = V[k1] + j * aggregate_size;
      for (k2 = 0; k2 < k1; k2++) {
        v_pt2 = V[k2] + j * aggregate_size;
        alpha1 = 0;
        alpha2 = 0;
        // V[k1] -= <V[k2],V[k1]> V[k2] | 2*j-th and 2*j+1-st aggregate
        for (i = 0; i < aggregate_size;) {
          for (k = 0; k < offset; k++, i++)
            alpha1 += conj(v_pt2[i]) * v_pt1[i];
          for (k = 0; k < offset; k++, i++)
            alpha2 += conj(v_pt2[i]) * v_pt1[i];
        }
        for (i = 0; i < aggregate_size;) {
          for (k = 0; k < offset; k++, i++)
            v_pt1[i] -= alpha1 * v_pt2[i];
          for (k = 0; k < offset; k++, i++)
            v_pt1[i] -= alpha2 * v_pt2[i];
        }
      }

      norm1 = 0;
      norm2 = 0;
      // V[k1] = V[k1]/norm(V[k1]) | 2*j-th and 2*j+1-st aggregate
      for (i = 0; i < aggregate_size;) {
        for (k = 0; k < offset; k++, i++)
          norm1 += NORM_SQUARE_double(v_pt1[i]);
        for (k = 0; k < offset; k++, i++)
          norm2 += NORM_SQUARE_double(v_pt1[i]);
      }

      norm1 = 1 / sqrt(norm1);
      norm2 = 1 / sqrt(norm2);

      for (i = 0; i < aggregate_size;) {
        for (k = 0; k < offset; k++, i++) {
          v_pt1[i] = norm1 * real(v_pt1[i]) + I * norm1 * imag(v_pt1[i]);
        }
        for (k = 0; k < offset; k++, i++) {
          v_pt1[i] = norm2 * real(v_pt1[i]) + I * norm2 * imag(v_pt1[i]);
        }
      }
    }
  }

  PROF_double_STOP(_GRAM_SCHMIDT_ON_AGGREGATES, 1);
}

void spinwise_double_skalarmultiply(vector_double eta1, vector_double eta2, vector_double phi,
                                    complex_double alpha, int start, int end, Level *l) {

  PROF_double_START(_LA6);
  for (int i = start; i < end;) {
    FOR6(eta1[i] = alpha * phi[i]; eta2[i] = std::complex<double>(0, 0); i++;)
    FOR6(eta2[i] = alpha * phi[i]; eta1[i] = std::complex<double>(0, 0); i++;)
  }
  PROF_double_STOP(_LA6, 1);
}

void set_boundary_double(vector_double phi, complex_double alpha, Level *l) {

  PROF_double_START(_SET);

//   for( int i=l->inner_vector_size; i<l->vector_size; i++)
//     phi[i] = alpha;
  VECTOR_FOR(int i=l->inner_vector_size, i < l->vector_size, phi[i] = alpha, i++, l);

  PROF_double_STOP(_SET,
                   (double)(l->vector_size - l->inner_vector_size) / (double)l->inner_vector_size);
}

void gram_schmidt_double(vector_double *V, complex_double *buffer, const int begin, const int n,
                         Level *l) {

  PROF_double_START(_LA);
  double beta;
  int i, j, start = 0, end = l->inner_vector_size;
  complex_double *tmp = new complex<double>[n];

  for (i = begin; i < n; i++) {
    process_multi_inner_product_double(i, tmp, V, V[i], 0, l->inner_vector_size, l);
    for (j = 0; j < i; j++) {
      buffer[j] = tmp[j];
    }
    if (i > 0) {
      PROF_double_START(_ALLR);
      MPI_Allreduce(buffer, buffer + n, i, MPI_COMPLEX_double, MPI_SUM,
                    (l->type == _FINE) ? l->global->comm_cart : l->gs_double.level_comm);
      PROF_double_STOP(_ALLR, 1);
    }
    for (j = 0; j < i; j++) {
      vector_double_saxpy(V[i], V[i], V[j], -(buffer + n)[j], start, end, l);
    }
    beta = global_norm_double(V[i], 0, l->inner_vector_size, l);
    vector_double_real_scale(V[i], V[i], real(1.0 / beta), start, end, l);
  }
  delete[] tmp;
  PROF_double_STOP(_LA, 1);
}

#if !defined(SSE) || !defined(GRAM_SCHMIDT_VECTORIZED_double)
void setup_gram_schmidt_double_compute_dots(complex_double *buffer, vector_double *V, int count,
                                            int offset, int start, int end, Level *l) {

  int cache_block_size = 12 * 64;
  complex_double *tmp = new complex_double[cache_block_size];
  for (int i = 0; i < 2 * offset; i++)
    buffer[i] = 0.0;

  for (int i = start; i < end; i += cache_block_size) {
    coarse_gamma5_double(tmp, V[count] + i, 0, cache_block_size, l);
    for (int j = 0; j < count; j++) {
      for (int k = 0; k < cache_block_size; k++) {
        buffer[j] += conj(V[j][i + k]) * V[count][i + k];
        buffer[j + offset] += conj(V[j][i + k]) * tmp[k];
      }
    }
  }

  delete[] tmp;
}
#endif

#if !defined(SSE) || !defined(GRAM_SCHMIDT_VECTORIZED_double)
void setup_gram_schmidt_double_axpys(complex_double *buffer, vector_double *V, int count,
                                     int offset, int start, int end, Level *l) {

  int cache_block_size = 12 * 64;
  complex_double *tmp = new complex_double[cache_block_size];

  for (int i = start; i < end; i += cache_block_size) {
    for (int j = 0; j < count; j++) {
      coarse_gamma5_double(tmp, V[j] + i, 0, cache_block_size, l);
      for (int k = 0; k < cache_block_size; k++) {
        V[count][i + k] -= buffer[2 * offset + j] * V[j][i + k];
        V[count][i + k] -= buffer[3 * offset + j] * tmp[k];
      }
    }
  }
  delete[] tmp;
}
#endif

void setup_gram_schmidt_double(vector_double *V, vector_double g5v, complex_double *buffer,
                               const int n, Level *l) {

  PROF_double_START(_GRAM_SCHMIDT);
  double beta;
  int i, j;
  int start = 0;
  int end = l->inner_vector_size;

  complex_double *buffer2 = new complex_double[4 * n];
  for (i = 0; i < 4 * n; i++)
    buffer2[i] = 0;

  for (i = 0; i < n; i++) {

    if (l->depth > 0) {
      coarse_gamma5_double(g5v, V[i], start, end, l);
      for (j = 0; j < i; j++) {
        buffer2[j] = process_inner_product_double(V[j], V[i], start, end, l);
        buffer2[j + n] = process_inner_product_double(V[j], g5v, start, end, l);
      }
    } else
      setup_gram_schmidt_double_compute_dots(buffer2, V, i, n, start, end, l);

    if (i > 0) {
      PROF_double_START(_ALLR);
      MPI_Allreduce(buffer2, buffer2 + 2 * n, 2 * n, MPI_COMPLEX_double, MPI_SUM,
                    (l->type == _FINE) ? l->global->comm_cart : l->gs_double.level_comm);
      PROF_double_STOP(_ALLR, 1);
    }

    if (l->depth > 0) {
      for (j = 0; j < i; j++) {
        vector_double_saxpy(V[i], V[i], V[j], -(buffer2 + 2 * n)[j], start, end, l);
        coarse_gamma5_double(g5v, V[j], start, end, l);
        vector_double_saxpy(V[i], V[i], g5v, -(buffer2 + 3 * n)[j], start, end, l);
      }
    } else {
      setup_gram_schmidt_double_axpys(buffer2, V, i, n, start, end, l);
    }

    beta = global_norm_double(V[i], start, end, l);
    vector_double_real_scale(V[i], V[i], 1.0 / beta, start, end, l);
  }
  delete[] buffer2;
  PROF_double_STOP(_GRAM_SCHMIDT, (double)(end - start) / (double)l->inner_vector_size);
}
