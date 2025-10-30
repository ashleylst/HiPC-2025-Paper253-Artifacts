/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel Simon Heybrock, Simone Bacchio, Bjoern
 * Leder.
 *
 * This file is part of the DDalphaAMG solver library.
 *
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The DDalphaAMG solver library is distributed in the hope that it will be useful
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 *
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 *
 */

/** @file linalg_generic.h
 * Linear Algebra functions .
 *
 * Contains the needed functions for linear algebra in DDalphaAMG such as inner products and norms.
 */

#ifndef LINALG_double_HEADER
#define LINALG_double_HEADER

#ifdef _M10TV
#define VECTOR_FOR(start, end, expression, update, l)                                              \
  do {                                                                                             \
    if (l->depth == 0) {                                                                           \
      for (start; end;)                                                                            \
        FOR12(expression; update;)                                                                 \
    } else {                                                                                       \
      for (start; end;)                                                                            \
        FOR20(expression; update;)                                                                 \
    }                                                                                              \
  } while (0)
#else
#define VECTOR_FOR(start, end, expression, update, l)                                              \
  do {                                                                                             \
    if (l->depth == 0) {                                                                           \
      for (start; end;)                                                                            \
        FOR12(expression; update;)                                                                 \
    } else {                                                                                       \
      for (start; end;)                                                                            \
        FOR2(expression; update;)                                                                  \
    }                                                                                              \
  } while (0)
#endif

#ifdef _M10TV
#define REAL_VECTOR_FOR(start, end, expression, update, l)                                         \
  do {                                                                                             \
    if (l->depth == 0) {                                                                           \
      for (start; end;)                                                                            \
        FOR24(expression; update;)                                                                 \
    } else {                                                                                       \
      for (start; end;)                                                                            \
        FOR40(expression; update;)                                                                 \
    }                                                                                              \
  } while (0)
#else
#define REAL_VECTOR_FOR(start, end, expression, update, l)                                         \
  do {                                                                                             \
    if (l->depth == 0) {                                                                           \
      for (start; end;)                                                                            \
        FOR24(expression; update;)                                                                 \
    } else {                                                                                       \
      for (start; end;)                                                                            \
        FOR4(expression; update;)                                                                  \
    }                                                                                              \
  } while (0)
#endif

#ifdef _M10TV
#define THREADED_VECTOR_FOR(i, start_index, end_index, expression, update, l, threading)           \
  do {                                                                                             \
    int thread_start, thread_end;                                                                  \
    if (l->depth == 0) {                                                                           \
      compute_core_start_end_custom(start_index, end_index, &thread_start, &thread_end,            \
                                    l threading, 12);                                              \
      for (i = thread_start; i < thread_end;)                                                      \
        FOR12(expression; update;)                                                                 \
    } else {                                                                                       \
      compute_core_start_end_custom(start_index, end_index, &thread_start, &thread_end,            \
                                    l threading, 20);                                              \
      for (i = thread_start; i < thread_end;)                                                      \
        FOR20(expression; update;)                                                                 \
    }                                                                                              \
  } while (0)
#else
#define THREADED_VECTOR_FOR(i, start_index, end_index, expression, update, l, threading)           \
  do {                                                                                             \
    int thread_start, thread_end;                                                                  \
    if (l->depth == 0) {                                                                           \
      compute_core_start_end_custom(start_index, end_index, &thread_start, &thread_end,            \
                                    l threading, 12);                                              \
      for (i = thread_start; i < thread_end;)                                                      \
        FOR12(expression; update;)                                                                 \
    } else {                                                                                       \
      compute_core_start_end_custom(start_index, end_index, &thread_start, &thread_end,            \
                                    l threading, 2);                                               \
      for (i = thread_start; i < thread_end;)                                                      \
        FOR2(expression; update;)                                                                  \
    }                                                                                              \
  } while (0)
#endif

/**
 * Computes the inner product of two globally distributed vectors.
 *
 * Computes the inner product of two globally distributed vectors.
 *
 * @param x the first vector in the product.
 * @param y the second vector in the product.
 * @param start the first position of the vectors.
 * @param end the last position of the vectors.
 * @param l the Level where the vectors are.
 * @param threading a Thread struct.
 * @return result of the inner product (complex valued).
 */
complex_double global_inner_product_double(vector_double x, vector_double y, int start, int end,
                                           Level *l);

complex_double process_inner_product_double(vector_double phi, vector_double psi, int start,
                                            int end, Level *l);

void process_multi_inner_product_double(int count, complex_double *results, vector_double *phi,
                                        vector_double psi, int start, int end, Level *l);

/**
 * Computes the norm of a globally distributed vector.
 *
 * Computes the norm of a globally distributed vector.
 *
 * @param phi the vector which norm is to be computed.
 * @param start the first position of the vector.
 * @param end the last position of the vector.
 * @param l the Level where the vector is.
 * @param threading a Thread struct.
 * @return result of the norm (double).
 */
double global_norm_double(vector_double phi, int start, int end, Level *l);

double process_norm_double(vector_double x, int start, int end, Level *l);

complex_double local_xy_over_xx_double(vector_double phi, vector_double psi, int start, int end,
                                       Level *l);

void vector_double_define(vector_double phi, complex_double value, int start, int end);

void vector_double_define_random(vector_double phi, int start, int end);

/**
 * Computes the sum of two globally distributed vectors.
 *
 *  z := x + y
 *
 * @param z the vector where the sum will be stored.
 * @param x the first vector in the sum.
 * @param y the second vector in the sum.
 * @param start the first position of the vectors.
 * @param end the last position of the vectors.
 * @param l the Level where the vectors are.
 * @param threading a Thread struct.
 */
void vector_double_plus(vector_double z, vector_double x, vector_double y, int start, int end,
                        Level *l); // z := x + y

/**
 * Computes the difference of two globally distributed vectors.
 *
 * z := x - y
 *
 * @param z the vector where the difference will be stored.
 * @param x the first vector in the difference.
 * @param y the second vector in the difference.
 * @param start the first position of the vectors.
 * @param end the last position of the vectors.
 * @param l the Level where the vectors are.
 */
void vector_double_minus(vector_double z, vector_double x, vector_double y, int start, int end,
                         Level *l); // z := x - y

/**
 * Computes the product of a vector with a complex scalar.
 *
 * z := alpha*x
 *
 * @param z the vector where the sum will be stored.
 * @param x the vector involved in the product.
 * @param alpha the scalar.
 * @param start the first position of the vectors.
 * @param end the last position of the vectors.
 * @param l the Level where the vectors are.
 * @param threading a Thread struct.
 */
void vector_double_scale(vector_double z, vector_double x, complex_double alpha, int start, int end,
                         Level *l); // z := alpha*x

/**
 * Computes the product of a vector with a real scalar.
 *
 * z := alpha*x
 *
 * @param z the vector where the sum will be stored.
 * @param x the vector involved in the product.
 * @param alpha the scalar.
 * @param start the first position of the vectors.
 * @param end the last position of the vectors.
 * @param l the Level where the vectors are.
 * @param threading a Thread struct.
 */
void vector_double_real_scale(vector_double z, vector_double x, complex_double alpha, int start,
                              int end, Level *l);

/**
 * Computes the saxpy operation in double.
 *
 *  z := x + alpha*y.
 *
 * @param z the vector where the result will be stored.
 * @param x the first vector in the operation.
 * @param y the second vector in the operation.
 * @param start the first position of the vectors.
 * @param alpha the involved scalar.
 * @param end the last position of the vectors.
 * @param l the Level where the vectors are.
 * @param threading a Thread struct.
 */
void vector_double_saxpy(vector_double z, vector_double x, vector_double y, complex_double alpha,
                         int start, int end, Level *l); // z := x + alpha*y

/**
 * Copies a vector into another.
 *
 * z := x
 *
 * @param z the vector where the copy will be stored.
 * @param x the vector to be copied.
 * @param start the first position of the vectors.
 * @param end the last position of the vectors.
 * @param l the Level where the vectors are.
 * @param threading a Thread struct.
 */
void vector_double_copy(vector_double z, vector_double x, int start, int end,
                        Level *l); // z := x

template<typename T>
void vector_copy(std::complex<T>* z, std::complex<T>* x, int start, int end){
  for(int i = start; i < end; i++){
    z[i-start] = x[i];
  }
}

template<typename T>
void vector_minus(std::complex<T>* z, std::complex<T> *x, std::complex<T> *y, int start, int end){
  for(int i = start; i < end; i++){
    z[i] = x[i] - y[i];
  }
}

void vector_double_copy_3(vector_double z, vector_double x, int start, int end, Level *l);

// TODO: Document
void vector_double_projection(vector_double z, complex_double *innerProducts, vector_double v,
                              int k, vector_double *W, int orthogonal, Level *l);

// TODO: Document
void gram_schmidt_on_aggregates_double(vector_double *V, const int num_vec, Level *l);

// TODO: Document
// Gram-Schmidt on a block of vectors, used by Block-Gram-Schmidt //TODO: Document
void aggregate_gram_schmidt_block_double(double *V, int num_vec, int leading_dimension, Level *l);
// used by Block-Gram-Schmidt //TODO: Document
void aggregate_orthogonalize_block_wrt_orthonormal_block_double(double *B, double *U, int num_vec,
                                                                Level *l);

// used by Block-Gram-Schmidt  //TODO: Document
void aggregate_block_dot_block_double(double *S, double *U, double *B, int num_vec,
                                      int leading_dimension, Level *l);

// used by Block-Gram-Schmidt  //TODO: Document
void aggregate_block_minus_block_times_dot_double(double *B, double *U, double *S, int num_vec,
                                                  int leading_dimension, Level *l);
// TODO: Document
void gram_schmidt_double(vector_double *V, complex_double *buffer, const int start, const int n,
                         Level *l);

// TODO: Document
void setup_gram_schmidt_double(vector_double *V, vector_double g5v, complex_double *buffer,
                               const int n, Level *l);

// TODO: Document
void spinwise_double_skalarmultiply(vector_double eta1, vector_double eta2, vector_double phi,
                                    complex_double alpha, int start, int end, Level *l);

/**
 * Sets all the components of a vector to a given value.
 *
 * Sets all the components of a vector to a given value.
 *
 * @param phi the vector which components are to be set.
 * @param alpha the complex-value to be used.
 * @param l the Level where the vector is.
 * @param threading a Thread struct.
 */
void set_boundary_double(vector_double phi, complex_double alpha, Level *l);

#endif
