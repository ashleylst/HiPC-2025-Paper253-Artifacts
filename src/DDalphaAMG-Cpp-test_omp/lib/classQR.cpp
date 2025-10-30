/*
#include "main.h"

void Orthogonalization::allocateMemory() {
  coefficients = new complex_double[2 * numberOfVectors];
  tmp = new complex_double[numberOfVectors];
  if (type == _DENSE_QR)
    work = new complex_double[lwork];
}

Orthogonalization::Orthogonalization()
    : Orthogonalization(0, 0, 0, 0, nullptr, nullptr) {} // delegating default constructor

Orthogonalization::Orthogonalization(int type, int vectorSize, int numberOfVectors,
                                     int beginOrthogonalization, Level *level, Thread *threading)
    : level(level)(threading), type(type), vectorSize(vectorSize),
      numberOfVectors(numberOfVectors), beginOrthogonalization(beginOrthogonalization),
lwork(vectorSize * numberOfVectors) { Orthogonalization::allocateMemory();
}

Orthogonalization::~Orthogonalization() {
  delete[] coefficients;
  delete[] tmp;
  if (type == _DENSE_QR)
    delete[] work;
}

void Orthogonalization::orthogonalize(complex_double **V) {
  switch (type) {
  case _DENSE_QR: denseQR(V); break;
  case _GRAM_SCHMID: gramSchmidt(V); break;
  case _MODIFIED_GRAM_SCHMIDT: modifiedGramSchmidt(V); break;
  case _HOUSEHOLDER_QR: householderQR(V); break;
  case _GRAM_SCHMID_ON_AGGREGATES: gramSchmidtOnAggregates(V); break;
  default: error0("Invalid Orthogonalization type\n");
  }
}

void Orthogonalization::denseQR(complex_double **V) {

  zgeqrf_(&vectorSize, &numberOfVectors, *V, &vectorSize, tmp, work, &lwork,
          &info); // compute QR factorization
  if (!info)
    zungqr_(&vectorSize, &numberOfVectors, &numberOfVectors, *V, &vectorSize, tmp, work, &lwork,
            &info); // recover Q from Householder reflections
  else
    error0("Error in zgeqrf_, Parameter %d had an illegal value\n", -1 * info);
  if (info)
    error0("Error in zungqr_, Parameter %d had an illegal value\n", -1 * info); // TODO: Try-catch?
}

void Orthogonalization::gramSchmidt(complex_double **V) {

  // NOTE: only thread safe, if "coefficients" is the same buffer for all threads belonging to a
common
  // MPI process

  compute_core_start_end_custom(0, vectorSize, &startThreading, &endThreading, level,
                                level->num_lattice_site_var);

  for (int i = beginOrthogonalization; i < numberOfVectors; i++) {
    process_multi_inner_product_double(i, tmp, V, V[i], 0, vectorSize, level);


    for (int j = 0; j < i; j++) {
      coefficients[j] = tmp[j];
    }


    if (i > 0) {

      MPI_Allreduce(coefficients, coefficients + numberOfVectors, i, MPI_COMPLEX_double, MPI_SUM,
                    (level->depth == 0) ? g.comm_cart : level->gs_double.level_comm);


    }
    for (int j = 0; j < i; j++) {
      vector_double_saxpy(V[i], V[i], V[j], -(coefficients + numberOfVectors)[j], startThreading,
                          endThreading, level);

    }

    norm = global_norm_double(V[i], 0, vectorSize, level);

    vector_double_scale(V[i], V[i], 1.0 / norm, startThreading, endThreading, level);

  }
}

void Orthogonalization::modifiedGramSchmidt(complex_double **V) {
  error0("modifiedGramSchmidt not yet implemented\n");
}

void Orthogonalization::householderQR(complex_double **V) {
  error0("householderQR not yet implemented\n");
}

void Orthogonalization::gramSchmidtOnAggregates(complex_double **V) {

  //   PROF_double_START(_GRAM_SCHMIDT_ON_AGGREGATES);


  int i, j, k, k1, k2, num_aggregates = level->s_double.num_aggregates,
                       aggregate_size = vectorSize / num_aggregates,
                       offset = level->num_lattice_site_var / 2;

  complex_double alpha1, alpha2;
  vector_double v_pt1, v_pt2;
  double norm1, norm2;

  for (j = threading->n_thread * threading->core + threading->thread; j < num_aggregates;
       j += threading->n_thread * threading->n_core) {
    for (k1 = 0; k1 < numberOfVectors; k1++) {
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
        for (k = 0; k < offset; k++, i++)
          v_pt1[i] = norm1 * real(v_pt1[i]) + I * norm1 * imag(v_pt1[i]);
        for (k = 0; k < offset; k++, i++)
          v_pt1[i] = norm2 * real(v_pt1[i]) + I * norm2 * imag(v_pt1[i]);
      }
    }
  }


  //   PROF_double_STOP(_GRAM_SCHMIDT_ON_AGGREGATES, 1);
}*/
