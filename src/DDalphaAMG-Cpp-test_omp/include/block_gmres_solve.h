//
// Created by shiting on 2024-11-18.
//

#ifndef DDALPHAAMG_BLOCK_GMRES_SOLVE_H
#define DDALPHAAMG_BLOCK_GMRES_SOLVE_H

#include "classBlockSolver.h"

static inline int get_H_idx(int row, int col, int bnum, int blen){
  return (row*row + 3*row)/2*blen + blen*col + bnum;
}

template<typename T>
static inline T max_in_block(T *z, int blen){
  T max = z[0];
  for (int i = 1; i < blen; i++) {
    if(z[i] > max){
      max = z[i];
    }
  }
  return max;
}

template<typename T>
void check_staganation(std::complex<T> *z, bool *marked, T tol, int blen){
  for (int i = 0; i < blen; ++i) {
    if(abs(z[i]) < tol){
      marked[i] = true;
    }
  }
}

/**
 * y = ||z||_2
 * @tparam T precision
 */
template<typename T>
void set_global_norm(std::complex<T> *z, std::complex<T> *y, int end, int blen, Level *l){
  auto *local_alpha = new T[blen];
  auto *global_alpha = new T[blen];
  memset(local_alpha, 0, blen * sizeof (T));
  memset(global_alpha, 0, blen * sizeof (T));

  int nv = l->num_lattice_site_var;
  int num_lattice_sites = end/nv;

  //std::cout << local_alpha[2] << std::endl;
  for (int i = 0; i < num_lattice_sites; i++) {
    for (int j = 0; j < blen; j++) {
      for (int k = 0; k < nv; k++) {
        local_alpha[j] += real(z[i*blen*nv + j*nv + k]) * real(z[i*blen*nv + j*nv + k]) +
                          imag(z[i*blen*nv + j*nv + k]) * imag(z[i*blen*nv + j*nv + k]);
      }
    }
  }

  //std::cout << std::endl;
  if (l->global->num_processes > 1) {
    //std::cout << local_alpha[2] << std::endl;
    MPI_Allreduce(local_alpha, global_alpha, blen, MPI_double, MPI_SUM,
                  (l->depth == 0) ? l->global->comm_cart : l->gs_double.level_comm);
    for (int i = 0; i < blen; i++) {
      y[i] = (T)sqrt((T)global_alpha[i]);
    }
  } else {
    for (int i = 0; i < blen; i++) {
      y[i] = (T)sqrt((T)local_alpha[i]);
    }
  }

  delete[] local_alpha;
  delete[] global_alpha;
}

template<typename T>
void set_global_norm_per_entry(std::complex<T> *z, std::complex<T> *y, int end, int blen, Level *l){
  auto *local_alpha = new T[blen];
  auto *global_alpha = new T[blen];
  memset(local_alpha, 0, blen * sizeof (T));
  memset(global_alpha, 0, blen * sizeof (T));


  //std::cout << local_alpha[2] << std::endl;
  for (int i = 0; i < end; i++) {
    for (int j = 0; j < blen; j++) {
      local_alpha[j] += real(z[i*blen + j]) * real(z[i*blen + j]) +
                        imag(z[i*blen + j]) * imag(z[i*blen + j]);
    }
  }
  if (l->global->num_processes > 1) {
    //std::cout << local_alpha[2] << std::endl;
    MPI_Allreduce(local_alpha, global_alpha, blen, MPI_double, MPI_SUM,
                  (l->depth == 0) ? l->global->comm_cart : l->gs_double.level_comm);
    for (int i = 0; i < blen; i++) {
      y[i] = (T)sqrt((T)global_alpha[i]);
    }
  } else {
    for (int i = 0; i < blen; i++) {
      y[i] = (T)sqrt((T)local_alpha[i]);
    }
  }

  delete[] local_alpha;
  delete[] global_alpha;
}

/** z = x-alpha*y
 * */
template<typename T>
void vector_saxpy_neg(std::complex<T> *z, std::complex<T> *x, std::complex<T> *y,
                      std::complex<T> *alpha, int end, int blen, Level *l){
  int nv = l->num_lattice_site_var;
  int num_lattice_sites = end/nv;

  for (int i = 0; i < num_lattice_sites; i++) {
    for (int j = 0; j < blen; j++) {
      for (int k = 0; k < nv; k++) {
        z[i*blen*nv + j*nv + k] = x[i*blen*nv + j*nv + k] - alpha[j] * y[i*blen*nv + j*nv + k];
      }
    }
  }
}

template<typename T>
void vector_saxpy_neg_per_entry(std::complex<T> *z, std::complex<T> *x, std::complex<T> *y,
                                std::complex<T> *alpha, int end, int blen, Level *l){
  for (int i = 0; i < end; i++) {
    for (int j = 0; j < blen; j++) {
      z[i*blen + j] = x[i*blen + j] - alpha[j] * y[i*blen + j];
    }
  }
}

/** z = x+alpha*y
 * */
template<typename T>
void vector_saxpy_with_mark(std::complex<T> *z, std::complex<T> *x, std::complex<T> *y,
                            std::complex<T> *alpha, bool *marked, int end, int blen, Level *l){
  int nv = l->num_lattice_site_var;
  int num_lattice_sites = end/nv;
  bool flag = false;

  for (int j = 0; j < blen; j++) {
    if (marked[j]){
      flag = true;
    }
  }

  if(flag){
    for (int i = 0; i < num_lattice_sites; i++) {
      for (int j = 0; j < blen; j++) {
        if (marked[j]){
          continue;
        }
        for (int k = 0; k < nv; k++) {
          z[i*blen*nv + j*nv + k] = x[i*blen*nv + j*nv + k] + alpha[j] * y[i*blen*nv + j*nv + k];
        }
      }
    }
  } else{
    for (int i = 0; i < num_lattice_sites; i++) {
      for (int j = 0; j < blen; j++) {
        for (int k = 0; k < nv; k++) {
          z[i*blen*nv + j*nv + k] = x[i*blen*nv + j*nv + k] + alpha[j] * y[i*blen*nv + j*nv + k];
        }
      }
    }
  }
}

template<typename T>
void vector_saxpy_mark_per_entry(std::complex<T> *z, std::complex<T> *x, std::complex<T> *y,
                                 std::complex<T> *alpha, bool *marked, int end, int blen, Level *l){
  bool flag = false;

  for (int j = 0; j < blen; j++) {
    if (marked[j]){
      flag = true;
    }
  }

  if (flag){
    for (int i = 0; i < end; i++) {
      for (int j = 0; j < blen; j++) {
        if (marked[j]){
          continue;
        }
        z[i*blen + j] = x[i*blen + j] + alpha[j] * y[i*blen + j];
      }
    }
  } else{
    for (int i = 0; i < end; i++) {
      for (int j = 0; j < blen; j++) {
        z[i*blen + j] = x[i*blen + j] + alpha[j] * y[i*blen + j];
      }
    }
  }
}

/** z = alpha*x
 * */
template<typename T>
void vector_scale(std::complex<T> *z, std::complex<T> *x, std::complex<T> *alpha, int end, int blen,
                  Level *l) {
  int nv = l->num_lattice_site_var;
  int num_lattice_sites = end/nv;

  for (int i = 0; i < num_lattice_sites; i++) {
    for (int j = 0; j < blen; j++) {
      for (int k = 0; k < nv; k++) {
        z[i*blen*nv + j*nv + k] = alpha[j] * x[i*blen*nv + j*nv + k];
      }
    }
  }
}

template<typename T>
void vector_scale_per_entry(std::complex<T> *z, std::complex<T> *x, std::complex<T> *alpha, int end, int blen,
                            Level *l){
  for (int i = 0; i < end; i++) {
    for (int j = 0; j < blen; j++) {
      z[i*blen + j] = alpha[j] * x[i*blen + j];
    }
  }
}

/** res = phi[i]_H * psi[i], i in 0...nv-1
 * */
template<typename T>
static inline std::complex<T> segment_dot_product(int nv, complex<T> *phi, complex<T> *psi){
  std::complex<T> res = 0;
  for (int i = 0; i < nv; i++) {
    res += conj(phi[i]) * psi[i];
  }
  return res;
}

template<typename T>
void multi_inner_product(int count, std::complex<T> *results, std::complex<T> **phi,
                         std::complex<T> *psi, int end, int blen, Level *l) {
  std::fill(results, results+blen*count, 0);
  int nv = l->num_lattice_site_var;
  int num_lattice_sites = end/nv;

  for (int i = 0; i < num_lattice_sites; i++) {
    for (int c = 0; c < count; c++) {
      for (int j = 0; j < blen; j++) {
        results[blen*c + j] +=
            segment_dot_product(nv, &phi[c][i * nv * blen + j * nv], &psi[i * nv * blen + j * nv]);
      }
    }
  }
}

template<typename T>
void multi_inner_product_per_entry(int count, std::complex<T> *results, std::complex<T> **phi,
                                   std::complex<T> *psi, int end, int blen, Level *l) {
  std::fill(results, results+blen*count, 0);

  for (int i = 0; i < end; i++) {
    for (int j = 0; j < count; j++) {
      for (int k = 0; k < blen; k++) {
        results[blen*j + k] += conj(phi[j][i*blen + k]) * psi[i*blen + k];
      }
    }
  }
}

template<typename T>
void BlockGmres<T>::arnoldiStep(complex<T> **V, complex<T> *H, complex<T> *w, complex<T> *y,
                                int curLoop,
                                void (*evalOp)(BlockVecPerSite<T>*, BlockVecPerSite<T>*,
                                               operator_struct<T>*, Level*, const int)) {
  Global *global = this->level->global;
  int vlen = this->vectorEnd;
  int blen = this->blockLength;

  auto *blockv = new BlockVecPerSite<T>(V[curLoop], 12);
  auto *blockw = new BlockVecPerSite<T>(w, 12);

  evalOp(blockw, blockv, this->matrix, this->level, blen);

  /// orthogonalization
  /// y[c] = conj(V[c][i]) * w[i], where c in 0...j and i in vectorStart..vectorEnd
  multi_inner_product(curLoop+1, y, V, w, vlen, blen, this->level);

  /// H[j][i] = <V[i]_H, w>
  if (global->num_processes > 1) {
    //std::cout << H[get_H_idx(curLoop,0,1,blen)] << std::endl;

    MPI_Allreduce(y, &H[get_H_idx(curLoop,0,0,blen)], blen*(curLoop+1), MPI_COMPLEX_double, MPI_SUM,
                  (this->level->depth == 0) ? global->comm_cart : this->level->gs_double.level_comm);
  } else {
    for (int i = 0; i < curLoop + 1; i++){
      for (int j = 0; j < blen; j++) {
        H[get_H_idx(curLoop, i, j, blen)] = y[blen*i + j];
      }
    }
  }

  /// w = w-H[j][i]*V[i]
  for (int i = 0; i <= curLoop; i++)
    vector_saxpy_neg(w, w, V[i], &H[get_H_idx(curLoop, i, 0, blen)], vlen, blen, this->level);

  /// H[j][j+1] = ||w||_2
  set_global_norm(w, &H[get_H_idx(curLoop, curLoop+1, 0, blen)], vlen, blen, this->level);

  auto *alpha = new std::complex<T>[blen];
  /// V[j+1] = w/H[j][j+1]
  for (int i = 0; i < blen; i++) {
    //std::cout << abs(H[get_H_idx(curLoop, curLoop+1, 0, blen)+i]) << std::endl;
    if (abs(H[get_H_idx(curLoop, curLoop+1, 0, blen)+i]) > 1e-15){
      alpha[i] = real(1.0/H[get_H_idx(curLoop, curLoop+1, 0, blen)+i]);
    } else{
      alpha[i] = 1.0;
    }
  }
  vector_scale(V[curLoop+1], w, alpha, vlen, blen, this->level);

  delete[] alpha;
}

template <typename T>
void BlockGmres<T>::arnoldiStep(complex<T> **V, complex<T> *H, complex<T> *w, complex<T> *y,
                                int curLoop,
                                void (*evalOp)(BlockVecPerEntry<T> *, BlockVecPerEntry<T> *,
                                               operator_struct<T> *, Level *, const int)) {
  Global *global = this->level->global;
  int vlen = this->vectorEnd;
  int blen = this->blockLength;

  auto *blockv = new BlockVecPerEntry<T>(V[curLoop]);
  auto *blockw = new BlockVecPerEntry<T>(w);

  evalOp(blockw, blockv, this->matrix, this->level, blen);

  /// orthogonalization
  /// y[c] = conj(V[c][i]) * w[i], where c in 0...j and i in vectorStart..vectorEnd
  multi_inner_product_per_entry(curLoop+1, y, V, w, vlen, blen, this->level);

  /// H[j][i] = <V[i]_H, w>
  if (global->num_processes > 1) {
    MPI_Allreduce(y, &H[get_H_idx(curLoop,0,0,blen)], blen*(curLoop+1), MPI_COMPLEX_double, MPI_SUM,
                  (this->level->depth == 0) ? global->comm_cart : this->level->gs_double.level_comm);
  } else {
    for (int i = 0; i < curLoop + 1; i++){
      for (int j = 0; j < blen; j++) {
        H[get_H_idx(curLoop, i, j, blen)] = y[blen*i + j];
      }
    }
  }

  /// w = w-H[j][i]*V[i]
  for (int i = 0; i <= curLoop; i++)
    vector_saxpy_neg_per_entry(w, w, V[i], &H[get_H_idx(curLoop, i, 0, blen)], vlen, blen, this->level);

  /// H[j][j+1] = ||w||_2
  set_global_norm_per_entry(w, &H[get_H_idx(curLoop, curLoop+1, 0, blen)], vlen, blen, this->level);

  auto *alpha = new std::complex<T>[blen];
  /// V[j+1] = w/H[j][j+1]
  for (int i = 0; i < blen; i++) {
    if (abs(H[get_H_idx(curLoop, curLoop+1, 0, blen)+i]) > 1e-15){
      alpha[i] = real(1.0/H[get_H_idx(curLoop, curLoop+1, 0, blen)+i]);
    } else{
      alpha[i] = 1.0;
    }
  }
  vector_scale_per_entry(V[curLoop+1], w, alpha, vlen, blen, this->level);

  delete[] alpha;
}

template<typename T>
void BlockGmres<T>::qrUpdate(std::complex<T> *H, std::complex<T> *s, std::complex<T> *c,
                             std::complex<T> *gamma, int curLoop){
  const int blen = this->blockLength;
  std::complex<T> beta[blen];
  for (int i = 0; i < curLoop; i++) {
    for (int j = 0; j < blen; j++) {
      beta[j] = (-s[i*blen + j]) * H[get_H_idx(curLoop, i, j, blen)] +
                (c[i*blen + j]) * H[get_H_idx(curLoop, i+1, j, blen)];
      H[get_H_idx(curLoop, i, j, blen)] = conj(c[i*blen + j]) * H[get_H_idx(curLoop, i, j, blen)] +
                                          conj(s[i*blen + j]) * H[get_H_idx(curLoop, i+1, j, blen)];
      H[get_H_idx(curLoop, i+1, j, blen)] = beta[j];
    }
  }

  for (int i = 0; i < blen; i++) {
    beta[i] = (std::complex<T>)sqrt(NORM_SQUARE_double(H[get_H_idx(curLoop, curLoop, i, blen)]) +
                                    NORM_SQUARE_double(H[get_H_idx(curLoop, curLoop+1, i, blen)]));
    s[curLoop*blen + i] = H[get_H_idx(curLoop, curLoop+1, i, blen)]/beta[i];
    c[curLoop*blen + i] = H[get_H_idx(curLoop, curLoop, i, blen)]/beta[i];
    /// update right column
    gamma[(curLoop+1)*blen + i] = (-s[curLoop*blen + i]) * gamma[curLoop*blen + i];
    gamma[curLoop*blen + i] = conj(c[curLoop*blen + i]) * gamma[curLoop*blen + i];
    /// apply current Givens rotation
    H[get_H_idx(curLoop, curLoop, i, blen)] = beta[i];
    H[get_H_idx(curLoop, curLoop+1, i, blen)] = 0;
  }

}

template<typename T>
void BlockGmres<T>::computeSolution(std::complex<T> **V, std::complex<T> *y,std::complex<T> *gamma,
                                    std::complex<T> *H, bool *marked, int curLoop, int outerLoop,
                                    int layout){
  const int blen = this->blockLength;
  int vlen = this->vectorEnd;

  /// backward substitution
  for (int i = curLoop; i >= 0 ; i--) {
    for (int j = 0; j < blen; j++) {
      y[i*blen + j] = gamma[i*blen + j];
    }

    for (int k = i + 1; k <= curLoop; k++) {
      for (int j = 0; j < blen; j++) {
        y[i*blen + j] -= H[get_H_idx(k, i, j, blen)] * y[k*blen + j];
      }
    }

    for (int j = 0; j < blen; j++) {
      y[i*blen + j] /= H[get_H_idx(i, i, j, blen)];
    }
  }

  /// x = x + V*y
  if (outerLoop) {
    if (layout == _PER_SITE){
      for (int i = 0; i <= curLoop; i++) {
        vector_saxpy_with_mark(this->x, this->x, V[i], &y[blen * i], marked, vlen, blen, this->level);
      }
    } else if (layout == _PER_ENTRY){
      for (int i = 0; i <= curLoop; i++) {
        vector_saxpy_mark_per_entry(this->x, this->x, V[i], &y[blen * i], marked, vlen, blen, this->level);
      }
    }
  } else {
    if (layout == _PER_SITE){
      vector_scale(this->x, V[0], &y[0], vlen, blen, this->level);
      for (int i = 1; i <= curLoop; i++) {
        vector_saxpy_with_mark(this->x, this->x, V[i], &y[blen * i], marked, vlen, blen, this->level);
      }
    } else if (layout == _PER_ENTRY){
      vector_scale_per_entry(this->x, V[0], &y[0], vlen, blen, this->level);
      for (int i = 1; i <= curLoop; i++) {
        vector_saxpy_mark_per_entry(this->x, this->x, V[i], &y[blen * i], marked, vlen, blen, this->level);
      }
    }

  }
}

template<typename T>
void BlockGmres<T>::solve(void (*evalOp)(BlockVecPerSite<T>*, BlockVecPerSite<T>*,
                                         operator_struct<T>*, Level*, const int)) {
  Global *global = this->level->global;
  int restartLen = this->restartLen;
  const int blen = this->blockLength;
  int vlen = this->vectorEnd;
  int iter = 0;
  int j = -1;

  auto *H =
      new std::complex<T>[(((restartLen + 1) * (restartLen + 1) + 3 * (restartLen + 1)) / 2) *
                          blen];
  std::complex<T> **V;
  NEW2D(V, std::complex<T>, blen * vlen, restartLen + 1);

  auto *w = new std::complex<T>[vlen * blen];
  auto *y = new std::complex<T>[(restartLen+1) * blen];
  auto *s = new std::complex<T>[blen * restartLen];
  auto *c = new std::complex<T>[blen * restartLen];
  auto *gamma = new std::complex<T>[(restartLen+1) * blen];

  auto *norm = new std::complex<T>[blen];
  auto *normR0 = new std::complex<T>[blen];
  auto *gammaJPlusOne = new std::complex<T>[blen];

  auto *alpha = new std::complex<T>[blen];
  auto *relnorm = new T[blen];
  auto *marked = new bool[blen];
  memset(marked, 0, blen * sizeof(marked[0]));

  /// wrapper for the data layout
  auto *blockx = new BlockVecPerSite<T>(this->x, 12);
  auto *blockw = new BlockVecPerSite<T>(w, 12);
  //auto *blockb = new BlockVecPerSite<T>(this->b, 12);
  //blockb->transform(vlen, blen);

  if (this->level->type == _FINE)
    elapsedTime = -MPI_Wtime();

  if (this->verbose && global->print > 0)
    printf0("\n+----------------------------------------------------------+\n");
  for (int outerLoop = 0; outerLoop < this->restartNum && !this->finished; outerLoop++) {
    if (outerLoop == 0 && this->initialGuessIsZero) {
      /// set residual flag to no residual
      res = _NO_RES;
      /// set r = b
      vector_copy(this->r, this->b, 0, blen * vlen);
    } else {
      res = _RES;
      /// w = Ax
      evalOp(blockw, blockx, this->matrix, this->level, blen);
      /// r = b - w
      vector_minus(this->r, this->b, w, 0, blen * vlen);
    }

    /// norm = ||r||_2 Euclidean norm of r
    set_global_norm(this->r, &norm[0], vlen, blen, this->level);
    vector_copy(&gamma[0], &norm[0], 0, blen);

    if (outerLoop == 0) {
      if (this->level->type == _FINE && !this->initialGuessIsZero) {
        ///normR0 = ||b||_2
        set_global_norm(this->b, &normR0[0], vlen, blen, this->level);
      } else {
        vector_copy(&normR0[0], &norm[0], 0, blen);
      }
    }

    /// V[0] = 1/gamma[0] * r
    for (int i = 0; i < blen; i++) {
      alpha[i] = 1.0 / gamma[i];
    }
    vector_scale(V[0], this->r, alpha, vlen, blen, this->level);
    //vector_copy(inject, V[0], 0, blen * vlen);

    for (int innerLoop = 0; innerLoop < restartLen && !this->finished; innerLoop++) {
      j = innerLoop;
      iter++;

      /// simply do one step of Arnoldi
      arnoldiStep(V, H, w, y, j, evalOp);

      /// check for every single rhs to see if it stagnates
      check_staganation(&H[get_H_idx(j,j+1, 0, blen)], marked, this->tolerance / 10, blen);

      /// do QR nevertheless
      qrUpdate(H, s, c, gamma, j);
      for (int i = 0; i < blen; ++i) {
        gammaJPlusOne[i] = abs(gamma[blen*(j+1)+i]);
      }

      for (int i = 0; i < blen; ++i) {
        relnorm[i] = real(gammaJPlusOne[i]) / real(normR0[i]);
      }

      /// printing
      if (iter % 10 == 0 || this->level->depth > 0) {
        if (this->verbose && global->print > 0)
          for (int i = 0; i < blen; ++i) {
            printf0("| approx. rel. res. after  %-6d iterations: %e |\n", iter,
                    relnorm[i]);
          }
      }

      /// here it can go on below tolerance
      T maxnorm = max_in_block(relnorm, blen);
      if (maxnorm < this->tolerance || maxnorm > 1E+5) { // if satisfied ... stop
        this->finished = 1;
        if (maxnorm > 1E+5)
          printf0("Divergence of fgmres_double, iter = %d, depth=%d\n", iter, this->level->depth);
      }

    } // end of a single restart
    //vector_copy(inject, H, 0, blen * ((((restartLen + 1) * (restartLen + 1) + 3 * (restartLen + 1)) / 2)));
    computeSolution(V, y, gamma, H, marked, j, (res == _NO_RES) ? outerLoop : 1, _PER_SITE);

  }

  if (this->level->type == _FINE) {
    elapsedTime += MPI_Wtime();
    global->total_time = elapsedTime;
    global->iter_count = iter;
    global->norm_res = max_in_block(relnorm, blen);
  }

  if (this->verbose) {
#ifdef FGMRES_RESTEST
    evalOp(blockw, blockx, this->matrix, this->level, blen);
    vector_minus(this->r, this->b, w, 0, blen * vlen);
    set_global_norm(this->r, &norm[0], vlen, blen, this->level);
#else
    vector_copy(&norm[0], &gammaJPlusOne[0], 0, blen);
#endif
    T normres = real(norm[0])/real(normR0[0]);
    for (int i = 1; i < blen; ++i) {
      if(real(norm[i])/real(normR0[i]) > normres){
        normres = real(norm[i])/real(normR0[i]);
      }
    }
    global->norm_res = normres;
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    if (global->print > 0)

      printf0("+----------------------------------------------------------+\n\n");
#endif
    printf0("+----------------------------------------------------------+\n");
    printf0("|       FGMRES iterations: %-6d coarse average: 0        |\n", iter); // TODO Fix coarse average
    printf0("| exact relative residual: ||r||/||b|| = %e      |\n", global->norm_res);
    printf0("| elapsed wall clock time: %-8.4lf seconds                |\n", elapsedTime);
    if (global->coarse_time > 0)
      printf0("|        coarse grid time: %-8.4lf seconds (%04.1lf%%)        |\n",
              global->coarse_time, 100 * (global->coarse_time / elapsedTime));
    printf0("|      coarsest grid time: %-8.4lf seconds (%04.1lf%%)        |\n",
            global->coarsest_time, 100 * (global->coarsest_time / elapsedTime));
    printf0("|  consumed core minutes: %-8.2le (solve only)            |\n",
            (elapsedTime * global->num_processes / 60.0));
    printf0("|    max used mem/MPIproc: %-8.2le GB                     |\n",
            0.0); // TODO: Fix mem usage(?)
    printf0("+----------------------------------------------------------+\n");
  }

  if (this->level->type == _COARSEST) {
    global->coarse_iter_count += iter;
  }

  /// return with x transformed back to sequential block store
  blockx->detransform(vlen, blen);

  /// clean up / reset relevant class members
  this->finished = 0;

  delete[] H;
  delete[] w;
  delete[] y;
  delete[] s;
  delete[] c;
  delete[] gamma;
  delete[] norm;
  delete[] normR0;
  delete[] alpha;
  delete[] relnorm;
  delete[] marked;
  DELETE2D(V);
}

template <typename T>
void BlockGmres<T>::solve(void (*evalOp)(BlockVecPerEntry<T> *, BlockVecPerEntry<T> *,
                                         operator_struct<T> *, Level *, const int)) {
  Global *global = this->level->global;
  int restartLen = this->restartLen;
  const int blen = this->blockLength;
  int vlen = this->vectorEnd;
  int iter = 0;
  int j = -1;

  auto *H =
      new std::complex<T>[(((restartLen + 1) * (restartLen + 1) + 3 * (restartLen + 1)) / 2) *
                          blen];
  std::complex<T> **V;
  NEW2D(V, std::complex<T>, blen * vlen, restartLen + 1);

  auto *w = new std::complex<T>[vlen * blen];
  auto *y = new std::complex<T>[(restartLen+1) * blen];
  auto *s = new std::complex<T>[blen * restartLen];
  auto *c = new std::complex<T>[blen * restartLen];
  auto *gamma = new std::complex<T>[(restartLen+1) * blen];

  auto *norm = new std::complex<T>[blen];
  auto *normR0 = new std::complex<T>[blen];
  auto *gammaJPlusOne = new std::complex<T>[blen];

  auto *alpha = new std::complex<T>[blen];
  auto *relnorm = new T[blen];
  auto *marked = new bool[blen];
  memset(marked, 0, blen * sizeof(marked[0]));

  /// wrapper for the data layout
  auto *blockx = new BlockVecPerEntry<T>(this->x);
  auto *blockw = new BlockVecPerEntry<T>(w);
  //auto *blockb = new BlockVecPerEntry<T>(this->b);
  //blockb->transform(vlen, blen);

  if (this->level->type == _FINE)
    elapsedTime = -MPI_Wtime();

  if (this->verbose && global->print > 0)
    printf0("\n+----------------------------------------------------------+\n");
  for (int outerLoop = 0; outerLoop < this->restartNum && !this->finished; outerLoop++) {
    if (outerLoop == 0 && this->initialGuessIsZero) {
      /// set residual flag to no residual
      res = _NO_RES;
      /// set r = b
      vector_copy(this->r, this->b, 0, blen * vlen);
    } else {
      res = _RES;
      /// w = Ax
      evalOp(blockw, blockx, this->matrix, this->level, blen);
      /// r = b - w
      vector_minus(this->r, this->b, w, 0, blen * vlen);
    }

    /// norm = ||r||_2 Euclidean norm of r
    set_global_norm_per_entry(this->r, &norm[0], vlen, blen, this->level);
    vector_copy(&gamma[0], &norm[0], 0, blen);

    if (outerLoop == 0) {
      if (this->level->type == _FINE && !this->initialGuessIsZero) {
        ///normR0 = ||b||_2
        set_global_norm_per_entry(this->b, &normR0[0], vlen, blen, this->level);
      } else {
        vector_copy(&normR0[0], &norm[0], 0, blen);
      }
    }

    /// V[0] = 1/gamma[0] * r
    for (int i = 0; i < blen; i++) {
      alpha[i] = 1.0 / gamma[i];
    }
    vector_scale_per_entry(V[0], this->r, alpha, vlen, blen, this->level);
    //vector_copy(inject, V[0], 0, blen * vlen);

    for (int innerLoop = 0; innerLoop < restartLen && !this->finished; innerLoop++) {
      j = innerLoop;
      iter++;

      /// simply do one step of Arnoldi
      arnoldiStep(V, H, w, y, j, evalOp);

      /// check for every single rhs to see if it stagnates
      check_staganation(&H[get_H_idx(j,j+1, 0, blen)], marked, this->tolerance / 10, blen);

      /// do QR nevertheless
      qrUpdate(H, s, c, gamma, j);
      for (int i = 0; i < blen; ++i) {
        gammaJPlusOne[i] = abs(gamma[blen*(j+1)+i]);
      }

      for (int i = 0; i < blen; ++i) {
        relnorm[i] = real(gammaJPlusOne[i]) / real(normR0[i]);
      }

      /// printing
      if (iter % 10 == 0 || this->level->depth > 0) {
        if (this->verbose && global->print > 0)
          for (int i = 0; i < blen; ++i) {
            printf0("| approx. rel. res. after  %-6d iterations: %e |\n", iter,
                    relnorm[i]);
          }
      }

      /// here it can go on below tolerance
      T maxnorm = max_in_block(relnorm, blen);
      if (maxnorm < this->tolerance || maxnorm > 1E+5) { // if satisfied ... stop
        this->finished = 1;
        if (maxnorm > 1E+5)
          printf0("Divergence of fgmres_double, iter = %d, depth=%d\n", iter, this->level->depth);
      }

    } // end of a single restart
    //vector_copy(inject, H, 0, blen * ((((restartLen + 1) * (restartLen + 1) + 3 * (restartLen + 1)) / 2)));
    computeSolution(V, y, gamma, H, marked, j, (res == _NO_RES) ? outerLoop : 1, _PER_ENTRY);

  }

  if (this->level->type == _FINE) {
    elapsedTime += MPI_Wtime();
    global->total_time = elapsedTime;
    global->iter_count = iter;
    global->norm_res = max_in_block(relnorm, blen);
  }

  if (this->verbose) {
#ifdef FGMRES_RESTEST
    evalOp(blockw, blockx, this->matrix, this->level, blen);
    vector_minus(this->r, this->b, w, 0, blen * vlen);
    set_global_norm_per_entry(this->r, &norm[0], vlen, blen, this->level);
#else
    vector_copy(&norm[0], &gammaJPlusOne[0], 0, blen);
#endif
    T normres = real(norm[0])/real(normR0[0]);
    for (int i = 1; i < blen; ++i) {
      if(real(norm[i])/real(normR0[i]) > normres){
        normres = real(norm[i])/real(normR0[i]);
      }
    }
    global->norm_res = normres;
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    if (global->print > 0)

      printf0("+----------------------------------------------------------+\n\n");
#endif
    printf0("+----------------------------------------------------------+\n");
    printf0("|       FGMRES iterations: %-6d coarse average: 0        |\n", iter); // TODO Fix coarse average
    printf0("| exact relative residual: ||r||/||b|| = %e      |\n", global->norm_res);
    printf0("| elapsed wall clock time: %-8.4lf seconds                |\n", elapsedTime);
    if (global->coarse_time > 0)
      printf0("|        coarse grid time: %-8.4lf seconds (%04.1lf%%)        |\n",
              global->coarse_time, 100 * (global->coarse_time / elapsedTime));
    printf0("|      coarsest grid time: %-8.4lf seconds (%04.1lf%%)        |\n",
            global->coarsest_time, 100 * (global->coarsest_time / elapsedTime));
    printf0("|  consumed core minutes: %-8.2le (solve only)            |\n",
            (elapsedTime * global->num_processes / 60.0));
    printf0("|    max used mem/MPIproc: %-8.2le GB                     |\n",
            0.0); // TODO: Fix mem usage(?)
    printf0("+----------------------------------------------------------+\n");
  }

  if (this->level->type == _COARSEST) {
    global->coarse_iter_count += iter;
  }

  /// return with x transformed back to sequential block store
  blockx->detransform(vlen, blen);

  /// clean up / reset relevant class members
  this->finished = 0;

  delete[] H;
  delete[] w;
  delete[] y;
  delete[] s;
  delete[] c;
  delete[] gamma;
  delete[] norm;
  delete[] normR0;
  delete[] alpha;
  delete[] relnorm;
  delete[] marked;
  DELETE2D(V);
}



#endif // DDALPHAAMG_BLOCK_GMRES_SOLVE_H
