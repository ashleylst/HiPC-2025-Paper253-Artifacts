#ifndef DDALPHAAMG_ODDEVEN_GENERIC_TEMPLATE_H
#define DDALPHAAMG_ODDEVEN_GENERIC_TEMPLATE_H

template<typename T>
static inline void LLH_perform_fwd_bwd_subs(std::complex<T>* x, std::complex<T>* b,
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

template<typename T>
static inline void LLH_multiply(std::complex<T>* y, std::complex<T>* x, std::complex<T>* L) {

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

template<typename T>
void diag_ee(std::complex<T>* y, std::complex<T>* x, operator_struct<T> *op, Level *l,
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
      LLH_multiply(y, x, sc);
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

template<typename T>
void diag_oo_inv(std::complex<T>* y, std::complex<T>* x, operator_struct<T> *op, Level *l,
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
      LLH_perform_fwd_bwd_subs(y, x, sc);
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


#endif // DDALPHAAMG_ODDEVEN_GENERIC_TEMPLATE_H
