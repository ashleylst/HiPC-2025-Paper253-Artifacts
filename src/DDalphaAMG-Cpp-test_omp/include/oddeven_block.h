#ifndef DDALPHAAMG_ODDEVEN_BLOCK_H
#define DDALPHAAMG_ODDEVEN_BLOCK_H

/** @Brief
 * Solves L*(L^H)*x = b for x, i.e., the clover coupling for a single lattice site.
 * @tparam T Precision
 * @param x Output vector.
 * @param b Input vector.
 * @param L Cholesky factor ( lower triangular matrix ).
 * @param blen block size.
 */
template<typename T>
static inline void LLH_perform_fwd_bwd_subs_block(std::complex<T>* x, std::complex<T>* b,
                                                  std::complex<T>* L, const int blen) {

  for (int n = 0; n < 2; n++) {
    /// forward substitution with L
    for (int i = 0; i < 6; i++) {
      for (int k = 0; k < blen; k++) {
        x[(i + 6*n)*blen + k] = b[(i + 6*n)*blen + k];
        for (int j = 0; j < i; j++) {
          x[(i + 6*n)*blen + k] -= L[(i*i+i)/2 + j + 21*n] * x[(j + 6*n)*blen + k];
        }
        x[(i + 6*n)*blen + k] /= L[(i*i+i)/2 + i + 21*n];
      }
    }

    /// backward substitution with L^H
    for (int i = 5; i >= 0; i--) {
      for (int k = 0; k < blen; k++) {
        for (int j = i + 1; j < 6; j++) {
          x[(i + 6*n)*blen + k] -= conj(L[(j*j+j)/2 + i + 21*n]) * x[(j + 6*n)*blen + k];
        }
        x[(i + 6*n)*blen + k] /= conj(L[(i*i+i)/2 + i + 21*n]);
      }
    }
  }
}

/** @Brief
 * Apply Clover to x, i.e., y = L(L^H)x.
 * @tparam T Precision
 * @param y Output vector.
 * @param x Input vector.
 * @param L Cholesky factor ( lower triangular matrix ).
 * @param blen block size.
 */
template<typename T>
static inline void LLH_multiply_block(std::complex<T>* y, std::complex<T>* x, std::complex<T>* L,
                                      const int blen) {
  auto *z = new std::complex<T>[6*blen];

  for (int n = 0; n < 2; n++) {
    /// z = L^H x
    for (int j = 0; j < 6; j++) {   // columns
      for (int k = 0; k < blen; k++) {
        for (int i = 0; i < j; i++) { // rows
          z[i*blen + k] += conj(L[(j*j+j)/2 + i + 21*n]) * x[(j + 6*n)*blen + k];
        }
        z[j*blen + k] = conj(L[(j*j+j)/2 + j + 21*n]) * x[(j + 6*n)*blen + k];
      }
    }
    /// y = L*z;
    for (int i = 0; i < 6; i++) { // rows
      for (int k = 0; k < blen; k++) {
        y[(i + 6*n)*blen + k] = L[(i*i+i)/2 + 21*n] * z[k];
        for (int j = 1; j <= i; j++) { // columns
          y[(i + 6*n)*blen + k] += L[(i*i+i)/2 + j + 21*n] * z[j*blen + k];
        }
      }
    }
  }

  delete[] z;
}

/** @Brief
 * Apply the inverse of the clover to in, i.e., out = D^-1(in).
 * @tparam T Precision.
 * @param out Output vector.
 * @param in Input vector.
 * @param op Operator.
 * @param l Level.
 * @param start Vector starting position.
 * @param end Vector ending position
 * @param blen Block size.
 */
template<typename T>
void diag_oo_inv_block(BlockVecPerEntry<T> *out, BlockVecPerEntry<T> *in, operator_struct<T> *op,
                       Level *l, int start, int end, const int blen){
  Global *global = l->global;
  std::complex<T>* sc = op->clover;

  std::complex<T> *x = in->v;
  std::complex<T> *y = out->v;

  if (global->csw){
    for (int i = start/12; i < end/12; i++) {
      LLH_perform_fwd_bwd_subs_block(&y[i*12*blen], &x[i*12*blen], &sc[i*42], blen);
    }
  } else {
    for (int i = start/12; i < end/12; i++) {
      for (int j = 0; j < 12; j++) {
        for (int k = 0; k < blen; k++){
          y[(i*12 + j)*blen + k] = x[(i*12 + j)*blen + k] / sc[i*12 + j];
        }
      }
    }
  }
}

template<typename T>
void diag_oo_inv_block(BlockVecPerSite<T> *out, BlockVecPerSite<T> *in, operator_struct<T> *op,
                       Level *l, int start, int end, const int blen){
  Global *global = l->global;
  std::complex<T>* sc = op->clover;

  std::complex<T> *x = in->v;
  std::complex<T> *y = out->v;

  if (global->csw){
    for (int i = start/12; i < end/12; i++) {
      for (int j = 0; j < blen; j++) {
        LLH_perform_fwd_bwd_subs(&y[i * 12 * blen + j * 12], &x[i * 12 * blen + j * 12],
                                 &sc[i * 42]);
      }
    }
  } else {
    for (int i = start/12; i < end/12; i++) {
      for (int j = 0; j < blen; j++) {
        for (int k = 0; k < 12; k++){
          y[i*12*blen + j*12 + k] = x[i*12*blen + j*12 + k] / sc[i*12 + k];
        }
      }
    }
  }
}

/** @Brief
 * Apply the clover to in, i.e., out = D(in)
 * @tparam T Precision.
 * @param out Output vector.
 * @param in Input vector.
 * @param op Operator.
 * @param l Level.
 * @param start Vector starting position.
 * @param end Vector ending position.
 * @param blen Block size.
 */
template<typename T>
void diag_ee_block(BlockVecPerEntry<T>* out, BlockVecPerEntry<T>* in, operator_struct<T> *op,
                   Level *l, int start, int end, const int blen) {
  Global *global = l->global;
  std::complex<T>* sc = op->clover;

  std::complex<T> *x = in->v;
  std::complex<T> *y = out->v;
  if (global->csw) {
    /// diagonal blocks applied to the even sites
    for (int i = start/12; i < end/12; i++) {
      LLH_multiply_block(&y[i*12*blen], &x[i*12*blen], &sc[i*42], blen);
    }
  } else {
    for (int i = start/12; i < end/12; i++) {
      for (int j = 0; j < 12; j++) {
        for (int k = 0; k < blen; k++){
          y[(i*12 + j)*blen + k] = x[(i*12 + j)*blen + k] * sc[i*12 + j];
        }
      }
    }
  }
}

template<typename T>
void diag_ee_block(BlockVecPerSite<T>* out, BlockVecPerSite<T>* in, operator_struct<T>* op,
                   Level *l, int start, int end, const int blen){
  Global *global = l->global;
  std::complex<T>* sc = op->clover;

  std::complex<T> *x = in->v;
  std::complex<T> *y = out->v;
  if (global->csw) {
    /// diagonal blocks applied to the even sites
    for (int i = start/12; i < end/12; i++) {
      for (int j = 0; j < blen; j++) {
        LLH_multiply(&y[i*12*blen + j*12], &x[i*12*blen + j*12], &sc[i*42]);
      }
    }
  } else {
    for (int i = start/12; i < end/12; i++) {
      for (int j = 0; j < blen; j++) {
        for (int k = 0; k < 12; k++){
          y[i*12*blen + j*12 + k] = x[i*12*blen + j*12 + k] * sc[i*12 + k];
        }
      }
    }
  }
}

template<typename T>
void hopping_term_block(BlockVecPerEntry<T> *bout, BlockVecPerEntry<T> *bv,
                        operator_struct<T> *op, const int amount, Level *l, const int blen) {

  int start_even, end_even, start_odd, end_odd,
      n = l->num_inner_lattice_sites, *neighbor = op->neighbor_table, start = 0,
      plus_dir_param = _FULL_SYSTEM, minus_dir_param = _FULL_SYSTEM;
  std::complex<T> *phi = bv->v;
  std::complex<T> *eta = bout->v;

  auto *prnT = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];
  auto *prnZ = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];
  auto *prnY = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];
  auto *prnX = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];

  auto *prpT = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];
  auto *prpZ = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];
  auto *prpY = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];
  auto *prpX = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];

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

  pi_minus_per_entry(phi, prnT, prnZ, prnY, prnX, start*12, n*12, blen);

  /// start communication in negative direction
  ghost_sendrecv_per_entry(prnT, _T, -1, &(op->c), minus_dir_param, l, blen);
  ghost_sendrecv_per_entry(prnZ, _Z, -1, &(op->c), minus_dir_param, l, blen);
  ghost_sendrecv_per_entry(prnY, _Y, -1, &(op->c), minus_dir_param, l, blen);
  ghost_sendrecv_per_entry(prnX, _X, -1, &(op->c), minus_dir_param, l, blen);

  /// project plus dir and multiply with U dagger
  pi_plus_per_entry(phi, op->D, prpT, prpZ, prpY, prpX, neighbor, start*12, n*12, blen);

  if (amount == _EVEN_SITES) {
    start = start_even, n = end_even;
  } else if (amount == _ODD_SITES) {
    start = start_odd, n = end_odd;
  }

  /// start communication in positive direction
  ghost_sendrecv_per_entry(prpT, _T, +1, &(op->c), plus_dir_param, l, blen);
  ghost_sendrecv_per_entry(prpZ, _Z, +1, &(op->c), plus_dir_param, l, blen);
  ghost_sendrecv_per_entry(prpY, _Y, +1, &(op->c), plus_dir_param, l, blen);
  ghost_sendrecv_per_entry(prpX, _X, +1, &(op->c), plus_dir_param, l, blen);

  /// wait for communication in negative direction
  ghost_wait_per_entry(prnT, _T, -1, &(op->c), minus_dir_param, l, blen);
  ghost_wait_per_entry(prnZ, _Z, -1, &(op->c), minus_dir_param, l, blen);
  ghost_wait_per_entry(prnY, _Y, -1, &(op->c), minus_dir_param, l, blen);
  ghost_wait_per_entry(prnX, _X, -1, &(op->c), minus_dir_param, l, blen);

  /// multiply with U and lift up minus dir
  apply_U_per_entry(eta, op->D, prnT, prnZ, prnY, prnX, neighbor, start*12, n*12, blen);

  /// wait for communication in positive direction
  ghost_wait_per_entry(prpT, _T, +1, &(op->c), plus_dir_param, l, blen);
  ghost_wait_per_entry(prpZ, _Z, +1, &(op->c), plus_dir_param, l, blen);
  ghost_wait_per_entry(prpY, _Y, +1, &(op->c), plus_dir_param, l, blen);
  ghost_wait_per_entry(prpX, _X, +1, &(op->c), plus_dir_param, l, blen);

  /// lift up plus dir
  lift_pi_plus_per_entry(eta, prpT, prpZ, prpY, prpX, start*12, n*12, blen);

  delete[] prnT;
  delete[] prnZ;
  delete[] prnY;
  delete[] prnX;

  delete[] prpT;
  delete[] prpZ;
  delete[] prpY;
  delete[] prpX;
}

template<typename T>
void hopping_term_block(BlockVecPerSite<T> *bout, BlockVecPerSite<T> *bv,
                        operator_struct<T> *op, const int amount, Level *l, const int blen) {

  int start_even, end_even, start_odd, end_odd,
      n = l->num_inner_lattice_sites, *neighbor = op->neighbor_table, start = 0,
      plus_dir_param = _FULL_SYSTEM, minus_dir_param = _FULL_SYSTEM;
  std::complex<T> *phi = bv->v;
  std::complex<T> *eta = bout->v;

  auto *prnT = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];
  auto *prnZ = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];
  auto *prnY = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];
  auto *prnX = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];

  auto *prpT = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];
  auto *prpZ = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];
  auto *prpY = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];
  auto *prpX = new std::complex<T>[(l->num_lattice_site_var)/2 * l->num_lattice_sites * blen];

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

  pi_minus(phi, prnT, prnZ, prnY, prnX, start*12*blen, n*12*blen);

  /// start communication in negative direction
  ghost_sendrecv_full(prnT, _T, -1, &(op->c), minus_dir_param, l, blen);
  ghost_sendrecv_full(prnZ, _Z, -1, &(op->c), minus_dir_param, l, blen);
  ghost_sendrecv_full(prnY, _Y, -1, &(op->c), minus_dir_param, l, blen);
  ghost_sendrecv_full(prnX, _X, -1, &(op->c), minus_dir_param, l, blen);

  /// project plus dir and multiply with U dagger
  pi_plus_kron_UH_psi(phi, op->D, prpT, prpZ, prpY, prpX, neighbor, start*12, n*12, blen);

  if (amount == _EVEN_SITES) {
    start = start_even, n = end_even;
  } else if (amount == _ODD_SITES) {
    start = start_odd, n = end_odd;
  }

  /// start communication in positive direction
  ghost_sendrecv_full(prpT, _T, +1, &(op->c), plus_dir_param, l, blen);
  ghost_sendrecv_full(prpZ, _Z, +1, &(op->c), plus_dir_param, l, blen);
  ghost_sendrecv_full(prpY, _Y, +1, &(op->c), plus_dir_param, l, blen);
  ghost_sendrecv_full(prpX, _X, +1, &(op->c), plus_dir_param, l, blen);

  /// wait for communication in negative direction
  ghost_wait_full(prnT, _T, -1, &(op->c), minus_dir_param, l, blen);
  ghost_wait_full(prnZ, _Z, -1, &(op->c), minus_dir_param, l, blen);
  ghost_wait_full(prnY, _Y, -1, &(op->c), minus_dir_param, l, blen);
  ghost_wait_full(prnX, _X, -1, &(op->c), minus_dir_param, l, blen);

  /// multiply with U and lift up minus dir
  apply_U_psi_plus_mu(eta, op->D, prnT, prnZ, prnY, prnX, neighbor, start*12, n*12, blen);

  /// wait for communication in positive direction
  ghost_wait_full(prpT, _T, +1, &(op->c), plus_dir_param, l, blen);
  ghost_wait_full(prpZ, _Z, +1, &(op->c), plus_dir_param, l, blen);
  ghost_wait_full(prpY, _Y, +1, &(op->c), plus_dir_param, l, blen);
  ghost_wait_full(prpX, _X, +1, &(op->c), plus_dir_param, l, blen);

  /// lift up plus dir
  lift_pi_plus_dir(eta, prpT, prpZ, prpY, prpX, start*12*blen, n*12*blen);

  delete[] prnT;
  delete[] prnZ;
  delete[] prnY;
  delete[] prnX;

  delete[] prpT;
  delete[] prpZ;
  delete[] prpY;
  delete[] prpX;
}

template<typename T>
void apply_schur_complement_block(BlockVecPerEntry<T>* vout, BlockVecPerEntry<T>* vin, operator_struct<T> *op,
                                  Level *l, const int blen) {

  /*********************************************************************************
   * Applies the Schur complement to a vector.
   *********************************************************************************/

  // start and end indices for vector functions
  int start_even = 0;
  int end_even = op->num_even_sites * l->num_lattice_site_var;
  int start_odd = op->num_even_sites * l->num_lattice_site_var;
  int end_odd = l->inner_vector_size;

  auto *tmp_odd = new std::complex<T>[blen * end_odd];
  auto *tmp_even = new std::complex<T>[blen * end_odd];
  std::complex<T>* out = vout->v;
  auto *vodd = new BlockVecPerEntry<T>(tmp_odd);
  auto *veven = new BlockVecPerEntry<T>(tmp_even);

  vector_double_define(tmp_odd, 0, start_odd * blen, end_odd * blen);
  vector_double_define(tmp_odd, 0, start_even * blen, end_even * blen);
  diag_ee_block(vout, vin, op, l, start_even, end_even, blen);
  hopping_term_block(vodd, vin, op, _ODD_SITES, l, blen);

  diag_oo_inv_block(veven, vodd, op, l, start_odd, end_odd, blen);
  hopping_term_block(vodd, veven, op, _EVEN_SITES, l, blen);
  vector_minus(out, out, tmp_odd, start_even * blen, end_even * blen);

  delete[] tmp_odd;
  delete[] tmp_even;
}

template<typename T>
void apply_schur_complement_block(BlockVecPerSite<T>* vout, BlockVecPerSite<T>* vin, operator_struct<T> *op,
                                  Level *l, const int blen) {

  /*********************************************************************************
   * Applies the Schur complement to a vector.
   *********************************************************************************/

  // start and end indices for vector functions
  int start_even = 0;
  int end_even = op->num_even_sites * l->num_lattice_site_var;
  int start_odd = op->num_even_sites * l->num_lattice_site_var;
  int end_odd = l->inner_vector_size;

  auto *tmp_odd = new std::complex<T>[blen * end_odd];
  auto *tmp_even = new std::complex<T>[blen * end_odd];
  std::complex<T>* out = vout->v;
  auto *vodd = new BlockVecPerSite<T>(tmp_odd, 12);
  auto *veven = new BlockVecPerSite<T>(tmp_even, 12);

  vector_double_define(tmp_odd, 0, start_odd * blen, end_odd * blen);
  vector_double_define(tmp_odd, 0, start_even * blen, end_even * blen);
  diag_ee_block(vout, vin, op, l, start_even, end_even, blen);
  hopping_term_block(vodd, vin, op, _ODD_SITES, l, blen);

  diag_oo_inv_block(veven, vodd, op, l, start_odd, end_odd, blen);
  hopping_term_block(vodd, veven, op, _EVEN_SITES, l, blen);
  vector_minus(out, out, tmp_odd, start_even * blen, end_even * blen);

  delete[] tmp_odd;
  delete[] tmp_even;
}

#endif // DDALPHAAMG_ODDEVEN_BLOCK_H
