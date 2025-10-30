#ifndef SCHWARZ_UTIL_HEADER
#define SCHWARZ_UTIL_HEADER

void Schwarz::blockBoundaryOperator(vector_double eta, vector_double phi, int direction, int k) {
  // k: number of current block
  int *bbl = block_boundary_length;

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
  double *Dplus = op.D_vectorized;
  double *Dminus = op.D_transformed_vectorized;

  for (int mu = 0; mu < 4; mu++) {
    if (direction == _POSITIVE) {
      boundary_plus_coupling_double((double *)eta, Dplus, (double *)phi, mu, bbl[2 * mu],
                                    bbl[2 * mu + 1], block[k].bt, nullptr);
      boundary_minus_coupling_double((double *)eta, Dminus, (double *)phi, mu, bbl[2 * mu + 1],
                                     bbl[2 * mu + 2], block[k].bt, nullptr);
    } else {
      boundary_nplus_coupling_double((double *)eta, Dplus, (double *)phi, mu, bbl[2 * mu],
                                     bbl[2 * mu + 1], block[k].bt, nullptr);
      boundary_nminus_coupling_double((double *)eta, Dminus, (double *)phi, mu, bbl[2 * mu + 1],
                                      bbl[2 * mu + 2], block[k].bt, nullptr);
    }
  }
#else
  int i, mu, index, neighbor_index;
  config_double D_pt, D = op.D;
  vector_double phi_pt, eta_pt;

#ifdef HAVE_TM1p1
  Global *global = level->global;
  if (global->n_flavours == 2) {
    complex_double buf1[24], *buf2 = buf1 + 12;
    mu = _T;
    // plus mu direction
    for (i = bbl[2 * mu]; i < bbl[2 * mu + 1]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * index + 9 * mu;
      phi_pt = phi + 24 * neighbor_index;
      eta_pt = eta + 24 * index;
      dprp_T_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvm_double(buf2, D_pt, buf1) : nmvm_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvm_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      (direction == _POSITIVE) ? mvm_double(buf2 + 6, D_pt, buf1 + 6)
                               : nmvm_double(buf2 + 6, D_pt, buf1 + 6);
      (direction == _POSITIVE) ? mvm_double(buf2 + 9, D_pt, buf1 + 9)
                               : nmvm_double(buf2 + 9, D_pt, buf1 + 9);
      dpbp_su3_T_double(buf2, eta_pt);
    }
    // minus mu direction
    for (i = bbl[2 * mu + 1]; i < bbl[2 * mu + 2]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * neighbor_index + 9 * mu;
      phi_pt = phi + 24 * neighbor_index;
      eta_pt = eta + 24 * index;
      dprn_T_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvmh_double(buf2, D_pt, buf1) : nmvmh_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 6, D_pt, buf1 + 6)
                               : nmvmh_double(buf2 + 6, D_pt, buf1 + 6);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 9, D_pt, buf1 + 9)
                               : nmvmh_double(buf2 + 9, D_pt, buf1 + 9);
      dpbn_su3_T_double(buf2, eta_pt);
    }

    mu = _Z;
    // plus mu direction
    for (i = bbl[2 * mu]; i < bbl[2 * mu + 1]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * index + 9 * mu;
      phi_pt = phi + 24 * neighbor_index;
      eta_pt = eta + 24 * index;
      dprp_Z_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvm_double(buf2, D_pt, buf1) : nmvm_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvm_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      (direction == _POSITIVE) ? mvm_double(buf2 + 6, D_pt, buf1 + 6)
                               : nmvm_double(buf2 + 6, D_pt, buf1 + 6);
      (direction == _POSITIVE) ? mvm_double(buf2 + 9, D_pt, buf1 + 9)
                               : nmvm_double(buf2 + 9, D_pt, buf1 + 9);
      dpbp_su3_Z_double(buf2, eta_pt);
    }
    // minus mu direction
    for (i = bbl[2 * mu + 1]; i < bbl[2 * mu + 2]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * neighbor_index + 9 * mu;
      phi_pt = phi + 24 * neighbor_index;
      eta_pt = eta + 24 * index;
      dprn_Z_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvmh_double(buf2, D_pt, buf1) : nmvmh_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 6, D_pt, buf1 + 6)
                               : nmvmh_double(buf2 + 6, D_pt, buf1 + 6);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 9, D_pt, buf1 + 9)
                               : nmvmh_double(buf2 + 9, D_pt, buf1 + 9);
      dpbn_su3_Z_double(buf2, eta_pt);
    }

    mu = _Y;
    // plus mu direction
    for (i = bbl[2 * mu]; i < bbl[2 * mu + 1]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * index + 9 * mu;
      phi_pt = phi + 24 * neighbor_index;
      eta_pt = eta + 24 * index;
      dprp_Y_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvm_double(buf2, D_pt, buf1) : nmvm_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvm_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      (direction == _POSITIVE) ? mvm_double(buf2 + 6, D_pt, buf1 + 6)
                               : nmvm_double(buf2 + 6, D_pt, buf1 + 6);
      (direction == _POSITIVE) ? mvm_double(buf2 + 9, D_pt, buf1 + 9)
                               : nmvm_double(buf2 + 9, D_pt, buf1 + 9);
      dpbp_su3_Y_double(buf2, eta_pt);
    }
    // minus mu direction
    for (i = bbl[2 * mu + 1]; i < bbl[2 * mu + 2]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * neighbor_index + 9 * mu;
      phi_pt = phi + 24 * neighbor_index;
      eta_pt = eta + 24 * index;
      dprn_Y_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvmh_double(buf2, D_pt, buf1) : nmvmh_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 6, D_pt, buf1 + 6)
                               : nmvmh_double(buf2 + 6, D_pt, buf1 + 6);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 9, D_pt, buf1 + 9)
                               : nmvmh_double(buf2 + 9, D_pt, buf1 + 9);
      dpbn_su3_Y_double(buf2, eta_pt);
    }

    mu = _X;
    // plus mu direction
    for (i = bbl[2 * mu]; i < bbl[2 * mu + 1]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * index + 9 * mu;
      phi_pt = phi + 24 * neighbor_index;
      eta_pt = eta + 24 * index;
      dprp_X_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvm_double(buf2, D_pt, buf1) : nmvm_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvm_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      (direction == _POSITIVE) ? mvm_double(buf2 + 6, D_pt, buf1 + 6)
                               : nmvm_double(buf2 + 6, D_pt, buf1 + 6);
      (direction == _POSITIVE) ? mvm_double(buf2 + 9, D_pt, buf1 + 9)
                               : nmvm_double(buf2 + 9, D_pt, buf1 + 9);
      dpbp_su3_X_double(buf2, eta_pt);
    }
    // minus mu direction
    for (i = bbl[2 * mu + 1]; i < bbl[2 * mu + 2]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * neighbor_index + 9 * mu;
      phi_pt = phi + 24 * neighbor_index;
      eta_pt = eta + 24 * index;
      dprn_X_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvmh_double(buf2, D_pt, buf1) : nmvmh_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 6, D_pt, buf1 + 6)
                               : nmvmh_double(buf2 + 6, D_pt, buf1 + 6);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 9, D_pt, buf1 + 9)
                               : nmvmh_double(buf2 + 9, D_pt, buf1 + 9);
      dpbn_su3_X_double(buf2, eta_pt);
    }
  } else {
#endif
    complex_double buf1[12], *buf2 = buf1 + 6;
    mu = _T;
    // plus mu direction
    for (i = bbl[2 * mu]; i < bbl[2 * mu + 1]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * index + 9 * mu;
      phi_pt = phi + 12 * neighbor_index;
      eta_pt = eta + 12 * index;
      prp_T_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvm_double(buf2, D_pt, buf1) : nmvm_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvm_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      pbp_su3_T_double(buf2, eta_pt);
    }
    // minus mu direction
    for (i = bbl[2 * mu + 1]; i < bbl[2 * mu + 2]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * neighbor_index + 9 * mu;
      phi_pt = phi + 12 * neighbor_index;
      eta_pt = eta + 12 * index;
      prn_T_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvmh_double(buf2, D_pt, buf1) : nmvmh_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      pbn_su3_T_double(buf2, eta_pt);
    }

    mu = _Z;
    // plus mu direction
    for (i = bbl[2 * mu]; i < bbl[2 * mu + 1]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * index + 9 * mu;
      phi_pt = phi + 12 * neighbor_index;
      eta_pt = eta + 12 * index;
      prp_Z_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvm_double(buf2, D_pt, buf1) : nmvm_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvm_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      pbp_su3_Z_double(buf2, eta_pt);
    }
    // minus mu direction
    for (i = bbl[2 * mu + 1]; i < bbl[2 * mu + 2]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * neighbor_index + 9 * mu;
      phi_pt = phi + 12 * neighbor_index;
      eta_pt = eta + 12 * index;
      prn_Z_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvmh_double(buf2, D_pt, buf1) : nmvmh_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      pbn_su3_Z_double(buf2, eta_pt);
    }

    mu = _Y;
    // plus mu direction
    for (i = bbl[2 * mu]; i < bbl[2 * mu + 1]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * index + 9 * mu;
      phi_pt = phi + 12 * neighbor_index;
      eta_pt = eta + 12 * index;
      prp_Y_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvm_double(buf2, D_pt, buf1) : nmvm_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvm_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      pbp_su3_Y_double(buf2, eta_pt);
    }
    // minus mu direction
    for (i = bbl[2 * mu + 1]; i < bbl[2 * mu + 2]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * neighbor_index + 9 * mu;
      phi_pt = phi + 12 * neighbor_index;
      eta_pt = eta + 12 * index;
      prn_Y_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvmh_double(buf2, D_pt, buf1) : nmvmh_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      pbn_su3_Y_double(buf2, eta_pt);
    }

    mu = _X;
    // plus mu direction
    for (i = bbl[2 * mu]; i < bbl[2 * mu + 1]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * index + 9 * mu;
      phi_pt = phi + 12 * neighbor_index;
      eta_pt = eta + 12 * index;
      prp_X_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvm_double(buf2, D_pt, buf1) : nmvm_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvm_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvm_double(buf2 + 3, D_pt, buf1 + 3);
      pbp_su3_X_double(buf2, eta_pt);
    }
    // minus mu direction
    for (i = bbl[2 * mu + 1]; i < bbl[2 * mu + 2]; i += 2) {
      index = block[k].bt[i];
      neighbor_index = block[k].bt[i + 1];
      D_pt = D + 36 * neighbor_index + 9 * mu;
      phi_pt = phi + 12 * neighbor_index;
      eta_pt = eta + 12 * index;
      prn_X_double(buf1, phi_pt);
      (direction == _POSITIVE) ? mvmh_double(buf2, D_pt, buf1) : nmvmh_double(buf2, D_pt, buf1);
      (direction == _POSITIVE) ? mvmh_double(buf2 + 3, D_pt, buf1 + 3)
                               : nmvmh_double(buf2 + 3, D_pt, buf1 + 3);
      pbn_su3_X_double(buf2, eta_pt);
    }
#ifdef HAVE_TM1p1
  }
#endif
#endif
}

void Schwarz::coarseBlockBoundaryOperator(vector_double eta, vector_double phi, int direction,
                                          int k) {
  // k: number of current block
  int *bbl = block_boundary_length, n = level->num_lattice_site_var;

#ifdef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_double
  int column_offset = 2 * SIMD_LENGTH_double *
                      ((level->num_parent_eig_vect + SIMD_LENGTH_double - 1) / SIMD_LENGTH_double);
  int vectorized_link_offset = 4 * level->num_parent_eig_vect * column_offset;

  for (int mu = 0; mu < 4; mu++) {
    OPERATOR_TYPE_double *Dplus = op.D_vectorized + mu * vectorized_link_offset;
    OPERATOR_TYPE_double *Dminus = op.D_transformed_vectorized + mu * vectorized_link_offset;
    // plus mu direction
    for (int i = bbl[2 * mu]; i < bbl[2 * mu + 1]; i += 2) {
      int index = block[k].bt[i];
      int neighbor_index = block[k].bt[i + 1];
      vector_double phi_pt = phi + n * neighbor_index;
      vector_double eta_pt = eta + n * index;
      (direction == _POSITIVE)
          ? coarse_hopp_double_vectorized(eta_pt, phi_pt,
                                          Dplus + 4 * vectorized_link_offset * index, level)
          : coarse_n_hopp_double_vectorized(eta_pt, phi_pt,
                                            Dplus + 4 * vectorized_link_offset * index, level);
    }
    // minus mu direction
    for (int i = bbl[2 * mu + 1]; i < bbl[2 * mu + 2]; i += 2) {
      int index = block[k].bt[i];
      int neighbor_index = block[k].bt[i + 1];
      vector_double phi_pt = phi + n * neighbor_index;
      vector_double eta_pt = eta + n * index;
      (direction == _POSITIVE)
          ? coarse_hopp_double_vectorized(
                eta_pt, phi_pt, Dminus + 4 * vectorized_link_offset * neighbor_index, level)
          : coarse_n_hopp_double_vectorized(
                eta_pt, phi_pt, Dminus + 4 * vectorized_link_offset * neighbor_index, level);
    }
  }
#else
  config_double D = op.D;
  int link_size = SQUARE(2 * level->num_parent_eig_vect), site_size = 4 * link_size;

  for (int mu = 0; mu < 4; mu++) {
    // plus mu direction
    for (int i = bbl[2 * mu]; i < bbl[2 * mu + 1]; i += 2) {
      int index = block[k].bt[i];
      int neighbor_index = block[k].bt[i + 1];
      vector_double phi_pt = phi + n * neighbor_index;
      vector_double eta_pt = eta + n * index;
      config_double D_pt = D + site_size * index + link_size * mu;
      (direction == _POSITIVE) ? coarse_hopp_double(eta_pt, phi_pt, D_pt, level)
                               : coarse_n_hopp_double(eta_pt, phi_pt, D_pt, level);
    }
    // minus mu direction
    for (int i = bbl[2 * mu + 1]; i < bbl[2 * mu + 2]; i += 2) {
      int index = block[k].bt[i];
      int neighbor_index = block[k].bt[i + 1];
      vector_double phi_pt = phi + n * neighbor_index;
      vector_double eta_pt = eta + n * index;
      config_double D_pt = D + site_size * neighbor_index + link_size * mu;
      (direction == _POSITIVE) ? coarse_daggered_hopp_double(eta_pt, phi_pt, D_pt, level)
                               : coarse_n_daggered_hopp_double(eta_pt, phi_pt, D_pt, level);
    }
  }
#endif
}

void Schwarz::localMinres(vector_double phi, vector_double eta, vector_double latest_iter,
                          int minresStart) {

  /*********************************************************************************
   * Minimal Residual iteration solver used to solve the block systems
   *     blockD phi = eta
   * within the Schwarz method, phi contains an initial guess and its updated version
   * is returned after the block solve has been performed.
   * eta is overwritten by the block residual r.
   * To calculate the missing contributions to r on the current Schwarz block
   * coming from outside of the block, an update "phi_new - phi_old" is returned in
   * latest_iter -> cheaper residual update in the Schwarz method
   *********************************************************************************/

  int minresEnd = (level->odd_even && level->type == _FINE)
                      ? (minresStart + level->num_lattice_site_var * number_of_even_block_sites)
                      : (minresStart + block_vector_size);
  vector_double Dr = local_minres_buffer[0];
  vector_double r = local_minres_buffer[1];
  vector_double lphi = local_minres_buffer[2];
  complex_double alpha;
  void (*minresOp)(vector_double, vector_double, int, Schwarz *, Level *) =
      level->type == _FINE
          ? (level->odd_even ? apply_block_schur_complement_double : block_d_plus_clover_double)
          : coarse_block_operator_double;

  vector_double_copy(r, eta, minresStart, minresEnd, level);
  vector_double_define(lphi, 0, minresStart, minresEnd);

  for (int i = 0; i < level->block_iter; i++) {
    minresOp(Dr, r, minresStart, this, level);                             // Dr = blockD*r
    alpha = local_xy_over_xx_double(Dr, r, minresStart, minresEnd, level); // alpha = <Dr,r>/<Dr,Dr>
    vector_double_saxpy(lphi, lphi, r, alpha, minresStart, minresEnd, level); // phi += alpha * r
    vector_double_saxpy(r, r, Dr, -alpha, minresStart, minresEnd, level);     // r -= alpha * Dr
  }

  if (latest_iter != nullptr)
    vector_double_copy(latest_iter, lphi, minresStart, minresEnd, level);
  if (phi != nullptr)
    vector_double_plus(phi, phi, lphi, minresStart, minresEnd, level);
  vector_double_copy(eta, r, minresStart, minresEnd, level);
}

void Schwarz::blockSolveOddeven(vector_double phi, vector_double r, vector_double latest_iter,
                                int minresStart) {
  int minresEnd = minresStart + block_vector_size;

  // odd to even
  vector_double_copy(odd_even_buffer[3], r, minresStart, minresEnd, level);
  block_diag_oo_inv_double(odd_even_buffer[2], odd_even_buffer[3], minresStart, this, level);
  block_n_hopping_term_double(odd_even_buffer[3], odd_even_buffer[2], minresStart, _EVEN_SITES,
                              this, level);

  Schwarz::localMinres(nullptr, odd_even_buffer[3], odd_even_buffer[2], minresStart);

  // even to odd
  block_n_hopping_term_double(odd_even_buffer[3], odd_even_buffer[2], minresStart, _ODD_SITES, this,
                              level);
  block_diag_oo_inv_double(odd_even_buffer[2], odd_even_buffer[3], minresStart, this, level);

  // update phi, latest_iter
  vector_double_copy(latest_iter, odd_even_buffer[2], minresStart, minresEnd, level);
  vector_double_plus(phi, phi, odd_even_buffer[2], minresStart, minresEnd, level);
  // update r
  vector_double_copy(r, odd_even_buffer[3], minresStart,
                     minresStart + level->num_lattice_site_var * number_of_even_block_sites, level);
  vector_double_define(r, 0, minresStart + level->num_lattice_site_var * number_of_even_block_sites,
                       minresEnd);
}

inline int Schwarz::connectLink(int t, int z, int y, int x, int mu, int dir, int *dt, int *it) {

  int coord[4];
  coord[_T] = t;
  coord[_Z] = z;
  coord[_Y] = y;
  coord[_X] = x;
  coord[mu] += dir;
  if (level->global_splitting[mu] > 1) {
    return site_mod_index(coord[_T], coord[_Z], coord[_Y], coord[_X], dt, it);
  } else {
    coord[mu] = (coord[mu] + level->local_lattice[mu]) % level->local_lattice[mu];
    return site_index(coord[_T], coord[_Z], coord[_Y], coord[_X], dt, it);
  }

  return 0; //To please the compiler gods
}

void Schwarz::defineLayout() {

  int a0, b0, c0, d0, a1, b1, c1, d1, block_split[4], block_size[4], agg_split[4], i, j, k, mu,
      index, x, y, z, t, ls[4], le[4], l_st[4], l_en[4],
      *dt = op.table_dim, *dt_mod = op.table_mod_dim, *it = op.index_table, *count[4];

  const int sigma[16] = {0, 1, 3, 2, 6, 4, 5, 7, 15, 14, 12, 13, 9, 11, 10, 8};
  const int color_to_comm[16][2] = {{_T, -1}, {_X, +1}, {_Y, +1}, {_X, -1}, {_Z, +1}, {_Y, -1},
                                    {_X, +1}, {_Y, +1}, {_T, +1}, {_X, -1}, {_Y, -1}, {_X, +1},
                                    {_Z, -1}, {_Y, +1}, {_X, -1}, {_Y, -1}};
  int color_counter[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  Global *global = level->global;

  // Define coloring
  switch (global->method) {
  case _ADDITIVE:
    number_of_colors = 1;
    break;
  case _TWO_COLOR:
    number_of_colors = 2;
    break;
  case _SIXTEEN_COLOR: {
    // check if 16 color Schwarz is possible
    for (mu = 0; mu < 4; mu++) {
      if ((level->local_lattice[mu] / level->block_lattice[mu]) % 2 == 1) {
        number_of_colors = 2;
        printf0("depth: %d, switching to red black schwarz as smoother\n", level->depth);
        break;
      }
    }
    number_of_colors = 16;
    break;
  }
  default:
    number_of_colors = 0;
    break;
  }

  number_of_block_sites = 1;
  block_odd_even_offset = 0;
  number_of_aggregates = 1;

  if (global->method == _TWO_COLOR) {
    for (i = 0; i < 8; i++)
      block_list_length[i] = 0;
  }

  for (mu = 0; mu < 4; mu++) {
    number_of_block_sites *= level->block_lattice[mu];
    block_odd_even_offset += ((level->local_lattice[mu] / level->block_lattice[mu]) *
                              (global->my_coords[mu] / level->comm_offset[mu])) %
                             2;
    ls[mu] = 0;
    le[mu] = ls[mu] + level->local_lattice[mu];
    // set the dimension according to its initialization in Schwarz::allocateMemory()
    if(level->type != _COARSEST)
      dt[mu] = level->local_lattice[mu] + 2;
    else
      dt[mu] = level->local_lattice[mu] + 1;
    dt_mod[mu] = level->local_lattice[mu] + 2;
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
    agg_split[mu] = level->local_lattice[mu] / level->coarsening[mu];
    block_split[mu] = level->coarsening[mu] / level->block_lattice[mu];
    block_size[mu] = level->block_lattice[mu];
    number_of_aggregates *= agg_split[mu];
  }
  block_odd_even_offset = block_odd_even_offset % 2;
  block_vector_size = number_of_block_sites * level->num_lattice_site_var;

  i = 0;
  j = 0;
  // inner hyper cuboid
  count[_T] = &d1;
  count[_Z] = &c1;
  count[_Y] = &b1;
  count[_X] = &a1;
  for (d0 = 0; d0 < agg_split[_T]; d0++)
    for (c0 = 0; c0 < agg_split[_Z]; c0++)
      for (b0 = 0; b0 < agg_split[_Y]; b0++)
        for (a0 = 0; a0 < agg_split[_X]; a0++) {

          for (d1 = d0 * block_split[_T]; d1 < (d0 + 1) * block_split[_T]; d1++)
            for (c1 = c0 * block_split[_Z]; c1 < (c0 + 1) * block_split[_Z]; c1++)
              for (b1 = b0 * block_split[_Y]; b1 < (b0 + 1) * block_split[_Y]; b1++)
                for (a1 = a0 * block_split[_X]; a1 < (a0 + 1) * block_split[_X]; a1++) {

                  block[j].start = i;
                  block[j].no_comm = 1;
                  if (number_of_colors == 1) {
                    block[j].color = 0;
                  } else if (number_of_colors == 2) {
                    block[j].color = (d1 + c1 + b1 + a1 + block_odd_even_offset) % 2;
                  } else if (number_of_colors == 16) {
                    for (k = 0; k < 16; k++)
                      if (sigma[k] == 8 * (d1 % 2) + 4 * (c1 % 2) + 2 * (b1 % 2) + 1 * (a1 % 2)) {
                        block[j].color = k;
                        block_list[k][color_counter[k]] = j;
                        color_counter[k]++;
                        break;
                      }
                  }

                  if (number_of_colors == 1 || number_of_colors == 2) {
                    for (mu = 0; mu < 4; mu++) {
                      if (((*count[mu]) == 0) || ((*count[mu] + 1) == le[mu] / block_size[mu]))
                        block[j].no_comm = 0;
                    }

                    if (number_of_colors == 2) {
                      // calculate boundary correspondence of the block
                      int count_plus = 0, count_minus = 0, count_inner = 0, index;
                      for (mu = 0; mu < 4; mu++) {
                        if ((*count[mu]) == 0)
                          count_minus++;
                        if ((*count[mu] + 1) == le[mu] / block_size[mu])
                          count_plus++;
                        if ((*count[mu]) != 0 && (*count[mu] + 1) != le[mu] / block_size[mu])
                          count_inner++;
                      }

                      if (count_inner == 4) {
                        index = 4 * block[j].color;
                      } else if (count_minus == 0) {
                        if (block[j].color == 0)
                          index = 1;
                        else
                          index = 7;
                      } else if (count_plus == 0) {
                        if (block[j].color == 0)
                          index = 3;
                        else
                          index = 5;
                      } else {
                        index = 2 + 4 * block[j].color;
                      }

                      block_list[index][block_list_length[index]] = j;
                      block_list_length[index]++;
                    }

                  } else if (number_of_colors == 16) {
                    k = block[j].color;
                    if (k == 0) {
                      for (mu = 0; mu < 4; mu++) {
                        if ((*count[mu]) == 0)
                          block[j].no_comm = 0;
                      }
                    } else {
                      mu = color_to_comm[k][0];
                      if ((color_to_comm[k][1] == +1 &&
                           (*count[mu] + 1) == le[mu] / block_size[mu]) ||
                          (color_to_comm[k][1] == -1 && (*count[mu]) == 0))
                        block[j].no_comm = 0;
                    }
                  }

                  j++;

                  // set up index table
                  if (level->type == _FINE && level->odd_even) {
                    // odd even on the blocks
                    // even sites
                    for (t = d1 * block_size[_T]; t < (d1 + 1) * block_size[_T]; t++)
                      for (z = c1 * block_size[_Z]; z < (c1 + 1) * block_size[_Z]; z++)
                        for (y = b1 * block_size[_Y]; y < (b1 + 1) * block_size[_Y]; y++)
                          for (x = a1 * block_size[_X]; x < (a1 + 1) * block_size[_X]; x++) {
                            if (((t - d1 * block_size[_T]) + (z - c1 * block_size[_Z]) +
                                 (y - b1 * block_size[_Y]) + (x - a1 * block_size[_X])) %
                                    2 ==
                                0) {
                              index = lex_index(t, z, y, x, dt);
                              it[index] = i;
                              i++;
                            }
                          }
                    // odd sites
                    for (t = d1 * block_size[_T]; t < (d1 + 1) * block_size[_T]; t++)
                      for (z = c1 * block_size[_Z]; z < (c1 + 1) * block_size[_Z]; z++)
                        for (y = b1 * block_size[_Y]; y < (b1 + 1) * block_size[_Y]; y++)
                          for (x = a1 * block_size[_X]; x < (a1 + 1) * block_size[_X]; x++) {
                            if (((t - d1 * block_size[_T]) + (z - c1 * block_size[_Z]) +
                                 (y - b1 * block_size[_Y]) + (x - a1 * block_size[_X])) %
                                    2 ==
                                1) {
                              index = lex_index(t, z, y, x, dt);
                              it[index] = i;
                              i++;
                            }
                          }
                  } else {
                    // no odd even
                    for (t = d1 * block_size[_T]; t < (d1 + 1) * block_size[_T]; t++)
                      for (z = c1 * block_size[_Z]; z < (c1 + 1) * block_size[_Z]; z++)
                        for (y = b1 * block_size[_Y]; y < (b1 + 1) * block_size[_Y]; y++)
                          for (x = a1 * block_size[_X]; x < (a1 + 1) * block_size[_X]; x++) {
                            index = lex_index(t, z, y, x, dt);
                            it[index] = i;
                            i++;
                          }
                  }
                }
        }

  // boundaries
  for (mu = 0; mu < 4; mu++) {
    l_st[mu] = le[mu];
    l_en[mu] = le[mu] + 1;
    for (t = l_st[_T]; t < l_en[_T]; t++)
      for (z = l_st[_Z]; z < l_en[_Z]; z++)
        for (y = l_st[_Y]; y < l_en[_Y]; y++)
          for (x = l_st[_X]; x < l_en[_X]; x++) {
            index = lex_index(t, z, y, x, dt);
            it[index] = i;
            i++;
          }
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }

  // define negative boundaries
  for (mu = 0; mu < 4; mu++) {
    l_st[mu] = ls[mu] - 1;
    l_en[mu] = ls[mu];
    for (t = l_st[_T]; t < l_en[_T]; t++)
      for (z = l_st[_Z]; z < l_en[_Z]; z++)
        for (y = l_st[_Y]; y < l_en[_Y]; y++)
          for (x = l_st[_X]; x < l_en[_X]; x++) {
            index = lex_mod_index(t, z, y, x, dt);
            it[index] = i;
            i++;
          }
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }

  i = 0;
  j = 0;
  // block boundary table
  for (d0 = 0; d0 < agg_split[_T]; d0++)
    for (c0 = 0; c0 < agg_split[_Z]; c0++)
      for (b0 = 0; b0 < agg_split[_Y]; b0++)
        for (a0 = 0; a0 < agg_split[_X]; a0++)

          for (d1 = d0 * block_split[_T]; d1 < (d0 + 1) * block_split[_T]; d1++)
            for (c1 = c0 * block_split[_Z]; c1 < (c0 + 1) * block_split[_Z]; c1++)
              for (b1 = b0 * block_split[_Y]; b1 < (b0 + 1) * block_split[_Y]; b1++)
                for (a1 = a0 * block_split[_X]; a1 < (a0 + 1) * block_split[_X]; a1++) {
                  // for all blocks
                  i = 0;
                  int block_start[4], block_end[4], tmp;

                  block_start[_T] = d1 * block_size[_T];
                  block_start[_Z] = c1 * block_size[_Z];
                  block_start[_Y] = b1 * block_size[_Y];
                  block_start[_X] = a1 * block_size[_X];
                  block_end[_T] = (d1 + 1) * block_size[_T];
                  block_end[_Z] = (c1 + 1) * block_size[_Z];
                  block_end[_Y] = (b1 + 1) * block_size[_Y];
                  block_end[_X] = (a1 + 1) * block_size[_X];

                  for (mu = 0; mu < 4; mu++) {
                    tmp = block_start[mu];

                    // minus dir
                    block_start[mu] = block_end[mu] - 1;
                    for (t = block_start[_T]; t < block_end[_T]; t++)
                      for (z = block_start[_Z]; z < block_end[_Z]; z++)
                        for (y = block_start[_Y]; y < block_end[_Y]; y++)
                          for (x = block_start[_X]; x < block_end[_X]; x++) {
                            block[j].bt[i] = site_index(t, z, y, x, dt, it);
                            i++;
                            block[j].bt[i] = connectLink(t, z, y, x, mu, +1, dt, it);
                            i++;
                          }

                    block_start[mu] = tmp;
                    tmp = block_end[mu];

                    // plus dir
                    block_end[mu] = block_start[mu] + 1;
                    for (t = block_start[_T]; t < block_end[_T]; t++)
                      for (z = block_start[_Z]; z < block_end[_Z]; z++)
                        for (y = block_start[_Y]; y < block_end[_Y]; y++)
                          for (x = block_start[_X]; x < block_end[_X]; x++) {
                            block[j].bt[i] = site_index(t, z, y, x, dt, it);
                            i++;
                            block[j].bt[i] = connectLink(t, z, y, x, mu, -1, dt, it);
                            i++;
                          }
                    block_end[mu] = tmp;
                  }
                  j++;
                }

  // index table for block dirac operator
  if (level->type == _FINE && level->odd_even) {
    count[_T] = &t;
    count[_Z] = &z;
    count[_Y] = &y;
    count[_X] = &x;
    i = 0;
    j = 0;
    for (t = 0; t < block_size[_T]; t++)
      for (z = 0; z < block_size[_Z]; z++)
        for (y = 0; y < block_size[_Y]; y++)
          for (x = 0; x < block_size[_X]; x++) {
            if ((t + z + y + x) % 2 == 0) {
              i++;
            } else {
              j++;
            }
          }
    number_of_even_block_sites = i;
    number_of_odd_block_sites = j;

    for (mu = 0; mu < 4; mu++) {
      // even sites, plus dir ( = odd sites, minus dir )
      i = 0;
      j = 0;
      for (t = 0; t < block_size[_T]; t++)
        for (z = 0; z < block_size[_Z]; z++)
          for (y = 0; y < block_size[_Y]; y++)
            for (x = 0; x < block_size[_X]; x++) {
              if ((t + z + y + x) % 2 == 0) {
                if (*(count[mu]) < block_size[mu] - 1) {
                  odd_even_index_table[mu][j] = i;
                  j++;
                }
                i++;
              }
            }
      direction_length_even[mu] = j;
      // odd sites, plus dir ( = even sites, minus dir )
      j = 0;
      for (t = 0; t < block_size[_T]; t++)
        for (z = 0; z < block_size[_Z]; z++)
          for (y = 0; y < block_size[_Y]; y++)
            for (x = 0; x < block_size[_X]; x++) {
              if ((t + z + y + x) % 2 == 1) {
                if (*(count[mu]) < block_size[mu] - 1) {
                  odd_even_index_table[mu][direction_length_even[mu] + j] = i;
                  j++;
                }
                i++;
              }
            }
      direction_length_odd[mu] = j;
    }
  }

  count[_T] = &t;
  count[_Z] = &z;
  count[_Y] = &y;
  count[_X] = &x;
  for (mu = 0; mu < 4; mu++) {
    j = 0;
    for (t = 0; t < block_size[_T]; t++)
      for (z = 0; z < block_size[_Z]; z++)
        for (y = 0; y < block_size[_Y]; y++)
          for (x = 0; x < block_size[_X]; x++) {
            if (*(count[mu]) < block_size[mu] - 1) {
              index_table[mu][j] = site_index(t, z, y, x, dt, it);
              j++;
            }
          }
  }

  // define neighbor table (for the application of the entire operator),
  // negative inner boundary table (for communication),
  // translation table (for translation to lexicographical site ordnering)
  define_nt_bt_tt(op.neighbor_table, op.backward_neighbor_table, op.c.boundary_table,
                  op.translation_table, it, dt, level);
}

void Schwarz::boundaryUpdate() {

  /*********************************************************************************
   * Updates the current level hopping term in "op.D" on the process boundaries
   * in all negative directions. This is necessary for enabling Schwarz to perform
   * local block residual updates on demand.
   *********************************************************************************/

  int i, t, z, y, x, mu, nu, index, *it = op.index_table, *dt = op.table_dim, ls[4], le[4],
                                    buf_length[4], link_size;
  vector_double buf[4] = {nullptr, nullptr, nullptr, nullptr},
                rbuf[4] = {nullptr, nullptr, nullptr, nullptr};
  config_double D = op.D;
  Global *global = level->global;

  for (mu = 0; mu < 4; mu++) {
    ls[mu] = 0;
    le[mu] = level->local_lattice[mu];
    buf_length[mu] = 0;
  }

  if (level->type == _FINE)
    link_size = 4 * 9;
  else
    link_size = 4 * SQUARE(level->num_lattice_site_var);

  // allocate buffers
  for (mu = 0; mu < 4; mu++) {
    if (level->global_splitting[mu] > 1) {
      buf_length[mu] = link_size;
      for (nu = 0; nu < 4; nu++) {
        if (nu != mu)
          buf_length[mu] *= le[nu];
      }
      buf[mu] = new complex_double[buf_length[mu]];
      rbuf[mu] = new complex_double[buf_length[mu]];
    }
  }

  // post recv for desired directions
  for (mu = 0; mu < 4; mu++) {
    if (level->global_splitting[mu] > 1) {
      MPI_Irecv(rbuf[mu], buf_length[mu], MPI_COMPLEX_double, level->neighbor_rank[2 * mu + 1],
                2 * mu + 1, global->comm_cart, &(op.c.rreqs[2 * mu + 1]));
    }
  }

  // buffer data for send and send it
  for (mu = 0; mu < 4; mu++) {
    if (level->global_splitting[mu] > 1) {
      ls[mu] = level->local_lattice[mu] - 1;
      i = 0;
      for (t = ls[_T]; t < le[_T]; t++)
        for (z = ls[_Z]; z < le[_Z]; z++)
          for (y = ls[_Y]; y < le[_Y]; y++)
            for (x = ls[_X]; x < le[_X]; x++) {
              index = site_index(t, z, y, x, dt, it);
              vector_double_copy(buf[mu] + i * link_size, D + index * link_size, 0, link_size,
                                 level);
              i++;
            }
      MPI_Isend(buf[mu], buf_length[mu], MPI_COMPLEX_double, level->neighbor_rank[2 * mu],
                2 * mu + 1, global->comm_cart, &(op.c.sreqs[2 * mu + 1]));
      ls[mu] = 0;
    }
  }

  // store links in desired ordering after recv
  for (mu = 0; mu < 4; mu++) {
    if (level->global_splitting[mu] > 1) {
      MPI_Wait(&(op.c.rreqs[2 * mu + 1]), MPI_STATUS_IGNORE);
      ls[mu] = -1;
      le[mu] = 0;
      i = 0;
      for (t = ls[_T]; t < le[_T]; t++)
        for (z = ls[_Z]; z < le[_Z]; z++)
          for (y = ls[_Y]; y < le[_Y]; y++)
            for (x = ls[_X]; x < le[_X]; x++) {
              index = site_mod_index(t, z, y, x, dt, it);
              vector_double_copy(D + index * link_size, rbuf[mu] + i * link_size, 0, link_size,
                                 level);
              i++;
            }
      ls[mu] = 0;
      le[mu] = level->local_lattice[mu];
    }
  }

  // free buffers
  for (mu = 0; mu < 4; mu++) {
    if (level->global_splitting[mu] > 1) {
      MPI_Wait(&(op.c.sreqs[2 * mu + 1]), MPI_STATUS_IGNORE);
      delete[] buf[mu];
      delete[] rbuf[mu];
    }
  }
}

void Schwarz::setupOddEvenOperator() {

  int mu, i, d0, c0, b0, a0, d1, c1, b1, a1, t, z, y, x, agg_split[4], block_split[4],
      block_size[4];
  int n1 = number_of_even_block_sites;
#ifdef HAVE_TM
  config_double tm_term_pt = op.tm_term;
#endif

  for (mu = 0; mu < 4; mu++) {
    agg_split[mu] = level->local_lattice[mu] / level->coarsening[mu];
    block_split[mu] = level->coarsening[mu] / level->block_lattice[mu];
    block_size[mu] = level->block_lattice[mu];
  }

  if (level->csw) {
#ifndef OPTIMIZED_SELF_COUPLING_double
    config_double clover_pt = op.clover, clover_oo_inv_pt = op.clover_oo_inv;
    complex_double buffer[42];
    int cs = 42;
#else
    double *clover_pt = op.clover_vectorized, *clover_oo_inv_pt = op.clover_oo_inv_vectorized;
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
  config_double clover_oo_inv_pt = op.clover_doublet_oo_inv, clover_pt = op.clover;
  int cs = global->csw ? 42 : 12;
#else
  double *clover_pt = global->csw ? op.clover_doublet_vectorized : (double *)op.clover,
         *clover_oo_inv_pt = op.clover_doublet_oo_inv_vectorized;
  int cs = level->csw ? 288 : 24;
#endif
  config_double eps_term_pt = op.epsbar_term;
#ifdef HAVE_TM
  tm_term_pt = op.tm_term;
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
                            if (level->csw) {
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

void Schwarz::setupFineOperator(operator_struct<double> *op_in) {

  /*********************************************************************************
   * Copies the Dirac operator and the clover term from op_in into the schwarz_double_struct
   * struct.
   * - operator_struct<double> *op_in: Input operator.
   *********************************************************************************/

  int n = level->num_inner_lattice_sites, *tt = op.translation_table;
  config_double D_out_pt, clover_out_pt, odd_proj_out_pt;
  config_double D_in_pt = op_in->D, clover_in_pt = op_in->clover, odd_proj_in_pt = op_in->odd_proj;

  Global *global = level->global;

  op.m0 = op_in->m0;
  for (int i = 0; i < n; i++) {
    D_out_pt = op.D + 36 * tt[i];
    FOR36(*D_out_pt = (complex_double)*D_in_pt; D_out_pt++; D_in_pt++;);
  }

  if (global->csw != 0) {
    for (int i = 0; i < n; i++) {
      clover_out_pt = op.clover + 42 * tt[i];
      FOR42(*clover_out_pt = (complex_double)*clover_in_pt; clover_out_pt++; clover_in_pt++;);
    }
  } else {
    for (int i = 0; i < n; i++) {
      clover_out_pt = op.clover + 12 * tt[i];
      FOR12(*clover_out_pt = (complex_double)*clover_in_pt; clover_out_pt++; clover_in_pt++;);
    }
  }

  for (int i = 0; i < n; i++) {
    odd_proj_out_pt = op.odd_proj + 12 * tt[i];
    FOR12(*odd_proj_out_pt = (complex_double)*odd_proj_in_pt; odd_proj_out_pt++; odd_proj_in_pt++;);
  }

#ifdef HAVE_TM
  tm_term_double_setup((double)(global->mu_factor[level->depth] * op_in->mu),
                       (double)(global->mu_factor[level->depth] * op_in->mu_even_shift),
                       (double)(global->mu_factor[level->depth] * op_in->mu_odd_shift), &op, level);
#endif

#ifdef HAVE_TM1p1
  epsbar_term_double_setup(
      (double)(global->epsbar_factor[level->depth] * op_in->epsbar),
      (double)(global->epsbar_factor[level->depth] * op_in->epsbar_ig5_even_shift),
      (double)(global->epsbar_factor[level->depth] * op_in->epsbar_ig5_odd_shift), &op, level);
#endif

  boundaryUpdate();

  operator_double_set_couplings(&op, level);
}

void Schwarz::setupIntermediateOperator() {
  if (!level->idle) {
    conf_double_gather(&op, &level->op_double, level);
    Schwarz::boundaryUpdate();
    if (level->global->method >= 4 && level->odd_even) {
      coarse_oddeven_alloc_double(level);
      coarse_oddeven_setup_double(&op, _REORDER, level);
    }
    coarse_operator_double_set_couplings(&op, level);
  }
}

void Schwarz::setupCoarseOperator() {
  conf_double_gather(&(op), &(level->op_double), level);
  if (level->global->method >= 4 && level->odd_even) {
    coarse_oddeven_alloc_double(level);
    coarse_oddeven_setup_double(&op, _NO_REORDERING, level);
  } else {
    coarse_operator_double_set_couplings(&op, level);
  }
}
#endif
