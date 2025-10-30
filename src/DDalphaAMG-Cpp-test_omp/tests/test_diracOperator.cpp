//
// Created by shiting on 2024-10-04.
//
#include <main.h>
#ifndef BLEN
#define BLEN = 2
#endif

template<typename T>
static inline void check_relative_diff(std::complex<T>* v1, std::complex<T>* v2, int len, Level* l){
  vector_double_minus(v2, v2, v1, 0, len, l);
  double numerator = global_norm_double(v2, 0, len, l);
  double denominator = global_norm_double(v1, 0, len, l);
  ASSERT_LT(numerator/denominator, 1e-15);
}

// The fixture for testing class Level.
class DiracTest : public ::testing::Test {
protected:
  // You can remove any or all of the following functions if their bodies would
  // be empty.

  DiracTest() {
    // You can do set-up work for each test here.

  }

  ~DiracTest() override {
    // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  void SetUp() override {
    // Code here will be called immediately after the constructor (right
    // before each test).

    level = EnvironmentDDalphaAMG::level;

    vector_double_define(level[0].gmres.b, 1, 0, level[0].inner_vector_size);
    vector_double_define(level[0].gmres.w, 1, 0, level[0].inner_vector_size);
  }

  void TearDown() override {
    // Code here will be called immediately after each test (right
    // before the destructor).
  }

  // Class members declared here can be used by all tests in the test suite

  Level *level;

};

TEST_F(DiracTest, CheckGhostSend){
  int nv = level[0].num_lattice_site_var;
  int mu = _T;

  int offset = nv/2 * level[0].num_lattice_sites;
  operator_struct<double>* op = level[0].gmres.matrix;

  /// prepare the arrays to be transferred
  auto *prnT = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnTb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  vector_double_define_random(prnT, 0, BLEN * offset);
  vector_double_copy(prnTb, prnT, 0, BLEN * offset, &level[0]);
  transform_vector_block(prnTb, offset, nv/2, BLEN);

  auto buffer = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto bufferb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];

  for (int k = 0; k < BLEN; k++) {
    ghost_sendrecv_double(&prnT[k*offset], mu, -1, &(op->c), _FULL_SYSTEM, &level[0]);
    ghost_wait_double(&prnT[k*offset], mu, -1, &(op->c), _FULL_SYSTEM, &level[0]);
  }

  ghost_sendrecv_full(prnTb, mu, -1, &(op->c), _FULL_SYSTEM, &level[0], BLEN);
  ghost_wait_full(prnTb, mu, -1, &(op->c), _FULL_SYSTEM, &level[0], BLEN);

  separate_vectors_in_block(prnTb, offset, nv/2, BLEN);
  check_relative_diff(prnT, prnTb, BLEN * offset, &level[0]);

  /*
   * for testing the reordering of data
  int * table = op->c.boundary_table[2 * _T + 1];
  int num_boundary_sites = op->c.num_boundary_sites[2 * _T];
  for (int k = 0; k < BLEN; ++k) {
    auto phi_pt_start = &prnT[k*offset];
    for (int j = 0; j < num_boundary_sites; j++) {
      auto phi_pt = phi_pt_start + table[j] * 6;
      for (int i = 0; i < 6; i++) {
        buffer[i] = phi_pt[i];
      }
      buffer += 6;
    }
  }
  buffer -= 6 * num_boundary_sites * BLEN;

  for (int j = 0; j < num_boundary_sites; j++) {
    for (int k = 0; k < BLEN; k++) {
      for (int i= 0; i < 6; i++) {
        bufferb[j*BLEN*6 + k*6 + i] = prnTb[table[j]*BLEN*6 + k*6 + i];
      }
    }
  }
  separate_vectors_in_block(bufferb, 6 * num_boundary_sites, nv/2);

  check_relative_diff(buffer, bufferb, 6 * num_boundary_sites * BLEN, &level[0]);
  */

  delete[] prnT;
  delete[] buffer;
  delete[] prnTb;
  delete[] bufferb;

  auto *prpT = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpTb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  vector_double_define_random(prpT, 0, BLEN * offset);
  vector_double_copy(prpTb, prpT, 0, BLEN * offset, &level[0]);
  transform_vector_block(prpTb, offset, nv/2, BLEN);

  for (int k = 0; k < BLEN; k++) {
    ghost_sendrecv_double(&prpT[k*offset], mu, +1, &(op->c), _FULL_SYSTEM, &level[0]);
    ghost_wait_double(&prpT[k*offset], mu, +1, &(op->c), _FULL_SYSTEM, &level[0]);
  }

  ghost_sendrecv_full(prpTb, mu, +1, &(op->c), _FULL_SYSTEM, &level[0], BLEN);
  ghost_wait_full(prpTb, mu, +1, &(op->c), _FULL_SYSTEM, &level[0], BLEN);

  separate_vectors_in_block(prpTb, offset, nv/2, BLEN);
  check_relative_diff(prpT, prpTb, BLEN * offset, &level[0]);

  delete[] prpT;
  delete[] prpTb;
}

TEST_F(DiracTest, CheckBlockComm){
  int nv = level[0].num_lattice_site_var;
  int mu = _T;
  int dir = 1;

  int offset = nv/2 * level[0].num_lattice_sites;
  operator_struct<double>* op = level[0].gmres.matrix;

  /// prepare the arrays to be transferred
  auto *prnT = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnTb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  vector_double_define_random(prnT, 0, BLEN * offset);
  vector_double_copy(prnTb, prnT, 0, BLEN * offset, &level[0]);

  auto *vsout = new BlockVecPerSite<double>(prnTb, nv/2);
  vsout->transform(offset, BLEN);
  ghost_sendrecv_full(prnTb, mu, dir, &(op->c), _FULL_SYSTEM, &level[0], BLEN);
  ghost_wait_full(prnTb, mu, dir, &(op->c), _FULL_SYSTEM, &level[0], BLEN);
  vsout->detransform(offset, BLEN);

  auto *veout = new BlockVecPerEntry<double>(prnT);
  veout->transform(offset, BLEN);
  ghost_sendrecv_per_entry(prnT, mu, dir, &(op->c), _FULL_SYSTEM, &level[0], BLEN);
  ghost_wait_per_entry(prnT, mu, dir, &(op->c), _FULL_SYSTEM, &level[0], BLEN);
  veout->detransform(offset, BLEN);

  check_relative_diff(prnT, prnTb, BLEN * offset, &level[0]);

  delete[] prnTb;
  delete[] prnT;
}

TEST_F(DiracTest, EvaluateDiracOperator) {
  /// len = (4*n^3+n^4) * 12, in practice the length of vector should be n^4 * 12
  int len = level[0].gmres.vectorLength;
  auto* res = new std::complex<double>[len];
  int nv = level[0].num_lattice_site_var;

  /// sanity check for x=0
  level[0].gmres.applyOperator(res, level[0].gmres.x, level[0].gmres.matrix, &level[0]);
  for (int i = 0; i < len; ++i) {
    EXPECT_EQ(res[i], level[0].gmres.x[i]) << "Vectors x and y differ at index " << i;
  }
  delete[] res;

  /// initialize phi and eta
  auto* block_phi = new std::complex<double>[BLEN*len];
  auto* block_eta = new std::complex<double>[BLEN*len];
  auto* phi = new std::complex<double>[BLEN*len];
  auto* eta = new std::complex<double>[BLEN*len];
  // TODO:change the phi to something meaningful, here is random
  vector_double_define_random(block_phi, 0, len*BLEN);
  vector_double_copy(phi, block_phi, 0, len*BLEN, &level[0]);

  /// do normal apply dirac operator
  for (int k = 0; k < BLEN; k++) {
    level[0].gmres.applyOperator(&eta[k*len], &phi[k*len], level[0].gmres.matrix, &level[0]);
  }

  /// do block apply dirac operator
  auto *blockphi = new BlockVecPerSite<double>(block_phi, nv);
  auto *blocketa = new BlockVecPerSite<double>(block_eta, nv);
  blockphi->transform(len, BLEN);
  d_plus_clover_block(blocketa, blockphi, level[0].gmres.matrix, &level[0], BLEN);
  blocketa->detransform(len, BLEN);

  /// compare
  check_relative_diff(&blocketa->v[0], &eta[0], BLEN * len, &level[0]);

  delete[] block_phi;
  delete[] block_eta;
  delete[] phi;
  delete[] eta;
}

TEST_F(DiracTest, CheckBlockDirac){
  /// len = (4*n^3+n^4) * 12, in practice the length of vector should be n^4 * 12
  int len = level[0].gmres.vectorLength;
  int nv = level[0].num_lattice_site_var;

  /// initialize phi and eta
  auto* block_phi = new std::complex<double>[BLEN*len];
  auto* block_eta = new std::complex<double>[BLEN*len];
  auto* phi = new std::complex<double>[BLEN*len];
  auto* eta = new std::complex<double>[BLEN*len];
  // TODO:change the phi to something meaningful, here is random
  vector_double_define_random(block_phi, 0, len*BLEN);
  vector_double_copy(phi, block_phi, 0, len*BLEN, &level[0]);

  /// do block per site
  auto *vsphi = new BlockVecPerSite<double>(block_phi, nv);
  auto *vseta = new BlockVecPerSite<double>(block_eta, nv);
  vsphi->transform(len, BLEN);
  d_plus_clover_block(vseta, vsphi, level[0].gmres.matrix, &level[0], BLEN);
  vseta->detransform(len, BLEN);

  /// do block per entry
  auto *vephi = new BlockVecPerEntry<double>(phi);
  auto *veeta = new BlockVecPerEntry<double>(eta);
  vephi->transform(len, BLEN);
  d_plus_clover_block(veeta, vephi, level[0].gmres.matrix, &level[0], BLEN);
  veeta->detransform(len, BLEN);

  /// compare
  check_relative_diff(veeta->v, vseta->v, len*BLEN, &level[0]);

  delete[] block_phi;
  delete[] block_eta;
  delete[] phi;
  delete[] eta;
}

TEST_F(DiracTest, CheckTransform){
  const int len = level[0].gmres.vectorLength;
  int nv = 12;
  std::complex<double> vs[BLEN][len];

  /// define random numbers to vectors vs
  for(auto & v : vs){
    vector_double_define_random(v, 0, len);
  }
  /// concatenate vs to v
  auto* v = new std::complex<double>[len*BLEN];
  for (int i = 0; i < BLEN; i++) {
    std::copy(vs[i], vs[i]+len, v+i*len);
  }

  /// test transformation
  transform_vector_block(v, len, nv, BLEN);
  for (int i = 0; i < len/nv; i++) {
    for (int j = 0; j < BLEN; j++) {
      for (int k = 0; k < nv; k++) {
        EXPECT_EQ(v[i*BLEN*nv+j*nv+k], vs[j][i*nv+k]) << "problem occurs at index " << i;
      }
    }
  }

  /// test de-transformation
  separate_vectors_in_block(v, len, nv, BLEN);
  for (int i = 0; i < BLEN; i++) {
    for (int j = 0; j < len/nv; j++) {
      for (int k = 0; k < nv; k++) {
        EXPECT_EQ(v[i*len+j*nv+k], vs[i][j*nv+k]) << "problem occurs at index " << i;
      }
    }
  }

  delete[] v;
}

TEST_F(DiracTest, CheckClover){
  /// len = (4*n^3+n^4) * 12, in practice the length of vector should be n^4 * 12
  int len = level[0].gmres.vectorLength;
  int nv = level[0].num_lattice_site_var;
  int n = level[0].num_inner_lattice_sites;
  int end =  nv * n;

  /// initialize phi and eta
  auto* block_phi = new std::complex<double>[BLEN*len];
  auto* block_eta = new std::complex<double>[BLEN*len];
  auto* phi = new std::complex<double>[BLEN*len];
  auto* eta = new std::complex<double>[BLEN*len];
  // TODO:change the phi to something meaningful, here is random
  vector_double_define_random(block_phi, 0, len*BLEN);
  vector_double_copy(phi, block_phi, 0, BLEN * len, &level[0]);

  /// do normal clover
  for (int i = 0; i < BLEN; ++i) {
    clover_double(&eta[i*len], &phi[i*len], level[0].gmres.matrix, 0, end, &level[0]);
  }

  /// do block clover
  transform_vector_block(block_phi, len, nv, BLEN);
  clover_block(block_eta, block_phi, level[0].gmres.matrix, end, &level[0], BLEN);
  separate_vectors_in_block(block_eta, len, nv, BLEN);

  /// compare
  check_relative_diff(block_eta, eta, len*BLEN, &level[0]);

  delete[] eta;
  delete[] phi;
  delete[] block_phi;
  delete[] block_eta;
}

TEST_F(DiracTest, CheckBlockClover){
  /// len = (4*n^3+n^4) * 12, in practice the length of vector should be n^4 * 12
  int len = level[0].gmres.vectorLength;
  int nv = level[0].num_lattice_site_var;
  int n = level[0].num_inner_lattice_sites;
  int end =  nv * n;

  /// initialize phi and eta
  auto* first = new std::complex<double>[BLEN*len];
  auto* second = new std::complex<double>[BLEN*len];
  auto* out1 = new std::complex<double>[BLEN*len];
  auto* out2 = new std::complex<double>[BLEN*len];

  vector_double_define_random(first, 0, BLEN*len);
  vector_double_copy(second, first, 0, BLEN*len, &level[0]);

  auto *vpsite = new BlockVecPerSite<double>(first, nv);
  vpsite->transform(len, BLEN);
  clover_block(out1, first, level[0].gmres.matrix, end, &level[0], BLEN);
  auto *vsout = new BlockVecPerSite<double>(out1, nv);
  vsout->detransform(len, BLEN);

  auto*vpentry = new BlockVecPerEntry<double>(second);
  vpentry->transform(len, BLEN);
  clover_per_entry(out2, second, level[0].gmres.matrix, end, &level[0], BLEN);
  auto* vecout = new BlockVecPerEntry<double>(out2);
  vecout->detransform(len, BLEN);

  check_relative_diff(out1, out2, len*BLEN, &level[0]);

  delete[] first;
  delete[] second;
  delete[] out2;
  delete[] out1;
  delete vpsite;
  delete vsout;
  delete vpentry;
  delete vecout;
}


TEST_F(DiracTest, CheckPiMinus){
  int n = level[0].num_inner_lattice_sites, nv = level[0].num_lattice_site_var;
  int end =  nv * n;
  int len = level[0].gmres.vectorLength;

  /// initialize phi
  auto* block_phi = new std::complex<double>[BLEN*len];
  auto* phi = new std::complex<double>[BLEN*len];
  vector_double_define_random(block_phi, 0, BLEN*len);
  for (int i = 0; i < BLEN*len; i++) {
    phi[i] = block_phi[i];
  }
  /// transform to block_phi
  transform_vector_block(block_phi, len, nv, BLEN);

  /// do normal pi minus
  auto *prnT = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnZ = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnY = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnX = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  for (int k = 0; k < BLEN; ++k) {
    auto phi_pt = &phi[k*len];
    for (int i = 0; i < end / 2; i += 6, phi_pt += 12) {
      prp_T_double(prnT + i + k * nv/2 * level[0].num_lattice_sites, phi_pt);
      prp_Z_double(prnZ + i + k * nv/2 * level[0].num_lattice_sites, phi_pt);
      prp_Y_double(prnY + i + k * nv/2 * level[0].num_lattice_sites, phi_pt);
      prp_X_double(prnX + i + k * nv/2 * level[0].num_lattice_sites, phi_pt);
    }
  }

  /// do block pi minus
  auto *prnTb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnZb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnYb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnXb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  pi_minus(block_phi, prnTb, prnZb, prnYb, prnXb, 0, end*BLEN);

  /// transform block back and compare
  separate_vectors_in_block(prnTb, nv/2 * level[0].num_lattice_sites, nv/2, BLEN);
  separate_vectors_in_block(prnZb, nv/2 * level[0].num_lattice_sites, nv/2, BLEN);
  separate_vectors_in_block(prnYb, nv/2 * level[0].num_lattice_sites, nv/2, BLEN);
  separate_vectors_in_block(prnXb, nv/2 * level[0].num_lattice_sites, nv/2, BLEN);

  /// compare
  int offset = nv/2 * level[0].num_lattice_sites;
  check_relative_diff(prnTb, prnT, BLEN * offset, &level[0]);
  check_relative_diff(prnZb, prnZ, BLEN * offset, &level[0]);
  check_relative_diff(prnYb, prnY, BLEN * offset, &level[0]);
  check_relative_diff(prnXb, prnX, BLEN * offset, &level[0]);

  delete[] prnX;
  delete[] prnY;
  delete[] prnZ;
  delete[] prnT;
  delete[] prnXb;
  delete[] prnYb;
  delete[] prnZb;
  delete[] prnTb;
  delete[] block_phi;
  delete[] phi;
}

TEST_F(DiracTest, CheckBlockPiMinus){
  int n = level[0].num_inner_lattice_sites, nv = level[0].num_lattice_site_var;
  int end =  nv * n;
  int len = level[0].gmres.vectorLength;
  /// initialize phi
  auto *phi1 = new std::complex<double>[BLEN*len];
  auto *phi2 = new std::complex<double>[BLEN*len];
  vector_double_define_random(phi1, 0, BLEN*len);
  vector_copy(phi2, phi1, 0, len*BLEN);

  auto *vs = new BlockVecPerSite<double>(phi1, nv);
  vs->transform(len, BLEN);
  auto *ve = new BlockVecPerEntry<double>(phi2);
  ve->transform(len, BLEN);

  auto *prnTb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnZb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnYb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnXb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  pi_minus(phi1, prnTb, prnZb, prnYb, prnXb, 0, end*BLEN);
  /// wrap and de-transform
  auto *vsout_T = new BlockVecPerSite<double>(prnTb, nv/2);
  vsout_T->detransform(len/2, BLEN);
  auto *vsout_Z = new BlockVecPerSite<double>(prnZb, nv/2);
  vsout_Z->detransform(len/2, BLEN);
  auto *vsout_Y = new BlockVecPerSite<double>(prnYb, nv/2);
  vsout_Y->detransform(len/2, BLEN);
  auto *vsout_X = new BlockVecPerSite<double>(prnXb, nv/2);
  vsout_X->detransform(len/2, BLEN);

  auto *prnT = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnZ = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnY = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnX = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  pi_minus_per_entry(phi2, prnT, prnZ, prnY, prnX, 0, end, BLEN);
  /// wrap and de-transform
  auto *veout_T = new BlockVecPerEntry<double>(prnT);
  veout_T->detransform(len/2, BLEN);
  auto *veout_Z = new BlockVecPerEntry<double>(prnZ);
  veout_Z->detransform(len/2, BLEN);
  auto *veout_Y = new BlockVecPerEntry<double>(prnY);
  veout_Y->detransform(len/2, BLEN);
  auto *veout_X = new BlockVecPerEntry<double>(prnX);
  veout_X->detransform(len/2, BLEN);

  check_relative_diff(vsout_T->v, veout_T->v, len/2*BLEN, &level[0]);
  check_relative_diff(vsout_Z->v, veout_Z->v, len/2*BLEN, &level[0]);
  check_relative_diff(vsout_Y->v, veout_Y->v, len/2*BLEN, &level[0]);
  check_relative_diff(vsout_X->v, veout_X->v, len/2*BLEN, &level[0]);

  delete[] prnXb;
  delete[] prnYb;
  delete[] prnZb;
  delete[] prnTb;
  delete[] prnX;
  delete[] prnY;
  delete[] prnZ;
  delete[] prnT;
  delete[] phi1;
  delete[] phi2;
}

TEST_F(DiracTest, CheckPiPlusProduct){
  int n = level[0].num_inner_lattice_sites, nv = level[0].num_lattice_site_var;
  int end =  nv * n;
  int len = level[0].gmres.vectorLength;

  /// initialize phi
  auto* block_phi = new std::complex<double>[BLEN*len];
  auto* phi = new std::complex<double>[BLEN*len];
  vector_double_define_random(block_phi, 0, BLEN*len);
  for (int i = 0; i < BLEN*len; i++) {
    phi[i] = block_phi[i];
  }
  /// transform to block_phi
  transform_vector_block(block_phi, len, nv, BLEN);

  /// do normal pi plus product
  auto *prpT = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpZ = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpY = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpX = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  std::complex<double> *end_pt, *phi_pt;
  std::complex<double> pbuf[6];
  int j;
  int offset = nv/2 * level[0].num_lattice_sites;

  for (int k = 0; k < BLEN; ++k) {
    auto *D_pt = level[0].gmres.matrix->D;
    auto *nb_pt = level[0].gmres.matrix->neighbor_table;

    for (phi_pt = &phi[k*len], end_pt = phi_pt + end; phi_pt < end_pt; phi_pt += 12) {
      // T dir
      j = 6 * (*nb_pt);
      nb_pt++;
      prn_T_double(pbuf, phi_pt); /// (I4 + Gamma_mu) kronecker I3 * psi(x)
      mvmh_double(k*offset + prpT + j, D_pt, pbuf); /// apply U_mu^H(x), first half
      mvmh_double(k*offset + prpT + j + 3, D_pt, pbuf + 3); /// second half
      D_pt += 9;
      // Z dir
      j = 6 * (*nb_pt);
      nb_pt++;
      prn_Z_double(pbuf, phi_pt);
      mvmh_double(k*offset + prpZ + j, D_pt, pbuf);
      mvmh_double(k*offset + prpZ + j + 3, D_pt, pbuf + 3);
      D_pt += 9;
      // Y dir
      j = 6 * (*nb_pt);
      nb_pt++;
      prn_Y_double(pbuf, phi_pt);
      mvmh_double(k*offset + prpY + j, D_pt, pbuf);
      mvmh_double(k*offset + prpY + j + 3, D_pt, pbuf + 3);
      D_pt += 9;
      // X dir
      j = 6 * (*nb_pt);
      nb_pt++;
      prn_X_double(pbuf, phi_pt);
      mvmh_double(k*offset + prpX + j, D_pt, pbuf);
      mvmh_double(k*offset + prpX + j + 3, D_pt, pbuf + 3);
      D_pt += 9;
    }
  }

  /// do block pi plus product
  auto *prpTb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpZb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpYb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpXb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  pi_plus_kron_UH_psi(block_phi, level[0].gmres.matrix->D, prpTb, prpZb, prpYb, prpXb,
                      level[0].gmres.matrix->neighbor_table, 0, end, BLEN);

  separate_vectors_in_block(prpTb, offset, nv/2, BLEN);
  separate_vectors_in_block(prpZb, offset, nv/2, BLEN);
  separate_vectors_in_block(prpYb, offset, nv/2, BLEN);
  separate_vectors_in_block(prpXb, offset, nv/2, BLEN);

  /// compare
  check_relative_diff(prpTb, prpT, BLEN * offset, &level[0]);
  check_relative_diff(prpZb, prpZ, BLEN * offset, &level[0]);
  check_relative_diff(prpYb, prpY, BLEN * offset, &level[0]);
  check_relative_diff(prpXb, prpX, BLEN * offset, &level[0]);

  delete[] prpT;
  delete[] prpZ;
  delete[] prpX;
  delete[] prpY;
  delete[] prpTb;
  delete[] prpZb;
  delete[] prpYb;
  delete[] prpXb;
}

TEST_F(DiracTest, CheckBlockPiPlus){
  int n = level[0].num_inner_lattice_sites, nv = level[0].num_lattice_site_var;
  int end =  nv * n;
  int len = level[0].gmres.vectorLength;

  /// initialize phi
  auto *phi1 = new std::complex<double>[BLEN*len];
  auto *phi2 = new std::complex<double>[BLEN*len];
  vector_double_define_random(phi1, 0, BLEN*len);
  vector_copy(phi2, phi1, 0, len*BLEN);

  auto *vs = new BlockVecPerSite<double>(phi1, nv);
  vs->transform(len, BLEN);
  auto *ve = new BlockVecPerEntry<double>(phi2);
  ve->transform(len, BLEN);

  auto *prpTb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpZb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpYb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpXb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  pi_plus_kron_UH_psi(phi1, level[0].gmres.matrix->D, prpTb, prpZb, prpYb, prpXb,
                      level[0].gmres.matrix->neighbor_table, 0, end, BLEN);
  /// wrap and de-transform
  auto *vsout_T = new BlockVecPerSite<double>(prpTb, nv/2);
  vsout_T->detransform(len/2, BLEN);
  auto *vsout_Z = new BlockVecPerSite<double>(prpZb, nv/2);
  vsout_Z->detransform(len/2, BLEN);
  auto *vsout_Y = new BlockVecPerSite<double>(prpYb, nv/2);
  vsout_Y->detransform(len/2, BLEN);
  auto *vsout_X = new BlockVecPerSite<double>(prpXb, nv/2);
  vsout_X->detransform(len/2, BLEN);

  auto *prpT = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpZ = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpY = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpX = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  pi_plus_per_entry(phi2, level[0].gmres.matrix->D, prpT, prpZ, prpY, prpX,
                    level[0].gmres.matrix->neighbor_table, 0, end, BLEN);
  /// wrap and de-transform
  auto *veout_T = new BlockVecPerEntry<double>(prpT);
  veout_T->detransform(len/2, BLEN);
  auto *veout_Z = new BlockVecPerEntry<double>(prpZ);
  veout_Z->detransform(len/2, BLEN);
  auto *veout_Y = new BlockVecPerEntry<double>(prpY);
  veout_Y->detransform(len/2, BLEN);
  auto *veout_X = new BlockVecPerEntry<double>(prpX);
  veout_X->detransform(len/2, BLEN);

  check_relative_diff(vsout_T->v, veout_T->v, len/2*BLEN, &level[0]);
  check_relative_diff(vsout_Z->v, veout_Z->v, len/2*BLEN, &level[0]);
  check_relative_diff(vsout_Y->v, veout_Y->v, len/2*BLEN, &level[0]);
  check_relative_diff(vsout_X->v, veout_X->v, len/2*BLEN, &level[0]);

  delete[] prpXb;
  delete[] prpYb;
  delete[] prpZb;
  delete[] prpTb;
  delete[] prpX;
  delete[] prpY;
  delete[] prpZ;
  delete[] prpT;
  delete[] phi1;
  delete[] phi2;
}

TEST_F(DiracTest, CheckApplyUPsi){
  int n = level[0].num_inner_lattice_sites, nv = level[0].num_lattice_site_var;
  int end =  nv * n;
  int len = level[0].gmres.vectorLength;

  auto *prnT = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnZ = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnY = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnX = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];

  auto *prnTb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnZb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnYb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnXb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];

  int offset = nv/2 * level[0].num_lattice_sites;

  /// initialize prn
  vector_double_define_random(prnT, 0, BLEN*offset);
  vector_double_define_random(prnZ, 0, BLEN*offset);
  vector_double_define_random(prnY, 0, BLEN*offset);
  vector_double_define_random(prnX, 0, BLEN*offset);
  for (int i = 0; i < BLEN*offset; i++) {
    prnTb[i] = prnT[i];
    prnZb[i] = prnZ[i];
    prnYb[i] = prnY[i];
    prnXb[i] = prnX[i];
  }
  transform_vector_block(prnTb, offset, nv/2, BLEN);
  transform_vector_block(prnZb, offset, nv/2, BLEN);
  transform_vector_block(prnYb, offset, nv/2, BLEN);
  transform_vector_block(prnXb, offset, nv/2, BLEN);

  /// apply U psi(x+mu) normally
  std::complex<double> *eta_pt, *end_pt;
  int j;
  std::complex<double> pbuf[6];
  auto *eta = new std::complex<double>[BLEN * len];
  for (int k = 0; k < BLEN; k++) {
    auto *D_pt = level[0].gmres.matrix->D;
    auto *nb_pt = level[0].gmres.matrix->neighbor_table;
    for (eta_pt = &eta[k*len], end_pt = eta_pt + end; eta_pt < end_pt; eta_pt += 12) {
      // T dir
      j = 6 * (*nb_pt);
      nb_pt++;
      mvm_double(pbuf, D_pt, k*offset + prnT + j); /// apply U_mu(x) to prn, first half of x
      mvm_double(pbuf + 3, D_pt, k*offset + prnT + j + 3); /// second half of x
      /// each line of second half of eta can be represented by lines of first half of eta,
      /// here we apply the effect of pi_mu_minus kronecker U_mu(x) * psi(x+mu_hat)
      pbp_su3_T_double(pbuf, eta_pt);
      D_pt += 9;
      // Z dir
      j = 6 * (*nb_pt);
      nb_pt++;
      mvm_double(pbuf, D_pt, k*offset + prnZ + j);
      mvm_double(pbuf + 3, D_pt, k*offset + prnZ + j + 3);
      pbp_su3_Z_double(pbuf, eta_pt);
      D_pt += 9;
      // Y dir
      j = 6 * (*nb_pt);
      nb_pt++;
      mvm_double(pbuf, D_pt, k*offset + prnY + j);
      mvm_double(pbuf + 3, D_pt, k*offset + prnY + j + 3);
      pbp_su3_Y_double(pbuf, eta_pt);
      D_pt += 9;
      // X dir
      j = 6 * (*nb_pt);
      nb_pt++;
      mvm_double(pbuf, D_pt, k*offset + prnX + j);
      mvm_double(pbuf + 3, D_pt, k*offset + prnX + j + 3);
      pbp_su3_X_double(pbuf, eta_pt);
      D_pt += 9;
    }
  }

  auto *etab = new std::complex<double>[BLEN * len];
  apply_U_psi_plus_mu(etab, level[0].gmres.matrix->D, prnTb, prnZb, prnYb, prnXb,
                      level[0].gmres.matrix->neighbor_table, 0, end, BLEN);

  separate_vectors_in_block(etab, len, nv, BLEN);

  check_relative_diff(etab, eta, BLEN * len, &level[0]);

  delete[] prnX;
  delete[] prnY;
  delete[] prnZ;
  delete[] prnT;
  delete[] prnXb;
  delete[] prnYb;
  delete[] prnZb;
  delete[] prnTb;
  delete[] eta;
  delete[] etab;
}

TEST_F(DiracTest, CheckBlockApplyU){
  int n = level[0].num_inner_lattice_sites, nv = level[0].num_lattice_site_var;
  int end =  nv * n;
  int len = level[0].gmres.vectorLength;

  /// initialize prns
  auto *prnTb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnZb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnYb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnXb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnT = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnZ = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnY = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prnX = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];

  vector_double_define_random(prnTb, 0, BLEN*len/2);
  vector_double_define_random(prnZb, 0, BLEN*len/2);
  vector_double_define_random(prnYb, 0, BLEN*len/2);
  vector_double_define_random(prnXb, 0, BLEN*len/2);
  vector_copy(prnT, prnTb, 0, BLEN*len/2);
  vector_copy(prnZ, prnZb, 0, BLEN*len/2);
  vector_copy(prnY, prnYb, 0, BLEN*len/2);
  vector_copy(prnX, prnXb, 0, BLEN*len/2);

  auto *vsout_T = new BlockVecPerSite<double>(prnTb, nv/2);
  vsout_T->transform(len/2, BLEN);
  auto *vsout_Z = new BlockVecPerSite<double>(prnZb, nv/2);
  vsout_Z->transform(len/2, BLEN);
  auto *vsout_Y = new BlockVecPerSite<double>(prnYb, nv/2);
  vsout_Y->transform(len/2, BLEN);
  auto *vsout_X = new BlockVecPerSite<double>(prnXb, nv/2);
  vsout_X->transform(len/2, BLEN);

  /// initialize etas
  auto *etab = new std::complex<double>[BLEN * len];
  auto *eta = new std::complex<double>[BLEN*len];
  vector_double_define_random(etab, 0, BLEN*len);
  vector_copy(eta, etab, 0, BLEN*len);
  auto *vseta = new BlockVecPerSite<double>(etab, nv);
  auto *veeta = new BlockVecPerEntry<double>(eta);

  /// do per site
  vseta->transform(len, BLEN);
  apply_U_psi_plus_mu(etab, level[0].gmres.matrix->D, prnTb, prnZb, prnYb, prnXb,
                      level[0].gmres.matrix->neighbor_table, 0, end, BLEN);
  vseta->detransform(len, BLEN);

  /// prepare prns for per entry
  auto *veout_T = new BlockVecPerEntry<double>(prnT);
  veout_T->transform(len/2, BLEN);
  auto *veout_Z = new BlockVecPerEntry<double>(prnZ);
  veout_Z->transform(len/2, BLEN);
  auto *veout_Y = new BlockVecPerEntry<double>(prnY);
  veout_Y->transform(len/2, BLEN);
  auto *veout_X = new BlockVecPerEntry<double>(prnX);
  veout_X->transform(len/2, BLEN);

  /// do per entry
  veeta->transform(len, BLEN);
  apply_U_per_entry(eta, level[0].gmres.matrix->D, prnT, prnZ, prnY, prnX,
                    level[0].gmres.matrix->neighbor_table, 0, end, BLEN);
  veeta->detransform(len, BLEN);

  check_relative_diff(vseta->v, veeta->v, len*BLEN, &level[0]);

  delete[] prnXb;
  delete[] prnYb;
  delete[] prnZb;
  delete[] prnTb;
  delete[] prnX;
  delete[] prnY;
  delete[] prnZ;
  delete[] prnT;
  delete[] eta;
  delete[] etab;
}

TEST_F(DiracTest, CheckLiftUpPlusDir){
  int n = level[0].num_inner_lattice_sites, nv = level[0].num_lattice_site_var;
  int end =  nv * n;
  int len = level[0].gmres.vectorLength;

  auto *prpT = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpZ = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpY = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpX = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];

  auto *prpTb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpZb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpYb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpXb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];

  int offset = nv/2 * level[0].num_lattice_sites;

  /// initialize prp
  vector_double_define_random(prpT, 0, BLEN*offset);
  vector_double_define_random(prpZ, 0, BLEN*offset);
  vector_double_define_random(prpY, 0, BLEN*offset);
  vector_double_define_random(prpX, 0, BLEN*offset);
  for (int i = 0; i < BLEN*offset; i++) {
    prpTb[i] = prpT[i];
    prpZb[i] = prpZ[i];
    prpYb[i] = prpY[i];
    prpXb[i] = prpX[i];
  }
  transform_vector_block(prpTb, offset, nv/2, BLEN);
  transform_vector_block(prpZb, offset, nv/2, BLEN);
  transform_vector_block(prpYb, offset, nv/2, BLEN);
  transform_vector_block(prpXb, offset, nv/2, BLEN);

  /// do normal lift up
  auto *eta = new std::complex<double>[BLEN * len];
  for (int k = 0; k < BLEN; k++) {
    std::complex<double> *eta_pt = &eta[k*len];
    for (int i = 0; i < end / 2; i += 6, eta_pt += 12) {
      /// each line of second half of eta can be represented by lines of first half of eta,
      /// here we apply the effect of pi_mu_plus kronecker U_mu(x-mu_hat) * psi(x-mu_hat)
      pbn_su3_T_double(k*offset + prpT + i, eta_pt);
      pbn_su3_Z_double(k*offset + prpZ + i, eta_pt);
      pbn_su3_Y_double(k*offset + prpY + i, eta_pt);
      pbn_su3_X_double(k*offset + prpX + i, eta_pt);
    }
  }

  auto *etab = new std::complex<double>[BLEN * len];
  lift_pi_plus_dir(etab, prpTb, prpZb, prpYb, prpXb, 0, end*BLEN);

  separate_vectors_in_block(etab, len, nv, BLEN);

  check_relative_diff(etab, eta, BLEN * len, &level[0]);

  delete[] prpT;
  delete[] prpZ;
  delete[] prpX;
  delete[] prpY;
  delete[] prpTb;
  delete[] prpZb;
  delete[] prpYb;
  delete[] prpXb;
}

TEST_F(DiracTest, CheckBlockLiftUp){
  int n = level[0].num_inner_lattice_sites, nv = level[0].num_lattice_site_var;
  int end =  nv * n;
  int len = level[0].gmres.vectorLength;

  /// initialize prns
  auto *prpTb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpZb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpYb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpXb = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpT = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpZ = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpY = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];
  auto *prpX = new std::complex<double>[BLEN * nv/2 * level[0].num_lattice_sites];

  vector_double_define_random(prpTb, 0, BLEN*len/2);
  vector_double_define_random(prpZb, 0, BLEN*len/2);
  vector_double_define_random(prpYb, 0, BLEN*len/2);
  vector_double_define_random(prpXb, 0, BLEN*len/2);
  vector_copy(prpT, prpTb, 0, BLEN*len/2);
  vector_copy(prpZ, prpZb, 0, BLEN*len/2);
  vector_copy(prpY, prpYb, 0, BLEN*len/2);
  vector_copy(prpX, prpXb, 0, BLEN*len/2);

  auto *vsout_T = new BlockVecPerSite<double>(prpTb, nv/2);
  vsout_T->transform(len/2, BLEN);
  auto *vsout_Z = new BlockVecPerSite<double>(prpZb, nv/2);
  vsout_Z->transform(len/2, BLEN);
  auto *vsout_Y = new BlockVecPerSite<double>(prpYb, nv/2);
  vsout_Y->transform(len/2, BLEN);
  auto *vsout_X = new BlockVecPerSite<double>(prpXb, nv/2);
  vsout_X->transform(len/2, BLEN);

  /// initialize etas
  auto *etab = new std::complex<double>[BLEN * len];
  auto *eta = new std::complex<double>[BLEN*len];
  vector_double_define_random(etab, 0, BLEN*len);
  vector_copy(eta, etab, 0, BLEN*len);
  auto *vseta = new BlockVecPerSite<double>(etab, nv);
  auto *veeta = new BlockVecPerEntry<double>(eta);

  /// do per site
  vseta->transform(len, BLEN);
  lift_pi_plus_dir(etab, prpTb, prpZb, prpYb, prpXb, 0, end*BLEN);
  vseta->detransform(len, BLEN);

  /// prepare prns for per entry
  auto *veout_T = new BlockVecPerEntry<double>(prpT);
  veout_T->transform(len/2, BLEN);
  auto *veout_Z = new BlockVecPerEntry<double>(prpZ);
  veout_Z->transform(len/2, BLEN);
  auto *veout_Y = new BlockVecPerEntry<double>(prpY);
  veout_Y->transform(len/2, BLEN);
  auto *veout_X = new BlockVecPerEntry<double>(prpX);
  veout_X->transform(len/2, BLEN);

  /// do per entry
  veeta->transform(len, BLEN);
  lift_pi_plus_per_entry(eta, prpT, prpZ, prpY, prpX, 0, end, BLEN);
  veeta->detransform(len, BLEN);

  check_relative_diff(vseta->v, veeta->v, len*BLEN, &level[0]);

  delete[] prpXb;
  delete[] prpYb;
  delete[] prpZb;
  delete[] prpTb;
  delete[] prpX;
  delete[] prpY;
  delete[] prpZ;
  delete[] prpT;
  delete[] eta;
  delete[] etab;
}

TEST_F(DiracTest, MeasurementTest) {
  /// len = (4*n^3+n^4) * 12, in practice the length of vector should be n^4 * 12
  int len = level[0].gmres.vectorLength;
  int nv = level[0].num_lattice_site_var;

  double t0;
  double t_nb = 0.0, t_b = 0.0, t_be = 0.0;
  const int measurement_round = 2;

  /// initialize phi and eta
  auto* phi2 = new std::complex<double>[BLEN*len];
  auto* eta2 = new std::complex<double>[BLEN*len];
  auto* block_phi = new std::complex<double>[BLEN*len];
  auto* block_eta = new std::complex<double>[BLEN*len];
  auto* phi = new std::complex<double>[BLEN*len];
  auto* eta = new std::complex<double>[BLEN*len];

  for (int i = 0; i < measurement_round; ++i) {
    memset(eta, 0, BLEN * len * sizeof(std::complex<double>));
    memset(block_eta, 0, BLEN * len * sizeof(std::complex<double>));
    memset(eta2, 0, BLEN * len * sizeof(std::complex<double>));
    vector_double_define_random(block_phi, 0, len*BLEN);
    vector_double_copy(phi, block_phi, 0, len*BLEN, &level[0]);
    vector_copy(phi2, block_phi, 0, len*BLEN);

    /// preparing the objects for block D_W
    BlockVecPerSite<double> *veta[BLEN];
    BlockVecPerSite<double> *vphi[BLEN];
    auto *vbeta = new BlockVecPerSite<double>(block_eta, 12);
    auto *vbphi = new BlockVecPerSite<double>(block_phi, 12);
    auto *veeta = new BlockVecPerEntry<double>(eta2);
    auto *vephi = new BlockVecPerEntry<double>(phi2);

    for (int k = 0; k < BLEN; ++k) {
      veta[k] = new BlockVecPerSite<double>(&eta[k*len], 12);
      vphi[k] = new BlockVecPerSite<double>(&phi[k*len], 12);
    }

    t0 = MPI_Wtime();
    for (int k = 0; k < BLEN; k++) {
      d_plus_clover_block(veta[k], vphi[k], level[0].gmres.matrix, &level[0], 1);
    }
    t_nb += MPI_Wtime() - t0;

    vbphi->transform(len, BLEN);
    t0 = MPI_Wtime();
    d_plus_clover_block(vbeta, vbphi, level[0].gmres.matrix, &level[0], BLEN);
    t_b += MPI_Wtime() - t0;
    vbeta->detransform(len, BLEN);

    vephi->transform(len, BLEN);
    t0 = MPI_Wtime();
    d_plus_clover_block(veeta, vephi, level[0].gmres.matrix, &level[0], BLEN);
    t_be += MPI_Wtime() - t0;
    veeta->detransform(len, BLEN);

    check_relative_diff(&vbeta->v[0], &veta[0]->v[0], BLEN * len, &level[0]);
    check_relative_diff(&vbeta->v[0], &veeta->v[0], BLEN * len, &level[0]);
  }

  printf0("No blocking for %d block length: %-8.4lf seconds, measurement rounds: %d\n",
          BLEN, t_nb, measurement_round);
  printf0("blocking for %d block length: %-8.4lf seconds, measurement rounds: %d\n",
          BLEN, t_b, measurement_round);
  printf0("blocking for %d block per entry length: %-8.4lf seconds, measurement rounds: %d\n",
          BLEN, t_be, measurement_round);

  delete[] block_phi;
  delete[] block_eta;
  delete[] phi;
  delete[] eta;
}