//
// Created by shiting on 2024-10-04.
//
#include <main.h>
#include "papi.h"
#define PAPI_CHECK(PAPI_CMD, MSG)                                              \
  do {                                                                         \
    int retval = (PAPI_CMD);                                                   \
    if ((retval) != PAPI_OK) {                                                 \
      PAPI_perror(MSG);                                                        \
    }                                                                          \
  } while (0);
const int blen = 2;

template<typename T>
static inline void check_relative_diff(std::complex<T>* v1, std::complex<T>* v2, int len, Level* l){
  vector_double_minus(v2, v2, v1, 0, len, l);
  double numerator = global_norm_double(v2, 0, len, l);
  double denominator = global_norm_double(v1, 0, len, l);
  ASSERT_LT(numerator/denominator, 1e-15);
}

int papi_start(const char** eventNames, int eventNum){
  int retval;

  static int eventSet = PAPI_NULL;
  int EventCode = PAPI_NULL;

  retval=PAPI_library_init(PAPI_VER_CURRENT);
  if (retval!=PAPI_VER_CURRENT) {
    fprintf(stderr,"Error initializing PAPI! %s\n",
            PAPI_strerror(retval));
    exit(1);
  }

  PAPI_CHECK(PAPI_create_eventset(&eventSet), "Error to create PAPI event set");

  for(int i = 0; i < eventNum; i++){
    PAPI_CHECK(PAPI_event_name_to_code(eventNames[i],&EventCode), "Error to convert name to code");

    PAPI_CHECK(PAPI_add_event(eventSet, EventCode),
               "Error adding PAPI event");
  }

  PAPI_CHECK(PAPI_start(eventSet), "Failed starting PAPI counters.");

  return eventSet;
}

void papi_stop(int eventSet){
  int num_events = PAPI_num_events(eventSet);
  if (num_events <= 0) {
    PAPI_perror("Failed determining the number of PAPI counters.");
  }
  long long *counters = (long long *)malloc(num_events * sizeof(long long));

  PAPI_stop(eventSet, counters);

  FILE *fh = stdout;
  int events[num_events];
  PAPI_list_events(eventSet, events, &num_events);

  for (int i = 0; i < num_events; i++) {
    char eventName[PAPI_MAX_STR_LEN];
    if(PAPI_event_code_to_name(events[i], eventName) != PAPI_OK){
      PAPI_perror("Error transforming code to name");
    }
    fprintf(fh, "%s\t", eventName);
  }
  fprintf(fh, "\n");

  for (int i = 0; i < num_events; i++) {
    fprintf(fh, "%lld\t", counters[i]);
  }
  fprintf(fh, "\n");

  free(counters);
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
  auto *prnT = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnTb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  vector_double_define_random(prnT, 0, blen * offset);
  vector_double_copy(prnTb, prnT, 0, blen * offset, &level[0]);
  transform_vector_block(prnTb, offset, nv/2, blen);

  auto buffer = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto bufferb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];

  for (int k = 0; k < blen; k++) {
    ghost_sendrecv_double(&prnT[k*offset], mu, -1, &(op->c), _FULL_SYSTEM, &level[0]);
    ghost_wait_double(&prnT[k*offset], mu, -1, &(op->c), _FULL_SYSTEM, &level[0]);
  }

  ghost_sendrecv_full(prnTb, mu, -1, &(op->c), &level[0], blen);
  ghost_wait_full(prnTb, mu, -1, &(op->c), &level[0], blen);

  separate_vectors_in_block(prnTb, offset, nv/2, blen);
  check_relative_diff(prnT, prnTb, blen * offset, &level[0]);

  /*
   * for testing the reordering of data
  int * table = op->c.boundary_table[2 * _T + 1];
  int num_boundary_sites = op->c.num_boundary_sites[2 * _T];
  for (int k = 0; k < blen; ++k) {
    auto phi_pt_start = &prnT[k*offset];
    for (int j = 0; j < num_boundary_sites; j++) {
      auto phi_pt = phi_pt_start + table[j] * 6;
      for (int i = 0; i < 6; i++) {
        buffer[i] = phi_pt[i];
      }
      buffer += 6;
    }
  }
  buffer -= 6 * num_boundary_sites * blen;

  for (int j = 0; j < num_boundary_sites; j++) {
    for (int k = 0; k < blen; k++) {
      for (int i= 0; i < 6; i++) {
        bufferb[j*blen*6 + k*6 + i] = prnTb[table[j]*blen*6 + k*6 + i];
      }
    }
  }
  separate_vectors_in_block(bufferb, 6 * num_boundary_sites, nv/2);

  check_relative_diff(buffer, bufferb, 6 * num_boundary_sites * blen, &level[0]);
  */

  delete[] prnT;
  delete[] buffer;
  delete[] prnTb;
  delete[] bufferb;

  auto *prpT = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpTb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  vector_double_define_random(prpT, 0, blen * offset);
  vector_double_copy(prpTb, prpT, 0, blen * offset, &level[0]);
  transform_vector_block(prpTb, offset, nv/2, blen);

  for (int k = 0; k < blen; k++) {
    ghost_sendrecv_double(&prpT[k*offset], mu, +1, &(op->c), _FULL_SYSTEM, &level[0]);
    ghost_wait_double(&prpT[k*offset], mu, +1, &(op->c), _FULL_SYSTEM, &level[0]);
  }

  ghost_sendrecv_full(prpTb, mu, +1, &(op->c), &level[0], blen);
  ghost_wait_full(prpTb, mu, +1, &(op->c), &level[0], blen);

  separate_vectors_in_block(prpTb, offset, nv/2, blen);
  check_relative_diff(prpT, prpTb, blen * offset, &level[0]);

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
  auto *prnT = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnTb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  vector_double_define_random(prnT, 0, blen * offset);
  vector_double_copy(prnTb, prnT, 0, blen * offset, &level[0]);

  auto *vsout = new BlockVecPerSite<double>(prnTb, nv/2);
  vsout->transform(offset, blen);
  ghost_sendrecv_full(prnTb, mu, dir, &(op->c), &level[0], blen);
  ghost_wait_full(prnTb, mu, dir, &(op->c), &level[0], blen);
  vsout->detransform(offset, blen);

  auto *veout = new BlockVecPerEntry<double>(prnT);
  veout->transform(offset, blen);
  ghost_sendrecv_per_entry(prnT, mu, dir, &(op->c), &level[0], blen);
  ghost_wait_per_entry(prnT, mu, dir, &(op->c), &level[0], blen);
  veout->detransform(offset, blen);

  check_relative_diff(prnT, prnTb, blen * offset, &level[0]);

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
  auto* block_phi = new std::complex<double>[blen*len];
  auto* block_eta = new std::complex<double>[blen*len];
  auto* phi = new std::complex<double>[blen*len];
  auto* eta = new std::complex<double>[blen*len];
  // TODO:change the phi to something meaningful, here is random
  vector_double_define_random(block_phi, 0, len*blen);
  vector_double_copy(phi, block_phi, 0, len*blen, &level[0]);

  /// do normal apply dirac operator
  for (int k = 0; k < blen; k++) {
    level[0].gmres.applyOperator(&eta[k*len], &phi[k*len], level[0].gmres.matrix, &level[0]);
  }

  /// do block apply dirac operator
  auto *blockphi = new BlockVecPerSite<double>(block_phi, nv);
  auto *blocketa = new BlockVecPerSite<double>(block_eta, nv);
  blockphi->transform(len, blen);
  d_plus_clover_block(blocketa, blockphi, level[0].gmres.matrix, &level[0], blen);
  blocketa->detransform(len, blen);

  /// compare
  check_relative_diff(&blocketa->v[0], &eta[0], blen * len, &level[0]);

  delete[] block_phi;
  delete[] block_eta;
  delete[] phi;
  delete[] eta;
}

TEST_F(DiracTest, CheckBlockDirac){
  /// len = (4*n^3+n^4) * 12, in practice the length of vector should be n^4 * 12
  int len = level[0].gmres.vectorLength;
  auto* res = new std::complex<double>[len];
  int nv = level[0].num_lattice_site_var;

  /// initialize phi and eta
  auto* block_phi = new std::complex<double>[blen*len];
  auto* block_eta = new std::complex<double>[blen*len];
  auto* phi = new std::complex<double>[blen*len];
  auto* eta = new std::complex<double>[blen*len];
  // TODO:change the phi to something meaningful, here is random
  vector_double_define_random(block_phi, 0, len*blen);
  vector_double_copy(phi, block_phi, 0, len*blen, &level[0]);

  /// do block per site
  auto *vsphi = new BlockVecPerSite<double>(block_phi, nv);
  auto *vseta = new BlockVecPerSite<double>(block_eta, nv);
  vsphi->transform(len, blen);
  d_plus_clover_block(vseta, vsphi, level[0].gmres.matrix, &level[0], blen);
  vseta->detransform(len, blen);

  /// do block per entry
  auto *vephi = new BlockVecPerEntry<double>(phi);
  auto *veeta = new BlockVecPerEntry<double>(eta);
  vephi->transform(len, blen);
  d_plus_clover_block(veeta, vephi, level[0].gmres.matrix, &level[0], blen);
  veeta->detransform(len, blen);

  /// compare
  check_relative_diff(veeta->v, vseta->v, len*blen, &level[0]);

  delete[] block_phi;
  delete[] block_eta;
  delete[] phi;
  delete[] eta;
}

TEST_F(DiracTest, CheckTransform){
  const int len = level[0].gmres.vectorLength;
  int nv = 12;
  std::complex<double> vs[blen][len];

  /// define random numbers to vectors vs
  for(auto & v : vs){
    vector_double_define_random(v, 0, len);
  }
  /// concatenate vs to v
  auto* v = new std::complex<double>[len*blen];
  for (int i = 0; i < blen; i++) {
    std::copy(vs[i], vs[i]+len, v+i*len);
  }

  /// test transformation
  transform_vector_block(v, len, nv, blen);
  for (int i = 0; i < len/nv; i++) {
    for (int j = 0; j < blen; j++) {
      for (int k = 0; k < nv; k++) {
        EXPECT_EQ(v[i*blen*nv+j*nv+k], vs[j][i*nv+k]) << "problem occurs at index " << i;
      }
    }
  }

  /// test de-transformation
  separate_vectors_in_block(v, len, nv, blen);
  for (int i = 0; i < blen; i++) {
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
  auto* block_phi = new std::complex<double>[blen*len];
  auto* block_eta = new std::complex<double>[blen*len];
  auto* phi = new std::complex<double>[blen*len];
  auto* eta = new std::complex<double>[blen*len];
  // TODO:change the phi to something meaningful, here is random
  vector_double_define_random(block_phi, 0, len*blen);
  vector_double_copy(phi, block_phi, 0, blen * len, &level[0]);

  /// do normal clover
  for (int i = 0; i < blen; ++i) {
    clover_double(&eta[i*len], &phi[i*len], level[0].gmres.matrix, 0, end, &level[0]);
  }

  /// do block clover
  transform_vector_block(block_phi, len, nv, blen);
  clover_block(block_eta, block_phi, level[0].gmres.matrix, end, &level[0], blen);
  separate_vectors_in_block(block_eta, len, nv, blen);

  /// compare
  check_relative_diff(block_eta, eta, len*blen, &level[0]);

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
  auto* first = new std::complex<double>[blen*len];
  auto* second = new std::complex<double>[blen*len];
  auto* out1 = new std::complex<double>[blen*len];
  auto* out2 = new std::complex<double>[blen*len];

  vector_double_define_random(first, 0, blen*len);
  vector_double_copy(second, first, 0, blen*len, &level[0]);

  auto *vpsite = new BlockVecPerSite<double>(first, nv);
  vpsite->transform(len, blen);
  clover_block(out1, first, level[0].gmres.matrix, end, &level[0], blen);
  auto *vsout = new BlockVecPerSite<double>(out1, nv);
  vsout->detransform(len, blen);

  auto*vpentry = new BlockVecPerEntry<double>(second);
  vpentry->transform(len, blen);
  clover_per_entry(out2, second, level[0].gmres.matrix, end, &level[0], blen);
  auto* vecout = new BlockVecPerEntry<double>(out2);
  vecout->detransform(len, blen);

  check_relative_diff(out1, out2, len*blen, &level[0]);

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
  auto* block_phi = new std::complex<double>[blen*len];
  auto* phi = new std::complex<double>[blen*len];
  vector_double_define_random(block_phi, 0, blen*len);
  for (int i = 0; i < blen*len; i++) {
    phi[i] = block_phi[i];
  }
  /// transform to block_phi
  transform_vector_block(block_phi, len, nv, blen);

  /// do normal pi minus
  auto *prnT = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnZ = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnY = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnX = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  for (int k = 0; k < blen; ++k) {
    auto phi_pt = &phi[k*len];
    for (int i = 0; i < end / 2; i += 6, phi_pt += 12) {
      prp_T_double(prnT + i + k * nv/2 * level[0].num_lattice_sites, phi_pt);
      prp_Z_double(prnZ + i + k * nv/2 * level[0].num_lattice_sites, phi_pt);
      prp_Y_double(prnY + i + k * nv/2 * level[0].num_lattice_sites, phi_pt);
      prp_X_double(prnX + i + k * nv/2 * level[0].num_lattice_sites, phi_pt);
    }
  }

  /// do block pi minus
  auto *prnTb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnZb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnYb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnXb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  pi_minus(block_phi, prnTb, prnZb, prnYb, prnXb, end*blen);

  /// transform block back and compare
  separate_vectors_in_block(prnTb, nv/2 * level[0].num_lattice_sites, nv/2, blen);
  separate_vectors_in_block(prnZb, nv/2 * level[0].num_lattice_sites, nv/2, blen);
  separate_vectors_in_block(prnYb, nv/2 * level[0].num_lattice_sites, nv/2, blen);
  separate_vectors_in_block(prnXb, nv/2 * level[0].num_lattice_sites, nv/2, blen);

  /// compare
  int offset = nv/2 * level[0].num_lattice_sites;
  check_relative_diff(prnTb, prnT, blen * offset, &level[0]);
  check_relative_diff(prnZb, prnZ, blen * offset, &level[0]);
  check_relative_diff(prnYb, prnY, blen * offset, &level[0]);
  check_relative_diff(prnXb, prnX, blen * offset, &level[0]);

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
  auto *phi1 = new std::complex<double>[blen*len];
  auto *phi2 = new std::complex<double>[blen*len];
  vector_double_define_random(phi1, 0, blen*len);
  vector_copy(phi2, phi1, 0, len*blen);

  auto *vs = new BlockVecPerSite<double>(phi1, nv);
  vs->transform(len, blen);
  auto *ve = new BlockVecPerEntry<double>(phi2);
  ve->transform(len, blen);

  auto *prnTb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnZb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnYb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnXb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  pi_minus(phi1, prnTb, prnZb, prnYb, prnXb, end*blen);
  /// wrap and de-transform
  auto *vsout_T = new BlockVecPerSite<double>(prnTb, nv/2);
  vsout_T->detransform(len/2, blen);
  auto *vsout_Z = new BlockVecPerSite<double>(prnZb, nv/2);
  vsout_Z->detransform(len/2, blen);
  auto *vsout_Y = new BlockVecPerSite<double>(prnYb, nv/2);
  vsout_Y->detransform(len/2, blen);
  auto *vsout_X = new BlockVecPerSite<double>(prnXb, nv/2);
  vsout_X->detransform(len/2, blen);

  auto *prnT = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnZ = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnY = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnX = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  pi_minus_per_entry(phi2, prnT, prnZ, prnY, prnX, end, blen);
  /// wrap and de-transform
  auto *veout_T = new BlockVecPerEntry<double>(prnT);
  veout_T->detransform(len/2, blen);
  auto *veout_Z = new BlockVecPerEntry<double>(prnZ);
  veout_Z->detransform(len/2, blen);
  auto *veout_Y = new BlockVecPerEntry<double>(prnY);
  veout_Y->detransform(len/2, blen);
  auto *veout_X = new BlockVecPerEntry<double>(prnX);
  veout_X->detransform(len/2, blen);

  check_relative_diff(vsout_T->v, veout_T->v, len/2*blen, &level[0]);
  check_relative_diff(vsout_Z->v, veout_Z->v, len/2*blen, &level[0]);
  check_relative_diff(vsout_Y->v, veout_Y->v, len/2*blen, &level[0]);
  check_relative_diff(vsout_X->v, veout_X->v, len/2*blen, &level[0]);

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
  auto* block_phi = new std::complex<double>[blen*len];
  auto* phi = new std::complex<double>[blen*len];
  vector_double_define_random(block_phi, 0, blen*len);
  for (int i = 0; i < blen*len; i++) {
    phi[i] = block_phi[i];
  }
  /// transform to block_phi
  transform_vector_block(block_phi, len, nv, blen);

  /// do normal pi plus product
  auto *prpT = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpZ = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpY = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpX = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  std::complex<double> *end_pt, *phi_pt;
  std::complex<double> pbuf[6];
  int j;
  int offset = nv/2 * level[0].num_lattice_sites;

  for (int k = 0; k < blen; ++k) {
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
  auto *prpTb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpZb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpYb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpXb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  pi_plus_kron_UH_psi(block_phi, level[0].gmres.matrix->D, prpTb, prpZb, prpYb, prpXb,
                      level[0].gmres.matrix->neighbor_table, end, blen);

  separate_vectors_in_block(prpTb, offset, nv/2, blen);
  separate_vectors_in_block(prpZb, offset, nv/2, blen);
  separate_vectors_in_block(prpYb, offset, nv/2, blen);
  separate_vectors_in_block(prpXb, offset, nv/2, blen);

  /// compare
  check_relative_diff(prpTb, prpT, blen * offset, &level[0]);
  check_relative_diff(prpZb, prpZ, blen * offset, &level[0]);
  check_relative_diff(prpYb, prpY, blen * offset, &level[0]);
  check_relative_diff(prpXb, prpX, blen * offset, &level[0]);

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
  auto *phi1 = new std::complex<double>[blen*len];
  auto *phi2 = new std::complex<double>[blen*len];
  vector_double_define_random(phi1, 0, blen*len);
  vector_copy(phi2, phi1, 0, len*blen);

  auto *vs = new BlockVecPerSite<double>(phi1, nv);
  vs->transform(len, blen);
  auto *ve = new BlockVecPerEntry<double>(phi2);
  ve->transform(len, blen);

  auto *prpTb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpZb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpYb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpXb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  pi_plus_kron_UH_psi(phi1, level[0].gmres.matrix->D, prpTb, prpZb, prpYb, prpXb,
                      level[0].gmres.matrix->neighbor_table, end, blen);
  /// wrap and de-transform
  auto *vsout_T = new BlockVecPerSite<double>(prpTb, nv/2);
  vsout_T->detransform(len/2, blen);
  auto *vsout_Z = new BlockVecPerSite<double>(prpZb, nv/2);
  vsout_Z->detransform(len/2, blen);
  auto *vsout_Y = new BlockVecPerSite<double>(prpYb, nv/2);
  vsout_Y->detransform(len/2, blen);
  auto *vsout_X = new BlockVecPerSite<double>(prpXb, nv/2);
  vsout_X->detransform(len/2, blen);

  auto *prpT = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpZ = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpY = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpX = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  pi_plus_per_entry(phi2, level[0].gmres.matrix->D, prpT, prpZ, prpY, prpX,
                    level[0].gmres.matrix->neighbor_table, end, blen);
  /// wrap and de-transform
  auto *veout_T = new BlockVecPerEntry<double>(prpT);
  veout_T->detransform(len/2, blen);
  auto *veout_Z = new BlockVecPerEntry<double>(prpZ);
  veout_Z->detransform(len/2, blen);
  auto *veout_Y = new BlockVecPerEntry<double>(prpY);
  veout_Y->detransform(len/2, blen);
  auto *veout_X = new BlockVecPerEntry<double>(prpX);
  veout_X->detransform(len/2, blen);

  check_relative_diff(vsout_T->v, veout_T->v, len/2*blen, &level[0]);
  check_relative_diff(vsout_Z->v, veout_Z->v, len/2*blen, &level[0]);
  check_relative_diff(vsout_Y->v, veout_Y->v, len/2*blen, &level[0]);
  check_relative_diff(vsout_X->v, veout_X->v, len/2*blen, &level[0]);

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

  auto *prnT = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnZ = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnY = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnX = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];

  auto *prnTb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnZb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnYb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnXb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];

  int offset = nv/2 * level[0].num_lattice_sites;

  /// initialize prn
  vector_double_define_random(prnT, 0, blen*offset);
  vector_double_define_random(prnZ, 0, blen*offset);
  vector_double_define_random(prnY, 0, blen*offset);
  vector_double_define_random(prnX, 0, blen*offset);
  for (int i = 0; i < blen*offset; i++) {
    prnTb[i] = prnT[i];
    prnZb[i] = prnZ[i];
    prnYb[i] = prnY[i];
    prnXb[i] = prnX[i];
  }
  transform_vector_block(prnTb, offset, nv/2, blen);
  transform_vector_block(prnZb, offset, nv/2, blen);
  transform_vector_block(prnYb, offset, nv/2, blen);
  transform_vector_block(prnXb, offset, nv/2, blen);

  /// apply U psi(x+mu) normally
  std::complex<double> *eta_pt, *end_pt;
  int j;
  std::complex<double> pbuf[6];
  auto *eta = new std::complex<double>[blen * len];
  for (int k = 0; k < blen; k++) {
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

  auto *etab = new std::complex<double>[blen * len];
  apply_U_psi_plus_mu(etab, level[0].gmres.matrix->D, prnTb, prnZb, prnYb, prnXb,
                      level[0].gmres.matrix->neighbor_table, end, blen);

  separate_vectors_in_block(etab, len, nv, blen);

  check_relative_diff(etab, eta, blen * len, &level[0]);

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
  auto *prnTb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnZb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnYb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnXb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnT = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnZ = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnY = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prnX = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];

  vector_double_define_random(prnTb, 0, blen*len/2);
  vector_double_define_random(prnZb, 0, blen*len/2);
  vector_double_define_random(prnYb, 0, blen*len/2);
  vector_double_define_random(prnXb, 0, blen*len/2);
  vector_copy(prnT, prnTb, 0, blen*len/2);
  vector_copy(prnZ, prnZb, 0, blen*len/2);
  vector_copy(prnY, prnYb, 0, blen*len/2);
  vector_copy(prnX, prnXb, 0, blen*len/2);

  auto *vsout_T = new BlockVecPerSite<double>(prnTb, nv/2);
  vsout_T->transform(len/2, blen);
  auto *vsout_Z = new BlockVecPerSite<double>(prnZb, nv/2);
  vsout_Z->transform(len/2, blen);
  auto *vsout_Y = new BlockVecPerSite<double>(prnYb, nv/2);
  vsout_Y->transform(len/2, blen);
  auto *vsout_X = new BlockVecPerSite<double>(prnXb, nv/2);
  vsout_X->transform(len/2, blen);

  /// initialize etas
  auto *etab = new std::complex<double>[blen * len];
  auto *eta = new std::complex<double>[blen*len];
  vector_double_define_random(etab, 0, blen*len);
  vector_copy(eta, etab, 0, blen*len);
  auto *vseta = new BlockVecPerSite<double>(etab, nv);
  auto *veeta = new BlockVecPerEntry<double>(eta);

  /// do per site
  vseta->transform(len, blen);
  apply_U_psi_plus_mu(etab, level[0].gmres.matrix->D, prnTb, prnZb, prnYb, prnXb,
                      level[0].gmres.matrix->neighbor_table, end, blen);
  vseta->detransform(len, blen);

  /// prepare prns for per entry
  auto *veout_T = new BlockVecPerEntry<double>(prnT);
  veout_T->transform(len/2, blen);
  auto *veout_Z = new BlockVecPerEntry<double>(prnZ);
  veout_Z->transform(len/2, blen);
  auto *veout_Y = new BlockVecPerEntry<double>(prnY);
  veout_Y->transform(len/2, blen);
  auto *veout_X = new BlockVecPerEntry<double>(prnX);
  veout_X->transform(len/2, blen);

  /// do per entry
  veeta->transform(len, blen);
  apply_U_per_entry(eta, level[0].gmres.matrix->D, prnT, prnZ, prnY, prnX,
                    level[0].gmres.matrix->neighbor_table, end, blen);
  veeta->detransform(len, blen);

  check_relative_diff(vseta->v, veeta->v, len*blen, &level[0]);

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

  auto *prpT = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpZ = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpY = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpX = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];

  auto *prpTb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpZb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpYb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpXb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];

  int offset = nv/2 * level[0].num_lattice_sites;

  /// initialize prp
  vector_double_define_random(prpT, 0, blen*offset);
  vector_double_define_random(prpZ, 0, blen*offset);
  vector_double_define_random(prpY, 0, blen*offset);
  vector_double_define_random(prpX, 0, blen*offset);
  for (int i = 0; i < blen*offset; i++) {
    prpTb[i] = prpT[i];
    prpZb[i] = prpZ[i];
    prpYb[i] = prpY[i];
    prpXb[i] = prpX[i];
  }
  transform_vector_block(prpTb, offset, nv/2, blen);
  transform_vector_block(prpZb, offset, nv/2, blen);
  transform_vector_block(prpYb, offset, nv/2, blen);
  transform_vector_block(prpXb, offset, nv/2, blen);

  /// do normal lift up
  auto *eta = new std::complex<double>[blen * len];
  for (int k = 0; k < blen; k++) {
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

  auto *etab = new std::complex<double>[blen * len];
  lift_pi_plus_dir(etab, prpTb, prpZb, prpYb, prpXb, end*blen);

  separate_vectors_in_block(etab, len, nv, blen);

  check_relative_diff(etab, eta, blen * len, &level[0]);

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
  auto *prpTb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpZb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpYb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpXb = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpT = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpZ = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpY = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];
  auto *prpX = new std::complex<double>[blen * nv/2 * level[0].num_lattice_sites];

  vector_double_define_random(prpTb, 0, blen*len/2);
  vector_double_define_random(prpZb, 0, blen*len/2);
  vector_double_define_random(prpYb, 0, blen*len/2);
  vector_double_define_random(prpXb, 0, blen*len/2);
  vector_copy(prpT, prpTb, 0, blen*len/2);
  vector_copy(prpZ, prpZb, 0, blen*len/2);
  vector_copy(prpY, prpYb, 0, blen*len/2);
  vector_copy(prpX, prpXb, 0, blen*len/2);

  auto *vsout_T = new BlockVecPerSite<double>(prpTb, nv/2);
  vsout_T->transform(len/2, blen);
  auto *vsout_Z = new BlockVecPerSite<double>(prpZb, nv/2);
  vsout_Z->transform(len/2, blen);
  auto *vsout_Y = new BlockVecPerSite<double>(prpYb, nv/2);
  vsout_Y->transform(len/2, blen);
  auto *vsout_X = new BlockVecPerSite<double>(prpXb, nv/2);
  vsout_X->transform(len/2, blen);

  /// initialize etas
  auto *etab = new std::complex<double>[blen * len];
  auto *eta = new std::complex<double>[blen*len];
  vector_double_define_random(etab, 0, blen*len);
  vector_copy(eta, etab, 0, blen*len);
  auto *vseta = new BlockVecPerSite<double>(etab, nv);
  auto *veeta = new BlockVecPerEntry<double>(eta);

  /// do per site
  vseta->transform(len, blen);
  lift_pi_plus_dir(etab, prpTb, prpZb, prpYb, prpXb, end*blen);
  vseta->detransform(len, blen);

  /// prepare prns for per entry
  auto *veout_T = new BlockVecPerEntry<double>(prpT);
  veout_T->transform(len/2, blen);
  auto *veout_Z = new BlockVecPerEntry<double>(prpZ);
  veout_Z->transform(len/2, blen);
  auto *veout_Y = new BlockVecPerEntry<double>(prpY);
  veout_Y->transform(len/2, blen);
  auto *veout_X = new BlockVecPerEntry<double>(prpX);
  veout_X->transform(len/2, blen);

  /// do per entry
  veeta->transform(len, blen);
  lift_pi_plus_per_entry(eta, prpT, prpZ, prpY, prpX, end, blen);
  veeta->detransform(len, blen);

  check_relative_diff(vseta->v, veeta->v, len*blen, &level[0]);

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
  const int measurement_round = 1;

  /// initialize phi and eta
  auto* phi2 = new std::complex<double>[blen*len];
  auto* eta2 = new std::complex<double>[blen*len];
  auto* block_phi = new std::complex<double>[blen*len];
  auto* block_eta = new std::complex<double>[blen*len];
  auto* phi = new std::complex<double>[blen*len];
  auto* eta = new std::complex<double>[blen*len];

  const char *eventNames[5];
  eventNames[0] = "PAPI_TOT_INS";
  eventNames[1] = "PAPI_TOT_CYC";


  for (int i = 0; i < measurement_round; ++i) {
    memset(eta, 0, blen * len * sizeof(std::complex<double>));
    memset(block_eta, 0, blen * len * sizeof(std::complex<double>));
    memset(eta2, 0, blen * len * sizeof(std::complex<double>));
    vector_double_define_random(block_phi, 0, len*blen);
    vector_double_copy(phi, block_phi, 0, len*blen, &level[0]);
    vector_copy(phi2, block_phi, 0, len*blen);

    /// preparing the objects for block D_W
    BlockVecPerSite<double> *veta[blen];
    BlockVecPerSite<double> *vphi[blen];
    auto *vbeta = new BlockVecPerSite<double>(block_eta, 12);
    auto *vbphi = new BlockVecPerSite<double>(block_phi, 12);
    auto *veeta = new BlockVecPerEntry<double>(eta2);
    auto *vephi = new BlockVecPerEntry<double>(phi2);

    for (int k = 0; k < blen; ++k) {
      veta[k] = new BlockVecPerSite<double>(&eta[k*len], 12);
      vphi[k] = new BlockVecPerSite<double>(&phi[k*len], 12);
    }

    t0 = MPI_Wtime();
    for (int k = 0; k < blen; k++) {
      d_plus_clover_block(veta[k], vphi[k], level[0].gmres.matrix, &level[0], 1);
    }
    t_nb += MPI_Wtime() - t0;

    //int eventSet = papi_start(eventNames, 2);
    vbphi->transform(len, blen);
    t0 = MPI_Wtime();
    d_plus_clover_block(vbeta, vbphi, level[0].gmres.matrix, &level[0], blen);
    t_b += MPI_Wtime() - t0;
    vbeta->detransform(len, blen);
    //papi_stop(eventSet);

    int eventSet2 = papi_start(eventNames, 2);
    vephi->transform(len, blen);
    t0 = MPI_Wtime();
    d_plus_clover_block(veeta, vephi, level[0].gmres.matrix, &level[0], blen);
    t_be += MPI_Wtime() - t0;
    veeta->detransform(len, blen);
    papi_stop(eventSet2);

    check_relative_diff(&vbeta->v[0], &veta[0]->v[0], blen * len, &level[0]);
    check_relative_diff(&vbeta->v[0], &veeta->v[0], blen * len, &level[0]);
  }

  printf0("No blocking for %d block length: %-8.4lf seconds, mesurement rounds: %d\n",
          blen, t_nb, measurement_round);
  printf0("blocking for %d block length: %-8.4lf seconds, mesurement rounds: %d\n",
          blen, t_b, measurement_round);
  printf0("blocking for %d block per entry length: %-8.4lf seconds, mesurement rounds: %d\n",
          blen, t_be, measurement_round);

  delete[] block_phi;
  delete[] block_eta;
  delete[] phi;
  delete[] eta;
}