#include "main.h"
#ifndef BLEN
#define BLEN 2
#endif

class OETest : public ::testing::Test {
protected:
  // You can remove any or all of the following functions if their bodies would
  // be empty.

  OETest() {
    // You can do set-up work for each test here.

  }

  ~OETest() override {
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

template<typename T>
static inline void check_relative_diff(std::complex<T>* v1, std::complex<T>* v2, int len,
                                       double prec, Level* l){
  vector_double_minus(v2, v2, v1, 0, len, l);
  double numerator = global_norm_double(v2, 0, len, l);
  double denominator = global_norm_double(v1, 0, len, l);
  ASSERT_LT(numerator/denominator, prec);
}

TEST_F(OETest, CheckInv){
  int vecLen = level[0].gmres.vectorLength;
  int vec_inner = level[0].gmres.vectorEnd;

  level[0].gmresSmoother.initialize(&level[0].oe_op_double, 20 , 20,
                                    vec_inner, 1e-8, _COARSE_GMRES, _NOTHING, 0, 0, vecLen,
                                    nullptr, apply_schur_complement_double, &level[0]);
  int start = level[0].inner_vector_size/2;
  int end = level[0].inner_vector_size;
  oddeven_setup_double(&level[0].global->op_double, &level[0]);

  auto* phi = new std::complex<double>[BLEN*vecLen];
  auto* eta = new std::complex<double>[BLEN*vecLen];
  auto* phi2 = new std::complex<double>[BLEN*vecLen];
  auto* eta2 = new std::complex<double>[BLEN*vecLen];

  auto* vephi = new BlockVecPerEntry<double>(phi2);
  auto* veeta = new BlockVecPerEntry<double>(eta2);

  vector_double_define_random(eta, 0, vecLen*BLEN);
  vector_copy(eta2, eta, 0, BLEN*vecLen);
  veeta->transform(vecLen, BLEN);

  for (int i = 0; i < BLEN; ++i) {
    //LLH_perform_fwd_bwd_subs_double(&phi[i*vecLen], &eta[i*vecLen], level[0].oe_op_double.clover);
    diag_oo_inv(&phi[i*vecLen], &eta[i*vecLen], &level[0].oe_op_double, &level[0], start, end);
  }
  //LLH_perform_fwd_bwd_subs_block(phi2, eta2, level[0].oe_op_double.clover, BLEN);
  diag_oo_inv_block(vephi, veeta, &level[0].oe_op_double, &level[0], start, end, BLEN);

  vephi->detransform(vecLen, BLEN);
  check_relative_diff(phi, phi2, BLEN*vecLen, 1e-15, &level[0]);

  delete[] phi;
  delete[] eta;
  delete[] phi2;
  delete[] eta2;
}

TEST_F(OETest, CheckSiteInv){
  int vecLen = level[0].gmres.vectorLength;
  int vec_inner = level[0].gmres.vectorEnd;

  level[0].gmresSmoother.initialize(&level[0].oe_op_double, 20 , 20,
                                    vec_inner, 1e-8, _COARSE_GMRES, _NOTHING, 0, 0, vecLen,
                                    nullptr, apply_schur_complement_double, &level[0]);
  int start = level[0].inner_vector_size/2;
  int end = level[0].inner_vector_size;
  oddeven_setup_double(&level[0].global->op_double, &level[0]);

  auto* phi = new std::complex<double>[BLEN*vecLen];
  auto* eta = new std::complex<double>[BLEN*vecLen];
  auto* phi2 = new std::complex<double>[BLEN*vecLen];
  auto* eta2 = new std::complex<double>[BLEN*vecLen];

  auto* vephi = new BlockVecPerEntry<double>(phi2);
  auto* veeta = new BlockVecPerEntry<double>(eta2);

  auto* vsphi = new BlockVecPerSite<double>(phi, 12);
  auto* vseta = new BlockVecPerSite<double>(eta, 12);

  vector_double_define_random(eta, 0, vecLen*BLEN);
  vector_copy(eta2, eta, 0, BLEN*vecLen);
  veeta->transform(vecLen, BLEN);
  vseta->transform(vecLen, BLEN);

  diag_oo_inv_block(vephi, veeta, &level[0].oe_op_double, &level[0], start, end, BLEN);

  diag_oo_inv_block(vsphi, vseta, &level[0].oe_op_double, &level[0], start, end, BLEN);

  vephi->detransform(vecLen, BLEN);
  vsphi->detransform(vecLen, BLEN);
  check_relative_diff(phi, phi2, BLEN*vecLen, 1e-15, &level[0]);

  delete[] phi;
  delete[] eta;
  delete[] phi2;
  delete[] eta2;
}

TEST_F(OETest, CheckClover){
  int vecLen = level[0].gmres.vectorLength;
  int vec_inner = level[0].gmres.vectorEnd;

  level[0].gmresSmoother.initialize(&level[0].oe_op_double, 20 , 20,
                                    vec_inner, 1e-8, _COARSE_GMRES, _NOTHING, 0, 0, vecLen,
                                    nullptr, apply_schur_complement_double, &level[0]);
  int start = 0;
  int end = level[0].inner_vector_size/2;
  oddeven_setup_double(&level[0].global->op_double, &level[0]);

  auto* phi = new std::complex<double>[BLEN*vecLen];
  auto* eta = new std::complex<double>[BLEN*vecLen];
  auto* phi2 = new std::complex<double>[BLEN*vecLen];
  auto* eta2 = new std::complex<double>[BLEN*vecLen];

  auto* vephi = new BlockVecPerEntry<double>(phi2);
  auto* veeta = new BlockVecPerEntry<double>(eta2);

  vector_double_define_random(eta, 0, vecLen*BLEN);
  vector_copy(eta2, eta, 0, BLEN*vecLen);
  veeta->transform(vecLen, BLEN);

  for (int i = 0; i < BLEN; ++i) {
    //LLH_multiply(&phi[i*vecLen], &eta[i*vecLen], level[0].oe_op_double.clover);
    diag_ee(&phi[i*vecLen], &eta[i*vecLen], &level[0].oe_op_double, &level[0], start, end);
  }
  //LLH_multiply_block(phi2, eta2, level[0].oe_op_double.clover, BLEN);
  diag_ee_block(vephi, veeta, &level[0].oe_op_double, &level[0], start, end, BLEN);

  vephi->detransform(vecLen, BLEN);
  check_relative_diff(phi, phi2, BLEN*vecLen, 1e-15, &level[0]);

  delete[] phi;
  delete[] eta;
  delete[] phi2;
  delete[] eta2;
}

TEST_F(OETest, CheckSiteClover){
  int vecLen = level[0].gmres.vectorLength;
  int vec_inner = level[0].gmres.vectorEnd;

  level[0].gmresSmoother.initialize(&level[0].oe_op_double, 20 , 20,
                                    vec_inner, 1e-8, _COARSE_GMRES, _NOTHING, 0, 0, vecLen,
                                    nullptr, apply_schur_complement_double, &level[0]);
  int start = 0;
  int end = level[0].inner_vector_size/2;
  oddeven_setup_double(&level[0].global->op_double, &level[0]);

  auto* phi = new std::complex<double>[BLEN*vecLen];
  auto* eta = new std::complex<double>[BLEN*vecLen];
  auto* phi2 = new std::complex<double>[BLEN*vecLen];
  auto* eta2 = new std::complex<double>[BLEN*vecLen];

  auto* vephi = new BlockVecPerEntry<double>(phi2);
  auto* veeta = new BlockVecPerEntry<double>(eta2);

  auto* vsphi = new BlockVecPerSite<double>(phi, 12);
  auto* vseta = new BlockVecPerSite<double>(eta, 12);

  vector_double_define_random(eta, 0, vecLen*BLEN);
  vector_copy(eta2, eta, 0, BLEN*vecLen);
  veeta->transform(vecLen, BLEN);
  vseta->transform(vecLen, BLEN);

  diag_ee_block(vephi, veeta, &level[0].oe_op_double, &level[0], start, end, BLEN);

  diag_ee_block(vsphi, vseta, &level[0].oe_op_double, &level[0], start, end, BLEN);

  vephi->detransform(vecLen, BLEN);
  vsphi->detransform(vecLen, BLEN);
  check_relative_diff(phi, phi2, BLEN*vecLen, 1e-15, &level[0]);

  delete[] phi;
  delete[] eta;
  delete[] phi2;
  delete[] eta2;
}

TEST_F(OETest, CheckHoppingTerm){
  int vecLen = level[0].gmres.vectorLength;
  int vec_inner = level[0].gmres.vectorEnd;
  int amount = _EVEN_SITES;
  auto* phi = new std::complex<double>[vecLen * BLEN];
  auto* phi2 = new std::complex<double>[vecLen * BLEN];
  auto* eta = new std::complex<double>[vecLen * BLEN];
  auto* eta2 = new std::complex<double>[vecLen * BLEN];

  level[0].gmresSmoother.initialize(&level[0].oe_op_double, 20 , 20,
                                    vec_inner, 1e-8, _COARSE_GMRES, _NOTHING, 0, 0, vecLen,
                                    nullptr, apply_schur_complement_double, &level[0]);
  oddeven_setup_double(&level[0].global->op_double, &level[0]);

  vector_double_define_random(phi, 0, vecLen * BLEN);
  vector_copy(phi2, phi, 0, vecLen * BLEN);
  auto* vephi = new BlockVecPerEntry<double>(phi2);
  auto* veeta = new BlockVecPerEntry<double>(eta2);

  for (int i = 0; i < BLEN; i++) {
    hopping_term_double(&eta[i*vecLen], &phi[i*vecLen], &level[0].oe_op_double, amount, &level[0]);
  }

  vephi->transform(vecLen, BLEN);
  hopping_term_block(veeta, vephi, &level[0].oe_op_double, amount, &level[0], BLEN);
  veeta->detransform(vecLen, BLEN);

  check_relative_diff(eta, eta2, vecLen*BLEN, 1e-15, &level[0]);

  delete[] phi;
  delete[] phi2;
  delete[] eta;
  delete[] eta2;
}

TEST_F(OETest, CheckSiteHoppingTerm){
  int vecLen = level[0].gmres.vectorLength;
  int vec_inner = level[0].gmres.vectorEnd;
  int amount = _EVEN_SITES;
  auto* phi = new std::complex<double>[vecLen * BLEN];
  auto* phi2 = new std::complex<double>[vecLen * BLEN];
  auto* eta = new std::complex<double>[vecLen * BLEN];
  auto* eta2 = new std::complex<double>[vecLen * BLEN];

  level[0].gmresSmoother.initialize(&level[0].oe_op_double, 20 , 20,
                                    vec_inner, 1e-8, _COARSE_GMRES, _NOTHING, 0, 0, vecLen,
                                    nullptr, apply_schur_complement_double, &level[0]);
  oddeven_setup_double(&level[0].global->op_double, &level[0]);

  vector_double_define_random(phi, 0, vecLen * BLEN);
  vector_copy(phi2, phi, 0, vecLen * BLEN);
  auto* vephi = new BlockVecPerEntry<double>(phi2);
  auto* veeta = new BlockVecPerEntry<double>(eta2);

  auto* vsphi = new BlockVecPerSite<double>(phi, 12);
  auto* vseta = new BlockVecPerSite<double>(eta, 12);

  vephi->transform(vecLen, BLEN);
  hopping_term_block(veeta, vephi, &level[0].oe_op_double, amount, &level[0], BLEN);
  veeta->detransform(vecLen, BLEN);

  vsphi->transform(vecLen, BLEN);
  hopping_term_block(vseta, vsphi, &level[0].oe_op_double, amount, &level[0], BLEN);
  vseta->detransform(vecLen, BLEN);

  check_relative_diff(eta, eta2, vecLen*BLEN, 1e-15, &level[0]);

  delete[] phi;
  delete[] phi2;
  delete[] eta;
  delete[] eta2;
}

TEST_F(OETest, CheckSchurComplement){
  int vecLen = level[0].gmres.vectorLength;
  int vec_inner = level[0].gmres.vectorEnd;

  auto *btmp = new std::complex<double>[vecLen*BLEN];
  auto *brhs = new std::complex<double>[vecLen*BLEN];
  //auto *bx = new std::complex<double>[vecLen*BLEN];

  auto *tmp = new std::complex<double>[vecLen*BLEN];
  auto *rhs = new std::complex<double>[vecLen*BLEN];
  //auto *x = new std::complex<double>[vecLen*BLEN];

  level[0].gmresSmoother.initialize(&level[0].oe_op_double, 20 , 20,
                                    vec_inner, 1e-8, _COARSE_GMRES, _NOTHING, 0, 0, vecLen,
                                    nullptr, apply_schur_complement_double, &level[0]);
  level[0].gmresSmoother.verbose = 1;
  int start = level[0].inner_vector_size/2;
  int end = level[0].inner_vector_size;
  int half = start;
  oddeven_setup_double(&level[0].global->op_double, &level[0]);

  vector_double_define_random(brhs, 0, vecLen * BLEN);

  vector_copy(rhs, brhs, 0, vecLen * BLEN);

  for (int i = 0; i < BLEN; ++i) {
    /// tmp = Doo^-1 eta_o
    diag_oo_inv(&tmp[i*vecLen], &rhs[i*vecLen], &level[0].oe_op_double, &level[0], start, end);
    /// tmp = -Doo^-1 eta_o
    vector_double_scale(&tmp[i*vecLen], &tmp[i*vecLen], -1, start, end, &level[0]);
    /// b_e = eta = eta_e - Deo Doo^-1 eta_o
    hopping_term_double(&rhs[i*vecLen], &tmp[i*vecLen], &level[0].oe_op_double, _EVEN_SITES, &level[0]);

    //apply_schur_complement_double(&x[i*vecLen], &rhs[i*vecLen], &level[0].oe_op_double, &level[0]);
    vector_copy(level[0].gmresSmoother.b, &rhs[i*vecLen], 0, vecLen);
    level[0].gmresSmoother.solve();
  }

  auto *vtmp = new BlockVecPerEntry<double>(btmp);
  auto *vrhs = new BlockVecPerEntry<double>(brhs);
  //auto *vx = new BlockVecPerEntry<double>(bx);
  vrhs->transform(vecLen, BLEN);

  /// set up b for odd even
  /// tmp = Doo^-1 eta_o
  diag_oo_inv_block(vtmp, vrhs, &level[0].oe_op_double, &level[0], start, end, BLEN);
  /// tmp = -Doo^-1 eta_o
  vector_scale(btmp, btmp, std::complex<double>(-1), start, end, BLEN);

  /// b_e = eta = eta_e - Deo Doo^-1 eta_o
  hopping_term_block(vrhs, vtmp, &level[0].oe_op_double, _EVEN_SITES, &level[0], BLEN);

  //apply_schur_complement_block(vx, vrhs, &level[0].oe_op_double, &level[0], BLEN);

  auto *bgmres = new BlockGmres<double>(level[0].oe_op_double.num_even_sites * 12, BLEN, level[0].gmresSmoother.matrix, &level[0],
                                        1e-8, 20, 20, true, true);
  void (*op)(BlockVecPerEntry<double>*, BlockVecPerEntry<double>*,
             operator_struct<double>*, Level*, const int) = &apply_schur_complement_block;
  vector_copy(bgmres->b, brhs, 0, half*BLEN);
  bgmres->solve(op);

  //auto *vr = new BlockVecPerEntry<double>(bgmres->r);
  //vr->detransform(vecLen, BLEN);
  //check_relative_diff(rtmp, rhs, half, 1e-15, &level[0]);
  //check_relative_diff(vr->v, rtmp, half, 1e-15, &level[0]);

  //vrhs->detransform(vecLen, BLEN);
  //check_relative_diff(rhs, brhs, vecLen*BLEN, 1e-15, &level[0]);
  //vx->detransform(vecLen, BLEN);
  //check_relative_diff(x, bx, vecLen*BLEN, 1e-15, &level[0]);

  delete[] brhs;
  delete[] btmp;
  //delete[] bx;
  //delete[] x;
  delete[] rhs;
  delete[] tmp;
}

TEST_F(OETest, CheckSiteSchurComplement){
  int vecLen = level[0].gmres.vectorLength;
  int vec_inner = level[0].gmres.vectorEnd;

  auto *btmp = new std::complex<double>[vecLen*BLEN];
  auto *brhs = new std::complex<double>[vecLen*BLEN];
  //auto *bx = new std::complex<double>[vecLen*BLEN];

  auto *tmp = new std::complex<double>[vecLen*BLEN];
  auto *rhs = new std::complex<double>[vecLen*BLEN];
  //auto *x = new std::complex<double>[vecLen*BLEN];

  level[0].gmresSmoother.initialize(&level[0].oe_op_double, 20 , 20,
                                    vec_inner, 1e-8, _COARSE_GMRES, _NOTHING, 0, 0, vecLen,
                                    nullptr, apply_schur_complement_double, &level[0]);
  level[0].gmresSmoother.verbose = 1;
  int start = level[0].inner_vector_size/2;
  int end = level[0].inner_vector_size;
  int half = start;
  oddeven_setup_double(&level[0].global->op_double, &level[0]);

  vector_double_define_random(brhs, 0, vecLen * BLEN);

  vector_copy(rhs, brhs, 0, vecLen * BLEN);

  for (int i = 0; i < BLEN; ++i) {
    /// tmp = Doo^-1 eta_o
    diag_oo_inv(&tmp[i*vecLen], &rhs[i*vecLen], &level[0].oe_op_double, &level[0], start, end);
    /// tmp = -Doo^-1 eta_o
    vector_double_scale(&tmp[i*vecLen], &tmp[i*vecLen], -1, start, end, &level[0]);
    /// b_e = eta = eta_e - Deo Doo^-1 eta_o
    hopping_term_double(&rhs[i*vecLen], &tmp[i*vecLen], &level[0].oe_op_double, _EVEN_SITES, &level[0]);

    //apply_schur_complement_double(&x[i*vecLen], &rhs[i*vecLen], &level[0].oe_op_double, &level[0]);
    vector_copy(level[0].gmresSmoother.b, &rhs[i*vecLen], 0, vecLen);
    level[0].gmresSmoother.solve();
  }

  auto *vtmp = new BlockVecPerSite<double>(btmp, 12);
  auto *vrhs = new BlockVecPerSite<double>(brhs, 12);
  //auto *vx = new BlockVecPerEntry<double>(bx);
  vrhs->transform(vecLen, BLEN);

  /// set up b for odd even
  /// tmp = Doo^-1 eta_o
  diag_oo_inv_block(vtmp, vrhs, &level[0].oe_op_double, &level[0], start, end, BLEN);
  /// tmp = -Doo^-1 eta_o
  vector_scale(btmp, btmp, std::complex<double>(-1), start, end, BLEN);

  /// b_e = eta = eta_e - Deo Doo^-1 eta_o
  hopping_term_block(vrhs, vtmp, &level[0].oe_op_double, _EVEN_SITES, &level[0], BLEN);

  //apply_schur_complement_block(vx, vrhs, &level[0].oe_op_double, &level[0], BLEN);

  auto *bgmres = new BlockGmres<double>(level[0].oe_op_double.num_even_sites * 12, BLEN, level[0].gmresSmoother.matrix, &level[0],
                                        1e-8, 20, 20, true, true);
  void (*op)(BlockVecPerSite<double>*, BlockVecPerSite<double>*,
             operator_struct<double>*, Level*, const int) = &apply_schur_complement_block;
  vector_copy(bgmres->b, brhs, 0, half*BLEN);
  bgmres->solve(op);

  //auto *vr = new BlockVecPerEntry<double>(bgmres->r);
  //vr->detransform(vecLen, BLEN);
  //check_relative_diff(rtmp, rhs, half, 1e-15, &level[0]);
  //check_relative_diff(vr->v, rtmp, half, 1e-15, &level[0]);

  //vrhs->detransform(vecLen, BLEN);
  //check_relative_diff(rhs, brhs, vecLen*BLEN, 1e-15, &level[0]);
  //vx->detransform(vecLen, BLEN);
  //check_relative_diff(x, bx, vecLen*BLEN, 1e-15, &level[0]);

  delete[] brhs;
  delete[] btmp;
  //delete[] bx;
  //delete[] x;
  delete[] rhs;
  delete[] tmp;
}

TEST_F(OETest, CheckOddEven){
  int vecLen = level[0].gmres.vectorLength;
  int vec_inner = level[0].gmres.vectorEnd;
  auto* tmp = new std::complex<double>[vecLen];

  level[0].gmresSmoother.initialize(&level[0].oe_op_double, 20 , 20,
                                    vec_inner, 1e-8, _COARSE_GMRES, _NOTHING, 0, 0, vecLen,
                                    nullptr, apply_schur_complement_double, &level[0]);
  int start = level[0].inner_vector_size/2;
  int end = level[0].inner_vector_size;
  oddeven_setup_double(&level[0].global->op_double, &level[0]);

  vector_double_define_random(level[0].gmres.b, 0, vecLen);
  //vector_copy(level[0].gmresSmoother.b, level[0].gmres.b, 0, vecLen);
  serial_to_oddeven_double(level[0].gmresSmoother.b, level[0].gmres.b, &level[0]);

  /// tmp = Doo^-1 eta_o
  diag_oo_inv(&tmp[0], level[0].gmresSmoother.b, &level[0].oe_op_double, &level[0], start, end);
  /// tmp = -Doo^-1 eta_o
  vector_double_scale(tmp, tmp, -1, start, end, &level[0]);
  /// b_e = eta = eta_e - Deo Doo^-1 eta_o
  hopping_term_double(level[0].gmresSmoother.b, tmp, &level[0].oe_op_double, _EVEN_SITES, &level[0]);

  level[0].gmresSmoother.verbose = 1;
  level[0].gmresSmoother.solve();

  /// psi_o = Doo^-1 eta_o
  diag_oo_inv(level[0].gmresSmoother.x, level[0].gmresSmoother.b, &level[0].oe_op_double, &level[0], start, end);
  /// tmp = Doe psi_e
  vector_double_define(tmp, 0, start, end);
  hopping_term_double(tmp, level[0].gmresSmoother.x, &level[0].oe_op_double, _ODD_SITES, &level[0]);
  /// b_o = Doo^-1 Doe psi_e
  diag_oo_inv(level[0].gmresSmoother.b, tmp, &level[0].oe_op_double, &level[0], start, end);
  /// psi_o = psi_o - b_o = Doo^-1 eta_o - Doo^-1 Doe psi_e
  vector_double_minus(level[0].gmresSmoother.x, level[0].gmresSmoother.x, level[0].gmresSmoother.b, start, end, &level[0]);

  oddeven_to_serial_double(tmp, level[0].gmresSmoother.x, &level[0]);
  level[0].gmres.solve();

  check_relative_diff(tmp, level[0].gmres.x, vecLen, 1e-7, &level[0]);
  delete[] tmp;
}