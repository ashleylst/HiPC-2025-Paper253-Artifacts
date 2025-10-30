#include <main.h>

// The fixture for testing class Level.
class BlockTest : public ::testing::Test {
protected:
  // You can remove any or all of the following functions if their bodies would
  // be empty.

  BlockTest() {
    // You can do set-up work for each test here.

  }

  ~BlockTest() override {
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

const int blen = 2;

template<typename T>
static inline void check_relative_diff(std::complex<T>* v1, std::complex<T>* v2, int len,
                                       double prec, Level* l){
  vector_double_minus(v2, v2, v1, 0, len, l);
  double numerator = global_norm_double(v2, 0, len, l);
  double denominator = global_norm_double(v1, 0, len, l);
  ASSERT_LT(numerator/denominator, prec);
}

template<typename T>
void transform_H(std::complex<T> *H, int veclen, int blen){
  auto *tmp = new std::complex<T>[blen*veclen];
  for(int i = 0; i < veclen; i++){
    for (int j = 0; j < blen; j++) {
      tmp[i*blen + j] = H[j*veclen + i];
    }
  }
  for (int i = 0; i < blen*veclen; ++i) {
    H[i] = tmp[i];
  }
  delete[] tmp;
}

template<typename T>
void detransform_H(std::complex<T> *H, int veclen, int blen){
  auto *tmp = new std::complex<T>[blen*veclen];
  for(int i = 0; i < veclen; i++){
    for (int j = 0; j < blen; j++) {
      tmp[j*veclen + i] = H[i*blen + j];
    }
  }
  for (int i = 0; i < blen*veclen; ++i) {
    H[i] = tmp[i];
  }
  delete[] tmp;
}



TEST_F(BlockTest, CheckQR){
  int restartLen = level[0].gmres.restartLength;
  int restartNum = level[0].gmres.restartNumber;
  int problemSize = level[0].gmres.vectorEnd;
  int Hsize = ((restartLen+1)*(restartLen+1) + 3*(restartLen+1))/2;
  int testnum = 19;

  auto *H = new complex_double[blen * Hsize];
  auto *gamma = new complex_double[blen * (restartLen+1)];
  auto *s = new complex_double[blen * restartLen];
  auto *c = new complex_double[blen * restartLen];

  auto *bH = new complex_double[blen * Hsize];
  auto *bgamma = new complex_double[blen * (restartLen+1)];
  auto *bs = new complex_double[blen * restartLen];
  auto *bc = new complex_double[blen * restartLen];

  vector_double_define_random(bH, 0, blen*Hsize);
  vector_double_copy(H, bH, 0, blen*Hsize, &level[0]);
  vector_double_define_random(bgamma, testnum*blen, (testnum+1)*blen);
  vector_double_define_random(bc, 0, testnum*blen);
  vector_double_define_random(bs, 0, testnum*blen);

  for (int i = 0; i < blen; ++i) {
    level[0].gmres.j = testnum;
    level[0].gmres.gamma[testnum] = bgamma[testnum*blen + i];
    for (int j = 0; j < testnum; ++j) {
      level[0].gmres.cosinus[j] = bc[j*blen + i];
      level[0].gmres.sinus[j] = bs[j*blen + i];
    }

    level[0].gmres.qrUpdate(&H[i * Hsize]);

    vector_double_copy_3(&s[i*restartLen], level[0].gmres.sinus, 0, restartLen, &level[0]);
    vector_double_copy_3(&c[i*restartLen], level[0].gmres.cosinus, 0, restartLen, &level[0]);
    vector_double_copy_3(&gamma[i*(restartLen+1)], level[0].gmres.gamma, 0, restartLen+1, &level[0]);
  }

  transform_H(bH, Hsize, blen);
  auto* b = new BlockGmres<double>(problemSize, blen, level[0].gmres.matrix, &level[0],
                                    1e-9, restartNum, restartLen, true, false);
  b->qrUpdate(bH, bs, bc, bgamma, testnum);
  detransform_H(bH, Hsize, blen);
  detransform_H(bs, restartLen, blen);
  detransform_H(bc, restartLen, blen);
  detransform_H(bgamma, restartLen+1, blen);

  check_relative_diff(bH, H, blen * Hsize, 1e-12, &level[0]);
  check_relative_diff(bs, s, blen * restartLen, 1e-12, &level[0]);
  check_relative_diff(bc, c, blen * restartLen, 1e-12, &level[0]);
  check_relative_diff(bgamma, gamma, blen * (restartLen+1), 1e-12, &level[0]);

  delete[] H;
  delete[] gamma;
  delete[] s;
  delete[] c;
  delete[] bH;
  delete[] bgamma;
  delete[] bs;
  delete[] bc;
}

TEST_F(BlockTest, CheckComputeSol){
  int vecLen = level[0].gmres.vectorEnd;
  int restartLen = level[0].gmres.restartLength;
  int restartNum = level[0].gmres.restartNumber;
  int Hsize = ((restartLen+1)*(restartLen+1) + 3*(restartLen+1))/2;
  int testnum = 19;
  int outerloop = 1;

  auto *H = new complex_double[blen * Hsize];
  auto *y = new std::complex<double>[blen * (restartLen+1)];
  auto *gamma = new complex_double[restartLen+1];
  auto *x = new std::complex<double>[blen * vecLen];

  auto *by = new std::complex<double>[blen * (restartLen+1)];
  auto *bH = new complex_double[blen * Hsize];
  auto *bgamma = new complex_double[blen * (restartLen+1)];

  std::complex<double> **bV;
  NEW2D(bV, complex_double, blen * vecLen, restartLen + 1);
  std::complex<double> **V;
  NEW2D(V, complex_double, vecLen, restartLen + 1);

  vector_double_define_random(bH, 0, blen*Hsize);
  vector_double_copy(H, bH, 0, blen*Hsize, &level[0]);
  vector_double_define_random(bgamma, 0, blen*(restartLen+1));

  for (int i = 0; i < restartLen + 1; ++i) {
    vector_double_define_random(bV[i], 0, blen*vecLen);
  }

  DELETE2D(level[0].gmres.V);
  delete[] level[0].gmres.gamma;
  level[0].gmres.gamma = gamma;
  level[0].gmres.V = V;
  for (int i = 0; i < blen; i++) {
    for (int j = 0; j < restartLen + 1; ++j) {
      vector_copy(V[j], bV[j], i * vecLen, (i + 1) * vecLen);
    }
    level[0].gmres.j = testnum;
    memset(level[0].gmres.x, 0, vecLen*sizeof(std::complex<double>));
    vector_copy(gamma, bgamma, i * (restartLen + 1), (i + 1) * (restartLen + 1));
    level[0].gmres.computeSolution(outerloop, &H[i * Hsize]);

    vector_double_copy_3(&y[i*(restartLen+1)], level[0].gmres.y, 0, restartLen+1, &level[0]);
    vector_double_copy_3(&x[i*vecLen], level[0].gmres.x, 0, vecLen, &level[0]);
  }

  bool *marked = new bool[blen];
  memset(marked, 0, blen*sizeof(marked[0]));
  transform_H(bH, Hsize, blen);
  transform_H(bgamma, restartLen+1, blen);
  for(int i = 0; i < restartLen+1; i++){
    transform_vector_block(bV[i], vecLen, 12, blen);
  }
  auto* b = new BlockGmres<double>(vecLen, blen, level[0].gmres.matrix, &level[0],
                                   1e-9, restartNum, restartLen, true, false);
  b->computeSolution(bV, by, bgamma, bH, marked, testnum, outerloop);
  detransform_H(by, restartLen+1, blen);
  separate_vectors_in_block(b->x, vecLen, 12, blen);

  check_relative_diff(by, y, blen*(restartLen+1), 1e-12, &level[0]);
  check_relative_diff(b->x, x, blen*vecLen, 1e-12,  &level[0]);

  delete[] H;
  delete[] x;
  delete[] by;
  delete[] bH;
  delete[] bgamma;
  DELETE2D(bV);
}

TEST_F(BlockTest, CheckArnoldi){
  int vecLen = level[0].gmres.vectorEnd;
  int restartLen = level[0].gmres.restartLength;
  int restartNum = level[0].gmres.restartNumber;
  int Hsize = ((restartLen+1)*(restartLen+1) + 3*(restartLen+1))/2;
  int testnum = 18;

  /// V has blen*vecLen columns, restartLen+1 rows -> it has restartLen+1 vectors
  std::complex<double> **bV;
  NEW2D(bV, complex_double, blen * vecLen, restartLen + 1);
  std::complex<double> **V;
  NEW2D(V, complex_double, vecLen, restartLen + 1);
  auto *v = new std::complex<double>[blen * vecLen];

  auto *by = new std::complex<double>[blen * (restartLen+1)];
  auto *bw = new std::complex<double>[blen * vecLen];
  auto *w = new std::complex<double>[blen * vecLen];

  /// H is stored H[0]_b1 H[0]_b2 ... H[1]_b1 H[1]_b2
  auto *H = new complex_double[blen * Hsize];
  auto *bH = new complex_double[blen * Hsize];

  for (int i = 0; i < restartLen + 1; ++i) {
    vector_double_define_random(bV[i], 0, blen*vecLen);
  }

  DELETE2D(level[0].gmres.V);
  level[0].gmres.V = V;
  //level[0].gmres.j = 0;
  for (int i = 0; i < blen; ++i) {
    for (int j = 0; j < restartLen + 1; ++j) {
      vector_copy(V[j], bV[j], i * vecLen, (i + 1) * vecLen);
    }

    for (int k = 0; k < 2; ++k) {
      level[0].gmres.j = testnum + k;
      level[0].gmres.arnoldiStep(&H[i * Hsize]);
      vector_double_copy_3(&w[i*vecLen], level[0].gmres.w, 0, vecLen, &level[0]);
      vector_double_copy_3(&v[i*vecLen], level[0].gmres.V[testnum+1], 0, vecLen, &level[0]);
    }
  }

  for(int i = 0; i < restartLen+1; i++){
    transform_vector_block(bV[i], vecLen, 12, blen);
  }
  auto* b = new BlockGmres<double>(vecLen, blen, level[0].gmres.matrix, &level[0],
                                   1e-9, restartNum, restartLen, true, false);
  for (int k = 0; k < 2; ++k) {
    b->arnoldiStep(bV, bH, bw, by, testnum + k, d_plus_clover_block);
  }

  separate_vectors_in_block(bw, vecLen, 12, blen);
  check_relative_diff(bw, w, blen*vecLen, 1e-12, &level[0]);

  for (int j = 0; j < blen; ++j) {
    for (int k = 0; k < 2; ++k) {
      for (int i = 0; i < testnum + k + 2; ++i) {
        complex_double a1 = bH[get_H_idx(testnum+k, i, j, blen)] -
                            H[j * Hsize + ((testnum+k) * (testnum+k) + 3 * (testnum+k)) / 2 + i];
        complex_double a2 = H[j * Hsize + ((testnum+k) * (testnum+k) + 3 * (testnum+k)) / 2 + i];

        ASSERT_LT(real(a1) / real(a2), 1e-11) << "i " << i << " k " << k << std::endl;
        if (i != testnum + k + 1) {
          ASSERT_LT(imag(a1) / imag(a2), 1e-11) << "i " << i << std::endl;
        }
      }
    }
  }

  separate_vectors_in_block(bV[testnum+1], vecLen, 12, blen);
  check_relative_diff(bV[testnum+1], v, blen*vecLen, 1e-12, &level[0]);

  DELETE2D(bV);
  delete[] by;
  delete[] bw;
  delete[] bH;
  delete[] H;
  delete[] w;
}

TEST_F(BlockTest, CheckInnerProduct){
  int vecLen = level[0].gmres.vectorEnd;
  int restartLen = level[0].gmres.restartLength;

  std::complex<double> **bV;
  NEW2D(bV, complex_double, blen * vecLen, restartLen + 1);
  std::complex<double> **V;
  NEW2D(V, complex_double, vecLen, restartLen + 1);

  auto *by = new std::complex<double>[blen * (restartLen+1)];
  auto *bw = new std::complex<double>[blen * vecLen];
  auto *y = new std::complex<double>[blen * (restartLen+1)];
  auto *w = new std::complex<double>[vecLen];

  /// initialize V and w
  for (int i = 0; i < restartLen + 1; ++i) {
    vector_double_define_random(bV[i], 0, blen*vecLen);
  }
  vector_double_define_random(bw, 0, blen * vecLen);

  for (int i = 0; i < blen; ++i) {
    for (int j = 0; j < restartLen+1; ++j) {
      vector_copy(V[j], bV[j], i * vecLen, (i + 1) * vecLen);
    }
    vector_copy(w, bw, i * vecLen, (i + 1) * vecLen);

    /// normal multi dot product y[i] = <V_H[i], w>
    process_multi_inner_product_double(restartLen + 1, &y[i*(restartLen+1)], V, w, 0, vecLen, &level[0]);
  }

  for(int i = 0; i < restartLen+1; i++){
    transform_vector_block(bV[i], vecLen, 12, blen);
  }
  transform_vector_block(bw, vecLen, 12, blen);
  /// block multi dot product
  multi_inner_product(restartLen+1, by, bV, bw, vecLen, blen, &level[0]);

  /// compare
  for (int i = 0; i < restartLen+1; ++i) {
    for (int j = 0; j < blen; ++j) {
      complex_double a = by[blen*i + j] - y[j*(restartLen+1) + i];
      complex_double b = y[j*(restartLen+1) + i];

      ASSERT_LT(real(a)/ real(b), 1e-11) << "i " << i << "j " << j << std::endl;
      ASSERT_LT(imag(a)/ imag(b), 1e-11) << "i " << i << "j " << j << std::endl;
    }
  }


  DELETE2D(bV);
  DELETE2D(V);
  delete[] by;
  delete[] bw;
  delete[] y;
  delete[] w;
}

TEST_F(BlockTest, CheckSolve){
  //level[0].gmres.solve();

  int vecLen = level[0].gmres.vectorEnd;
  int restartLen = level[0].gmres.restartLength;
  int restartNum = level[0].gmres.restartNumber;
  int Hsize = ((restartLen+1)*(restartLen+1) + 3*(restartLen+1))/2;

  auto *rhs = new std::complex<double>[blen * vecLen];

  //auto *norm = new std::complex<double>[blen];
  //auto *tmp = new std::complex<double>[blen * vecLen];
  //auto *inject = new std::complex<double>[blen * vecLen];

  vector_double_define_random(rhs, 0, blen*vecLen);

  level[0].gmres.restartLength = restartLen;
  level[0].gmres.restartNumber = restartNum;
  level[0].gmres.tolerance = 1e-8;
  for (int i = 0; i < blen; ++i) {
    memset(level[0].gmres.x, 0, vecLen*sizeof(std::complex<double>));
    memset(level[0].gmres.r, 0, vecLen*sizeof(std::complex<double>));
    vector_copy(level[0].gmres.b, rhs, i*vecLen, (i+1)*vecLen);
    level[0].gmres.solve();
    //norm[i] = level[0].gmres.norm;
    //vector_double_copy_3(&tmp[i*vecLen], level[0].gmres.x, 0, vecLen, &level[0]);
  }

  auto *bgmres = new BlockGmres<double>(vecLen, blen, level[0].gmres.matrix, &level[0],
                                        1e-8, restartNum, restartLen, true, true);
  vector_copy(bgmres->b, rhs, 0, blen*vecLen);
  void (*op)(BlockVecPerSite<double>*, BlockVecPerSite<double>*,
      operator_struct<double>*, Level*, const int) = &d_plus_clover_block;
  bgmres->solve(op);

  //detransform_H(inject, Hsize, blen);
  //check_relative_diff(bgmres->x, tmp,  blen*vecLen, 1e-11, &level[0]);
  /*separate_vectors_in_block(inject, vecLen, 12, blen);
  for (int i = 0; i < blen*vecLen; ++i) {
    ASSERT_LT(real(inject[i])- real(tmp[i]), 1e-11) << "i " << i << std::endl;
    ASSERT_LT(imag(inject[i])- imag(tmp[i]), 1e-11) << "i " << i << std::endl;
  }
  check_relative_diff(inject, tmp,  blen*vecLen, 1e-11, &level[0]);*/


  delete[] rhs;
  //delete[] norm;
  //delete[] tmp;
  //delete[] inject;
}