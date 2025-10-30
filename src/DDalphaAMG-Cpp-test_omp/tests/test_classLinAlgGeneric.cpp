#include <main.h>


// The fixture for testing Linalggeneric.
class LinAlgGenericTest : public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if their bodies would
  // be empty.
  LinAlgGenericTest() {

    // Vectors: 1, -1, i, -i
    vector_double_define(level[0].vbuf_double[8],  1.0, 0, level[0].inner_vector_size);
    vector_double_define(level[0].vbuf_double[7], -1.0, 0, level[0].inner_vector_size);
    vector_double_define(level[0].vbuf_double[6],    I, 0, level[0].inner_vector_size);
    vector_double_define(level[0].vbuf_double[5],   -I, 0, level[0].inner_vector_size);
  }

  ~LinAlgGenericTest() override {
     // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  void SetUp() override {
     // Code here will be called immediately after the constructor (right
     // before each test).
     level = EnvironmentDDalphaAMG::level;
     norm_vector_one = sqrt(level[0].inner_vector_size*level[0].inner_vector_size);
  }

  void TearDown() override {
     // Code here will be called immediately after each test (right
     // before the destructor).
  }

  // Class members declared here can be used by all tests in the test suite
  Level *level = EnvironmentDDalphaAMG::level;
  complex_double product;
  double norm_vector_one;
};


//Syntax: Suite
TEST_F(LinAlgGenericTest, InnerProductTest) {


    // Product of one vector i with vector i.
    product = global_inner_product_double(level[0].vbuf_double[6], level[0].vbuf_double[6], 0, level[0].inner_vector_size, &level[0]);
    ASSERT_DOUBLE_EQ(real(product), norm_vector_one);
    ASSERT_DOUBLE_EQ(imag(product), 0.0);

   // Product of one vector i with vector -i.
    product = global_inner_product_double(level[0].vbuf_double[6], level[0].vbuf_double[5], 0, level[0].inner_vector_size, &level[0]);
    ASSERT_DOUBLE_EQ(real(product), -norm_vector_one);
    ASSERT_DOUBLE_EQ(imag(product), 0.0);

    // Product of one vector 1 with vector 1.
    product = global_inner_product_double(level[0].vbuf_double[8], level[0].vbuf_double[8], 0, level[0].inner_vector_size, &level[0]);
    ASSERT_DOUBLE_EQ(real(product), norm_vector_one);
    ASSERT_DOUBLE_EQ(imag(product), 0.0);

    // Product of one vector 1 with vector -1.
    product = global_inner_product_double(level[0].vbuf_double[8], level[0].vbuf_double[8], 0, level[0].inner_vector_size, &level[0]);
    ASSERT_DOUBLE_EQ(real(product), norm_vector_one);
    ASSERT_DOUBLE_EQ(imag(product), 0.0);

   // Product of one vector i with vector 1.
    product = global_inner_product_double(level[0].vbuf_double[6], level[0].vbuf_double[8], 0, level[0].inner_vector_size, &level[0]);
    ASSERT_DOUBLE_EQ(real(product), 0.0);
    ASSERT_DOUBLE_EQ(imag(product), -norm_vector_one);


}


