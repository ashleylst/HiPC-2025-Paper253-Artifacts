#include <main.h>

// The fixture for testing class Level.
class LevelTest : public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if their bodies would
  // be empty.

  LevelTest() {
     // You can do set-up work for each test here.

  }

  ~LevelTest() override {
     // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  void SetUp() override {
     // Code here will be called immediately after the constructor (right
     // before each test).

    foo_double=1.0, bar_double=2.0;


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

  double foo_double, bar_double;
};


//Syntax: Suite
TEST_F(LevelTest, RunSolve) {

    level[0].gmres.applyOperator(level[0].gmres.x, level[0].gmres.b, level[0].gmres.matrix, &level[0]);
    double testnorm = global_norm_double(level[0].gmres.x, 0, level[0].inner_vector_size, &level[0]);
    printf0("matvec test: %lf\n", testnorm);

   /* vector_double_define(level[1].s_double.b, 1, 0, level[1].inner_vector_size, &level[1]);
    apply_coarse_operator_double(level[1].s_double.x, level[1].s_double.b, &level[1].s_double.op, &level[1]);
    testnorm = global_norm_double(level[1].s_double.x, 0, level[1].inner_vector_size, &level[1]);
    printf0("matvec test coarse: %lf\n", testnorm);*/

    level[0].gmres.solve();
}


TEST_F(LevelTest, CheckResidual) {
//1. Generate rhs vector (b) done in set-up
//2. Solve 
  level[0].gmres.solve();
//3. Apply operator on solution (Dx)
  level[0].gmres.applyOperator(level[0].gmres.w, level[0].gmres.x, level[0].gmres.matrix, &level[0]);

//4. compute residual(Dx - b)
  int start = 0, end = level[0].inner_vector_size; //TODO: fix once threading is working

  vector_double_minus(level[0].gmres.w, level[0].gmres.w, level[0].gmres.b, start, end, &level[0]);

//5. check relative residual	(||Ax -b||/||b|| < global_tol)
  double w_norm = global_norm_double(level[0].gmres.w, 0, level[0].inner_vector_size, &level[0]);
  double b_norm = global_norm_double(level[0].gmres.b, 0, level[0].inner_vector_size, &level[0]);

  ASSERT_LT((w_norm)/b_norm, EnvironmentDDalphaAMG::global.tol);
}
