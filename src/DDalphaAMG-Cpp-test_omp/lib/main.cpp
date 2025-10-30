/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern
 * Leder.
 *
 * This file is part of the DDalphaAMG solver library.
 *
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 *
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 *
 */

#include "main.h"

int test_argc;
char **test_argv;

Level *EnvironmentDDalphaAMG::level;
Global EnvironmentDDalphaAMG::global;

int main(int argc, char **argv) {

  test_argc = argc;
  test_argv = argv;

  /*
   * Here are the commands for the filters one can apply to  run ONLY a test suite or a test
   *    ::testing::GTEST_FLAG(filter) = "FooSuite.*"   Runs everything in test suite FooSuite
   *    ::testing::GTEST_FLAG(filter) = "FooSuite.FooTest"  Runs test FooTest in suite FooTest
   *    ::testing::GTEST_FLAG(filter) = "~FooSuite.FooTest"  Skips test FooTest in suite FooTest
   * If not specified, all tests will be run.
   *
   */
  ::testing::GTEST_FLAG(filter) = "DiracTest.*";

  ::testing::InitGoogleTest(&argc, argv);

  /*
   * Definition of the Gtest environment
   * its setUp Phase is called before all tests are run
   */
  auto my_env = new EnvironmentDDalphaAMG;
  ::testing::AddGlobalTestEnvironment(my_env);

  return RUN_ALL_TESTS();
  /* The tearDown phase of the environment is called after
   * all tests are run*/
}
