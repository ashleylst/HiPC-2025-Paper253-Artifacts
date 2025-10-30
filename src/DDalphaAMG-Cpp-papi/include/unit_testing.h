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

#ifndef UNIT_TESTING_HEADER
#define UNIT_TESTING_HEADER

#include "main.h"
#include "stdio.h"

using namespace std;

class EnvironmentDDalphaAMG : public ::testing::Environment {
public:
  // Assume there's only going to be a single instance of this class

  static Level *level;
  static Global global;

  // remove braces to use SetUp() in .cpp file
  void SetUp(); // setting up env...

  // remove braces to use TearDown() in .cpp file
  void TearDown();
};

#endif
