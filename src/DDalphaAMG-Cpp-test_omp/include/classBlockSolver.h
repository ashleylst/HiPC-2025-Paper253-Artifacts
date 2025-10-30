/*
* Copyright (C) 2024, Shiting Long.
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
 */
#ifndef DDALPHAAMG_CLASSBLOCKSOLVER_H
#define DDALPHAAMG_CLASSBLOCKSOLVER_H

#include "block.h"
#include "dirac_block_vpsite.h"

template<typename T>
class BlockSolver{

public:
  BlockSolver(int problemSize, int blen, operator_struct<T> *op, Level *l,
              T tol, int rn, int rl, bool InitGuessZero, bool v){
    this->vectorEnd = problemSize;
    //this->vectorLength = vecLen;
    this->blockLength = blen;
    this->matrix = op;
    this->level = l;
    //this->applyOperator = evalOp;
    this->restartNum = rn;
    this->restartLen = rl;
    this->initialGuessIsZero = InitGuessZero;
    this->verbose = v;
    this->tolerance = tol;

    x = new std::complex<T>[blockLength * vectorEnd];
    b = new std::complex<T>[blockLength * vectorEnd];
    r = new std::complex<T>[blockLength * vectorEnd];
  }; // default constructor

  BlockSolver(const BlockSolver &) = delete;            // copy constructor
  BlockSolver &operator=(const BlockSolver &) = delete; // assignment operator

  BlockSolver(const BlockSolver &&) = delete;            // move constructor
  BlockSolver &operator=(const BlockSolver &&) = delete; // move assignment operator

  ~BlockSolver(){
    delete[] x;
    delete[] b;
    delete[] r;
  }; // default destructor

  //void (*applyOperator)(std::complex<T>*, std::complex<T>*, operator_struct<T>*,
                     //   Level*, const int);

  T tolerance;
  int vectorEnd;
  //int vectorLength;
  int blockLength;
  int restartLen;
  int restartNum;
  bool verbose;
  bool finished = false;
  bool initialGuessIsZero = true;

  /// result
  std::complex<T> *b; /// rhs
  std::complex<T> *r; /// residual

  operator_struct<T> *matrix; /// Dirac operator info
  Level *level;
  std::complex<T> *x;
};

template<typename T>
class BlockGmres : public BlockSolver<T>{
public:
  BlockGmres<T>(int probSize, int blen, operator_struct<T>* op, Level *l,
                T tol, int rn, int rl, bool init, bool v)
      : BlockSolver<T>(probSize, blen, op, l, tol, rn, rl, init, v){
    //this->applyOperator = evalOp;
    this->res = 0;
    this->elapsedTime = 0.0;
  };

  BlockGmres(const BlockGmres &) = delete;            // copy constructor
  BlockGmres &operator=(const BlockGmres &) = delete; // assignment operator

  BlockGmres(const BlockGmres &&) = delete;            // move constructor
  BlockGmres &operator=(const BlockGmres &&) = delete; // move assignment operator

  ~BlockGmres()= default;

  void solve(void (*evalOp)(BlockVecPerSite<T>*, BlockVecPerSite<T>*,
      operator_struct<T>*, Level*, const int));

  void solve(void (*evalOp)(BlockVecPerEntry<T>*, BlockVecPerEntry<T>*,
                            operator_struct<T>*, Level*, const int));

  void arnoldiStep(complex<T> **V, complex<T> *H, complex<T> *w, complex<T> *y, int curLoop,
                   void (*evalOp)(BlockVecPerSite<T>*, BlockVecPerSite<T>*,
                       operator_struct<T>*, Level*, const int));

  void arnoldiStep(complex<T> **V, complex<T> *H, complex<T> *w, complex<T> *y, int curLoop,
                   void (*evalOp)(BlockVecPerEntry<T>*, BlockVecPerEntry<T>*,
                                  operator_struct<T>*, Level*, const int));

  void qrUpdate(std::complex<T> *H, std::complex<T> *s, std::complex<T> *c, std::complex<T> *gamma,
                int curLoop);

  void computeSolution(std::complex<T> **V, std::complex<T> *y,std::complex<T> *gamma,
                       std::complex<T> *H, bool *marked, int curLoop, int outerLoop, int layout);

private:
  int res;
  double elapsedTime;
  //void (*applyOperator)(std::complex<T>*, std::complex<T>*, operator_struct<T>*,
                        //Level*, const int);
};



#endif // DDALPHAAMG_CLASSBLOCKSOLVER_H
