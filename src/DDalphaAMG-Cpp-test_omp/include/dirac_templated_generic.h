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

#ifndef DDALPHAAMG_DIRAC_TEMPLATED_GENERIC_H
#define DDALPHAAMG_DIRAC_TEMPLATED_GENERIC_H

template<typename T>
static inline void prp_T(std::complex<T> *prp_pt, std::complex<T> *l_pt) {
  /// for each line, complex - complex plus complex * complex plus real * real,
  /// floating point operation is therefore 2+6+1
  prp_pt[0] = l_pt[0] - (std::complex<T>)GAMMA_T_SPIN0_VAL * l_pt[3 * GAMMA_T_SPIN0_CO];
  prp_pt[1] = l_pt[1] - (std::complex<T>)GAMMA_T_SPIN0_VAL * l_pt[3 * GAMMA_T_SPIN0_CO + 1];
  prp_pt[2] = l_pt[2] - (std::complex<T>)GAMMA_T_SPIN0_VAL * l_pt[3 * GAMMA_T_SPIN0_CO + 2];
  prp_pt[3] = l_pt[3] - (std::complex<T>)GAMMA_T_SPIN1_VAL * l_pt[3 * GAMMA_T_SPIN1_CO];
  prp_pt[4] = l_pt[4] - (std::complex<T>)GAMMA_T_SPIN1_VAL * l_pt[3 * GAMMA_T_SPIN1_CO + 1];
  prp_pt[5] = l_pt[5] - (std::complex<T>)GAMMA_T_SPIN1_VAL * l_pt[3 * GAMMA_T_SPIN1_CO + 2];
}

template<typename T>
static inline void prp_T_block(std::complex<T> *__restrict__ prp, std::complex<T> *__restrict__ phi, const int blen){
  for (int j = 0; j < blen; j++) {
    prp[j] = phi[j] - (std::complex<T>)GAMMA_T_SPIN0_VAL * phi[3*GAMMA_T_SPIN0_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[blen + j] = phi[blen + j] - (std::complex<T>)GAMMA_T_SPIN0_VAL *
                                        phi[(3*GAMMA_T_SPIN0_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[2*blen + j] = phi[2*blen + j] - (std::complex<T>)GAMMA_T_SPIN0_VAL *
                                                phi[(3*GAMMA_T_SPIN0_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    prp[3*blen + j] = phi[3*blen + j] - (std::complex<T>)GAMMA_T_SPIN1_VAL *
                                                phi[3*GAMMA_T_SPIN1_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[4*blen + j] = phi[4*blen + j] - (std::complex<T>)GAMMA_T_SPIN1_VAL *
                                                phi[(3*GAMMA_T_SPIN1_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[5*blen + j] = phi[5*blen + j] - (std::complex<T>)GAMMA_T_SPIN1_VAL *
                                                phi[(3*GAMMA_T_SPIN1_CO+2)*blen + j];
  }
}

template<typename T>
static inline void prp_Z(std::complex<T> *prp_pt, std::complex<T> *l_pt) {
  prp_pt[0] = l_pt[0] - (std::complex<T>)GAMMA_Z_SPIN0_VAL * l_pt[3 * GAMMA_Z_SPIN0_CO];
  prp_pt[1] = l_pt[1] - (std::complex<T>)GAMMA_Z_SPIN0_VAL * l_pt[3 * GAMMA_Z_SPIN0_CO + 1];
  prp_pt[2] = l_pt[2] - (std::complex<T>)GAMMA_Z_SPIN0_VAL * l_pt[3 * GAMMA_Z_SPIN0_CO + 2];
  prp_pt[3] = l_pt[3] - (std::complex<T>)GAMMA_Z_SPIN1_VAL * l_pt[3 * GAMMA_Z_SPIN1_CO];
  prp_pt[4] = l_pt[4] - (std::complex<T>)GAMMA_Z_SPIN1_VAL * l_pt[3 * GAMMA_Z_SPIN1_CO + 1];
  prp_pt[5] = l_pt[5] - (std::complex<T>)GAMMA_Z_SPIN1_VAL * l_pt[3 * GAMMA_Z_SPIN1_CO + 2];
}

template<typename T>
static inline void prp_Z_block(std::complex<T> *__restrict__ prp, std::complex<T> *__restrict__ phi, const int blen){
  for (int j = 0; j < blen; j++) {
    prp[j] = phi[j] - (std::complex<T>)GAMMA_Z_SPIN0_VAL * phi[3*GAMMA_Z_SPIN0_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[blen + j] = phi[blen + j] - (std::complex<T>)GAMMA_Z_SPIN0_VAL *
                                                phi[(3*GAMMA_Z_SPIN0_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[2*blen + j] = phi[2*blen + j] - (std::complex<T>)GAMMA_Z_SPIN0_VAL *
                                                phi[(3*GAMMA_Z_SPIN0_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    prp[3*blen + j] = phi[3*blen + j] - (std::complex<T>)GAMMA_Z_SPIN1_VAL *
                                                phi[3*GAMMA_Z_SPIN1_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[4*blen + j] = phi[4*blen + j] - (std::complex<T>)GAMMA_Z_SPIN1_VAL *
                                                phi[(3*GAMMA_Z_SPIN1_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[5*blen + j] = phi[5*blen + j] - (std::complex<T>)GAMMA_Z_SPIN1_VAL *
                                                phi[(3*GAMMA_Z_SPIN1_CO+2)*blen + j];
  }
}

template<typename T>
static inline void prp_Y(std::complex<T> *prp_pt, std::complex<T> *l_pt) {
  prp_pt[0] = l_pt[0] - (std::complex<T>)GAMMA_Y_SPIN0_VAL * l_pt[3 * GAMMA_Y_SPIN0_CO];
  prp_pt[1] = l_pt[1] - (std::complex<T>)GAMMA_Y_SPIN0_VAL * l_pt[3 * GAMMA_Y_SPIN0_CO + 1];
  prp_pt[2] = l_pt[2] - (std::complex<T>)GAMMA_Y_SPIN0_VAL * l_pt[3 * GAMMA_Y_SPIN0_CO + 2];
  prp_pt[3] = l_pt[3] - (std::complex<T>)GAMMA_Y_SPIN1_VAL * l_pt[3 * GAMMA_Y_SPIN1_CO];
  prp_pt[4] = l_pt[4] - (std::complex<T>)GAMMA_Y_SPIN1_VAL * l_pt[3 * GAMMA_Y_SPIN1_CO + 1];
  prp_pt[5] = l_pt[5] - (std::complex<T>)GAMMA_Y_SPIN1_VAL * l_pt[3 * GAMMA_Y_SPIN1_CO + 2];
}

template<typename T>
static inline void prp_Y_block(std::complex<T> *__restrict__ prp, std::complex<T> *__restrict__ phi, const int blen){
  for (int j = 0; j < blen; j++) {
    prp[j] = phi[j] - (std::complex<T>)GAMMA_Y_SPIN0_VAL * phi[3*GAMMA_Y_SPIN0_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[blen + j] = phi[blen + j] - (std::complex<T>)GAMMA_Y_SPIN0_VAL *
                                                phi[(3*GAMMA_Y_SPIN0_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[2*blen + j] = phi[2*blen + j] - (std::complex<T>)GAMMA_Y_SPIN0_VAL *
                                                phi[(3*GAMMA_Y_SPIN0_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    prp[3*blen + j] = phi[3*blen + j] - (std::complex<T>)GAMMA_Y_SPIN1_VAL *
                                                phi[3*GAMMA_Y_SPIN1_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[4*blen + j] = phi[4*blen + j] - (std::complex<T>)GAMMA_Y_SPIN1_VAL *
                                                phi[(3*GAMMA_Y_SPIN1_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[5*blen + j] = phi[5*blen + j] - (std::complex<T>)GAMMA_Y_SPIN1_VAL *
                                                phi[(3*GAMMA_Y_SPIN1_CO+2)*blen + j];
  }
}

template<typename T>
static inline void prp_X(std::complex<T> *prp_pt, std::complex<T> *l_pt) {
  prp_pt[0] = l_pt[0] - (std::complex<T>)GAMMA_X_SPIN0_VAL * l_pt[3 * GAMMA_X_SPIN0_CO];
  prp_pt[1] = l_pt[1] - (std::complex<T>)GAMMA_X_SPIN0_VAL * l_pt[3 * GAMMA_X_SPIN0_CO + 1];
  prp_pt[2] = l_pt[2] - (std::complex<T>)GAMMA_X_SPIN0_VAL * l_pt[3 * GAMMA_X_SPIN0_CO + 2];
  prp_pt[3] = l_pt[3] - (std::complex<T>)GAMMA_X_SPIN1_VAL * l_pt[3 * GAMMA_X_SPIN1_CO];
  prp_pt[4] = l_pt[4] - (std::complex<T>)GAMMA_X_SPIN1_VAL * l_pt[3 * GAMMA_X_SPIN1_CO + 1];
  prp_pt[5] = l_pt[5] - (std::complex<T>)GAMMA_X_SPIN1_VAL * l_pt[3 * GAMMA_X_SPIN1_CO + 2];
}

template<typename T>
static inline void prp_X_block(std::complex<T> *__restrict__ prp, std::complex<T> *__restrict__ phi, const int blen){
  for (int j = 0; j < blen; j++) {
    prp[j] = phi[j] - (std::complex<T>)GAMMA_X_SPIN0_VAL * phi[3*GAMMA_X_SPIN0_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[blen + j] = phi[blen + j] - (std::complex<T>)GAMMA_X_SPIN0_VAL *
                                                phi[(3*GAMMA_X_SPIN0_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[2*blen + j] = phi[2*blen + j] - (std::complex<T>)GAMMA_X_SPIN0_VAL *
                                                phi[(3*GAMMA_X_SPIN0_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    prp[3*blen + j] = phi[3*blen + j] - (std::complex<T>)GAMMA_X_SPIN1_VAL *
                                                phi[3*GAMMA_X_SPIN1_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[4*blen + j] = phi[4*blen + j] - (std::complex<T>)GAMMA_X_SPIN1_VAL *
                                                phi[(3*GAMMA_X_SPIN1_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prp[5*blen + j] = phi[5*blen + j] - (std::complex<T>)GAMMA_X_SPIN1_VAL *
                                                phi[(3*GAMMA_X_SPIN1_CO+2)*blen + j];
  }
}

template<typename T>
static inline void prn_T(std::complex<T> *prn_pt, std::complex<T> *l_pt) {
  prn_pt[0] = l_pt[0] + (std::complex<T>)GAMMA_T_SPIN0_VAL * l_pt[3 * GAMMA_T_SPIN0_CO];
  prn_pt[1] = l_pt[1] + (std::complex<T>)GAMMA_T_SPIN0_VAL * l_pt[3 * GAMMA_T_SPIN0_CO + 1];
  prn_pt[2] = l_pt[2] + (std::complex<T>)GAMMA_T_SPIN0_VAL * l_pt[3 * GAMMA_T_SPIN0_CO + 2];
  prn_pt[3] = l_pt[3] + (std::complex<T>)GAMMA_T_SPIN1_VAL * l_pt[3 * GAMMA_T_SPIN1_CO];
  prn_pt[4] = l_pt[4] + (std::complex<T>)GAMMA_T_SPIN1_VAL * l_pt[3 * GAMMA_T_SPIN1_CO + 1];
  prn_pt[5] = l_pt[5] + (std::complex<T>)GAMMA_T_SPIN1_VAL * l_pt[3 * GAMMA_T_SPIN1_CO + 2];
}

template<typename T>
static inline void prn_T_block(std::complex<T> *__restrict__ prn, std::complex<T> *__restrict__ phi, const int blen){
  for (int j = 0; j < blen; j++) {
    prn[j] = phi[j] + (std::complex<T>)GAMMA_T_SPIN0_VAL * phi[3*GAMMA_T_SPIN0_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[blen + j] = phi[blen + j] + (std::complex<T>)GAMMA_T_SPIN0_VAL *
                                                phi[(3*GAMMA_T_SPIN0_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[2*blen + j] = phi[2*blen + j] + (std::complex<T>)GAMMA_T_SPIN0_VAL *
                                                phi[(3*GAMMA_T_SPIN0_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    prn[3*blen + j] = phi[3*blen + j] + (std::complex<T>)GAMMA_T_SPIN1_VAL *
                                                phi[3*GAMMA_T_SPIN1_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[4*blen + j] = phi[4*blen + j] + (std::complex<T>)GAMMA_T_SPIN1_VAL *
                                                phi[(3*GAMMA_T_SPIN1_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[5*blen + j] = phi[5*blen + j] + (std::complex<T>)GAMMA_T_SPIN1_VAL *
                                                phi[(3*GAMMA_T_SPIN1_CO+2)*blen + j];
  }
}

template<typename T>
static inline void prn_Z(std::complex<T> *prn_pt, std::complex<T> *l_pt) {
  prn_pt[0] = l_pt[0] + (std::complex<T>)GAMMA_Z_SPIN0_VAL * l_pt[3 * GAMMA_Z_SPIN0_CO];
  prn_pt[1] = l_pt[1] + (std::complex<T>)GAMMA_Z_SPIN0_VAL * l_pt[3 * GAMMA_Z_SPIN0_CO + 1];
  prn_pt[2] = l_pt[2] + (std::complex<T>)GAMMA_Z_SPIN0_VAL * l_pt[3 * GAMMA_Z_SPIN0_CO + 2];
  prn_pt[3] = l_pt[3] + (std::complex<T>)GAMMA_Z_SPIN1_VAL * l_pt[3 * GAMMA_Z_SPIN1_CO];
  prn_pt[4] = l_pt[4] + (std::complex<T>)GAMMA_Z_SPIN1_VAL * l_pt[3 * GAMMA_Z_SPIN1_CO + 1];
  prn_pt[5] = l_pt[5] + (std::complex<T>)GAMMA_Z_SPIN1_VAL * l_pt[3 * GAMMA_Z_SPIN1_CO + 2];
}

template<typename T>
static inline void prn_Z_block(std::complex<T> *__restrict__ prn, std::complex<T> *__restrict__ phi, const int blen){
  for (int j = 0; j < blen; j++) {
    prn[j] = phi[j] + (std::complex<T>)GAMMA_Z_SPIN0_VAL * phi[3*GAMMA_Z_SPIN0_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[blen + j] = phi[blen + j] + (std::complex<T>)GAMMA_Z_SPIN0_VAL *
                                                phi[(3*GAMMA_Z_SPIN0_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[2*blen + j] = phi[2*blen + j] + (std::complex<T>)GAMMA_Z_SPIN0_VAL *
                                                phi[(3*GAMMA_Z_SPIN0_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    prn[3*blen + j] = phi[3*blen + j] + (std::complex<T>)GAMMA_Z_SPIN1_VAL *
                                                phi[3*GAMMA_Z_SPIN1_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[4*blen + j] = phi[4*blen + j] + (std::complex<T>)GAMMA_Z_SPIN1_VAL *
                                                phi[(3*GAMMA_Z_SPIN1_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[5*blen + j] = phi[5*blen + j] + (std::complex<T>)GAMMA_Z_SPIN1_VAL *
                                                phi[(3*GAMMA_Z_SPIN1_CO+2)*blen + j];
  }
}

template<typename T>
static inline void prn_Y(std::complex<T> *prn_pt, std::complex<T> *l_pt) {
  prn_pt[0] = l_pt[0] + (std::complex<T>)GAMMA_Y_SPIN0_VAL * l_pt[3 * GAMMA_Y_SPIN0_CO];
  prn_pt[1] = l_pt[1] + (std::complex<T>)GAMMA_Y_SPIN0_VAL * l_pt[3 * GAMMA_Y_SPIN0_CO + 1];
  prn_pt[2] = l_pt[2] + (std::complex<T>)GAMMA_Y_SPIN0_VAL * l_pt[3 * GAMMA_Y_SPIN0_CO + 2];
  prn_pt[3] = l_pt[3] + (std::complex<T>)GAMMA_Y_SPIN1_VAL * l_pt[3 * GAMMA_Y_SPIN1_CO];
  prn_pt[4] = l_pt[4] + (std::complex<T>)GAMMA_Y_SPIN1_VAL * l_pt[3 * GAMMA_Y_SPIN1_CO + 1];
  prn_pt[5] = l_pt[5] + (std::complex<T>)GAMMA_Y_SPIN1_VAL * l_pt[3 * GAMMA_Y_SPIN1_CO + 2];
}

template<typename T>
static inline void prn_Y_block(std::complex<T> *__restrict__ prn, std::complex<T> *__restrict__ phi, const int blen){
  for (int j = 0; j < blen; j++) {
    prn[j] = phi[j] + (std::complex<T>)GAMMA_Y_SPIN0_VAL * phi[3*GAMMA_Y_SPIN0_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[blen + j] = phi[blen + j] + (std::complex<T>)GAMMA_Y_SPIN0_VAL *
                                                phi[(3*GAMMA_Y_SPIN0_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[2*blen + j] = phi[2*blen + j] + (std::complex<T>)GAMMA_Y_SPIN0_VAL *
                                                phi[(3*GAMMA_Y_SPIN0_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    prn[3*blen + j] = phi[3*blen + j] + (std::complex<T>)GAMMA_Y_SPIN1_VAL *
                                                phi[3*GAMMA_Y_SPIN1_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[4*blen + j] = phi[4*blen + j] + (std::complex<T>)GAMMA_Y_SPIN1_VAL *
                                                phi[(3*GAMMA_Y_SPIN1_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[5*blen + j] = phi[5*blen + j] + (std::complex<T>)GAMMA_Y_SPIN1_VAL *
                                                phi[(3*GAMMA_Y_SPIN1_CO+2)*blen + j];
  }
}

template<typename T>
static inline void prn_X(std::complex<T> *prn_pt, std::complex<T> *l_pt) {
  prn_pt[0] = l_pt[0] + (std::complex<T>)GAMMA_X_SPIN0_VAL * l_pt[3 * GAMMA_X_SPIN0_CO];
  prn_pt[1] = l_pt[1] + (std::complex<T>)GAMMA_X_SPIN0_VAL * l_pt[3 * GAMMA_X_SPIN0_CO + 1];
  prn_pt[2] = l_pt[2] + (std::complex<T>)GAMMA_X_SPIN0_VAL * l_pt[3 * GAMMA_X_SPIN0_CO + 2];
  prn_pt[3] = l_pt[3] + (std::complex<T>)GAMMA_X_SPIN1_VAL * l_pt[3 * GAMMA_X_SPIN1_CO];
  prn_pt[4] = l_pt[4] + (std::complex<T>)GAMMA_X_SPIN1_VAL * l_pt[3 * GAMMA_X_SPIN1_CO + 1];
  prn_pt[5] = l_pt[5] + (std::complex<T>)GAMMA_X_SPIN1_VAL * l_pt[3 * GAMMA_X_SPIN1_CO + 2];
}

template<typename T>
static inline void prn_X_block(std::complex<T> *__restrict__ prn, std::complex<T> *__restrict__ phi, const int blen){
  for (int j = 0; j < blen; j++) {
    prn[j] = phi[j] + (std::complex<T>)GAMMA_X_SPIN0_VAL * phi[3*GAMMA_X_SPIN0_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[blen + j] = phi[blen + j] + (std::complex<T>)GAMMA_X_SPIN0_VAL *
                                                phi[(3*GAMMA_X_SPIN0_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[2*blen + j] = phi[2*blen + j] + (std::complex<T>)GAMMA_X_SPIN0_VAL *
                                                phi[(3*GAMMA_X_SPIN0_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    prn[3*blen + j] = phi[3*blen + j] + (std::complex<T>)GAMMA_X_SPIN1_VAL *
                                                phi[3*GAMMA_X_SPIN1_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[4*blen + j] = phi[4*blen + j] + (std::complex<T>)GAMMA_X_SPIN1_VAL *
                                                phi[(3*GAMMA_X_SPIN1_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    prn[5*blen + j] = phi[5*blen + j] + (std::complex<T>)GAMMA_X_SPIN1_VAL *
                                                phi[(3*GAMMA_X_SPIN1_CO+2)*blen + j];
  }
}

template<typename T>
static inline void mvmh(std::complex<T> *eta, std::complex<T> *D,
                               std::complex<T> *phi) {
  /// for each line below, we ignore the conj operation because it's flipping a sign,
  /// we do complex * complex, floating point operation is 6 (2*2+2)
  eta[0] = conj(D[0]) * phi[0];
  eta[1] = conj(D[1]) * phi[0];
  eta[2] = conj(D[2]) * phi[0];
  /// for each line below, we do complex * complex and then complex + complex, FLOP is 6+2
  eta[0] += conj(D[3]) * phi[1];
  eta[1] += conj(D[4]) * phi[1];
  eta[2] += conj(D[5]) * phi[1];
  eta[0] += conj(D[6]) * phi[2];
  eta[1] += conj(D[7]) * phi[2];
  eta[2] += conj(D[8]) * phi[2];
}

template<typename T>
static inline void mvmh_block(std::complex<T> *__restrict__ eta, std::complex<T> *__restrict__ D, std::complex<T> *__restrict__ phi,
                              const int blen){
  for (int j = 0; j < blen; j++) {
    eta[j] = conj(D[0]) * phi[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] = conj(D[1]) * phi[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] = conj(D[2]) * phi[j];
  }

  for (int j = 0; j < blen; j++) {
    eta[j] += conj(D[3]) * phi[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] += conj(D[4]) * phi[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] += conj(D[5]) * phi[blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[j] += conj(D[6]) * phi[2*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] += conj(D[7]) * phi[2*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] += conj(D[8]) * phi[2*blen + j];
  }
}

template<typename T>
static inline void mvm(std::complex<T> *eta, std::complex<T> *D,
                              std::complex<T> *phi) {
  /// for each block below, we do complex * complex 3 times, complex + complex 2 times
  /// floating point operation is 3*6+2*2 = 22
  eta[0] = D[0] * phi[0];
  eta[0] += D[1] * phi[1];
  eta[0] += D[2] * phi[2];

  eta[1] = D[3] * phi[0];
  eta[1] += D[4] * phi[1];
  eta[1] += D[5] * phi[2];

  eta[2] = D[6] * phi[0];
  eta[2] += D[7] * phi[1];
  eta[2] += D[8] * phi[2];
}

template<typename T>
static inline void mvm_block(std::complex<T> *eta, std::complex<T> *D, std::complex<T> *phi,
                              const int blen){
  for (int j = 0; j < blen; j++) {
    eta[j] = D[0] * phi[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[j] += D[1] * phi[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[j] += D[2] * phi[2*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[blen + j] = D[3] * phi[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] += D[4] * phi[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] += D[5] * phi[2*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] = D[6] * phi[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] += D[7] * phi[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] += D[8] * phi[2*blen + j];
  }
}

template<typename T>
static inline void pbp_su3_T(std::complex<T> *prp_su3_pt, std::complex<T> *l_pt) {
  /// for each line below, we do complex - complex, FLOP is 2
  l_pt[0] -= prp_su3_pt[0];
  l_pt[1] -= prp_su3_pt[1];
  l_pt[2] -= prp_su3_pt[2];
  l_pt[3] -= prp_su3_pt[3];
  l_pt[4] -= prp_su3_pt[4];
  l_pt[5] -= prp_su3_pt[5];
  /// for each line below, we do real * real plus complex * complex plus complex + complex,
  /// FLOP is 1+6+2
  l_pt[6] += (std::complex<T>)GAMMA_T_SPIN2_VAL * prp_su3_pt[3 * GAMMA_T_SPIN2_CO];
  l_pt[7] += (std::complex<T>)GAMMA_T_SPIN2_VAL * prp_su3_pt[3 * GAMMA_T_SPIN2_CO + 1];
  l_pt[8] += (std::complex<T>)GAMMA_T_SPIN2_VAL * prp_su3_pt[3 * GAMMA_T_SPIN2_CO + 2];
  l_pt[9] += (std::complex<T>)GAMMA_T_SPIN3_VAL * prp_su3_pt[3 * GAMMA_T_SPIN3_CO];
  l_pt[10] += (std::complex<T>)GAMMA_T_SPIN3_VAL * prp_su3_pt[3 * GAMMA_T_SPIN3_CO + 1];
  l_pt[11] += (std::complex<T>)GAMMA_T_SPIN3_VAL * prp_su3_pt[3 * GAMMA_T_SPIN3_CO + 2];
}

template<typename T>
static inline void pbp_su3_T_block(std::complex<T> *__restrict__ prp, std::complex<T> *__restrict__ eta, const int blen){
  for (int j = 0; j < blen; j++) {
    eta[j] -= prp[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] -= prp[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] -= prp[2*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] -= prp[3*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] -= prp[4*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] -= prp[5*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] += (std::complex<T>)GAMMA_T_SPIN2_VAL * prp[3*GAMMA_T_SPIN2_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] += (std::complex<T>)GAMMA_T_SPIN2_VAL * prp[(3*GAMMA_T_SPIN2_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] += (std::complex<T>)GAMMA_T_SPIN2_VAL * prp[(3*GAMMA_T_SPIN2_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] += (std::complex<T>)GAMMA_T_SPIN3_VAL * prp[3*GAMMA_T_SPIN3_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] += (std::complex<T>)GAMMA_T_SPIN3_VAL * prp[(3*GAMMA_T_SPIN3_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] += (std::complex<T>)GAMMA_T_SPIN3_VAL * prp[(3*GAMMA_T_SPIN3_CO+2)*blen + j];
  }
}

template<typename T>
static inline void pbp_su3_Z(std::complex<T> *prp_su3_pt, std::complex<T> *l_pt) {
  l_pt[0] -= prp_su3_pt[0];
  l_pt[1] -= prp_su3_pt[1];
  l_pt[2] -= prp_su3_pt[2];
  l_pt[3] -= prp_su3_pt[3];
  l_pt[4] -= prp_su3_pt[4];
  l_pt[5] -= prp_su3_pt[5];
  l_pt[6] += (std::complex<T>)GAMMA_Z_SPIN2_VAL * prp_su3_pt[3 * GAMMA_Z_SPIN2_CO];
  l_pt[7] += (std::complex<T>)GAMMA_Z_SPIN2_VAL * prp_su3_pt[3 * GAMMA_Z_SPIN2_CO + 1];
  l_pt[8] += (std::complex<T>)GAMMA_Z_SPIN2_VAL * prp_su3_pt[3 * GAMMA_Z_SPIN2_CO + 2];
  l_pt[9] += (std::complex<T>)GAMMA_Z_SPIN3_VAL * prp_su3_pt[3 * GAMMA_Z_SPIN3_CO];
  l_pt[10] += (std::complex<T>)GAMMA_Z_SPIN3_VAL * prp_su3_pt[3 * GAMMA_Z_SPIN3_CO + 1];
  l_pt[11] += (std::complex<T>)GAMMA_Z_SPIN3_VAL * prp_su3_pt[3 * GAMMA_Z_SPIN3_CO + 2];
}

template<typename T>
static inline void pbp_su3_Z_block(std::complex<T> *__restrict__ prp, std::complex<T> *__restrict__ eta, const int blen){
  for (int j = 0; j < blen; j++) {
    eta[j] -= prp[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] -= prp[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] -= prp[2*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] -= prp[3*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] -= prp[4*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] -= prp[5*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] += (std::complex<T>)GAMMA_Z_SPIN2_VAL * prp[3*GAMMA_Z_SPIN2_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] += (std::complex<T>)GAMMA_Z_SPIN2_VAL * prp[(3*GAMMA_Z_SPIN2_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] += (std::complex<T>)GAMMA_Z_SPIN2_VAL * prp[(3*GAMMA_Z_SPIN2_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] += (std::complex<T>)GAMMA_Z_SPIN3_VAL * prp[3*GAMMA_Z_SPIN3_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] += (std::complex<T>)GAMMA_Z_SPIN3_VAL * prp[(3*GAMMA_Z_SPIN3_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] += (std::complex<T>)GAMMA_Z_SPIN3_VAL * prp[(3*GAMMA_Z_SPIN3_CO+2)*blen + j];
  }
}

template<typename T>
static inline void pbp_su3_Y(std::complex<T> *prp_su3_pt, std::complex<T> *l_pt) {
  l_pt[0] -= prp_su3_pt[0];
  l_pt[1] -= prp_su3_pt[1];
  l_pt[2] -= prp_su3_pt[2];
  l_pt[3] -= prp_su3_pt[3];
  l_pt[4] -= prp_su3_pt[4];
  l_pt[5] -= prp_su3_pt[5];
  l_pt[6] += (std::complex<T>)GAMMA_Y_SPIN2_VAL * prp_su3_pt[3 * GAMMA_Y_SPIN2_CO];
  l_pt[7] += (std::complex<T>)GAMMA_Y_SPIN2_VAL * prp_su3_pt[3 * GAMMA_Y_SPIN2_CO + 1];
  l_pt[8] += (std::complex<T>)GAMMA_Y_SPIN2_VAL * prp_su3_pt[3 * GAMMA_Y_SPIN2_CO + 2];
  l_pt[9] += (std::complex<T>)GAMMA_Y_SPIN3_VAL * prp_su3_pt[3 * GAMMA_Y_SPIN3_CO];
  l_pt[10] += (std::complex<T>)GAMMA_Y_SPIN3_VAL * prp_su3_pt[3 * GAMMA_Y_SPIN3_CO + 1];
  l_pt[11] += (std::complex<T>)GAMMA_Y_SPIN3_VAL * prp_su3_pt[3 * GAMMA_Y_SPIN3_CO + 2];
}

template<typename T>
static inline void pbp_su3_Y_block(std::complex<T> *__restrict__ prp, std::complex<T> *__restrict__ eta, const int blen){
  for (int j = 0; j < blen; j++) {
    eta[j] -= prp[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] -= prp[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] -= prp[2*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] -= prp[3*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] -= prp[4*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] -= prp[5*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] += (std::complex<T>)GAMMA_Y_SPIN2_VAL * prp[3*GAMMA_Y_SPIN2_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] += (std::complex<T>)GAMMA_Y_SPIN2_VAL * prp[(3*GAMMA_Y_SPIN2_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] += (std::complex<T>)GAMMA_Y_SPIN2_VAL * prp[(3*GAMMA_Y_SPIN2_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] += (std::complex<T>)GAMMA_Y_SPIN3_VAL * prp[3*GAMMA_Y_SPIN3_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] += (std::complex<T>)GAMMA_Y_SPIN3_VAL * prp[(3*GAMMA_Y_SPIN3_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] += (std::complex<T>)GAMMA_Y_SPIN3_VAL * prp[(3*GAMMA_Y_SPIN3_CO+2)*blen + j];
  }
};

template<typename T>
static inline void pbp_su3_X(std::complex<T> *prp_su3_pt, std::complex<T> *l_pt) {
  l_pt[0] -= prp_su3_pt[0];
  l_pt[1] -= prp_su3_pt[1];
  l_pt[2] -= prp_su3_pt[2];
  l_pt[3] -= prp_su3_pt[3];
  l_pt[4] -= prp_su3_pt[4];
  l_pt[5] -= prp_su3_pt[5];
  l_pt[6] += (std::complex<T>)GAMMA_X_SPIN2_VAL * prp_su3_pt[3 * GAMMA_X_SPIN2_CO];
  l_pt[7] += (std::complex<T>)GAMMA_X_SPIN2_VAL * prp_su3_pt[3 * GAMMA_X_SPIN2_CO + 1];
  l_pt[8] += (std::complex<T>)GAMMA_X_SPIN2_VAL * prp_su3_pt[3 * GAMMA_X_SPIN2_CO + 2];
  l_pt[9] += (std::complex<T>)GAMMA_X_SPIN3_VAL * prp_su3_pt[3 * GAMMA_X_SPIN3_CO];
  l_pt[10] += (std::complex<T>)GAMMA_X_SPIN3_VAL * prp_su3_pt[3 * GAMMA_X_SPIN3_CO + 1];
  l_pt[11] += (std::complex<T>)GAMMA_X_SPIN3_VAL * prp_su3_pt[3 * GAMMA_X_SPIN3_CO + 2];
}

template<typename T>
static inline void pbp_su3_X_block(std::complex<T> *__restrict__ prp, std::complex<T> *__restrict__ eta, const int blen){
  for (int j = 0; j < blen; j++) {
    eta[j] -= prp[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] -= prp[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] -= prp[2*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] -= prp[3*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] -= prp[4*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] -= prp[5*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] += (std::complex<T>)GAMMA_X_SPIN2_VAL * prp[3*GAMMA_X_SPIN2_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] += (std::complex<T>)GAMMA_X_SPIN2_VAL * prp[(3*GAMMA_X_SPIN2_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] += (std::complex<T>)GAMMA_X_SPIN2_VAL * prp[(3*GAMMA_X_SPIN2_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] += (std::complex<T>)GAMMA_X_SPIN3_VAL * prp[3*GAMMA_X_SPIN3_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] += (std::complex<T>)GAMMA_X_SPIN3_VAL * prp[(3*GAMMA_X_SPIN3_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] += (std::complex<T>)GAMMA_X_SPIN3_VAL * prp[(3*GAMMA_X_SPIN3_CO+2)*blen + j];
  }
}

template<typename T>
static inline void pbn_su3_T(std::complex<T> *prn_su3_pt, std::complex<T> *l_pt) {
  l_pt[0] -= prn_su3_pt[0];
  l_pt[1] -= prn_su3_pt[1];
  l_pt[2] -= prn_su3_pt[2];
  l_pt[3] -= prn_su3_pt[3];
  l_pt[4] -= prn_su3_pt[4];
  l_pt[5] -= prn_su3_pt[5];
  l_pt[6] -= (std::complex<T>)GAMMA_T_SPIN2_VAL * prn_su3_pt[3 * GAMMA_T_SPIN2_CO];
  l_pt[7] -= (std::complex<T>)GAMMA_T_SPIN2_VAL * prn_su3_pt[3 * GAMMA_T_SPIN2_CO + 1];
  l_pt[8] -= (std::complex<T>)GAMMA_T_SPIN2_VAL * prn_su3_pt[3 * GAMMA_T_SPIN2_CO + 2];
  l_pt[9] -= (std::complex<T>)GAMMA_T_SPIN3_VAL * prn_su3_pt[3 * GAMMA_T_SPIN3_CO];
  l_pt[10] -= (std::complex<T>)GAMMA_T_SPIN3_VAL * prn_su3_pt[3 * GAMMA_T_SPIN3_CO + 1];
  l_pt[11] -= (std::complex<T>)GAMMA_T_SPIN3_VAL * prn_su3_pt[3 * GAMMA_T_SPIN3_CO + 2];
}

template<typename T>
static inline void pbn_su3_T_block(std::complex<T> *__restrict__ prn, std::complex<T> *__restrict__ eta, const int blen){
  for (int j = 0; j < blen; j++) {
    eta[j] -= prn[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] -= prn[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] -= prn[2*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] -= prn[3*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] -= prn[4*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] -= prn[5*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] -= (std::complex<T>)GAMMA_T_SPIN2_VAL *
                         prn[3*GAMMA_T_SPIN2_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] -= (std::complex<T>)GAMMA_T_SPIN2_VAL *
                         prn[(3*GAMMA_T_SPIN2_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] -= (std::complex<T>)GAMMA_T_SPIN2_VAL *
                         prn[(3*GAMMA_T_SPIN2_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] -= (std::complex<T>)GAMMA_T_SPIN3_VAL *
                         prn[3*GAMMA_T_SPIN3_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] -= (std::complex<T>)GAMMA_T_SPIN3_VAL *
                          prn[(3*GAMMA_T_SPIN3_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] -= (std::complex<T>)GAMMA_T_SPIN3_VAL *
                          prn[(3*GAMMA_T_SPIN3_CO+2)*blen + j];
  }
}

template<typename T>
static inline void pbn_su3_Z(std::complex<T> *prn_su3_pt, std::complex<T> *l_pt) {
  l_pt[0] -= prn_su3_pt[0];
  l_pt[1] -= prn_su3_pt[1];
  l_pt[2] -= prn_su3_pt[2];
  l_pt[3] -= prn_su3_pt[3];
  l_pt[4] -= prn_su3_pt[4];
  l_pt[5] -= prn_su3_pt[5];
  l_pt[6] -= (std::complex<T>)GAMMA_Z_SPIN2_VAL * prn_su3_pt[3 * GAMMA_Z_SPIN2_CO];
  l_pt[7] -= (std::complex<T>)GAMMA_Z_SPIN2_VAL * prn_su3_pt[3 * GAMMA_Z_SPIN2_CO + 1];
  l_pt[8] -= (std::complex<T>)GAMMA_Z_SPIN2_VAL * prn_su3_pt[3 * GAMMA_Z_SPIN2_CO + 2];
  l_pt[9] -= (std::complex<T>)GAMMA_Z_SPIN3_VAL * prn_su3_pt[3 * GAMMA_Z_SPIN3_CO];
  l_pt[10] -= (std::complex<T>)GAMMA_Z_SPIN3_VAL * prn_su3_pt[3 * GAMMA_Z_SPIN3_CO + 1];
  l_pt[11] -= (std::complex<T>)GAMMA_Z_SPIN3_VAL * prn_su3_pt[3 * GAMMA_Z_SPIN3_CO + 2];
}

template<typename T>
static inline void pbn_su3_Z_block(std::complex<T> *__restrict__ prn, std::complex<T> *__restrict__ eta, const int blen){
  for (int j = 0; j < blen; j++) {
    eta[j] -= prn[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] -= prn[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] -= prn[2*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] -= prn[3*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] -= prn[4*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] -= prn[5*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] -= (std::complex<T>)GAMMA_Z_SPIN2_VAL *
                         prn[3*GAMMA_Z_SPIN2_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] -= (std::complex<T>)GAMMA_Z_SPIN2_VAL *
                         prn[(3*GAMMA_Z_SPIN2_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] -= (std::complex<T>)GAMMA_Z_SPIN2_VAL *
                         prn[(3*GAMMA_Z_SPIN2_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] -= (std::complex<T>)GAMMA_Z_SPIN3_VAL *
                         prn[3*GAMMA_Z_SPIN3_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] -= (std::complex<T>)GAMMA_Z_SPIN3_VAL *
                          prn[(3*GAMMA_Z_SPIN3_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] -= (std::complex<T>)GAMMA_Z_SPIN3_VAL *
                          prn[(3*GAMMA_Z_SPIN3_CO+2)*blen + j];
  }
}

template<typename T>
static inline void pbn_su3_Y(std::complex<T> *prn_su3_pt, std::complex<T> *l_pt) {
  l_pt[0] -= prn_su3_pt[0];
  l_pt[1] -= prn_su3_pt[1];
  l_pt[2] -= prn_su3_pt[2];
  l_pt[3] -= prn_su3_pt[3];
  l_pt[4] -= prn_su3_pt[4];
  l_pt[5] -= prn_su3_pt[5];
  l_pt[6] -= (std::complex<T>)GAMMA_Y_SPIN2_VAL * prn_su3_pt[3 * GAMMA_Y_SPIN2_CO];
  l_pt[7] -= (std::complex<T>)GAMMA_Y_SPIN2_VAL * prn_su3_pt[3 * GAMMA_Y_SPIN2_CO + 1];
  l_pt[8] -= (std::complex<T>)GAMMA_Y_SPIN2_VAL * prn_su3_pt[3 * GAMMA_Y_SPIN2_CO + 2];
  l_pt[9] -= (std::complex<T>)GAMMA_Y_SPIN3_VAL * prn_su3_pt[3 * GAMMA_Y_SPIN3_CO];
  l_pt[10] -= (std::complex<T>)GAMMA_Y_SPIN3_VAL * prn_su3_pt[3 * GAMMA_Y_SPIN3_CO + 1];
  l_pt[11] -= (std::complex<T>)GAMMA_Y_SPIN3_VAL * prn_su3_pt[3 * GAMMA_Y_SPIN3_CO + 2];
}

template<typename T>
static inline void pbn_su3_Y_block(std::complex<T> *__restrict__ prn, std::complex<T> *__restrict__ eta, const int blen){
  for (int j = 0; j < blen; j++) {
    eta[j] -= prn[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] -= prn[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] -= prn[2*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] -= prn[3*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] -= prn[4*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] -= prn[5*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] -= (std::complex<T>)GAMMA_Y_SPIN2_VAL *
                         prn[3*GAMMA_Y_SPIN2_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] -= (std::complex<T>)GAMMA_Y_SPIN2_VAL *
                         prn[(3*GAMMA_Y_SPIN2_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] -= (std::complex<T>)GAMMA_Y_SPIN2_VAL *
                         prn[(3*GAMMA_Y_SPIN2_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] -= (std::complex<T>)GAMMA_Y_SPIN3_VAL *
                         prn[3*GAMMA_Y_SPIN3_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] -= (std::complex<T>)GAMMA_Y_SPIN3_VAL *
                          prn[(3*GAMMA_Y_SPIN3_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] -= (std::complex<T>)GAMMA_Y_SPIN3_VAL *
                          prn[(3*GAMMA_Y_SPIN3_CO+2)*blen + j];
  }
}

template<typename T>
static inline void pbn_su3_X(std::complex<T> *prn_su3_pt, std::complex<T> *l_pt) {
  l_pt[0] -= prn_su3_pt[0];
  l_pt[1] -= prn_su3_pt[1];
  l_pt[2] -= prn_su3_pt[2];
  l_pt[3] -= prn_su3_pt[3];
  l_pt[4] -= prn_su3_pt[4];
  l_pt[5] -= prn_su3_pt[5];
  l_pt[6] -= (std::complex<T>)GAMMA_X_SPIN2_VAL * prn_su3_pt[3 * GAMMA_X_SPIN2_CO];
  l_pt[7] -= (std::complex<T>)GAMMA_X_SPIN2_VAL * prn_su3_pt[3 * GAMMA_X_SPIN2_CO + 1];
  l_pt[8] -= (std::complex<T>)GAMMA_X_SPIN2_VAL * prn_su3_pt[3 * GAMMA_X_SPIN2_CO + 2];
  l_pt[9] -= (std::complex<T>)GAMMA_X_SPIN3_VAL * prn_su3_pt[3 * GAMMA_X_SPIN3_CO];
  l_pt[10] -= (std::complex<T>)GAMMA_X_SPIN3_VAL * prn_su3_pt[3 * GAMMA_X_SPIN3_CO + 1];
  l_pt[11] -= (std::complex<T>)GAMMA_X_SPIN3_VAL * prn_su3_pt[3 * GAMMA_X_SPIN3_CO + 2];
}

template<typename T>
static inline void pbn_su3_X_block(std::complex<T> *__restrict__ prn, std::complex<T> *__restrict__ eta, const int blen){
  for (int j = 0; j < blen; j++) {
    eta[j] -= prn[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[blen + j] -= prn[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] -= prn[2*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] -= prn[3*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] -= prn[4*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] -= prn[5*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] -= (std::complex<T>)GAMMA_X_SPIN2_VAL *
                               prn[3*GAMMA_X_SPIN2_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] -= (std::complex<T>)GAMMA_X_SPIN2_VAL *
                               prn[(3*GAMMA_X_SPIN2_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] -= (std::complex<T>)GAMMA_X_SPIN2_VAL *
                               prn[(3*GAMMA_X_SPIN2_CO+2)*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] -= (std::complex<T>)GAMMA_X_SPIN3_VAL *
                               prn[3*GAMMA_X_SPIN3_CO*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] -= (std::complex<T>)GAMMA_X_SPIN3_VAL *
                               prn[(3*GAMMA_X_SPIN3_CO+1)*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] -= (std::complex<T>)GAMMA_X_SPIN3_VAL *
                               prn[(3*GAMMA_X_SPIN3_CO+2)*blen + j];
  }
}

template<typename T>
static inline void site_clover(std::complex<T> *eta, std::complex<T> *phi,
                                      std::complex<T> *clover) {
  // diagonal
  eta[0] = clover[0] * phi[0];
  eta[1] = clover[1] * phi[1];
  eta[2] = clover[2] * phi[2];
  eta[3] = clover[3] * phi[3];
  eta[4] = clover[4] * phi[4];
  eta[5] = clover[5] * phi[5];
  eta[6] = clover[6] * phi[6];
  eta[7] = clover[7] * phi[7];
  eta[8] = clover[8] * phi[8];
  eta[9] = clover[9] * phi[9];
  eta[10] = clover[10] * phi[10];
  eta[11] = clover[11] * phi[11];
  // spin 0 and 1, row major
  eta[0] += clover[12] * phi[1];
  eta[0] += clover[13] * phi[2];
  eta[0] += clover[14] * phi[3];
  eta[0] += clover[15] * phi[4];
  eta[0] += clover[16] * phi[5];
  eta[1] += clover[17] * phi[2];
  eta[1] += clover[18] * phi[3];
  eta[1] += clover[19] * phi[4];
  eta[1] += clover[20] * phi[5];
  eta[2] += clover[21] * phi[3];
  eta[2] += clover[22] * phi[4];
  eta[2] += clover[23] * phi[5];
  eta[3] += clover[24] * phi[4];
  eta[3] += clover[25] * phi[5];
  eta[4] += clover[26] * phi[5];
  eta[1] += conj(clover[12]) * phi[0];
  eta[2] += conj(clover[13]) * phi[0];
  eta[3] += conj(clover[14]) * phi[0];
  eta[4] += conj(clover[15]) * phi[0];
  eta[5] += conj(clover[16]) * phi[0];
  eta[2] += conj(clover[17]) * phi[1];
  eta[3] += conj(clover[18]) * phi[1];
  eta[4] += conj(clover[19]) * phi[1];
  eta[5] += conj(clover[20]) * phi[1];
  eta[3] += conj(clover[21]) * phi[2];
  eta[4] += conj(clover[22]) * phi[2];
  eta[5] += conj(clover[23]) * phi[2];
  eta[4] += conj(clover[24]) * phi[3];
  eta[5] += conj(clover[25]) * phi[3];
  eta[5] += conj(clover[26]) * phi[4];
  // spin 2 and 3, row major
  eta[6] += clover[27] * phi[7];
  eta[6] += clover[28] * phi[8];
  eta[6] += clover[29] * phi[9];
  eta[6] += clover[30] * phi[10];
  eta[6] += clover[31] * phi[11];
  eta[7] += clover[32] * phi[8];
  eta[7] += clover[33] * phi[9];
  eta[7] += clover[34] * phi[10];
  eta[7] += clover[35] * phi[11];
  eta[8] += clover[36] * phi[9];
  eta[8] += clover[37] * phi[10];
  eta[8] += clover[38] * phi[11];
  eta[9] += clover[39] * phi[10];
  eta[9] += clover[40] * phi[11];
  eta[10] += clover[41] * phi[11];
  eta[7] += conj(clover[27]) * phi[6];
  eta[8] += conj(clover[28]) * phi[6];
  eta[9] += conj(clover[29]) * phi[6];
  eta[10] += conj(clover[30]) * phi[6];
  eta[11] += conj(clover[31]) * phi[6];
  eta[8] += conj(clover[32]) * phi[7];
  eta[9] += conj(clover[33]) * phi[7];
  eta[10] += conj(clover[34]) * phi[7];
  eta[11] += conj(clover[35]) * phi[7];
  eta[9] += conj(clover[36]) * phi[8];
  eta[10] += conj(clover[37]) * phi[8];
  eta[11] += conj(clover[38]) * phi[8];
  eta[10] += conj(clover[39]) * phi[9];
  eta[11] += conj(clover[40]) * phi[9];
  eta[11] += conj(clover[41]) * phi[10];
}


template<typename T>
static inline void site_clover_block(std::complex<T> *__restrict__ eta, std::complex<T> *__restrict__ phi,
                               std::complex<T> *__restrict__ clover, const int blen) {
  /// diagonal
  for (int j = 0; j < blen; ++j) {
    eta[j] = clover[0] * phi[j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[blen + j] = clover[1] * phi[blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[2*blen + j] = clover[2] * phi[2*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[3*blen + j] = clover[3] * phi[3*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[4*blen + j] = clover[4] * phi[4*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[5*blen + j] = clover[5] * phi[5*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[6*blen + j] = clover[6] * phi[6*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[7*blen + j] = clover[7] * phi[7*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[8*blen + j] = clover[8] * phi[8*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[9*blen + j] = clover[9] * phi[9*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[10*blen + j] = clover[10] * phi[10*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[11*blen + j] = clover[11] * phi[11*blen + j];
  }

  /// spin 0 and 1, row major
  for (int j = 0; j < blen; ++j) {
    eta[j] += clover[12] * phi[blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[j] += clover[13] * phi[2*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[j] += clover[14] * phi[3*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[j] += clover[15] * phi[4*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[j] += clover[16] * phi[5*blen + j];
  }

  for (int j = 0; j < blen; ++j) {
    eta[blen + j] += clover[17] * phi[2*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[blen + j] += clover[18] * phi[3*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[blen + j] += clover[19] * phi[4*blen + j];
  }
  for (int j = 0; j < blen; ++j) {
    eta[blen + j] += clover[20] * phi[5*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] += clover[21] * phi[3*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] += clover[22] * phi[4*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] += clover[23] * phi[5*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] += clover[24] * phi[4*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] += clover[25] * phi[5*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] += clover[26] * phi[5*blen + j];
  }
  /// lower triangular matrix
  for (int j = 0; j < blen; j++) {
    eta[blen + j] += conj(clover[12]) * phi[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] += conj(clover[13]) * phi[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] += conj(clover[14]) * phi[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] += conj(clover[15]) * phi[j];
  }
  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] += conj(clover[16]) * phi[j];
  }

  for (int j = 0; j < blen; j++) {
    eta[2*blen + j] += conj(clover[17]) * phi[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] += conj(clover[18]) * phi[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] += conj(clover[19]) * phi[blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] += conj(clover[20]) * phi[blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[3*blen + j] += conj(clover[21]) * phi[2*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] += conj(clover[22]) * phi[2*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] += conj(clover[23]) * phi[2*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[4*blen + j] += conj(clover[24]) * phi[3*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] += conj(clover[25]) * phi[3*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[5*blen + j] += conj(clover[26]) * phi[4*blen + j];
  }

  /// spin 2 and 3, row major
  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] += clover[27] * phi[7*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] += clover[28] * phi[8*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] += clover[29] * phi[9*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] += clover[30] * phi[10*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[6*blen + j] += clover[31] * phi[11*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] += clover[32] * phi[8*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] += clover[33] * phi[9*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] += clover[34] * phi[10*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] += clover[35] * phi[11*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] += clover[36] * phi[9*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] += clover[37] * phi[10*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] += clover[38] * phi[11*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] += clover[39] * phi[10*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] += clover[40] * phi[11*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] += clover[41] * phi[11*blen + j];
  }
  /// lower triangular matrix
  for (int j = 0; j < blen; j++) {
    eta[7*blen + j] += conj(clover[27]) * phi[6*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] += conj(clover[28]) * phi[6*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] += conj(clover[29]) * phi[6*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] += conj(clover[30]) * phi[6*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] += conj(clover[31]) * phi[6*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[8*blen + j] += conj(clover[32]) * phi[7*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] += conj(clover[33]) * phi[7*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] += conj(clover[34]) * phi[7*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] += conj(clover[35]) * phi[7*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[9*blen + j] += conj(clover[36]) * phi[8*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] += conj(clover[37]) * phi[8*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] += conj(clover[38]) * phi[8*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[10*blen + j] += conj(clover[39]) * phi[9*blen + j];
  }
  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] += conj(clover[40]) * phi[9*blen + j];
  }

  for (int j = 0; j < blen; j++) {
    eta[11*blen + j] += conj(clover[41]) * phi[10*blen + j];
  }
}

#endif // DDALPHAAMG_DIRAC_TEMPLATED_GENERIC_H
