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
#ifndef DDALPHAAMG_DIRAC_BLOCK_VPSITE_H
#define DDALPHAAMG_DIRAC_BLOCK_VPSITE_H
#include "omp.h"

#include "dirac_templated_generic.h"
#include "block.h"
#include "mpi.h"

struct mpi_type_handler {
  static MPI_Datatype get_datatype(std::complex<double>* x){
    return MPI_COMPLEX_double;
  }

  static MPI_Datatype get_datatype(std::complex<float>* x){
    return MPI_COMPLEX_float;
  }

  template<typename T>
  MPI_Datatype get_datatype(T x){
    throw "MPI type not defined";
  }
};


template<typename T>
void ghost_sendrecv_full(std::complex<T> *phi, const int mu, const int dir, comm_struct<T> *c,
                         const int amount, Level *l, const int blen) {

  Global *global = l->global;

  // does not allow sending in both directions at the same time
  if (l->global_splitting[mu] > 1) {

    int *table = nullptr, mu_dir = 2 * mu - MIN(dir, 0), offset = c->offset,
              length[2] = {0, 0}, comm_start = 0;
    complex<T> *buffer, *phi_pt;

    if (amount == _FULL_SYSTEM) {
      length[0] = (c->num_boundary_sites[2 * mu]) * offset * blen;
      // num_boundary_sites[2 * mu]:num_boundary in -mu dir
      length[1] = (c->num_boundary_sites[2 * mu + 1]) * offset * blen;
      // num_boundary_sites[2 * mu + 1]: num_boundary in +mu dir
      comm_start = c->comm_start[mu_dir] * blen;
      // comm_start: start index of +mu boundary
    } else if (amount == _EVEN_SITES) {
      length[0] = c->num_even_boundary_sites[2 * mu] * offset * blen;
      length[1] = c->num_even_boundary_sites[2 * mu + 1] * offset * blen;
      comm_start = c->comm_start[mu_dir] * blen;
    } else if (amount == _ODD_SITES) {
      length[0] = c->num_odd_boundary_sites[2 * mu] * offset * blen;
      length[1] = c->num_odd_boundary_sites[2 * mu + 1] * offset * blen;
      comm_start = (c->comm_start[mu_dir] + c->num_even_boundary_sites[mu_dir] * offset) * blen;
    }

    ASSERT(c->in_use[mu_dir] == 0);
    c->in_use[mu_dir] = 1;

    if (MAX(length[0], length[1]) > c->max_length[mu]) {
      printf("CAUTION: my_rank: %d, not enough comm buffer\n", global->my_rank);
      fflush(0);
      ghost_free(c, l);
      ghost_alloc(MAX(length[0], length[1]), c, l);
    }

    buffer = (complex<T>*)c->buffer[mu_dir];

    // dir = senddir
    if (dir == 1) {
      // data to be communicated is stored serially in the vector phi
      // recv target is a buffer
      // afterwards (in ghost_wait) the data has to be distributed onto the correct sites
      // touching the respective boundary in -mu direction

      phi_pt = phi + comm_start;
      if (length[1] > 0) {
        MPI_Irecv(buffer, length[1], mpi_type_handler().get_datatype(phi), l->neighbor_rank[2 * mu + 1], 2 * mu,
                  global->comm_cart, &(c->rreqs[2 * mu]));
      }
      if (length[0] > 0) {
        MPI_Isend(phi_pt, length[0], mpi_type_handler().get_datatype(phi), l->neighbor_rank[2 * mu], 2 * mu,
                  global->comm_cart, &(c->sreqs[2 * mu]));
      }

    } else if (dir == -1) {
      // data to be communicated is stored on the sites touching the boundary in -mu direction
      // this data is gathered in a buffer in the correct ordering
      // which is required on the boundary of the vector phi
      int num_boundary_sites = length[1] / offset / blen;
      //table stores the -mu boundary
      if (amount == _ODD_SITES){
        table = c->boundary_table[2 * mu + 1] + c->num_even_boundary_sites[mu_dir];
      } else {
        table = c->boundary_table[2 * mu + 1];
      }
#pragma omp parallel for shared(num_boundary_sites, offset, buffer, phi, table)
      for (int j = 0; j < num_boundary_sites; j++) {
        for (int k = 0; k < blen; k++) {
          for (int i = 0; i < offset; i++) {
            buffer[j*blen*offset + k*offset + i] = phi[table[j]*blen*offset + k*offset + i];
          }
        }
      }

      buffer = (std::complex<T>*)c->buffer[mu_dir];
      phi_pt = phi + comm_start;

      if (length[0] > 0) {
        MPI_Irecv(phi_pt, length[0], mpi_type_handler().get_datatype(phi), l->neighbor_rank[2 * mu], 2 * mu + 1,
                  global->comm_cart, &(c->rreqs[2 * mu + 1]));

      }
      if (length[1] > 0) {
        MPI_Isend(buffer, length[1], mpi_type_handler().get_datatype(phi), l->neighbor_rank[2 * mu + 1], 2 * mu + 1,
                  global->comm_cart, &(c->sreqs[2 * mu + 1]));
      }

    } else
      ASSERT(dir == 1 || dir == -1);
  }
}

template<typename T>
void ghost_wait_full(std::complex<T> *phi, const int mu, const int dir, comm_struct<T> *c,
                     const int amount, Level *l, const int blen) {

  if (l->global_splitting[mu] > 1) {
    int mu_dir = 2 * mu - MIN(dir, 0);
    int i, j, *table, offset = c->offset, length[2] = {0, 0};
    std::complex<T> *buffer;

    if (amount == _FULL_SYSTEM) {
      length[0] = (c->num_boundary_sites[2 * mu]) * offset * blen;
      // num_boundary_sites[2 * mu]:num_boundary in -mu dir
      length[1] = (c->num_boundary_sites[2 * mu + 1]) * offset * blen;
      // num_boundary_sites[2 * mu + 1]: num_boundary in +mu dir
    } else if (amount == _EVEN_SITES) {
      length[0] = c->num_even_boundary_sites[2 * mu] * offset * blen;
      length[1] = c->num_even_boundary_sites[2 * mu + 1] * offset * blen;
    } else if (amount == _ODD_SITES) {
      length[0] = c->num_odd_boundary_sites[2 * mu] * offset * blen;
      length[1] = c->num_odd_boundary_sites[2 * mu + 1] * offset * blen;
    }

    ASSERT(c->in_use[mu_dir] == 1);

    if (dir == 1) {

      int num_boundary_sites = length[0] / offset / blen;

      buffer = (std::complex<T>*)c->buffer[mu_dir];
      if (amount == _ODD_SITES){
        table = c->boundary_table[2 * mu + 1] + c->num_even_boundary_sites[mu_dir];
      } else{
        table = c->boundary_table[2 * mu + 1];
      }

      if (length[0] > 0) {
        MPI_Wait(&(c->sreqs[2 * mu]), MPI_STATUS_IGNORE);
      }
      if (length[1] > 0) {
        MPI_Wait(&(c->rreqs[2 * mu]), MPI_STATUS_IGNORE);
      }

      if (l->depth == 0) {
#pragma omp parallel for shared(num_boundary_sites, offset, buffer, phi, table)
        for (j = 0; j < num_boundary_sites; j++) {
          for (int k = 0; k < blen; k++) {
            for (i = 0; i < offset; i++) {
              phi[table[j]*blen*offset + k*offset + i] = buffer[j*blen*offset + k*offset + i];
            }
          }
        }
      } else {
#pragma omp parallel for shared(num_boundary_sites, offset, buffer, phi, table)
        for (j = 0; j < num_boundary_sites; j++) {
          for (int k = 0; k < blen; k++) {
            for (i = 0; i < offset; i++) {
              phi[table[j]*blen*offset + k*offset + i] += buffer[j*blen*offset + k*offset + i];
            }
          }
        }
      }
    } else if (dir == -1) {
      if (length[1] > 0) {
        MPI_Wait(&(c->sreqs[2 * mu + 1]), MPI_STATUS_IGNORE);
      }
      if (length[0] > 0) {
        MPI_Wait(&(c->rreqs[2 * mu + 1]), MPI_STATUS_IGNORE);
      }

    } else
      ASSERT(dir == 1 || dir == -1);
    c->in_use[mu_dir] = 0;
  }
}

template<typename T>
static inline void clover_block(std::complex<T> *eta, std::complex<T> *phi, operator_struct<T> *op,
                   int end, Level *l, const int blen) {

  if (l->global->csw == 0.0) {
    std::complex<T> *clover = op->clover;
#pragma omp parallel for shared(end, eta, phi, clover)
    for (int i = 0; i < end/12; i++) {
      for (int j = 0; j < blen; j++){
        for(int k = 0; k < 12; k++){
          eta[i*12*blen+j*12+k] = phi[i*12*blen+j*12+k] * clover[i*12+k];
        }
      }
    }


  } else {
    std::complex<T> *clover = op->clover;
#pragma omp parallel for shared(end, eta, phi, clover)
    for (int i = 0; i < end/12; i++) {
      for (int j = 0; j < blen; j++) {
        site_clover(&eta[i*blen*12 + j*12], &phi[i*blen*12 + j*12], &clover[i*42]);
      }
    }

  }
}

template<typename T>
static inline void pi_minus(std::complex<T>* phi, std::complex<T>* prnT, std::complex<T>* prnZ,
                            std::complex<T>* prnY, std::complex<T>* prnX, int start, int end){
#pragma omp parallel for shared(start, end, phi, prnT, prnZ, prnY, prnX)
  for (int i = start/2; i < end/2; i += 6) {
        prp_T(&prnT[i], &phi[i*2]);
        prp_Z(&prnZ[i], &phi[i*2]);
        prp_Y(&prnY[i], &phi[i*2]);
        prp_X(&prnX[i], &phi[i*2]);
  }
}


template<typename T>
static inline void pi_plus_kron_UH_psi(std::complex<T>* phi, std::complex<T>* D, std::complex<T>* prpT,
                                       std::complex<T>* prpZ, std::complex<T>* prpY, std::complex<T>* prpX,
                                       const int* neighbor, int start, int end, const int blen){
  std::complex<T> pbuf[6];
  int nt_id;

#pragma omp parallel for shared(start, end, phi, D, neighbor, prpT, prpZ, prpY, prpX) private(pbuf, nt_id)
  for (int i = start; i < end; i+=12) {
    nt_id = i/3;  /// nt_id is in range [site_nr * 4, site_nr * 4 + 3] in each iter, site_nr = i/12 here
    for (int j = 0; j < blen; j++) {
      /// T dir
      prn_T(pbuf, &phi[i*blen+j*12]); /// (I4 + Gamma_mu) kronecker I3 * psi(x)
      mvmh(&prpT[blen*6*neighbor[nt_id]+j*6], &D[9*nt_id], pbuf); /// apply U_mu^H(x), first half
      mvmh(&prpT[blen*6*neighbor[nt_id]+j*6 + 3], &D[9*nt_id], &pbuf[3]); /// second half

      /// Z dir
      prn_Z(pbuf, &phi[i*blen+j*12]);
      mvmh(&prpZ[blen*6*neighbor[nt_id+1]+j*6], &D[9*(nt_id+1)], pbuf);
      mvmh(&prpZ[blen*6*neighbor[nt_id+1]+j*6 + 3], &D[9*(nt_id+1)], &pbuf[3]);

      /// Y dir
      prn_Y(pbuf, &phi[i*blen+j*12]);
      mvmh(&prpY[blen*6*neighbor[nt_id+2]+j*6], &D[9*(nt_id+2)], pbuf);
      mvmh(&prpY[blen*6*neighbor[nt_id+2]+j*6 + 3], &D[9*(nt_id+2)], &pbuf[3]);

      /// X dir
      prn_X(pbuf, &phi[i*blen+j*12]);
      mvmh(&prpX[blen*6*neighbor[nt_id+3]+j*6], &D[9*(nt_id+3)], pbuf);
      mvmh(&prpX[blen*6*neighbor[nt_id+3]+j*6 + 3], &D[9*(nt_id+3)], &pbuf[3]);
    }
  }
}

template<typename T>
static inline void apply_U_psi_plus_mu(std::complex<T>* eta, std::complex<T>* D, std::complex<T>* prnT,
                                       std::complex<T>* prnZ, std::complex<T>* prnY, std::complex<T>* prnX,
                                       const int* neighbor, int start, int end, const int blen){
  std::complex<T> pbuf[6];
  int nt_id;

#pragma omp parallel for shared(start, end, eta, D, neighbor, prnT, prnZ, prnY, prnX) private(pbuf, nt_id)
  for (int i = start; i < end; i+=12) {
    nt_id = i/3;
    for (int j = 0; j < blen; j++) {
      /// T dir
      mvm(pbuf, &D[9*nt_id], &prnT[blen*6*neighbor[nt_id]+j*6]); /// apply U_mu(x) to prn, first half of x
      mvm(&pbuf[3], &D[9*nt_id], &prnT[blen*6*neighbor[nt_id]+j*6 + 3]); /// second half of x
      /// each line of second half of eta can be represented by lines of first half of eta,
      /// here we apply the effect of pi_mu_minus kronecker U_mu(x) * psi(x+mu_hat)
      pbp_su3_T(pbuf, &eta[i*blen+j*12]);

      /// Z dir
      mvm(pbuf, &D[9*(nt_id+1)], &prnZ[blen*6*neighbor[nt_id+1]+j*6]);
      mvm(&pbuf[3], &D[9*(nt_id+1)], &prnZ[blen*6*neighbor[nt_id+1]+j*6] + 3);
      pbp_su3_Z(pbuf, &eta[i*blen+j*12]);

      /// Y dir
      mvm(pbuf, &D[9*(nt_id+2)], &prnY[blen*6*neighbor[nt_id+2]+j*6]);
      mvm(&pbuf[3], &D[9*(nt_id+2)], &prnY[blen*6*neighbor[nt_id+2]+j*6 + 3]);
      pbp_su3_Y(pbuf, &eta[i*blen+j*12]);

      /// X dir
      mvm(pbuf, &D[9*(nt_id+3)], &prnX[blen*6*neighbor[nt_id+3]+j*6]);
      mvm(&pbuf[3], &D[9*(nt_id+3)], &prnX[blen*6*neighbor[nt_id+3]+j*6 + 3]);
      pbp_su3_X(pbuf, &eta[i*blen+j*12]);
    }
  }
}

template<typename T>
static inline void lift_pi_plus_dir(std::complex<T>* eta, std::complex<T>* prpT, std::complex<T>* prpZ,
                                    std::complex<T>* prpY, std::complex<T>* prpX, int start, int end){
#pragma omp parallel for shared(start, end, prpT, prpZ, prpY, prpX, eta)
  for (int i = start / 2; i < end / 2; i += 6) {
    /// each line of second half of eta can be represented by lines of first half of eta,
    /// here we apply the effect of pi_mu_plus kronecker U_mu(x-mu_hat) * psi(x-mu_hat)
    pbn_su3_T(&prpT[i], &eta[2*i]);
    pbn_su3_Z(&prpZ[i], &eta[2*i]);
    pbn_su3_Y(&prpY[i], &eta[2*i]);
    pbn_su3_X(&prpX[i], &eta[2*i]);
  }
}

template<typename T>
void d_plus_clover_block(BlockVecPerSite<T> *bout, BlockVecPerSite<T> *bv, operator_struct<T> *op,
                         Level *l, const int blen){
  int n = l->num_inner_lattice_sites, *neighbor = op->neighbor_table, nv = l->num_lattice_site_var;
  int end =  nv * n;
  std::complex<T> *phi = bv->v;
  std::complex<T> *eta = bout->v;

  auto *prnT = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];
  auto *prnZ = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];
  auto *prnY = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];
  auto *prnX = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];

  auto *prpT = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];
  auto *prpZ = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];
  auto *prpY = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];
  auto *prpX = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];

  clover_block(eta, phi, op, end, l, blen);

  pi_minus(phi, prnT, prnZ, prnY, prnX, 0, end*blen);

  ///send prn to mu-1 direction, receive corresponding prn from mu+1 neighbor
  ghost_sendrecv_full(prnT, _T, -1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_sendrecv_full(prnZ, _Z, -1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_sendrecv_full(prnY, _Y, -1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_sendrecv_full(prnX, _X, -1, &(op->c), _FULL_SYSTEM, l, blen);

  pi_plus_kron_UH_psi(phi, op->D, prpT, prpZ, prpY, prpX, neighbor, 0, end, blen);

  /// send prp to mu+1 neighbour, receive corresponding prp from mu-1 neighbor
  ghost_sendrecv_full(prpT, _T, +1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_sendrecv_full(prpZ, _Z, +1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_sendrecv_full(prpY, _Y, +1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_sendrecv_full(prpX, _X, +1, &(op->c), _FULL_SYSTEM, l, blen);

  /// wait for communication of prn
  ghost_wait_full(prnT, _T, -1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_wait_full(prnZ, _Z, -1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_wait_full(prnY, _Y, -1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_wait_full(prnX, _X, -1, &(op->c), _FULL_SYSTEM, l, blen);

  apply_U_psi_plus_mu(eta, op->D, prnT, prnZ, prnY, prnX, neighbor, 0, end, blen);

  /// wait for communication of prp
  ghost_wait_full(prpT, _T, +1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_wait_full(prpZ, _Z, +1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_wait_full(prpY, _Y, +1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_wait_full(prpX, _X, +1, &(op->c), _FULL_SYSTEM, l, blen);

  lift_pi_plus_dir(eta, prpT, prpZ, prpY, prpX, 0, end*blen);

  delete[] prnT;
  delete[] prnZ;
  delete[] prnY;
  delete[] prnX;

  delete[] prpT;
  delete[] prpZ;
  delete[] prpY;
  delete[] prpX;
}


// below are auxiliary functions for unit testing
/**
 * @brief Transform the rhs (containing multiple vectors sequentially) to block form
 * such that the layout is nv of first vector -> nv of second vector -> nv of first -> ...
 * every nv elements should contain info about 1 lattice site.
 * */
template<typename T>
void transform_vector_block(std::complex<T>* v, int len, int nv, const int blen){
  auto *tmp = new std::complex<T>[blen * len];

  for(int i = 0; i < len/nv; i++){
    for(int j = 0; j < blen; j++){
      for (int k = 0; k < nv; ++k) {
        tmp[i*blen*nv+j*nv+k] = v[j*len+i*nv+k];
      }
    }
  }

  for(int i = 0; i < blen*len; i++){
    v[i]=tmp[i];
  }

  delete[] tmp;
}

/**
 * @brief Transform the block layout back to normal rhs layout
 * */
template<typename T>
void separate_vectors_in_block(std::complex<T>* v, int len, int nv, const int blen){
  auto *tmp = new std::complex<T>[blen*len];

  for (int i = 0; i < len/nv; i++) {
    for(int j = 0; j < blen; j++){
      for (int k = 0; k < nv; k++) {
        tmp[j*len+i*nv+k] = v[i*blen*nv+j*nv+k];
      }
    }
  }

  for(int i = 0; i < blen*len; i++){
    v[i]=tmp[i];
  }

  delete[] tmp;
}

#endif // DDALPHAAMG_DIRAC_BLOCK_VPSITE_H
