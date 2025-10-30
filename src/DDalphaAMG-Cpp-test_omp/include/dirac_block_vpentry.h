//
// Created by shiting on 2024-11-07.
//

#ifndef DDALPHAAMG_DIRAC_BLOCK_VPENTRY_H
#define DDALPHAAMG_DIRAC_BLOCK_VPENTRY_H
#include "omp.h"

#include "block.h"
#include "mpi.h"
#include "dirac_templated_generic.h"

static inline int get_neighbour(int t, int z, int y, int x, int dir, int *dim, int *proc){
  switch (dir) {
  case _T:
    return lex_index((proc[_T]>1)? t+1 : (t+1)%(dim[_T]-1), z, y, x, dim);
  case _Z:
    return lex_index(t, (proc[_Z]>1)? z+1 : (z+1)%(dim[_Z]-1), y, x, dim);
  case _Y:
    return lex_index(t, z, (proc[_Y]>1)? y+1 : (y+1)%(dim[_Y]-1), x, dim);
  case _X:
    return lex_index(t, z, y, (proc[_X]>1)? x+1 : (x+1)%(dim[_X]-1), dim);
  default:
    error0("invalid direction used.");
    exit(-1);
  }
}

template<typename T>
static inline void clover_per_entry(std::complex<T> *eta, std::complex<T> *phi, operator_struct<T> *op,
                                int end, Level *l, const int blen) {

  if (l->global->csw == 0.0) {
    std::complex<T> *clover = op->clover;
#pragma omp parallel for shared(end, phi, eta, clover)
    for (int i = 0; i < end; i++) {
      for (int j = 0; j < blen; j++){
        eta[i*blen + j] = phi[i*blen + j] * clover[i];
      }
    }

  } else {
    std::complex<T> *clover = op->clover;
#pragma omp parallel for shared(end, phi, eta, clover)
    for (int i = 0; i < end/12; i++) {
      site_clover_block(&eta[i*blen*12], &phi[i*blen*12], &clover[i*42], blen);
    }

  }
}

template<typename T>
static inline void pi_minus_per_entry(std::complex<T>* phi, std::complex<T>* prnT, std::complex<T>* prnZ,
                                      std::complex<T>* prnY, std::complex<T>* prnX,
                                      int start, int end, const int blen){
#pragma omp parallel for shared(start, end, phi, prnT, prnZ, prnY, prnX)
  for (int i = start/12; i < end/12; i++) {
    prp_T_block(&prnT[i*blen*6], &phi[i*blen*12], blen);
    prp_Z_block(&prnZ[i*blen*6], &phi[i*blen*12], blen);
    prp_Y_block(&prnY[i*blen*6], &phi[i*blen*12], blen);
    prp_X_block(&prnX[i*blen*6], &phi[i*blen*12], blen);
  }
}

template<typename T>
static inline void pi_plus_per_entry(std::complex<T>* phi, std::complex<T>* D, std::complex<T>* prpT,
                                     std::complex<T>* prpZ, std::complex<T>* prpY, std::complex<T>* prpX,
                                     const int* neighbor, int start, int end, const int blen){
  //auto *pbuf = new std::complex<T>[6*blen];

#pragma omp parallel for shared(start, end, phi, prpT, prpZ, prpY, prpX, D, neighbor)
  for (int i = start/12; i < end/12; i++) {
    /// T dir
    auto *pbuf = new std::complex<T>[6*blen];
    prn_T_block(pbuf, &phi[i*blen*12], blen); /// (I4 + Gamma_mu) kronecker I3 * psi(x)
    mvmh_block(&prpT[blen*6*neighbor[4*i]], &D[36*i], pbuf, blen); /// apply U_mu^H(x), first half
    mvmh_block(&prpT[blen*6*neighbor[4*i] + 3*blen], &D[36*i], &pbuf[3*blen], blen); /// second half

    /// Z dir
    prn_Z_block(pbuf, &phi[i * blen * 12], blen);
    mvmh_block(&prpZ[blen * 6 * neighbor[4 * i + 1]], &D[36 * i + 9], pbuf, blen);
    mvmh_block(&prpZ[blen * 6 * neighbor[4 * i + 1] + 3 * blen], &D[36 * i + 9],
               &pbuf[3 * blen], blen);

    /// Y dir
    prn_Y_block(pbuf, &phi[i * blen * 12], blen);
    mvmh_block(&prpY[blen * 6 * neighbor[4 * i + 2]], &D[36 * i + 18], pbuf, blen);
    mvmh_block(&prpY[blen * 6 * neighbor[4 * i + 2] + 3 * blen], &D[36 * i + 18],
               &pbuf[3 * blen], blen);

    /// X dir
    prn_X_block(pbuf, &phi[i*blen*12], blen);
    mvmh_block(&prpX[blen*6*neighbor[4*i+3]], &D[36*i+27], pbuf, blen);
    mvmh_block(&prpX[blen*6*neighbor[4*i+3] + 3*blen], &D[36*i+27], &pbuf[3*blen], blen);

    delete[] pbuf;
  }

 // delete[] pbuf;
}

template<typename T>
static inline void apply_U_per_entry(std::complex<T>* eta, std::complex<T>* D, std::complex<T>* prnT,
                                     std::complex<T>* prnZ, std::complex<T>* prnY, std::complex<T>* prnX,
                                     const int* neighbor, int start, int end, const int blen){

#pragma omp parallel for shared(start, end, eta, prnT, prnZ, prnY, prnX, D, neighbor)
  for (int i = start/12; i < end/12; i++) {
    auto *pbuf = new std::complex<T>[6*blen];
    /// T dir
    mvm_block(pbuf, &D[36 * i], &prnT[blen * 6 * neighbor[4 * i]],
              blen); /// apply U_mu(x) to prn, first half of x
    mvm_block(&pbuf[3 * blen], &D[36 * i], &prnT[blen * 6 * neighbor[4 * i] + 3 * blen],
              blen); /// second half of x
    pbp_su3_T_block(pbuf, &eta[i * blen * 12], blen);

    /// Z dir
    mvm_block(pbuf, &D[36 * i + 9], &prnZ[blen * 6 * neighbor[4 * i + 1]], blen);
    mvm_block(&pbuf[3 * blen], &D[36 * i + 9], &prnZ[blen * 6 * neighbor[4 * i + 1] + 3 * blen],
              blen);
    pbp_su3_Z_block(pbuf, &eta[i * blen * 12], blen);

    /// Y dir
    mvm_block(pbuf, &D[36 * i + 18], &prnY[blen * 6 * neighbor[4 * i + 2]], blen);
    mvm_block(&pbuf[3 * blen], &D[36 * i + 18],
              &prnY[blen * 6 * neighbor[4 * i + 2] + 3 * blen], blen);
    pbp_su3_Y_block(pbuf, &eta[i * blen * 12], blen);

    /// X dir
    mvm_block(pbuf, &D[36 * i + 27], &prnX[blen * 6 * neighbor[4 * i + 3]], blen);
    mvm_block(&pbuf[3 * blen], &D[36 * i + 27],
              &prnX[blen * 6 * neighbor[4 * i + 3] + 3 * blen], blen);
    pbp_su3_X_block(pbuf, &eta[i * blen * 12], blen);

    delete[] pbuf;
  }

}

template<typename T>
static inline void lift_pi_plus_per_entry(std::complex<T>* eta, std::complex<T>* prpT, std::complex<T>* prpZ,
                                          std::complex<T>* prpY, std::complex<T>* prpX,  int start,
                                          int end, const int blen){
#pragma omp parallel for shared(start, end, prpT, prpZ, prpY, prpX, eta)
  for (int i = start/12; i < end/12; i++) {
    pbn_su3_T_block(&prpT[i*blen*6], &eta[i*blen*12], blen);
    pbn_su3_Z_block(&prpZ[i*blen*6], &eta[i*blen*12], blen);
    pbn_su3_Y_block(&prpY[i*blen*6], &eta[i*blen*12], blen);
    pbn_su3_X_block(&prpX[i*blen*6], &eta[i*blen*12], blen);
  }
}

template<typename T>
void ghost_sendrecv_per_entry(std::complex<T> *phi, const int mu, const int dir, comm_struct<T> *c,
                              const int amount, Level *l, const int blen) {

  Global *global = l->global;

  // does not allow sending in both directions at the same time
  if (l->global_splitting[mu] > 1) {

    int *table = nullptr, mu_dir = 2 * mu - MIN(dir, 0), offset = c->offset,// offset = 6
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
        /// receive from -mu neighbour
        MPI_Irecv(buffer, length[1], mpi_type_handler().get_datatype(phi), l->neighbor_rank[2 * mu + 1], 2 * mu,
                  global->comm_cart, &(c->rreqs[2 * mu]));
      }
      if (length[0] > 0) {
        /// send to +mu neighbour
        MPI_Isend(phi_pt, length[0], mpi_type_handler().get_datatype(phi), l->neighbor_rank[2 * mu], 2 * mu,
                  global->comm_cart, &(c->sreqs[2 * mu]));
      }

    } else if (dir == -1) {
      // data to be communicated is stored on the sites touching the boundary in -mu direction
      // this data is gathered in a buffer in the correct ordering
      // which is required on the boundary of the vector phi
      int num_boundary_sites = length[1] / offset/ blen;
      //table stores the -mu boundary
      if (amount == _ODD_SITES){
        table = c->boundary_table[2 * mu + 1] + c->num_even_boundary_sites[mu_dir];
      } else {
        table = c->boundary_table[2 * mu + 1];
      }
#pragma omp parallel for shared(num_boundary_sites, offset, buffer, phi, table)
      for (int j = 0; j < num_boundary_sites; j++) {
        for (int i = 0; i < offset; i++) {
          for (int k = 0; k < blen; k++) {
            buffer[j*blen*offset + i*blen + k] = phi[table[j]*blen*offset + i*blen + k];
          }
        }
      }

      buffer = (std::complex<T>*)c->buffer[mu_dir];
      phi_pt = phi + comm_start;

      if (length[0] > 0) {
        /// receive from +mu neighbour
        MPI_Irecv(phi_pt, length[0], mpi_type_handler().get_datatype(phi), l->neighbor_rank[2 * mu], 2 * mu + 1,
                  global->comm_cart, &(c->rreqs[2 * mu + 1]));

      }
      if (length[1] > 0) {
        /// send to -mu neighbour
        MPI_Isend(buffer, length[1], mpi_type_handler().get_datatype(phi), l->neighbor_rank[2 * mu + 1], 2 * mu + 1,
                  global->comm_cart, &(c->sreqs[2 * mu + 1]));
      }

    } else
      ASSERT(dir == 1 || dir == -1);
  }
}

template<typename T>
void ghost_wait_per_entry(std::complex<T> *phi, const int mu, const int dir, comm_struct<T> *c,
                          const int amount, Level *l, const int blen) {

  if (l->global_splitting[mu] > 1) {
    int mu_dir = 2 * mu - MIN(dir, 0);
    int *table, offset = c->offset, length[2] = {0, 0};
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
#pragma omp parallel for shared(num_boundary_sites, offset, phi, table, buffer)
        for (int j = 0; j < num_boundary_sites; j++) {
          for (int i = 0; i < offset; i++) {
            for (int k = 0; k < blen; k++) {
              phi[table[j]*blen*offset + i*blen + k] = buffer[j*blen*offset + i*blen + k];
            }
          }
        }
      } else {
#pragma omp parallel for shared(num_boundary_sites, offset, phi, table, buffer)
        for (int j = 0; j < num_boundary_sites; j++) {
          for (int i = 0; i < offset; i++) {
            for (int k = 0; k < blen; k++) {
              phi[table[j]*blen*offset + i*blen + k] += buffer[j*blen*offset + i*blen + k];
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
void d_plus_clover_block(BlockVecPerEntry<T> *bout, BlockVecPerEntry<T> *bv, operator_struct<T> *op,
                         Level *l, const int blen){
  int n = l->num_inner_lattice_sites, *neighbor = op->neighbor_table, nv = l->num_lattice_site_var;
  int end =  nv * n;
  std::complex<T> *phi = bv->v;
  std::complex<T> *eta = bout->v;
  //omp_set_num_threads(2);

  auto *prnT = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];
  auto *prnZ = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];
  auto *prnY = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];
  auto *prnX = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];

  auto *prpT = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];
  auto *prpZ = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];
  auto *prpY = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];
  auto *prpX = new std::complex<T>[blen * nv/2 * l->num_lattice_sites];

  clover_per_entry(eta, phi, op, end, l, blen);

  pi_minus_per_entry(phi, prnT, prnZ, prnY, prnX, 0, end, blen);

  ///send prn to mu-1 direction, receive corresponding prn from mu+1 neighbor
  ghost_sendrecv_per_entry(prnT, _T, -1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_sendrecv_per_entry(prnZ, _Z, -1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_sendrecv_per_entry(prnY, _Y, -1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_sendrecv_per_entry(prnX, _X, -1, &(op->c), _FULL_SYSTEM, l, blen);

  pi_plus_per_entry(phi, op->D, prpT, prpZ, prpY, prpX, neighbor, 0, end, blen);

  /// send prp to mu+1 neighbour, receive corresponding prp from mu-1 neighbor
  ghost_sendrecv_per_entry(prpT, _T, +1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_sendrecv_per_entry(prpZ, _Z, +1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_sendrecv_per_entry(prpY, _Y, +1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_sendrecv_per_entry(prpX, _X, +1, &(op->c), _FULL_SYSTEM, l, blen);

  /// wait for communication of prn
  ghost_wait_per_entry(prnT, _T, -1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_wait_per_entry(prnZ, _Z, -1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_wait_per_entry(prnY, _Y, -1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_wait_per_entry(prnX, _X, -1, &(op->c), _FULL_SYSTEM, l, blen);

  apply_U_per_entry(eta, op->D, prnT, prnZ, prnY, prnX, neighbor, 0, end, blen);

  /// wait for communication of prp
  ghost_wait_per_entry(prpT, _T, +1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_wait_per_entry(prpZ, _Z, +1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_wait_per_entry(prpY, _Y, +1, &(op->c), _FULL_SYSTEM, l, blen);
  ghost_wait_per_entry(prpX, _X, +1, &(op->c), _FULL_SYSTEM, l, blen);

  lift_pi_plus_per_entry(eta, prpT, prpZ, prpY, prpX, 0, end, blen);

  delete[] prnT;
  delete[] prnZ;
  delete[] prnY;
  delete[] prnX;

  delete[] prpT;
  delete[] prpZ;
  delete[] prpY;
  delete[] prpX;
}

#endif // DDALPHAAMG_DIRAC_BLOCK_VPENTRY_H
