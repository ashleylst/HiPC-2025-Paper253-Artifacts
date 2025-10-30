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

#ifndef GHOST_double_HEADER
#define GHOST_double_HEADER

void negative_sendrecv_double(vector_double phi, const int mu, comm_struct<double> *c, Level *l);

// as negative_sendrecv_double, but for count vectors stored in phi in vector-fused data layout
// buffer must be big enough to hold the surface data for count vectors (in one direction)
void negative_sendrecv_double_vectorized(complex_double *phi, const int mu, comm_struct<double> *c,
                                         Level *l, int count, complex_double *buffer);
void negative_wait_double(const int mu, comm_struct<double> *c, Level *l);

void ghost_alloc_double(int buffer_size, comm_struct<double> *c, Level *l);
void ghost_alloc_double(int buffer_size, comm_struct<double> *c, Global *g);
void ghost_free_double(comm_struct<double> *c, Level *l);
void ghost_sendrecv_init_double(const int type, comm_struct<double> *c, Level *l);
void ghost_sendrecv_double(vector_double phi, const int mu, const int dir, comm_struct<double> *c,
                           const int amount, Level *l);
void ghost_wait_double(vector_double phi, const int mu, const int dir, comm_struct<double> *c,
                       const int amount, Level *l);

void ghost_update_double(vector_double phi, const int mu, const int dir, comm_struct<double> *c,
                         Level *l);
void ghost_update_wait_double(vector_double phi, const int mu, const int dir, comm_struct<double> *c,
                              Level *l);

template<typename T>
void ghost_free(comm_struct<T> *c, Level *l){
  int mu;

  for (mu = 0; mu < 4; mu++) {
    delete[] c->buffer[2 * mu];
    delete[] c->buffer[2 * mu + 1];
  }
}

template<typename T>
void ghost_alloc(int buffer_size, comm_struct<T> *c, Level *l) {

  int mu, nu, factor = 1;

  if (l->depth > 0) {
    c->offset = l->num_lattice_site_var;
  } else {
    c->offset = l->num_lattice_site_var / 2;
    if (l->global->method < 5)
      factor = 2;
  }

#ifdef HAVE_TM1p1
  factor *= 2;
#endif

  if (buffer_size <= 0) {
    c->comm_start[0] = c->offset * l->num_inner_lattice_sites;
    c->comm_start[1] = c->offset * l->num_inner_lattice_sites;
    for (mu = 0; mu < 4; mu++) {
      if (mu > 0) {
        c->comm_start[2 * mu] = c->comm_start[2 * (mu - 1)] + buffer_size;
        c->comm_start[2 * mu + 1] = c->comm_start[2 * (mu - 1) + 1] + buffer_size;
      }
      buffer_size = c->offset;
      for (nu = 0; nu < 4; nu++) {
        if (nu != mu) {
          buffer_size *= l->local_lattice[nu];
        }
      }
      c->length[2 * mu] = buffer_size;
      c->length[2 * mu + 1] = buffer_size;
      c->max_length[mu] = factor * buffer_size;
      c->buffer[2 * mu] = new std::complex<T>[factor * buffer_size];
      c->buffer[2 * mu + 1] = new std::complex<T>[factor * buffer_size];
      c->in_use[2 * mu] = 0;
      c->in_use[2 * mu + 1] = 0;
    }
  } else {
    for (mu = 0; mu < 4; mu++) {
      c->max_length[mu] = buffer_size;
#ifdef HAVE_TM1p1
      c->buffer[2 * mu] = new std::complex<T>[2 * buffer_size];
      c->buffer[2 * mu + 1] = new std::complex<T>[2 * buffer_size];
#else
      c->buffer[2 * mu] = new std::complex<T>[buffer_size];
      c->buffer[2 * mu + 1] = new std::complex<T>[buffer_size];
#endif
    }
  }

  if (l->vbuf_double[8] == nullptr) {
#ifdef HAVE_TM1p1
    MALLOC(l->vbuf_double[8], complex_double, 2 * l->vector_size);
#else
    MALLOC(l->vbuf_double[8], complex_double, l->vector_size);
#endif
  }
}

#ifdef PERS_COMMS
void pers_comms_init_double(vector_double phi, const int mu, const int dir,
                            const int pers_comms_id1, const int pers_comms_id2,
                            comm_struct<double> *c, const int amount, Level *l);
void pers_comms_open_double(Level *l);
void pers_comms_free_double(vector_double phi, const int mu, const int dir,
                            const int pers_comms_id1, const int pers_comms_id2,
                            comm_struct<double> *c, const int amount, Level *l);
void pers_comms_close_double(Level *l);
#endif

#endif
