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

#ifndef IO_HEADER
#define IO_HEADER
#ifdef HAVE_HDF5
#include "hdf5.h"
#define ROOTGROUPNAME "/rqcd"
#define EIGENMODEGROUPNAME "eigenmodes"

typedef struct Hdf5_fileinfo {
  char *filename;
  hid_t file_id, rootgroup_id, configgroup_id, eigenmodegroup_id, thiseigenmodegroup_id;
  unsigned int isOpen;
  int mode;
  double ioTime;
} Hdf5_fileinfo;

extern Hdf5_fileinfo h5info;
#endif
void byteswap(char *in);
void byteswap8(char *in);
void read_conf(double *input_data, char *input_name, double *conf_plaq, Level *l);
void read_conf(double *input_data, char *input_name, double *conf_plaq, Global *global);
void vector_io(double *phi, char *filename, const int mode, Level *l);
void vector_io_single_file(double *psi, double *lambda, char *filename, const int mode, int n,
                           char *vector_type, Level *l);
void d_dump(config_double D, Level *l);

// static inline int process_index(int t, int z, int y, int x, int ll[4]) {
//
//   if (g.num_processes > 1) {
//     int pcoords[4];
//     int rank;
//     pcoords[_T] = t / ll[_T];
//     pcoords[_Z] = z / ll[_Z];
//     pcoords[_Y] = y / ll[_Y];
//     pcoords[_X] = x / ll[_X];
//     g.Cart_rank(g.comm_cart, pcoords, &rank);
//     return rank;
//   } else {
//     return 0;
//   }
// }

static inline int process_index(int t, int z, int y, int x, int ll[4], Global *g) {

  if (g->num_processes == 1)
    return 0;

  int pcoords[4];
  int rank;

  pcoords[_T] = t / ll[_T];
  pcoords[_Z] = z / ll[_Z];
  pcoords[_Y] = y / ll[_Y];
  pcoords[_X] = x / ll[_X];
  g->Cart_rank(g->comm_cart, pcoords, &rank);

  return rank;
}

#ifdef HAVE_HDF5
unsigned int stepIntoEigenmode(int index);
unsigned int initFile(char *filename, const int mode, Level *l);
void closeFile();

#endif
#endif
