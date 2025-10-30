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

#ifndef DATA_LAYOUT_HEADER
#define DATA_LAYOUT_HEADER

void data_layout_init(Level *l);
void data_layout_n_flavours(int n, Level *l);
void define_eot(int *eot, int *N, Level *l);
void define_eo_bt(int **bt, int *eot, int *n_ebs, int *n_obs, int *n_bs, int *N, Level *l);
void define_nt_bt_tt(int *nt, int *backward_nt, int **bt, int *tt, int *it, int *dt, Level *l);
void define_nt_bt_tt(int *nt, int *backward_nt, int **bt, int *tt, int *it, int *dt, Global *g);

static inline int lex_index(int t, int z, int y, int x, int N[3]) {
  return x + N[_X] * (y + N[_Y] * (z + N[_Z] * t));
}

static inline int lex_mod_index(int t, int z, int y, int x, int N[4]) {
  return (x + N[_X]) % N[_X] +
         N[_X] *
             ((y + N[_Y]) % N[_Y] + N[_Y] * ((z + N[_Z]) % N[_Z] + N[_Z] * ((t + N[_T]) % N[_T])));
}

static inline int site_index(int t, int z, int y, int x, int N[3], int *index_table) {
  return index_table[lex_index(t, z, y, x, N)];
}

static inline int site_mod_index(int t, int z, int y, int x, int N[4], int *index_table) {
  return index_table[lex_mod_index(t, z, y, x, N)];
}

#endif
