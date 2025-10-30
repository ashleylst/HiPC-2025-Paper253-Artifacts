/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel Simon Heybrock, Simone Bacchio, Bjoern
 * Leder.
 *
 * This file is part of the DDalphaAMG solver library.
 *
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The DDalphaAMG solver library is distributed in the hope that it will be useful
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 *
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 *
 */

#ifndef ODDEVEN_double_HEADER
#define ODDEVEN_double_HEADER

void selfcoupling_cholesky_decomposition_double(const config_double output, config_double input);

void hopping_term_double(vector_double eta, vector_double phi, operator_struct<double> *op,
                         const int amount, Level *l);

void oddeven_setup_double(operator_struct<double> *in, Level *l);
void oddeven_free_double(Level *l);

void oddeven_to_serial_double(vector_double out, vector_double in, Level *l);
void serial_to_oddeven_double(vector_double out, vector_double in, Level *l);

void oddeven_to_block_double(vector_double out, vector_double in, Level *l);
void block_to_oddeven_double(vector_double out, vector_double in, Level *l);

void block_hopping_term_double(vector_double eta, vector_double phi, int start, int amount,
                               Schwarz *s, Level *l);
void block_n_hopping_term_double(vector_double eta, vector_double phi, int start, int amount,
                                 Schwarz *s, Level *l);
void block_diag_oo_inv_double(vector_double eta, vector_double phi, int start, Schwarz *s,
                              Level *l);
void block_diag_oo_double(vector_double eta, vector_double phi, int start, Schwarz *s, Level *l);
void block_diag_ee_double(vector_double eta, vector_double phi, int start, Schwarz *s, Level *l);

void apply_schur_complement_double(vector_double out, vector_double in, operator_struct<double> *op,
                                   Level *l);
void solve_oddeven_double(Gmres *p, operator_struct<double> *op, Level *l);
void g5D_apply_schur_complement_double(vector_double out, vector_double in,
                                       operator_struct<double> *op, Level *l);
void g5D_solve_oddeven_double(Gmres *p, operator_struct<double> *op, Level *l);

void schwarz_double_oddeven_setup(Schwarz *s, Level *l);

void apply_block_schur_complement_double(vector_double out, vector_double in, int start, Schwarz *s,
                                         Level *l);
void block_solve_oddeven_double(vector_double phi, vector_double r, vector_double latest_iter,
                                int start, Schwarz *s, Level *l);

void oddeven_double_test(Level *l);
void block_oddeven_double_test(Level *l);

#endif
