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

#ifndef COARSE_ODDEVEN_double_HEADER
#define COARSE_ODDEVEN_double_HEADER

void coarse_oddeven_alloc_double(Level *l);

void coarse_oddeven_setup_double(operator_struct<double> *in, int reorder, Level *l);
void coarse_oddeven_double_set_self_couplings(Level *l);

void coarse_oddeven_free_double(Level *l);

void coarse_solve_odd_even_double(Gmres *p, operator_struct<double> *op, Level *l);
void coarse_apply_schur_complement_double(vector_double out, vector_double in,
                                          operator_struct<double> *op, Level *l);
void g5D_coarse_solve_odd_even_double(Gmres *p, operator_struct<double> *op, Level *l);
void g5D_coarse_apply_schur_complement_double(vector_double out, vector_double in,
                                              operator_struct<double> *op, Level *l);

void coarse_hopping_term_double(vector_double out, vector_double in, operator_struct<double> *op,
                                const int amount, Level *l);
void coarse_n_hopping_term_double(vector_double out, vector_double in, operator_struct<double> *op,
                                  const int amount, Level *l);
void coarse_hopping_term_double_vectorized(vector_double out, vector_double in,
                                           operator_struct<double> *op, const int amount, Level *l);
void coarse_pn_hopping_term_double_vectorized(vector_double out, vector_double in,
                                              operator_struct<double> *op, const int amount,
                                              Level *l, int sign);
void coarse_n_hopping_term_double_vectorized(vector_double out, vector_double in,
                                             operator_struct<double> *op, const int amount,
                                             Level *l);

void coarse_odd_even_double_test(vector_double c4, vector_double c1, Level *l);

void coarse_diag_ee_double(vector_double y, vector_double x, operator_struct<double> *op, Level *l);

void coarse_diag_oo_double(vector_double y, vector_double x, operator_struct<double> *op, Level *l);

void coarse_diag_oo_inv_double(vector_double y, vector_double x, operator_struct<double> *op,
                               Level *l);

void coarse_perform_fwd_bwd_subs_double(vector_double x, vector_double b, config_double A,
                                        Level *l);

#endif
