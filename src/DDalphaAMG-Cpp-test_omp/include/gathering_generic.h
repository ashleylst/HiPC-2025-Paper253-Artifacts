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

#ifndef GATHERING_double_HEADER
#define GATHERING_double_HEADER

void gathering_double_next_level_init(gathering_double_struct *gs, Level *l);
void gathering_double_setup(gathering_double_struct *gs, Level *l);

void gathering_double_init(gathering_double_struct *gs, Level *l);
void gathering_double_setup_new(gathering_double_struct *gs, Level *l);

void gathering_double_free(gathering_double_struct *gs, Level *l);

void conf_double_gather(operator_struct<double> *out, operator_struct<double> *in, Level *l);
void vector_double_gather(vector_double gath, vector_double dist, Level *l);
void vector_double_distribute(vector_double dist, vector_double gath, Level *l);

void distribution_double_next_level_test(Level *l);

#endif
