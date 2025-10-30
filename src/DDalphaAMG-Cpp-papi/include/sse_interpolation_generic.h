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

#ifndef SSE_INTERPOLATION_double_HEADER
#define SSE_INTERPOLATION_double_HEADER

#ifdef SSE
void interpolation_double_alloc(Level *l);
void interpolation_double_free(Level *l);
void interpolation_double_dummy_alloc(Level *l);
void interpolation_double_dummy_free(Level *l);

void interpolate_double(vector_double phi, vector_double phi_c, Level *l, Thread *threading);
void interpolate3_double(vector_double phi, vector_double phi_c, Level *l, Thread *threading);
void restrict_double(vector_double phi_c, vector_double phi, Level *l, Thread *threading);
#endif

#endif