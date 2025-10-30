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

/** @file top_Level.h
 * Top Level functions and drivers.
 *
 * Contains the otermost functions and drivers of DDalphaAMG.
 */

#ifndef TOP_LEVEL_HEADER
#define TOP_LEVEL_HEADER

struct Thread;

/**
 * Defines the right hand side (rhs) of the system.
 *
 * Sets the rhs vector with the information contained in the object Level.
 * such as the lenght of the vector and type of entries (ones, random,...)
 *
 * @param rhs a vector for the rhs.
 * @param l the Level of the rhs.
 * @param threading a Thread struct.
 */
void rhs_define(vector_double rhs, Level *l);

/**
 * Driver for Wilson system solver.
 *
 * Solves the system using the set solver (fgmres,CGN,...)
 * If the WILSON_BENCHMARK is activated, it solves the system 100 times
 * then prints the average and the minimum solve times.
 *
 * @param solution a vector where the solution will be stored.
 * @param source a vector with the righ hand side (rhs) of the system.
 * @param l the Level where the solves are made.
 * @param threading a Thread struct.
 * @return The number of iterations when solving.
 */
int wilson_driver(vector_double solution, vector_double source, Level *l);

/**
 * Solves the system.
 *
 * Solves the system calling the wilson_driver.
 * TODO: if evaluation.... the rhs is set with vector_double_define_random.

 * @param solution a vector where the solution will be stored.
 * @param source a vector with the righ hand side (rhs) of the system.
 * @param l the Level_struct of the Level where the solves are made.
 * @param threading a Thread struct.
 * @return void.
 */
void solve(vector_double solution, vector_double source, Level *l);

// TODO: Document
void solve_driver(Level *l);

#endif
