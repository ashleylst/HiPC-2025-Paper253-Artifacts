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

#include "main.h"

void prof_double_init(Level *l) {

  /*********************************************************************************
   * Initializes the profiling struct by specifying the name of every entry.
   *********************************************************************************/

  if (l != nullptr) {
    for (int i = 0; i < _NUM_PROF; i++) {
      l->prof_double.time[i] = 0.0;
      l->prof_double.count[i] = 0.0;
      l->prof_double.flop[i] = 0.0;
    }

    double level_ratio = 1;
    for (int mu = 0; mu < 4; mu++)
      level_ratio *= (double)l->global->global_lattice[l->depth][mu] /
                     (double)l->global->global_lattice[0][mu];

    sprintf(l->prof_double.name[_GIP], "global inner product, double");
    l->prof_double.flop[_GIP] = level_ratio * l->num_lattice_site_var * 8.0;
    sprintf(l->prof_double.name[_PIP], "process inner product, double");
    l->prof_double.flop[_PIP] = level_ratio * l->num_lattice_site_var * 8.0;
    sprintf(l->prof_double.name[_LA2], "2 flop vector operations, double");
    l->prof_double.flop[_LA2] = level_ratio * l->num_lattice_site_var * 2.0;
    sprintf(l->prof_double.name[_LA6], "6 flop vector operations, double");
    l->prof_double.flop[_LA6] = level_ratio * l->num_lattice_site_var * 6.0;
    sprintf(l->prof_double.name[_LA8], "8 flop vector operations, double");
    l->prof_double.flop[_LA8] = level_ratio * l->num_lattice_site_var * 8.0;
    sprintf(l->prof_double.name[_LA], "other vector operations, double");
    sprintf(l->prof_double.name[_GRAM_SCHMIDT], "Gram-Schmidt, double");
    sprintf(l->prof_double.name[_GRAM_SCHMIDT_ON_AGGREGATES], "Gram-Schmidt on aggregates, double");
    sprintf(l->prof_double.name[_CPY], "copy operations, double");
    sprintf(l->prof_double.name[_SET], "set value operations, double");
    sprintf(l->prof_double.name[_PR], "interpolation and restriction, double");
    l->prof_double.flop[_PR] =
        level_ratio * l->num_lattice_site_var * 8.0 * (l->num_lattice_site_var / 2);
    sprintf(l->prof_double.name[_SC], "self coupling, double");
    l->prof_double.flop[_SC] =
        (l->depth == 0) ? 552.0 : level_ratio * SQUARE(l->num_lattice_site_var) * 8.0;
    sprintf(l->prof_double.name[_NC], "neighbor coupling, double");
    l->prof_double.flop[_NC] =
        (l->depth == 0) ? 1368.0 : level_ratio * 8.0 * SQUARE(l->num_lattice_site_var) * 8.0;
    sprintf(l->prof_double.name[_SM], "smoother, double");
    double ncflops = l->prof_double.flop[_SC];
    for (int mu = 0; mu < 4; mu++)
      ncflops += (l->prof_double.flop[_NC] / 4.0) *
                 ((double)(l->block_lattice[mu] - 1) / (double)l->block_lattice[mu]);
    l->prof_double.flop[_SM] =
        ncflops * (double)(l->global->odd_even ? l->block_iter + 1 : l->block_iter);
    l->prof_double.flop[_SM] += (l->prof_double.flop[_NC] + l->prof_double.flop[_SC]);
    sprintf(l->prof_double.name[_OP_COMM], "operator comm init, double");
    sprintf(l->prof_double.name[_OP_IDLE], "operator comm wait, double");
    sprintf(l->prof_double.name[_ALLR], "allreduces, double");
    sprintf(l->prof_double.name[_GD_COMM], "data re-distribution comm init, double");
    sprintf(l->prof_double.name[_GD_IDLE], "data re-distribution comm wait, double");
    sprintf(l->prof_double.name[_SM1], "smoother - pt 1, res no comm, double");
    sprintf(l->prof_double.name[_SM2], "smoother - pt 2, solve no comm, double");
    sprintf(l->prof_double.name[_SM3], "smoother - pt 3, res comm, double");
    sprintf(l->prof_double.name[_SM4], "smoother - pt 4, solve comm, double");

    sprintf(l->prof_double.name[_SMALL1], "Hessenberg: qr update double");
    sprintf(l->prof_double.name[_SMALL2], "Hessenberg: bkwd subst double");
  }
}

double prof_double_print(Level *l) {
  double flop = 0;
  for (int i = 0; i < _NUM_PROF; i++)
    if (l->prof_double.count[i] > 0) {
      if (l->prof_double.count[i] > 9999999)
        printf0("| %37s: %8.2le(%7.1le) |\n", l->prof_double.name[i], l->prof_double.time[i],
                l->prof_double.count[i]);
      else
        printf0("| %37s: %8.2le(%7d) |\n", l->prof_double.name[i], l->prof_double.time[i],
                (int)l->prof_double.count[i]);
      flop += (double)l->prof_double.count[i] * l->prof_double.flop[i];
    }
  return flop;
}
