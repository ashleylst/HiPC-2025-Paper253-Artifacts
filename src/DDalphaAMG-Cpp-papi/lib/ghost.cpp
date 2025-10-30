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

void neighbor_define(Level *l) {

  /*********************************************************************************
   * Determines the rank of the (cartesian) neighbors of a process and stores them
   * in l->neighbor_rank.
   *********************************************************************************/

  int neighbor_coords[4];
  Global *global = l->global;

  for (int mu = 0; mu < 4; mu++) {
    neighbor_coords[mu] = global->my_coords[mu];
  }

  for (int mu = 0; mu < 4; mu++) {
    neighbor_coords[mu] += l->comm_offset[mu];
    global->Cart_rank(global->comm_cart, neighbor_coords, &(l->neighbor_rank[2 * mu]));
    neighbor_coords[mu] -= 2 * l->comm_offset[mu];

    global->Cart_rank(global->comm_cart, neighbor_coords, &(l->neighbor_rank[2 * mu + 1]));
    neighbor_coords[mu] += l->comm_offset[mu];
  }
}

void neighbor_define(Global *g) {

  /*********************************************************************************
   * Determines the rank of the (cartesian) neighbors of a process and stores them
   * in g->neighbor_rank.
   *********************************************************************************/

  int neighbor_coords[4];

  for (int mu = 0; mu < 4; mu++) {
    neighbor_coords[mu] = g->my_coords[mu];
  }

  for (int mu = 0; mu < 4; mu++) {
    neighbor_coords[mu] += g->comm_offset[mu];
    g->Cart_rank(g->comm_cart, neighbor_coords, &(g->neighbor_rank[2 * mu]));
    neighbor_coords[mu] -= 2 * g->comm_offset[mu];

    g->Cart_rank(g->comm_cart, neighbor_coords, &(g->neighbor_rank[2 * mu + 1]));
    neighbor_coords[mu] += g->comm_offset[mu];
  }
}

void cart_define(MPI_Comm comm, Level *l) {

  /*********************************************************************************
   * Creates a cartesian MPI Communicator and assignes 4D Coordinates to every
   * process. Also defines the neighbor processes for every process.
   *********************************************************************************/

  int num_processes;

  Global *global = l->global;

  MPI_Comm_size(comm, &num_processes);
  if (num_processes != global->num_processes) {
    error0("Number of processes has to be %d, but is %d\n", global->num_processes, num_processes);
  }
  MPI_Cart_create(comm, 4, l->global_splitting, l->periodic_bc, 1, &(global->comm_cart));
  MPI_Comm_rank(global->comm_cart, &(global->my_rank));
  global->Cart_coords(global->comm_cart, global->my_rank, 4, global->my_coords);

  for (int mu = 0; mu < 4; mu++) {
    l->num_processes_dir[mu] = l->global_lattice[mu] / l->local_lattice[mu];
    l->comm_offset[mu] = 1;
  }

  neighbor_define(l);
  MPI_Comm_group(global->comm_cart, &(global->global_comm_group));
}

void cart_define(MPI_Comm comm, Global *g) {

  /*********************************************************************************
   * Creates a cartesian MPI Communicator and assignes 4D Coordinates to every
   * process. Also defines the neighbor processes for every process.
   *********************************************************************************/

  int num_processes;

  MPI_Comm_size(comm, &num_processes);
  if (num_processes != g->num_processes) {
    error0("Number of processes has to be %d, but is %d\n", g->num_processes, num_processes);
  }

  MPI_Cart_create(comm, 4, g->process_grid, g->periodic_bc, 1, &(g->comm_cart));
  MPI_Comm_rank(g->comm_cart, &(g->my_rank));
  g->Cart_coords(g->comm_cart, g->my_rank, 4, g->my_coords);

  for (int mu = 0; mu < 4; mu++) {
    g->num_processes_dir[mu] = g->global_lattice[0][mu] / g->local_lattice[0][mu];
    g->comm_offset[mu] = 1;
  }

  neighbor_define(g);
  MPI_Comm_group(g->comm_cart, &(g->global_comm_group));
}

void cart_validate(MPI_Comm comm, Level *l) {

  /*********************************************************************************
   * Validates a cartesian MPI Communicator and assignes 4D Coordinates to every
   * process. Also defines the neighbor processes for every process.
   *********************************************************************************/

  int num_processes;
  Global *global = l->global;

  MPI_Comm_size(comm, &num_processes);
  if (num_processes != global->num_processes) {
    error0("Number of processes has to be %d, but is %d\n", global->num_processes, num_processes);
  }
  // MPI_Cart_create( comm, 4, l->global_splitting, l->periodic_bc, 1, &(global->comm_cart) );
  // TODO: Validation procedure as input of line above
  MPI_Comm_dup(comm, &(global->comm_cart));
  MPI_Comm_rank(global->comm_cart, &(global->my_rank));
  global->Cart_coords(global->comm_cart, global->my_rank, 4, global->my_coords);

  for (int mu = 0; mu < 4; mu++) {
    l->num_processes_dir[mu] = l->global_lattice[mu] / l->local_lattice[mu];
    l->comm_offset[mu] = 1;
  }

  neighbor_define(l);
  MPI_Comm_group(global->comm_cart, &(global->global_comm_group));
}

void cart_free(Level *l) {

  /*********************************************************************************
   * Deletes the communicator and the according group.
   *********************************************************************************/

  Global *global = l->global;
  MPI_Group_free(&(global->global_comm_group));
  MPI_Comm_free(&(global->comm_cart));
}
