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

#include "DDalphaAMG.h"
#include "main.h"

void method_update(int setup_iter, Level *l) {

  Global *global = l->global;

  if (global->method > 0 && global->interpolation && global->num_levels > 1 && setup_iter > 0) {

    double t0 = 0, t1 = 0;

    global->in_setup = 1;
    t0 = MPI_Wtime();

    if (global->setup_m0 != global->m0) {
      //       m0_update((complex_double)global->setup_m0, l); // Why the cast to
      //       complex? (happend again several times below in this file)
      m0_update(global->setup_m0, l);
#ifdef HAVE_TM
    }
    if (global->setup_mu != global->mu) {
      tm_term_update(global->setup_mu, l);
      finalize_operator_update(l);
    } else if (global->setup_m0 != global->m0) {
#endif
      finalize_operator_update(l);
    }

    if (global->mixed_precision)
      iterative_double_setup(setup_iter, l);
    else
      iterative_double_setup(setup_iter, l);
    if (global->setup_m0 != global->m0) {
      m0_update(global->m0, l);
#ifdef HAVE_TM
    }
    if (global->setup_mu != global->mu) {
      tm_term_update(global->mu, l);
      finalize_operator_update(l);
    } else if (global->setup_m0 != global->m0) {
#endif
      finalize_operator_update(l);
    }

    t1 = MPI_Wtime();
    global->total_time = t1 - t0;
    printf0("\nperformed %d iterative setup steps\n", setup_iter);
    printf0("elapsed time: %lf seconds (%lf seconds on coarse grid)\n\n", t1 - t0,
            global->coarse_time);
  }

  global->in_setup = 0;

#ifdef DEBUG
  test_routine(l);
#endif
}

int read_parameter(void **save_at, const char *search_pattern, const char *read_format, int number,
                   FILE *read_from, int set_default) {

  /*********************************************************************************
   * Reads input parameters from provided inputfile.
   * - void **save_at: Points to variable, where parameters are stored.
   * - const char *search_pattern: Gives the name of parameter to be read.
   * - const char *read_format: Specifies datatype of parameter.
   * - int number: Specifies how many numbers need to be read.
   * - FILE *read_from: The inputfile.
   * - int set_default: Either _NO_DEFAULT_SET for no default value, or DEFAULT_SET,
   *   if a default value is assigned to this parameter.
   *   If set_default = _NO_DEFAULT_SET, then this parameter MUST be set in the
   *   inputfile.
   *
   * For examples, see lg_in(...) further below.
   *********************************************************************************/

  int i = 0, j, k, n = strlen(search_pattern), match = 0;
  char read_pattern[100000], *read_pattern_pt, buffer[50];
  var_table_entry e;

  if (read_from == nullptr) {
    if (!set_default)
      error0("FILE nullptr, unable to find string \"%s\" --- fatal error\n", search_pattern);
    else
      return match;
  }

  fseek(read_from, 0L, SEEK_SET);

  while (!match && fgets(read_pattern, 100000, read_from)) {

    k = strlen(read_pattern);
    /*
    j = 0;
    for ( i=0; i<k && !match; i++ ) {
      if ( search_pattern[j] == read_pattern[i] )
        j++;
      else
        j = 0;
      if ( j == n ) {
        match = 1;
      }
    }
    */ // replace it with a search just at the beginning of the line.
    if (k > n) {
      match = 1;
      i = 0;
      while (i < n && match) {
        if (search_pattern[i] != read_pattern[i])
          match = 0;
        i++;
      }
    }
  }

  read_pattern_pt = read_pattern + i;
  while (*read_pattern_pt == ' ')
    read_pattern_pt++;

  if (match) {
    if (strcmp(read_format, "%s") != 0) {
      e.pt = *save_at;
      for (i = 0; i < n - 1; i++) {
        e.name[i] = search_pattern[i];
      }
      e.name[n - 1] = '\0';
      if (strcmp(read_format, "%d") == 0) {
        // int
        for (j = 0; j < number; j++) {
          sscanf(read_pattern_pt, read_format, &(((int *)*save_at)[j]));
          sscanf(read_pattern_pt, "%s", buffer);
          read_pattern_pt += strlen(buffer);
          while (*read_pattern_pt == ' ')
            read_pattern_pt++;
        }
        sprintf(e.datatype, "int");
      } else {
        // double
        for (j = 0; j < number; j++) {
          sscanf(read_pattern_pt, read_format, &(((double *)*save_at)[j]));
          sscanf(read_pattern_pt, "%s", buffer);
          read_pattern_pt += strlen(buffer);
          while (*read_pattern_pt == ' ')
            read_pattern_pt++;
        }
        sprintf(e.datatype, "double");
      }

      //       if (number == 1) {
      //         var_table_insert(&(l->global->vt), e);
      //       }
    } else {
      // string
      sprintf(((char *)*save_at), "%s", read_pattern_pt);
      ((char *)*save_at)[strlen(read_pattern_pt) - 1] = '\0';
    }
  } else {
    if (!set_default)
      error0("unable to find string \"%s\" --- fatal error\n", search_pattern);
  }

  return match;
}

void set_DDalphaAMG_parameters(struct init *params, Level *l) {

  //   FILE *in = nullptr;
  //
  //   if (params->init_file != nullptr)
  //     ASSERT((in = fopen(params->init_file, "r")) != nullptr);
  //
  //   l->global->num_levels = params->number_of_levels;
  //   l->global->num_desired_levels = l->global->num_levels;
  //
  //   int ls = MAX(l->global->num_levels, 2);
  //   allocate_for_global_struct_after_read_global_info(ls);
  //
  //   set_global_info(params, l);
  //   read_solver_parameters(in, l);
  //   read_geometry_data(in, ls);
  //
  //   ls = MAX(l->global->num_levels, 2);
  //
  //   set_level_and_global_structs_according_to_global_struct(l);
  //
  //   read_testvector_io_data_if_necessary(in);
  //   read_evaluation_parameters_if_necessary(in);
  //   read_kcycle_data(in);
  //
  //   validate_parameters(ls, l);
  //
  //   if (params->init_file != nullptr)
  //     fclose(in);
}
