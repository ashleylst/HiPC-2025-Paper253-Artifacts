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

void read_tv_from_file_double(Level *l) {

  char filename[STRINGLENGTH + 1], tv[] = "test vectors";
  vector_double tmp = new complex_double[l->inner_vector_size];
  Global *global = l->global;

  if (l->type == _FINE) {
    if (global->tv_io_single_file) {
      vector_io_single_file(nullptr, nullptr, global->tv_io_file_name, _READ, l->num_eig_vect, tv, l);
      re_setup_double(l);
    } else {
      for (int i = 0; i < l->num_eig_vect; i++) {
        error0("Fix -Wformat-overflow warning for next two lines before uncommenting!"); // TODO
        //         sprintf(filename, "%s.%02d", global->tv_io_file_name, i);
        //         printf0("%s.%02d\n", global->tv_io_file_name, i);
        vector_io((double *)tmp, filename, _READ, l);
        trans_double(l->is_double.test_vector[i], tmp, l->s_double.op.translation_table, l);
      }
      re_setup_double(l);
    }
  }
  delete[] tmp;
}

void re_setup_double(Level *l) {

  if (l->type != _COARSEST) {
    if (!l->idle) {
#ifdef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_double
      define_interpolation_double_operator(l->is_double.test_vector, l);
      gram_schmidt_on_aggregates_double_vectorized(l->is_double.operator, l->num_eig_vect, l);
      if (l->type != _FINE)
        gram_schmidt_on_aggregates_double_vectorized(l->is_double.operator, l->num_eig_vect, l);
      coarse_operator_double_setup_vectorized(l->is_double.operator, l);

#else
      for (int i = 0; i < l->num_eig_vect; i++) {
        vector_double_copy(l->is_double.interpolation[i], l->is_double.test_vector[i], 0,
                           l->inner_vector_size, l);
      }
      gram_schmidt_on_aggregates_double(l->is_double.interpolation, l->num_eig_vect, l);
      if (l->type != _FINE)
        gram_schmidt_on_aggregates_double(l->is_double.interpolation, l->num_eig_vect, l);
      define_interpolation_double_operator(l->is_double.interpolation, l);

      coarse_operator_double_setup(l->is_double.interpolation, l);
#endif
      conf_double_gather(&(l->next_level->s_double.op), &(l->next_level->op_double), l->next_level);

      if (!l->next_level->idle && l->next_level->type != _COARSEST) {

        schwarz_double_boundary_update(&(l->next_level->s_double), l->next_level);

        if (l->global->method >= 4 && l->odd_even) {
          coarse_oddeven_setup_double(&(l->next_level->s_double.op), _REORDER, l->next_level);
        } else {
          coarse_operator_double_set_couplings(&(l->next_level->s_double.op), l->next_level);
        }
      }
      if (!l->next_level->idle && l->next_level->type == _COARSEST && l->odd_even) {
        coarse_oddeven_setup_double(&(l->next_level->s_double.op), _NO_REORDERING, l->next_level);
      } else if (!l->next_level->idle && l->next_level->type == _COARSEST) {
        coarse_operator_double_set_couplings(&(l->next_level->s_double.op), l->next_level);
      }
      re_setup_double(l->next_level);
    }
  }
#if defined(POLYPREC) || defined(GCRODR) || defined(BLOCK_JACOBI)
  else {
    // this runs on currentLevel 0 only
#ifdef POLYPREC
    l->p_double.polyprec_double.update_lejas = 1;
    l->p_double.polyprec_double.preconditioner = nullptr;
#endif
#ifdef GCRODR
    l->p_double.gcrodr_double.update_CU = 1;
    l->p_double.gcrodr_double.upd_ctr = 0;
#endif
#ifdef BLOCK_JACOBI
    l->p_double.block_jacobi_double.local_p.polyprec_double.update_lejas = 1;
    l->p_double.block_jacobi_double.BJ_usable = 0;
#endif
  }
#endif
}

void set_kcycle_tol_double(double tol, Level *l) {

  if (!l->idle) {
    l->gmresKCycle.tolerance = tol;
#ifdef BLOCK_JACOBI
    if (l->type == _COARSEST)
      l->p_double.block_jacobi_double.local_p.tol = tol;
#endif
  }

  if (l->next_level->type != _COARSEST)
    set_kcycle_tol_double(tol, l->next_level);
}

void test_vector_double_update(int i, Level *l) {

  double norm;
  if (l->next_level->type != _COARSEST)
    test_vector_double_update(i, l->next_level);

  if (!l->idle && i<l->num_eig_vect) {
    norm = global_norm_double(l->vbuf_double[1], 0, l->inner_vector_size, l);
    vector_double_real_scale(l->is_double.test_vector[i], l->vbuf_double[1], 1.0 / norm, 0,
                             l->inner_vector_size, l);
  }
}

void inv_iter_inv_fcycle_double(int setup_iter, Level *l) {

  Global *global = l->global;

  complex_double *buffer = new complex_double[2 * l->num_eig_vect];

  if (l->type == _FINE)
    set_kcycle_tol_double(global->coarse_tol, l);

  if (!l->idle) {
    for (int j = 0; j < setup_iter; j++) {
      int pc = 0;
#ifdef DEBUG
      int pi = 1, pn = l->num_eig_vect * l->post_smooth_iter;
#endif

      if (global->print > 0)
        printf0("depth: %d, bootstrap step number %d...\n", l->depth, j + 1);
#ifdef DEBUG
      if (global->print > 0) {
        printf0("\033[0;42m\033[1;37m|");
        if (global->my_rank == 0)
          fflush(0);
      }
#endif

      gram_schmidt_double(l->is_double.test_vector, buffer, 0, l->num_eig_vect,
                          l); // NOTE: Can we drop this if we use .interpolation vectors instead?
      //       gram_schmidt_double(l->is_double.test_vector, buffer, 0, l->num_eig_vect, l);


      //NOTE: IMPORTANT: The vbuf in below for loop HAS TO BE l->vbuf_double[1(!)] in order to match with how the vcycle stores results on each level, i.e., also in l->vbuf_double[1]
      //TODO: Untangle this somehow, this is unnecessarily hard to maintain
      for (int i = 0; i < l->num_eig_vect; i++) {
        vcycle_double(l->vbuf_double[1], nullptr, l->is_double.test_vector[i], _NO_RES, l);
        double norm = global_norm_double(l->vbuf_double[1], 0, l->inner_vector_size, l);
        printf0("norm iterative TV[%d] = %.15lf\n", i, norm);
        test_vector_double_update(i, l);

        pc += l->post_smooth_iter;
#ifdef DEBUG

        if (pc >= (int)((0.2 * pi) * pn)) {
          if (global->print > 0) {
            printf0("%4d%% |", 20 * pi);
            if (global->my_rank == 0)
              fflush(0);
          }
          pi++;
        }

#endif
      }

#ifdef DEBUG

      if (global->print > 0)
        printf0("\033[0m\n");

#endif
      re_setup_double(l);
      if (l->type == _FINE && l->next_level->type != _COARSEST) {
        inv_iter_inv_fcycle_double(
            MAX(1, round(((double)(j + 1) * l->next_level->setup_iter) / ((double)setup_iter))),
            l->next_level);
      }
    }
    if (l->type == _INTERMEDIATE && l->next_level->type != _COARSEST) {
      inv_iter_inv_fcycle_double(
          MAX(1, round((double)(l->next_level->setup_iter * setup_iter) / ((double)l->setup_iter))),
          l->next_level);
    }
  }

  delete[] buffer;

  if (l->type == _FINE) {

    set_kcycle_tol_double(global->kcycle_tol, l);
  }
}

void testvector_analysis_double(vector_double *test_vectors, Level *l) {
  // #ifdef TESTVECTOR_ANALYSIS
  //   double mu;
  //   complex_double lambda;
  //
  //   //
  //   if (l->depth == 0) {
  //     printf0("--------------------------------------- depth: %d "
  //             "----------------------------------------\n",
  //             l->depth);
  //     for (int i = 0; i < l->num_eig_vect; i++) {
  //       printf0("vector #%02d: ", i + 1);
  //       apply_operator_double(l->vbuf_double[3], test_vectors[i], &(l->p_double), l);
  //       coarse_gamma5_double(l->vbuf_double[0], l->vbuf_double[3], 0, l->inner_vector_size, l);
  //       lambda = global_inner_product_double(test_vectors[i], l->vbuf_double[0], 0,
  //                                            l->inner_vector_size, l);
  //       lambda /= global_inner_product_double(test_vectors[i], test_vectors[i], 0,
  //                                             l->inner_vector_size, l);
  //       vector_double_saxpy(l->vbuf_double[1], l->vbuf_double[0], test_vectors[i], -lambda, 0,
  //                           l->inner_vector_size, l);
  //       mu = global_norm_double(l->vbuf_double[1], 0, l->inner_vector_size, l) /
  //            global_norm_double(test_vectors[i], 0, l->inner_vector_size, l);
  //       printf0("singular value: %+lf%+lfi, singular vector precision: %le\n", real(lambda),
  //               imag(lambda), mu);
  //     }
  //     printf0("--------------------------------------- depth: %d "
  //             "----------------------------------------\n",
  //             l->depth);
  //   }
  // //
  // #endif
}

void iterative_double_setup(int setup_iter, Level *l) {
  if (l->type == _FINE) {
    switch (l->global->interpolation) {
    case 1:
      inv_iter_inv_fcycle_double(setup_iter, l);
      break;
    case 2:
      read_tv_from_file_double(l);
      break;
    default:
      error0("invalid interpolation type\n");
      break;
    }
  }

  Level *lp = l;
  while (lp->currentLevel > 0) {
    testvector_analysis_double(lp->is_double.test_vector, lp);
    lp = lp->next_level;
    if (lp == nullptr)
      break;
  }
}
