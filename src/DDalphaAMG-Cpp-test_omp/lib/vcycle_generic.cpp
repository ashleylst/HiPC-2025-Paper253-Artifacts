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
#include "vcycle_generic.h"

void smoother_double(vector_double phi, vector_double Dphi, vector_double eta, int n, int res,
                     Level *l) {

  Global *global = l->global;
  ASSERT(phi != eta);

  if (global->method > 0 && global->method < 4) {
    l->s_double.has_residual = res;
    l->s_double.number_of_cycles = n;
    vector_double_copy(l->s_double.x, phi, 0, l->inner_vector_size, l);
    vector_double_copy(l->s_double.b, eta, 0, l->inner_vector_size, l);
    l->s_double.solve();
    vector_double_copy(phi, l->s_double.x, 0, l->inner_vector_size, l); // TODO: pretty ugly
  } else {
    l->gmresSmoother.initialGuessIsZero = res;
    l->gmresSmoother.restartNumber = n;
    if (global->method == 4 || global->method == 6) {
      if (global->odd_even) {
        if (res == _RES) {
          d_plus_clover_double(l->gmresSmoother.x, phi, &(l->s_double.op), l);
          vector_double_minus(l->gmresSmoother.x, eta, l->gmresSmoother.x, 0, l->inner_vector_size,
                              l);
        }
        block_to_oddeven_double(l->gmresSmoother.b, res == _RES ? l->gmresSmoother.x : eta, l);
        l->gmresSmoother.initialGuessIsZero = _NO_RES;
        if (global->method == 6) {
          if (l->type == _FINE)
            g5D_solve_oddeven_double(&(l->gmresSmoother), &(l->oe_op_double), l);
          else
            g5D_coarse_solve_odd_even_double(&(l->gmresSmoother), &(l->oe_op_double), l);
        } else {
          if (l->type == _FINE)
            solve_oddeven_double(&(l->gmresSmoother), &(l->oe_op_double), l);
          else
            coarse_solve_odd_even_double(&(l->gmresSmoother), &(l->oe_op_double), l);
        }
        if (res == _NO_RES) {
          oddeven_to_block_double(phi, l->gmresSmoother.x, l);
        } else {
          oddeven_to_block_double(l->gmresSmoother.b, l->gmresSmoother.x, l);
          vector_double_plus(phi, phi, l->gmresSmoother.b, 0, l->inner_vector_size, l);
        }
      } else {
        vector_double_copy(l->gmresSmoother.x, phi, 0, l->inner_vector_size, l);
        vector_double_copy(l->gmresSmoother.b, eta, 0, l->inner_vector_size, l);
        l->gmresSmoother.solve();
        vector_double_copy(phi, l->gmresSmoother.x, 0, l->inner_vector_size, l);
      }
    }
    ASSERT(Dphi == nullptr);
  }
}

void vcycle_double(vector_double phi, vector_double Dphi, vector_double eta, int res, Level *l) {

  Global *global = l->global;
  if (global->interpolation && l->type != _COARSEST) {
    for (int i = 0; i < l->n_cy; i++) {
      if (i == 0 && res == _NO_RES) {
        restrict_double(l->next_level->vbuf_double[0], eta, l);
      } else {
        l->gmresSmoother.applyOperator(l->vbuf_double[2], phi, &l->op_double,
                                       l); // NOTE: Not sure about
        vector_double_minus(l->vbuf_double[3], eta, l->vbuf_double[2], 0, l->inner_vector_size, l);
        restrict_double(l->next_level->vbuf_double[0], l->vbuf_double[3], l);
      }
      if (!l->next_level->idle) {

        if (l->depth == 0)
          global->coarse_time -= MPI_Wtime();

        if (l->next_level->type != _COARSEST) {
          if (global->kcycle) {
            vector_double_copy(l->next_level->gmresKCycle.b, l->next_level->vbuf_double[0], 0,
                               l->next_level->inner_vector_size, l->next_level);
            l->next_level->gmresKCycle.solve();
            vector_double_copy(l->next_level->vbuf_double[1], l->next_level->gmresKCycle.x, 0,
                               l->next_level->inner_vector_size, l->next_level);
          } else
            vcycle_double(l->next_level->vbuf_double[1], nullptr, l->next_level->vbuf_double[0],
                          _NO_RES, l->next_level);
        } else {
          vector_double_copy(l->next_level->gmresCoarseSolve.b, l->next_level->vbuf_double[0], 0,
                             l->next_level->inner_vector_size, l->next_level);
          if (l->odd_even) {
            if (global->method == 6)
              g5D_coarse_solve_odd_even_double(&(l->next_level->gmresCoarseSolve),
                                               &(l->next_level->oe_op_double), l->next_level);
            else
              coarse_solve_odd_even_double(&(l->next_level->gmresCoarseSolve),
                                           &(l->next_level->oe_op_double), l->next_level);
          } else {
            l->next_level->gmresCoarseSolve.solve();
          }
          vector_double_copy(l->next_level->vbuf_double[1], l->next_level->gmresCoarseSolve.x, 0,
                             l->next_level->inner_vector_size, l->next_level);
        }
        if (l->depth == 0)
          global->coarse_time += MPI_Wtime();
      }
      if (i == 0 && res == _NO_RES)
        interpolate3_double(phi, l->next_level->vbuf_double[1], l);
      else
        interpolate_double(phi, l->next_level->vbuf_double[1], l);
      smoother_double(phi, Dphi, eta, l->post_smooth_iter, _RES, l);
      res = _RES;
    }
  } else {
    smoother_double(phi, Dphi, eta, (l->depth == 0) ? l->n_cy : l->post_smooth_iter, res, l);
  }
}
