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
#include "preconditioner.h"

void preconditioner(vector_double phi, vector_double Dphi, vector_double eta, const int res,
                    Level *l) {

  Global *global = l->global;
  if (global->method == 0)
    vector_double_copy(phi, eta, 0, l->inner_vector_size, l);
  else if (global->method < 5 || global->method == 6 || !global->odd_even) {
    if (global->mixed_precision) {
      trans_double(l->sbuf_double[0], eta, l->s_double.op.translation_table, l);
      vcycle_double(l->sbuf_double[1], nullptr, l->sbuf_double[0], res, l);
      trans_back_double(phi, l->sbuf_double[1], l->s_double.op.translation_table, l);
    } else {
      trans_double(l->sbuf_double[0], eta, l->s_double.op.translation_table, l);
      vcycle_double(l->sbuf_double[1], nullptr, l->sbuf_double[0], res, l);
      trans_back_double(phi, l->sbuf_double[1], l->s_double.op.translation_table, l);
    }
  } else {
    if (global->mixed_precision) {
      //
      //       l->gmresSmoother.restartNumber = l->n_cy;
      //       l->gmresSmoother.initialGuessIsZero = res;
      //
      //       serial_to_oddeven_double(l->gmresSmoother.b, eta, l);
      //       if (global->method == 6) {
      //         g5D_solve_oddeven_double(&(l->gmresSmoother), &(l->oe_op_double), l);
      //       } else {
      //         solve_oddeven_double(&(l->gmresSmoother), &(l->oe_op_double), l);
      //       }
      //       oddeven_to_serial_double(phi, l->gmresSmoother.x, l);
    } else {

      l->gmresSmoother.restartNumber = l->n_cy;
      l->gmresSmoother.initialGuessIsZero = res;

      serial_to_oddeven_double(l->gmresSmoother.b, eta, l);
      if (global->method == 6) {
        g5D_solve_oddeven_double(&(l->gmresSmoother), &(l->oe_op_double), l);
      } else {
        solve_oddeven_double(&(l->gmresSmoother), &(l->oe_op_double), l);
      }
      oddeven_to_serial_double(phi, l->gmresSmoother.x, l);
    }
  }
  ASSERT(global->mixed_precision != 2);
  ASSERT(Dphi == nullptr);
}
