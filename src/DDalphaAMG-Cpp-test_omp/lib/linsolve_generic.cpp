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

void print_Matrix_double(complex_double **A, int mv, int mh) {
  int i, j;

  // printf("\n\n");
  // for (i=0; i < mv; i++)
  // {
  //     for(j=0; j < mh; j++)
  //     {
  //             fprintf(stdout, "%6.6f +i%6.6f\t", real(A[i*mh + j]), imag(A[i*mh+j]));
  //     }
  //     fprintf(stdout, "\n");
  // }
  // printf("--\n");

  printf("\n\n");
  for (i = 0; i < mh; i++) {
    for (j = 0; j < mv; j++) {
      // fprintf(stdout, "%6.6f +i%6.6f\t", real(A[j*mh + i]), imag(A[j*mh+i]));
      fprintf(stdout, "%6.6f +i%6.6f\t", real(A[j][i]), imag(A[j][i]));
    }
    fprintf(stdout, "\n");
  }
  printf("--\n");
  printf("\n\n");
}

void print_double_eigenvalues(char *desc, int n, complex_double *w) {
  int j;
  printf("\n %s\n", desc);

  for (j = 0; j < n; j++) {
    printf(" (%6.2f,%6.2f)", real(w[j]), imag(w[j]));
  }
  printf("\n");
}

void fgmres_double_struct_init(gmres_double_struct *p) {

  /*********************************************************************************
   * Initializes all declared pointers with nullptr.
   *********************************************************************************/

  p->Z = nullptr;
  p->V = nullptr;
  p->H = nullptr;
  p->x = nullptr;
  p->b = nullptr;
  p->r = nullptr;
  p->w = nullptr;
  p->y = nullptr;
  p->gamma = nullptr;
  p->c = nullptr;
  p->s = nullptr;
  p->preconditioner = nullptr;
  p->eval_operator = nullptr;

  // copy of Hessenberg matrix
#if defined(GCRODR) && defined(POLYPREC)
  p->gcrodr_double.eigslvr.Hc = nullptr;
  p->polyprec_double.eigslvr.Hc = nullptr;
#elif defined(GCRODR)
  p->gcrodr_double.eigslvr.Hc = nullptr;
#elif defined(POLYPREC)
  p->polyprec_double.eigslvr.Hc = nullptr;
#endif

#ifdef POLYPREC
  p->polyprec_double.Hcc = nullptr;
  p->polyprec_double.L = nullptr;
  p->polyprec_double.col_prods = nullptr;
  p->polyprec_double.accum_prod = nullptr;
  p->polyprec_double.product = nullptr;
  p->polyprec_double.temp = nullptr;
  p->polyprec_double.h_ritz = nullptr;
  p->polyprec_double.lejas = nullptr;
  p->polyprec_double.random_rhs = nullptr;
  p->polyprec_double.xtmp = nullptr;

  p->polyprec_double.eigslvr.vl = nullptr;
  p->polyprec_double.eigslvr.vr = nullptr;
  p->polyprec_double.dirctslvr.ipiv = nullptr;
  p->polyprec_double.dirctslvr.x = nullptr;
  p->polyprec_double.dirctslvr.b = nullptr;
#endif

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  p->Va = nullptr;
  p->Za = nullptr;
#endif

#ifdef BLOCK_JACOBI
  p->block_jacobi_double.b_backup = nullptr;
  local_fgmres_double_struct_init(&(p->block_jacobi_double.local_p));
#endif
}

void fgmres_double_struct_alloc(
    int m, int n, long int vl, double tol, const int type, const int prec_kind,
    void (*precond)(vector_double, vector_double, vector_double, const int, Level *),
    void (*eval_op)(vector_double, vector_double, operator_struct<double> *, Level *),
    gmres_double_struct *p, Level *l) {

  /*********************************************************************************
   * Allocates memory for the fgmres struct and sets its values.
   * int m: Restart length
   * int n: Number of restarts
   * int vl: System size
   * double tol: Tolerance for relative residual
   * const int type: Specifies the problem for which fgmres should be applied
   *                 (_GLOBAL_FGMRES, _K_CYCLE, _COARSE_GMRES)
   * const int prec_kind: type of preconditioning: _RIGHT (flexible preconditioner),
   *                                               _LEFT (stationary preconditioner)
   *                                               or _NOTHING
   * void (*precond): Function pointer to the preconditioner
   *********************************************************************************/

  long int total = 0;
  int i, k = 0;

  p->restart_length = m;
  p->num_restart = n;

  p->preconditioner = precond;

  p->eval_operator = eval_op;
  p->tol = tol;
  p->kind = prec_kind;

#ifdef HAVE_TM1p1
  vl *= 2;
#endif

  if (m > 0) {
    total += (m + 1) * m; // Hessenberg matrix
    MALLOC(p->H, complex_double *, m);

    total += (5 + m) * vl; // x, r, b, w, V
    MALLOC(p->V, complex_double *, m + 1);

    if (precond != nullptr) {
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
      if (l->currentLevel == 0 && l->depth > 0) {
        total += (m + 2) * vl;
        k = m + 2;
        MALLOC(p->Z, complex_double *, k);
      } else {
#endif
        if (prec_kind == _RIGHT) {
          total += (m + 1) * vl; // Z
          k = m + 1;
        } else {
          total += vl;
          k = 1;
        }
        MALLOC(p->Z, complex_double *, k);
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
      }
#endif
    } else {
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
      if (l->currentLevel == 0 && l->depth > 0) {
        total += (m + 2) * vl;
        k = m + 2;
        MALLOC(p->Z, complex_double *, k);
      }
#else
      k = 0;
#endif
    }

#ifdef PERS_COMMS
    l->global->pers_comms_nrZs = k;
#endif

    total += 4 * (m + 1); // y, gamma, c, s

    p->H[0] = nullptr; // allocate connected memory
    MALLOC(p->H[0], complex_double, total);

    p->total_storage = total;
    total = 0;

    // ordering: H, y, gamma, c, s, w, V, Z, x, r, b
    // H
    for (i = 1; i < m; i++)
      p->H[i] = p->H[0] + i * (m + 1);
    total += m * (m + 1);

    // y
    p->y = p->H[0] + total;
    total += m + 1;
    // gamma
    p->gamma = p->H[0] + total;
    total += m + 1;
    // c
    p->c = p->H[0] + total;
    total += m + 1;
    // s
    p->s = p->H[0] + total;
    total += m + 1;
    // w
    p->w = p->H[0] + total;
    total += vl;
    // V
    for (i = 0; i < m + 1; i++) {
      p->V[i] = p->H[0] + total;
      total += vl;
    }
    // Z
    for (i = 0; i < k; i++) {
      p->Z[i] = p->H[0] + total;
      total += vl;
    }

    // x
    p->x = p->H[0] + total;
    total += vl;
    // r
    p->r = p->H[0] + total;
    total += vl;
    // b
    p->b = p->H[0] + total;
    total += vl;

    ASSERT(p->total_storage == total);
  }

  if (type == _GLOBAL_FGMRES) {
    p->timing = 1;
    p->print = l->global->vt.evaluation ? 0 : 1;
    p->initial_guess_zero = 1;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    p->op = &(l->global->op_double);
  } else if (type == _K_CYCLE) {
    // these settings also work for GMRES as a smoother
    p->timing = 0;
    p->print = 0;
    p->initial_guess_zero = 1;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    p->op = &(l->s_double.op);
  } else if (type == _COARSE_GMRES) {
    p->timing = 0;
    p->print = 0;
    p->initial_guess_zero = 1;
    p->layout = -1;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    if (l->global->odd_even)
      p->op = &(l->oe_op_double);
    else
      p->op = &(l->s_double.op);
  } else {
    ASSERT(type < 3);
  }

#if defined(GCRODR) || defined(POLYPREC)
  if (l->currentLevel == 0) {
#endif

    // FIXME : is this function-pointer-assignment really necessary ?
#if defined(GCRODR) || defined(POLYPREC)
    // p->polyprec_double.eigslvr.eigslvr_double = eigslvr_double;
#endif

    // copy of Hessenberg matrix
#if defined(GCRODR) && defined(POLYPREC)
    MALLOC(p->gcrodr_double.eigslvr.Hc, complex_double *, m);
    p->polyprec_double.eigslvr.Hc = p->gcrodr_double.eigslvr.Hc;
    p->gcrodr_double.eigslvr.Hc[0] = nullptr; // allocate connected memory
    MALLOC(p->gcrodr_double.eigslvr.Hc[0], complex_double, m * (m + 1));
    for (i = 1; i < m; i++)
      p->gcrodr_double.eigslvr.Hc[i] = p->gcrodr_double.eigslvr.Hc[0] + i * (m + 1);
    p->polyprec_double.eigslvr.Hc[0] = p->gcrodr_double.eigslvr.Hc[0];
#elif defined(GCRODR)
  MALLOC(p->gcrodr_double.eigslvr.Hc, complex_double *, m);
  p->gcrodr_double.eigslvr.Hc[0] = nullptr; // allocate connected memory
  MALLOC(p->gcrodr_double.eigslvr.Hc[0], complex_double, m * (m + 1));
  for (i = 1; i < m; i++)
    p->gcrodr_double.eigslvr.Hc[i] = p->gcrodr_double.eigslvr.Hc[0] + i * (m + 1);
#elif defined(POLYPREC)
  MALLOC(p->polyprec_double.eigslvr.Hc, complex_double *, m);
  p->polyprec_double.eigslvr.Hc[0] = nullptr; // allocate connected memory
  MALLOC(p->polyprec_double.eigslvr.Hc[0], complex_double, m * (m + 1));
  for (i = 1; i < m; i++)
    p->polyprec_double.eigslvr.Hc[i] = p->polyprec_double.eigslvr.Hc[0] + i * (m + 1);
#endif

#ifdef POLYPREC
    p->polyprec_double.d_poly = l->global->polyprec_d;
    int d_poly = p->polyprec_double.d_poly;

    MALLOC(p->polyprec_double.col_prods, complex_double, d_poly);
    MALLOC(p->polyprec_double.h_ritz, complex_double, d_poly);
    MALLOC(p->polyprec_double.lejas, complex_double, d_poly);
    MALLOC(p->polyprec_double.random_rhs, complex_double, vl);
    MALLOC(p->polyprec_double.accum_prod, complex_double, vl);
    MALLOC(p->polyprec_double.product, complex_double, vl);
    MALLOC(p->polyprec_double.temp, complex_double, vl);

    MALLOC(p->polyprec_double.xtmp, complex_double, vl);

    MALLOC(p->polyprec_double.Hcc, complex_double, d_poly * d_poly);
    MALLOC(p->polyprec_double.L, complex_double *, d_poly + 1);

    p->polyprec_double.L[0] = nullptr;

    MALLOC(p->polyprec_double.L[0], complex_double, (d_poly + 1) * d_poly);

    for (i = 1; i < d_poly + 1; i++) {
      p->polyprec_double.L[i] = p->polyprec_double.L[0] + i * d_poly;
    }

    MALLOC(p->polyprec_double.dirctslvr.ipiv, int, d_poly);
    MALLOC(p->polyprec_double.dirctslvr.x, complex_double, d_poly);
    MALLOC(p->polyprec_double.dirctslvr.b, complex_double, d_poly);

    p->polyprec_double.dirctslvr.N = d_poly;
    p->polyprec_double.dirctslvr.lda = d_poly; // m here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    p->polyprec_double.dirctslvr.ldb = d_poly;
    p->polyprec_double.dirctslvr.nrhs = 1;
    p->polyprec_double.dirctslvr.Hcc = p->polyprec_double.Hcc;
    p->polyprec_double.dirctslvr.dirctslvr_double = dirctslvr_double;

    MALLOC(p->polyprec_double.eigslvr.vl, complex_double, d_poly * d_poly);
    MALLOC(p->polyprec_double.eigslvr.vr, complex_double, d_poly * d_poly);

    p->polyprec_double.eigslvr.jobvl = 'N';
    p->polyprec_double.eigslvr.jobvr = 'N';

    p->polyprec_double.eigslvr.N = d_poly;
    p->polyprec_double.eigslvr.lda = p->restart_length + 1;
    p->polyprec_double.eigslvr.ldvl = d_poly;
    p->polyprec_double.eigslvr.ldvr = d_poly;
    p->polyprec_double.eigslvr.w = p->polyprec_double.h_ritz;
    p->polyprec_double.Hc = p->polyprec_double.eigslvr.Hc;
    p->polyprec_double.eigslvr.eigslvr_double = eigslvr_double;

    p->polyprec_double.update_lejas = 1;
    p->polyprec_double.preconditioner = nullptr;
    p->polyprec_double.preconditioner_bare = p->preconditioner;
    p->polyprec_double.syst_size = vl;

    p->polyprec_double.eigslvr.A = p->polyprec_double.Hc[0];
#endif

#if defined(GCRODR) || defined(POLYPREC)
  }
#endif

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  p->syst_size = vl;
  MALLOC(p->Va, complex_double *, m + 2);
  MALLOC(p->Za, complex_double *, m + 2);
  p->Va[0] = nullptr;
  p->Za[0] = nullptr;
  MALLOC(p->Va[0], complex_double, (m + 2) * vl);
  MALLOC(p->Za[0], complex_double, (m + 2) * vl);

  for (i = 1; i < m + 2; i++) {
    p->Va[i] = p->Va[0] + i * vl;
    p->Za[i] = p->Za[0] + i * vl;
  }

#ifdef PERS_COMMS
  l->global->pers_comms_nrZas = m + 2;
#endif
#endif

#ifdef BLOCK_JACOBI
  p->block_jacobi_double.syst_size = vl;

  if (l->currentLevel == 0) {
    // these two always go together
    p->block_jacobi_double.BJ_usable = 0;
    p->block_jacobi_double.local_p.polyprec_double.update_lejas = 1;

    MALLOC(p->block_jacobi_double.b_backup, complex_double, vl);
    MALLOC(p->block_jacobi_double.xtmp, complex_double, vl);

    p->block_jacobi_double.local_p.polyprec_double.d_poly = l->global->local_polyprec_d;

    local_fgmres_double_struct_alloc(
        l->global->local_polyprec_d, 1, vl, l->global->coarse_tol, _COARSE_GMRES, _NOTHING, nullptr,
        coarse_local_apply_schur_complement_double, &(p->block_jacobi_double.local_p), l);
  }
#endif
}

void fgcr_double(gmres_double_struct *p, Level *l) {

  /*********************************************************************************
   * Uses FGCR to solve the system D x = b, where b is taken from p->b and x is
   * stored in p->x.
   *********************************************************************************/

  int i, j = -1, finish = 0, iter = 0, il, ol;
  complex_double beta = 0, alpha;
  double r0_norm = 0, t0 = 0, t1 = 0;

  if (p->timing || p->print)
    t0 = MPI_Wtime();
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  if (p->print)
    printf0("+----------------------------------------------------------+\n");
#endif
  for (ol = 0; ol < p->num_restart && finish == 0; ol++) {

    if (ol == 0 && p->initial_guess_zero) {
      vector_double_copy(p->r, p->b, p->v_start, p->v_end, l);

    } else {
      apply_operator_double(p->w, p->x, p, l);                        // compute w = D*x
      vector_double_minus(p->r, p->b, p->w, p->v_start, p->v_end, l); // compute r = b - w
    }

    if (ol == 0) {
      r0_norm = global_norm_double(p->r, p->v_start, p->v_end, l);
    }

    for (il = 0; il < p->restart_length && finish == 0; il++) {

      j = il;
      iter++;

      p->preconditioner(p->V[j], nullptr, p->r, _NO_RES, l);
      apply_operator_double(p->Z[j], p->V[j], p, l);

      for (i = 0; i < j; i++) {
        beta = global_inner_product_double(p->Z[i], p->Z[j], p->v_start, p->v_end, l) / p->gamma[i];
        vector_double_saxpy(p->V[j], p->V[j], p->V[i], -beta, p->v_start, p->v_end, l);
        vector_double_saxpy(p->Z[j], p->Z[j], p->Z[i], -beta, p->v_start, p->v_end, l);
      }

      p->gamma[j] = global_inner_product_double(p->Z[j], p->Z[j], p->v_start, p->v_end, l);
      alpha = global_inner_product_double(p->Z[j], p->r, p->v_start, p->v_end, l) / p->gamma[j];
      vector_double_saxpy(p->x, p->x, p->V[j], alpha, p->v_start, p->v_end, l);
      vector_double_saxpy(p->r, p->r, p->Z[j], -alpha, p->v_start, p->v_end, l);

      alpha = global_norm_double(p->r, p->v_start, p->v_end, l) / r0_norm;
      if (real(alpha) < p->tol) {
        finish = 1;
        break;
      } else {
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
        if (iter % 10 == 0 && p->print)
          printf0("| approx. rel. res. after  %-6d iterations: %e |\n", iter, alpha);
#endif
      }
    } // end of restart
  }   // end of fgcr

  if (p->timing || p->print)
    t1 = MPI_Wtime();
  if (p->print) {
    apply_operator_double(p->w, p->x, p, l);
    vector_double_minus(p->r, p->b, p->w, p->v_start, p->v_end, l);
    beta = global_norm_double(p->r, p->v_start, p->v_end, l);
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    printf0("+----------------------------------------------------------+\n");
    printf0("\n");
#endif
    printf0("+----------------------------------------------------------+\n");
    printf0("|         FGCR iterations: %-6d                          |\n", iter);
    printf0("| exact relative residual: ||r||/||b|| = %e      |\n", real(beta) / r0_norm);
    printf0("| elapsed wall clock time: %-7lf seconds                |\n", t1 - t0);
    if (l->global->coarse_time > 0)
      printf0("|        coarse grid time: %-7lf seconds (%04.1lf%%)        |\n",
              l->global->coarse_time, 100 * (l->global->coarse_time / (t1 - t0)));
    printf0("+----------------------------------------------------------+\n\n");
  }
}
