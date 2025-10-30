
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
#include "unit_testing.h"

// Global g;
int my_rank;
#ifdef HAVE_HDF5
Hdf5_fileinfo h5info;
#endif

void EnvironmentDDalphaAMG::SetUp() {

  const char *copyright = "\n\n+----------------------------------------------------------+\n"
                          "| The DDalphaAMG solver library.                           |\n"
                          "| Copyright (C) 2016, Matthias Rottmann, Artur Strebel,    |\n"
                          "|       Simon Heybrock, Simone Bacchio, Bjoern Leder.      |\n"
                          "|                                                          |\n"
                          "| This program comes with ABSOLUTELY NO WARRANTY.          |\n"
                          "+----------------------------------------------------------+\n";

  char input_file_name[STRINGLENGTH];


#ifdef HAVE_HDF5
  h5info.filename = NULL;
  h5info.file_id = -1;
  h5info.rootgroup_id = -1;
  h5info.configgroup_id = -1;
  h5info.eigenmodegroup_id = -1;
  h5info.thiseigenmodegroup_id = -1;
  h5info.isOpen = 0;
  h5info.mode = -1;
#endif

#ifdef WRITE_LOGFILE
  EnvironmentDDalphaAMG::global.logfile = fopen("output.log", "w");
  fprintf(EnvironmentDDalphaAMG::global.logfile, "---------- log file -----------\n\n");
  fflush(EnvironmentDDalphaAMG::global.logfile);
#endif

  // initialize MPI
  MPI_Init(&test_argc, &test_argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &(EnvironmentDDalphaAMG::global.my_rank));
  my_rank = EnvironmentDDalphaAMG::global
                .my_rank; // NOTE: EnvironmentDDalphaAMG::global rank info needed for printf0(...)
                          // to not call MPI_Comm_rank each time, TODO: Check if it matters

  // So that only one process prints g-test
  ::testing::TestEventListeners& listeners =
  ::testing::UnitTest::GetInstance()->listeners();
  if (EnvironmentDDalphaAMG::global.my_rank != 0) {
    delete listeners.Release(listeners.default_result_printer());
  }

  // open input file
  //   std::string input_file_name((test_argc > 1) ? test_argv[1] : "sample.ini");
  strcpy(input_file_name, (test_argc > 1) ? test_argv[1] : "sample.ini");
  ASSERT((EnvironmentDDalphaAMG::global.inputfile = fopen(input_file_name, "r")) != nullptr);

  // read input file and store parameters into EnvironmentDDalphaAMG::global
  EnvironmentDDalphaAMG::global.readGlobalInfo();
  EnvironmentDDalphaAMG::global.allocateDataAfterReadingInputFile();
  EnvironmentDDalphaAMG::global.readNoDefaultInfo();
  EnvironmentDDalphaAMG::global.readSolverInfo();
  EnvironmentDDalphaAMG::global.readGeometryInfo();
  EnvironmentDDalphaAMG::global.readTestvectorIOData();
  EnvironmentDDalphaAMG::global.readEvaluationParameters();
  EnvironmentDDalphaAMG::global.readKcycleData();
  EnvironmentDDalphaAMG::global.validateParameters();
  EnvironmentDDalphaAMG::global.initialize();
  EnvironmentDDalphaAMG::global.setupFineDiracOperator();
  fclose(EnvironmentDDalphaAMG::global.inputfile);

  printf0("%s\nconfiguration: %s\n", copyright, EnvironmentDDalphaAMG::global.in);
  if (EnvironmentDDalphaAMG::global.rhs == 4)
    printf0("source list: %s\n", EnvironmentDDalphaAMG::global.source_list);

  if(global.method == 0){
    global.num_levels = 1;
  }
  // initialize multigrid level hierarchy
  EnvironmentDDalphaAMG::level = new Level[EnvironmentDDalphaAMG::global.num_levels];
  Level *level = EnvironmentDDalphaAMG::level;
  for (int i = 0; i < EnvironmentDDalphaAMG::global.num_levels; i++) {
    level[i].depth = i;
    level[i].type = (i == 0)                                              ? _FINE
                    : (i != EnvironmentDDalphaAMG::global.num_levels - 1) ? _INTERMEDIATE
                                                                          : _COARSEST;
    level[i].global = &EnvironmentDDalphaAMG::global;
    level[i].previous_level = (i == 0) ? nullptr : &level[i - 1];
    level[i].next_level =
        (i != EnvironmentDDalphaAMG::global.num_levels - 1) ? &level[i + 1] : nullptr;
    level[i].initialize();
  }

  method_update(level[0].setup_iter, &level[0]);

  // TODO : remove this solve, this shouldn't happen here
  //vector_double_define(level[0].gmres.b, 1, 0, level[0].inner_vector_size);
  //level[0].gmres.solve();

  // NOTE: Some code snippets I like to keep while developing/debugging!
  //     vector_double_define(level[0].gmres.b, 1, 0, level[0].inner_vector_size);
  //     level[0].gmres.applyOperator(level[0].gmres.x, level[0].gmres.b, level[0].gmres.matrix,
  //     &level[0] ); double testnorm = global_norm_double(level[0].gmres.x, 0,
  //     level[0].inner_vector_size, &level[0] ); printf0("matvec test: %lf\n", testnorm);
  //
  //     vector_double_define(level[0].s_double.b, 1, 0, level[0].inner_vector_size);
  //     level[0].s_double.applyOperator(level[0].s_double.x, level[0].s_double.b,
  //     level[0].s_double.matrix, &level[0] ); testnorm = global_norm_double(level[0].s_double.x,
  //     0, level[0].inner_vector_size, &level[0] ); printf0("schwarz matvec test: %lf\n",
  //     testnorm);
  //
  //     vector_double_define(level[1].s_double.b, 1, 0, level[1].inner_vector_size);
  //     apply_coarse_operator_double(level[1].s_double.x, level[1].s_double.b, &level[1].op_double,
  //     &level[1] ); testnorm = global_norm_double(level[1].s_double.x, 0,
  //     level[1].inner_vector_size, &level[1] ); printf0("matvec test coarse: %lf\n", testnorm);

  //     vector_double_define(level[1].s_double.b, 1, 0, level[1].inner_vector_size);
  //     apply_coarse_operator_double(level[1].s_double.x, level[1].s_double.b,
  //     &level[1].s_double.op, &level[1] ); testnorm = global_norm_double(level[1].s_double.x, 0,
  //     level[1].inner_vector_size, &level[1] ); printf0("schwarz matvec test coarse: %lf\n",
  //     testnorm);
  //
  //     vector_double_define(level[2].s_double.b, 1, 0, level[2].inner_vector_size);
  //     apply_coarse_operator_double(level[2].s_double.x, level[2].s_double.b,
  //     &level[2].s_double.op, &level[2] ); testnorm = global_norm_double(level[2].s_double.x, 0,
  //     level[2].inner_vector_size, &level[2] ); printf0("matvec test coarsest: %lf\n", testnorm);

  //   FILE* fin = NULL;
  //   ASSERT( (fin = fopen( "Dcpp", "w" )) != NULL );
  //   int size = 36*level[0].num_inner_lattice_sites;
  //   fprintf(fin, "size = %d\noe op = [\n", size);
  //   for( int kk=0; kk< size; kk++) {
  //     fprintf(fin, "%+lf%+lfi\n", real(level[0].oe_op_double.D[kk]),
  //     imag(level[0].oe_op_double.D[kk]));
  // //     fprintf(fin, "%lf\n", fabs(creal(l->next_level->p_PRECISION.b[kk])));
  // //     fprintf(fin, "%lf\n", fabs(creal(eta[kk])));
  //   }
  //   fprintf(fin, "];\n");
  //   fclose(fin);
  //   error0("");
}

void EnvironmentDDalphaAMG::TearDown() {

  delete[] EnvironmentDDalphaAMG::level;
  MPI_Finalize();
}
