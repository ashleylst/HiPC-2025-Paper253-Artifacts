#!/bin/bash -x
# budget account where contingent is taken from
#SBATCH --nodes=1

# can be omitted if --nodes and --ntasks-per-node
# are given
#SBATCH --ntasks-per-node=32
# if keyword omitted: Max. 96 tasks per node
# (SMT enabled, see comment below)
#SBATCH --cpus-per-task=4
# for OpenMP/hybrid jobs only

# if keyword omitted: Default is slurm-%j.out in
# the submission directory (%j is replaced by
# the job ID).

# if keyword omitted: Default is slurm-%j.out in
# the submission directory.
#SBATCH --time=10:00
#SBATCH --partition=cn-eth

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

export OMP_DISPLAY_ENV=true
export OMP_DISPLAY_AFFINITY=true
export OMP_AFFINITY_FORMAT="Thread Affinity: %0.3L %.8n %.15{thread_affinity} %.12H"
export OMP_PROC_BIND=true

cmake --build ./cmake-build-debug/ -j 14

srun ./cmake-build-debug/ddalphaamg --gtest_filter=DiracTest.MeasurementTest ./sample64x16to3.ini

