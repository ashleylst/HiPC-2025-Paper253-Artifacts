#!/usr/bin/env bash

#SBATCH --job-name=ddalphaamg
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=128
#SBATCH --cpus-per-task=3
#SBATCH --time=20:00
#SBATCH -p all-nodes

# unload any modules currently loaded
module purge

omp_threads=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$omp_threads

# enable thread binding and print out info on thread affinity
export OMP_DISPLAY_ENV=true
export OMP_DISPLAY_AFFINITY=true
export OMP_AFFINITY_FORMAT="Thread Affinity: %0.3L %.8n %.15{thread_affinity} %.12H"
export OMP_PROC_BIND=true

# load the gcc 11 module
#module load CPE
#module load cray-mvapich2_nogpu_sve/2.3.6
module load openmpi/gcc13.1.0/4.1.5 gcc/13.1.0 openblas cmake/3.25.2
#module load arm-modules/23.10
#module load openmpi/arm23.10/5.0.0
module load slurm

cmake --build /lustre/home/shlong/test-omp/DDalphaAMG-Cpp/cmake-build-debug --target ddalphaamg -j 32

export OMP_PLACES=sockets

mpirun -n 2048 /lustre/home/shlong/test-omp/DDalphaAMG-Cpp/cmake-build-debug/ddalphaamg --gtest_filter=DiracTest.MeasurementTest /lustre/home/shlong/test-omp/DDalphaAMG-Cpp/sample128.ini