# HiPC-2025-Paper253-Artifacts
This repository contains the code, scripts, and plots associated with the paper titled:

"Performance-Portable Optimization and Analysis of Multiple Right-Hand Sides in a Lattice QCD Solver"

*Accepted for presentation at the 2025 International Conference on High Performance Computing (HiPC).*

## GitLab repository

The code relating to the paper is hosted here: https://git.uni-wuppertal.de/strebel/DDalphaAMG-Cpp

## Repository Structure

- `batchscripts/`: Contains batch scripts and input files to DDalphaAMG for running experiments on JUWELS, Ookami and HAICGU.
- `fig/`: Includes figure files used in the paper.
- `plot_script/`: Provides scripts for generating plots from the output data of experiments.
- `src/`
  - `DDalphaAMG-Cpp-papi`: source code of DDalphaAMG for PAPI readings
  - `DDalphaAMG-Cpp-test_omp`: source code of DDalphaAMG with the multiple rhs optimizations

## License


This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.


