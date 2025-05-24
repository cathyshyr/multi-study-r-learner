# Multi-Study $R$-Learner for Estimating Heterogeneous Treatment Effects Across Studies Using Statistical Machine Learning


This repository contains the code to reproduce the numerical results from the multi-study $R$-learner paper.

### Files:
1. **Utils.R** – Utility functions
2. **edat_orig.R** – Data preprocessing code for simulations using the `curatedOvarianData` R package
3. **ovarianDataSimulation_ScenarioA.R** – Ovarian cancer simulations (Scenario A: linear)
5. **ovarianDataSimulation_ScenarioB.R** – Ovarian cancer simulations (Scenario B: nonlinear)
6. **run_SimA_wrapper.R** - Wrapping function for ovarianDataSimulation_ScenarioA.R
7. **run_SimB_wrapper.R** - Wrapping function for ovarianDataSimulation_ScenarioB.R
8. **breastDataApplication_Revision.R** – Breast cancer data application with the `curatedBreastData` R package. The proc_curatedBreastDataExprSetList.RData can be downloaded [here](https://drive.google.com/file/d/1hpY4yuBgrZRRg8r_dn_P3jbRZkIORnvl/view?usp=sharing). Alternatively, please refer to the [documentation](https://www.bioconductor.org/packages/devel/data/experiment/manuals/curatedBreastData/man/curatedBreastData.pdf) on the `curatedBreastData` R package on how to obtain the curated breast data expression set. 

**Note:** The simulation code is designed for parallelization and is best run on a computing cluster. Alternatively, it can be executed in a standard loop, but the computational overhead may be significant.
