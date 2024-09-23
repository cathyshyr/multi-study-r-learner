# Multi-Study $R$-Learner for Estimating Heterogeneous Treatment Effects Across Studies Using Statistical Machine Learning


This repository contains the code to reproduce the numerical results from the multi-study $R$-learner paper.

### Files:
1. **Utils.R** – Utility functions
2. **edat_orig.R** – Data preprocessing code for simulations using the `curatedOvarianData` R package
3. **ovarianDataSimulation_ScenarioA.R** – Ovarian cancer simulations (Scenario A: linear)
4. **ovarianDataSimulation_ScenarioB.R** – Ovarian cancer simulations (Scenario B: nonlinear)
5. **breastDataApplication.R** – Breast cancer data application with the `curatedBreastData` R package

**Note:** The code is designed for parallelization and is best run on a computing cluster. Alternatively, it can be executed in a standard loop, but the computational overhead may be significant.
