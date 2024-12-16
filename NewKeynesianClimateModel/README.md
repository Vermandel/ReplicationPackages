# New Keynesian Climate Model Replication

This repository replicates the **New Keynesian Climate Model** as presented in the CEPR discussion paper [DP19745](https://cepr.org/publications/dp19745). The simulations, figures, and results are implemented using **MATLAB** with the **Dynare package** (version 6.x).

## Authors
- **Jean-Guillaume Sahuc**
- **Frank Smets**
- **Gauthier Vermandel**

## Prerequisites

1. **MATLAB** (version compatible with Dynare 6.x).
2. **Dynare 6.x**:
   - Download Dynare 6.x from the official website: [https://www.dynare.org/download/](https://www.dynare.org/download/).
   - Follow the installation instructions for MATLAB.

## File Descriptions

### Core Files
- `model_file.mod` :
  - This file contains the **core implementation** of the New Keynesian Climate Model.

- `run_simulations.mod` :
  - Runs simulations to generate:
    - Out-of-sample predictions up to 2100.
    - Six trends discussed in the paper.
    - Decomposition of aggregate demand and supply under different scenarios.
    - Changes in monetary policy rules.

- `run_IRFs.mod` :
  - Simulates **Impulse Response Functions (IRFs)** for the model.

- `get_Z.m` :
  - MATLAB script to compute the **asymptotic TFP value**.

- `SSV_sims0.m` :
  - Provides a **guess value** for the path of endogenous variables. This step is not compulsory as simulations can work without it, but it may improve convergence.

### Estimation Data
- `run_estimation_filtered_data.mat` :
  - Contains the filtered shocks and variables from the estimation.

- `run_estimation_mle_temp.mat` or `run_estimation_mode.mat` :
  - Files that store the **estimated parameters** of the model.

## Steps to Run the Model

1. Install Dynare 6.x and ensure it is correctly added to your MATLAB path.
2. Open MATLAB and navigate to the repository folder.
3. Run the following files in sequence depending on your needs:
   - Simulations:
     ```
     dynare run_simulations.mod
     ```
   - Impulse Response Functions:
     ```
     dynare run_IRFs.mod
     ```

## Notes
- Ensure all `.mod` and `.m` files are in the same directory.
- Dynare must be initialized before running `.mod` files.
- Results include trends, out-of-sample predictions, decompositions, and monetary policy variations as detailed in the paper.

## Contact
For questions or issues, please refer to the CEPR paper or contact the authors:
- **Gauthier Vermandel**: [gauthier@vermandel.fr](mailto:gauthier@vermandel.fr)
