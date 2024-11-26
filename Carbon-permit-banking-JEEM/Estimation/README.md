# Replication Package for "A General Equilibrium Approach to Carbon Permit Banking"

This repository contains the MATLAB and Dynare code required to replicate the results presented in the paper:

- **Title**: *A General Equilibrium Approach to Carbon Permit Banking*  
- **Authors**: Loick Dubois, Jean-Guillaume Sahuc, Gauthier Vermandel 
- **Published in**: *Journal of Environmental Economics and Management (JEEM)*, 2024  
- **DOI**: [https://doi.org/10.1016/j.jeem.2024.103076](https://www.sciencedirect.com/science/article/pii/S0095069624001505)

---

## Table of Contents
1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Project Structure](#project-structure)
6. [Contact](#contact)

---

## Overview

This repository provides the MATLAB and Dynare code used to replicate the figures presented in the paper. The paper analyzes the effects of various carbon cap policies, frontloading mechanisms, and the market stability reserve on the European Union Emissions Trading System (EU-ETS) using a general equilibrium framework.

---

## Requirements

To run this replication package, you need:

- **MATLAB**: Version R2021a or later.
- **Dynare**: Version 5.x only
- **Operating System**: Compatible with Windows, macOS, or Linux.

### MATLAB Toolboxes
No specific toolboxes are required, but the Optimization Toolbox is recommended for advanced performance.

---

## Installation

1. Clone this repository to your local machine:
   ```bash
   git clone https://github.com/Vermandel/ReplicationPackages/tree/main/Carbon-permit-banking-JEEM/Estimation
   cd carbon-permit-banking
   ```

2. Add the repository to your MATLAB path:
   ```matlab
   addpath(genpath('path/to/your/repository'));
   ```

3. Verify that Dynare is correctly installed:
   ```matlab
   dynare_version
   ```

---

## Usage

### Running the Scripts

Each `.mod` file corresponds to specific figures in the paper. Run the files using Dynare as follows:

1. **Table 2: "Estimated parameters"**, **Figure 3: "Impulse response functions"**
   ```matlab
   dynare DSV.mod
   ```


### Output
The scripts will generate:
- **Figures**: Saved as `.png` or `.pdf` in the `/results/figures` folder.
- **Simulation data**: Stored in `/results/simulations`.

---

## Project Structure

```
/carbon-permit-banking
├── README.md                			# Documentation
├── LICENSE                  			# License file
├── /estimation/             			# folder for the nonlinear estimation routine
├── DSV.mod 				 			# Dynare file for the model estimation
├── DSV_mle_estimates_temp.mat			# Matlab file containing estimated mode
├── DSV_mode.mat      					# Matlab file containing estimated mode
├── Data_Estimation_062009_122019.xlsx 	# Excel file containing the set of observable variables
```

---

## Contact

For questions or issues regarding this replication package, please contact:

- **Loick Dubois**: [ldubois@london.edu](mailto:ldubois@london.edu)  
- **Jean-Guillaume Sahuc**: [jgsahuc@gmail.com](mailto:jgsahuc@gmail.com)  
- **Gauthier Vermandel**: [gauthier@vermandel.fr](mailto:gauthier@vermandel.fr)

---
