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
7. [License](#license)

---

## Overview

This repository provides the MATLAB and Dynare code used to replicate the figures presented in the paper. The paper analyzes the effects of various carbon cap policies, frontloading mechanisms, and the market stability reserve on the European Union Emissions Trading System (EU-ETS) using a general equilibrium framework.

---

## Requirements

To run this replication package, you need:

- **MATLAB**: Version R2021a or later.
- **Dynare**: Version 5.x or 6.x.  
- **Operating System**: Compatible with Windows, macOS, or Linux.

### MATLAB Toolboxes
No specific toolboxes are required, but the Optimization Toolbox is recommended for advanced performance.

---

## Installation

1. Clone this repository to your local machine:
   ```bash
   git clone https://github.com/username/carbon-permit-banking.git
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

1. **Figure 6: "European-Union Emission Trading System Cap"**
   ```matlab
   run Figure_Cap.m
   ```

2. **Figure 7: "Projections under pre and post EU-ETS 2023 cap reform"**
   ```matlab
   dynare run_policy1_baseline.mod
   ```

3. **Figure 10: "The impacts of frontloading allowances"**
   ```matlab
   dynare run_policy2_frontloading.mod
   ```

4. **Figure 11: "Permit banking and the market stability reserve"**
   ```matlab
   dynare run_policy3_msr.mod
   ```

5. **Figure 12: "Cap policy versus carbon tax"**
   ```matlab
   dynare run_policy_comparison.mod
   ```

6. **Figure A.1: "Carbon emissions with limited borrowing (in million tons of CO2)"**
   ```matlab
   run_policy1_baseline_borro_graph.m
   ```

### Output
The scripts will generate:
- **Figures**: Saved as `.png` or `.pdf` in the `/results/figures` folder.
- **Simulation data**: Stored in `/results/simulations`.

---

## Project Structure

```
/carbon-permit-banking
├── README.md                # Documentation
├── LICENSE                  # License file
├── Figure_Cap.m             # Code for Figure 6
├── run_policy1_baseline.mod # Code for Figure 7
├── run_policy2_frontloading.mod # Code for Figure 10
├── run_policy3_msr.mod      # Code for Figure 11
├── run_policy_comparison.mod # Code for Figure 12
├── run_policy1_baseline_borro_graph.m # Code for Figure A.1
```

---

## Contact

For questions or issues regarding this replication package, please contact:

- **Loick Dubois**: [ldubois@london.edu](mailto:ldubois@london.edu)  
- **Jean-Guillaume Sahuc**: [jgsahuc@gmail.com](mailto:jgsahuc@gmail.com)  
- **Gauthier Vermandel**: [gauthier@vermandel.fr](mailto:gauthier@vermandel.fr)

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
