# Green Asset Pricing: Code and Replication Files

This repository contains the code and data necessary to replicate the results presented in the paper **"Green Asset Pricing"** ([ECB Working Paper No. 2477](https://www.ecb.europa.eu/pub/pdf/scpwps/ecb.wp2477~e636f9c496.en.pdf)). Specifically, it provides scripts and models to replicate **Table II** and **Table III** from the paper, which compare the moments under different preferences and policy regimes.

---

## Repository Structure

### Core Files
- `/core_files/EZWH_ss.mod`: Dynare file embedding the macroeconomic model.
- `/core_files/DataFinance.xlsx`: Excel file containing the financial data used for estimation.

### Simulation and Estimation Files
- `/SMM/RUN_table_II.m`: MATLAB script that:
  - Calls Dynare.
  - Loads parameters under Epstein-Zin with habits (EZWH) and Constant Relative Risk Aversion (CRRA) preferences.
  - Computes the moments for **Table II**.

- `/SMM/RUN_table_III.m`: MATLAB script that:
  - Calls Dynare.
  - Loads parameters under EZWH and CRRA preferences for optimal vs. laissez-faire policies.
  - Computes the moments for **Table III**.

### Supporting Files
- `/SMM/smm_moment.m`: MATLAB function that generates the simulated moments.
- `/SMM/data_CARA_smm_temp.mat`: MATLAB file containing estimated parameters under CRRA preferences.
- `/SMM/data_ez_smm_temp.mat`: MATLAB file containing estimated parameters under EZWH preferences.
- `/SMM/call_EZ.mod`: Dynare file including Epstein-Zin Preferences with habits (EZWH).
- `/SMM/call_CRRA.mod`: Dynare file including Constant Risk Aversion (CRRA).
- `/SMM/empirical_moments.m`: MATLAB file that:
  - Reads data from `DataFinance.xlsx`.
  - Generates empirical moments for the analysis.

---

## How to Run the Code

### Prerequisites
1. **MATLAB** and **Dynare 6.x** installed.
2. Ensure that the required toolbox for reading Excel files in MATLAB is installed.
3. Clone this repository and set the base directory in MATLAB to the repository's root folder.

### Replicating Table II
1. Navigate to `/SMM/` in MATLAB.
2. Run the script:
   ```matlab
   RUN_table_II
   ```
3. The script will:
   - Call Dynare with the necessary models.
   - Load parameters for EZWH and CRRA preferences.
   - Output the moments for **Table II**.

### Replicating Table III
1. Navigate to `/SMM/` in MATLAB.
2. Run the script:
   ```matlab
   RUN_table_III
   ```
3. The script will:
   - Call Dynare with the necessary models.
   - Load parameters for EZWH and CRRA preferences, considering optimal vs. laissez-faire policies.
   - Output the moments for **Table III**.

---

## Notes
- Ensure `DataFinance.xlsx` is present in `/core_files/` before running the scripts.
- The empirical moments are generated via `empirical_moments.m`, which reads data from the Excel file and processes it for further analysis.
- The Dynare `.mod` files are pre-configured to run the macroeconomic models for both preference specifications (EZWH and CRRA).

---

## Citation
If you use this repository or any part of the code in your work, please cite:

> Benmir, G., Jaccard, I., & G. Vermandel "Green Asset Pricing." ECB Working Paper No. 2477.

---

For any issues or questions, please open an issue in this repository or contact the authors.

