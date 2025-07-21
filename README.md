# SEIR-WW-EKF: Epidemic Modeling with Wastewater Data using an Extended Kalman Filter

This project implements an extended **SEIR** epidemiological model using **wastewater data (WW)** and **clinical case data (CC)**. An **Extended Kalman Filter (EKF)** is used to estimate the epidemic dynamics and the effective reproduction number $R_\mathrm{eff}$ in real time.

## 📌 Objectives

- Estimate $R_\mathrm{eff}$ in real time using:
  - clinical case data only
  - wastewater data only
  - both signals jointly
- Generate future case projections
- Evaluate the correlation between wastewater signals and clinical data

---

## 📁 Project Structure


---

## ▶️ How to Run the Project

### 1. Prerequisites

- MATLAB (recommended version R2020b or later)
- Required toolboxes:
  - Optimization Toolbox
  - Statistics and Machine Learning Toolbox
- [matlab2tikz](https://github.com/matlab2tikz/matlab2tikz) to export plots to LaTeX

### 2. Configuration

Edit the `main.m` script to customize the analysis:

```matlab
sensitivity_analysis = true;               % Enable/disable sensitivity analysis
sens_analysis_R_matrix = true;             % Sensitivity w.r.t. observation covariance matrix R
sens_analysis_dark_number = false;         % Sensitivity w.r.t. dark number 
regionName = 'wwtp1';                      % Region to analyze
new_dark_number = 1.5;                     % Value of the dark number to use
```

## 👥 Authors
Enea Dragoni - GitHub: https://github.com/EneaDrago
Simone Cirelli - GitHub: https://github.com/Cirocirotondo
Duccio Petreni - GitHub: https://github.com/ducciopet

