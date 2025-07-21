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
```
.
├── .gitignore
├── main.asv
├── main.m
├── README.md
├── sens_analysis.m
├── setup.m
│
├── data
│   ├── Luxembourg.xlsx
│   ├── wwtp1.xlsx
│   ├── wwtp2.xlsx
│   ├── wwtp3.xlsx
│   └── wwtp4.xlsx
│
├── img
│   ├── ...
│
├── parameters
│   ├── .gitkeep
│   ├── params_wwtp1.mat
│   ├── params_wwtp2.mat
│   ├── params_wwtp3.mat
│   └── params_wwtp4.mat
│
└── SEIRWWfiles
    ├── .gitkeep
    ├── paramFit.m
    ├── prediction.m
    ├── RWest.m
    ├── SEIRreaction.m
    ├── SEIRWWcalibrate.m
    ├── SEIRWWinit.m
    ├── SEIR_WW.asv
    ├── SEIR_WW.m
    ├── SEIR_WW_FWD.m
    ├── SEIR_WW_sens.m
    ├── SEIR_WW_sens_R.asv
    ├── SEIR_WW_sens_R.m
    └── WWinterpol.m
```
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
```

### 3. Run
Simply run the main.m file.
All necessary settings and configurations will be automatically adjusted based on the selected parameters.

## 👥 Authors

- Enea Dragoni - GitHub: [https://github.com/EneaDrago](https://github.com/EneaDrago) 
- Simone Cirelli - GitHub: [https://github.com/Cirocirotondo](https://github.com/Cirocirotondo) 
- Duccio Petreni - GitHub: [https://github.com/ducciopet](https://github.com/ducciopet) 
 
