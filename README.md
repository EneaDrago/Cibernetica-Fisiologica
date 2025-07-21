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

C:.
│   .gitignore
│   A_data_correlation.m
│   A_plot_data.m
│   main.asv
│   main.m
│   README.md
│   sens_analysis.m
│   setup.m
│
├───data
│       Luxembourg.xlsx
│       wwtp1.xlsx
│       wwtp2.xlsx
│       wwtp3.xlsx
│       wwtp4.xlsx
│
├───img
│   ├───Luxembourg
│   │       img_casi_corr_casi_sens_Luxembourg.tex
│   │       img_casi_corr_ww_sens_Luxembourg.tex
│   │       img_cumsum_casi_corr_casi_sens_Luxembourg.tex
│   │       img_cumsum_casi_corr_ww_sens_Luxembourg.tex
│   │       img_epsilon_Luxembourg.tex
│   │       img_gamma_Luxembourg.tex
│   │       img_J_Luxembourg.tex
│   │       img_nu_Luxembourg.tex
│   │       img_ww_corr_casi_sens_Luxembourg.tex
│   │       img_ww_corr_ww_sens_Luxembourg.tex
│   │
│   ├───wwtp1
│   │       grafico_casi_cumulativi_wwtp1.tex
│   │       grafico_casi_giornalieri_wwtp1.tex
│   │       grafico_Reff_wwtp1.tex
│   │       img_casi_corr_casi_sens_wwtp1.tex
│   │       img_casi_corr_casi_wwtp1.tex
│   │       img_casi_corr_ww_sens_wwtp1.tex
│   │       img_casi_corr_ww_wwtp1.tex
│   │       img_cumsum_casi_corr_casi_sens_wwtp1.tex
│   │       img_cumsum_casi_corr_casi_wwtp1.tex
│   │       img_cumsum_casi_corr_ww_sens_wwtp1.tex
│   │       img_cumsum_casi_corr_ww_wwtp1.tex
│   │       img_data_wwtp1.tex
│   │       img_epsilon_wwtp1.tex
│   │       img_gamma_wwtp1.tex
│   │       img_J_wwtp1.tex
│   │       img_nu_wwtp1.tex
│   │       img_sensR_casi_wwtp1.tex
│   │       img_sensR_cumsum_wwtp1.tex
│   │       img_sensR_ww_wwtp1.tex
│   │       img_ww_corr_casi_sens_wwtp1.tex
│   │       img_ww_corr_casi_wwtp1.tex
│   │       img_ww_corr_ww_sens_wwtp1.tex
│   │       img_ww_corr_ww_wwtp1.tex
│   │
│   ├───wwtp2
│   │       grafico_casi_cumulativi_wwtp2.tex
│   │       grafico_casi_giornalieri_wwtp2.tex
│   │       grafico_Reff_wwtp2.tex
│   │       img_casi_corr_casi_wwtp2.tex
│   │       img_casi_corr_ww_wwtp2.tex
│   │       img_cumsum_casi_corr_casi_wwtp2.tex
│   │       img_cumsum_casi_corr_ww_wwtp2.tex
│   │       img_data_wwtp2.tex
│   │       img_epsilon_wwtp2.tex
│   │       img_gamma_wwtp2.tex
│   │       img_J_wwtp2.tex
│   │       img_nu_wwtp2.tex
│   │       img_sensR_casi_wwtp2.tex
│   │       img_sensR_cumsum_wwtp2.tex
│   │       img_sensR_ww_wwtp2.tex
│   │       img_ww_corr_casi_wwtp2.tex
│   │       img_ww_corr_ww_wwtp2.tex
│   │
│   ├───wwtp3
│   │       grafico_casi_cumulativi_wwtp3.tex
│   │       grafico_casi_giornalieri_wwtp3.tex
│   │       grafico_Reff_wwtp3.tex
│   │       img_casi_corr_casi_wwtp3.tex
│   │       img_casi_corr_ww_wwtp3.tex
│   │       img_cumsum_casi_corr_casi_wwtp3.tex
│   │       img_cumsum_casi_corr_ww_wwtp3.tex
│   │       img_data_wwtp3.tex
│   │       img_epsilon_wwtp3.tex
│   │       img_gamma_wwtp3.tex
│   │       img_J_wwtp3.tex
│   │       img_nu_wwtp3.tex
│   │       img_sensR_casi_wwtp3.tex
│   │       img_sensR_cumsum_wwtp3.tex
│   │       img_sensR_ww_wwtp3.tex
│   │       img_ww_corr_casi_wwtp3.tex
│   │       img_ww_corr_ww_wwtp3.tex
│   │
│   └───wwtp4
│           grafico_casi_cumulativi_wwtp4.tex
│           grafico_casi_giornalieri_wwtp4.tex
│           grafico_Reff_wwtp4.tex
│           img_casi_corr_casi_wwtp4.tex
│           img_casi_corr_ww_wwtp4.tex
│           img_cumsum_casi_corr_casi_wwtp4.tex
│           img_cumsum_casi_corr_ww_wwtp4.tex
│           img_data_wwtp4.tex
│           img_epsilon_wwtp4.tex
│           img_gamma_wwtp4.tex
│           img_J_wwtp4.tex
│           img_nu_wwtp4.tex
│           img_sensR_casi_wwtp4.tex
│           img_sensR_cumsum_wwtp4.tex
│           img_sensR_ww_wwtp4.tex
│           img_ww_corr_casi_wwtp4.tex
│           img_ww_corr_ww_wwtp4.tex
│
├───parameters
│       .gitkeep
│       params_wwtp1.mat
│       params_wwtp2.mat
│       params_wwtp3.mat
│       params_wwtp4.mat
│
└───SEIRWWfiles
        .gitkeep
        paramFit.m
        prediction.m
        RWest.m
        SEIRreaction.m
        SEIRWWcalibrate.m
        SEIRWWinit.m
        SEIR_WW.asv
        SEIR_WW.m
        SEIR_WW_FWD.m
        SEIR_WW_sens.m
        SEIR_WW_sens_R.asv
        SEIR_WW_sens_R.m
        WWinterpol.m
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

- Enea Dragoni - GitHub: [https://github.com/EneaDrago](https://github.com/EneaDrago) 
- Simone Cirelli - GitHub: [https://github.com/Cirocirotondo](https://github.com/Cirocirotondo) 
- Duccio Petreni - GitHub: [https://github.com/ducciopet](https://github.com/ducciopet) 
 
