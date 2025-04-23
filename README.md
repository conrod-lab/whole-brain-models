# ROI Covariance Modelling

A simulation framework for generating and analyzing longitudinal cortical thickness data with realistic ROI-by-ROI covariance structures. Compare per-ROI (massive univariate, but multilevel at the subject level) vs pooled multilevel mixed-effects models.

## Table of Contents
- [Features](#features)
- [Dependencies](#dependencies)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Scripts](#scripts)
- [Hyperparameter Tuning](#hyperparameter-tuning)
- [Project Structure](#project-structure)
- [Contributing](#contributing)

## Dependencies

This project requires:

- R (>= 4.0)
- CmdStanR (and CmdStan)
- R packages: cmdstanr, posterior, reshape2, dplyr, tidyr, ggplot2, lme4, lmerTest, purrr, scales, cowplot, shiny

Refer to [Quick Start](#quick-start) for installation commands.

## Features
This repository contains code to simulate longitudinal cortical‐thickness data under a drug intervention, with realistic ROI-by-ROI covariance structures, and to compare two analysis approaches:

1. **Massive univariate** (per‐ROI) mixed-effects models with FDR correction  
2. **Multilevel (pooled) mixed-effects** model using ROI-level partial pooling  

---

##  Quick Start
1. Install Dependencies
This project uses R (≥ 4.0) + CmdStanR. From an R console:
```
install.packages(c(
  "cmdstanr", "posterior", "reshape2", "dplyr",
  "tidyr",    "ggplot2",  "lme4",     "lmerTest",
  "purrr",    "scales",   "cowplot", "shiny",
))
cmdstanr::install_cmdstan()
```
Running roi_cov_modelling_v3.R will:

1. Simulate one fixed-param draw of Y[s,v,r] using your chosen Stan program.

2. Melt the 3D array into a long data frame.

3. Fit a pooled multilevel model (lmer() with (1 + visit:drug || roi)).

4. Fit unpooled per-ROI models and apply FDR correction.

5. Summarize detection rates & slope-recovery (RMSE, bias).

6. Produce diagnostic plots.


## Hyperparameter Tuning

Edit the `stan_data` list in each script to adjust:

n_subj, n_roi, n_visit in stan_data

Covariance parameters: rho_intra, rho_inter, rho_visit

Variance components: mu_i_sd, tau_i_sd, v_r_sd, subj_roi_sd, sigma_eps

Fixed effects: gamma_time, gamma_drug_int, gamma_global

These let you stress-test your analysis pipeline under more or less signal, varying ROI covariance structures, or sparser time-points.

## Usage

Run the scripts interactively in R or from the command line:

```bash
Rscript roi_cov_modelling_v3.R
Rscript roi_cov_plot.R
```

## Scripts

- `sim_cortical_thickness.stan`: Stan model for simulating longitudinal cortical thickness data.
- `roi_cov_modelling.R`: Simulate data, fit both per-ROI (unpooled) and pooled multilevel models, compare detection rates and slope recovery, and produce diagnostic plots.
- `roi_cov_plot.R`: Generate quick summary plots of global aging and baseline drug effects.

## Project Structure

```
.
├── sim_cortical_thickness.stan  # Stan simulation model
├── roi_cov_modelling.R       # Main analysis and comparison script
└── roi_cov_plot.R               # Quick visualization script
```

## Contributing

Contributions, issues, and feature requests are welcome. Please open an issue or submit a pull request.
