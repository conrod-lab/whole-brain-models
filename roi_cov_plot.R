# ===============================================================
# 0)  Packages ---------------------------------------------------
library(cmdstanr)
library(posterior)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(purrr)
library(scales)
library(cowplot)

# ===============================================================
# 1)  Compile Stan generator ------------------------------------
mod <- cmdstan_model("sim_cortical_thickness2.stan")   # your updated Stan

# ===============================================================
# 2)  Hyper-parameters (low-power scenario) ---------------------
stan_data <- list(
  n_subj         = 20,      # fewer subjects → low power
  n_roi          = 200,     # many ROIs → harsh multiple testing
  n_visit        = 3,       # three repeated visits
  
  # spatial/temporal smoothing
  n_net          = 5,       # assume 5 functional networks
  roi_net        = sample(1:5, 200, replace=TRUE),  # random assignment of each ROI to a network
  rho_intra      = 0.6,     # within-network ROI corr
  rho_inter      = 0.2,     # between-network ROI corr
  rho_visit      = 0.3,     # AR(1) temporal corr
  
  # variance components
  mu_i_sd        = 0.20,    # subj intercept SD (mm)
  tau_i_sd       = 0.05,    # subj slope SD
  v_r_sd         = 0.08,    # ROI intercept SD
  subj_roi_sd    = 0.05,    # subj×ROI jitter SD
  sigma_eps      = 0.40,    # measurement noise SD
  
  # fixed effects
  gamma_time     = -0.10,   # global thinning per unit time
  gamma_drug_int =  0.20,   # drug main‐effect on intercept
  gamma_global   = -0.05    # Drug×Time interaction
)


# ===============================================================
# 3)  Simulate one draw (fixed-param) ----------------------------
fit <- mod$sample(
  data           = stan_data,
  chains         = 1,
  iter_sampling  = 1,
  iter_warmup    = 0,
  fixed_param    = TRUE,
  refresh        = 0
)

# ===============================================================
# 3a)  Extract & reshape  ---------------------------------------
# pull out Y and drug assignment
Y_vec    <- as.numeric(fit$draws("Y",    format = "draws_matrix")[1, ])
drug_vec <- as.integer(fit$draws("drug", format = "draws_matrix")[1, ])

# build the 3d array and long data.frame
Y_array <- array(Y_vec,
                 dim = c(stan_data$n_subj,
                         stan_data$n_visit,
                         stan_data$n_roi))

df <- reshape2::melt(Y_array,
                     varnames   = c("subj", "visit", "roi"),
                     value.name = "thk") %>%
  mutate(
    subj  = as.integer(subj),
    visit = as.integer(visit),
    roi   = as.integer(roi),
    drug  = factor(drug_vec[subj],
                   levels = c(0,1),
                   labels = c("Placebo","Drug"))
  )

# ===============================================================
# 4)  Quick visualization of the TWO NEW GLOBAL PATTERNS --------

# A) Global ageing effect: average across ALL ROIs & subjects
df_age <- df %>%
  group_by(visit) %>%
  summarize(mean_th = mean(thk), .groups="drop")

p_age <- ggplot(df_age, aes(x = visit, y = mean_th)) +
  geom_line(size=1.2, color="steelblue") +
  geom_point(size=2, color="steelblue") +
  labs(
    title = "Overall average cortical thickness ↓ over time",
    subtitle = sprintf("True γ_time = %.2f", stan_data$gamma_time),
    x = "Visit (scaled)",
    y = "Mean thickness [mm]"
  ) +
  theme_minimal(base_size=14)

# B) Drug main‐effect on intercept: baseline difference
df_base <- df %>%
  filter(visit == 1) %>%
  group_by(drug) %>%
  summarize(
    mean_th = mean(thk),
    se_th   = sd(thk)/sqrt(n()),
    .groups = "drop"
  )

p_base <- ggplot(df_base, aes(x=drug, y=mean_th, fill=drug)) +
  geom_col(width=.6) +
  geom_errorbar(aes(ymin=mean_th - 1.96*se_th,
                    ymax=mean_th + 1.96*se_th),
                width=.2) +
  scale_fill_manual(values=c("grey70","firebrick")) +
  labs(
    title = "Baseline (Visit 1) thickness by Drug",
    subtitle = sprintf("True γ_drug_int = %.2f", stan_data$gamma_drug_int),
    x = "Drug group",
    y = "Mean thickness [mm]"
  ) +
  theme_minimal(base_size=14) +
  theme(legend.position="none")

# stitch them side‐by‐side
plot_grid(p_base, p_age, ncol=2, rel_widths=c(1,1.2))
