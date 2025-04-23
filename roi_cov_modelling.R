# ===============================================================
# 0)  Packages ---------------------------------------------------
library(cmdstanr)      # cmdstanr for simulation
library(posterior)     # for handling draws
library(reshape2)      # melt()
library(dplyr)         # data wrangling
library(tidyr)         # pivot_*
library(ggplot2)       # plotting
library(lme4)          # lmer()
library(lmerTest)      # p-values for lmer
library(purrr)         # map_dfr()
library(scales)        # percent_format()
library(cowplot)       # plot_grid()

# ===============================================================
# 1)  Compile Stan generator ------------------------------------
mod <- cmdstan_model("sim_cortical_thickness.stan")

# ===============================================================
# 2)  Hyper-parameters (low-power scenario) ---------------------
stan_data <- list(
  n_subj         = 150,      # fewer subjects → low power
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
# 3)  Simulate one draw (fixed_param) ----------------------------
fit_sim <- mod$sample(
  data           = stan_data,
  chains         = 1,
  iter_sampling  = 1,
  iter_warmup    = 0,
  fixed_param    = TRUE,
  refresh        = 0
)

# pull into a posterior::draws_array
draws <- fit_sim$draws(format="draws_array")

# ===============================================================
# 4)  Extract & reshape Y  --------------------------------------
# 4a) drug assignment per subject
drug_vec <- as.integer(draws[1,1, grep("^drug\\[",   dimnames(draws)$variable)])

# 4b) thickness Y: one long vector, then array
Y_vec <- as.numeric(draws[1,1, grep("^Y\\[",     dimnames(draws)$variable)])
Y_array <- array(Y_vec,
                 dim = c(stan_data$n_subj,
                         stan_data$n_visit,
                         stan_data$n_roi))

# 4c) build a long data.frame
df <- reshape2::melt(Y_array,
                     varnames   = c("subj","visit","roi"),
                     value.name = "thk") %>%
  transmute(
    subj   = as.integer(subj),
    visit  = as.integer(visit),
    roi    = as.integer(roi),
    thk,
    drug   = drug_vec[subj]
  ) %>%
  mutate(
    visit_c  = visit - mean(visit),  # centered visit
    drug_num = drug                  # 0/1 numeric
  )

# ===============================================================
# 5)  Compute “true” ROI slopes ----------------------------------
beta_vec <- as.numeric(draws[1,1, grep("^beta_vr\\[", dimnames(draws)$variable)])
beta_mat <- array(beta_vec,
                  dim = c(stan_data$n_visit, stan_data$n_roi))
true_df <- tibble(
  ROI        = 1:stan_data$n_roi,
  true_slope = stan_data$gamma_global + beta_mat[stan_data$n_visit, ]
)

# ===============================================================
# 6)  Joint multilevel model (pooled) ----------------------------
fit_hier <- lmer(
  thk ~ visit_c * drug_num +
    (1 | subj) +
    (1 + visit_c:drug_num || roi),
  data = df,
  REML = FALSE
)
beta_fixed    <- fixef(fit_hier)["visit_c:drug_num"]
pooled_df <- ranef(fit_hier)$roi %>%
  as_tibble(rownames="ROI") %>%
  transmute(
    ROI          = as.integer(ROI),
    slope_pooled = `visit_c:drug_num` + beta_fixed
  )

# ===============================================================
# 7) Massive-univariate (unpooled) + FDR -------------------------
# You might get singular fit warnings. Ask yourself why this might
# be the case and how can you erase this warning by modifying 
# the generated data
unpooled <- map_dfr(
  1:stan_data$n_roi,
  function(r) {
    d <- filter(df, roi==r)
    m <- tryCatch(
      lmer(thk ~ visit_c * drug_num + (1|subj), data=d),
      error = function(e) NULL
    )
    if (is.null(m)) {
      return(tibble(ROI=r, slope_unp=NA_real_, se_unp=NA_real_, p_unp=NA_real_))
    }
    co <- summary(m)$coefficients
    if (!"visit_c:drug_num" %in% rownames(co))
      return(tibble(ROI=r, slope_unp=NA_real_, se_unp=NA_real_, p_unp=NA_real_))
    tibble(
      ROI       = r,
      slope_unp = co["visit_c:drug_num","Estimate"],
      se_unp    = co["visit_c:drug_num","Std. Error"],
      p_unp     = co["visit_c:drug_num","Pr(>|t|)"]
    )
  }
) %>%
  mutate(
    fdr     = p.adjust(p_unp, method="fdr"),
    sig_fdr = fdr < 0.05
  )

# ===============================================================
# 8)  Compare & summarise ----------------------------------------
compare_df <- pooled_df %>%
  left_join(unpooled, by="ROI") %>%
  left_join(true_df,   by="ROI")

# detection rates
det_raw <- mean(compare_df$p_unp < 0.05, na.rm=TRUE)
det_fdr <- mean(compare_df$sig_fdr,   na.rm=TRUE)
det_pooled <- 1  # by design we detect the global effect

## pooled: did we detect the global Drug×Time effect?  ----
coef_tbl <- summary(fit_hier)$coefficients
p_global <- coef_tbl["visit_c:drug_num", "Pr(>|t|)"]
sig_pooled <- (p_global < 0.05)

cat(sprintf("Pooled model global test p = %.3f → %s\n",
            p_global,
            if (sig_pooled) "SIGNIFICANT" else "not significant"))


## unpooled: % of ROIs detected (before & after FDR)  ----
det_raw  <- mean(compare_df$p_unp  < 0.05, na.rm=TRUE)
det_fdr  <- mean(compare_df$sig_fdr,        na.rm=TRUE)

cat(sprintf("Unpooled detection:\n • raw p<.05: %.1f%%\n • FDR<.05: %.1f%%\n",
            det_raw*100, det_fdr*100))

cat(sprintf(
  "Detection rates:\n • raw per‐ROI p<.05: %.1f%%\n • FDR<.05: %.1f%%\n • pooled: 100%%\n",
  det_raw*100, det_fdr*100
))

# RMSE
rmse <- function(x,y) sqrt(mean((x-y)^2, na.rm=TRUE))
rmse_unp <- rmse(compare_df$slope_unp, compare_df$true_slope)
rmse_pld <- rmse(compare_df$slope_pooled,compare_df$true_slope)
cat(sprintf("RMSE:\n • unpooled: %.3f\n • pooled:   %.3f\n",
            rmse_unp, rmse_pld))

# ===============================================================
# 9)  Plots -------------------------------------------------------

p1 <- ggplot(compare_df,
             aes(x=true_slope, y=slope_unp, color=sig_fdr)) +
  geom_abline(slope=1,intercept=0, linetype="dashed",col="grey60") +
  geom_point() +
  scale_color_manual(values=c("grey80","forestgreen")) +
  labs(title="Unpooled vs True", color="sig_FDR") +
  theme_minimal()

det_df <- tibble(
  Method = c("Pooled global", "Unpooled p<.05", "Unpooled FDR"),
  Detected = c(
    as.numeric(sig_pooled),
    det_raw,
    det_fdr
  )
)


p2 <- ggplot(det_df, aes(Method, Detected, fill=Method))+
  geom_col(width=0.6)+
  scale_y_continuous(labels = scales::percent_format(1), limits=c(0,1))+
  labs(title="% Effects Detected",
       y="% tests significant", x=NULL) +
  theme_minimal() +
  theme(legend.position="none")

p3 <- tibble(
  Method=c("Unpooled","Pooled"),
  RMSE=c(rmse_unp,rmse_pld)
) %>% ggplot(aes(Method,RMSE,fill=Method))+
  geom_col()+
  labs(title="RMSE of slope estimates")

plot_grid(p1,p2,p3,ncol=1,rel_heights=c(2,1,1))

# ===============================================================
# 8b)  Global effect recovery -----------------------------------
# — true vs estimated fixed-effect interaction -----------------

# 1) True parameter
true_gamma <- stan_data$gamma_global

# 2) Pooled estimate (from the joint multilevel model)
est_pooled <- beta_fixed

# 3) Unpooled “average” estimate (mean of all per-ROI slopes)
est_unp_mean <- mean(compare_df$slope_unp, na.rm=TRUE)

# 4) Bias and RMSE for the global effect estimates
bias_pooled   <- est_pooled   - true_gamma
bias_unpooled <- est_unp_mean - true_gamma
rmse_global <- function(est) sqrt((est - true_gamma)^2)
rmse_pooled   <- rmse_global(est_pooled)
rmse_unpooled <- rmse_global(est_unp_mean)

# 5) Print to console
cat(sprintf(
  "Global effect (%s) recovery:\n", "γ_global"
))
cat(sprintf(
  "  true:     %.3f\n  pooled:   %.3f  (bias %+ .3f, RMSE %.3f)\n",
  true_gamma, est_pooled,   bias_pooled,   rmse_pooled
))
cat(sprintf(
  "  unpooled: %.3f  (bias %+ .3f, RMSE %.3f)\n\n",
  est_unp_mean, bias_unpooled, rmse_unpooled
))

# 6) Visual comparison
eff_df <- tibble(
  Method   = c("Pooled", "Mean Unpooled", "True"),
  Estimate = c(est_pooled, est_unp_mean, true_gamma)
)

p_global <- ggplot(eff_df, aes(x = Method, y = Estimate, fill = Method)) +
  geom_col(data = filter(eff_df, Method != "True"), width = .6) +
  geom_point(data = filter(eff_df, Method == "True"),
             aes(x = Method, y = Estimate),
             color = "black", size = 4, shape = 18) +
  labs(
    title = "Recovery of Global Drug×Time Effect",
    y     = "Estimate of γ_global",
    x     = NULL
  ) +
  theme_minimal(base_size = 13) +
  scale_fill_manual(values = c("Pooled" = "tomato", "Mean Unpooled" = "steelblue")) +
  theme(legend.position = "none")


p_global

