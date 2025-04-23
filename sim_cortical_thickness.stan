// ------------------------------------------------------------------
// Simulate longitudinal cortical‐thickness under:
//  • realistic block‐diagonal ROI covariance (functional networks)
//  • subject×ROI intercept jitter
//  • global ageing effect + drug main‐effect on intercept
//  • Drug×Time (global + ROI×Time deviations)
// ------------------------------------------------------------------
data {
  int<lower=1> n_subj;             // # subjects
  int<lower=1> n_roi;              // # cortical parcels
  int<lower=1> n_visit;            // # repeated visits per subject
  int<lower=1> n_net;              // # functional networks

  // which network each ROI belongs to
  array[n_roi] int<lower=1,upper=n_net> roi_net;

  // within‐ vs between‐network ROI correlations
  real<lower=0,upper=1> rho_intra;
  real<lower=0,upper=1> rho_inter;

  // AR(1) corr‐over‐time
  real<lower=-1,upper=1> rho_visit;

  // variance components
  real<lower=0> mu_i_sd;        // σ_μ: subj‐baseline SD
  real<lower=0> tau_i_sd;       // σ_τ: subj ageing‐slope SD
  real<lower=0> v_r_sd;         // σ_ν: ROI‐baseline SD
  real<lower=0> subj_roi_sd;    // σ_sr: subj×ROI intercept jitter
  real<lower=0> sigma_eps;      // σ_ε: measurement‐error SD

  // fixed‐effect parameters
  real gamma_time;              // global thinning per unit time
  real gamma_drug_int;          // drug main‐effect on intercept
  real gamma_global;            // Drug×Time interaction
}
transformed data {
  // 1) build ROI×ROI cov matrix
  matrix[n_roi,n_roi] K_roi;
  for (i in 1:n_roi) {
    for (j in 1:n_roi) {
      if (i == j) {
        K_roi[i,j] = 1;
      } else if (roi_net[i] == roi_net[j]) {
        K_roi[i,j] = rho_intra;
      } else {
        K_roi[i,j] = rho_inter;
      }
    }
  }

  // 2) build AR(1) over visits
  matrix[n_visit,n_visit] K_visit;
  for (i in 1:n_visit)
    for (j in 1:n_visit)
      K_visit[i,j] = pow(rho_visit, abs(i - j));
}
generated quantities {
  // 0) per‐subject drug assignment
  array[n_subj] int<lower=0,upper=1> drug;
  for (s in 1:n_subj)
    drug[s] = bernoulli_rng(0.5);

  // 1) subject‐level intercepts & ageing deviations
  vector[n_subj] mu_i;
  vector[n_subj] tau_i;
  for (s in 1:n_subj) {
    mu_i[s]  = normal_rng(2.5, mu_i_sd);
    tau_i[s] = normal_rng(0,   tau_i_sd);
  }

  // 2) ROI‐level intercepts via network‐block cov
  vector[n_roi] v_r;
  {
    // a) network covariance
    matrix[n_net,n_net] K_net;
    for (m in 1:n_net)
      for (n in 1:n_net)
        K_net[m,n] = (m==n ? 1 : rho_inter);
    matrix[n_net,n_net] L_net = cholesky_decompose(K_net);

    // b) draw one intercept per network
    vector[n_net] z_net;
    for (m in 1:n_net) z_net[m] = normal_rng(0,1);
    vector[n_net] w_net = v_r_sd * (L_net * z_net);

    // c) assign to each ROI
    for (r in 1:n_roi)
      v_r[r] = w_net[ roi_net[r] ];
  }

  // 3) ROI×visit Drug×Time deviations
  matrix[n_visit,n_roi] Z;
  for (v in 1:n_visit)
    for (r in 1:n_roi)
      Z[v,r] = normal_rng(0,1);
  matrix[n_visit,n_roi] beta_vr =
    cholesky_decompose(K_visit) * Z * cholesky_decompose(K_roi)';
  beta_vr *= 0.05;  // scale

  // 4) unstructured subj×ROI intercept jitter
  array[n_subj,n_roi] real u_sr;
  for (s in 1:n_subj)
    for (r in 1:n_roi)
      u_sr[s,r] = normal_rng(0, subj_roi_sd);

  // 5) simulate final thickness array
  array[n_subj,n_visit,n_roi] real Y;
  for (s in 1:n_subj) {
    for (v in 1:n_visit) {
      // *** floating‐point time scaling! ***
      real t = (v - 1) * 1.0 / (n_visit - 1);
      for (r in 1:n_roi) {
        real m =
            mu_i[s]                             // subj baseline
          + gamma_time   * t                   // global ageing
          + tau_i[s]     * t                   // subj‐specific slope
          + v_r[r]                             // ROI intercept
          + u_sr[s,r]                          // subj×ROI jitter
          + gamma_drug_int * drug[s]           // drug main‐effect
          + (gamma_global + beta_vr[v,r])      // Drug×Time
            * (drug[s] * t);

        Y[s,v,r] = normal_rng(m, sigma_eps);
      }
    }
  }
}
