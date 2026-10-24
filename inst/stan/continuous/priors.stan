log_s0 ~ normal(prior_log_s0[1], prior_log_s0[2]);
eta0 ~ normal(prior_eta0[1], prior_eta0[2]);
if (has_trend == 1)
  trend_slope[1] ~ normal(0, prior_trend_sd);
if (estimate_df == 1)
  eta_obs[1] ~ normal(prior_eta_obs[1], prior_eta_obs[2]);
to_vector(gamma) ~ normal(0, prior_variance_sd);
for (c in 1:C) {
  if (location_intercept[c] == 1)
    beta[c] ~ normal(prior_intercept[1], prior_intercept[2]);
  else
    beta[c] ~ normal(prior_coefficient[1], prior_coefficient[2]);
}
for (g in 1:G) {
  // Scale with df to reduce the variance-hyperparameter funnel; retain its Jacobian.
  target += continuous_standardized_variance_prior_lpdf(log_variance_raw[g] | nu0);
}
if (K1 > 0) {
  mu_omega_1 ~ normal(prior_mu_omega[1], prior_mu_omega[2]);
  tau_omega_1 ~ normal(0, prior_tau_omega_sd);
  for (g in 1:G) {
    z_omega_1[g] ~ std_normal();
    target += random_prior_1[g];
    L_1[g] ~ lkj_corr_cholesky(prior_lkj);
  }
}
if (K2 > 0) {
  mu_omega_2 ~ normal(prior_mu_omega[1], prior_mu_omega[2]);
  tau_omega_2 ~ normal(0, prior_tau_omega_sd);
  for (g in 1:G) {
    z_omega_2[g] ~ std_normal();
    target += random_prior_2[g];
    L_2[g] ~ lkj_corr_cholesky(prior_lkj);
  }
}
if (K3 > 0) {
  mu_omega_3 ~ normal(prior_mu_omega[1], prior_mu_omega[2]);
  tau_omega_3 ~ normal(0, prior_tau_omega_sd);
  for (g in 1:G) {
    z_omega_3[g] ~ std_normal();
    target += random_prior_3[g];
    L_3[g] ~ lkj_corr_cholesky(prior_lkj);
  }
}
if (K4 > 0) {
  mu_omega_4 ~ normal(prior_mu_omega[1], prior_mu_omega[2]);
  tau_omega_4 ~ normal(0, prior_tau_omega_sd);
  for (g in 1:G) {
    z_omega_4[g] ~ std_normal();
    target += random_prior_4[g];
    L_4[g] ~ lkj_corr_cholesky(prior_lkj);
  }
}
