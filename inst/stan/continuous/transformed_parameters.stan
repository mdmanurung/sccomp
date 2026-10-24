real<lower=2> nu0 = 2 + exp(eta0);
real<lower=0> variance_prior_scale = sqrt(2 / nu0);
vector[G] log_sigma_sq = rep_vector(2 * log_s0, G) +
  variance_prior_scale * log_variance_raw;
if (has_trend == 1)
  log_sigma_sq += 2 * trend_slope[1] * trend_covariate;
vector<lower=0>[G] sigma_sq = exp(log_sigma_sq);
real<lower=0> nu_obs = student_t == 0 ? 0 :
  (estimate_df == 1 ? 2 + exp(eta_obs[1]) : fixed_df);
matrix<lower=0>[G, K1] omega_1;
array[G] matrix[K1, J1] random_effect_1;
array[G * (K1 > 0)] real random_prior_1;
if (K1 > 0) {
  for (g in 1:G) {
    omega_1[g] = exp(mu_omega_1 + tau_omega_1 .* z_omega_1[g])';
    matrix[K1 + 1, J1] coordinates = continuous_random_coordinates(
      random_raw_1[g], omega_1[g]', L_1[g], y[, g] - X * beta[, g],
      log_sigma_sq[g] + W * gamma[, g], Z1, group1, use_data);
    random_effect_1[g] = continuous_random_effect(
      coordinates[1:K1, ], omega_1[g]', L_1[g]);
    random_prior_1[g] = std_normal_lpdf(to_vector(coordinates[1:K1, ])) +
      sum(coordinates[K1 + 1, ]);
  }
}
matrix<lower=0>[G, K2] omega_2;
array[G] matrix[K2, J2] random_effect_2;
array[G * (K2 > 0)] real random_prior_2;
if (K2 > 0) {
  for (g in 1:G) {
    omega_2[g] = exp(mu_omega_2 + tau_omega_2 .* z_omega_2[g])';
    matrix[K2 + 1, J2] coordinates = continuous_random_coordinates(
      random_raw_2[g], omega_2[g]', L_2[g], y[, g] - X * beta[, g],
      log_sigma_sq[g] + W * gamma[, g], Z2, group2, use_data);
    random_effect_2[g] = continuous_random_effect(
      coordinates[1:K2, ], omega_2[g]', L_2[g]);
    random_prior_2[g] = std_normal_lpdf(to_vector(coordinates[1:K2, ])) +
      sum(coordinates[K2 + 1, ]);
  }
}
matrix<lower=0>[G, K3] omega_3;
array[G] matrix[K3, J3] random_effect_3;
array[G * (K3 > 0)] real random_prior_3;
if (K3 > 0) {
  for (g in 1:G) {
    omega_3[g] = exp(mu_omega_3 + tau_omega_3 .* z_omega_3[g])';
    matrix[K3 + 1, J3] coordinates = continuous_random_coordinates(
      random_raw_3[g], omega_3[g]', L_3[g], y[, g] - X * beta[, g],
      log_sigma_sq[g] + W * gamma[, g], Z3, group3, use_data);
    random_effect_3[g] = continuous_random_effect(
      coordinates[1:K3, ], omega_3[g]', L_3[g]);
    random_prior_3[g] = std_normal_lpdf(to_vector(coordinates[1:K3, ])) +
      sum(coordinates[K3 + 1, ]);
  }
}
matrix<lower=0>[G, K4] omega_4;
array[G] matrix[K4, J4] random_effect_4;
array[G * (K4 > 0)] real random_prior_4;
if (K4 > 0) {
  for (g in 1:G) {
    omega_4[g] = exp(mu_omega_4 + tau_omega_4 .* z_omega_4[g])';
    matrix[K4 + 1, J4] coordinates = continuous_random_coordinates(
      random_raw_4[g], omega_4[g]', L_4[g], y[, g] - X * beta[, g],
      log_sigma_sq[g] + W * gamma[, g], Z4, group4, use_data);
    random_effect_4[g] = continuous_random_effect(
      coordinates[1:K4, ], omega_4[g]', L_4[g]);
    random_prior_4[g] = std_normal_lpdf(to_vector(coordinates[1:K4, ])) +
      sum(coordinates[K4 + 1, ]);
  }
}
