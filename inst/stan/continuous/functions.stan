real continuous_student_sd_multiplier(real nu) {
  return sqrt((nu - 2) / nu);
}

real continuous_student_scale(real variance, real nu) {
  return sqrt(variance) * continuous_student_sd_multiplier(nu);
}

real continuous_variance_prior_lpdf(real log_variance, real nu0,
                                    real log_s0) {
  real shape = nu0 / 2;
  real log_rate = log(shape) + 2 * log_s0;
  return shape * log_rate - lgamma(shape) - shape * log_variance -
    exp(log_rate - log_variance);
}

real continuous_standardized_variance_prior_lpdf(real raw, real nu0) {
  real scale = sqrt(2 / nu0);
  return continuous_variance_prior_lpdf(raw * scale | nu0, 0) + log(scale);
}

matrix continuous_random_effect(matrix raw, vector omega, matrix cholesky) {
  return diag_pre_multiply(omega, cholesky) * raw;
}

// Gaussian conditional coordinates; the final row holds the log Jacobian in z-space.
matrix continuous_random_coordinates(matrix raw, vector omega, matrix correlation,
                                      vector residual, vector log_variance,
                                      matrix design, array[] int group, int use_data) {
  int K = rows(raw);
  int J = cols(raw);
  matrix[K + 1, J] result = rep_matrix(0, K + 1, J);
  if (use_data == 0) {
    result[1:K, ] = raw;
    return result;
  }
  matrix[K, K] scale = diag_pre_multiply(omega, correlation);
  array[J] matrix[K, K] precision;
  matrix[K, J] information = rep_matrix(0, K, J);
  for (j in 1:J) precision[j] = identity_matrix(K);
  for (n in 1:rows(design)) {
    vector[K] a = scale' * design[n]';
    real weight = exp(-log_variance[n]);
    precision[group[n]] += weight * (a * a');
    information[, group[n]] += weight * a * residual[n];
  }
  for (j in 1:J) {
    matrix[K, K] precision_cholesky = cholesky_decompose(precision[j]);
    vector[K] h = mdivide_left_tri_low(precision_cholesky, information[, j]);
    result[1:K, j] = mdivide_right_tri_low((h + raw[, j])', precision_cholesky)';
    result[K + 1, j] = -sum(log(diagonal(precision_cholesky)));
  }
  return result;
}

real continuous_partial_sum(array[] int features, int start, int end,
                                matrix y, matrix mu, matrix log_variance,
                                int student_t, real nu) {
  real lp = 0;
  for (s in 1:size(features)) {
    int g = features[s];
    vector[rows(y)] sd = exp(0.5 * log_variance[, g]);
    if (student_t == 1)
      lp += student_t_lpdf(y[, g] | nu, mu[, g],
                          sd * continuous_student_sd_multiplier(nu));
    else
      lp += normal_lpdf(y[, g] | mu[, g], sd);
  }
  return lp;
}
