functions {
  #include common_functions.stan

  real student_t_scale_from_log_precision(real student_t_df,
                                           real log_precision) {
    return sqrt((student_t_df - 2) / student_t_df) *
      exp(-0.5 * log_precision);
  }

  int pair_is_included(int n, int m, array[,] int excluded_pairs) {
    int n_excluded = dims(excluded_pairs)[1];
    for (k in 1:n_excluded) {
      if (excluded_pairs[k, 1] == n && excluded_pairs[k, 2] == m)
        return 0;
    }
    return 1;
  }

  real abundance_variability_regression_student_t(
      row_vector log_precision,
      row_vector abundance,
      real prec_intercept_1,
      real prec_slope_1,
      real prec_slope_2,
      real prec_intercept_2,
      real prec_sd,
      int bimodal_mean_variability_association,
      real mix_p,
      int exclude_mean_variability_association) {
    real slope_1 = exclude_mean_variability_association == 1 ? 0 : prec_slope_1;
    real slope_2 = exclude_mean_variability_association == 1 ? 0 : prec_slope_2;
    real lp = 0;

    if (bimodal_mean_variability_association == 1) {
      for (m in 1:cols(log_precision)) {
        lp += log_mix(
          mix_p,
          normal_lpdf(log_precision[m] |
                      abundance[m] * slope_1 + prec_intercept_1, prec_sd),
          normal_lpdf(log_precision[m] |
                      abundance[m] * slope_2 + prec_intercept_2, prec_sd)
        );
      }
    } else {
      lp = normal_lpdf(
        log_precision |
        abundance * slope_1 + prec_intercept_1,
        prec_sd
      );
    }
    return lp;
  }

  real partial_sum_student_t_lpmf(
      array[] int idx_y,
      int start,
      int end,
      matrix y,
      matrix X,
      matrix beta,
      matrix Xa,
      matrix alpha,
      real student_t_df,
      int M,
      array[] int ncol_X_random_eff,
      matrix X_random_effect_1,
      matrix X_random_effect_2,
      matrix X_random_effect_3,
      matrix X_random_effect_4,
      matrix random_effect_1,
      matrix random_effect_2,
      matrix random_effect_3,
      matrix random_effect_4,
      array[,] int excluded_pairs) {
    int N_slice = size(idx_y);
    matrix[N_slice, M] mu = X[idx_y, ] * beta;
    matrix[N_slice, M] log_precision = Xa[idx_y, ] * alpha;
    real lp = 0;

    if (ncol_X_random_eff[1] > 0)
      mu += X_random_effect_1[idx_y, ] * random_effect_1;
    if (ncol_X_random_eff[2] > 0)
      mu += X_random_effect_2[idx_y, ] * random_effect_2;
    if (ncol_X_random_eff[3] > 0)
      mu += X_random_effect_3[idx_y, ] * random_effect_3;
    if (ncol_X_random_eff[4] > 0)
      mu += X_random_effect_4[idx_y, ] * random_effect_4;

    for (i in 1:N_slice) {
      for (m in 1:M) {
        int n = idx_y[i];
        if (pair_is_included(n, m, excluded_pairs)) {
          lp += student_t_lpdf(
            y[n, m] |
            student_t_df,
            mu[i, m],
            student_t_scale_from_log_precision(
              student_t_df, log_precision[i, m]
            )
          );
        }
      }
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> C;
  int<lower=1> A;
  int<lower=0> A_intercept_columns;
  int<lower=0> B_intercept_columns;
  matrix[N, M] y;
  matrix[N, C] X;
  matrix[N, A] Xa;
  array[A] int<lower=1, upper=C> variability_to_composition_map;

  int<lower=0> TNS;
  int<lower=0> TNIM;
  array[TNIM, 2] int<lower=1> truncation_not_idx_minimal;

  array[2] real prior_prec_intercept;
  array[2] real prior_prec_slope;
  array[2] real prior_prec_sd;
  array[2] real prior_mean_intercept;
  array[2] real prior_mean_coefficients;

  int<lower=0, upper=1> exclude_mean_variability_association;
  int<lower=0, upper=1> bimodal_mean_variability_association;
  int<lower=0, upper=1> use_data;
  int<lower=1> grainsize;
  int<lower=0, upper=1> intercept_in_design;
  int<lower=0, upper=1> enable_loo;
  real<lower=2> student_t_df;

  int<lower=0, upper=1> is_random_effect;
  array[4] int<lower=0> ncol_X_random_eff;
  matrix[N, ncol_X_random_eff[1]] X_random_effect_1;
  matrix[N, ncol_X_random_eff[2]] X_random_effect_2;
  matrix[N, ncol_X_random_eff[3]] X_random_effect_3;
  matrix[N, ncol_X_random_eff[4]] X_random_effect_4;

  array[4] int<lower=0> n_groups;
  array[4] int<lower=0> how_many_factors_in_random_design;
  array[how_many_factors_in_random_design[1], n_groups[1]] int group_factor_indexes_for_covariance_1;
  array[how_many_factors_in_random_design[2], n_groups[2]] int group_factor_indexes_for_covariance_2;
  array[how_many_factors_in_random_design[3], n_groups[3]] int group_factor_indexes_for_covariance_3;
  array[how_many_factors_in_random_design[4], n_groups[4]] int group_factor_indexes_for_covariance_4;
}

transformed data {
  array[N] int array_N;
  int ncol_X_random_eff_safe_1 = max(ncol_X_random_eff[1], 1);
  int ncol_X_random_eff_safe_2 = max(ncol_X_random_eff[2], 1);
  int ncol_X_random_eff_safe_3 = max(ncol_X_random_eff[3], 1);
  int ncol_X_random_eff_safe_4 = max(ncol_X_random_eff[4], 1);
  for (n in 1:N) array_N[n] = n;
}

parameters {
  matrix[C, M] beta;
  matrix[A, M] alpha;
  array[A] ordered[1 + bimodal_mean_variability_association] prec_intercept;
  array[A] real prec_slope_1;
  array[A * bimodal_mean_variability_association] real prec_slope_2;
  array[A] real log_prec_sd;
  real<lower=0, upper=1> mix_p;

  array[ncol_X_random_eff[1] * (ncol_X_random_eff[1] > 0)] vector[M] random_effect_raw_1;
  array[M * (ncol_X_random_eff[1] > 0)] vector[how_many_factors_in_random_design[1]] random_effect_sigma_raw_1;
  array[M * (ncol_X_random_eff[1] > 0)] cholesky_factor_corr[how_many_factors_in_random_design[1] * (ncol_X_random_eff[1] > 0)] sigma_correlation_factor_1;

  array[ncol_X_random_eff[2] * (ncol_X_random_eff[2] > 0)] vector[M] random_effect_raw_2;
  array[M * (ncol_X_random_eff[2] > 0)] vector[how_many_factors_in_random_design[2]] random_effect_sigma_raw_2;
  array[M * (ncol_X_random_eff[2] > 0)] cholesky_factor_corr[how_many_factors_in_random_design[2] * (ncol_X_random_eff[2] > 0)] sigma_correlation_factor_2;

  array[ncol_X_random_eff[3] * (ncol_X_random_eff[3] > 0)] vector[M] random_effect_raw_3;
  array[M * (ncol_X_random_eff[3] > 0)] vector[how_many_factors_in_random_design[3]] random_effect_sigma_raw_3;
  array[M * (ncol_X_random_eff[3] > 0)] cholesky_factor_corr[how_many_factors_in_random_design[3] * (ncol_X_random_eff[3] > 0)] sigma_correlation_factor_3;

  array[ncol_X_random_eff[4] * (ncol_X_random_eff[4] > 0)] vector[M] random_effect_raw_4;
  array[M * (ncol_X_random_eff[4] > 0)] vector[how_many_factors_in_random_design[4]] random_effect_sigma_raw_4;
  array[M * (ncol_X_random_eff[4] > 0)] cholesky_factor_corr[how_many_factors_in_random_design[4] * (ncol_X_random_eff[4] > 0)] sigma_correlation_factor_4;

  array[4 * (is_random_effect > 0)] real random_effect_sigma_mu;
  array[4 * (is_random_effect > 0)] real random_effect_sigma_sigma;
  array[is_random_effect > 0] real zero_random_effect;
}

transformed parameters {
  array[A] real prec_intercept_1;
  array[A * bimodal_mean_variability_association] real prec_intercept_2;
  array[A] real<lower=0> prec_sd;
  real mix_p_scalar = bimodal_mean_variability_association == 1 ? mix_p : 0.5;

  matrix[ncol_X_random_eff_safe_1 * (is_random_effect > 0), M] random_effect_1;
  matrix[ncol_X_random_eff_safe_2 * (is_random_effect > 0), M] random_effect_2;
  matrix[ncol_X_random_eff_safe_3 * (is_random_effect > 0), M] random_effect_3;
  matrix[ncol_X_random_eff_safe_4 * (is_random_effect > 0), M] random_effect_4;

  for (a in 1:A) {
    prec_intercept_1[a] = prec_intercept[a][1];
    prec_sd[a] = exp(log_prec_sd[a]);
    if (bimodal_mean_variability_association == 1)
      prec_intercept_2[a] = prec_intercept[a][2];
  }

  if (ncol_X_random_eff[1] > 0) {
    array[ncol_X_random_eff[1]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[1]) raw_vec[i] = random_effect_raw_1[i];
    random_effect_1 = build_re_block(
      M, n_groups[1], how_many_factors_in_random_design[1], ncol_X_random_eff[1],
      group_factor_indexes_for_covariance_1, raw_vec,
      random_effect_sigma_mu[1], random_effect_sigma_sigma[1],
      random_effect_sigma_raw_1, sigma_correlation_factor_1
    );
  }
  if (ncol_X_random_eff[2] > 0) {
    array[ncol_X_random_eff[2]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[2]) raw_vec[i] = random_effect_raw_2[i];
    random_effect_2 = build_re_block(
      M, n_groups[2], how_many_factors_in_random_design[2], ncol_X_random_eff[2],
      group_factor_indexes_for_covariance_2, raw_vec,
      random_effect_sigma_mu[2], random_effect_sigma_sigma[2],
      random_effect_sigma_raw_2, sigma_correlation_factor_2
    );
  }
  if (ncol_X_random_eff[3] > 0) {
    array[ncol_X_random_eff[3]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[3]) raw_vec[i] = random_effect_raw_3[i];
    random_effect_3 = build_re_block(
      M, n_groups[3], how_many_factors_in_random_design[3], ncol_X_random_eff[3],
      group_factor_indexes_for_covariance_3, raw_vec,
      random_effect_sigma_mu[3], random_effect_sigma_sigma[3],
      random_effect_sigma_raw_3, sigma_correlation_factor_3
    );
  }
  if (ncol_X_random_eff[4] > 0) {
    array[ncol_X_random_eff[4]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[4]) raw_vec[i] = random_effect_raw_4[i];
    random_effect_4 = build_re_block(
      M, n_groups[4], how_many_factors_in_random_design[4], ncol_X_random_eff[4],
      group_factor_indexes_for_covariance_4, raw_vec,
      random_effect_sigma_mu[4], random_effect_sigma_sigma[4],
      random_effect_sigma_raw_4, sigma_correlation_factor_4
    );
  }
}

model {
  if (use_data == 1) {
    target += reduce_sum(
      partial_sum_student_t_lupmf,
      array_N,
      grainsize,
      y,
      X,
      beta,
      Xa,
      alpha,
      student_t_df,
      M,
      ncol_X_random_eff,
      X_random_effect_1,
      X_random_effect_2,
      X_random_effect_3,
      X_random_effect_4,
      random_effect_1,
      random_effect_2,
      random_effect_3,
      random_effect_4,
      truncation_not_idx_minimal
    );
  }

  for (a in 1:A) {
    target += abundance_variability_regression_student_t(
      alpha[a],
      beta[variability_to_composition_map[a]],
      prec_intercept_1[a],
      prec_slope_1[a],
      bimodal_mean_variability_association == 1 ? prec_slope_2[a] : 0,
      bimodal_mean_variability_association == 1 ? prec_intercept_2[a] : 0,
      prec_sd[a],
      bimodal_mean_variability_association,
      mix_p_scalar,
      exclude_mean_variability_association
    );
  }

  if (bimodal_mean_variability_association == 1)
    mix_p ~ beta(1, 5);
  else
    mix_p ~ beta(1, 1);

  for (a in 1:A) {
    if (intercept_in_design == 1 && a == 1) {
      prec_intercept[a][1] ~ student_t(3, prior_prec_intercept[1], prior_prec_intercept[2]);
      if (bimodal_mean_variability_association == 1)
        prec_intercept[a][2] ~ student_t(3, prior_prec_intercept[1], prior_prec_intercept[2]);
    } else {
      prec_intercept[a][1] ~ student_t(3, 0, prior_prec_intercept[2]);
      if (bimodal_mean_variability_association == 1)
        prec_intercept[a][2] ~ student_t(3, 0, prior_prec_intercept[2]);
    }
    prec_slope_1[a] ~ student_t(3, prior_prec_slope[1], prior_prec_slope[2]);
    if (bimodal_mean_variability_association == 1)
      prec_slope_2[a] ~ student_t(3, prior_prec_slope[1], prior_prec_slope[2]);
    log_prec_sd[a] ~ normal(prior_prec_sd[1], prior_prec_sd[2]);
  }

  for (c in 1:B_intercept_columns)
    beta[c] ~ normal(prior_mean_intercept[1], prior_mean_intercept[2]);
  if (C > B_intercept_columns) {
    for (c in (B_intercept_columns + 1):C)
      beta[c] ~ normal(prior_mean_coefficients[1], prior_mean_coefficients[2]);
  }

  if (is_random_effect > 0) {
    random_effect_sigma_mu ~ std_normal();
    random_effect_sigma_sigma ~ std_normal();
    zero_random_effect ~ std_normal();
  }

  if (ncol_X_random_eff[1] > 0) {
    for (i in 1:ncol_X_random_eff[1]) random_effect_raw_1[i] ~ std_normal();
    for (m in 1:M) random_effect_sigma_raw_1[m] ~ std_normal();
    for (m in 1:M) sigma_correlation_factor_1[m] ~ lkj_corr_cholesky(2);
  }
  if (ncol_X_random_eff[2] > 0) {
    for (i in 1:ncol_X_random_eff[2]) random_effect_raw_2[i] ~ std_normal();
    for (m in 1:M) random_effect_sigma_raw_2[m] ~ std_normal();
    for (m in 1:M) sigma_correlation_factor_2[m] ~ lkj_corr_cholesky(2);
  }
  if (ncol_X_random_eff[3] > 0) {
    for (i in 1:ncol_X_random_eff[3]) random_effect_raw_3[i] ~ std_normal();
    for (m in 1:M) random_effect_sigma_raw_3[m] ~ std_normal();
    for (m in 1:M) sigma_correlation_factor_3[m] ~ lkj_corr_cholesky(2);
  }
  if (ncol_X_random_eff[4] > 0) {
    for (i in 1:ncol_X_random_eff[4]) random_effect_raw_4[i] ~ std_normal();
    for (m in 1:M) random_effect_sigma_raw_4[m] ~ std_normal();
    for (m in 1:M) sigma_correlation_factor_4[m] ~ lkj_corr_cholesky(2);
  }
}

generated quantities {
  vector[TNS] log_lik = rep_vector(0, TNS);

  if (enable_loo == 1) {
    matrix[N, M] mu = X * beta;
    matrix[N, M] log_precision = Xa * alpha;
    int k = 1;

    if (ncol_X_random_eff[1] > 0)
      mu += X_random_effect_1 * random_effect_1;
    if (ncol_X_random_eff[2] > 0)
      mu += X_random_effect_2 * random_effect_2;
    if (ncol_X_random_eff[3] > 0)
      mu += X_random_effect_3 * random_effect_3;
    if (ncol_X_random_eff[4] > 0)
      mu += X_random_effect_4 * random_effect_4;

    for (n in 1:N) {
      for (m in 1:M) {
        if (pair_is_included(n, m, truncation_not_idx_minimal)) {
          log_lik[k] = student_t_lpdf(
            y[n, m] |
            student_t_df,
            mu[n, m],
            student_t_scale_from_log_precision(
              student_t_df, log_precision[n, m]
            )
          );
          k += 1;
        }
      }
    }
  }
}
