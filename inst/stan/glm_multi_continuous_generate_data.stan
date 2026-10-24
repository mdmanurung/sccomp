functions {
  #include continuous/functions.stan
}
data {
  #include continuous/data.stan
}
transformed data {
  array[G] int features;
  for (g in 1:G) features[g] = g;
  if (estimate_df == 1 && student_t == 0)
    reject("Gaussian family cannot estimate Student-t df");
}
parameters {
  #include continuous/parameters.stan
}
transformed parameters {
  #include continuous/transformed_parameters.stan
}
model {}
generated quantities {
  matrix[N, G] mu = X * beta;
  matrix[N, G] residual_variance =
    exp(W * gamma + rep_matrix(log_sigma_sq', N));
  matrix[N, G] y_rep;
  matrix[N, G] log_lik;
  #include continuous/location.stan
  for (n in 1:N) {
    for (g in 1:G) {
      if (student_t == 1) {
        real scale = continuous_student_scale(residual_variance[n, g], nu_obs);
        y_rep[n, g] = student_t_rng(nu_obs, mu[n, g], scale);
        log_lik[n, g] = student_t_lpdf(y[n, g] | nu_obs, mu[n, g], scale);
      } else {
        real sd = sqrt(residual_variance[n, g]);
        y_rep[n, g] = normal_rng(mu[n, g], sd);
        log_lik[n, g] = normal_lpdf(y[n, g] | mu[n, g], sd);
      }
    }
  }
}
