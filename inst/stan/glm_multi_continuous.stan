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
model {
  matrix[N, G] mu = X * beta;
  matrix[N, G] log_variance = W * gamma + rep_matrix(log_sigma_sq', N);
  #include continuous/location.stan
  #include continuous/priors.stan
  if (use_data == 1)
    target += reduce_sum(continuous_partial_sum, features, grainsize,
                         y, mu, log_variance, student_t, nu_obs);
}
