functions {
  #include continuous/functions.stan
}
data {
  matrix[5, 3] y;
  matrix[5, 3] mu;
  matrix[5, 3] log_variance;
  matrix[2, 2] raw;
  vector[2] omega;
  matrix[2, 2] correlation;
  matrix[5, 2] design;
  array[5] int group;
}
generated quantities {
  array[3] int features = {1, 2, 3};
  real gaussian = continuous_partial_sum(features, 1, 3, y, mu,
                                        log_variance, 0, 0);
  real student = continuous_partial_sum(features, 1, 3, y, mu,
                                       log_variance, 1, 7);
  real gaussian_limit = continuous_partial_sum(features, 1, 3, y, mu,
                                              log_variance, 1, 1e7);
  real variance_prior = continuous_variance_prior_lpdf(-0.7 | 9, 0.2);
  real standardized_variance_prior = continuous_standardized_variance_prior_lpdf(-0.7 | 9);
  matrix[3, 2] coordinates = continuous_random_coordinates(raw, omega, correlation,
    y[, 1] - mu[, 1], log_variance[, 1], design, group, 1);
  matrix[2, 2] effects = continuous_random_effect(coordinates[1:2, ], omega, correlation);
  matrix[3, 2] prior_coordinates = continuous_random_coordinates(raw, omega, correlation,
    y[, 1] - mu[, 1], log_variance[, 1], design, group, 0);
  real coordinate_prior = std_normal_lpdf(to_vector(coordinates[1:2, ])) + sum(coordinates[3, ]);
  vector[5] residual = y[, 1] - mu[, 1] - rows_dot_product(design, effects[, group]');
  real coordinate_gaussian_joint = coordinate_prior + normal_lpdf(
    residual | 0, exp(.5 * log_variance[, 1]));
  real coordinate_student_joint = coordinate_prior + student_t_lpdf(
    residual | 7, 0, exp(.5 * log_variance[, 1]) * sqrt(5.0 / 7));
  real prior_coordinate_error = sum(square(to_vector(prior_coordinates[1:2, ] - raw))) +
    sum(square(prior_coordinates[3, ]));
}
