root <- normalizePath(".")
.libPaths(c(file.path(root, ".local/continuous/library"), .libPaths()))
run_id <- commandArgs(trailingOnly = TRUE)
if (!length(run_id)) run_id <- format(Sys.time(), "densities-%Y%m%dT%H%M%S", tz = "UTC")
stopifnot(length(run_id) == 1L, grepl("^[a-zA-Z0-9_-]+$", run_id))
run <- file.path(root, "dev/continuous/runs", run_id)
if (dir.exists(run)) stop("Run exists; never overwrite a density receipt.")
dir.create(run, recursive = TRUE)
model_dir <- file.path(root, ".local/continuous/density-model")
dir.create(file.path(model_dir, "continuous"), recursive = TRUE,
           showWarnings = FALSE)
stopifnot(file.copy("inst/stan/continuous/functions.stan",
                    file.path(model_dir, "continuous"), overwrite = TRUE))
stopifnot(file.copy("dev/continuous/density-reference.stan", model_dir,
                    overwrite = TRUE))
model <- cmdstanr::cmdstan_model(
  file.path(model_dir, "density-reference.stan"), include_paths = model_dir,
  cpp_options = list(stan_threads = TRUE), force_recompile = TRUE
)
y <- matrix(seq(-3, 4, length.out = 15), 5, 3)
mu <- matrix(seq(-1, 1, length.out = 15), 5, 3)
log_variance <- matrix(seq(-0.5, 0.8, length.out = 15), 5, 3)
raw <- matrix(c(-.3, .4, .2, -.7), 2, 2)
omega <- c(1.3, .4)
correlation <- matrix(c(1, .3, 0, sqrt(1 - .3^2)), 2, 2)
design <- cbind(1, c(0, 1, -.5, .5, 2))
group <- c(1L, 1L, 1L, 2L, 2L)
fit <- model$sample(data = list(y = y, mu = mu, log_variance = log_variance,
                              raw = raw, omega = omega, correlation = correlation,
                              design = design, group = group),
                    chains = 1, threads_per_chain = 1,
                    fixed_param = TRUE, iter_sampling = 1,
                    seed = 10102026L, output_dir = run, sig_figs = 17)
observed_matrix <- fit$draws(c("gaussian", "student", "gaussian_limit",
                         "variance_prior", "standardized_variance_prior",
                         "coordinate_prior", "coordinate_gaussian_joint",
                         "coordinate_student_joint", "prior_coordinate_error"), format = "matrix")
observed <- stats::setNames(as.numeric(observed_matrix[1L, ]), colnames(observed_matrix))
sd <- exp(log_variance / 2)
t_scale <- sd * sqrt(5 / 7)
shape <- 9 / 2
rate <- 9 * exp(0.4) / 2
variance <- exp(-0.7)
expected <- c(
  gaussian = sum(stats::dnorm(y, mu, sd, log = TRUE)),
  student = sum(stats::dt((y - mu) / t_scale, 7, log = TRUE) - log(t_scale)),
  gaussian_limit = sum(stats::dnorm(y, mu, sd, log = TRUE)),
  variance_prior = shape * log(rate) - lgamma(shape) -
    (shape + 1) * log(variance) - rate / variance + log(variance),
  standardized_variance_prior = shape * log(shape) - lgamma(shape) -
    shape * (-.7 * sqrt(2 / 9)) - shape * exp(.7 * sqrt(2 / 9)) + log(sqrt(2 / 9))
)
multivariate_normal <- function(x, covariance) {
  -.5 * length(x) * log(2 * pi) - .5 * as.numeric(determinant(covariance, logarithm = TRUE)$modulus) -
    .5 * as.numeric(crossprod(x, solve(covariance, x)))
}
scale_matrix <- diag(omega) %*% correlation
covariance <- tcrossprod(scale_matrix)
coordinate_prior <- gaussian_marginal <- gaussian_likelihood <- student_likelihood <- 0
for (j in 1:2) {
  indices <- which(group == j)
  z <- design[indices, , drop = FALSE]
  variance <- exp(log_variance[indices, 1])
  residual <- y[indices, 1] - mu[indices, 1]
  precision <- diag(2) + crossprod(z %*% scale_matrix, (z %*% scale_matrix) / variance)
  information <- crossprod(z %*% scale_matrix, residual / variance)
  effect <- scale_matrix %*% (solve(precision, information) + solve(chol(precision), raw[, j]))
  jacobian <- as.numeric(determinant(scale_matrix, logarithm = TRUE)$modulus) -
    sum(log(diag(chol(precision))))
  coordinate_prior <- coordinate_prior + multivariate_normal(effect, covariance) + jacobian
  adjusted <- residual - as.vector(z %*% effect)
  gaussian_likelihood <- gaussian_likelihood + sum(stats::dnorm(adjusted, 0, sqrt(variance), log = TRUE))
  student_scale <- sqrt(variance * 5 / 7)
  student_likelihood <- student_likelihood + sum(stats::dt(adjusted / student_scale, 7, log = TRUE) - log(student_scale))
  gaussian_marginal <- gaussian_marginal + multivariate_normal(
    residual, diag(variance) + z %*% covariance %*% t(z))
}
expected <- c(expected, coordinate_prior = coordinate_prior,
  coordinate_gaussian_joint = gaussian_marginal + sum(stats::dnorm(raw, log = TRUE)),
  coordinate_student_joint = coordinate_prior + student_likelihood,
  prior_coordinate_error = 0)
stopifnot(abs(coordinate_prior + gaussian_likelihood - expected["coordinate_gaussian_joint"]) < 1e-10)
stopifnot(setequal(names(observed), names(expected)))
observed <- observed[names(expected)]
result <- data.frame(parameter = names(expected), observed = unname(observed),
                      expected = unname(expected),
                      absolute_error = unname(abs(observed - expected)))
utils::write.csv(result, file.path(run, "density-comparison.csv"),
                 row.names = FALSE)
writeLines(capture.output(tools::md5sum(c(
  "inst/stan/continuous/functions.stan", "dev/continuous/density-reference.stan"
))), file.path(run, "source-md5.txt"))
stopifnot(is.numeric(result$absolute_error), length(result$absolute_error) == 9L,
          all(is.finite(result$absolute_error)),
          all(result$absolute_error[-3L] < 1e-10),
          result$absolute_error[3] < 1e-5)
cat("Production likelihood and log-variance-prior densities match independent R calculations.\n")
