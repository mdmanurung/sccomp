# Independent generative reference; do not call sccomp or Stan RNG helpers here.
simulate_continuous_reference <- function(seed, family = c("student_t", "gaussian"),
                                          design = c("fixed", "intercept", "slope"),
                                          n_features = 20L, n_donors = 20L,
                                          variance_effect = FALSE) {
  family <- match.arg(family)
  design <- match.arg(design)
  stopifnot(n_features >= 2L, n_donors >= 2L)
  set.seed(seed)
  repeats <- if (design == "slope") 4L else 2L
  n <- repeats * n_donors
  donor <- rep(seq_len(n_donors), each = repeats)
  condition <- rep(c(0, 1), n / 2L)
  x <- cbind(`(Intercept)` = 1, conditiontreated = condition)
  w <- if (variance_effect) matrix(1 - 2 * condition, ncol = 1L) else {
    matrix(numeric(), nrow = n, ncol = 0L)
  }
  if (variance_effect) colnames(w) <- "condition1"
  beta <- rbind(stats::rnorm(n_features, 0, 5),
                stats::rnorm(n_features, 0, 2.5))
  gamma <- matrix(stats::rnorm(ncol(w) * n_features),
                   nrow = ncol(w), ncol = n_features)
  log_s0 <- stats::rnorm(1, 0, 1)
  eta0 <- stats::rnorm(1, log(8), 1)
  nu0 <- 2 + exp(eta0)
  sigma_sq <- 1 / stats::rgamma(
    n_features, shape = nu0 / 2, rate = nu0 * exp(2 * log_s0) / 2
  )
  mu <- x %*% beta
  random <- NULL
  if (design != "fixed") {
    k <- if (design == "slope") 2L else 1L
    z <- if (k == 2L) cbind(1, condition) else matrix(1, n, 1L)
    mu_omega <- stats::rnorm(k, log(0.5), 1)
    tau_omega <- abs(stats::rnorm(k, 0, 0.5))
    omega <- exp(sweep(sweep(matrix(stats::rnorm(k * n_features),
                                    nrow = n_features),
                              2L, tau_omega, "*"),
                        2L, mu_omega, "+"))
    correlation <- if (k == 2L) 2 * stats::rbeta(n_features, 2, 2) - 1 else {
      rep(0, n_features)
    }
    b <- array(0, c(n_donors, k, n_features))
    for (g in seq_len(n_features)) {
      corr <- if (k == 2L) matrix(c(1, correlation[g], correlation[g], 1), 2) else {
        matrix(1, 1, 1)
      }
      covariance <- diag(omega[g, ], nrow = k) %*% corr %*%
        diag(omega[g, ], nrow = k)
      b[, , g] <- matrix(stats::rnorm(n_donors * k), nrow = n_donors) %*%
        chol(covariance)
      mu[, g] <- mu[, g] + rowSums(z * matrix(b[donor, , g], nrow = n))
    }
    random <- list(mu_omega = mu_omega, tau_omega = tau_omega, omega = omega,
                   correlation = correlation, b = b)
  }
  residual_variance <- exp(w %*% gamma +
                             matrix(log(sigma_sq), n, n_features, byrow = TRUE))
  eta_obs <- if (family == "student_t") stats::rnorm(1, log(8), 1) else NULL
  nu_obs <- if (family == "student_t") 2 + exp(eta_obs) else NULL
  noise <- if (family == "student_t") {
    stats::rt(n * n_features, nu_obs) * sqrt((nu_obs - 2) / nu_obs)
  } else stats::rnorm(n * n_features)
  y <- mu + matrix(noise, n, n_features) * sqrt(residual_variance)
  samples <- sprintf("sample_%02d", seq_len(n))
  features <- sprintf("feature_%02d", seq_len(n_features))
  data <- data.frame(
    sample = rep(samples, times = n_features),
    feature = rep(features, each = n),
    value = as.vector(y),
    condition = factor(rep(ifelse(condition == 1, "treated", "control"),
                           times = n_features), levels = c("control", "treated")),
    donor = factor(rep(donor, times = n_features))
  )
  list(data = data, seed = seed, family = family, design = design,
       truth = list(beta = beta, gamma = gamma, sigma_sq = sigma_sq,
                    log_s0 = log_s0, eta0 = eta0, nu0 = nu0,
                    eta_obs = eta_obs, nu_obs = nu_obs, random = random,
                    mu = mu, residual_variance = residual_variance),
       x = x, w = w)
}

if ("--self-check" %in% commandArgs(trailingOnly = TRUE)) {
  for (family in c("gaussian", "student_t")) {
    for (design in c("fixed", "intercept", "slope")) {
      for (variance_effect in c(FALSE, TRUE)) {
        a <- simulate_continuous_reference(10102026L, family, design,
                                           variance_effect = variance_effect)
        b <- simulate_continuous_reference(10102026L, family, design,
                                           variance_effect = variance_effect)
        stopifnot(identical(a, b), nrow(a$data) == if (design == "slope") 1600L else 800L,
                  all(is.finite(a$data$value)),
                  all(a$truth$residual_variance > 0), all(colMeans(a$w) == 0))
        if (design != "fixed") {
          stopifnot(all(a$truth$random$tau_omega >= 0),
                    all(a$truth$random$omega > 0),
                    all(abs(a$truth$random$correlation) < 1))
        }
      }
    }
  }
  set.seed(10102026L)
  corrected <- stats::rt(200000L, 8) * sqrt(6 / 8)
  stopifnot(abs(mean(corrected)) < 0.015, abs(stats::var(corrected) - 1) < 0.025)
  cat("Independent simulator determinism, shape, hierarchy and t variance checks passed.\n")
}
