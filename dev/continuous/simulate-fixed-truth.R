# Independent repeated-sampling generator; truths stay fixed across replicates.
simulate_continuous_fixed_truth <- function(
    seed, scenario = "gaussian", design = "fixed", n_features = 20L,
    n_donors = 20L, variance_effect = FALSE,
    hypervariable_fraction = NULL, noise_family = NULL, batch_layout = "none") {
  scenarios <- c("null", "gaussian", "student_t", "contamination", "correlated",
                 "hypervariable", "unbalanced", "few_donors", "non_ig", "trend",
                 "variance_null", "variance_alternative")
  stopifnot(scenario %in% scenarios, design %in% c("fixed", "intercept", "slope"),
            n_features >= 2L, n_donors >= 2L, is.finite(seed))
  stopifnot(batch_layout %in% c("none", "balanced", "partial", "confounded"))
  if (!is.null(hypervariable_fraction)) stopifnot(
    length(hypervariable_fraction) == 1L, is.finite(hypervariable_fraction),
    hypervariable_fraction >= 0, hypervariable_fraction <= 1)
  if (is.null(noise_family)) noise_family <- if (scenario == "student_t") "student_t" else "gaussian"
  stopifnot(length(noise_family) == 1L, noise_family %in% c("gaussian", "student_t"))
  if (scenario == "few_donors") n_donors <- 5L
  repeats <- if (design == "slope") 4L else 2L
  n <- repeats * n_donors
  donor <- rep(seq_len(n_donors), each = repeats)
  condition <- rep(c(0, 1), n / 2L)
  if (scenario == "unbalanced") condition[donor > floor(n_donors / 2)] <- 0
  x <- cbind(`(Intercept)` = 1, conditiontreated = condition)
  w <- matrix(1 - 2 * condition, ncol = 1L)
  w <- sweep(w, 2L, colMeans(w), "-")
  colnames(w) <- "condition1"
  g <- seq_len(n_features)
  set.seed(50102026L)
  beta <- rbind(stats::rnorm(n_features, 3, 2),
                ifelse(g %% 4L == 2L, 1.5, ifelse(g %% 4L == 3L, -1.5, 0)))
  if (scenario %in% c("null", "variance_null", "variance_alternative")) beta[2, ] <- 0
  set.seed(60102026L)
  sigma_sq <- 1 / stats::rgamma(n_features, shape = 5, rate = 5)
  if (scenario == "non_ig") {
    set.seed(60102027L)
    sigma_sq <- exp(stats::rnorm(n_features, 0, .7))
  }
  if (scenario == "hypervariable") sigma_sq[g %% 10L == 0L] <- 9 * sigma_sq[g %% 10L == 0L]
  hypervariable <- rep(FALSE, n_features)
  if (!is.null(hypervariable_fraction)) {
    # Even spacing includes feature 10 in the 20-feature, 20% stress panel.
    count <- floor(n_features * hypervariable_fraction)
    if (count) hypervariable[floor(seq_len(count) * n_features / count)] <- TRUE
    sigma_sq[hypervariable] <- 10 * sigma_sq[hypervariable]
  }
  gamma <- matrix(0, 1L, n_features)
  if (variance_effect || scenario == "variance_alternative") {
    gamma[1, ] <- ifelse(g %% 4L == 2L, -.5 * log(2), ifelse(g %% 4L == 3L, .5 * log(2), 0))
  }
  if (scenario == "trend") {
    true_grand_mean <- as.numeric(colMeans(x %*% beta))
    sigma_sq <- sigma_sq * exp(.5 * as.numeric(scale(true_grand_mean)))
  }
  mu <- x %*% beta
  random <- NULL
  if (design != "fixed") {
    k <- if (design == "slope") 2L else 1L
    omega <- vapply(seq_len(k), function(coefficient) {
      set.seed(70102026L + coefficient)
      exp(stats::rnorm(n_features, log(.5), .25))
    }, numeric(n_features))
    if (k == 2L) omega[, 2L] <- .6 * omega[, 2L]
    rho <- if (k == 2L) rep(.25, n_features) else rep(0, n_features)
    z <- if (k == 2L) cbind(1, condition) else matrix(1, n, 1L)
    b <- array(0, c(n_donors, k, n_features))
    set.seed(seed)
    for (feature in g) {
      corr <- if (k == 2L) matrix(c(1, rho[feature], rho[feature], 1), 2) else matrix(1, 1, 1)
      covariance <- diag(omega[feature, ], nrow = k) %*% corr %*% diag(omega[feature, ], nrow = k)
      b[, , feature] <- matrix(stats::rnorm(n_donors * k), n_donors, k) %*% chol(covariance)
      mu[, feature] <- mu[, feature] + rowSums(z * matrix(b[donor, , feature], n))
    }
    random <- list(omega = omega, correlation = rho, b = b)
  }
  set.seed(seed + 100000L)
  noise <- matrix(if (noise_family == "student_t") {
    stats::rt(n * n_features, 5) * sqrt(3 / 5)
  } else stats::rnorm(n * n_features), n, n_features)
  if (scenario == "correlated") {
    for (feature in 2:n_features) {
      noise[, feature] <- .6 * noise[, feature - 1L] + sqrt(1 - .6^2) * noise[, feature]
    }
  }
  contamination <- matrix(FALSE, n, n_features)
  mixture_variance <- 1
  if (scenario == "contamination") {
    contamination[] <- stats::runif(length(noise)) < .03
    noise <- noise + contamination * sample(c(-8, 8), length(noise), replace = TRUE)
    mixture_variance <- 1 + .03 * 8^2
  }
  base_variance <- exp(w %*% gamma + matrix(log(sigma_sq), n, n_features, byrow = TRUE))
  batch <- batch_beta <- batch_gamma <- NULL
  if (batch_layout != "none") {
    batch <- switch(batch_layout,
      balanced = rep(c(0, 0, 1, 1), length.out = n),
      partial = ifelse(donor %% 4L == 0L, 1 - condition, condition),
      confounded = condition)
    batch_beta <- ifelse(g %% 2L == 0L, .8, -.8)
    batch_gamma <- rep(-.5 * log(2), n_features)
    mu <- mu + outer(batch, batch_beta)
    batch_w <- (1 - 2 * batch) - mean(1 - 2 * batch)
    base_variance <- base_variance * exp(outer(batch_w, batch_gamma))
  }
  y <- mu + sqrt(base_variance) * noise
  data <- data.frame(sample = rep(sprintf("sample_%02d", seq_len(n)), n_features),
    feature = rep(sprintf("feature_%04d", g), each = n), value = as.vector(y),
    condition = factor(rep(ifelse(condition == 1, "treated", "control"), n_features),
                         levels = c("control", "treated")),
    donor = factor(rep(donor, n_features)))
  if (!is.null(batch)) data$batch <- factor(rep(ifelse(batch == 1, "B", "A"), n_features),
                                            levels = c("A", "B"))
  list(data = data, scenario = scenario, design = design, seed = seed,
    truth = list(beta = beta, gamma = gamma, sigma_sq = sigma_sq,
      mu = mu, base_variance = base_variance,
      observation_variance = base_variance * mixture_variance, random = random,
      contamination = contamination, hypervariable = hypervariable,
      batch_beta = batch_beta, batch_gamma = batch_gamma), x = x, w = w)
}

if ("--self-check" %in% commandArgs(trailingOnly = TRUE)) {
  for (scenario in c("null", "gaussian", "student_t", "contamination", "correlated",
                      "hypervariable", "unbalanced", "few_donors", "non_ig", "trend",
                      "variance_null", "variance_alternative")) {
    for (design in c("fixed", "intercept", "slope")) {
      a <- simulate_continuous_fixed_truth(10102026L, scenario, design)
      b <- simulate_continuous_fixed_truth(10102027L, scenario, design)
      stopifnot(all(is.finite(a$data$value)), all(a$truth$observation_variance > 0),
                identical(a$truth$beta, b$truth$beta), identical(a$truth$sigma_sq, b$truth$sigma_sq),
                identical(a$truth$gamma, b$truth$gamma), !identical(a$data$value, b$data$value),
                all(abs(colMeans(a$w)) < 1e-12))
      if (scenario == "null") stopifnot(all(a$truth$beta[2, ] == 0))
      if (scenario == "few_donors") stopifnot(length(unique(a$data$donor)) == 5L)
      if (scenario == "unbalanced") {
        metadata <- a$data[!duplicated(a$data$sample), ]
        counts <- table(metadata$condition)
        stopifnot(counts["control"] == 3L * counts["treated"],
          any(vapply(split(metadata$condition, metadata$donor), function(x) {
            length(unique(x)) == 2L
          }, logical(1))))
      }
      if (scenario == "variance_alternative") {
        stopifnot(any(abs(-2 * a$truth$gamma) > log(1.5)))
      }
    }
  }
  for (design in c("fixed", "intercept", "slope")) {
    a <- simulate_continuous_fixed_truth(10102026L, design = design, n_features = 100L)
    b <- simulate_continuous_fixed_truth(10102026L, design = design, n_features = 500L)
    stopifnot(identical(a$data$value, b$data$value[seq_len(nrow(a$data))]))
  }
  null <- simulate_continuous_fixed_truth(10102026L, "null", n_features = 1000L, n_donors = 100L)
  standardized <- (matrix(null$data$value, 200L) - null$truth$mu) / sqrt(null$truth$base_variance)
  stopifnot(abs(mean(standardized)) < .01, abs(stats::var(as.vector(standardized)) - 1) < .015)
  cat("Fixed-truth invariance, misspecification panels, nulls and Gaussian moments passed.\n")
}
