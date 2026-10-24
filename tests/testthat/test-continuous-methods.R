# Exact-draw fixtures test API arithmetic; they are not calibration evidence.
continuous_draw_fixture <- function(family = "gaussian", random = FALSE) {
  n <- if (random) 24L else 12L
  data <- expand.grid(sample = sprintf("sample_%02d", seq_len(n)),
                       feature = sprintf("feature_%02d", seq_len(12)),
                       stringsAsFactors = FALSE)
  data$value <- sin(seq_len(nrow(data)))
  data$condition <- factor(rep(c("control", "treated"), 6L * n))
  data$donor <- factor(rep(rep(seq_len(6), each = n / 6L), 12L))
  input <- continuous_model_input(
    data, if (random) ~ condition + (1 + condition | donor) else ~ condition,
    ~ condition, "sample", "feature", "value", family, "none", "none",
    "estimate", sccomp_continuous_prior()
  )
  input$global_center <- 10
  input$global_scale <- 2
  input$coefficient_offset <- c(1, 0)
  names <- c(outer(seq_len(2), seq_len(12), function(c, g) sprintf("beta[%d,%d]", c, g)),
             sprintf("gamma[1,%d]", seq_len(12)), sprintf("sigma_sq[%d]", seq_len(12)))
  if (family == "student_t") names <- c(names, "nu_obs")
  if (random) {
    names <- c(names, as.vector(outer(seq_len(12), seq_len(2), function(g, q) sprintf("omega_1[%d,%d]", g, q))))
    for (g in seq_len(12)) {
      for (q in seq_len(2)) {
        names <- c(names, sprintf("random_effect_1[%d,%d,%d]", g, q, seq_len(6)))
        names <- c(names, sprintf("L_1[%d,%d,%d]", g, q, seq_len(2)))
      }
    }
  }
  array <- array(0, c(1000, 4, length(names)), dimnames = list(NULL, NULL, names))
  for (g in seq_len(12)) {
    array[, , sprintf("beta[1,%d]", g)] <- -3
    array[, , sprintf("beta[2,%d]", g)] <- if (g == 1) 2 else 0
    array[, , sprintf("gamma[1,%d]", g)] <- if (g == 1) -.4 else 0
    array[, , sprintf("sigma_sq[%d]", g)] <- 3
    if (random) {
      array[, , sprintf("omega_1[%d,1]", g)] <- .5
      array[, , sprintf("omega_1[%d,2]", g)] <- .25
      array[, , sprintf("L_1[%d,1,1]", g)] <- 1
      array[, , sprintf("L_1[%d,2,1]", g)] <- .6
      array[, , sprintf("L_1[%d,2,2]", g)] <- .8
      for (j in seq_len(6)) {
        array[, , sprintf("random_effect_1[%d,1,%d]", g, j)] <- j / 10
        array[, , sprintf("random_effect_1[%d,2,%d]", g, j)] <- j / 20
      }
    }
  }
  if (family == "student_t") array[, , "nu_obs"] <- 8
  draws <- posterior::as_draws_array(array)
  fit <- continuous_estimate_table(draws, input)
  class(fit) <- c("sccomp_continuous_tbl", "sccomp_tbl", class(fit))
  attr(fit, "fit") <- draws
  attr(fit, "continuous") <- list(input = input, family = family, seed = 1001,
    inference_method = "hmc", inference_status = "hmc", formula_location = ~ condition,
    formula_variance = ~ condition, diagnostics = list(hmc_qualified = TRUE))
  attr(fit, "inference_method") <- "hmc"
  fit
}

test_that("continuous tables back-transform location, variance and random effects", {
  fit <- continuous_draw_fixture(random = TRUE)
  expect_equal(fit$estimate[fit$component == "location" & fit$parameter == "(Intercept)"], rep(4, 12))
  expect_equal(fit$estimate[fit$component == "residual_variance"], rep(12, 12))
  expect_equal(fit$estimate[fit$component == "random_sd"], rep(c(1, .5), 12))
  expect_equal(fit$estimate[fit$component == "random_correlation"], rep(.6, 12))
  expect_named(fit, c("feature", "component", "parameter", "estimate", "lower", "upper", "rhat", "ess_bulk", "ess_tail"))
})

test_that("hyperparameter log SDs return to assay units without rescaling df or dispersion", {
  variables <- c("log_s0", "nu0", "mu_omega_1[1]", "tau_omega_1[1]")
  set.seed(20261010L)
  values <- array(stats::rnorm(1600), c(100, 4, 4),
                  dimnames = list(NULL, NULL, variables))
  draws <- posterior::as_draws_array(values)
  original <- unclass(draws)
  a <- continuous_hyperparameter_table(draws, variables, 1)
  b <- continuous_hyperparameter_table(draws, variables, 3)
  shifted <- a$variable %in% c("log_s0", "mu_omega_1[1]")
  expect_equal(b$mean - a$mean, as.numeric(shifted) * log(3), tolerance = 1e-12)
  expect_equal(b$sd, a$sd, tolerance = 1e-12)
  expect_equal(b$rhat, a$rhat, tolerance = 1e-12)
  expect_identical(unclass(draws), original)
})

test_that("contrast arithmetic rejects nonlinear expressions and code execution", {
  columns <- c("a", "b:c")
  expect_equal(continuous_contrast_weights("(a + 2 * `b:c`) / 2", columns), c(a = .5, `b:c` = 1))
  expect_equal(continuous_contrast_weights("-a", columns), c(a = -1, `b:c` = 0))
  for (expression in c("a * a", "a / a", "a / 0", "a + 1", "system('false')", "a; a", "missing", "a-a")) {
    expect_error(continuous_contrast_weights(expression, columns))
  }
})

test_that("combined contrasts cannot inherit passing coefficient diagnostics", {
  fit <- continuous_draw_fixture()
  set.seed(20261010L)
  common <- matrix(stats::rnorm(4000, sd = 50), 1000, 4)
  slow <- matrix(rep(c(-.1, .1, -.1, .1), each = 1000), 1000, 4) +
    matrix(stats::rnorm(4000, sd = .01), 1000, 4)
  draws <- attr(fit, "fit")
  draws[, , "beta[1,1]"] <- common + slow
  draws[, , "beta[2,1]"] <- common - slow
  attr(fit, "fit") <- draws
  for (values in list(common + slow, common - slow)) {
    expect_lt(posterior::rhat(values), 1.01)
    expect_gte(posterior::ess_bulk(values), 400)
  }
  expect_error(sccomp_test(fit, contrasts = "`(Intercept)` - conditiontreated",
                          effect_threshold = 1), "diagnostics failed for requested contrast")
  expect_error(sccomp_test(fit, contrasts = "1e308 * conditiontreated",
                          effect_threshold = 1), "Contrast values must be finite")
})

test_that("continuous discovery requires units and valid inference", {
  fit <- continuous_draw_fixture()
  expect_error(sccomp_test(fit), "effect_threshold")
  expect_error(sccomp_test(fit, effect_threshold = -1), "positive")
  expect_error(sccomp_test(fit, effect_threshold = 1, posterior_fdr = 5), "zero and one")
  expect_error(sccomp_test(fit, percent_false_positive = 1), "effect_threshold")
  expect_error(sccomp_test(fit, test_composition_above_logit_fold_change = 1), "effect_threshold")
  expect_equal(sccomp_test(fit, NULL, effect_threshold = 1)$selected,
               sccomp_test(fit, effect_threshold = 1)$selected)
  expect_error(sccomp_test(fit, effect_threshold = 1, variance_ratio_threshold = 1), "greater than one")
  meta <- attr(fit, "continuous")
  meta$diagnostics$hmc_qualified <- FALSE
  attr(fit, "continuous") <- meta
  expect_error(sccomp_test(fit, effect_threshold = 1), "diagnostics")
  meta$inference_method <- "variational"
  meta$inference_status <- "exploratory"
  attr(fit, "continuous") <- meta
  expect_error(sccomp_test(fit, effect_threshold = 1), "allow_exploratory")
  result <- sccomp_test(fit, effect_threshold = 1, allow_exploratory = TRUE)
  expect_true(all(result$inference_status == "exploratory"))
  expect_true(all(is.na(result$probability_mcse)))
  meta$inference_status <- "prior_only"
  attr(fit, "continuous") <- meta
  expect_error(sccomp_test(fit, effect_threshold = 1, allow_exploratory = TRUE), "Prior-only")
})

test_that("precise small effects do not pass a biological effect threshold", {
  fit <- continuous_draw_fixture()
  draws <- attr(fit, "fit")
  set.seed(20261013L)
  draws[, , "beta[2,1]"] <- matrix(stats::rnorm(4000, .1, .005), 1000, 4)
  draws[, , "beta[2,2]"] <- matrix(stats::rnorm(4000, 1, .05), 1000, 4)
  attr(fit, "fit") <- draws
  tested <- sccomp_test(fit, effect_threshold = 1)
  small <- tested[tested$feature == "feature_01", ]
  large <- tested[tested$feature == "feature_02", ]
  expect_lt(small$upper - small$lower, .1)
  expect_gt(small$lower, 0)
  expect_equal(small$p_null, 1)
  expect_equal(small$p_false_sign, 0)
  expect_equal(small$p_error, 1)
  expect_false(small$selected)
  expect_true(large$selected)
})

test_that("posterior error families separate components and contrasts", {
  fit <- continuous_draw_fixture()
  result <- sccomp_test(fit, contrasts = list(location = c(forward = "conditiontreated", reverse = "-conditiontreated"),
      log_variance = c(treated_control = "-2 * condition1")),
    effect_threshold = 1, variance_ratio_threshold = 1.5)
  expect_equal(nrow(result), 36)
  expect_equal(result$estimate[result$feature == "feature_01"], c(4, -4, .8))
  expect_equal(result$posterior_fdr[result$feature == "feature_01"], c(0, 0, 0))
  expect_equal(sum(result$selected), 3)
  expect_true(all(result$p_error[result$feature != "feature_01"] == 1))
  expect_equal(result$posterior_fdr[result$parameter == "forward"], c(0, rep(11 / 12, 11)))
  expect_error(sccomp_test(fit, contrasts = list(log_variance = "condition1"), effect_threshold = 1), "Variance contrasts require")
  expect_error(sccomp_test(fit, contrasts = list(log_variance = "conditiontreated"), effect_threshold = 1,
                          variance_ratio_threshold = 2), "Unknown contrast coefficient")
  stripped <- sccomp_test(fit, effect_threshold = 1, pass_fit = FALSE)
  expect_null(attr(stripped, "fit"))
  expect_error(sccomp_predict(stripped), "no retained")
})

test_that("tied posterior error probabilities cannot select arbitrary features", {
  fit <- continuous_draw_fixture()
  draws <- attr(fit, "fit")
  set.seed(20261010L)
  for (g in 2:3) {
    values <- rep(0, 4000L)
    values[sample.int(4000L, 3640L)] <- 2
    draws[, , sprintf("beta[2,%d]", g)] <- matrix(values, 1000L)
  }
  attr(fit, "fit") <- draws
  result <- sccomp_test(fit, effect_threshold = 1)
  expect_equal(result$p_error[2:3], c(.09, .09))
  expect_equal(result$posterior_fdr[2:3], c(.06, .06))
  expect_equal(which(result$selected), 1L)
})

test_that("continuous predictions preserve original units, coding, draws and RNG", {
  fit <- continuous_draw_fixture()
  set.seed(481)
  before <- .Random.seed
  predictions <- sccomp_predict(fit, number_of_draws = 10, mcmc_seed = 2001,
                                summary_instead_of_draws = FALSE)
  expect_identical(.Random.seed, before)
  expect_identical(predictions, sccomp_replicate(fit, number_of_draws = 10, mcmc_seed = 2001))
  control <- predictions$sample == "sample_01" & predictions$feature == "feature_01"
  treated <- predictions$sample == "sample_02" & predictions$feature == "feature_01"
  expect_equal(predictions$location[control], rep(4, 10))
  expect_equal(predictions$location[treated], rep(8, 10))
  expect_equal(predictions$residual_variance[control], rep(12 * exp(-.4), 10))
  expect_equal(predictions$residual_variance[treated], rep(12 * exp(.4), 10))
  new <- attr(fit, "continuous")$input$sample_data[2, , drop = FALSE]
  single <- sccomp_predict(fit, new_data = new, number_of_draws = 10,
                           mcmc_seed = 2001, summary_instead_of_draws = FALSE)
  expect_equal(single$residual_variance[single$feature == "feature_01"], rep(12 * exp(.4), 10))
  new$condition <- "unknown"
  expect_error(sccomp_predict(fit, new_data = new), "new level")
  expect_error(sccomp_predict(fit, formula_composition = ~ condition), "formula_location")
  expect_error(sccomp_predict(fit, robust = TRUE), "robust")
  expect_error(sccomp_replicate(fit, formula_variability = ~ condition), "formula_location")
  residuals <- sccomp_calculate_residuals(fit)
  expect_equal(residuals$residual, residuals$observed_value - residuals$location)
  expect_equal(residuals$standardized_residual, residuals$residual / sqrt(residuals$residual_variance))
})

test_that("Student-t prediction uses residual variance rather than scale squared", {
  fit <- continuous_draw_fixture("student_t")
  result <- sccomp_replicate(fit, number_of_draws = 4000, mcmc_seed = 7821)
  control <- result$feature != "feature_01"
  residuals <- result$generated_value[control] - result$location[control]
  expect_equal(stats::var(residuals), 12, tolerance = .03)
})

test_that("continuous random prediction and unsupported APIs dispatch safely", {
  fit <- continuous_draw_fixture(random = TRUE)
  result <- sccomp_predict(fit, number_of_draws = 1, mcmc_seed = 101,
                           summary_instead_of_draws = FALSE)
  expect_equal(result$location[result$feature == "feature_01" & result$sample == "sample_08"], 8.6)
  new <- attr(fit, "continuous")$input$sample_data[1, , drop = FALSE]
  new$donor <- 999
  expect_error(sccomp_predict(fit, new_data = new), "observed random-effect")
  expect_error(sccomp_remove_outliers(fit), "unavailable")
  expect_error(sccomp_proportional_fold_change(fit), "unavailable")
  expect_error(sccomp_remove_unwanted_effects(fit), "unavailable")
  expect_output(print(fit), "sccomp continuous")
  expect_s3_class(plot(fit), "ggplot")
})

test_that("continuous prediction formulas remove fitted terms in input units", {
  fit <- continuous_draw_fixture(random = TRUE)
  population <- sccomp_predict(fit, formula_location = ~ condition,
    number_of_draws = 1, mcmc_seed = 1001, summary_instead_of_draws = FALSE)
  expect_equal(population$location[population$feature == "feature_01" &
                                   population$sample == "sample_04"], 8)
  baseline <- sccomp_predict(fit, formula_location = ~ 1,
    number_of_draws = 1, mcmc_seed = 1001, summary_instead_of_draws = FALSE)
  expect_equal(baseline$location, rep(4, nrow(baseline)))
  no_intercept <- sccomp_predict(fit, formula_location = ~ 0 + condition,
    number_of_draws = 1, mcmc_seed = 1001, summary_instead_of_draws = FALSE)
  expect_equal(no_intercept$location[no_intercept$feature == "feature_01" &
                                     no_intercept$sample == "sample_04"], 4)
  constant_variance <- sccomp_replicate(fit, formula_variance = ~ 1,
                                        number_of_draws = 1, mcmc_seed = 1001)
  expect_equal(constant_variance$residual_variance, rep(12, nrow(constant_variance)))
  expect_error(sccomp_predict(fit, formula_location = ~ unknown), "fitted terms")
  expect_error(sccomp_predict(fit, formula_location = ~ condition + (1 | donor)), "match fitted clauses")
  full <- sccomp_predict(fit, formula_location = ~ condition + (1 + condition | donor),
    number_of_draws = 1, mcmc_seed = 1001, summary_instead_of_draws = FALSE)
  expect_equal(full, sccomp_predict(fit, number_of_draws = 1, mcmc_seed = 1001,
                                    summary_instead_of_draws = FALSE))
})
