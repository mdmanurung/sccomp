# Synthetic posterior arithmetic only; no sampler or calibration claim.
root <- normalizePath(".")
.libPaths(c(file.path(root, ".local/continuous/library"), .libPaths()))
pkgload::load_all(root, quiet = TRUE, helpers = FALSE)
source("dev/continuous/simulate-fixed-truth.R")
source("dev/continuous/discovery-metrics.R")
source("dev/continuous/fixed-truth-campaign.R")
reference <- simulate_continuous_fixed_truth(10102026L, "variance_alternative", n_features = 12L)
input <- sccomp:::continuous_model_input(reference$data, ~ condition, ~ condition,
  "sample", "feature", "value", "gaussian", "none", "global",
  "estimate", sccomp_continuous_prior())
input$coefficient_offset <- c(1, 0)
variables <- c(as.vector(outer(1:2, 1:12, function(q, g) sprintf("beta[%d,%d]", q, g))),
               sprintf("gamma[1,%d]", 1:12), sprintf("sigma_sq[%d]", 1:12))
set.seed(95102026L)
values <- array(0, c(1000L, 4L, length(variables)), dimnames = list(NULL, NULL, variables))
for (g in 1:12) {
  for (q in 1:2) {
    values[, , sprintf("beta[%d,%d]", q, g)] <- matrix(stats::rnorm(4000, sd = .1), 1000L) +
      (reference$truth$beta[q, g] - if (q == 1L) input$global_center else 0) / input$global_scale
  }
  values[, , sprintf("gamma[1,%d]", g)] <- matrix(stats::rnorm(4000, sd = .03), 1000L) + reference$truth$gamma[1L, g]
  values[, , sprintf("sigma_sq[%d]", g)] <- exp(matrix(stats::rnorm(4000, sd = .01), 1000L))
}
draws <- posterior::as_draws_array(values)
fit <- sccomp:::continuous_estimate_table(draws, input)
class(fit) <- c("sccomp_continuous_tbl", "sccomp_tbl", class(fit))
attr(fit, "fit") <- draws
attr(fit, "continuous") <- list(input = input, inference_method = "hmc", inference_status = "hmc",
                                diagnostics = list(hmc_qualified = TRUE))
tested <- sccomp_test(fit, effect_threshold = 1, variance_ratio_threshold = 1.5,
  contrasts = list(location = c(`treated-control` = "conditiontreated"),
                   log_variance = c(`treated/control` = "-2 * condition1")))
features <- fixed_truth_feature_metrics(fit, tested, reference)
variance <- features[features$component == "log_variance", ]
location <- features[features$component == "location", ]
stopifnot(nrow(features) == 24L,
  identical(variance$truth, as.numeric(-2 * reference$truth$gamma[1L, ])),
  identical(location$truth, as.numeric(reference$truth$beta[2L, ])),
  all(abs(variance$estimate - variance$truth) < .01),
  all(abs(location$estimate - location$truth) < .02),
  all(features$lower_90 >= features$lower), all(features$upper_90 <= features$upper),
  all(features$covered_90), all(features$covered_95))
metrics <- fixed_truth_dataset_metrics(features)
stopifnot(nrow(metrics) == 2L,
  metrics$discoveries[metrics$component == "location"] == 0L,
  metrics$discoveries[metrics$component == "log_variance"] == 6L,
  all(metrics$fdp == 0), all(metrics$false_sign_proportion == 0))
cat("Original-unit posterior extraction, sum-coded variance ratio, 90/95 intervals and discovery arithmetic passed.\n")
