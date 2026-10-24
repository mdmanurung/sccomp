# Matching continuous values/designs; no count-specific weights or transformations.
fit_continuous_comparator <- function(data, design, effect_threshold = 1, fdr = .05,
                                      trend = FALSE, robust = FALSE,
                                      formula_location = NULL) {
  stopifnot(design %in% c("fixed", "intercept", "slope"),
            is.finite(effect_threshold), effect_threshold > 0, fdr > 0, fdr < 1,
            is.logical(trend), is.logical(robust),
            !anyDuplicated(data[c("sample", "feature")]))
  samples <- unique(as.character(data$sample))
  features <- unique(as.character(data$feature))
  stopifnot(nrow(data) == length(samples) * length(features))
  if (is.null(formula_location)) formula_location <- switch(design,
    fixed = ~ condition, intercept = ~ condition + (1 | donor),
    slope = ~ condition + (1 + condition | donor))
  metadata <- data[match(samples, data$sample),
                   unique(c("sample", "condition", "donor", all.vars(formula_location)))]
  rownames(metadata) <- samples
  y <- matrix(NA_real_, length(features), length(samples),
              dimnames = list(features, samples))
  y[cbind(match(data$feature, features), match(data$sample, samples))] <- data$value
  stopifnot(all(is.finite(y)))
  coefficient <- "conditiontreated"
  if (design == "fixed") {
    x <- stats::model.matrix(formula_location, metadata)
    raw_fit <- limma::lmFit(y, x)
    fit <- limma::eBayes(raw_fit, trend = trend, robust = robust)
    intervals <- lapply(c(.90, .95), function(level) limma::topTable(
      fit, coef = coefficient, number = Inf, sort.by = "none", confint = level))
    decision <- limma::topTreat(limma::treat(raw_fit, lfc = effect_threshold, trend = trend, robust = robust),
                               coef = coefficient, number = Inf, sort.by = "none")
  } else {
    formula <- formula_location
    # dream sets options(warn = 2), so lmerTest's Satterthwaite Hessian warning on boundary
    # (singular) REML fits drops the feature. Muffle only that warning; estimates are kept and
    # flagged in `boundary`. Every other warning and dream check keeps its default.
    muffle_hessian <- function(w) {
      if (grepl("^Model failed to converge with [0-9]+ negative eigenvalue", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
    raw_fit <- withCallingHandlers(variancePartition::dream(y, formula, metadata,
      ddf = "Satterthwaite", hideErrorsInBackend = FALSE,
      BPPARAM = BiocParallel::SerialParam()), warning = muffle_hessian)
    boundary <- vapply(features, function(feature) {
      lme4::isSingular(suppressMessages(suppressWarnings(
        lme4::lmer(update(formula, y ~ .), transform(metadata, y = y[feature, ]), REML = TRUE))))
    }, logical(1))
    # limma's trend covariate is the feature mean (what lmFit stores for fixed designs).
    if (is.null(raw_fit$Amean)) raw_fit$Amean <- rowMeans(y)[rownames(raw_fit$coefficients)]
    fit <- variancePartition::eBayes(raw_fit, trend = trend, robust = robust)
    intervals <- lapply(c(.90, .95), function(level) variancePartition::topTable(
      fit, coef = coefficient, number = Inf, sort.by = "none", confint = level))
    decision <- variancePartition::getTreat(fit, lfc = effect_threshold,
      coef = coefficient, number = Inf, sort.by = "none")
  }
  if (design == "fixed") boundary <- rep(NA, length(features))
  stopifnot(all(features %in% rownames(decision)),
            all(vapply(intervals, function(x) all(features %in% rownames(x)), logical(1))))
  intervals <- lapply(intervals, function(x) x[features, , drop = FALSE])
  decision <- decision[features, , drop = FALSE]
  result <- data.frame(feature = features, component = "location",
    parameter = coefficient, method = paste0(if (design == "fixed") "limma" else "dream",
                    if (trend) "_trend", if (robust) "_robust"),
    estimate = intervals[[2L]]$logFC,
    lower_90 = intervals[[1L]]$CI.L, upper_90 = intervals[[1L]]$CI.R,
    lower_95 = intervals[[2L]]$CI.L, upper_95 = intervals[[2L]]$CI.R,
    threshold = effect_threshold, p_value = decision$P.Value,
    adjusted_p_value = decision$adj.P.Val, selected = decision$adj.P.Val <= fdr,
    direction = ifelse(decision$logFC >= 0, "positive", "negative"),
    boundary = unname(boundary))
  result$valid <- is.finite(result$estimate) & is.finite(result$lower_95) &
    is.finite(result$upper_95) & is.finite(result$adjusted_p_value)
  list(results = result, fit = fit, raw_fit = raw_fit,
       inputs = list(samples = samples, features = features, metadata = metadata, y = y),
       decision_rule = "TREAT minimum-effect test with BH adjustment; distinct from posterior directional error.",
       settings = list(effect_threshold = effect_threshold, fdr = fdr, trend = trend, robust = robust,
                       formula_location = formula_location,
                       ddf = if (design == "fixed") NULL else "Satterthwaite",
                       hessian_warning = if (design == "fixed") NULL else "negative-eigenvalue warning muffled; boundary fits flagged",
                       variancePartition = as.character(utils::packageVersion("variancePartition"))))
}

# This tests spread association at zero, not a minimum variance ratio.
fit_continuous_variance_comparator <- function(comparison, fdr = .05) {
  stopifnot(is.finite(fdr), fdr > 0, fdr < 1)
  raw <- comparison$raw_fit
  if (is.null(raw$formula)) {
    raw <- variancePartition::dream(comparison$inputs$y,
      comparison$settings$formula_location, comparison$inputs$metadata,
      BPPARAM = BiocParallel::SerialParam(), hideErrorsInBackend = FALSE)
  }
  fit <- variancePartition::diffVar(raw, method = "AD", scale = "leverage",
    BPPARAM = BiocParallel::SerialParam(), hideErrorsInBackend = FALSE)
  table <- variancePartition::topTable(fit, coef = "conditiontreated",
                                        number = Inf, sort.by = "none")
  features <- comparison$inputs$features
  stopifnot(all(features %in% rownames(table)))
  table <- table[features, , drop = FALSE]
  data.frame(feature = features, component = "spread_association",
    method = "dream_diffVar_AD", estimate = table$logFC,
    p_value = table$P.Value, adjusted_p_value = table$adj.P.Val,
    selected = table$adj.P.Val <= fdr,
    direction = ifelse(table$logFC >= 0, "positive", "negative"),
    valid = is.finite(table$logFC) & is.finite(table$adj.P.Val))
}

# Per-dataset comparator metrics against the simulation truth (95% bounds plus 90% coverage).
comparator_dataset_metrics <- function(results, reference, n_features) {
  stopifnot(all(results$valid), nrow(results) == n_features)
  results$lower <- results$lower_95
  results$upper <- results$upper_95
  truth <- stats::setNames(reference$truth$beta[2L, ], unique(reference$data$feature))
  metrics <- continuous_discovery_metrics(results, truth, 1)
  metrics$coverage_90 <- mean(truth[results$feature] >= results$lower_90 &
                              truth[results$feature] <= results$upper_90)
  metrics
}

if (identical(commandArgs(trailingOnly = TRUE), "--self-check")) {
  source("dev/continuous/simulate-fixed-truth.R")
  reference <- simulate_continuous_fixed_truth(10102026L, n_features = 20L)
  result <- fit_continuous_comparator(reference$data, "fixed")
  y <- result$inputs$y
  treated <- result$inputs$metadata$condition == "treated"
  difference <- rowMeans(y[, treated]) - rowMeans(y[, !treated])
  stopifnot(max(abs(result$results$estimate - difference)) < 1e-12,
            all(result$results$valid), all(result$results$lower_95 < result$results$lower_90),
            all(result$results$upper_95 > result$results$upper_90))
  reordered <- fit_continuous_comparator(reference$data[nrow(reference$data):1L, ], "fixed")
  stopifnot(isTRUE(all.equal(result$results,
    reordered$results[match(result$results$feature, reordered$results$feature), ],
    check.attributes = FALSE, tolerance = 1e-12)))
  reference <- simulate_continuous_fixed_truth(10102026L, design = "intercept")
  result <- fit_continuous_comparator(reference$data, "intercept")
  stopifnot(all(result$results$valid), nrow(result$results) == 20L,
            all(result$results$lower_95 < result$results$upper_95),
            identical(result$inputs$y[, 1L],
                      stats::setNames(reference$data$value[reference$data$sample == "sample_01"],
                                      result$inputs$features)))
  # Frozen v2 slope dataset: dream dropped two boundary-fit features before this was accepted.
  reference <- simulate_continuous_fixed_truth(13102029L, "gaussian", "slope", 20L, 20L, FALSE)
  result <- fit_continuous_comparator(reference$data, "slope")
  stopifnot(all(result$results$valid), nrow(result$results) == 20L, any(result$results$boundary))
  source("dev/continuous/discovery-metrics.R")
  stopifnot(is.finite(comparator_dataset_metrics(result$results, reference, 20L)$coverage_90))
  # Robust/trend variants keep the estimates and change only the moderated variance.
  base <- fit_continuous_comparator(reference$data, "slope")$results
  variant <- fit_continuous_comparator(reference$data, "slope", trend = TRUE, robust = TRUE)$results
  stopifnot(all(variant$valid), identical(unique(variant$method), "dream_trend_robust"),
            isTRUE(all.equal(base$estimate, variant$estimate)),
            !isTRUE(all.equal(base$lower_95, variant$lower_95)))
  for (design in c("fixed", "intercept", "slope")) {
    reference <- simulate_continuous_fixed_truth(15102027L, "variance_alternative", design)
    comparison <- fit_continuous_comparator(reference$data, design)
    spread <- fit_continuous_variance_comparator(comparison)
    stopifnot(all(spread$valid), identical(spread$feature, comparison$inputs$features),
              all(spread$component == "spread_association"),
              identical(spread$selected, spread$adjusted_p_value <= .05))
  }
  reference <- simulate_continuous_fixed_truth(15102027L, batch_layout = "partial")
  comparison <- fit_continuous_comparator(reference$data, "fixed", formula_location = ~ condition + batch)
  x <- stats::model.matrix(~ condition + batch, comparison$inputs$metadata)
  direct <- stats::lm.fit(x, t(comparison$inputs$y))$coefficients["conditiontreated", ]
  stopifnot(max(abs(direct - comparison$results$estimate)) < 1e-12)
  cat("Matched limma/dream inputs, OLS contrasts, intervals and reordering checks passed.\n")
}
