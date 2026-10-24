# Review extensions reuse the independent fixed-truth simulator and frozen workers.
review_manifest <- function(stage) {
  stopifnot(stage %in% c("review_pilot", "review_truth"))
  pilot <- stage == "review_pilot"
  template <- fixed_truth_manifest("fixed_pilot")[1L, ]
  panels <- expand.grid(noise_family = c("gaussian", "student_t"),
    hypervariable_fraction = c(0, .05, .20), stringsAsFactors = FALSE)
  panels$batch_layout <- "none"
  batches <- data.frame(noise_family = "gaussian", hypervariable_fraction = 0,
                        batch_layout = c("balanced", "partial", "confounded"))
  panels <- rbind(panels, batches)
  rows <- list()
  index <- 0L
  for (p in seq_len(nrow(panels))) for (design in c("fixed", "intercept", "slope")) {
    panel <- panels[p, ]
    scenario <- if (panel$batch_layout == "none") {
      sprintf("review_variance_%s_%02d", panel$noise_family, 100 * panel$hypervariable_fraction)
    } else paste0("review_batch_", panel$batch_layout)
    adjustments <- if (panel$batch_layout %in% c("none", "confounded")) "joint" else {
      c("joint", "omit", "precorrect")
    }
    for (replication in seq_len(if (pilot) 1L else 100L)) {
      index <- index + 1L
      dataset <- sprintf("%s-%s-%s-%04d", stage, scenario, design, replication)
      for (adjustment in adjustments) for (family in c("gaussian", "student_t")) {
        for (trend in if (panel$batch_layout == "none") c("none", "linear") else "none") {
          row <- template
          row$id <- paste(dataset, adjustment, family, trend, sep = "-")
          row$dataset_id <- dataset
          row$stage <- stage
          row$scenario <- scenario
          row$design <- design
          row$replication <- replication
          row$family <- family
          row$variance_trend <- trend
          row$seed <- 10102026L + index + if (pilot) 5000000L else 6000000L
          row$fit_seed <- 20102026L + length(rows) + 1L + if (pilot) 5000000L else 6000000L
          row$observations_per_donor <- if (design == "slope") 4L else 2L
          row$adjustment <- adjustment
          row$noise_family <- panel$noise_family
          row$hypervariable_fraction <- panel$hypervariable_fraction
          row$batch_layout <- panel$batch_layout
          row$expected_rejection <- panel$batch_layout == "confounded"
          row$protocol <- "review-v1"
          rows[[length(rows) + 1L]] <- row
        }
      }
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  stopifnot(!anyDuplicated(result$id), !anyDuplicated(result$fit_seed))
  result
}

review_simulate <- function(row) {
  simulate_continuous_fixed_truth(row$seed, "gaussian", row$design,
    row$n_features, row$n_donors, row$variance_effect,
    hypervariable_fraction = row$hypervariable_fraction,
    noise_family = row$noise_family, batch_layout = row$batch_layout)
}

review_analysis <- function(reference, adjustment) {
  stopifnot(adjustment %in% c("joint", "omit", "precorrect"))
  data <- reference$data
  if (adjustment == "precorrect") {
    stopifnot("batch" %in% names(data))
    # Estimate batch alone, then treat its subtraction as fixed: a diagnostic comparator.
    for (indices in split(seq_len(nrow(data)), data$feature)) {
      x <- stats::model.matrix(~ batch, data[indices, ])
      coefficient <- stats::lm.fit(x, data$value[indices])$coefficients[2L]
      stopifnot(is.finite(coefficient))
      data$value[indices] <- data$value[indices] - x[, 2L] * coefficient
    }
  }
  location <- switch(reference$design, fixed = ~ condition,
    intercept = ~ condition + (1 | donor), slope = ~ condition + (1 + condition | donor))
  variance <- ~ condition
  if (adjustment == "joint" && "batch" %in% names(data)) {
    location <- update(location, ~ . + batch)
    variance <- ~ condition + batch
  }
  list(data = data, formula_location = location, formula_variance = variance)
}

review_validate_design <- function(analysis, family = "gaussian") {
  getFromNamespace("continuous_model_input", "sccomp")(
    analysis$data, analysis$formula_location, analysis$formula_variance,
    "sample", "feature", "value", family, "none", "global", "estimate",
    sccomp::sccomp_continuous_prior())
  invisible(TRUE)
}

review_nuisance_metrics <- function(fit, reference) {
  if (is.null(reference$truth$batch_beta)) return(NULL)
  rows <- list()
  for (component in c("location", "log_variance")) {
    parameter <- if (component == "location") "batchB" else "batch1"
    x <- fit[fit$component == component & fit$parameter == parameter, ]
    if (!nrow(x)) next
    truth <- if (component == "location") reference$truth$batch_beta else reference$truth$batch_gamma
    names(truth) <- unique(reference$data$feature)
    truth <- truth[x$feature]
    rows[[length(rows) + 1L]] <- data.frame(component = component, parameter = parameter,
      rmse = sqrt(mean((x$estimate - truth)^2)), bias = mean(x$estimate - truth),
      coverage_95 = mean(truth >= x$lower & truth <= x$upper))
  }
  do.call(rbind, rows)
}

if (identical(commandArgs(trailingOnly = TRUE), "--self-check") && sys.nframe() == 0L) {
  pkgload::load_all(".", quiet = TRUE, helpers = FALSE)
  source("dev/continuous/fixed-truth-campaign.R")
  source("dev/continuous/simulate-fixed-truth.R")
  pilot <- review_manifest("review_pilot")
  full <- review_manifest("review_truth")
  stopifnot(nrow(pilot) == 114L, nrow(full) == 11400L,
    length(unique(full$dataset_id)) == 2700L,
    !any(pilot$seed %in% full$seed), all(full$protocol == "review-v1"))
  for (design in c("fixed", "intercept", "slope")) {
    baseline <- simulate_continuous_fixed_truth(15102027L, design = design,
                                               hypervariable_fraction = 0)
    for (fraction in c(.05, .20)) {
      stressed <- simulate_continuous_fixed_truth(15102027L, design = design,
                                                  hypervariable_fraction = fraction)
      stopifnot(sum(stressed$truth$hypervariable) == 20L * fraction,
        identical(stressed$truth$beta, baseline$truth$beta),
        all(abs(stressed$truth$sigma_sq / baseline$truth$sigma_sq -
                  ifelse(stressed$truth$hypervariable, 10, 1)) < 1e-12))
    }
    for (layout in c("balanced", "partial", "confounded")) {
      reference <- simulate_continuous_fixed_truth(15102027L, design = design, batch_layout = layout)
      analysis <- review_analysis(reference, "joint")
      error <- tryCatch({ review_validate_design(analysis); NULL }, error = conditionMessage)
      if (layout == "confounded") stopifnot(!is.null(error), grepl("full rank", error)) else {
        stopifnot(is.null(error))
        corrected <- review_analysis(reference, "precorrect")
        stopifnot(!identical(corrected$data$value, reference$data$value),
                  identical(corrected$data$condition, reference$data$condition))
      }
    }
  }
  bounds <- fixed_truth_bounds(full)
  stopifnot(all(bounds$datasets == 100L))
  cat("Review manifest, shared inputs, tenfold variance stress and nuisance rank checks passed.\n")
}
