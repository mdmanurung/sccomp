# Resample independent datasets jointly for paired method comparisons.
fixed_truth_bootstrap <- function(values, statistic = mean, repetitions = 2000L) {
  stopifnot(nrow(as.matrix(values)) >= 1L, all(is.finite(values)))
  n <- nrow(as.matrix(values))
  observed <- statistic(values)
  sampled <- replicate(repetitions, {
    indices <- sample.int(n, n, replace = TRUE)
    statistic(if (is.matrix(values)) values[indices, , drop = FALSE] else values[indices])
  })
  c(estimate = observed, lower = unname(stats::quantile(sampled, .025)),
    upper = unname(stats::quantile(sampled, .975)))
}

assess_fixed_truth <- function(report, bounds, repetitions = 2000L, seed = 90102026L) {
  manifest <- report$manifest
  status <- report$status
  stopifnot(identical(manifest$id, status$id), !anyDuplicated(manifest$id),
            nrow(manifest) > 0L, !anyDuplicated(report$comparator_status$id))
  set.seed(seed)
  metrics <- report$metrics
  features <- report$features
  if (!is.null(metrics)) stopifnot(all(c("id", "dataset_id", "component",
    "max_probability_mcse", "mc_unstable_selections") %in% names(metrics)),
    !anyDuplicated(metrics[c("id", "component")]), all(metrics$id %in% manifest$id))
  passed <- status$id[status$status == "hmc_passed" & !nzchar(status$error)]
  if (!is.null(metrics)) metrics <- metrics[metrics$id %in% passed, ]
  if (!is.null(features)) features <- features[features$id %in% passed, ]
  comparison <- report$comparator_metrics
  comparator_ids <- report$comparator_status$id[
    report$comparator_status$status == "comparator_passed" & !nzchar(report$comparator_status$error)]
  if (!is.null(comparison)) comparison <- comparison[comparison$dataset_id %in% comparator_ids, ]
  cell_columns <- c("scenario", "design", "family", "variance_trend", "student_t_df")
  if ("adjustment" %in% names(manifest)) cell_columns <- c(cell_columns, "adjustment")
  cells <- unique(manifest[cell_columns])
  match_cell <- function(table, cell) {
    if (is.null(table) || !nrow(table)) return(logical())
    Reduce(`&`, lapply(cell_columns, function(name) table[[name]] == cell[[name]]))
  }
  operating <- list()
  coverage <- list()
  accuracy <- list()
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]
    planned <- manifest[match_cell(manifest, cell), ]
    bound <- bounds[match_cell(bounds, cell), ]
    stopifnot(nrow(bound) == 2L, all(bound$datasets == nrow(planned)))
    for (component in c("location", "log_variance")) {
      x <- if (is.null(metrics)) NULL else metrics[
        match_cell(metrics, cell) & metrics$component == component, ]
      n <- if (is.null(x)) 0L else nrow(x)
      if (n) stopifnot(!anyDuplicated(x$dataset_id), all(x$id %in% planned$id))
      summary <- cbind(cell, data.frame(component = component,
        attempted = nrow(planned), valid = n, failed = nrow(planned) - n,
        discoveries = if (n) sum(x$discoveries) else NA_integer_,
        no_discovery_datasets = if (n) sum(x$no_discoveries) else NA_integer_,
        mc_unstable_selections = if (n) sum(x$mc_unstable_selections) else NA_integer_))
      for (metric in c("fdp", "false_sign_proportion", "directional_fdp",
                       "coverage_90", "coverage_95", "mean_width_95")) {
        estimate <- if (n) fixed_truth_bootstrap(x[[metric]], repetitions = repetitions) else {
          c(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
        }
        for (field in names(estimate)) summary[[paste(metric, field, sep = "_")]] <- unname(estimate[field])
      }
      global_null <- if (component == "location") {
        cell$scenario %in% c("null", "variance_null", "variance_alternative")
      } else cell$scenario != "variance_alternative"
      summary$global_null <- global_null
      summary$selection_mc_qualified <- n == nrow(planned) && n > 0L &&
        all(is.finite(x$max_probability_mcse)) && all(x$mc_unstable_selections == 0L)
      summary$any_discovery_datasets <- if (n) sum(x$any_discovery) else NA_integer_
      null_ci <- if (n && global_null) stats::binom.test(sum(x$any_discovery), n)$conf.int else c(NA, NA)
      summary$null_any_discovery_lower <- null_ci[1L]
      summary$null_any_discovery_upper <- null_ci[2L]
      summary$null_compatibility_upper_count <- bound$null_any_discovery_upper_count[1L]
      summary$decision_compatible <- n == nrow(planned) && n >= 100L && if (global_null) {
        sum(x$any_discovery) <= summary$null_compatibility_upper_count
      } else summary$directional_fdp_lower <= .05
      # A wide interval can be compatible without being precise enough to qualify.
      summary$decision_precise <- n >= 100L && if (global_null) {
        summary$null_any_discovery_upper - summary$null_any_discovery_lower <= .05
      } else summary$directional_fdp_upper - summary$directional_fdp_lower <= .05
      operating[[length(operating) + 1L]] <- summary
      for (feature in c("feature_0001", "feature_0010")) {
        f <- if (is.null(features)) NULL else features[
          match_cell(features, cell) & features$component == component & features$feature == feature, ]
        if (!is.null(f)) stopifnot(!anyDuplicated(f$dataset_id))
        valid <- if (is.null(f)) 0L else nrow(f)
        for (level in c(.90, .95)) {
          b <- bound[bound$nominal == level, ]
          count <- if (valid) sum(f[[paste0("covered_", round(100 * level))]]) else 0L
          ci <- if (valid) stats::binom.test(count, valid)$conf.int else c(NA, NA)
          coverage[[length(coverage) + 1L]] <- cbind(cell, data.frame(component = component,
            feature = feature, nominal = level, attempted = nrow(planned), valid = valid,
            covered = count, attempted_coverage = count / nrow(planned),
            conditional_coverage = if (valid) count / valid else NA_real_,
            mc_lower = ci[1L], mc_upper = ci[2L],
            acceptance_lower_count = b$coverage_lower_count,
            acceptance_upper_count = b$coverage_upper_count,
            passed = valid == nrow(planned) && valid >= 100L &&
              count >= b$coverage_lower_count && count <= b$coverage_upper_count))
        }
      }
      if (component == "location" && cell$scenario %in% c("gaussian", "review_variance_gaussian_00") && cell$family == "gaussian" &&
          cell$variance_trend == "none" && n) {
        comparator <- if (is.null(comparison)) NULL else comparison[
          match(x$dataset_id, comparison$dataset_id), ]
        valid_comparator <- !is.null(comparator) && nrow(comparator) == n &&
          all(!is.na(comparator$dataset_id)) && identical(x$dataset_id, comparator$dataset_id)
        ratio <- if (valid_comparator) fixed_truth_bootstrap(cbind(x$rmse, comparator$rmse),
          statistic = function(z) sqrt(mean(z[, 1L]^2) / mean(z[, 2L]^2)),
          repetitions = repetitions) else c(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
        accuracy[[length(accuracy) + 1L]] <- cbind(cell, data.frame(attempted = nrow(planned),
          paired = if (valid_comparator) n else 0L, rmse_ratio = unname(ratio["estimate"]),
          mc_lower = unname(ratio["lower"]), mc_upper = unname(ratio["upper"]),
          passed = valid_comparator && n == nrow(planned) && n >= 100L && ratio["upper"] <= 1.05))
      }
    }
  }
  widths <- list()
  if (!is.null(metrics)) {
    student <- metrics[metrics$family == "student_t" & metrics$student_t_df == "estimate", ]
    gaussian <- metrics[metrics$family == "gaussian", ]
    keys <- c("dataset_id", "variance_trend", "component")
    if ("adjustment" %in% names(metrics)) keys <- c(keys, "adjustment")
    paired <- merge(student, gaussian, by = keys, suffixes = c("_student", "_gaussian"))
    if (nrow(paired)) for (indices in split(seq_len(nrow(paired)), interaction(
      paired$scenario_student, paired$design_student, paired$variance_trend, paired$component,
      if ("adjustment" %in% names(paired)) paired$adjustment else "joint", drop = TRUE))) {
      x <- paired[indices, ]
      stopifnot(!anyDuplicated(x$dataset_id))
      ratio <- fixed_truth_bootstrap(cbind(x$mean_width_95_student, x$mean_width_95_gaussian),
        statistic = function(z) mean(z[, 1L]) / mean(z[, 2L]), repetitions = repetitions)
      improvement <- fixed_truth_bootstrap(x$coverage_95_student - x$coverage_95_gaussian,
                                            repetitions = repetitions)
      attempted <- sum(manifest$scenario == x$scenario_student[1L] & manifest$design == x$design_student[1L] &
        manifest$family == "student_t" & manifest$variance_trend == x$variance_trend[1L] &
        manifest$student_t_df == "estimate")
      if ("adjustment" %in% names(manifest)) attempted <- sum(
        manifest$scenario == x$scenario_student[1L] & manifest$design == x$design_student[1L] &
        manifest$family == "student_t" & manifest$variance_trend == x$variance_trend[1L] &
        manifest$student_t_df == "estimate" & manifest$adjustment == x$adjustment[1L])
      widths[[length(widths) + 1L]] <- data.frame(scenario = x$scenario_student[1L], design = x$design_student[1L],
        variance_trend = x$variance_trend[1L], component = x$component[1L], attempted = attempted,
        paired = nrow(x), width_ratio = unname(ratio["estimate"]),
        width_mc_lower = unname(ratio["lower"]), width_mc_upper = unname(ratio["upper"]),
        coverage_improvement = unname(improvement["estimate"]),
        coverage_mc_lower = unname(improvement["lower"]), coverage_mc_upper = unname(improvement["upper"]),
        passed = nrow(x) == attempted && attempted >= 100L &&
          (ratio["upper"] <= 1.10 || improvement["lower"] > 0))
      if ("adjustment" %in% names(x)) widths[[length(widths)]]$adjustment <- x$adjustment[1L]
    }
  }
  list(all_hmc = length(passed) == nrow(manifest),
    all_comparators = length(comparator_ids) == length(unique(manifest$dataset_id)),
    operating_characteristics = do.call(rbind, operating), coverage = do.call(rbind, coverage),
    gaussian_accuracy = do.call(rbind, accuracy), student_width = do.call(rbind, widths),
    bootstrap_seed = seed, bootstrap_repetitions = repetitions,
    scientific_qualified = FALSE,
    reason = "Panel compatibility, Monte Carlo precision, SBC, performance and application review remain separate gates.")
}

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--self-check") && sys.nframe() == 0L) {
  source("dev/continuous/fixed-truth-campaign.R")
  manifest <- fixed_truth_manifest("fixed_pilot")
  report <- list(manifest = manifest,
    status = data.frame(id = manifest$id, status = "hmc_passed", error = ""),
    comparator_status = data.frame(id = unique(manifest$dataset_id), status = "comparator_passed", error = ""))
  report$metrics <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
    cbind(manifest[i, ], data.frame(component = c("location", "log_variance"),
      fdp = 0, false_sign_proportion = 0, directional_fdp = 0, coverage_90 = .9,
      coverage_95 = .95, mean_width_95 = 2, discoveries = 1L, no_discoveries = FALSE,
      any_discovery = TRUE, mc_unstable_selections = 0L, rmse = .2))
  }))
  report$metrics$max_probability_mcse <- .001
  report$comparator_metrics <- data.frame(dataset_id = unique(manifest$dataset_id), rmse = .2)
  report$features <- NULL
  result <- assess_fixed_truth(report, fixed_truth_bounds(manifest), repetitions = 20L)
  stopifnot(result$all_hmc, result$all_comparators, !result$scientific_qualified,
    !any(result$coverage$passed), all(result$gaussian_accuracy$rmse_ratio == 1),
    all(result$student_width$width_ratio == 1), !any(result$student_width$passed))
  report$status$status[1L] <- "analysis_failed"
  result <- assess_fixed_truth(report, fixed_truth_bounds(manifest), repetitions = 20L)
  stopifnot(!result$all_hmc, sum(result$operating_characteristics$failed) == 2L,
            !any(result$coverage$passed), !result$scientific_qualified)
  values <- cbind(c(1, 2, 3), c(1, 2, 3))
  paired <- fixed_truth_bootstrap(values, function(x) mean(x[, 1L] - x[, 2L]), repetitions = 20L)
  stopifnot(all(paired == 0))
  cat("Independent-dataset bootstrap, paired comparisons, failure denominators and closed pilot gates passed.\n")
} else if (length(args) && sys.nframe() == 0L) {
  stopifnot(length(args) == 2L)
  campaign <- normalizePath(args[1L])
  output <- args[2L]
  if (file.exists(output)) stop("Assessment exists; never overwrite.")
  report <- readRDS(file.path(campaign, "report.rds"))
  bounds <- utils::read.csv(file.path(campaign, "precomputed-fixed-truth-bounds.csv"),
                            colClasses = c(student_t_df = "character"))
  result <- assess_fixed_truth(report, bounds)
  saveRDS(result, output)
  for (name in c("operating_characteristics", "coverage", "gaussian_accuracy", "student_width")) {
    if (!is.null(result[[name]])) utils::write.csv(result[[name]],
      sub("\\.rds$", paste0("-", name, ".csv"), output), row.names = FALSE)
  }
  cat("All HMC:", result$all_hmc, "| Overall scientific qualification:", result$scientific_qualified, "\n")
}
