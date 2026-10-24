# Freeze shared dataset identities before any fit; models differ only by settings.
fixed_truth_manifest <- function(stage) {
  stopifnot(stage %in% c("fixed_pilot", "fixed_truth"))
  pilot <- stage == "fixed_pilot"
  scenarios <- if (pilot) "gaussian" else c(
    "null", "gaussian", "student_t", "contamination", "correlated",
    "hypervariable", "unbalanced", "few_donors", "non_ig", "trend",
    "variance_null", "variance_alternative")
  rows <- list()
  dataset_index <- 0L
  for (scenario in scenarios) for (design in c("fixed", "intercept", "slope")) {
    for (replication in seq_len(if (pilot) 1L else 100L)) {
      dataset_index <- dataset_index + 1L
      dataset <- sprintf("%s-%s-%s-%04d", stage, scenario, design, replication)
      settings <- data.frame(family = c("gaussian", "student_t"),
        variance_trend = "none", student_t_df = "estimate", stringsAsFactors = FALSE)
      if (!pilot && scenario == "trend") {
        extra <- settings
        extra$variance_trend <- "linear"
        settings <- rbind(settings, extra)
      }
      if (!pilot && scenario == "gaussian") settings <- rbind(settings,
        data.frame(family = "student_t", variance_trend = "none", student_t_df = "1000000"))
      for (j in seq_len(nrow(settings))) {
        index <- length(rows) + 1L
        rows[[index]] <- data.frame(
          id = paste(dataset, settings$family[j], settings$variance_trend[j],
                     settings$student_t_df[j], sep = "-"),
          dataset_id = dataset, stage = stage, scenario = scenario, design = design,
          replication = replication, family = settings$family[j],
          variance_trend = settings$variance_trend[j], student_t_df = settings$student_t_df[j],
          seed = 10102026L + dataset_index + if (pilot) 3000000L else 4000000L,
          fit_seed = 20102026L + index + if (pilot) 3000000L else 4000000L,
          n_features = 20L, n_donors = if (scenario == "few_donors") 5L else 20L,
          observations_per_donor = if (design == "slope") 4L else 2L,
          variance_effect = scenario == "variance_alternative", scaling = "global",
          chains = 4L, warmup = 1500L, sampling = 1500L,
          adapt_delta = .995, max_treedepth = 12L,
          effect_threshold = 1, variance_ratio_threshold = 1.5,
          posterior_fdr = .05, protocol = "v2-identifiable-slopes", stringsAsFactors = FALSE)
      }
    }
  }
  manifest <- do.call(rbind, rows)
  stopifnot(!anyDuplicated(manifest$id), !anyDuplicated(manifest$fit_seed))
  manifest
}

fixed_truth_bounds <- function(manifest) {
  columns <- c("scenario", "design", "family", "variance_trend", "student_t_df")
  if ("adjustment" %in% names(manifest)) columns <- c(columns, "adjustment")
  cells <- unique(manifest[columns])
  cells$datasets <- vapply(seq_len(nrow(cells)), function(i) {
    sum(Reduce(`&`, lapply(columns, function(name) manifest[[name]] == cells[[name]][i])))
  }, integer(1))
  # Two features, two components, two coverage levels, with simultaneous tails.
  strata <- nrow(cells) * 2L * 2L
  bounds <- cells[rep(seq_len(nrow(cells)), each = 2L), ]
  bounds$nominal <- rep(c(.90, .95), nrow(cells))
  bounds$coverage_lower_count <- stats::qbinom(.05 / (4 * strata), bounds$datasets, bounds$nominal)
  bounds$coverage_upper_count <- stats::qbinom(1 - .05 / (4 * strata), bounds$datasets, bounds$nominal)
  bounds$diagnostic_strata <- strata
  bounds$null_any_discovery_upper_count <- stats::qbinom(
    1 - .05 / (2 * nrow(cells)), bounds$datasets, .05)
  bounds
}

fixed_truth_feature_metrics <- function(fit, tested, reference) {
  input <- attr(fit, "continuous")$input
  features <- input$features
  generated_features <- unique(as.character(reference$data$feature))
  stopifnot(identical(features, generated_features),
            setequal(unique(tested$component), c("location", "log_variance")),
            nrow(tested) == 2L * length(features))
  draws <- unclass(attr(fit, "fit"))
  result <- lapply(c("location", "log_variance"), function(component) {
    x <- as.data.frame(lapply(tested[tested$component == component, ], identity))
    stopifnot(identical(as.character(x$feature), features))
    truth <- if (component == "location") reference$truth$beta[2L, ] else {
      -2 * reference$truth$gamma[1L, ]
    }
    coefficient <- if (component == "location") {
      match("conditiontreated", colnames(input$stan_data$X))
    } else match("condition1", colnames(input$stan_data$W))
    prefix <- if (component == "location") "beta" else "gamma"
    stopifnot(!is.na(coefficient))
    intervals <- vapply(seq_along(features), function(g) {
      values <- draws[, , sprintf("%s[%d,%d]", prefix, coefficient, g)]
      values <- if (component == "location") values * input$global_scale else -2 * values
      stopifnot(all(is.finite(values)), abs(mean(values) - x$estimate[g]) < 1e-8)
      as.numeric(stats::quantile(values, c(.05, .95)))
    }, numeric(2L))
    x$truth <- truth
    x$lower_90 <- intervals[1L, ]
    x$upper_90 <- intervals[2L, ]
    x$covered_90 <- truth >= x$lower_90 & truth <= x$upper_90
    x$covered_95 <- truth >= x$lower & truth <= x$upper
    x
  })
  do.call(rbind, result)
}

fixed_truth_dataset_metrics <- function(features, fdr = .05) {
  source_discovery <- continuous_discovery_metrics
  rows <- lapply(split(features, features$component), function(x) {
    result <- source_discovery(x, stats::setNames(x$truth, x$feature), unique(x$threshold))
    result$coverage_90 <- mean(x$covered_90)
    result$component <- unique(x$component)
    result$parameter <- unique(x$parameter)
    # Perturb directional error probabilities using their autocorrelation-aware MCSE.
    uncertainty <- 2 * x$probability_mcse
    uncertainty[!is.finite(uncertainty)] <- 1
    selection <- function(probability) {
      order <- order(probability)
      selected <- logical(length(probability))
      sorted <- probability[order]
      ends <- which(!duplicated(sorted, fromLast = TRUE))
      cumulative <- cumsum(sorted) / seq_along(order)
      selected[order] <- cumulative[ends][match(sorted, sorted[ends])] <= fdr
      selected
    }
    low <- selection(pmax(0, x$p_error - uncertainty))
    high <- selection(pmin(1, x$p_error + uncertainty))
    result$mc_unstable_selections <- sum(low != high)
    result$max_probability_mcse <- if (any(!is.finite(x$probability_mcse))) NA_real_ else {
      max(x$probability_mcse)
    }
    result
  })
  do.call(rbind, rows)
}

if (identical(commandArgs(trailingOnly = TRUE), "--self-check") && sys.nframe() == 0L) {
  pilot <- fixed_truth_manifest("fixed_pilot")
  full <- fixed_truth_manifest("fixed_truth")
  stopifnot(nrow(pilot) == 6L, nrow(full) == 8100L,
    length(unique(full$dataset_id)) == 3600L,
    !any(pilot$seed %in% full$seed), all(full$scaling == "global"),
    all(full$effect_threshold == 1), all(full$variance_ratio_threshold == 1.5))
  datasets <- split(full, full$dataset_id)
  stopifnot(all(vapply(datasets, function(x) length(unique(x$seed)) == 1L, logical(1))))
  bounds <- fixed_truth_bounds(full)
  stopifnot(all(bounds$datasets == 100L), all(bounds$coverage_lower_count <= bounds$coverage_upper_count),
    nrow(bounds) == 162L, all(bounds$diagnostic_strata == 324L))
  source("dev/continuous/discovery-metrics.R")
  x <- data.frame(feature = c("a", "b"), component = "location", parameter = "treated-control",
    estimate = c(0, 2), truth = c(0, 2), lower = c(-.5, 1.5), upper = c(.5, 2.5),
    threshold = 1, selected = c(FALSE, TRUE), direction = "positive",
    covered_90 = TRUE, covered_95 = TRUE, p_error = c(.99, .01), probability_mcse = .001)
  result <- fixed_truth_dataset_metrics(x)
  stopifnot(result$discoveries == 1L, result$fdp == 0, result$mc_unstable_selections == 0L)
  x$p_error[2L] <- .05
  x$probability_mcse[2L] <- .01
  stopifnot(fixed_truth_dataset_metrics(x)$mc_unstable_selections == 1L)
  x$probability_mcse[2L] <- NA_real_
  stopifnot(is.na(fixed_truth_dataset_metrics(x)$max_probability_mcse))
  cat("Fixed-truth shared seeds, full-panel manifest, frozen bounds and selection-MC checks passed.\n")
}
