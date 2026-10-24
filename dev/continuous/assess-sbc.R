# Across independent datasets, never across dependent features within a dataset.
assess_sbc <- function(manifest, status, metrics, alpha = .05, strata = 91L) {
  stopifnot(!anyDuplicated(manifest$id), identical(manifest$id, status$id))
  stopifnot(is.character(status$error), !anyNA(status$error))
  all_hmc <- all(status$status == "hmc_passed") && all(!nzchar(status$error))
  if (is.null(metrics) || !nrow(metrics)) {
    return(list(all_hmc = all_hmc, coverage = NULL, ranks = NULL,
                sbc_campaign_passed = FALSE, scientific_qualified = FALSE))
  }
  stopifnot(!anyDuplicated(metrics[c("id", "feature", "parameter")]))
  metrics <- metrics[metrics$feature %in% c("global", "feature_01", "feature_10") &
                       metrics$id %in% status$id[status$status == "hmc_passed"], ]
  keys <- interaction(metrics$family, metrics$design, metrics$feature,
                       metrics$parameter, drop = TRUE)
  coverage <- list()
  ranks <- list()
  for (indices in split(seq_len(nrow(metrics)), keys)) {
    x <- metrics[indices, ]
    n <- sum(manifest$family == x$family[1L] & manifest$design == x$design[1L])
    identifying <- x[1L, c("family", "design", "feature", "parameter")]
    required <- if (x$design[1L] == "fixed") 200L else 100L
    for (level in c(.90, .95)) {
      count <- sum(x[[paste0("covered_", round(100 * level))]])
      acceptance <- stats::qbinom(c(alpha / (4 * strata), 1 - alpha / (4 * strata)), n, level)
      observed_ci <- if (nrow(x)) stats::binom.test(count, nrow(x))$conf.int else c(NA, NA)
      coverage[[length(coverage) + 1L]] <- cbind(identifying, data.frame(
        nominal = level, attempted = n, valid_intervals = nrow(x), covered = count,
        attempted_coverage = count / n,
        conditional_coverage = if (nrow(x)) count / nrow(x) else NA_real_,
        mc_lower = observed_ci[1L], mc_upper = observed_ci[2L],
        acceptance_lower = acceptance[1L] / n, acceptance_upper = acceptance[2L] / n,
        passed = n >= required && nrow(x) == n && count >= acceptance[1L] && count <= acceptance[2L]
      ))
    }
    stopifnot(length(unique(x$rank_draws)) == 1L, all(x$rank >= 0),
              all(x$rank <= x$rank_draws))
    grid <- seq.int(0L, x$rank_draws[1L])
    expected <- (grid + 1L) / (x$rank_draws[1L] + 1L)
    ecdf <- vapply(grid, function(k) mean(x$rank <= k), numeric(1))
    epsilon <- sqrt(log(2 * strata / alpha) / (2 * n))
    ranks[[length(ranks) + 1L]] <- cbind(identifying, data.frame(
      attempted = n, valid_ranks = nrow(x), grid = grid, expected = expected,
      observed = ecdf, lower = pmax(0, expected - epsilon), upper = pmin(1, expected + epsilon),
      max_deviation = max(abs(ecdf - expected)), envelope = epsilon,
      passed = n >= required && nrow(x) == n && max(abs(ecdf - expected)) <= epsilon
    ))
  }
  coverage <- if (length(coverage)) do.call(rbind, coverage) else NULL
  ranks <- if (length(ranks)) do.call(rbind, ranks) else NULL
  cells <- unique(manifest[c("family", "design")])
  complete_designs <- nrow(cells) == 6L && all(c("gaussian", "student_t") %in% cells$family) &&
    all(c("fixed", "intercept", "slope") %in% cells$design)
  sbc_pass <- all_hmc && complete_designs && !is.null(coverage) &&
    nrow(unique(coverage[c("family", "design", "feature", "parameter")])) == strata &&
    all(coverage$passed) && all(ranks$passed) && all(manifest$stage == "sbc")
  thinning <- if (all(c("rank_ess_bulk", "max_abs_rank_acf") %in% names(metrics))) {
    metrics[c("id", "family", "design", "feature", "parameter", "rank_draws",
              "rank_ess_bulk", "max_abs_rank_acf")]
  } else NULL
  list(all_hmc = all_hmc, coverage = coverage, ranks = ranks,
       rank_thinning_diagnostics = thinning,
       sbc_campaign_passed = sbc_pass, scientific_qualified = FALSE,
       rank_inspection = "Required separately; no automated human-inspection claim.",
       uncertainty = "Bonferroni binomial coverage bounds and DKW ECDF envelopes over 91 preselected strata; broad bounds are not precision evidence.")
}

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--self-check")) {
  manifest <- data.frame(id = sprintf("sbc-gaussian-fixed-%03d", 1:200),
                         family = "gaussian", design = "fixed", stage = "sbc")
  status <- data.frame(id = manifest$id, status = "hmc_passed", error = "")
  status_file <- tempfile(fileext = ".csv")
  utils::write.csv(status, status_file, row.names = FALSE)
  restored <- utils::read.csv(status_file, na.strings = character(),
                              colClasses = c(error = "character"))
  unlink(status_file)
  stopifnot(identical(status, restored))
  metrics <- data.frame(id = manifest$id, family = "gaussian", design = "fixed",
    feature = "feature_01", parameter = "beta[1]", covered_90 = seq_len(200) <= 180,
    covered_95 = seq_len(200) <= 190, rank = rep(0:100, length.out = 200), rank_draws = 100)
  result <- assess_sbc(manifest, status, metrics)
  stopifnot(all(result$coverage$passed), all(result$ranks$passed),
            !result$sbc_campaign_passed, !result$scientific_qualified)
  bad <- metrics
  bad$covered_95 <- FALSE
  bad$rank <- 0L
  result <- assess_sbc(manifest, status, bad)
  stopifnot(!all(result$coverage$passed), !all(result$ranks$passed))
  status$status[1L] <- "fit_failed"
  result <- assess_sbc(manifest, status, metrics)
  stopifnot(!result$all_hmc, all(result$coverage$attempted == 200),
            all(result$coverage$valid_intervals == 199), !any(result$coverage$passed))
  cat("Independent-dataset coverage/rank assessment self-checks passed.\n")
} else if (length(args)) {
  stopifnot(length(args) == 2L)
  campaign <- normalizePath(args[1L])
  output <- args[2L]
  if (file.exists(output)) stop("Assessment receipt exists; never overwrite.")
  manifest <- utils::read.csv(file.path(campaign, "manifest.csv"), stringsAsFactors = FALSE)
  status <- utils::read.csv(file.path(campaign, "fit-status.csv"), stringsAsFactors = FALSE,
                            na.strings = character(), colClasses = c(error = "character"))
  metrics_file <- file.path(campaign, "sbc-metrics.csv")
  metrics <- if (file.exists(metrics_file)) utils::read.csv(metrics_file, stringsAsFactors = FALSE) else NULL
  result <- assess_sbc(manifest, status, metrics)
  saveRDS(result, output)
  if (!is.null(result$coverage)) utils::write.csv(result$coverage,
    sub("\\.rds$", "-coverage.csv", output), row.names = FALSE)
  if (!is.null(result$ranks)) utils::write.csv(result$ranks,
    sub("\\.rds$", "-rank-ecdf.csv", output), row.names = FALSE)
  if (!is.null(result$rank_thinning_diagnostics)) utils::write.csv(
    result$rank_thinning_diagnostics, sub("\\.rds$", "-rank-thinning.csv", output),
    row.names = FALSE)
  cat("SBC campaign automated gates:", result$sbc_campaign_passed,
      "| Overall scientific qualification:", result$scientific_qualified, "\n")
}
