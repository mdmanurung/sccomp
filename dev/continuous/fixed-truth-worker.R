args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L, args[1L] %in% c("simulate", "fit", "comparator", "summarize"))
action <- args[1L]
campaign <- normalizePath(args[2L])
id <- args[3L]
output <- args[4L]
paths <- readLines(file.path(campaign, "paths.txt"))
source_root <- paths[2L]
.libPaths(c(file.path(campaign, "library"), paths[4L], .libPaths()))
library(sccomp)
hashes <- readRDS(file.path(campaign, "source-md5.rds"))
stopifnot(identical(unname(tools::md5sum(file.path(source_root, names(hashes)))), unname(hashes)))
for (file in c("simulate-fixed-truth.R", "comparators.R", "discovery-metrics.R", "fixed-truth-campaign.R")) {
  source(file.path(source_root, "dev/continuous", file))
}
source(file.path(source_root, "dev/continuous/review-study.R"))
manifest <- utils::read.csv(file.path(campaign, "manifest.csv"), stringsAsFactors = FALSE,
                           colClasses = c(student_t_df = "character"))
stopifnot(all(manifest$stage %in% c("fixed_pilot", "fixed_truth", "review_pilot", "review_truth")))
review <- all(startsWith(manifest$stage, "review_"))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
if (file.exists(output)) stop("Frozen outputs may not be overwritten.")
dataset_columns <- c("dataset_id", "scenario", "design", "replication", "seed",
                     "n_features", "n_donors", "observations_per_donor", "variance_effect")
if (review) dataset_columns <- c(dataset_columns, "noise_family", "hypervariable_fraction", "batch_layout")
if (action %in% c("simulate", "comparator")) {
  row <- unique(manifest[manifest$dataset_id == id, dataset_columns])
  stopifnot(nrow(row) == 1L)
} else if (action == "fit") {
  row <- manifest[manifest$id == id, ]
  stopifnot(nrow(row) == 1L)
}
if (action != "summarize") rownames(row) <- NULL
if (action == "simulate") {
  reference <- if (review) review_simulate(row) else {
    simulate_continuous_fixed_truth(row$seed, row$scenario, row$design,
      row$n_features, row$n_donors, row$variance_effect)
  }
  stopifnot(length(unique(reference$data$sample)) == row$n_donors * row$observations_per_donor)
  saveRDS(list(dataset = row, reference = reference), output)
  quit(status = 0L)
}
if (action %in% c("fit", "comparator")) {
  dataset_id <- if (action == "fit") row$dataset_id else id
  input_file <- file.path(campaign, "inputs", paste0(dataset_id, ".rds"))
  input <- readRDS(input_file)
  stopifnot(identical(row[dataset_columns], input$dataset))
  reference <- input$reference
  started <- Sys.time()
  receipt <- list(manifest = row, input_md5 = unname(tools::md5sum(input_file)),
    source_md5 = hashes, started = started, runtime = sessionInfo(),
    status = "fit_failed")
  attempted <- tryCatch({
    expected_rejection <- review && row$batch_layout == "confounded"
    if (expected_rejection) {
      error <- tryCatch({ review_validate_design(review_analysis(reference, "joint")); NULL },
                        error = conditionMessage)
      if (is.null(error) || !grepl("full rank", error)) stop("Expected confounded-design rank rejection was not observed.")
      receipt$rejection <- error
      receipt$status <- "design_rejected"
    } else if (action == "comparator") {
      adjustments <- if (review) unique(manifest$adjustment[manifest$dataset_id == id]) else "joint"
      variants <- list(plain = list(), robust = list(robust = TRUE),
                       trend_robust = list(trend = TRUE, robust = TRUE))
      comparisons <- list()
      for (adjustment in adjustments) {
        analysis <- review_analysis(reference, adjustment)
        comparison_base <- NULL
        for (variant in names(variants)) {
          key <- paste(adjustment, variant, sep = "-")
          started_variant <- Sys.time()
          comparisons[[key]] <- tryCatch({
            comparison <- do.call(fit_continuous_comparator,
              c(list(data = analysis$data, design = row$design,
                     formula_location = analysis$formula_location), variants[[variant]]))
            fit_file <- file.path(dirname(output), paste0(id,
              if (adjustment == "joint" && variant == "plain") "" else paste0("-", key), "-fit.rds"))
            if (file.exists(fit_file)) stop("Comparator fit exists; never overwrite.")
            saveRDS(comparison, fit_file)
            if (variant == "plain") comparison_base <- comparison
            metrics <- comparator_dataset_metrics(comparison$results, reference, row$n_features)
            metrics$method <- unique(comparison$results$method)
            metrics$adjustment <- adjustment
            list(status = "comparator_passed", metrics = metrics, features = comparison$results,
                 settings = comparison$settings, fit_file = fit_file)
          }, error = function(e) list(status = "fit_failed", error = conditionMessage(e)))
          comparisons[[key]]$elapsed_seconds <- as.numeric(difftime(Sys.time(), started_variant, units = "secs"))
        }
        key <- paste(adjustment, "diffVar", sep = "-")
        comparisons[[key]] <- tryCatch({
          if (is.null(comparison_base)) stop("Plain comparison failed; no spread reference is available.")
          spread <- fit_continuous_variance_comparator(comparison_base)
          stopifnot(all(spread$valid))
          list(status = "comparator_passed", features = spread,
               adjustment = adjustment,
               decision_rule = "BH spread association at zero; not the minimum variance-ratio rule.")
        }, error = function(e) list(status = "fit_failed", error = conditionMessage(e)))
      }
      receipt$comparisons <- comparisons
      primary <- comparisons[[paste(adjustments[1L], "plain", sep = "-")]]
      receipt$features <- primary$features
      receipt$metrics <- primary$metrics
      receipt$comparator_variants <- do.call(rbind, lapply(comparisons, `[[`, "metrics"))
      receipt$status <- if (all(vapply(comparisons, function(x) x$status == "comparator_passed", logical(1)))) {
        "comparator_passed"
      } else "comparator_failed"
      if (receipt$status == "comparator_failed") receipt$error <- paste(
        vapply(comparisons[vapply(comparisons, function(x) x$status == "fit_failed", logical(1))],
               `[[`, character(1), "error"), collapse = "; ")
    } else {
      analysis <- review_analysis(reference, if (review) row$adjustment else "joint")
      directory <- file.path(campaign, "draws", id)
      df <- if (row$student_t_df == "estimate") "estimate" else as.numeric(row$student_t_df)
      fit <- sccomp_estimate_continuous(analysis$data, analysis$formula_location, analysis$formula_variance,
        "sample", "feature", "value", family = row$family,
        variance_trend = row$variance_trend, student_t_df = df, scaling = row$scaling,
        cores = 4L, chains = row$chains, iter_warmup = row$warmup, iter_sampling = row$sampling,
        mcmc_seed = row$fit_seed, output_directory = directory,
        cache_stan_model = file.path(campaign, "compiled"), adapt_delta = row$adapt_delta,
        max_treedepth = row$max_treedepth, refresh = 0L, sig_figs = 16)
      saveRDS(fit, file.path(directory, "fit.rds"))
      receipt$diagnostics <- attr(fit, "continuous")$diagnostics
      receipt$model_load_elapsed_seconds <- attr(fit, "continuous")$model_load_elapsed_seconds
      if (review && row$adjustment == "joint") receipt$nuisance_metrics <- review_nuisance_metrics(fit, reference)
      receipt$status <- if (isTRUE(receipt$diagnostics$hmc_qualified)) "hmc_passed" else "diagnostics_failed"
      if (receipt$status == "hmc_passed") {
        tested <- sccomp_test(fit,
          contrasts = list(location = c(`treated-control` = "conditiontreated"),
                           log_variance = c(`treated/control` = "-2 * condition1")),
          effect_threshold = row$effect_threshold, variance_ratio_threshold = row$variance_ratio_threshold,
          posterior_fdr = row$posterior_fdr)
        saveRDS(tested, file.path(directory, "tested.rds"))
        receipt$features <- fixed_truth_feature_metrics(fit, tested, reference)
        receipt$metrics <- fixed_truth_dataset_metrics(receipt$features, row$posterior_fdr)
      }
    }
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(attempted)) {
    receipt$error <- attempted
    if (!is.null(receipt$diagnostics)) receipt$status <- "analysis_failed"
  }
  receipt$finished <- Sys.time()
  receipt$elapsed_seconds <- as.numeric(difftime(receipt$finished, started, units = "secs"))
  saveRDS(receipt, output)
  cat(id, receipt$status, "elapsed", receipt$elapsed_seconds, "seconds\n")
  quit(status = 0L)
}
receipts <- lapply(manifest$id, function(id) readRDS(file.path(campaign, "receipts", paste0(id, ".rds"))))
comparators <- lapply(unique(manifest$dataset_id), function(id) {
  readRDS(file.path(campaign, "comparators", paste0(id, ".rds")))
})
status <- function(receipts, ids) data.frame(id = ids,
  status = vapply(receipts, `[[`, character(1), "status"),
  elapsed_seconds = vapply(receipts, `[[`, numeric(1), "elapsed_seconds"),
  error = vapply(receipts, function(x) if (is.null(x$error)) "" else x$error, character(1)))
fit_status <- status(receipts, manifest$id)
comparator_status <- status(comparators, unique(manifest$dataset_id))
utils::write.csv(fit_status, file.path(campaign, "fit-status.csv"), row.names = FALSE)
utils::write.csv(comparator_status, file.path(campaign, "comparator-status.csv"), row.names = FALSE)
bind_metrics <- function(receipts, field) {
  do.call(rbind, lapply(receipts, function(x) {
    if (is.null(x[[field]])) return(NULL)
    cbind(x$manifest, x[[field]])
  }))
}
metrics <- bind_metrics(receipts, "metrics")
features <- bind_metrics(receipts, "features")
comparator_metrics <- bind_metrics(comparators, "metrics")
comparator_variants <- bind_metrics(comparators, "comparator_variants")
nuisance_metrics <- bind_metrics(receipts, "nuisance_metrics")
for (name in c("metrics", "features", "comparator_metrics")) {
  table <- get(name)
  if (!is.null(table)) utils::write.csv(table, file.path(campaign, paste0(name, ".csv")), row.names = FALSE)
}
saveRDS(list(stage = unique(manifest$stage), manifest = manifest,
  status = fit_status, comparator_status = comparator_status, metrics = metrics,
  features = features, comparator_metrics = comparator_metrics, comparator_variants = comparator_variants,
  nuisance_metrics = nuisance_metrics,
  attempted = nrow(manifest), hmc_passed = sum(fit_status$status == "hmc_passed"),
  scientific_qualified = FALSE,
  reason = "Execution is separate from dataset-level uncertainty, calibration and performance qualification."), output)
