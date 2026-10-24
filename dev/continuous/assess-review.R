# Assess the review panels without treating expected rank rejection as sampler failure.
assess_review <- function(report, bounds, repetitions = 2000L) {
  source("dev/continuous/assess-fixed-truth.R", local = TRUE)
  expected <- report$manifest$expected_rejection
  stopifnot(is.logical(expected), !anyNA(expected),
            identical(report$manifest$id, report$status$id))
  expected_status <- report$status[expected, ]
  rejection_passed <- nrow(expected_status) > 0L &&
    all(expected_status$status == "design_rejected" & !nzchar(expected_status$error))
  ids <- report$manifest$dataset_id[!expected]
  comparison <- report$comparator_status[!report$comparator_status$id %in% ids, ]
  rejection_passed <- rejection_passed && nrow(comparison) > 0L &&
    all(comparison$status == "design_rejected" & !nzchar(comparison$error))
  regular <- report
  regular$manifest <- report$manifest[!expected, ]
  regular$status <- report$status[!expected, ]
  regular$comparator_status <- report$comparator_status[report$comparator_status$id %in% ids, ]
  assessment <- assess_fixed_truth(regular, bounds, repetitions = repetitions)
  assessment$expected_rejections <- expected_status
  assessment$expected_rejections_passed <- rejection_passed
  assessment$comparator_variants <- report$comparator_variants
  assessment$nuisance_metrics <- report$nuisance_metrics
  assessment$scientific_qualified <- FALSE
  assessment$reason <- paste("Review panels supplement v1; full SBC, software, performance and application gates remain mandatory.",
    "Omitted/precorrected nuisance fits are diagnostic misspecification controls, not promotion candidates.")
  assessment
}

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--self-check") && sys.nframe() == 0L) {
  source("dev/continuous/fixed-truth-campaign.R")
  source("dev/continuous/review-study.R")
  manifest <- review_manifest("review_pilot")
  expected <- manifest$expected_rejection
  datasets <- unique(manifest[c("dataset_id", "expected_rejection")])
  report <- list(manifest = manifest,
    status = data.frame(id = manifest$id,
      status = ifelse(expected, "design_rejected", "diagnostics_failed"), error = ""),
    comparator_status = data.frame(id = datasets$dataset_id,
      status = ifelse(datasets$expected_rejection, "design_rejected", "comparator_passed"), error = ""),
    metrics = NULL, features = NULL, comparator_metrics = NULL)
  result <- assess_review(report, fixed_truth_bounds(manifest), repetitions = 20L)
  stopifnot(result$expected_rejections_passed, !result$all_hmc, !result$scientific_qualified,
            sum(result$operating_characteristics$attempted) == 2L * sum(!expected),
            !any(result$coverage$passed))
  report$status$status[which(expected)[1L]] <- "hmc_passed"
  result <- assess_review(report, fixed_truth_bounds(manifest), repetitions = 20L)
  stopifnot(!result$expected_rejections_passed, !result$scientific_qualified)
  report$status$status <- ifelse(expected, "design_rejected", "hmc_passed")
  report$metrics <- do.call(rbind, lapply(which(!expected), function(i) {
    cbind(manifest[rep(i, 2L), ], data.frame(component = c("location", "log_variance"),
      fdp = 0, false_sign_proportion = 0, directional_fdp = 0,
      coverage_90 = .9, coverage_95 = .95, mean_width_95 = 2,
      discoveries = 0L, no_discoveries = TRUE, any_discovery = FALSE,
      mc_unstable_selections = 0L, max_probability_mcse = .001, rmse = .2))
  }))
  report$comparator_metrics <- data.frame(dataset_id = datasets$dataset_id, rmse = .2)
  result <- assess_review(report, fixed_truth_bounds(manifest), repetitions = 20L)
  stopifnot(result$all_hmc, result$expected_rejections_passed,
    all(result$student_width$width_ratio == 1), !any(result$student_width$passed),
    all(c("joint", "omit", "precorrect") %in% result$student_width$adjustment),
    !result$scientific_qualified)
  cat("Expected rank rejection, failure denominators and closed review qualification gates passed.\n")
} else if (length(args) && sys.nframe() == 0L) {
  stopifnot(length(args) == 2L)
  campaign <- normalizePath(args[1L])
  output <- args[2L]
  if (file.exists(output)) stop("Assessment exists; never overwrite.")
  report <- readRDS(file.path(campaign, "report.rds"))
  bounds <- utils::read.csv(file.path(campaign, "precomputed-fixed-truth-bounds.csv"),
                            colClasses = c(student_t_df = "character"))
  result <- assess_review(report, bounds)
  saveRDS(result, output)
  for (name in c("operating_characteristics", "coverage", "gaussian_accuracy", "student_width",
                 "comparator_variants", "nuisance_metrics", "expected_rejections")) {
    if (!is.null(result[[name]])) utils::write.csv(result[[name]],
      sub("\\.rds$", paste0("-", name, ".csv"), output), row.names = FALSE)
  }
  cat("Expected rejections:", result$expected_rejections_passed,
      "| All HMC:", result$all_hmc, "| Scientific qualification:", result$scientific_qualified, "\n")
}
