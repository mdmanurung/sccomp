args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L, args[1L] %in% c("simulate", "fit", "summarize"))
action <- args[1L]
campaign <- normalizePath(args[2L])
id <- args[3L]
output <- args[4L]
paths <- readLines(file.path(campaign, "paths.txt"))
source_root <- paths[2L]
.libPaths(c(file.path(campaign, "library"), paths[4L], .libPaths()))
library(sccomp)
source(file.path(source_root, "dev/continuous/simulate-reference.R"))
hashes <- readRDS(file.path(campaign, "source-md5.rds"))
stopifnot(identical(unname(tools::md5sum(file.path(source_root, names(hashes)))), unname(hashes)))
manifest <- utils::read.csv(file.path(campaign, "manifest.csv"), stringsAsFactors = FALSE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
if (file.exists(output)) stop("Frozen outputs may not be overwritten.")
if (action != "summarize") {
  row <- manifest[manifest$id == id, ]
  stopifnot(nrow(row) == 1L)
}
if (action == "simulate") {
  reference <- simulate_continuous_reference(row$seed, row$family, row$design,
    row$n_features, row$n_donors, row$variance_effect)
  saveRDS(list(manifest = row, reference = reference), output)
  quit(status = 0L)
}
if (action == "fit") {
  input_file <- file.path(campaign, "inputs", paste0(id, ".rds"))
  input <- readRDS(input_file)
  stopifnot(identical(row, input$manifest))
  reference <- input$reference
  formula <- switch(row$design, fixed = ~ condition,
                    intercept = ~ condition + (1 | donor),
                    slope = ~ condition + (1 + condition | donor))
  started <- Sys.time()
  fit_directory <- file.path(campaign, "draws", id)
  receipt <- list(manifest = row, input_md5 = tools::md5sum(input_file),
                   started = started, source_md5 = hashes,
                   runtime = sessionInfo(), status = "fit_failed")
  attempted <- tryCatch({
    fit <- sccomp_estimate_continuous(
      reference$data, formula, ~ condition, "sample", "feature", "value",
      family = row$family, scaling = "none", prior = sccomp_continuous_prior(),
      cores = 4L, chains = row$chains, iter_warmup = row$warmup,
      iter_sampling = row$sampling, mcmc_seed = row$fit_seed,
      output_directory = fit_directory,
      cache_stan_model = file.path(campaign, "compiled"),
      adapt_delta = row$adapt_delta, max_treedepth = row$max_treedepth,
      refresh = 0L, sig_figs = 16
    )
    saveRDS(fit, file.path(fit_directory, "fit.rds"))
    receipt$diagnostics <- attr(fit, "continuous")$diagnostics
    receipt$status <- if (isTRUE(receipt$diagnostics$hmc_qualified)) "hmc_passed" else "diagnostics_failed"
    variables <- character()
    truth <- numeric()
    features <- character()
    add <- function(name, value, feature = "global") {
      variables <<- c(variables, name)
      truth <<- c(truth, value)
      features <<- c(features, feature)
    }
    for (g in seq_len(row$n_features)) {
      feature <- sprintf("feature_%02d", g)
      for (q in 1:2) add(sprintf("beta[%d,%d]", q, g), reference$truth$beta[q, g], feature)
      add(sprintf("gamma[1,%d]", g), reference$truth$gamma[1, g], feature)
      add(sprintf("sigma_sq[%d]", g), reference$truth$sigma_sq[g], feature)
      if (row$design != "fixed") {
        for (q in seq_len(ncol(reference$truth$random$omega))) {
          add(sprintf("omega_1[%d,%d]", g, q), reference$truth$random$omega[g, q], feature)
        }
        if (row$design == "slope") add(sprintf("L_1[%d,2,1]", g), reference$truth$random$correlation[g], feature)
      }
    }
    add("log_s0", reference$truth$log_s0)
    add("nu0", reference$truth$nu0)
    if (row$family == "student_t") add("nu_obs", reference$truth$nu_obs)
    if (row$design != "fixed") {
      for (q in seq_along(reference$truth$random$mu_omega)) {
        add(sprintf("mu_omega_1[%d]", q), reference$truth$random$mu_omega[q])
        add(sprintf("tau_omega_1[%d]", q), reference$truth$random$tau_omega[q])
      }
    }
    draws <- unclass(attr(fit, "fit"))[, , variables, drop = FALSE]
    set.seed(row$fit_seed + 100000L)
    stride <- floor(dim(draws)[1L] / row$rank_draws_per_chain)
    rank_indices <- lapply(seq_len(dim(draws)[2L]), function(chain) {
      sample.int(stride, 1L) + stride * (seq_len(row$rank_draws_per_chain) - 1L)
    })
    metrics <- lapply(seq_along(variables), function(j) {
      full <- draws[, , j]
      thin_chains <- vapply(seq_along(rank_indices), function(chain) {
        full[rank_indices[[chain]], chain]
      }, numeric(row$rank_draws_per_chain))
      thin <- as.vector(thin_chains)
      autocorrelation <- apply(thin_chains, 2L, function(x) {
        as.numeric(stats::acf(x, lag.max = 3L, plot = FALSE)$acf)[-1L]
      })
      base_name <- sub("\\[.*", "", variables[j])
      # Keep coefficient slots distinct without treating features as independent runs.
      parameter <- if (base_name %in% c("beta", "gamma")) sub(
        ",[0-9]+\\]$", "]", variables[j]) else if (grepl("^omega", base_name)) {
        sub("\\[[0-9]+,", "[", variables[j])
      } else if (grepl("^L_", base_name)) base_name else variables[j]
      data.frame(variable = variables[j], feature = features[j], parameter = parameter,
        truth = truth[j], estimate = mean(full),
        covered_90 = truth[j] >= stats::quantile(full, .05) && truth[j] <= stats::quantile(full, .95),
        covered_95 = truth[j] >= stats::quantile(full, .025) && truth[j] <= stats::quantile(full, .975),
        rank = sum(thin < truth[j]), rank_draws = length(thin),
        rank_ess_bulk = posterior::ess_bulk(thin_chains),
        max_abs_rank_acf = max(abs(autocorrelation)),
        width_95 = diff(stats::quantile(full, c(.025, .975))),
        rhat = posterior::rhat(full), ess_bulk = posterior::ess_bulk(full),
        ess_tail = posterior::ess_tail(full))
    })
    receipt$metrics <- do.call(rbind, metrics)
    receipt$rank_indices <- rank_indices
    receipt$rank_thinning_review <- "Inspect per-chain lag 1-3 autocorrelation and rank ESS across independent datasets before interpreting ranks; short-series diagnostics are noisy."
    if (receipt$status == "hmc_passed") {
      tested <- sccomp_test(fit, effect_threshold = row$effect_threshold,
                           variance_ratio_threshold = row$variance_ratio_threshold)
      saveRDS(tested, file.path(fit_directory, "tested.rds"))
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
status <- data.frame(id = manifest$id, family = manifest$family, design = manifest$design,
  status = vapply(receipts, `[[`, character(1), "status"),
  elapsed_seconds = vapply(receipts, `[[`, numeric(1), "elapsed_seconds"),
  error = vapply(receipts, function(x) if (is.null(x$error)) "" else x$error, character(1)))
utils::write.csv(status, file.path(campaign, "fit-status.csv"), row.names = FALSE)
metrics <- lapply(seq_along(receipts), function(i) {
  if (is.null(receipts[[i]]$metrics)) return(NULL)
  cbind(manifest[i, c("id", "stage", "family", "design", "replication")], receipts[[i]]$metrics)
})
metrics <- do.call(rbind, metrics)
if (!is.null(metrics)) utils::write.csv(metrics, file.path(campaign, "sbc-metrics.csv"), row.names = FALSE)
report <- list(stage = unique(manifest$stage), status = status, metrics = metrics,
  attempted = nrow(manifest), hmc_passed = sum(status$status == "hmc_passed"),
  scientific_qualified = FALSE,
  reason = "Pilot/campaign execution alone does not establish coverage, SBC rank, discovery or performance gates.")
saveRDS(report, output)
