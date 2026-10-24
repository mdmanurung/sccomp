# Both likelihoods and a correlated mixed design, with independent GQ comparison.
root <- normalizePath(".")
.libPaths(c(file.path(root, ".local/continuous/library"), .libPaths()))
pkgload::load_all(root, quiet = TRUE, helpers = FALSE)
source("dev/continuous/simulate-reference.R")
run_id <- commandArgs(trailingOnly = TRUE)
if (!length(run_id)) run_id <- format(Sys.time(), "integration-%Y%m%dT%H%M%S", tz = "UTC")
stopifnot(length(run_id) == 1L, grepl("^[a-zA-Z0-9_-]+$", run_id))
run <- file.path(root, "dev/continuous/runs", run_id)
if (dir.exists(run)) stop("Run directory exists; choose a new identity.")
dir.create(run, recursive = TRUE)
sources <- c(list.files("R", "continuous", full.names = TRUE),
             list.files("inst/stan/continuous", full.names = TRUE),
             file.path("inst/stan", c("glm_multi_continuous.stan",
                                      "glm_multi_continuous_generate_data.stan")))
writeLines(capture.output(tools::md5sum(sources)), file.path(run, "source-md5.txt"))
writeLines(capture.output(sessionInfo()), file.path(run, "session.txt"))
for (family in c("gaussian", "student_t")) {
  design <- if (family == "gaussian") "fixed" else "slope"
  reference <- simulate_continuous_reference(10102029L, family, design,
                                             n_features = 12L, n_donors = 12L,
                                             variance_effect = TRUE)
  formula <- if (design == "fixed") ~ condition else ~ condition + (1 + condition | donor)
  family_run <- file.path(run, family)
  dir.create(family_run)
  gq_model <- sccomp:::continuous_load_model("glm_multi_continuous_generate_data", family,
    file.path(root, ".local/continuous/api-models"))$model
  saveRDS(reference, file.path(family_run, "independent-input.rds"))
  elapsed <- system.time(fit <- sccomp_estimate_continuous(
    reference$data, formula_location = formula, formula_variance = ~ condition,
    sample = sample, feature = feature, value = value, family = family,
    scaling = "none", prior = sccomp_continuous_prior(),
    cores = 4, chains = 4, iter_warmup = 2000, iter_sampling = 4000,
    mcmc_seed = 10102030L, output_directory = family_run,
    cache_stan_model = file.path(root, ".local/continuous/api-models"),
    adapt_delta = .995, max_treedepth = 12L, refresh = 0L, sig_figs = 16
  ))
  saveRDS(fit, file.path(family_run, "fit.rds"))
  stopifnot(identical(attr(fit, "continuous")$sampler_arguments$extra$init, .5),
            all(dirname(attr(fit, "continuous")$sampler_files) != family_run))
  writeLines(capture.output(elapsed), file.path(family_run, "elapsed.txt"))
  gq <- gq_model$generate_quantities(
    data = attr(fit, "model_input"), fitted_params = attr(fit, "continuous")$sampler_files,
    seed = 10102031L, threads_per_chain = 1L, output_dir = family_run, sig_figs = 16
  )
  predictions <- sccomp_predict(fit, number_of_draws = prod(dim(attr(fit, "fit"))[1:2]),
                                 mcmc_seed = 2001L)
  gq_draws <- gq$draws(c("mu", "residual_variance"), format = "draws_matrix")
  stopifnot(nrow(gq_draws) == prod(dim(attr(fit, "fit"))[1:2]))
  input <- attr(fit, "continuous")$input
  comparisons <- list()
  for (g in seq_len(input$stan_data$G)) {
    p <- predictions[predictions$feature == input$features[g], ]
    p <- p[match(input$samples, p$sample), ]
    stopifnot(nrow(p) == input$stan_data$N,
              length(p$location) == input$stan_data$N,
              length(p$residual_variance) == input$stan_data$N,
              all(is.finite(p$location)), all(is.finite(p$residual_variance)))
    stan_mu <- colMeans(gq_draws[, sprintf("mu[%d,%d]", seq_along(input$samples), g), drop = FALSE])
    stan_v <- colMeans(gq_draws[, sprintf("residual_variance[%d,%d]", seq_along(input$samples), g), drop = FALSE])
    stopifnot(max(abs(p$location - stan_mu)) < 1e-9,
              max(abs(p$residual_variance - stan_v)) < 1e-8)
    comparisons[[g]] <- data.frame(feature = input$features[g],
      max_location_difference = max(abs(p$location - stan_mu)),
      max_variance_difference = max(abs(p$residual_variance - stan_v)))
  }
  utils::write.csv(do.call(rbind, comparisons), file.path(family_run, "gq-comparison.csv"), row.names = FALSE)
  residuals <- sccomp_calculate_residuals(fit)
  stopifnot(nrow(residuals) == nrow(reference$data), all(is.finite(residuals$standardized_residual)))
  if (isTRUE(attr(fit, "continuous")$diagnostics$hmc_qualified)) {
    tested <- sccomp_test(fit, effect_threshold = 1, variance_ratio_threshold = 1.5)
    saveRDS(tested, file.path(family_run, "tested.rds"))
  } else {
    stopifnot(inherits(try(sccomp_test(fit, effect_threshold = 1), silent = TRUE), "try-error"))
  }
  cat(family, design, "API/GQ checks passed; diagnostic qualification:",
      attr(fit, "continuous")$diagnostics$hmc_qualified, "\n")
}
writeLines("Technical integration checks completed; NOT SBC or scientific qualification.",
           file.path(run, "technical-integration.txt"))
