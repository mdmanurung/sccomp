# Preserve the failed fit and change numerical settings on the same observed data.
root <- normalizePath(".")
.libPaths(c(file.path(root, ".local/continuous/library"), .libPaths()))
pkgload::load_all(root, quiet = TRUE, helpers = FALSE)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, grepl("^[a-zA-Z0-9_-]+$", args[2L]))
previous_file <- normalizePath(args[1L])
previous <- readRDS(previous_file)
meta <- attr(previous, "continuous")
run <- file.path(root, "dev/continuous/runs", args[2L])
if (dir.exists(run)) stop("Run exists; never replace an earlier fit.")
dir.create(run, recursive = TRUE)
saveRDS(list(previous_file = previous_file, previous_md5 = tools::md5sum(previous_file),
             previous_diagnostics = meta$diagnostics,
             amendment = "Same physical model/data/priors; current numerical coordinates and initialization, 2000 warmup/4000 samples."),
          file.path(run, "reason.rds"))
writeLines(capture.output(tools::md5sum(c(list.files("R", "continuous", full.names = TRUE),
  list.files("inst/stan/continuous", full.names = TRUE),
  file.path("inst/stan", c("glm_multi_continuous.stan",
                           "glm_multi_continuous_generate_data.stan"))))), file.path(run, "source-md5.txt"))
arguments <- list(
  .data = meta$input$original_data, formula_location = meta$formula_location,
  formula_variance = meta$formula_variance, sample = meta$input$sample,
  feature = meta$input$feature, value = meta$input$value, family = meta$family,
  scaling = "none", prior = meta$input$prior, student_t_df = meta$student_t_df,
  cores = 4, chains = 4, iter_warmup = 2000, iter_sampling = 4000,
  mcmc_seed = meta$seed, adapt_delta = .995, max_treedepth = 12L,
  output_directory = run, cache_stan_model = file.path(root, ".local/continuous/api-models"),
  refresh = 0L, sig_figs = 16
)
elapsed <- system.time(fit <- do.call(sccomp_estimate_continuous, arguments))
saveRDS(fit, file.path(run, "fit.rds"))
writeLines(capture.output(elapsed), file.path(run, "elapsed.txt"))
diagnostics <- attr(fit, "continuous")$diagnostics
saveRDS(diagnostics, file.path(run, "diagnostics.rds"))
cat("Reference refit diagnostic gate:", diagnostics$hmc_qualified, "\n")
print(diagnostics$chains)
failed <- with(diagnostics$parameters, !is.finite(rhat) | !is.finite(ess_bulk) |
  !is.finite(ess_tail) | rhat >= 1.01 | ess_bulk < 400 | ess_tail < 400)
print(diagnostics$parameters[failed, ])
