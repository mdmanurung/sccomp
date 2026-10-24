# Compilation and routing checks; short chains are not scientific qualification.
root <- normalizePath(".")
.libPaths(c(file.path(root, ".local/continuous/baseline-library"),
            file.path(root, ".local/continuous/library"), .libPaths()))
library(sccomp)
source("R/sccomp_continuous_input.R")
formula_to_random_effect_formulae <- getFromNamespace(
  "formula_to_random_effect_formulae", "sccomp"
)
source("dev/continuous/simulate-reference.R")
run_id <- commandArgs(trailingOnly = TRUE)
if (!length(run_id)) run_id <- format(Sys.time(), "core-%Y%m%dT%H%M%S", tz = "UTC")
stopifnot(length(run_id) == 1L, grepl("^[a-zA-Z0-9_-]+$", run_id))
run <- file.path(root, "dev/continuous/runs", run_id)
if (dir.exists(run)) stop("Run directory already exists; use a new run identity.")
dir.create(run, recursive = TRUE, showWarnings = FALSE)
model_dir <- file.path(root, ".local/continuous/core-models")
dir.create(file.path(model_dir, "continuous"), recursive = TRUE,
           showWarnings = FALSE)
source_files <- c(
  list.files("inst/stan/continuous", full.names = TRUE),
  file.path("inst/stan", c("glm_multi_continuous.stan",
                           "glm_multi_continuous_generate_data.stan"))
)
writeLines(capture.output(tools::md5sum(source_files)),
           file.path(run, "source-md5.txt"))
for (source in source_files) {
  destination <- file.path(model_dir,
                           if (grepl("/continuous/", source)) "continuous" else "",
                           basename(source))
  if (!file.exists(destination) ||
      unname(tools::md5sum(source)) != unname(tools::md5sum(destination))) {
    stopifnot(file.copy(source, destination, overwrite = TRUE))
  }
}
writeLines(capture.output(sessionInfo()), file.path(run, "session.txt"))
writeLines(c(cmdstanr::cmdstan_path(), as.character(cmdstanr::cmdstan_version()),
             paste("MAKEFLAGS", Sys.getenv("MAKEFLAGS"))),
           file.path(run, "toolchain.txt"))
hashes <- tools::md5sum(source_files)
hash_receipt <- file.path(model_dir, "compiled-source-md5.rds")
recompile <- !file.exists(hash_receipt) || !identical(readRDS(hash_receipt), hashes)
models <- lapply(c("glm_multi_continuous", "glm_multi_continuous_generate_data"),
                  function(name) {
                    cmdstanr::cmdstan_model(
                      file.path(model_dir, paste0(name, ".stan")),
                      include_paths = model_dir,
                      cpp_options = list(stan_threads = TRUE),
                      force_recompile = recompile
                    )
                  })
saveRDS(hashes, hash_receipt)
saveRDS(models, file.path(run, "models.rds"))
for (family in c("gaussian", "student_t")) {
 for (design in c("fixed", "slope")) {
  case_run <- file.path(run, paste(family, design, sep = "-"))
  dir.create(case_run)
  reference <- simulate_continuous_reference(10102026L, family,
                                             design = design, n_features = 12L, n_donors = 6L)
  formula <- if (design == "fixed") ~ condition else ~ condition + (1 + condition | donor)
  input <- continuous_model_input(
    reference$data, formula, ~ 1, "sample", "feature", "value", family,
    "none", "none", "estimate", sccomp_continuous_prior(), use_data = FALSE
  )
  data <- input$stan_data
  saveRDS(list(input = input, reference = reference),
          file.path(case_run, "input.rds"))
  fit <- models[[1]]$sample(
    data = data, seed = 10102026L, chains = 4L, parallel_chains = 4L,
    threads_per_chain = 1L,
    iter_warmup = 500L, iter_sampling = 500L, adapt_delta = 0.99,
    refresh = 0L, output_dir = case_run, sig_figs = 16
  )
  diagnostics <- fit$diagnostic_summary()
  saveRDS(diagnostics, file.path(case_run, "prior-diagnostics.rds"))
  gq <- models[[2]]$generate_quantities(
    data = data, fitted_params = fit$output_files(), seed = 10102027L,
    threads_per_chain = 1L, output_dir = case_run
  )
  draws <- gq$draws(c("mu", "residual_variance", "y_rep", "log_lik"),
                     format = "matrix")
  stopifnot(all(is.finite(draws)),
            all(draws[, grepl("^residual_variance\\[", colnames(draws))] > 0))
  saveRDS(gq$summary(c("mu", "residual_variance", "y_rep", "log_lik")),
          file.path(case_run, "gq-summary.rds"))
  cat(family, design, "compiled and executed prior-only/GQ routing.\n")
 }
}
writeLines("Technical routing checks completed; NOT scientific validation.",
           file.path(run, "technical-routing.txt"))
