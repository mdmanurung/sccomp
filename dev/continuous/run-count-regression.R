# Refit immutable inherited fixtures through the current released count API.
root <- normalizePath(".")
.libPaths(c(file.path(root, ".local/continuous/library"), .libPaths()))
pkgload::load_all(root, quiet = TRUE, helpers = FALSE)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, grepl("^[a-zA-Z0-9_-]+$", args[2L]))
baseline <- normalizePath(args[1L])
run <- file.path(root, "dev/continuous/runs", args[2L])
if (dir.exists(run)) stop("Run exists; preserve earlier receipts.")
dir.create(run, recursive = TRUE)
assignInNamespace("sccomp_stan_models_cache_dir",
                  file.path(root, ".local/continuous/baseline-cache"), ns = "sccomp")
writeLines(capture.output(sessionInfo()), file.path(run, "session.txt"))
files <- c(list.files("R", full.names = TRUE),
           list.files("inst/stan", recursive = TRUE, full.names = TRUE))
files <- files[!dir.exists(files)]
writeLines(capture.output(tools::md5sum(files)), file.path(run, "source-md5.txt"))
for (name in c("my_estimate", "my_estimate_with_variance", "my_estimate_intercept_only")) {
  previous_file <- file.path(baseline, paste0(name, ".rds"))
  previous <- readRDS(previous_file)
  seed <- attr(previous, "fit")$metadata()$seed
  arguments <- list(.data = attr(previous, "count_data"),
    formula_composition = attr(previous, "formula_composition"),
    formula_variability = attr(previous, "formula_variability"),
    sample = rlang::as_name(attr(previous, ".sample")),
    cell_group = rlang::as_name(attr(previous, ".cell_group")),
    abundance = rlang::as_name(attr(previous, ".count")),
    inference_method = "pathfinder", cores = 1L, mcmc_seed = seed,
    max_sampling_iterations = 1000L, verbose = FALSE,
    output_directory = file.path(run, name))
  current <- do.call(sccomp_estimate, arguments)
  saveRDS(current, file.path(run, paste0(name, ".rds")))
  stopifnot(identical(attr(previous, "model_input"), attr(current, "model_input")),
            identical(lapply(previous, identity), lapply(current, identity)),
            identical(attr(previous, "fit")$draws(format = "matrix"),
                      attr(current, "fit")$draws(format = "matrix")))
  old_test <- sccomp_test(previous)
  new_test <- sccomp_test(current)
  stopifnot(identical(lapply(old_test, identity), lapply(new_test, identity)))
  saveRDS(list(previous_file = previous_file, previous_md5 = tools::md5sum(previous_file),
    seed = seed, model_input_identical = TRUE, estimates_identical = TRUE,
    draws_identical = TRUE, test_statistics_identical = TRUE),
    file.path(run, paste0(name, "-comparison.rds")))
  cat(name, "count input, draws, estimates and test statistics identical.\n")
}
writeLines("Three inherited packaged-example count fixtures reproduced exactly; source provenance is undocumented.",
           file.path(run, "count-regression.txt"))
