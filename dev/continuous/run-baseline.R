# Run the full inherited suite against the immutable source snapshot.
root <- normalizePath(".")
baseline <- file.path(root, "dev/continuous/runs/baseline-2026-10-10")
run_id <- commandArgs(trailingOnly = TRUE)
if (!length(run_id)) run_id <- format(Sys.time(), "baseline-%Y%m%dT%H%M%S", tz = "UTC")
stopifnot(length(run_id) == 1L, grepl("^[a-zA-Z0-9_-]+$", run_id))
run <- file.path(root, "dev/continuous/runs", run_id)
if (dir.exists(run)) stop("Run directory already exists; use a new run identity.")
dir.create(run, recursive = TRUE)
local <- file.path(root, ".local/continuous")
lib <- file.path(local, "baseline-library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, file.path(local, "library"), .libPaths()))
Sys.setenv(NOT_CRAN = "true", R_LIBS = paste(.libPaths(), collapse = ":"))
stopifnot(packageVersion("instantiate") >= "0.2.3")
stopifnot(packageVersion("cmdstanr") >= "0.9.0")
source_root <- file.path(local, "baseline-source")
if (!dir.exists(source_root)) {
  dir.create(source_root)
  utils::untar(file.path(baseline, "source.tar.gz"), exdir = source_root)
}
source <- file.path(source_root, "sccomp")
if (!file.exists(file.path(lib, "sccomp/DESCRIPTION"))) {
  install_status <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-multiarch", paste0("--library=", shQuote(lib)),
      shQuote(source)),
    stdout = file.path(run, "package-install.log"),
    stderr = file.path(run, "package-install.log")
  )
  stopifnot(install_status == 0L)
}
library(sccomp)
cache <- file.path(local, "baseline-cache")
assignInNamespace("sccomp_stan_models_cache_dir", cache, ns = "sccomp")
writeLines(capture.output(sessionInfo()), file.path(run, "session.txt"))
writeLines(c(cmdstanr::cmdstan_path(), as.character(cmdstanr::cmdstan_version())),
           file.path(run, "cmdstan.txt"))
start <- Sys.time()
testthat::set_max_fails(Inf)
test_env <- new.env(parent = asNamespace("sccomp"))
result <- testthat::test_dir(
  file.path(source, "tests/testthat"),
  env = test_env, reporter = "summary", stop_on_failure = FALSE,
  stop_on_warning = FALSE
)
saveRDS(result, file.path(run, "test-results.rds"))
summary <- as.data.frame(result)
utils::write.csv(summary[!vapply(summary, is.list, logical(1))],
                 file.path(run, "test-results.csv"), row.names = FALSE)
writeLines(as.character(difftime(Sys.time(), start, units = "secs")),
           file.path(run, "elapsed-seconds.txt"))
for (object in c("my_estimate", "my_estimate_with_variance",
                 "my_estimate_intercept_only")) {
  if (exists(object, envir = test_env, inherits = FALSE)) {
    saveRDS(get(object, envir = test_env),
            file.path(run, paste0(object, ".rds")))
  }
}
stopifnot(all(summary$failed == 0L), !any(summary$error))
