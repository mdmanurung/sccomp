# Run all tests against a verified, installed immutable source snapshot.
root <- normalizePath(".")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, grepl("^[a-zA-Z0-9_-]+$", args[2L]))
snapshot <- normalizePath(args[1L])
paths <- readLines(file.path(snapshot, "paths.txt"))
source_root <- paths[2L]
hashes <- readRDS(file.path(snapshot, "source-md5.rds"))
stopifnot(identical(unname(tools::md5sum(file.path(source_root, names(hashes)))), unname(hashes)))
.libPaths(c(file.path(snapshot, "library"), paths[4L], .libPaths()))
Sys.setenv(NOT_CRAN = "true", R_LIBS = paste(.libPaths(), collapse = ":"))
library(sccomp)
run <- file.path(root, "dev/continuous/runs", args[2L])
if (dir.exists(run)) stop("Run exists; preserve the earlier receipt.")
dir.create(run, recursive = TRUE)
cache <- file.path(root, ".local/continuous/baseline-cache")
assignInNamespace("sccomp_stan_models_cache_dir", cache, ns = "sccomp")
saveRDS(list(source_root = source_root, source_md5 = hashes,
             library = find.package("sccomp")), file.path(run, "source.rds"))
writeLines(capture.output(sessionInfo()), file.path(run, "session.txt"))
testthat::set_max_fails(Inf)
started <- Sys.time()
result <- testthat::test_dir(file.path(source_root, "tests/testthat"),
  env = new.env(parent = asNamespace("sccomp")), reporter = "summary",
  stop_on_failure = FALSE, stop_on_warning = FALSE)
saveRDS(result, file.path(run, "test-results.rds"))
summary <- as.data.frame(result)
utils::write.csv(summary[!vapply(summary, is.list, logical(1))],
                 file.path(run, "test-results.csv"), row.names = FALSE)
writeLines(as.character(difftime(Sys.time(), started, units = "secs")),
           file.path(run, "elapsed-seconds.txt"))
stopifnot(all(summary$failed == 0L), !any(summary$error))
