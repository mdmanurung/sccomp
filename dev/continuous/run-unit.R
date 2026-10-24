root <- normalizePath(".")
.libPaths(c(file.path(root, ".local/continuous/library"), .libPaths()))
pkgload::load_all(root, quiet = TRUE, helpers = FALSE)
run_id <- commandArgs(trailingOnly = TRUE)
if (!length(run_id)) run_id <- format(Sys.time(), "unit-%Y%m%dT%H%M%S", tz = "UTC")
stopifnot(length(run_id) == 1L, grepl("^[a-zA-Z0-9_-]+$", run_id))
run <- file.path(root, "dev/continuous/runs", run_id)
if (dir.exists(run)) stop("Run exists; use a new identity.")
dir.create(run, recursive = TRUE)
files <- c(list.files("R", "continuous", full.names = TRUE),
           list.files("tests/testthat", "continuous", full.names = TRUE))
writeLines(capture.output(tools::md5sum(files)), file.path(run, "source-md5.txt"))
writeLines(capture.output(sessionInfo()), file.path(run, "session.txt"))
result <- testthat::test_dir("tests/testthat", filter = "continuous", load_helpers = FALSE,
                            env = new.env(parent = asNamespace("sccomp")),
                            reporter = "summary", stop_on_failure = FALSE)
saveRDS(result, file.path(run, "test-results.rds"))
summary <- as.data.frame(result)
utils::write.csv(summary[!vapply(summary, is.list, logical(1))],
                 file.path(run, "test-results.csv"), row.names = FALSE)
stopifnot(!any(summary$failed > 0), !any(summary$error))
