# Build vignettes and check the complete package from a verified frozen source.
root <- normalizePath(".")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, grepl("^[a-zA-Z0-9_-]+$", args[2L]))
snapshot <- normalizePath(args[1L])
paths <- readLines(file.path(snapshot, "paths.txt"))
source_root <- paths[2L]
hashes <- readRDS(file.path(snapshot, "source-md5.rds"))
stopifnot(identical(unname(tools::md5sum(file.path(source_root, names(hashes)))), unname(hashes)))
run <- file.path(root, "dev/continuous/runs", args[2L])
if (dir.exists(run)) stop("Check run exists; preserve its receipt.")
dir.create(run, recursive = TRUE)
.libPaths(c(file.path(snapshot, "library"), paths[4L], .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = ":"), NOT_CRAN = "true",
  SCCOMP_STAN_CACHE_DIR = file.path(root, ".local/continuous/baseline-cache"),
  MAKEFLAGS = "PRECOMPILED_HEADERS=false")
saveRDS(list(source_root = source_root, source_md5 = hashes, runtime = sessionInfo(),
  environment = Sys.getenv(c("R_LIBS", "NOT_CRAN", "SCCOMP_STAN_CACHE_DIR", "MAKEFLAGS", "PATH", "RSTUDIO_PANDOC")),
  pandoc = rmarkdown::find_pandoc()),
  file.path(run, "source-runtime.rds"))
writeLines(c("R CMD build SOURCE", "R CMD check --no-multiarch --timings sccomp_2.5.0.tar.gz",
             "Vignettes, tests, examples and manual are enabled."), file.path(run, "commands.txt"))
setwd(run)
started <- Sys.time()
build <- system2(file.path(R.home("bin"), "R"), c("CMD", "build", shQuote(source_root)),
                 stdout = "build.log", stderr = "build.log")
saveRDS(list(build_exit = build, started = started, build_finished = Sys.time()), "build-status.rds")
stopifnot(build == 0L)
archive <- list.files(pattern = "^sccomp_.*[.]tar[.]gz$")
stopifnot(length(archive) == 1L)
check <- system2(file.path(R.home("bin"), "R"),
  c("CMD", "check", "--no-multiarch", "--timings", shQuote(archive)),
  stdout = "check.log", stderr = "check.log")
saveRDS(list(build_exit = build, check_exit = check, archive_md5 = tools::md5sum(archive),
  started = started, finished = Sys.time(),
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))), "check-status.rds")
stopifnot(check == 0L)
cat("Package check command completed; inspect check.log for warnings/notes before qualification.\n")
