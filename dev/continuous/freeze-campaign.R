# Freeze package/scripts/seeds before running; never load the mutable checkout.
root <- normalizePath(".")
.libPaths(c(file.path(root, ".local/continuous/library"), .libPaths()))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, grepl("^[a-zA-Z0-9_-]+$", args[1L]),
          args[2L] %in% c("pilot", "sbc", "fixed_pilot", "fixed_truth", "review_pilot", "review_truth"))
campaign <- file.path(root, ".local/continuous/campaigns", args[1L])
if (dir.exists(campaign)) stop("Campaign exists; never overwrite a frozen run.")
dir.create(campaign, recursive = TRUE)
package_source <- file.path(campaign, "source/sccomp")
files <- system2("git", c("ls-files", "-co", "--exclude-standard"), stdout = TRUE)
files <- files[!dir.exists(files) & !grepl("^dev/continuous/runs/", files)]
for (file in files) {
  target <- file.path(package_source, file)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file, target))
}
source_hashes <- tools::md5sum(file.path(package_source, files))
names(source_hashes) <- files
saveRDS(source_hashes, file.path(campaign, "source-md5.rds"))
utils::write.csv(data.frame(file = files, md5 = unname(source_hashes)),
                 file.path(campaign, "source-md5.csv"), row.names = FALSE)
setwd(dirname(package_source))
utils::tar(file.path(campaign, "source.tar.gz"), files = "sccomp",
            compression = "gzip", tar = "internal")
setwd(root)
library_path <- file.path(campaign, "library")
dir.create(library_path)
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = ":"))
status <- system2(file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "--no-multiarch", paste0("--library=", shQuote(library_path)),
    shQuote(package_source)), stdout = file.path(campaign, "install.log"),
  stderr = file.path(campaign, "install.log"))
stopifnot(status == 0L)
rows <- list()
fixed_truth <- args[2L] %in% c("fixed_pilot", "fixed_truth", "review_pilot", "review_truth")
seed_offset <- if (args[2L] == "pilot") 1000000L else 2000000L
if (!fixed_truth) for (family in c("gaussian", "student_t")) {
  for (design in c("fixed", "intercept", "slope")) {
    repetitions <- if (args[2L] == "pilot") 1L else if (design == "fixed") 200L else 100L
    for (replication in seq_len(repetitions)) {
      index <- length(rows) + 1L
      rows[[index]] <- data.frame(
        id = sprintf("%s-%s-%s-%04d", args[2L], family, design, replication),
        stage = args[2L], family = family, design = design, replication = replication,
        seed = 10102026L + index + seed_offset,
        fit_seed = 20102026L + index + seed_offset,
        n_features = 20L, n_donors = 20L, variance_effect = TRUE,
        observations_per_donor = if (design == "slope") 4L else 2L,
        chains = 4L, warmup = 1500L, sampling = 1500L,
        adapt_delta = .995, max_treedepth = 12L,
        effect_threshold = 1, variance_ratio_threshold = 1.5,
        rank_draws_per_chain = 25L, protocol = "v2-identifiable-slopes"
      )
    }
  }
}
if (fixed_truth) {
  source("dev/continuous/fixed-truth-campaign.R")
  if (startsWith(args[2L], "review_")) {
    source("dev/continuous/review-study.R")
    manifest <- review_manifest(args[2L])
  } else manifest <- fixed_truth_manifest(args[2L])
} else manifest <- do.call(rbind, rows)
utils::write.csv(manifest, file.path(campaign, "manifest.csv"), row.names = FALSE)
if (!fixed_truth) {
bounds <- unique(manifest[c("family", "design")])
bounds$datasets <- vapply(seq_len(nrow(bounds)), function(i) {
  sum(manifest$family == bounds$family[i] & manifest$design == bounds$design[i])
}, integer(1))
bounds <- bounds[rep(seq_len(nrow(bounds)), each = 2L), ]
bounds$nominal_coverage <- rep(c(.90, .95), nrow(bounds) / 2L)
bounds$coverage_lower_count <- stats::qbinom(.05 / (4 * 91), bounds$datasets, bounds$nominal_coverage)
bounds$coverage_upper_count <- stats::qbinom(1 - .05 / (4 * 91), bounds$datasets, bounds$nominal_coverage)
bounds$rank_ecdf_envelope <- sqrt(log(2 * 91 / .05) / (2 * bounds$datasets))
bounds$diagnostic_strata <- 91L
}
bound_file <- if (fixed_truth) "precomputed-fixed-truth-bounds.csv" else "precomputed-sbc-bounds.csv"
if (fixed_truth) bounds <- fixed_truth_bounds(manifest)
utils::write.csv(bounds, file.path(campaign, bound_file), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(campaign, "freeze-session.txt"))
writeLines(capture.output(utils::installed.packages()[, c("Package", "Version", "LibPath")]),
           file.path(campaign, "dependency-versions.txt"))
writeLines(c("source_root", package_source, "dependency_library",
             file.path(root, ".local/continuous/library"), "Rscript", file.path(R.home("bin"), "Rscript")),
           file.path(campaign, "paths.txt"))
Sys.chmod(file.path(package_source, files), "0444")
Sys.chmod(file.path(campaign, c("manifest.csv", "source-md5.rds", "source-md5.csv", "source.tar.gz",
                                bound_file)), "0444")
cat("Frozen", campaign, "with", nrow(manifest), "planned fits; no jobs launched.\n")
