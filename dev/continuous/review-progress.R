# A receipt snapshot reports pending work without creating the campaign's final report.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, grepl("^[a-zA-Z0-9_-]+$", args[2L]))
campaign <- normalizePath(args[1L])
output <- file.path("dev/continuous/runs", args[2L])
if (dir.exists(output)) stop("Progress run exists; preserve it.")
manifest <- utils::read.csv(file.path(campaign, "manifest.csv"),
                            colClasses = c(student_t_df = "character"))
stopifnot(all(startsWith(manifest$stage, "review_")), !anyDuplicated(manifest$id))
hashes <- readRDS(file.path(campaign, "source-md5.rds"))
source_root <- readLines(file.path(campaign, "paths.txt"))[2L]
stopifnot(identical(unname(tools::md5sum(file.path(source_root, names(hashes)))), unname(hashes)))
read_receipts <- function(ids, folder) lapply(ids, function(id) {
  file <- file.path(campaign, folder, paste0(id, ".rds"))
  if (!file.exists(file)) return(list(status = "pending", elapsed_seconds = NA_real_))
  receipt <- readRDS(file)
  stopifnot(identical(receipt$source_md5, hashes),
            id %in% c(receipt$manifest$id, receipt$manifest$dataset_id))
  receipt
})
fits <- read_receipts(manifest$id, "receipts")
datasets <- unique(manifest$dataset_id)
comparisons <- read_receipts(datasets, "comparators")
status <- function(receipts, ids) data.frame(id = ids,
  status = vapply(receipts, `[[`, character(1), "status"),
  elapsed_seconds = vapply(receipts, `[[`, numeric(1), "elapsed_seconds"),
  error = vapply(receipts, function(x) if (is.null(x$error)) "" else x$error, character(1)))
bind <- function(receipts, field) do.call(rbind, lapply(receipts, function(x) {
  if (is.null(x[[field]])) return(NULL)
  cbind(x$manifest[rep(1L, nrow(x[[field]])), ], x[[field]])
}))
report <- list(manifest = manifest, status = status(fits, manifest$id),
  comparator_status = status(comparisons, datasets),
  metrics = bind(fits, "metrics"), features = bind(fits, "features"),
  comparator_metrics = bind(comparisons, "metrics"),
  comparator_variants = bind(comparisons, "comparator_variants"),
  nuisance_metrics = bind(fits, "nuisance_metrics"),
  partial = any(vapply(c(fits, comparisons), function(x) x$status == "pending", logical(1))),
  scientific_qualified = FALSE, source_md5 = hashes, captured = Sys.time())
stopifnot(identical(report$status$id, manifest$id),
          nrow(report$comparator_status) == length(datasets))
source("dev/continuous/assess-review.R")
bounds <- utils::read.csv(file.path(campaign, "precomputed-fixed-truth-bounds.csv"),
                          colClasses = c(student_t_df = "character"))
assessment <- assess_review(report, bounds)
dir.create(output, recursive = TRUE)
saveRDS(report, file.path(output, "report.rds"))
saveRDS(assessment, file.path(output, "assessment.rds"))
utils::write.csv(report$status, file.path(output, "fit-status.csv"), row.names = FALSE)
utils::write.csv(report$comparator_status, file.path(output, "comparator-status.csv"), row.names = FALSE)
for (name in c("metrics", "comparator_variants", "nuisance_metrics")) {
  if (!is.null(report[[name]])) utils::write.csv(report[[name]],
    file.path(output, paste0(name, ".csv")), row.names = FALSE)
}
saveRDS(list(script_md5 = tools::md5sum(c("dev/continuous/review-progress.R",
  "dev/continuous/assess-review.R", "dev/continuous/assess-fixed-truth.R")),
  campaign = campaign, runtime = sessionInfo()), file.path(output, "provenance.rds"))
print(table(report$status$status))
cat("Partial:", report$partial, "| Scientific qualification:", assessment$scientific_qualified, "\n")
