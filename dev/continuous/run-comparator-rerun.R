# Rerun comparators on a frozen campaign's inputs with the current comparators.R.
# Never writes into the campaign; failed original receipts stay in its denominator.
# Usage: run-comparator-rerun.R CAMPAIGN RUN-ID [plain|robust|trend_robust]
args <- commandArgs(trailingOnly = TRUE)
variants <- list(plain = list(), robust = list(robust = TRUE), trend_robust = list(trend = TRUE, robust = TRUE))
variant <- if (length(args) == 3L) args[3L] else "plain"
stopifnot(length(args) %in% 2:3, grepl("^[a-zA-Z0-9_-]+$", args[2L]), variant %in% names(variants))
campaign <- normalizePath(args[1L])
out <- file.path("dev/continuous/runs", args[2L])
if (dir.exists(out)) stop("Run exists; never overwrite.")
dir.create(out, recursive = TRUE)
.libPaths(c(file.path(campaign, "library"), readLines(file.path(campaign, "paths.txt"))[4L], .libPaths()))
sources <- file.path("dev/continuous", c("simulate-fixed-truth.R", "discovery-metrics.R", "comparators.R"))
for (file in sources) source(file)
utils::write.csv(data.frame(file = sources, md5 = unname(tools::md5sum(sources))),
                 file.path(out, "source-md5.csv"), row.names = FALSE)
rows <- lapply(list.files(file.path(campaign, "inputs"), full.names = TRUE), function(file) {
  input <- readRDS(file)
  id <- input$dataset$dataset_id
  started <- Sys.time()
  receipt <- list(manifest = input$dataset, input_md5 = unname(tools::md5sum(file)), status = "fit_failed")
  error <- tryCatch({
    comparison <- do.call(fit_continuous_comparator,
      c(list(input$reference$data, input$dataset$design), variants[[variant]]))
    receipt$features <- comparison$results
    receipt$variant <- variant
    receipt$metrics <- comparator_dataset_metrics(comparison$results, input$reference, input$dataset$n_features)
    receipt$status <- "comparator_passed"
    NULL
  }, error = conditionMessage)
  receipt$error <- error
  receipt$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  saveRDS(receipt, file.path(out, paste0(id, ".rds")))
  original <- file.path(campaign, "comparators", paste0(id, ".rds"))
  data.frame(dataset_id = id, status = receipt$status, error = if (is.null(error)) "" else error,
             boundary_features = sum(receipt$features$boundary, na.rm = TRUE), original_receipt = original,
             original_status = readRDS(original)$status)
})
links <- do.call(rbind, rows)
utils::write.csv(links, file.path(out, "links.csv"), row.names = FALSE)
print(links[c("dataset_id", "status", "boundary_features", "original_status")])
