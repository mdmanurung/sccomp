args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
archive <- normalizePath(args[1L])
output <- args[2L]
if (file.exists(output)) stop("The frozen application input already exists.")
root <- normalizePath(".")
.libPaths(c(file.path(root, ".local/continuous/library"), .libPaths()))
suppressPackageStartupMessages(library(airway))
data("airway", package = "airway")
stopifnot(utils::packageVersion("airway") == "1.30.0")
counts <- as.matrix(SummarizedExperiment::assay(airway))
metadata <- as.data.frame(SummarizedExperiment::colData(airway))
stopifnot(!anyNA(counts), all(is.finite(counts)), all(counts >= 0),
          identical(colnames(counts), rownames(metadata)))
sample_data <- data.frame(sample = rownames(metadata),
  condition = factor(ifelse(metadata$dex == "trt", "treated", "control"),
                       levels = c("control", "treated")),
  donor = factor(metadata$cell))
stopifnot(all(table(sample_data$donor, sample_data$condition) == 1L))
nonzero <- rowSums(counts) > 0
dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts[nonzero, ]), method = "TMM")
cpm <- edgeR::cpm(dge, log = FALSE)
eligible <- rownames(cpm)[rowSums(cpm >= 1 & dge$counts >= 10) == ncol(cpm)]
features <- head(sort(eligible), 100L)
stopifnot(length(features) == 100L)
values <- edgeR::cpm(dge, log = TRUE, prior.count = 2)[features, , drop = FALSE]
stopifnot(all(is.finite(values)))
long <- data.frame(sample = rep(sample_data$sample, each = length(features)),
                     feature = rep(features, nrow(sample_data)), value = as.vector(values))
long <- merge(long, sample_data, by = "sample", sort = FALSE)
protocol <- readLines("dev/continuous/APPLICATION.md")
provenance <- list(
  source = "https://bioconductor.org/packages/3.22/data/experiment/src/contrib/airway_1.30.0.tar.gz",
  source_sha256 = system2("sha256sum", shQuote(archive), stdout = TRUE),
  source_md5 = tools::md5sum(archive), package_version = utils::packageVersion("airway"),
  citation = utils::citation("airway"), geo = "GSE52778", package_license = "LGPL",
  protocol = protocol, runtime = sessionInfo(), prepared = Sys.time(),
  assay_units = "TMM-normalized log2 CPM, prior.count=2",
  independent_unit = "primary cell line", feature_rule = "first 100 lexical eligible Ensembl IDs",
  feature_filter = "CPM >=1 and raw count >=10 in all eight samples",
  excluded_samples = character(), imputation = "none"
)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(data = long, values = values, sample_data = sample_data,
             original_metadata = metadata, normalization = dge$samples,
             features = features, provenance = provenance), output)
utils::write.csv(sample_data, file.path(dirname(output), "sample-metadata.csv"), row.names = FALSE)
writeLines(features, file.path(dirname(output), "selected-features.txt"))
writeLines(capture.output(provenance[c("source", "source_sha256", "package_version", "assay_units", "independent_unit")]),
           file.path(dirname(output), "provenance.txt"))
Sys.chmod(output, "0444")
cat("Prepared", nrow(sample_data), "samples,", nlevels(sample_data$donor),
    "independent cell lines and", length(features), "features. No fitting performed.\n")
