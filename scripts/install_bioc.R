#!/usr/bin/env Rscript
# Install the packages required for the core GSE196006 DESeq2 analysis.
# GO enrichment is optional: the analysis script deliberately continues with
# an empty enrichment table if the optional clusterProfiler stack is unavailable.
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", ask = FALSE, update = FALSE)
}
options(repos = BiocManager::repositories())

required_cran <- c("ggplot2")
required_bioc <- c("DESeq2")
missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, ask = FALSE, update = FALSE)
}
missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, dependencies = TRUE, ask = FALSE, update = FALSE)
}

required <- c(required_cran, required_bioc)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Required package installation incomplete: ", paste(missing, collapse = ", "))
}
message("Core DESeq2 analysis dependencies are installed.")

# Optional pathway-enrichment stack. Failure here must not block differential
# expression, plots, metadata checks, or provenance outputs.
optional_bioc <- c("AnnotationDbi", "org.Hs.eg.db", "clusterProfiler")
missing_optional <- optional_bioc[!vapply(optional_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_optional)) {
  message("Attempting optional GO-enrichment dependencies: ", paste(missing_optional, collapse = ", "))
  tryCatch(
    BiocManager::install(missing_optional, ask = FALSE, update = FALSE),
    error = function(e) message("Optional GO-enrichment dependencies could not be installed: ", conditionMessage(e))
  )
}
still_missing <- optional_bioc[!vapply(optional_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) {
  message("WARNING: optional GO enrichment will be skipped; unavailable package(s): ",
          paste(still_missing, collapse = ", "))
} else {
  message("Optional GO-enrichment dependencies are installed.")
}
