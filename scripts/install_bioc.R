#!/usr/bin/env Rscript
# Install R/Bioconductor dependencies for scripts/analyze_gse196006.R.
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", ask = FALSE, update = FALSE)
}
cran_packages <- c("ggplot2")
bioc_packages <- c("DESeq2", "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler")
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran, ask = FALSE, update = FALSE)
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
required <- c(cran_packages, bioc_packages)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Package installation incomplete: ", paste(missing, collapse = ", "))
message("All analysis R/Bioconductor dependencies are installed.")
