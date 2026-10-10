#!/usr/bin/env Rscript
# Reproducible installer for the GSE196006 analysis.
# Pin Bioconductor 3.21 to match R 4.5 in .github/workflows/ci.yml.
# This avoids mixing Bioconductor releases when repositories advance.
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", ask = FALSE, update = FALSE)
}
BiocManager::install(version = "3.21", ask = FALSE, update = FALSE)
options(repos = BiocManager::repositories(version = "3.21"))

required_cran <- c("ggplot2")
required_bioc <- c("DESeq2")
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, ask = FALSE, update = FALSE)
  }
}
missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, version = "3.21", dependencies = TRUE,
                       ask = FALSE, update = FALSE)
}
required <- c(required_cran, required_bioc)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Required package installation incomplete: ", paste(missing, collapse = ", "))
}
message("Core dependencies installed using R ", getRversion(),
        " and Bioconductor ", as.character(BiocManager::version()), ".")

# GO enrichment is optional. It must never block the core DESeq2 analysis.
optional_bioc <- c("AnnotationDbi", "org.Hs.eg.db", "clusterProfiler")
missing_optional <- optional_bioc[!vapply(optional_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_optional)) {
  message("Attempting optional GO-enrichment dependencies: ", paste(missing_optional, collapse = ", "))
  tryCatch(
    BiocManager::install(missing_optional, version = "3.21", dependencies = TRUE,
                         ask = FALSE, update = FALSE),
    error = function(e) message("Optional GO-enrichment install failed: ", conditionMessage(e))
  )
}
still_missing <- optional_bioc[!vapply(optional_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) {
  message("WARNING: GO enrichment will be skipped; unavailable package(s): ",
          paste(still_missing, collapse = ", "))
} else {
  message("Optional GO-enrichment dependencies are installed.")
}
