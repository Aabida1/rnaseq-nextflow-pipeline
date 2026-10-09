if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c("DESeq2", "org.Hs.eg.db", "AnnotationDbi"), ask = FALSE, update = FALSE)
install.packages("ggplot2", repos = "https://cloud.r-project.org")
# Enrichment is optional: clusterProfiler has a larger dependency tree and can fail to compile on hosted runners.
tryCatch(
  BiocManager::install("clusterProfiler", ask = FALSE, update = FALSE),
  error = function(e) message("Optional clusterProfiler install failed: ", conditionMessage(e))
)
