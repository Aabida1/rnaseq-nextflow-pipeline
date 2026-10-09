if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi"), ask = FALSE, update = FALSE)
install.packages("ggplot2", repos = "https://cloud.r-project.org")
