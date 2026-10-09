#!/usr/bin/env Rscript
# Post-analysis QC and annotation for GSE196006.
# Run after analyze_gse196006.R in the same output directory:
#   Rscript scripts/qc_gse196006.R --indir results
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
indir <- "results"
if ("--indir" %in% args) {
  i <- match("--indir", args)
  if (i < length(args)) indir <- args[[i + 1]]
}
need <- c("raw_counts.tsv", "normalized_counts.tsv", "sample_metadata.tsv", "differential_expression.tsv")
missing <- need[!file.exists(file.path(indir, need))]
if (length(missing)) stop("Missing required input(s): ", paste(missing, collapse = ", "), ". Run analyze_gse196006.R first.")
suppressPackageStartupMessages({
  library(ggplot2)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
out <- file.path(indir, "qc")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
counts_df <- read.delim(file.path(indir, "raw_counts.tsv"), check.names = FALSE, row.names = 1)
norm_df <- read.delim(file.path(indir, "normalized_counts.tsv"), check.names = FALSE)
meta <- read.delim(file.path(indir, "sample_metadata.tsv"), check.names = FALSE, colClasses = "character")
de <- read.delim(file.path(indir, "differential_expression.tsv"), check.names = FALSE)
if (!all(c("sample", "patient", "condition") %in% names(meta))) stop("Metadata must contain sample, patient, condition.")
meta$patient <- factor(meta$patient)
meta$condition <- factor(tolower(meta$condition), levels = c("normal", "tumour"))
if (anyNA(meta$condition)) stop("Unexpected condition labels; expected normal and tumour.")
if (!setequal(colnames(counts_df), meta$sample)) stop("Sample names in raw counts and metadata do not match.")
counts_df <- counts_df[, meta$sample, drop = FALSE]
lib <- colSums(counts_df)
detected <- colSums(counts_df >= 10)
qc <- data.frame(sample = meta$sample, patient = meta$patient, condition = meta$condition,
                 raw_library_size = as.numeric(lib), genes_count_ge_10 = as.numeric(detected),
                 normalized_library_sum = NA_real_)
if (all(meta$sample %in% names(norm_df))) {
  nmat <- as.matrix(norm_df[, meta$sample, drop = FALSE])
  storage.mode(nmat) <- "numeric"
  qc$normalized_library_sum <- colSums(nmat)
}
write.table(qc, file.path(out, "sample_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
p <- ggplot(qc, aes(x = reorder(sample, raw_library_size), y = raw_library_size, fill = condition)) +
  geom_col() + coord_flip() + theme_bw(base_size = 10) +
  labs(title = "Raw library sizes by sample", x = "Sample", y = "Raw total counts")
ggsave(file.path(out, "library_sizes.png"), p, width = 9, height = 10, dpi = 160)
cts <- as.matrix(counts_df)
storage.mode(cts) <- "integer"
rownames(meta) <- meta$sample
dds <- DESeqDataSetFromMatrix(countData = cts, colData = meta, design = ~ patient + condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
vsd <- vst(dds, blind = FALSE)
pc <- plotPCA(vsd, intgroup = c("condition", "patient"), returnData = TRUE)
pv <- round(100 * attr(pc, "percentVar"))
p <- ggplot(pc, aes(PC1, PC2, color = condition, label = patient)) +
  geom_point(size = 3, alpha = .9) + geom_text(check_overlap = TRUE, vjust = -0.7, size = 2.7) +
  theme_bw() + labs(title = "PCA labelled by patient", x = paste0("PC1: ", pv[1], "% variance"),
                    y = paste0("PC2: ", pv[2], "% variance"))
ggsave(file.path(out, "pca_patient_labels.png"), p, width = 9, height = 7, dpi = 160)
mat <- assay(vsd)
cor_mat <- cor(mat, method = "pearson")
write.table(cor_mat, file.path(out, "sample_correlation.tsv"), sep = "\t", quote = FALSE, col.names = NA)
cor_long <- as.data.frame(as.table(cor_mat))
names(cor_long) <- c("sample_x", "sample_y", "correlation")
p <- ggplot(cor_long, aes(sample_x, sample_y, fill = correlation)) +
  geom_tile() + coord_fixed() + theme_minimal(base_size = 7) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  labs(title = "Sample correlation (VST counts)", x = NULL, y = NULL)
ggsave(file.path(out, "sample_correlation.png"), p, width = 11, height = 10, dpi = 160)
ids <- sub("\\.[0-9]+$", "", as.character(de$gene_id))
valid <- grepl("^ENSG[0-9]+$", ids)
symbols <- rep(NA_character_, length(ids))
if (any(valid)) {
  mapped <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(ids[valid]), keytype = "ENSEMBL",
                                  column = "SYMBOL", multiVals = "first")
  symbols[valid] <- unname(mapped[ids[valid]])
}
de$ensembl_id_no_version <- ids
de$gene_symbol <- symbols
de$symbol_mapping_status <- ifelse(!valid, "not_standard_Ensembl_gene_id",
                                   ifelse(is.na(symbols) | symbols == "", "unmapped", "mapped"))
write.table(de, file.path(out, "differential_expression_annotated.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
padj <- de$padj
lfc <- de$log2FoldChange
summary <- c(
  paste("Samples:", nrow(meta)),
  paste("Unique patients:", length(unique(meta$patient))),
  paste("Raw count genes:", nrow(counts_df)),
  paste("DE result rows:", nrow(de)),
  paste("FDR < 0.05:", sum(!is.na(padj) & padj < 0.05)),
  paste("FDR < 0.05 and |log2FC| >= 1:", sum(!is.na(padj) & padj < 0.05 & !is.na(lfc) & abs(lfc) >= 1)),
  paste("Upregulated (same thresholds):", sum(!is.na(padj) & padj < 0.05 & !is.na(lfc) & lfc >= 1)),
  paste("Downregulated (same thresholds):", sum(!is.na(padj) & padj < 0.05 & !is.na(lfc) & lfc <= -1)),
  paste("Gene-symbol mapping success:", sum(de$symbol_mapping_status == "mapped")),
  paste("Gene-symbol mapping missing/unmapped:", sum(de$symbol_mapping_status != "mapped")),
  "No sample is automatically excluded. Review QC plots and sample metadata before making exclusions."
)
writeLines(summary, file.path(out, "qc_summary.txt"))
message("QC outputs written to: ", out)
