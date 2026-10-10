#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
outdir <- "results"
if (length(args) >= 2 && args[[1]] == "--outdir") outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
has_annotation <- requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
                  requireNamespace("AnnotationDbi", quietly = TRUE)
has_clusterprofiler <- requireNamespace("clusterProfiler", quietly = TRUE)
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE196006&format=file&file=GSE196006_raw_counts.csv.gz"
dest <- file.path(outdir, "GSE196006_raw_counts.csv.gz")
if (!file.exists(dest) || file.info(dest)$size < 1000) download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
if (!file.exists(dest) || file.info(dest)$size < 1000) stop("GEO count matrix download failed or file is unexpectedly small.")
writeLines(c("GEO accession: GSE196006", paste("Source URL:", url), paste("File:", basename(dest)), paste("MD5:", unname(tools::md5sum(dest)))), file.path(outdir, "input_provenance.txt"))
raw <- read.csv(gzfile(dest), check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(raw) < 3 || nrow(raw) < 100) stop("Downloaded matrix has an unexpected shape.")
rownames(raw) <- as.character(raw[[1]])
raw[[1]] <- NULL
counts <- as.matrix(raw)
storage.mode(counts) <- "numeric"
if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) stop("Count matrix contains invalid/non-integer values.")
counts <- round(counts[!grepl("^__", rownames(counts)), , drop = FALSE])
samples <- colnames(counts)
sample_match <- regexec("^X?([0-9]+\\.[0-9]+)_[A-Za-z]([07])_G821_htseq\\.out$", samples, ignore.case = TRUE)
sample_parts <- regmatches(samples, sample_match)
if (any(lengths(sample_parts) != 3L)) {
  writeLines(samples, file.path(outdir, "unparsed_sample_columns.txt"))
  stop("Could not parse patient/condition labels from GEO matrix column names; see unparsed_sample_columns.txt.")
}
patient <- vapply(sample_parts, `[[`, character(1), 2)
condition_digit <- vapply(sample_parts, `[[`, character(1), 3)
condition <- ifelse(condition_digit == "0", "normal", "tumour")
metadata <- data.frame(sample = samples, patient = patient, condition = factor(condition, levels = c("normal", "tumour")))
if (nrow(metadata) != 42 || length(unique(metadata$patient)) != 21) stop("Expected 42 samples from 21 patients.")
pair_table <- table(metadata$patient, metadata$condition)
if (ncol(pair_table) != 2 || any(pair_table != 1)) stop("Each patient must have exactly one normal and one tumour sample.")
rownames(metadata) <- metadata$sample
write.table(metadata, file.path(outdir, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(counts, file.path(outdir, "raw_counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ patient + condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "tumour", "normal"), alpha = 0.05)
res_df <- as.data.frame(res); res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id"))]
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
write.table(res_df, file.path(outdir, "differential_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
norm <- counts(dds, normalized = TRUE)
write.table(data.frame(gene_id = rownames(norm), norm, check.names = FALSE), file.path(outdir, "normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Sample-level QC: library size and number of detected genes.
qc <- data.frame(
  sample = colnames(counts),
  patient = metadata[colnames(counts), "patient"],
  condition = metadata[colnames(counts), "condition"],
  raw_library_size = colSums(counts),
  detected_genes = colSums(counts > 0),
  normalized_library_size = colSums(norm),
  stringsAsFactors = FALSE
)
qc$matched_pair_correlation <- NA_real_
for (id in unique(qc$patient)) {
  pair_samples <- qc$sample[qc$patient == id]
  if (length(pair_samples) == 2L) {
    pair_cor <- suppressWarnings(cor(
      assay(vst(dds[, pair_samples], blind = TRUE)),
      method = "pearson"
    )[1, 2])
    qc$matched_pair_correlation[qc$sample %in% pair_samples] <- pair_cor
  }
}
write.table(qc, file.path(outdir, "sample_qc_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
sample_cor <- cor(vst_mat, method = "pearson", use = "pairwise.complete.obs")
write.table(data.frame(sample = rownames(sample_cor), sample_cor, check.names = FALSE),
            file.path(outdir, "sample_correlation_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cor_long <- as.data.frame(as.table(sample_cor), stringsAsFactors = FALSE)
names(cor_long) <- c("sample_x", "sample_y", "correlation")
cor_long$sample_x <- factor(cor_long$sample_x, levels = colnames(sample_cor))
cor_long$sample_y <- factor(cor_long$sample_y, levels = rev(colnames(sample_cor)))
p_cor <- ggplot(cor_long, aes(sample_x, sample_y, fill = correlation)) +
  geom_tile() + scale_fill_gradient2(limits = c(-1, 1), midpoint = 0) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 5),
        axis.title = element_blank(), panel.grid = element_blank()) +
  labs(fill = "Pearson r", title = "Sample correlation (VST counts)")
ggsave(file.path(outdir, "sample_correlation_heatmap.png"), p_cor,
       width = 12, height = 11, dpi = 180)

pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"))
pca$sample <- if ("name" %in% names(pca)) as.character(pca$name) else rownames(pca)
p <- ggplot(pca, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 3) + geom_text(check_overlap = TRUE, nudge_y = 1.2, size = 2.4, show.legend = FALSE) +
  xlab(paste0("PC1: ", pv[1], "% variance")) + ylab(paste0("PC2: ", pv[2], "% variance")) +
  theme_bw() + labs(title = "PCA of variance-stabilized counts (sample labels)")
ggsave(file.path(outdir, "pca.png"), p, width = 10, height = 7, dpi = 180)

pair_ids <- unique(metadata$patient)
pair_qc <- do.call(rbind, lapply(pair_ids, function(id) {
  pair_samples <- rownames(metadata)[metadata$patient == id]
  if (length(pair_samples) != 2L) return(NULL)
  data.frame(patient = id,
             normal_sample = pair_samples[metadata[pair_samples, "condition"] == "normal"],
             tumour_sample = pair_samples[metadata[pair_samples, "condition"] == "tumour"],
             vst_pearson_correlation = unname(cor(vst_mat[, pair_samples[1]], vst_mat[, pair_samples[2]], method = "pearson")),
             stringsAsFactors = FALSE)
}))
write.table(pair_qc, file.path(outdir, "matched_pair_correlations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
plot_df <- res_df
plot_df$neglog10padj <- -log10(pmax(plot_df$padj, 1e-300))
plot_df$significance <- ifelse(!is.na(plot_df$padj) & plot_df$padj < 0.05 & abs(plot_df$log2FoldChange) >= 1, "FDR < 0.05 and |log2FC| >= 1", "Other")
p <- ggplot(plot_df, aes(log2FoldChange, neglog10padj, color = significance)) + geom_point(alpha = 0.55, size = 1) + theme_bw() + labs(title = "Tumour versus adjacent normal", x = "log2 fold change", y = "-log10 adjusted p-value")
ggsave(file.path(outdir, "volcano.png"), p, width = 7, height = 5, dpi = 160)
png(file.path(outdir, "ma_plot.png"), width = 1200, height = 900, res = 150); plotMA(res, ylim = c(-5, 5)); dev.off()
sig <- sub("\\.[0-9]+$", "", res_df$gene_id[!is.na(res_df$padj) & res_df$padj < 0.05])
go_result <- data.frame()
if (has_clusterprofiler && has_annotation && length(sig) > 0 && all(grepl("^ENSG[0-9]+$", sig))) {
  # Attach the annotation package so its OrgDb object is available by name.
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
  })
  go_result <- tryCatch({
    mapped <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(sig), keytype = "ENSEMBL", column = "ENTREZID", multiVals = "first")
    mapped <- unique(na.omit(unname(mapped)))
    if (length(mapped) > 0) {
      ego <- clusterProfiler::enrichGO(gene = mapped, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = TRUE)
      as.data.frame(ego)
    } else {
      data.frame()
    }
  }, error = function(e) {
    message("Optional GO enrichment skipped: ", conditionMessage(e))
    data.frame()
  })
} else {
  message("Optional GO enrichment skipped because its packages or eligible gene IDs are unavailable.")
}
write.table(go_result, file.path(outdir, "go_bp_enrichment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c("Dataset: NCBI GEO GSE196006", "Contrast: tumour versus matched adjacent normal", "Model: ~ patient + condition", "FDR threshold: 0.05", "GO enrichment uses human Ensembl IDs; empty results can mean no significant genes or unmapped IDs.", "Adjacent normal tissue may show field effects."), file.path(outdir, "analysis_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("Analysis completed; inspect outputs and QC before interpreting biology.")
