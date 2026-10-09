#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: run_deseq2.R counts.txt samplesheet.csv output_dir control_condition")
counts_file <- args[[1]]; samplesheet_file <- args[[2]]; outdir <- args[[3]]; control <- args[[4]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
tab <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
if (!all(c("Geneid", "Length") %in% names(tab)) || ncol(tab) < 9) stop("Unexpected featureCounts table.")
samples <- read.csv(samplesheet_file, check.names = FALSE)
required <- c("sample", "patient", "condition", "fastq_1", "fastq_2")
if (!all(required %in% names(samples))) stop("Samplesheet lacks required columns.")
bam_cols <- names(tab)[7:ncol(tab)]
sample_ids <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(bam_cols))
if (!all(samples$sample %in% sample_ids)) stop("Could not match BAM columns to sample IDs.")
idx <- match(samples$sample, sample_ids)
counts <- as.matrix(tab[, 6 + idx, drop = FALSE]); storage.mode(counts) <- "integer"
rownames(counts) <- tab$Geneid; colnames(counts) <- samples$sample
rownames(samples) <- samples$sample
samples$condition <- relevel(factor(samples$condition), ref = control)
if (length(unique(samples$condition)) != 2) stop("Exactly two conditions are supported.")
if (any(table(samples$patient, samples$condition) != 1)) stop("Each patient must have one sample per condition.")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~ patient + condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]; dds <- DESeq(dds)
other <- setdiff(levels(samples$condition), control)[[1]]
res <- results(dds, contrast = c("condition", other, control))
out <- as.data.frame(res); out$gene_id <- rownames(out)
write.table(out, file.path(outdir, "differential_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
