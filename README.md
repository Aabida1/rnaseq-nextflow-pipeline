# Early-onset colorectal cancer RNA-seq

Reproducible RNA-seq portfolio project using real public human bulk RNA-seq data: **NCBI GEO GSE196006**, with 21 early-onset colorectal cancer patients and paired tumour/adjacent-normal tissue (42 samples). Sources: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196006) · [publication](https://doi.org/10.3389/fonc.2024.1365762).

## Real-data analysis
Install dependencies with `Rscript env/install_bioc.R`, then run:

```bash
Rscript scripts/analyze_gse196006.R --outdir results
```

The script downloads the public raw-count matrix, checks the expected 21 matched pairs, runs paired DESeq2 (`~ patient + condition`), and generates differential-expression tables, PCA/volcano/MA plots, GO enrichment, provenance and R session information.

## Raw-read workflow
Requires Nextflow 23.10+, Docker, compatible FASTA/GTF and paired-end FASTQs. The sample sheet columns are `sample,patient,condition,fastq_1,fastq_2`. Validate with `python scripts/validate_samplesheet.py samplesheet.csv`, build with `docker build -t rnaseq-pipeline:1.0.0 .`, then run:

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --fasta reference/GRCh38.fa --gtf reference/GRCh38.gtf --outdir results --control_condition normal
```

### FASTQ quality-control stages
1. **Raw-read QC:** FastQC runs on both original FASTQ files before trimming. Reports are published under `results/qc/raw_fastqc/`.
2. **Read trimming and filtering:** fastp removes adapters and applies its default quality/length filters; its HTML and JSON reports are saved under `results/trimmed/`.
3. **Post-trimming QC:** FastQC runs on the trimmed FASTQs. Reports are published under `results/qc/trimmed_fastqc/`.
4. **Combined report:** MultiQC aggregates raw FastQC, fastp and post-trimming FastQC reports into `results/multiqc/multiqc_report.html`.

Review the MultiQC report before interpreting downstream alignment and counts. QC reports are diagnostic; the pipeline does not automatically discard samples based on a single metric. Confirm SRA read layout before using GSE196006 raw reads; the GEO count-matrix analysis does not assume read layout.

## Tools
FastQC, fastp, STAR, featureCounts, MultiQC, Nextflow DSL2, Docker, DESeq2, clusterProfiler, Python validation and GitHub Actions.

## Validation
GitHub Actions is configured for Python tests and a Nextflow stub-mode workflow check. Stub mode checks workflow wiring only; it does not validate real alignment or biological results. Do not claim the real-data analysis completed until the GEO download, R script, and CI have actually run successfully.

Expected outputs: `sample_metadata.tsv`, `differential_expression.tsv`, `normalized_counts.tsv`, `pca.png`, `volcano.png`, `ma_plot.png`, `go_bp_enrichment.tsv`, `input_provenance.txt`, `analysis_notes.txt`, and `sessionInfo.txt`.

See `docs/dataset.md`, `docs/analysis_plan.md`, and `docs/provenance.md`. Adjacent-normal tissue can show field effects; enrichment is exploratory and association is not causation.
