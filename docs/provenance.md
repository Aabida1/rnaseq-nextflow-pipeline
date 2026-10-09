# Provenance

The GEO analysis records accession, source URL, downloaded-file checksum, parsed sample metadata, contrast/model, and R session information.

For raw-read runs, record SRA accessions and FASTQ checksums; reference FASTA/GTF source, release and SHA-256 checksums; Nextflow version; container tag/digest; command line; Git commit; trace; report and timeline.

```bash
sha256sum samplesheet.csv reference/GRCh38.fa reference/GRCh38.gtf > input_checksums.sha256
nextflow run main.nf -profile docker --input samplesheet.csv --fasta reference/GRCh38.fa --gtf reference/GRCh38.gtf --outdir results -with-trace results/trace.tsv -with-report results/nextflow-report.html -with-timeline results/timeline.html
```

Avoid committing large data files or patient-level FASTQs. The GEO script downloads the public count matrix at runtime.
