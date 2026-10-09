# Dataset card: GSE196006

- Repository: NCBI Gene Expression Omnibus (GEO)
- Title: Transcriptome analysis implicates MYC as a driver of early onset colorectal cancer
- Organism/assay: Homo sapiens, bulk RNA-seq
- Design: 21 surgically resected early-onset colorectal tumours and patient-matched adjacent colonic tissue; 42 samples total.
- GEO sample records specify hg38 and provide raw gene counts; raw sequencing reads are linked through SRA.

Links: [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196006), [raw count matrix](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE196006&format=file&file=GSE196006_raw_counts.csv.gz), [publication](https://doi.org/10.3389/fonc.2024.1365762), [SRA project PRJNA802883](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA802883).

The analysis script parses patient/condition labels from GEO matrix column names and stops unless it finds 21 complete pairs. Inspect metadata before interpreting results.
