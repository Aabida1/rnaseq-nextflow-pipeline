# Analysis plan

Question: which genes and pathways differ between early-onset colorectal tumour tissue and patient-matched adjacent-normal tissue?

Use DESeq2 with paired design `~ patient + condition`, reference `normal`, contrast `tumour - normal`. DESeq2 receives raw integer counts; Benjamini-Hochberg adjusted p-values control false discovery rate.

Outputs: differential-expression table, normalized counts for visualization, PCA, volcano and MA plots, GO Biological Process enrichment where mapping succeeds, parsed sample metadata, checksum and R session information.

Safeguards: remove HTSeq special rows beginning with `__`; require unique samples and exactly one normal and tumour sample per patient; never use normalized counts as DESeq2 input. Adjacent normal tissue may show field effects and is not equivalent to healthy volunteer tissue. Enrichment is exploratory.
