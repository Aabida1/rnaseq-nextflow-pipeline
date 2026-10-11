#!/usr/bin/env python3
"""Generate a deterministic, tiny paired-end RNA-seq dataset for real-tool smoke tests."""
import gzip
import random
import sys
from pathlib import Path

out = Path(sys.argv[1] if len(sys.argv) > 1 else "smoke-data")
out.mkdir(parents=True, exist_ok=True)
rng = random.Random(196006)
alphabet = "ACGT"
complement = str.maketrans("ACGT", "TGCA")
genome = "".join(rng.choice(alphabet) for _ in range(1_000_000))
(out / "reference.fa").write_text(">chrSmoke\n" + "\n".join(genome[i:i+80] for i in range(0, len(genome), 80)) + "\n")

genes = []
for i in range(40):
    start0 = 10_000 + i * 20_000
    end0 = start0 + 5_000
    gene_id = f"ENSGSMOKE{i+1:011d}"
    genes.append((gene_id, start0, end0))
with (out / "genes.gtf").open("w") as handle:
    for gene_id, start0, end0 in genes:
        handle.write(f'chrSmoke\tsmoke\texon\t{start0+1}\t{end0}\t.\t+\t.\tgene_id "{gene_id}"; transcript_id "TX_{gene_id}";\n')

samples = [
    ("P01_N", "P01", "normal"),
    ("P01_T", "P01", "tumour"),
    ("P02_N", "P02", "normal"),
    ("P02_T", "P02", "tumour"),
]
with (out / "samplesheet.csv").open("w") as sheet:
    sheet.write("sample,patient,condition,fastq_1,fastq_2\n")
    for sample, patient, condition in samples:
        r1_path = out / f"{sample}_R1.fastq.gz"
        r2_path = out / f"{sample}_R2.fastq.gz"
        with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
            read_number = 0
            for gi, (gene_id, start0, end0) in enumerate(genes):
                base = 24 + (gi * 13) % 37
                condition_effect = 2.0 if condition == "tumour" and gi < 12 else (0.5 if condition == "tumour" and 12 <= gi < 24 else 1.0)
                patient_effect = 1.15 if patient == "P02" and gi % 3 == 0 else (0.85 if patient == "P02" and gi % 3 == 1 else 1.0)
                sample_jitter = rng.uniform(0.65, 1.35)
                n_pairs = max(12, round(base * condition_effect * patient_effect * sample_jitter))
                for j in range(n_pairs):
                    max_start = end0 - start0 - 225
                    offset = rng.randrange(0, max_start)
                    pos = start0 + offset
                    read1_seq = genome[pos:pos+75]
                    read2_seq = genome[pos+125:pos+200].translate(complement)[::-1]
                    read_number += 1
                    name = f"@{sample}:{read_number}:{gene_id}:{pos}"
                    r1.write(f"{name}/1\n{read1_seq}\n+\n{'I'*75}\n")
                    r2.write(f"{name}/2\n{read2_seq}\n+\n{'I'*75}\n")
        sheet.write(f"{sample},{patient},{condition},{r1_path.resolve()},{r2_path.resolve()}\n")
print(f"Generated {len(samples)} paired-end samples and {len(genes)} genes in {out.resolve()}")
