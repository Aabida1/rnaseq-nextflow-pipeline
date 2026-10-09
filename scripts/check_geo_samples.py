#!/usr/bin/env python3
"""Fast, dependency-free validation of GSE196006 GEO matrix sample labels."""
import csv
import gzip
import re
import urllib.request

URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE196006&format=file&file=GSE196006_raw_counts.csv.gz"
PATTERN = re.compile(r"^X?([0-9]+\.[0-9]+)_[A-Za-z]([07])_G821_htseq\.out$", re.IGNORECASE)


def main():
    request = urllib.request.Request(URL, headers={"User-Agent": "rnaseq-nextflow-pipeline/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        with gzip.GzipFile(fileobj=response) as zipped:
            header = next(csv.reader((line.decode("utf-8") for line in zipped), delimiter=","))
    samples = header[1:]
    if len(samples) != 42:
        raise SystemExit(f"Expected 42 sample columns, found {len(samples)}")
    pairs = {}
    bad = []
    for sample in samples:
        match = PATTERN.fullmatch(sample)
        if not match:
            bad.append(sample)
            continue
        patient, condition_digit = match.groups()
        condition = "normal" if condition_digit == "0" else "tumour"
        pairs.setdefault(patient, []).append(condition)
    if bad:
        raise SystemExit("Unparsed GEO sample columns: " + ", ".join(bad))
    if len(pairs) != 21:
        raise SystemExit(f"Expected 21 patients, found {len(pairs)}")
    incomplete = {patient: labels for patient, labels in pairs.items()
                  if sorted(labels) != ["normal", "tumour"]}
    if incomplete:
        raise SystemExit(f"Incomplete or duplicate pairs: {incomplete}")
    print("PASS: parsed 42 GEO sample columns into 21 matched normal/tumour pairs.")
    for sample in samples:
        match = PATTERN.fullmatch(sample)
        print(f"{sample}\t{match.group(1)}\t{match.group(2)}")


if __name__ == "__main__":
    main()
