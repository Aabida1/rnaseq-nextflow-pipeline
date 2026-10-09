#!/usr/bin/env python3
import argparse, csv, sys
from pathlib import Path
REQUIRED=("sample","patient","condition","fastq_1","fastq_2")
def validate_samplesheet(path, check_files=True):
    errors=[]
    try:
        with Path(path).open(newline="",encoding="utf-8-sig") as handle:
            reader=csv.DictReader(handle)
            if not reader.fieldnames: return ["Samplesheet is empty or has no header."]
            missing=[c for c in REQUIRED if c not in reader.fieldnames]
            if missing: return ["Missing required column(s): "+", ".join(missing)]
            rows=list(reader)
    except OSError as exc: return [f"Cannot read samplesheet: {exc}"]
    if not rows: errors.append("Samplesheet must contain at least one row.")
    seen=set(); pairs={}
    for line,row in enumerate(rows,2):
        for col in REQUIRED:
            if not (row.get(col) or "").strip(): errors.append(f"Line {line}: {col} is empty.")
        sample=(row.get("sample") or "").strip()
        if sample in seen: errors.append(f"Line {line}: duplicate sample ID {sample}.")
        seen.add(sample)
        if any(ch.isspace() for ch in sample): errors.append(f"Line {line}: sample ID contains whitespace.")
        if check_files:
            for col in ("fastq_1","fastq_2"):
                val=(row.get(col) or "").strip()
                if val and not Path(val).is_file(): errors.append(f"Line {line}: file does not exist: {val}")
        pairs.setdefault((row.get("patient") or "").strip(),[]).append((row.get("condition") or "").strip())
    if len({(r.get("condition") or "").strip() for r in rows})<2: errors.append("At least two conditions are required.")
    for patient,conds in pairs.items():
        if len(conds)!=2 or len(set(conds))!=2: errors.append(f"Patient {patient} must have exactly two samples in different conditions.")
    return errors
def main():
    p=argparse.ArgumentParser(); p.add_argument("samplesheet"); p.add_argument("--no-file-check",action="store_true"); a=p.parse_args()
    errors=validate_samplesheet(a.samplesheet,not a.no_file_check)
    if errors:
        print("\n".join("ERROR: "+e for e in errors),file=sys.stderr); return 1
    print("Samplesheet OK:",a.samplesheet); return 0
if __name__=="__main__": raise SystemExit(main())
