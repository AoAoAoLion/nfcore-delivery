#!/usr/bin/env python3
"""Build and validate an nf-core/rnaseq samplesheet from a directory of FASTQ files.

The most expensive failure mode in a contract RNA-seq job is launching a 6-hour
nf-core run on a samplesheet that is silently wrong: mismatched pairs, a sample
named twice, a strandedness typo, or an empty file that only fails at the
quantification step. This script catches those in two seconds.

Expected input layout (any nesting depth):

    fastq/
      ctrl_1_R1.fastq.gz  ctrl_1_R2.fastq.gz
      ctrl_2_R1.fastq.gz  ctrl_2_R2.fastq.gz

Usage:
    python scripts/make_samplesheet.py --fastq-dir fastq --out samplesheet.csv \
        --strandedness reverse
    python scripts/make_samplesheet.py --fastq-dir fastq --out samplesheet.csv \
        --conditions conditions.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

R1_PATTERNS = [r"_R1$", r"[._]1$"]
R2_PATTERNS = [r"_R2$", r"[._]2$"]
VALID_STRANDEDNESS = {"forward", "reverse", "unstranded", "auto"}
VALID_SUFFIX = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
MIN_BYTES = 1024


def classify(path: Path) -> tuple[str | None, str]:
    """Return (read, sample_stem) where read is 'R1', 'R2' or None."""
    name = path.name
    for suffix in VALID_SUFFIX:
        if name.endswith(suffix):
            stem = name[: -len(suffix)]
            break
    else:
        return None, name

    for pat in R2_PATTERNS:
        if re.search(pat, stem):
            return "R2", re.sub(pat, "", stem)
    for pat in R1_PATTERNS:
        if re.search(pat, stem):
            return "R1", re.sub(pat, "", stem)
    return None, stem


def collect(fastq_dir: Path) -> tuple[dict[str, dict[str, Path]], list[str]]:
    """Group files by sample. Returns (samples, problems)."""
    samples: dict[str, dict[str, Path]] = defaultdict(dict)
    problems: list[str] = []

    for path in sorted(fastq_dir.rglob("*")):
        if not path.is_file():
            continue
        if not path.name.endswith(VALID_SUFFIX):
            if path.name.lower().endswith((".bam", ".cram", ".sam")):
                problems.append(f"{path.name}: alignment file in the FASTQ directory")
            continue
        read, stem = classify(path)
        if read is None:
            problems.append(f"{path.name}: cannot tell R1 from R2 in the filename")
            continue
        if read in samples[stem]:
            problems.append(f"{stem}: two files resolve to {read} "
                            f"({samples[stem][read].name} and {path.name})")
            continue
        samples[stem][read] = path

    for stem, reads in samples.items():
        if "R1" not in reads:
            problems.append(f"{stem}: no R1 file")
        if "R2" not in reads:
            problems.append(f"{stem}: no R2 file (single-end runs need --single-end)")
        for read, path in reads.items():
            if path.stat().st_size < MIN_BYTES:
                problems.append(f"{path.name}: {path.stat().st_size} bytes, looks truncated or empty")
    return dict(samples), problems


def read_conditions(path: Path) -> dict[str, str]:
    """Two-column TSV/CSV: sample<TAB>condition."""
    out: dict[str, str] = {}
    text = path.read_text()
    delim = "\t" if "\t" in text.splitlines()[0] else ","
    for i, line in enumerate(text.splitlines()):
        if not line.strip() or (i == 0 and line.lower().startswith("sample")):
            continue
        parts = [p.strip() for p in line.split(delim)]
        if len(parts) < 2:
            raise SystemExit(f"{path}: line {i + 1} needs two columns (sample, condition)")
        out[parts[0]] = parts[1]
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fastq-dir", required=True)
    ap.add_argument("--out", default="samplesheet.csv")
    ap.add_argument("--strandedness", default="auto",
                    help="forward | reverse | unstranded | auto (auto runs the pipeline's inference)")
    ap.add_argument("--conditions", default=None,
                    help="optional TSV of sample<TAB>condition for the downstream report")
    ap.add_argument("--single-end", action="store_true")
    ap.add_argument("--allow-absolute", action="store_true",
                    help="write absolute paths instead of paths relative to the samplesheet")
    args = ap.parse_args()

    if args.strandedness not in VALID_STRANDEDNESS:
        raise SystemExit(f"strandedness must be one of {sorted(VALID_STRANDEDNESS)}")

    fastq_dir = Path(args.fastq_dir)
    if not fastq_dir.is_dir():
        raise SystemExit(f"not a directory: {fastq_dir}")

    samples, problems = collect(fastq_dir)

    if args.single_end:
        problems = [p for p in problems if "single-end runs need" not in p]

    if not samples:
        raise SystemExit(f"no FASTQ files found under {fastq_dir}")

    conditions = read_conditions(Path(args.conditions)) if args.conditions else {}
    unknown = set(conditions) - set(samples)
    if unknown:
        problems.append(f"conditions file names samples with no FASTQ: {sorted(unknown)}")

    out_path = Path(args.out)
    out_dir = out_path.parent if str(out_path.parent) else Path(".")
    rows = []
    for sample in sorted(samples):
        reads = samples[sample]
        fq1 = reads["R1"]
        fq2 = reads.get("R2")
        fmt = (lambda p: str(p.resolve()) if args.allow_absolute
               else str(p.resolve().relative_to(out_dir.resolve()))
               if _is_relative(p.resolve(), out_dir.resolve()) else str(p.resolve()))
        row = {
            "sample": sample,
            "fastq_1": fmt(fq1),
            "fastq_2": fmt(fq2) if fq2 else "",
            "strandedness": args.strandedness,
        }
        if conditions:
            row["condition"] = conditions.get(sample, "NA")
        rows.append(row)

    if problems:
        print("samplesheet NOT written, fix these first:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        raise SystemExit(1)

    fields = list(rows[0].keys())
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {out_path}: {len(rows)} samples, strandedness={args.strandedness}")
    for row in rows:
        cond = f"  [{row['condition']}]" if "condition" in row else ""
        print(f"  {row['sample']}{cond}")
    if conditions:
        missing = set(samples) - set(conditions)
        if missing:
            print(f"note: no condition given for {sorted(missing)}", file=sys.stderr)


def _is_relative(path: Path, base: Path) -> bool:
    try:
        path.relative_to(base)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    main()
