#!/usr/bin/env python3
"""Gate nf-core/rnaseq output on QC metrics before anything is reported to a client.

Why this exists: a pipeline that exits 0 has not necessarily produced usable data. Low
alignment rate, heavy duplication or a sample that is 40% rRNA will still generate figures,
and those figures are what gets a deliverable rejected. This script reads the MultiQC data
that nf-core already writes, applies explicit thresholds, and says which samples need a
conversation before the report is written.

Metric sources are matched per module (STAR, salmon, Picard, FastQC, nf-core custom content),
because MultiQC stores each tool's numbers under its own keys and a naive "any key containing
percent" match silently picks the wrong column.

Usage:
    python scripts/qc_gate.py --results results --table qc_flags.csv [--strict]

Exit code is 1 with --strict when any sample fails a hard threshold.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

# module -> (key, scale) in priority order
ALIGNED_SOURCES = [
    ("star", "uniquely_mapped_percent", 1.0),
    ("salmon", "percent_mapped", 1.0),
    ("rseqc", "unique_percent", 1.0),
    ("samtools_3", "mapped_passed_pct", 1.0),
    ("custom_content", "% Aligned", 1.0),
]
DUP_SOURCES = [
    ("picard", "PERCENT_DUPLICATION", 100.0),          # Picard reports a fraction
    ("fastqc_trimmed", "percent_duplicates", 1.0),
    ("fastqc_raw", "percent_duplicates", 1.0),
]
READS_SOURCES = [
    ("star", "total_reads", 1.0),
    ("rseqc", "total_records", 1.0),
    ("samtools", "reads_mapped", 1.0),
]
RRNA_SOURCES = [
    ("custom_content", "percent_rRNA", 1.0),
    ("sortmerna", "percent", 1.0),
]

DEFAULTS = {"min_aligned": 60.0, "max_duplication": 70.0, "min_reads_m": 5.0, "max_rrna": 25.0}


def normalize(sample: str) -> str:
    """FastQC reports one row per read file (SAMPLE_1, SAMPLE_2); merge them."""
    return re.sub(r"_[12]$", "", sample)


def load_modules(results: Path) -> dict[str, dict[str, dict]]:
    """Return {module: {sample: metrics}} from whatever MultiQC wrote."""
    candidates = sorted(results.glob("**/multiqc_general_stats.json")) + \
                 sorted(results.glob("**/multiqc_data.json"))
    if not candidates:
        raise SystemExit(f"no MultiQC JSON under {results} — did the pipeline run with --skip_multiqc?")

    data = json.loads(candidates[0].read_text())
    modules: dict[str, dict[str, dict]] = {}

    if "report_general_stats_data" in data:                 # multiqc_data.json
        sections = data["report_general_stats_data"]
        items = sections.items() if isinstance(sections, dict) else enumerate(sections)
        for module, samples in items:
            if isinstance(samples, dict):
                modules[str(module)] = samples
    else:                                                   # multiqc_general_stats.json
        for module, samples in data.items():
            if isinstance(samples, dict) and any(isinstance(v, dict) for v in samples.values()):
                modules[module] = samples

    # merge read-pair rows within each module
    merged: dict[str, dict[str, dict]] = {}
    for module, samples in modules.items():
        bucket: dict[str, dict] = {}
        for sample, metrics in samples.items():
            if not isinstance(metrics, dict):
                continue
            bucket.setdefault(normalize(sample), {}).update(metrics)
        merged[module] = bucket
    return merged


def lookup(modules: dict, sources, sample: str) -> float | None:
    for module, key, scale in sources:
        value = modules.get(module, {}).get(sample, {}).get(key)
        if isinstance(value, (int, float)):
            return float(value) * scale
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default="results")
    ap.add_argument("--table", default="qc_flags.csv")
    ap.add_argument("--min-aligned", type=float, default=DEFAULTS["min_aligned"])
    ap.add_argument("--max-duplication", type=float, default=DEFAULTS["max_duplication"])
    ap.add_argument("--min-reads-m", type=float, default=DEFAULTS["min_reads_m"])
    ap.add_argument("--max-rrna", type=float, default=DEFAULTS["max_rrna"])
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    modules = load_modules(Path(args.results))
    samples = sorted({s for mod in modules.values() for s in mod})
    # keep only samples that appear in an alignment or quantification module
    quantified = set(modules.get("star", {})) | set(modules.get("salmon", {})) | \
                 set(modules.get("rseqc", {}))
    samples = [s for s in samples if s in quantified] or samples
    if not samples:
        raise SystemExit("MultiQC data contains no per-sample entries")

    rows, failures = [], []
    for sample in samples:
        aligned = lookup(modules, ALIGNED_SOURCES, sample)
        dup = lookup(modules, DUP_SOURCES, sample)
        reads = lookup(modules, READS_SOURCES, sample)
        rrna = lookup(modules, RRNA_SOURCES, sample)
        reads_m = None if reads is None else reads / 1e6   # STAR total_reads is a count

        flags = []
        if aligned is not None and aligned < args.min_aligned:
            flags.append(f"aligned {aligned:.1f}% < {args.min_aligned}%")
        if dup is not None and dup > args.max_duplication:
            flags.append(f"duplication {dup:.1f}% > {args.max_duplication}%")
        if reads_m is not None and reads_m < args.min_reads_m:
            flags.append(f"reads {reads_m:.1f}M < {args.min_reads_m}M")
        if rrna is not None and rrna > args.max_rrna:
            flags.append(f"rRNA {rrna:.1f}% > {args.max_rrna}%")

        if flags:
            failures.append((sample, flags))
        rows.append({
            "sample": sample,
            "aligned_pct": "" if aligned is None else round(aligned, 2),
            "duplication_pct": "" if dup is None else round(dup, 2),
            "reads_millions": "" if reads_m is None else round(reads_m, 3),
            "rrna_pct": "" if rrna is None else round(rrna, 2),
            "status": "FAIL" if flags else "PASS",
            "flags": "; ".join(flags),
        })

    with open(args.table, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"{'sample':<26}{'aligned%':>10}{'dup%':>9}{'reads(M)':>11}{'rRNA%':>8}  status")
    for r in rows:
        print(f"{r['sample']:<26}{str(r['aligned_pct']):>10}{str(r['duplication_pct']):>9}"
              f"{str(r['reads_millions']):>11}{str(r['rrna_pct']):>8}  {r['status']}")

    if failures:
        print("\nsamples needing a decision before reporting:", file=sys.stderr)
        for sample, flags in failures:
            print(f"  {sample}: {'; '.join(flags)}", file=sys.stderr)
        if args.strict:
            raise SystemExit(1)
    print(f"\nwrote {args.table}")


if __name__ == "__main__":
    main()
