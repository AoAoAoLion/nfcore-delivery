# nf-core delivery layer

nf-core gives you a pipeline. It does not give you a deliverable.

After `nextflow run nf-core/rnaseq` finishes you have alignments, counts, and a MultiQC
report — plus a samplesheet that may be quietly wrong and QC numbers nobody has looked at.
This repo is the small layer between "pipeline exited 0" and "client can read the result":

- **`make_samplesheet.py`** — build and validate an nf-core samplesheet before a six-hour run
- **`qc_gate.py`** — read the MultiQC data nf-core already wrote and apply explicit thresholds
- **`build_report.py`** — fit DESeq2 on the pipeline's own counts and write an HTML report

Everything runs on the pipeline's output. Nothing re-quantifies, nothing hides samples.

```bash
python scripts/make_samplesheet.py --fastq-dir fastq --out samplesheet.csv \
    --strandedness reverse --conditions conditions.tsv

bash scripts/run_rnaseq.sh --input samplesheet.csv --outdir results --genome GRCh38

python scripts/qc_gate.py --results results --table qc_flags.csv --strict
python scripts/build_report.py --results results --samplesheet samplesheet.csv \
    --ref control --test treated --outdir report --qc-table qc_flags.csv
```

## Why each piece exists

**Samplesheet validation.** The expensive failure is not a crash, it is a run that completes on
mismatched pairs or a duplicated sample name and produces figures nobody can defend. The
validator refuses to write a samplesheet when R1/R2 are unpaired, a sample appears twice, a
file is truncated, or the conditions file names samples with no FASTQ.

**QC gate.** A pipeline that exits 0 has not necessarily produced usable data. The gate reads
MultiQC's per-sample metrics (alignment rate from STAR or salmon, duplication from Picard or
FastQC, rRNA fraction from nf-core's custom content) and compares them to thresholds you can
see and change:

| Metric | Default | Flag when |
|---|---|---|
| alignment | 60% | below |
| duplication | 70% | above |
| reads | 5M | below |
| rRNA | 25% | above |

`--strict` exits non-zero so a wrapper script can stop before writing a report. Failing samples
are still listed in the report with their reason — they are not dropped silently.

**Report.** DESeq2 on the pipeline's gene counts with the contrast named explicitly, volcano and
PCA figures, the QC table, and a methods paragraph that matches what actually ran.

## Verified against a real nf-core run

Tested on `nf-core/rnaseq 3.22.2` (Nextflow 25.04.6, Docker) on AWS EC2, both the pipeline's own
test profile and a 4-sample human run (GM12878 vs K562, GENCODE GRCh38, salmon). Two things that
test caught and are now handled:

- nf-core 3.22 writes MultiQC data under `multiqc/<aligner>/multiqc_report_data/`, not
  `multiqc/multiqc_data/` — the loader searches instead of hard-coding a path
- MultiQC reports FastQC per read file (`SAMPLE_1`, `SAMPLE_2`) while STAR reports per sample,
  so names are normalised before merging, and per-module metric keys are used instead of a
  fuzzy "any key containing percent" match

Repeated sample rows (common when one library is split across lanes) are collapsed when their
conditions agree and rejected when they conflict.

## Requirements

`pandas`, `pydeseq2`, `matplotlib`, `seaborn`, `jinja2` (see `requirements.txt`).
`nf-core/rnaseq` itself needs Nextflow and Docker; this layer assumes the run already happened.

## What this layer does not do

- No alignment, quantification or trimming. That is nf-core's job.
- No automatic sample exclusion. Dropping a sample changes the answer, so it is a documented
  decision, not a preprocessing default.
- No batch correction unless it is in the samplesheet's condition column.
- No publication figure polish. The figures are honest defaults, not journal artwork.

## Tests

```bash
pytest tests -q
```

Covers samplesheet rejection paths, QC gating on both MultiQC layouts, report generation from a
directory shaped like nf-core output, and refusal to run with a missing condition column.

## License

MIT.
