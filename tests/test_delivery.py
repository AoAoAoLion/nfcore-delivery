"""Tests for the nf-core delivery layer.

These cover the three things that break contract jobs: a wrong samplesheet, QC that was
looked at but not acted on, and a report built from the wrong table. The report test uses a
small simulated count matrix laid out exactly like nf-core/rnaseq output, so the parsing and
the statistics are exercised without needing a full pipeline run.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
PY = sys.executable


def make_fastq(directory: Path, sample: str, size: int = 4096) -> None:
    for read in ("R1", "R2"):
        path = directory / f"{sample}_{read}.fastq.gz"
        path.write_bytes(b"\x1f\x8b\x08\x00" + b"0" * size)


@pytest.fixture()
def fastq_dir(tmp_path):
    d = tmp_path / "fastq"
    d.mkdir()
    for s in ("ctrl_1", "ctrl_2", "trt_1", "trt_2"):
        make_fastq(d, s)
    return d


def run(args, **kw):
    return subprocess.run([PY, *args], capture_output=True, text=True, **kw)


def test_samplesheet_happy_path(fastq_dir, tmp_path):
    cond = tmp_path / "conditions.tsv"
    cond.write_text("sample\tcondition\nctrl_1\tcontrol\nctrl_2\tcontrol\ntrt_1\ttreated\ntrt_2\ttreated\n")
    out = tmp_path / "samplesheet.csv"
    proc = run([str(REPO / "scripts" / "make_samplesheet.py"), "--fastq-dir", str(fastq_dir),
                "--out", str(out), "--strandedness", "reverse", "--conditions", str(cond)])
    assert proc.returncode == 0, proc.stderr
    sheet = pd.read_csv(out)
    assert list(sheet.columns) == ["sample", "fastq_1", "fastq_2", "strandedness", "condition"]
    assert len(sheet) == 4
    assert (sheet["strandedness"] == "reverse").all()
    assert sheet.loc[sheet["sample"] == "trt_1", "condition"].iloc[0] == "treated"


def test_samplesheet_missing_pair(tmp_path):
    d = tmp_path / "fastq"
    d.mkdir()
    make_fastq(d, "good")
    (d / "broken_R1.fastq.gz").write_bytes(b"\x1f\x8b\x08\x00" + b"0" * 4096)
    proc = run([str(REPO / "scripts" / "make_samplesheet.py"), "--fastq-dir", str(d),
                "--out", str(tmp_path / "s.csv")])
    assert proc.returncode == 1
    assert "no R2" in proc.stderr
    assert not (tmp_path / "s.csv").exists()


def test_samplesheet_truncated_file(tmp_path):
    d = tmp_path / "fastq"
    d.mkdir()
    make_fastq(d, "good")
    for read in ("R1", "R2"):
        (d / f"tiny_{read}.fastq.gz").write_bytes(b"x")
    proc = run([str(REPO / "scripts" / "make_samplesheet.py"), "--fastq-dir", str(d),
                "--out", str(tmp_path / "s.csv")])
    assert proc.returncode == 1
    assert "truncated or empty" in proc.stderr


def test_qc_gate_flags_low_alignment(tmp_path):
    mq = tmp_path / "results" / "multiqc" / "multiqc_data"
    mq.mkdir(parents=True)
    (mq / "multiqc_general_stats.json").write_text(json.dumps({
        "STAR_mqc-generalstats-star": {
            "good_sample": {"percent_aligned": 92.5, "total_reads": 30e6},
            "bad_sample": {"percent_aligned": 31.0, "total_reads": 30e6},
        },
        "FastQC_mqc-generalstats-fastqc": {
            "good_sample": {"percent_duplicates": 22.0},
            "bad_sample": {"percent_duplicates": 88.0},
        },
    }))
    table = tmp_path / "qc_flags.csv"
    proc = run([str(REPO / "scripts" / "qc_gate.py"), "--results", str(tmp_path / "results"),
                "--table", str(table), "--strict"])
    assert proc.returncode == 1
    flags = pd.read_csv(table).set_index("sample")
    assert flags.loc["good_sample", "status"] == "PASS"
    assert flags.loc["bad_sample", "status"] == "FAIL"


def test_qc_gate_passes_clean_run(tmp_path):
    mq = tmp_path / "results" / "multiqc" / "multiqc_data"
    mq.mkdir(parents=True)
    (mq / "multiqc_general_stats.json").write_text(json.dumps({
        "STAR_mqc-generalstats-star": {"s1": {"percent_aligned": 95.0, "total_reads": 40e6}},
        "FastQC_mqc-generalstats-fastqc": {"s1": {"percent_duplicates": 15.0}},
    }))
    table = tmp_path / "qc_flags.csv"
    proc = run([str(REPO / "scripts" / "qc_gate.py"), "--results", str(tmp_path / "results"),
                "--table", str(table), "--strict"])
    assert proc.returncode == 0, proc.stderr


@pytest.fixture()
def nfcore_like(tmp_path):
    """A directory tree shaped like nf-core/rnaseq output, with simulated counts."""
    rng = np.random.default_rng(3)
    samples = [f"ctrl_{i}" for i in (1, 2, 3)] + [f"trt_{i}" for i in (1, 2, 3)]
    genes = [f"ENSG{i:08d}" for i in range(1, 501)]
    counts = rng.negative_binomial(5, 0.02, size=(500, 6)).astype(int)
    counts[0:60, 3:] *= 4          # spiked signal so the test has something to find
    df = pd.DataFrame(counts, index=genes, columns=samples)
    df.insert(0, "gene_name", [f"GENE{i}" for i in range(1, 501)])

    results = tmp_path / "results" / "star_salmon"
    results.mkdir(parents=True)
    df.to_csv(results / "salmon.merged.gene_counts.tsv", sep="\t")

    sheet = tmp_path / "samplesheet.csv"
    rows = []
    for s in samples:
        rows.append({"sample": s, "fastq_1": f"{s}_R1.fastq.gz", "fastq_2": f"{s}_R2.fastq.gz",
                     "strandedness": "reverse", "condition": "control" if s.startswith("ctrl") else "treated"})
    pd.DataFrame(rows).to_csv(sheet, index=False)

    mq = tmp_path / "results" / "multiqc" / "multiqc_data"
    mq.mkdir(parents=True)
    (mq / "multiqc_general_stats.json").write_text(json.dumps({
        "STAR_mqc-generalstats-star": {s: {"percent_aligned": 93.0, "total_reads": 25e6} for s in samples},
    }))
    return tmp_path


def test_report_builds_from_nfcore_output(nfcore_like):
    qc = nfcore_like / "qc_flags.csv"
    run([str(REPO / "scripts" / "qc_gate.py"), "--results", str(nfcore_like / "results"),
         "--table", str(qc)])
    out = nfcore_like / "report"
    proc = run([str(REPO / "scripts" / "build_report.py"),
                "--results", str(nfcore_like / "results"),
                "--samplesheet", str(nfcore_like / "samplesheet.csv"),
                "--outdir", str(out), "--qc-table", str(qc), "--cpus", "2"])
    assert proc.returncode == 0, proc.stderr + proc.stdout
    html = (out / "report.html").read_text()
    assert "Differential expression" in html
    assert "nf-core/rnaseq" in html
    res = pd.read_csv(out / "tables" / "deseq2_results.csv", index_col=0)
    assert "padj" in res.columns and len(res) > 100
    assert (out / "figures" / "volcano.png").exists()
    summary = json.loads((out / "summary.json").read_text())
    assert summary["n_samples"] == 6
    assert summary["n_sig"] > 0, summary


def test_report_refuses_missing_condition_column(nfcore_like):
    sheet = pd.read_csv(nfcore_like / "samplesheet.csv").drop(columns=["condition"])
    bad = nfcore_like / "bad.csv"
    sheet.to_csv(bad, index=False)
    proc = run([str(REPO / "scripts" / "build_report.py"),
                "--results", str(nfcore_like / "results"),
                "--samplesheet", str(bad), "--outdir", str(nfcore_like / "r2")])
    assert proc.returncode != 0
    assert "condition" in (proc.stderr + proc.stdout)
