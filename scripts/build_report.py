#!/usr/bin/env python3
"""Turn nf-core/rnaseq output into a client-facing report.

nf-core gives you alignments, quantification and a MultiQC page. What it does not give you
is a short document a biologist can read: which samples passed QC, what changed between
groups, which figures go in the paper, and a methods paragraph that matches what was run.

This script reads the pipeline's own output (no re-quantification), fits DESeq2 on the
gene counts, and writes an HTML report that cites the MultiQC report for the raw QC.

Usage:
    python scripts/build_report.py --results results --samplesheet samplesheet.csv \
        --condition-col condition --ref control --test treated --outdir report
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from jinja2 import Template
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

sns.set_theme(style="whitegrid", context="notebook")

# nf-core/rnaseq writes counts here depending on the aligner
COUNT_CANDIDATES = [
    ("star_salmon", "salmon.merged.gene_counts.tsv"),
    ("salmon", "salmon.merged.gene_counts.tsv"),
    ("star_rsem", "rsem.merged.gene_counts.tsv"),
    ("featurecounts", "merged_gene_counts.tsv"),
]


def log(msg: str) -> None:
    print(f"[report] {msg}", flush=True)


def find_counts(results: Path) -> tuple[Path, str]:
    for folder, name in COUNT_CANDIDATES:
        path = results / folder / name
        if path.exists():
            return path, folder
    hint = [str(p.relative_to(results)) for p in results.rglob("*gene_counts*.tsv")][:5]
    raise SystemExit(f"no gene count table found under {results}. Looked for: "
                     f"{[f'{f}/{n}' for f, n in COUNT_CANDIDATES]}. Found: {hint}")


def load_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    # nf-core writes gene_id as index and gene_name as first column
    if "gene_name" in df.columns:
        df = df.drop(columns=["gene_name"])
    df = df.drop(columns=[c for c in ("gene_id",) if c in df.columns], errors="ignore")
    df = df.loc[:, [c for c in df.columns if not c.startswith("tpm")]]
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    return df


def load_metadata(samplesheet: Path, condition_col: str) -> pd.DataFrame:
    sheet = pd.read_csv(samplesheet)
    if "sample" not in sheet.columns:
        raise SystemExit(f"{samplesheet} has no 'sample' column")
    if condition_col not in sheet.columns:
        raise SystemExit(f"{samplesheet} has no '{condition_col}' column; "
                         f"add it (sample,fastq_1,fastq_2,strandedness,{condition_col})")
    # nf-core collapses repeated sample names (e.g. multiple lanes of the same library)
    # into one column. A sample sheet that keeps both rows would silently mis-assign
    # conditions, so collapse identical conditions and refuse conflicting ones.
    dup = sheet["sample"][sheet["sample"].duplicated()].unique().tolist()
    if dup:
        grouped = sheet.groupby("sample")[condition_col].nunique()
        conflicting = sorted(grouped[grouped > 1].index.tolist())
        if conflicting:
            raise SystemExit(f"samplesheet has duplicate rows for {conflicting} with different "
                             f"{condition_col} values; keep one row per sample")
        log(f"collapsed duplicate rows for {dup} (same {condition_col}, treated as one library)")
    meta = sheet.drop_duplicates("sample").set_index("sample")[[condition_col]]
    return meta


def qc_table(path: Path) -> list[dict] | None:
    if not path.exists():
        return None
    return pd.read_csv(path).to_dict("records")


def plots(res: pd.DataFrame, counts: pd.DataFrame, meta: pd.DataFrame,
          condition_col: str, outdir: Path) -> pd.DataFrame:
    outdir.mkdir(parents=True, exist_ok=True)
    sig = res[(res["padj"] < 0.05) & (res["log2FoldChange"].abs() >= 1)]

    fig, ax = plt.subplots(figsize=(9, 6.5))
    ns = res.drop(sig.index)
    ax.scatter(ns["log2FoldChange"], -np.log10(ns["padj"]), s=6, c="#B0B7C3", label="not significant")
    up, down = sig[sig["log2FoldChange"] > 0], sig[sig["log2FoldChange"] < 0]
    ax.scatter(up["log2FoldChange"], -np.log10(up["padj"]), s=9, c="#C44E52", label=f"up ({len(up)})")
    ax.scatter(down["log2FoldChange"], -np.log10(down["padj"]), s=9, c="#4C72B0", label=f"down ({len(down)})")
    for x in (-1, 1):
        ax.axvline(x, ls="--", c="k", lw=1)
    ax.axhline(-np.log10(0.05), ls="--", c="k", lw=1)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("Volcano plot")
    ax.legend(frameon=False, fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "volcano.png", dpi=150)
    plt.close(fig)

    logcpm = np.log2(counts.div(counts.sum(axis=0), axis=1) * 1e6 + 1)
    x = logcpm.sub(logcpm.mean(axis=1), axis=0)
    u, s, vt = np.linalg.svd(x.values, full_matrices=False)
    var = (s ** 2) / np.sum(s ** 2) * 100
    pcs = vt[:2].T
    fig, ax = plt.subplots(figsize=(7, 5.5))
    groups = meta.loc[counts.columns, condition_col].astype(str)
    palette = dict(zip(groups.unique(), sns.color_palette("deep", groups.nunique())))
    for grp in groups.unique():
        m = (groups == grp).values
        ax.scatter(pcs[m, 0], pcs[m, 1], s=70, label=grp, color=palette[grp], edgecolor="k")
    dist = np.hypot(pcs[:, 0], pcs[:, 1])
    for grp in groups.unique():
        idx = np.where((groups == grp).values)[0]
        for i in idx[np.argsort(dist[idx])[-2:]]:
            ax.annotate(counts.columns[i], (pcs[i, 0], pcs[i, 1]), fontsize=8,
                        xytext=(5, 5), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.set_title("PCA on log2 CPM")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "pca.png", dpi=150)
    plt.close(fig)
    return sig


TEMPLATE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>RNA-seq report</title>
<style>
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;max-width:960px;margin:40px auto;padding:0 20px;color:#1b2430;line-height:1.6}
h1{border-bottom:3px solid #2ec4b6;padding-bottom:10px}
h2{margin-top:34px;border-bottom:1px solid #e3e8ef;padding-bottom:6px}
table{border-collapse:collapse;width:100%;font-size:14px;margin:12px 0}
th,td{border:1px solid #e3e8ef;padding:8px 10px;text-align:left}
th{background:#f6f8fb}
img{max-width:100%;border:1px solid #e3e8ef;border-radius:6px;margin:8px 0}
code{background:#f2f4f7;padding:2px 5px;border-radius:3px;font-size:13px}
.meta{color:#5b6b7c;font-size:14px}
.fail{color:#b42318;font-weight:600}
</style></head><body>
<h1>RNA-seq differential expression report</h1>
<p class="meta">Contrast {{ test }} vs {{ ref }} &middot; {{ n_samples }} samples &middot;
quantified with nf-core/rnaseq ({{ aligner }}) &middot; {{ date }}</p>

<h2>1. Samples and QC</h2>
{% if qc %}
<table><tr><th>Sample</th><th>Aligned %</th><th>Duplication %</th><th>Reads (M)</th><th>rRNA %</th><th>Status</th></tr>
{% for r in qc %}<tr><td>{{ r.sample }}</td><td>{{ r.aligned_pct }}</td><td>{{ r.duplication_pct }}</td><td>{{ r.reads_millions }}</td><td>{{ r.rrna_pct }}</td><td class="{{ 'fail' if r.status == 'FAIL' else '' }}">{{ r.status }}</td></tr>{% endfor %}
</table>
<p>Full per-tool QC is in <code>results/multiqc/multiqc_report.html</code>. Samples flagged FAIL
are reported as they are, with the reason; they were not silently dropped.</p>
{% else %}
<p>No QC table was found next to this report. Run <code>scripts/qc_gate.py</code> first; the
raw MultiQC report is at <code>results/multiqc/multiqc_report.html</code>.</p>
{% endif %}

<h2>2. Differential expression</h2>
<table><tr><th>Metric</th><th>Value</th></tr>
<tr><td>Genes tested</td><td>{{ n_tested }}</td></tr>
<tr><td>padj &lt; 0.05</td><td>{{ n_padj }}</td></tr>
<tr><td>padj &lt; 0.05 and |log2FC| &ge; 1</td><td>{{ n_sig }}</td></tr>
<tr><td>Up / down</td><td>{{ n_up }} / {{ n_down }}</td></tr>
</table>
<img src="figures/volcano.png" alt="volcano">
<img src="figures/pca.png" alt="PCA">
<p>Model <code>~{{ condition }}</code>, DESeq2 with median-of-ratios normalisation
(PyDESeq2 {{ version }}), Wald test, Benjamini-Hochberg adjustment.</p>

<h2>3. Methods paragraph</h2>
<p>Reads were processed with nf-core/rnaseq ({{ aligner }}), including adapter trimming, rRNA
removal and quantification against the supplied reference. Gene-level counts were modelled with
DESeq2 using the design <code>~{{ condition }}</code>; the contrast {{ test }} versus {{ ref }} was
tested with the Wald test and p-values were adjusted with the Benjamini-Hochberg procedure.
Genes with adjusted p &lt; 0.05 were considered differentially expressed; figure highlighting
additionally requires |log2 fold change| &ge; 1. Per-sample QC metrics are reported above and in
the accompanying MultiQC report.</p>

<h2>4. Files</h2>
<ul>
<li><code>report/tables/deseq2_results.csv</code> &mdash; full result table</li>
<li><code>report/figures/</code> &mdash; volcano and PCA</li>
<li><code>results/multiqc/multiqc_report.html</code> &mdash; per-tool QC</li>
<li><code>results/</code> &mdash; nf-core/rnaseq output, including the Nextflow execution report</li>
</ul>
</body></html>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default="results")
    ap.add_argument("--samplesheet", default="samplesheet.csv")
    ap.add_argument("--condition-col", default="condition")
    ap.add_argument("--ref", default="control")
    ap.add_argument("--test", default="treated")
    ap.add_argument("--outdir", default="report")
    ap.add_argument("--qc-table", default="qc_flags.csv")
    ap.add_argument("--cpus", type=int, default=2)
    args = ap.parse_args()

    results = Path(args.results)
    outdir = Path(args.outdir)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    (outdir / "tables").mkdir(parents=True, exist_ok=True)

    counts_path, aligner = find_counts(results)
    counts = load_counts(counts_path)
    meta = load_metadata(Path(args.samplesheet), args.condition_col)
    missing = set(counts.columns) - set(meta.index)
    if missing:
        raise SystemExit(f"samples in counts but not in samplesheet: {sorted(missing)}")
    meta = meta.loc[counts.columns]

    log(f"counts: {counts.shape[0]} genes x {counts.shape[1]} samples ({counts_path.name}, aligner={aligner})")
    counts = counts.loc[counts.sum(axis=1) >= 10]

    dds = DeseqDataSet(counts=counts.T, metadata=meta, design=f"~{args.condition_col}",
                       refit_cooks=True, n_cpus=args.cpus, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=[args.condition_col, args.test, args.ref],
                       n_cpus=args.cpus, quiet=True)
    stats.summary()
    res = stats.results_df.copy()
    res.index.name = "gene"
    res.sort_values("padj", na_position="last").to_csv(outdir / "tables" / "deseq2_results.csv")

    sig = plots(res, counts, meta, args.condition_col, outdir / "figures")
    qc = qc_table(Path(args.qc_table))

    import pydeseq2

    ctx = dict(
        test=args.test, ref=args.ref, n_samples=int(counts.shape[1]), aligner=aligner,
        n_tested=int(res["padj"].notna().sum()), n_padj=int((res["padj"] < 0.05).sum()),
        n_sig=int(len(sig)), n_up=int((sig["log2FoldChange"] > 0).sum()),
        n_down=int((sig["log2FoldChange"] < 0).sum()), condition=args.condition_col,
        version=pydeseq2.__version__, qc=qc,
        date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"),
    )
    (outdir / "report.html").write_text(Template(TEMPLATE).render(**ctx))
    (outdir / "summary.json").write_text(json.dumps(
        {k: ctx[k] for k in ("test", "ref", "n_samples", "n_tested", "n_padj", "n_sig", "n_up", "n_down")},
        indent=2))
    log(f"done: {ctx['n_padj']} genes padj<0.05 -> {outdir/'report.html'}")


if __name__ == "__main__":
    main()
