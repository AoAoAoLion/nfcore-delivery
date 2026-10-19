#!/usr/bin/env bash
# Launch nf-core/rnaseq with the settings we actually want for contract work.
#
# Everything a client cares about is pinned and visible: pipeline version, profile,
# genome, strandedness and the resume behaviour. No hidden defaults.
#
# Usage:
#   bash scripts/run_rnaseq.sh --input samplesheet.csv --outdir results \
#        --genome GRCh38 --profile docker --strandedness reverse
set -euo pipefail

VERSION="${NF_VERSION:-3.22.2}"
INPUT="samplesheet.csv"
OUTDIR="results"
GENOME="GRCh38"
PROFILE="docker"
STRANDEDNESS="reverse"
MAX_CPUS="${MAX_CPUS:-8}"
MAX_MEMORY="${MAX_MEMORY:-32.GB}"
EXTRA="${EXTRA:-}"

while [ $# -gt 0 ]; do
  case "$1" in
    --input) INPUT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --genome) GENOME="$2"; shift 2 ;;
    --profile) PROFILE="$2"; shift 2 ;;
    --strandedness) STRANDEDNESS="$2"; shift 2 ;;
    --version) VERSION="$2"; shift 2 ;;
    *) EXTRA="$EXTRA $1"; shift ;;
  esac
done

if ! command -v nextflow >/dev/null 2>&1; then
  echo "nextflow not found. Install: curl -s https://get.nextflow.io | bash && mv nextflow ~/bin/" >&2
  exit 1
fi
if [ ! -s "$INPUT" ]; then
  echo "samplesheet not found: $INPUT (build it with scripts/make_samplesheet.py)" >&2
  exit 1
fi

echo "[run] nf-core/rnaseq ${VERSION} | genome ${GENOME} | strandedness ${STRANDEDNESS} | profile ${PROFILE}"
echo "[run] input ${INPUT} -> ${OUTDIR}"

nextflow run nf-core/rnaseq \
  -r "${VERSION}" \
  -profile "${PROFILE}" \
  --input "$INPUT" \
  --outdir "$OUTDIR" \
  --genome "$GENOME" \
  --aligner star_salmon \
  --remove_ribo_rna \
  --skip_bbsplit true \
  --skip_deseq2_qc false \
  --max_cpus "$MAX_CPUS" \
  --max_memory "$MAX_MEMORY" \
  -resume \
  $EXTRA

echo "[run] done. Next steps:"
echo "  python scripts/qc_gate.py --results ${OUTDIR} --table qc_flags.csv --strict"
echo "  python scripts/build_report.py --results ${OUTDIR} --samplesheet ${INPUT} --outdir report"
