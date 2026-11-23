#!/bin/bash
# =============================================================================
# krakenuniq_pipeline.sh  –  Unified KrakenUniq pipeline (Nanopore + Illumina)
#
# Usage (called by run_pipeline.sh via sbatch, not directly):
#   Nanopore:  krakenuniq_pipeline.sh <reads.fastq.gz>
#   Illumina:  krakenuniq_pipeline.sh <R1.fastq.gz> <R2.fastq.gz>
#
# Runs with the partition/cores/RAM submitted by run_pipeline.sh (the sbatch
# flags there override the #SBATCH defaults below, which only apply to a
# direct `sbatch scripts/krakenuniq_pipeline.sh ...` call).
# --preload-size loads the DB in chunks with sufficient headroom to avoid OOM.
# Intermediate files (read_merger.pl) are written to ./tmp.
#
# All paths are relative and assume that the job runs in the package root
# directory (run_pipeline.sh calls sbatch from there without changing the
# directory first -- SLURM sets the job's working directory to this submit
# directory by default). R1/R2 are therefore also expected relative to the
# package root (e.g. "fastq/sample_R1_001.fastq.gz"), not as absolute paths.
# =============================================================================
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=zen3
#SBATCH --time=1:00:00
#SBATCH --mem=140G
#SBATCH --job-name=krakenuniq
#SBATCH --mail-type=ALL
#SBATCH --error=./log/%x_%j.err.txt
#SBATCH --output=./log/%x_%j.out.txt

set -euo pipefail

# ── Deployment-specific values ───────────────────────────────────────────────
# The only paths that point outside this package. Can be overridden via
# environment variables without editing this file, e.g.:
#   KRAKENUNIQ_DB=/other/path bash run_pipeline.sh
: "${KRAKENUNIQ_DB:=/scratch/tmp/thomachr/references/krakenuniq/microbial_db}"
: "${KRAKENUNIQ_BIN_DIR:=/scratch/tmp/thomachr/software/krakenuniq}"
: "${EXTRA_BIN_DIR:=/home/t/thomachr/bin}"

export PATH="$KRAKENUNIQ_BIN_DIR:$EXTRA_BIN_DIR:$PATH"

DATABASE="$KRAKENUNIQ_DB"
REPORT_SCRIPT="scripts/generate_report_v3.py"

# ── Input & mode detection ────────────────────────────────────────────────────
R1="$1"
R1_BASE=$(basename "$R1" .fastq.gz)

if [[ -n "${2:-}" ]]; then
    MODE="illumina"
    R2="$2"
    KRAKEN_INPUT="--paired $R1 $R2"
else
    MODE="nanopore"
    KRAKEN_INPUT="$R1"
fi

THREADS=${SLURM_CPUS_PER_TASK:-24}

echo "=== KrakenUniq Pipeline ==="
echo "Mode        : $MODE"
echo "Input       : $R1${R2:+ / $R2}"
echo "Threads     : $THREADS"
echo "Preload-size: 64G"
echo "Database    : $DATABASE"

# ── Directories ───────────────────────────────────────────────────────────────
mkdir -p log tmp output

OUT_DIR="output/${R1_BASE}"
mkdir -p "$OUT_DIR"

REPORT_TXT="$OUT_DIR/${R1_BASE}.krakenuniq.report.txt"

# ── KrakenUniq ────────────────────────────────────────────────────────────────
module purge
ml palma/2024a GCC/13.3.0 Jellyfish/2.3.1 bzip2/1.0.8

# Change into ./tmp so that read_merger.pl intermediate files end up there
# instead of in the project root. Since tmp/ is a direct subdirectory of the
# package root, the input/report paths point back there with "../" from here.
cd tmp

if [[ "$MODE" == "illumina" ]]; then
    KRAKEN_INPUT="--paired ../$R1 ../$R2"
else
    KRAKEN_INPUT="../$R1"
fi

krakenuniq \
    --preload-size 64G \
    --report-file "../$REPORT_TXT" \
    --db "$DATABASE" \
    --threads "$THREADS" \
    --output - \
    $KRAKEN_INPUT

cd ..

echo "KrakenUniq finished → $REPORT_TXT"

# ── PDF report ────────────────────────────────────────────────────────────────
if [[ -f "$REPORT_SCRIPT" ]]; then
    echo "Generating PDF report..."
    python3 "$REPORT_SCRIPT" \
        --input  "$REPORT_TXT" \
        --output "$OUT_DIR/${R1_BASE}.metagenomics_report.pdf" \
        --sample "$R1_BASE"
    echo "PDF report → $OUT_DIR/${R1_BASE}.metagenomics_report.pdf"
else
    echo "WARNING: Report script not found at $REPORT_SCRIPT – skipping PDF."
fi

echo "=== Done ==="
