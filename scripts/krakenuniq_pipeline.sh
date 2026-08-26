#!/bin/bash
# =============================================================================
# krakenuniq_pipeline.sh  –  Unified KrakenUniq pipeline (Nanopore + Illumina)
#
# Usage (called by run_pipeline.sh via sbatch, not directly):
#   Nanopore:  krakenuniq_pipeline.sh <reads.fastq.gz>
#   Illumina:  krakenuniq_pipeline.sh <R1.fastq.gz> <R2.fastq.gz>
#
# Runs on zen3 with cores/RAM as submitted by run_pipeline.sh (sbatch flags
# there override the #SBATCH defaults below, which only apply for a direct
# `sbatch scripts/krakenuniq_pipeline.sh ...` call).
# --preload-size loads the DB in chunks with sufficient headroom to avoid OOM.
# Intermediate files (read_merger.pl) are written to ./tmp.
#
# Alle Pfade sind relativ und gehen davon aus, dass der Job im Paket-
# Wurzelverzeichnis läuft (run_pipeline.sh ruft sbatch von dort aus auf, ohne
# vorher das Verzeichnis zu wechseln -- SLURM setzt das Arbeitsverzeichnis des
# Jobs standardmäßig auf dieses Submit-Verzeichnis). R1/R2 werden daher
# ebenfalls relativ zum Paket-Root erwartet (z.B. "fastq/probe_R1_001.fastq.gz"),
# nicht absolut.
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

# ── Deployment-spezifische Werte ─────────────────────────────────────────────
# Vorher in scripts/config.sh, jetzt direkt hier. Weiterhin per Umgebungs-
# variable überschreibbar, ohne diese Datei zu ändern, z.B.:
#   KRAKENUNIQ_DB=/anderer/pfad bash run_pipeline.sh
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

# In ./tmp wechseln, damit read_merger.pl-Zwischendateien dort landen statt im
# Projekt-Root. Da tmp/ eine direkte Unterebene des Paket-Roots ist, zeigen
# die Pfade für Input/Report ab hier mit "../" zurück dorthin.
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
