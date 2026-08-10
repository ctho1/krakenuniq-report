#!/bin/bash
# =============================================================================
# krakenuniq_pipeline.sh  –  Unified KrakenUniq pipeline (Nanopore + Illumina)
#
# Usage (called by run_pipeline.sh via sbatch, not directly):
#   Nanopore:  krakenuniq_pipeline.sh <reads.fastq.gz>
#   Illumina:  krakenuniq_pipeline.sh <R1.fastq.gz> <R2.fastq.gz>
#
# Runs on zen4 with cores/RAM as submitted by run_pipeline.sh (sbatch flags
# there override the #SBATCH defaults below, which only apply for a direct
# `sbatch scripts/krakenuniq_pipeline.sh ...` call).
# --preload-size loads the DB in chunks with sufficient headroom to avoid OOM.
# Intermediate files (read_merger.pl) are written to ./tmp.
#
# Alle paketinternen Pfade (dieses Skript, generate_report_v3.py, analysis/)
# werden relativ zum Paket-Wurzelverzeichnis aufgelöst -- das Paket kann also
# komplett an eine andere Stelle/ein anderes Konto kopiert werden. Die
# einzigen externen, deployment-spezifischen Pfade (Referenz-DB, krakenuniq-
# Installation) stehen in config.sh (im selben scripts/-Ordner).
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

# ── Paket-Wurzelverzeichnis + Konfiguration ──────────────────────────────────
# Robust aus dem eigenen Skriptpfad abgeleitet (scripts/ -> eine Ebene hoch),
# nicht aus pwd, damit das Skript unabhängig vom SLURM-Arbeitsverzeichnis
# funktioniert (run_pipeline.sh setzt zusätzlich --chdir auf denselben Pfad).
PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=config.sh
source "$PACKAGE_ROOT/scripts/config.sh"

export PATH="$KRAKENUNIQ_BIN_DIR:$EXTRA_BIN_DIR:$PATH"

DATABASE="$KRAKENUNIQ_DB"
REPORT_SCRIPT="$PACKAGE_ROOT/scripts/generate_report_v3.py"

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

PROJECT_ROOT="$PACKAGE_ROOT"

echo "=== KrakenUniq Pipeline ==="
echo "Mode        : $MODE"
echo "Input       : $R1${R2:+ / $R2}"
echo "Threads     : $THREADS"
echo "Preload-size: 64G"
echo "Database    : $DATABASE"
echo "Project root: $PROJECT_ROOT"

# ── Directories ───────────────────────────────────────────────────────────────
mkdir -p "$PROJECT_ROOT"/log \
         "$PROJECT_ROOT"/tmp \
         "$PROJECT_ROOT"/output

OUT_DIR="$PROJECT_ROOT/output/${R1_BASE}"
mkdir -p "$OUT_DIR"

REPORT_TXT="$OUT_DIR/${R1_BASE}.krakenuniq.report.txt"

# ── KrakenUniq ────────────────────────────────────────────────────────────────
module purge
ml palma/2024a GCC/13.3.0 Jellyfish/2.3.1 bzip2/1.0.8

# cd into ./tmp so read_merger.pl intermediate files land there, not project root
cd "$PROJECT_ROOT/tmp"

krakenuniq \
    --preload-size 64G \
    --report-file "$REPORT_TXT" \
    --db "$DATABASE" \
    --threads "$THREADS" \
    --output - \
    $KRAKEN_INPUT

cd "$PROJECT_ROOT"

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
