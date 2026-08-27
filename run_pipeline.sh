#!/bin/bash
# =============================================================================
# run_pipeline.sh  –  Submit KrakenUniq jobs for all FASTQ files in ./fastq/
#
# Automatically distinguishes Illumina (paired *_R1_001.fastq.gz / *_R2_001.fastq.gz)
# from Nanopore (all other *.fastq.gz) and submits scripts/krakenuniq_pipeline.sh
# as an sbatch job for each. Both modes use zen3 / 48 cores / 140G RAM.
#
# Verwendet ausschließlich relative Pfade und muss daher aus diesem Ordner
# heraus aufgerufen werden (siehe Usage unten). Nutzer-/deployment-spezifische
# Werte (Referenz-DB, krakenuniq-Installation, Mail-Adresse) stehen direkt in
# diesem Skript bzw. in scripts/krakenuniq_pipeline.sh.
#
# Usage (aus diesem Ordner heraus aufrufen):
#   bash run_pipeline.sh [--mail user@uni-muenster.de]
# =============================================================================

set -euo pipefail

PIPELINE="scripts/krakenuniq_pipeline.sh"

# E-Mail für SLURM-Job-Benachrichtigungen. Standardmäßig leer (keine Mail-
# Benachrichtigung) -- per `--mail user@uni-muenster.de` oder als
# Umgebungsvariable setzen.
: "${MAIL_USER:=}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mail) MAIL_USER="$2"; shift 2 ;;
        *)      echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# --mail-user nur an sbatch übergeben, wenn eine Adresse gesetzt ist.
MAIL_ARGS=()
if [[ -n "$MAIL_USER" ]]; then
    MAIL_ARGS=(--mail-user="$MAIL_USER")
fi

mkdir -p log output fastq

submitted=0
skipped=0

# ── Illumina: paired R1 / R2 ─────────────────────────────────────────────────
echo "=== Illumina samples ==="

while IFS= read -r -d '' R1; do
    R2="${R1/R1_001/R2_001}"

    if [[ ! -f "$R2" ]]; then
        echo "  SKIP $R1  (no matching R2 found: $R2)"
        (( skipped++ )) || true
        continue
    fi

    BASE=$(basename "$R1" .fastq.gz)
    echo "  Submitting $BASE"

    sbatch \
        "${MAIL_ARGS[@]+"${MAIL_ARGS[@]}"}" \
        --job-name="$BASE" \
        --partition=requeue \
        --cpus-per-task=36 \
        --mem=140G \
        --time=1:00:00 \
        "$PIPELINE" "$R1" "$R2"

    (( submitted++ )) || true

done < <(find fastq -type f -name "*_R1_001*.fastq.gz" -print0 | sort -z)

# ── Nanopore: single-end ──────────────────────────────────────────────────────
echo ""
echo "=== Nanopore samples ==="

while IFS= read -r -d '' FASTQ; do
    BASE=$(basename "$FASTQ" .fastq.gz)
    echo "  Submitting $BASE"

    sbatch \
        "${MAIL_ARGS[@]+"${MAIL_ARGS[@]}"}" \
        --job-name="$BASE" \
        --partition=requeue \
        --cpus-per-task=36 \
        --mem=140G \
        --time=1:00:00 \
        "$PIPELINE" "$FASTQ"

    (( submitted++ )) || true

done < <(find fastq -type f -name "*.fastq.gz" \
         ! -name "*_R1_001*.fastq.gz" \
         ! -name "*_R2_001*.fastq.gz" \
         -print0 | sort -z)

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Submitted: $submitted job(s)  |  Skipped: $skipped ==="
