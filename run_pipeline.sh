#!/bin/bash
# =============================================================================
# run_pipeline.sh  –  Submit KrakenUniq jobs for all FASTQ files in ./fastq/
#
# Automatically distinguishes Illumina (paired *_R1_001.fastq.gz / *_R2_001.fastq.gz)
# from Nanopore (all other *.fastq.gz) and submits scripts/krakenuniq_pipeline.sh
# as an sbatch job for each. Both modes use partition requeue / 36 cores /
# 140G RAM / 1 h time limit.
#
# Uses relative paths only and must therefore be called from this folder
# (see Usage below). User/deployment-specific values (reference DB,
# krakenuniq installation) are set in scripts/krakenuniq_pipeline.sh; the
# mail address is set here (--mail or MAIL_USER).
#
# Usage (call from this folder):
#   bash run_pipeline.sh [--mail user@uni-muenster.de]
# =============================================================================

set -euo pipefail

PIPELINE="scripts/krakenuniq_pipeline.sh"

# E-mail address for SLURM job notifications. Empty by default (no mail
# notifications) -- set via `--mail user@uni-muenster.de` or as an
# environment variable.
: "${MAIL_USER:=}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mail) MAIL_USER="$2"; shift 2 ;;
        *)      echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# Only pass --mail-user to sbatch if an address is set.
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
