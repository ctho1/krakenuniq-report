# fastq/

Place FASTQ files here (`.gz`-compressed). `run_pipeline.sh` detects the
sequencing type automatically from the file name:

- **Illumina** (paired-end): `*_R1_001.fastq.gz` + `*_R2_001.fastq.gz` (same
  prefix). Only R1 files with a matching R2 file are submitted.
- **Nanopore** (single-end): any other `*.fastq.gz` file.

Subdirectories are searched as well (`find` is called without a depth limit).
