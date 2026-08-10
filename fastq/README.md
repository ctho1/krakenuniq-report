# fastq/

FASTQ-Dateien hier ablegen (`.gz`-komprimiert). `run_pipeline.sh` erkennt den
Sequenziertyp automatisch am Dateinamen:

- **Illumina** (paired-end): `*_R1_001.fastq.gz` + `*_R2_001.fastq.gz` (gleicher
  Präfix). Nur R1-Dateien mit passender R2-Datei werden eingereicht.
- **Nanopore** (single-end): jede andere `*.fastq.gz`-Datei.

Unterordner werden nicht durchsucht (keine Rekursion).
