# KrakenUniq Pipeline – PALMA Deployment

KrakenUniq classification and report pipeline for mNGS analyses, adapted for the PALMA HPC cluster (University of Münster). All paths **within** this folder are relative – the complete folder can be copied to any location and works unchanged. The only exceptions (reference database, krakenuniq installation) are set at the top of [`scripts/krakenuniq_pipeline.sh`](scripts/krakenuniq_pipeline.sh).

## Directory structure

```
krakenuniq-report/
├── run_pipeline.sh            # entry point: submits sbatch jobs for all files in fastq/
├── scripts/
│   ├── requirements.txt        # Python dependencies (reportlab only)
│   ├── krakenuniq_pipeline.sh  # sbatch job script: KrakenUniq + PDF report per sample
│   ├── generate_report_v3.py   # PDF report generator
│   ├── kraken_tree.py          # shared parsing logic (used by generate_report_v3.py)
│   └── build_reference_db.py   # rebuilds analysis/*.csv from a corpus of historical reports
├── analysis/                   # reference/background data for z-score & prevalence (see below)
│   ├── genus_reference_stats.csv
│   ├── taxon_prevalence_genus_min10kb.csv
│   └── reference_meta.json
├── fastq/                      # place FASTQ files here (see fastq/README.md)
└── log/, output/, tmp/         # created automatically on the first run
```

## Requirements

- Access to PALMA with the modules `palma/2024a GCC/13.3.0 Jellyfish/2.3.1 bzip2/1.0.8`
  (loaded by `scripts/krakenuniq_pipeline.sh` via `module purge && ml ...`)
- A working krakenuniq installation and the reference database
  (both external, see the top of `scripts/krakenuniq_pipeline.sh`)
- Python 3 with `reportlab` (`pip install -r scripts/requirements.txt`, or make
  it available as a module)

## Setup

1. Check/adjust the paths at the top of `scripts/krakenuniq_pipeline.sh`
   (`KRAKENUNIQ_DB`, `KRAKENUNIQ_BIN_DIR`, `EXTRA_BIN_DIR`). Not necessary
   if account and setup are unchanged. All three can also be overridden via
   environment variables, e.g. `KRAKENUNIQ_DB=/other/path bash run_pipeline.sh`.
2. `pip install -r scripts/requirements.txt` (or make reportlab available as a
   module).
3. Copy the FASTQ files into `fastq/`.

## Usage

```bash
cd krakenuniq-report
bash run_pipeline.sh                               # no mail notifications by default
bash run_pipeline.sh --mail user@uni-muenster.de   # enable mail notifications
```

`run_pipeline.sh` searches `fastq/`, automatically distinguishes Illumina
(`*_R1_001.fastq.gz` + `*_R2_001.fastq.gz`) from Nanopore files (all other
`*.fastq.gz`) and submits a separate `sbatch` job for each sample (partition
`requeue`, 36 cores, 140G RAM, 1 h time limit). Each job runs:

1. **KrakenUniq classification** → `output/<sample>/<sample>.krakenuniq.report.txt`
2. **PDF report** (`generate_report_v3.py`) → `output/<sample>/<sample>.metagenomics_report.pdf`,
   including top-hits/supplement tables with z-scores relative to the bundled
   reference cohort (`analysis/`)

Logs are written to `log/`, intermediate files to `tmp/`.

A PDF report can also be generated locally from an existing KrakenUniq report:

```bash
python3 scripts/generate_report_v3.py --input <sample>.krakenuniq.report.txt \
    --output <sample>.metagenomics_report.pdf --sample "<sample name>"
python3 scripts/generate_report_v3.py --output demo.pdf   # built-in demo data
```

## Reference database for z-score & prevalence (`analysis/`)

**The "database of 437 reference cases" is already fully included** – not as
437 individual raw files, but already aggregated into the three files in
`analysis/` (mean/spread per genus for the z-score, percentiles for the
boxplots, prevalence percentages). This is all that `generate_report_v3.py`
reads at runtime – the sample size (437, >=10kB cohort) and the creation date
are stored in `reference_meta.json`. The raw report files of the 437 cases
themselves are deliberately **not** included (only needed to generate these
three files, not at runtime; see below).

These three files are **not** updated automatically. To recompute them with
newly added reports (e.g. those generated on PALMA itself):

```bash
python3 scripts/build_reference_db.py \
    --reports-dir /path/to/collected/raw/reports \
    --out-dir analysis
```

(Defaults without arguments: `reports/` and `analysis/` in the package root,
i.e. relative to the `scripts/` folder – not to the current working
directory. Besides the three runtime files, the script also writes further
statistics files to `--out-dir` that are not needed for the report.)
