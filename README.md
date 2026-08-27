# KrakenUniq-Pipeline – PALMA-Deployment

KrakenUniq-Klassifikations- und Report-Pipeline für mNGS-Analysen. Angepasst für den den PALMA-HPC-Cluster (Uni Münster). Alle Pfade **innerhalb** dieses Ordners sind relativ – der komplette Ordner kann an eine beliebige Stelle kopiert werden und funktioniert unverändert. Die einzigen Ausnahmen (Referenzdatenbank, krakenuniq-Installation, Mail-Adresse) stehen zentral in [`scripts/config.sh`](scripts/config.sh).

## Verzeichnisstruktur

```
palma_deploy/
├── run_pipeline.sh            # Einstiegspunkt: reicht sbatch-Jobs für alle fastq/-Dateien ein
├── scripts/
│   ├── config.sh               # einzige Stelle mit externen/nutzerspezifischen Pfaden
│   ├── requirements.txt        # Python-Abhängigkeiten (nur reportlab)
│   ├── krakenuniq_pipeline.sh  # sbatch-Jobskript: KrakenUniq + PDF-Report je Sample
│   ├── generate_report_v3.py   # PDF-Report-Generator
│   ├── kraken_tree.py          # gemeinsame Parsing-Logik (von generate_report_v3.py verwendet)
│   └── build_reference_db.py   # baut analysis/*.csv aus einem Korpus historischer Reports neu
├── analysis/                   # Referenz-/Hintergrunddaten für z-Score & Prävalenz (siehe unten)
│   ├── genus_reference_stats.csv
│   ├── taxon_prevalence_genus_min10kb.csv
│   └── reference_meta.json
├── fastq/                      # hier FASTQ-Dateien ablegen (siehe fastq/README.md)
└── log/, output/, tmp/         # werden beim ersten Lauf automatisch angelegt
```

## Voraussetzungen

- Zugang zu PALMA mit den Modulen `palma/2024a GCC/13.3.0 Jellyfish/2.3.1 bzip2/1.0.8`
  (werden von `scripts/krakenuniq_pipeline.sh` per `module purge && ml ...` geladen)
- Eine funktionierende krakenuniq-Installation sowie die Referenzdatenbank
  (beides extern, siehe `scripts/config.sh`)
- Python 3 mit `reportlab` (`pip install -r scripts/requirements.txt`, oder als
  Modul verfügbar machen)

## Einrichtung

1. Ordner nach PALMA kopieren (z.B. `scp -r palma_deploy/ palma:~/krakenuniq_run/`
   oder `rsync -av`).
2. `scripts/config.sh` prüfen/anpassen (Referenz-DB-Pfad, krakenuniq-
   Installationspfad, ggf. Mail-Adresse) – bei unverändertem Konto/Setup meist
   nicht nötig.
3. `pip install -r scripts/requirements.txt` (bzw. reportlab per Modul
   verfügbar machen).
4. FASTQ-Dateien nach `fastq/` legen.

## Ausführung

```bash
cd palma_deploy
bash run_pipeline.sh                              # standardmäßig ohne Mail-Benachrichtigung
bash run_pipeline.sh --mail user@uni-muenster.de   # Mail-Benachrichtigung aktivieren
```

Es wird **standardmäßig keine Mail-Adresse** gesetzt (SLURM schickt dann keine
Job-Benachrichtigungen). Adresse per `--mail`, per Umgebungsvariable
(`MAIL_USER=... bash run_pipeline.sh`) oder dauerhaft in `scripts/config.sh`
setzen.

`run_pipeline.sh` durchsucht `fastq/`, erkennt Illumina- (`*_R1_001.fastq.gz` +
`*_R2_001.fastq.gz`) und Nanopore-Dateien (alle übrigen `*.fastq.gz`) automatisch
und reicht für jedes Sample einen eigenen `sbatch`-Job ein (Partition `zen4`,
48 Cores, 140G RAM, 1h Zeitlimit). Jeder Job führt aus:

1. **KrakenUniq-Klassifikation** → `output/<sample>/<sample>.krakenuniq.report.txt`
2. **PDF-Report** (`generate_report_v3.py`) → `output/<sample>/<sample>.metagenomics_report.pdf`,
   inkl. Top-Hits/Supplement-Tabellen mit z-Scores gegenüber der mitgelieferten
   Referenzkohorte (`analysis/`)

Logs landen in `log/`, Zwischendateien in `tmp/`.

## Referenzdatenbank für z-Score & Prävalenz (`analysis/`)

**Die "Datenbank mit den 437 Referenzfällen" ist bereits vollständig
enthalten** – allerdings nicht als 437 einzelne Rohdateien, sondern bereits zu
den drei Dateien in `analysis/` aggregiert (Mittelwert/Streuung je Gattung für
den z-Score, Perzentile für die Boxplots, Prävalenz-Prozentzahlen). Das ist
alles, was `generate_report_v3.py` zur Laufzeit liest – die Stichprobengröße
(437, >=10kB-Kohorte) und das Erstellungsdatum stehen in `reference_meta.json`.
Die rohen Report-Dateien der 437 Fälle selbst sind bewusst **nicht** enthalten
(nur zur Erzeugung dieser drei Dateien nötig, nicht zur Laufzeit; siehe unten).

Diese drei Dateien werden **nicht** automatisch aktualisiert. Um sie mit neu
hinzugekommenen Reports (z.B. den auf PALMA selbst erzeugten) neu zu
berechnen:

```bash
python3 scripts/build_reference_db.py \
    --reports-dir /pfad/zu/gesammelten/rohen/reports \
    --out-dir analysis
```

(Standardwerte ohne Argumente: `./reports` bzw. `./analysis`, relativ zum
Ausführungsverzeichnis – nicht relativ zum Skript.)

## Bewusst nicht enthalten

Dieses Paket enthält nur, was für den laufenden Klassifikations-/Report-Betrieb
auf PALMA nötig ist. Nicht enthalten sind die lokalen Forschungs-/Validierungs-
Bestände des Hauptprojekts (roher Report-Korpus der 437 Referenzfälle und aller
übrigen historischen Reports, Referenzkohorten-Rohdaten, Validierungs-/
Analyseskripte, Manuskript-Tabellen) – diese bleiben im Hauptarbeitsverzeichnis
und sind für den Pipeline-Betrieb nicht erforderlich (siehe oben: die daraus
abgeleiteten `analysis/`-Dateien reichen zur Laufzeit).
