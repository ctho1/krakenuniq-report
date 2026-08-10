#!/bin/bash
# =============================================================================
# config.sh – Einzige Stelle mit deployment-/nutzerspezifischen Werten.
#
# Alle anderen Pfade in diesem Paket sind relativ zum Paket-Wurzelverzeichnis
# (also portabel: der ganze Ordner kann an eine beliebige Stelle auf PALMA
# oder in ein anderes Nutzerkonto kopiert werden). Nur die hier gelisteten
# Werte verweisen auf Ressourcen AUSSERHALB dieses Pakets (Referenzdatenbank,
# krakenuniq-Installation, Mail-Adresse) und müssen bei einem neuen
# Nutzer/Konto ggf. angepasst werden.
#
# Wird von run_pipeline.sh und scripts/krakenuniq_pipeline.sh eingelesen.
# Werte können auch per Umgebungsvariable überschrieben werden, ohne diese
# Datei zu ändern, z.B.:
#   KRAKENUNIQ_DB=/anderer/pfad bash run_pipeline.sh
# =============================================================================

# E-Mail für SLURM-Job-Benachrichtigungen. Standardmäßig leer (keine Mail-
# Benachrichtigung) -- per `--mail user@uni-muenster.de` an run_pipeline.sh
# oder als Umgebungsvariable setzen.
: "${MAIL_USER:=}"

# KrakenUniq-Referenzdatenbank (groß, extern -- nicht Teil dieses Pakets)
: "${KRAKENUNIQ_DB:=/scratch/tmp/thomachr/references/krakenuniq/microbial_db}"

# Zusätzliche PATH-Einträge für die krakenuniq-Installation bzw. persönliche
# Hilfsskripte (extern -- nicht Teil dieses Pakets)
: "${KRAKENUNIQ_BIN_DIR:=/scratch/tmp/thomachr/software/krakenuniq}"
: "${EXTRA_BIN_DIR:=/home/t/thomachr/bin}"

export MAIL_USER KRAKENUNIQ_DB KRAKENUNIQ_BIN_DIR EXTRA_BIN_DIR
