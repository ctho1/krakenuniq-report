#!/usr/bin/env python3
"""
Baut die Referenz-Datenbank fuer die KrakenUniq-Filterstrategie.

Liest alle KrakenUniq-Reports in ./reports ein und erzeugt:
  1. sample_metrics.csv           - QC-Metriken je Sample (alle Samples)
  2. taxon_prevalence_all.csv     - Taxon-Praevalenz (Species/Strain-Ebene)
                                     ueber ALLE Samples
  3. taxon_prevalence_min10kb.csv - wie 2., nur Samples mit Report-Groesse
                                     >= 10 KB
  4. taxon_reference_stats.csv    - Verteilungskennzahlen (Perzentile) je
                                     Taxon fuer die >=10kb-Kohorte
  5. taxon_prevalence_genus_all.csv - Praevalenz auf GENUS-Ebene (ueber den
                                     Taxonomiebaum aggregiert) ueber ALLE
                                     Samples. Ein Sample zaehlt einmal pro
                                     Genus, unabhaengig davon wie viele
                                     Species/Strains dieser Gattung darin
                                     vorkommen. Wird von generate_report_v3.py
                                     fuer die Haupt-/Supplement-Priorisierung
                                     verwendet, da Genus-Ebene robuster gegen
                                     "Avalanchen" einzelner Referenzstaemme
                                     ist (z.B. viele Bradyrhizobium sp. XY).
  6. genus_reference_stats.csv     - wie taxon_reference_stats.csv, aber je
                                     Genus statt je Species/Strain: pro Sample
                                     werden Reads/k-mers aller Species/Strains
                                     einer Gattung summiert (Coverage als Max,
                                     dup Read-gewichtet), dann ueber die
                                     >=10kb-Kohorte aggregiert. Enthaelt auch
                                     adj_kmers_log_mean/_sd (coverage-
                                     adjustierte k-mers, kraken_tree.
                                     adjusted_kmers(): k-mers je Spezies durch
                                     deren dup-Wert gedaempft, ausser bei
                                     Viren) -- Basis des Genus-Ebene z-Scores
                                     in generate_report_v3.py.

Homo sapiens (taxID 9606) wird von vornherein ausgeschlossen.
"""

import argparse
import csv
import json
import math
import os
import re
import glob
import sys
from collections import defaultdict
from datetime import date

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kraken_tree

REPORTS_DIR = os.path.join(os.path.dirname(__file__), "..", "reports")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "analysis")
os.makedirs(OUT_DIR, exist_ok=True)

SIZE_THRESHOLD_BYTES = 10 * 1024  # 10 KB

# Taxonomic ranks that represent a concrete organism-level call
# (mirrors RELEVANT_RANKS in generate_report_v3.py)
RELEVANT_RANKS = {
    "species", "subspecies", "strain", "isolate",
    "no rank", "varietas", "serotype",
    "species group", "species subgroup",
}

# Container / aggregate taxIDs that must never be counted as a "hit"
# even if some row of that name happened to carry taxReads > 0.
EXCLUDED_TAXIDS = {"0", "1", "131567", "9606"}
EXCLUDED_NAMES = {"Homo sapiens", "unclassified", "root", "cellular organisms"}

PLATE_RE = re.compile(r"SIL-\d+")


def parse_report(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()

    total_classified = total_unclassified = total_human = 0
    pct_classified = pct_unclassified = pct_human = 0.0
    db = None
    date = None

    for line in text.split("\n"):
        if not line.startswith("#"):
            continue
        if db is None:
            m = re.search(r"DB:(\S+)", line)
            if m:
                db = m.group(1)
        if date is None:
            m = re.search(r"DATE:(\S+)", line)
            if m:
                date = m.group(1)

    entries = kraken_tree.parse_kraken_rows(text)

    for e in entries:
        if e["taxID"] == "0":
            total_unclassified = e["reads"]
            pct_unclassified = e["pct"]
        if e["taxID"] == "1":
            total_classified = e["reads"]
            pct_classified = e["pct"]
        if e["taxID"] == "9606":
            total_human = e["reads"]
            pct_human = e["pct"]

    hits = [
        e for e in entries
        if e["taxReads"] > 0
        and e["rank"] in RELEVANT_RANKS
        and e["taxID"] not in EXCLUDED_TAXIDS
        and e["name"] not in EXCLUDED_NAMES
    ]

    total_all = total_classified + total_unclassified

    return dict(
        db=db, date=date,
        total_all=total_all, total_classified=total_classified,
        total_unclassified=total_unclassified, total_human=total_human,
        pct_classified=pct_classified, pct_unclassified=pct_unclassified,
        pct_human=pct_human, hits=hits,
    )


def percentile(sorted_vals, p):
    """Nearest-rank percentile, 0<=p<=100, sorted_vals ascending, non-empty."""
    if not sorted_vals:
        return None
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    k = (len(sorted_vals) - 1) * (p / 100)
    f = int(k)
    c = min(f + 1, len(sorted_vals) - 1)
    if f == c:
        return sorted_vals[f]
    d0 = sorted_vals[f] * (c - k)
    d1 = sorted_vals[c] * (k - f)
    return d0 + d1


def main():
    global REPORTS_DIR, OUT_DIR
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reports-dir", default=REPORTS_DIR,
                     help="Verzeichnis mit *.krakenuniq.report.txt (Default: ./reports)")
    ap.add_argument("--out-dir", default=OUT_DIR,
                     help="Ausgabeverzeichnis fuer die CSV/JSON-Referenzdateien (Default: ./analysis)")
    args = ap.parse_args()
    REPORTS_DIR = args.reports_dir
    OUT_DIR = args.out_dir
    os.makedirs(OUT_DIR, exist_ok=True)

    paths = sorted(glob.glob(os.path.join(REPORTS_DIR, "*.krakenuniq.report.txt")))
    print(f"Gefundene Reports: {len(paths)}")

    sample_rows = []
    # taxon_key -> {"all": [sample_records...], "big": [...]}
    taxon_hits_all = defaultdict(list)
    taxon_hits_big = defaultdict(list)
    # genus_key -> set of sample names in which the genus occurred at least once
    genus_samples_all = defaultdict(set)
    genus_samples_big = defaultdict(set)
    genus_names = {}  # genus_taxid -> genus_name
    # genus_key -> [per-sample aggregated dict, ...] (>=10kb cohort only) -
    # basis for genus_reference_stats.csv (z-score/boxplot reference at
    # Genus-Ebene). Ein Sample liefert genau EINEN aggregierten Datenpunkt
    # je Genus (Reads/k-mers summiert über alle Species/Strains dieser
    # Gattung im Sample, Coverage als Max, dup Read-gewichtet).
    genus_hits_big = defaultdict(list)

    n_all = 0
    n_big = 0

    for path in paths:
        fname = os.path.basename(path)
        sample = fname.replace(".krakenuniq.report.txt", "")
        size_bytes = os.path.getsize(path)
        is_big = size_bytes >= SIZE_THRESHOLD_BYTES

        n_all += 1
        if is_big:
            n_big += 1

        parsed = parse_report(path)
        hits = parsed["hits"]

        plate_m = PLATE_RE.search(sample)
        plate_id = plate_m.group(0) if plate_m else ""

        total_nonhuman_hit_reads = sum(h["taxReads"] for h in hits)
        n_distinct_taxa = len(hits)
        if hits:
            top_hit = max(hits, key=lambda h: h["taxReads"])
            top_taxon_name = top_hit["name"]
            top_taxon_reads = top_hit["taxReads"]
        else:
            top_taxon_name = ""
            top_taxon_reads = 0

        sample_rows.append(dict(
            sample=sample,
            plate_id=plate_id,
            report_date=parsed["date"] or "",
            file_size_bytes=size_bytes,
            pass_min10kb=is_big,
            total_reads=parsed["total_all"],
            classified_reads=parsed["total_classified"],
            pct_classified=round(parsed["pct_classified"], 3),
            unclassified_reads=parsed["total_unclassified"],
            pct_unclassified=round(parsed["pct_unclassified"], 3),
            human_reads=parsed["total_human"],
            pct_human=round(parsed["pct_human"], 3),
            nonhuman_classified_reads=max(
                parsed["total_classified"] - parsed["total_human"], 0
            ),
            n_distinct_taxa_hits=n_distinct_taxa,
            total_nonhuman_hit_reads=total_nonhuman_hit_reads,
            top_taxon_name=top_taxon_name,
            top_taxon_reads=top_taxon_reads,
        ))

        sample_genus_agg = defaultdict(
            lambda: dict(reads=0, kmers=0, cov=0.0, dup_wsum=0.0, adj_kmers=0.0, kingdom=None))

        for h in hits:
            key = (h["taxID"], h["name"], h["rank"])
            record = dict(
                sample=sample,
                taxReads=h["taxReads"],
                kmers=h["kmers"],
                dup=h["dup"],
                cov=h["cov"],
            )
            taxon_hits_all[key].append(record)
            if is_big:
                taxon_hits_big[key].append(record)

            genus_key = h["genus_taxid"]
            genus_names[genus_key] = h["genus_name"]
            genus_samples_all[genus_key].add(sample)
            if is_big:
                genus_samples_big[genus_key].add(sample)
                agg = sample_genus_agg[genus_key]
                # Kingdom is inferred once per genus/sample, from the first
                # member encountered -- mirrors group_by_genus() in
                # generate_report_v3.py (kingdom=members[0]["kingdom"]), so
                # the reference distribution and live scoring use the same
                # adjusted_kmers() basis for a given genus.
                if agg["kingdom"] is None:
                    agg["kingdom"] = kraken_tree.infer_kingdom(h["name"])
                agg["reads"] += h["taxReads"]
                agg["kmers"] += h["kmers"]
                agg["cov"] = max(agg["cov"], h["cov"])
                agg["dup_wsum"] += h["dup"] * h["taxReads"]
                agg["adj_kmers"] += kraken_tree.adjusted_kmers(h["kmers"], h["dup"], agg["kingdom"])

        if is_big:
            for genus_key, agg in sample_genus_agg.items():
                dup_avg = agg["dup_wsum"] / agg["reads"] if agg["reads"] else 0.0
                genus_hits_big[genus_key].append(dict(
                    taxReads=agg["reads"], kmers=agg["kmers"],
                    cov=agg["cov"], dup=dup_avg, adj_kmers=agg["adj_kmers"],
                ))

    # ── sample_metrics.csv ──────────────────────────────────────────────
    sample_fields = list(sample_rows[0].keys()) if sample_rows else []
    with open(os.path.join(OUT_DIR, "sample_metrics.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=sample_fields)
        w.writeheader()
        w.writerows(sample_rows)

    # ── taxon prevalence tables ─────────────────────────────────────────
    def write_prevalence(taxon_hits, n_total, out_name):
        rows = []
        for (taxID, name, rank), records in taxon_hits.items():
            n_present = len(records)
            reads_vals = sorted(r["taxReads"] for r in records)
            rows.append(dict(
                taxID=taxID,
                name=name,
                rank=rank,
                n_samples_present=n_present,
                pct_samples_present=round(100 * n_present / n_total, 4) if n_total else 0,
                median_taxReads=round(percentile(reads_vals, 50), 2),
                max_taxReads=reads_vals[-1],
            ))
        rows.sort(key=lambda r: -r["n_samples_present"])
        with open(os.path.join(OUT_DIR, out_name), "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "taxID", "name", "rank", "n_samples_present",
                "pct_samples_present", "median_taxReads", "max_taxReads",
            ])
            w.writeheader()
            w.writerows(rows)
        return rows

    write_prevalence(taxon_hits_all, n_all, "taxon_prevalence_all.csv")
    write_prevalence(taxon_hits_big, n_big, "taxon_prevalence_min10kb.csv")

    # ── Genus-Ebene-Prävalenz (ein Sample zählt je Genus nur einmal) ────
    def write_genus_prevalence(genus_samples, n_total, out_name):
        rows = []
        for genus_taxid, samples in genus_samples.items():
            n_present = len(samples)
            rows.append(dict(
                genus_taxid=genus_taxid,
                genus_name=genus_names[genus_taxid],
                n_samples_present=n_present,
                pct_samples_present=round(100 * n_present / n_total, 4) if n_total else 0,
            ))
        rows.sort(key=lambda r: -r["n_samples_present"])
        with open(os.path.join(OUT_DIR, out_name), "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "genus_taxid", "genus_name", "n_samples_present", "pct_samples_present",
            ])
            w.writeheader()
            w.writerows(rows)
        return rows

    write_genus_prevalence(genus_samples_all, n_all, "taxon_prevalence_genus_all.csv")
    write_genus_prevalence(genus_samples_big, n_big, "taxon_prevalence_genus_min10kb.csv")

    # ── taxon_reference_stats.csv (percentile distributions, >=10kb cohort) ──
    # log10(reads+1) mean/SD per taxon -> z-score basis in generate_report_v3.py
    # ("how many SDs above typical background is this finding's read count").
    # Reads are strongly right-skewed, hence the log transform. A global pooled
    # SD (across ALL hit records, not per-taxon) is written to reference_meta.json
    # as a fallback for taxa with too few reference occurrences (n<3) to fit a
    # stable taxon-specific SD.
    import statistics

    def log_mean_sd(vals):
        logs = [math.log10(v + 1) for v in vals]
        n = len(logs)
        mean = sum(logs) / n
        sd = statistics.stdev(logs) if n >= 2 else 0.0
        return mean, sd

    def log_mean_sd_direct(vals, floor=1e-9):
        # For quantities that are always small fractions (coverage): log10(v+1)
        # would collapse to ~0 for every value and destroy the signal, so the
        # untransformed log10(v) is used instead (floored to avoid log(0)).
        logs = [math.log10(max(v, floor)) for v in vals]
        n = len(logs)
        mean = sum(logs) / n
        sd = statistics.stdev(logs) if n >= 2 else 0.0
        return mean, sd

    def pooled_log_sd(all_records, key="taxReads", direct=False):
        if direct:
            logs = [math.log10(max(r[key], 1e-9)) for r in all_records]
        else:
            logs = [math.log10(r[key] + 1) for r in all_records]
        return statistics.stdev(logs) if len(logs) >= 2 else 1.0

    def ref_stats_row(records, n_total):
        n_present = len(records)
        n_absent = max(n_total - n_present, 0)
        reads_vals = sorted(r["taxReads"] for r in records)
        kmers_vals = sorted(r["kmers"] for r in records)
        dup_vals = sorted(r["dup"] for r in records)
        cov_vals = sorted(r["cov"] for r in records)
        # Both the z-score log-mean/SD AND the boxplot percentiles must
        # reflect the FULL cohort, not just the samples where this taxon/
        # genus happened to be observed: a genus seen in only 4/437 samples
        # is ABSENT (no signal) in the other 433, and "typical background"
        # has to include that near-universal absence. Computing these only
        # from the n_present positive records silently conditions on
        # "detected at all" and badly understates how rare a genuine finding
        # is (reported case: Cytomegalovirus, 4/437 samples -- z-score far
        # too low for how seldom it appears, and its boxplot only showed the
        # 4 positive detections with no visual indication of the other 433
        # absences). Fixed by zero-padding reads/kmers/cov with
        # (n_total - n_present) explicit absences before both the log
        # transform and the percentile calculation (log10(0+1)=0 for counts;
        # log_mean_sd_direct's floor handles cov's 0 the same way). An
        # intermediate shrinkage-estimator version (partial correction,
        # tuned for best TP/FP separation) was tried and reverted on
        # request in favour of this full, statistically unweighted
        # population accounting -- see reference_cohort_filter_strategy.md
        # for that comparison (full zero-padding measured AUC 0.579 vs.
        # shrinkage's 0.817 on the reference cohort; the drop is expected
        # and accepted here, not a bug: most genera in this cohort are
        # individually rare (median 10/437 samples), so an unweighted
        # population mean pushes the majority of background/contaminant
        # hits to a high z-score too, same as it correctly does for true
        # rare pathogens -- z-score alone is no longer a reliable Top-Hits
        # filter under this accounting, which is why Kontaminant-tagged
        # genus rows are now excluded from Top-Hits outright regardless of
        # z-score, see select_top_hits() in generate_report_v3.py).
        #
        # dup has no natural "absent" value (it is a ratio defined only
        # when reads exist at all), so its percentiles/boxplot intentionally
        # stay presence-conditional -- unlike reads/kmers/cov, "0 dup" would
        # not mean anything.
        reads_vals_full = sorted(reads_vals + [0] * n_absent)
        kmers_vals_full = sorted(kmers_vals + [0] * n_absent)
        cov_vals_full = sorted(cov_vals + [0.0] * n_absent)

        reads_log_mean, reads_log_sd = log_mean_sd(reads_vals_full)
        kmers_log_mean, kmers_log_sd = log_mean_sd(kmers_vals_full)
        cov_log_mean, cov_log_sd = log_mean_sd_direct(cov_vals_full)
        row = dict(
            n_samples_present=n_present,
            pct_samples_present=round(100 * n_present / n_total, 4) if n_total else 0,
            reads_log_mean=round(reads_log_mean, 4),
            reads_log_sd=round(reads_log_sd, 4),
            kmers_log_mean=round(kmers_log_mean, 4),
            kmers_log_sd=round(kmers_log_sd, 4),
            cov_log_mean=round(cov_log_mean, 4),
            cov_log_sd=round(cov_log_sd, 4),
            reads_p50=round(percentile(reads_vals_full, 50), 2),
            reads_p90=round(percentile(reads_vals_full, 90), 2),
            reads_p95=round(percentile(reads_vals_full, 95), 2),
            reads_p99=round(percentile(reads_vals_full, 99), 2),
            reads_max=reads_vals_full[-1],
            # Full 5-number summary for kmers/cov -> mini-boxplots in
            # generate_report_v3.py, zero-padded like the log-mean above.
            kmers_min=kmers_vals_full[0],
            kmers_p25=round(percentile(kmers_vals_full, 25), 2),
            kmers_p50=round(percentile(kmers_vals_full, 50), 2),
            kmers_p75=round(percentile(kmers_vals_full, 75), 2),
            kmers_p95=round(percentile(kmers_vals_full, 95), 2),
            kmers_max=kmers_vals_full[-1],
            # Full 5-number summary for dup -> mini-boxplot in generate_report_v3.py.
            # Presence-conditional (NOT zero-padded, see comment above).
            dup_min=dup_vals[0],
            dup_p25=round(percentile(dup_vals, 25), 3),
            dup_p50=round(percentile(dup_vals, 50), 3),
            dup_p75=round(percentile(dup_vals, 75), 3),
            dup_p95=round(percentile(dup_vals, 95), 3),
            dup_max=dup_vals[-1],
            cov_min=cov_vals_full[0],
            cov_p25=round(percentile(cov_vals_full, 25), 10),
            cov_p50=round(percentile(cov_vals_full, 50), 10),
            cov_p75=round(percentile(cov_vals_full, 75), 10),
            cov_p95=round(percentile(cov_vals_full, 95), 10),
            cov_max=cov_vals_full[-1],
        )
        # adj_kmers (coverage-adjusted k-mers, kraken_tree.adjusted_kmers()) is
        # only computed for genus-level records (see sample_genus_agg above) --
        # z-score basis in generate_report_v3.py, replacing the earlier
        # reads-based z-score (k-mers separate TP/FP better than reads, see
        # analysis/threshold_analysis_report.md; dup is folded in to weight
        # broad, low-duplication coverage over narrow, repetitive matches).
        if records and "adj_kmers" in records[0]:
            adj_vals_full = sorted([r["adj_kmers"] for r in records] + [0] * n_absent)
            adj_log_mean, adj_log_sd = log_mean_sd(adj_vals_full)
            row["adj_kmers_log_mean"] = round(adj_log_mean, 4)
            row["adj_kmers_log_sd"] = round(adj_log_sd, 4)
        return row

    pooled_reads_log_sd = pooled_log_sd(
        [r for records in taxon_hits_big.values() for r in records])
    pooled_reads_log_sd_genus = pooled_log_sd(
        [r for records in genus_hits_big.values() for r in records])
    pooled_adj_kmers_log_sd_genus = pooled_log_sd(
        [r for records in genus_hits_big.values() for r in records], key="adj_kmers")
    pooled_kmers_log_sd_genus = pooled_log_sd(
        [r for records in genus_hits_big.values() for r in records], key="kmers")
    pooled_cov_log_sd_genus = pooled_log_sd(
        [r for records in genus_hits_big.values() for r in records], key="cov", direct=True)

    with open(os.path.join(OUT_DIR, "reference_meta.json"), "w", encoding="utf-8") as fh:
        json.dump(dict(
            generated=date.today().isoformat(),
            n_samples_all=n_all,
            n_samples_min10kb=n_big,
            size_threshold_bytes=SIZE_THRESHOLD_BYTES,
            human_excluded=True,
            zero_padded_absent_samples=True,
            pooled_reads_log_sd=round(pooled_reads_log_sd, 4),
            pooled_reads_log_sd_genus=round(pooled_reads_log_sd_genus, 4),
            pooled_adj_kmers_log_sd_genus=round(pooled_adj_kmers_log_sd_genus, 4),
            pooled_kmers_log_sd_genus=round(pooled_kmers_log_sd_genus, 4),
            pooled_cov_log_sd_genus=round(pooled_cov_log_sd_genus, 4),
        ), fh, indent=2)

    ref_rows = []
    for (taxID, name, rank), records in taxon_hits_big.items():
        ref_rows.append(dict(taxID=taxID, name=name, rank=rank,
                              **ref_stats_row(records, n_big)))
    ref_rows.sort(key=lambda r: -r["n_samples_present"])
    with open(os.path.join(OUT_DIR, "taxon_reference_stats.csv"), "w", newline="", encoding="utf-8") as fh:
        fieldnames = list(ref_rows[0].keys()) if ref_rows else []
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(ref_rows)

    # ── genus_reference_stats.csv (percentile distributions, >=10kb cohort,
    # aggregiert über alle Species/Strains je Gattung und Sample) ──────────
    genus_ref_rows = []
    for genus_taxid, records in genus_hits_big.items():
        genus_ref_rows.append(dict(
            genus_taxid=genus_taxid, genus_name=genus_names[genus_taxid],
            **ref_stats_row(records, n_big)))
    genus_ref_rows.sort(key=lambda r: -r["n_samples_present"])
    with open(os.path.join(OUT_DIR, "genus_reference_stats.csv"), "w", newline="", encoding="utf-8") as fh:
        fieldnames = list(genus_ref_rows[0].keys()) if genus_ref_rows else []
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(genus_ref_rows)

    # ── sample-level reference percentiles (for QC placement of new runs) ──
    def col_percentiles(rows, key):
        vals = sorted(r[key] for r in rows)
        return {p: percentile(vals, p) for p in (5, 10, 25, 50, 75, 90, 95)}

    big_rows = [r for r in sample_rows if r["pass_min10kb"]]
    qc_cols = [
        "total_reads", "pct_classified", "pct_human",
        "nonhuman_classified_reads", "n_distinct_taxa_hits",
        "total_nonhuman_hit_reads",
    ]
    with open(os.path.join(OUT_DIR, "sample_qc_reference_percentiles.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "n", "p5", "p10", "p25", "p50", "p75", "p90", "p95"])
        for col in qc_cols:
            pc = col_percentiles(big_rows, col)
            w.writerow([col, len(big_rows)] + [round(pc[p], 4) for p in (5, 10, 25, 50, 75, 90, 95)])

    print(f"Samples gesamt: {n_all}")
    print(f"Samples >= 10 KB: {n_big}")
    print(f"Distinkte Taxa (alle Samples): {len(taxon_hits_all)}")
    print(f"Distinkte Taxa (>=10kb Kohorte): {len(taxon_hits_big)}")
    print(f"Distinkte Genera (alle Samples): {len(genus_samples_all)}")
    print(f"Distinkte Genera (>=10kb Kohorte, mit Referenzstatistik): {len(genus_hits_big)}")
    print(f"Output-Verzeichnis: {os.path.abspath(OUT_DIR)}")


if __name__ == "__main__":
    main()
