"""
Gemeinsame Parsing-Logik für KrakenUniq-Reports.

KrakenUniq-Reports sind ein eingerückter Taxonomiebaum (root -> ... ->
genus -> species/strain/...). Die Einrückung der taxName-Spalte codiert die
Tiefe im Baum. parse_kraken_rows() nutzt das, um jeder Zeile ihren nächsten
Genus-Vorfahren zuzuordnen (für eine Genus-Ebene-Prävalenzfilterung, robust
gegenüber der Vielzahl an Stamm-/Isolat-Namen innerhalb einer Gattung).

Wird sowohl von build_reference_db.py als auch generate_report_v3.py
verwendet, damit die Genus-Zuordnung an beiden Stellen identisch ist.
"""

import math

COLUMN_HEADER_PREFIX = "%\treads\ttaxReads"


def adjusted_kmers(kmers, dup, kingdom):
    """
    Coverage-adjusted k-mer count, used as the basis for the genus-level
    z-score (build_reference_db.py + generate_report_v3.py). k-mers alone
    (distinct k-mers matched) separates true pathogens from background
    better than raw reads (see analysis/threshold_analysis_report.md), but
    does not distinguish "many distinct genome positions covered once each"
    (broad, credible coverage) from "few positions hit over and over"
    (narrow/repetitive, less credible) -- dup (mean multiplicity per
    distinct k-mer) captures exactly that, so it is used here to discount
    the k-mer count for high duplication.

    Viruses are exempted (kmers returned unchanged): their genomes are
    small enough that very high dup is routinely produced by genuinely
    deep, broad sequencing coverage rather than repetitive/contaminant
    cross-mapping -- the same reasoning already used to drop the dup
    evidence-filter cap for Viren (see generate_report_v3.assess()). A
    sqrt (rather than linear) discount is used for the other kingdoms so a
    single very high dup value cannot collapse the score to near zero.
    """
    if kingdom == "Viren":
        return kmers
    return kmers / math.sqrt(max(dup, 1.0))


def infer_kingdom(name):
    n = name.lower()
    virus_kw = [
        "virus", "phage", "viridae", "viricota", "viricetes", "virales",
        "herpesvirus", "betaherpes", "alphaherpes", "gammaherpes",
        "polyomavirus", "circovirus", "lyssavirus", "morbillivirus",
        "coronavirus", "gemycircularvirus",
    ]
    fungus_kw = [
        "ascomycota", "basidiomycota", "saccharomycetes", "tremellomycetes",
        "trichosporon", "candida", "aspergillus", "cryptococcus",
        "histoplasma", "microsporidia", "tubulinosema", "anncaliia",
        "spraguea", "mucoromycota", "mucorales", "rhizopus",
    ]
    parasite_kw = [
        "plasmodium", "toxoplasma", "cryptosporidium", "trypanosoma",
        "leishmania", "naegleria", "acanthamoeba", "apicomplexa",
        "kinetoplastea", "coccidia", "eimeria", "giardia", "entamoeba",
    ]
    if any(k in n for k in virus_kw):    return "Viren"
    if any(k in n for k in parasite_kw): return "Parasiten"
    if any(k in n for k in fungus_kw):   return "Pilze"
    return "Bakterien"


def parse_kraken_rows(text):
    """
    Parst den Tabellenteil eines KrakenUniq-Reports.
    Gibt eine Liste von dicts zurück: pct, reads, taxReads, kmers, dup, cov,
    taxID, rank, name, indent, genus_taxid, genus_name.

    genus_taxid/genus_name: taxID/Name des nächsten Vorfahren mit rank=="genus"
    im Taxonomiebaum. Falls kein Genus-Vorfahre existiert (z.B. Viren ohne
    Genus-Rang, oder Knoten oberhalb der Genus-Ebene), fällt es auf die
    Zeile selbst zurück (genus_taxid = eigene taxID).

    Manche Reports enthalten die komplette Tabelle versehentlich zweimal
    (Konkatenations-Artefakt) – ab dem zweiten Auftreten der Spaltenüber-
    schrift wird das Parsen abgebrochen.
    """
    rows = []
    stack = []  # dicts: indent, taxID, name, rank
    seen_col_header = False

    for line in text.split("\n"):
        if line.startswith("#") or not line.strip():
            continue
        if line.startswith(COLUMN_HEADER_PREFIX):
            if seen_col_header:
                break
            seen_col_header = True
            continue

        cols = line.split("\t")
        if len(cols) != 9:
            # Rare KrakenUniq report corruption: stray extra tab(s) or two
            # numbers merged into one field around very large dup/cov values
            # (e.g. "8.3e+03" rows). taxID/rank/name are reliably the last
            # three tab-separated fields regardless; the six numeric fields
            # in front are recovered if dropping empty fields yields exactly
            # six, otherwise zeroed. Zeroed taxReads can never register as a
            # hit (taxReads>0 required), but the node stays in the ancestor
            # stack so genus resolution for its children is not broken.
            if len(cols) < 3:
                continue
            id_col, rank_col, name_col = cols[-3], cols[-2], cols[-1]
            numeric = [c for c in cols[:-3] if c != ""]
            if len(numeric) == 6:
                cols = numeric + [id_col, rank_col, name_col]
            else:
                cols = ["0", "0", "0", "0", "0", "0", id_col, rank_col, name_col]
        try:
            pct      = float(cols[0])
            reads    = int(cols[1])
            taxReads = int(cols[2])
            kmers    = int(cols[3])
            dup      = float(cols[4])
            cov      = float(cols[5])
            taxID    = cols[6].strip()
            rank     = cols[7].strip()
            raw_name = cols[8]
        except (ValueError, IndexError):
            continue

        indent = len(raw_name) - len(raw_name.lstrip(" "))
        name = raw_name.strip()

        while stack and stack[-1]["indent"] >= indent:
            stack.pop()

        genus_taxid, genus_name = None, None
        for anc in reversed(stack):
            if anc["rank"] == "genus":
                genus_taxid, genus_name = anc["taxID"], anc["name"]
                break
        if genus_taxid is None:
            genus_taxid, genus_name = taxID, name

        rows.append(dict(
            pct=pct, reads=reads, taxReads=taxReads, kmers=kmers,
            dup=dup, cov=cov, taxID=taxID, rank=rank, name=name,
            indent=indent, genus_taxid=genus_taxid, genus_name=genus_name,
        ))
        stack.append(dict(indent=indent, taxID=taxID, name=name, rank=rank))

    return rows
