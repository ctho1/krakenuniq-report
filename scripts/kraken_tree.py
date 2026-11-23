"""
Shared parsing logic for KrakenUniq reports.

KrakenUniq reports are an indented taxonomy tree (root -> ... ->
genus -> species/strain/...). The indentation of the taxName column encodes
the depth in the tree. parse_kraken_rows() uses this to assign each row its
nearest genus ancestor (for genus-level prevalence filtering, robust
against the large number of strain/isolate names within a genus).

Used by both build_reference_db.py and generate_report_v3.py, so that the
genus assignment is identical in both places.
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
    evidence-filter cap for viruses (see generate_report_v3.assess()).
    kingdom is the category from category_from_lineage(). A
    sqrt (rather than linear) discount is used for the other kingdoms so a
    single very high dup value cannot collapse the score to near zero.
    """
    if kingdom == VIRUSES:
        return kmers
    return kmers / math.sqrt(max(dup, 1.0))


# ── Category (kingdom) from the taxonomy tree ────────────────────────────────
# Each hit is assigned to a category via its ancestors in the KrakenUniq
# taxonomy tree (superkingdom, and kingdom "Fungi" within Eukaryota) -- not
# via its name.
VIRUSES = "Viruses"
BACTERIA = "Bacteria"
ARCHAEA = "Archaea"
FUNGI = "Fungi"
OTHER_EUKARYOTES = "Parasites and other eukaryotes"
UNASSIGNED = "Unassigned"
CATEGORY_ORDER = [VIRUSES, BACTERIA, ARCHAEA, FUNGI, OTHER_EUKARYOTES, UNASSIGNED]

# Top-level ranks: "superkingdom" in the NCBI taxonomy used by the current
# database; "domain"/"acellular root" in newer NCBI taxonomy releases.
_TOP_RANKS = {"superkingdom", "domain", "acellular root"}
_TOP_CATEGORY = {"Viruses": VIRUSES, "Bacteria": BACTERIA, "Archaea": ARCHAEA}


def category_from_lineage(lineage):
    """
    lineage: list of (rank, name) from the root down to the taxon itself.
    Returns one of CATEGORY_ORDER:
      superkingdom Viruses / Bacteria / Archaea -> that category
      superkingdom Eukaryota -> Fungi if below kingdom "Fungi",
                                otherwise OTHER_EUKARYOTES
      no superkingdom ancestor (e.g. "other entries") -> UNASSIGNED
    """
    top = next((name for rank, name in lineage if rank in _TOP_RANKS), None)
    if top in _TOP_CATEGORY:
        return _TOP_CATEGORY[top]
    if top == "Eukaryota":
        if any(rank == "kingdom" and name == "Fungi" for rank, name in lineage):
            return FUNGI
        return OTHER_EUKARYOTES
    return UNASSIGNED


def parse_kraken_rows(text):
    """
    Parses the table part of a KrakenUniq report.
    Returns a list of dicts: pct, reads, taxReads, kmers, dup, cov,
    taxID, rank, name, indent, genus_taxid, genus_name, kingdom.

    genus_taxid/genus_name: taxID/name of the nearest ancestor with
    rank=="genus" in the taxonomy tree. If no genus ancestor exists (e.g.
    viruses without a genus rank, or nodes above genus level), it falls back
    to the row itself (genus_taxid = own taxID).

    kingdom: category from the ancestors in the taxonomy tree, see
    category_from_lineage().

    Some reports accidentally contain the complete table twice
    (concatenation artefact) -- parsing stops at the second occurrence of
    the column header.
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
        taxID    = cols[6].strip()
        rank     = cols[7].strip()
        raw_name = cols[8]
        try:
            pct      = float(cols[0])
            reads    = int(cols[1])
            taxReads = int(cols[2])
            kmers    = int(cols[3])
            dup      = float(cols[4])
            cov      = float(cols[5])
        except ValueError:
            # Same corruption with exactly 9 fields (e.g. an empty field plus
            # two numbers merged into one, seen on a "Viruses" superkingdom
            # row): numbers are zeroed, but the node is kept in the ancestor
            # stack -- dropping it would attach all of its descendants to
            # the wrong branch (wrong genus/category).
            pct, reads, taxReads, kmers, dup, cov = 0.0, 0, 0, 0, 0.0, 0.0

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

        lineage = [(anc["rank"], anc["name"]) for anc in stack] + [(rank, name)]

        rows.append(dict(
            pct=pct, reads=reads, taxReads=taxReads, kmers=kmers,
            dup=dup, cov=cov, taxID=taxID, rank=rank, name=name,
            indent=indent, genus_taxid=genus_taxid, genus_name=genus_name,
            kingdom=category_from_lineage(lineage),
        ))
        stack.append(dict(indent=indent, taxID=taxID, name=name, rank=rank))

    return rows
