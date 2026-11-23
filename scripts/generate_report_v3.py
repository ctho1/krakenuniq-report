"""
KrakenUniq PDF Report Generator – v5
Changes vs. v4:
  - "Top hits" table: a fixed number of prioritised hits (rarest first,
    topped up with more common ones if needed) reliably fills page 1
    instead of leaving a lot of white space when there are few rare hits
  - Mini boxplots (reference distribution, >=10kb cohort) for k-mers and
    coverage per hit, with a red dot for the current value
  - Remaining hits (rare ones beyond the top N + common ones) in a compact
    supplement, split into "additional rare" and "common" findings

Changes vs. v3:
  - Slimmer header (title/meta/stat tiles more compact)
  - Automatically generated narrative findings text removed
  - Reference prevalence (n / % of the reference database, genus level) per
    hit instead of a plausibility badge ("Detected" / "Notable" removed)
  - Comment column: Pathogen / Commensal / Contaminant
"""

import csv
import json
import math
import os
import re
import sys
import tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kraken_tree
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor, white
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    HRFlowable, KeepTogether, PageBreak, CondPageBreak,
)
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT
from reportlab.pdfbase import pdfmetrics
from reportlab.graphics.shapes import Drawing, Line, Rect, Circle

# ── Palette ──────────────────────────────────────────────────────────────────
C_WHITE       = HexColor("#FFFFFF")
C_NAVY        = HexColor("#1D1D1F")
C_MID         = HexColor("#6E6E73")
C_LIGHT       = HexColor("#AEAEB2")
C_RULE        = HexColor("#E5E5EA")
C_RULE_DARK   = HexColor("#C7C7CC")

C_TILE_BLUE   = HexColor("#0071E3")
C_TILE_GREEN  = HexColor("#34C759")
C_TILE_GREY   = HexColor("#636366")

C_GREEN_FILL  = HexColor("#F0FDF4")
C_GREEN_TEXT  = HexColor("#1A7F3C")
C_ORANGE_FILL = HexColor("#FFF7ED")
C_ORANGE_TEXT = HexColor("#B45309")
C_GREY_FILL   = HexColor("#F5F5F7")
C_GREY_TEXT   = HexColor("#636366")
C_RED_FILL    = HexColor("#FEF2F2")
C_RED_TEXT    = HexColor("#B91C1C")
C_PURPLE_TEXT = HexColor("#7C3AED")

C_TABLE_HEAD  = HexColor("#1D1D1F")
C_ROW_EVEN    = HexColor("#FFFFFF")
C_ROW_ODD     = HexColor("#F9F9FB")

PAGE_W, PAGE_H = A4
ML = 20*mm
MR = 20*mm
MT = 16*mm      # slimmer top margin (was 22mm)
MB = 20*mm
CW = PAGE_W - ML - MR  # ≈ 170 mm

# Column widths (mm) – top-hits table (page 1, incl. mini boxplots)
# The top-3 species are shown as a grey sub-line in the organism cell, not
# as a separate column. k-mers/dup/Cov. show the number and the boxplot
# stacked in ONE column (saves width compared to separate value/chart
# columns).
# Total must equal CW (170 mm)
COL_W = {
    "organism":  35*mm,
    "reads":     14*mm,
    "kmers":     24*mm,
    "dup":       20*mm,
    "cov":       26*mm,
    "zscore":    22*mm,
    "comment":   29*mm,
}
assert abs(sum(COL_W.values()) - CW) < 1, \
    f"Col widths sum {sum(COL_W.values())/mm:.1f}mm ≠ CW {CW/mm:.1f}mm"
# Horizontal boxplots: width = full cell width (KMERS_INNER_W/COV_INNER_W/
# DUP_INNER_W), small height to keep the row height of the top-hits table
# low.
BOX_H = 4*mm
KMERS_INNER_W = COL_W["kmers"] - 4*mm
DUP_INNER_W   = COL_W["dup"]   - 4*mm
COV_INNER_W   = COL_W["cov"] - 4*mm


# Prevalence thresholds for the tier structure (relative to the complete
# reference database, not only the >=10kB cohort)
RARE_THRESHOLD_1PCT = 1.0
RARE_THRESHOLD_5PCT = 5.0

# z-score thresholds (standard deviations from the taxon's log10(reads+1)
# background noise in the >=10kB cohort) for the colour coding in the
# top-hits table.
Z_THRESHOLD_HIGH = 3.0
Z_THRESHOLD_ELEVATED = 1.5

# Number of prioritised hits in the top-hits table on page 1. Topped up
# with the next most common hits when there are only a few rare hits, so
# that page 1 is reliably filled instead of showing a lot of white space.
TOP_N_HITS = 15

# Minimum reads above which a curated Pathogen/Commensal hit is always
# included in the top hits regardless of its z-score (see
# select_top_hits()). Set lower for "Pathogen" (established CNS pathogens,
# small curated list) than for "Commensal" (human resident flora, higher
# probability of contamination -> more evidence required). Since switching
# to pure z-score sorting (without the genus-prevalence tier as a pre-sort)
# this specifically rescues cases with a very low but genuine signal: the
# reference-cohort validation otherwise placed HIV-1 (9 reads), BoDV-1
# (5/8 reads) and SARS-CoV-2 (5 reads) outside the top 15, because other
# hits in their (partly highly polymicrobial) samples reached a higher
# z-score.
FORCED_MIN_READS = {"Pathogen": 5, "Commensal": 20}
# Cap on the number of forced hits, so that polymicrobially colonised
# samples do not fill the entire top-hits table with commensals.
FORCED_MAX_N = 10
# Minimum number of slots in the top-hits table that are always reserved
# for pure genus-rarity prioritisation (cannot be displaced by forced
# hits).
MIN_RARITY_SLOTS = 5

# Additional inclusion criterion ONLY for the top-hits table (page 1), on
# request: a genus row (read-weighted mean dup value, see group_by_genus())
# with dup >= TOP_HITS_MAX_DUP is no longer included in the top hits
# (neither forced nor via z-score), but still appears in the supplement as
# long as it passes the evidence filter -- the evidence filter itself
# (assess()) is unchanged, so a hit does not disappear from the report
# entirely. Kingdom "Viruses" is exempt, for the same reason as the dup cap
# in the evidence filter (VIRAL_MAX_DUP): small viral genomes structurally
# produce high dup values under genuine deep coverage (e.g. BoDV-1 dup
# 50-400, TBEV dup 2920, JC virus dup 215 in the reference cohort) -- an
# unfiltered dup<2 would have displaced practically all confirmed viral top
# hits.
TOP_HITS_MAX_DUP = 2.0

# ── Minimum criteria for inclusion in the report (noise suppression) ───────
THRESHOLDS = {
    "CONVINCING": dict(minReads=10, minKmers=500, maxDup=5,  minCov=1e-4),
    "NOTABLE":    dict(minReads=5,  minKmers=50,  maxDup=30, minCov=1e-6),
}

# dup ("mean k-mer duplication rate") scales with sequencing depth per
# genome size. Small viral genomes structurally produce very high dup values
# (>1000) under genuinely very high coverage (e.g. TBEV with 750,000+ reads
# on a ~11 kb genome) without this indicating contamination/cross-
# reactivity. Reference-cohort validation (63 cases with a confirmed
# pathogen, Supplementary Table 2): 6 of 6 viral cases with dup>30 were
# genuine, highly covered pathogens (TBEV dup=2920, BoDV-1 dup=50-400, JCV
# dup=215) -- the dup<=30 cutoff would have discarded them despite tens of
# thousands of reads. The dup cutoff is therefore lifted for kingdom
# "Viruses"; for bacteria/fungi/parasites (much larger genomes) dup remains
# a meaningful contamination/cross-reactivity indicator and is still
# applied.
VIRAL_MAX_DUP = math.inf

RELEVANT_RANKS = {
    "species", "subspecies", "strain", "isolate",
    "no rank", "varietas", "serotype",
    "species group", "species subgroup",
}

SKIP_TAXA = {
    "Homo sapiens", "synthetic construct", "unclassified",
    "root", "cellular organisms", "other entries",
    "other sequences", "artificial sequences",
}

# Contaminant reference:
# Laurence M et al. (2014) PLoS ONE 9(5):e97876. doi:10.1371/journal.pone.0097876
# Supplemented with well-established skin flora and environmental contaminants
# documented in clinical metagenomics of tissue samples.
#
# THREE-TIER CLASSIFICATION:
# LIKELY_CONTAMINANTS  – no human pathogen potential; pure env./reagent/skin contaminants
# OPPORTUNISTIC_FLORA  – human commensals with opportunistic CNS pathogen potential
#                        (e.g. oral streptococci in brain abscess); labelled separately
# KNOWN_PATHOGENS      – established CNS pathogens (defined below)

LIKELY_CONTAMINANTS = {
    # ── Ultrapure water / reagent contaminants (Laurence et al. 2014) ─────────
    "Bradyrhizobium japonicum",
    "Bradyrhizobium elkanii",
    "Bradyrhizobium sp. DFCI-1",
    "Methylobacterium extorquens",
    "Methylobacterium radiotolerans",
    "Methylobacterium mesophilicum",
    "Methylobacterium tardum",
    "Methylorubrum populi",
    "Agrobacterium fabrum",
    "Agrobacterium tumefaciens",
    "Rhizobium etli",
    "Sphingomonas aerolata",
    "Sphingomonas aliaeris",
    "Sphingomonas morindae",
    "Sphingopyxis sp. DBS4",
    "Ralstonia pickettii",
    "Ralstonia mannitolilytica",
    "Cupriavidus metallidurans",
    "Xanthomonas campestris",
    "Stenotrophomonas maltophilia",
    "Schlegelella aquatica",
    "Limnobacter sp. SAORIC-580",
    "Moraxella osloensis",             # environmental Moraxella, reagent contaminant
    # ── Skin flora / tissue processing contaminants ───────────────────────────
    "Cutibacterium acnes",
    "Cutibacterium granulosum",
    "Cutibacterium avidum",
    "Cutibacterium modestum",
    "Cutibacterium acnes subsp. acnes",
    "Cutibacterium acnes subsp. elongatum",
    "unclassified Plantactinospora",
    "Staphylococcus epidermidis",
    "Staphylococcus capitis",
    "Staphylococcus caprae",
    "Staphylococcus warneri",
    "Staphylococcus hominis",
    "Staphylococcus haemolyticus",
    "Staphylococcus pasteuri",
    "Mammaliicoccus sciuri",
    "Corynebacterium tuberculostearicum",
    "Corynebacterium kefirresidentii",
    "Corynebacterium sp. SCR221107",
    "Corynebacterium ureicelerivorans",
    "Corynebacterium appendicis",
    "Brevibacterium casei",
    "Micrococcus luteus",
    "Micrococcus yunnanensis",
    "Kocuria rhizophila",
    "Winkia neuii",
    "Malassezia restricta",
    "Malassezia globosa",
    "Malassezia restricta CBS 7877",
    "Malassezia globosa CBS 7966",
    # ── Food / fermentation organisms (no human CNS pathogen potential) ──────
    "Streptococcus thermophilus",
    "Lactococcus cremoris",
    "Lactococcus lactis",
    "Ligilactobacillus salivarius",
    "Lactiplantibacillus plantarum",
    "Lacticaseibacillus rhamnosus",
    "Lactobacillus crispatus",
    "Leuconostoc carnosum",
    # ── Environmental / soil organisms ───────────────────────────────────────
    "Wolbachia pipientis",
    "Bdellovibrio bacteriovorus",
    "Candidatus Carsonella ruddii",
    "Hanseniaspora guilliermondii",
    "Plantactinospora sp. BB1",
    "Plantactinospora sp. BC1",
    "Cloacibacterium caeni",
    "Cloacibacterium normanense",
    "Streptomyces sp. A10(2020)",
    "Streptomyces albidoflavus",
    "Lysobacter capsici",
    "Rothia mucilaginosa",
    # Rothia dentocariosa: removed from LIKELY_CONTAMINANTS – classified as
    # OPPORTUNISTIC_FLORA due to rare brain abscess reports. assess() prioritises
    # is_opportunist over is_contaminant via OPPORTUNISTIC_FLORA set.
    "Rothia aeria",
    "Novosphingobium humi",
    "Paracoccus sp. MC1862",
    "Paracoccus yeei",
    "Paracoccus sanguinis",
    "Haematobacter massiliensis",
    "Paraburkholderia aromaticivorans",
    "Paraburkholderia fungorum",
    # ── Non-human parasites / animal pathogens ────────────────────────────────
    "Anncaliia algerae", "Spraguea lophii",
    "Eimeria acervulina", "Eimeria tenella",
    "Trypanosoma evansi", "Cryptosporidium baileyi",
    # ── Non-human viruses / phages ────────────────────────────────────────────
    "Bat gemycircularvirus", "Pleurochrysis carterae circular virus",
    "Emiliania huxleyi virus 86", "Gemykibivirus hipla1",
    "Gallid alphaherpesvirus 2", "Cyprinid herpesvirus 2",
    "Vibrio phage ValB1MD-2", "Colobine gammaherpesvirus 1",
    # ── Environmental / soil / water organisms (no human pathogen potential) ──
    # Identified from reference-cohort false-positive frequency analysis
    # (analysis/threshold_dataset.csv, high-occurrence uncommented genera);
    # each species checked individually for CNS/human pathogen case-report
    # history before inclusion -- mixed-pathogenicity genera (e.g. Bacillus,
    # Vibrio, Burkholderia, Actinomyces, Mycobacterium, Clostridium,
    # Bordetella, Brucella, Salmonella) were deliberately NOT blanket-tagged,
    # since group_by_genus() already lets a real Pathogen/Commensal member
    # species override a Contaminant one within the same genus row.
    "Edwardsiella ictaluri",            # fish pathogen only, not human-pathogenic
    "Curtobacterium citreum",
    "Curtobacterium flaccumfaciens",
    "unclassified Curtobacterium",
    "Acidovorax temperans",
    "Deinococcus geothermalis DSM 11300",
    "Delftia acidovorans",
    "Delftia tsuruhatensis",
    "Stutzerimonas stutzeri",
    "Stutzerimonas stutzeri group",
    "Brevundimonas mediterranea",
    "unclassified Brevundimonas",
    "Epilithonimonas vandammei",
    "Rhodopseudomonas palustris",
    "Tepidimonas taiwanensis",
    "Microbacterium aurum",
    "unclassified Microbacterium",
    "Sorangium cellulosum",
    "Variovorax paradoxus",
    "unclassified Variovorax",
    "Herbaspirillum huttiense",
    # ── Insect endosymbionts (obligate, cannot infect humans) ──────────────────
    "Buchnera aphidicola",
    "Candidatus Karelsulcia muelleri",
    "Candidatus Purcelliella pentastirinorum",
    # ── Gut / oral commensal flora (no CNS pathogen potential) ─────────────────
    "Oxalobacter aliiformigenes",
    "Dolosigranulum pigrum",
    "Bifidobacterium thermophilum",
    "Bifidobacterium adolescentis",
    "Bifidobacterium dentium",
    "Phocaeicola dorei",
    "Phocaeicola vulgatus",
    "Lautropia mirabilis",
    # ── Skin flora (sebaceous-gland-associated, cf. Cutibacterium) ─────────────
    "Lawsonella clevelandensis",
}

# Opportunistic flora: human commensals that can cause CNS infections
# (brain abscess, meningitis in immunocompromised) but are also frequent
# biopsy/FFPE contaminants. Labelled "Commensal" in the comment column to
# distinguish them from pure environmental contaminants and primary pathogens.
OPPORTUNISTIC_FLORA = {
    # ── Oral streptococci (brain abscess, endocarditis-related meningitis) ────
    "Streptococcus mitis",
    "Streptococcus oralis",
    "Streptococcus salivarius",
    "Streptococcus sanguinis",
    "Streptococcus sp. LPB0220",
    "Streptococcus intermedius",       # most common oral Strep in brain abscess
    "Streptococcus anginosus",
    "Streptococcus constellatus",
    "Gemella sanguinis",
    "Granulicatella adiacens",
    # ── Anaerobes (brain abscess, polymicrobial) ──────────────────────────────
    "Finegoldia magna",
    "Veillonella rogosae",
    "Veillonella parvula",
    "Lachnoanaerobaculum gingivalis",
    "Fusobacterium nucleatum",
    "Fusobacterium necrophorum",
    "Prevotella melaninogenica",
    "Prevotella intermedia",
    # ── Oral anaerobes / mixed flora ──────────────────────────────────────────
    "Schaalia odontolytica",
    "Rothia dentocariosa",             # rare brain abscess reports
    "Haemophilus parainfluenzae",      # rare CNS infections
    # ── Coagulase-negative staphylococci (device-related, post-neurosurgical) ─
    "Staphylococcus aureus",           # primary pathogen but also contaminant
    "Staphylococcus lugdunensis",
    # ── Gram-negative opportunists ────────────────────────────────────────────
    "Klebsiella pneumoniae",
    "Escherichia coli",
    "Acinetobacter baumannii",
    "Pseudomonas aeruginosa",
    "Enterococcus faecalis",
    "Enterococcus faecium",
    "Enterococcus cecorum",
}

KNOWN_PATHOGENS = {
    # CNS-relevant primary pathogens (curated list for signal enhancement)
    # Reference: CZ ID pathogen list (https://czid.org/pathogen_list, Chan Zuckerberg Initiative)
    # ── Viruses ───────────────────────────────────────────────────────────────
    "Cytomegalovirus humanbeta5",   # CMV / HHV-5
    "Human betaherpesvirus 5",
    "Roseolovirus humanbeta6a",     # HHV-6A
    "Human betaherpesvirus 6A",
    "Roseolovirus humanbeta6b",     # HHV-6B
    "Human betaherpesvirus 6B",
    "Human betaherpesvirus 6",
    "Roseolovirus humanbeta7",      # HHV-7
    "Simplexvirus humanalpha1",     # HSV-1
    "Human alphaherpesvirus 1",
    "Simplexvirus humanalpha2",     # HSV-2
    "Human alphaherpesvirus 2",
    "Varicellovirus humanalpha3",   # VZV
    "Human alphaherpesvirus 3",
    "Lymphocryptovirus humangamma4", # EBV
    "Human gammaherpesvirus 4",
    "Rhadinovirus humangamma8",     # KSHV
    "Human gammaherpesvirus 8",
    "Alphapolyomavirus quintihominis",  # JC virus
    "JC polyomavirus",
    "Betapolyomavirus hominis",         # BK virus
    "BK polyomavirus",
    "Lyssavirus rabies",            # rabies
    "Morbillivirus hominis",        # measles
    "Orthobornavirus bornaense",    # Borna virus encephalitis
    "Borna disease virus 1",        # NCBI taxonomy name for BoDV-1 in KrakenUniq
    "Human immunodeficiency virus 1",  # HIV-1 (CNS involvement/HIV encephalopathy)
    "Orthoflavivirus nilense",      # West Nile virus
    "West Nile virus",
    "Orthoflavivirus japonicum",    # Japanese encephalitis
    "Orthoflavivirus encephalitidis", # tick-borne encephalitis (TBE)
    "Tick-borne encephalitis virus",
    "Severe acute respiratory syndrome-related coronavirus",
    "Severe acute respiratory syndrome coronavirus 2",   # NCBI taxon name in KrakenUniq
    "SARS-CoV-2",
    # ── Bacteria ─────────────────────────────────────────────────────────────
    "Listeria monocytogenes",
    "Neisseria meningitidis",
    "Streptococcus pneumoniae",
    "Mycobacterium tuberculosis",
    "Mycobacterium tuberculosis complex",  # species-group node in KrakenUniq
    "Treponema pallidum",
    "Borreliella burgdorferi",
    "Borrelia burgdorferi",
    "Rickettsia prowazekii",
    "Rickettsia rickettsii",
    "Ehrlichia chaffeensis",
    "Bartonella henselae",
    "Brucella melitensis",
    "Tropheryma whipplei",
    "Orientia tsutsugamushi",
    "Coxiella burnetii",
    "Francisella tularensis",
    # ── Fungi ────────────────────────────────────────────────────────────────
    "Cryptococcus neoformans",
    "Cryptococcus gattii",
    "Aspergillus fumigatus",
    "Candida albicans",
    "Candida auris",
    "[Candida] auris",
    "Trichosporon asahii",
    "Histoplasma capsulatum",
    "Coccidioides immitis",
    "Coccidioides posadasii",
    "Pneumocystis jirovecii",
    # ── Parasites ────────────────────────────────────────────────────────────
    "Toxoplasma gondii",
    "Cryptosporidium hominis",
    "Cryptosporidium parvum",
    "Naegleria fowleri",
    "Balamuthia mandrillaris",
    "Acanthamoeba castellanii",
    "Acanthamoeba culbertsoni",
    "Trypanosoma brucei",
    "Trypanosoma cruzi",
    "Plasmodium falciparum",
    "Plasmodium vivax",
    "Plasmodium knowlesi",
    "Taenia solium",
    "Echinococcus granulosus",
    "Echinococcus multilocularis",
    "Angiostrongylus cantonensis",
}


# ─────────────────────────────────────────────────────────────────────────────
# Parser
# ─────────────────────────────────────────────────────────────────────────────
def detect_platform(header_lines):
    """
    Try to infer sequencing platform from the KrakenUniq CL header line.
    KrakenUniq itself is platform-agnostic, but the input filename sometimes
    contains hints (e.g. 'illumina', 'ont', 'nanopore').
    Returns 'Nanopore Metagenomics' or 'Illumina Metagenomics'.
    """
    for line in header_lines:
        lower = line.lower()
        if any(k in lower for k in ["illumina", "ilmn", "miseq", "nextseq",
                                     "hiseq", "novaseq", "_r1_", "_r2_"]):
            return "Illumina Metagenomics"
        if any(k in lower for k in ["nanopore", "ont", "fast5", "pod5",
                                     "pass_barcode", "guppy", "dorado"]):
            return "Nanopore Metagenomics"
    return "Nanopore Metagenomics"   # default


def assess(entry):
    name  = entry["name"]

    if name in SKIP_TAXA:                   return None
    if entry["rank"] not in RELEVANT_RANKS: return None
    # Require directly assigned reads (taxReads > 0).
    # taxReads=0 means reads belong entirely to child nodes → duplicate if included.
    if entry["taxReads"] == 0:              return None

    reads = entry["taxReads"]
    kmers = entry["kmers"]
    dup   = entry["dup"]
    cov   = entry["cov"]
    # Category from the taxonomy tree (kraken_tree.category_from_lineage()).
    kingdom = entry["kingdom"]

    is_contaminant    = name in LIKELY_CONTAMINANTS
    is_opportunist    = name in OPPORTUNISTIC_FLORA
    is_known_pathogen = name in KNOWN_PATHOGENS

    t_c = THRESHOLDS["CONVINCING"]
    t_n = THRESHOLDS["NOTABLE"]
    # dup cutoff lifted for viruses (see VIRAL_MAX_DUP comment above).
    max_dup_c = VIRAL_MAX_DUP if kingdom == kraken_tree.VIRUSES else t_c["maxDup"]
    max_dup_n = VIRAL_MAX_DUP if kingdom == kraken_tree.VIRUSES else t_n["maxDup"]

    if (reads >= t_c["minReads"] and kmers >= t_c["minKmers"]
            and dup <= max_dup_c and cov >= t_c["minCov"]):
        level = "CONVINCING"
    elif (reads >= t_n["minReads"] and kmers >= t_n["minKmers"]
          and dup <= max_dup_n and cov >= t_n["minCov"]):
        level = "NOTABLE"
    else:
        level = "LOW"

    # Pathogen bonus A: NOTABLE → CONVINCING only if dup is also within the
    # CONVINCING threshold. A known pathogen with a dirty dup signal (e.g.
    # Trichosporon with dup=1020) should not be silently promoted.
    if is_known_pathogen and level == "NOTABLE" and dup <= max_dup_c:
        level = "CONVINCING"

    # Pathogen bonus B: LOW → NOTABLE only with a minimum of reads, k-Mer
    # diversity, AND an acceptable dup value. reads>=3 (instead of >=5)
    # captures borderline cases with a very low but diagnostically confirmed
    # pathogen load (reference-cohort validation: A. fumigatus with 4 reads/
    # 30 k-mers in ITS-confirmed aspergillosis). The dup guard remains to
    # avoid false positives such as Trypanosoma (dup=234) or Cryptosporidium
    # (dup=864) -- for viruses the lifted max_dup_n applies.
    if (is_known_pathogen and level == "LOW"
            and reads >= 3 and kmers >= 10 and dup <= max_dup_n):
        level = "NOTABLE"

    # Contaminant cap: never promote a pure contaminant to CONVINCING.
    # Opportunistic flora can legitimately reach CONVINCING (e.g. S. intermedius
    # in a brain abscess with strong signal).
    if is_contaminant and level == "CONVINCING":
        level = "NOTABLE"

    if level == "LOW":
        return None

    if is_known_pathogen:
        comment = "Pathogen"
    elif is_opportunist:
        comment = "Commensal"
    elif is_contaminant:
        comment = "Contaminant"
    else:
        comment = ""

    return {
        **entry,
        "reads_display":  reads,
        "level":          level,
        "is_contaminant": is_contaminant,
        "is_opportunist": is_opportunist,
        "is_pathogen":    is_known_pathogen,
        "comment":        comment,
        "kingdom":        kingdom,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Reference database (prevalence across all previous KrakenUniq reports)
# ─────────────────────────────────────────────────────────────────────────────
_REF_DB_CACHE = None


def load_reference_db():
    """
    Loads the genus-level prevalence reference database
    (analysis/taxon_prevalence_genus_min10kb.csv, produced by
    scripts/build_reference_db.py). Deliberately based on the >=10kB cohort
    (n=437) rather than the complete database (n=1194): the same cohort the
    boxplot/z-score reference distributions (taxon_reference_stats.csv) come
    from -- very small/nearly empty reports (median report size < 10kB)
    would otherwise bias the prevalence estimate without adding resolution.

    Genus level rather than species/strain level: a single genus (e.g.
    Bradyrhizobium) can be spread across hundreds of individual strain/
    isolate names in the reference database, each of which is rare on its
    own -- even though the genus as a whole is a common background
    contaminant. Prioritisation in the report is therefore done at genus
    level; the exact species/strain name is still displayed.

    The result is cached per process, since the table has several thousand
    rows.
    """
    global _REF_DB_CACHE
    if _REF_DB_CACHE is not None:
        return _REF_DB_CACHE

    base = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(base, "..", "analysis", "taxon_prevalence_genus_min10kb.csv")
    meta_path = os.path.join(base, "..", "analysis", "reference_meta.json")

    n_total = 437
    if os.path.exists(meta_path):
        try:
            with open(meta_path, encoding="utf-8") as fh:
                n_total = json.load(fh).get("n_samples_min10kb", n_total)
        except (OSError, ValueError):
            pass

    db = {}
    if os.path.exists(csv_path):
        with open(csv_path, encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh):
                try:
                    db[row["genus_taxid"]] = (
                        int(row["n_samples_present"]),
                        float(row["pct_samples_present"]),
                    )
                except (KeyError, ValueError):
                    continue

    _REF_DB_CACHE = (db, n_total)
    return _REF_DB_CACHE


def attach_reference(hit, ref_db, n_total):
    """Enriches a hit with its genus-level reference prevalence (n / % of n_total)."""
    n_present, pct_present = ref_db.get(hit["genus_taxid"], (0, 0.0))
    hit["ref_n"] = n_present
    hit["ref_pct"] = pct_present
    hit["ref_total"] = n_total
    return hit


# Shared field set for taxon_reference_stats.csv AND genus_reference_stats.csv
# (identical schema, see build_reference_db.py: ref_stats_row()).
_TAXON_STATS_FIELDS = [
    "kmers_min", "kmers_p25", "kmers_p50", "kmers_p75", "kmers_max",
    "cov_min", "cov_p25", "cov_p50", "cov_p75", "cov_max",
    "dup_min", "dup_p25", "dup_p50", "dup_p75", "dup_max",
    "reads_log_mean", "reads_log_sd",
    "adj_kmers_log_mean", "adj_kmers_log_sd",
    "kmers_log_mean", "kmers_log_sd",
    "cov_log_mean", "cov_log_sd",
    "n_samples_present",
]


_POOLED_SD_CACHE = None


def load_pooled_reads_log_sd():
    """
    Loads the pooled standard deviation of log10(reads+1) across ALL hits
    of the >=10kB reference cohort (analysis/reference_meta.json). Serves as
    the fallback denominator for the z-score when a taxon occurs too rarely
    in the reference (n<3) to estimate its own stable SD.
    """
    global _POOLED_SD_CACHE
    if _POOLED_SD_CACHE is not None:
        return _POOLED_SD_CACHE
    base = os.path.dirname(os.path.abspath(__file__))
    meta_path = os.path.join(base, "..", "analysis", "reference_meta.json")
    sd = 1.0
    if os.path.exists(meta_path):
        try:
            with open(meta_path, encoding="utf-8") as fh:
                sd = json.load(fh).get("pooled_reads_log_sd", sd)
        except (OSError, ValueError):
            pass
    _POOLED_SD_CACHE = sd
    return sd


def compute_z_score(value, box_stats, pooled_sd, min_n_for_own_sd=1,
                     mean_key="adj_kmers_log_mean", sd_key="adj_kmers_log_sd",
                     transform="add1"):
    """
    z-score: how many standard deviations (log10 scale) is this hit above
    the typical background noise for this taxon in the >=10kB reference
    cohort?

    By default (transform="add1", value = coverage-adjusted k-mer count, see
    kraken_tree.adjusted_kmers()): k-mers separate true pathogens from
    background noise better than raw reads (see
    analysis/threshold_analysis_report.md, AUC 0.873 vs. 0.755), and the
    dup adjustment weights broad, credible genome coverage (many distinct
    positions hit once each) higher than a narrow, highly duplicated
    signature (few positions, hit often) -- except for viruses, where high
    dup values arise from genuine deep coverage of small genomes and are
    therefore not penalised. log10(value+1), because value as count data
    (k-mers, reads) can be 0.

    transform="direct": for quantities that are always tiny fractions
    (coverage fraction) -- log10(value+1) would collapse to ~0 for every
    value and destroy the signal, so the unshifted log10(max(value,1e-9)) is
    used instead (identical to the reference distribution, see
    build_reference_db.py log_mean_sd_direct()).

    The reference mean/SD (box_stats[mean_key]/[sd_key], see
    build_reference_db.py: ref_stats_row()) are computed over the ENTIRE
    >=10kB cohort, not only over the cases with a detection -- samples
    without a detection explicitly enter with 0. As a result, the
    taxon-specific SD is already well conditioned from a single reference
    occurrence onwards (for rare taxa it is dominated by the many "absence"
    zeros), which is why the taxon's own SD is used by default from
    n_ref>=1:

      - Taxon with >=1 reference occurrence AND SD>0: own mean/own SD.
      - Taxon never seen in the reference (n=0): background is assumed to
        be "silence" (mean=0), with the pooled SD (across all taxa) as the
        denominator.
    """
    log_val = math.log10(value + 1) if transform == "add1" else math.log10(max(value, 1e-9))
    if pooled_sd <= 0:
        pooled_sd = 1.0

    n_ref = int(box_stats["n_samples_present"]) if box_stats else 0
    if box_stats and n_ref >= min_n_for_own_sd and box_stats[sd_key] > 0:
        mean, sd = box_stats[mean_key], box_stats[sd_key]
    elif box_stats and n_ref >= 1:
        mean, sd = box_stats[mean_key], pooled_sd
    else:
        mean, sd = 0.0, pooled_sd

    return (log_val - mean) / sd


_GENUS_STATS_CACHE = None


def load_genus_stats():
    """
    Loads the genus-level distribution statistics
    (analysis/genus_reference_stats.csv, >=10kB cohort): reads/k-mers summed
    per sample over all species/strains of a genus (coverage as max, dup
    read-weighted), then aggregated across the cohort. Basis for the
    genus-level boxplots/z-scores in the report (see group_by_genus()).
    """
    global _GENUS_STATS_CACHE
    if _GENUS_STATS_CACHE is not None:
        return _GENUS_STATS_CACHE

    base = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(base, "..", "analysis", "genus_reference_stats.csv")

    db = {}
    if os.path.exists(csv_path):
        with open(csv_path, encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh):
                try:
                    db[row["genus_taxid"]] = {f: float(row[f]) for f in _TAXON_STATS_FIELDS}
                except (KeyError, ValueError):
                    continue

    _GENUS_STATS_CACHE = db
    return db


_POOLED_SD_META_CACHE = {}


def _load_pooled_sd_genus(meta_key):
    """
    Pooled SD (log scale) of a genus-level metric across the >=10kB
    reference cohort (analysis/reference_meta.json) -- fallback denominator
    for compute_z_score() when a genus occurs too rarely in the reference to
    estimate its own stable SD.
    """
    if meta_key in _POOLED_SD_META_CACHE:
        return _POOLED_SD_META_CACHE[meta_key]
    base = os.path.dirname(os.path.abspath(__file__))
    meta_path = os.path.join(base, "..", "analysis", "reference_meta.json")
    sd = 1.0
    if os.path.exists(meta_path):
        try:
            with open(meta_path, encoding="utf-8") as fh:
                sd = json.load(fh).get(meta_key, sd)
        except (OSError, ValueError):
            pass
    _POOLED_SD_META_CACHE[meta_key] = sd
    return sd


def load_pooled_adj_kmers_log_sd_genus():
    return _load_pooled_sd_genus("pooled_adj_kmers_log_sd_genus")


def load_pooled_kmers_log_sd_genus():
    return _load_pooled_sd_genus("pooled_kmers_log_sd_genus")


def load_pooled_cov_log_sd_genus():
    return _load_pooled_sd_genus("pooled_cov_log_sd_genus")


def group_by_genus(hits):
    """
    Aggregates species/strain-level hits (already filtered by assess())
    into one row per genus:
      - reads_display / kmers: sum over all species/strains of the genus in
        this report
      - cov: maximum (coverage depends on genome size and therefore cannot
        be meaningfully summed; the max value shows the strongest single
        piece of evidence)
      - dup: read-weighted mean
      - adj_kmers: sum of the members' coverage-adjusted k-mer counts
        (kraken_tree.adjusted_kmers(), computed per member with its own dup
        value, then summed) -- basis of the z-score, see
        attach_genus_stats()/compute_z_score()
      - comment: highest priority among the members
        (Pathogen > Commensal > Contaminant > "")
      - top3: the 3 species/strains with the most k-mers (for the organism
        sub-line in the report)
    z-score/boxplot reference data (box_stats, z_score) are NOT set here --
    call attach_genus_stats() after grouping for that.
    """
    groups = defaultdict(list)
    for h in hits:
        groups[h["genus_taxid"]].append(h)

    comment_priority = {"Pathogen": 0, "Commensal": 1, "Contaminant": 2, "": 3}

    genus_rows = []
    for genus_taxid, members in groups.items():
        top3 = sorted(members, key=lambda h: -h["kmers"])[:3]
        reads_sum = sum(h["taxReads"] for h in members)
        kmers_sum = sum(h["kmers"] for h in members)
        cov_max = max(h["cov"] for h in members)
        dup_avg = (sum(h["dup"] * h["taxReads"] for h in members) / reads_sum
                   if reads_sum else 0.0)
        kingdom = members[0]["kingdom"]
        adj_kmers_sum = sum(
            kraken_tree.adjusted_kmers(h["kmers"], h["dup"], kingdom) for h in members)
        comment = min((h["comment"] for h in members), key=lambda c: comment_priority[c])
        rep = members[0]
        genus_rows.append(dict(
            genus_taxid=genus_taxid,
            name=rep["genus_name"],
            kingdom=rep["kingdom"],
            reads_display=reads_sum,
            kmers=kmers_sum,
            adj_kmers=adj_kmers_sum,
            cov=cov_max,
            dup=dup_avg,
            comment=comment,
            ref_n=rep["ref_n"], ref_pct=rep["ref_pct"], ref_total=rep["ref_total"],
            top3=top3,
            n_species=len(members),
            members=members,
        ))
    return genus_rows


def attach_genus_stats(genus_row, genus_stats, pooled_sd_genus):
    """
    Sets z_score (main score, coverage-adjusted k-mers, see
    compute_z_score()) as well as the separate kmers_z_score/cov_z_score --
    same methodology, but on the raw k-mers or the coverage of this genus,
    shown below the respective boxplot in the top-hits table
    (make_hit_table()), so that it becomes visible how strongly k-mer count
    and coverage EACH deviate from the background, not only the combined
    dup-adjusted metric.
    """
    genus_row["box_stats"] = genus_stats.get(genus_row["genus_taxid"])
    genus_row["z_score"] = compute_z_score(
        genus_row["adj_kmers"], genus_row["box_stats"], pooled_sd_genus)
    genus_row["kmers_z_score"] = compute_z_score(
        genus_row["kmers"], genus_row["box_stats"], load_pooled_kmers_log_sd_genus(),
        mean_key="kmers_log_mean", sd_key="kmers_log_sd")
    genus_row["cov_z_score"] = compute_z_score(
        genus_row["cov"], genus_row["box_stats"], load_pooled_cov_log_sd_genus(),
        mean_key="cov_log_mean", sd_key="cov_log_sd", transform="direct")
    return genus_row


def parse_report(text, ref_db=None, ref_total=None):
    header_lines = [l for l in text.split("\n") if l.startswith("#")]
    platform = detect_platform(header_lines)

    entries = kraken_tree.parse_kraken_rows(text)

    total_classified = total_unclassified = total_human = 0
    pct_classified   = pct_unclassified   = pct_human   = 0.0
    for e in entries:
        if e["taxID"] == "0":
            total_unclassified = e["reads"]
            pct_unclassified   = e["pct"]
        if e["taxID"] == "1":
            total_classified   = e["reads"]
            pct_classified     = e["pct"]
        if e["taxID"] == "9606":
            total_human        = e["reads"]
            pct_human          = e["pct"]

    hits = [h for h in (assess(e) for e in entries) if h is not None]

    if ref_db is None or ref_total is None:
        ref_db, ref_total = load_reference_db()
    hits = [attach_reference(h, ref_db, ref_total) for h in hits]
    # z-score/boxplot reference data are needed at genus level (see
    # group_by_genus()/attach_genus_stats() in select_top_hits()), no longer
    # per individual species/strain hit.

    total_all = total_classified + total_unclassified

    return {
        "platform":           platform,
        "total_all":          total_all,
        "total_classified":   total_classified,
        "total_unclassified": total_unclassified,
        "total_human":        total_human,
        "pct_classified":     pct_classified,
        "pct_unclassified":   pct_unclassified,
        "pct_human":          pct_human,
        "ref_total":          ref_total,
        "hits":               hits,
    }


# ─────────────────────────────────────────────────────────────────────────────
# PDF helpers
# ─────────────────────────────────────────────────────────────────────────────
def fmt_num(n):
    return f"{int(n):,}"

def fmt_cov(cov):
    """Return coverage as percentage string.
    >= 0.1%  → 1 decimal place (e.g. 1.1%)
    < 0.1%   → scientific notation (e.g. 1.8×10⁻²%)
    == 0     → "0%"
    """
    if cov == 0:
        return "0%"
    pct = cov * 100
    if pct >= 0.1:
        return f"{pct:.1f}%"
    # Below 0.1%: use scientific notation to avoid "0.0%"
    s = f"{pct:.1e}"
    m, e = s.split("e")
    return f"{m}×10<super>{int(e)}</super>%"

def fmt_dup(dup):
    """1 decimal place; from >= 100 rounded to an integer (column width)."""
    if dup >= 100:
        return f"{dup:.0f}"
    return f"{dup:.1f}"


def comment_text_color(comment):
    """Text colour of the comment column per category."""
    return {
        "Pathogen":    C_RED_TEXT,
        "Commensal":   C_PURPLE_TEXT,
        "Contaminant": C_GREY_TEXT,
    }.get(comment, C_MID)




def zscore_tier_props(z):
    """
    Background colour of the z-score cell in the top-hits table.
    >= Z_THRESHOLD_HIGH     -> strong deviation from background noise (reddish)
    >= Z_THRESHOLD_ELEVATED -> elevated deviation (orange)
    otherwise               -> within the usual noise (neutral)
    """
    if z >= Z_THRESHOLD_HIGH:
        return C_RED_FILL, C_RED_TEXT
    if z >= Z_THRESHOLD_ELEVATED:
        return C_ORANGE_FILL, C_ORANGE_TEXT
    return None, C_MID


def fmt_ref(hit):
    return f"{hit['ref_n']} ({hit['ref_pct']:.1f}%)"


def fmt_zscore(z):
    return f"{z:+.1f}"


def fmt_num_short(n):
    """Compact number for the top-3 species sub-line, e.g. 8212 -> '8.2k'."""
    if n >= 1000:
        s = f"{n/1000:.1f}k"
        return s.replace(".0k", "k")
    return str(n)


def fmt_top3_species(genus_row):
    """
    Formats the 3 species/strains of a genus row with the most k-mers as a
    compact, '·'-separated list for the organism sub-line in the report,
    e.g. "coli (8.2k) · albertii (540) · fergusonii (12)". The genus name is
    removed from the species name if it occurs as a prefix (saves space,
    since the genus is already shown in the main line).
    """
    genus = genus_row["name"]
    parts = []
    for h in genus_row["top3"]:
        name = h["name"]
        if name.lower().startswith(genus.lower() + " "):
            label = name[len(genus) + 1:]
        else:
            label = name
        parts.append(f"{label} ({fmt_num_short(h['kmers'])})")
    return " · ".join(parts)


# ── Mini boxplots (reference distribution k-mers/coverage, log scale) ─────
C_BOX_FILL = HexColor("#E5E5EA")   # pale grey (reference distribution)
C_BOX_LINE = HexColor("#AEAEB2")
C_BOX_DOT  = HexColor("#FF3B30")   # red dot: current case


def _log_frac(v, lo, hi):
    """Position of v on a log10 scale between lo and hi, as 0..1."""
    eps = 1e-12
    v, lo, hi = max(v, eps), max(lo, eps), max(hi, eps)
    if hi <= lo:
        return 0.5
    frac = (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
    return min(max(frac, 0.0), 1.0)


def make_boxplot(current_value, stats, prefix, width, height):
    """
    Horizontal mini boxplot on a log10 scale (saves row height compared to
    the vertical variant):
      - grey box: reference interquartile range (p25-p75), >=10kB cohort
      - grey line: reference median (p50)
      - grey whiskers: reference min/max
      - red dot: value of the current case
    Returns None if no reference data exist for this taxon (e.g. new/only
    observed outside the >=10kB cohort).
    """
    if not stats:
        return None
    vmin = stats[f"{prefix}_min"]
    p25  = stats[f"{prefix}_p25"]
    p50  = stats[f"{prefix}_p50"]
    p75  = stats[f"{prefix}_p75"]
    vmax = stats[f"{prefix}_max"]

    lo = min(vmin, current_value)
    hi = max(vmax, current_value)

    pad = 1.5 * mm
    plot_w = width - 2 * pad

    def x(v):
        return pad + _log_frac(v, lo, hi) * plot_w

    d = Drawing(width, height)
    cy = height / 2.0
    box_h = height * 0.55

    d.add(Line(x(vmin), cy, x(vmax), cy, strokeColor=C_BOX_LINE, strokeWidth=0.6))
    x0, x1 = x(p25), x(p75)
    d.add(Rect(min(x0, x1), cy - box_h / 2, max(abs(x1 - x0), 0.3), box_h,
               fillColor=C_BOX_FILL, strokeColor=C_BOX_LINE, strokeWidth=0.5))
    xm = x(p50)
    d.add(Line(xm, cy - box_h / 2, xm, cy + box_h / 2,
               strokeColor=C_BOX_LINE, strokeWidth=0.9))
    xv = x(current_value)
    d.add(Circle(xv, cy, 1.15, fillColor=C_BOX_DOT, strokeColor=None))
    return d


def S(name, **kw):
    return ParagraphStyle(name, **kw)

def P(text, style, **kw):
    if kw:
        style = ParagraphStyle("_", parent=style, **kw)
    return Paragraph(text, style)

def HR(sb=3*mm, sa=3*mm):
    return HRFlowable(width=CW, thickness=0.5, color=C_RULE,
                      spaceBefore=sb, spaceAfter=sa)


# ─────────────────────────────────────────────────────────────────────────────
# Header / Footer
# ─────────────────────────────────────────────────────────────────────────────
def draw_chrome(c, doc, sample_name, total_pages):
    w, h = A4
    page_num = doc.page

    # Top rule
    c.setStrokeColor(C_RULE_DARK)
    c.setLineWidth(0.5)
    c.line(ML, h - MT + 5*mm, w - MR, h - MT + 5*mm)

    c.setFillColor(C_NAVY)
    c.setFont("Helvetica-Bold", 8.5)
    c.drawString(ML, h - MT + 8*mm, "Metagenomic Analysis")

    c.setFillColor(C_MID)
    c.setFont("Helvetica", 8)
    c.drawRightString(w - MR, h - MT + 8*mm, sample_name)

    # Bottom rule + footer
    c.setStrokeColor(C_RULE_DARK)
    c.setLineWidth(0.5)
    c.line(ML, MB - 3*mm, w - MR, MB - 3*mm)

    c.setFillColor(C_LIGHT)
    c.setFont("Helvetica", 7)
    c.drawString(ML, MB - 7*mm,
                 "Institute of Neuropathology · University Hospital Münster")
    c.drawRightString(w - MR, MB - 7*mm,
                      f"Page {page_num} of {total_pages}")


# ─────────────────────────────────────────────────────────────────────────────
# Story builder
# ─────────────────────────────────────────────────────────────────────────────
def priority_key(h):
    """
    Sort descending by z-score (standard deviations from the taxon-specific
    background noise, >=10kB reference cohort).

    Note: an earlier version additionally sorted by genus-prevalence tier
    (< 1% / < 5% / >= 5%) BEFORE the z-score, because a pure z-score tends
    to rank even harmless commensals high in heavily polymicrobial (e.g.
    post-mortem) samples (practically every species there is strongly
    elevated relative to its own low background). The tier pre-sort was
    removed on request; the only remaining safeguard against exactly this
    scenario is the forced inclusion of curated pathogens/commensals in
    select_top_hits() (with "Pathogen" taking precedence over "Commensal"
    when the budget is limited). If needed, re-validate against the 63-case
    cohort with scripts/build_reference_cohort_summary.py.
    """
    return -h["z_score"]


def select_top_hits(hits):
    """
    Groups the species/strain-level hits of a report into genus rows
    (group_by_genus() + attach_genus_stats(), see there) and selects the
    top-hits table (page 1) from them. Single implementation of this
    selection logic -- used both by the PDF report (build_story) and by the
    reference-cohort validation (scripts/build_reference_cohort_summary.py),
    so that the two never diverge.

    Prioritised overall list: descending by z-score (see priority_key()),
    now computed at genus level (reads/k-mers summed over all species/
    strains of the genus). The top N reliably fill page 1 even when there
    are only a few notable hits (topped up with the next-highest z-scores)
    -- the rest appear with the same metrics/boxplots in the supplement,
    split into "additional rare" and "common" findings (there still based
    on the genus-prevalence tier, < 5% vs. >= 5%).

    Forced inclusion regardless of z-score: genus rows containing a curated
    pathogen/commensal (comment column) with non-trivial read support are
    never pushed into the supplement solely because of a lower z-score.
    Reference-cohort validation (63 cases with a confirmed pathogen):
    without this rule, e.g. S. aureus (43 reads), K. pneumoniae (318) and
    A. fumigatus (133,475) ended up outside the top 15 despite a dominant
    signal.

    The forced list is capped at FORCED_MAX_N (strongest reads first):
    heavily polymicrobial samples (e.g. autopsy/post-mortem brain samples
    with extensive bacterial overgrowth) can contain 20+ commensal genera
    with >= FORCED_MIN_READS. Without a cap they would fill the entire
    top-hits table and displace the actually rare, diagnostically decisive
    finding. At least MIN_RARITY_SLOTS slots are therefore always reserved
    for pure z-score prioritisation.

    Returns: (top_hits, rest_rare, rest_common) -- lists of genus rows.
    """
    genus_stats = load_genus_stats()
    pooled_sd_genus = load_pooled_adj_kmers_log_sd_genus()
    rows = [attach_genus_stats(r, genus_stats, pooled_sd_genus)
            for r in group_by_genus(hits)]

    # Inclusion criteria ONLY for the top hits (dup cutoff, see
    # TOP_HITS_MAX_DUP, kingdom "Viruses" exempt; curated contaminants, see
    # LIKELY_CONTAMINANTS). Hits that fail them remain in the pool for
    # rest_rare/rest_common (supplement), but can neither be forced nor
    # selected into top_hits via z-score. The contaminant exclusion is
    # deliberately hard (not just a bonus/penalty on the z-score): since the
    # reference distribution counts absence as 0 (see change history in
    # reference_cohort_filter_strategy.md), "rare in the reference" is the
    # normal case for practically every genus, not specific to true
    # pathogens -- curated contaminants would otherwise reach a high z-score
    # just as often as rare genuine findings and dilute page 1. A
    # Contaminant tag at genus level only arises if NO member of that genus
    # is classified as Pathogen/Commensal (group_by_genus() takes the
    # highest priority among the members), so a "contaminated" finding with
    # a genuine pathogen species in the same genus is not excluded.
    eligible = [r for r in rows
                if (r["kingdom"] == kraken_tree.VIRUSES or r["dup"] < TOP_HITS_MAX_DUP)
                and r["comment"] != "Contaminant"]

    # Curated pathogens take precedence over commensals within the forced
    # budget: in heavily polymicrobial samples, commensal genera (e.g.
    # post-mortem E. coli overgrowth, 13,529 reads) can occupy the entire
    # cap purely through raw read count and displace a genuine but
    # read-weaker curated genus (e.g. BoDV-1, 92 reads).
    forced_pathogens = sorted(
        [r for r in eligible if r["comment"] == "Pathogen"
         and r["reads_display"] >= FORCED_MIN_READS["Pathogen"]],
        key=lambda r: -r["reads_display"],
    )
    forced_commensals = sorted(
        [r for r in eligible if r["comment"] == "Commensal"
         and r["reads_display"] >= FORCED_MIN_READS["Commensal"]],
        key=lambda r: -r["reads_display"],
    )
    forced_cap = min(FORCED_MAX_N, TOP_N_HITS - MIN_RARITY_SLOTS)
    forced = forced_pathogens[:forced_cap]
    forced += forced_commensals[:max(forced_cap - len(forced), 0)]
    forced_ids = {id(r) for r in forced}
    eligible_sorted = sorted(eligible, key=priority_key)
    fill_n     = max(TOP_N_HITS - len(forced), 0)
    filler     = [r for r in eligible_sorted if id(r) not in forced_ids][:fill_n]
    top_hits   = sorted(forced + filler, key=priority_key)
    top_ids    = {id(r) for r in top_hits}
    all_sorted = sorted(rows, key=priority_key)
    rest       = [r for r in all_sorted if id(r) not in top_ids]
    rest_rare   = [r for r in rest if r["ref_pct"] < RARE_THRESHOLD_5PCT]
    rest_common = [r for r in rest if r["ref_pct"] >= RARE_THRESHOLD_5PCT]
    return top_hits, rest_rare, rest_common


def build_story(sample_name, date, parsed, platform_override=None):
    platform = platform_override or parsed["platform"]
    hits     = parsed["hits"]

    # Styles
    st_title     = S("title",   fontName="Helvetica-Bold", fontSize=17,
                     textColor=C_NAVY, spaceAfter=0.5*mm, leading=20, leftIndent=-6)
    st_subtitle  = S("sub",     fontName="Helvetica",      fontSize=9,
                     textColor=C_MID,  spaceAfter=2*mm, leading=12, leftIndent=-6)
    st_section   = S("sec",     fontName="Helvetica-Bold", fontSize=12,
                     textColor=C_NAVY, spaceBefore=5*mm, spaceAfter=2.5*mm, leading=15, leftIndent=-6)
    st_meta_lbl  = S("mlbl",    fontName="Helvetica-Bold", fontSize=6.5,
                     textColor=C_LIGHT, spaceAfter=0.5*mm, leading=8, letterSpacing=0.4)
    st_meta_val  = S("mval",    fontName="Helvetica",      fontSize=8,
                     textColor=C_NAVY,  spaceAfter=1.5*mm, leading=10)
    st_body      = S("body",    fontName="Helvetica",      fontSize=9.5,
                     textColor=C_MID,   leading=14, spaceAfter=3*mm)
    st_tbl_hdr   = S("th",      fontName="Helvetica-Bold", fontSize=7.5,
                     textColor=white,   alignment=TA_LEFT,   leading=9)
    st_tbl_hdr_r = S("thr",     fontName="Helvetica-Bold", fontSize=7.5,
                     textColor=white,   alignment=TA_RIGHT,  leading=9)
    st_tbl_hdr_c = S("thc",     fontName="Helvetica-Bold", fontSize=7.5,
                     textColor=white,   alignment=TA_CENTER, leading=9)
    st_cell      = S("tc",      fontName="Helvetica",      fontSize=8,
                     textColor=C_NAVY,  alignment=TA_LEFT,   leading=10)
    st_cell_i    = S("tci",     fontName="Helvetica-Oblique", fontSize=8,
                     textColor=C_NAVY,  alignment=TA_LEFT,   leading=10)
    st_cell_r    = S("tcr",     fontName="Helvetica",      fontSize=8,
                     textColor=C_NAVY,  alignment=TA_RIGHT,  leading=10)
    st_cell_mid  = S("tcm",     fontName="Helvetica",      fontSize=8,
                     textColor=C_MID,   alignment=TA_LEFT,   leading=10)
    st_leg_term  = S("lt",      fontName="Helvetica-Bold", fontSize=8,
                     textColor=C_NAVY,  leading=11)
    st_leg_desc  = S("ld",      fontName="Helvetica",      fontSize=8,
                     textColor=C_MID,   leading=11)

    story = []

    # ── TITLE (slim) ─────────────────────────────────────────────────────
    story.append(P("Metagenomic Analysis", st_title))
    story.append(P(f"{sample_name} · {date}", st_subtitle))
    story.append(HR(0.5*mm, 2*mm))

    # ── META (single compact row) ─────────────────────────────────────────
    # Truncate very long sample names in the meta row to avoid line breaks
    def trunc(s, max_len=34):
        return s if len(s) <= max_len else s[:max_len - 1] + "…"

    meta_cols = [
        ("SAMPLE",    trunc(sample_name), 70*mm),
        ("DATE",      date,               26*mm),
        ("METHOD",    platform,           42*mm),
        ("REFERENCE", f"n = {parsed.get('ref_total', 437)} cases", 32*mm),
    ]
    meta_cells = []
    meta_widths = []
    for lbl, val, w in meta_cols:
        meta_cells.append([P(lbl, st_meta_lbl), P(val, st_meta_val)])
        meta_widths.append(w)

    mt = Table([meta_cells], colWidths=meta_widths)
    mt.setStyle(TableStyle([
        ("VALIGN",       (0,0), (-1,-1), "TOP"),
        ("LEFTPADDING",  (0,0), (-1,-1), 0),
        ("RIGHTPADDING", (0,0), (-1,-1), 4*mm),
        ("TOPPADDING",   (0,0), (-1,-1), 0),
        ("BOTTOMPADDING",(0,0), (-1,-1), 0),
    ]))
    story.append(mt)
    story.append(HR(1.5*mm, 2*mm))

    # ── STAT TILES (slimmer) ────────────────────────────────────────────
    tile_data = [
        (C_TILE_BLUE,
         "TOTAL READS",
         fmt_num(parsed["total_all"]),
         ""),
        (C_TILE_GREEN,
         "CLASSIFIED",
         f"{parsed['pct_classified']:.1f}%",
         fmt_num(parsed["total_classified"]) + " reads"),
        (HexColor("#30A46C"),
         "HUMAN",
         f"{parsed['pct_human']:.1f}%",
         fmt_num(parsed["total_human"]) + " reads"),
        (C_TILE_GREY,
         "UNCLASSIFIED",
         f"{parsed['pct_unclassified']:.1f}%",
         fmt_num(parsed["total_unclassified"]) + " reads"),
    ]

    def tile_cell(colour, label, big, sub):
        white_lbl  = S("_tl", fontName="Helvetica-Bold",  fontSize=6,
                        textColor=HexColor("#FFFFFFBB"), leading=7,
                        letterSpacing=0.5, spaceAfter=1*mm)
        white_big  = S("_tb", fontName="Helvetica-Bold",  fontSize=14,
                        textColor=white, leading=17, spaceAfter=0)
        white_sub  = S("_ts", fontName="Helvetica",       fontSize=6.5,
                        textColor=HexColor("#FFFFFFAA"), leading=8, spaceAfter=0)
        items = [P(label, white_lbl), P(big, white_big)]
        if sub:
            items.append(P(sub, white_sub))
        return items

    tile_cells  = [[tile_cell(*td) for td in tile_data]]
    tile_col_w  = CW / 4 - 1*mm
    tile_tbl    = Table(tile_cells, colWidths=[tile_col_w]*4, rowHeights=[13*mm])
    ts_cmds = [
        ("VALIGN",       (0,0), (-1,-1), "MIDDLE"),
        ("LEFTPADDING",  (0,0), (-1,-1), 3.5*mm),
        ("RIGHTPADDING", (0,0), (-1,-1), 2.5*mm),
        ("TOPPADDING",   (0,0), (-1,-1), 2*mm),
        ("BOTTOMPADDING",(0,0), (-1,-1), 2*mm),
        ("ROUNDEDCORNERS", [4]),
    ]
    for i, (colour, *_) in enumerate(tile_data):
        ts_cmds.append(("BACKGROUND", (i,0), (i,0), colour))
    tile_tbl.setStyle(TableStyle(ts_cmds))
    story.append(tile_tbl)
    story.append(HR(3*mm, 1*mm))

    # ── FINDINGS TABLES ───────────────────────────────────────────────────
    kingdoms = kraken_tree.CATEGORY_ORDER
    cw = list(COL_W.values())

    top_hits, rest_rare, rest_common = select_top_hits(hits)

    # Column indices in COL_W (order matters!) for SPAN/ALIGN rules
    IDX = {name: i for i, name in enumerate(COL_W.keys())}

    def stacked_metric_cell(value_para, chart, inner_w, z_para=None):
        """Number stacked above a mini boxplot (above an optional z-score
        line), as one compact cell."""
        rows = [[value_para], [chart]]
        if z_para is not None:
            rows.append([z_para])
        t = Table(rows, colWidths=[inner_w])
        t.setStyle(TableStyle([
            ("ALIGN",        (0,0), (-1,-1), "CENTER"),
            ("VALIGN",       (0,0), (-1,-1), "MIDDLE"),
            ("TOPPADDING",   (0,0), (-1,-1), 0.1*mm),
            ("BOTTOMPADDING",(0,0), (-1,-1), 0.1*mm),
            ("LEFTPADDING",  (0,0), (-1,-1), 0),
            ("RIGHTPADDING", (0,0), (-1,-1), 0),
        ]))
        return t

    def make_hit_table(group):
        """Top-hits table: full metrics + mini boxplots + reference + comment."""
        tp = bp = 1.5
        lp = rp = 2.5
        fs = 8.0
        fs_hdr = 7.5

        s_hdr   = S("_sh",  fontName="Helvetica-Bold",    fontSize=fs_hdr,
                    textColor=white,  alignment=TA_LEFT,   leading=fs_hdr+1.5)
        s_hdr_r = S("_shr", fontName="Helvetica-Bold",    fontSize=fs_hdr,
                    textColor=white,  alignment=TA_RIGHT,  leading=fs_hdr+1.5)
        s_hdr_c = S("_shc", fontName="Helvetica-Bold",    fontSize=fs_hdr,
                    textColor=white,  alignment=TA_CENTER, leading=fs_hdr+1.5)
        s_ci    = S("_sci", fontName="Helvetica-Oblique", fontSize=fs,
                    textColor=C_NAVY, alignment=TA_LEFT,   leading=fs+1.5)
        s_cr    = S("_scr", fontName="Helvetica",         fontSize=fs,
                    textColor=C_NAVY, alignment=TA_RIGHT,  leading=fs+1.5)
        s_cc    = S("_scc", fontName="Helvetica",         fontSize=fs,
                    textColor=C_NAVY, alignment=TA_CENTER, leading=fs+0.5)
        s_dash  = S("_sda", fontName="Helvetica",         fontSize=fs,
                    textColor=C_LIGHT, alignment=TA_CENTER, leading=fs+0.5)
        s_metz  = S("_smz", fontName="Helvetica",         fontSize=6,
                    textColor=HexColor("#AEAEB2"), alignment=TA_CENTER, leading=6.5)

        hdr = [
            P("Organism",    s_hdr),
            P("Reads",       s_hdr_r),
            P("k-mers",      s_hdr_c),
            P("dup",         s_hdr_c),
            P("Cov. (%)",    s_hdr_c),
            P("z-score",     s_hdr_c),
            P("Comment",     s_hdr),
        ]
        rows = [hdr]
        for hit in group:
            _, z_txt_col = zscore_tier_props(hit["z_score"])
            z_st = S("_rf", fontName="Helvetica-Bold", fontSize=fs,
                     textColor=z_txt_col, alignment=TA_CENTER, leading=fs+0.5)
            cm_st = S("_cm", fontName="Helvetica", fontSize=fs,
                      textColor=comment_text_color(hit["comment"]),
                      alignment=TA_LEFT, leading=fs+1.5)

            kbox = make_boxplot(hit["kmers"], hit.get("box_stats"), "kmers", KMERS_INNER_W, BOX_H)
            dbox = make_boxplot(hit["dup"],   hit.get("box_stats"), "dup",   DUP_INNER_W,   BOX_H)
            cbox = make_boxplot(hit["cov"],   hit.get("box_stats"), "cov",   COV_INNER_W,   BOX_H)
            kmers_cell = stacked_metric_cell(
                P(str(hit["kmers"]), s_cc), kbox or P("–", s_dash), KMERS_INNER_W,
                P(f"z {fmt_zscore(hit['kmers_z_score'])}", s_metz))
            # dup has no z-score of its own (see legend), but gets an empty
            # third line with the same line height as kmers_z/cov_z
            # (instead of None) -- otherwise the 2-line dup cell would be
            # vertically centred within the row height (set to 3 lines by
            # kmers/Cov.), and its boxplot would sit visibly lower than the
            # other two boxplots.
            dup_cell = stacked_metric_cell(
                P(fmt_dup(hit["dup"]), s_cc), dbox or P("–", s_dash), DUP_INNER_W,
                P("&nbsp;", s_metz))
            cov_cell = stacked_metric_cell(
                P(fmt_cov(hit["cov"]), s_cc), cbox or P("–", s_dash), COV_INNER_W,
                P(f"z {fmt_zscore(hit['cov_z_score'])}", s_metz))
            z_cell = P(
                f"{fmt_zscore(hit['z_score'])}<br/>"
                f"<font size=6 color='#AEAEB2'>{fmt_ref(hit)}</font>", z_st)

            rows.append([
                P(f"<i>{hit['name']}</i><br/>"
                  f"<font size=6 color='#AEAEB2'>{fmt_top3_species(hit)}</font>", s_ci),
                P(str(hit["reads_display"]), s_cr),
                kmers_cell,
                dup_cell,
                cov_cell,
                z_cell,
                P(hit["comment"],            cm_st),
            ])

        from reportlab.lib.units import mm as _mm
        tbl = Table(rows, colWidths=cw, repeatRows=1)
        style_cmds = [
            ("BACKGROUND",    (0,0), (-1,0),  C_TABLE_HEAD),
            ("TOPPADDING",    (0,0), (-1,0),  tp*_mm),
            ("BOTTOMPADDING", (0,0), (-1,0),  bp*_mm),
            ("TOPPADDING",    (0,1), (-1,-1), tp*_mm),
            ("BOTTOMPADDING", (0,1), (-1,-1), bp*_mm),
            ("LEFTPADDING",   (0,0), (-1,-1), lp*_mm),
            ("RIGHTPADDING",  (0,0), (-1,-1), rp*_mm),
            ("ROWBACKGROUNDS",(0,1), (-1,-1), [C_ROW_EVEN, C_ROW_ODD]),
            ("LINEBELOW",     (0,0), (-1,-1), 0.4, C_RULE),
            ("VALIGN",        (0,0), (-1,-1), "MIDDLE"),
            ("ALIGN",         (IDX["reads"],0),   (IDX["reads"],-1),   "RIGHT"),
            ("ALIGN",         (IDX["dup"],0),     (IDX["dup"],-1),     "CENTER"),
            ("ALIGN",         (IDX["kmers"],0),   (IDX["kmers"],-1),   "CENTER"),
            ("ALIGN",         (IDX["cov"],0),     (IDX["cov"],-1),     "CENTER"),
            ("ALIGN",         (IDX["zscore"],0),  (IDX["zscore"],-1),  "CENTER"),
        ]
        for ri, hit in enumerate(group, start=1):
            z_fill, _ = zscore_tier_props(hit["z_score"])
            if z_fill is not None:
                style_cmds.append(("BACKGROUND", (IDX["zscore"], ri), (IDX["zscore"], ri), z_fill))
        tbl.setStyle(TableStyle(style_cmds))
        return tbl

    st_kg_lbl = S("_kgl", fontName="Helvetica-Bold", fontSize=8,
                  textColor=C_MID, spaceBefore=2*mm, spaceAfter=1.5*mm,
                  leading=10, leftIndent=-6)

    def kingdom_grouped_tables(group_all, sort_key):
        """Renders group_all as one make_hit_table() per kingdom (order see
        `kingdoms`), with a kingdom sub-heading -- shared by top hits and
        supplement, so that both are structured and illustrated
        (boxplots/z-scores) identically."""
        for kingdom in kingdoms:
            kg_group = sorted(
                [h for h in group_all if h["kingdom"] == kingdom],
                key=sort_key,
            )
            if not kg_group:
                continue
            story.append(CondPageBreak(20*mm))
            story.append(P(kingdom, st_kg_lbl))
            story.append(make_hit_table(kg_group))
            story.append(Spacer(1, 1.5*mm))

    # Top hits: fixed number (TOP_N_HITS), prioritised by rarity, reliably
    # fills page 1 instead of leaving a lot of white space when there are
    # few rare hits. Grouped by category (kraken_tree.CATEGORY_ORDER),
    # still descending by z-score within each group.
    if top_hits:
        story.append(HR(5*mm, 2.5*mm))
        story.append(P("Top Hits", st_section, spaceBefore=0))
        kingdom_grouped_tables(top_hits, priority_key)
    else:
        story.append(HR(5*mm, 2.5*mm))
        story.append(P("No genus meets the top-hits criteria.", st_body))

    # Supplement: everything outside the top hits, with the same metrics/
    # boxplots/z-scores as the top hits, split into "additional rare"
    # (< 5% genus prevalence, but no longer among the top N) and "common"
    # findings (>= 5% genus prevalence, background).
    has_supplement = bool(rest_rare or rest_common)

    if has_supplement:
        story.append(PageBreak())

        st_supp_head = S("_sph", fontName="Helvetica-Bold", fontSize=12,
                         textColor=C_MID, spaceBefore=0, spaceAfter=2*mm,
                         leading=15, leftIndent=-6)
        st_supp_note = S("_spn", fontName="Helvetica", fontSize=7,
                         textColor=C_MID, leading=10, spaceAfter=2*mm,
                         leftIndent=-6)

        def supplement_block(title, note, group_all, sort_key):
            if not group_all:
                return
            story.append(P(title, st_supp_head))
            story.append(P(note, st_supp_note))
            kingdom_grouped_tables(group_all, sort_key)

        supplement_block(
            "Supplement – Rare Findings",
            f"Genera not among the top hits, detected in &lt; "
            f"{RARE_THRESHOLD_5PCT:g}% of the reference cases. Sorted by z-score.",
            rest_rare,
            priority_key,
        )
        if rest_rare and rest_common:
            story.append(HR(3*mm, 2*mm))
        supplement_block(
            "Supplement – Common Findings",
            f"Genera not among the top hits, detected in &ge; "
            f"{RARE_THRESHOLD_5PCT:g}% of the reference cases. Sorted by reads.",
            rest_common,
            lambda h: -h["reads_display"],
        )

    # ── LEGEND ────────────────────────────────────────────────────────────
    # All thresholds are taken from the constants above, so the text always
    # matches the logic.
    story.append(HR(5*mm, 3*mm))
    story.append(P("Methodological Notes", st_section, spaceBefore=0))

    n_ref = parsed.get('ref_total', 437)
    t_n = THRESHOLDS["NOTABLE"]
    forced_cap = min(FORCED_MAX_N, TOP_N_HITS - MIN_RARITY_SLOTS)
    legend = [
        ("Reference cohort",
         f"{n_ref} NGS cases from University Hospital Münster."),
        ("Inclusion",
         "Taxa of rank species or below, species groups and taxa without "
         f"rank, with &ge; {t_n['minReads']} directly assigned reads, "
         f"&ge; {t_n['minKmers']} distinct k-mers, dup &le; {t_n['maxDup']} "
         f"(no limit for viruses) and Cov. &ge; {t_n['minCov'] * 100:g}%. "
         "Curated pathogens are also included with &ge; 3 reads and &ge; 10 "
         "k-mers (same dup limit, no Cov. minimum)."),
        ("Genus rows",
         "Included taxa are combined per genus (taxa without a genus form "
         "their own row): reads and k-mers summed, dup as read-weighted mean, "
         "Cov. as maximum. Grey line: the 3 species/"
         "strains with the most k-mers (k-mer count in brackets)."),
        ("Categories",
         "From the KrakenUniq taxonomy: superkingdom Viruses, Bacteria or "
         "Archaea; within Eukaryota, kingdom Fungi; all other eukaryotes "
         "(protozoa, helminths, oomycetes, algae) under \"Parasites and other "
         "eukaryotes\"."),
        ("Top hits",
         f"Up to {TOP_N_HITS} genera, sorted by z-score. Not eligible: genera "
         f"tagged Contaminant and non-viral genera with dup &ge; "
         f"{TOP_HITS_MAX_DUP:g}. Eligible genera tagged Pathogen (&ge; "
         f"{FORCED_MIN_READS['Pathogen']} reads) or Commensal (&ge; "
         f"{FORCED_MIN_READS['Commensal']} reads) are included regardless of "
         f"z-score (max. {forced_cap}; pathogens before commensals, each by "
         "read count); the remaining places go to the highest z-scores. All "
         "other genera are listed in the supplement."),
        ("z-score",
         "log10(adjusted k-mers + 1) of the genus, compared with its mean and "
         f"SD in the {n_ref} reference cases (cases without the genus count as "
         "0). Adjusted k-mers: k-mers / &radic;dup per species/strain (not for "
         "viruses), summed per genus. Genera absent from the reference: mean 0, "
         f"SD pooled across all genera. Red: z &ge; {Z_THRESHOLD_HIGH:g}; "
         f"orange: z &ge; {Z_THRESHOLD_ELEVATED:g}. Small number below: "
         "reference cases (n, %) in which the genus was detected."),
        ("Boxplots",
         f"Distribution of the genus in the {n_ref} reference cases, log scale "
         "(box: 25th–75th percentile, line: median, whiskers: min–max; red "
         "dot: this sample). k-mers and Cov. count cases without the genus as "
         "0, dup uses only cases with the genus. \"–\": genus absent from the "
         "reference. Small \"z\": z-score of k-mers or Cov. alone, without dup "
         "adjustment."),
        ("Columns",
         "Reads: reads assigned directly to the included taxa (KrakenUniq "
         "taxReads). k-mers: distinct k-mers. dup: mean number of times each "
         "distinct k-mer was found. Cov.: proportion of the taxon's k-mers in "
         "the database that were found."),
        ("Comment",
         "From curated lists. Pathogen: established CNS pathogen. Commensal: "
         "human flora with opportunistic CNS pathogenic potential. "
         "Contaminant: typical reagent, skin or environmental contaminant. "
         "Within a genus the highest category applies (Pathogen &gt; "
         "Commensal &gt; Contaminant)."),
    ]
    for term, desc in legend:
        row = Table(
            [[P(f"<b>{term}</b>", st_leg_term), P(desc, st_leg_desc)]],
            colWidths=[30*mm, CW - 30*mm]
        )
        row.setStyle(TableStyle([
            ("VALIGN",       (0,0), (-1,-1), "TOP"),
            ("LEFTPADDING",  (0,0), (-1,-1), 0),
            ("RIGHTPADDING", (1,0), (1,-1),  0),
            ("TOPPADDING",   (0,0), (-1,-1), 2*mm),
            ("BOTTOMPADDING",(0,0), (-1,-1), 2*mm),
            ("LINEBELOW",    (0,0), (-1,-1), 0.3, C_RULE),
        ]))
        story.append(row)

    # ── REFERENCES ───────────────────────────────────────────────────────
    story.append(HR(3*mm, 3*mm))
    st_ref = S("_ref", fontName="Helvetica", fontSize=7, textColor=C_LIGHT,
               leading=10, leftIndent=-6, spaceAfter=2*mm)
    story.append(P(
        "<b>Classification:</b> Breitwieser FP, Baker DN, Salzberg SL. "
        "KrakenUniq: confident and fast metagenomics classification using unique "
        "k-mer counts. <i>Genome Biol.</i> 2018;19(1):198. "
        "doi:10.1186/s13059-018-1568-0",
        st_ref))
    story.append(P(
        "<b>Contaminant list:</b> based on Laurence M, Hatzis C, Brash DE. "
        "Common contaminants in next-generation sequencing that hinder discovery "
        "of low-abundance microbes. <i>PLoS One.</i> 2014;9(5):e97876. "
        "doi:10.1371/journal.pone.0097876, supplemented with further curated "
        "species.",
        st_ref))

    story.append(Spacer(1, 3*mm))
    return story


# ─────────────────────────────────────────────────────────────────────────────
# Two-pass PDF build
# ─────────────────────────────────────────────────────────────────────────────
class _Counter:
    def __init__(self):
        self.total = 1

def generate_pdf(output_path, report_text=None, sample_name=None, date=None):
    """
    Main entry point.
    report_text: raw KrakenUniq report string (if None, uses hardcoded demo data)
    """
    from datetime import date as dt_date

    if report_text is None:
        # Demo data mirroring Barcode 40 analysis
        parsed = _demo_parsed()
        s_name = "PBK06944 · Barcode 40"
        s_date = "10 April 2026"
    else:
        parsed = parse_report(report_text)
        s_name = sample_name or "Unknown sample"
        s_date = date or dt_date.today().strftime("%d %B %Y").lstrip("0")

    counter = _Counter()  # instance variable, safe for batch usage

    def _build(path, draw_fn):
        doc = SimpleDocTemplate(
            path,
            pagesize=A4,
            leftMargin=ML, rightMargin=MR,
            topMargin=MT + 8*mm,   # space for header chrome
            bottomMargin=MB + 8*mm,
            title="Metagenomic Analysis",
            author="Institute of Neuropathology, University Hospital Münster",
        )
        doc.build(
            build_story(s_name, s_date, parsed),
            onFirstPage=draw_fn,
            onLaterPages=draw_fn,
        )
        return doc.page

    # Pass 1: count pages
    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
        tmp_path = tmp.name
    counter.total = _build(tmp_path, lambda c, d: None)
    try:
        os.unlink(tmp_path)
    except Exception:
        pass

    # Pass 2: render with correct page numbers
    def draw(c, doc):
        draw_chrome(c, doc, s_name, counter.total)

    _build(output_path, draw)
    print(f"PDF created: {output_path}  ({counter.total} page(s))")


# ─────────────────────────────────────────────────────────────────────────────
# Demo data (Barcode 40)
# ─────────────────────────────────────────────────────────────────────────────
def _demo_parsed():
    """
    Demo data (species/strain level, as a real report provides them).
    Genus grouping, boxplots and z-scores are computed at display time by
    select_top_hits()/group_by_genus()/attach_genus_stats() -- as in real
    operation -- not hard-coded here. Without a real
    analysis/genus_reference_stats.csv row for the (made-up) demo genus
    taxIDs, the z-score falls back to the "silence" baseline (mean=0,
    pooled SD) -- sufficient for demo purposes.
    """
    demo_ref_total = 437

    def _hit(taxid, name, genus_taxid, genus_name, kingdom, reads, kmers, dup, cov,
             is_contaminant, is_opportunist, is_pathogen, ref_n, ref_pct):
        comment = ("Pathogen" if is_pathogen else
                   "Commensal" if is_opportunist else
                   "Contaminant" if is_contaminant else "")
        return dict(
            taxID=taxid, name=name, genus_taxid=genus_taxid, genus_name=genus_name,
            kingdom=kingdom, rank="species",
            taxReads=reads, reads_display=reads, kmers=kmers, dup=dup, cov=cov,
            is_contaminant=is_contaminant, is_opportunist=is_opportunist,
            is_pathogen=is_pathogen, comment=comment,
            ref_n=ref_n, ref_pct=ref_pct, ref_total=demo_ref_total,
        )

    return {
        "platform":           "Nanopore Metagenomics",
        "total_all":          958615,
        "total_classified":   814549,
        "total_unclassified": 144066,
        "total_human":        680000,
        "pct_classified":     85.0,
        "pct_unclassified":   15.0,
        "pct_human":          70.9,
        "ref_total":          demo_ref_total,
        "hits": [
            _hit("d1", "Human betaherpesvirus 5", "dg1", "Cytomegalovirus", "Viruses",
                 28, 2525, 1.05, 0.01108, False, False, True, ref_n=3, ref_pct=0.25),
            _hit("d2", "Human betaherpesvirus 6B", "dg2", "Roseolovirus", "Viruses",
                 2, 4, 13.0, 3.32e-5, False, False, True, ref_n=5, ref_pct=0.42),
            _hit("d3", "Staphylococcus epidermidis", "dg3", "Staphylococcus", "Bacteria",
                 150, 16593, 1.04, 1.75e-3, True, False, False, ref_n=192, ref_pct=16.1),
            _hit("d4", "Staphylococcus hominis", "dg3", "Staphylococcus", "Bacteria",
                 12, 210, 1.4, 6.0e-5, True, False, False, ref_n=192, ref_pct=16.1),
            _hit("d5", "Burkholderia cenocepacia", "dg5", "Burkholderia", "Bacteria",
                 2, 99, 1.0, 4.82e-6, False, False, False, ref_n=283, ref_pct=23.7),
            _hit("d6", "Escherichia coli", "dg6", "Escherichia", "Bacteria",
                 18, 43, 1.63, 6.08e-7, False, True, False, ref_n=416, ref_pct=34.8),
            _hit("d7", "Trichosporon asahii", "dg7", "Trichosporon", "Fungi",
                 59, 42, 1020.0, 1.78e-6, False, False, True, ref_n=9, ref_pct=0.75),
        ],
    }


if __name__ == "__main__":
    import argparse
    from datetime import date as _date

    parser = argparse.ArgumentParser(
        description="KrakenUniq PDF Report Generator"
    )
    parser.add_argument(
        "--input", "-i",
        help="KrakenUniq report .txt file. Omit to use built-in demo data.",
        default=None,
    )
    parser.add_argument(
        "--output", "-o",
        help="Output PDF path (default: <input_basename>.metagenomics_report.pdf)",
        default=None,
    )
    parser.add_argument(
        "--sample", "-s",
        help="Sample name shown in the report header (default: input filename)",
        default=None,
    )
    parser.add_argument(
        "--date", "-d",
        help="Date string for the report (default: today)",
        default=None,
    )
    args = parser.parse_args()

    # Resolve output path
    if args.output is None:
        if args.input:
            base = os.path.splitext(os.path.splitext(args.input)[0])[0]
            args.output = base + ".metagenomics_report.pdf"
        else:
            args.output = "KrakenUniq_Report_demo.pdf"

    # Sample name fallback
    sample_name = args.sample
    if sample_name is None and args.input:
        sample_name = os.path.basename(args.input).replace(".krakenuniq.report.txt", "")

    # Date fallback
    report_date = args.date or _date.today().strftime("%d %B %Y").lstrip("0")

    # Read report text
    report_text = None
    if args.input:
        with open(args.input, "r", encoding="utf-8") as fh:
            report_text = fh.read()

    generate_pdf(
        output_path=args.output,
        report_text=report_text,
        sample_name=sample_name,
        date=report_date,
    )
