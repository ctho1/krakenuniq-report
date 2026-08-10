"""
KrakenUniq PDF Report Generator – v5
Änderungen gegenüber v4:
  - "Top-Hits"-Tabelle: feste Anzahl der priorisierten Treffer (seltenste
    zuerst, bei Bedarf mit häufigeren aufgefüllt) füllt Seite 1 zuverlässig,
    statt bei wenigen seltenen Treffern viel Weißraum zu lassen
  - Vertikale Mini-Boxplots (Referenzverteilung, >=10kb-Kohorte) für k-mers
    und Coverage je Treffer, mit rotem Punkt für den aktuellen Wert
  - Restliche Treffer (seltene über Top-N hinaus + häufige) im kompakten
    Supplement, getrennt nach "weitere seltene" und "häufige" Befunde

Änderungen gegenüber v3:
  - Schlankerer Kopfbereich (Titel/Meta/Stat-Kacheln kompakter)
  - Automatisch generierter Befund-Fließtext entfernt
  - Referenz-Prävalenz (n / % der Referenzdatenbank, Genus-Ebene) je Treffer
    statt Plausibilitäts-Badge ("Detektiert" / "Auffällig" entfernt)
  - Kommentarspalte: Pathogen / Kommensale / Kontaminant
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

# Column widths (mm) – Top-Hits-Tabelle (Seite 1, inkl. Mini-Boxplots)
# Rang wird als graue Unterzeile im Organismus-Feld geführt, keine eigene
# Spalte mehr. k-mers/dup/Cov. zeigen Zahl + Boxplot gestapelt in EINER
# Spalte (spart Breite gegenüber getrennten Wert-/Chart-Spalten).
# Total must equal CW (170 mm)
COL_W = {
    "organismus": 35*mm,
    "reads":      14*mm,
    "kmers":      24*mm,
    "dup":        20*mm,
    "cov":        26*mm,
    "referenz":   22*mm,
    "kommentar":  29*mm,
}
assert abs(sum(COL_W.values()) - CW) < 1, \
    f"Col widths sum {sum(COL_W.values())/mm:.1f}mm ≠ CW {CW/mm:.1f}mm"
# Horizontale Boxplots: Breite = volle Zellbreite (KMERS_INNER_W/COV_INNER_W/
# DUP_INNER_W), Höhe klein, um die Zeilenhöhe der Top-Hits-Tabelle gering zu
# halten.
BOX_H = 4*mm
KMERS_INNER_W = COL_W["kmers"] - 4*mm
DUP_INNER_W   = COL_W["dup"]   - 4*mm
COV_INNER_W   = COL_W["cov"] - 4*mm


# Prävalenz-Schwellenwerte für die Tier-Struktur (bezogen auf die komplette
# Referenzdatenbank, nicht nur die >=10kB-Kohorte)
RARE_THRESHOLD_1PCT = 1.0
RARE_THRESHOLD_5PCT = 5.0

# z-Score-Schwellenwerte (Standardabweichungen vom log10(Reads+1)-Hinter-
# grundrauschen des Taxons in der >=10kB-Kohorte) für die Farbmarkierung in
# der Top-Hits-Tabelle.
Z_THRESHOLD_HIGH = 3.0
Z_THRESHOLD_ELEVATED = 1.5

# Anzahl der priorisierten Treffer in der Top-Hits-Tabelle auf Seite 1.
# Wird bei wenigen seltenen Treffern mit den nächst-häufigeren aufgefüllt,
# damit Seite 1 zuverlässig gefüllt ist statt viel Weißraum zu zeigen.
TOP_N_HITS = 15

# Mindest-Reads, ab denen ein kuratierter Pathogen/Kommensalen-Treffer
# unabhängig von seinem z-Score zwingend in die Top-Hits aufgenommen wird
# (siehe select_top_hits()). Für "Pathogen" (etablierte ZNS-Erreger, kleine
# kuratierte Liste) niedriger angesetzt als für "Kommensale" (humane
# Standortflora, höhere Kontaminationswahrscheinlichkeit → mehr Evidenz
# nötig). Seit dem Wechsel auf reine z-Score-Sortierung (ohne Genus-
# Prävalenz-Tier als Vorsortierung) senkt dies gezielt Fälle mit sehr
# niedrigem, aber echtem Signal ab: Referenzkohorten-Validierung zeigte
# HIV-1 (9 Reads), BoDV-1 (5/8 Reads) und SARS-CoV-2 (5 Reads) sonst
# außerhalb der Top-15, weil in ihren (z.T. stark polymikrobiellen) Proben
# andere Treffer einen höheren z-Score erreichten.
FORCED_MIN_READS = {"Pathogen": 5, "Kommensale": 20}
# Deckel für die Anzahl forcierter Treffer, damit polymikrobiell besiedelte
# Proben nicht die komplette Top-Hits-Tabelle mit Kommensalen füllen.
FORCED_MAX_N = 10
# Mindestanzahl Plätze, die in der Top-Hits-Tabelle immer der reinen
# Genus-Seltenheits-Priorisierung vorbehalten bleiben (nicht von forcierten
# Treffern verdrängbar).
MIN_RARITY_SLOTS = 5

# Zusätzliches Aufnahmekriterium NUR für die Top-Hits-Tabelle (Seite 1),
# auf Wunsch: eine Genus-Zeile (read-gewichteter mittlerer dup-Wert, siehe
# group_by_genus()) mit dup >= TOP_HITS_MAX_DUP wird nicht mehr in die
# Top-Hits aufgenommen (weder forciert noch per z-Score), erscheint aber
# weiterhin im Supplement, sofern sie den Evidenzfilter besteht -- der
# Evidenzfilter selbst (assess()) bleibt unverändert, ein Treffer verschwindet
# dadurch also nicht komplett aus dem Report. Kingdom "Viren" ist davon
# ausgenommen, aus demselben Grund wie beim dup-Cap im Evidenzfilter
# (VIRAL_MAX_DUP): kleine virale Genome erzeugen bei echter, tiefer Coverage
# strukturell hohe dup-Werte (z.B. BoDV-1 dup 50-400, TBEV dup 2920, JC-Virus
# dup 215 in der Referenzkohorte) -- ein ungefiltertes dup<2 hätte praktisch
# alle bestätigten viralen Top-Hits verdrängt.
TOP_HITS_MAX_DUP = 2.0

# ── Mindestkriterien für Aufnahme in den Report (Rauschunterdrückung) ──────
THRESHOLDS = {
    "CONVINCING": dict(minReads=10, minKmers=500, maxDup=5,  minCov=1e-4),
    "NOTABLE":    dict(minReads=5,  minKmers=50,  maxDup=30, minCov=1e-6),
}

# dup ("mittlere k-mer-Duplikationsrate") skaliert mit Sequenziertiefe pro
# Genomgröße. Kleine virale Genome erzeugen bei genuine sehr hoher Coverage
# (z.B. TBEV mit 750.000+ Reads auf einem ~11kb-Genom) strukturell sehr hohe
# dup-Werte (>1000), ohne dass das Kontamination/Kreuzreaktivität bedeutet.
# Referenzkohorten-Validierung (63 Fälle mit bestätigtem Erreger, Supplementary
# Table 2): 6 von 6 viralen Fällen mit dup>30 waren echte, stark abgedeckte
# Erreger (TBEV dup=2920, BoDV-1 dup=50-400, JCV dup=215) – der dup<=30-Cutoff
# hätte sie trotz zehntausender Reads verworfen. Für die Kingdom "Viren" wird
# der dup-Cutoff daher aufgehoben; bei Bakterien/Pilzen/Parasiten (deutlich
# größere Genome) bleibt dup ein sinnvoller Kontaminations-/Kreuzreaktivitäts-
# Indikator und wird weiter angewendet.
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
# OPPORTUNISTIC_FLORA  – human commensals with opportunistic ZNS pathogen potential
#                        (e.g. oral streptococci in brain abscess); labelled separately
# KNOWN_PATHOGENS      – established ZNS pathogens (defined below)

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
    # ── Food / fermentation organisms (kein humanpathogenes ZNS-Potenzial) ────
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
    # since group_by_genus() already lets a real Pathogen/Kommensale member
    # species override a Kontaminant one within the same genus row.
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
    # ── Gut / oral commensal flora (no ZNS pathogen potential) ─────────────────
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

# Opportunistic flora: human commensals that can cause ZNS infections
# (brain abscess, meningitis in immunocompromised) but are also frequent
# biopsy/FFPE contaminants. Labelled "Auffällig – Opportunist" to distinguish
# from pure environmental contaminants and primary pathogens.
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
    "Haemophilus parainfluenzae",      # rare ZNS infections
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
    # ZNS-relevante Primärpathogene (kuratierte Liste für Signalverstärkung)
    # Referenz: CZ ID pathogen list (https://czid.org/pathogen_list, Chan Zuckerberg Initiative)
    # ── Viren ─────────────────────────────────────────────────────────────────
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
    "Alphapolyomavirus quintihominis",  # JC-Virus
    "JC polyomavirus",
    "Betapolyomavirus hominis",         # BK-Virus
    "BK polyomavirus",
    "Lyssavirus rabies",            # Tollwut
    "Morbillivirus hominis",        # Masern
    "Orthobornavirus bornaense",    # Bornavirus-Enzephalitis
    "Borna disease virus 1",        # NCBI-Taxonomiename für BoDV-1 in KrakenUniq
    "Human immunodeficiency virus 1",  # HIV-1 (ZNS-Beteiligung/HIV-Enzephalopathie)
    "Orthoflavivirus nilense",      # West-Nil-Virus
    "West Nile virus",
    "Orthoflavivirus japonicum",    # Japanische Enzephalitis
    "Orthoflavivirus encephalitidis", # FSME
    "Tick-borne encephalitis virus",
    "Severe acute respiratory syndrome-related coronavirus",
    "Severe acute respiratory syndrome coronavirus 2",   # NCBI taxon name in KrakenUniq
    "SARS-CoV-2",
    # ── Bakterien ────────────────────────────────────────────────────────────
    "Listeria monocytogenes",
    "Neisseria meningitidis",
    "Streptococcus pneumoniae",
    "Mycobacterium tuberculosis",
    "Mycobacterium tuberculosis complex",  # species-group-Knoten in KrakenUniq
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
    # ── Pilze ────────────────────────────────────────────────────────────────
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
    # ── Parasiten ────────────────────────────────────────────────────────────
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


infer_kingdom = kraken_tree.infer_kingdom


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
    kingdom = infer_kingdom(name)

    is_contaminant    = name in LIKELY_CONTAMINANTS
    is_opportunist    = name in OPPORTUNISTIC_FLORA
    is_known_pathogen = name in KNOWN_PATHOGENS

    t_c = THRESHOLDS["CONVINCING"]
    t_n = THRESHOLDS["NOTABLE"]
    # dup-Cutoff für Viren aufgehoben (siehe VIRAL_MAX_DUP-Kommentar oben).
    max_dup_c = VIRAL_MAX_DUP if kingdom == "Viren" else t_c["maxDup"]
    max_dup_n = VIRAL_MAX_DUP if kingdom == "Viren" else t_n["maxDup"]

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
    # diversity, AND an acceptable dup value. reads>=3 (statt >=5) erfasst
    # Grenzfälle mit sehr niedriger, aber diagnostisch bestätigter Erreger-
    # last (Referenzkohorten-Validierung: A. fumigatus mit 4 Reads/30 k-mers
    # bei ITS-bestätigter Aspergillose). Der dup-Guard bleibt, um False
    # Positives wie Trypanosoma (dup=234), Cryptosporidium (dup=864) zu
    # vermeiden – für Viren gilt dabei der aufgehobene max_dup_n.
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
        comment = "Kommensale"
    elif is_contaminant:
        comment = "Kontaminant"
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
# Referenzdatenbank (Prävalenz über alle bisherigen KrakenUniq-Reports)
# ─────────────────────────────────────────────────────────────────────────────
_REF_DB_CACHE = None


def load_reference_db():
    """
    Lädt die Genus-Ebene-Prävalenz-Referenzdatenbank
    (analysis/taxon_prevalence_genus_min10kb.csv, erzeugt von
    scripts/build_reference_db.py). Bezieht sich bewusst auf die >=10kB-
    Kohorte (n=437) statt der kompletten Datenbank (n=1194): dieselbe
    Kohorte, aus der auch die Boxplot-/z-Score-Referenzverteilungen
    (taxon_reference_stats.csv) stammen – sehr kleine/nahezu leere Reports
    (Median-Reportgröße < 10kB) verzerren sonst die Prävalenzschätzung, ohne
    zur Auflösung beizutragen.

    Genus-Ebene statt Species/Strain-Ebene: eine einzelne Gattung (z.B.
    Bradyrhizobium) kann in der Referenzdatenbank über hunderte einzelne
    Stamm-/Isolat-Namen verteilt sein, von denen jeder für sich genommen
    selten ist – obwohl die Gattung insgesamt ein häufiger Hintergrund-
    kontaminant ist. Die Priorisierung im Report erfolgt daher auf
    Genus-Ebene; angezeigt wird weiterhin der genaue Species/Strain-Name.

    Ergebnis wird pro Prozess gecached, da die Tabelle mehrere tausend
    Zeilen hat.
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
    """Reichert einen Treffer um Genus-Ebene-Referenz-Prävalenz (n / % von n_total) an."""
    n_present, pct_present = ref_db.get(hit["genus_taxid"], (0, 0.0))
    hit["ref_n"] = n_present
    hit["ref_pct"] = pct_present
    hit["ref_total"] = n_total
    return hit


# Gemeinsames Feld-Set für taxon_reference_stats.csv UND genus_reference_stats.csv
# (identisches Schema, siehe build_reference_db.py: ref_stats_row()).
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
    Lädt die gepoolte Standardabweichung von log10(Reads+1) über ALLE
    Treffer der >=10kB-Referenzkohorte (analysis/reference_meta.json).
    Dient als Fallback-Nenner für den z-Score, wenn ein Taxon zu selten in
    der Referenz vorkommt (n<3), um eine eigene, stabile SD zu schätzen.
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
    z-Score: wie viele Standardabweichungen liegt dieser Treffer (log10-Skala)
    über dem typischen Hintergrundrauschen für dieses Taxon in der
    >=10kB-Referenzkohorte?

    Standardmäßig (transform="add1", value = coverage-adjustierte k-mer-Zahl,
    siehe kraken_tree.adjusted_kmers()): k-mers trennen echte Erreger von
    Hintergrundrauschen besser als rohe Reads (siehe
    analysis/threshold_analysis_report.md, AUC 0.873 vs. 0.755), und die
    Adjustierung durch dup gewichtet breite, glaubwürdige Genom-Coverage
    (viele verschiedene Positionen je einmal getroffen) höher als eine
    schmale, stark duplizierte Signatur (wenige Positionen, oft getroffen)
    -- außer bei Viren, wo hohe dup-Werte durch echte tiefe Coverage kleiner
    Genome entstehen und daher nicht abgewertet werden. log10(value+1), da
    value als Zähldaten (k-mers, Reads) auch 0 sein kann.

    transform="direct": für Größen, die immer winzige Brüche sind (Coverage-
    Anteil) -- log10(value+1) würde für jeden Wert auf ~0 kollabieren und das
    Signal zerstören, daher wird das unverschobene log10(max(value,1e-9))
    verwendet (identisch zur Referenzverteilung, siehe build_reference_db.py
    log_mean_sd_direct()).

    Referenz-Mittelwert/SD (box_stats[mean_key]/[sd_key], siehe
    build_reference_db.py: ref_stats_row()) werden über die GESAMTE
    >=10kB-Kohorte berechnet, nicht nur über die Fälle mit Nachweis -- Proben
    ohne Nachweis fließen explizit mit 0 ein. Dadurch ist die taxon-eigene
    SD bereits ab einem einzigen Referenz-Vorkommen gut konditioniert (bei
    seltenen Taxa von den vielen "Abwesenheits"-Nullen dominiert), weshalb
    hier standardmäßig ab n_ref>=1 die eigene SD verwendet wird:

      - Taxon mit >=1 Referenz-Vorkommen UND SD>0: eigener Mittelwert/eigene SD.
      - Taxon noch nie in der Referenz gesehen (n=0): Hintergrund wird als
        "Stille" angenommen (mean=0), gepoolte SD (über alle Taxa) als Nenner.
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
    Lädt die Genus-Ebene-Verteilungskennzahlen (analysis/genus_reference_stats.csv,
    >=10kB-Kohorte): pro Sample über alle Species/Strains einer Gattung
    summierte Reads/k-mers (Coverage als Max, dup Read-gewichtet), dann über
    die Kohorte aggregiert. Grundlage für die Genus-Ebene-Boxplots/z-Scores
    im Report (siehe group_by_genus()).
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
    Gepoolte SD (log-Skala) einer Genus-Ebene-Metrik über die >=10kB-
    Referenzkohorte (analysis/reference_meta.json) -- Fallback-Nenner für
    compute_z_score(), wenn eine Gattung zu selten in der Referenz vorkommt,
    um eine eigene, stabile SD zu schätzen.
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
    Aggregiert Species/Strain-Ebene-Treffer (bereits durch assess() gefiltert)
    zu einer Zeile je Gattung:
      - reads_display / kmers: Summe über alle Species/Strains der Gattung
        in diesem Report
      - cov: Maximum (Coverage ist genomgrößenabhängig und daher nicht
        sinnvoll summierbar; der Max-Wert zeigt die stärkste Einzelevidenz)
      - dup: Read-gewichteter Mittelwert
      - adj_kmers: Summe der coverage-adjustierten k-mer-Zahlen der
        Mitglieder (kraken_tree.adjusted_kmers(), je Mitglied mit dessen
        eigenem dup-Wert berechnet, dann summiert) -- Basis des z-Scores,
        siehe attach_genus_stats()/compute_z_score()
      - comment: höchste Priorität unter den Mitgliedern
        (Pathogen > Kommensale > Kontaminant > "")
      - top3: die 3 Species/Strains mit den meisten k-mers (für die
        Organismus-Unterzeile im Report)
    z-Score/Boxplot-Referenzdaten (box_stats, z_score) werden hier NICHT
    gesetzt – dafür attach_genus_stats() nach dem Gruppieren aufrufen.
    """
    groups = defaultdict(list)
    for h in hits:
        groups[h["genus_taxid"]].append(h)

    comment_priority = {"Pathogen": 0, "Kommensale": 1, "Kontaminant": 2, "": 3}

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
    Setzt z_score (Haupt-Score, coverage-adjustierte k-mers, siehe
    compute_z_score()) sowie die separaten kmers_z_score/cov_z_score --
    dieselbe Methodik, aber auf den rohen k-mers bzw. der Coverage dieser
    Gattung, unter dem jeweiligen Boxplot in der Top-Hits-Tabelle angezeigt
    (make_hit_table()), damit sichtbar wird, wie stark k-mer-Zahl und
    Coverage EINZELN vom Hintergrund abweichen, nicht nur die kombinierte
    dup-adjustierte Metrik.
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
    # z-Score/Boxplot-Referenzdaten werden auf Genus-Ebene benötigt (siehe
    # group_by_genus()/attach_genus_stats() in select_top_hits()), nicht
    # mehr je einzelnem Species/Strain-Treffer.

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
    return f"{int(n):,}".replace(",", ".")

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
    """1 Nachkommastelle, ab >= 100 gerundet auf Ganzzahl (Spaltenbreite)."""
    if dup >= 100:
        return f"{dup:.0f}"
    return f"{dup:.1f}"


def comment_text_color(comment):
    """Textfarbe der Kommentarspalte je Kategorie."""
    return {
        "Pathogen":    C_RED_TEXT,
        "Kommensale":  C_PURPLE_TEXT,
        "Kontaminant": C_GREY_TEXT,
    }.get(comment, C_MID)




def zscore_tier_props(z):
    """
    Hintergrundfarbe der z-Score-Zelle in der Top-Hits-Tabelle.
    >= Z_THRESHOLD_HIGH     -> starke Abweichung vom Hintergrundrauschen (rötlich)
    >= Z_THRESHOLD_ELEVATED -> erhöhte Abweichung (orange)
    sonst                   -> im Rahmen des üblichen Rauschens (neutral)
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
    """Kompakte Zahl für die Top-3-Species-Unterzeile, z.B. 8212 -> '8.2k'."""
    if n >= 1000:
        s = f"{n/1000:.1f}k"
        return s.replace(".0k", "k")
    return str(n)


def fmt_top3_species(genus_row):
    """
    Formatiert die 3 k-mer-stärksten Species/Strains einer Genus-Zeile als
    kompakte, mit '·' getrennte Liste für die Organismus-Unterzeile im
    Report, z.B. "coli (8.2k) · albertii (540) · fergusonii (12)". Der
    Gattungsname wird aus dem Species-Namen entfernt, falls er als Präfix
    vorkommt (spart Platz, da die Gattung bereits in der Hauptzeile steht).
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


# ── Mini-Boxplots (Referenzverteilung k-mers/Coverage, log-Skala) ─────────
C_BOX_FILL = HexColor("#E5E5EA")   # blasses Grau (Referenzverteilung)
C_BOX_LINE = HexColor("#AEAEB2")
C_BOX_DOT  = HexColor("#FF3B30")   # roter Punkt: aktueller Fall


def _log_frac(v, lo, hi):
    """Position von v auf einer Log10-Skala zwischen lo und hi, als 0..1."""
    eps = 1e-12
    v, lo, hi = max(v, eps), max(lo, eps), max(hi, eps)
    if hi <= lo:
        return 0.5
    frac = (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
    return min(max(frac, 0.0), 1.0)


def make_boxplot(current_value, stats, prefix, width, height):
    """
    Horizontaler Mini-Boxplot auf Log10-Skala (spart Zeilenhöhe gegenüber
    der vertikalen Variante):
      - graue Box: Referenz-Interquartilsbereich (p25-p75), >=10kB-Kohorte
      - graue Linie: Referenz-Median (p50)
      - graue Whisker: Referenz-Min/Max
      - roter Punkt: Wert des aktuellen Falls
    Gibt None zurück, wenn keine Referenzdaten für diesen Taxon vorliegen
    (z.B. neu/nur außerhalb der >=10kB-Kohorte beobachtet).
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
    c.drawString(ML, h - MT + 8*mm, "Metagenomics Untersuchung")

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
                 "Institut für Neuropathologie · Universitätsklinikum Münster")
    c.drawRightString(w - MR, MB - 7*mm,
                      f"Seite {page_num} von {total_pages}")


# ─────────────────────────────────────────────────────────────────────────────
# Story builder
# ─────────────────────────────────────────────────────────────────────────────
def priority_key(h):
    """
    Sortierung absteigend nach z-Score (Standardabweichungen vom taxon-
    spezifischen Hintergrundrauschen, >=10kB-Referenzkohorte).

    Hinweis: eine frühere Version sortierte zusätzlich nach Genus-Prävalenz-
    Tier (< 1% / < 5% / >= 5%) VOR dem z-Score, weil reiner z-Score in stark
    polymikrobiell besiedelten (z.B. postmortalen) Proben dazu neigt, auch
    harmlose Kommensalen weit oben einzusortieren (praktisch jede Spezies ist
    dort gegenüber ihrem eigenen, niedrigen Hintergrund extrem erhöht). Die
    Tier-Vorsortierung wurde auf Wunsch entfernt; als Absicherung gegen genau
    dieses Szenario bleibt ausschließlich die forcierte Aufnahme kuratierter
    Pathogene/Kommensalen in select_top_hits() bestehen (mit Vorrang für
    "Pathogen" vor "Kommensale" bei begrenztem Budget). Bei Bedarf erneut mit
    scripts/build_reference_cohort_summary.py gegen die 63-Fälle-Kohorte
    validieren.
    """
    return -h["z_score"]


def select_top_hits(hits):
    """
    Gruppiert die Species/Strain-Ebene-Treffer eines Reports zu Genus-Zeilen
    (group_by_genus() + attach_genus_stats(), siehe dort) und wählt daraus
    die Top-Hits-Tabelle (Seite 1). Einzige Implementierung dieser Auswahl-
    logik – wird sowohl vom PDF-Report (build_story) als auch von der
    Referenzkohorten-Validierung (scripts/build_reference_cohort_summary.py)
    verwendet, damit beide nie auseinanderlaufen.

    Priorisierte Gesamtliste: absteigend nach z-Score (siehe priority_key()),
    jetzt auf Genus-Ebene berechnet (Reads/k-mers über alle Species/Strains
    der Gattung summiert). Die Top-N füllen Seite 1 zuverlässig auch dann,
    wenn es nur wenige auffällige Treffer gibt (aufgefüllt mit den nächst-
    höchsten z-Scores) – der Rest erscheint mit denselben Metriken/Boxplots
    im Supplement, getrennt nach "weitere seltene" und "häufige" Befunde
    (dort weiterhin anhand der Genus-Prävalenz-Tier, < 5% vs. >= 5%).

    Forcierte Aufnahme unabhängig vom z-Score: Genus-Zeilen mit kuratiertem
    Pathogen/Kommensalen-Anteil (Kommentar-Spalte) und nicht-trivialem
    Read-Support werden nie allein wegen niedrigerem z-Score ins Supplement
    verdrängt. Referenzkohorten-Validierung (63 Fälle mit bestätigtem
    Erreger): ohne diese Regel landeten z.B. S. aureus (43 Reads),
    K. pneumoniae (318), A. fumigatus (133.475) trotz dominanten Signals
    außerhalb der Top-15.

    Die forcierte Liste wird auf FORCED_MAX_N gedeckelt (stärkste Reads
    zuerst): stark polymikrobiell besiedelte Proben (z.B. Autopsie-/
    postmortale Hirnproben mit ausgedehnter bakterieller Überwucherung)
    können 20+ Kommensal-Gattungen mit >= FORCED_MIN_READS enthalten. Ohne
    Deckel würden diese die komplette Top-Hits-Tabelle füllen und den
    eigentlich seltenen, diagnostisch entscheidenden Befund verdrängen.
    Mindestens MIN_RARITY_SLOTS Plätze bleiben daher immer für die reine
    z-Score-Priorisierung reserviert.

    Returns: (top_hits, rest_rare, rest_common) – Listen von Genus-Zeilen.
    """
    genus_stats = load_genus_stats()
    pooled_sd_genus = load_pooled_adj_kmers_log_sd_genus()
    rows = [attach_genus_stats(r, genus_stats, pooled_sd_genus)
            for r in group_by_genus(hits)]

    # Aufnahmekriterien NUR für Top-Hits (dup-Cutoff, siehe TOP_HITS_MAX_DUP,
    # Kingdom "Viren" ausgenommen; kuratierte Kontaminanten, siehe
    # LIKELY_CONTAMINANTS). Treffer, die daran scheitern, bleiben im Pool
    # für rest_rare/rest_common (Supplement) erhalten, sind aber weder
    # forcierbar noch per z-Score in top_hits wählbar. Die Kontaminanten-
    # Ausnahme ist bewusst hart (nicht nur ein Bonus/Malus im z-Score): seit
    # die Referenzverteilung Abwesenheit als 0 mitzählt (s. Änderungshistorie
    # in reference_cohort_filter_strategy.md), ist "selten in der Referenz"
    # der Normalfall für praktisch jede Gattung, nicht spezifisch für echte
    # Erreger – kuratierte Kontaminanten würden sonst genauso oft einen hohen
    # z-Score erreichen wie seltene echte Funde und Seite 1 verwässern. Ein
    # Kontaminant-Tag auf Genus-Ebene entsteht nur, wenn KEIN Mitglied dieser
    # Gattung als Pathogen/Kommensale eingestuft ist (group_by_genus() nimmt
    # die höchste Priorität unter den Mitgliedern), ein "verunreinigter"
    # Fund mit einer echten Erreger-Spezies in derselben Gattung wird also
    # nicht mit ausgeschlossen.
    eligible = [r for r in rows
                if (r["kingdom"] == "Viren" or r["dup"] < TOP_HITS_MAX_DUP)
                and r["comment"] != "Kontaminant"]

    # Kuratierte Pathogene erhalten innerhalb des forcierten Budgets Vorrang
    # vor Kommensalen: in stark polymikrobiellen Proben können Kommensal-
    # Gattungen (z.B. postmortale E. coli-Überwucherung, 13.529 Reads) den
    # kompletten Cap allein durch rohe Read-Zahl belegen und eine echte, aber
    # read-schwächere kuratierte Gattung (z.B. BoDV-1, 92 Reads) verdrängen.
    forced_pathogens = sorted(
        [r for r in eligible if r["comment"] == "Pathogen"
         and r["reads_display"] >= FORCED_MIN_READS["Pathogen"]],
        key=lambda r: -r["reads_display"],
    )
    forced_commensals = sorted(
        [r for r in eligible if r["comment"] == "Kommensale"
         and r["reads_display"] >= FORCED_MIN_READS["Kommensale"]],
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

    # ── TITLE (schmal) ───────────────────────────────────────────────────
    story.append(P("Metagenomics Untersuchung", st_title))
    story.append(P(f"{sample_name} · {date}", st_subtitle))
    story.append(HR(0.5*mm, 2*mm))

    # ── META (single compact row) ─────────────────────────────────────────
    # Truncate very long sample names in the meta row to avoid line breaks
    def trunc(s, max_len=34):
        return s if len(s) <= max_len else s[:max_len - 1] + "…"

    meta_cols = [
        ("PROBE",    trunc(sample_name), 70*mm),
        ("DATUM",    date,               26*mm),
        ("METHODE",  platform,           42*mm),
        ("REFERENZ", f"n = {parsed.get('ref_total', 437)} Fälle", 32*mm),
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

    # ── STAT TILES (schmaler) ───────────────────────────────────────────
    tile_data = [
        (C_TILE_BLUE,
         "READS GESAMT",
         fmt_num(parsed["total_all"]),
         ""),
        (C_TILE_GREEN,
         "KLASSIFIZIERT",
         f"{parsed['pct_classified']:.1f}%",
         fmt_num(parsed["total_classified"]) + " Reads"),
        (HexColor("#30A46C"),
         "HUMAN",
         f"{parsed['pct_human']:.1f}%",
         fmt_num(parsed["total_human"]) + " Reads"),
        (C_TILE_GREY,
         "UNKLASSIFIZIERT",
         f"{parsed['pct_unclassified']:.1f}%",
         fmt_num(parsed["total_unclassified"]) + " Reads"),
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

    # ── BEFUNDTABELLEN ────────────────────────────────────────────────────
    kingdoms = ["Viren", "Bakterien", "Pilze", "Parasiten"]
    cw = list(COL_W.values())

    top_hits, rest_rare, rest_common = select_top_hits(hits)

    # Spaltenindizes in COL_W (Reihenfolge!) für SPAN/ALIGN-Regeln
    IDX = {name: i for i, name in enumerate(COL_W.keys())}

    def stacked_metric_cell(value_para, chart, inner_w, z_para=None):
        """Zahl über Mini-Boxplot (über optionaler z-Score-Zeile) gestapelt,
        als eine kompakte Zelle."""
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
        """Top-Hits-Tabelle: volle Metriken + Mini-Boxplots + Referenz + Kommentar."""
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
            P("Organismus",  s_hdr),
            P("Reads",       s_hdr_r),
            P("k-mers",      s_hdr_c),
            P("dup",         s_hdr_c),
            P("Cov. (%)",    s_hdr_c),
            P("z-Score",     s_hdr_c),
            P("Kommentar",   s_hdr),
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
            # dup hat keinen eigenen z-Score (siehe Legende), bekommt aber
            # eine leere dritte Zeile in derselben Zeilenhöhe wie kmers_z/
            # cov_z (statt None) -- sonst würde die 2-zeilige dup-Zelle
            # innerhalb der (durch kmers/Cov. auf 3 Zeilen gesetzten)
            # Zeilenhöhe vertikal zentriert und der Boxplot läge dadurch
            # sichtbar tiefer als die anderen beiden Boxplots.
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
            ("ALIGN",         (IDX["referenz"],0),(IDX["referenz"],-1),"CENTER"),
        ]
        for ri, hit in enumerate(group, start=1):
            z_fill, _ = zscore_tier_props(hit["z_score"])
            if z_fill is not None:
                style_cmds.append(("BACKGROUND", (IDX["referenz"], ri), (IDX["referenz"], ri), z_fill))
        tbl.setStyle(TableStyle(style_cmds))
        return tbl

    st_kg_lbl = S("_kgl", fontName="Helvetica-Bold", fontSize=8,
                  textColor=C_MID, spaceBefore=2*mm, spaceAfter=1.5*mm,
                  leading=10, leftIndent=-6)

    def kingdom_grouped_tables(group_all, sort_key):
        """Rendert group_all als eine make_hit_table() je Kingdom (Reihenfolge
        s. `kingdoms`), mit Kingdom-Zwischenüberschrift -- gemeinsam von
        Top-Hits und Supplement verwendet, damit beide identisch gegliedert
        und identisch bebildert (Boxplots/z-Scores) sind."""
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

    # Top-Hits: feste Anzahl (TOP_N_HITS), priorisiert nach Seltenheit, füllt
    # Seite 1 zuverlässig statt bei wenigen seltenen Treffern viel Weißraum
    # zu lassen. Gegliedert nach Kingdom (Viren/Bakterien/Pilze/Parasiten),
    # innerhalb jeder Gruppe weiterhin absteigend nach z-Score.
    if top_hits:
        story.append(HR(5*mm, 2.5*mm))
        story.append(P("Top-Hits", st_section, spaceBefore=0))
        kingdom_grouped_tables(top_hits, priority_key)
    else:
        story.append(HR(5*mm, 2.5*mm))
        story.append(P(
            "Keine Befunde oberhalb der Nachweisgrenze "
            "(&ge; 5 Reads &middot; &ge; 50 distinkte k-mers).",
            st_body))

    # Supplement: alles außerhalb der Top-Hits, mit denselben Metriken/
    # Boxplots/z-Scores wie Top-Hits, getrennt nach "weitere seltene"
    # (< 5% Genus-Prävalenz, aber nicht mehr unter den Top-N) und "häufige"
    # Befunde (>= 5% Genus-Prävalenz, Hintergrund).
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
            "Supplement – Weitere seltene Befunde",
            "Treffer mit &lt; 5% Genus-Prävalenz, die aufgrund des Platzlimits "
            f"nicht mehr unter den Top {TOP_N_HITS} auf Seite 1 erscheinen. "
            "Sortiert nach Prävalenz (aufsteigend), dann Reads (absteigend).",
            rest_rare,
            priority_key,
        )
        if rest_rare and rest_common:
            story.append(HR(3*mm, 2*mm))
        supplement_block(
            "Supplement – Häufige Befunde (Hintergrund)",
            f"Treffer mit &ge; 5% Genus-Prävalenz in der Referenzdatenbank (n = "
            f"{parsed.get('ref_total', 437)} Fälle) – mit hoher Wahrscheinlichkeit "
            "Hintergrund oder Kontamination. Sortiert nach Reads (absteigend).",
            rest_common,
            lambda h: -h["reads_display"],
        )

    # ── LEGENDE ───────────────────────────────────────────────────────────
    story.append(HR(5*mm, 3*mm))
    story.append(P("Methodische Hinweise", st_section, spaceBefore=0))

    legend = [
        ("Aufnahmekriterium",
         "≥ 5 Reads · ≥ 50 distinkte k-mers · dup ≤ 30 · Genomabdeckung ≥ 0,0001%. "
         "Bekannte Pathogene werden bei der Aufnahme bevorzugt bewertet. Treffer "
         "unterhalb dieser Kriterien erscheinen nicht im Report."),
        ("Organismus / Top-3-Species",
         "Jede Zeile fasst ALLE Treffer einer Gattung (Genus) im Report "
         "zusammen (Reads/k-mers summiert, s.u.). Die graue Unterzeile zeigt "
         "die 3 Species/Strains dieser Gattung mit den meisten k-mers samt "
         "k-mer-Zahl, z.B. \"coli (8.2k) · albertii (540)\" – der "
         "Gattungsname wird dabei weggelassen, da er bereits in der "
         "Hauptzeile steht."),
        ("Top-Hits",
         f"Die {TOP_N_HITS} Genus-Zeilen mit der höchsten Priorität erscheinen "
         "mit voller Metrik auf Seite 1, absteigend nach z-Score sortiert. "
         "Reicht die Anzahl auffälliger Gattungen nicht aus, wird mit den "
         "nächst-höchsten z-Scores aufgefüllt. Gattungen mit kuratiertem "
         "Pathogen/Kommensalen-Anteil und nicht-trivialem Read-Support werden "
         "zusätzlich unabhängig vom z-Score forciert aufgenommen (s. z-Score "
         "unten), damit ein starkes Signal nie allein wegen häufiger "
         "Begleitflora im Supplement landet. Zusätzliche Aufnahmebedingung "
         f"nur für Seite 1: dup &lt; {TOP_HITS_MAX_DUP:.0f} (read-gewichteter "
         "Mittelwert über die Gattung) – bei Bakterien/Pilzen/Parasiten, nicht "
         "bei Viren (dort strukturell hohe dup-Werte auch bei echter, tiefer "
         "Coverage kleiner Genome, s. z-Score-Erklärung unten). Gattungen mit "
         "kuratiertem Kontaminanten-Kommentar (s. Spalte \"Kommentar\") sind "
         "unabhängig vom z-Score grundsätzlich von Seite 1 ausgeschlossen und "
         "erscheinen ausschließlich im Supplement – seit die Referenz "
         "Nichtvorkommen als 0 mitzählt (s. z-Score-Erklärung unten), ist "
         "\"selten in der Referenz\" der Normalfall für praktisch jede "
         "Gattung, nicht spezifisch für echte Erreger, sodass auch kuratierte "
         "Kontaminanten regelmäßig hohe z-Scores erreichen können. Enthält "
         "eine Gattung daneben auch eine als Pathogen/Kommensale eingestufte "
         "Spezies, hat diese Vorrang und der Ausschluss greift nicht. "
         "Gattungen, die an einem der beiden Kriterien scheitern, bestehen "
         "weiterhin den Evidenzfilter und erscheinen im Supplement, nur eben "
         "nicht auf Seite 1. Sowohl "
         "Top-Hits als auch Supplement sind zusätzlich nach Kingdom "
         "gegliedert (Viren / Bakterien / Pilze / Parasiten, jeweils mit "
         "Zwischenüberschrift), innerhalb jeder Gruppe weiterhin absteigend "
         "nach z-Score. Alle übrigen Gattungen erscheinen mit denselben "
         "Metriken/Boxplots wie Top-Hits im Supplement, getrennt nach "
         "weiteren seltenen und häufigen (Hintergrund-)Befunden."),
        ("Boxplot (k-mers / dup / Cov.)",
         "Horizontale Referenzverteilung (blassgrau, Log-Skala) aus der >=10kB-"
         f"Kohorte (n = {parsed.get('ref_total', 437)} Reports) auf GENUS-Ebene: pro Referenz-Report werden "
         "Reads/k-mers aller Species/Strains derselben Gattung summiert (Cov. "
         "als Max, dup als Read-gewichteter Mittelwert), daraus Box = "
         "Interquartilsbereich (25.–75. Perzentil), graue Linie = Median, "
         "Whisker = Minimum/Maximum. Roter Punkt = aggregierter Wert dieser "
         "Gattung im aktuellen Fall. \"–\" bedeutet: keine Referenzdaten für "
         "diese Gattung in der >=10kB-Kohorte vorhanden. Die kleine graue "
         "Zahl unter k-mers/Cov. (\"z ...\") ist der z-Score dieser einen "
         "Metrik für sich genommen – Standardabweichungen von log10(k-mers+1) "
         "bzw. log10(Cov.) gegenüber dem Gattungs-Hintergrund, nach derselben "
         "Methodik wie der z-Score in der Hauptspalte (siehe dort), aber ohne "
         "dup-Adjustierung bzw. Kombination der Metriken. So wird sichtbar, "
         "ob eine auffällige Gesamt-Abweichung eher aus breiter Coverage, aus "
         "roher k-mer-Zahl oder aus beidem stammt. Der dup-Boxplot zeigt nur "
         "die Referenzverteilung, ohne eigenen z-Score (dup fließt bereits "
         "in den k-mers-z-Score der Hauptspalte ein, s.u.)."),
        ("z-Score",
         "Standardabweichungen von log10(coverage-adjustierte k-mers + 1) "
         "gegenüber dem üblichen Hintergrundrauschen dieser Gattung in der "
         f"Referenzdatenbank (>=10kB-Kohorte, n = {parsed.get('ref_total', 437)} Reports). "
         "\"Coverage-adjustiert\" heißt: je Species/Strain wird die k-mer-"
         "Zahl durch die Wurzel der Duplikationsrate (dup) geteilt, dann "
         "über die Gattung summiert – k-mers allein trennen echte Erreger "
         "besser von Hintergrundrauschen als Reads (Breite der Genom-"
         "Coverage statt reiner Read-Menge), und die dup-Gewichtung wertet "
         "eine breite, gering-duplizierte Signatur (viele verschiedene "
         "Genompositionen je einmal getroffen) höher als eine schmale, "
         "stark duplizierte (wenige Positionen, oft getroffen). Ausnahme: "
         "bei Viren wird k-mers unverändert übernommen (keine dup-"
         "Abwertung), da hohe dup-Werte dort durch echte tiefe Coverage "
         "kleiner Genome entstehen, nicht durch Kontamination. Positive "
         "Werte = mehr Signal als typischerweise für diese Gattung "
         "beobachtet; z.B. z = +4,6 bedeutet 4,6 Standardabweichungen über "
         "dem üblichen Rauschpegel. Proben, in denen die Gattung NICHT "
         "nachgewiesen wurde, gehen dabei explizit mit dem Wert 0 in "
         "Mittelwert und Streuung ein (nicht nur die Proben mit Nachweis) – "
         "sonst würde bei sehr seltenen Gattungen (z.B. 4 von "
         f"{parsed.get('ref_total', 437)} Fällen) der Mittelwert allein aus "
         "diesen wenigen positiven Nachweisen berechnet und dadurch so hoch "
         "angesetzt, dass ein weiterer, eigentlich seltener Fund fälschlich "
         "unauffällig erschiene. Da praktisch jede Gattung in der Referenz "
         "selten ist (die meisten kommen nur in einer Handvoll der "
         f"{parsed.get('ref_total', 437)} Fälle vor), erreichen dadurch auch "
         "viele kuratierte Kontaminanten regelmäßig einen hohen z-Score – "
         "dieser allein ist deshalb kein verlässliches Ausschlusskriterium "
         "mehr; kuratierte Kontaminanten werden stattdessen unabhängig vom "
         "z-Score generell von Seite 1 ausgeschlossen (s. \"Top-Hits\" "
         "oben). Bei Gattungen ohne jedes Referenz-Vorkommen gilt \"Stille\" "
         "(0 in allen Referenzfällen) als Hintergrund, mit einer über alle "
         "Gattungen gepoolten Streuung als Nenner. Rot hinterlegt: "
         f"z &ge; {Z_THRESHOLD_HIGH:.0f} (starke Abweichung). Orange "
         f"hinterlegt: z &ge; {Z_THRESHOLD_ELEVATED:.1f} (erhöhte Abweichung). "
         "Bestimmt absteigend die Sortierung der Top-Hits-Tabelle."),
        ("Referenz (kleine Zahl unter dem z-Score)",
         f"Anzahl und Anteil der Referenzfälle (>=10kB-Kohorte, n = {parsed.get('ref_total', 437)} "
         "Reports, H. sapiens ausgeschlossen), in denen mindestens ein "
         "Treffer derselben Gattung auftrat. Rein informativ in den Top-Hits "
         "(kein Sortierkriterium); bestimmt im Supplement, ob eine Gattung "
         "als \"weitere seltene\" (< 5%) oder \"häufige\" (&ge; 5%) Befunde "
         "einsortiert wird, dort weiterhin als eigene Spalte \"Referenz\", "
         "rot/orange nach Genus-Prävalenz < 1% / < 5% hinterlegt."),
        ("Kommentar",
         "Pathogen: mindestens ein Species/Strain dieser Gattung ist ein "
         "etablierter ZNS-Erreger (kuratierte Liste). Kommensale: mindestens "
         "einer ist humane Standortflora mit opportunistischem ZNS-"
         "Pathogenpotenzial – erfordert klinische Korrelation. Kontaminant: "
         "typischer Reagenz-/Haut-/Laborkontaminant. Bei mehreren Kategorien "
         "innerhalb derselben Gattung gilt die höchste Priorität (Pathogen > "
         "Kommensale > Kontaminant). Ohne Eintrag: keine der drei kuratierten "
         "Listen trifft zu."),
        ("k-mers",
         "Summe der distinkten 31-mere aller Species/Strains dieser Gattung "
         "im Report. Hohe Diversität bei dup ≈ 1 spricht für genuine "
         "genomische Abdeckung."),
        ("dup",
         "Read-gewichteter Mittelwert der k-mer-Duplikationsrate über alle "
         "Species/Strains dieser Gattung, mit demselben horizontalen "
         "Referenz-Boxplot wie k-mers/Cov. (s.o.). Werte >> 1 deuten auf "
         "repetitive Matches oder Kreuzreaktivität hin – außer bei Viren, "
         "wo hohe dup-Werte durch echte tiefe Coverage kleiner Genome "
         "entstehen (s. z-Score-Erklärung unten)."),
        ("Cov. (%)",
         "Höchster Wert unter den Species/Strains dieser Gattung: Prozentanteil "
         "des jeweiligen Referenzgenoms, der durch gematchte k-mers abgedeckt "
         "wird (1 Nachkommastelle). Nicht summiert, da Coverage genomgrößen-"
         "abhängig ist und zwischen verschiedenen Species nicht additiv ist."),
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

    # ── REFERENZEN ───────────────────────────────────────────────────────
    story.append(HR(3*mm, 3*mm))
    st_ref = S("_ref", fontName="Helvetica", fontSize=7, textColor=C_LIGHT,
               leading=10, leftIndent=-6, spaceAfter=2*mm)
    story.append(P(
        "<b>Klassifikationsmethode:</b> Breitwieser FP, Baker DN, Salzberg SL. "
        "KrakenUniq: confident and fast metagenomics classification using unique "
        "k-mer counts. <i>Genome Biol.</i> 2018;19(1):198. "
        "doi:10.1186/s13059-018-1568-0. PMID: 30445993",
        st_ref))
    story.append(P(
        "<b>Kontaminantenliste:</b> Laurence M, Hatzis C, Brash DE. "
        "Common contaminants in next-generation sequencing that hinder discovery "
        "of low-abundance microbes. <i>PLoS One.</i> 2014;9(5):e97876. "
        "doi:10.1371/journal.pone.0097876. PMID: 24837716 — ergänzt um "
        "etablierte Hautflora- und Umweltkontaminanten.",
        st_ref))

    story.append(P(
        f"<b>Referenzdatenbank:</b> Institutsinterne Prävalenzdatenbank aus "
        f"{parsed.get('ref_total', 437)} KrakenUniq-Reports (>=10kB-Kohorte "
        "kuratierter mNGS-Kontrollfälle ohne Infektionsverdacht, z.B. Tumor-, "
        "Referenzpathologie- und entzündliches Läsionsmaterial gemäß "
        "mNGS-controls.xlsx; H. sapiens ausgeschlossen). Wird bei Vorliegen "
        "neuer Kontrollfälle aktualisiert (scripts/build_reference_db.py "
        "--reports-dir control_cohort).",
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
        s_date = "10. April 2026"
    else:
        parsed = parse_report(report_text)
        s_name = sample_name or "Unbekannte Probe"
        s_date = date or dt_date.today().strftime("%d. %B %Y").lstrip("0")

    counter = _Counter()  # instance variable, safe for batch usage

    def _build(path, draw_fn):
        doc = SimpleDocTemplate(
            path,
            pagesize=A4,
            leftMargin=ML, rightMargin=MR,
            topMargin=MT + 8*mm,   # space for header chrome
            bottomMargin=MB + 8*mm,
            title="Metagenomics Untersuchung",
            author="Institut für Neuropathologie, UKM",
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
    print(f"PDF erstellt: {output_path}  ({counter.total} Seite(n))")


# ─────────────────────────────────────────────────────────────────────────────
# Demo data (Barcode 40)
# ─────────────────────────────────────────────────────────────────────────────
def _demo_parsed():
    """
    Demo-Daten (Species/Strain-Ebene, wie ein echter Report sie liefert).
    Die Genus-Gruppierung, Boxplots und z-Scores werden – wie im echten
    Betrieb – von select_top_hits()/group_by_genus()/attach_genus_stats()
    zur Anzeigezeit berechnet, nicht hier hartkodiert. Ohne eine echte
    analysis/genus_reference_stats.csv-Zeile für die (frei erfundenen)
    Demo-Genus-taxIDs fällt der z-Score auf die "Stille"-Baseline zurück
    (mean=0, gepoolte SD) – für Demozwecke ausreichend.
    """
    demo_ref_total = 437

    def _hit(taxid, name, genus_taxid, genus_name, kingdom, reads, kmers, dup, cov,
             is_contaminant, is_opportunist, is_pathogen, ref_n, ref_pct):
        comment = ("Pathogen" if is_pathogen else
                   "Kommensale" if is_opportunist else
                   "Kontaminant" if is_contaminant else "")
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
            _hit("d1", "Human betaherpesvirus 5", "dg1", "Cytomegalovirus", "Viren",
                 28, 2525, 1.05, 0.01108, False, False, True, ref_n=3, ref_pct=0.25),
            _hit("d2", "Human betaherpesvirus 6B", "dg2", "Roseolovirus", "Viren",
                 2, 4, 13.0, 3.32e-5, False, False, True, ref_n=5, ref_pct=0.42),
            _hit("d3", "Staphylococcus epidermidis", "dg3", "Staphylococcus", "Bakterien",
                 150, 16593, 1.04, 1.75e-3, True, False, False, ref_n=192, ref_pct=16.1),
            _hit("d4", "Staphylococcus hominis", "dg3", "Staphylococcus", "Bakterien",
                 12, 210, 1.4, 6.0e-5, True, False, False, ref_n=192, ref_pct=16.1),
            _hit("d5", "Burkholderia cenocepacia", "dg5", "Burkholderia", "Bakterien",
                 2, 99, 1.0, 4.82e-6, False, False, False, ref_n=283, ref_pct=23.7),
            _hit("d6", "Escherichia coli", "dg6", "Escherichia", "Bakterien",
                 18, 43, 1.63, 6.08e-7, False, True, False, ref_n=416, ref_pct=34.8),
            _hit("d7", "Trichosporon asahii", "dg7", "Trichosporon", "Pilze",
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
    report_date = args.date or _date.today().strftime("%d. %B %Y").lstrip("0")

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
