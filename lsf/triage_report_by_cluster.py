#!/usr/bin/env python
"""
Cluster-deduplicated triage report (see README "Raw gene count vs. distinct
clusters" section for the full writeup).

The built-in 2+/3+ score (and this fork's flat triage_report.py) count each
flagged core gene independently. But ARTS's flagged core genes often occur as
whole multi-gene operons -- most dramatically flagellar biosynthesis genes,
which are duplicated/phylogenetically discordant as a block, for reasons
having nothing to do with antibiotic self-resistance. When several flagged
genes from the same operon land next to the same BGC, the raw score counts
them as independent evidence when they are really one signal.

Confirmed on real data (2026-09-10): in ALL 10 of the batch's top-10 strains
by raw 3+ score, 7-8 of the ~9-10 flagged genes are a single flagellar operon
(either the rod/hook genes FlgB/C/F/G/K/L/M in one BGC, or the export-
apparatus genes FliE/F/FlhB/FliP/FliR/FliQ/FliN/fliI_yscN in another) sitting
next to one BGC -- and 3 of those 10 strains (S2394, S4051, S1193) have ZERO
resistance-annotated gene anywhere in their entire "3+" set once the
flagellar operon and two unrelated housekeeping hits are accounted for.

This script re-ranks by DISTINCT FLAGGED BGC CLUSTERS instead of raw flagged-
gene count, and separately flags whether any of a strain's flagged genes
carries a resistance-plausible TIGRFAM annotation (see README's "Cross-
checking ARTS's candidates against independent ARG calls" section: across
the whole batch, only TIGR00710/efflux_Bcr_CflA does).

Usage:
    python triage_report_by_cluster.py
Writes: /work3/josne/Projects/Vibrio_Galathea3/arts_results/triage_report_by_cluster.tsv
"""
import ast
import csv
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from triage_report import BATCH_BASE, PHYLA, parse_genome_hits  # reuse the same log-parsing logic/bug fixes

ARTS_DIR = "/work3/josne/github/arts"
OUTPUT_TSV = "/work3/josne/Projects/Vibrio_Galathea3/arts_results/triage_report_by_cluster.tsv"

# Keyword scan for resistance-plausible TIGRFAM annotations (see README).
# Confirmed manually 2026-09-09: only 3 total keyword matches across all 347
# genes ever flagged in this batch; 2 are false positives, kept here as an
# explicit, documented exclusion list (same pattern as
# build_phylum_samplesheet.py's KNOWN_GENUS_FALLBACK) rather than a smarter
# but unverified heuristic:
#   TIGR01388 ("rnd")     = ribonuclease D, not an RND efflux pump
#   TIGR02314 (ABC_MetN)  = D-methionine ABC IMPORTER, not a drug-efflux pump
RESISTANCE_KEYWORDS = ["resist", "efflux", "beta-lactam", "betalactam", "multidrug",
                        "mdr", "abc transporter", "rnd", "mfs", "antibiotic", "drug",
                        "vanc", "tetracycl", "aminoglycoside"]
KNOWN_KEYWORD_FALSE_POSITIVES = {"TIGR01388", "TIGR02314"}


def load_model_metadata():
    """phylum -> {tigr_id: [short, desc, category, ...]}"""
    out = {}
    for phylum in PHYLA:
        path = os.path.join(ARTS_DIR, "reference", phylum, "model_metadata.json")
        with open(path) as fh:
            out[phylum] = json.load(fh)
    return out


def is_resistance_plausible(tigr_id, entry):
    if tigr_id in KNOWN_KEYWORD_FALSE_POSITIVES or not entry:
        return False
    text = " ".join(str(x) for x in entry[:3]).lower()
    return any(k in text for k in RESISTANCE_KEYWORDS)


def load_bgc_gene_map(phylum, strain):
    """tigr_id -> list of (cluster_id, cluster_type) it appears as a Core hit in"""
    path = os.path.join(BATCH_BASE, phylum, strain, "tables", "bgctable.tsv")
    gene_map = {}
    if not os.path.isfile(path):
        return gene_map
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 7:
                continue
            cluster_id, ctype = row[0], row[1]
            try:
                genelist = ast.literal_eval(row[6])
            except (ValueError, SyntaxError):
                continue
            for entry in genelist:
                if len(entry) < 2:
                    continue
                gene_map.setdefault(entry[1], []).append((cluster_id, ctype))
    return gene_map


def main():
    metadata = load_model_metadata()
    out_rows = []

    for phylum in PHYLA:
        for resultdir in sorted(glob.glob(os.path.join(BATCH_BASE, phylum, "*/"))):
            strain = os.path.basename(resultdir.rstrip("/"))
            log_path = os.path.join(resultdir, "arts-query.log")
            two_n, three_n, two_genes, three_genes, complete = parse_genome_hits(log_path)
            if not complete:
                continue

            gene_map = load_bgc_gene_map(phylum, strain)

            def cluster_set(genes):
                s = set()
                for g in genes:
                    for cid, ctype in gene_map.get(g, []):
                        s.add((cid, ctype))
                return s

            clusters_2 = cluster_set(two_genes)
            clusters_3 = cluster_set(three_genes)

            resistance_genes = sorted(
                g for g in two_genes
                if is_resistance_plausible(g, metadata[phylum].get(g))
            )

            out_rows.append((
                strain, phylum, two_n, three_n,
                len(clusters_2), len(clusters_3),
                bool(resistance_genes), ",".join(resistance_genes),
                ",".join(sorted(set(ct for _, ct in clusters_3))),
            ))

    # Rank by distinct 3+ clusters desc, then distinct 2+ clusters desc --
    # the cluster-aware analogue of triage_report.py's (-three, -two) sort.
    out_rows.sort(key=lambda r: (-r[5], -r[4]))

    with open(OUTPUT_TSV, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "strain", "phylum", "raw_2plus_genes", "raw_3plus_genes",
            "distinct_2plus_clusters", "distinct_3plus_clusters",
            "has_resistance_annotated_gene", "resistance_annotated_genes",
            "3plus_cluster_types",
        ])
        writer.writerows(out_rows)

    print("Ranked %d truly-completed genome(s) by distinct flagged BGC clusters" % len(out_rows))
    print("Wrote %s" % OUTPUT_TSV)
    print()
    print("Top 15 by distinct 3+ CLUSTERS (not raw gene count):")
    print("%-12s %-22s %6s %6s %8s %8s %-5s %s" % (
        "strain", "phylum", "raw2+", "raw3+", "clust2+", "clust3+", "res?", "resistance genes"))
    for r in out_rows[:15]:
        strain, phylum, two_n, three_n, c2, c3, has_res, res_genes, ctypes = r
        print("%-12s %-22s %6d %6d %8d %8d %-5s %s" % (
            strain, phylum, two_n, three_n, c2, c3, "YES" if has_res else "no", res_genes))


if __name__ == "__main__":
    main()
