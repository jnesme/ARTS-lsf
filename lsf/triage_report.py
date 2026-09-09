#!/usr/bin/env python
"""
"Easy way" batch triage report (see README "Using ARTS output at batch scale").

Ranks every completed genome by ARTS's own built-in composite score (the
"Hits with two or more / three or more criteria" set already computed per
genome in arts-query.log -- no re-derivation needed), then for each flagged
gene looks up its cross-genome recurrence in that phylum's
combined_core_table.tsv (len(Dup_orgs)/len(Phyl_orgs)/len(BGC_orgs)/
len(Known_orgs) against len(Core_orgs) -- NOT the raw fraction column, see
README "Interpreting the combined tables" for why).

This is deliberately the simple, flat-TSV version: one row per
(strain, flagged gene). It does not yet resolve which specific BGC/product
each flagged gene sits next to (that needs parsing bgctable.tsv's nested
Genelist string) -- a natural "harder way" follow-up once this is useful.

Usage:
    python triage_report.py
Writes: /work3/josne/Projects/Vibrio_Galathea3/arts_results/triage_report.tsv
"""
import ast
import csv
import glob
import os
import re
import sys

# combined_core_table.tsv's [Hits_listed] column concatenates every sharing
# organism's Genelist string -- for a widely-shared core gene this grows with
# genome count and exceeds Python's default 128KB csv field limit once enough
# genomes are combined (hit in practice at n=194; not at n=2/22/65). We never
# read that column here (see load_combined_core_lookup), so it's safe to just
# raise the limit rather than change what's parsed. sys.maxsize can overflow
# the platform's C long for this call on some systems -- back off until it fits.
_limit = sys.maxsize
while True:
    try:
        csv.field_size_limit(_limit)
        break
    except OverflowError:
        _limit = int(_limit / 10)

BATCH_BASE = "/work3/josne/Projects/Vibrio_Galathea3/arts_results/batch"
COMBINED_BASE = "/work3/josne/Projects/Vibrio_Galathea3/arts_results/batch_combined"
PHYLA = ["gammaproteobacteria", "alphaproteobacteria"]
OUTPUT_TSV = "/work3/josne/Projects/Vibrio_Galathea3/arts_results/triage_report.tsv"

HITS_RE = re.compile(r"Hits with (two|three) or more criteria: (\d+) : (\{.*\}|set\(\))")


def load_combined_core_lookup(phylum):
    """core_gene_id -> {Core_orgs, Dup_orgs, BGC_orgs, Phyl_orgs, Known_orgs} lengths"""
    path = os.path.join(COMBINED_BASE, phylum, "tables", "combined_core_table.tsv")
    lookup = {}
    if not os.path.isfile(path):
        return lookup
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            if len(row) < 13:
                continue
            gene_id = row[0]
            try:
                core_orgs = ast.literal_eval(row[8])
                dup_orgs = ast.literal_eval(row[9])
                bgc_orgs = ast.literal_eval(row[10])
                phyl_orgs = ast.literal_eval(row[11])
                known_orgs = ast.literal_eval(row[12])
            except (ValueError, SyntaxError):
                continue
            lookup[gene_id] = {
                "core_n": len(core_orgs),
                "dup_n": len(dup_orgs),
                "bgc_n": len(bgc_orgs),
                "phyl_n": len(phyl_orgs),
                "known_n": len(known_orgs),
            }
    return lookup


def parse_genome_hits(arts_query_log):
    """
    Returns (two_plus_n, three_plus_n, two_plus_genes:set, three_plus_genes:set, complete:bool).
    complete=False means "Hits with two or more criteria" was never logged --
    writecoretable() writes an early N/A-placeholder version of coretable.tsv
    before the phylogeny/RangerDTL step, and only the final write (after this
    log line) reflects real data. A genome with no such line yet is still
    mid-run, not a genuine zero-hit result -- must not be silently treated
    as 0/0 (confirmed as a real issue on real batch data: 20/22 "completed"
    coretable.tsv files were actually this early, incomplete write).

    three_plus_genes is a SUBSET of two_plus_genes (ARTS logs both sets
    independently) -- without tracking it separately, every 2+ gene looks
    equally strong even though only the 3-criteria ones are the true
    standouts (confirmed on real data: S4733 had 48 genes at 2+ but only 3
    of those specific genes were actually at 3+).
    """
    two_plus_n, three_plus_n = 0, 0
    two_plus_genes, three_plus_genes, complete = set(), set(), False
    if not os.path.isfile(arts_query_log):
        return two_plus_n, three_plus_n, two_plus_genes, three_plus_genes, complete
    with open(arts_query_log) as fh:
        for line in fh:
            m = HITS_RE.search(line)
            if not m:
                continue
            kind, n, genes_repr = m.groups()
            try:
                genes = ast.literal_eval(genes_repr)
            except (ValueError, SyntaxError):
                # Reading a log file that a still-running genome is actively
                # appending to can catch a line mid-flush -- rare but real,
                # confirmed once against a live batch. Treat as "not this
                # line" rather than crashing the whole report; the next run
                # (this script is meant to be re-run repeatedly against a
                # live batch) will see the completed line instead.
                continue
            if kind == "two":
                two_plus_n, two_plus_genes, complete = int(n), genes, True
            else:
                three_plus_n, three_plus_genes = int(n), genes
    return two_plus_n, three_plus_n, two_plus_genes, three_plus_genes, complete


def main():
    strain_rows = []  # (strain, phylum, two_plus_n, three_plus_n, two_plus_genes, three_plus_genes)
    skipped_incomplete = []
    for phylum in PHYLA:
        for resultdir in sorted(glob.glob(os.path.join(BATCH_BASE, phylum, "*/"))):
            strain = os.path.basename(resultdir.rstrip("/"))
            coretable = os.path.join(resultdir, "tables", "coretable.tsv")
            if not os.path.isfile(coretable) or os.path.getsize(coretable) == 0:
                continue
            log_path = os.path.join(resultdir, "arts-query.log")
            two_n, three_n, two_genes, three_genes, complete = parse_genome_hits(log_path)
            if not complete:
                skipped_incomplete.append((strain, phylum))
                continue
            strain_rows.append((strain, phylum, two_n, three_n, two_genes, three_genes))

    strain_rows.sort(key=lambda r: (-r[3], -r[2]))

    combined_lookup_cache = {p: load_combined_core_lookup(p) for p in PHYLA}

    with open(OUTPUT_TSV, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "strain", "phylum", "strain_2plus_total", "strain_3plus_total",
            "flagged_gene", "criteria_tier", "recur_core_n", "recur_dup_n", "recur_bgc_n",
            "recur_phyl_n", "recur_known_n",
        ])
        for strain, phylum, two_n, three_n, two_genes, three_genes in strain_rows:
            lookup = combined_lookup_cache[phylum]
            if not two_genes:
                writer.writerow([strain, phylum, two_n, three_n, "", "", "", "", "", "", ""])
                continue
            for gene in sorted(two_genes):
                rec = lookup.get(gene, {})
                tier = 3 if gene in three_genes else 2
                writer.writerow([
                    strain, phylum, two_n, three_n, gene, tier,
                    rec.get("core_n", ""), rec.get("dup_n", ""),
                    rec.get("bgc_n", ""), rec.get("phyl_n", ""), rec.get("known_n", ""),
                ])

    print("Ranked %d truly-completed genome(s) (%d Gammaproteobacteria, %d Alphaproteobacteria)"
          % (len(strain_rows),
             sum(1 for r in strain_rows if r[1] == "gammaproteobacteria"),
             sum(1 for r in strain_rows if r[1] == "alphaproteobacteria")))
    if skipped_incomplete:
        print("Skipped %d still-mid-run genome(s) (coretable.tsv exists but phylogeny step "
              "not finished yet): %s"
              % (len(skipped_incomplete), ", ".join(s for s, _ in skipped_incomplete)))
    print("Wrote %s" % OUTPUT_TSV)
    print("")
    print("Top 10 strains by 3+/2+ criteria:")
    print("%-12s %-22s %6s %6s" % ("strain", "phylum", "2+", "3+"))
    for strain, phylum, two_n, three_n, _, _ in strain_rows[:10]:
        print("%-12s %-22s %6d %6d" % (strain, phylum, two_n, three_n))


if __name__ == "__main__":
    main()
