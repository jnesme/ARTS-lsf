#!/usr/bin/env python
"""
Combines ARTS results across every completed genome within ONE phylum group
of the full Galathea3 batch run (arts_phylum_samplesheet.csv +
submit_batch_array.sh), producing the same combined_core_table/
combined_known_table/combined_dup_table/summary_table outputs ARTS's own
native multi-genome mode would produce -- see lsf/combine_smoketest_results.py
for the 2-genome proof-of-concept this generalizes, and the README
"Added tooling" / "Interpreting the combined tables" sections for why this
works and how to read the output.

Phylum groups are combined SEPARATELY and must never be mixed: ARTS's
core-gene comparison is refdir-specific (different marker HMMs/gene matrix
per phylum), so a genome run against gammaproteobacteria and one run
against alphaproteobacteria are not comparable in one combined table.

Usage:
    python combine_batch_results.py gammaproteobacteria
    python combine_batch_results.py alphaproteobacteria
"""
import csv
import os
import sys

ARTS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ARTS_DIR)

from artspipeline1 import write_paths_to_file
from combine_results import (
    generate_summary,
    combine_core_results,
    combine_known_results,
    combine_dup_results,
)

SAMPLESHEET = "/work3/josne/Projects/Vibrio_Galathea3/arts_phylum_samplesheet.csv"
BATCH_RESULTS_BASE = "/work3/josne/Projects/Vibrio_Galathea3/arts_results/batch"
COMBINED_OUTPUT_BASE = "/work3/josne/Projects/Vibrio_Galathea3/arts_results/batch_combined"


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: %s <phylum_refdir>  (e.g. gammaproteobacteria, alphaproteobacteria)" % sys.argv[0])
    phylum = sys.argv[1]

    result_dict = {}
    not_yet_done = []
    missing_entirely = []

    with open(SAMPLESHEET, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row["phylum_refdir"] != phylum:
                continue
            strain = row["strain"]
            gbk_path = row["gbk_path"]
            resultdir = os.path.join(BATCH_RESULTS_BASE, phylum, strain)
            coretable = os.path.join(resultdir, "tables", "coretable.tsv")

            if not os.path.isdir(resultdir):
                missing_entirely.append(strain)
                continue
            if not os.path.isfile(coretable) or os.path.getsize(coretable) == 0:
                not_yet_done.append(strain)
                continue

            result_dict[resultdir] = gbk_path

    if not result_dict:
        sys.exit("No completed genomes found for phylum '%s' under %s -- nothing to combine."
                  % (phylum, os.path.join(BATCH_RESULTS_BASE, phylum)))

    outdir = os.path.join(COMBINED_OUTPUT_BASE, phylum)
    os.makedirs(outdir, exist_ok=True)

    print("Combining %d completed genome(s) for phylum '%s'" % (len(result_dict), phylum))
    if not_yet_done:
        print("Skipping %d genome(s) not yet complete (no/empty coretable.tsv): %s"
              % (len(not_yet_done), ", ".join(sorted(not_yet_done))))
    if missing_entirely:
        print("Skipping %d genome(s) with no result directory at all (not submitted/started yet): %s"
              % (len(missing_entirely), ", ".join(sorted(missing_entirely))))

    write_paths_to_file(result_dict, outdir)
    generate_summary(result_dict, outdir)
    combine_core_results(result_dict, outdir)
    combine_known_results(result_dict, outdir)
    combine_dup_results(result_dict, outdir)

    print("")
    print("Done. Combined tables written under: %s" % os.path.join(outdir, "tables"))
    for fname in ("summary_table.tsv", "combined_core_table.tsv", "combined_known_table.tsv", "combined_dup_table.tsv"):
        fpath = os.path.join(outdir, "tables", fname)
        print("  %-28s %s" % (fname, "OK" if os.path.isfile(fpath) else "MISSING"))


if __name__ == "__main__":
    main()
