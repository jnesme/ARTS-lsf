#!/usr/bin/env python
"""
Combine already-completed, independent single-genome ARTS result directories into
the same combined_core_table/combined_known_table/combined_dup_table/summary_table
outputs that ARTS's own native multi-genome mode would produce.

Why this works: ARTS's multi-genome mode (artspipeline1.py call_startquery(), the
len(input_list) > 1 branch) runs each genome through startquery() completely
independently (same static refdir passed to every genome, no cross-genome state),
then calls exactly the functions below, once, on a result_dict mapping each
genome's already-completed result directory to its original input file. Since our
solo runs already produced those same per-genome result directories using the same
refdir/options, calling those functions ourselves reproduces the native multi-genome
output without redundantly re-running the expensive per-genome MAFFT/TrimAl/RAxML
analysis a second time.

Deliberately NOT called here (both require combine_bigscape_results() to have run
first, which only happens when ARTS is invoked with -rbsc/BiG-SCAPE):
  - parse_json()      -- native ARTS itself only calls this wrapped in try/except,
                          because it unconditionally opens combined_bgc_table.tsv,
                          which doesn't exist without BiG-SCAPE.
  - generate_plots()  -- native ARTS calls this UNGUARDED; it also unconditionally
                          opens combined_bgc_table.tsv, so a real multi-genome ARTS
                          run without -rbsc would crash here too. Skipped here since
                          we don't need the plots and don't want a spurious crash.

Usage:
    python combine_smoketest_results.py
"""
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

SMOKETEST_DIR = "/work3/josne/Projects/Vibrio_Galathea3/arts_results/smoketest"

# result_dict maps each genome's already-completed ARTS result directory to the
# original input file used for that genome's run -- exactly what ARTS's own
# call_startquery() builds internally for its multi-genome branch.
result_dict = {
    os.path.join(SMOKETEST_DIR, "vibrio_only"):
        "/work3/josne/Projects/Vibrio_Galathea3/vibrio_seq/funcscan_results/bgc/antismash/S0204/S0204.gbk",
    os.path.join(SMOKETEST_DIR, "pseudoalteromonas_only"):
        "/work3/josne/Projects/Vibrio_Galathea3/pseudoalteromonas_seq/funcscan_results/bgc/antismash/S1608/S1608.gbk",
}

OUTDIR = os.path.join(SMOKETEST_DIR, "merged")

if __name__ == "__main__":
    for res_dir in result_dict:
        if not os.path.isdir(res_dir):
            sys.exit("ERROR: expected completed result directory not found: %s" % res_dir)
        if not os.path.isfile(os.path.join(res_dir, "tables", "coretable.tsv")):
            sys.exit("ERROR: %s has no tables/coretable.tsv -- run did not complete" % res_dir)

    os.makedirs(OUTDIR, exist_ok=True)

    print("Writing combined results to: %s" % OUTDIR)
    write_paths_to_file(result_dict, OUTDIR)
    generate_summary(result_dict, OUTDIR)
    combine_core_results(result_dict, OUTDIR)
    combine_known_results(result_dict, OUTDIR)
    combine_dup_results(result_dict, OUTDIR)

    print("Done. Combined tables written under: %s" % os.path.join(OUTDIR, "tables"))
    for fname in ("summary_table.tsv", "combined_core_table.tsv", "combined_known_table.tsv", "combined_dup_table.tsv"):
        fpath = os.path.join(OUTDIR, "tables", fname)
        print("  %-28s %s" % (fname, "OK" if os.path.isfile(fpath) else "MISSING"))
