#!/usr/bin/env python
"""
Cross-checks ARTS's top candidate resistance-gene family (TIGR00710,
efflux_Bcr_CflA -- present in nearly every genome per ARTS's broad HMM, see
README "Using ARTS output at batch scale") against funcscan's independent,
per-tool ARG calls for the specific Bcr/CflA-family gene in this locus:
tet(35) ("tetracycline efflux Na+/H+ antiporter family transporter Tet(35)"),
which abricate/AMRFinderPlus/RGI call in near-perfect 3-way agreement at the
same genomic coordinates.

An earlier, narrower version of this check searched hamronization output for
the literal strings "bcr"/"cfla" and found only DeepARG's low-confidence
"potential ARG" pass flagging a hit under that name -- a single-tool, weakly
supported result. tet(35) is the real, strongly cross-validated candidate at
this locus (confirmed 2026-09-09 by checking genomic coordinates: DeepARG's
"Bcr" locus and the abricate/amrfinderplus/rgi "tet(35)" locus are DIFFERENT
loci in the same genome -- this script checks tet(35), not the DeepARG name).

Data source notes (each a real gotcha hit while building this):
  - vibrio_seq has THREE hamronization_combined_report.tsv files on disk
    (funcscan_results, funcscan_results_merged, funcscan_results_test5/batch2).
    funcscan_results_merged is the current, complete one (532 distinct
    strains, dated 2026-08-24) -- it's also what
    build_phylum_samplesheet.py's gbk_template points ARTS at for vibrio_seq.
    funcscan_results (370 strains, dated 2026-03-04) is a stale partial run;
    using it silently drops >160 strains from this cross-check.
  - pseudoalteromonas_seq's funcscan_results (478 strains, 2026-09-04) is its
    own separate, correct report -- same hamronization column schema,
    confirmed via diff before combining.
  - 6 strain IDs appear in both collections (same physical strain assembled
    twice -- see build_phylum_samplesheet.py's COLLECTION_PREFERENCE). ARTS
    was run on the vibrio_seq assembly for these, so this script must use
    vibrio_seq's ARG calls for them too, not pseudoalteromonas_seq's --
    otherwise the ARG evidence and the ARTS result being compared could come
    from two different physical assemblies of the same strain.

Usage:
    python argcheck_crosscheck.py
Writes: /work3/josne/Projects/Vibrio_Galathea3/arts_results/tet35_argcheck_crosscheck.tsv
"""
import csv
import os

HAMRONIZATION_SOURCES = [
    ("vibrio_seq", "/work3/josne/Projects/Vibrio_Galathea3/vibrio_seq/funcscan_results_merged/reports/hamronization_summarize/hamronization_combined_report.tsv"),
    ("pseudoalteromonas_seq", "/work3/josne/Projects/Vibrio_Galathea3/pseudoalteromonas_seq/funcscan_results/reports/hamronization_summarize/hamronization_combined_report.tsv"),
]
COLLECTION_PREFERENCE = ["vibrio_seq", "pseudoalteromonas_seq"]

TRIAGE_REPORT = "/work3/josne/Projects/Vibrio_Galathea3/arts_results/triage_report.tsv"
OUTPUT_TSV = "/work3/josne/Projects/Vibrio_Galathea3/arts_results/tet35_argcheck_crosscheck.tsv"

TARGET_GENE_SYMBOL = "tet(35)"


def strain_from_input_file_name(name):
    # hamronization's input_file_name bakes in the calling tool's own file
    # suffix (e.g. "S2687.tsv.amrfinderplus", "S2687.txt.rgi", "S2687" for
    # abricate) -- the strain id is always the leading dot-separated token.
    return name.split(".")[0]


def load_tet35_calls():
    """strain -> {"S1234"} sorted set of analysis_software_name values that called tet(35)"""
    calls_by_collection = {}
    for coll_name, path in HAMRONIZATION_SOURCES:
        calls = {}
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("gene_symbol", "") != TARGET_GENE_SYMBOL:
                    continue
                strain = strain_from_input_file_name(row["input_file_name"])
                calls.setdefault(strain, set()).add(row["analysis_software_name"])
        calls_by_collection[coll_name] = calls

    merged = {}
    for strain in set().union(*calls_by_collection.values()):
        by_coll = {c: calls_by_collection[c][strain] for c in calls_by_collection if strain in calls_by_collection[c]}
        for pref in COLLECTION_PREFERENCE:
            if pref in by_coll:
                merged[strain] = by_coll[pref]
                break
        else:
            # strain present under a collection name not in COLLECTION_PREFERENCE
            merged[strain] = next(iter(by_coll.values()))
    return merged


def load_arts_ranking():
    """strain -> (phylum, two_plus, three_plus); also returns the full rank order list"""
    seen = {}
    with open(TRIAGE_REPORT, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            strain = row["strain"]
            if strain in seen:
                continue
            seen[strain] = (row["phylum"], int(row["strain_2plus_total"]), int(row["strain_3plus_total"]))
    ranked = sorted(seen.items(), key=lambda kv: (-kv[1][2], -kv[1][1]))
    rank_of = {strain: i + 1 for i, (strain, _) in enumerate(ranked)}
    return seen, rank_of, len(ranked)


def main():
    tet35_calls = load_tet35_calls()
    arts_by_strain, rank_of, total_ranked = load_arts_ranking()

    rows = []
    for strain in sorted(tet35_calls):
        tools = sorted(tet35_calls[strain])
        if strain in arts_by_strain:
            phylum, two_n, three_n = arts_by_strain[strain]
            rank = rank_of[strain]
            rows.append((strain, len(tools), ",".join(tools), phylum, rank, total_ranked, two_n, three_n))
        else:
            rows.append((strain, len(tools), ",".join(tools), "", "", total_ranked, "", ""))

    rows.sort(key=lambda r: (r[4] if r[4] != "" else 10**9, r[0]))

    with open(OUTPUT_TSV, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "strain", "n_tools_confirming_tet35", "tools_confirming_tet35",
            "arts_phylum", "arts_rank", "arts_rank_out_of",
            "arts_2plus_total", "arts_3plus_total",
        ])
        writer.writerows(rows)

    resolved = [r for r in rows if r[4] != ""]
    unresolved = [r for r in rows if r[4] == ""]

    print("tet(35)-positive strains (>=1 tool, either collection): %d" % len(rows))
    print("  with a completed ARTS rank: %d" % len(resolved))
    print("  not yet in the ARTS batch / not completed: %d" % len(unresolved))
    if unresolved:
        print("    %s" % ", ".join(r[0] for r in unresolved))
    if resolved:
        ranks = [r[4] for r in resolved]
        print("")
        print("ARTS rank of tet(35)-confirmed strains (out of %d ranked genomes):" % total_ranked)
        print("  best (lowest/top) rank : %d" % min(ranks))
        print("  worst rank             : %d" % max(ranks))
        print("  mean rank              : %.1f" % (sum(ranks) / len(ranks)))
        n_3tool = sum(1 for r in resolved if r[1] == 3)
        print("  3-tool agreement (abricate+amrfinderplus+rgi): %d/%d" % (n_3tool, len(resolved)))
    print("")
    print("Wrote %s" % OUTPUT_TSV)


if __name__ == "__main__":
    main()
