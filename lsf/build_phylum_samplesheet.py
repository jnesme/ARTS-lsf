#!/usr/bin/env python
"""
Builds the strain,phylum,gbk_path samplesheet needed to drive ARTS batch
submission (lsf/submit_one_genome.sh requires phylum as a mandatory,
validated argument -- see README "Added tooling").

Taxonomy source: KmerFinder, via each collection's kmerfinder_summary.csv
(bacass repo output). ARTS itself never uses/reports phylum -- it's an
external, pre-invocation decision only (which reference/<phylum>/ set to
run a genome against). See '07-kmerfinder_best_hit_Taxonomy', a
semicolon-delimited lineage string, e.g.:
    cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria; ...
For Proteobacteria, ARTS's reference sets are organized at the CLASS level
(alphaproteobacteria/betaproteobacteria/gammaproteobacteria/
delta_epsilon-proteobacteria), so the class (lineage index 3) is used, not
the phylum (index 2). For every other phylum, ARTS's reference sets are at
the phylum level, so index 2 is used directly.

A handful of samples get zero KmerFinder hits (Total_hits_07_kmerfinder=0,
empty Taxonomy) -- confirmed via manual cross-check against
PRJNA242743_AssemblyDetails.txt's genus-level taxonomy for the specific
cases seen in this project (see KNOWN_GENUS_FALLBACK below). This is a
fallback for these specific known cases, not a general genus->phylum
lookup -- any *new* KmerFinder failure not in this dict is reported as
UNRESOLVED rather than guessed.

Usage:
    python build_phylum_samplesheet.py
Writes: <SMOKETEST_DIR's parent>/arts_phylum_samplesheet.csv
"""
import csv
import os
import sys

ARTS_DIR = "/work3/josne/github/arts"

# Available ARTS reference sets (checked at runtime against what's actually
# on disk, not hardcoded as a static assumption).
def available_refdirs():
    refdirs = set()
    refroot = os.path.join(ARTS_DIR, "reference")
    for name in os.listdir(refroot):
        if os.path.isfile(os.path.join(refroot, name, "coremodels.hmm")):
            refdirs.add(name)
    return refdirs

# Lineage-string phylum/class name -> ARTS reference/<name>/ directory name.
TAXONOMY_TO_REFDIR = {
    "actinobacteria": "actinobacteria",
    "bacteroidetes": "bacteroidetes",
    "chlamydiae": "chlamydiae",
    "cyanobacteria": "cyanobacteria",
    "deinococcus-thermus": "deinococcus-thermus",
    "firmicutes": "firmicutes",
    "fusobacteria": "fusobacteria",
    "spirochaetes": "spirochaetes",
    "tenericutes": "tenericutes",
    # Proteobacteria is split at the class level in ARTS's reference sets:
    "alphaproteobacteria": "alphaproteobacteria",
    "betaproteobacteria": "betaproteobacteria",
    "gammaproteobacteria": "gammaproteobacteria",
    "deltaproteobacteria": "delta_epsilon-proteobacteria",
    "epsilonproteobacteria": "delta_epsilon-proteobacteria",
}

# Known KmerFinder-classification failures (Total_hits_07_kmerfinder=0),
# resolved via PRJNA242743_AssemblyDetails.txt genus-level taxonomy.
# collection -> {strain: (genus, refdir)}
KNOWN_GENUS_FALLBACK = {
    "pseudoalteromonas_seq": {
        "S3726": ("Marinomonas", "gammaproteobacteria"),      # Oceanospirillales
        "S4079": ("Loktanella", "alphaproteobacteria"),        # Rhodobacterales
        "S4388": ("Pseudoalteromonas", "gammaproteobacteria"), # Alteromonadales
        "S2756": ("Pseudoalteromonas", "gammaproteobacteria"),
        "S4382": ("Pseudoalteromonas", "gammaproteobacteria"),
    },
}

COLLECTIONS = [
    {
        "name": "vibrio_seq",
        "kmerfinder_csv": "/work3/josne/Projects/Vibrio_Galathea3/vibrio_seq/Bacass_results_merged/Kmerfinder/kmerfinder_summary.csv",
        "gbk_template": "/work3/josne/Projects/Vibrio_Galathea3/vibrio_seq/funcscan_results_merged/bgc/antismash/{strain}/{strain}.gbk",
    },
    {
        "name": "pseudoalteromonas_seq",
        "kmerfinder_csv": "/work3/josne/Projects/Vibrio_Galathea3/pseudoalteromonas_seq/bacass_results/Kmerfinder/kmerfinder_summary.csv",
        "gbk_template": "/work3/josne/Projects/Vibrio_Galathea3/pseudoalteromonas_seq/funcscan_results/bgc/antismash/{strain}/{strain}.gbk",
    },
]

OUTPUT_CSV = "/work3/josne/Projects/Vibrio_Galathea3/arts_phylum_samplesheet.csv"


def taxonomy_to_refdir(taxonomy_str):
    parts = [p.strip() for p in taxonomy_str.split(";") if p.strip()]
    try:
        bidx = parts.index("Bacteria")
    except ValueError:
        return None, None
    phylum = parts[bidx + 1] if len(parts) > bidx + 1 else None
    klass = parts[bidx + 2] if len(parts) > bidx + 2 else None
    if phylum and phylum.lower() == "proteobacteria" and klass:
        return TAXONOMY_TO_REFDIR.get(klass.lower()), "class:%s" % klass
    elif phylum:
        return TAXONOMY_TO_REFDIR.get(phylum.lower()), "phylum:%s" % phylum
    return None, None


# Preference order when the same strain ID appears in more than one
# collection (confirmed 2026-09-08: 6 such collisions, e.g. S2043/S2757 --
# not the same file, but almost certainly the same physical strain
# assembled twice: once via vibrio_seq's raw-read bacass pipeline, once via
# pseudoalteromonas_seq's pre-assembled/NCBI-sourced entry point. Same
# KmerFinder top-hit species/reference accession in both, but very
# different depth/coverage values and different .gbk checksums. User
# decision: keep the vibrio_seq assembly (full raw-read pipeline), drop
# the pseudoalteromonas_seq duplicate, rather than double-counting the
# strain or arbitrarily disambiguating both into the batch.
COLLECTION_PREFERENCE = ["vibrio_seq", "pseudoalteromonas_seq"]


def main():
    valid_refdirs = available_refdirs()
    rows = []
    unresolved = []
    missing_gbk = []
    unsupported_refdir = []
    duplicate_strains_dropped = []

    for coll in COLLECTIONS:
        name = coll["name"]
        with open(coll["kmerfinder_csv"], newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                strain = row["sample_name"]
                taxonomy = row.get("07-kmerfinder_best_hit_Taxonomy", "")
                refdir, source = taxonomy_to_refdir(taxonomy) if taxonomy.strip() else (None, None)

                if refdir is None:
                    fallback = KNOWN_GENUS_FALLBACK.get(name, {}).get(strain)
                    if fallback:
                        genus, refdir = fallback
                        source = "fallback:AssemblyDetails genus=%s" % genus
                    else:
                        unresolved.append((name, strain))
                        continue

                if refdir not in valid_refdirs:
                    unsupported_refdir.append((name, strain, refdir))
                    continue

                gbk_path = coll["gbk_template"].format(strain=strain)
                if not os.path.isfile(gbk_path):
                    missing_gbk.append((name, strain, gbk_path))
                    continue

                rows.append((name, strain, refdir, source, gbk_path))

    # Resolve cross-collection strain-ID collisions per COLLECTION_PREFERENCE.
    by_strain = {}
    for row in rows:
        by_strain.setdefault(row[1], []).append(row)
    deduped_rows = []
    for strain, dups in by_strain.items():
        if len(dups) == 1:
            deduped_rows.append(dups[0])
            continue
        dups_by_coll = {r[0]: r for r in dups}
        for pref in COLLECTION_PREFERENCE:
            if pref in dups_by_coll:
                kept = dups_by_coll.pop(pref)
                deduped_rows.append(kept)
                break
        for name, r in dups_by_coll.items():
            duplicate_strains_dropped.append((r[0], r[1]))
    rows = deduped_rows

    with open(OUTPUT_CSV, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["collection", "strain", "phylum_refdir", "taxonomy_source", "gbk_path"])
        writer.writerows(rows)

    print("Wrote %d rows to %s" % (len(rows), OUTPUT_CSV))
    print("")
    print("=== Refdir distribution ===")
    counts = {}
    for _, _, refdir, _, _ in rows:
        counts[refdir] = counts.get(refdir, 0) + 1
    for refdir, n in sorted(counts.items(), key=lambda x: -x[1]):
        print("  %-25s %d" % (refdir, n))

    if unresolved:
        print("")
        print("=== UNRESOLVED (no KmerFinder hit and no known fallback) - excluded, needs manual review ===")
        for name, strain in unresolved:
            print("  %s / %s" % (name, strain))

    if unsupported_refdir:
        print("")
        print("=== UNSUPPORTED refdir (taxonomy resolved but no matching reference/<phylum>/ on disk) - excluded ===")
        for name, strain, refdir in unsupported_refdir:
            print("  %s / %s -> %s" % (name, strain, refdir))

    if missing_gbk:
        print("")
        print("=== MISSING antiSMASH .gbk (taxonomy resolved but no gbk file found) - excluded ===")
        for name, strain, path in missing_gbk:
            print("  %s / %s -> %s" % (name, strain, path))

    if duplicate_strains_dropped:
        print("")
        print("=== DUPLICATE strain ID across collections - dropped per COLLECTION_PREFERENCE ===")
        for name, strain in duplicate_strains_dropped:
            print("  dropped %s / %s (kept the other collection's assembly instead)" % (name, strain))

    print("")
    print("Total excluded: %d (%d unresolved, %d unsupported refdir, %d missing gbk, %d duplicate)"
          % (len(unresolved) + len(unsupported_refdir) + len(missing_gbk) + len(duplicate_strains_dropped),
             len(unresolved), len(unsupported_refdir), len(missing_gbk), len(duplicate_strains_dropped)))


if __name__ == "__main__":
    main()
