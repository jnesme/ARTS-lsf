# Antibiotic Resistant Target Seeker (ARTS) Overview

ARTS is a webserver and analysis pipeline for screening for known and putative antibiotic resistance markers in order to identify and prioritize their corresponding biosynthetic gene clusters. 
ARTS allows for specific and efficient genome mining for antibiotics with interesting and novel targets by rapidly linking housekeeping and known resistance genes to BGC proximity, duplication and horizontal gene transfer (HGT) events.

ARTS can be installed locally, or you can use the free public webserver located at https://arts.ziemertlab.com

See https://github.com/ziemertlab/artswebapp for a guide on installing the webserver independently.


# Installation of ARTS

There are three options for installing ARTS:

- Using Docker Images 
- Using Anaconda/Miniconda
- Manual Installation for Linux/Ubuntu

## 1- Using Docker Image:

- Firstly, if you don't have Docker, you should install the Docker engine on your computer. Please check out the latest version of Docker 
[on the official website.](https://docs.docker.com/get-docker/)

- To run ARTS Image, you should download the "docker_run_arts.py" file from the command line or from the repository using a web browser.
```bash
    mkdir ARTSdocker && cd ARTSdocker
    wget https://github.com/ziemertlab/arts/raw/master/docker_run_arts.py
```
**Note:** Python 3.x is needed to run "docker_run_arts.py".

- ARTS Image include only Actinobacteria reference set. If you need other reference sets, please download (~2.2GB) and unzip all of them.
```bash
    mkdir ARTSdocker && cd ARTSdocker
    wget https://arts.ziemertlab.com/static/zip_refsets/all_references.zip
    unzip all_references.zip 
```
- Enter the required arguments and run the script
```bash
    python docker_run_arts.py [-h] [input] [resultdir] [-optional_arguments]
```
- You can see the other details [on Docker Hub](https://hub.docker.com/r/ziemertlab/arts-beta)

## 2- Using Anaconda/Miniconda:
We recommend [Anaconda3/Miniconda3](https://docs.anaconda.com/free/anaconda/install/index.html) (with python >=3.8) and 
it is necessery for the [conda](https://docs.conda.io/en/latest/index.html) package manager.

- Clone/Download the repository (root / sudo required):
```bash
    git clone https://github.com/ziemertlab/arts
```
- Enter the arts folder:
```bash
    cd arts
```
- Prepare Actinobacteria reference set and related HMMs:
```bash
    unzip reference/'*.zip' -d reference/ 
```
- ARTS GitHub repository include only Actinobacteria reference set. If you need other reference sets, please download (~2.2GB) and unzip all of them.
```bash
    wget https://arts.ziemertlab.com/static/zip_refsets/all_references.zip
    unzip all_references.zip 
```
- Create a new environment and install all the packages using the environment.yml file with conda:
```bash
    conda env create -f environment.yml
```
- Activate arts environment:
```bash
    conda activate arts
```
- Install required binary or use pre-compiled linux64bit bin (root / sudo required):
  - Dependency:
    - Ranger-DTL : ranger-dtl-U => http://compbio.mit.edu/ranger-dtl/
  - Pre-compiled bin:
    ```bash
        tar -zxvf linux_64bins.tar.gz -C /usr/local/bin/ ranger-dtl-U
    ```
- Run ARTS (See [Usage](https://github.com/ZiemertLab/ARTS/tree/master#usage) for more):
```bash
    python artspipeline1.py [-h] [input] [refdir] [-optional_arguments]
```

## 3- Manual Installation for Linux/Ubuntu:
The analysis server will start a local antiSMASH job if cluster annotation is not already provided as input. We recommend antiSMASH version >= 6.0.1.
See [antiSMASH](https://docs.antismash.secondarymetabolites.org/install/) for installation instructions.

**Note:** Python version 3.8 or higher is recommended.

- Clone/Download the repository (root / sudo required):
```bash
    git clone https://github.com/ziemertlab/arts
```
- Enter the arts folder:
```bash
    cd arts
```
- Prepare Actinobacteria reference set and related HMMs:
```bash
    unzip reference/'*.zip' -d reference/ 
```
- ARTS GitHub repository include only Actinobacteria reference set. If you need other reference sets, please download (~2.2GB) and unzip all of them.
```bash
    wget https://arts.ziemertlab.com/static/zip_refsets/all_references.zip
    unzip all_references.zip 
```
- Install required libraries and applications (root / sudo required):
```bash
    apt-get update
    apt-get install -y python3-dev liblzma-dev default-jdk hmmer2 hmmer diamond-aligner fasttree prodigal ncbi-blast+ muscle mafft
    pip install -r requirements.txt
```
- Install required binaries or use pre-compiled linux64bit bins (root / sudo required):
  - Dependencies:
    - TrimAl : trimal => https://github.com/inab/trimal
    - RaxML : raxmlHPC-SSE3 => https://github.com/stamatak/standard-RAxML
    - Ranger-DTL : ranger-dtl-U => http://compbio.mit.edu/ranger-dtl/
    - Glimmer : glimmer3 => https://ccb.jhu.edu/software/glimmer/index.shtml
    - GlimmerHMM : glimmerhmm => https://ccb.jhu.edu/software/glimmerhmm/
  - Pre-compiled bins:
    ```bash
        tar -zxvf linux_64bins.tar.gz -C /usr/local/bin/
    ```

- Run ARTS (See [Usage](https://github.com/ZiemertLab/ARTS/tree/master#usage) for more):
```bash
    python artspipeline1.py [-h] [input] [refdir] [-optional_arguments]
```


## Optional: For comparing the results of multi-genome analysis:

The BiG-SCAPE algorithm is used to compare the results of multi-genome analysis. 
All clustered BGCs from antiSMASH results are analyzed to determine BGC similarity. 
The BiG-SCAPE algorithm generates sequence similarity networks of BGCs and classifies them into gene cluster families (GCFs).

To install [the BiG-SCAPE](https://bigscape-corason.secondarymetabolites.org/index.html), please see https://github.com/medema-group/BiG-SCAPE/wiki/installation

**Note:** Make sure that the Pfam database is in the same folder as bigscape.py

# Running ARTS
ARTS uses a webserver to queue jobs to the analysis pipeline. Details on webserver usage can be found at: https://arts.ziemertlab.com/help 

Alternatively jobs can be run directly using the artspipeline1.py script (see -h for options).

````
usage: artspipeline1.py [-h] [-hmms HMMDBLIST] [-khmms KNOWNHMMS] [-duf DUFHMMS] [-cchmms CUSTCOREHMMS] [-chmms CUSTOMHMMS] [-rhmm RNAHMMDB] [-t THRESH]
                        [-td TEMPDIR] [-rd RESULTDIR] [-ast ASTRAL] [-cpu MULTICPU] [-opt OPTIONS] [-org ORGNAME] [-pbt PREBUILTTREES] [-ras]
                        [-asp ANTISMASHPATH] [-bcp BIGSCAPEPATH] [-rbsc]
                        input refdir

Start from genbank file and compare with pre-computed reference for Duplication and Transfers

positional arguments:
  input                 gbk file to start query
  refdir                Directory of precomputed reference files

optional arguments:
  -h, --help            show this help message and exit
  -hmms HMMDBLIST, --hmmdblist HMMDBLIST
                        hmm file, directory, or list of hmm models for core gene id
  -khmms KNOWNHMMS, --knownhmms KNOWNHMMS
                        Resistance models hmm file
  -duf DUFHMMS, --dufhmms DUFHMMS
                        Domains of unknown function hmm file
  -cchmms CUSTCOREHMMS, --custcorehmms CUSTCOREHMMS
                        User supplied core models. hmm file
  -chmms CUSTOMHMMS, --customhmms CUSTOMHMMS
                        User supplied resistance models. hmm file
  -rhmm RNAHMMDB, --rnahmmdb RNAHMMDB
                        RNA hmm models to run (default: None)
  -t THRESH, --thresh THRESH
                        Hmm reporting threshold. Use global bitscore value or Model specific options: gathering= GA, trusted= TC, noise= NC(default: none)
  -td TEMPDIR, --tempdir TEMPDIR
                        Directory to create unique results folder
  -rd RESULTDIR, --resultdir RESULTDIR
                        Directory to store results
  -ast ASTRAL, --astral ASTRAL
                        Location of Astral jar executable default: Value of environment var 'ASTRALJAR'
  -cpu MULTICPU, --multicpu MULTICPU
                        Turn on Multi processing set # Cpus (default: Off, 1)
  -opt OPTIONS, --options OPTIONS
                        Analysis to run. phyl=phylogeny, kres=known resistance, duf=Domain of unknown function, expert=Exploration mode (default: phyl,kres,duf)
  -org ORGNAME, --orgname ORGNAME
                        Explicitly specify organism name
  -pbt PREBUILTTREES, --prebuilttrees PREBUILTTREES
                        Directory of prebuilt trees
  -ras, --runantismash  Run input file through antismash first
  -asp ANTISMASHPATH, --antismashpath ANTISMASHPATH
                        Location of the executable file of antismash or location of antismash 'run_antismash.py' script
  -bcp BIGSCAPEPATH, --bigscapepath BIGSCAPEPATH
                        location of bigscape 'bigscape.py' script
  -rbsc, --runbigscape  Run antismash results through bigscape
````

# Usage 

- For basic run with positional arguments;
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria
````

- To save all output data files: `-rd`, `--resultdir`
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria -rd /PATH/result_folder
````

- To use antiSMASH: `-asp`, `--antismashpath` and to run antiSMASH: `-ras`, `--runantismash`
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria -asp /PATH/antismash -ras -rd /PATH/result_folder
````

- If there is an exsiting antiSMASH job, .json files of antiSMASH results are available fo ARTS: `-asp`, `--antismashpath`
````
    python artspipeline1.py /PATH/antismash_result.json /PATH/arts/reference/actinobacteria -asp /PATH/antismash -rd /PATH/result_folder
````

- To run ARTS with exploration mode, please use `-opt`, `--options` parameter;
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria -asp /PATH/antismash -ras -opt 'expert' 
````

- To identify known resistance, please use `-khmms`, `--knownhmms` and `-opt`, `--options` parameters;
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria -asp /PATH/antismash -ras -khmms /PATH/arts/reference/knownresistance.hmm -opt 'kres'
````

- To identify domain of unknown function(DUF), please use `-duf`, `--dufhmms` and `-opt`, `--options` parameters;
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria -asp /PATH/antismash -ras -khmms /PATH/arts/reference/dufmodels.hmm -opt 'duf'
````

- To run ARTS with phylogeny screening, please use `-ast`, `--astral` and `-opt`, `--options` parameter ;
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria -asp /PATH/antismash -ras -ast /PATH/arts/astral/astral.5.7.7.jar -opt 'phly' 
````

- For multi-genome input, it is enough to put commas without any space between the paths of genome files;
````
    python artspipeline1.py /PATH/input_genome1.gbk,/PATH/input_genome2.gbk,/PATH/input_genome3.gbk /PATH/arts/reference/actinobacteria -rd /PATH/result_folder
````

- To run the BiG-SCAPE algorithms, please use `-bcp`, `--bigscapepath` and `-rbsc`, `--runbigscape`
````
    python artspipeline1.py /PATH/input_genome1.gbk,/PATH/input_genome2.gbk /PATH/arts/reference/actinobacteria -bcp /PATH/BiG-SCAPE_1.1.5/bigscape.py -rbsc -rd /PATH/result_folder
````

# Support
If you have any issues please feel free to contact us at arts-support@ziemertlab.com

# Licence
This software is licenced under the GPLv3. See LICENCE.txt for details.

# Publication
If you found ARTS to be helpful, please [cite us](https://doi.org/10.1093/nar/gkaa374):

Mungan,M.D., Alanjary,M., Blin,K., Weber,T., Medema,M.H. and Ziemert,N. (2020) ARTS 2.0: feature updates and expansion 
of the Antibiotic Resistant Target Seeker for comparative genome mining. 
Nucleic Acids Res.,[10.1093/nar/gkaa374](https://doi.org/10.1093/nar/gkaa374)

Alanjary,M., Kronmiller,B., Adamek,M., Blin,K., Weber,T., Huson,D., Philmus,B. and Ziemert,N. (2017) The Antibiotic 
Resistant Target Seeker (ARTS), an exploration engine for antibiotic cluster prioritization and novel drug target discovery. 
Nucleic Acids Res.,[10.1093/nar/gkx360](https://doi.org/10.1093/nar/gkx360)

# Fork notes (jnesme/ARTS-lsf)

This fork (`git@github.com:jnesme/ARTS-lsf.git`) adds an LSF (DTU HPC) deployment of ARTS for
analyzing the Galathea3 Vibrio/Pseudoalteromonas strain collections, reusing antiSMASH results
already computed by an upstream bacass+funcscan Nextflow pipeline rather than re-running
antiSMASH. Remote convention: `origin` points at this fork (push target), `upstream` points at
[ZiemertLab/ARTS](https://github.com/ziemertlab/arts) (pull future upstream releases via
`git fetch upstream && git merge upstream/master`).

## Added tooling

- `lsf/setup.sh`, `lsf/submit_vibrio_only.sh`, `lsf/submit_pseudoalteromonas_only.sh`,
  `lsf/submit_merged.sh` — bsub submission scripts for the original 2-genome ARTS smoke test on
  the `hpc` queue, feeding antiSMASH `.gbk` output directly into `artspipeline1.py` (no
  `-ras`/antiSMASH re-run). `submit_merged.sh` is kept only as a reference example — its native
  2-genome run is provably redundant with running the two solo scripts and combining afterward
  (see `combine_smoketest_results.py` below), so it should not normally be resubmitted.
- `lsf/submit_one_genome.sh` — the generic, reusable per-genome submission template (`<STRAIN>
  <GBK_PATH> <PHYLUM> <RESULTDIR>`), used for batch-scale runs going forward. **`PHYLUM` is a
  mandatory, validated argument, never defaulted**: ARTS's core-gene comparison is phylum-specific
  (marker HMMs + gene matrix built per phylum reference set), and running a genome against the
  wrong phylum's refdir produces no error — just silently meaningless results. The script checks
  `reference/<phylum>/` actually exists and is complete (`coremodels.hmm`/`genematrix.txt`/
  `model_metadata.json` present), printing the list of valid phyla and exiting 1 on any mismatch.
  This matters concretely for Galathea3: the `pseudoalteromonas_seq` collection is taxonomically
  mixed (most strains Gammaproteobacteria, but ~48/143 are Alphaproteobacteria — Tritonibacter/
  Paracoccus/Loktanella), so phylum can never be safely assumed from source directory alone.
- `lsf/combine_smoketest_results.py` — combines independently-completed single-genome ARTS
  result directories into the same `combined_core_table`/`combined_known_table`/
  `combined_dup_table`/`summary_table` outputs that ARTS's own native multi-genome mode would
  produce. This works because ARTS's multi-genome branch runs every genome through `startquery()`
  completely independently (same static `refdir`, no cross-genome state) and only combines
  already-written per-genome tables at the end — so running genomes as separate jobs and
  combining them afterward reproduces identical results without redundantly repeating each
  genome's expensive MAFFT/TrimAl/RAxML analysis inside one long-running multi-genome job.
- `lsf/build_phylum_samplesheet.py` — builds the `strain,phylum,gbk_path` samplesheet
  `submit_batch_array.sh` (below) consumes, sourced from KmerFinder taxonomy
  (`kmerfinder_summary.csv`, from the sibling `bacass` repo — ARTS itself never uses/reports
  phylum, it's purely an external, pre-invocation choice of which `reference/<phylum>/` set to
  run a genome against). Reduces each genome's semicolon-delimited `Taxonomy` lineage string to
  the class-level token ARTS's reference sets actually use for Proteobacteria
  (alphaproteobacteria/betaproteobacteria/gammaproteobacteria/delta_epsilon-proteobacteria), or
  the phylum-level token directly for every other phylum. Validates every row against what's
  actually on disk (both the `reference/<phylum>/` set and the antiSMASH `.gbk` file) rather than
  assuming, and includes a documented fallback (cross-checked against
  `PRJNA242743_AssemblyDetails.txt`) for the handful of samples KmerFinder fails to classify.
  Also detects and resolves strain-ID collisions across collections (confirmed cases: the same
  physical strain assembled twice via two different pipelines, one raw-reads/bacass, one
  pre-assembled/NCBI-sourced — kept per a fixed collection-preference order, never silently
  double-counted).
- `lsf/submit_batch_array.sh` — the full-batch LSF job-array driver. One array task per
  samplesheet row (looked up via `$LSB_JOBINDEX`), delegating to `submit_one_genome.sh` for the
  actual phylum-validated ARTS invocation. Submit with `bsub < submit_batch_array.sh`; the
  `[1-N]%20` array spec throttles concurrency to be considerate of the shared `hpc` queue (the
  array size `N` must match the samplesheet's row count — update both together if regenerated).
- `lsf/combine_batch_results.py <phylum>` — the batch-scale generalization of
  `combine_smoketest_results.py`: reads the samplesheet, finds every genome for the given phylum
  that has *truly* completed (see the completion-check pitfall documented below — a non-empty
  `coretable.tsv` alone is not sufficient), reports which genomes are skipped and why (not yet
  submitted vs. still running), and combines only the completed ones. Safe to re-run at any point
  while the batch is still in progress — it just combines whatever has finished so far. **Run
  once per phylum group, never combining across phyla** (see "Using ARTS output at batch scale"
  below).
- `lsf/triage_report.py` — the "easy way" batch triage report (flat TSV, one row per
  flagged gene per genome): ranks every truly-completed genome by ARTS's own built-in
  `2+`/`3+` composite score, and for every gene in a genome's own `2+` set, looks up its
  cross-genome recurrence in that phylum's `combined_core_table.tsv`
  (`recur_dup_n`/`recur_bgc_n`/`recur_phyl_n`/`recur_known_n` against `recur_core_n` — see
  "Interpreting the combined tables"). Also tags each gene with a `criteria_tier` column (2 or 3)
  since the `3+` set is a genuine subset of the `2+` set, not a separate list — conflating them
  makes every `2+` gene look equally strong even though only the `3+` ones are the true
  standouts. Does not yet resolve which specific BGC/product a flagged gene sits next to (needs
  parsing `bgctable.tsv`'s nested `Genelist` string) — a natural "harder way" follow-up.

## Bug fixes made in this fork

1. **`NameError: makeantismashresults` on any direct `.gbk` input.** `artspipeline1.py`'s
   non-`-ras` code path unconditionally called `makeantismashresults()` for `.gbk` input, but
   that function is commented out upstream while the call site was left active — crashing
   instantly on exactly the "I already have antiSMASH results" usage this README documents. The
   dead function only ever regenerated a legacy antiSMASH-3.0.5-style HTML page via an internal
   API that no longer exists in antiSMASH ≥5, so the call is now skipped rather than restored.
2. **TrimAl errors on lowercase ambiguity characters MAFFT inserts via `--add`.** When adding a
   new query sequence to a reference alignment, MAFFT can fill newly-created columns in the
   *reference* sequences with a lowercase ambiguity character (e.g. `n`) where it has no data.
   TrimAl's `-automated1` scoring matrix only recognizes uppercase symbols and errors repeatedly
   on lowercase ones (confirmed via direct testing: uppercase `N` and any-case `acgt` are fine,
   only lowercase ambiguity codes trigger it). Sequence lines are now uppercased immediately after
   MAFFT writes its output, before TrimAl reads it.
3. **Silent worker-exception swallowing in parallel tree building.** With `-cpu > 1`,
   `startquery()` dispatches `buildtrees()` across a `multiprocessing.Pool` via `apply_async()`
   without ever collecting the results. Any exception raised inside a worker process (as opposed
   to `buildtrees()`'s own handled False-return path, which *is* logged as `"BuildTree Failed"`)
   was silently discarded — the pool finishes, the pipeline exits 0, and the affected marker's
   core gene tree is simply missing with zero trace anywhere in the log. Fixed by collecting each
   `(marker, AsyncResult)` pair and calling `.get()` on every one after `pool.join()`, logging any
   exception with the marker name attached, plus an explicit
   `"Tree building summary: X/Y markers succeeded"` log line (both the parallel and sequential
   code paths) so a shortfall is always visible rather than requiring a manual file-count
   cross-check.

   Verified three ways: (a) neither original smoke-test run (S0204, S1608) actually hit this —
   `coregenes/*.fna` counts matched `BuildTree`-finished counts exactly in both; (b) a
   fault-injection test against this real (unmodified) `buildtrees()` — 3 genuine markers
   processed normally via real mafft/trimal/raxml calls, alongside one marker wired through a
   wrapper that deliberately raises — confirmed the exception is caught, logged with the correct
   marker name and full traceback, and correctly excluded from the success count (`3/4`, not a
   silently-wrong `4/4`); (c) a 10-genome heavier batch test (5 Vibrio + 5 species-diverse
   Pseudoalteromonas strains, run as independent LSF jobs) produced zero worker exceptions across
   ~900 real marker-tree builds, alongside the same expected ~20-25-per-genome benign
   coverage-threshold filter messages (see below) seen in the original two genomes.

## Incidents in this fork's own batch tooling (not ARTS bugs)

**CRLF line endings silently broke every file-existence check in the first full-batch launch.**
The first 438-genome `submit_batch_array.sh` run failed 438/438 within seconds, every task
reporting `Input GenBank file not found` for a file that demonstrably existed. Root cause:
Python's `csv.writer` defaults to `\r\n` line terminators (per the CSV spec), so
`build_phylum_samplesheet.py`'s output had every row's last field (`gbk_path`) carrying an
invisible trailing carriage return once read back with plain shell `sed`/`cut` in
`submit_batch_array.sh` — the compared path string never matched a real file on disk. Confirmed
this wasn't a filesystem/symlink issue first: both collections failed identically, including
`pseudoalteromonas_seq`'s non-symlinked real files, which ruled that out immediately. Fixed at the
source (`csv.writer(fh, lineterminator="\n")`) and defensively in the array script (`tr -d '\r'`
on the extracted row, so a future hand-edited or differently-generated CSV can't reintroduce
this), then verified with a local `LSB_JOBINDEX=1` dry run reaching the real ARTS invocation
banner, and confirmed for real via a 3-genome small-scale test spanning both collections and both
phyla before relaunching the full batch. No compute was wasted — each failure took ~5 seconds
before any ARTS analysis started, and no result directories were left behind (the file-existence
check runs before `mkdir -p` in `submit_one_genome.sh`).

**A non-empty `coretable.tsv` does not mean a genome has finished.** `writecoretable()`
(`artspipeline1.py`) is called *twice* per genome: once early, before the phylogeny/RangerDTL
step, writing `N/A` placeholders into the `Phylogeny` column, and again at the very end with real
data. `combine_batch_results.py`'s and `triage_report.py`'s first versions both checked only
"`coretable.tsv` exists and is non-empty" as their completion signal — which the early write
already satisfies. Worse, `combine_core_results()` silently treats `N/A` as `"No"`, so combining
a still-running genome doesn't just include incomplete data, it *actively misrepresents* its
phylogeny status as a confirmed negative in the combined table. Caught by direct inspection of
real batch data: of 22 genomes whose `coretable.tsv` looked "done," only 2 actually were — the
other 20 had 100% `N/A` in their `Phylogeny` column and were still mid-tree-building. Fixed in
both scripts by requiring `"Hits with two or more criteria"` to appear in `arts-query.log`
instead — that line is only written after the second, final `writecoretable()` call, making it
the true completion signal for any run using `-opt ...,phyl`.

## Expected (benign) log messages — not bugs

`extractdbgenes.py`'s initial HMM search uses a deliberately loose e-value cutoff (`evalue<=0.1`)
so it doesn't miss weakly-scoring true orthologs, then applies a stricter coverage check
(`genecov>=0.5` or `hmmcov>=0.5`) before accepting a hit as usable in the core-gene alignment. Any
marker whose only hits fail *both* coverage checks gets its placeholder file removed, logged as:

```
INFO - extractdbgenes - None found passing coverage thresholds in potential core: <marker>
```

This is expected, not an error — typically caused by gene fragments truncated at contig edges,
pseudogenes, or divergent paralogs that only weakly cross-hit a profile. Consistently seen at
~20-25 markers per genome (out of ~580-600 initial candidates) across every genome tested in this
fork's smoke test and heavier-batch test, including the same specific marker names recurring
across independent Vibrio genomes (e.g. `TIGR00399`, `TIGR00706`, `TIGR01954` were filtered in
both S0204 and S0276) — reflecting a stable property of specific reference profiles' fit to a
given genus, not a per-genome anomaly. `summary_table.tsv`'s "Core Genes" count (from
`combine_results.py`) reflects the pre-filter candidate count, not the post-filter
`coretable.tsv` row count — the two are expected to differ by roughly this amount.

## Interpreting the combined tables

Two things worth knowing before prioritizing hits from `combined_core_table.tsv`/
`combined_known_table.tsv` (native or via `combine_smoketest_results.py`):

- **`combined_core_table.tsv`'s per-gene fraction columns (Duplication/BGC_Proximity/Phylogeny/
  Known_target) divide by the total organism count in the run, not by how many organisms actually
  have that gene.** A gene found in only 1 of N genomes with a "Yes" flag there shows the same
  fraction as a gene found in all N genomes with only 1 "Yes" — very different signals, identical
  number. Compare `len(Dup_orgs)`/`len(Phyl_orgs)`/etc. against `len(Core_orgs)` (the per-row
  organism-lists also provided) rather than trusting the fraction column alone, especially at
  larger sample sizes where this dilution effect grows.
- **`combined_known_table.tsv`/`combined_dup_table.tsv` row counts are not directly comparable to
  `summary_table.tsv`'s raw per-genome counts.** `combine_known_results()` intentionally collapses
  multiple paralog hits against the same resistance-gene model within one organism to a single
  organism-presence flag before merging across organisms — so its row count can be *smaller* than
  either individual genome's raw `knownhits.tsv` row count. This is correct, organism-presence
  semantics, just different from `combined_core_table.tsv`'s per-gene-fraction approach.

## Using ARTS output at batch scale

ARTS's premise (see its own README intro): a housekeeping gene sitting next to a BGC, especially
if it's also duplicated and/or phylogenetically discordant from the species tree, is a strong
self-resistance-gene candidate — and a strong hint about the BGC's product's mechanism of action.
**The flagged marker genes are ordinary housekeeping gene families (e.g. `ackA` acetate kinase,
`ribB` riboflavin biosynthesis, `sigpep_I_bact` signal peptidase) — not resistance genes by
name.** They're candidates *because of where/how they appear* (duplicated, BGC-adjacent,
phylogenetically odd), not because the gene itself is known to confer resistance. (One partial
exception seen in practice: `TIGR00710`/`efflux_Bcr_CflA` is itself TIGRFAM-classified as a
drug-resistance transporter family, making it a more directly compelling hit than the others.)

### The built-in per-genome shortlist — and what it actually scores

ARTS already computes a composite score for you — no extra tooling needed for this part. Every
completed genome's `arts-query.log` (and, once combined, `summary_table.tsv`'s `2+`/`3+` columns)
contains lines like:

```
INFO - artspipeline1 - Hits with two or more criteria: 16 : {'TIGR01534', 'TIGR01892', ...}
INFO - artspipeline1 - Hits with three or more criteria: 0 : set()
```

(real output from strain F3329, this fork's first Alphaproteobacteria run). **This composite only
combines three of the four `coretable.tsv` columns: Duplication, BGC_Proximity, and Phylogeny**
(`artspipeline1.py:958-967` — `twoplus`/`threeplus` are set intersections over exactly
`rslt["phylogeny"]`, `rslt["proximity"]`, `rslt["duplicates"]`). **`Known_target` (homology to a
characterized resistance/target gene family) is a separate, independent screen and never
contributes to this score.** Confirmed on real data: `TIGR01534` (GAPDH-I) had
Duplication=Yes/Phylogeny=Yes/Known_target=Yes but BGC_Proximity=No — it landed in the `2+` set
(duplicates∩phylogeny) but never the `3+` set, which strictly requires all of
phylogeny∩proximity∩duplicates together, Known_target regardless. Don't assume "3+" means "hit
on any 3 of the 4 columns" — check which 3 specifically.

ARTS's "Known Resistance Hits" (`knownhits.tsv`) is worth checking independently of the `2+`/
`3+` score, not as a component of it — it's matched against ARTS's own curated database of gene
families *documented elsewhere in the literature as natural-product targets* (e.g. DNA gyrase B
as the known target of aminocoumarins), not a general clinical AMR-surveillance database like
CARD/ResFinder — don't expect to see classic named resistance genes (`vanA`, `tetM`, `blaTEM`)
here.

### Batch workflow

1. **Combine per phylum group, separately.** Once genomes in one phylum group have (at least
   partially) completed, run:
   ```
   python lsf/combine_batch_results.py gammaproteobacteria
   python lsf/combine_batch_results.py alphaproteobacteria
   ```
   Each run reports which genomes it included vs. skipped (not yet submitted, or still running)
   and is safe to re-run at any point — it just combines whatever has finished so far. **Never**
   combine across phyla: ARTS's core-gene comparison is refdir-specific (different marker
   HMMs/gene matrix per phylum), so a Gammaproteobacteria genome and an Alphaproteobacteria genome
   are not comparable in one combined table, and combining them would silently produce
   meaningless results — the same underlying principle as `submit_one_genome.sh`'s mandatory
   phylum validation (above): ARTS never checks this for you, so combining across phyla is just
   as silently wrong as running a genome against the wrong `refdir` in the first place.
2. **Triage strains first, not genes.** Sort the combined `summary_table.tsv` by `3+` then `2+`
   descending. This immediately surfaces which *strains* across the whole collection carry the
   richest self-resistance signal, before drilling into any specific gene or BGC.
3. **Drill into a promising strain's own `bgctable.tsv`** to map its `2+`/`3+`-flagged gene IDs
   back to the actual BGC region and predicted product they sit next to.
4. **Cross-check recurrence in `combined_core_table.tsv`.** Does the same gene show up flagged in
   other strains too? Compare `len(Dup_orgs)`/`len(Phyl_orgs)`/etc. against `len(Core_orgs)` for
   that row (see "Interpreting the combined tables" above for why the raw fraction column
   shouldn't be trusted directly) — a gene independently flagged across multiple, phylogenetically
   distinct strains is a far stronger, population-validated candidate than a single-genome hit,
   and is precisely the kind of signal that running many genomes together (rather than one at a
   time) is meant to surface.

## Known unfixed issue (not hit by this fork's usage, documented for awareness)

`combine_results.py`'s `generate_plots()` is called *unguarded* by `artspipeline1.py`'s native
multi-genome mode, but unconditionally opens `combined_bgc_table.tsv` — a file only created by
`combine_bigscape_results()`, which only runs when ARTS is invoked with `-rbsc`/BiG-SCAPE. A
native multi-genome ARTS run without `-rbsc` will therefore crash at this final step (after all
the actually-useful combined tables have already been written). `parse_json()` has the same
unconditional dependency but is at least wrapped in a try/except upstream. Neither function is
called by `lsf/combine_smoketest_results.py` in this fork.
