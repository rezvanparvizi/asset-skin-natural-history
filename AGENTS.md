# AGENTS.md — orientation

**Read this first.** It exists so that an AI agent or a new collaborator can
understand this repository in two minutes instead of reverse-engineering it
from the filesystem.

---

## What this project is

Single-cell RNA-seq analysis of serial skin biopsies from the **ASSET** trial
(NCT02161406), a phase 2 randomized placebo-controlled trial of abatacept in
early diffuse cutaneous systemic sclerosis (dcSSc).

**Project 1 — the primary question.** mRSS (skin score) improved substantially
in the ASSET **placebo** arm, in patients receiving no treatment at all. Why?
What are the molecular and cellular mechanisms of spontaneous regression, and
what distinguishes untreated patients who improve from those who do not?

**Why it matters.** No background immunomodulatory therapy was permitted in
ASSET, so the placebo arm is a treatment-naive, immunosuppression-free early
dcSSc cohort with serial skin biopsies — effectively unobtainable today.
Natural history is the single largest reason mRSS-endpoint trials in early
dcSSc fail, and the field has almost no direct molecular description of it.

Project 2 (a quantitative systems pharmacology model of SSc skin) is a
background thread and lives in a **separate repository**, `asset-qsp-skin`.

Owner: postdoc in Dinesh Khanna's lab, University of Michigan Scleroderma
Program. Target venue class: Arthritis & Rheumatology / JCI Insight / Annals.

---

## Data

| Modality | Scope |
|---|---|
| scRNA-seq skin (**primary**) | baseline / month 3 / month 6, both arms, plus healthy controls collected simultaneously. ~80 patients, ~240 samples, >1.3M cells. Likely 10x Flex / FFPE. |
| Bulk RNA-seq skin | same biopsies, ~240 samples. Almost certainly GEO **GSE217067**, so intrinsic-subset calls and CD28 module scores already exist. |
| Bulk RNA-seq PBMC | baseline, ~70 samples |
| Clinical | U-M Data Coordinating Center: mRSS 0/1/3/6/9/12 mo, HAQ-DI, ACR CRISS, FVC%, DLCO, autoantibodies, escape-therapy dates, arm, disease-duration stratum |

**Patient-level data never enters git.** Not the clinical table, not the
manifest, not barcode-level label tables. See `.gitignore` and the "Data and
privacy" section of `README.md`.

---

## The three rules that matter most

### 1. Cluster on everything, infer on subsets

Three layers, with different scoping:

- **Layer 1 — Reference (cohort-agnostic).** QC, ambient handling,
  integration, lineage assignment, compartment subclustering, annotation. Runs
  on **all** cells: both arms, all timepoints, healthy controls. Product is a
  frozen barcode → label table.
- **Layer 2 — Inference (cohort-scoped).** Subset by cohort, then composition
  testing, pseudobulk DE, module scores, clinical correlation. **Never
  re-cluster here.**
- **Layer 3 — Exploratory.** Cohort-scoped re-clustering is allowed, must live
  under `results/{freeze}/exploratory/`, must be flagged exploratory in
  `docs/runs.csv`, and must never be the sole support for a claim.

Why: more cells give better cluster definitions and rare-state detection, and a
single common label space is required for any cross-cohort comparison. Cluster
placebo and abatacept separately and their labels are not matched, so their
proportions are not comparable.

### 2. Integrate on `batch_id`, never on sample, subject, or timepoint

Within-patient change across timepoints **is the signal** of Project 1.
`RunHarmony(obj, "orig.ident")` would erase it. Patient effects are handled
downstream by a random effect in the pseudobulk model, not by removing them
from the embedding.

### 3. Do not let the outcome define the clusters

Tuning clustering resolution until a cluster separates improvers from
non-improvers, then testing that cluster for improver association, uses the
outcome twice and makes the p-value meaningless. Choose resolution on
marker-driven criteria, record it in `docs/decisions.md` **before** looking at
improver status, and do not revisit it.

---

## Three indexing axes

Every result is identified by three things. Conflating them is how a project
ends up with `Tcellss_ASSET_byBprop_pseudobulk_2026-05-29`.

| Axis | Meaning | Defined in |
|---|---|---|
| **freeze** | which libraries are in the dataset | `config/freezes/*.yml` |
| **cohort** | which samples an analysis uses | `config/cohorts/*.yml` |
| **labelset** | which annotation version | `metadata/labels/` |

Freeze and cohort are directory levels. Labelset is a **field** in every
`run_config.yml` and a column in `docs/runs.csv`.

```
results/{freeze}/{reference|cohort|exploratory}/{stage}/{run_name}/
├── run_config.yml     resolved parameters, copied in at runtime
├── session_info.txt   R + package versions
├── git_sha.txt        code version
├── objects/ plots/ tables/
```

Every run is self-describing. There is never a need to reconstruct what
`res0p3_from3b` meant.

**Freezes exist because some libraries are being resequenced.** When they
arrive, do not edit `freeze01.yml` — write `freeze02.yml` and re-run. Nothing
is overwritten, and freeze01 vs freeze02 becomes a free robustness check.

---

## Repository map

```
AGENTS.md              this file
README.md              layout, conventions, naming rules
QUICKSTART.md          first-time setup, in order
config/
  paths.R              ALL filesystem paths, in one place
  freezes/             which libraries are in each freeze
  cohorts/             sample membership rules
  de_runs/             one config per DE contrast
metadata/
  library_manifest.csv       [gitignored] one row per library
  celltype_dictionary.csv    valid labels + defining markers
  labels/                    [gitignored] barcode -> label tables
  gene_sets/                 CD28 module, SASP, subtype centroids
R/                     FUNCTIONS only, no top-level execution
  provenance.R           init_run() / finalize_run() — the run contract
  io.R                   config/label loading, BPCells path repair
  paths_check.R          environment sanity check
pipelines/             config-driven, reusable (pseudobulk_de.R)
analysis/              numbered stages, run in order:
  00_manifest_and_batch    metadata, manifest validation, BATCH HAZARD CHECK
  01_bulk_skin             bulk skin RNA-seq (GSE217067): subsets, CD28, trajectory
  02_bulk_pbmc             bulk baseline PBMC RNA-seq
  03_sc_qc  04_sc_ambient  05_sc_integration  06_sc_lineage
  07..11_sc_subcluster_*   fibroblast / immune (all) / lymphoid (T-NK, B-plasma) /
                           myeloid / vascular
  12_sc_assemble_labels    -> the frozen barcode->label table (Layer 1 output)
  13_sc_composition  14_sc_pseudobulk_de
  15_deconvolution         sc reference -> bulk skin (needs 12)
  16_multimodal  17_clinical_models  18_figures
jobs/
  run.sh                 launch a script in a DETACHED screen session
  status.sh              running jobs, load, recent log verdicts
scripts/
  bootstrap_server.sh    R + OpenBLAS + smoke test
  bootstrap_renv.R       project package library
  bootstrap_laptop.sh    macOS dev environment
docs/
  PROGRESS.md            WHERE THE PROJECT STANDS — read this first
  ENVIRONMENT.md         every environment decision and its reasoning
  runs.csv               append-only registry of every run + verdict
  decisions.md           dated log of judgment calls
  inherited_objects.md   provenance of anything from a colleague
  data_dictionary.md     field definitions
  figure_manifest.md     figure -> script that makes it
  server_notes.md        screen, Positron, SSH, VPN
data/ results/ figures/ logs/    [symlinks to large storage, gitignored]
```

---

## Conventions

**Script naming:** `{stage}_{seq}_{verb_phrase}.R` →
`07_3_label_fibroblast_subtypes.R`. Lexical sort equals execution order.

**Run naming:** `{compartment}_{contrast}_{scope}` →
`fib_improver_vs_non_m0`. Describes the *question*, not the date or the
parameters.

**Never in a filename:** dates, resolutions (`res0p3`), `v2`, `new`, `final`,
`old`, `test`, `superseded`, initials. Each of those is a job that git,
`run_config.yml`, or `docs/runs.csv` does better.

**Superseding a run:** do not move, rename, or delete it. Set `status` in
`docs/runs.csv` to `superseded` and put the replacement's `run_id` in
`verdict`.

**Every analysis script** starts with

```r
source("config/paths.R"); source("R/provenance.R"); source("R/io.R")
run <- init_run(stage = "...", run_name = "...", notes = "...")
```

and ends with

```r
finalize_run(run, status = "ok", verdict = "one sentence on what you found")
```

**Statistics defaults:** pseudobulk per sample per cell type with a patient
random effect for DE (never per-cell tests, which treat cells as independent
and inflate significance); proportion models (propeller, scCODA,
Dirichlet-multinomial) for composition. `dream`/`limma-voom` rather than
DESeq2 for any contrast where a subject contributes multiple timepoints —
DESeq2 has no random-effects term.

---

## Environment

Do **not** use the system R. See `docs/ENVIRONMENT.md` for the full reasoning.

```
R            ~/R/R-4.6.1  (source-built; system R is 4.3.2 / Bioc 3.18, stale)
BLAS         ~/opt/openblas  OpenBLAS 0.3.34 (system has only reference BLAS)
Packages     renv, per project, shared cache at ~/.cache/R/renv
             292 installed, 277 in renv.lock
CRAN         pinned to P3M snapshot focal/2026-06-01
Bioconductor 3.23
Pinned       Seurat 5.5.0 · SeuratObject 5.4.0 · Matrix 1.7-5 · harmony 2.0.3
             BPCells 0.3.1  (r-universe; SHA adc4a3c3)
Ceiling      gcc 9.4 -> C++17 only, no C++20
Not present  scDblFinder (C++20 ceiling), MuSiC (TOAST gone from Bioc 3.23)
```

**Launch pattern matters.** Use `bash -c`, NOT `bash -lc` — the `-l` sources
`.bashrc`, which can re-activate conda and shadow system libraries. Keep
`/sbin` on PATH and `USE_BUNDLED_LIBUV=1` set. These are workarounds 1, 3 and
4 in `docs/ENVIRONMENT.md`; read that file's workaround section BEFORE
debugging any package build failure.

Server is `mininubio.ddns.med.umich.edu`: Ubuntu 20.04, 32 cores, 503 GB RAM,
**no root, no job scheduler**. Jobs run under `screen` via `jobs/run.sh`.
Shared box with sustained load 16-18, so thread caps default to 4.

---

## Known hazards, in priority order

These are recorded in `README.md` too, but they govern how any analysis here
should be read.

1. **BATCH.** If prep/capture batch correlates with timepoint or arm, the
   longitudinal comparison is confounded and integration will not rescue it.
   `analysis/00_manifest_and_batch/` runs first, before any biology, and
   writes a `VERDICT.txt`.
2. **ESCAPE THERAPY / INFORMATIVE CENSORING.** 16/44 placebo participants
   started immunomodulators from month 6, *because they were worsening*.
   Month-6 placebo samples are partly drug-exposed and systematically depleted
   of non-improvers. Cohort `placebo` censors them; `placebo_all` does not and
   exists for sensitivity analysis.
   **Resolved for tissue (decisions.md, 2026-10-01):** no Placebo patient
   escaped before the M06 biopsy, so no tissue sample is censored.
3. **TISSUE ENDS AT MONTH 6, OUTCOME AT MONTH 12.** Frame as prediction, not
   as observing tissue at maximal improvement.
4. **FALLING mRSS ≠ REVERSED FIBROSIS.** Dermal atrophy and appendage/adipose
   loss mimic improvement. Single-cell data can discriminate (adipocyte-lineage
   recovery vs loss) — test it, do not assume. Disease duration is the
   separating covariate.
5. **REPEAT-BIOPSY ARTIFACT.** Serial biopsies were taken within ~1 cm of the
   prior site, so months 3 and 6 may carry wound-healing/repair biology.
6. **CONFOUNDERS.** Anti-RNAP3, tendon friction rubs, disease duration,
   baseline mRSS (inclusion caps of 10-35 / 15-45 truncate the range),
   regression to the mean, mRSS inter-rater variability.
7. **POWER.** ~19 placebo improvers vs ~19 non-improvers. Subset-stratified
   analyses within placebo drop to single digits. Say so; do not let an
   exploratory split read as a result.

---

## Prior work on this cohort

Seven papers in project knowledge, cited by number in conversation.

1. **Khanna 2020** (Arthritis Rheumatol 72:125-136) — primary trial. Primary
   endpoint not met: mRSS −6.24 (abatacept) vs −4.49 (placebo), P=0.28.
   Baseline intrinsic subsets (n=84): 33 inflammatory, 33 normal-like, 18
   fibroproliferative.
2. **Chung 2020** (Lancet Rheumatol 2:e743) — open-label extension to month 18.
   Both arms improved further; authors note natural history as a likely
   contributor.
3. **Mehta 2022** (JCI Insight 7:e155282) — bulk skin, GSE217067. CD28
   costimulation elevated at baseline in the inflammatory subset, falls on
   abatacept in improvers. **Directly relevant:** within PLACEBO, normal-like
   patients showed the *opposite* correlation — higher baseline CD28 tracked
   with worsening (r=0.67, P=0.016).
4. **Chakravarty 2015** (Arthritis Res Ther 17:159) — 10-patient pilot,
   GSE66321.
5. **Domsic 2023** (Rheumatology 62:1543) — RNAP3 and TFR predict mRSS
   trajectory. In placebo, RNAP3+ patients declined markedly less than RNAP3−
   (4.52, P=0.03 at 6 months).
6. **Gurrea-Rubio 2026** (Arthritis Rheumatol art.70301) — ASSET PBMC flow
   cytometry; CD319/SLAMF7+ cytotoxic T cells. **Its skin scRNA-seq is a
   separate 8-patient FFPE cohort (GSE343804), NOT this dataset.**
7. **Talia 2021** (J Scleroderma Relat Disord 6:194) — morphea case + CIBERSORT.

**Anchors for Project 1** (all placebo-arm results in this cohort): the Mehta
CD28 finding in (3); the observation in (3) that placebo patients largely stay
in their baseline subset while the field reports transcriptomic
"normalisation"; and the Domsic serology predictors in (5).

---

## Epistemic rules

- Label every mechanistic claim as one of: **(a)** established in human SSc
  skin, **(b)** established in liver/lung fibrosis or mouse models and
  imported as hypothesis, **(c)** observed in this data. Never let (b) be
  stated as (a). **Most of the fibrosis-regression literature is (b).**
- Never produce a citation that has not been verified. If unsure whether a
  paper, statistic, or GEO accession is real, say so.
- Prior ASSET secondary and subset analyses were explicitly
  hypothesis-generating and not multiplicity-corrected. Hold work here to a
  higher standard than that precedent.
- When a question cannot be answered with the available data, say that first,
  then say what would be needed. Do not construct an analysis that technically
  runs but cannot support the claim.
- The owner's knowledge of recent literature exceeds any model's training
  cutoff. For anything recent — including papers from this lab — search or ask
  rather than assert.

---

## Working practices

### Read these first, in this order

1. **`docs/PROGRESS.md`** — where the project stands right now: done, in
   flight, blocked, next, and the open questions that must not be defaulted
   on. This is the only file recording *temporal* state.
2. **`AGENTS.md`** (this file) — structure and rules.
3. **`docs/ENVIRONMENT.md`** — before running or installing anything. Its
   workaround section lists six server-specific failures with symptoms; check
   it before debugging any build error.
4. **`docs/runs.csv`** — before proposing an analysis. It may already be done,
   or already marked `invalid`.
5. **`docs/decisions.md`** — before changing a parameter. The current value
   may be deliberate and justified.

### One task per session

Keep a session scoped to one analysis stage or one bounded piece of work.
Stage directories under `analysis/` are the natural units.

The reason is context degradation, not tidiness. In a long session the risk
shifts from "does not understand the task" to "optimises for consistency with
earlier turns rather than correctness" — which produces confident, wrong,
internally-coherent output. A concrete instance: during the environment build,
a migration script was shipped that had rewritten its own rename table with
`sed`, and the error was not caught because the session had been running for
hours.

When a session is done, update `docs/PROGRESS.md`. Two minutes.

### Handshake before implementing

For anything involving an analytical choice — a statistical design, a
clustering decision, a cohort definition, an outcome variable — state the plan
and get agreement before writing code:

1. Restate the goal in one or two sentences.
2. Say what will be computed, on which cohort, with which model.
3. Name the assumption most likely to be wrong.
4. Wait for confirmation.

Skip the handshake for mechanical work (installing a package, fixing a typo,
renaming a file).

The asymmetry that justifies this: a build error announces itself, whereas a
wrong statistical framing produces code that runs cleanly and answers a
different question than the one asked. In this project the second failure mode
is the expensive one.

### Decisions get recorded, not just made

`docs/decisions.md` is the decision log — the equivalent of Architecture
Decision Records, kept as one dated append-only file rather than one file per
decision. Every entry names the alternative rejected and why.

Add an entry when: a threshold or resolution is chosen; a covariate is
included or excluded; an outcome is defined; a package version is pinned; an
inherited object is adopted or rebuilt.

### Correctness checks, not unit tests

Analysis code does not benefit from unit tests the way application code does.
What it does benefit from is assertions on the properties that matter, run
inside the pipeline:

- pseudobulk aggregation preserves total counts per sample
- a label join achieves complete barcode coverage — fail loudly, never
  silently produce `NA`s
- the design matrix is full rank before any model is fitted
- one pseudobulk unit maps to exactly one subject, timepoint and batch
- cell counts per unit exceed the declared minimum

`pipelines/pseudobulk_de.R` already implements several of these. Add to them
rather than importing a testing framework.

### If you are an AI agent

- Never write patient-level data into the repository.
- Never run heavy computation in the foreground on the server; use
  `jobs/run.sh`, which detaches into `screen`.
- Use `init_run()` / `finalize_run()` for anything that produces output, so
  provenance is recorded automatically.
- Respect the three-layer rule: do not re-cluster in a cohort-scoped analysis.
- Do not invent citations, statistics or GEO accessions. If unsure whether
  something is real, say so.
- Label mechanistic claims (a) established in human SSc skin, (b) imported
  from liver/lung fibrosis or mouse models as hypothesis, or (c) observed in
  this data. Most of the fibrosis-regression literature is (b).
- If a script you produced has a bug, say so plainly rather than letting the
  owner assume operator error.


## Rules for AI agents working in this repository

- Never read, print, summarise or search files under data/,
  metadata/clinical/, metadata/labels/, metadata/inherited_labels/, or
  metadata/library_manifest.csv, by any route, including shell commands
  such as cat, head, less, grep, or R/Python code that loads them.
- If a task seems to require patient-level data, stop and ask. Work from
  the column names or the TEMPLATE files instead.
- Do not propose moving patient-level tables into results/ or any other
  readable folder. The one exception is metadata/design/ (owner's
  decision, docs/decisions.md 2026-10-01): AI agents may read the subject
  design table there (Subject_ID, arm, mRSS_category, Ever_escaped,
  Escape_month, timepoints and modalities per subject), the ID-match
  report, and the bulk file-name inventory. Nothing else from the
  clinical table goes there without a new decision.
- Heavy computation never runs in the foreground: use jobs/run.sh.
- Follow the three-layer rule and the integration rule stated above.

