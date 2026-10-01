# ASSET skin scRNA-seq — natural history of early dcSSc

Analysis repository for Project 1: **why does untreated dcSSc skin improve?**

Single-cell RNA-seq of serial skin biopsies (baseline / month 3 / month 6) from
the ASSET trial (NCT02161406), plus simultaneously collected healthy controls,
matched bulk skin RNA-seq, and baseline bulk PBMC RNA-seq.

Maintainer: <your name> · Khanna Lab, University of Michigan Scleroderma Program

---

## Contents

- [Data and privacy](#data-and-privacy)
- [Three-layer analysis structure](#three-layer-analysis-structure)
- [The three indexing axes](#the-three-indexing-axes)
- [Directory map](#directory-map)
- [Naming rules](#naming-rules)
- [Running an analysis](#running-an-analysis)
- [Documentation you must keep current](#documentation-you-must-keep-current)
- [Setup](#setup)

---

## Data and privacy

**No patient-level data belongs in this repository.** Not in a private repo,
not in a branch, not in a notebook output. This repo is hosted on personal
GitHub, which is third-party storage outside U-M's institutional agreements.

| Artifact | In git? | Where it lives |
|---|---|---|
| Code, configs, docs | Yes | this repo |
| `metadata/library_manifest.csv` (technical + design fields) | **No** — gitignored | server only; a `_TEMPLATE` version is tracked |
| Clinical tables (mRSS, autoantibodies, escape dates) | **Never** | server only, under `data/clinical/` |
| Cell-level label tables | No (large) | server only, under `metadata/labels/` |
| Seurat / BPCells objects, counts | No | server large storage, via symlink |
| Small derived summary tables | Yes (opt-in) | `results/**/*.csv` allow-listed in `.gitignore` |

Before committing anything new, ask: *could a row in this file be traced to a
trial participant?* If maybe, it does not go in git.

---

## Three-layer analysis structure

The single most important convention in this repo. Cluster on everything,
infer on subsets.

### Layer 1 — Reference (cohort-agnostic)

QC, ambient handling, integration, lineage assignment, compartment
subclustering, annotation. Runs on **all cells**: both arms, all timepoints,
healthy controls.

Built once per freeze. Its product is a frozen barcode → label table
(`metadata/labels/{freeze}_{labelset}_celltype_labels.csv.gz`).

Rationale: more cells give better cluster definitions and rare-state
detection; and a single common label space is required for any cross-cohort
comparison. If you cluster placebo and abatacept separately, their labels are
not matched and their proportions are not comparable.

**Rules for Layer 1**
- Integrate/correct on **technical batch** (`batch_id`), never on `sample_id`,
  `Subject_ID`, or `timepoint`. Within-patient change across timepoints is the
  signal; patient effects are handled downstream by a random effect, not by
  erasing them from the embedding.
- Choose clustering resolution on marker-driven criteria **before** looking at
  improver status. Record the decision in `docs/decisions.md` with a
  timestamp, then do not revisit it.
- After integration, verify that (a) no cluster is dominated by one batch and
  (b) HC and SSc still separate where biology says they should. Under-
  correction is safer than over-correction here.

### Layer 2 — Inference (cohort-scoped)

Subset cells by cohort, then: composition testing, pseudobulk DE, module
scores, clinical correlation, multimodal integration.

**Never re-cluster in Layer 2.** Use the frozen labels.

### Layer 3 — Exploratory (cohort-scoped re-clustering)

Permitted, e.g. focused re-clustering of placebo fibroblasts to surface a
state diluted in the global model. Must live under `results/{freeze}/exploratory/`,
must be labeled exploratory in `docs/runs.csv`, and must never be the sole
support for a confirmatory claim.

**Circularity warning.** Do not tune clustering until a cluster separates
improvers from non-improvers and then test that cluster for improver
association. That uses the outcome twice and the p-value is meaningless.
Testing composition between patient groups with outcome-blind clusters is
fine — but the blindness has to be real.

---

## The three indexing axes

Every result is identified by three things. Conflating them is how you end up
with `Tcellss_ASSET_byBprop_pseudobulk_2026-05-29`.

| Axis | Meaning | Defined in | Example |
|---|---|---|---|
| **freeze** | which libraries are in the dataset | `config/freezes/*.yml` | `freeze01` (pre-resequencing) |
| **cohort** | which samples this analysis uses | `config/cohorts/*.yml` | `placebo`, `placebo_hc`, `allarms` |
| **labelset** | which annotation version | `metadata/labels/` | `labelset01` |

Freeze and cohort are directory levels. Labelset is **not** — it is a field in
every `run_config.yml` and a column in `docs/runs.csv`, so that
"which runs used the old fibroblast labels?" is one `grep`.

Results path:

```
results/{freeze}/{reference|cohort|exploratory}/{stage}/{run_name}/
```

Example:

```
results/freeze01/placebo/13_sc_pseudobulk_de/fib_improver_vs_non_m0/
├── run_config.yml        # the RESOLVED config, copied in at runtime
├── session_info.txt      # R + package versions
├── git_sha.txt           # code version
├── objects/
├── plots/
└── tables/
```

Every run carries its own parameters, code version, and environment. You will
never again have to reconstruct what `res0p3_from3b` meant.

---

## Directory map

```
config/
  paths.R              all filesystem paths, in ONE place
  freezes/             which libraries are in each freeze
  cohorts/             sample membership rules
  de_runs/             one config per DE contrast
metadata/
  library_manifest.csv       [gitignored] technical + design fields per library
  library_manifest_TEMPLATE.csv
  celltype_dictionary.csv    committed annotation labels + defining markers
  labels/                    [gitignored] barcode -> label tables per labelset
  inherited_labels/          [gitignored] label tables received from colleagues
  gene_sets/                 CD28 module, SASP, Poon T-cell, Ma fibroblast, etc.
R/
  paths_check.R        sanity-check that paths resolve
  provenance.R         init_run() / finalize_run() — the run contract
  io.R                 object loading, BPCells path repair
  qc.R  ambient.R  subcluster.R  pseudobulk.R  composition.R  modules.R
  theme_asset.R        shared ggplot theme
pipelines/
  pseudobulk_de.R      config-driven, standardized output contract
  composition_test.R
analysis/
  00_manifest_and_batch/ ... 17_figures/    numbered stages, run in order
jobs/
  run.sh               launch a script in a detached screen session
  status.sh            list running jobs and tail their logs
docs/
  PROGRESS.md          where the project stands — read first each session
  runs.csv             append-only registry of every run + verdict
  decisions.md         dated log of judgment calls
  inherited_objects.md provenance of anything received from a colleague
  data_dictionary.md   field definitions
  figure_manifest.md   figure number -> script that makes it
  server_notes.md      screen, Positron, git, VPN — practical how-to
data/ results/ figures/ logs/     [symlinks to large storage, gitignored]
```

---

## Naming rules

**Scripts:** `{stage}_{seq}_{verb_phrase}.R`

```
14_1_subset_fibroblasts.R
14_2_sweep_resolution.R
14_3_label_fibroblast_subtypes.R
```

Zero-padded, snake_case, lexical sort equals execution order.

**Run names:** `{compartment}_{contrast}_{scope}` — describe the *question*,
not the date or the parameters.

```
fib_improver_vs_non_m0
myeloid_m0_vs_m3_placebo
adipo_ssc_vs_hc_baseline
```

**Figures:** `fig{N}{panel}_{description}.pdf`

```
fig2c_fibroblast_composition_by_improver.pdf
```

**Never put in a filename:** dates, resolutions (`res0p3`), `v2`, `new`,
`final`, `old`, `test`, `superseded`, your initials. Every one of those is a
job that git, `run_config.yml`, or `docs/runs.csv` does better.

**When a run is superseded:** do not move, rename, or delete it. Set its
`status` in `docs/runs.csv` to `superseded` and put the replacement's `run_id`
in the `verdict` column. Encoding obsolescence in the filesystem produces
archival trees nobody can safely clean up.

---

## Running an analysis

```bash
# on the server, from the repo root
jobs/run.sh analysis/07_sc_subcluster_fibroblast/14_1_subset_fibroblasts.R \
  --freeze freeze01 --cohort reference --labelset labelset01
```

This launches the script in a **detached** `screen` session, logs stdout and
stderr to `logs/`, and returns your prompt immediately. See
`docs/server_notes.md` for reattaching and monitoring.

Inside every analysis script, the first lines are always:

```r
source("config/paths.R")
source("R/provenance.R")

run <- init_run(
  stage    = "07_sc_subcluster_fibroblast",
  run_name = "fib_subcluster_res_sweep",
  notes    = "Resolution sweep 0.1-0.8 on full-cohort fibroblasts"
)
# run$dir, run$objects, run$plots, run$tables are created and ready
# run_config.yml, session_info.txt, git_sha.txt are already written
```

and the last line is:

```r
finalize_run(run, status = "ok",
             verdict = "res=0.3 chosen; 12 clusters, all with distinguishing markers")
```

`finalize_run()` appends the row to `docs/runs.csv`.

---

## Documentation you must keep current

In descending order of how much they will save you:

1. **`docs/runs.csv`** — one row per run, including failures and invalid runs.
   Written automatically by `init_run()`/`finalize_run()`. This is your lab
   notebook and your defense against re-deriving conclusions.
2. **`docs/decisions.md`** — one dated entry per judgment call, always naming
   the alternative you rejected and why. Six months from now a reviewer asks
   why resolution 0.3; this is the answer, already written.
3. **`docs/inherited_objects.md`** — write the entry the day you copy anything
   from a colleague's directory. Especially note BPCells paths.
4. **`config/freezes/*.yml`, `config/cohorts/*.yml`** — these *are*
   documentation, executable rather than prose, which is why they do not rot.

---

## Setup

**Environment is already built** on `mininubio`: R 4.6.1, OpenBLAS 0.3.34,
Bioconductor 3.23, 277 packages locked in `renv.lock`. Read
`docs/ENVIRONMENT.md` before installing anything or debugging a build — it
records six installation workarounds specific to this server.

See `docs/server_notes.md` for the full walkthrough (SSH keys, Positron remote,
VPN, `screen`). Short version, on the server:

```bash
git clone git@github.com:<username>/asset-skin-natural-history.git
cd asset-skin-natural-history

# point the heavy directories at large storage
BIG=/path/to/large/storage/asset
mkdir -p $BIG/{data,results,figures,logs}
ln -s $BIG/data data && ln -s $BIG/results results
ln -s $BIG/figures figures && ln -s $BIG/logs logs

# edit config/paths.R so it matches your server layout, then verify
Rscript R/paths_check.R

# restore the package environment
Rscript -e 'renv::restore()'
```

---

## Related repositories

- `asset-qsp-skin` — Project 2, quantitative systems pharmacology model.
  Deliberately separate: different formalism, dependencies, and timeline. It
  will consume this project's natural-history results as an input.
