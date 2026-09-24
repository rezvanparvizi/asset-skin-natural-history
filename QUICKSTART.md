# Quickstart

Ten minutes, in order. Full detail in `docs/server_notes.md`.

## 1. Create the repo on GitHub

Private repo named `asset-skin-natural-history`. Do **not** initialize it
with a README, .gitignore, or license — this scaffold provides them.

## 2. Push this scaffold from your laptop

```bash
cd asset-skin-natural-history
git init -b main
git add -A
git commit -m "Scaffold: three-layer structure, freeze/cohort config, provenance"
git remote add origin git@github.com:<username>/asset-skin-natural-history.git
git push -u origin main
```

(If you have not set up an SSH key on the laptop yet, see
`docs/server_notes.md` → "One-time: SSH key for GitHub".)

## 3. Clone onto the server

VPN first, then:

```bash
ssh <uniqname>@mininubio.ddns.med.umich.edu
mkdir -p ~/projects && cd ~/projects
git clone git@github.com:<username>/asset-skin-natural-history.git
cd asset-skin-natural-history
```

The server needs its own SSH key for GitHub — again see
`docs/server_notes.md`.

## 4. Wire up storage

```bash
df -h                              # find large storage
./setup.sh /path/to/large/storage/$USER/asset
```

This creates the symlinks, the `data/` subdirectories, locks down
`data/clinical/`, makes the scripts executable, and runs the path check.

## 5. Install packages

```bash
Rscript -e 'install.packages("renv"); renv::init()'
Rscript -e 'install.packages(c("yaml","Matrix","ggplot2","dplyr"))'
# then Seurat, BPCells, DESeq2, limma, edgeR, variancePartition, speckle,
# harmony, SingleR as you need them
Rscript -e 'renv::snapshot()'
git add renv.lock && git commit -m "renv: initial snapshot" && git push
```

## 6. Fill in the two things only you can fill in

- `metadata/library_manifest.csv` — from the sample manifest.
  `setup.sh` created it from the template; replace the example rows.
  It is gitignored on purpose.
- `config/freezes/freeze01.yml` — the `TODO` fields: platform,
  chemistry, probe set, CellRanger version, reference. Get these from
  the sequencing core. If this is probe-based Flex, the panel is not
  whole-transcriptome, which constrains which genes you can ask about
  at all.

Put the clinical table at `data/clinical/asset_clinical.csv`. Never
commit it.

## 7. Run the first two analyses

```bash
jobs/run.sh analysis/00_manifest_and_batch/00_1_build_manifest.R
jobs/status.sh

jobs/run.sh analysis/00_manifest_and_batch/00_2_batch_vs_design_crosstab.R
cat results/freeze01/reference/00_manifest_and_batch/batch_vs_design_crosstab/VERDICT.txt
```

The second is hazard #1. **Read its verdict before doing any biology.**
If batch tracks timepoint, the longitudinal comparison at the centre of
Project 1 is confounded and integration will not rescue it.

## 8. Then

Write `analysis/01_sc_qc/01_1_qc.R` and onward. Every script starts:

```r
source("config/paths.R")
source("R/provenance.R")
source("R/io.R")

run <- init_run(stage = "01_sc_qc", run_name = "initial_qc",
                notes = "per-library QC thresholds")
```

and ends:

```r
finalize_run(run, status = "ok", verdict = "one sentence on what you found")
```

---

## Where things are

| I want to... | Go to |
|---|---|
| understand the layout and conventions | `README.md` |
| set up SSH, Positron, learn `screen` | `docs/server_notes.md` |
| know what a manifest column means | `docs/data_dictionary.md` |
| see what I have already run | `docs/runs.csv` |
| remember why I chose a parameter | `docs/decisions.md` |
| record an object from a colleague | `docs/inherited_objects.md` |
| define a new sample subset | `config/cohorts/` |
| add libraries after resequencing | new `config/freezes/freeze02.yml` |
| set up a DE contrast | copy `config/de_runs/EXAMPLE_*.yml` |

## The three rules that matter most

1. **Cluster on everything, infer on subsets.** Layer 1 builds one
   label space from all cells. Layer 2 subsets by cohort and never
   re-clusters.
2. **Integrate on `batch_id`, never on sample, subject, or timepoint.**
   Within-patient change across timepoints is the signal.
3. **Commit before a run you intend to keep.** A run whose `git_sha`
   does not describe its code is not reproducible, and `init_run()`
   will warn you.
