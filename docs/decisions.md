# Decisions log

One dated entry per judgment call. Always name the alternative you
rejected and why. Six months from now a reviewer asks why resolution
0.3; this file is the answer, already written.

**Format.** Newest at the bottom, so the file reads as a chronology.
Keep entries short. Link to the run directory or run_id that supports
the decision.

**What belongs here**
- clustering resolutions and the criteria used to pick them
- QC thresholds, and what they removed
- integration variable choices
- cell-type label assignments and merges
- outcome variable definitions
- covariates included or excluded from a model
- anything you would struggle to justify from the code alone

**What does not belong here** — parameters already captured in
`run_config.yml`, and results (those go in `docs/runs.csv`).

---

## 2026-09-23 — Repository structure: three layers

Adopted cluster-on-everything / infer-on-subsets.

- Layer 1 (reference): QC, ambient, integration, subclustering,
  annotation on ALL cells (both arms, all timepoints, healthy
  controls). Product is a frozen barcode -> label table.
- Layer 2 (inference): cohort-scoped composition, DE, module scores,
  clinical correlation. No re-clustering.
- Layer 3 (exploratory): cohort-scoped re-clustering, labelled
  exploratory, never the sole support for a claim.

Rejected clustering each cohort separately: labels would not be matched
across cohorts, so proportions would not be comparable, and the
abatacept-arm-as-validation option would be foreclosed. Since the final
cohort for the paper is not yet decided, a single common label space is
the only choice that keeps options open.

## 2026-09-23 — Indexing axes: freeze / cohort / labelset

freeze = which libraries are in the dataset; cohort = which samples an
analysis uses; labelset = which annotation version. Freeze and cohort
are directory levels; labelset is a field in `run_config.yml` and a
column in `docs/runs.csv`.

Rejected date-stamped output directories (the pattern in the
colleague's tree): four near-identical copies of one analysis
distinguishable only by date, with no record of what changed between
them. The date was doing a job it cannot do.

## 2026-09-23 — Pseudobulk engine depends on the contrast

DESeq2 for cross-sectional contrasts (one pseudobulk unit per subject).
limma-voom with `duplicateCorrelation`, or `variancePartition::dream`,
for any contrast where a subject contributes multiple timepoints.

Reason: DESeq2 has no random-effects term, so repeated measures are
treated as independent and significance is inflated. The colleague's
pipeline uses DESeq2 throughout, which is appropriate for her mostly
cross-sectional contrasts but not for the longitudinal placebo
comparisons at the centre of this project. Output contract kept
identical so results remain comparable.

---

<!-- TEMPLATE — copy below this line

## YYYY-MM-DD — <short title>

<Decision, in one or two sentences.>

Chose: <what>
Rejected: <alternative> because <reason>
Evidence: <run_id or results path>
Revisit when: <condition, e.g. "freeze02 adds libraries">

-->

## 2026-09-29 — AI coding assistant: data-access rules

Decision: Claude Code (personal subscription, VS Code extension in Positron)
is permitted on the server for code, logs and results. It is blocked from
patient-level data.

Convention: patient-level tables (clinical, mRSS, autoantibodies, escape
therapy, the real library manifest, barcode-level labels) live ONLY in
data/, metadata/clinical/, metadata/labels/ or metadata/inherited_labels/.
Never in results/ or anywhere else.

Enforcement: .claude/settings.json denies Claude's file tools on those
folders, their real paths under /home/parvizi/asset-data/, and the /hits
protection copy. Shell commands require manual approval; refuse any that
touch those paths.

Limitation: deny rules are a guardrail, not a security boundary. The
folder convention is what makes them work.



## 2026-09-29 — Raw data: copy, not symlink

Decision: every library's CellRanger output is copied into data/raw/
by analysis/00_manifest_and_batch/00_0_import_raw_libraries.R, md5-
verified against the source, and made read-only. Checksums in
data/raw/checksums.csv. Per library: sample_filtered and sample_raw
feature matrices (.h5) and metrics_summary.csv; per pool: the pool raw
matrix and config.csv. BAMs and molecule_info are not copied.

Rejected: symlinking to /hits. The upstream files belong to other
people and have already changed once (wasikowr's seurat.RDS rewritten
2026-09-02, after the colleague's objects were built from it). A
symlinked freeze can change silently. Cost of copying is ~12 GB.

## 2026-09-29 — Batch candidates for hazard #1

The core has not confirmed a prep/capture batch, so every technical
grouping that can be reconstructed is tested against timepoint, arm
(SSc only) and group:

- batch_id = sequencing-core run (e.g. 13452-JF). Each is one 10x Flex
  pool of up to 16 probe-barcoded samples = one GEM capture, so it is
  the most plausible capture batch and shares one ambient profile.
- sample_name_date = the YYMMDD prefix of the sample name (9 dates,
  2025-04-22 .. 05-30, each spanning 1-3 pools). Meaning unconfirmed
  (fixation? hybridisation?), so stored under a neutral name.
- cellranger_run.
- site_id (recruiting hospital; biopsy collection and handling) — added
  from the clinical table once it is in data/clinical/. Site cannot
  confound timepoint (every patient has all timepoints) but can
  confound arm, pool, and later improver status.

Rejected: the colleague's `batch` column, which is unique per sample
and therefore not a batch.

Known before running: all 12 healthy controls were captured in their
own pools (12952-JF, 8911-JF, 9260-JF), so SSc-vs-HC differences are
partly technical by construction.

## 2026-09-29 — Ambient RNA: redo SoupX; keep raw counts as the base layer

Decision: raw integer counts stay the canonical `counts` layer, because
a correction can be redone from raw but never undone. SoupX is re-run
here from each library's sample_raw_feature_bc_matrix.h5, output with
`roundToInt = TRUE` (stochastic rounding, preserves expected totals) as
a separate layer. Clustering/annotation use the corrected layer;
pseudobulk DE uses corrected counts with raw counts as a sensitivity
analysis, because ambient fraction may differ by timepoint.

Rejected: (a) round() on wasikowr's non-integer corrected matrix —
deterministic rounding is biased for the many small corrected values,
and her per-sample rho is undocumented; (b) raw counts only, as the
colleague did — reintroduces the contamination SoupX removed.

Open, decide at stage 04 BEFORE looking at outcome: how rho is set.
The colleague tuned 0.1-0.25 by eye per sample. Whatever rule is used
must be written here first, judged on marker leakage (e.g. keratin /
collagen / HBB in immune cells), identical for every library, and never
revisited with improver status in view.

## 2026-09-29 — Patient-level run outputs go to data/patient_level/

Any run output carrying subject IDs (the generated manifest, per-subject
tables) is written by save_patient_table() to
data/patient_level/<freeze>/<scope>/<stage>/<run_name>/, mode 600; the
run's tables/ gets a pointer file. results/ holds only library-level
technical data, counts and summaries.

Why: the data-access rule above says patient-level tables never live in
results/, and the scaffold's 00_1 / 00_2 had been writing
subject-level tables there.

## 2026-09-29 — jobs/run.sh launches with bash -c, not bash -lc

The screen wrapper used `bash -lc`, contradicting AGENTS.md: a login
shell sources .bashrc and can re-activate conda (ENVIRONMENT.md
workaround 1). Changed to `bash -c`; screen inherits PATH from the
launching shell.
