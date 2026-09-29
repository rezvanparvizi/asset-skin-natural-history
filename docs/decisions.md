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


