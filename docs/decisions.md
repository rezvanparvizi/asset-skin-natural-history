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

## 2026-10-01 — Publication nomenclature; Subject_ID as the participant key

Values: arm `Placebo` / `Abatacept`, group `SSc` / `HC`, never all caps.
Participant column: `Subject_ID` everywhere, matching the clinical
master table. Applied to the manifest generator (00_0 used toupper()
on the treatment column, which produced PLACEBO / ABATACEPT), cohort
and freeze configs, palettes in R/theme_asset.R, io.R, the pseudobulk
pipeline and the data dictionary. The same convention applies to bulk
skin and bulk PBMC.

Why: labels flow unchanged into plot legends and figure panels; one
spelling at the source avoids per-figure relabelling and mismatched
factor levels across modalities.

Rejected: keeping internal codes (SSC, subject_id) and relabelling at
plot time — every figure script would need a lookup, and joins with
the clinical table would need a rename step.

Escape-therapy censoring now reads `Escape_month` from the clinical
master table (was the placeholder `escape_start_month`).

## 2026-10-01 — mRSS source columns: `0`, `1`, `3`, `6`, `9`, `12`

The clinical master table has two near-identical mRSS series: columns
named `0`…`12` and `MRSSTOTALCALCULATION_0`…`_12`. The owner's decision:
use `0`…`12`. The calculation columns are not used; their derivation is
unknown.

Rejected: `MRSSTOTALCALCULATION_*` — undocumented provenance. A later
check of how often the two series disagree, and by how much, would be
worth doing once; if they differ for some subjects, that is a question
for the DCC.

`mRSS_category` levels: Improver / Worsened / Stable / Set_aside.
Definitions still to be recorded in data_dictionary.md.

## 2026-10-01 — Subject design table readable by AI agents

The owner decided that AI agents may see Subject_ID, treatment arm,
mRSS_category, Ever_escaped, Escape_month, and which timepoints and
modalities each subject has (these are published for this trial). They
are extracted by 00_3_build_subject_design_table.R into metadata/design/
(gitignored, mode 600), outside the data/ tree the agent deny rules
cover. The clinical master table itself, and every other clinical
column (age, sex, site, visit days, serology, lung function), stay
unreadable.

Rejected: lifting the deny rule on data/clinical/ — exposes
quasi-identifiers that no analysis step needs an agent to see; and
pasting values into the chat by hand — slow, and not reproducible.

Open for the owner: confirm the data use agreement permits sending
these fields to the AI service.

mRSS_category calls may be revised by the owner/PI later. Every stage
reads them from the design table at run time, so a revision means
re-running 00_3 and the downstream plots, not editing scripts.

## 2026-10-01 — Subject ID corrections for the scRNA-seq sample map

Two subject IDs in the scRNA-seq sample map are typos, confirmed by the
owner against the clinical master table (found by 00_3, R0007). The
corrections live in metadata/design/subject_id_corrections.csv
(gitignored: sample-map ID, correct Subject_ID, who confirmed, date) and
are applied by 00_0 as the sample map is read, so the manifest and every
downstream table carry the corrected Subject_ID. Any re-run of 00_0, and
any later freeze built from the same sample map, inherits the fix.

Rejected: correcting the IDs downstream (in 00_3 or at each join) —
every consumer of the manifest would have to remember to do it.

To pass on: the colleague's Seurat objects carry the uncorrected IDs in
their patient metadata. Added to the colleague discussion list.

## 2026-10-01 — No escape-therapy censoring of tissue samples

The owner checked the Placebo patients: none started escape therapy
before the month-6 biopsy, and no immunosuppressive treatment was
allowed in Placebo at months 0, 3 or 6. So every Placebo M00/M03/M06
biopsy is treatment-free, and `exclude_post_escape` is false in every
cohort. `placebo` and `placebo_all` now resolve to the same libraries;
both are kept so existing run names stay valid.

Medication review over months 0-6 is being done separately by the owner
with a colleague; escape and treatment are not tracked in this
repository's processing stages.

Rejected: censoring by Escape_month (the scaffold's default) — it would
have dropped every M06 biopsy from the 10 Placebo patients whose escape
month is 6, although those biopsies precede escape therapy. Those
patients are concentrated in the Worsened group, so the censoring
itself would have biased the comparison.

## 2026-10-01 — Cell QC rule (stage 03): per-library adaptive + hard floor

Chosen before any cell-level data are examined, applied identically to
every library, never revisited with mRSS_category in view.

Per library, a cell fails QC if ANY of:
- log10 UMI   < median − 3 MAD (that library)      [low quality / empty]
- log10 genes < median − 3 MAD (that library)
- genes < 200 or UMI < 500                          [hard floor]
- log10 UMI   > median + 5 MAD (that library)      [crude doublet guard
  until a doublet method is chosen; lenient so large cells survive]
- % mitochondrial > median + 3 MAD AND > 10%, ONLY if the Flex probe
  panel contains MT- genes (the script reports whether it does)

Failing cells are kept in the per-cell table with qc_pass = FALSE and a
reason, not dropped, so the rule can be revised without recomputing.
Libraries are flagged for review (not removed) if < 500 cells pass or
> 50% fail. The owner reviews the QC outputs before integration.

Rejected: (A) one fixed threshold for all libraries — depth differs
widely between pools (6 CellRanger outliers, forced-cell HC re-runs),
so a fixed cut makes pool decide which cells survive; (B) MAD alone —
a uniformly poor library has a poor median and keeps its junk.

## 2026-10-01 — freeze01 library exclusions after cell QC (R0013)

Excluded (freeze01.yml excluded_libraries): 13718-JF-14 (M00),
13719-JF-10 (M06), 13596-JF-9 (M03). Cell calling failed (median 3, 7
and 107 UMI per barcode); 0, 0 and 127 cells pass QC. One per timepoint,
so not design-correlated. freeze01 = 218 libraries.

Kept, listed under qc_review_libraries (owner decision): 13481-JF-8,
13719-JF-7 and HC 8911-JF-14, 9260-JF-16. Their passing cells are used
everywhere; any result can be re-checked without them.

Rejected: keeping 13596-JF-9's 127 passing cells — too few to represent
a biopsy, and drawn from a library whose cell calling failed.

## 2026-10-01 — Doublets and ambient genes: follow the lab's routine (jarnagin)

Owner decision: doublet/contamination handling follows the routine in
/home/jarnagin/ASSET_Flex (skin-specific experience of the lab), made
reproducible here:

1. Ambient genes do not drive IMMUNE-compartment clustering. In immune
   (T/NK, myeloid, B/plasma) subcluster embeddings only — as in her
   4_Tcell_Common.R — variable genes matching her ambient set are
   removed before PCA. NOT applied to whole-tissue, fibroblast,
   keratinocyte or vascular embeddings: collagens, DCN, LUM, ACTA2 and
   keratins define those cell types (her whole-tissue and fibroblast
   scripts do not exclude them either). Set — regex
   `^(KRT[0-9AP]|COL[0-9]|MT-|RP[LS][0-9]|MTRNR|HB[ABDGQZ][0-9]?$)` plus
   her explicit list (KRTDAP, SBSN, ..., IGKC, IGHG1, IGHM; copied
   verbatim from 4_Tcell_Common.R into metadata/gene_sets/). Genes stay
   in the data; they are only excluded from HVG selection.
2. Doublet/contamination pockets are called per compartment at
   subclustering, on clusters, with TWO conditions: foreign-lineage
   panel score > 0.40 AND own-lineage identity < 0.30 (her OFF_ABS /
   ON_MIN). Foreign signal alone is ambient, not grounds to drop a cell.
3. Every dropped cluster gets her doublet check: median UMI relative to
   the kept cells; > 1.5x is doublet-like, otherwise it is dropped as
   low-confidence identity, NOT as a doublet. Reported every time.
4. Cells labelled fibroblast/keratinocyte inside the immune compartment
   after subclustering are removed (her 2b step).

Cells are flagged in the label table, not deleted from the counts.
No algorithmic doublet caller in the primary route (scDblFinder is
blocked here and the lab does not use one in this pipeline).

Rejected: scDblFinder / scds as the primary route — not the lab's
routine; may be added later as a sensitivity check only.

## 2026-10-01 — SoupX rule (stage 04), approved by the owner

Per library, from sample_raw_feature_bc_matrix.h5 (that library's own
empty droplets):
- clusters for SoupX: fixed in advance, identical for every library
  (log-normalise, 2,000 HVG — ambient set NOT excluded, since keratin
  and collagen are what separate the clusters SoupX relies on — 30 PCs,
  Louvain resolution 0.8);
- rho from autoEstCont; accepted if 0.01 <= rho <= 0.30, otherwise the
  median rho of the same pool (one capture = one ambient environment);
- output adjustCounts(roundToInt = TRUE): integer corrected counts as a
  separate layer; raw counts remain the base layer.

Skin needs more than one check (owner): every library is judged on
several criteria, reported side by side, before the corrected counts
are used —
  (a) rho per library vs pool, timepoint and arm (must not track design);
  (b) marker leakage before/after: keratin (KRT1/5/10/14) in immune and
      fibroblast cells, collagen (COL1A1/COL1A2/COL3A1) in immune and
      keratinocyte cells, HBB outside erythroid, IGKC/JCHAIN outside
      plasma cells;
  (c) on-lineage markers must NOT drop (COL1A1 in fibroblasts, KRT14 in
      basal keratinocytes, PTPRC in immune) — evidence of over-correction;
  (d) the soup profile's top genes per library, which should be
      keratinocyte/fibroblast-dominated in skin.
A library failing (b)-(c) is reviewed with the owner, not silently
re-tuned. Never judged with mRSS_category in view.

Rejected: per-sample rho tuned by eye (0.1-0.25, upstream object) —
not reproducible and can track timepoint.

## 2026-10-01 — Integration (stage 05): Harmony on batch_id only, with timepoint-preservation checks

Owner requirement: the integration must keep M00 / M03 / M06 differences.
Approved plan:
- Harmony on batch_id (Flex pool) only — pool is independent of
  timepoint (00_2: p = 1.0), so there is no timepoint signal for Harmony
  to remove. Never subject, sample, run date or timepoint.
- Integration changes only the embedding and clusters; composition and
  DE across timepoints are computed from counts + labels with a patient
  random effect, so measured timepoint change cannot be removed by it.
- Checks reported with every integration run: unintegrated vs
  integrated UMAP side by side; timepoint composition per cluster before
  vs after (a timepoint-dominated cluster must not dissolve); mixing
  scores (LISI) for pool (should rise) and timepoint (must not be forced).

Rejected: Harmony on subject (the colleague's usual choice) — it aligns
patients to each other and so removes between-patient differences such
as Improver vs Worsened; on sample — additionally removes within-patient
timepoint change.

## 2026-10-01 — Cell-typing follows the lab routine; compartment stages renamed

Owner decision: adopt the lab's cell-typing routine and gene lists
(docs/lab_routine_celltyping.md; metadata/gene_sets/lab_markers_jarnagin.yml).
Top level uses wasikowr's 16-label vocabulary; compartments follow
jarnagin: fibroblast, immune (all immune incl. mast), lymphoid second
level (T/NK, B/plasma), plus vascular. Stage folders renamed accordingly:
08_sc_subcluster_myeloid -> 08_sc_subcluster_immune,
09_sc_subcluster_tnk -> 09_sc_subcluster_lymphoid,
10_sc_subcluster_adipo_vascular -> 10_sc_subcluster_vascular.
Adipocytes (~0.5% of cells; PLIN1/FABP4/ADIPOQ+) keep a top-level label
but get no compartment.

Deviations from the lab code, each for a stated reason: Harmony on
batch_id not subject; neighbours from PCs not the UMAP; contaminated
clusters flagged rather than deleted; SFPR4 read as SFRP4.

Rejected: a separate myeloid compartment — the lab annotates all immune
cells together and then goes deeper on T/NK and B/plasma; splitting
myeloid out would break comparability with her labels.
