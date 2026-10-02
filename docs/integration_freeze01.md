# Integration, freeze01 — what was done, what was checked, where it stands

Stage 05, Layer 1 (all 1,191,247 QC-pass cells: both arms, all timepoints,
HC). Plan approved 2026-10-01 (decisions.md, "Integration (stage 05)").
Reviewed and extended 2026-10-02. Aggregate numbers only; no IDs.

## 1. What was run

| Run | Script | What | Status |
|---|---|---|---|
| R0022 | `05_1_harmony_batch_first_pass.R` (sha 029d560) | merge SoupX counts, PCA, Harmony, UMAP; object saved | `failed` in runs.csv — only the plotting failed (png needs X11); the object is complete |
| R0023 | `05_2_plot_first_pass_integration.R` | checks + plots | superseded (HC timepoint "" -> NA mixing) |
| R0024 | same | checks + plots | superseded by R0025 |
| R0025 | `05_2_assess_first_pass_integration.R` (renamed 2026-10-02) | R0024 + SSc-only per-lineage and within-patient checks, PCA elbow, panel-label fix | ok (git_dirty: edits not yet committed at run time) |

**Object of record:** `data/objects/freeze01_reference_harmony_first_pass.rds`
(from R0022). The compute part of 05_1 is byte-for-byte the code that ran in
R0022: `git diff 029d560 HEAD` on 05_1 shows only removed plotting and one
added `fwrite` of the unintegrated UMAP, after the compute. No re-run of 05_1
is needed to make the object match the script.

**Pipeline (05_1):** per-library SoupX integer counts (04_1, R0020) merged on
disk (BPCells) -> LogNormalize -> 2,000 HVG (vst, all cells, ambient genes
not excluded) -> ScaleData -> PCA 50 -> PCs by the lab's sd rule
(first PC whose sd is > 99.5% of the previous; floor 15) = **49** ->
`RunHarmony(group.by.vars = "batch_id")`, harmony 2.0.3 defaults (theta 2)
-> UMAP (uwot, 30 neighbours, min_dist 0.3) on the Harmony PCs. No
clustering. `batch_id` = Flex pool (17 pools; 14 carry Placebo libraries).

## 2. Checks and results

Neighbour score = among a cell's 30 nearest neighbours, the fraction sharing
its label, divided by the fraction expected by chance (1 = fully mixed).
Computed in PCA space vs Harmony space on the same cells.

### 2a. Global (R0024/R0025; 100k random cells)

| | PCA | Harmony |
|---|---|---|
| pool | 2.00 | 1.50 |
| timepoint (HC counted as a level) | 1.17 | 1.08 |
| patient | 4.72 | 3.42 |

**Re-assessment:** these numbers mix three things — cell-type composition
(a keratinocyte-rich library has keratinocyte neighbours), the HC pools
(HC-only, so HC separation is pool-driven), and patient nesting in pools.
The drop in "timepoint" from 1.17 to 1.08 was mostly HC being aligned, not
M00/M03/M06 being merged. They are kept for continuity; 2b-2d are the
checks to read.

### 2b. Within patient: are M00 / M03 / M06 still apart? (R0025)

74 SSc patients with >= 2 timepoints (>= 200 cells each; up to 1,000 per
timepoint), kNN among that patient's own cells.

| patients | timepoint score, PCA -> Harmony (median) | per-patient Harmony/PCA ratio (p10 / p50 / p90) |
|---|---|---|
| one pool | 1.35 -> 1.35 | 0.99 / 1.00 / 1.00 |
| multi-pool | 1.44 -> 1.41 | 0.96 / 0.99 / 1.01 |

No patient below a ratio of 0.8. **Within-patient timepoint structure is
untouched** — the owner's requirement for this stage holds. The small drop
in multi-pool patients is expected: there part of the timepoint difference
is pool.

### 2c. SSc only, within each provisional lineage (R0025; up to 30k cells per lineage)

| lineage | pool | timepoint | patient |
|---|---|---|---|
| Keratinocyte | 1.69 -> 1.43 | 1.08 -> 1.07 | 3.76 -> 3.01 |
| Fibroblast | 1.88 -> 1.65 | 1.12 -> 1.10 | 4.55 -> 3.95 |
| Endothelial | 1.43 -> 1.35 | 1.06 -> 1.05 | 2.98 -> 2.77 |
| Immune | 1.72 -> 1.71 | 1.16 -> 1.15 | 3.58 -> 3.51 |
| Melanocyte | 1.81 -> 1.72 | 1.13 -> 1.12 | 4.46 -> 4.10 |
| Gland | 1.54 -> 1.50 | 1.11 -> 1.10 | 3.12 -> 2.98 |
| Mural | 2.08 -> 2.09 | 1.23 -> 1.20 | 4.90 -> 4.54 |
| Plasma (2,539) | 1.62 -> 3.05 | 1.10 -> 1.16 | 1.93 -> 3.47 |

Timepoint is unchanged in every lineage. Patient drops 5-20% in the large
lineages: because pools hold only a few patients each, correcting pool
removes some patient-level variation too (see 2e). The plain pool score
stays well above 1, which looked like weak correction — 2d shows it is
patient nesting, not batch.

### 2d. Pool among neighbours from OTHER patients (scratch check, 2026-10-02)

Same cells as 2c. For each cell, among neighbours from a different patient,
the fraction from the same pool, over the fraction expected (that pool's
share of other patients' cells). This is the clean batch measure, since
same-patient neighbours are automatically same-pool. Not a registered run
(script was in the session scratchpad; aggregate output only) — to be built
into the 05_3 checks.

| lineage | PCA | Harmony |
|---|---|---|
| Keratinocyte | 1.18 | 1.04 |
| Fibroblast | 1.09 | 0.94 |
| Endothelial | 1.05 | 1.00 |
| Immune | 1.06 | 1.06 |
| Melanocyte | 1.12 | 1.09 |
| Gland | 1.00 | 0.99 |
| **Mural** | 1.07 | **1.24** |
| **Plasma** | 1.22 | **2.68** |

- The true pool effect within lineages was modest before integration (1.0-1.18)
  and is at chance after Harmony in the large lineages. Fibroblasts are
  slightly below 1 (mild over-correction).
- **Mural and plasma get worse after Harmony.** Small populations are
  softly assigned to Harmony clusters dominated by large lineages and
  receive those clusters' pool corrections — a known failure mode for rare
  types. Plasma neighbours also become more same-patient (fraction from
  other patients 0.85 -> 0.74). Caveat: "Plasma" is the noisy provisional
  label (2.5k cells, may include Ig-ambient cells).

### 2e. Pool vs mRSS_category (Placebo; design table; 2026-10-02)

Pool taken from the library-ID prefix, 103 Placebo libraries (Set_aside and
excluded libraries removed).
- Per library: Fisher p < 1e-4, Cramér's V 0.53 — **inflated**: a
  patient's libraries mostly share a pool (28 of 38 Placebo patients in one
  pool), so libraries are not independent.
- Per subject, M00 only (one library each): Fisher p = 0.59. No evidence
  that outcome groups were assigned to different pools.
- But several pools are one-sided (one pool 0 Improver / 5 Stable / 6
  Worsened libraries; two pools Improver-only). Because pools ~ small sets
  of patients, Harmony on pool can remove some outcome-associated
  between-patient structure from the embedding. 2c quantifies the patient
  part: 5-20% in large lineages.
- 10 of 38 Placebo patients span > 1 pool: pool must stay a covariate in
  within-subject models (planned since 00_2).

### 2f. Other observations

- **HC:** in the unintegrated UMAP the HC cells form their own streak in the
  keratinocytes; after Harmony it is gone. HC pools are HC-only (00_2,
  V = 1), so Harmony on pool necessarily aligns HC to SSc. This embedding
  cannot show SSc-vs-HC differences; any such comparison must come from
  counts with that confound stated.
- **PCs:** sd falls slowly (PC1 8.26, PC20 2.35, PC30 1.81, PC50 1.34). The
  0.995 rule first triggers at PC48->49, so it chose 49 of 50: it did not
  find an elbow. Acceptable for a heterogeneous top level; noted that the
  rule is not discriminating at this scale.
- **Harmony convergence: unknown.** Run with `verbose = FALSE`; harmony
  2.0.3 stores nothing in the reduction's misc slot. Not evidence of a
  problem; record it next time (verbose log) rather than re-run for it.
- **Provisional `coarse_lineage` is unreliable for immune and gland cells**
  (much of the immune island is labelled Keratinocyte): per-library SoupX
  clusters, argmax of panel means, keratin ambient wins in small libraries.
  Used only for plotting and stratifying checks; replaced at 05_3/06.
- In the Placebo panels, one small dense keratinocyte island appears only
  in Worsened M03 — most likely a single sample. Noted, not interpreted; it
  must not influence the 05_3 resolution.
- The owner has now seen Placebo UMAPs by mRSS_category. The clustering
  resolution must therefore be chosen by a criterion written down before
  05_3 is run (rule 3).

## 3. Assessment

**On the right path for the top level.** The approved design (Harmony on pool
only) does what it should: technical pool structure in the large lineages
goes to chance, and within-patient timepoint change — the Project 1 signal —
is preserved exactly. Nothing found argues for re-running 05_1 or changing
the integration variable.

Its limits, which constrain what this embedding is used for:
1. **Top-level clustering and broad labels only.** Rare and small
   populations (plasma, mural; by extension B cells) are distorted by the
   global Harmony; their subclusters must come from compartment-level
   PCA + Harmony (stages 07-11), with checks 2b and 2d repeated there.
2. **Some between-patient structure is lost from the embedding** (pools are
   small patient groups). Inference does not use the embedding: composition
   and DE come from counts + labels with a patient random effect and pool
   covariate. But a patient-specific state could be merged into neighbours
   at clustering time — 05_3 should flag clusters dominated by one patient
   or one pool, and compare with clusters on the unintegrated PCA.
3. **No SSc-vs-HC reading of the embedding.**

## 4. Recommendations

1. Accept the R0022 object as the Layer-1 top-level embedding (owner's call).
2. 05_3 (not started): neighbours from the 49 Harmony PCs; resolution chosen
   by a marker criterion written in decisions.md first; per cluster report
   timepoint composition (PCA vs Harmony clusters — decision requirement),
   patient and pool dominance, and the across-patient pool score (2d).
3. Compartment stages: re-run PCA + Harmony on batch_id inside each
   compartment; repeat 2b/2d; B/plasma and mural especially.
4. Next time 05_1 runs (freeze02): Harmony verbose (convergence in the log),
   save PCA sd table.
5. Commit today's script changes so the next run is not `git_dirty`.

## 5. Dependency on stage 04 (added 2026-10-02, later)

The SoupX re-assessment (docs/ambient_soupx_freeze01.md) suggests our rho
(median 0.054) under-corrects, possibly by ~2x. If the owner changes the
rho rule, the counts feeding 05_1 change and 05_1 + 05_2 must be re-run.
The integration design and checks stay the same; acceptance of the
embedding should wait for that decision.
