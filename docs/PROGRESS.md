# PROGRESS — where this project stands

**The handoff file.** Read this first in any new session, human or AI. It is
the only file that records *temporal* state: what is done, what is in flight,
what is blocked, and what comes next. `AGENTS.md` says how things are
organised; this says where we are.

**Update it at the end of each working session.** Two minutes. Move finished
items to Done, revise In flight, add anything newly blocked. If it goes stale
it becomes worse than nothing, because it will be believed.

Last updated: **2026-10-03** (evening) · Updated by: Rezvan Parvizi (with Claude Code)

---

## Current focus

**Project 1** — molecular basis of spontaneous mRSS improvement in the ASSET
placebo arm. freeze01 (218 libraries): stage 00, cell QC (03_1) and the
upstream comparison (06_0) are done. Ambient correction re-done with a
fixed rho = 0.13 (04_2, R0026); integration of record R0028 (05_1) checked
by R0034 (05_2). 05_3 top-level clustering ran (R0039); neither rule A
(R0039) nor rule B (06_1, R0042) gave broad labels. **Broad labels are on
hold**: the owner wants wasikowr's own labels for all cells and has asked
both colleagues for them (2026-10-03).
Bulk metadata cleaned (docs/clinical_metadata.md). freeze02 is ON
HOLD (owner decision).

---

## Done

### 2026-10-03 (later) — rule B result; colleague objects surveyed

- **06_1 R0042, rule B (clusters named by wasikowr's labels): no resolution
  qualified.** At res 0.2-0.5 only pericyte vs smooth muscle failed (mixed
  3.0%). kNN cross-check recovers her labels 88-99% (leave-one-out), with
  the same confidence at M03/M06 as at baseline. Rules A and B are set
  aside for now (owner).
- **Her labels exist for all cells only in the June-August version of
  wasikowr's seurat.RDS**, which the 2 Sep rewrite replaced (no `celltype`
  column). jarnagin has no all-cells object; her readable all-timepoint
  objects cover immune cells (104,675, 12 labels) and fibroblasts (143,963,
  Round 1 only, no subtypes). No exported table covers all cells.
  `/hits/home/jarnagin/{ASSET,freezes}` are unreadable to the owner too.
  All of it, plus how she subset immune cells and fibroblasts, is in
  docs/lab_sources.md.

### 2026-10-03 — 05_3 result; 06_0 hardened and re-run

- **05_3 R0039: no resolution met rule A** (owner review). Ambiguous
  clusters hold 29-41% of cells at every resolution (limit 5%); required
  labels covered 7/16 (res 0.05) to 15/16 (0.5, 0.8). B Cells are never
  covered: B and plasma share one cluster (margin 0.01-0.03). Diagnosis
  (aggregate tables only, no outcome): (1) sibling pairs (T/NK, B/plasma,
  Langerhans/myeloid, endothelial/mural at 0.1) give small margins by
  construction; (2) large keratinocyte clusters fail the 1.0 margin
  because keratinocyte clusters are most of the clusters, so z-scoring
  across clusters compresses their own panel's z (e.g. res 0.5 cluster 0,
  212k cells, margin 0.70). Clustering itself looks sound (ARI 0.77 vs
  jarnagin's immune labels at res 0.1). Cluster IDs at all six
  resolutions are saved in R0039 `objects/top_level_clusters.csv.gz`, so a
  revised naming rule can be applied without rebuilding the graph (~4 h).
- 06_0 fixed (commit a1b5c23): fails loudly on a missing SoupX table, an
  unmapped celltype, a non-unique orig.ident map, or cell names that do not
  reduce to a 24-nt barcode. Re-run R0041 reproduces R0021 exactly; R0021
  superseded. 05_3 will also compare with wasikowr's res-0.1 clusters
  (all timepoints) on its next run; R0039 not re-run.
- **04_1 `coarse_lineage` is unreliable for immune cells** (new table
  `wasikowr_vs_our_coarse_lineage`, R0041): her T, Langerhans, mast, B and
  myeloid cells sit mostly in our per-library "Keratinocyte" clusters
  (0.54-0.87). Counts are unaffected (04_2 uses fixed rho 0.13), but do not
  use 04_1 coarse_lineage downstream.

### 2026-10-02 — integration review, SoupX re-decided, colleague findings

- Integration R0022 reviewed (docs/integration_freeze01.md): Harmony on
  pool keeps within-patient timepoint structure exactly; pool effect at
  chance in large lineages (across-patient check); rare types (plasma,
  mural) distorted -> compartment-level re-integration later; HC aligned to
  SSc by construction. Pool vs mRSS_category: no subject-level confounding
  (p = 0.59) but pools are small patient groups. 05_2 extended (R0025).
- SoupX: autoEstCont rho (median 0.054) judged too low; owner chose a fixed
  rho; 0.13 recorded (decisions.md). 04 split: 04_1 estimate (R0020),
  04_2 apply fixed rho. wasikowr uses rho 0.1; her counts are ~30% below
  CellRanger h5 for a non-ambient reason (asked).
- jarnagin's `==` bug quantified: base_control holds ~half the
  Baseline+Control cells; all Baseline_Analysis downstream affected; B-cell
  concentration in few patients is real (data_issues_freeze01.md).
- New docs/lab_sources.md (what is where in colleagues' directories);
  AGENTS.md rules: record colleague findings, no pseudocode.
- Integration of record on rho 0.13 counts: R0028 (30 PCs, Harmony
  converged in 6 iterations) + R0034 checks; R0022/R0025 superseded.
- 05_3 rule A recorded before clustering (decisions.md); top-level marker
  panels approved (metadata/gene_sets/top_level_panels.yml). TPSAB1/TPSB2
  are not in the Flex panel (lab mast panels affected).
- QC by cell type (03_2, R0027) and by timepoint (03_3, R0029): M00
  libraries are shallower; libraries listed in docs/qc_freeze01.md.
- Bulk/clinical metadata (owner's three workbooks, scripts-only access):
  IDs checked against all matrices; skin coding and mRSS_category fixed;
  skin QC with the PBMC metrics (01_1, R0036; all 9 outliers already in
  the owner's Exclude); `Use_downstream` in skin (219) and PBMC (67);
  master `Bulk_baseline_data`; misspelt skin matrix header corrected.
  One-time scripts removed. docs/clinical_metadata.md.

### 2026-10-01 — nomenclature, design table, cell QC

- Nomenclature everywhere: `Placebo`/`Abatacept`, `SSc`/`HC`, `Subject_ID`
  (00_0 had used toupper()).
- Subject design table `metadata/design/subject_design.csv` (00_3, R0016;
  agent-readable by owner decision): Subject_ID, arm, mRSS_category,
  escape fields, scRNA-seq timepoints + library IDs. 88 subjects (44/44),
  79 SSc with scRNA-seq. Placebo mRSS_category: Improver 19, Stable 13,
  Worsened 8, Set_aside 4 (only 1 Set_aside has scRNA-seq).
- Two sample-map ID typos corrected via
  `metadata/design/subject_id_corrections.csv`, applied in 00_0.
- No escape censoring: no Placebo patient escaped before the M06 biopsy.
  Escape/medications are out of scope for processing (owner).
- Cell QC 03_1 (R0017): per-library 3-MAD + floor (>=200 genes, >=500
  UMI) + MT; 1,191,247 of 1,213,176 cells pass (98.2%). freeze01 now
  excludes 3 failed libraries (13718-JF-14, 13719-JF-10, 13596-JF-9);
  4 weak ones kept and listed under `qc_review_libraries`.
- Decisions logged: QC rule, SoupX rule with skin-specific checks, lab
  doublet/ambient routine (jarnagin's pocket screen; ambient-gene HVG
  exclusion in immune embeddings only).
- Reviewed the colleague's pipeline (see Open questions and
  data_issues_freeze01.md).

### Stage 00 — freeze01 (2026-09-29)

- **R0001 `import_raw_libraries`.** 221 libraries (209 SSc from 79
  subjects, 12 HC) in 17 pools. 697 files (12.2 GB) copied into
  `data/raw/`, md5-verified and read-only. Manifest generated.
- **R0002 `build_manifest`.** Structural checks pass. `placebo_all` = 107
  libraries from 39 subjects; `basectrl` = 72.
- **R0003 `batch_vs_design_crosstab`.** Read its `VERDICT.txt`.
  - Pool vs **timepoint: independent** (p = 1.0, V = 0.15). The
    longitudinal design is clean.
  - Pool vs **arm**: p = 0.02, V = 0.34, driven by 13678-DP / 13679-DP.
    Irrelevant within placebo; matters for any arm comparison.
  - Pool vs **group (HC)**: V = 1. The controls are in their own pools
    (CellRanger 9.0.0 and forced cell counts for 5 of them). SSc-vs-HC
    differences are partly technical.
  - 23 of 74 multi-sample subjects (31%) span more than one pool, so pool
    enters within-subject contrasts. Plan: pool as a covariate in
    pseudobulk models. Not yet recorded as a decision.
- Library flags: `13639-JF-11` and `-12` are `timepoint_unresolved` (same
  subject, both labelled M03, no M00). Kept for Layer 1, excluded from all
  cohorts.
- **Sample sheets for building the clinical table** (patient-level):
  `data/patient_level/freeze01/reference/00_manifest_and_batch/import_raw_libraries/`
  `sample_sheet_per_library.csv` and `sample_sheet_per_subject.csv`
- **Discussion list for colleagues:** `docs/data_issues_freeze01.md`
  (controls, pools, 16 unmapped libraries, the timepoint flag, 6 metric
  outliers, the upstream object, code points to pass on).

### Environment (complete, verified)

- R 4.6.1 built from source at `~/R/R-4.6.1` (`--enable-R-shlib`, no X11)
- OpenBLAS 0.3.34 at `~/opt/openblas` (`DYNAMIC_ARCH`, pthread)
- Bioconductor 3.23 — escaped the system R's 3.18
- 292 packages installed; **277 recorded in `renv.lock`**, committed
- CRAN pinned to P3M snapshot `focal/2026-06-01`
- Object-critical pins verified matching the handoff objects:
  Seurat 5.5.0 · SeuratObject 5.4.0 · Matrix 1.7-5 · harmony 2.0.3 ·
  BPCells 0.3.1
- Functional tests passed: Seurat 5 object from an on-disk BPCells matrix with
  PCA; `dream()` fitting `~ timepoint + (1|subject)` on 12 subjects × 2
  timepoints
- Six installation workarounds recorded in `docs/ENVIRONMENT.md`
- `conda config --set auto_activate_base false` — deliberate, see workaround 1

### Repository

- Private GitHub repo `rezvanparvizi/asset-skin-natural-history`, pushed
- Scaffold: three-layer structure, freeze/cohort/labelset axes, provenance
  contract (`init_run()` / `finalize_run()`)
- Analysis stages reordered: metadata → bulk skin → bulk blood → single cell
  (commit `0fb2f17`)
- Bootstrap scripts for server and laptop; `jobs/run.sh` detaching into screen
- `AGENTS.md` + `CLAUDE.md` symlink; `docs/ENVIRONMENT.md`;
  `docs/project_instructions_sections_9-12.md`

### Machines

- U-M MacBook: R 4.6.1 arm64, gfortran 14.2, Xcode CLT licensed, SSH keys to
  GitHub and to `mininubio`, Positron connected via Remote SSH
- Server: SSH key to GitHub, ssh-agent block in `.bashrc` guarded for
  interactive shells only

---

## In flight

- Nothing running. Waiting on wasikowr / jarnagin for the per-cell
  broad labels (see Blocked).

---

## Blocked, waiting on someone else

| What | Waiting on | Why it matters |
|---|---|---|
| **wasikowr's broad labels for all cells** (old seurat.RDS, or a per-cell barcode -> celltype table) | wasikowr / jarnagin (asked 2026-10-03; not before Monday) | Stage 06 broad labels wait for it; rules A/B set aside. |
| **Meaning of pool and sample-name date** | sequencing core / Jarnagin | batch_id is currently the pool (the best available stand-in). Confirm it is the capture batch. |
| **Site ID per patient** | owner, from the clinical table | Adds site to the batch check (00_2 picks up `site_id` automatically). |
| **Which of 13639-JF-11/-12 is Baseline** | Jarnagin / DCC | Both are excluded from cohorts until resolved. |
| **Chemistry / fixation** confirmation | sequencing core | Probe set (v1.1.0), CellRanger (9.0.1; 9.0.0 for 5 controls) and reference are confirmed. Chemistry and fixation are still TODO in freeze01.yml. |
| **Colleague's label tables** (barcode → lineage/celltype/subtype) | Jarnagin | Preferred over inheriting her Seurat objects; decouples us from her ongoing iteration. Questions listed in `docs/inherited_objects.md`. |
| **Which of her objects is current** | Jarnagin | `res0p2` / `res0p3` / `from3b` / `_v2` / `Round1` / `Round2` cannot be resolved from filenames. |
| **`/hits/home/parvizi` ownership** | wasikowr | Currently owned by her, mode 777. Cannot `chmod` it; working in a `700` subdirectory instead, but the parent stays world-writable so the directory could be deleted wholesale. |
| **Lab Dropbox folder convention** | Jarnagin / lab | Needed before writing the sharing side of `scripts/backup.sh`. |
| **Whether `/hits` is IT-backed** | whoever administers it | Determines whether Dropbox is needed as a third copy. |

---

## Next up, in order

0. **When the colleagues reply**: bring in wasikowr's per-cell labels
   (record in inherited_objects.md), compare with R0042, then decide how
   stage 06 assigns broad labels. Item 1 below is on hold until then.
1. **(On hold) Revise the 05_3 naming rule (owner decision, then decisions.md)**:
   R0039 met rule A at no resolution. Decide, on marker criteria only and
   without mRSS_category in view: sibling pairs (T/NK, B/plasma,
   Langerhans/myeloid, pericyte/SMC, endothelial/lymphatic) named at the
   parent level instead of "ambiguous"; B vs plasma left to stage 09; a
   panel score that is not compressed when one lineage dominates the
   cluster count. Apply it to R0039's saved cluster IDs (no graph
   rebuild). Then review timepoint composition and patient/pool-dominated
   clusters at the chosen resolution.
2. **Stage 06, broad labels**: freeze the top-level labels; lab pocket
   screen / doublet routine (flag, not delete). Plan for owner review first.
3. **Compartment stages 07-11**: PCA + Harmony per compartment; repeat the
   within-patient and across-patient pool checks (plasma, mural, B).
   Per-compartment QC loss by timepoint (qc_freeze01.md).
4. **Colleague meeting**: docs/data_issues_freeze01.md (subset bug,
   tryptase genes, wasikowr's counts ~30% below h5, freeze02 samples).
5. **Bulk skin stage 01 analysis**: choose the canonical count matrix
   (proposal: Ensembl raw counts, our own symbol mapping); filter on
   `Use_downstream`.
6. freeze02: 23 newer libraries ALREADY exist (wasikowr, processed Aug
   2026; see data_issues_freeze01.md). 21 are repeats of low-yield
   freeze01 libraries, 2 are unidentified. Waiting on wasikowr's sample
   sheet; then `freeze02.yml` (00_0 must handle the flat CellRanger
   layout) and decide merge_strategy per repeated biopsy.
7. Deferred from before: backfill environment decisions; backup script.

---

## Open questions — decide deliberately, do not default

These are recorded because defaulting on them silently would bias the
analysis.

**Outcome definition.** Global mRSS change versus forearm-site skin score. The
biopsy is forearm; local score tracks local biology better (THBS1 r = 0.76
local versus weaker against global). And binary improver (≥5 points or >20% at
12 months, the definition in papers 01/03) versus continuous trajectory —
binary is comparable to prior work, continuous is better powered with ~19 vs
~19. **Owner's grouping: `mRSS_category` (Improver / Stable / Worsened /
Set_aside) from the clinical master table**; calls may be revised by the
PI. Formal outcome for inference still to be stated before any test.

**Library quality covariate (new 2026-10-02).** M00 libraries are
shallower than M06 (median UMI 1,238 vs 1,401, p = 0.001; a few poor M00
libraries). Decide before stage 14 whether longitudinal pseudobulk models
carry a library-quality covariate (data_issues_freeze01.md).

**Bulk metadata.** Cleaned 2026-10-02; read the workbooks directly and
filter on `Use_downstream` (docs/clinical_metadata.md). Open: autoantibody
and escape columns (owner confirming), canonical skin count matrix.

**Primary timepoint.** Escape is no longer a concern for tissue (no Placebo
escape before the M06 biopsy). The owner wants M00 alone and all
timepoints. **Which is primary for inference: not yet decided.**

**Starting point: wasikowr's object vs raw CellRanger.** Currently raw
(00_0). The colleague kept wasikowr's cell set but re-read raw integer
counts herself (1b_RawCountsReading.R). Owner is weighing lab convention;
a barcode-level comparison of the two cell sets would make the choice
concrete.

**Inherit or rebuild Layer 1.** Taking the colleague's annotation saves months
but inherits her filtering history — note the `prefilter_bcell_2026-08-26`
invalidation in her tree, where an annotation change retroactively invalidated
downstream results. Building fresh gives a clean documented reference layer
both could share. **Worth agreeing a single versioned canonical label set with
her before either of us has a figure.**

**SoupX rho** — decided 2026-10-01 (decisions.md). **Doublets** — decided
2026-10-01: lab routine (pocket screen). Note: no scDblFinder results were
found in the upstream objects, contrary to the 2026-09-29 note.

---

## Deferred TODOs

- macOS Keychain for SSH passphrases (`ssh-add --apple-use-keychain`) — does
  not need an Apple ID; the local login keychain suffices
- Homebrew + HDF5 on the laptop, required before BPCells works locally
- `MuSiC` / `TOAST` — `TOAST` no longer downloads from Bioc 3.23. Needed only
  at stage 14. Options: Bioc archive, vendor from GitHub, Scaden, CIBERSORTx
- Request Apptainer from whoever administers the servers — it is the cleanest
  path to identical environments across the group's 6-7 boxes
- Fill in `config/freezes/freeze01.yml` TODOs once the core confirms them
- Project 2 (QSP model) — separate repository `asset-qsp-skin`, not started,
  deliberately not competing with Project 1

---

## Session log

Newest last. One or two lines per session: what was attempted, what landed.
Detail belongs in `docs/decisions.md` (judgment) and `docs/runs.csv` (runs).

```
2026-09-22  Reviewed the 7 prior ASSET papers; drafted project framing.
            Surveyed the colleague's server directory (1,600-line ls -R):
            hierarchical subclustering pipeline, BPCells on-disk objects,
            MuSiC/Scaden deconvolution, date-stamped output duplication.
2026-09-23  Designed the repository: three layers, freeze/cohort/labelset,
            provenance contract. Scaffold created.
2026-09-24  Laptop setup: R 4.6.1 arm64, gfortran, Xcode license. SSH keys to
            GitHub and mininubio. Repo created and pushed. Surveyed the server
            environment: system R 4.3.2 / Bioc 3.18, no OpenBLAS, gcc 9,
            no containers, no P3M focal binaries. Chose home-built R + renv.
            OpenBLAS 0.3.34 built. R build FAILED — conda ICU 73 vs system 66.
2026-09-25  R 4.6.1 built successfully in a conda-free shell. Bootstrap
            packages installed. BPCells and variancePartition smoke tests
            passed (dream() fits a repeated-measures design).
            conda auto_activate_base set to false.
2026-09-28  renv install across several passes. Workarounds: USE_BUNDLED_LIBUV
            for fs; /sbin on PATH for curl's static-libcurl fallback; presto
            is GitHub-only not CRAN. scDblFinder abandoned (scrapper needs
            C++20; conda-forge gcc route also needs gfortran_linux-64).
            renv.lock initially recorded 1 package — snapshot.type had to be
            "all", not "explicit". 277 packages now locked and committed.
            Analysis stages reordered to metadata -> bulk -> single cell.
            docs/ENVIRONMENT.md rewritten with all six workarounds.
2026-09-29  Surveyed colleague's ASSET_Flex tree (read-only). Stage 00
            built: raw copied + md5 (R0001), manifest (R0002), batch
            check (R0003: timepoint clean, HC fully confounded with
            pool). Patient-level outputs moved to data/patient_level/.
            Discussion list for colleagues written.
2026-10-01  Nomenclature fix; design table (00_3) with 2 ID corrections;
            escape censoring dropped; cell QC 03_1 (98.2% pass, 3 failed
            libraries excluded). Reviewed jarnagin/wasikowr pipeline;
            adopted lab doublet routine; SoupX rule logged. SoupX run and
            accepted (rho median 0.054). Lab cell-typing routine and gene
            lists copied; stages renumbered (08 immune, 09 lymphoid,
            10 myeloid, 11 vascular, 12-18). Upstream match: 98.3% of
            wasikowr's cells pass our QC; her broad labels check out.
            Found 23 newer libraries in her object (freeze02, on hold).
            Run-id race fixed in provenance.R. Unpushed history scrubbed
            of Subject_IDs/specimen codes before push.
2026-10-02  Integration R0022 reviewed and checks extended (R0025).
            jarnagin's subset bug quantified. SoupX re-decided: fixed rho
            0.13 (04_2); chain 04_2 -> 05_1 -> 05_2 launched. wasikowr's
            counts ~30% below h5 (not ambient). docs/lab_sources.md.
            Integration of record R0028/R0034 (rho 0.13). Rule A + panels
            for 05_3; 05_3 running. QC by timepoint (qc_freeze01.md).
            Bulk/clinical workbooks cleaned via scripts; Use_downstream
            in skin (219) and PBMC (67); one-time scripts removed.
2026-10-03  05_3 R0039 finished: no resolution met rule A (ambiguous
            29-41%; B never separates from plasma). 06_0 hardened and
            re-run (R0041 = R0021; adds coarse-lineage table, which shows
            04_1 coarse_lineage misplaces immune cells). R0021 superseded.
2026-10-03  06_1 rule B (R0042): no resolution qualified (pericyte vs SMC).
            Surveyed all jarnagin objects/exports: no all-cells labels
            readable; labels requested from colleagues. lab_sources.md
            records her immune/fibroblast subsetting steps and outcomes.
```
