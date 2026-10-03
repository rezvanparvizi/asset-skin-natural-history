# Lab sources: what is where in colleagues' directories

**Read this before opening anything in a colleague's directory.** It
records what has already been found there: paths, what each object or
file contains, and what was learned from it, with the date it was checked.
Add to it in the same session as the check; if a check produces
barcode- or patient-level output, ask the owner whether and where to save
it (AGENTS.md, "If you are an AI agent").

Everything here was read in place, read-only; nothing was copied unless
stated. Copies of colleague files get an entry in `inherited_objects.md`
instead. Discussion points and questions for colleagues live in
`data_issues_freeze01.md`; the lab's cell-typing routine in
`lab_routine_celltyping.md`. No Subject_IDs or specimen codes in this file.

---

## wasikowr — upstream processing (CellRanger, QC, SoupX, top-level labels)

Root: `/hits/home/wasikowr/ffpe/`

| Path (under root) | What it is | Checked |
|---|---|---|
| `Asset/seurat.RDS` | Main object, rewritten 2026-09-02, 1.54 GB. 1,302,390 cells, 18,123 genes, **242 samples** = our 219 + 23 newer libraries. Counts layer non-integer (SoupX-corrected). Labels: `RNA_snn_res.0.1` / `seurat_clusters` only — no cell-type column. | 2026-10-01, 10-02 |
| `Asset/seurat_amb.RDS` | 65 MB, 2026-08-31. Same 1,302,390 cells, 18,123 genes; metadata `orig.ident, nCount_RNA, nFeature_RNA, batch, orig.nb`; counts identical to seurat.RDS (non-integer, same totals) — not an uncorrected copy. No rho stored. | 2026-10-02 |
| `Asset/firstrun/demo.txt` | Sample map; identical (md5) to jarnagin's `results/demo.txt`; the source of our freeze01 sample map. | 2026-10-01 |
| `Asset/firstrun/input.txt` | Our 209 SSc libraries. | 2026-10-01 |
| `Asset/firstrun/asset_batch5.csv` | The 23 newer libraries: `sampleID, path, batch` (= original sample number), `condition`. | 2026-10-01 |
| `Asset/bpcells_obj1a`, `2a`, `3a`, `4a`, `5a` | **BPCells count matrices (genes x cells) behind `seurat.RDS`**: 297,717 + 482,819 + 212,626 + 217,305 + 91,923 (5a, Aug 31 = the 23 newer libraries) = 1,302,390 cells, exactly seurat.RDS. 18,058-18,123 genes. **Counts only, no labels.** jarnagin's all-sample objects point at 1a-4a. `bpcells_obj1`..`4` (no "a", Apr 30) are older versions of 296,225 / 479,876 / 211,960 / 215,977 cells. | 2026-10-03 (dimensions only) |
| `Asset/firstrun/all cells ... .csv` (36 files, May 5) | Her DE results (autoantibody groups at M3+M6 vs negative, Baseline vs M3/M6, SSc vs Control), plain and `_filtered`. Not matrices, not labels. | 2026-10-03 (listing only) |
| `Asset/Batch1`..`Batch5`, `Batch*a` | Per-batch working folders. `Batch1` holds `ffpe_report.pdf` (95 MB), `harmony.RData` (47 GB), `markers.csv`. Reports **not read** (no PDF tools on the server) — may document her settings. | 2026-10-02 (listing only) |
| `Asset/Batch{1,1a,2a,3,3a,4,4a,5}/markers.csv` | Seurat `FindAllMarkers` output per batch (`p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene`; Batch1 12,757 rows). Cluster NUMBERS only — no cluster -> cell-type map, so it does not give her 16-label marker list. | 2026-10-02 |
| `{14794,15195,15231,15242}-JF_v1/` | **CellRanger outputs for the 23 newer libraries (freeze02 material)**, re-run by wasikowr Aug 2026. Flat `per_sample_outs` layout (no `count/` subfolder — 00_0 must handle it). All forced cells; probe set v1.1.0; probe barcodes of a different series than freeze01. `15231-JF_v1` is a combined pool with 15230-JF. A v2-probe run of 15242-JF was superseded by `15242-JF_v1`. | 2026-10-01 |
| other `ffpe/*-JF*` folders | Other projects' runs (controls, keloids, scalp, ...). Not ASSET; not examined. | 2026-10-02 (listing only) |

Findings:
- **No SoupX script** anywhere under `Asset/`. Owner asked her (2026-10-02):
  **her pipeline uses SoupX rho = 0.1.**
- Her counts hold ~30% fewer UMIs than the CellRanger h5 for the same
  cells, near-uniform across genes (soup genes lose less than average) —
  **not** the ambient correction; source unknown, asked
  (`ambient_soupx_freeze01.md`).
- 98.3% of her cells pass our QC; her broad labels (16 types, via
  jarnagin's base_control object) check out against markers (06_0, R0021).
- Her clusters were computed on the UMAP, not PCs (`lab_routine_celltyping.md`).
- **Where her 16 broad labels for ALL cells went (2026-10-03).** jarnagin's
  full-dataset scripts (`FullDataset_Analysis/1a`, `2a`, `3a`) all
  `readRDS("/hits/home/wasikowr/ffpe/Asset/seurat.RDS")` and subset it on
  `celltype`, so the June-August version of seurat.RDS carried `celltype`
  for every cell, all timepoints. The 2 Sep rewrite dropped that column. No
  copy of the old version, and no exported per-cell label table for all
  cells, exists anywhere readable (jarnagin's `ASSET_Flex` tree, the
  `bpcells_obj*` matrices, `firstrun/`). Readable labels at all timepoints
  exist only for immune cells and fibroblasts (jarnagin section below).
  Asked both colleagues for the old object or a per-cell table (owner,
  2026-10-03); waiting.
- freeze02: 21 of the 23 newer libraries repeat freeze01 libraries (same
  biopsy; matched by sample number / specimen code); #32 and #68 have no
  original in the sample map (inferred repeats of 13454-JF-7 and
  13484-JF-11, unconfirmed). On hold (decisions.md 2026-10-01). Detail:
  `data_issues_freeze01.md`, "23 newer libraries" and "Cross-check".

## jarnagin — compartment subclustering, labels, downstream analysis

Scripts: `/home/jarnagin/ASSET_Flex/A_expr_scRNA_ASSET/`
(`Baseline_Analysis/`, `FullDataset_Analysis/`, `BulkRNAseq_Analysis/`,
`TCAT/`, `4_Tcell_Common.R`; not under git).
Objects: `/home/jarnagin/ASSET_Flex/objects/`. Also `results/`, `reference/`.
Not readable: `/hits/home/jarnagin/ASSET` (mode 770) and
`/hits/home/jarnagin/freezes` (mode 750); the owner cannot open them
either (2026-10-03). `A_expr_scRNA_ASSET/objects/` is a set of symlinks
into `freezes/ASSET_Flex_freeze_v2/objects/` with the same file names as
`ASSET_Flex/objects/` plus `ASSET_sim.h5ad`; whether those targets are
newer versions is unknown.

### All objects in `ASSET_Flex/objects/` (inspected 2026-10-03)

Every Seurat object was loaded read-only and summarised (aggregate output
only; log `logs/20261003_162633_inspect_jarnagin_objects_freeze01_reference.log`,
script in the session scratchpad, not kept). **jarnagin has no all-cells
object of her own**: every object is either Baseline + HC only, or one
compartment. Counts are wasikowr's SoupX-corrected, non-integer counts
unless marked "integer" (raw CellRanger h5, her `1b`). Integration in
every compartment object is Harmony on `subject` (= one value per library)
or `run_date` (T cells), never pool.

| Object | Cells | Libraries / timepoints | Counts | Labels | Notes |
|---|---|---|---|---|---|
| `ASSET_base_control.RDS` (+ `base_control_BPCells/`) | 181,766 | 83; Baseline + HC | corrected | `celltype` (16) | halved by the `==` bug |
| `ASSET_base_control_integer_BPCells.RDS` | 181,766 | same | integer | same | |
| `ASSET_base_control_subclustered_integers_BPCells.RDS`, `..._subclustered_v2_BPCells.RDS` | 180,459 | same | integer | `celltype` | Jun 15; immune cleaned |
| `ASSET_base_control_Subtyped_res0p3_BPCells.RDS`, `..._Subtyped_LineageCalls.RDS` | 179,358 | same | integer | + fibroblast subtypes, lineage calls | Aug 18 / 21 |
| `ASSET_base_control_BcellPlasma_Clean.RDS` | 179,282 | same | integer | + B/plasma | Aug 26 |
| `ASSET_MuSiC_reference_sce.rds` | 82,369 | Baseline + HC | — | 18 types incl. plasma, pDC | subsample for MuSiC/Scaden |
| `ASSET_base_control_immunecells.RDS` | 15,979 | Baseline + HC | corrected | 4 immune types | Jun 11 |
| `ASSET_basecontrol_Fibs_*` (6 objects) | 21,944 unfiltered; 20,843 filtered | Baseline + HC | corrected | fibroblast subtypes (res 0.2 `from3b` / 0.3, Steele-aligned) | |
| `ASSET_baseline_BcellPlasma_Subclustered.rds` | 649 | 55; Baseline | integer | B cells | |
| **`ASSET_AllSamples_Samples_immunecells.RDS`** | 118,087 | **219; M00/M03/M06/HC** | corrected | `celltype`: Myeloid 71,709, T 28,603, B 8,984, Langerhans 8,791 | Jun 23, earlier all-sample attempt; re-embedded (harmony 19 dims) |
| **`ASSET_Round1_imm.RDS`** | 107,818 | 219; all | corrected | `celltype` (4 immune), `temp_subclusters` | Jun 25 |
| **`ASSET_Round2_immFilt.RDS`** | 104,675 | 219; all | corrected | `celltype`, **`subclusters` (12 immune labels)** | Jun 25; source of all later immune work |
| **`ASSET_Round1_fibs.RDS`** | 143,963 | 219; all | corrected | `celltype` = Fibroblasts | Jun 29; no Round 2 |
| `ASSET_FullDataset_Tcells_Triaged.rds` / `_Subtyping_Subclustered.rds` / `_Subtyped.rds` | 25,322 | 219; all | corrected | T/NK subtypes (Oct 2 version) | |
| `ASSET_FullDataset_CD4CD8_Reclustered.rds` | 18,186 | 217 | corrected | | Aug 21, superseded |
| `ASSET_FullDataset_CD4_Reclustered.rds` / `_CD4_Labelled.rds` | 15,660 / 15,730 | 218 | corrected | CD4 subtypes; TCAT assay | Sep 23 / Oct 2 |
| `ASSET_FullDataset_CD8_Reclustered.rds` / `_CD8_Labelled.rds` | 6,069 | 207 | corrected | CD8 subtypes; TCAT assay | Oct 2 |
| `ASSET_FullDataset_BcellPlasma_Subclustered.rds` | 7,058 | 189; all | corrected | B / plasma | Aug 25 |
| `ASSET_FullDataset_BcellPlasma_Clean.rds` | 6,582 | 181; all | corrected | B / plasma clean | Aug 25 |
| `*_SingleR_Monaco*.rds`, `*_ClusterDE_*.rds`, `ASSET_bulk_*.rds` | — | — | — | SingleR score lists, DE tables, bulk matrices | not cell objects |
| `ASSET_base_control_cohort_barcodes.txt` | — | — | — | barcodes of the base_control cohort (used by her 4a) | |

Counts of the all-sample objects (`AllSamples_Samples_immunecells`,
`Round1_imm`, `Round2_immFilt`, `Round1_fibs`) live in wasikowr's
`/hits/home/wasikowr/ffpe/Asset/bpcells_obj1a..4a`; the labels are in the
`.RDS` metadata and survive without them.

### Exported tables and matrices (not R objects; checked 2026-10-03)

None covers all cells with `celltype`.

| File (under `ASSET_Flex/`) | Cells | Scope | Content |
|---|---|---|---|
| `results/BulkRNAseq_Analysis/Scaden/data/ref_counts.mtx`, `ref_barcodes.tsv`, `ref_celltypes.tsv` | 82,369 | Baseline + HC (= MuSiC reference) | integer counts, 17,959 genes; 18 cell types |
| `results/Baseline_Analysis/2026-08-21_ASSET_basecontrol_BaselineLineageCalls.csv` | — | Baseline + HC | fine lineage calls |
| `results/FullDataset_Analysis/TCAT/2026-08-20_ASSET_fulldataset_TNK_forTCAT/` (`counts.mtx`, `barcodes.tsv`, `metadata.csv`) | ~25k | all timepoints, T/NK | counts + metadata |
| `results/FullDataset_Analysis/2026-08-25_ASSET_fulldataset_BcellPlasma_CleanLabels.csv` | 6,582 | all timepoints, B/plasma | per-cell labels (her "transfer table") |
| `objects/Steele_LabelTransfer/..._query_mtx/` (`counts.mtx`, `metadata.csv`) | 20,843 | Baseline fibroblasts | counts + fibroblast subtypes |
| `results/BulkRNAseq_Analysis/Scaden/data/*.h5ad` | — | — | Scaden inputs and simulated bulk |

### How she subset immune cells and fibroblasts on the full dataset (read 2026-10-03)

Scripts: `FullDataset_Analysis/2a_ImmuneCell_Subclustering.R` (immune),
`3a_Fibroblast_Subclustering.R` (fibroblasts), then `4a-4c` (T/NK) and
`6a-6b` (B/plasma) from the immune object. Outcome numbers come from her
saved objects and her `results/FullDataset_Analysis/*ClusterSummary.csv`,
`*ContaminationScreen.csv`, `*SubtypeCounts.csv` and `*CleanComposition.csv`.

**Shared functions (2a and 3a are identical up to the marker lists).**
- `reembed()`: NormalizeData -> 2,000 HVG -> ScaleData -> PCA 50 -> PCs up
  to where sd stops falling (next sd > 0.995 x previous) ->
  `RunHarmony(group.by.vars = "subject")` -> UMAP on the Harmony dims.
  `subject` is one value per library, so within-patient timepoint
  differences are removed.
- `cluster_and_clean()`: Louvain on the Harmony dims at resolution 0.2;
  drop every cluster where > 30% of cells express any listed foreign marker
  (log-normalised > 0.5); re-embed; repeat up to 5 times or until clean.
- Per round, for review: FindAllMarkers, SingleR (Human Primary Cell
  Atlas, main labels), CellCycleScoring, marker dot plots and heatmaps,
  and a cluster summary CSV. **Labels are then typed in by hand** as a
  cluster -> label vector in the script.

**Immune (2a), all 219 libraries, all timepoints.**
1. Subset the old seurat.RDS on `celltype %in% c("Myeloid Cells", "T Cells",
   "B Cells", "Langerhans Cells")`. Mast cells are **not** included. The
   input count was not saved. The earlier all-sample attempt
   (`AllSamples_Samples_immunecells.RDS`, Jun 23) holds 118,087 cells with
   these four labels, plausibly the same subset, but this is not confirmed.
2. Round 1: re-embed, then the clean loop with fibro (COL1A1, FN1) and kerat
   (KRT1/2/10/14/15/17) markers -> **107,818 cells, 18 clusters**
   (`ASSET_Round1_imm.RDS`). Hand labels (`temp_subclusters`) include one
   "Keratinocytes" cluster (cluster 8, 3,143 cells; top genes KRT2, DSP,
   KRT1).
3. Round 2: drop the "Keratinocytes"/"Fibroblasts"-labelled clusters
   (3,143 cells), re-embed, clean loop with fibro markers only (removed
   nothing more) -> **104,675 cells, 17 clusters** -> hand labels in
   `subclusters` (`ASSET_Round2_immFilt.RDS`):

   | subclusters | cells | | subclusters | cells |
   |---|---|---|---|---|
   | Macrophages (5 clusters) | 36,316 | | cDC1 | 3,881 |
   | T Cells | 20,356 | | NK Cells | 3,181 |
   | cDC2B (2 clusters) | 16,672 | | Plasma Cells | 3,114 |
   | Langerhans Cells | 7,435 | | cDC2A | 1,377 |
   | Proliferating | 6,324 | | pDCs | 1,363 |
   | B Cells | 3,944 | | Neutrophils | 712 |

   Two hand labels to question: Round 2 cluster 9 (2,336 cells; VSTM1,
   S100A12, VCAN, LILRA5, MCEMP1; SingleR Monocyte) is labelled cDC2B, and
   cluster 14 (809; NR4A1, CXCL8, FOSB, ATF3; SingleR Monocyte) is labelled
   Macrophages. Both read as monocyte-like. Her broad `celltype` "B Cells"
   (8,113 here) is a B + plasma + pDC mixture; `subclusters` separates them.
4. T/NK (4a-4c): subset `subclusters %in% c("T Cells", "NK Cells",
   "Proliferating")` = 29,861 cells; re-embed with ambient genes removed
   from the HVGs, Harmony on `run_date` (`4_Tcell_Common.R`), resolution
   0.4; lineage triage by absolute detection rates moved 4 clusters out
   (4,539 cells: proliferating myeloid 1,457 + 781, myeloid 1,362,
   proliferating DC 939) -> **25,322 T/NK cells**. 4b: cell cycle
   regressed, resolution 0.6, 14 clusters, four callers (Poon panels,
   canonical markers, SingleR Monaco, TCAT); 4c: hand labels, NK split on
   2026-09-17. Labels on 2026-10-02: CD4 T cells c0-c5 (13,397), CD8 T
   cells c0-c1 (4,762), Tregs 2,263, NK cells 2,085, Cytotoxic NK cells
   696, Gamma-delta T cells 648 (her own caveat: could be ILC),
   Proliferating T cells 1,471. CD4 and CD8 were then re-clustered
   separately (`Tcell_Lineage_Subtyping/CD4`, `CD8`; Oct 2 objects).
5. B/plasma (6a-6b): B Cells + Plasma Cells from `subclusters` = 7,058;
   Harmony on `subject`, resolution 0.3. 6b drops cells only when foreign
   signal is present AND own identity is lost (pocket screen; median UMI
   showed the flagged clusters are not doublets) -> **6,582 cells**
   (476 dropped). 3 small clusters are B/plasma-mixed ("bridge
   candidates", 45-132 cells each); pseudotime was run on Aug 27.
6. Not covered at all timepoints: mast cells (left out at step 1).

**Fibroblasts (3a), all 219 libraries, all timepoints.**
1. Subset the old seurat.RDS on `celltype == "Fibroblasts"`. The input
   count was not saved.
2. Round 1: re-embed, clean loop with nerve (NRXN1, CDH19) and kerat
   (KRT1/2/10/14/15/17) markers -> **143,963 cells, 14 clusters**
   (`ASSET_Round1_fibs.RDS`, Jun 29). Marker panels for subtyping: SFRP2,
   CCL19, CLDN1, COL8A1 (with the `SFPR4` typo), FMO1, FMO2, TNN.
3. **Round 2 is entirely commented out** (it is still the immune template,
   unedited). So the full-dataset fibroblasts have **no hand labels, no
   subtypes, and no second cleaning pass**. Fibroblast subtypes exist only
   for Baseline + HC (`Baseline_Analysis/3a-3d`, 20,843 cells, 7 subtypes,
   plus Steele scVI label transfer), which is downstream of the `==` bug.
4. Residual non-fibroblast clusters remain in Round 1 (her cluster
   summary): cluster 5 (5,304 cells; CCL14, PLVAP, VWF, PECAM1 =
   endothelial) and cluster 10 (2,574; KRT2, DMKN, ABCA12 = keratinocyte);
   clusters 1, 7, 9 are called smooth muscle by SingleR. Endothelial and
   mural markers were never in the fibroblast contamination list.

**Observations to keep in mind when reusing these labels** (my reading of
her tables, not her statements):
- The `cycling` column is not usable as biology: CellCycleScoring calls
  25-50% of immune cells and 54-70% of every fibroblast cluster "cycling"
  (S or G2M), whereas true proliferating clusters are at ~100%. Her T-cell
  scripts use it as a nuisance covariate only.
- `2a_ImmuneCell_Subclustering.R:452` ends in a stray `ƒ`
  (`imm_filt$subclusters <- Idents(imm_filt)ƒ`), so the script no longer
  runs past that line as written. The saved Jun 25 object does carry
  `subclusters`, so the character was probably added after the run.
- Harmony on `subject` / `run_date` means none of her embeddings or
  clusters can be reused for within-patient timepoint questions; her
  **labels** can still be transferred by barcode.

Findings:
- **Recycling-subset bug**, `Baseline_Analysis/1a_Baseline.R:90`
  (`time == c("Baseline","Control")`): base_control holds ~47-49% of the
  Baseline+Control cells (~200,600 lost). Every `Baseline_Analysis/` script,
  the 4a T-cell cohort, TCAT `basecontrol`, and the MuSiC/Scaden reference
  descend from it; `FullDataset_Analysis/` does not. B cells 1,648 -> 836;
  11 samples pushed under her "> 3 B cells" cutoff; B-cell concentration in
  few patients is real (top-5 samples hold 57% either way). Detail:
  `data_issues_freeze01.md`, "Recycling-subset bug".
- She did not run SoupX: she uses wasikowr's corrected counts (immune, 4a)
  or raw h5 counts (`1b`, TCAT). Her 4a comment documents raw-count ambient
  in T cells (KRT14 in 46% of T cells raw vs 0% corrected).
- `3a_Fibroblast_Subclustering.R`: `SFPR4` typo (SFRP4 silently dropped).
- `4_Tcell_Common.R:88,113`: mast panels use TPSAB1 / TPSB2, which are not in
  the Flex probe panel (score nothing; checked 2026-10-02).
- Harmony on `subject` (= per-biopsy specimen code) in compartment scripts.
- Her metadata `patient` = our Subject_ID + "B"; carries the same two
  sample-map typos we corrected.
- Which object is current per compartment (`res0p2`/`res0p3`/`from3b`/
  `_v2`/`Round1`/`Round2`) is still unanswered (`inherited_objects.md`).

## What we derived from these sources

| Output | Where | Run |
|---|---|---|
| Barcode -> label table: our cell IDs joined to labels from the 7 objects above (wasikowr `seurat.RDS` + six jarnagin objects) | `metadata/inherited_labels/freeze01_upstream_labels.csv.gz` (barcode-level, gitignored; agents do not open it) | R0021 (06_0), re-run as R0041 |
| Match and coverage summaries (aggregate) | `results/freeze01/reference/06_sc_lineage/match_upstream_cells_and_labels/tables/` — `upstream_sources.csv`, `cellset_overlap_per_library.csv`, `inherited_label_coverage.csv`, `inherited_celltype_by_timepoint.csv`, `wasikowr_vs_jarnagin_labels.csv`, ... | R0021 |
| Bug impact, wasikowr counts vs raw, per-gene check | numbers in `data_issues_freeze01.md` and `ambient_soupx_freeze01.md`; scratch scripts, aggregate output only, not kept | 2026-10-02 |
