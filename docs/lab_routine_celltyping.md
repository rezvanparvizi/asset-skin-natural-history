# The lab's cell-typing routine (wasikowr → jarnagin), and how we adopt it

Reviewed 2026-10-01 from `/hits/home/wasikowr/ffpe/Asset/` and
`/home/jarnagin/ASSET_Flex/A_expr_scRNA_ASSET/` (read-only). Gene lists are
copied verbatim into `metadata/gene_sets/lab_markers_jarnagin.yml` and
`metadata/gene_sets/ambient_genes_lab.yml`. Owner decision (2026-10-01):
follow this routine.

## Who did what

| Level | Done by | Where |
|---|---|---|
| Cell calling, >200-gene filter, SoupX, Harmony, UMAP, **top-level cell type** | wasikowr | `seurat.RDS` (no script in the folder) |
| Compartment subclustering, contamination cleaning, subtypes | jarnagin | `FullDataset_Analysis/2a, 3a, 4a-c, 6a-b` |
| Raw integer counts re-read for wasikowr's cells | jarnagin | `Baseline_Analysis/1b_RawCountsReading.R` |

## 1. Top level — wasikowr's 16 cell types

jarnagin's `1a_FullCellTyping.R` only plots the `celltype` column already in
wasikowr's object; every compartment script subsets on it. Labels and
counts in the baseline + control object (`ASSET_base_control.RDS`, n = 181,766):

| celltype | cells | | celltype | cells |
|---|---|---|---|---|
| Keratinocytes | 93,232 | | Melanocytes | 3,397 |
| Fibroblasts | 21,944 | | Pericytes | 2,560 |
| Follicle Cells | 20,859 | | Sebocytes | 2,282 |
| Myeloid Cells | 9,711 | | L Endothelial Cells | 1,617 |
| Endothelial Cells | 8,736 | | Smooth Muscle Cells | 1,568 |
| Eccrine Cells | 7,397 | | Langerhans Cells | 1,507 |
| T Cells | 3,831 | | Mast Cells | 1,235 |
| | | | Adipocytes | 960 |
| | | | B Cells | 930 |

- **How these were assigned is not documented.** Each batch folder has a
  `markers.csv` (FindAllMarkers per cluster) and wasikowr made marker PDFs,
  but the marker list and the cluster → label map are not there. Clusters in
  the current object come from neighbours on the 2-D UMAP at resolution 0.1.
  **Ask wasikowr for her lineage marker list.**
- **Adipocytes are real but rare:** ~0.5% of cells; of the 960, PLIN1 90%,
  FABP4 91%, ADIPOQ 71% detected. 3 mm biopsies occasionally reach fat.
- The `remove` column (32% "remove") is **not a QC flag**: it is exactly
  "all non-epithelial cell types" (keratinocyte, follicle, eccrine,
  sebocyte = keep). jarnagin drops the column, not the cells.

## 2. Compartments — jarnagin's routine

Same pattern in every compartment (`reembed()` + `cluster_and_clean()`):

1. Subset by wasikowr's `celltype`.
   - Immune: Myeloid + T + B + Langerhans (Mast is NOT included).
   - Fibroblast: Fibroblasts.
2. Re-embed: LogNormalize, 2,000 HVG, PCA, PCs chosen where sd stops
   dropping (> 0.995 x previous), Harmony on `subject`, UMAP.
   Immune T-cell passes additionally exclude ambient genes from the HVGs
   and (for subtyping) regress cell cycle; Harmony there is on `run_date`.
3. **Contamination loop** at resolution 0.2: a cluster is dropped if
   > 30% of its cells express any foreign marker (log-normalised > 0.5);
   re-embed and re-cluster; up to 5 iterations.
   - Immune: fibro (COL1A1, FN1), kerat (KRT1/2/10/14/15/17).
   - Fibroblast: nerve (NRXN1, CDH19), kerat (same six keratins).
4. Annotate with marker dot plots and heatmaps (lists in the YAML),
   SingleR against the Human Primary Cell Atlas, cell-cycle scoring, and
   FindAllMarkers; labels are assigned **by hand** per cluster.
5. Round 2: drop residual keratinocyte/fibroblast clusters inside the
   immune compartment, repeat 2-4.
6. Second level inside immune:
   - **T/NK** (`4a-4c`): lineage triage by detection-rate evidence panels,
     SingleR Monaco, CD4/CD8 re-clustering, TCAT programs.
   - **B/plasma** (`6a-6b`): two-criteria pocket screen (foreign > 0.40
     AND own identity < 0.30) + doublet check (median UMI > 1.5x).
7. Fibroblast subtypes (`3a-3c`): 7 marker-defined groups (SFRP2, CCL19,
   CLDN1, COL8A1, FMO1, FMO2, TNN) plus scVI label transfer from the
   Steele skin atlas for baseline fibroblasts.

Final immune labels (2a round 2): Macrophages, T Cells, NK Cells,
B Cells, Plasma Cells, cDC1, cDC2A, cDC2B, pDCs, Langerhans Cells,
Neutrophils, Proliferating.

## 3. How this repository adopts it

| Step | Adopted as-is | Deviation, and why |
|---|---|---|
| Top-level lineages | wasikowr's 16 labels as the vocabulary | Assigned here by marker panels (lineage + evidence panels), because her assignment rule is not documented. Compared to her labels by barcode where cells overlap. |
| Compartments | fibroblast; immune (all); lymphoid second level; vascular | Mast cells go into the immune compartment (her subset leaves them out). Vascular (endothelial, lymphatic, pericyte, smooth muscle) has no lab script yet. Keratinocytes/epidermal appendages and adipocytes keep top-level labels only. |
| Contamination loop | thresholds 0.5 / 30% / 5 iterations, her marker lists | Clusters are flagged in the label table, not deleted from the counts. |
| Pocket screen + doublet check | as-is (0.40 / 0.30 / 1.5x) | — |
| Ambient genes out of HVGs | immune embeddings only, her list | — |
| Integration | — | Harmony on `batch_id` (pool), never `subject`/`run_date`: integrating on subject removes between-patient differences (Improver vs Worsened). decisions.md 2026-10-01. |
| Clustering | resolution 0.2 inside compartments | Neighbours from the Harmony PCs, never from the UMAP. Resolution recorded before mRSS_category is looked at. |
| Counts | — | SoupX re-run with a documented rule; raw counts kept alongside. |

Stage folders: `07_sc_subcluster_fibroblast`, `08_sc_subcluster_immune`
(all immune cells, first pass), then two subsets of it:
`09_sc_subcluster_lymphoid` (T/NK and B/plasma) and
`10_sc_subcluster_myeloid` (macrophages, DCs, Langerhans, mast);
`11_sc_subcluster_vascular`.

## 4. Issues found in the lab code (to pass on)

- `Baseline_Analysis/1a_Baseline.R:90` — `time == c("Baseline", "Control")`
  recycles and keeps about half the cells of each group (her baseline object
  has 0.47 / 0.49 of the expected cells). Should be `%in%`.
- `FullDataset_Analysis/3a_Fibroblast_Subclustering.R` — `SFPR4` (typo for
  SFRP4), silently dropped from the COL8A1 panel.
- wasikowr's clusters were computed on the UMAP, not on PCs.
- Harmony on `subject` in most compartment scripts (see table above).

## Questions for wasikowr / jarnagin

1. wasikowr: the marker list and cluster → label map for the 16 cell types.
2. wasikowr: SoupX settings (rho per sample) and the Harmony variable.
3. jarnagin: the final label tables (barcode → label) for immune, T, B/plasma,
   fibroblast, to compare with ours.
