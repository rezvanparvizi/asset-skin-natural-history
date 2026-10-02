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
| `Asset/Batch1`..`Batch5`, `Batch*a`, `bpcells_obj*` | Per-batch working folders. `Batch1` holds `ffpe_report.pdf` (95 MB), `harmony.RData` (47 GB), `markers.csv`. Reports **not read** (no PDF tools on the server) — may document her settings. | 2026-10-02 (listing only) |
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
Not readable: `/hits/home/jarnagin/ASSET` (permission denied),
`/hits/home/jarnagin/freezes` (mode 750).

| Object (under objects/) | Cells | Counts | Labels | Checked |
|---|---|---|---|---|
| `ASSET_base_control.RDS` (+ `base_control_BPCells/`) | 181,766 | wasikowr's (corrected) | `celltype` (wasikowr's 16 broad labels) | 2026-10-01, 10-02 |
| `ASSET_base_control_integer_BPCells.RDS` | same 181,766 | raw CellRanger h5 (her `1b`) | as above | 2026-10-02 |
| `ASSET_Round2_immFilt.RDS` | 104,675 | wasikowr's | full-dataset immune `celltype`, `subclusters` | 2026-10-01, 10-02 |
| `ASSET_FullDataset_Tcells_Subtyped.rds` | 25,322 | — | T subtype / lineage calls | 2026-10-01 |
| `ASSET_FullDataset_BcellPlasma_Clean.rds` | 6,582 | — | B/plasma `subclusters` | 2026-10-01 |
| `ASSET_Round1_fibs.RDS` | 143,963 | — | fibroblast `celltype` | 2026-10-01 |
| `ASSET_MuSiC_reference_sce.rds` | 82,369 | — | many (subtype, lineage, compartment) | 2026-10-01 |
| `ASSET_base_control_cohort_barcodes.txt` | — | — | barcodes of the base_control cohort (used by her 4a) | 2026-10-02 |

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
- Harmony on `subject` (= per-biopsy specimen code) in compartment scripts.
- Her metadata `patient` = our Subject_ID + "B"; carries the same two
  sample-map typos we corrected.
- Which object is current per compartment (`res0p2`/`res0p3`/`from3b`/
  `_v2`/`Round1`/`Round2`) is still unanswered (`inherited_objects.md`).

## What we derived from these sources

| Output | Where | Run |
|---|---|---|
| Barcode -> label table: our cell IDs joined to labels from the 7 objects above (wasikowr `seurat.RDS` + six jarnagin objects) | `metadata/inherited_labels/freeze01_upstream_labels.csv.gz` (barcode-level, gitignored; agents do not open it) | R0021 (06_0) |
| Match and coverage summaries (aggregate) | `results/freeze01/reference/06_sc_lineage/match_upstream_cells_and_labels/tables/` — `upstream_sources.csv`, `cellset_overlap_per_library.csv`, `inherited_label_coverage.csv`, `inherited_celltype_by_timepoint.csv`, `wasikowr_vs_jarnagin_labels.csv`, ... | R0021 |
| Bug impact, wasikowr counts vs raw, per-gene check | numbers in `data_issues_freeze01.md` and `ambient_soupx_freeze01.md`; scratch scripts, aggregate output only, not kept | 2026-10-02 |
