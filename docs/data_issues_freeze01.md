# freeze01 — data issues to resolve with colleagues

Prepared 2026-09-29 from the import of the 221 libraries in the
colleague's sample map (`00_0_import_raw_libraries`). Library IDs only;
patient-level detail is in
`data/patient_level/freeze01/reference/00_manifest_and_batch/import_raw_libraries/`.

Each section ends with the questions to ask. Record answers in
`docs/decisions.md` and, where they change metadata, in
`config/freezes/`.

---

## 1. Healthy controls were processed apart from the SSc samples

All 12 controls sit in pools that contain no SSc sample, so every
SSc-vs-control difference is partly technical by construction.

| Pool | Controls | CellRanger | Cell calling | Processed by |
|---|---|---|---|---|
| 12952-JF | 7 (`12952-JF-10` .. `-16`) | 9.0.1 | automatic | sequencing core |
| 8911-JF | 4 (`8911-JF-13` .. `-16`) | 9.0.0 | **forced** (`force_cells` 5000 / 500 / 3500 / …) | wasikowr |
| 9260-JF | 1 (`9260-JF-16`) | 9.0.0 | **forced** (1000) | wasikowr |

Same probe panel everywhere (Human Transcriptome v1.1.0, identical
file) and same reference (GRCh38-2024-A).

Questions:
- Were the controls fixed, stored and hybridised with the same protocol
  and in the same period as the ASSET biopsies? Where do they come from?
- Why were 8911-JF and 9260-JF re-run with `force_cells`, and how were
  the numbers chosen? (`8911-JF-14` forced to 500 cells.)
- Pool 12952-JF has 16 libraries; only 7 are in the sample map. What are
  `12952-JF-1` .. `-9`?
- Are there more controls anywhere, ideally captured alongside SSc
  samples? Even two or three would let the control/SSc batch effect be
  estimated.

## 2. Runs, pools and batches

- 15 SSc pools (one 10x Flex pool of up to 16 probe-barcoded samples =
  one GEM capture). Every pool mixes baseline, month 3 and month 6 and
  both arms — good for the longitudinal design.
- Two pools are lopsided by arm: **13678-DP** (14 placebo / 2
  abatacept) and **13679-DP** (3 placebo / 13 abatacept).
- Sample names start with a date (YYMMDD): 9 dates from 2025-04-22 to
  2025-05-30, each covering 1–3 pools.
- The colleague's `batch` column is unique per sample, so it is not a
  batch.

Questions:
- Is the pool the capture batch? Were the pools prepared on different
  days from when they were sequenced?
- What does the date at the start of the sample name mean — fixation,
  sectioning, hybridisation, capture?
- How were samples assigned to pools? Was it randomised, or by patient
  or site?
- **Site:** which hospital collected and handled each biopsy? (Coming
  from the clinical table; it will be added to the batch check.)

## 3. Libraries sequenced but absent from the sample map

Present in the ASSET pools' `config.csv`, not in the sample map, so not
in freeze01:

- `13452-JF-1` .. `-7` (7 libraries; the only pool with 9 of 16 in the map)
- `13454-JF-7`, `13484-JF-11`, `13596-JF-4`, `13634-JF-13`,
  `13639-JF-4`, `13639-JF-10`, `13719-JF-8`, `13719-JF-15`

Questions:
- Are these ASSET samples that failed, pilot or technical samples, or
  another study?
- If they are failed ASSET samples: which arm and timepoint? Libraries
  missing for a reason linked to arm or timepoint are informative
  missingness, not just QC.
- Are any of them among the ~20 samples expected in freeze02?

## 4. Timepoint labelling

`13639-JF-11` and `13639-JF-12` belong to the same (abatacept) subject.
Both are labelled Month 3, and that subject has no Baseline. One is
probably a mislabelled Baseline.

Current handling: both stay in freeze01 for Layer 1, flagged
`timepoint_unresolved` in `config/freezes/freeze01.yml`, and are
excluded from every inference cohort until resolved.

Question: which one is the Baseline biopsy? (Check the DCC sample log.)

## 5. Libraries with CellRanger metrics outliers (preliminary)

These are flags only. The formal decisions come at QC (stage 03).

| Library | Cells called | Median genes/cell | Median UMI | Mean reads/cell | % reads mapped in cells | Likely issue |
|---|---|---|---|---|---|---|
| 13718-JF-14 | 87,685 | 3 | 3 | 103 | 98.2 | cell calling failed (empty droplets called as cells) |
| 13719-JF-10 | 1,144 | 7 | 7 | 256 | 15.9 | library failed |
| 13596-JF-9 | 17,557 | 102 | 107 | 1,570 | 80.3 | very low depth / many low-quality barcodes |
| 9260-JF-16 | 1,000 (forced) | 283 | 338 | 3,250 | 48.7 | low quality control library |
| 13481-JF-8 | 1,025 | 408 | 504 | 10,691 | 33.2 | low fraction of reads in cells |
| 13452-JF-8 | 2,160 | 546 | 689 | 4,103 | 39.3 | low fraction of reads in cells |

Across all 221 libraries: median 4,815 cells, 952 genes/cell, 1,283 UMI/cell.

Questions:
- Are any of these being re-sequenced or re-made?
- Were they already known to the colleague, and excluded in her
  analysis?

## 6. The upstream object and ambient correction

- `wasikowr/ffpe/Asset/seurat.RDS` was rewritten on 2026-09-02, after
  the colleague's objects were built from it (May–June). A `Batch5`
  folder appeared on 2026-08-26.
- Its counts layer is SoupX-corrected and non-integer. We are re-running
  SoupX from the raw per-sample matrices.

Questions for wasikowr:
- What changed in the September rewrite? Which samples does `Batch5` add?
- What SoupX settings were used, and which rho for each sample? We will
  set our own rho by a written rule, but hers are a useful comparison.
- Was doublet detection run? If so, with what?

## 7. Points on the colleague's code (to pass on)

- `Baseline_Analysis/1a_Baseline.R`:
  `subset(orig_obj, time == c("Baseline","Control"))` compares against a
  two-item vector, so R recycles it and cells are tested alternately
  against "Baseline" and "Control". It probably drops about half the
  cells. It should be `time %in% c(...)`. A quick check:
  `table(base_control$time)` against `table(orig_obj$time)`.
- Harmony in `FullDataset_Analysis/2a` and `3a` integrates on `subject`,
  which in the sample map is unique per sample. That removes the
  differences between timepoints within a patient.

## 8. freeze02 (~20 new samples, about a month away)

Questions:
- Which subjects, arms and timepoints? Are they new pools, and were
  they processed with the same protocol and CellRanger version?
- Can their sample map come in the same format, with the pool and the
  sample-name date?

## Subject ID typos in the sample map (found 2026-10-01)

Two `patient` values in the scRNA-seq sample map do not match the
clinical master table; the owner confirmed both as typos. Corrected in
this repository via metadata/design/subject_id_corrections.csv. The
colleague's objects (and anything built from the sample map) still carry
the wrong IDs — worth fixing at the source.

## Lab code issues found while reviewing the cell-typing routine (2026-10-01)

Details in docs/lab_routine_celltyping.md, section 4: the recycling
subset bug in jarnagin's 1a_Baseline.R (baseline object holds about half
the cells), the SFPR4 typo, wasikowr's clustering on the UMAP, and
Harmony on `subject`. Also ask wasikowr for her top-level marker list
and SoupX settings.

## 23 newer libraries already processed by wasikowr (found 2026-10-01)

wasikowr's seurat.RDS (2 Sep) has 242 samples: our 219 plus 23 newer
libraries, listed in /hits/home/wasikowr/ffpe/Asset/firstrun/asset_batch5.csv
(sampleID, path, batch = original sample number, condition). CellRanger
outputs (re-run by wasikowr, Aug 2026; flat per_sample_outs layout, no
count/ subfolder): /hits/home/wasikowr/ffpe/{14794,15195,15231,15242}-JF_v1/.

- 21 are repeats of freeze01 libraries (matched by sample number), mostly
  low-yield originals (~1,000-2,000 cells). Same biopsy -> a second
  library_id for one sample_id; merge_strategy decision needed.
- 2 (sample numbers 32 and 68) have no original in our sample map.
  Subject/timepoint unknown — ask wasikowr.
- None repeats our 3 failed libraries.
These are freeze02 material. Asked wasikowr (Slack, 2026-10-01) for a
sample sheet for all 242 samples, CellRanger paths, per-cell labels and
how many samples are still waiting to be sequenced.

## Cross-check of sample records: ours vs wasikowr vs jarnagin (2026-10-01)

- Sample map: wasikowr's firstrun/demo.txt and jarnagin's results/demo.txt
  are identical (same md5); it is our freeze01 source. wasikowr's
  firstrun/input.txt = our 209 SSc libraries.
- Samples: ours 221 imported / 218 in freeze01; jarnagin 219 (all
  timepoints; = ours + 13596-JF-9, minus our 2 zero-cell failed
  libraries); wasikowr 2 Sep object 242 = those 219 + 23 newer.
- Subject IDs: jarnagin's `patient` column = our Subject_ID + trailing
  "B" for all 207 SSc libraries, once our two corrections are applied —
  her metadata carries the same two typos.
  Timepoint and arm agree for every library.
- jarnagin's `subject` column is a per-biopsy SPECIMEN code (site
  prefix + number), unique per library; the prefix is the site. Specimen numbers
  are NOT chronological within patient (4 subjects out of order), so they
  cannot resolve 13639-JF-11/-12 (two different specimens, both "Month 3";
  the patient's M06 specimen has a higher number than both).
- Sample-map gaps: 14 sample numbers are absent (32, 68, 77, 118, 125,
  131, 198, 203, 208, 209, 211, 214, 215, 222). The first six fall
  between consecutive libraries of one pool, i.e. most likely the pool
  libraries missing from the sample map: 13454-JF-7, 13484-JF-11,
  13596-JF-4, 13634-JF-13, 13639-JF-4, 13639-JF-10 (inferred, not
  confirmed). Their per-sample outputs were not copied into data/raw.
- Repeats (asset_batch5.csv): `batch` = original sample number.
  15 of 15 repeats that carry a specimen code match the original
  library's specimen -> same biopsy, confirmed. 6 repeats (15195-JF-11..16)
  carry no code -> same biopsy by sample number only. 2 repeats (#32,
  #68) have no original in the sample map; by the gap
  analysis they repeat 13454-JF-7 and 13484-JF-11. Patient/timepoint
  for these two unknown.
- New runs: every library forced (force_cells), unlike freeze01 SSc
  pools; probe set v1.1.0 (same as freeze01; a v2-probe run of 15242-JF
  was superseded by 15242-JF_v1); probe barcodes of a different series
  (BC0xx / A-, B-, D- plate wells) than freeze01 — chemistry version to
  confirm with the core. 15231-JF_v1 is a combined pool with 15230-JF.

## Agenda for the conversation with wasikowr (2026-10-02)

Her reply (2026-10-01): she is finishing the latest build; was puzzled
that the cell type metadata is missing; objects need library(BPCells).

1. **Cell types.** The 2 Sep seurat.RDS metadata has only orig.ident,
   nCount_RNA, nFeature_RNA, batch, orig.nb, RNA_snn_res.0.1,
   seurat_clusters. The 16 `celltype` labels exist in the earlier version
   (jarnagin's objects built from it, June-Aug). Ask: will the new build
   carry celltype for all cells, and can she export a per-cell table
   (barcode, orig.ident, celltype) plus the marker genes / cluster ->
   celltype step?
2. **Samples.** Confirm 242 = 219 + 23 newer (asset_batch5.csv). Patient
   and timepoint for #32 and #68; are they repeats of
   13454-JF-7 and 13484-JF-11? Are 15195-JF-11..16 repeats of the same
   biopsies as their sample numbers say? Will the new build add more?
3. **Missing from the sample map.** Why were 13454-JF-7, 13484-JF-11,
   13596-JF-4, 13634-JF-13, 13639-JF-4, 13639-JF-10 (and 8 other sample
   numbers) left out — failed, or withdrawn?
4. **Processing of the newer runs.** force_cells for every library — why;
   the probe-barcode series differs from freeze01 (chemistry version?).
5. **Her QC / SoupX.** Filter used (we see > 200 genes); SoupX rho per
   sample; Harmony variable; why ~2,800 high-UMI cells present in the raw
   outputs are absent from her object.
6. **How many samples are still waiting to be sequenced.**
7. **13639-JF-11 / -12** (one patient, both labelled Month 3): which is
   Baseline?
