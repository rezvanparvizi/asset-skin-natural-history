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
