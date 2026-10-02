# Clinical and bulk sample metadata — where it is, how to filter it

**The master table is the core of all downstream analysis.** Patient-level
files: never in git, never opened by AI agents directly — agents work
through scripts (decisions.md 2026-10-02, "scripts-only access"). This file
holds structure, rules and counts only; no IDs, no values.

## Where the files are

Directory: `data/clinical/` = `/home/parvizi/asset-data/data/clinical/`
(`CLINICAL` in config/paths.R; the master is `CLINICAL_MASTER`).

| File | Sheet | Rows x columns | Unit |
|---|---|---|---|
| `ASSET_clinical_data_master_subject_based.xlsx` | `Sheet1` | 88 x 76 | one row per subject (44 Placebo / 44 Abatacept) |
| same | `notes_extra_data` | 10 x 2 (`topic`, `note`) | free-text notes |
| `ASSET_clinical_data_skin_bulk_metadata.xlsx` | `ASSET_skin_bulk` | 234 x 46 | one row per skin bulk sample |
| `ASSET_clinical_data_PBMC_bulk_metadata.xlsx` | `ASSET_blood_combined_filtered_7` | 70 x 18 | one row per PBMC bulk sample (baseline) |

Backups made by scripts before any in-place edit: `data/clinical/backup/<run_id>/`.

Bulk matrices: `data/bulk/skin/` (6 matrices, 234 sample columns; Ensembl
TPM has one extra annotation column) and `data/bulk/blood/` (4 matrices,
70 sample columns + blank-named gene column; QC scripts 1.3/1.4, README).
Column inventories (names, types, missing/distinct counts) are regenerated
by 00_4 into `metadata/design/clinical_column_inventory.csv`,
`clinical_file_inventory.csv` and `bulk_matrix_inventory.csv`.

## Keys

- `Subject_ID` (shape `99-9999`) — joins all three files and the scRNA-seq
  design table.
- Skin `Sample_ID` (`99_9999_<Timepoint>`, `_A` suffix on a replicate) =
  the skin matrix column names.
- PBMC `Sample_ID` (`A_99_9999`, one with a letter suffix) = the blood
  matrix column names. `Matched_Skin_Sample_ID` = the skin `Sample_ID` of the
  same subject's baseline biopsy (69 of 70; one row has no match).
- Skin `Timepoint` (`Baseline` / `Month3` / `Month6`) and `Time` (0 / 3 / 6)
  are consistent. Clinical scores in the skin sheet (`MRSS`, `HAQ_DI`, `PGA`,
  `PtGA`, `FVC`, `DLCO`) are the master's value at that timepoint
  (all 234 rows agree; screening values for Time 0).
- Master time-indexed columns: mRSS `0 1 3 6 9 12` (also
  `MRSSTOTALCALCULATION_*`), `FVCPTP_*`, `DLCO_corrected_*`, `HAQ_DI_Month*`,
  `PGA_Physician_Month*`, `PtGA_Patient_Month*`, `Days_*`.

## Codings

Harmonised in the skin workbook on 2026-10-02 (00_7, R0033) to the
master's coding: Sex `Female`/`Male` (was `F`/`M`); Improver and
True_improver `Non-Improver` (was `Non_Improver`). Still different (escape
left out by the owner): Ever_escaped `True`/`False` (master) vs
`TRUE`/`FALSE` (skin).

**Blank cells are often the text `NA`** in these workbooks. Any script
must treat `"NA"` as missing before testing a flag (R0035 crashed on this).

## How to filter the skin samples (owner's rule)

1. Drop `Exclude == "Exclude"` (excluded / low quality; owner-curated):
   11 samples.
2. Drop `Repeat == "Repeat"`: 4 more (8 samples carry the flag; 4 of them are
   already excluded at step 1).
3. Result: **219 samples = exactly one per Subject x Timepoint**, 84
   subjects: Baseline 80, Month3 69, Month6 70.

| | samples | subject x timepoint units | units with > 1 sample | Baseline / M3 / M6 subjects |
|---|---|---|---|---|
| all | 234 | 226 | 8 | 84 / 70 / 72 |
| after Exclude | 223 | 219 | 4 | 80 / 69 / 70 |
| after Exclude + Repeat | 219 | 219 | 0 | 80 / 69 / 70 |

`replicate` (1 / 2; 8 samples = 2) and the `_A` Sample_ID suffix mark the
same 8 second samples. `ASSETpaper_2022` marks the 140 samples used in the
2022 paper. **Use `Use_downstream == "Yes"`** (added 2026-10-02 by 01_1, R0036): it
applies Exclude, Repeat and the QC rule in one column, with the reason in
`Use_reason`. Same 219 samples as the owner's two filters.

PBMC: `QC_flag` = `ok` 62, `check` 6, `exclude_suggested` 2 (owner's QC).
One subject has two PBMC samples.

## Coverage (from 00_5, R0031)

- Master subjects with >= 1 skin bulk sample: 84 / 88; with PBMC: 69 / 88.
- All skin and PBMC Subject_IDs exist in the master.
- Master column **`Bulk_baseline_data`** (01_2, R0037): usable BASELINE bulk
  data per patient. Skin = a Baseline sample with `Use_downstream == "Yes"`;
  PBMC = a PBMC sample whose `QC_flag` is not `exclude_suggested`.

| Bulk_baseline_data | Abatacept | Placebo | Total |
|---|---|---|---|
| skin+PBMC | 29 | 34 | 63 |
| skin only | 14 | 3 | 17 |
| PBMC only | 0 | 4 | 4 |
| none | 1 | 3 | 4 |

Placebo by mRSS_category (Improver / Stable / Worsened / Set_aside):
skin+PBMC 16 / 10 / 6 / 2; skin only 0 / 1 / 2 / 0; PBMC only 1 / 2 / 0 / 1;
none 2 / 0 / 0 / 1. QC filters change: Baseline skin 84 -> 80 patients
(owner's Exclude), PBMC 69 -> 67 (exclude_suggested).

## Inconsistencies found (2026-10-02; 00_5 R0031, 00_6 R0032) — details for the owner

ID-level lists: `data/patient_level/freeze01/reference/00_manifest_and_batch/`
`check_bulk_sample_metadata/consistency_details.csv` and
`resolve_bulk_metadata_mismatches/disagreement_details.csv`.

| What | Rows | Subjects |
|---|---|---|
| skin Sample_ID missing from all 6 skin matrices; the matrices carry one unexplained column; edit distance 2, same ID shape — very likely a typo of the same sample | 1 | 1 |
| skin `mRSS_category` differs from master | 2 | 1 |
| skin `Antibody_category_lab` differs from master | 12 | 4 |
| skin `Antibody_category_lab` blank, master filled | 6 | 2 |
| skin `Ever_escaped` blank, master filled | 12 | 3 |
| skin `Escape_6mo` differs from master | 3 | 1 |
| PBMC: one row with no matched skin sample; the duplicated subject's two PBMC samples share one `Matched_Skin_Sample_ID` | 1 / 1 | |
| PBMC rows whose matched skin sample is Exclude- or Repeat-flagged | 3 | |

Everything else agrees: arm, age, autoantibodies, escape 3/9/12 mo, all
timepoint scores, PBMC site.

## Bulk skin QC (01_1, R0036) — added to the skin workbook

Columns added (same definitions and rules as the owner's PBMC QC,
`analysis/02_bulk_pbmc/1.4_QC_70_samples.R`, gitignored because it holds
IDs): `raw_library_size`, `corrected_library_size` (blank: no ComBat-seq for
skin), `genes_detected_TPM1`, `hemoglobin_TPM_pct`, `mito_TPM_pct`,
`median_cor_to_others`, `cor_before_vs_after_combat` (blank), `PC1`, `PC2`,
`PCA_distance_PC1to5`, `Note`, `QC_flag`, `Use_downstream`, `Use_reason`.
Expressed genes: TPM >= 1 in >= 20% of samples (17,151). Library size from
the Ensembl raw counts; the rest from the symbol TPM matrix.

Result: QC ok 215, check 10, exclude_suggested 9. **All 9 PCA outliers are
already in the owner's Exclude**; the other 2 owner exclusions are QC
"check". No owner-kept sample is an outlier; 5 kept samples are "check"
(kept; reason in Note). Rules triggered: low similarity 9, PCA outlier 9,
few genes 16, low depth (< 3M) 5, many genes 1. The 2022-paper set explains
<= 3.4% of any of PC1-5 (no batch split).

## Matrix column fix

One skin matrix column is misspelt (the metadata is right: it matches the
ID built from Subject_ID + Timepoint). Rename map:
`data/clinical/bulk_skin_matrix_column_fixes.csv` — every script that
reads the skin matrices must apply it (01_1 does).

## Edit log of the owner's workbooks (backups in data/clinical/backup/<run_id>/)

| Run | Script | Workbook | Change |
|---|---|---|---|
| R0033 | 00_7 | skin | Sex, Improver, True_improver coding; mRSS_category from the master (1 subject, 2 rows: Set_aside -> Worsened) |
| R0036 | 01_1 | skin | 14 QC / use columns added |
| R0037 | 01_2 | master | `Bulk_baseline_data` added |

Not changed, pending the owner: autoantibody columns (being confirmed
externally; many statuses wrong), escape columns.
