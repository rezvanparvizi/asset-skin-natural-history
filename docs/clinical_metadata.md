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

## Codings that differ between files (harmonise before joining)

| Field | master | skin / PBMC sheets |
|---|---|---|
| Sex | `Female` / `Male` | `F` / `M` |
| Improver, True_improver | `Non-Improver` | `Non_Improver` |
| Ever_escaped | `True` / `False` | `TRUE` / `FALSE` |

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
2022 paper. Final-use column after the skin QC: **pending** (item 4 of the
2026-10-02 request).

PBMC: `QC_flag` = `ok` 62, `check` 6, `exclude_suggested` 2 (owner's QC).
One subject has two PBMC samples.

## Coverage (from 00_5, R0031)

- Master subjects with >= 1 skin bulk sample: 84 / 88; with PBMC: 69 / 88.
- All skin and PBMC Subject_IDs exist in the master.
- Per-patient availability column in the master: **pending** (item 5).

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
