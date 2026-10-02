# Clinical and bulk sample metadata — read this before using them

**The master table is the core of all downstream analysis.** The three
workbooks below were cleaned on 2026-10-02 and are the files to use: read
them directly, with no rename maps or correction steps. Patient-level: never in
git; AI agents work through scripts only and see counts, not values
(decisions.md 2026-10-02). This file holds structure, rules and counts only.

## Files

Directory `data/clinical/` = `/home/parvizi/asset-data/data/clinical/`
(`CLINICAL` in config/paths.R; the master is `CLINICAL_MASTER`).

| File | Sheet | Rows | Unit |
|---|---|---|---|
| `ASSET_clinical_data_master_subject_based.xlsx` | `Sheet1` | 88 | one row per subject (44 Placebo / 44 Abatacept) |
| same | `notes_extra_data` | 10 | free-text notes (`topic`, `note`) |
| `ASSET_clinical_data_skin_bulk_metadata.xlsx` | `ASSET_skin_bulk` | 234 | one row per skin bulk sample |
| `ASSET_clinical_data_PBMC_bulk_metadata.xlsx` | `ASSET_blood_combined_filtered_7` | 70 | one row per PBMC bulk sample (all baseline) |

`data/clinical/backup/R0033/` holds the three workbooks as they were before
the 2026-10-02 edits (read-only), plus the record of the one matrix column
rename. Nothing else is needed from that folder.

Bulk matrices: `data/bulk/skin/` — 6 matrices, columns = skin `Sample_ID`
(all 234 match exactly); the Ensembl TPM file has one extra annotation
column. `data/bulk/blood/` — 4 matrices, columns = PBMC `Sample_ID` (all 70
match), first column blank-named (genes). Column names, types and
missing/distinct counts of the workbooks: re-run 00_4 ->
`metadata/design/clinical_column_inventory.csv`.

## How to filter — one column

**`Use_downstream == "Yes"`** in both bulk workbooks (reason for "No" in
`Use_reason`):

| | Use_downstream = Yes | What it removes |
|---|---|---|
| skin | **219 samples**, exactly one per Subject x Timepoint; 84 subjects: Baseline 80, Month3 69, Month6 70 | owner's `Exclude` (11) and `Repeat` (4 more; 8 flagged), and QC PCA outliers (all 9 already in Exclude) |
| PBMC | **67 samples**, one per subject | QC_flag `exclude_suggested` (2); the second sample of the one subject with two (technical replicate; kept: QC `ok` over `check`, then deeper library) |

**Baseline availability per patient**: master column `Bulk_baseline_data`
(`skin+PBMC`, `skin only`, `PBMC only`, `none`), from the usable samples above.

| Bulk_baseline_data | Abatacept | Placebo | Total |
|---|---|---|---|
| skin+PBMC | 29 | 34 | 63 |
| skin only | 14 | 3 | 17 |
| PBMC only | 0 | 4 | 4 |
| none | 1 | 3 | 4 |

Placebo by mRSS_category (Improver / Stable / Worsened / Set_aside):
skin+PBMC 16 / 10 / 6 / 2; skin only 0 / 1 / 2 / 0; PBMC only 1 / 2 / 0 / 1;
none 2 / 0 / 0 / 1.

## Keys and codings

- `Subject_ID` (`99-9999`) joins all three workbooks and the scRNA-seq design
  table. Every skin and PBMC subject is in the master.
- Skin `Sample_ID` = `<Subject_ID with _>_<Timepoint>`; the 8 second samples
  carry an `_A` suffix and `replicate = 2`. `Timepoint` Baseline / Month3 /
  Month6 = `Time` 0 / 3 / 6. `ASSETpaper_2022` marks the 140 samples of the
  2022 paper.
- PBMC `Sample_ID` = `A_99_9999` (one with a letter suffix: the replicate).
  `Matched_Skin_Sample_ID` = that subject's baseline skin `Sample_ID` (one
  PBMC row has none).
- Clinical scores in the skin sheet (`MRSS`, `HAQ_DI`, `PGA_Physician`,
  `PtGA_Patient`, `FVC`, `DLCO`) equal the master's value at that timepoint
  (screening values for Time 0). Master time-indexed columns: mRSS `0 1 3 6
  9 12` (use these, not `MRSSTOTALCALCULATION_*`), `FVCPTP_*`,
  `DLCO_corrected_*`, `HAQ_DI_Month*`, `PGA_Physician_Month*`,
  `PtGA_Patient_Month*`, `Days_*`.
- Codings match the master: Sex `Female`/`Male`; Improver / True_improver
  `Improver` / `Non-Improver` / `No_data`. Exception: `Ever_escaped` is
  `True`/`False` in the master and `TRUE`/`FALSE` in the skin sheet.
- **Blank cells are often the text `NA`** — treat `"NA"` as missing before
  testing any column.

## Bulk skin QC columns (01_1, R0036)

Same definitions and rules as the owner's PBMC QC (`1.4_QC_70_samples.R`,
with the blood data): `raw_library_size` (Ensembl raw counts),
`corrected_library_size` and `cor_before_vs_after_combat` (blank — no
ComBat-seq for skin), `genes_detected_TPM1`, `hemoglobin_TPM_pct`,
`mito_TPM_pct`, `median_cor_to_others`, `PC1`, `PC2`, `PCA_distance_PC1to5`
(log2 TPM+1, 17,151 genes with TPM >= 1 in >= 20% of samples), `Note`,
`QC_flag`. Result: ok 215, check 10, exclude_suggested 9 — all 9 already in
the owner's Exclude; 5 kept samples are "check" (reason in `Note`). The
2022-paper set explains <= 3.4% of any of PC1-5 (no batch split). Re-running
01_1 rewrites these columns.

## Open (owner)

- **Autoantibody columns**: many statuses being confirmed externally; the
  skin sheet also differs from the master for `Antibody_category_lab`
  (4 subjects differ, 2 blank). Not edited.
- **Escape columns**: skin `Ever_escaped` blank for 3 subjects that the master
  fills; `Escape_6mo` differs for 1 subject. Not edited.
  ID-level list for both: `data/patient_level/freeze01/reference/`
  `00_manifest_and_batch/resolve_bulk_metadata_mismatches/disagreement_details.csv`.
- **Which symbol-level raw count matrix is canonical**: two exist
  (39,921 and 40,279 genes), no README for skin. Proposal: start from the
  Ensembl raw counts and map to symbols in our own code.

## What was changed on 2026-10-02 (one-time; scripts removed afterwards)

Done by one-time scripts (recoverable from git at the sha in docs/runs.csv):
skin coding harmonised (Sex, Improver, True_improver); skin `mRSS_category`
set to the master's value (1 subject, 2 rows: Set_aside -> Worsened); QC
and `Use_downstream` columns added to skin (01_1) and `Use_downstream` to
PBMC; `Bulk_baseline_data` added to the master; one misspelt skin `Sample_ID`
column header corrected in all 6 skin matrices (the metadata was right).
Runs: R0031-R0033, R0036, R0037, R0040.
