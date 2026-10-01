# Data dictionary

Field definitions for the manifest, the label tables, and the clinical
join. Update whenever a column is added.

## metadata/library_manifest.csv  [gitignored]

One row per sequencing library. Technical and design fields only — no
clinical values.

| Column | Type | Meaning |
|---|---|---|
| `library_id` | chr | unique per sequencing library. **Primary key.** |
| `sample_id` | chr | unique per biopsy. A resequenced biopsy has ONE `sample_id` and TWO `library_id`s. |
| `Subject_ID` | chr | unique per participant. Multiple samples per subject (M00/M03/M06). |
| `group` | chr | `SSc` or `HC` |
| `arm` | chr | `Placebo`, `Abatacept`, or `NA` for healthy controls |
| `timepoint` | chr | `M00`, `M03`, `M06`, or `NA` for HC |
| `batch_id` | chr | library prep / capture batch. **Load-bearing: hazard #1.** |
| `capture_date` | date | |
| `run_id` | chr | sequencing run |
| `cellranger_run` | chr | CellRanger version + reference used |
| `qc_status` | chr | `pass`, `fail_low_genes`, `fail_ambient`, ... |
| `resequenced_of` | chr | the `library_id` this one replaces, if any |
| `merge_strategy` | chr | `keep`, `exclude`, `replace`, `merge_reads` |
| `notes` | chr | free text |

### Why `library_id` and `sample_id` are separate

A resequenced sample yields a second library. You then have to decide,
per sample, whether to replace the original library or pool reads from
both. `merge_strategy` records that decision so it lives in the data
rather than in your head. Resequenced libraries will also differ from
the originals in depth *and* batch — and if the original failures were
not random with respect to arm or timepoint, freeze02 introduces a
confounder freeze01 did not have. Check it explicitly.

## metadata/labels/{freeze}_{labelset}_celltype_labels.csv.gz  [gitignored]

The Layer 1 -> Layer 2 interface. One row per cell.

| Column | Meaning |
|---|---|
| `barcode` | cell barcode, with library suffix as it appears in the object |
| `library_id` | joins to the manifest |
| `sample_id` | joins to the manifest |
| `lineage` | coarsest level, e.g. Stromal / Myeloid / TNK / BPlasma / Vascular / Keratinocyte |
| `celltype` | mid level, e.g. Fibroblast / Macrophage / CD8_T |
| `subtype` | finest level, e.g. SFRP4_FNDC1 / MMP_restorative |
| `qc_pass` | logical; cells excluded at QC are retained with FALSE rather than dropped |

Valid values for the three label columns are enumerated in
`metadata/celltype_dictionary.csv`. Do not introduce a label that is
not in the dictionary — add it there first.

## Nomenclature — used everywhere (tables, configs, plots, figures)

Values are spelled the way they appear in the publication. Never all caps.

| Concept | Column | Values |
|---|---|---|
| participant | `Subject_ID` | as in the clinical master table. Same name in every table. |
| disease group | `group` | `SSc`, `HC` |
| treatment arm | `arm` | `Placebo`, `Abatacept` (`NA` for HC) |
| timepoint | `timepoint` | `M00`, `M03`, `M06` |

The same nomenclature applies to bulk skin and bulk PBMC tables.

## ASSET cohort sizes

| | Subjects |
|---|---|
| Randomised | **88** (44 Placebo, 44 Abatacept) |
| Bulk skin RNA-seq | **84** subjects (no skin sequencing for 4 of 88) |
| Bulk PBMC RNA-seq, baseline | **69** subjects from **70** sequenced samples (one sample is a technical duplicate) |
| scRNA-seq skin (freeze01) | 79 SSc subjects + 12 HC (see `config/freezes/freeze01.yml`) |

Any stage that sees a different count stops and says so.

## Clinical master table  [NEVER in git]

`data/clinical/ASSET_clinical_data_master_subject_based.xlsx`
(i.e. `/home/parvizi/asset-data/data/clinical/`), one row per subject, all
88 randomised patients. Joined on `Subject_ID`. Column names only below;
values are never written into the repository.

**`mRSS_category` is the owner's primary outcome grouping.** Levels:
`Improver`, `Worsened`, `Stable`, `Set_aside`. The cut-offs defining
each level, and what puts a subject in `Set_aside`, are to be recorded
here before any analysis uses it.

| Columns | Meaning |
|---|---|
| `Number`, `Subject_ID` | row number; participant ID (join key) |
| `Treatment_arm` | `Placebo` / `Abatacept` |
| `Site_ID`, `Site_name`, `PI_name` | enrolling site |
| `True_improver`, `Improver` | improver calls. **The difference between the two is to confirm.** |
| `MRSS_Visit1`, `0`, `1`, `3`, `6`, `9`, `12`, `MRSS_MonthNA` | mRSS by visit month. **`0`, `1`, `3`, `6`, `9`, `12` are THE mRSS columns used in this project** (decisions.md, 2026-10-01). The bare numeric names become `X0`…`X12` under R's default `check.names`; read with `readxl` / `check.names = FALSE` and rename explicitly. |
| `MRSSTOTALCALCULATION_0` … `_12` | similar to the mRSS columns; meaning unknown. **Not used.** |
| `FVCPTP_screening`, `FVCPTP_Month3` … `_Month12`, `FVCPTP_MonthNA` | FVC % predicted |
| `Sex`, `Age` | |
| `Autoantibody_group`, `Autoantibody`, `Scl70`, `RNApol1`, `RNApol3`, `CENPB`, `Antibody_category`, `Antibody_category_lab` | serology. RNAP3 predicts placebo mRSS trajectory (paper 5). |
| `Note_on_Escape`, `Ever_escaped`, `Escape_month`, `Escape_3mo` … `Escape_12mo` | escape therapy. **`Escape_month` drives cohort censoring** (`cohort_samples()` in `R/io.R`); must be numeric months, empty if never escaped. |
| `Days_Screening`, `Days_Month0` … `Days_Month12`, `Days_MonthNA` | visit day |
| `PGA_Physician_Month0` … `_Month12`, `_MonthNA` | physician global assessment |
| `PtGA_Patient_Month0` … `_Month12`, `_MonthNA` | patient global assessment |
| `HAQ_DI_Month0` … `_Month12`, `_MonthNA` | HAQ-DI |
| `DLCO_corrected_Screening`, `_Month6`, `_Month12` | DLCO, corrected |
| `mRSS_change_6mo`, `mRSS_pctchange_6mo`, `mRSS_change_12mo` | mRSS change from baseline |
| `mRSS_category`, `mRSS_category_note` | **primary outcome grouping**: `Improver` / `Worsened` / `Stable` / `Set_aside` — see above |

Not in this table, and still needed: disease duration (separates matrix
resolution from late atrophy), tendon friction rubs, intrinsic subset and
CD28 score (from bulk skin, GSE217067), forearm site-level skin score if
the DCC has it.

### Outcome definition is an open decision

Not yet settled, and it matters:

- **global mRSS vs forearm-site score** — the biopsy is forearm. THBS1
  correlates with local skin score at r = 0.76, more weakly with global
  mRSS.
- **binary improver vs continuous delta** — binary is comparable to
  papers 01 and 03; continuous is better powered with ~19 vs ~19.

Record the choice in `docs/decisions.md` when made.
