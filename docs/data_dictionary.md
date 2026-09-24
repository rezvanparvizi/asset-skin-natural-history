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
| `subject_id` | chr | unique per participant. Multiple samples per subject (M00/M03/M06). |
| `group` | chr | `SSC` or `HC` |
| `arm` | chr | `PLACEBO`, `ABATACEPT`, or `NA` for healthy controls |
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

## data/clinical/asset_clinical.csv  [NEVER in git]

From the U-M Data Coordinating Center. Joined on `subject_id`.

| Column | Meaning |
|---|---|
| `subject_id` | |
| `arm` | |
| `disease_duration_yr` | at baseline; separates matrix resolution from late atrophy |
| `mrss_m00` ... `mrss_m12` | global mRSS at 0/1/3/6/9/12 months |
| `mrss_forearm_m00` ... | **site-level forearm score if available.** The biopsy is forearm; local score tracks local biology better than global mRSS. Ask the DCC. |
| `improver_12m` | `improver` / `non_improver`: >=5-point OR >20% mRSS reduction at M12 (definition used in Khanna 2020, Mehta 2022) |
| `delta_mrss_12m` | continuous; better powered than the binary version |
| `escape_start_month` | month escape therapy began; `NA` if never. Drives cohort censoring. |
| `rnap3` | anti-RNA polymerase III status. Reported inconsistently across the ASSET papers (40%/51% by arm; 39/85; 32/64) — **verify against the DCC before using as a covariate.** |
| `scl70`, `aca` | other autoantibodies |
| `tfr` | tendon friction rubs at baseline |
| `intrinsic_subset` | inflammatory / fibroproliferative / normal-like, from bulk skin (GSE217067) |
| `cd28_score_m00` | CD28 costimulation module score from bulk skin |
| `hhaq_di`, `criss`, `fvc_pct`, `dlco_pct`, `ptga`, `phga` | secondary outcomes |

### Outcome definition is an open decision

Not yet settled, and it matters:

- **global mRSS vs forearm-site score** — the biopsy is forearm. THBS1
  correlates with local skin score at r = 0.76, more weakly with global
  mRSS.
- **binary improver vs continuous delta** — binary is comparable to
  papers 01 and 03; continuous is better powered with ~19 vs ~19.

Record the choice in `docs/decisions.md` when made.
