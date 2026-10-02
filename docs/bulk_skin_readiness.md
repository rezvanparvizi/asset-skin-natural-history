# Bulk skin RNA-seq (stage 01) — what is present, what is missing

Checked 2026-10-02 from `metadata/design/bulk_file_inventory.csv` (file
names and sizes only; no file under data/ was opened). Numbers in file
names are read as genes x samples; that reading is an inference from the
names, not verified against the contents.

## Present in data/bulk/skin/

| File | Reads as | Size |
|---|---|---|
| `raw_counts_mtx_ASSET_skin_EMSEMBL_64252_234.txt` | raw counts, Ensembl IDs, 64,252 genes x 234 samples | 41 MB |
| `raw_counts_mtx_ASSET_skin_geneSymvol_39921_234.txt` | raw counts, gene symbols, 39,921 genes | 28 MB |
| `raw_counts_mtx_ASSET_skin_geneSymvol_40279_234.txt` | raw counts, gene symbols, 40,279 genes | 28 MB |
| `TPM_exp_mtx_ASSET_skin_ENSEMBL_64252_234.txt` | TPM, Ensembl | 167 MB |
| `TPM_exp_mtx_ASSET_skin_geneSymbol_39921_234.txt` | TPM, gene symbols | 118 MB |
| `logTPM_exp_mtx_ASSET_skin_geneSymbol_39921_234.txt` | log TPM, gene symbols | 118 MB |

("geneSymvol" is the file names' spelling.)

## Missing — compare with data/bulk/blood/, which has all three

1. **Sample sheet**: column name -> Subject_ID, timepoint, arm, and any
   sequencing batch / site. Blood has
   `ASSET_blood_combined_filtered_70_sample_info.csv` and
   `batch_info_filtered_70_samples.txt`; skin has neither. Without it
   nothing in stage 01 can start, and `bulk_skin` in the design table stays NA.
2. **Provenance README**: aligner / quantifier, genome and annotation
   version, how Ensembl IDs were collapsed to symbols, any filtering.
   Blood has `README_ASSET_blood_70_samples.md`.
3. **Which symbol-level count matrix is canonical**: two exist (39,921
   and 40,279 genes; 358 more in the second). Most likely two different
   Ensembl -> symbol collapses (duplicate symbols summed vs dropped, or
   different annotation versions) — not checked.

## Questions (owner / DCC / whoever processed it)

1. Is there a sample sheet for the 234 skin columns (Subject_ID,
   timepoint, arm, batch/site)? If the column names already encode them,
   what is the format?
2. 234 samples: do they include healthy controls, technical repeats, or
   samples outside the 88 trial subjects? (AGENTS.md expects ~240.)
3. Is this the matrix behind GEO GSE217067 / Mehta 2022 (paper 3), so that
   the published intrinsic-subset calls and CD28 module scores map onto
   these columns?
4. Which of the two symbol-level raw count matrices was used for the
   published analyses?

## Proposal for stage 01 import (owner to decide)

Start from the **Ensembl raw counts** (64,252 genes) and do the
Ensembl -> symbol mapping in our own code, recorded in run_config; then
check it against both provided symbol matrices. Reason: reproducible and
independent of an undocumented collapse; the comparison tells us which of
the two provided matrices matches. Alternative: adopt whichever symbol
matrix the published analyses used (answer to question 4), for direct
comparability with papers 1 and 3.
