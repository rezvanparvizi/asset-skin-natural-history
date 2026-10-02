# Cell and library QC, freeze01 — report

Runs: 03_1 cell QC (R0017), 03_2 loss by cell type (R0027), 03_3 loss by
cell type x timepoint (R0029). Library IDs only; no Subject_IDs.
**Re-run 03_1 -> 03_3 for every new freeze and compare with this file**
(freeze02 includes re-sequenced repeats of low-yield freeze01 libraries,
which may fix some of the libraries below).

## Rule (decisions.md 2026-10-01)

Per library: fail if log10 UMI or log10 genes < median - 3 MAD, genes < 200
or UMI < 500 (floor), log10 UMI > median + 5 MAD, or %MT > median + 3 MAD
and > 10%. 1,191,247 of 1,213,176 cells pass (98.2%). Three libraries
excluded from freeze01 at this step (13718-JF-14, 13719-JF-10, 13596-JF-9).

## Finding: library quality differs by timepoint

SSc libraries (03_1 library table):

| | M00 (71) | M03 (69) | M06 (66) |
|---|---|---|---|
| median UMI per cell, library median | 1,238 | 1,242 | 1,401 |
| cells called per library, median | 3,791 | 4,290 | 6,114 |
| floor failures, median library | 1.1% | 1.1% | 0.85% |
| floor failures, mean | 3.0% | 1.4% | 0.9% |
| floor-failing cells, total | 8,280 | 4,225 | 3,917 |
| libraries losing > 10% at the floor | 4 | 1 | 0 |

Median UMI by timepoint: Kruskal p = 0.001; floor fraction: p = 0.0007.
Pool is independent of timepoint (00_2), so this is not pool. Every cell
type (fibroblasts included) fails ~2x more at M00 than at M06 (03_3):
library quality, not cell-type biology. Cause unknown — observation (c).

## Libraries driving it

SSc libraries with > 5% floor failures or median UMI in the lowest 5%
(03_1 values, 2026-10-02):

| library_id | pool | timepoint | cells called | pass | floor fail | median UMI | median genes | 03_1 flag |
|---|---|---|---|---|---|---|---|---|
| 13719-JF-7 | 13719-JF | M00 | 6,256 | 3,123 | 50.0% | 500 | 433 | most_cells_fail |
| 13481-JF-8 | 13481-JF | M00 | 1,025 | 492 | 49.1% | 504 | 408 | few_passing_cells, most_cells_fail |
| 13454-JF-6 | 13454-JF | M00 | 3,249 | 2,617 | 19.3% | 717 | 609 | |
| 13452-JF-8 | 13452-JF | M03 | 2,160 | 1,817 | 15.2% | 689 | 546 | |
| 13680-DP-16 | 13680-DP | M00 | 3,489 | 2,920 | 13.4% | 752 | 567 | |
| 13680-DP-11 | 13680-DP | M00 | 3,709 | 3,483 | 6.1% | 1,070 | 827 | |
| 13634-JF-12 | 13634-JF | M03 | 4,769 | 4,637 | 2.3% | 850 | 670 | |
| 13453-JF-5 | 13453-JF | M03 | 1,634 | 1,600 | 1.8% | 836 | 704 | |
| 13481-JF-6 | 13481-JF | M03 | 2,632 | 2,563 | 1.7% | 837 | 665 | |
| 13718-JF-3 | 13718-JF | M00 | 2,071 | 2,029 | 1.6% | 804 | 657 | |
| 13597-JF-6 | 13597-JF | M03 | 1,652 | 1,624 | 1.5% | 828 | 607 | |

The first five carry most of it: the top five M00 libraries hold 60% of
all M00 floor failures. 13719-JF-7 and 13481-JF-8 were already on 03_1's
review list. No M06 library appears.

## Consequences and open decision

- QC loss per cell type is small (<= 6%), so composition is only mildly
  affected; the loss ratio M00:M06 is similar across types.
- Technical depth that differs by timepoint can look like M00 -> M06
  change in pseudobulk DE. **Owner to decide before stage 14:** a
  library-quality covariate (e.g. log median UMI) in longitudinal models,
  and/or a sensitivity analysis without the libraries above.
- Ask the core / jarnagin whether baseline blocks were cut or processed
  differently (block age, order).

## Low-UMI cell types (03_2)

Fail rate vs the 3.0% average (wasikowr's labels, Baseline + Control):
mast 8.5% (2.9x), T 5.8% (2.0x), smooth muscle 4.7%, B 4.6%, eccrine 4.6%
(MT cut part of it), almost all at the 500-UMI floor. Kept for Layer 1;
check per compartment.
