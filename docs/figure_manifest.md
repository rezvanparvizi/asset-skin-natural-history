# Figure manifest

Figure number -> the script and run that produced it. Fill in as
figures stabilize. The point is that any figure in the manuscript can
be regenerated without guesswork.

| Figure | File | Script | Run ID | Cohort | Caption gist |
|---|---|---|---|---|---|
| 1a | `fig1a_cohort_schematic.pdf` | `17_1_cohort_schematic.R` | | placebo | trial design, biopsy timepoints, escape-therapy censoring |
| 1b | `fig1b_umap_all_lineages.pdf` | `17_2_overview_umaps.R` | | reference | all cells, lineage labels |
| 2a | | | | placebo | composition by improver status, baseline |
| 2b | | | | placebo | composition trajectory M0 -> M3 -> M6 |
| 2c | | | | placebo_hc | distance to healthy-control centroid over time |
| 3a | | | | placebo | fibroblast subtype proportions, improver vs non |
| 3b | | | | placebo | adipocyte-lineage recovery vs loss (atrophy test) |
| 4a | | | | placebo | pseudobulk DE, fibroblast, improver vs non |
| 4b | | | | placebo | CD28 module score, cell-type resolved |
| 5a | | | | placebo | bulk deconvolution vs single-cell proportions |
| S1 | | | | reference | batch vs design crosstab (hazard #1 evidence) |
| S2 | | | | placebo_all | escape-therapy sensitivity analysis |
| S3 | | | | placebo | freeze01 vs freeze02 robustness |

Naming: `fig{N}{panel}_{description}.pdf`, e.g.
`fig2c_fibroblast_composition_by_improver.pdf`.

Every figure script reads from `results/` and writes to `figures/`. No
figure script should compute a statistic — if a number appears in a
panel, it came from a run in `docs/runs.csv`.
