# Cohort definitions

A **cohort** is a membership rule over libraries: which samples a given
Layer 2 analysis uses. Cohorts are cheap — add one the day you need it.
You do not have to know now which cohort the paper will end up using.

Layer 1 (QC, integration, clustering, annotation) uses cohort
`reference`, i.e. all cells. Only Layer 2 subsets.

| Cohort | Purpose |
|---|---|
| `reference` | all cells; used for Layer 1 only |
| `placebo` | **primary Project 1 cohort** — placebo arm, post-escape samples censored |
| `placebo_all` | placebo arm, uncensored; sensitivity analysis against `placebo` |
| `placebo_hc` | placebo arm + healthy controls; for the "normalization toward healthy" anchor |
| `allarms` | both arms; contrast and validation |
| `basectrl` | baseline SSc + healthy controls; cross-sectional disease signature |

Rules of use:

- Cohort selection never changes clusters or labels. Those come from the
  frozen labelset built on `reference`.
- Any cohort that censors post-escape samples needs the clinical table,
  which is not in git. `cohort_samples()` will refuse to run without it.
- When you add a cohort, add a row to the table above.
