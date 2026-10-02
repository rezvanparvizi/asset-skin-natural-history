# Ambient RNA (stage 04), freeze01 — re-assessment 2026-10-02

Run of record: R0020 (`04_1_soupx_per_library.R`). Rule: decisions.md
2026-10-01. Aggregate numbers only.

## How cell types were obtained in 04_1

No external labels, no integration: 04 runs before 05. For each library
separately: QC-pass cells -> LogNormalize -> 2,000 HVG -> 30 PCs -> Louvain
0.8 (same recipe everywhere). SoupX uses only these cluster IDs
(`setClusters`). For the checks, each cluster got a provisional
`coarse_lineage` = the panel with the highest mean normalised expression
(Keratinocyte, Fibroblast, Immune, Plasma, Endothelial, Mural, Melanocyte,
Gland; genes in the script). wasikowr's labels were not used.

Known weakness (found 2026-10-02, integration_freeze01.md 2f): this
label is unreliable for immune and gland cells — keratin ambient can win
the argmax in small libraries. The leakage checks are therefore computed on
noisy groups (a few keratinocytes or doublets in "Immune" inflate keratin in
immune cells; immune cells labelled Keratinocyte are missing from it).

## Our rho vs a marker-based estimate (from R0020 tables)

SoupX's manual logic: in cells that cannot express a gene set, its share of
UMIs = rho x its share of the soup. Soup share (median): KRT1/5/10/14 4.2%,
COL1A1/1A2/3A1 1.3%; all 4 keratins are in the top-20 soup genes in all 218
libraries.

| estimate | p10 | p25 | median | p75 | p90 |
|---|---|---|---|---|---|
| autoEstCont (used) | 0.026 | 0.038 | **0.054** | 0.073 | 0.090 |
| keratin in immune-labelled cells | 0.046 | 0.074 | **0.121** | 0.174 | 0.226 |
| keratin in fibroblast-labelled cells | 0.062 | 0.087 | **0.127** | 0.195 | 0.276 |
| collagen in immune-labelled cells | 0.062 | 0.107 | **0.146** | 0.233 | 0.348 |

Spearman(autoEstCont, keratin-in-fibroblast estimate) = 0.59: the ranking
of libraries agrees, the level does not. The marker estimates are likely
upper-biased (noisy labels, doublets), so the true value is probably
between 0.05 and ~0.13. autoEstCont's default prior is centred at 0.05
(SoupX default as I recall it; to verify), which may be pulling it down.
**Reading: our correction is probably too weak, by up to ~2x.** Consistent
with leakage only halving (keratin in immune 0.52% -> 0.29% of UMIs).

## What the lab used

- No SoupX script exists in /hits/home/wasikowr/ffpe/Asset; no rho is
  stored in seurat.RDS or seurat_amb.RDS (both hold the same corrected,
  non-integer counts). Batch*/ffpe_report.pdf not read (no PDF tools on
  the server) — may document it.
- **Owner asked wasikowr (2026-10-02): her pipeline uses rho = 0.1.**
- Her counts vs raw h5 for the same cells (jarnagin's paired objects,
  181,766 cells; ASSET_base_control.RDS vs
  ASSET_base_control_integer_BPCells.RDS): 29.6% fewer UMIs per library
  (median; 80% of libraries 0.29-0.30). **This is NOT an ambient
  correction** (per-gene check, 2026-10-02): soup genes lose LESS than
  average (KRT5/KRT14 kept 0.75, KRT1/10 0.77, COL1A1 0.82, vs 0.705
  overall); genes in abundance deciles 1-6 all lose ~31-32%; the rarest
  decile loses 7%. SoupX would remove most from the soup genes. So a rho
  0.1 correction is compatible, on top of an unexplained, near-uniform
  ~25-30% reduction from some other step (different CellRanger run /
  settings? a transformation?). **An earlier note in this file read the
  30% as "fixed rho ~0.3" — that was wrong and is withdrawn.**
- The note "tuned 0.1-0.25 by eye" (decisions.md 2026-09-29, 2026-10-01)
  has no traceable source in her folder; superseded by her answer (0.1).
- jarnagin did not run SoupX; she used wasikowr's corrected counts (immune,
  4a) or raw h5 counts (1b, TCAT).

## Options (owner to decide; not yet changed)

1. **Marker-based rho per library** (SoupX `estimateNonExpressingCells` +
   `calculateContaminationFraction`) with fixed gene sets (keratins,
   collagens, HB) on reliable clusters — the method SoupX recommends when
   the automatic estimate is doubtful. Needs better clusters than the
   per-library argmax: the 05_3 broad lineages (projected back per library)
   or wasikowr's broad labels (cover only ~182k cells).
2. autoEstCont with a weaker/higher prior — still opaque; not preferred.
3. Cross-check with DecontX (celda; per-cell contamination from cluster
   labels, no empty droplets needed) on a subset of libraries.
   CellBender is not proposed: GPU-oriented and Flex support unverified.
Any change: written rule first, same checks (a)-(d), labels not outcome,
then 04 -> 05_1 -> 05_2 re-run (~4-5 h of jobs).

## Decision (owner, 2026-10-02)

Fixed rho for every library; no marker-based estimation (no reliable
clusters yet). Value 0.13 — see decisions.md 2026-10-02. Implemented as a
separate step, `04_2_correct_ambient_fixed_rho.R`, reusing 04_1's
clusters and soup; output in data/bpcells/freeze01/counts_corrected/.
