# Gene sets

Plain-text (one gene per line) or CSV files. Keep the source of every
set in the table below, because a module score is only interpretable if
you know where the genes came from.

| File | Source | Used for |
|---|---|---|
| `cd28_costimulation_reactome.txt` | Reactome, via Mehta 2022 (paper 03) core enrichment genes | anchoring single-cell findings to the published bulk result |
| `sasp_senescence.txt` | TODO — pick and cite one | senescence/SASP decline during regression |
| `tgfb_response_thbs1_comp.txt` | TODO — Sargent 2010 TGF-beta responsive signature | THBS1/COMP/CTGF fall as skin improves |
| `poon_tcell_subtypes.csv` | Poon et al., via colleague's `reference/` | T-cell subtype labeling |
| `ma2024_fibroblast_subtypes.csv` | Ma 2024, via colleague's `reference/` | fibroblast subtype labeling |
| `intrinsic_subset_centroids.csv` | Milano 2008 / Franks 2019 | mapping single-cell pseudobulk onto intrinsic subsets |

Label every set in `docs/decisions.md` as one of:
(a) established in human SSc skin,
(b) established in liver/lung fibrosis or mouse models and imported as
    hypothesis,
(c) derived from this dataset.

Most published fibrosis-regression gene sets are (b). Never let (b) be
reported as (a).
