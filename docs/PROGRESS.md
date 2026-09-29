# PROGRESS — where this project stands

**The handoff file.** Read this first in any new session, human or AI. It is
the only file that records *temporal* state: what is done, what is in flight,
what is blocked, and what comes next. `AGENTS.md` says how things are
organised; this says where we are.

**Update it at the end of each working session.** Two minutes. Move finished
items to Done, revise In flight, add anything newly blocked. If it goes stale
it becomes worse than nothing, because it will be believed.

Last updated: **2026-09-29** · Updated by: Rezvan Parvizi (with Claude Code)

---

## Current focus

**Project 1** — molecular basis of spontaneous mRSS improvement in the ASSET
placebo arm. Stage 00 is done for freeze01. Next are the clinical table and
the meeting with colleagues about the data issues, then stage 03 (QC).

---

## Done

### Stage 00 — freeze01 (2026-09-29)

- **R0001 `import_raw_libraries`.** 221 libraries (209 SSc from 79
  subjects, 12 HC) in 17 pools. 697 files (12.2 GB) copied into
  `data/raw/`, md5-verified and read-only. Manifest generated.
- **R0002 `build_manifest`.** Structural checks pass. `placebo_all` = 107
  libraries from 39 subjects; `basectrl` = 72.
- **R0003 `batch_vs_design_crosstab`.** Read its `VERDICT.txt`.
  - Pool vs **timepoint: independent** (p = 1.0, V = 0.15). The
    longitudinal design is clean.
  - Pool vs **arm**: p = 0.02, V = 0.34, driven by 13678-DP / 13679-DP.
    Irrelevant within placebo; matters for any arm comparison.
  - Pool vs **group (HC)**: V = 1. The controls are in their own pools
    (CellRanger 9.0.0 and forced cell counts for 5 of them). SSc-vs-HC
    differences are partly technical.
  - 23 of 74 multi-sample subjects (31%) span more than one pool, so pool
    enters within-subject contrasts. Plan: pool as a covariate in
    pseudobulk models. Not yet recorded as a decision.
- Library flags: `13639-JF-11` and `-12` are `timepoint_unresolved` (same
  subject, both labelled M03, no M00). Kept for Layer 1, excluded from all
  cohorts.
- **Sample sheets for building the clinical table** (patient-level):
  `data/patient_level/freeze01/reference/00_manifest_and_batch/import_raw_libraries/`
  `sample_sheet_per_library.csv` and `sample_sheet_per_subject.csv`
- **Discussion list for colleagues:** `docs/data_issues_freeze01.md`
  (controls, pools, 16 unmapped libraries, the timepoint flag, 6 metric
  outliers, the upstream object, code points to pass on).

### Environment (complete, verified)

- R 4.6.1 built from source at `~/R/R-4.6.1` (`--enable-R-shlib`, no X11)
- OpenBLAS 0.3.34 at `~/opt/openblas` (`DYNAMIC_ARCH`, pthread)
- Bioconductor 3.23 — escaped the system R's 3.18
- 292 packages installed; **277 recorded in `renv.lock`**, committed
- CRAN pinned to P3M snapshot `focal/2026-06-01`
- Object-critical pins verified matching the handoff objects:
  Seurat 5.5.0 · SeuratObject 5.4.0 · Matrix 1.7-5 · harmony 2.0.3 ·
  BPCells 0.3.1
- Functional tests passed: Seurat 5 object from an on-disk BPCells matrix with
  PCA; `dream()` fitting `~ timepoint + (1|subject)` on 12 subjects × 2
  timepoints
- Six installation workarounds recorded in `docs/ENVIRONMENT.md`
- `conda config --set auto_activate_base false` — deliberate, see workaround 1

### Repository

- Private GitHub repo `rezvanparvizi/asset-skin-natural-history`, pushed
- Scaffold: three-layer structure, freeze/cohort/labelset axes, provenance
  contract (`init_run()` / `finalize_run()`)
- Analysis stages reordered: metadata → bulk skin → bulk blood → single cell
  (commit `0fb2f17`)
- Bootstrap scripts for server and laptop; `jobs/run.sh` detaching into screen
- `AGENTS.md` + `CLAUDE.md` symlink; `docs/ENVIRONMENT.md`;
  `docs/project_instructions_sections_9-12.md`

### Machines

- U-M MacBook: R 4.6.1 arm64, gfortran 14.2, Xcode CLT licensed, SSH keys to
  GitHub and to `mininubio`, Positron connected via Remote SSH
- Server: SSH key to GitHub, ssh-agent block in `.bashrc` guarded for
  interactive shells only

---

## In flight

Nothing. Clean stopping point.

---

## Blocked, waiting on someone else

| What | Waiting on | Why it matters |
|---|---|---|
| **Meaning of pool and sample-name date** | sequencing core / Jarnagin | batch_id is currently the pool (the best available stand-in). Confirm it is the capture batch. |
| **Site ID per patient** | owner, from the clinical table | Adds site to the batch check (00_2 picks up `site_id` automatically). |
| **Which of 13639-JF-11/-12 is Baseline** | Jarnagin / DCC | Both are excluded from cohorts until resolved. |
| **Chemistry / fixation** confirmation | sequencing core | Probe set (v1.1.0), CellRanger (9.0.1; 9.0.0 for 5 controls) and reference are confirmed. Chemistry and fixation are still TODO in freeze01.yml. |
| **Clinical table** from the U-M DCC | DCC | Needed for cohort censoring (`placebo` drops post-escape samples), and for every model's covariates. |
| **Colleague's label tables** (barcode → lineage/celltype/subtype) | Jarnagin | Preferred over inheriting her Seurat objects; decouples us from her ongoing iteration. Questions listed in `docs/inherited_objects.md`. |
| **Which of her objects is current** | Jarnagin | `res0p2` / `res0p3` / `from3b` / `_v2` / `Round1` / `Round2` cannot be resolved from filenames. |
| **`/hits/home/parvizi` ownership** | wasikowr | Currently owned by her, mode 777. Cannot `chmod` it; working in a `700` subdirectory instead, but the parent stays world-writable so the directory could be deleted wholesale. |
| **Lab Dropbox folder convention** | Jarnagin / lab | Needed before writing the sharing side of `scripts/backup.sh`. |
| **Whether `/hits` is IT-backed** | whoever administers it | Determines whether Dropbox is needed as a third copy. |

---

## Next up, in order

1. **Clinical table.** The owner builds it against the sample sheets and
   puts it in `data/clinical/`. It must include site ID and the
   escape-therapy start month. Then re-run 00_2 to add site.
2. **Colleague meeting** on `docs/data_issues_freeze01.md`. Record the
   answers in `docs/decisions.md`.
3. **Stage 03: QC** per library from `data/raw/*/*/sample_filtered_feature_bc_matrix.h5`.
   The 6 CellRanger outliers are already flagged.
4. **Stage 04: SoupX** from the sample_raw matrices. Write the rho rule
   in decisions.md BEFORE running it (see the 2026-09-29 entry).
5. Decide the doublet route (scDblFinder is blocked).
6. freeze02 (~20 new samples, about a month away): add them with a new
   `freeze02.yml` and re-run everything. Hold the expensive manual
   subtype annotation (stages 07–10) until then.
7. Deferred from before: backfill environment decisions; backup script;
   bulk skin stage 01 once the owner has placed the inputs.

---

## Open questions — decide deliberately, do not default

These are recorded because defaulting on them silently would bias the
analysis.

**Outcome definition.** Global mRSS change versus forearm-site skin score. The
biopsy is forearm; local score tracks local biology better (THBS1 r = 0.76
local versus weaker against global). And binary improver (≥5 points or >20% at
12 months, the definition in papers 01/03) versus continuous trajectory —
binary is comparable to prior work, continuous is better powered with ~19 vs
~19. **Not yet decided.** Ask the DCC whether site-level scores exist.

**Primary timepoint.** Month 6 is contaminated by escape therapy (16/44
placebo, starting *because* they were worsening). Baseline-predicts-12-month
is the cleaner framing. **Not yet decided.**

**Inherit or rebuild Layer 1.** Taking the colleague's annotation saves months
but inherits her filtering history — note the `prefilter_bcell_2026-08-26`
invalidation in her tree, where an annotation change retroactively invalidated
downstream results. Building fresh gives a clean documented reference layer
both could share. **Worth agreeing a single versioned canonical label set with
her before either of us has a figure.**

**How rho is set for SoupX.** The colleague tuned 0.1–0.25 by eye per
sample. Write a rule first (marker-leakage criteria, the same for every
library, never judged with improver status in view). **Not yet decided.**

**Doublet detection route.** `scDblFinder` is blocked by the gcc 9 / C++20
ceiling. Options: `DoubletFinder`, `scds`, or inherit her existing calls (she
has `scDblFinder 1.16.0` results on this dataset). **Not yet decided.**

---

## Deferred TODOs

- macOS Keychain for SSH passphrases (`ssh-add --apple-use-keychain`) — does
  not need an Apple ID; the local login keychain suffices
- Homebrew + HDF5 on the laptop, required before BPCells works locally
- `MuSiC` / `TOAST` — `TOAST` no longer downloads from Bioc 3.23. Needed only
  at stage 14. Options: Bioc archive, vendor from GitHub, Scaden, CIBERSORTx
- Request Apptainer from whoever administers the servers — it is the cleanest
  path to identical environments across the group's 6-7 boxes
- Fill in `config/freezes/freeze01.yml` TODOs once the core confirms them
- Project 2 (QSP model) — separate repository `asset-qsp-skin`, not started,
  deliberately not competing with Project 1

---

## Session log

Newest last. One or two lines per session: what was attempted, what landed.
Detail belongs in `docs/decisions.md` (judgment) and `docs/runs.csv` (runs).

```
2026-09-22  Reviewed the 7 prior ASSET papers; drafted project framing.
            Surveyed the colleague's server directory (1,600-line ls -R):
            hierarchical subclustering pipeline, BPCells on-disk objects,
            MuSiC/Scaden deconvolution, date-stamped output duplication.
2026-09-23  Designed the repository: three layers, freeze/cohort/labelset,
            provenance contract. Scaffold created.
2026-09-24  Laptop setup: R 4.6.1 arm64, gfortran, Xcode license. SSH keys to
            GitHub and mininubio. Repo created and pushed. Surveyed the server
            environment: system R 4.3.2 / Bioc 3.18, no OpenBLAS, gcc 9,
            no containers, no P3M focal binaries. Chose home-built R + renv.
            OpenBLAS 0.3.34 built. R build FAILED — conda ICU 73 vs system 66.
2026-09-25  R 4.6.1 built successfully in a conda-free shell. Bootstrap
            packages installed. BPCells and variancePartition smoke tests
            passed (dream() fits a repeated-measures design).
            conda auto_activate_base set to false.
2026-09-28  renv install across several passes. Workarounds: USE_BUNDLED_LIBUV
            for fs; /sbin on PATH for curl's static-libcurl fallback; presto
            is GitHub-only not CRAN. scDblFinder abandoned (scrapper needs
            C++20; conda-forge gcc route also needs gfortran_linux-64).
            renv.lock initially recorded 1 package — snapshot.type had to be
            "all", not "explicit". 277 packages now locked and committed.
            Analysis stages reordered to metadata -> bulk -> single cell.
            docs/ENVIRONMENT.md rewritten with all six workarounds.
2026-09-29  Surveyed colleague's ASSET_Flex tree (read-only). Stage 00
            built: raw copied + md5 (R0001), manifest (R0002), batch
            check (R0003: timepoint clean, HC fully confounded with
            pool). Patient-level outputs moved to data/patient_level/.
            Discussion list for colleagues written.
```
