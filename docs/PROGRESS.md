# PROGRESS — where this project stands

**The handoff file.** Read this first in any new session, human or AI. It is
the only file that records *temporal* state: what is done, what is in flight,
what is blocked, and what comes next. `AGENTS.md` says how things are
organised; this says where we are.

**Update it at the end of each working session.** Two minutes. Move finished
items to Done, revise In flight, add anything newly blocked. If it goes stale
it becomes worse than nothing, because it will be believed.

Last updated: **2026-09-28** · Updated by: Rezvan Parvizi

---

## Current focus

**Project 1** — molecular basis of spontaneous mRSS improvement in the ASSET
placebo arm. Infrastructure is complete; the next work is data-facing.

Nothing is in flight right now. The environment build finished 2026-09-28 and
the repository is committed and pushed.

---

## Done

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
| **Library manifest** with `batch_id` | whoever generated the libraries | Hazard #1. `analysis/00_manifest_and_batch/` cannot run without it, and no longitudinal result is trustworthy until batch is checked against timepoint and arm. |
| **Platform / chemistry / probe set** confirmation | sequencing core | If it is probe-based 10x Flex, the panel is not whole-transcriptome — that constrains which genes can be asked about at all. `config/freezes/freeze01.yml` has these as `TODO`. |
| **Clinical table** from the U-M DCC | DCC | Needed for cohort censoring (`placebo` drops post-escape samples), and for every model's covariates. |
| **Colleague's label tables** (barcode → lineage/celltype/subtype) | Jarnagin | Preferred over inheriting her Seurat objects; decouples us from her ongoing iteration. Questions listed in `docs/inherited_objects.md`. |
| **Which of her objects is current** | Jarnagin | `res0p2` / `res0p3` / `from3b` / `_v2` / `Round1` / `Round2` cannot be resolved from filenames. |
| **`/hits/home/parvizi` ownership** | wasikowr | Currently owned by her, mode 777. Cannot `chmod` it; working in a `700` subdirectory instead, but the parent stays world-writable so the directory could be deleted wholesale. |
| **Lab Dropbox folder convention** | Jarnagin / lab | Needed before writing the sharing side of `scripts/backup.sh`. |
| **Whether `/hits` is IT-backed** | whoever administers it | Determines whether Dropbox is needed as a third copy. |

---

## Next up, in order

1. **Storage setup.** `./setup.sh /home/parvizi/asset-data`, then
   `Rscript R/paths_check.R`. Creates the `data/`, `results/`, `figures/`,
   `logs/` symlinks and locks `data/clinical/` to mode 700.
2. **Backfill `docs/decisions.md`** with the environment decisions. The
   reasoning currently lives in `ENVIRONMENT.md` and in chat history, not in
   the dated decision log.
3. **Bulk skin — stage 01.** The only modality available today. GSE217067
   with the intrinsic-subset calls and CD28 module scores from Mehta 2022.
   Reproduce the subset assignments, then the placebo-arm trajectory. This is
   also the material the owner knows best, so it is the natural first real
   analysis.
4. **Backup script.** `rsync` tiers to `/hits/home/parvizi/asset`. The Dropbox
   leg waits on the folder convention.
5. **Manifest and batch check — stage 00.** The moment the manifest arrives,
   this runs before any biology. Read its `VERDICT.txt`.
6. **Python envs.** `scripts/bootstrap_python.sh` for scVI, TCAT, Scaden.
   Deferred until the R side has been used in anger.

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
```
