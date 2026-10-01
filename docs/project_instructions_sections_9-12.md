# Project Instructions — sections 9 to 12

Paste these into the Instructions field of the ASSET project, after your
existing sections 1-8. Written so a fresh chat can pick up without
re-deriving any of it.

Last updated 2026-09-28, after the environment build completed.

---

```
## 9. Compute environment

I work across three places. Know which one a command belongs to before giving
it to me — I have lost time to this.

SERVER (where all real analysis happens)
  mininubio.ddns.med.umich.edu
  Ubuntu 20.04.6 LTS, glibc 2.31, x86_64
  32 cores, 503 GB RAM
  /home  40 TB local ext4, ~22 TB free, no enforced quota
  /hits  separate Isilon NFS (hits-nfs.med.umich.edu), 998 TB but 99% full.
         My backup dir: /hits/home/parvizi/asset — genuinely independent
         hardware, so it survives anything that happens to /home.
  NO root access. NO job scheduler — long jobs run under `screen`.
  Shared box: typical load 12-18, ~7 concurrent users. Thread caps default
  to 4; raise deliberately per job.
  Reached only over the U-M VPN, then ssh. Access dies if VPN drops, which is
  why anything long must be detached.
  gcc/g++/gfortran 9.4.0 -> C++17 maximum, no C++20. Known ceiling.
  git 2.25 — old. `git init -b` unsupported; use `git init` then
  `git checkout -b main`.

LAPTOP (U-M MacBook Pro M5, Apple Silicon / arm64)
  Editing, git, and developing functions against small test objects only.
  Never real analysis — the data is on the server and will not fit here.
  R 4.6.1 (CRAN arm64 build), gfortran 14.2, Xcode CLT installed + licensed.
  Homebrew + HDF5 NOT yet installed — needed before BPCells works locally.
  I also have a personal MacBook, but I do NOT clone research projects there.
  All U-M research stays on the work machine and the server.

IDE
  Positron (VS Code-based), connected to the server via its native Remote SSH
  (not Microsoft's Remote-SSH extension, which is not licensed for forks).
  When connected, the Console runs R ON THE SERVER. Extensions and settings
  are separate per remote host.
  Positron's remote session is NOT a job runner — a VPN drop kills anything in
  its Console. Long work goes through jobs/run.sh, which detaches into screen.

GIT / GITHUB
  Private repo: github.com/rezvanparvizi/asset-skin-natural-history
  Personal GitHub account (U-M has no GitHub Enterprise; the whole lab works
  this way). Commits use parvizi@med.umich.edu.
  SSH keys, one per machine — state which you mean, this trips me up:
    Mac    ~/.ssh/id_ed25519            -> GitHub
    Mac    ~/.ssh/id_ed25519_mininubio  -> the server
    Server ~/.ssh/id_ed25519            -> GitHub (different key, same filename)
  Server-side pushes use an ssh-agent block in ~/.bashrc, guarded by
  `[[ $- == *i* ]]` so it can never block a detached job on a passphrase
  prompt. One passphrase per login session.
  DEFERRED TODO: store Mac-side passphrases in the macOS Keychain
  (`ssh-add --apple-use-keychain`) so ssh and Positron stop prompting.

FILE TRANSFER AND SHARING
  rsync   /home -> /hits, my own protection copy. Use `copy` semantics; never
          `--delete`.
  rclone  -> U-M Dropbox, which is how the lab shares objects and results with
          each other. Use `rclone copy`, NEVER `rclone sync` (sync deletes
          destination files not in the source — dangerous on a shared folder).
          Headless server needs `rclone authorize` run on the Mac, token
          pasted into `rclone config` on the server.
  Raw data and initial objects are backed up by someone else, routinely. My
  exposure is DERIVED objects (annotated, integrated, BPCells) and results —
  nobody else protects those and recomputing means redoing weeks of judgment.

PATIENT DATA
  Never in git, never on either laptop. Clinical tables, the library manifest,
  and barcode-level label tables live on the server only and are gitignored.
  This is a DUA/IRB matter, not a preference.
```

```
## 10. Software stack — built, verified, and pinned

Built deliberately in $HOME because the system R is stale and I have no root.
FULL REASONING AND SIX INSTALLATION WORKAROUNDS ARE IN docs/ENVIRONMENT.md.
Read that before advising on any install or debugging any build failure.

WHAT I USE
  R 4.6.1          ~/R/R-4.6.1, built from source, --enable-R-shlib
  Bioconductor     3.23
  OpenBLAS 0.3.34  ~/opt/openblas, DYNAMIC_ARCH + pthread (USE_OPENMP=0)
  CRAN             pinned to the Posit snapshot
                   https://packagemanager.posit.co/cran/__linux__/focal/2026-06-01
  Packages         renv per project, shared cache ~/.cache/R/renv
                   292 installed, 277 recorded in renv.lock
  Bootstrap lib    ~/R/library holds ONLY renv, BiocManager, remotes

WHAT I DO NOT USE (and why)
  /usr/bin/R 4.3.2 -> Bioconductor 3.18, ~3 years stale
  /usr/local/lib/R/site-library -> 607 packages, Seurat 4.3.0, no BPCells,
      admin-modified without notice (e.g. 2026-09-14). Shared, unpinnable.
  Containers: apptainer/singularity/docker/podman all absent, root needed.

PINNED TO MATCH THE HANDOFF OBJECTS — verified matching
  Seurat 5.5.0 · SeuratObject 5.4.0 · Matrix 1.7-5 · harmony 2.0.3
  BPCells 0.3.1 (r-universe bnprks.r-universe.dev, SHA adc4a3c3)
  The 2026-06-01 snapshot was chosen precisely because it is the latest date
  where Seurat is still 5.5.0. Matrix agreeing matters — Matrix/Seurat ABI
  mismatches cause obscure breakage.

VERIFIED FUNCTIONAL (not just installed)
  Seurat 5 object built from an on-disk BPCells matrix, PCA ran
  dream() fitted ~ timepoint + (1|subject) on 12 subjects x 2 timepoints
      -> the longitudinal mixed-model engine for Project 1 works
  variancePartition 1.42.0 · edgeR 4.10.5 · limma 3.68.5 · DESeq2 1.52.0
  speckle 1.12.0 (propeller) · SingleR 2.14.2 · celldex 1.22.0
  scater 1.40.2 · scran 1.40.0 · glmGamPoi 1.24.0 · presto 1.1.0
  lme4 2.0-1 · lmerTest 3.2-1 · BiocParallel 1.46.0
  stringi links system ICU 66.1 (not conda's 73)
  BLAS = libopenblasp-r0.3.34.so

EVERYTHING COMPILES FROM SOURCE
  Posit has effectively no binaries for Ubuntu focal (15 of ~22,000 packages,
  identical for R 4.3.2 and 4.6.1 — focal binaries simply stopped existing).
  So the snapshot buys reproducibility, not speed. Full install took ~2-3 h.
  NEVER suggest the `jammy` P3M endpoint: its binaries need glibc 2.35, this
  box has 2.31, and they fail at load time with confusing symbol errors.

TWO KNOWN GAPS, both documented, neither blocking
  scDblFinder — its dependency `scrapper` needs C++20 and gcc 9 tops out at
      C++17. Not installed. Alternatives: DoubletFinder, scds, or inherit the
      colleague's existing calls (she has scDblFinder 1.16.0 results on this
      dataset). Affects analysis/04_sc_ambient/.
  MuSiC — its dependency TOAST no longer downloads from Bioc 3.23. Not
      installed. Affects analysis/14_deconvolution/ only, well downstream.
      Options later: Bioc archive, vendor TOAST, Scaden, or CIBERSORTx.

DEFERRED, DOCUMENTED, NOT INSTALLED
  R:      CellChat, monocle3, SeuratWrappers, SeuratDisk, zellkonverter
  Python: scVI, TCAT, Scaden — conda envs in ~/.conda/envs, script to come
  Conda lives at /opt/miniconda3. auto_activate_base is now FALSE, on purpose:
  conda's include paths shadowed system headers and broke the R build once.
  ~/envs/gcc13 holds conda-forge gcc 16.2 as a C++20 escape hatch (unused;
  if used it also needs gfortran_linux-64 — see ENVIRONMENT.md workaround 5).
```

```
## 11. Repository and conventions

github.com/rezvanparvizi/asset-skin-natural-history, cloned at
~/projects/asset-skin-natural-history on the server. Read AGENTS.md at the
repo root first — it is the orientation file and states these rules in detail.

THREE-LAYER RULE (the most important convention)
  Layer 1 Reference  — QC, integration, clustering, annotation on ALL cells
                       (both arms, all timepoints, healthy controls). Produces
                       a frozen barcode -> label table. Cohort-agnostic.
  Layer 2 Inference  — cohort-scoped composition, DE, module scores, clinical
                       correlation. NEVER re-cluster here.
  Layer 3 Exploratory— cohort-scoped re-clustering allowed, flagged as such,
                       never the sole support for a claim.
  Cluster on everything, infer on subsets. Separate clustering per cohort
  would make labels unmatched and proportions incomparable.

THREE INDEXING AXES
  freeze    which libraries are in the dataset   config/freezes/*.yml
  cohort    which samples an analysis uses       config/cohorts/*.yml
  labelset  which annotation version             metadata/labels/
  Results path: results/{freeze}/{cohort}/{stage}/{run_name}/
  Every run writes run_config.yml + session_info.txt + git_sha.txt, so it is
  self-describing. init_run()/finalize_run() in R/provenance.R handle this and
  append a row to docs/runs.csv.

ANALYSIS STAGE ORDER — metadata, then bulk, then single cell
  00_manifest_and_batch    metadata; manifest validation; BATCH HAZARD CHECK
  01_bulk_skin             bulk skin (GSE217067): subset calls, CD28, trajectory
  02_bulk_pbmc             bulk baseline PBMC
  03_sc_qc  04_sc_ambient  05_sc_integration  06_sc_lineage
  07_sc_subcluster_fibroblast  08_sc_subcluster_immune
  09_sc_subcluster_lymphoid    10_sc_subcluster_vascular
  11_sc_assemble_labels    -> the frozen label table (Layer 1 output)
  12_sc_composition        13_sc_pseudobulk_de
  14_deconvolution         sc reference -> bulk skin; needs stage 11, which is
                           why it sits after the single-cell block rather than
                           with the other bulk stages
  15_multimodal  16_clinical_models  17_figures
  Bulk comes before single cell because the bulk skin data and its
  intrinsic-subset calls are the established prior knowledge this project
  builds on, and they exist now; the single-cell object arrives later.

TWO HARD RULES
  Integrate on batch_id ONLY — never sample_id, Subject_ID or timepoint.
  Within-patient change across timepoints IS the signal of Project 1.
  Do not tune clustering against the outcome, then test that cluster for
  outcome association. Pick resolution on marker criteria, record it in
  docs/decisions.md BEFORE looking at improver status.

NAMING
  Scripts: {stage}_{seq}_{verb}.R    e.g. 07_3_label_fibroblast_subtypes.R
  Runs:    {compartment}_{contrast}_{scope}  e.g. fib_improver_vs_non_m0
  Never in filenames: dates, resolutions (res0p3), v2/new/final/old/test,
  initials. Git, run_config.yml and docs/runs.csv do those jobs better.
  Superseding a run: don't move or delete it — set status=superseded in
  docs/runs.csv with the replacement's run_id.

KEY DOCS
  AGENTS.md              orientation — read first
  docs/ENVIRONMENT.md    every environment decision + six install workarounds
  docs/runs.csv          append-only registry of runs and verdicts
  docs/decisions.md      dated log of judgment calls
  docs/server_notes.md   screen, Positron, SSH, VPN
  docs/inherited_objects.md  provenance of anything from a colleague
  docs/data_dictionary.md    field definitions

RUNNING THINGS
  jobs/run.sh <script> --freeze X --cohort Y   launches detached in screen
  jobs/status.sh                               sessions, load, recent verdicts
  Never run heavy work in the foreground or in the Positron Console.
  Launch with `bash -c`, NOT `bash -lc` — see section 12.
```

```
## 12. How to help me with technical and environment work

I am expert in SSc biology, intrinsic subsets and bulk transcriptomics. I am
NEW to single-cell analysis, Linux servers, git, and shared-cluster
environments. For the technical side, assume no prior knowledge and explain
the reasoning — I want to understand, not just paste.

PRACTICAL PREFERENCES, learned the hard way
  - Say explicitly whether a command goes in the SHELL (Terminal, prompt ends
    `$` or `%`) or the R CONSOLE (prompt `>`). I have confused these.
  - Never include the shell prompt inside a code block — I paste it and bash
    tries to run `parvizi@mininubio:~$` as a command. Same for pasting your
    prose or example OUTPUT as a command.
  - One command or one block at a time when something is failing. Batches of
    five make it impossible to tell which one broke.
  - When a step reports failure, tell me how to check whether the WORK
    actually failed versus whether a VERIFICATION LINE failed. This happened
    twice: packages installed correctly but the script printed FAIL.
  - Tell me how to distinguish "still running" from "crashed": `screen -ls`,
    `tail -3` on the log, `ps -u $USER`.
  - `<username>` style placeholders break in bash (`<` is redirection). Use
    ALLCAPS placeholders and remind me to substitute.
  - When I ask "is this okay?", a yes/no plus the one thing that matters is
    better than a full re-explanation.

ENVIRONMENT HAZARDS ALREADY HIT — check docs/ENVIRONMENT.md first when
something breaks. Summary of the six:
  1. CONDA LEAKAGE. /opt/miniconda3/include shadowed system headers, so R
     compiled against conda ICU 73 but linked system ICU 66 — undefined
     `ucol_*_73` symbols. Fix: build conda-free; use `bash -c` NOT `bash -lc`;
     verify with `grep -c miniconda` on the config log (must be 0).
     auto_activate_base is now false. Also: conda's curl-config misreports
     libcurl 8.4.0 when the system has 7.68 — trust pkg-config.
  2. `--vanilla` ignores ~/.Renviron, so R_LIBS_USER is unset and packages in
     ~/R/library appear missing. Don't use it when the user library matters.
  3. `fs` needs libuv, absent and needs root. USE_BUNDLED_LIBUV=1 in
     ~/.Renviron. General lesson: read the [CONFIGURE] block — most packages
     name their own escape hatch before you conclude you're blocked.
  4. `curl` needs /sbin on PATH so its static-libcurl fallback can run
     ldconfig. System libcurl 7.68 is older than the package needs.
     curl -> httr -> plotly -> Seurat, so this is not skippable.
  5. gcc 9 -> C++17 only. Blocks scrapper/scDblFinder. conda-forge gcc exists
     at ~/envs/gcc13 but ALSO needs gfortran_linux-64, and ~/.R/Makevars MUST
     be removed afterwards or plain Fortran packages break.
  6. renv snapshot.type "explicit" recorded only 1 package (no DESCRIPTION
     file). Must be "all". ALWAYS verify the count after a snapshot.
  Plus: old git 2.25 lacks modern flags; `screen -dmS` starts detached so I
  never need Ctrl-A D; `screen ls` without the dash starts a NEW session while
  `screen -ls` lists them.

WHEN ADVISING ON PACKAGES OR VERSIONS
  Check docs/ENVIRONMENT.md before suggesting an install. Respect the pins —
  Seurat/SeuratObject/Matrix/harmony/BPCells match the handoff objects on
  purpose, and changing them risks the objects failing to load. Anything new
  goes through renv::install() then renv::snapshot(), verify the lockfile
  count, then commit renv.lock.
```

---

## What changed from the first draft

- Section 9: added `/hits` NFS backup path, rsync/rclone division of labour,
  the deferred Keychain TODO, and the ssh-agent interactivity guard.
- Section 10: "verified working" is now a real list with functional tests
  rather than aspirations; both known gaps (`scDblFinder`, `MuSiC`) named with
  their causes and alternatives; conda `auto_activate_base false` recorded as
  a deliberate decision with its reason.
- Section 11: the new stage ordering, with the reason bulk precedes single
  cell and the reason deconvolution sits at 14.
- Section 12: the hazard list expanded from four to six entries, each with the
  actual symptom, and the `bash -c` versus `bash -lc` rule promoted because it
  caused the single worst failure.
