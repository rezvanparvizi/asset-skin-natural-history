# Computing environment — what, where, and why

Written so someone arriving cold — a collaborator, a future you, an AI agent —
understands this setup without reverse-engineering it from the filesystem.
Every non-obvious choice has its reason recorded next to it, and every
workaround discovered during installation is written down.

Built and verified 2026-09-24 → 2026-09-28 on `mininubio.ddns.med.umich.edu`.

---

## Contents

- [Summary](#summary)
- [The machines](#the-machines)
- [Decisions and reasoning](#decisions-and-reasoning)
- [Installation workarounds — READ THIS FIRST when something fails](#installation-workarounds--read-this-first-when-something-fails)
- [What is installed](#what-is-installed)
- [Known gaps](#known-gaps)
- [Reproducing this environment](#reproducing-this-environment)
- [Routine operations](#routine-operations)
- [What was surveyed before deciding](#what-was-surveyed-before-deciding)

---

## Summary

Three layers, each pinning a different thing:

| Layer | Pinned by | Where |
|---|---|---|
| R itself | source build | `~/R/R-4.6.1` |
| BLAS / system libs | OpenBLAS build + the OS | `~/opt/openblas`, `/usr/lib` |
| R packages | **renv**, per project | `<project>/renv/library` |

The system R (4.3.2) and the shared site-library are **not used**.

**Verified working:** R 4.6.1 · Bioconductor 3.23 · OpenBLAS 0.3.34 ·
Seurat 5.5.0 + BPCells 0.3.1 on-disk matrices · `dream()` mixed models ·
292 packages installed, 277 recorded in `renv.lock`.

---

## The machines

### Server: mininubio

```
Ubuntu 20.04.6 LTS   glibc 2.31   x86_64
32 cores, 503 GB RAM
/home  on /dev/md0 (local ext4), 40 TB, ~22 TB free, no enforced quota
/hits  on hits-nfs.med.umich.edu:/ifs/nas11/hits-nfs/MM-DERM-TSOI-Lab
       998 TB, 99% used — a separate Isilon NAS, genuinely independent storage
gcc / g++ / gfortran 9.4.0     cmake 3.31.8     git 2.25
No root access. No job scheduler — jobs run under screen.
Typical load 12-18, ~7 concurrent users. Reached only over the U-M VPN.
```

Pre-existing R installations, none of which we use:

- `/usr/bin/R` — 4.3.2 (Oct 2023), Bioconductor **3.18**
- `/usr/local/lib/R/site-library` — 607 packages, Seurat **4.3.0**, no
  BPCells, admin-modified without notice (e.g. 2026-09-14)
- `/opt/miniconda3/envs/scRNAseqCelltypeAnnotation` — root-owned, contains R
- A colleague's `~/R/R-4.6.0` + `~/R/library` (636 packages) — the precedent
  this setup follows

### Laptop: U-M MacBook Pro M5 (arm64)

```
R 4.6.1 from the CRAN arm64 pkg      gfortran 14.2 from mac.r-project.org
Xcode CLT installed, license accepted (needed `sudo xcodebuild -license accept`)
Homebrew + HDF5: NOT YET INSTALLED — required before BPCells works locally
```

Editing, git, and developing functions against small test objects. Never real
analysis.

**Architecture note:** laptop is arm64, server is x86_64. `renv.lock` is
portable between them; `renv/library/` is **not** — compiled packages can never
be copied across. Parity means *same versions, independently built*.

**Policy: keep both machines on R 4.6.x and upgrade them together,
deliberately.** If macOS R auto-updates to 4.7 while the server stays on 4.6,
`renv.lock` stops restoring and parity is silently gone.

---

## Decisions and reasoning

### R built from source in `$HOME`

System R is 4.3.2, which pins Bioconductor to **3.18 (October 2023)** — nearly
three years stale, with no override possible.

Rejected alternatives:

- **System R + a user library.** Bioc stays 3.18, and the shared site-library
  holds Seurat 4.3.0 / SeuratObject 4.1.3 while we need Seurat 5. The search
  path would carry two incompatible major versions plus 607 packages compiled
  against the old one. Mixed-ABI stacks fail in ways that are miserable to
  diagnose.
- **Asking an admin to upgrade.** Asked; not happening on a useful timescale.
- **Conda R.** Works, proven on this box by a colleague, but compiles from
  source anyway, ships duplicate system libraries, and — as it turned out —
  actively caused the ICU failure below. Kept only for Python tools.
- **Containers.** Technically best for cross-server reproducibility. Rejected:
  `apptainer`, `singularity`, `docker`, `podman` all absent, root needed.
  **Worth requesting** — Apptainer is the standard academic-HPC ask.

### R 4.6.1 specifically

Matches the laptop. R keys its library directory to `major.minor`, so 4.6.0 and
4.6.1 share `library/4.6` and are package-compatible — the colleague's
4.6.0-built objects load fine.

### PATH is *prepended*, unlike the local precedent

The colleague's `.bashrc` appends her R:

```bash
path_append "$HOME/R/R-4.6.0/bin"   # comment: "NOT prepended — conda env R must win"
```

Appending leaves `/usr/bin` earlier in the path, so a bare `R` runs the system
**4.3.2**; her 4.6.0 only runs when invoked by full path. Her stated goal —
letting an active conda env take precedence — is achieved anyway by
prepending, because `conda activate` prepends *later* in the shell's lifetime
and therefore still wins.

So we prepend. Verify with `which R`.

### OpenBLAS built in `$HOME`

The system has only the **reference** BLAS — no OpenBLAS, ATLAS or MKL
anywhere, which is why the colleague's R shows `BLAS_LIBS = -lblas`.

Honest scope: most of the sparse single-cell path (irlba PCA on sparse
matrices, neighbour graphs, UMAP, marker tests) is **not** dense-BLAS-bound, so
this is not a blanket speedup. It matters for Harmony (dense operations on the
embedding, repeated k-means), `ScaleData`, and dense correlation matrices —
steps run many times over a project.

Build flags and why:

- `DYNAMIC_ARCH=1` — one binary selecting the right kernel at runtime, so this
  install is portable to the group's other servers.
- `USE_OPENMP=0` — pthread threading. An OpenMP-threaded BLAS fights with R
  packages that use OpenMP themselves (`data.table`, `glmGamPoi`), causing
  thread oversubscription. Capped at runtime instead.
- `NUM_THREADS=64` — compile-time ceiling, not the runtime default.

R links with `-Wl,-rpath,$HOME/opt/openblas/lib`, so it finds the library
without `LD_LIBRARY_PATH`.

Confirmed: `sessionInfo()$BLAS` →
`/home/parvizi/opt/openblas/lib/libopenblasp-r0.3.34.so` (the `p` confirms the
pthread build).

### CRAN pinned to the P3M snapshot `focal/2026-06-01`

Chosen because it is the **latest date at which Seurat is still 5.5.0** — the
version the handoff objects were written with. Verified:

```
2026-05-01  Seurat 5.5.0
2026-06-01  Seurat 5.5.0   <- chosen
2026-07-01  Seurat 5.5.1
```

At that snapshot all four object-critical packages match the colleague's
library exactly, with no hand-pinning: Seurat 5.5.0, SeuratObject 5.4.0,
Matrix 1.7-5, harmony 2.0.3.

`Matrix` agreeing matters more than it looks — it is a recommended package
bundled with R, and a Matrix/Seurat ABI mismatch is a classic source of
obscure breakage.

### Everything compiles from source — there are no binaries

Posit Package Manager has effectively no binaries for Ubuntu focal:

```
R 4.3.2 -> 15 binaries    of ~22,000 packages
R 4.6.1 -> 15 binaries
```

Identical across R versions, so this is **not** an R-version problem — focal
binaries simply no longer exist for anyone. Focal reached
end-of-standard-support in May 2025.

So the snapshot buys **reproducibility, not speed**. Full stack install was
~2-3 hours of wall clock across several passes.

**Do not use the `jammy` endpoint**, which does have binaries. They link
against glibc 2.35; this box has 2.31. They would install and then fail at load
time with confusing symbol errors.

### renv per project, with a shared cache

`renv` pins R *packages*; the source build pins *R*. Different layers, both
needed.

1. **It makes the freeze concept real.** `git tag freeze01-analysis` records
   the code; `renv.lock` records the packages. Without the lockfile the tag is
   half a promise.
2. **It insulates against the shared library changing.** Something modified
   `/usr/local/lib/R/site-library` on 2026-09-14 with no notice.
3. **It is the reproducibility artifact** reviewers and journals ask for.
4. **It makes cross-checking with collaborators tractable** — exchange
   lockfiles instead of doing filesystem archaeology.

Shared cache at `~/.cache/R/renv`: packages are hard-linked into each project's
library, so a second project costs almost no extra disk and no recompilation.
Given that everything compiles from source here, this is worth a lot.

**`snapshot.type` must be `"all"`, not `"explicit"`.** See workaround 6 below —
this was a real bug in the first attempt.

### Conservative thread caps

`~/.Renviron` sets BLAS and OpenMP threads to **4**. Shared 32-core box, no
scheduler, sustained load 12-18, ~7 users. Override per job:

```bash
OMP_NUM_THREADS=16 jobs/run.sh <script>
```

### Bootstrap library kept minimal

`~/R/library` holds only `renv`, `BiocManager`, `remotes`. Every project's real
packages live in its own `renv/library/`. This is what lets Project 1 and a
future project use different Seurat versions without fighting.

---

## Installation workarounds — READ THIS FIRST when something fails

Six problems hit during installation, in order. Each is the kind of thing that
costs an hour to rediscover. **When a package fails to build on this box, check
this list before anything else.**

### 1. Conda leakage broke the R build (ICU version clash)

**Symptom.** `make` fails at the final link with:

```
/usr/bin/ld: ../../lib/libR.so: undefined reference to `ucol_open_73'
/usr/bin/ld: ../../lib/libR.so: undefined reference to `u_getVersion_73'
collect2: error: ld returned 1 exit status
```

**Cause.** Conda's `base` env was active, so `configure` picked up
`-I/opt/miniconda3/include` and compiled against **conda's ICU 73** headers
(symbols suffixed `_73`). At link time `-licuuc` resolved to the **system ICU
66**. Headers and library disagreed.

```
/opt/miniconda3/include/unicode/uvernum.h   U_ICU_VERSION_MAJOR_NUM 73
ldconfig -p | grep libicuuc                 libicuuc.so.66
```

**Fix.** Build in a conda-free environment:

```bash
conda deactivate; conda deactivate
unset CONDA_PREFIX CONDA_DEFAULT_ENV
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
export CPPFLAGS=""
export PKG_CONFIG_PATH="/usr/lib/x86_64-linux-gnu/pkgconfig:/usr/share/pkgconfig"
```

**Verify before spending 30 minutes on `make`:**

```bash
grep -c miniconda ~/src/r_configure.log              # MUST be 0
grep -c miniconda ~/R/R-4.6.1/lib/R/etc/Makeconf     # MUST be 0
```

**Prevention, now in place.** `conda config --set auto_activate_base false`.
Conda still works via `conda activate <env>`; it just no longer shadows system
paths in every shell.

**Launch jobs with `bash -c`, NOT `bash -lc`.** The `-l` sources `.bashrc`,
which can re-activate conda inside a detached screen session. This is why
`jobs/run.sh` and the bootstrap launch commands use `bash -c`.

Related trap: `curl-config --version` reports **8.4.0** (conda's) while
`pkg-config --modversion libcurl` reports **7.68.0** (the system's). Trust
pkg-config.

### 2. `--vanilla` hides the user library

**Symptom.** Packages install successfully, then a verification line fails:

```
* DONE (BiocManager)
Error in loadNamespace(x) : there is no package called 'BiocManager'
```

**Cause.** `Rscript --vanilla` ignores `~/.Renviron`, so `R_LIBS_USER` is never
set and `~/R/library` is absent from `.libPaths()`. `install.packages(lib=...)`
worked because it names the library explicitly; `BiocManager::version()`
searched `.libPaths()` and found nothing.

**Fix.** Do not use `--vanilla` when the user library matters. Fixed in
`scripts/bootstrap_server.sh` (commit `2b6602b`).

**General lesson.** When a script reports failure, check whether the *work*
failed or only a *verification line* failed. This happened twice.

### 3. `fs` needs libuv, which is not installed

**Symptom.**

```
Configuration failed because libuv was not found.
<stdin>:1:10: fatal error: uv.h: No such file or directory
ERROR: configuration failed for package 'fs'
```

**Cause.** `libuv1-dev` is absent and needs root.

**Fix.** The package bundles a static libuv. Now permanent in `~/.Renviron`:

```
USE_BUNDLED_LIBUV=1
```

**General lesson.** Most packages handle a missing system library by bundling a
fallback behind an environment variable. **Read the `[CONFIGURE]` block** — it
usually names the escape hatch — before concluding you are blocked.

### 4. `curl` needs `/sbin` on PATH

**Symptom.**

```
Local libcurl is too old. Downloading a static libcurl for legacy Linux...
./configure: ./get-curl-linux.sh: ldconfig: not found
No suitable OpenSSL or NSS found
handle.c:489:51: error: 'CURLINFO_EFFECTIVE_METHOD' undeclared
ERROR: compilation failed for package 'curl'
```

**Cause.** Two-stage. System libcurl is **7.68**; `CURLINFO_EFFECTIVE_METHOD`
arrived in **7.72**. The package detected this and tried its fallback —
downloading a static modern libcurl — but that script needs `ldconfig`, which
lives in `/sbin`, and the sanitized PATH from workaround 1 excluded it.

**Fix.** Include `/sbin` and `/usr/sbin`:

```bash
export PATH="$HOME/R/R-4.6.1/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"
```

`curl` 7.1.0 then built in 7 seconds. **This is why `curl` matters:**
`curl` → `httr` → `plotly` → **Seurat**. Not skippable.

### 5. gcc 9.4 → C++17 ceiling (blocks `scDblFinder`)

**Symptom.**

```
run_pca.cpp:97:48: error: missing template arguments before 'opt'
   scran_pca::SubsetPcaBlockedOptions opt;
ERROR: compilation failed for package 'scrapper'
Error: failed to install "scDblFinder"
```

**Cause.** `scrapper` uses class template argument deduction for aggregates, a
**C++20** feature. gcc 9 implements complete C++17 and no C++20. Not a missing
library — no environment variable fixes it.

**What was tried and did not work:**

- **Pinning an older `scrapper`** (`bioc::scrapper@1.4.0`) → "failed to
  download". Bioconductor 3.23 serves only the current version; archives are
  not reachable this way.
- **conda-forge gcc 16.2 via `~/.R/Makevars`** → compiled everything, then
  failed at link:
  ```
  x86_64-conda-linux-gnu-ld: cannot find -lgfortran
  ```
  conda's `ld` searches only conda's library paths, and the `gcc_linux-64`
  package does not ship `libgfortran`.

**If you retry the conda-toolchain route**, two things are needed:

```bash
conda create -p ~/envs/gcc13 -c conda-forge \
  gcc_linux-64 gxx_linux-64 gfortran_linux-64
# and in ~/.R/Makevars:
#   LDFLAGS = -L/usr/lib/gcc/x86_64-linux-gnu/9 -L$HOME/envs/gcc13/lib \
#             -Wl,-rpath,$HOME/envs/gcc13/lib
```

`~/envs/gcc13` (gcc 16.2.0, 168 MB) is **kept** as the starting point for any
future C++20 requirement.

**CRITICAL: remove `~/.R/Makevars` when finished.** Leaving it in place made
`nnls` and `SparseM` — plain Fortran packages — fail with the same
`-lgfortran` error, breaking the MuSiC dependency chain. Verify with:

```bash
R CMD config CXX    # expect: g++ -std=gnu++17
R CMD config FC     # expect: gfortran
```

**Current resolution.** `scDblFinder` is **not installed**. It is commented out
in `scripts/bootstrap_renv.R` with a pointer here. Three alternatives:
`DoubletFinder` (GitHub, Seurat-native), `scds` (Bioconductor), or inherit the
colleague's existing doublet calls — she has `scDblFinder 1.16.0` results on
this very dataset.

### 6. `renv` `snapshot.type = "explicit"` recorded only 1 package

**Symptom.** Install completes, then:

```
[ ok ] renv.lock written
packages recorded: 1
```

**Cause.** `"explicit"` mode records only packages declared as dependencies in
a `DESCRIPTION` file. This repo has no `DESCRIPTION`, so almost nothing was
recorded. `renv::restore()` elsewhere would have installed nothing — defeating
the entire purpose.

**Fix.**

```r
renv::settings$snapshot.type("all")
renv::snapshot(prompt = FALSE)
```

277 packages recorded. **Always verify the count after a snapshot:**

```r
lk <- jsonlite::fromJSON("renv.lock"); length(lk$Packages)
```

**Minor known quirk.** The lockfile records CRAN as
`https://packagemanager.posit.co/cran/2026-06-01`, dropping the
`__linux__/focal/` segment. Harmless here (focal has no binaries either way),
but on a distro *with* binaries, `renv::restore()` would fetch source instead.

---

## What is installed

```
~/R/R-4.6.1/          R, source-built                      ~306 MB
~/R/library/          bootstrap only: renv, BiocManager, remotes
~/opt/openblas/       OpenBLAS 0.3.34 (DYNAMIC_ARCH, pthread)
~/envs/gcc13/         conda-forge gcc 16.2 — C++20 escape hatch, unused
~/.cache/R/renv/      shared renv cache; all projects hard-link here
~/src/                build trees + logs (safe to delete)
~/projects/asset-skin-natural-history/    the repo, with its own renv/library
~/.conda/envs/        Python envs (scVI, TCAT, Scaden) — not yet created
```

### Configuration files

`~/.Renviron`

```
OPENBLAS_NUM_THREADS=4
OMP_NUM_THREADS=4
MKL_NUM_THREADS=4
RENV_PATHS_CACHE=/home/parvizi/.cache/R/renv
R_LIBS_USER=/home/parvizi/R/library
USE_BUNDLED_LIBUV=1          # workaround 3
```

`~/.bashrc` — two managed blocks: the R PATH block (prepended), and an
`ssh-agent` block guarded by `[[ $- == *i* ]]` so it never blocks a detached
job on a passphrase prompt.

`~/.R/Makevars` — **must not exist** in normal operation. See workaround 5.

### Verified functional

| Check | Result |
|---|---|
| R version | 4.6.1 (2026-06-24) |
| Bioconductor | 3.23 |
| BLAS | `libopenblasp-r0.3.34.so` |
| `stringi` ICU | 66.1 (system, not conda's 73) |
| Seurat 5 + BPCells on-disk | object created, PCA ran |
| `dream()` mixed model | `~ timepoint + (1|subject)` fitted, 12 subjects × 2 timepoints |
| Packages installed | 292 |
| Packages in `renv.lock` | 277 |

### Pinned to match the handoff objects

Seurat 5.5.0 · SeuratObject 5.4.0 · Matrix 1.7-5 · harmony 2.0.3 ·
BPCells 0.3.1 (r-universe, SHA `adc4a3c30f60a03522f58947d733d7d77a6eb2cf`)

### Core analysis set

**Bioconductor 3.23:** DESeq2 1.52.0 · limma 3.68.5 · edgeR 4.10.5 ·
variancePartition 1.42.0 · speckle 1.12.0 · SingleR 2.14.2 · celldex 1.22.0 ·
SingleCellExperiment 1.34.0 · scater 1.40.2 · scran 1.40.0 ·
glmGamPoi 1.24.0 · BiocParallel 1.46.0 · HDF5Array · DelayedArray

**CRAN:** presto 1.1.0 · lme4 2.0-1 · lmerTest 3.2-1 · pbkrtest ·
data.table 1.18.4 · dplyr · tidyr · readr · ggplot2 4.0.3 · patchwork ·
ggrepel · viridis · pheatmap · cowplot · yaml · jsonlite · future ·
future.apply · furrr · reticulate

### Deliberately NOT installed

Each is heavy, none is needed for Project 1, and all can live in their own
renv project so their dependency trees cannot perturb this library.

| Package | Install with | Note |
|---|---|---|
| CellChat | `renv::install("Jin-s-Lab/CellChat")` | R, not Python. Sprawling deps. |
| monocle3 | `renv::install("cole-trapnell-lab/monocle3")` | pseudotime; heavy |
| SeuratWrappers | `renv::install("satijalab/seurat-wrappers")` | bridge to external tools; Harmony has native Seurat support |
| SeuratDisk | `renv::install("mojaveazure/seurat-disk")` | h5ad; thinly maintained, imperfect Seurat 5 support |
| zellkonverter | `renv::install("bioc::zellkonverter")` | **better h5ad route** — via SingleCellExperiment, actively maintained |

When the R ↔ Python round-trip is needed (scVI, TCAT), install SeuratDisk
**and** zellkonverter, test both against the real objects, and record which one
works in `docs/decisions.md`. Conversion is where these pipelines break.

---

## Known gaps

### `scDblFinder` — blocked by the gcc 9 / C++20 ceiling

See workaround 5. Alternatives: `DoubletFinder`, `scds`, or the colleague's
existing calls. Affects `analysis/04_sc_ambient/`.

### `MuSiC` — blocked by `TOAST` disappearing from Bioconductor 3.23

`MuSiC` downloaded fine; its dependency `TOAST` returns "failed to download",
presumably deprecated in 3.23. Affects `analysis/15_deconvolution/`.

Options when needed: install `TOAST` from the Bioconductor archive, vendor it
from GitHub, use `Scaden` (Python, already in the plan), or `CIBERSORTx`. The
colleague has MuSiC 1.0.0 working as a reference point.

Not urgent — deconvolution is well downstream.

### Toolchain ceiling will bite again

gcc 9 gives complete C++17, no C++20. Packages published from 2026 onward
increasingly want C++20. When one is load-bearing, the route is workaround 5
*with* `gfortran_linux-64`, applied only for that package, then Makevars
removed.

### No containers, no optimized system BLAS, focal is EOL

All need root. Worth requesting Apptainer.

### Missing system libraries

`libglpk`, `libgit2`, `libmagick` absent. None blocks the current stack —
igraph bundles its own GLPK, git2r only affects devtools plumbing, magick is
optional plotting.

---

## Reproducing this environment

On a new server:

```bash
git clone git@github.com:rezvanparvizi/asset-skin-natural-history.git
cd asset-skin-natural-history
chmod +x scripts/*.sh

# conda must NOT be active — see workaround 1
conda config --set auto_activate_base false
exit   # then log in again

# 1. R + OpenBLAS + bootstrap packages + smoke test (~1 hour)
screen -dmS bootstrap bash -c 'scripts/bootstrap_server.sh 2>&1 | tee ~/bootstrap.log'
tail -f ~/bootstrap.log

# 2. new shell so PATH takes effect
exec bash -l
which R && R --version | head -1
Rscript -e 'cat(sessionInfo()$BLAS, "\n")'

# 3. project library from the lockfile (~2-3 hours; source compilation)
screen -dmS renv bash -c 'cd ~/projects/asset-skin-natural-history; \
  unset CONDA_PREFIX CONDA_DEFAULT_ENV; \
  export PATH="$HOME/R/R-4.6.1/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"; \
  export ASSET_NCPUS=8 USE_BUNDLED_LIBUV=1; \
  Rscript -e "renv::restore()" 2>&1 | tee ~/renv_restore.log'
```

Note the launch pattern: `bash -c` not `bash -lc`, conda variables unset,
`/sbin` on PATH, `USE_BUNDLED_LIBUV=1` set. Those are workarounds 1, 3 and 4.

On the laptop:

```bash
./scripts/bootstrap_laptop.sh            # level (b): dev + single-cell
./scripts/bootstrap_laptop.sh --minimal  # level (a): editing tools only
```

---

## Routine operations

**Add a package**

```r
renv::install("somePackage")        # or "bioc::x", or "user/repo"
renv::snapshot()
```
```bash
# verify the count did not collapse — see workaround 6
Rscript -e 'cat(length(jsonlite::fromJSON("renv.lock")$Packages), "\n")'
git add renv.lock && git commit -m "renv: add somePackage" && git push
```

**Snapshot at a data freeze**

```bash
Rscript -e 'renv::snapshot(prompt = FALSE)'
git add renv.lock && git commit -m "renv: snapshot for freeze01"
git tag -a freeze01-env -m "Package state used for freeze01 results"
git push --tags
```

**Check which R is active**

```bash
which R                                       # expect ~/R/R-4.6.1/bin/R
Rscript -e 'cat(sessionInfo()$BLAS, "\n")'    # expect .../openblas
R CMD config CXX                              # expect g++ -std=gnu++17
```

**Raise threads for one heavy job**

```bash
OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 jobs/run.sh <script>
```

**Recover a broken library**

```r
renv::restore()
```

**Is it running or did it crash?**

```bash
screen -ls                     # a session listed = running
tail -3 ~/<log>                # ends in FAIL / Done = over
ps -u $USER | grep -E "R$|Rscript|make"
```

---

## What was surveyed before deciding

Recorded so nobody repeats the investigation.

| Question | Finding |
|---|---|
| System R | 4.3.2 → Bioc 3.18, ~3 years stale |
| Shared site-library | 607 packages, Seurat **4.3.0**, no BPCells, modified 2026-09-14 |
| Writable? | `/usr/lib/R/site-library` not writable by us |
| Colleague's setup | own R 4.6.0 in `~/R`, 636 pkgs / 5.9 GB; per-tool conda envs in `~/.conda/envs`; a `conda_create_r()` wrapper writing `R_LIBS_USER` into each env's `Renviron.site`; **no renv** |
| Colleague's PATH | appends her R → a bare `R` runs system 4.3.2 (we prepend) |
| BLAS available | reference only; no OpenBLAS/ATLAS/MKL |
| Compilers | gcc/g++/gfortran 9.4.0, cmake 3.31.8 → C++17 ceiling |
| R build headers | zlib, bzlib, lzma, pcre2, curl, readline — all present |
| HDF5 | 1.10.4 serial; working BPCells links `libhdf5_serial.so.103` |
| GSL | `libgsl.so.23` + headers present |
| libcurl | system 7.68 (pkg-config); conda's `curl-config` misreports 8.4.0 |
| ICU | system 66; conda 73 — cause of workaround 1 |
| Missing libs | libglpk, libgit2, libmagick |
| P3M focal binaries | 15 of ~22,000 — effectively none, for any R version |
| Snapshot matching handoff | `focal/2026-06-01` → Seurat 5.5.0, SeuratObject 5.4.0, Matrix 1.7-5, harmony 2.0.3 |
| Containers | none, root needed |
| Root access | none; admin upgrade not coming soon |
| Disk | 22 TB free on `/home`, no enforced quota; `/hits` separate NFS at 99% |
| Hardware | 32 cores, 503 GB RAM, load 12-18, ~7 users |
| Pre-existing R config | none — clean slate before this setup |
| git on server | 2.25 — `git init -b` unsupported; use `git init` then `git checkout -b main` |
