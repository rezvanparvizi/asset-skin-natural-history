# Computing environment — what, where, and why

Written so that someone arriving cold (a collaborator, a future you, an AI
agent) understands this setup without reverse-engineering it from the
filesystem. Every non-obvious choice below has a reason recorded next to it.

Last reviewed: 2026-09-24 · Servers: `mininubio.ddns.med.umich.edu`

---

## Contents

- [Summary](#summary)
- [The machines](#the-machines)
- [Decisions and reasoning](#decisions-and-reasoning)
- [What is installed where](#what-is-installed-where)
- [Known limitations](#known-limitations)
- [How to reproduce this environment](#how-to-reproduce-this-environment)
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

The system R (4.3.2) and the shared site-library are **not used**. Reasons
below.

---

## The machines

### Server: mininubio

```
Ubuntu 20.04.6 LTS   glibc 2.31   x86_64
32 cores, 503 GB RAM
/home on /dev/md0, 40 TB, no enforced quota
gcc / g++ / gfortran 9.4.0     cmake 3.31.8
No root access. No job scheduler — jobs run under screen.
Typical load: 16-18 sustained, ~7 concurrent users.
```

Pre-existing R installations, none of which we use:

- `/usr/bin/R` — 4.3.2 (Oct 2023), Bioconductor **3.18**
- `/usr/local/lib/R/site-library` — 607 packages, admin-maintained, shared by
  all users, **modified 2026-09-14 without notice**
- `/opt/miniconda3/envs/scRNAseqCelltypeAnnotation` — root-owned, contains R
- A colleague's `~/R/R-4.6.0` + `~/R/library` (636 packages) — the precedent
  this setup follows

### Laptop: U-M MacBook Pro (Apple Silicon)

```
macOS, arm64
R 4.6.1 from the CRAN arm64 pkg
gfortran 14.2 from mac.r-project.org
Xcode CLT installed, license accepted
Homebrew for HDF5 and friends
```

Used for editing, git, and developing functions against small test objects.
Never for real analysis — the data is on the server.

---

## Decisions and reasoning

### R built from source in `$HOME`, not system R, not conda

System R is 4.3.2, which pins Bioconductor to **3.18 (October 2023)**. Nearly
three years stale in a field that moves fast, and the ceiling is hard: R
version determines Bioconductor version with no override.

Rejected alternatives:

- **System R + a user library.** Fastest, but (a) Bioc stays 3.18, and (b) the
  shared site-library holds Seurat 4.3.0 / SeuratObject 4.1.3 while we need
  Seurat 5, so the search path would carry two incompatible major versions and
  607 packages compiled against the old one. Mixed-ABI stacks fail in ways that
  are miserable to diagnose.
- **Asking an admin to upgrade.** Asked; not happening on a useful timescale.
- **Conda R.** Works and is proven on this box, but always compiles from
  source, ships its own duplicate system libraries, and buys nothing over a
  source build here — because everything compiles from source anyway (see
  below). Kept for Python tools, where it is the right instrument.
- **Containers.** Technically the best answer for reproducibility across the
  group's 6-7 servers. Rejected because `apptainer`, `singularity`, `docker`
  and `podman` are all absent and installing them needs root. **Worth
  requesting** — Apptainer is the standard academic-HPC ask and does not
  require users to have root.

Source build gives current R (4.6.1 → Bioc 3.23), needs no permission, links
against the system's existing libraries, and is already proven on this exact
machine by a colleague's successful 4.6.0 build.

### R 4.6.1 specifically

Matches the laptop exactly. R keys its library directory to `major.minor`, so
4.6.0 and 4.6.1 share `library/4.6` and are package-compatible — the
colleague's 4.6.0-built objects load without complaint.

**Policy: keep both machines on 4.6.x and upgrade them together, deliberately.**
If macOS R auto-updates to 4.7 while the server stays on 4.6, `renv.lock`
stops restoring and laptop/server parity is silently gone.

### PATH is *prepended*, unlike the local precedent

The colleague's `.bashrc` appends her R to `PATH`:

```bash
path_append "$HOME/R/R-4.6.0/bin"   # comment: "NOT prepended — conda env R must win"
```

Appending leaves `/usr/bin` earlier in the path, so a bare `R` runs the system
**4.3.2**; her 4.6.0 only runs when invoked by full path. Her stated goal —
letting an active conda env take precedence — is achieved anyway by
prepending, because `conda activate` prepends the env's `bin` *later* in the
shell's lifetime and therefore still wins.

So we prepend. A bare `R` is 4.6.1; conda envs override when active. Verify
with `which R`.

### OpenBLAS built in `$HOME`

The system has only the **reference** BLAS:

```
libblas.so.3 -> /usr/lib/x86_64-linux-gnu/blas/libblas.so.3   (priority 10, sole option)
```

No OpenBLAS, ATLAS or MKL anywhere, which is why the colleague's R shows
`BLAS_LIBS = -lblas`.

Honest scope of the benefit: most of the sparse single-cell path (irlba PCA on
sparse matrices, neighbour graphs, UMAP, marker tests) is **not**
dense-BLAS-bound, so this is not a blanket speedup. It matters concretely for
Harmony (dense operations on the embedding, repeated k-means), `ScaleData`,
dense correlation matrices, and any dense linear algebra — steps run many
times over the project.

Build flags and why:

- `DYNAMIC_ARCH=1` — one binary that selects the right kernel at runtime, so
  this install is portable to the group's other servers with different CPUs.
- `USE_OPENMP=0` — pthread threading. An OpenMP-threaded BLAS fights with R
  packages that use OpenMP themselves (`data.table`, `glmGamPoi`), causing
  thread oversubscription. Threads are capped at runtime instead.
- `NUM_THREADS=64` — compile-time ceiling, not the runtime default.

R is linked with `-Wl,-rpath,$HOME/opt/openblas/lib`, so it finds the library
without needing `LD_LIBRARY_PATH` set.

### CRAN pinned to the P3M snapshot `focal/2026-06-01`

```
https://packagemanager.posit.co/cran/__linux__/focal/2026-06-01
```

Chosen because it is the **latest date at which Seurat is still 5.5.0** — the
version the handoff objects were written with. Verified:

```
2026-05-01  Seurat 5.5.0
2026-06-01  Seurat 5.5.0   <- chosen
2026-07-01  Seurat 5.5.1
```

At that snapshot, all four object-critical packages match the colleague's
library exactly, with no hand-pinning needed:

| Package | Snapshot | Her library |
|---|---|---|
| Seurat | 5.5.0 | 5.5.0 |
| SeuratObject | 5.4.0 | 5.4.0 |
| Matrix | 1.7-5 | 1.7-5 |
| harmony | 2.0.3 | 2.0.3 |

`Matrix` agreeing matters more than it looks — it is a recommended package
bundled with R, and a Matrix/Seurat ABI mismatch is a classic source of
obscure breakage.

Pinning a dated snapshot pins the **entire dependency graph** reproducibly,
rather than four packages by hand with everything else floating.

### Everything compiles from source — there are no binaries

Posit Package Manager has effectively no binaries for Ubuntu focal:

```
R 4.3.2 -> 15 binaries    of ~22,000 packages
R 4.6.1 -> 15 binaries
```

Identical across R versions, so this is not an R-version problem: focal
binaries simply no longer exist for anyone. Focal reached end-of-standard-
support in May 2025.

Consequences, all accounted for:

- Full stack install is ~2-3 hours, once, on 32 cores.
- The smoke-test-first ordering in the bootstrap scripts is **load-bearing**,
  not decorative.
- The shared renv cache matters: compile once, hard-link into every future
  project.

**Do not be tempted by the `jammy` endpoint**, which does have binaries. Jammy
binaries link against glibc 2.35; this box has 2.31. They would install and
then fail at load time with confusing symbol errors.

### Bioconductor 3.23

Follows automatically from R 4.6.x. No action needed, and it matches the
colleague's environment.

### BPCells from r-universe

Not on CRAN. Installed from `https://bnprks.r-universe.dev`, the same source
as the 0.3.1 build already working on this machine:

```
Version:      0.3.1
Repository:   https://bnprks.r-universe.dev
RemoteUrl:    https://github.com/bnprks/BPCells
RemoteSha:    adc4a3c30f60a03522f58947d733d7d77a6eb2cf
RemoteSubdir: r
Built:        R 4.6.0; x86_64-pc-linux-gnu
```

The SHA is recorded here and in `renv.lock` as a fallback, because r-universe
builds can rotate. BPCells links against the system HDF5 **1.10.4**
(`libhdf5_serial.so.103`) — confirmed by `ldd` on the working build.

### renv per project, with a shared cache

`renv` pins R *packages*; conda or the source build pins *R*. Different layers,
both needed.

What renv delivers here:

1. **It makes the freeze concept real.** `git tag freeze01-analysis` records
   the code; `renv.lock` records the packages. Together, freeze01 results stay
   reproducible. Without the lockfile the tag is half a promise.
2. **It insulates against the shared library changing.** Something modified
   `/usr/local/lib/R/site-library` on 2026-09-14 with no notice. A
   project-local library cannot be reached by that.
3. **It is the reproducibility artifact** reviewers and journals ask for — a
   small text file in the repo that fully answers "what versions did you use".
4. **It makes cross-checking with collaborators tractable** — exchange
   lockfiles instead of doing filesystem archaeology.

Shared cache at `~/.cache/R/renv`: packages are hard-linked into each
project's library, so a second project costs almost no extra disk and no
recompilation. Given that everything compiles from source here, this is worth
a lot.

`snapshot.type = "explicit"` — `renv.lock` records only the packages we
declare, not everything the code happens to load. Keeps the library a
deliberate decision.

### Conservative thread caps

`~/.Renviron` sets BLAS and OpenMP threads to **4**. This is a shared 32-core
box with no scheduler, sustained load 16-18, and ~7 concurrent users; one R
process taking all 32 cores harms everyone. Override deliberately per job:

```bash
OMP_NUM_THREADS=16 jobs/run.sh <script>
```

No lab-wide norm exists on this. If one is agreed later, update here.

### Bootstrap library kept minimal

`~/R/library` holds only `renv`, `BiocManager`, `remotes` — the tooling needed
to drive renv. Every project's real packages live in its own
`renv/library/`. This is what lets Project 1 and a future project use
different Seurat versions without fighting.

---

## What is installed where

```
~/R/R-4.6.1/                R installation (source-built)   ~306 MB
~/R/library/                bootstrap only: renv, BiocManager, remotes
~/opt/openblas/             OpenBLAS (DYNAMIC_ARCH, pthread)
~/.cache/R/renv/            shared renv cache — all projects hard-link here
~/src/                      build trees and logs (safe to delete after install)
~/projects/<repo>/          each repo, each with its own renv/library
~/asset-data/               data; the repo's data/ symlink points here
~/.conda/envs/              Python envs (scVI, TCAT, Scaden) — set up later
```

### Package set

**Pinned to match the handoff objects:** Seurat 5.5.0 · SeuratObject 5.4.0 ·
Matrix 1.7-5 · harmony 2.0.3 · BPCells 0.3.1

**Bioconductor 3.23:** DESeq2 · limma · edgeR · variancePartition · speckle ·
scDblFinder · SingleR · celldex · SingleCellExperiment · SummarizedExperiment ·
scater · scran · glmGamPoi · BiocParallel

**CRAN:** presto · lme4 · data.table · dplyr · tidyr · tibble · stringr ·
purrr · readr · ggplot2 · patchwork · ggrepel · RColorBrewer · viridis ·
scales · pheatmap · cowplot · yaml · jsonlite · R.utils · future ·
future.apply · furrr

**GitHub:** MuSiC (`xuranw/MuSiC`) — bulk deconvolution for
`analysis/13_bulk_skin_deconvolution/`

### Deliberately NOT installed

Each is heavy, none is needed for Project 1, and all can live in their own
renv project so their dependency trees cannot perturb this library.

| Package | Install with | Note |
|---|---|---|
| CellChat | `renv::install("Jin-s-Lab/CellChat")` | R package, not Python. Sprawling deps, historically awkward. |
| monocle3 | `renv::install("cole-trapnell-lab/monocle3")` | pseudotime; heavy |
| SeuratWrappers | `renv::install("satijalab/seurat-wrappers")` | bridge to external tools; Harmony has native Seurat support, so mainly useful for monocle3 conversion |
| SeuratDisk | `renv::install("mojaveazure/seurat-disk")` | h5ad conversion; thinly maintained, imperfect Seurat 5 support |
| zellkonverter | `renv::install("bioc::zellkonverter")` | **better h5ad route than SeuratDisk** — via SingleCellExperiment, actively maintained |

When the R ↔ Python round-trip is needed (scVI, TCAT), install SeuratDisk
**and** zellkonverter, test both against the real objects, and record which one
handles them in `docs/decisions.md`. Conversion is where these pipelines break.

---

## Known limitations

### Toolchain ceiling: gcc 9.4 → C++17

gcc 9 provides complete C++17, partial C++20, no C++23. The colleague's R
`Makeconf` shows `CC23 =` empty, confirming no C23 compiler was detected.

Current stack is unaffected — BPCells needs C++17 and works. But packages
published from 2026 onward increasingly want C++20, and this is the most
likely way the environment will eventually bite.

**When a package demands C++20:**

1. Check whether a binary exists elsewhere first (a newer distro codename, or
   r-universe), since those are built with modern toolchains.
2. Otherwise install a newer compiler into `$HOME` via conda-forge:
   ```bash
   conda create -p ~/envs/gcc13 -c conda-forge gcc_linux-64 gxx_linux-64
   ```
   Then point R at it for that package only, via `~/.R/Makevars`. Mixing
   toolchains has ABI caveats, so treat it as a targeted fix, never a default.
3. Or request that an admin upgrade the OS. Ubuntu 20.04 is past
   end-of-standard-support, so this may already be planned.

### No containers

`apptainer`, `singularity`, `docker`, `podman` all absent; installing them
needs root. Worth requesting — it is the cleanest path to identical
environments across the group's servers.

### No optimized BLAS system-wide

Addressed by the `$HOME` OpenBLAS build, but note it is per-user: anyone else
on this box is still on reference BLAS.

### Missing system libraries

`libglpk`, `libgit2`, `libmagick` are absent. None blocks the current stack —
igraph bundles its own GLPK, git2r only affects devtools plumbing, magick is
optional plotting. Installing them needs root.

### Focal is EOL

No security updates for the base OS, and no package-manager binaries. Both
are consequences of the same fact. Nothing to do about it without admin.

---

## How to reproduce this environment

On a new server (any of the group's boxes):

```bash
git clone git@github.com:<username>/asset-skin-natural-history.git
cd asset-skin-natural-history
chmod +x scripts/*.sh

# 1. R + OpenBLAS + bootstrap packages + smoke test  (~1 hour)
screen -dmS bootstrap bash -lc 'scripts/bootstrap_server.sh 2>&1 | tee ~/bootstrap.log'
tail -f ~/bootstrap.log

# 2. new shell so PATH takes effect
exec bash -l
which R && R --version | head -1

# 3. project library from the lockfile  (~2-3 hours; source compilation)
screen -dmS renv bash -lc 'Rscript scripts/bootstrap_renv.R 2>&1 | tee ~/renv_install.log'
tail -f ~/renv_install.log
```

Once `renv.lock` exists and is committed, step 3 on any new machine becomes:

```r
renv::restore()
```

which reproduces the exact package set.

On the laptop:

```bash
./scripts/bootstrap_laptop.sh          # level (b): dev + single-cell
./scripts/bootstrap_laptop.sh --minimal  # level (a): editing tools only
```

---

## Routine operations

**Add a package to the project**

```r
renv::install("somePackage")      # or "bioc::x", or "user/repo"
renv::snapshot()
# then: git add renv.lock && git commit -m "renv: add somePackage"
```

**Take a package snapshot at a data freeze**

```r
renv::snapshot()
```
```bash
git add renv.lock
git commit -m "renv: snapshot for freeze01"
git tag -a freeze01-env -m "Package state used for freeze01 results"
```

**Check which R is active**

```bash
which R          # expect ~/R/R-4.6.1/bin/R
R --version | head -1
Rscript -e 'cat(sessionInfo()$BLAS, "\n")'   # expect .../openblas
```

**Raise threads for one heavy job**

```bash
OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 jobs/run.sh <script>
```

**Recover from a broken project library**

```r
renv::restore()   # rebuild from renv.lock
```

**Check machine load before a big run**

```bash
jobs/status.sh
```

---

## What was surveyed before deciding

Recorded so nobody repeats the investigation.

| Question | Finding |
|---|---|
| System R version | 4.3.2 → Bioc 3.18, ~3 years stale |
| Shared site-library | 607 packages, Seurat **4.3.0**, no BPCells, modified 2026-09-14 |
| Writable? | `/usr/lib/R/site-library` not writable by us |
| Colleague's setup | own R 4.6.0 in `~/R`, 636 pkgs, 5.9 GB; per-tool conda envs in `~/.conda/envs`; a `conda_create_r()` wrapper writing `R_LIBS_USER` into each env's `Renviron.site`; **no renv** |
| Colleague's PATH | appends her R → a bare `R` runs system 4.3.2 (we prepend instead) |
| BLAS available | reference only; no OpenBLAS/ATLAS/MKL |
| Compilers | gcc/g++/gfortran 9.4.0, cmake 3.31.8 → C++17 ceiling |
| R build headers | zlib, bzlib, lzma, pcre2, curl, readline — all present |
| HDF5 | 1.10.4 serial; working BPCells links `libhdf5_serial.so.103` |
| GSL | `libgsl.so.23` + headers present |
| Missing libs | libglpk, libgit2, libmagick |
| P3M focal binaries | 15 of ~22,000 — effectively none, for any R version |
| Snapshot matching handoff | `focal/2026-06-01` → Seurat 5.5.0, SeuratObject 5.4.0, Matrix 1.7-5, harmony 2.0.3 |
| Containers | none, and root needed to install |
| Root access | none on the server; admin upgrade not coming soon |
| Disk | 22 TB free on `/home`, no enforced quota |
| Hardware | 32 cores, 503 GB RAM, load 16-18, ~7 users |
| Pre-existing R config | no `~/.Renviron`, `~/.Rprofile` or `~/.R/Makevars` — clean slate |
