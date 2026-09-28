#!/usr/bin/env Rscript
# =============================================================================
# bootstrap_renv.R
#
# Builds this project's renv library and writes renv.lock.
#
# Run AFTER scripts/bootstrap_server.sh, from the repository root:
#   cd ~/projects/asset-skin-natural-history
#   screen -dmS renv bash -lc 'Rscript scripts/bootstrap_renv.R 2>&1 | tee ~/renv_install.log'
#   tail -f ~/renv_install.log
#
# Expect 2-3 hours. Posit Package Manager has no binaries for Ubuntu focal
# (verified: 15 of ~22,000 packages), so every package compiles from source.
# The shared renv cache means you pay this once, not per project.
#
# Version strategy
#   CRAN        pinned to the P3M snapshot 2026-06-01, chosen because it is the
#               latest date where Seurat is 5.5.0 — matching the objects being
#               handed over. Pinning the snapshot pins the WHOLE dependency
#               graph, not just a few packages by hand.
#   Bioconductor 3.23, determined automatically by R 4.6.x.
#   BPCells     r-universe (bnprks.r-universe.dev), the same source used for
#               the 0.3.1 build already working on this machine.
#
# Rationale in docs/ENVIRONMENT.md.
# =============================================================================

options(warn = 1)

CRAN_SNAPSHOT   <- "https://packagemanager.posit.co/cran/__linux__/focal/2026-06-01"
BPCELLS_UNIVERSE <- "https://bnprks.r-universe.dev"
BPCELLS_SHA     <- "adc4a3c30f60a03522f58947d733d7d77a6eb2cf"  # provenance only
NCPUS           <- as.integer(Sys.getenv("ASSET_NCPUS", "8"))

hr   <- function() cat(strrep("-", 72), "\n")
step <- function(...) { cat("\n"); hr(); cat("==>", ..., "\n"); hr() }
ok   <- function(...) cat("  [ ok ]", ..., "\n")
warn <- function(...) cat("  [warn]", ..., "\n")

# ----- sanity ----------------------------------------------------------------
if (!file.exists("config/paths.R")) {
  stop("Run this from the repository root (config/paths.R not found).",
       call. = FALSE)
}
if (getRversion() < "4.6.0") {
  stop("Expected R >= 4.6.0, found ", getRversion(),
       ". Did you open a new shell after bootstrap_server.sh? ",
       "Check `which R`.", call. = FALSE)
}

step("Environment")
cat("  R          ", R.version.string, "\n")
cat("  platform   ", R.version$platform, "\n")
cat("  BLAS       ", sessionInfo()$BLAS %||% "unknown", "\n")
cat("  cpus       ", NCPUS, "\n")
cat("  snapshot   ", CRAN_SNAPSHOT, "\n")

`%||%` <- function(a, b) if (is.null(a)) b else a

options(repos = c(CRAN = CRAN_SNAPSHOT), Ncpus = NCPUS)

if (!requireNamespace("renv", quietly = TRUE))
  stop("renv not found. Run scripts/bootstrap_server.sh first.", call. = FALSE)

# ----- initialise renv -------------------------------------------------------
step("renv::init")

# bare = TRUE: create the project library but do not scan the code and install
# whatever it finds. We want an explicit, curated package list.
if (!file.exists("renv.lock") && !dir.exists("renv")) {
  renv::init(bare = TRUE, restart = FALSE, settings = list(
    snapshot.type = "explicit"   # renv.lock records only what we declare
  ))
  ok("renv initialised")
} else {
  ok("renv already present")
}

# Make the snapshot repo sticky for this project, so every future R session in
# it resolves to the same package versions.
prof <- ".Rprofile"
marker <- "# ASSET: pinned CRAN snapshot"
if (!any(grepl(marker, readLines(prof, warn = FALSE), fixed = TRUE))) {
  cat(sprintf('\n%s\noptions(repos = c(CRAN = "%s",\n                  BPCELLS = "%s"))\n',
              marker, CRAN_SNAPSHOT, BPCELLS_UNIVERSE),
      file = prof, append = TRUE)
  ok("pinned repos written to .Rprofile")
}
options(repos = c(CRAN = CRAN_SNAPSHOT, BPCELLS = BPCELLS_UNIVERSE))

# ----- package sets ----------------------------------------------------------
# Kept explicit rather than inferred from the code, so the library is a
# deliberate decision that can be reviewed.

cran_core <- c(
  # single cell — versions come from the snapshot and match the handoff objects
  "Seurat",           # 5.5.0
  "SeuratObject",     # 5.4.0
  "Matrix",           # 1.7-5  (bundled/recommended; snapshot agrees)
  "harmony",          # 2.0.3
  # modelling
  "lme4",             # variancePartition dependency, also used directly
  # data wrangling
  "data.table", "dplyr", "tidyr", "tibble", "stringr", "purrr", "readr",
  # plotting
  "ggplot2", "patchwork", "ggrepel", "RColorBrewer", "viridis", "scales",
  "pheatmap", "cowplot",
  # infrastructure
  "yaml", "jsonlite", "R.utils", "future", "future.apply", "furrr"
)

bioc_core <- c(
  "SingleCellExperiment", "SummarizedExperiment",
  "DESeq2", "limma", "edgeR",
  "variancePartition",     # dream() — longitudinal pseudobulk
  "speckle",               # propeller — composition testing
  # "scDblFinder",   # BLOCKED: dep scrapper needs C++20 (gcc 9 ceiling).
                     # See docs/ENVIRONMENT.md. Alternatives: DoubletFinder,
                     # scds; or inherit QC calls from the colleague.
  "SingleR", "celldex",
  "scater", "scran",
  "glmGamPoi",
  "BiocParallel"
)

github_core <- c(
  MuSiC  = "xuranw/MuSiC",             # bulk deconvolution, analysis/13_*
  presto = "immunogenomics/presto"     # fast Wilcoxon markers; NOT on CRAN
)

deferred <- c(
  "CellChat (Jin-s-Lab/CellChat) — heavy; install in its own project when needed",
  "monocle3 (cole-trapnell-lab/monocle3) — heavy; pseudotime, not needed for Project 1",
  "SeuratWrappers (satijalab/seurat-wrappers) — bridge to external tools",
  "SeuratDisk (mojaveazure/seurat-disk) — h5ad conversion; consider zellkonverter instead",
  "zellkonverter (Bioc) — better-maintained h5ad route than SeuratDisk"
)

# ----- install ---------------------------------------------------------------
# Order matters: BPCells and variancePartition first. They are the two most
# likely to fail on this machine, and both are load-bearing. Failing at minute
# ten beats failing at hour two.

step("1/5  BPCells  (C++17 + system HDF5 1.10.4)")
tryCatch({
  renv::install("BPCells")
  library(BPCells)
  ok("BPCells", as.character(packageVersion("BPCells")))
  cat("       expected 0.3.1 from r-universe; recorded SHA", substr(BPCELLS_SHA, 1, 8), "\n")
}, error = function(e) {
  warn("BPCells FAILED:", conditionMessage(e))
  warn("This blocks the on-disk single-cell workflow. Stopping.")
  quit(status = 1)
})

step("2/5  variancePartition  (no precedent on this machine)")
vp_ok <- tryCatch({
  renv::install("bioc::variancePartition")
  suppressPackageStartupMessages(library(variancePartition))
  ok("variancePartition", as.character(packageVersion("variancePartition")))
  TRUE
}, error = function(e) {
  warn("variancePartition FAILED:", conditionMessage(e))
  warn("Fallback: limma::duplicateCorrelation() for repeated measures.")
  warn("Set engine: limma_voom in config/de_runs/*.yml and note it in")
  warn("docs/decisions.md. Continuing with the rest of the install.")
  FALSE
})

step("3/5  CRAN core  (longest stage)")
renv::install(cran_core)
ok(length(cran_core), "CRAN packages")

step("4/5  Bioconductor core")
renv::install(paste0("bioc::", setdiff(bioc_core, "variancePartition")))
ok(length(bioc_core) - 1, "Bioconductor packages")

step("5/5  GitHub")
for (nm in names(github_core)) {
  tryCatch({
    renv::install(github_core[[nm]])
    ok(nm, as.character(packageVersion(nm)))
  }, error = function(e) warn(nm, "failed:", conditionMessage(e)))
}

# ----- verify ----------------------------------------------------------------
step("Verification")

check <- c("Seurat", "SeuratObject", "BPCells", "Matrix", "harmony",
           "DESeq2", "limma", "edgeR", "speckle", "SingleR",
           "SingleCellExperiment", "glmGamPoi", "MuSiC",
           if (vp_ok) "variancePartition")

expected <- c(Seurat = "5.5.0", SeuratObject = "5.4.0",
              BPCells = "0.3.1", Matrix = "1.7-5", harmony = "2.0.3")

for (p in check) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
  tag <- ""
  if (p %in% names(expected)) {
    tag <- if (identical(v, expected[[p]])) "  <- matches handoff objects" else
      sprintf("  <- EXPECTED %s, CHECK THIS", expected[[p]])
  }
  cat(sprintf("  %-22s %-12s%s\n", p, v %||% "MISSING", tag))
}

step("Functional test: Seurat 5 + BPCells on-disk")
tryCatch({
  suppressPackageStartupMessages({ library(Seurat); library(BPCells) })
  set.seed(1)
  m <- matrix(rpois(20000, 2), nrow = 500,
              dimnames = list(paste0("gene", 1:500), paste0("cell", 1:40)))
  d <- file.path(tempdir(), "bpc_verify")
  if (dir.exists(d)) unlink(d, recursive = TRUE)
  write_matrix_dir(as(m, "dgCMatrix"), dir = d)
  obj <- CreateSeuratObject(counts = open_matrix_dir(d))
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 100, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 5, verbose = FALSE)
  ok(sprintf("Seurat %s object from on-disk matrix, PCA ran (%d cells)",
             packageVersion("Seurat"), ncol(obj)))
  unlink(d, recursive = TRUE)
}, error = function(e) warn("functional test failed:", conditionMessage(e)))

# ----- snapshot --------------------------------------------------------------
step("renv::snapshot")
renv::snapshot(prompt = FALSE)
ok("renv.lock written")

lk <- jsonlite::fromJSON("renv.lock")
cat("  packages recorded:", length(lk$Packages), "\n")
cat("  R version        :", lk$R$Version, "\n")
cat("  repositories     :\n")
print(lk$R$Repositories)

# ----- deferred --------------------------------------------------------------
step("Deliberately NOT installed")
for (d in deferred) cat("  -", d, "\n")
cat("\n  Install these in their own project when needed, so their dependency\n")
cat("  trees cannot perturb this library. See docs/ENVIRONMENT.md.\n")

step("Done")
cat("
  NEXT

  1. Commit the lockfile — this is what makes freeze01 reproducible:
         git add renv.lock .Rprofile renv/settings.json
         git commit -m 'renv: initial snapshot (CRAN 2026-06-01, Bioc 3.23)'
         git push

  2. Point Positron at this R. In the remote window, the interpreter picker
     (top right of the Console pane) should offer R", as.character(getRversion()), "

  3. Sanity-check the environment from the repo:
         Rscript R/paths_check.R

  4. Then the first real analysis:
         jobs/run.sh analysis/00_manifest_and_batch/00_1_build_manifest.R

  If variancePartition failed, record the fallback in docs/decisions.md now.
")
