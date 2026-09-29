# ==============================================================
# config/paths.R
#
# The ONE place where filesystem paths are defined. Every analysis
# script starts with source("config/paths.R"). No script should ever
# contain a hard-coded absolute path.
#
# Usage:
#   source("config/paths.R")
#   readRDS(file.path(OBJ, "ASSET_full_lineage.rds"))
#
# Override any path with an environment variable of the same name,
# which is what jobs/run.sh does. Example:
#   ASSET_DATA=/other/mount Rscript analysis/03_sc_qc/03_1_qc.R
# ==============================================================

# ---- repo root -----------------------------------------------
# Resolved from this file's own location so it works no matter where
# R was launched from (repo root, a subdirectory, or Positron's console).
.find_repo_root <- function() {
  # 1. explicit override
  env <- Sys.getenv("ASSET_REPO", unset = NA_character_)
  if (!is.na(env) && dir.exists(env)) return(normalizePath(env))

  # 2. walk up from the working directory looking for a marker
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (file.exists(file.path(d, "config", "paths.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }

  stop("Could not locate the repo root. Either run R from inside the ",
       "repository, or set ASSET_REPO=/path/to/asset-skin-natural-history",
       call. = FALSE)
}

REPO <- .find_repo_root()

# ---- helper: env var with a default --------------------------
.p <- function(var, default) {
  v <- Sys.getenv(var, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) default else v
}

# ==============================================================
# EDIT THIS BLOCK ONCE, ON FIRST SETUP
#
# data/ results/ figures/ logs/ are symlinks in the repo pointing at
# large storage. If you would rather not use symlinks, set the four
# paths below to absolute locations instead.
# ==============================================================

DATA    <- .p("ASSET_DATA",    file.path(REPO, "data"))
RESULTS <- .p("ASSET_RESULTS", file.path(REPO, "results"))
FIGURES <- .p("ASSET_FIGURES", file.path(REPO, "figures"))
LOGS    <- .p("ASSET_LOGS",    file.path(REPO, "logs"))

# ---- inside data/ --------------------------------------------
RAW      <- file.path(DATA, "raw")          # CellRanger outs, per library
OBJ      <- file.path(DATA, "objects")      # your Seurat / BPCells objects
BPCELLS  <- file.path(DATA, "bpcells")      # on-disk count matrices
BULK     <- file.path(DATA, "bulk")         # bulk skin + PBMC counts
CLINICAL <- file.path(DATA, "clinical")     # DCC clinical tables — NEVER in git
EXTERNAL <- file.path(DATA, "external")     # reference atlases for label transfer
INHERIT  <- file.path(DATA, "inherited")    # objects copied from colleagues
PATIENT  <- file.path(DATA, "patient_level") # run outputs carrying patient IDs — NEVER results/

# ---- inside the repo (tracked) -------------------------------
CONFIG   <- file.path(REPO, "config")
METADATA <- file.path(REPO, "metadata")
GENESETS <- file.path(METADATA, "gene_sets")
LABELS   <- file.path(METADATA, "labels")
DOCS     <- file.path(REPO, "docs")
RFUN     <- file.path(REPO, "R")
PIPELINE <- file.path(REPO, "pipelines")

# ---- run context ---------------------------------------------
# Set by jobs/run.sh, or set by hand in an interactive session.
# init_run() reads these and refuses to run if a required one is missing.
FREEZE   <- .p("ASSET_FREEZE",   "freeze01")
COHORT   <- .p("ASSET_COHORT",   "reference")
LABELSET <- .p("ASSET_LABELSET", "labelset01")

# ---- convenience constructors --------------------------------

#' Path to a freeze definition file
freeze_file <- function(freeze = FREEZE) {
  file.path(CONFIG, "freezes", paste0(freeze, ".yml"))
}

#' Path to a cohort definition file
cohort_file <- function(cohort = COHORT) {
  file.path(CONFIG, "cohorts", paste0(cohort, ".yml"))
}

#' Path to a labelset's barcode -> label table
labels_file <- function(freeze = FREEZE, labelset = LABELSET) {
  file.path(LABELS, sprintf("%s_%s_celltype_labels.csv.gz", freeze, labelset))
}

#' Results directory for a stage, optionally a named run inside it
#'
#' results/{freeze}/{cohort}/{stage}/{run_name}/
results_dir <- function(stage, run_name = NULL,
                        freeze = FREEZE, cohort = COHORT) {
  p <- file.path(RESULTS, freeze, cohort, stage)
  if (!is.null(run_name)) p <- file.path(p, run_name)
  p
}

# ---- print a summary when sourced interactively --------------
if (interactive() && !isTRUE(getOption("asset.quiet"))) {
  message(sprintf(
    "ASSET paths loaded\n  repo:   %s\n  data:   %s\n  results:%s\n  context: freeze=%s cohort=%s labelset=%s",
    REPO, DATA, RESULTS, FREEZE, COHORT, LABELSET))
}
