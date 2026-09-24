#!/usr/bin/env Rscript
# ==============================================================
# R/paths_check.R
#
# Run this once after cloning, and any time something mysterious
# breaks. It verifies that paths resolve, symlinks point somewhere
# real, required packages exist, and the config files parse.
#
#   Rscript R/paths_check.R
# ==============================================================

options(asset.quiet = TRUE)

ok   <- function(m) cat("  [ ok ]", m, "\n")
warn <- function(m) cat("  [warn]", m, "\n")
bad  <- function(m) cat("  [FAIL]", m, "\n")

n_fail <- 0L
n_warn <- 0L

cat("\n=== paths.R ===\n")
res <- tryCatch({ source("config/paths.R"); TRUE },
                error = function(e) { bad(conditionMessage(e)); FALSE })
if (!isTRUE(res)) {
  cat("\nCannot continue without config/paths.R. Are you in the repo root?\n")
  quit(status = 1L)
}
ok(paste("repo root:", REPO))

cat("\n=== directories ===\n")
dirs <- c(DATA = DATA, RESULTS = RESULTS, FIGURES = FIGURES, LOGS = LOGS,
          RAW = RAW, OBJ = OBJ, BPCELLS = BPCELLS, BULK = BULK,
          CLINICAL = CLINICAL, EXTERNAL = EXTERNAL, INHERIT = INHERIT,
          CONFIG = CONFIG, METADATA = METADATA, GENESETS = GENESETS,
          LABELS = LABELS, DOCS = DOCS)

for (nm in names(dirs)) {
  p <- dirs[[nm]]
  if (!file.exists(p)) {
    if (nm %in% c("DATA", "RESULTS", "FIGURES", "LOGS")) {
      bad(sprintf("%-9s %s  <- missing; create it or fix the symlink", nm, p))
      n_fail <- n_fail + 1L
    } else {
      warn(sprintf("%-9s %s  <- missing (fine if not needed yet)", nm, p))
      n_warn <- n_warn + 1L
    }
    next
  }
  extra <- ""
  li <- Sys.readlink(p)
  if (nzchar(li)) {
    extra <- paste0("  -> ", li)
    if (!file.exists(li)) {
      bad(sprintf("%-9s %s%s  <- BROKEN SYMLINK", nm, p, extra))
      n_fail <- n_fail + 1L
      next
    }
  }
  ok(sprintf("%-9s %s%s", nm, p, extra))
}

cat("\n=== writability ===\n")
for (nm in c("RESULTS", "FIGURES", "LOGS")) {
  p <- dirs[[nm]]
  if (!file.exists(p)) next
  tf <- file.path(p, paste0(".write_test_", Sys.getpid()))
  wrote <- tryCatch({ writeLines("x", tf); file.remove(tf); TRUE },
                    error = function(e) FALSE, warning = function(w) FALSE)
  if (wrote) ok(paste(nm, "writable")) else {
    bad(paste(nm, "NOT writable")); n_fail <- n_fail + 1L
  }
}

cat("\n=== clinical directory permissions ===\n")
if (dir.exists(CLINICAL)) {
  m <- file.mode(CLINICAL)
  if (as.integer(m) %% 64L != 0L) {
    warn(sprintf("%s is mode %s. chmod 700 it: patient-level data.",
                 CLINICAL, as.character(m)))
    n_warn <- n_warn + 1L
  } else {
    ok("clinical/ is owner-only")
  }
} else {
  warn("clinical/ does not exist yet")
  n_warn <- n_warn + 1L
}

cat("\n=== config files parse ===\n")
if (!requireNamespace("yaml", quietly = TRUE)) {
  bad("package 'yaml' not installed — provenance.R needs it")
  n_fail <- n_fail + 1L
} else {
  for (f in list.files(file.path(CONFIG, "freezes"), "\\.ya?ml$",
                       full.names = TRUE)) {
    r <- tryCatch({ yaml::read_yaml(f); TRUE },
                  error = function(e) { bad(paste(basename(f), "-",
                                                  conditionMessage(e))); FALSE })
    if (isTRUE(r)) ok(paste("freeze:", basename(f))) else n_fail <- n_fail + 1L
  }
  for (f in list.files(file.path(CONFIG, "cohorts"), "\\.ya?ml$",
                       full.names = TRUE)) {
    r <- tryCatch({ yaml::read_yaml(f); TRUE },
                  error = function(e) { bad(paste(basename(f), "-",
                                                  conditionMessage(e))); FALSE })
    if (isTRUE(r)) ok(paste("cohort:", basename(f))) else n_fail <- n_fail + 1L
  }
}

cat("\n=== metadata ===\n")
mf <- file.path(METADATA, "library_manifest.csv")
if (file.exists(mf)) {
  ok("library_manifest.csv present")
  src <- tryCatch({ source("R/io.R"); read_manifest(); TRUE },
                  error = function(e) { bad(conditionMessage(e)); FALSE })
  if (!isTRUE(src)) n_fail <- n_fail + 1L
} else {
  warn(paste("library_manifest.csv not present yet. Copy",
             "library_manifest_TEMPLATE.csv and fill it in."))
  n_warn <- n_warn + 1L
}
if (file.exists(file.path(METADATA, "celltype_dictionary.csv"))) {
  ok("celltype_dictionary.csv present")
}

cat("\n=== packages ===\n")
core <- c("yaml", "Matrix")
soft <- c("Seurat", "SeuratObject", "BPCells", "DESeq2", "limma", "edgeR",
          "variancePartition", "speckle", "ggplot2", "dplyr", "harmony",
          "SingleR", "apeglm")
for (p in core) {
  if (requireNamespace(p, quietly = TRUE)) ok(paste("core:", p)) else {
    bad(paste("core package missing:", p)); n_fail <- n_fail + 1L
  }
}
missing_soft <- soft[!vapply(soft, requireNamespace, logical(1),
                             quietly = TRUE)]
if (length(missing_soft)) {
  warn(paste("not installed (install as needed):",
             paste(missing_soft, collapse = ", ")))
} else {
  ok("all analysis packages present")
}

cat("\n=== git ===\n")
sha <- suppressWarnings(tryCatch(
  system2("git", c("rev-parse", "--short", "HEAD"),
          stdout = TRUE, stderr = FALSE),
  error = function(e) character(0)))
if (length(sha)) {
  ok(paste("HEAD:", sha[1]))
  dirty <- suppressWarnings(system2("git", c("status", "--porcelain"),
                                    stdout = TRUE, stderr = FALSE))
  if (length(dirty)) {
    warn(sprintf("%d uncommitted change(s). Commit before runs you keep.",
                 length(dirty)))
    n_warn <- n_warn + 1L
  } else ok("working tree clean")
} else {
  warn("not a git repository, or git unavailable")
  n_warn <- n_warn + 1L
}

cat("\n=== run context ===\n")
cat(sprintf("  freeze=%s  cohort=%s  labelset=%s\n", FREEZE, COHORT, LABELSET))
cat(sprintf("  example results dir: %s\n",
            results_dir("11_sc_pseudobulk_de", "example_run")))

cat("\n=== summary ===\n")
cat(sprintf("  %d failure(s), %d warning(s)\n\n", n_fail, n_warn))
if (n_fail > 0L) {
  cat("Fix the failures before running any analysis.\n\n")
  quit(status = 1L)
}
cat("Ready.\n\n")
