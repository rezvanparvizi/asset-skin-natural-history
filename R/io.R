# ==============================================================
# R/io.R
#
# Object loading, config loading, cohort subsetting, and BPCells
# path repair.
#
# The BPCells section is the most practically important part of this
# file. See docs/inherited_objects.md for why.
# ==============================================================

# ---- config loading ------------------------------------------

#' Read a freeze definition
read_freeze <- function(freeze = FREEZE) {
  f <- freeze_file(freeze)
  if (!file.exists(f)) stop("No such freeze: ", f, call. = FALSE)
  yaml::read_yaml(f)
}

#' Read a cohort definition
read_cohort <- function(cohort = COHORT) {
  f <- cohort_file(cohort)
  if (!file.exists(f)) stop("No such cohort: ", f, call. = FALSE)
  yaml::read_yaml(f)
}

#' Read the library manifest
#'
#' Technical + design fields, one row per library. Gitignored — lives on
#' the server only.
read_manifest <- function(path = file.path(METADATA, "library_manifest.csv")) {
  if (!file.exists(path)) {
    stop("Manifest not found at ", path, "\n",
         "Copy metadata/library_manifest_TEMPLATE.csv and fill it in. ",
         "The real manifest is gitignored on purpose.", call. = FALSE)
  }
  m <- utils::read.csv(path, stringsAsFactors = FALSE)
  req <- c("library_id", "sample_id", "subject_id", "group", "arm",
           "timepoint", "batch_id", "qc_status")
  miss <- setdiff(req, names(m))
  if (length(miss)) {
    stop("Manifest is missing required columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(m$library_id)) {
    stop("Duplicate library_id in manifest.", call. = FALSE)
  }
  m
}

#' Libraries included in a freeze
#'
#' Resolves the freeze's exclusion list against the manifest.
freeze_libraries <- function(freeze = FREEZE, manifest = read_manifest()) {
  fz <- read_freeze(freeze)
  excl <- vapply(fz$excluded_libraries %||% list(),
                 function(x) x$library_id, character(1))
  keep <- manifest[!manifest$library_id %in% excl, ]

  # honour merge_strategy for resequenced libraries
  if ("merge_strategy" %in% names(keep)) {
    keep <- keep[!keep$merge_strategy %in% c("exclude"), ]
  }
  if (!is.null(fz$n_libraries) && nrow(keep) != fz$n_libraries) {
    warning(sprintf(
      "Freeze %s declares n_libraries = %d but the manifest resolves to %d. ",
      freeze, fz$n_libraries, nrow(keep)), call. = FALSE)
  }
  keep
}

#' Samples belonging to a cohort
#'
#' Applies the cohort's include/exclude rules to the freeze's libraries.
#' Escape-therapy censoring needs the clinical table, which is NOT in git;
#' pass it explicitly when the cohort requires it.
cohort_samples <- function(cohort   = COHORT,
                           freeze   = FREEZE,
                           manifest = read_manifest(),
                           clinical = NULL) {

  ch  <- read_cohort(cohort)
  lib <- freeze_libraries(freeze, manifest)

  # Flagged libraries (freeze library_flags) stay in Layer 1 but are
  # never part of an inference cohort until the flag is resolved.
  if (isTRUE(ch$layer == 2) && "design_flag" %in% names(lib)) {
    flagged <- !is.na(lib$design_flag) & nzchar(lib$design_flag)
    if (any(flagged)) {
      message(sprintf("Cohort '%s': dropping %d flagged libraries (%s).",
                      cohort, sum(flagged),
                      paste(unique(lib$design_flag[flagged]), collapse = ", ")))
    }
    lib <- lib[!flagged, ]
  }

  inc <- ch$include %||% list()
  for (field in names(inc)) {
    if (!field %in% names(lib)) {
      stop("Cohort '", cohort, "' filters on '", field,
           "' which is not a manifest column.", call. = FALSE)
    }
    lib <- lib[lib[[field]] %in% unlist(inc[[field]]), ]
  }

  if (isTRUE(ch$exclude_post_escape)) {
    if (is.null(clinical)) {
      stop("Cohort '", cohort, "' censors post-escape-therapy samples, so it ",
           "needs the clinical table. Pass clinical = <data.frame with ",
           "subject_id, escape_start_month>.", call. = FALSE)
    }
    tp_month <- as.integer(sub("^M", "", lib$timepoint))
    esc <- clinical$escape_start_month[match(lib$subject_id,
                                             clinical$subject_id)]
    drop <- !is.na(esc) & tp_month >= esc
    if (any(drop)) {
      message(sprintf("Cohort '%s': censoring %d post-escape libraries.",
                      cohort, sum(drop)))
    }
    lib <- lib[!drop, ]
  }

  lib
}

# ---- labels --------------------------------------------------

#' Read the frozen barcode -> label table for a labelset
#'
#' This is the interface between Layer 1 (reference) and Layer 2
#' (inference). Layer 2 scripts should load labels this way rather than
#' relying on whatever metadata happens to be inside an .rds.
read_labels <- function(freeze = FREEZE, labelset = LABELSET) {
  f <- labels_file(freeze, labelset)
  if (!file.exists(f)) {
    stop("Label table not found: ", f, "\n",
         "Layer 1 (reference) must be run for this freeze before any ",
         "cohort-scoped analysis.", call. = FALSE)
  }
  utils::read.csv(f, stringsAsFactors = FALSE)
}

#' Attach frozen labels to a Seurat object by barcode
#'
#' Fails loudly on incomplete coverage rather than silently producing NAs.
attach_labels <- function(obj,
                          freeze   = FREEZE,
                          labelset = LABELSET,
                          cols     = c("lineage", "celltype", "subtype"),
                          require_complete = TRUE) {

  lab <- read_labels(freeze, labelset)
  if (!"barcode" %in% names(lab)) {
    stop("Label table has no 'barcode' column.", call. = FALSE)
  }
  i <- match(colnames(obj), lab$barcode)
  n_miss <- sum(is.na(i))
  if (n_miss > 0) {
    msg <- sprintf("%d of %d cells (%.1f%%) have no entry in labelset %s.",
                   n_miss, ncol(obj), 100 * n_miss / ncol(obj), labelset)
    if (require_complete) {
      stop(msg, "\nEither the object and the labelset come from different ",
           "freezes, or barcode suffixes differ. Check before proceeding.",
           call. = FALSE)
    }
    warning(msg, call. = FALSE)
  }
  for (cl in intersect(cols, names(lab))) {
    obj[[cl]] <- lab[[cl]][i]
  }
  obj[["labelset"]] <- labelset
  obj
}

# ==============================================================
# BPCells path repair
#
# A BPCells-backed Seurat object is a thin .rds holding an ABSOLUTE
# PATH pointer to an on-disk counts directory:
#
#   objects/ASSET_full_BPCells.rds            <- the object
#   objects/ASSET_full_BPCells_counts/        <- the actual data
#       col_names idxptr index_data index_starts row_names shape val version
#
# Consequences:
#   - Copy the .rds alone and it loads fine, then fails the moment you
#     touch the expression matrix.
#   - Copy both and the pointer still refers to the ORIGINAL absolute
#     path, i.e. your colleague's directory.
#   - If that directory is reorganized or purged, your object breaks
#     silently.
#
# So: always copy .rds AND _counts/ together, record the original path
# in docs/inherited_objects.md, and re-point with repoint_bpcells().
# ==============================================================

#' Report the on-disk paths a BPCells-backed object points to
bpcells_paths <- function(obj) {
  out <- list()
  for (a in names(obj@assays)) {
    for (lyr in c("counts", "data")) {
      m <- tryCatch(SeuratObject::LayerData(obj, assay = a, layer = lyr),
                    error = function(e) NULL)
      if (is.null(m)) next
      p <- tryCatch(attr(m, "dir"), error = function(e) NULL)
      if (is.null(p)) {
        # BPCells stores the path inside the iterable description
        p <- tryCatch(m@dir, error = function(e) NULL)
      }
      if (!is.null(p)) {
        out[[paste(a, lyr, sep = "/")]] <-
          list(path = p, exists = dir.exists(p))
      }
    }
  }
  if (!length(out)) message("No BPCells-backed layers detected.")
  out
}

#' Check that every BPCells path a object references actually exists
check_bpcells <- function(obj, stop_on_missing = TRUE) {
  p <- bpcells_paths(obj)
  if (!length(p)) return(invisible(TRUE))
  bad <- names(p)[!vapply(p, `[[`, logical(1), "exists")]
  for (nm in names(p)) {
    message(sprintf("  %-20s %s  [%s]", nm, p[[nm]]$path,
                    if (p[[nm]]$exists) "ok" else "MISSING"))
  }
  if (length(bad)) {
    msg <- paste0("BPCells directories missing for layer(s): ",
                  paste(bad, collapse = ", "),
                  "\nUse repoint_bpcells() after copying the _counts/ ",
                  "directories to your own storage.")
    if (stop_on_missing) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(length(bad) == 0)
}

#' Re-point a BPCells-backed object at counts directories under your own storage
#'
#' @param obj A Seurat object loaded from a colleague's .rds.
#' @param new_dir Directory containing the copied *_counts/ folders.
#' @param assay,layer Which layer to re-point.
#'
#' Requires the BPCells package. Read the counts back with
#' BPCells::open_matrix_dir() and reassign, rather than editing slots.
repoint_bpcells <- function(obj, new_dir, assay = "RNA", layer = "counts") {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    stop("Package 'BPCells' is required to re-point on-disk matrices.",
         call. = FALSE)
  }
  if (!dir.exists(new_dir)) {
    stop("new_dir does not exist: ", new_dir, call. = FALSE)
  }
  m <- BPCells::open_matrix_dir(dir = new_dir)
  old_cells <- colnames(obj)
  if (!all(old_cells %in% colnames(m))) {
    stop("The matrix at ", new_dir, " does not contain all cells in the ",
         "object (", sum(!old_cells %in% colnames(m)), " missing). ",
         "Wrong counts directory?", call. = FALSE)
  }
  m <- m[, old_cells]
  SeuratObject::LayerData(obj, assay = assay, layer = layer) <- m
  message("Re-pointed ", assay, "/", layer, " -> ", new_dir)
  obj
}

#' Load an object and immediately verify its on-disk backing
#'
#' Use this instead of readRDS() for anything BPCells-backed.
load_object <- function(path, check = TRUE) {
  if (!file.exists(path)) stop("No such object: ", path, call. = FALSE)
  obj <- readRDS(path)
  message("Loaded ", basename(path), ": ",
          ncol(obj), " cells x ", nrow(obj), " features")
  if (check) check_bpcells(obj, stop_on_missing = FALSE)
  obj
}

# ---- small utilities -----------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Defined in R/provenance.R. Fallback here so io.R also works standalone.
if (!exists(".rel_to_repo", mode = "function")) {
  .rel_to_repo <- function(p, repo = REPO) {
    p <- tryCatch(normalizePath(p, mustWork = FALSE), error = function(e) p)
    if (startsWith(p, repo)) sub(paste0("^", repo, "/?"), "", p) else p
  }
}

#' Write a table into a run's tables/ directory
save_table <- function(x, run, name) {
  stopifnot(inherits(run, "asset_run"))
  f <- file.path(run$tables, paste0(name, ".csv"))
  utils::write.csv(x, f, row.names = FALSE)
  message("  wrote ", .rel_to_repo(f))
  invisible(f)
}

#' Write a table that carries patient-level identifiers
#'
#' Anything with subject IDs, clinical values or barcode-level labels goes
#' under data/patient_level/, never results/ (docs/decisions.md,
#' 2026-09-29). The run directory gets a pointer file instead.
save_patient_table <- function(x, run, name) {
  stopifnot(inherits(run, "asset_run"))
  d <- file.path(PATIENT, run$context$freeze, run$context$scope,
                 run$stage, run$run_name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  f <- file.path(d, paste0(name, ".csv"))
  utils::write.csv(x, f, row.names = FALSE, na = "")
  Sys.chmod(f, "0600")
  cat(.rel_to_repo(f), "\n", file = file.path(run$tables, "PATIENT_LEVEL_TABLES.txt"),
      append = TRUE)
  message("  wrote ", .rel_to_repo(f), " (patient-level)")
  invisible(f)
}

#' Save a ggplot into a run's plots/ (or plots/QC/) directory
save_plot <- function(p, run, name, qc = FALSE,
                      width = 7, height = 5, device = "pdf") {
  stopifnot(inherits(run, "asset_run"))
  dest <- if (qc) run$qc else run$plots
  f <- file.path(dest, paste0(name, ".", device))
  ggplot2::ggsave(f, p, width = width, height = height, device = device)
  message("  wrote ", .rel_to_repo(f))
  invisible(f)
}
