# ==============================================================
# 00_0_import_raw_libraries.R
#
# Copy every library's CellRanger output into data/raw/ and generate
# metadata/library_manifest.csv from the source sample map.
#
# Why copy rather than symlink: the upstream files belong to other
# people and have already changed once (wasikowr's seurat.RDS was
# rewritten 2026-09-02, after the colleague's objects were built from
# it). A freeze must not move under us. Copies are verified by md5 and
# made read-only. See docs/decisions.md, 2026-09-29.
#
# What is copied, per library (from per_sample_outs/<library>/):
#   sample_filtered_feature_bc_matrix.h5   cells called by CellRanger
#   sample_raw_feature_bc_matrix.h5        incl. empty droplets -> SoupX
#   metrics_summary.csv                    CellRanger QC metrics
# and per pool (one Flex capture = one GEM reaction):
#   raw_feature_bc_matrix.h5               all barcodes, pool level
#   config.csv                             probe set, reference, BC map
#
# Layout:  data/raw/<run_id>/<library_id>/...   data/raw/<run_id>/_pool/...
# Raw data is freeze-independent: freeze02 adds directories, never
# replaces them.
#
# Idempotent: a file already present with a matching md5 is skipped. A
# file present with a DIFFERENT md5 stops the script — raw data is never
# silently overwritten.
#
#   [ASSET_REPLACE_MANIFEST=1] jobs/run.sh analysis/00_manifest_and_batch/00_0_import_raw_libraries.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")

fz  <- read_freeze()
src <- fz$sources
if (is.null(src$sample_map) || !file.exists(src$sample_map)) {
  stop("freeze ", FREEZE, " has no readable sources$sample_map", call. = FALSE)
}

run <- init_run(
  stage    = "00_manifest_and_batch",
  run_name = "import_raw_libraries",
  params   = list(sample_map = src$sample_map),
  notes    = "Copy CellRanger per-sample and pool outputs; generate manifest"
)

# ---- 1. sample map -------------------------------------------

map <- utils::read.delim(src$sample_map, stringsAsFactors = FALSE,
                         check.names = FALSE)
need <- c("sampleID", "path", "condition", "time", "patient", "treatment")
miss <- setdiff(need, names(map))
if (length(miss)) stop("sample map lacks columns: ", paste(miss, collapse = ", "))
if (anyDuplicated(map$sampleID)) stop("duplicate sampleID in sample map")
message(sprintf("Sample map: %d libraries", nrow(map)))

# Subject ID corrections, confirmed by the owner, applied before anything
# else uses map$patient (metadata/design/subject_id_corrections.csv,
# gitignored; docs/decisions.md 2026-10-01). A correction whose
# sample-map ID is absent is reported, not fatal: a later sample map may
# already be fixed upstream.
ID_FIX <- file.path(DESIGN, "subject_id_corrections.csv")
if (file.exists(ID_FIX)) {
  fix <- utils::read.csv(ID_FIX, stringsAsFactors = FALSE, colClasses = "character")
  stopifnot(all(c("sample_map_id", "Subject_ID") %in% names(fix)),
            !anyDuplicated(fix$sample_map_id))
  hit <- match(map$patient, fix$sample_map_id)
  map$patient[!is.na(hit)] <- fix$Subject_ID[hit[!is.na(hit)]]
  message(sprintf("Subject ID corrections: %d of %d applied (%d libraries)",
                  length(unique(stats::na.omit(hit))), nrow(fix), sum(!is.na(hit))))
} else {
  message("No subject ID corrections file at ", .rel_to_repo(ID_FIX))
}

# per_sample_outs/<lib>/count/sample_filtered_feature_bc_matrix.h5
map$sample_dir <- dirname(dirname(map$path))
map$pool_dir   <- sub("/per_sample_outs/.*$", "", map$path)
stopifnot(basename(map$sample_dir) == map$sampleID)

# run_id: the sequencing-core submission, e.g. 13452-JF (8911_JF -> 8911-JF)
map$run_id <- sub("_", "-", regmatches(map$pool_dir,
                                       regexpr("[0-9]{4,5}[-_][A-Z]{2}",
                                               map$pool_dir)))
stopifnot(length(map$run_id) == nrow(map),
          startsWith(map$sampleID, map$run_id))

# ---- 2. copy with md5 verification ---------------------------

RAW_CHECKSUMS <- file.path(RAW, "checksums.csv")
done <- if (file.exists(RAW_CHECKSUMS)) {
  utils::read.csv(RAW_CHECKSUMS, stringsAsFactors = FALSE)
} else NULL

copy_verified <- function(from, to) {
  if (!file.exists(from)) stop("source missing: ", from, call. = FALSE)
  md5_src <- unname(tools::md5sum(from))
  # Verified copies are made read-only. A writable file at the
  # destination is an interrupted copy from an earlier attempt: redo it.
  if (file.exists(to) && file.access(to, 2) == 0) file.remove(to)
  if (file.exists(to)) {
    md5_dst <- unname(tools::md5sum(to))
    if (md5_dst != md5_src) {
      stop("destination exists with a different md5; refusing to overwrite:\n  ",
           to, "\nsource: ", from, call. = FALSE)
    }
    action <- "present"
  } else {
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    part <- paste0(to, ".partial")
    if (!file.copy(from, part, overwrite = TRUE, copy.date = TRUE)) {
      stop("copy failed: ", from, call. = FALSE)
    }
    if (unname(tools::md5sum(part)) != md5_src) {
      file.remove(part)
      stop("md5 mismatch after copy: ", from, call. = FALSE)
    }
    file.rename(part, to)
    Sys.chmod(to, "0444")
    action <- "copied"
  }
  data.frame(dest = .rel_to_repo(to), source = from,
             bytes = file.size(to), md5 = md5_src, action = action,
             stringsAsFactors = FALSE)
}

lib_files <- c(filtered = "count/sample_filtered_feature_bc_matrix.h5",
               raw      = "count/sample_raw_feature_bc_matrix.h5",
               metrics  = "metrics_summary.csv")

log_rows <- list()
for (i in seq_len(nrow(map))) {
  for (f in lib_files) {
    log_rows[[length(log_rows) + 1]] <- copy_verified(
      file.path(map$sample_dir[i], f),
      file.path(RAW, map$run_id[i], map$sampleID[i], basename(f)))
  }
  if (i %% 20 == 0) message(sprintf("  libraries: %d / %d", i, nrow(map)))
}

pools <- unique(map[, c("run_id", "pool_dir")])
stopifnot(!anyDuplicated(pools$run_id))
for (i in seq_len(nrow(pools))) {
  for (f in c("multi/count/raw_feature_bc_matrix.h5", "config.csv")) {
    log_rows[[length(log_rows) + 1]] <- copy_verified(
      file.path(pools$pool_dir[i], f),
      file.path(RAW, pools$run_id[i], "_pool", basename(f)))
  }
}
copy_log <- do.call(rbind, log_rows)
save_table(copy_log, run, "copy_log")

ck_cols <- c("dest", "source", "bytes", "md5")
checksums <- unique(rbind(if (!is.null(done)) done[, ck_cols],
                          copy_log[, ck_cols]))
utils::write.csv(checksums, RAW_CHECKSUMS, row.names = FALSE)
message(sprintf("Copied %d files, %d already present, %.1f GB total",
                sum(copy_log$action == "copied"),
                sum(copy_log$action == "present"),
                sum(copy_log$bytes) / 1e9))

# a copy of the sample map itself, as received
map_copy <- copy_verified(src$sample_map,
                          file.path(RAW, "_sample_maps",
                                    paste0(FREEZE, "_", basename(src$sample_map))))

# ---- 3. pool configs: probe set, reference, unmapped libraries ----

# config.csv comes in two dialects: the core's, and wasikowr's control
# re-runs (trailing commas, a force_cells column). Normalise both.
read_pool_config <- function(f) {
  x <- trimws(sub(",+\\s*$", "", readLines(f, warn = FALSE)))
  sec <- cumsum(grepl("^\\[", x))
  s_start <- which(x == "[samples]")
  rows <- x[sec == sec[s_start]][-1]
  rows <- rows[nzchar(rows)]
  hdr <- trimws(strsplit(rows[1], ",")[[1]])
  samples <- sub(",.*$", "", rows[-1])
  ref <- grep("^reference,", x, value = TRUE)
  list(reference = if (length(ref)) basename(sub("^[^,]*,", "", ref[1])) else NA,
       samples = samples, force_cells = "force_cells" %in% hdr)
}

# Panel name from the probe_set.csv CellRanger copies into each library's
# outs, not from the config path (the control re-runs used a renamed copy).
probe_panel <- function(sample_dir) {
  f <- file.path(sample_dir, "count", "probe_set.csv")
  if (!file.exists(f)) return(NA_character_)
  h <- grep("^#panel_name=", readLines(f, n = 10), value = TRUE)
  if (length(h)) sub("^#panel_name=", "", h[1]) else NA_character_
}

cellranger_version <- function(sample_dir) {
  h <- file.path(sample_dir, "web_summary.html")
  if (!file.exists(h)) return(NA_character_)
  txt <- readLines(h, warn = FALSE)
  v <- unique(regmatches(txt, regexpr("cellranger-[0-9]+\\.[0-9]+\\.[0-9]+", txt)))
  if (length(v)) paste(v, collapse = ";") else NA_character_
}

pool_info <- do.call(rbind, lapply(seq_len(nrow(pools)), function(i) {
  cfg <- read_pool_config(file.path(RAW, pools$run_id[i], "_pool", "config.csv"))
  first_lib <- map$sample_dir[map$run_id == pools$run_id[i]][1]
  unmapped <- setdiff(cfg$samples, map$sampleID)
  data.frame(run_id = pools$run_id[i],
             n_in_pool = length(cfg$samples),
             n_in_sample_map = sum(map$run_id == pools$run_id[i]),
             not_in_sample_map = paste(unmapped, collapse = ";"),
             probe_set = probe_panel(first_lib), reference = cfg$reference,
             force_cells = cfg$force_cells,
             cellranger = cellranger_version(first_lib),
             stringsAsFactors = FALSE)
}))
save_table(pool_info, run, "pool_info")

# ---- 4. CellRanger metrics per library (technical, no patient data) ----

read_metrics <- function(f) {
  m <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  m <- m[m$Category == "Cells" & m[["Library Type"]] == "Gene Expression", ]
  num <- function(name) {
    v <- m[["Metric Value"]][m[["Metric Name"]] == name]
    if (!length(v)) return(NA_real_)
    as.numeric(gsub("[,%]", "", v[1]))
  }
  c(cells = num("Cells"), median_umi = num("Median UMI counts per cell"),
    median_genes = num("Median genes per cell"),
    mean_reads = num("Mean reads per cell"),
    pct_mapped_in_cells = num("Confidently mapped reads in cells"))
}
metrics <- do.call(rbind, lapply(seq_len(nrow(map)), function(i) {
  f <- file.path(RAW, map$run_id[i], map$sampleID[i], "metrics_summary.csv")
  data.frame(library_id = map$sampleID[i], t(read_metrics(f)))
}))
save_table(metrics, run, "cellranger_metrics")

# ---- 5. manifest ---------------------------------------------

tp <- c("Baseline" = "M00", "Month 3" = "M03", "Month 6" = "M06",
        "Control" = NA)
if (!all(map$time %in% names(tp))) {
  stop("unexpected time values: ",
       paste(setdiff(unique(map$time), names(tp)), collapse = ", "))
}
is_hc <- map$time == "Control"
stopifnot(all(is_hc == (map$treatment == "Control")),
          all(is_hc == (map$patient == "Control")))

# The leading YYMMDD of the sample name (e.g. 250422_ASSE_Skin_12).
# Its meaning (fixation? hybridisation? capture?) is NOT confirmed, so it
# is stored under a neutral name and tested as a batch candidate only.
name_date <- ifelse(grepl("^[0-9]{6}_", map$condition),
                    format(as.Date(substr(map$condition, 1, 6), "%y%m%d")),
                    NA_character_)

pinfo <- pool_info[match(map$run_id, pool_info$run_id), ]

# Arm and group values are spelled as they appear in publication figures:
# Placebo / Abatacept, SSc / HC. Never all caps. Anything else maps to NA
# and fails the assertion below.
arm_label <- function(x) {
  lab <- c(placebo = "Placebo", abatacept = "Abatacept")
  unname(lab[tolower(trimws(x))])
}

man <- data.frame(
  library_id        = map$sampleID,
  sample_id         = map$sampleID,     # one library per biopsy in freeze01
  Subject_ID        = ifelse(is_hc, map$sampleID, map$patient),
  group             = ifelse(is_hc, "HC", "SSc"),
  arm               = ifelse(is_hc, NA, arm_label(map$treatment)),
  timepoint         = unname(tp[map$time]),
  batch_id          = map$run_id,       # Flex pool = one GEM capture
  capture_date      = NA,               # unknown; see sample_name_date
  run_id            = map$run_id,
  cellranger_run    = pinfo$cellranger,
  qc_status         = "pending",
  resequenced_of    = NA,
  merge_strategy    = NA,
  notes             = NA,
  sample_name_date  = name_date,
  probe_set         = pinfo$probe_set,
  design_flag       = NA_character_,
  source_orig_ident = map$condition,    # matches orig.ident in wasikowr's object
  raw_dir           = file.path(map$run_id, map$sampleID),
  stringsAsFactors  = FALSE
)

# flags declared in the freeze
for (fl in fz$library_flags %||% list()) {
  i <- match(fl$library_id, man$library_id)
  if (is.na(i)) stop("library_flags names an unknown library: ", fl$library_id)
  man$design_flag[i] <- fl$flag
}

fail <- function(msg) {
  finalize_run(run, status = "failed", verdict = substr(msg, 1, 200))
  stop(msg, call. = FALSE)
}

# assertions
stopifnot(!anyDuplicated(man$library_id),
          all(man$arm[!is_hc] %in% c("Placebo", "Abatacept")))
# a subject may appear twice at one timepoint only if every library
# involved is flagged in the freeze
visit <- paste(man$Subject_ID, man$timepoint)[!is_hc]
dup <- visit %in% visit[duplicated(visit)]
unflagged_dup <- dup & is.na(man$design_flag[!is_hc])
if (any(unflagged_dup)) {
  fail(paste("subject x timepoint not unique and not flagged for libraries:",
             paste(man$library_id[!is_hc][unflagged_dup], collapse = ", ")))
}
arm_per_subject <- tapply(man$arm[!is_hc], man$Subject_ID[!is_hc],
                          function(a) length(unique(a)))
if (any(arm_per_subject > 1)) fail("subject assigned to >1 arm")

# The manifest carries subject IDs: it is patient-level and never goes
# in results/.
mf <- file.path(METADATA, "library_manifest.csv")
tmpl <- file.path(METADATA, "library_manifest_TEMPLATE.csv")
new_file <- save_patient_table(man, run, "library_manifest_generated")

replaceable <- !file.exists(mf) ||
  identical(unname(tools::md5sum(mf)), unname(tools::md5sum(tmpl))) ||
  identical(unname(tools::md5sum(mf)), unname(tools::md5sum(new_file)))
# A deliberate regeneration (e.g. a nomenclature change) is requested
# with ASSET_REPLACE_MANIFEST=1. The old manifest is kept next to the new
# generated one, never discarded.
if (!replaceable && identical(Sys.getenv("ASSET_REPLACE_MANIFEST"), "1")) {
  prev <- file.path(dirname(new_file), "library_manifest_previous.csv")
  file.copy(mf, prev, overwrite = FALSE)
  Sys.chmod(prev, "0600")
  message("ASSET_REPLACE_MANIFEST=1: previous manifest kept at ", .rel_to_repo(prev))
  replaceable <- TRUE
}
if (!replaceable) {
  finalize_run(run, status = "failed",
               verdict = "existing manifest differs from the generated one; not overwritten")
  stop("metadata/library_manifest.csv exists and differs from the generated ",
       "manifest. Compare with\n  ", new_file,
       "\nand replace it by hand if the generated one is right.", call. = FALSE)
}
file.copy(new_file, mf, overwrite = TRUE)
Sys.chmod(mf, "0600")
message("wrote ", .rel_to_repo(mf))

# ---- 6. sample sheets for clinical-data preparation -------------
# One row per library and one per subject, for the owner to build the
# clinical table against. Patient-level, so data/patient_level/ only.
# The QC note is PRELIMINARY — CellRanger metrics only; formal QC is
# stage 03.

qc_note <- with(metrics, trimws(paste(
  ifelse(median_genes < 300, sprintf("median genes/cell %g;", median_genes), ""),
  ifelse(cells > 30000, sprintf("%g cells called (cell calling likely failed);", cells), ""),
  ifelse(pct_mapped_in_cells < 40, sprintf("%g%% reads mapped in cells;", pct_mapped_in_cells), ""))))
qc_note[!nzchar(qc_note)] <- NA

sheet <- merge(man[, c("library_id", "Subject_ID", "group", "arm", "timepoint",
                       "design_flag", "batch_id", "sample_name_date",
                       "cellranger_run", "source_orig_ident")],
               cbind(metrics, qc_note_preliminary = qc_note),
               by = "library_id", all.x = TRUE, sort = FALSE)
sheet$force_cells <- pool_info$force_cells[match(sheet$batch_id, pool_info$run_id)]
sheet <- sheet[order(sheet$group, sheet$Subject_ID, sheet$timepoint), ]
stopifnot(nrow(sheet) == nrow(man))
save_patient_table(sheet, run, "sample_sheet_per_library")

ssc_sheet <- sheet[sheet$group == "SSc", ]
subj <- do.call(rbind, lapply(split(ssc_sheet, ssc_sheet$Subject_ID), function(d)
  data.frame(Subject_ID = d$Subject_ID[1], arm = d$arm[1],
             n_libraries = nrow(d),
             has_M00 = "M00" %in% d$timepoint, has_M03 = "M03" %in% d$timepoint,
             has_M06 = "M06" %in% d$timepoint,
             libraries = paste(d$library_id, collapse = ";"),
             pools = paste(sort(unique(d$batch_id)), collapse = ";"),
             flags = paste(stats::na.omit(unique(d$design_flag)), collapse = ";"),
             qc_notes = paste(stats::na.omit(d$qc_note_preliminary), collapse = " | "),
             stringsAsFactors = FALSE)))
save_patient_table(subj, run, "sample_sheet_per_subject")

# ---- 7. verdict ------------------------------------------------

n_unmapped <- sum(nzchar(pool_info$not_in_sample_map) &
                  pool_info$n_in_sample_map < pool_info$n_in_pool)
versions <- unique(c(pool_info$cellranger, pool_info$probe_set))
summary_lines <- c(
  sprintf("Libraries: %d (%d SSc from %d subjects, %d HC)", nrow(man),
          sum(!is_hc), length(unique(man$Subject_ID[!is_hc])), sum(is_hc)),
  sprintf("Pools: %d", nrow(pools)),
  sprintf("Pools with libraries absent from the sample map: %d", n_unmapped),
  paste("CellRanger / probe set:", paste(versions, collapse = " | ")))
writeLines(summary_lines, file.path(run$dir, "SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

finalize_run(run, status = "ok", verdict = paste(
  sprintf("%d libraries (%d SSc / %d subjects, %d HC) in %d pools copied and md5-verified;",
          nrow(man), sum(!is_hc), length(unique(man$Subject_ID[!is_hc])),
          sum(is_hc), nrow(pools)),
  "manifest generated; batch_id = pool run_id"))
