# ==============================================================
# 00_7_clean_bulk_metadata_workbooks.R
#
# Corrects the bulk sample workbooks IN PLACE, after a backup (owner
# request 2026-10-02; decisions.md "scripts-only access ... in-place
# Excel edits with backup"). Inconsistencies come from 00_5 / 00_6.
#
# Changes (skin workbook, unless stated):
#   1. Sex "F"/"M" -> "Female"/"Male"                    (master coding)
#   2. Improver, True_improver "Non_Improver" -> "Non-Improver"
#   3. mRSS_category <- the master's value for the subject (the master is
#      the authoritative source of the owner's outcome grouping)
#   4. Sample_ID typo: the one skin Sample_ID missing from the matrices vs
#      the one unexplained matrix column. Whichever equals the ID built
#      from the row's own Subject_ID and Timepoint
#      ("<Subject_ID with _>_<Timepoint>", "_A" for replicate 2) is right.
#      Metadata wrong -> metadata corrected (and the PBMC
#      Matched_Skin_Sample_ID that points to it). Matrix column wrong ->
#      metadata left alone; the column rename is written to
#      data/clinical/bulk_skin_matrix_column_fixes.csv for the import step.
#      Neither -> stop for the owner.
# NOT touched (owner): autoantibody columns, escape columns.
#
# Only the changed columns are rewritten (openxlsx, formatting kept);
# afterwards every other cell is re-read and must be unchanged.
# Counts to tables/; ID-level change log to data/patient_level/.
#
#   jobs/run.sh analysis/00_manifest_and_batch/00_7_clean_bulk_metadata_workbooks.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(openxlsx)
})

P <- list(skin_file = "ASSET_clinical_data_skin_bulk_metadata.xlsx",
          pbmc_file = "ASSET_clinical_data_PBMC_bulk_metadata.xlsx",
          sex_map = c(F = "Female", M = "Male"),
          improver_map = c(Non_Improver = "Non-Improver"))

run <- init_run(
  stage    = "00_manifest_and_batch",
  run_name = "clean_bulk_metadata_workbooks",
  params   = P,
  notes    = "In-place corrections of the bulk workbooks (coding, mRSS_category from master, Sample_ID typo)"
)

# ---- backup ------------------------------------------------------
bk <- file.path(CLINICAL, "backup", run$run_id)
dir.create(bk, recursive = TRUE, showWarnings = FALSE, mode = "0700")
for (f in c(P$skin_file, P$pbmc_file, basename(CLINICAL_MASTER))) {
  stopifnot(file.copy(file.path(CLINICAL, f), file.path(bk, f), copy.date = TRUE))
  Sys.chmod(file.path(bk, f), "0400")
  stopifnot(unname(tools::md5sum(file.path(bk, f))) == unname(tools::md5sum(file.path(CLINICAL, f))))
}
message("  backup: ", .rel_to_repo(bk))

rd <- function(f, sheet = 1) as.data.table(read_excel(file.path(CLINICAL, f), sheet = sheet, col_types = "text"))
master <- rd(basename(CLINICAL_MASTER)); skin0 <- rd(P$skin_file); pbmc0 <- rd(P$pbmc_file)
skin <- copy(skin0); pbmc <- copy(pbmc0)
changes <- list(); counts <- list()
log_change <- function(tab, id, col, old, new) {
  ch <- (is.na(old) != is.na(new)) | (!is.na(old) & !is.na(new) & old != new)
  if (any(ch)) changes[[length(changes) + 1]] <<- data.table(table = tab, row_id = id[ch], column = col,
                                                             old = old[ch], new = new[ch])
  counts[[length(counts) + 1]] <<- data.table(table = tab, column = col, rows_changed = sum(ch))
}

# ---- 1-2. coding -------------------------------------------------
stopifnot(all(na.omit(skin$Sex) %in% c(names(P$sex_map), P$sex_map)))
new <- ifelse(skin$Sex %in% names(P$sex_map), P$sex_map[skin$Sex], skin$Sex)
log_change("skin", skin$Sample_ID, "Sex", skin$Sex, new); skin$Sex <- new
for (col in c("Improver", "True_improver")) {
  new <- ifelse(skin[[col]] %in% names(P$improver_map), P$improver_map[skin[[col]]], skin[[col]])
  log_change("skin", skin$Sample_ID, col, skin[[col]], new); skin[[col]] <- new
}

# ---- 3. mRSS_category from the master -----------------------------
mc <- master$mRSS_category[match(skin$Subject_ID, master$Subject_ID)]
stopifnot(!anyNA(mc))
transitions <- data.table(skin_value = skin$mRSS_category, master_value = mc)[
  is.na(skin_value) | skin_value != master_value, .N, by = .(skin_value, master_value)]
log_change("skin", skin$Sample_ID, "mRSS_category", skin$mRSS_category, mc); skin$mRSS_category <- mc

# ---- 4. Sample_ID typo ---------------------------------------------
mat <- list.files(file.path(BULK, "skin"), pattern = "raw_counts.*EMSEMBL", full.names = TRUE)[1]
cols <- names(fread(mat, nrows = 0))[-1]
expected <- paste0(gsub("-", "_", skin$Subject_ID), "_", skin$Timepoint,
                   ifelse(skin$replicate == "2", "_A", ""))
add_count <- function(item, n) counts[[length(counts) + 1]] <<- data.table(table = "skin", column = item, rows_changed = n)
add_count("check: Sample_ID equals ID built from Subject_ID + Timepoint (+_A), all rows", sum(skin$Sample_ID == expected))
if (mean(skin$Sample_ID == expected) < 0.95) {
  finalize_run(run, status = "failed", verdict = "Sample_ID construction rule holds for < 95% of rows; cannot arbitrate the typo")
  stop("Sample_ID rule does not hold generally (", sum(skin$Sample_ID == expected), "/", nrow(skin), ")", call. = FALSE)
}
miss <- which(!skin$Sample_ID %in% cols); extra <- setdiff(cols, skin$Sample_ID)
id_fix <- data.table()
if (length(miss) == 1 && length(extra) == 1) {
  meta_ok <- skin$Sample_ID[miss] == expected[miss]; col_ok <- extra == expected[miss]
  if (!meta_ok && col_ok) {
    old <- skin$Sample_ID[miss]
    log_change("skin", skin$Sample_ID, "Sample_ID", skin$Sample_ID, replace(skin$Sample_ID, miss, extra))
    skin$Sample_ID[miss] <- extra
    i <- which(pbmc$Matched_Skin_Sample_ID == old)
    new <- replace(pbmc$Matched_Skin_Sample_ID, i, extra)
    log_change("pbmc", pbmc$Sample_ID, "Matched_Skin_Sample_ID", pbmc$Matched_Skin_Sample_ID, new)
    pbmc$Matched_Skin_Sample_ID <- new
    add_count("Sample_ID typo: metadata corrected to the matrix column", 1L)
  } else if (meta_ok && !col_ok) {
    id_fix <- data.table(matrix_column = extra, correct_Sample_ID = skin$Sample_ID[miss],
                         note = "matrix column misspelt; rename at import")
    f <- file.path(CLINICAL, "bulk_skin_matrix_column_fixes.csv")
    fwrite(id_fix, f); Sys.chmod(f, "0600")
    add_count("Sample_ID typo: metadata correct, matrix column misspelt (rename file written)", 1L)
  } else {
    finalize_run(run, status = "failed", verdict = "Sample_ID typo: neither spelling matches Subject_ID + Timepoint; owner review")
    stop("Sample_ID typo unresolved: neither spelling matches the expected ID", call. = FALSE)
  }
} else add_count("Sample_ID typo: not the expected 1 vs 1 pattern; nothing changed", 0L)

# ---- write back only the changed columns ---------------------------
write_cols <- function(f, sheet, x0, x) {
  path <- file.path(CLINICAL, f)
  wb <- loadWorkbook(path)
  sh <- if (is.numeric(sheet)) names(wb)[sheet] else sheet
  changed <- names(x)[vapply(names(x), function(c) !identical(x[[c]], x0[[c]]), TRUE)]
  for (c in changed) writeData(wb, sh, x = x[[c]], startCol = match(c, names(x0)), startRow = 2, colNames = FALSE)
  tmp <- file.path(dirname(path), paste0(".tmp_", run$run_id, "_", f))
  saveWorkbook(wb, tmp, overwrite = TRUE)
  # verification on the temporary file: every cell outside the changed columns is identical
  y <- as.data.table(read_excel(tmp, sheet = sh, col_types = "text"))
  stopifnot(identical(names(y), names(x0)), nrow(y) == nrow(x0))
  for (c in setdiff(names(x0), changed)) if (!identical(y[[c]], x0[[c]])) stop(f, ": column ", c, " changed unexpectedly")
  for (c in changed) if (!identical(y[[c]], x[[c]])) stop(f, ": column ", c, " not written as intended")
  if (length(changed)) stopifnot(file.rename(tmp, path)) else unlink(tmp)
  changed
}
ch_skin <- write_cols(P$skin_file, 1, skin0, skin)
ch_pbmc <- write_cols(P$pbmc_file, 1, pbmc0, pbmc)

counts <- rbindlist(counts)
save_table(counts, run, "changes_by_column")
save_table(transitions, run, "mRSS_category_transitions")           # category labels only
chg <- rbindlist(changes, fill = TRUE)
save_patient_table(chg, run, "change_log")
if (nrow(id_fix)) save_patient_table(id_fix, run, "matrix_column_fixes")

writeLines(c(sprintf("Backup: %s", .rel_to_repo(bk)),
             sprintf("Skin workbook, columns rewritten: %s", paste(ch_skin, collapse = ", ")),
             sprintf("PBMC workbook, columns rewritten: %s", if (length(ch_pbmc)) paste(ch_pbmc, collapse = ", ") else "none"),
             "Rows changed per column:", sprintf("  %-5s %-75s %d", counts$table, counts$column, counts$rows_changed),
             "mRSS_category transitions (skin value -> master value, rows):",
             if (nrow(transitions)) sprintf("  %s -> %s: %d", transitions$skin_value, transitions$master_value, transitions$N) else "  none",
             "Not touched: autoantibody and escape columns (owner)."),
           file.path(run$dir, "CLEAN_SUMMARY.txt"))
finalize_run(run, status = "ok", verdict = sprintf("skin: %s rewritten; pbmc: %s; backup %s",
  paste(ch_skin, collapse = "/"), if (length(ch_pbmc)) paste(ch_pbmc, collapse = "/") else "none", .rel_to_repo(bk)))
