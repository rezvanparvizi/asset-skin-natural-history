# ==============================================================
# 00_8_finalize_bulk_metadata.R   (one-time fix; removed from the repo
# after it ran — recoverable from git at the sha in docs/runs.csv)
#
# 1. PBMC workbook: add Use_downstream / Use_reason (as the skin
#    workbook). Rule: drop QC_flag "exclude_suggested"; one sample per
#    subject — prefer QC_flag "ok" over "check", then the larger
#    raw_library_size. Must give 67 samples (owner).
#    The PBMC workbook is unchanged since R0033, whose backup holds it
#    (md5 checked), so no new backup copy is made.
# 2. Skin matrices: the one misspelt Sample_ID column header (found by
#    00_7) is corrected in all 6 files under data/bulk/skin/, so every
#    script reads the matrices without a rename map. Header line only;
#    each file is rewritten via a temp file and verified (header = the
#    metadata IDs; line count unchanged). The rename map is moved into
#    data/clinical/backup/R0033/ as the record.
# Counts only to tables/.
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(openxlsx)
})

P <- list(pbmc_file = "ASSET_clinical_data_PBMC_bulk_metadata.xlsx",
          skin_file = "ASSET_clinical_data_skin_bulk_metadata.xlsx",
          backup_of_originals = "R0033", expected_pbmc_use = 67L)

run <- init_run(stage = "00_manifest_and_batch", run_name = "finalize_bulk_metadata", params = P,
                notes = "One-time: PBMC Use_downstream; skin matrix header typo corrected")
na_txt <- function(v) { v[!is.na(v) & toupper(trimws(v)) %in% c("NA", "")] <- NA; v }
counts <- list(); add <- function(i, n) counts[[length(counts) + 1]] <<- data.table(item = i, n = n)

# ---- 1. PBMC Use_downstream ----------------------------------------
path <- file.path(CLINICAL, P$pbmc_file)
bk <- file.path(CLINICAL, "backup", P$backup_of_originals, P$pbmc_file)
stopifnot(file.exists(bk), unname(tools::md5sum(path)) == unname(tools::md5sum(bk)))
pbmc <- as.data.table(read_excel(path, col_types = "text"))
q <- pbmc[, .(Sample_ID, Subject_ID, flag = na_txt(QC_flag), depth = as.numeric(raw_library_size))]
q[, rank := frank(list(fifelse(flag == "ok", 0L, 1L), -depth), ties.method = "first"), by = Subject_ID]
q[, Use_reason := fifelse(flag == "exclude_suggested", "QC: exclude_suggested", NA_character_)]
q[is.na(Use_reason), rank2 := frank(rank, ties.method = "first"), by = Subject_ID]
q[is.na(Use_reason) & rank2 > 1, Use_reason := "Technical replicate (other sample of the subject kept)"]
q[, Use_downstream := fifelse(is.na(Use_reason), "Yes", "No")]
stopifnot(q[Use_downstream == "Yes", .N] == P$expected_pbmc_use, q[Use_downstream == "Yes", !anyDuplicated(Subject_ID)])
add("PBMC samples", nrow(q)); add("PBMC Use_downstream = Yes", q[Use_downstream == "Yes", .N])
for (r in na.omit(unique(q$Use_reason))) add(paste("PBMC No:", r), q[Use_reason == r, .N])

wb <- loadWorkbook(path); sh <- names(wb)[1]; hdr <- names(pbmc)
for (c in c("Use_downstream", "Use_reason")) {
  j <- match(c, hdr); if (is.na(j)) { hdr <- c(hdr, c); j <- length(hdr) }
  writeData(wb, sh, x = c, startCol = j, startRow = 1)
  writeData(wb, sh, x = q[[c]], startCol = j, startRow = 2, colNames = FALSE)
}
tmp <- file.path(CLINICAL, paste0(".tmp_", run$run_id, "_", P$pbmc_file))
saveWorkbook(wb, tmp, overwrite = TRUE)
y <- as.data.table(read_excel(tmp, col_types = "text"))
stopifnot(identical(names(y), hdr))
for (c in names(pbmc)) if (!identical(y[[c]], pbmc[[c]])) stop("PBMC column ", c, " changed unexpectedly")
stopifnot(identical(y$Use_downstream, q$Use_downstream), file.rename(tmp, path))

# ---- 2. skin matrix header typo --------------------------------------
fix_file <- file.path(CLINICAL, "bulk_skin_matrix_column_fixes.csv")
fx <- fread(fix_file); stopifnot(nrow(fx) == 1)
skin_ids <- as.data.table(read_excel(file.path(CLINICAL, P$skin_file), col_types = "text"))$Sample_ID
mats <- list.files(file.path(BULK, "skin"), pattern = "\\.txt$", full.names = TRUE)
for (f in mats) {
  first <- readLines(f, n = 1)
  sep <- if (grepl("\t", first)) "\t" else ","
  fields <- strsplit(first, sep, fixed = TRUE)[[1]]
  unq <- function(x) gsub('^"|"$', "", x)
  hit <- which(unq(fields) == fx$matrix_column)
  if (length(hit) == 0 && fx$correct_Sample_ID %in% unq(fields)) { add(paste(basename(f), "already correct"), 1L); next }
  stopifnot(length(hit) == 1)
  fields[hit] <- sub(fx$matrix_column, fx$correct_Sample_ID, fields[hit], fixed = TRUE)
  tmp <- paste0(f, ".tmp_", run$run_id)
  writeLines(paste(fields, collapse = sep), tmp)
  stopifnot(system(sprintf("tail -n +2 %s >> %s", shQuote(f), shQuote(tmp))) == 0)
  n_old <- as.integer(sub(" .*", "", system2("wc", c("-l", shQuote(f)), stdout = TRUE)))
  n_new <- as.integer(sub(" .*", "", system2("wc", c("-l", shQuote(tmp)), stdout = TRUE)))
  cols <- names(fread(tmp, nrows = 0))
  stopifnot(n_old == n_new, all(skin_ids %in% cols), !fx$matrix_column %in% cols)
  stopifnot(file.rename(tmp, f))
  add(paste(basename(f), "header corrected; all 234 metadata IDs present"), 1L)
}
stopifnot(file.rename(fix_file, file.path(CLINICAL, "backup", P$backup_of_originals, basename(fix_file))))

counts <- rbindlist(counts)
save_table(counts, run, "finalize_counts")
finalize_run(run, status = "ok", verdict = sprintf("PBMC Use_downstream Yes = %d; skin matrix header corrected in %d files",
  q[Use_downstream == "Yes", .N], sum(grepl("header corrected", counts$item))))
