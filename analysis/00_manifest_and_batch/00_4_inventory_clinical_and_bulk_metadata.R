# ==============================================================
# 00_4_inventory_clinical_and_bulk_metadata.R
#
# Structure of the clinical / bulk-sample workbooks and the bulk matrices,
# WITHOUT any values (decisions.md 2026-10-02, "scripts-only access"):
#
#   metadata/design/clinical_file_inventory.csv
#       file, sheet, rows, columns — every .xlsx/.csv under data/clinical/
#   metadata/design/clinical_column_inventory.csv
#       file, sheet, position, column name, R class, n missing, n distinct
#   metadata/design/bulk_matrix_inventory.csv
#       every matrix under data/bulk/: rows, sample columns, name of the
#       first (gene) column — header line and row count only
#
# These three are agent-readable. No IDs, no values. Read-only on data/.
#
#   jobs/run.sh analysis/00_manifest_and_batch/00_4_inventory_clinical_and_bulk_metadata.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

run <- init_run(
  stage    = "00_manifest_and_batch",
  run_name = "inventory_clinical_and_bulk_metadata",
  params   = list(clinical_dir = .rel_to_repo(CLINICAL), bulk_dir = .rel_to_repo(BULK)),
  notes    = "File / sheet / column structure of clinical and bulk metadata; no values"
)

write_design <- function(x, name) {
  f <- file.path(DESIGN, paste0(name, ".csv"))
  utils::write.csv(x, f, row.names = FALSE, na = "")
  Sys.chmod(f, "0600")
  cat(.rel_to_repo(f), "\n", file = file.path(run$tables, "DESIGN_TABLES.txt"), append = TRUE)
  message("  wrote ", .rel_to_repo(f))
  invisible(f)
}

# ---- 1. clinical workbooks -----------------------------------

files <- list.files(CLINICAL, pattern = "\\.(xlsx|xls|csv)$", recursive = TRUE, ignore.case = TRUE)
files <- files[!grepl("^backup/", files)]
if (!length(files)) stop("No .xlsx/.csv files under ", CLINICAL, call. = FALSE)

file_inv <- list(); col_inv <- list()
for (f in files) {
  path <- file.path(CLINICAL, f)
  sheets <- if (grepl("\\.csv$", f, ignore.case = TRUE)) NA_character_ else excel_sheets(path)
  for (s in sheets) {
    x <- if (is.na(s)) fread(path, colClasses = "character", na.strings = c("", "NA"))
         else as.data.table(read_excel(path, sheet = s, guess_max = 10000, .name_repair = "minimal"))
    file_inv[[length(file_inv) + 1]] <- data.table(file = f, sheet = s, rows = nrow(x), columns = ncol(x),
                                                   duplicated_column_names = sum(duplicated(names(x))))
    if (ncol(x)) col_inv[[length(col_inv) + 1]] <- data.table(
      file = f, sheet = s, position = seq_len(ncol(x)), column = names(x),
      class = vapply(x, function(v) class(v)[1], ""),
      n_missing = vapply(x, function(v) sum(is.na(v) | (is.character(v) & !nzchar(trimws(v)))), 0L),
      n_distinct = vapply(x, function(v) uniqueN(v[!is.na(v)]), 0L))
  }
}
file_inv <- rbindlist(file_inv); col_inv <- rbindlist(col_inv)
write_design(file_inv, "clinical_file_inventory")
write_design(col_inv, "clinical_column_inventory")

# ---- 2. bulk matrices: header and row count only -------------

mats <- list.files(BULK, pattern = "\\.(txt|tsv|csv)$", recursive = TRUE)
mats <- mats[!grepl("sample_info|batch_info|README", mats, ignore.case = TRUE)]
mat_inv <- rbindlist(lapply(mats, function(m) {
  path <- file.path(BULK, m)
  hdr <- names(fread(path, nrows = 0))
  rows <- as.integer(system2("wc", c("-l", shQuote(path)), stdout = TRUE) |> sub(pattern = "^\\s*(\\d+).*", replacement = "\\1"))
  data.table(file = m, rows_excluding_header = rows - 1L, columns = length(hdr),
             first_column_name = hdr[1], header_blank_first = !nzchar(hdr[1]))
}))
write_design(mat_inv, "bulk_matrix_inventory")

summary_lines <- c(
  sprintf("Clinical files: %d (%s)", length(files), paste(files, collapse = "; ")),
  sprintf("  %s [%s]: %d rows x %d columns", file_inv$file, file_inv$sheet, file_inv$rows, file_inv$columns),
  sprintf("Bulk matrices: %d", nrow(mat_inv)),
  sprintf("  %s: %d rows x %d columns (first column '%s')", mat_inv$file, mat_inv$rows_excluding_header,
          mat_inv$columns, mat_inv$first_column_name))
writeLines(summary_lines, file.path(run$dir, "INVENTORY_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

finalize_run(run, status = "ok", verdict = sprintf("%d clinical files, %d sheets, %d columns inventoried; %d bulk matrices",
                                                   length(files), nrow(file_inv), nrow(col_inv), nrow(mat_inv)))
