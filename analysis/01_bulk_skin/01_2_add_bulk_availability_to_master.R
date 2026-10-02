# ==============================================================
# 01_2_add_bulk_availability_to_master.R
#
# Adds `Bulk_baseline_data` to the clinical master (in place, after a
# backup; decisions.md 2026-10-02): whether each patient has usable
# BASELINE bulk data — "skin+PBMC", "skin only", "PBMC only" or "none".
#   skin : a Baseline skin sample with Use_downstream = "Yes" (01_1)
#   PBMC : a PBMC sample whose QC_flag is not "exclude_suggested" (owner's
#          PBMC QC; all PBMC samples are baseline)
# Counts to tables/ (also: how many patients the QC filters change).
#
#   jobs/run.sh analysis/01_bulk_skin/01_2_add_bulk_availability_to_master.R --freeze freeze01 --cohort reference
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
          column = "Bulk_baseline_data")

run <- init_run(
  stage    = "01_bulk_skin",
  run_name = "add_bulk_availability_to_master",
  params   = P,
  notes    = "Master column: usable baseline bulk skin / PBMC per patient"
)

rd <- function(f, sheet = 1) as.data.table(read_excel(file.path(CLINICAL, f), sheet = sheet, col_types = "text"))
na_txt <- function(v) { v[!is.na(v) & toupper(trimws(v)) %in% c("NA", "")] <- NA; v }
master <- rd(basename(CLINICAL_MASTER)); skin <- rd(P$skin_file); pbmc <- rd(P$pbmc_file)
stopifnot(nrow(master) == 88, "Use_downstream" %in% names(skin), "QC_flag" %in% names(pbmc))

skin_any  <- unique(skin[Timepoint == "Baseline", Subject_ID])
skin_use  <- unique(skin[Timepoint == "Baseline" & Use_downstream == "Yes", Subject_ID])
pbmc_any  <- unique(pbmc$Subject_ID)
pbmc_use  <- unique(pbmc[na_txt(QC_flag) != "exclude_suggested" | is.na(na_txt(QC_flag)), Subject_ID])
s <- master$Subject_ID %in% skin_use; b <- master$Subject_ID %in% pbmc_use
val <- fifelse(s & b, "skin+PBMC", fifelse(s, "skin only", fifelse(b, "PBMC only", "none")))

tab <- as.data.table(table(Bulk_baseline_data = val, Treatment_arm = master$Treatment_arm))
tab <- dcast(tab, Bulk_baseline_data ~ Treatment_arm, value.var = "N")
tab[, Total := rowSums(.SD), .SDcols = setdiff(names(tab), "Bulk_baseline_data")]
pl <- master$Treatment_arm == "Placebo"
by_cat <- as.data.table(table(Bulk_baseline_data = val[pl], mRSS_category = master$mRSS_category[pl]))
save_table(tab, run, "bulk_baseline_by_arm")
save_table(by_cat, run, "bulk_baseline_placebo_by_mRSS_category")
filt <- data.table(item = c("patients with any Baseline skin sample", "  ... usable (Use_downstream = Yes)",
                            "patients with any PBMC sample", "  ... usable (not exclude_suggested)"),
                   n = c(sum(master$Subject_ID %in% skin_any), sum(s), sum(master$Subject_ID %in% pbmc_any), sum(b)))
save_table(filt, run, "effect_of_qc_filters")

# ---- write (backup, temp file, verify) ----------------------------
bk <- file.path(CLINICAL, "backup", run$run_id)
dir.create(bk, recursive = TRUE, showWarnings = FALSE, mode = "0700")
src <- CLINICAL_MASTER
stopifnot(file.copy(src, file.path(bk, basename(src)), copy.date = TRUE)); Sys.chmod(file.path(bk, basename(src)), "0400")
wb <- loadWorkbook(src); sh <- names(wb)[1]
j <- match(P$column, names(master)); if (is.na(j)) j <- ncol(master) + 1
writeData(wb, sh, x = P$column, startCol = j, startRow = 1)
writeData(wb, sh, x = val, startCol = j, startRow = 2, colNames = FALSE)
tmp <- file.path(CLINICAL, paste0(".tmp_", run$run_id, "_", basename(src)))
saveWorkbook(wb, tmp, overwrite = TRUE)
y <- as.data.table(read_excel(tmp, sheet = sh, col_types = "text"))
for (c in setdiff(names(master), P$column)) if (!identical(y[[c]], master[[c]])) stop("column ", c, " changed unexpectedly")
stopifnot(identical(y[[P$column]], val), identical(excel_sheets(tmp), excel_sheets(src)))
stopifnot(file.rename(tmp, src))

writeLines(c("Bulk_baseline_data by arm:", capture.output(print(tab)), "",
             "Placebo by mRSS_category:", capture.output(print(dcast(by_cat, Bulk_baseline_data ~ mRSS_category, value.var = "N"))), "",
             sprintf("%-45s %d", filt$item, filt$n), "", sprintf("Master updated; backup %s", .rel_to_repo(bk))),
           file.path(run$dir, "AVAILABILITY_SUMMARY.txt"))
finalize_run(run, status = "ok", verdict = paste(sprintf("%s %d", tab$Bulk_baseline_data, tab$Total), collapse = "; "))
