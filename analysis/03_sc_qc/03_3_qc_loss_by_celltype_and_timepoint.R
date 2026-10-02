# ==============================================================
# 03_3_qc_loss_by_celltype_and_timepoint.R
#
# Follow-up to 03_2 (R0027): low-UMI types (T, B, mast) fail the 03_1
# floor more often. What matters for inference is whether that loss
# DIFFERS BY TIMEPOINT — a type lost more at M06 than at M00 would bias
# every longitudinal comparison of it.
#
# Labels with all timepoints (06_0 table, R0021):
#   jf_immune__celltype     B, T, myeloid, Langerhans (jarnagin, full dataset)
#   jf_fibroblast__celltype fibroblasts (jarnagin Round 1)
# SSc libraries only, descriptive (no test): per type x timepoint, cells,
# pooled fail fraction, median of per-library fail fractions (libraries
# with >= 20 cells of the type), and the fail ratio to M00.
# Their cells passed wasikowr's >200-gene filter (03_2 caveat applies).
# Nothing here looks at mRSS_category.
#
#   jobs/run.sh analysis/03_sc_qc/03_3_qc_loss_by_celltype_and_timepoint.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")
source("R/theme_asset.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

P <- list(label_columns = c("jf_immune__celltype", "jf_fibroblast__celltype"),
          min_cells_per_library = 20, ratio_flag = 1.5)

run <- init_run(
  stage    = "03_sc_qc",
  run_name = "qc_loss_by_celltype_and_timepoint",
  params   = P,
  notes    = "Does QC loss of a cell type differ by timepoint? (follow-up to 03_2)"
)

qc <- fread(file.path(results_dir("03_sc_qc", "cell_qc_per_library", cohort = "reference"),
                      "objects", "cell_qc.csv.gz"),
            select = c("cell_id", "library_id", "qc_pass", "qc_fail_reason"))
lab <- fread(file.path(METADATA, "inherited_labels", sprintf("%s_upstream_labels.csv.gz", FREEZE)),
             select = c("cell_id", P$label_columns))
man <- as.data.table(freeze_libraries())[group == "SSc" & timepoint %in% c("M00", "M03", "M06"),
                                         .(library_id, timepoint)]
d <- merge(merge(qc, man, by = "library_id"), lab, by = "cell_id")
d[is.na(qc_fail_reason), qc_fail_reason := ""]

tab <- rbindlist(lapply(P$label_columns, function(col) {
  x <- d[!is.na(get(col)) & nzchar(get(col)), .(label = get(col), library_id, timepoint, qc_pass, qc_fail_reason)]
  per_lib <- x[, .(n = .N, fail = mean(!qc_pass)), by = .(label, timepoint, library_id)][n >= P$min_cells_per_library]
  pooled <- x[, .(cells = .N, libraries = uniqueN(library_id), fail_pooled = mean(!qc_pass),
                  fail_floor = mean(grepl("floor", qc_fail_reason, fixed = TRUE))), by = .(label, timepoint)]
  med <- per_lib[, .(libraries_ge_min = .N, fail_median_per_library = median(fail)), by = .(label, timepoint)]
  out <- merge(pooled, med, by = c("label", "timepoint"), all.x = TRUE)
  out[, source := col]
  out[, ratio_to_M00 := fail_pooled / fail_pooled[timepoint == "M00"], by = label]
  out[]
}))
setcolorder(tab, c("source", "label", "timepoint"))
setorder(tab, source, label, timepoint)
tab[, flag := ifelse(timepoint != "M00" & (ratio_to_M00 >= P$ratio_flag | ratio_to_M00 <= 1 / P$ratio_flag),
                     sprintf("fail rate differs >= %.1fx from M00", P$ratio_flag), "")]
save_table(tab, run, "qc_loss_by_celltype_timepoint")

p <- ggplot(tab, aes(timepoint, fail_pooled, group = label, colour = label)) +
  geom_line() + geom_point() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Cells failing 03_1 QC", colour = NULL,
       subtitle = "SSc cells, colleagues' labels (all timepoints)") +
  theme_asset()
save_plot(p, run, "qc_loss_by_celltype_timepoint", width = 5, height = 3.5)

flagged <- tab[nzchar(flag)]
fmt <- function(x) sprintf("%5.1f%%", 100 * x)
summary_lines <- c(
  "QC failure by cell type and timepoint (SSc; colleagues' labels):",
  sprintf("  %-14s %s  %7s cells  %2d libs  fail %s (floor %s)  median/library %s  x%.2f vs M00",
          tab$label, tab$timepoint, format(tab$cells, big.mark = ","), tab$libraries,
          fmt(tab$fail_pooled), fmt(tab$fail_floor), fmt(tab$fail_median_per_library), tab$ratio_to_M00),
  sprintf("Type x timepoint differing >= %.1fx from M00: %d", P$ratio_flag, nrow(flagged)),
  if (nrow(flagged)) paste0("    ", flagged$label, " ", flagged$timepoint))
writeLines(summary_lines, file.path(run$dir, "QC_LOSS_TIMEPOINT_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

finalize_run(run, status = "ok", verdict = sprintf(
  "%d type x timepoint cells; %d differ >= %.1fx from M00 (%s)", nrow(tab), nrow(flagged), P$ratio_flag,
  if (nrow(flagged)) paste(flagged$label, flagged$timepoint, collapse = ", ") else "none"))
