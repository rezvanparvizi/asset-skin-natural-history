# ==============================================================
# 03_2_qc_loss_by_inherited_celltype.R
#
# Does the 03_1 cell QC (R0017) remove some cell types more than others?
# Low-UMI types (T cells, mast, smooth muscle run at ~0.8x the median UMI)
# could be depleted by the 500-UMI / 200-gene floor or the MT cut, which
# would bias every rare-population analysis downstream.
#
# Labels: the colleagues' labels joined to our cells by 06_0 (R0021),
#   wasikowr_basectrl__celltype  16 broad types; Baseline + Control cells of
#                                jarnagin's base_control object (about half
#                                of those cells: her subset bug drops every
#                                other cell, independently of type, so the
#                                sample is not biased by type)
#   jf_immune__celltype          immune types (B, T, myeloid, Langerhans),
#                                all timepoints
# Their cells passed wasikowr's >200-gene filter, so cells below that are
# absent from both sources: loss rates are for cells she kept.
#
# Output (aggregate only): per label, cells, fraction failing our QC by
# reason, fraction under the UMI floor alone, UMI distribution.
# Nothing here looks at mRSS_category.
#
#   jobs/run.sh analysis/03_sc_qc/03_2_qc_loss_by_inherited_celltype.R \
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

P <- list(label_columns = c("wasikowr_basectrl__celltype", "jf_immune__celltype"),
          min_umi = 500, min_genes = 200, excess_flag = 1.5)

run <- init_run(
  stage    = "03_sc_qc",
  run_name = "qc_loss_by_inherited_celltype",
  params   = P,
  notes    = "QC failure rate per inherited cell type (does the floor / MT cut deplete a lineage?)"
)

qc <- fread(file.path(results_dir("03_sc_qc", "cell_qc_per_library", cohort = "reference"),
                      "objects", "cell_qc.csv.gz"),
            select = c("cell_id", "library_id", "n_umi", "n_genes", "pct_mt", "qc_pass", "qc_fail_reason"))
lab_file <- file.path(METADATA, "inherited_labels", sprintf("%s_upstream_labels.csv.gz", FREEZE))
if (!file.exists(lab_file)) stop("06_0 label table not found: ", lab_file, call. = FALSE)
lab <- fread(lab_file, select = c("cell_id", P$label_columns))
stopifnot(!anyDuplicated(qc$cell_id), !anyDuplicated(lab$cell_id),
          all(lab$cell_id %in% qc$cell_id))

# restrict to the 218-library freeze, as R0017 did for the pass/fail calls
qc <- qc[library_id %in% freeze_libraries()$library_id]
d <- merge(qc, lab, by = "cell_id", all.x = TRUE)
d[, under_floor := n_umi < P$min_umi | n_genes < P$min_genes]
d[is.na(qc_fail_reason), qc_fail_reason := ""]

reasons <- c("floor", "low_genes", "low_umi", "high_umi", "high_mt")
per_label <- rbindlist(lapply(P$label_columns, function(col) {
  x <- d[!is.na(get(col)) & nzchar(get(col))]
  overall_fail <- mean(!x$qc_pass)
  out <- x[, c(list(cells = .N,
                    median_umi = as.numeric(median(n_umi)),
                    fail = mean(!qc_pass),
                    under_floor = mean(under_floor)),
               setNames(lapply(reasons, function(r) mean(grepl(r, qc_fail_reason, fixed = TRUE))),
                        paste0("fail_", reasons))),
           by = .(label = get(col))]
  out[, `:=`(source = col, fail_vs_source_overall = fail / overall_fail)]
  setcolorder(out, c("source", "label"))
  out[order(-fail_vs_source_overall)]
}))
per_label[, flag := ifelse(fail_vs_source_overall >= P$excess_flag & cells >= 100,
                           sprintf("fails >= %.1fx the source average", P$excess_flag), "")]
save_table(per_label, run, "qc_loss_by_label")

# UMI distribution by broad type, with the floor
col1 <- P$label_columns[1]
pd <- d[!is.na(get(col1)) & nzchar(get(col1)), .(label = get(col1), n_umi, qc_pass)]
ord <- pd[, .(m = median(n_umi)), by = label][order(m)]$label
pd[, label := factor(label, ord)]
p <- ggplot(pd, aes(n_umi, label)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_vline(xintercept = P$min_umi, linetype = 2, linewidth = 0.3) +
  scale_x_log10() +
  labs(x = "UMI per cell (raw, log scale)", y = NULL,
       subtitle = "wasikowr's broad labels (Baseline + Control); dashed = 500-UMI floor") +
  theme_asset()
save_plot(p, run, "umi_by_inherited_celltype", width = 6, height = 5)

flagged <- per_label[nzchar(flag)]
fmt <- function(x) sprintf("%5.1f%%", 100 * x)
summary_lines <- c(
  "QC failure by inherited label (cells she kept; our 03_1 calls):",
  unlist(lapply(P$label_columns, function(col) {
    x <- per_label[source == col]
    c(sprintf("  %s  (overall fail %s)", col, fmt(sum(x$fail * x$cells) / sum(x$cells))),
      sprintf("    %-20s %7s cells  fail %s  (x%.2f)  floor %s  high_mt %s  median UMI %6.0f",
              x$label, format(x$cells, big.mark = ","), fmt(x$fail), x$fail_vs_source_overall,
              fmt(x$fail_floor), fmt(x$fail_high_mt), x$median_umi))
  })),
  sprintf("Labels failing >= %.1fx their source average (>= 100 cells): %d", P$excess_flag, nrow(flagged)),
  if (nrow(flagged)) paste0("    ", flagged$source, ": ", flagged$label))
writeLines(summary_lines, file.path(run$dir, "QC_LOSS_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

finalize_run(run, status = "ok", verdict = sprintf(
  "%d labels checked; %d fail >= %.1fx their source average (%s)", nrow(per_label), nrow(flagged),
  P$excess_flag, if (nrow(flagged)) paste(flagged$label, collapse = ", ") else "none"))
