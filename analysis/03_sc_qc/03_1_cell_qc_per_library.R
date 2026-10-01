# ==============================================================
# 03_1_cell_qc_per_library.R
#
# Per-cell QC metrics and pass/fail calls for every library in the
# freeze (Layer 1: all arms, all timepoints, HC).
#
# Rule (docs/decisions.md, 2026-10-01, fixed before looking at data):
# per library, a cell fails if ANY of
#   log10 UMI   < median - 3 MAD        log10 genes < median - 3 MAD
#   genes < 200 or UMI < 500            (hard floor)
#   log10 UMI   > median + 5 MAD        (crude doublet guard)
#   % MT > median + 3 MAD AND > 10%     (only if the panel has MT- genes)
# Failing cells are kept with qc_pass = FALSE and a reason. Libraries
# with < 500 passing cells or > 50% failing are FLAGGED for review, not
# removed.
#
# Input : data/raw/<run_id>/<library_id>/sample_filtered_feature_bc_matrix.h5
# Output: objects/cell_qc.csv.gz            one row per cell (no subject IDs)
#         tables/library_qc_summary.csv     one row per library
#         tables/qc_pass_by_design.csv      pass fraction by pool/arm/timepoint/group
#         plots/                            per-library distributions with thresholds
#         QC_SUMMARY.txt                    what to review
#
# Nothing in this script looks at mRSS_category.
#
#   jobs/run.sh analysis/03_sc_qc/03_1_cell_qc_per_library.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")
source("R/theme_asset.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(BPCells)
})

QC <- list(nmads_low = 3, nmads_high_umi = 5, min_genes = 200, min_umi = 500,
           mt_nmads = 3, mt_floor_pct = 10,
           flag_min_pass_cells = 500, flag_max_fail_frac = 0.5)

run <- init_run(
  stage    = "03_sc_qc",
  run_name = "cell_qc_per_library",
  params   = QC,
  notes    = "Per-library adaptive (MAD) cell QC with hard floor; flags, does not drop"
)

THREADS <- as.integer(Sys.getenv("OMP_NUM_THREADS", "4"))

man <- freeze_libraries()
message(sprintf("Freeze %s: %d libraries", FREEZE, nrow(man)))

# ---- 1. per-library metrics -----------------------------------

read_gene_names <- function(f) {
  nm <- as.character(rhdf5::h5read(f, "matrix/features/name"))
  ft <- as.character(rhdf5::h5read(f, "matrix/features/feature_type"))
  rhdf5::h5closeAll()
  nm[ft == "Gene Expression"]
}

mad_bounds <- function(x, nmads) {
  lx <- log10(x)
  m  <- stats::median(lx)
  d  <- stats::mad(lx)
  c(low = 10^(m - nmads * d), high = 10^(m + nmads * d))
}

qc_library <- function(i) {
  lib <- man$library_id[i]
  f <- file.path(RAW, man$raw_dir[i], "sample_filtered_feature_bc_matrix.h5")
  if (!file.exists(f)) stop("missing: ", f)

  mat   <- BPCells::open_matrix_10x_hdf5(f, feature_type = "Gene Expression")
  genes <- read_gene_names(f)
  if (length(genes) != nrow(mat)) {
    stop(lib, ": ", length(genes), " gene names vs ", nrow(mat), " matrix rows")
  }

  n_umi   <- colSums(mat)
  n_genes <- BPCells::matrix_stats(mat, col_stats = "nonzero")$col_stats["nonzero", ]
  pct_of  <- function(idx) {
    if (!length(idx)) return(rep(NA_real_, ncol(mat)))
    100 * colSums(mat[idx, ]) / pmax(n_umi, 1)
  }
  mt   <- which(grepl("^MT-", genes))
  dt <- data.table(
    library_id = lib,
    barcode    = colnames(mat),
    n_umi      = as.numeric(n_umi),
    n_genes    = as.numeric(n_genes),
    pct_mt     = pct_of(mt),
    pct_ribo   = pct_of(which(grepl("^RP[SL][0-9]", genes))),
    pct_hb     = pct_of(which(genes %in% c("HBB", "HBA1", "HBA2"))),
    pct_krt    = pct_of(which(grepl("^KRT[0-9]", genes))),
    pct_col    = pct_of(which(genes %in% c("COL1A1", "COL1A2", "COL3A1")))
  )

  # thresholds, this library only
  b_umi   <- mad_bounds(dt$n_umi,   QC$nmads_low)
  b_genes <- mad_bounds(dt$n_genes, QC$nmads_low)
  hi_umi  <- mad_bounds(dt$n_umi,   QC$nmads_high_umi)["high"]
  has_mt  <- length(mt) > 0
  mt_hi   <- if (has_mt) {
    max(stats::median(dt$pct_mt) + QC$mt_nmads * stats::mad(dt$pct_mt), QC$mt_floor_pct)
  } else NA_real_

  r <- data.table(
    low_umi   = dt$n_umi   < b_umi["low"],
    low_genes = dt$n_genes < b_genes["low"],
    floor     = dt$n_genes < QC$min_genes | dt$n_umi < QC$min_umi,
    high_umi  = dt$n_umi   > hi_umi,
    high_mt   = if (has_mt) dt$pct_mt > mt_hi else rep(FALSE, nrow(dt))
  )
  dt[, qc_fail_reason := apply(r, 1, function(z) paste(names(r)[z], collapse = ";"))]
  dt[, qc_pass := !nzchar(qc_fail_reason)]

  thr <- data.table(library_id = lib, n_mt_genes = length(mt),
                    thr_umi_low = unname(b_umi["low"]),
                    thr_genes_low = unname(b_genes["low"]),
                    thr_umi_high = unname(hi_umi), thr_mt_high = mt_hi)
  list(cells = dt, thr = thr)
}

res <- parallel::mclapply(seq_len(nrow(man)), function(i) {
  tryCatch(qc_library(i), error = function(e) e)
}, mc.cores = THREADS, mc.preschedule = FALSE)

err <- vapply(res, inherits, logical(1), "error")
if (any(err)) {
  msg <- paste(man$library_id[err], vapply(res[err], conditionMessage, ""),
               sep = ": ", collapse = "\n")
  finalize_run(run, status = "failed",
               verdict = sprintf("%d libraries failed to read", sum(err)))
  stop("Libraries failed:\n", msg, call. = FALSE)
}

cells <- rbindlist(lapply(res, `[[`, "cells"))
thr   <- rbindlist(lapply(res, `[[`, "thr"))
cells[, cell_id := paste(library_id, barcode, sep = "_")]
stopifnot(!anyDuplicated(cells$cell_id),
          setequal(unique(cells$library_id), man$library_id))

# one row per cell; library-level IDs only, no subject IDs
fwrite(cells[, .(cell_id, library_id, barcode, n_umi, n_genes, pct_mt,
                 pct_ribo, pct_hb, pct_krt, pct_col, qc_pass, qc_fail_reason)],
       file.path(run$objects, "cell_qc.csv.gz"))
message("  wrote ", .rel_to_repo(file.path(run$objects, "cell_qc.csv.gz")))

# ---- 2. per-library summary ----------------------------------

design_cols <- c("library_id", "batch_id", "group", "arm", "timepoint",
                 "cellranger_run", "design_flag")
reasons <- c("low_umi", "low_genes", "floor", "high_umi", "high_mt")
lib_sum <- cells[, c(list(cells_called = .N, cells_pass = sum(qc_pass),
                          median_umi = stats::median(n_umi),
                          median_genes = stats::median(n_genes),
                          median_pct_mt = stats::median(pct_mt),
                          median_pct_hb = stats::median(pct_hb),
                          median_pct_krt = stats::median(pct_krt),
                          median_pct_col = stats::median(pct_col)),
                     lapply(setNames(reasons, paste0("n_", reasons)),
                            function(r) sum(grepl(r, qc_fail_reason, fixed = TRUE)))),
                 by = library_id]
lib_sum[, frac_pass := cells_pass / cells_called]
lib_sum <- merge(lib_sum, thr, by = "library_id")
lib_sum <- merge(as.data.table(man[, design_cols]), lib_sum, by = "library_id")
lib_sum[, review_flag := trimws(paste(
  ifelse(cells_pass < QC$flag_min_pass_cells, "few_passing_cells", ""),
  ifelse(1 - frac_pass > QC$flag_max_fail_frac, "most_cells_fail", ""),
  ifelse(!is.na(design_flag) & nzchar(design_flag), design_flag, "")))]
setorder(lib_sum, batch_id, library_id)
save_table(lib_sum, run, "library_qc_summary")

# ---- 3. does QC depend on design? ----------------------------
# QC removal should not track timepoint or arm; if it does, QC itself
# becomes a confounder. Tested within SSc, libraries as units.

ssc <- lib_sum[group == "SSc"]
by_design <- rbindlist(lapply(c("batch_id", "arm", "timepoint", "group"), function(v) {
  d <- if (v == "group") lib_sum else ssc
  s <- d[, .(libraries = .N, median_frac_pass = stats::median(frac_pass),
             min_frac_pass = min(frac_pass), cells_pass = sum(cells_pass)),
         by = c(v)]
  setnames(s, v, "level")
  p <- stats::kruskal.test(d$frac_pass, factor(d[[v]]))$p.value
  cbind(variable = v, s, kruskal_p_frac_pass = signif(p, 3))
}))
save_table(by_design, run, "qc_pass_by_design")

# ---- 4. plots ------------------------------------------------

set.seed(1)
plot_cells <- cells[, .SD[sample(.N, min(.N, 3000))], by = library_id]
plot_cells <- merge(plot_cells, lib_sum[, .(library_id, batch_id, group, timepoint)],
                    by = "library_id")
lib_order <- lib_sum$library_id
plot_cells[, library_id := factor(library_id, levels = lib_order)]
thr_long <- merge(thr, lib_sum[, .(library_id, batch_id)], by = "library_id")
thr_long[, library_id := factor(library_id, levels = lib_order)]

per_lib_violin <- function(y, ylab, low = NULL, high = NULL, log = TRUE) {
  p <- ggplot(plot_cells, aes(library_id, .data[[y]], fill = group)) +
    geom_violin(scale = "width", linewidth = 0.1) +
    facet_grid(~ batch_id, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = PAL_GROUP) +
    labs(x = NULL, y = ylab, fill = NULL,
         subtitle = "Up to 3,000 cells per library; red ticks = that library's thresholds") +
    theme_asset() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
          strip.text = element_text(size = 5))
  if (log) p <- p + scale_y_log10()
  for (t in c(low, high)) {
    p <- p + geom_point(data = thr_long, aes(library_id, .data[[t]]),
                        inherit.aes = FALSE, shape = 95, size = 3, colour = "red")
  }
  p
}

W <- 24
save_plot(per_lib_violin("n_umi", "UMI per cell", "thr_umi_low", "thr_umi_high"),
          run, "qc_umi_per_library", width = W, height = 5)
save_plot(per_lib_violin("n_genes", "Genes per cell", "thr_genes_low"),
          run, "qc_genes_per_library", width = W, height = 5)
if (any(thr$n_mt_genes > 0)) {
  save_plot(per_lib_violin("pct_mt", "% mitochondrial", high = "thr_mt_high", log = FALSE),
            run, "qc_pct_mt_per_library", width = W, height = 5)
}
for (v in c("pct_hb", "pct_krt", "pct_col")) {
  save_plot(per_lib_violin(v, paste0(v, " (informational, not filtered)"), log = FALSE),
            run, paste0("qc_", v, "_per_library"), qc = TRUE, width = W, height = 5)
}

p_frac <- ggplot(lib_sum, aes(factor(library_id, levels = lib_order), frac_pass,
                              fill = group)) +
  geom_col() +
  geom_hline(yintercept = 1 - QC$flag_max_fail_frac, linetype = 2, linewidth = 0.2) +
  facet_grid(~ batch_id, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = PAL_GROUP) +
  labs(x = NULL, y = "Fraction of called cells passing QC", fill = NULL) +
  theme_asset() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
        strip.text = element_text(size = 5))
save_plot(p_frac, run, "qc_fraction_pass_per_library", width = W, height = 4)

p_scatter <- ggplot(plot_cells, aes(n_umi, n_genes)) +
  geom_bin2d(bins = 120) +
  scale_x_log10() + scale_y_log10() +
  scale_fill_viridis_c(trans = "log10") +
  facet_wrap(~ qc_pass, labeller = label_both) +
  labs(x = "UMI per cell", y = "Genes per cell") +
  theme_asset()
save_plot(p_scatter, run, "qc_umi_vs_genes", width = 7, height = 3.5)

p_design <- ggplot(ssc, aes(timepoint, frac_pass, colour = arm)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 0.6) +
  scale_colour_manual(values = PAL_ARM) +
  labs(x = NULL, y = "Fraction passing QC (per library)", colour = NULL,
       subtitle = "SSc libraries; QC loss should not depend on timepoint or arm") +
  theme_asset()
save_plot(p_design, run, "qc_fraction_pass_by_timepoint_arm", width = 4, height = 3)

p_cells <- ggplot(lib_sum, aes(cells_called, cells_pass, colour = group)) +
  geom_abline(linetype = 2, linewidth = 0.2) +
  geom_point(size = 0.8) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = PAL_GROUP) +
  labs(x = "Cells called by CellRanger", y = "Cells passing QC", colour = NULL) +
  theme_asset()
save_plot(p_cells, run, "qc_cells_called_vs_pass", width = 4, height = 3.5)

# ---- 5. summary for review ------------------------------------

flagged <- lib_sum[nzchar(review_flag)]
reason_tot <- vapply(reasons, function(r) sum(grepl(r, cells$qc_fail_reason, fixed = TRUE)), 0)
p_tp  <- by_design[variable == "timepoint", kruskal_p_frac_pass][1]
p_arm <- by_design[variable == "arm", kruskal_p_frac_pass][1]
p_bat <- by_design[variable == "batch_id", kruskal_p_frac_pass][1]

summary_lines <- c(
  sprintf("Libraries: %d   cells called: %s   pass: %s (%.1f%%)",
          nrow(lib_sum), format(nrow(cells), big.mark = ","),
          format(sum(cells$qc_pass), big.mark = ","), 100 * mean(cells$qc_pass)),
  sprintf("MT- genes in probe panel: %s (n = %s)",
          if (any(thr$n_mt_genes > 0)) "yes" else "NO - % mito not used",
          paste(unique(thr$n_mt_genes), collapse = "/")),
  "Cells failing each criterion (a cell can fail several):",
  paste0("  ", names(reason_tot), ": ", format(reason_tot, big.mark = ",")),
  sprintf("Fraction passing vs timepoint (SSc): Kruskal p = %s", p_tp),
  sprintf("Fraction passing vs arm (SSc):       Kruskal p = %s", p_arm),
  sprintf("Fraction passing vs pool (SSc):      Kruskal p = %s", p_bat),
  sprintf("Libraries flagged for review: %d", nrow(flagged)),
  if (nrow(flagged)) paste0("  ", flagged$library_id, "  [", flagged$batch_id, "]  ",
                            flagged$cells_pass, "/", flagged$cells_called, " pass  ",
                            flagged$review_flag),
  "",
  "Review before integration: plots/qc_*_per_library.pdf (thresholds in red),",
  "plots/qc_fraction_pass_by_timepoint_arm.pdf, tables/library_qc_summary.csv.")
writeLines(summary_lines, file.path(run$dir, "QC_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

finalize_run(run, status = "ok", verdict = sprintf(
  "%s of %s cells pass (%.1f%%); %d libraries flagged for review; pass fraction vs timepoint p=%s, vs arm p=%s",
  format(sum(cells$qc_pass), big.mark = ","), format(nrow(cells), big.mark = ","),
  100 * mean(cells$qc_pass), nrow(flagged), p_tp, p_arm))
