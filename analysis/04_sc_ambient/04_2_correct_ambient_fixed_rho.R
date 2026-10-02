# ==============================================================
# 04_2_correct_ambient_fixed_rho.R
#
# Ambient RNA correction with ONE fixed contamination fraction for every
# library (owner decision, docs/decisions.md 2026-10-02, "SoupX rho fixed").
#
# Step 1 of stage 04 (04_1, run R0020) already did per library: the soup
# profile from the library's own empty droplets, the fixed clustering SoupX
# needs, and the automatic rho (kept for comparison only). This script is
# step 2: it reuses those clusters and applies rho = SX$rho to all
# libraries. The soup profile is recomputed from the same raw h5 (cheap and
# deterministic). Nothing is re-clustered.
#
#   soup     : that library's own empty droplets (sample_raw h5)
#   clusters : 04_1's per-library clusters (objects/soupx_cells.csv.gz)
#   rho      : SX$rho, identical for every library
#   output   : adjustCounts(roundToInt = TRUE) -> integer corrected counts
#              data/bpcells/<freeze>/counts_corrected/<lib>/
#              (04_1's counts_soupx/ and counts_raw/ are left untouched)
#
# Same checks as 04_1, reported side by side: (b) marker leakage before ->
# after, (c) retention of each lineage's own marker, and the comparison
# with 04_1's automatic-rho correction. Nothing here looks at
# mRSS_category.
#
#   OMP_NUM_THREADS=8 jobs/run.sh analysis/04_sc_ambient/04_2_correct_ambient_fixed_rho.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")
source("R/ambient.R")
source("R/theme_asset.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(BPCells)
  library(SoupX)
})

SX <- list(rho = 0.13, round_to_int = TRUE, retention_flag = 0.70,
           step1_run = "soupx_per_library", out_subdir = "counts_corrected", seed = 1)

run <- init_run(
  stage    = "04_sc_ambient",
  run_name = "correct_ambient_fixed_rho",
  params   = c(SX, list(leakage_genes = AMBIENT_LEAK, own_markers = as.list(AMBIENT_OWN))),
  notes    = "SoupX with one fixed rho for all libraries; clusters and soup from 04_1 (R0020)"
)
THREADS <- as.integer(Sys.getenv("OMP_NUM_THREADS", "4"))
man <- freeze_libraries()

step1 <- results_dir("04_sc_ambient", SX$step1_run, cohort = "reference")
cells <- fread(file.path(step1, "objects", "soupx_cells.csv.gz"))
auto  <- fread(file.path(step1, "tables", "soupx_per_library.csv"))
auto_leak <- fread(file.path(step1, "tables", "leakage_per_library.csv"))
stopifnot(setequal(unique(cells$library_id), man$library_id))

OUT <- file.path(BPCELLS, FREEZE, SX$out_subdir)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- correct one library -------------------------------------

correct_library <- function(i) {
  lib <- man$library_id[i]
  set.seed(SX$seed + i)                        # stochastic rounding, reproducible per library
  d   <- file.path(RAW, man$raw_dir[i])
  toc <- read_h5_gex(file.path(d, "sample_filtered_feature_bc_matrix.h5"))
  tod <- read_h5_gex(file.path(d, "sample_raw_feature_bc_matrix.h5"))
  cl  <- cells[library_id == lib]
  stopifnot(all(cl$barcode %in% colnames(toc)))
  toc <- toc[, cl$barcode, drop = FALSE]
  al  <- align_cell_droplet_genes(toc, tod, lib); toc <- al$toc; tod <- al$tod

  sc  <- SoupChannel(tod, toc, calcSoupProfile = TRUE)
  sc  <- setClusters(sc, setNames(as.character(cl$soupx_cluster), cl$barcode))
  sc  <- setContaminationFraction(sc, SX$rho, forceAccept = TRUE)
  out <- adjustCounts(sc, roundToInt = SX$round_to_int, verbose = 0)
  stopifnot(identical(dim(out), dim(toc)), all(out@x == round(out@x)),
            sum(out) <= sum(toc))

  colnames(out) <- paste(lib, colnames(out), sep = "_")
  dir <- file.path(OUT, lib)
  if (dir.exists(dir)) unlink(dir, recursive = TRUE)
  BPCells::write_matrix_dir(BPCells::convert_matrix_type(as(out, "IterableMatrix"), "uint32_t"),
                            dir, overwrite = TRUE)
  colnames(out) <- cl$barcode
  chk <- ambient_checks(toc, out, cl, lib)
  chk[, `:=`(total_raw = sum(toc), total_corrected = sum(out))]
  message(sprintf("  done: %s (%d cells, %.1f%% removed)", lib, ncol(toc), 100 * (1 - sum(out) / sum(toc))))
  chk
}

message(sprintf("Correcting %d libraries at rho = %.3f, %d threads", nrow(man), SX$rho, THREADS))
res <- parallel::mclapply(seq_len(nrow(man)), function(i)
  tryCatch(correct_library(i), error = function(e) e),
  mc.cores = THREADS, mc.preschedule = FALSE)
err <- vapply(res, inherits, logical(1), "error")
if (any(err)) {
  msg <- paste(man$library_id[err], vapply(res[err], conditionMessage, ""), sep = ": ", collapse = "\n")
  finalize_run(run, status = "failed", verdict = sprintf("correction failed for %d libraries", sum(err)))
  stop("Correction failed:\n", msg, call. = FALSE)
}
leak <- rbindlist(res)

# ---- tables --------------------------------------------------

tot <- unique(leak[, .(library_id, total_raw, total_corrected)])
lib_tab <- merge(tot, auto[, .(library_id, batch_id, group, arm, timepoint, n_cells,
                               rho_auto = rho_used, frac_removed_auto = frac_counts_removed)],
                 by = "library_id")
lib_tab[, `:=`(rho = SX$rho, frac_removed = 1 - total_corrected / total_raw)]
ret <- dcast(leak[grepl("^retention_", metric)], library_id ~ metric, value.var = "ratio")
lib_tab <- merge(lib_tab, ret, by = "library_id")
ret_cols <- grep("^retention_", names(lib_tab), value = TRUE)
lib_tab[, review_flag := ifelse(apply(as.matrix(.SD) < SX$retention_flag, 1, any, na.rm = TRUE),
                                "own_marker_loss", ""), .SDcols = ret_cols]
setorder(lib_tab, batch_id, library_id)
stopifnot(nrow(lib_tab) == nrow(man))
save_table(lib_tab, run, "ambient_per_library")
save_table(leak[, !c("total_raw", "total_corrected")], run, "leakage_per_library")

# leakage: raw vs automatic rho (04_1) vs fixed rho (here)
lk <- leak[!grepl("^retention_", metric), .(library_id, metric, raw, fixed = corrected)]
lk <- merge(lk, auto_leak[, .(library_id, metric, auto = soupx)], by = c("library_id", "metric"))
med <- lk[, .(raw = median(raw, na.rm = TRUE), auto_rho = median(auto, na.rm = TRUE),
              fixed_rho = median(fixed, na.rm = TRUE)), by = metric]
save_table(med, run, "leakage_median_raw_auto_fixed")

# ---- plots ---------------------------------------------------

lk_long <- melt(lk, id.vars = c("library_id", "metric"), measure.vars = c("raw", "auto", "fixed"),
                variable.name = "counts", value.name = "fraction")
lk_long[, counts := factor(counts, c("raw", "auto", "fixed"),
                           c("raw", "auto rho (04_1)", sprintf("rho %.2f", SX$rho)))]
p_leak <- ggplot(lk_long, aes(counts, fraction)) +
  geom_line(aes(group = library_id), alpha = 0.15, linewidth = 0.2) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.4) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "Fraction of UMIs",
       subtitle = "Marker leakage: raw -> automatic rho -> fixed rho (one line per library)") +
  theme_asset() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_leak, run, "leakage_raw_auto_fixed", width = 8, height = 5)

rt <- melt(lib_tab[, c("library_id", ret_cols), with = FALSE], id.vars = "library_id",
           variable.name = "marker", value.name = "retention")
p_ret <- ggplot(rt, aes(marker, retention)) +
  geom_hline(yintercept = SX$retention_flag, linetype = 2, linewidth = 0.2) +
  geom_jitter(width = 0.15, size = 0.6) +
  labs(x = NULL, y = "Corrected / raw mean expression",
       subtitle = "Own-lineage markers should be retained (dashed = review threshold)") +
  theme_asset() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_ret, run, "own_marker_retention", width = 5, height = 3.5)

# ---- summary -------------------------------------------------

flagged <- lib_tab[nzchar(review_flag)]
summary_lines <- c(
  sprintf("Libraries: %d   cells: %s   rho: %.3f (fixed, all libraries)",
          nrow(lib_tab), format(sum(lib_tab$n_cells), big.mark = ","), SX$rho),
  sprintf("Counts removed per library: median %.1f%% (range %.1f-%.1f%%); automatic rho (04_1) removed %.1f%%",
          100 * median(lib_tab$frac_removed), 100 * min(lib_tab$frac_removed),
          100 * max(lib_tab$frac_removed), 100 * median(lib_tab$frac_removed_auto)),
  "Leakage, median over libraries (raw -> auto rho -> fixed rho):",
  sprintf("    %-20s %.4f -> %.4f -> %.4f", med$metric, med$raw, med$auto_rho, med$fixed_rho),
  "Own-marker retention, median (min):",
  sprintf("    %-40s %.3f (%.3f)", ret_cols,
          sapply(ret_cols, function(r) median(lib_tab[[r]], na.rm = TRUE)),
          sapply(ret_cols, function(r) min(lib_tab[[r]], na.rm = TRUE))),
  sprintf("Libraries flagged (own marker < %.2f): %d", SX$retention_flag, nrow(flagged)),
  if (nrow(flagged)) paste0("    ", flagged$library_id, " [", flagged$batch_id, "] ",
                            format(flagged$n_cells, big.mark = ","), " cells"),
  "",
  sprintf("Corrected counts: %s/", .rel_to_repo(OUT)))
writeLines(summary_lines, file.path(run$dir, "AMBIENT_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

finalize_run(run, status = "ok", verdict = sprintf(
  "rho %.2f fixed; %.1f%% of counts removed (median); keratin in immune %.4f -> %.4f; %d libraries flagged",
  SX$rho, 100 * median(lib_tab$frac_removed),
  med[metric == "keratin_in_immune"]$raw, med[metric == "keratin_in_immune"]$fixed_rho, nrow(flagged)))
