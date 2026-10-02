# ==============================================================
# 05_2_assess_integration_all_cells.R
#
# Integration checks and the owner's first figures, from the object saved
# by 05_1_integrate_all_cells.R (no recomputation of PCA/Harmony).
#
#   - neighbour mixing: for a random 100k cells, how much more often a
#     cell's 30 nearest neighbours share its pool / timepoint / patient
#     than expected by chance, in PCA vs Harmony space. Harmony should
#     pull pool towards 1 and leave timepoint and patient largely where
#     PCA had them;
#   - the same, SSc cells only, within each provisional lineage; and
#     within each SSc patient, whether their own timepoints stay apart
#     (one-pool vs multi-pool patients); pool among neighbours from
#     OTHER patients only (patients are nested in pools, so same-patient
#     neighbours inflate the plain pool score). PCA elbow;
#   - UMAP without vs with Harmony, same cells, by pool and timepoint;
#   - all cells coloured by provisional lineage, pool, timepoint, group,
#     and wasikowr's labels where inherited;
#   - Placebo cells (cohort `placebo`) by mRSS_category (Improver |
#     Stable | Worsened | Set_aside), M00 only and all timepoints, equal
#     cells per panel for the three main groups.
# (Per-cluster timepoint composition needs clusters: 05_3.)
#
# Colours are PROVISIONAL lineages, not the final annotation.
#
# PNGs use the cairo device: R on this server has no X11.
#
#   OMP_NUM_THREADS=8 jobs/run.sh analysis/05_sc_integration/05_2_assess_integration_all_cells.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")
source("R/theme_asset.R")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(BPCells)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
options(bitmapType = "cairo")

INT <- list(umap_neighbors = 30, umap_min_dist = 0.3, unintegrated_cells = 300000,
            plot_cells = 300000, mixing_cells = 100000, mixing_k = 30,
            panel_max_cells = 20000, strat_cells = 30000,
            within_pt_cells = 1000, within_pt_min_cells = 200, seed = 1)
CAT_LEVELS <- c("Improver", "Stable", "Worsened", "Set_aside")

run <- init_run(
  stage    = "05_sc_integration",
  run_name = "assess_integration_all_cells",
  params   = INT,
  notes    = "Integration checks and Placebo UMAPs by mRSS_category from the 05_1 object"
)
THREADS <- as.integer(Sys.getenv("OMP_NUM_THREADS", "4"))
set.seed(INT$seed)

obj_file <- file.path(OBJ, sprintf("%s_reference_integrated_all_cells.rds", FREEZE))
if (!file.exists(obj_file)) stop("05_1 object not found: ", obj_file, call. = FALSE)
obj <- readRDS(obj_file)
pc_dim <- ncol(Embeddings(obj, "harmony"))
for (v in c("timepoint", "batch_id", "Subject_ID", "group", "coarse_lineage",
            "wasikowr_celltype", "mRSS_category", "library_id")) {
  x <- as.character(obj@meta.data[[v]])
  x[!is.na(x) & !nzchar(x)] <- NA            # the manifest stores missing as "" (e.g. HC timepoint)
  obj@meta.data[[v]] <- x
}
message(sprintf("Loaded %s cells, %d Harmony dims", format(ncol(obj), big.mark = ","), pc_dim))

# unintegrated UMAP subset: from 05_1 if it saved one, otherwise computed here
# (only from THIS object's 05_1 run: a subset from an earlier embedding would be wrong)
un_file <- c(file.path(results_dir("05_sc_integration", "integrate_all_cells"),
                       "objects", "umap_unintegrated_subset.csv.gz"),
             file.path(run$objects, "umap_unintegrated_subset.csv.gz"))
un_file <- un_file[file.exists(un_file)][1]
if (!is.na(un_file)) {
  und <- fread(un_file)
  sub_cells <- und$cell_id; un <- as.matrix(und[, .(un_1, un_2)])
} else {
  message("No saved unintegrated UMAP; computing on a random subset")
  sub_cells <- sample(colnames(obj), min(INT$unintegrated_cells, ncol(obj)))
  un <- uwot::umap(Embeddings(obj, "pca")[sub_cells, seq_len(pc_dim)],
                   n_neighbors = INT$umap_neighbors, min_dist = INT$umap_min_dist,
                   n_threads = THREADS, seed = INT$seed)
  fwrite(data.table(cell_id = sub_cells, un_1 = un[, 1], un_2 = un[, 2]),
         file.path(run$objects, "umap_unintegrated_subset.csv.gz"))
}

# ---- 1. neighbour mixing, PCA vs Harmony ---------------------

mix_cells <- sample(colnames(obj), min(INT$mixing_cells, ncol(obj)))
mix_md <- obj@meta.data[mix_cells, c("batch_id", "timepoint", "Subject_ID", "group")]
mix_md$timepoint[is.na(mix_md$timepoint)] <- "HC"
knn <- function(X, k) {
  a <- RcppAnnoy::AnnoyEuclidean$new(ncol(X))
  for (i in seq_len(nrow(X))) a$addItem(i - 1, X[i, ])
  a$build(50)
  t(vapply(seq_len(nrow(X)), function(i) a$getNNsByItem(i - 1, k + 1)[-1] + 1, numeric(k)))
}
enrichment <- function(nn, lab) {
  p   <- prop.table(table(lab))
  obs <- rowMeans(matrix(lab[nn] == rep(lab, ncol(nn)), nrow = nrow(nn)))
  mean(obs) / mean(p[lab])
}
mix <- rbindlist(lapply(c("pca", "harmony"), function(r) {
  nn <- knn(Embeddings(obj, r)[mix_cells, seq_len(pc_dim)], INT$mixing_k)
  data.table(space = r,
             variable = c("pool (batch_id)", "timepoint", "patient (Subject_ID)"),
             same_label_enrichment = c(enrichment(nn, mix_md$batch_id),
                                       enrichment(nn, mix_md$timepoint),
                                       enrichment(nn, mix_md$Subject_ID)))
}))
save_table(mix, run, "neighbour_mixing_pca_vs_harmony")

# ---- 1b. stratified checks (added 2026-10-02) -----------------
# The global score above mixes three things: cell-type composition (a
# keratinocyte-rich library has keratinocyte neighbours), the HC pools
# (HC-only, so their separation is pool-driven), and real biology. These
# two checks remove the first two:
#   (i)  SSc cells only, WITHIN each provisional lineage: pool / patient /
#        timepoint enrichment, PCA vs Harmony;
#   (ii) WITHIN each SSc patient with >= 2 timepoints: do that patient's
#        own M00 / M03 / M06 cells stay distinguishable? Split by whether
#        the patient's libraries sit in one pool or span several. For
#        one-pool patients Harmony applies the same correction to every
#        timepoint, so this score should not move; for multi-pool patients
#        part of the timepoint difference is pool, and some drop is the
#        correction working, not signal lost.

# among neighbours from a different patient, the share from the same pool,
# over the share expected (that pool's share of other patients' cells)
pool_other_patients <- function(nn, pool, pat) {
  N <- length(pool); k <- ncol(nn)
  pp <- paste(pool, pat)
  expct <- (as.numeric(table(pool)[pool]) - as.numeric(table(pp)[pp])) /
           (N - as.numeric(table(pat)[pat]))
  other <- matrix(pat[nn] != rep(pat, k), nrow = N)
  same  <- matrix(pool[nn] == rep(pool, k), nrow = N)
  ok <- rowSums(other) > 0
  mean(rowSums(same & other)[ok] / rowSums(other)[ok]) / mean(expct[ok])
}

ssc <- obj@meta.data$group == "SSc" & !is.na(obj@meta.data$timepoint)
md_all <- obj@meta.data
strat <- rbindlist(lapply(sort(unique(md_all$coarse_lineage[ssc])), function(lin) {
  cells <- rownames(md_all)[ssc & md_all$coarse_lineage == lin]
  if (length(cells) < 2000) return(NULL)
  cells <- sample(cells, min(length(cells), INT$strat_cells))
  m <- md_all[cells, ]
  rbindlist(lapply(c("pca", "harmony"), function(r) {
    nn <- knn(Embeddings(obj, r)[cells, seq_len(pc_dim)], INT$mixing_k)
    data.table(lineage = lin, n_cells = length(cells), space = r,
               variable = c("pool (batch_id)", "timepoint", "patient (Subject_ID)",
                            "pool, other-patient neighbours"),
               same_label_enrichment = c(enrichment(nn, m$batch_id), enrichment(nn, m$timepoint),
                                         enrichment(nn, m$Subject_ID),
                                         pool_other_patients(nn, m$batch_id, m$Subject_ID)))
  }))
}))
save_table(strat, run, "neighbour_mixing_by_lineage_ssc")

pt <- as.data.table(md_all[ssc, c("Subject_ID", "timepoint", "batch_id")], keep.rownames = "cell")
pt_info <- pt[, .(n_tp = uniqueN(timepoint), n_pool = uniqueN(batch_id),
                  min_tp_cells = min(table(timepoint))), by = Subject_ID]
pt_use <- pt_info[n_tp >= 2 & min_tp_cells >= INT$within_pt_min_cells]
within <- rbindlist(lapply(pt_use$Subject_ID, function(s) {
  x <- pt[Subject_ID == s]
  cells <- x[, .(cell = sample(cell, min(.N, INT$within_pt_cells))), by = timepoint]$cell
  tp <- md_all[cells, "timepoint"]
  rbindlist(lapply(c("pca", "harmony"), function(r) {
    nn <- knn(Embeddings(obj, r)[cells, seq_len(pc_dim)], INT$mixing_k)
    data.table(subject_idx = match(s, pt_use$Subject_ID),   # index only, no IDs
               pools = ifelse(pt_use[Subject_ID == s]$n_pool > 1, "multi-pool", "one pool"),
               space = r, timepoint_enrichment = enrichment(nn, tp))
  }))
}))
save_table(within, run, "within_patient_timepoint_separation")
within_w <- dcast(within, subject_idx + pools ~ space, value.var = "timepoint_enrichment")
within_sum <- within_w[, .(patients = .N, pca_median = median(pca), harmony_median = median(harmony),
                           median_ratio = median(harmony / pca)), by = pools]

# PCA elbow (the lab's sd rule) and Harmony convergence, from the saved object
sv <- Stdev(obj, "pca")
save_table(data.table(pc = seq_along(sv), stdev = sv, ratio_to_previous = c(NA, sv[-1] / sv[-length(sv)])),
           run, "pca_stdev")
p_elb <- ggplot(data.table(pc = seq_along(sv), sd = sv), aes(pc, sd)) +
  geom_point(size = 0.8) + geom_vline(xintercept = pc_dim, linetype = 2, linewidth = 0.2) +
  labs(x = "PC", y = "Standard deviation", subtitle = sprintf("PCs used: %d (dashed)", pc_dim)) +
  theme_asset()
save_plot(p_elb, run, "pca_elbow", width = 4.5, height = 3)
hm_misc <- obj[["harmony"]]@misc
message("Harmony misc slots: ", paste(names(hm_misc), collapse = ", "))

# ---- 2. plots: integration checks -----------------------------

pc <- sample(colnames(obj), min(INT$plot_cells, ncol(obj)))
pd <- cbind(obj@meta.data[pc, ], Embeddings(obj, "umap")[pc, ])
pd$timepoint[is.na(pd$timepoint)] <- "HC"
umap_plot <- function(d, colour, title, pal = NULL, x = "umap_1", y = "umap_2") {
  p <- ggplot(d, aes(.data[[x]], .data[[y]], colour = .data[[colour]])) +
    geom_point(size = 0.05, alpha = 0.4, shape = 16) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(title = title, x = NULL, y = NULL, colour = NULL) +
    theme_asset() + theme(axis.text = element_blank(), axis.ticks = element_blank())
  if (!is.null(pal)) p <- p + scale_colour_manual(values = pal, na.value = "grey85")
  p
}
names(pd)[names(pd) %in% c("umap_1", "umap_2", "UMAP_1", "UMAP_2")] <- c("umap_1", "umap_2")
pal_time <- c(PAL_TIME, HC = "#7F7F7F")

ud <- data.frame(obj@meta.data[sub_cells, c("batch_id", "timepoint", "coarse_lineage")],
                 un_1 = un[, 1], un_2 = un[, 2],
                 umap_1 = Embeddings(obj, "umap")[sub_cells, 1],
                 umap_2 = Embeddings(obj, "umap")[sub_cells, 2])
ud$timepoint[is.na(ud$timepoint)] <- "HC"
p_side <- (umap_plot(ud, "batch_id", "PCA (no integration): pool", x = "un_1", y = "un_2") +
           umap_plot(ud, "batch_id", "Harmony on pool: pool")) /
          (umap_plot(ud, "timepoint", "PCA (no integration): timepoint", pal_time, "un_1", "un_2") +
           umap_plot(ud, "timepoint", "Harmony on pool: timepoint", pal_time)) +
          plot_layout(guides = "collect")
ggsave(file.path(run$plots, "integration_check_unintegrated_vs_harmony.png"), p_side,
       width = 12, height = 11, dpi = 200)

for (v in c("coarse_lineage", "batch_id", "timepoint", "group", "wasikowr_celltype")) {
  pal <- switch(v, timepoint = pal_time, group = PAL_GROUP, NULL)
  ggsave(file.path(run$plots, paste0("umap_all_cells_by_", v, ".png")),
         umap_plot(pd, v, paste("All cells, Harmony on pool:", v), pal),
         width = 8, height = 6.5, dpi = 200)
}

p_mix <- ggplot(mix, aes(variable, same_label_enrichment, fill = space)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.2) +
  labs(x = NULL, y = "Neighbours sharing the label / expected by chance",
       subtitle = "Pool should drop towards 1 after Harmony; timepoint and patient should not be forced to 1") +
  theme_asset()
save_plot(p_mix, run, "integration_check_neighbour_mixing", width = 5.5, height = 3.5)

# ---- 3. owner's figures: Placebo by mRSS_category ------------

pl_libs <- cohort_samples("placebo", manifest = read_manifest())$library_id
pl <- obj@meta.data[obj$library_id %in% pl_libs & !is.na(obj$mRSS_category), ]
pl$mRSS_category <- factor(pl$mRSS_category, levels = CAT_LEVELS)
emb_all <- Embeddings(obj, "umap")
bg <- data.frame(emb_all[sample(rownames(emb_all), 100000), ])
names(bg) <- c("umap_1", "umap_2")

panel_plot <- function(d, colour_by, title, rows = NULL) {
  # equal cells per panel (rows x category), capped
  # equalised across Improver / Stable / Worsened; Set_aside (very few
  # patients) shows up to the same number, i.e. all of its cells if fewer
  grp <- if (is.null(rows)) as.character(d$mRSS_category) else paste(d[[rows]], d$mRSS_category)
  main <- table(grp[d$mRSS_category != "Set_aside"])
  n_eq <- min(min(main), INT$panel_max_cells)
  keep <- unlist(lapply(split(rownames(d), grp), function(x) sample(x, min(length(x), n_eq))))
  e <- data.frame(d[keep, ], emb_all[keep, ])
  names(e)[(ncol(e) - 1):ncol(e)] <- c("umap_1", "umap_2")
  n_sub <- if (is.null(rows)) tapply(d$Subject_ID, d$mRSS_category, function(x) length(unique(x)))
           else NULL
  n_shown <- table(factor(d[keep, "mRSS_category"], levels = names(n_sub)))
  facet_lab <- if (is.null(rows)) {
    setNames(sprintf("%s\n(%d patients; %s cells shown)", names(n_sub), n_sub,
                     format(as.integer(n_shown), big.mark = ",")), names(n_sub))
  } else waiver()
  p <- ggplot(e, aes(umap_1, umap_2)) +
    geom_point(data = bg, colour = "grey92", size = 0.05, shape = 16) +
    geom_point(aes(colour = .data[[colour_by]]), size = 0.08, alpha = 0.6, shape = 16) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(title = title, x = NULL, y = NULL, colour = NULL,
         caption = "Provisional lineage colours (not final annotation). Equal cells per panel.") +
    theme_asset() + theme(axis.text = element_blank(), axis.ticks = element_blank())
  if (is.null(rows)) p + facet_wrap(~ mRSS_category, nrow = 1, drop = FALSE,
                                    labeller = labeller(mRSS_category = facet_lab))
  else p + facet_grid(reformulate("mRSS_category", rows), drop = FALSE)
}

m00 <- pl[pl$timepoint == "M00", ]
ggsave(file.path(run$plots, "placebo_M00_umap_by_mrss_category.png"),
       panel_plot(m00, "coarse_lineage", "Placebo, baseline (M00), by mRSS_category"),
       width = 14, height = 4.5, dpi = 250)
ggsave(file.path(run$plots, "placebo_all_timepoints_umap_by_mrss_category.png"),
       panel_plot(pl, "coarse_lineage", "Placebo, M00 / M03 / M06, by mRSS_category", rows = "timepoint"),
       width = 14, height = 11, dpi = 250)
if (any(!is.na(m00$wasikowr_celltype))) {
  ggsave(file.path(run$plots, "placebo_M00_umap_by_mrss_category_wasikowr_labels.png"),
         panel_plot(m00[!is.na(m00$wasikowr_celltype), ], "wasikowr_celltype",
                    "Placebo M00, wasikowr's broad labels (cells she labelled)"),
         width = 14, height = 4.5, dpi = 250)
}

# counts behind the panels (aggregate, no IDs)
panel_counts <- as.data.table(pl)[, .(cells = .N, patients = uniqueN(Subject_ID)),
                                  by = .(timepoint, mRSS_category)][order(timepoint, mRSS_category)]
save_table(panel_counts, run, "placebo_panel_counts")

# ---- summary ------------------------------------------------

mw <- dcast(mix, variable ~ space, value.var = "same_label_enrichment")
e_of <- function(v, sp) mw[[sp]][mw$variable == v]
summary_lines <- c(
  sprintf("Cells: %s   genes: %d   PCs: %d   Harmony on: batch_id",
          format(ncol(obj), big.mark = ","), nrow(obj), pc_dim),
  "Neighbour enrichment (same label among 30 neighbours / chance), PCA -> Harmony:",
  sprintf("  %-22s %.2f -> %.2f", mw$variable, mw$pca, mw$harmony),
  "Same, SSc cells only, within each provisional lineage (PCA -> Harmony):",
  {
    sw <- dcast(strat, lineage + n_cells + variable ~ space, value.var = "same_label_enrichment")
    sprintf("  %-13s %-22s %.2f -> %.2f", sw$lineage, sw$variable, sw$pca, sw$harmony)
  },
  "Within-patient timepoint separation (SSc patients with >= 2 timepoints; median enrichment):",
  sprintf("  %-10s %2d patients  PCA %.2f -> Harmony %.2f  (median per-patient ratio %.2f)",
          within_sum$pools, within_sum$patients, within_sum$pca_median,
          within_sum$harmony_median, within_sum$median_ratio),
  sprintf("PCA: %d PCs computed, %d used; sd of PC%d / PC1 = %.3f", length(sv), pc_dim, pc_dim, sv[pc_dim] / sv[1]),
  sprintf("Harmony misc slots: %s", if (length(hm_misc)) paste(names(hm_misc), collapse = ", ") else "none"),
  "Placebo panels (cells / patients):",
  sprintf("  %s %-10s %8s cells  %2d patients", panel_counts$timepoint, panel_counts$mRSS_category,
          format(panel_counts$cells, big.mark = ","), panel_counts$patients),
  "",
  "Review: plots/integration_check_*.png before reading the Placebo panels.")
writeLines(summary_lines, file.path(run$dir, "INTEGRATION_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

finalize_run(run, status = "ok", verdict = sprintf(
  "%s cells, %d PCs; pool enrichment %.2f -> %.2f, timepoint %.2f -> %.2f, patient %.2f -> %.2f",
  format(ncol(obj), big.mark = ","), pc_dim,
  e_of("pool (batch_id)", "pca"), e_of("pool (batch_id)", "harmony"),
  e_of("timepoint", "pca"), e_of("timepoint", "harmony"),
  e_of("patient (Subject_ID)", "pca"), e_of("patient (Subject_ID)", "harmony")))
