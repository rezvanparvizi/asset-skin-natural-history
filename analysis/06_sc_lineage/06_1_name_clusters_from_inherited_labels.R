# ==============================================================
# 06_1_name_clusters_from_inherited_labels.R
#
# Broad labels for ALL cells from wasikowr's broad labels, by rule B
# (docs/decisions.md 2026-10-03, recorded before this script first ran;
# replaces rule A's panel naming, which R0039 met at no resolution):
#
#   R0039's saved Louvain clusters (05_3, Harmony SNN, fixed grid) are
#   named by the reference cells they contain: our QC-pass cells carrying
#   `wasikowr_basectrl__celltype` (baseline + HC only). Per cluster, label =
#   most common reference label if n_ref >= 50 and purity >= 0.75; "mixed"
#   if n_ref >= 50 and purity < 0.75; "unlabelled" if n_ref < 50. Chosen =
#   the LOWEST resolution at which each of her 16 labels names >= 1 cluster
#   and mixed + unlabelled clusters hold <= 5% of cells. If none
#   qualifies, nothing is chosen (owner review).
#
# Cross-check, never used to assign labels: kNN transfer on the Harmony
# dims (k = 30, Annoy; a reference cell never counts itself). Reported:
# leave-one-out confusion on reference cells, cluster label vs kNN label,
# and kNN confidence by timepoint -- the check on the main assumption, that
# baseline labels describe M03/M06 cells.
#
# Barcode-level output goes to data/patient_level/ (save_patient_table).
# Nothing here looks at mRSS_category.
#
#   jobs/run.sh analysis/06_sc_lineage/06_1_name_clusters_from_inherited_labels.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")
source("R/theme_asset.R")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(BiocNeighbors)
  library(ggplot2)
})

CL <- list(clusters_run = "cluster_top_level", grid = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8),
           label_column = "wasikowr_basectrl__celltype", n_required = 16L,
           min_ref = 50, min_purity = 0.75, max_unresolved_frac = 0.05,
           knn_k = 30, low_conf = 0.5, threads = 4L, seed = 1)

run <- init_run(
  stage    = "06_sc_lineage",
  run_name = "name_clusters_from_inherited_labels",
  params   = CL,
  notes    = "Rule B (decisions.md 2026-10-03): R0039 clusters named by wasikowr's broad labels; kNN transfer as cross-check"
)
set.seed(CL$seed)

# ---- 1. inputs -----------------------------------------------

obj <- readRDS(file.path(OBJ, sprintf("%s_reference_integrated_all_cells.rds", FREEZE)))
emb <- Embeddings(obj, "harmony")
md <- data.table(cell_id = colnames(obj))
for (v in c("timepoint", "group")) {
  x <- as.character(obj@meta.data[[v]]); x[!is.na(x) & !nzchar(x)] <- NA
  md[[v]] <- x
}
rm(obj); invisible(gc())
stopifnot(all(stats::na.omit(md$timepoint) %in% c("M00", "M03", "M06")),
          all(md$group %in% c("SSc", "HC")))
md[, tp := fifelse(group == "HC", "HC", timepoint)]
message(sprintf("Loaded %s cells, %d Harmony dims", format(nrow(md), big.mark = ","), ncol(emb)))

clu_file <- file.path(results_dir("05_sc_integration", CL$clusters_run, cohort = "reference"),
                      "objects", "top_level_clusters.csv.gz")
res_cols <- paste0("res_", CL$grid)
clu <- fread(clu_file, select = c("cell_id", res_cols), colClasses = "character")
m <- match(md$cell_id, clu$cell_id)
if (anyNA(m) || nrow(clu) != nrow(md))
  stop("R0039 clusters do not cover the object's cells exactly: ", sum(is.na(m)), " missing", call. = FALSE)
clu <- clu[m]

lab_file <- file.path(METADATA, "inherited_labels", sprintf("%s_upstream_labels.csv.gz", FREEZE))
inh <- fread(lab_file, select = c("cell_id", CL$label_column))
m <- match(md$cell_id, inh$cell_id)
if (anyNA(m)) stop(sum(is.na(m)), " object cells missing from the inherited label table", call. = FALSE)
ref_label <- as.character(inh[[CL$label_column]][m])
ref_label[!is.na(ref_label) & !nzchar(ref_label)] <- NA
rm(inh); invisible(gc())
is_ref <- !is.na(ref_label)
required <- sort(unique(ref_label[is_ref]))
if (length(required) != CL$n_required)
  stop("expected ", CL$n_required, " reference labels, found ", length(required), call. = FALSE)
message(sprintf("Reference cells: %s (%s labels)", format(sum(is_ref), big.mark = ","), length(required)))

# ---- 2. rule B over the grid ---------------------------------

ssc <- md$group == "SSc" & !is.na(md$timepoint)

name_by_ref <- function(cl) {
  d <- data.table(cluster = cl, lab = ref_label, tp = md$tp, ssc = ssc)
  tot <- d[, .(cells = .N, n_ref = sum(!is.na(lab)),
               ssc_cells = sum(ssc),
               share_M00 = sum(tp == "M00", na.rm = TRUE) / max(sum(ssc), 1),
               share_M03 = sum(tp == "M03", na.rm = TRUE) / max(sum(ssc), 1),
               share_M06 = sum(tp == "M06", na.rm = TRUE) / max(sum(ssc), 1),
               share_HC  = sum(tp == "HC", na.rm = TRUE) / .N), by = cluster]
  cnt <- d[!is.na(lab), .N, by = .(cluster, lab)][order(cluster, -N)]
  top <- cnt[, .(top = lab[1], purity = N[1] / sum(N),
                 second = if (.N > 1) lab[2] else NA_character_,
                 second_share = if (.N > 1) N[2] / sum(N) else 0), by = cluster]
  tab <- merge(tot, top, by = "cluster", all.x = TRUE)
  tab[, label := fifelse(n_ref < CL$min_ref, "unlabelled",
                         fifelse(purity >= CL$min_purity, top, "mixed"))]
  tab[, ref_per_1000 := round(1000 * n_ref / cells, 1)]
  tab[order(-cells)]
}

sweep <- list(); per_cluster <- list(); labs <- list()
for (r in CL$grid) {
  cl <- clu[[paste0("res_", r)]]
  tab <- name_by_ref(cl)
  covered <- intersect(required, tab$label)
  unres <- sum(tab$cells[tab$label %in% c("mixed", "unlabelled")]) / sum(tab$cells)
  sweep[[as.character(r)]] <- data.table(
    resolution = r, n_clusters = nrow(tab), required_covered = length(covered),
    required_missing = paste(setdiff(required, covered), collapse = "; "),
    mixed_cells_frac = sum(tab$cells[tab$label == "mixed"]) / sum(tab$cells),
    unlabelled_cells_frac = sum(tab$cells[tab$label == "unlabelled"]) / sum(tab$cells),
    meets_rule_B = length(covered) == length(required) && unres <= CL$max_unresolved_frac)
  per_cluster[[as.character(r)]] <- cbind(resolution = r, tab)
  labs[[paste0("label_res_", r)]] <- tab$label[match(cl, tab$cluster)]
  message(sprintf("  resolution %.2f: %d clusters, %d/%d labels, mixed %.1f%%, unlabelled %.1f%%",
                  r, nrow(tab), length(covered), length(required),
                  100 * sweep[[as.character(r)]]$mixed_cells_frac,
                  100 * sweep[[as.character(r)]]$unlabelled_cells_frac))
}
sweep <- rbindlist(sweep); per_cluster <- rbindlist(per_cluster)
save_table(sweep, run, "rule_B_sweep")
save_table(per_cluster, run, "clusters_by_resolution")
chosen <- if (any(sweep$meets_rule_B)) min(sweep$resolution[sweep$meets_rule_B]) else NA_real_
message("Chosen resolution: ", ifelse(is.na(chosen), "NONE (owner review)", chosen))

# ---- 3. kNN transfer (cross-check) ---------------------------

k <- CL$knn_k
q <- queryKNN(X = emb[is_ref, , drop = FALSE], query = emb, k = k + 1,
              BNPARAM = AnnoyParam(), num.threads = CL$threads, get.distance = FALSE)
idx <- q$index; rm(q)
ref_row <- integer(nrow(md)); ref_row[is_ref] <- seq_len(sum(is_ref))
self_hit <- idx == ref_row                      # column-wise recycling: row i vs ref_row[i]
has_self <- rowSums(self_hit) > 0
nb <- idx[, seq_len(k), drop = FALSE]
if (any(has_self)) {
  x <- idx[has_self, , drop = FALSE]; x[self_hit[has_self, , drop = FALSE]] <- NA
  nb[has_self, ] <- t(apply(x, 1, function(v) v[!is.na(v)][seq_len(k)]))
}
stopifnot(!anyNA(nb))
rm(idx, self_hit); invisible(gc())
nb_lab <- matrix(ref_label[is_ref][nb], nrow = nrow(nb))
cnt <- vapply(required, function(l) rowSums(nb_lab == l), numeric(nrow(nb_lab)))
best <- max.col(cnt, ties.method = "first")
knn_label <- required[best]
knn_conf <- cnt[cbind(seq_len(nrow(cnt)), best)] / k
rm(nb, nb_lab, cnt); invisible(gc())
message(sprintf("kNN: self excluded for %s of %s reference cells",
                format(sum(has_self & is_ref), big.mark = ","), format(sum(is_ref), big.mark = ",")))

loo <- data.table(truth = ref_label[is_ref], knn = knn_label[is_ref])[, .N, by = .(truth, knn)]
loo[, share_of_truth := N / sum(N), by = truth]
save_table(loo[order(truth, -N)], run, "knn_leave_one_out_confusion")
recall <- loo[truth == knn, .(truth, recall = share_of_truth)][order(recall)]

conf_tp <- data.table(tp = md$tp, knn_label = knn_label, knn_conf = knn_conf)[
  , .(cells = .N, median_conf = stats::median(knn_conf), low_conf_share = mean(knn_conf < CL$low_conf)),
  by = .(tp, knn_label)][order(knn_label, tp)]
save_table(conf_tp, run, "knn_confidence_by_timepoint_and_label")
conf_tp_all <- data.table(tp = md$tp, knn_conf = knn_conf)[
  , .(cells = .N, median_conf = stats::median(knn_conf), low_conf_share = mean(knn_conf < CL$low_conf)), by = tp][order(tp)]
save_table(conf_tp_all, run, "knn_confidence_by_timepoint")

agree <- rbindlist(lapply(CL$grid, function(r) {
  l <- labs[[paste0("label_res_", r)]]
  data.table(cluster_label = l, knn_label = knn_label)[
    , .(cells = .N, agree_with_knn = mean(cluster_label == knn_label)), by = cluster_label][
    , resolution := r][]
}))
setcolorder(agree, c("resolution", "cluster_label", "cells", "agree_with_knn"))
save_table(agree[order(resolution, -cells)], run, "cluster_label_vs_knn_label")

# ---- 4. per-cell output (barcode-level) ----------------------

cells_out <- data.table(cell_id = md$cell_id, as.data.table(labs),
                        knn_label = knn_label, knn_conf = round(knn_conf, 3))
if (!is.na(chosen)) {
  cells_out[, cluster_chosen := clu[[paste0("res_", chosen)]]]
  cells_out[, broad_label := labs[[paste0("label_res_", chosen)]]]
}
save_patient_table(cells_out, run, "broad_labels_per_cell")

# ---- 5. plots ------------------------------------------------

lp <- copy(loo)[, `:=`(truth = factor(truth, levels = required), knn = factor(knn, levels = required))]
p_loo <- ggplot(lp, aes(knn, truth, fill = share_of_truth)) + geom_tile() +
  geom_text(aes(label = ifelse(share_of_truth >= 0.01, sprintf("%.2f", share_of_truth), "")), size = 1.8) +
  scale_fill_gradient(low = "white", high = "#2F4F7F", limits = c(0, 1)) +
  labs(x = "kNN label (k = 30, self excluded)", y = "wasikowr label", fill = "share of\nher label",
       subtitle = "Leave-one-out kNN on Harmony dims, reference cells (baseline + HC)") +
  theme_asset() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_loo, run, "knn_leave_one_out_confusion", width = 8, height = 6)

if (!is.na(chosen)) {
  pc <- per_cluster[resolution == chosen]
  pc[, late_share := share_M03 + share_M06]
  p_pc <- ggplot(pc, aes(late_share, purity, size = cells, colour = label)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = CL$min_purity, linetype = 2) +
    labs(x = "M03 + M06 share of SSc cells in cluster", y = "purity of reference labels",
         size = "cells", colour = NULL,
         title = sprintf("Clusters at resolution %s (rule B)", chosen),
         subtitle = "Clusters to the right are named by few baseline cells; check them first") +
    theme_asset()
  save_plot(p_pc, run, "cluster_purity_vs_late_share", width = 9, height = 6)
}

# ---- summary -------------------------------------------------

summary_lines <- c(
  sprintf("Cells: %s   reference cells: %s   labels: %d   k = %d",
          format(nrow(md), big.mark = ","), format(sum(is_ref), big.mark = ","), length(required), k),
  "Rule B sweep (n_ref >= 50, purity >= 0.75; mixed + unlabelled <= 5% of cells):",
  sprintf("  res %-5s %3d clusters  labels %2d/%d  mixed %5.1f%%  unlabelled %5.1f%%  meets: %s  missing: %s",
          sweep$resolution, sweep$n_clusters, sweep$required_covered, length(required),
          100 * sweep$mixed_cells_frac, 100 * sweep$unlabelled_cells_frac, sweep$meets_rule_B,
          sweep$required_missing),
  sprintf("CHOSEN resolution: %s", ifelse(is.na(chosen), "NONE — no resolution met rule B; owner review", chosen)),
  "",
  "kNN leave-one-out recall per label (reference cells):",
  sprintf("  %-22s %.3f", recall$truth, recall$recall),
  "",
  "kNN confidence by timepoint (median; share < 0.5):",
  sprintf("  %-4s %9s cells  median %.2f  low %.3f", conf_tp_all$tp,
          format(conf_tp_all$cells, big.mark = ","), conf_tp_all$median_conf, conf_tp_all$low_conf_share))
if (!is.na(chosen)) {
  a <- agree[resolution == chosen]
  summary_lines <- c(summary_lines, "",
    sprintf("Cluster label vs kNN label at resolution %s (share of cells agreeing):", chosen),
    sprintf("  %-22s %9s cells  %.3f", a$cluster_label, format(a$cells, big.mark = ","), a$agree_with_knn))
}
writeLines(summary_lines, file.path(run$dir, "BROAD_LABEL_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

finalize_run(run, status = "ok", verdict = if (is.na(chosen))
  "no resolution met rule B — owner review (see rule_B_sweep)" else
  sprintf("rule B chose resolution %s (%d clusters, %d/%d labels); kNN agrees with cluster labels for %.1f%% of cells",
          chosen, sweep[resolution == chosen]$n_clusters, sweep[resolution == chosen]$required_covered,
          length(required),
          100 * mean(labs[[paste0("label_res_", chosen)]] == knn_label)))
