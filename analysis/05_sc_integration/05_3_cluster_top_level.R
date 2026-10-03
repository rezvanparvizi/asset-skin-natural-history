# ==============================================================
# 05_3_cluster_top_level.R
#
# Top-level clusters of ALL cells and the resolution chosen by rule A
# (docs/decisions.md 2026-10-02, recorded before this script first ran):
#
#   SNN graph on the Harmony PCs of 05_1 `integrate_all_cells`; Louvain at
#   a fixed grid. Each cluster is named by marker panels
#   (metadata/gene_sets/top_level_panels.yml): per gene, the cluster mean of
#   log-normalised expression, z-scored across clusters; panel score = mean
#   z of its genes; the top panel is the name if it beats the second by
#   >= 1.0, else "ambiguous". Chosen = the LOWEST resolution at which every
#   required label is the top call of >= 1 cluster and ambiguous clusters
#   hold <= 5% of cells. If none qualifies, nothing is chosen.
#
# Reported with the choice (never part of it):
#   - timepoint composition per cluster, Harmony clusters vs clusters on
#     the unintegrated PCA at the same resolution (SSc cells);
#   - clusters dominated by one patient or one pool;
#   - comparison with wasikowr's broad labels and jarnagin's immune labels
#     (06_0 table), at the chosen resolution and at 0.1 (her top level).
#
# Nothing here looks at mRSS_category.
#
#   OMP_NUM_THREADS=8 jobs/run.sh analysis/05_sc_integration/05_3_cluster_top_level.R \
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
})
options(bitmapType = "cairo",
        future.globals.maxSize = 8 * 1024^3)   # FindNeighbors ships ~1 GB of globals at 1.2M cells (R0038 failed on the 500 MB default)

CL <- list(grid = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), k_param = 20, louvain_n_start = 3,
           margin = 1.0, max_ambiguous_frac = 0.05, compare_resolution = 0.1,
           tp_enrich_flag = 2, tp_flag_min_cells = 500,
           dominance_flag = 0.5, dominance_min_cells = 200,
           panel_file = "metadata/gene_sets/top_level_panels.yml",
           label_columns = c("wasikowr_basectrl__celltype", "jf_immune__celltype",
                             "wasikowr_object__RNA_snn_res.0.1"),   # her res-0.1 clusters: all timepoints, 1.19M of our cells
           plot_cells = 300000, seed = 1)

run <- init_run(
  stage    = "05_sc_integration",
  run_name = "cluster_top_level",
  params   = CL,
  notes    = "Top-level Louvain sweep on Harmony PCs; resolution by rule A (decisions.md 2026-10-02)"
)
set.seed(CL$seed)

# ---- 1. object, panels ---------------------------------------

obj <- readRDS(file.path(OBJ, sprintf("%s_reference_integrated_all_cells.rds", FREEZE)))
pc_dim <- ncol(Embeddings(obj, "harmony"))
for (v in c("timepoint", "batch_id", "Subject_ID", "group")) {
  x <- as.character(obj@meta.data[[v]]); x[!is.na(x) & !nzchar(x)] <- NA
  obj@meta.data[[v]] <- x
}
message(sprintf("Loaded %s cells, %d Harmony dims", format(ncol(obj), big.mark = ","), pc_dim))

pan <- yaml::read_yaml(CL$panel_file)
panels <- c(lapply(pan$required, `[[`, "genes"), lapply(pan$reported, `[[`, "genes"))
required <- names(pan$required)
genes <- unique(unlist(panels))
missing <- setdiff(genes, rownames(obj))
if (length(missing)) {
  finalize_run(run, status = "failed", verdict = paste("panel genes missing:", paste(missing, collapse = ", ")))
  stop("Panel genes not in the object: ", paste(missing, collapse = ", "), call. = FALSE)
}
M <- as(LayerData(obj, assay = "RNA", layer = "data")[genes, ], "dgCMatrix")   # genes x cells

name_clusters <- function(cl) {
  f <- factor(cl)
  X <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)))
  n <- as.numeric(table(f))
  mu <- as.matrix(M %*% X) %*% diag(1 / n, nrow = length(n))     # genes x clusters
  z  <- t(apply(mu, 1, function(r) { s <- sd(r); if (!is.finite(s) || s == 0) r * 0 else (r - mean(r)) / s }))
  sc <- t(sapply(panels, function(g) colMeans(z[g, , drop = FALSE])))   # panels x clusters
  colnames(sc) <- levels(f)
  o <- apply(sc, 2, order, decreasing = TRUE)
  top <- rownames(sc)[o[1, ]]; second <- rownames(sc)[o[2, ]]
  margin <- sc[cbind(o[1, ], seq_len(ncol(sc)))] - sc[cbind(o[2, ], seq_len(ncol(sc)))]
  list(table = data.table(cluster = levels(f), cells = n, top = top, second = second,
                          margin = margin, label = ifelse(margin >= CL$margin, top, "ambiguous")),
       scores = sc)
}

# ---- 2. graph, sweep, rule A ---------------------------------

obj <- FindNeighbors(obj, reduction = "harmony", dims = seq_len(pc_dim), k.param = CL$k_param,
                     graph.name = c("harmony_nn", "harmony_snn"), verbose = FALSE)
sweep <- list(); per_cluster <- list(); cl_cols <- list()
for (r in CL$grid) {
  t0 <- Sys.time()
  obj <- FindClusters(obj, graph.name = "harmony_snn", resolution = r, algorithm = 1,
                      n.start = CL$louvain_n_start, random.seed = CL$seed, verbose = FALSE)
  cl <- as.character(obj$seurat_clusters)
  nm <- name_clusters(cl)
  tab <- nm$table
  covered <- intersect(required, tab$label)
  amb <- sum(tab$cells[tab$label == "ambiguous"]) / sum(tab$cells)
  sweep[[as.character(r)]] <- data.table(
    resolution = r, n_clusters = nrow(tab), required_covered = length(covered),
    required_missing = paste(setdiff(required, covered), collapse = "; "),
    ambiguous_cells_frac = amb,
    reported_found = paste(intersect(names(pan$reported), tab$label), collapse = "; "),
    meets_rule_A = length(covered) == length(required) && amb <= CL$max_ambiguous_frac)
  per_cluster[[as.character(r)]] <- cbind(resolution = r, tab)
  cl_cols[[paste0("res_", r)]] <- cl
  message(sprintf("  resolution %.2f: %d clusters, %d/%d required, ambiguous %.1f%% (%.1f min)",
                  r, nrow(tab), length(covered), length(required), 100 * amb,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
sweep <- rbindlist(sweep); per_cluster <- rbindlist(per_cluster)
save_table(sweep, run, "rule_A_sweep")
save_table(per_cluster, run, "clusters_by_resolution")
chosen <- if (any(sweep$meets_rule_A)) min(sweep$resolution[sweep$meets_rule_A]) else NA_real_
message("Chosen resolution: ", ifelse(is.na(chosen), "NONE (owner review)", chosen))

lab_at <- function(r) {
  cl <- cl_cols[[paste0("res_", r)]]
  tab <- per_cluster[resolution == r]
  tab$label[match(cl, tab$cluster)]
}
cells_out <- data.table(cell_id = colnames(obj), as.data.table(cl_cols))
for (r in unique(c(chosen, CL$compare_resolution))) if (!is.na(r)) cells_out[[paste0("label_res_", r)]] <- lab_at(r)
fwrite(cells_out, file.path(run$objects, "top_level_clusters.csv.gz"))   # barcode-level, no subject IDs

# ---- 3. checks at the chosen resolution ----------------------

md <- obj@meta.data
ssc <- md$group == "SSc" & !is.na(md$timepoint)
tp_overall <- prop.table(table(md$timepoint[ssc]))

tp_comp <- function(cl, space) {
  d <- data.table(cluster = cl[ssc], timepoint = md$timepoint[ssc])
  w <- dcast(d[, .N, by = .(cluster, timepoint)], cluster ~ timepoint, value.var = "N", fill = 0)
  tps <- setdiff(names(w), "cluster")
  w[, ssc_cells := rowSums(.SD), .SDcols = tps]
  for (t in tps) w[[paste0("enrich_", t)]] <- (w[[t]] / w$ssc_cells) / as.numeric(tp_overall[t])
  ecols <- paste0("enrich_", tps)
  w[, max_enrich := do.call(pmax, .SD), .SDcols = ecols]
  w[, flag := ssc_cells >= CL$tp_flag_min_cells & max_enrich >= CL$tp_enrich_flag]
  cbind(space = space, w)
}

dominance <- function(cl) {
  d <- data.table(cluster = cl, patient = md$Subject_ID, pool = md$batch_id)
  d[, .(cells = .N,
        max_patient_share = max(table(patient, useNA = "no")) / sum(!is.na(patient)),
        n_patients = uniqueN(patient[!is.na(patient)]),
        max_pool_share = max(table(pool)) / .N, n_pools = uniqueN(pool)), by = cluster][
    , flag := cells >= CL$dominance_min_cells &
        (max_patient_share > CL$dominance_flag | max_pool_share > CL$dominance_flag)][order(cluster)]
}

ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  s <- sum(choose(tab, 2)); sa <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
  e <- sa * sb / choose(n, 2)
  (s - e) / ((sa + sb) / 2 - e)
}

if (!is.na(chosen)) {
  cl_h <- cl_cols[[paste0("res_", chosen)]]
  obj <- FindNeighbors(obj, reduction = "pca", dims = seq_len(pc_dim), k.param = CL$k_param,
                       graph.name = c("pca_nn", "pca_snn"), verbose = FALSE)
  obj <- FindClusters(obj, graph.name = "pca_snn", resolution = chosen, algorithm = 1,
                      n.start = CL$louvain_n_start, random.seed = CL$seed, verbose = FALSE)
  cl_p <- as.character(obj$seurat_clusters)
  fwrite(data.table(cell_id = colnames(obj), pca_cluster = cl_p),
         file.path(run$objects, "pca_clusters_at_chosen.csv.gz"))
  tpc <- rbind(tp_comp(cl_h, "harmony"), tp_comp(cl_p, "pca"), fill = TRUE)
  save_table(tpc, run, "timepoint_composition_per_cluster")
  # where did timepoint-enriched PCA clusters go after Harmony?
  fl <- tpc[space == "pca" & flag == TRUE, cluster]
  went <- rbindlist(lapply(fl, function(k) {
    idx <- which(cl_p == k & ssc)
    dest <- sort(table(cl_h[idx]), decreasing = TRUE)
    hk <- names(dest)[1]
    data.table(pca_cluster = k, ssc_cells = length(idx), top_harmony_cluster = hk,
               share_in_top = as.numeric(dest[1]) / length(idx),
               pca_max_enrich = tpc[space == "pca" & cluster == k, max_enrich],
               harmony_max_enrich = tpc[space == "harmony" & cluster == hk, max_enrich])
  }))
  save_table(went, run, "timepoint_enriched_pca_clusters_after_harmony")
  save_table(dominance(cl_h), run, "cluster_patient_pool_dominance")
  ari_hp <- ari(cl_h, cl_p)
} else {
  tpc <- data.table(); went <- data.table(); ari_hp <- NA_real_
}

# ---- 4. comparison with the colleagues' labels (report only) --

lab_file <- file.path(METADATA, "inherited_labels", sprintf("%s_upstream_labels.csv.gz", FREEZE))
inh <- fread(lab_file, select = c("cell_id", CL$label_columns))
inh <- inh[match(colnames(obj), cell_id)]
cmp <- list(); xt <- list()
for (r in unique(c(chosen, CL$compare_resolution))) {
  if (is.na(r)) next
  ours <- lab_at(r); clu <- cl_cols[[paste0("res_", r)]]
  for (col in CL$label_columns) {
    theirs <- as.character(inh[[col]]); ok <- !is.na(theirs) & nzchar(theirs)
    is_cluster_id <- grepl("RNA_snn_res|seurat_clusters", col)   # unnamed clusters: ARI and top_our_label only
    x <- data.table(theirs = theirs[ok], ours = ours[ok])
    a <- ari(clu[ok], theirs[ok])
    cmp[[paste(r, col)]] <- x[, .(cells = .N, agree = if (is_cluster_id) NA_real_ else mean(ours == theirs),
                                  top_our_label = names(sort(table(ours), decreasing = TRUE))[1]),
                              by = theirs][, `:=`(resolution = r, source = col, ari_clusters_vs_theirs = a)][]
    xt[[paste(r, col)]] <- x[, .N, by = .(theirs, ours)][, `:=`(resolution = r, source = col)][]
  }
}
cmp <- rbindlist(cmp); xt <- rbindlist(xt)
setcolorder(cmp, c("resolution", "source", "theirs", "cells", "agree", "top_our_label"))
save_table(cmp, run, "agreement_with_colleague_labels")
save_table(xt, run, "crosstab_with_colleague_labels")

# ---- 5. plots ------------------------------------------------

if (!is.na(chosen)) {
  sc <- name_clusters(cl_cols[[paste0("res_", chosen)]])$scores
  hm <- as.data.table(as.table(sc)); setnames(hm, c("panel", "cluster", "z"))
  hm[, cluster := factor(cluster, levels = colnames(sc)[order(as.integer(colnames(sc)))])]
  p_hm <- ggplot(hm, aes(cluster, panel, fill = z)) + geom_tile() +
    scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#B63B3B") +
    labs(x = sprintf("cluster (resolution %s)", chosen), y = NULL, fill = "panel z") +
    theme_asset() + theme(axis.text.x = element_text(angle = 90, size = 5))
  save_plot(p_hm, run, "panel_scores_chosen", width = 10, height = 5)

  pc <- sample(colnames(obj), min(CL$plot_cells, ncol(obj)))
  pd <- data.frame(Embeddings(obj, "umap")[pc, ], label = lab_at(chosen)[match(pc, colnames(obj))])
  names(pd)[1:2] <- c("umap_1", "umap_2")
  p_um <- ggplot(pd, aes(umap_1, umap_2, colour = label)) +
    geom_point(size = 0.05, alpha = 0.4, shape = 16) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(title = sprintf("Top-level labels, resolution %s (rule A)", chosen), x = NULL, y = NULL, colour = NULL) +
    theme_asset() + theme(axis.text = element_blank(), axis.ticks = element_blank())
  ggsave(file.path(run$plots, "umap_top_level_labels_chosen.png"), p_um, width = 9, height = 7, dpi = 200)
}
if (nrow(xt)) {
  xp <- xt[source == CL$label_columns[1]]
  xp[, frac := N / sum(N), by = .(resolution, theirs)]
  p_xt <- ggplot(xp, aes(ours, theirs, fill = frac)) + geom_tile() +
    facet_wrap(~ paste("resolution", resolution)) +
    scale_fill_gradient(low = "white", high = "#2F4F7F") +
    labs(x = "ours", y = "wasikowr", fill = "share of\nher label") +
    theme_asset() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_xt, run, "crosstab_wasikowr_labels", width = 11, height = 5)
}

# ---- summary -------------------------------------------------

fmtp <- function(x) sprintf("%.1f%%", 100 * x)
summary_lines <- c(
  sprintf("Cells: %s   Harmony dims: %d   k = %d   Louvain n.start = %d",
          format(ncol(obj), big.mark = ","), pc_dim, CL$k_param, CL$louvain_n_start),
  "Rule A sweep (required = wasikowr's 16 labels):",
  sprintf("  res %-5s %3d clusters  required %2d/%d  ambiguous %6s  meets: %-5s missing: %s  reported found: %s",
          sweep$resolution, sweep$n_clusters, sweep$required_covered, length(required),
          fmtp(sweep$ambiguous_cells_frac), sweep$meets_rule_A, sweep$required_missing, sweep$reported_found),
  sprintf("CHOSEN resolution: %s", ifelse(is.na(chosen), "NONE — no resolution met rule A; owner review", chosen)),
  if (!is.na(chosen)) c(
    sprintf("Harmony vs unintegrated-PCA clusters at %s: ARI %.3f", chosen, ari_hp),
    sprintf("Timepoint-enriched clusters (SSc, >= %dx overall, >= %d cells): PCA %d, Harmony %d",
            CL$tp_enrich_flag, CL$tp_flag_min_cells, sum(tpc[space == "pca"]$flag), sum(tpc[space == "harmony"]$flag)),
    if (nrow(went)) sprintf("  PCA cluster %s (%d cells, x%.1f) -> Harmony cluster %s (%.0f%% of its cells, x%.1f)",
                            went$pca_cluster, went$ssc_cells, went$pca_max_enrich, went$top_harmony_cluster,
                            100 * went$share_in_top, went$harmony_max_enrich)),
  "Agreement with colleague labels (share of their cells that we call the same; report only):",
  if (nrow(cmp)) sprintf("  res %-4s %-28s %-20s %7s cells  agree %6s  (ours mostly: %s)",
                         cmp$resolution, cmp$source, cmp$theirs, format(cmp$cells, big.mark = ","),
                         fmtp(cmp$agree), cmp$top_our_label),
  if (nrow(cmp)) sprintf("  ARI clusters vs labels: %s",
                         paste(unique(sprintf("res %s %s %.3f", cmp$resolution, cmp$source, cmp$ari_clusters_vs_theirs)),
                               collapse = "; ")))
writeLines(summary_lines, file.path(run$dir, "CLUSTER_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

finalize_run(run, status = "ok", verdict = if (is.na(chosen))
  "no resolution met rule A — owner review (see rule_A_sweep)" else
  sprintf("rule A chose resolution %s (%d clusters); timepoint-enriched clusters PCA %d / Harmony %d",
          chosen, sweep[resolution == chosen]$n_clusters,
          sum(tpc[space == "pca"]$flag), sum(tpc[space == "harmony"]$flag)))
