# ==============================================================
# 04_1_soupx_per_library.R
#
# Ambient RNA correction with SoupX, one library at a time, on the cells
# that passed QC (03_1). Rule fixed before running
# (docs/decisions.md, 2026-10-01, "SoupX rule"):
#
#   soup     : that library's own empty droplets (sample_raw h5)
#   clusters : log-normalise, 2,000 HVG (ambient genes NOT excluded),
#              30 PCs, Louvain resolution 0.8 — identical for every library
#   rho      : autoEstCont; accepted if 0.01 <= rho <= 0.30, otherwise the
#              median accepted rho of the same pool
#   output   : adjustCounts(roundToInt = TRUE) -> integer corrected counts
#
# Skin-specific checks, reported side by side for every library:
#   (a) rho vs pool / timepoint / arm
#   (b) marker leakage before -> after (keratin in immune and fibroblast
#       clusters; collagen in immune and keratinocyte clusters; HBB;
#       immunoglobulin outside immune clusters)
#   (c) retention of each lineage's own markers (over-correction check)
#   (d) the soup's top genes per library
# Clusters are given a coarse lineage from marker panels (cluster level,
# so ambient RNA in single cells does not decide it).
#
# Input : data/raw/<run>/<lib>/sample_{filtered,raw}_feature_bc_matrix.h5
#         results/<freeze>/reference/03_sc_qc/cell_qc_per_library/objects/cell_qc.csv.gz
# Output: data/bpcells/<freeze>/counts_raw/<lib>/     raw integer counts, QC-pass cells
#         data/bpcells/<freeze>/counts_soupx/<lib>/   corrected integer counts
#         objects/soupx_cells.csv.gz                  cell -> SoupX cluster, coarse lineage
#         tables/ plots/ SOUPX_SUMMARY.txt
#
# Nothing here looks at mRSS_category.
#
#   jobs/run.sh analysis/04_sc_ambient/04_1_soupx_per_library.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")
source("R/theme_asset.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(BPCells)
  library(Seurat)
  library(SoupX)
})

SX <- list(n_hvg = 2000, n_pcs = 30, resolution = 0.8,
           rho_min = 0.01, rho_max = 0.30, round_to_int = TRUE,
           retention_flag = 0.70, seed = 1)

PANELS <- list(
  Keratinocyte = c("KRT5", "KRT14", "KRT1", "KRT10", "KRT15", "PERP"),
  Fibroblast   = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
  Immune       = c("PTPRC", "CD3E", "CD2", "CD74", "LYZ", "CD68", "CD14", "MS4A1"),
  Plasma       = c("JCHAIN", "MZB1", "XBP1"),
  Endothelial  = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  Mural        = c("ACTA2", "TAGLN", "RGS5", "MYH11"),
  Melanocyte   = c("PMEL", "MLANA", "TYRP1"),
  Gland        = c("DCD", "SCGB2A2", "MUCL1", "PIP")
)
LEAK <- list(keratin = c("KRT1", "KRT5", "KRT10", "KRT14"),
             collagen = c("COL1A1", "COL1A2", "COL3A1"),
             hb = c("HBB", "HBA1", "HBA2"),
             ig = c("IGKC", "JCHAIN", "IGHG1", "IGHA1", "IGHM"))
OWN <- c(Fibroblast = "COL1A1", Keratinocyte = "KRT14", Immune = "PTPRC")

run <- init_run(
  stage    = "04_sc_ambient",
  run_name = "soupx_per_library",
  params   = c(SX, list(panels = PANELS, leakage_genes = LEAK, own_markers = as.list(OWN))),
  notes    = "SoupX per library: fixed clustering, autoEstCont with bounds and pool fallback"
)

THREADS <- as.integer(Sys.getenv("OMP_NUM_THREADS", "4"))
man <- freeze_libraries()
qc_file <- file.path(results_dir("03_sc_qc", "cell_qc_per_library", cohort = "reference"),
                     "objects", "cell_qc.csv.gz")
if (!file.exists(qc_file)) stop("cell QC table not found: ", qc_file, call. = FALSE)
qc <- fread(qc_file, select = c("library_id", "barcode", "qc_pass"))
qc <- qc[qc_pass == TRUE]

OUT_RAW   <- file.path(BPCELLS, FREEZE, "counts_raw")
OUT_SOUPX <- file.path(BPCELLS, FREEZE, "counts_soupx")
dir.create(OUT_RAW,   recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_SOUPX, recursive = TRUE, showWarnings = FALSE)

# ---- helpers ------------------------------------------------

read_h5 <- function(f) {
  m  <- BPCells::open_matrix_10x_hdf5(f, feature_type = "Gene Expression")
  nm <- as.character(rhdf5::h5read(f, "matrix/features/name"))
  ft <- as.character(rhdf5::h5read(f, "matrix/features/feature_type"))
  rhdf5::h5closeAll()
  nm <- nm[ft == "Gene Expression"]
  stopifnot(length(nm) == nrow(m))
  m <- as(m, "dgCMatrix")
  rownames(m) <- make.unique(nm)
  m
}

# Flex: the raw (droplet) matrix lists every reference gene (~39k), the
# filtered (cell) matrix only the probe-panel genes (~18k). Every
# cell-matrix gene must exist in the raw matrix; the raw-only genes are
# outside the panel and are dropped. Differences are reported.
align_genes <- function(toc, tod, lib) {
  common <- intersect(rownames(toc), rownames(tod))
  if (length(common) < nrow(toc)) {
    stop(lib, ": ", nrow(toc) - length(common), " cell-matrix genes absent from ",
         "the raw matrix (filtered ", nrow(toc), ", raw ", nrow(tod), ")")
  }
  diff <- data.table(library_id = lib, n_genes_filtered = nrow(toc),
                     n_genes_raw = nrow(tod), n_genes_shared = length(common),
                     only_in_filtered = paste(setdiff(rownames(toc), common), collapse = ";"),
                     n_only_in_raw = nrow(tod) - length(common),
                     raw_only_counts = sum(tod[setdiff(rownames(tod), common), , drop = FALSE]))
  list(toc = toc[common, , drop = FALSE], tod = tod[common, , drop = FALSE], diff = diff)
}

frac_of <- function(m, genes, cells) {
  g <- intersect(genes, rownames(m))
  if (!length(g) || !length(cells)) return(NA_real_)
  sum(m[g, cells, drop = FALSE]) / max(sum(m[, cells, drop = FALSE]), 1)
}
mean_of <- function(m, gene, cells) {
  if (!gene %in% rownames(m) || !length(cells)) return(NA_real_)
  mean(m[gene, cells])
}

# Fixed per-library clustering and coarse lineage per cluster.
cluster_library <- function(toc) {
  set.seed(SX$seed)
  s <- CreateSeuratObject(toc, min.cells = 0, min.features = 0)
  s <- NormalizeData(s, verbose = FALSE)
  s <- FindVariableFeatures(s, nfeatures = SX$n_hvg, verbose = FALSE)
  s <- ScaleData(s, verbose = FALSE)
  npc <- min(SX$n_pcs, ncol(s) - 1, length(VariableFeatures(s)) - 1)
  s <- RunPCA(s, npcs = npc, verbose = FALSE)
  s <- FindNeighbors(s, dims = seq_len(npc), verbose = FALSE)
  s <- FindClusters(s, resolution = SX$resolution, random.seed = SX$seed, verbose = FALSE)
  cl <- as.character(Idents(s))
  # coarse lineage: mean normalised panel expression per cluster, argmax
  dat <- LayerData(s, layer = "data")
  sc  <- sapply(PANELS, function(g) {
    g <- intersect(g, rownames(dat))
    if (!length(g)) return(rep(0, ncol(dat)))
    Matrix::colMeans(dat[g, , drop = FALSE])
  })
  agg <- apply(sc, 2, function(v) tapply(v, cl, mean))
  if (is.null(dim(agg))) agg <- matrix(agg, nrow = 1, dimnames = list(unique(cl), names(PANELS)))
  lin <- setNames(colnames(agg)[max.col(agg, ties.method = "first")], rownames(agg))
  list(cluster = cl, lineage = unname(lin[cl]),
       pcs = Embeddings(s, "pca"))
}

# ---- pass 1: clusters + automatic rho per library ------------

estimate_library <- function(i) {
  lib <- man$library_id[i]
  d   <- file.path(RAW, man$raw_dir[i])
  toc <- read_h5(file.path(d, "sample_filtered_feature_bc_matrix.h5"))
  tod <- read_h5(file.path(d, "sample_raw_feature_bc_matrix.h5"))
  keep <- qc$barcode[qc$library_id == lib]
  toc <- toc[, colnames(toc) %in% keep, drop = FALSE]
  stopifnot(ncol(toc) == length(keep))
  al <- align_genes(toc, tod, lib); toc <- al$toc; tod <- al$tod

  cl <- cluster_library(toc)
  sc <- SoupChannel(tod, toc, calcSoupProfile = TRUE)
  sc <- setClusters(sc, setNames(cl$cluster, colnames(toc)))
  est <- tryCatch({
    s2 <- autoEstCont(sc, doPlot = FALSE, verbose = FALSE, forceAccept = TRUE)
    s2$metaData$rho[1]
  }, error = function(e) NA_real_)

  soup <- sc$soupProfile
  top  <- head(soup[order(-soup$est), , drop = FALSE], 20)
  message(sprintf("  pass 1 done: %s (%d cells, rho_auto %s)", lib, ncol(toc),
                  format(round(est, 3))))
  list(
    gene_diff = al$diff,
    lib = lib,
    cells = data.table(library_id = lib, barcode = colnames(toc),
                       soupx_cluster = cl$cluster, coarse_lineage = cl$lineage),
    est = data.table(library_id = lib, n_cells = ncol(toc),
                     n_clusters = length(unique(cl$cluster)),
                     n_empty_droplets = sum(Matrix::colSums(tod) > 0) - ncol(toc),
                     rho_auto = est),
    soup_top = data.table(library_id = lib, rank = seq_len(nrow(top)),
                          gene = rownames(top), soup_frac = top$est)
  )
}

message(sprintf("Pass 1: %d libraries, %d threads", nrow(man), THREADS))
p1 <- parallel::mclapply(seq_len(nrow(man)), function(i)
  tryCatch(estimate_library(i), error = function(e) e),
  mc.cores = THREADS, mc.preschedule = FALSE)
err <- vapply(p1, inherits, logical(1), "error")
if (any(err)) {
  msg <- paste(man$library_id[err], vapply(p1[err], conditionMessage, ""), sep = ": ", collapse = "\n")
  finalize_run(run, status = "failed", verdict = sprintf("pass 1 failed for %d libraries", sum(err)))
  stop("Pass 1 failed:\n", msg, call. = FALSE)
}

est <- rbindlist(lapply(p1, `[[`, "est"))
est <- merge(est, as.data.table(man[, c("library_id", "batch_id", "group", "arm", "timepoint")]),
             by = "library_id")
est[, rho_accepted := !is.na(rho_auto) & rho_auto >= SX$rho_min & rho_auto <= SX$rho_max]
pool_med <- est[rho_accepted == TRUE, .(pool_rho = stats::median(rho_auto)), by = batch_id]
est <- merge(est, pool_med, by = "batch_id", all.x = TRUE)
if (any(!est$rho_accepted & is.na(est$pool_rho))) {
  finalize_run(run, status = "failed", verdict = "a pool has no library with an accepted rho")
  stop("No accepted rho in pool(s): ",
       paste(unique(est$batch_id[!est$rho_accepted & is.na(est$pool_rho)]), collapse = ", "))
}
est[, rho_used := ifelse(rho_accepted, rho_auto, pool_rho)]
est[, rho_source := ifelse(rho_accepted, "autoEstCont",
                           ifelse(is.na(rho_auto), "pool_median (auto failed)",
                                  "pool_median (auto out of bounds)"))]
cells <- rbindlist(lapply(p1, `[[`, "cells"))
soup_top <- rbindlist(lapply(p1, `[[`, "soup_top"))
gene_diff <- rbindlist(lapply(p1, `[[`, "gene_diff"))
rm(p1); invisible(gc())

# ---- pass 2: correct, write matrices, leakage metrics --------

correct_library <- function(i) {
  lib <- man$library_id[i]
  d   <- file.path(RAW, man$raw_dir[i])
  toc <- read_h5(file.path(d, "sample_filtered_feature_bc_matrix.h5"))
  tod <- read_h5(file.path(d, "sample_raw_feature_bc_matrix.h5"))
  cl  <- cells[library_id == lib]
  toc <- toc[, cl$barcode, drop = FALSE]
  al  <- align_genes(toc, tod, lib); toc <- al$toc; tod <- al$tod
  rho <- est$rho_used[est$library_id == lib]

  sc <- SoupChannel(tod, toc, calcSoupProfile = TRUE)
  sc <- setClusters(sc, setNames(cl$soupx_cluster, cl$barcode))
  sc <- setContaminationFraction(sc, rho, forceAccept = TRUE)
  out <- adjustCounts(sc, roundToInt = SX$round_to_int, verbose = 0)
  stopifnot(all(out@x == round(out@x)), identical(dim(out), dim(toc)))
  message(sprintf("  pass 2 done: %s (rho %.3f)", lib, rho))

  cid <- paste(lib, colnames(toc), sep = "_")
  write_one <- function(m, root) {
    colnames(m) <- cid
    dir <- file.path(root, lib)
    if (dir.exists(dir)) unlink(dir, recursive = TRUE)
    BPCells::write_matrix_dir(BPCells::convert_matrix_type(as(m, "IterableMatrix"), "uint32_t"),
                              dir, overwrite = TRUE)
  }
  write_one(toc, OUT_RAW)
  write_one(out, OUT_SOUPX)

  by_lin <- function(l) cl$barcode[cl$coarse_lineage %in% l]
  imm <- by_lin(c("Immune", "Plasma")); fib <- by_lin("Fibroblast")
  ker <- by_lin("Keratinocyte"); nonimm <- setdiff(cl$barcode, imm)
  m <- function(x, f) setNames(c(f(toc, x), f(out, x)), c("raw", "soupx"))
  leak <- rbind(
    keratin_in_immune   = m(imm, function(z, x) frac_of(z, LEAK$keratin, x)),
    keratin_in_fibro    = m(fib, function(z, x) frac_of(z, LEAK$keratin, x)),
    collagen_in_immune  = m(imm, function(z, x) frac_of(z, LEAK$collagen, x)),
    collagen_in_kerat   = m(ker, function(z, x) frac_of(z, LEAK$collagen, x)),
    hb_all_cells        = m(cl$barcode, function(z, x) frac_of(z, LEAK$hb, x)),
    ig_outside_immune   = m(nonimm, function(z, x) frac_of(z, LEAK$ig, x)))
  own <- sapply(names(OWN), function(lin) {
    x <- by_lin(lin); r <- mean_of(toc, OWN[[lin]], x); c <- mean_of(out, OWN[[lin]], x)
    if (is.na(r) || r == 0) NA_real_ else c / r
  })
  data.table(library_id = lib, metric = c(rownames(leak), paste0("retention_", OWN, "_in_", names(OWN))),
             raw = c(leak[, "raw"], rep(NA, length(OWN))),
             soupx = c(leak[, "soupx"], rep(NA, length(OWN))),
             ratio = c(leak[, "soupx"] / leak[, "raw"], own),
             n_cells_group = c(length(imm), length(fib), length(imm), length(ker),
                               nrow(cl), length(nonimm), lengths(lapply(names(OWN), by_lin))),
             total_raw = sum(toc), total_soupx = sum(out))
}

message("Pass 2: correcting and writing matrices")
p2 <- parallel::mclapply(seq_len(nrow(man)), function(i)
  tryCatch(correct_library(i), error = function(e) e),
  mc.cores = THREADS, mc.preschedule = FALSE)
err <- vapply(p2, inherits, logical(1), "error")
if (any(err)) {
  msg <- paste(man$library_id[err], vapply(p2[err], conditionMessage, ""), sep = ": ", collapse = "\n")
  finalize_run(run, status = "failed", verdict = sprintf("pass 2 failed for %d libraries", sum(err)))
  stop("Pass 2 failed:\n", msg, call. = FALSE)
}
leak <- rbindlist(p2)

# ---- tables -------------------------------------------------

cells[, cell_id := paste(library_id, barcode, sep = "_")]
fwrite(cells[, .(cell_id, library_id, barcode, soupx_cluster, coarse_lineage)],
       file.path(run$objects, "soupx_cells.csv.gz"))

tot <- unique(leak[, .(library_id, total_raw, total_soupx)])
est <- merge(est, tot, by = "library_id")
est[, frac_counts_removed := 1 - total_soupx / total_raw]
ret <- dcast(leak[grepl("^retention_", metric)], library_id ~ metric, value.var = "ratio")
est <- merge(est, ret, by = "library_id")
ret_cols <- grep("^retention_", names(est), value = TRUE)
est[, review_flag := trimws(paste(
  ifelse(rho_source != "autoEstCont", "rho_from_pool", ""),
  ifelse(apply(as.matrix(.SD) < SX$retention_flag, 1, any, na.rm = TRUE),
         "own_marker_loss", ""))), .SDcols = ret_cols]
setorder(est, batch_id, library_id)
save_table(est, run, "soupx_per_library")
save_table(leak[, !c("total_raw", "total_soupx")], run, "leakage_per_library")
save_table(soup_top, run, "soup_top_genes")
save_table(gene_diff, run, "gene_lists_filtered_vs_raw")

# (a) rho vs design, SSc libraries
ssc <- est[group == "SSc"]
design_p <- rbindlist(lapply(c("batch_id", "timepoint", "arm"), function(v)
  data.table(variable = v,
             kruskal_p_rho = signif(stats::kruskal.test(ssc$rho_used, factor(ssc[[v]]))$p.value, 3))))
save_table(design_p, run, "rho_vs_design")

# (d) soup composition by panel
soup_panel <- soup_top[, .(frac = sum(soup_frac)), by = .(library_id,
  panel = vapply(gene, function(g) {
    hit <- names(PANELS)[vapply(PANELS, function(p) g %in% p, logical(1))]
    if (grepl("^KRT", g)) "Keratinocyte" else if (grepl("^COL[0-9]", g)) "Fibroblast"
    else if (length(hit)) hit[1] else "other"
  }, ""))]
save_table(dcast(soup_panel, library_id ~ panel, value.var = "frac", fill = 0),
           run, "soup_top20_by_lineage")

# ---- plots --------------------------------------------------

p_rho <- ggplot(ssc, aes(timepoint, rho_used, colour = arm)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(aes(shape = rho_source),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 0.8) +
  scale_colour_manual(values = PAL_ARM) +
  labs(x = NULL, y = "Contamination fraction (rho)", colour = NULL, shape = NULL,
       subtitle = "(a) rho must not track timepoint or arm") +
  theme_asset()
save_plot(p_rho, run, "soupx_rho_by_timepoint_arm", width = 5, height = 3.5)

p_pool <- ggplot(est, aes(batch_id, rho_used, colour = group)) +
  geom_jitter(width = 0.15, size = 0.8) +
  scale_colour_manual(values = PAL_GROUP) +
  labs(x = NULL, y = "rho", colour = NULL, subtitle = "(a) rho per library, by pool") +
  theme_asset() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_pool, run, "soupx_rho_by_pool", width = 6, height = 3.5)

lk <- leak[!grepl("^retention_", metric)]
lk_long <- melt(lk, id.vars = c("library_id", "metric"), measure.vars = c("raw", "soupx"),
                variable.name = "counts", value.name = "fraction")
p_leak <- ggplot(lk_long, aes(counts, fraction)) +
  geom_line(aes(group = library_id), alpha = 0.2, linewidth = 0.2) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.4) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "Fraction of UMIs", subtitle = "(b) marker leakage, before -> after SoupX (one line per library)") +
  theme_asset()
save_plot(p_leak, run, "soupx_leakage_before_after", width = 7, height = 4.5)

rt <- melt(est[, c("library_id", ret_cols), with = FALSE], id.vars = "library_id",
           variable.name = "marker", value.name = "retention")
p_ret <- ggplot(rt, aes(marker, retention)) +
  geom_hline(yintercept = SX$retention_flag, linetype = 2, linewidth = 0.2) +
  geom_jitter(width = 0.15, size = 0.6) +
  labs(x = NULL, y = "Corrected / raw mean expression",
       subtitle = "(c) own-lineage markers should be retained (dashed = review threshold)") +
  theme_asset() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_ret, run, "soupx_own_marker_retention", width = 5, height = 3.5)

# ---- summary ------------------------------------------------

med_leak <- lk[, .(raw = stats::median(raw, na.rm = TRUE),
                   soupx = stats::median(soupx, na.rm = TRUE)), by = metric]
flagged <- est[nzchar(review_flag)]
top_soup <- head(soup_top[rank <= 5, .N, by = gene][order(-N)], 10)
summary_lines <- c(
  sprintf("Libraries: %d   cells: %s", nrow(est), format(nrow(cells), big.mark = ",")),
  sprintf("Genes: %s probe-panel genes used; raw matrices carry %s extra non-panel genes with %s total counts",
          paste(unique(gene_diff$n_genes_shared), collapse = "/"),
          paste(unique(gene_diff$n_only_in_raw), collapse = "/"),
          format(sum(gene_diff$raw_only_counts), big.mark = ",")),
  sprintf("rho: median %.3f, range %.3f-%.3f; from autoEstCont %d, pool median %d",
          stats::median(est$rho_used), min(est$rho_used), max(est$rho_used),
          sum(est$rho_source == "autoEstCont"), sum(est$rho_source != "autoEstCont")),
  sprintf("Counts removed: median %.1f%% per library",
          100 * stats::median(est$frac_counts_removed)),
  "(a) rho vs design (SSc, Kruskal):",
  paste0("    ", design_p$variable, ": p = ", design_p$kruskal_p_rho),
  "(b) leakage, median over libraries (raw -> SoupX):",
  sprintf("    %-20s %.4f -> %.4f", med_leak$metric, med_leak$raw, med_leak$soupx),
  "(c) own-marker retention, median:",
  sprintf("    %-40s %.3f", ret_cols, sapply(ret_cols, function(r) stats::median(est[[r]], na.rm = TRUE))),
  "(d) genes most often in a library's top-5 soup genes:",
  paste0("    ", top_soup$gene, " (", top_soup$N, " libraries)"),
  sprintf("Libraries flagged for review: %d", nrow(flagged)),
  if (nrow(flagged)) paste0("    ", flagged$library_id, " [", flagged$batch_id, "] rho ",
                            sprintf("%.3f", flagged$rho_used), "  ", flagged$review_flag),
  "",
  "Review with the owner before the corrected counts are used (stage 05).")
writeLines(summary_lines, file.path(run$dir, "SOUPX_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

finalize_run(run, status = "ok", verdict = sprintf(
  "rho median %.3f (%.3f-%.3f); %d of %d auto-accepted; rho vs timepoint p=%s; %d libraries flagged for review",
  stats::median(est$rho_used), min(est$rho_used), max(est$rho_used),
  sum(est$rho_source == "autoEstCont"), nrow(est),
  design_p[variable == "timepoint", kruskal_p_rho], nrow(flagged)))
