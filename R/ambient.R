# ==============================================================
# R/ambient.R — helpers for ambient-RNA correction (stage 04)
#
# Functions only. Used by 04_2_correct_ambient_fixed_rho.R. The same
# helpers appear inline in 04_1_soupx_per_library.R, which is left as it
# ran for R0020.
# ==============================================================

# Marker panels for the provisional lineage of a SoupX cluster, and the
# genes used for the leakage / retention checks (decisions.md 2026-10-01).
AMBIENT_PANELS <- list(
  Keratinocyte = c("KRT5", "KRT14", "KRT1", "KRT10", "KRT15", "PERP"),
  Fibroblast   = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
  Immune       = c("PTPRC", "CD3E", "CD2", "CD74", "LYZ", "CD68", "CD14", "MS4A1"),
  Plasma       = c("JCHAIN", "MZB1", "XBP1"),
  Endothelial  = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  Mural        = c("ACTA2", "TAGLN", "RGS5", "MYH11"),
  Melanocyte   = c("PMEL", "MLANA", "TYRP1"),
  Gland        = c("DCD", "SCGB2A2", "MUCL1", "PIP")
)
AMBIENT_LEAK <- list(keratin  = c("KRT1", "KRT5", "KRT10", "KRT14"),
                     collagen = c("COL1A1", "COL1A2", "COL3A1"),
                     hb       = c("HBB", "HBA1", "HBA2"),
                     ig       = c("IGKC", "JCHAIN", "IGHG1", "IGHA1", "IGHM"))
AMBIENT_OWN <- c(Fibroblast = "COL1A1", Keratinocyte = "KRT14", Immune = "PTPRC")

# Gene-expression matrix from a CellRanger h5, gene symbols as rownames.
read_h5_gex <- function(f) {
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

# Flex: the raw (droplet) matrix lists every reference gene, the filtered
# (cell) matrix only the probe-panel genes. Keep the panel genes; every one
# must exist in the raw matrix.
align_cell_droplet_genes <- function(toc, tod, lib) {
  common <- intersect(rownames(toc), rownames(tod))
  if (length(common) < nrow(toc)) {
    stop(lib, ": ", nrow(toc) - length(common), " cell-matrix genes absent from ",
         "the raw matrix (filtered ", nrow(toc), ", raw ", nrow(tod), ")")
  }
  list(toc = toc[common, , drop = FALSE], tod = tod[common, , drop = FALSE])
}

# Share of a cell group's UMIs that fall in a gene set.
umi_frac_of <- function(m, genes, cells) {
  g <- intersect(genes, rownames(m))
  if (!length(g) || !length(cells)) return(NA_real_)
  sum(m[g, cells, drop = FALSE]) / max(sum(m[, cells, drop = FALSE]), 1)
}

# Mean count of one gene over a cell group.
gene_mean_of <- function(m, gene, cells) {
  if (!gene %in% rownames(m) || !length(cells)) return(NA_real_)
  mean(m[gene, cells])
}

# Leakage (raw vs corrected) and own-marker retention for one library.
# `cl` has barcode and coarse_lineage for the library's cells.
ambient_checks <- function(raw, corr, cl, lib) {
  by_lin <- function(l) cl$barcode[cl$coarse_lineage %in% l]
  imm <- by_lin(c("Immune", "Plasma")); fib <- by_lin("Fibroblast")
  ker <- by_lin("Keratinocyte"); nonimm <- setdiff(cl$barcode, imm)
  pair <- function(x, f) c(raw = f(raw, x), corrected = f(corr, x))
  leak <- rbind(
    keratin_in_immune  = pair(imm, function(z, x) umi_frac_of(z, AMBIENT_LEAK$keratin, x)),
    keratin_in_fibro   = pair(fib, function(z, x) umi_frac_of(z, AMBIENT_LEAK$keratin, x)),
    collagen_in_immune = pair(imm, function(z, x) umi_frac_of(z, AMBIENT_LEAK$collagen, x)),
    collagen_in_kerat  = pair(ker, function(z, x) umi_frac_of(z, AMBIENT_LEAK$collagen, x)),
    hb_all_cells       = pair(cl$barcode, function(z, x) umi_frac_of(z, AMBIENT_LEAK$hb, x)),
    ig_outside_immune  = pair(nonimm, function(z, x) umi_frac_of(z, AMBIENT_LEAK$ig, x)))
  own <- sapply(names(AMBIENT_OWN), function(lin) {
    x <- by_lin(lin); r <- gene_mean_of(raw, AMBIENT_OWN[[lin]], x)
    c <- gene_mean_of(corr, AMBIENT_OWN[[lin]], x)
    if (is.na(r) || r == 0) NA_real_ else c / r
  })
  data.table::data.table(
    library_id = lib,
    metric = c(rownames(leak), paste0("retention_", AMBIENT_OWN, "_in_", names(AMBIENT_OWN))),
    raw = c(leak[, "raw"], rep(NA, length(AMBIENT_OWN))),
    corrected = c(leak[, "corrected"], rep(NA, length(AMBIENT_OWN))),
    ratio = c(leak[, "corrected"] / leak[, "raw"], own),
    n_cells_group = c(length(imm), length(fib), length(imm), length(ker), nrow(cl),
                      length(nonimm), lengths(lapply(names(AMBIENT_OWN), by_lin))))
}
