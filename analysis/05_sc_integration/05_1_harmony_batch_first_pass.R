# ==============================================================
# 05_1_harmony_batch_first_pass.R
#
# First-pass Layer-1 embedding of ALL cells (all arms, timepoints, HC):
# SoupX-corrected counts -> LogNormalize -> 2,000 HVG -> PCA -> Harmony
# on batch_id (pool) ONLY -> UMAP. Approved plan: docs/decisions.md,
# 2026-10-01, "Integration (stage 05)".
#
# This script computes and saves only: the object (PCA, Harmony, UMAP)
# to data/objects/, the Harmony UMAP and an unintegrated UMAP of a random
# subset to objects/. The integration checks and the owner's Placebo
# figures are made by 05_2_plot_first_pass_integration.R from these.
#
# No clustering here, so no resolution is chosen with mRSS_category in
# view (rule 3).
#
#   OMP_NUM_THREADS=8 jobs/run.sh analysis/05_sc_integration/05_1_harmony_batch_first_pass.R \
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
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

INT <- list(n_hvg = 2000, max_pcs = 50, pc_sd_ratio = 0.995, min_pcs = 15,
            harmony_by = "batch_id", umap_neighbors = 30, umap_min_dist = 0.3,
            unintegrated_cells = 300000, plot_cells = 300000,
            mixing_cells = 100000, mixing_k = 30, panel_max_cells = 20000, seed = 1)
CAT_LEVELS <- c("Improver", "Stable", "Worsened", "Set_aside")

run <- init_run(
  stage    = "05_sc_integration",
  run_name = "harmony_batch_first_pass",
  params   = INT,
  notes    = "All cells; SoupX counts; Harmony on batch_id only; mixing checks; Placebo UMAPs by mRSS_category"
)
THREADS <- as.integer(Sys.getenv("OMP_NUM_THREADS", "4"))
set.seed(INT$seed)

# ---- 1. counts: merge per-library SoupX matrices on disk -----

man <- freeze_libraries()
sx_cells <- fread(file.path(results_dir("04_sc_ambient", "soupx_per_library", cohort = "reference"),
                            "objects", "soupx_cells.csv.gz"))
sx_dir <- file.path(BPCELLS, FREEZE, "counts_soupx")
miss <- man$library_id[!dir.exists(file.path(sx_dir, man$library_id))]
if (length(miss)) stop("SoupX matrices missing for: ", paste(miss, collapse = ", "), call. = FALSE)

mats  <- lapply(man$library_id, function(l) open_matrix_dir(file.path(sx_dir, l)))
genes <- rownames(mats[[1]])
stopifnot(all(vapply(mats, function(m) identical(rownames(m), genes), logical(1))))
merged_dir <- file.path(BPCELLS, FREEZE, "merged_soupx")
if (dir.exists(merged_dir)) unlink(merged_dir, recursive = TRUE)
counts <- write_matrix_dir(do.call(cbind, mats), merged_dir)
rm(mats)
stopifnot(!anyDuplicated(colnames(counts)), setequal(colnames(counts), sx_cells$cell_id))
message(sprintf("Merged: %d genes x %s cells", nrow(counts), format(ncol(counts), big.mark = ",")))

# ---- 2. cell metadata ---------------------------------------

md <- sx_cells[match(colnames(counts), cell_id), .(cell_id, library_id, coarse_lineage)]
md <- cbind(md, as.data.table(man[match(md$library_id, man$library_id),
                                  c("Subject_ID", "group", "arm", "timepoint", "batch_id", "design_flag")]))
des <- read_design()
md[, mRSS_category := des$mRSS_category[match(Subject_ID, des$Subject_ID)]]
up_file <- file.path(METADATA, "inherited_labels", sprintf("%s_upstream_labels.csv.gz", FREEZE))
if (file.exists(up_file)) {
  up <- fread(up_file, select = c("cell_id", "wasikowr_basectrl__celltype"))
  md[, wasikowr_celltype := up$wasikowr_basectrl__celltype[match(cell_id, up$cell_id)]]
} else md[, wasikowr_celltype := NA_character_]
md <- as.data.frame(md); rownames(md) <- md$cell_id

# ---- 3. normalise, HVG, PCA ---------------------------------

obj <- CreateSeuratObject(counts = counts, meta.data = md)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = INT$n_hvg, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = INT$max_pcs, seed.use = INT$seed, verbose = FALSE)
sv  <- Stdev(obj, "pca")
idx <- which(sv[-1] > INT$pc_sd_ratio * sv[-length(sv)])          # lab rule (reembed())
pc_dim <- max(INT$min_pcs, if (length(idx)) min(idx) + 1 else length(sv))
message("PCs used: ", pc_dim)

# ---- 4. Harmony on pool, UMAP --------------------------------

obj <- RunHarmony(obj, group.by.vars = INT$harmony_by, reduction.use = "pca",
                  dims.use = seq_len(pc_dim), reduction.save = "harmony", verbose = FALSE)
um <- uwot::umap(Embeddings(obj, "harmony")[, seq_len(pc_dim)],
                 n_neighbors = INT$umap_neighbors, min_dist = INT$umap_min_dist,
                 n_threads = THREADS, seed = INT$seed)
dimnames(um) <- list(colnames(obj), c("umap_1", "umap_2"))
obj[["umap"]] <- CreateDimReducObject(um, key = "umap_", assay = "RNA")
rm(um)

# unintegrated UMAP on a random subset, for the side-by-side in 05_2
sub_cells <- sample(colnames(obj), min(INT$unintegrated_cells, ncol(obj)))
un <- uwot::umap(Embeddings(obj, "pca")[sub_cells, seq_len(pc_dim)],
                 n_neighbors = INT$umap_neighbors, min_dist = INT$umap_min_dist,
                 n_threads = THREADS, seed = INT$seed)
fwrite(data.table(cell_id = sub_cells, un_1 = un[, 1], un_2 = un[, 2]),
       file.path(run$objects, "umap_unintegrated_subset.csv.gz"))

# scale.data (dense, ~19 GB at this size) is not kept; counts stay on disk
saveRDS(DietSeurat(obj, layers = c("counts", "data"), dimreducs = c("pca", "harmony", "umap")),
        file.path(OBJ, sprintf("%s_reference_harmony_first_pass.rds", FREEZE)))
emb <- data.table(cell_id = colnames(obj),
                  UMAP_1 = Embeddings(obj, "umap")[, 1], UMAP_2 = Embeddings(obj, "umap")[, 2])
fwrite(emb, file.path(run$objects, "umap_harmony.csv.gz"))

# ---- done: checks and figures are 05_2 ----------------------

finalize_run(run, status = "ok", verdict = sprintf(
  "%s cells, %d PCs, Harmony on batch_id; object and UMAPs saved (checks/plots: 05_2)",
  format(ncol(obj), big.mark = ","), pc_dim))
