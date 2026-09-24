# ==============================================================
# pipelines/pseudobulk_de.R
#
# Config-driven pseudobulk differential expression, producing the same
# output contract your colleague's pipeline uses so results are directly
# comparable:
#
#   <run_dir>/
#     run_config.yml  session_info.txt  git_sha.txt
#     objects/  colData.rds  dds.rds  pseudo_filtered.rds  res.rds
#               res_df.rds   vsd.rds
#     plots/    volcano.pdf  heatmap_top_up_genes.pdf
#               top_gene_norm_counts.pdf
#       QC/     PCA, dispersion, sample correlation, cells-per-donor,
#               libsize-vs-genes, target-cell-proportion, MA,
#               top500 variable-gene heatmap
#     tables/   DESeq2_all_results.csv  DESeq2_results_annotated.csv
#               pseudobulk_counts.csv   colData.csv
#               donor_sanity_check.csv
#
# Run it:
#   jobs/run.sh pipelines/pseudobulk_de.R \
#     --freeze freeze01 --cohort placebo --labelset labelset01 \
#     --config config/de_runs/fib_improver_vs_non_b.yml
#
# ---------------------------------------------------------------
# STATUS: SKELETON. The scaffolding (config handling, aggregation,
# filtering, output contract, provenance) is complete and functional.
# The three model engines are laid out but NOT validated against your
# data — read them, decide on the design, and check the model matrix
# before trusting any p-value. Marked TODO(review) below.
# ==============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(yaml)
})

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")

# ---- 1. config ------------------------------------------------

cfg_path <- Sys.getenv("ASSET_CONFIG", unset = NA_character_)
if (is.na(cfg_path) || !nzchar(cfg_path)) {
  stop("Set ASSET_CONFIG, or pass --config to jobs/run.sh", call. = FALSE)
}
cfg <- yaml::read_yaml(cfg_path)

need <- c("run_name", "cell_selection", "aggregation", "model")
miss <- setdiff(need, names(cfg))
if (length(miss)) {
  stop("Config is missing: ", paste(miss, collapse = ", "), call. = FALSE)
}

run <- init_run(
  stage    = "11_sc_pseudobulk_de",
  run_name = cfg$run_name,
  params   = cfg,
  notes    = cfg$question %||% ""
)

on.exit({
  if (!isTRUE(getOption("asset.run_finalized"))) {
    finalize_run(run, status = "failed",
                 verdict = "script exited before finalize_run()")
  }
}, add = TRUE)

# ---- 2. load cells --------------------------------------------

message("\n== Loading cells ==")

# TODO(you): point this at the object your Layer 1 stages produce.
obj_path <- file.path(OBJ, sprintf("ASSET_%s_annotated.rds", FREEZE))
obj <- load_object(obj_path)

obj <- attach_labels(obj, freeze = FREEZE, labelset = LABELSET)

sel <- cfg$cell_selection
lab_col <- sel$label_column
if (!lab_col %in% colnames(obj[[]])) {
  stop("Label column '", lab_col, "' not found on the object.", call. = FALSE)
}

keep_cell <- obj[[lab_col]][, 1] %in% unlist(sel$labels)
if (!is.null(sel$timepoints)) {
  keep_cell <- keep_cell & obj$timepoint %in% unlist(sel$timepoints)
}

# cohort membership — requires clinical for post-escape censoring
clin <- NULL
ch <- read_cohort(COHORT)
if (isTRUE(ch$requires_clinical)) {
  clin_path <- file.path(CLINICAL, "asset_clinical.csv")   # never in git
  if (!file.exists(clin_path)) {
    stop("Cohort '", COHORT, "' needs the clinical table at ", clin_path,
         call. = FALSE)
  }
  clin <- utils::read.csv(clin_path, stringsAsFactors = FALSE)
}
cohort_libs <- cohort_samples(COHORT, FREEZE, clinical = clin)
keep_cell <- keep_cell & obj$library_id %in% cohort_libs$library_id

message(sprintf("  %d of %d cells selected (%.1f%%)",
                sum(keep_cell), length(keep_cell),
                100 * sum(keep_cell) / length(keep_cell)))
if (sum(keep_cell) < 500) {
  warning("Fewer than 500 cells selected. Check the label column and values.",
          call. = FALSE)
}

obj <- obj[, keep_cell]

# ---- 3. aggregate to pseudobulk -------------------------------

message("\n== Aggregating to pseudobulk ==")

agg  <- cfg$aggregation
unit <- agg$unit %||% "sample_id"

# Aggregate RAW COUNTS by summing. Do not average normalized values:
# the count-based models downstream assume counts, and averaging
# normalized expression discards the library-size information they need.
cells_by_unit <- split(colnames(obj), obj[[unit]][, 1])

counts_layer <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")

pb <- vapply(cells_by_unit, function(cc) {
  Matrix::rowSums(counts_layer[, cc, drop = FALSE])
}, FUN.VALUE = numeric(nrow(counts_layer)))
pb <- as.matrix(pb)
storage.mode(pb) <- "integer"

n_cells <- vapply(cells_by_unit, length, integer(1))

# per-unit metadata, taken from the first cell of each unit and checked
# for internal consistency
meta_cols <- intersect(
  c("sample_id", "subject_id", "library_id", "group", "arm", "timepoint",
    "batch_id", "run_id"),
  colnames(obj[[]]))

col_data <- do.call(rbind, lapply(names(cells_by_unit), function(u) {
  cc <- cells_by_unit[[u]]
  m  <- obj[[]][cc, meta_cols, drop = FALSE]
  n_unique <- vapply(m, function(x) length(unique(x)), integer(1))
  if (any(n_unique > 1)) {
    stop("Unit '", u, "' spans multiple values of: ",
         paste(names(n_unique)[n_unique > 1], collapse = ", "),
         ". Wrong aggregation unit?", call. = FALSE)
  }
  m[1, , drop = FALSE]
}))
rownames(col_data) <- names(cells_by_unit)
col_data$n_cells <- n_cells[rownames(col_data)]

# fraction of the sample's cells that are the target cell type — the
# "target cell proportion" QC your colleague plots, and a genuine
# confounder: composition shifts can masquerade as expression changes
tot <- table(obj[[unit]][, 1])
col_data$target_cell_prop <- as.numeric(col_data$n_cells / tot[rownames(col_data)])

message(sprintf("  %d pseudobulk units x %d genes",
                ncol(pb), nrow(pb)))

# drop low-cell units
min_cells <- agg$min_cells_per_unit %||% 20
drop_u <- col_data$n_cells < min_cells
if (any(drop_u)) {
  message(sprintf("  dropping %d units with < %d cells: %s",
                  sum(drop_u), min_cells,
                  paste(rownames(col_data)[drop_u], collapse = ", ")))
}
pb <- pb[, !drop_u, drop = FALSE]
col_data <- col_data[!drop_u, , drop = FALSE]

# join clinical variables needed by the design
if (!is.null(clin)) {
  j <- match(col_data$subject_id, clin$subject_id)
  for (v in setdiff(names(clin), "subject_id")) col_data[[v]] <- clin[[v]][j]
}

# gene filter
frac  <- agg$min_samples_expressing %||% 0.2
mincn <- agg$min_count %||% 10
keep_g <- rowSums(pb >= mincn) >= ceiling(frac * ncol(pb))
message(sprintf("  gene filter: %d of %d genes retained",
                sum(keep_g), length(keep_g)))
pb <- pb[keep_g, , drop = FALSE]

# ---- 4. donor sanity check ------------------------------------
# One row per subject: how many units, which timepoints, cell counts.
# The single most useful table for catching a mis-specified design
# before you fit anything.

donor_check <- do.call(rbind, lapply(
  split(col_data, col_data$subject_id), function(d) {
    data.frame(
      subject_id = d$subject_id[1],
      group      = d$group[1],
      arm        = if ("arm" %in% names(d)) d$arm[1] else NA,
      n_units    = nrow(d),
      timepoints = paste(sort(unique(d$timepoint)), collapse = ","),
      batches    = paste(sort(unique(d$batch_id)), collapse = ","),
      n_cells    = sum(d$n_cells),
      min_cells  = min(d$n_cells),
      stringsAsFactors = FALSE)
  }))

save_table(donor_check, run, "donor_sanity_check")
save_table(cbind(unit = rownames(col_data), col_data), run, "colData")
saveRDS(col_data, file.path(run$objects, "colData.rds"))
saveRDS(pb,       file.path(run$objects, "pseudo_filtered.rds"))
utils::write.csv(cbind(gene = rownames(pb), as.data.frame(pb)),
                 file.path(run$tables, "pseudobulk_counts.csv"),
                 row.names = FALSE)

# If any subject contributes more than one unit, a fixed-effects model
# is wrong — the units are not independent.
if (any(donor_check$n_units > 1) &&
    identical(cfg$model$engine, "deseq2") &&
    is.null(cfg$model$random_effect)) {
  warning("Subjects contribute multiple pseudobulk units, but engine is ",
          "deseq2 with no random effect. Repeated measures are being ",
          "treated as independent, which inflates significance. Use ",
          "engine: limma_voom (with duplicateCorrelation) or dream.",
          call. = FALSE)
}

# ---- 5. fit ---------------------------------------------------
# TODO(review): confirm the design and inspect the model matrix before
# trusting output. Check that every design term is (a) present in
# col_data, (b) not collinear with the contrast, (c) has enough levels.

message("\n== Fitting ==")
engine   <- cfg$model$engine
design_f <- stats::as.formula(cfg$model$design)
alpha    <- cfg$thresholds$alpha %||% 0.05

mm <- stats::model.matrix(design_f, data = col_data)
message("  model matrix: ", nrow(mm), " x ", ncol(mm),
        " (rank ", qr(mm)$rank, ")")
if (qr(mm)$rank < ncol(mm)) {
  stop("Design matrix is rank-deficient — terms are collinear. ",
       "Check tables/colData.csv for a covariate that is constant or ",
       "perfectly confounded with the contrast.", call. = FALSE)
}

res_df <- NULL

if (engine == "deseq2") {
  # Appropriate for CROSS-SECTIONAL contrasts only (one unit per subject).
  stopifnot(requireNamespace("DESeq2", quietly = TRUE))
  dds <- DESeq2::DESeqDataSetFromMatrix(pb, col_data, design_f)
  dds <- DESeq2::DESeq(dds)
  ct  <- cfg$model$contrast
  res <- DESeq2::results(
    dds, contrast = c(ct$variable, ct$numerator, ct$denominator),
    alpha = alpha)
  if (isTRUE(cfg$thresholds$lfc_shrink) &&
      requireNamespace("apeglm", quietly = TRUE)) {
    res <- DESeq2::lfcShrink(dds, res = res, type = "apeglm",
                             coef = DESeq2::resultsNames(dds)[
                               length(DESeq2::resultsNames(dds))])
  }
  vsd <- DESeq2::vst(dds, blind = TRUE)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)

  saveRDS(dds, file.path(run$objects, "dds.rds"))
  saveRDS(res, file.path(run$objects, "res.rds"))
  saveRDS(vsd, file.path(run$objects, "vsd.rds"))

} else if (engine == "limma_voom") {
  # Repeated measures via duplicateCorrelation on subject_id.
  stopifnot(requireNamespace("limma", quietly = TRUE),
            requireNamespace("edgeR", quietly = TRUE))
  d <- edgeR::DGEList(pb); d <- edgeR::calcNormFactors(d)
  v <- limma::voom(d, mm, plot = FALSE)
  dc <- limma::duplicateCorrelation(v, mm, block = col_data$subject_id)
  message("  consensus within-subject correlation: ",
          round(dc$consensus.correlation, 3))
  fit <- limma::lmFit(v, mm, block = col_data$subject_id,
                      correlation = dc$consensus.correlation)
  fit <- limma::eBayes(fit)
  coef_name <- tail(colnames(mm), 1)   # TODO(review): name the coefficient explicitly
  res_df <- limma::topTable(fit, coef = coef_name, number = Inf,
                            sort.by = "P")
  res_df$gene <- rownames(res_df)
  saveRDS(fit, file.path(run$objects, "dds.rds"))
  saveRDS(v,   file.path(run$objects, "vsd.rds"))

} else if (engine == "dream") {
  # Proper mixed model. Preferred for longitudinal contrasts.
  stopifnot(requireNamespace("variancePartition", quietly = TRUE),
            requireNamespace("edgeR", quietly = TRUE))
  d <- edgeR::DGEList(pb); d <- edgeR::calcNormFactors(d)
  form <- stats::as.formula(paste(cfg$model$design,
                                  cfg$model$random_effect, sep = " + "))
  vobj <- variancePartition::voomWithDreamWeights(d, form, col_data)
  fit  <- variancePartition::dream(vobj, form, col_data)
  fit  <- variancePartition::eBayes(fit)
  coef_name <- tail(colnames(fit$coefficients), 1)  # TODO(review)
  res_df <- variancePartition::topTable(fit, coef = coef_name, number = Inf)
  res_df$gene <- rownames(res_df)
  saveRDS(fit,  file.path(run$objects, "dds.rds"))
  saveRDS(vobj, file.path(run$objects, "vsd.rds"))

} else {
  stop("Unknown engine: ", engine,
       " (expected deseq2, limma_voom, or dream)", call. = FALSE)
}

saveRDS(res_df, file.path(run$objects, "res_df.rds"))
save_table(res_df, run, "DESeq2_all_results")

# ---- 6. annotate and summarise --------------------------------

p_col  <- intersect(c("padj", "adj.P.Val"), names(res_df))[1]
lfc_col <- intersect(c("log2FoldChange", "logFC"), names(res_df))[1]
lfc_min <- cfg$thresholds$report_lfc_min %||% 0

sig <- res_df[!is.na(res_df[[p_col]]) &
              res_df[[p_col]] < alpha &
              abs(res_df[[lfc_col]]) >= lfc_min, ]
sig <- sig[order(sig[[p_col]]), ]
save_table(sig, run, "DESeq2_results_annotated")

message(sprintf("\n  %d genes at %s < %.3f and |%s| >= %.2f",
                nrow(sig), p_col, alpha, lfc_col, lfc_min))

# ---- 7. plots -------------------------------------------------
# TODO(you): fill these in using R/theme_asset.R. Keep the filenames
# exactly as they are so this run directory is interchangeable with
# your colleague's.
#
#   plots/volcano.pdf
#   plots/heatmap_top_up_genes.pdf
#   plots/top_gene_norm_counts.pdf
#   plots/QC/QC_PCA_plots.pdf
#   plots/QC/QC_dispersion.pdf
#   plots/QC/QC_sample_correlation.pdf
#   plots/QC/QC_cells_per_donor.pdf
#   plots/QC/QC_libsize_vs_genes.pdf
#   plots/QC/QC_pseudobulk_library_sizes.pdf
#   plots/QC/QC_genes_detected_per_sample.pdf
#   plots/QC/QC_target_cell_proportion.pdf
#   plots/QC/QC_MA_plot.pdf
#   plots/QC/QC_top500_variable_genes_heatmap.pdf
#
# Two of these carry real interpretive weight for this project:
#   QC_PCA_plots            colour by batch_id AND by timepoint. If PC1-2
#                           separate by batch, stop and fix that first.
#   QC_target_cell_proportion  vs the contrast variable. If improvers have
#                           systematically more of the target cell type,
#                           a composition shift can masquerade as an
#                           expression change.

message("\n  TODO: plotting block not yet implemented")

# ---- 8. close -------------------------------------------------

options(asset.run_finalized = TRUE)
finalize_run(
  run,
  status  = "ok",
  verdict = sprintf("%d DE genes (%s<%.2f, |%s|>=%.2f) in %s; engine=%s",
                    nrow(sig), p_col, alpha, lfc_col, lfc_min,
                    paste(unlist(cfg$cell_selection$labels), collapse = "/"),
                    engine)
)
