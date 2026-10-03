# ==============================================================
# 06_0_match_upstream_cells_and_labels.R
#
# Match our cells to the lab's upstream objects (wasikowr, jarnagin) by
# barcode, and bring in every cell label that still exists. Three jobs:
#
#   1. Cell-set comparison: of wasikowr's 1.30M cells, how many do we
#      call and keep, and why do the others differ (our QC reason)?
#   2. Inherited labels: wasikowr's 16 broad cell types and jarnagin's
#      compartment labels, joined to our cell IDs. Written to
#      metadata/inherited_labels/ (barcode-level, gitignored).
#   3. Quality of wasikowr's broad labels: for each label, the fraction of
#      cells detecting each lab marker panel (a good label is high on its
#      own panel and low on the others), median UMI (doublet-like if well
#      above the rest), and agreement with jarnagin's later refinements.
#
# Matching key: library + 24-nt barcode. Upstream cell names carry other
# suffixes ("-3_1"); orig.ident is the sample-map `condition`, which our
# manifest keeps as source_orig_ident.
#
# Read-only on the upstream files. Provenance is recorded in
# docs/inherited_objects.md. Nothing here looks at mRSS_category.
#
#   jobs/run.sh analysis/06_sc_lineage/06_0_match_upstream_cells_and_labels.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")
source("R/theme_asset.R")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(Seurat)
  library(BPCells)
  library(SingleCellExperiment)
})

JF <- "/home/jarnagin/ASSET_Flex/objects"
SOURCES <- list(
  wasikowr_object = "/hits/home/wasikowr/ffpe/Asset/seurat.RDS",
  wasikowr_basectrl = file.path(JF, "ASSET_base_control.RDS"),       # older wasikowr version, celltype present
  music_reference = file.path(JF, "ASSET_MuSiC_reference_sce.rds"),
  jf_immune       = file.path(JF, "ASSET_Round2_immFilt.RDS"),
  jf_tcells       = file.path(JF, "ASSET_FullDataset_Tcells_Subtyped.rds"),
  jf_bplasma      = file.path(JF, "ASSET_FullDataset_BcellPlasma_Clean.rds"),
  jf_fibroblast   = file.path(JF, "ASSET_Round1_fibs.RDS")
)
LABEL_RE <- "^(celltype|subclusters|compartment|seurat_clusters|RNA_snn_res\\..*|.*subtype.*|.*lineage.*|.*label.*|.*_call)$"
NOT_LABEL <- c("subject", "patient", "treatment", "time", "antibody", "batch", "temp")

run <- init_run(
  stage    = "06_sc_lineage",
  run_name = "match_upstream_cells_and_labels",
  params   = list(sources = SOURCES, label_regex = LABEL_RE),
  notes    = "Barcode match to wasikowr/jarnagin objects; inherited labels; checks on wasikowr's broad labels"
)

man <- read_manifest()
fz  <- freeze_libraries(manifest = man)
stopifnot(!anyDuplicated(stats::na.omit(man$source_orig_ident)))   # lib_of below must be one-to-one
n_excluded <- length(setdiff(man$library_id, fz$library_id))

# Cell name -> 24-nt barcode. Assumes suffix-style names ("ACGT...-3_1");
# a prefix-style name ("S12_ACGT...") would reduce to the prefix, so fail
# loudly rather than let every key collide and be dropped as a duplicate.
core <- function(x, what) {
  b <- sub("[-_].*$", "", x)
  bad <- !grepl("^[ACGT]{24}$", b)
  if (any(bad)) stop(sprintf("%s: %s of %s cell names do not reduce to a 24-nt barcode (e.g. '%s')",
                             what, format(sum(bad), big.mark = ","), format(length(b), big.mark = ","),
                             x[which(bad)[1]]), call. = FALSE)
  b
}

qc <- fread(file.path(results_dir("03_sc_qc", "cell_qc_per_library", cohort = "reference"),
                      "objects", "cell_qc.csv.gz"),
            select = c("cell_id", "library_id", "barcode", "n_umi", "n_genes", "qc_pass", "qc_fail_reason"))
qc <- qc[library_id %in% fz$library_id]
qc[, key := paste(library_id, core(barcode, "our cell_qc"), sep = "|")]
stopifnot(!anyDuplicated(qc$key))

# ---- read metadata from each upstream object ------------------

read_meta <- function(path, keep_object = FALSE) {
  o  <- readRDS(path)
  md <- if (inherits(o, "SingleCellExperiment")) as.data.frame(colData(o)) else o@meta.data
  md <- as.data.table(md, keep.rownames = "upstream_cell")
  lab <- setdiff(grep(LABEL_RE, names(md), value = TRUE), NOT_LABEL)
  lab <- lab[!vapply(md[, ..lab], is.numeric, logical(1))]
  keep <- intersect(c("upstream_cell", "orig.ident", "nCount_RNA", "nFeature_RNA", "time", lab), names(md))
  if (!keep_object) { rm(o); o <- NULL; invisible(gc()) }
  list(meta = md[, ..keep], labels = lab, object = o)
}

lib_of <- setNames(man$library_id, man$source_orig_ident)
meta <- list(); lab_cols <- list(); src_summary <- list(); bc <- NULL
for (s in names(SOURCES)) {
  if (!file.exists(SOURCES[[s]])) { message("  missing, skipped: ", SOURCES[[s]]); next }
  message("Reading ", s)
  r <- read_meta(SOURCES[[s]], keep_object = (s == "wasikowr_basectrl"))   # section 3 needs its counts
  if (s == "wasikowr_basectrl") bc <- r$object
  r$object <- NULL
  m <- r$meta
  m[, library_id := unname(lib_of[as.character(orig.ident)])]
  m[!is.na(library_id), key := paste(library_id, core(upstream_cell, s), sep = "|")]   # unmapped libraries: key NA
  dup <- !is.na(m$key) & (duplicated(m$key) | duplicated(m$key, fromLast = TRUE))
  src_summary[[s]] <- data.table(
    source = s, path = SOURCES[[s]], cells = nrow(m),
    orig_ident_not_in_manifest = sum(is.na(m$library_id)),
    in_excluded_libraries = sum(!is.na(m$library_id) & !m$library_id %in% fz$library_id),
    duplicated_keys = sum(dup),
    matched_our_called = sum(m$key %in% qc$key),
    matched_our_pass = sum(m$key %in% qc$key[qc$qc_pass]),
    label_columns = paste(r$labels, collapse = ";"))
  meta[[s]] <- m[!dup & !is.na(library_id)]
  lab_cols[[s]] <- r$labels
}
src_summary <- rbindlist(src_summary)
save_table(src_summary, run, "upstream_sources")

# ---- 1. cell-set comparison against wasikowr's object ---------

w <- meta$wasikowr_object
qc[, in_wasikowr := key %in% w$key]
per_lib <- qc[, .(our_called = .N, our_pass = sum(qc_pass),
                  wasikowr = sum(in_wasikowr),
                  both_pass = sum(in_wasikowr & qc_pass),
                  ours_only_pass = sum(!in_wasikowr & qc_pass),
                  hers_failed_ours = sum(in_wasikowr & !qc_pass)), by = library_id]
hers_not_called <- w[library_id %in% fz$library_id & !key %in% qc$key, .N, by = library_id]
per_lib <- merge(per_lib, hers_not_called[, .(library_id, hers_not_called_by_us = N)],
                 by = "library_id", all.x = TRUE)
per_lib[is.na(hers_not_called_by_us), hers_not_called_by_us := 0L]
per_lib <- merge(as.data.table(fz[, c("library_id", "batch_id", "group", "arm", "timepoint")]),
                 per_lib, by = "library_id")
save_table(per_lib, run, "cellset_overlap_per_library")

why <- qc[in_wasikowr & !qc_pass, .N, by = qc_fail_reason][order(-N)]
save_table(why, run, "wasikowr_cells_failing_our_qc_by_reason")
ours_only <- qc[!in_wasikowr & qc_pass, .(cells = .N, median_umi = stats::median(n_umi),
                                          median_genes = stats::median(n_genes))]
save_table(ours_only, run, "our_cells_absent_from_wasikowr")

# ---- 2. inherited label table ---------------------------------

lab <- qc[, .(cell_id, key, qc_pass)]
for (s in names(meta)) {
  cols <- lab_cols[[s]]
  if (!length(cols)) next
  m <- meta[[s]][, c("key", cols), with = FALSE]
  setnames(m, cols, paste(s, cols, sep = "__"))
  lab <- merge(lab, m, by = "key", all.x = TRUE)
}
lab[, key := NULL]
dir.create(file.path(METADATA, "inherited_labels"), showWarnings = FALSE)
lab_file <- file.path(METADATA, "inherited_labels",
                      sprintf("%s_upstream_labels.csv.gz", FREEZE))
fwrite(lab, lab_file)
Sys.chmod(lab_file, "0600")
cat(.rel_to_repo(lab_file), "\n", file = file.path(run$tables, "BARCODE_LEVEL_TABLES.txt"))
message("  wrote ", .rel_to_repo(lab_file), " (barcode-level)")

# coverage per label column (aggregate only)
cov <- rbindlist(lapply(setdiff(names(lab), c("cell_id", "qc_pass")), function(cn)
  data.table(column = cn, cells_labelled = sum(!is.na(lab[[cn]])),
             of_our_pass = sum(!is.na(lab[[cn]]) & lab$qc_pass),
             n_levels = length(unique(stats::na.omit(lab[[cn]]))))))
save_table(cov, run, "inherited_label_coverage")

# broad label by timepoint (from whichever source carries `celltype`)
bt <- melt(lab[, c("cell_id", grep("__celltype$", names(lab), value = TRUE)), with = FALSE],
           id.vars = "cell_id", variable.name = "source", value.name = "celltype", na.rm = TRUE)
bt[, library_id := qc$library_id[match(cell_id, qc$cell_id)]]
bt[, timepoint := fz$timepoint[match(library_id, fz$library_id)]]
bt[, group := fz$group[match(library_id, fz$library_id)]]
save_table(dcast(bt, source + celltype ~ group + timepoint, fun.aggregate = length,
                 value.var = "cell_id"), run, "inherited_celltype_by_timepoint")

# ---- 3. quality of wasikowr's broad labels --------------------
# On her own labelled object (baseline + control): detection of each lab
# panel per label, and median UMI per label.

mk <- yaml::read_yaml(file.path(GENESETS, "lab_markers_jarnagin.yml"))
panels <- c(mk$evidence_panels,
            list(Sebocyte = c("MGST1", "FADS2", "AWAT2", "DGAT2"),       # appendage panels added here;
                 Eccrine  = c("DCD", "SCGB2A2", "MUCL1", "KRT19"),        # the lab panels have none
                 Follicle = c("KRT75", "SOX9", "KRT6B", "KRT17"),
                 Pericyte_SMC = c("RGS5", "ACTA2", "MYH11", "TAGLN"),
                 Lymphatic = c("PROX1", "LYVE1", "CCL21", "TFF3"),
                 Adipocyte = c("ADIPOQ", "PLIN1", "FABP4", "PLIN4")))
if (is.null(bc)) stop("wasikowr_basectrl object not read: ", SOURCES$wasikowr_basectrl, call. = FALSE)
genes <- intersect(unique(unlist(panels)), rownames(bc))
cnt <- as(LayerData(bc, assay = "RNA", layer = "counts")[genes, ], "dgCMatrix")
ct  <- as.character(bc$celltype)
det <- rbindlist(lapply(names(panels), function(p) {
  g <- intersect(panels[[p]], genes)
  if (!length(g)) return(NULL)
  pos <- Matrix::colSums(cnt[g, , drop = FALSE] > 0) > 0
  data.table(panel = p, celltype = names(tapply(pos, ct, mean)), frac_cells = as.numeric(tapply(pos, ct, mean)))
}))
umi <- data.table(celltype = ct, n_umi = bc$nCount_RNA)[, .(cells = .N, median_umi = stats::median(n_umi)), by = celltype]
umi[, umi_vs_median_label := round(median_umi / stats::median(median_umi), 2)]
rm(bc, cnt); invisible(gc())
save_table(dcast(det, celltype ~ panel, value.var = "frac_cells"), run, "wasikowr_label_panel_detection")
save_table(umi[order(-cells)], run, "wasikowr_label_umi")

p_det <- ggplot(det, aes(panel, celltype, fill = frac_cells)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", frac_cells)), size = 1.6) +
  scale_fill_gradient(low = "white", high = "firebrick3", limits = c(0, 1)) +
  labs(x = "Marker panel (lab evidence panels + appendage panels)", y = "wasikowr celltype",
       fill = "Fraction of\ncells detecting",
       subtitle = "A good label is high on its own panel and low elsewhere (baseline + control cells)") +
  theme_asset() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_det, run, "wasikowr_label_panel_detection", width = 8, height = 5)

# agreement: wasikowr broad label vs jarnagin's refinements (counts only)
agree <- list()
for (s in c("jf_immune", "jf_tcells", "jf_bplasma", "jf_fibroblast")) {
  if (is.null(meta[[s]]) || !"celltype" %in% names(meta[[s]])) next
  ref <- grep("^(subclusters|lineage_fine_clean|.*subtype.*)$", lab_cols[[s]], value = TRUE)[1]
  if (is.na(ref)) next
  a <- meta[[s]][, .N, by = c("celltype", ref)]
  setnames(a, c("wasikowr_celltype", "jarnagin_label", "cells"))
  agree[[s]] <- cbind(source = s, jarnagin_column = ref, a)
}
if (length(agree)) save_table(rbindlist(agree)[order(source, -cells)], run,
                              "wasikowr_vs_jarnagin_labels")

# against our own per-library coarse lineage (SoupX stage). Required: R0021
# ran while SoupX was still running and skipped this table without notice.
sx <- file.path(results_dir("04_sc_ambient", "soupx_per_library", cohort = "reference"),
                "objects", "soupx_cells.csv.gz")
if (!file.exists(sx)) stop("SoupX cell table missing (run 04_1 first): ", sx, call. = FALSE)
s4 <- fread(sx, select = c("cell_id", "coarse_lineage"))
wc <- grep("__celltype$", names(lab), value = TRUE)[1]
x <- merge(lab[, c("cell_id", wc), with = FALSE], s4, by = "cell_id")
setnames(x, wc, "wasikowr_celltype")
save_table(x[!is.na(wasikowr_celltype), .N, by = .(wasikowr_celltype, coarse_lineage)][order(wasikowr_celltype, -N)],
           run, "wasikowr_vs_our_coarse_lineage")

# ---- summary ------------------------------------------------

wo <- src_summary[source == "wasikowr_object"]
own_panel <- c(Keratinocytes = "Keratinocyte", Fibroblasts = "Fibroblast", `T Cells` = "panT",
               `Myeloid Cells` = "Myeloid", `Endothelial Cells` = "Endothelial",
               Melanocytes = "Melanocyte", `Mast Cells` = "Mast", `B Cells` = "B_Plasma",
               `Langerhans Cells` = "DC_LC", Sebocytes = "Sebocyte", `Eccrine Cells` = "Eccrine",
               `Follicle Cells` = "Follicle", Pericytes = "Pericyte_SMC",
               `Smooth Muscle Cells` = "Pericyte_SMC", `L Endothelial Cells` = "Lymphatic",
               Adipocytes = "Adipocyte")
unmapped <- setdiff(unique(det$celltype), names(own_panel))
if (length(unmapped)) stop("wasikowr celltype(s) with no own panel: ", paste(unmapped, collapse = ", "), call. = FALSE)
stopifnot(all(own_panel %in% det$panel))
own <- det[, .(own = frac_cells[panel == own_panel[celltype[1]]],
               best_other = max(frac_cells[panel != own_panel[celltype[1]]]),
               best_other_panel = panel[panel != own_panel[celltype[1]]][which.max(frac_cells[panel != own_panel[celltype[1]]])]),
           by = celltype]
summary_lines <- c(
  sprintf("wasikowr object: %s cells; matched to our called cells: %s; to our QC-pass cells: %s",
          format(wo$cells, big.mark = ","), format(wo$matched_our_called, big.mark = ","),
          format(wo$matched_our_pass, big.mark = ",")),
  sprintf("  her cells from libraries not in our sample map: %s; from our %d excluded libraries: %s",
          format(wo$orig_ident_not_in_manifest, big.mark = ","), n_excluded,
          format(wo$in_excluded_libraries, big.mark = ",")),
  sprintf("  our QC-pass cells absent from her object: %s (median %s UMI)",
          format(ours_only$cells, big.mark = ","), ours_only$median_umi),
  "  her cells failing our QC, by reason:",
  paste0("    ", why$qc_fail_reason, ": ", format(why$N, big.mark = ",")),
  "",
  "Inherited labels (cells labelled / of our QC-pass cells):",
  sprintf("  %-55s %9s / %9s", cov$column, format(cov$cells_labelled, big.mark = ","),
          format(cov$of_our_pass, big.mark = ",")),
  "",
  "wasikowr's broad labels: fraction detecting own panel vs best other panel:",
  sprintf("  %-22s own %.2f   best other %.2f (%s)", own$celltype, own$own, own$best_other, own$best_other_panel))
writeLines(summary_lines, file.path(run$dir, "UPSTREAM_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

finalize_run(run, status = "ok", verdict = sprintf(
  "%s of %s wasikowr cells match our QC-pass cells; %s of ours absent from hers; broad labels exist for %s of our cells",
  format(wo$matched_our_pass, big.mark = ","), format(wo$cells, big.mark = ","),
  format(ours_only$cells, big.mark = ","),
  format(max(cov$cells_labelled[grepl("__celltype$", cov$column)], 0), big.mark = ",")))
