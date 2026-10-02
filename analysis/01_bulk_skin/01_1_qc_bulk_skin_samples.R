# ==============================================================
# 01_1_qc_bulk_skin_samples.R
#
# Per-sample QC of the 234 bulk skin samples, with the SAME metrics and
# rules as the owner's PBMC QC (analysis/02_bulk_pbmc/1.4_QC_70_samples.R,
# gitignored), and the final "use downstream" column. Results are added to
# the skin workbook in place, after a backup (decisions.md 2026-10-02).
#
# Metrics (as blood): raw_library_size (column sums of the Ensembl raw
# counts), genes_detected_TPM1, hemoglobin_TPM_pct (HBA1/HBA2/HBB/HBD),
# mito_TPM_pct (MT-*), median_cor_to_others (Pearson, log2(TPM+1), expressed
# genes = TPM >= 1 in >= 20% of samples), PC1, PC2, PCA_distance_PC1to5
# (PCs scaled by their sd). corrected_library_size and
# cor_before_vs_after_combat are blood-only (ComBat-seq); left blank.
# Rules (as blood): Low depth < 3M raw counts; median correlation robust
# z < -3; PCA distance robust z > 3; genes detected robust z < -3 or > 3.
# QC_flag: "exclude_suggested" if a PCA outlier, "check" for the others,
# else "ok". Hemoglobin and mito are reported, not flagged (as blood).
#
# Use_downstream = "Yes" if Exclude and Repeat are blank AND QC_flag is not
# "exclude_suggested"; otherwise "No", with Use_reason.
#
# Matrix column fixes from 00_7 (data/clinical/bulk_skin_matrix_column_fixes.csv)
# are applied when reading the matrices.
# Counts and plots to results/; per-sample values to data/patient_level/.
#
#   jobs/run.sh analysis/01_bulk_skin/01_1_qc_bulk_skin_samples.R --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")
source("R/theme_asset.R")

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(openxlsx)
  library(ggplot2)
})

P <- list(skin_file = "ASSET_clinical_data_skin_bulk_metadata.xlsx",
          counts_file = "skin/raw_counts_mtx_ASSET_skin_EMSEMBL_64252_234.txt",
          tpm_file = "skin/TPM_exp_mtx_ASSET_skin_geneSymbol_39921_234.txt",
          expr_tpm = 1, expr_frac = 0.2, depth_cut = 3e6, rz_cut = 3, n_pcs = 5,
          hb_genes = c("HBA1", "HBA2", "HBB", "HBD"))

run <- init_run(
  stage    = "01_bulk_skin",
  run_name = "qc_bulk_skin_samples",
  params   = P,
  notes    = "Bulk skin per-sample QC with the PBMC metrics; Use_downstream column; written into the skin workbook"
)

# ---- inputs ------------------------------------------------------
skin <- as.data.table(read_excel(file.path(CLINICAL, P$skin_file), col_types = "text"))
fix_file <- file.path(CLINICAL, "bulk_skin_matrix_column_fixes.csv")
fixes <- if (file.exists(fix_file)) fread(fix_file) else data.table(matrix_column = character(), correct_Sample_ID = character())
read_mat <- function(f) {
  m <- fread(file.path(BULK, f))
  ids <- names(m)[-1]
  ids[ids %in% fixes$matrix_column] <- fixes$correct_Sample_ID[match(ids[ids %in% fixes$matrix_column], fixes$matrix_column)]
  keep <- ids %in% skin$Sample_ID                                  # drops annotation columns
  x <- as.matrix(m[, -1][, ..keep]); colnames(x) <- ids[keep]; rownames(x) <- m[[1]]
  x[, skin$Sample_ID]
}
cnt <- read_mat(P$counts_file)
tpm <- read_mat(P$tpm_file)
stopifnot(identical(colnames(cnt), skin$Sample_ID), identical(colnames(tpm), skin$Sample_ID))

# ---- metrics (definitions as the PBMC QC) --------------------------
rz <- function(x) (x - median(x)) / mad(x)
logt <- log2(tpm + 1)
expr <- rowSums(tpm >= P$expr_tpm) >= P$expr_frac * ncol(tpm)
pca <- prcomp(t(logt[expr, ]), center = TRUE, scale. = FALSE)
cc <- cor(logt[expr, ])
hb <- intersect(P$hb_genes, rownames(tpm))
qc <- data.table(
  Sample_ID = skin$Sample_ID,
  raw_library_size = colSums(cnt),
  corrected_library_size = NA_real_,
  genes_detected_TPM1 = colSums(tpm >= 1),
  hemoglobin_TPM_pct = round(colSums(tpm[hb, , drop = FALSE]) / 1e4, 2),
  mito_TPM_pct = round(colSums(tpm[grep("^MT-", rownames(tpm)), , drop = FALSE]) / 1e4, 2),
  median_cor_to_others = round(apply(cc, 2, function(x) median(x[x < 1])), 3),
  cor_before_vs_after_combat = NA_real_,
  PC1 = round(pca$x[, 1], 2), PC2 = round(pca$x[, 2], 2),
  PCA_distance_PC1to5 = round(sqrt(rowSums(scale(pca$x[, 1:P$n_pcs], center = TRUE,
                                                  scale = pca$sdev[1:P$n_pcs])^2)), 2))
z_cor <- rz(qc$median_cor_to_others); z_pca <- rz(qc$PCA_distance_PC1to5); z_gen <- rz(qc$genes_detected_TPM1)
note <- vapply(seq_len(nrow(qc)), function(i) {
  n <- c(if (qc$raw_library_size[i] < P$depth_cut) sprintf("Low depth: %.1fM raw counts (<3M)", qc$raw_library_size[i] / 1e6),
         if (z_cor[i] < -P$rz_cut) sprintf("Low similarity to other samples (median r=%.3f; robust z<-3)", qc$median_cor_to_others[i]),
         if (z_pca[i] > P$rz_cut) "PCA outlier (PC1-5 distance robust z>3)",
         if (z_gen[i] < -P$rz_cut) sprintf("Few genes detected (%d with TPM>=1; robust z<-3)", qc$genes_detected_TPM1[i]),
         if (z_gen[i] > P$rz_cut) sprintf("Many genes detected (%d with TPM>=1; robust z>3)", qc$genes_detected_TPM1[i]))
  if (is.null(n)) NA_character_ else paste(n, collapse = "; ")
}, "")
qc[, Note := note]
qc[, QC_flag := fifelse(grepl("PCA outlier", Note), "exclude_suggested",
                fifelse(!is.na(Note), "check", "ok"))]

# ---- final column ---------------------------------------------------
qc[, `:=`(Exclude = skin$Exclude, Repeat = skin$Repeat, Subject_ID = skin$Subject_ID, Timepoint = skin$Timepoint)]
qc[, Use_reason := trimws(paste(fifelse(!is.na(Exclude), "Exclude (owner)", ""),
                                fifelse(!is.na(Repeat), "Repeat", ""),
                                fifelse(QC_flag == "exclude_suggested", "QC: PCA outlier", "")))]
qc[, Use_reason := gsub("\\s+", "; ", Use_reason)]
qc[Use_reason == "", Use_reason := NA_character_]
qc[, Use_downstream := fifelse(is.na(Use_reason), "Yes", "No")]
stopifnot(qc[Use_downstream == "Yes", !anyDuplicated(paste(Subject_ID, Timepoint))])

# ---- comparisons (counts only) -------------------------------------
cmp <- qc[, .N, by = .(Exclude = fifelse(is.na(Exclude), "kept", "Exclude"),
                      Repeat = fifelse(is.na(Repeat), "kept", "Repeat"), QC_flag)][order(Exclude, Repeat, QC_flag)]
save_table(cmp, run, "qc_flag_vs_owner_flags")
# units where the owner's kept sample is a QC outlier but a Repeat-flagged twin passes
tw <- qc[, .(n = .N, kept_bad = any(is.na(Exclude) & is.na(Repeat) & QC_flag == "exclude_suggested"),
             repeat_ok = any(!is.na(Repeat) & is.na(Exclude) & QC_flag != "exclude_suggested")), by = .(Subject_ID, Timepoint)]
units <- qc[Use_downstream == "Yes", .N, by = Timepoint]
lost <- qc[is.na(Exclude) & is.na(Repeat) & QC_flag == "exclude_suggested", .N]
batch_r2 <- if ("ASSETpaper_2022" %in% names(skin)) {
  g <- factor(is.na(skin$ASSETpaper_2022))
  sapply(1:5, function(k) summary(lm(pca$x[, k] ~ g))$r.squared)
} else rep(NA_real_, 5)
summ <- data.table(item = c("samples", "expressed genes used", "QC ok", "QC check", "QC exclude_suggested",
                            "owner-kept samples that QC would drop (PCA outlier)",
                            "units where kept sample fails QC but its Repeat twin passes",
                            "Use_downstream = Yes", paste("Use_downstream = Yes at", units$Timepoint),
                            "subjects with >= 1 usable sample", "subjects with a usable Baseline sample",
                            sprintf("R2 of PC%d with the 2022-paper set", 1:5)),
                   n = c(nrow(qc), sum(expr), sum(qc$QC_flag == "ok"), sum(qc$QC_flag == "check"),
                         sum(qc$QC_flag == "exclude_suggested"), lost, sum(tw$kept_bad & tw$repeat_ok),
                         sum(qc$Use_downstream == "Yes"), units$N,
                         qc[Use_downstream == "Yes", uniqueN(Subject_ID)],
                         qc[Use_downstream == "Yes" & Timepoint == "Baseline", uniqueN(Subject_ID)],
                         round(batch_r2, 3)))
save_table(summ, run, "qc_summary")
note_kinds <- qc[!is.na(Note), .(rule = unlist(strsplit(Note, "; ")))][, .(rule = sub(":.*| \\(.*", "", rule))][, .N, by = rule]
save_table(note_kinds, run, "qc_rules_triggered")
save_patient_table(qc, run, "skin_sample_qc")

p <- ggplot(qc, aes(PC1, PC2, colour = QC_flag, shape = Use_downstream)) + geom_point(size = 1.6) +
  labs(subtitle = sprintf("Bulk skin, log2(TPM+1), %d expressed genes", sum(expr))) + theme_asset()
save_plot(p, run, "pca_qc_flag", width = 6, height = 4.5)
if (!all(is.na(batch_r2))) {
  p2 <- ggplot(cbind(qc, set = ifelse(is.na(skin$ASSETpaper_2022), "not in 2022 paper", "2022 paper")),
               aes(PC1, PC2, colour = set)) + geom_point(size = 1.4) + theme_asset()
  save_plot(p2, run, "pca_by_2022_paper_set", width = 6, height = 4.5)
}

# ---- write into the skin workbook (backup, temp file, verify) -------
bk <- file.path(CLINICAL, "backup", run$run_id)
dir.create(bk, recursive = TRUE, showWarnings = FALSE, mode = "0700")
stopifnot(file.copy(file.path(CLINICAL, P$skin_file), file.path(bk, P$skin_file), copy.date = TRUE))
Sys.chmod(file.path(bk, P$skin_file), "0400")
new_cols <- c("raw_library_size", "corrected_library_size", "genes_detected_TPM1", "hemoglobin_TPM_pct",
              "mito_TPM_pct", "median_cor_to_others", "cor_before_vs_after_combat", "PC1", "PC2",
              "PCA_distance_PC1to5", "Note", "QC_flag", "Use_downstream", "Use_reason")
path <- file.path(CLINICAL, P$skin_file)
wb <- loadWorkbook(path); sh <- names(wb)[1]
hdr <- names(skin)
for (c in new_cols) {
  j <- match(c, hdr); if (is.na(j)) { hdr <- c(hdr, c); j <- length(hdr) }
  writeData(wb, sh, x = qc[[c]], startCol = j, startRow = 2, colNames = FALSE)   # NA -> blank cell
  writeData(wb, sh, x = c, startCol = j, startRow = 1)
}
tmp <- file.path(CLINICAL, paste0(".tmp_", run$run_id, "_", P$skin_file))
saveWorkbook(wb, tmp, overwrite = TRUE)
y <- as.data.table(read_excel(tmp, col_types = "text"))
stopifnot(identical(names(y), hdr), nrow(y) == nrow(skin))
for (c in names(skin)[!names(skin) %in% new_cols]) if (!identical(y[[c]], skin[[c]])) stop("column ", c, " changed unexpectedly")
stopifnot(identical(y$Use_downstream, qc$Use_downstream), identical(y$QC_flag, qc$QC_flag))
stopifnot(file.rename(tmp, path))

writeLines(c(sprintf("%-62s %s", summ$item, summ$n), "", "Owner flags x QC_flag:",
             sprintf("  Exclude=%-8s Repeat=%-7s QC=%-18s %d", cmp$Exclude, cmp$Repeat, cmp$QC_flag, cmp$N),
             "", "Rules triggered:", sprintf("  %-45s %d", note_kinds$rule, note_kinds$N),
             "", sprintf("Workbook updated (backup %s); %d columns added/updated.", .rel_to_repo(bk), length(new_cols))),
           file.path(run$dir, "SKIN_QC_SUMMARY.txt"))
finalize_run(run, status = "ok", verdict = sprintf("%d ok / %d check / %d exclude_suggested; Use_downstream Yes = %d (%d subjects)",
  sum(qc$QC_flag == "ok"), sum(qc$QC_flag == "check"), sum(qc$QC_flag == "exclude_suggested"),
  sum(qc$Use_downstream == "Yes"), qc[Use_downstream == "Yes", uniqueN(Subject_ID)]))
