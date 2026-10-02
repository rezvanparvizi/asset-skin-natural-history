# ==============================================================
# 00_6_resolve_bulk_metadata_mismatches.R
#
# Follow-up to 00_5 (R0030). 00_5 compared raw values, so pure coding
# differences counted as mismatches: Sex "F"/"M" (bulk) vs "Female"/"Male"
# (master); "Non_Improver" (bulk) vs "Non-Improver" (master);
# "TRUE"/"FALSE" vs "True"/"False". This run harmonises those codings
# and reports only the disagreements that remain, by kind:
#   bulk missing / master missing / both present and different.
# It also measures how far the one skin Sample_ID missing from the
# matrices is from the unexplained extra matrix column (edit distance).
#
# Scripts-only access (decisions.md 2026-10-02): counts to tables/,
# ID-level rows to data/patient_level/. Read-only on data/.
#
#   jobs/run.sh analysis/00_manifest_and_batch/00_6_resolve_bulk_metadata_mismatches.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

P <- list(skin_file = "ASSET_clinical_data_skin_bulk_metadata.xlsx",
          pbmc_file = "ASSET_clinical_data_PBMC_bulk_metadata.xlsx",
          coding = list(Sex = c(F = "Female", M = "Male"),
                        any = c(Non_Improver = "Non-Improver", TRUE. = "True", FALSE. = "False")))

run <- init_run(
  stage    = "00_manifest_and_batch",
  run_name = "resolve_bulk_metadata_mismatches",
  params   = P,
  notes    = "Field disagreements bulk metadata vs master after harmonising coding; residual Sample_ID mismatch"
)

rd <- function(f) {
  x <- as.data.table(read_excel(file.path(CLINICAL, f), col_types = "text"))
  x[, names(x) := lapply(.SD, function(v) { v <- trimws(v); v[!is.na(v) & (!nzchar(v) | toupper(v) == "NA")] <- NA; v })]
  x
}
master <- rd(basename(CLINICAL_MASTER)); skin <- rd(P$skin_file); pbmc <- rd(P$pbmc_file)

harmonise <- function(v, field) {
  if (field == "Sex") v <- ifelse(v %in% names(P$coding$Sex), P$coding$Sex[v], v)
  v <- ifelse(v == "Non_Improver", "Non-Improver", v)
  v <- ifelse(toupper(v) == "TRUE", "True", ifelse(toupper(v) == "FALSE", "False", v))
  v
}
same <- function(a, b) {
  na <- suppressWarnings(as.numeric(a)); nb <- suppressWarnings(as.numeric(b))
  ifelse(!is.na(na) & !is.na(nb), abs(na - nb) < 1e-6, a == b)
}

fields <- c("Treatment_arm", "Sex", "Age", "True_improver", "Improver", "mRSS_category", "Autoantibody_group",
            "Autoantibody", "Scl70", "RNApol1", "RNApol3", "CENPB", "Antibody_category",
            "Antibody_category_lab", "Ever_escaped", "Escape_3mo", "Escape_6mo", "Escape_9mo", "Escape_12mo")
counts <- list(); details <- list()
for (tg in c("skin", "pbmc")) {
  x <- get(tg)
  j <- merge(x, master, by = "Subject_ID", suffixes = c("", ".master"))
  for (f in intersect(fields, names(x))) {
    a <- harmonise(j[[f]], f); b <- harmonise(j[[paste0(f, ".master")]], f)
    kind <- fifelse(is.na(a) & is.na(b), "both missing",
            fifelse(is.na(a), "bulk missing, master has value",
            fifelse(is.na(b), "master missing, bulk has value",
            fifelse(same(a, b), "agree", "both present, differ"))))
    t <- as.data.table(table(kind))
    counts[[paste(tg, f)]] <- data.table(table = tg, field = f, kind = t$kind, rows = t$N,
                                         subjects = sapply(t$kind, function(k) uniqueN(j$Subject_ID[kind == k])))
    bad <- !kind %in% c("agree", "both missing")
    if (any(bad)) details[[paste(tg, f)]] <- data.table(table = tg, Sample_ID = j$Sample_ID[bad],
      Subject_ID = j$Subject_ID[bad], field = f, kind = kind[bad], bulk_value = j[[f]][bad],
      master_value = j[[paste0(f, ".master")]][bad])
  }
}
counts <- rbindlist(counts)
save_table(counts, run, "field_disagreements_after_harmonising")

# residual Sample_ID mismatch: metadata IDs absent from a skin matrix vs its extra columns
mat <- list.files(file.path(BULK, "skin"), pattern = "raw_counts.*EMSEMBL", full.names = TRUE)[1]
cols <- names(fread(mat, nrows = 0))
miss <- setdiff(skin$Sample_ID, cols)
extra <- setdiff(cols, skin$Sample_ID)[-1]                          # first column is the gene ID
d <- if (length(miss) && length(extra)) adist(miss, extra) else matrix(NA_integer_, 0, 0)
id_counts <- data.table(item = c("metadata IDs missing from matrix", "extra matrix columns (excl. gene column)",
                                 "minimum edit distance between them", "same ID shape"),
                        n = c(length(miss), length(extra), if (length(d)) min(d) else NA_integer_,
                              if (length(miss) && length(extra))
                                sum(outer(gsub("[0-9]", "9", gsub("[A-Za-z]", "A", miss)),
                                          gsub("[0-9]", "9", gsub("[A-Za-z]", "A", extra)), "==")) else NA_integer_))
save_table(id_counts, run, "residual_sample_id_mismatch")
details[["ids"]] <- data.table(table = "skin matrix", Sample_ID = miss, field = "Sample_ID",
                               kind = "in metadata, not in matrix",
                               bulk_value = paste(extra, collapse = "; "), master_value = NA_character_)

det <- rbindlist(details, fill = TRUE)
save_patient_table(det, run, "disagreement_details")

real <- counts[!kind %in% c("agree", "both missing")]
writeLines(c("Disagreements after harmonising coding (rows / subjects):",
             sprintf("  %-5s %-22s %-32s %4d rows  %3d subjects", real$table, real$field, real$kind, real$rows, real$subjects),
             "", "Residual Sample_ID mismatch (skin matrix):",
             sprintf("  %-45s %s", id_counts$item, id_counts$n),
             "", sprintf("ID-level details for the owner: %d rows", nrow(det))),
           file.path(run$dir, "MISMATCH_SUMMARY.txt"))
finalize_run(run, status = "ok", verdict = sprintf("%d field x kind disagreements remain after harmonising coding; %d ID rows for owner",
                                                   nrow(real), nrow(det)))
