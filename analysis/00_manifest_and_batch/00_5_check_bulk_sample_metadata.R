# ==============================================================
# 00_5_check_bulk_sample_metadata.R
#
# Consistency of the bulk sample metadata with the bulk matrices and with
# the clinical master (owner request 2026-10-02, items 2 and 6). Read-only.
# Scripts-only access (decisions.md 2026-10-02): this run reports COUNTS
# and category levels to tables/ (agent-readable); every row with an ID is
# written by save_patient_table() to data/patient_level/ for the owner.
#
# Checks
#   A. Sample_ID: skin metadata vs the columns of each skin matrix; PBMC
#      metadata vs each blood matrix. Exact match, and match after
#      normalising (trim, case, "-" "." "_" " " removed) — the second
#      catches typos of separator/case kind.
#   B. Subject_ID: skin and PBMC metadata vs the master (88).
#   C. PBMC Matched_Skin_Sample_ID: exists in the skin metadata, same
#      subject, and which skin timepoint it points to.
#   D. Shared clinical fields: skin / PBMC rows vs the master row of the
#      same subject (static fields; mRSS, HAQ-DI, PGA, PtGA, FVC, DLCO at
#      the sample's timepoint).
#   E. ID shapes (letters -> A, digits -> 9) and stray whitespace.
#   F. Skin design: samples per Subject x Timepoint before / after the
#      owner's Exclude and Repeat filters.
#   G. Levels of every low-cardinality column (<= 12 values), counts only.
#
#   jobs/run.sh analysis/00_manifest_and_batch/00_5_check_bulk_sample_metadata.R \
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
          level_max = 12, num_tol = 1e-6)

run <- init_run(
  stage    = "00_manifest_and_batch",
  run_name = "check_bulk_sample_metadata",
  params   = P,
  notes    = "Sample/Subject IDs across bulk metadata, matrices and the clinical master; field consistency"
)

rd <- function(f, sheet = 1) {
  x <- as.data.table(read_excel(file.path(CLINICAL, f), sheet = sheet, col_types = "text"))
  x[, names(x) := lapply(.SD, function(v) { v[!is.na(v) & !nzchar(trimws(v))] <- NA; v })]
  x
}
master <- rd(basename(CLINICAL_MASTER))
skin   <- rd(P$skin_file)
pbmc   <- rd(P$pbmc_file)
stopifnot(nrow(master) == 88, !anyDuplicated(master$Subject_ID))

norm_id <- function(x) toupper(gsub("[-._ ]", "", trimws(x)))
shape   <- function(x) gsub("[0-9]", "9", gsub("[A-Za-z]", "A", x))
counts <- list(); patient <- list()
add <- function(check, item, n, of = NA_integer_) {
  counts[[length(counts) + 1]] <<- data.table(check = check, item = item, n = n, of = of)
}

# ---- A. Sample_IDs vs matrices --------------------------------

mat_cols <- function(sub) {
  fs <- list.files(file.path(BULK, sub), pattern = "\\.txt$", full.names = TRUE)
  fs <- fs[!grepl("sample_info|batch_info", fs)]
  setNames(lapply(fs, function(f) names(fread(f, nrows = 0))), basename(fs))
}
check_matrix_ids <- function(meta_ids, mats, tag) {
  for (m in names(mats)) {
    cols <- mats[[m]]
    sample_cols <- cols[cols %in% meta_ids | norm_id(cols) %in% norm_id(meta_ids)]
    other <- setdiff(cols, sample_cols)
    add(paste0("A_", tag), paste(m, "metadata IDs found exactly"), sum(meta_ids %in% cols), length(meta_ids))
    add(paste0("A_", tag), paste(m, "metadata IDs found only after normalising"),
        sum(!meta_ids %in% cols & norm_id(meta_ids) %in% norm_id(cols)), length(meta_ids))
    add(paste0("A_", tag), paste(m, "metadata IDs not found"), sum(!norm_id(meta_ids) %in% norm_id(cols)), length(meta_ids))
    add(paste0("A_", tag), paste(m, "matrix columns not in metadata (incl. gene/annotation columns)"),
        length(other), length(cols))
    add(paste0("A_", tag), paste(m, "duplicated column names"), sum(duplicated(cols)), length(cols))
    patient[[paste("A", tag, m)]] <<- rbind(
      data.table(matrix = m, side = "metadata ID not exact in matrix",
                 id = meta_ids[!meta_ids %in% cols],
                 normalised_match = cols[match(norm_id(meta_ids[!meta_ids %in% cols]), norm_id(cols))]),
      data.table(matrix = m, side = "matrix column not in metadata", id = other, normalised_match = NA_character_))
  }
}
check_matrix_ids(skin$Sample_ID, mat_cols("skin"), "skin")
check_matrix_ids(pbmc$Sample_ID, mat_cols("blood"), "blood")

# ---- B. Subject_IDs vs master ---------------------------------

for (tg in c("skin", "pbmc")) {
  x <- get(tg)
  s <- unique(x$Subject_ID)
  add("B_subject", paste(tg, "subjects in master (exact)"), sum(s %in% master$Subject_ID), length(s))
  add("B_subject", paste(tg, "subjects in master only after normalising"),
      sum(!s %in% master$Subject_ID & norm_id(s) %in% norm_id(master$Subject_ID)), length(s))
  bad <- s[!s %in% master$Subject_ID]
  if (length(bad)) patient[[paste("B", tg)]] <- data.table(table = tg, Subject_ID = bad,
    normalised_match = master$Subject_ID[match(norm_id(bad), norm_id(master$Subject_ID))])
}
add("B_subject", "master subjects with >= 1 skin sample", sum(master$Subject_ID %in% skin$Subject_ID), 88L)
add("B_subject", "master subjects with a PBMC sample", sum(master$Subject_ID %in% pbmc$Subject_ID), 88L)
dup_pbmc <- pbmc[, .N, by = Subject_ID][N > 1]
add("B_subject", "PBMC subjects with > 1 sample", nrow(dup_pbmc), uniqueN(pbmc$Subject_ID))
if (nrow(dup_pbmc)) patient[["B pbmc dup"]] <- pbmc[Subject_ID %in% dup_pbmc$Subject_ID,
                                                    .(table = "pbmc duplicate subject", Subject_ID, Sample_ID)]

# ---- C. PBMC -> matched skin sample ---------------------------

m <- merge(pbmc[, .(pbmc_Sample_ID = Sample_ID, Subject_ID, Matched_Skin_Sample_ID)],
           skin[, .(Matched_Skin_Sample_ID = Sample_ID, skin_Subject_ID = Subject_ID, skin_Timepoint = Timepoint,
                    skin_Exclude = Exclude, skin_Repeat = Repeat)],
           by = "Matched_Skin_Sample_ID", all.x = TRUE)
add("C_match", "PBMC rows whose matched skin sample exists", sum(!is.na(m$skin_Subject_ID)), nrow(m))
add("C_match", "matched skin sample from the same subject", sum(m$skin_Subject_ID == m$Subject_ID, na.rm = TRUE), nrow(m))
add("C_match", "duplicated Matched_Skin_Sample_ID", sum(duplicated(m$Matched_Skin_Sample_ID)), nrow(m))
for (tp in sort(unique(na.omit(m$skin_Timepoint)))) add("C_match", paste("matched skin timepoint =", tp), sum(m$skin_Timepoint == tp, na.rm = TRUE), nrow(m))
add("C_match", "matched skin sample is Exclude/Repeat-flagged (see G for levels)",
    sum(!is.na(m$skin_Exclude) & (m$skin_Exclude != names(which.max(table(skin$Exclude))) |
                                  m$skin_Repeat != names(which.max(table(skin$Repeat))))), nrow(m))
patient[["C"]] <- m[is.na(skin_Subject_ID) | skin_Subject_ID != Subject_ID | duplicated(Matched_Skin_Sample_ID) |
                    duplicated(Matched_Skin_Sample_ID, fromLast = TRUE)]

# ---- D. shared fields vs master -------------------------------

static <- c("Treatment_arm", "Sex", "Age", "True_improver", "Improver", "mRSS_category", "Autoantibody_group",
            "Autoantibody", "Scl70", "RNApol1", "RNApol3", "CENPB", "Antibody_category",
            "Antibody_category_lab", "Ever_escaped", "Escape_3mo", "Escape_6mo", "Escape_9mo", "Escape_12mo")
same_value <- function(a, b) {
  na <- suppressWarnings(as.numeric(a)); nb <- suppressWarnings(as.numeric(b))
  ifelse(!is.na(na) & !is.na(nb), abs(na - nb) <= P$num_tol,
         ifelse(is.na(a) & is.na(b), TRUE, toupper(trimws(a)) == toupper(trimws(b))))
}
compare_fields <- function(x, tag, pairs) {
  j <- merge(x, master, by = "Subject_ID", suffixes = c("", ".master"))
  for (p in pairs) {
    a <- j[[p[1]]]; b <- j[[p[2]]]
    if (is.null(a) || is.null(b)) { add(paste0("D_", tag), paste(p[1], "vs master", p[2], ": column missing"), NA_integer_); next }
    ok <- same_value(a, b); ok[is.na(ok)] <- FALSE
    add(paste0("D_", tag), paste(p[1], "vs master", p[2], "— rows differing"), sum(!ok), length(ok))
    if (any(!ok)) patient[[paste("D", tag, p[1], p[2])]] <<- data.table(
      table = tag, Sample_ID = j$Sample_ID[!ok], Subject_ID = j$Subject_ID[!ok],
      field = p[1], value = a[!ok], master_field = p[2], master_value = b[!ok])
  }
}
pairs_static <- lapply(intersect(static, names(skin)), function(f) c(f, paste0(f, ".master")))
compare_fields(skin, "skin_static", pairs_static)

# timepoint-specific: skin Time (0/3/6) -> master column
tp_cols <- list(MRSS = c(`0` = "0", `3` = "3", `6` = "6"),
                HAQ_DI = c(`0` = "HAQ_DI_Month0", `3` = "HAQ_DI_Month3", `6` = "HAQ_DI_Month6"),
                PGA_Physician = c(`0` = "PGA_Physician_Month0", `3` = "PGA_Physician_Month3", `6` = "PGA_Physician_Month6"),
                PtGA_Patient = c(`0` = "PtGA_Patient_Month0", `3` = "PtGA_Patient_Month3", `6` = "PtGA_Patient_Month6"),
                FVC = c(`0` = "FVCPTP_screening", `3` = "FVCPTP_Month3", `6` = "FVCPTP_Month6"),
                DLCO = c(`0` = "DLCO_corrected_Screening", `6` = "DLCO_corrected_Month6"))
add("D_skin_time", "skin Time values", uniqueN(skin$Time))
for (f in names(tp_cols)) {
  map <- tp_cols[[f]]
  x <- skin[Time %in% names(map), .(Sample_ID, Subject_ID, Time, value = get(f))]
  x[, master_col := map[as.character(Time)]]
  x[, master_value := mapply(function(s, mc) { r <- master[[mc]][master$Subject_ID == s]
                                                if (length(r)) r[1] else NA_character_ },
                              Subject_ID, master_col, USE.NAMES = FALSE)]
  ok <- same_value(x$value, x$master_value); ok[is.na(ok)] <- FALSE
  add("D_skin_time", paste(f, "at sample timepoint vs master — rows differing"), sum(!ok), nrow(x))
  if (any(!ok)) patient[[paste("D time", f)]] <- x[!ok, .(table = "skin", Sample_ID, Subject_ID, field = f,
                                                          value, master_field = master_col, master_value)]
}
# PBMC: site vs master Site_ID
compare_fields(pbmc, "pbmc", list(c("site", "Site_ID")))

# ---- E. ID shapes and whitespace -------------------------------

for (tg in c("master", "skin", "pbmc")) {
  x <- get(tg)
  for (col in intersect(c("Subject_ID", "Sample_ID", "Matched_Skin_Sample_ID"), names(x))) {
    v <- x[[col]]
    sh <- as.data.table(table(shape(v)))
    for (i in seq_len(nrow(sh))) add(paste("E_shape", tg, col), sh$V1[i], sh$N[i], length(v))
    add(paste("E_space", tg, col), "values with leading/trailing/inner spaces", sum(grepl("^\\s|\\s$|\\s", v)), length(v))
  }
}

# ---- F. skin design: one sample per subject x timepoint? -------

lv_main <- function(v) names(which.max(table(v)))      # majority level = "keep"
keep_ex <- lv_main(skin$Exclude); keep_rp <- lv_main(skin$Repeat)
f0 <- skin[, .N, by = .(Subject_ID, Timepoint)]
f1 <- skin[Exclude == keep_ex, .N, by = .(Subject_ID, Timepoint)]
f2 <- skin[Exclude == keep_ex & Repeat == keep_rp, .N, by = .(Subject_ID, Timepoint)]
add("F_design", sprintf("majority (kept) levels: Exclude = '%s', Repeat = '%s'", keep_ex, keep_rp), NA_integer_)
for (st in list(list("all samples", f0), list("after Exclude", f1), list("after Exclude + Repeat", f2))) {
  d <- st[[2]]
  add("F_design", paste(st[[1]], ": samples"), sum(d$N))
  add("F_design", paste(st[[1]], ": subject x timepoint units"), nrow(d))
  add("F_design", paste(st[[1]], ": units with > 1 sample"), sum(d$N > 1))
  add("F_design", paste(st[[1]], ": subjects"), uniqueN(d$Subject_ID))
  for (tp in sort(unique(d$Timepoint))) add("F_design", paste(st[[1]], ": subjects at", tp), sum(d$Timepoint == tp))
}
add("F_design", "Timepoint vs Time combinations", nrow(unique(skin[, .(Timepoint, Time)])))
add("F_design", "Repeat-flagged samples whose unit has no other sample",
    nrow(merge(skin[Repeat != keep_rp, .(Subject_ID, Timepoint)], f0[N == 1], by = c("Subject_ID", "Timepoint"))))
patient[["F units >1 after filters"]] <- merge(skin[, .(Sample_ID, Subject_ID, Timepoint, Exclude, Repeat, replicate)],
                                              f2[N > 1, .(Subject_ID, Timepoint)], by = c("Subject_ID", "Timepoint"))

# ---- G. category levels (counts only) --------------------------

lv <- rbindlist(lapply(c("master", "skin", "pbmc"), function(tg) {
  x <- get(tg)
  rbindlist(lapply(names(x), function(col) {
    v <- x[[col]]
    if (uniqueN(v) > P$level_max || grepl("ID$|note", col, ignore.case = TRUE)) return(NULL)   # no IDs, no free text
    t <- as.data.table(table(v, useNA = "ifany"))
    data.table(table = tg, column = col, level = as.character(t$v), n = t$N)
  }))
}))
save_table(lv, run, "category_levels")

# ---- outputs ---------------------------------------------------

counts <- rbindlist(counts)
save_table(counts, run, "consistency_counts")
pt <- rbindlist(patient, fill = TRUE, idcol = "check")
save_patient_table(pt, run, "consistency_details")

writeLines(c(sprintf("%-14s %-90s %6s / %s", counts$check, counts$item, counts$n, counts$of),
             "", sprintf("ID-level details for the owner: %d rows (data/patient_level/.../consistency_details.csv)", nrow(pt))),
           file.path(run$dir, "CONSISTENCY_SUMMARY.txt"))
finalize_run(run, status = "ok", verdict = sprintf("%d checks; %d ID-level rows for owner review", nrow(counts), nrow(pt)))
