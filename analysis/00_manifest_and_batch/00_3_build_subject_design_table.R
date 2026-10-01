# ==============================================================
# 00_3_build_subject_design_table.R
#
# One row per subject: the design fields every later stage joins on.
#
#   Subject_ID, group, arm, mRSS_category, Ever_escaped, Escape_month,
#   scRNA-seq timepoints present, bulk skin / bulk PBMC availability
#
# Written to metadata/design/ (gitignored). By the owner's decision
# (docs/decisions.md, 2026-10-01) this table — and only this table — is
# readable by AI agents. The clinical master table itself stays in
# data/clinical/ and is read only by this script.
#
# Also writes metadata/design/bulk_file_inventory.csv: the file names
# and sizes under data/bulk/ (names only, no contents), so the bulk
# import stages can be written against the real layout.
#
# Fails loudly if the clinical table does not have 88 subjects
# (44 Placebo, 44 Abatacept), if an mRSS_category value is unexpected,
# or if any scRNA-seq subject is missing from the clinical table or
# carries a different arm there. ID mismatches are written to
# metadata/design/id_match_report.csv before stopping.
#
#   jobs/run.sh analysis/00_manifest_and_batch/00_3_build_subject_design_table.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")

N_SUBJECTS   <- 88
N_PER_ARM    <- c(Placebo = 44, Abatacept = 44)
MRSS_LEVELS  <- c("Improver", "Stable", "Worsened", "Set_aside")
KEEP_COLS    <- c("Subject_ID", "Treatment_arm", "mRSS_category",
                  "Ever_escaped", "Escape_month")

run <- init_run(
  stage    = "00_manifest_and_batch",
  run_name = "build_subject_design_table",
  params   = list(clinical_master = .rel_to_repo(CLINICAL_MASTER),
                  n_subjects = N_SUBJECTS, n_per_arm = as.list(N_PER_ARM),
                  mrss_levels = MRSS_LEVELS),
  notes    = "Subject design table (arm, mRSS_category, escape, data availability)"
)

fail <- function(msg) {
  finalize_run(run, status = "failed", verdict = substr(msg, 1, 200))
  stop(msg, call. = FALSE)
}

dir.create(DESIGN, showWarnings = FALSE, mode = "0700")
write_design <- function(x, name) {
  f <- file.path(DESIGN, paste0(name, ".csv"))
  utils::write.csv(x, f, row.names = FALSE, na = "")
  Sys.chmod(f, "0600")
  cat(.rel_to_repo(f), "\n", file = file.path(run$tables, "DESIGN_TABLES.txt"),
      append = TRUE)
  message("  wrote ", .rel_to_repo(f))
  invisible(f)
}

# ---- 1. clinical master table --------------------------------
# All columns as text: keeps IDs exactly as typed and avoids type
# guessing. Only the five design columns are kept.

if (!file.exists(CLINICAL_MASTER)) fail(paste("not found:", CLINICAL_MASTER))
clin <- readxl::read_excel(CLINICAL_MASTER, sheet = 1, col_types = "text",
                           .name_repair = "minimal")
clin <- as.data.frame(clin, stringsAsFactors = FALSE)
miss <- setdiff(KEEP_COLS, names(clin))
if (length(miss)) fail(paste("clinical table lacks columns:", paste(miss, collapse = ", ")))
clin <- clin[, KEEP_COLS]
clin <- clin[!is.na(clin$Subject_ID) & nzchar(trimws(clin$Subject_ID)), ]
clin$Subject_ID <- trimws(clin$Subject_ID)

if (anyDuplicated(clin$Subject_ID)) fail("duplicate Subject_ID in clinical table")
if (nrow(clin) != N_SUBJECTS) {
  fail(sprintf("clinical table has %d subjects, expected %d", nrow(clin), N_SUBJECTS))
}

clin$arm <- arm_label(clin$Treatment_arm)
if (anyNA(clin$arm)) {
  fail(paste("unrecognised Treatment_arm values:",
             paste(unique(clin$Treatment_arm[is.na(clin$arm)]), collapse = ", ")))
}
n_arm <- table(factor(clin$arm, levels = names(N_PER_ARM)))
if (!all(n_arm == N_PER_ARM)) {
  fail(sprintf("arm counts %s, expected %s",
               paste(names(n_arm), n_arm, collapse = " / "),
               paste(names(N_PER_ARM), N_PER_ARM, collapse = " / ")))
}

clin$mRSS_category <- trimws(clin$mRSS_category)
bad_cat <- setdiff(unique(stats::na.omit(clin$mRSS_category)), MRSS_LEVELS)
if (length(bad_cat)) {
  fail(paste("unexpected mRSS_category values:", paste(bad_cat, collapse = ", ")))
}

esc_txt <- trimws(clin$Escape_month)
esc_txt[esc_txt %in% c("", "NA", "N/A", "na")] <- NA
esc_num <- suppressWarnings(as.numeric(esc_txt))
if (any(!is.na(esc_txt) & is.na(esc_num))) {
  fail(paste("non-numeric Escape_month values:",
             paste(unique(esc_txt[!is.na(esc_txt) & is.na(esc_num)]), collapse = ", ")))
}
clin$Escape_month <- esc_num

# ---- 2. scRNA-seq availability, from the manifest ------------

man <- read_manifest()
ssc_lib <- man[man$group == "SSc", ]
hc_lib  <- man[man$group == "HC", ]

sc_avail <- do.call(rbind, lapply(split(ssc_lib, ssc_lib$Subject_ID), function(d) {
  ok <- is.na(d$design_flag) | !nzchar(d$design_flag)
  data.frame(Subject_ID      = d$Subject_ID[1],
             sc_arm_manifest = d$arm[1],
             sc_M00 = any(d$timepoint == "M00" & ok),
             sc_M03 = any(d$timepoint == "M03" & ok),
             sc_M06 = any(d$timepoint == "M06" & ok),
             sc_n_libraries = nrow(d),
             sc_flagged     = paste(unique(d$design_flag[!ok]), collapse = ";"),
             stringsAsFactors = FALSE)
}))

# ID matching report: written before any failure so the mismatch can be
# inspected.
issue_rows <- function(ids, issue) {
  data.frame(Subject_ID = as.character(ids), issue = rep(issue, length(ids)),
             stringsAsFactors = FALSE)
}
j <- match(sc_avail$Subject_ID, clin$Subject_ID)
arm_mismatch <- !is.na(j) & sc_avail$sc_arm_manifest != clin$arm[j]
id_report <- rbind(
  issue_rows(setdiff(sc_avail$Subject_ID, clin$Subject_ID),
             "in scRNA-seq manifest, not in clinical table"),
  issue_rows(setdiff(clin$Subject_ID, sc_avail$Subject_ID),
             "in clinical table, no scRNA-seq library in this freeze"),
  issue_rows(sc_avail$Subject_ID[arm_mismatch],
             "arm differs between scRNA-seq manifest and clinical table"))
write_design(id_report, "id_match_report")

n_sc_unmatched <- sum(!sc_avail$Subject_ID %in% clin$Subject_ID)
if (n_sc_unmatched > 0) {
  fail(sprintf(paste("%d of %d scRNA-seq subjects not found in the clinical",
                     "table (ID format?); see metadata/design/id_match_report.csv"),
               n_sc_unmatched, nrow(sc_avail)))
}
if (any(arm_mismatch)) {
  fail(sprintf("%d subjects have a different arm in the manifest and the clinical table",
               sum(arm_mismatch)))
}

# ---- 3. bulk file inventory (names and sizes only) -----------
# Per-subject bulk availability is added once the bulk layout is known
# (stages 01 / 02); until then the columns are NA.

bulk_files <- list.files(BULK, recursive = TRUE, all.files = FALSE)
bulk_inv <- data.frame(path  = bulk_files,
                       bytes = file.size(file.path(BULK, bulk_files)),
                       stringsAsFactors = FALSE)
write_design(bulk_inv, "bulk_file_inventory")

# ---- 4. design table -----------------------------------------

design <- data.frame(Subject_ID    = clin$Subject_ID,
                     group         = "SSc",
                     arm           = clin$arm,
                     mRSS_category = clin$mRSS_category,
                     Ever_escaped  = clin$Ever_escaped,
                     Escape_month  = clin$Escape_month,
                     stringsAsFactors = FALSE)
k <- match(design$Subject_ID, sc_avail$Subject_ID)
for (cl in c("sc_M00", "sc_M03", "sc_M06")) design[[cl]] <- !is.na(k) & sc_avail[[cl]][k] %in% TRUE
design$sc_n_libraries <- ifelse(is.na(k), 0L, sc_avail$sc_n_libraries[k])
design$sc_flagged     <- ifelse(is.na(k), "", sc_avail$sc_flagged[k])
design$bulk_skin      <- NA
design$bulk_pbmc_M00  <- NA

hc <- data.frame(Subject_ID = unique(hc_lib$Subject_ID), group = "HC",
                 arm = NA, mRSS_category = NA, Ever_escaped = NA,
                 Escape_month = NA, sc_M00 = FALSE, sc_M03 = FALSE,
                 sc_M06 = FALSE, sc_n_libraries = NA, sc_flagged = "",
                 bulk_skin = NA, bulk_pbmc_M00 = NA, stringsAsFactors = FALSE)
hc$sc_n_libraries <- as.integer(table(hc_lib$Subject_ID)[hc$Subject_ID])

design <- rbind(design, hc)
stopifnot(!anyDuplicated(design$Subject_ID))
write_design(design, "subject_design")

# ---- 5. aggregate summaries (no IDs) into results/ -----------

ssc <- design[design$group == "SSc", ]
ssc$mRSS_category[is.na(ssc$mRSS_category)] <- "(missing)"
cat_by_arm <- as.data.frame.matrix(table(ssc$mRSS_category, ssc$arm))
save_table(cbind(mRSS_category = rownames(cat_by_arm), cat_by_arm), run,
           "mrss_category_by_arm")
sc_counts <- do.call(rbind, lapply(split(ssc, ssc$arm), function(d)
  data.frame(arm = d$arm[1], subjects = nrow(d),
             sc_any = sum(d$sc_n_libraries > 0), sc_M00 = sum(d$sc_M00),
             sc_M03 = sum(d$sc_M03), sc_M06 = sum(d$sc_M06),
             escaped = sum(!is.na(d$Escape_month)))))
save_table(sc_counts, run, "sc_availability_by_arm")

pl <- ssc[ssc$arm == "Placebo", ]
verdict <- sprintf(paste("%d subjects (44/44) + %d HC; %d SSc with scRNA-seq;",
                         "Placebo mRSS_category: %s; %d bulk files inventoried"),
                   nrow(ssc), nrow(hc), sum(ssc$sc_n_libraries > 0),
                   paste(names(table(pl$mRSS_category)), table(pl$mRSS_category),
                         sep = "=", collapse = " "),
                   nrow(bulk_inv))
finalize_run(run, status = "ok", verdict = verdict)
