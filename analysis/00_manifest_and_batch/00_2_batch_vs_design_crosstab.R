# ==============================================================
# 00_2_batch_vs_design_crosstab.R
#
# HAZARD #1. This is the first analysis in the project, not a QC
# footnote.
#
# If library prep / capture batch is correlated with timepoint or arm,
# the longitudinal comparison at the centre of Project 1 is confounded
# and no amount of integration rescues it. Harmony can remove a batch
# effect; it cannot tell you whether the difference it removed was
# batch or biology when the two are collinear.
#
# What this script does:
#   1. crosstabs batch_id against timepoint, arm, and group
#   2. tests each association (chi-square / Fisher)
#   3. computes Cramer's V as an effect size
#   4. checks whether each SUBJECT's serial samples sit in one batch
#      or are split across batches
#   5. writes a verdict you can act on
#
# Interpretation guide:
#   - batch vs timepoint independent        -> longitudinal design is clean
#   - each subject's samples in ONE batch   -> batch is nested in subject;
#                                              subject random effect absorbs it
#   - batch tracks timepoint                -> STOP. The M0 -> M6 comparison
#                                              cannot be separated from batch.
#   - batch tracks arm                      -> arm comparisons compromised
#   - HC in their own batches               -> SSc-vs-HC differences are
#                                              partly technical
#
# Run:
#   jobs/run.sh analysis/00_manifest_and_batch/00_2_batch_vs_design_crosstab.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")
source("R/theme_asset.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

run <- init_run(
  stage    = "00_manifest_and_batch",
  run_name = "batch_vs_design_crosstab",
  notes    = paste("Hazard #1: is batch_id confounded with timepoint,",
                   "arm, or group?")
)

# ---- 1. manifest, restricted to the freeze --------------------

man <- read_manifest()
lib <- freeze_libraries(FREEZE, man)

message(sprintf("\nFreeze %s: %d libraries, %d subjects",
                FREEZE, nrow(lib), length(unique(lib$subject_id))))

if (!"batch_id" %in% names(lib) || all(is.na(lib$batch_id))) {
  finalize_run(run, status = "failed",
               verdict = "batch_id absent from the manifest — get it from the core")
  stop("No batch_id. Ask the sequencing core for the prep/capture batch ",
       "assignment per library. Without it, hazard #1 cannot be assessed ",
       "and no longitudinal result is trustworthy.", call. = FALSE)
}

# ---- 2. crosstabs and association tests -----------------------

cramers_v <- function(tbl) {
  chi <- suppressWarnings(stats::chisq.test(tbl)$statistic)
  n   <- sum(tbl)
  k   <- min(dim(tbl)) - 1
  if (k <= 0 || n == 0) return(NA_real_)
  unname(sqrt(chi / (n * k)))
}

test_assoc <- function(tbl, name) {
  tbl <- tbl[rowSums(tbl) > 0, colSums(tbl) > 0, drop = FALSE]
  if (any(dim(tbl) < 2)) {
    return(data.frame(comparison = name, test = "none",
                      p_value = NA_real_, cramers_v = NA_real_,
                      note = "only one level present",
                      stringsAsFactors = FALSE))
  }
  small <- any(suppressWarnings(stats::chisq.test(tbl)$expected) < 5)
  if (small) {
    p <- tryCatch(
      stats::fisher.test(tbl, simulate.p.value = TRUE, B = 20000)$p.value,
      error = function(e) NA_real_)
    tst <- "Fisher (simulated)"
  } else {
    p <- suppressWarnings(stats::chisq.test(tbl)$p.value)
    tst <- "chi-square"
  }
  data.frame(comparison = name, test = tst, p_value = p,
             cramers_v = cramers_v(tbl),
             note = if (small) "sparse table" else "",
             stringsAsFactors = FALSE)
}

design_vars <- intersect(c("timepoint", "arm", "group", "run_id",
                           "cellranger_run"), names(lib))
design_vars <- setdiff(design_vars, "batch_id")

tests <- do.call(rbind, lapply(design_vars, function(v) {
  tbl <- table(lib$batch_id, lib[[v]], useNA = "ifany")
  utils::write.csv(as.data.frame.matrix(tbl),
                   file.path(run$tables,
                             sprintf("crosstab_batch_by_%s.csv", v)))
  test_assoc(tbl, paste("batch_id vs", v))
}))

save_table(tests, run, "batch_association_tests")
print(tests, row.names = FALSE)

# ---- 3. is batch nested within subject? -----------------------
# If every subject's serial samples were prepped in one batch, batch is
# nested in subject and the subject random effect absorbs it. That is
# the good case. If subjects are split across batches, batch enters the
# within-subject comparison directly.

per_subject <- do.call(rbind, lapply(
  split(lib, lib$subject_id), function(d) {
    data.frame(subject_id   = d$subject_id[1],
               group        = d$group[1],
               arm          = if ("arm" %in% names(d)) d$arm[1] else NA,
               n_libraries  = nrow(d),
               n_batches    = length(unique(d$batch_id)),
               batches      = paste(sort(unique(d$batch_id)),
                                    collapse = ";"),
               timepoints   = paste(sort(unique(d$timepoint)),
                                    collapse = ";"),
               stringsAsFactors = FALSE)
  }))
save_table(per_subject, run, "batch_nesting_per_subject")

multi <- per_subject[per_subject$n_libraries > 1, ]
n_split <- sum(multi$n_batches > 1)
pct_split <- if (nrow(multi)) 100 * n_split / nrow(multi) else NA

message(sprintf(
  "\nSubjects with >1 library: %d; of those, %d (%.0f%%) span multiple batches",
  nrow(multi), n_split, pct_split))

# ---- 4. plots -------------------------------------------------

for (v in design_vars) {
  df <- as.data.frame(table(batch = lib$batch_id, level = lib[[v]]))
  p <- ggplot(df, aes(batch, Freq, fill = level)) +
    geom_col(position = "stack", width = 0.8) +
    scale_fill_manual(values = pal_asset(length(unique(df$level))),
                      name = v) +
    labs(title = paste("Libraries per batch, coloured by", v),
         subtitle = sprintf("%s: p = %s, Cramer's V = %s",
                            tests$test[tests$comparison ==
                                         paste("batch_id vs", v)],
                            signif(tests$p_value[tests$comparison ==
                                                   paste("batch_id vs", v)], 3),
                            signif(tests$cramers_v[tests$comparison ==
                                                     paste("batch_id vs", v)], 2)),
         x = "batch_id", y = "libraries") +
    theme_asset() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, run, paste0("QC_batch_by_", v), qc = TRUE,
            width = 7, height = 3.5)
}

p_nest <- ggplot(multi, aes(factor(n_batches))) +
  geom_bar(fill = "#4C72B0", width = 0.7) +
  labs(title = "Batches per subject (subjects with >1 library)",
       subtitle = "1 = batch nested in subject (good); >1 = batch enters the within-subject comparison",
       x = "distinct batches", y = "subjects") +
  theme_asset()
save_plot(p_nest, run, "QC_batch_nesting_per_subject", qc = TRUE,
          width = 4, height = 3)

# ---- 5. verdict -----------------------------------------------

alpha <- 0.05
v_thresh <- 0.3   # Cramer's V above this = substantively associated

flagged <- tests[!is.na(tests$p_value) &
                 (tests$p_value < alpha |
                  (!is.na(tests$cramers_v) & tests$cramers_v > v_thresh)), ]

verdict_lines <- c(
  sprintf("Freeze: %s | %d libraries | %d subjects",
          FREEZE, nrow(lib), length(unique(lib$subject_id))),
  sprintf("Subjects spanning multiple batches: %d of %d (%.0f%%)",
          n_split, nrow(multi), pct_split),
  "",
  "Associations tested:",
  utils::capture.output(print(tests, row.names = FALSE)),
  ""
)

if (nrow(flagged) == 0) {
  status  <- "ok"
  verdict <- "No batch-design confounding detected; longitudinal design is clean"
  verdict_lines <- c(verdict_lines,
    "VERDICT: no significant or substantive association between batch and",
    "any design variable. Proceed. Still integrate on batch_id only.")
} else {
  status  <- "ok"
  verdict <- paste0("BATCH CONFOUNDING FLAGGED: ",
                    paste(flagged$comparison, collapse = "; "),
                    " — see run directory before any longitudinal analysis")
  verdict_lines <- c(verdict_lines,
    "VERDICT: batch is associated with design variable(s):",
    paste("  -", flagged$comparison,
          sprintf("(p = %s, V = %s)",
                  signif(flagged$p_value, 3),
                  signif(flagged$cramers_v, 2))),
    "",
    "Do NOT proceed to biological interpretation of the affected contrast",
    "until this is resolved. Options, in order of preference:",
    "  1. If batch is nested within subject, a subject random effect",
    "     absorbs it — verify with tables/batch_nesting_per_subject.csv.",
    "  2. Include batch_id as a covariate in every pseudobulk model, and",
    "     check the design matrix is not rank-deficient.",
    "  3. Restrict to the subset of batches that are balanced across the",
    "     contrast, and report the reduced n.",
    "  4. If batch is collinear with timepoint, the longitudinal claim",
    "     cannot be made from this freeze. Say so explicitly.",
    "",
    "Record the chosen option in docs/decisions.md.")
}

writeLines(verdict_lines, file.path(run$dir, "VERDICT.txt"))
cat("\n", paste(verdict_lines, collapse = "\n"), "\n", sep = "")

finalize_run(run, status = status, verdict = verdict)
