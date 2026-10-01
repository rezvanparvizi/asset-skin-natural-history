# ==============================================================
# 00_1_build_manifest.R
#
# Validate the library manifest and freeze the resolved library list,
# so the freeze stays reproducible even if the manifest changes later.
#
# Run FIRST, before anything else in the project.
#
#   jobs/run.sh analysis/00_manifest_and_batch/00_1_build_manifest.R \
#     --freeze freeze01 --cohort reference
# ==============================================================

source("config/paths.R")
source("R/provenance.R")
source("R/io.R")

run <- init_run(
  stage    = "00_manifest_and_batch",
  run_name = "build_manifest",
  notes    = "Validate manifest, resolve freeze membership, report cohort sizes"
)

man <- read_manifest()
message(sprintf("\nManifest: %d libraries", nrow(man)))

# ---- structural checks ---------------------------------------

issues <- character(0)

# one sample_id should map to one Subject_ID
bad_map <- names(which(vapply(split(man$Subject_ID, man$sample_id),
                              function(x) length(unique(x)), integer(1)) > 1))
if (length(bad_map)) {
  issues <- c(issues, paste("sample_id mapping to >1 Subject_ID:",
                            paste(bad_map, collapse = ", ")))
}

# resequenced libraries must share their parent's sample_id
res <- man[!is.na(man$resequenced_of) & nzchar(man$resequenced_of), ]
if (nrow(res)) {
  parent_sample <- man$sample_id[match(res$resequenced_of, man$library_id)]
  mism <- res$library_id[!is.na(parent_sample) &
                         parent_sample != res$sample_id]
  if (length(mism)) {
    issues <- c(issues, paste("resequenced library with a different",
                              "sample_id than its parent:",
                              paste(mism, collapse = ", ")))
  }
  orphan <- res$library_id[is.na(parent_sample)]
  if (length(orphan)) {
    issues <- c(issues, paste("resequenced_of points at an unknown",
                              "library_id:", paste(orphan, collapse = ", ")))
  }
  message(sprintf("  %d resequenced libraries present", nrow(res)))
}

# every resequenced library needs an explicit merge_strategy
if ("merge_strategy" %in% names(man) && nrow(res)) {
  nostrat <- res$library_id[is.na(res$merge_strategy) |
                            !nzchar(res$merge_strategy)]
  if (length(nostrat)) {
    issues <- c(issues,
      paste("resequenced library with no merge_strategy (keep/exclude/",
            "replace/merge_reads):", paste(nostrat, collapse = ", ")))
  }
}

# healthy controls should have no arm or timepoint
hc <- man[man$group == "HC", ]
if (nrow(hc)) {
  odd <- hc$library_id[!(is.na(hc$arm) | hc$arm %in% c("NA", ""))]
  if (length(odd)) {
    issues <- c(issues, paste("HC library with a treatment arm:",
                              paste(odd, collapse = ", ")))
  }
  message(sprintf("  %d healthy control libraries", nrow(hc)))
}

if (length(issues)) {
  for (i in issues) message("  ISSUE: ", i)
} else {
  message("  structural checks passed")
}

# ---- resolve freeze membership --------------------------------

lib <- freeze_libraries(FREEZE, man)
out <- file.path(METADATA, sprintf("%s_libraries.txt", FREEZE))
writeLines(sort(lib$library_id), out)
message(sprintf("\nFreeze %s resolves to %d libraries -> %s",
                FREEZE, nrow(lib), .rel_to_repo(out)))

# ---- cohort sizes ---------------------------------------------
# Reported without clinical data, so post-escape censoring is NOT
# applied here. Cohorts needing it will differ; that is expected.

coh_files <- list.files(file.path(CONFIG, "cohorts"), "\\.ya?ml$",
                        full.names = TRUE)
sizes <- do.call(rbind, lapply(coh_files, function(f) {
  cid <- tools::file_path_sans_ext(basename(f))
  n <- tryCatch(nrow(cohort_samples(cid, FREEZE, man, clinical = NULL)),
                error = function(e) NA_integer_)
  ns <- tryCatch(length(unique(cohort_samples(cid, FREEZE, man,
                                              clinical = NULL)$Subject_ID)),
                 error = function(e) NA_integer_)
  ch <- yaml::read_yaml(f)
  data.frame(cohort = cid, layer = ch$layer %||% NA,
             n_libraries = n, n_subjects = ns,
             needs_clinical = isTRUE(ch$requires_clinical),
             stringsAsFactors = FALSE)
}))
save_table(sizes, run, "cohort_sizes")
print(sizes, row.names = FALSE)

# ---- design summary -------------------------------------------

design <- table(lib$group, lib$arm, lib$timepoint, useNA = "ifany")
utils::write.csv(as.data.frame(design),
                 file.path(run$tables, "design_summary.csv"),
                 row.names = FALSE)

save_patient_table(lib, run, "freeze_libraries")   # has Subject_ID

status  <- if (length(issues)) "invalid" else "ok"
verdict <- if (length(issues)) {
  paste0(length(issues), " manifest issue(s): ",
         paste(substr(issues, 1, 60), collapse = " | "))
} else {
  sprintf("Manifest clean; freeze %s = %d libraries, %d subjects (%d HC)",
          FREEZE, nrow(lib), length(unique(lib$Subject_ID)), nrow(hc))
}

finalize_run(run, status = status, verdict = verdict)

if (length(issues)) {
  stop("Fix the manifest issues above before proceeding.", call. = FALSE)
}
