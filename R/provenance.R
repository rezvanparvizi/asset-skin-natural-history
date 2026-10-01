# ==============================================================
# R/provenance.R
#
# The run contract. Every analysis script calls init_run() first and
# finalize_run() last. Between them, the run directory is created and
# is self-describing: it carries its own resolved parameters, package
# versions, and code version.
#
# This is the single highest-value convention in the repository. It is
# what lets you answer, two years from now, "what exactly produced
# this figure?" without archaeology.
#
# Requires: yaml
# ==============================================================

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required. install.packages('yaml')", call. = FALSE)
}

# ---- internal helpers ----------------------------------------

.git_sha <- function(repo = REPO) {
  out <- suppressWarnings(tryCatch(
    system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_))
  if (length(out) == 0 || is.na(out[1])) return(NA_character_)
  out[1]
}

.git_dirty <- function(repo = REPO) {
  out <- suppressWarnings(tryCatch(
    system2("git", c("-C", shQuote(repo), "status", "--porcelain"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_))
  if (length(out) == 0) return(FALSE)
  if (length(out) == 1 && is.na(out[1])) return(NA)
  # The run registry is written by every run; it does not describe code.
  out <- out[!grepl("docs/runs\\.csv$", out)]
  length(out) > 0
}

.git_branch <- function(repo = REPO) {
  out <- suppressWarnings(tryCatch(
    system2("git", c("-C", shQuote(repo), "rev-parse", "--abbrev-ref", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_))
  if (length(out) == 0) NA_character_ else out[1]
}

.calling_script <- function() {
  # Rscript path, when run non-interactively
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(sub("^--file=", "", f[1]))
  # sys.frames fallback for source()
  for (i in rev(seq_len(sys.nframe()))) {
    fn <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(fn)) return(fn)
  }
  "interactive"
}

.rel_to_repo <- function(p, repo = REPO) {
  p <- tryCatch(normalizePath(p, mustWork = FALSE), error = function(e) p)
  if (startsWith(p, repo)) sub(paste0("^", repo, "/?"), "", p) else p
}

.next_run_id <- function(registry) {
  if (!file.exists(registry)) return("R0001")
  reg <- utils::read.csv(registry, stringsAsFactors = FALSE,
                         colClasses = "character")
  if (nrow(reg) == 0) return("R0001")
  nums <- suppressWarnings(as.integer(sub("^R", "", reg$run_id)))
  nums <- nums[!is.na(nums)]
  sprintf("R%04d", if (length(nums)) max(nums) + 1L else 1L)
}

REGISTRY_COLS <- c("run_id", "date", "freeze", "cohort", "labelset",
                   "stage", "run_name", "script", "git_sha", "git_dirty",
                   "status", "verdict")

# Runs can start concurrently (several screen jobs), so a run_id is
# RESERVED at init_run() by appending a "running" row under a lock, and
# finalize_run() updates that row. Before 2026-10-01 the row was written
# only at the end, so two runs started together got the same run_id.
.with_registry_lock <- function(registry, expr, timeout = 120) {
  lock <- paste0(registry, ".lock")
  t0 <- Sys.time()
  while (!dir.create(lock, showWarnings = FALSE)) {
    if (difftime(Sys.time(), t0, units = "secs") > timeout) {
      stop("Could not lock ", registry, " (stale lock? remove ", lock, ")", call. = FALSE)
    }
    Sys.sleep(0.2)
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  force(expr)
}

.registry_row <- function(run_id, context, stage, run_name, code, status, verdict) {
  data.frame(
    run_id    = run_id,
    date      = format(Sys.Date(), "%Y-%m-%d"),
    freeze    = context$freeze,
    cohort    = context$scope,
    labelset  = context$labelset,
    stage     = stage,
    run_name  = run_name,
    script    = code$script,
    git_sha   = substr(ifelse(is.na(code$git_sha), "NA", code$git_sha), 1, 7),
    git_dirty = as.character(code$git_dirty),
    status    = status,
    verdict   = verdict,
    stringsAsFactors = FALSE
  )[, REGISTRY_COLS]
}

.ensure_registry <- function(registry) {
  if (!file.exists(registry)) {
    dir.create(dirname(registry), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      as.data.frame(setNames(rep(list(character()), length(REGISTRY_COLS)),
                             REGISTRY_COLS)),
      registry, sep = ",", row.names = FALSE, qmethod = "double")
  }
  invisible(registry)
}

# ---- init_run ------------------------------------------------

#' Start a run: create its directory and write its provenance
#'
#' @param stage Stage directory name, e.g. "14_sc_pseudobulk_de".
#' @param run_name Descriptive run name, e.g. "fib_improver_vs_non_m0".
#'   Describes the QUESTION, not the date or the parameters.
#' @param params Named list of analysis parameters to record verbatim
#'   (resolutions, thresholds, contrast specifications, covariates...).
#' @param notes One-line human description of intent.
#' @param freeze,cohort,labelset Run context. Default to the values in
#'   config/paths.R, which jobs/run.sh sets from the command line.
#' @param overwrite If FALSE (default; TRUE when ASSET_OVERWRITE=1) and the run directory already
#'   exists and is non-empty, stop. Prevents silently clobbering a
#'   previous result.
#' @param exploratory Mark this run as Layer 3. Forces the cohort
#'   directory to "exploratory" so it can never be mistaken for a
#'   confirmatory result.
#'
#' @return A list with: dir, objects, plots, tables, run_id, and the
#'   full context. Pass it to finalize_run().
init_run <- function(stage,
                     run_name,
                     params      = list(),
                     notes       = "",
                     freeze      = FREEZE,
                     cohort      = COHORT,
                     labelset    = LABELSET,
                     overwrite   = identical(Sys.getenv("ASSET_OVERWRITE"), "1"),
                     exploratory = FALSE) {

  stopifnot(is.character(stage), length(stage) == 1,
            is.character(run_name), length(run_name) == 1)

  bad <- grepl("[^A-Za-z0-9_]", run_name)
  if (bad) {
    stop("run_name must be snake_case alphanumeric only (got '", run_name,
         "'). No dates, no dots, no spaces, no 'v2'/'final'.", call. = FALSE)
  }
  if (grepl("(?i)(final|v[0-9]+|new|old|test|tmp|[0-9]{4}-[0-9]{2}-[0-9]{2})",
            run_name, perl = TRUE)) {
    warning("run_name '", run_name, "' contains a version/date token. ",
            "Those belong in run_config.yml and docs/runs.csv, not the name.",
            call. = FALSE)
  }

  scope <- if (exploratory) "exploratory" else cohort

  dir  <- file.path(RESULTS, freeze, scope, stage, run_name)
  subs <- c(objects = file.path(dir, "objects"),
            plots   = file.path(dir, "plots"),
            qc      = file.path(dir, "plots", "QC"),
            tables  = file.path(dir, "tables"))

  if (dir.exists(dir) && length(list.files(dir)) > 0 && !overwrite) {
    stop("Run directory already exists and is not empty:\n  ", dir,
         "\nEither choose a new run_name, or pass overwrite = TRUE if you ",
         "really mean to replace it (and mark the old run 'superseded' in ",
         "docs/runs.csv).", call. = FALSE)
  }

  for (p in c(dir, subs)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

  script <- .rel_to_repo(.calling_script())
  sha    <- .git_sha()
  dirty  <- .git_dirty()

  registry <- file.path(DOCS, "runs.csv")
  .ensure_registry(registry)
  run_id <- .with_registry_lock(registry, {
    id <- .next_run_id(registry)
    utils::write.table(
      .registry_row(id, list(freeze = freeze, scope = scope, labelset = labelset),
                    stage, run_name, list(script = script, git_sha = sha, git_dirty = dirty),
                    "running", ""),
      registry, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE,
      qmethod = "double")
    id
  })

  cfg <- list(
    run_id      = run_id,
    run_name    = run_name,
    stage       = stage,
    started     = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    notes       = notes,
    exploratory = exploratory,
    context = list(
      freeze   = freeze,
      cohort   = cohort,
      labelset = labelset,
      scope    = scope
    ),
    code = list(
      script     = script,
      git_branch = .git_branch(),
      git_sha    = sha,
      git_dirty  = dirty
    ),
    inputs = list(
      freeze_file = .rel_to_repo(freeze_file(freeze)),
      cohort_file = .rel_to_repo(cohort_file(cohort)),
      labels_file = .rel_to_repo(labels_file(freeze, labelset))
    ),
    params = params,
    system = list(
      user     = unname(Sys.info()[["user"]]),
      host     = unname(Sys.info()[["nodename"]]),
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      wd       = getwd()
    )
  )

  yaml::write_yaml(cfg, file.path(dir, "run_config.yml"))
  writeLines(capture.output(utils::sessionInfo()),
             file.path(dir, "session_info.txt"))
  writeLines(c(paste("sha   :", sha),
               paste("branch:", .git_branch()),
               paste("dirty :", dirty)),
             file.path(dir, "git_sha.txt"))

  if (isTRUE(dirty)) {
    warning("Uncommitted changes present. This run's git_sha does not fully ",
            "describe the code that produced it. Commit before a run you ",
            "intend to keep.", call. = FALSE)
  }

  message(sprintf("[%s] %s / %s / %s\n  -> %s",
                  run_id, freeze, scope, run_name, dir))

  structure(
    c(cfg, list(dir = dir, objects = subs[["objects"]],
                plots = subs[["plots"]], qc = subs[["qc"]],
                tables = subs[["tables"]])),
    class = "asset_run"
  )
}

# ---- finalize_run --------------------------------------------

#' Close a run: record status and a one-line verdict in the registry
#'
#' @param run The object returned by init_run().
#' @param status One of "ok", "failed", "invalid", "superseded",
#'   "exploratory".
#'   - ok:          completed, result is usable
#'   - failed:      crashed or did not finish
#'   - invalid:     completed, but the result is wrong (bug, wrong
#'                  filter, wrong contrast). KEEP the row. Say why.
#'   - superseded:  replaced by a later run; put that run_id in verdict
#'   - exploratory: hypothesis-generating only
#' @param verdict One sentence on what you concluded. This is the field
#'   you will actually grep in a year. Write it for your future self.
finalize_run <- function(run,
                         status  = c("ok", "failed", "invalid",
                                     "superseded", "exploratory"),
                         verdict = "") {

  status <- match.arg(status)
  if (!inherits(run, "asset_run")) {
    stop("finalize_run() expects the object returned by init_run().",
         call. = FALSE)
  }
  if (!nzchar(verdict)) {
    warning("Empty verdict. A row in runs.csv with no verdict is nearly ",
            "useless later. One sentence is enough.", call. = FALSE)
  }

  registry <- file.path(DOCS, "runs.csv")
  .ensure_registry(registry)

  row <- .registry_row(run$run_id, run$context, run$stage, run$run_name,
                       run$code, status, verdict)

  .with_registry_lock(registry, {
    reg <- utils::read.csv(registry, stringsAsFactors = FALSE, colClasses = "character")
    i <- which(reg$run_id == run$run_id & reg$status == "running" &
               reg$run_name == run$run_name)
    if (length(i) == 1) {
      reg[i, ] <- row
      utils::write.table(reg, registry, sep = ",", row.names = FALSE,
                         col.names = TRUE, qmethod = "double")
    } else {
      utils::write.table(row, registry, sep = ",", append = TRUE,
                         col.names = FALSE, row.names = FALSE, qmethod = "double")
    }
  })

  # also stamp the run directory so it is self-contained
  cfg <- yaml::read_yaml(file.path(run$dir, "run_config.yml"))
  cfg$finished <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  cfg$status   <- status
  cfg$verdict  <- verdict
  yaml::write_yaml(cfg, file.path(run$dir, "run_config.yml"))

  message(sprintf("[%s] %s — %s", run$run_id, status, verdict))
  invisible(row)
}

# ---- convenience ---------------------------------------------

#' Mark an earlier run as superseded (edits docs/runs.csv in place)
mark_superseded <- function(run_id, by_run_id, reason = "") {
  registry <- file.path(DOCS, "runs.csv")
  reg <- utils::read.csv(registry, stringsAsFactors = FALSE,
                         colClasses = "character")
  i <- which(reg$run_id == run_id)
  if (!length(i)) stop("run_id not found in registry: ", run_id, call. = FALSE)
  reg$status[i]  <- "superseded"
  reg$verdict[i] <- paste0("superseded by ", by_run_id,
                           if (nzchar(reason)) paste0(" — ", reason) else "")
  utils::write.csv(reg, registry, row.names = FALSE, qmethod = "double")
  message("Marked ", run_id, " superseded by ", by_run_id)
  invisible(reg[i, ])
}

#' Read the run registry as a data frame
read_runs <- function() {
  registry <- file.path(DOCS, "runs.csv")
  if (!file.exists(registry)) return(NULL)
  utils::read.csv(registry, stringsAsFactors = FALSE, colClasses = "character")
}

#' Which runs used a given labelset? (the question you will ask after
#' re-annotating)
runs_using <- function(labelset = NULL, freeze = NULL, status = NULL) {
  reg <- read_runs()
  if (is.null(reg)) return(NULL)
  if (!is.null(labelset)) reg <- reg[reg$labelset == labelset, ]
  if (!is.null(freeze))   reg <- reg[reg$freeze == freeze, ]
  if (!is.null(status))   reg <- reg[reg$status %in% status, ]
  reg
}
