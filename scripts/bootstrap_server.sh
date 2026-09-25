#!/usr/bin/env bash
# =============================================================================
# bootstrap_server.sh
#
# Builds a self-contained R environment in $HOME on an Ubuntu server where you
# have NO root access. Designed for mininubio.ddns.med.umich.edu but written to
# be portable across the group's other servers.
#
# What it does, in order:
#   1. Preflight checks (compilers, headers, disk, load)
#   2. OpenBLAS -> ~/opt/openblas          (system has only reference BLAS)
#   3. R 4.6.1 from source -> ~/R/R-4.6.1  (system R is 4.3.2 / Bioc 3.18)
#   4. PATH + ~/.Renviron configuration
#   5. Bootstrap packages (renv, BiocManager, remotes) -> ~/R/library
#   6. SMOKE TEST: BPCells, then variancePartition  <-- risk front-loaded here
#
# It does NOT install the analysis stack. That happens per project via renv,
# using scripts/bootstrap_renv.R. Rationale in docs/ENVIRONMENT.md.
#
# Idempotent: safe to re-run. Completed stages are skipped unless --force.
#
# USAGE
#   chmod +x scripts/bootstrap_server.sh
#   screen -dmS bootstrap bash -lc 'scripts/bootstrap_server.sh 2>&1 | tee ~/bootstrap.log'
#   tail -f ~/bootstrap.log
#
# OPTIONS
#   --skip-openblas     use the system reference BLAS instead
#   --skip-r            only do config + packages (R already built)
#   --jobs N            parallel make jobs (default 8)
#   --force             rebuild even if a stage looks complete
#   --dry-run           print the plan and exit
# =============================================================================

set -euo pipefail

# ----- configuration ---------------------------------------------------------
R_VERSION="4.6.1"                  # matches the MacBook, library-compatible with 4.6.0
R_PREFIX="$HOME/R/R-${R_VERSION}"
OPENBLAS_PREFIX="$HOME/opt/openblas"
BOOTSTRAP_LIB="$HOME/R/library"
RENV_CACHE="$HOME/.cache/R/renv"
SRC_DIR="$HOME/src"
CRAN_SNAPSHOT="https://packagemanager.posit.co/cran/__linux__/focal/2026-06-01"
BPCELLS_UNIVERSE="https://bnprks.r-universe.dev"
BPCELLS_SHA="adc4a3c30f60a03522f58947d733d7d77a6eb2cf"   # recorded for provenance
RUNTIME_THREADS=4                  # conservative on a shared 32-core box

JOBS=8
SKIP_OPENBLAS=0
SKIP_R=0
FORCE=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-openblas) SKIP_OPENBLAS=1; shift ;;
    --skip-r)        SKIP_R=1; shift ;;
    --jobs)          JOBS="$2"; shift 2 ;;
    --force)         FORCE=1; shift ;;
    --dry-run)       DRY_RUN=1; shift ;;
    -h|--help)       sed -n '3,40p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# ----- output helpers --------------------------------------------------------
BOLD=$'\033[1m'; RESET=$'\033[0m'; RED=$'\033[31m'; GRN=$'\033[32m'; YLW=$'\033[33m'
step() { echo; echo "${BOLD}==> $*${RESET}"; }
ok()   { echo "  ${GRN}[ ok ]${RESET} $*"; }
warn() { echo "  ${YLW}[warn]${RESET} $*"; }
die()  { echo "  ${RED}[FAIL]${RESET} $*" >&2; exit 1; }
note() { echo "         $*"; }

START_TS=$(date +%s)
echo "${BOLD}ASSET server bootstrap${RESET}"
echo "  host        $(hostname)"
echo "  user        $USER"
echo "  date        $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "  R version   ${R_VERSION} -> ${R_PREFIX}"
echo "  OpenBLAS    $( [[ $SKIP_OPENBLAS -eq 1 ]] && echo 'SKIPPED (system BLAS)' || echo "${OPENBLAS_PREFIX}" )"
echo "  make jobs   ${JOBS}"
echo "  CRAN        ${CRAN_SNAPSHOT}"

if [[ $DRY_RUN -eq 1 ]]; then
  echo; echo "Dry run. Nothing was changed."; exit 0
fi

# =============================================================================
step "1. Preflight"
# =============================================================================

[[ "$(uname -s)" == "Linux" ]] || die "This script is for Linux. Use bootstrap_laptop.sh on macOS."

for c in gcc g++ gfortran make curl tar sed awk; do
  command -v "$c" >/dev/null || die "missing required tool: $c"
done
ok "compilers and tools present ($(gcc -dumpversion))"

# gcc 9 caps us at C++17. Documented limitation, not a blocker for the current
# stack (BPCells needs C++17 and works). See docs/ENVIRONMENT.md.
GCC_MAJOR=$(gcc -dumpversion | cut -d. -f1)
if [[ "$GCC_MAJOR" -lt 10 ]]; then
  warn "gcc ${GCC_MAJOR} -> C++17 maximum, no C++20"
  note "Fine for Seurat/BPCells/Bioc 3.23. If a future package needs C++20,"
  note "see the 'Toolchain ceiling' section of docs/ENVIRONMENT.md."
fi

MISSING_HDR=()
for h in zlib.h bzlib.h lzma.h pcre2.h curl/curl.h readline/readline.h; do
  find /usr/include -name "$(basename "$h")" 2>/dev/null | grep -q . || MISSING_HDR+=("$h")
done
if [[ ${#MISSING_HDR[@]} -gt 0 ]]; then
  die "missing R build headers: ${MISSING_HDR[*]} (needs root to install)"
fi
ok "R build headers present"

AVAIL_GB=$(df -BG --output=avail "$HOME" 2>/dev/null | tail -1 | tr -dc '0-9')
[[ -n "$AVAIL_GB" && "$AVAIL_GB" -ge 20 ]] || warn "less than 20 GB free in \$HOME (have ${AVAIL_GB:-?} GB)"
ok "disk: ${AVAIL_GB:-?} GB available"

NCPU=$(nproc)
LOAD=$(awk '{printf "%.1f", $1}' /proc/loadavg)
ok "cpus: ${NCPU}, load: ${LOAD}"
if awk "BEGIN{exit !($LOAD > $NCPU * 0.75)}"; then
  warn "machine is heavily loaded. Consider running this later, or --jobs 4."
fi

mkdir -p "$SRC_DIR" "$BOOTSTRAP_LIB" "$RENV_CACHE"

# =============================================================================
step "2. OpenBLAS"
# =============================================================================

if [[ $SKIP_OPENBLAS -eq 1 ]]; then
  warn "skipped by request; R will use the system reference BLAS"
  BLAS_CONFIGURE=(--with-blas --with-lapack)
  BLAS_LDFLAGS=""
elif [[ -f "$OPENBLAS_PREFIX/lib/libopenblas.so" && $FORCE -eq 0 ]]; then
  ok "already built at $OPENBLAS_PREFIX"
  BLAS_CONFIGURE=(--with-blas="-L$OPENBLAS_PREFIX/lib -lopenblas" --with-lapack)
  BLAS_LDFLAGS="-Wl,-rpath,$OPENBLAS_PREFIX/lib"
else
  # Resolve the newest release; fall back to a known-good tag if the API is
  # unreachable. The org moved from xianyi/ to OpenMathLib/.
  OB_VER="$(curl -fsSL https://api.github.com/repos/OpenMathLib/OpenBLAS/releases/latest 2>/dev/null \
            | sed -n 's/.*"tag_name": *"v\([0-9.]*\)".*/\1/p' | head -1)"
  if [[ -z "$OB_VER" ]]; then
    OB_VER="0.3.28"
    warn "GitHub API unreachable; falling back to OpenBLAS ${OB_VER}"
  fi
  note "building OpenBLAS ${OB_VER} (10-20 min)"

  cd "$SRC_DIR"
  OB_TAR="OpenBLAS-${OB_VER}.tar.gz"
  [[ -f "$OB_TAR" ]] || curl -fsSL -o "$OB_TAR" \
    "https://github.com/OpenMathLib/OpenBLAS/releases/download/v${OB_VER}/${OB_TAR}" \
    || die "could not download OpenBLAS ${OB_VER}"
  rm -rf "OpenBLAS-${OB_VER}"
  tar xzf "$OB_TAR"
  cd "OpenBLAS-${OB_VER}"

  # DYNAMIC_ARCH=1  -> one binary that picks the right kernel at runtime, so
  #                    this install is portable to the group's other servers.
  # USE_OPENMP=0    -> pthread threading. OpenMP-threaded BLAS fights with
  #                    R packages that use OpenMP themselves (data.table,
  #                    glmGamPoi). Threads are capped at runtime instead.
  # NUM_THREADS=64  -> compile-time ceiling, not the runtime default.
  make -j"$JOBS" DYNAMIC_ARCH=1 USE_OPENMP=0 NUM_THREADS=64 \
    >"$SRC_DIR/openblas_build.log" 2>&1 || die "OpenBLAS build failed; see $SRC_DIR/openblas_build.log"
  make PREFIX="$OPENBLAS_PREFIX" install >>"$SRC_DIR/openblas_build.log" 2>&1 \
    || die "OpenBLAS install failed"

  [[ -f "$OPENBLAS_PREFIX/lib/libopenblas.so" ]] || die "libopenblas.so not found after install"
  ok "OpenBLAS ${OB_VER} -> $OPENBLAS_PREFIX"
  BLAS_CONFIGURE=(--with-blas="-L$OPENBLAS_PREFIX/lib -lopenblas" --with-lapack)
  BLAS_LDFLAGS="-Wl,-rpath,$OPENBLAS_PREFIX/lib"
fi

# =============================================================================
step "3. R ${R_VERSION} from source"
# =============================================================================

if [[ $SKIP_R -eq 1 ]]; then
  warn "skipped by request"
elif [[ -x "$R_PREFIX/bin/R" && $FORCE -eq 0 ]]; then
  ok "already installed: $("$R_PREFIX/bin/R" --version | head -1)"
else
  cd "$SRC_DIR"
  R_TAR="R-${R_VERSION}.tar.gz"
  R_MAJOR="${R_VERSION%%.*}"
  if [[ ! -f "$R_TAR" ]]; then
    note "downloading R ${R_VERSION} source"
    curl -fsSL -o "$R_TAR" \
      "https://cran.r-project.org/src/base/R-${R_MAJOR}/${R_TAR}" \
      || die "could not download R ${R_VERSION}"
  fi

  # Record the hash for provenance. CRAN's per-directory MD5SUM file is not
  # reliably present, so verify against it when available and otherwise just
  # log what we got. The hash lands in docs/ENVIRONMENT.md.
  R_SHA256=$(sha256sum "$R_TAR" | awk '{print $1}')
  note "sha256: $R_SHA256"
  if curl -fsSL "https://cran.r-project.org/src/base/R-${R_MAJOR}/MD5SUM" -o /tmp/cran_md5 2>/dev/null; then
    EXPECT=$(grep -m1 " ${R_TAR}\$" /tmp/cran_md5 | awk '{print $1}' || true)
    GOT=$(md5sum "$R_TAR" | awk '{print $1}')
    if [[ -n "$EXPECT" ]]; then
      [[ "$EXPECT" == "$GOT" ]] && ok "md5 verified against CRAN" \
        || die "MD5 MISMATCH: expected $EXPECT got $GOT — delete $SRC_DIR/$R_TAR and retry"
    else
      warn "no CRAN md5 entry for ${R_TAR}; hash logged but unverified"
    fi
  else
    warn "CRAN MD5SUM unavailable; hash logged but unverified"
  fi

  rm -rf "R-${R_VERSION}"
  tar xzf "$R_TAR"
  cd "R-${R_VERSION}"

  # --enable-R-shlib       required by Positron's R backend (libR.so)
  # --without-x            headless server; cairo still provides png/pdf
  # --with-cairo           good raster output without X11
  # --enable-memory-profiling  Rprofmem(), useful on 1.3M-cell objects
  note "configure"
  LDFLAGS="$BLAS_LDFLAGS" ./configure \
    --prefix="$R_PREFIX" \
    --enable-R-shlib \
    --enable-memory-profiling \
    --without-x \
    --with-cairo \
    --with-libpng \
    --with-jpeglib \
    --with-libtiff \
    --with-recommended-packages \
    "${BLAS_CONFIGURE[@]}" \
    >"$SRC_DIR/r_configure.log" 2>&1 || die "configure failed; see $SRC_DIR/r_configure.log"

  note "make -j${JOBS}  (20-40 min)"
  make -j"$JOBS" >"$SRC_DIR/r_make.log" 2>&1 || die "make failed; see $SRC_DIR/r_make.log"

  note "make install"
  make install >>"$SRC_DIR/r_make.log" 2>&1 || die "make install failed"

  [[ -x "$R_PREFIX/bin/R" ]] || die "R binary missing after install"
  ok "$("$R_PREFIX/bin/R" --version | head -1)"

  # Confirm the BLAS that actually got linked, not the one we asked for.
  BLAS_LINKED=$("$R_PREFIX/bin/Rscript" -e 'cat(sessionInfo()$BLAS)' 2>/dev/null || echo unknown)
  ok "BLAS linked: $BLAS_LINKED"
  if [[ $SKIP_OPENBLAS -eq 0 && "$BLAS_LINKED" != *openblas* ]]; then
    warn "expected OpenBLAS but got: $BLAS_LINKED — check $SRC_DIR/r_configure.log"
  fi
fi

# =============================================================================
step "4. Shell and R configuration"
# =============================================================================

# PATH: PREPEND, unlike the append pattern used elsewhere on this box.
# Appending leaves /usr/bin/R (4.3.2) winning for a bare `R`. Prepending makes
# 4.6.1 the default, and an active conda env still overrides because
# `conda activate` prepends later in the shell's lifetime.
MARK="# >>> ASSET R environment >>>"
if grep -qF "$MARK" "$HOME/.bashrc" 2>/dev/null && [[ $FORCE -eq 0 ]]; then
  ok ".bashrc block already present"
else
  # remove any previous block before rewriting
  if grep -qF "$MARK" "$HOME/.bashrc" 2>/dev/null; then
    sed -i "/$MARK/,/# <<< ASSET R environment <<</d" "$HOME/.bashrc"
  fi
  cat >> "$HOME/.bashrc" <<EOF

$MARK
# Managed by scripts/bootstrap_server.sh — edit there, not here.
# Prepended so a bare \`R\` is ${R_VERSION}, not the system 4.3.2.
# An active conda env still wins, because conda prepends after this runs.
export PATH="$R_PREFIX/bin:\$PATH"
export LD_LIBRARY_PATH="$OPENBLAS_PREFIX/lib:\${LD_LIBRARY_PATH:-}"
# <<< ASSET R environment <<<
EOF
  ok "PATH block written to ~/.bashrc"
fi

# ~/.Renviron: thread caps and the shared renv cache.
if [[ -f "$HOME/.Renviron" && $FORCE -eq 0 ]] && grep -q "RENV_PATHS_CACHE" "$HOME/.Renviron"; then
  ok "~/.Renviron already configured"
else
  [[ -f "$HOME/.Renviron" ]] && cp "$HOME/.Renviron" "$HOME/.Renviron.bak.$(date +%Y%m%d%H%M%S)"
  cat > "$HOME/.Renviron" <<EOF
# Managed by scripts/bootstrap_server.sh

# Conservative thread caps. This is a shared box with no scheduler; a single
# R process taking all 32 cores harms everyone. Override per job:
#   OMP_NUM_THREADS=16 jobs/run.sh <script>
OPENBLAS_NUM_THREADS=${RUNTIME_THREADS}
OMP_NUM_THREADS=${RUNTIME_THREADS}
MKL_NUM_THREADS=${RUNTIME_THREADS}

# Shared renv cache: packages are hard-linked into each project's library, so
# a second project costs almost no extra disk and no recompilation.
RENV_PATHS_CACHE=${RENV_CACHE}

# Bootstrap library. Project packages live in <project>/renv/library, NOT here.
R_LIBS_USER=${BOOTSTRAP_LIB}
EOF
  ok "~/.Renviron written (threads=${RUNTIME_THREADS}, renv cache shared)"
fi

export PATH="$R_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$OPENBLAS_PREFIX/lib:${LD_LIBRARY_PATH:-}"
RSCRIPT="$R_PREFIX/bin/Rscript"

# =============================================================================
step "5. Bootstrap packages"
# =============================================================================

# Deliberately minimal: only what is needed to drive renv. The analysis stack
# lives in per-project renv libraries.
"$RSCRIPT" -e "
options(repos = c(CRAN = '${CRAN_SNAPSHOT}'),
        Ncpus = ${JOBS}, warn = 1)
lib <- '${BOOTSTRAP_LIB}'
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
need <- setdiff(c('renv','BiocManager','remotes'),
                rownames(installed.packages(lib.loc = lib)))
if (length(need)) install.packages(need, lib = lib) else cat('already installed\n')
for (p in c('renv','BiocManager','remotes'))
  cat(sprintf('  %-14s %s\n', p, as.character(packageVersion(p, lib.loc = lib))))
cat('Bioconductor for this R:', as.character(BiocManager::version()), '\n')
" || die "bootstrap package installation failed"
ok "renv, BiocManager, remotes installed"

# =============================================================================
step "6. Smoke test (risk front-loaded)"
# =============================================================================

# These two are installed FIRST, before any project renv work, because they are
# the two most likely to fail on this machine:
#
#   BPCells           needs C++17 + HDF5. Her build proves it works here, but
#                     ours must link against the same system HDF5 1.10.4.
#   variancePartition never installed on this box by anyone. It is the engine
#                     for the longitudinal pseudobulk contrasts (repeated
#                     measures across M0/M3/M6), so a failure changes the
#                     analysis plan, not just the environment.
#
# Fifteen minutes here beats discovering a problem in hour three of the full
# install. If variancePartition fails, the fallback is
# limma::duplicateCorrelation() — workable, less good. Recorded either way.

SMOKE_LIB="$HOME/R/smoke-test-lib"
mkdir -p "$SMOKE_LIB"

echo; note "6a. BPCells (C++17 + HDF5)"
"$RSCRIPT" -e "
options(repos = c(BPCELLS = '${BPCELLS_UNIVERSE}', CRAN = '${CRAN_SNAPSHOT}'),
        Ncpus = ${JOBS}, warn = 1)
.libPaths(c('${SMOKE_LIB}', .libPaths()))
if (!requireNamespace('BPCells', quietly = TRUE))
  install.packages('BPCells', lib = '${SMOKE_LIB}')
library(BPCells)
cat('  BPCells', as.character(packageVersion('BPCells')), '\n')
set.seed(1)
m <- matrix(rpois(1000, 3), nrow = 100,
            dimnames = list(paste0('g', 1:100), paste0('c', 1:10)))
d <- file.path(tempdir(), 'bpcells_smoke')
if (dir.exists(d)) unlink(d, recursive = TRUE)
write_matrix_dir(as(m, 'dgCMatrix'), dir = d)
back <- open_matrix_dir(d)
stopifnot(identical(dim(back), dim(m)))
cat('  on-disk write + read: OK\n')
unlink(d, recursive = TRUE)
" && ok "BPCells works" || {
  warn "BPCells FAILED — this blocks the single-cell workflow"
  note "check HDF5: pkg-config --modversion hdf5  (expect 1.10.4)"
  note "compare with the working build: ldd /home/jarnagin/R/library/BPCells/libs/BPCells.so"
  BPCELLS_OK=0
}

echo; note "6b. variancePartition / dream (no precedent on this machine)"
"$RSCRIPT" -e "
options(repos = c(CRAN = '${CRAN_SNAPSHOT}'), Ncpus = ${JOBS}, warn = 1)
.libPaths(c('${SMOKE_LIB}', .libPaths()))
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager', lib = '${SMOKE_LIB}')
if (!requireNamespace('variancePartition', quietly = TRUE))
  BiocManager::install('variancePartition', lib = '${SMOKE_LIB}',
                       ask = FALSE, update = FALSE)
suppressPackageStartupMessages(library(variancePartition))
cat('  variancePartition', as.character(packageVersion('variancePartition')), '\n')
# Minimal mixed-model fit mirroring the real design: 12 subjects x 2 timepoints
set.seed(1)
n_sub <- 12; tp <- c('M00','M03')
meta <- expand.grid(subject = factor(paste0('s', 1:n_sub)), timepoint = factor(tp))
cnt <- matrix(rnbinom(200 * nrow(meta), mu = 50, size = 5), nrow = 200,
              dimnames = list(paste0('g', 1:200), paste0('u', 1:nrow(meta))))
suppressPackageStartupMessages({library(edgeR)})
dge <- calcNormFactors(DGEList(cnt))
form <- ~ timepoint + (1 | subject)
vobj <- voomWithDreamWeights(dge, form, meta, quiet = TRUE)
fit  <- dream(vobj, form, meta, quiet = TRUE)
fit  <- eBayes(fit)
res  <- topTable(fit, coef = 'timepointM03', number = 3)
cat('  dream() fit on a repeated-measures design: OK\n')
cat('  (', nrow(res), 'rows returned )\n')
" && ok "variancePartition / dream works" || {
  warn "variancePartition FAILED"
  note "Fallback: limma::duplicateCorrelation() for repeated measures."
  note "Record the decision in docs/decisions.md and set"
  note "engine: limma_voom in config/de_runs/*.yml"
}

# =============================================================================
step "Done"
# =============================================================================

ELAPSED=$(( ($(date +%s) - START_TS) / 60 ))
cat <<EOF

  elapsed         ${ELAPSED} min
  R               $R_PREFIX/bin/R
  OpenBLAS        $( [[ $SKIP_OPENBLAS -eq 1 ]] && echo '(system BLAS)' || echo "$OPENBLAS_PREFIX" )
  bootstrap lib   $BOOTSTRAP_LIB
  renv cache      $RENV_CACHE
  smoke-test lib  $SMOKE_LIB   (safe to delete: rm -rf $SMOKE_LIB)
  logs            $SRC_DIR/{openblas_build,r_configure,r_make}.log

  NEXT
  1. Open a new shell so the PATH change takes effect:
         exec bash -l
     Then confirm:
         which R && R --version | head -1
     Expect $R_PREFIX/bin/R and ${R_VERSION}.

  2. Set up the project library (~2-3 hours, everything compiles from source
     because Posit has no focal binaries — see docs/ENVIRONMENT.md):
         cd ~/projects/asset-skin-natural-history
         screen -dmS renv bash -lc 'Rscript scripts/bootstrap_renv.R 2>&1 | tee ~/renv_install.log'
         tail -f ~/renv_install.log

  3. Commit the lockfile when it finishes:
         git add renv.lock && git commit -m 'renv: initial snapshot (freeze01)' && git push

  Record anything surprising in docs/decisions.md while it is fresh.
EOF
