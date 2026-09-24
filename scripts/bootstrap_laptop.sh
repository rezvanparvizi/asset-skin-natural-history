#!/usr/bin/env bash
# =============================================================================
# bootstrap_laptop.sh   (macOS, Apple Silicon)
#
# Sets up the U-M MacBook for "level (b)" development: enough of the stack to
# write and debug functions locally against small test objects, then run them
# unchanged on the server. Not for real analysis — the data lives on the server
# and a 1.3M-cell object will not fit here.
#
# What it does:
#   1. Preflight (R version, Xcode tools, gfortran)
#   2. Homebrew, if absent
#   3. HDF5 + friends via brew  (BPCells needs HDF5)
#   4. ~/.Renviron with the same thread caps and shared renv cache pattern
#   5. Dev + single-cell package set, pinned to the SAME CRAN snapshot as the
#      server so a function written here behaves identically there
#
# ARCHITECTURE NOTE
#   Laptop is arm64, server is x86_64. renv.lock is portable between them;
#   renv/library/ is NOT — compiled packages can never be copied across.
#   Parity means "same versions, independently built".
#
# USAGE
#   chmod +x scripts/bootstrap_laptop.sh
#   ./scripts/bootstrap_laptop.sh
#
# OPTIONS
#   --minimal     editing tools only (level a): yaml, jsonlite, languageserver
#   --skip-brew   assume HDF5 is already present
# =============================================================================

set -euo pipefail

R_EXPECTED="4.6.1"
CRAN_SNAPSHOT="https://packagemanager.posit.co/cran/2026-06-01"
RENV_CACHE="$HOME/Library/Caches/org.R-project.R/R/renv"
RUNTIME_THREADS=4

MINIMAL=0
SKIP_BREW=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --minimal)   MINIMAL=1; shift ;;
    --skip-brew) SKIP_BREW=1; shift ;;
    -h|--help)   sed -n '3,30p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

BOLD=$'\033[1m'; RESET=$'\033[0m'; RED=$'\033[31m'; GRN=$'\033[32m'; YLW=$'\033[33m'
step() { echo; echo "${BOLD}==> $*${RESET}"; }
ok()   { echo "  ${GRN}[ ok ]${RESET} $*"; }
warn() { echo "  ${YLW}[warn]${RESET} $*"; }
die()  { echo "  ${RED}[FAIL]${RESET} $*" >&2; exit 1; }
note() { echo "         $*"; }

echo "${BOLD}ASSET laptop bootstrap${RESET}"
echo "  host      $(hostname)"
echo "  arch      $(uname -m)"
echo "  macOS     $(sw_vers -productVersion 2>/dev/null || echo '?')"
echo "  mode      $( [[ $MINIMAL -eq 1 ]] && echo 'minimal (level a)' || echo 'development (level b)' )"

# =============================================================================
step "1. Preflight"
# =============================================================================

[[ "$(uname -s)" == "Darwin" ]] || die "macOS only. Use bootstrap_server.sh on Linux."
[[ "$(uname -m)" == "arm64" ]] || warn "not arm64 — this script assumes Apple Silicon"

command -v R >/dev/null || die "R not found. Install R-${R_EXPECTED}-arm64.pkg from https://cran.r-project.org/bin/macosx/"
R_VER=$(R --version | head -1 | sed 's/R version \([0-9.]*\).*/\1/')
R_ARCH=$(R -q -e 'cat(R.version$arch)' --vanilla 2>/dev/null | tail -1)
ok "R ${R_VER} (${R_ARCH})"
[[ "$R_ARCH" == "aarch64" ]] || warn "R is ${R_ARCH}, not aarch64 — you may have the Intel build"
if [[ "$R_VER" != "$R_EXPECTED" ]]; then
  warn "R ${R_VER} here vs ${R_EXPECTED} expected on the server"
  note "R keys its library to major.minor, so 4.6.x are mutually compatible."
  note "A different MAJOR.MINOR breaks renv.lock parity — keep both machines"
  note "on 4.6.x and upgrade them together, deliberately."
fi

xcode-select -p >/dev/null 2>&1 || die "Xcode CLT missing: xcode-select --install"
ok "Xcode command line tools"

FC=$(R CMD config FC 2>/dev/null || echo "")
[[ -n "$FC" ]] && ok "Fortran: $FC" || warn "no Fortran compiler; install gfortran from https://mac.r-project.org/tools/"

# =============================================================================
step "2. Homebrew"
# =============================================================================

if [[ $SKIP_BREW -eq 1 ]]; then
  warn "skipped by request"
elif command -v brew >/dev/null 2>&1; then
  ok "already installed ($(brew --version | head -1))"
else
  note "installing Homebrew (will prompt for your Mac password)"
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" \
    || die "Homebrew install failed"
  # The installer adds the shellenv line to ~/.zprofile itself; add it only if
  # it is genuinely absent (an earlier manual attempt may have been cleaned up).
  if ! grep -q 'brew shellenv' "$HOME/.zprofile" 2>/dev/null; then
    echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> "$HOME/.zprofile"
  fi
  eval "$(/opt/homebrew/bin/brew shellenv)"
  ok "Homebrew installed"
fi

# =============================================================================
step "3. System libraries"
# =============================================================================

if [[ $SKIP_BREW -eq 1 || $MINIMAL -eq 1 ]]; then
  warn "skipped"
else
  # hdf5 is the one that matters — BPCells links against it. The rest cover
  # common compile failures for spatial/plotting packages.
  BREW_PKGS=(hdf5 gsl glpk libgit2 pkg-config cmake)
  for p in "${BREW_PKGS[@]}"; do
    if brew list --formula "$p" >/dev/null 2>&1; then
      ok "$p (already)"
    else
      note "installing $p"
      brew install "$p" >/dev/null || warn "$p failed — continuing"
    fi
  done
  HDF5_PREFIX="$(brew --prefix hdf5 2>/dev/null || true)"
  [[ -n "$HDF5_PREFIX" ]] && ok "hdf5 at $HDF5_PREFIX"
fi

# =============================================================================
step "4. ~/.Renviron"
# =============================================================================

if [[ -f "$HOME/.Renviron" ]] && grep -q RENV_PATHS_CACHE "$HOME/.Renviron"; then
  ok "already configured"
else
  [[ -f "$HOME/.Renviron" ]] && cp "$HOME/.Renviron" "$HOME/.Renviron.bak.$(date +%Y%m%d%H%M%S)"
  mkdir -p "$RENV_CACHE"
  cat > "$HOME/.Renviron" <<EOF
# Managed by scripts/bootstrap_laptop.sh

# Same conservative thread caps as the server, so timing behaviour is
# comparable when you test a function in both places.
OPENBLAS_NUM_THREADS=${RUNTIME_THREADS}
OMP_NUM_THREADS=${RUNTIME_THREADS}
VECLIB_MAXIMUM_THREADS=${RUNTIME_THREADS}

# Shared renv cache across local projects.
RENV_PATHS_CACHE=${RENV_CACHE}
EOF
  ok "written (threads=${RUNTIME_THREADS})"
fi

# Help compilers find Homebrew headers when building from source.
if [[ $MINIMAL -eq 0 ]] && command -v brew >/dev/null 2>&1; then
  BP="$(brew --prefix)"
  mkdir -p "$HOME/.R"
  if ! grep -q "ASSET" "$HOME/.R/Makevars" 2>/dev/null; then
    [[ -f "$HOME/.R/Makevars" ]] && cp "$HOME/.R/Makevars" "$HOME/.R/Makevars.bak.$(date +%Y%m%d%H%M%S)"
    cat >> "$HOME/.R/Makevars" <<EOF

# ASSET: Homebrew include/lib paths so source installs find hdf5, gsl, glpk
CPPFLAGS += -I${BP}/include
LDFLAGS  += -L${BP}/lib
EOF
    ok "~/.R/Makevars updated for Homebrew paths"
  else
    ok "~/.R/Makevars already configured"
  fi
fi

# =============================================================================
step "5. R packages"
# =============================================================================

# macOS note: the P3M snapshot serves SOURCE for macOS (their binaries are
# Linux-only), so these compile. That is deliberate: matching the server's
# versions matters more here than install speed, and this is a one-time cost
# on a machine that only ever holds small test objects.

if [[ $MINIMAL -eq 1 ]]; then
  Rscript --vanilla -e "
options(repos = c(CRAN = '${CRAN_SNAPSHOT}'), Ncpus = 4, warn = 1)
pkgs <- c('yaml','jsonlite','languageserver','renv','ggplot2','dplyr','data.table')
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need) else cat('all present\n')
for (p in pkgs) cat(sprintf('  %-16s %s\n', p, tryCatch(as.character(packageVersion(p)), error=function(e) 'MISSING')))
" || die "package install failed"
  ok "minimal (level a) set installed"
else
  note "development set — expect 30-60 min, compiling from source"
  Rscript --vanilla -e "
options(repos = c(CRAN = '${CRAN_SNAPSHOT}',
                  BPCELLS = 'https://bnprks.r-universe.dev'),
        Ncpus = 4, warn = 1)

# editing and infrastructure
dev <- c('languageserver','renv','BiocManager','remotes','yaml','jsonlite',
         'devtools','testthat','usethis')
# analysis packages needed to write and debug real functions locally
sc  <- c('Seurat','SeuratObject','Matrix','harmony','data.table','dplyr',
         'tidyr','ggplot2','patchwork','ggrepel','RColorBrewer','viridis')

need <- setdiff(c(dev, sc), rownames(installed.packages()))
if (length(need)) install.packages(need) else cat('CRAN set already present\n')

# BPCells — the local install most likely to fail (needs Homebrew hdf5)
if (!requireNamespace('BPCells', quietly = TRUE)) {
  cat('\ninstalling BPCells\n')
  tryCatch(install.packages('BPCells'),
           error = function(e) message('BPCells failed: ', conditionMessage(e)))
}

cat('\nInstalled:\n')
for (p in c(dev, sc, 'BPCells'))
  cat(sprintf('  %-16s %s\n', p,
      tryCatch(as.character(packageVersion(p)), error = function(e) 'MISSING')))
" || warn "some packages failed — see output above"
  ok "development (level b) set attempted"
fi

# =============================================================================
step "Done"
# =============================================================================

cat <<EOF

  NEXT

  1. Restart Positron so it picks up ~/.Renviron and any new packages.

  2. Verify the local stack matches the server's versions:
         Rscript -e 'for (p in c("Seurat","SeuratObject","Matrix","harmony","BPCells")) cat(sprintf("%-14s %s\\n", p, tryCatch(as.character(packageVersion(p)), error=function(e) "MISSING")))'
     Expect Seurat 5.5.0, SeuratObject 5.4.0, Matrix 1.7-5, harmony 2.0.3.

  3. Make a small test object on the SERVER and copy it here, so you can
     develop against real structures:
         # on mininubio, once you have the real object
         # obj_small <- obj[, sample(colnames(obj), 20000)]
         # saveRDS(obj_small, "~/asset-data/objects/test_20k.rds")
         # then, on the laptop:
         scp mininubio:~/asset-data/objects/test_20k.rds ~/asset-test/

  4. Keep both machines on R 4.6.x and upgrade them together, on purpose.
     If macOS R jumps to 4.7 while the server stays on 4.6, renv.lock stops
     restoring and parity is silently gone.

  WHAT THIS MACHINE IS FOR
     Writing and debugging functions in R/ against small test objects, editing
     scripts, git. NOT real analysis — the data stays on the server, and
     anything heavy goes through jobs/run.sh in a screen session there.
EOF
