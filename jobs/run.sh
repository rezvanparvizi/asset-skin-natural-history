#!/usr/bin/env bash
# ==============================================================
# jobs/run.sh — launch an R script in a DETACHED screen session
#
# This server has no scheduler, so `screen` is the job manager. The key
# trick is `screen -dmS`, which starts the session already detached:
# you never have to remember Ctrl-A D, and closing your laptop or
# dropping VPN cannot kill the job.
#
# Usage:
#   jobs/run.sh analysis/03_sc_integration/03_1_harmony.R
#   jobs/run.sh analysis/05_sc_subcluster_fibroblast/05_2_sweep_resolution.R \
#       --freeze freeze01 --cohort reference --labelset labelset01
#   jobs/run.sh pipelines/pseudobulk_de.R \
#       --cohort placebo --config config/de_runs/fib_improver_vs_non_m0.yml
#
# Options:
#   --freeze    FREEZE      default freeze01
#   --cohort    COHORT      default reference
#   --labelset  LABELSET    default labelset01
#   --config    FILE        passed through to the script as ASSET_CONFIG
#   --name      NAME        screen session name (default derived from script)
#   --fg                    run in the foreground instead (for quick tests)
#
# Then:
#   jobs/status.sh                 what is running, and recent log tails
#   screen -ls                     list sessions
#   screen -r <name>               attach (Ctrl-A then D to detach again)
#   tail -f logs/<logfile>         watch output without attaching
# ==============================================================

set -euo pipefail

# ---- locate repo root ---------------------------------------
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

# ---- defaults ------------------------------------------------
FREEZE="freeze01"
COHORT="reference"
LABELSET="labelset01"
CONFIG=""
NAME=""
FOREGROUND=0

if [[ $# -lt 1 ]]; then
  sed -n '3,30p' "$0" | sed 's/^# \{0,1\}//'
  exit 1
fi

SCRIPT="$1"; shift

while [[ $# -gt 0 ]]; do
  case "$1" in
    --freeze)   FREEZE="$2";   shift 2 ;;
    --cohort)   COHORT="$2";   shift 2 ;;
    --labelset) LABELSET="$2"; shift 2 ;;
    --config)   CONFIG="$2";   shift 2 ;;
    --name)     NAME="$2";     shift 2 ;;
    --fg)       FOREGROUND=1;  shift   ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# ---- checks --------------------------------------------------
if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: script not found: $SCRIPT" >&2
  exit 1
fi
if [[ -n "$CONFIG" && ! -f "$CONFIG" ]]; then
  echo "ERROR: config not found: $CONFIG" >&2
  exit 1
fi
if [[ ! -f "config/freezes/${FREEZE}.yml" ]]; then
  echo "ERROR: no such freeze: config/freezes/${FREEZE}.yml" >&2
  exit 1
fi
if [[ ! -f "config/cohorts/${COHORT}.yml" ]]; then
  echo "ERROR: no such cohort: config/cohorts/${COHORT}.yml" >&2
  exit 1
fi
if ! command -v screen >/dev/null 2>&1; then
  echo "ERROR: screen is not installed on this host." >&2
  exit 1
fi

mkdir -p logs

# Warn if the working tree is dirty. The run will record a git SHA that
# does not fully describe the code that produced it.
if command -v git >/dev/null 2>&1 && git rev-parse --git-dir >/dev/null 2>&1; then
  if [[ -n "$(git status --porcelain)" ]]; then
    echo "WARNING: uncommitted changes in the working tree." >&2
    echo "         This run's git_sha will not fully describe the code." >&2
  fi
fi

# ---- naming --------------------------------------------------
STAMP="$(date +%Y%m%d_%H%M%S)"
BASE="$(basename "$SCRIPT")"; BASE="${BASE%.*}"
if [[ -z "$NAME" ]]; then
  NAME="${BASE}_${COHORT}"
fi
LOG="logs/${STAMP}_${BASE}_${FREEZE}_${COHORT}.log"

# Refuse to start a second session with the same name — otherwise two
# jobs write into the same run directory and race.
if screen -ls 2>/dev/null | grep -q "[.]${NAME}[[:space:]]"; then
  echo "ERROR: a screen session named '${NAME}' is already running." >&2
  echo "       screen -r ${NAME}   to attach, or pass --name to differentiate." >&2
  exit 1
fi

# ---- environment passed to R --------------------------------
export ASSET_REPO="$REPO"
export ASSET_FREEZE="$FREEZE"
export ASSET_COHORT="$COHORT"
export ASSET_LABELSET="$LABELSET"
[[ -n "$CONFIG" ]] && export ASSET_CONFIG="$CONFIG"

# Keep BLAS/OpenMP from oversubscribing a shared machine. Raise
# deliberately for a job you know needs it.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-4}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-4}"

# ---- header written into the log ----------------------------
{
  echo "=============================================================="
  echo " script   : $SCRIPT"
  echo " freeze   : $FREEZE"
  echo " cohort   : $COHORT"
  echo " labelset : $LABELSET"
  [[ -n "$CONFIG" ]] && echo " config   : $CONFIG"
  echo " session  : $NAME"
  echo " started  : $(date '+%Y-%m-%d %H:%M:%S %Z')"
  echo " host     : $(hostname)"
  echo " user     : $USER"
  echo " git sha  : $(git rev-parse --short HEAD 2>/dev/null || echo NA)"
  echo " git dirty: $([[ -n "$(git status --porcelain 2>/dev/null)" ]] && echo yes || echo no)"
  echo " threads  : $OMP_NUM_THREADS"
  echo "=============================================================="
} > "$LOG"

# ---- launch --------------------------------------------------
RCMD="Rscript --no-save --no-restore '$SCRIPT'"

if [[ "$FOREGROUND" -eq 1 ]]; then
  echo "Running in foreground; output -> $LOG"
  eval "$RCMD" 2>&1 | tee -a "$LOG"
  exit "${PIPESTATUS[0]}"
fi

# -d -m : start detached immediately
# -S    : name the session
# The wrapper appends an exit-status line so a failed job is obvious in
# the log even if you were not watching.
screen -dmS "$NAME" bash -lc "
  cd '$REPO'
  { eval $RCMD ; } >> '$LOG' 2>&1
  code=\$?
  {
    echo '--------------------------------------------------------------'
    echo \" finished : \$(date '+%Y-%m-%d %H:%M:%S %Z')\"
    echo \" exit code: \$code\"
    if [ \$code -ne 0 ]; then echo ' STATUS   : FAILED'; else echo ' STATUS   : OK'; fi
    echo '--------------------------------------------------------------'
  } >> '$LOG'
"

sleep 1
cat <<EOF

Launched in detached screen session: ${NAME}

  watch output     tail -f ${LOG}
  attach           screen -r ${NAME}
  detach again     Ctrl-A then D
  list sessions    screen -ls
  overview         jobs/status.sh
  kill it          screen -S ${NAME} -X quit

Your prompt is free. The job survives logout and VPN drops.
EOF
