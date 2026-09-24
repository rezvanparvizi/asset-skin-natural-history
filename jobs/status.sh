#!/usr/bin/env bash
# ==============================================================
# jobs/status.sh — what is running, and how the recent runs ended
#
# Usage:
#   jobs/status.sh            sessions + last 5 logs
#   jobs/status.sh -n 20      last 20 logs
#   jobs/status.sh -f         follow the newest log (Ctrl-C to stop)
# ==============================================================

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

N=5
FOLLOW=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    -n) N="$2"; shift 2 ;;
    -f) FOLLOW=1; shift ;;
    *)  echo "Usage: jobs/status.sh [-n N] [-f]" >&2; exit 1 ;;
  esac
done

bold() { printf '\033[1m%s\033[0m\n' "$1"; }

# ---- running sessions ---------------------------------------
bold "Screen sessions"
if command -v screen >/dev/null 2>&1; then
  if screen -ls 2>/dev/null | grep -q "Socket\|Sockets\|Detached\|Attached"; then
    screen -ls || true
  else
    echo "  (none)"
  fi
else
  echo "  screen not available on this host"
fi
echo

# ---- resource use -------------------------------------------
bold "Your R processes"
# shellcheck disable=SC2009
ps -u "$USER" -o pid,etime,%cpu,%mem,rss,comm 2>/dev/null \
  | awk 'NR==1 || /R$|Rscript/' \
  | awk '{ if (NR==1) print "  "$0; else printf "  %s  %sMB\n", $0, int($5/1024) }' \
  || echo "  (none)"
echo

bold "Memory on this host"
free -h 2>/dev/null | sed 's/^/  /' || echo "  free(1) unavailable"
echo

# ---- recent logs --------------------------------------------
if [[ ! -d logs ]]; then
  echo "No logs/ directory yet."
  exit 0
fi

mapfile -t LOGS < <(ls -1t logs/*.log 2>/dev/null | head -n "$N")

if [[ ${#LOGS[@]} -eq 0 ]]; then
  echo "No logs yet."
  exit 0
fi

if [[ "$FOLLOW" -eq 1 ]]; then
  echo "Following ${LOGS[0]}  (Ctrl-C to stop)"
  echo
  tail -f "${LOGS[0]}"
  exit 0
fi

bold "Last $N runs"
for L in "${LOGS[@]}"; do
  STATUS="$(grep -m1 -o 'STATUS   : [A-Z]*' "$L" 2>/dev/null | awk '{print $3}')"
  [[ -z "$STATUS" ]] && STATUS="RUNNING?"
  SIZE="$(du -h "$L" | cut -f1)"
  printf '  %-10s %-6s %s\n' "$STATUS" "$SIZE" "$L"
done
echo

bold "Tail of newest log: ${LOGS[0]}"
tail -n 25 "${LOGS[0]}" | sed 's/^/  /'
echo

# ---- run registry -------------------------------------------
if [[ -f docs/runs.csv ]]; then
  bold "Last 8 registry entries (docs/runs.csv)"
  { head -1 docs/runs.csv; tail -n 8 docs/runs.csv; } \
    | cut -d, -f1,2,3,4,11,12 \
    | column -t -s, 2>/dev/null | sed 's/^/  /' \
    || tail -n 8 docs/runs.csv | sed 's/^/  /'
fi
