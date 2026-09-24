#!/usr/bin/env bash
# ==============================================================
# setup.sh — one-time setup on the server, after cloning
#
#   ./setup.sh /path/to/large/storage/asset
#
# Creates the data/results/figures/logs symlinks, the subdirectories
# under data/, sets permissions on the clinical directory, makes the
# job scripts executable, and runs the path check.
# ==============================================================

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: ./setup.sh /path/to/large/storage/asset" >&2
  echo >&2
  echo "Pick a location on LARGE storage. Confirm with your systems" >&2
  echo "contact which paths are backed up and which are purged before" >&2
  echo "putting anything irreplaceable there." >&2
  exit 1
fi

BIG="$1"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO"

echo "repo: $REPO"
echo "data: $BIG"
echo

mkdir -p "$BIG"/{data,results,figures,logs}
mkdir -p "$BIG"/data/{raw,objects,bpcells,bulk,clinical,external,inherited}
chmod 700 "$BIG/data/clinical"
echo "created directory tree under $BIG"

for d in data results figures logs; do
  if [[ -L "$d" ]]; then
    echo "  symlink $d already exists -> $(readlink "$d")"
  elif [[ -e "$d" ]]; then
    echo "  WARNING: $d exists and is not a symlink; leaving it alone" >&2
  else
    ln -s "$BIG/$d" "$d"
    echo "  linked $d -> $BIG/$d"
  fi
done

chmod +x jobs/run.sh jobs/status.sh R/paths_check.R setup.sh
echo "made scripts executable"

if [[ ! -f metadata/library_manifest.csv ]]; then
  cp metadata/library_manifest_TEMPLATE.csv metadata/library_manifest.csv
  echo "created metadata/library_manifest.csv from the template (gitignored)"
  echo "  -> replace the example rows with your real manifest"
fi

echo
echo "=== path check ==="
Rscript R/paths_check.R || true

cat <<EOF

Next:
  1. Fill in metadata/library_manifest.csv from the sample manifest.
  2. Fill in the TODOs in config/freezes/freeze01.yml (platform, probe
     set, CellRanger version, reference) from the sequencing core.
  3. Put the clinical table at data/clinical/asset_clinical.csv
     (never committed).
  4. Rscript -e 'renv::restore()'
  5. jobs/run.sh analysis/00_manifest_and_batch/00_1_build_manifest.R
  6. jobs/run.sh analysis/00_manifest_and_batch/00_2_batch_vs_design_crosstab.R

Step 6 is hazard #1. Read its VERDICT.txt before any biology.
EOF
