#!/bin/bash
# One-time conda environment setup. Run from the repo root on a LOGIN node.
set -euo pipefail

module purge
module load conda

# Keep conda's package cache and envs off $HOME, which has a small quota.
: "${DENGUE_FL_CONDA_ROOT:?set DENGUE_FL_CONDA_ROOT, e.g. /blue/<group>/$USER/conda}"
export CONDA_PKGS_DIRS="${DENGUE_FL_CONDA_ROOT}/pkgs"
export CONDA_ENVS_PATH="${DENGUE_FL_CONDA_ROOT}/envs"
mkdir -p "$CONDA_PKGS_DIRS" "$CONDA_ENVS_PATH"

conda env create -f envs/dengue-fl.yaml
# shellcheck disable=SC1091
source activate dengue-fl

missing=0
for binary in augur snakemake nextclade sqlite3 seqkit mafft iqtree \
              csvtk tsv-select tsv-filter tsv-join tsv-append \
              datasets dataformat unzip awk python3; do
  if ! command -v "$binary" >/dev/null; then
    echo "MISSING: $binary"
    missing=1
  fi
done
[[ $missing -eq 0 ]] || exit 1

augur --version
nextclade --version
snakemake --version

# augur merge shells out to sqlite3 and seqkit, so exercise it for real rather
# than trusting that the binaries are on PATH.
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
printf 'accession\tdate\nAB1\t2020-01-01\n' > "$tmp/full.tsv"
printf 'accession\tdate\n' > "$tmp/empty.tsv"
augur merge --metadata a="$tmp/full.tsv" b="$tmp/empty.tsv" \
    --metadata-id-columns accession --output-metadata "$tmp/out.tsv"
echo "augur merge tolerates a header-only table:"
cat "$tmp/out.tsv"

echo "Environment OK."
