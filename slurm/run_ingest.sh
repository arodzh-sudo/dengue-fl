#!/bin/bash
#SBATCH --job-name=dengue-ingest
#SBATCH --account=CHANGEME
#SBATCH --qos=CHANGEME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=08:00:00
#SBATCH --output=logs/slurm/ingest_%j.out
#SBATCH --error=logs/slurm/ingest_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=CHANGEME@ufl.edu

# Downloads dengue sequences from NCBI, curates them, and assigns v-gen-lab
# lineages. Needs outbound network for `datasets download` and
# `nextclade dataset get`.
#
# If compute nodes on your partition cannot reach the internet, prefetch on a
# LOGIN node first, then submit this script:
#
#   cd ingest
#   snakemake --cores 2 --notemp \
#       data/ncbi_dataset.zip \
#       data/nextclade_data/v-gen-lab/denv{1,2,3,4}.zip
#
# --notemp matters: fetch_ncbi_dataset_package marks data/ncbi_dataset.zip as
# temp(), so without it Snakemake deletes the file you just went to the login
# node to download.

set -euo pipefail

module purge
module load conda
# shellcheck disable=SC1091
source activate dengue-fl

mkdir -p logs/slurm
export TMPDIR="${SLURM_SUBMIT_DIR}/.tmp/${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR}/ingest"

snakemake \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --rerun-incomplete \
    --printshellcmds \
    --show-failed-logs
