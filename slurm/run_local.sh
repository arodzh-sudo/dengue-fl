#!/bin/bash
#SBATCH --job-name=dengue-local
#SBATCH --account=CHANGEME
#SBATCH --qos=CHANGEME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=01:00:00
#SBATCH --output=logs/slurm/local_%j.out
#SBATCH --error=logs/slurm/local_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=CHANGEME@ufl.edu

# Validates the hand-filled metadata, assigns lineages to the local consensus
# genomes, and writes the per-serotype files the phylogenetic build merges in.
#
# Set nextclade.dataset_zip in local/defaults/config.yaml to reuse the datasets
# ingest already downloaded, and this job needs no network at all.

set -euo pipefail

module purge
module load conda
# shellcheck disable=SC1091
source activate dengue-fl

mkdir -p logs/slurm
export TMPDIR="${SLURM_SUBMIT_DIR}/.tmp/${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

cd "${SLURM_SUBMIT_DIR}/local"

snakemake \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --rerun-incomplete \
    --printshellcmds \
    --show-failed-logs

# Put the report in the job output so a failure explains itself in the email.
echo "=============== LOCAL METADATA VALIDATION REPORT ==============="
cat results/validation_report.txt
