#!/bin/bash
#SBATCH --job-name=dengue-phylo
#SBATCH --account=CHANGEME
#SBATCH --qos=CHANGEME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128gb
#SBATCH --time=72:00:00
#SBATCH --output=logs/slurm/phylo_%j.out
#SBATCH --error=logs/slurm/phylo_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=CHANGEME@ufl.edu

# Builds all ten Auspice datasets: 5 serotypes x 2 genes.
#
# Wall clock is driven by augur tree (iqtree) and augur refine (treetime) over
# 4000-tip alignments. Time and memory here are a starting point; tune them from
# phylogenetic/benchmarks/ after the first run.
#
# If this exceeds your QoS time limit, submit one job per serotype instead, each
# targeting its two JSONs:
#   snakemake ... auspice/dengue_denv1_genome.json auspice/dengue_denv1_E.json

set -euo pipefail

module purge
module load conda
# shellcheck disable=SC1091
source activate dengue-fl

mkdir -p logs/slurm
export TMPDIR="${SLURM_SUBMIT_DIR}/.tmp/${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

# augur align and augur tree already get 8 threads each from Snakemake. Stop
# numpy inside concurrent treetime runs from oversubscribing the allocation.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR}/phylogenetic"

snakemake \
    --configfile build-configs/florida/config.yaml \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --rerun-incomplete \
    --printshellcmds \
    --keep-going \
    --show-failed-logs
