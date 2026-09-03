# Running on UF HiPerGator

Four scripts. Edit the `#SBATCH --account`, `--qos`, and `--mail-user` lines in
all three job scripts before the first submission.

Work from a `/blue` directory rather than `$HOME`. The workflows write tens of
gigabytes of intermediates and `$HOME` quotas are small.

## One-time setup

```sh
cd /blue/<group>/$USER/dengue-fl
export DENGUE_FL_CONDA_ROOT=/blue/<group>/$USER/conda
bash slurm/setup_env.sh
```

Run this on a login node, not through sbatch. Building a conda environment needs
network access, and the script finishes in minutes.

`setup_env.sh` does more than create the environment. It checks that every
binary the rules shell out to is on `PATH`, and it runs a real `augur merge` on a
header-only table. That second check matters twice over: `augur merge` depends on
`sqlite3` and `seqkit`, which are not `augur` dependencies, and the `local`
workflow relies on merging tables that are empty whenever a run contains no
samples of some serotype.

## Each run

```sh
cd /blue/<group>/$USER/dengue-fl
mkdir -p logs/slurm

# Stage the private inputs
cat /path/to/daytona_output/assemblies_qc_pass/*.fasta > local/input/sequences.fasta
cp local/defaults/metadata_template.tsv local/input/metadata.tsv   # then fill it in

# ingest and local are independent, so run them concurrently
JOB_INGEST=$(sbatch --parsable slurm/run_ingest.sh)
JOB_LOCAL=$(sbatch  --parsable slurm/run_local.sh)

# phylogenetic needs both
sbatch --dependency=afterok:${JOB_INGEST}:${JOB_LOCAL} slurm/run_phylo.sh
```

`local` reads nothing from `ingest`, so the two really are independent. The
dependency exists only because the phylogenetic build merges their outputs.

Once `run_phylo.sh` finishes, copy `phylogenetic/auspice/*.json` off the cluster
and hand them to colleagues. They drop them on
[auspice.us](https://auspice.us), which needs no server and nothing installed.

## Network access

`run_ingest.sh` is the only job that needs the internet, for the NCBI download
and the Nextclade dataset fetch. If your partition's compute nodes cannot reach
out, prefetch on a login node first:

```sh
cd ingest
snakemake --cores 2 --notemp \
    data/ncbi_dataset.zip \
    data/nextclade_data/v-gen-lab/denv{1,2,3,4}.zip
```

`--notemp` is not optional here. The fetch rule marks `data/ncbi_dataset.zip` as
a temporary file, so without it Snakemake deletes the download as soon as the
rules that consume it finish, and the sbatch job would try to fetch it again.

To skip the network in `run_local.sh` as well, point it at the datasets `ingest`
already downloaded, in `local/defaults/config.yaml`:

```yaml
nextclade:
  dataset_zip: "../ingest/data/nextclade_data/v-gen-lab/{serotype}.zip"
```

## Resources

The numbers in the scripts are a starting point, not measurements.

| Job | Cores | Memory | Time |
|---|---|---|---|
| `run_ingest.sh` | 8 | 32 GB | 8 h |
| `run_local.sh` | 4 | 8 GB | 1 h |
| `run_phylo.sh` | 32 | 128 GB | 72 h |

Tune them after the first successful run using the per-rule timings in
`phylogenetic/benchmarks/` and `ingest/benchmarks/`.

The phylogenetic job is the one that will need attention. Ten builds run
independently, and `augur align` and `augur tree` each declare eight threads, so
32 cores keeps four of them busy. Wall clock is dominated by `augur tree`
(iqtree) and `augur refine` (treetime, single-threaded) over 4000-tip
alignments.

If 72 hours is more than your QoS allows, submit one job per serotype and name
its two targets explicitly:

```sh
snakemake --configfile build-configs/florida/config.yaml --cores 8 \
    auspice/dengue_denv1_genome.json auspice/dengue_denv1_E.json
```

The `all` build is the long pole, since it carries every serotype.

## After a failure

Every rule writes to `logs/` and `benchmarks/` inside its workflow directory.
The job scripts pass `--show-failed-logs`, so the failing rule's log is already
in the SLURM output.

If a job is killed rather than failing cleanly, Snakemake leaves the outputs it
was mid-way through marked incomplete. The scripts pass `--rerun-incomplete`, so
resubmitting picks up from there.
