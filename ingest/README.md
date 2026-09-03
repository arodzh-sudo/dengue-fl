# dengue-fl/ingest

Fetches dengue sequences and metadata from NCBI GenBank, curates them, and
assigns dengue lineages with the `community/v-gen-lab/dengue` Nextclade
datasets. Output is the public backbone the phylogenetic build places local
sequences against.

## Software requirements

Create the conda environment described in [`envs/dengue-fl.yaml`](../envs/dengue-fl.yaml):

```sh
bash slurm/setup_env.sh          # on HiPerGator, from the repo root
# or
conda env create -f envs/dengue-fl.yaml
```

## Usage

Run from the `ingest` directory:

```sh
cd ingest
snakemake --cores 8
```

Fetch the raw NDJSON without running the rest:

```sh
snakemake --cores 2 data/ncbi.ndjson
```

On HiPerGator, submit [`slurm/run_ingest.sh`](../slurm/run_ingest.sh) instead.
That script also documents what to do when compute nodes cannot reach NCBI.

This produces 10 files (within the `ingest` directory):

A pair of files with all the dengue sequences:

- `ingest/results/metadata_all.tsv`
- `ingest/results/sequences_all.fasta`

A pair of files for each dengue serotype (denv1 - denv4)

- `ingest/results/metadata_denv1.tsv`
- `ingest/results/sequences_denv1.fasta`
- `ingest/results/metadata_denv2.tsv`
- `ingest/results/sequences_denv2.fasta`
- `ingest/results/metadata_denv3.tsv`
- `ingest/results/sequences_denv3.fasta`
- `ingest/results/metadata_denv4.tsv`
- `ingest/results/sequences_denv4.fasta`


## Configuration

Configuration takes place in `defaults/config.yaml` by default.

## Input data

### GenBank data

GenBank sequences and metadata are fetched via [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).
