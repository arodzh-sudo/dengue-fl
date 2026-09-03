# dengue-fl

Dengue virus phylogenetics for the Florida Bureau of Public Health Laboratories,
forked from [nextstrain/dengue](https://github.com/nextstrain/dengue).

This fork runs entirely on local infrastructure. It fetches public GenBank data
for global context, adds consensus genomes from the
[Daytona_dengue](https://github.com/BPHL-Molecular/Daytona_dengue) pipeline, and
writes Auspice JSON files that colleagues open by drag and drop at
[auspice.us](https://auspice.us). Nothing reads from or writes to AWS.

## Workflows

Run them in this order.

1. [`ingest/`](./ingest) downloads dengue sequences from NCBI GenBank, curates
   them, and assigns dengue lineages with the `community/v-gen-lab/dengue`
   Nextclade datasets. Re-run it when you want fresher public context.
2. [`local/`](./local) validates your hand-filled metadata, assigns lineages to
   your consensus genomes with the same Nextclade datasets, and writes
   per-serotype files. Re-run it whenever you have new samples.
3. [`phylogenetic/`](./phylogenetic) merges the two and builds the trees.

[`nextclade/`](./nextclade) builds the Nextstrain genotype datasets and is
inherited from upstream. It is not part of the Florida build.

## Output

The Florida build produces ten Auspice datasets, the full cross product of five
serotype groupings and two genes:

| | Whole genome | E gene |
|---|---|---|
| All serotypes together | `dengue_all_genome.json` | `dengue_all_E.json` |
| DENV1 | `dengue_denv1_genome.json` | `dengue_denv1_E.json` |
| DENV2 | `dengue_denv2_genome.json` | `dengue_denv2_E.json` |
| DENV3 | `dengue_denv3_genome.json` | `dengue_denv3_E.json` |
| DENV4 | `dengue_denv4_genome.json` | `dengue_denv4_E.json` |

Each comes with a `_tip-frequencies.json` sidecar. All twenty land in
`phylogenetic/auspice/`.

Every tree carries the dengue lineage system of Hill et al. 2024 at three levels
of detail, from the Nextclade datasets at
[community/v-gen-lab/dengue](https://github.com/nextstrain/nextclade_data/tree/master/data/community/v-gen-lab/dengue),
plus a "Data source" coloring that separates local sequences from public ones.

## Quick start

```sh
conda env create -f envs/dengue-fl.yaml
conda activate dengue-fl

cat /path/to/daytona_output/assemblies_qc_pass/*.fasta > local/input/sequences.fasta
cp local/defaults/metadata_template.tsv local/input/metadata.tsv   # then fill it in

(cd ingest && snakemake --cores 8)
(cd local  && snakemake --cores 4)
(cd phylogenetic && snakemake --configfile build-configs/florida/config.yaml --cores 32)
```

On UF HiPerGator, use the SLURM scripts in [`slurm/`](./slurm) instead. See
[`slurm/README.md`](./slurm/README.md).

## Private data

`local/input/`, `local/data/`, and `local/results/` are all gitignored, so
consensus sequences and patient metadata cannot be committed by accident. Keep
it that way.

## Documentation

- [Contributor documentation](./CONTRIBUTING.md)
