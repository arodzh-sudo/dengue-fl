# local

Prepares privately sequenced dengue consensus genomes so the phylogenetic
workflow can place them on a tree alongside public GenBank data.

This workflow validates a hand-filled metadata table, assigns dengue lineages
with the same Nextclade datasets `ingest` uses, and writes one metadata file,
one sequence file, and one include list per serotype.

## What you provide

Two files, both in `local/input/`, which is gitignored so consensus sequences
cannot be committed by accident.

### `local/input/sequences.fasta`

Every consensus genome for this run, concatenated into one FASTA. From a
`Daytona_dengue` run that is:

```sh
cat /path/to/daytona_output/assemblies_qc_pass/*.fasta > local/input/sequences.fasta
```

Each record header must be the sample identifier, either alone or as the first
whitespace-delimited token. The workflow rewrites headers to the bare
identifier, matching the bare accessions `ingest` writes.

### `local/input/metadata.tsv`

Copy `defaults/metadata_template.tsv` and fill one row per sample.

| Column | Required | Notes |
|---|---|---|
| `sample_id` | yes | must match a FASTA header, and must be unique |
| `serotype` | yes | `denv1` through `denv4`; `DENV1` and `1` are accepted and normalized |
| `collection_date` | yes | `YYYY-MM-DD`, or `YYYY-MM-XX` / `YYYY-XX-XX` when partial |
| `location` | no | county |
| `case_origin` | no | `local`, `travel-associated`, or `unknown` |
| `travel_country` | no | country of exposure for travel-associated cases |
| `host` | no | defaults to `Homo sapiens` |
| `strain` | no | overrides the derived Auspice display name |
| `country`, `region`, `division` | no | default to `USA`, `North America`, `Florida` |
| `authors`, `institution`, `notes` | no | default to the values in `defaults/config.yaml` |

Everything else is filled in automatically. `serotype_genbank`, `is_lab_host`,
`host_genus`, `host_type`, `length`, `data_source`, and the four lineage columns
are all derived, and they use the same values and spellings as the public
metadata so that colorings do not split into duplicate categories.

## Sample identifiers

```
FL-BPHL-<YY>-<NNNN>
```

For example `FL-BPHL-26-0042`. The identifier lands in the `accession` column,
which is what the phylogenetic workflow keys tips on, and it must equal the
FASTA header.

Three properties matter, and the validator enforces all three:

- The `FL-` prefix makes a collision with a GenBank accession structurally
  impossible. On a collision, `augur merge` would silently keep one sequence and
  drop the other.
- Only letters, digits, and hyphens. Spaces, slashes, colons, commas, and pipes
  break Newick trees and Auspice node names.
- Stable across runs, so a sample keeps the same tip every time you rebuild.

Change `validate.id_regex` in `defaults/config.yaml` if your lab uses a different
scheme, but keep the three properties.

The Auspice display name is derived as `DENV1/USA/FL-BPHL-26-0042/2026`, matching
the pattern the public GenBank records use, so local and public tips read the
same way in the tree.

## Running it

```sh
cd local
snakemake --cores 4
```

Outputs land in `local/results/`:

- `metadata_{all,denv1..denv4}.tsv`
- `sequences_{all,denv1..denv4}.fasta`
- `include_{all,denv1..denv4}.txt`
- `validation_report.txt`

All five serotypes always get a file, even when a run has no samples for one of
them, because the phylogenetic workflow expands its input paths over every
serotype on every run.

Read `validation_report.txt` before releasing a build. It lists per-sample
length, unambiguous base count, and ambiguity fraction, along with every error
and warning.

When validation fails, Snakemake deletes the outputs of the failed rule, so
`results/validation_report.txt` will not be there. The same errors and warnings
are in `logs/validate_local_metadata.txt`, which Snakemake leaves alone. The
SLURM script passes `--show-failed-logs`, so on the cluster they also appear in
the job output.

## What the validator rejects

The run stops, and nothing is written, on any of these:

- a required column missing, or an unrecognized column present
- an empty, duplicated, or malformed `sample_id`
- a serotype outside `denv1` to `denv4`
- an unparseable, future, or pre-1950 collection date
- `country` or `region` resolving to empty or `?`, which `augur filter` silently
  drops later
- a sample in the metadata with no FASTA record, or the reverse
- a duplicate FASTA header
- a character outside the IUPAC nucleotide alphabet

Warnings do not stop the run unless you set `validate.strict: true`. The ones
worth reading are short-sequence and high-ambiguity warnings, because those
samples enter the genome tree only because `include.txt` forces them past
`augur filter --min-length 5000`.

## Reusing the Nextclade datasets from ingest

By default this workflow downloads the v-gen-lab datasets itself. To skip the
network and reuse what `ingest` already fetched, set in `defaults/config.yaml`:

```yaml
nextclade:
  dataset_zip: "../ingest/data/nextclade_data/v-gen-lab/{serotype}.zip"
```

If you instead let both workflows download, pin `nextclade.dataset_tag` to the
same value here and in `ingest/defaults/config.yaml`. Two different dataset
versions give local and public sequences lineage calls that are not comparable,
and nothing in the pipeline will warn you.

## A note on disclosure

The exported Auspice JSONs carry collection dates at day resolution and
county-level `location` values for local samples. Whether that is releasable
outside your organization is a policy question, not a technical one. Coarsen
`collection_date` to the month and drop `location` in the input file if the
answer is no.
