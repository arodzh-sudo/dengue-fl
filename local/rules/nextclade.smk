"""
This part of the workflow assigns dengue lineages to the local sequences using
the same community/v-gen-lab datasets that ingest uses for public sequences.

REQUIRED INPUTS:

    sequences = data/sequences_{denv1..denv4}.fasta

OUTPUTS:

    nextclade = data/nextclade_metadata.tsv
"""


def _nextclade_dataset(wildcards):
    override = config["nextclade"].get("dataset_zip") or ""
    if override:
        return override.format(serotype=wildcards.serotype)
    return f"data/nextclade_data/{wildcards.serotype}.zip"


rule get_nextclade_dataset:
    """Download the same v-gen-lab dataset that ingest uses."""
    output:
        dataset="data/nextclade_data/{serotype}.zip",
    params:
        dataset_name=lambda wildcards: config["nextclade"]["dataset_name"].format(
            serotype=wildcards.serotype
        ),
    wildcard_constraints:
        serotype="|".join(NEXTCLADE_SEROTYPES),
    retries: 3
    log:
        "logs/get_nextclade_dataset_{serotype}.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        nextclade dataset get \
            --name={params.dataset_name:q} \
            --output-zip={output.dataset:q}
        """


rule run_nextclade_local:
    """
    Genotype local sequences.

    A run can contain no samples of a given serotype, leaving an empty input
    FASTA. Rather than depend on how Nextclade handles that, the empty case
    writes the header-only TSV that the concat rule expects.
    """
    input:
        sequences="data/sequences_{serotype}.fasta",
        dataset=_nextclade_dataset,
    output:
        nextclade="data/nextclade/{serotype}/nextclade.tsv",
    threads: 2
    params:
        columns="\\t".join(config["nextclade"]["field_map"].keys()),
    wildcard_constraints:
        serotype="|".join(NEXTCLADE_SEROTYPES),
    log:
        "logs/run_nextclade_{serotype}.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        if [[ -s {input.sequences:q} ]]; then
            nextclade run \
                --input-dataset {input.dataset:q} \
                -j {threads} \
                --output-tsv {output.nextclade:q} \
                --silent \
                {input.sequences:q}
        else
            echo "No local {wildcards.serotype} sequences this run." >&2
            printf '{params.columns}\n' > {output.nextclade:q}
        fi
        """


rule concat_nextclade_results:
    """Mirrors ingest/rules/nextclade.smk rule concat_genotype_nextclade_results."""
    input:
        nextclade_files=expand(
            "data/nextclade/{serotype}/nextclade.tsv", serotype=NEXTCLADE_SEROTYPES
        ),
    output:
        genotype_nextclade="data/nextclade_metadata.tsv",
    params:
        input_nextclade_fields=",".join(config["nextclade"]["field_map"].keys()),
        output_nextclade_fields=",".join(config["nextclade"]["field_map"].values()),
    log:
        "logs/concat_nextclade_results.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        echo "{params.output_nextclade_fields}" \
        | tr ',' '\t' \
        > {output.genotype_nextclade:q}

        tsv-select -H -f "{params.input_nextclade_fields}" {input.nextclade_files} \
        | awk 'NR>1 {{print}}' \
        >> {output.genotype_nextclade:q}
        """
