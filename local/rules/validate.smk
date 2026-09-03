"""
This part of the workflow validates the hand-filled metadata against the
consensus FASTA and prepares per-serotype sequence files for Nextclade.

REQUIRED INPUTS:

    metadata  = config.local_metadata
    sequences = config.local_sequences

OUTPUTS:

    metadata  = data/metadata_validated.tsv
    sequences = data/sequences_all.fasta
                data/sequences_{denv1..denv4}.fasta
    report    = results/validation_report.txt
"""

import json


rule validate_local_metadata:
    """
    Check the hand-filled TSV against the FASTA and emit a metadata table with
    the columns ingest produces. Fails the run on any structural problem.
    """
    input:
        metadata=config["local_metadata"],
        sequences=config["local_sequences"],
    output:
        metadata="data/metadata_validated.tsv",
        report="results/validation_report.txt",
    params:
        # These are lambdas because Snakemake expands wildcards in plain param
        # strings, and both the JSON blobs and the strain template contain
        # braces that are not wildcards.
        constants=lambda wildcards: json.dumps(config["constants"]),
        rules=lambda wildcards: json.dumps(config["validate"]),
        strain_template=lambda wildcards: config["strain_template"],
        strict="--strict" if config["validate"]["strict"] else "",
    log:
        "logs/validate_local_metadata.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/validate-local-metadata.py \
            --metadata {input.metadata:q} \
            --sequences {input.sequences:q} \
            --constants {params.constants:q} \
            --rules {params.rules:q} \
            --strain-template {params.strain_template:q} \
            --output-metadata {output.metadata:q} \
            --output-report {output.report:q} \
            {params.strict}
        """


rule normalize_local_fasta:
    """
    Rewrite FASTA headers to the bare sample identifier, matching the bare
    accession headers ingest writes.
    """
    input:
        sequences=config["local_sequences"],
        metadata="data/metadata_validated.tsv",
    output:
        sequences="data/sequences_all.fasta",
    log:
        "logs/normalize_local_fasta.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/normalize-local-fasta.py \
            --sequences {input.sequences:q} \
            --metadata {input.metadata:q} \
            --id-column accession \
            --output {output.sequences:q}
        """


rule split_sequences_for_nextclade:
    """
    Each serotype is genotyped against its own v-gen-lab dataset, so sequences
    are split before Nextclade runs.
    """
    input:
        sequences="data/sequences_all.fasta",
        metadata="data/metadata_validated.tsv",
    output:
        sequences=expand("data/sequences_{serotype}.fasta", serotype=NEXTCLADE_SEROTYPES),
    log:
        "logs/split_sequences_for_nextclade.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/split-by-serotype.py \
            --metadata {input.metadata:q} \
            --sequences {input.sequences:q} \
            --id-column accession \
            --serotype-column serotype_genbank \
            --serotypes denv1 denv2 denv3 denv4 \
            --output-dir data \
            --sequences-only
        """
