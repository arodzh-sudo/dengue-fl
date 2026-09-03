"""
This part of the workflow merges the Nextclade calls into the validated
metadata, derives the lineage levels, and writes the per-serotype files the
phylogenetic workflow consumes as additional_inputs.

REQUIRED INPUTS:

    metadata  = data/metadata_validated.tsv
    nextclade = data/nextclade_metadata.tsv
    sequences = data/sequences_all.fasta

OUTPUTS:

    metadata  = results/metadata_{serotype}.tsv
    sequences = results/sequences_{serotype}.fasta
    include   = results/include_{serotype}.txt
"""


rule append_nextclade_columns:
    """Mirrors ingest rule append_nextclade_and_gene_coverage_columns, without the
    gene coverage columns, which the phylogenetic workflow does not read."""
    input:
        metadata="data/metadata_validated.tsv",
        nextclade="data/nextclade_metadata.tsv",
    output:
        metadata="data/metadata_nextclade.tsv",
    params:
        id_field="accession",
    log:
        "logs/append_nextclade_columns.txt",
    shell:
        r"""
        augur merge \
            --metadata \
                metadata={input.metadata:q} \
                nextclade={input.nextclade:q} \
            --metadata-id-columns \
                metadata={params.id_field:q} \
                nextclade={params.id_field:q} \
            --output-metadata {output.metadata:q} \
            --no-source-columns \
        &> {log:q}
        """


rule infer_major_lineage:
    """
    Derive the three lineage levels from the Nextclade call.

    Keep this identical to rule infer_major_lineage in
    ingest/rules/nextclade.smk. If the two diverge, local and public rows stop
    agreeing on what major_lineage means and the Auspice coloring becomes
    meaningless.

      1I_A.2.3 -> genotype 1I, major_lineage 1I_A, minor_lineage 1I_A.2.3
    """
    input:
        metadata="data/metadata_nextclade.tsv",
    output:
        metadata="data/metadata_lineages.tsv",
    params:
        nextclade_field="genotype_nextclade",
    log:
        "logs/infer_major_lineage.txt",
    shell:
        """
        cat {input.metadata:q} \
        | csvtk -tl mutate \
          -f {params.nextclade_field} \
          -n genotype \
          -p "^([0-9][A-Z]+)" \
        | csvtk -tl mutate \
          -f {params.nextclade_field} \
          -n major_lineage \
          -p "^([0-9][A-Z]+(?:_[A-Z])?)" \
        | csvtk -tl mutate \
          -f {params.nextclade_field} \
          -n minor_lineage \
        > {output.metadata:q}
        """


rule split_outputs_by_serotype:
    """
    Write all fifteen outputs in one rule so that "no samples for this serotype"
    is an ordinary path rather than an error. The phylogenetic workflow expands
    additional_inputs over all five serotypes on every run, so every file has to
    exist even when it is empty.
    """
    input:
        metadata="data/metadata_lineages.tsv",
        sequences="data/sequences_all.fasta",
    output:
        metadata=expand("results/metadata_{serotype}.tsv", serotype=SEROTYPES),
        sequences=expand("results/sequences_{serotype}.fasta", serotype=SEROTYPES),
        include=expand("results/include_{serotype}.txt", serotype=SEROTYPES),
    params:
        serotypes=" ".join(SEROTYPES),
        columns=",".join(config["metadata_columns"]),
    log:
        "logs/split_outputs_by_serotype.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/split-by-serotype.py \
            --metadata {input.metadata:q} \
            --sequences {input.sequences:q} \
            --id-column accession \
            --serotype-column serotype_genbank \
            --serotypes {params.serotypes} \
            --columns {params.columns:q} \
            --output-dir results \
            --write-include
        """
