"""
This part of the workflow collects the phylogenetic tree and annotations to
export a Nextstrain dataset.
REQUIRED INPUTS:
    metadata        = results/{serotype}/metadata.tsv
    tree            = results/{serotype}/{gene}/tree.nwk
    branch_lengths  = results/{serotype}/{gene}/branch-lengths.json
    node_data       = results/{serotype}/{gene}/*.json
OUTPUTS:
    auspice_json      = auspice/dengue_{serotype}_{gene}.json
    frequencies_json  = auspice/dengue_{serotype}_{gene}_tip-frequencies.json
    There are optional sidecar JSON files that can be exported as part of the dataset.
    See Nextstrain's data format docs for more details on sidecar files:
    https://docs.nextstrain.org/page/reference/data-formats.html
This part of the workflow usually includes the following steps:
    - augur export v2
    - augur frequencies
See Augur's usage docs for these commands for more details.
"""

import json


rule colors:
    input:
        color_schemes = "defaults/color_schemes.tsv",
        color_orderings = "defaults/color_orderings.tsv",
        metadata = "results/{serotype}/metadata.tsv",
        manual_colors = "defaults/colors.tsv"
    output:
        colors = "results/{serotype}/colors.tsv"
    benchmark:
        "benchmarks/{serotype}/colors.txt"
    params:
        # Columns with too many values to tell apart by colour. Every observed
        # value is muted, then defaults/colors.tsv pins the handful that matter.
        muted_columns = ["country", "country_exposure", "location"],
        muted_color = "#dcdcdc",
    shell:
        """
        python3 scripts/assign-colors.py \
            --color-schemes {input.color_schemes} \
            --ordering {input.color_orderings} \
            --metadata {input.metadata} \
            --output {output.colors}.ramp

        # A 200-colour ramp is a smooth interpolation, so neighbouring countries
        # are indistinguishable. Mute them all, then let the manual file give
        # strong colours to the few that matter for this build.
        #
        # augur passes every row of the colours file through to the exported
        # scale without deduplicating, and Auspice then uses the first entry it
        # finds for a value. So the sections are concatenated in priority order,
        # manual pins first and the generated ramp last, and the whole thing is
        # deduplicated on trait and value keeping the first row. That leaves one
        # colour per value and relies on no undocumented precedence.
        #
        # The mute list comes from the metadata rather than a fixed list so that
        # no value, however it is spelled, falls through to Auspice's own scale.
        (
        cat {input.manual_colors}

        for column in {params.muted_columns}; do
            awk -F'\t' -v col="$column" -v grey='{params.muted_color}' -v OFS='\t' '
                NR == 1 {{ for (i = 1; i <= NF; i++) if ($i == col) c = i; next }}
                c {{
                    value = $c
                    # augur merge quotes any field that is not a bare single word,
                    # so the raw column holds "Sri Lanka" rather than Sri Lanka.
                    # augur reads its own quoting back off, so the colours have to
                    # be keyed on the unquoted value or they never match.
                    gsub(/^"|"$/, "", value)
                    if (value != "" && value != "?") seen[value] = 1
                }}
                END {{ for (value in seen) print col, value, grey }}
            ' {input.metadata}
        done

        cat {output.colors}.ramp
        ) | awk -F'\t' '$0 == "" || !seen[$1, $2]++' > {output.colors}

        rm {output.colors}.ramp
        """


rule prepare_auspice_config:
    """Prepare the auspice config file for each serotypes"""
    output:
        auspice_config="results/defaults/{serotype}/{gene}/auspice_config.json",
    benchmark:
        "benchmarks/{serotype}/{gene}/prepare_auspice_config.txt"
    params:
        replace_clade_key=lambda wildcard: r"clade_membership" if wildcard.gene in ['genome'] else r"major_lineage",
        replace_clade_title=lambda wildcard: r"Serotype" if wildcard.serotype in ['all'] else r"Genotype (Nextclade)",
    run:
        export_config = config.get("export", {})
        data = {
            "title": export_config.get("title", "Real-time tracking of dengue virus evolution"),
            "maintainers": export_config.get("maintainers", [
              {"name": "the Nextstrain team", "url": "https://nextstrain.org/team"}
            ]),
            "data_provenance": [
              {
                "name": "GenBank",
                "url": "https://www.ncbi.nlm.nih.gov/genbank/"
              }
            ],
            "build_url": export_config.get("build_url", "https://github.com/nextstrain/dengue"),
            "colorings": [
              {
                "key": "gt",
                "title": "Genotype",
                "type": "categorical"
              },
              {
                "key": "data_source",
                "title": "Data source",
                "type": "categorical"
              },
              {
                "key": "num_date",
                "title": "Date",
                "type": "continuous"
              },
              {
                "key": "country",
                "title": "Country",
                "type": "categorical"
              },
              {
                "key": "region",
                "title": "Region",
                "type": "categorical"
              },
              {
                "key": "serotype_genbank",
                "title": "Serotype (Genbank metadata)",
                "type": "categorical"
              },
              {
                "key": "genotype",
                "title": "Genotype (Nextclade)",
                "type": "categorical"
              },
              {
                "key": "major_lineage",
                "title": "Major lineage (Nextclade)",
                "type": "categorical"
              },
              {
                "key": "minor_lineage",
                "title": "Minor lineage (Nextclade)",
                "type": "categorical"
              },
              {
                "key": "host",
                "title": "Host Species",
                "type": "categorical"
              },
              {
                "key": "host_genus",
                "title": "Host Genus",
                "type": "categorical"
              },
              {
                "key": "host_type",
                "title": "Host Type",
                "type": "categorical"
              },
              {
                "key": "division",
                "title": "State or division",
                "type": "categorical"
              },
              {
                "key": "location",
                "title": "County",
                "type": "categorical"
              },
              {
                "key": "case_origin",
                "title": "Case origin",
                "type": "categorical"
              },
              {
                "key": "country_exposure",
                "title": "Country of exposure",
                "type": "categorical"
              },
              {
                "key": "region_exposure",
                "title": "Region of exposure",
                "type": "categorical"
              }
            ],
            "geo_resolutions": [
              "country",
              "region"
            ],
            "display_defaults": {
              "map_triplicate": True,
              "color_by": params.replace_clade_key,
              "tip_label": "strain"
            },
            "filters": [
              "data_source",
              "case_origin",
              "location",
              "division",
              "country",
              "region",
              "author"
            ],
            "panels": [
              "tree",
              "map",
              "entropy",
              "frequencies"
            ],
            "metadata_columns": [
              "accession",
              "strain",
              "url",
              "data_source",
              "division",
              "location"
            ]
          }

        # During genome/dengue_all workflows, clade membership represents Serotype
        # While genome/dengue_denvX workflows, clade_membership represents the more detailed Genotype
        if params.replace_clade_key == 'clade_membership':
            if wildcards.gene in ['genome'] and wildcards.serotype in ['all']:
                clade_membership_title="Serotype (Nextstrain)"
            else:
                clade_membership_title="Genotype (Nextstrain)"

            data["colorings"].append({
                "key": "clade_membership",
                "title": clade_membership_title,
                "type": "categorical"
            })
        else:
            # During E/dengue_all workflows, default color by Serotype
            if wildcards.serotype in ['all']:
                data["display_defaults"]["color_by"]="serotype_genbank"

        with open(output.auspice_config, 'w') as fh:
            json.dump(data, fh, indent=2)


rule export:
    """Exporting data files for auspice"""
    input:
        tree = "results/{serotype}/{gene}/tree.nwk",
        metadata = "results/{serotype}/metadata.tsv",
        branch_lengths = "results/{serotype}/{gene}/branch-lengths.json",
        traits = "results/{serotype}/{gene}/traits.json",
        clades = lambda wildcard: "results/{serotype}/{gene}/clades.json" if wildcard.gene in ['genome'] else [],
        nt_muts = "results/{serotype}/{gene}/nt-muts.json",
        aa_muts = "results/{serotype}/{gene}/aa-muts.json",
        description = config["export"]["description"],
        auspice_config = "results/defaults/{serotype}/{gene}/auspice_config.json",
        colors = "results/{serotype}/colors.tsv",
    output:
        auspice_json = "auspice/dengue_{serotype}_{gene}.json"
    benchmark:
        "benchmarks/{serotype}/{gene}/export.txt"
    params:
        strain_id = config.get("strain_id_field", "strain"),
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --node-data {input.branch_lengths} {input.traits} {input.clades} {input.nt_muts} {input.aa_muts} \
            --colors {input.colors} \
            --description {input.description} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence-inline \
            --output {output.auspice_json}
        """

rule tip_frequencies:
    """
    Estimating KDE frequencies for tips
    """
    input:
        tree = "results/{serotype}/{gene}/tree.nwk",
        metadata = "results/{serotype}/metadata.tsv",
    output:
        tip_freq = "auspice/dengue_{serotype}_{gene}_tip-frequencies.json"
    benchmark:
        "benchmarks/{serotype}/{gene}/tip_frequencies.txt"
    params:
        strain_id = config["strain_id_field"],
        min_date = config["tip_frequencies"]["min_date"],
        max_date = config["tip_frequencies"]["max_date"],
        narrow_bandwidth = config["tip_frequencies"]["narrow_bandwidth"],
        wide_bandwidth = config["tip_frequencies"]["wide_bandwidth"]
    shell:
        r"""
        augur frequencies \
            --method kde \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --wide-bandwidth {params.wide_bandwidth} \
            --output {output.tip_freq}
        """
