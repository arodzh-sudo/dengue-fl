#!/usr/bin/env python3
"""
Split local metadata and sequences into per-serotype files.

`augur filter` is deliberately not used here: its --empty-output-reporting
defaults to error, so a run containing no samples for one serotype would abort.
Every requested serotype always gets a file, empty or not, because the
phylogenetic workflow's additional_inputs paths must resolve for all of them.
"""

import argparse
import csv
import sys
from datetime import date

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--id-column", default="accession")
    parser.add_argument("--serotype-column", default="serotype_genbank")
    parser.add_argument("--serotypes", nargs="+", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--columns", help="comma-separated output column order")
    parser.add_argument("--sequences-only", action="store_true")
    parser.add_argument("--write-include", action="store_true")
    parser.add_argument(
        "--empty-sentinel",
        action="store_true",
        help="give empty serotypes one placeholder row that augur filter's "
        "--exclude-where country=? region=? date=? removes",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.metadata, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        present = list(reader.fieldnames or [])

    columns = [c.strip() for c in args.columns.split(",")] if args.columns else present
    for column in present:
        if column not in columns:
            columns.append(column)

    sequences = {}
    for record in SeqIO.parse(args.sequences, "fasta"):
        sequences[record.id.split()[0].strip()] = str(record.seq)

    for serotype in args.serotypes:
        if serotype == "all":
            selected = rows
        else:
            selected = [row for row in rows if row.get(args.serotype_column) == serotype]
        ids = [row[args.id_column] for row in selected]

        fasta_path = f"{args.output_dir}/sequences_{serotype}.fasta"
        with open(fasta_path, "w", encoding="utf-8") as out:
            for sample_id in ids:
                if sample_id in sequences:
                    out.write(f">{sample_id}\n{sequences[sample_id]}\n")

        if args.sequences_only:
            print(f"{serotype}: {len(ids)} sequences", file=sys.stderr)
            continue

        tsv_path = f"{args.output_dir}/metadata_{serotype}.tsv"
        with open(tsv_path, "w", newline="", encoding="utf-8") as out:
            writer = csv.DictWriter(
                out, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore"
            )
            writer.writeheader()
            for row in selected:
                writer.writerow({column: row.get(column, "") for column in columns})
            if not selected and args.empty_sentinel:
                sentinel = {column: "" for column in columns}
                sentinel[args.id_column] = f"__no_local_{serotype}__"
                sentinel.update({"date": "?", "country": "?", "region": "?"})
                writer.writerow(sentinel)

        if args.write_include:
            include_path = f"{args.output_dir}/include_{serotype}.txt"
            with open(include_path, "w", encoding="utf-8") as out:
                out.write(
                    f"# Local sequences force-included past subsampling. "
                    f"Generated {date.today().isoformat()} by local/.\n"
                )
                for sample_id in ids:
                    out.write(f"{sample_id}\n")

        print(f"{serotype}: {len(selected)} rows, {len(ids)} sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
