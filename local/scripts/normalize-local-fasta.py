#!/usr/bin/env python3
"""
Rewrite consensus FASTA headers to the bare sample identifier, matching the bare
accession headers that ingest produces, and drop records absent from the
validated metadata.
"""

import argparse
import csv
import sys

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--id-column", default="accession")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.metadata, newline="", encoding="utf-8") as handle:
        wanted = {row[args.id_column] for row in csv.DictReader(handle, delimiter="\t")}

    written = set()
    with open(args.output, "w", encoding="utf-8") as out:
        for record in SeqIO.parse(args.sequences, "fasta"):
            seq_id = record.id.split()[0].strip()
            if seq_id not in wanted or seq_id in written:
                continue
            written.add(seq_id)
            out.write(f">{seq_id}\n{str(record.seq).upper()}\n")

    missing = wanted - written
    if missing:
        raise SystemExit(f"no sequence found for: {', '.join(sorted(missing))}")

    print(f"Wrote {len(written)} sequences to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
