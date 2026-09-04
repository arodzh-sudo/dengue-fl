#!/usr/bin/env python3
"""
Turn a Daytona_dengue summary_report.txt into a pre-filled local metadata table.

Fills in sample_id, serotype, and the v-gen-lab lineage call the pipeline already
made, and leaves collection_date and the epidemiological columns blank for a
human. Only samples
whose vadr_flag is accepted are carried through, so this is also where the
decision about what enters the trees gets made.
"""

import argparse
import csv
import re
import sys

TEMPLATE_COLUMNS = [
    "sample_id",
    "serotype",
    "nextclade_clade",
    "collection_date",
    "location",
    "case_origin",
    "travel_country",
    "host",
    "strain",
    "notes",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-report", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--vadr-flags",
        default="PASS",
        help="comma-separated vadr_flag values to keep, or 'any' to keep every sample",
    )
    parser.add_argument(
        "--vector-pattern",
        default=r"(?i)mosquito",
        help="sample_id regex marking a vector pool rather than a human case",
    )
    parser.add_argument(
        "--vector-host",
        default="Aedes aegypti",
        help="host for samples matching --vector-pattern; must appear in ingest's host map",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    keep_any = args.vadr_flags.strip().lower() == "any"
    keep = {flag.strip().upper() for flag in args.vadr_flags.split(",")}

    with open(args.summary_report, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    for column in ("sample_id", "serotype", "vadr_flag", "nextclade_clade"):
        if not rows or column not in rows[0]:
            raise SystemExit(f"{args.summary_report} has no {column!r} column")

    vector = re.compile(args.vector_pattern)

    kept = []
    vectors = []
    dropped_qc = []
    dropped_serotype = []

    for row in rows:
        sample_id = (row.get("sample_id") or "").strip()
        serotype = (row.get("serotype") or "").strip()
        vadr = (row.get("vadr_flag") or "").strip().upper()

        if not keep_any and vadr not in keep:
            dropped_qc.append((sample_id, vadr or "NA"))
            continue
        # "unclassified" and "NA" carry no serotype digit, so there is no
        # v-gen-lab dataset to place the sample against.
        if not re.search(r"[1-4]", serotype):
            dropped_serotype.append((sample_id, serotype or "NA"))
            continue

        row_out = {
            "sample_id": sample_id,
            "serotype": serotype,
            "nextclade_clade": (row.get("nextclade_clade") or "").strip(),
        }
        # A mosquito pool is a vector sample, not a human case. Leaving host at
        # the workflow default would label it Homo sapiens, and host_genus and
        # host_type are derived from host, so the error would propagate.
        if vector.search(sample_id):
            row_out["host"] = args.vector_host
            vectors.append(sample_id)
        kept.append(row_out)

    with open(args.output, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(
            out, fieldnames=TEMPLATE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in kept:
            writer.writerow({column: row.get(column, "") for column in TEMPLATE_COLUMNS})

    print(f"Wrote {len(kept)} samples to {args.output}", file=sys.stderr)
    if vectors:
        print(f"\nTreated {len(vectors)} as {args.vector_host} vector pools:", file=sys.stderr)
        for sample_id in vectors:
            print(f"  {sample_id}", file=sys.stderr)
    if dropped_qc:
        print(f"\nDropped {len(dropped_qc)} on vadr_flag (keeping {sorted(keep)}):", file=sys.stderr)
        for sample_id, flag in dropped_qc:
            print(f"  {sample_id}\t{flag}", file=sys.stderr)
    if dropped_serotype:
        print(f"\nDropped {len(dropped_serotype)} with no assignable serotype:", file=sys.stderr)
        for sample_id, serotype in dropped_serotype:
            print(f"  {sample_id}\t{serotype}", file=sys.stderr)
    print(f"\nFill in collection_date for all {len(kept)} rows before running the workflow.", file=sys.stderr)


if __name__ == "__main__":
    main()
