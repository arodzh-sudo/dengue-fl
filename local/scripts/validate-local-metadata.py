#!/usr/bin/env python3
"""
Validate the hand-filled local metadata TSV against the consensus FASTA and emit
a metadata table with the same columns ingest produces.
"""

import argparse
import csv
import json
import re
import sys
from datetime import date, datetime

from Bio import SeqIO

UNAMBIGUOUS = set("ACGTU")
VALID_NT = set("ACGTURYSWKMBDHVN-.")
CASE_ORIGINS = {"local", "travel-associated", "unknown", ""}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--constants", required=True, help="JSON object of constant column values")
    parser.add_argument("--rules", required=True, help="JSON object of validation rules")
    parser.add_argument("--strain-template", required=True)
    parser.add_argument("--output-metadata", required=True)
    parser.add_argument("--output-report", required=True)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def read_metadata(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        lines = [line for line in handle if not line.lstrip().startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit(f"{path} is empty")
    columns = [name.strip() for name in reader.fieldnames]
    rows = []
    for row in reader:
        rows.append({key.strip(): (value or "").strip() for key, value in row.items() if key is not None})
    return columns, rows


def read_sequences(path):
    records = {}
    duplicates = []
    for record in SeqIO.parse(path, "fasta"):
        seq_id = record.id.split()[0].strip()
        if seq_id in records:
            duplicates.append(seq_id)
        records[seq_id] = str(record.seq).upper()
    return records, duplicates


def normalise_serotype(value):
    match = re.search(r"([1-4])", value or "")
    if not match:
        return ""
    return f"denv{match.group(1)}"


def parse_date(value):
    """Accept a full ISO date or a partial date with XX for unknown parts."""
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", value):
        return datetime.strptime(value, "%Y-%m-%d").date(), int(value[:4])
    if re.fullmatch(r"\d{4}-\d{2}-XX", value):
        return None, int(value[:4])
    if re.fullmatch(r"\d{4}-XX-XX", value):
        return None, int(value[:4])
    raise ValueError(value)


def sequence_stats(sequence):
    ungapped = sequence.replace("-", "").replace(".", "")
    unambiguous = sum(1 for base in ungapped if base in UNAMBIGUOUS)
    ambiguous_fraction = 0.0 if not ungapped else 1 - (unambiguous / len(ungapped))
    return len(ungapped), unambiguous, ambiguous_fraction


def main():
    args = parse_args()
    constants = json.loads(args.constants)
    rules = json.loads(args.rules)
    strict = args.strict or rules.get("strict", False)

    errors = []
    warnings = []

    columns, rows = read_metadata(args.metadata)
    sequences, duplicate_headers = read_sequences(args.sequences)

    known = set(rules["required_columns"]) | set(rules["optional_columns"])
    for column in columns:
        if column not in known:
            errors.append(f"unknown column {column!r}; expected one of {sorted(known)}")
    for column in rules["required_columns"]:
        if column not in columns:
            errors.append(f"required column {column!r} is missing")

    for seq_id in sorted(set(duplicate_headers)):
        errors.append(f"duplicate FASTA record {seq_id!r}")

    if errors:
        report_and_exit(args, errors, warnings, [])

    id_pattern = re.compile(rules["id_regex"])
    allowed = set(rules["allowed_serotypes"])
    today = date.today()

    seen = set()
    records = []
    for index, row in enumerate(rows, start=2):
        where = f"row {index}"
        sample_id = row.get("sample_id", "")

        if not sample_id:
            errors.append(f"{where}: empty sample_id")
            continue
        if sample_id in seen:
            errors.append(f"{where}: duplicate sample_id {sample_id!r}")
            continue
        seen.add(sample_id)
        if not id_pattern.fullmatch(sample_id):
            errors.append(f"{where}: sample_id {sample_id!r} does not match {rules['id_regex']}")
            continue

        serotype = normalise_serotype(row.get("serotype", ""))
        if serotype not in allowed:
            errors.append(f"{sample_id}: serotype {row.get('serotype', '')!r} is not one of {sorted(allowed)}")
            continue

        raw_date = row.get("collection_date", "")
        try:
            parsed, year = parse_date(raw_date)
        except ValueError:
            errors.append(
                f"{sample_id}: collection_date {raw_date!r} is not YYYY-MM-DD, YYYY-MM-XX, or YYYY-XX-XX"
            )
            continue
        if parsed and parsed > today:
            errors.append(f"{sample_id}: collection_date {raw_date} is in the future")
            continue
        if year < 1950:
            errors.append(f"{sample_id}: collection_date {raw_date} predates 1950")
            continue

        record = dict(constants)
        record.update({key: value for key, value in row.items() if value})

        # augur filter drops rows where any of these is a literal "?".
        for field in ("country", "region"):
            if record.get(field, "") in ("", "?"):
                errors.append(f"{sample_id}: {field} resolved to {record.get(field, '')!r}")

        host = record.get("host", "")
        record.update(
            {
                "accession": sample_id,
                "accession_version": sample_id,
                "date": raw_date,
                "serotype_genbank": serotype,
                "host": host,
                "date_released": today.isoformat(),
                "date_updated": today.isoformat(),
            }
        )
        for consumed in ("sample_id", "collection_date", "serotype"):
            record.pop(consumed, None)

        if not record.get("strain"):
            record["strain"] = args.strain_template.format(
                serotype_number=serotype[-1],
                country=record["country"].replace(" ", "_").upper(),
                accession=sample_id,
                year=year,
            )

        case_origin = record.get("case_origin", "")
        if case_origin not in CASE_ORIGINS:
            warnings.append(f"{sample_id}: case_origin {case_origin!r} is not one of {sorted(CASE_ORIGINS)}")
        if record.get("travel_country") and case_origin != "travel-associated":
            warnings.append(f"{sample_id}: travel_country is set but case_origin is {case_origin!r}")

        records.append(record)

    metadata_ids = {record["accession"] for record in records}
    missing_sequence = sorted(metadata_ids - set(sequences))
    missing_metadata = sorted(set(sequences) - metadata_ids)
    for sample_id in missing_sequence:
        errors.append(f"{sample_id}: present in metadata but has no FASTA record")
    for seq_id in missing_metadata:
        errors.append(f"{seq_id}: present in the FASTA but has no metadata row")

    stats = []
    for record in records:
        sequence = sequences.get(record["accession"])
        if sequence is None:
            continue
        invalid = set(sequence) - VALID_NT
        if invalid:
            errors.append(f"{record['accession']}: sequence contains {sorted(invalid)}")
        length, unambiguous, ambiguous_fraction = sequence_stats(sequence)
        record["length"] = str(length)
        stats.append((record["accession"], record["serotype_genbank"], length, unambiguous, ambiguous_fraction))

        if unambiguous < rules["min_ungapped_length"]:
            warnings.append(
                f"{record['accession']}: {unambiguous} unambiguous bases is below "
                f"{rules['min_ungapped_length']}; it will only enter the genome tree because "
                "include.txt forces it past --min-length"
            )
        if unambiguous < 1000:
            warnings.append(
                f"{record['accession']}: {unambiguous} unambiguous bases is below 1000; "
                "E gene extraction will drop it"
            )
        if ambiguous_fraction > rules["max_ambiguous_fraction"]:
            warnings.append(
                f"{record['accession']}: {ambiguous_fraction:.1%} ambiguous exceeds "
                f"{rules['max_ambiguous_fraction']:.0%}"
            )

    if errors or (strict and warnings):
        report_and_exit(args, errors, warnings, stats)

    write_report(args.output_report, errors, warnings, stats)

    fieldnames = sorted({key for record in records for key in record})
    with open(args.output_metadata, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow(record)

    print(f"Validated {len(records)} local samples.", file=sys.stderr)


def write_report(path, errors, warnings, stats):
    lines = ["# Local metadata validation report", ""]

    by_serotype = {}
    for _, serotype, _, _, _ in stats:
        by_serotype[serotype] = by_serotype.get(serotype, 0) + 1
    lines.append(f"Samples: {len(stats)}")
    for serotype in sorted(by_serotype):
        lines.append(f"  {serotype}: {by_serotype[serotype]}")
    lines.append("")

    lines.append("sample\tserotype\tlength\tunambiguous\tambiguous_fraction")
    for accession, serotype, length, unambiguous, fraction in sorted(stats):
        lines.append(f"{accession}\t{serotype}\t{length}\t{unambiguous}\t{fraction:.4f}")
    lines.append("")

    lines.append(f"Errors: {len(errors)}")
    lines.extend(f"  ERROR {message}" for message in errors)
    lines.append("")
    lines.append(f"Warnings: {len(warnings)}")
    lines.extend(f"  WARNING {message}" for message in warnings)
    lines.append("")

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def report_and_exit(args, errors, warnings, stats):
    write_report(args.output_report, errors, warnings, stats)
    for message in errors:
        print(f"ERROR {message}", file=sys.stderr)
    for message in warnings:
        print(f"WARNING {message}", file=sys.stderr)
    print(f"\nValidation failed. Full report: {args.output_report}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
