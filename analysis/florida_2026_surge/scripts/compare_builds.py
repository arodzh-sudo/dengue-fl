#!/usr/bin/env python3
"""Show what moved between two builds, which is the list of sentences to revisit.

Usage: compare_builds.py <older build dir> <newer build dir>
"""

import csv
import os
import sys


def facts(build):
    path = os.path.join(build, "results", "report_facts.tsv")
    if not os.path.exists(path):
        raise SystemExit("{} has no results/report_facts.tsv, run analyze_v2.py on it".format(build))
    with open(path, newline="", encoding="utf-8") as handle:
        return {row["fact"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__.strip())
    older, newer = facts(sys.argv[1]), facts(sys.argv[2])
    width = max(len(name) for name in set(older) | set(newer))
    changed = 0
    print("{:<{}}  {:<22} {}".format("fact", width, os.path.basename(sys.argv[1].rstrip("/\\")),
                                     os.path.basename(sys.argv[2].rstrip("/\\"))))
    for name in list(older) + [n for n in newer if n not in older]:
        was, now = older.get(name, "absent"), newer.get(name, "absent")
        if was != now:
            changed += 1
            print("{:<{}}  {:<22} {}".format(name, width, was, now))
    print()
    print("{} of {} facts moved".format(changed, len(set(older) | set(newer))))
    if changed:
        print("every sentence quoting one of those needs rereading before the report goes out")


if __name__ == "__main__":
    main()
