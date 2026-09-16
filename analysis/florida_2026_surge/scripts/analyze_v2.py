#!/usr/bin/env python3
"""Gather every number the Florida dengue report needs from the 2026 build."""

import argparse
import csv
import datetime as dt
import json
import os
import re
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "..", "2026")
SEROTYPES = ["denv1", "denv2", "denv3", "denv4"]


def set_base(path):
    """Point every script at one build directory, holding json/ and the run reports."""
    global BASE
    BASE = os.path.abspath(path)


def build_dir_argument(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--build-dir", default=BASE,
                        help="folder holding json/, summary_report.txt, input_report.txt, "
                             "validation_report.txt and metadata.txt")
    arguments = parser.parse_args()
    set_base(arguments.build_dir)
    return arguments


def jsons():
    return os.path.join(BASE, "json")


def results():
    return os.path.join(BASE, "results")


def load(serotype):
    with open(os.path.join(jsons(), "dengue_{}_genome.json".format(serotype)), encoding="utf-8") as handle:
        return json.load(handle)


def walk(node, parent=None, parents=None):
    if parents is not None:
        parents[node["name"]] = parent
    yield node
    for child in node.get("children", []):
        yield from walk(child, node, parents)


def attr(node, key):
    return node.get("node_attrs", {}).get(key, {}).get("value")


def is_tip(node):
    return not node.get("children")


def decimal_to_date(value):
    if value is None:
        return ""
    year = int(value)
    start = dt.date(year, 1, 1).toordinal()
    length = dt.date(year + 1, 1, 1).toordinal() - start
    return dt.date.fromordinal(start + round((value - year) * length)).isoformat()


def parse_date(text):
    for form in ("%m/%d/%Y", "%Y-%m-%d", "%m/%d/%y"):
        try:
            return dt.datetime.strptime(text.strip(), form).date()
        except ValueError:
            continue
    return None


def read_tsv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(name, columns, rows):
    with open(os.path.join(results(), name), "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})


def sequences(tree, root_sequence, wanted):
    """Nucleotide sequence of each wanted tip, rebuilt from the root and the branch mutations."""
    out = {}
    seq = list(root_sequence)

    def descend(node):
        changed = []
        for mutation in node.get("branch_attrs", {}).get("mutations", {}).get("nuc", []):
            match = re.fullmatch(r"([A-Z-])(\d+)([A-Z-])", mutation)
            if match:
                site = int(match.group(2)) - 1
                changed.append((site, seq[site]))
                seq[site] = match.group(3)
        if node["name"] in wanted:
            out[node["name"]] = "".join(seq)
        for child in node.get("children", []):
            descend(child)
        for site, base in reversed(changed):
            seq[site] = base

    descend(tree)
    return out


def distance(one, two):
    """Differing sites and sites compared, ignoring anything that is not an unambiguous base."""
    calls = set("ACGT")
    differ = shared = 0
    for a, b in zip(one, two):
        if a in calls and b in calls:
            shared += 1
            if a != b:
                differ += 1
    return differ, shared


def mrca(parents, names):
    lineages = []
    for name in names:
        chain, cursor = [], name
        while cursor is not None:
            chain.append(cursor)
            cursor = parents[cursor]["name"] if parents[cursor] else None
        lineages.append(chain)
    common = set(lineages[0])
    for chain in lineages[1:]:
        common &= set(chain)
    for name in lineages[0]:
        if name in common:
            return name
    return None


def main():
    os.makedirs(results(), exist_ok=True)
    findings = []

    def say(text=""):
        findings.append(text)
        print(text)

    summary = read_tsv(os.path.join(BASE, "summary_report.txt"))
    metadata = {row["sample_id"]: row for row in read_tsv(os.path.join(BASE, "metadata.txt"))}
    validation = {}
    with open(os.path.join(BASE, "validation_report.txt"), encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) == 5 and parts[1].startswith("denv"):
                validation[parts[0]] = {"length": parts[2], "unambiguous": parts[3],
                                        "ambiguous_fraction": parts[4].strip()}

    def collected(name, node):
        """Collection date as the laboratory recorded it, since the tree's decimal dates round."""
        recorded = parse_date(metadata.get(name, {}).get("collection_date", ""))
        return recorded.isoformat() if recorded else decimal_to_date(attr(node, "num_date"))

    trees, parents_by, tips_by = {}, {}, {}
    for serotype in SEROTYPES:
        data = load(serotype)
        parents = {}
        nodes = list(walk(data["tree"], None, parents))
        trees[serotype] = data
        parents_by[serotype] = parents
        tips_by[serotype] = [n for n in nodes if is_tip(n)]

    # 1. the funnel
    specimens = {}
    for row in summary:
        base_id = re.sub(r"_r$", "", row["sample_id"])
        specimens.setdefault(base_id, []).append(row)
    flags = Counter(rows[0]["vadr_flag"] for rows in specimens.values())
    serotyped = sum(1 for rows in specimens.values() if re.search(r"[1-4]", rows[0]["serotype"]))
    in_build = {}
    for serotype in SEROTYPES:
        for tip in tips_by[serotype]:
            if attr(tip, "data_source") == "Florida BPHL":
                in_build[tip["name"]] = serotype
    funnel = [
        {"step": "rows in summary_report.txt", "n": len(summary)},
        {"step": "specimens after removing the _r repeat rows", "n": len(specimens)},
        {"step": "serotype assigned", "n": serotyped},
        {"step": "VADR PASS", "n": flags["PASS"]},
        {"step": "VADR REVIEW", "n": flags["REVIEW"]},
        {"step": "VADR FAIL", "n": flags["FAIL"]},
        {"step": "genomes in the build", "n": len(in_build)},
    ]
    write_tsv("funnel.tsv", ["step", "n"], funnel)
    say("FUNNEL")
    for row in funnel:
        say("  {:>4}  {}".format(row["n"], row["step"]))

    # 2. every Florida tip
    rows = []
    for serotype in SEROTYPES:
        for tip in tips_by[serotype]:
            if attr(tip, "data_source") != "Florida BPHL":
                continue
            run = specimens.get(tip["name"], [{}])[0]
            rows.append({
                "sample": tip["name"],
                "serotype": serotype,
                "lineage": attr(tip, "minor_lineage"),
                "date": collected(tip["name"], tip),
                "county": attr(tip, "location"),
                "case_origin": attr(tip, "case_origin") or "(missing)",
                "country_exposure": attr(tip, "country_exposure"),
                "host_type": attr(tip, "host_type"),
                "vadr_flag": attr(tip, "vadr_flag"),
                "percent_assembled": run.get("percent_genome_cov_assembled", ""),
                "ambiguous_fraction": validation.get(tip["name"], {}).get("ambiguous_fraction", ""),
                "private_mutations": len(tip.get("branch_attrs", {}).get("mutations", {}).get("nuc", [])),
            })
    rows.sort(key=lambda r: (r["serotype"], r["date"]))
    write_tsv("tips.tsv", list(rows[0]), rows)
    say("")
    say("Florida tips written: {}".format(len(rows)))
    say("  missing case origin: " + ", ".join(r["sample"] for r in rows if r["case_origin"] == "(missing)"))

    # 3. the Hillsborough and Pinellas DENV2 clade
    denv2 = trees["denv2"]
    parents = parents_by["denv2"]
    nodes = {n["name"]: n for n in walk(denv2["tree"])}
    seed = [t["name"] for t in tips_by["denv2"]
            if attr(t, "data_source") == "Florida BPHL"
            and attr(t, "location") in ("Hillsborough", "Pinellas")
            and (attr(t, "case_origin") == "local" or attr(t, "host_type") == "Mosquito")]
    say("")
    say("CLUSTER SEED (local Hillsborough and Pinellas plus the pool): " + ", ".join(sorted(seed)))
    ancestor = mrca(parents, seed)
    clade = [n["name"] for n in walk(nodes[ancestor]) if is_tip(n)]
    confidence = nodes[ancestor]["node_attrs"]["num_date"]["confidence"]
    say("MRCA of the seed: {} {} CI {} to {}".format(
        ancestor, decimal_to_date(attr(nodes[ancestor], "num_date")),
        decimal_to_date(confidence[0]), decimal_to_date(confidence[1])))
    say("tips under it: {}".format(len(clade)))

    cluster_rows = []
    for name in clade:
        tip = nodes[name]
        cluster_rows.append({
            "sample": name,
            "source": attr(tip, "data_source"),
            "county": attr(tip, "location"),
            "country_exposure": attr(tip, "country_exposure"),
            "case_origin": attr(tip, "case_origin") or "(missing)",
            "host_type": attr(tip, "host_type"),
            "date": collected(name, tip),
            "vadr_flag": attr(tip, "vadr_flag"),
            "private_mutations": len(tip.get("branch_attrs", {}).get("mutations", {}).get("nuc", [])),
            "ambiguous_fraction": validation.get(name, {}).get("ambiguous_fraction", ""),
        })
    cluster_rows.sort(key=lambda r: r["date"])
    write_tsv("denv2_cluster.tsv", list(cluster_rows[0]), cluster_rows)
    for row in cluster_rows:
        say("  " + "\t".join(str(row[k]) for k in ("date", "sample", "county", "case_origin",
                                                   "vadr_flag", "private_mutations",
                                                   "ambiguous_fraction")))
    say("  VADR flags in the clade: {}".format(dict(Counter(r["vadr_flag"] for r in cluster_rows))))

    # 3b. how the clade is built inside
    say("")
    say("STRUCTURE INSIDE THE CLADE")
    inside_rows = []
    for node in walk(nodes[ancestor]):
        if is_tip(node):
            continue
        members = [n["name"] for n in walk(node) if is_tip(n)]
        muts = node.get("branch_attrs", {}).get("mutations", {})
        row = {
            "node": node["name"],
            "date": decimal_to_date(attr(node, "num_date")),
            "ci": "{} to {}".format(*[decimal_to_date(v) for v in node["node_attrs"]["num_date"]["confidence"]]),
            "n_tips": len(members),
            "counties": ",".join(sorted({attr(nodes[m], "location") for m in members})),
            "nuc_changes_on_branch": len(muts.get("nuc", [])),
            "aa_changes_on_branch": "; ".join("{} {}".format(g, ",".join(m))
                                              for g, m in muts.items() if g != "nuc"),
            "tips": " ".join(sorted(members)),
        }
        inside_rows.append(row)
        say("  {}\t{}\t{} tips\t{}\t{} nt\t{}".format(row["node"], row["date"], row["n_tips"],
                                                      row["counties"], row["nuc_changes_on_branch"],
                                                      row["aa_changes_on_branch"]))
        say("      " + row["tips"])
    write_tsv("denv2_cluster_structure.tsv", list(inside_rows[0]), inside_rows)

    # the path from the clade up to the root, with the reconstructed exposure at each step
    say("")
    say("PATH FROM THE CLUSTER ANCESTOR TO THE ROOT")
    cursor, path_rows = ancestor, []
    seen, outside = set(clade), []
    while cursor is not None:
        node = nodes[cursor]
        conf = node["node_attrs"].get("country_exposure", {}).get("confidence", {})
        top = sorted(conf.items(), key=lambda kv: -kv[1])[:3]
        muts = node.get("branch_attrs", {}).get("mutations", {})
        changes = "; ".join("{} {}".format(gene, ",".join(m)) for gene, m in muts.items() if gene != "nuc")
        row = {
            "node": cursor,
            "date": decimal_to_date(attr(node, "num_date")),
            "tips": sum(1 for n in walk(node) if is_tip(n)),
            "florida_tips": sum(1 for n in walk(node) if is_tip(n) and attr(n, "data_source") == "Florida BPHL"),
            "nuc_changes_on_branch": len(muts.get("nuc", [])),
            "aa_changes_on_branch": changes,
            "country_exposure": " ".join("{}:{:.2f}".format(k, v) for k, v in top),
        }
        path_rows.append(row)
        say("  " + "\t".join(str(row[k]) for k in ("node", "date", "tips", "florida_tips",
                                                   "nuc_changes_on_branch", "country_exposure")))
        if row["aa_changes_on_branch"]:
            say("      " + row["aa_changes_on_branch"])
        joining = [n["name"] for n in walk(node) if is_tip(n) and n["name"] not in seen]
        for name in joining:
            seen.add(name)
            tip = nodes[name]
            outside.append({"joins_at": cursor, "sample": name,
                            "source": attr(tip, "data_source"),
                            "county": attr(tip, "location"),
                            "country": attr(tip, "country"),
                            "country_exposure": attr(tip, "country_exposure"),
                            "case_origin": attr(tip, "case_origin") or "(missing)",
                            "date": collected(name, tip)})
            if len(outside) <= 25:
                say("      joins here: {} {} {} {} {}".format(
                    name, outside[-1]["date"], outside[-1]["county"] or outside[-1]["country"],
                    outside[-1]["country_exposure"], outside[-1]["case_origin"]))
        cursor = parents[cursor]["name"] if parents[cursor] else None
    write_tsv("denv2_cluster_path.tsv", list(path_rows[0]), path_rows)
    write_tsv("denv2_cluster_relatives.tsv", list(outside[0]), outside[:40])

    # 4. pairwise distances over the focus set
    focus = list(clade) + [row["sample"] for row in outside[:12]]
    for extra in ("TVU26000552", "TVU26000553", "TVU26000530"):
        if extra in nodes and extra not in focus:
            focus.append(extra)
    seqs = sequences(denv2["tree"], denv2["root_sequence"]["nuc"], set(focus))
    pairs = []
    for i, one in enumerate(focus):
        for two in focus[i + 1:]:
            differ, shared = distance(seqs[one], seqs[two])
            called = [int(validation.get(n, {}).get("unambiguous") or 0) for n in (one, two)]
            pairs.append({"a": one, "b": two, "snps": differ, "sites_compared": shared,
                          "unambiguous_a": called[0], "unambiguous_b": called[1],
                          "assembled_sites_in_common_at_most": min(called) if all(called) else ""})
    write_tsv("denv2_distances.tsv", ["a", "b", "snps", "sites_compared", "unambiguous_a",
                                      "unambiguous_b", "assembled_sites_in_common_at_most"], pairs)
    say("")
    say("PAIRWISE DISTANCES over {} genomes, {} pairs".format(len(focus), len(pairs)))
    for row in sorted(pairs, key=lambda r: r["snps"])[:12]:
        say("  {:>3} SNP  {} vs {}  ({} sites)".format(row["snps"], row["a"], row["b"], row["sites_compared"]))
    inside = set(clade)
    across = [p for p in pairs if (p["a"] in inside) != (p["b"] in inside)]
    within = [p for p in pairs if p["a"] in inside and p["b"] in inside]
    say("  within the clade: {} to {} SNP".format(min(p["snps"] for p in within),
                                                  max(p["snps"] for p in within)))
    closest = min(across, key=lambda p: p["snps"])
    say("  clade to anything outside it: closest {} SNP, {} vs {}".format(
        closest["snps"], closest["a"], closest["b"]))
    for target in ("TVU26000552", "MosquitoPool_K26-10948"):
        near = sorted((p for p in pairs if target in (p["a"], p["b"])), key=lambda r: r["snps"])[:5]
        say("  nearest to {}: ".format(target) + ", ".join(
            "{} {}".format(p["b"] if p["a"] == target else p["a"], p["snps"]) for p in near))

    # 4b. what sits nearest to each local case that is not in the cluster
    say("")
    say("NEAREST RELATIVES OF EVERY LOCAL CASE OUTSIDE THE HILLSBOROUGH CLADE")
    nearest_rows = []
    for serotype in SEROTYPES:
        queries = [t["name"] for t in tips_by[serotype]
                   if attr(t, "data_source") == "Florida BPHL"
                   and t["name"] not in clade
                   and (attr(t, "case_origin") == "local"
                        or attr(t, "host_type") == "Mosquito"
                        or attr(t, "case_origin") is None)]
        if not queries:
            continue
        tree = trees[serotype]
        every = {n["name"] for n in tips_by[serotype]}
        seqs_all = sequences(tree["tree"], tree["root_sequence"]["nuc"], every)
        index = {n["name"]: n for n in walk(tree["tree"])}
        for query in queries:
            scored = []
            for other in every:
                if other == query:
                    continue
                differ, shared = distance(seqs_all[query], seqs_all[other])
                scored.append((differ, other, shared))
            scored.sort()
            for differ, other, shared in scored[:5]:
                tip = index[other]
                row = {"serotype": serotype, "sample": query, "neighbour": other, "snps": differ,
                       "source": attr(tip, "data_source"),
                       "place": attr(tip, "location") or attr(tip, "country"),
                       "country_exposure": attr(tip, "country_exposure"),
                       "case_origin": attr(tip, "case_origin") or "",
                       "date": collected(other, tip)}
                nearest_rows.append(row)
            best = nearest_rows[-5:]
            say("  {} {} ({}, {}): ".format(serotype, query,
                                            attr(index[query], "location"),
                                            collected(query, index[query]))
                + ", ".join("{} {} SNP [{} {} {}]".format(r["neighbour"], r["snps"], r["source"],
                                                          r["place"], r["date"]) for r in best))
    write_tsv("nearest_relatives.tsv", list(nearest_rows[0]), nearest_rows)

    # 4d. the two clades behind the Florida DENV2 genomes
    def marked_by(gene, change):
        """The clade whose own branch carries this amino acid change."""
        for node in walk(denv2["tree"]):
            if change in node.get("branch_attrs", {}).get("mutations", {}).get(gene, []):
                return node
        raise SystemExit("no branch carries {} {}".format(gene, change))

    clades = [("outbreak", marked_by("E", "S7A")), ("travelers", marked_by("NS2A", "I33L"))]
    sublineage_rows, member_of, context_rows = [], {}, []
    say("")
    say("THE TWO CLADES")
    for label, node in clades:
        members = [n for n in walk(node) if is_tip(n) and attr(n, "data_source") == "Florida BPHL"]
        muts = node.get("branch_attrs", {}).get("mutations", {})
        say("  {} {} {} | {} Florida genomes | {} nt on its branch | {}".format(
            label, node["name"], decimal_to_date(attr(node, "num_date")), len(members),
            len(muts.get("nuc", [])),
            "; ".join("{} {}".format(g, ",".join(m)) for g, m in muts.items() if g != "nuc")))
        say("     counties: {}".format(dict(Counter(attr(t, "location") for t in members))))
        say("     origin: {}".format(dict(Counter(attr(t, "case_origin") or "(missing)" for t in members))))
        for tip in members:
            member_of[tip["name"]] = label
            sublineage_rows.append({
                "sample": tip["name"], "sublineage": label, "node": node["name"],
                "county": attr(tip, "location"), "country_exposure": attr(tip, "country_exposure"),
                "case_origin": attr(tip, "case_origin") or "(missing)",
                "date": collected(tip["name"], tip), "vadr_flag": attr(tip, "vadr_flag"),
            })
    sublineage_rows.sort(key=lambda r: (r["sublineage"], r["date"]))
    write_tsv("denv2_sublineages.tsv", list(sublineage_rows[0]), sublineage_rows)

    # what each clade is related to, level by level up the tree
    say("")
    say("WHAT EACH CLADE IS RELATED TO")
    for label, node in clades:
        seen = {n["name"] for n in walk(node) if is_tip(n)}
        inside = [n for n in walk(node) if is_tip(n) and attr(n, "data_source") != "Florida BPHL"]
        say("  {}: {} tips, public inside: {}".format(
            label, len(seen), ", ".join("{} {} {}".format(
                n["name"], attr(n, "country"), collected(n["name"], n)) for n in inside) or "none"))
        cursor = parents[node["name"]]
        for level in range(1, 4):
            if cursor is None:
                break
            joining = [n for n in walk(cursor) if is_tip(n) and n["name"] not in seen]
            for tip in joining:
                seen.add(tip["name"])
                context_rows.append({
                    "clade": label, "level": level, "node": cursor["name"],
                    "sample": tip["name"], "source": attr(tip, "data_source"),
                    "country": attr(tip, "country"), "division": attr(tip, "division"),
                    "date": collected(tip["name"], tip)})
            say("    level {} {} {}: {} tips joining, {}".format(
                level, cursor["name"], decimal_to_date(attr(cursor, "num_date")), len(joining),
                dict(Counter(attr(n, "country") for n in joining))))
            cursor = parents[cursor["name"]]
    write_tsv("denv2_clade_context.tsv",
              ["clade", "level", "node", "sample", "source", "country", "division", "date"],
              context_rows)

    # every close pair of Florida genomes outside the outbreak group
    florida = [t["name"] for t in tips_by["denv2"] if attr(t, "data_source") == "Florida BPHL"]
    all_seqs = sequences(denv2["tree"], denv2["root_sequence"]["nuc"], set(florida))
    close = []
    for i, one in enumerate(florida):
        for two in florida[i + 1:]:
            if one in inside and two in inside:
                continue
            differ, _ = distance(all_seqs[one], all_seqs[two])
            if differ > 3:
                continue
            close.append({
                "snps": differ, "a": one, "b": two,
                "county_a": attr(nodes[one], "location"), "county_b": attr(nodes[two], "location"),
                "exposure_a": attr(nodes[one], "country_exposure"),
                "exposure_b": attr(nodes[two], "country_exposure"),
                "date_a": collected(one, nodes[one]), "date_b": collected(two, nodes[two]),
                "sublineage": member_of.get(one, "") if member_of.get(one) == member_of.get(two) else "",
            })
    close.sort(key=lambda r: (r["snps"], r["date_a"]))
    write_tsv("denv2_close_pairs.tsv",
              ["snps", "a", "b", "county_a", "county_b", "exposure_a", "exposure_b",
               "date_a", "date_b", "sublineage"], close)
    say("")
    say("FLORIDA DENV2 PAIRS WITHIN 3 SNP, OUTSIDE THE OUTBREAK GROUP")
    for row in close:
        say("  {:>2} SNP  {} {} {} {}  vs  {} {} {} {}  [{}]".format(
            row["snps"], row["a"], row["county_a"], row["exposure_a"], row["date_a"],
            row["b"], row["county_b"], row["exposure_b"], row["date_b"], row["sublineage"]))

    # 4c. what each serotype contributes
    say("")
    say("BY SEROTYPE")
    for serotype in SEROTYPES:
        florida = [t for t in tips_by[serotype] if attr(t, "data_source") == "Florida BPHL"]
        say("  {}: {} Florida genomes, lineages {}, counties {}, origins {}".format(
            serotype, len(florida),
            dict(Counter(attr(t, "minor_lineage") for t in florida)),
            dict(Counter(attr(t, "location") for t in florida)),
            dict(Counter(attr(t, "case_origin") or "(missing)" for t in florida))))

    # 5. the lineage in context
    context = []
    for tip in tips_by["denv2"]:
        if attr(tip, "minor_lineage") == "2II_F.1.1.2":
            context.append({"country": attr(tip, "country"),
                            "year": int(attr(tip, "num_date")) if attr(tip, "num_date") else "",
                            "source": attr(tip, "data_source")})
    counts = Counter((c["country"], c["year"]) for c in context)
    write_tsv("lineage_2II_F_1_1_2.tsv", ["country", "year", "n"],
              [{"country": k[0], "year": k[1], "n": v} for k, v in sorted(counts.items())])
    say("")
    say("2II_F.1.1.2 BY COUNTRY: " + ", ".join(
        "{} {}".format(k, v) for k, v in Counter(c["country"] for c in context).most_common()))
    cuba = [int(attr(t, "num_date")) for t in tips_by["denv2"]
            if attr(t, "country") == "Cuba" and attr(t, "num_date")]
    say("Cuban DENV2 genomes in the build: {}, latest {}, of lineage 2II_F.1.1.2: {}".format(
        len(cuba), max(cuba), sum(1 for c in context if c["country"] == "Cuba")))

    # 6. epidemiological denominators
    weeks, months, counties, origins = Counter(), Counter(), Counter(), Counter()
    specimen_rows = []
    for sample, runs in sorted(specimens.items()):
        row = metadata.get(sample, {})
        date = parse_date(row.get("collection_date", ""))
        origin = row.get("case_origin") or "blank"
        serotype = in_build.get(sample) or runs[0]["serotype"].lower()
        week = (date - dt.timedelta(days=date.weekday())).isoformat() if date else ""
        specimen_rows.append({
            "sample": sample,
            "date": date.isoformat() if date else "",
            "week": week,
            "county": row.get("location", ""),
            "case_origin": origin,
            "travel_country": row.get("travel_country", ""),
            "serotype": serotype,
            "vadr_flag": runs[0]["vadr_flag"],
            "percent_assembled": runs[0].get("percent_genome_cov_assembled", ""),
            "in_build": "yes" if sample in in_build else "no",
        })
        if date and date.year == 2026:
            weeks[(week, origin)] += 1
            months[(date.strftime("%Y-%m"), origin)] += 1
            counties[(row.get("location", ""), origin)] += 1
            origins[origin] += 1
    write_tsv("sequenced_specimens.tsv", list(specimen_rows[0]), specimen_rows)
    write_tsv("sequenced_weeks_2026.tsv", ["week", "case_origin", "n"],
              [{"week": k[0], "case_origin": k[1], "n": v} for k, v in sorted(weeks.items())])
    write_tsv("sequenced_months_2026.tsv", ["month", "case_origin", "n"],
              [{"month": k[0], "case_origin": k[1], "n": v} for k, v in sorted(months.items())])
    write_tsv("sequenced_counties_2026.tsv", ["county", "case_origin", "n"],
              [{"county": k[0], "case_origin": k[1], "n": v} for k, v in sorted(counties.items())])
    dated = [r for r in specimen_rows if r["date"]]
    say("")
    say("SEQUENCED SPECIMENS: {}, with a collection date {}, in the build {}".format(
        len(specimen_rows), len(dated), sum(1 for r in specimen_rows if r["in_build"] == "yes")))
    say("  collected in 2026: {} ({})".format(sum(origins.values()), dict(origins)))
    say("  2026 local by county: " + ", ".join(
        "{} {}".format(k[0], v) for k, v in sorted(counties.items(), key=lambda kv: -kv[1]) if k[1] == "local"))
    say("  2026 by week: " + ", ".join("{} {} {}".format(k[0], k[1][:5], v) for k, v in sorted(weeks.items())))
    say("  date range: {} to {}".format(min(r["date"] for r in dated), max(r["date"] for r in dated)))
    say("  specimens with no metadata row: " + ", ".join(
        r["sample"] for r in specimen_rows if not r["date"]))

    with open(os.path.join(results(), "findings.txt"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(findings) + "\n")


if __name__ == "__main__":
    build_dir_argument(__doc__)
    main()
