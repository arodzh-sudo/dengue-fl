#!/usr/bin/env python3
"""Draw every figure in the Florida dengue report from the 2026 build."""

import csv
import datetime as dt
import os
from collections import Counter, OrderedDict

import analyze_v2
from analyze_v2 import (SEROTYPES, attr, build_dir_argument, decimal_to_date, distance, is_tip,
                        load, mrca, results, sequences, walk)


def figures():
    return os.path.join(analyze_v2.BASE, "report", "figures")

SURFACE = "#ffffff"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
INK_MUTED = "#8a8a85"
GRID = "#e4e3df"
HILLSBOROUGH = "#2a78d6"
PINELLAS = "#eb6834"
DADE = "#1baf7a"
OTHER_COUNTY = "#6f6e69"
PUBLIC = "#c9c8c2"
LOCAL = "#1c5cab"
IMPORTED = "#86b6ef"
UNDETERMINED = "#8a8a85"
WEDGE = "#4a3aa7"
BLUES = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
         "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
FONT = "Inter, 'Helvetica Neue', Arial, sans-serif"
MONO = "'SF Mono', Menlo, Consolas, monospace"

COUNTY_COLOR = {"Hillsborough": HILLSBOROUGH, "Pinellas": PINELLAS, "Dade": DADE}
ORIGIN_COLOR = {"local": LOCAL, "travel-associated": IMPORTED, "unknown": UNDETERMINED,
                 "(missing)": UNDETERMINED, "blank": UNDETERMINED}


class Svg:
    def __init__(self, width, height):
        self.width, self.height = width, height
        self.parts = []

    def add(self, markup):
        self.parts.append(markup)

    def rect(self, x, y, w, h, fill, **kw):
        self.add('<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}" fill="{}"{}/>'.format(
            x, y, max(w, 0), max(h, 0), fill, attrs(kw)))

    def line(self, x1, y1, x2, y2, stroke=GRID, width=1, **kw):
        self.add('<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="{}" '
                 'stroke-width="{}"{}/>'.format(x1, y1, x2, y2, stroke, width, attrs(kw)))

    def text(self, x, y, value, size=11, fill=INK, anchor="start", weight=400, family=FONT, **kw):
        self.add('<text x="{:.1f}" y="{:.1f}" font-family="{}" font-size="{}" fill="{}" '
                 'text-anchor="{}" font-weight="{}"{}>{}</text>'.format(
                     x, y, family, size, fill, anchor, weight, attrs(kw), escape(value)))

    def circle(self, x, y, r, fill, **kw):
        self.add('<circle cx="{:.1f}" cy="{:.1f}" r="{}" fill="{}"{}/>'.format(x, y, r, fill, attrs(kw)))

    def diamond(self, x, y, r, fill, **kw):
        points = " ".join("{:.1f},{:.1f}".format(px, py) for px, py in
                          ((x, y - r), (x + r, y), (x, y + r), (x - r, y)))
        self.add('<polygon points="{}" fill="{}"{}/>'.format(points, fill, attrs(kw)))

    def save(self, name):
        os.makedirs(figures(), exist_ok=True)
        head = ('<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
                'viewBox="0 0 {w} {h}" role="img">'.format(w=self.width, h=self.height))
        body = '<rect width="{}" height="{}" fill="{}"/>'.format(self.width, self.height, SURFACE)
        with open(os.path.join(figures(), name), "w", encoding="utf-8") as handle:
            handle.write(head + body + "".join(self.parts) + "</svg>\n")
        print("wrote", name)


def attrs(kw):
    return "".join(' {}="{}"'.format(k.replace("_", "-"), v) for k, v in kw.items())


def escape(value):
    return (str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def recorded_dates():
    """Collection dates as the laboratory recorded them, keyed by sample."""
    return {row["sample"]: row["date"] for row in read("tips.tsv")}


def read(name):
    with open(os.path.join(results(), name), newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def legend(svg, x, y, entries, size=9, gap=118, shape="circle"):
    for index, (label, color) in enumerate(entries):
        cx = x + index * gap
        if shape == "square":
            svg.rect(cx, y - 7, 9, 9, color, rx=2)
        else:
            svg.circle(cx + 4, y - 3, 4.5, color)
        svg.text(cx + 14, y, label, size=size, fill=INK_SOFT)


def year_fraction(date):
    start = dt.date(date.year, 1, 1).toordinal()
    length = dt.date(date.year + 1, 1, 1).toordinal() - start
    return date.year + (date.toordinal() - start) / length


def month_ticks(low, high):
    ticks, year, month = [], int(low), 1
    while year <= int(high) + 1:
        point = year_fraction(dt.date(year, month, 1))
        if low <= point <= high:
            ticks.append((point, dt.date(year, month, 1)))
        month += 1
        if month > 12:
            month, year = 1, year + 1
    return ticks


# ---------------------------------------------------------------- figure 1
def figure_one():
    weeks = read("sequenced_weeks_2026.tsv")
    svg = Svg(900, 300)
    base = -10
    counts = {}
    for row_ in weeks:
        counts.setdefault(row_["week"], Counter())[row_["case_origin"]] = int(row_["n"])
    ordered = sorted(counts)
    first = dt.date.fromisoformat(ordered[0])
    last = dt.date.fromisoformat(ordered[-1])
    span = [(first + dt.timedelta(days=7 * i)).isoformat()
            for i in range((last - first).days // 7 + 1)]
    plot_top, plot_bottom, left, right = base + 30, base + 210, 66, 880
    tallest = max(sum(c.values()) for c in counts.values())
    step = (right - left) / max(len(span), 1)
    for level in range(0, tallest + 1, 2):
        y = plot_bottom - level / tallest * (plot_bottom - plot_top)
        svg.line(left, y, right, y, stroke=GRID)
        svg.text(left - 8, y + 4, str(level), size=9, fill=INK_MUTED, anchor="end")
    for index, week in enumerate(span):
        bucket = counts.get(week, Counter())
        x = left + index * step
        stacked = 0
        for origin in ("travel-associated", "unknown", "local"):
            value = bucket.get(origin, 0)
            if not value:
                continue
            height = value / tallest * (plot_bottom - plot_top)
            y = plot_bottom - stacked - height
            svg.rect(x + 2, y, step - 6, height - 2, ORIGIN_COLOR[origin], rx=2)
            stacked += height
        if sum(bucket.values()):
            svg.text(x + (step - 4) / 2, plot_bottom - stacked - 6, str(sum(bucket.values())),
                     size=9, fill=INK_SOFT, anchor="middle")
        if dt.date.fromisoformat(week).day <= 7 or index == 0:
            svg.text(x + (step - 4) / 2, plot_bottom + 16,
                     dt.date.fromisoformat(week).strftime("%b %d"), size=9,
                     fill=INK_MUTED, anchor="middle")
    svg.line(left, plot_bottom, right, plot_bottom, stroke=INK_MUTED)
    svg.add('<text x="{:.1f}" y="{:.1f}" font-family="{}" font-size="11" fill="{}" '
            'text-anchor="middle" transform="rotate(-90 {:.1f} {:.1f})">{}</text>'.format(
                18, (plot_top + plot_bottom) / 2, FONT, INK_SOFT, 18,
                (plot_top + plot_bottom) / 2, "Sequenced specimens"))
    svg.text((left + right) / 2, plot_bottom + 42, "Week of collection, weeks starting Monday",
             size=11, fill=INK_SOFT, anchor="middle")
    legend(svg, 66, plot_bottom + 74,
           [("Acquired in Florida", LOCAL), ("Travel-associated", IMPORTED),
            ("Origin not recorded", UNDETERMINED)], shape="square", gap=150)
    svg.height = int(plot_bottom) + 92
    svg.save("fig1_weeks.svg")


# ---------------------------------------------------------------- trees
def ladder(node):
    """Order children so the smaller clade sits on top, which keeps the ladder readable."""
    children = node.get("children", [])
    for child in children:
        ladder(child)
    children.sort(key=lambda c: (sum(1 for n in walk(c) if is_tip(n)), attr(c, "num_date") or 0))


def layout(root, collapse=()):
    y_of, order = {}, []

    def assign(node):
        if is_tip(node) or node["name"] in collapse:
            order.append(node["name"])
            y_of[node["name"]] = len(order) - 1
            return y_of[node["name"]]
        values = [assign(child) for child in node["children"]]
        y_of[node["name"]] = (min(values) + max(values)) / 2
        return y_of[node["name"]]

    assign(root)
    return y_of, order


def draw_branches(svg, root, y_of, to_x, to_y, color_of, width_of, collapse=()):
    for node in walk(root):
        if not node.get("children") or node["name"] in collapse:
            continue
        if any(child["name"] not in y_of for child in node["children"]):
            continue
        x = to_x(attr(node, "num_date"))
        ys = [to_y(y_of[child["name"]]) for child in node["children"]]
        svg.line(x, min(ys), x, max(ys), stroke=GRID, width=1.5)
        for child in node["children"]:
            cx = to_x(attr(child, "num_date"))
            cy = to_y(y_of[child["name"]])
            svg.line(x, cy, cx, cy, stroke=color_of(child), width=width_of(child))


def figure_two():
    data = load("denv2")
    parents = {}
    nodes = list(walk(data["tree"], None, parents))
    index = {n["name"]: n for n in nodes}
    florida = [n["name"] for n in nodes if is_tip(n) and attr(n, "data_source") == "Florida BPHL"]
    root = index[mrca(parents, florida)]
    ladder(root)
    y_of, order = layout(root)

    top, bottom, left, right = 40, 40 + len(order) * 6.3, 60, 700
    svg = Svg(1040, int(bottom) + 76)

    low = attr(root, "num_date")
    high = max(attr(index[name], "num_date") for name in order)
    pad = (high - low) * 0.02

    def to_x(value):
        return left + (value - low) / (high - low + pad) * (right - left)

    def to_y(slot):
        return top + slot * 6.3

    def color_of(node):
        if attr(node, "data_source") != "Florida BPHL":
            return PUBLIC
        return COUNTY_COLOR.get(attr(node, "location"), OTHER_COUNTY)

    for point, when in month_ticks(low, high):
        if when.month == 1:
            svg.line(to_x(point), top - 10, to_x(point), bottom + 6, stroke=GRID)
            svg.text(to_x(point), top - 16, str(when.year), size=10, fill=INK_MUTED, anchor="middle")

    draw_branches(svg, root, y_of, to_x, to_y, color_of,
                  lambda n: 1.8 if attr(n, "data_source") == "Florida BPHL" else 0.8)

    for name in order:
        node = index[name]
        x, y = to_x(attr(node, "num_date")), to_y(y_of[name])
        local = attr(node, "data_source") == "Florida BPHL"
        if not local:
            svg.circle(x, y, 1.8, PUBLIC)
            continue
        color = color_of(node)
        if attr(node, "host_type") == "Mosquito":
            svg.diamond(x, y, 4, color, stroke=SURFACE, stroke_width=1)
        elif attr(node, "case_origin") == "local":
            svg.circle(x, y, 3.4, color, stroke=SURFACE, stroke_width=1)
        else:
            svg.circle(x, y, 3.4, SURFACE, stroke=color, stroke_width=1.6)
        svg.text(x + 7, y + 3.2, name, size=7.4, fill=INK, family=MONO)

    cluster = [n["name"] for n in walk(index["NODE_0003278"]) if is_tip(n)]
    top_y = min(to_y(y_of[n]) for n in cluster) - 5
    bottom_y = max(to_y(y_of[n]) for n in cluster) + 5
    svg.line(866, top_y, 866, bottom_y, stroke=HILLSBOROUGH, width=2)
    svg.text(874, (top_y + bottom_y) / 2 - 4, "Hillsborough and", size=10, fill=HILLSBOROUGH, weight=600)
    svg.text(874, (top_y + bottom_y) / 2 + 9, "Pinellas group", size=10, fill=HILLSBOROUGH, weight=600)

    legend(svg, 0, bottom + 40, [("Hillsborough", HILLSBOROUGH), ("Pinellas", PINELLAS),
                                 ("Miami-Dade", DADE), ("Other county", OTHER_COUNTY),
                                 ("GenBank", PUBLIC)], gap=120)
    svg.text(0, bottom + 62,
             "Filled circle: acquired in Florida.  Open circle: travel-associated.  "
             "Diamond: mosquito pool.", size=10, fill=INK_MUTED)
    svg.save("fig2_denv2_florida_tree.svg")


def tree_figure(root, filename, collapse=(), band=None, collapse_label="",
                row=26, right=520, scale=1.0):
    """A time-scaled tree with county colors, SNP counts on the branches and one text column."""
    dates = recorded_dates()
    parents = {}
    index = {n["name"]: n for n in walk(root, None, parents)}
    ladder(root)
    y_of, order = layout(root, collapse)

    top, left = 54, 70
    bottom = top + (len(order) - 1) * row
    svg = Svg(1120, int(bottom) + 120)
    tips = [attr(index[n], "num_date") for n in order]
    # anything older than the tips themselves is drawn at the left edge, so one deep
    # relative cannot squeeze every recent branch into a sliver
    low = max(attr(root, "num_date") - 0.1, min(tips) - 1.2)
    high = max(tips) + 0.05

    def to_x(value):
        return max(left, left + (value - low) / (high - low) * (right - left))

    def to_y(slot):
        return top + slot * row

    span = high - low
    months = (1, 4, 7, 10) if span <= 2.2 else (1, 7) if span <= 4 else (1,)
    for point, when in month_ticks(low, high):
        if when.month in months:
            x = to_x(point)
            svg.line(x, top - 14, x, bottom + 14, stroke=GRID)
            svg.text(x, top - 20, when.strftime("%b %Y") if len(months) > 1 else str(when.year),
                     size=9 * scale, fill=INK_MUTED, anchor="middle")

    if band is not None:
        confidence = band["node_attrs"]["num_date"]["confidence"]
        svg.rect(to_x(confidence[0]), top - 16, to_x(confidence[1]) - to_x(confidence[0]),
                 bottom - top + 30, HILLSBOROUGH, opacity="0.10", rx=3)
        svg.text((to_x(confidence[0]) + to_x(confidence[1])) / 2, bottom + 32,
                 "common ancestor of the cluster {}, range {} to {}".format(
                     decimal_to_date(attr(band, "num_date")),
                     decimal_to_date(confidence[0]), decimal_to_date(confidence[1])),
                 size=9, fill=HILLSBOROUGH, anchor="middle")

    def color_of(node):
        if attr(node, "data_source") != "Florida BPHL":
            return PUBLIC
        return COUNTY_COLOR.get(attr(node, "location"), OTHER_COUNTY)

    draw_branches(svg, root, y_of, to_x, to_y, color_of, lambda n: 2.2, collapse)

    for node in walk(root):
        if node["name"] not in y_of and node["name"] not in collapse:
            continue
        muts = node.get("branch_attrs", {}).get("mutations", {})
        count = len(muts.get("nuc", []))
        changes = "; ".join("{}:{}".format(gene, ",".join(m)) for gene, m in muts.items()
                            if gene != "nuc")
        parent = parents[node["name"]]
        if not count or parent is None or parent["name"] not in y_of:
            continue
        x = (to_x(attr(parent, "num_date")) + to_x(attr(node, "num_date"))) / 2
        y = to_y(y_of[node["name"]])
        listed = sum(len(m) for gene, m in muts.items() if gene != "nuc")
        label = "{} [{}]".format(changes, count) if changes and listed <= 3 else "[{}]".format(count)
        svg.text(x, y - 6, label, size=8.5 * scale, fill=INK_SOFT if listed <= 3 and changes else INK_MUTED,
                 anchor="middle", family=MONO)

    for name in order:
        node = index[name]
        x, y = to_x(attr(node, "num_date")), to_y(y_of[name])
        color = color_of(node)
        if name in collapse:
            members = [n for n in walk(node) if is_tip(n)]
            far = to_x(max(attr(n, "num_date") for n in members))
            svg.add('<polygon points="{:.1f},{:.1f} {:.1f},{:.1f} {:.1f},{:.1f}" fill="{}" '
                    'opacity="0.8"/>'.format(x, y, far, y - 9, far, y + 9, WEDGE))
            svg.line(far + 8, y, 552, y, stroke=GRID, width=0.8, stroke_dasharray="1 3")
            svg.text(560, y + 4, collapse_label, size=10 * scale, weight=700, fill=INK)
            svg.text(760, y + 4, "{} genomes".format(len(members)), size=9 * scale, fill=INK_MUTED)
            continue
        if attr(node, "host_type") == "Mosquito":
            svg.diamond(x, y, 5.5, color, stroke=SURFACE, stroke_width=1.2)
        elif attr(node, "case_origin") == "local":
            svg.circle(x, y, 4.5, color, stroke=SURFACE, stroke_width=1.2)
        else:
            svg.circle(x, y, 4.5, SURFACE, stroke=color, stroke_width=2)
        svg.line(x + 8, y, 552, y, stroke=GRID, width=0.8, stroke_dasharray="1 3")
        local = attr(node, "data_source") == "Florida BPHL"
        svg.text(560, y + 4, name, size=10 * scale, family=MONO, fill=INK)
        when = dates.get(name) or decimal_to_date(attr(node, "num_date"))
        if not local:
            note = "{}  {}".format(when, attr(node, "country"))
        elif attr(node, "case_origin") == "travel-associated":
            note = "{}  travel to {}".format(when, attr(node, "country_exposure"))
        else:
            note = when
        svg.text(760, y + 4, note, size=9 * scale, fill=INK_MUTED)

    used = {color_of(index[n]) for n in order if n not in collapse}
    entries = [("Hillsborough", HILLSBOROUGH), ("Pinellas", PINELLAS), ("Miami-Dade", DADE),
               ("Other county", OTHER_COUNTY), ("GenBank", PUBLIC)]
    entries = [e for e in entries if e[1] in used]
    if collapse:
        entries.append((collapse_label, WEDGE))
    legend(svg, 0, bottom + 54, entries, gap=150)
    note = ("Filled: acquired in Florida.  Open: travel-associated.  Diamond: mosquito pool.  "
            "Branch labels: amino acid changes and, in brackets, SNPs.")
    if band is not None:
        note += "  Shaded band: the range of dates for the cluster's common ancestor."
    svg.text(0, bottom + 76, note, size=10, fill=INK_MUTED)
    svg.save(filename)


def denv2_tree(filename, root_marker=None, root_name=None, levels_up=0, collapse_marker=None,
               band_name=None, collapse_label="", row=26, scale=1.0):
    data = load("denv2")
    parents = {}
    nodes = list(walk(data["tree"], None, parents))
    index = {n["name"]: n for n in nodes}
    root = index[root_name] if root_name else marked_by(data["tree"], *root_marker)
    for _ in range(levels_up):
        root = parents[root["name"]]
    collapse = set()
    if collapse_marker:
        collapse.add(marked_by(root, *collapse_marker)["name"])
    tree_figure(root, filename, collapse, index[band_name] if band_name else None,
                collapse_label, row, scale=scale)


def prune(node, keep, carried=None):
    """A copy of the tree holding only the kept tips, with dropped branches folded into the next one."""
    carried = carried or {}
    muts = node.get("branch_attrs", {}).get("mutations", {})
    merged = {gene: list(carried.get(gene, [])) + list(values) for gene, values in muts.items()}
    for gene, values in carried.items():
        merged.setdefault(gene, list(values))
    if is_tip(node):
        if node["name"] not in keep:
            return None
        return {"name": node["name"], "node_attrs": node["node_attrs"],
                "branch_attrs": {"mutations": merged}}
    kept = [child for child in (prune(c, keep) for c in node["children"]) if child]
    if not kept:
        return None
    if len(kept) == 1:
        only = kept[0]
        for gene, values in merged.items():
            only["branch_attrs"]["mutations"][gene] = values + only["branch_attrs"]["mutations"].get(gene, [])
        return only
    return {"name": node["name"], "node_attrs": node["node_attrs"],
            "branch_attrs": {"mutations": merged}, "children": kept}


def nearest_public(name, parents, limit=2):
    """The public genomes in the smallest clade above this tip that holds any."""
    cursor = parents[name]
    while cursor is not None:
        public = [n["name"] for n in walk(cursor)
                  if is_tip(n) and attr(n, "data_source") != "Florida BPHL"]
        if public:
            return public[:limit]
        cursor = parents[cursor["name"]]
    return []


def serotype_tree(serotype, filename, row=26, scale=1.0):
    """Every Florida genome of one serotype with the public genomes nearest to each of them."""
    data = load(serotype)
    parents = {}
    nodes = list(walk(data["tree"], None, parents))
    florida = [n["name"] for n in nodes if is_tip(n) and attr(n, "data_source") == "Florida BPHL"]
    keep = set(florida)
    for name in florida:
        keep.update(nearest_public(name, parents))
    root = prune(data["tree"], keep)
    tree_figure(root, filename, row=row, scale=scale)


# ---------------------------------------------------------------- figure 4
def figure_four():
    rows = read("denv2_distances.tsv")
    cluster = [r["sample"] for r in read("denv2_cluster.tsv")]
    order = cluster + ["TVU26000552", "TVU26000553", "TVU26000530"]
    tips = read("tips.tsv")
    ambiguity = {r["sample"]: float(r["ambiguous_fraction"] or 0) for r in tips}
    county_of = {r["sample"]: r["county"] for r in tips}
    pairs = {}
    for row in rows:
        pairs[(row["a"], row["b"])] = int(row["snps"])
        pairs[(row["b"], row["a"])] = int(row["snps"])

    cell, left, top = 30, 190, 120
    size = len(order)
    svg = Svg(900, top + size * cell + 70)

    biggest = max(v for k, v in pairs.items() if k[0] in order and k[1] in order)
    for column, name in enumerate(order):
        x = left + column * cell
        svg.add('<text x="{:.1f}" y="{:.1f}" font-family="{}" font-size="9" fill="{}" '
                'text-anchor="start" transform="rotate(-60 {:.1f} {:.1f})">{}</text>'.format(
                    x + cell / 2, top - 8, MONO, INK_SOFT, x + cell / 2, top - 8, name))
    for row_index, name in enumerate(order):
        y = top + row_index * cell
        color = COUNTY_COLOR.get(county_of.get(name, ""), OTHER_COUNTY)
        svg.circle(left - 12, y + cell / 2, 4, color)
        svg.text(left - 20, y + cell / 2 + 3.5, name, size=9, family=MONO, fill=INK, anchor="end")
        for column, other in enumerate(order):
            x = left + column * cell
            if name == other:
                svg.rect(x + 1, y + 1, cell - 2, cell - 2, "#f0efec", rx=3)
                continue
            value = pairs.get((name, other))
            step = min(len(BLUES) - 1, int(value / biggest * (len(BLUES) - 1)))
            shade = BLUES[len(BLUES) - 1 - step]
            svg.rect(x + 1, y + 1, cell - 2, cell - 2, shade, rx=3)
            svg.text(x + cell / 2, y + cell / 2 + 3.5, str(value), size=8.5,
                     fill="#ffffff" if step < (len(BLUES) - 1) * 0.55 else INK,
                     anchor="middle", family=MONO)
    right = left + size * cell + 24
    svg.text(right, top - 8, "unresolved bases", size=9, fill=INK_SOFT)
    for row_index, name in enumerate(order):
        y = top + row_index * cell
        share = ambiguity.get(name, 0)
        svg.rect(right, y + cell / 2 - 5, max(share * 600, 1), 10, "#eb6834", rx=2)
        svg.text(right + max(share * 600, 1) + 6, y + cell / 2 + 3.5,
                 "{:.1f}%".format(share * 100), size=8.5, fill=INK_MUTED, family=MONO)
    counties = [("Hillsborough", HILLSBOROUGH), ("Pinellas", PINELLAS), ("Miami-Dade", DADE)]
    used = {COUNTY_COLOR.get(county_of.get(name, ""), OTHER_COUNTY) for name in order}
    legend(svg, 0, top + size * cell + 34, [c for c in counties if c[1] in used], gap=120)
    svg.text(0, top + size * cell + 58,
             "The dot beside each name gives the county. Darker cells are more similar. The bars "
             "on the right are the share of each assembly that was unresolved before the tree "
             "filled it in.", size=10, fill=INK_MUTED)
    svg.save("fig4_distances.svg")


# ---------------------------------------------------------------- figure 5
def figure_five():
    rows = read("lineage_2II_F_1_1_2.tsv")
    years = sorted({int(r["year"]) for r in rows})
    countries = [c for c, _ in Counter({r["country"]: 0 for r in rows}).items()]
    totals = Counter()
    grid = {}
    for row in rows:
        grid[(row["country"], int(row["year"]))] = int(row["n"])
        totals[row["country"]] += int(row["n"])
    countries = [c for c, _ in totals.most_common()]

    cell_w, cell_h, left, top = 64, 30, 150, 46
    svg = Svg(900, top + len(countries) * cell_h + 70)
    biggest = max(grid.values())
    for column, year in enumerate(years):
        svg.text(left + column * cell_w + cell_w / 2, top - 10, str(year), size=10,
                 fill=INK_SOFT, anchor="middle")
    for row_index, country in enumerate(countries):
        y = top + row_index * cell_h
        svg.text(left - 12, y + cell_h / 2 + 4, country, size=11, fill=INK, anchor="end")
        for column, year in enumerate(years):
            x = left + column * cell_w
            value = grid.get((country, year))
            if not value:
                svg.rect(x + 1, y + 1, cell_w - 3, cell_h - 3, "#f4f3f0", rx=3)
                continue
            shade = BLUES[min(len(BLUES) - 1, int((value / biggest) ** 0.45 * (len(BLUES) - 1)))]
            svg.rect(x + 1, y + 1, cell_w - 3, cell_h - 3, shade, rx=3)
            svg.text(x + cell_w / 2, y + cell_h / 2 + 4, str(value), size=10, family=MONO,
                     fill="#ffffff" if value > biggest * 0.3 else INK, anchor="middle")
        svg.text(left + len(years) * cell_w + 12, y + cell_h / 2 + 4, str(totals[country]),
                 size=10, weight=700, fill=INK_SOFT, family=MONO)
    svg.text(left + len(years) * cell_w + 12, top - 10, "total", size=10, fill=INK_SOFT)
    svg.text(150, top + len(countries) * cell_h + 34,
             "USA includes the 54 Florida genomes of this build.", size=10, fill=INK_MUTED)
    svg.save("fig5_lineage_context.svg")


# ---------------------------------------------------------------- figure 6
def figure_six():
    tips = read("tips.tsv")
    svg = Svg(900, 420)

    dates = [dt.date.fromisoformat(r["date"]) for r in tips]
    low, high = year_fraction(min(dates)) - 0.02, year_fraction(max(dates)) + 0.02
    left, right, top = 150, 870, 62

    def to_x(value):
        return left + (value - low) / (high - low) * (right - left)

    placements, heights = {}, {}
    for serotype in SEROTYPES:
        rows = sorted((r for r in tips if r["serotype"] == serotype), key=lambda r: r["date"])
        placed = []
        for record in rows:
            x = to_x(year_fraction(dt.date.fromisoformat(record["date"])))
            level = 0
            while any(abs(x - px) < 9 and level == pl for px, pl in placed):
                level += 1
            placed.append((x, level))
            placements.setdefault(serotype, []).append((record, x, level))
        heights[serotype] = max((p[1] for p in placed), default=0)

    bands, cursor = {}, top
    for serotype in SEROTYPES:
        bands[serotype] = cursor
        cursor += 46 + heights[serotype] * 10
    svg.height = int(cursor) + 70

    for point, when in month_ticks(low, high):
        x = to_x(point)
        svg.line(x, top - 16, x, cursor - 30, stroke=GRID)
        svg.text(x, top - 22, when.strftime("%b"), size=9, fill=INK_MUTED, anchor="middle")
        if when.month == 1:
            svg.text(x, top - 38, str(when.year), size=10, fill=INK_SOFT, anchor="middle")

    for serotype in SEROTYPES:
        y = bands[serotype]
        rows = [r for r in tips if r["serotype"] == serotype]
        svg.text(0, y + 4, serotype.upper(), size=12, weight=700)
        svg.text(0, y + 20, "{} genomes, {}".format(
            len(rows), ", ".join(sorted(Counter(r["lineage"] for r in rows)))),
            size=9, fill=INK_MUTED)
        for record, x, level in placements.get(serotype, []):
            color = COUNTY_COLOR.get(record["county"], OTHER_COUNTY)
            cy = y + level * 10
            if record["host_type"] == "Mosquito":
                svg.diamond(x, cy, 4.5, color, stroke=SURFACE, stroke_width=1)
            elif record["case_origin"] == "local":
                svg.circle(x, cy, 3.8, color, stroke=SURFACE, stroke_width=1)
            else:
                svg.circle(x, cy, 3.8, SURFACE, stroke=color, stroke_width=1.5)

    legend(svg, 0, cursor - 6, [("Hillsborough", HILLSBOROUGH), ("Pinellas", PINELLAS),
                                ("Miami-Dade", DADE), ("Other county", OTHER_COUNTY)], gap=120)
    svg.text(0, cursor + 18,
             "Filled: acquired in Florida.  Open: travel-associated.  Diamond: mosquito pool.",
             size=10, fill=INK_MUTED)
    svg.save("fig6_serotypes.svg")


CUBA = "#eb6834"
FLORIDA = "#1c5cab"
ELSEWHERE = "#8a8a85"


def exposure_color(node):
    if attr(node, "data_source") != "Florida BPHL":
        return PUBLIC
    if attr(node, "case_origin") is None and attr(node, "country_exposure") == "USA":
        return ELSEWHERE
    return {"Cuba": CUBA, "USA": FLORIDA}.get(attr(node, "country_exposure"), ELSEWHERE)


def marked_by(tree, gene, change):
    """The clade whose own branch carries this amino acid change."""
    for node in walk(tree):
        if change in node.get("branch_attrs", {}).get("mutations", {}).get(gene, []):
            return node
    raise SystemExit("no branch carries {} {}".format(gene, change))


def tip_mark(svg, node, x, y, color, size=3.6):
    if attr(node, "host_type") == "Mosquito":
        svg.diamond(x, y, size + 1.4, color, stroke=SURFACE, stroke_width=1)
    elif attr(node, "case_origin") == "local":
        svg.circle(x, y, size, color, stroke=SURFACE, stroke_width=1)
    else:
        svg.circle(x, y, size, SURFACE, stroke=color, stroke_width=1.6)


def exposure_legend(svg, x, y, used=None):
    entries = [("Acquired in Cuba", CUBA), ("Acquired in Florida", FLORIDA),
               ("Elsewhere or not recorded", ELSEWHERE), ("GenBank", PUBLIC)]
    if used is not None:
        entries = [e for e in entries if e[1] in used]
    legend(svg, x, y, entries, gap=136)
    svg.text(x, y + 22, "Filled: acquired in Florida.  Open: travel-associated.  "
                        "Diamond: mosquito pool.", size=10, fill=INK_MUTED)


# ---------------------------------------------------------------- figure 5a
def figure_five_a():
    data = load("denv2")
    roots = [("Sub-lineage carrying the outbreak", marked_by(data["tree"], "E", "S7A")),
             ("Imported-only sub-lineage", marked_by(data["tree"], "NS2A", "I33L"))]
    for _, root in roots:
        ladder(root)
    layouts = [(label, root) + layout(root) for label, root in roots]

    row, left, right, top, gap = 8.0, 40, 700, 56, 34
    low = min(attr(root, "num_date") for _, root in roots) - 0.3
    high = max(attr(index, "num_date") for _, root, y_of, order in layouts
               for index in [max((n for n in walk(root) if is_tip(n)),
                                 key=lambda n: attr(n, "num_date"))]) + 0.1

    def to_x(value):
        return left + (value - low) / (high - low) * (right - left)

    height = top
    bands = []
    for label, root, y_of, order in layouts:
        bands.append((label, root, y_of, order, height))
        height += len(order) * row + gap
    svg = Svg(1000, int(height) + 76)

    for point, when in month_ticks(low, high):
        if when.month == 1:
            x = to_x(point)
            svg.line(x, top - 12, x, height - gap, stroke=GRID)
            svg.text(x, top - 18, str(when.year), size=10, fill=INK_MUTED, anchor="middle")

    for label, root, y_of, order, offset in bands:
        def to_y(slot, offset=offset):
            return offset + slot * row

        draw_branches(svg, root, y_of, to_x, to_y, exposure_color,
                      lambda n: 1.8 if attr(n, "data_source") == "Florida BPHL" else 0.8)
        for name in order:
            node = [n for n in walk(root) if n["name"] == name][0]
            tip_mark(svg, node, to_x(attr(node, "num_date")), to_y(y_of[name]),
                     exposure_color(node), size=3.4)
        first, last = to_y(0) - 6, to_y(len(order) - 1) + 6
        svg.line(right + 30, first, right + 30, last, stroke=INK_MUTED, width=2)
        svg.text(right + 40, (first + last) / 2 - 5, label, size=11, weight=600, fill=INK)
        florida = sum(1 for n in walk(root) if is_tip(n) and attr(n, "data_source") == "Florida BPHL")
        svg.text(right + 40, (first + last) / 2 + 10,
                 "{} Florida genomes".format(florida), size=10, fill=INK_MUTED)

    exposure_legend(svg, 0, height + 20)
    svg.save("fig5a_exposure_tree.svg")


# ---------------------------------------------------------------- supplementary
def figure_supplementary():
    data = load("denv2")
    parents = {}
    nodes = list(walk(data["tree"], None, parents))
    index = {n["name"]: n for n in nodes}
    florida = [n["name"] for n in nodes if is_tip(n) and attr(n, "data_source") == "Florida BPHL"]
    root = index[mrca(parents, florida)]
    ladder(root)
    y_of, order = layout(root)

    top, bottom, left, right = 46, 46 + len(order) * 7.0, 60, 860
    svg = Svg(1000, int(bottom) + 76)
    low = attr(root, "num_date")
    high = max(attr(index[n], "num_date") for n in order)

    def to_x(value):
        return left + (value - low) / (high - low + 0.2) * (right - left)

    def to_y(slot):
        return top + slot * 6.3

    for point, when in month_ticks(low, high):
        if when.month == 1:
            svg.line(to_x(point), top - 12, to_x(point), bottom + 6, stroke=GRID)
            svg.text(to_x(point), top - 18, str(when.year), size=10, fill=INK_MUTED, anchor="middle")

    draw_branches(svg, root, y_of, to_x, to_y, exposure_color,
                  lambda n: 1.6 if attr(n, "data_source") == "Florida BPHL" else 0.8)
    for name in order:
        node = index[name]
        tip_mark(svg, node, to_x(attr(node, "num_date")), to_y(y_of[name]),
                 exposure_color(node), size=3.0)
    exposure_legend(svg, 0, bottom + 32)
    svg.save("supplementary_lineage_exposure_tree.svg")


def main():
    figure_one()
    figure_two()
    figure_four()
    figure_five_a()
    figure_five()
    figure_six()
    denv2_tree("fig3_cluster.svg", root_name="NODE_0003277", band_name="NODE_0003278",
               row=31, scale=1.18)
    denv2_tree("fig8_outbreak_clade.svg", root_marker=("E", "S7A"), levels_up=2,
               collapse_marker=("NS5", "T363I"),
               collapse_label="Hillsborough and Pinellas cluster")
    denv2_tree("fig9_traveler_clade.svg", root_marker=("NS2A", "I33L"), levels_up=2)
    serotype_tree("denv4", "fig11_denv4_tree.svg")
    serotype_tree("denv3", "fig12_denv3_tree.svg")
    figure_supplementary()


if __name__ == "__main__":
    build_dir_argument(__doc__)
    main()
