# Report scripts

These build the Florida dengue report from a finished Nextstrain build. They know nothing about
any particular build: the clusters, the clades, the counties and the colors are all discovered
from the files you point them at. Nothing here needs a package beyond the Python standard library.

## Layout

```
analysis/florida_2026_surge/
  scripts/                    these files, the only part that is tracked in git
  builds/
    2026-09-15/               one folder per build, named for the day it was built
      json/                   dengue_denv{1,2,3,4}_genome.json from phylogenetic/auspice/
      summary_report.txt      the Daytona_dengue run summary
      input_report.txt        local/results/input_report.txt
      validation_report.txt   local/results/validation_report.txt
      metadata.txt            the sample table, the same one local/ was given
      report.yaml             what this report features
      content.py              the narrative for this build
      results/                written by analyze_v2.py
      report/                 figures, HTML and Word
```

A bare command uses the newest dated folder. `--build-dir` points at any other.

## Building a report for a new build

1. Run the Nextstrain workflows as usual. They are entirely independent of anything here.
2. Make a dated folder and copy the five inputs into it.
3. Run the analysis and read what it found:

   ```sh
   python3 scripts/analyze_v2.py --build-dir builds/2026-11-20
   ```

   Look at `results/clusters.tsv` for every cluster of genomes acquired in Florida,
   `results/unclustered.tsv` for the genomes no cluster claimed, which are findings in their own
   right, and `results/clades.tsv` for the parts of the tree the Florida genomes sit in.
4. Write `report.yaml`, copying the previous build's as a starting point. It says which cluster
   the report features and which extra cases belong in the distance matrix.
5. Copy the previous `content.py` next to the new build and edit the narrative.
6. See what moved:

   ```sh
   python3 scripts/compare_builds.py builds/2026-09-15 builds/2026-11-20
   ```

   Every fact it lists is a sentence to reread before the report goes out.
7. Build it:

   ```sh
   bash scripts/build_report.sh builds/2026-11-20
   ```

## What each script does

| Script | Output |
|---|---|
| `analyze_v2.py` | `results/`: the specimen funnel, every Florida tip, the clusters and clades it found, the genomes in no cluster with their nearest relatives, pairwise SNPs, the lineage in context, per week and per county counts, `landmarks.tsv` naming the tree nodes the figures need, `report_facts.tsv` for checking a draft, and `findings.txt` as a log of the run |
| `make_figures.py` | `report/figures/*.svg`, twelve numbered figures plus one supplementary tree |
| `rasterize.sh` | `report/figures/*.png` at twice the nominal size |
| `make_report.py` | `report/florida_dengue_2026_v3.html` and `.docx`, both from the build's `content.py` |
| `compare_builds.py` | the facts that moved between two builds |
| `build_report.sh` | the four steps in order |

On HiPerGator `rsvg-convert` is on the path and is used automatically. The script also accepts
inkscape, cairosvg or a chromium build, and falls back to headless Edge on Windows. Without any
of them the HTML still renders, since it references the SVG files, and the Word file is skipped
with a message.

## How things are found rather than named

- **A cluster** is a clade where every genome was acquired in Florida, mosquito pools included,
  with at least `min_local` of them. The report features the largest, or the one holding the
  sample named in `report.yaml`.
- **A clade** is a part of the tree holding at least `min_florida` Florida genomes with at most
  `max_public` public genomes inside. Each is described by what surrounds it, which is written to
  `results/clade_context.tsv`.
- **The nearest relative** of a genome is found among the tips in the smallest clade above it
  that holds enough of them to compare, then measured exactly in SNPs.
- **Counties** take the fixed hues in order of how many genomes each has, so a county that has
  never appeared before still gets a color and a legend entry.

Nothing in `scripts/` names a node, a sample or a county. If a build no longer contains something
`report.yaml` names, the run stops and says which one.

## Auspice screenshots

Figure 7 is a screenshot, not something the scripts draw. Drop it into `report/figures/` as
`fig7a_map_exposure.png` before the final run. Without it the HTML builds with a gap and the Word
file is skipped.

## Editing the narrative

The prose lives in `content()` in the build's `content.py`. The house rules: no semicolons, no
comma before "and", SNPs rather than nucleotide differences, US spelling, and figures carry no
title inside the image because the caption does that job.
