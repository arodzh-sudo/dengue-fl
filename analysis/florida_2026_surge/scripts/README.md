# Report scripts

These build the Florida dengue report from one finished Nextstrain build. They read only the
build directory you point them at and write back into it, so the same command reproduces the
report anywhere. Nothing here needs a package beyond the Python standard library.

## What the build directory has to hold

```
<build-dir>/
  json/                     dengue_denv{1,2,3,4}_genome.json from phylogenetic/auspice/
  summary_report.txt        the Daytona_dengue run summary
  input_report.txt          local/results/input_report.txt
  validation_report.txt     local/results/validation_report.txt
  metadata.txt              the sample table, the same one local/ was given
```

`results/` and `report/` are created inside it.

## Running it

```sh
bash scripts/build_report.sh /path/to/build-dir
```

That runs the four steps in order. Each also runs on its own:

```sh
python3 scripts/analyze_v2.py   --build-dir /path/to/build-dir   # counts, clades, distances
python3 scripts/make_figures.py --build-dir /path/to/build-dir   # twelve SVG figures
bash    scripts/rasterize.sh    /path/to/build-dir               # the same figures as PNG
python3 scripts/make_report.py  --build-dir /path/to/build-dir   # HTML and Word
```

On HiPerGator `rsvg-convert` is on the path and is used automatically. The script also accepts
inkscape, cairosvg or a chromium build, and falls back to headless Edge on Windows. Without any
of them the HTML report still renders, since it references the SVG files, and the Word file is
skipped with a message.

## What each script does

| Script | Output |
|---|---|
| `analyze_v2.py` | `results/`: the specimen funnel, every Florida tip, the outbreak cluster and its structure, pairwise SNPs, the two DENV2 clades and what they are related to, close pairs, the lineage in context, per-week and per-county counts, and `findings.txt` as a log of the run |
| `make_figures.py` | `report/figures/*.svg`, twelve numbered figures plus one supplementary tree |
| `rasterize.sh` | `report/figures/*.png` at twice the nominal size |
| `make_report.py` | `report/florida_dengue_2026_v3.html` and `.docx`, both written from one copy of the text |

Auspice map screenshots are dropped into `report/figures/` by hand as `fig7a_map_exposure.png`.
Every other figure is drawn from the build.

## Editing the report

All of the prose lives in `content()` in `make_report.py`. The house rules for it: no semicolons,
no comma before "and", SNPs rather than nucleotide differences, US spelling, and figures carry no
title inside the image because the caption does that job.
