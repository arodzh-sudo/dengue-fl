#!/usr/bin/env bash
# Build the Florida dengue report from one Nextstrain build directory.
# Usage: build_report.sh [build-dir]      default: ../2026
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build="$(cd "${1:-$here/../2026}" && pwd)"

python3 "$here/analyze_v2.py" --build-dir "$build"
python3 "$here/make_figures.py" --build-dir "$build"
bash "$here/rasterize.sh" "$build" || echo "skipping PNG conversion, the Word file will be short of figures"
python3 "$here/make_report.py" --build-dir "$build"

echo
echo "results: $build/results"
echo "report:  $build/report"
