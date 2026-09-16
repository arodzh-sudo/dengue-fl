#!/usr/bin/env bash
# Rasterize the report SVGs to PNG, which the Word file embeds.
# Usage: rasterize.sh [build-dir]        default: ../2026
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build="${1:-$here/../2026}"
dir="$(cd "$build/report/figures" && pwd)"

edge=""
for candidate in "/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"                  "/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"; do
  [ -x "$candidate" ] && edge="$candidate" && break
done
if command -v rsvg-convert >/dev/null 2>&1; then
  engine=rsvg
elif command -v inkscape >/dev/null 2>&1; then
  engine=inkscape
elif command -v chromium >/dev/null 2>&1; then
  engine=chromium; browser=chromium
elif command -v chromium-browser >/dev/null 2>&1; then
  engine=chromium; browser=chromium-browser
elif command -v google-chrome >/dev/null 2>&1; then
  engine=chromium; browser=google-chrome
elif python3 -c "import cairosvg" >/dev/null 2>&1; then
  engine=cairosvg
elif [ -n "$edge" ]; then
  engine=chromium; browser="$edge"
else
  echo "no SVG converter found: install rsvg-convert, inkscape, cairosvg or a chromium build" >&2
  echo "the HTML report renders from the SVGs and needs no conversion" >&2
  exit 1
fi

for svg in "$dir"/*.svg; do
  name="$(basename "$svg" .svg)"
  png="$dir/$name.png"
  read -r width height < <(sed -n 's/.*width="\([0-9]*\)" height="\([0-9]*\)".*/\1 \2/p' "$svg" | head -1)
  case "$engine" in
    rsvg)     rsvg-convert -z 2 -o "$png" "$svg" ;;
    inkscape) inkscape --export-type=png --export-dpi=192 --export-filename="$png" "$svg" >/dev/null 2>&1 ;;
    cairosvg) python3 -c "import cairosvg,sys; cairosvg.svg2png(url=sys.argv[1], write_to=sys.argv[2], scale=2)" "$svg" "$png" ;;
    chromium)
      target="$png"; source="$svg"
      if command -v cygpath >/dev/null 2>&1; then
        target="$(cygpath -w "$png")"; source="file:///$(cygpath -w "$svg" | tr '\\' '/')"
      else
        source="file://$svg"
      fi
      MSYS_NO_PATHCONV=1 "$browser" --headless --disable-gpu --force-device-scale-factor=2 \
        --screenshot="$target" --window-size="$width,$height" --hide-scrollbars \
        "$source" >/dev/null 2>&1 ;;
  esac
  echo "$name.png"
done
