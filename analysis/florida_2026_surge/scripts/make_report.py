#!/usr/bin/env python3
"""Write the Florida dengue report once, as HTML and as Word."""

import os
import struct
import zipfile
from xml.sax.saxutils import escape

import importlib.util

import analyze_v2
from analyze_v2 import build_dir_argument


def narrative():
    """The blocks for this build, from content.py beside it."""
    path = os.path.join(analyze_v2.BASE, "content.py")
    if not os.path.exists(path):
        raise SystemExit("no content.py in {}, copy one from an earlier build".format(
            analyze_v2.BASE))
    spec = importlib.util.spec_from_file_location("content", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def report_dir():
    return os.path.join(analyze_v2.BASE, "report")


def figures_dir():
    return os.path.join(report_dir(), "figures")

# ---------------------------------------------------------------- HTML
STYLE = """
:root { --surface: #ffffff; --ink: #0b0b0b; --soft: #52514e; --muted: #8a8a85;
        --rule: #e4e3df; --accent: #2a78d6; }
* { box-sizing: border-box; }
body { margin: 0; background: #f2f1ee; color: var(--ink);
       font-family: Inter, 'Helvetica Neue', Arial, sans-serif; line-height: 1.55; }
main { max-width: 980px; margin: 0 auto; padding: 48px 32px 96px; background: var(--surface);
       min-height: 100vh; }
h1 { font-size: 30px; line-height: 1.2; margin: 0 0 8px; }
.dek { font-size: 16px; color: var(--soft); margin: 0 0 24px; max-width: 70ch; }
h2 { font-size: 19px; margin: 44px 0 12px; padding-bottom: 6px; border-bottom: 2px solid var(--rule); }
p { max-width: 78ch; margin: 0 0 14px; }
ul { max-width: 78ch; padding-left: 20px; }
li { margin-bottom: 9px; }
table { border-collapse: collapse; width: 100%; margin: 8px 0 18px; font-size: 14px; }
th { text-align: left; font-weight: 600; color: var(--soft); border-bottom: 2px solid var(--rule);
     padding: 8px 10px; vertical-align: top; }
td { border-bottom: 1px solid var(--rule); padding: 8px 10px; vertical-align: top; }
td:first-child, th:first-child { padding-left: 0; }
figure { margin: 26px 0; }
figure img { width: 100%; height: auto; border: 1px solid var(--rule); border-radius: 6px;
             background: var(--surface); }
figcaption { font-size: 12.5px; color: var(--muted); margin-top: 8px; max-width: 86ch; }
figcaption b { color: var(--soft); }
.meta { display: flex; flex-wrap: wrap; gap: 26px; border-top: 1px solid var(--rule);
        border-bottom: 1px solid var(--rule); padding: 14px 0; margin-bottom: 8px; }
.meta div span { display: block; font-size: 11px; text-transform: uppercase;
                 letter-spacing: .04em; color: var(--muted); }
.meta div b { font-size: 15px; font-weight: 600; }
code { font-family: 'SF Mono', Menlo, Consolas, monospace; font-size: 13px; }
.footer { margin-top: 48px; padding-top: 14px; border-top: 1px solid var(--rule);
          font-size: 12.5px; color: var(--muted); max-width: 86ch; }
@media (max-width: 720px) { main { padding: 28px 16px 64px; } h1 { font-size: 24px; } }
"""


def render_html(blocks, path):
    out = ['<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">',
           '<meta name="viewport" content="width=device-width, initial-scale=1">',
           "<title>{}</title><style>{}</style></head><body><main>".format(escape(TITLE), STYLE),
           "<h1>{}</h1>".format(escape(TITLE)),
           '<p class="dek">{}</p>'.format(escape(SUBTITLE))]
    for item in blocks:
        kind = item["kind"]
        if kind == "meta":
            out.append('<div class="meta">' + "".join(
                "<div><span>{}</span><b>{}</b></div>".format(escape(k), escape(v))
                for k, v in item["rows"]) + "</div>")
        elif kind == "h2":
            out.append("<h2>{}</h2>".format(escape(item["text"])))
        elif kind == "p":
            out.append("<p>{}</p>".format(mark(item["text"])))
        elif kind == "bullets":
            out.append("<ul>" + "".join("<li>{}</li>".format(mark(i)) for i in item["items"]) + "</ul>")
        elif kind == "table":
            head = "".join("<th>{}</th>".format(escape(c)) for c in item["header"])
            body = "".join("<tr>{}</tr>".format(
                "".join("<td>{}</td>".format(mark(c)) for c in row)) for row in item["rows"])
            out.append("<table><thead><tr>{}</tr></thead><tbody>{}</tbody></table>".format(head, body))
        elif kind == "figure":
            source = "figures/{}.{}".format(item["name"], "png" if item.get("screenshot") else "svg")
            out.append('<figure><img src="{}" alt="{}"><figcaption><b>Figure {}.</b> {}'
                       "</figcaption></figure>".format(source, escape(item["caption"][:120]),
                                                       item["number"], mark(item["caption"])))
        elif kind == "footer":
            out.append('<div class="footer">{}</div>'.format(mark(item["text"])))
    out.append("</main></body></html>")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(out))
    print("wrote", os.path.basename(path))


def mark(text):
    """Escape, then set identifiers and file names in the monospace face."""
    escaped = escape(text)
    words = []
    for word in escaped.split(" "):
        bare = word.strip(".,;:()")
        if bare and (bare.startswith(("TVU", "TVA", "JVV", "MosquitoPool", "NODE_", "OQ", "PQ", "PP",
                                      "PV", "PX", "OR"))
                     or bare.endswith((".json", ".tsv", ".txt", ".py"))
                     or bare in ("case_origin", "metadata.txt", "augur", "2II_F.1.1.2",
                                 "3III_B.3.2", "4II_B.1.3")):
            word = word.replace(bare, "<code>{}</code>".format(bare))
        words.append(word)
    return " ".join(words)


# ---------------------------------------------------------------- Word
EMU_PER_PX = 9525
PAGE_WIDTH_EMU = 5943600


def png_size(path):
    with open(path, "rb") as handle:
        header = handle.read(24)
    return struct.unpack(">II", header[16:24])


def w_paragraph(runs, style=None, spacing=(0, 120)):
    properties = "<w:pPr>"
    if style:
        properties += '<w:pStyle w:val="{}"/>'.format(style)
    properties += '<w:spacing w:before="{}" w:after="{}"/></w:pPr>'.format(*spacing)
    return "<w:p>" + properties + runs + "</w:p>"


def w_run(text, bold=False, size=20, color="0B0B0B", mono=False):
    fonts = '<w:rFonts w:ascii="{0}" w:hAnsi="{0}"/>'.format("Consolas" if mono else "Calibri")
    return ('<w:r><w:rPr>{}{}<w:sz w:val="{}"/><w:color w:val="{}"/></w:rPr>'
            '<w:t xml:space="preserve">{}</w:t></w:r>'.format(
                fonts, "<w:b/>" if bold else "", size, color, escape(text)))


def w_image(index, width_emu, height_emu):
    return ('<w:p><w:pPr><w:spacing w:before="160" w:after="80"/></w:pPr><w:r><w:drawing>'
            '<wp:inline distT="0" distB="0" distL="0" distR="0">'
            '<wp:extent cx="{w}" cy="{h}"/><wp:docPr id="{i}" name="Figure {i}"/>'
            '<a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
            '<a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
            '<pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
            '<pic:nvPicPr><pic:cNvPr id="{i}" name="Figure {i}"/><pic:cNvPicPr/></pic:nvPicPr>'
            '<pic:blipFill><a:blip r:embed="rId{i}"/><a:stretch><a:fillRect/></a:stretch>'
            "</pic:blipFill>"
            '<pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{w}" cy="{h}"/></a:xfrm>'
            '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic>'
            "</a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>".format(
                w=width_emu, h=height_emu, i=index))


def w_table(header, rows):
    widths = [int(9000 / len(header))] * len(header)

    def cell(text, bold=False, shade=None):
        properties = '<w:tcPr><w:tcW w:w="{}" w:type="dxa"/>'.format(widths[0])
        if shade:
            properties += '<w:shd w:val="clear" w:fill="{}"/>'.format(shade)
        properties += "</w:tcPr>"
        return "<w:tc>" + properties + w_paragraph(w_run(text, bold=bold, size=18),
                                                   spacing=(40, 40)) + "</w:tc>"

    out = ['<w:tbl><w:tblPr><w:tblW w:w="9000" w:type="dxa"/>'
           '<w:tblBorders><w:top w:val="single" w:sz="4" w:color="E4E3DF"/>'
           '<w:bottom w:val="single" w:sz="4" w:color="E4E3DF"/>'
           '<w:insideH w:val="single" w:sz="4" w:color="E4E3DF"/></w:tblBorders></w:tblPr>']
    out.append("<w:tr>" + "".join(cell(c, bold=True, shade="F4F3F0") for c in header) + "</w:tr>")
    for row in rows:
        out.append("<w:tr>" + "".join(cell(c) for c in row) + "</w:tr>")
    out.append("</w:tbl>" + w_paragraph("", spacing=(0, 120)))
    return "".join(out)


def render_docx(blocks, path):
    missing = [item["name"] for item in blocks if item["kind"] == "figure"
               and not os.path.exists(os.path.join(figures_dir(), item["name"] + ".png"))]
    if missing:
        print("no PNG for {}, skipping the Word file".format(", ".join(missing)))
        print("run rasterize.sh once an SVG converter is available")
        return
    body, images = [], []
    body.append(w_paragraph(w_run(TITLE, bold=True, size=36)))
    body.append(w_paragraph(w_run(SUBTITLE, size=22, color="52514E")))

    for item in blocks:
        kind = item["kind"]
        if kind == "meta":
            body.append(w_table(["", ""], [[k, v] for k, v in item["rows"]]))
        elif kind == "h2":
            body.append(w_paragraph(w_run(item["text"], bold=True, size=26), spacing=(280, 100)))
        elif kind == "p":
            body.append(w_paragraph(w_run(item["text"])))
        elif kind == "bullets":
            for entry in item["items"]:
                body.append(w_paragraph(w_run("•  " + entry), spacing=(0, 80)))
        elif kind == "table":
            body.append(w_table(item["header"], item["rows"]))
        elif kind == "figure":
            name = item["name"] + ".png"
            source = os.path.join(figures_dir(), name)
            width, height = png_size(source)
            scale = min(PAGE_WIDTH_EMU / (width * EMU_PER_PX), 1.0)
            images.append(source)
            body.append(w_image(len(images), int(width * EMU_PER_PX * scale),
                                int(height * EMU_PER_PX * scale)))
            body.append(w_paragraph(
                w_run("Figure {}. ".format(item["number"]), bold=True, size=17, color="52514E")
                + w_run(item["caption"], size=17, color="8A8A85"), spacing=(0, 200)))
        elif kind == "footer":
            body.append(w_paragraph(w_run(item["text"], size=17, color="8A8A85"),
                                    spacing=(240, 0)))

    document = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
        'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing">'
        "<w:body>" + "".join(body) +
        '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/>'
        '<w:pgMar w:top="1080" w:right="1080" w:bottom="1080" w:left="1080"/></w:sectPr>'
        "</w:body></w:document>")

    relationships = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                     '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']
    overrides = []
    for index, source in enumerate(images, start=1):
        target = "media/{}".format(os.path.basename(source))
        relationships.append('<Relationship Id="rId{}" Type="http://schemas.openxmlformats.org/'
                             'officeDocument/2006/relationships/image" Target="{}"/>'.format(
                                 index, target))
        overrides.append('<Default Extension="png" ContentType="image/png"/>')
    relationships.append("</Relationships>")

    types = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
             '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
             '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.'
             'relationships+xml"/><Default Extension="xml" ContentType="application/xml"/>'
             '<Default Extension="png" ContentType="image/png"/>'
             '<Override PartName="/word/document.xml" ContentType="application/vnd.'
             'openxmlformats-officedocument.wordprocessingml.document.main+xml"/></Types>')

    root_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                 '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                 '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/'
                 '2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>')

    stem = path
    for attempt in range(1, 6):
        try:
            open(path, "ab").close()
            break
        except PermissionError:
            path = stem.replace(".docx", "-rebuilt{}.docx".format(attempt))
    if path != stem:
        print("the Word file is open, writing", os.path.basename(path), "instead")
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", types)
        archive.writestr("_rels/.rels", root_rels)
        archive.writestr("word/document.xml", document)
        archive.writestr("word/_rels/document.xml.rels", "".join(relationships))
        for source in images:
            archive.write(source, "word/media/" + os.path.basename(source))
    print("wrote", os.path.basename(path))


def main():
    module = narrative()
    global TITLE, SUBTITLE
    TITLE, SUBTITLE = module.TITLE, module.SUBTITLE
    blocks = module.content()
    render_html(blocks, os.path.join(report_dir(), "florida_dengue_2026_v3.html"))
    render_docx(blocks, os.path.join(report_dir(), "florida_dengue_2026_v3.docx"))


if __name__ == "__main__":
    build_dir_argument(__doc__)
    main()
