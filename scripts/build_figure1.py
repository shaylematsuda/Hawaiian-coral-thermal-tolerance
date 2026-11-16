"""Rebuild Figure 1 panels from source assets; never substitute schematic data for S1.
Official R notebook: ACE21_Figure1.Rmd.
CLI: python scripts/build_figure1.py (openpyxl, reportlab, pypdf; Rscript and pdftoppm).
Use --assemble-only to reuse the freshly exported sensor panel without rerunning it.
The approved sensor workflow exports Figures/assets/Fig1/Fig1B_temperature.pdf.
"""
import csv
import json
import math
from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.lib.utils import ImageReader
from reportlab.platypus import Table, TableStyle, Paragraph
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader, PdfWriter, Transformation

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / 'Figures/assets/Fig1'
FIGURES = ROOT / 'Figures'
SITE = (21.463086, -157.812518)


def build_map():
    meta = json.loads((ASSETS / 'kaneohe_satellite_metadata.json').read_text())
    extent = meta['extent']
    width, height = 432, 360
    path = FIGURES / 'Fig1A_collection_map.pdf'
    c = canvas.Canvas(str(path), pagesize=(width, height))
    c.drawImage(str(ASSETS / 'kaneohe_satellite.jpg'), 0, 0, width, height)
    radius = 6378137
    lat, lon = SITE
    east = radius * math.radians(lon)
    north = radius * math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
    x = (east - extent['xmin']) / (extent['xmax'] - extent['xmin']) * width
    y = (north - extent['ymin']) / (extent['ymax'] - extent['ymin']) * height
    c.setStrokeColor(colors.HexColor('#FFD700'))
    c.setLineWidth(1.5)
    c.circle(x, y, 4.5, stroke=1, fill=0)
    c.line(x + 5, y + 4, x + 44, y + 34)
    c.setFillColor(colors.HexColor('#FFD700'))
    c.setFont('Helvetica-Bold', 10)
    c.drawString(x + 47, y + 34, 'Collection reef')
    c.setFillColor(colors.white)
    c.setFont('Helvetica-Bold', 12)
    c.drawString(10, height - 20, "Kane'ohe Bay, O'ahu, Hawai'i")
    c.setFont('Helvetica', 6)
    c.drawString(8, 9, 'Imagery: Esri, Maxar, Earthstar Geographics and GIS User Community')
    c.setFont('Helvetica', 7)
    c.drawString(10, height - 33, '21.463086 N, 157.812518 W')
    c.save()
    return path


def build_traits():
    # Embed all font faces so PDF viewers preserve scientific-name italics.
    font_dir = Path('/System/Library/Fonts/Supplemental')
    faces = {'Fig1Arial': 'Arial.ttf', 'Fig1ArialBold': 'Arial Bold.ttf',
             'Fig1ArialItalic': 'Arial Italic.ttf',
             'Fig1ArialBoldItalic': 'Arial Bold Italic.ttf'}
    for name, filename in faces.items():
        pdfmetrics.registerFont(TTFont(name, str(font_dir / filename)))
    pdfmetrics.registerFontFamily('Fig1Arial', normal='Fig1Arial',
        bold='Fig1ArialBold', italic='Fig1ArialItalic', boldItalic='Fig1ArialBoldItalic')
    # Exact transcription of the supplied figure; scientific values are not re-estimated.
    # TODO (author): replace * and ** with actual reference numbers when supplied.
    # Preserve the placeholders until then; do not infer or fabricate citations.
    # Reproductive mode is retained pending user clarification; it repeats the
    # transmission labels in the supplied PDF and is documented in the source notes.
    rows = [
        ['Tissue thickness (µm)**', '1678', '1684', '1313', '491'],
        ['Skeletal morphology', 'Porous', 'Porous', 'Non-porous, ridges<br/>and valleys', 'Non-porous'],
        ['Growth form', 'Branching, plating', 'Branching', 'Encrusting', 'Branching'],
        ['Symbiodiniaceae<br/>transmission mode', 'Vertical', 'Vertical', 'Horizontal', 'Vertical'],
        ["Bleaching susceptibility<br/>(2014/2015 events<br/>in Kane'ohe Bay)*", '36.7/16.6%', '43.7/19.7%', '62.9/24.9%', '90.6/81.4%*'],
        ['Reproductive mode', 'Vertical', 'Vertical', 'Horizontal', 'Vertical'],
        ['Symbiodiniaceae genera', '<i>Cladocopium,<br/>Durusdinium</i>', '<i>Cladocopium</i>', '<i>Cladocopium</i>', '<i>Cladocopium</i>'],
    ]
    with (ASSETS / 'coral_traits_source.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Trait', 'Montipora capitata', 'Porites compressa', 'Pavona varians', 'Pocillopora acuta'])
        writer.writerows(rows)
    path = FIGURES / 'Fig1C_traits.pdf'
    width, height = 540, 300
    c = canvas.Canvas(str(path), pagesize=(width, height))
    style = ParagraphStyle('cell', fontName='Fig1Arial', fontSize=8.2, leading=9.6, alignment=1)
    left = ParagraphStyle('left', parent=style, alignment=0)
    headers = ['', '<b><i>M. capitata</i></b>', '<b><i>P. compressa</i></b>', '<b><i>P. varians</i></b>', '<b><i>P. acuta</i></b>']
    data = [[Paragraph(v, style) for v in headers]]
    data += [[Paragraph(v, left if i == 0 else style) for i, v in enumerate(row)] for row in rows]
    table = Table(data, colWidths=[132, 102, 102, 102, 102], rowHeights=[22, 22, 34, 26, 30, 40, 22, 36])
    table.setStyle(TableStyle([
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('GRID', (0, 1), (-1, -1), .35, colors.grey),
        ('LEFTPADDING', (0, 0), (-1, -1), 5),
        ('RIGHTPADDING', (0, 0), (-1, -1), 5),
    ]))
    table.wrapOn(c, width, height)
    table.drawOn(c, 0, 66)
    for i in range(4):
        image = ASSETS / f'reference-{i+1:03d}.jpg'
        c.drawImage(str(image), 132+i*102, 21, 102, 44, preserveAspectRatio=False)
    c.save()
    return path


def assemble(map_path, traits_path, temperature_path):
    # Same A/B top-row and full-width C arrangement as the supplied figure.
    writer = PdfWriter()
    page = writer.add_blank_page(width=612, height=792)
    for path, x, y, w, h in [
        (map_path, 25, 490, 255, 212.5),
        (temperature_path, 303, 490, 284, 212.5),
        (traits_path, 25, 125, 562, 312.2),
    ]:
        source = PdfReader(str(path)).pages[0]
        scale = min(w / float(source.mediabox.width), h / float(source.mediabox.height))
        page.merge_transformed_page(source, Transformation().scale(scale).translate(x, y))
    # Vector labels and explanatory notes remain editable in the assembled PDF.
    from io import BytesIO
    overlay = BytesIO()
    labels = canvas.Canvas(overlay, pagesize=(612, 792))
    labels.setFont('Helvetica-Bold', 14)
    for label, x, y in [('A', 25, 719), ('B', 303, 719), ('C', 25, 451)]:
        labels.drawString(x, y, label)
    note_style = ParagraphStyle('notes', fontName='Helvetica', fontSize=7.3, leading=10)
    notes = Paragraph(
        '<b>B.</b> Mean temperature ± 1 SE across available tanks at each timestamp, '
        'March 20-May 17, 2017. Logger-removal data for tanks 1-4 during May 1-15 '
        'and six cleaning-interruption readings on April 17-19 were excluded. '
        'SE is not shown where fewer than two tanks remain.',
        note_style)
    _, note_height = notes.wrap(562, 80)
    notes.drawOn(labels, 25, 90-note_height)
    labels.save()
    overlay.seek(0)
    page.merge_page(PdfReader(overlay).pages[0])
    with (FIGURES / 'Fig1.pdf').open('wb') as f:
        writer.write(f)


if __name__ == '__main__':
    import argparse
    import shutil
    import subprocess
    import sys
    parser = argparse.ArgumentParser(description='Build the official Figure 1 from source panels.')
    parser.add_argument('--assemble-only', action='store_true',
        help='Import the existing Fig1B export; do not rerun temperature preparation.')
    args = parser.parse_args()
    FIGURES.mkdir(exist_ok=True)
    if not args.assemble_only:
        # The default CLI rebuilds all temperature inputs before importing panel B.
        subprocess.run([sys.executable, str(ROOT / 'scripts/import_apex_records.py')],
                       cwd=ROOT, check=True)
        subprocess.run(['Rscript', str(ROOT / 'scripts/plot_apex_temperature.R')],
                       cwd=ROOT, check=True)
    temperature = ASSETS / 'Fig1B_temperature.pdf'
    if not temperature.exists():
        raise FileNotFoundError('Run the Figure 1B chunk in ACE21_Sensor data.Rmd first: ' + str(temperature))
    a = build_map()
    c = build_traits()
    assemble(a, c, temperature)
    renderer = shutil.which('pdftoppm')
    if renderer:
        subprocess.run([renderer, '-scale-to', '1800', '-png', '-singlefile',
                        str(FIGURES / 'Fig1.pdf'), str(FIGURES / 'Fig1')], check=True)
    else:
        raise RuntimeError('Fig1.pdf saved; install Poppler (pdftoppm) to render the required Fig1.png preview.')
    print('Saved Figures/Fig1.pdf and Figures/Fig1.png from the exported Figure 1B panel.')
