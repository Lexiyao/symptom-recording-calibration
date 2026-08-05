"""
Extract two panels from Barclay et al. (2024) Figure 2 and add a 3% reference
line at a position calibrated from the figure's own gridlines.

Source:  Barclay ME, Renzi C, Harrison H, et al. Cancer incidence and competing
         mortality risk following 15 presenting symptoms in primary care.
         BMJ Oncology 2024;3:e000500.  doi:10.1136/bmjonc-2024-000500
Licence: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

Modifications made to the original figure, declared here and on the slide:
  1. Two of the fifteen panels (Haematuria, Fatigue) are shown; the rest are cropped away.
  2. A horizontal 3% reference line is added.
  3. Two short text annotations are added, reporting values measured off the
     figure by measure_bands() below.

Nothing in the plotted data is altered. The 3% line is not placed by eye: the
panel's y-axis is calibrated by locating its own gridlines and the black
baseline, fitting rows against known axis values, and asserting the fit
residual is under 3 pixels before anything is drawn.

Usage:  python code/barclay_panels.py path/to/bmjonc-2024-000500.pdf

The PDF is not redistributed in this repository. Download it from the
publisher (open access) and pass the path.
"""

import collections
import pathlib
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "figures"

FIG2_PAGE = 8          # 1-indexed page of Figure 2 in the published PDF
RENDER_SCALE = 6.0     # pdfium render scale; all pixel constants assume this
NAVY = (31, 78, 121)

# stacked-area fill colours, sampled from the paper's own legend swatches
SITE_COLOURS = {
    "Haematological": (76, 147, 195),
    "Lower GI": (252, 174, 173),
    "Prostate": (193, 229, 161),
    "Lung": (253, 204, 140),
    "Upper GI": (233, 72, 73),
    "Urological": (92, 179, 86),
    "Other cancer": (200, 224, 240),
}


def render_page(pdf_path):
    import pypdfium2 as pdfium
    doc = pdfium.PdfDocument(str(pdf_path))
    return doc[FIG2_PAGE - 1].render(scale=RENDER_SCALE).to_pil()


def find_title_bar_rows(page):
    """Panel rows are separated by wide solid-black title bars."""
    dark = np.array(page.convert("L")) < 70
    frac = dark.mean(axis=1)
    bands, inside, start = [], False, 0
    for y, v in enumerate(frac):
        if v > 0.25 and not inside:
            inside, start = True, y
        elif v <= 0.25 and inside:
            inside = False
            if y - start > 10:
                bands.append((start, y))
    return bands


def find_panel_columns(page, bar_top, bar_bottom):
    """x-extent of each title bar in one panel row."""
    dark = np.array(page.convert("L")) < 70
    frac = dark[bar_top:bar_bottom, :].mean(axis=0)
    segs, inside, start = [], False, 0
    for x, v in enumerate(frac):
        if v > 0.5 and not inside:
            inside, start = True, x
        elif v <= 0.5 and inside:
            inside = False
            if x - start > 80:
                segs.append((start, x))
    return segs


def calibrate_axis(page, top, x0, x1):
    """Return (row_of_3pc, baseline_row) for one panel, from its own gridlines.

    Raises if the linear fit of gridline row against axis value does not hold
    to within 3 pixels — a silently miscalibrated reference line would be a
    data-fidelity error, not a cosmetic one.
    """
    grey = np.array(page.convert("L"))
    # the 0.00 axis is a near-black line spanning the full plot width
    dark_frac = (grey[top:top + 1000, x0 + 20:x1 - 12] < 90).mean(axis=1)
    baseline = next(y for y, v in enumerate(dark_frac) if v > 0.85 and y > 600)

    soft = (grey[top:top + int(baseline), x0 + 20:x1 - 12] < 232).mean(axis=1)
    rows, group = [], []
    for y, v in enumerate(soft):
        if v > 0.60:
            group.append(y)
        elif group:
            rows.append(int(np.mean(group)))
            group = []
    if group:
        rows.append(int(np.mean(group)))

    # gridlines sit at 0.05 intervals above the baseline; keep the ones that
    # land within a pixel of a multiple and label them from the baseline up
    step = None
    for r in rows:
        d = baseline - r
        if d > 100:
            step = d / round(d / 133.0) if round(d / 133.0) else None
            if step and 125 < step < 140:
                break
    if step is None:
        raise AssertionError("could not establish gridline spacing")

    known = {(baseline - r) / step * 0.05: r for r in rows if baseline - r > 50}
    pts = {round(v * 20) / 20: r for v, r in known.items() if abs(v * 20 - round(v * 20)) < 0.25}
    pts[0.00] = baseline
    vals = np.array(sorted(pts))
    obs = np.array([pts[v] for v in vals], dtype=float)
    slope, intercept = np.polyfit(vals, obs, 1)
    resid = np.abs(obs - (slope * vals + intercept)).max()
    if resid > 3.0:
        raise AssertionError(f"axis calibration residual {resid:.1f}px exceeds 3px")
    return slope * 0.03 + intercept, baseline


def find_plot_right_edge(page, top, baseline, x_from, x_to):
    """Rightmost x carrying stacked-area fill.

    The black title bar is WIDER than the plotting area, so its right edge must
    not be used as a sampling column — doing so samples empty margin and
    silently returns near-zero band thicknesses.
    """
    hsv = np.array(page.convert("HSV"))
    saturated = (hsv[:, :, 1].astype(int) > 60) & (hsv[:, :, 2].astype(int) > 90)
    strip = saturated[top:top + int(baseline), x_from:x_to]
    filled = np.nonzero(strip.mean(axis=0) > 0.01)[0]
    if not len(filled):
        raise AssertionError("no stacked-area fill found in panel")
    return x_from + int(filled.max())


def measure_bands(page, top, baseline, plot_right, rows_per_005):
    """Thickness of each cancer site's band at the oldest age, in risk units.

    Barclay's claim concerns risk at INDIVIDUAL sites, but the panels are
    stacked, so the stack total is not the quantity of interest. This measures
    each band separately so the annotation states a single-site value.
    """
    rgb = np.array(page.convert("RGB"))
    column = rgb[top:top + int(baseline), plot_right - 6:plot_right - 2].reshape(-1, 3)
    counts = collections.Counter()
    for pixel in column:
        dist, site = min(
            (sum((int(a) - int(b)) ** 2 for a, b in zip(colour, pixel)), name)
            for name, colour in SITE_COLOURS.items())
        if dist < 2500:
            counts[site] += 1
    per_px = 0.05 / rows_per_005
    return {site: (n / 4) * per_px for site, n in counts.items()}


def compose(page, pdf_path, scale=0.70):
    bars = find_title_bar_rows(page)
    row3_bar = bars[2]                       # third panel row holds Haematuria and Fatigue
    top = row3_bar[0] - 8
    cols = find_panel_columns(page, *row3_bar)
    haem, fatigue = cols[-2], cols[-1]

    row_3pc, baseline = calibrate_axis(page, top, *haem)
    rows_per_005 = (baseline - row_3pc) / 0.03 * 0.05

    right_h = find_plot_right_edge(page, top, baseline, haem[0] - 360, haem[1] + 40)
    right_f = find_plot_right_edge(page, top, baseline, fatigue[0] - 360, fatigue[1] + 40)
    measured = {"Haematuria": measure_bands(page, top, baseline, right_h, rows_per_005),
                "Fatigue": measure_bands(page, top, baseline, right_f, rows_per_005)}
    for panel, bands in measured.items():
        top_site, top_val = max(bands.items(), key=lambda kv: kv[1])
        print(f"{panel}: largest single site = {top_site} at {top_val * 100:.1f}% "
              f"({'exceeds' if top_val > 0.03 else 'below'} the 3% threshold)")
    # the whole point of the annotation: haematuria has a single site above 3%,
    # fatigue does not. Assert it rather than trusting the eye.
    assert max(measured["Haematuria"].values()) > 0.03
    assert max(measured["Fatigue"].values()) < 0.03

    def left_gutter(x_start):
        """Walk left to the white gutter so no neighbouring panel bleeds in."""
        ink = (np.array(page.convert("L"))[top:baseline + top, :] < 245).mean(axis=0)
        run, x = 0, x_start - 1
        while x > 0:
            if ink[x] < 0.005:
                run += 1
                if run >= 22:
                    return x + run
            else:
                run = 0
            x -= 1
        return x_start - 320

    e_h, e_f = left_gutter(haem[0]), left_gutter(fatigue[0])
    axis_y = top + int(baseline)
    strips = [(page.crop((e_h, top, haem[1] + 16, axis_y + 2)), e_h, haem[0],
               f"Urological alone\nreaches {measured['Haematuria']['Urological'] * 100:.1f}%", NAVY),
              (page.crop((e_f, top, fatigue[1] + 16, axis_y + 2)), e_f, fatigue[0],
               "No single site\nreaches 3%", (110, 110, 110))]
    labels = page.crop((e_f, axis_y + 2, fatigue[1] + 16, axis_y + 88))

    pw, ph = strips[0][0].size
    tw, th = int(pw * scale), int(ph * scale)
    lh = int(labels.size[1] * scale)
    gap, pad_l, pad_t = 34, 10, 10
    try:
        f_thr = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 23)
        f_ann = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 21)
        f_note = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 20)
        f_src = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 17)
    except OSError:
        f_thr = f_ann = f_note = f_src = ImageFont.load_default()

    note = ("Cumulative risk to 12 months by age (years), stacked by cancer site "
            "\u00b7 male non-smokers")
    src = ("Barclay et al. (2024) BMJ Oncology 3:e000500, Fig 2 (CC BY 4.0) "
           "\u2014 2 of 15 panels; annotation added")
    width = max(pad_l + tw * 2 + gap + 10, pad_l + int(f_src.getlength(src)) + 14)
    canvas = Image.new("RGB", (width, pad_t + th + lh + 54), "white")
    draw = ImageDraw.Draw(canvas)
    y3 = pad_t + int(row_3pc * scale)

    for i, (strip, edge, plot_x0, annotation, colour) in enumerate(strips):
        x = pad_l + i * (tw + gap)
        w = int(strip.size[0] * scale)
        canvas.paste(strip.resize((w, th), Image.LANCZOS), (x, pad_t))
        canvas.paste(labels.resize((w, lh), Image.LANCZOS), (x, pad_t + th))
        plot_left = x + int((plot_x0 - edge) * scale)
        draw.line([(plot_left, y3), (x + w - 3, y3)], fill=NAVY, width=3)
        draw.multiline_text((plot_left + 10, pad_t + int(110 * scale)),
                            annotation, font=f_ann, fill=colour, spacing=5)

    draw.text((pad_l + int((strips[0][2] - strips[0][1]) * scale) + 10, y3 + 7),
              "3% referral threshold", font=f_thr, fill=NAVY)
    y_note = pad_t + th + lh + 4
    draw.text((pad_l + 2, y_note), note, font=f_note, fill=(70, 70, 70))
    draw.text((pad_l + 2, y_note + 25), src, font=f_src, fill=(125, 125, 125))

    OUT.mkdir(exist_ok=True)
    out = OUT / "fig_barclay_panels.png"
    canvas.save(out)
    print(f"wrote {out}")
    return measured


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__.strip().splitlines()[-3])
    compose(render_page(pathlib.Path(sys.argv[1])), sys.argv[1])
