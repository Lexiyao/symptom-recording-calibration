"""
Figures for the DPhil interview presentation:
"Differential symptom recording and the calibration of cancer risk
prediction models in UK primary care"

Reproduces three figures from the CSVs in ../data:

  fig_forest.png: d'Elia et al. (2025) Table 3: odds of a symptom being
                       coded, by ethnicity and deprivation (CPRD Aurum)
  fig_threshold.png: Price et al. (2016): positive predictive values before
                       and after recovering free-text records
  fig_pooling.png: illustration of subgroup calibration slopes pooling to
                       a value that describes neither subgroup

Usage:  python code/figures.py            (writes into figures/)

Every number plotted here is read from ../data/*.csv, which were transcribed
from the published tables named in each file. No values are hard-coded in the
plotting code.
"""

import csv
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

NAVY = "#1F4E79"
GREY = "#6E6E6E"
RULE = "#C8CDD4"
FONT = "Arial"

GROUP_ORDER = ["Black", "South Asian", "Other ethnicity", "Most deprived (IMD 1)"]


def house_style():
    """Single font, three-step size ladder, open frame."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [FONT, "DejaVu Sans"],
        "font.size": 11.5,
        "axes.titlesize": 13,
        "axes.labelsize": 12.5,
        "xtick.labelsize": 12,
        "ytick.labelsize": 11.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#8A8A8A",
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def read_rows(name):
    with open(DATA / name, newline="") as fh:
        return list(csv.DictReader(fh))


def forest():
    """d'Elia Table 3: odds that a symptom was coded, by patient group.

    The paper reports significant differences for twelve symptoms. This shows 15
    rows because three symptoms are significant in more than one patient group,
    so the unit here is the symptom-by-group association, not the symptom.
    """
    raw = read_rows("delia2025_table3.csv")
    grouped = {g: [r for r in raw if r["group"] == g] for g in GROUP_ORDER}
    n_total = len(raw)
    n_above = sum(1 for r in raw if float(r["odds_ratio"]) > 1)

    fig, ax = plt.subplots(figsize=(7.7, 4.55))
    y = 0.0
    ticks, labels, uppers = [], [], []
    headers = []
    for g in GROUP_ORDER:
        y -= 0.95
        headers.append((y, g))
        for r in grouped[g]:
            y -= 1.0
            lo, hi = float(r["ci_low_99_76"]), float(r["ci_high_99_76"])
            orr = float(r["odds_ratio"])
            # the one association that runs the other way is drawn in grey
            colour = GREY if hi < 1 else NAVY
            ax.plot([lo, hi], [y, y], color=colour, lw=2.7, solid_capstyle="round", zorder=2)
            ax.plot([orr], [y], "o", ms=7.0, color=colour, zorder=3)
            ticks.append(y)
            labels.append(r["symptom"])
            uppers.append(hi)
    ybot = y - 0.55

    for yy, g in headers:
        ax.text(-0.415, yy, g, transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=12.5, fontweight="bold",
                color="#2D2D2D", clip_on=False)

    ax.axvline(1, color="#8A8A8A", ls=(0, (4, 3)), lw=1.3, zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(0.42, 11.6)
    ax.set_ylim(ybot, 0.15)
    ax.set_xticks([0.5, 1, 2, 4, 8])
    ax.set_xticklabels(["0.5", "1", "2", "4", "8"])
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Odds ratio that the symptom was coded  (99.76% CI)", labelpad=3)

    # annotate the exception by looking its row up, never by a guessed coordinate
    hip = ticks[labels.index("Hip pain")]
    ax.annotate("the one that runs\nthe other way", xy=(0.92, hip), xytext=(1.45, hip),
                fontsize=11, color=GREY, va="center", ha="left", linespacing=1.35,
                arrowprops=dict(arrowstyle="-", color=GREY, lw=0.9, shrinkA=3, shrinkB=2))
    ax.text(2.05, ticks[1], f"{n_above} of {n_total} symptom-by-group\nassociations lie above 1",
            fontsize=11.5, color=NAVY, fontweight="bold", ha="left", va="center",
            linespacing=1.40)

    fig.subplots_adjust(left=0.295, right=0.985, top=0.99, bottom=0.115)
    check_no_overlap(fig, ax)
    fig.savefig(OUT / "fig_forest.png", facecolor="white")
    plt.close(fig)
    return n_above, n_total



def check_line_clear_of_text(fig, ax, text_artist, x_data):
    """Assert no drawn line passes through a text label's box.

    check_no_overlap compares text against text. A dashed reference line is a
    Line2D, not a Text, so it passes that check while striking through the
    label. This closes that gap for one label and one vertical line.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    box = text_artist.get_window_extent(renderer)
    x_px = ax.transData.transform((x_data, 0))[0]
    if not (box.x0 <= x_px <= box.x1):
        return
    for line in ax.lines:
        xs, ys = line.get_xdata(), line.get_ydata()
        if len(xs) < 2 or not all(abs(v - x_data) < 1e-9 for v in xs):
            continue
        y_top_px = ax.transData.transform((x_data, max(ys)))[1]
        assert y_top_px <= box.y0 + 1, (
            f"the vertical line at x={x_data} reaches y={max(ys)}, which is "
            f"inside the label {text_artist.get_text()!r}"
        )

def threshold():
    """Price 2016: PPV with coded records only, then with recovered free text.

    All four symptom/cancer pairs the paper reports are shown, grouped as the
    paper groups them, into alarm features and one "low-risk but not no-risk"
    feature. Showing three of four would look like selection.

    The x axis is logarithmic because the four values span 0.14% to 12.8%. On a
    linear 0-15% axis the jaundice movement fills 43% of the width while
    haematuria, the only pair that crosses the referral threshold, fills 7%,
    so the visual emphasis lands on the secondary finding.
    """
    raw = read_rows("price2016_ppv.csv")
    groups = []
    for r in raw:
        if not groups or groups[-1][0] != r["group"]:
            groups.append((r["group"], []))
        groups[-1][1].append(r)

    fig, ax = plt.subplots(figsize=(8.2, 3.05))
    y = 0.0
    ticks, labels, headers = [], [], []
    for gname, items in groups:
        y -= 1.05
        headers.append((y, gname))
        for r in items:
            y -= 1.0
            a = float(r["ppv_coded_only"])
            b = float(r["ppv_codes_plus_text"])
            lo = min(float(r["ci_low_coded"]), float(r["ci_low_both"]))
            hi = max(float(r["ci_high_coded"]), float(r["ci_high_both"]))
            colour = NAVY if r["crosses_3pc"] == "yes" else GREY
            # pale band spans both intervals, so the reader sees that the
            # corrected haematuria CI (2.6-3.2) still includes 3
            ax.plot([lo, hi], [y, y], color=colour, lw=6.5, alpha=0.18,
                    solid_capstyle="butt", zorder=1)
            ax.annotate("", xy=(b, y), xytext=(a, y),
                        arrowprops=dict(arrowstyle="-|>", color=colour, lw=2.4,
                                        shrinkA=6, shrinkB=2, mutation_scale=15), zorder=2)
            ax.plot([a], [y], "o", ms=9, color=colour, mfc="white", mew=2.2, zorder=3)
            ax.plot([b], [y], "o", ms=9, color=colour, zorder=3)
            if abs(a - b) < 1e-9:
                ax.text(hi * 1.16, y, f"{a}%, unchanged", color=colour,
                        ha="left", va="center", fontsize=11.5)
            else:
                ax.text(lo * 0.84, y, f"{b}%", color=colour, ha="right", va="center",
                        fontsize=11.5, fontweight="bold")
                ax.text(hi * 1.16, y, f"{a}%", color=colour, ha="left", va="center", fontsize=11.5)
            ticks.append(y)
            labels.append(f"{r['symptom']}, {r['cancer']}")
    ybot = y - 0.72

    for yy, gname in headers:
        ax.text(-0.475, yy, gname, transform=ax.get_yaxis_transform(), ha="left",
                va="center", fontsize=11.5, fontweight="bold", color="#2D2D2D", clip_on=False)
    # the dashed line must stop below the label, or it strikes through the text
    label_y = 0.34
    ax.plot([3.0, 3.0], [ybot, label_y - 0.06], color=NAVY, ls=(0, (4, 3)),
            lw=1.7, zorder=1)
    thresh_label = ax.text(3.0, label_y, "3% referral threshold", color=NAVY,
                           fontsize=11.5, fontweight="bold", ha="center",
                           va="bottom")
    ax.set_xscale("log")
    ax.set_xlim(0.085, 34)
    ax.set_ylim(ybot, 1.18)
    ax.set_xticks([0.1, 0.3, 1, 3, 10, 30])
    ax.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30"], fontsize=11)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=11.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Positive predictive value for cancer (%), log scale", fontsize=11.5, labelpad=3)
    fig.text(0.5, 0.032,
             "open circle: coded records only          filled circle: codes plus recovered text"
             "          pale band: span of both 95% CIs",
             ha="center", fontsize=10, color="#6A6A6A")
    fig.subplots_adjust(left=0.325, right=0.885, top=0.975, bottom=0.255)
    check_no_overlap(fig, ax)
    check_line_clear_of_text(fig, ax, thresh_label, x_data=3.0)
    fig.savefig(OUT / "fig_threshold.png", facecolor="white")
    plt.close(fig)
    return len(raw)


def pooling():
    """Subgroup calibration slopes and the pooled value that describes none of them.

    Values are read from data/calibration_slopes.csv, transcribed from the
    results table of Lexiyao/model-evaluation-from-scratch. All THREE subgroups
    are shown: with only the two extremes on screen a reader who averages them
    gets 1.18, not the reported pooled 1.60, and the figure looks wrong.
    """
    rows_ = read_rows("calibration_slopes.csv")
    subs = [r for r in rows_ if r["is_pooled"] == "no"]
    pooled = next(r for r in rows_ if r["is_pooled"] == "yes")
    ordered = subs + [pooled]

    fig, ax = plt.subplots(figsize=(3.46, 2.28))
    ys = list(range(len(ordered) - 1, -1, -1))
    for r, yv in zip(ordered, ys):
        xv = float(r["calibration_slope"])
        is_pooled = r["is_pooled"] == "yes"
        colour = NAVY if is_pooled else "#2D2D2D"
        ax.plot([1.0, xv], [yv, yv], color=colour, lw=1.6, alpha=0.40, zorder=1)
        ax.plot([xv], [yv], "D" if is_pooled else "o", ms=9 if is_pooled else 8,
                color=colour, zorder=3)
        ax.text(xv, yv + 0.26, f"{xv:.2f}", fontsize=10, color=colour, ha="center",
                va="bottom", fontweight="bold" if is_pooled else "normal")
    ax.axvline(1.0, color=GREY, ls=(0, (3, 2.2)), lw=1.2, zorder=1)
    ax.text(1.0, max(ys) + 0.72, "1.00 = perfect", fontsize=9.5, color=GREY,
            ha="center", va="center")
    ax.axhline(0.55, color=RULE, lw=1.0, zorder=1)   # separates pooled from subgroups
    ax.set_xlim(0.30, 2.16)
    ax.set_ylim(-0.75, max(ys) + 1.05)
    ax.set_yticks(ys)
    ax.set_yticklabels([r["group"] for r in ordered], fontsize=9.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xticks([0.5, 1.0, 1.5, 2.0])
    ax.set_xticklabels(["0.50", "1.00", "1.50", "2.00"], fontsize=9.5)
    ax.set_xlabel("Calibration slope", fontsize=9.5, labelpad=1)
    ax.text(float(pooled["calibration_slope"]), -0.42, "describes no one", fontsize=9.5,
            color=NAVY, fontweight="bold", ha="center", va="center")
    fig.subplots_adjust(left=0.395, right=0.985, top=0.985, bottom=0.215)
    check_no_overlap(fig, ax)
    fig.savefig(OUT / "fig_pooling.png", facecolor="white")
    plt.close(fig)
    return [float(r["calibration_slope"]) for r in subs], float(pooled["calibration_slope"])


def check_no_overlap(fig, ax):
    """Fail loudly if any two visible text boxes collide, or text sits on a spine.

    A figure that silently overlaps its own labels is a wrong figure. This is a
    geometric check only: it cannot see a leader line pointing at the wrong
    row, so anchor annotations to looked-up coordinates instead of eyeballed
    ones (see forest(), which looks up the 'Hip pain' row).
    """
    import matplotlib as mpl
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(renderer)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    canvas = fig.bbox
    for t, box in texts:
        if box.x0 < canvas.x0 - 1 or box.x1 > canvas.x1 + 1:
            problems_outside = f"{t.get_text()[:24]!r} runs outside the canvas"
            raise AssertionError(problems_outside)
    ticklabels = set(ax.get_xticklabels() + ax.get_yticklabels())
    problems = []
    for i, (a, box_a) in enumerate(texts):
        for b, box_b in texts[i + 1:]:
            if box_a.overlaps(box_b):
                problems.append(f"{a.get_text()[:24]!r} overlaps {b.get_text()[:24]!r}")
    for t, box in texts:
        if t in ticklabels:
            continue
        for spine in ax.spines.values():
            if spine.get_visible() and box.overlaps(spine.get_window_extent(renderer)):
                problems.append(f"{t.get_text()[:24]!r} sits on a spine")
    if problems:
        raise AssertionError("figure has colliding labels:\n  " + "\n  ".join(problems))


if __name__ == "__main__":
    house_style()
    n_above, n_total = forest()
    n_pairs = threshold()
    slopes, pooled = pooling()
    print(f"fig_forest.png      {n_above} of {n_total} associations above 1")
    print(f"fig_threshold.png   {n_pairs} symptom/cancer pairs, log scale")
    print(f"fig_pooling.png     subgroups {slopes} pool to {pooled}")
    print(f"all figures in {OUT}")
