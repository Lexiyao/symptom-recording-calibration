"""
Figures for the DPhil interview presentation:
"Differential symptom recording and the calibration of cancer risk
prediction models in UK primary care"

Reproduces three figures from the CSVs in ../data:

  fig_forest.png     — d'Elia et al. (2025) Table 3: odds of a symptom being
                       coded, by ethnicity and deprivation (CPRD Aurum)
  fig_threshold.png  — Price et al. (2016): positive predictive values before
                       and after recovering free-text records
  fig_pooling.png    — illustration of subgroup calibration slopes pooling to
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
    """d'Elia Table 3 — odds that a symptom was coded, by patient group."""
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
    ax.text(3.30, ticks[1], f"{n_above} of {n_total} associations\nlie above 1",
            fontsize=12, color=NAVY, fontweight="bold", ha="left", va="center",
            linespacing=1.40)

    fig.subplots_adjust(left=0.295, right=0.985, top=0.99, bottom=0.115)
    check_no_overlap(fig, ax)
    fig.savefig(OUT / "fig_forest.png", facecolor="white")
    plt.close(fig)
    return n_above, n_total


def threshold():
    """Price 2016 — PPV with coded records only vs codes plus recovered text.

    Drawn as small multiples with a per-panel scale, deliberately. On a shared
    0-15% axis the jaundice movement (12.8 to 6.3) occupies 43% of the axis
    while haematuria (4.0 to 2.9) occupies 7% - so the visually dominant series
    is the secondary finding and the one that actually crosses the referral
    threshold is nearly invisible. Per-panel scales give each movement the same
    visual weight; the axis labels carry the magnitude.
    """
    raw = read_rows("price2016_ppv.csv")
    # per-panel x-range and ticks, chosen so each movement fills its panel
    VIEW = {"Haematuria": ((2.42, 4.58), [3, 4], "in"),
            "Jaundice": ((4.15, 15.1), [6, 9, 12], "below"),
            "Abdominal pain": ((0.1775, 0.2625), [0.20, 0.24], "above")}
    NOTE = {"in": "3% referral threshold (dashed)",
            "below": "3% threshold lies below this range",
            "above": "3% threshold lies above this range"}
    VERDICT = {"in": "crosses the threshold", "below": "halves, but stays above",
               "above": "does not move"}

    fig, axes = plt.subplots(1, 3, figsize=(8.6, 2.72))
    Y, L0, L1 = 0.50, 0.30, 0.70
    for ax, r in zip(axes, raw):
        a = float(r["ppv_coded_only"])
        b = float(r["ppv_codes_plus_text"])
        xr, ticks, thr = VIEW[r["symptom"]]
        colour = NAVY if r["crosses_3pc"] == "yes" else GREY
        ax.set_xlim(*xr)
        ax.set_ylim(0, 1)
        if thr == "in":
            # a short vertical segment, not axvline: a full-height line runs
            # through the value labels and the verdict text
            ax.vlines(3.0, L0, L1, color=NAVY, ls=(0, (4, 3)), lw=1.7, zorder=1)
        ax.annotate("", xy=(b, Y), xytext=(a, Y),
                    arrowprops=dict(arrowstyle="-|>", color=colour, lw=2.6,
                                    shrinkA=5, shrinkB=1, mutation_scale=17))
        ax.plot([a], [Y], "o", ms=10, color=colour, mfc="white", mew=2.4, zorder=3)
        ax.plot([b], [Y], "o", ms=10, color=colour, zorder=3)
        span = xr[1] - xr[0]
        ax.text(b - span * 0.055, Y, f"{b}%", color=colour, ha="right", va="center",
                fontsize=12, fontweight="bold")
        ax.text(a + span * 0.055, Y, f"{a}%", color=colour, ha="left", va="center", fontsize=12)
        ax.set_title(f"{r['symptom']}\n{r['cancer'].lower()} cancer", fontsize=12.5,
                     color="#2D2D2D", pad=8, linespacing=1.35)
        ax.set_yticks([])
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t:g}" for t in ticks], fontsize=11)
        ax.tick_params(axis="x", length=3, colors="#6A6A6A")
        ax.text(0.5, 0.955, NOTE[thr], transform=ax.transAxes, color=NAVY, fontsize=10.5,
                ha="center", va="top", fontweight="bold" if thr == "in" else "normal")
        ax.text(0.5, 0.085, VERDICT[thr], transform=ax.transAxes, color=colour, fontsize=11.5,
                ha="center", va="center", fontweight="bold" if thr == "in" else "normal")

    fig.text(0.5, 0.075, "Positive predictive value for cancer (%) \u2014 each panel on its own scale",
             ha="center", fontsize=11.5, color="#4A4A4A")
    fig.text(0.5, 0.012,
             "open circle: coded records only          filled circle: codes plus recovered free text",
             ha="center", fontsize=10.5, color="#6A6A6A")
    fig.subplots_adjust(left=0.055, right=0.978, top=0.775, bottom=0.235, wspace=0.20)
    check_panels_clean(fig, axes, dashed_at=(axes[0], 3.0, L0, L1))
    fig.savefig(OUT / "fig_threshold.png", facecolor="white")
    plt.close(fig)


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


def check_panels_clean(fig, axes, dashed_at=None):
    """Multi-panel version of check_no_overlap, plus a dashed-line check.

    dashed_at=(ax, x, y0, y1) asserts that the short dashed segment drawn at
    data-x on that axis does not pass through any text. The plain text-vs-text
    check cannot catch this: a Line2D is not a Text, so a line running straight
    through four labels passes silently.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    items = []
    for ax in axes:
        for t in list(ax.texts) + [ax.title] + ax.get_xticklabels():
            if t.get_text().strip():
                items.append((ax, t, t.get_window_extent(r)))
    problems = []
    for i, (_, a, ba) in enumerate(items):
        for _, b, bb in items[i + 1:]:
            if ba.overlaps(bb):
                problems.append(f"{a.get_text()[:22]!r} overlaps {b.get_text()[:22]!r}")
    fb = fig.bbox
    for _, t, bb in items:
        if bb.x0 < fb.x0 or bb.x1 > fb.x1 or bb.y0 < fb.y0 or bb.y1 > fb.y1:
            problems.append(f"{t.get_text()[:22]!r} runs outside the canvas")
    if dashed_at:
        ax, xdata, y0, y1 = dashed_at
        xpix = ax.transData.transform((xdata, 0))[0]
        ylo = ax.transAxes.transform((0, y0))[1]
        yhi = ax.transAxes.transform((0, y1))[1]
        for a_, t, bb in items:
            if a_ is ax and bb.x0 <= xpix <= bb.x1 and bb.y1 >= ylo and bb.y0 <= yhi:
                problems.append(f"dashed line crosses {t.get_text()[:22]!r}")
    if problems:
        raise AssertionError("figure has layout faults:\n  " + "\n  ".join(problems))


def check_no_overlap(fig, ax):
    """Fail loudly if any two visible text boxes collide, or text sits on a spine.

    A figure that silently overlaps its own labels is a wrong figure. This is a
    geometric check only — it cannot see a leader line pointing at the wrong
    row, so anchor annotations to looked-up coordinates rather than eyeballed
    ones (see forest(), which looks up the 'Hip pain' row).
    """
    import matplotlib as mpl
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(renderer)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
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
    threshold()
    slopes, pooled = pooling()
    print(f"fig_forest.png      {n_above} of {n_total} associations above 1")
    print("fig_threshold.png   3 panels, each on its own scale")
    print(f"fig_pooling.png     subgroups {slopes} pool to {pooled}")
    print(f"all figures in {OUT}")
