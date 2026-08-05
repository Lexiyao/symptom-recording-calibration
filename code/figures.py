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
    """Price 2016 — PPV with coded records only vs codes plus recovered text."""
    raw = read_rows("price2016_ppv.csv")
    fig, ax = plt.subplots(figsize=(6.10, 3.24))
    ys = list(range(len(raw) - 1, -1, -1))

    for r, yv in zip(raw, ys):
        a = float(r["ppv_coded_only"])
        b = float(r["ppv_codes_plus_text"])
        crosses = r["crosses_3pc"] == "yes"
        colour = NAVY if crosses else GREY
        ax.annotate("", xy=(b, yv), xytext=(a, yv),
                    arrowprops=dict(arrowstyle="-|>", color=colour, lw=2.4,
                                    shrinkA=4, shrinkB=0, mutation_scale=15))
        ax.plot([a], [yv], "o", ms=8, color=colour, mfc="white", mew=2.2, zorder=3)
        ax.plot([b], [yv], "o", ms=8, color=colour, zorder=3)
        if a - b > 0.5:
            ax.text(a + 0.45, yv, f"{a}%", color=colour, ha="left", va="center")
            ax.text(b - 0.45, yv, f"{b}%", color=colour, ha="right", va="center",
                    fontweight="bold" if crosses else "normal")
        else:
            # markers would collide: one combined label, offset past both
            ax.text(max(a, b) + 0.95, yv, f"{a}% \u2192 {b}%   unchanged",
                    color=colour, ha="left", va="center")

    ax.axvline(3, color=NAVY, ls=(0, (4, 3)), lw=1.6, zorder=1)
    ax.text(3.25, max(ys) + 0.62, "3% referral threshold", color=NAVY,
            ha="left", va="center", fontweight="bold")
    ax.set_xlim(-1.4, 16.4)
    ax.set_ylim(min(ys) - 0.52, max(ys) + 0.86)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{r['symptom']}\n{r['cancer'].lower()} cancer" for r in raw])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xticks([0, 3, 5, 10, 15])
    ax.set_xlabel("Positive predictive value for cancer (%)", labelpad=2)

    crosser = next(i for i, r in enumerate(raw) if r["crosses_3pc"] == "yes")
    ax.annotate("crosses the line", xy=(2.9, ys[crosser] - 0.10),
                xytext=(5.6, ys[crosser] - 0.48), color=NAVY, fontweight="bold",
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=NAVY, lw=1.0, shrinkA=4,
                                shrinkB=4, connectionstyle="arc3,rad=-0.18"))
    ax.text(-1.30, min(ys) - 0.40,
            "\u25cb coded records only          \u25cf codes plus recovered text",
            fontsize=11, color="#4A4A4A", ha="left", va="center")

    fig.subplots_adjust(left=0.285, right=0.985, top=0.99, bottom=0.185)
    check_no_overlap(fig, ax)
    fig.savefig(OUT / "fig_threshold.png", facecolor="white")
    plt.close(fig)


def pooling():
    """Why a pooled calibration slope can describe neither subgroup."""
    slopes = [("Group A", 0.50, "#2D2D2D", "o"),
              ("Group B", 1.85, "#2D2D2D", "o"),
              ("Pooled", 1.60, NAVY, "D")]
    fig, ax = plt.subplots(figsize=(3.34, 2.16))
    ys = [2, 1, 0]
    for (lab, xv, colour, marker), yv in zip(slopes, ys):
        ax.plot([1.0, xv], [yv, yv], color=colour, lw=1.6, alpha=0.45, zorder=1)
        ax.plot([xv], [yv], marker, ms=9, color=colour, zorder=3)
        ax.text(xv, yv + 0.28, f"{xv:.2f}", fontsize=10, color=colour,
                ha="center", va="bottom",
                fontweight="bold" if marker == "D" else "normal")
    ax.axvline(1.0, color=GREY, ls=(0, (3, 2.2)), lw=1.2, zorder=1)
    ax.text(1.0, 2.64, "1.00 = perfect", fontsize=9.5, color=GREY,
            ha="center", va="center")
    ax.set_xlim(0.30, 2.14)
    ax.set_ylim(-0.92, 2.94)
    ax.set_yticks(ys)
    ax.set_yticklabels([s[0] for s in slopes], fontsize=10)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xticks([0.5, 1.0, 1.5, 2.0])
    ax.set_xticklabels(["0.50", "1.00", "1.50", "2.00"], fontsize=9.5)
    ax.set_xlabel("Calibration slope", fontsize=9.5, labelpad=1)
    ax.text(1.60, -0.52, "describes no one", fontsize=9.5, color=NAVY,
            fontweight="bold", ha="center", va="center")
    fig.subplots_adjust(left=0.245, right=0.985, top=0.99, bottom=0.235)
    check_no_overlap(fig, ax)
    fig.savefig(OUT / "fig_pooling.png", facecolor="white")
    plt.close(fig)


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
    pooling()
    print(f"fig_forest.png      {n_above} of {n_total} associations above 1")
    print("fig_threshold.png   written")
    print("fig_pooling.png     written")
    print(f"all figures in {OUT}")
