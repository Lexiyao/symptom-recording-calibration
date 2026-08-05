"""
Consistency checks on the published values used in the presentation.

Run:  python verification/check_sources.py

Two kinds of check:

  A. Internal consistency of each published confidence interval — a 99.76% CI on
     an odds ratio should be roughly symmetric on the log scale, so its
     geometric centre sqrt(lower*upper) should sit close to the point estimate.
     An interval that fails this is either a typo or a different quantity.

  B. Threshold arithmetic — which predictive values cross the 3% line that NICE
     used when building the urgent referral criteria.

Check A found three discrepancies between the ABSTRACT of d'Elia et al. (2025)
and Table 3 of the same paper. The presentation uses the Table 3 values
throughout. The discrepancies are reported below, and on an appendix slide, not
because they matter to the argument but because checking that an interval is
consistent with its own point estimate is the kind of check the symptom code
lists in Aim 1 will need.
"""

import csv
import math
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

# values as printed in the ABSTRACT of d'Elia et al. (2025), for comparison
# against Table 3 (which is what data/delia2025_table3.csv holds)
ABSTRACT_AS_PRINTED = [
    ("Muscle cramps", "South Asian", 1.71, 1.44, 2.57),
    ("Falls", "Most deprived (IMD 1)", 1.37, 2.03, 1.82),
    ("Shoulder pain", "South Asian", 1.44, 1.25, 1.66),
]

NICE_THRESHOLD = 3.0


def log_centred(point, lower, upper, tolerance=0.02):
    """Is the interval's geometric centre within tolerance of the point estimate?"""
    if lower <= 0 or upper <= 0:
        return False, float("nan")
    centre = math.sqrt(lower * upper)
    return abs(centre - point) / point <= tolerance, centre


def check_table3():
    print("A. Internal consistency of d'Elia Table 3 (as used in the deck)")
    rows = list(csv.DictReader(open(DATA / "delia2025_table3.csv", newline="")))
    failures = 0
    for r in rows:
        point = float(r["odds_ratio"])
        lo, hi = float(r["ci_low_99_76"]), float(r["ci_high_99_76"])
        ok, centre = log_centred(point, lo, hi)
        if hi < lo:
            print(f"   FAIL  {r['symptom']:22s} upper limit below lower limit")
            failures += 1
        elif not ok:
            print(f"   FAIL  {r['symptom']:22s} OR {point} vs geometric centre {centre:.2f}")
            failures += 1
    print(f"   {len(rows) - failures}/{len(rows)} intervals internally consistent")
    assert failures == 0, "Table 3 values should all be self-consistent"
    return rows


def check_abstract():
    print("\nB. The same three associations as printed in the paper's ABSTRACT")
    for symptom, group, point, lo, hi in ABSTRACT_AS_PRINTED:
        ok, centre = log_centred(point, lo, hi)
        if hi < lo:
            verdict = f"limits reversed (lower {lo} > upper {hi})"
        elif not ok:
            verdict = f"implies OR {centre:.2f}, not {point}"
        else:
            verdict = "consistent, but attributed to a different group in Table 3"
        print(f"   {symptom:16s} {point} ({lo}\u2013{hi})  \u2192 {verdict}")
    print("   The deck uses the Table 3 values in every case.")


def check_thresholds():
    print(f"\nC. Which predictive values cross the {NICE_THRESHOLD}% referral threshold")
    rows = list(csv.DictReader(open(DATA / "price2016_ppv.csv", newline="")))
    for r in rows:
        before = float(r["ppv_coded_only"])
        after = float(r["ppv_codes_plus_text"])
        crosses = (before > NICE_THRESHOLD) != (after > NICE_THRESHOLD)
        claimed = r["crosses_3pc"] == "yes"
        assert crosses == claimed, f"{r['symptom']}: crossing flag disagrees with the values"
        fold = before / after if after else float("nan")
        note = "CROSSES the threshold" if crosses else f"stays on one side ({fold:.2f}-fold move)"
        print(f"   {r['symptom']:16s} {before}% \u2192 {after}%   {note}")
    n = sum(1 for r in rows if r["crosses_3pc"] == "yes")
    print(f"   {n} of {len(rows)} symptoms cross the line when text records are recovered")


def check_barclay_bands():
    print("\nD. Single-site risk measured off Barclay Fig 2 (see code/barclay_panels.py)")
    rows = list(csv.DictReader(open(DATA / "barclay2024_measured_bands.csv", newline="")))
    for panel in ("Haematuria", "Fatigue"):
        bands = [r for r in rows if r["panel"] == panel]
        top = max(bands, key=lambda r: float(r["band_thickness_pc_at_oldest_age"]))
        total = sum(float(r["band_thickness_pc_at_oldest_age"]) for r in bands)
        above = [r["cancer_site"] for r in bands
                 if float(r["band_thickness_pc_at_oldest_age"]) > NICE_THRESHOLD]
        print(f"   {panel:12s} largest single site {top['cancer_site']} "
              f"{float(top['band_thickness_pc_at_oldest_age']):.2f}%  "
              f"| stacked total {total:.2f}%  | sites above 3%: {above or 'none'}")
    print("   Barclay's claim is about risk at INDIVIDUAL sites, so the stacked")
    print("   total is not the quantity of interest; the annotation reports a single site.")


def check_subgroup_power():
    print("\nE. Feasibility floor for estimating a subgroup calibration slope")
    print("   Rule of thumb: ~100 events minimum per subgroup, 200 preferred.")
    for share in (0.02, 0.05, 0.10, 0.20):
        print(f"   subgroup = {share * 100:4.0f}% of cohort \u2192 need "
              f"{100 / share:7,.0f} total events (floor), {200 / share:7,.0f} (preferred)")
    print("   This is why Aim 2 prespecifies which subgroups can be powered, and")
    print("   reports subgroup sizes before any calibration comparison.")


if __name__ == "__main__":
    check_table3()
    check_abstract()
    check_thresholds()
    check_barclay_bands()
    check_subgroup_power()
    print("\nAll assertions passed.")
