"""
What the simulation in Lexiyao/symptom-recording-simulation does and does not
establish.

Run:  python verification/check_simulation.py

The simulation uses synthetic data with every parameter declared, so it cannot
support any claim about the size of an effect in real records. This script
separates the three things it can support, and shows the arithmetic for each.

  1. A code check. At mnar_strength = 0 the coding indicator is independent of
     true symptom status, which is ordinary non-differential misclassification
     of a binary exposure with sensitivity < 1 and specificity = 1. That has a
     closed form. If the simulation did not reproduce it, the simulation would
     be wrong. It does reproduce it, to within Monte Carlo error, and that is
     all this agreement means.

  2. The content. Sweeping mnar_strength has no simple closed form, and the two
     quantities move apart: attenuation of the odds ratio improves markedly
     while calibration among true symptom carriers barely improves and never
     reaches 1. A study reporting only the first would conclude the problem is
     shrinking.

  3. A falsified expectation. The proposal assumed the damage came mainly from
     recording that tracks unobserved severity. Most of the attenuation is
     already present when coding is at random, so the proposal's target moved
     from the odds ratio to subgroup calibration.

Declared parameters, from sim_params() in R/simulate_recording_bias.R:
     n = 20000, prev_symptom = 0.08, baseline_risk = 0.004,
     log_or_symptom = 2.34, base_coding_prob = 0.85 (alarm) or 0.45 (vague),
     threshold = 0.03
"""

import csv
import math
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

PREV_SYMPTOM = 0.08
BASELINE_RISK = 0.004
LOG_OR_SYMPTOM = 2.34
U_COEF = 0.35          # coefficient on unobserved severity in the outcome model
MC_TOLERANCE = 0.6     # percentage points; Monte Carlo SE is ~1 point at nsim=200


def marginal_risk(symptom, n_quad=20001):
    """P(Y=1 | S=symptom), integrating the unobserved severity U ~ N(0,1)."""
    intercept = math.log(BASELINE_RISK / (1 - BASELINE_RISK))
    total_w = 0.0
    total = 0.0
    for i in range(n_quad):
        u = (i / (n_quad - 1)) * 12 - 6
        w = math.exp(-u * u / 2)
        eta = intercept + LOG_OR_SYMPTOM * symptom + U_COEF * u
        total += w / (1 + math.exp(-eta))
        total_w += w
    return total / total_w


def closed_form(sensitivity):
    """Observed OR and O/E under non-differential misclassification.

    An uncoded symptom reads as absent, so sensitivity = P(coded | present).
    A symptom that was never present is never coded, so specificity = 1.
    """
    p1, p0 = marginal_risk(1), marginal_risk(0)
    or_true = (p1 / (1 - p1)) / (p0 / (1 - p0))

    a = PREV_SYMPTOM * sensitivity * p1
    b = PREV_SYMPTOM * sensitivity * (1 - p1)
    c = PREV_SYMPTOM * (1 - sensitivity) * p1 + (1 - PREV_SYMPTOM) * p0
    d = PREV_SYMPTOM * (1 - sensitivity) * (1 - p1) + (1 - PREV_SYMPTOM) * (1 - p0)
    or_obs = (a / b) / (c / d)

    # among TRUE symptom carriers, the model's prediction is a mixture: those
    # coded get the exposed-group risk, those uncoded get the unexposed-group risk
    w_sym = PREV_SYMPTOM * (1 - sensitivity)
    w_asym = 1 - PREV_SYMPTOM
    r_unexposed = (w_sym * p1 + w_asym * p0) / (w_sym + w_asym)
    predicted = sensitivity * p1 + (1 - sensitivity) * r_unexposed
    return or_true, or_obs, (1 - or_obs / or_true) * 100, p1 / predicted


def load_sweep():
    path = DATA / "sweep_results.csv"
    if not path.exists():
        raise SystemExit(f"missing {path}. Copy figures/sweep_results.csv from "
                         "Lexiyao/symptom-recording-simulation into data/.")
    rows = list(csv.DictReader(open(path, newline="")))
    for r in rows:
        for k, v in r.items():
            if k != "regime":
                r[k] = float(v)
    return rows


def check_code_against_theory(rows):
    print("1. CODE CHECK: simulated attenuation at random coding vs closed form")
    ok = True
    for label, coding_prob in (("vague", 0.45), ("alarm", 0.85)):
        at_zero = next(r for r in rows
                       if r["base_coding_prob"] == coding_prob and r["mnar_strength"] == 0)
        simulated = (1 - at_zero["or_standard"] / at_zero["or_oracle"]) * 100
        _, _, analytic, oe_analytic = closed_form(coding_prob)
        gap = abs(simulated - analytic)
        print(f"   {label:6s} sensitivity {coding_prob}:  simulated {simulated:5.1f}%   "
              f"closed form {analytic:5.1f}%   difference {gap:.1f} points")
        print(f"          observed/predicted:  simulated {at_zero['cal_oe']:.3f}   "
              f"closed form {oe_analytic:.3f}")
        if gap > MC_TOLERANCE:
            ok = False
    assert ok, "simulation departs from the closed-form value by more than Monte Carlo error"
    print("   Agreement means the code is correct. It is not a result.")


def check_divergence(rows):
    print("\n2. CONTENT: the two quantities move apart across the sweep")
    vague = sorted((r for r in rows if r["base_coding_prob"] == 0.45),
                   key=lambda r: r["mnar_strength"])
    first, last = vague[0], vague[-1]
    att_first = (1 - first["or_standard"] / first["or_oracle"]) * 100
    att_last = (1 - last["or_standard"] / last["or_oracle"]) * 100
    print(f"   attenuation:          {att_first:5.1f}%  ->  {att_last:5.1f}%   "
          f"(improves by {att_first - att_last:.1f} points)")
    print(f"   observed/predicted:   {first['cal_oe']:5.3f}  ->  {last['cal_oe']:5.3f}   "
          f"(improves by {first['cal_oe'] - last['cal_oe']:.3f})")
    print(f"   best observed/predicted is still {(last['cal_oe'] - 1) * 100:.0f}% above 1")
    assert att_first - att_last > 15, "expected the odds ratio to recover across the sweep"
    assert last["cal_oe"] > 1.4, "expected calibration to stay far from 1"
    print("   A study reporting only the odds ratio would call the problem shrinking.")


def check_falsified_expectation(rows):
    print("\n3. FALSIFIED EXPECTATION: where the damage actually comes from")
    vague = [r for r in rows if r["base_coding_prob"] == 0.45]
    at_zero = next(r for r in vague if r["mnar_strength"] == 0)
    worst = max(vague, key=lambda r: (1 - r["or_standard"] / r["or_oracle"]))
    att_zero = (1 - at_zero["or_standard"] / at_zero["or_oracle"]) * 100
    att_worst = (1 - worst["or_standard"] / worst["or_oracle"]) * 100
    print(f"   attenuation at random coding:  {att_zero:.1f}%")
    print(f"   worst attenuation anywhere:    {att_worst:.1f}% "
          f"(at mnar_strength = {worst['mnar_strength']})")
    assert worst["mnar_strength"] == 0, \
        "the proposal assumed informative recording would be worse; it is not"
    print("   The proposal assumed informative recording was the main mechanism.")
    print("   It is not, so the target moved to subgroup calibration.")


if __name__ == "__main__":
    rows = load_sweep()
    print(f"loaded {len(rows)} sweep rows "
          f"({int(rows[0]['nsim'])} replicates per point)\n")
    check_code_against_theory(rows)
    check_divergence(rows)
    check_falsified_expectation(rows)
    print("\nAll assertions passed. Nothing here is a claim about real patients.")
