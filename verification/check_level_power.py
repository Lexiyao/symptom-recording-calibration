"""Precision of calibration-in-the-large, the PRIMARY estimand of Aim 3.

Appendix Z of the interview deck quotes the numbers this script prints, and
slide 5 quotes the 20,000-per-group row. verification/check_null_power.py does
the same job for the calibration SLOPE, which is the secondary bound; the two
scripts share their design assumptions and differ only in the estimand.

Why a separate script. The two-group simulation (symptom-recording-simulation,
R/simulate_two_group.R) found the recording gap lands on the calibration LEVEL
and not on the gradient: across nine grid points the slope gap ran -0.018 to
+0.003 against a Monte Carlo SE of 0.008. The estimand moved accordingly, so
the power calculation had to move with it.

Estimand. Calibration-in-the-large is the recalibration intercept from
    logit(P(Y=1)) = a + offset(lp)
with the model's linear predictor held at coefficient one. Exponentiating a
gives the observed-to-expected risk ratio, which is what the deck reports as
an O:E gap.

Run:  python verification/check_level_power.py
"""

import math

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

EVENT_RATE = 0.015     # two-year lung + bowel cancer incidence, order of magnitude
LP_SD = 1.25           # SD of the linear predictor; same declared value as check_null_power.py
ALPHA = 0.05
POWER = 0.80
SIZES = (5_000, 20_000, 100_000)
THRESHOLD = 0.03       # NICE NG12 referral threshold
N_SIM = 400            # replicates for the simulation cross-check
SEED = 20260726

# what Appendix Z and slide 5 claim, asserted below
CLAIMED_SE = {5_000: 0.119, 20_000: 0.059, 100_000: 0.027}
CLAIMED_OE = {5_000: 1.60, 20_000: 1.27, 100_000: 1.11}
CLAIMED_PCT = (3.37, 2.67)   # the 20,000 row, expressed at the referral line

Z = norm.ppf(1 - ALPHA / 2) + norm.ppf(POWER)


def _grid(n=200_001):
    g = np.linspace(-9, 9, n)
    w = np.exp(-0.5 * g ** 2)
    return g, w / w.sum()


def intercept_for_rate(rate, sd):
    """Intercept that makes E[expit(mu + sd*Z)] equal the target event rate."""
    g, w = _grid()
    f = lambda mu: float((w / (1 + np.exp(-(mu + sd * g)))).sum()) - rate
    return brentq(f, -15, 5, xtol=1e-14)


def fisher_se_intercept(n, sd, rate):
    """Closed-form SE of the recalibration intercept with the slope fixed at one.

    The offset model has a single parameter, so the Fisher information is
    n * E[p(1-p)] and the SE is its inverse square root. Using pbar(1-pbar)
    instead of E[p(1-p)] overstates the information, because p(1-p) is concave:
    at a 1.5% rate with SD 1.25 the naive version gives 0.0577 rather than
    0.0594, which is the difference between an O:E gap of 1.26 and 1.27.
    """
    mu = intercept_for_rate(rate, sd)
    g, w = _grid()
    p = 1 / (1 + np.exp(-(mu + sd * g)))
    return 1.0 / math.sqrt(n * float((w * p * (1 - p)).sum()))


def simulated_se_intercept(n, sd, rate, n_sim=N_SIM, seed=SEED):
    """Same quantity by simulation: fit the offset model, keep the intercept."""
    rng = np.random.default_rng(seed)
    mu = intercept_for_rate(rate, sd)
    out = []
    for _ in range(n_sim):
        lp = mu + sd * rng.normal(size=n)
        y = rng.binomial(1, 1 / (1 + np.exp(-lp)))
        a = 0.0
        for _ in range(60):                       # Newton-Raphson on the intercept
            p = 1 / (1 + np.exp(-(a + lp)))
            grad = float((y - p).sum())
            hess = float((p * (1 - p)).sum())
            if hess <= 0:
                break
            step = grad / hess
            a += step
            if abs(step) < 1e-12:
                break
        out.append(a)
    return float(np.std(out, ddof=1))


def at_threshold(ratio, threshold=THRESHOLD):
    """Split an O:E ratio symmetrically about the referral threshold."""
    return threshold * math.sqrt(ratio) * 100, threshold / math.sqrt(ratio) * 100


def main():
    print("1. PRECISION OF CALIBRATION-IN-THE-LARGE (the primary estimand)")
    print(f"   linear predictor ~ N({intercept_for_rate(EVENT_RATE, LP_SD):.3f}, {LP_SD}^2)"
          f" gives a {EVENT_RATE:.1%} event rate")
    print(f"   recalibration intercept, slope held at 1; {N_SIM} replicates per size\n")
    print(f"{'n':>11} {'events':>8} {'SE (sim)':>10} {'SE (Fisher)':>12} {'claimed':>8}")
    ses = {}
    for n in SIZES:
        fisher = fisher_se_intercept(n, LP_SD, EVENT_RATE)
        sim = simulated_se_intercept(n, LP_SD, EVENT_RATE)
        ses[n] = fisher
        print(f"{n:>11,} {n * EVENT_RATE:>8.0f} {sim:>10.3f} {fisher:>12.3f} {CLAIMED_SE[n]:>8.3f}")
        assert abs(sim - fisher) < 0.004, f"simulation {sim:.4f} vs closed form {fisher:.4f} at n={n}"
        assert abs(fisher - CLAIMED_SE[n]) < 0.001, (
            f"appendix Z claims SE {CLAIMED_SE[n]} at n={n}, closed form gives {fisher:.4f}")
    print("\n   Simulation and closed form agree, so the standard errors are a")
    print("   property of the design and not an artefact of one random draw.")

    print("\n2. WHAT THAT DETECTS, BETWEEN TWO SUBGROUPS")
    print(f"   two-sided alpha {ALPHA}, power {POWER:.0%}:"
          f" MDE = {Z:.3f} x SE x sqrt(2) = {Z * math.sqrt(2):.3f} x SE\n")
    print(f"{'n per group':>13} {'SE':>8} {'MDE (log-odds)':>16} {'as O:E ratio':>14} {'claimed':>8}")
    for n in SIZES:
        mde = Z * math.sqrt(2) * ses[n]
        ratio = math.exp(mde)
        print(f"{n:>13,} {ses[n]:>8.3f} {mde:>16.3f} {ratio:>14.3f} {CLAIMED_OE[n]:>8.2f}")
        assert abs(round(ratio, 2) - CLAIMED_OE[n]) < 0.005, (
            f"appendix Z claims O:E {CLAIMED_OE[n]} at n={n}, this gives {ratio:.4f}")

    print("\n3. THE 20,000 ROW AT THE REFERRAL LINE")
    ratio = math.exp(Z * math.sqrt(2) * ses[20_000])
    hi, lo = at_threshold(ratio)
    print(f"   an O:E gap of {ratio:.4f} ({ratio:.2f} to two places), split about the {THRESHOLD:.0%} line,")
    print(f"   is {hi:.2f}% against {lo:.2f}% -- the difference between referring and not.")
    print(f"   claimed on slide 5: {CLAIMED_PCT[0]}% against {CLAIMED_PCT[1]}%")
    assert abs(hi - CLAIMED_PCT[0]) < 0.005 and abs(lo - CLAIMED_PCT[1]) < 0.005, (
        f"slide 5 claims {CLAIMED_PCT}, this gives ({hi:.3f}, {lo:.3f})")
    print(f"\n   Note on rounding: {hi:.4f} and {lo:.4f} each round to two places as")
    print(f"   quoted, but {CLAIMED_PCT[0]}/{CLAIMED_PCT[1]} = {CLAIMED_PCT[0] / CLAIMED_PCT[1]:.4f}."
          f" The ratio {ratio:.2f} is computed before")
    print("   rounding, not from the rounded percentages.")

    print("\n4. HOW MUCH THE ANSWER DEPENDS ON THE DECLARED SD")
    print(f"{'LP SD':>8} {'SE @20,000':>12} {'O:E detectable':>16}")
    for sd in (1.00, 1.25, 1.50):
        se = fisher_se_intercept(20_000, sd, EVENT_RATE)
        print(f"{sd:>8.2f} {se:>12.4f} {math.exp(Z * math.sqrt(2) * se):>16.3f}")
    print("\n   A flatter risk distribution carries slightly MORE information about")
    print("   the level (unlike the slope, where it carries less), because p(1-p)")
    print("   is larger when risk is concentrated near the mean. The SD is a")
    print("   declared assumption either way: Aim 3 reports the observed SD of")
    print("   the linear predictor next to the power.")

    print("\nAll assertions passed.")


if __name__ == "__main__":
    main()
