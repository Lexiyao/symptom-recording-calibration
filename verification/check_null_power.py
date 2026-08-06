"""
What a null result in Aim 3 would and would not establish.

Run:  python verification/check_null_power.py

Aim 3 compares the calibration slope between patient subgroups. If that
comparison comes back null, the result is only interpretable if the design could
have detected a difference worth detecting. This script computes what it could
have detected, so that a null is reported as a bound and never as an absence.

Three steps, in the order the appendix argues them:

  1. The precision. Simulate logistic recalibration at the event rate Aim 3
     expects and report the standard error of the fitted calibration slope at
     three subgroup sizes. The simulation is cross-checked against the Fisher
     information for the slope, which has a closed form once the intercept is
     profiled out. If the two disagreed, the simulation would be wrong.

  2. What that detects. Convert each standard error into the smallest
     between-group difference in slope detectable at 80% power and a two-sided
     5% level. Two independent groups, so the difference has variance 2 * se^2.

  3. Why 0.20 is the yardstick. Translate a slope back into the quantity a
     referral decision uses: the risk assigned to a patient the model places at
     the 3% threshold. This is where the parameterisation matters, and step 3
     shows both versions because only one of them answers the clinical question.

Every number here is a property of the design under declared assumptions. None
is an estimate from patient data, and no subgroup size below is a claim about
what CPRD will yield.

Declared assumptions:
     event rate 1.5%, linear predictor SD 1.25 on the log-odds scale,
     calibration slope 1.0 under the null, two-sided alpha 0.05, power 0.80

The linear predictor SD is a declared assumption and it drives the answer, so it
is stated rather than buried: at SD 1.25 a 0.20 slope gap needs about 20,000 per
group, and at SD 1.00 the same gap needs about 29,000. A flatter risk
distribution carries less information about the slope. Aim 3 reports the observed
SD of the linear predictor alongside the power calculation for this reason.
"""

import statistics

import numpy as np

SEED = 20260806
N_REPLICATES = 400
EVENT_RATE = 0.015
LP_SD = 1.25
LP_SD_ALTERNATIVE = 1.00   # reported as a sensitivity, see check_sd_sensitivity
ALPHA = 0.05
POWER = 0.80

# subgroup sizes the appendix quotes, with the standard errors it claims
SUBGROUP_SIZES = (5_000, 20_000, 100_000)
CLAIMED_SE = {5_000: 0.101, 20_000: 0.050, 100_000: 0.022}
CLAIMED_MDE = {5_000: 0.40, 20_000: 0.20, 100_000: 0.09}

REFERRAL_THRESHOLD = 0.03
YARDSTICK_SLOPE = 0.80


def _logit(p):
    return np.log(p / (1 - p))


def _expit(x):
    return 1.0 / (1.0 + np.exp(-x))


def _intercept_for_event_rate(target_rate, lp_sd, n_quad=20001):
    """Mean of the linear predictor giving a marginal event rate of target_rate.

    The linear predictor is N(mu, lp_sd^2); the marginal risk is the integral of
    expit over that distribution. Solved by bisection on mu.
    """
    grid = np.linspace(-8.0, 8.0, n_quad)
    weights = np.exp(-0.5 * (grid / 1.0) ** 2)
    weights /= weights.sum()

    def marginal(mu):
        return float((_expit(mu + lp_sd * grid) * weights).sum())

    low, high = -12.0, 4.0
    for _ in range(200):
        mid = 0.5 * (low + high)
        if marginal(mid) < target_rate:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def fit_recalibration(linear_predictor, outcome, max_iter=50, tol=1e-10):
    """Logistic regression of outcome on the model's own linear predictor.

    Returns (intercept, slope). The slope is the calibration slope: 1.0 means the
    predictions are correctly scaled, below 1.0 means they are too extreme.

    Newton-Raphson by hand, so the estimator is visible rather than delegated.
    """
    design = np.column_stack([np.ones_like(linear_predictor), linear_predictor])
    beta = np.zeros(2)
    for _ in range(max_iter):
        eta = design @ beta
        mu = _expit(eta)
        weight = np.clip(mu * (1 - mu), 1e-12, None)
        gradient = design.T @ (outcome - mu)
        hessian = design.T @ (design * weight[:, None])
        step = np.linalg.solve(hessian, gradient)
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    return beta[0], beta[1]


def analytic_slope_se(linear_predictor):
    """Fisher-information standard error for the slope, intercept profiled out.

    For a two-column design the inverse information gives
        var(slope) = 1 / sum(w_i * (lp_i - lp_bar_w)^2)
    where w_i = p_i(1 - p_i) and lp_bar_w is the w-weighted mean of lp. This is
    the closed form the simulation is checked against.
    """
    p = _expit(linear_predictor)
    w = p * (1 - p)
    centred = linear_predictor - (w * linear_predictor).sum() / w.sum()
    return float(1.0 / np.sqrt((w * centred**2).sum()))


def check_precision():
    """Step 1: the standard error of the calibration slope at each subgroup size."""
    print("1. PRECISION OF THE CALIBRATION SLOPE")
    mu = _intercept_for_event_rate(EVENT_RATE, LP_SD)
    print(f"   linear predictor ~ N({mu:.3f}, {LP_SD}^2) gives a "
          f"{EVENT_RATE:.1%} event rate")
    print(f"   {N_REPLICATES} replicates per size, slope fitted by "
          f"Newton-Raphson\n")
    print(f"   {'n':>8} {'events':>8} {'SE (sim)':>10} {'SE (Fisher)':>12} "
          f"{'claimed':>8}")

    standard_errors = {}
    rng = np.random.default_rng(SEED)
    for n in SUBGROUP_SIZES:
        slopes = []
        analytic = []
        for _ in range(N_REPLICATES):
            lp = mu + LP_SD * rng.standard_normal(n)
            y = (rng.random(n) < _expit(lp)).astype(float)
            if y.sum() < 5:
                continue
            slopes.append(fit_recalibration(lp, y)[1])
            analytic.append(analytic_slope_se(lp))
        simulated = statistics.stdev(slopes)
        fisher = float(np.mean(analytic))
        standard_errors[n] = simulated
        print(f"   {n:>8,} {n * EVENT_RATE:>8.0f} {simulated:>10.3f} "
              f"{fisher:>12.3f} {CLAIMED_SE[n]:>8.3f}")

        assert abs(simulated - fisher) / fisher < 0.10, (
            f"simulation and closed form disagree at n={n}: "
            f"{simulated:.4f} vs {fisher:.4f}")
        assert abs(simulated - CLAIMED_SE[n]) < 0.006, (
            f"appendix Z claims SE {CLAIMED_SE[n]} at n={n}, "
            f"simulation gives {simulated:.4f}")

    print("\n   Simulation and closed form agree, so the standard errors are a")
    print("   property of the design and not an artefact of one random draw.")
    return standard_errors


def check_detectable_difference(standard_errors):
    """Step 2: smallest between-group slope difference detectable at 80% power."""
    print("\n2. WHAT THAT DETECTS, BETWEEN TWO SUBGROUPS")
    z_alpha = statistics.NormalDist().inv_cdf(1 - ALPHA / 2)
    z_power = statistics.NormalDist().inv_cdf(POWER)
    multiplier = (z_alpha + z_power) * np.sqrt(2.0)
    print(f"   two-sided alpha {ALPHA}, power {POWER:.0%}: "
          f"MDE = {z_alpha:.3f} + {z_power:.3f} times SE times sqrt(2)")
    print(f"        = {multiplier:.3f} * SE, the sqrt(2) because two "
          f"independent groups\n")
    print(f"   {'n per group':>12} {'SE':>8} {'MDE (slope units)':>20} "
          f"{'claimed':>8}")

    detectable = {}
    for n in SUBGROUP_SIZES:
        mde = multiplier * standard_errors[n]
        detectable[n] = mde
        print(f"   {n:>12,} {standard_errors[n]:>8.3f} {mde:>20.3f} "
              f"{CLAIMED_MDE[n]:>8.2f}")
        assert abs(mde - CLAIMED_MDE[n]) < 0.02, (
            f"appendix Z claims MDE {CLAIMED_MDE[n]} at n={n}, "
            f"computed {mde:.3f}")

    print("\n   So at 20,000 per group a slope difference of 0.20 is detectable,")
    print("   and a null result bounds any real difference below roughly that.")
    print("   Nothing here licenses a claim about subgroups smaller than this.")
    return detectable


def check_sd_sensitivity():
    """The power calculation depends on the linear predictor SD, so state it."""
    print("\n3. HOW MUCH THE ANSWER DEPENDS ON THE DECLARED SD")
    z_alpha = statistics.NormalDist().inv_cdf(1 - ALPHA / 2)
    z_power = statistics.NormalDist().inv_cdf(POWER)
    multiplier = (z_alpha + z_power) * np.sqrt(2.0)
    target_gap = 0.20

    print(f"   n per group needed to detect a {target_gap} slope gap:\n")
    print(f"   {'LP SD':>7} {'info per obs':>14} {'n per group':>13}")
    required = {}
    for sd in (LP_SD, LP_SD_ALTERNATIVE):
        mu = _intercept_for_event_rate(EVENT_RATE, sd)
        grid = np.linspace(-8.0, 8.0, 200001)
        density = np.exp(-0.5 * grid**2)
        density /= density.sum()
        lp = mu + sd * grid
        p = _expit(lp)
        w = p * (1 - p)
        e_w = (w * density).sum()
        e_wlp = (w * lp * density).sum()
        e_wlp2 = (w * lp**2 * density).sum()
        info = e_wlp2 - e_wlp**2 / e_w
        n_needed = (multiplier / target_gap) ** 2 / info
        required[sd] = n_needed
        print(f"   {sd:>7.2f} {info:>14.6f} {n_needed:>13,.0f}")

    assert 19_000 < required[LP_SD] < 21_000, (
        f"at SD {LP_SD} expected about 20,000 per group, "
        f"got {required[LP_SD]:,.0f}")
    assert required[LP_SD_ALTERNATIVE] > required[LP_SD] * 1.3, (
        "a flatter linear predictor should demand a materially larger sample")

    print(f"\n   A flatter risk distribution carries less information about the")
    print(f"   slope, so the same gap costs about "
          f"{required[LP_SD_ALTERNATIVE] / required[LP_SD]:.1f} times the sample.")
    print(f"   The SD is therefore a declared assumption, not a detail: Aim 3")
    print(f"   reports the observed SD of the linear predictor next to the power.")
    return required


def check_clinical_yardstick():
    """Step 4: what a slope of 0.80 does to a patient the model puts at 3%."""
    print("\n4. WHY 0.20 IS THE YARDSTICK, NOT AN ARBITRARY NUMBER")
    lp_at_threshold = _logit(REFERRAL_THRESHOLD)
    mean_lp = _logit(EVENT_RATE)

    # The recalibration model is centred on the mean linear predictor, which is
    # what a fitted slope-and-intercept pair means in practice: the slope rotates
    # predictions about the cohort average, it does not pivot them about zero.
    centred = _expit(mean_lp + YARDSTICK_SLOPE * (lp_at_threshold - mean_lp))
    # The same slope applied with the intercept pinned at zero is a different
    # model and gives a different, larger number. It is shown so that the
    # distinction is on the record rather than discovered in a viva.
    uncentred = _expit(YARDSTICK_SLOPE * lp_at_threshold)

    print(f"   a patient the model places at {REFERRAL_THRESHOLD:.0%}, "
          f"cohort mean risk {EVENT_RATE:.1%}")
    print(f"   slope 1.00 (calibrated):              "
          f"{REFERRAL_THRESHOLD * 100:.2f}%")
    print(f"   slope {YARDSTICK_SLOPE:.2f}, recalibration centred on the "
          f"cohort mean: {centred * 100:.2f}%")
    print(f"   slope {YARDSTICK_SLOPE:.2f}, intercept pinned at zero "
          f"(NOT the estimand): {uncentred * 100:.2f}%")

    assert abs(centred - 0.0261) < 0.0005, (
        f"appendix Z quotes 2.61%, computed {centred * 100:.3f}%")
    assert uncentred > 0.05, "the uncentred version should be visibly different"

    print(f"\n   The centred figure is the estimand: the recalibration model")
    print(f"   regresses the outcome on the linear predictor, so the slope")
    print(f"   rotates predictions about the cohort average.")
    print(f"   {REFERRAL_THRESHOLD * 100:.1f}% against {centred * 100:.2f}% is "
          f"the difference between referring and not referring,")
    print(f"   which is why a slope gap of 0.20 is the smallest gap worth")
    print(f"   calling a finding, and why it is also what 20,000 per group buys.")
    return centred


def check_null_is_a_bound(detectable):
    """The reporting rule that follows from the steps above."""
    print("\n5. WHAT THE REPORTING RULE HAS TO BE")
    powered = [n for n in SUBGROUP_SIZES if detectable[n] <= 0.20 + 1e-9]
    print(f"   sizes at which a 0.20 gap is detectable: "
          f"{', '.join(f'{n:,}' for n in powered)}")
    print(f"   smallest such size: {min(powered):,} per group "
          f"({min(powered) * EVENT_RATE:.0f} events)")
    assert min(powered) == 20_000, "expected 20,000 to be the powered size"
    print("   Below that, a null bounds nothing worth reporting, so subgroup")
    print("   sizes are prespecified and reported BEFORE any comparison and an")
    print("   underpowered null is never written up as a null.")
    print("   Finding no difference at 20,000 per group is itself the")
    print("   contribution: no one has measured that gap for symptom codes.")


if __name__ == "__main__":
    standard_errors = check_precision()
    detectable = check_detectable_difference(standard_errors)
    check_sd_sensitivity()
    check_clinical_yardstick()
    check_null_is_a_bound(detectable)
    print("\nAll assertions passed. Every quantity is a design property under")
    print("declared assumptions; none is an estimate from patient data.")
