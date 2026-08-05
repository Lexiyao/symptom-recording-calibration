"""
Two arithmetic checks on the Aim 2 design.

Run:  python verification/check_identification.py

1. Practice fixed effects cannot be used. Coding propensity is defined at
   practice level, so it is a linear combination of the practice indicators and
   its coefficient is not identified. This is shown by constructing the design
   matrix and reporting its rank deficiency, not by assertion.

2. How much unmeasured case-mix would be needed to reproduce an observed
   effect. Standard omitted-variable arithmetic for a practice-level exposure:
   the bias equals beta_U * corr(U, P) * sd(U) / sd(P). Setting the bias equal
   to the observed effect gives the confounder strength required to explain it
   away. This is the number to report alongside the estimate.

Neither check uses patient data. Both are properties of the design.
"""

import numpy as np

N_PRACTICES = 40
PATIENTS_PER_PRACTICE = 200
TRUE_PRACTICE_EFFECT = 1.2
SEED = 20260805


def check_fixed_effects_collinearity():
    """A practice-level exposure is absorbed by practice fixed effects."""
    print("1. CAN PRACTICE FIXED EFFECTS BE USED?")
    rng = np.random.default_rng(SEED)
    n = N_PRACTICES * PATIENTS_PER_PRACTICE
    propensity = rng.uniform(0.3, 0.9, N_PRACTICES)   # one value per practice
    practice = np.repeat(np.arange(N_PRACTICES), PATIENTS_PER_PRACTICE)
    covariate = rng.normal(size=n)
    outcome = (TRUE_PRACTICE_EFFECT * propensity[practice]
               + 0.4 * covariate
               + rng.normal(scale=0.5, size=n))

    dummies = np.zeros((n, N_PRACTICES))
    dummies[np.arange(n), practice] = 1
    with_fe = np.column_stack([propensity[practice], covariate, dummies])
    rank = np.linalg.matrix_rank(with_fe)
    deficiency = with_fe.shape[1] - rank
    print(f"   with practice fixed effects: {with_fe.shape[1]} columns, "
          f"rank {rank}, deficiency {deficiency}")
    assert deficiency >= 1, \
        "expected the practice-level exposure to be collinear with the dummies"
    print("   The exposure coefficient is not identified. Fixed effects are out.")

    without_fe = np.column_stack([np.ones(n), propensity[practice], covariate])
    beta = np.linalg.lstsq(without_fe, outcome, rcond=None)[0]
    print(f"   without fixed effects: coefficient {beta[1]:.3f} "
          f"(data generated with {TRUE_PRACTICE_EFFECT})")
    assert abs(beta[1] - TRUE_PRACTICE_EFFECT) < 0.15
    print("   Recovered, which is why the model uses practice RANDOM effects")
    print("   with patient-level standardisation instead.")


def check_confounding_bound(observed_effect=0.30, sd_propensity=0.15):
    """Confounder strength required to reproduce an observed effect.

    bias = beta_U * corr(U, P) * sd(U) / sd(P)
    """
    print("\n2. HOW STRONG WOULD UNMEASURED CASE-MIX HAVE TO BE?")
    print(f"   effect to explain away: {observed_effect} calibration-slope units")
    print(f"   between-practice SD of coding propensity: {sd_propensity}")
    print(f"\n   {'corr(U,P)':>10} {'sd(U)':>7} {'required slope effect per SD of U':>34}")
    results = {}
    for correlation in (0.3, 0.5, 0.7):
        for sd_u in (0.5, 1.0):
            required = observed_effect / (correlation * sd_u / sd_propensity)
            results[(correlation, sd_u)] = required
            print(f"   {correlation:>10.1f} {sd_u:>7.1f} {required:>34.3f}")
    headline = results[(0.5, 1.0)]
    print(f"\n   Headline: a confounder correlated 0.5 with coding propensity and")
    print(f"   SD 1.0 must move the calibration slope by {headline:.2f} per SD to")
    print(f"   reproduce the finding. That magnitude is reportable, and it is what")
    print(f"   the negative-control analyses are designed to detect.")
    assert abs(headline - 0.09) < 0.005, "headline bound changed unexpectedly"
    return results


if __name__ == "__main__":
    check_fixed_effects_collinearity()
    check_confounding_bound()
    print("\nAll assertions passed. Neither check uses patient data.")
