# Differential symptom recording and the calibration of cancer risk prediction models in UK primary care

Figures, extracted source values, and verification code for a five-minute DPhil
interview presentation (Nuffield Department of Primary Care Health Sciences,
University of Oxford).

Everything plotted in the presentation can be regenerated from this repository.
No value is typed into the plotting code: each figure reads from a CSV in
[`data/`](data), and each CSV records which published table it was transcribed
from.

---

## Why this repository exists

The presentation makes three empirical claims. Each is checkable here:

| Claim | Evidence | Where |
|---|---|---|
| Symptom coding differs by patient ethnicity and deprivation, in the database the project will use | d'Elia et al. (2025), CPRD Aurum, 70,115 adults, Table 3 | [`figures/fig_forest.png`](figures/fig_forest.png) |
| Recording completeness alone moves a predictive value across the 3% referral threshold | Price et al. (2016), CPRD, 20,958 symptom records | [`figures/fig_threshold.png`](figures/fig_threshold.png) |
| The symptoms nearest that threshold are the non-specific ones | Barclay et al. (2024), CPRD, 1,622,419 patients, Figure 2 | [`figures/fig_barclay_panels.png`](figures/fig_barclay_panels.png) |

The two symptoms that appear in **both** d'Elia's coding-disparity list and
Barclay's fifteen cancer presenting symptoms are **fatigue** and **weight loss**.
These are the non-specific ones, and they are also the ones whose single-site
cancer risk never crosses 3%. That overlap is what closes the argument, and it is
computed in `verification/check_sources.py`, not asserted.

---

## Reproducing

```bash
pip install matplotlib numpy pillow pypdfium2

# three figures, from the CSVs in data/
python code/figures.py

# consistency and threshold checks on every published value used
python verification/check_sources.py

# panels extracted from Barclay Fig 2 (download the open-access PDF first)
python code/barclay_panels.py /path/to/bmjonc-2024-000500.pdf
```

`code/figures.py` fails with an `AssertionError` if any two text labels collide,
if a label runs off the canvas, or if the dashed threshold segment passes through
any text. A figure that silently overlaps its own annotations cannot be
committed.

Those checks exist because earlier drafts of these figures failed in ways a
text-vs-text test could not see:

| Fault | Why the check missed it | Fix |
|---|---|---|
| A dashed threshold line ran through four labels | a `Line2D` is not a `Text` | explicit line-vs-text geometry (`check_panels_clean`) |
| A leader line terminated on the wrong row | geometry was consistent, just wrong | anchor annotations to looked-up row coordinates |
| An arrow pointed right while its label read "never crosses" | no check on semantic direction | the schematic was replaced by the source figure |

The last one is the reason `fig_barclay_panels.png` reuses the published figure
instead of redrawing it. Where a source figure exists and its licence permits
reuse, reusing it removes a whole class of error.

---

## What the verification script found

`verification/check_sources.py` tests each published confidence interval for
internal consistency: a 99.76% CI on an odds ratio should be approximately
symmetric on the log scale, so `sqrt(lower × upper)` should sit close to the
point estimate.

All fifteen associations in **Table 3** of d'Elia et al. (2025) pass. Three
values as printed in the **abstract** of the same paper do not:

| As printed in the abstract | Problem | Table 3 gives |
|---|---|---|
| Muscle cramps 1.71 (1.44–2.57) | geometric centre is 1.92, not 1.71 | 1.71 (1.14–2.57) |
| Falls, most deprived 1.37 (2.03–1.82) | lower limit exceeds the upper limit | 1.37 (1.03–1.82) |
| Shoulder pain 1.44 (1.25–1.66), South Asian | Table 3 attributes this estimate to Black ethnicity | South Asian is 1.33 (1.07–1.65) |

**The presentation uses the Table 3 values throughout.** This is recorded
because verifying that an interval is consistent with its own point estimate is
the kind of check the symptom code lists in Aim 1 will require. It does not affect
the argument, which holds on either set of numbers.

---

## A note on the Barclay figure

Barclay et al. report cumulative cancer risk **stacked by cancer site**, and
their claim concerns risk at *individual* sites. The stacked total is therefore
not the quantity of interest, and a 3% line drawn across a stacked axis invites
a misreading.

`code/barclay_panels.py` measures each site's band separately, by matching pixel
colours against the paper's own legend swatches:

- **Haematuria**: urological cancer alone reaches **11.6%** at the oldest age, and two sites exceed 3%.
- **Fatigue**: the largest single site reaches **1.8%**. No site exceeds 3%, although the stack sums to 6.1%.

The annotation on the figure reports the single-site value, and the script
asserts this contrast holds before drawing anything.

The 3% reference line is not placed by eye. The script locates the panel's own
gridlines and black baseline, fits pixel row against axis value, and raises an
`AssertionError` if the fit residual exceeds 3 pixels.

---

## Why these four symptom/cancer pairs

All four pairs that Price et al. report are shown. The selection is theirs, not
mine: they chose three clinical features as representative of alarm features and
of "low-risk but not no-risk" features, namely visible haematuria for bladder
cancer, jaundice for pancreatic cancer, and abdominal pain as a low-risk symptom
common to both cancers. Abdominal pain therefore yields two pairs, one per
cancer, giving four in total.

Only haematuria crosses the 3% threshold when free-text records are recovered,
and one detail matters for how strongly that can be claimed: the corrected
estimate is 2.9% with a 95% CI of 2.6 to 3.2, which still includes 3. The figure
shows the span of both intervals as a pale band, so this is visible instead of
buried in the caption. The claim the project rests on is that recording
completeness moves a predictive value across a decision threshold, not that this
particular point estimate is definitively below it.

---

## How this sits against existing work in the department

The Nuffield Department of Primary Care Health Sciences has already asked a
closely related question of a different predictor class:

- **Virdee, Oke, Bankhead, et al. (2025)**, *JMIR Cancer* 11:e70275. Reviewed 7 cancer risk models built on blood-test trends. Across 14 external validations, **one** assessed calibration, and all but one study scored high risk of bias in the analysis domain.
- **Nicholson, Bankhead, et al. (2019)**, *BMC Medicine* 17:222. Characterised weight recording in 5 million adults, and found the **white** ethnic group *less* likely to have weight recorded.

Two consequences for this proposal:

1. The Year 1 scoping review is not a duplicate. Theirs covers blood-test trends, where the predictor is ordered and returned as a number; this one covers symptom codes, where the predictor exists only if someone entered it. A Europe PMC search for reviews combining symptom-recording completeness with prediction-model calibration in primary care returns nothing.
2. The direction of any disparity is an open question, not an assumption. Weight recording and symptom coding run in *opposite* directions by ethnicity, plausibly because one is driven by chronic-disease review and the other by what is raised in a consultation. The fairness argument does not need a direction: if recording differs between groups, model trustworthiness differs between groups.

---

## Repository contents

```
code/figures.py                       three figures, plotted from data/*.csv
code/barclay_panels.py                panel extraction + per-site band measurement
data/delia2025_table3.csv             15 associations, transcribed from Table 3
data/price2016_ppv.csv                predictive values before/after text recovery
data/calibration_slopes.csv           subgroup and pooled calibration slopes
data/barclay2024_measured_bands.csv   per-site bands, written by barclay_panels.py
verification/check_sources.py         consistency, threshold, bridge and power checks
figures/                              generated output
```

---

## Sources

- **d'Elia A, Baranskaya A, Haroon S, et al.** Prodromal symptoms of rheumatoid arthritis in a primary care database: variation by ethnicity and socioeconomic status. *Rheumatology* 2025;64(3):1029–1035. [doi:10.1093/rheumatology/keae157](https://doi.org/10.1093/rheumatology/keae157). CC BY 4.0
- **Price SJ, Stapley SA, Shephard E, et al.** Is omission of free text records a possible source of data loss and bias in CPRD studies? *BMJ Open* 2016;6:e011664. [doi:10.1136/bmjopen-2016-011664](https://doi.org/10.1136/bmjopen-2016-011664)
- **Virdee PS, Collins KK, Smith CF, et al.** Clinical prediction models incorporating blood test trend for cancer detection: systematic review, meta-analysis, and critical appraisal. *JMIR Cancer* 2025;11:e70275. [doi:10.2196/70275](https://doi.org/10.2196/70275)
- **Nicholson BD, Aveyard P, Bankhead C, et al.** Determinants and extent of weight recording in UK primary care. *BMC Medicine* 2019;17:222. [doi:10.1186/s12916-019-1446-y](https://doi.org/10.1186/s12916-019-1446-y)
- **Barclay ME, Renzi C, Harrison H, et al.** Cancer incidence and competing mortality risk following 15 presenting symptoms in primary care. *BMJ Oncology* 2024;3:e000500. [doi:10.1136/bmjonc-2024-000500](https://doi.org/10.1136/bmjonc-2024-000500). CC BY 4.0
- **Luijken K, Groenwold RHH, Van Calster B, et al.** Impact of predictor measurement heterogeneity across settings on the performance of prediction models. *Statistics in Medicine* 2019;38(18):3444–3459. [doi:10.1002/sim.8183](https://doi.org/10.1002/sim.8183)
- **Van Calster B, McLernon DJ, van Smeden M, et al.** Calibration: the Achilles heel of predictive analytics. *BMC Medicine* 2019;17:230. [doi:10.1186/s12916-019-1466-7](https://doi.org/10.1186/s12916-019-1466-7)

### Figure reuse

`figures/fig_barclay_panels.png` reproduces part of Figure 2 of Barclay et al.
(2024), *BMJ Oncology* 3:e000500, under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Modifications: two of
fifteen panels shown; a 3% reference line added; two text annotations added
reporting values measured from the figure. The underlying plotted data are
unaltered. The source PDF is not redistributed here.

---

## Licence

Code and generated figures: [MIT](LICENSE), except
`figures/fig_barclay_panels.png`, which is a derivative of a CC BY 4.0 figure and
carries that licence with attribution as above.

Transcribed values in `data/` are facts from the published papers cited, held
here for verification.
