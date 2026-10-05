# Validation and discrepancies

Numerical consistency, reproduction of a source model, and agreement with
experimental measurements are distinct checks. The records below describe
specific configurations, not universal forecast accuracy.

## TF-QKD

- An independent implementation of spectra, integration, window, and duty gives a maximum relative discrepancy of **1.978 × 10⁻⁷**. [Record](../results/stage3_d/independent.json)
- Table I fast/reference grids use **4,097 / 65,537** points. Maximum variance discrepancy is **5.101 × 10⁻⁵**, below the **10⁻⁴** criterion. This grid check applies to the tested configurations. [Record](../results/stage3_a/grid_regression.json)
- Held-out Zhou rates, after fitting one intrinsic-error parameter at one length, differ by **8–11%** at the other shorter length. Equivalent length differences are **3.98–4.23 km** under the stated loss interpolation. These compare the protocol/loss calculation with a fitted parameter; they do not validate a spectral phase forecast. [Rates](../results/out_of_sample/validation.json), [length conversion](../results/stage3_d/validation_distance.json)
- At **615.59 km**, the calibrated model yields zero against the published **0.32 bit/s**. No empirical correction is applied.
- The parametric fiber spectrum does not reproduce the full shape of the measured 114 km curve. Classical compensation omits phase-measurement detection noise.

For the checked Bertaina inputs, classical compensation can lose its calculated
benefit around **30–70 km per arm**, while dual-band stabilization approaches
the ideal fiber-noise ceiling. This is configuration-dependent.
[Calculation record](../results/release_v1_0/conclusion_checks.json)

## Phase-observation comparisons

Pittaluga/Zhou examples use explicit phase brackets rather than same-window
quantum-channel RMS measurements at every key-rate length. Pittaluga's
0.104 rad comes from 605 km, while the CAL comparison is at 368.702 km.
Zhou's 0.099 rad reference-channel quantity is not a quantum-channel RMS
measurement at each of the three key-rate lengths. The comparisons assume
nondecreasing residual RMS and, for Zhou, comparable reference/quantum RMS.
Their effective 500 MHz clocks and zero extra stabilization downtime are timing
idealizations. Missing Pittaluga CAL QBER and numerical error-correction
inefficiency leave its calibration conditional on ideal limits.

## Mode-pairing

The spectrum-to-drift short-interval limit is checked analytically and
numerically. Transferring another line's fiber coefficients overpredicts the
published spool drift by **1.6–2.8×**. The four-point length fit has exponent
**1.168** with a conditional residual-based 95% interval **[0.743, 1.592]**.
The independent-section exponent 0.5 lies outside that interval, but published
measurement uncertainties are unavailable. This does not identify a spatial
correlation mechanism or calibrate field-route coefficients.
[Phase record](../results/mp_phase/details.json),
[length fit](../results/mp_spatial_analysis/details.json)

Key comparisons use observed finite counts as asymptotic expectation proxies,
without finite-session confidence bounds. Unreported vacuum X error gains use
an explicit random-error assumption of 1/2. Zhang's 101 km table prints
X2mu,0 = 145543908 and X0,2mu = 13488265; the resulting 10.7904 ratio is retained
literally rather than corrected as a presumed typo. Zhang held-out comparisons are
conditional on measured gains at each length. A fitted error absorbs estimator,
timing, tracking, and finite-key differences; it is not a measurement of
intrinsic misalignment. One calibration gives held-out ratios **0.780 / 1.265**;
a different calibration gives **1.655** at 403 km.
[Key comparison record](../results/mp_key/details.json)

For Zhu at **202 km**, the asymptotic analytical-decoy result is **1569.78 bit/s**,
**24.53%** below the reported **2080 bit/s**. Separately, replaying the published
finite components gives **1826.80 bit/s**, a **12.17%** discrepancy. Even zero
error-correction expense reaches only **1892.94 bit/s** with the printed privacy
components. Rounding and exact sent totals do not close the discrepancy; the
source does not identify a unique erroneous table entry.

At Zhu's **304 and 407 km**, the analytical X/privacy bounds establish no
positive key. Zhu documents decoy estimation and Chernoff-Hoeffding bounds but
does not identify its actual solver. Its solver must not be labelled LP on that
basis. No replacement estimator or empirical correction is applied. Zhu phase
forecasts are absent because drift and tracking residuals are unreported.
[Diagnostic record](../results/mp_zhu_diagnostics_v1_1/details.json),
[checks](../results/mp_zhu_diagnostics_v1_1/verification.json)

## Reproduce numerical checks

From the repository root:

```bash
python -m unittest discover -s tests
TFQKD_FULL_TESTS=1 python -m unittest discover -s tests
python scripts/validate/validate_stage3_a.py --calibrate
python scripts/validate/validate_mp_phase.py
```

These commands run without article PDFs. Source-dependent reproduction checks
require the papers listed in [SOURCES.md](../SOURCES.md). Unit tests and grid
convergence do not establish experimental forecast accuracy.
