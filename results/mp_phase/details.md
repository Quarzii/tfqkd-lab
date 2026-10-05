# Mode-pairing phase block — executed evidence

Date: 2026-10-03. Scope: phase increments and published raw X-error comparisons only.
No mode-pairing key-rate expression, pairing optimizer, or frequency estimator is implemented.
The existing TF-QKD modules, web files and calculation configurations remain byte-identical
across the 40-file SHA-256 baseline (`before_hashes.json`). Existing uncommitted web changes
were preserved. This task adds an independent module rather than changing TF physics.

## Reproduction

From the repository root:

```bash
python -m tfqkd.mp_phase examples/mode_pairing/phase.toml --output results/mp_phase/example.json
python tools/digitize_mp_curves.py sources/data/mp_phase/digitization_config.json --output sources/data/mp_phase/digitized
python -m scripts.validate.validate_mp_phase
python -m unittest discover -s tests
TFQKD_FULL_TESTS=1 python -m unittest discover -s tests
```

`details.json` includes executed unit-test, CLI, tracking sensitivity and TF regression records.
The phase validator writes `analytic.json` **before** computing experimental comparisons.
Re-running that validator regenerates the physical comparison evidence; external test logs
and this narrative remain separate files. No numerical coefficient is fitted in this task.
The command-line demonstration uses explicit research parameters, not apparatus defaults.

## Implemented relations and their provenance

| Relation | Source / status |
|---|---|
| Free laser r3/f³ + r2/f² [fc/(f+fc)]² | Bertaina F1; separate summands keep separate lower limits |
| Stabilized laser, residual controller and cavity terms | Bertaina F2–F4; existing TF controller reused without modification |
| Single-pass fiber lL/f² [fc1/(f+fc1)]² | Bertaina Eq.6; existing TF spectrum reused |
| Independent sum of both lasers and both arms | Bertaina Eq.7; one-pass geometry in Zeng Fig.2 and Box1 |
| Increment filter 4 sin²(πfΔt) | Derived phase-increment variance from Di Domenico Eq.1 (one-sided convention) |
| White frequency level h0=linewidth/π and S_phase=h0/f² | Di Domenico white-frequency paragraph after Eq.2; Zhang Eq.3/Appendix B Lorentzian model |
| Raw X error 1/2 − M sin(2π/M)/(8π) exp(−D/2) | Zhang Eq.B10; Gaussian increment model, weak symmetric channels, perfect-detector approximation in Appendix B |
| Raw X error 1/2 − V2/2 exp(−D/2) cos(2πΔfΔt) | Zhou Eq.2; Gaussian stochastic drift plus explicitly supplied residual beat |
| Count-weighted discrete average | Zhang Eq.D5 and following discrete-sum recommendation; actual Zhu TableIV counts |
| Fiber short-interval expansion and equivalent drift | Derived below from Bertaina Eq.6 and Di Domenico Eq.1; not asserted to be a separately numbered published formula |

Primary sources: [Bertaina](https://arxiv.org/abs/2310.08621),
[Di Domenico](https://doi.org/10.1364/AO.49.004801),
[Zeng](https://doi.org/10.1038/s41467-022-31534-7),
[Zhou async](https://arxiv.org/abs/2212.14190v4),
[Zhang PRX](https://doi.org/10.1103/PhysRevX.15.021037),
[Zhu](https://arxiv.org/abs/2208.05649v2).
PDFs are in `sources/papers/`; extracted texts in `sources/data/mode_pairing/`.

### Geometry and normalization

Alice and Bob send their pulses once to Charlie (Zeng Fig.2/Box1). For independent
arm disturbances, the differential phase PSD is the sum of their one-pass PSDs;
the sign of the arm phase does not change its PSD. The common-laser round-trip K
of Bertaina Eq.5 is therefore absent. Two temporal samples introduce the distinct
filter |exp(i2πfΔt)−1|²=4 sin²(πfΔt). This temporal factor is retained exactly once.
The two lasers are explicit independent inputs; their PSDs need not be identical.
Lengths in drift tables are L_A+L_B, not the length of each arm.

All PSDs supplied to this module are one-sided phase spectra in rad²/Hz.
Di Domenico Eq.1 writes the field correlation as exp[−2∫S_frequency sin²(πfΔt)/f² df].
For its Gaussian phase model this exponent is −D/2, giving the implemented
increment integral and the frequency-to-phase conversion S_phase=S_frequency/f².
A PSD alone does not prove Gaussian statistics; the X-error conversion remains
conditional on the Gaussian source model used by Zhou Eq.2 and Zhang Appendix B.
The 25% floor concerns raw coherent-state X errors, not a single-photon privacy error.
For M=32, Zhang B10 gives 0.25160328721394853 at D=0; Zhou's V2=0.46 gives 0.27.

### Component-specific lower limits and tracking

For S∝f^(−p), the increment integrand scales as f^(2−p) near zero. It is integrable
only when p<3. Consequently free-laser r3 and nonzero cavity C3/C4 require an
explicit T_track; white-frequency, free-laser r2, cavity C2, servo residual and
Eq.6 fiber start at zero. The controller residual's zero-frequency behavior follows
Bertaina F4. A missing T_track with any divergent component produces an explanatory error.
No apparatus tracking model is inferred from the time window.

Approved convention: only divergent components have f_min=1/T_track.
Each output lists cut_components and every component's actual lower bound.
The same inputs are computed at T_track/10, T_track and 10*T_track.
The numerical positive log-grid floor is not a physical low-frequency cutoff:
convergent terms include the zero endpoint through the analytic limit of f²S.
A divergent cutoff below that floor is covered by a newly generated log grid,
rather than a single wide trapezoid; the regression test covers this case.

## Analytical derivation, checked before experimental comparisons

Let t=|Δt|, A=l(L_A+L_B), c_f=fc1, x=f/c_f, q=2πc_f t. From Eq.6 and the increment filter:

\[
D_F=\frac{2A}{c_f}J(q),\qquad
J(q)=\int_0^\infty\frac{1-\cos(qx)}{x^2(1+x)^2}\,dx.
\]

The two derivatives used here are integrable:

\[
J''(q)=\int_0^\infty\frac{\cos(qx)}{(1+x)^2}\,dx
=1-q\int_0^\infty\frac{\sin(qx)}{1+x}\,dx.
\]

Setting y=qx in the last integral gives ∫sin(y)/(q+y)dy. Its q→0+ limit is π/2
(the Dirichlet integral); separating 0<y<1 from the tail bounds the correction by
O(q|log q|). With J(0)=J'(0)=0, integration twice gives:

\[
J(q)=\frac{q^2}{2}-\frac{\pi q^3}{12}+O(q^4|\log q|),
\]

and therefore

\[
D_F=4\pi^2 A c_f t^2\left[1-\frac{\pi^2}{3}c_f t
+O((c_f t)^2|\log(c_f t)|)\right],\qquad
\sigma_L=2\pi\sqrt{A c_f}.
\]

The next term is negative and cubic in |Δt|. Expanding sin² to fourth order
inside the infinite integral would be invalid: Eq.6 has an f^(−4) tail,
so the fourth spectral moment diverges. The expansion requires fc1|Δt|≪1,
and its coefficients describe the zero-lower-limit, infinite-upper-limit model.
They are not a long-interval drift formula or an extra fitted noise contribution.

Independent adaptive quadrature uses R(a)=∫sinc(ax)²/(1+x)² dx, a=fc1 t;
0..1 directly, the middle interval with log(x), and the tail with z=ax.
The reference tail ends at z=200; its omitted ratio is bounded by
 a/(3π²·200³). This bound, rather than a fitted uncertainty, is recorded in analytic.json.
The implemented finite-band trapezoid is compared against that independent reference.

### Analytic results

| fc1 Δt | Relative leading overestimate | Relative cubic error | Observed (1−R)/a |
|---:|---:|---:|---:|
| 1e-06 | 3.28979355e-06 | -8.54063487e-11 | 3.289782728 |
| 3e-06 | 9.8689979e-06 | -7.03905489e-10 | 3.289633501 |
| 1e-05 | 3.28927352e-05 | -7.02822556e-09 | 3.289165334 |
| 3e-05 | 9.86490286e-05 | -5.67516801e-08 | 3.287976598 |
| 0.0001 | 0.000328543411 | -5.51488901e-07 | 3.284355056 |
| 0.0003 | 0.00098361544 | -4.31578981e-06 | 3.275496304 |
| 0.001 | 0.00326046641 | -4.01282258e-05 | 3.249870319 |
| 0.003 | 0.00966689149 | -0.000298121305 | 3.191445803 |
| 0.01 | 0.0313501123 | -0.00257994641 | 3.039715796 |

The limiting cubic coefficient is π²/3=3.2898681337.
Maximum fiber quadrature discrepancy: 5.26913974674e-08; tolerance 2e-06.

### White-frequency normalization and mesh convergence

For one laser, S_frequency=h0, linewidth=πh0 (Di Domenico), the exact finite-band result is

\[
D=4h_0\left[\pi t\,\mathrm{Si}(2\pi f_{max}t)
-\frac{\sin^2(\pi f_{max}t)}{f_{max}}\right].
\]

For f_max→∞ this tends to 2π·linewidth·t per laser, hence
4π·linewidth·t for two identical independent lasers, matching Zhang Eq.3.
This explicitly checks the one-sided normalization rather than hiding a factor of two.
The initial unit assertion used the infinite-band law against a finite-band integral:
at t=10 µs the discrepancy was approximately 1.0411e−4 versus tolerance 1e−4.
The finite upper-band deficit alone is about 1.0132e−4. The test reference was
corrected to Di Domenico's exact finite-band expression; no coefficient or PSD was tuned.

Chosen validation log grid: positive floor 1e−8 Hz, upper endpoint 1e8 Hz,
65537 points, plus zero for convergent terms. These are numerical decisions.
They do not modify the existing TF fast/reference grids. Convergence is checked
against independent references rather than against an expected experimental curve.

| Points | Fiber maximum relative error vs adaptive quadrature | White maximum relative error vs exact finite band |
|---:|---:|---:|
| 16385 | 8.42735913942e-07 | 2.79061643138e-06 |
| 65537 | 5.26913974674e-08 | 1.22165613403e-06 |
| 131073 | 1.31891757515e-08 | 1.6189097507e-07 |

The 16385-point numerical test grid does **not** pass the stricter 2e−6 white reference gate (2.79062e−6). The 65537-point validation grid passes both gates. Curve refinement against 131073 points is independently gated at 1e−4.

## Explicit T_track example (not an apparatus prediction)

Free lasers and fiber coefficients are Table III Bertaina; arm lengths 100.93 km each are Zhou Table S1 geometry only. T_track=0.01 s is a demonstration input. Only laser_A.r3 and laser_B.r3 are cut. Both r2 and both fiber terms retain f_min=0.

| Δt [s] | D at T_track/10 [rad²] | D nominal [rad²] | D at 10 T_track [rad²] | Short / long change [%] |
|---:|---:|---:|---:|---:|
| 1e-06 | 0.0101537284955 | 0.0106991426337 | 0.0112445571538 | -5.0977369 / 5.0977404 |
| 1e-05 | 0.20344808321 | 0.257985678462 | 0.312527092283 | -21.139776 / 21.141256 |
| 0.0001 | 4.84094802585 | 10.2567744092 | 15.7105339344 | -52.802432 / 53.172267 |

## Drift-rate comparison without fitting

l=44 rad² Hz/km and fc1=100 Hz are from Bertaina Table III, a different line. This is the requested order comparison; no apparatus transfer accuracy is assumed. Zhou Fig.3 publishes empirical drift rates; Zhang Eq.3/Fig.1 adopts 4000 rad/s in its source model.

| Source | Total L [km] | Published σ_L [rad/s] | Predicted [rad/s] | Predicted / published | Difference [%] |
|---|---:|---:|---:|---:|---:|
| zhou2023_async | 201.86 | 2100 | 5921.494648 | 2.81975936 | 181.97594 |
| zhou2023_async | 306.31 | 3400 | 7294.353297 | 2.14539803 | 114.5398 |
| zhou2023_async | 413.73 | 5300 | 8477.439776 | 1.59951694 | 59.951694 |
| zhou2023_async | 508.16 | 5900 | 9395.208983 | 1.5924083 | 59.24083 |
| zhang2025 | 403 | 4000 | 8366.787322 | 2.09169683 | 109.16968 |

The Eq.6 model scales as sqrt(L); the exponent inferred from Zhou's **two endpoints only** is 1.1189237679, not 0.5. No coefficient was estimated from this diagnostic. Published σ_L/sqrt(L) values are 147.806713, 194.266685, 260.56577, 261.728956. The transferred coefficients overpredict every published drift value checked.

## Published X-error curves: digitization and conditional comparison

Tool: Poppler `pdftocairo -svg` plus the new `tools/digitize_mp_curves.py` using
standard-library XML and numpy. This is an explicitly adopted digitization method;
no curve fitting, smoothing, statistical covariance or hidden timestamps are recovered.
`digitization_config.json` records axis calibration in PDF points, color/path selection,
figure identities, exact selected path indices, sampling stride, and bin boundaries.
`digitized/digitization.json` records source/SVG/CSV hashes and conversion stderr/status.

- Zhou Fig.3(a,b), PDF p.4: eight experimental connected traces, every fifteenth
  vector vertex plus the last. The Exp.[Sim.] legend identifies bracketed symbols
  as simulation; the connected experimental traces are used for numerical comparison.
- Zhang Fig.1(b), PDF p.3: all selected blue post-tracking experimental marker centers.
  The left and right segments of the broken logarithmic axis have separate calibrations.
- Zhu Fig.6(a), Appendix B3, PDF p.12: six strong-reference bin curves, four length
  categories each. Their actual pulse-gap bounds are retained, not replaced by bin midpoints.

A Poppler conversion for Zhou returned −6 (assertion) after writing a complete,
parseable SVG. Its status is retained, not treated as a clean conversion. The plotted
coordinates/axis selection were visually checked against an independently rendered PDF
page; source-model residuals below give a separate numerical control. Zhang/Zhu conversions
returned zero. Digitization remains a graphical approximation, not raw experimental data.

| Curve | Selected points | Source |
|---|---:|---|
| zhou_3a_0 | 101 | Fig. 3(a), zhou2023_async |
| zhou_3a_1 | 101 | Fig. 3(a), zhou2023_async |
| zhou_3a_2 | 101 | Fig. 3(a), zhou2023_async |
| zhou_3a_3 | 101 | Fig. 3(a), zhou2023_async |
| zhou_3b_0 | 101 | Fig. 3(b), zhou2023_async |
| zhou_3b_1 | 101 | Fig. 3(b), zhou2023_async |
| zhou_3b_2 | 101 | Fig. 3(b), zhou2023_async |
| zhou_3b_3 | 101 | Fig. 3(b), zhou2023_async |
| zhang_1b_left | 9 | Fig. 1(b), zhang2025 |
| zhang_1b_right | 240 | Fig. 1(b), zhang2025 |
| zhu_6a_bin_0 | 4 | Fig. 6(a), Appendix B3, zhu2023 |
| zhu_6a_bin_1 | 4 | Fig. 6(a), Appendix B3, zhu2023 |
| zhu_6a_bin_2 | 4 | Fig. 6(a), Appendix B3, zhu2023 |
| zhu_6a_bin_3 | 4 | Fig. 6(a), Appendix B3, zhu2023 |
| zhu_6a_bin_4 | 4 | Fig. 6(a), Appendix B3, zhu2023 |
| zhu_6a_bin_5 | 4 | Fig. 6(a), Appendix B3, zhu2023 |

Total: 1081 points on 16 curves.

### Parameter provenance and unfilled inputs

Zhou: arm lengths from Table S1; linewidth 1 Hz per laser from the main setup;
V2=0.46, drift rates and nominal frequency offsets from Eq.2/Fig.3.
The green curve is published as a negligible offset, nominal 0.01 kHz / <10 Hz;
10 Hz is used as the paper's nominal/boundary comparison, not claimed as the exact
mean of that trace. Supplementary Fig.S5 gives mean and SD for the larger beats;
it does not establish their full noise PSD. The phase calculation models both lasers
as Lorentzian white-frequency sources, an **explicit approximation**, rather than
claiming that their linewidth specifies a measured colored spectrum.

Zhang: linewidth 100 Hz, M=32, 403-km geometry, source σ_L=4000 rad/s from
Eq.3/Appendix B/Fig.1/III. Equal arms split the symmetric total; identical fiber
coefficients make their sum the only length input to D. The source simplified
comparison omits residual frequency-estimator error. That idealization is stated,
not turned into a measured zero. The specific differential fiber PSD is unavailable.

Thus the following curves are **conditional source-model comparisons** under their
Gaussian/Lorentzian approximations. They are not a certified apparatus prediction or
an independent validation of a complete tracking model. The source formula controls
use the published source drift, separately from the transferred Eq.6 PSD. No physical
parameter is adjusted to minimize the resulting discrepancies.

| Curve | RMS source formula − graphic [probability] | RMS transferred PSD − graphic [probability] | Mean transferred residual | Maximum relative mesh change |
|---|---:|---:|---:|---:|
| zhou_3a_0 | 0.002330797183 | 0.07578353615 | 0.05822832308 | 5.501239642e-08 |
| zhou_3a_1 | 0.005263284272 | 0.05512564438 | -0.01572222967 | 5.501239286e-08 |
| zhou_3a_2 | 0.002578259666 | 0.05146754771 | -0.002934326365 | 5.501238509e-08 |
| zhou_3a_3 | 0.003064317695 | 0.0508226387 | 0.0002733542487 | 5.50123358e-08 |
| zhou_3b_0 | 0.002368255372 | 0.0756351879 | 0.05817208189 | 5.011498061e-08 |
| zhou_3b_1 | 0.001088886323 | 0.04757152666 | 0.02890092976 | 4.651483954e-08 |
| zhou_3b_2 | 0.00172626032 | 0.02591132294 | 0.01302137611 | 4.470109194e-08 |
| zhou_3b_3 | 0.001700614292 | 0.02454184871 | 0.01171941555 | 4.373802343e-08 |
| zhang_1b_left | 0.002787551491 | 0.002787620284 | -0.002352453712 | 5.060790698e-08 |
| zhang_1b_right | 0.009775607348 | 0.04515832581 | 0.02962670641 | 5.692424537e-07 |

RMS is an unweighted diagnostic over the selected published graphic points. It is
not an interval-distribution average and has no established statistical significance
without event covariance/confidence metadata. For length-dependent Zhou Fig.3(b),
transferred-noise RMS is 2.45–7.56 percentage points; on Zhang's right-axis segment
it is 4.52 percentage points. The transferred line systematically increases the
phase decoherence rate in those comparisons, consistent with its drift overprediction.
The oscillatory beat curves can have either sign of residual, so no all-frequency
sign is claimed. Plot: comparison (`x_error_comparison.png`; research archive not included in the public snapshot); all points in comparison CSVs.

### Zhu: digitized, prediction intentionally not calculated

Published 2-kHz laser linewidth and 625-MHz clock do not specify residual MLE
frequency-estimate noise, the distribution of pair intervals within the Fig.6 bins,
or a complete residual phase/selection model for those strong-reference bin averages.
All six curves therefore carry NO_PREDICTION with named missing inputs in JSON.
No estimator residual, midpoint, pairing distribution, or theoretical bin average is invented.
Plot: digitized curves (`zhu_digitized.png`; research archive not included in the public snapshot). This is a validation limit, not a failed
numeric implementation of a fully specified case.

The separately published four-bin counts for the 202-km run in TableIV allow an
exact **observed** discrete average: ∑errors / ∑pairs, equivalent to the count-weighted
sum recommended after Zhang D5. This does not predict Fig.6 or recover within-bin intervals.

| Data | Pairs | Error pairs | Count-weighted raw X error |
|---|---:|---:|---:|
| strong | 87829739 | 25281666 | 0.287848583952 |
| qkd | 506700 | 155979 | 0.3078330373 |

## Executed checks

| Check | Result | Evidence |
|---|---|---|
| New phase unit checks | 12 passed (zero interval, no K, required tracking band, cutoff isolation/full coverage, PSD decomposition, normalization, cubic limit, raw error floor, discrete counts, input rejection) | tests/test_mp_phase.py; full_tests.log |
| Quick unit group | 102 executed and passed, 6 full-only skipped; process 4.450 s | quick_tests.json/log |
| Full unit group | 108 executed and passed; process 9.130 s | full_tests.json/log |
| Independent fiber / exact white reference | PASS at configured 2e−6 gate | analytic.json |
| MP experimental-curve mesh refinement | PASS at 1e−4; max 5.692424537e-07 | details.json |
| Existing TF files | 40 byte-identical | before_hashes.json + details.json tf_unchanged |
| Existing TableI fast/reference | PASS; variance 5.100273135e-05 / 0.0001; time 1.483966661e-05 / 0.0001; key 1.15959765e-05 / 0.001 | tf_regression.json |
| Example CLI | exit0; process 0.516 s | example_run.json/example.json |

The analytic-and-curve validation process took 20.419 s, including independent
adaptive quadrature, 10 comparisons, mesh refinement and plot generation. This is not
a phase-only user-run timing. Existing reference T1–T7 acceptance remains archived;
this task reran the unit groups and TableI fast/reference regression, not the entire
historical plotting pipeline. No claim of a fresh full T1–T7 pipeline run is made.

## Decisions and remaining limits

- Approved band convention 1/T_track for divergent increment components only; not a tracking model.
- Numerical log grid, explicit removable zero endpoint, finite upper bound, interval batching and tolerances are calculation choices, not apparatus coefficients.
- Zhou linewidth is represented by the explicitly declared Lorentzian white-frequency approximation; short-term linewidth does not determine a complete PSD.
- Zhang comparison uses the source simplified model without estimator-residual noise; its absence is conditional idealization, not a measured zero.
- Vector digitization samples every fifteenth Zhou trace vertex and all selected Zhang/Zhu marker centers; RMS uses selected published coordinates, not event probabilities.
- Zhu bin averages are not assigned artificial midpoint intervals; predictions remain blocked.


Missing open-data requirements are recorded in docs/UNKNOWNS.md as
MP_TRACKING_RESPONSE, MP_TRANSFERRED_FIBER, MP_LASER_PSD and
MP_ZHU_TRACKING_AND_BINS, under research-stage assumptions. They do not alter TF
laboratory reports. MP_GRAPHIC_COMPARISON is a validation limit.
No MP key-rate computation or actual frequency-tracking transfer function is implemented.
