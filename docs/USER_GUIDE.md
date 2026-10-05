# TF-QKD Laboratory User Guide — v1.1

This tool estimates Twin-Field QKD key rate, phase-noise impact and operating
reach from inputs that a laboratory can state explicitly. It has no equipment
preset library. Published examples in `examples/` are reproducibility fixtures,
not defaults for another installation.

The calculation uses [Bertaina2024](https://doi.org/10.1002/qute.202400032), the
Di Domenico linewidth conversion, the Williams round-trip residual model, and
the author protocol code where applicable. The wrapper adds laboratory input
handling, measured-input modes, conditional requirements for unknown quantities,
ranked sensitivity, HTML/Markdown reporting, and a local one-page web console.

## Quick Start

Run commands from the repository root.

```bash
python -m tfqkd.lab_run --help
python -m tfqkd.lab_run examples/bertaina2024_table3.toml --output results/my_run
python -m tfqkd.lab_web --port 8765
```

The CLI writes `report.md`, `report.html`, `result.json`, `spectrum.png`, and
`cumulative_variance.png` for spectral runs. Measured-phase runs do not produce
spectrum plots because no spectral chain is evaluated.

The web command starts an offline local server at `http://127.0.0.1:8765/`.
It binds only to loopback. CSV uploads are stored in a temporary directory and
retained for at most one hour (up to eight reports) and deleted when the server exits. The web console uses `workers=1`. Local reports use one worker by default to reuse cached transfer functions; set `numerics.workers` explicitly for independent parallel input cases.

In restricted sandboxes, process creation can be blocked. If a CLI run fails
with `Operation not permitted` while spawning workers, set:

```toml
[numerics]
workers = 1
```


## What The Tool Calculates

For a spectral prediction, the tool calculates:

- one-sided phase-noise spectra for laser, fiber and detection terms;
- phase variance by integrating above `1/tau_Q`;
- `tau_Q`, duty cycle and phase error;
- protocol key rate for SNS-AOPP or CAL;
- four fiber-stabilization rates: free, perfect fiber cancellation, classical
  round-trip compensation, and dual-band stabilization;
- ranked sensitivity for approved factor-two changes;
- optional numerical reach if `[reach] enabled = true`.

For a measured-phase run, the spectral chain is bypassed. You supply residual
phase RMS, an operating window, and the phase-estimation overhead. The tool reports
`e_phi = sigma_phi^2 / 4` (Eq.1, Bertaina2024) and retains the author protocol
Gaussian phase-error function for the key calculation. It reports this
as an estimate from an observation, not as a forecast of another line length or
another compensation scheme.

## Required Input Style

Every apparatus value must come from one of these routes:

- direct model coefficients;
- passport or familiar measured quantities with a documented conversion;
- a measured CSV spectrum or actuator response with mandatory normalization
  metadata;
- a measured residual phase RMS and time window.

If a required geometric, protocol, cutoff, wavelength or actuator parameter is
missing, the resolver returns a missing-input list. It does not fill Bertaina
values. The only reference exception is `physics.n`: when omitted, the tool uses
`n = 1.4682`, the typical effective group index at 1550 nm from Corning SMF-28
Ultra PI-1424-AEN, July 2025, page 2. The report marks that reference as a
default and records the previous n-sensitivity result: up to 1.818719% variance
change over n=1.44..1.48 in the Table I check.

## How To Get Each Input

Use explicit units. Keep source comments in your TOML, especially for values
taken from publications or datasheets.

| Node | Accepted input | How to obtain it |
| --- | --- | --- |
| Laser | `r3`, `r2`, `fc_hz`; plus `C4`, `C3`, `C2`, `B_hz`, `gamma`, `delta` for cavity mode | Fit a free-running laser PSD with F1, or enter published/model coefficients. |
| Laser linewidth | `laser.lorentz_width_hz` | Use a Lorentzian linewidth in Hz; the tool converts it to `r2` with Di Domenico2010 Eqs. 1 and 5. It does not infer `r3` or `fc_hz`. |
| Laser CSV | `[laser.spectrum] mode="fit"` or `mode="direct"` | Use a measured frequency or phase PSD. Fit mode identifies F1 coefficients; direct mode interpolates the PSD on the grid. |
| Line geometry | `arm_a_km` and `arm_b_km`, or `length_km` plus signed `imbalance_km` | Use physical single-pass arm lengths. If `dL` crosses zero in a range, the known zero-laser-noise extremum is included. |
| Fiber noise | `l`, `fc1_hz`, or `[line.spectrum]` | Fit or supply Eq. 6 coefficients, or use a direct measured arm PSD. A direct arm PSD describes only the measured compensation state. |
| Loss | `attenuation_db_per_km`, or complete `loss_a_db` and `loss_b_db` | Use per-km attenuation for length scans. Use complete arm losses, including insertions, for a point measurement. Complete losses do not identify reach versus length. |
| Detectors | scalar `efficiency`, `dark_count_rate_hz`, `error`, or two detector `channels` | Use calibrated receiver data or datasheet values. If two detectors are supplied, the unchanged scalar protocol is run twice, once per detector projection. No averaging is performed. |
| Actuator | `omega_a_rad_s`, `bandwidth_hz`, or `[actuator.response]` CSV | Use a measured pole, a nominal bandwidth converted by `omega_a = 2*pi*f_a`, or fit a one-pole response. This is an engineering one-pole proxy. |
| Classical loop | automatic `g`, or manual `loop.g_per_s` | Automatic mode maximizes key rate subject to `g <= safety_fraction * g_crit`. Manual `g >= g_crit` is rejected. |
| Dual-band detection noise | `s0`, `fc2_hz`, `lambda_s_nm`, `lambda_q_nm` | Supply the measured or published phase-detection noise and wavelengths for the dual-band measurement path. |
| Operation | `sigma_limit_rad`, `tau_max_s`, `tau_ps_s` | `tau_max_s` is an upper cap. The working `tau_Q` is solved from the phase threshold. |
| Protocol | SNS-AOPP or CAL key-rate fields | Use the published or measured protocol parameters for your run. Unused protocol fields are neutral storage only. |

Detector PSD files are not accepted. Actuator CSV is not a PSD; it must contain
frequency, magnitude and phase.

## CSV Metadata

Laser and line PSD CSV files need positive frequencies and linear PSD values.
The metadata must state quantity, units, sidedness and pass:

```toml
[line.spectrum]
file = "my_line.csv"
mode = "fit"              # or "direct"
quantity = "phase"        # "phase" or "frequency"
frequency_unit = "Hz"
psd_unit = "rad^2/Hz"
sidedness = "one-sided"   # or "two-sided"
pass = "single"           # or "round-trip"
measurement_length_km = 86.0
```

Two-sided positive-frequency bins are converted to one-sided PSD by multiplying
by 2. Frequency PSD is converted to phase PSD with the Di Domenico relation.
Round-trip PSD needs an explicit calibrated factor; the tool does not infer a
factor of 4 from the word `round-trip`. Uncalibrated SSB dBc/Hz is rejected.

Fit mode uses F1 for lasers or Eq. 6 for lines. It checks identifiability and
rejects large residuals. The default residual gate is RMS <= 0.5 dex, an
approved engineering criterion, not a statistical confidence level.

Direct mode interpolates the measured PSD in log-log coordinates and inserts it
directly into Eq. 5 or Eq. 7. Outside the measured frequency band, the default is
rejection. Explicit endpoint power-law extrapolation is available, and the
report flags it as an engineering assumption. Without a parametric fit, inverse
requirements and factor-two sensitivity are omitted for that node.

Actuator response CSV uses:

```toml
[actuator.response]
file = "response.csv"
frequency_unit = "Hz"
magnitude_unit = "linear" # or "dB"
phase_unit = "deg"
maximum_rms_phase_deg = 5.0
```

The fitted response is the unit-DC one-pole model `1/(1+s/omega_a)`. The phase
gate defaults to 5 degrees and the magnitude gate to 0.5 dex. Both are
user-approved engineering gates, not publication-derived hardware tolerances.

## Missing Inputs And Conditional Requirements

Noise amplitudes may be left unknown: `r3`, `r2`, `l`, `s0`, `C2`, `C3`, `C4`,
and detector dark counts. The tool sets each missing amplitude to zero only to
produce an optimistic upper key-rate estimate. It then searches for the maximum
value that would reduce that optimistic rate by the configured loss fraction
(`requirements.loss_fraction`, default 0.10).

These requirements are one-at-a-time conditional limits. They assume every other
unknown contribution remains zero. Do not combine separate maxima as a joint
equipment budget.

If `keyrate.detector_error` or `keyrate.f_error` is missing, the first report
line says the result is an upper model estimate. Ideal limiting values are used
only for that estimate. If `requirements.target_key_bps` is supplied, the tool
also finds conditional limits for the missing protocol quantity. The `e_d`
requirement fixes `f_EC = 1.16`, a literature reference from Liu et al.,
Quantum Frontiers 2, 16 (2023), Section 2.1 discussion of Eq. 1. This value is
not used as a silent apparatus default in the main calculation.

## Reach

Reach scans are optional and off by default:

```toml
[reach]
enabled = true
geometry = "fixed_imbalance"
points = 5 # numerical bracketing; roots are refined near crossings
minimum_total_km = 100.0
maximum_total_km = 1000.0
```

If explicit bounds are omitted, the tool brackets from the input length. The
scan keeps the signed arm imbalance fixed and re-evaluates the actual key-rate
curve. A factor-of-two rate uncertainty is converted into kilometers by finding
the numerical crossings at `R_target`, `R_target/2`, and `2*R_target`. It is not
derived from a closed-form loss law.

Measured-phase runs and point runs with complete arm losses do not identify a
length dependence. They report reach as undefined.

Stage D timing evidence is in `results/stage3_d/cli_performance.json`: fresh
CLI processes, including Markdown/HTML/JSON and plots, using the unchanged
4097-point fast grid. The optional length scan uses five bracketing points and
refines actual crossings. It re-optimizes gain at every evaluated length; it does
not replace the length dependence by a fitted rate law. Exact transfer caching
and reuse of already evaluated gains affect scheduling only. Reference validation
retains 65537 points. Timings apply to the documented input fixtures, not arbitrary
CSV sizes or arbitrarily many nonmonotone input intervals.

Independent threshold crossings use `reach.workers = 4` by default. Set it to
one for a restricted environment. Gain/frequency arrays within each crossing
remain vectorized. The local web console uses one worker throughout.

For a range input, positive lower-rate crossings use the highest-noise endpoint
within each fixed nonmonotone corner, under the already reported B4 monotonicity
assumption. Nonpositive signed rates still evaluate both endpoints: their raw
values need not be monotone. This preserves the original signed zero envelope.

The seven fresh CLI checks, including point and range inputs in both schemes,
measured at most 4.877 s without a length scan (target 5 s) and 8.674 s with
one (target 10 s). The range fixture is an explicitly labelled factor-two
stress on `l`, not published uncertainty for an apparatus. The largest
recalculated relative target-rate error at the reported roots was `2.48e-6`
(tolerance `1e-3`).

Six expensive report/reference checks belong to the full test group. The
historical slow quick-run investigation was closed at the user's request; no
cause is asserted and no further investigation is scheduled. Archived timings
remain in `results/stage3_d/`.

## Reading The Report

The Markdown and HTML reports are ordered for laboratory decisions:

1. Fiber-stabilization ceiling and realized fractions.
2. Selected key rate, `tau_Q`, duty cycle and estimate label.
3. Ranked sensitivity.
4. Conditional requirements for unmeasured amplitudes or protocol quantities.
5. Recommended classical loop gain, `g_crit`, and gain margin.
6. Phase-variance contributions by laser, fiber and detection terms.
7. Tool limitations and input provenance.

The ceiling message answers whether fiber stabilization is worth improving at
that length. `H = R_perfect / R_free` is the ideal fiber-stabilization ceiling.
The realized fraction for a scheme is `(R_scheme - R_free) / (R_perfect -
R_free)`. `R_classical` is optimistic because the model has no phase-detection
noise for the classical round-trip beat measurement.

Sensitivity changes one approved input at a time. Noise amplitudes, dark
counts, detector error, line attenuation and absolute imbalance are reduced by
two. Detector efficiency halves the shortfall to one. Clock rate and actuator
pole frequency are doubled. Frequency cutoffs, wavelengths, refractive index,
operation thresholds and protocol probabilities are not ranked. Unknown
amplitudes are not ranked; their conditional requirements are shown instead.

The recommended classical `g` is the optimization result under the selected
safety fraction. Its margin is `1 - g/g_crit`, a gain margin fraction, not a
phase margin in degrees or a hardware stability certificate.

## Comparison Mode

A comparison file lists two or three complete configurations:

```toml
[comparison]
configurations = ["cavity_classical.toml", "cavity_dual.toml"]
labels = ["Classical round-trip", "Dual-band"]
working_total_lengths_km = [50.0, 100.0, 150.0, 200.0, 250.0]
workers = 1
```

Run:

```bash
python -m tfqkd.lab_run examples/part_c/compare_two.toml --output results/my_compare
```

Each variant receives its own report under `variant_1/`, `variant_2/`, etc. The
working-length table reports key rate, `tau_Q`, duty cycle, and ratios to the
first variant. Measured-phase inputs, direct arm residual PSDs and complete-loss
point inputs are rejected for working-length comparisons because they do not
identify another length or compensation state.

Part D laboratory examples are expected at:

- `examples/laboratory/diode_classical.toml`
- `examples/laboratory/stabilized_dual.toml`
- `examples/laboratory/actuator_comparison_10km_arms.toml`
- `examples/laboratory/actuator_comparison_100km_arms.toml`
- reports under `examples/laboratory/reports/`

## Validation Status

Run checks from the project root:

```bash
python -m unittest discover -s tests
python scripts/validate/validate_stage3_d.py --regression
python scripts/validate/validate_stage3_d.py --full
```

The full group includes the long T1–T7 source checks and the six expensive
reference/report tests omitted from quick development runs.

The spectral/integration/time/duty calculation was independently rewritten in
`scripts/validate/audit_independent.py` without importing the core. It agrees
with the core below `2e-7`; evidence is in
`results/stage3_d/independent.json` (maximum relative difference
`1.978e-7` on the reference grid).

The direct PSD regression feeds spectra generated from Table III through direct
CSV mode and reproduces the parametric path within grid tolerance; evidence is
in `results/out_of_sample/regression/direct_regression.json` and
`results/stage3_d/regression/direct_regression.json`.

The held-out Zhou2023 check calibrates one common intrinsic misalignment `e_d`
per detector projection at one shorter length with fixed `sigma = 0.05 rad`,
then predicts the other shorter length. The predicted-to-measured ratios are
0.8902..1.0786, an 8-11% rate error, within the requested factor 1.5. Mapping
that rate error through the published Table S5 loss curve gives 3.984..4.227 km
of equivalent distance error. At 615.59 km, the calibrated model predicts no
secure key while the publication reports 0.32 bit/s, so no multiplicative
correction is introduced.

These experiment checks primarily test protocol and loss bookkeeping under the
stated phase assumptions. They do not validate the spectral phase-noise forecast
for those experiments. The fitted `e_d` also absorbs finite-key differences
because the model is asymptotic.

Changing the approved phase endpoint changes positive calibrated Zhou rates by
4.58..7.27% and Pittaluga rates by 5.50..5.72%. The longest Zhou point has zero
calibrated rate at both endpoints, so its relative change is undefined.

Pittaluga2021 has one CAL point. The tool can fit the `e_d` needed to reproduce
that rate, but no held-out CAL point or own CAL QBER was located in the supplied
material. This is a diagnostic, not validation.

## What the tool does not account for

The report includes only tool limitations that affect user runs. The main ones
are:

- asymptotic key rates; finite-size effects are not included;
- the classical round-trip calculation omits classical phase-detection noise,
  so `R_classical` is optimistic relative to dual-band comparisons;
- the actuator model is one pole plus the adopted Williams delayed loop, not a
  full apparatus transfer function;
- two unequal detectors are represented as two scalar projections, not as a
  strict two-detector receiver model;
- direct PSD mode cannot infer another compensation state, inverse coefficient
  requirement, or coefficient sensitivity for that direct node;
- measured phase mode is a conditional key estimate from an observation, not a
  spectrum or reach forecast;
- separate inverse requirements for missing parameters are conditional and not a
  joint tolerance budget;
- T4 line-spectrum validation has unresolved level and shape disagreement;
  the project records this as a validation limit, not a fitted correction.

Detailed limitation text is maintained in `docs/UNKNOWNS.md`. User reports
include the "Tool limitations" section only, not the research-stage assumptions.

## Sources

- [Bertaina2024](https://doi.org/10.1002/qute.202400032): spectra, phase variance, dual-band noise
  separation, protocols and Table III example coefficients.
- [Author notebook](../sources/data/QKD.ipynb): author implementation and protocol reference.
- [Di Domenico2010](https://doi.org/10.1364/AO.49.004801): frequency-noise PSD and linewidth
  conversion.
- [Williams2008](https://doi.org/10.1364/JOSAB.25.001284): round-trip delayed-loop residual model.
- [Clivati2022](https://doi.org/10.1038/s41467-021-27808-1): measured 114 km line-noise validation input.
- [Pittaluga2021](https://arxiv.org/abs/2012.15099) and [Zhou2023](https://arxiv.org/abs/2208.09347):
  measured-phase and held-out protocol/loss checks.
- Corning SMF-28 Ultra PI-1424-AEN, July 2025, page 2: reference group index
  `n = 1.4682`: [datasheet](https://www.corning.com/content/dam/corning/media/worldwide/coc/documents/Fiber/product-information-sheets/PI-1424-AEN.pdf).
- Liu et al., Quantum Frontiers 2, 16 (2023), Section 2.1: conditional
  reference `f_EC = 1.16` for inverse `e_d` requirements: [publication](https://link.springer.com/article/10.1007/s44214-023-00039-9).

Use `docs/UNIVERSAL_INPUTS.md`, `docs/MEASURED_INPUTS.md`, and
`docs/PART_C_OUTPUTS.md` for lower-level field details and provenance notes.

## Web input validation and published examples

The local form checks supplied values before sending a calculation. The server
repeats these checks for Advanced TOML and comparison variants. Errors appear
beside the relevant field; errors in a comparison configuration identify the
parameter beside that variant's editor. Help buttons work on hover and keyboard
focus. The small caption retains the internal parameter name.

Web bounds (interface policy requested by the user, not new physical equations):

- Efficiencies and basis probabilities: `(0, 1]`.
- SNS sending probability epsilon: `(0, 1)`, since both sending and not-sending
  events are needed by the author AOPP expression (Bertaina Appendix A / notebook
  Cell 21). At 0 or 1 the field explains this requirement.
- Intrinsic misalignment error: `[0, 0.5)`; observed phase RMS: `[0, pi]` rad.
- Error-correction inefficiency: `f_EC >= 1`.
- Lengths, complete arm losses and pulse rate: strictly positive; noise-count
  rates: nonnegative. Operating windows are positive; stabilization downtime may
  be zero.
- SNS decoys: `big > medium > mini >= 0`.

An RMS above 1 rad, efficiency below 0.05, f_EC above 2, or epsilon outside
0.05–0.3 gives a warning and permits calculation while still inside its valid
bounds. These warning bands are explicitly user-selected screening policies;
they do not define equipment uncertainty or certify the validity of all other
model approximations. CLI behavior and the scientific kernel are unchanged.

### Intensities: form versus Advanced TOML

Enter **SNS** signal and decoy intensities **per user, as usually published**.
The form multiplies each by two for `keyrate.decoy_big`, `decoy_medium` and
`decoy_mini`, following Bertaina Appendix D and QKD.ipynb Cell 21. For example,
Zhou's 0.493 / 0.105 / 0.0002 per-user values become
0.986 / 0.210 / 0.0004 in the kernel.

**CAL** `u_cal` already denotes the one-user intensity alpha squared (Bertaina
Appendix B / QKD.ipynb Cell 19), so it is **not doubled**. The web report and its
JSON list the per-user input, core value, conversion factor and source for each
active intensity. These values are displayed also when entered as Advanced TOML.

The collapsed **Advanced — Configuration editor** uses the **core convention**:
SNS values there are already two-user values. Editing that text switches the
run to the edited configuration; it is not doubled again. Editing form inputs
or choosing **Build configuration** rebuilds the configuration from the form.
Additional TOML is also collapsed under Advanced.

When both complete arm losses are present, per-km attenuation is dimmed and not
used. Advanced input follows the same priority, including attenuation supplied
under `keyrate`; the ignored representation is recorded in the report. Supplying
only one complete loss is an error beside the missing arm field.

### Load example

The button fills explicit examples, never installation defaults:

- **Zhou 2023, 403 km:** geometry and losses from Supplementary Tables S1/S2,
  receiver inputs from Note 3 / Table S3, per-user intensities and probabilities
  from Table S4, f_EC from Methods. The quantum-channel residual RMS at that
  length is **not published**: the loaded sigma=0 is an explicit ideal phase
  bound, not a measurement; intrinsic e_d remains blank and gives the existing
  labelled upper estimate. The 200 ns window and zero extra stabilization time
  preserve the previously accepted timing idealization and effective 500 MHz
  clock convention. These qualifications are shown beside the loader and in
  the generated report.
- **Bertaina 2024:** explicit noise/protocol coefficients from Table III,
  Table II/Appendix D and QKD.ipynb. The actuator pole remains an engineering
  demonstration value, not a measured published specification.

Sources are the local [Zhou paper](https://arxiv.org/abs/2208.09347),
[Bertaina paper](https://doi.org/10.1002/qute.202400032) and
[author notebook](../sources/data/QKD.ipynb). The held-out check is preserved;
see web regression evidence (`../results/web_input_update/details.md`; research archive not included in the public snapshot).

## Mode-pairing asymptotic key extension

The separate `tfqkd.mp_key` entry point implements the base IID mode-pairing
protocol, [Zeng2022](https://doi.org/10.1038/s41467-022-31534-7), Eqs.4 and7.
It does not change the TF-QKD calculation or implement Zhou's async vacuum-key
term or click filtering. Run the explicit published-count example with:

```bash
python -m tfqkd.mp_key examples/mode_pairing/key_zhang202_given_gains.toml --output results/my_mp_key
python -m scripts.validate.validate_mp_key
```

The example is a calculation **conditioned on observed detection gains**, not a
forecast from detector efficiency, loss and phase noise alone. Required inputs
are click probability, maximum gap in rounds, signal Z-pair fraction `r_s`,
per-user per-pulse `mu` and `nu`, `e_z`, `f_ec`, and 3x3 Z/X gain and X error-gain
matrices in vacuum/weak/signal order. These gains must use a common pair
normalization within each basis (Zeng Supplement Eqs.53,54,56,69,70).
The example's original counts, source probabilities, phase-sifting factors and
missing-data handling are recorded in
[comparison inputs](../sources/data/mode_pairing/key_comparison_inputs.json)
and [detailed validation](../results/mp_key/details.md).
No missing apparatus or protocol parameter is silently supplied.

### Phase input and privacy inference

`phase.mode = "published_drift"` requires `sigma_L_rad_s`, the Lorentzian
`linewidth_hz` of **each** of two equal independent lasers, `phase_slices`, and
explicit `intrinsic_x_error`. Zhang Eq.3 describes the drift/linewidth model
under ideal frequency compensation; it does not specify real tracking errors.
`phase.mode = "spectral"` uses the accepted `mp_phase` components, numerical
grid and an explicit `T_track_s` if a component's increment integral diverges.
Only divergent components receive the approved 1/T_track cutoff, which is a
band convention rather than a frequency-estimator model. The phase module can
separately evaluate tenfold tracking-window sensitivity.

`intervals.mode = "iid_zeng"` uses the truncated geometric law (Zeng Eqs.15/16),
with durations at the **active** round clock. It omits minimum pairing gaps and
reference-frame gaps. Alternatively, `mode = "measured"` requires arrays
`delta_t_s` and `counts` and uses their discrete sum (Zhang Appendix D5).
Never substitute bin midpoints when the interval distribution is unknown.

Raw coherent X errors follow Zhang B10 and its D2 misalignment factor. In
`decoy.error_mode = "phase_weak_equal"`, this replaces only the balanced weak
X-setting error gain; all other gains/error gains remain explicit. In
`"observed"` mode the phase calculation is a diagnostic and does not replace
measured error gains. The single-photon privacy bound follows Zhang C4–C7,
with vanishing finite-key sampling term. **The coherent 25% floor is not
subtracted to obtain a single-photon error.** A nonpositive yield lower bound
or an error upper interval containing 1/2 supplies no positive privacy bound.
Negative error-yield upper bounds are rejected as inconsistent inputs.

### Units and limits of validation

Output is bits per shared quantum emission round. Conversion to bit/s uses
`active_round_clock_hz * (1-reference_fraction-recovery_fraction)`; all three
inputs are explicit. The denominator is neither a successful pair nor one
user's pulse. Zhu's `K/N_pair`, with `N_pair=N_rounds/2`, is twice the per-round
rate (Zhu Eq.1 and its following paragraph). Reference slots must not be counted
again when an already effective quantum clock is supplied.

Published drift parameters in [Zhang2025](https://doi.org/10.1103/PhysRevX.15.021037),
[Zhu2023](https://arxiv.org/abs/2208.05649) and
[Zhou async](https://arxiv.org/abs/2212.14190) come from laboratory
fiber spools. **Do not transfer them to field routes.** The four published Zhou
drifts give alpha=1.168 with a conditional residual-based 95% interval
[0.743,1.592]; alpha=0.5 is excluded within that regression, but published
measurement uncertainties are unavailable. This supports rejecting the
independent-section scaling for this spool dataset under the stated fit
assumptions; it does not prove spatial correlation or identify a universal
coefficient. See [spatial analysis](../results/mp_spatial_analysis/details.md).

The new checks use published finite counts as asymptotic expectation proxies;
they are **not finite-session security certificates**. Analytical decoy bounds
can be looser than a publication's LP bounds, so removing finite-key penalties
does not guarantee a higher rate across different estimators. At Zhu's 304 and
407 km the present analytical bounds establish no positive key. Zhu phase
forecasts are deliberately absent: drift and tracking residuals are unreported.
Zhang's one-parameter held-out checks remain conditional on measured detection
gains at each length. Their fitted `e_d` absorbs estimator, timing, tracking and
finite-key differences and is not a measurement of intrinsic misalignment.
All discrepancies and the literal published-table anomaly are retained in the
[MP key validation report](../results/mp_key/details.md); no correction factor
is applied.

### Mode-pairing phase results and equipment-transfer limit

For the single-pass fiber Eq.6 of Bertaina, the phase-increment integral has the
analytically derived short-interval limit

\[
D_F(\Delta t)=4\pi^2 l(L_A+L_B)f_{c1}\Delta t^2
\left[1-\frac{\pi^2}{3}f_{c1}|\Delta t|
+O((f_{c1}\Delta t)^2|\log(f_{c1}|\Delta t|)|)\right].
\]

Thus the leading drift-model parameter is
`sigma_L = 2*pi*sqrt(l*(L_A+L_B)*fc1)`, when `fc1*|delta_t| << 1` and the
convergent fiber increment is integrated from zero. Together with white
frequency noise from two equal independent lasers, the leading total variance
becomes Zhang Eq.3, `D = 4*pi*linewidth*|delta_t| + sigma_L^2*delta_t^2`.
This is an asymptotic connection, not a new long-interval model. The cubic
correction is negative; the derivation and numerical checks are in the
[phase report](../results/mp_phase/details.md), based on Bertaina Eq.6 and
Di Domenico Eq.1.

Using another line's published `l=44` and `fc1=100 Hz` without fitting predicts
drifts **1.592–2.820 times above** Zhou async's published empirical drift inputs;
Zhang's 403 km comparison has ratio2.092. This establishes order of magnitude
for these conditional comparisons, not calibrated coefficients for a spool or
field route. The four Zhou values have log-slope alpha1.168 with conditional
95% interval[0.743,1.592]: the independent-section alpha0.5 is outside that
interval. Measurement uncertainties are not published; the residual-based
regression does not identify a spatial-correlation mechanism. See the
[original source/fit analysis](../results/mp_spatial_analysis/details.md).

### Zhu diagnostics: distinct rate discrepancies

The accepted **asymptotic analytical-decoy** result at202 km is1569.78bit/s,
−24.53% against2080bit/s. The separate −12.17% discrepancy belongs to a replay
of Zhu's **published finite components**, which produces1826.80bit/s from
TableVI/VII and Eq.1/A1. The two comparisons have different inputs and must
not be confused.

At202 km the printed single-photon lower count and phase upper give8.93637
million privacy bits, minus0.31223million EC bits. The published2080bit/s
requires9.81948million final bits at the stated normalization. Even zero EC
expense gives only1892.94bit/s, so error correction or an additional
nonnegative finite-key penalty cannot close the difference. Exact sent counts
and printed-value rounding also fail to close it. The inconsistency is in
the relation between the published privacy components and final rate at that
normalization; the source does not identify a unique erroneous scalar. No
entry is corrected to make the rates agree.

Zhu uses a minimum pairing gap63 and Eq.B7 (PDFp11), whereas the base module
uses Zeng Eq.4. B7 matches Zhu TableVI `r_p` within0.53%. In these observed-count
checks `r_p*r_s = M_signal/N_rounds`, so this prefactor difference cancels and
does not explain the finite-component replay discrepancy. B7 is evaluated
only by the diagnostic script; no new protocol option is added.

Our analytical `q_(1,1)` lower is0.75–4.47% **higher** than the published finite
fractions across all four Zhu lengths. At304 km the analytical X yield lower
is negative, leaving only a trivial error upper1. At407 km the error upper is
0.555958 versus the published0.3468. Both include the worst privacy entropy
at1/2 and establish zero key. The looser X/privacy bound, rather than the Z
single-photon fraction, is the limiting estimate.

The source documents decoy-state estimation and Chernoff-Hoeffding bounds
(Zhu AppendixA step9/p6, A2/p7, A3/p8), using weak-weak X pairs, and cites
Zeng2022. Zeng Supplement Note3 Eqs.64/65 (PDFp33) formulates LP problems and
also describes an analytical solution onp34. **Zhu does not identify which
solver it actually used.** References to "published finite-key LP estimates"
in the earlier MP comparison must not be read as verification of Zhu's solver.
The numerical intermediate X single-photon/error bounds are also unreported.
No alternative estimator is introduced in v1.1.

Full term comparisons, all four decoy comparisons, diagnostic one-term
replacements and their source locations are in
[v1.1 diagnostics](../results/mp_zhu_diagnostics_v1_1/details.md) and JSON.
