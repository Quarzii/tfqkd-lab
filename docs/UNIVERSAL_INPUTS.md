# Universal laboratory inputs

The instrument has no equipment or route preset library. Every apparatus value
comes from your explicit input or an accepted CSV fit. Published equipment
examples are demonstration inputs, not values inherited by another laboratory.

Measured-input modes are documented in
[MEASURED_INPUTS.md](MEASURED_INPUTS.md): residual sigma with a measured window,
direct CSV without fitting, complete arm losses and two detector projections.
These paths have mode-specific mandatory inputs and omit predictions not
identified by the supplied observation. The remaining sections below describe
the parametric/fitted path unless explicitly stated otherwise.

```bash
python -m tfqkd.lab_run examples/bertaina2024_table3.toml --output results/my_run
python -m tfqkd.lab_run examples/missing_r3.toml --output results/my_requirement
python -m tfqkd.lab_inputs SETUP.toml
python -m tfqkd.measured_inputs examples/b6/jiang2008_urban86.json --node line
python -m tfqkd.lab_requirements SETUP.toml
python -m tfqkd.lab_ceiling SETUP.toml
python -m tfqkd.lab_uncertainty SETUP.toml --workers 1
```

The first command produces JSON, Markdown, standalone HTML, and scientific
plots. Reports include sensitivity, conditional equipment requirements,
compensation comparisons, and noise contributions. See [report details](OUTPUTS.md)
and the [user guide](USER_GUIDE.md) for the browser interface.

## What must be supplied

Start from `examples/bertaina2024_table3.toml` to understand the complete shape,
then replace its values with your own. Copying that example explicitly selects
that published calculation; it does not identify your apparatus. Do not retain
example values for quantities you do not know.

| Node | Direct input | Familiar units / alternative | Missing values |
| --- | --- | --- | --- |
| Laser | `laser.r3`, `r2`, `fc_hz`, `model="free"` or `"cavity"` | `laser.lorentz_width_hz`, or `laser.spectrum` | Missing r3/r2 gives optimistic zero noise and conditional equipment requirements; fc must be explicit. |
| Cavity laser | Additional `C4`, `C3`, `C2`, `B_hz`, `gamma`, `delta` | A Lorentz input refers to its **free plant**, not its locked output | Missing C2–C4 may give inverse requirements; loop shape parameters are required. |
| Line | `line.arm_a_km` and `arm_b_km`; `l`, `fc1_hz`, `attenuation_db_per_km` | Instead of arms: `length_km` **and** signed `imbalance_km`; or `line.spectrum` for noise | Geometry, attenuation and cutoff are required. Missing l may give a requirement. |
| Detectors | `detector.efficiency`, `dark_count_rate_hz`, `error` | Same calibrated/passport quantities; efficiency is a fraction, not percent | Efficiency is required. Missing error uses a labelled ideal upper estimate and a conditional target-rate requirement; missing dark counts may give a noise requirement. PSD files are not accepted for detectors. |
| Actuator | `actuator.omega_a_rad_s` | `bandwidth_hz`, or `actuator.response` CSV | A pole is required for classical compensation and for the four-rate comparison. No equipment type is selected. |
| Dual compensation | `physics.lambda_s_nm`, `lambda_q_nm`, `fc2_hz`, `s0` | Explicit wavelengths and phase-measurement noise | Wavelengths/cutoff are required; missing s0 is optimistic zero detection noise. |
| Operation | `operation.sigma_limit_rad`, `tau_max_s`, `tau_ps_s` | radians, seconds, seconds | All explicit. tau_max is a cap; tau_Q is solved, not fixed at the cap. |
| Protocol | `protocol.name="SNS-AOPP"` or `"CAL"`, and `keyrate` fields below | Clock [Hz], probabilities, mean photon numbers, error correction factor | Selected intensities/probabilities/clock are explicit. Missing fEC is labelled as an ideal upper estimate; see target requirements in USER_GUIDE. |

The direct core aliases `[physics]` and `[keyrate]` are also accepted. Choose one
representation per parameter: a duplicate core/node value, even an equal one,
is rejected. Both explicit arm lengths cannot be combined with total length or
imbalance. No unpublished parameter is filled from Bertaina.

Both protocols require `clockrate_hz`, a loss representation and
`detector_efficiency`. Unreported `detector_error` and `f_error` give only a
first-line labelled upper model estimate; conditional limits require an explicit
`requirements.target_key_bps`. Missing dark counts are handled as an unknown
noise amplitude. SNS-AOPP additionally requires `decoy_big`,
`decoy_medium`, `decoy_mini`, `pz_sns`, `eps_sns_aopp`. CAL additionally requires
`pz_cal`, `u_cal`, `nmin_cal`, `nmax_cal`. Unused protocol fields are neutral
dataclass storage; they are never used as apparatus values by the selected
protocol. Overhead is supplied by `operation.tau_ps_s`; a duplicate
`keyrate.stab_overhead_s` must agree.

`laser.model`, arrangement/compensation and protocol-name defaults are model
choices, not measured equipment properties. `configs/lab_defaults.toml` contains only
these choices, numerical grid/search settings, the configured fit gate and inverse
loss target. The original `configs/config.toml` serves historical core validation only;
laboratory calculations do not import its apparatus/protocol values. The optional Python
`core_config` argument transfers only computational settings and validation
scenario labels, never physical values.

### The reference index exception

If `physics.n` is omitted, the reference is **1.4682**, the typical
effective **group index at 1550 nm** in Corning SMF-28 Ultra PI-1424-AEN,
July 2025, **page 2, Performance Characterizations**. It is appropriate as a
reference propagation-delay input, not an identification of every fiber or a
phase refractive index. Override it if your wavelength/fiber warrants another
value. Output records the value, source and `default_used=true`.

The existing Table I calculation in `results/universal_inputs/core_regression/full/stage2_checks.json`
checks n=1.44–1.48 relative to n=1.45: its largest variance change is about
1.818719%. This recorded sensitivity is **not** a guaranteed bound for arbitrary
geometry, controller or spectrum. The example retains explicit n=1.45 to
reproduce the authors. c=299792.458 km/s is the author code's speed-of-light
constant; K=4 defaults to correlated passes, with K=2 available explicitly
(Bertaina Eq.5 discussion). Neither is a fitted noise coefficient.

## PSD CSV: laser and line only

CSV columns are named `frequency,psd`. Data must be finite, strictly positive
linear PSD, with positive distinct Fourier frequencies in Hz. Rows are sorted.
A file labelled two-sided supplies its positive-frequency bins; negative bins
and DC are not used by this fitting interface. The metadata is mandatory:

```toml
[line.spectrum]
file = "my_line.csv" # relative to this configuration, not the project directory
quantity = "phase" # phase or frequency
frequency_unit = "Hz"
psd_unit = "rad^2/Hz" # frequency PSD: Hz^2/Hz
sidedness = "one-sided" # or two-sided
pass = "single" # or round-trip
measurement_length_km = 86.0 # EXAMPLE measurement geometry, not a default
# round_trip_psd_factor = 4.0 # must be calibrated and explicitly supplied IF round-trip
```

A laser uses `[laser.spectrum]` without measurement length. It must represent a
free laser plant for F1 fitting; a locked output is not silently interpreted as
the free plant in F2. A line's measurement length is the single-pass path length
underlying that measurement, and may differ from the target link length.

Normalization is explicit: two-sided → one-sided multiplies the positive-bin
PSD by 2; frequency → phase divides by f² (Di Domenico Eq.1); round-trip divides
by the **user-supplied calibrated PSD factor**. Four is not inferred from the
word “round-trip”: the correlation assumptions matter (Bertaina Eq.5).
Uncalibrated SSB dBc/Hz is rejected, rather than treated as phase PSD.

Fits use exactly F1 or Eq.6 with positive coefficients and equal log-residual
weights, a log-space fitting method. The local
finite-difference Jacobian SVD checks identifiability. Rank deficiency or
unidentified coefficient uncertainty rejects the fit. RMS log residual must
not exceed **0.5 dex**, the engineering gate, configurable through
`[fit].maximum_rms_log10_residual` or that spectrum's setting of the same name.
This gate is not a source-derived statistical confidence test. Fit coefficients,
conditional local SE, residuals and normalization are retained in provenance.
SE describes scatter under the adopted local residual model; it does not become
an invented equipment population range or calibration error.

The five B6 curves and exact digitization provenance are in `examples/b6/`.
`*.json` contains the CSV specification; `*.csv` contains normalized positive
one-sided single-pass phase PSD. Pixel coordinates, image anchors and source
normalization evidence remain in `results/stage3_b_amendments/b6/` and
`examples/b6/digitization_provenance.json`. Jiang passes the configured gate;
NF6700/Droste fail residual checks and Orbits/Snigirev fail identifiability
(Orbits also exceeds the residual gate). None is an equipment selector.

## Actuator response CSV: not a PSD

Columns: `frequency,magnitude,phase`. Frequencies are positive Hz, magnitude
is positive linear amplitude or amplitude dB, phase is degrees:

```toml
[actuator.response]
file = "response.csv"
frequency_unit = "Hz"
magnitude_unit = "linear" # or dB; amplitude dB uses 20 log10
phase_unit = "deg"
maximum_rms_phase_deg = 5.0 # configurable engineering default; change as needed
maximum_rms_log_magnitude_dex = 0.5 # configurable magnitude residual gate, configurable
```

The user-adopted model has unit DC gain `1/(1+s/omega_a)`. The fit minimizes
the concatenated log amplitude residual and wrapped phase residual in radians,
with equal dimensionless weights. This is an explicitly adopted engineering
fit, not a Williams equipment-identification formula. It reports omega_a,
conditional local SE and separate magnitude/phase residuals. Either excessive
residual or unresolved pole causes refusal. The phase gate defaults to **5 degrees**, specified as an engineering
criterion, not a source value or apparatus calibration tolerance. Override it
through `[fit].maximum_rms_phase_deg` or the response metadata field of that name. Magnitude calibration
and the unit-DC convention must be appropriate to the measured tract.

Alternatively `bandwidth_hz` uses omega_a=2πf_a as the accepted **engineering
one-pole proxy**. Nominal bandwidth does not establish a full loaded actuator
response or a hardware stability guarantee. Direct omega_a remains an option.

## Missing noise: upper estimates and conditional requirements

Missing amplitudes r3, r2, l, s0, C2–C4 or detector dark counts contribute zero
to the optimistic rate. They are prominently labelled **unmeasured**, not
replaced by equipment defaults. Geometry, pole, cutoffs and protocol parameters
cannot use this rule. Missing mandatory inputs produce a collected list.

For each missing amplitude, a one-dimensional monotone search solves
`R(parameter) = (1-loss_fraction) R_upper`, reusing the selected core model.
`requirements.loss_fraction` defaults to the user-requested **0.10**. A numerical
bracket expands/contracts by factors of two; it is not an apparatus prior.
Observed monotonicity is checked with the existing numerical rate tolerance,
and the returned coefficient is recalculated to verify the achieved reduction.
A zero optimistic rate, inactive coefficient, exhausted search or model failure
does not produce an invented finite equipment bound.

**Every requirement is conditional on the other unmeasured contributions being
zero.** Do not set all missing parameters to their separate maxima and interpret
that as a guaranteed combined 10% loss. The report lists the other zeroed
parameters for each requirement. At common-laser imbalance zero the laser term
is exactly absent by Eq.5, so no finite laser-amplitude bound is inferred.

## Four-rate stabilization ceiling

Every fully specified parametric/fitted `lab_run` point (including each evaluated
range case) calculates:

1. R_free: no fiber stabilization.
2. R_perfect: the same free-fiber baseline, but **only its fiber term is zero**.
3. R_classical: the adopted Williams residual with automatically optimized g.
4. R_dual: the separated fiber and detection contributions of Eq.8/Appendix G.

Laser, geometry, protocol and operation inputs stay fixed. As explicitly agreed,
S_det appears only in R_dual; free/classical/perfect retain zero S_det. This
preserves the core's physical detection convention. The ceiling is
`H=R_perfect/R_free`; a scheme's realized fraction is
`(R_scheme-R_free)/(R_perfect-R_free)`. R_free=0 or zero headroom is reported
as infinite/undefined where appropriate, never assigned an arbitrary factor.
R_perfect below either scheme raises an **error**, not a warning or rate repair.
For H<1.1 the report states that even ideal fiber stabilization offers less
than a 10% gain under those inputs. With unknown amplitudes, this comparison is
conditional on optimistic inputs, not a guaranteed prediction for the apparatus.

Classical g obeys the configured safety fraction (default 0.5) times the numerically
found g_crit. Manual g is optional for the selected scheme; g≥g_crit is rejected.
The four-rate ceiling always optimizes its classical comparison, even when the
main selected calculation uses a manual g. Margin is `1-g/g_crit`, not degrees
of phase margin. Numerical scan/refinement is resolution checked, not a proof
of the global optimum for arbitrary apparatus.

The existing integral of free fiber noise above f_b remains **secondary**. It
is not the A8 residual or a fixed-window reach test. Eq.4 can reach the phase
threshold by shortening tau_Q, even when that diagnostic exceeds one threshold.
An optimized g=0 reports no resolved key benefit at the search tolerance; it
does not declare that no operating window exists.

## User-supplied ranges

Input ranges are supported. A range uses a supported dotted
parameter and explicit `minimum`, `maximum`, `sources`. The program records
your source references without claiming it independently verifies them. No
percentage range is added by the instrument or derived from a conditional fit SE.

Additive noise parameters passing the Table I monotonicity checks use two
endpoint calculations. Other parameters use corners; an imbalance interval
containing zero additionally evaluates that known extremum. Other potential
interior extrema retain `NONMONOTONE_INTERIOR`: the output is a sampled envelope,
not certified continuum bounds. Independent endpoint/corner cases use processes;
g/frequency work remains vectorized. Set `numerics.workers=1` if process creation
is unavailable; Python callers using processes need a normal `__main__` guard.

`lab_run` reports the selected-scheme range envelope and stores four-rate
comparisons for every evaluated case in JSON. The headline H is the value at
the explicitly supplied point; it is not presented as a guaranteed interval for
H over a nonmonotone continuum. Counts include actual input cases, comparison
cases, inverse-search evaluations and spectral gain candidates separately.

## Applicability and sources

The TF-QKD calculation uses the same PSD for both independent lasers; arbitrary
different laser spectra require a separately specified extension. Protocol
loss equalizes the shorter arm to the longer arm: `2 max(L_A,L_B) alpha`
(Bertaina Sec.IV, authors' `calc_sigma_tau_loss`). Actual physical lengths enter
fiber noise and delay. The classical residual assumes uniform, uncorrelated
spatial noise (Williams A8/A11); applying it to the common-laser geometry with
the original K is the already documented engineering composition.

See [model limitations](UNKNOWNS.md) and [validation](VALIDATION.md)
for assumptions and discrepancies. Classical measurement-detection noise
is absent: R_classical and its realized ceiling share are optimistic, and the
comparison is biased in favor of classical compensation.

- [Bertaina 2024](https://arxiv.org/abs/2310.08621): Eqs.1,4–8, F1–F4,
  Appendix G; Table II/III, Appendix D; protocol Appendices A/B.
- [Author notebook/data](https://zenodo.org/records/10911121): `QKD.ipynb`,
  `calc_spectra`, `stabfibnoise`, `calc_sigma_tau_loss`, `G1`.
- [Di Domenico 2010](https://doi.org/10.1364/AO.49.004801): Eqs.1,5,
  one-sided frequency noise and pure-white Lorentz width.
- `williams2008.pdf`: Appendix A, especially A6/A8/A11, JOSA B 25,1284.
- [Corning SMF-28 Ultra](https://www.corning.com/content/dam/corning/media/worldwide/coc/documents/Fiber/product-information-sheets/PI-1424-AEN.pdf):
  July 2025, page 2, effective group index at 1550 nm.
- [Jiang 2008](https://arxiv.org/abs/0807.1882): Sec.3, Figs.1/3,86-km Paris line.
- [Li 2013](https://thesis.caltech.edu/7799/): Appendix B, Eq.B.1 and Fig.B.1,
  PSD normalization and ECDL/fiber-laser curves.
- [Snigirev 2023](https://doi.org/10.1038/s41586-023-05724-2): Fig.2(h),Methods,
  one-sided noise of the free DFB; original [data](https://zenodo.org/records/7371066).
- [Droste 2013](https://doi.org/10.1103/PhysRevLett.111.110801): Fig.4,Eq.2,
  single-pass intercity-line spectrum.

## SNS decoy normalization

`decoy_big`, `decoy_medium`, `decoy_mini` are **total two-user** mean photon numbers, following bertaina2024 Appendix D and QKD.ipynb Cell 23. For symmetric setups, multiply a publication’s per-user decoy intensities by two before input. The core fixes the per-user signal and not-send intensities to `decoy_big/2` and `decoy_mini/2`. Two asymmetric source sets are not supported.
