# Measured phase, spectra, losses, and detectors

The tool accepts an observed residual phase, a measured spectrum without a
fit, complete arm losses, and two separate detector descriptions. It does not
import Bertaina apparatus parameters into an unspecified installation.
Key rates use the asymptotic protocol model.

## Residual phase observation

Use `[phase] mode="measured"` and explicitly supply all three quantities:

```toml
[phase]
mode = "measured"
sigma_phi_rad = 0.104 # rad; Pittaluga Fig.2f at 605 km, NOT a measurement at another length
tau_s = 0.00002504 # s; Pittaluga Note VI frame duration; conditional comparison window
tau_ps_s = 0.0 # s; explicit timing idealization, not an apparatus default
source = "Pittaluga2021 Fig.2f and Note VI; explicit phase-bound/window assumptions"
description = "Phase bound from another length, not a same-window measurement."
```

The complete executable example is
[pittaluga2021_368.702_phase_bound.toml](../examples/measured_inputs/pittaluga2021_368.702_phase_bound.toml).
It supplies complete losses, both detectors and CAL parameters as well.

```bash
python -m tfqkd.lab_run examples/measured_inputs/pittaluga2021_368.702_phase_bound.toml --output results/my_phase_run
```

The spectral chain is bypassed. The reported variance is sigma squared;
`e_phi = sigma_phi²/4` follows Eq.1, bertaina2024. Duty is
`tau/(tau+tau_PS)`, using the author duty formula with your observed window.
There is no threshold-time solve, phase cap, required spectral cutoff, laser
lock shape, fiber model or actuator pole. No propagation-index default is used
when no propagation calculation is made. Complete measured losses make physical
lengths optional for this point calculation.

The protocol still uses its original Gaussian phase-error function from the
notebook; reporting Eq.1 does not replace that function. This preserves the
validated core. Missing intrinsic detector error and error correction factor
are handled by the user-authorized ideal values zero and one, with explicit
provenance and an **upper model estimate** label. Missing scalar dark counts
retain the existing conditional requirement mechanism. Other selected-protocol
parameters remain mandatory. Missing intrinsic error or fEC is labelled in the
first report line. Set `requirements.target_key_bps` for separate conditional
protocol requirements; details and the published fEC reference are in
[USER_GUIDE](USER_GUIDE.md#unknown-protocol-errors-and-operational-reach).

This is a key estimate conditional on an observation or labelled phase bound.
It does not predict another length, another compensation scheme, a gain, a
stabilization ceiling, or separate noise components. Those sections are absent
from its Markdown/HTML report. Measured-phase range input is currently rejected:
supply and execute separately labelled endpoint files.

## Direct spectrum, without F1/Eq.6 fitting

```toml
[laser.spectrum]
file = "laser_output.csv" # resolved relative to this TOML file
mode = "direct" # mode="fit" retains the existing coefficient fitting path
quantity = "phase" # or frequency
frequency_unit = "Hz"
psd_unit = "rad^2/Hz" # frequency alternative: Hz^2/Hz
sidedness = "one-sided" # or two-sided, explicitly normalized
pass = "single" # round-trip requires an explicit calibrated PSD factor
extrapolation = "reject" # default; explicit alternative: power-law
```

The two-column CSV has a header and positive, ordered frequency/PSD samples.
Normalization metadata is mandatory in both fit and direct modes. The existing
normalizer converts frequency PSD to phase PSD and two-sided to one-sided PSD;
uncalibrated round-trip and SSB dBc/Hz inputs remain rejected.

Direct PSD is interpolated in log-log coordinates onto the calculation grid.
The default rejects a grid extending beyond the measured band. Explicit
`extrapolation="power-law"` uses each endpoint's two-point logarithmic slope;
the report records its use, slopes and measured band. This is an engineering
input representation, not a claim about unmeasured noise.

A direct laser file is the **actual output** of the laser, including any lock
already used in the measurement. F2 is not applied again. Its r3/r2/fc and
cavity coefficients/loop shape are unnecessary. Eq.5 or Eq.7 combines that
output with the arm terms. The existing core still assumes identical spectra
for two independent lasers.

A `[line.spectrum]` input has the same metadata plus
`measurement_length_km`, the length represented by its normalized single-pass
arm PSD. It must be a **fiber-only contribution for the selected compensation
state**, excluding laser and central detection noise. No fc1 or l is required
for this direct node. Linear scaling with length follows the explicit Eq.6/8
spatial assumption; it does not identify another route's noise. Dual-band
central detection remains a separate, once-only contribution (Appendix G).

An already measured arm residual cannot identify another compensation state,
another controller gain or a four-scheme stabilization ceiling. Those forecasts
are omitted or rejected with an explanation. Positive gain queries are rejected
for such a residual. For prediction of another state, supply its measurement
separately or use `mode="fit"` with an applicable source model.

Coefficient sensitivity and inverse coefficient requirements are omitted for
each direct node, with a reason. Other explicit parametric nodes retain their
analysis. No model coefficient is inferred silently. A synthetic executable
demonstration is [direct_table3.toml](../examples/measured_inputs/direct_table3.toml);
its CSV is generated from F1, **not an experimental measurement**.

## Complete losses and separate detectors

```toml
[line]
loss_a_db = 34.21625623318439 # Pittaluga Table II + Table I insertion transmission
loss_b_db = 34.33392836595395 # complete dB, not just attenuation times length

[detector]
channels = [
  { name="D0", efficiency=0.73, dark_count_rate_hz=14.0 },
  { name="D1", efficiency=0.77, dark_count_rate_hz=14.0 }
] # Pittaluga Table I; fractions and Hz
```

Both complete losses are required together. Do not also supply a per-km
attenuation coefficient. Physical arm lengths are still required in spectral
runs for noise and propagation delay.

The author loss convention attenuates the stronger arm to the weaker arm:
effective total loss is twice the larger complete arm loss (Bertaina Sec.IV /
`calc_sigma_tau_loss`). This is an explicit equalization approximation, not a
general asymmetric protocol. A complete-loss point cannot be extrapolated in
the working-length comparison by inventing a per-km attenuation coefficient.

The author protocol accepts a scalar receiver. Two actual detectors are kept
as two scalar **projections**, each assuming that both ports have that detector's
characteristics. The report shows both raw signed bounds and nonnegative secure
rates. It does not average detectors. The smaller projection is the displayed
reference; their interval is **not a certified physical bound** for an unequal
two-detector receiver. Per-detector coefficient ranking is omitted, rather than
ranking an invented scalar average.

Optional `background_count_rate_hz` preserves a separately supplied background.
The explicitly supplied dark and background rates are added for the scalar
spurious-count input (Zhou Supplementary Note3/Table S3). Absence of a background
field is not a claim that a measured apparatus has zero background.

See [validation](VALIDATION.md) for numerical checks and discrepancies.
