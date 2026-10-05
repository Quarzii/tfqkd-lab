# User guide

Use this tool to estimate phase noise and asymptotic key rates from explicit
installation data. Supplied examples reproduce particular calculation inputs;
they do not describe your equipment until you replace their values.

## Installation

Use Python 3.11 or newer. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, use `.venv\Scripts\activate` instead of the activation command above.

## TF-QKD in the browser

```bash
python -m tfqkd.lab_web --port 8765
```

Open <http://127.0.0.1:8765>. Choose SNS-AOPP or CAL, enter your installation
parameters, and calculate a report. Loaded examples are explicitly identified
calculation cases. The Zhou example uses an ideal phase bound rather than a
published quantum-channel RMS measurement at that length.

The server listens on your computer's loopback address. Uploaded CSV files are
kept temporarily, for at most one hour with up to eight reports, and deleted on
exit. Save any report you need before closing the server with Ctrl+C.

In the form, SNS signal/decoy intensities are **per user**. In TOML,
`decoy_big`, `decoy_medium`, and `decoy_mini` are **total two-user intensities**:
for symmetric sources, multiply published per-user values by two. The form
performs this conversion; editing the advanced TOML does not apply it again.
CAL `u_cal` is already a one-user intensity and is not doubled.

## TF-QKD from a file

Choose a [laboratory example](../examples/laboratory/README.md), replace its
values, then run:

```bash
python -m tfqkd.lab_run examples/laboratory/stabilized_dual.toml --output results/my_run
python -m tfqkd.lab_run --help
```

For a spectral prediction, supply laser and fiber coefficients or measured
PSDs, geometry, losses, detector data, timing, and the selected protocol inputs.
Classical compensation also needs an actuator pole or measured response.
The [input reference](UNIVERSAL_INPUTS.md) lists required fields and units.

For a measured-phase estimate, supply residual RMS in radians, its operating
window and stabilization overhead in seconds, and receiver/loss/protocol data.
This estimates key rate at that observation; it cannot infer another line's
noise or compensation performance. See [measured inputs](MEASURED_INPUTS.md).

PSD CSVs require positive frequencies, linear PSD values, and explicit quantity,
units, sidedness, and single/round-trip normalization. A linewidth alone does
not specify low-frequency colored laser noise. Direct PSD interpolation rejects
out-of-band evaluation unless extrapolation is explicitly selected.

Missing required geometry or spectral-shape parameters produce an input error.
Certain unknown noise amplitudes, intrinsic error, or error-correction factors
instead give a labelled upper model estimate. Conditional equipment limits
assume the other unknown contributions are zero; separate maxima are not a
joint equipment budget. If omitted, `physics.n` uses the documented reference
group index 1.4682 at 1550 nm, rather than a measured value for your fiber.

## Read the TF-QKD report

The output directory contains `report.md`, `report.html`, and `result.json`.
Spectral runs also write `spectrum.png` and `cumulative_variance.png`;
measured-phase runs have no inferred spectrum plots.

- **Key rate:** asymptotic estimate in bit/s for the selected protocol.
- **Operating window (`tau_Q`):** usable transmission time in seconds, solved from the phase threshold; `tau_max_s` is a cap.
- **Duty cycle:** fraction of the cycle available for key transmission.
- **Noise contributions:** laser, fiber, and detection variance at the same selected window.
- **Compensation comparison:** free fiber, ideal fiber cancellation, classical round-trip compensation, and dual-band stabilization.
- **Sensitivity and requirements:** conditional effects of changing specified inputs, rather than certified equipment tolerances.

Two supplied detectors produce two scalar receiver projections. Their range
is not a certified bound for a receiver with unequal detectors. Complete arm
losses describe a point condition and do not establish attenuation for a length
scan. Optional reach scanning is disabled by default and requires suitable
length-dependent inputs. See [report details](OUTPUTS.md).

## Mode-pairing QKD

Run phase increments separately from key estimation:

```bash
python -m tfqkd.mp_phase examples/mode_pairing/phase.toml --output results/my_mp_phase.json
python -m tfqkd.mp_key examples/mode_pairing/key_zhang202_given_gains.toml --output results/my_mp_key
```

The phase example evaluates independent single-pass laser/fiber components at
specified pairing intervals. Its transferred fiber coefficients are demonstration
inputs, not calibrated coefficients for the cited mode-pairing experiments.
A tracking window supplies the `1/T_track` lower band edge only for divergent
terms; it does not describe a frequency-tracking response.

The key example implements the base IID Zeng protocol. Supply click probability,
maximum gap in rounds, signal Z-pair fraction, per-user per-pulse intensities,
Z error, error-correction factor, and 3×3 Z/X gain and X error-gain matrices in
vacuum/weak/signal order. Gains require consistent pair normalization within
each basis. Original counts and normalization are recorded in
[comparison inputs](../sources/data/mode_pairing/key_comparison_inputs.json).

Phase inputs can use explicit published drift and linewidth or a spectral
calculation. Observed error gains can remain unchanged, with phase used only as
a diagnostic; `phase_weak_equal` instead replaces the balanced weak X error
gain. Measured interval counts can replace the IID interval law. A coherent
X-error floor is not subtracted to invent a single-photon privacy error.

Rates are bits per shared quantum emission round. Conversion to bit/s uses
`active_round_clock_hz * (1-reference_fraction-recovery_fraction)`.
Do not count reference time again if the supplied clock is already effective.
The base calculation omits experimental minimum gaps, frame gaps, and async
vacuum-key terms. Laboratory-spool drift coefficients are not field-route
defaults. See [mode-pairing examples](../examples/mode_pairing/README.md).

## Check results and assumptions

All key models are asymptotic and do not certify finite-session security.
Classical compensation omits phase-measurement detection noise and uses a
one-pole actuator approximation. Consult [model limitations](UNKNOWNS.md) and
[validation results](VALIDATION.md) before interpreting a calculated rate as an
experimental forecast.

```bash
python -m unittest discover -s tests
TFQKD_FULL_TESTS=1 python -m unittest discover -s tests
python scripts/validate/validate_stage3_a.py --calibrate
```

These checks do not require PDFs. Source-dependent reproduction checks require
the separately obtained papers listed in [SOURCES.md](../SOURCES.md).
