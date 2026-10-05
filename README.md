# Fiber QKD noise and key-rate tool — v1.1

[Русское описание / Russian](README.ru.md) · [User guide](docs/USER_GUIDE.md) · [Sources](SOURCES.md)

Estimate how laser and fiber phase noise limits quantum key distribution over
optical fiber. This Python tool connects measured or explicitly specified
installation parameters to an operating window, duty cycle and asymptotic key
rate. It helps identify the dominant noise contribution and compare possible
equipment improvements before changing an experiment.

## Idea and scope

The starting point is [Bertaina et al. (2024)](https://doi.org/10.1002/qute.202400032)
and their [open reference implementation](https://doi.org/10.5281/zenodo.10911121).
The TF-QKD calculation follows:

**Laser and fiber spectra → interference phase noise → usable key window →
duty cycle → SNS-AOPP or CAL key rate.**

The project reproduces that model, adds a delay-dependent classical round-trip
compensation calculation based on Williams, and provides a laboratory input and
reporting layer. The actuator pole and gain-selection choices are explicitly
identified engineering assumptions. This is a calculation tool; it does not
operate an optical link or generate cryptographic keys.

Mode-pairing QKD is a separate research module: phase increments from spectra,
their short-interval drift limit, and the base asymptotic Zeng key expression
using explicit decoy gains. Its documented limitations differ from TF-QKD.

## What you can do

- Supply model coefficients, supported datasheet quantities, or measured laser/fiber PSD files.
- Use a measured residual phase RMS instead of predicting its spectrum.
- Specify complete losses for each arm and separate detector descriptions; the scalar protocol's two detector projections are reported explicitly.
- Compare free fiber, classical round-trip compensation, dual-band stabilization and the ideal fiber-noise ceiling.
- Rank equipment changes by calculated gain, obtain conditional requirements for missing noise amplitudes, and optionally scan target-rate reach.
- Produce Markdown, HTML, JSON and plots through a CLI or local web page.

Apparatus parameters are not silently borrowed from a reference experiment.
If intrinsic misalignment or error-correction inefficiency is missing, the result
is labelled an upper model estimate. Published configurations are examples,
not a library of defaults for another laboratory.

## Quick start

Clone this repository, enter its directory, and use Python **3.11 or newer**:

```bash
python -m pip install -r requirements.txt
python -m tfqkd.lab_run examples/laboratory/stabilized_dual.toml --output results/my_run
python -m tfqkd.lab_web --port 8765
```

The web page is served locally at http://127.0.0.1:8765. Stop it with Ctrl+C.
Replace the example inputs with your installation's values; source citations
and units are included in each configuration.

Separate mode-pairing example:

```bash
python -m tfqkd.mp_key examples/mode_pairing/key_zhang202_given_gains.toml --output results/my_mp_key
```

## Validation and its scope

- Independent reimplementation of TF spectra, integration, window and duty: maximum relative discrepancy **1.978 × 10⁻⁷**. [Evidence](results/stage3_d/independent.json)
- Table I fast/reference regression: **4,097 / 65,537** points; maximum variance error **5.101 × 10⁻⁵**, below the **10⁻⁴** criterion. [Evidence](results/stage3_a/grid_regression.json)
- Held-out Zhou TF comparison: fit one intrinsic-error parameter at one length, predict the other; rate error **8–11%**, equivalent to **3.98–4.23 km** under the stated loss interpolation. This checks the protocol/loss calculation with a calibrated parameter, not the spectral phase forecast. [Rate evidence](results/out_of_sample/details.md), [distance evidence](results/stage3_d/validation_distance.json)
- At Zhou's **615.59 km** edge, the model yields zero against a published **0.32 bit/s**; no empirical correction is added.
- Mode-pairing spectrum-to-drift limit is checked analytically and numerically. Transferring another line's fiber coefficients overpredicts published spool drift by **1.6–2.8×**; the independent-section length law is outside the **conditional** fitted 95% interval. [Phase evidence](results/mp_phase/details.md), [length analysis](results/mp_spatial_analysis/details.md)

For the checked Bertaina configurations, classical compensation can lose its
benefit around **30–70 km per arm**, while dual-band stabilization captures most
of the calculated ideal fiber benefit. This is a configuration-dependent result,
not a general fiber-length limit. [Calculation evidence](results/release_v1_0/conclusion_checks.json)

## Run checks

```bash
python -m unittest discover -s tests
TFQKD_FULL_TESTS=1 python -m unittest discover -s tests
python scripts/validate/validate_stage3_a.py --calibrate
```

These commands work without article PDFs. Historical source-based T1–T7 scripts
also require local article files; see [SOURCES.md](SOURCES.md). Retained evidence
records the accepted T4 discrepancy rather than declaring that measurement reproduced.

## Known limits

The key models are asymptotic; finite-size effects are not included. The classical
scheme omits phase-measurement detection noise and its actuator is a one-pole
model. The fiber spectrum does not reproduce the measured 114 km curve's full
shape. Scalar protocol functions do not constitute a two-detector asymmetric
model. Mode-pairing spool drift is not a field-line default, and its analytical
decoy bounds can give zero where published experiments report a positive key.

See [tool and validation limits](docs/UNKNOWNS.md) and
[mode-pairing diagnostics](results/mp_zhu_diagnostics_v1_1/details.md).

## Repository contents and licensing

Code, tests, explicit examples, documentation, small validation records and
attributed CC BY reference data are included. Article PDFs, copied figures,
large scan outputs, private installation data, logs, caches and the original
research Git history are excluded.

[MIT](LICENSE) covers original contributions. Reference software/data and adapted
author functions retain their attribution and licence; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
