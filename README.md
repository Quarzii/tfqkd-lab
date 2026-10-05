# TF-QKD and mode-pairing noise model

[Русский](README.ru.md) · [User guide](docs/USER_GUIDE.md) · [Validation](docs/VALIDATION.md) · [Sources](SOURCES.md)

Estimate laser and fiber phase noise in quantum key distribution, compare
stabilization schemes, and calculate asymptotic key rates from specified or
measured inputs.

**TF-QKD:** calculate phase-noise spectra, phase variance, the usable transmission
window, duty cycle, and SNS-AOPP or CAL key rate. Compare uncompensated fiber,
classical round-trip compensation, dual-band stabilization, and ideal fiber-noise
cancellation. Measured spectra and residual phase RMS are also supported.

**Mode-pairing QKD:** calculate phase increments between paired pulses and a
base asymptotic key bound from explicit detection gains. The phase calculation
and key-rate calculation are separate commands. A key estimate based on observed
gains is conditional on those observations.

## Get started

Requires Python **3.11+**. Run these commands in a terminal:

```bash
git clone https://github.com/Quarzii/tfqkd-lab.git
cd tfqkd-lab
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m tfqkd.lab_web --port 8765
```

On Windows, activate the environment with `.venv\Scripts\activate`.
Open <http://127.0.0.1:8765> in your browser. Stop the server with Ctrl+C.
The interface runs locally and offers TF-QKD input forms and reports.

For a TF-QKD calculation from a file:

```bash
python -m tfqkd.lab_run examples/laboratory/stabilized_dual.toml --output results/my_run
```

This writes Markdown, HTML, JSON, and, for spectral calculations, plots.
Replace example parameters with your installation's values; example coefficients
are calculation inputs, not equipment defaults.

For mode-pairing phase noise and key rate:

```bash
python -m tfqkd.mp_phase examples/mode_pairing/phase.toml --output results/my_mp_phase.json
python -m tfqkd.mp_key examples/mode_pairing/key_zhang202_given_gains.toml --output results/my_mp_key
```

## Choose your inputs

- **Noise prediction:** supply explicit laser/fiber coefficients or measured PSD files.
- **Measured phase:** supply residual phase RMS and its operating window to estimate the key rate at that measured condition.
- **Mode-pairing key rate:** supply detection gains, error gains, pairing parameters, and phase or drift inputs.

[The user guide](docs/USER_GUIDE.md) explains units, outputs, and how to select an
example. [Input reference](docs/UNIVERSAL_INPUTS.md) covers TF-QKD configuration
fields and CSV formats.

## Scope and accuracy

The TF-QKD model follows Bertaina et al.; classical compensation uses the
Williams round-trip model. The mode-pairing calculation uses the base Zeng
protocol and documented phase/decoy assumptions. See [sources](SOURCES.md).

Key rates are asymptotic: finite-size security analysis is outside the model.
Classical compensation omits phase-measurement detection noise, so its rate is
optimistic. Results depend on the supplied spectra, timing, and apparatus values;
coefficients from one route are not calibrated predictions for another.
Missing noise amplitudes or intrinsic-error parameters can produce a labelled
upper model estimate rather than an installation forecast.

Numerical consistency checks and published-experiment comparisons have different
scopes. The tool does not reproduce every published rate or measured spectrum.
See [validation and discrepancies](docs/VALIDATION.md) and
[model limitations](docs/UNKNOWNS.md).

## Check the installation

```bash
python -m unittest discover -s tests
TFQKD_FULL_TESTS=1 python -m unittest discover -s tests
python scripts/validate/validate_stage3_a.py --calibrate
```

These checks run without article PDFs. Paper-dependent reproduction checks are
explained in [SOURCES.md](SOURCES.md).

## License

[MIT](LICENSE) covers original contributions. Reference software, adapted author
functions, and measurement files retain their own licenses and attribution:
[third-party notices](THIRD_PARTY_NOTICES.md).
