# Status — v1.1 public snapshot — 2026-10-05

## Available
- TF-QKD spectral forecasts and measured-phase estimates; SNS-AOPP and CAL.
- Common/independent lasers; free, classical and dual-band fiber schemes.
- Explicit coefficients, datasheet conversions, fitted/direct PSD, actuator CSV.
- Noise attribution, upgrade ranking, conditional requirements and optional reach.
- Local web form, CLI, Markdown/HTML/JSON reports and plots.
- Separate mode-pairing phase and asymptotic Zeng key modules with explicit decoy gains.

## Verified
- Independent TF implementation: maximum relative difference 1.978e-7.
- Table I fast/reference: 4097/65537 points; variance error 5.101e-5 < 1e-4.
- Zhou TF held-out: one calibrated e_d, 8–11% rate error; 3.98–4.23 km under stated loss interpolation.
- Clean snapshot: 115 quick tests passed (6 full-only skips); all 121 full unit tests passed.
- Example CLI and fast/reference calibration run successfully without PDFs.

## Limits
- Asymptotic keys; calibrated e_d also absorbs finite-key differences.
- Zhou TF edge 615.59 km: calculated zero, published 0.32 bit/s; no correction applied.
- Unknown e_d/f_EC: upper model estimate and conditional target-rate requirements.
- T4 shape mismatch; one-pole actuator; classical detection noise omitted; scalar detector projections.
- MP spool drift is not a field default; analytical decoy bounds and Zhu table-closure limits remain.

## Publication / next
- Code, tests, explicit examples, documentation and small evidence only; fresh Git history.
- Article PDFs, copied figures, large outputs, logs, credentials and caches are excluded.
- Original contributions: MIT; Bertaina reference/adaptations: attributed CC BY 4.0.
- Published: https://github.com/Quarzii/tfqkd-lab ; no further model development is scheduled.

## Run
```bash
python -m tfqkd.lab_run examples/laboratory/stabilized_dual.toml --output results/my_run
python -m tfqkd.lab_web --port 8765
python -m unittest discover -s tests
```
